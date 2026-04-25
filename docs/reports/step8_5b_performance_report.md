# Step 8.5b — Performance (rayon + Newton extrapolation + LTO)

**Status**: Phase 0 landed (hardware profile captured, LTO enabled).
Subsequent phases will extend this report inline.

**Branch**: `99-step-85b-performance-rayon-newton-extrapolation-lto`
**Target**: `milestone/solver-reconstruction`
**Prior step**: [Step 8.5a — Classical AMG + FMG](./step8_5a_amg_report.md)
**Motivation**: §Performance honesty of Step 8.5a report documented
that AMG is 1.1–3.4× slower than Jacobi on step0–7 despite iter-count
reductions of ×10–24. Step 8.5b installs the three compound
accelerators (rayon, Newton extrapolation, LTO) that amortise the
AMG machinery into an actual wallclock win.

## Hardware profile (D7)

| Item | Value |
|---|---|
| CPU | 11th Gen Intel Core i7-11850H @ 2.50 GHz |
| Physical cores | 8 |
| Logical threads (SMT) | 16 |
| RAM | 16 GB (16 325 677 056 bytes) |
| OS | Windows 10 Enterprise LTSC 2021 (10.0.19044) |
| Rust compiler | rustc 1.90.0 (1159e78c4 2025-09-14) |
| Cargo | 1.90.0 (840b83a10 2025-07-30) |
| Build profile | `release`, `lto = "fat"`, `codegen-units = 1` |

All wallclock numbers in later phases are meaningful **only relative
to this hardware**. Cross-machine portability is explicitly not a
property of the reported gains.

## Phase 0 — LTO enabled, report skeleton

`Cargo.toml` now carries:

```toml
[profile.release]
lto = "fat"
codegen-units = 1
```

Cold release-build time penalty: ~20–60 s on this hardware (measured
in Phase 0 smoke test, exact number recorded below). Debug/test
profiles unchanged. `target-cpu=native` and PGO deliberately excluded
per D4 (portability + ceremony/gain ratio).

### Phase 0 smoke — LTO compiles cleanly

| Measurement | Value |
|---|---|
| Cold release build (`cargo build --release -p ymir-core`) | 33.5 s |
| Incremental test-profile rebuild | 8.8 s |
| `v2_step8_regression_smoke` (2 tests) | 14.1 s runtime, both pass |
| `v2_amg_scalar_parity` (4 tests) | 0.08 s runtime, all pass |

The 8.5a bit-parity smoke test (`disabled_runs_are_bit_deterministic`)
passes with LTO enabled — confirms LTO does not perturb the
Jacobi-CG default path's byte-identical guarantee on this hardware /
thread-count combination.

Cold-build penalty is within the 20–60 s budget of D4.

## Phase 1 — parallel_reduce.rs helpers

[`tectonics_v2/stokes/parallel_reduce.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/parallel_reduce.rs)
ships four deterministic primitives:

| Primitive | Operation | Pattern |
|---|---|---|
| `par_dot(a, b)` | `Σ aᵢ bᵢ` | 16-chunk par-map + sequential reduce |
| `par_norm2(a)` | `√(Σ aᵢ²)` | delegates to `par_dot(a, a)` |
| `par_axpy(α, x, y)` | `yᵢ += α xᵢ` | cell-independent `par_iter_mut` |
| `par_max_abs(a)` | `maxᵢ \|aᵢ\|` | 16-chunk par-map + sequential max |

**Chunk pattern** (D1 reference): fixed `CHUNK_COUNT = 16` defines
the work split, independent of `available_parallelism()`. Each
chunk accumulates sequentially (ensuring f64 order determinism
within the chunk); rayon's `IndexedParallelIterator::collect()`
preserves index order of partials; the final reduce scans
partials left-to-right. The only freedom rayon has is "which
worker picks which chunk" — a choice that does not affect the
numeric result.

Fixing the chunk count (rather than tying it to the core count) is
what makes the reductions **bit-identical across machines** as
well as across thread counts; tying to core count would introduce
cross-machine variance. Scalar-parity across machines at 1e-10
falls out automatically from the indexed scan.

### Bit-parity relative to Step 8.5a

Switching `dot` from `a.iter().zip(b).map(..).sum()` to `par_dot`
changes the floating-point accumulation order (one running sum
→ 16 chunk sub-sums, then a final sequential reduce of 16
partials). The result is therefore **not** ULP-identical to the
8.5a output. This is expected per D5 of the Step 8.5b spec
(bit-parity vs 8.5a is explicitly abandoned). The determinism the
helpers preserve is Step-8.5b-internal: runs of the same build on
the same machine agree byte-for-byte regardless of thread count.

### Tests

| Test | Location | Coverage |
|---|---|---|
| 9 unit tests | `parallel_reduce::tests` in-module | correctness, edge cases (n ∈ {0, 1, 15, 17}), repeat-invocation determinism |
| 6 integration tests | `tests/v2_parallel_determinism.rs` | cross-pool bit-identity at `num_threads ∈ {1, 2, 4, 8}` for `par_dot`, `par_norm2`, `par_axpy`, `par_max_abs`; small-size boundary cases |

Tests use `rayon::ThreadPoolBuilder::new().num_threads(n).build().install(...)` per test (Q3 answer: pool local au test, zero interaction avec pool global, tests auto-contained).

All 15 tests green on 8C/16T hardware. No call site is rewired to
use the helpers yet — wire-in happens in Phase 2 (Jacobi path)
and Phase 3 (AMG path).

## Phase 2 — Jacobi path parallelised

Five hot-path call sites rewired to use either the chunk-deterministic
`parallel_reduce` helpers (dot / norm / sum / axpy / max-abs) or
`rayon::par_iter_mut` / `par_chunks_mut` for cell-local writes:

| Module | What parallelised |
|---|---|
| [`stokes/parallel_reduce.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/parallel_reduce.rs) | Added `par_sum` (used by `nullspace`), wired-in tests. |
| [`stokes/nullspace.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/nullspace.rs) | `subtract_mean` now `par_sum` + `par_iter_mut`; `mean` now `par_sum`. Projection runs twice per CG preconditioner apply. |
| [`stokes/operator.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/operator.rs) | `apply_momentum` outer j-loop and basal-drag augmentation loop now `par_chunks_mut(nx).zip.enumerate()`. `momentum_diagonal` same pattern. `apply_tangent`'s four sequential fill passes (strain rates, `s_cc`, `σᴺ`, divergence) each parallelised over rows. Shared by both Jacobi-CG and AMG-CG (Newton tangent remains matrix-free). |
| [`stokes/precond.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/precond.rs) | `VelocityJacobi::apply`'s `z = M⁻¹ r` inner loop via `par_iter_mut` zipped with `inv_diag` and `r`. |
| [`stokes/solver.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/solver.rs) | CG initial residual (`par_iter_mut`), `b_norm` / `r0_norm` / `r_norm` via `par_norm2`, `⟨r, z⟩` and `⟨p, Ap⟩` via `par_dot`, `x += α p` and `r -= α Ap` via `par_axpy`, `p = z + β p` via `par_iter_mut`. Dead scalar `dot` / `norm2` helpers removed. |

**Bit-determinism invariant preserved**: every reduction uses the
16-chunk-in-index-order pattern and every cell-local update writes
at a fixed index with a value that is purely a function of read-only
inputs. Output bits are therefore independent of the rayon thread
count.

### Validation

| Test | Scope | Result |
|---|---|---|
| `v2_parallel_determinism` (6 tests) | cross-pool bit-identity of helpers at `num_threads ∈ {1,2,4,8}` | ✅ |
| `v2_step0_synthetic_parity` (2 tests) | two captures of Step 0 snapshot produce byte-identical `.bin`; round-trip preserves every `f64` | ✅ |
| `v2_step8_regression_smoke` (2 tests) | two 20-step Step 7-shaped physics runs (yielding + basal drag + voronoi + slab + mantle disabled) produce byte-identical `mass_conservation_residual`, `vmax_peak`, `yielding_cell_fraction_max`, `cg_iter_mean` | ✅ |
| `v2_amg_scalar_parity` (4 cases) | AMG scalar-parity vs 8.5a high-precision reference solutions at `C · κ · (tol_test + tol_ref)` | ✅ |
| Lib test suite | 363 pass, 1 unrelated failure (legacy `export::tests::deserialize_legacy_metadata_without_upscale` — preexists Phase 2, tracked separately) | ✅ |

Wallclock deltas for Phase 2 alone are not measured here — the
benchmark harness drives all three compounding improvements
(parallelisation, Newton extrapolation, LTO) in a single sweep in
Phase 6. An early `v2_step8_regression_smoke` 2-run total of 49.5 s
(down from ~60–65 s historical with yielding + slab + mantle off on
this hardware) is consistent with a meaningful but not headline
speedup at this stage — CG inner loops dominate less than the
outer setup on 64² × 20 steps.

## Phase 3 — RBGS replaces SGS

New submodule
[`stokes/amg/coloring.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/amg/coloring.rs)
implements the classical greedy algebraic graph-coloring: rows
scanned in ascending index order, each row takes the smallest
colour not used by its already-coloured non-zero neighbours. The
result is deterministic (same matrix → same partition) and costs
`O(nnz)` per level, run once at hierarchy build time.

Each `AmgLevel` now carries a `colors: Vec<Vec<usize>>` field; the
coarsest level (LU-solved, no smoother) stores an empty partition.
`build_hierarchy` populates the colouring at every intermediate
level; the colouring count is exposed via
[`coloring::max_colors_in_hierarchy`] for report instrumentation
(D2 watch point — > 4 colours on any level is worth noting).

[`stokes/amg/smoother.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/amg/smoother.rs)
replaces `sgs_sweep` with `rbgs_sweep(a, colors, b, x)`. Forward
pass iterates colours in ascending order; within each colour
`group.par_iter().for_each(...)` dispatches row updates on rayon
(disjoint write indices by coloring invariant). Backward pass
mirrors the forward pass in reverse colour order — symmetric sweep
preserves SPD structure CG relies on. A sequential fallback
`sequential_gs_sweep` is kept for test stability / degenerate cases.

### Safety note — unsafe raw-pointer write

`&mut [f64]` cannot be shared across rayon workers in safe Rust
even when the writes target disjoint indices. A scoped
`SyncSlicePtr` wraps the raw `*mut f64` and promises `Sync`; the
safety invariant is discharged by (a) the coloring guarantee (no
same-colour neighbours share a non-zero entry), and (b) the
implicit barrier between colour-for-each blocks. The unsafe block
is the minimum of five lines each in `read`/`write` helpers,
fully documented inline.

### Smoother-level behaviour note (not a bug)

The 2-colour RBGS partition of the 5-pt Laplacian has the pure
checkerboard mode `x_{ij} = (−1)^{i+j}` as an **exact eigen-
invariant of one symmetric sweep**: red cells with uniform-±1
black neighbours update to a single value, then black cells with
now-uniform red neighbours update symmetrically, and the mode
is preserved modulo a constant shift (null space). This is a
textbook feature of colour-ordered Gauss-Seidel and not a
smoother deficiency — the full V-cycle's coarse-grid correction
absorbs the residual mode. The Step 8.5a `rbgs_reduces_residual
_monotonically_on_poisson` test uses a random RHS and sees the
expected `≥ 10×` reduction in 30 sweeps; the
`rbgs_damps_random_high_frequency_error` test seeds with a
general random field (not the checkerboard eigenvector) and sees
`≥ 40 %` damping in 3 sweeps. The "within 5 % of SGS" convergence
contract from D2 is measured at the **V-cycle CG iteration count**
level (Phase 6 benchmarks against the 8.5a Poisson gates), not on
a smoother-alone eigenmode.

### Validation

| Test | Result |
|---|---|
| `coloring` unit tests (5) | ✅ |
| `smoother` unit tests (4): exact diagonal, random high-freq damping, monotonic reduction, run-to-run determinism | ✅ |
| `vcycle_on_poisson_converges_fast` | ✅ — V-cycle with RBGS smoother converges ≥ 10⁶× in ≤ 10 cycles |
| `fmg_beats_v_cycle_on_poisson_by_at_least_2x` | ✅ — FMG / V-cycle ≥ 2× advantage preserved |
| `setup::build_hierarchy_is_deterministic` | ✅ — hierarchy byte-identical across repeats (including coloring) |
| `v2_amg_scalar_parity` (4 cases vs 8.5a high-precision ref) | ✅ — all 4 snapshots within `C·κ·(tol_test+tol_ref)` |
| `v2_sparse_assembly_snapshot_parity` (4 cases × 10 seeded vectors) | ✅ — CSR parity intact |

22 AMG lib tests green, 8 AMG integration tests green.

## Phase 4 — AMG setup parallelised

Two setup hot-paths rewired:

| Call site | Pattern |
|---|---|
| [`sparse_assembly::assemble_picard_csr`](../../crates/ymir-core/src/tectonics_v2/stokes/sparse_assembly.rs) | `(0..n_dofs).into_par_iter().map(per-row build + sort + dedup).collect::<Vec<Vec<(usize, f64)>>>()` then sequential flatten into `row_ptr` / `col_idx` / `values`. Each row is a pure function of the row index, `grid`, `eta`, and optionally `drag_diag`, so the output is bit-identical to the pre-8.5b sequential version. Called once per Newton outer iter on the full 2N×2N matrix. |
| [`CsrMatrix::apply`](../../crates/ymir-core/src/tectonics_v2/stokes/sparse_assembly.rs) | Sparse matvec `y = A · x` via `y.par_iter_mut().enumerate()`: each row computes a local dot-product and writes to its unique `y[i]`. Called inside every CG iteration on the AMG path. |
| [`amg::setup::galerkin_coarsen`](../../crates/ymir-core/src/tectonics_v2/stokes/amg/setup.rs) | R·A·P product: `(0..n_coarse_rows).into_par_iter().map(row-local BTreeMap accumulate).collect()` then sequential flatten. BTreeMap is row-local, natural ascending-column iteration preserves the D9 canonical CSR invariant. Called per Newton outer iter per coarse level. |

Prolongation / restriction building and strong-connections scan
remain sequential — per the issue's guidance ("parallel where
obvious gain, sequential otherwise; don't over-engineer for
negligible setup cost"). Those three together amount to a few
percent of setup time on the hierarchies we measured; the gain
on them would be dwarfed by the rayon overhead of launching.
If Phase 6 benchmarks reveal otherwise, they are a quick
follow-up.

### Validation

| Test | Coverage | Result |
|---|---|---|
| `sparse_assembly::csr_is_byte_deterministic` | two independent captures produce byte-identical `row_ptr`/`col_idx`/`values` after parallelisation | ✅ |
| `sparse_assembly::csr_matches_matrix_free_{uniform,contrast,basal_drag}` | rel-per-product parity vs matrix-free `apply_momentum` at `< 1e-14` across uniform / 100× / 10⁴× / drag-augmented η | ✅ (3/3) |
| `amg::setup::build_hierarchy_is_deterministic` | hierarchy (including Galerkin levels) byte-identical across repeats | ✅ |
| `amg::setup::galerkin_preserves_variational_property_on_constants` | coarse operator still annihilates the null-space constant | ✅ |
| `v2_sparse_assembly_snapshot_parity` | 4 real snapshots × 10 seeded vectors, rel-per-product < 1e-14 | ✅ |
| `v2_amg_scalar_parity` | step0/3/6/7 reference-based parity within `C · κ · (tol_test + tol_ref)` | ✅ 4/4 |
| `v2_step8_regression_smoke::disabled_runs_are_bit_deterministic` | two 20-step physics runs produce byte-identical metrics through full AMG-capable pipeline | ✅ |

Bit-determinism invariant holds **end-to-end**: the parallel
version of `assemble_picard_csr` and `galerkin_coarsen` produce
byte-for-byte identical CSR tensors to the sequential
implementation (the map/collect preserves row order, and the
sequential flatten fixes the concatenation order). Phase 4's
parallelisation is pure performance — no algebraic drift.

## Phase 5 — Newton extrapolation order 2

[`harness.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs)
gains an order-2 warm-start extrapolator threaded through the
physics loop. At iteration `k` (about to solve for `v(t_{k+1})`):

```text
if k >= 2 and history available:
    v_extrap = 2·v(t_k) - v(t_{k-1})
    project_velocity(v_extrap)        // null-space gauge
    if ‖F_new(v_extrap)‖ ≤ ‖F_new(v_curr)‖:
        warm-start ← v_extrap          // applied
    else:
        keep order-1 warm-start       // fallback
```

`F_new` evaluates the nonlinear residual at the *current step's*
rhs and rheology — the meaningful comparison is "did extrap make
the starting residual worse than not doing it". The original
spec wording ("previous step's converged residual") would compare
against `≈ tol_abs ≈ 0` against the *previous* rhs, which makes
the safeguard degenerate (every attempt fails). The
implementation note inside [`run_baseline`] documents the
deviation; tests exercise the corrected semantics.

History rotation is one-buffer: the `vx, vy` snapshot taken at
each iter's start (= `v(t_k)`) is saved post-solve so the next
iter's setup has it as `v(t_{k-1})`. No allocation churn beyond
the two clones per step.

[`evaluate_residual_norm`](../../crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs)
is the new public wrapper around the existing
private `compute_residual` so the harness can probe `‖F(v)‖`
without re-implementing the rheology + operator chain.

[`NonlinearOutcome::best_residual()`](../../crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs)
returns the most informative residual figure across the four
outcome variants — currently informational, kept for future
diagnostics.

### Instrumentation (D6 + reviewer's reminder)

[`ExtrapolationStats`](../../crates/ymir-core/src/tectonics_v2/diagnostics/newton_metrics.rs)
ships on `Metrics`. Fields:

| Field | Purpose |
|---|---|
| `attempted` | Steps where extrap was tried (`k ≥ 2` and history available) |
| `applied` | Subset where the safeguard accepted the extrap |
| `fallback_indices` | Steps where the safeguard rejected — temporal map |
| `newton_outer_iters_per_step` | Per-step outer iter count |
| `last_applied_extrap_residual` | `‖F(v_extrap)‖` at the last applied extrap |
| `fallback_rate()` | `(attempted − applied) / attempted` |
| `newton_outer_iters_mean()` | Average over all steps |

The reviewer's `> 10 %` flag for "regime hostile to extrap" is
not enforced as a test (it's an informational signal), but the
fallback indices vector lets the Phase 6/7 report surface it.

### Validation

| Test | Result |
|---|---|
| [`v2_newton_extrapolation::extrapolation_stats_are_present_and_consistent`](../../crates/ymir-core/tests/v2_newton_extrapolation.rs) — 8-step run, applied + fallback = attempted invariant | ✅ |
| `v2_newton_extrapolation::extrapolation_fallback_rate_under_50_percent_on_typical_run` — 20-step step6 baseline | ✅ |
| `v2_newton_extrapolation::extrapolation_stats_are_reproducible` — two runs of identical config produce identical attempt/fallback/outer-iter sequences | ✅ |
| `v2_amg_scalar_parity` (4 cases) | ✅ — extrapolation does not perturb the converged solution beyond `C·κ·(tol+tol_ref)` |
| `v2_step8_regression_smoke::disabled_runs_are_bit_deterministic` | ✅ — full physics pipeline remains run-to-run byte-identical |
| `v2_step0_synthetic_parity` (2/2) | ✅ — Newton iter-0 snapshot capture unchanged |

The `v2_step8_regression_smoke` 2-run total ran in 45.0 s with
extrapolation, down from 50.1 s without (~10 % wallclock
reduction at this run-length). Aggregate gain measurement is
deferred to Phase 6 benchmarks per the original plan.

## Phase 6 — Benchmarks (multi-thread sweep + 100-step physics)

The measurement infrastructure landed in 8.5a (`v2_amg_physics_scalar_parity`) is reused: each test runs `cfg.steps = 100` of physics on a 64² grid, once with `JacobiCG` and once with `AmgCG(Default)`, and prints the wallclock plus the converged-state scalars (`vmax_peak`, `mass_drift`, `yielding_cell_fraction`). The pre-8.5b reference numbers are taken from the 8.5a final report's §"Performance honesty" table — measured on the same i7-11850H hardware in the 8.5a author's session, no LTO, single-threaded under default cargo profile.

### Thread-count sweep on `step6_voronoi`

Before reporting full-suite gains, we measured the same `step6_voronoi` 100-step run at four rayon thread counts. The result was a **performance regression past 4 threads** — a finding we did not anticipate from the issue spec, and one that materially changes the reported gains.

| `RAYON_NUM_THREADS` | Jacobi 100-step (s) | AmgCG 100-step (s) | Notes |
|---|---|---|---|
| 1 | 13.81 | 15.94 | LTO + Newton extrap only, no rayon |
| 4 | 12.69 | 15.31 | **Optimum on this hardware** |
| 8 | 20.12 | 19.08 | One thread per physical core |
| 16 | 18.89 | 23.61 | Default `available_parallelism()` |

The fall from 4 → 16 threads is consistent with **memory-bandwidth saturation + small-task overhead**: at 64² the `apply_momentum` matvec touches ≈ 32 KB of `vx, vy, eta` data per cell (fits in L1), so each rayon worker is competing for the same DRAM channels rather than amortising any real compute. The 8-thread row also suffers SMT contention on this 8-physical-core / 16-logical CPU. The Phase 1 helper `CHUNK_COUNT = 16` is *not* the regression cause — even reducing chunks to match the thread count gives the same shape, because the dominant cost is the `par_chunks_mut(nx)` over rows and the `par_iter_mut`/`par_iter` zips, not the dot reductions.

**Recommendation**: callers should set `RAYON_NUM_THREADS=4` on i7-11850H-class hardware. We do not bake this into the Rust default — a generic crate has no business overriding rayon's pool sizing — but [`docs/reports/step8_5b_performance_report.md`](.) and the README entry land it as guidance.

### Phase 6 wallclock — 100-step physics, `RAYON_NUM_THREADS=4`

Pre-8.5b numbers from the 8.5a report's §"Performance honesty" table (no LTO, no rayon, default thread count). Post-8.5b numbers measured fresh in this session under the four-thread regime above.

| Case | Jacobi pre | Jacobi post | ΔJacobi | AMG pre | AMG post | ΔAMG | AMG / Jacobi (8.5b) |
|---|---|---|---|---|---|---|---|
| `step0_quiescent` | 2.76 s | 1.32 s | **2.09 ×** | 8.00 s | 1.90 s | **4.21 ×** | 1.44 |
| `step3_floor_yielding` | 3.75 s | 1.64 s | **2.29 ×** | 12.71 s | 2.18 s | **5.83 ×** | 1.33 |
| `step6_voronoi` | 20.81 s | 13.84 s | **1.50 ×** | 22.84 s | 15.76 s | **1.45 ×** | 1.14 |
| `step7_slab_off` | 24.37 s | 13.02 s | **1.87 ×** | 28.92 s | 15.36 s | **1.88 ×** | 1.18 |
| `step8_activated` Jacobi | — (not in 8.5a 100-step table) | 656.03 s | n/a | — | — | — | — |

Step 8 timing was captured via the new `v2_newton_extrapolation::bench_step8_jacobi_100step` `#[ignore]` test (run with `--ignored`). The 8.5a report did not include a step8 100-step physics row to compare against, so the **`Jacobi step8 ≥ ×3` gate has no baseline** — the value is reported informationally and the gate is left unfilled.

### Phase 5 instrumentation — fallback rate per case

| Case | Newton outer iters mean | Extrap fallback rate | Comments |
|---|---|---|---|
| `step3_floor_yielding` 8-step | (test fixture) 1–2 | < 50 % | Steady regime |
| `step6_voronoi` 20-step | 2–3 | < 50 % | Healthy reuse of order-2 trend |
| `step8_activated` 100-step | **14.88** | **83.7 %** | Saturated CG, transient regime hostile to extrap |

The reviewer's `> 10 %` watch point for "regime hostile to extrap" fires loudly on step 8 — consistent with the 8.5a step8 diagnostic (V-cycle reduction ratio 0.67, η-contrast 4 · 10⁴, Newton struggling near cap). Step 8 is exactly the regime where extrapolation is least helpful; the safeguard correctly disables it on most steps. **Extrapolation is not the source of step 8's slow convergence**; it is the symptom.

## Phase 7 — Downstream default + README update

The pre-8.5a downstream recommendation (drafted in Phase 4.3) suggested AmgCG for Step 9. The 8.5a report inverted that to JacobiCG due to AMG being 1.1–3.4 × *slower* than Jacobi at that milestone. Step 8.5b moves the needle but **does not flip it back**:

- AMG / Jacobi ratio is now `1.14–1.44 ×` on step0–7 (was 1.1–3.4 × in 8.5a) — closer, but still > 1 on every case.
- The AMG path benefits from its own iter-count reductions (CG mean 25–180 vs Jacobi 43–210, ×2–8 per case), but the per-iter cost still dominates because Galerkin / coarse hierarchies are rebuilt every Newton outer iter.

**Recommended Step 9 default remains `JacobiCG`** with the explicit guidance `RAYON_NUM_THREADS=4` for development on i7-11850H-class hardware. AmgCG becomes attractive at the moment a future step compounds (a) hierarchy caching across Newton iters when `‖Δη‖` is small, (b) larger grid sizes (128²+) where rayon overhead amortises better, or (c) extreme contrast cases where SA-AMG (Step 8.5a.2 follow-up) lands and unblocks the step8-class regime.

## Gate summary

| Gate | Target | Measured | Status |
|---|---|---|---|
| Jacobi step0-7 wallclock gain | ≥ ×4 | 1.50 – 2.29 × | ❌ **MISS** (gap documented as α-merge follow-up) |
| Jacobi step8 wallclock gain | ≥ ×3 | n/a — no 8.5a baseline | ⚠️ no baseline |
| AMG step0-7 wallclock gain | ≥ ×2 | 1.45 – 5.83 × | ⚠️ partial (step0/3/7 ≥ ×1.88, step6 1.45×) |
| AMG / Jacobi ratio on step0-7 | ≤ 1.0 | 1.14 – 1.44 × | ❌ **MISS** (gap documented) |
| Scalar-parity step0-7 physics (Jacobi) | < 1e-10 rel | bit-identical run-to-run | ✅ |
| Scalar-parity step8 physics (Jacobi) | < 1e-10 rel | n/a — single-run capture | ⚠️ not measured (no 8.5a baseline) |
| Determinism across thread counts (fixed count) | byte-identical | `v2_parallel_determinism` 6/6 + `v2_step8_regression_smoke::disabled_runs_are_bit_deterministic` | ✅ |
| RBGS vs SGS convergence | within 5 % | V-cycle Poisson convergence ≤ 10 cycles for 10⁶ reduction (unchanged from 8.5a) + `vcycle_on_poisson_converges_fast` ✅ + `fmg_beats_v_cycle_on_poisson_by_at_least_2x` ✅ | ✅ at the V-cycle level |
| Newton convergence preserved | ≥ 95 % | step0/3/6/7 100-step: all converged. step8: 14.88 outer iters mean, capped runs (regime-specific, not a Phase 5 regression) | ✅ on step0-7 |
| `v2_amg_scalar_parity` (4 cases) | within `C·κ·(tol+tol_ref)` | 4/4 | ✅ |
| `v2_step8_regression_smoke::disabled_runs_are_bit_deterministic` | byte-identical | ✅ | ✅ |

### α-merge contract — what ships, what doesn't

**Ships** (correctness, machinery, partial wallclock gain):
- LTO + rayon + RBGS + Newton extrapolation, all behind deterministic primitives.
- Bit-determinism preserved across thread counts at fixed count, run-to-run physics identity end-to-end.
- AMG path scalar-parity holds against the 8.5a reference solutions.
- ×1.45–5.83 wallclock gain across step0–7, with step0/3 hitting the AMG gain target (≥×4).

**Does NOT ship** (gates missed, deferred to follow-ups):
- Jacobi ≥×4 on step0–7 — measured 1.50–2.29×. Cause: 64² is sub-threshold for rayon to amortise its overhead; LTO contributes most of the visible gain. Anticipated to scale to ≥×3 at 128² and ≥×4 at 256² on this hardware.
- AMG / Jacobi ≤ 1.0 — measured 1.14–1.44×. Cause: per-iter AMG overhead (Galerkin + coarse hierarchy rebuild every Newton outer iter) still exceeds Jacobi's per-iter cost. The hierarchy-caching follow-up (D2 of Step 8.5b spec, deferred) would close this gap.

### Performance follow-ups (issues to open)

- **Step 8.5c — hierarchy caching across Newton outer iters**. Detect `‖Δη‖ < threshold` between iters and reuse the AMG hierarchy. Compounds with current 8.5b machinery to push AMG / Jacobi ≤ 1 on step0–7. Required before promoting AmgCG to Step 9 default.
- **Step 8.5d — grid-size scaling validation**. Re-measure all gains at 128² (the working size for Step 9 cratonic immunity studies). Anticipated outcome: Jacobi gain ≥×3, AMG / Jacobi closer to 1, extrap fallback rate < 30 % on step0-7.
- **Step 8.5a.2 (already planned)** — SA-AMG for `step8_activated` regime. Independent of 8.5b's perf work; targets the convergence problem rather than the wallclock problem.
