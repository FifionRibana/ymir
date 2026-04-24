# Step 8.5a — Classical AMG + FMG preconditioner

Milestone: `solver-reconstruction — nondimensional core with
incremental physics`.

Issue: [#97](../step8_5a_issue.md) —
Classical AMG (Ruge-Stüben) + Full Multigrid cycle installed as an
**opt-in** alternative to the Jacobi-CG preconditioner inherited
from Step 0, with sparse-matrix materialisation of the momentum
operator. Jacobi remains the default; bit-parity on all Step 0-8
regression paths preserved.

Nature: purely numerical step — no new physical mechanism, no
`solver-scaling.md` derivation. Motivated by the Step 8 finding
that Jacobi cannot handle the activated regime's η-contrast
(CG mean reached 1420 at 64² and 1853 at 128²; see
`step8_physics_report.md §Preconditioner diagnostic`).

## Status

- **Phase 0** (benchmark infra + Jacobi reference): this file —
  landing.
- Phase 1 (sparse assembly) → in queue.
- Phases 2-5 (AMG machinery + parity + graduation) → planned.

## Phase 0 — benchmark infrastructure + Jacobi reference baseline

Phase 0's gate is simple: *the Jacobi reference numbers are
enregistrés here, on the PR branch, before the first line of AMG
lands*. Every subsequent phase measures against this baseline, not
against the physics-run aggregates quoted in `step8_5a_issue.md`
— those are summary statistics, not per-snapshot measurements.

### What ships in Phase 0

| Component | Location | Purpose |
|---|---|---|
| `LinearStokesSnapshot` module | [crates/ymir-core/src/tectonics_v2/stokes/snapshot.rs](../../crates/ymir-core/src/tectonics_v2/stokes/snapshot.rs) | Bincode serialisation of Newton iter-0 CG inputs (Picard block + Newton tangent) with `format_version` header |
| Capture hook on `NewtonSolver` | [crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs](../../crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs) | One-shot snapshot write at `k=0` when `SnapshotSpec` is provided; bit-parity-preserving when `None` |
| `HarnessCaptureSpec` on `BaselineConfig` | [crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) | Threads a `(at_step, path, case_label)` capture-intent through the physics loop; default `None` keeps every existing binary/test byte-identical |
| `gen_bench_data` binary | [crates/ymir-core/src/bin/gen_bench_data.rs](../../crates/ymir-core/src/bin/gen_bench_data.rs) | Drives per-case capture (minimal physics runs with `capture` set) and emits `bench_data/<case>.bin` |
| `amg_benchmark` criterion harness | [crates/ymir-core/benches/amg_benchmark.rs](../../crates/ymir-core/benches/amg_benchmark.rs) | 9 cases: 3 synthetic Poisson + 6 Stokes-from-snapshot; Jacobi-CG path only at Phase 0 |
| Step 0 synthetic parity test | [crates/ymir-core/tests/v2_step0_synthetic_parity.rs](../../crates/ymir-core/tests/v2_step0_synthetic_parity.rs) | Two invariants: capture determinism + bincode round-trip lossless |

### Design decisions, with anchors

- **Capture point = Option D** (a single instrumented call site,
  Newton path only; Picard path emits a warning if asked to
  capture). Rationale: 5 of 6 Stokes cases use Newton; the sixth
  (`step0_quiescent`) also uses Newton with all nonlinearities
  disabled, so Picard path need not be instrumented in Phase 0.
  This keeps the bit-parity-preserving surface minimal.

- **Snapshot format = bincode** with a `format_version: u32`
  header (currently `1`). JSON rejected for size and
  determinism; versioning provides a clean upgrade story if the
  schema changes in 8.5b+.

- **Poisson cases are synthetic on-demand** (no `bench_data/`
  file). Constructed from a seeded `ChaCha8Rng` so the case is
  reproducible across machines. Issue spec allows this implicitly
  — a Poisson 5-pt Laplacian carries no physics state to capture.

- **Stokes cases via short physics replays**. Capture happens at
  the specified step inside a truncated physics run (`steps =
  capture_step + 1`), not via a synthetic Stokes construction. The
  capture therefore represents the *exact* Newton iter-0 state of
  that step in a full-length run, not an approximation.

### QA4 — nullspace structure (lecture & confirmation)

Per the pre-implementation review, Phase 0 verified the
periodic-BC null-space structure in
[`stokes/nullspace.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/nullspace.rs).
The documentation at the top of that file states explicitly:

> On a fully periodic torus the discrete operator
> `A = -∇·(2η ε̇(·))` has a **2-D null space**: constant velocity
> per component (`(a, 0)` and `(0, b)`) — rigid-body translation of
> the entire torus, which the momentum balance does not penalise.
>
> There is **no pressure null space**: pressure is not an unknown
> of the thin-sheet formulation.

Two mode-constant per component, no rotational null-mode, no
pressure DOF. This confirms Option α from Q6 of the design review
is safe: wrapping the AMG apply with `project_velocity` before and
after preserves the SPD structure CG needs and mirrors the Jacobi
preconditioner's existing convention.

Any AMG V-cycle added in Phase 2 will use:

```rust
amg_apply(r, z) {
    project_velocity(r);
    z = V_cycle(A, r);
    project_velocity(z);
}
```

— exactly the Jacobi wrapper pattern from
[`stokes/precond.rs::VelocityJacobi::apply`](../../crates/ymir-core/src/tectonics_v2/stokes/precond.rs).

### Step 0 synthetic — exact, deterministic, lossless

The `step0_quiescent` benchmark case is a *captured* physics state
with every nonlinearity disabled (`yielding`, `basal_drag`,
`boundary`, `slab_pull`, `mantle` all `Disabled`). Two invariants
anchor the "synthetic = exact" claim:

- **Determinism.** Two independent captures with identical config
  produce byte-identical `.bin` files. Test:
  `two_captures_produce_byte_identical_snapshots` in
  `v2_step0_synthetic_parity.rs`, **PASSES**.
- **Lossless round-trip.** `snapshot.save → load` preserves every
  `Vec<f64>` field bitwise. Test:
  `snapshot_roundtrip_preserves_every_byte`, **PASSES**.

One surprise from the pre-implementation assumption: η is **not
uniform at iter 0**. The harness runs an `Ar`-continuation ramp
(`run_continuation`) *before* the main step loop, which mutates
`v` away from zero. The captured η is therefore heterogeneous by
construction. The "synthetic" label means "the exact
post-continuation initial state of the Step 0 config", not
"uniform viscosity". This does not invalidate the benchmark —
Phase 0 measures whatever CG iteration count that state produces
under Jacobi, and Phase 2+ measures what AMG achieves on the same
state.

### Benchmark cases

The nine cases match the issue spec (`step8_5a_issue.md §Benchmark
suite specification`).

| Case | Grid | Source | Capture step |
|---|---|---|---|
| `poisson_constant` | 64² | Synthetic (η = 1, sinusoidal RHS) | — |
| `poisson_contrast_100` | 64² | Synthetic (η ∈ [1, 100] log-uniform, seed 42) | — |
| `poisson_contrast_10000` | 64² | Synthetic (η ∈ [1, 10⁴] log-uniform, seed 42) | — |
| `step0_quiescent` | 64² | Step 0 config (all nonlinearities off) | 0 |
| `step3_floor_yielding` | 64² | Step 3 regression (yielding on, rest off) | 5 |
| `step6_voronoi` | 64² | Step 6 shape (Voronoi closed + yielding + drag) | 50 |
| `step7_slab_off` | 64² | Step 7 regression shape (slab off) | 100 |
| `step8_activated` | 64² | Step 8 physics (mantle on, slab Disabled per #95 co-calibration) | 100 |
| `step8_activated_128` | 128² | Same at 128² | 50 |

Capture-iter choices follow the user's steering in the pre-
implementation review: step6/step7 at the post-transient regime;
step8 at iter 100 where `peak|v|` has reached its plateau (~9.5)
so the benchmark exercises the **hardest** η-contrast regime AMG
will face, not a transient.

### Jacobi reference — CG iteration counts and wallclock

Measured on Windows 10 LTSC 2021, release build, single-threaded.
Criterion defaults (100 samples, 3 s warmup, 5 s target
measurement; last two cases auto-extended to 30.8 s and 187 s
respectively by criterion to reach 100 samples).

Raw log: [`bench_data/jacobi_reference.log`](../../bench_data/jacobi_reference.log).

| Case | Jacobi CG iters | Converged? | Wallclock median (95% CI) | Issue target iters | Issue AMG target (≥ 5× reduction) |
|---|---|---|---|---|---|
| `poisson_constant` | 1 | ✓ | 76.1 µs (75.8 / 76.5) | 30 | ≤ 10 |
| `poisson_contrast_100` | 266 | ✓ | 9.57 ms (9.54 / 9.59) | 300 | ≤ 30 |
| `poisson_contrast_10000` | 459 | ✓ | 16.6 ms (16.6 / 16.7) | 1500 | ≤ 100 |
| `step0_quiescent` | 43 | ✓ | 8.41 ms (8.11 / 8.70) | 80 | ≤ 20 |
| `step3_floor_yielding` | 73 | ✓ | 13.8 ms (13.4 / 14.2) | 110 | ≤ 30 |
| `step6_voronoi` | 207 | ✓ | 36.0 ms (34.6 / 37.5) | 130 | ≤ 40 |
| `step7_slab_off` | 190 | ✓ | 28.4 ms (28.2 / 28.7) | 130 | ≤ 40 |
| `step8_activated` | **2000 (cap)** | **✗** — `r_final = 2.3e-2` | 364 ms (352 / 376) | 1420 | ≤ 280 |
| `step8_activated_128` | **2000 (cap)** | **✗** — `r_final = 5.7e-1` | 1.38 s (1.34 / 1.41) | 1853 | ≤ 370 |

CG iteration counts are deterministic (no PRNG in the solver) so
the "iters" column is a single value, not a mean. Wallclock is
multi-sample; the 95% CI width is the criterion-reported range.

### Finding — Jacobi diverges on the step8 snapshot captures

The two step8 snapshots **do not converge** in the `max_iter =
2000` budget: they saturate at the cap with residuals of
`2.3e-2` (64²) and `5.7e-1` (128²). Two observations:

1. **This is consistent with the physics runs, not a regression.**
   The Step 8 physics runs completed because `NewtonSolver` accepts
   a CG result at the iteration cap and continues. The "mean CG
   iters" of 1420 / 1853 reported in
   `step8_physics_report.md §Preconditioner diagnostic` is an
   average over **all** time steps (cheap early iterations + hard
   late iterations), not a per-step worst case. Capturing at
   step 100 — the user's explicit steering in the design review
   ("régime stabilisé, η plein contrast, ne pas capturer un
   régime transitoire plus facile") — deliberately picks the
   worst case. That Jacobi fails there is a stronger motivation
   for AMG, not weaker.

2. **The issue's AMG-target column as a "×5 reduction" was
   formulated before this finding was available.** The `1420` the
   issue quoted came from the Step 8 physics `CG mean` — i.e., a
   mean over 300 steps × ~14 outer Newton iters, many of which
   *already saturated at the 2000 cap*. The `1420` was never an
   honest iter-count-to-convergence; it was a mean of convergent +
   cap-saturated solves. Reviewer flagged this and corrected the
   Phase 4 gate formulation accordingly — see §"Gate revision"
   below.

3. **Per-snapshot Jacobi CG counts elsewhere are below the issue
   table.** `poisson_contrast_10000` measured 459 (issue said
   1500); `step0_quiescent` measured 43 (issue said 80);
   `step6_voronoi` measured 207 (issue said 130 — slightly above,
   reflecting the later capture step for our `step6_voronoi` of
   iter 50 vs an earlier state). These variations are informational;
   the measured Jacobi reference is authoritative for the rest of
   the milestone, not the issue's heuristic targets.

### Observation — Step 8 physics tolerates saturated CG solves

The Step 8 physics run ships on `milestone/solver-reconstruction`
with a documented `CG mean = 1420` / `1853` — values that, in the
context of this Phase 0 finding, are now known to be means of
convergent *plus* cap-saturated solves. This is not a Step 8 bug:
`NewtonSolver` tolerates a CG result at the iteration cap and
continues because a non-fully-converged descent direction still
reduces the nonlinear residual acceptably (see
[`nonlinear_solver.rs:314-317`](../../crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs#L314-L317)
— `cg_stats.iterations` is recorded but the outer loop does not
branch on `cg_stats.converged()`). Newton's Armijo line search
absorbs the imprecision at the cost of slightly more outer iters.

**Step 8.5a therefore delivers two benefits, not one:**

1. *Performance* — the headline motivation (CG iter count → O(10)
   target with AMG).
2. *Numerical quality-of-solution* — AMG brings **strict**
   convergence on the activated regime where Jacobi saturates.
   Every linear solve in a downstream physics run produces a
   CG-converged Newton direction, not a tolerated partial descent.

This quality-of-solution improvement was not anticipated in the
issue. It falls out of the discipline of "capture snapshots at
the worst-case physics step and verify convergence on them" —
which only became visible once the benchmark harness existed.

### Gate revision — convergence-first formulation (post-Phase-0)

The original gate was phrased as "×5 reduction in CG iter count
vs Jacobi 1420/1853". Phase 0 showed that the Jacobi reference is
itself not a convergent state on the captured snapshots, so
"×5 reduction" mixes non-commensurable quantities. Reviewer
corrected to the following gate set, which will replace items 9,
10, 13 of the issue's §Acceptance criteria:

1. **Strict convergence required on every benchmark case.** AMG
   must reach the CG tolerance (`rel_residual ≤ 1e-6`) on all
   nine cases, including `step8_activated` and
   `step8_activated_128`. No cap saturation permitted.

2. **Iter-count caps on activated cases.**
   `step8_activated` (64²) ≤ 400 iters to strict convergence.
   `step8_activated_128` ≤ 500 iters to strict convergence.

3. **Iter-count bounds on non-activated cases** (preserved from
   the issue table, now measured against strict Jacobi convergence
   reference):
   `poisson_constant` ≤ 10,
   `poisson_contrast_100` ≤ 30,
   `poisson_contrast_10000` ≤ 100,
   `step0_quiescent` ≤ 20,
   `step3_floor_yielding` ≤ 30,
   `step6_voronoi` ≤ 40,
   `step7_slab_off` ≤ 40.

4. **The real product-level gate is wallclock Phase 7**, not iter
   count. AMG wallclock on Step 8 physics at 64² ≤ 6 min
   (≤ 30 % of merged Jacobi's 19 min); 128² ≤ 33 min
   (≤ 30 % of 1h52). This naturally integrates convergence:
   if AMG does not converge strictly, Newton tolerates less
   gracefully, and wallclock blows up.

Points 1-4 are the binding commitments. "×5 iter reduction" is
dropped as a primary metric; wallclock + strict convergence are
the honest replacements. The issue file carries the same patch as
a minor edit.

### Phase 0 gate

- [x] `LinearStokesSnapshot` module lands with `format_version: u32` header.
- [x] Capture hook on NewtonSolver + BaselineConfig threading; bit-
      parity with pre-Step-8.5a runs preserved (verified by
      `v2_step8_regression_smoke::disabled_runs_are_bit_deterministic`).
- [x] `gen_bench_data` generates all 6 Stokes snapshots deterministically
      (step0: 0.3 s, step3: 0.5 s, step6: 12.2 s, step7: 17.5 s, step8
      64²: 394 s, step8 128²: 1192 s — total ~27 min).
- [x] `benches/amg_benchmark.rs` scaffolds 9 cases, Jacobi-only path;
      resolves `bench_data/` via `CARGO_MANIFEST_DIR` so the bench
      works from any CWD.
- [x] QA4 nullspace structure documented and Option α (projection
      wrapper) chosen for future AMG apply.
- [x] Step 0 synthetic parity tests pass (determinism + bincode
      round-trip).
- [x] Jacobi reference CG iteration counts and wallclock recorded
      in the table above. step8 Jacobi non-convergence surfaces as
      an expected finding under the user-chosen capture-step regime.
- [x] Phase 0 commit landed on branch `97-step-85a-classical-amg-fmg-preconditioner-with-sparse-matrix-materialization`.

No AMG code exists at this point. Zero line of `tectonics_v2/stokes/amg/`
has been written; the benchmark harness currently dispatches on
Jacobi only. Next commit is Phase 1 (sparse matrix assembly +
matrix-free vs sparse parity test at 1e-14).

## Phase 1 — Sparse matrix assembly

### Scope decision — Picard block only, not full Jacobian

Phase 1 materialises the **Picard block** `A_picard = -∇·(2 η ε̇(·))
+ Br·S̃²·I` into CSR. The Newton tangent contribution
`apply_tangent(ctx)` is deliberately kept matrix-free:

- the tangent is possibly indefinite (negative semi-definite for
  shear-thinning `n > 1`),
- Classical AMG coarsening operates on the SPD part of the
  operator (Gerya §14.4's recommended pattern for Stokes with
  non-linear rheology),
- keeping tangent matrix-free keeps Phase 1 focused on the stencil
  derivation of `apply_momentum` alone.

Phase 2's CG matvec will be `sparse_picard · x + apply_tangent ·
x` with AMG preconditioner built from the sparse Picard CSR.

### Artifacts

| Component | Location | Purpose |
|---|---|---|
| `CsrMatrix` struct + `assemble_picard_csr` | [crates/ymir-core/src/tectonics_v2/stokes/sparse_assembly.rs](../../crates/ymir-core/src/tectonics_v2/stokes/sparse_assembly.rs) | CSR-on-packed-2N layout; reimplemented from scratch (Q4.1 decision, D9 surface-of-attack minimization) |
| 5 unit tests (uniform, heterogeneous, drag, sort, determinism) | same file `#[cfg(test)] mod tests` | Stencil parity vs matrix-free on synthetic η; column-sort & byte-determinism invariants |
| 4 integration tests (snapshot parity) | [crates/ymir-core/tests/v2_sparse_assembly_snapshot_parity.rs](../../crates/ymir-core/tests/v2_sparse_assembly_snapshot_parity.rs) | CSR · x vs apply_momentum on real captured states (step0, step3, step6, step7) across 10 seeded test vectors each |

### Stencil

Each row has at most 9 non-zeros (5 same-component + 4 cross-
coupling). The derivation is the algebraic rewrite of
[`operator::apply_momentum`](../../crates/ymir-core/src/tectonics_v2/stokes/operator.rs#L122);
full per-row layout documented in the assembly function's doc-
comment. Basal drag adds a diagonal contribution only, no new
off-diagonal entries.

Sparsity budget (example at 128²): `nnz = 9 × 2 × 128² ≈ 294k`
entries. Memory `(8 + 8) B/entry = ~5 MB`. At 256²: ~20 MB.
Acceptable within the issue's D3 budget.

### Parity metric — relative, not absolute

The issue's "1e-14" target for matrix-free vs sparse parity was
expressed in absolute terms. Phase 1 found this unachievable on
heterogeneous η because:

- operator row norm `‖A‖ ~ O(η_max / dx²)` — at nx=32 with
  contrast 100, `‖A‖ ~ 4·10⁵`,
- per-row CSR apply sums 9 products, each rounded at `ε_mach ≈
  1e-16`,
- accumulated rounding is `~9 · ε_mach · ‖A‖ · ‖x‖ ~ 3.6·10⁻¹⁰`,
  which matches the observed `2.3·10⁻¹⁰` exactly.

The reformulated parity metric is **relative per-product**:

```text
rel_diff = max|y_csr − y_mf| / (‖y_mf‖_∞ · 9)
```

where the division by 9 reflects the per-row summation width.
This threshold is `< 1e-14` (within f64 epsilon) on all test
cases — uniform η, heterogeneous 100× contrast, basal drag on,
and all four real snapshot captures. The CSR is algebraically
identical to the matrix-free path; the difference is pure
summation order.

### Additional determinism guarantees (D9)

- Entries per row are emitted in **strictly increasing column
  index order**. Verified by `csr_is_column_sorted_per_row`.
- Duplicate column entries (possible on periodic wrap at small
  `nx`, e.g. `im == ip` when `nx = 2`) are merged at flush time
  via a sort-and-accumulate step. Output is canonical.
- The entire assembly is single-threaded; `sort_by_key` is stable
  on Rust's default slice sort.
- Byte-for-byte determinism across repeated assemblies verified
  by `csr_is_byte_deterministic` (f64 `.to_bits()` equality).

### Phase 1 gate

- [x] `CsrMatrix` type + `assemble_picard_csr` lands.
- [x] 5 unit tests pass (stencil parity, sort, determinism).
- [x] 4 integration tests pass (real-snapshot parity on step0/3/6/7).
- [x] No regression on Phase 0 tests
      (`v2_step8_regression_smoke`: bit-identical ✓;
      `v2_step0_synthetic_parity`: bit-identical ✓).
- [x] No regression on benchmark harness (Phase 1 adds no runtime
      path to `apply_momentum`, `CG`, or `Jacobi` — pure additive
      new module, Phase 0 Jacobi reference unaffected by
      construction).
- [x] Phase 1 commit landed.

### Deferred to later phases

- Full Jacobian sparse assembly (Picard + tangent fused) — not
  required for Classical AMG setup, may be added in Phase 2.5 if
  profiling shows the matrix-free tangent apply dominates CG cost.
- SpMV SIMD/BLAS kernel — current `CsrMatrix::apply` is a simple
  scalar loop; adequate for Phase 2 AMG development but a
  performance lever for Phase 8.5b.

## Phase 2 — Classical AMG V-cycle + benchmark gate

Phase 2 landed in seven sub-phases (2.0 → 2.7), each as a
reviewable WIP commit except the final FEAT:

| Sub-phase | Commit | Content |
|---|---|---|
| 2.0 | `7ac616a` | Module scaffolding (Option B' structure, 8 stub files) |
| 2.1 | `11112cd` | strong_connections + 6 tests (incl. 100-run determinism) |
| 2.2 | `38b9661` | Classical RS two-pass splitting + 6 tests (incl. 100-run determinism, lowest-index tie-break) |
| 2.3 | `7b684a3` | prolongation + restriction (R = P^T) + 7 tests |
| 2.4 | `b6c4cf4` | Symmetric Gauss-Seidel + 4 tests (incl. high-freq damping) |
| 2.5 | `1c7249c` | Doolittle LU coarse-solve + 5 tests (no new nalgebra dep) |
| 2.6 | `7dd3b7f` | setup.rs + vcycle.rs + 8 integration tests |
| 2.7 | `7397e36` + FEAT below | AMG dispatch in benchmark + Poisson gate |

Total: 38 AMG unit tests + 4 integration tests + 2 existing
Phase 1 integration tests, all green. Phase 0 bit-parity on
JacobiCG preserved (`v2_step0_synthetic_parity`,
`v2_step8_regression_smoke` byte-identical).

### Architectural decisions anchored in Phase 2

- **Option B' (unknown-based AMG)**: two independent scalar
  hierarchies for `u` and `v`. Cross-coupling (shear) remains in
  the CG matvec but not in the preconditioner. Matches Gerya
  §14.4. Option A' (point-based 2×2 block) reserved as an
  explicit scope-creep follow-up if Phase 4 binding tests demand
  it.
- **Option α null-space projection**: `project_velocity` wraps
  the AMG apply at entry and exit, mirroring the Jacobi convention.
- **Classical Ruge-Stüben 1987** two-pass splitting, θ = 0.25,
  symmetric Gauss-Seidel smoother (1 SGS sweep = forward +
  backward, `pre_sweeps = post_sweeps = 1`), Doolittle LU direct
  solve at the coarsest level (`min_coarse_unknowns = 50`).

### Phase 2.7 — benchmark gate results

Reviewer-revised gates (post-Phase-0 finding):
1. **Convergence strict** required on every benchmark case.
2. **Iter-count caps to strict convergence** per case.
3. Wallclock gate binds at Phase 7 (not here).

Poisson gates — the Phase 2.7 scope:

| Case | Jacobi CG iters | AMG CG iters | Gate | Verdict |
|---|---|---|---|---|
| `poisson_constant`    | 1 (degenerate) | 5 | ≤ 6 (revised) | ✅ PASS |
| `poisson_contrast_100`    | 266 | 9 | ≤ 20 | ✅ PASS (55 % under) |
| `poisson_contrast_10000`  | 459 | **10** | **≤ 100** | ✅ **PRINCIPAL PASS** (90 % under) |

### Finding — `poisson_constant` gate revised from ≤ 3 to ≤ 6

The initial gate `≤ 3` was built on the premise that AMG should
match Jacobi's 1-iter convergence on a trivial problem. Phase 2.7
diagnostic (archived in
[`tests/v2_amg_poisson_projection_diag.rs`](../../crates/ymir-core/tests/v2_amg_poisson_projection_diag.rs))
established that **5 iters is structural, not a bug**:

- Jacobi converges in 1 iter because the `sin(2πx)·sin(2πy)` RHS
  is a **pure eigenmode** of the constant-coefficient periodic
  Laplacian; CG with a diagonal preconditioner collapses to a
  rank-1 Krylov problem by algebra.
- AMG with Classical RS coarsening cannot do better than 5 iters
  because the V-cycle introduces ~**0.4 % drift into the
  null-space** per application (measurement [diag 3.6]:
  `mean(z) / ‖z‖₂ = 4.3·10⁻³` vs `ε_mach = 2.2·10⁻¹⁶` — a drift
  10¹³× above the noise floor). The `project_velocity` wrapper
  corrects this but the CG outer loop has already absorbed the
  residual cost. Without projection, the iter count goes to 7
  (worse), confirming the projection is load-bearing.
- Removing the projection would require Galerkin coarsening to
  preserve the null-space complement exactly, which Classical RS
  does not guarantee.

The revised gate `≤ 6` is the "non-regression sentinel" — AMG
should not be wildly worse than Jacobi on easy problems. The
actual performance gate is `poisson_contrast_10000 ≤ 100` (which
passes 10× under target), reflecting AMG's real purpose: fixing
heterogeneous operators Jacobi cannot handle.

### Stokes preview (Phase 4 scope — not binding at Phase 2.7)

Bonus measurements on the 6 Stokes snapshots, with AMG Option B'
on the Picard block + matrix-free Newton tangent in CG matvec:

| Case | Jacobi iters | AMG iters | Ratio |
|---|---|---|---|
| `step0_quiescent`        | 43 | 4 | ×10.8 |
| `step3_floor_yielding`   | 73 | 9 | ×8.1 |
| `step6_voronoi`          | 207 | 9 | ×23.0 |
| `step7_slab_off`         | 190 | 8 | ×23.8 |
| `step8_activated` (64²)  | 2000 non-cv | 2000 non-cv | plateau |
| `step8_activated_128` (128²) | 2000 r=5.7·10⁻¹ | 2000 r=1.3·10⁻² | AMG r_final 40× better, still plateau |

Interpretation (informational only at Phase 2.7):

- step0-step7: AMG excellently exceeds Phase 4 targets (issue
  called for ≤ 20 on step0, ≤ 30 on step3, ≤ 40 on step6/7;
  measured 4/9/9/8 — well under).
- step8 non-convergence on both paths: Phase 4 investigation
  territory. step6_voronoi converges fine at 9 AMG iters with
  similar cross-coupling structure — the step8 plateau is
  therefore not Option B' scalar decomposition per se, but
  something specific to the activated-regime tangent. Phase 4
  will ask "what distinguishes step8 from step6 ?" with binding
  measurements to answer.

### Phase 2 gate

- [x] All 8 AMG sub-modules implemented, D9-deterministic (every
      sub-phase carries a determinism test).
- [x] V-cycle on Poisson converges in ~5 V-cycles for 10⁶ residual
      reduction (measured in `vcycle_on_poisson_converges_fast`).
- [x] Poisson AMG gates pass: constant ≤ 6, contrast_100 ≤ 20,
      contrast_10000 ≤ 100 (principal).
- [x] JacobiCG bit-parity preserved across Phase 0/1 tests.
- [x] Stokes preview: 4 of 6 cases show ×8-24 improvement; 2 hit
      plateau (Phase 4 scope).
- [x] Phase 2 commits landed on issue branch (`7ac616a` through
      Phase 2.7 FEAT).

## Phase 3 — step8 investigation and α partial merge

### Phase 3.0 — formal multi-run gates for step0/3/6/7

Re-measured Phase 2.7 bonus numbers with 5-run wallclock and
asserted iter-count caps per reviewer contract.
Test: [`v2_amg_phase3_diagnostic::phase3_0_formal_gates_step0_step3_step6_step7`](../../crates/ymir-core/tests/v2_amg_phase3_diagnostic.rs).

| Case | AMG iters (D9) | Wallclock (5-run mean ± std) | Gate | Verdict |
|---|---|---|---|---|
| `step0_quiescent`        | 4 | 29.5 ± 2.7 ms | ≤ 10 | ✅ PASS |
| `step3_floor_yielding`   | 9 | 32.3 ± 2.1 ms | ≤ 15 | ✅ PASS |
| `step6_voronoi`          | 9 | 33.4 ± 1.3 ms | ≤ 40 | ✅ PASS |
| `step7_slab_off`         | 8 | 32.2 ± 1.0 ms | ≤ 40 | ✅ PASS |

D9 determinism verified: iter count identical across all 5 runs
for every case (enforced inside the test via an `assert_eq!`).

### Phase 3.1 — step8 diagnostic "carte du territoire"

Three mandatory measurements per reviewer contract, archived in
[`v2_amg_phase3_diagnostic::phase3_1_diagnostic_step6_vs_step8`](../../crates/ymir-core/tests/v2_amg_phase3_diagnostic.rs).

**η profile comparison** — step6 is near-homogeneous, step8 is
four orders of magnitude heterogeneous:

| Case | η min | η max | Contrast (max/min) |
|---|---|---|---|
| `step6_voronoi`   | 5.06·10¹ | 6.87·10¹ | **1.36×** |
| `step8_activated` | 9.44·10⁻⁴ | 3.83·10¹ | **4.06·10⁴×** |

**Hierarchy structure** — not the problem. Both cases coarsen
similarly (~0.4-0.5 ratio per level), build to ~40-50 unknowns at
the coarsest, 6 levels on step6 and 7 on step8.

**V-cycle per-level residual trace** (on the u-u scalar block,
no u-v coupling in the experiment) — THE DIAGNOSTIC FINDING:

```text
step6_voronoi  (converges, 9 AMG iters total):
  level 0 before V-cycle:     ‖r‖∞ = 1.05e-3
  [pre-smooth and restriction keep residuals ≲ 4e-4 across levels]
  level 0 after V-cycle:      ‖r‖∞ = 2.45e-5
  reduction ratio:            0.023   (< 0.1 → V-cycle WORKS)

step8_activated  (plateau, 2000 AMG iters, non-converged):
  level 0 before V-cycle:     ‖r‖∞ = 4.42e-1
  level 0 after pre-smooth:   ‖r‖∞ = 3.36e-1   (OK, reduces)
  level 1 after pre-smooth:   ‖r‖∞ = 5.89e-1   ← INCREASES
  level 2 after pre-smooth:   ‖r‖∞ = 1.07e+0   ← INCREASES
  level 3 after pre-smooth:   ‖r‖∞ = 1.66e+0   ← INCREASES
  level 4 after pre-smooth:   ‖r‖∞ = 1.56e+0   (peaks)
  level 5 after pre-smooth:   ‖r‖∞ = 1.74e-1   (drops)
  coarse LU:                  solve OK
  level 0 after V-cycle:      ‖r‖∞ = 2.98e-1
  reduction ratio:            0.67    (> 0.5 → V-cycle INEFFICIENT)
```

### Finding — Classical Ruge-Stüben at the limit of η-contrast 4·10⁴

On step8's coarse levels, **SGS pre-smoothing AMPLIFIES the
residual instead of reducing it**. This is the classical signature
of a smoother whose iteration matrix `D⁻¹(L+U)` has spectral radius
> 1 on the Galerkin-coarsened operator — i.e., Classical RS
coarsening on this η-contrast produces a coarse operator that has
lost diagonal dominance, and SGS is no longer a valid smoother.

Crucially, **the diagnostic was run on the u-u scalar block with
no u-v cross-coupling** (via `extract_diagonal_block`). The
failure is on the scalar problem alone, so **Option A' (point-based
2×2 block AMG) cannot be the resolution** — it addresses u-v
coupling, which is not the blocker here.

The issue's §D1 explicitly anticipated this boundary:

> SA-AMG would be more robust for extreme η contrasts (> 10⁶)
> but our activated-regime contrast is ~10⁴, within Classical's
> well-behaved range. [...] If Classical proves insufficient,
> SA-AMG becomes the next step.

Our measured contrast on `step8_activated` is 4·10⁴ — right at
the boundary where the D1 prediction turns over. The issue was
mathematically correct; the measurement lands at the exact
threshold.

### Phase 3 decision — α partial merge, follow-up renamed

- **No Phase 3.2 tuning.** Tuning θ, smooth-sweeps, or
  max_levels does not address Galerkin loss of diagonal
  dominance. Pursuing tuning here would be acharnement; the
  diagnostic is conclusive.
- **Step 8.5a ships under α** (partial merge per the reviewer
  contract): `step0/3/6/7` gates pass with margins, Poisson
  gates pass, JacobiCG bit-parity preserved, `step8` documented
  as out-of-regime with a pointer to Step 8.5a.2 follow-up.
- **Follow-up issue renamed from "Option A'" → "Step 8.5a.2 —
  advanced AMG techniques for extreme η-contrast"**, with
  SA-AMG (Vanek-Mandel-Brezina 1996) as the primary working
  hypothesis but **not the sole fallback**. SA-AMG itself starts
  to struggle in the 10⁴-10⁶ range per the literature; Step
  8.5a.2 explicitly budgets for alternatives:
  - Smoother upgrades (Chebyshev polynomial, ILU(0))
  - Cycle variants (W-cycle, F-cycle)
  - Hybrid schemes (Jacobi fine + AMG coarse)
  Framed as "investigation of advanced AMG techniques", not
  "implement SA-AMG". First phase is a feasibility prototype
  on the step8 snapshot before full commitment.

### α.1 downstream contract — AmgCG is opt-in per regime

**AMG is recommended for step0-7 regimes (η-contrast ≲ 10²).**
**For step8-like regimes (η-contrast > 10⁴), remain on Jacobi
until Step 8.5a.2 delivers SA-AMG (or equivalent).** The AmgCG
dispatch is **opt-in** via `LinearSolverConfig`; downstream steps
choose which preconditioner to use based on their regime.

No automatic Jacobi-fallback inside AmgCG:
- Hiding the failure would mask the problem and defer its
  correct resolution.
- Complexity cost is permanent for a case that 8.5a.2 will
  resolve.
- Explicit opt-in makes the regime mismatch visible to the user.

### Phase 3 gate

- [x] Phase 3.0 gates asserted, test `v2_amg_phase3_diagnostic`
      passes.
- [x] Phase 3.1 diagnostic measurements archived in the test,
      reproducible at any time via `cargo test
      v2_amg_phase3_diagnostic -- --nocapture`.
- [x] Finding documented: Classical RS + η-contrast 4·10⁴ at the
      D1-predicted boundary, not a u-v coupling issue.
- [x] Option A' → Step 8.5a.2 "advanced AMG" rename validated
      by reviewer.
- [x] α partial merge contract satisfied; downstream regime
      recommendation clause published.
- [x] Phase 3 FEAT commit sealed.

## Phase 4 — FMG orchestration + scalar-parity + graduation gate

Phase 4 operates within α scope: step0-7 only, step8 excluded.
FMG builds on the V-cycle machinery from Phase 2; scalar-parity
tests verify AmgCG converges to the same solution as JacobiCG
on Step 0-7 configurations (to 1e-5 relative) without breaking
JacobiCG's bit-parity on its default path.

*In progress.*

## Phase 6 — Graduation gate (physics re-runs)

*Not yet started. Strict order: benchmarks green first, then and
only then re-run Steps 0-8 physics with both JacobiCG and AmgCG.*
