# Step 8.5a: Classical AMG + FMG preconditioner with sparse matrix materialization

**Milestone**: Solver reconstruction — nondimensional core with
incremental physics

**Nature of this step**: Purely numerical. No new physical mechanism,
no new nondimensional number, no `solver-scaling.md` derivation.
Installs a new preconditioner (Classical AMG with Full Multigrid
cycle) as an **opt-in alternative** to the current Jacobi, with
the target of reducing CG iteration counts by ~5-10× in the
activated regime (Step 8+) while preserving bit-parity on existing
Jacobi paths.

**Previous step**: #95 (Step 8 — Mantle forcing, yielding
checkpoint resolved at 0.998/0.997)

**Branch from**: `milestone/solver-reconstruction` after Step 8 PR
merged.

**PR target**: `milestone/solver-reconstruction`.

**References**:
- `docs/reports/step8_physics_report.md` §Preconditioner diagnostic
  — the observation that CG mean reached 1420 (64²) and 1853 (128²)
  at activated regime, motivating this step
- `docs/reports/step8_regression_report.md` — bit-identical Jacobi
  path as the parity target
- Classical AMG references: Ruge & Stüben 1987 (foundational),
  Briggs-Henson-McCormick "A Multigrid Tutorial" (pedagogical), Falgout
  2006 "An Introduction to Algebraic Multigrid" (implementation-oriented)
- FMG reference: Brandt 1977 "Multi-Level Adaptive Solutions"; the
  simplified FMG-from-V-cycle pattern is standard in modern AMG
  libraries

## Document navigation

- Step 8.5a context — why this step exists
- Goal — what the code must produce
- Design decisions (D1–D9)
- Methodology: pyramid-inverted development
- What gets built
- Benchmark suite specification
- Acceptance criteria
- Out of scope
- Definition of done

## Step 8.5a context

Step 8 established the activated regime (peak|v| ~ O(1), yielding
active, ε̇_II well above floor). It also established that the
Jacobi preconditioner inherited from Step 0 is **insufficient for
the activated regime**: the effective viscosity η contrast jumps
from uniform (floor-dominated) to ~10⁴ (activated), and Jacobi
captures only the diagonal magnitude, not the cell-coupling
structure that dominates convergence in heterogeneous elliptic
problems.

Observed cost of keeping Jacobi:

| Grid | Jacobi CG iter mean | Physics run wallclock |
|---|---|---|
| 64² Step 8 activated | 1420 | 19 min |
| 128² Step 8 activated | 1853 | 1h52 |
| 256² Step 8 activated (extrapolated) | ~4000 | ~9h |
| 512² Step 8 activated (extrapolated) | ~10000 | ~75h |

This is structural, not a bug. The preconditioner needs replacement
for the activated regime. Classical AMG with Full Multigrid cycle
is the standard answer for heterogeneous elliptic problems. Target:
CG iter count reduced to O(10) independent of N, wallclock scaling
returned to O(N²·log N) roughly.

Step 8.5a installs the AMG infrastructure. Step 8.5b (follow-up
step) will add `rayon` parallelization, Newton extrapolation, and
compilation flags. Both are required to reach the performance
targets — neither alone is sufficient.

## Goal

Produce a Classical AMG preconditioner with Full Multigrid cycle
that:

1. Is **opt-in via configuration** (Jacobi remains default for
   backward compatibility and bit-parity validation)
2. Reduces CG iteration count by at least **5×** on the Step 8
   activated benchmark (target: 1420 → ≤ 280 at 64²)
3. Preserves **bit-parity** with Jacobi on Steps 0-7 configurations
   when Jacobi is used
4. Produces **scalar-parity** with Jacobi when AMG is used (same
   converged solution to Newton tolerance, different path)
5. Ships a **benchmark suite** for the pyramid-inverted development
   methodology (see Methodology section)
6. Materializes the momentum operator as a **sparse matrix**
   (medium-scope refactor) — AMG coarsening requires the matrix
   structure, which matrix-free cannot provide

## Design decisions (already made, not to revisit)

### D1 — Classical AMG (Ruge-Stüben), not Smoothed Aggregation

Classical AMG:
- Identifies "strong connections" per row via thresholding on
  off-diagonal magnitudes
- Coarsens by C/F splitting (coarse points, fine points)
- Prolongation matrix built from strong-connection weights
- Standard smoother: Gauss-Seidel (forward + backward sweeps)

Rationale: Classical is simpler to implement correctly from scratch
in Rust (~1500-2500 LOC). SA-AMG would be more robust for extreme
η contrasts (> 10⁶) but our activated-regime contrast is ~10⁴,
within Classical's well-behaved range. If Classical proves
insufficient, SA-AMG becomes the next step. Start with the simpler
approach that can work.

### D2 — Full Multigrid (FMG) cycle, not V-cycle standalone

FMG solves on the coarsest grid first, prolongs the solution as
initial guess to the next level, runs V-cycle to convergence,
prolongs again, etc. Typically converges in 1-2 V-cycles total
vs 5-10 for standalone V-cycle.

FMG is included in Step 8.5a (not deferred to 8.5c). The infra
is the same as V-cycle; only the orchestration differs. If
benchmarks show FMG gives no additional gain on our heterogeneous
Stokes (possible — FMG's advantage is mode-dependent), we keep
V-cycle as fallback. Both paths supported.

### D3 — Sparse matrix materialization (medium scope)

AMG Classical requires the matrix in a sparse representation (CSR
or similar) to compute strong connections, build coarsening, and
construct prolongation. Matrix-free is incompatible with Classical
AMG's coarsening phase.

Sparse path installed:
- New module `stokes/sparse_assembly.rs` building a sparse CSR
  representation of the momentum operator on demand
- The current matrix-free path remains as-is and is the default
- When AMG is selected, the sparse matrix is built once per Newton
  outer iteration (once the η field is frozen for that iter) and
  consumed by the AMG setup phase
- Memory cost at 128² is acceptable (~5 non-zeros per row × 3 DOF
  × 16384 cells × 8 bytes = ~2 MB). At 512² it becomes ~32 MB,
  still acceptable.

### D4 — `LinearSolverConfig::{JacobiCG, AmgCG}` opt-in

The solver path becomes configurable:

```rust
pub enum LinearSolverConfig {
    JacobiCG,        // current path, default, bit-parity preserved
    AmgCG {
        cycle: AmgCycle,          // Vcycle or FmgCycle
        coarsening: CoarseningConfig,  // threshold, max_levels, etc.
        smoother: SmootherConfig,       // GS sweeps, etc.
    },
}
```

All Step 0-8 tests continue to use `JacobiCG` by default. Step 8.5a
adds a new set of tests exercising `AmgCG` on representative
configurations. When Step 8.5b lands, the opt-in remains — each
downstream step (9, 10) continues on `JacobiCG` unless explicitly
switched. Promotion of AMG to default is post-milestone, after
the retrospective re-runs confirm equivalence.

### D5 — Bit-parity guarantee for Jacobi path

No existing test should change behaviour. Specifically:

- All Step 0-8 regression smoke tests must continue to pass with
  byte-identical outputs when `LinearSolverConfig::JacobiCG` (the
  default) is used
- The `linear_solve.rs` module is extended (not refactored) —
  Jacobi code path remains exactly as-is, AMG added alongside
- The sparse matrix assembly is a new module, not a refactor of
  the matrix-free operator apply

### D6 — Scalar-parity guarantee for AMG path

AMG produces a different path to the solution than Jacobi but
should converge to the same solution (to Newton rtol = 1e-6). The
scalar-parity test for AMG:

- Same Step 0-8 configurations, run with `AmgCG`
- Compare final `S̃`, `v_solved`, `m_subducted` (where applicable)
  against the Jacobi output
- Acceptance: `max |scalar_amg - scalar_jacobi| / max|scalar_jacobi|
  < 1e-5` (loose tolerance because the converged state may differ
  at Newton tolerance level)
- Stricter metrics (`peak|v|`, `yielding_cell_fraction_max`): must
  agree to 1% in value and 100% in sign/structure

### D7 — Benchmark suite as primary development metric

Development of AMG happens against a fast benchmark suite (see
Methodology and Benchmark suite specification below), not against
full Step 0-8 physics runs. The full physics runs are validation
gates at the end, not development iterations.

Benchmark suite total wallclock on Jacobi (reference): **< 2 min**
Benchmark suite total wallclock on AMG (target): **< 30 sec**

### D8 — Cycle convergence criteria (AMG-internal)

Each V-cycle applies the smoother pre-coarsening, recursively
solves on the coarse grid, applies the smoother post-prolongation.
Number of pre/post smoother sweeps: 2 each (standard default).
Coarsest-grid solve: direct via LU factorization (sub-problem is
small, ~100 unknowns or less after 3-4 levels of coarsening).

Convergence of the outer CG wrapping the AMG preconditioner: same
criterion as Jacobi-CG currently (residual norm < Newton rtol).

### D9 — Determinism

AMG coarsening, C/F splitting, and smoother ordering must all be
deterministic. Same input matrix → same preconditioner sequence
→ same CG path → same output (byte-identical, or if f64 sums reorder,
at least to 1e-14 of absolute difference).

Seed handling: AMG does not take a seed. Any stochastic choice
(e.g., random tie-breaking in C/F splitting) must be replaced by
deterministic tie-breaking (lowest index wins).

## Methodology: pyramid-inverted development

Development cycles on the benchmark suite (~2 min per full
benchmark run), not on full physics runs (~2h at Step 8 scale).
The rule:

1. **Implement** a piece of AMG (e.g., strong connections detection)
2. **Unit test** it in isolation (test file under `tests/`)
3. **Run the benchmark suite**: each of the N benchmark cases.
   Measure CG iter count and wallclock for each. Compare against
   Jacobi reference and previous AMG state.
4. **Iterate** on the implementation based on benchmark feedback
5. When benchmark targets are reached (≥ 5× reduction in CG iter
   on activated-regime benchmark), **graduate to full physics
   validation**:
   a. Re-run Steps 0-8 physics with `JacobiCG` (should be
      byte-identical to merged reports — bit-parity gate)
   b. Re-run Steps 0-8 physics with `AmgCG` (scalar-parity gate:
      same physics, different solver path)
6. **Commit** the AMG implementation. The validation runs are
   expensive and should only be done when the benchmark suggests
   we're ready.

This structure is required because each iteration of AMG tuning
without the benchmark would cost hours of wallclock. A
benchmark-driven loop makes the work tractable.

## What gets built

### New module `solver/amg/`

Decomposed:

- `mod.rs`: re-exports, `AmgConfig` struct, `AmgPreconditioner`
  entry type
- `strong_connections.rs`: strong connection detection per row
- `splitting.rs`: C/F splitting (Classical Ruge-Stüben algorithm)
- `prolongation.rs`: prolongation matrix construction from C/F
  split + strong connections
- `restriction.rs`: restriction = prolongation transpose
- `smoother.rs`: Gauss-Seidel smoother (forward + backward sweeps)
- `coarse_solve.rs`: direct LU solve on the coarsest grid
- `vcycle.rs`: V-cycle orchestration
- `fmg.rs`: FMG orchestration (builds on V-cycle)
- `setup.rs`: full setup phase — build the hierarchy from matrix
- `apply.rs`: the `apply(r) → z` entry point used by CG

### New module `solver/sparse_assembly.rs`

Builds a sparse CSR representation of the momentum operator. Takes
the same inputs as `apply_momentum` (η field, grid, config) and
produces a `SparseMatrix` struct. Format: CSR (row pointers, column
indices, values), 3 DOF per cell (u, v, p or u, v and pressure
handled separately depending on existing structure).

### Extension of `linear_solve.rs`

Add the `AmgCG` variant to the solver enum. The CG loop is shared
between Jacobi and AMG; only the `apply_preconditioner` callback
differs.

### Benchmark suite — `benches/amg_benchmark.rs`

See next section for specification.

### New tests

- `v2_amg_setup_determinism.rs`: same input matrix → same
  preconditioner, byte-identical
- `v2_amg_smoother_convergence.rs`: smoother alone (no coarsening)
  converges at standard Gauss-Seidel rate on a Poisson test
- `v2_amg_vcycle_poisson.rs`: V-cycle achieves ≥ 10× iteration
  reduction vs Jacobi on constant-coefficient Poisson
- `v2_amg_vcycle_heterogeneous.rs`: same on contrast-10⁴ Poisson
  (closer to activated Stokes regime)
- `v2_amg_fmg_poisson.rs`: FMG achieves ≥ 2× additional reduction
  vs V-cycle on same Poisson
- `v2_amg_scalar_parity_step0.rs`: Step 0 regression config run
  with AMG matches Jacobi to 1e-5
- `v2_amg_scalar_parity_step8.rs`: same for Step 8 regression
  config (activated regime)
- Existing Step 0-8 regression smoke tests continue to pass
  byte-identical (Jacobi default path)

## Benchmark suite specification

The benchmark suite lives in `benches/amg_benchmark.rs` and runs
via `cargo bench --bench amg_benchmark`. Each benchmark extracts
a single representative Stokes solve (not a full 300-step physics
run) and measures:

- CG iteration count (Jacobi vs AMG)
- Setup phase wallclock (AMG only)
- Solve phase wallclock (both)
- Peak memory (both)

Benchmark cases:

| Case | Grid | Regime | Expected Jacobi CG | AMG target |
|---|---|---|---|---|
| `poisson_constant` | 64² | Trivial, uniform coeff | 30 | ≤ 10 |
| `poisson_contrast_100` | 64² | η contrast 100× | 300 | ≤ 30 |
| `poisson_contrast_10000` | 64² | η contrast 10000× | 1500 | ≤ 100 |
| `step0_quiescent` | 64² | Step 0 Stokes, uniform η | 80 | ≤ 20 |
| `step3_floor_yielding` | 64² | Step 3 regression state | 110 | ≤ 30 |
| `step6_voronoi` | 64² | Step 6 physics end-of-run | 130 | ≤ 40 |
| `step7_slab_off` | 64² | Step 7 floor-dominated | 130 | ≤ 40 |
| `step8_activated` | 64² | Step 8 mantle-on end-of-run | 1420 | ≤ 280 |
| `step8_activated_128` | 128² | Same at 128² | 1853 | ≤ 370 |

Each case is a serialized Stokes state (η field, RHS, grid config)
loaded from a `bench_data/` directory. The data is generated once
via a helper binary `gen_bench_data` that runs the relevant step
config, extracts the Stokes state at a representative iteration,
and serializes to disk.

Total wallclock:
- Jacobi reference (all cases): ~90 sec
- AMG target (all cases): ~25 sec (gain ×3.5 on wallclock
  thanks to iter reduction dominating setup cost)

## Acceptance criteria

### Numerical correctness

1. **AMG setup determinism**: same sparse matrix → byte-identical
   preconditioner. Test `v2_amg_setup_determinism`.
2. **Smoother convergence**: Gauss-Seidel alone converges at
   expected rate on Poisson. Test `v2_amg_smoother_convergence`.
3. **V-cycle on Poisson constant**: ≥ 10× CG iter reduction vs
   Jacobi. Test `v2_amg_vcycle_poisson`.
4. **V-cycle on Poisson heterogeneous** (contrast 10⁴): ≥ 5×
   reduction. Test `v2_amg_vcycle_heterogeneous`.
5. **FMG vs V-cycle**: FMG achieves ≥ 2× additional iter reduction
   on Poisson. Test `v2_amg_fmg_poisson`.
6. **Scalar-parity on Step 0 config**: AMG converged state matches
   Jacobi to 1e-5 relative. Test `v2_amg_scalar_parity_step0`.
7. **Scalar-parity on Step 8 config**: same at activated regime.
   Test `v2_amg_scalar_parity_step8`. peak|v| and
   yielding_cell_fraction_max within 1%.

### Regression (bit-parity)

8. **Jacobi default path unchanged**: all Step 0-8 regression
   smoke tests pass byte-identical. Wallclock within [0.95, 1.05]
   of merged reports. No existing numerical result changes.

### Benchmark performance

> **Patch applied after Phase 0 finding** (Phase 0 measured Jacobi
> saturating at the 2000-iter cap on `step8_activated` /
> `step8_activated_128` snapshots, invalidating "×5 reduction from
> 1420" as a commensurable target — see
> [`docs/reports/step8_5a_amg_report.md §Gate revision`](reports/step8_5a_amg_report.md)).
> Items 9-10 are reformulated convergence-first. Item 11 preserved.

9. **Strict convergence required on every benchmark case.** AMG
   must reach CG tolerance `rel_residual ≤ 1e-6` on all nine
   cases, including `step8_activated` and `step8_activated_128`.
   No cap saturation permitted.
10. **Iter-count caps to strict convergence**:
    - `step8_activated` ≤ 400 iters, `step8_activated_128` ≤ 500
      iters (binding).
    - `poisson_constant` ≤ 10, `poisson_contrast_100` ≤ 30,
      `poisson_contrast_10000` ≤ 100, `step0_quiescent` ≤ 20,
      `step3_floor_yielding` ≤ 30, `step6_voronoi` ≤ 40,
      `step7_slab_off` ≤ 40 (binding; measured against strict
      Jacobi convergence reference, not issue heuristics).
11. **Benchmark total wallclock < 30 sec** on AMG path.

### Physics validation (graduation gate, run only when
benchmark targets met)

12. **Re-run Step 0-8 physics with JacobiCG**: byte-identical to
    merged reports. This is the final bit-parity gate. If any
    numerical value differs, the Jacobi path has been accidentally
    modified — remontée.
13. **Re-run Step 0-8 physics with AmgCG**: scalar-parity with
    merged Jacobi reports. peak|v|, mass_conservation_residual,
    yielding_cell_fraction_max all within 1%. Wallclock targets
    (post-Phase-0 revision): AMG wallclock on Step 8 physics
    ≤ 30 % of merged Jacobi wallclock — i.e., ≤ 6 min at 64²
    (vs 19 min Jacobi) and ≤ 33 min at 128² (vs 1h52 Jacobi).
    Steps 0-7 ≥ 3× reduction preserved as before. This wallclock
    gate is the binding product-level commitment; it integrates
    convergence naturally (AMG saturation blows up Newton
    iterations and wallclock explodes).

### Solver health

14. **Newton convergence preserved**: both paths maintain ≥ 95%
    Newton convergence on Step 0-8 configs.
15. **AMG setup phase cost amortized**: setup phase wallclock per
    Newton outer iter < 10% of total solve time. If setup dominates,
    coarsening is too aggressive or levels too many.

### Performance

16. **Wallclock indicative targets** (with AMG, Jacobi-CG reference
    in parentheses):
    - Step 0-5 regression smoke tests: unchanged (Jacobi default)
    - Step 6 physics 64² with AMG: ≤ 10s (vs 34s Jacobi)
    - Step 7 physics 64² with AMG: ≤ 10s (vs 36s Jacobi)
    - Step 8 physics 64² with AMG: ≤ 4 min (vs 19 min Jacobi)
    - Step 8 physics 128² with AMG: ≤ 30 min (vs 1h52 Jacobi)

These are Step 8.5a alone targets. Step 8.5b will compound with
additional ~3× gain (parallelization + Newton extrapolation),
approaching the final cibles (64² at ~40s total physics, 128² at
~5 min total physics).

### Reporting

17. Single report `docs/reports/step8_5a_amg_report.md` covering:
    - Benchmark suite results: Jacobi vs AMG, each case
    - Setup phase profiling (coarsening stats, level count, memory)
    - Scalar-parity diagnostics (max relative error vs Jacobi on
      Step 0-8)
    - Full physics re-run results (Step 0-8 wallclock with AMG)
    - Observations on where AMG helps most vs least
    - Provisional performance note: what remains for Step 8.5b

## Out of scope

- `rayon` parallelization (Step 8.5b)
- Newton extrapolation (Step 8.5b)
- Compilation flags LTO/PGO/target-cpu (Step 8.5b)
- Smoothed Aggregation AMG (fallback if Classical proves insufficient)
- AMG-specific Stokes treatment (block AMG, Vanka smoothers) —
  start with scalar-level Classical and upgrade only if needed
- Promotion of AMG to default (post-milestone task)
- Visual inspection of physics results (Step 9, 10, 10.5)

## Definition of done

- `solver/amg/` module complete with 9 submodules listed
- `solver/sparse_assembly.rs` materializes the momentum operator
- `LinearSolverConfig::{JacobiCG, AmgCG}` implemented
- Benchmark suite `benches/amg_benchmark.rs` with the 9 cases
- 7 new tests for AMG (setup, smoother, V-cycle, FMG, scalar-parity)
- All 9 benchmark cases meet their AMG target
- Benchmark total wallclock < 30 sec with AMG
- Full physics re-run Step 0-8 with JacobiCG: byte-identical to
  merged reports
- Full physics re-run Step 0-8 with AmgCG: scalar-parity within
  1% on peak|v|, yielding, mass conservation
- Wallclock indicative targets met for Step 6-8 with AMG
- `docs/reports/step8_5a_amg_report.md` published
- All Step 0-8 tests still pass in their respective default
  (Jacobi) modes
- PR opened against `milestone/solver-reconstruction`

## Labels

`domain::solver`, `domain::performance`, `type::feature`,
`prio::critical`, `epic::reconstruction`
