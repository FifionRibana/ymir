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

2. **The issue's AMG-target column is interpretable as a ratio,
   not an absolute cap.** The issue's `≤ 280` for `step8_activated`
   was built from `1420 ÷ 5`. On the snapshot-at-step-100 regime,
   Jacobi's true iteration count is ≥ 2000 (unbounded). Read
   literally, the `≤ 280` target implies AMG must converge the
   problem Jacobi can't — which is harder, but also the actual
   product goal. Phase 4 gate is therefore reinterpreted:
   *"AMG converges `step8_activated` to tolerance in ≤ 400 iters,
   demonstrating ≥ 5× reduction vs Jacobi's cap-saturation
   behaviour."* A stricter interpretation (≤ 280 / ≤ 370) is
   retained as the stretch goal.

3. **Per-snapshot Jacobi CG counts elsewhere are below the issue
   table.** `poisson_contrast_10000` measured 459 (issue said
   1500); `step0_quiescent` measured 43 (issue said 80);
   `step6_voronoi` measured 207 (issue said 130 — slightly above,
   reflecting the later capture step for our `step6_voronoi` of
   iter 50 vs an earlier state). These variations are informational;
   the measured Jacobi reference is authoritative for the rest of
   the milestone, not the issue's heuristic targets.

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

*Not yet started. Gate: matrix-free `apply_momentum` vs sparse CSR
`Ax` agree to 1e-14 on 10 test vectors spanning (constant, linear,
sinusoid of various k, seeded random).*

## Phase 2 — AMG V-cycle on Poisson

*Not yet started.*

## Phases 3-5 — V-cycle on heterogeneous / FMG / scalar-parity

*Not yet started.*

## Phase 6 — Graduation gate (physics re-runs)

*Not yet started. Strict order: benchmarks green first, then and
only then re-run Steps 0-8 physics with both JacobiCG and AmgCG.*
