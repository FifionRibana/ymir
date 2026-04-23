# Solver reconstruction — implementation roadmap

This document tracks the incremental rebuild of the tectonic solver from
a nondimensional core. It is the implementation companion to
`solver-scaling.md` (physics reference).

## Status

| Step | Description | Status | PR | Report |
|---|---|---|---|---|
| 0 | Nondim Stokes core + S advection (incl. null-space precond) | shipped | #79 | `docs/reports/step0_report.md` |
| 1 | Power-law rheology + Newton solver + continuation | shipped | #81 | `docs/reports/step1_report.md` |
| 2 | GPE spreading (Ar·∇(½S²), staggered flux form) | shipped | #83 | `docs/reports/step2_physics_report.md` (physics) + `docs/reports/step2_regression_report.md` (Step-1 mirror). **#78 remains open** — the GPE gradient spike concerns sharp material interfaces introduced at Steps 5/6, not this step. |
| 3 | Plastic yielding (Von Mises / Bingham, stateless) | shipped | #85 | `docs/reports/step3_physics_report.md` + `docs/reports/step3_regression_report.md`. **Baseline `Bi = 0.15` with GPE-only forcing runs floor-dominated** — `yielding_cell_fraction = 0` at Bi=0.15 is the expected physical outcome (analytic criterion: yielding dominance requires `Bi < ε̇_min^(1/3) = 0.1` when ε̇_II is floor-dominated, cf. report's "Strain-rate regime diagnostic" section). The yielding mechanism itself is validated through the Bi sweep (`yielding_cell_fraction = 1.0` at `Bi ≤ 0.10`), MMS slopes (2.001), Jacobian symmetry (3e-14), and zero-cost regression (ratio 1.00/1.04 vs Step 2). **Checkpoint flagged**: revisit `yielding_cell_fraction` at the Step 4 / 5 / 7 / 8 physics reports — if it remains 0 by Step 7, the source-mechanism coupling to ε̇ is under-dimensioned and warrants remontée. |
| 4 | Basal drag (velocity-damping via operator diagonal) | shipped | #87 | `docs/reports/step4_physics_report.md` + `docs/reports/step4_regression_report.md`. Stateless `f_drag = -Br · S̃² · ṽ` as a positive diagonal augmentation to the Picard block. Preconditioner follows case (B) — diagonal reconstructed analytically in `stokes/operator.rs::momentum_diagonal`, consistency with `apply_momentum`'s stencil guarded by `tests/v2_precond_drag_diagonal`. **Zero-cost when disabled** confirmed by regression ratio ≈ 1.00 vs Step 3 (wallclock and CG iters). **Baseline drag/visc ratio ≈ 10⁻⁷** (corrected vs prompt: η ≈ 100 at the floor-dominated power-law baseline, not η ≈ 1), so peak \|v\| damping is below 4-digit precision; the Br sweep's strict `peak \|v\|` monotonicity (at f64) is the physical signal. `yielding_cell_fraction = 0` at baseline (yielding Disabled here; full checkpoint deferred to Step 5/7 per roadmap). |
| 5 | Boundary sources/sinks | shipped | #89 | `docs/reports/step5_physics_report.md` + `docs/reports/step5_regression_report.md` + `docs/reports/step5_reference_variant_report.md`. Five source/sink terms (Q_sub, Q_arc, Q_spread, Q_coll-v, Q_rift-v) on a prescribed static layout (`horizontal_oceanic_strip`). `BoundaryConfig::Disabled` → structural bypass (advection-only, no Q, no clamp, no tracking). Regression convention installed: step N activates all mechanisms through N-1 in their canonical configuration. Because Step 4 physics was non-canonical (yielding `Disabled` for Br isolation), Step 5 produces a "reference variant" on this branch (Step 4 config with yielding Enabled) as the regression's comparison target. **Baseline s_oceanic_mean = 0.2158 (64²) / 0.2085 (128²)** both in `[0.18, 0.22]` post-calibration (calibrated `k_spread = 0.050`, first-probe convergence at the bracket's low end — at Step 5 baseline's tiny convergent-velocity regime, Q_sub barely fires and any sizable `k_spread` grows the oceanic strip; the balance appears at Steps 7/8). **mass_balance_residual < 10⁻¹²** at both grids (spec acceptance `< 1%`). **Clamp never fires** (clamp_activation_fraction = 0). **k_sub sweep monotonicity strict at f64 precision**. **yielding_cell_fraction = 0** — yielding checkpoint still deferred to Step 7 as expected; the Step 5 baseline remains floor-dominated. **CG ratio vs Step 4 physics ≈ 2.2×** — the heterogeneous S̃² (oceanic 0.04 vs continental 1.0 adjacent) stresses the Jacobi preconditioner; above the spec's 1.3× target, but Newton still converges 100% at both grids. **boundary_type_diversity = 2** (subduction + spread active). **#78 monitoring**: peak f_GPE on oceanic/continental interfaces ≈ 3.58 at 64², recorded as trajectory telemetry through Steps 5–8. |
| 6 | Conservative recycling | blocked by 5 | — | — |
| 7 | Slab-pull regularized RHS | blocked by 6 | — | — |
| 8 | Mantle flow | blocked by 7 | — | — |
| 9 | Cratonic plastic immunity | blocked by 8 | — | — |
| 10 | Geological age field | blocked by 9 | — | — |

| Transverse | Description | Status |
|---|---|---|
| ~~T1~~ | ~~Null-space-aware preconditioner~~ | **absorbed into Step 0** — shipped as part of `stokes/precond.rs` + `stokes/nullspace.rs` with verification in `tests/v2_nullspace.rs` |
| T2 | Diagnostics framework extension | MVP introduced step 0; extended each subsequent step |
| T3 | Stochastic validation harness | not started (introduced after step 4) |

## Conventions

- Each step lives on its own commit on `reconstruction/solver-from-scratch`
- Each step has exactly one PR against the reconstruction branch
- Each step produces a report in `logs/stepN_report.md` with metrics
  delta vs step N-1
- The `tectonics_v2/` module tree is built incrementally; no code from
  `tectonics/` is imported except for utility types (Field2D,
  PeriodicIndex)

## Reusable infrastructure

The reporting framework introduced in step 0 is the contract. Each
step extends it with new metrics but does not modify the core
comparison logic.

Reference scenario for cross-step comparison:
- Grid: 64² and 128²
- Seed: 42
- Duration: T = 6·τ*
- Number of macro steps: 300
- Active extensions: those introduced up to and including the current
  step

## Final integration

Once step 10 lands, T1 and T3 must be done before merging the
reconstruction branch back to main. T2 is built incrementally
throughout.

The final merge is a single PR with full diff review. The new tree
replaces the existing `tectonics/` module wholesale; UI integration in
ymir-viz is updated in the same merge.