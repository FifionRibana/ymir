# Solver reconstruction — implementation roadmap

This document tracks the incremental rebuild of the tectonic solver from
a nondimensional core. It is the implementation companion to
`solver-scaling.md` (physics reference).

## Status

| Step | Description | Status | PR | Report |
|---|---|---|---|---|
| 0 | Nondim Stokes core + S advection (incl. null-space precond) | shipped | #79 | `docs/reports/step0_report.md` |
| 1 | Power-law rheology | blocked by 0 | — | — |
| 2 | GPE spreading | blocked by 1 | — | — |
| 3 | Plastic yielding | blocked by 2 | — | — |
| 4 | Basal drag | blocked by 3 | — | — |
| 5 | Boundary sources/sinks | blocked by 4 | — | — |
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