# R5b D1 — amg_cg

32² × 10 steps × mf=1.0 × workflow OFF, seed=42.

## Solver health (aggregate)

- Linear solver: `amg_cg`
- Total runtime: **47.0s** (4.70s/step)
- CG iter mean: **302.2**  (max: 2000)
- Kappa estimate: 9.99e2
- Peak |v|: 2.993e-3
- Outcomes (Converged / Stalled / Diverged / Capped): **15 / 0 / 0 / 0**

## Per-step Newton outer iterations

| step | Newton outer iter |
|---|---|
| 1 | 1 |
| 2 | 6 |
| 3 | 5 |
| 4 | 5 |
| 5 | 4 |
| 6 | 10 |
| 7 | 7 |
| 8 | 4 |
| 9 | 4 |
| 10 | 4 |
| 11 | 4 |
| 12 | 4 |
| 13 | 4 |
| 14 | 4 |
| 15 | 4 |

## CG iterations per Newton inner solve

- Total Newton inner solves: **70**
- CG iter min/mean/max: **3 / 302.2 / 2000**

First 20 inner-CG iter counts: `[3, 52, 52, 48, 42, 38, 32, 45, 47, 41, 35, 27, 43, 46, 40, 34, 22, 41, 45, 38]`

## Histogram (CG iter per inner solve)

- bin edges (≤): `[402, 801, 1201, 1600, 2000]`
- counts: `[60, 0, 0, 0, 10]`

## Cost breakdown (AMG-specific, not instrumented)

AMG setup cost (hierarchy build) and per-V-cycle time are **not** exposed by `Metrics` in this iteration. Inferring from total wallclock − Newton iter count × estimated CG iter time is approximate. If decisive, instrument in D2.

## Regime

**CONVERGENT** (CG mean ≤ 500)
