# Step 1 — Power-law rheology + Newton solver (baseline)

> **Step 1 reference run for milestone "Solver reconstruction".**
> Compared against Step 0 (`docs/reports/step0_report.md`).
> Subsequent steps' reports will diff against this one.

- Seed: `42`
- Formulation: thin viscous sheet + power-law rheology (n > 1). Linear internal system per Newton iteration is symmetric (variational structure); CG suffices.

## Physical scales

```
Scales:
- L* = 3.500e5 m (350 km)
- S* = 3.500e4 m (35 km)
- τ* = 9.467e14 s (30.00 Myr)
- ρ* = 3300.0 kg/m³
- v* = 3.697e-10 m/s (1.167 cm/yr)
- ε̇* = 1.056e-15 1/s
- η* = 1.073e24 Pa·s
- σ* = 1.133e9 Pa
- p* = 1.133e9 Pa
- f* = 3.237e3 N/m³
- Ar = 0.100
- 2π check (not used, informational): 6.283185
```

## Discretisation validation (MMS convergence at report time)

The baseline run at Step 1 does not fully exercise power-law behaviour (the placeholder forcing saturates against `ε̇_min`). The following MMS convergence checks are run **every time the report is generated** so the discretisation remains visibly verified.

### Constant η (Picard path, Step 0 operator)

Manufactured solution `v = (sin(2πx), sin(2πy))`, `η = 1`.

| N | v_err RMS | slope to next |
|---|---|---|
| 16 | 9.158e-3 | 2.008 |
| 32 | 2.276e-3 | 2.002 |
| 64 | 5.682e-4 | 2.001 |
| 128 | 1.420e-4 | — |

Final slope: `2.001` (expected ≥ 1.7; quadratic target = 2.0).

### Variable η (linear, prescribed η field)

Manufactured solution `v = (sin(2πx)cos(2πy), -cos(2πx)sin(2πy))`, `η(x,y) = 1 + 0.5·sin(2πx)·cos(2πy)`. Validates the η-variable Picard path used under Newton.

| N | v_err RMS | slope to next |
|---|---|---|
| 32 | 1.545e-3 | 2.002 |
| 64 | 3.857e-4 | 2.001 |
| 128 | 9.640e-5 | — |

Final slope: `2.001` (expected ≥ 1.7).

### Nonlinear Newton tail (n = 3)

Target-generated RHS on 32² grid. Newton outer iterations: `8`.

Residual trail:

```
  8.695e2
  6.847e2
  4.629e2
  2.082e2
  4.584e1
  2.598e0
  1.338e-2
  4.310e-5
  1.759e-7
```

Tail reductions: `310.4×` then `244.9×` (super-linear target: both ≥ 100×; strict quadratic requires an exact inner solve).

## Grid 64×64

### Solver configuration

| field | value |
|---|---|
| formulation | thin viscous sheet (elliptic, no pressure) with power-law rheology |
| discretization | MAC staggered (v face / η S cell-centre / ε̇_xy corner) |
| η averaging to corners | arithmetic 4-point at corners (see operator.rs) |
| preconditioner | velocity Jacobi (Picard-block diagonal), null-space wrapped |
| gauge fixing | mean(vx), mean(vy) projected before & after every M⁻¹ + post-solve |
| preset | `dynamic-accidented` |
| nonlinear solver | `newton` |
| rheology `n` (after continuation) | 3.00 |
| rheology `ε̇_min` | 1.0e-3 |
| rheology `η_max` (soft cap) | 1.0e3 |
| continuation schedule | `[1.0, 1.5, 2.0, 2.5, 3.0]` |
| Newton rel tol | 1.0e-6 |
| Newton max outer iters | 20 |
| CG tolerance | 1.0e-8 |
| CG max iter | 2000 |
| CFL factor | 0.30 |
| grid spacing (nondim) | 0.015625 |
| body force | SinusoidalForce(ε=10) |
| seed | 42 |

### Timing

- wallclock total: `0.501 s`
- wallclock per step (mean): `1.669 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `5.79e0`
- CG iterations per Newton step — mean: `22.8`, max: `38`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 8 | 2 |
  | 15 | 1 |
  | 23 | 18 |
  | 30 | 8 |
  | 38 | 9 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `0.1`, max: `9`
- effective η_max/η_min over run — mean: `26.37`, max: `26.37`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%` (spec target < 1%)
- continuation ramp: ✅ all 5 sub-solves converged

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `1.665e-15`

### Null-space health

- max |mean(vx)| across solves: `2.435e-20`
- max |mean(vy)|: `0.000e0`

### Velocity magnitude

- peak |v|: `2.738e-2`

### Heightmaps of S

- `docs/reports/step1_heightmaps/s_64x64_t0000.png`
- `docs/reports/step1_heightmaps/s_64x64_t0150.png`
- `docs/reports/step1_heightmaps/s_64x64_t0300.png`

### Comparison with Step 0

#### Grid 64×64 — comparison vs Step 0

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 0.152 | 0.501 | ×3.29 |
| CG iters / linear solve (mean) | 0.0 (solver-trivial) | 22.8 | N/A — no denominator; report absolute |
| S mass drift (relative) | 2.331e-15 | 1.665e-15 | gate 1e-10 |
| max \|mean(vx)\| | 2.647e-23 | 2.435e-20 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 0.000e0 | bruit machine |

### Dormant metrics (inactive at Step 1)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| age field stats | Step 10 |

## Grid 128×128

### Solver configuration

| field | value |
|---|---|
| formulation | thin viscous sheet (elliptic, no pressure) with power-law rheology |
| discretization | MAC staggered (v face / η S cell-centre / ε̇_xy corner) |
| η averaging to corners | arithmetic 4-point at corners (see operator.rs) |
| preconditioner | velocity Jacobi (Picard-block diagonal), null-space wrapped |
| gauge fixing | mean(vx), mean(vy) projected before & after every M⁻¹ + post-solve |
| preset | `dynamic-accidented` |
| nonlinear solver | `newton` |
| rheology `n` (after continuation) | 3.00 |
| rheology `ε̇_min` | 1.0e-3 |
| rheology `η_max` (soft cap) | 1.0e3 |
| continuation schedule | `[1.0, 1.5, 2.0, 2.5, 3.0]` |
| Newton rel tol | 1.0e-6 |
| Newton max outer iters | 20 |
| CG tolerance | 1.0e-8 |
| CG max iter | 2000 |
| CFL factor | 0.30 |
| grid spacing (nondim) | 0.007812 |
| body force | SinusoidalForce(ε=10) |
| seed | 42 |

### Timing

- wallclock total: `2.442 s`
- wallclock per step (mean): `8.140 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `2.52e1`
- CG iterations per Newton step — mean: `48.2`, max: `77`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 16 | 1 |
  | 31 | 0 |
  | 46 | 19 |
  | 61 | 11 |
  | 77 | 9 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `0.1`, max: `9`
- effective η_max/η_min over run — mean: `29.13`, max: `29.13`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%` (spec target < 1%)
- continuation ramp: ✅ all 5 sub-solves converged

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `-1.998e-15`

### Null-space health

- max |mean(vx)| across solves: `1.784e-20`
- max |mean(vy)|: `0.000e0`

### Velocity magnitude

- peak |v|: `2.736e-2`

### Heightmaps of S

- `docs/reports/step1_heightmaps/s_128x128_t0000.png`
- `docs/reports/step1_heightmaps/s_128x128_t0150.png`
- `docs/reports/step1_heightmaps/s_128x128_t0300.png`

### Comparison with Step 0

#### Grid 128×128 — comparison vs Step 0

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 0.527 | 2.442 | ×4.63 |
| CG iters / linear solve (mean) | 0.0 (solver-trivial) | 48.2 | N/A — no denominator; report absolute |
| S mass drift (relative) | -3.553e-15 | -1.998e-15 | gate 1e-10 |
| max \|mean(vx)\| | 1.423e-22 | 1.784e-20 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 0.000e0 | bruit machine |

### Dormant metrics (inactive at Step 1)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| age field stats | Step 10 |

---
*Generated by `cargo run --release --bin step_baseline`.*
