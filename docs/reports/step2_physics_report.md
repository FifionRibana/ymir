# Step 2 — GPE spreading (physics)

> **Step 2 physics run for milestone "Solver reconstruction".**
> This run uses `GpeForce` — the first **physical** term in the milestone. The placeholder sinusoidal force is retained for the companion regression report.
> Compared against Step 1 only on physical quantities (peak |v|, S range, variance, max |∇S|); numerical solver regression lives in the companion regression report.

- Seed: `42`
- Ar (Argand) = `0.100` — **derived** from the 4 primary scales; never a direct knob. See `scales::Scales::argand_number` for the `solver-scaling.md` §5.1 range inconsistency note.

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

### GPE force (staggered `-Ar·∇(½S²)`, smooth S)

Smooth manufactured `S = 1 + 0.1·sin(2πx)·cos(2πy)` at `Ar = 2`. Validates the GPE discretisation introduced at Step 2 against the analytic `-Ar·S·∇S`.

| N | v_err RMS | slope to next |
|---|---|---|
| 32 | 1.024e-3 | 1.999 |
| 64 | 2.561e-4 | 2.000 |
| 128 | 6.402e-5 | — |

Final slope: `2.000` (expected ≥ 1.7).

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
| body force | ForceSum [gpe]: GpeForce (Ar = 0.100 from scales) |
| seed | 42 |

### Timing

- wallclock total: `2.062 s`
- wallclock per step (mean): `6.874 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `1.19e1`
- CG iterations per Newton step — mean: `33.5`, max: `90`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 22 | 1 |
  | 39 | 300 |
  | 56 | 2 |
  | 73 | 1 |
  | 90 | 9 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `1.0`, max: `4`
- effective η_max/η_min over run — mean: `1.00`, max: `1.00`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### S field evolution

- Var(S̃) timeline: initial `1.000e-4`, middle `9.970e-5`, final `9.940e-5` (Δ = `-0.60%` vs initial)
- max|∇S̃| timeline: initial `1.255e-1`, peak `1.255e-1`, final `1.251e-1`

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `-4.885e-15`

### Null-space health

- max |mean(vx)| across solves: `7.005e-24`
- max |mean(vy)|: `1.004e-23`

### Velocity magnitude

- peak |v|: `7.969e-7`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step2_physics_heightmaps/s_64x64_t0000.png` | `9.800e-1` | `1.020e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_64x64_t0000_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_64x64_t0150.png` | `9.801e-1` | `1.020e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_64x64_t0150_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_64x64_t0300.png` | `9.801e-1` | `1.020e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_64x64_t0300_colorbar.png` |

### Comparison vs Step 1 (advisory — physics changed, not a regression test)

#### Grid 64×64 — comparison vs Step 1

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 0.501 | 2.062 | ×4.12 |
| CG iters / linear solve (mean) | 22.8 | 33.5 | ×1.47 [acceptable] |
| S mass drift (relative) | 1.665e-15 | -4.885e-15 | gate 1e-10 |
| max \|mean(vx)\| | 2.435e-20 | 7.005e-24 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 1.004e-23 | bruit machine |

### Dormant metrics (inactive at Step 2)

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
| body force | ForceSum [gpe]: GpeForce (Ar = 0.100 from scales) |
| seed | 42 |

### Timing

- wallclock total: `14.754 s`
- wallclock per step (mean): `49.179 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `5.21e1`
- CG iterations per Newton step — mean: `69.5`, max: `185`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 45 | 1 |
  | 80 | 299 |
  | 115 | 3 |
  | 150 | 1 |
  | 185 | 9 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `1.0`, max: `4`
- effective η_max/η_min over run — mean: `1.00`, max: `1.00`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### S field evolution

- Var(S̃) timeline: initial `1.000e-4`, middle `9.970e-5`, final `9.940e-5` (Δ = `-0.60%` vs initial)
- max|∇S̃| timeline: initial `1.256e-1`, peak `1.256e-1`, final `1.252e-1`

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `7.772e-15`

### Null-space health

- max |mean(vx)| across solves: `6.294e-24`
- max |mean(vy)|: `9.819e-24`

### Velocity magnitude

- peak |v|: `7.974e-7`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step2_physics_heightmaps/s_128x128_t0000.png` | `9.800e-1` | `1.020e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_128x128_t0000_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_128x128_t0150.png` | `9.800e-1` | `1.020e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_128x128_t0150_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_128x128_t0300.png` | `9.801e-1` | `1.020e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_128x128_t0300_colorbar.png` |

### Comparison vs Step 1 (advisory — physics changed, not a regression test)

#### Grid 128×128 — comparison vs Step 1

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 2.442 | 14.754 | ×6.04 |
| CG iters / linear solve (mean) | 48.2 | 69.5 | ×1.44 [acceptable] |
| S mass drift (relative) | -1.998e-15 | 7.772e-15 | gate 1e-10 |
| max \|mean(vx)\| | 1.784e-20 | 6.294e-24 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 9.819e-24 | bruit machine |

### Dormant metrics (inactive at Step 2)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| age field stats | Step 10 |

---
*Generated by `cargo run --release --bin step_baseline`.*
