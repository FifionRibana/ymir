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

## Ar sweep (diagnostic)

Honest thin-sheet value from the default scales: **`Ar = S*/L* = 0.100`** (used in the baseline above).

The design note's historical target `Ar ∈ [1, 5]` is mathematically incompatible with `S* ≪ L*` — the thin-sheet assumption forces `Ar ≪ 1`. The sweep below tabulates the GPE-only response at 64²·300 steps for a range that brackets both the honest value and the historical band, so the discretisation's behaviour across `Ar` is visible quantitatively.

| Ar | Var(S̃) init | Var(S̃) final | ratio | peak \|∇S̃\| | peak \|v\| | Newton conv | CG mean | mass drift | wallclock (s) |
|---|---|---|---|---|---|---|---|---|---|
| `0.10` | `1.000e-2` | `9.938e-3` | `0.994` | `1.255e0` | `8.155e-6` | `100%` | `51.5` | `2.00e-15` | `5.692` |
| `0.50` | `1.000e-2` | `9.660e-3` | `0.966` | `1.255e0` | `4.482e-5` | `100%` | `58.4` | `2.22e-16` | `14.851` |
| `1.00` | `1.000e-2` | `9.238e-3` | `0.924` | `1.255e0` | `1.016e-4` | `100%` | `62.2` | `-3.33e-16` | `22.096` |
| `2.00` | `1.000e-2` | `8.119e-3` | `0.812` | `1.255e0` | `2.674e-4` | `100%` | `68.1` | `-2.22e-16` | `34.289` |
| `5.00` | `1.000e-2` | `3.799e-3` | `0.380` | `1.255e0` | `1.631e-3` | `100%` | `83.1` | `3.55e-15` | `35.851` |

**Interpretation** — GPE dissipation scales as `Ar/τ*`, so the characteristic spreading time is `τ*/Ar`. At `Ar = 0.1` this is `~10·τ*`, ten times the tectonic time scale and well beyond the 300-step run (`6·τ*`). The variance ratio across the sweep confirms the expected monotonic response: lower `Ar` → slower spreading → larger `Var(S̃)_final / Var(S̃)_initial`. Narrative-level dynamics (continents building and breaking on the run window) must therefore come from the mechanisms being added at Steps 3–10, not from GPE alone — see `solver-scaling.md` §5.1.bis for the characteristic-time ordering.

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

- wallclock total: `5.439 s`
- wallclock per step (mean): `18.129 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `2.85e1`
- CG iterations per Newton step — mean: `51.5`, max: `109`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 25 | 1 |
  | 46 | 5 |
  | 67 | 597 |
  | 88 | 9 |
  | 109 | 8 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `2.0`, max: `6`
- effective η_max/η_min over run — mean: `1.03`, max: `1.03`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### S field evolution

- Var(S̃) timeline: initial `1.000e-2`, middle `9.969e-3`, final `9.938e-3` (Δ = `-0.62%` vs initial)
- max|∇S̃| timeline: initial `1.255e0`, peak `1.255e0`, final `1.251e0`

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `1.998e-15`

### Null-space health

- max |mean(vx)| across solves: `7.858e-23`
- max |mean(vy)|: `1.380e-22`

### Velocity magnitude

- peak |v|: `8.155e-6`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step2_physics_heightmaps/s_64x64_t0000.png` | `8.005e-1` | `1.200e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_64x64_t0000_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_64x64_t0150.png` | `8.007e-1` | `1.199e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_64x64_t0150_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_64x64_t0300.png` | `8.009e-1` | `1.199e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_64x64_t0300_colorbar.png` |

### Comparison vs Step 1 (advisory — physics changed, not a regression test)

#### Grid 64×64 — comparison vs Step 1

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 0.501 | 5.439 | ×10.86 |
| CG iters / linear solve (mean) | 22.8 | 51.5 | ×2.26 [acceptable] |
| S mass drift (relative) | 1.665e-15 | 1.998e-15 | gate 1e-10 |
| max \|mean(vx)\| | 2.435e-20 | 7.858e-23 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 1.380e-22 | bruit machine |

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

- wallclock total: `49.381 s`
- wallclock per step (mean): `164.602 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `1.50e2`
- CG iterations per Newton step — mean: `117.2`, max: `222`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 50 | 1 |
  | 93 | 0 |
  | 136 | 601 |
  | 179 | 9 |
  | 222 | 9 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `2.0`, max: `6`
- effective η_max/η_min over run — mean: `1.03`, max: `1.03`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### S field evolution

- Var(S̃) timeline: initial `1.000e-2`, middle `9.969e-3`, final `9.938e-3` (Δ = `-0.62%` vs initial)
- max|∇S̃| timeline: initial `1.256e0`, peak `1.256e0`, final `1.252e0`

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `4.663e-15`

### Null-space health

- max |mean(vx)| across solves: `8.261e-23`
- max |mean(vy)|: `2.057e-22`

### Velocity magnitude

- peak |v|: `8.160e-6`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step2_physics_heightmaps/s_128x128_t0000.png` | `8.001e-1` | `1.200e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_128x128_t0000_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_128x128_t0150.png` | `8.003e-1` | `1.199e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_128x128_t0150_colorbar.png` |
| `docs/reports/step2_physics_heightmaps/s_128x128_t0300.png` | `8.006e-1` | `1.199e0` | `1.000e0` | `docs/reports/step2_physics_heightmaps/s_128x128_t0300_colorbar.png` |

### Comparison vs Step 1 (advisory — physics changed, not a regression test)

#### Grid 128×128 — comparison vs Step 1

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 2.442 | 49.381 | ×20.22 (>20×, flag) |
| CG iters / linear solve (mean) | 48.2 | 117.2 | ×2.43 [acceptable] |
| S mass drift (relative) | -1.998e-15 | 4.663e-15 | gate 1e-10 |
| max \|mean(vx)\| | 1.784e-20 | 8.261e-23 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 2.057e-22 | bruit machine |

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
