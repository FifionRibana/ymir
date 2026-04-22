# Step 3 — GPE spreading + plastic yielding (physics)

> **Step 3 physics run for milestone "Solver reconstruction".**
> `GpeForce` (Ar = 0.1) + `YieldingConfig::Enabled` with `Bi = 0.15`. Von Mises / Bingham, stateless — no plastic memory, no healing, no cratonic immunity. Power-law + smooth cap unchanged since Step 1; only the effective-viscosity blend and the Jacobian chain rule differ.
> Solver unchanged: CG (the tangent Jacobian remains symmetric under arithmetic corner averaging, whether or not yielding is active — see `stokes/operator.rs` doc-comment).

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

## Bi sweep (diagnostic)

Baseline `Bi = 0.15` (preset `dynamic-accidented`, design note §5.1 centre of range). The sweep below covers `Bi ∈ {0.05, 0.10, 0.15, 0.30, 0.50}` at 64²·N steps with `GpeForce` + yielding Enabled. Expected qualitative behaviour: yielding cells widespread at low Bi (plasticity takes over early), rare at high Bi (plastic branch rarely wins the soft-min). Strict monotonic decrease of `yielding_cell_fraction` with Bi is the acceptance invariant (issue #85).

| Bi | yielding_cell_fraction | yielding_intensity | S min | S max | S mean | Newton conv | Newton iter mean | CG mean | peak \|v\| | mass drift | wallclock (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `0.05` | `1.000` | `3.206` | `NaN` | `NaN` | `NaN` | `100%` | `5.0` | `53.2` | `3.701e-5` | `-1.11e-15` | `14.66` |
| `0.10` | `1.000` | `1.077` | `NaN` | `NaN` | `NaN` | `100%` | `3.0` | `55.0` | `1.737e-5` | `8.88e-16` | `9.31` |
| `0.15` | `0.000` | `0.446` | `NaN` | `NaN` | `NaN` | `100%` | `3.0` | `50.7` | `1.192e-5` | `-4.44e-16` | `8.73` |
| `0.30` | `0.000` | `0.000` | `NaN` | `NaN` | `NaN` | `100%` | `2.0` | `52.1` | `8.558e-6` | `-3.77e-15` | `6.14` |
| `0.50` | `0.000` | `0.000` | `NaN` | `NaN` | `NaN` | `100%` | `2.0` | `51.5` | `8.211e-6` | `2.22e-15` | `6.09` |

**Monotonicity of `yielding_cell_fraction` vs Bi** : ✅ strictly non-increasing across the 5 points, as required.

**Interpretation** — low-Bi yielding is pervasive (every cell tastes the plastic branch once its strain is above the floor), high-Bi yielding is confined to the occasional active zone. Newton iteration budget grows modestly at low Bi (more cells to linearise through the blend), in line with expectations.

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

- wallclock total: `8.960 s`
- wallclock per step (mean): `29.868 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `2.85e1`
- CG iterations per Newton step — mean: `50.7`, max: `109`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 25 | 1 |
  | 46 | 300 |
  | 67 | 602 |
  | 88 | 8 |
  | 109 | 10 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `3.0`, max: `6`
- effective η_max/η_min over run — mean: `1.07`, max: `1.07`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Plastic yielding

- Bi = `0.150`
- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `0.000`
- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `0.446`
- Definition notes: the `< 0.5·η_visc` criterion captures "yielding dominant", not "yielding present anywhere" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).

### Strain-rate regime diagnostic (floor-domination)

- `ε̇_min` (regularisation floor): `1.000e-3`
- `ε̇_II` at final timestep: mean = `5.078e-5`, max = `8.087e-5`
- Fraction of cells with `ε̇_II < 10·ε̇_min = 1.000e-2` at final timestep: `1.000` (1.0 = everywhere in the floor-dominated band)
- `max(ε̇_II) / ε̇_min` = `0.08` — ratio of the strongest strain-rate cell to the regularisation floor.

**Verdict: floor-dominated.** `ε̇_II` lies below `ε̇_min` over most of the domain; the viscous and plastic branches both saturate at their floor values. The analytic criterion for yielding dominance in this regime is

```
  Bi < ε̇_min^(1/3)   (n = 3)
  ⟺ Bi < 0.100
```

with the default scales. The baseline `Bi = 0.150` sits above this threshold, so `yielding_cell_fraction = 0` is the **expected** diagnostic outcome for Ar = 0.1 + GPE-only forcing. The Bi sweep at `Bi ≤ 0.10` crosses the threshold and shows `yielding_cell_fraction = 1.0`, confirming the yielding mechanism is wired correctly — it simply is not activated by the weak GPE regime at this baseline.

**Anticipated cross-over at later steps.** Mechanisms introduced in Steps 4 (basal drag), 5 (boundary sources), 7 (slab pull), and 8 (mantle flow) inject energy at faster time scales than GPE and should push `ε̇_II` into the O(1) range in active zones. As soon as active-zone `ε̇_II > ε̇_min`, the `Bi = 0.15` criterion for yielding dominance will hold locally and `yielding_cell_fraction > 0` should appear naturally. If `yielding_cell_fraction` is still 0 after Step 7 (slab pull, which acts at τ*/Sp ≈ 10–60 Myr), the coupling between source mechanisms and ε̇ is under-dimensioned and warrants a remontée — **this is flagged as a checkpoint** for the Step 4, 5, 7, 8 physics reports.

Basal drag (Step 4) is dissipative and may *not* raise `ε̇_II`; the threshold check starts in earnest at Step 5 (boundary sources inject mass and create strain) and carries through Steps 7–8.

### S field evolution

- Var(S̃) timeline: initial `1.000e-2`, middle `9.954e-3`, final `9.909e-3` (Δ = `-0.91%` vs initial)
- max|∇S̃| timeline: initial `1.255e0`, peak `1.255e0`, final `1.249e0`

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `-4.441e-16`

### Null-space health

- max |mean(vx)| across solves: `1.067e-22`
- max |mean(vy)|: `1.733e-22`

### Velocity magnitude

- peak |v|: `1.192e-5`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step3_physics_heightmaps/s_64x64_t0000.png` | `8.005e-1` | `1.200e0` | `1.000e0` | `docs/reports/step3_physics_heightmaps/s_64x64_t0000_colorbar.png` |
| `docs/reports/step3_physics_heightmaps/s_64x64_t0150.png` | `8.008e-1` | `1.199e0` | `1.000e0` | `docs/reports/step3_physics_heightmaps/s_64x64_t0150_colorbar.png` |
| `docs/reports/step3_physics_heightmaps/s_64x64_t0300.png` | `8.012e-1` | `1.198e0` | `1.000e0` | `docs/reports/step3_physics_heightmaps/s_64x64_t0300_colorbar.png` |

### Comparison vs Step 2 (advisory — yielding added, not a regression test)

#### Grid 64×64 — comparison vs Step 1

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 5.439 | 8.960 | ×1.65 |
| CG iters / linear solve (mean) | 51.5 | 50.7 | ×0.99 [idéal] |
| S mass drift (relative) | 1.998e-15 | -4.441e-16 | gate 1e-10 |
| max \|mean(vx)\| | 7.858e-23 | 1.067e-22 | bruit machine |
| max \|mean(vy)\| | 1.380e-22 | 1.733e-22 | bruit machine |

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

- wallclock total: `69.898 s`
- wallclock per step (mean): `232.995 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `1.35e2`
- CG iterations per Newton step — mean: `111.4`, max: `222`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 50 | 1 |
  | 93 | 300 |
  | 136 | 601 |
  | 179 | 8 |
  | 222 | 11 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `3.0`, max: `6`
- effective η_max/η_min over run — mean: `1.07`, max: `1.07`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Plastic yielding

- Bi = `0.150`
- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `0.000`
- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `0.446`
- Definition notes: the `< 0.5·η_visc` criterion captures "yielding dominant", not "yielding present anywhere" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).

### Strain-rate regime diagnostic (floor-domination)

- `ε̇_min` (regularisation floor): `1.000e-3`
- `ε̇_II` at final timestep: mean = `5.075e-5`, max = `8.104e-5`
- Fraction of cells with `ε̇_II < 10·ε̇_min = 1.000e-2` at final timestep: `1.000` (1.0 = everywhere in the floor-dominated band)
- `max(ε̇_II) / ε̇_min` = `0.08` — ratio of the strongest strain-rate cell to the regularisation floor.

**Verdict: floor-dominated.** `ε̇_II` lies below `ε̇_min` over most of the domain; the viscous and plastic branches both saturate at their floor values. The analytic criterion for yielding dominance in this regime is

```
  Bi < ε̇_min^(1/3)   (n = 3)
  ⟺ Bi < 0.100
```

with the default scales. The baseline `Bi = 0.150` sits above this threshold, so `yielding_cell_fraction = 0` is the **expected** diagnostic outcome for Ar = 0.1 + GPE-only forcing. The Bi sweep at `Bi ≤ 0.10` crosses the threshold and shows `yielding_cell_fraction = 1.0`, confirming the yielding mechanism is wired correctly — it simply is not activated by the weak GPE regime at this baseline.

**Anticipated cross-over at later steps.** Mechanisms introduced in Steps 4 (basal drag), 5 (boundary sources), 7 (slab pull), and 8 (mantle flow) inject energy at faster time scales than GPE and should push `ε̇_II` into the O(1) range in active zones. As soon as active-zone `ε̇_II > ε̇_min`, the `Bi = 0.15` criterion for yielding dominance will hold locally and `yielding_cell_fraction > 0` should appear naturally. If `yielding_cell_fraction` is still 0 after Step 7 (slab pull, which acts at τ*/Sp ≈ 10–60 Myr), the coupling between source mechanisms and ε̇ is under-dimensioned and warrants a remontée — **this is flagged as a checkpoint** for the Step 4, 5, 7, 8 physics reports.

Basal drag (Step 4) is dissipative and may *not* raise `ε̇_II`; the threshold check starts in earnest at Step 5 (boundary sources inject mass and create strain) and carries through Steps 7–8.

### S field evolution

- Var(S̃) timeline: initial `1.000e-2`, middle `9.955e-3`, final `9.910e-3` (Δ = `-0.90%` vs initial)
- max|∇S̃| timeline: initial `1.256e0`, peak `1.256e0`, final `1.251e0`

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `-2.220e-15`

### Null-space health

- max |mean(vx)| across solves: `8.758e-23`
- max |mean(vy)|: `2.708e-22`

### Velocity magnitude

- peak |v|: `1.193e-5`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step3_physics_heightmaps/s_128x128_t0000.png` | `8.001e-1` | `1.200e0` | `1.000e0` | `docs/reports/step3_physics_heightmaps/s_128x128_t0000_colorbar.png` |
| `docs/reports/step3_physics_heightmaps/s_128x128_t0150.png` | `8.005e-1` | `1.199e0` | `1.000e0` | `docs/reports/step3_physics_heightmaps/s_128x128_t0150_colorbar.png` |
| `docs/reports/step3_physics_heightmaps/s_128x128_t0300.png` | `8.008e-1` | `1.199e0` | `1.000e0` | `docs/reports/step3_physics_heightmaps/s_128x128_t0300_colorbar.png` |

### Comparison vs Step 2 (advisory — yielding added, not a regression test)

#### Grid 128×128 — comparison vs Step 1

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 49.381 | 69.898 | ×1.42 |
| CG iters / linear solve (mean) | 117.2 | 111.4 | ×0.95 [idéal] |
| S mass drift (relative) | 4.663e-15 | -2.220e-15 | gate 1e-10 |
| max \|mean(vx)\| | 8.261e-23 | 8.758e-23 | bruit machine |
| max \|mean(vy)\| | 2.057e-22 | 2.708e-22 | bruit machine |

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
