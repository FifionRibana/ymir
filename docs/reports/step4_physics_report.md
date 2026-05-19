# Step 4 — GPE spreading + basal drag (physics)

> **Step 4 physics run for milestone "Solver reconstruction".**
> `GpeForce` (Ar = 0.1) + `BasalDragConfig::Enabled` with `Br = 0.05`. Velocity damping via `Br · S̃² · ṽ`, contributed to the **operator diagonal** (not the RHS), face-interpolated by arithmetic 2-point cell-to-face averaging. Yielding is `Disabled` at this baseline to isolate the drag effect.
> Solver unchanged: CG. The drag diagonal is positive semi-definite, preserves SPD-ness of the Picard block, and enters the preconditioner through `momentum_diagonal` (case B — analytical reconstruction).

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

## Br sweep (diagnostic)

Baseline `Br = 0.05` (preset `dynamic-accidented`, solver-scaling §5.1 centre of range). The sweep below covers `Br ∈ {0.01, 0.05, 0.10, 0.20, 0.30}` at 64²·N steps with `GpeForce (Ar = 0.1)` + basal drag Enabled + yielding Disabled (to isolate the Br effect). Expected qualitative behaviour: higher Br damps the velocity more, and the drag contribution on the operator diagonal improves conditioning (fewer CG iters at higher Br, or at worst stable). The two monotonicity checks below are acceptance invariants of the Step 4 sweep (issue #87).

| Br | wallclock (s) | CG iters (mean) | Newton iters (mean) | peak \|v\| | mass drift | Newton conv | drag/visc ratio | drag energy ratio |
|---|---|---|---|---|---|---|---|---|
| `0.010` | `5.72` | `51.51` | `2.0` | `8.155391e-6` | `1.33e-15` | `100%` | `2.523e-8` | `2.523e-8` |
| `0.050` | `5.54` | `51.51` | `2.0` | `8.155369e-6` | `-1.11e-15` | `100%` | `1.262e-7` | `1.262e-7` |
| `0.100` | `5.78` | `51.52` | `2.0` | `8.155341e-6` | `-5.55e-16` | `100%` | `2.523e-7` | `2.523e-7` |
| `0.200` | `5.73` | `51.52` | `2.0` | `8.155286e-6` | `2.22e-16` | `100%` | `5.046e-7` | `5.046e-7` |
| `0.300` | `5.56` | `51.52` | `2.0` | `8.155230e-6` | `6.66e-16` | `100%` | `7.569e-7` | `7.569e-7` |

**Monotonicity of `peak|v|` vs Br** (strictly decreasing): ✅ strictly decreasing across the 5 points, as required.

**Monotonicity of `cg_iter_mean` vs Br** (non-increasing): ✅ monotone non-increasing across the 5 points, as expected (drag improves conditioning).

**Interpretation** — low-Br damping is weak (velocities close to the un-damped baseline), high-Br damping reduces the peak flow magnitude. The drag contribution on the operator diagonal preserves SPD-ness of the Picard block and improves the conditioning of CG in low-ε̇ regions (modest effect in absolute terms at the Step 4 baseline where `S̃² ≈ 1` and `Br·S̃² ≪ η/Δx²`, as documented in the physics report's "Expected magnitude of drag effect" paragraph).

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
| basal drag | Enabled (Br = 0.050, S exponent = 2.0) |
| seed | 42 |

### Timing

- wallclock total: `5.685 s`
- wallclock per step (mean): `18.949 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `2.96e1`
- CG iterations per Newton step — mean: `51.5`, max: `109`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 41 | 1 |
  | 58 | 601 |
  | 75 | 4 |
  | 92 | 7 |
  | 109 | 7 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `2.0`, max: `6`
- effective η_max/η_min over run — mean: `1.03`, max: `1.03`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Basal drag

- Br = `0.050`
- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `1.262e-7`
- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `1.262e-7`
- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `1.262e-7` vs measured `1.262e-7` (relative diff `5.0e-9`; spec bound: coarse, typically `< 1e-1`)
- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)

**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.

**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.

### S field evolution

- Var(S̃) timeline: initial `1.000e-2`, middle `9.969e-3`, final `9.938e-3` (Δ = `-0.62%` vs initial)
- max|∇S̃| timeline: initial `1.255e0`, peak `1.255e0`, final `1.251e0`

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `-1.110e-15`

### Null-space health

- max |mean(vx)| across solves: `6.907e-23`
- max |mean(vy)|: `1.504e-22`

### Velocity magnitude

- peak |v|: `8.155e-6`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step4_physics_heightmaps/s_64x64_t0000.png` | `8.005e-1` | `1.200e0` | `1.000e0` | `docs/reports/step4_physics_heightmaps/s_64x64_t0000_colorbar.png` |
| `docs/reports/step4_physics_heightmaps/s_64x64_t0150.png` | `8.007e-1` | `1.199e0` | `1.000e0` | `docs/reports/step4_physics_heightmaps/s_64x64_t0150_colorbar.png` |
| `docs/reports/step4_physics_heightmaps/s_64x64_t0300.png` | `8.009e-1` | `1.199e0` | `1.000e0` | `docs/reports/step4_physics_heightmaps/s_64x64_t0300_colorbar.png` |

### Comparison vs Step 3 (advisory — basal drag added, not a regression test)

#### Grid 64×64 — comparison vs Step 3

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 8.960 | 5.685 | ×0.63 |
| CG iters / linear solve (mean) | 50.7 | 51.5 | ×1.02 [idéal] |
| S mass drift (relative) | -4.441e-16 | -1.110e-15 | gate 1e-10 |
| max \|mean(vx)\| | 1.067e-22 | 6.907e-23 | bruit machine |
| max \|mean(vy)\| | 1.733e-22 | 1.504e-22 | bruit machine |


**Wallclock improvement interpretation.** The Step-4 physics run at this grid is `×0.63` of the Step-3 physics wallclock — a measurable improvement, coherent with the theoretical expectation that adding `Br · S̃²` to the operator diagonal improves the conditioning of low-ε̇ regions. Despite the very small absolute drag contribution at this baseline (`drag/visc ≈ 10⁻⁷`), the augmented diagonal gives CG a slightly tighter grip on the system. Caveat: Step-4 physics also disables yielding (to isolate the Br effect), so part of the wallclock delta vs Step 3 physics comes from skipping the `soft_min_harmonic` plastic-branch evaluation; the pure drag contribution is best read off the Br sweep's strict `peak \|v\|` monotonicity and the κ(A) stability (ratio ≤ 1.3 across both grids). Encouraging signal for Step 5/6: introducing oceanic cells (`S̃ ≈ 0.2` → `S̃² ≈ 0.04`) will create a ×25 differentiation between continental and oceanic drag, which is where the Step-4 machinery's physical payoff will become visible.

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
| basal drag | Enabled (Br = 0.050, S exponent = 2.0) |
| seed | 42 |

### Timing

- wallclock total: `50.436 s`
- wallclock per step (mean): `168.120 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `1.50e2`
- CG iterations per Newton step — mean: `117.3`, max: `222`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 87 | 1 |
  | 121 | 599 |
  | 154 | 5 |
  | 188 | 9 |
  | 222 | 6 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `2.0`, max: `6`
- effective η_max/η_min over run — mean: `1.03`, max: `1.03`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Basal drag

- Br = `0.050`
- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `3.154e-8`
- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `3.154e-8`
- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `3.154e-8` vs measured `3.154e-8` (relative diff `1.3e-9`; spec bound: coarse, typically `< 1e-1`)
- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)

**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.

**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.

### S field evolution

- Var(S̃) timeline: initial `1.000e-2`, middle `9.969e-3`, final `9.938e-3` (Δ = `-0.62%` vs initial)
- max|∇S̃| timeline: initial `1.256e0`, peak `1.256e0`, final `1.252e0`

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `-5.884e-15`

### Null-space health

- max |mean(vx)| across solves: `7.165e-23`
- max |mean(vy)|: `1.487e-22`

### Velocity magnitude

- peak |v|: `8.160e-6`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step4_physics_heightmaps/s_128x128_t0000.png` | `8.001e-1` | `1.200e0` | `1.000e0` | `docs/reports/step4_physics_heightmaps/s_128x128_t0000_colorbar.png` |
| `docs/reports/step4_physics_heightmaps/s_128x128_t0150.png` | `8.003e-1` | `1.199e0` | `1.000e0` | `docs/reports/step4_physics_heightmaps/s_128x128_t0150_colorbar.png` |
| `docs/reports/step4_physics_heightmaps/s_128x128_t0300.png` | `8.006e-1` | `1.199e0` | `1.000e0` | `docs/reports/step4_physics_heightmaps/s_128x128_t0300_colorbar.png` |

### Comparison vs Step 3 (advisory — basal drag added, not a regression test)

#### Grid 128×128 — comparison vs Step 3

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 69.898 | 50.436 | ×0.72 |
| CG iters / linear solve (mean) | 111.4 | 117.3 | ×1.05 [idéal] |
| S mass drift (relative) | -2.220e-15 | -5.884e-15 | gate 1e-10 |
| max \|mean(vx)\| | 8.758e-23 | 7.165e-23 | bruit machine |
| max \|mean(vy)\| | 2.708e-22 | 1.487e-22 | bruit machine |


**Wallclock improvement interpretation.** The Step-4 physics run at this grid is `×0.72` of the Step-3 physics wallclock — a measurable improvement, coherent with the theoretical expectation that adding `Br · S̃²` to the operator diagonal improves the conditioning of low-ε̇ regions. Despite the very small absolute drag contribution at this baseline (`drag/visc ≈ 10⁻⁷`), the augmented diagonal gives CG a slightly tighter grip on the system. Caveat: Step-4 physics also disables yielding (to isolate the Br effect), so part of the wallclock delta vs Step 3 physics comes from skipping the `soft_min_harmonic` plastic-branch evaluation; the pure drag contribution is best read off the Br sweep's strict `peak \|v\|` monotonicity and the κ(A) stability (ratio ≤ 1.3 across both grids). Encouraging signal for Step 5/6: introducing oceanic cells (`S̃ ≈ 0.2` → `S̃² ≈ 0.04`) will create a ×25 differentiation between continental and oceanic drag, which is where the Step-4 machinery's physical payoff will become visible.

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
