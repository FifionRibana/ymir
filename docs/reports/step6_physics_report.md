# Step 6 — Dynamic boundaries, Voronoi plates, conservative recycling (physics)

> **Step 6 physics run for milestone "Solver reconstruction".**
> `GpeForce (Ar = 0.1)` + `YieldingConfig::Enabled (Bi = 0.15)` + `BasalDragConfig::Enabled (Br = 0.05)` + `BoundaryConfig::Enabled` with **Voronoi tessellation** (num_plates=8, continental_ratio=0.3) and **Closed-mode recycling** (arc=0.15, coll_v=0.03, rift_v=0.02, spread=0.80, mantle_loss=0.00, mantle_delay=20 steps).
> Boundary flags are **detected dynamically per step** from the velocity divergence (threshold=1e-4). Subducted mass feeds arc/collision/rift volcanism immediately and mid-ocean spreading through the delayed buffer — `k_spread` disappears as a rate in favour of the spread fraction of the recycled budget. No taper at oceanic/continental interfaces: #78 monitored as trajectory telemetry at t ∈ {1, 10, 50, 150, 300}·Δt.
> Solver unchanged: CG on the Picard block, Newton outer. Conservation invariant (Step 6): `|Δmass_obs + mantle_loss + buffer_fill + pending_immediate − clamp_flux| / initial_mass < 1e-6` with `mantle_loss_fraction = 0`.

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
| basal drag | Enabled (Br = 0.050, S exponent = 2.0) |
| boundary config | Enabled [Closed (arc=0.150, coll_v=0.030, rift_v=0.020, spread=0.800, mantle_loss=0.000, delay=20 steps)] (layout='voronoi_seed42_n8', k_sub=0.500, k_arc=0.000, k_spread=0.000, k_coll-v=0.000, k_rift-v=0.000) |
| boundary layout | voronoi_seed42_n8 |
| seed | 42 |

### Timing

- wallclock total: `34.448 s`
- wallclock per step (mean): `114.826 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `1.85e2`
- CG iterations per Newton step — mean: `129.6`, max: `2000`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 430 | 1514 |
  | 822 | 0 |
  | 1215 | 0 |
  | 1607 | 0 |
  | 2000 | 15 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `5.0`, max: `8`
- effective η_max/η_min over run — mean: `1.33`, max: `1.33`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Plastic yielding

- Bi = `0.150`
- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `0.000`
- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `0.477`
- Definition notes: the `< 0.5·η_visc` criterion captures "yielding dominant", not "yielding present anywhere" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).

### Strain-rate regime diagnostic (floor-domination)

- `ε̇_min` (regularisation floor): `1.000e-3`
- `ε̇_II` at final timestep: mean = `1.383e-4`, max = `3.979e-4`
- Fraction of cells with `ε̇_II < 10·ε̇_min = 1.000e-2` at final timestep: `1.000` (1.0 = everywhere in the floor-dominated band)
- `max(ε̇_II) / ε̇_min` = `0.40` — ratio of the strongest strain-rate cell to the regularisation floor.

**Verdict: floor-dominated.** `ε̇_II` lies below `ε̇_min` over most of the domain; the viscous and plastic branches both saturate at their floor values. The analytic criterion for yielding dominance in this regime is

```
  Bi < ε̇_min^(1/3)   (n = 3)
  ⟺ Bi < 0.100
```

with the default scales. The baseline `Bi = 0.150` sits above this threshold, so `yielding_cell_fraction = 0` is the **expected** diagnostic outcome for Ar = 0.1 + GPE-only forcing. The Bi sweep at `Bi ≤ 0.10` crosses the threshold and shows `yielding_cell_fraction = 1.0`, confirming the yielding mechanism is wired correctly — it simply is not activated by the weak GPE regime at this baseline.

**Anticipated cross-over at later steps.** Mechanisms introduced in Steps 4 (basal drag), 5 (boundary sources), 7 (slab pull), and 8 (mantle flow) inject energy at faster time scales than GPE and should push `ε̇_II` into the O(1) range in active zones. As soon as active-zone `ε̇_II > ε̇_min`, the `Bi = 0.15` criterion for yielding dominance will hold locally and `yielding_cell_fraction > 0` should appear naturally. If `yielding_cell_fraction` is still 0 after Step 7 (slab pull, which acts at τ*/Sp ≈ 10–60 Myr), the coupling between source mechanisms and ε̇ is under-dimensioned and warrants a remontée — **this is flagged as a checkpoint** for the Step 4, 5, 7, 8 physics reports.

Basal drag (Step 4) is dissipative and may *not* raise `ε̇_II`; the threshold check starts in earnest at Step 5 (boundary sources inject mass and create strain) and carries through Steps 7–8.

### Basal drag

- Br = `0.050`
- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `9.406e-8`
- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `9.406e-8`
- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `9.406e-8` vs measured `9.406e-8` (relative diff `1.1e-7`; spec bound: coarse, typically `< 1e-1`)
- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)

**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.

**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.

### Boundary source/sink diagnostics

- `s_oceanic_mean` = `0.2003` (std `0.0207`)  — target `[0.18, 0.22]` post-calibration
- `s_continental_interior_mean` = `0.8279` (std `0.0163`)  — target `[0.9, 1.1]`
- `boundary_type_diversity` = `1` (number of distinct mechanisms active on the run)
- `clamp_activation_fraction` — mean `0.000e0`, max `0.000e0` (healthy: mean < 1%, max < 5%)
- `∫Q dt dA` = `-2.446e-4`; `∫clamp_flux dt dA` = `0.000e0`
- `mass_balance_residual` = `1.150e-13` (issue #89 D5; acceptance `< 1%`)

### Voronoi plate geometry

- distinct plate_count = `8` (expected 8 for `num_plates=8`)
- plate_type_distribution (oceanic, continental) = `(0.554, 0.446)` — target continental ∈ [0.15, 0.45]

### Boundary dynamics (dynamic detection per step)

- `boundary_flag_transition_rate` — mean `0.000e0`, max `0.000e0`
  - Fraction of cells whose `boundary_flag` changed vs the previous step. Telemetry only — no acceptance. Expected transient spike early in the run (flags emerging from `None` as the first Stokes solves produce non-trivial divergence), then stabilisation.

### Recycling health (Closed mode)

- `recycling_buffer_fill` — mean `1.184e-4`, max `2.358e-4`, final `2.358e-4`
- `immediate_pending_max` over run = `8.842e-6`, final sum = `8.842e-6`
- `clamp_activation_during_spinup_max` = `0.000e0` (target 0 — clamp should not fire during the buffer fill-up)

### Mass balance (Step 6 closed recycling, 5 components)

- Δmass_observed (dimensionless, S̃ sum): initial `2.281932e3`, final `2.280930e3`, Δ = `-1.002e0`
- `buffer_fill_final` (cell-area units) = `2.358e-4`
- `pending_immediate_final` (cell-area units) = `8.842e-6`
- `clamp_flux_integral` (cell-area units) = `0.000e0`
- `mantle_loss_integral` (cell-area units) = `0.000e0` (zero when mantle_loss_fraction=0)
- **`mass_conservation_residual` = `2.800e-17`** (target `< 1e-6`)

Formula: `|Δmass_obs + mantle_loss + buffer_fill + pending − clamp_flux| / initial_mass`. All five components are tracked; the residual is the absolute sum divided by `initial_mass`. A `< 1e-6` residual means the pipeline is mass-exact at machine precision; all deviations from exact conservation are accounted for by the known components (loss + in-transit buffer mass + rollover pending + clamp artificial flux).

### Issue #78 trajectory (5 instants: t ∈ {1, 10, 50, 150, 300}·Δt)

| step | max\|∇S̃\|_interface | max\|∇S̃\|_global | peak\|f_GPE\|_interface | peak\|f_GPE\|_global | buffer_fill |
|---|---|---|---|---|---|
| `1` | `8.669e1` | `8.669e1` | `6.210e0` | `6.210e0` | `7.873e-7` |
| `10` | `8.664e1` | `8.664e1` | `6.208e0` | `6.208e0` | `7.873e-6` |
| `50` | `8.642e1` | `8.642e1` | `6.199e0` | `6.199e0` | `3.936e-5` |
| `150` | `8.585e1` | `8.585e1` | `6.175e0` | `6.175e0` | `1.180e-4` |
| `300` | `8.502e1` | `8.502e1` | `6.138e0` | `6.138e0` | `2.358e-4` |

**Interpretation.** No taper was applied at the Voronoi oceanic/continental interfaces (per Step 6 D5 — #78 is tested, not contoured). A spike that appears at step 1 and damps by step 50 is a transient artefact of the raw contrast; a spike that grows monotonically across the 5 instants is a real signal that #78 has activated and must be addressed before Step 7. **Absolute critical threshold**: `peak|f_GPE| > 100` at any instant = red-flag bug.

### Note on OceanicSubduction drain symmetry

When two oceanic cells meet at a convergent boundary, both are flagged `OceanicSubduction` and both contribute to `Q_sub`. This effectively doubles the local drain compared to Oceanic/Continental subduction (where only the oceanic cell drains). This is an assumed approximation in the absence of an age field (Step 10) that would resolve which cell actually subducts. The mass budget stays correct because the combined drain feeds the same recycling pool: total mass conservation is satisfied independently of which side is drained. To be refined at Step 10.

### Yielding activation checkpoint (Step 6)

- Bi = `0.150`, `yielding_cell_fraction_max` = `0.000`

**Checkpoint status: still 0 at Step 6.** Step 6 was the last step before slab-pull forcing that could plausibly activate yielding without an external mechanism. `yielding_cell_fraction = 0` here means the checkpoint migrates to Step 7 — slab-pull at `τ*/Sp ≈ 10–60 Myr` is the expected activation trigger. If still 0 at Step 7, remontée required.

### Preconditioner surveillance (continued from Step 5)

Step 5 physics: CG mean = 108.5 (64²) / 205.0 (128²), ≈ 2× Step 4. Step 6 adds Voronoi interfaces (sharper contrasts, more heterogeneity). If the CG ratio vs Step 5 is ≤ 2× (i.e., vs Step 4 ≤ 4×), continue surveillance. If > 10× Step 4, the preconditioner has reached its usable limit and the maintenance task (block-Jacobi / ILU(0)) should be scheduled before Step 7.

### S field evolution

- Var(S̃) timeline: initial `1.633e-1`, middle `1.629e-1`, final `1.626e-1` (Δ = `-0.46%` vs initial)
- max|∇S̃| timeline: initial `8.670e1`, peak `8.670e1`, final `8.502e1`

### Mass conservation of S

- initial mass: `2.281931793e3`
- final mass: `2.280929763e3`
- relative drift: `-4.391e-4`

### Null-space health

- max |mean(vx)| across solves: `3.220e-22`
- max |mean(vy)|: `5.954e-21`

### Velocity magnitude

- peak |v|: `3.602e-5`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step6_physics_heightmaps/s_64x64_t0000.png` | `1.601e-1` | `1.198e0` | `5.571e-1` | `docs/reports/step6_physics_heightmaps/s_64x64_t0000_colorbar.png` |
| `docs/reports/step6_physics_heightmaps/s_64x64_t0150.png` | `1.599e-1` | `1.197e0` | `5.570e-1` | `docs/reports/step6_physics_heightmaps/s_64x64_t0150_colorbar.png` |
| `docs/reports/step6_physics_heightmaps/s_64x64_t0300.png` | `1.597e-1` | `1.196e0` | `5.569e-1` | `docs/reports/step6_physics_heightmaps/s_64x64_t0300_colorbar.png` |

### Comparison vs Step 5 physics (advisory — Voronoi + dynamic detection + Closed recycling added)

#### Grid 64×64 — comparison vs Step 5 physics

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 35.058 | 34.448 | ×0.98 |
| CG iters / linear solve (mean) | 108.5 | 129.6 | ×1.19 [idéal] |
| S mass drift (relative) | 6.342e-3 | -4.391e-4 | gate 1e-10 |
| max \|mean(vx)\| | 8.148e-23 | 3.220e-22 | bruit machine |
| max \|mean(vy)\| | 1.754e-20 | 5.954e-21 | bruit machine |

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
| boundary config | Enabled [Closed (arc=0.150, coll_v=0.030, rift_v=0.020, spread=0.800, mantle_loss=0.000, delay=20 steps)] (layout='voronoi_seed42_n8', k_sub=0.500, k_arc=0.000, k_spread=0.000, k_coll-v=0.000, k_rift-v=0.000) |
| boundary layout | voronoi_seed42_n8 |
| seed | 42 |

### Timing

- wallclock total: `337.492 s`
- wallclock per step (mean): `1124.972 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `6.31e2`
- CG iterations per Newton step — mean: `240.4`, max: `2000`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 445 | 1808 |
  | 834 | 6 |
  | 1222 | 0 |
  | 1611 | 0 |
  | 2000 | 15 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `6.0`, max: `8`
- effective η_max/η_min over run — mean: `1.38`, max: `1.38`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Plastic yielding

- Bi = `0.150`
- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `0.000`
- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `0.477`
- Definition notes: the `< 0.5·η_visc` criterion captures "yielding dominant", not "yielding present anywhere" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).

### Strain-rate regime diagnostic (floor-domination)

- `ε̇_min` (regularisation floor): `1.000e-3`
- `ε̇_II` at final timestep: mean = `1.383e-4`, max = `4.563e-4`
- Fraction of cells with `ε̇_II < 10·ε̇_min = 1.000e-2` at final timestep: `1.000` (1.0 = everywhere in the floor-dominated band)
- `max(ε̇_II) / ε̇_min` = `0.46` — ratio of the strongest strain-rate cell to the regularisation floor.

**Verdict: floor-dominated.** `ε̇_II` lies below `ε̇_min` over most of the domain; the viscous and plastic branches both saturate at their floor values. The analytic criterion for yielding dominance in this regime is

```
  Bi < ε̇_min^(1/3)   (n = 3)
  ⟺ Bi < 0.100
```

with the default scales. The baseline `Bi = 0.150` sits above this threshold, so `yielding_cell_fraction = 0` is the **expected** diagnostic outcome for Ar = 0.1 + GPE-only forcing. The Bi sweep at `Bi ≤ 0.10` crosses the threshold and shows `yielding_cell_fraction = 1.0`, confirming the yielding mechanism is wired correctly — it simply is not activated by the weak GPE regime at this baseline.

**Anticipated cross-over at later steps.** Mechanisms introduced in Steps 4 (basal drag), 5 (boundary sources), 7 (slab pull), and 8 (mantle flow) inject energy at faster time scales than GPE and should push `ε̇_II` into the O(1) range in active zones. As soon as active-zone `ε̇_II > ε̇_min`, the `Bi = 0.15` criterion for yielding dominance will hold locally and `yielding_cell_fraction > 0` should appear naturally. If `yielding_cell_fraction` is still 0 after Step 7 (slab pull, which acts at τ*/Sp ≈ 10–60 Myr), the coupling between source mechanisms and ε̇ is under-dimensioned and warrants a remontée — **this is flagged as a checkpoint** for the Step 4, 5, 7, 8 physics reports.

Basal drag (Step 4) is dissipative and may *not* raise `ε̇_II`; the threshold check starts in earnest at Step 5 (boundary sources inject mass and create strain) and carries through Steps 7–8.

### Basal drag

- Br = `0.050`
- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `2.349e-8`
- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `2.349e-8`
- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `2.349e-8` vs measured `2.349e-8` (relative diff `2.7e-8`; spec bound: coarse, typically `< 1e-1`)
- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)

**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.

**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.

### Boundary source/sink diagnostics

- `s_oceanic_mean` = `0.2003` (std `0.0209`)  — target `[0.18, 0.22]` post-calibration
- `s_continental_interior_mean` = `0.8272` (std `0.0161`)  — target `[0.9, 1.1]`
- `boundary_type_diversity` = `1` (number of distinct mechanisms active on the run)
- `clamp_activation_fraction` — mean `0.000e0`, max `0.000e0` (healthy: mean < 1%, max < 5%)
- `∫Q dt dA` = `-2.447e-4`; `∫clamp_flux dt dA` = `0.000e0`
- `mass_balance_residual` = `1.602e-11` (issue #89 D5; acceptance `< 1%`)

### Voronoi plate geometry

- distinct plate_count = `8` (expected 8 for `num_plates=8`)
- plate_type_distribution (oceanic, continental) = `(0.555, 0.445)` — target continental ∈ [0.15, 0.45]

### Boundary dynamics (dynamic detection per step)

- `boundary_flag_transition_rate` — mean `4.083e-7`, max `6.104e-5`
  - Fraction of cells whose `boundary_flag` changed vs the previous step. Telemetry only — no acceptance. Expected transient spike early in the run (flags emerging from `None` as the first Stokes solves produce non-trivial divergence), then stabilisation.

### Recycling health (Closed mode)

- `recycling_buffer_fill` — mean `1.184e-4`, max `2.359e-4`, final `2.359e-4`
- `immediate_pending_max` over run = `8.845e-6`, final sum = `8.845e-6`
- `clamp_activation_during_spinup_max` = `0.000e0` (target 0 — clamp should not fire during the buffer fill-up)

### Mass balance (Step 6 closed recycling, 5 components)

- Δmass_observed (dimensionless, S̃ sum): initial `9.124226e3`, final `9.120217e3`, Δ = `-4.009e0`
- `buffer_fill_final` (cell-area units) = `2.359e-4`
- `pending_immediate_final` (cell-area units) = `8.845e-6`
- `clamp_flux_integral` (cell-area units) = `0.000e0`
- `mantle_loss_integral` (cell-area units) = `0.000e0` (zero when mantle_loss_fraction=0)
- **`mass_conservation_residual` = `3.921e-15`** (target `< 1e-6`)

Formula: `|Δmass_obs + mantle_loss + buffer_fill + pending − clamp_flux| / initial_mass`. All five components are tracked; the residual is the absolute sum divided by `initial_mass`. A `< 1e-6` residual means the pipeline is mass-exact at machine precision; all deviations from exact conservation are accounted for by the known components (loss + in-transit buffer mass + rollover pending + clamp artificial flux).

### Issue #78 trajectory (5 instants: t ∈ {1, 10, 50, 150, 300}·Δt)

| step | max\|∇S̃\|_interface | max\|∇S̃\|_global | peak\|f_GPE\|_interface | peak\|f_GPE\|_global | buffer_fill |
|---|---|---|---|---|---|
| `1` | `1.734e2` | `1.734e2` | `1.246e1` | `1.246e1` | `7.876e-7` |
| `10` | `1.732e2` | `1.732e2` | `1.245e1` | `1.245e1` | `7.876e-6` |
| `50` | `1.724e2` | `1.724e2` | `1.242e1` | `1.242e1` | `3.937e-5` |
| `150` | `1.704e2` | `1.704e2` | `1.235e1` | `1.235e1` | `1.180e-4` |
| `300` | `1.673e2` | `1.673e2` | `1.223e1` | `1.223e1` | `2.359e-4` |

**Interpretation.** No taper was applied at the Voronoi oceanic/continental interfaces (per Step 6 D5 — #78 is tested, not contoured). A spike that appears at step 1 and damps by step 50 is a transient artefact of the raw contrast; a spike that grows monotonically across the 5 instants is a real signal that #78 has activated and must be addressed before Step 7. **Absolute critical threshold**: `peak|f_GPE| > 100` at any instant = red-flag bug.

### Note on OceanicSubduction drain symmetry

When two oceanic cells meet at a convergent boundary, both are flagged `OceanicSubduction` and both contribute to `Q_sub`. This effectively doubles the local drain compared to Oceanic/Continental subduction (where only the oceanic cell drains). This is an assumed approximation in the absence of an age field (Step 10) that would resolve which cell actually subducts. The mass budget stays correct because the combined drain feeds the same recycling pool: total mass conservation is satisfied independently of which side is drained. To be refined at Step 10.

### Yielding activation checkpoint (Step 6)

- Bi = `0.150`, `yielding_cell_fraction_max` = `0.000`

**Checkpoint status: still 0 at Step 6.** Step 6 was the last step before slab-pull forcing that could plausibly activate yielding without an external mechanism. `yielding_cell_fraction = 0` here means the checkpoint migrates to Step 7 — slab-pull at `τ*/Sp ≈ 10–60 Myr` is the expected activation trigger. If still 0 at Step 7, remontée required.

### Preconditioner surveillance (continued from Step 5)

Step 5 physics: CG mean = 108.5 (64²) / 205.0 (128²), ≈ 2× Step 4. Step 6 adds Voronoi interfaces (sharper contrasts, more heterogeneity). If the CG ratio vs Step 5 is ≤ 2× (i.e., vs Step 4 ≤ 4×), continue surveillance. If > 10× Step 4, the preconditioner has reached its usable limit and the maintenance task (block-Jacobi / ILU(0)) should be scheduled before Step 7.

### S field evolution

- Var(S̃) timeline: initial `1.633e-1`, middle `1.630e-1`, final `1.626e-1` (Δ = `-0.46%` vs initial)
- max|∇S̃| timeline: initial `1.735e2`, peak `1.735e2`, final `1.673e2`

### Mass conservation of S

- initial mass: `9.124226001e3`
- final mass: `9.120216559e3`
- relative drift: `-4.394e-4`

### Null-space health

- max |mean(vx)| across solves: `4.205e-22`
- max |mean(vy)|: `7.775e-21`

### Velocity magnitude

- peak |v|: `3.625e-5`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step6_physics_heightmaps/s_128x128_t0000.png` | `1.600e-1` | `1.198e0` | `5.569e-1` | `docs/reports/step6_physics_heightmaps/s_128x128_t0000_colorbar.png` |
| `docs/reports/step6_physics_heightmaps/s_128x128_t0150.png` | `1.598e-1` | `1.197e0` | `5.568e-1` | `docs/reports/step6_physics_heightmaps/s_128x128_t0150_colorbar.png` |
| `docs/reports/step6_physics_heightmaps/s_128x128_t0300.png` | `1.597e-1` | `1.197e0` | `5.567e-1` | `docs/reports/step6_physics_heightmaps/s_128x128_t0300_colorbar.png` |

### Comparison vs Step 5 physics (advisory — Voronoi + dynamic detection + Closed recycling added)

#### Grid 128×128 — comparison vs Step 5 physics

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 298.899 | 337.492 | ×1.13 |
| CG iters / linear solve (mean) | 205.0 | 240.4 | ×1.17 [idéal] |
| S mass drift (relative) | 3.198e-3 | -4.394e-4 | gate 1e-10 |
| max \|mean(vx)\| | 9.213e-23 | 4.205e-22 | bruit machine |
| max \|mean(vy)\| | 2.354e-20 | 7.775e-21 | bruit machine |

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
