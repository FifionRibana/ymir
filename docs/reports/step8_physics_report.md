# Step 8 — Mantle bootstrap validation (slab-pull held disabled pending co-calibration)

> **Step 8 physics run for milestone "Solver reconstruction".**
> Setup: **Step 6 physics base** (`GpeForce` + yielding Enabled + basal drag Enabled + Voronoi + dynamic detection + Closed recycling) plus `MantleConfig::Enabled` with baseline `(Mf = 1.0, coupling = 1.0, num_modes = 6, seed = 42, evolution_rate = 0)`. **Slab-pull is held Disabled** for this step — see §Slab+Mantle interaction instability finding below and the regression-convention exception in `tectonics_v2/README.md`.
> Formulation: `f_mantle = coupling · S̃ · (Mf · v_pattern − v_solved)`. The `-coupling · S̃ · v_solved` part is folded into the momentum-operator diagonal (same as basal drag Step 4) for exact self-consistency at every Newton outer iteration; the constant RHS part `coupling · S̃ · Mf · v_pattern` is assembled as a body force. Pattern is div-free by construction (staggered curl of a nodal Fourier stream function) and static at Step 8 (time evolution deferred per D6).
> **Yielding checkpoint STRICT — last chance.** Per the amplifier-vs-initiator revision at Step 7, mantle forcing is the INITIATOR of the mechanism hierarchy. Mantle-alone is the configuration that resolves the checkpoint; the slab+mantle interaction requires co-calibration deferred to a dedicated follow-up issue (see the finding section below).

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
| slab-pull | Disabled |
| mantle forcing | Enabled (Mf = 1.000, coupling = 1.000, num_modes = 6, seed = 42, evolution_rate = 0.000) |
| seed | 42 |

### Timing

- wallclock total: `1145.551 s`
- wallclock per step (mean): `3818.504 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `2.21e4`
- CG iterations per Newton step — mean: `1420.7`, max: `2000`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 419 | 657 |
  | 814 | 658 |
  | 1209 | 472 |
  | 1604 | 264 |
  | 2000 | 2667 |

### Newton (nonlinear) health

- outcome distribution — Converged: `99.3%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.7%`
- Newton outer iters per timestep — mean: `15.5`, max: `20`
- effective η_max/η_min over run — mean: `65803.49`, max: `106706.93`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Plastic yielding

- Bi = `0.150`
- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `0.998`
- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `27.650`
- Definition notes: the `< 0.5·η_visc` criterion captures "yielding dominant", not "yielding present anywhere" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).

### Strain-rate regime diagnostic (floor-domination)

- `ε̇_min` (regularisation floor): `1.000e-3`
- `ε̇_II` at final timestep: mean = `1.733e1`, max = `1.206e2`
- Fraction of cells with `ε̇_II < 10·ε̇_min = 1.000e-2` at final timestep: `0.026` (1.0 = everywhere in the floor-dominated band)
- `max(ε̇_II) / ε̇_min` = `120606.09` — ratio of the strongest strain-rate cell to the regularisation floor.

**Verdict:** partial floor-domination — `97.4%` of cells are above the `10·ε̇_min` threshold. `yielding_cell_fraction` should be roughly consistent with that active fraction.

### Basal drag

- Br = `0.050`
- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `1.158e-3`
- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `1.164e-3`
- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `1.163e-3` vs measured `1.158e-3` (relative diff `4.0e-3`; spec bound: coarse, typically `< 1e-1`)
- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)

**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.

**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.

### Boundary source/sink diagnostics

- `s_oceanic_mean` = `0.5235` (std `0.1775`)  — target `[0.18, 0.22]` post-calibration
- `s_continental_interior_mean` = `0.5109` (std `0.1205`)  — target `[0.9, 1.1]`
- `s_continental_collision_mean` = `0.6095` — telemetry only (orogen thickening, tracked through Steps 5-10)
- `boundary_type_diversity` = `4` (number of distinct mechanisms active on the run)
- `clamp_activation_fraction` — mean `3.196e-3`, max `5.859e-3` (healthy: mean < 1%, max < 5%)
- `∫Q dt dA` = `-1.195e-2`; `∫clamp_flux dt dA` = `4.401e-3`
- `mass_balance_residual` = `6.679e-14` (issue #89 D5; acceptance `< 1%`)

### Voronoi plate geometry

- distinct plate_count = `8` (expected 8 for `num_plates=8`)
- plate_type_distribution (oceanic, continental) = `(0.554, 0.446)` — target continental ∈ [0.15, 0.45]

### Boundary dynamics (dynamic detection per step)

- `boundary_flag_transition_rate` — mean `3.935e-3`, max `7.251e-2`
  - Fraction of cells whose `boundary_flag` changed vs the previous step. Telemetry only — no acceptance. Expected transient spike early in the run (flags emerging from `None` as the first Stokes solves produce non-trivial divergence), then stabilisation.
- flag counts **at step 1** (proving detection fired): None=`7`, Subduction=`63`, OceanicSubduction=`1249`, Rift=`1928`, ContinentalCollision=`849`
- flag counts **at final step**: None=`31`, Subduction=`39`, OceanicSubduction=`1072`, Rift=`2108`, ContinentalCollision=`846`

### Recycling health (Closed mode)

- `recycling_buffer_fill` — mean `1.087e-2`, max `1.238e-2`, final `1.195e-2`
- `immediate_pending_max` over run = `0.000e0`, final sum = `0.000e0`
- `clamp_activation_during_spinup_max` = `2.441e-4` (target 0 — clamp should not fire during the buffer fill-up)

### Mass balance (Step 6 closed recycling, 5 components)

- Δmass_observed (dimensionless, S̃ sum): initial `2.281932e3`, final `2.251013e3`, Δ = `-3.092e1`
- `buffer_fill_final` (cell-area units) = `1.195e-2`
- `pending_immediate_final` (cell-area units) = `0.000e0`
- `clamp_flux_integral` (cell-area units) = `4.401e-3`
- `mantle_loss_integral` (cell-area units) = `0.000e0` (zero when mantle_loss_fraction=0)
- **`mass_conservation_residual` = `1.080e-15`** (target `< 1e-6`)

Formula: `|Δmass_obs + mantle_loss + buffer_fill + pending − clamp_flux| / initial_mass`. All five components are tracked; the residual is the absolute sum divided by `initial_mass`. A `< 1e-6` residual means the pipeline is mass-exact at machine precision; all deviations from exact conservation are accounted for by the known components (loss + in-transit buffer mass + rollover pending + clamp artificial flux).

### Issue #78 trajectory (5 instants: t ∈ {1, 10, 50, 150, 300}·Δt)

| step | max\|∇S̃\|_interface | max\|∇S̃\|_global | peak\|f_GPE\|_interface | peak\|f_GPE\|_global | buffer_fill |
|---|---|---|---|---|---|
| `1` | `8.658e1` | `8.658e1` | `6.250e0` | `6.250e0` | `2.961e-4` |
| `10` | `8.402e1` | `8.402e1` | `6.643e0` | `6.643e0` | `3.804e-3` |
| `50` | `8.130e1` | `8.130e1` | `6.409e0` | `6.409e0` | `1.085e-2` |
| `150` | `4.314e1` | `4.692e1` | `4.269e0` | `4.269e0` | `1.154e-2` |
| `300` | `4.362e1` | `4.684e1` | `5.344e0` | `5.344e0` | `1.195e-2` |

**Interpretation.** No taper was applied at the Voronoi oceanic/continental interfaces (per Step 6 D5 — #78 is tested, not contoured). A spike that appears at step 1 and damps by step 50 is a transient artefact of the raw contrast; a spike that grows monotonically across the 5 instants is a real signal that #78 has activated and must be addressed before Step 7. **Absolute critical threshold**: `peak|f_GPE| > 100` at any instant = red-flag bug.

### Continental mass balance (Closed mode)

Continental cells cannot drain via Q_sub (Step 5 invariant: Q_sub fires only on `(Oceanic, is_subduction())` cells). Continental thickness changes come from three sources: (1) **immediate recycling returns** (`Q_arc + Q_coll_v + Q_rift_v`, all applied to continental eligible cells), (2) **advection** across the continental/oceanic boundary, driven by GPE spreading, and (3) **no other Q contribution**.

- `M_sub_total` (integrated drain, all oceanic subducting cells): `2.216e-1`
- `∫Q_arc dt dA` (continental return, arc volcanism): `3.324e-2` — fraction `0.150` of M_sub
- `∫Q_coll_v dt dA` (continental return, collision volcanism): `6.648e-3` — fraction `0.030`
- `∫Q_rift_v dt dA` (continental return, rift volcanism): `4.432e-3` — fraction `0.020`
- Total continental return: `4.432e-2` — fraction `0.200` of M_sub
- `∫Q_spread dt dA` (oceanic return, mid-ocean ridges): `1.653e-1` — fraction `0.746` of M_sub

`s_continental_interior_mean = 0.5109` at end of run (target `[0.9, 1.1]`).

**Interpretation** — with default fractions `(arc 0.15, coll_v 0.03, rift_v 0.02, spread 0.80)` the immediate continental return is **20% of M_sub** while 80% is routed through the delayed buffer to OCEANIC ridges. Net continental balance depends on (a) how much mass the Voronoi advection pushes across the continental/oceanic boundary, and (b) how evenly the 20% immediate return is distributed over the continental cell population.

If `s_continental_interior_mean < 0.9`, the interpretation is that the **continental set is a net mass exporter** to the oceanic set via advection — GPE drives flow away from high-S continental cells toward the thinner oceanic strip, and only 20% of the subducted mass returns to continental via arc + collision + rift volcanism. Global mass is conserved (the spread_fraction=0.80 returns to oceanic cells via the delayed buffer), but the continental/oceanic **partition** is not invariant.

This is expected physics, not a bug. The `[0.9, 1.1]` target band from issue #90 was set against the Step 5 static layout (where continental cells sat in spatial isolation from subduction). With a Voronoi tessellation where continental patches are surrounded by advecting oceanic zones, mass redistribution over 300 steps is larger — the continental mean drifts toward a new Voronoi-specific equilibrium that is not 1.0. Adjusting the acceptance band to reflect Voronoi dynamics is follow-up work; the mass budget itself (`mass_conservation_residual < 1e-6`) holds unambiguously.

### Note on OceanicSubduction drain symmetry

When two oceanic cells meet at a convergent boundary, both are flagged `OceanicSubduction` and both contribute to `Q_sub`. This effectively doubles the local drain compared to Oceanic/Continental subduction (where only the oceanic cell drains). This is an assumed approximation in the absence of an age field (Step 10) that would resolve which cell actually subducts. The mass budget stays correct because the combined drain feeds the same recycling pool: total mass conservation is satisfied independently of which side is drained. To be refined at Step 10.

### Yielding activation checkpoint (Step 6)

- Bi = `0.150`, `yielding_cell_fraction_max` = `0.998`

**Checkpoint status: ✅ activated at Step 6.** Dynamic boundary geometry + closed recycling produced enough convergent strain at some cells to push `ε̇_II > ε̇_min` locally, crossing the Bi threshold. The mechanism is wired and active; expect further growth at Steps 7 (slab pull) and 8 (mantle forcing).

### Preconditioner surveillance (continued from Step 5)

Step 5 physics: CG mean = 108.5 (64²) / 205.0 (128²), ≈ 2× Step 4. Step 6 adds Voronoi interfaces (sharper contrasts, more heterogeneity). If the CG ratio vs Step 5 is ≤ 2× (i.e., vs Step 4 ≤ 4×), continue surveillance. If > 10× Step 4, the preconditioner has reached its usable limit and the maintenance task (block-Jacobi / ILU(0)) should be scheduled before Step 7.

### Mantle bootstrap (Step 8)

- Mf = `1.000` (target band [0.3, 2.0] per §4.9)
- coupling = `1.000` (target band [0.1, 10.0])
- num_modes = `6`, seed = `42`

- `peak|v_mantle|` (= Mf · peak|v_pattern|) = `1.653e1`
- `peak|v_solved|` (max over run) = `9.552e0`
- `v_solved_to_v_mantle_alignment` (mean of `<v, Mf·v_m>/|Mf·v_m|²`) = `0.240`
- `div_v_mantle_max` = `1.137e-13` (strict acceptance `< 1e-10`)

**Bootstrap: ✅ system escaped floor-domination.** `peak|v_solved|` exceeds 0.1 — three or more orders of magnitude above the Step 7 baseline (3.6e-5). Mantle forcing is performing its role as the mechanism-hierarchy initiator (see §4.8 activation-regime note).

### Force hierarchy (Step 8)

- `peak|f_GPE|` = `7.350e0`
- `peak|f_slab|` = `0.000e0`
- `peak|f_mantle|` = `2.012e1`
- `f_mantle / f_GPE` (mean per step) = `2.604e0`

**Interpretation bands** (telemetry, not acceptance — except the pathological case):
- `f_mantle ≫ f_GPE` (ratio ≥ 10): mantle bootstrapped. Success.
- `f_mantle ~ f_slab` (ratio 0.1–10): healthy coupling.
- `f_slab ≫ f_mantle` (ratio < 0.1): non-pathological, document.
- `f_mantle ≪ f_GPE` (ratio < 0.1): PATHOLOGICAL — correlates with bootstrap failed, remontée required.

### Yielding activation (Step 8 — STRICT, last chance)

- Bi = `0.150`, `yielding_cell_fraction_max` = `9.976e-1`
- `max(ε̇_II) / ε̇_min` = `1.586e5` (floor-dominated if ≤ 1)

**Yielding activation: ✅ RESOLVED.** The checkpoint transported since Step 3 (and strictly enforced here as last-chance per the Step 7 revision) is met: yielding fires in a non-marginal fraction of cells. Mantle forcing has bootstrapped `ε̇_II` above the regularisation floor locally, and the Bingham criterion (`η_eff < 0.5 · η_visc`) captures the resulting yielding-dominated regime. The mechanism hierarchy is confirmed.

### Slab+Mantle interaction instability finding (Step 8)

The Step 8 baseline above holds **slab-pull Disabled** by deliberate choice. During Step 8 development, running the nominal spec configuration (Step 7 physics + mantle Enabled) produced catastrophic numerical divergence within 15–20 timesteps at 64² × Mf=1.0, `coupling=1.0`, slab-pull at Step 7's `(Sp=1.5, τ_slab=0.5, k_slab_accum=1.0)`. The runaway is physically real (captured in the `v2_mantle_runaway_diagnostic` ignored test); it is not a bug in the mantle or slab implementations individually.

**Trajectory** (20 steps at 64², mantle+slab, baseline parameters):

| steps | peak\|v_solved\| | peak\|f_slab\| | alignment |
|---|---|---|---|
| 5 | `9.6e0` | `9.8e0` | `+0.22` |
| 10 | `3.3e1` | `5.5e1` | `+0.23` |
| 15 | `1.5e7` | `1.0e6` | `−48` |
| 20 | `7.9e14` | `4.0e13` | `−1.9e9` |

**Closed-loop gain analysis (Step 8 regime, bootstrapped).** Once mantle forcing pulls `v ~ O(Mf) = O(1)`, the power-law rheology exits the floor-dominated band: `ε̇_II ~ v/L = O(1)` → `η_newton ≈ ε̇^{1/n−1} ≈ 1`, so the viscous diagonal `2·η·k² ≈ 80` at `k=1` on a 64² grid. In the same regime the discrete divergence operator in `Q_sub_conv = k_slab · max(0,−div v)` amplifies `|div v|_max ≈ 2·|v|/dx = 128·|v|` at grid spacing `dx = 1/64`. Then `m_subducted ≈ Q · τ_slab = 64·v`, and `f_slab = Sp · m ≈ 1.5 · 64 · v = 96·v`. The slab contribution to the momentum balance scales as `96·v` while the viscous dissipation scales as `80·v` — closed-loop gain

```
G_activated = (Sp · k_slab_accum · τ_slab · (2/dx)) / (2·η_op·k²)
            ≈ (1.5 · 1 · 0.5 · 128) / 80
            ≈ 96 / 80
            ≈ 1.2  > 1
```

— linear instability in the activated regime. The §4.8 target band `Sp ∈ [0.5, 3]` was calibrated against quiescent-regime balance assumptions and is **not co-calibrated** with §4.9's `Mf ∈ [0.3, 2]` in the mantle-activated regime.

**This is the second §4.x refutation this milestone.** Step 7 established that slab-pull alone cannot bootstrap out of floor-domination. Step 8 establishes that slab-pull + mantle together in the activated regime produce unbounded positive feedback at the §4.8 baseline parameters. Both findings are revisions of implicit assumptions in `solver-scaling.md`, not implementation bugs.

**Three resolution paths, none selected at this step:**
- **(a) Recalibrate `Sp` in the activated regime.** Stability condition: `Sp · k_slab_accum · τ_slab · (2/dx) / (2·η_op · k²) < 1`. At 64² baseline, this reduces to `Sp < 80/128 ≈ 0.6` — below the §4.8 band's lower edge. A full recalibration would reset the band based on the activated-regime operator.
- **(b) Modify the discrete divergence operator used in `Q_sub_conv`.** The `1/dx` amplification is a discretisation choice; a smoothed or gradient-bounded variant would reduce the gain without altering the §4.8 `Sp` band.
- **(c) Physical saturation of `m_subducted`.** Introduce an upper bound or nonlinear growth law that prevents `m_steady = Q·τ` from scaling linearly with `|div v|` when `|div v|` is already large. Changes slab-pull's contract and is the most invasive path.

**Follow-up issue:** a dedicated slab+mantle co-calibration issue is drafted in `docs/followup_slab_mantle_cocalibration.md` for opening post-Step 8. It does not block Step 9 (cratonic immunity), which can proceed on the mantle-only base.

**Permanent oracle:** the `v2_mantle_runaway_diagnostic` test (currently `#[ignore]`-d) reproduces the runaway with the offending parameter combination. After the co-calibration issue is resolved, that test will be switched to a non-ignored regression guard — any future change that re-introduces the instability will trip it.

### S field evolution

- Var(S̃) timeline: initial `1.633e-1`, middle `6.889e-2`, final `3.119e-2` (Δ = `-80.90%` vs initial)
- max|∇S̃| timeline: initial `8.670e1`, peak `8.670e1`, final `4.684e1`

### Mass conservation of S

- initial mass: `2.281931793e3`
- final mass: `2.251012804e3`
- relative drift: `-1.355e-2`

### Null-space health

- max |mean(vx)| across solves: `1.554e-16`
- max |mean(vy)|: `1.412e-16`

### Velocity magnitude

- peak |v|: `9.528e0`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step8_physics_heightmaps/s_64x64_t0000.png` | `1.601e-1` | `1.198e0` | `5.571e-1` | `docs/reports/step8_physics_heightmaps/s_64x64_t0000_colorbar.png` |
| `docs/reports/step8_physics_heightmaps/s_64x64_t0150.png` | `5.000e-2` | `1.518e0` | `5.481e-1` | `docs/reports/step8_physics_heightmaps/s_64x64_t0150_colorbar.png` |
| `docs/reports/step8_physics_heightmaps/s_64x64_t0300.png` | `5.000e-2` | `1.396e0` | `5.496e-1` | `docs/reports/step8_physics_heightmaps/s_64x64_t0300_colorbar.png` |

### Comparison vs Step 6 physics (advisory — mantle added, slab still off)

The Step 8 physics baseline sits on the Step 6 setup (GPE + yielding + basal drag + Voronoi + Closed recycling) with mantle forcing added on top and slab-pull held Disabled. The mantle contribution bootstraps the system out of floor-domination, so large deltas vs Step 6 are expected in `peak|v|`, `yielding_cell_fraction_max`, strain-rate distribution, and CG iteration counts. This is an advisory comparison only, not a regression test.

#### Grid 64×64 — comparison vs Step 6 physics

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 34.402 | 1145.551 | ×33.30 (>20×, flag) |
| CG iters / linear solve (mean) | 129.6 | 1420.7 | ×10.96 [fail] |
| S mass drift (relative) | -4.391e-4 | -1.355e-2 | gate 1e-10 |
| max \|mean(vx)\| | 3.220e-22 | 1.554e-16 | bruit machine |
| max \|mean(vy)\| | 5.954e-21 | 1.412e-16 | bruit machine |

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
| slab-pull | Disabled |
| mantle forcing | Enabled (Mf = 1.000, coupling = 1.000, num_modes = 6, seed = 42, evolution_rate = 0.000) |
| seed | 42 |

### Timing

- wallclock total: `6704.293 s`
- wallclock per step (mean): `22347.643 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `3.76e4`
- CG iterations per Newton step — mean: `1853.5`, max: `2000`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 520 | 99 |
  | 890 | 173 |
  | 1260 | 219 |
  | 1630 | 187 |
  | 2000 | 4256 |

### Newton (nonlinear) health

- outcome distribution — Converged: `99.0%`, Stalled: `0.3%`, Diverged: `0.0%`, CappedIters: `0.7%`
- Newton outer iters per timestep — mean: `16.2`, max: `20`
- effective η_max/η_min over run — mean: `139449.30`, max: `265505.28`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%`
- continuation ramp: ✅ all 5 sub-solves converged

### Plastic yielding

- Bi = `0.150`
- yielding_cell_fraction (max over run, criterion `η_eff < 0.5·η_visc`): `0.985`
- yielding_intensity (mean of `η_visc/η_eff − 1` where `η_eff < 0.9·η_visc`, max over run): `25.900`
- Definition notes: the `< 0.5·η_visc` criterion captures "yielding dominant", not "yielding present anywhere" (the legacy `η_p < η_v` metric saturated near 1.0 — cf. issue #75).

### Strain-rate regime diagnostic (floor-domination)

- `ε̇_min` (regularisation floor): `1.000e-3`
- `ε̇_II` at final timestep: mean = `1.631e1`, max = `1.360e2`
- Fraction of cells with `ε̇_II < 10·ε̇_min = 1.000e-2` at final timestep: `0.058` (1.0 = everywhere in the floor-dominated band)
- `max(ε̇_II) / ε̇_min` = `135971.73` — ratio of the strongest strain-rate cell to the regularisation floor.

**Verdict:** partial floor-domination — `94.2%` of cells are above the `10·ε̇_min` threshold. `yielding_cell_fraction` should be roughly consistent with that active fraction.

### Basal drag

- Br = `0.050`
- `basal_drag_energy_ratio` (mean over run of `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`): `3.710e-4`
- `drag_vs_visc_diagonal_ratio` (mean over run of `mean_cells(Br·S̃² / (η/Δx²))`): `3.716e-4`
- Algebraic identity check `r/(1+r) ≈ energy_ratio`: predicted `3.715e-4` vs measured `3.710e-4` (relative diff `1.3e-3`; spec bound: coarse, typically `< 1e-1`)
- `peak_v_damping_ratio`: — (requires regression run; use `--forcing both` or `--forcing sinusoidal`)

**Expected magnitude of drag effect at Step 4 (corrected vs spec).** The baseline has `S̃ ≈ 1` uniformly and `Br = 0.05`, so `Br · S̃² ≈ 0.05` per cell. The viscous diagonal per cell is `η · N²` (approximately — see `momentum_diagonal` for the exact stencil): at 64² that's `N² = 4096`, at 128² it's `16 384`. The Step-4 spec's sketched band `[10⁻⁶, 10⁻⁴]` assumed `η ≈ 1`, but the power-law rheology at this baseline is floor-dominated — `ε̇_II` lies below `ε̇_min = 10⁻³` everywhere, so `η_newton = ε̇_min^(1/n-1) = (10⁻³)^{-2/3} ≈ 100` in the bulk, which the soft cap (η_max = 10³) barely attenuates. With the corrected `η ≈ 100`, the drag/viscous ratio sits in `[10⁻⁸, 10⁻⁷]`: ≈ `1.2×10⁻⁷` at 64² and ≈ `3×10⁻⁸` at 128². Both measured values above fall in this corrected band — the smallness is **by construction of the Step 4 baseline** (no oceanic cells yet; those arrive at Step 5/6 with `S̃ ≈ 0.2` so `S̃² ≈ 0.04` creates ×25 differentiation between continental and oceanic drag; Step 9 will raise the cratonic η). Step 4 installs the machinery; its full physical effect shows up later.

**Yielding checkpoint.** Basal drag is dissipative — it removes kinetic energy rather than injecting strain. `yielding_cell_fraction` is expected to stay at 0 at the Step 4 baseline (yielding is Disabled here and would anyway remain floor-dominated under Br alone). The yielding activation threshold will be re-checked at Step 5 (boundary sources inject mass) and Step 7 (slab pull operates at τ*/Sp ≈ 10–60 Myr), not at this step.

### Boundary source/sink diagnostics

- `s_oceanic_mean` = `0.4837` (std `0.2814`)  — target `[0.18, 0.22]` post-calibration
- `s_continental_interior_mean` = `0.6059` (std `0.2309`)  — target `[0.9, 1.1]`
- `s_continental_collision_mean` = `0.6656` — telemetry only (orogen thickening, tracked through Steps 5-10)
- `boundary_type_diversity` = `4` (number of distinct mechanisms active on the run)
- `clamp_activation_fraction` — mean `3.490e-3`, max `5.737e-3` (healthy: mean < 1%, max < 5%)
- `∫Q dt dA` = `-5.894e-3`; `∫clamp_flux dt dA` = `4.492e-3`
- `mass_balance_residual` = `1.628e-13` (issue #89 D5; acceptance `< 1%`)

### Voronoi plate geometry

- distinct plate_count = `8` (expected 8 for `num_plates=8`)
- plate_type_distribution (oceanic, continental) = `(0.555, 0.445)` — target continental ∈ [0.15, 0.45]

### Boundary dynamics (dynamic detection per step)

- `boundary_flag_transition_rate` — mean `3.950e-3`, max `2.122e-1`
  - Fraction of cells whose `boundary_flag` changed vs the previous step. Telemetry only — no acceptance. Expected transient spike early in the run (flags emerging from `None` as the first Stokes solves produce non-trivial divergence), then stabilisation.
- flag counts **at step 1** (proving detection fired): None=`59`, Subduction=`105`, OceanicSubduction=`5446`, Rift=`7183`, ContinentalCollision=`3591`
- flag counts **at final step**: None=`89`, Subduction=`66`, OceanicSubduction=`4905`, Rift=`8203`, ContinentalCollision=`3121`

### Recycling health (Closed mode)

- `recycling_buffer_fill` — mean `5.188e-3`, max `5.926e-3`, final `5.894e-3`
- `immediate_pending_max` over run = `0.000e0`, final sum = `0.000e0`
- `clamp_activation_during_spinup_max` = `0.000e0` (target 0 — clamp should not fire during the buffer fill-up)

### Mass balance (Step 6 closed recycling, 5 components)

- Δmass_observed (dimensionless, S̃ sum): initial `9.124226e3`, final `9.101252e3`, Δ = `-2.297e1`
- `buffer_fill_final` (cell-area units) = `5.894e-3`
- `pending_immediate_final` (cell-area units) = `0.000e0`
- `clamp_flux_integral` (cell-area units) = `4.492e-3`
- `mantle_loss_integral` (cell-area units) = `0.000e0` (zero when mantle_loss_fraction=0)
- **`mass_conservation_residual` = `1.508e-15`** (target `< 1e-6`)

Formula: `|Δmass_obs + mantle_loss + buffer_fill + pending − clamp_flux| / initial_mass`. All five components are tracked; the residual is the absolute sum divided by `initial_mass`. A `< 1e-6` residual means the pipeline is mass-exact at machine precision; all deviations from exact conservation are accounted for by the known components (loss + in-transit buffer mass + rollover pending + clamp artificial flux).

### Issue #78 trajectory (5 instants: t ∈ {1, 10, 50, 150, 300}·Δt)

| step | max\|∇S̃\|_interface | max\|∇S̃\|_global | peak\|f_GPE\|_interface | peak\|f_GPE\|_global | buffer_fill |
|---|---|---|---|---|---|
| `1` | `1.716e2` | `1.716e2` | `1.229e1` | `1.229e1` | `6.467e-5` |
| `10` | `1.694e2` | `1.694e2` | `1.279e1` | `1.279e1` | `1.720e-3` |
| `50` | `1.757e2` | `1.757e2` | `1.177e1` | `1.177e1` | `4.440e-3` |
| `150` | `1.499e2` | `1.499e2` | `1.467e1` | `1.467e1` | `5.912e-3` |
| `300` | `1.113e2` | `1.113e2` | `1.188e1` | `1.188e1` | `5.894e-3` |

**Interpretation.** No taper was applied at the Voronoi oceanic/continental interfaces (per Step 6 D5 — #78 is tested, not contoured). A spike that appears at step 1 and damps by step 50 is a transient artefact of the raw contrast; a spike that grows monotonically across the 5 instants is a real signal that #78 has activated and must be addressed before Step 7. **Absolute critical threshold**: `peak|f_GPE| > 100` at any instant = red-flag bug.

### Continental mass balance (Closed mode)

Continental cells cannot drain via Q_sub (Step 5 invariant: Q_sub fires only on `(Oceanic, is_subduction())` cells). Continental thickness changes come from three sources: (1) **immediate recycling returns** (`Q_arc + Q_coll_v + Q_rift_v`, all applied to continental eligible cells), (2) **advection** across the continental/oceanic boundary, driven by GPE spreading, and (3) **no other Q contribution**.

- `M_sub_total` (integrated drain, all oceanic subducting cells): `1.059e-1`
- `∫Q_arc dt dA` (continental return, arc volcanism): `1.588e-2` — fraction `0.150` of M_sub
- `∫Q_coll_v dt dA` (continental return, collision volcanism): `3.177e-3` — fraction `0.030`
- `∫Q_rift_v dt dA` (continental return, rift volcanism): `2.118e-3` — fraction `0.020`
- Total continental return: `2.118e-2` — fraction `0.200` of M_sub
- `∫Q_spread dt dA` (oceanic return, mid-ocean ridges): `7.882e-2` — fraction `0.744` of M_sub

`s_continental_interior_mean = 0.6059` at end of run (target `[0.9, 1.1]`).

**Interpretation** — with default fractions `(arc 0.15, coll_v 0.03, rift_v 0.02, spread 0.80)` the immediate continental return is **20% of M_sub** while 80% is routed through the delayed buffer to OCEANIC ridges. Net continental balance depends on (a) how much mass the Voronoi advection pushes across the continental/oceanic boundary, and (b) how evenly the 20% immediate return is distributed over the continental cell population.

If `s_continental_interior_mean < 0.9`, the interpretation is that the **continental set is a net mass exporter** to the oceanic set via advection — GPE drives flow away from high-S continental cells toward the thinner oceanic strip, and only 20% of the subducted mass returns to continental via arc + collision + rift volcanism. Global mass is conserved (the spread_fraction=0.80 returns to oceanic cells via the delayed buffer), but the continental/oceanic **partition** is not invariant.

This is expected physics, not a bug. The `[0.9, 1.1]` target band from issue #90 was set against the Step 5 static layout (where continental cells sat in spatial isolation from subduction). With a Voronoi tessellation where continental patches are surrounded by advecting oceanic zones, mass redistribution over 300 steps is larger — the continental mean drifts toward a new Voronoi-specific equilibrium that is not 1.0. Adjusting the acceptance band to reflect Voronoi dynamics is follow-up work; the mass budget itself (`mass_conservation_residual < 1e-6`) holds unambiguously.

### Note on OceanicSubduction drain symmetry

When two oceanic cells meet at a convergent boundary, both are flagged `OceanicSubduction` and both contribute to `Q_sub`. This effectively doubles the local drain compared to Oceanic/Continental subduction (where only the oceanic cell drains). This is an assumed approximation in the absence of an age field (Step 10) that would resolve which cell actually subducts. The mass budget stays correct because the combined drain feeds the same recycling pool: total mass conservation is satisfied independently of which side is drained. To be refined at Step 10.

### Yielding activation checkpoint (Step 6)

- Bi = `0.150`, `yielding_cell_fraction_max` = `0.985`

**Checkpoint status: ✅ activated at Step 6.** Dynamic boundary geometry + closed recycling produced enough convergent strain at some cells to push `ε̇_II > ε̇_min` locally, crossing the Bi threshold. The mechanism is wired and active; expect further growth at Steps 7 (slab pull) and 8 (mantle forcing).

### Preconditioner surveillance (continued from Step 5)

Step 5 physics: CG mean = 108.5 (64²) / 205.0 (128²), ≈ 2× Step 4. Step 6 adds Voronoi interfaces (sharper contrasts, more heterogeneity). If the CG ratio vs Step 5 is ≤ 2× (i.e., vs Step 4 ≤ 4×), continue surveillance. If > 10× Step 4, the preconditioner has reached its usable limit and the maintenance task (block-Jacobi / ILU(0)) should be scheduled before Step 7.

### Mantle bootstrap (Step 8)

- Mf = `1.000` (target band [0.3, 2.0] per §4.9)
- coupling = `1.000` (target band [0.1, 10.0])
- num_modes = `6`, seed = `42`

- `peak|v_mantle|` (= Mf · peak|v_pattern|) = `1.661e1`
- `peak|v_solved|` (max over run) = `9.622e0`
- `v_solved_to_v_mantle_alignment` (mean of `<v, Mf·v_m>/|Mf·v_m|²`) = `0.235`
- `div_v_mantle_max` = `2.274e-13` (strict acceptance `< 1e-10`)

**Bootstrap: ✅ system escaped floor-domination.** `peak|v_solved|` exceeds 0.1 — three or more orders of magnitude above the Step 7 baseline (3.6e-5). Mantle forcing is performing its role as the mechanism-hierarchy initiator (see §4.8 activation-regime note).

### Force hierarchy (Step 8)

- `peak|f_GPE|` = `1.637e1`
- `peak|f_slab|` = `0.000e0`
- `peak|f_mantle|` = `2.034e1`
- `f_mantle / f_GPE` (mean per step) = `1.350e0`

**Interpretation bands** (telemetry, not acceptance — except the pathological case):
- `f_mantle ≫ f_GPE` (ratio ≥ 10): mantle bootstrapped. Success.
- `f_mantle ~ f_slab` (ratio 0.1–10): healthy coupling.
- `f_slab ≫ f_mantle` (ratio < 0.1): non-pathological, document.
- `f_mantle ≪ f_GPE` (ratio < 0.1): PATHOLOGICAL — correlates with bootstrap failed, remontée required.

### Yielding activation (Step 8 — STRICT, last chance)

- Bi = `0.150`, `yielding_cell_fraction_max` = `9.849e-1`
- `max(ε̇_II) / ε̇_min` = `3.140e5` (floor-dominated if ≤ 1)

**Yielding activation: ✅ RESOLVED.** The checkpoint transported since Step 3 (and strictly enforced here as last-chance per the Step 7 revision) is met: yielding fires in a non-marginal fraction of cells. Mantle forcing has bootstrapped `ε̇_II` above the regularisation floor locally, and the Bingham criterion (`η_eff < 0.5 · η_visc`) captures the resulting yielding-dominated regime. The mechanism hierarchy is confirmed.

### Slab+Mantle interaction instability finding (Step 8)

The Step 8 baseline above holds **slab-pull Disabled** by deliberate choice. During Step 8 development, running the nominal spec configuration (Step 7 physics + mantle Enabled) produced catastrophic numerical divergence within 15–20 timesteps at 64² × Mf=1.0, `coupling=1.0`, slab-pull at Step 7's `(Sp=1.5, τ_slab=0.5, k_slab_accum=1.0)`. The runaway is physically real (captured in the `v2_mantle_runaway_diagnostic` ignored test); it is not a bug in the mantle or slab implementations individually.

**Trajectory** (20 steps at 64², mantle+slab, baseline parameters):

| steps | peak\|v_solved\| | peak\|f_slab\| | alignment |
|---|---|---|---|
| 5 | `9.6e0` | `9.8e0` | `+0.22` |
| 10 | `3.3e1` | `5.5e1` | `+0.23` |
| 15 | `1.5e7` | `1.0e6` | `−48` |
| 20 | `7.9e14` | `4.0e13` | `−1.9e9` |

**Closed-loop gain analysis (Step 8 regime, bootstrapped).** Once mantle forcing pulls `v ~ O(Mf) = O(1)`, the power-law rheology exits the floor-dominated band: `ε̇_II ~ v/L = O(1)` → `η_newton ≈ ε̇^{1/n−1} ≈ 1`, so the viscous diagonal `2·η·k² ≈ 80` at `k=1` on a 64² grid. In the same regime the discrete divergence operator in `Q_sub_conv = k_slab · max(0,−div v)` amplifies `|div v|_max ≈ 2·|v|/dx = 128·|v|` at grid spacing `dx = 1/64`. Then `m_subducted ≈ Q · τ_slab = 64·v`, and `f_slab = Sp · m ≈ 1.5 · 64 · v = 96·v`. The slab contribution to the momentum balance scales as `96·v` while the viscous dissipation scales as `80·v` — closed-loop gain

```
G_activated = (Sp · k_slab_accum · τ_slab · (2/dx)) / (2·η_op·k²)
            ≈ (1.5 · 1 · 0.5 · 128) / 80
            ≈ 96 / 80
            ≈ 1.2  > 1
```

— linear instability in the activated regime. The §4.8 target band `Sp ∈ [0.5, 3]` was calibrated against quiescent-regime balance assumptions and is **not co-calibrated** with §4.9's `Mf ∈ [0.3, 2]` in the mantle-activated regime.

**This is the second §4.x refutation this milestone.** Step 7 established that slab-pull alone cannot bootstrap out of floor-domination. Step 8 establishes that slab-pull + mantle together in the activated regime produce unbounded positive feedback at the §4.8 baseline parameters. Both findings are revisions of implicit assumptions in `solver-scaling.md`, not implementation bugs.

**Three resolution paths, none selected at this step:**
- **(a) Recalibrate `Sp` in the activated regime.** Stability condition: `Sp · k_slab_accum · τ_slab · (2/dx) / (2·η_op · k²) < 1`. At 64² baseline, this reduces to `Sp < 80/128 ≈ 0.6` — below the §4.8 band's lower edge. A full recalibration would reset the band based on the activated-regime operator.
- **(b) Modify the discrete divergence operator used in `Q_sub_conv`.** The `1/dx` amplification is a discretisation choice; a smoothed or gradient-bounded variant would reduce the gain without altering the §4.8 `Sp` band.
- **(c) Physical saturation of `m_subducted`.** Introduce an upper bound or nonlinear growth law that prevents `m_steady = Q·τ` from scaling linearly with `|div v|` when `|div v|` is already large. Changes slab-pull's contract and is the most invasive path.

**Follow-up issue:** a dedicated slab+mantle co-calibration issue is drafted in `docs/followup_slab_mantle_cocalibration.md` for opening post-Step 8. It does not block Step 9 (cratonic immunity), which can proceed on the mantle-only base.

**Permanent oracle:** the `v2_mantle_runaway_diagnostic` test (currently `#[ignore]`-d) reproduces the runaway with the offending parameter combination. After the co-calibration issue is resolved, that test will be switched to a non-ignored regression guard — any future change that re-introduces the instability will trip it.

### S field evolution

- Var(S̃) timeline: initial `1.633e-1`, middle `1.192e-1`, final `7.745e-2` (Δ = `-52.59%` vs initial)
- max|∇S̃| timeline: initial `1.735e2`, peak `2.005e2`, final `1.113e2`

### Mass conservation of S

- initial mass: `9.124226001e3`
- final mass: `9.101252054e3`
- relative drift: `-2.518e-3`

### Null-space health

- max |mean(vx)| across solves: `1.819e-16`
- max |mean(vy)|: `3.010e-16`

### Velocity magnitude

- peak |v|: `9.586e0`

### Heightmaps of S (dynamic remap with bounds)

| snapshot | min | max | mean | colour-bar |
|---|---|---|---|---|
| `docs/reports/step8_physics_heightmaps/s_128x128_t0000.png` | `1.600e-1` | `1.198e0` | `5.569e-1` | `docs/reports/step8_physics_heightmaps/s_128x128_t0000_colorbar.png` |
| `docs/reports/step8_physics_heightmaps/s_128x128_t0150.png` | `5.000e-2` | `1.665e0` | `5.535e-1` | `docs/reports/step8_physics_heightmaps/s_128x128_t0150_colorbar.png` |
| `docs/reports/step8_physics_heightmaps/s_128x128_t0300.png` | `5.000e-2` | `1.606e0` | `5.555e-1` | `docs/reports/step8_physics_heightmaps/s_128x128_t0300_colorbar.png` |

### Comparison vs Step 6 physics (advisory — mantle added, slab still off)

The Step 8 physics baseline sits on the Step 6 setup (GPE + yielding + basal drag + Voronoi + Closed recycling) with mantle forcing added on top and slab-pull held Disabled. The mantle contribution bootstraps the system out of floor-domination, so large deltas vs Step 6 are expected in `peak|v|`, `yielding_cell_fraction_max`, strain-rate distribution, and CG iteration counts. This is an advisory comparison only, not a regression test.

#### Grid 128×128 — comparison vs Step 6 physics

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 312.293 | 6704.293 | ×21.47 (>20×, flag) |
| CG iters / linear solve (mean) | 240.4 | 1853.5 | ×7.71 [suspect] — suspect (no justification on file) |
| S mass drift (relative) | -4.394e-4 | -2.518e-3 | gate 1e-10 |
| max \|mean(vx)\| | 4.205e-22 | 1.819e-16 | bruit machine |
| max \|mean(vy)\| | 7.775e-21 | 3.010e-16 | bruit machine |

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
