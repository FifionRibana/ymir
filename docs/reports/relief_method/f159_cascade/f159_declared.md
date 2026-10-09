# F159 — declared before any measurement (2026-10-09)

The multi-scale cascade prototype: 64² → 128² → 256² → 512² on the témoin. Each level gets tectonic uplift, erosion and
diffusion, and each level is visible in the viz. **Production is unchanged.** The cascade is a separate chain, gated
(nothing calls it except the viz's Cascade window and the `f159_*` benches), and read-only on the tectonic state.

Written after reading:
- the solver (`stream_power.rs`);
- the production chain (`production_upscale.rs`, `upscale.rs`);
- Perron 2009, re-read for this document. L_c = (D/K)^(1/(2m+1)), eq. 3. The first-order valley spacing is
  6.4 L_c ≤ λ ≤ 12.7 L_c (Fig. 2). « L_c² is approximately the drainage area at which the topography transitions
  from a concave-down, creep-dominated hillslope to a concave-up, stream-incision-dominated valley ».

Nothing has been run on the cascade.

## D1 — the levels

| level | cells | cell (km) | cell area (km²) | role |
|---|---|---|---|---|
| 64² | 4 096 | 6.25 | 39.1 | the C1 tectonics, **unchanged**: `c1_coarse_normalized_altitude(state, IsostasyConfig::c1_default(), ss, None)` (`target_land_fraction` = production's `None`) |
| 128² | 16 384 | 3.125 | 9.77 | cascade level 1 |
| 256² | 65 536 | 1.5625 | 2.44 | cascade level 2 |
| 512² | 262 144 | 0.78125 | 0.610 | cascade level 3 |

**The framing is the témoin's**: the 400 km domain and the roll (6, 37) of 64² (`CANONICAL_ORIGIN`).
- Pixel i of a level of N cells sits at coarse coordinate 6 + i·64 / N, as in `upscale_with_fbm` (`sx = origin +
  i·scale`).
- So the even pixels of level 2N coincide with level N's pixels.
- The domain is the torus, read as a bounded map, as in production: the incision's routing ends at the map edge.

## D2 — the uplift U(x)

**U(x) = U₀ · max(0, h_iso(x)) / 1 000 m**, in m/yr. h_iso is the coarse C1 altitude in metres: (norm − 0.5) · 11
300, the `c1_altitude_norm_to_metres` contract.

**What it derives from**: h_iso is C1's isostatic altitude (`c1_production_altitude_craton`). It integrates:
- the crustal thickness s̃, which carries the convergent thickening;
- the thermal age, ocean side;
- the cratonic buoyancy (the craton mask).

It does **not** read the plate boundaries directly.

**The reading**: the isostatic altitude is the height the tectonics holds up. The uplift needed to hold it against
erosion is taken proportional to it. This is the author's « la tectonique actuelle est la base du soulèvement ».

**The sea**: cells with h_iso ≤ 0 get no uplift.

**Constant in time this round.** Baldwin 2003 and Braun & Robert 2005 (post-orogenic decay) are NOT applied.
- A constant U has a steady state, and the steady state is the stop criterion (D7).
- A decaying U has none, so the duration would again be a proxy (F44).
- The decay is a later round's question.

**U₀ is calibrated once, at 128²**:
- **The target**: the steady state's mean land altitude equals h_iso's, both at 128² (h_iso bicubic-upscaled).
- **The method**:
  1. one trial run at U₀ = 1e-4 m/yr to equilibrium;
  2. then U₀ = 1e-4 · target / measured.
- **Why one rescale suffices**: for n = 1 with linear diffusion, the steady state is linear in U.
- **The caveat**: the talus is not linear. Its activity at 128² is reported; if it acts on more than 1 % of the land
  the linearity is broken, and this is said.
- **Cost**: the trial run's is reported inside the cascade's cost.

**Upscaling of U**: bicubic ×2 from level to level, like z (D5).

## D3 — the processes per level

Each level is a loop of steps, with the operators split in production's order:
1. **uplift**: z += U · dt on land cells (z > sea);
2. **erosion**: the existing implicit stream power, `incise_with_floor` with `iterations = 1`. It is `relief_v3` at
   the level's cell:
   - MFD p = 2 for the area;
   - D8 receivers and stack on the filled surface;
   - the base-level floor at 0.5 m;
   - lateral erosion K_lat = 4.0 (physical half-width in metres);
   - A_c = 0.1 km² (physical).
   - **Inside this call the talus and the diffusion are switched off**; they run in step 3.
3. **diffusion**, the hillslope closures, in production's order:
   - the talus (tan 33°, 4 passes, factor 0.5);
   - then the linear diffusion, everywhere on land (`diffuse_channels = true`).
   - **Both are the SAME code** as in `incise_with_floor`, extracted unchanged into two functions that
     `incise_with_floor` now calls. A pure refactor; the guard 6 / 6 checks it.

**Declared absences**:
- No deposition.
- No basin graph: depressions are filled by the priority flood at each step, as now. A1 is off.
- No lithology / fracture K field (uniform K).
- No FBM.
- No edifices.
- No coast warp. So the cascade's coastline is the bicubic of the 64² coast, not the témoin's warped coast; the
  instruments read each field on its own land.
- No bathymetry: the sea cells get neither uplift nor erosion.
  - **During a level they are held at sea level**, the LEM's base level. Otherwise the diffusion would pull the coast
    down toward a −2 000 m neighbour, and the talus could shed onto the sea.
  - **The land mask is fixed for the level.**
  - The upscaled bathymetry is restored at the end of the level.

## D4 — the resolution dependence (the central point)

**K is physical and the same at every level.** K = 2.0e-3 in the code's units (E [m/yr] = K · A_km²^0.5 · S), i.e.
K_SI = 2e-6 yr⁻¹ for m = 0.5, n = 1.
- With n = 1 the steady state depends on U / K only, so K sets the clock and U₀ sets the height.
- **Kwang & Parker 2017**: with m / n = 0.5 and no diffusion, the stream power law is scale-free. The valleys then
  collapse to the grid. **So K alone cannot carry a length scale, and scaling K between levels is not the fix.**

**The length scale comes from D, through Perron 2009**:
- L_c = √(D / K_SI) (m = 0.5);
- the first-order valley spacing is 6.4 L_c ≤ λ ≤ 12.7 L_c.

**D per level**:

  **D_L = max(D_phys, α · K_SI · cell_L²)**, with α = 0.16.

- **The sub-grid term** α K cell² gives L_c = 0.4 cell, so λ = 2.6–5.1 cells. These are the finest valleys a level
  can carry; Nyquist is 2 cells.
  - Each level adds the valleys its resolution can carry: the LOD principle. The coarse level's diffusion stands for
    the dissection it cannot resolve.
  - This is the sub-grid closure Armitage 2019 and Perron 2008 point at, done here by the diffusivity, not by the
    routing (declared: the routing remains MFD p = 2 on the grid).
- **The physical term** D_phys = K_SI · A_c = 2e-6 × 1e5 m² = **0.2 m²/yr**. It ties D to the code's existing
  channel head A_c = 0.1 km², by Perron's L_c² ≈ A_c, so L_c,phys = 316 m.
  - **PROXY**: it inherits A_c's calibration (F1), not a measured D.
  - **Its consequence for later rounds, stated now**: from 1 024² down, D_phys governs. The finest valleys of the
    final world are then λ = 2.0–4.0 km.
- **At the prototype's levels**:

  | level | α K cell² (m²/yr) | D_phys | D_L | L_c (cells) | λ Perron (cells) |
  |---|---|---|---|---|---|
  | 128² | 3.13 | 0.2 | 3.13 | 0.40 | 2.6–5.1 |
  | 256² | 0.78 | 0.2 | 0.78 | 0.40 | 2.6–5.1 |
  | 512² | 0.195 | 0.2 | **0.2** | 0.405 | 2.6–5.2 |

  The two terms cross at 512² (cell 790 m).

**How the code reads it**: the existing linear diffusion's weight is dimensionless per step: κΔt / dx² (F45; F158
corrected). It becomes w_L = D_L · dt / cell_L² = 0.032 at every level (4 sub-steps of 0.008; explicit, stable
≤ 0.25).
- This answers F45's dimensional debt **for the cascade only**: the weight is derived from a D in m²/yr.
- Production's 0.08 is untouched.

**The time step**: dt = 1e5 yr at every level. The scheme is implicit, so dt only sets the clock:
- for n = 1 the split step's fixed point satisfies U = K A^m S (dt cancels);
- the diffusion's weight is linear in dt.

**The head threshold**: A_c = 0.1 km², physical, from production. It is below one cell at 128²–512², so every land
cell is a channel cell at these levels. The hillslope is then D_L's sub-grid closure. Reported, not hidden.

**The talus**: a physical slope. tan 33° is a drop in metres per metre of run (`base_drop = S_c · cell_m`), so it is
resolution-invariant by definition. Its activity is reported per level.

**Lateral erosion**: a physical half-width K_lat · A^0.5 in metres. It is sub-cell at these levels for A below about
38 000 km² (at 512²), so expected inert. Reported.

## D5 — the upscaling

- **The method**: periodic Catmull-Rom bicubic ×2, on the cell-corner mapping of D1. Even pixels are the coarse values
  exactly; odd pixels are interpolated.
- **What is passed on**: z and U. Nothing else: no network, no lakes, no K field, no sea mask. A cell is sea at the
  new level if its bicubic z ≤ sea.

## D6 — keeping the coarse level

**Measure only this round, no constraint.** Reported per level:
- level n+1 brought back to level n, against level n: the mean bias and mean |Δ| (m), on land.
  - The restriction is the centred full weighting (1/4, 1/2, 1/4, separable), because level n's pixel I sits on n+1's
    pixel 2I.
  - A plain 2×2 block would be off by half a cell. Corrected before any measurement.
- the trunk holding (D8).

## D7 — the stop per level

- **Equilibrium**: the mean |Δz| per step over land < **1 %** of the mean uplift per step (mean U · dt over land), for
  5 consecutive steps.
- **Caps**: 400 steps (128²), 300 (256²), 200 (512²).
- **The calibration run** (128²) uses the same rule.
- A level that stops at its cap is reported as « not at equilibrium », with its last ratio.

## D8 — the instruments (per level, cascade and reference alike)

Slopes are central differences in metres, as in F157-B7.

1. **Planar walls (F155)**: % of land at 28° ± 0.5°.
   - **Added**: % of land at 33° ± 0.5°, the cascade's talus angle. A talus could make its own planar band, and
     measuring only the construction's angle would miss it.
2. **Sharp crests**:
   - **crest cells**: land cells higher than both neighbours along at least one of the four axes (x, y, two diagonals);
   - **reported**: their share of land, and p50 / p90 of −∇²z at crests (5-point Laplacian, m / km²).
3. **Crest facets** (new) and **their control**:
   - **The definition**: a crest cell is a facet crest when the median |∇²z| of the non-crest land cells of its ring
     at Chebyshev distance 2 is < 0.25 × its own |∇²z|. That is a sharp edge between flanks with no curvature, the
     author's « crêtes en facettes ». The facet share = facet crests / crest cells.
   - **The control** (rule 13), synthetic, at each level's grid:
     - a field of pyramids (planar faces): facet share ≥ 60 % required;
     - a field of paraboloid hills (uniform curvature, a diffusive hilltop): ≤ 20 % required.
     - If the control fails at a level, the instrument is declared blind there and cannot feed the stop rule.
   - **For information**: S1 and the construction (8192², block-averaged as the reference).
4. **Comb teeth**: F124's teeth need the construction's wall mask, so **they do not apply to the cascade**. The
   declared stand-in is F84's coherent-window share (`aniso`, windows of 8 cells with coherence > 0.7): parallel
   valleys or ridges. It is a stand-in, not F124's number.
5. **Axis alignment (F126)**: the terrain R8 (`aniso`, windows of 8 cells at every level, land-only windows). The
   noise floor ≈ 1 / √(windows) is printed beside it.
6. **Trunk holding**:
   - **trunks**: land cells with D8 area ≥ 100 km² (`compute_flow` on the level's field, sea 0.5);
   - **per transition n → n+1**, the distance from each trunk cell of n+1 to the nearest trunk cell of n:
     p50 / p90 / max, in km and in cells of level n. Also the reverse direction.
   - **Also**: the cascade's trunks against the reference's at the same level.
7. **Valley spacing against Perron**:
   - **the method**: along every row and column, in land runs ≥ 16 cells, count the local minima of z with a
     prominence ≥ 1 m within the run. λ = run length / minima.
   - **the comparison**: against 6.4–12.7 L_c in cells. A transect crosses valleys obliquely, so it reads λ high by
     an orientation factor (≈ π/2 for an isotropic network). Declared, not corrected.
8. **Relief**:
   - local relief: max − min over non-overlapping 12.5 km all-land blocks (4 / 8 / 16 cells), p50;
   - the mean land altitude and the land p90.
9. **Cost per level**:
   - the wall time per sub-step (upscale, uplift, erosion, diffusion) and the total;
   - the steps and the calibration;
   - the memory, estimated from the retained buffers (cells × grids × 4 bytes). Declared an estimate, not a
     measured peak.
   - **The extrapolation to 8192² (ESTIMATE)**: T = Σ over N = 1 024 … 8 192 of t_step(512) · (N² log N²) /
     (512² log 512²) · steps(512), plus the measured 128²–512². It assumes each finer level needs as many steps as
     512². Compared with 15 min and 1 h.

## D9 — the reference

**The production témoin**: C2 /10 col (`Knobs { valley: ValleyConstruction::new(F121_AGE_K, Some(0.1)),
slope_floor_abs: 0.024, ..Knobs::passes(2) }`, PSEED, CANONICAL_ORIGIN), built at 8192².
- **Brought to each level by centred block means**: for level pixel I, the témoin pixels f·I − f/2 … f·I + f/2 − 1,
  wrapped, with f = 64 / 32 / 16. Centred so that the block and the level's pixel share a centre.
- **The caveat**: averaging lowers the slopes and curvatures. The reference at a level is « what the témoin carries
  at that resolution ».

## R — the stop rule (the brief's, verbatim, plus how it reads)

« s'arrêter au niveau où l'un est vrai :
- murs plans ou facettes de crête au-dessus du témoin au même niveau ;
- un tronc déplacé de plus de 2 cellules du niveau précédent (p90) ;
- le coût extrapolé à 8192² dépasse 1 h. »

**How it reads, declared**:
- **« au-dessus »**: strictly greater share, for either wall instrument (28° or 33°) or the crest-facet share. The
  facets count only where their control passes.
- **« un tronc déplacé »**: the p90 of the n+1 → n trunk distance, in cells of level n.
- **The cost**: the D8.9 ESTIMATE.

## I — the images (declared crops)

Three 100 km windows, chosen by rule on the témoin's 512² block mean, coordinates printed by the bench:
- **A, « chaîne »** (range): the window with the highest local relief p50.
- **B, « plateau »**: among windows ≥ 90 % land with mean altitude above the land p75, the lowest local relief.
- **C, « bassin + lac »**: the window centred on F157's lake 4, (1913, 3277) at 8192².

**What is shown**:
- per level, and for the last level's sub-steps: Ombrage (the viz's hillshade), north up, the same window;
- beside each: the témoin at the same level (block mean) and at 8192².
