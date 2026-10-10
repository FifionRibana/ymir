# F165 — declared before any measurement of our worlds (2026-10-10)

**The continent's physiography: the physics level follows C1's history.**

**Production is unchanged**: the guard 6 / 6 (field and lakes, a8d2d538d692c2f0) runs before and after. C1's outputs
stay bit-identical (the history is recorded read-only, under a permanent test). The cascade is frozen at F164's
setting: p6 · W-tous W1, R4+π, F163's roughness, the budgets and k of F164's best chain (512² 5.623e2,
1 024² 1.000e3, 2 048² 3.162e2), with no recalibration downstream.

**Not blind**:
- F159–F164's readings;
- the code inventory of A1 below, read before writing this (for example, the rigid continental crust);
- the témoin's physics level at 256² (λ 17.6 km and the « filins » the author saw).

Nothing about the four seeds' C1 fields, the history, Europe or the classifier on any terrain has been measured.

## A1 — the inventory (read from the code, before any measurement)

- **C1's time loop** (`tectonics_c1/time_loop.rs::run_with_closures_observed`):
  - **300 steps** (`C1TimeLoopConfig::n_steps`, the benches' and the viz's value);
  - Δt = 0.5 · dx / max|v| in non-dimensional units, ≈ 0.69 per step at 64²;
  - `docs/c1_lightweight_dynamic_tectonics.md` §11 maps it to ~30–100 ka per step, ~10–30 Ma for 300 steps.
    `SteinSteinParams::age_to_ma = 0.667` maps the same 300 steps to ~200 Ma for the ocean ages;
  - **C1 has no pinned physical duration.** Hence T is the calibrated constant (C2).
- **Where the fields change**:
  - `s` is advected (step 1), then:
    - Davis-Suppe's orogenic source on upper-plate wedge cells (2);
    - the equilibrium-height sink (3);
    - the stream-power sink (4);
    - subduction (5a);
    - the rift thinning (5b);
  - `age` is advected and reset at rift splits;
  - `plate_type` changes at subduction's floor trigger (oceanic → continental);
  - `plate_id` changes at subduction, accretion merges and rift splits;
  - `cratonic_mask` is built at init and never recomputed;
  - **h_iso is not a state field**: it is computed from `s`, `age`, `plate_type` and `cratonic_mask` by
    `production_upscale::c1_coarse_raw_altitude` (Airy, plus Stein-Stein on oceanic cells, plus the craton
    treatment), then normalised by `c1_normalize_coarse`.
- **The transport**: Eulerian first-order upwind (`step_upwind_masked`) on the fixed 64² grid, with one velocity per
  plate (`fill_velocity_field`).
  - **With `rigid_continental_crust = true`** (the benches' and production's setting), continental cells get zero
    velocity (`apply_continental_rigidity`), and a no-flux boundary stops oceanic crust from flowing into them.
  - Only oceanic crust moves.
- **Following the matter is clean for the land.** A continental column stays in its cell, so the Eulerian h_iso at a
  continental cell IS the material history. There are no fictitious uplift waves: the advected material, the ocean,
  is never land.
  - **The exception**: cells that change type during the run (oceanic → continental by subduction, or accreted). Their
    record before the change belongs to moving oceanic crust. They are counted and reported, not corrected (an ocean
    column becoming land is what happened at that cell).
- **The parameter that chooses the continent type** (the author: « on a déjà ça »): `Phase2InitParams`:
  - `num_plates` (default 8);
  - `cluster.seed_cluster_count` (default 1; the M1 island preset uses 16 plates and 3 clusters);
  - `cluster.continental_fraction` (0.29);
  - `craton_shield_fraction` (Some(0.15): 15 % of the cratonic area as high shield, the rest low platform);
  - `craton_thickness_ratio` (1.25).
  - The benches use the defaults.
- **The gorges' age selector**: `valley_construction::F121_AGE_K` = 0.07183, « a calibration on a TARGET, not a
  measurement of an age » (its docstring). **No link to C1's duration.**
- **The current physics level**:
  - `cascade.rs::uplift_field`: U = U₀ · max(0, h_iso) / 1 000 m, with U₀ calibrated on the peak (one trial, then a
    linear rescale);
  - `run_level` stops at equilibrium: mean |Δz| / mean U·dt < 0.01 for 5 consecutive steps, or at the step cap
    (300 at 256²).

## Seeds (A2, D)

- **The témoin's seed is the author's seed**: `PSEED` = 10481999410520546993 is the very number of the author's
  uncommitted default (`config.rs`). They are one world for C1.
  - The témoin's knobs (valley construction, slope floor) act after C1 and do not touch it.
- **The four seeds**: 10481999410520546993 (the témoin and the author's), **42** (the old default), **1** and **9**.
- The C1 run is the benches' (`build_world`'s): 64², 300 steps, rigid continental crust, `Phase2InitParams::default()`,
  `C1Closures::default()`, `IsostasyConfig::c1_default()`.
- The coarse field is `c1_coarse_normalized_altitude` with the benches' `target_land_fraction`.
- **If none of the four has collision or active-margin cells, that is said** (no seed is swapped after measuring).

## A2 — the diagnostic measures (each seed)

- **The final plates** (distinct `plate_id`), and the cells of:
  - **collision**: Convergent cells whose two sides are continental;
  - **subduction upper plate**: `upper_plate_mask` after `retarget_upper_plate_continental`;
  - **craton**: `cratonic_mask`;
  - **rift**: Divergent cells on continental crust.
  - All from `classify_boundaries` on the final state.
- **The sutures**: every cell classified Convergent at ANY step of the run (recorded read-only through the observed
  callback; F151's fossil-belt notion).
- **h_iso at 64²** (m, the coarse field, (norm − 0.5) · n2m):
  - the land deciles;
  - for h_iso > 500, 1 000 and 2 000 m, the 8-connected components, with:
    - the count;
    - the width, 2 × the median Euclidean distance transform over the component cells (cells and km, 6.25 km per
      cell);
    - the elongation, √(λ₁ / λ₂) of the cells' coordinate covariance, area-weighted over the components;
    - the share of the area within ≤ 2 cells (Chebyshev) of a current convergent cell or a suture.
- **The same measures on the current physics level** (256², F164's regime: p = 6 with roughness, to equilibrium), in
  its own cells (1.5625 km).
- **Verdict A**:
  - **« inherited »** if h_iso > 1 000 m already has a median width ≤ 3 cells (19 km) and ≥ 70 % near a convergent
    boundary or suture;
  - **« regime »** if h_iso is wide (> 3 cells) but the physics level's land above 1 000 m is narrower by ≥ ×2;
  - **« both »** otherwise, with the width ratio given as the parts.

## B — the references (declared before any measurement of our worlds)

### B1 — data

- **The macro reference: Europe**, from **ETOPO 2022 (NOAA NCEI), surface elevation, GeoTIFF**.
  - **Amended before any measurement**: the 30″ product exists only as one 1.5 GB global file, so the **15″ tiles**
    (15° × 15°, 12 tiles, 330 MB) are used instead.
  - Each 1.5625 km cell is the mean of 4 × 4 bilinear samples, i.e. an area mean over ~3 × 3 to 5 × 5 source pixels.
  - The 20 360 cells (0.55 %) beyond the tiles, in the north-west Atlantic, are set to −3 000 m (sea).
  - It is in the public domain (US government work).
  - It is downloaded into its own directory outside the repo (untrusted data; Python run with `-I`), and the derived
    grids go to `data/europe/` (git-ignored).
  - **If the NOAA server refuses, I stop and ask** (no other source is tried without the author).
- **The window**: a square of **3 000 km** on a **Lambert azimuthal equal-area projection centred at 50.0° N, 10.0° E**
  (ETRS89-LAEA's centre, spherical R = 6 371 km), resampled bilinearly to **1 920 × 1 920 cells of exactly 1.5625 km**,
  our 256² cell. It spans roughly 36°–63° N and 12° W–33° E.
- **The alpine control window**: 400 km square, LAEA centred at 46.3° N, 9.5° E, at 1.5625 km.
- **Corsica**: `data/corsica/corse_128.bin` (200 km at 1.5625 km) for the classes, and `corse_1024` (195 m) and
  `corse_4096` (49 m) for the texture and the coast.
- **The qualitative regions** (inside the Europe window, lat/lon boxes):
  - the Paris Basin: 47.8–49.3° N, 1.0–3.5° E;
  - the Massif Central: 44.7–46.0° N, 2.3–4.0° E;
  - the Meseta: 39.0–41.5° N, 5.5–2.5° W;
  - Corsica: 41.4–43.0° N, 8.5–9.6° E.

### B2 — the macro classifier (real metric on both sides; 1.5625 km cells)

Inputs: z (m), the slope (central differences, degrees), and **R7 = max − min of z within a disc of radius 7 km** (the
cells within 4.48 cells). The land is z > 0 m, 8-connected to nothing in particular (inland basins below 0 stay
land).

The classes, in order:
1. **Mountain**: **Kapos et al. 2000** (UNEP-WCMC; the GMBA definition, Körner et al. 2011), built on 1 km DEM cells,
   so it transposes to 1.5 km:
   - z ≥ 2 500 m;
   - or 1 500 ≤ z < 2 500 and slope ≥ 2°;
   - or 1 000 ≤ z < 1 500 and (slope ≥ 5° or R7 ≥ 300 m);
   - or 300 ≤ z < 1 000 and R7 ≥ 300 m.
2. **Hill**: not mountain, **R7 ≥ 150 m**.
3. **Plateau**: not mountain or hill, **z ≥ 500 m** (Meybeck et al. 2001's plateau logic: high and smooth).
4. **Plain**: the rest (z < 500 m, R7 < 150 m).

**The thresholds**:
- Kapos's are the literature's.
- The two others (150 m, 500 m) are declared here from Meybeck's class logic. Meybeck's own roughness thresholds are
  defined on 30′ cells and are **not** transposed (the brief's trap).
- **One adjustment pass on Europe is allowed** before any world of ours is classified: if the qualitative control
  fails on a region for a threshold reason, one of the two may move once, with the reason written. **After that they
  are frozen.**

**The qualitative control** (declared): Corsica mostly mountain; the Paris Basin plain; the Massif Central and the
Meseta plateau. What holds and what does not is listed (the class with the largest share in each box).

**The one adjustment, made on Europe before any world of ours was classified** (`f165_europe_output.txt` as first
run, `f165_europe_output_adjusted.txt` after):
- **As declared**, Europe read plateau 1.6 %, and the Meseta read hill (34 %) over plateau (31 %): the 150 m hill
  threshold caught its moderate dissection at 500–900 m.
- **The adjustment** (a rule order, not a new number): **hills are low**. In Meybeck et al. 2001's typology the hills
  sit below the plateaus, so a non-mountain cell at z ≥ 500 m is a plateau, and the hill test (R7 ≥ 150 m) applies
  only below 500 m.
- **After it**: Europe plateau 4.7 %; the Meseta plateau (61 %).
- The Massif Central still reads mountain (74 %): Kapos's definition counts it as mountain (R7 ≥ 300 m at
  300–1 000 m). **It does not hold, and it is not adjusted.**
- **The thresholds are frozen from here.**
- **Meybeck 2001's global figures could not be read** (BioOne serves a page, not the paper). Nothing is quoted.

### B3 — the macro quantities

- the four classes' fractions of the land;
- the land deciles;
- for the mountain class: the 8-connected components (≥ 4 cells), their count, and **the width / the continent's
  size**:
  - the width is 2 × the median distance transform over the mountain cells, in km;
  - the size is √(the land area);
  - this is the author's « macro-filiform » test.
- **Meybeck et al. 2001's global figures**, if the paper can be read (open access or not); otherwise that is said and
  nothing is quoted from memory.

### B4 — the coast's grain (Corsica the reference)

At 195 m: our worlds at 2 048², Corsica `corse_1024`. Corsica `corse_4096` (49 m) is reported for information. The
coastline is the land/sea boundary of the land mask (z > 0 m, the ocean as `ocean_mask`). Four measures:
1. **The box-counting dimension** of the coastline cells over box sizes 2, 4, 8, 16, 32, 64 cells (0.39–12.5 km): the
   least-squares slope of log N against log(1 / size).
2. **The length ratio**: L(0.39 km) / L(6.25 km), with L(ε) = N(ε) · ε.
3. **The aligned share**:
   - the coastline is traced cell by cell (Moore neighbour boundary of each land component of ≥ 64 cells);
   - it is cut into 8-step pieces;
   - the chord share within ±5° of 0 / 45 / 90 / 135° is taken, as A_dir.
4. **The 64² block index**: the share of the land/sea cell edges lying on the 64² block lines (every 32 cells at
   2 048²), over the 1 / 32 a random coast would give. Corsica is read with a 32-cell lattice too (≈ 1 expected).

- **Degraded** means |m_history − m_Corsica| > |m_F164 − m_Corsica| + 0.10 · |m_F164|.
- **The coast's displacement**: the land gained and lost between F164's final field and the history's, at 256² and at
  2 048², as % of F164's land.

## C — the mechanism

### C1 — the record

- **N = 31 snapshots** of the raw coarse altitude (`c1_coarse_raw_altitude`): the state before the first step, then
  after every 10th step (10, 20 … 300). Each is normalised exactly like the final field.
- **Recorded through `run_with_closures_observed`'s read-only callback.** **A permanent test**: s, age, plate_id,
  plate_type and cratonic_mask are bit-identical with and without the recorder.
- **The material frame**: Eulerian = material on continental cells (A1). The cells that changed type are reported.

### C2 — the historical physics level (256²)

- **The initial state**: h_iso(t₀) upsampled as today (`upsample_to`, bicubic), plus F163's roughness (the « π » of the
  brief, read as the frozen physics-level perturbation: 0.5 × Corsica's 3.1–6.25 km octave at 1 563 m, on the land).
- **The time**: **N_phys = 300 steps over T years** (dt = T / 300). The physics step's snapshot time maps linearly to
  C1's 300 steps.
- **At each step**:
  - **U(x, t) = Δh_iso / Δt** between the two bracketing snapshots (piecewise constant per interval, 30 intervals),
    upsampled as h_iso today, in m / yr, positive or negative, applied to every cell, land and sea;
  - **the erosion**: the implicit stream power (`incise_with_floor`, K, m, n, dt per year as the current level, MFD
    p = 6 as F164); cells ≤ sea level are its base level;
  - the talus and the linear diffusion unchanged; their change is kept on the land only (a cell at or below sea level
    keeps its pre-step value);
  - **the sea is not held at sea level**: it follows U, and the land / sea mask is the current z;
  - **the isostatic rebound**: Δz = (ρc / ρm) · (G ∗ E), applied to every cell.
    - E is the step's stream-power lowering (≥ 0, land).
    - **G is the thin-elastic-plate filter**: Ĝ(k) = 1 / (1 + (k α)⁴ / 4), FFT on the periodic 256² grid, isotropic,
      in real km. **Never in 64² blocks.**
    - ρc = 2 700 kg/m³, ρm = 3 300 kg/m³.
    - **α = (4 D / (ρm g))^¼**, with D = E_y Te³ / (12 (1 − ν²)), E_y = 100 GPa, ν = 0.25, g = 9.81 m/s², and
      **Te = 25 km** (the brief's 20–30 km, the mid continental-interior value of Watts 2001). So **α = 64.4 km**.
- **No equilibrium criterion**: it stops at the end of C1's history. **No deposition** (a known limit). The subsided
  basins below sea level (interior and connected) are measured and reported.
- **The one calibrated constant, T**:
  - calibrated **once, on the témoin** (rebound on, α as declared), so its land peak at 256² falls in **[2 700,
    3 000] m** (target 2 850 m);
  - the method: trials at log₁₀ T = 5, 6, 7, 8, 9, then bisection in log T on the bracketing pair, ≤ 12 trials;
  - **T is then fixed** for the other seeds and every control. K, ρ and α are not touched.
- **The chain**: the frozen cascade (no recalibration) from each physics level to 2 048².

### C3 — the controls

- **Negative 1** (a permanent test): the historical level given a « history » that is U ∝ the final h_iso, constant,
  run to equilibrium, with no rebound, **is the current `physics_level`, bit for bit**. It is the same loop, with F164's
  output hashed before the refactor.
- **Control 2**: the rebound set to 0.
- **Control 3**: the whole Δh_iso (final − initial) concentrated in the last 10 % of T (U = 0 before 0.9 T, the same
  total). The expectation is a little-dissected relief.
- **Sensitivity**: α × 0.5 and α × 2. Measured only.

## D — the measures

- **Macro (B3)**: F164's physics level, the history, and controls 2 and 3, for each seed, beside Europe (and Meybeck
  if read).
- **Texture**: F163–F164's Corsica battery at 2 048² against Corsica at 195 m.
  - It is applied **to the mountain class** (the 256² class map, ×8 nearest, both sides, Corsica classified at 1.5625 km
    then ×8). The non-mountain cells are excluded like D5 cells.
  - It is applied to the whole world as well.
  - F is reported per class.
  - **P5's « elements that passed at F164 »**: the instruments within ×1.5 for p6 · W-tous W1 at F164, i.e. all but
    A_dir and F.
- **The coast (B4)**: F164 against the history, against Corsica; the displacement.
- **The new depressions**: a component after a warp step is new if none of its cells was in a depression before it
  (matched by footprint). The net count is reported beside it.
- **The general checks**: the guard 6 / 6, the cascade's drift, and the s / world (reference 249.8 s).

## R — the stop rules (the brief's, verbatim in substance)

- **Success**, on ≥ 3 of 4 seeds with T only:
  - the mountain, plateau and plain fractions each within ×1.5 of Europe;
  - the mountain width / continent size within ×2 of Europe;
  - P5 holds;
  - and no B4 measure degraded by more than 10 % against F164.
- **« C1 » stop**: P1 holds and the plateau fraction stays under 50 % of Europe's. Then nothing more is tuned, and the
  report describes what C1 lacks.
- **« Régime » stop**: P1 refuted (h_iso wide), but the history does not create plateaus (the fraction under 50 % of
  Europe's on ≥ 3 seeds). Then I report, with no second mechanism.
- **Forbidden**:
  - any per-zone or per-class knob;
  - noise, terraces, clipping;
  - any downstream recalibration;
  - any constant calibrated on our worlds other than T.
