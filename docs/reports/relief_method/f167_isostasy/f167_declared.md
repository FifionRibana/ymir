# F167 — declared before any measurement (2026-10-10)

**The thickness → altitude mapping: physical isostasy and an absolute sea level, in the cascade only.**

**Production, C1 and the guard 6 / 6 are unchanged.** The cascade is frozen (p6 · W-tous W1, F164's budgets and k).
**T stays at 5.62 Myr** (F166), with no recalibration. **No constant is calibrated on our worlds in this finding**:
every constant of the mapping comes from the literature, cited below by file and page.

**Not blind**:
- F166's readings;
- the code read for A1–A2;
- `f167_inventory_output.txt`: s on the oceanic and continental cells of the four seeds at t₀ and at the end, read
  before the declaration to decide B2 (below). No world was mapped or measured with the new mapping.

## The sources, read (`docs/refs/`)

- **Whitehead & Clift 2009** (`2009_Whitehead_ContinentElevation_JGR.pdf`, JGR 114, B05410; the reviewer's
  « Whitehead 2009 »):
  - **p. 3, §2.1**:
    - « average stable, continental crust has a thickness of ~38 km [Christensen and Mooney, 1995] and that crust now
      at sea level averages ~35 km thickness »;
    - « Using an average continental crustal density of 2800 kg/m³ and mantle density of 3300 kg/m³ (… textbook values
      … from Turcotte and Schubert [2002]) we estimate that 42.5-km-thick crust stands at ~1.5 km elevation above
      modern sea level, while 58-km-thick crust rises to 3.8 km »;
  - **p. 7, eq. (1)**: ρo = 1 030 kg m⁻³ (seawater), ρc = 2 800 (continent and sediment), ρm = 3 300;
  - **p. 8, eq. (2)–(3)**: e = (3300 − 2800)/3300 · (dc − ds) + 1030/3300 · do + C, with C = −2.2580 km, the mean
    continental elevation 0.835 km and the mean crust dc = 38.037 km.
- **Lachenbruch & Morgan 1990** (`lachenbruch1990.pdf`, Tectonophysics 174, 39–62):
  - **p. 42, eq. (4)–(7)**: the reference column is the mid-ocean ridge, with ρa = 3.2 g/cm³ (asthenosphere),
    ρl = 2.8, L = 5.5 km, ρw = 1.0, and **E = −2.5 km**, giving A₀ ≈ 3.5 km and H₀ ≈ 2.4 km;
  - **p. 42–43, Fig. 2**: the curve ρl = 2.8 passes through E = −2.5 km at L = 5.5 km (the ridge).
- **Talling et al. 1997** (`1997_Tallingetal.pdf`, Basin Research 9, 275–302):
  - **p. 276**: « For 11 orogens, the spacing ratio only varied between 1.91 and 2.23 » (Hovius 1996);
  - **p. 277, eq. (1)**: R = W/S, W the topographic half-width (mountain front to the ridge crest), S the outlet
    spacing.

**Where the reviewer's citations differ from the papers** (the papers prevail):
1. **Whitehead's three numbers are not one Airy line.**
   - With the slope (3300 − 2800)/3300 = 0.1515, 42.5 km → 1.5 km puts sea level at 32.6 km, and 58 km gives 3.85 km
     (the paper's 3.8).
   - So the paper's own isostatic numbers put **35 km at +0.36 km and 38 km at +0.82 km**, not 35 km at sea level.
   - The « ~35 km at sea level » is the paper's seismic statement (Christensen & Mooney), not its Airy line.
   - The brief's « 38 km → ~+0.45 km » follows neither the seismic statement (+0.45 needs sea at 35 km) nor the Airy
     line (+0.82). **The paper cannot settle the anchor** between these readings.
2. **Lachenbruch & Morgan's ridge column is a LITHOSPHERE column** (mean density 2.8 over an asthenosphere of 3.2,
   ρw = 1.0), not crust over a 3.3 mantle.
   - Their −2.5 km ridge depends on the hot, light ridge mantle, which a crust-only Airy (ρm 3 300) cannot represent:
     with the brief's densities a 5.5 km crust lands at −5.8 km.
   - So the ridge cannot anchor the same function as the continents' freeboard.
3. « Whitehead 2009 » is **Whitehead and Clift 2009**.
4. Talling / Hovius: **confirmed as cited** (1.91–2.23 over 11 orogens).

## A1 — C1's current mapping

- **The land ramp** (`tectonics/isostasy.rs::compute_isostasy_inner`, lines 268–360):
  - h_raw = s · b, with b = 1 − ρc/ρm: ρc 2 750 (`IsostasyConfig` default, line 113), ρm 3 300; the cratons
    ρc = 2 900 (`c1_default`, `craton_rho_crust`);
  - **the datum** h_sea = h_min + 0.4 · (p92(h_raw) − h_min), a percentile of the WHOLE field (lines 281–292);
  - the land: norm = sea_norm + (h − h_sea)/(land_ceiling − h_sea) · (1 − sea_norm), with land_ceiling = 2.0 · b (S̃ =
    2.0) ↔ 5 650 m (`max_elevation_m`, `land_ref_thickness`, lines 314–322, 177–183);
  - then a Gaussian blur σ = 0.5 (line 367);
  - then to metres (`production_upscale.rs:111–130`).
  - **Oceanic cells are overwritten** by Stein-Stein 1992: ridge 2 600 m → asymptote 5 651 m, by age
    (`oceanic_bathymetry`, `age_to_ma` 0.667).
  - The result is normalised by `c1_normalize_coarse` ((raw + 1.13)/2.26, clamped to [0, 1]; **5 650 m is the ceiling**
    of the normalised field).
- **The convention of s**: S̃ = 1 ↔ ~35 km, normal continental crust; S̃ = 0.2 ↔ ~7 km oceanic; S̃ = 2 ↔ ~70 km
  (`docs/c1_lightweight_dynamic_tectonics.md` §11, table).
- **The unit of `age`**: a « step », converted by `age_to_ma` = 0.667 Ma for Stein-Stein.
  - **But `age` is not a clock**: it is initialised (ridge-aligned), advected and reset to 0 at rift splits; **it is
    never incremented in time** (`time_loop.rs`, `init_r7/age_init.rs:174`, `rifting/split.rs:277`).
  - So **no thermal subsidence is added** (B2).

## A2 — what reads h in C1

- **The erosion closure** (term 4) reads the altitude (isostasy + Stein-Stein) for its slope, and drains on an S̃-space
  sea level (`compute_sea_level_ref_s_space`).
- **Davis-Suppe reads s, not h**: its h_crit is a thickness target in S̃. So do the equilibrium sink, the subduction,
  the rift and the transport.
- **So correcting the mapping in the cascade only leaves out exactly one thing**: C1's erosion closure keeps eroding
  with C1's slopes. Its Δs is the term the sources' drive already excludes (F166).
- The approximation leaves C1's own erosion field, and therefore C1's s (which the steady F164 level reads), as they
  were. **The sources' drive itself is not approximated**: s_src carries no term that reads h.

## B — the mapping (cascade only)

- **B1, the continental side**:
  - Airy with ρc = 2 800 and ρm = 3 300 kg/m³ (Whitehead & Clift p. 3), so a slope of 0.1515 km of altitude per km of
    crust, ×35 km per unit S̃ (C1's convention) = **5.303 km per unit S̃**;
  - **the anchor**: the paper cannot settle between 35 and 38 km (above), so I take **the brief's central value:
    S̃ = 1 → +0.45 km** (the 0.3–0.6 km freeboard band). That fixes the constant: h = 5.303 · S̃ − 4.853 km above sea
    level;
  - **the high check**: S̃ = 2 (70 km) → +5.75 km, inside 5–6.5 km and beside C1's 5.65 km and Whitehead's line
    (70 km → 5.67 km).
    - The normalised field's ceiling clamps altitudes at 5 650 m (S̃ > 1.98), as C1 does. Reported, not corrected.
  - **One density for all the crust**: C1's craton density 2 900 is not in these sources and is not used. The cratons
    (×1.25 thickness) map to +1.78 km.
- **B2, the oceanic side**:
  - the same function, water-loaded below sea level: e = (5.303 · S̃ − 4.853 km) · ρm/(ρm − ρw), ρw = 1 030 (Whitehead
    & Clift p. 7);
  - continuous at 0 and monotone in S̃: S̃ = 0.2 (7 km) → −5.51 km;
  - **the ridge column**: Lachenbruch & Morgan's −2.5 km cannot be the reference of this crust-only function (above).
    It is reported as a check: our function puts a 5.5 km crust at −5.8 km;
  - **no thermal subsidence** (`age` is not a clock);
  - **the oceanic cells' s, declared from the inventory**: at the end of C1, 76–186 oceanic cells per seed have
    S̃ > 0.915 (up to 2.18). These are C1's documented « phantom oceanic advective spike » (`production_upscale.rs:90–99`),
    which C1 itself never maps to altitude (Stein-Stein overwrites it). Under one function of s they would stand as
    islands up to 5.7 km.
    - **Main**: the oceanic-plate cells are read at C1's oceanic reference thickness, **S̃ = 0.2** (`init_s_field`).
      The function is the same; only the input of the cells whose s C1 does not trust is substituted.
    - **Variant « own s »**: every cell is read with its own s. Its phantom islands are counted.
- **B3, the sea level: 0 m absolute.** No percentile, so the ocean–land coupling of the datum is gone.
- **The grid filter**: C1's Gaussian blur σ = 0.5 is applied to the mapped altitude, so that only the mapping changes.
- **B4, the drive** (F166's, unchanged):
  - h_src = the new mapping of s_src (C1's erosion excluded);
  - U = Δh_src/Δt;
  - the initial state is the new mapping of s(t₀), plus the roughness;
  - **T = 5.62 Myr** fixed; the peaks reported.
- **B5, the variants and the controls**:
  - **the anchor ±0.2 km** (S̃ = 1 → +0.25 and +0.65 km), measured only, at 256²;
  - **the decoupling test**: « sources » against « sources − EH » (the EH sink excluded too), under the new mapping
    (main and « own s »). **The land must be identical within 0.1 %**;
  - **the negative control** (permanent): with C1's mapping, the record is F166's bit for bit;
  - « own s » as above.
- **B6, the viz**: « Correspondance : C1 / Airy » in the Cascade window.

## C — the measures (seeds: the témoin, 42, 1, 9; physics 256²; the frozen cascade to 2 048²)

- **Macro**: the classes against Europe; the cratons' altitude (p50); the land lost and gained against F164; the
  interior at ≤ 0 m; the lakes.
- **The peaks and the slope guard** (mountain-class slope p50 / p90 against Corsica).
- **The texture**: the Corsica battery, whole and mountain class; F per class.
- **The coast**: F165's B4 against Corsica, and « degraded » against F166.
- **The regularity baseline** (for F169), on F164, F166 and F167 at 1 024² (391 m) and 2 048² (195 m), against
  Corsica:
  - CV_λ, CV_L, σ_θ and the mean confluence angle, R2_g (`cascade::predict`);
  - **Hovius's ratio R = W/S on the mountain fronts** (the mountain class: the 256² classes ×8):
    - an outlet is a mountain-class cell whose D8 receiver lies outside the class;
    - every mountain cell belongs to the outlet it drains to;
    - **amended at its control, before any measurement of our worlds: R is read per transverse basin**, R = the mean
      of L²/A, where A is the basin's area inside the class and L the largest distance from the outlet to a cell of the
      basin (W for a transverse basin; S = A/L its width along the front);
      - **why**: the declared « front length / outlets » estimator read **1.24** on the synthetic range built at
        W/S = 2. The front's ends and its cell staircase inflate S.
      - The per-basin form reads it within ±35 % (`hovius_reads_a_synthetic_range`);
    - A_min = 25 km² (main), 10 and 100 km² reported;
    - the same on Corsica (its mountain class, its fronts mostly coastal);
    - the reference 1.91–2.23.
  - For each instrument: where we are more regular than Corsica, and at which cell.
- **The images**: the four seeds F166 | F167 (shaded and classes, 256² and 2 048²); the mapping h(s), C1 against
  Airy; 25 km coast crops Corsica | F166 | F167.
- **The cost** in CPU seconds per world.

## R — the stop rules (the brief's)

- **Macro success**, on ≥ 3 of 4 seeds:
  - the mountains in 16.5–37.2 %;
  - plain + plateau in 37.5–84.5 %;
  - the land lost under 15 %;
  - no coast degraded by more than 10 % against F166.
  - On success, F168 follows.
- **If the mountains still exceed 37.2 %**: decompose them (cratons, Davis-Suppe zones, arcs) and stop. No new
  mechanism.
- **If the land lost exceeds 15 %**: report the anchor and its ±0.2 km sensitivity, and stop. The anchor is not
  calibrated on our worlds.
- **Forbidden**: as F165–F166. No calibrated constant.
