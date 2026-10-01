# F134 — instruments, declared BEFORE the measurements (2026-10-01)

**Two worlds, said so.**
- **The viz world**, `ValleyConstruction::new(F121_AGE_K, Some(0.1))`, slope floor 0.024, canonical framing. It
  is the guarded "C2 /10 col (défaut)" (field `a8d2d538d692c2f0`), and it is where Part O's 44.3 M cells and
  Part D's p90 1 177 m were read (F133v). **Used for O and D.**
- **The benches' témoin**, the same plus `wall_sea_floor_m = Some(0.5)` (F126-B's mur ↔ mer clause, carried
  by every témoin since Finding 127). The viz never sets it. It is where F133's T, U and Δ numbers were read
  (θ 0.493 → 0.436, under-law 32.1 → 49.6 %, residual 12 509). **Used for T, U and Δ**, so that each
  number is reproduced before it is split.
- "ON" = the extended lake base (`InputLakesAndBasins`, the viz toggle) unless a line says F132
  (`InputLakes`).

## 0.2 / 0.3 — the lake guard
- **Key**: `hd_drainage_key(ekey, dcfg_viz, ss, 45, PrecipParams::default(), Some(40), 7.5, 400)`.
- **Fingerprint**, taken AFTER the C-2 crater pass: FNV-1a over `lake_map` (u32 LE, row-major) `/` FNV-1a
  over the lake list (id, `{:?}` of the type, `level_m` bits, `area_km2` bits), in list order.
- **The bench's tail** (`common::viz_hd_lakes_on`) calls the functions `run_hd` calls:
  - pre-breach `c1_drainage_windowed(None)`;
  - `breach_monotone_protected` with `crater_protect_mask` (moved to core);
  - `c1_climate_placed(45, 40)`;
  - `build_hd_infiltration` on `run_coarse_tectonics`' state;
  - `assemble_hd_drainage(.., infil, ..)`;
  - `crater_lake_pass` (moved to core).
- **Identity, lake by lake**: `bench_guard::lake_listing` (id, type, level, area with bits, cell count,
  FNV-1a of the cell indices). Written by the bench and by `run_hd` (`f134v_lake_listing`), OFF and ON,
  then diffed.

## O — the ocean (viz world)
- **Stages**, each built by the pipeline with a `Knobs` switch, OFF vs ON:
  - (a) FBM + craters + construction (+ rims): `no_incision, erosion_off, bathymetry_off`;
  - (b) + the light pass: `erosion_off, bathymetry_off`;
  - (c) + the droplet erosion: `bathymetry_off`;
  - (d) + the bathymetry, i.e. the eroded world.
- **Sea cell at a stage**: `z ≤ 0.5` (norm) in both OFF and ON, split by `water_class` into OCEAN (edge
  connected) and INLAND. The "first difference" is the first stage where an OCEAN cell differs bitwise.
- **Mean depth**: exactly `apply_bathymetry_profile`'s, on the stage (c) field: the mean over cells with
  `z ≤ 0.5` of `(0.5 − z)·depth_per_norm`, in f64. `depth_per_norm = 2·1.13·depth_scale_m`.
- **Ratio**: `z_ON / z_OFF` in metres on the stage (d) cells that are sea in both.
  - Its median, and the share within 10⁻³ of the median.
  - With `texture = 1` the profile is `depth = env(dist)·c/m`, so where the pre-bathymetry depth `c` is
    bitwise equal, the predicted ratio is `m_OFF / m_ON` exactly, clamps aside. The share within 10⁻⁵ of that
    prediction is reported, split by: `c` equal / `c` different / clamped (depth at `shelf_min` or at the
    floor in either world).

## D — lakes 2, 11, 1 (viz world, F133v's ids, found by centre)
- **Field**: F133v's (`HdResult.eroded` = the CONDITIONED field, no water). The pre-breach eroded field is
  read beside it.
- **Δz** = z_ON − z_OFF (signed).
- **Populations** (F133v's attribution: the OFF drainage's D8 path to the first ON final lake):
  - (i) the lake's ON footprint;
  - (ii) outside every ON lake and attributed to it ("upstream").
  - Each is split into **carved** (lowered by the construction in OFF or ON: `carve` on the construction
    input S1, with `mk.carved`) and **kept at the terrain** (carved in neither).
- **Per population**: count, p50 / p90 of |Δz|, the shares Δz > 0 and < 0, and Σ|Δz|·A.
- **Beside it**: the lake's ON level and depth, and the col (`ocean_flood` spill on S1's breached field, at
  the lake's cells, max).
- **Map**: PNG, centre ± 307 cells, Δz coloured, the ON and OFF lake outlines.

## T, U (benches' témoin; OFF skeleton's trunk cells ≥ 10 km², as in F133)
- **Based set** (per variant): F132 = input lakes; extended = input lakes ∪ `basin_base` closed
  depressions.
- **Populations of a link** (by its upstream cell):
  - upstream = the OFF skeleton's D8 path reaches the based set;
  - downstream = some ancestor on the OFF skeleton is in the based set (propagated in increasing-area
    order);
  - untouched = neither (cells inside the based set are excluded).
- **θ** = `theta_laid_ci` with the carved mask restricted to the population (both ends in it).
- **"Off their law > 1 m"**: carved trunk cells with |z_built − floor_m| > 1 m (the variant's own
  skeleton law).
  - Distance to the based set by D8 (Chebyshev) BFS in cells, capped at 1 000 cells (48.8 km).
  - Shares within 2 km; p50 / p90 of the distance.
- **U**: the cells under their law in ON but not in OFF.
  - Distance to the based shore (as above): p50, p90, share < 2 km.
  - `law − terrain` (floor_m − z_S1): p50, p90.
  - The under-law run each sits in, walked downstream on the OFF skeleton while under the law: its length
    (p50, p90) and whether it ends at the based set's shore (the "delta" signature) or elsewhere.

## Δ (benches' témoin, extended; F133-F's A/B flat-pointer patches)
- **Residual** = cells outside the present lakes with z_A ≠ z_B, with p50 / p90 / max of |Δz|, overall and
  per F133-F class.
- **The laying sample**: `carve_diag`'s `who` in A and in B, compared by position (x, y), floor `zf` and
  half-width. The classes are: moved, same position but another floor, no sample on one side, and the
  confluence band.
- **The function** is named from the code.
