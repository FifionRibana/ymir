# Finding 150 — predictions, written 2026-10-06 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–149 and the reviewer's predictions.

**Read in the code before writing these (part I's code reading, no measurement):**
- **`C1State` (`tectonics_c1/state.rs`)**, at the 64² tectonic grid:
  - `s` (the crust thickness, non-dimensional), `age`, `plate_id`, `plate_type` (continental / oceanic);
  - `cratonic_mask` (a BFS at init).
  - The production upscale reads `s`, `age`, `plate_type` and `cratonic_mask` for the isostasy
    (`production_upscale.rs:307–310`).
- **`debug_labels::CoarseTectonicLabels`** (continental, craton, rift, subduction upper / slab, collision, divergent,
  age_norm):
  - built for the viz's "Tectonique" overlay and carried in `HdResult.tectonic`;
  - **not** in production; nothing exported.
- **The viz cell panel hard-codes the three fields** (`workspace.rs:1922–1924`): `kv(…, "—")`, with the comment
  "coarse — deferred".
- **C-3 lithology** is an ERODIBILITY multiplier with three classes: hard basement ×1, rift-soft ×10 (continental
  `age < 1`), volcaniclastic ×3 (edifice discs).
  - **C-3b fracture** is the density from the convergent and transform contacts, giving K = 1 + 6·density.
  - **C-2 volcanism** places the edifices and craters.
  - **H-1 infiltration** gives a permeability.
- **In the viz the four closures are opt-in, OFF by default. In the benches and in the 6/6 guard's states, C-2, C-3
  and C-3b are ON.**
- **The export writes no geological layer**: height, coastline, cliffs, temperature, precipitation, biome,
  lake_mask, flow_accumulation, rivers, rivers_ll, lakes, water_class.
- C-3's own comment: production erosion is detachment-limited, so there is NO sedimentary-basin signal and no
  deposit field.

## Predictions

- **P-I**:
  - "Craton" and "Type de plaque" ARE computed for a C1 continent: `cratonic_mask` and `plate_type` at 64², both
    filled, and in `HdResult.tectonic` whenever the viz requests labels. **The panel is empty because it is
    hard-coded**, not because nothing is computed. "Épaisseur crustale" (`s`) is computed but never carried to the
    viz.
    - Against the reviewer's "not computed for a C1 continent".
  - **No "type de sol" with classes exists** (with the reviewer). The nearest are:
    - the C-3 lithology classes (hard / rift-soft / volcaniclastic: three, erodibility only);
    - the Whittaker biomes;
    - `water_class`.
  - **The age is degenerate** (C-3's comment): its map will show little structure.
- **P-D** (C-3 ON against OFF, C-2 and C-3b ON in both, the témoin):
  - the land relief p50 changes by **< 2 %** (with the reviewer's < 5 %);
  - the land |Δz| p50 is < 1 m, p90 < 50 m, max several hundred metres;
  - **≥ 80 % of the |Δz| > 10 m cells lie on or downstream within ~5 km of the rift-soft class**, a few on the
    volcanic discs (the reviewer: "rift / volcanoes", agreeing);
  - R8 changes by < 0.01.
- **P-C**: with the existing fields, **8 of the 11 starting resources can be zoned as a possibility**:
  - gold (craton + collision + placers along rivers), iron (craton; bog iron in wetlands), copper (subduction upper
    plate), tin (collision), salt (endorheic lakes + arid climate), building stone (by class), clay (floodplains,
    lake flats), sulphur and obsidian (edifices);
  - **not with the existing fields**: coal and the carbonate-hosted Ag-Pb-Zn (they need sedimentary basins, and Ymir
    has no deposit field), and gemstones (pegmatite and metamorphic-grade proxies are missing; only their placers
    follow the rivers).
  - With the reviewer (≥ half), larger.

## Meta

My disagreement with the reviewer: P-I (computed, not wired). **At least one of my predictions is wrong.**
