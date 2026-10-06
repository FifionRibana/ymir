# F152 — the build's parameters and the gates' instruments, declared BEFORE the measurements (2026-10-06)

**Before / after**: the guard 6/6 (`ymir-viz f123_viz_guard`), lib and viz tests, `cargo check --workspace`.
"Before" = F151's after-checks (`5a124bc`).

## The build (spec F151, with these concrete choices)

- **The history**: `run_with_closures_observed` in a geology pass of its own (0.2 s at 64²). The terrain's cached
  tectonic run is untouched, so the history cannot change the field.
- **The fossil belt** (F151's PROXY): collision ≥ 10 steps, or a continental suture cell, and not convergent at the
  end. **The active belt**: convergent at the end, dilated by 2 coarse cells.
- **The structural density** (the zoning only): `max(C-3b density, exp(−d_km / 25))`, where d is the distance to the
  sutures and the past collision cells.
  - At 64², sampled bilinearly at the cell.
  - It is never passed to the incision. A test pins that C-3b's own density is unchanged by building it.
- **The contacts** (craton, belt, rift):
  - the coarse mask is sampled bilinearly onto a 1024² grid with the terrain's own mapping;
  - a Gaussian of σ = 8 cells at 1024² (= 0.5 coarse cell ≈ 3 km, the spec's PROXY);
  - then bilinear to the cell, and > 0.5.
- **Volcanic**: within 1.5 basal radii of an edifice, the cells whose altitude stands above the base level (the median
  of the 1.3–1.6 radii ring) by > 10 % of the edifice's height. Its class follows its setting: arc → arc volcanic,
  else basaltic.
- **Loose deposits** (PROXY):
  - floodplain: slope < 0.02 and A ≥ 1 km²;
  - lake margin: a land cell within 2 cells (≈ 100 m) of a lake;
  - slope foot: slope < 0.05 within 3 cells of a slope > 0.3.
- **Evaporites**: an ENDORHEIC lake's footprint (the playa), plus its land within 1 km where precipitation < 400 mm.
  - **A deviation from the spec, declared**: the spec put lake water in class 0. An endorheic lake's bed is its
    evaporite flat, so it is class 7; exorheic lakes stay 0.
- **`valley_floor`**: the construction's floor mask is not in `run_hd`. So `valley_floor` = the floodplain criterion
  above (PROXY, declared).
- **The context at ¼ (2048²)**:
  - rock = the majority of the 4×4 block's land cells;
  - slope = the block mean;
  - A = the block max;
  - lake, wetland and endorheic = any cell of the block;
  - climate = the block mean;
  - the distances on the ¼ grid (a chamfer);
  - edifices: the distance to the centre minus the basal radius.
- **The placers**: on the HD rivers (A ≥ the rule's area), upstream to downstream in increasing accumulation.
  - The source = the HD cells whose ¼ non-placer favourability is ≥ `source_min_chance`.
  - The value carried loses `decay_per_km` per km of D8 path, and is set where the slope is < `max_slope`.
  - Then the max into the ¼ grid.
- **`modulate = "structural"`**: favourability × g(d), with g = g_min + (1 − g_min)·min(d / d_ref, 1).
  - **g_min = 0.2** (DECISION, declared in the predictions with its consequence).
  - **d_ref = the p90 of d over the land cells of the belts (fossil and active) and of the arc zone (within 30 km of
    an upper-plate cell)** on the témoin. MEASURED first, then written into the rules file's `[modulation]` and fixed.
  - Set on the vein rules only: gold (belt, arc), silver (arc, belt), tin (belt), copper (arc).
- **The export**:
  - `geologie_roches.png` (u8 class codes, full resolution, north at the bottom like every raster: y = 0 = south);
  - `favorabilite_<id>.png` (u8 0–100, ¼);
  - `geologie.json` (the rules' SHA-256, the classes, the resources);
  - documented in `docs/geology_format.md`.
- **The zoning cache**: the context (rocks + fields) is kept with the world. A rules reload reruns only the zoning.

## The gates' instruments

- **G-field**:
  - the guard 6/6;
  - the test of the structural density above;
  - the geology stage reads the world and writes only its own buffers.
- **G-blocks**:
  - the land cells on a rock contact (a 4-neighbour class change, both cells on land), and the ¼ favourability edges
    (a step of ≥ 10 points between 4-neighbours);
  - the share of each within 1 HD cell (¼: 1 cell) of a coarse grid line, that is where the coarse coordinate's
    fractional part is within 1/128 of 0 or 0.5, either axis;
  - the control is the same share over all the land cells;
  - **pass if the ratio ≤ 1.5**.
- **G-circles**:
  - each volcanic patch (4-connected, class 3 or 4, ≥ 20 cells): Q = 4πA/P², P = the count of its 4-neighbour edges;
  - against Q of a rasterised disc of the same area, computed the same way;
  - **pass if Q_patch / Q_disc ≤ 0.8 for every patch**; the patches and their ratios are reported.
- **G-nonempty**: per resource, the ¼ land cells with favourability ≥ 10. Pass if > 0, or the absence explained.
- **G-placers**:
  - the HD cells set by a placer are all on a river (A ≥ the rule's area);
  - each has a source upstream (by construction, asserted by a test);
  - the ¼ cells they reach are reported.
- **G-rules**: a permanent test. A synthetic context, two rules files differing in ONE rule's chance; only that
  resource's grid changes.
- **G-veins**:
  - over the land cells of the belts, the modulated gold favourability: p10 / p50 / p90, and the share above half the
    belt maximum;
  - **the control**: the same with d = C-3b alone;
  - the change in the mean favourability over the FOSSIL belts.
- **The cost**: the geology stage's seconds (history, rocks, context, zoning) and the re-zoning alone.

## Amendment after run 1 (2026-10-06, NOT blind: run 1 read)

**Run 1** (`f152_geology_raw_run1.txt`):
- **G-blocks FAILED on the coarse-sourced contacts** (craton / belt / rift / basement): 4.01 × the control. All contacts
  read 1.29 and the favourability edges 0.98.
- **The cause is geometric**: a long straight edge of the coarse mask (a coarse row or column) stays straight, and ON
  the grid line, after any symmetric smoothing. A larger σ cannot fix it.

**The amendment: the spec's optional "calage sur le terrain"** (§1.2), declared before run 2:
- The coarse sources' smooth field gets a terrain term:
  - `value + β·a` for the hard classes (craton, belt);
  - `value − β·a` for the soft rift fill.
- `a = clamp((z − z̄) / 100 m, −1, 1)`, where z̄ is the altitude smoothed at σ = 4 cells of the 1024² grid (≈ 1.6 km).
  β = 0.25.
- **Its rationale (PROXY)**: differential erosion. Hard rock holds the ridges, soft rock the hollows, so a contact
  runs along the valleys.
- β = 0.25 lets a contact move within the Gaussian's transition band, ≈ ±2 km, the spec's "within 2 km".
- **Run 2 re-measures every gate**, with the file's d_ref set to run 1's measurement (0.994, already the value run 1
  used in memory).

