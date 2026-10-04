# Finding 146 — the rivers for Living Landz, round 1 of 2: `rivers_ll.json` (main stems, Chaikin in the valley, attributes, and the hex edges on the author's 40 m flat-top axial grid) and the viz layer are built, at ~0.6 s per world; the polylines leave the stairs (turns ≥ 45° 19 % → 4.5 %) within 0.1 cell of the D8 trace; two things remain open — the 1 m bed instrument fails its own negative control, and the D8 network already crosses itself (574 X-crossings, 662 after smoothing)

**Status: built and measured. The selection is the author's (nothing fixed). Round 2 is the author's look in Living
Landz.**

**Raw outputs, in this folder:** `f146_rivers_raw.txt`, `f146_hex_raw.txt`, `checks_before.txt`,
`checks_after.txt`. The predictions
were written before any measurement (`f146_predictions.md`). The methods and instruments were declared before
(`f146_declared.md`). The format is **`docs/rivers_ll_format.md`**.

## The author's criterion (recorded as given, at the head of the work)

- Rendu : « Polylignes lissées en SDF (comme les routes) ».
- Rôle en jeu : « Obstacle, gué, bonus (bords d'hex) ».
- Sélection : « Pour l'instant pas décidé. Certainement les grandes et celles qui ont un intérêt. Long cours d'eau
  même petit peut avoir un intérêt. À voir une fois la carte rendue dans LL. »
- **Critère d'acceptation** (the reviewer's, to be validated by the author at round 2): in Living Landz the rivers
  shown follow their valley floor, with no staircase blocks and no bundles of parallel traces, and the selection is
  tuned without remaking the world.
- **Limit**: 2 rounds. This one builds and measures; round 2 is the author's look in Living Landz, then the
  adjustments.
- Cost principle: seconds per world, against 249.8 s.

## Partie 0

1. **F145 committed** (`5132caf`).
2. **The cones are closed as an ACCEPTED defect.** `docs/adr/accepted_defects.md` is created (§ 1), and ADR F145
   records it at F146's commit, with the reasons:
   - the defect is measured against the pre-carve divides;
   - 64 % of its cells lie within 2 cells of the old divide, where the ridge is the walls' intersection;
   - the captures are under the noise;
   - C1 adds a seam.
   - The exception is the drained bowls (the paused gorge).
3. **The author's decisions (2026-10-04)** are recorded in ADR F146:
   - the order of work with its round limits: cones, rivers, geology v1, lakes / gorges / falls resumed, resources;
   - the three answers on the rivers, above.
4. **Checks**: before = F145's after-checks (`5132caf`; no file under `src` changed since); after = `checks_after.txt`.

## Q — the code, before anything

1. **`rivers.json` today** (`export/hydro.rs`): per D8 segment,
   - `points` (cell indices), `strahler_order`, `avg_flow` / `max_flow`, `basin_id`, and the links `upstream` /
     `downstream`;
   - plus `catchment_km2`, `drainage_km2`, `navigability`, `discharge_m3s`, `width_m`, `profile_m` (the bed per
     point), `kind` (Watercourse / Spillway) and `source_lake_id`, and `features` (the gorge's falls, omitted when
     empty).
   - **Living Landz's side in Ymir**: only consumer notes in the ADR (F23: LL renders rivers at constant width; F47 /
     F86 on `strahler_order` and spillways). There is no LL code here.
2. **The polylines are D8, cell by cell** (`flow.rs`, `extract_rivers` / `trace_segment`): cell centres, cut at the
   junctions. **The staircase is born there**: each step is a multiple of 45°.
3. **Ymir does NOT know Living Landz's hex grid.** There is only `km_per_cell` "≈ one hex" (`container.rs`) and a
   bench-side assumption of a 1 km-edge hex "for the author to rescale" (`h1_infiltration_sweep.rs`). No size,
   orientation, origin or coordinate system. **Part H stops here** (questions at the end).

## Built

- **`ymir_core::export::rivers_ll`** (new):
  - the main stems (at each confluence the largest-catchment upstream reach continues; every segment in exactly one
    river);
  - Chaikin ×3 with both endpoints and every confluence point fixed (a tributary's last point is a vertex of the
    river it joins);
  - the valley constraint (a vertex ≤ 1 m above its nearest trace point's bed, else pulled back by bisection);
  - Ramer–Douglas–Peucker ε 0.1 cell;
  - the attributes (A).
- **`rivers_ll.json`**, format 0.1.0, a new container layer (`rivers_ll`), documented in `docs/rivers_ll_format.md`.
  **`rivers.json` is byte-identical** (a test, and measured on the témoin).
- **The viz**: `HdResult.rivers_ll` (timed in the log), a "Rivières LL" toggle with a filters menu (min length km,
  min catchment km², min Strahler, logarithmic sliders) showing the rivers kept and their km. The layer draws the
  kept rivers with a stroke by `width_m`.
- **Tests**:
  - `chaikin_keeps_fixed_points_and_smooths_a_staircase` (negative control: the staircase turns at every step);
  - `rdp_keeps_fixed_points`;
  - `rivers_ll_is_a_new_layer_and_leaves_rivers_json_byte_identical`.
  - `hex_of_flat_top_axial` (Part H);
  - lib 597 → 601.

## L — smoothed in the valley (`f146_rivers_raw.txt`; the témoin, 16 226 segments → 9 786 rivers)

| | terrain − bed (m) p50 / p99 / max | length above 1 m | ridge crossings (samples) | stair index (turns ≥ 45°) |
|---|---|---|---|---|
| **raw D8 trace (negative control)** | 0.00 / 12.28 / 559.4 | **17.27 %** | 640 of 1 854 451 | 0.194 |
| smoothed, before the constraint | 0.04 / 14.38 / 561.8 | 23.33 % | — | — |
| **smoothed, after the constraint** | 0.13 / 12.05 / 576.0 | **10.69 %** | **1** of 3 576 377 | **0.045** |

- Lateral deviation from the D8 trace: p50 0.000, **p99 0.094 cell** (≈ 4.6 m). 23.1 % of the vertices were pulled
  back by the constraint.
- **The bed instrument fails its own negative control.** The raw D8 trace, in the valley by construction, already has
  17.3 % of its length more than 1 m above the bed. Bilinear terrain between cell centres mixes in the neighbouring
  cells, i.e. the valley walls wherever the floor is one cell wide.
  - **So the 1 m tolerance cannot judge "in the floor" at sub-cell positions.** The smoothed line does better than the
    raw trace on it (10.7 % against 17.3 %), and crosses a ridge once against 640 times.
- **What locates the line is the lateral deviation**: 99 % of it lies within 0.094 cell of the D8 trace. The smoothing
  takes the stairs off (0.194 → 0.045) without leaving the trace's cells.
- **Crossings between two rivers away from a confluence**: raw **574**, smoothed **662**. Self-intersections: raw 613,
  smoothed 271.
  - **The D8 network already crosses itself**: two neighbouring cells flowing diagonally across each other (an "X")
    is allowed by D8.
  - The smoothing does not remove those, and adds about 15 %. **The acceptance's "no crossing outside a confluence"
    FAILS as built**: open, for round 2. A repair (separating the two strokes at an X) is not made.
  - The raw self-intersections (613) are not attributed (a chained trace should not revisit a cell); the instrument
    may count the torus seam or the chain's joins.

## P — the parallel bundles (k = 2 cells, L = 2 km; no hydrology touched)

- **1 599 pairs**, runs totalling 5 354 km.
- By the smaller river's catchment:
  - **10–100 km²: 1 465 pairs (4 952 km)**;
  - 1–10 km²: 84 (187 km);
  - ≥ 100 km²: 50 (215 km).
  - The rivers start at the drainage's head threshold, p10 catchment 6.7 km², so "under 10 km²" is rare by
    construction.
- Where the runs lie:
  - **88.1 % "elsewhere"**, i.e. not on the construction's carved cells: the uncarved terrain;
  - 5.3 % valley floors;
  - **3.9 % planar walls (F124)**;
  - 2.7 % volcanic edifices.
  - Planarity of the uncarved slopes was not measured.
- **Some pairs are a river and a spillway on the same course.** Rivers 6206 and 9779 share a 4 673.6 km² catchment
  (46.7 against 2 591 m³/s, 19.6 km); likewise 8472 and 9785. The spillways are appended reaches that retrace a
  watercourse.
- **The proposed export rules (not applied)**:
  - (a) merge each pair into its higher-discharge river: removes **921 rivers, 3 431 km**;
  - (b) drop the rivers under 10 km²: removes 3 015 rivers (1 581 km) but leaves 1 515 of the 1 599 pairs. Under 1 or
    3 km² it removes nothing.
  - **(a) acts on the bundles; (b) does not.**

## A — the selection attributes (per river; no threshold fixed)

- **9 786 rivers.**
  - Length: p10 0.18, p50 1.19, p90 5.18, max 59.2 km.
  - Catchment: p10 6.7, p50 15.6, p90 304, max 72 522 km².
  - Strahler (the max over its segments): 1 → 9 269, 2 → 419, 3 → 70, 4 → 13, spillways 15.
- Ends:
  - confluence 7 524;
  - lake 1 767;
  - endorheic lake 409;
  - **sea 82**;
  - terminal 4.
  - 25 rivers touch a lake on their course.
  - The rivers are clipped at the lakes (F37), so a course through a lake is two rivers.
- **"Long, even small"**: no river ≥ the p90 length (5.2 km) has a catchment < 10 km². Long rivers are not small in
  this world.
- For the author's choice, what a few filters keep:
  - length ≥ 5 km: 1 036 rivers;
  - length ≥ 10 km: 275;
  - catchment ≥ 10 km²: 6 771;
  - catchment ≥ 100 km²: 1 653.
- Falls: empty in every river (the gorge is paused). The field exists, in the gorge's format.

## H — the hex edges (`f146_hex_raw.txt`; the grid given by the author mid-round)

**Q.3 stopped H. The author then gave the grid (2026-10-05)**: « 40 m de rayon environ ; flat top ; l'origine est en bas
à gauche ; axial ». Read as:
- the radius = centre to corner, 40 m;
- flat-top, axial (q, r);
- **hex (0, 0) centred on the map's bottom-left corner** (`HexGrid::living_landz()`, `origin_m` = (0, 0)).
- If the origin is the hex's corner rather than its centre, `origin_m` shifts and nothing else changes.

**Built**: `river_hex_edges` walks each river every 0.1 hex size and records every change of hex as a crossing, with:
- the two hexes;
- the river;
- the catchment at the crossing (the flow accumulation's 3×3 maximum, capped at the mouth);
- the discharge (the mouth's, scaled by the catchment ratio);
- the local slope along the river (±1 cell of arc);
- W(A) = 200 m·(A/10 km²)^0.3 (the témoin's construction law).

The edges are written columnar in `rivers_ll.json` (format 0.2.0). Each river carries its `hex_edges` count, and the
viz's filter menu shows the edges kept. A test, `hex_of_flat_top_axial`, places the six neighbours' centres (and
off-centre points) in their hexes.

**On the témoin**:
- **395 399 crossings** by the 9 786 rivers, **all between neighbouring hexes** (0 jumps);
- 18.1 edges per km of river; per river p50 21, p90 95, max 1 038.
- At the crossings:
  - catchment p10 / p50 / p90: 0.1 / 0.6 / 36 km² (the head extensions carry most of the length);
  - discharge: 0.00 / 0.01 / 0.5 m³/s;
  - slope: 0.007 / 0.070 / 0.361 (2.0 % negative: bilinear noise and lake flats);
  - W(A): 55 / 87 / 295 m.
- **What a selection keeps**:
  - length ≥ 5 km: 1 036 rivers, 177 118 edges;
  - length ≥ 10 km: 275 rivers, 83 317 edges;
  - catchment ≥ 100 km²: 1 653 rivers, 192 843 edges.
- **Cost: 0.15 s per world.** `rivers_ll.json` grows from 5.1 MB to **31.9 MB** with the edges.

## E — the export, the preview, the cost

- **`rivers_ll.json` on the témoin: 5.1 MB** (112 698 points) without the hex edges, **31.9 MB with them**, against
  `rivers.json`'s 35.6 MB.
- The viz layer and its filters are built: the toggle "Rivières LL" and the "filtres ▾" menu.
- **The cost: 0.43 s build + 0.15 s hex edges + the serialization, about 0.6 s per world**, against 249.8 s (+0.2 %).
  In the viz it is timed as `[HD timing] rivers_ll`.

## Predictions

**The reviewer's (hypotheses to check):**
- **L** (> 99 % of the smoothed length within the valley-floor tolerance; no ridge crossing): **refuted as measured**
  (10.7 % above the 1 m tolerance after the constraint; 1 ridge crossing). The instrument itself fails on the raw
  trace (17.3 %).
- **P**:
  - "mostly under 10 km²": **refuted** (92 % at 10–100 km²; the head threshold sits near 7 km²);
  - "more than half on planar slopes or cones": **refuted on the construction's walls** (3.9 % walls + 2.7 %
    edifices). Uncarved planar slopes were not measured.
- **Q.3** (Ymir does not know the hex grid): **held**. Read, not blind; the author then gave it.
- **Cost < 10 s**: **held** (~0.6 s with the hex edges).
- **Meta**: **held**.

**Mine** (`f146_predictions.md`):
- **P-L**:
  - "1–5 % above before the constraint": **refuted** (23.3 %);
  - "< 0.5 % after": **refuted** (10.7 %; the instrument);
  - ridge crossings ~0 after: **held** (1);
  - lateral deviation p50 0.15 / p99 0.5 cell: **refuted, smaller** (0.000 / 0.094);
  - "≥ 1 crossing between rivers": **held** (662).
- **P-P**: "mostly < 10 km²" **refuted**; "> half on walls / edifices" **refuted**.
- **P-A**: "long-and-small rivers exist" **refuted** (0).
- **P-cost < 5 s**: **held**.
- **Meta**: **held**.

## What round 2 needs (for the author)

1. **Confirm the hex grid's reading**: is the origin the centre of hex (0, 0), or its corner?
2. **The look in Living Landz**: does the smoothing suffice visually? It stays within 0.1 cell of the D8 trace.
3. **Decisions**: rule (a) against the bundles, or not; and the spillways that retrace a watercourse.
4. **Open on Ymir's side, to fix in round 2 if the author wants them gone**:
   - the X-crossings of the D8 network (574 raw, 662 smoothed);
   - the bed instrument, to replace by the cell-under-sample terrain or the lateral deviation.

## Limitations, stated

1. **One world, one seed** (the témoin).
2. **The bed instrument's tolerance** (1 m, bilinear) is too strict for one-cell floors (above).
3. **The location classes** of P use the construction's masks on S1, not the eroded world's slopes.
4. **The crossing instrument** does not separate X-crossings from possible torus-seam or join artefacts.

## State

**Uncommitted**, awaiting the go-ahead:
- `export/rivers_ll.rs` (new), `export/mod.rs`;
- `export/container.rs` (the `rivers_ll` layer);
- `export/hydro.rs` (the test);
- the viz `bridge/c1/hd.rs` and `ui/workspace.rs` (the layer and its filters);
- `f126_coast.rs` (`f146_rivers`);
- `docs/rivers_ll_format.md`, `docs/adr/accepted_defects.md`;
- ADR (F145's closing, Finding 146);
- this folder.
