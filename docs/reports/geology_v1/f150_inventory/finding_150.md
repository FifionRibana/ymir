# Finding 150 — geology v1, round 1 of 3: what Ymir already knows about the rock, and the catalogue (nothing built)

**Status: inventory and catalogue; no production code; nothing committed.**

**Files:**
- predictions: `f150_predictions.md`, written before any measurement;
- instruments: `f150_declared.md`;
- raw output: `f150_geology_raw.txt` (bench `crates/ymir-core/tests/f150_geology.rs`);
- the maps: `map_*.png`;
- the catalogue: `catalogue_geologie_v1.md`;
- checks: `checks_before.txt` and `checks_after.txt`.

## The author's criterion, recorded as given

- Usage : « Chutes + ressources. Ce qui est visible dans LL, on verra plus tard, car ce sera peut-être un
  sous-ensemble. Là, on va faire du zoning : ici, il y a de l'or possible ; ici, du fer, etc. »
- Détail : « On va lister ce qu'on intègre avant. »
- Relief : « C'est peut-être déjà le cas. Mais pour moi le type de sol est déjà existant (craton, meuble, etc.), donc
  la roche en découle. À voir. Mais non : on fait la carte, puis on voit après. »
- **The limit: 3 rounds.**
  1. Inventory and catalogue (this round).
  2. The author picks the list; the specification is written.
  3. The map, the zones, the gates, the author's look.
- The cost principle: seconds per world against 249.8 s.

## Part 0

- F149 committed (`9e710be`, not pushed).
- **Recorded**:
  - **the author on the Strahler ≥ 2 filter**: « Filtrage avec Strahler semble très bien » (714 rivers kept of 8 857
    on the témoin);
  - **the reviewer's remark for Living Landz**: the viz's filters combine with AND. Keeping a long first-order river
    needs an OR rule on Living Landz's side (Strahler ≥ 2, or length ≥ X).
- Checks before = F149's after-checks.
- **Checks after** (`checks_after.txt`):
  - `cargo check --workspace` clean;
  - lib 610 passed;
  - viz 30 passed;
  - **the guard 6 / 6, field and lakes = banc** (C2 /10 col `a8d2d538d692c2f0`).
  - Production code is untouched this round.

## I — the inventory

**The grep** (rules 11 / 11c), over `ymir-core/src` and `ymir-viz/src`. Files hit per term:

| term | files |
|---|---|
| litholog | 13 |
| C-3 / C-3b | 10 / 8 |
| craton | 56 |
| plate type | 62 |
| crustal thickness | 15 |
| rift | 81 |
| arc | 40 |
| orogen | 26 |
| volcanism | 23 |
| sediment | 16 |
| fracture | 10 |
| hotspot | 4 |
| inherited structure | 4 |
| basalt | 3 |
| soil | 2 |
| granite | 1 |
| crust age | 2 |
| meuble, regolith, substrat | **0** |

**The fields** (measured on the témoin: `build_world(C2 /10 col)`, C-2, C-3 and C-3b ON as in the benches and the
6/6 guard):

| field | type, unit, grid | computed from, stage | filled on the témoin | exported | read by |
|---|---|---|---|---|---|
| `C1State.s` | f64, crust thickness non-dim., 64² | advected by the C1 time loop | yes: p10 / p50 / p90 0.06 / 0.24 / 0.71, max 2.18; continental p50 0.39, oceanic 0.22 | no | the isostasy → production altitude (`production_upscale.rs:307`) |
| `C1State.age` | f64, non-dim., 64² | ridge-aligned at init, advected, set to 0 on rift cells | yes: p50 0.07, p90 7.0, max 156; 3 294 distinct values | no | isostasy; C-3's rift class; the debug labels |
| `C1State.plate_type` | continental / oceanic, 64² | Voronoi init; mutated by subduction / accretion / rifting | yes: 832 continental, 3 264 oceanic | no | isostasy; C-3; volcanism (arc seeds); labels |
| `C1State.plate_id` | u16, 64² | Voronoi init; merges and splits | yes, but **2 plates left of 8 at init** | no | boundary classification, kinematics |
| `C1State.cratonic_mask` | bool, 64² | BFS at init over continental plates; static | yes: **57 cells** (6.9 % of the continental) | no | the time loop (rigid transport); isostasy; labels |
| boundary classes (`classify_boundaries`) | enum, 64² | `plate_id` + kinematics | yes: convergent 154, divergent 142, transform 74 | no | volcanism (arc sites), C-3b, labels |
| debug labels (`CoarseTectonicLabels`) | bool / f32, 64² | derived from the above, read-only | yes: craton 57, rift 66, subduction upper **7**, slab 147, **collision 0**, divergent 142 | no | **the viz's « Tectonique » overlay only**, in `HdResult.tectonic` |
| C-3 lithology K | f32 multiplier, HD | `build_coarse_k` (rift = continental `age < 1` → ×10) + edifice discs (×3) | yes: rift-soft 4.9 % of the land, volcaniclastic 4.3 %, hard 90.8 % | no | the incision only |
| C-3b fracture density | f32 [0,1], 64² → HD | distance to convergent / transform contacts, decay 25 km | yes: land p10 / p50 / p90 0.01 / 0.04 / 0.43 | no | the incision (K = 1 + 6·density); H-1 |
| C-2 edifices / craters | list | tectonic setting (arc / hotspot / rift) | 15 edifices, 15 craters | **crater lakes only** (`lakes.json`) | the upscale (stamped), C-3, the crater lakes, the viz markers |
| H-1 permeability | f32, HD | C-3 classes + C-3b density | **off** in the benches | no | the water balance when on |
| surface fields | HD | climate, drainage | yes | yes: biome, temperature, precipitation, water_class, lake_mask, flow_accumulation | Living Landz |

**In the viz, C-2, C-3, C-3b and H-1 are opt-in boxes, unchecked at startup** (`workspace.rs:458–461`). The benches
and the 6/6 guard's states run C-2, C-3 and C-3b ON (`f123_viz_guard`, the workspace's literal with them checked).
**Which of the two is "the production world" is the author's to confirm.**

**The three empty panel fields** ("Épaisseur crustale", "Type de plaque", "Craton"):
- they are hard-coded `"—"` in `workspace.rs:1922–1924` (`// coarse — deferred`);
- **"Type de plaque" and "Craton" ARE computed for a C1 continent**: `plate_type` and `cratonic_mask` at 64², and
  both already in `HdResult.tectonic` (`continental`, `craton`) when the viz requests the labels (« Générer HD » does).
  **Computed but not wired**;
- **"Épaisseur crustale" (`s`) is computed** and consumed by the isostasy, but **never carried to the viz**. It is
  not in the labels.

**The author's « type de sol » (craton, meuble, …)**: **it does not exist in Ymir as a classed field.**
- The nearest are:
  - the C-3 classes, three: hard basement / rift-soft / volcaniclastic. They are an ERODIBILITY, not a map of
    soils;
  - `cratonic_mask`;
  - the biomes (surface vegetation).
- There is no "meuble" (loose) class, no regolith, no sediment or deposit field. Production erosion is
  detachment-limited, so nothing is deposited (C-3's module doc).

**The maps** (north up, sea greyed; the coarse fields at the témoin's framing):

| map | what it shows |
|---|---|
| `map_crust_thickness_s.png` | `s` |
| `map_age.png` | `age` |
| `map_plate_type.png` | `plate_type` |
| `map_plate_id.png` | `plate_id` (2 plates) |
| `map_craton.png` | the craton |
| `map_tectonic_settings.png` | the labels: subduction upper orange, slab blue, divergent cyan, rift teal, craton gold, continental grey |
| `map_c3_lithology_class.png` | the C-3 classes: rift teal, volcaniclastic violet, hard tan |
| `map_c3b_fracture_density.png` | the C-3b density |
| `map_d_dz_c3_on_minus_off.png` | part D |

They show the C-3 lesson again: **the craton and the rifts are 6.25 km blocks, the volcanic class perfect discs.**

## D — does the rock already change the relief?

C-3 ON (the témoin) against C-3 OFF (`lithology_off`). C-2 and C-3b stay ON in both: this is the viz box
« Lithologie (C-3) » alone.

| measure | value |
|---|---|
| \|Δz\| p50 / p90 / p99 / max (land, 11.3 M cells) | **0.00 / 0.00 / 148 / 1 040 m** |
| signed p1 / p99 | −148 / 0.0 m: **C-3 only lowers** |
| cells \|Δz\| > 10 m | 372 787 (**3.3 % of the land**), of which 163 012 > 100 m |
| where the \|Δz\| > 10 m cells are | **61.7 % on a volcanic footprint, 35.4 % on the rift-soft class**, 2.8 % within 5 km of either, 0.0 % elsewhere |
| land altitude p50 | 515.3 m against 520.8 m (**−1.05 %**) |
| R8 | 0.0451 against 0.0459 (−0.0008) |
| local σ 3×3 p50 | 5.58 against 5.65 m |

**The answer**: C-3 changes the relief **only inside its two soft classes**.
- It deepens the volcanic edifices' flanks and the rifts by up to ~1 km.
- **Nothing on the 91 % hard basement moves.**
- So the geology already shapes the relief, but only through erodibility on 9 % of the land, and through the
  isostasy (`s`, `age`, `plate_type`, `cratonic_mask` set the altitude before any erosion: not measured here).

## C — the catalogue

`catalogue_geologie_v1.md`:
- C1, seven candidate rock classes, with their hardness for the falls, their placing fields and their artefact risks;
- C2, the reviewer's 11 resources plus 5 proposed, with their real contexts (cited from memory, to check: Robb 2005,
  Evans 1993, Cox & Singer 1986, Clifford 1966) and the existing fields;
- C3, three forms of zoning.

- **By the catalogue's grading** (a judgement on the code, not a measurement): **9 of the 11 starting resources are
  zonable with existing fields or a field derived from them.** Coal and the carbonate-hosted Ag-Pb-Zn need a
  sedimentary-basin field Ymir does not have.
- **But on the témoin the collision and active-margin settings are almost absent** (0 collision cells, 7 upper-plate
  cells): tin, orogenic gold and metamorphic gems would have empty zones there.

## Predictions

**The reviewer's (hypotheses to check):**
- **I**:
  - "Craton and Type de plaque exist in the schema but are not computed for a C1 continent": **REFUTED**. Both are
    computed and filled (57 and 832 cells) and sit in `HdResult.tectonic`; the panel is hard-coded.
  - "No classed « type de sol »": **HELD**.
- **D, "C-3 changes the relief p50 by < 5 %, mostly at the rift / volcano contact"**: **HELD** on the size (−1.05 %).
  On the place: ON the soft classes themselves (97 %), not at their contacts.
- **C, "at least half of the starting resources zone with existing fields"**: **HELD** (9 of 11 by the catalogue's
  grading), with the témoin's missing belts as a caveat.
- **Meta**: HELD.

**Mine:**
- **P-I**:
  - "computed, not wired": **HELD**;
  - "no « type de sol »": **HELD**;
  - "the age is degenerate, little structure": **REFUTED** (p50 0.07, p90 7, max 156, 3 294 distinct values).
- **P-D**:
  - "p50 < 2 %": HELD (1.05 %);
  - "|Δz| p50 < 1 m, p90 < 50 m": HELD (0, 0);
  - "max several hundred metres": REFUTED (1 040 m);
  - "≥ 80 % on or within 5 km of the soft classes": HELD (100 %);
  - "R8 < 0.01": HELD.
- **P-C, "8 of 11"**: REFUTED by one (the catalogue grades gemstones partly zonable).
- **Meta**: HELD.

**Not foreseen by anyone**: the témoin ends with 2 plates and no collision belt.

## For the author (round 2)

1. **Pick the rock classes and the resources** from `catalogue_geologie_v1.md`.
2. **The zoning's form** (C3): a possibility raster per resource, polygons, or rock classes plus rules.
3. **Confirm which world is "production"**: the viz starts with C-2 / C-3 / C-3b unchecked, while the benches and the
   guard have them on.
4. Wiring the three panel fields is a small viz change (two of them are already in `HdResult`), for whenever you want
   it.
