# Finding 154 — lakes / gorges / falls resumed, round 1 of 4: the geology split into substratum and surface (the final grid bit-identical, the substratum available at the construction stage for 0.35 s); every lip on hard rock with no contrast downstream, so rock does not discriminate the témoin's falls; P4 (no light pass on the kept bodies' catchments) cuts G-pits to 2–18 but not to 0, so it is NOT retained; the breach's ramp anchored at the pit's spill removes ON's 292 below-sea cells and the sea's advance at lake 4

**Status: built (the split, the bench options), measured, nothing promoted. The gorge stays gated, its toggle
hidden. Nothing committed.**

**Raw outputs, in this folder:**
- `f154_s_head_raw.txt`, `f154_s_raw.txt` (S);
- `f154_hlbr_raw.txt` (H, L, Br);
- `f154_br_attr_raw.txt`, `f154_br_spill_raw.txt` (the two Br amendments);
- `checks_before.txt`, `checks_after.txt`.

**Process:**
- The predictions were written before any measurement (`f154_predictions.md`, with the non-blind items declared).
- The instruments, the φ values (DECISION) and the selection rule were declared before the benches
  (`f154_declared.md`).
- **Two Br amendments came after the run. They are declared there as NOT blind.**

## The author's criterion, recorded as given

- « Je veux réduire la taille des lacs avec l'âge […]. » (2026-09-30)
- « Un lac présent est un niveau de base […] ; un lac vidé ne l'est plus. » (2026-09-29)
- « Les lacs des bassins sous la mer sont des lacs présents. » (2026-09-30)
- « Chute ou gorge, les deux sont valides. » « […] faire apparaître des chutes d'eau (et donc de les tagger comme tel)
  là où c'est pertinent. Ça donnera de la diversité au monde. » (2026-10-02)
- « C'est bien d'avoir l'âge qui fixe jusqu'où la gorge a reculé. » (2026-10-02)

**The decisions in force** (spec v5, `etat_au_F144.md`, reread):
- 18 bodies;
- the retreat per lake, p ∈ {0, 0.25, 0.5}, default 0.5;
- map (a); m = 10, the « raidir » rule;
- the fall capped by the river's size; the minimal plain;
- the below-sea lakes are out of v1.

## Partie 0

- **F153 committed** (`aa8ca60`). No wrong behaviour was reported from the author's look at the menus and the
  inspector, so nothing was corrected first.
- **The gorge stays gated, its toggle hidden.**
- **Checks before** (HEAD `aa8ca60`): `cargo check` clean; lib 617; viz 31; **the guard 6 / 6, field and lakes =
  banc**.
- **Checks after**: `cargo check` clean; lib **619** (+2, the two permanent tests below); viz 31; **the guard 6 / 6 =
  banc** (C2 /10 col `a8d2d538d692c2f0`). The geology stage logs 1.9–2.0 s per world.

## S — the geology split into substratum and surface

**Built** (`geology/rocks.rs`, `geology/mod.rs`):
- **`rocks::build_substratum`** holds what the tectonics and the edifices give, on every cell (water included): craton,
  belt, basaltic or arc volcanic, rift fill, basement.
  - It reads a relief: the contacts' terrain snapping and the cones' relief test.
- **`rocks::apply_surface`** lays over it what the end of the chain gives: water, evaporites, loose deposits.
- `build_rocks` = the surface over the substratum.
- **`geology::tectonic_sources`** (the coarse history, its masks, the edifices) and
  **`TectonicSources::substratum`** give the substratum from C1 alone, without the final drainage. `build_geology`
  goes through them.

**The gate — bit-identical, HELD.**
- The témoin's final rock grid, FNV-1a 64, was **`dbd91290b6c6148c` at HEAD before the split**, and the same after it.
- The class counts are identical.
- **Permanent test** (rule 13), `the_rock_grid_is_the_surface_over_the_substratum`:
  - F152's monolithic per-cell order is restated, and both paths give it;
  - negative control: a lake-margin cell is loose deposits at the surface over a rock substratum; the substratum holds
    no surface class.

**The stage and its cost:**
- The substratum needs **C1's history and the edifices**: 0.27 s, the geology's own coarse pass.
- It needs **a relief to read**: at the construction stage, the construction's input S1.
- **Its build on S1 takes 0.35 s**, or **0.62 s with its sources** if the history is not shared with the end-of-chain
  geology stage.
- **The substratum read on S1 differs from the one read on the final field on 5.04 % of the final land cells.** The
  terrain snapping follows the relief the erosion changes.
  - For the 18 lips the S1, final and final-substratum classes agree (H below).

## H — falls by rock (×1, p 0.5, the 18 bodies)

**The lip** is the col. **The downstream samples** lie on the construction's outlet path at 1 km and 3 km. The class
is the substratum on S1. The φ values are DECISION (mine, declared for the tabulation; the author's to set):
- (i) hard 0.5, medium 0.25, soft 0;
- (ii) a soft lip 0; a lip harder than its downstream 0.5; otherwise 0.1;
- (iii) the current draw.

**The head fall is replayed** with the construction's arithmetic (`GorgeBody::head_for_phi`, H_cap(A) unchanged).
- **Self-check: the replay at each body's own φ gives its built head fall to 0.0000 m.**
- Permanent test: `p4_excludes_the_kept_bodies_catchments_and_the_head_replays`.

**The lips:**
- **All 18 lips are on HARD rock**: 17 basement, 1 basaltic volcanic (lake (1871, 3666)).
- **The class 1 km and 3 km downstream is the lip's own in every case** (the 3 km sample "ends" for 3 bodies).
- No lip is on medium or soft rock.
- The final rock grid at the lip is the surface's class for 7 of 18 (loose deposits 6, water 1). **This is why the
  split was needed**: the surface hides the lip's rock.

| rule | tagged / 18 | median (m) | max (m) | by lip class (tagged / lips) |
|---|---|---|---|---|
| (i) hardness | **16** | 74.0 | 239.9 | basement 15/17 · basaltic 1/1 |
| (ii) contrast | **14** | 34.0 | 127.2 | basement 13/17 · basaltic 1/1 |
| (iii) the current draw (reference) | **11** | 60.4 | 239.9 | basement 10/17 · basaltic 1/1 |

- **On the témoin, rock does not discriminate.**
  - Rule (i) is φ = 0.5 everywhere. It tags every body except two, whose D_g is ≤ 6 m: lake 13 (2439, 5070), bounded,
    D_g −5.0, and body 15 (2541, 5265).
  - **Rule (ii) finds no contrast**. Its weak φ of 0.1 still tags 14: H_f is 10 m, and D_g is 86–1 208 m on 14 of
    the 18 outlets. **The "cascades" of a weak φ are falls by the tag.**
  - The per-lake table is in `f154_hlbr_raw.txt`.
- **Why**: the rock grid's contacts are smoothed at ≈ 3 km (σ = 8 on the 1 024 grid). At these outlets the coarse
  masks give basement over tens of km. **A rock rule can only diversify the falls where a contact crosses an outlet,
  and on the témoin none does.**

## L — P4: the light pass not applied on the kept bodies' catchments

**E**, as declared: the cells whose D8 path on the light pass's input (its own `compute_flow`) reaches a kept body's
footprint, footprints included. After the light pass, E is restored (bench option `light_mode = 4`,
`gorge_catchment_mask`).
- **E is 28.2–28.3 % of the land** (7 618–7 636 km²).
- 4–5 of the 18 cols fall in E: they drain into another kept body.

| world | cand. | held / 14 | **G-pits** | G-levels | G-rim | G-ring | G-slope100 | G-drained (constr.) | G-tag | G-sea | θ | edge (mask) | **divide edge** | build |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ×1 p .5 | P0 | 6 | 43 190 | 7/18 | 5 | 99 | 1/14 (31.5°) | 0 (0) | 194 | +0 | 0.444 IN | 573 | 851 | 362 s* |
| | **P4** | 6 | **2** | 10/18 | 4 | 136 | 2/15 (36.2°) | 0 (0) | 147 | +0 | 0.443 IN | 364 | **1 309** | 153 s |
| ×1.4 p .5 | P0 | 3 | 6 461 | 2/9 | 1 | 103 | 1/5 (32.6°) | 9 781 (2 190) | 110 | **+112** | 0.440 IN | 500 | 813 | 145 s |
| | **P4** | 3 | **4** | 2/9 | 1 | 138 | 1/5 (32.3°) | 6 513 (2 190) | 78 | **+122** | 0.445 IN | 354 | **1 360** | 160 s |
| ×1.4 p 0 | P0 | 7 | 6 336 | — | 0 | 105 | 0 | 14 366 (2 318) | 107 | +0 | 0.442 IN | 523 | 851 | 147 s |
| | **P4** | 7 | **18** | — | 0 | 143 | 0 | 7 001 (2 318) | 77 | +0 | 0.445 IN | 372 | **1 423** | 152 s |

\* The first world of the run, cold; the other builds are 145–160 s.

**Common to all six worlds:**
- F38 never fires; canyons 0; 18 bodies kept.
- G-bodies is 11–12 / 18 at ×1, 11 / 18 at ×1.4 p .5 and 15 / 18 at ×1.4 p 0.
- Drained absent: 8 / 9 and 15 / 18.
- **P0 reproduces F144's numbers** (G-pits 43 190 / 6 461 / 6 336; edge 573 / 500 / 523; G-sea +112).
- **The coast**: P4 against P0 is +0.04 % / +0.00 % / +0.00 %, within ±2 %.
- **θ** (corridors excluded): IN everywhere with the declared tolerance [0.423, 0.449]. Under F144's strict CI
  [0.428, 0.444], P4's 0.445 would be OUT at ×1.4: the gate still decides on noise.
- **G-area** is non-increasing from ×1 to ×1.4 (p .5) for both: P0 1 168.6 → 70.3 km², P4 1 179.5 → 68.3 km².

**Rule 14** (km³ against S1, in all / inside E):

| world | P0 removal | P4 removal | P0 raise | P4 raise | the light pass's removal in E under P0 (construction − eroded) |
|---|---|---|---|---|---|
| ×1 p .5 | 4 624.4 / 942.6 | 4 411.9 / 730.1 | 0.71 / 0.157 | 0.62 / 0.060 | 212.8 |
| ×1.4 p .5 | 4 723.0 / 1 176.1 | 4 452.3 / 905.4 | 0.74 / 0.178 | 0.65 / 0.094 | 271.1 |
| ×1.4 p 0 | 4 744.2 / 1 196.6 | 4 467.1 / 919.5 | 0.74 / 0.182 | 0.66 / 0.100 | 277.5 |

- P4 removes 4.6–5.8 % less in all; the difference is the light pass's removal in E.
- Under P4, construction − eroded in E is 0.00 km³: nothing after the pass lowers E measurably.

**What P4 does:**
- **G-pits falls by > 99.7 %, to 2 / 4 / 18, not to 0.** C2 (F143), which froze the rings and the corridors as well,
  reached 0.
  - The remaining pits were **not located**.
  - **Hypothesis to check**: pits in body cells under a rim cell outside E that the light pass raised.
- **G-drained falls by a third to a half** (9 781 → 6 513; 14 366 → 7 001), still about 3 × the construction's
  (2 190 / 2 318).
  - **The light pass makes part of G-drained.** F144 had attributed it to the construction's cones (Z); the cones
    stay the larger share's source.
- G-levels at ×1 rises from 7 to 10 of 18. **The off bodies have their col outside E**: the light pass still lowers
  the outlet.
- G-tag falls (194 → 147, 110 → 78, 107 → 77).
- **But G-ring rises 31–38 %** (99 → 136, 103 → 138, 105 → 143), and the steepest outlet rises from 31.5° to 36.2°
  at ×1. The edge on the design mask falls (573 → 364).
- **The divide edge** (a > 28° step across E's boundary) is **1 309 / 1 360 / 1 423 under P4**, against 851 / 813 /
  851 under P0 (×1.5–1.7). It is born at 93–131 on the construction.
- **G-sea at ×1.4 p .5 is unchanged** (+112 → +122). It comes from the breach (Br below).

**The selection rule (declared)**: « P4 est retenu s'il tient G-pits sur les trois mondes, et si ses cellules de bord
restent sous le tiers de celles de P1 au F144. »
- The divide edge stays under the third: 1 309 < 1 645; 1 360 < 2 080; 1 423 < 2 122.
- **G-pits is 2 / 4 / 18, not 0, in all three worlds.**
- **→ P4 is NOT retained. Report; nothing further is built.**

## Br — the breach's ramp, a bench copy (no production change)

**As declared** (the anchor `z[ci]`, the outlet cell's level when the flood reaches the pit). The copy with the
change off is bit-identical to production.
- **On ON:**
  - **286 of the 292** below-sea land cells remain (all under lake 1000016);
  - 12 954 conditioned cells change; **no lake changes** (34 / 34);
  - the river segments go 13 477 → 13 456 (−0.16 %).
- **On the gorge P4 ×1.4 p .5:**
  - G-sea goes from +122 / −0 to **+117 / −6**;
  - G-drained on the conditioned field goes 6 583 → 6 571; on the eroded field it is 6 513 (read before the breach).
- **The declared copy did not test the brief's idea.** That is the first amendment.

**Amendment 1 (NOT blind), the attribution (`f154_br_attr`):**
- **All 292 cells are taken by ramps**: 102 ramps from a row of pits along lake 1000016's fringe (y ≈ 5 205–5 255,
  x ≈ 2 710–2 925).
- The pits' floors are 7.0–12.7 m and their spills on the pre-breach flood 13.5–16.9 m. The ramps run 65–154 steps
  at 0.113 m per step.
- The cells' original heights are 0.34 / 11.3 / 16.9 m above the sea (p0 / p50 / p100).
- **When the flood reaches each pit, its outlet z[ci] has already been lowered by the neighbouring pits' ramps.** It
  lies 0.0–0.5 m above the floor. So "z[ci]" is not the overflow level: the ramps stay nearly as deep (73 ramps,
  286 cells).

**Amendment 2 (NOT blind): the anchor at the pit's spill**, `max(height[nb], filled[nb])` on the pre-breach flood
(`f154_br_spill`). The copy with the change off is bit-identical to production.
- **On ON:**
  - **the 292 below-sea land cells → 0**;
  - **22 863 conditioned cells change, all RAISED** (the pits are filled to their spill instead of trenched to the
    sea), none lowered;
  - **no lake changes** (34 / 34, levels and areas within 1 m / 5 %);
  - the river segments go 13 477 → 13 578 (+0.75 %); 16 941 river cells are in one network only.
- **On the gorge P4 ×1.4 p .5:**
  - below-sea land cells 414 → **0**;
  - G-sea against ON's production set: +122 / −0 → **+0 / −292**. The −292 are ON's own 292, which the same breach
    also removes, so **against ON with the same breach, G-sea is exact**.
  - **The sea's advance at lake 4 disappears.**
- **The decision is round 2's.** Nothing in production changed. The spill copy was not timed (the change is one
  more comparison per pit; the declared copy took 119 s against production's 117 s in this bench).

## The cost (against the reference run_hd 249.8 s)

- **S**: the split adds no work at the end of the chain (bit-identical; the guard logs geology 1.9–2.0 s). **A
  substratum at the construction stage would add 0.35 s**, plus 0.27 s for its sources if the history is not shared
  (+0.1–0.25 %).
- **P4**: **+8.4–8.6 s per world** (17.1 s on the cold first world): E's flow routing and walk, timed alone (+3.4 %).
  The light pass still runs on the whole grid, so **P4 does not reduce the time**. The build times (153 / 160 / 152 s
  against 145 / 147 s) are single timings within the noise.
- **Br (spill)**: not timed. The declared copy took 119 s against production's 117 s in the bench.

## Predictions

**The reviewer's (hypotheses to check; judged here):**
- **S**:
  - "bit-identical": **held**;
  - "< 0.5 s": **held for the substratum alone** (0.35 s); 0.62 s with its sources if not shared.
- **H**:
  - "rule (i) tags ≥ 80 % of the lakes": **held** (16 / 18, 89 %);
  - "rule (ii) < half": **refuted** (14 / 18). No contrast exists, and the weak φ passes H_f.
- **L**:
  - "P4 holds G-pits on all three with fewer edge cells than P1": **refuted on G-pits** (2 / 4 / 18); held on the
    edge (1 309–1 423 < 4 936–6 367);
  - "G-sea at ×1.4 and G-drained still fail": **held**.
- **Br** ("it removes the 292 below-sea cells and the sea's advance at lake 4"):
  - **refuted as declared** (the anchor z[ci]: 286 remain; +117);
  - **held under amendment 2** (the anchor at the spill: 0 and +0).
- **Cost** ("P4 reduces the time"): **refuted** (+8.4 s).
- **Meta**: held. Several refuted.

**Mine (`f154_predictions.md`):**
- **S**:
  - P-S1 bit-identical: **held**;
  - P-S2 (< 1 s; the history 0.5–2 s more): **held** for the substratum (0.35 s); the history is 0.27 s, under my
    range;
  - P-S3 (< 3 % differ): **refuted** (5.04 %).
- **H**:
  - P-H1 (≥ 14 hard, ≤ 2 rift): **held** (18 hard, 0 rift);
  - P-H2 (the same class downstream for ≥ 14): **held** (all 18);
  - P-H3 (rule (i) ≥ 14, median above (iii)'s): **held** (16; 74.0 > 60.4);
  - P-H4 (≤ 3 contrasts, 4–12 tagged): **held on the contrasts** (0); **refuted on the count** (14);
  - P-H5 (11 tagged; the replay < 1 cm): **held**.
- **L**:
  - P-L1 (G-pits 0 everywhere): **refuted**;
  - P-L2: under a third of P1 **held**; ≥ 1.5 × P0 **held** (1.54 / 1.67 / 1.67); "100–800" **refuted** (1 309–1 423);
  - P-L3 (retained): **refuted**;
  - P-L4: G-levels fails at ×1 **held** (10 / 18); "G-ring near P0's" **refuted** (+31–38 %);
  - P-L5: G-drained fails **held**, "at about the construction's count" **refuted** (P4 halves P0's, 3 × the
    construction's); G-sea +100 to +125 **held** (+122);
  - P-L6 (E 5–15 % of the land; the removal falling by that share): **refuted** (E 28 %; the removal −4.6 to −5.8 %);
  - P-L7 (θ in tolerance): **held**;
  - P-L8 (no reduction, +3 to +15 s): **held** (+8.4–8.6 s).
- **Br**:
  - P-Br1 (the 292 → 0): **refuted as declared**, **held under amendment 2**;
  - P-Br2 (> 100 000 cells, ≥ 5 lakes, > 1 % of the segments): **refuted** (12 954 / 22 863 cells; 0 lakes; −0.16 % /
    +0.75 %). My "anticipated" mechanism was wrong for z[ci] (amendment 1).
  - P-Br3 (G-sea's inlet gone): **refuted as declared**, **held under amendment 2**.
- **Meta**: **held**.

## What the round says (no decision taken; nothing promoted)

1. **The geology has two layers.** The substratum is available from C1 and a relief, so at the construction stage,
   for 0.35 s. The final grid is unchanged.
2. **Falls by rock: on the témoin every lip is on hard basement, with the same rock 1 and 3 km downstream.**
   - A hardness rule (i) makes almost every outlet a fall (16 / 18).
   - A contrast rule (ii) finds no contrast. Its weak φ still tags 14, since H_f (10 m) is small against D_g.
   - **Rock gives no diversity here. It would need contacts crossing outlets, or a finer structure.**
   - Named for the author: if rock is to drive the falls, the weak φ must fall under H_f, or H_f must rise for weak
     contrasts.
3. **P4 is not retained** (G-pits 2–18).
   - It is the best light-pass candidate so far on G-pits (> 99.7 % cut) and G-drained (−33 to −51 %).
   - It costs +8.5 s.
   - It makes a divide edge of 1 309–1 423 cells and a G-ring of 136–143.
   - **The declared rule stops the round here: nothing further is built.**
4. **The breach's ramp anchored at the pit's spill removes every below-sea land cell** (ON's 292; the gorge's 414 and
   the inlet at lake 4). It changes no lake, raises 22 863 cells and changes the river segments by +0.75 %.
   - **The declared z[ci] anchor did not**: earlier ramps lower the outlets.
   - **The decision is round 2's.**

**Named for the author:**
- **Br, the spill anchor**: a production change to decide. The breach is a cached step, so an `ALGO_*` bump would
  follow, and the guard's lakes and field would be re-read (ON's conditioned field changes by 22 863 cells).
- **P4's remaining pits** (2–18): to locate before any P4 variant, for example E plus the ring of the kept bodies.
- **Falls by rock**: no rule discriminates on the témoin. Keep the draw (iii), or change H_f for weak contrasts.
- The weak "cascade" φ of rule (ii) is a fall by the current tag.

## Limitations, stated

1. **One world, one seed.** The rock at the 18 lips is the témoin's: another seed may put contacts on outlets.
2. **H's replay** holds every other body's raises as built; a different φ on one body could move another's law (bodies
   are built in sequence).
3. **P4's pits were not located.**
4. **Br's second anchor** came after the run (amendment 2, not blind). The declared one was an incorrect reading of
   "overflow level".
5. **Single timings** (±30 s noise); the first world was built cold.

## State

**Uncommitted**, awaiting the feu vert.
- `geology/rocks.rs`: `SubstratumInputs`, `build_substratum`, `SurfaceInputs`, `apply_surface`, `build_rocks`
  delegating, and the test.
- `geology/mod.rs`: `TectonicSources`, `tectonic_sources`, `TectonicSources::substratum`, `build_geology` through
  them.
- `tectonics_c1/valley_construction.rs`:
  - `GorgeBody::{outlet_path, outlet_s_m, outlet_law_m, slope_law}`, `GorgeBody::head_for_phi`;
  - `gorge_catchment_mask`;
  - `light_mode` 4 documented;
  - the test.
- `tectonics_c1/production_upscale.rs`: `light_mode` 4 (P4, a bench option).
- `tests/common/mod.rs`: `viz_hd_lakes_with`, `BreachFn` (`viz_hd_lakes_on` delegates, unchanged).
- `tests/f126_coast.rs`: `f154_hlbr`. `tests/f154_resume.rs`: `f154_s`, `f154_br_attr`, `f154_br_spill`.
- ADR Finding 154; this folder.
