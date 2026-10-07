# Finding 155 — lakes / gorges / falls resumed, round 2 of 4: the spill-anchored breach is NOT confirmed blind (the stop rule fires on all six guard states: the lakes' footprints change), so it does not ship; the soft-lip falls rule is built gated and is exactly the current draw on the témoin; the light pass removed: the gorge holds G-pits and G-rim everywhere, but the world shows the construction's naked facets, θ leaves its CI and the cost falls by 22–32 s

**Status: measured; the falls rule built gated; nothing promoted; the breach unchanged. Nothing committed.**

**Raw outputs, in this folder:**
- `f155_br_raw.txt`, `f155_br_loc_raw.txt` (Br);
- `f155_f_raw.txt` (F);
- `f155_p_world_raw.txt`, `f155_p_gorge_raw.txt` (P);
- the six `f155_ombrage_*.png`;
- `checks_before.txt`, `checks_after.txt`.

**Process:**
- The predictions were written before any measurement (`f155_predictions.md`, with the non-blind items declared).
- The instruments and the stop rule were declared before the benches (`f155_declared.md`).
- **One amendment came after the Br run (the localisation). It is declared there as NOT blind.**

## The author's criterion and decisions (2026-10-07), recorded as given

- The breach's ramp anchored at the overflow level: « Confirmer sur les 5 autres états, puis production ».
- The falls: « Tirage actuel, sauf lèvre tendre = pas de chute ».
- The light pass: « Mesurer ce qu'elle apporte (option : la retirer) ».
- F154's criterion stands (the lakes shrink with age; a present lake is a base level; fall or gorge, both valid; the
  age sets how far the gorge has retreated).

## Partie 0

- **F154 committed** (`297972a`), with three notes added to its ADR entry ("Recorded at commit (F155)"):
  - P4 is not retained; it was also worse on the walls, G-ring and the cost;
  - rock does not discriminate the témoin's falls (17 of 18 lips in the undifferentiated basement), which is v1's
    limit;
  - Br's amendment 2 is NOT blind, to be confirmed here.
- **Checks before** = F154's after-checks: no file under `src` changed between them and this round's start.
- **Checks after**: `cargo check` clean; lib **620** (+1, `a_soft_lip_makes_no_fall`); viz 31; **the guard 6 / 6, field
  and lakes = banc** (C2 /10 col `a8d2d538d692c2f0`).

## Br — the breach: the blind confirmation fails, nothing ships

**F154's amendment 2, untouched** (the ramp anchored at `max(height[nb], filled[nb])`), against production on the
guard's six states, through the same run_hd tail. **Production's lake fingerprint equals the guard's stored hash in
all six**, so the bench is the guard.

| state | below-sea land cells: production → spill (new) | lakes | footprint (mask hash) | list hash | raised / lowered cells | river segments |
|---|---|---|---|---|---|---|
| livré | 10 375 → **6 294** (0 new) | 60 → 60 | **differs** | differs | 337 506 / 892 | 7 301 → 8 068 (+10.5 %) |
| A1+B2 | 365 → **137** (0) | 38 → 38 | **differs** | differs | 176 819 / 86 | 11 806 → 12 054 (+2.1 %) |
| C1 nue | 414 → 0 (0) | 27 → 27 | **differs** | differs | 17 380 / 0 | 15 054 → 15 218 (+1.1 %) |
| C2 /10 col (témoin) | 415 → 0 (0) | 26 → 26 | **differs** | differs | 27 302 / 0 | 16 226 → 16 389 (+1.0 %) |
| C2 /3 | 416 → 0 (0) | 26 → 26 | **differs** | differs | 57 140 / 0 | 15 327 → 15 567 (+1.6 %) |
| C2 /10 niveau mer | 272 → **2** (0) | 24 → 24 | **differs** | differs | 22 373 / 0 | 12 796 → 12 841 (+0.3 %) |

**The declared stop rule FIRES on all six states**: « si un état change un lac (compte ou empreinte) […] ne pas passer
en production. Rapporter et localiser. » The count never changes. **The footprint changes everywhere.** No new
below-sea land cell appears anywhere. **Nothing ships.**

**Localised** (`f155_br_loc`, livré, the témoin and C2 /10 niveau mer):
- **The témoin**: 4 lakes **lose** cells, none gains: 1000016 −553, 1000018 −23, 1000019 −14, 1000017 −3.
  - These are the cells the floor anchor's ramps carved below the sea, which the merged below-sea lakes covered.
  - The spill anchor leaves them land, so the lakes lose them.
  - The production ramps: anchors 7.0–22.2 m above the sea, 113–196 steps.
- **C2 /10 niveau mer**: 3 lakes lose 4–10 cells.
  - **2 below-sea cells remain**: one ramp anchored 23.7 m above the sea, 209 steps long. A spill high enough still
    reaches the sea at 0.113 m per step.
- **Livré** (no construction): 14 of the 60 listing lines differ.
  - **The ids are permuted**: production's 1000002 (endorheic, 173.2 km²) is the spill world's 1000006 (172.7 km²),
    and 1000005 / 1000006 swap the same way.
  - Beyond the permutation, the lakes lose cells (1000001 −1 604, 1000007 −1 046, …).
  - **6 294 below-sea cells remain.** Their ramps start from coastal pits whose spill lies only 0.50–2.73 m above
    the sea (p50 0.67 m), and run 13–56 steps. **Anchoring at the spill does not stop a ramp that starts near the
    sea.**
- **Read**: in the five constructed states the footprint change is the same object in every case. The lakes shrink
  by exactly the cells the old ramps dug under them. Whether that is a lake "changing" in the rule's sense is the
  author's to say. The rule as written fires.

**The cache, found by reading the code** (it was never exercised, since nothing ships):
- `ALGO_BREACH` alone would NOT suffice. The climate (`climate_key`) and the HD drainage bundle (`hd_drainage_key`)
  read the CONDITIONED field but are keyed on the ERODED key. A breach change therefore needs `ALGO_BREACH`,
  `ALGO_CLIMATE` and `ALGO_HD_DRAINAGE`; otherwise run_hd serves a stale climate, lakes and rivers.
- **The skeleton's own breach** (`valley_construction.rs:642`, `breach_monotone`) uses the same ramp. Changing the
  function would change the construction, hence the eroded field. A production change must keep it apart (named for
  the author).

**The defect stays in the queue** (F142, `etat_au_F144.md` defect 3). The spill anchor was tested and did not pass
the declared rule.

## F — the falls' rule: the current draw, unless the lip is soft (built, gated)

- **Built**:
  - `GorgeRetreat::soft_lip` (absent from the key when off);
  - `GorgeRetreat::v6` = v5 + `soft_lip`;
  - `skeleton_with_soft_lips`.
- **The soft cells** in `production_upscale` behind the gate:
  - F154's substratum, read on the construction's input, at hardness 0;
  - the rift mask and the edifices only (the craton and belt masks are hard and tested after, so they are passed
    empty).
- **Permanent test** (rule 13), `a_soft_lip_makes_no_fall`:
  - a soft col gives φ = 0;
  - negative controls: a hard col, and the rule off, keep the draw;
  - the key moves.
- **The témoin** (`f155_f`, v5 against v6, ×1 p .5):
  - **0 soft lips** of 18; every φ identical;
  - the 11 tagged falls identical;
  - the constructed field identical (`ec60bde885e04bf4`);
  - **the production path's field identical** (`a219359ba8317d66`, 11 / 11 falls).
- The soft substratum (rift fill) is 3.36 % of the land on S1; none of it lies at a kept lip.
- **The cost when the gate is on**: the soft mask (the tectonic labels plus the substratum) ≈ 0.5 s per world. The
  bench's 0.48 s also reran C1's history, which production does not need.

## P — what the light pass brings to the world

**Grep (rules 11 / 11c), the origin, quoted before any number:**
- **F121** (ADR L17955): « the hybrid at 8192²: the shape constructed along the pre-incision trunks, **the texture
  left to a light incision** ». The construction replaced the shipped incision « by ONE light pass at `f × k_time` ».
  Its measured reasons (`f121_hybrid/finding_121.md`):
  - « **C1 [bare]: neither too smooth nor combed, but STRIPED at the scale of the valleys** … the walls are planes
    whose orientation follows trunks »;
  - « **C2 lowers the anisotropy instead of raising it**: R8 terrain 0.0798 → 0.0608 »;
  - « the light incision breaks the flat » of the floors;
  - and the C1 design doc's reason 3: *"'looks invented' is a real failure mode"*.
- **But the light pass also carries A1+B2** (F109's absolute slope floor with F104's depression floor:
  `light.slope_floor_uk = Some(s_eq); light.depression_floor = true`). That is a drainage closure. At F102 it took the
  canyon class from 29.6 % to 3.4 % and recovered the coast (Δ +3).
- `light incision` 1 hit · `texture` 32 · `A1+B2` 39 · `k_time` 65 · `light_k_time_fraction` 0 (code only).

**The production world** (the témoin, gorge off), with and without the pass:

| | WITH (production) | WITHOUT |
|---|---|---|
| relief p50 (paired) | 515.4 m | **565.8 m (+9.8 %)** |
| σ local p50 / p90 | 5.57 / 21.46 m | 4.78 / 21.20 m (−14 %) |
| R8 terrain | 0.0451 | **0.0577 (+28 %)** |
| R8 network | 0.4086 | 0.5154 (+26 %) |
| comb tile R8 | 0.0883 | 0.0724 (−18 %) |
| θ (trunks, carved, 32 450 / 34 264 links) | 0.436 [0.428, 0.444] | **0.493 [0.491, 0.495] OUT** |
| canyons | 0 | 0 |
| river segments · rivers_ll · Strahler ≥ 2 | 16 226 · 10 246 · 952 | 15 054 · 9 978 · 965 |
| parallel bundles (F146) | 1 626 | **4 133 (× 2.5)** |
| lakes (area) | 26 (3 492.8 km²) | 27 (3 332.9 km²) |
| biomes changed | — | 2.97 % of the land |
| Δz without − with, p10 / p50 / p90 | — | −0.0 / +0.1 / +113.9 m |
| planar walls (28° ± 0.5°) | 9.51 % of the land | **18.24 % (× 1.9)** |
| sharp crests (both sides > 28°) | **2 739** | 613 |
| build (twice) | 154 / 155 s | **123 / 123 s (−31.5 s)** |

- **The light pass is the world's texture.** Without it the construction shows bare: the images show cones and planar
  facets meeting in straight seams (`f155_ombrage_{montagne,plateau,vallee}_{with,without}.png`, 512² cells, north up;
  the crops were chosen on the WITH world first).
- **The sharp crests and the comb tile go the other way from my prediction.** The pass's gullies make knife-edged
  spurs, which the bare construction does not have. Without the pass the R8 terrain rises and the comb tile falls.
- The planar walls double, θ goes to the construction's 0.493, and the parallel bundles multiply by 2.5.
- **Rule 14, on the gorge's ×1 world** (P0 against NP, below): the pass removes **1 034 km³** of 4 624 km³ (3 590 km³
  without it). The raise falls from 0.71 to 0.06 km³.

**The gorge without the pass** (v5; F144's gates, F154's code; P0 reproduces F144 in all three worlds):

| world | cand. | held / 14 | **G-pits** | G-levels | **G-rim** | G-ring | G-slope100 (max) | **G-drained** (constr.) | G-tag | G-sea | θ | coast vs P0 | G-area | build |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ×1 p .5 | P0 | 6 | 43 190 | 7/18 | 5 | 99 | 1/14 (31.5°) | 0 (0) | 194 | +0 | 0.444 IN | — | 1 168.6 | 143 s |
| | **NP** | 6 | **0** | **15/18** | **0** | 112 | 3/16 (49.1°) | 0 (0) | 82 | +28/−17 | **0.495 OUT** | **−32.8 %** | 1 162.9 | **120 s** |
| ×1.4 p .5 | P0 | 3 | 6 461 | 2/9 | 1 | 103 | 1/5 (32.6°) | 9 781 (2 190) | 110 | +112 | 0.440 IN | — | 70.3 | 144 s |
| | **NP** | 4 | **0** | 3/9 | **0** | 102 | 1/6 (40.5°) | **2 190 (2 190)** | 18 | +139/−17 | **0.495 OUT** | +0.7 % | 78.6 | **121 s** |
| ×1.4 p 0 | P0 | 7 | 6 336 | — | 0 | 105 | 0 | 14 366 (2 318) | 107 | +0 | 0.442 IN | — | 0.0 | 140 s |
| | **NP** | 6 | **0** | — | **0** | 89 | 0 | **2 318 (2 318)** | 18 | +28/−17 | **0.495 OUT** | +0.7 % | 0.0 | **121 s** |

- **Without the pass the gorge holds G-pits = 0 and G-rim = 0 in all three worlds.** G-levels rises to 15 / 18 at ×1,
  but is 3 / 9 at ×1.4 p .5.
- **G-drained falls exactly to the construction's count** (2 190 / 2 318): it is the cones' (F144-Z), and it fails.
- G-tag falls to the construction's (82 / 18 / 18).
- **G-sea is not readable for NP.** Its reference is ON's set with the pass. ON without the pass was not built, so
  +28 / −17 measures the pass, not the gorge.
- **θ is OUT in all three** (0.495, the construction's own).
- At ×1, the coast falls by 32.8 % and the steepest outlet reaches 49.1°.
- G-area is non-increasing for both. Canyons 0, F38 never fires.

**The two options, costed (no choice made):**
- **(a) no light pass when the gorge is on**: **−20 to −23 s per gorge world** (−8 to −9 % of 249.8 s).
  - The gorge holds G-pits and G-rim everywhere, and G-levels at ×1 15 / 18.
  - G-drained (the cones) and G-ring / G-slope100 (the head falls) still fail.
  - The price is the world's texture: θ 0.495 OUT, the naked facets, ×2.5 bundles, the coast −33 % at ×1.
- **(b) no light pass at all**: **−31.5 s per world** (−12.6 %).
  - The production world changes in all the ways in the table above: relief +9.8 %, R8 +28 %, θ out, bundles ×2.5,
    planar walls ×1.9, biomes 3 %.
  - The A1+B2 closure is gone (the canyons stay at 0 on the témoin).
  - **The F121 reason the pass was added (the construction "looks invented") comes back.** The images show it.

## Predictions

**The reviewer's (hypotheses to check; judged here):**
- Br, "the amendment holds on the 5 states, no lake changes, the below-sea cells fall to 0 wherever there were
  some": **refuted**. Every footprint changes; livré keeps 6 294, A1+B2 137 and niveau mer 2.
- P grep, "added to soften the construction's geometry, not for a hydrological reason": **held in part**. F121's
  reason is the texture and the anisotropy, but the pass also carries the A1+B2 drainage closure.
- P without the pass:
  - relief p50 < 2 %: **refuted** (+9.8 %);
  - R8 terrain > +10 %: **held** (+28 %);
  - θ toward 0.5: **held** (0.493);
  - the cost falls by > 5 s: **held** (−31.5 s).
- The gorge without the pass:
  - G-pits and G-rim hold on all three: **held**;
  - G-levels holds: **refuted** (15/18, 3/9);
  - G-drained still fails at ×1.4: **held**.
- F, "exactly the current draw on the témoin": **held**.
- Meta: **held**.

**Mine (`f155_predictions.md`):**
- **Br**:
  - P-Br1 (the footprint identical, the rule does not fire): **refuted**;
  - P-Br2 (the list changes): **held**;
  - P-Br3 (below-sea → 0 everywhere): **refuted** (livré, A1+B2, niveau mer); "none new" held;
  - P-Br4 (5 000–60 000 raised, 0 lowered, < 2 % segments): **refuted for livré and A1+B2** (337 506 / 176 819
    raised, 892 / 86 lowered, +10.5 % / +2.1 %); **held for the four C states**;
  - P-Br5 (the three ALGO bumps): **held by reading**, not exercised.
- **F**: P-F1 **held**.
- **P**:
  - P-P1 **held**;
  - P-P2:
    - relief, R8, σ, θ, canyons, segments, Strahler, bundles, biomes and the cost: **held**;
    - "lakes +10 to +40 %": **refuted** (26 → 27, area −4.6 %);
  - P-P3:
    - planar walls ≥ 3×: **refuted** (×1.9);
    - "the sharp crests rise": **refuted** (2 739 → 613);
    - "the comb tile's R8 rises": **refuted** (0.0883 → 0.0724);
  - P-P4:
    - G-pits 0, G-rim, G-levels failing at ×1, G-drained ≈ the construction's, G-ring failing: **held**;
    - "G-sea holds": **refuted**, and not readable.
- **Meta**: **held**.

## What the round says (no decision taken)

1. **The spill-anchored breach does not pass the declared rule.**
   - In the five constructed states, the lakes lose exactly the cells the old ramps dug under them.
   - In livré the ids are permuted, and 6 294 coastal cells stay below the sea (pits whose spill lies < 3 m above it).
   - Nothing ships. The cache needs three bumps, not one. The skeleton's breach must stay apart.
2. **The falls' rule is built gated**, and on the témoin it is the current draw exactly (no soft lip).
3. **The light pass is the world's texture and its drainage closure (A1+B2).**
   - Without it, the gorge holds G-pits and G-rim.
   - But the construction's facets show, θ leaves its CI, the bundles multiply by 2.5 and the coast drops at ×1.
   - It costs 22–32 s per world.

**Named for the author:**
- **Br**: is "a lake loses the cells the old ramps dug under it" a change in the rule's sense? A second rule, for
  example no lake gains a cell and no id changes, would pass the four C states and fail livré on its ids. A ramp
  limit near the sea (livré's coastal pits) would be another round.
- **P**: option (a) or (b), or neither. With neither, the gorge's G-pits stay at the light pass's mercy (F142–F154).
- **G-drained (the cones)** is now the gorge's remaining construction-born failure, with G-ring / G-slope100 (the
  head falls).

## Limitations, stated

1. **One world, one seed** for P and F. Br used the guard's six states on one seed.
2. **G-sea for NP** has no matching reference (ON without the pass was not built).
3. **The images** are three 25 km crops: the reader's eye, not a measure.
4. **Single timings**, ±30 s, but each P world was built twice: 154 / 155 against 123 / 123 s.
5. **The localisation** came after the run (amendment, not blind). It covers 3 of the 6 states.

## State

**Uncommitted**, awaiting the feu vert.
- `tectonics_c1/valley_construction.rs`:
  - `GorgeRetreat::soft_lip`, `GorgeRetreat::v6`;
  - `skeleton_with_soft_lips` (`skeleton` and `skeleton_patched` delegate to one `skeleton_impl`);
  - φ = 0 on a soft col;
  - the test.
- `tectonics_c1/production_upscale.rs`: the soft cells behind the gate.
- `tests/f155_round2.rs` (`f155_br`), `tests/f155_soft.rs` (`f155_f`, `f155_br_loc`), `tests/f126_coast.rs`
  (`f155_p_gorge`, `f155_p_world`).
- The breach, the cache versions and the guard are unchanged.
- ADR Finding 155; this folder.
