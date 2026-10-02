# Finding 138 — specification v2 of the gorge's retreat by age, and the tables for the author's decisions (nothing built)

**Raw outputs, in this folder:**
- the tables bench `f138_t.txt`;
- the first, defective run, kept: `f138_t_first_run_defective.txt`;
- `checks_before.txt` and `checks_after.txt`.

The predictions were written before any measurement (`f138_predictions.md`). The instruments were declared
before the measurements (`f138_declared.md`), with two amendments after the first run, said so below.
**The specification is `spec_gorge_age_v2.md`.**

**World**: the témoin (viz state C2 /10 col), ON extended.
- 14 lakes with a D8 outlet and 9 with a spillway path: 23 descents.
- T2 also reads the delivered world.

## Partie 0

- **F137 committed** (`c0852b1`): the f137_q bench, the folder, spec v1 as a draft, Finding 137.
- Added to ADR F137:
  - the eight-point review of v1;
  - the author's decisions of 2026-10-02, as given: the retreat follows the existing age selector; the gorge's
    slope by the river's size; « les deux » falls, with the reviewer's answer to « Qu'est-ce qui serait le plus
    cohérent ? ».
- **Checks**: before = F137's after-checks (no code changed since, declared). After = `checks_after.txt`. F138
  changes no production code (one bench).

## The amendments to the declared instrument (after the first run)

1. **T1 / T4's law below** was declared as "the extended skeleton's `floor_m` on the path".
   - The true col is a cell of the input lake's own bed (F136-K7), where `floor_m` is the LAKE's base (≥ the
     level). Every gorge "fit in 0.00 km", and T4's head falls came out negative.
   - **Amended**: the law below = B + k·χ computed **along the path from the next base B**, with the world's
     accumulation clamped at the outlet area.
2. **The outlet area** of the spillway lakes was 0.1–3.8 km²: the D8 accumulation stops in the basin.
   - **Amended**: A = max(the accumulation on the path's first km, the lake's inflow).
3. **φ's hash**: FNV-1a's high bits barely moved with the lake id (first run: φ ≈ 0.43–0.46 on 12 of 14 lakes).
   - **Replaced by splitmix64**. The printout's label read "FNV": a leftover string. It is *corrected in F139* in
     `f138_t.txt` (marked) and in the bench's source.
4. **T2 and T3 are unchanged between the two runs.** T3 also prints the totals of the 14 D8 lakes alone, because
   T3's instrument does not fit the below-sea lakes (below).

## T1 — the gorge by m (23 descents)

`S_g = min(m·S_loi(A), tan 28°)`, `S_loi = k·A^−0.5`, L = today's level, B = the next base on the ON path.

| | m = 3 | m = 10 | m = 30 |
|---|---|---|---|
| S_g range | 0.2–6.6° | 0.7–20.9° | 2.1–28° |
| fit / **no room** | 4 / **19** | 13 / **10** | 21 / **2** |

- At m = 10, the 10 with no room and their shortage falls:
  - D8 lakes: 1 (153 m), 2 (159 m), 3 (129 m), 9 (14 m), 10 (581 m);
  - spillway lakes: 1000011 (26 m), 1000014 (40 m), 1000016 (72 m), 1000018 (155 m), 1000019 (210 m).
- The large outlets (A 240–1 440 km²) have S_loi 0.11–0.27°. At m = 10 their gorges are 1.1–3.3° and need
  2.8–36 km.
- Lake 10 (D = 898 m to lake 1000008 over 12.6 km) cannot fit even at m = 10 (needs 35.7 km).

## T2 — the max slope over 100 m on today's steps (drop > 50 m)

| world | steps > 50 m | > 28° / > 40° / > 60° | among those > 200 m |
|---|---|---|---|
| ON extended | 8 | 3 / **2** / 0 | 6: 3 / **2** / 0 |
| livré | 9 | 4 / **4** / 0 | 6: 4 / **4** / 0 |

- ON: lake 10 **47.1°**, lake 2 **43.6°**, lake 7 31.3°; the others 22–28°.
- livré: 55.6°, 55.6°, 53.8°, 45.5°.
- **Short walls exist** behind mean slopes of 8.5–15° (F137-Q3 read only means). None passes 60° over 100 m.

## T3 — lakes by r (the input bowls' hypsometry)

| r | 0 | 0.5 | 1 | 1.5 | 2 |
|---|---|---|---|---|---|
| total area, the 14 D8 lakes | **1 308.5** | 1 223.6 | **1 145.7** | 367.3 | 0 km² |
| total area, 22 lakes with a body | 2 823.8 | 2 738.9 | 2 661.0 | 1 215.8 | 0 km² |
| lakes ≥ 1 km² | 22 | 22 | 22 | 22 | 0 |

- **The ratio r = 0 / r = 1 is 1.14 (D8 lakes) and 1.06 (all)**: from L_in to the bed sill the lakes lose little
  (3–48 m of level).
- The loss is beyond r = 1: −68 % at r = 1.5 (D8), and every lake is empty at r = 2 by definition.
- Per lake:
  - lake 2: 57.8 → 35.5 → 4.2 → 0 km² (r = 0, 1, 1.5, 2);
  - lake 10: 287.1 → 249.0 → 58.0 → 0.
- **The below-sea lakes do not fit T3's instrument**: their final level (the merge's) can exceed their body's
  ring (1000011: 459.3 against L_in 325.8), so L(r) stays flat up to r = 1. 1000001 has no body. Named, open
  (spec § 7).
- **The selector age each r would propose** (DECISION options): (a) ×0.7 → 0, ×1 → 1, ×1.4 → 2;
  (b) ×0.7 → 0.5, ×1 → 1, ×1.4 → 1.5.

## T4 — falls at m = 10

φ is per lake (splitmix64 of seed and id): 0 with probability 1/3, otherwise uniform in [0.1, 0.5].
H_f = max(1.5·S_g·100 m, 10 m).

| | count > 0 | tagged | median | max |
|---|---|---|---|---|
| head falls | 17 | 11 | **26.2 m** | 199.5 m (lake 2) |
| shortage falls | 10 | 10 | **129.1 m** | 581.4 m (lake 10) |

- **No gorge cell passes H_f** (0 two-cell drops over H_f).
- Every H_f but one (1000001's 57.4 m) is at the 10 m floor: 1.5·S_g·100 m is 2–6 m at m = 10.

## Predictions

**Mine:**
- **P-T1**: "at m = 10 at least half have no room" **refuted** (10 of 23, 43 %); "m = 3: ≤ 3 fit" refuted (4);
  "m = 30: most fit" held (21 of 23).
- **P-T2**: "ON: ≥ 2 of 6 > 28°" held (3); "**at most 1 > 40°**" **refuted** (2); "livré ≥ 2 > 40°" held (4).
- **P-T3**: held (1.06 / 1.14, not > 1.5).
- **P-T4**: held (shortages 129 m against heads 26 m in median; no gorge cell over H_f).
- **Meta** ("at least one of my two disagreements wrong"): **held** (T2).

**The reviewer's:**
- **T1** ("at m = 10, ≥ 1/3 have no room"): **held** (43 %).
- **T2** ("≥ 2 of 6 > 40°"): **held** (2).
- **T3** ("r = 0 > 1.5 × r = 1"): **refuted** (1.06 / 1.14).
- **T4**: **held**.
- **Meta**: held.

## Limitations, stated

1. **T1 / T4 use today's levels** (r = 1) and B + k·χ along the ON path, with the world's accumulation. The built
   gorge would read the skeleton's own areas.
2. **T3's area** grows 8-connected cells under L from the body's lowest cell, on the input field. Sub-bowls under
   a level are joined only if connected.
3. **The below-sea lakes**: T3 does not fit them, and three have no spillway (open).
4. **φ's law** is a proposal; its draw is one per lake per world.
5. **The amendments** (law below, area, hash) came after the first run. The first run is kept as raw output.

## State

**Uncommitted**, awaiting the go-ahead: `f126_coast.rs` (`f138_t`), this folder (with `spec_gorge_age_v2.md`), and
ADR Finding 138. **No production code. Nothing built.**
