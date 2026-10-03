# Finding 142 — the gorge's retreat v4 (14 D8 lakes, the minimal plain), re-tabulated and gated: the construction holds (1:1 bodies, no pits, no new sea, F38 never fires, a negligible deposit), but the light pass erodes the tagged head falls and reopens the lakes and plains (G-pits, G-rim, G-levels fail); θ falls with the decided steepness; four gates are not located

**Status: built gated (off by default), measured, NOT promoted. The viz toggle stays hidden. Nothing committed beyond Partie 0.1.**

**Raw outputs, in this folder:**
- `f142_r_raw.txt` (the re-tabulated tables);
- `f142_g_raw.txt` (the gates) and `f142_g_first_run_defective.txt` (stopped, G-pits amended);
- `f142_theta_raw.txt`, `f142_sea_raw.txt` (the attributions);
- `checks_before.txt`, `checks_after.txt`.

The predictions were written before any measurement (`f142_predictions.md`). The instruments were declared before
each bench, with three amendments, each declared before its run (`f142_declared.md`). The specification is
**`spec_gorge_age_v4.md`**.

## Partie 0

1. **F141 committed** (`f4a6a56`): `f141_attr`, `f141_b4`, `f141_attribution/`, ADR Finding 141.
2. **ADR F140 corrected in place, flagged "[CORRECTED in F142, after F141.]"** in three places. The same note heads
   "Where it comes from" in `finding_140.md`.
   - "lake 4: 1 571 against 298 km²" was an index pairing. For the true lake 4 (body 3), A agrees under all three
     definitions (300.7 / 298.1 / 305.2 km²).
   - "The cause is deviation 1" is replaced by F141's α / β / γ.
   - The reviewer took the figure up without checking it.
   - The method rule (re-tabulate before building) stands.
3. **The defect, named and queued, out of this work**: the breach's ramp starts at each pit's floor and descends
   0.113 m per cell towards the flood's root, without reading the spill (`flow.rs:1066–1075`).
   - In production it leaves 292 land cells under the sea, all under lake 1000016.
   - **Finding 38's invariant holds by covering.**
   - Recorded in ADR Finding 142's queue.
4. **The author's decisions (2026-10-02), recorded as given** (ADR F142, spec v4 § 0):
   - « Seulement les 14 lacs D8 (sous la mer inchangés) »;
   - the minimal plain: the question « Les deux étant possibles, qu'est-ce qui les discerne ? », the reviewer's answer
     (two stages of one history), and « Ok go ».
5. **`spec_gorge_age_v4.md`** = v3 plus:
   - the scope (a construction-side proxy, deviation 6);
   - the minimal plain;
   - the disposition of F141's eight differences, with their reasons;
   - F141's rule (pairing by footprint, never by index).
6. **The viz toggle stays HIDDEN** (`GORGE_TOGGLE_VISIBLE = false`). It now builds v4.
7. **Checks**:
   - before = F141's after-checks (`f4a6a56`; no file under `src` changed between them and this round's start);
   - after = `checks_after.txt`.

## R — the tables re-tabulated as built, the input alone (`f142_r_raw.txt`)

**The pairing** (by footprint, F141's rule): the 14 D8 lakes go 1:1 to 14 distinct bodies. The body's share is
60–98 % (lakes 2 and 7: 61 % and 60 %). The lake's share is 99–100 %.

### T3 — level and area by r

| | r 0 | 0.5 | 1 | 1.5 | 2 |
|---|---|---|---|---|---|
| **as built** (km², the 14) | 1 308.6 | **1 306.5** | **1 304.5** | 403.4 | 0 |
| F138-T3 | 1 308.5 | 1 223.6 | 1 145.7 | 367.3 | 0 |

- **As built, L_bed is within 0.0–2.0 m of L_in, and 2.8–46.6 m ABOVE today's level**:
  - lake 2 +46.6;
  - lakes 7 and 11 +28.3;
  - lake 10 +26.6;
  - lake 19 +2.8.
- The construction's minimax from the lowest cell to the outflow, on the input, does not find the bed sill that
  today's level shows. **Deviation 2 is not "the same object"**: F138-T3 took L_bed = today's level.
- **So r ≤ 1 barely shrinks the lakes** (−0.3 %, against F138's −12 %). The retreat lives in r > 1.

### R — the lakes left and the total area by r_world and p (A as built; A_ref 418.7)

The as-built median of A over the 14 is 425.1 km² (+1.5 % on A_ref). A_ref stays the declared constant.

| p | r_world 0 | 0.5 | 1 | 1.5 | 2 |
|---|---|---|---|---|---|
| 0 | 14, 1 308.6 | 14, 1 306.5 | 14, 1 304.5 | 14, 403.4 | 0, 0 |
| 0.25 | 14, 1 308.6 | 14, 1 306.3 | 14, 986.6 | 13, 164.9 | 3, 5.9 |
| 0.5 | 14, 1 308.6 | 14, 1 306.2 | 14, **653.4** | 9, 121.1 | 4, 22.1 |
| F139-R, p 0.5 | 14, 1 308.5 | 14, 1 202.7 | 14, 521.3 | 9, 103.0 | 4, 19.0 |

**At ×1, p = 0.5 the 14 lose 50 %** (653.4 km²), against F139's 55 % (521.3). The author's −55 % becomes −50 %.

**The emptying ages, p = 0.5** (map (a)):

| lake | 1 | 2 | 3 | 10 | 11 | 12 | 13 | 15 | 19 |
|---|---|---|---|---|---|---|---|---|---|
| F139-R | ×1.05 | never | ×1.22 | ×1.17 | ×1.14 | ×1.39 | ×1.12 | ×1.40 | ×1.03 |
| as built | ×1.04 | **×1.38** | **×1.27** | ×1.17 | ×1.13 | ×1.39 | ×1.12 | **never (2.01)** | **×1.21** |
| Δ | −0.01 | enters | **+0.05** | 0 | −0.01 | 0 | 0 | leaves | **+0.18** |

- **8 lakes empty within the range** (×1.04–×1.39, spread 50 %), as in F139, but not the same 8: lake 2 enters and
  lake 15 leaves.
- 5 of F139's 8 move by ≤ 0.01. **Three move by more than 0.03**: lake 3 (+0.05), lake 19 (+0.18; A as built is half
  of F139's) and lake 15 (to never).
- Never within the range: lakes 4, 7, 8, 9, 15 and 20.
- p = 0.25: 8 empty, ×1.19–×1.40 (29 %). p = 0: all 14 at ×1.40 (v2's defect, unchanged).

### T4 — the head falls as built (steepen, m = 10; the construction's own values)

| | head falls > 0 | tagged (> H_f) | median | max | shortage falls |
|---|---|---|---|---|---|
| r_world 1, p 0 (r_lake 1: F139-C4's level) | 14 | **10** | 34.2 m | 463.6 m (lake 10) | 0 |
| r_world 1, p 0.5 | 14 | **9** | 34.2 m | 418.0 m (lake 10) | 0 |

- The four lakes at φ = 0 (2, 9, 11, 20) carry the one-cell step S·s₁ (1.6–4.7 m), below H_f.
- **φ by the construction's draw (deviation 3)**: 10 of the 14 lakes have φ > 0, as in F139-C4, but not the same 10.
  - As built: 1, 3, 4, 7, 8, 10, 12, 13, 15, 19.
  - F139: 2, 3, 4, 7, 9, 11, 12, 13, 15, 20.
  - So the falls change place. Lake 10's head fall goes 0 → 464 m; lake 2's goes 199 → 4 m.
- **At p = 0.5, lake 13** (r 1.53) has D_g = −6 m: its level L(r) = 129.5 m stands under the law below at the first
  fit, so there is no fall.
- **No shortage fall** in either (steepen: F139-C3).

## B — what is built (gated: `gorge_retreat: None` by default; v4 = `GorgeRetreat::v4`)

All in `valley_construction.rs` unless stated.
- **`GorgeRetreat`** gains `d8_scope` and `plain` (`#[serde(default)]`, serialised, so v4 moves the key against v3).
  `v4(r, p)` = v3's values with both on. v3 is unchanged (F140's setting stays reproducible).
- **The scope** in `gorge_bodies`: under `d8_scope`, a body is kept only if it is an input lake whose ring touches
  no other body's closed-depression cell and no cell under the sea.
  - The other bodies are not created. Their cells keep the extended lake base: no retreat, rim or gorge.
  - `GorgeBody` gains `input_lake`, `touches_depression` and `touches_below_sea`, as diagnostics.
- **The minimal plain**: `pub fn gorge_plain(field, sk, ss, min_r) -> Vec<GorgePlain>`.
  - A priority flood restricted to each kept body's footprint with r_lake > `min_r`, seeded on its ring, in
    normalised units.
  - A cell whose fill level is the footprint's spill is the lake's bed and is kept. A higher fill level is a
    residual hollow, and is raised to it.
  - **The construction's first term that RAISES the terrain** (a deposit). `carve` still only lowers.
  - **Placement**: `production_upscale.rs`, right after `carve` and before the stream-power block (the light pass),
    under `vc.gorge_retreat.is_some_and(|g| g.plain)` with `min_r = 1`.
- **The viz** builds v4. The toggle stays hidden (`GORGE_TOGGLE_VISIBLE = false`).
- **Permanent tests** (rule 13, negative controls first; lib 593 → 595):
  - `the_minimal_plain_fills_a_drained_bowls_residual_hollows_and_only_raises`:
    - control: a dimple in a drained bowl is a hollow without the plain;
    - with it: only up, only inside the footprint, by the reported volume, idempotent;
    - nothing when r ≤ 1.
  - `v4_keeps_only_the_input_lakes_standing_alone`:
    - control: v3 treats a body merged with a below-sea basin;
    - v4 keeps only the lone bowl; a dropped body's cells keep the extended base bit for bit;
    - v4's key differs from v3's.

## G — the gates (`f142_g_raw.txt`; the first run, stopped and amended: `f142_g_first_run_defective.txt`)

**Each world ran isolated. No world panicked, and F38's invariant fired in none (8 of 8).**

**G-inert HOLDS**: with the gate off, the after-checks' guard reads 6 / 6, field and lakes (C2 /10 col = `a8d2d538d692c2f0`). Lib 595, viz 30, `cargo check` clean (`checks_after.txt`).

| world (age, p) | ×0.7, .5 | ×0.85, .5 | ×1, .5 | ×1.2, .5 | ×1.4, .5 | ×1, 0 | ×1, .25 | ×1.4, 0 |
|---|---|---|---|---|---|---|---|---|
| **G-bodies** (1:1, no split) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **G-pits** at (a) (construction + plain) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **G-pits** at (b) = (c) = (d) | — | — | **43 078** (188 m) | **8 788** (113 m) | **9 628** (122 m) | — | **80 424** (257 m) | **9 215** (122 m) |
| **G-sea** (ON's 292 exactly) | ✓ | ✓ | ✓ | ✓ | **+112** | ✓ | ✓ | ✓ |
| **F38** fires | no | no | no | no | no | no | no | no |
| **G-deposit** (km³) | 0 | 0 | 0.060 | 0.066 | 0.118 | 0 | 0.0001 | 0.133 |
| **G-levels** (not drained, ±1 m) | **2/18** | **7/18** | **7/18** | **3/13** | **2/9** | **8/18** | **9/18** | — |
| drained: absent at the end | — | — | — | 3/5 | 7/9 | — | — | **14/18** |
| **G-rim**, construction / eroded world | 0 / **12** of 18 | 0 / **12** | 0 / **5** of 9 | 0 / **2** of 3 | 0 / **1** of 2 | 0 / **12** | 0 / **5** | — |
| **G-ring** (> 28° beside a clamped cell) | **201** | **214** | **162** | **127** | **113** | **205** | **156** | **95** |
| **G-slope100** (> 28°; ON: 4 of 14) | 1 | 1 | 2 | 2 | 1 | 2 | 2 | 0 (no path) |
| **G-drained** (< L_floor − 1 m) | 0 | 0 | 0 | **53** | **2 223** | 0 | 0 | **2 396** |
| **G-tag**: tagged head falls | 12 | 12 | **11** | 8 | 7 | 12 | 12 | 7 |
| G-tag: 2-cell drops > H_f beyond the lip | **227** | **225** | **145** | **106** | **58** | **207** | **183** | **51** |
| the 14 D8 bodies' final lakes (km²) | 1 246.3 | 1 256.9 | 658.1 | 141.7 | 67.9 | 1 251.8 | 972.1 | 0 |
| F142-R's R (p 0.5; km²) | 1 308.6 | 1 306.2 | 653.4 | 121.1 | 22.1 | — | — | — |
| cost added (s, on ~115–160 s) | −3.8 | +9.5 | −3.1 | +2.2 | +0.9 | +7.5 | +0.9 | +1.6 |

"—": not applicable. No body has r_lake > 1, so there is no plain and nothing to measure at (b)–(c), or there is
no undrained lake.

**Shortage falls: 0 in every world** (steepen). **Canyons: 0** (×1 p .5, ×1.4 p 0, ON, OFF).

### The failures, attributed (rule 12) where measured

1. **G-pits fails from the light pass on**, never at the construction.
   - The pipeline's (a) equals the bench's construction + plain within 0.01 m (9 640 cells differ by ≤ 1 cm, in
     bodies 13–15), so **the production wiring of the plain holds**.
   - At (b), the ring's minimum falls under L(r) in step with the body's head fall:

     | lake | head fall (F142-R T4, p 0.5) | ring minimum at (b) − L(r) | pit cells at (b) |
     |---|---|---|---|
     | 10 | 418 m | **−283.6 m** | 41 236 (to 188 m) |
     | 3 | 406 m | **−169.4 m** | 43 (to 45.6 m) |
     | 19 | 60 m | −18.8 m | 0 |
     | 1 | 37 m | −7.6 m | 1 773 (to 0.7 m) |
     | 12 | 14 m | −1.8 m | 0 |
     | 2, 11 (φ = 0) and 13 (head 0.1 m) | ≤ 5 m | −5.6, 0, 0 | 1, 16, 8 |

   - **The light pass (one stream-power pass, A1 + B2) erodes the tagged head falls.** It does not know the tag. It
     incises the one-cell fall and cuts a canyon back from the lip into the lake or the plain (lake 10: 284 m). That
     opens closed hollows at many levels.
   - (b) = (c) = (d): the droplets and the bathymetry add nothing.
   - **Read in rank over 9 bodies, not tested separately** (no run with the falls removed).
2. **G-levels and G-rim (on the eroded world)** fail at every age, ×0.7 included. The construction's rim holds
   (0 cuts > 10 m at every age). The eroded world's rim is cut on 12 of 18 bodies at r_world 0, where every lake
   keeps L_in and carries its head fall (12 tagged, to 414 m). **The same light-pass erosion of the lip is the
   reading; not measured body by body.**
   - **Lake 13** (L_floor 13.9 m) and **the extra body 15**, drained at ×1.2 and ×1.4, are covered at 60 % and 73 %
     by a final lake at 116.7–116.9 m, 1000016's merged level (115.3 m in ON). **Once lowered below ~116 m, lake
     13's bowl joins the below-sea basin's merged water.** Its body touches no depression cell on the input.
3. **G-ring fails at every age** (95–214 cells, up to 83°). The listed cells sit beside body 0's col, (3382–3386,
   2598–2602): the clamp holds the ring at L(r) next to the gorge's cut. **Not attributed further.**
4. **G-drained fails at ×1.2 and ×1.4** (53; 2 223; 2 396 cells upstream of drained lakes, under L_floor − 1 m on
   the construction). **Not located.** F141-Z's mechanism (a neighbouring line's cones) is a candidate, not a
   finding.
5. **G-tag**: 51–227 two-cell drops > H_f along the gorge beyond the lip. **Not located.** At ×1, p 0.5 the world
   carries 11 tagged head falls; F142-R's T4 counts 9 on the 14 for the same setting, so 2 lie on the 4 extra bodies
   (by difference).
6. **G-area is not monotone** at p 0.5: 1 246.3 → 1 256.9 km² from ×0.7 to ×0.85 (+10.6), then 658.1 / 141.7 /
   67.9. At ×1.2 and ×1.4 the final areas exceed R's (141.7 vs 121.1; 67.9 vs 22.1) because of the lakes that refill
   (lake 13 at 116.7 m).
7. **G-sea fails at ×1.4, p 0.5 only** (+112 cells), attributed by `f142_sea` (below). F38 does not fire there.

### θ (non-regression) fails, and is attributed (`f142_theta_raw.txt`)

| θ on the carved links, construction mask | value [95 % CI] (links) |
|---|---|
| OFF ×1 | 0.493 [0.491, 0.495] (34 264) |
| ON ×1 | 0.493 [0.491, 0.495] (36 943) |
| **GORGE v4 ×1, p 0.5** | **0.264** [0.253, 0.274] (36 746) |
| GORGE v4 ×1.4, p 0 | **0.383** (f142_g) |
| GORGE ×1 without the 1 278 links touching a cell whose base the gorge RAISES | **0.478** [0.474, 0.483] |
| GORGE ×1 without the links the rebase LOWERS | 0.225 [0.212, 0.238] |
| GORGE ×1 without both | **0.493** [0.491, 0.496] (= ON on the same links: 0.493) |
| GORGE ×1 on the raised links only / the lowered only | 0.482 [0.451, 0.514] / 0.435 [0.419, 0.450] |

- Each class keeps a concave θ on its own. **The pooled fit falls because the gorge paths are a steeper population
  at the largest areas.** S = max(m·S_loi, S_req) with m = 10 is ten times the law's slope for the same A, and that
  flattens log S against log A across the mix.
- Removing those 1 278 links (3.5 %) brings θ back to 0.478, just under OFF's CI. The rest of the gap is the
  lowered catchments'.
- **The θ gate measures the decided steepness (m = 10), not a defect of the law.** It still fails as written.

### Rule 14 and rule 18 (×1; the extras)

| world | lakes (non-crater km²) | removal S1 − world, land (km³) | inside the kept bodies | raise world − S1 (km³) | the plain (km³) |
|---|---|---|---|---|---|
| OFF ×1 | 26 (3 491.6) | 5 321.2 | — | 0.64 | — |
| ON ×1 | 34 (4 445.0) | 4 515.9 | — | 0.66 | — |
| GORGE v4 ×1, p .5 | 34 (3 961.7) | 4 634.3 | 18.7 | 0.71 (0.08 inside) | 0.060 |
| GORGE v4 ×1.4, p 0 | 27 (3 623.4) | 4 798.4 | 119.3 | 0.77 (0.14 inside) | 0.133 |

- **The construction's removal against ON at the same age** (G-deposit's control) is 150.4 km³ at ×1, p .5 and
  434.8 km³ at ×1.4, p .5. Most of it is the rebase: every tributary of a lowered lake is laid from a lower base.
- **The deposit is 0.03–0.04 % of that, and 0.11 % of the removal inside the kept bodies at ×1.4, p 0.** The deepest
  fill is 37 m (lake 1 at ×1).
- **The coast**: spurs near a coastal wall 0.0142 /km in ON and GORGE ×1 (unchanged), 0.0050 at GORGE ×1.4, p 0.
  The age also changes k there; ON was not run at ×1.4, so this is not attributed.
- **The cost**: +0.9 to +9.5 s and −3.1 / −3.8 s on constructions of 114–160 s. The median added is +1.25 s. The
  timing noise (± ~10 s, the bench sharing the machine) does not resolve "< 5 s".

### G-sea at ×1.4, p 0.5, attributed (`f142_sea_raw.txt`)

- **The 112 added cells form one line, all ocean-connected (0 inland), covered by no lake.**
  - They lie 30–129 cells (3–13 km) from the nearest kept body, about 10 cells per 10-cell band, in x 2035–2128,
    y 3344–3429.
  - That is downstream of body 3 (lake 4, col (2160, 3454), L_floor 13.4 m), drained at this age (r 1.69, no final
    lake).
  - They stood at 0.8–31.7 m before the breach (p50 18.0).
- **The reading is F141's γ**: the breach's ramp, anchored on a low pit's floor in the emptied bowl, descends 0.113 m
  per cell along the outlet towards the open sea. It reaches the sea before its root.
- Ocean-connected, these cells are not "enclosed below-sea components", so **F38 does not fire**.
- **At ×1.4, p 0 the same lake is fully drained and no cell is added.** Why not was not measured.

## Predictions

**The reviewer's:**
- **R** (the emptying ages within ±0.03 of F139-R for the 8): **refuted**: lake 3 +0.05, lake 19 +0.18, lake 15 goes to
  "never", lake 2 enters.
- **G-bodies 14 / 14**: **held** (1:1 in all eight worlds).
- **G-pits and G-sea hold at every age, F38 fires nowhere (×1.4 p 0 included)**:
  - G-pits **refuted** (from the light pass, in the five worlds with a drained bowl);
  - G-sea **refuted** at ×1.4 p 0.5 (+112, ocean-connected);
  - F38 **held** (8 of 8).
- **G-ring fails on at least one lake**: **held** (95–214 cells).
- **G-deposit < 5 % of the gorge's removal at ×1.4**: **held** (0.03 %).
- **G-tag 5–10 tagged falls at ×1, p 0.5**: **refuted** (11; 9 on the 14 by F142-R).
- **Cost < 5 s per world**: **not resolved** (median +1.25 s; −3.8 to +9.5 s in noise of ~±10 s).
- **Meta**: **held**.

**Mine** (`f142_predictions.md`):
- **P-R** (not blind): lakes 19 (≈ ×1.21), 15 (→ never) and 2 (≈ ×1.38) **held**. "The others within ±0.03"
  **refuted** by lake 3 (+0.05).
- **P-T3** (within 10 % of F138's 1 145.7 at r = 1; L_bed ≈ today's level): **refuted** (1 304.5, +14 %; L_bed is
  2.8–46.6 m above today's level).
- **P-T4** (the count kept, the set moved): **held** (10 of 14 with φ > 0, a different 10).
- **P-G-bodies** (1:1, but the proxy keeps extras): **held** (4 extras).
- **P-G-pits** (holds at the construction, fails later): **held** (from the light pass).
- **P-G-sea and F38 hold everywhere**: G-sea **refuted** (×1.4 p .5); F38 **held**.
- **P-G-ring fails**: **held**.
- **P-G-deposit > 5 %** (against the reviewer): **refuted** (0.03 %).
- **P-G-tag 5–10**: **refuted** (11).
- **P-G-levels ≥ 12 / 14 at ×1, p 0.5**: **refuted** (7 of the 14).
- **P-G-area non-increasing**: **refuted** (+10.6 km² from ×0.7 to ×0.85).
- **P-cost**: not resolved.
- **Meta**: **held**.

## What the table says (no decision taken; nothing promoted)

1. **The construction itself holds what F140 lacked.**
   - No lake is split, F38 never fires, G-sea holds in 7 of 8, and there are no pits at the construction.
   - The rim is uncut at the construction, and G-drained holds to ×1.
   - The plain deposits a negligible volume (≤ 0.13 km³). **It fills what it was asked to fill, and the production
     wiring is verified.**
2. **The light pass undoes it.** One stream-power pass erodes the tagged head falls the construction lays (up to
   418 m on one cell) and cuts canyons back into the lakes and plains (lake 10: −284 m at the lip).
   - That alone accounts for G-pits, and reads as the cause of G-rim on the eroded world and of most of G-levels.
   - **"Chute ou gorge, les deux sont valides"** meets a finish that does not know a fall is meant to stay.
3. **θ falls because the decided steepness (m = 10) puts a ten-times-steeper population at the largest areas.**
   Without those 3.5 % of links θ is 0.478.
4. **The scope proxy holds the 14, with 4 small extras.**
   - L_bed as built is not today's level, so r ≤ 1 is almost inert: the author's −55 % at ×1 becomes −50 %.
5. **Lake 13 joins 1000016's merged water once lowered below ~116 m.** A D8 lake by the input, it becomes part of
   a below-sea water body with age.
6. **Not located**: G-ring (beside body 0's col), G-drained at ×1.2 / ×1.4, G-tag's drops beyond the lip.

**Options, named for the reviewer, none applied:**
- (i) exempt the tagged falls from the light pass, or tag them as non-erodible;
- (ii) lay no head fall (φ = 0) under the light pass;
- (iii) measure θ without the gorge paths, or retune m;
- (iv) locate G-ring / G-drained / G-tag before any v5.

## Limitations, stated

1. **One world, one seed** (the témoin).
2. **The light pass as G-pits' cause** is read from 9 bodies in rank (head fall against the ring's drop), not tested
   by a run without falls.
3. **The cost** is not resolved below the timing noise.
4. **The coast at ×1.4** is compared with ON ×1 only.
5. **G-pits measures only bodies with r_lake > 1** (the plain's bodies). At r ≤ 1 the light pass's cuts show in G-rim
   and G-levels instead.
6. **The first G run was stopped and its G-pits measure amended after seeing its first world** (declared,
   non-blind; the output is kept).

## State

**Uncommitted**, awaiting the go-ahead.
- `valley_construction.rs`: `GorgeRetreat::{d8_scope, plain, v4}`, the scope in `gorge_bodies`, the `GorgeBody`
  diagnostics, `gorge_plain` / `GorgePlain`, two tests.
- `production_upscale.rs`: the plain's call.
- `workspace.rs`: v4 (the toggle still hidden).
- `f126_coast.rs`: `f142_r`, `f142_g`, `f142_theta`, `f142_sea`.
- ADR: F140 corrected in place (flagged); Finding 142.
- This folder.
