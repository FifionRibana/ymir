# Finding 143 — F142's failures attributed by control: freezing the design through the light pass removes G-pits, G-rim and G-levels at ×1 (the light pass is the cause), φ = 0 does not (the falls are not); the head falls make G-ring and G-slope100; G-drained, G-sea at ×1.4, G-tag and θ are born at the construction; the scope is now exactly 14

**Status: diagnosis plus one fix (the exact scope, gated). Nothing promoted. The viz toggle stays hidden.**

**Raw outputs, in this folder:**
- `f143_l_raw.txt` (S, L, θ, H);
- `f143_c_raw.txt` (the controls);
- `checks_before.txt`, `checks_after.txt`.

The predictions were written before any measurement (`f143_predictions.md`, with the non-blind items declared). The
instruments were declared before both benches (`f143_declared.md`).

## Partie 0

1. **F142 committed** (`7a91ea0`): the gated code with the toggle hidden. Added to ADR F142 ("Recorded at commit
   (F143)"):
   - the light pass's cause is correlated over 9 bodies, not attributed;
   - the light pass is a second age process (grid celerity, F115 Courant 3 699): the age is counted twice;
   - head falls up to 418 m on large outlets, F139's objection; the anchors (Victoria, Iguazú, Niagara under
     ~100 m) are to be sourced;
   - θ: the gorges are designed breaks, outside the graded trunk's θ; the rebased links' effect unexplained;
   - R: r = 1 is no longer ON's level; −55 % → −50 %.
2. **Checks**:
   - before = F142's after-checks (`7a91ea0`; no file under `src` changed between them and this round's start);
   - after = `checks_after.txt`.

## S — the exact scope (the round's only fix; `f143_l_raw.txt`)

**Built (gated, a per-world list):** `GorgeRetreat::scope_lows: Option<[u32; 16]>`, the lowest cells of the bodies
to keep (`u32::MAX` pads).
- It replaces the proxy test in `gorge_bodies`.
- It is skipped in the key when `None`, so v4's key and caches are unchanged.
- **The construction knows no final lake** (F132-P4's circularity). So the exact set is a measurement: ON ×1's 14 D8
  lakes matched by footprint to the v4-proxy bodies. In the témoin: 14 lowest cells (listed in the raw output).
- **Another world needs its own list.** Without one, the proxy stays.

**The 4 bodies that were in excess, and why the proxy kept them.** All four are input lakes whose ring touches no
other body's depression cell and no sea cell: the proxy's definition. None is a D8 lake of the final world.

| body | at | km² | L_in | A (km²) | in ON: ring min on the construction / at (d) | ON's final lakes over it |
|---|---|---|---|---|---|---|
| 5 | (3549, 3970) | 5.2 | 1 148.2 | 76.3 | 277.0 (−871) / 274.7 (−874) | none |
| 9 | (1668, 4125) | 5.3 | 723.4 | 287.2 | 400.8 (−323) / 458.4 (−265) | none |
| 11 | (1897, 4546) | 8.2 | 1 257.4 | 78.6 | 131.7 (−1 126) / 320.8 (−937) | none |
| 15 | (2541, 5265) | 6.1 | 163.7 | 1 108.3 | 117.8 (−46) / 119.9 (−44) | 1000016 (507 cells) |

- **Bodies 5, 9 and 11 are drained in ON by its own construction.** A valley cuts their rim by 265–1 126 m, so no lake
  is left.
- **Body 15 lies under the merged below-sea lake 1000016**, which holds it at 115.3 m.

**Permanent test** (rule 13; lib 595 → 596): `the_exact_scope_keeps_only_the_listed_bodies`.
- Negative control first: the proxy keeps both lone bowls. The list keeps only the listed one, and the other's cells
  keep the extended base bit for bit.
- The list, φ = 0 and the frozen design each move the key. v4 serialises as before.
- **The bench options** `phi_zero` (C1) and `freeze_design` (C2: `gorge_design_mask` + the restore after the light pass
  in `production_upscale.rs`) are built alongside, gated, off and absent from the key.

## C — the decisive controls (`f143_c_raw.txt`; the exact scope; G-bodies **exactly 14, 1:1** in all twelve worlds)

**No world panicked and F38 fired in none (12 of 12).** **G-inert holds**: with the gate off the after-checks read the guard 6 / 6, field and lakes (C2 /10 col = `a8d2d538d692c2f0`), lib 596, viz 30, `cargo check` clean. "con." = on the construction (C2 = C0 and C3 = C1 there, by
construction); "fin." = on the final world (G-ring, G-drained and G-tag on the eroded world).

### ×1, p 0.5

| gate | C0 con. | C0 fin. | C1 (φ = 0) con. | C1 fin. | C2 (frozen) fin. | C3 fin. |
|---|---|---|---|---|---|---|
| G-pits | 0 | **43 070** | 0 | **43 078** | **0** | **0** |
| G-levels (±1 m) | 12/14 | 7/14 | 12/14 | 8/14 | **12/14** | **12/14** |
| G-rim (> 10 m) | 0/6 | 3/6 | 0/6 | 3/6 | **0/6** | **0/6** |
| G-ring | **68** | 70 | **0** | 12 | **305** | 208 |
| G-slope100 > 28° (max) | 2/13 (68.4°) | 2/13 (35.6°) | **0/13** (26.5°) | 2/13 (32.8°) | 2/13 (**57.1°**) | 3/13 (**52.9°**) |
| G-drained | 0 | 0 | 0 | 0 | 0 | 0 |
| G-tag drops | 60 | 249 | 159 | 342 | 58 | 159 |
| tagged falls | 9 | 9 | 0 | 0 | 9 | 0 |
| G-area, the 14 (km²) | 702.3 | 694.4 | 702.3 | 695.8 | 703.5 | 703.5 |
| G-sea | +0/−0 | +0/−0 | +0/−0 | +0/−0 | +0/−0 | +0/−0 |
| θ | 0.284 | 0.238 | 0.259 | 0.222 | 0.242 | 0.223 |
| removal (kept bodies), km³ | — | 4 654.6 (18.2) | — | 4 643.9 (18.1) | 4 649.1 (13.9) | 4 637.9 (13.9) |

### ×1.4, p 0.5

| gate | C0 con. | C0 fin. | C1 con. | C1 fin. | C2 fin. | C3 fin. |
|---|---|---|---|---|---|---|
| G-pits | 0 | **9 360** | 0 | **8 245** | **0** | **0** |
| G-levels / drained absent | **1/6** · 7/8 | 1/6 · 7/8 | 1/6 · 7/8 | 1/6 · 7/8 | 1/6 · 7/8 | 1/6 · 7/8 |
| G-ring | 37 | 47 | 0 | 7 | 520 | 429 |
| G-slope100 > 28° (max) | 1/4 (40.5°) | 1/4 (32.6°) | 0/4 (9.8°) | 0/4 (15.7°) | 1/4 (37.0°) | 0/4 (14.5°) |
| G-drained | **2 205** | **10 653** | **2 177** | 9 647 | 7 567 | 7 178 |
| G-tag drops | 23 | 108 | 55 | 133 | 32 | 55 |
| G-sea | **+111 (vs ON's con.)** | **+112** | +111 | +112 | **+122** | +122 |
| θ | 0.397 | 0.339 | 0.369 | 0.323 | 0.340 | 0.321 |

### ×1.4, p 0 (every lake drained)

| gate | C0 con. | C0 fin. | C1 con. | C1 fin. | C2 fin. | C3 fin. |
|---|---|---|---|---|---|---|
| G-pits | 0 | **8 965** | 0 | **7 844** | **0** | **0** |
| drained absent | 11/14 | 11/14 | 11/14 | 11/14 | 11/14 | 11/14 |
| G-ring | 37 | 42 | 0 | 7 | 516 | 417 |
| G-drained | **2 207** | **13 831** | **2 179** | 12 549 | 10 375 | 9 731 |
| G-tag drops | 20 | 101 | 54 | 127 | 29 | 54 |
| G-sea | +0 | +0 | +0 | +0 | +0 | +0 |
| θ | 0.400 | 0.342 | 0.372 | 0.326 | 0.343 | 0.325 |

**The controls ×1** (θ on the construction / on the eroded world):
- ON: 34 lakes, θ 0.493 / **0.423**;
- OFF: 26 lakes, θ 0.493 / **0.436**.

**Canyons 0 everywhere. The coast is unchanged within an age**: 0.0142 /km at ×1 for every control (= ON / OFF), and
0.0050 at ×1.4.

### The answers

**What disappears under C1 (φ = 0)**:
- **G-ring at the construction** (68 → 0; 37 → 0 at ×1.4);
- **G-slope100** (max 68.4° → 26.5° at ×1);
- the tagged falls (9 → 0).
- **G-ring is the head fall's lip** next to the clamped ring.
- **NOT G-pits** (43 070 → 43 078 at ×1; −12 % at ×1.4), **NOT G-levels** (7 → 8 of 14), **NOT G-rim** (3 of 6).
  **The head falls are not what the light pass erodes into the lakes**: F142's correlation is refuted by this
  control.
- G-tag's drops rise (60 → 159): without the fall, the gorge descends over more cells at S > H_f / 2 cells, the
  instrument defect.

**What disappears under C2 (the design frozen through the light pass)**:
- **G-pits** (0 in all three worlds);
- **G-rim** (0/6);
- **G-levels at ×1** (12 of 14, the construction's value).
- **The light pass is the cause, attributed by the control.** It erodes the designed geometry: the lips and the
  gorges, not the falls.
- **C2 makes its own failures**: walls at the frozen mask's edge (G-ring on the eroded world 305–520; G-slope100 up
  to 57°).
- C2 does not remove what is born at the construction (below).

**Under C3**: as C2, with C1's lip gone. G-ring at the mask's edge is still 208–429; G-slope100 reaches 52.9° at ×1.

**What already exists BEFORE the light pass** (on the construction):
- **G-ring** (68 / 37 / 37: the head fall's lip);
- **G-tag's drops** (60 / 23 / 20: mostly the instrument, a few hanging junctions);
- **G-drained at ×1.4** (2 205 / 2 207: other lines' cones). The light pass then multiplies it (×5–6).
- **G-sea at ×1.4, p 0.5** (+111: the breach on the construction already opens the inlet below lake 4);
- **G-levels at ×1.4** (1 of 6: partly drained lakes do not stand at L(r) even on the construction);
- **θ** (0.284 at ×1 against ON's 0.493: the gorge corridors).

**Passing at the construction**:
- G-pits (0 in every world), G-rim (0), G-levels at ×1 (12 of 14), G-sea at ×1 and at ×1.4 p 0;
- G-drained at ×1 (0).

**θ on the eroded world**: ON 0.423 and OFF 0.436, so **the light pass lowers θ by ~0.06–0.07 even without the
gorge**. The non-regression reference for an eroded-world θ is OFF's eroded 0.436 [0.428, 0.444], not 0.493. The
gorge stays below it (0.222–0.343).

## L — what was not located (`f143_l_raw.txt`)

### G-ring: born at the construction (the clamp beside cut cells)

At ×1, p 0.5 with the exact scope:
- **68 clamped ring cells with a > 28° drop to an outside neighbour on the construction; the same cells give 70 on the
  light pass's world.** Not the light pass.
- Every neighbour was cut by the construction:
  - **33 by the OUTLET line**: the lip. On body 0 the col's flanks fall from L(r) = 182.6 to 145.9 m, its 36.8 m
    head fall, and the clamp holds the ring at 182.6 beside them;
  - **35 by ANOTHER line's** cones.
- **G-ring is the clamp meeting the construction's own cuts**: the head fall's lip, and other valleys' walls.

### G-drained: born at the construction (another line's cone)

| world | under L_floor − 1 m upstream of a drained body | laid by a line draining ELSEWHERE, zf under the floor | laid in the same body (the corridor beside the lip) |
|---|---|---|---|
| ×1.2, p 0.5 | 35 (at `carve` and after the plain) | 22 | 13 |
| ×1.4, p 0 | 2 207 | **2 192** | 14 (+1 upstream) |

- **`carve_diag`'s nearest-sample cone from another catchment's line lays these cells.** For example body 1 (L_floor
  934.5 m): cells at 884–917 m laid from a sample at zf 119–121 m, about 32 cells away.
- This is F141-Z's mechanism (the cone of a neighbouring line with a lower base), here within the D8 scope.
- **The plain does not touch them** (they lie outside the footprints). Present at `carve`.

### G-tag's drops: mostly the instrument; a few real hanging junctions

At ×1, p 0.5 on the construction: **60 drops**, in bodies 1 (18), 7 (18), 6 (11), 2 (8), 9 (4) and 3 (1).
- **Most are the gorge's own designed slope over two cells.** Body 1: S = 3.87° and H_f = 10.1 m. Two cells span
  195–276 m, so the designed slope alone drops 13–19 m > H_f; the body-1 examples are steady 10.5–12.4 m drops.
  - **F140's G-tag compares S × (two cells) with 1.5·S·100 m, so a straight gorge fails it by construction.**
    Declared as an instrument defect, not corrected here.
- **A few are real steps.** Body 1 at 10.8 km from the col: 591.1 → 202.4 m in one cell.
  - The gorge's path cell is laid by ANOTHER line's sample: a hanging junction where the raised gorge meets a lower
    trunk, at the construction.
- **24 of the 60 lie wholly off the carved floor**, 17 of them with no laying sample at all: the path leaves the valley
  onto S1. 18 lie on another line's floor (the confluences). **Only 7 are on the outlet line's floor.**

### G-sea at ×1.4, p 0.5: the breach makes a 4.3 km sea inlet

- **112 cells, all ocean-connected, from 11 pits, all in body 3 (lake 4)'s bowl.** The pit floors are 13.8–25.5 m,
  their spill 61.8 m (L(r)). On the construction those cells stood at 1.1–39.6 m.
- The ramps run 124–226 steps and cross the sea at steps 122–225. They stop at the sea (2034, 3343) or at a cell an
  earlier ramp lowered.
- **They lie up to 4.30 km from the coast of (d) (p50 2.64 km): the sea advances inland along lake 4's outlet.** The
  stage is the breach (F141's γ), once the light pass has drained the lake (lip −76 m at ×1, no final lake at
  ×1.4).

## θ — the rebased links (`f143_l_raw.txt`)

S / S_law on OFF's carved trunk links, by distance upstream of the nearest kept body:

| links | stage | 0–1 km | 1–2 | 2–5 | 5–10 | > 10 km |
|---|---|---|---|---|---|---|
| **rebased** | construction | 1.00 (4 %) | 1.00 (4 %) | 1.00 (5 %) | 1.00 (4 %) | 1.00 (3 %) |
| rebased | light pass (b) | 1.00 (10 %) | 1.00 (10 %) | 1.00 (12 %) | 1.00 (15 %) | 1.00 (14 %) |
| rebased | ON construction | 1.00 (4 %) | 1.00 (5 %) | 1.00 (7 %) | 1.00 (4 %) | 1.00 (3 %) |
| **unchanged** | construction | **10.00** (61 %) | 1.23 (47 %) | 1.00 (27 %) | 1.00 (8 %) | 1.00 (3 %) |
| unchanged | light pass (b) | 4.61 (63 %) | 1.13 (51 %) | 1.03 (35 %) | 1.01 (27 %) | 1.00 (27 %) |
| unchanged | ON construction | 1.00 (7 %) | 1.00 (4 %) | 1.00 (4 %) | 1.00 (4 %) | 1.00 (3 %) |

(median; in brackets the share outside [0.5, 2].)

- **At the construction the rebased links lie on their law exactly as ON's do** (median 1.00, 3–5 % outside).
  **The light pass spreads them** (10–15 % outside).
- **The departure from the law in the construction is on the UNCHANGED links within 2 km of the bodies**: the gorge
  corridors, at ten times the law (m = 10). They are "unchanged" because the gorge raises their base, it does not
  lower it.
- **What stays unexplained**: F142-θ read θ = 0.435 for the rebased links alone on the construction, while their ratio
  is 1.00 with 4 % tails. The median does not move the fit; the tails (or the area the θ fit uses) would have to.
  Not decomposed.

## H — the head falls at ×1, p 0.5 (the exact scope; a table for the author)

| lake (body) | A km² | D_g m | φ | head fall m | lip drop in C0 m | H_cap(A) m |
|---|---|---|---|---|---|---|
| 10 (6) | 819.4 | 695.8 | 0.50 | **348.6** | **−283.5** | 71 |
| 3 (2) | 598.7 | 917.6 | 0.32 | **299.7** | **−169.4** | 84 |
| 7 (4) | 244.7 | 403.7 | 0.44 | **178.8** | −86.5 | 131 |
| 19 (12) | 727.2 | 123.9 | 0.47 | 60.4 | −18.8 | 76 |
| 1 (0) | 1 371.7 | 182.0 | 0.19 | 36.8 | −7.6 | 55 |
| 4 (3) | 300.7 | 150.3 | 0.21 | 34.2 | **−76.2** | 118 |
| 8 (5) | 156.4 | 86.3 | 0.30 | 29.9 | −4.3 | 164 |
| 12 (9) | 425.1 | 108.9 | 0.11 | 14.0 | −1.9 | 99 |
| 15 (11) | 416.6 | 52.0 | 0.20 | 13.0 | −5.6 | 100 |
| 2 (1) | 442.7 | 1 144.4 | 0 | 4.7 | −5.5 | 97 |
| 20 (13) | 227.2 | 53.9 | 0 | 3.3 | **−45.1** | 136 |
| 9 (7) | 393.3 | 322.3 | 0 | 1.8 | −6.0 | 103 |
| 11 (8) | 947.1 | 140.2 | 0 | 1.6 | 0.0 | 66 |
| 13 (10) | 983.3 | −46.7 | 0.24 | −9.4 (no fall) | 0.0 | 65 |

| cap | tagged | median (m) | max (m) |
|---|---|---|---|
| 50 m | 9 | 36.8 | 50.0 |
| 100 m | 9 | 36.8 | 100.0 |
| 200 m | 9 | 36.8 | 200.0 |
| none | 9 | 36.8 | 348.6 |
| H_cap(A) = 100 m·(A/418.7)^−0.5 (PROXY) | 9 | 36.8 | 130.8 |

- **No cap changes the tagged count** (9): every cap ≥ H_f keeps a tagged fall tagged. A cap acts on the heights only.
- **Three lakes carry a fall > 100 m today** (10, 3, 7). H_cap(A) cuts lakes 10 and 3 to 71 and 84 m.
- **The lip drop against the head fall.** The three large falls lose 48–81 % of their height at the lip.
  - **Two lakes with small falls lose a lot too**: lake 4 (34 m → −76 m) and lake 20 (3.3 m → −45 m). Both are low
    lakes (L_floor 13.4 and 21.5 m).
  - **So the head fall is not the light pass's only lever.** C1 says how much it is.
- (The exact scope changed some heads from F142-R's: lake 10 418 → 349, lake 3 406 → 300. Removing the 4 extra
  bodies changes the outlet paths' ends, so their D_g.)

## Predictions

**The reviewer's:**
- **C1**:
  - "at ×1 G-pits, G-levels, G-rim and G-area pass": **refuted**: G-pits 43 078, G-levels 8/14 and G-rim 3/6. G-area
    is close to C0's.
  - "at ×1.4 G-pits still fails": **held** (8 245; 7 844).
- **C2** ("every failed gate of F142 passes except G-ring"): **refuted**.
  - G-pits, G-rim and G-levels at ×1 pass.
  - **G-drained (7 567–10 375), G-sea at ×1.4 p .5 (+122), G-levels at ×1.4 (1/6), θ and G-slope100 (57°) still
    fail.** Several of them are born at the construction.
- **G-ring is already at the construction (the clamp)**: **held** (68). It is the clamp beside the head fall's lip;
  φ = 0 removes it.
- **G-drained is born in the light pass**: **refuted**. 2 205 / 2 207 cells on the construction, laid by other lines'
  cones; the light pass multiplies them.
- **θ: the rebased links' departure is born in the light pass**: **held** for the slope against the law (construction
  1.00 like ON; the light pass spreads them, 10–15 % outside [0.5, 2]). F142's θ = 0.435 for these links stays
  unexplained.
- **H**: "under a 100 m cap the tagged count stays within ±2 of 11": **held** (9, at every cap). "At least 3 lakes over
  100 m": **held** (3).
- **Meta**: **held**.

**Mine** (`f143_predictions.md`):
- **P-S**: **held** (4 dropped; the scope is a per-world list).
- **P-C1**:
  - "G-levels still fails": **held** (8/14);
  - "at ×1.4 G-pits remains": **held**;
  - "G-pits at ×1 falls by > 90 %": **refuted** (unchanged);
  - "the lip drops under 10 m": not measured in C1.
- **P-C2**:
  - "G-pits = 0 in ≥ 2 of 3": **held** (3 of 3);
  - "G-rim on the eroded world passes": **held**;
  - "G-levels improves, not to all": **held**;
  - "a new failure at the mask's edge": **held** (G-ring 305–520, G-slope100 57°);
  - "G-ring, G-drained and G-tag exactly as C0": **refuted** (G-ring rises, G-drained and G-tag fall).
- **P-C3** (0 tagged): **held**.
- **Before the light pass**:
  - G-ring, G-drained, G-tag fail: **held**;
  - G-pits and G-rim pass: **held**;
  - G-levels ≥ 12/14 at ×1: **held**;
  - **G-sea passes: refuted** (+111 at ×1.4 p .5 on the construction).
- **P-L-ring** (the outlet line, not another): **refuted** (33 / 35).
- **P-L-drained** (another line's cone): **held**.
- **P-L-tag** (where the path leaves the carved floor): **mostly refuted**. The main cause is the instrument's
  threshold against the designed slope; 24 of 60 lie off the floor.
- **P-L-sea** (a ramp from lake 4's bowl to the sea, a sea inlet): **held** (4.30 km).
- **P-θ** (the rebased links depart at the construction): **refuted** (1.00, as ON).
- **P-H**: the counts **held**. The lip-drop rule is **refuted**: lake 7 loses 48 % (< half), and lakes 4 and 20
  lose 76 and 45 m with falls of 34 and 3 m.
- **Meta**: **held**.

## What the controls say (no decision taken; only the scope is fixed)

1. **The light pass, not the falls, reopens the lakes and the plains.**
   - φ = 0 changes nothing to G-pits, G-levels or G-rim at ×1.
   - Freezing the designed geometry through the light pass makes all three pass.
   - The light pass erodes what the construction designs: the lips and the gorges at m = 10, ten times the graded
     slope.
   - **This is the "second age process"** recorded at F142's commit: r already sets the retreat, the light pass adds
     one at a grid celerity.
2. **The head falls make G-ring and G-slope100**: the lip beside the clamped ring. φ = 0 removes both.
3. **Born at the construction, independent of the light pass and the falls**:
   - G-drained at ×1.4 (other lines' cones: F141-Z's mechanism);
   - G-sea at ×1.4 p .5 (the breach's ramp from lake 4's emptied bowl);
   - G-levels at ×1.4 (the partly drained lakes do not stand at L(r));
   - G-tag (the instrument's threshold, and a few hanging junctions);
   - θ (the steepened gorge corridors).
4. **Freezing the design is not a fix**: it makes walls at the mask's edge (G-ring 305–520, G-slope100 57°).
5. **H**: a cap changes the heights, not the tagged count. Three falls exceed 100 m. H_cap(A) would bring them to
   71–131 m.

**Options, named for the author (none applied):**
- (i) no light pass over the designed geometry, with a transition instead of a hard mask;
- (ii) a gorge steepness that the light pass leaves graded (m smaller), against the decision m = 10;
- (iii) a cap on the head falls (the H table);
- (iv) for G-drained: the cones of a line limited to its own catchment;
- (v) for G-sea at ×1.4: F141's queued breach defect;
- (vi) correct G-tag's threshold (H_f over two cells) and measure θ on the eroded world against OFF's eroded θ.

## Limitations, stated

1. **One world, one seed.**
2. **C2's mask** (footprints, rings, gorge corridor of radius 2) is one declared choice. Another radius would move its
   edge walls.
3. **C2 / C3's construction gates are C0 / C1's** (identical constructions by design, not re-run).
4. **The lip drop under C1** was not measured directly (G-rim and G-pits stand for it).
5. **F142's θ = 0.435 for the rebased links** is not explained by their slope ratios.

## State

**Uncommitted**, awaiting the go-ahead.
- `valley_construction.rs`:
  - `GorgeRetreat::{scope_lows, phi_zero, freeze_design}` (the last two are bench options);
  - the exact scope in `gorge_bodies`;
  - `Skeleton::gorge_path`, `gorge_design_mask`;
  - the test `the_exact_scope_keeps_only_the_listed_bodies`.
- `production_upscale.rs`: the frozen-design restore around the light pass (bench option).
- `f126_coast.rs`: `f143_l`, `f143_c`.
- ADR: F142's additions (committed with F142), Finding 143.
- This folder.

**No other fix.**
