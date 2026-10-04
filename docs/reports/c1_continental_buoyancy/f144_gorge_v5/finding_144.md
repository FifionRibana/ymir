# Finding 144 — spec v5 (the input scope, the outlet-base bound, the capped falls) and phase 1 of the light pass's candidates: the input scope keeps 18; the bound stops four lakes at a merged below-sea level but not lake 4; no candidate (floor, transition, order) holds G-pits in the three worlds, all make walls at the design's edge, the order (P3) costs a second construction; REPORT AND STOP, no phase 2

**Status: built gated (v5, off by default), measured, NOT promoted. The viz toggle stays hidden. The declared rule
stopped the round after phase 1.**

**Raw outputs, in this folder:** `f144_p_raw.txt`, `checks_before.txt`, `checks_after.txt`.

The predictions were written before any measurement (`f144_predictions.md`, with the non-blind items declared). The
instruments and the selection rule were declared before the bench (`f144_declared.md`). The specification is
**`spec_gorge_age_v5.md`**.

## Partie 0

1. **F143 committed** (`2adba34`). Added to ADR F143 ("Recorded at commit (F144)"):
   - C1 refutes the "falls" correlation (the reviewer had taken it up); C2 establishes the light pass on the lips and
     gorges as the cause;
   - a hard mask is not a fix (C2's walls);
   - the pattern "one invariant, two paths";
   - the list of 14 depended on ON's defect (bodies 5, 9, 11);
   - L(r) can go below the outlet's base;
   - the instruments: H_f probably 10·S_loi; θ against OFF's eroded θ.
2. **The author's decisions (2026-10-02), recorded as given**:
   - « le recul de la gorge s'applique à tous les lacs d'entrée hors des bassins sous la mer (17) »;
   - « la hauteur des chutes de sortie : plafond selon la taille de la rivière (~70–130 m) ».
3. **`spec_gorge_age_v5.md`** = v4 + these decisions + S, B, H and P.
4. **The toggle stays hidden.** Checks:
   - before = F143's after-checks (`2adba34`; no file under `src` changed between them and this round's start);
   - after = `checks_after.txt`.

## S — the scope computed on the input (`f144_p_raw.txt`)

**Built** (gated, `GorgeRetreat::input_scope`). A body is kept iff it is an input lake, none of its cells lies in a
`basin_base` closed depression, and none lies under the sea.

- **It keeps 18 bodies, not 17**: F143's 14, plus bodies 5 (3549, 3970), 9 (1668, 4125), 11 (1897, 4546) **and 15
  (2541, 5265)**.
- **The 7 dropped** are the input lakes inside a merged below-sea water body: bodies at (4344, 3579) (F141's body 4),
  (4735, 3564), (3721, 4958), (1796, 5239), (2154, 5284), (1523, 5500) and (3112, 6002). For each, in a depression =
  true.
- **Body 15 is not in a depression on the input.** In ON it is 20 % under 1000016, because ON's construction lowers it
  (its ring −46 m, F143). **The input criterion cannot drop it without the final world.**
  - So the author's "17" is 18 by the criterion he described. **Named for the author.**
- **Permanent test** (rule 13; lib 596 → 597): `v5_keeps_the_input_lakes_outside_the_basins`.
  - Negative control first: v3 treats the basin's land fringe, a body in a depression. v5 drops it and keeps the
    lone bowl; every kept body passes the criterion.
  - The scope moves the key and is absent when off.
  - **An input lake NESTED in a depression could not be produced by the synthetic pre-drainage** (three geometries
    tried: the nested pit was absorbed into the depression component). The témoin's 7 dropped bodies are that case.

## B — L(r) bounded by the outlet's base (`f144_p_raw.txt`)

| age, p | bodies raised (level ← B; their L_floor) |
|---|---|
| ×0.7, ×0.85 (p .5); ×1, p 0 | none |
| ×1, p .5 | lake 13 (2439, 5070), r 1.53: 123.1 (B 123.1; L_floor 13.9) |
| ×1.2, p .5 | lake 11 (1463, 4246) 49.3 (= 1000014's level); lake 13 115.3 and body 15 115.3 (= 1000016's level) |
| ×1.4, p .5 | the same three + lake 15 (2652, 5041) 115.3 (r 1.99, L_floor 84.2) |

- **The bound binds the lakes whose outlet enters a merged below-sea water body**: 1000014 at 49.3 m, 1000016 at
  115.3 m. They stop at that water's level instead of draining to their floor.
- At ×1, lake 13 is first held at 123.1 m, the next base on its path.
- **Lake 4**: L_floor 13.4 m and **B = 0.5 m (the sea's `base_m`)** at every age. **The bound does not bind it** (its
  floor lies above the sea). Its level is L(r): 172.1 → 129.0 → 61.8 m from ×0.7 to ×1.4.
- **The 4.3 km sea inlet stays**: G-sea +112 at ×1.4, p .5 under P0, as in F142 / F143. **It is not a level error.**
  It comes from the breach's ramps from the residual pits of lake 4's bowl (F141's γ), which the pre-breach drainage
  does not hold as a lake.
- **No lake is bounded by the sea with its floor below it** in the témoin: no case of "a lake at the sea's level".

## H — the capped head falls (×1, p 0.5)

- **11 tagged on the 18; median 60.4 m; max 239.9 m.**
- H_cap(A) = 100 m·(A/418.7)^−0.5 rises to ~240 m on the smallest outlets (body 5: A 76 km²). **The ~70–130 m the
  author named holds for the 14's outlets (A 150–1 400 km²), not for the three small extras.**
- **H_f already used the real slope**: the code is `(1.5 * sg * 100.0).max(10.0)` with `sg = slope_for(top)`. **The
  F143-commit item "H_f probably 10·S_loi" was wrong**; read here, as declared non-blind in the predictions.
- **G-tag's defect is its span**: the old instrument counts 203 drops on the construction. Corrected (the excess over
  the designed slope over the same two cells), 82 remain: the hanging junctions and the path cells laid by other
  lines (F143-L).

## Z — the laying cones' catchment (measured, nothing changed)

| ×1.4 | G-drained's violators | the drained bowls' footprints | upstream of the drained bowls | every carved cell |
|---|---|---|---|---|
| p 0.5 | **2 174 / 2 190 (99.3 %)** | 40 128 / 221 699 (18.1 %) | 226 858 / 810 276 (28.0 %) | 917 177 / 2 516 993 (36.4 %) |
| p 0 | **2 179 / 2 318 (94.0 %)** | 44 609 / 263 204 (16.9 %) | 251 402 / 958 767 (26.2 %) | 921 422 / 2 541 143 (36.3 %) |

(A cell is "outside" when its laying cone's line is neither the first line its D8 path meets nor a line that one joins
downstream. "Without a laying sample" cells are excluded from the shares.)

- **G-drained is almost entirely laid from outside the catchment.**
- **The same happens to 36 % of the world's carved cells.** Part of that is legitimate: a wall cell near a divide is
  laid by the nearest valley, even if its D8 path leaves on the other side. **The instrument does not separate the
  two**; the share is a ceiling, not a defect count.

## The cost (the red flag on world generation; recorded at the author's request)

The author's reference, a viz `run_hd` of the témoin (2026-10-04): **run_hd TOTAL 249.8 s** (eroded 149.2 s, breach
46.3 s, drainage 30.8 s, prebreach drainage 10.9 s). It already sits in the high part of acceptable durations.

| candidate | ×1 p .5 | ×1.4 p .5 | ×1.4 p 0 | added against P0 |
|---|---|---|---|---|
| P0 (v5, light pass free; v5 itself ≈ +1 s against ON, F142) | 140 s | 215 s | 139 s | — |
| P1 (floor) | 156 s | 207 s | 140 s | +16 / −8 / +1 s (within noise) |
| P2 (transition; = P1 here) | 153 s | 188 s | 147 s | +13 / −27 / +8 s (within noise) |
| P3 (order: a second, bare construction) | **329 s** | **269 s** | **259 s** | **+189 / +54 / +120 s: RED FLAG** (+22–76 % of run_hd) |

`build_world` stands for the viz's "eroded" stage (149.2 s in the author's log). These are single bench timings,
±30 s of noise (×1.4 p .5's P0 at 215 s is itself high). The machine slept during P3 ×1's later stages, not during
its build.

**P1 / P2 add nothing measurable. P3 adds a full extra skeleton + carve.**

## P — phase 1: the light pass and the designed geometry (`f144_p_raw.txt`)

**d_t, measured before P2 by the declared rule: 1 cell = 49 m.**
- C2's steep share minus P0's is −2.46 pp in band 1 and exactly +0.00 pp in bands 2–20.
- **The instrument was blind to C2's walls.** It counted the steep drops *from* cells outside the mask, while C2's
  walls are drops from the frozen mask cells *into* the eroded cells around them.
- Applied as declared, d_t = 49 m makes **P2 identical to P1**: no lowering on the mask, full erosion from the first
  cell outside. The two give the same numbers in all three worlds. **P2 as intended (a real transition) was not
  tested.**

| world | candidate | held / 14 | G-pits | G-levels | G-rim | G-ring | G-slope100 (max) | G-drained | G-tag (corr.) | G-sea | θ (corridors excl.) | edge | build (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ×1, p .5 | P0 | 5 | 43 190 | 7/18 | 5 | 99 | 1 (31.5°) | 0 | 194 | +0 | 0.444 OUT | 573 | 140 |
| | P1 | 6 | 1 792 | 15/18 | 0 | 308 | 4 (52.4°) | 0 | 94 | +0 | 0.445 OUT | 4 936 | 156 |
| | P2 (= P1) | 6 | 1 792 | 15/18 | 0 | 308 | 4 (52.4°) | 0 | 94 | +0 | 0.445 OUT | 4 936 | 153 |
| | P3 | 6 | **187** | 15/18 | 0 | 561 | 7 (**80.1°**) | 0 | 84 | +0 | 0.442 IN | **8 369** | **329** |
| ×1.4, p .5 | P0 | 3 | 6 461 | 2/9 | 1 | 103 | 1 (32.6°) | 9 781 | 110 | **+112** | 0.440 IN | 500 | 215 |
| | P1 / P2 | 4 | 2 330 | 3/9 | 0 | 550 | 2 (48.1°) | 7 241 | 31 | **+118** | 0.441 IN | 6 242 | 207 / 188 |
| | P3 | 2 | 2 037 | 3/9 | 0 | 790 | 3 (76.2°) | 9 186 | 22 | **+118** | 0.444 OUT | 8 152 | 269 |
| ×1.4, p 0 | P0 | 7 | 6 336 | — | 0 | 105 | 0 | 14 366 | 107 | +0 | 0.442 IN | 523 | 139 |
| | P1 / P2 | 7 | 2 487 | — | 0 | 586 | 0 | 11 272 | 36 | +0 | 0.443 IN | 6 367 | 140 / 147 |
| | P3 | 6 | 2 234 | — | 0 | 850 | 0 | 15 920 | 22 | +0 | 0.446 OUT | 8 230 | 259 |

**G-inert holds**: with the gate off the after-checks read the guard 6 / 6, field and lakes (C2 /10 col = `a8d2d538d692c2f0`), lib 597, viz 30, `cargo check` clean.

**Common to all twelve worlds:**
- F38 never fires; canyons 0, except **P3 at ×1 and ×1.4 p .5 (canyons fail)**;
- the kept bodies are 18.
- **G-bodies** against the gorge world's own lakes is 11–15 of 18, never 18. The drained or merged bodies have a lake,
  and some undrained bodies none.
- **drained absent**: 8/9 at ×1.4 p .5, 15/18 at ×1.4 p 0.
- **The edge cells on the construction** are 209–221 (born there). The light pass multiplies them under P1–P3.

**Phase 1's totals (14 gates × 3 worlds = 42):**

| candidate | held | G-ring (sum) | G-pits in every world |
|---|---|---|---|
| P0 (reference) | 15 | 307 | no |
| P1 | **17** | 1 444 | no |
| P2 (= P1 at d_t = 49 m) | **17** | 1 444 | no |
| P3 | 14 | 2 201 | no |

**The declared rule: no candidate holds G-pits in all three worlds, so REPORT AND STOP. No phase 2 was run.**

**Read across the candidates:**
- **Every candidate cuts G-pits** by 64–99 % against P0, and none reaches 0. **P1** brings G-levels at ×1 from 7 to
  15 of 18 and G-rim to 0. **P3** comes closest on G-pits at ×1 (187).
- **Every candidate makes walls**: edge cells ×9–16 and G-ring ×3–8 against P0, slopes to 48–80°.
  - P1 / P2's floor is C2's hard edge with deposition allowed.
  - P3 is the worst: the design laid back above the eroded surroundings, up to 80°.
- **P3 costs +54 to +189 s per world** (a second, bare construction), against a run_hd of ~250 s: **a red flag**
  (the cost section).
- **G-sea at ×1.4, p .5 stays** under every candidate (+112 to +118): F141's γ, untouched by P or B.
- **G-drained at ×1.4 stays** (7 241–15 920): the cones from outside the catchment (Z), untouched by P.
- **θ (corridors excluded) sits on OFF's eroded CI's upper edge**: 0.440–0.446 against [0.428, 0.444]. IN or OUT by
  a few thousandths, so the gate decides on noise.

## Predictions

**The reviewer's:**
- **S** ("exactly 17: the 14 + 5, 9, 11"): **refuted**, 18 (body 15 kept). "5, 9 and 11 are present lakes in the
  gorge world": not tabulated per body. G-bodies is 11–15 / 18 overall.
- **B** ("lake 4's floor is under its outlet base, and the inlet disappears"): **refuted**. Lake 4's floor (13.4 m)
  lies above its base (the sea, 0.5 m), the bound does not bind it, and the inlet stays (+112).
- **H**:
  - "at ×1, 10–13 tagged on the 17": **held** (11 on 18);
  - "none above 131 m": **refuted** (239.9 m; H_cap reaches ~240 m on the small outlets);
  - "with H_f corrected, no gorge cell passes": **refuted** (corrected G-tag 82 on the construction). H_f was already
    the real slope.
- **P**:
  - "P3 holds G-pits, G-levels and G-rim on the three worlds and is selected": **refuted** (G-pits 187 / 2 037 /
    2 234; none selected);
  - "P1 reproduces C2-type walls (G-ring > 100)": **held** (308–586);
  - "P2 halves G-ring against C2 but fails G-pits at ×1.4": **G-pits held, G-ring not testable** (P2 = P1 at d_t =
    49 m).
- **θ in OFF's eroded CI at ×1, corridors excluded**: **refuted at the edge**. P0 0.444 and P1 0.445 OUT, P3 0.442 IN,
  on a 0.444 bound.
- **Z** (> half of G-drained's violators laid from outside their catchment): **held** (99.3 % / 94.0 %).
- **Meta**: **held**.

**Mine** (`f144_predictions.md`):
- **P-S** (18, body 15 kept): **held**. "G-bodies fails for ≥ 1 body at ×1.4": **held**.
- **P-B** (lake 4 not bound, the inlet stays; lake 13 bound from ×1.2): **held**. Lake 13 is bound from ×1, and lakes
  11 and 15 and body 15 too.
- **P-H**:
  - 9–12 tagged: **held** (11);
  - "none above its H_cap": **held**;
  - H_f not the cause: **held**;
  - ≥ 1 corrected drop: **held**.
- **P-P**:
  - "P1, P2, P3 all hold G-pits": **refuted** (none);
  - "P1 reproduces C2's walls": **held**;
  - "P2 is selected": **refuted** (no selection; P2 = P1);
  - "P3 makes edge walls": **held**;
  - "no candidate holds G-sea at ×1.4 p .5": **held**.
- **P-θ** (below OFF's eroded CI at ×1): **refuted in part**. P0 and P1 are 0.001 above, P3 is inside.
- **P-Z**: "> 90 % of the violators": **held**. "< 10 % over the bowls and the world": **refuted** (17–36 %).
- **Meta**: **held**.

## What the round says (no decision taken; nothing promoted; no phase 2)

1. **The input scope is 18, not 17.** Body 15 lies in no depression on the input. The author's "17" counts it out only
   through ON's world.
2. **The outlet-base bound works where it should**: four lakes stop at a merged below-sea water level (49.3 / 115.3 m)
   instead of draining to their floor. **It does not touch lake 4**: the inlet is a breach effect, not a level error.
3. **The capped falls** keep 11 tagged. The cap is ~70–130 m on the large outlets and up to ~240 m on the smallest.
   **H_f was never the G-tag defect; the span was.**
4. **No light-pass candidate holds G-pits, and all make walls.**
   - The floor (P1, = P2 here) is the best by count (17 / 42) and P3 the closest on G-pits at ×1.
   - **P3 costs a second construction**: +54 to +189 s per world, a red flag against ~250 s.
   - **The pattern holds**: whatever keeps the design out of the light pass's reach leaves a seam at the design's edge.
5. **Z**: G-drained is laid by cones from outside their catchment (94–99 %). The same holds for 36 % of the world's
   carved cells, an upper bound mixing defects and legitimate walls.
6. **Instruments to fix before any rerun**:
   - d_t (count the drops from the mask outwards);
   - G-bodies' rule for drained / merged bodies;
   - θ's gate on a CI edge.

**Options, named for the author (none applied):**
- (i) the light pass's own invariant: a stream-power floor at the design inside the light pass, not after it;
- (ii) P2 with a d_t measured from the mask outwards (a real transition, untested here);
- (iii) no light pass on the kept bodies' catchments (cost-free; it changes the finish);
- (iv) cones limited to their catchment (Z);
- (v) the breach's queued defect (G-sea);
- (vi) body 15 by decision.

## Limitations, stated

1. **One world, one seed.**
2. **d_t's instrument was blind** to the walls it was meant to measure; P2 was tested only at 49 m.
3. **The cost** is single bench timings (±30 s noise at ×1.4; the machine slept during P3 ×1's extras, not its build).
4. **The construction-stage gates were not re-measured** in phase 1. The edge cells on the construction (209–221) and
   corrected G-tag (82 at ×1, 18 at ×1.4) are given.
5. **Z's catchment test** uses the raw lines' owners and the parent chain. It does not separate defects from legitimate
   divide walls.

## State

**Uncommitted**, awaiting the go-ahead.
- `valley_construction.rs`:
  - `GorgeRetreat::{input_scope, bound_b, head_cap, light_mode, light_dt_m, bare, v5}`;
  - the v5 scope;
  - the outlet-base bound;
  - the capped head fall;
  - `GorgeBody::{in_depression, in_sea, b_out, bounded, head_cap_m}`;
  - `Skeleton::line_owner`;
  - `gorge_design_distance_m`;
  - the test `v5_keeps_the_input_lakes_outside_the_basins`.
- `production_upscale.rs`: the P1–P3 bench options around the light pass.
- `f126_coast.rs`: `f144_p`.
- ADR: F143's additions (committed with F143), Finding 144.
- This folder (`spec_gorge_age_v5.md`).
