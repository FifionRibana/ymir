# Finding 156 — lakes / gorges / falls resumed, round 3 of 4: the spill breach fails its blind criterion on three new seeds (the lakes lose the pits the old ramps drained, not only the ramps' cells; one seed's production world crashes in the drainage), so nothing ships; d_t = 3 cells (146 m, the fallback); the true transition T fails G-pits in all three worlds (exactly F144-P1's counts: the light pass's deposition on the design) and its walls exceed 1.5 × P0, so T is NOT retained; the diagnostic T-all (deposition weighted too) holds G-pits = 0 and G-rim = 0 everywhere

**Status: measured; the transition built gated; nothing promoted; the breach unchanged. Nothing committed.**

**Raw outputs, in this folder:**
- `f156_br_raw_run1.txt`, `f156_br_raw.txt`, `f156_br_loc_raw.txt` (Br);
- `f156_t_raw_run1.txt`, `f156_t_raw_run2_oom.txt`, `f156_t_raw.txt` (T);
- the 12 `f156_ombrage_lac{1,2,11}_{p0,t,ta,np}.png`;
- `checks_before.txt`, `checks_after.txt`.

**Process:**
- The predictions were written before any measurement (`f156_predictions.md`, with the non-blind items declared).
- The instruments, the criterion and the selection rule were declared before the benches (`f156_declared.md`).
- **Four amendments came after a run, declared there as NOT blind**: the Br localisation; the Br bench's robustness
  (`catch_unwind`); T's memory and resume; and the stopping of two IDE processes (below).

## The author's criterion and decisions (2026-10-08), recorded as given

- The breach: « Nouvelle confirmation à l'aveugle sur 2–3 graines, critère précisé ».
- The light pass, after the hillshades with and without it: « Inacceptable : un dernier candidat, la vraie
  transition ». Options (a) and (b) of F155 are set aside.
- Round 4 is the author's look, then the decision: ship the gorge, or pause it for good.

## Partie 0

- **F155 committed** (`42fffe2`). Its ADR notes were committed apart (`01530f1`):
  - the Br stop rule was badly written (not rewritten; a new blind confirmation follows);
  - what the light pass brings, and the author's verdict;
  - the breach's second defect queued (livré's coastal pits).
- **Checks before** = F155's after-checks: no file under `src` changed between them and this round's start.
- **Checks after**: `cargo check` clean; lib **621** (+1, `the_transition_spares_the_design_and_fades_out`); viz 31;
  **the guard 6 / 6, field and lakes = banc** (C2 /10 col `a8d2d538d692c2f0`).

## Br — the blind confirmation on new seeds fails; nothing ships

- **Seeds**: 20261008001 / 002 / 003, never used before (grep-checked).
- **State**: production C2 /10 col, framed by the viz's automatic roll.
- **Copy**: F154's amendment 2, untouched.
- **R**: the cells production's ramps lower, recorded by an instrumented copy checked bit-identical to production.

| seed | lakes | cells gained | cells lost (not an old ramp's) | R under a lake, still under it | new below-sea land | below-sea land: production → spill | raised / lowered | segments |
|---|---|---|---|---|---|---|---|---|
| 20261008001 | 26 → 26 | **0** | 534 (**204**) | 1 097 (767) | **0** | 267 → 17 | 33 587 / 9 | +0.81 % |
| 20261008002 | — | — | — | — | — | — | — | — |
| 20261008003 | 30 → 30 | **0** | 1 828 (**1 380**) | 2 812 (2 364) | **0** | 101 → 13 | 72 005 / 8 | +0.56 % |

**Seed 20261008002**:
- **the run_hd tail panics in BOTH chains, production and spill alike**;
- the error is F85's fixed-point assertion: « the level/outlet/balance fixed point did not converge in 16 passes (20
  class merges) » (`drainage.rs:2831`);
- **this is a production defect, independent of the breach's anchor**: the production world does not generate on
  this seed.

**The stop rule fires on all three seeds** (item 3 on 001 and 003; a chain's panic on 002). **Nothing ships.**

**Localised** (seed 001, amendment), the lost cells outside R per lake:
- 1000001: 80 cells; 1000015: 123; 1000016: 1.
- They lie **1–16 cells from an old ramp** (p50 2–3).
- In the spill world they are **raised by +1 to +73 m** (p50 +17 to +48 m).
- None lies in a pre-breach lake.
- **They are the pits themselves.** The old ramps drained them, leaving them below the lake's level. The spill anchor
  leaves them to the mop-up fill, which raises them to their outlet, out of the lake.
- **The criterion counted only the cells the ramps LOWERED, not the pits they DRAINED.** It is not rewritten here.
  Seed 003 was not localised.

**The reverse gap**: 70 % (001) and 84 % (003) of the old ramps' cells under a production lake stay under a spill
lake. The ramps' trenches inside a lake mostly stay lake.

**The production change prepared but NOT applied**: `RampAnchor`, the skeleton's explicit floor anchor, and the three
ALGO bumps.

**The queue**:
- the first defect (the ramp anchored on the floor) stays;
- so does the second (livré's coastal pits);
- **a third is added: seed 20261008002's F85 non-convergence in production.**

## T — the true transition

**d_t, measured first** (×1 p .5, P0 against C2):
- **The profile** D(j) = median(z_C2 − z_P0) is +0.74 m on the mask's edge cells and **+0.00 m from j = 1 to 20**.
- **Its d_t is 1 cell.** Outside the mask C2 and P0 share the light pass, so their difference there is the droplets'
  and the breach's, not the pass's (F144's blind instrument, again).
- **The declared fallback**: C2's edge step is p50 0.7 m, **p90 69.9 m** and max 421.6 m, over 22 651 edge cells
  (5 337 steeper than 28°). max(2, ceil(69.9 / (tan 28° · 48.8 m))) = **3 cells**.
- **d_t = 3 cells = 146 m.** F144's wall width read 1 cell.

**The rule** (`light_mode` 5, gated):
- the light pass's LOWERING is weighted by smoothstep(d / d_t), with 0 on `gorge_design_mask` (the rings, the
  corridors, the plains, the lake beds);
- deposition is not weighted.
- **T-all** (`light_mode` 6) weights the whole change; it is a declared diagnostic, not eligible.
- `valley_construction::{transition_weight, blend_light_pass}`.
- Permanent test `the_transition_spares_the_design_and_fades_out` (with its negative control).

| world | cand. | held / 14 | **G-pits** | G-levels | **G-rim** | G-ring (in / beyond the fall) | G-slope100 (max; in) | G-drained | G-tag | G-sea | θ | **edge** | build |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ×1 p .5 | P0 | 6 | 43 190 | 7/18 | 5 | 99 (31 / 68) | 1/14 (31.5°; 1) | 0 | 194 | +0 | 0.444 IN | **573** | 191 s |
| | C2 | 8 | 0 | 15/18 | 0 | 310 (47 / 263) | 4/16 (49.1°; 2) | 0 | 84 | +0 | 0.445 IN | 4 858 | 194 s |
| | **T** | 7 | **1 792** | 15/18 | **0** | 168 (41 / 127) | 4/16 (44.2°; 2) | 0 | 94 | +0 | 0.445 IN | **1 195** | 233 s |
| | T-all | 8 | **0** | 15/18 | **0** | 171 (43 / 128) | 3/16 (49.1°; 2) | 0 | 84 | +0 | 0.445 IN | 1 120 | 166 s |
| ×1.4 p .5 | P0 | 3 | 6 461 | 2/9 | 1 | 103 (28 / 75) | 1/5 (32.6°; 0) | 9 781 | 110 | +112 | 0.440 IN | **500** | 192 s |
| | C2 | 5 | 0 | 3/9 | 0 | 552 (45 / 507) | 2/6 (48.1°; 1) | 7 241 | 22 | +122 | 0.441 IN | 6 155 | 154 s |
| | **T** | 4 | **2 330** | 3/9 | **0** | 243 (38 / 205) | 2/6 (36.1°; 1) | 7 176 | 32 | +117 | 0.439 IN | **1 210** | 227 s |
| | T-all | 5 | **0** | 3/9 | **0** | 245 (40 / 205) | 2/6 (43.1°; 1) | 7 176 | 23 | +122 | 0.440 IN | 1 127 | 175 s |
| ×1.4 p 0 | P0 | 7 | 6 336 | — | 0 | 105 (29 / 76) | 0 | 14 366 | 107 | +0 | 0.442 IN | **523** | 160 s |
| | C2 | 8 | 0 | — | 0 | 588 (45 / 543) | 0 | 11 272 | 22 | +0 | 0.443 IN | 6 277 | 155 s |
| | **T** | 7 | **2 487** | — | **0** | 234 (38 / 196) | 0 | 11 118 | 37 | +0 | 0.442 IN | **1 346** | 161 s |
| | T-all | 8 | **0** | — | **0** | 236 (40 / 196) | 0 | 11 118 | 23 | +0 | 0.442 IN | 1 259 | 161 s |

**Common to all worlds:**
- G-drained at ×1 is 0;
- the canyons are 0, **in the zone ≤ d_t too**;
- the coast is within ±2 % of P0;
- G-area is non-increasing ×1 → ×1.4 for every candidate (T: 1 187.8 → 76.8 km²);
- F38 never fires.

**The selection rule (declared)**: « T est retenu s'il tient G-pits = 0 et G-rim = 0 sur les trois mondes, avec des
murs au bord du masque au plus 1,5 fois ceux de P0. »
- **G-rim = 0** in all three: holds.
- **G-pits is 1 792 / 2 330 / 2 487**: fails.
- **The walls are 1 195 / 1 210 / 1 346** against 1.5 × P0 = 860 / 750 / 785: fail (× 2.1–2.6 P0).
- **→ T is NOT retained. Round 4 decides the pause.**

**Read:**
- **T's G-pits are exactly F144's P1** (1 792 / 2 330 / 2 487).
  - On the design, T (w = 0, deposition free) and P1 (a floor at the construction, deposition allowed) are the same
    operation.
  - **The pits on the design come from the light pass's DEPOSITION, not from its erosion.**
- **T-all**, which weights the deposition too, **holds G-pits = 0 and G-rim = 0 in all three worlds**, as C2 does.
  - Its walls are 1 120 / 1 127 / 1 259: 4–5 × fewer than C2's (4 858–6 277), but 2.0–2.4 × P0's.
  - It is the best transition measured, and it was declared a diagnostic: it cannot be selected this round.
- G-levels rises from 7 to 15 / 18 at ×1 under every design-protecting candidate (C2, T, T-all). It stays 3 / 9 at
  ×1.4 p .5.
- **G-drained** falls from P0's 9 781 / 14 366 to 7 176 / 11 118 (T, T-all), still 3–5 × the construction's
  (2 190 / 2 318): the cones.
- **G-sea** at ×1.4 p .5 stays (+117 / +122): the breach.
- **The attribution** (G-ring and G-slope100):
  - **24–31 % of G-ring's violators lie in the head falls' footprint; the majority lies beyond it**, a real defect
    of the ring's edge. C2's beyond count (263–543) shows the frozen ring's walls.
  - G-slope100's steepest window lies in the footprint for half the violators (P0 1 / 1; T 2 / 4, 1 / 2).
- **The cost of T is not isolated.**
  - T's builds are 233 / 227 / 161 s against P0's 191 / 192 / 160 s.
  - The two first T worlds ran under memory pressure (run 1, the IDE's clangd growing).
  - At ×1.4 p 0, T takes +1 s. The distance transform was not timed alone.

**The images for round 4**:
- `f156_ombrage_lac{1,2,11}_{p0,t,ta,np}.png`: Relief → Ombrage, north up, 512² cells (25 km);
- centred on ON's final lakes 1 (3 490, 2 500), 2 (4 019, 2 838) and 11 (1 428, 4 287);
- at ×1 p .5: P0, T, T-all and F155's NP (no light pass).
- On lake 1, P0, T and T-all differ only along the design's edge (a 146 m band). NP shows the bare construction.

## Predictions

**The reviewer's (hypotheses to check; judged):**
- Br, "the criterion holds on every new seed": **refuted**. The lakes also lose the pits the old ramps drained, and
  one seed's world crashes in production.
- T, d_t "3–8 cells": **held** (3, by the fallback; the profile gave 1).
- T, gates:
  - G-pits = 0 on the three: **refuted** (1 792–2 487);
  - G-rim = 0: **held**;
  - walls ≤ 1.5 × P0: **refuted** (× 2.1–2.6);
  - G-drained fails at ×1.4: **held**.
- G-ring, "the majority of its failures lies in the head falls' footprint": **refuted** (24–31 %).
- The cost, "T adds < 2 s": **not determined** (+1 s at ×1.4 p 0; the others under memory pressure).
- Meta: **held**.

**Mine (`f156_predictions.md`):**
- **Br**:
  - P-Br1 (≥ 2 of 3 hold): **refuted**;
  - P-Br2 (no new below-sea cell): **held**;
  - P-Br3 (lost = R under the lakes): **refuted**;
  - P-Br4 (5 000–80 000 raised, segments ±2 %): **held**.
- **T**:
  - P-T1 (the profile < 2, the fallback 2–5): **held** (1; 3);
  - **P-T2 (T fails G-pits through the unweighted deposition, of P1's order; not retained): held**, and the counts
    are P1's exactly;
  - P-T3 (G-rim 0): **held**;
  - P-T4 (T's walls ≤ P0's): **refuted** (× 2.1–2.6);
  - P-T5 (G-drained fails): **held**;
  - P-T6 (T-all G-pits 0, walls ≤ P0's): **held for G-pits**, **refuted for the walls** (× 2.0–2.4);
  - P-T7 (+1 to +5 s): **not determined**;
  - P-T8 (≥ 60 % of G-ring in the footprint): **refuted** (24–31 %); "G-slope100 mostly in it": half;
  - P-T9 (canyons 0 in the zone): **held**.
- **Meta**: **held**.

## What the round says (no decision taken)

1. **The spill breach does not ship.**
   - The lakes lose not only the old ramps' cells but the pits those ramps drained, which the spill anchor now fills.
     The criterion did not foresee that.
   - **One new seed's production world crashes** (F85, in both chains): a defect of the production drainage,
     separate from the breach.
2. **The true transition T fails**, for the reason F144-P1 already showed: the light pass's deposition on the design
   makes the pits. Its walls are × 2.1–2.6 P0.
3. **T-all holds G-pits and G-rim everywhere**, with walls × 2.0–2.4 P0. It is the best candidate measured, but it is
   not eligible: it was declared a diagnostic before the run.
4. **G-drained (the cones) and the ring's edge walls beyond the falls** remain the gorge's construction-born failures.

**Named for round 4 (the author's look, then the decision):**
- the images;
- T-all as a possible candidate, with its walls at × 2.0–2.4 P0 against the rule's 1.5;
- the pause.

## Limitations, stated

1. **Br**: 3 seeds; one crashed; seed 003 not localised. The criterion's R is the lowered cells only.
2. **T**: one seed (PSEED).
   - d_t comes from the fallback; the profile is blind outside the mask.
   - The two runs' worlds were assembled (run 1: seven worlds; run 2: five).
   - T's cost is not isolated.
3. **The memory**: run 1 and the first resume ran out of memory beside the IDE's language servers. I stopped clangd
   (PID 39916, 7.8 GB) and the older rust-analyzer (PID 24132, 4.5 GB), at the author's choice.
4. **The images** are 25 km crops for the eye, not a measure.

## State

**Uncommitted**, awaiting the feu vert.
- `tectonics_c1/valley_construction.rs`: `transition_weight`, `blend_light_pass`, `light_mode` 5 / 6 documented, the
  test.
- `tectonics_c1/production_upscale.rs`: `light_mode` 5 / 6 (gated).
- `tests/f156_round3.rs` (`f156_br`), `tests/f156_loc.rs` (`f156_br_loc`), `tests/f126_coast.rs` (`f156_t`).
- The breach, the cache versions and the guard are unchanged.
- ADR Finding 156; this folder.
