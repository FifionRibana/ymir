# F155 — my predictions, written BEFORE any measurement (2026-10-07)

**Not blind, declared.** Before writing these I had read:
- **F121's table**, C1 bare (no light pass) against C2 /10 (with it), measured at F121 and since moved by later rounds:
  - relief p50 519.0 against 457.1 m;
  - R8 terrain 0.0798 against 0.0608;
  - σ p50 5.36 against 6.34 m;
  - lakes 16 and 16;
  - build 78 s against 105–110 s.
- F121's title: « the texture left to a light incision ».
- **F143**:
  - C2, the design frozen through the light pass, holds G-pits = 0 in three worlds and G-levels 12 / 14 at ×1;
  - the light pass lowers θ (0.493 on the construction against 0.423–0.436 eroded).
- **F154**:
  - P4 shows the light pass makes part of G-drained;
  - the spill-anchored breach on ON: 292 → 0 below-sea cells, 22 863 raised, no lake moved by > 1 m / 5 %.
- **The code**, read this round:
  - the light pass carries the absolute slope floor A1+B2 (`production_upscale.rs`);
  - the climate and the HD drainage bundle are keyed on the ERODED key, not the conditioned one;
  - the skeleton's own breach (`valley_construction.rs:642`) is `breach_monotone`, the same function.

## Br — blind confirmation on the guard's other five states

- **P-Br1**: no state changes its lake COUNT. The lake mask's hash (the first half of the guard's fingerprint) is
  identical in all six. **The stop rule does not fire.**
- **P-Br2**: the lake LIST's hash (levels and areas, f32 bits; the second half) changes in ≥ 1 state, through the
  climate read on the raised cells.
- **P-Br3**: the below-sea land cells the breach leaves go to 0 in every state that has any, and states at 0 stay
  at 0.
- **P-Br4**: 5 000–60 000 cells raised per state, 0 lowered. The river segments move by < 2 %.
- **P-Br5 (production)**: `ALGO_BREACH` alone does not suffice. The climate and the HD bundle are keyed on the eroded
  field, so `ALGO_CLIMATE` and `ALGO_HD_DRAINAGE` must move too. The lake guard's digests then change, while its
  hashes stay equal where P-Br1 holds. The field guard is unchanged.

## F — the soft-lip rule

- **P-F1**: on the témoin no kept lip is soft (F154: 17 basement, 1 basaltic). The gorge is bit-identical with the
  rule on, and every φ is the current draw.

## P — what the light pass brings

- **P-P1 (grep)**:
  - F121 added it as a texture: it breaks the construction's flat floors and planar walls and lowers the
    anisotropy. That is a geometric reason, not a hydrological one.
  - **But it also carries the absolute slope floor A1+B2 (F109)**, a drainage closure. Without the pass, that
    closure goes too.
- **P-P2 (the production world without the pass)**:
  - relief p50 +8 to +15 %;
  - R8 terrain +20 to +40 %; σ p50 −10 to −20 %;
  - θ rises to 0.47–0.50, outside OFF's eroded CI;
  - canyons 0;
  - lakes +10 to +40 % in number;
  - the river segments within ±10 %; Strahler ≥ 2 within ±10 %;
  - the parallel bundles rise;
  - the biomes change on 1–5 % of the land;
  - **the cost falls by 20–40 s per world**.
- **P-P3 (what appears without it)**:
  - the planar-wall share (slopes at 28° ± 0.5°) rises ≥ 3×;
  - the sharp crests rise;
  - the comb tile's R8 rises.
- **P-P4 (the gorge without the pass, ×1 p .5, ×1.4 p .5, ×1.4 p 0)**:
  - G-pits = 0 in all three;
  - G-rim holds;
  - G-levels fails for ≥ 1 body at ×1 (C2 held 12 of 14);
  - G-drained still fails at ×1.4, near the construction's count (2 190 / 2 318, within ×1.5);
  - G-ring fails (the head falls);
  - G-sea holds where the production breach is the spill one.

## Meta

- At least one of these predictions is refuted.
