# F161 — my predictions, written BEFORE any measurement (2026-10-09)

**Not blind, declared.** Before writing these I had seen:
- everything F160 measured: the transposed Schott erosion inert, the band RMS of the témoin and of the cascade, k ×100
  lowering the land;
- the Corsica grids' land area, peak and mean at each cell size, printed by the preparation script (8 700 km²; peak
  2 661 → 2 367 m from 49 m to 1 563 m cells; mean ≈ 572 m).

None of the Corsica instruments and nothing of N1 has been run. The design (`f161_declared.md`) is mine.

## The Corsica reference

- **P-C1**: the slope p50 / p90 at 391 m cells is ≈ 0.20 / 0.45 (11° / 24°), and it rises by ≈ 15–25 % per halving of
  the cell (a self-affine surface with H ≈ 0.75–0.85).
- **P-C2**: the band RMS (2–4 cells of wavelength) at a fixed wavelength in km is independent of the cell size within
  ±10 %, so the « spectrum in km » is one curve.
- **P-C3**: Corsica's crest facets are 10–20 % of its crests at every cell size; its walls at 28° are ≤ 1 %.
- **P-C4**: Corsica's terrain R8 is at its noise floor (no axis preference).

## N1 (n = 1, the area capped, the smooth retargeting, k calibrated)

- **P-N1**: k_L calibrates to ±20 % at every level within ≤ 8 trials, and **k_L rises** at the finer levels (×2–5 per
  level).
- **P-N2**: the valley spacing in km falls from level to level (a new generation of valleys), and it stays within ×1.5
  of Corsica's at 512² and 1 024². At 2 048² it falls outside: our continent's hillslopes are longer.
- **P-N3**: the restricted drift is < 0.5 % per level with the retargeting; the peak stays at 2 850 ± 5 %.
- **P-N4**: the slope p90 stays within ×1.5 of Corsica's at 512² and 1 024²; the p50 is below Corsica's by more than
  ×1.5 at 512² (our lowlands are flatter).
- **P-N5**: the inherited facets fall under 1.5 × Corsica from 1 024².
- **P-N6**: R fires first on the slope p50 (at 512²).

## The variants

- **P-V1**: without the retargeting (same k) the peak falls by > 20 % by 2 048², and the drift exceeds 2 % at the
  first level.
- **P-V2**: without the area cap the trunks move more, but the p90 stays ≤ 1.5 cells (the implicit solver grades the
  big rivers once, it does not wander them).
- **P-V3**: the talus at a quarter changes the slope p90 by < 5 % and divides its time by ≈ 4.

## Cost

- **P-T1**: 512² → 2 048² for one N1 chain in < 3 min without the calibration; the calibration itself (≤ 8 trials per
  level) takes 5–15 min.
- **P-T2**: the ESTIMATE to 8 192² exceeds 30 min, with the erosion the largest part (the implicit solver's priority
  flood per iteration).

## The reviewer's (hypotheses to check), my reading

- « k calibrates to ±20 %, k rising at the fine levels »: held.
- « the valley spacing falls level to level and stays within ×1.5 of Corsica to 1 024² »: held.
- « the drift < 0.5 % with the retargeting; without it, at the same k, the peak falls by > 20 % »: held.
- « without the cap, the trunks' p90 exceeds 1.5 cells at one level at least »: **refuted**.
- « the inherited facets under 1.5 × Corsica from 1 024² »: held.
- « the talus at a quarter changes the slopes by < 10 % and divides its time by 4 »: held.
- « cost to 2 048² under 3 min »: held without the calibration, refuted with it.
- meta: held.

## Meta

- At least one of my predictions is refuted.
