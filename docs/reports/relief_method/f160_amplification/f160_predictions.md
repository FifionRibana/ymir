# F160 — my predictions, written BEFORE any measurement (2026-10-09)

**Not blind, declared.** Before writing these I had read:
- Schott 2024 (§3–§6, Tables 1, 2, 4) and its published code (shaders and the preset sequence);
- F159's measurements: 128² at equilibrium in 36 steps, the talus inert at 128², the trunks' p90 1.41 / 0.71 cells;
- `f160_declared.md`, which I wrote myself, so these predictions are about my own design.

Nothing of the amplification had been run.

## The physics level (P128, P256)

- **P-P1**: the U₀ calibrated on the peak (2 850 m) is 0.6–0.8 × F159's (which matched the mean). Both levels are at
  equilibrium, P128 in < 100 steps and P256 in < 150.
- **P-P2**: P256's peak is reached with its own U₀ (calibrated at 256²). Its mean land altitude is within 15 % of
  P128's.

## The amplification (S, moyen budget, from P128)

- **P-A1**: before retargeting, the peak never exceeds the physics level's (the erosion lowers, the deposition fills
  valleys and pits). After retargeting it is within ±5 % of 2 850 m.
- **P-A2**: the restricted drift is negative at every level (−1 % to −6 % of the mean land altitude), under the 10 %
  stop.
- **P-A3**: the carved volume per level is 5–40 % of its theoretical bound (the slopes are mostly under s_ref), and it
  halves with each level (± 30 %).
- **P-A4**: the trunks ≥ 100 km² move by a p90 ≤ 1.0 cell at every transition (the area cap keeps the big rivers from
  re-carving).
- **P-A5**: the planar walls (28°) stay under 20 cells or under the témoin at 256² and 512². The facets stay under the
  témoin at 256² and 512² and go above it at 1024².
- **P-A6**: the detail per band (the band-pass RMS) is below the témoin's in the two finest bands of every level, by a
  factor 1.5–4.
- **P-A7**: S+ρ lowers the terrain R8 and the coherent-window share by 10–20 % against S: **less than the reviewer's
  quarter**.
- **P-A8**: the faible / fort budgets move the finest-band RMS by ×0.6 / ×1.6, not ×0.5 / ×2. The deposition and the
  talus do not scale linearly.

## R and cost

- **P-R1**: the best variant reaches 1 024² without firing R. If 2 048² is measured, R fires there on the facets.
- **P-C1**: on the CPU, the talus dominates the cost at 1 024² and 2 048² (it has the most iterations). The
  per-process extrapolation to 8 192² exceeds 15 min.

## The viz

- **P-V1**: the guard 6 / 6 holds before and after (field and lakes).

## The reviewer's (hypotheses to check), my reading

- « the peak never exceeds the physics level's before retargeting »: held;
- « the drift negative and under 10 % per level »: held;
- « the carved volume under the bound and decreasing with resolution »: held;
- « S+ρ lowers R8 and the comb by at least a quarter »: **refuted**;
- « the facets stay under the témoin with the noisy talus »: refuted at 1 024²;
- « the trunks under 1.5 cells (p90) »: held;
- « on the CPU the talus dominates at the fine levels; the extrapolation exceeds 15 min »: held;
- meta: held.

## Meta

- At least one of my predictions is refuted.
