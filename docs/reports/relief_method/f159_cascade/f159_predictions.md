# F159 — my predictions, written BEFORE any measurement (2026-10-09)

**Not blind, declared.** Before writing these I had read:
- the solver's code: the implicit stack, MFD p = 2, the linear diffusion's per-step weight (not scaled; F45), the
  talus, the bank planing;
- F155's cost of one 8192² light pass (≈ 30 s);
- F157-B7's numbers on the témoin at 8192²: planar walls 9 % of the land, R8 terrain 0.088;
- Perron 2009's 6.4 L_c ≤ λ ≤ 12.7 L_c;
- Schott 2024's Table 1 (hundreds to thousands of iterations per level).
- I chose the cascade's design (`f159_declared.md`) myself, so the predictions below are about my own design.

## M — per level (128², 256², 512²)

- **P-M1 (planar walls)**: the cascade's share of land at 28° ± 0.5° is ≤ the reference's at every level. The
  reference (the témoin block-averaged) loses most of its 28° walls in the averaging: both are under 1 % at 512².
- **P-M2 (crest facets)**:
  - the cascade's facet share is below the reference's at 512²;
  - the declared control holds at every level: pyramids ≥ 60 %, paraboloids ≤ 20 % (edited before any measurement,
    when the declaration replaced the construction-vs-S1 control with a synthetic one);
  - for information, the construction's share at 512² (block mean) is above S1's.
- **P-M3 (trunks)**: the trunks > 100 km² move by p90 ≤ 1 cell of the previous level between 128² and 256², and
  between 256² and 512². The max is larger (> 3 cells), from captures.
- **P-M4 (valley spacing)**: λ measured at 2–6 cells at each level, near Perron's band with L_c = 0.4 cell (2.6–5.1
  cells). It scales with the cell: the LOD principle holds by construction of D.
- **P-M5 (R8 terrain)**: the cascade's R8 is ABOVE the reference's at 512². The cell-scale valleys follow the D8 /
  MFD lattice: Schott's « axis-aligned artifacts », with no noise to break them.
- **P-M6 (relief)**: the land relief p50 moves by 5–20 % between consecutive levels. It stays within ±25 % of the
  reference's at the same level, after the 128² calibration of U0.
- **P-M7 (equilibrium)**: 128² reaches the declared equilibrium in 100–400 steps; 256² and 512² often stop at their
  step cap.
- **P-M8 (cost)**:
  - 128² → 512² in **< 60 s** in all;
  - the ESTIMATE extrapolated to 8192² is **10–30 min**, above the 15 min target but under the 1 h stop.

## R — the stop rule

- **P-R1**: the rule does NOT fire at 128² or 256². At 512² it may fire on R8 (not a stop criterion) but not on the
  three stop criteria. **The cascade reaches 512².**

## B — the build

- **P-B1**: the production world and the guard are untouched: the guard 6 / 6, field and lakes, before and after.

## The reviewer's (hypotheses to check)

- "at 512², fewer than half the reference's planar walls, and no crest facet beyond the reference": I expect held,
  but nearly vacuous (both near 0 at 512²);
- "trunks > 100 km² move < 1 cell (p90) between 256² and 512²": I expect held, at the edge;
- "without scaling, the relief changes by > 20 % between two levels": not tested (my design scales D);
- "128² → 512² under 60 s; extrapolation under 15 min": the first held, **the second refuted** (10–30 min);
- meta: held.

## Meta

- At least one of my predictions is refuted.
