# Finding 145 — predictions, written 2026-10-04 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–144 and the reviewer's hypotheses, in particular:
- **F144-Z**:
  - 94–99 % of G-drained's violators at ×1.4 are laid by a cone outside their LINE catchment (the first line met down
    the D8 path and its parent chain);
  - 36 % of the world's carved cells too (an upper bound mixing defects and legitimate divide walls).
- **The code of `carve_diag`, read for Part Q before predicting**:
  - a cell's sample is the NEAREST, by a Dijkstra-like propagation from each sample's own cell;
  - a front only stops where the cone rises above the terrain;
  - a minimum is taken against the 8 neighbours' samples of OTHER lines;
  - the témoin has no confluence band, no wall profile, no rim.
- **F127-C's four canyons live in B2 → A_c (a_min = 0.1 km²), not in the témoin** (canyons 0 in the témoin since
  F142). Their dams are tributary FLOORS laid higher at the col, i.e. RAISED cells.

## Predictions

- **P-M1** (the basin is the D8 outlet basin, coarser than F144-Z's line catchment):
  - **the foreign-laid cells are 8–20 % of the carved cells** (below F144-Z's 36 %, because two tributaries of one
    river share a basin);
  - **most of them are deepened (Δz < −1 m)**: the nearest sample is the foreign line's, so its cone lies lower than
    the own lines' at that cell;
  - **the real defect is 5–15 % of the carved cells** (against the reviewer's 3–10 %, overlapping);
  - **less than 40 %** of the foreign-laid cells are within 2 cells of a divide with |Δz| ≤ 1 m (against the
    reviewer's > 60 %).
- **P-M2**: **at most 1 of F127's 4 dams** is a real-defect cell. A dam is a raised col, and the defect is a deepened
  cell (against the reviewer's ≥ 3).
- **P-M3**: divides are lowered by foreign cones (hundreds to thousands of cells, p50 > 5 m), and **at least one
  artificial capture** exists (with the reviewer). **The stop rule does not trigger.**
- **P-C1** (measured on the construction; no production code, so no full world):
  - canyons stay 0 (with the reviewer);
  - F127's dams do not change, because they are same-basin tributary floors (against the reviewer's "decrease");
  - **walls across the divides rise by more than ×2** (with the reviewer);
  - captures fall to ~0;
  - θ within ±0.01 of the témoin;
  - the bench's filtered carve costs < +5 s per world. The labels are free in production, since `compute_flow`
    already returns basins.

## Meta

My disagreements with the reviewer:
- M1's legit share;
- M2;
- C1's dams.

**At least one of my predictions is wrong.**
