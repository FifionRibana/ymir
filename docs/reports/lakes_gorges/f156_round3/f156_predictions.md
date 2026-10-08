# F156 — my predictions, written BEFORE any measurement (2026-10-08)

**Not blind, declared.** Before writing these I had read:
- **F155's Br** on seed 1's six states:
  - the lakes LOSE exactly the old ramps' cells in the four C states (témoin: 1000016 −553), and none gains;
  - livré permutes ids and keeps 6 294 coastal cells below the sea.
- **F143's C2** (the design frozen exactly): G-pits 0 in three worlds.
- **F144**:
  - P1 (a floor on the design, deposition allowed): **G-pits 1 792 / 2 330 / 2 487**, G-rim 0;
  - P2 at d_t = 49 m gave P1's numbers;
  - the d_t instrument was blind (C2 and P0 share the light pass outside the mask).
- **F154's P4** and **F155's NP**: G-pits 2–18 and 0; G-rim 0 under NP.
- The light pass's code (`production_upscale.rs`): P2 weights only the LOWERING; deposition passes unweighted.

## Br — new seeds, the production state

- **P-Br1**: the criterion holds on ≥ 2 of the 3 seeds. My best guess: all 3 hold, and the breach goes to production.
- **P-Br2**: no new below-sea land cell on any seed. The spill anchor never starts deeper than the floor.
- **P-Br3**: the lost cells equal the old ramps' cells under the lakes exactly on ≥ 2 seeds; elsewhere the gap is
  < 5 % of them.
- **P-Br4**: raised cells 5 000–80 000 per seed; the river segments within ±2 %.

## T — the true transition

- **P-T1 (d_t)**: the declared profile gives d_t < 2 cells again. Outside the mask C2 and P0 share the light pass, so
  their difference is ≈ 0 from the first cell. The declared fallback (the p90 edge step at 28°) gives **2–5 cells**
  (100–250 m).
- **P-T2**: **T FAILS G-pits** in ≥ 1 world. Its deposition is unweighted (the brief weights the erosion), and
  deposition on the rings is what left P1 its pits. T's G-pits are of P1's order (10²–10³). **T is not retained.**
- **P-T3**: G-rim = 0 in the three worlds (the rings are in the mask, w = 0).
- **P-T4**: T's walls at the mask edge ≤ P0's in all three, well under 1.5×.
- **P-T5**: G-drained still fails at ×1.4 (the cones).
- **P-T6**: the diagnostic T-all (the whole change weighted, deposition included) holds G-pits = 0 in the three worlds,
  as C2 does, with walls ≤ P0's.
- **P-T7 (cost)**: T adds +1 to +5 s per world (the distance transform on 8192²).
- **P-T8 (attribution)**: ≥ 60 % of G-ring's violators lie within the head falls' footprint (Chebyshev ≤ 2 of the
  col or of the fall's foot). The steepest 100 m of G-slope100 starts within it for most outlets.
- **P-T9**: canyons stay 0 in the mask zone widened by d_t.

## Meta

- At least one of these predictions is refuted.
