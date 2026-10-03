# Finding 141 — predictions, written 2026-10-03 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–140 and the reviewer's predictions, and in particular
on:
- **F140's diagnosis**: body 4, A (col, construction) 1 571 km² against F139's 298; 0.6 of 52.4 km² held
  pre-breach; 2 821 land cells ≤ sea after the breach in GORGE against 292 in ON; a cell at 91.1 m, 78 m under body 4's
  input floor, at (4148, 3674).
- **The breach's code, read in F140** (`flow.rs`, `breach_monotone_protected`): when a pit `nb` is reached from a
  higher outlet `ci`, the loop lowers `ci → backlink → …` to `height[nb] − EPS`, then `− 2·EPS`, …. It stops at a base
  or at the first cell already lower than the ramp. **Br's central question is therefore half-read, not blind.**

## Predictions
- **P-D**: deviations 1 (A) and 3 (φ's seed) touch tabulated quantities.
  - Deviation 1: F139-R's r_lake and emptying ages; F138-T1's S_loi; F139-C2's m_fit; T4's H_f.
  - Deviation 3: T4's φ per lake, so the falls' counts and heights.
  - Deviation 2 (L_bed from the input's minimax instead of today's level) touches F138-T3 / F139-R's L(r); it
    should agree on the D8 lakes (F136-K7: today's level is the bed sill).
  - With the reviewer ("at least one other than 1"); I expect two others (2 and 3).
- **P-A**:
  - (iii), the inflow into the construction's input body, and (ii), F139's instrument, agree within ×1.5 for ≥ 15 of
    the 22 lakes (with the reviewer);
  - (i), the col's area, differs from (ii) by > ×2 for ≥ 4 (with the reviewer), lake 4 among them;
  - the cause for lake 4: the input footprint intercepts a large river that the final lake does not (the river runs
    along the bowl's edge).
- **P-Br**:
  - **the ramp starts at the PIT's floor, not at its spill**, and is carved towards the base along the
    priority flood's backlinks. It stops at the first cell lower than the ramp, which is not the pit's spill side;
  - the ON cells cross the sea because their pits' floors are low (≤ ~10 m), near the coast or the below-sea basins;
  - the ON 292 are covered by final lakes at ≥ 90 % (with the reviewer), mostly the below-sea basins' merged lakes;
  - in GORGE the 2 821 cross the sea from the drained bowls' pits, which are not held.
- **P-Z** (against the reviewer, half): the cell is laid by **the light pass**, not by the construction. The
  construction's lake cells are based on L(r) = 178.7 ≥ L_floor and their cones hang from that. Once the bowl is open
  and its pits are not held as depressions, the stream power's A1 no longer guards them, and the incision relaxes
  cells onto low receivers. **Second choice**: a cone of the outlet line's downstream samples reaching into the bowl.

## Meta
- Disagreements: P-Z (the light pass against the construction's cone). **At least one of my predictions is wrong.**
