# F144 — instruments, declared BEFORE the measurements (2026-10-04)

**The worlds**: the témoin, `GorgeRetreat::v5` (+ `light_mode` per candidate), k = F121_AGE_K × age. The
construction is `skeleton` + `carve` + `gorge_plain` on S1. The final world is `build_world` (d) + the run_hd tail,
isolated twice in `catch_unwind` as in F142 / F143.

## S, B, H, Z (bench `f144_p`, on the constructions)

- **S** (×1, p 0.5): the bodies v5 keeps against F143's list of 14 (their lowest cells), and v4's proxy (18). For
  every input-lake body of v3: in_depression, in_sea, touches_depression.
- **B**: at ×0.7 / ×0.85 / ×1 / ×1.2 / ×1.4 (p 0.5), ×1 (p 0, 0.25) and ×1.4 (p 0), every body the bound raises:
  its B, its spec level and its level. **Lake 4** (lowest cell (2237, 3545)) at each age: L_floor, B and the base's
  kind.
- **H** (×1, p 0.5): per body:
  - A, D_g, φ;
  - the capped head fall and the cap;
  - S and H_f;
  - tagged or not.
  - G-tag on the construction, the old instrument against the corrected one.
- **Z** (×1.4, p 0.5 and p 0): **a line's catchment** = the cells whose skeleton D8 path first meets that line or a
  line it joins downstream:
  - the first line met along the D8 path (`Skeleton::line_owner`);
  - then its chain of parents (`line_parent`).
  - A cell is laid **outside its catchment** when the line of `carve_diag`'s laying sample (the band's sample when
    the confluence clause decided) is not in that chain.
  - The shares among: G-drained's violators; the drained bowls' footprints; the cells upstream of the drained bowls;
    every carved cell of the world. No change is made.

## P — phase 1 (bench `f144_p`)

- **The worlds**: ×1 (p 0.5), ×1.4 (p 0.5), ×1.4 (p 0); P0, P1, P2, P3 each.
- **d_t (P2)**, measured at ×1 before P2 runs. On the final worlds of P0 and C2 (v5 + `freeze_design`), for each
  distance band j = 1…20 cells outside the design mask: the share of cells whose steepest drop to an 8-neighbour
  exceeds 28°. **d_t = the first band where C2's share exceeds P0's by less than 0.5 percentage point**, times the
  cell size (provenance: the width of C2's walls). It is used at all three ages.
- **The gates, per world, on the final world** (F142's, with the declared changes):
  1. **G-bodies**: each kept body with r_lake < 2 is matched 1:1 to a final lake of the gorge world (the body's
     share > 50 %, distinct lakes); each drained body has no final lake over > 50 % of it.
  2. G-pits = 0.
  3. G-levels: every undrained kept body within ±1 m.
  4. Drained bodies absent at the end (all).
  5. G-rim = 0.
  6. G-ring = 0 (on the eroded world).
  7. G-slope100 = 0.
  8. G-drained = 0 (on the eroded world).
  9. **G-tag = 0, corrected**: a two-cell drop on the gorge path beyond the lip exceeding S × the span by more than
     H_f.
  10. G-sea: exactly ON ×1's set.
  11. F38 does not fire.
  12. Canyons = 0.
  13. **θ on the eroded world, without the links touching a gorge-corridor cell** (the path cells, radius 2),
      inside **OFF ×1's eroded CI** (measured in the bench, same links rule).
  14. **Edge cells = 0**: a design-mask cell with an outside 8-neighbour more than 28° lower, on the eroded world.
  - Also reported, not gated: the removal and the raise (km³; inside the kept bodies), the plain's volume, the
    tagged falls, G-area, the coast, the time of each world's build.
- **The selection rule, written before the measurement**:
  - the candidate holding the most of these 14 gates summed over the three worlds;
  - at a tie, the one with the fewest G-ring cells (summed);
  - **if no candidate holds G-pits in all three worlds: report and stop (no phase 2).**
- **Phase 2** (bench `f144_p2`, written after phase 1): the selected candidate on F142's 8 worlds, with the same
  gates.
