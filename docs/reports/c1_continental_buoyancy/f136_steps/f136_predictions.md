# Finding 136 — predictions, written 2026-10-02 BEFORE any grep, reading or measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–135 and on the reviewer's predictions,
and in particular on:
- the K2 profiles (OFF skeleton path; the steps of ~300 / ~280 / ~900 m; lake 11's col 1.69 km off that
  path, lake 2's 0.53 km, lake 1's 0.07 km);
- K1's cuts (ring and at-col);
- F135-M2's link counts: 31 593 based, 26 418 also union, so 5 175 in the based set only;
- the code seen in F134 / F135: `LithologyConfig { soft_multiplier, volcanic_multiplier, rift_age_threshold }`,
  `production_k_field` (the C-3 erodibility multiplier), `debug_labels::LithoClass {Hard, RiftSoft}`, the edifice
  basal discs (the viz's "volcaniclastic" overlay), `RiverSegment.points: Vec<(u32, u32)>`.

## G
- **P-G1**: nothing in the code joins two base levels with a slope. The construction lays each reach from its own
  base, and the breach only makes the profile monotone. The ADR names outlet incision / knickpoint retreat only as
  literature or as an unbuilt item (H-2 / F114 among them, I do not remember F114's content: blind there).
- **P-G2** (non-blind, from the code): **a lithology / erodibility field EXISTS**. It is the C-3 K multiplier (hard
  1, rift-soft, volcaniclastic), from the coarse tectonics and the edifice footprints. It is not a geological map
  of rock types. A volcanic mask exists (the edifices' basal discs). Against the reviewer's "no lithology field".
- **P-G3**: the rivers' export has no tagged point list. Segments carry points and per-segment scalars; a
  per-point tag needs a new field.

## K5
- **P-K5**: on the ON water's own path (the ON world's D8 from the lake's real outlet):
  - the step exists for lakes 1 and 2, within ± 30 % of K2's heights;
  - for **lake 11**, whose col lay 1.69 km off K2's path, the ON water leaves by another way and the step is
    less than half of K2's ~280 m.

  This half-disagrees with the reviewer (three lakes, ± 20 %).

## K6
- **P-K6**:
  - 6–12 steps > 200 m and ≥ 3 > 500 m among the present lakes;
  - most lakes ≥ 5 km² have a step > 50 m;
  - the canyons instrument sees none (with the reviewer on "none"; ≥ 8 > 200 m is in my range, not certain).

## K7
- **P-K7** (against the reviewer):
  - the lake's lowering (15–48 m) is NOT the notch's depth;
  - the notch is the outlet valley's floor laid from the downstream base, by the construction's nearest-sample
    propagation (not the cross-line minimum, not the breach), and it is hundreds of metres deep;
  - the lowering is set by the **input lake's own bed**: the new sill is a bed cell between the deep part and
    the notch. The construction does not carve a based bed, so lowering < notch, by far.

## L
- **P-L**:
  - the 5 absent lakes are emptied by the construction's notch where the bed drains to it (≥ 3 of 5 have no bowl
    left in the eroded world). The others are emptied by the balance or the below-sea cleanup;
  - the floors of the ones the terrain empties are the S1 bed (FBM), not a plain: < 30 % of the footprint under
    0.1 % slope;
  - **they carry < 50 % of the 5 175 links**; the rest are the margins of shrunk lakes and closed depressions
    that are not final lakes. Against the reviewer's > 80 %, and against its "flat plain".

## K8
- **P-K8**:
  - hanging junctions > 50 m are FEWER than outlet steps > 200 m. The construction keeps junctions continuous
    (the cross-line minimum and the shared χ), so ≤ 5 hang by > 50 m. Against the reviewer;
  - few are F127-type dams; most of the few are clean breaks;
  - coastal-wall mouths above the sea: ≤ 10, under 20 m;
  - in the volcanic mask, the max 2 km slope of the trunk reaches is higher than outside.

## Meta
- Disagreements with the reviewer: P-G2, P-K5 (lake 11), P-K7, P-L, P-K8. **At least two of these five are
  wrong.**
