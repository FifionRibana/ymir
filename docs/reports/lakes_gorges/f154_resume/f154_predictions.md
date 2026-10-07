# F154 — my predictions, written BEFORE any measurement (2026-10-07)

**Not blind, declared.** Before writing these I had read:
- F144's numbers: the P0–P3 table, the edge cells (P1: 4 936 / 6 242 / 6 367), the 11 tagged falls at ×1, and G-sea
  +112 at ×1.4 p .5;
- F143's C2 (the design frozen through the light pass): G-pits 0 in all three worlds;
- the geology code: 4 of the 6 substratum classes are hard (craton, belt, basaltic volcanic, basement), arc volcanic is
  medium and rift fill soft. The surface layer's lake margin (2 cells) covers every lip;
- the current φ draw: 0 with probability 1/3, else U[0.1, 0.5];
- F152's geology cost (≈ 3 s per world in all);
- the production breach (`flow.rs` 1018–1080) and F141's account of the 292 cells.

## S — the substratum / surface split

- **P-S1**: the final rock grid is bit-identical to F152's. `build_rocks` becomes the surface applied over the
  substratum, with the same arithmetic per cell.
- **P-S2**: the substratum is available at the construction stage once the C1 history and the edifices exist. Its own
  build (the three smoothed masks, the terrain anomaly, the cones) takes < 1 s. The history (a C1 rerun) costs 0.5–2 s
  more if it is not shared with the geology stage.
- **P-S3**: the substratum read on the construction's input differs from the one read on the final field on < 3 % of
  the land cells. The snapping and the cones' relief move little.

## H — falls by rock (×1, p 0.5, 18 bodies)

- **P-H1**: ≥ 14 of the 18 lips are on a hard class, and ≤ 2 on rift fill.
- **P-H2**: the class 1 km and 3 km downstream equals the lip's class for ≥ 14 of 18. The contacts are smoothed at
  ≈ 3 km.
- **P-H3**: rule (i) tags ≥ 14 of 18, with a median height above rule (iii)'s.
- **P-H4**: rule (ii) gives the strong φ (a contrast) to ≤ 3 lakes and tags 4–12 in all. The weak φ of 0.1 still
  passes H_f where D_g is large.
- **P-H5**: rule (iii) reproduces F144: 11 of 18. The re-tabulation reproduces each body's built head fall at its
  own φ to < 1 cm.

## L — P4 (the light pass not applied on the kept bodies' catchments)

- **P-L1**: P4 holds G-pits = 0 in all three worlds. C2 did, and the excluded set contains the footprints.
- **P-L2**: P4's edge cells at the divides are < a third of P1's at F144 in all three worlds (I expect 100–800 each).
  They are ≥ 1.5 × P0's on the same instrument.
- **P-L3**: **P4 is retained** by the declared rule.
- **P-L4**: the cols lie outside the catchments, so the light pass still lowers them. G-levels fails for ≥ 1 body at
  ×1, and G-ring / G-slope100 stay near P0's.
- **P-L5**: G-drained still fails at ×1.4 (p .5 and p 0), at about the construction's count (the cones, F144-Z).
  G-sea at ×1.4 p .5 still fails (+100 to +125).
- **P-L6**: the land removed from the pass is 5–15 % of the land. The light pass's removal falls by about that share.
- **P-L7**: θ (corridors excluded) is within the declared tolerance in all three worlds.
- **P-L8 (cost)**: P4 does NOT reduce the time. The light pass still runs on the whole grid and the catchments are
  restored after it; the extra flow routing adds +3 to +15 s per world.

## Br — the breach's ramp aimed at the overflow level (bench copy)

- **P-Br1**: on ON, the 292 below-sea land cells disappear (0 left).
- **P-Br2**: the change is large. With the ramp anchored at the spill, each residual pit is filled to its spill by the
  mop-up fill instead of being trenched:
  - > 100 000 conditioned cells change;
  - ≥ 5 lakes change (footprint or level beyond 1 m);
  - the number of river segments moves by > 1 %.
- **P-Br3**: on the gorge P4 at ×1.4 p .5, G-sea's inlet disappears (≤ 5 cells left). G-drained on the eroded field is
  unchanged by definition, since it is read before the breach.

## Meta

- At least one of these predictions is refuted.
