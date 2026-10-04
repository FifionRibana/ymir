# Finding 146 — predictions, written 2026-10-04 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–145 and the reviewer's hypotheses, in particular:
- **Q, read before predicting**:
  - `rivers.json`'s segments are D8 cell chains (`flow.rs`, `extract_rivers`), split at confluences, with
    `upstream` / `downstream` links, `profile_m`, `catchment_km2`, `discharge_m3s`, `width_m` and `strahler_order`;
  - **the staircase is born there**: points at cell centres, steps of 45° multiples;
  - **Ymir knows no Living Landz hex grid**: only `km_per_cell` "≈ one hex" in `container.rs` and a bench-side
    assumption of a 1 km-edge hex "for the author to rescale".
  - **The reviewer's Q.3 is therefore read, not predicted.**
- F124's planar walls: the construction's 28° walls, carved and not floor.

## Predictions

- **P-L** (Chaikin with fixed endpoints, then the valley constraint):
  - **before the constraint**, 1–5 % of the smoothed length stands more than the tolerance (1 m) above the bed, on the
    small rivers whose floor is one cell wide (the corner cut climbs a 28° wall);
  - **after the constraint, < 0.5 %** (the reviewer says > 99 % within; I agree after the constraint, not before);
  - ridge crossings: 0 after the constraint, a handful before;
  - lateral deviation from the D8 trace: p50 ≈ 0.15 cell, p99 ≈ 0.5 cell.
  - **Crossings between two smoothed rivers away from a confluence: ≥ 1 before any repair**, where two rivers run one
    cell apart.
- **P-P** (k = 2 cells, L = 2 km): parallel pairs exist (tens to hundreds):
  - **mostly < 10 km² of drained area** (with the reviewer);
  - **more than half of their length on the construction's planar walls** (F124) or volcanic edifices (with the
    reviewer).
- **P-A**: most rivers by count are short headwater stems. Length and drained area are correlated, but long-and-small
  rivers exist; they are the author's "long, even small" case.
- **P-cost**: building `rivers_ll` adds < 5 s per world (the reviewer says < 10 s).

## Meta

My disagreement with the reviewer: L before the constraint. **At least one of my predictions is wrong.**
