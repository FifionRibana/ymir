# Finding 151 — predictions, written 2026-10-06 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–150 and the reviewer's predictions, in particular:
- F150 on the témoin: 8 plates at init and 2 at the end (so ~6 merges, plus any rift splits), 154 convergent and 0
  collision cells at the end, a craton of 57 cells, 832 continental cells;
- C-3: hard 90.8 %, rift-soft 4.9 %, volcaniclastic 4.3 % of the land;
- the code: `run_with_closures` gives its callback the state but not the kinematics, which change at each merge and
  split (`time_loop.rs:848–851`). Classifying a past boundary needs them, so F adds a read-only observer.

## Predictions

- **P-F (the fossil belts)**:
  - over the 300 steps, the cells that were EVER convergent are **600–1 500** of 4 096 (the boundaries move and merge);
  - the cells that were ever **collision** (C-C convergent) are **50–300**, against 0 at the end;
  - the merges leave **≥ 4 sutures**, at least half on continental crust;
  - the field is bit-identical, at < 1 s.
  - With the reviewer (≥ 3 sutures, ≥ 100 cells), larger.
- **P-S (the rock classes, a proposal measured as a share of the land)**:
  - the default class (undifferentiated basement) covers **35–55 %** of the land. Fossil belts plus loose recent
    deposits (valley floors, lake beds, slope feet) will take a large share;
  - against the reviewer's > 50 %, uncertain.
- **P-export** (12 chance grids, u8):
  - raw at full resolution 12 × 67 MB = 805 MB;
  - at ¼, **12 × 4.2 MB = 50 MB raw**, and **≤ 10 MB compressed** (PNG / deflate: the chances are mostly 0 with
    patches);
  - at full resolution compressed, 50–150 MB.
  - With the reviewer (< 50 MB at ¼).
- **P-viz**: the guard reads 6/6 with the new default boxes, since they are exactly the guard's literal. With the
  reviewer.
- **P-cost**: the observer and the recording cost < 0.5 s over the 300 steps.

## Meta

My disagreement with the reviewer: the default class's share. **At least one of my predictions is wrong.**
