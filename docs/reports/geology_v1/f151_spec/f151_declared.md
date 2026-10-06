# F151 — instruments and methods, declared BEFORE the measurements (2026-10-06)

**Before / after**: the guard 6/6 (`ymir-viz f123_viz_guard`), lib and viz tests, `cargo check --workspace`.
"Before" = F150's after-checks (`267f868`).

## 0 — the viz's two small changes

- **The defaults**: `volcanism`, `lithology` and `fracture` checked at startup; `infiltration` stays off. A permanent
  test asserts that the default `WorkspaceState`'s closures are the guard's literal: C-2 / C-3 / C-3b on, H-1 off.
- **The panel**:
  - « Type de plaque » (continental / océanique) and « Craton » (oui / non) from `HdResult.tectonic`, sampled at the
    cell with the overlay's own mapping (round to the nearest coarse cell);
  - « Épaisseur crustale / S̃ » from a new `crust_s` in `CoarseTectonicLabels` (the debug labels, read-only).
    Non-dimensional: no km scale exists in the code.
- **The gate**: guard 6/6 unchanged.

## F — the fossil-belt prototype

**The observer**: `run_with_closures_observed`, whose callback also receives the kinematics. `run_with_closures`
delegates to it.
- **This is a production-code change of the signature only**, read-only: the callback cannot mutate.
- Declared because the brief asks for no production code in F. The kinematics change at merges and splits, so without
  them a past boundary cannot be classified.
- Its gates:
  - the final state is bit-identical with and without the recording (`s`, `age`, `plate_id`, `plate_type`,
    `cratonic_mask`);
  - the guard 6/6 holds.

**Recorded per coarse cell over the 300 steps** (the témoin's run, the benches' config), with `derive_tectonic_labels`
at each step:
- the number of steps convergent (any), collision (C-C), subduction upper plate, and divergent;
- the last step of each.

**The sutures**: at each step where a plate id disappears (a merge), the previous step's cells of the vanished plate
that touch the plate that absorbed it. Recorded: the step, the two plates, the number of cells, and the share on
continental crust.

**Reported**:
- the counts;
- the maps (north up): ever-convergent, ever-collision, sutures, and the proposed fossil belt;
- the overlap with the craton and with the current boundaries;
- the recording's cost (the observed run against the plain run).

**The time scale**: ≈ 0.67 Ma per step, from `state.rs`'s "lag is ~0.67 Ma" comment. A PROXY, to check.

## S — the measurements the spec needs

**The default class's share**: the spec's classes are not built this round. The bench measures the land shares of
their PRIMARY sources:
- craton;
- the proposed fossil belt and the current belts (convergent within 2 coarse cells);
- volcanic discs;
- rift-soft;
- loose deposits: the construction's valley floor, or a lake footprint, or a slope < 0.02 with accumulation ≥ 1 km²;
- evaporites: an endorheic lake footprint.
- The default = the land in none of them.
- These thresholds are this measurement's, not the spec's.

**The export sizes**: 12 PROXY chance grids (the rules are not built), each a declared field mask × a distance-decayed
chance (100 inside, −10 per km outside to 0), quantised to u8. For each, at full, ½ and ¼ resolution:
- the raw size;
- PNG (8-bit grey, the `image` crate's deflate).
- The PROXY is labelled: the real grids come in round 3. Their entropy may differ.
