# F154 — instruments, declared BEFORE the measurements (2026-10-07)

**The worlds**:
- **The témoin**: C2 /10 col, seed PSEED.
- **ON**: `on_at(1.0)`, the témoin with `LakeBase::InputLakesAndBasins`, no gorge (F144's reference for G-sea).
- **The gorge**: `GorgeRetreat::v5(map_a(age), p)` on `on_at(age)`, with `light_mode` 0 (P0) or 4 (P4).
- The construction is `skeleton` + `carve` + `gorge_plain` on S1, as in F144. The final world is `build_world` plus
  the run_hd tail, inside `catch_unwind`.
- Bench `f154_resume` in `crates/ymir-core/tests/f154_resume.rs`.

## S — the substratum / surface split

- **The split**:
  - `rocks::build_substratum` holds what the tectonics and the edifices give: craton, belt, basaltic or arc volcanic,
    rift fill, basement. It is defined on every cell, water included.
  - `rocks::apply_surface` lays over it what the end of the chain gives: water, evaporites, loose deposits.
  - `build_rocks` = `apply_surface(build_substratum(…))`, unchanged in signature.
  - **The substratum reads a relief** (the contacts' terrain snapping, the cones' relief test). At the construction
    stage that relief is the construction's input.
- **The gate**:
  - the témoin's final rock grid, from `build_geology`, hashed (FNV-1a 64 over the bytes) and counted per class;
  - **first at HEAD before the split** (rocks.rs and geology/mod.rs are F152's, unchanged since `f85f431`), then
    after the split;
  - bit-identical, or the gate fails.
- **A permanent test**:
  - on a synthetic world, `build_rocks` equals the surface over the substratum cell for cell;
  - negative control: a lake-margin cell's substratum is not its final class.
- **The stage and the cost**:
  - `record_history` timed alone (it reruns C1);
  - `build_substratum` timed on S1 (the construction's input) and on the final field.
- **The agreement**: the share of the final land cells (not sea, not lake) whose substratum read on S1 differs from
  the one read on the final field.

## H — falls by rock (×1, p 0.5, the 18 bodies of v5)

- **The lip**: the body's col (`GorgeBody::col`).
- **Downstream**: the construction's outlet path from the col (the skeleton's D8, the path the head fall is computed
  on), with the cumulative distance s along it.
  - The samples are the first path cell with s ≥ 1 000 m and the first with s ≥ 3 000 m.
  - "ends" when the path stops first (at the sea, or at another body).
- **The class read**: the substratum on S1, the stage where the falls are designed.
  - Also reported: the substratum on the final P0 ×1 p .5 field at the same cells, and that field's final (surface)
    class.
  - The hardness is `ROCK_CLASSES[..].hardness` (2 hard, 1 medium, 0 soft).
- **The φ values, DECISION** (mine, for the tabulation; the author's to set):
  - **rule (i)**: hard → φ_d = 0.5, medium → φ_m = 0.25, soft → 0;
  - **rule (ii)**:
    - a soft lip → 0;
    - a lip harder than the softest available downstream sample (1 km, 3 km) → φ_strong = 0.5;
    - otherwise (equal hardness, a harder downstream, or no sample) → φ_weak = 0.1;
  - **rule (iii)**: the body's own φ (`gorge_phi`, the reference).
- **The head fall per rule**: `GorgeBody::head_for_phi(φ)`, the construction's arithmetic on the body's frozen path
  and law, with H_cap(A) unchanged.
  - top = level − min(φ·D_g, H_cap); the gorge's slope is slope_for(top); head = level − (top − S·s₁).
  - D_g does not depend on φ.
  - **Approximation, declared**: the other bodies' raises are not replayed for a different φ (bodies are built in
    sequence).
  - **Self-check**: at the body's own φ, `head_for_phi` must return the built `head_fall_m`.
- **Tagged**: head > H_f, the code's H_f = max(1.5·S·100 m, 10 m).
- **Per rule**: the count tagged, the median and the max of the tagged heights, the split by lip class, and the
  per-lake table.

## L — P4 (bench option `light_mode = 4`)

- **The excluded set E**:
  - **The field**: the light pass's input, the construction plus the minimal plain.
  - **The routing**: the light pass's own, `compute_flow(field, FlowConfig { sea_level: 0.5, ..default })`
    (priority-flood fill + D8).
  - **The definition**: E = the cells whose D8 path reaches a cell of a kept body's footprint (`gorge_body_of`),
    footprints included. A path leaving the map leaves E.
  - **Applied**: after the light pass, E's cells are reset to the pass's input (F143-C2's restore, on E instead of the
    design mask). The droplets, the active rims, the bathymetry and the breach run as under P0.
  - **The cols are not in E** (a col drains away from its body).
- **The worlds**: ×1 p .5, ×1.4 p .5, ×1.4 p 0; P0 and P4 each.
- **The gates**: F144's 14, with F144's definitions, plus:
  - **G-area**: the final area of the undrained kept bodies' lakes does not rise from ×1 to ×1.4 (p .5), per
    candidate;
  - **coast**: F144's near-coast ratio, within ±2 % of P0's on the same world;
  - **θ**: on the eroded world without the corridors' links, IN if within OFF ×1's eroded CI widened by **0.005** on
    each side (the tolerance; F144's values sat within 0.002 of the bound).
- **The divide edge** (new):
  - a land cell c ∈ E with a land 8-neighbour m ∉ E where |z_c − z_m| / d > tan 28°, on the eroded world (m);
  - counted on P0 and P4 (same construction, same E);
  - F144's design-mask edge instrument is also reported on both.
- **The land removed from the pass**: E's land cells (S1 above the sea), in km² and as a share of the land.
- **Rule 14**: the removal and the raise against S1 (km³) as F144's extras, in all and inside E.
- **The cost**:
  - `build_world`'s time per world;
  - the E computation timed alone on the construction (`compute_flow` plus the walk), since single timings carry
    ±30 s of noise.
- **The selection rule, written before measuring**:
  - « P4 est retenu s'il tient G-pits sur les trois mondes, et si ses cellules de bord restent sous le tiers de celles
    de P1 au F144. »
  - The thresholds, with P4's divide edge cells: 4 936 / 3 = **1 645** (×1 p .5), 6 242 / 3 = **2 080** (×1.4 p .5),
    6 367 / 3 = **2 122** (×1.4 p 0).
  - Otherwise: report, and build nothing further.

## Br — the breach's ramp aimed at the overflow level (bench copy, no production change)

- **The copy**: `breach_monotone_protected` copied in the bench, with one change:
  - the ramp's anchor becomes `target = z[ci] − EPS`, the level of the cell through which the pit's flood overflows,
    when the priority flood reaches the pit;
  - it was `height[nb] − EPS`, the pit's floor.
  - The mop-up fill is unchanged.
  - **Anticipated, declared**: the ramp then lowers its path by EPS steps only, and the mop-up fill raises the pit to
    its overflow level.
  - **Check**: with the change off, the copy is bit-identical to production on ON.
- **On ON**, the run_hd tail with the production breach, then with the copy (the same tail otherwise):
  - **the below-sea land cells**: S > sea and conditioned ≤ sea, and those covered by a lake;
  - **the conditioned cells changed** by more than 1 cm;
  - **the lakes**, paired by footprint (the best overlap, > 50 % of the production lake). A lake changes when it is
    unmatched, or its level moves > 1 m, or its area > 5 %. Also new lakes;
  - **the rivers**: the segment count, and the river cells (the union of the segments' points) that differ.
- **On the gorge P4 at ×1.4, p .5**:
  - G-sea, added / missing against ON's production set, with the production breach and with the copy;
  - G-drained on the eroded field (read before the breach, so unchanged by definition) and on the conditioned field,
    with each breach.

## Amendment after the run (NOT blind): Br, ON's below-sea cells attributed

- **Why**: the run showed the overflow copy keeps 286 of ON's 292 below-sea land cells. My declared anticipation
  ("the ramp then lowers its path by EPS steps only") does not explain that.
- **The instrument** (bench `f154_br_attr` in `f154_resume.rs`): an instrumented copy of the breach's ramp stage
  records, for every land cell a ramp takes to ≤ sea, its ramp (the pit, the outlet cell ci, the anchor) and its step.
  - The anchors: the floor (production) and the overflow.
  - Reported: the count, the share taken by a ramp, the cells' original height above the sea (p0 / p50 / p100), and
    the largest ramps.
  - **Check**: the ramp stage's below-sea set equals production's (the mop-up fill raises no base).

## Amendment 2 after the run (NOT blind): Br, the ramp anchored at the pit's spill

- **Why**: `f154_br_attr` showed that the declared reading, `z[ci]`, is not the pit's overflow level here. When the
  priority flood reaches a pit of the row along lake 1000016's fringe, its outlet cell `ci` has already been lowered
  by the neighbouring pits' ramps: z[ci] lies 0.1–0.5 m above the floor, while the pit's spill on the pre-breach flood
  is 13.5–16.9 m. **The declared copy did not test the brief's idea.**
- **The second copy** (bench `f154_br_spill`, `f154_resume.rs`): the anchor is `max(height[nb], filled[nb])`, the
  pit's spill on the pre-breach priority-flood fill (the `filled` the breach already receives). The rest is
  unchanged. Check: with the change off it is bit-identical to production.
- **Measured**, as in the declared Br:
  - on ON: the below-sea land cells (covered or not), the conditioned cells changed, the lakes paired by footprint,
    the rivers;
  - on the gorge P4 ×1.4 p .5: the below-sea cells and G-sea against ON's production set (a rebuilt world, same
    config).

