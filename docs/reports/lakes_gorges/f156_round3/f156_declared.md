# F156 — instruments, declared BEFORE the measurements (2026-10-08)

## Br — new seeds, blind

- **The seeds** (never used in this work, nor anywhere in the tests; grep-checked): **20261008001, 20261008002,
  20261008003**.
- **The state**: production, C2 /10 col (`ValleyConstruction::new(F121_AGE_K, Some(0.1))`, S_eq 0.024, volcanism,
  lithology and fracturation on).
- **The framing**: the viz's automatic roll per seed (`c1_coarse_land_report`'s largest-mass centre, rolled to the
  domain's centre, as `run_hd` does without a manual pan); 45° / 40°.
- **The copy**: F154's amendment 2, untouched. The ramp is anchored at `max(height[nb], filled[nb])`.
- **The pairing**, by footprint overlap and never by id:
  - each production lake P is paired with the spill lake S sharing most of its cells;
  - each spill lake S is paired with the production lake sharing most of its cells.
- **The old ramps' cells R**: the cells production's ramp stage (the floor anchor) lowers, recorded by an instrumented
  copy checked bit-identical to production. **R_lake** = R ∩ the production lakes' cells.
- **The criterion (each item, each seed)**:
  1. the lake count is unchanged;
  2. **no lake gains a cell**: no cell of a spill lake lies outside its paired production lake;
  3. **the cells the lakes lose are the old ramps' cells**: lost = the production lake cells not in their paired
     spill lake. **Violated if any lost cell lies outside R.** The equality is reported both ways: |lost \ R|, and
     |R_lake \ lost| (R cells still under a lake, reported, not a violation);
  4. **no land cell is created below the sea**: none ≤ sea under the copy that is not ≤ sea under production.
  - Also reported: the cells raised / lowered (> 1 cm), the river segments' variation, the below-sea land cells under
    each anchor.
- **The stop rule**: one violation, on one seed, and no production. Localise it.
- **Production, otherwise** (in the HD assembly only):
  - `terrain::flow::breach_monotone_anchored(.., RampAnchor)`; `run_hd`'s conditioning and its bench mirror
    (`common::viz_hd_lakes_on`) take `RampAnchor::Spill`;
  - the skeleton's breach becomes an explicit `breach_monotone_anchored(.., SKELETON_RAMP_ANCHOR)` with
    `SKELETON_RAMP_ANCHOR = RampAnchor::Floor`, under a test (the constant, and the construction bit-identical to
    the floor breach);
  - `ALGO_BREACH` 1 → 2 (the conditioned field), `ALGO_CLIMATE` 1 → 2 and `ALGO_HD_DRAINAGE` 12 → 13 (both read the
    conditioned field, keyed on the eroded key);
  - the permanent test `a_spill_anchored_ramp_stays_above_the_sea`, with its negative control (the floor anchor);
  - **the guard**:
    - `bench_field_hashes.json` must not change (the eroded field);
    - `bench_lake_hashes.json` is regenerated (`f134_lake_guard`): its digests (ALGO_HD_DRAINAGE) and its hashes (the
      footprints) change. The diff is reported;
  - the ADR queue: the first defect removed, the second (livré's coastal pits) kept.

## T — the true transition (the gorge, gated)

**d_t, measured BEFORE T is built** (bench `f156_t`, at ×1 p .5, on P0 and C2 = v5 + `freeze_design`):
- **The profile**: the distance from the design mask (`gorge_design_distance_m`, in cells).
  - For bands j = 1…20 outside the mask (and j = 0 on its edge cells), D(j) = the median, over the band's land
    cells, of (z_C2 − z_P0) on the eroded worlds (m).
  - **d_t = the first j ≥ 1 with D(j) < 1 m.**
- **If d_t < 2 cells, the declared fallback**:
  - the step at C2's mask edge: S = the p90, over C2's edge cells (mask cells with an outside 8-neighbour), of the
    largest drop to an outside neighbour (m);
  - **d_t = max(2, ceil(S / (tan 28° · cell)))** cells: the distance a 28° ramp needs to absorb that step;
  - reported with why the profile gives < 2.
- **F144's wall width is reported beside it** (its instrument gave 1 cell) and C2's edge-wall count.

**The rule** (`light_mode = 5`, `light_dt_m = d_t`):
- after the light pass, a cell's LOWERING is weighted by w(d) = smoothstep(min(1, d / d_t)) = x²(3 − 2x), with
  d its distance from the design mask;
- **w = 0 on the design**: the rings, the gorge corridors (radius 2), the plains and the lake beds, i.e.
  `gorge_design_mask`;
- **deposition is not weighted** (the brief weights the erosion).
- `valley_construction::{transition_weight, blend_light_pass}`.
- **A diagnostic, declared, NOT eligible for selection**: **T-all** (`light_mode = 6`), the whole change weighted,
  deposition included.
- **Permanent test**: w(0) = 0, w(d_t) = 1, monotone; on the design the lowering is undone, beyond d_t it is
  kept, and deposition passes (T) or not (T-all). **Negative control**: without the blend the design cell is lowered.

**The worlds**: ×1 p .5, ×1.4 p .5, ×1.4 p 0. P0, C2, T and T-all each (C2 also feeds d_t and the wall comparison).

**The gates**: F144's 14 (F154 / F155's code and tolerances), plus:
- **the walls at the mask's edge**: F144's `edge` (a design-mask cell with an outside 8-neighbour more than 28°
  lower), against P0 and against C2;
- **the canyons in the zone**: the over-dug bodies whose floor cell lies within d_t of the mask;
- **the cost**: `build_world`'s time.

**The attribution of G-ring and G-slope100**:
- **the head fall's footprint** = the cells within Chebyshev 2 of the body's col (`outlet_path[0]`) or of the
  fall's foot (`outlet_path[1]`);
- **G-ring**: each violating ring cell is counted IN or BEYOND that footprint;
- **G-slope100**: the window holding the outlet's steepest 100 m is IN the footprint when it starts within 0.15 km
  (3 cells) of the outlet's col, BEYOND otherwise.

**The selection rule, written before measuring**: « T est retenu s'il tient G-pits = 0 et G-rim = 0 sur les trois
mondes, avec des murs au bord du masque au plus 1,5 fois ceux de P0. Sinon, rapporter, et le tour 4 décide de la
pause. » The walls are F144's `edge` count, per world.

**The images** (round 4):
- the hillshade of `run_hd`'s Relief → Ombrage, north up, 512² cells;
- **centred on lakes 1, 2 and 11**: ON ×1's final lake ids (F139's numbering), their cells' centroid;
- at ×1 p .5, for P0, T, T-all and the gorge without the light pass (F155's NP, rebuilt).

## Amendment after the Br run (NOT blind): the localisation

- **Why**: seed 20261008001 violates item 3. 204 of the 534 lost cells are not an old ramp's cell, so the stop rule
  fires. « Localiser. »
- **The instrument** (bench `f156_br_loc` in `tests/f156_loc.rs`, on each violating seed), per production lake with
  such cells:
  - its paired spill lake, both lakes' type, level and area;
  - the lost cells outside R: their Chebyshev distance to the nearest old-ramp cell, their eroded height, their
    conditioned Δz (spill − production), how many lay in a pre-breach lake, their bounding box.
- **No production change**: the breach stays as it is. T is measured on the current production.

## Amendment after the Br run (NOT blind): the bench's robustness

- **Why**: on seed 20261008002 the run_hd tail panicked in the drainage (F85's fixed-point assertion,
  `drainage.rs:2831`), and the bench stopped before seed 20261008003.
- **The change**: each chain (production, spill) runs in `catch_unwind`. A panic is reported with its chain and counted
  as a violation for that seed, and the bench goes on. The criterion is unchanged.
- The bench is rerun in full (the three seeds), then the localisation on the violating seeds.


## Amendment after T's run 1 (NOT blind): the memory and the resume

- **Why**: run 1 (`f156_t_raw_run1.txt`) aborted at its eighth world, « TA ×1.4 p 0.5 », on a failed 512 MB
  allocation. It had held the ×1 buffers (P0's and C2's eroded fields, the distance, the mask) from d_t's
  measurement to the end.
- **The change**: the buffers are freed once d_t is measured. A resume (`F156_DT_M`, `F156_SKIP`) injects run 1's
  measured d_t (3 cells = 146.484375 m) and skips the seven worlds run 1 measured.
- Run 2 (`f156_t_raw.txt`) measures the five others and the NP crops. The selection is read over both runs' lines.

## Note after T's first resume (NOT blind): the machine's memory

- The first resume (`f156_t_raw_run2_oom.txt`) also failed at « TA ×1.4 p 0.5 », its first measured world. The commit
  charge was full (11 GB free of 94 GB): the IDE's clangd (PID 39916, 7.8 GB) and two rust-analyzers (5.7 and
  4.5 GB) held it.
- At the author's choice, I stopped clangd 39916 and the older rust-analyzer 24132 (18 GB free afterwards), then
  resumed. Nothing in the bench changed.

