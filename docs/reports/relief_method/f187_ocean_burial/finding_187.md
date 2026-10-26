# Finding 187 — the ocean by connectivity and the burial bounded by the 256² surface: the panel's worlds carry unconnected below-sea depressions on almost every step; the flag changes the world without the deposition little; the bounded burial removes F186's 1.8–3.0 km bumps (now 52–128 m, no filled cell > 50 m above the 256² surface); family (iii) does not move with A

**Status: Parts 0, A, B and C done, on mains. Nothing committed.**
- Declared: `f187_declared.md`. Predictions: `f187_predictions.md`.
- Bench output: `f187_part0_output.txt`, `f187_tier1_output.txt`, `f187_tier2_output.txt`, `checks_f187.txt`.

## A correction first

The F186 commit is **5e40634**, not 6a1d0f2. I gave the wrong hash in my message after the commit, and the F187 brief took it
up. The declaration and the bench name 5e40634.

## Part 0 — the inventory, and the depressions without the deposition

### The definitions (read at 5e40634, `f187_declared.md`)

- **The physics level** tests the sea per cell, with no connectivity:
  - ≤ for the incision's base nodes (`stream_power.rs:670`), the routing (`flow.rs:198`), the MFD area (`flow.rs:902`) and
    the fluvial pass (`deposition.rs:537`, `:595`, `:628`);
  - > for the land mask (`cascade.rs:642`), so the hillslope holds every `!land` cell (`cascade.rs:841-846`);
  - < for the marine domain (`deposition.rs:186`) and the walk (`:145`).
- **The chain** takes the largest 4-connected periodic component of z ≤ 0 (`amplify.rs:491-526`). It puts the ocean at 0 m
  (`amplify.rs:885-890`) and lifts every other cell ≤ 0.5 m to 0.5 m (`amplify.rs:903` → `:442`).
- **A1 differs from the chain only by the strict comparison.** A4 aligns the chain behind its own flag,
  `AmpConfig::sea_strict`, which tier 2 turns on.

### The depressions (the panel, intermediate, no deposition, flag off)

| seed | at the end: count, area, largest, depth | steps with ≥ 1 (of 300) | first / first ≥ 50 cells | the maxima: count, area, largest, depth |
|---|---|---|---|---|
| 8 | 10, 598 km², 83 cells, 1 018 m | 285 | 14 / 28 | 11, 627 km², 137 cells, 1 675 m |
| 18 | 4, 688 km², 181 cells, 1 096 m | 293 | 8 / 14 | 17, 1 636 km², 344 cells, 1 601 m |
| 19 | 4, 1 619 km², 271 cells, 3 069 m | 291 | 2 / 13 | 20, 3 125 km², 562 cells, 3 069 m |
| 23 | 7, 2 898 km², 421 cells, 4 045 m | 292 | 4 / 12 | 23, 7 239 km², 598 cells, 4 243 m |
| 29 | 14, 2 566 km², 252 cells, 3 691 m | 294 | 1 / 13 | 17, 5 999 km², 985 cells, 4 385 m |

- **X1 held, 5 / 5.**
- Every world carries them on 95–98 % of its steps. The large ones (≥ 50 cells) appear from steps 12–28, early in the history,
  not late as I predicted (Y1). They reach 1.6–4.4 km below the sea.
- Without the flag the physics level treats them as the sea and as base levels; the chain lifts them to land. This is the
  incoherence the reviewer read.

## Part A — the ocean by connectivity (`CascadeConfig::sea: Sea { connected }`, off by default)

- **A1**: `amplify::ocean_connected` / `largest_component`, the largest 4-connected periodic component of h < h_sea (the ties
  as in the chain). It is re-read after each uplift (each fluvial sub-step included) and at each marine sub-step. The labelling
  is redone only when a cell crossed the sea, which gives exactly the labelling of every sub-step.
- **A2**: the ocean mask is passed as an argument:
  - to the incision (`incise_with_floor_torus_ocean`): the base nodes, the routing (`compute_flow_ocean`), the MFD area
    (`mfd_accumulation_ocean`), the banks;
  - to the fluvial pass (`pass_routed`): the flood's seeds, the fixed receivers, the lakes, the outlets;
  - to the live marine step per cell (`marine_step_live_sea`): the domain, and the walk fills only the ocean;
  - the land mask is the complement of the ocean, so the hillslope holds only the ocean's cells.
- **A3**: the other below-sea cells are closed depressions. The routing fills them to their sill, the pass fills them as lakes,
  and they receive their catchment's flux.
- **Cost**: `StreamPowerConfig` gets no field, so the cache digests do not move.
- **Declared at the controls**: with the marine step, a cell filled to the sea exactly is the sea per cell (≤) and land for A1
  (<). So the « everything connected » control runs without the marine step.
- **A5's controls**:
  - **the flag off = F186 bit for bit**. A worktree at 5e40634 and the current tree give the same hashes for seed 23: without
    the deposition 6b3a0aa02911cd27, with F186's configuration A fb9a0ba78e86a1bd, its cover 4c13633b7778294b.
  - **a connected world** (a full-width sea strip on the torus) gives the same field with and without the flag, bit for bit:
    without the deposition, and with the lakes and K_eff (`the_f187_controls`);
  - **an interior basin** with the deposition and the live marine step: both mass controls < 10⁻⁹; the basin takes no marine
    deposit and is filled as a lake. Without the flag the marine step fills it.
  - **production's C1 hash, the guard, lib and viz**: see the checks.

## Part B — the burial bounded by the 256² surface (`Transmit::bury_bounded`)

- On a covered cell the volume is max(0, z_256 − b̄) × n, b̄ the mean of the fine basement the chain left. L follows as at F186.
- **Since Σ max(0, L − b) ≥ n (L − b̄), L ≤ z_256 on each cell.**
- **The controls** (`the_burial_closure`):
  - a chain that leaves the coarse basement (b̄ = z_256 − H) gives F186's burial exactly;
  - H = 0 changes nothing;
  - a basement lifted to 0.5 m under a 3 km cover: F186 fills to > 3 000 m, F187 to the surface, 56 m.

## Part C — the measures

### Tier 1: the panel at 256², intermediate (off / off + A / deposition + A, the marine step per cell)

| seed | the sea | circumnavigable | (v) | W₁₀ | plains | mountain | median | large lakes' area | depressions at the end |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 49.3 / 49.3 / 51.4 % | yes / yes / yes | 1 / 1 / 1 | 0.150 / 0.148 / 0.149 | 23.6 / 24.6 / 23.8 % | 47.7 / 47.0 / 46.2 % | 377 / 387 / 502 m | 906 / 1 099 / 2 249 km² | 10 / 8 / 0 |
| 18 | 51.3 / 51.3 / 53.6 % | yes / yes / yes | 1 / 1 / 1 | 0.144 / 0.143 / 0.173 | 58.2 / 58.2 / 54.7 % | 20.2 / 20.1 / 23.7 % | 62 / 64 / 67 m | 410 / 164 / 2 126 km² | 4 / 4 / 0 |
| 19 | 47.3 / 47.4 / 52.2 % | yes / yes / yes | 1 / 1 / 1 | 0.033 / 0.033 / 0.036 | 63.6 / 63.9 / 55.2 % | 12.9 / 13.0 / 19.3 % | 42 / 42 / 54 m | 654 / 659 / 469 km² | 4 / 4 / 4 |
| 23 | 38.1 / 38.3 / 36.6 % | no / no / no | 1 / 1 / 1 | 0.210 / 0.211 / 0.159 | 60.3 / 61.0 / 62.9 % | 18.7 / 18.7 / 21.0 % | 79 / 85 / 69 m | 1 033 / 947 / 1 470 km² | 7 / 7 / 1 |
| 29 | 44.3 / 44.4 / 46.0 % | yes / yes / yes | 1 / 1 / 1 | 0.193 / 0.193 / 0.179 | 43.7 / 45.0 / 46.0 % | 26.7 / 26.5 / 29.0 % | 180 / 184 / 144 m | 735 / 886 / 857 km² | 14 / 14 / 0 |

- **No stop fires**:
  - the mass controls' worst is 7.5 × 10⁻¹¹;
  - no solver;
  - the circumnavigation is lost to the deposition on no seed (depA against offA);
  - A moves the plains by at most +1.3 points (the stop is 5).
- **A without the deposition changes little: X4 held, 5 / 5** (the median within +10 m, the plains within +1.3 points).
  - Nothing fills the depressions without the deposition, so they stay: at the end 4–14 per world, the same as without the
    flag.
- **With the deposition and A the depressions are filled**: at the end 0 on 3 seeds, 1 cell on seed 23, 4 lakes on seed 19.
  - Of the cells that were in a depression at some step, 79–96 % end as land and the rest join the ocean.
  - F186's deposition without A also fills them, as the sea's (seed 23: 3 left, 7 km²). A changes how they are filled (as lakes,
    not by the marine step), not whether.
- **The cost of A**: the physics level's wall median is 89.8 s with the deposition and A, against ~84 s for F186 per cell
  (seed 23: 92.8 against 83.3 s). That is +6 to +10 s. Without the deposition A adds < 1 s.
- **E3**: 2 / 5.
- **The shelves are back** with the per-cell marine step (passive medians 4.7 / 6.0 / 8.2 / 3.1 km on seeds 18, 19, 23, 29,
  against 2.2 km without the deposition).

**Family (iii)**, seed 23: the rock (z − H) at F186's cells:

| cell | off | off + A | deposition | deposition + A | the gap without A → with A |
|---|---|---|---|---|---|
| 41567 | 1 043 | 1 042 | 1 833 | 1 838 | +790 → +796 m |
| 44269 | 217 | 212 | 838 | 705 | +621 → +493 m |
| 45401 | 653 | 645 | 1 002 | 1 009 | +349 → +365 m |
| 60574 | 31 | 90 | 203 | 162 | +172 → +71 m |
| 39851 | 459 | 459 | 1 496 | 1 473 | +1 036 → +1 013 m |
| 39603 | 730 | 718 | 1 712 | 1 742 | +983 → +1 024 m |
| 41684 | 649 | 647 | 1 537 | 1 540 | +888 → +893 m |
| 42442 | 528 | 524 | 1 630 | 1 653 | +1 101 → +1 129 m |
| 39102 | 676 | 673 | 1 627 | 1 644 | +951 → +972 m |

- **X3 refuted.** The gap halves on 1 cell of 9 (60574), shrinks by ~20 % on one (44269), and does not move on 7.
- The raised rock is not the sea's definition. Both runs fill the basins, by the marine step without A and as lakes with it, so
  A does not test the hypothesis that the base level moved. The cause of (iii) stays open.

### Tier 2: the chain to 2 048², seed 23, with A (the physics level and the chain's ocean) and B

- **The stages**:
  - C1 1.3 s;
  - the physics level with the deposition and A 98.6 s;
  - the amplifications 5.3 / 8.7 / 40.1 s;
  - the burial 0.13 s.
- **X2 held**:
  - of 2 458 545 filled fine cells, 0 are more than 50 m above the bilinear 256² surface; the filled cells' greatest excess is
    40 m;
  - with every cell counted, the basement's own amplification reaches 1 852 m above that surface.
- **F186's spots** (the fine maximum around each 256² cell, buried / the basement chain):
  - the lake crop's oval (cell 4980): 52 / 5 m (F186: 1 761 m);
  - the islet (17360): 72 / 3 m (F186: 2 516 m);
  - the bottom spot (56482): 128 / 4 m (F186: 3 008 m).
- **The filled basins** (52 groups of 256² cells whose basement is below the sea, covered, not the ocean):
  - their fine maximum stays near their 256² surface (e.g. 666 against 713 m, 147 against 142 m);
  - where it exceeds the surface (1 277 against 1 131 m), the excess is the basement's amplification, not the fill.
- **The volume, fine against coarse**: every cell +0.0 % (139 053 against 139 027 km³); the land −6.1 %.
- **At 2 048²**, without the deposition (with A) → with: plains 49.8 → 54.1 %, hills 31.2 → 24.7 %, mountain 19.0 → 21.2 %,
  median 78 → 68 m.
- **The plates** (for the author):
  - `images/t2_map2048_off_A_seed23.png`, the whole map;
  - `images/t2_crops_lake_delta_plateau_off_A_seed23.png`, lake, delta and plateau, without | with.
  - « Without » is the world without the deposition, with A.
  - What they show: no white bump is left; the shelves and the filled bays read as smooth shallows and plains; the delta crop
    keeps the basement's straight coast (F186's diagnosis).

## The predictions

| | the claim | the result |
|---|---|---|
| X1 | every world has a depression at some time | **held** 5 / 5 (on 95–98 % of the steps) |
| X2 | no fill > 50 m above the 256² surface | **held** (0 cells; worst 40 m) |
| X3 | (iii) comes down by half with A | **refuted** (1 cell of 9) |
| X4 | A changes the world without deposition little | **held** 5 / 5 |
| X5 | the circumnavigation and (v) hold on 5 / 5 | **held** ((v) 5 / 5; the circumnavigation on the 4 that have it, none lost) |
| Y1 | most small, the large ones late | **refuted**: the large ones come at steps 12–28 |
| Y2 | X2 holds against the bilinear surface | **held** |
| Y3 | X3 fails | **held** |
| Y4 | X4 holds on 4 / 5 | refuted in the good direction (5 / 5) |
| Y5 | the circumnavigation on the 4 seeds that have it, (v) 5 / 5 | **held** |
| Y6 | the labelling costs 15–40 s per world | **refuted**: +6 to +10 s |
| Y7 | the strict chain changes < 1 m in median | not isolated: tier 2 ran the chain strict only |

## For the author

1. **The adoption's conditions**: the mass, the circumnavigation and (v) hold, and no burial bump is left (X2). What remains is
   your look at tier 2's plates.
2. **Family (iii)** is not the sea's definition. The rock rises 0.3–1.1 km under the deposition, with or without A. Its cause
   is open: the base level moved by the filled basins is not tested by A, since both runs fill them.
3. **The cost**: the per-cell marine step with A takes ~90 s at the physics level; the chain's budget is F188.

## The checks (`checks_f187.txt`)

- **The guard 6 / 6 Match** (field and lakes). No digest moved, because `StreamPowerConfig` got no field.
- **Production's C1 hash** is unchanged.
- **The tests**: 673 lib tests pass (F186's 672 plus `the_f187_controls`; `the_burial_closure` extended); 33 viz tests pass.
- `cargo check --workspace --tests`: clean.

## Family (iii)'s diagnosis (after the report; no correction)

**The setup**:
- seed 23, 256², intermediate, the deposition with A, on mains (636 s);
- each run cuts one element:
  - (a) the lakes' and depressions' fill;
  - (b) the marine step (what reaches the sea is exported);
  - (c) the load;
  - (d) the compaction;
  - (e) the recycling;
  - (f) K_eff, read as no fluvial deposition: K_eff off and G = 0, so the explicit pass does not take its place.
- Every run's mass controls are ≤ 4.4 × 10⁻¹¹.
- Output: `f187_family3_output.txt`.

**The rock (z − H) minus the world without the deposition (with A)**, at F186's 9 cells (m):

| cell | whole | (a) no fill | (b) no marine | (c) no load | (d) no compaction | (e) no recycling | (f) no K_eff |
|---|---|---|---|---|---|---|---|
| 41567 | +796 | +798 | +884 | +934 | +899 | +560 | +22 |
| 44269 | +493 | +246 | +190 | +307 | +731 | +332 | +530 |
| 45401 | +365 | +549 | +61 | +164 | +423 | +229 | +397 |
| 60574 | +71 | −6 | +63 | +202 | −56 | +101 | −13 |
| 39851 | +1 013 | +1 018 | +1 167 | +1 220 | +792 | +904 | −21 |
| 39603 | +1 024 | +1 014 | +1 126 | +1 190 | +863 | +885 | +62 |
| 41684 | +893 | +826 | +897 | +943 | +863 | +612 | +327 |
| 42442 | +1 129 | +1 130 | +1 135 | +1 226 | +1 082 | +931 | +122 |
| 39102 | +972 | +988 | +956 | +1 058 | +1 001 | +782 | +102 |
| **mean** | **+751** | +729 | +720 | +805 | +733 | +593 | **+170** |

**The share of the whole run's mean gap that each cut removes**:

| (a) | (b) | (c) | (d) | (e) | (f) |
|---|---|---|---|---|---|
| 3 % | 4 % | −7 % | 2 % | 21 % | **77 %** |

- **K_eff carries the gap.** Without it the mean falls from +751 to +170 m; on 6 of the 9 cells the rock comes back within
  ±130 m.
  - This is the steepening Yuan's equation 24 predicts (JGR p. 9: K p̃^m A^m S^n = U (1 + G/p̃)). The divided K keeps the
    channels steeper, so the rock stands higher above the same base level.
  - The two cells where K_eff does not carry it (44269, 45401) move with (a) and (b) instead.
- **The recycling carries the next 21 %.** That figure is measured; its mechanism is not.
- **The reviewer's lead is refuted on these cells** (the coast advances, the rivers lengthen, the profile rises upstream):
  - with the whole deposition the coast moves in both directions (land gained 3 445 km², lost 3 870 km²);
  - the distance to the nearest ocean does not change at 8 of the 9 cells; it shortens 4 km at 44269;
  - the rivers to the sea are as long or shorter (−78 to +4 km);
  - without the marine step (b) the coast barely moves (+200 / −398 km²), and the gap stays (+720 m against +751 m).
- **Corrected nothing.** The decision on K_eff's effect goes to the author.

## The author's words, the recycling with the load cut, the provenance of K

**The author's words (F187)**:
- On the marine step: « Oui marin en cellule. »
  - The per-cell marine step is the reference. F188 will try to accelerate it: the exact solution of the diffusion by FFT, then
    the projection of the constraints.
  - If the acceleration fails, the fallback is « no marine step » (the diagnosis' run (b)), not the blocks.
- On F187's plates: « Tel quelle oui les planches semble ok. »
- Nothing is adopted until the author decides on K_eff.

**The recycling with the load cut** (seed 23, 256², intermediate, the deposition with A; the test the reviewer asked for, 88 s
on mains):
- **The reviewer's hypothesis**: the recycling lightens the load offshore, and the flexure lifts the neighbouring coasts.
- **With the load**, the recycling moves the mean rock gap at F186's 9 cells by 158 m (+751 m with it, +593 m without).
- **With the load cut**, it moves it by 14 m (+805 m with it, +819 m without). Per cell the two runs differ by 11–35 m.
- So without the load the recycling no longer changes the rock: the hypothesis' test, as written, holds. Reported, not
  concluded.

**The provenance of the cascade's K** (no bench; read in the code, the declarations and the sources in `docs/refs`):
- **The cascade's K is not relief_v3's.**
  - `CascadeConfig::solver` overwrites it (`cascade.rs:173`, `sp.k = self.k`) with `CascadeConfig::declared`'s
    `k = 2.0e-3` (`cascade.rs:126`), i.e. K_SI = 2 × 10⁻⁶ yr⁻¹ for m = 0.5, n = 1.
  - It was declared at F159 (`f159_declared.md`, D4, lines 99–100: « K is physical and the same at every level »).
  - **No source and no page are cited for this value**, and no measurement derives it. F170's reading calls it « an erodible
    value ».
  - So it is documented as neither (a) nor (b).
- **relief_v3's own K** (`RELIEF_V1_K × RELIEF_V3_K_MULT` = 1 500 × 3, `stream_power.rs:526-530`, `:465`) is unused by the
  cascade. It came from a morphometric target, « half the incision-reproducing K=3000 », ADR §sculpt and Finding 6, with no
  published source.
- **K2, the rock** (F171, `f171_declared.md` lines 57 and 92):
  - ρ_K = 10 and σ_lnK = 0.5 are read from Stock & Montgomery 1999 (JGR 104(B3)), p. 4983 (abstract) and Table 2, p. 4986.
  - Only ratios are used: ρ_K multiplies the unsourced K_couv above.
- **How Stock & Montgomery obtained K** (p. 4983, the abstract): « We address these uncertainties by simulating the lowering
  of ancient river profiles using the finite difference method. We vary m, n, and K to match the evolved profile as closely as
  possible to the corresponding modern river profile over a time period constrained by the age of the mapped paleoprofiles. »
  - Their law is KA^m S^n with no deposition term, fitted on natural rivers (Australia, Kauai, Japan, California).
  - So their K is of kind (a), effective.
  - The ADR noted this in the C-3 lithology section (line 2520): « the detachment-limited production regime is exactly the domain Stock &
    Montgomery calibrated K for ».
- **Yuan et al. 2019** (JGR manuscript, p. 9, eq. 24–25): « S = [(1 + G/p̃) U / K_f]^(1/n) A^(−m/n) p̃^(−m/n) », and « a new
  steepness index can be defined: k′s = (1 + G/p̃)^(1/n) (U/K_f)^(1/n) p̃^(−m/n) ».
  - p. 9: « For typical values of G/p̃ = 1 (Davy and Lague, 2009) and n = 1 (Whipple and Tucker, 1999), the effect of sediment
    deposition is to increase the steepness index by a factor of 2 ».
  - p. 18, §5.1: « the steepness index, usually used to derive information about uplift and landscape properties, is also
    sensitive to the amount of deposition within the fluvial network ».
  - There K_f is the detachment coefficient of their model with explicit deposition (kind b).
- **In short, reported without concluding**:
  - the absolute K the cascade divides at K_eff has no documented provenance;
  - the only sourced part, K2's ratios, comes from a K of kind (a).
  - The decision on K_eff is the author's.

## The adoption's measures, without K_eff

The author's decision on K_eff, verbatim: « oui on supprime k_eff ».

**The measures for the adoption, without K_eff** (`a187`: F186's A with `fluvial_keff` off and G = 0, with A and B; on mains;
outputs `f187_adopt_panel_output.txt`, `f187_adopt_chain_output.txt`):
- **The panel at 256², intermediate** (the world without the deposition, with A → with the deposition):

| seed | the sea | circumnavigable | (v) | W₁₀ (active / passive) | plains | mountain | median | shelf passive / active | large lakes (area) | all lakes |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 49.3 → 49.8 % | yes → yes | 1 | 0.148 → 0.168 (0.351 → 0.196 / 0.123 → 0.164) | 24.6 → 23.3 % | 47.0 → 43.7 % | 387 → 429 m | 2.2 / 2.2 → 3.8 / 6.0 km | 3 → 4 (1 099 → 7 739 km²) | 1 746 → 8 445 km² |
| 18 | 51.3 → 53.6 % | yes → yes | 1 | 0.143 → 0.144 | 58.2 → 53.9 % | 20.1 → 20.6 % | 64 → 74 m | 2.2 / 2.2 → 5.3 / 7.5 km | 1 → 2 (164 → 6 289 km²) | 637 → 6 641 km² |
| 19 | 47.4 → 54.3 % | yes → yes | 1 | 0.033 → 0.031 | 63.9 → 61.9 % | 13.0 → 12.7 % | 42 → 47 m | 2.2 / 2.2 → 6.2 / 56.2 km | 2 → 2 (659 → 332 km²) | 1 064 → 898 km² |
| 23 | 38.3 → 34.9 % | no → no | 1 | 0.211 → 0.142 (0.410 → 0.268 / 0.091 → 0.098) | 61.0 → 63.2 % | 18.7 → 17.1 % | 85 → 60 m | 2.2 / 3.8 → 6.9 / 6.0 km | 3 → 5 (947 → 2 876 km²) | 1 487 → 4 006 km² |
| 29 | 44.4 → 44.8 % | yes → yes | 1 | 0.193 → 0.167 (0.450 → 0.350 / 0.087 → 0.093) | 45.0 → 49.3 % | 26.5 → 24.9 % | 184 → 129 m | 2.2 / 3.1 → 3.1 / 4.7 km | 4 → 3 (886 → 4 568 km²) | 1 921 → 5 679 km² |

- **The conditions** (the panel):
  - the mass controls' worst is 6.7 × 10⁻¹¹;
  - the circumnavigation is lost to the deposition on no seed;
  - (v) holds 5 / 5;
  - the physics level's wall median is 89.5 s.
- **What changes from the run with K_eff**:
  - the lakes grow a lot on 4 seeds (all lakes ×2.7 to ×13; seeds 8, 18 and 29 gain large lakes of 4 600–7 700 km²);
  - unconnected depressions are left at the end on seed 19 (30, 1 208 km²), where K_eff's run filled them.
- **Family (iii)**: this configuration is exactly the diagnosis' run (f) (K_eff off, G = 0, A on). At F186's 9 cells the rock
  minus the world without the deposition is +22, +530, +397, −13, −21, +62, +327, +122, +102 m, mean **+170 m**, against
  +751 m with K_eff.
- **The chain to 2 048², seed 23** (A, B, the chain aligned):
  - the stages: the physics level with the deposition 98.5 s; the amplifications 6.6 / 9.8 / 46.9 s; the burial 0.10 s;
  - **X2 now fails narrowly**: 7 of 2 397 347 filled fine cells are more than 50 m above the bilinear 256² surface (the worst
    110.5 m). With K_eff: 0.
  - F186's spots are at 57 / 47 / 137 m;
  - the volume, fine against coarse: −0.1 % on every cell, −3.5 % on land;
  - at 2 048²: plains 49.8 → 55.4 %, hills 31.2 → 26.9 %, mountain 19.0 → 17.7 %, median 78 → 59 m.
- **The plates** (for the author):
  - `images/ad_map2048_off_A_seed23.png`, the whole map;
  - `images/ad_crops_lake_delta_plateau_off_A_seed23.png`, lake, delta, plateau; without | with; « without » is the world without
    the deposition, with A.
  - What I see, for the author's eye:
    - **small brown specks along the south coast** of the map (not located; they are possibly the 7 cells above 50 m);
    - in the delta crop, **the whole bay is filled to land**, and where the fill meets the mountain flank (top left) its edge
      shows small steps;
    - the lake and plateau crops are smooth.
- **Not done before the author's look**: the adoption into the continent profile, and the commit.
- The checks are to be rerun before the commit; the code changed only by doc comments since F187's checks.

## Before the adoption: the aspect changes and two corrections

**Before the adoption: the aspect changes and two corrections** (on mains; outputs `f187_load_lakes_output.txt`,
`f187_x2_specks_output.txt`, `f187_plates_output.txt`):

- **1. The lakes and seed 19's drowned coast: the load** (the panel without K_eff, the load cut, the rest unchanged). Each cell
  reads the world without the deposition / the full configuration / the load cut:

| seed | the sea | the coast's length ratio | large lakes (area) | all lakes | unconnected depressions | circumnavigable |
|---|---|---|---|---|---|---|
| 8 | 49.3 / 49.8 / 45.2 % | 1.01 / 1.02 / 1.07 | 1 099 / 7 739 / 3 379 km² | 1 746 / 8 445 / 4 556 km² | 8 / 1 / 1 | yes / yes / yes |
| 18 | 51.3 / 53.6 / 50.3 % | 1.20 / 1.11 / 1.16 | 164 / 6 289 / 1 023 km² | 637 / 6 641 / 1 208 km² | 4 / 0 / 0 | yes / yes / yes |
| 19 | 47.4 / 54.3 / 44.7 % | 1.16 / 1.68 / 1.04 | 659 / 332 / 803 km² | 1 064 / 898 / 1 042 km² | 4 / 30 / 3 | yes / yes / yes |
| 23 | 38.3 / 34.9 / **14.8 %** | 1.12 / 1.05 / 1.11 | 947 / 2 876 / 20 959 km² | 1 487 / 4 006 / 22 073 km² | 7 / 5 / 9 | no / no / no |
| 29 | 44.4 / 44.8 / **0.0 %** | 1.23 / 0.93 / — | 886 / 4 568 / 0 km² | 1 921 / 5 679 / 0 km² | 14 / 3 / 0 | yes / yes / **no** |

  - **The load makes the new lakes and the drowned coast.** Without it:
    - the large lakes shrink by 56 % (seed 8) and 84 % (seed 18);
    - seed 19's coast returns: the sea 54.3 → 44.7 %, the length ratio 1.68 → 1.04, 30 → 3 depressions.
  - **But without the load the deposition fills the ocean again** (F182's problem): seed 29's sea goes, and its
    circumnavigation with it; seed 23's falls to 14.8 %.
  - **Where the new large lakes are**: near the coast (their nearest cell 3–13 km from the ocean; one at 70 km), their outlets
    at 0–40 m.
  - **The load's deflection at their outlets is negative everywhere**: −61 to −1 041 m. It is the full run's cumulated rebound
    minus the run without the load, an approximation (the erosion's rebound differs a little between the two).
  - **The reviewer's hypothesis is refuted**: there is no flexural bulge barring the rivers. The outlets are lowered, so the
    lakes are coastal lowlands the load subsided toward the sea level.
  - Maps: `images/ll_lakes_full_seed*.png`, `ll_lakes_noload_seed*.png`, and `ll_deflection_seed*.png` (red up, blue down, the
    lakes in grey).
- **2. The 7 cells above the surface and the south coast's brown specks** (seed 23, 2 048²): **they do not coincide** (the
  specks are 864–1 132 fine cells from the 7).
  - **The 7 cells** are in the deep sea (around −3 050 m, map pixel 192–193, 287), in one covered cell (36704, H 1 198 m, no
    lake), 70–110 m above the bilinear surface.
    - The origin is in `amplify.rs`, `Transmit::bury_with`: on a covered cell whose bounded volume is 0 (z_256 ≤ b̄), L was
      the minimum of the fine basement.
    - That minimum can stand above z_256 where the chain raised the basement. The interpolation carried it into the
      neighbours, so B2's bound held per cell but not across them.
  - **The south coast's specks** are the basement chain's own peaks on uncovered cells (H = 0; buried = the basement chain,
    712–970 m). The world without the deposition has them too (613–878 m); the moved coast now leaves them at the shore.
  - Separately, sediment piles of ~1 km stand within 2 km of the coast on the north centre (z_256 0.65–1.12 km, H 0.67–1.47 km,
    against ~0 m without the deposition). Their origin is not diagnosed; a hypothesis to check: lakes filled to a high spill.
- **3. The corrected burial, `Transmit::bury_edge`** (no parameter):
  - on an uncovered cell, L is the minimum of its fine basement (the author's step 3);
  - on a covered cell with no bounded volume, the same;
  - both capped by z_256 (declared: the cap keeps B2's bound across the cells, and is step 2's correction).
  - **Controls**: H = 0 changes nothing; a covered cell with no room fills nothing.

  | seed | X2 (cells > 50 m), worst | the volume filled outside the covered cells | the total fill against the coarse |
  |---|---|---|---|
  | 8 | 2 → 0 cells (53.1 → 13.6 m) | 242 → 38 km³ | +0.2 → +0.1 % |
  | 18 | 12 → 0 cells (141.6 → 2.7 m) | 181 → 19 km³ | +0.5 → +0.3 % |
  | 19 | 66 → **18** cells (156.3 → 74.4 m) | 100 → 7 km³ | −0.1 → −0.2 % |
  | 23 | 7 → 0 cells (110.5 → 12.6 m) | 263 → 44 km³ | −0.1 → −0.3 % |

  - Seed 19 keeps 18 cells above 50 m (the worst 74 m). Not located this round.
- **4. The plates** (for the author; the chain to 2 048², without K_eff, A, B corrected; « without » is the world without the
  deposition, with A):
  - seeds 8, 18, 19, the whole map without | with: `images/pl_map2048_off_A_seed8.png`, `_seed18.png`, `_seed19.png`;
  - seed 23, the delta crop without | the bounded burial | corrected: `images/pl_crop2048_delta_off_old_new_seed23.png`.
  - **What I see, for the author's eye**:
    - seed 23's delta: the steps at the fill's edge are smoothed (a small notch remains at the top);
    - seed 8: the east coast breaks into fringes and islets;
    - seed 18: a wide shallow shelf with a ragged edge on the west coast;
    - seed 19: the south coast steps back behind a wide shallow sea (the load's drowned coast).
- **Nothing adopted, nothing committed**, before the author's look at these plates.

## The targeted round before the adoption

**The targeted round before the adoption** (on mains; outputs `f187_round2_output.txt`, `f187_deflection_rel_output.txt`; seeds 8,
18, 19, 23 at 2 048², without K_eff, A, B):

- **My error, corrected**: last round I read the load's deflection at the new lakes' outlets as negative everywhere and wrote
  « the reviewer's hypothesis is refuted ». That deflection held the torus' uniform part: the flexure of the mean load at
  wavenumber 0, from −457 to −936 m depending on the seed, which dominated it.
  - **With the uniform part removed** (`ll_deflection_rel_seed*.png`), the land stands raised against the loaded sea, at
    +100 to +570 m around the lakes.
  - But each outlet sits level with the land within 30 km: the outlet against that land's median is −27 to +36 m on 10 of 11
    lakes, and −73 m on one (seed 29). There is no local bulge on the outlets.
  - **The lake floor against its outlet**, the comparison that would decide whether the load digs the lakes, is not measured.
    The reviewer's hypothesis is therefore neither confirmed nor refuted.
- **The reviewer's reading, recorded** (a hypothesis to check): the large coastal lakes come from the subsidence under the
  sediment load, amplified by the absence of fluvial aggradation on the coastal plains (the loss declared at F186, A3). On
  Earth, the rivers aggrade these sinking lowlands.
- **1. The islets at the filled coasts** (the fine land cells whose 256² parent is the sea):

  | seed | islets with `bury_edge` | the basement chain itself | under a covered parent | with `bury_smooth` |
  |---|---|---|---|---|
  | 8 | 11 688 | 87 % | 49 % | 11 744 |
  | 18 | 12 187 | 94 % | 46 % | 12 213 |
  | 19 | 47 721 | 71 % | 75 % | 36 491 (−24 %) |
  | 23 | 11 355 | 91 % | 47 % | 11 294 |

  - **The diagnosis: yes, the islets are the basement's amplification piercing the surface** (71–94 % of them are the
    basement chain itself).
  - **But the correction (`Transmit::bury_smooth`) barely touches them**. Half of them sit under an uncovered 256² sea cell,
    where the correction does not apply. Under the covered ones the cover is thinner than the fine relief, so the condition
    H > max − mean does not fire. The islets change on seed 19 only.
  - What the correction does do:
    - X2 is 0 on the 4 seeds (seed 19: 18 → 0; the worst 74 → 46 m);
    - the net fill against the coarse stays −0.2 to +0.4 %;
    - H = 0 changes nothing (the unit control).
- **2. Seed 18's straight shelf edges are not the 64² recycling boundaries** (diagnosis only):
  - of the −150 m contour's 2 262 map pixels, 1 % lie within 1 pixel of a boundary of the cells C1 ever consumed, against 3 % of
    all map pixels;
  - map `r2_seed18_recycling_boundaries.png`. Their origin is not found.
- **3. The origins** (diagnosis only):
  - **Seed 19's 18 cells above the surface** are all covered lakes at the coast. Each is raised to its parent cell's flat bed:
    `amplify.rs:1530`, `if self.lake[p] && self.cover_m[p] > 0.0 { z = z.max(surface[p]) }`.
    - That reads the parent cell, against the author's rule. Its sea-side neighbours pull the bilinear surface below it by
      50–64 m.
    - `bury_smooth` removes them where the cover exceeds the relief.
  - **Seed 23's piles near the north coast** (cells 23007, 26553, 26560, 27053, 27050, 26800, 26803):
    - z 653–1 122 m against ~0 m without the deposition; H 674–1 470 m;
    - the cover is marine (2.1–2.9 km cumulated), with no lake and no fluvial part;
    - the basement is at −16 to −408 m; the cells are now neither the ocean nor a lake.
    - A marine deposit above the sea can come only from the walk's stranding (`deposition.rs:199-211`: the remainder put on
      the cell where the descent ends), or from a later rise of the cell.
    - **A hypothesis, not instrumented**: the second, through the rock law, which counts the deposit as basement once E_cum
      passes the initial cover (`cascade.rs:684`: K / ρ_K, ten times less erodible). The uplifted marine fill would then
      erode 10 × slower than in the world without the deposition. This is the queued item « sediment vs rock erodibility ».
- **4. For the author**:
  - the lakes: `images/r2_lakes_full_noload_deflection_seeds8_18_19_29.png` (with the load | without | the deflection with its
    uniform part, so read `ll_deflection_rel_seed*.png` for the deflection);
  - the 2 048² maps redone with `bury_smooth`: `images/r2_map2048_off_A_seed{8,18,19,23}.png`.
- **Nothing adopted, nothing committed.**
