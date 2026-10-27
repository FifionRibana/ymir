# Finding 189 — the stranding on dead ends (Part 1, the diagnosis)

**The step commit before F189**: `88c089f` (B adopted into `Deposition::continent_default()`, A and A′ removed; the guard
6 / 6 Match without re-indexing, production's C1 hash unchanged, 674 lib and 33 viz tests pass, checked in a worktree at
HEAD with only these files, the author's `config.rs` and `workspace.rs` left out).

**Part 1, the diagnosis** (one instrumented run; seed 8, intermediate, the committed configuration: `continent_default()`,
`Sea { connected: true }`, C1's consumption; on mains; output `f189_part1_output.txt`, 138 s):
- **The instrument**: observers only (`deposition::StrandObs`, `physics_history_strand`): the fluvial pass at each sub-step
  (the surface it routed, the sea mask, the pass), each stranding of the marine walk (`MarineStrand`). The control: the field
  with and without them is the same bit for bit (68 s with, 67 s without).
- **Where the masks are built**:
  - the sea mask, both for the routing and for the lakes: `cascade.rs:609` `land_conn` → `amplify::ocean_connected`
    (`amplify.rs:540`) → `largest_component` (`amplify.rs:547`), **4-connected** (`amplify.rs:563`); the fluvial stage's
    `sea_mask` is `!land` (`cascade.rs:760`);
  - the routing: `deposition::pass_routed` → `flow::compute_flow_ocean` (`flow.rs:204`). The flood `pit_fill`
    (`flow.rs:250`) seeds the **8**-neighbours of the ocean at their own height (`flow.rs:256`). D8 (`compute_d8`) takes a
    strictly lower neighbour (`flow.rs:331`), else an equal flat neighbour with a strictly smaller Garbrecht-Martz gradient
    (`flow.rs:352`). `flat_resolution` counts any ocean neighbour as lower (`flow.rs:628`);
  - the lakes: `is_lake` (`deposition.rs:723`) reads the same `flow.filled` as the routing. **The routing and the lake mask
    cannot diverge**: they come from one flood of one surface. The divergence is elsewhere: the ocean mask is 4-connected, but
    the flood and D8 are 8-connected and take the ocean as lower than any land.
- **Fluvial** (`deposition.rs:761`, `p.stranded[k] += tot`): 154 strandings on 95 sub-steps, **3 956.64 km³ of solid**. The
  record's `fluvial_m` (G = 0) gives the same total.
  - **(iii) the ocean mask against the flood: 114 strandings, 1 527.66 km³ (38.6 %)**:
    - every one is below the sea, out of the ocean, and touches the ocean only diagonally;
    - such a cell is a below-sea cell that the 4-connected labelling leaves out of the ocean, and that the 8-connected flood
      seeds at its own height. When it stands lower than every ocean neighbour, it has no lower neighbour, no lake (filled =
      surface) and no receiver;
    - the dead ends have 1 cell at the median (at most 17).
  - **(ii) a flat whose exit Garbrecht-Martz leaves unrouted: 40 strandings, 2 428.98 km³ (61.4 %)**:
    - every one is a flat cell (`needs`) of a filled depression with flat_grad 0, whose only exit is an equal-height land cell
      outside the flat, also at flat_grad 0. D8's flat pass asks for a strictly smaller gradient (`flow.rs:352`), so the cell
      has no receiver (40 / 40; 0 next to a higher ocean, 0 other);
    - the gradient is tl·(fh_max + 1) + (fh_max − fh). A spill cell (tl = 0) gets 0 when fh = fh_max, the grid's maximum, or
      when fh is unreached. A hypothesis, not checked: on these sub-steps fh_max = 0 (only flats one cell wide), so every
      spill ties;
    - 11 of the 40 are below the sea, 7 touch the ocean diagonally; the dead ends have 3 cells at the median (at most 27).
  - **(i) a depression missing from the lake mask: 0.** (iv) other: 0.
- **Marine** (the walk, `deposition.rs:226–242`): 211 958 strandings, **13 020.51 km³ in thickness** (4 817.59 km³ of solid
  at φ 0.63; the record's `marine_stranded_m` gives the same thickness).
  - **The cause**: the walk goes down the raw surface to a strictly lower neighbour (`deposition.rs:242`) and fills nothing on
    its way. **Every closed depression and every flat strands it.**
  - (a) below the sea, out of the ocean, touching it only diagonally: 7 577; 250.35 km³ (1.9 %); in pits.
  - (b) below the sea, out of the ocean, not touching it (a closed basin): 49 510; 2 124.33 km³ (16.3 %); 2 100.33 km³ in
    pits.
  - (c) at or above the sea: 154 185; **10 600.68 km³ (81.4 %)**; the median h − sea is 0.3 m. Of these:
    - 6 399.87 km³ in pits (every neighbour higher): land depressions the walk does not fill;
    - 4 200.81 km³ on flats, of which **2 439.98 km³ at the sea's level exactly**. These are the ocean cells the walk itself
      filled to the sea: at h = h_sea they leave the domain (h < h_sea), and they form a flat at the sea's level with no
      descent.
  - (d) in the ocean, filled to the sea, no lower neighbour: 686; 45.15 km³ (0.3 %).
  - (e) the walk's step limit: 0.
- **The 10 largest** (cumulated per cell): all fluvial. In order: (x 169, y 108) 1 479.85 km³ (ii); (138, 104) 1 319.90 km³
  (iii); (159, 114) 208.80 (ii); (173, 201) 197.45 (ii); (166, 119) 123.75 (ii); (168, 98) 104.40 (ii); (178, 96) 77.16 (iii);
  (146, 108) 48.41 (ii); (165, 191) 39.10 (iii); (166, 189) 38.58 (iii). The two largest make 71 % of the fluvial volume, in
  2 strandings each.
- **What Part 1 says for Part 2** (nothing corrected):
  - the fluvial causes are two inconsistencies of the routing, (iii) and (ii), not a missing lake. By the brief's order, the
    correction is to fix them, not the generic « every dead end a one-cell lake ».
    - (iii): the ocean mask and the flood must use one connectivity, or the flood must seed at the ocean's heights.
    - (ii): D8's flat pass must accept the equal non-flat exit.
  - **Both sit in `flow.rs`, shared with production** (`compute_flow` calls `compute_flow_ocean` with the per-cell mask). A
    correction there moves the guard unless it is gated to the deposition's routing. That decision is the reviewer's and the
    author's.
  - The marine stranding is 1.2 × the fluvial in solid (3.3 × in thickness). Its cause is the walk itself, which routes on the raw surface.
    Applying the principle (« toute cuvette fermée … est comblée par volume jusqu'à son déversoir ») means a walk on a flooded
    surface, plus a way out of the sea-level flats it creates.
- **F189 stops here**, before Part 2, for the reviewer's look.

**Part 2: the brief and the predictions** (written before any computation of Part 2):
- **The author's decision** (relayed by the reviewer): the corrections of `flow.rs` apply everywhere, production included;
  the guard will move.
- **The reviewer's predictions**:
  - **V2 bis redone** (V2's setup, the lakes on, no fluvial deposit; fluvial steps 160 / 40 / 10 kyr):
    - no land cell above 5 × the land's median, at any step;
    - on the quantities above E3's noise floor, r = change(40, 10) / change(160, 40) between 0.1 and 0.5, with the same sign.
  - **If V2 bis passes, E3's triplet 20 / 10 / 5 kyr redone** on the corrected adopted configuration:
    r = change(10, 5) / change(20, 10) between 0.25 and 1, with the same sign, on the quantities above the floor. Production
    is added against the 5 kyr step.
  - **The reading**:
    - if both pass, the step dependence came from the stranding. Production's change against the fine step is documented
      as an accepted limit, and the hillslope's splitting in sub-steps is not needed;
    - if V2 bis fails, stop: the hillslope's splitting is decided next.
- **CC's predictions**:
  - the controls:
    - the stranded volume (fluvial and marine) is 0 by construction;
    - no land cell without a receiver outside the lakes at any sub-step on seed 8;
    - the mass < 10⁻⁹ with the « exported » term;
  - **the « exported » volume** (the ocean filled with no deeper cell): below 1 % of the marine sources on every panel
    seed. A hypothesis: the open water of the panel's worlds stays deep enough at the intermediate age;
  - **the guard**: the six states' field and lake digests all move. `run_hd`'s routing goes through `compute_flow`; the
    per-cell sea makes (iii) impossible there, but the flats' tie (ii) touches every flat;
  - **V2 bis**: a hypothesis. G1 bis traced V2's spikes to the fluvial stranding (the talus shed 93 km in one step), so the
    spikes go and the 5 × criterion passes; r is not predicted beyond the reviewer's range;
  - **the step dependence (E3)**: a hypothesis, no number. Most of the fluvial stranding sat on two cells (71 %), so part of
    E3's dependence may go with it.

**Part 2.0, the readings** (one instrumented run each; output `results/p20.txt`):
- **(1) Where the marine walk starts.** The sources are built at `cascade.rs` (`let mut src`, the marine step):
  - on a sea cell: the rivers' outflow (`ps.to_sea`) plus the hillslope's net loss of the land, spread over the coastal sea
    cells by `deposition::coast_weights` (`deposition.rs`, a sea cell weighted by the height its land 4-neighbours stand above
    it);
  - on a land cell, only with the live marine step: a river's outflow on a cell that a later fluvial sub-step raised (k > 1).
  - **Measured on seed 8, the committed configuration**: 190 050 km³ of sources, **100.0 % in the ocean** (the rivers
    182 541 km³, the hillslope 7 510 km³). None on land, none on a sea cell raised above the sea before the marine step, none
    below the sea out of the ocean. The adopted profile runs k = 1.
  - So Part 1's class (c) (81 % of the marine stranding) did not start on land: a walk filled its ocean cells to the sea,
    they left the domain (h = h_sea is not h < h_sea), and the next parts stranded on that flat at the sea's level. Each
    stranding then raised the cell, and the next parts piled up on it (the 6 400 km³ of « pits » up to 1 910 m: 17 538
    strandings on their own source cell).
- **(2) The incision's routing** (`incise_with_floor_torus_ocean` routes with `compute_flow(_ocean)`, `stream_power.rs:672`):
  - a land cell without a receiver is a **fixed base node** (`receiver[k] == k`, skipped, `stream_power.rs:758`): it neither
    erodes nor rises, and its catchment incises toward it as toward the sea;
  - **the committed configuration** (seed 8): 87 of the 300 routed fields hold such cells, 112 cell-fields in all, every one
    of class (iii), none in a lake;
  - **the production cascade** (no deposition, the per-cell sea): 14 of 300 fields, 213 cell-fields, every one the flats'
    tie (ii), 211 of them in a filled depression. Class (iii) cannot occur there (every below-sea cell is the sea).

**Part 2.2, the hypothesis checked** (output `results/p22check.txt`): at the 42 fluvial strandings of the tie on seed 8,
fh_max = 0 on 41 (every flat cell next to higher terrain); the 42nd has fh = fh_max = 1, the formula's other branch, also 0.
The hypothesis holds. (Part 1 counted 40: two strandings whose flat touched the ocean were put in class (iii) there; here the
class is read per cell.)

**The reviewer's citation (recorded)**: « Barnes, Lehman & Mulla 2014 (Computers & Geosciences 62, 117–127) » is their
Priority-Flood paper (doi 10.1016/j.cageo.2013.04.024). The flats paper is « An efficient assignment of drainage direction
over flat surfaces in raster digital elevation models », Computers & Geosciences 62 (2014) 128–135, doi
10.1016/j.cageo.2013.01.009. Read in its open version, arXiv:1511.04433; the pages below are the arXiv version's, not the
journal's (the PDF is not deposited in `docs/refs`).

**Part 2.1 and 2.2, the routing's corrections** (`terrain/flow.rs`, everywhere, production included):
- **2.1, one connectivity for every link with the ocean** (the ocean stays the largest 4-connected component):
  - the flood `pit_fill` seeds only the ocean's 4-neighbours;
  - D8 (`compute_d8`, and `compute_dinf`'s slope pass) takes an ocean receiver only through a side, and a side ocean
    neighbour always counts as lower (its slope floored at the smallest positive value);
  - `flat_distances` (the classification of `flat_resolution`) counts an ocean neighbour only through a side; a corner ocean
    neighbour is no neighbour;
  - between two land cells D8 stays 8-connected.
  - Test `the_ocean_is_linked_through_sides_only`: a below-sea pocket touching the ocean by a corner, deeper than the ocean
    there, becomes a lake filled to its sill, has a receiver, and every land cell reaches the ocean. It fails on the old code.
- **2.2, the flats' exit** (Barnes et al. 2014):
  - the towards-lower distance starts at 1 on the flat cells next to the exit (§2.4, arXiv p. 6: « An iteration counter I is
    set to 1 »), so every flat cell's gradient stays above its exit's 0;
  - D8's flat pass (two passes now: the slopes first, then the flats) takes a flat neighbour with a smaller gradient, or an
    equal neighbour outside the flat that already has a receiver from the slope pass. That is Barnes's « low edge » (§2.2,
    arXiv p. 3: « c has a defined flow direction, … at least one neighbor n at the same elevation, and … this n does not have
    a defined flow direction »);
  - no cycle: along a route the flat cells' gradient strictly decreases, and a low edge leaves by a strict descent;
  - the weights are unchanged (tl × (fh_max + 1) + (fh_max − fh)), not Barnes's 2 L + (H_flat − H): they also make the
    towards-lower term dominate, and they keep #drainage fix A's noise bound;
  - test `a_one_cell_flat_drains_to_its_low_edge`: a flat one cell wide drains to its sill and every land cell reaches the
    ocean. It fails on the old code.
- `resolve_flats` went (`flat_resolution` returns the flat cells D8 now reads). The cached steps that route with `flow.rs`
  are bumped: `ALGO_UPSCALE_EROSION` 4 → 5, `ALGO_DRAINAGE` 9 → 10, `ALGO_BREACH` 1 → 2, `ALGO_HD_DRAINAGE` 12 → 13.

**Part 2.3, the marine walk** (`deposition::marine_step_live_sea` with `exported_out`; the cascade passes it):
- **The walk circulates only in the sea, on a filled surface**:
  - it moves on the traversable sea: the 4-connected component of h ≤ h_sea holding the ocean, the cells filled to the sea
    included. It is fixed for the call, since the marine step only raises cells and never above the sea;
  - on an ocean cell it deposits at most the room h_sea − h (the cell set to h_sea exactly when filled), and the rest goes
    on down the steepest descent among the traversable cells;
  - with no lower traversable neighbour, it goes to the nearest ocean cell with room through the traversable cells
    (breadth-first, 4-connected): the progradation over the filled cells;
  - with no such cell, the rest leaves the map (the open window): `SedRecord::marine_exported`, the budget's « exported »
    term, in both mass controls.
- **The sediment that arrives on land** (a source out of the traversable sea) is routed by a fluvial pass on the current
  surface (`pass_routed`: the lakes filled by volume, then the sea) before the walk: `SedRecord::marine_rerouted`. In the
  adopted profile there is none (2.0 (1)).
- The closed submarine depressions of Part 1's classes (a) and (b) are out of the ocean, so the walk never enters them. They
  are filled by the fluvial pass, as lakes up to their sill (F187).
- No new parameter. A source on a cell out of the traversable sea that reaches the walk anyway is stranded in place and
  observed (none measured).

**Part 2.4, the controls** (the mass, the strandings, the receivers):
- **Seed 8, the corrected adopted configuration** (`f189_seed8_controls`, 94 s; `results/p24_seed8.txt`):
  - land cells without a receiver: 0 on the 300 fields the pass routed and 0 on the 300 the incision routed (in lakes or
    out);
  - strandings: fluvial 0, marine 0;
  - the mass: 2.0 × 10⁻¹⁴ per step, 3.3 × 10⁻¹¹ on the stock;
  - exported: 0.000 km³ (0.0000 % of 173 843 km³ of marine sources); land sources rerouted: 0.
- **The permanent test `the_f189_no_stranding`**: on F187's control world with its interior basin below the sea, the adopted
  profile strands nothing, the mass holds at 10⁻⁹, nothing is rerouted, and the walk deposits.
  - Its first version used F187's narrow sea (8 rows of 64) with the rebound the load requires. The sea filled and rose
    above the sea level by step 20, so the flood had no seed and every pit kept no receiver (1 108 norm × cells
    stranded). **A world with no ocean left has no base level**, a case the correction does not cover (declared). The test
    uses a deep, wide sea (24 rows at 3 km).
- **The tests that changed, and why**:
  - `the_steady_physics_level_is_pinned` and `the_f159_extraction_changes_no_output` pin a field by hash: re-pinned, since
    the routing moves the field (the old values kept in the comments);
  - `arid_sink_enclosed_below_sea_stays_endorheic_below_col` (Finding 39's drainage test). With the old code, the flat land
    at 0.51 kept whole rows without a receiver (the tie), so their runoff never reached the pit, and the endorheic verdict
    leaned on that lost water. Routed, the whole flat drains to the pit and outruns the 9-cell floor's evaporation (exorheic,
    level 29 m). The cool catchment is now the pit's one-cell rim (16 cells), the rest hot, so the inflow stays below the
    evaporation by construction, which is the test's intent;
  - `f123_guard`'s control `F109_DELIVERED` re-pinned (0x316ca8f727ca6f62 → 0xa3bdb6b4e98810c0);
  - the F186–F188 controls and `f188_profile_is_b` pass unchanged.
- **The guard, old reference → new, per state** (the viz's eroded fields and final lakes from the cache, `results/
  guard_diff.txt`; the cold run 1 525 s):

| state | cells changed | |Δz| median / p99 / max | lakes | lake area | cells in / out of a lake changed |
|---|---|---|---|---|---|
| livré | 77.2 % | 0.00 / 0.0 / 50.7 m | 58 → 58 | 6 660 → 6 660 km² | 70 |
| A1+B2 | 80.6 % | 0.01 / 0.1 / 131.8 m | 36 → 36 | 5 217 → 5 217 km² | 259 |
| C1 nue | 1.6 % | 0.07 / 0.5 / 44.5 m | 25 → 25 | 3 330 → 3 330 km² | 9 |
| C2 /10 col (défaut) | 80.0 % | 0.00 / 0.1 / 131.9 m | 25 → 25 | 3 492 → 3 492 km² | 64 |
| C2 /3 | 80.6 % | 0.00 / 0.2 / 173.1 m | 25 → 25 | 3 683 → 3 683 km² | 114 |
| C2 /10 niveau mer | 80.5 % | 0.00 / 0.1 / 102.2 m | 23 → 23 | 4 314 → 4 314 km² | 202 |

  - Most cells move by less than a centimetre, a few by tens of metres; the lakes keep their number and area. (The lake count
    is the final drainage bundle's; the lake reference file counts after the crater pass.)

**The guard re-indexed** (in a worktree at 88c089f with only F189's files; outputs in the commit's checks):
- `f123_guard` rewrote `bench_field_hashes.json` (952 s, its own cache). Its six field hashes equal the viz's computed ones,
  so the identity between bench and viz holds.
- `f134_lake_guard` rewrote `bench_lake_hashes.json` (1 823 s). The lake counts after the crater pass are unchanged:
  60 / 38 / 27 / 26 / 26 / 24.
- **The viz guard after: 6 / 6 Match, field and lakes** (133 s warm).

**Part 2.5, V2 bis redone** (the same setup as F188's; 667 s; output `f189_v2bis_output.txt`): **it passes.**
- **The land's max over its median, worst over the steps**: 2.92 / 2.80 / 2.82 at 160 / 40 / 10 kyr, against
  48.2 / 22.3 / 7.7 in F188. The 5 × criterion holds at every step.
- **The quantities**: between 40 and 10 kyr every one is under E3's noise floor, so r is not read:
  - density −1.54 % → −0.99 %;
  - λ +5.33 % → +0.31 %;
  - S₁₀ −0.53 % → −0.18 %;
  - the classes' largest move 0.38 point.
- **The predictions**:
  - the reviewer's spike criterion holds;
  - the convergence holds. The quantities are all under the floor at the fine pair, so the range 0.1–0.5 does not apply,
    the rule written before;
  - CC's hypothesis (the spikes were the fluvial stranding) holds.
- **The reading**: V2's step dependence came from the stranding. E3's triplet follows.

**Part 2.5, E3's triplet redone** (seed 8, intermediate, the corrected adopted profile; production, 20, 10 and 5 kyr: 83 / 292
/ 489 / 855 s; 1 723 s in all; output `f189_e3_output.txt`): **it passes on all the land.**
- Every run strands nothing; the mass holds at ≤ 2.0 × 10⁻¹⁴ per step and ≤ 3.8 × 10⁻¹¹ on the stock.
- **All the land** (28 979 cells): between 10 and 5 kyr every quantity is under the noise floor:

| quantity | 20 vs 10 | 10 vs 5 | production vs 5 kyr |
|---|---|---|---|
| density | −1.48 % | +0.57 % | +1.22 % |
| λ | −1.62 % | +0.03 % | +0.49 % |
| S₁₀ | +0.50 % | −0.03 % | +0.10 % |
| plain | +0.38 pt | +0.03 pt | −1.10 pt |
| plateau | −0.28 pt | +0.33 pt | +3.40 pt |
| hill | −0.41 pt | −0.41 pt | −2.31 pt |
| mountain | +0.31 pt | +0.06 pt | +0.01 pt |

- **For reading only, the sub-zones by Courant number**:
  - C > 1 has 4 494 cells: density r +0.41 passes; plateau and mountain fail (r −0.62, −55.9);
  - C ≤ 1: λ fails (r −0.57);
  - the other quantities are under the floor.
  - These zones are small and the rule written before reads all the land. On C > 1, production's density and λ against 5 kyr
    (+44 %, +100 %) are computed on few channels and are not read.
- **The reading, by the rule written before**:
  - V2 bis and E3 both pass: **the step dependence came from the stranding**;
  - **production's change against the fine step is documented as an accepted limit**: +1.2 % density, +0.5 % λ, +0.1 % S₁₀;
    on the classes −1.1 / +3.4 / −2.3 / +0.0 points (plain / plateau / hill / mountain). The plateau's +3.4 points is the
    largest;
  - F188 measured +4 % density, +6 % λ and 2–3 points on the classes;
  - **the hillslope's splitting in sub-steps is not needed.**

**Part 2.5, the registered plates and the per-seed numbers** (seeds 8, 18, 19 and 23, intermediate; `f189_plates`, 650 s;
88c089f's numbers by the same function on a worktree at 88c089f, `results/before.txt`; output `f189_plates_output.txt`):
- **The plates**: `images/plates_88c089f_F189_seed{8,18,19,23}.png`. The left panel is F188's registered F187 + B panel (the
  code 88c089f committed); the right one is F189's. The amended registration check passes on every F189 panel: 100.00 % of the
  river pixels on land; mouths 38 / 38, 58 / 58, 58 / 58, 56 / 56.
- **88c089f → F189**:

| seed | all the lakes | stranded fluvial | stranded marine (solid) | exported | circumnavigable | (v) | the sea | the physics level |
|---|---|---|---|---|---|---|---|---|
| 8 | 3 147 → 1 455 km² | 3 956.64 → 0 km³ | 4 817.59 → 0 km³ | — → 0.000 km³ | yes → yes | yes → yes | 50.0 → 55.5 % | 81 → 86 s |
| 18 | 1 406 → 190 km² | 2 731.41 → 0 km³ | 3 893.26 → 0 km³ | — → 0.000 km³ | yes → yes | yes → yes | 54.4 → 55.2 % | 99 → 112 s |
| 19 | 842 → 303 km² | 2 049.46 → 0 km³ | 2 866.77 → 0 km³ | — → 0.000 km³ | yes → yes | yes → yes | 55.1 → 57.5 % | 91 → 111 s |
| 23 | 6 565 → 1 814 km² | 3 706.43 → 0 km³ | 8 366.69 → 0 km³ | — → 0.000 km³ | no → no | yes → yes | 34.6 → 37.1 % | 108 → 127 s |

- **What the plates show, for the author's eye**:
  - the flat filled surfaces with their 45° lines and solid patches recede, but not everywhere: they leave seed 18's
    north-east, and some remain on seed 23's south-east;
  - the islet mosaic of seed 19's south bay goes;
  - seed 8's east-coast steps (`bury_smooth`'s blocks) go;
  - the circled lakes of seed 23 (the raised sill, the merge) go;
  - **in their place come wide light-blue platforms at the sea's level**, notably seed 8's east coast and seed 19's south bay.
    These are the cells the walk now fills up to the sea before it progrades. On seed 8, 2 013 cells (4 914 km²) stand less
    than 1 m below the sea and 753 cells less than 1 m above it (`results/e3_prod.bin`). The platform does not stay at 0 m
    exactly: the re-compaction, the diffusion and the flexure move it.
  - A hypothesis, not checked: most of 88c089f's lakes and islets were dams and pits built by the stranded piles (the marine
    piles up to 1.9 km on the coast). The lake area falls by 54–86 %, and the sea grows by 0.8–5.5 points.
- **The cost**: the physics level takes +5 to +20 s per world (+6 to +22 %). That is the progradation's breadth-first search
  and the slope pass now separated from the flats' pass. **A red flag by the project's rule on generation time**; not
  optimised in this finding.
- **The checks of the code as it stands** (a worktree at 88c089f with F189's files): 677 lib and 33 viz tests pass;
  production's C1 hash unchanged; the guard 6 / 6 Match after re-indexing.
- **F189 stops here, nothing committed**, for the author's look at the plates. The order stays: F190 the torus' unfolding,
  F191 the hillslopes' unloading, F192 the fill by compaction with `bury_smooth`'s parent read removed, then the budget.

**F189's closure: the author's answers (verbatim)**:
1. On the sea-level platforms: « Elles me semblent réaliste. Il existe des zones completes avec des fond peu profonds et
   d'autres très profonds. Je pense que la continuité du plateau continental dans l'océan est très intéressant. On le voit sur
   la seed 8 par exemple. Le seul petit hic si s'en est un, serait peut-être des fines bandes de terre qui restent émergées je
   ne sais pas si c'est réaliste. Mais dans l'état c'est clairement le mieux qu'on ait eu depuis longtemps. »
2. On the lakes that went: « Oui ça va. Les zones humides étaient intéressantes mais si c'est des artefacts c'est tant mieux
   qu'elles soient plus là. Et d'ailleurs avec les biomes on aura certainement d'autres zones humides apparaitre. Peut-être
   avec les delta, etc. »
3. On the commit: « Oui commit. »
- **Not queued**: a shelf bounded by the wave base instead of the fill up to the sea level. The author keeps the fill up to
  the sea level.

**The check before the commit, seed 23** (the circumnavigation rule: a loss caused by the deposition against the same world
without it; `results/seed23_circ.txt`):
- the same world without the deposition: **not circumnavigable**, with the connected sea (the sea 38.3 %) and with the
  per-cell sea (38.1 %);
- the adopted profile: not circumnavigable (37.1 %).
- So the deposition loses nothing there. The commit goes ahead.

**Queued**:
- the adopted profile in the viz's « Dépôt (niveau de physique) » box, which still runs F182's declaration;
- **for the budget's finding**: the progradation's breadth-first search once per step rather than at each source part.
- **D, before F190**: the thin land strips left emerged on the platforms (the author's « fines bandes de terre »). The
  fields are read without correcting anything.

**What the commit holds** (the reviewer's rule):
- the corrections, the « exported » term, the permanent tests, the re-pinned tests, the re-indexed guard and this record;
- **not** Part 1's instrumentation (the stranding's observers, `physics_history_strand`) nor the F189 benches. The bench
  text needs the instrumentation to compile and is kept untracked at `bench/f189_cascade.rs.txt` (the report's folder). Its
  outputs are in `results/`.
- Seed 8's per-sub-step control (no land cell without a receiver) ran with that instrumentation. The permanent tests cover
  the property on synthetic worlds: `the_ocean_is_linked_through_sides_only`, `a_one_cell_flat_drains_to_its_low_edge`,
  `the_f189_no_stranding`.
