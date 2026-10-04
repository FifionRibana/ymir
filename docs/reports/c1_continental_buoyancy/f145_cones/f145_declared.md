# F145 — instruments, declared BEFORE the measurements (2026-10-04)

**The worlds** (construction stage only: S1 → `skeleton` → `carve`; no production code is added, so the candidate is
not run through the light pass, the droplets or the breach):
- **the témoin**: the viz state C2 /10 col, `ValleyConstruction::new(F121_AGE_K, Some(0.1))`;
- **ON extended**: + `lake_base: InputLakesAndBasins`;
- for M2: **B2 → A_c** (F127's canyons' world): `a_min_km2: 0.1` on the témoin.

## The basin label (declared)

- **A cell's basin** = the cell at which its D8 path, on the skeleton's own `direction`, reaches a terminal (DIR_NONE,
  which includes every sea cell). Memoised.
  - It is the outlet basin of the field the skeleton drains, the breached input: the same object as `compute_flow`'s
    `basins`, which the skeleton computes and drops.
- **A line's basin** = the majority basin of its polyline's sample cells.
- **A sample is foreign to a cell** when its line's basin differs from the cell's basin.
- **A divide cell** = a land cell with an 8-neighbour of another basin. The distance to a divide is in cells
  (`dist_from`).

## M (bench `f145_m`)

- **The own-basin carve** (the candidate C1, also M1's reference): a test-only copy of `carve_diag`'s path for this
  configuration (no confluence band, no wall profile, no wall-sea floor, no rim; asserted).
  - A sample may seed or propagate into a cell only if it is not foreign to that cell.
  - The cross-line minimum takes only non-foreign neighbours' samples.
  - **A cell no own-basin sample reaches keeps its terrain (uncarved)**: the restriction must not invent a floor
    from nowhere.
  - **Control**: with the filter off, the copy equals production `carve_diag` bit for bit (else nothing is read).
- **M1**: for every carved cell whose laying line (`who`) is foreign: Δz = z_prod − z_own (m).
  - The distribution (p5 / p25 / p50 / p75 / p95), and by the distance to a divide (0–2, 3–10, > 10 cells).
  - **The real defect** = Δz < −1 m: its count and share of the carved cells.
  - **The legitimate** = within 2 cells of a divide with |Δz| ≤ 1 m.
- **M2**: in B2 → A_c, F127's four true cols, (3852, 2337), (4551, 3280), (2253, 4495) and (1805, 4580):
  - foreign-laid or not;
  - Δz;
  - in the real defect or not.
  - F144-Z's shares in G-drained are recalled (99.3 % / 94.0 %).
- **M3**:
  - **lowered divides**: divide cells with z_prod < z_own − 1 m (count; the lowering's p50 / p90);
  - **artificial captures**: each field's drainage is the skeleton's flow chain (the pre-drainage, `breach_monotone`,
    `compute_flow`) on the input, on z_prod and on z_own. A cell's route is labelled by the INPUT basin of the last
    land cell on its D8 path.
  - A cell is captured artificially when its route under z_prod differs from the input's, while under z_own it equals
    the input's. Count, km², and where (the largest pairs of input basins).
- **Stop rule** (the brief's): if the real defect is < 1 % of the carved cells and there is no capture, report and do
  not run C.

## C — the candidate C1 = the own-basin carve (bench `f145_m`), against the témoin (rule 18)

On the témoin and on ON extended, **at the construction stage**:
- canyons (F95's criteria on the construction, as the extras);
- **F127's dams**: the four cols' heights in B2 → A_c under the témoin's carve and under C1, and its canyons;
- **walls across the divides**: adjacent cell pairs of different input basins with a slope > 28°. **Negative
  control: the témoin's carve.**
- captures (M3's measure);
- **θ** on the carved links (the construction's mask), **tolerance ±0.01** (declared);
- R8 terrain (Finding 125's windows, the témoin's classes);
- R8 network (chord 8, on the field's own rivers);
- the relief (F95's paired p50);
- the coast (spurs near a coastal wall);
- **the cost**: the filtered carve's seconds against production `carve`'s, on the same skeleton, plus the labelling.
  Against 249.8 s per world.
- No gate without its negative control. No promotion.

## Amendment to M3's captures, declared NON-BLIND after f145_m's témoin (before the refined bench runs)

f145_m's capture count on the témoin is **2 635 697 cells (6 284 km²) in 6 334 basin pairs**, and 717 905 cells
change their route even under the own-basin carve.
- The label (the INPUT route of the last land cell before the sea) separates the coast's many tiny basins: a carved
  mouth shifted by a few cells moves a cell to a neighbouring mini-basin.
- **The count is read as an upper bound, not as captures.**

**Refined** (bench `f145_cap`, the témoin and ON extended):
- an input basin's area = the land cells whose input route ends at its outlet;
- **a capture counts only when both the input basin and the basin it is routed to under the témoin's carve are
  ≥ 10 km²**, and the own-basin carve keeps the cell in its input basin;
- **the noise floor**: the same count for the own-basin carve against the input (cells it moves between basins
  ≥ 10 km²);
- the largest pairs (origin, destination, km²), with their location.
