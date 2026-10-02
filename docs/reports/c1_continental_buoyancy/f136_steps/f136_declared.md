# F136 — instruments, declared BEFORE the measurements (2026-10-02)

**World**: the realigned témoin = the viz state "C2 /10 col (défaut)".
- OFF = `new(k, Some(0.1))`, slope floor 0.024; ON = the same with `lake_base: InputLakesAndBasins`.
- Both run with `run_hd`'s tail (`common::viz_hd_lakes_on`). Field = the CONDITIONED field in metres (what the
  viz renders); drainage = the tail's `drainage` (lakes after the crater pass, D8 `flow.direction`).
- **Present lakes**: the ON final lakes ≥ 1 km² that are not crater lakes.

## K5 — the ON water's path
- **Start**: `Lake::outlet`, the lake's last water cell. The true col is its receiver (Finding 119).
- **Walk**: downstream along the ON drainage's `flow.direction`. It stops at the sea (z ≤ 0 m), at the entry into
  another lake, at the grid edge, or after 60 km.
- **Profiles** along those cells: z_ON, z_OFF (conditioned) and z_S1. Beside them, K2's OFF-skeleton path
  (Finding 135), with its distance counted from its point nearest the S1 col.
- SVG, rendered to PNG.
- **"The step exists on the ON path"**: K6's step height on that path is > 50 m.

## K6 — the step below each present lake (on its ON path)
- **Slope over 1 km** at a path cell: (z_i − z_j) / 1 km, j the first path cell ≥ 1 km downstream.
- **The step**: the first STEEP ZONE starting within 5 km of the col.
  - It starts at the first cell whose 1 km slope is ≥ 5 %.
  - It ends at the first later cell whose 1 km slope is < 2 % (the first graded reach).
- **Height** H = lake level − z(end of the zone). **Length** L = the path distance from the col to the end.
- **Max slope over 2 km**: the max of (z_i − z_j) / 2 km over the first 20 km.
- The same quantities in OFF, on the same cells, with OFF's z at the col as the start level.
- **Counts**: H > 50, > 200 and > 500 m.
- **The canyons instrument**: the count of W3's canyons (0 in ON extended, F135) whose floor lies on a step zone.

## K7 — the notch and the lowering
- **The input body**: the input lake (S1's `c1_drainage_windowed`, default config) or the closed depression
  (`basin_base` rule on S1's breached field) with the largest overlap with the final lake.
  - **L_in** = the lowest cell of its outer 8-ring (on S1 for a lake, on breached S1 for a depression).
  - **lowering** = L_in − the final lake's `level_m`.
- **Where**: the final lake's true col (the receiver of `Lake::outlet` in the ON drainage).
- **Which stage**: the true col's z in S1, then after the construction (`carve(S1)`, ON), the eroded world
  (pre-breach) and the conditioned world. The first stage whose z is < L_in − 1 m is the one that lowered the
  sill. "S1" means the true col is the input body's own bed (already under L_in before any construction).
- If it is the construction, the cone rebuild of Finding 135-C says nearest sample or cross-line minimum, and
  gives the laying line's area.
- **The notch** = L_in − the lowest constructed cell of the input body's ring (Finding 135-K1's ring cut, ON).
  Read against the lowering.

## L — the bases on lakes absent at the end
- **The 5 input bodies** with < 50 % of their cells under the ON final lakes (Finding 135-K4).
- **The stage that empties them**: the share of the body under water at:
  - S1 (by definition, 100 %);
  - the construction: fill depth > 0.5 m in `compute_flow(carve(S1))`;
  - the eroded world: the pre-breach `c1_drainage_windowed(None, viz dcfg)` lake map;
  - the final world: the tail's lake map.

  The first stage under 50 % empties them.
- **The floor in the final world** (conditioned):
  - the median slope over the footprint (central differences, m/m);
  - the share under 0.1 % (0.001);
  - the median z minus the base laid (L_in).
- **Their share of F135-M's links**: the OFF skeleton's trunk links ≥ 10 km² that touch the based set but no
  final lake of OFF / F132 / extended (5 175). Split into: those touching these 5 bodies; those touching another
  input lake's cells (the shrunk margins); those touching a closed depression only.

## K8 — contexts available without geology (counts only)
- **Hanging junctions**, on the ON drainage. A junction is a cell with ≥ 2 donors of ≥ 10 km² (accumulation in
  cells × the cell area).
  - The trunk is the largest donor; every other ≥ 10 km² donor is a tributary.
  - Hanging height = z(tributary's last cell) − z(junction cell), conditioned.
  - Counts > 10, > 50 and > 200 m.
  - **Dam-type** (Finding 127): the junction's constructed floor (`carve_diag(S1)`, ON) is laid by a sample whose
    line drains ≤ the tributary's area. Otherwise a **clean break**.
- **Mouths on a coastal wall**: Watercourse segments whose last point's D8 receiver is sea (z ≤ 0).
  - The mouth bed = the conditioned z at the last point.
  - Counted if it is > 1 m **and** within 2 cells of a coastal-wall cell. The wall mask is F126's: carved, not
    floor, above the sea, within 3 cells of the sea, on `carve(S1)` ON.
  - Counts and heights (p50, max).
- **Volcanic terrain**: the mask is the edifices' basal discs (`stamp_volcanic_k`, cells it changes from 1).
  - Trunk cells (accumulation ≥ 10 km²) inside / outside.
  - The slope over 2 km downstream: p50 / p90 / max, inside against outside.
- **K6's steps by context**: outlet of a lake (all of them by construction). Also, within 1 km of the zone: a
  hanging junction > 50 m, a coastal-wall mouth, the volcanic mask.
