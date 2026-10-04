# `rivers_ll.json` — rivers for Living Landz (format 0.2.0)

A vector layer of the `.ymir` container (layer id `rivers_ll`), written next to `rivers.json`, which is unchanged.
It holds the same drainage as `rivers.json`, regrouped into **rivers** (main stems), each a **smoothed polyline that
stays in its valley floor**, with the attributes a selection needs. **No selection is applied**: every river is
exported, and the consumer picks.

Produced by `ymir_core::export::rivers_ll` (ADR Finding 146).

## Coordinates

- `coordinate_space`: `"erosion_grid_cells_continuous"`, i.e. continuous cell coordinates on the container's rasters
  (`manifest.continent.grid`).
- **A cell's centre is at `(x + 0.5, y + 0.5)`.** `x` grows eastward and **`y = 0` is the SOUTH edge**, the
  container's orientation invariant.
- `km_per_cell` converts to kilometres.

## Document

```json
{
  "format_version": "0.2.0",
  "coordinate_space": "erosion_grid_cells_continuous",
  "km_per_cell": 0.0488,
  "rivers": [ RiverLl, ... ],
  "hex_grid": HexGrid,
  "hex_edges": HexEdges
}
```

A reader MUST refuse an unknown *major* version.

## `RiverLl`

| field | type | meaning |
|---|---|---|
| `id` | u32 | the river's id (its index in `rivers`) |
| `kind` | `"Watercourse"` \| `"Spillway"` | a spillway is the outflow of a closed below-sea basin over its col: no hierarchy |
| `segments` | [usize] | the `rivers.json` segments it chains, upstream → downstream |
| `points` | [[f32, f32]] | the smoothed polyline, upstream → downstream |
| `length_km` | f32 | the polyline's length |
| `catchment_km2` | f32 | geometric catchment at the mouth |
| `strahler` | u8 \| null | the largest Strahler order among its segments; `null` for a spillway |
| `discharge_m3s` | f32 | mean discharge at the mouth |
| `width_m` | f32 | bankfull channel width at the mouth (w = 5·Q^0.5); a stroke hint |
| `end` | object | how it ends, below |
| `lakes_crossed` | [u32] | the `lakes.json` ids touching its course, other than the one it ends in |
| `falls` | [Fall] | tagged falls (`kind`, `x`, `y`, `height_m`); **empty while the gorge work is paused** |
| `length_rank_in_basin` | u32 | 1 = the longest river of its mouth basin |
| `hex_edges` | u32 | the number of hex edges it crosses (its entries in `hex_edges`) |
| `mouth_basin` | u32 | the mouth basin's key: the last land cell (row-major index) on the D8 path from the river's end |

`end` is one of:
- `{"type": "sea"}`;
- `{"type": "lake", "lake_id": n}`;
- `{"type": "endorheic_lake", "lake_id": n}`;
- `{"type": "confluence", "river_id": n}`: the river joins a larger one. **Its last point is exactly a vertex of
  river `n`.**
- `{"type": "terminal"}`.

## What "a river" is

**A main stem (Hack's convention).** At each confluence, the upstream reach with the largest catchment continues the
river downstream; every other reach ends its own river there, with `end = confluence`. Every `rivers.json` segment
belongs to exactly one river.

## How the polyline is made

1. **The D8 trace**: the cell centres of the chained segments, with the receiving river's cell appended at a
   confluence.
2. **Chaikin corner cutting**, 3 iterations. Both endpoints and every point where a tributary joins stay fixed.
3. **The valley constraint**: a smoothed vertex may stand at most **1 m above the bed** of its nearest trace point
   (bilinear terrain on the conditioned field). A vertex that fails is pulled back towards that trace point.
4. **Ramer–Douglas–Peucker**, ε = 0.1 cell, keeping the fixed points.

Rendering as SDF strokes (like roads) is the intended use; `width_m` sizes the stroke.

## The hex grid and its edges (since 0.2.0)

`hex_grid` describes Living Landz's grid (the author, 2026-10-05):

```json
{ "size_m": 40.0, "orientation": "flat_top", "coordinates": "axial", "origin_m": [0.0, 0.0] }
```

- `size_m` is the hex's radius, centre to corner.
- **Hex (0, 0) is centred on the map's bottom-left corner** (`origin_m`, metres from that corner).
- `y` points north (the container's `y = 0` = south).
- A point at metres `(x, y)` lies in the axial hex obtained by `q = (2/3·x) / size`, `r = (−x/3 + √3/3·y) / size`,
  cube-rounded.

`hex_edges` lists every crossing of a hex edge by a river, **columnar** (equal-length arrays, one entry per crossing,
in river order then along the river):

| array | type | meaning |
|---|---|---|
| `q1`, `r1` | i32 | the hex the river leaves (axial) |
| `q2`, `r2` | i32 | the hex it enters (a neighbour of the first) |
| `river_id` | u32 | the river |
| `catchment_km2` | f32 | the catchment at the crossing (the flow accumulation, capped at the river's mouth) |
| `discharge_m3s` | f32 | the mouth's discharge scaled by the catchment ratio |
| `slope` | f32 | the local slope along the river over ±1 cell of arc (m/m, positive downstream) |
| `width_m` | f32 | the construction's floor width W(A) = 200 m·(A / 10 km²)^0.3 |

Ford, obstacle and bonus rules are Living Landz's to decide; Ymir only exports the edges.

