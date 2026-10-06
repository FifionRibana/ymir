# `rivers_ll.json` — rivers for Living Landz (format 0.4.0)

A vector layer of the `.ymir` container (layer id `rivers_ll`), written next to `rivers.json`, which is unchanged.
It holds the same drainage as `rivers.json`, regrouped into **rivers** (main stems), each a **smoothed polyline that
stays in its valley floor**, with the attributes a selection needs. **No selection is applied**: every river is
exported, and the consumer picks.

Produced by `ymir_core::export::rivers_ll` (ADR Findings 146–147).

## Changes

- **0.4.0** (Finding 148): **Hack's convention holds at every confluence.**
  - Where a reach joins another below its head (the drainage's T3b confluence), the larger stream continues and the
    receiving reach is cut there.
  - So a `rivers.json` segment may now belong to two rivers, and `segments` lists the segments a river covers wholly
    or in part (the 1 : 1 match is dropped).
- **0.3.0** (Finding 147):
  - **removed**: the hex edges (`hex_grid`, `hex_edges` and each river's `hex_edges` count). Living Landz computes them
    from the polylines, so the grid lives in one place. A reference for that computation is
    `docs/rivers_ll_hex_reference.json`.
  - **added**: per-vertex attributes (`vertex`); the width laws and the discharge's source in the header; the
    annotations `parallel_of`, `fusion_candidate`, `retraces_spillway`.
  - The polyline no longer crosses another river away from a confluence (see "How the polyline is made").
- **0.2.0** (Finding 146): the hex edges.
- **0.1.0**: the first layout.

## Coordinates

- `coordinate_space`: `"erosion_grid_cells_continuous"`, i.e. continuous cell coordinates on the container's rasters
  (`manifest.continent.grid`).
- **A cell's centre is at `(x + 0.5, y + 0.5)`.** `x` grows eastward and **`y = 0` is the SOUTH edge**, the
  container's orientation invariant.
- `km_per_cell` converts to kilometres.

## Document

```json
{
  "format_version": "0.4.0",
  "coordinate_space": "erosion_grid_cells_continuous",
  "km_per_cell": 0.0488,
  "discharge_source": "accumulated: the climate runoff … per trace point (…)",
  "bed_width_law": { "form": "a * Q^b", "a": 5.0, "b": 0.5, "a_status": "decision: …", "b_status": "anchored: …" },
  "valley_width_law": { "form": "coef_m * A_km2^exp", "coef_m": 100.24, "exp": 0.3, "status": "…" },
  "rivers": [ RiverLl, ... ]
}
```

A reader MUST refuse an unknown *major* version.

### The header's laws

- **`bed_width_law`**: the bed width is w = a·Q^b, with Q the mean discharge.
  - **b = 0.5 is ANCHORED**: downstream hydraulic geometry, Leopold & Maddock 1953 (USGS Professional Paper 252).
  - **a is a DECISION.** The default 5.0 is Ymir's own channel width (`rivers.json`'s `width_m`).
  - Every exported `bed_width_m` uses the header's `a`. A consumer who wants another coefficient rescales by
    `a_new / a`; no world needs regenerating.
- **`valley_width_law`**: the valley construction's floor width W(A) = 200 m·(A / 10 km²)^0.3 (a PROXY, Finding 120).
- **`discharge_source`**: the per-vertex discharge is the climate runoff (precipitation − evapotranspiration)
  accumulated down the drainage. It is not a prorata of the mouth's discharge.

## `RiverLl`

| field | type | meaning |
|---|---|---|
| `id` | u32 | the river's id (its index in `rivers`) |
| `kind` | `"Watercourse"` \| `"Spillway"` | a spillway is the outflow of a closed below-sea basin over its col: no hierarchy |
| `segments` | [usize] | the `rivers.json` segments it covers, wholly or in part, upstream → downstream |
| `points` | [[f32, f32]] | the smoothed polyline, upstream → downstream |
| `vertex` | object | per-vertex attributes, parallel to `points` (below) |
| `length_km` | f32 | the polyline's length |
| `catchment_km2` | f32 | geometric catchment at the mouth |
| `strahler` | u8 \| null | the largest Strahler order among its segments; `null` for a spillway |
| `discharge_m3s` | f32 | mean discharge at the mouth |
| `width_m` | f32 | `rivers.json`'s channel width at the mouth (5·Q^0.5) |
| `end` | object | how it ends, below |
| `lakes_crossed` | [u32] | the `lakes.json` ids touching its course, other than the one it ends in |
| `falls` | [Fall] | tagged falls (`kind`, `x`, `y`, `height_m`); **empty while the gorge work is paused** |
| `length_rank_in_basin` | u32 | 1 = the longest river of its mouth basin |
| `mouth_basin` | u32 | the mouth basin's key: the last land cell (row-major index) on the D8 path from the river's end |
| `parallel_of` | u32 \| null | the higher-discharge river of this river's longest parallel pair (its own id when it is that river); `null` without a pair |
| `fusion_candidate` | bool | the lower-discharge river of at least one parallel pair |
| `retraces_spillway` | bool | a spillway with ≥ 50 % of its D8 trace on watercourse cells; those stretches are drawn on the watercourse's vertices |

### `vertex`

Columnar: equal-length arrays, one entry per point of `points`. Each value is read at the D8 trace point the vertex
comes from, and rounded to 4 significant digits.

| array | type | meaning |
|---|---|---|
| `catchment_km2` | f32 | drained area: the flow accumulation × the cell area |
| `discharge_m3s` | f32 | mean discharge, accumulated (see `discharge_source`) |
| `slope` | f32 | the bed's slope along the river over the neighbouring trace points (m/m, positive downstream; a spillway climbing to its col reads negative) |
| `valley_width_m` | f32 | the construction's valley-floor width W(A) (`valley_width_law`) |
| `bed_width_m` | f32 | the bed width a·Q^b (`bed_width_law`) |

### `end`

`end` is one of:

- `{"type": "sea"}`;
- `{"type": "lake", "lake_id": n}`;
- `{"type": "endorheic_lake", "lake_id": n}`;
- `{"type": "confluence", "river_id": n}`: the river joins a larger one. **Its last point is exactly a vertex of
  river `n`.**
- `{"type": "terminal"}`.

### Parallel pairs

Two rivers form a pair when:
- their polylines, sampled every ≤ 1 cell, run within 2 cells of each other (Chebyshev, on the samples' cells);
- over more than 2 km of consecutive samples;
- they do not share a confluence;
- and they are not a retracing spillway with its watercourse.

Nothing is removed. The flags are for the consumer's selection.

## What "a river" is

**A main stem (Hack's convention).** At each confluence, the upstream reach with the largest catchment continues the
river downstream. Every other reach ends its own river there, with `end = confluence`.

The drainage lets a reach join another below its head (its T3b confluence). The receiving reach is then cut at the
junction, so the larger of the two streams continues and the other ends there. A `rivers.json` segment may thus be
shared by two rivers.

## How the polyline is made

1. **The D8 trace**: the cell centres of the chained segments. At a confluence, the receiving river's cell is
   appended.
2. **Fixed vertices**:
   - both endpoints;
   - every confluence point;
   - the receiver's trace points either side of a confluence point;
   - the tributary's second-to-last point;
   - every cell a spillway shares with a watercourse.

   A step between two fixed vertices stays the straight D8 step. So a tributary meets its receiver only at the
   confluence vertex.
3. **Chaikin corner cutting**, 3 iterations, keeping the fixed vertices.
4. **The valley bound**: a smoothed vertex stands at most **0.5 cell** from the D8 trace. A vertex beyond it is pulled
   back towards its nearest trace point.
5. **Spillways over watercourses**: each stretch of ≥ 2 trace cells a spillway shares with a watercourse is replaced by
   that watercourse's own vertices. The two lines then overlap exactly.
6. **Ramer–Douglas–Peucker**, ε = 0.1 cell, keeping the fixed and the shared vertices.

Rendering as SDF strokes (like roads) is the intended use. `vertex.bed_width_m` sizes the stroke along the river.

## Hex grids

The export carries no hex data since 0.3.0. Living Landz derives hex edges, fords, obstacles and bonuses from the
polylines on its own grid.

`docs/rivers_ll_hex_reference.json` gives 5 rivers with their polylines and the hex edges Ymir's reference computation
finds. It states its grid assumptions (40 m flat-top, axial, hex (0, 0) centred on the bottom-left corner), and the
author is asked to confirm them.
