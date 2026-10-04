# F146 — instruments and methods, declared BEFORE the measurements (2026-10-04)

**The world**: the témoin (C2 /10 col), the viz's HD drainage: the run_hd tail's assembly on the conditioned (breached)
field. These are the rivers `rivers.json` exports.

## The rivers (the unit of `rivers_ll.json`)

**A river is a main stem.** At each confluence, the upstream segment with the largest catchment continues the
downstream river; every other upstream segment ends its own river there (Hack's convention). Spillways are rivers too
(`kind` kept). Every segment belongs to exactly one river.

## L — the smoothing in the valley

- **The method: Chaikin corner cutting**, 3 iterations, on the river's chained cell-centre points (x + 0.5, y + 0.5).
  **Both endpoints and every confluence point stay fixed.** A tributary's last point is its confluence on the main
  stem, which stays a vertex.
- **The valley constraint**: each smoothed vertex must satisfy terrain(bilinear, conditioned field) − bed(the nearest
  trace point's `profile_m`) ≤ **1 m**.
  - Tolerance provenance: the ±1 m of the level gates used since F140.
  - A vertex that fails is moved towards its nearest trace point by bisection (8 steps); the trace point itself always
    satisfies it.
- **The instruments**, on the smoothed polyline sampled every 0.25 cell, before and after the constraint:
  - terrain − bed (p50 / p99 / max);
  - the share of length above the tolerance;
  - **ridge crossings**: a sample whose cell drains (the drainage's D8) to another outlet than the river's own trace
    cell;
  - **the lateral deviation** from the D8 trace (p50 / p99, cells);
  - crossings between two rivers away from a shared confluence, and self-intersections (a segment-intersection test
    on a 4-cell hash).
- **Negative control**: the raw D8 trace, sampled the same way (in the valley by construction, in steps). Its "stair
  index" is the share of vertices where the direction turns by ≥ 45°.

## P — the parallel bundles

- **A pair**: two rivers without a shared confluence, whose sampled points stay within **k = 2 cells** of each other
  over a run longer than **L = 2 km**.
  - Provenance: k ≈ the trunk floor's half-width at 10 km² (W/2 = 100 m ≈ 2 cells); L = 2 km excludes incidental
    approaches near junctions.
- **Per pair**: its length, the smaller river's catchment, and where the run lies:
  - on the construction's planar walls (carved, not floor);
  - on a valley floor;
  - on a volcanic edifice (within the craters' records' radius);
  - elsewhere.
- **Proposed export rules (not applied)**: (a) merge each pair into the higher-discharge river; (b) drop the rivers
  below an area threshold. Each with what it removes (count, km).

## A — the selection attributes (per river)

- the length (km);
- the catchment at the mouth (km²);
- the Strahler order (the maximum over its segments);
- the discharge at the mouth (m³/s);
- the downstream end: `sea` / `lake` (with its id; `endorheic_lake` when the lake is endorheic) / `confluence` (into
  a larger river) / `terminal`;
- the lakes crossed: the lakes touching its points, other than the end lake;
- the falls (the gorge's format; empty while it is paused);
- the length rank within its mouth basin.
- **No threshold is fixed.**

## H — the hex edges: blocked

Ymir knows no Living Landz hex grid (Q.3). The parameters are asked of the author.

## E — the export and the preview

- **`rivers_ll.json`**, `format_version` "0.1.0", documented in `docs/rivers_ll_format.md`. **`rivers.json` stays
  byte-identical** (a test).
- The viz layer "Rivières LL" with filters on length, area and Strahler, and counts.
- **The cost**: the build's seconds per world, against 249.8 s.
