# F147 — methods and instruments, declared BEFORE the measurements (2026-10-05)

**The world**: the témoin (C2 /10 col), the viz's HD drainage, as F146.

**Read before declaring** (`f147_xdiag_raw.txt`, a diagnosis of F146's build): the 662 crossings between two smoothed
rivers away from an end are:
- **583 a spillway retracing a watercourse** (they share trace cells);
- **76 two watercourses sharing trace cells**: a confluence's junction cell;
- **3 other**;
- **0 true D8 X.**

## X — the crossings, resolved in the exported geometry only (the hydrology untouched)

1. **The retracing spillways are spliced**: each run of a spillway's vertices whose trace cells belong to one
   watercourse is replaced by that watercourse's own smoothed vertices between the two nearest points. The overlap is
   then exact.
2. **The confluences are frozen**:
   - the receiving river's trace points within one index of each confluence point;
   - and the tributary's second-to-last trace point.
   - These stay vertices (Chaikin keeps them), so the last steps into a junction are the D8 steps, which meet only
     there.
3. **The smoothing constraint becomes lateral**: a vertex may stand at most **0.5 cell** from the D8 trace (pulled back
   towards its nearest trace point by bisection otherwise). Provenance: half a cell keeps the line within the trace
   cells' band. It replaces "1 m above the bed".
- **Gates**:
  - 0 crossings between two rivers away from a shared end (F146's instrument);
  - ridge crossings ≤ 1;
  - lateral deviation from the D8 trace p99 < 0.5 cell.
  - Self-intersections are reported.

## V — the per-vertex attributes (format 0.3.0; the hex edges leave the export)

Each final vertex keeps the trace point it comes from. At that point:
- **the catchment**: the flow accumulation × the cell area (km²);
- **the discharge**: `segment_discharge_profile_m3s`, the climate runoff accumulated per point. **It exists**, so it is
  used.
  - The prorata (the mouth's discharge × catchment ratio, F146's PROXY) is computed for the gap only.
  - **The gap is the share of vertices where they differ by > 20 %.**
- **the slope**: the bed (`profile_m`) over the neighbouring trace points, m/m;
- **the valley width** W(A) = 200 m·(A/10 km²)^0.3: the construction's floor, renamed;
- **the bed width** w = a·Q^b, with b = 0.5.
  - The form is ANCHORED: hydraulic geometry, Leopold & Maddock 1953, USGS Professional Paper 252.
  - a is the author's DECISION. **Default a = 5**, Ymir's own `width_m` law since F22. It is written in the header, so
    a consumer rescales without remaking the world.
- **The table**: for a ∈ {2.5, 3.5, 5, 7}, the share of river km whose bed is wider than:
  - **69.3 m** (one hex flat to flat, the radius read centre to corner);
  - **80 m** (the radius read centre to edge).
- The file's size.

## F — the bundles and spillways, annotated (nothing removed)

- **The parallel pairs**: F146's rule (k = 2 cells, L = 2 km, no shared confluence), computed in the export, excluding a
  retracing spillway's pair with its watercourse.
- **`parallel_of`** = the id of the higher-discharge river of the river's longest pair (its own id when it is that
  river), null without a pair.
- **`fusion_candidate`** = true for the lower-discharge river of at least one pair.
- **`retraces_spillway`** = a spillway whose trace cells are at least 50 % cells of watercourses.
- The viz filters "masquer les candidats à la fusion" and "masquer les déversoirs retracés", with the counts.

## T — the valley-floor instrument

The final polyline's lateral deviation from the D8 trace (p50 / p99 / max, cells), and the ridge crossings (a sample
whose cell drains to another outlet than its trace cell).

## R — the hex reference for Living Landz

- `docs/rivers_ll_hex_reference.json`: 5 varied rivers of the témoin, with their polylines and the hex edges F146's
  computation expects.
- The grid's hypotheses are written in it: radius 40 m centre to corner; hex (0, 0) centred on the bottom-left corner;
  flat-top; axial.
- **The author must confirm them**: the radius at a corner or at an edge's middle; the origin at the centre or at the
  corner of hex (0, 0).

## The cost

The build's and the serialization's seconds per world, against 249.8 s.

## Amendment after run 1 (2026-10-05, NOT blind: run 1 read)

**Run 1's output is kept as defective** (`f147_rivers_raw_run1_DEFECTIVE_instrument.txt`).

**The defect**: on the final polyline, run 1 gave:
- a lateral p99 of **42.7 cells** (max 136);
- **23 ridge crossings**;
- 36 % of samples above the 1 m line.

Meanwhile Chaikin's output gave p99 0.094 and the 0.5-cell bound pulled 0 vertices.

**The cause**: the instruments sample each final edge and look up the trace near ONE vertex's origin (± 3 trace steps).
- After the simplification, an edge joins vertices whose origins are many trace steps apart. The window then misses
  the trace under the edge's middle.
- F146 measured BEFORE the simplification, so the final polyline had never been measured.

**The amendment**: the window becomes the edge's two origins' span ± 3 (lateral) and ± 2 (nearest trace point: bed,
ridge).
- On the raw trace and Chaikin's output, consecutive origins differ by ≤ 1, so the window widens there by at most one
  step.
- Nothing in the build changes.

**Added to the bench**: a local dump (trace cells and vertices within 3 cells) of the first 6 crossings left, so the
71 crossings between confluence pairs can be read before any fix.

## Amendment 2 after run 2 (2026-10-05, NOT blind: run 2 read)

**Run 2** (`f147_rivers_raw_run2.txt`, instruments amended): lateral p99 0.133 (the gate holds). But it left:
- 74 crossings, 71 of them between confluence pairs;
- 88 self-intersections;
- 23 ridge crossings.

**The dump shows the cause, in F146's export itself**: a D8 trace that JUMPS. River 123 goes (4097,1974) -> (4100,1975)
and walks back. River 8781 passes next to itself again 110 steps later.
- `cont_up` chains a reach u onto its receiver d by the largest catchment.
- But the drainage's T3b confluence lets u join d BELOW d's head.
- The concatenated trace then jumps from u's end back to d's head, and walks down d's head reach again: a loop that
  crosses the tributaries joining there.

**The amendment (the export only)**:
- a reach continues its receiver only when it joins it AT ITS HEAD (its last cell is d's first, or drains into it);
- a reach joining below the head ends in a confluence on d's river, at the junction cell.
- Counted and reported: the mid-joins, and those larger than d's reach above the junction (Hack's convention is broken
  there; the alternative would split d, which the segments' 1:1 partition forbids).
- A permanent test with a negative control (the mid-joining reach is the larger one).

