# F148 — instruments and methods, declared BEFORE the measurements (2026-10-06)

**The world**: the témoin (C2 /10 col), the viz's HD drainage (`viz_hd_lakes_on`), as F146–F147.

**Before / after**:
- the guard 6/6 is `ymir-viz f123_viz_guard` (field and lakes through `run_hd`);
- C2 /10 col field `a8d2d538d692c2f0`;
- lib and viz tests, `cargo check --workspace`.
- "Before" = F147's after-checks at `254c66f` (the same tree).

## J — Hack and the ridges (the export only)

**J1 — the receiving reach is cut at each T3b junction.**
- `rivers_ll` builds on PIECES: each segment is cut at every index where another reach joins it below its head.
- A piece's upstream pieces all join at its head. The continuing one is the largest catchment:
  - an inner piece's catchment = the accumulation at its last cell, scaled to the segment's;
  - a last piece's = the segment's.
- So Hack's convention holds at every junction by construction.
- A river's `segments` lists the segments it covers wholly or in part. The 1 : 1 match with `rivers.json` is dropped
  (decided).
- **Gate**: 0 Hack violations.
  - The instrument: for each river ending in a confluence on R, its mouth catchment against R's per-vertex catchment
    at the vertex before the junction vertex.
  - A violation is when the tributary's is larger by more than 0.1 %.
  - Plus 0 crossings and 0 self-intersections (F147's instruments).

**J2 — the ridge instrument amended**: a crossing is a polyline stretch of length **> 0.001 cell** (5 cm) inside a
cell that drains to another outlet than the river's nearest trace cell.
- Each final edge is cut exactly at the cell boundaries.
- A zero-length touch at a cell corner is not a crossing.
- The count is per (edge, cell) stretch.
- **F147's sample instrument is reported beside it, on the same build.**

## M — the convergence measurement (the reviewer's hypothesis to check)

**The field**: the conditioned field the network is traced on, in metres.

**The index**: **C = div(∇z / |∇z|)**, the contour (planform) curvature, in km⁻¹.
- C > 0 is convergent (a hollow: flow lines meet); C < 0 is divergent (a nose); C ≈ 0 is a planar slope.
- Computed on the field smoothed by a Gaussian of **σ = 2 cells (98 m)**.
  - The reason: the cell-scale curvature is noise.
  - 2 cells is the narrowest constructed valley floor (W(A_c) ≈ 50 m is one cell; W(1 km²) ≈ 100 m is two).
- Sensitivity at σ = 1 and 4 cells is reported. The decision uses σ = 2.
- |∇z| < 1e-4 (flat) gives C = 0.
- A unit test with a negative control: a cone pit C > 0, a cone peak C < 0, a plane C ≈ 0.

**The A·S² criterion** (Montgomery & Dietrich, the reviewer's citation, from memory, source to check):
- A in m² (the accumulation × the cell area);
- S the D8 slope to the receiver on the same field.

**Populations** (rivers of the J build; a river in (c) is removed from (a) and (b), overlap reported):
- **(a)** the fusion candidates (F147: the lower-discharge river of a parallel pair). Sensitivity: every pair member.
- **(b)** the micro-rivers into a lake: `end` = lake or endorheic lake, and length **< 1 km** (below F146's median river
  length of 1.2 km; ~20 cells).
- **(c) the control**: length **≥ 10 km**, OR **≥ 50 % of its trace cells on a construction valley floor** (`carve`'s
  floor mask, the témoin's construction).

**Scores per river**, over its own trace cells (its last cell, which belongs to the receiving body, excluded):
- **max** — the rule's statistic: a river with one eligible cell is kept;
- **head** — the mean over its first 10 trace cells (blind to length);
- **median** — along the trace.

**Separation**: AUC = P(score(c) > score(a ∪ b)), ties counted half, for each score, for C and A·S².
- The max-score AUC is confounded by length (declared in the predictions). The head AUC is the honest separation.

**The rule simulated** (`build_rivers_ll` with an optional, gated head rule; off by default; not wired):
- a network cell is **eligible** if A ≥ A_min = A_c = 0.1 km² (production's head) and index ≥ θ;
- a segment is channel from its first eligible cell, or from the junction index of its first channel upstream reach,
  whichever comes first;
- **a reach starting on a lake's shore (an outlet) or a spillway is channel from its start**: water leaving a lake is a
  channel;
- downstream of a channel cell, everything is channel (the D8 network's closure).
- `dr.rivers` (rivers.json, lakes, outlets) is NOT touched: G shows the lakes read the segments.

**The table**, for each θ:
- θ = 0, plus the quantiles 0.5 / 1 / 2 / 5 / 10 / 20 % of (c)'s max score;
- reported: the rivers and km remaining; the parallel pairs remaining; the micro-lake rivers remaining; (a) and (b)
  removed; (c) lost.
- A river of a population is **removed / lost** when < 50 % of its trace cells remain channel cells.

**Stop rule (the brief's)**: if no θ removes ≥ 50 % of (a) + (b) while losing ≤ 1 % of (c), report and do not build B.

## K — the lake 1000011

On the témoin (`viz_hd_lakes_on`, the same assembly as `run_hd`):
- the lake's type and unresolved reason, its provenance (id ≥ 1 000 001: a below-sea basin);
- whether a spillway was traced for it;
- the segments whose first point lies within 1, 2 and 3 cells (Chebyshev) of its footprint, their kind and length;
- whether the lake is Exorheic before `resolve_exorheic_without_outlet`.

**Stage attribution**: the first stage at which a reach starting on its shore stops existing. **The viz panel's two
rules** (closed = Endorheic only; a 2-cell search) are compared with the core's (1 cell).

## The cost

The rule's and the convergence field's seconds per world, against 249.8 s.
