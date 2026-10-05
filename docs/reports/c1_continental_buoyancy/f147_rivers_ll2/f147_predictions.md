# Finding 147 — predictions, written 2026-10-05 BEFORE any measurement of this round (the crossing diagnosis was running, unread)

**Reading declaration, unfavourable.** Non-blind on Findings 73–146 and the reviewer's hypotheses, in particular:
- **The geometry of a D8 "X"**: cells P = (x, y) → Q = (x+1, y+1) for one river, R = (x+1, y) → S = (x, y+1) for the
  other. **The orthogonal neighbours of either diagonal are the OTHER river's cells**, so rerouting one river through
  an orthogonal neighbour makes it pass through a vertex of the other, on both sides of it: still a crossing.
  - When R and S lie on opposite sides of the other river's course, the crossing is **topological**. No local geometry
    removes it without changing which river joins which.
- **C1DrainageResult carries `segment_discharge_profile_m3s`**: the discharge per point of each segment, accumulated
  from the climate's runoff (F22 / F42). It is a per-point accumulated discharge, to be confirmed in the code before
  use.
- F146: 9 786 rivers, 112 698 points; 5.1 MB without the hex edges; spillways retracing watercourses exist.

## Predictions

- **P-X**:
  - the 662 crossings are mostly **not** true D8 X's. **Spillways retracing a watercourse, two rivers sharing trace
    cells, and crossings near a river's end make more than half.** True X's are 20–40 %;
  - the non-X classes are removed by the export's geometry (snapping a retrace onto its watercourse, near-end
    repairs);
  - **the true X's reach 0 only with a declared topological choice** (the smaller river ends at the larger at the X and
    its remainder becomes its own river). Against the reviewer's "0 without a new ridge crossing" by a local reroute;
  - ridge crossings stay ≤ 1, and the lateral p99 stays < 0.5 cell.
- **P-V**:
  - the file without hex edges but with the per-vertex attributes is **< 10 MB** (with the reviewer);
  - **the accumulated discharge exists** (`segment_discharge_profile_m3s`) and departs from the prorata by > 20 % on
    **> 30 %** of the vertices (the reviewer says ≥ 10 %; agreeing, larger).
- **P-bed width** (w = a·Q^0.5):
  - at a = 5, Ymir's own `width_m` law, **< 2 %** of the river km are wider than 69 m (it needs Q > 190 m³/s);
  - < 5 % for a up to 7 (with the reviewer).
- **P-cost**: **0.6–1.5 s** per world in total. Moving the parallel-pair detection and the crossing repair into the
  export may pass 1 s (against the reviewer's < 1 s).
- **P-T**: the final polyline's lateral deviation p99 < 0.1 cell; ridge crossings ≤ 1.

## Meta

My disagreements with the reviewer: X (topological X's), cost. **At least one of my predictions is wrong.**
