# Finding 147 — the rivers for Living Landz, round 2 of 2 (Ymir's side): `rivers_ll.json` 0.3.0

**Status: built and measured. Nothing committed, nothing selected (the author's).**

**Files:**
- predictions: `f147_predictions.md`, written before any measurement;
- instruments: `f147_declared.md`, with two non-blind amendments;
- final raw output: `f147_rivers_raw.txt` (run 5) and `f147_hexref_raw.txt`;
- defective or superseded runs, kept: `f147_rivers_raw_run1_DEFECTIVE_instrument.txt`, `…_run2_T3b_jumps.txt`,
  `…_run3.txt`, `…_run4.txt`;
- the crossing diagnosis: `f147_xdiag_raw.txt`;
- the format: `docs/rivers_ll_format.md` (0.3.0);
- the hex reference: `docs/rivers_ll_hex_reference.json`.

## Part 0

- F146 committed (`1026b63`) with the ADR's F146 additions.
- The guard and the tests before: `checks_before.txt`.
- **After** (`checks_after.txt`):
  - `cargo check --workspace` clean;
  - lib 605 passed (601 + 4), 0 failed;
  - viz 30 passed;
  - **the guard 6 / 6, field and lakes = banc, through `run_hd`**; C2 /10 col field = `a8d2d538d692c2f0`.
  - In production, `rivers_ll` takes 0.23–0.47 s per state.

## The run history (what the bench showed, in order)

1. **Run 1: the instrument was defective.**
   - It gave a final lateral p99 of 42.7 cells, while Chaikin gave 0.094 and the bound pulled no vertex.
   - The cause: the instruments looked up the trace near ONE vertex's origin. After the simplification, an edge spans
     origins far apart, so the lookup missed.
   - F146 measured before the simplification, so the final polyline had never been measured.
   - Amended (the window = the edge's origin span); declared non-blind.
2. **Run 2: the trace jumps.** The lateral gate held (p99 0.133), but 74 crossings, 88 self-intersections and 23
   ridge crossings were left. A local dump showed traces that **jump** back.
   - **A defect of F146's export**: `cont_up` chained a reach onto its receiver even when the reach joins it BELOW its
     head (the drainage's T3b confluence).
   - The trace then jumped back to the receiver's head and walked its head reach again: a loop crossing the
     tributaries there.
   - Amended: only a reach joining at the head continues it; declared non-blind; a permanent test with a negative
     control.
3. **Runs 3–5: the crossings are gone.** Runs 4 and 5 add the diagnosis of the ridge samples; the build is the same as
   run 3.

## X — the crossings

**The diagnosis before any fix** (F146's build, `f147_xdiag_raw.txt`). The 662 crossings were:
- 583 a spillway retracing a watercourse;
- 76 two watercourses sharing trace cells;
- 3 other;
- **0 true D8 X.** F146's reading ("the D8 network crosses itself, diagonal X's") is wrong.

**The method**, in the exported geometry only (the drainage untouched):
- the spillways' runs over watercourse cells are spliced onto the watercourses' own vertices;
- the confluence neighbourhoods are frozen (a step between two fixed vertices stays the straight D8 step);
- the T3b chaining is fixed (amendment 2);
- a 0.5-cell lateral bound.

**Results (run 5):**

| measure | F146 | F147 |
|---|---|---|
| crossings between two rivers away from an end, not a touch at a shared vertex | 662 | **0** |
| crossings through a shared run of vertices (arriving on one side, leaving on the other) | — | **0** (8 566 shared runs) |
| self-intersections | 271 | **0** |
| any intersection away from an end, touches included (F146's raw instrument) | 662 | 1 335 (the touches at shared vertices, overlaps and confluences) |

**Gate "0 crossings outside a confluence": HELD.**

- Spillway trace points spliced onto watercourse vertices: 850 of 1 175 (72.3 %).
- **T3b**: 2 261 reaches join their receiver below its head and now end in a confluence on it.
  - **831 of them are larger than the receiver's reach above the junction.** Hack's convention is broken there: the
    larger stream ends on the smaller one.
  - The alternative splits the receiver's reach at the junction. That breaks the 1 : 1 partition of `rivers.json`'s
    segments, so it is the author's decision.
- The river count goes from 9 786 to 10 246.

## T — the valley floor, on the final polyline

| measure | value |
|---|---|
| lateral deviation from the D8 trace, p50 / p99 / max (samples every 0.25 cell) | 0.000 / **0.133** / 0.221 cell |
| (Chaikin before the bound, p99 / max) | 0.094 / 0.133 |
| vertices pulled back by the 0.5-cell bound | 0 of 3 313 939 |
| ridge crossings (a sample's cell draining to another outlet than its trace cell) | raw trace 640 · **final 23** (F146: 1) |
| stair index (turns ≥ 45°) | raw 0.193 → final **0.089** (F146: 0.045) |
| F146's retired "1 m above the bed" | raw 16.96 % · Chaikin 22.47 % · final 23.54 % |

**Gate "lateral p99 < 0.5 cell": HELD.**

**Gate "ridge crossings ≤ 1": FAILED as declared (23).**
- The diagnosis: all 23 samples are OFF the trace, and **21 of 23 lie within 0.02 cell of a cell corner** (max
  0.0365).
- 16 sit exactly on an integer corner.
- The frozen D8 steps (new in F147) include diagonal steps. A straight diagonal step passes EXACTLY through the corner
  shared by four cells, and the instrument's `floor` assigns that corner point to a cell off the trace.
- So the line touches an off-valley cell at a point, over zero length.
- **This is an instrument at a boundary of measure zero, not a line leaving its valley.** But the gate was declared on
  this instrument, so it fails as declared.
- Amending it (e.g. ignoring samples within 0.05 cell of a corner) is for the reviewer.

**The price of the fixed vertices**: the stair index doubles against F146 (0.045 → 0.089). The frozen steps are D8
steps, so at each confluence and each shared cell the line keeps a D8 turn.

## V — the per-vertex attributes; format 0.3.0

**Removed:** the hex edges (and each river's count).

**Added**, per vertex: catchment, discharge, slope, valley width W(A), bed width a·Q^0.5. The header states the laws
and the discharge's source.

- **The discharge is the accumulated one**: `C1DrainageResult::segment_discharge_profile_m3s`, the climate runoff
  accumulated per trace point (drainage.rs:363). It exists, so it is used and cited.
- **Gap to F146's prorata** (mouth Q × A/A_mouth): **18.1 % of the vertices differ by > 20 %.**
  - prorata / accumulated p10 / p50 / p90: 0.64 / 0.96 / 1.00;
  - watercourses 17.8 %; spillways 73.1 % (672 vertices: a spillway's catchment says nothing of its lake's outflow);
  - 2 883 vertices carry no accumulated discharge (Q = 0).
- Per vertex, p10 / p50 / p90 / max:
  - catchment 6.84 / 54.05 / 2 851 / 72 520 km²;
  - discharge 0.075 / 1.05 / 49.1 / 5 300 m³/s;
  - slope 0.0064 / 0.051 / 0.50 (negative 0.0 %; the T3b loops made 0.7 % in run 1);
  - valley width 178 / 332 / 1 090 / 2 879 m;
  - bed width (a = 5) 1.4 / 5.1 / 35 / 364 m.
- **File size: 11.05 MB** (F146 without edges: 5.1 MB; with: 31.9 MB).
  - The five per-vertex arrays (4 significant digits, ~5 characters each over 122 868 vertices, ≈ 3.5 MB) make most
    of the increase.
  - Reducing it is a format choice (e.g. fewer digits on the points), not made: the gated geometry is the one
    exported.

**The bed width w = a·Q^0.5:**
- **b = 0.5 is ANCHORED** (Leopold & Maddock 1953, USGS PP 252).
- **a is the author's DECISION.** The default is a = 5 (`BED_WIDTH_A_DEFAULT`, Ymir's own channel width), in the
  header.

The share of river km whose bed is wider than one hex:

| a | Q needed (69.3 m / 80 m) | all rivers (21 541 km): > 69.3 m / > 80 m | watercourses only (21 476 km) |
|---|---|---|---|
| 2.5 | 768 / 1 024 m³/s | 0.35 % / 0.13 % | 0.14 % / 0.00 % |
| 3.5 | 392 / 522 m³/s | 0.63 % / 0.48 % | 0.39 % / 0.24 % |
| **5 (default)** | 192 / 256 m³/s | **1.11 % / 0.82 %** | 0.87 % / 0.59 % |
| 7 | 98 / 131 m³/s | 2.26 % / 1.68 % | 2.01 % / 1.42 % |

The two hex readings:
- **69.3 m** is flat to flat with the 40 m radius read centre to corner;
- **80 m** is flat to flat with the radius read centre to edge.

## F — the annotations (nothing removed)

- **1 631 parallel pairs** (≤ 2 cells over > 2 km, no shared confluence, no spillway with the watercourse it
  retraces). F146: 1 599.
- 1 423 rivers carry a `parallel_of`; **938 are fusion candidates** (3 448 km).
- **5 of the 15 spillways retrace watercourses** (≥ 50 % of their trace; 53 km).
- With both viz filters on, 9 303 of 10 246 rivers are left (18 041 of 21 541 km).
- The viz layer "Rivières LL" gets the filters « masquer les candidats à la fusion (n) » and « masquer les déversoirs
  retracés (n) ». The hex-edge count leaves its label.

## R — the hex reference

`docs/rivers_ll_hex_reference.json` holds 5 rivers of the témoin. Each pick is the first river by id in its class, not
already picked:
- a watercourse ending at the sea;
- a tributary at a confluence;
- a river ending in a lake;
- a spillway;
- a large river.

**For each river**, the file gives:
- its polyline in cells and in metres;
- the hexes of its two ends;
- the expected hex edges in order (F146-H's computation: walk every 4 m, each hex change is an edge).

**The assumptions, written in the file, which the author must confirm:**
- **the 40 m radius is centre to CORNER** (else, centre to edge, every edge changes);
- **hex (0, 0) is CENTRED on the map's bottom-left corner** (else every index shifts);
- flat-top; axial; y northward.

The picks (`f147_hexref_raw.txt`; the file is 47 kB):

| case | river | length | vertices | expected hex edges |
|---|---|---|---|---|
| a watercourse ending at the sea, 2–5 km | 5 | 2.95 km | 7 | 84 |
| a tributary ending at a confluence, 1–3 km | 0 | 2.57 km | 22 | 45 |
| a river ending in a lake, 1–5 km | 2049 | 2.25 km | 4 | 39 |
| a spillway, 0.5–5 km | 10237 | 2.14 km | 8 | 40 |
| a large river (catchment ≥ 50 km²) under 8 km | 2 | 2.61 km | 3 | 38 |

The first run picked river 0 twice (it is both a tributary and ≥ 50 km²). The picks now exclude a river already picked.

A negative r is normal: in flat-top axial, r = −x/3 + …

## The cost

**Build 0.37 s + serialization 0.03 s = 0.40 s per world** (+0.16 % of 249.8 s). F146: 0.43 s + 0.15 s for the hex
edges, now gone.

## Predictions

**The reviewer's (hypotheses to check):**
- **X, "0 crossings, no new ridge"**: 0 crossings HELD. "No new ridge" FAILED on the declared instrument (1 → 23, all
  corner touches).
- **V, "< 10 MB"**: REFUTED (11.05 MB).
- **V, "the accumulated discharge departs from the prorata by > 20 % on ≥ 10 % of the vertices"**: HELD (18.1 %).
- **Bed width, "< 5 % of km > 69 m for usual a"**: HELD (≤ 2.26 % up to a = 7).
- **Cost, "< 1 s"**: HELD (0.40 s).

**Mine:**
- **P-X**:
  - "mostly not true X's": HELD (0 true X).
  - "true X's 20–40 %": REFUTED (0 %).
  - "they reach 0 only by a topological choice": moot, there were none.
  - The real topological defect was another one (T3b), unforeseen.
- **P-V**: "< 10 MB" REFUTED; "the accumulated discharge exists" HELD; "> 30 % of vertices" REFUTED (18.1 %).
- **P-bed width**: < 2 % at a = 5 HELD (1.11 %); < 5 % up to a = 7 HELD.
- **P-cost, 0.6–1.5 s**: REFUTED (0.40 s, below).
- **P-T**: lateral p99 < 0.1 REFUTED (0.133); ridges ≤ 1 REFUTED (23).
- **Meta, "at least one of mine is wrong"**: HELD.

## For the author

**Decisions:**
- the coefficient `a`;
- the radius reading and the hex origin (the reference file);
- the selection, by the flags and the length / catchment / Strahler filters;
- the 831 T3b junctions where the larger stream ends on the smaller one: split the receiver's reach, or keep it;
- whether the corner-touch ridge instrument is amended.

**What remains is your look in Living Landz.**
