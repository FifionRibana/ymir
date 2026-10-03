# F141 — instruments, declared BEFORE the measurements (2026-10-03)

**Worlds**:
- **ON ×1** = the témoin with `lake_base: InputLakesAndBasins`, the gorge off.
- **GORGE ×1** = ON + `gorge_retreat: v3(1, 0.5)` (F140's crashing state).

Both are built through `build_world` and run_hd's pre-assembly steps. The protected breach is replayed by an
instrumented copy (below).

## A — the three definitions of A (the 22 lakes with a body: F139's 14 D8 + 8 spillway)
- **The lakes**: ON ×1's final lakes ≥ 1 km², not crater lakes, with a D8 outlet or a spillway path (F139's set).
  The body is the construction's (GORGE ×1's `gorge_bodies`) with the largest overlap of the final footprint.
- **(i)** the body's `a_out_km2`: the col's drained area in the construction (F140).
- **(ii)** F139's instrument: max(the ON accumulation on the path's first km, the final lake's inflow + its area).
- **(iii)** the inflow into the construction's input body: Σ of the skeleton's `area_km2` over the cells outside the
  body whose skeleton receiver is in it, plus the body's area.
- For each, at p = 0.5 and A_ref = 418.7 km²:
  - r_lake at r_world 0 / 1 / 1.5 / 2 (×0.7 / ×1 / ×1.2 / ×1.4, map (a));
  - the emptying age, map (a) of r_world = 2 / (A/A_ref)^0.5 ("never" above 2).
- **The divergence** (i) / (ii) > ×2 and lake 4: the body's outflow cell and col; its largest inflow (the outside
  cell with the largest area entering the body), its area, and whether it touches the FINAL lake's footprint
  (within 1 cell).

## Br — the breach that leaves land under the sea
- **An instrumented copy of `breach_monotone_protected`** (`flow.rs`), test-only, line for line. For every carve it
  records the pit that triggered it (`nb`) and the ramp's start (`height[nb] − EPS`).
  - **Control**: its output must equal the production function's **bit for bit**, or nothing is read from it.
- **The population**: the land cells (eroded > sea) that end ≤ sea after the breach.
- **For each**:
  - its pit: the pit's floor, and its fill level from the pre-breach `flow.filled` (its spill);
  - the base its ramp reaches (the backlink chain): the sea (z ≤ sea) or a pre-breach lake (id);
  - the number of ramp steps from the pit;
  - whether a final lake covers it (ON: run_hd's tail; GORGE: the tail panics, so the pre-breach lake map and the
    panic are said).
- Clustered by pit: the top pits by the number of cells they lower below the sea, with their location.
- **The central question** is answered from the code (lines cited) and checked on the records: whether the ramp
  starts at the pit's floor or at its spill, and where it stops.

## Z — the cell 78 m under body 4's drained floor
- The cell (4148, 3674) of GORGE ×1. Its z at:
  - S1;
  - the construction (`carve_diag(S1)`);
  - (a) the pipeline's construction + rims;
  - (b) + the light pass;
  - (c) + the droplets;
  - (d) + the bathymetry;
  - the protected breach.
- The construction's laying sample (`who`): its position, floor zf and half-width, its line, the skeleton's base and
  χ at its cell, and whether that cell is in body 4.
- The first stage under L_floor (168.9 m) names the function.

## D — read from the spec v3 and F138 / F139's tables (no measurement)

## Amendment, declared after f141_attr's Part A and BEFORE measuring it (2026-10-03)
f141_attr maps every one of the 22 lakes to a construction body. **Body 4 is none of them**: F139's lake 4 is body 3
(A 300.7 / 298.1 / 305.2). F140's "1 571 against F139's 298 for the same lake" paired body 4 with lake 4 by index.
A complementary bench, `f141_b4`, therefore measures:
- the ON final lake(s) over body 4's footprint: their id, area, type, level and outlet;
- why that lake is not in F139's set: whether its outlet has a D8 receiver, and whether it has spillway segments;
- body 4's (i), (iii) and its col's donors, as in Part A;
- (ii), computed on that lake with F139's instrument where it applies.
