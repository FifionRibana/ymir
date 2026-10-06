# Accepted defects

A measured behaviour that a Finding named as a defect, that the author then decided to **accept rather than fix**.
Each entry gives:
- what it is;
- the measurement;
- the reason it is accepted;
- its exception, if any;
- the date and the Finding.

An entry is removed only by a later Finding that reopens it.

## 1. `carve_diag`'s cones lay cells across the pre-carve divides (accepted 2026-10-04, F146 closing F145)

- **What**: the valley construction gives a cell its nearest sample whose 28° cone lies below the terrain
  (`valley_construction.rs`, `carve_diag`). Nothing bounds a cone but the terrain, so a valley's wall lays cells on the
  far side of the divide measured before carving.
- **Measured (F145, the témoin C2 /10 col, the construction stage)**:
  - 18.2 % of the carved cells (23.0 % in ON extended, 9.4 % in B2 → A_c) are laid by a line of another basin and
    deepened by > 1 m against the own-basin carve (Δz p50 −568 m);
  - 64 % of them lie within 2 cells of the old divide;
  - divides lowered: p50 425 m;
  - real captures between distant outlets: ~70 km² (~0.3 % of the land), below the measure's noise floor (~140–160
    km²).
- **Why it is accepted**:
  - the "defect" is measured against the divides of before the carving;
  - at the old divide the ridge becomes the intersection of the two valleys' 28° walls: the normal geometry of a
    dissected range;
  - the captures are below the measure's noise;
  - the own-basin restriction (C1) adds a seam: walls across the divides ×1.49 (×3.85 in B2 → A_c).
- **Exception, not accepted**: the drained bowls (F144-Z: 94–99 % of G-drained's violators are laid from outside
  their catchment). That belongs to the gorge's retreat, which is paused.

## 2. The smoothed river polylines clip the corners of off-valley cells (accepted 2026-10-06, F149 closing F148-J2)

- **What**: `rivers_ll.json`'s Chaikin-smoothed polylines (`export/rivers_ll.rs`) stay within 0.22 cell of the D8
  trace. In places they cut across the corner of a neighbouring cell that drains to another outlet.
- **Measured (F148-J2, the témoin)**: each final edge is cut exactly at the cell boundaries.
  - **134 stretches** run inside an off-outlet cell (zero-length corner touches excluded), all shorter than 0.1 cell:
    p50 0.012, max 0.098 cell = 4.8 m.
  - **122 m in all over 21 541 km of rivers.**
  - The 0.25-cell sample instrument of F146–F147 saw 23 of them.
- **Why it is accepted**:
  - the stretches are below the scale of one Living Landz hex (~70–80 m across) by an order of magnitude;
  - the river never leaves its valley floor (lateral p99 0.133 cell);
  - removing them would need either a second amendment of the instrument (a length tolerance), which the reviewer
    declined, or pinning more vertices to the D8 steps, which brings back the staircase.
- **No exception.**

## 3. Parallel bundles and micro-rivers into lakes stay in the network (accepted 2026-10-06, F149 closing F146–F148)

- **What**: the river designation is by area alone (`flow.rs`, `extract_rivers`: ≥ 20 km², plus a main-stem tail up
  to A_c = 0.1 km²). It keeps:
  - **parallel bundles**: two rivers within 2 cells of each other over > 2 km;
  - **micro-rivers into lakes**: shorter than 1 km, ending in a lake, which read as runoff on a shore.
- **Measured (F146–F149, the témoin)**:
  - 1 631–1 660 parallel pairs, 942 fusion candidates (3 430 km);
  - 726 micro-rivers into a lake (351 km).
  - **The heads of both are planar** (the mean contour convergence over their first 10 cells: median −0.12 and
    −0.07 km⁻¹, against +1.26 for a control of 2 828 rivers; AUC 0.784 at σ = 2 cells, 0.857 at σ = 1).
  - **But no head rule passes the stop rule** (≥ 50 % of them removed for ≤ 1 % of the control lost):
    - F148: one convergent cell anywhere, 2.4 % at 0.99 %;
    - F149: the head score, 3.0 % at 0.60 % (σ = 2), 7.2 % at 0.64 % (σ = 1);
    - A·S² separates the wrong way.
  - The control rivers a rule would lose first are real rivers (p50 8.2 km, 8 of 17 ≥ 10 km) whose heads rise on
    planar uncarved ground.
- **Why it is accepted**: the separation is real at the head but the control's tail overlaps it, so dropping them
  drops real rivers too. The F148 stop rule was not relaxed after the fact (rule 15).
- **Mitigation**: `rivers_ll.json` flags them (`fusion_candidate`, `micro_lake_inflow`). The viz box « Masquer
  faisceaux et micro-rivières » hides 1 668 rivers (3 781 km) on the témoin; Living Landz can do the same.
- **No exception.**
