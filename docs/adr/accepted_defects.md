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
