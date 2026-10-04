# Finding 144 — predictions, written 2026-10-04 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–143 and the reviewer's predictions, in particular:
- **The code of H_f** (`valley_construction.rs`, the gorge loop): `let hf = (1.5 * sg * 100.0).max(10.0);` with
  `sg = slope_for(top)`, which is the real slope max(m·S_loi, S_req), capped.
  - **H_f already uses the real slope, not 10·S_loi.** The ADR item recorded at F143's commit ("probably 10·S_loi")
    is read here, not measured.
  - F143-L located G-tag's defect in the instrument's span: two cells, 195–276 m, against H_f's 100 m basis.
- **F143-S**: the proxy's 18 = the 14 + bodies 5, 9, 11 and 15, all "lone input lakes" (their ring touches no
  depression cell, no sea cell). Body 15 is 20 % under 1000016 in ON.
- **F143-L / F143-C (G-sea)**:
  - lake 4's outlet runs to the open sea at (2034, 3343);
  - its L_floor is 13.4 m;
  - the 4.3 km inlet comes from the breach's ramps from pits of its bowl (13.8–25.5 m, spill 61.8 m);
  - it is already present on the construction's breach (+111), and under C2 (+122).
- F143-H: the 14's head falls at ×1, p 0.5 (349, 300, 179, 60, 37, 34, 30, 14, 13 m; 9 tagged).

## Predictions

- **P-S**: the input criterion (no cell in a `basin_base` closed depression, none under the sea) keeps **18** bodies:
  the 14 plus 5, 9, 11 **and 15**. Body 15 is not inside a depression on the input.
  - Against the reviewer's 17.
  - In the gorge world 5, 9 and 11 are present lakes (with the reviewer).
  - **G-bodies against the gorge world's own lakes fails for at least one body at ×1.4** (a drained or merged body:
    lake 13 joins 1000016's water, F142).
- **P-B**:
  - **Lake 4's B is the sea (0 m)**, below its L_floor (13.4 m). **The bound does not bind it**, and the inlet stays
    at ×1.4, p 0.5 (against the reviewer).
  - The bound binds the lakes whose outlet enters another water body above their floor: **lake 13**, whose next base
    is 1000016 at ~115 m, from ×1.2 on.
- **P-H**:
  - at ×1, p 0.5, 9–12 tagged head falls on the 17, none above its H_cap(A) (≤ 131 m);
  - **H_f was not the cause of G-tag's defect** (the code already uses the real slope);
  - with G-tag's span corrected, **at least one gorge drop still passes**: a hanging junction, as F143's 591 → 202 m.
- **P-P** (phase 1):
  - **P1, P2 and P3 all hold G-pits on the three worlds**: each keeps the footprints at or above the construction
    (P1, P3) or does not erode them (P2).
  - **P1 reproduces C2's walls** (G-ring > 100; with the reviewer).
  - **P2 has the fewest edge cells and is selected** (against the reviewer's P3).
  - P3 makes edge walls where the design is laid back above the eroded surroundings.
  - **No candidate holds G-sea at ×1.4, p 0.5** (the inlet is born at the construction's breach).
- **P-θ** (the eroded world, gorge corridors excluded): **below OFF's eroded CI at ×1** (against the reviewer). The
  light pass spreads the rebased links (F143: 10–15 % outside [0.5, 2] against 4 %).
- **P-Z** (not blind on the count; F143-L: 2 192 of 2 207 by a line "draining elsewhere"):
  - **more than 90 %** of G-drained's violators at ×1.4 are laid by a cone from outside their catchment;
  - over the drained bowls' cells and over the whole world the share is small (< 10 %).

## Meta

My disagreements with the reviewer:
- S (18, not 17);
- B (the inlet stays);
- G-tag (one drop stays);
- P (P2, not P3);
- θ.

**At least one of my predictions is wrong.**
