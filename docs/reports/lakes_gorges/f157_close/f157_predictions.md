# F157 — my predictions, written BEFORE any measurement (2026-10-08)

**Not blind, declared.** Before writing these I had read:
- **F133's table and F135-W's realigned replay** (OFF / F132 / extended):
  - canyons 0 / 0 / 0; spurs near a coastal wall 0.0142;
  - R8 terrain 0.0452 / 0.0395 / 0.0411; relief p50 515.4 / 551.9 / 546.6 m;
  - Δz between flat resolutions 1 150 588 / 235 602 / 12 509;
  - lakes 26 / 34 / 34 on the viz path;
  - θ on the carved links 0.493 → 0.436 [0.427, 0.445] (0.439 extended).
- **F133v's list**: the 11 lakes new in ON (1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 19), their crops, 615 cells each.
- **F136's K6b**: steps > 50 / > 200 / > 500 m = 8 / 6 / 2; lake 10 811 m; lake 2 862 m; lakes 7, 1, 9, 11 at
  260–317 m; max slope over 2 km 0.35.
- **F155's production world** (OFF): relief 515.4 m, R8 terrain 0.0451, planar walls 9.5 % of the land, sharp crests
  2 739, comb tile R8 0.0883.
- **The history since F135**: no change to the construction, the light pass, the droplets or the breach in
  production. F146–F149 changed the HD drainage bundle (types, rivers, an outlet rule: `ALGO_HD_DRAINAGE` 9 → 12).

## C — the gorge's closing summary

- Nothing to predict. A summary, not a measurement.

## V — the lakes' outline

- **P-V1**: the guard stays 6 / 6, field and lakes (the outline is a view).

## B — the lake base alone

- **P-B1**: F133's terrain columns are reproduced **to the digit**: canyons, coast, R8, relief p50, Δz between
  resolutions and θ, for OFF and extended. The lake counts are **26 / 34**. If anything moves, it is a lake count or
  area, through F146–F149's bundle changes.
- **P-B2**: the 11 new lakes are the same 11, at the same centres (±2 cells).
  - A control crop with zero changed cells: **none exists** at 615 × 615. ON moves the conditioned field over most
    of the land (44 M cells at F133v).
  - With the declared |Δz| ≤ 1 m threshold, a control exists far from the based lakes.
- **P-B4 (the steps)**: K6b's numbers are reproduced: 8 / 6 / 2 over 50 / 200 / 500 m.
  - ≥ 4 steps > 200 m;
  - ≥ 1 with a max slope over 100 m above 40°;
  - mean slopes of 2–15° over the step.
- **P-B5**: the two exports are written as the production world's `.ymir` containers (format 1.1.0, with rivers_ll
  and the geology).
- **P-B6**: promoting the extended lake base:
  - changes the production field's hash (the construction's χ changes), hence the field guard of every state with
    the construction (4 of 6) and their lake guard;
  - adds ≈ 0–10 s per world (the χ walk's extra stops), within the noise.
- **P-B7 (the production defects, OFF crops)**:
  - comb teeth (F124's definition) on the walls in ≥ 6 of the 12 crops;
  - planar walls 5–20 % of the land per crop;
  - the terrain R8 per crop 0.03–0.10.
  - Stages: the teeth and the planar walls are born at the **construction** (the walls, the D8-lattice trunks); the
    axis alignment comes partly from the **construction** (R8 rises from S1) and is cut by the light pass.

## Meta

- At least one of these predictions is refuted.
