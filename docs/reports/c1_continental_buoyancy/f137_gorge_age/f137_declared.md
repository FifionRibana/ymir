# F137 — instruments, declared BEFORE the measurements (2026-10-02)

Every world is built at the canonical framing and run through `run_hd`'s tail (`common::viz_hd_lakes_on`,
climate 45° / 40°). The field is the CONDITIONED one, in metres.

## Q1 — the steps below the lakes in the delivered worlds
- **Worlds**:
  - "livré" (`Knobs::passes(2)`, the delivered stream power);
  - "A1+B2" (+ `slope_floor_abs 0.024`, the world F114 measured, at today's framing);
  - ON extended (the construction), for comparison.
- **Lakes**: the final lakes ≥ 1 km², not crater lakes, with a D8 receiver at `Lake::outlet`. The number
  without one is given.
- **Instrument: F136-K6b exactly.**
  - Path: from the true col down the world's D8, to the sea, another lake, the edge, or 30 km.
  - Drop: from the lake's level to the first cell, within 10 km, whose slope over the next 2 km is < 2 %.
  - The distance to that cell is given, and the drop over the first 5 km.
  - Counts > 50 / > 200 / > 500 m.
- **F114's world and stage** are read from its bench (`tests/f114_sill.rs`): `Knobs { slope_floor_abs: Some(0.024),
  ..passes(2) }` at the framing of its time, read on the eroded and breached fields.

## Q2 — C-3's `production_k_field`
- **The field**: `production_k_field(state, kin, cfg, edifices, volc)` on the viz state's world, exactly the one the
  incision reads.
- **Components rebuilt beside it**:
  - the lithology classes (`build_coarse_k` → `upscale_k_to_hd`);
  - the volcanic stamp (`stamp_volcanic_k` on a field of ones, cells it changes);
  - the fracture factor (`build_hd_density_k`).
- **Over land (conditioned z > 0)**:
  - quantiles of each component;
  - a histogram of K;
  - the shares volcanic / rift-soft / fracture ≠ 1.
- **Correlation length**: the autocorrelation along x and along y (pairs on a 1-in-8 lattice, land pairs only), at
  lags 1, 2, 4, …, 2 048 cells. The length is the first lag below 1/e, for K, log K, each component, and the
  altitude.
- **Correlations**: the Pearson r of log K on land (1 cell in 7) with:
  - the altitude and the local slope (central differences);
  - the precipitation (`c1_climate_placed` 45 / 40);
  - the volcanic mask;
  - the cratonic mask (nearest-sampled at the altitude's mapping);
  - the rift-soft class;
  - the fracture factor.
- **The stages that read it**, from the code: `incise_with_floor` (`stream_power.rs:751`), called by the delivered
  incision and by the construction's light pass (`production_upscale.rs:509–553`); `carve` does not read it.
- **Map**: a PNG at 1/8 resolution, north up, log K on land, volcanic discs white.

## Q3 — F136's steps: lengths and classes
From F136's raw outputs (`f136_k.txt`, `f136_k6b.txt`; non-blind):
- **length**: the declared instrument's L (col → end of the first steep zone), and the amended instrument's
  distance (col → the first graded reach ≥ 2 km);
- **mean slope** = drop / length;
- **class** on the mean slope: wall > 45°, steep 10–45°, gorge < 10°.

## Q4 — the below-sea basins' spillways (ON extended)
- **Lakes**: the final lakes ≥ 1 km² with id in 1 000 000–1 999 999 (the below-sea merge) and no D8 receiver at
  their outlet.
- **Path**: their Spillway segments (`segment_kind == Spillway`, `segment_source_lake == id`), the first one and
  its downstream chain, then the D8 from its last point to the sea or 30 km.
- **The same K6b instrument** from the lake's level, the total drop to the path's end, and Q3's class.

## A — read, not measured
- The grep (rules 11 / 11c) is in `grep_rule11_age.txt`.
- Each answer cites the code line that carries the age.
