# F133v — the viz path, declared BEFORE the run (2026-09-30)

**Reading declaration, unfavourable.** Non-blind: the author's viz read 25 lakes OFF (badge Match, key
cf1d4539c7afe588, field a8d2d538d692c2f0) and 33 ON (NoReference, key 834fa66b41d6e05d, field a620a9aa0d5882a0).

- **Path**: `run_hd` itself (the viz's worker: `breach_monotone_protected`, its climate, its drainage), called from a
  `#[cfg(test)]` module of `ymir-viz` (as `f123_viz_guard`): no production code changes.
- **Parameters**: `f123_viz_guard`'s literal ("the workspace's own literal"): 8192², domain 400 km, latitude 45°,
  span 40°, manual offset [6/64, 37/64], stream power + closures + MFD p 2, slope floor 0.024, base level on,
  geo ratio 7.5, volcanism / lithology (10, 3, 1) / fracture (6, 25 km) / infiltration ON; the valley as the
  workspace builds it for "C2 /10 col (défaut)": `ValleyConstruction::new(k, Some(0.1))`, γ None, and
  `lake_base = Some(InputLakesAndBasins)` for ON. The author's own latitude sliders are NOT in the logs (climate
  HIT): the guard's literal is used and said so.
- **1 — Counts**: `drainage.lakes.len()` OFF and ON. Beside it, the BENCH assembly on the same eroded field
  (`breach_monotone`, climate placed at 45 / 40, the benches' `dcfg`): if the counts differ, a table of the lakes
  that differ, with the cause isolated by swapping ONE stage at a time (the breach: protected vs plain; the
  drainage config: the viz's vs the benches').
- **2 — Crops** (615 cells = 30 km, viz cell coordinates, the data row): (a) each ON lake with < 50 % of its cells
  under any OFF lake; (b) the 3 based lakes with the largest upstream |Δz|: each cell where the eroded field
  differs OFF vs ON is attributed to the first ON final lake its OFF D8 path (the run's flow direction) enters;
  **statistic: Σ |Δz| × cell area (km³) over the cells attributed to the lake** (p90 |Δz| beside); (c) a control
  crop: the 615 × 615 window (64-cell steps) with the MOST land among those with ZERO cells whose eroded z differs
  OFF vs ON (certified by that count).
- **3 — Lakes per crop**: the distinct ids of `drainage.lake_map` with >= 1 cell in the crop, OFF and ON.
- **4 — Badges**: `result.bench_guard` OFF and ON, with the key and the field hash.
- **5 — The guard's lake debt**: named from the code, not built.

## Amended BEFORE the run
- `HdResult.eroded` is the CONDITIONED field (post-breach, what the viz renders): Δz OFF vs ON and the control crop are
  read on it; the attribution walks the run's own `drainage.flow.direction` (OFF) to the first ON final lake.
- The bench-assembly comparison is run ONLY if the counts differ from 25 / 33 (the author's item 1), as a follow-up.
- "Lacs en base" for 2(b) = the ON world's final lakes (every present lake is a base under the extended variant).
