# F150 — instruments, declared BEFORE the measurements (2026-10-06)

**No production code this round.** Benches and documents only.

**The world**: the témoin as the benches build it.
- `build_world(Knobs { valley: C2 /10 col, slope_floor_abs: S_EQ, ..Knobs::passes(2) })`, seed `PSEED`, the
  canonical framing.
- **C-2 volcanism, C-3 lithology and C-3b fracture are ON**: the benches' default, and the 6/6 guard's.
- The coarse state: `run_coarse_tectonics(PSEED, 64, Phase2InitParams::default(), the benches' 300-step run,
  C1Closures::default())`, the same derivation production erodes.

**Before / after**: the guard 6/6 (`ymir-viz f123_viz_guard`), lib and viz tests, `cargo check --workspace`.
"Before" = F149's after-checks (`9e710be`).

## I — the inventory

- **One table row per field**, from the code:
  - name, type, unit;
  - computed from what, and at which stage;
  - filled on the témoin;
  - exported (which file);
  - read by whom.
- **"Filled" is measured** on the coarse state and the HD fields:
  - the share of non-default cells;
  - min / p50 / max;
  - the classes' counts.
- **The maps**:
  - one PNG per filled field, 1024², north up (internal y = 0 is south, so the rows are flipped);
  - the sea in grey;
  - the coarse fields sampled at the témoin's framing (the canonical origin, size 1): nearest for classes, bilinear
    for continuous fields.
  - The fields: crust thickness `s`, `age`, `plate_type`, `plate_id`, `cratonic_mask`, the boundary classes (labels),
    the C-3 class (hard / rift-soft / volcaniclastic), the C-3b density, the edifices.
- **The three empty panel fields**: answered from the code (where they would come from, why they are empty).

## D — does the rock already change the relief?

- **A** = the témoin (C-3 ON); **B** = the same with `lithology_off` (C-3 OFF; C-2 and C-3b stay ON). This is the
  viz box « Lithologie (C-3) » alone.
- **Δz = z_A − z_B** (m, the eroded field), over the land cells of either.
  - |Δz| p50 / p90 / max, and the signed p1 / p99.
- **Where**: the share of the |Δz| > 10 m cells:
  - on the rift-soft class;
  - on a volcanic footprint (within a crater record's radius × 2);
  - within 5 km of either;
  - elsewhere.
- **The relief** (F116's instruments): the land altitude p50; R8 = `aniso(f, land, 16).r8`; the 3×3 local σ p50.
  A against B.
- The world-time cost is not part of this round.

## C — the catalogue

`catalogue_geologie_v1.md`. The real deposit contexts are cited from textbook-level sources:
- Robb 2005, *Introduction to Ore-Forming Processes*;
- Evans 1993, *Ore Geology and Industrial Minerals*;
- the USGS deposit models (Cox & Singer 1986).
- These are cited from memory and **labelled to check**, since the PDFs are not in `docs/refs`.

**Each resource's zoning** is graded "existing fields", "existing + a derived field", or "missing field". The
fields are named.
