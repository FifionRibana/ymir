# Piste 4 — gallery-dynamics production + P95-cap render classification + permanent morphology acceptance

Retire the reclassify-driven workflow as the C1 production path; make gallery dynamics + the P95-cap sea level on the **render** path the production state, with a **permanent morphology acceptance gate**. Successor to Issue #141 (which delivered the `SeaLevelMode` mechanism but whose workflow calibration produced filaments — see `docs/reports/phase_1_5/README.md` post-merge correction).

## Branch + discipline

Branch `TBD-piste-4-gallery-production-morphology` off `milestone/c1-lightweight-dynamic-tectonics` **after this spec's branch (with the #141 corrective) merges**. Piste 4 **reuses** #141's `SeaLevelMode`/`c1_default` on the render path — does NOT re-implement or alter it.

- **Production C1 = gallery dynamics** (Issue #137 contract): standalone `run_with_closures`, **NO per-cycle `apply_post_tectonic` reclassify**. (Stage S confirmed: reclassify is invoked ONLY by the viz workflow toggle + tests; there is no headless C1 production entry. Gallery is already the de-facto production default.)
- **KEEP all #139 code** (hover, hypsometric, workflow worker, sliders, E4 continuation). Workflow toggle → **default OFF + relabeled "experimental — reclassify-driven migration is filamentous, non-production."** It is the instrument that proved the dead end + the **scaffold for the eustatic follow-up**. Do NOT remove.
- **Morphology metric is PERMANENT acceptance** (perim/area + largest-component + n_components), NOT fraction alone. The fraction-only gate let 30%-filaments through.
- **NO C1 export work here.** There is no C1 export — `export/mod.rs` is the v2 pipeline (unsuited for C1, a known separate chantier). The export land/sea classification at the low threshold comes WITH the future C1 export, NOT in piste 4. (This lightens Stage E2.)
- **Eustatic-coast is a SEPARATE follow-up** — register only; do NOT build here. And it requires its own first measurement (smooth glide vs block-jumps).
- v2 byte-identical + 9th-bit-identical guards hold (piste 4 = render/classification + acceptance + relabel; does not touch v2 or the core sea-level mechanism beyond reuse).

## Methodological context (memory)

- `reclassify-coast-dead-end-use-eustasy` — reclassify-coast is a 3×-proven paradigm dead end; production = gallery + P95-cap render; dynamic coast = eustasy (later, measured).
- `morphology-must-gate-acceptance` — the new gate; a scalar fraction does not validate spatial structure.
- `gallery-vs-workflow-reference` — gallery is the production reference.
- `visual-validation-supersedes-scalar-metrics` — the visual MUST actually run, not just be listed.

## Per-stage workflow

### Stage S — audit (~0.25d, mostly pre-confirmed)
Pre-confirmed by inspection: (1) reclassify only in viz workflow toggle + tests, no headless C1 prod entry → gallery is de-facto production; (2) NO C1 export (export/mod.rs is v2). Remaining to confirm + record: viz default pipeline = Gallery (`ActivePipeline` default); render altitude already uses `c1_default` (c1_viz `derive_altitude_field`, #141 E2); gallery worker uses `c1_default` (#141 E2). Establish baseline gallery-production morphology (seed 42 + 2 seeds). Commit `DOC Stage S`.

### Stage E1 — permanent morphology metric (~1d) — THE DURABLE ADD
Promote the throwaway morphology fn to a tested helper: `land_morphology(mask: &[bool], nx, ny) -> { area_frac, perimeter_over_area, n_components, largest_component_frac }` (4-neighbour). This is the net that would have caught the filaments. Unit-test on canonical shapes (single blob → low p/a, high largest, 1 comp; scattered/checkerboard → high p/a, many comp). Place where both acceptance tests and a future C1-export validator can use it.
W7: deterministic; tested on canonical shapes; connectivity documented. Commit `FEAT Stage E1`.

### Stage E2 — confirm production = gallery + render reuse (~0.5d, lightened)
Confirm production path = gallery dynamics (no reclassify) — already the default. Confirm `SeaLevelMode`/`c1_default` is on the render altitude (#141 E2) — **do NOT re-implement**. Register the C1 export (+ its land/sea classification at the low threshold) as a SEPARATE chantier (comes with the C1 export, which needs the v2 export pipeline reworked for C1). No code wiring of export here.
W7: production = standalone `run_with_closures`; render uses `c1_default` (reused, not re-implemented); C1-export registered separate, nothing wired. Commit `FEAT Stage E2` (or fold into E3 if trivially small).

### Stage E3 — demote workflow toggle (~0.5d)
Default pipeline = Gallery; workflow toggle relabeled experimental/dead-end; keep worker + sliders + E4 continuation (eustatic scaffold); hover + hypsometric untouched.
W7: default Gallery; workflow reachable + labeled + non-default; scaffold + instruments preserved. Commit `FEAT Stage E3`.

### Stage V — morphology + regression (~0.5d)
Headless morphology test, gallery seed 42: **perim/area < ~0.6, largest-component > ~0.6, n_components < ~20** (gallery ref 0.515 / 0.625 / 11, with per-seed margin) + fraction ~30%. v2 byte-identical (enumerate the FULL baseline — `enumerate-dont-tail-regression-baseline`; the 2 pre-existing failures: #47 fixed on its branch, `rectangular` tracked separately). 9th-bit-identical: gallery production does not exercise the `apply_post_tectonic` decomposition; the contract test (under `c1_default`, from #141) still passes — confirm unaffected.
W7: morphology gate per-seed margin; v2 full-baseline enumeration; 9th-bit unaffected. Commit `TEST Stage V`.

### Stage A — acceptance + VISUAL THAT RUNS (~1d)
Headless multi-seed: morphology + fraction band across ~5 seeds (natural variation, all masses). **Run the visual** (Viz-0.5 Gallery, seed 42 + 2 seeds): continental MASSES (not filaments), ~30%, static coast — the visual MUST execute (this is the correction of the filament gap; a listed-but-unrun visual is what let it through). The fraction-only acceptance is REPLACED by morphology + fraction.
W7: multi-seed morphology band; the visual is actually executed + recorded; fraction-only acceptance retired. Commit `TEST Stage A`.

### Stage Final — docs + memory + PR (~0.5d)
README `docs/reports/piste_4/`: reclassify-coast dead end (3×), eustatic follow-up **with its smooth-glide measurement prereq**, gallery-production + morphology gate, `SeaLevelMode` reused on render, workflow demoted, C1-export registered separate. Memory: create `project_piste_4_outcomes`; link the captured lessons (`reclassify-coast-dead-end-use-eustasy`, `morphology-must-gate-acceptance`, `isolate-subcomponent-not-block`, the corrected `c1-phase-1-5-outcomes`). Register follow-ups: **eustatic-coast (measure smooth-glide vs block-jumps FIRST)**, **C1 export + low-threshold land/sea classification** (needs v2-export reworked for C1), Phase 3 prerequisite (workflow quasi-static topology — now moot for production; relevant only if eustatic reuses the worker), gallery PNG regeneration, `rectangular` v2-solver issue. v2 + 9th-bit guards. PR to milestone.

## Anti-patterns
NE PAS reclassify per-cycle in production · NE PAS remove #139 workflow code (keep OFF/experimental, scaffold) · NE PAS accept on fraction alone (morphology gates + the visual runs) · NE PAS build eustatic-coast here · NE PAS create a partial C1 export here (separate chantier) · NE PAS re-implement/alter `SeaLevelMode` (reuse on render) · NE PAS skip the visual (the gap that let filaments through).

## Watchpoints global
W1 production = gallery, no reclassify · W2 morphology gates acceptance + the visual runs · W3 workflow kept OFF/experimental (instrument + eustatic scaffold) · W4 `SeaLevelMode`/`c1_default` reused on render, not re-implemented, not on reclassify · W5 v2 byte-identical + 9th-bit unaffected (full-baseline enumeration) · W6 eustatic-coast + C1-export are SEPARATE follow-ups (register, don't build) · W7 surface before code.

## Effort
~3 days: S(0.25) + E1(1) + E2(0.5) + E3(0.5) + V(0.5) + A(1) + Final(0.5). (E2 lightened — no export wiring.)
