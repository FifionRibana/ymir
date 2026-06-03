# Issue #141 — Phase 1.5 Stage S audit

Robust (P95-capped) sea-level definition for C1 land/sea classification. Core change (shared `tectonics` + `tectonics_v2/workflow`); **v2 byte-identical** is the hard guard. Branch `141-phase-15-robust-p95-capped-sea-level-definition-for-c1-landsea-classification` off `milestone/c1-lightweight-dynamic-tectonics` (head `659079f` — **#139/Viz-0.5 is merged in**, PR #140).

This audit **corrects the issue's predicted site list** and maps the real architecture before any code.

## Finding 1 — the issue's named sites are ALL test sites (the directory/line trap)

The issue predicted production sites `time_loop.rs:811,839` + `phase_a_c1.rs:192`. Confirmed against the code:

- `time_loop.rs:811` → inside `#[cfg(test)] mod tests`, test `uniform_field_advection_conserves_mass_exactly`.
- `time_loop.rs:839` → test `callback_fires_once_per_step`.
- `phase_a_c1.rs:192` → test `c1_phase_a_cycle_completes_300_steps` (`let iso_config = IsostasyConfig::default();`).

**None are production.** The issue's `IsostasyConfig::default()` grep caught test fixtures. The directory trap is real but inverted: the concern isn't "is `phase_a_c1.rs` C1?" (it is) — it's that the named lines aren't where the sea level is actually decided.

## Finding 2 — C1 sea level is CALLER-THREADED, not hardcoded

The C1 production path reads the sea level from a passed config, in two functions (the real change points):

- **`compute_isostasy` h_sea** ([isostasy.rs:88](../../../crates/ymir-core/src/tectonics/isostasy.rs#L88)): `h_sea = h_min + config.sea_level_fraction · h_range` (h = S̃·buoyancy). Drives the altitude heightmap + `land_ratio`.
- **`compute_sea_level_ref_s_space`** ([phase_a_common.rs:271](../../../crates/ymir-core/src/tectonics_v2/workflow/phase_a_common.rs#L271)): `s_min + iso_cfg.sea_level_fraction · (s_max − s_min)`. Drives reclassify + macro-redistribution + (per-step) drainage.

C1 production consumers, all reading a threaded config:

- `time_loop.rs:592` — `compute_isostasy(&state.s, &config.iso_config)` (per-step altitude).
- `time_loop.rs:614` — `compute_sea_level_ref_s_space(&state.s, &config.iso_config)` (per-step drainage sea level → erosion). **Runs PER STEP.**
- `apply_post_tectonic` — `compute_sea_level_ref_s_space(input.s_field, input.iso_cfg)` (per-cycle reclassify + macro).

`config.iso_config` is `C1TimeLoopConfig.iso_config`; `input.iso_cfg` is `PhaseACycleInputC1.iso_config` (`run_phase_a_cycle_c1` threads the caller's config into both — pass-through, not a construction site). **So the E1 design is valid** (add `sea_level_mode` to `IsostasyConfig`; branch both functions), but **E2's "wire 3 core sites" is wrong** — the wiring is at whoever *builds* the C1 config.

## Finding 3 — the REAL C1 production wiring sites are in ymir-viz (#139 merged)

Since #139 is in the milestone, the C1 configs are built in the viz worker:

- `thread.rs:1057` — gallery `RunBaseline` `C1TimeLoopConfig { iso_config: IsostasyConfig::default() }` → **C1, → c1_default**.
- `thread.rs:1228` — workflow `run_workflow_cycles` `let iso_config = IsostasyConfig::default();`, threaded to both the time loop (`:1236`) AND `apply_post_tectonic` (`:1267`) → **C1, → c1_default** (covers drainage + reclassify with one change).
- `c1_viz.rs:270` — `derive_altitude_field`: `compute_isostasy(&s_field, &IsostasyConfig::default())` for the rendered Altitude view → **C1, → c1_default**. **COHERENCE-CRITICAL (W2)**: if the render uses MinMaxFraction while reclassify uses P95-cap, the altitude=0 coast ≠ the plate_type coast — exactly the Viz-0 Stage A divergence. Both must use the same mode.

Plus ymir-core C1 **tests / gallery generators / acceptance** that build `C1TimeLoopConfig`/`PhaseACycleInputC1` (re-validation surface, Stage Cal).

## Finding 4 — v2 + export sites: UNTOUCHED (the byte-identical guard)

Keep `Default` (= `MinMaxFraction`) at every v2/export/legacy site:

- `phase_a_v2.rs:123/131` (v2 workflow), `phase_b.rs:78` (HD finalization, paradigm-agnostic but v2-default).
- `v2_viz.rs:715/956`, `phases/isostasy.rs:36/187` (viz v2 + standalone isostasy phase).
- `export/*`, and the `isostasy.rs` unit tests (146/165/180/197/217/235/237/257/258) — these pin MinMaxFraction behavior.

W1 verification: a diff of v2 test outputs must be unchanged.

## Finding 5 — serde (W4 satisfied)

`IsostasyConfig` is `#[derive(Debug, Clone, Serialize, Deserialize)]`. `SeaLevelMode` gets the same derives; the new field carries `#[serde(default)]` so legacy configs (no `sea_level_mode`) deserialize to `MinMaxFraction` — avoids the Issue #47-style missing-field break.

## Finding 6 — percentile cost + the Q4 convergence risk

`<[f64]>::select_nth_unstable_by(|a,b| a.partial_cmp(b).unwrap())` — O(N) average, std-stable. Needs a mutable copy of the field per call (~32 KB at 64²). The drainage consumer (`time_loop.rs:614`) runs **per step**, so P95 is recomputed per step. At 64² that's ~4 K ops + one clone per step — negligible runtime. **But** a per-step percentile threshold could make the drainage sink set jitter step-to-step → the **Q4 convergence risk** (Stage V gate). If it oscillates, the sub-fix is a per-cycle (stable-within-cycle) drainage threshold.

## Finding 7 — 9th bit-identical: 3 test files

`c1_phase_1_3_workflow.rs` (`c1_phase_a_decomposes_into_closures_then_post_tectonic`), `c1_phase_2_track_d_acceptance.rs` (`ninth_bit_identical_preservation_phase_2_r7`), `c1_phase_2_boundary_evolution.rs`. The decomposition CONTRACT (wrapper == steps) must stay exact under P95-cap; only reference VALUES change (Stage Final protocol).

## Corrected plan delta (vs the issue)

- **E2 is not "3 core sites → c1_default".** It's: 3 viz production sites (`thread.rs:1057`, `thread.rs:1228`, `c1_viz.rs:270`) + the C1 core test/gallery callers for re-validation. The core formula change is E1 (two functions + the enum).
- **Scope decision needed (surfaced below)**: how far to propagate `c1_default` into the ymir-core C1 **gallery generators** (Track A/B/D `*_visual_gallery.rs`, `#[ignore]`) and acceptance tests — i.e., do the committed-convention gallery PNGs + acceptance numbers move to P95-cap, or stay MinMaxFraction as a reference?

## W7 surfaces resolved

- Two formula instances confirmed (isostasy.rs:88 h-space, phase_a_common.rs:271 S̃-space); both read passed config → branch on `sea_level_mode`.
- C1 sites = caller-threaded; real production wiring is in viz (#139 merged), incl. the coherence-critical render altitude.
- v2/export sites enumerated for the untouched guard.
- serde default feasible; percentile O(N) per-step with a Q4 convergence risk; 9th bit-identical = 3 files, contract-exact / values-redefined.
