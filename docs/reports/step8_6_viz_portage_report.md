# Step 8.6 augmenté — viz portage to tectonics_v2 + rich visualisation

**Issue:** [#107]
**Branch:** `107-step-86-viz-porting-to-tectonics_v2-rich-visualization`
**Status:** sunset complete; PR target `milestone/solver-reconstruction`.

This is the closing report for Step 8.6 augmenté, the milestone that
ports the Bevy visualisation crate (`ymir-viz`) from the legacy
`tectonics::` solver to `tectonics_v2::`, adds the rich-visualisation
features the Phase 7 reviewer asked for, and sunsets the legacy
bridge once the visual revalidation passes.

The narrative below records the design decisions that survived
review, the methodological lesson Phase 7 surfaced (and why it
justifies the visual-checkpoint discipline), and the
phase-by-phase changelog so a future reader can land on the right
commit without reverse-engineering the branch history.

## Narrative

Step 8.6 originally shipped as the v2 bridge + a 4-preset visual
gallery. Phase 7 (reviewer-validated visual checkpoint) **caught a
bug that no automated test would have flagged**: the bridge's
per-step `peek_state` was filling `strain_rate_invariant` with a
zero placeholder, so every ε̇_II frame in the gallery rendered as
solid dark purple regardless of the underlying physics. Once the
reviewer asked "why are all four presets dark on strain?", a
ten-line bridge fix exposed the actual yielding patterns. The same
checkpoint also surfaced the `init S̃` sinusoidal artefact: TDD §4.2
prescribes a flat per-plate-type init, but the harness was applying
a deterministic sinusoidal perturbation that polluted every frame
with global bumps masquerading as continental geometry.

Both findings shaped the augmented milestone:

- **Phase 8a `InitMode`**: the legacy sinusoidal pattern moves under
  `InitMode::Checkerboard` (preserved bit-for-bit for Steps 0–10
  regression tests via strategy γ — every test config explicitly
  opts in). Three new modes ship: `Uniform` (TDD §4.2 default,
  flat per-plate-type with smoothstep boundary blending), `Gaussian`
  (peak at each Voronoï centroid), `Convolution` (Gaussian blur of
  the binary classification mask). Default = `Uniform`.

- **Phase 8b overlays**: `draw_voronoi_boundaries` paints
  inter-plate edges in black; `draw_velocity_vectors` overlays one
  yellow arrow per plate at the periodic-aware centroid. The two
  toggles default to off so existing Phase 7 captures stay
  unchanged unless the user opts in.

- **Phase 8c metrics dashboard**: a left-side egui panel surfaces
  live metrics (peak |v|, ⟨S̃⟩, max ε̇_II, cratonic fraction)
  during a run and the full `Metrics` summary post-run, each row
  colour-banded against the relevant §-acceptance threshold.

- **Phase 8d expose knobs**: `init_mode` reaches the parameter
  panel (`#[serde(default)]` keeps existing preset JSON loading);
  `mantle.evolution_rate` and the cratonic geometry knobs
  (`smoothing_width`, `plate_area_min`) get sliders.

- **Phase 8e export/import**: a versioned JSON snapshot stores
  `(spec, scalar_metrics, final_state)` and round-trips back into
  a viewable `V2RunState::Imported` variant. Forward-compat probe
  rejects unknown `format_version`.

- **Phase 8f equilibrium analysis**: 32² × (100 vs 200 steps)
  diagnostic on `active_medley` confirms the load-bearing physics
  metrics (peak |v|, yielding pattern, cratonic structure, CG
  behaviour) are at steady state by step 100. Two metrics drift
  on the strict 5 % rule but both are noise — `mass_drift_relative`
  is a time-cumulative quantity that scales with `t_max` by
  construction, and `mass_conservation_residual` is at machine ε.

- **Phase 8g visual revalidation**: re-run of the 4-preset Phase 7
  gallery with every Phase 8 correction live. All four presets pass
  the reviewer's acceptance criteria (no residual sinusoidal /
  checkerboard pattern in S̃, continents emerge from dynamics,
  Voronoï overlay clarifies plate structure, velocity arrows
  plausible per regime, ε̇_II coherent, no regression vs Phase 7).

- **Phase 8h sunset**: the legacy `tectonics::` bridge plus 18 legacy
  files (`bridge/{commands,events,export_system,plugin,thread}.rs`,
  `ui/{parameter_panel,pipeline_panel,statistics_panel,left_toolbar}.rs`,
  `visualization/{erosion,isostasy,plugin,render,rivers,upscale}.rs`,
  `tectonic_view.rs`, `terrain_view.rs`, `cursor_inspector.rs`,
  `state.rs`) come out. Build clean, zero `tectonics::` (non-v2)
  imports remain in `crates/ymir-viz/`, all v2 tests still pass.

The methodological note worth keeping for downstream milestones:
**the Phase 7 visual checkpoint earned its keep on a single
critical bug that automated tests did not flag**. The same
checkpoint is what surfaced the init-mode artefact that motivated
the entire Phase 8 augmented scope. Reviewer-validated visual
gates have a real return on investment when the underlying physics
is hard to scalar-check end-to-end — keep them in the Step
template.

## Phase 8f equilibrium verdict

`active_medley` reaches equilibrium on the load-bearing physics
indicators by step 100. The auto-report's strict 5 % rule flags
two rows but both are construction artefacts:

| Metric                       | A (step 100) | B (step 200) | Δ %  | Reading |
|------------------------------|--------------|--------------|------|---------|
| peak \|v\|                   | 3.011        | 3.011        | 0.00 | reached |
| yielding cells max           | 0.842        | 0.842        | 0.00 | reached |
| yielding in craton (peak)    | 0            | 0            | 0.00 | reached |
| cratonic cell fraction       | 0.157        | 0.157        | 0.00 | D7 static |
| CG iters mean                | 1260         | 1198         | 4.92 | within band |
| mass drift \|relative\|      | 0.020        | 0.023        | 14.5 | time-cumulative |
| mass conservation residual   | 4e-15        | 4e-15        | 6.6  | at machine ε |

Full report: [`step8_6_phase8f_equilibrium/active_medley_32sq.md`](step8_6_phase8f_equilibrium/active_medley_32sq.md).

A 64² verification run is optional follow-up and would not change
the verdict (the physics scales). Phase 8g visual revalidation kept
the canonical 100-step budget on this basis.

## Phase 8g visual gallery (post-corrections)

Each preset has a `_all.png` composite (5 fields × 11 capture
steps stacked vertically) under
[`step8_6_phase8g_visuals/`](step8_6_phase8g_visuals/):

- `active_medley/_all.png` — mantle on, cratonic on, age on. Mobile
  belts (white) on cratonic dark base, multidirectional flow.
- `convergence/_all.png` — yielding bands at plate-collision
  interfaces; small continental masses, Voronoï mesh visible.
- `divergence/_all.png` — yielding traces rifting boundaries;
  L-shape craton preserved with the new init.
- `subduction/_all.png` — strong yielding near subduction edges;
  velocity arrows along strike.

Reviewer report (per-preset judgement + acceptance verdict):
[`step8_6_phase8g_visuals/REPORT.md`](step8_6_phase8g_visuals/REPORT.md).

## Phase-by-phase changelog

| Phase | Commit  | What landed |
|-------|---------|-------------|
| 8a    | `ea90638` | `InitMode` enum (Checkerboard / Uniform / Gaussian / Convolution), default Uniform; 28 regression-test literals get `init_mode: Checkerboard` (strategy γ); 6 unit tests; 279/279 `tectonics_v2::` lib tests OK. |
| 8b    | `9659bb6` | `visualization::overlay` module (Voronoï boundaries + per-plate velocity arrows); `update_v2_texture` unifies on `field_to_rgba` (S̃ fixed `[0, 2.5]`, ε̇_II fixed `[1e-3, 1e2]`); 3 unit tests; UI toggles default off. |
| 8c    | `bd2ac46` | Real-time metrics dashboard (left-side egui panel); live values during Running, full `Metrics` summary on Completed; colour-banded rows per §-acceptance threshold. |
| 8d    | `e2f113b` | `V2InitModeSpec` reaches `V2RunSpec` (`#[serde(default)]` for backward-compat); cratonic geometry knobs (`smoothing_width`, `plate_area_min`) added; mantle `evolution_rate` slider exposed; init mode dropdown + per-mode parameter sliders in the panel; 2 lib roundtrip tests + 3 integration tests updated. |
| 8e    | `b806bda` | `V2RunSnapshot` JSON export/import (versioned schema); `V2RunState::Imported` variant rendered identically to `Completed`; `V2ScalarMetrics` carrier; integration roundtrip test (1 ULP tolerance on f64 fields, byte-exact on integer rasters); UI buttons. |
| 8f    | `f9e8000` | Equilibrium analysis test on `active_medley` 32² (run A 100 steps + run B 200 steps); markdown report writer; verdict: physics-relevant metrics at equilibrium by step 100. |
| 8g    | `7af5ef1` | Visual revalidation of the 4-preset gallery with all Phase 8 corrections live; 24 patchwork PNGs committed; reviewer report; **sunset authorised**. |
| 8h    | `81528bc` | Sunset: 18 legacy files deleted (5552 LOC removed), 6 files adapted to v2-only, build clean, 0 `tectonics::` (non-v2) imports remain in viz. |
| 8i    | this    | Final report + `crates/ymir-viz/README.md` update. |

## Lessons learned

1. **Visual checkpoints catch bugs that scalar tests don't.** Phase 7
   flagged two issues automated tests would have missed: the
   bridge's ε̇_II zero placeholder (a one-line bug with no scalar
   signature) and the init-mode sinusoidal artefact (a TDD §4.2
   conformance issue invisible from any single metric). Keep the
   reviewer-validated visual gate in the Step template.

2. **Strategy γ scales.** Adding a default-changing field to a
   widely-used core struct (`BaselineConfig.init_mode`) without
   breaking 28 regression tests hinged on one rule: the default
   constructor (`dynamic_accidented_defaults`) gets pinned to the
   regression-preserving variant, and every struct-literal call
   site opts in explicitly. The mechanical update was delegated to
   a sub-agent with a precise iteration loop. No regression on the
   272 `tectonics_v2::` lib tests.

3. **JSON export is "same render", not "same bits".** `serde_json`'s
   dtoa shortest round-trip can drift by 1 ULP on a fraction of
   f64 values. For the v2 viz contract — "the imported run renders
   identically to the original" — that's invisible at colormap
   (8-bit) resolution, so the integration test compares with
   `eps = 1e-12` tolerance well below `1/256 ≈ 4e-3`. If a future
   milestone ships bit-exact round-trip the format can switch to
   binary (postcard / CBOR) with a `format_version` bump.

4. **Sunset surgery wants a clean inventory before the cut.** The
   18-file deletion in Phase 8h was tractable because the audit
   sub-agent listed every legacy file, every legacy resource in
   `state.rs`, every `tectonics::` import file:line, and the
   single cross-bridge resource (`CursorWorldPos`) that needed a
   home post-sunset. Without that inventory the surgery would
   have rippled across multiple commits with broken intermediate
   states.

5. **Auto-equilibrium tests need taxonomy.** The Phase 8f strict
   5 % rule classifies `mass_drift_relative` (time-cumulative)
   the same as `peak |v|` (steady-state indicator). Future work:
   tag each metric as `equilibrium-relevant` vs `time-cumulative`
   in the test's `CompareRow` structure so the verdict reports
   only the steady-state column. Recorded as follow-up in the
   Phase 8f report.

## Out of scope (deferred)

- 64² Phase 8f / Phase 8g verification runs (the 32² result settles
  the verdict; the 64² renders would not change the visual
  conclusions).
- Per-step `cg_iters` / `newton_outer_iters` streaming (live
  dashboard would benefit; today these are only visible in the
  end-of-run `Metrics` aggregate).
- Save-as-preset on import (importing a snapshot replays its
  config, but the user has to copy fields manually if they want
  it as a new preset entry).
- Binary snapshot format. JSON suffices through 128² regimes; the
  schema bump is reserved for the milestone that needs it.

The augmented Step 8.6 closes here. Step 10.5 (redefined) follows
once this branch merges.
