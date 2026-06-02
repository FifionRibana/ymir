# Issue #139 — Viz-0.5: hover-to-inspect + workflow mode + hypsometric lens + continuation-ready

Builds on Viz-0 (Issue #137). Adds per-cell inspection, the calibrated Phase A **workflow** pipeline alongside the gallery path, a hypsometric meters lens, and continuation-ready worker structure. **Viz-only** — `tectonics_c1` untouched (9th bit-identical preserved).

## What ships

| Surface | Detail |
|---|---|
| **Hover-to-inspect** | Bottom-right "C1 Cell Inspector": `(i,j)`, S̃, age, altitude. Reads `CursorWorldPos` (world-space, camera-corrected), inverts the C1 sprite transform (`world_to_cell`, Y-flip-correct, resize-correct). Altitude derived lazily (`derive_altitude_field`, cached by snapshot step) → works in **all** views. |
| **Non-dim first, meters second** | Non-dim altitude is the verification value (shown first); meters are the cosmetic hypsometric lens (shown second, "after hypsometric curve"). |
| **Workflow pipeline** | `ActivePipeline::Workflow` → `C1Command::RunWorkflow`. The calibrated Phase A loop: `n_cycles` cycles of `run_with_closures(k_cycle)` + per-cycle `apply_post_tectonic` (sea-level → macro-redistribution → reclassify). Coast migrates per cycle. Gallery `RunBaseline` path left untouched. |
| **Calibrated cadence** | `PhaseAParams::default()` = `n_cycles 5 × k_cycle 20` = 100 tectonic steps, 5 macro passes. Total = `n_cycles × k_cycle`, **never** `n_steps` (the A1-c over-erosion guard). |
| **Hypsometric lens** | `HypsometricScale` resource: land scale/gamma tunable, ocean anchored to Stein-Stein `depth_scale_m = 5000`. Drives ONLY the hover meters — never the non-dim value or the model. |
| **Continuation-ready** | Worker retains the last completed run's `(C1State, PlateKinematics)`. `run_workflow_cycles` is the reusable continuation core (accepts `base_step`). A future `ContinueRun` would resume from the retained state. Capability only — no command/UI yet. |
| **9th bit-identical** | Preserved (viz-only; zero `tectonics_c1` changes). |

## Hover: non-dim = verification, meters = cosmetic

The non-dimensional altitude is the value to verify against — it comes straight from the Architecture C derivation (`compute_isostasy` + `apply_stein_stein_bathymetry`). The meters readout is a presentation lens (`hypsometric_meters`): land `= land_scale·alt^γ`, ocean `= −(ocean_scale·|alt|^γ)`. Tuning the lens never touches the non-dim value or the simulation. This separation is load-bearing — meters are for human intuition, non-dim is for validation.

## Workflow mode = the A1-c fix, calibrated

Viz-0 (Issue #137) reverted the A1-c worker (which ran `apply_post_tectonic` 6× the calibrated cadence — 50/300 — and over-eroded the continent) and shipped the static-coast gallery path, documenting that a *calibrated* workflow was a Viz-0-bis follow-up. Viz-0.5 ships it: the workflow worker reproduces `run_phase_a_cycle_c1(Enabled)` at the calibrated `PhaseAParams::default()` cadence, inlining `run_with_closures` only for the per-step animation hook. The cadence is `n_cycles × k_cycle`, with a UI calibration warning when the user moves off the default (raising the cadence without lowering `alpha` over-erodes — the A1-c mechanism).

## The continental-fraction finding (Stage V diagnostic)

The calibrated workflow converges to a continental/emergent fraction of **~0.058** at seed 42 (64²), not the ~0.20–0.45 originally estimated. The Stage V diagnostic established this is **correct**, on evidence:

1. **The ~0.27 init "continent" is a geometric LABEL, not emergent land.** Init `plate_type` = Voronoï continental seeds (0.27). The above-sea-level land (`compute_isostasy(s).land_ratio`, S̃-only) is ~0.045 even in the pure-tectonic gallery. The [0.20,0.45] target compared the geometric label to a final emergent fraction — apples/oranges.
2. **Not macro erosion.** `macro_redistribution` is mass-conserving (`Δmass ≈ 0`/cycle, `isostatic_rebound_ratio = 0.80`). The continental loss is **declassification**, not mass transport.
3. **Adaptive `sea_level_ref` drift is the mechanism.** `sea_level_ref = s_min + 0.4·(s_max − s_min)` jumps 0.52 → 0.98 in cycle 1 as `s_max` doubles (1.0 → 2.18) via Davis-Suppe orographic peaks — declassifying 753 cells by threshold rise. The workflow converges (last-cycle Δ 0.0007) and sits ABOVE the gallery isostatic floor (0.0588 > 0.0449).
4. **Pre-existing C1 debt, revealed not created.** The `min/max`-based `sea_level_ref` is fragile to peak outliers. A robust **percentile** sea level is a Phase 1.5 follow-up (core change, out of viz-only scope). This also explains "little emergent land at 256²" from the Viz-0 thread: C1 intrinsically produces little emergent land at seed 42 under the current sea-level formula — model behavior, not a viz bug. **For Phase 3 validation**: "little emergent land" is the normal seed-42 state.

## The A1-c regression guard (qualitative signature)

`workflow_mode_continent_preserved` does not assert a magic band; it asserts the signature that distinguishes the calibrated-converged workflow from the A1-c runaway:

1. mass-conserving — `|Δmass|/mass < 1e-3` per cycle;
2. converges — last-cycle `|Δfrac| < 0.01`;
3. above the isostatic floor — workflow `iso_land ≥ 0.9 × gallery iso_land` at the same 100 steps (apples-to-apples on the S̃-only emergent land, NOT the plate_type-gated rendered land).

## Stage outcomes

| Stage | Deliverable | Commit |
|---|---|---|
| S | Exploration (cursor→cell, split-borrow, calibration anchor, altitude cache) | `6ba2a73` |
| E1 | Hover-to-inspect per-cell readout | `a2c1a46` |
| E2 | Workflow-mode worker (calibrated PhaseAParams, per-cycle post-tectonic) | `7127e72` |
| E3 | Pipeline toggle + cadence sliders + hypsometric scale | `a0b6065` |
| E4 | Continuation-ready worker state retention (capability only) | `58e3902` |
| V | Workflow worker tests + A1-c regression guard | `097a35b` |
| A | Acceptance (migrating-coast product test + manual checklist) | `2b21cb6` |
| Final | This README + memory + PR | — |

## Tests

8 `bridge::c1::thread` tests PASS (4 gallery + 4 workflow) + 1 `#[ignore]` diagnostic; 6 `c1_plugin` tests (4 `world_to_cell` incl. Y-flip + resize, 2 `hypsometric_meters`); 7 `c1_viz` render tests (behavior-preserving after `derive_altitude_field` extraction). 9th bit-identical PASS (viz-only).

## Viz-0-bis / Phase 1.5 follow-ups

- **Percentile `sea_level_ref`** (Phase 1.5, core): robust to Davis-Suppe `s_max` outliers; would raise the emergent-land fraction. The reason emergent land is so low at seed 42.
- **`ContinueRun` command + UI**: resume from the retained state via `run_workflow_cycles` (structure already in place).
- Carried from Viz-0: `last_plate_velocities` (live mid-run kinematics), JSON presets, Track C kinematics, custom `IsostasyConfig`, true mid-run cancel.

## Files

- `crates/ymir-viz/src/bridge/c1/{commands,events,mod,plugin,thread}.rs`
- `crates/ymir-viz/src/visualization/{c1_plugin,c1_viz}.rs`
- `docs/reports/viz_0_5/{README,stage_s_exploration,acceptance_checklist}.md`
