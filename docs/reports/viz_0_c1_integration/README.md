# Issue #137 — Viz-0 C1 Integration

Wires the C1 lightweight dynamic tectonics solver (Track A/B/D foundations) into the existing ymir-viz Bevy/wgpu binary in parallel with the legacy v2 visualization. Live animation of `run_with_closures` from `bridge::c1` worker thread to a Bevy sprite, with 6 field views, overlays, and Track D evidence in egui control panels.

## CAVEAT FRANC — coast is static in the displayed path, not in C1

The visible coastline in Altitude view is **quasi-static** in Viz-0 MVP.

This is **NOT** a property of the C1 model. The C1 model's full Phase A pipeline (`run_phase_a_cycle_c1` with `Enabled` workflow) re-derives `plate_type` end-of-cycle via `apply_post_tectonic`'s reclassification, and **would** migrate the visible coast.

Viz-0 reproduces the **gallery path** (standalone `run_with_closures` — same code path as `c1_phase_2_track_d_visual_gallery.rs` which produces the reference PNGs at [`docs/reports/c1_phase_2_track_d_boundary_evolution/cycle_*.png`](../c1_phase_2_track_d_boundary_evolution/)). The gallery path does NOT call `apply_post_tectonic`. `plate_type` and `cratonic_mask` stay init-time. The visible coast lives at the static `plate_type` boundary and therefore does not migrate.

The "workflow mode" toggle that would migrate the coast is a documented Viz-0-bis follow-up (see §Backlog) and requires a Phase 1.5 multi-cycle calibration study for `macro_redistribution` before it can ship cleanly.

## What Viz-0 delivers

| Surface                            | Status                                                                                     |
|------------------------------------|--------------------------------------------------------------------------------------------|
| `bridge::c1` worker thread         | Spawned at app start, crossbeam `bounded(4)` commands / `bounded(2)` events                |
| Per-step animation                 | `StepCompleted { snapshot }` emitted per step from inside `run_with_closures`'s `on_step`  |
| 6 field views                      | S̃, Age, PlateId (12-hue hash-mod HSV), PlateType, Altitude (Architecture C), Cratonic    |
| Architecture C altitude            | `compute_isostasy` + `apply_stein_stein_bathymetry` re-applied per-snapshot in render path |
| Hypsometric bipolar palette        | `[-1.13, +1.13]` symmetric, sea_norm = 0.5 — matches gallery PNG conventions               |
| Overlays                           | Voronoi boundaries, velocity arrows (init-time only — see backlog item 1)                  |
| Track D evidence in UI             | Cumulative live stats: subduction cells, accretion merges, rifting splits, thinning cells  |
| Engine switcher                    | `ActiveEngine::{C1, V2}` resource flips `Visibility` on C1 sprite — v2 untouched           |
| Closure run-locking                | egui checkboxes disabled during `Running`, hint visible on hover                           |
| Cancel button (MVP Option C)       | Sets `AtomicBool`; effective between runs, not mid-run (Viz-0-bis hook backlog item 5)     |
| 9th bit-identical preserved        | `c1_phase_a_decomposes_into_closures_then_post_tectonic` PASS EXACT (zero ymir-core touch) |
| 4 worker tests                     | `spawns_and_runs`, `event_ordering_no_loss_under_backpressure`, `snapshot_carries_stats`, `acceptance_full_run_seed_42_pangaea_collapse` |

Production-scale measurement (seed 42 / 64² / 300 steps, default Track D closures):

```text
init num_plates       = 8
final live_plate_count = 2          ← Pangaea collapse via Track D accretion
cumulative subduction cells = 11,700
cumulative accretion merges = 6      ← matches Track D Stage A reference
cumulative rifting splits   = 0      ← expected at seed 42 (Track D Stage V evidence)
cumulative thinning cells   = 8,300
```

## A1-c detour — tried, reverted (commit `a2f5eec`)

A Stage A revision attempted to make the visible coast migrate by running the **full Phase A pipeline** in the worker (`run_with_closures` + `apply_post_tectonic` per cycle, `steps_per_cycle = 50`, 6 cycles for the default 300-step run).

**Failure mode**: `macro_redistribution` is calibrated by Phase 1.x workflow tests for **ONE** post-tectonic pass per cycle (single 300-step cycle). A1-c ran it **6× the calibrated cadence** → over-erosion → continental fraction collapsed to near-zero by step 300. PlateType view also showed rendering artefacts (white panels at later cycles) consistent with per-cycle `plate_type` mutation patterns.

**Diagnostic chain**: the failure was surfaced via user visual verification (test 4 in animation, not static PNG comparison). The diagnostic itself surfaced a learning — see the `feedback_visual_intuition_must_anchor_on_reference` memory entry.

**Reverted in this PR** to restore the standalone `run_with_closures` path matching the gallery contract.

## Architectural findings

1. **Coast migration in C1 requires reclassify.** The standalone gallery path (Viz-0) produces visually-static coast despite live S̃ dynamics. The full Phase A workflow re-derives `plate_type` via `apply_post_tectonic` end-of-cycle, and that is what migrates the visible coast. Property of the **displayed path**, not the C1 model.
2. **A1-c multi-cycle macro_redistribution failure mode.** Calibrated for 1 pass/cycle in Phase 1.x; 6× cadence over-erodes. A future workflow-mode toggle requires a Phase 1.5 calibration study for multi-cycle use.
3. **`plate_type` + `cratonic_mask` are init-only in C1.** Neither is recomputed inside `run_with_closures`. Only `apply_post_tectonic` mutates `plate_type` via reclassification; only `apply_post_tectonic` step 4 recomputes the cratonic factor (not the BoolField mask). This is by C1 design (gallery contract).

## Viz-0-bis backlog

1. **`last_plate_velocities` sibling on `C1State`** — live mid-run kinematics in snapshots. Currently snapshots carry init-time velocities only (Stage E2 W7 borrow-checker trade-off). Track D's accretion (mass-weighted merge) and rifting splits mutate kinematics mid-run; this is NOT reflected in overlays.
2. **JSON preset loader** — mirror v2 `presets.rs` for reproducible run specs sharable across runs.
3. **Track C kinematics presets** — current C1 viz wires Phase 1.1 only. Track C presets (when ready) would feed alternate plate motion patterns.
4. **Custom `IsostasyConfig` UI control** — currently hardcoded to `IsostasyConfig::default()` for gallery match.
5. **True mid-run cancel hook in `run_with_closures`** — ~1 day in ymir-core. Current MVP cancel is Option C between-runs only; current run continues to completion.
6. **Workflow mode toggle** — "Apply Phase A post-tectonic per cycle" UI checkbox. **This is the feature that would migrate the visible coast.** Requires a Phase 1.5 calibration study for `macro_redistribution` under multi-cycle use before it can ship without the A1-c-style over-erosion.

## Stage outcomes by sub-stage

| Stage | Deliverable                                                                                       |
|-------|---------------------------------------------------------------------------------------------------|
| S     | Exploration ([stage_s_exploration.md](stage_s_exploration.md)) — bridge/v2 template + workspace map |
| E1    | `C1StepStats` diagnostic field (Option B) on `C1State`. 9th bit-identical preserved.              |
| E2    | `bridge::c1` worker thread (Stage E2 W7 surfaces resolved)                                        |
| E3    | `C1Field` render + snapshot-cached view-switch during pause                                       |
| E4    | Architecture C altitude derivation + gallery palette continuity                                   |
| E5    | egui control panel + overlays + Track D live stats                                                |
| V     | 4 bridge worker tests (event ordering, snapshot stats, acceptance)                                |
| A     | Manual acceptance checklist + headless production-scale run assertion                             |
| A1-c  | Attempted Phase A workflow worker — reverted (this PR documents the detour)                       |
| Final | This README + 3 memory entries + revert commit + PR                                               |

## Files of record

- [acceptance_checklist.md](acceptance_checklist.md) — 12-section manual UI checklist + test 4 verdict + 3 findings + backlog
- [stage_s_exploration.md](stage_s_exploration.md) — Stage S workspace exploration notes
- This README — top-level Issue #137 outcome summary
