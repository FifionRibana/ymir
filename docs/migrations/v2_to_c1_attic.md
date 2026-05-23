# Issue #117 Phase 1 — Pre-refactor audit for v2 → `_attic` migration

**Status:** Phase 1 audit only. Read-only. No file moves, no source modifications, no `cfg` gates added.

**Scope:** Inventory of couplings between the modules retired per §4.7 of `docs/design/c1_lightweight_dynamic_tectonics.md` and the rest of the codebase, so Phase 2 (file moves) and Phase 3 (`cfg` gating + regression verification) execute without rediscovery.

**Audit cut-off commit:** `0969ee4` (Step 12 R7 Section 7), branch `112-step-12-interleaved-tectonic-erosion-workflow` post-merge into `milestone/solver-reconstruction`.

---

## Summary of findings

- **Six modules listed in §4.7** retire as planned. No surprise on those.
- **Hidden coupling 1 — `tectonics_v2/slab/` is NOT in §4.7** but is the partner of `forcing/slab_pull.rs` (which IS retired). They are inseparable. **Recommendation: add `tectonics_v2/slab/` to §4.7.**
- **Hidden coupling 2 — `tectonics_v2/presets.rs`** is in neither §4.7 nor §4.8. Contains `RheologyParams`, `YieldingConfig`, `Preset` — rheology-typed. **Recommendation: retire with rheology.**
- **Major surface to gate:** `diagnostics/harness.rs` is the orchestrator that ties stokes + mantle + slab + basal_drag + rheology + forcing into the v2 run path. It cannot be partially gated cleanly — recommend gating the whole file. Same for `diagnostics/{mms_bench, k_sub_sweep, num_plates_sweep, bi_sweep, br_sweep, ar_sweep}.rs`.
- **Bridge v2 build path** (`crates/ymir-viz/src/bridge/v2/build_config.rs`) is fully Stokes-coupled and must be feature-gated. The V2 spec layer (`spec.rs`) stays compilable as pure data.
- **Forcing module** stays in `tectonics_v2/forcing/` but the two retired implementations (`MantleForce`, `SlabPullForce`) and their re-exports in `forcing/mod.rs` need gating.
- **~30 integration tests + 6 binaries + 1 bench** need gating or relocation to `_attic`.
- **All 8 V2 JSON presets** contain retired-mapped fields. They are data (no code), so they stay in place; only the *consumers* gate.

Total file-touch count for Phase 3 Commit 2 (gating wiring), excluding `git mv` targets : **~ 50 files** (production: 12, tests: 32, bins+benches: 7, recommend audit per file before commit).

**Recommendation: pause before Phase 2 for user decision** on the two hidden couplings (`slab/`, `presets.rs`) plus the harness-gating strategy. Once those are resolved, Phase 2 + 3 should proceed cleanly.

---

## D1 — Direct imports (per retired module)

### D1.a — `tectonics_v2::stokes::*`

**Production code (gateable, not movable):**

| File | Lines |
|---|---|
| [`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) | 47–55 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/mms_bench.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/mms_bench.rs) | 19–24 |
| [`crates/ymir-viz/src/bridge/v2/build_config.rs`](../../crates/ymir-viz/src/bridge/v2/build_config.rs) | 18–20 |

Internal to retired subtree (move with parent):

| File | Lines |
|---|---|
| [`crates/ymir-core/src/tectonics_v2/stokes/sparse_assembly.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/sparse_assembly.rs) | 294 |
| [`crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/nonlinear_solver.rs) | 723 |
| [`crates/ymir-core/src/tectonics_v2/stokes/amg/mod.rs`](../../crates/ymir-core/src/tectonics_v2/stokes/amg/mod.rs) | 237–238 |

**Binaries:**

- [`crates/ymir-core/src/bin/gen_reference_solutions.rs`](../../crates/ymir-core/src/bin/gen_reference_solutions.rs) (lines 37–47)
- [`crates/ymir-core/src/bin/gen_bench_data.rs`](../../crates/ymir-core/src/bin/gen_bench_data.rs) (lines 44–45)

**Benches:**

- [`crates/ymir-core/benches/amg_benchmark.rs`](../../crates/ymir-core/benches/amg_benchmark.rs) (lines 28–38)

**Tests** (full list, both crates):

`v2_amg_physics_scalar_parity`, `v2_amg_phase3_diagnostic`, `v2_amg_scalar_parity`, `v2_amg_poisson_projection_diag`, `v2_basal_drag_mms`, `v2_baseline_32sq`, `v2_basal_drag_operator_symmetric`, `v2_continuation`, `v2_yielding_newton`, `v2_yielding_mms`, `v2_precond_drag_diagonal`, `v2_stokes_mms_variable_eta`, `v2_plate_kinematic_scenarios`, `v2_stokes_mms`, `v2_plate_kinematic`, `v2_step9_smoothing_probe`, `v2_picard_parity`, `v2_step9_physics_and_sweep`, `v2_step10_physics_and_regression`, `v2_parallel_determinism`, `v2_step0_synthetic_parity`, `v2_nullspace`, `v2_sparse_assembly_snapshot_parity`, `v2_newton_extrapolation`, `v2_step13_cg_ratio`, `v2_newton_convergence`, `v2_step13_5_cg_ratio`, `v2_slab_null_space_preservation`, `v2_mantle_relaxation`, `v2_mantle_null_space_preservation`, `v2_workflow_r4_visual_checkpoint` (lines 1596, 1603–1604).

### D1.b — `tectonics_v2::forcing::{slab_pull, mantle_force}::*`

**Direct module-path imports:** none found. The retired forcing implementations are consumed via the parent `forcing::*` re-export (`MantleForce`, `SlabPullForce`). See D1.g and the forcing/mod.rs note in §"Special cases".

### D1.c — `tectonics_v2::mantle::*`

**Production code:**

| File | Lines |
|---|---|
| [`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) | 36 |
| [`crates/ymir-viz/src/bridge/v2/build_config.rs`](../../crates/ymir-viz/src/bridge/v2/build_config.rs) | 12 |

**Binaries:** `step8_baseline.rs:40`, `gen_bench_data.rs:38`.

**Tests:** `v2_amg_physics_scalar_parity:35`, `v2_baseline_32sq:28`, `v2_mantle_force_mms:26`, `v2_mantle_evolution_rate:27–28`, `v2_mantle_div_free:14,17`, `v2_mantle_runaway_diagnostic:18`, `v2_mantle_relaxation:28–29`, `v2_mantle_null_space_preservation:21`, `v2_newton_extrapolation:24,163`, `v2_step10_physics_and_regression:31`, `v2_step0_synthetic_parity:23`, `v2_step9_smoothing_probe:25`, `v2_step9_physics_and_sweep:37`, `v2_step8_regression_smoke:29`, `v2_workflow_phase_b:30`, `v2_workflow_phase_a_loop:27`, `v2_step13_cg_ratio:41`, `v2_workflow_phase_a_cycle:31`, `v2_step13_5_cg_ratio:41`, `v2_workflow_phase_6_visual_checkpoint:50`, `v2_workflow_r4_visual_checkpoint:1448` (test-mod scope).

### D1.d — `tectonics_v2::rheology::*`

**Production code:**

| File | Lines |
|---|---|
| [`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) | 41 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/mms_bench.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/mms_bench.rs) | 18 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/k_sub_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/k_sub_sweep.rs) | 24 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/num_plates_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/num_plates_sweep.rs) | 24 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/bi_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/bi_sweep.rs) | 14 |
| [`crates/ymir-viz/src/bridge/v2/build_config.rs`](../../crates/ymir-viz/src/bridge/v2/build_config.rs) | 15 |

**Binaries:** `step5_baseline.rs:39`, `step6_baseline.rs:39`, `step7_baseline.rs:38`, `step8_baseline.rs:45`, `step_baseline.rs:35`, `gen_bench_data.rs:41`.

**Tests** (count: 23): `v2_plate_kinematic_scenarios`, `v2_plate_kinematic`, `v2_continuation`, `v2_picard_parity`, `v2_closed_recycling_conservation`, `v2_boundary_regression_smoke`, `v2_newton_extrapolation`, `v2_baseline_32sq`, `v2_newton_convergence`, `v2_step13_cg_ratio`, `v2_mass_balance_residual`, `v2_mantle_runaway_diagnostic`, `v2_step13_5_cg_ratio`, `v2_step10_physics_and_regression`, `v2_yielding_newton`, `v2_amg_physics_scalar_parity`, `v2_yielding_mms`, `v2_viscosity`, `v2_step9_physics_and_sweep`, `v2_step8_regression_smoke`, `v2_step9_smoothing_probe`, `v2_step7_regression_smoke`, `v2_step6_refactor_parity`.

### D1.e — `tectonics_v2::basal_drag::*`

**Production code:**

| File | Lines |
|---|---|
| [`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) | 20 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/num_plates_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/num_plates_sweep.rs) | 19 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/k_sub_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/k_sub_sweep.rs) | 20 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/br_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/br_sweep.rs) | 26 |
| [`crates/ymir-core/src/tectonics_v2/diagnostics/bi_sweep.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/bi_sweep.rs) | 11 |
| [`crates/ymir-viz/src/bridge/v2/build_config.rs`](../../crates/ymir-viz/src/bridge/v2/build_config.rs) | 5 |

**Binaries:** `step5_baseline.rs:24`, `step6_baseline.rs:24`, `step7_baseline.rs:26`, `step8_baseline.rs:30`, `step_baseline.rs:33`, `gen_bench_data.rs:32`.

**Tests** (count: 25): `v2_workflow_phase_b`, `v2_workflow_phase_a_loop`, `v2_workflow_phase_a_cycle`, `v2_workflow_phase_6_visual_checkpoint`, `v2_step9_smoothing_probe`, `v2_step9_physics_and_sweep`, `v2_step8_regression_smoke`, `v2_step7_regression_smoke`, `v2_step6_refactor_parity`, `v2_step13_cg_ratio`, `v2_step13_5_cg_ratio`, `v2_step10_physics_and_regression`, `v2_step0_synthetic_parity`, `v2_precond_drag_diagonal`, `v2_newton_extrapolation`, `v2_mass_balance_residual`, `v2_mantle_runaway_diagnostic`, `v2_closed_recycling_conservation`, `v2_boundary_regression_smoke`, `v2_baseline_32sq`, `v2_basal_drag_unit`, `v2_basal_drag_regression_smoke`, `v2_basal_drag_operator_symmetric`, `v2_basal_drag_mms`, `v2_amg_physics_scalar_parity`.

### D1.f — `tectonics_v2::slab::*` (hidden coupling, NOT in §4.7)

This is the slab-pull infrastructure that pairs with the retired `forcing/slab_pull.rs`. `forcing/slab_pull.rs` documents that it wraps `super::super::slab::state::SlabState` and `super::super::slab::convergence_direction`. **The two modules are inseparable.** §4.7 mentions `forcing/slab_pull.rs` but not `tectonics_v2/slab/`. Recommend treating them together.

**Production:** `diagnostics/harness.rs:43`, `bridge/v2/build_config.rs:17`, `step7_baseline.rs:40`, `step8_baseline.rs:47`, `gen_bench_data.rs:43`.

**Tests** (count: 17): `v2_amg_physics_scalar_parity:40`, `v2_baseline_32sq:35`, `v2_mantle_runaway_diagnostic:23`, `v2_newton_extrapolation:29`, `v2_slab_decay_mms:17`, `v2_slab_advection_mms:18`, `v2_workflow_phase_b:35`, `v2_workflow_phase_a_loop:32`, `v2_workflow_phase_a_cycle:36`, `v2_workflow_phase_6_visual_checkpoint:55`, `v2_step9_physics_and_sweep:44`, `v2_step13_cg_ratio:49`, `v2_step0_synthetic_parity:26`, `v2_step8_regression_smoke:34`, `v2_step9_smoothing_probe:30`, `v2_step13_5_cg_ratio:49`, `v2_step7_regression_smoke:29`, `v2_step10_physics_and_regression:38`.

### D1.g — `tectonics_v2::forcing::*` re-exports

`forcing/mod.rs` re-exports `MantleForce` and `SlabPullForce` via `pub use mantle_force::MantleForce;` (line 20) and `pub use slab_pull::SlabPullForce;` (line 22). The `pub mod mantle_force;` and `pub mod slab_pull;` declarations (lines 13, 15) make the modules visible. These four lines need gating.

**Consumers of `MantleForce` / `SlabPullForce` via the parent re-export:**

- [`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs:32`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) — `use crate::tectonics_v2::forcing::{MantleForce, SlabPullForce};`
- Tests: `v2_mantle_force_mms`, `v2_mantle_null_space_preservation`, `v2_mantle_relaxation`, `v2_slab_force_mms`, `v2_slab_null_space_preservation`.

`GpeForce`, `SinusoidalForce`, `ZeroForce`, `BodyForce`, `SimulationState`, `VectorField`, `ForceSum` from the same module **stay** (Gpe is the dome-flow term, preserved; the rest are abstractions).

---

## D2 — Indirect call-sites in preserved modules

### D2.a — `diagnostics::harness::run_baseline_with_progress`

[`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs:605`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs#L605)

The full body of `run_baseline_with_progress` orchestrates Newton + CG + Picard + continuation + mantle + slab + basal_drag + rheology + forcing. It cannot be partially gated cleanly because every internal call touches retired modules. The `BaselineConfig` struct (line 135) also carries fields typed by retired types (`mantle: MantleConfig`, `slab_pull: SlabPullConfig`, `basal_drag: BasalDragConfig`, `yielding: YieldingConfig`).

**Recommended gating:** the entire file `diagnostics/harness.rs` with `#![cfg(feature = "v2_legacy")]`. Consumers (bridge, tests) gate accordingly.

### D2.b — `diagnostics/mms_bench.rs`, `k_sub_sweep.rs`, `num_plates_sweep.rs`, `bi_sweep.rs`, `br_sweep.rs`, `ar_sweep.rs`

All import retired types directly (stokes, rheology, basal_drag, forcing). Same gating strategy: entire file `#![cfg(feature = "v2_legacy")]`.

`ar_sweep.rs` only imports `forcing::{ForceSum, GpeForce}` (preserved) plus the `diagnostics::harness` glue — but the glue (BaselineConfig) is retired. So `ar_sweep` indirectly depends on retired surface and must gate.

### D2.c — `workflow/phase_a.rs`

(Not in the D1 grep results since it doesn't import the retired modules directly.) The Phase A orchestrator calls `harness::run_baseline_with_progress` to execute the tectonic sub-cycle. If `harness.rs` gates, `workflow/phase_a.rs` either:

- gates with `#![cfg(feature = "v2_legacy")]` (loses Phase A entirely without feature), OR
- splits into a paradigm-agnostic shell + v2-specific implementation. C1 will need its own `run_baseline` analogue; the Phase A shell could be made tectonic-impl-agnostic.

**Recommendation:** for the migration commit, gate `workflow/phase_a.rs` entirely. The shell refactor is a C1 Phase 1.3 concern (per §7 C1.md), not this migration.

### D2.d — `bridge/v2/build_config.rs::build`

Function `build(spec: &V2RunSpec) -> BaselineConfig` translates the V2 spec into a v2 harness config. Every branch of the function uses a retired type. The V2RunSpec struct itself stays (data-only); only the `build()` call cannot run without `v2_legacy`. Gate the whole file with `#![cfg(feature = "v2_legacy")]`.

### D2.e — `bridge/v2/thread.rs` `spawn_v2_thread` and `V2Command::RunBaseline`

(Not grepped explicitly above but called from `build_config::build()` and downstream.) The thread harness wires the v2 solver loop. **All gating-dependent.** `spawn_v2_thread`, `V2Command::RunBaseline`, `V2Event::Progress { peek_state }`, `V2Event::Completed { metrics }` all touch the v2 surface and must gate.

The bridge's `V2BridgePlugin`, `V2SolverBridge` (Bevy resource), and the UI hooks consuming `V2RunState` need to either gate (if the UI runs only with `v2_legacy`) or get a paradigm-agnostic shim. Recommend: **gate the bridge layer entirely; C1 will introduce its own bridge.**

---

## D3 — Type consumers in the spec layer

| Type (in `crates/ymir-viz/src/bridge/v2/spec.rs`) | Status | Consumers |
|---|---|---|
| `V2MantleSpec` (line 23) | retired-mapped | `parameter_panel_v2.rs:253-271`, `build_config.rs:72-79`, all R4/R5b/R6/R7 tests (33+ matches) |
| `V2ForceKind` (line 133) | mixed — `Gpe` preserved, `Sinusoidal` is test-only | `parameter_panel_v2.rs:519-535`, `build_config.rs:66-67`, bridge tests |
| `slab_enabled: bool` (line 727) | preserved (data) | `build_config.rs:142-150` (currently hard-codes `SlabPullConfig::Disabled` regardless) |
| `V2LinearSolverSpec` | retired-mapped | spec.rs + build_config.rs |
| `V2YieldingLaw` | (not found — yielding is encoded inline in `V2RunSpec` fields `bi`, `br`) | spec.rs (data) + build_config.rs:132-135 |
| `V2CratonicSpec` | **preserved** (cratonic stays per §4.8; `b_factor` smoothstep amplification retires) | spec.rs + UI + build_config.rs |
| `V2AgeFieldSpec` | preserved | spec.rs + UI |
| `V2InitModeSpec` | preserved | spec.rs + UI |
| `V2PlateKinematicSpec` | preserved | spec.rs + UI |

**Gating strategy:**

- All `V2*Spec` enum *declarations* stay in `spec.rs` (data, no v2 imports needed at type level — they're plain serde structs).
- The `into_core()` methods that map V2 spec types to retired core types (`MantleConfig`, `LinearSolverConfig`, etc.) must gate with `#[cfg(feature = "v2_legacy")]`.
- `build_config.rs` calls those `into_core()` methods, so its `build()` function gates implicitly.

This keeps the JSON round-trip working (data loads fine without feature), but the spec only "builds" into a runnable config under `v2_legacy`.

---

## D4 — UI consumers

`crates/ymir-viz/src/ui/parameter_panel_v2.rs`:

| Lines | Section | Retired dependency |
|---|---|---|
| 226–248 | "Yielding & drag" CollapsingHeader | `V2RunSpec.bi`, `V2RunSpec.br` (yielding-rheology fields) + drag sliders |
| 249–298 | "Mantle" CollapsingHeader | `V2MantleSpec::On { mf, coupling, num_modes, seed, evolution_rate }` |
| 519–535 | Force toggle (Gpe / Sinusoidal) | `V2ForceKind` (Sinusoidal variant is test-mode only) |
| 17–18 | Imports | `V2ForceKind, V2InitModeSpec, V2LinearSolverSpec, V2MantleSpec` |

**Gating strategy:** the entire `parameter_panel_v2.rs` does not import `tectonics_v2::*` directly; it only sees `V2*Spec` enum variants. If the enum variants stay declared (D3), the UI compiles without `v2_legacy`. **But** the resulting `V2RunSpec` cannot be `build()`-ed without `v2_legacy` — i.e. running the simulation from UI requires the feature.

Decision needed: do we want the UI to render the v2 panel without the feature (greyed-out "Run" button), or hide the entire v2 panel when feature is off? Recommend: **gate the UI panel itself** (entire `parameter_panel_v2.rs` with `#![cfg(feature = "v2_legacy")]` or its mount point in the panels manager), and prepare a C1 panel for the default build.

---

## D5 — Test consumers (full inventory)

### D5.a — Tests gated by retired-module import

From D1.a–D1.f, **~ 36 unique test files** depend on at least one retired module. Listing them by primary coupling for clarity (a single test typically imports several retired modules):

**`crates/ymir-core/tests/`:**

- v2_amg_phase3_diagnostic, v2_amg_physics_scalar_parity, v2_amg_poisson_projection_diag, v2_amg_scalar_parity
- v2_basal_drag_mms, v2_basal_drag_operator_symmetric, v2_basal_drag_regression_smoke, v2_basal_drag_unit
- v2_baseline_32sq, v2_boundary_regression_smoke, v2_closed_recycling_conservation
- v2_continuation
- v2_mantle_div_free, v2_mantle_evolution_rate, v2_mantle_force_mms, v2_mantle_null_space_preservation, v2_mantle_relaxation, v2_mantle_runaway_diagnostic
- v2_mass_balance_residual
- v2_newton_convergence, v2_newton_extrapolation, v2_nullspace
- v2_parallel_determinism, v2_picard_parity, v2_plate_kinematic, v2_plate_kinematic_scenarios, v2_precond_drag_diagonal
- v2_slab_advection_mms, v2_slab_decay_mms, v2_slab_force_mms, v2_slab_null_space_preservation
- v2_sparse_assembly_snapshot_parity
- v2_step0_synthetic_parity, v2_step6_refactor_parity, v2_step7_regression_smoke, v2_step8_regression_smoke
- v2_step9_physics_and_sweep, v2_step9_smoothing_probe, v2_step10_physics_and_regression
- v2_step13_cg_ratio, v2_step13_5_cg_ratio
- v2_stokes_mms, v2_stokes_mms_variable_eta
- v2_viscosity, v2_workflow_phase_a_cycle, v2_workflow_phase_a_loop, v2_workflow_phase_b, v2_workflow_phase_6_visual_checkpoint
- v2_yielding_mms, v2_yielding_newton, v2_yielding_regression_smoke

**`crates/ymir-viz/tests/`:**

- v2_bridge_export_import_roundtrip (imports `V2ForceKind`, `V2MantleSpec`)
- v2_bridge_field_extraction (imports `V2ForceKind`, `V2MantleSpec`)
- v2_bridge_lifecycle (imports `V2ForceKind`, `V2MantleSpec`)
- v2_bridge_screenshot (imports `V2ForceKind`, `V2MantleSpec`)
- v2_workflow_r4_visual_checkpoint (full R4 + R5b + R6 + R7.A.2.4 sweep)
- v2_phase8g_visuals (imports `V2Field::VelocityMagnitude` etc — referenced earlier in session)
- v2_r7_omega_3_gradient_diagnostic (R7 ω.3 D — **init-only, does NOT import retired modules**; verify before gating)

### D5.b — Tests independent of retired modules

Checked: `v2_r7_omega_3_gradient_diagnostic` only imports `init`, `voronoi`, and the `V2InitModeSpec` enum (preserved). **Does not need gating.** This is the only R7 ω.3 era test that stands alone.

### D5.c — Recommended gating placement

Per-test `#![cfg(feature = "v2_legacy")]` at the top of each `.rs` file in `tests/`. Alternative: move all gated tests into a subdirectory `tests/_attic/` and add `[[test]]` entries with `required-features = ["v2_legacy"]` in `Cargo.toml`. The latter scales better with ~ 50 tests and signals intent more strongly.

**Decision needed** before Phase 2.

---

## D6 — JSON presets

All 8 V2 preset files in `crates/ymir-viz/presets/v2/` contain retired-mapped fields:

| File | `mantle` | `slab_enabled` | `force` | `linear_solver` | `cratonic.b_factor` |
|---|---|---|---|---|---|
| `active_medley.json` | On | false | `gpe` | `jacobi` | yes |
| `active_medley_composite.json` | On | false | `gpe` | `jacobi` | yes |
| `active_medley_orogenic.json` | On | false | `gpe` | `jacobi` | yes |
| `convergence.json` | On | false | `gpe` | `jacobi` | yes |
| `divergence.json` | On | false | `gpe` | `jacobi` | yes |
| `quiescent.json` | Off | false | `gpe` | `jacobi` | yes |
| `single_continent.json` | Off | false | `gpe` | `jacobi` | yes |
| `subduction.json` | On | false | `gpe` | `jacobi` | yes |

All presets are JSON data files, no Rust code. They are loaded by `presets.rs` (in `bridge/v2/`) and by `parameter_panel_v2.rs`. **They survive Phase 3 unchanged** because the deserialisation of `V2RunSpec` from JSON does not require the retired modules — only the subsequent `build()` call does.

**Recommendation:** keep `presets/v2/*.json` in place. Document in the migration doc that they only `build()` under `--features v2_legacy`. No move to `presets/v2_legacy/` needed.

The `b_factor` field on cratonic — explicitly noted in §4.8 as the "smoothstep amplification dropped" element. The field stays in the JSON (forward-compat with `--features v2_legacy`) but C1 ignores it.

---

## D7 — Documentation references

Files referencing retired-module paths (informational, no migration action; these become "historical" docs post-migration):

- [`docs/reports/step12_r7_final_report.md`](../reports/step12_r7_final_report.md) — Section 7 already filled with the exact §4.7/§4.8 lists.
- [`docs/design/c1_lightweight_dynamic_tectonics.md`](../design/c1_lightweight_dynamic_tectonics.md) — the source of truth for the retire/preserve split.
- [`docs/reports/step12_r6_mantle_evolution/R6_1_physics_decision.md`](../reports/step12_r6_mantle_evolution/) — R6 mantle evolution rationale.
- [`docs/reports/step12_solver_audit.md`](../reports/step12_solver_audit.md) — R5b D0 audit, source of the runtime motivation for C1.
- [`docs/reports/step8_5b_performance_report.md`](../reports/step8_5b_performance_report.md), [`docs/reports/step8_5a_amg_report.md`](../reports/step8_5a_amg_report.md), [`docs/PR_body_step8_5a.md`](../PR_body_step8_5a.md) — pre-C1 historical reports.

No documentation rewrite is needed for Phase 3; the docs accurately describe what's now in `_attic`.

---

## Per-file gating plan (Phase 3 Commit 2)

Files that need `#![cfg(feature = "v2_legacy")]` at the **file top** (because the entire file's content is Stokes-coupled):

**`crates/ymir-core/src/`:**

1. `tectonics_v2/diagnostics/harness.rs` — full file gate
2. `tectonics_v2/diagnostics/mms_bench.rs` — full file gate
3. `tectonics_v2/diagnostics/k_sub_sweep.rs` — full file gate
4. `tectonics_v2/diagnostics/num_plates_sweep.rs` — full file gate
5. `tectonics_v2/diagnostics/bi_sweep.rs` — full file gate
6. `tectonics_v2/diagnostics/br_sweep.rs` — full file gate
7. `tectonics_v2/diagnostics/ar_sweep.rs` — full file gate (uses preserved forcing types but pipes them into the retired harness)
8. `tectonics_v2/workflow/phase_a.rs` — full file gate (calls retired harness)
9. `tectonics_v2/workflow/phase_b.rs` — verify; if it consumes Phase A output as data only, may stay
10. `tectonics_v2/workflow/cycle.rs` and other sub-files of workflow — verify per-file

**Inside `tectonics_v2/diagnostics/mod.rs`** — gate the `pub mod harness;`, `pub mod mms_bench;`, `pub mod *_sweep;` declarations:

```rust
#[cfg(feature = "v2_legacy")]
pub mod harness;
#[cfg(feature = "v2_legacy")]
pub mod mms_bench;
// etc.
```

**Inside `tectonics_v2/forcing/mod.rs`** — gate four lines:

```rust
#[cfg(feature = "v2_legacy")]
pub mod mantle_force;
#[cfg(feature = "v2_legacy")]
pub mod slab_pull;
#[cfg(feature = "v2_legacy")]
pub use mantle_force::MantleForce;
#[cfg(feature = "v2_legacy")]
pub use slab_pull::SlabPullForce;
```

**`crates/ymir-viz/src/bridge/v2/`:**

11. `build_config.rs` — full file gate
12. `thread.rs` — gate the v2 worker thread loop body or full file
13. `events.rs` — gate the `V2Event::Progress`, `V2Event::Completed`, `V2Event::WorkflowPhaseACompleted`, `V2Event::WorkflowPhaseBCompleted` variants if their fields touch retired types; **partial gate** (the enum stays for downstream consumers)
14. `commands.rs` — gate the `V2Command::RunBaseline`, `V2Command::RunWorkflowPhaseA`, `V2Command::RunWorkflowPhaseB` variants
15. `plugin.rs` — gate the Bevy plugin v2 panels' mount points
16. `mod.rs` — gate the `pub mod` declarations + selective `pub use` of gated items

**`crates/ymir-viz/src/bridge/v2/spec.rs`:** stays compilable as plain serde data (no retired imports needed at type-declaration level). Gate only the `into_core()` impls that target retired types.

**`crates/ymir-viz/src/ui/parameter_panel_v2.rs`:** full file gate (the panel only makes sense if you can run a v2 simulation).

**Binaries** (in `crates/ymir-core/src/bin/`):

17. `step5_baseline.rs`, `step6_baseline.rs`, `step7_baseline.rs`, `step8_baseline.rs`, `step_baseline.rs` — gate the binary target in `Cargo.toml` with `required-features = ["v2_legacy"]`
18. `gen_bench_data.rs`, `gen_reference_solutions.rs` — same

**Benches** (in `crates/ymir-core/benches/`):

19. `amg_benchmark.rs` — gate in `Cargo.toml`

**Tests** (in both crates):

20. ~ 50 test files — gate via `#![cfg(feature = "v2_legacy")]` at file top, OR relocate to `tests/_attic/` with `required-features` in `Cargo.toml`. **Recommend the latter** for clearer signal and easier discoverability.

---

## Hidden couplings (W3 surfacing)

### HC1 — `tectonics_v2/slab/` is missing from §4.7

**Evidence:** `forcing/slab_pull.rs` docstring (line 8) cites `super::super::slab::state::SlabState` as the cell-centred field this body force depends on. `slab/mod.rs` is the slab-mass ODE + advection + convergence-direction machinery. The two modules cannot operate independently.

**Files that would orphan without retiring `slab/`:** the binaries `step7_baseline.rs`, `step8_baseline.rs`, `gen_bench_data.rs` all import `tectonics_v2::slab::SlabPullConfig` (the config enum) plus the rest. Tests `v2_slab_*` are slab-mass-mechanism-specific.

**Recommendation:** **add `tectonics_v2/slab/` to §4.7 retired list**. Treat slab subtree like mantle subtree — both are Stokes-coupled tectonic mechanisms with their own state + config + body-force wrapper.

### HC2 — `tectonics_v2/presets.rs` is in neither §4.7 nor §4.8

Contains `RheologyParams`, `YieldingConfig`, `Preset`. The presets defined here (`dynamic-accidented`, `stable-shield`, `soft-planet`) carry rheology-typed payloads that retire with `rheology.rs`. The file is imported by `bridge/v2/build_config.rs:13` for `Preset` and `YieldingConfig`.

**Recommendation:** **retire `tectonics_v2/presets.rs`** alongside `rheology.rs`. The C1 presets will live in a separate file.

### HC3 — `boundary_detection.rs`, `cancel.rs`, `scales.rs` are in neither list

Audited: these are small utility modules.

- `cancel.rs` — cancellation token framework, paradigm-agnostic, **preserved**.
- `scales.rs` — non-dim scale references, paradigm-agnostic, **preserved**.
- `boundary_detection.rs` — used by Voronoï classification; checked imports, paradigm-agnostic, **preserved**.

No issue; these were just absent from both lists. Implicit "preserved".

### HC4 — `workflow/` sub-files (RESOLVED by Phase 1.3 H2, Issue #125)

**Resolution status: RESOLVED.** Phase 1.3 H2 (Issue #125) shipped
the paradigm-agnostic workflow refactor on 2026-05-22 via 5 commits
on branch `125-c1-phase-13-equilibrium-height-closure-harness-paradigm-agnostic-refactor`.
Two-separate-functions approach (Option B), not trait abstraction.
See [`harness_paradigm_agnostic.md`](./harness_paradigm_agnostic.md)
for the full audit and § 12 of that doc for the implementation
record (SHAs `11c53b6` → `6ef32a0`).

Key outcomes:

- `workflow/phase_a.rs` was renamed to `phase_a_v2.rs` and stays
  gated under `v2_legacy`. A sibling `phase_a_c1.rs`
  (default-features-on) drives the C1 paradigm via
  `tectonics_c1::time_loop::run_with_closures`.
- `workflow/phase_b.rs` had its `v2_legacy` gate dropped — Phase B
  was already paradigm-agnostic at runtime (consumes `Field2D +
  sea_level + iso_config`, no `BaselineResult`).
- A shared `workflow/phase_a_common.rs` carries the paradigm-
  agnostic post-tectonic pass (sea-level → macro-redistribution
  → reclassification → cratonic recompute) consumed by both
  paradigm-specific Phase A entries.
- 58 / 58 regression tests PASS across both feature flags,
  including a bit-identical decomposition contract on the C1
  path that matches v2's bit-identical Disabled contract.

### HC4 (original audit — historical record)

§4.8 lists `tectonics_v2/workflow/*` as preserved. But the workflow's `phase_a.rs` calls `harness::run_baseline_with_progress`, which is retired. **The workflow orchestrator is conserve-architecturally but cannot execute without `v2_legacy`** until C1 provides its own per-cycle tectonic runner.

**Recommendation (since acted on):** add to the migration doc a note that the workflow shell is paradigm-agnostic *by design* but currently couples to v2 at run-time. C1 Phase 1.3+ (§7 C1.md) is where the shell becomes truly paradigm-portable. For the migration itself, gate `workflow/phase_a.rs` (and likely `phase_b.rs` if it consumes via harness too) and let C1 reintroduce a paradigm-agnostic form when the tectonic_c1 runner exists.

---

## Special cases

### SC1 — `harness.rs` cannot be partially gated

`BaselineConfig` carries fields of types `MantleConfig`, `SlabPullConfig`, `BasalDragConfig`, `YieldingConfig`. Even the data struct couples to retired surface. Splitting into "paradigm-agnostic shell + v2 implementation" is a substantial refactor that belongs in C1 Phase 1.3, not in the migration.

**Decision:** entire-file gate for the migration commit. Accept that the bridge's V2 run path is gated as a consequence.

### SC2 — Test relocation strategy: in-place gate vs `tests/_attic/`

Two viable approaches, each with trade-offs:

- **In-place gate** (`#![cfg(feature = "v2_legacy")]` at top of each test file):
  - Pro: minimal `git mv` churn; tests stay in their canonical location
  - Con: 50+ file headers, every future contributor sees v2 tests "side by side" with C1 tests
- **Move to `tests/_attic/`**:
  - Pro: clear visual signal "this is legacy code"; `Cargo.toml` `required-features` enforces gating
  - Con: 50+ file moves; tooling pointing at `tests/` may need to recursively look

**Recommendation: move to `tests/_attic/`** — aligned with the file-level `_attic/` pattern used for production code. Surface the decision to user for confirmation.

### SC3 — JSON presets default vs gated

Already addressed in D6: JSON presets stay in `presets/v2/`. Only consumer code gates.

### SC4 — `V2RunSpec` round-trip without `v2_legacy`

The V2 spec enum declarations and serde derives compile fine without `v2_legacy` (they don't import retired modules at type level — they translate via `into_core()` methods that ARE gated). So:

```bash
# Works without feature: serde round-trip + JSON validation
cargo test --workspace -p ymir-viz spec_roundtrip

# Requires feature: actually building a runnable BaselineConfig
cargo test --features v2_legacy --workspace
```

This preserves an "audit-without-running" capability that may be useful for the C1 development.

---

## Recommendations summary

1. **§4.7 amendment requested:**
   - Add `tectonics_v2/slab/` (HC1)
   - Add `tectonics_v2/presets.rs` (HC2)
2. **§4.8 amendment requested:**
   - Note that `workflow/*` is paradigm-agnostic *by design* but requires C1's tectonic runner to operate without `v2_legacy` (HC4)
   - Note that `diagnostics/` mod still hosts harness/mms_bench/sweeps, all of which gate (the umbrella stays preserved but its v2-coupled members are migrated to gated content)
3. **Test relocation:** move ~ 50 v2 tests from `crates/{ymir-core,ymir-viz}/tests/*.rs` to `crates/{ymir-core,ymir-viz}/tests/_attic/*.rs` with `required-features = ["v2_legacy"]` in each `[[test]]` Cargo entry (SC2)
4. **Bins + benches:** add `required-features = ["v2_legacy"]` to the v2-baseline / gen / bench targets in Cargo.toml
5. **JSON presets:** keep in place; document the implicit `--features v2_legacy` requirement in `presets/v2/README.md` (creating if absent)
6. **Harness gating:** entire file (SC1), accepting the bridge consequence
7. **Bridge layer gating:** full gate of `build_config.rs`, `thread.rs`, `parameter_panel_v2.rs`, plus per-variant gates in `events.rs`, `commands.rs`, `plugin.rs`, `mod.rs`. Spec layer stays compilable.

## Phase 2/3 readiness

- Phase 2 (`git mv` to `_attic/`) is **blocked** pending user decisions on:
  - HC1 — retire `slab/`? (Recommend yes)
  - HC2 — retire `presets.rs`? (Recommend yes)
  - SC2 — relocate tests to `tests/_attic/` or in-place gate? (Recommend relocate)
- Phase 3 (gating) inherits Phase 2 decisions. ~ 50 files touched; budget ~ 2-3 hours focused work for the gating + Cargo.toml edits.
- Phase 4 (regression verification) — run `cargo test --features v2_legacy --workspace`, capture summary.md outputs, diff against committed baselines. Cargo.lock + rustc toolchain version invariant must be documented.

**Estimated total migration effort:** 4-6 hours focused work across 3 commits. No surprises beyond HC1 and HC2.

---

*End of Phase 1 audit. Awaiting user decisions on §4.7/§4.8 amendments and SC2 test-relocation strategy before authorising Phase 2.*
