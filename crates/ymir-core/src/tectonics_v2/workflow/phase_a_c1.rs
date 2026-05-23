//! Phase A — low-res loop orchestration, **C1 path**.
//!
//! Per Phase 1.3 H2 (Issue #125), this module is the C1-paradigm
//! Phase A entry point — **default-features-on, no `v2_legacy`
//! gate**. The structurally-identical v2-paradigm path lives at
//! [`super::phase_a_v2`] (gated under `v2_legacy`).
//!
//! ## Asymmetric API by design
//!
//! Per the H1 audit's R1 mitigation, the C1 path API is **not**
//! force-symmetric with v2. The two paradigms have structurally
//! different state and progress shapes:
//!
//! | Aspect              | v2 path                          | C1 path (this module)              |
//! |---------------------|----------------------------------|------------------------------------|
//! | Tectonic input      | `&BaselineConfig` (23-field)     | `&mut C1State` + `&PlateKinematics` + `&C1Closures` + `&C1TimeLoopConfig` |
//! | Per-step progress   | `FnMut(&StepProgress<'_>) -> bool` (Stokes residuals, nonlinear iter count, …) | not exposed — see "deferred" note below |
//! | Result envelope     | `CycleOutputV2 { baseline: BaselineResult, common: CycleOutputCommon }` | `PhaseACycleOutputC1 { common, new_cratonic_factor }` — no `baseline` |
//! | State threading     | `final_state_to_continuation_v2` (returns `ContinuationState` for cycle N+1 warm-start) | implicit: caller owns `&mut C1State` across cycles |
//!
//! The shared post-tectonic-step pass
//! ([`super::phase_a_common::apply_post_tectonic`]) handles the
//! paradigm-agnostic part (sea-level → macro-redistribution →
//! reclassification → cratonic recompute) and is called identically
//! by both paths.
//!
//! ## On_progress variant deferred to Phase 4
//!
//! The v2 path ships both `run_phase_a_cycle_v2` and
//! `run_phase_a_cycle_with_progress_v2` so the viz bridge thread
//! can stream `V2Event::Progress` events to the metrics dashboard.
//! The C1 path ships only the no-progress variant for Phase 1.3 —
//! the C1 viz bridge will be added in a later phase (Phase 4 UI
//! integration per the milestone roadmap), and at that point a
//! `run_phase_a_cycle_with_progress_c1` callback signature can be
//! designed in concert with the bridge consumers. The C1
//! `time_loop::run_with_closures` already takes a per-step
//! callback (`FnMut(usize, &C1State)`), so the lift is trivial
//! once the consumer side is ready.

use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_v2::cratonic::CratonicConfigEnabled;
use crate::tectonics_v2::field::Field2D;

use crate::tectonics_c1::kinematics::PlateKinematics;
use crate::tectonics_c1::state::C1State;
use crate::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};

use super::phase_a_common::{apply_post_tectonic, extract_per_plate_type, PostTectonicInput};
use super::{CycleOutputCommon, WorkflowConfig};

/// Inputs to [`run_phase_a_cycle_c1`].
///
/// `state` is mutated in place by the tectonic step (advection
/// of `S̃` and `age`, plus per-cell closures from `closures`).
/// When `wf == Enabled(_)`, the shared post-tectonic pass mutates
/// `state.s` again (macro-redistribution) and `state.plate_type`
/// (reclassification).
///
/// `cratonic_config` is optional: C1 Phase 1.3 has no cratonic-
/// factor field on `C1State` (only a `BoolField cratonic_mask`),
/// so the recomputed factor returned in
/// [`PhaseACycleOutputC1::new_cratonic_factor`] would be discarded
/// by the caller. Set to `None` to skip the recompute step
/// entirely; Phase 1.4+ may add a `Field2D` cratonic factor to
/// `C1State` and consume the output.
pub struct PhaseACycleInputC1<'a> {
    pub state: &'a mut C1State,
    pub kinematics: &'a PlateKinematics,
    pub closures: &'a C1Closures,
    pub time_loop_config: &'a C1TimeLoopConfig,
    pub iso_config: &'a IsostasyConfig,
    pub cratonic_config: Option<&'a CratonicConfigEnabled>,
}

/// Output of [`run_phase_a_cycle_c1`].
///
/// Asymmetric vs v2's `CycleOutputV2`: no `baseline` field because
/// C1 mutates `&mut C1State` in place — the caller already holds
/// the post-cycle state.
pub struct PhaseACycleOutputC1 {
    pub common: CycleOutputCommon,
    /// Recomputed cratonic factor, or `None` when
    /// `input.cratonic_config` was `None` (or `wf == Disabled`).
    /// Discarded by C1 Phase 1.3 callers (no consuming field on
    /// `C1State`); reserved for Phase 1.4+.
    pub new_cratonic_factor: Option<Field2D>,
}

/// Run a single Phase A cycle on the C1 path.
///
/// `Disabled` → run the C1 tectonic sub-cycle (advection + the
/// per-cell closures in `closures`) and return an empty
/// `CycleOutputCommon` — no macro-redistribution, no
/// reclassification, no cratonic recompute. This is the parallel
/// of v2's `WorkflowConfig::Disabled` bit-identical-to-run-baseline
/// contract: `run_phase_a_cycle_c1(input, Disabled)` ≡ direct
/// `run_with_closures(state, ..., |_, _| {})`. See H3 Commit 5 for
/// the regression test.
///
/// `Enabled(params)` → 2-step pipeline:
///   1. Tectonic — `run_with_closures` (C1 advection + closures).
///   2. Shared post-tectonic pass — delegates to
///      [`apply_post_tectonic`] for sea-level, macro-redistribution,
///      reclassification, and cratonic recompute. The freshly-
///      recomputed cratonic factor (if any) is returned in
///      [`PhaseACycleOutputC1::new_cratonic_factor`] for the
///      caller to install or discard.
pub fn run_phase_a_cycle_c1(
    input: PhaseACycleInputC1<'_>,
    wf: &WorkflowConfig,
) -> PhaseACycleOutputC1 {
    let PhaseACycleInputC1 {
        state,
        kinematics,
        closures,
        time_loop_config,
        iso_config,
        cratonic_config,
    } = input;

    // Step 1 — Tectonic (C1 forward-Euler advection + closures).
    // C1 has no continuation/warm-start concept (the closures are
    // per-cell additive and the kinematics are static); the loop
    // mutates `state` in place.
    run_with_closures(state, kinematics, time_loop_config, closures, |_, _| {});

    match wf {
        WorkflowConfig::Disabled => PhaseACycleOutputC1 {
            common: CycleOutputCommon::default(),
            new_cratonic_factor: None,
        },
        WorkflowConfig::Enabled(params) => {
            // Capture the pre-reclassification per-plate type for
            // the D4 "was continental at init" gate. Must happen
            // BEFORE `apply_post_tectonic` mutates `plate_type`.
            let original_per_plate_type =
                extract_per_plate_type(&state.plate_id, &state.plate_type);

            // Step 2 — Shared post-tectonic pass. The struct
            // literal holds three disjoint borrows of `state` (mut
            // `s`, shared `plate_id`, mut `plate_type`), which the
            // borrow checker accepts under the splitting-borrow
            // rule.
            let post = apply_post_tectonic(PostTectonicInput {
                s_field: &mut state.s,
                plate_id: Some(&state.plate_id),
                plate_type: Some(&mut state.plate_type),
                previous_cratonic_factor: None,
                initial_per_plate_type: Some(&original_per_plate_type),
                params: &params.phase_a,
                iso_cfg: iso_config,
                cratonic_cfg: cratonic_config,
            });

            PhaseACycleOutputC1 {
                common: post.common,
                new_cratonic_factor: post.new_cratonic_factor,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::tectonics_c1::init::init_c1_state_phase_1_1;
    use crate::tectonics_v2::workflow::WorkflowParams;

    /// Smoke test: a 300-step C1 cycle through the full
    /// Phase A pipeline (`Enabled`) completes without producing
    /// NaN/Inf, lands within macro_redistribution's mass-
    /// conservation budget, and populates the
    /// `CycleOutputCommon` scalars.
    ///
    /// This is a sanity check, **not** a Phase 1.3 acceptance
    /// test — those live in Stage E3.
    #[test]
    fn c1_phase_a_cycle_completes_300_steps() {
        let grid = 64;
        let mut state = init_c1_state_phase_1_1(grid, 42);
        let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let time_loop_config = C1TimeLoopConfig {
            n_steps: 300,
            dx: 1.0 / grid as f64,
            dy: 1.0 / grid as f64,
        };
        let iso_config = IsostasyConfig::default();
        let wf = WorkflowConfig::Enabled(WorkflowParams::default());

        let output = run_phase_a_cycle_c1(
            PhaseACycleInputC1 {
                state: &mut state,
                kinematics: &kinematics,
                closures: &closures,
                time_loop_config: &time_loop_config,
                iso_config: &iso_config,
                cratonic_config: None,
            },
            &wf,
        );

        // Sanity 1 — no NaN / Inf in `S̃`.
        for &v in state.s.data() {
            assert!(
                v.is_finite(),
                "non-finite S̃ in state after 300 steps + post-pass"
            );
        }

        // Sanity 2 — sea-level was computed by the Enabled
        // post-pass (Disabled would have left this at 0.0).
        assert!(
            output.common.sea_level_normalized > 0.0,
            "sea_level_normalized should be set by Enabled post-pass; got {}",
            output.common.sea_level_normalized
        );

        // Sanity 3 — macro_redistribution is mass-conserving by
        // construction (Step 12 R3 drainage + isostatic rebound).
        // The `mass_drift` field measures the pre-vs-post-macro
        // delta inside `apply_post_tectonic` and should be at
        // machine-precision floor relative to the integrated
        // mass.
        let total_mass: f64 = state.s.data().iter().sum();
        let drift_budget = total_mass.abs() * 1e-9;
        assert!(
            output.common.mass_drift.abs() < drift_budget,
            "macro_redistribution mass drift {} exceeds {:.3e} budget (1e-9 × |mass| = {:.3e}); macro_redistribution must conserve mass",
            output.common.mass_drift, drift_budget, drift_budget
        );

        // Sanity 4 — cratonic_config was None → no recompute.
        assert!(
            output.new_cratonic_factor.is_none(),
            "C1 Phase 1.3 default-config call must produce no cratonic factor"
        );
        assert!(
            output.common.craton_recomputation_change.is_none(),
            "C1 Phase 1.3 default-config call must produce no craton change"
        );
    }
}
