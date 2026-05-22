//! Phase A — low-res loop orchestration, **v2 path**.
//!
//! Per Phase 1.3 H2 (Issue #125), this module is the v2-paradigm
//! Phase A entry point — gated under `v2_legacy`. The structurally-
//! identical C1 paradigm path lives at [`super::phase_a_c1`]
//! (default-features-on, Commit 4).
//!
//! `run_phase_a_cycle_v2` chains the 5-step single-cycle pipeline:
//!
//! 1. **Tectonic** — `run_baseline(cfg)`; the cfg may carry a
//!    [`ContinuationState`] for cycle-to-cycle warm-start (D3).
//!    **v2-specific** (Stokes + nonlinear continuation).
//! 2. **Post-tectonic pass (shared)** — delegates to
//!    [`super::phase_a_common::apply_post_tectonic`]: sea-level →
//!    macro-redistribution → reclassification → cratonic
//!    recompute. Paradigm-agnostic; the C1 path calls the same
//!    function in Commit 4.
//! 3. **Velocity reset** — sets `vx = vy = 0` after macro-
//!    redistribution. **v2-specific** (D1-ter empirical finding,
//!    see in-function comment for the rationale).
//! 4. **Cratonic factor install** — assigns the freshly-recomputed
//!    factor to `baseline.final_state.cratonic_factor`. **v2-
//!    specific** (data shape: v2 stores `Field2D`, C1 has only
//!    `BoolField cratonic_mask`).
//!
//! `WorkflowConfig::Disabled` short-circuits the entire pipeline:
//! the cycle is exactly `run_baseline(cfg)` with all extra scalars
//! at zero/`None`. The regression
//! `v2_workflow_disabled_regression::workflow_disabled_run_phase_a_cycle_is_bit_identical_to_run_baseline`
//! pins this contract byte-for-byte.

use super::phase_a_common::{apply_post_tectonic, extract_per_plate_type, PostTectonicInput};
use super::{CycleOutputCommon, CycleOutputV2, PhaseAOutputV2, WorkflowConfig};
use crate::tectonics::isostasy::IsostasyConfig;
use crate::tectonics_v2::cratonic::CratonicConfig;
use crate::tectonics_v2::diagnostics::harness::{
    run_baseline_with_progress, BaselineConfig, ContinuationState, FinalState, StepProgress,
};

/// Run a single Phase A cycle on the v2 path. Thin wrapper over
/// [`run_phase_a_cycle_with_progress_v2`] with a no-op callback that
/// never aborts, preserving the bit-identical regression contract
/// (acceptance #15) byte-for-byte: `run_baseline_with_progress(cfg,
/// |_| true)` is itself a wrapper over `run_baseline` from Step 8.6
/// follow-up, so the call chain reduces to the same primitive.
///
/// `Disabled` → direct `run_baseline(cfg)` passthrough wrapped in a
/// [`CycleOutputV2`] with all extra scalars at zero/`None`.
///
/// `Enabled(params)` → 4-step pipeline (tectonic → shared post-
/// tectonic pass → velocity reset → cratonic factor install).
/// Returns the post-cycle state suitable for cycle-to-cycle
/// continuation via [`final_state_to_continuation_v2`].
pub fn run_phase_a_cycle_v2(cfg: &BaselineConfig, wf: &WorkflowConfig) -> CycleOutputV2 {
    run_phase_a_cycle_with_progress_v2(cfg, wf, |_| true)
}

/// Streaming variant of [`run_phase_a_cycle_v2`]. The callback fires
/// once per completed harness step inside the cycle's tectonic
/// sub-phase (step 1 of the 4-step pipeline); returning `false`
/// requests a graceful abort of the harness step loop. Same
/// callback shape as
/// [`crate::tectonics_v2::diagnostics::harness::run_baseline_with_progress`].
///
/// Added in Step 12 follow-up so the v2 bridge can stream per-step
/// `V2Event::Progress` to the metrics dashboard during Phase A
/// (the dashboard previously froze between `WorkflowCycleCompleted`
/// events because `run_phase_a_cycle_v2` invoked `run_baseline` —
/// the `|_| true` callback wrapper — with no streaming hook). The
/// post-tectonic substeps (isostasy, erosion, reclassify, craton
/// recompute) are not currently streamed; they're sub-second on
/// 64² mantle-on, so a single "cycle progress" tick is the
/// pragmatic granularity.
pub fn run_phase_a_cycle_with_progress_v2<F>(
    cfg: &BaselineConfig,
    wf: &WorkflowConfig,
    on_progress: F,
) -> CycleOutputV2
where
    F: FnMut(&StepProgress<'_>) -> bool,
{
    match wf {
        WorkflowConfig::Disabled => {
            let baseline = run_baseline_with_progress(cfg, on_progress);
            CycleOutputV2 {
                baseline,
                common: CycleOutputCommon::default(),
            }
        }
        WorkflowConfig::Enabled(params) => {
            // Step 1 — Tectonic (v2-specific).
            let mut baseline = run_baseline_with_progress(cfg, on_progress);

            // Capture the pre-reclassification per-plate type for
            // the D4 "was continental at init" gate. Must happen
            // BEFORE the shared post-tectonic pass mutates
            // `plate_type` via reclassification.
            let original_per_plate_type = baseline
                .final_state
                .plate_type
                .as_ref()
                .zip(baseline.final_state.plate_id.as_ref())
                .map(|(pt, pid)| extract_per_plate_type(pid, pt));

            let cratonic_cfg_enabled = match &cfg.cratonic {
                CratonicConfig::Enabled(c) => Some(c),
                CratonicConfig::Disabled => None,
            };

            // Step 2 — Shared post-tectonic pass (paradigm-
            // agnostic). Sea-level → macro-redistribution →
            // reclassification → cratonic recompute. The
            // `final_state` fields are destructured so the borrow
            // checker accepts the disjoint mut borrows of
            // `s_field` and `plate_type` simultaneously.
            let FinalState {
                ref mut s_field,
                ref plate_id,
                ref mut plate_type,
                ref cratonic_factor,
                ..
            } = baseline.final_state;
            let iso_cfg = IsostasyConfig::default();
            let pt_result = apply_post_tectonic(PostTectonicInput {
                s_field,
                plate_id: plate_id.as_ref(),
                plate_type: plate_type.as_mut(),
                previous_cratonic_factor: cratonic_factor.as_ref(),
                initial_per_plate_type: original_per_plate_type.as_deref(),
                params: &params.phase_a,
                iso_cfg: &iso_cfg,
                cratonic_cfg: cratonic_cfg_enabled,
            });

            // Step 3 — Velocity reset (v2-specific, D1-ter).
            //
            // EMPIRICAL FINDING (counter-intuitive vs classical
            // Stokes wisdom): the warm-start `v = v_final_previous_
            // cycle` is not just sub-optimal post-macro, it is
            // **actively harmful**. `macro_redistribution::apply`
            // shifts S̃ enough that the GPE driver direction
            // changes; `v_warm_start` points in a direction now
            // anti-useful for the next tectonic step, and Newton
            // oscillates trying to correct it (sub-case C amplified
            // in D2-bis classification).
            //
            // The 3-variant D1-ter benchmark (commit 4969de9)
            // showed `v = 0` gives:
            //   - cycle 2 Converged 45/45 vs 14/41 with warm-start
            //   - 0 Oscillating vs 26 with warm-start
            //   - CG iter total over first 5 cycle-2 steps:
            //     28k vs 75k
            //   - ‖Δv‖/‖v‖ max: 1.00 vs 11.43
            // Variant C (Gaussian smoothing of v) is WORSE than
            // warm-start (90k CG iter, peak |v| explosion to 53) —
            // smoothing preserves the wrong direction.
            //
            // Gated by `WorkflowConfig::Enabled` (this match arm).
            // The Disabled branch is untouched and the bit-
            // identical regression
            // `v2_workflow_disabled_regression` continues to hold.
            //
            // Counter-intuitive — leave the comment block intact;
            // a future dev tempted to "re-enable warm-start because
            // it's faster on Stokes" would re-break the system.
            // See `docs/reports/step12_solver_audit.md` § F and
            // `docs/reports/step12_r5b_d1_ter_init_variants/` for
            // the full empirical record.
            for v in baseline.final_state.vx.iter_mut() {
                *v = 0.0;
            }
            for v in baseline.final_state.vy.iter_mut() {
                *v = 0.0;
            }

            // Step 4 — Install the new cratonic factor (v2-
            // specific data shape).
            if let Some(new_factor) = pt_result.new_cratonic_factor {
                baseline.final_state.cratonic_factor = Some(new_factor);
            }

            CycleOutputV2 {
                baseline,
                common: pt_result.common,
            }
        }
    }
}

/// Run the Phase A multi-cycle loop on the v2 path.
///
/// `Disabled` → exactly one cycle (single `run_baseline` passthrough).
/// The `&mut` requirement is preserved on this branch even though no
/// mutation actually fires, because the [`WorkflowConfig::Enabled`]
/// branch must mutate `cfg.continuation` between cycles to wire the
/// D3 warm-start contract.
///
/// `Enabled(params)` → loop `params.phase_a.n_cycles` cycles. After
/// each cycle (except the last) the loop sets
/// `cfg.continuation = Some(final_state_to_continuation_v2(...))` so
/// the next cycle's `run_baseline` warm-starts from the prior
/// cycle's post-erosion state. The S̃ field, velocity, age and
/// cratonic factor all thread through (D3 contract pinned by
/// `v2_workflow_continuation_no_transient`).
///
/// `cfg.steps` is consumed as the number of tectonic steps per cycle.
/// The convention is to set `cfg.steps = params.phase_a.k_cycle`
/// before calling, but the loop does not enforce this — the two are
/// independently configurable.
pub fn run_phase_a_loop_v2(cfg: &mut BaselineConfig, wf: &WorkflowConfig) -> PhaseAOutputV2 {
    match wf {
        WorkflowConfig::Disabled => {
            let cycle = run_phase_a_cycle_v2(cfg, wf);
            PhaseAOutputV2 { cycles: vec![cycle] }
        }
        WorkflowConfig::Enabled(params) => {
            let n_cycles = params.phase_a.n_cycles.max(1);
            let mut cycles: Vec<CycleOutputV2> = Vec::with_capacity(n_cycles);
            for cycle_idx in 0..n_cycles {
                let cycle = run_phase_a_cycle_v2(cfg, wf);
                // Set up the next cycle's warm-start *before* moving
                // `cycle` into the output vec. Skip the last cycle:
                // there is no next cycle to warm-start.
                if cycle_idx + 1 < n_cycles {
                    cfg.continuation =
                        Some(final_state_to_continuation_v2(&cycle.baseline.final_state));
                }
                cycles.push(cycle);
            }
            PhaseAOutputV2 { cycles }
        }
    }
}

/// Build a [`ContinuationState`] from a [`FinalState`].
///
/// The orchestrator (Phase 4) calls this at the end of cycle `N` to
/// build the input for cycle `N+1`'s `BaselineConfig.continuation`.
/// Step 8.6's `ContinuationState` carries everything `run_baseline`
/// needs to short-circuit re-init: `s, vx, vy, age, cratonic_factor`.
/// The Voronoï tessellation is implicitly preserved by the run's
/// `BoundaryConfig` (static for the run lifetime).
///
/// D3 contract: cycle `N+1` step 1 should produce a peak|v| within
/// 10 % of cycle `N` step `k_cycle` — pinned by the
/// `v2_workflow_continuation_no_transient` test.
///
/// **v2-specific.** C1 has no Stokes velocity field; its
/// equivalent threading function (if needed in Phase 1.4+) will
/// live in `phase_a_c1.rs`.
pub fn final_state_to_continuation_v2(fs: &FinalState) -> ContinuationState {
    ContinuationState {
        s: fs.s_field.clone(),
        vx: fs.vx.clone(),
        vy: fs.vy.clone(),
        age: fs.age_field.clone(),
        cratonic_factor: fs.cratonic_factor.clone(),
    }
}
