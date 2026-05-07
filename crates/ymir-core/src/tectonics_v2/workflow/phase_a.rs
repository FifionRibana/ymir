//! Phase A — low-res loop orchestration.
//!
//! Single-cycle and multi-cycle entry points. Phase 1 ships only the
//! `Disabled` passthrough: [`run_phase_a_cycle`] with
//! [`super::WorkflowConfig::Disabled`] is exactly one
//! [`crate::tectonics_v2::diagnostics::harness::run_baseline`] call,
//! [`run_phase_a_loop`] is exactly one cycle. The `Enabled(_)` branch
//! is wired in Phase 3 (single-cycle orchestration) and Phase 4
//! (multi-cycle loop with [`crate::tectonics_v2::diagnostics::harness::ContinuationState`]
//! warm-start between cycles).

use super::{CycleOutput, PhaseAOutput, WorkflowConfig};
use crate::tectonics_v2::diagnostics::harness::{run_baseline, BaselineConfig};

/// Run a single Phase A cycle.
///
/// `Disabled` → direct `run_baseline(cfg)` passthrough wrapped in a
/// [`CycleOutput`] with `erosion_volume_removed = 0.0`. The bit-
/// identical regression contract: every byte of the returned
/// `baseline.final_state` matches a parallel `run_baseline(cfg)` call.
///
/// `Enabled(_)` → Phase 3 work. Currently `unimplemented!`.
pub fn run_phase_a_cycle(cfg: &BaselineConfig, wf: &WorkflowConfig) -> CycleOutput {
    match wf {
        WorkflowConfig::Disabled => {
            let baseline = run_baseline(cfg);
            CycleOutput { baseline, erosion_volume_removed: 0.0 }
        }
        WorkflowConfig::Enabled(_) => {
            unimplemented!(
                "Phase A single-cycle orchestration (tectonic + isostasy + \
                 low_res_erosion + reclassify + recompute craton) lands in \
                 Step 12 Phase 3"
            );
        }
    }
}

/// Run the Phase A multi-cycle loop.
///
/// `Disabled` → exactly one cycle (single `run_baseline` passthrough).
/// `output.cycles.len() == 1` and `output.cycles[0]` is the direct
/// passthrough `CycleOutput`.
///
/// `Enabled(params)` → loop `params.phase_a.n_cycles` cycles, each
/// running `params.phase_a.k_cycle` tectonic steps before the cycle's
/// erosion pass. Continuation between cycles uses
/// [`crate::tectonics_v2::diagnostics::harness::ContinuationState`]
/// (Step 8.6 infrastructure). Phase 4 work, currently `unimplemented!`.
pub fn run_phase_a_loop(cfg: &BaselineConfig, wf: &WorkflowConfig) -> PhaseAOutput {
    match wf {
        WorkflowConfig::Disabled => {
            let cycle = run_phase_a_cycle(cfg, wf);
            PhaseAOutput { cycles: vec![cycle] }
        }
        WorkflowConfig::Enabled(_) => {
            unimplemented!("Phase A multi-cycle loop lands in Step 12 Phase 4");
        }
    }
}
