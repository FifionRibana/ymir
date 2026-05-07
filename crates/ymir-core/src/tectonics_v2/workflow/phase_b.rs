//! Phase B — HD finalization (one-shot upscale + erosion).
//!
//! Takes the final-cycle output of Phase A (a low-res
//! [`crate::tectonics_v2::diagnostics::harness::BaselineResult`]),
//! converts the `Field2D` S̃ field to `GridF32`, runs
//! [`crate::terrain::upscale::upscale_with_fbm`] to the configured HD
//! resolution, then [`crate::erosion::hydraulic::run_erosion`] on the
//! upscaled heightmap. Records the D5 grand-scale deviation
//! `‖HD_after - upscale(low_res)‖_∞`.
//!
//! Phase 1 (this commit) ships only the `Disabled` passthrough
//! (`None`); the actual HD pipeline lands in Step 12 Phase 5.

use super::{PhaseBOutput, WorkflowConfig};
use crate::tectonics_v2::diagnostics::harness::BaselineResult;

/// Run Phase B HD finalization on a Phase A output.
///
/// `Disabled` → returns `None`. The user is expected to consume the
/// low-res `BaselineResult` directly (Step 11 contract).
///
/// `Enabled(_)` → returns `Some(PhaseBOutput)`. Phase 5 work,
/// currently `unimplemented!`.
pub fn run_phase_b(_input: &BaselineResult, wf: &WorkflowConfig) -> Option<PhaseBOutput> {
    match wf {
        WorkflowConfig::Disabled => None,
        WorkflowConfig::Enabled(_) => {
            unimplemented!("Phase B HD finalization lands in Step 12 Phase 5");
        }
    }
}
