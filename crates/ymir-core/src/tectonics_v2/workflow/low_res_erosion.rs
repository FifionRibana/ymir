//! Phase A low-res parametric erosion (D2 of `step12_issue.md`).
//!
//! Per-cycle diffusive erosion algorithm applied to all continental
//! cells:
//!
//! ```text
//! For each cell i with S̃[i] > sea_level_reference:
//!     slope_i = max gradient of S̃ in 4-neighborhood
//!     Δh = α · slope_i · (S̃[i] - sea_level_reference)
//!     S̃[i] -= Δh
//!     if β > 0:
//!         downslope_neighbor.S̃ += β · Δh
//! ```
//!
//! The `sea_level_reference` is **not** hardcoded to `0.5` — it comes
//! from `compute_isostasy(s).sea_level_normalized` (Option 2 of the
//! Phase 0 architectural finding E). This makes the threshold
//! adaptive to mass drift cycle-after-cycle, which matters because
//! low-res erosion is non-conservative when `β = 0` and the
//! continental fraction can drift over a multi-cycle run.
//!
//! Phase 1 (this commit) ships only the module skeleton; the
//! algorithm itself lands in Step 12 Phase 2 along with three
//! acceptance tests:
//! - `v2_workflow_erosion_mass_balanced` (β=1.0 → conservation 1e-6)
//! - `v2_workflow_erosion_diffusive` (β=0.0 → mass monotonically
//!   decreases)
//! - `v2_workflow_erosion_applied_everywhere` (interior craton cells
//!   lose mass, not just coastlines — counter-isostasy contract).

use crate::tectonics_v2::field::Field2D;

/// Apply one pass of low-res parametric erosion in-place on `s`.
///
/// `sea_level_reference` is the continental/oceanic threshold used
/// both for the cell-selection mask (`S̃[i] > sea_level_reference`)
/// and for the depth-weighted erosion magnitude `(S̃[i] -
/// sea_level_reference)`. Pass the value returned by
/// `compute_isostasy(s).sea_level_normalized` for the adaptive
/// threshold contract.
///
/// Returns the integrated `Δh` over all eroded cells (`erosion_volume
/// _removed`), used by the workflow metrics dashboard.
///
/// Phase 1 stub: panics on call. Phase 2 ships the implementation.
pub fn apply(_s: &mut Field2D, _alpha: f64, _beta: f64, _sea_level_reference: f64) -> f64 {
    unimplemented!(
        "Step 12 Phase 2 — diffusive erosion algorithm + 3 acceptance tests"
    );
}
