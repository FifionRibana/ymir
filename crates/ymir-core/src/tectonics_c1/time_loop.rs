//! C1 time loop — forward-Euler advection of `S̃` and `age`.
//!
//! ## Phase 1.1 contract (filled in Stage 2)
//!
//! Advection only, no closures. At each step:
//!
//! 1. Compute `Δt = 0.5 · min(dx, dy) / max(|v|, 1e-12)` (CFL,
//!    safety factor 0.5 from §4.6 of the design doc).
//! 2. Build per-cell velocity field `(vx, vy)` from per-plate
//!    kinematics — every cell of plate `p` gets `velocities[p]`
//!    (cratonic and non-cratonic cells alike at this stage).
//! 3. Advect `S̃` and `age` independently using v2's
//!    [`crate::tectonics_v2::advection::step_upwind`] scheme.
//! 4. Invoke the caller-supplied `on_step` callback for
//!    diagnostics.
//!
//! Closures join in Phase 1.2+ between steps 3 and 4.
//!
//! ## Stage 1 status
//!
//! Signature only. Body lands in Stage 2 together with the
//! kinematics preset, since the two must be co-tested.

use super::kinematics::PlateKinematics;
use super::state::C1State;

/// Tunables for [`run_advection_only`].
///
/// `dx` and `dy` are the cell-size in non-dimensional length units.
/// For the typical unit-domain `1×1` non-dim setup they both equal
/// `1.0 / grid_size`.
#[derive(Clone, Copy, Debug)]
pub struct C1TimeLoopConfig {
    pub n_steps: usize,
    pub dx: f64,
    pub dy: f64,
}

/// Run advection-only forward in time. The callback fires once
/// per step **after** the state has been updated, with the
/// 0-based step index already incremented (so first invocation
/// is `on_step(0, …)` representing "post-step-1 state").
///
/// Stage 1: signature pinned. Stage 2 fills the body.
pub fn run_advection_only<F>(
    _state: &mut C1State,
    _kinematics: &PlateKinematics,
    _config: &C1TimeLoopConfig,
    _on_step: F,
) where
    F: FnMut(usize, &C1State),
{
    todo!("Stage 2 — CFL Δt, velocity field assembly, two step_upwind calls, callback fire")
}
