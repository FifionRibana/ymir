//! Per-plate kinematic field for C1.
//!
//! ## Phase 1.1 contract
//!
//! Constant-per-plate translation velocity. Every cell belonging
//! to plate `p` sees velocity `velocities[p]` exactly, with no
//! smoothing across plate boundaries (the discontinuity at
//! boundaries is the convergence/divergence signal that future
//! source terms will key off — see §4.4 of the design doc).
//!
//! ## Out of scope this phase
//!
//! - Boundary smoothing (§4.4): deferred to a later phase once
//!   a closure needs a continuous gradient at boundaries.
//! - R7-generalised kinematics sampling (§6.3): boundary
//!   displacement, continental clustering, constrained kinematics
//!   land in Phase 2.
//!
//! ## Stage 1 status (Issue #120)
//!
//! This file ships in Stage 1 as a signature-only skeleton.
//! [`PlateKinematics::preset_phase_1_1`] and any callers in the
//! time loop are filled in Stage 2. [`PlateKinematics::max_velocity`]
//! is implemented now because the time-loop signature in Stage 1
//! cannot reference an undefined helper.

/// Bundle of constant per-plate translation velocities. Index
/// into the vec is the plate id (`u16` upstream, `usize` here for
/// ergonomic indexing).
#[derive(Clone, Debug)]
pub struct PlateKinematics {
    pub velocities: Vec<(f64, f64)>,
}

impl PlateKinematics {
    /// Hand-tuned preset for Phase 1.1.
    ///
    /// Stage 2 fills this in. The Stage 2 implementation must
    /// produce visible convergent **and** divergent boundaries
    /// at the default 8-plate Voronoï layout (W4 of Issue #120).
    pub fn preset_phase_1_1(_num_plates: usize) -> Self {
        todo!("Stage 2 — hand-tuned preset producing convergence + divergence at 8-plate default")
    }

    /// Largest velocity magnitude `sqrt(vx² + vy²)` across all
    /// plates. Used for the CFL Δt calculation in the time loop.
    /// Returns `0.0` for an empty plate set; the caller guards
    /// against division by zero with a small floor.
    pub fn max_velocity(&self) -> f64 {
        self.velocities
            .iter()
            .map(|(vx, vy)| (vx * vx + vy * vy).sqrt())
            .fold(0.0, f64::max)
    }

    pub fn num_plates(&self) -> usize {
        self.velocities.len()
    }
}
