//! Step 11 — Plate kinematic drift.
//!
//! Steps 0-10 left the velocity field at rest; motion could only emerge
//! from forcing mechanisms (mantle convection, GPE). Step 8.6 Phase 7
//! visual review surfaced two limitations of that contract:
//!
//! 1. Without mantle, the system stayed quiescent regardless of plate
//!    geometry — the user could not configure scenarios driven by
//!    explicit plate motion.
//! 2. Activating mantle to obtain motion forced the system into the
//!    saturated regime documented in Step 8 (CG cap, η-contrast
//!    4×10⁴). There was no way to obtain moderate motion.
//!
//! ## Semantics: drift, not initial condition
//!
//! Step 11 was originally framed as an "initial velocity per plate"
//! mechanism. That framing is wrong for this codebase: the
//! thin-viscous-sheet Stokes solver is **quasi-static** (no inertia
//! term), so any `v(t=0)` is overwritten in full by the first solve
//! — `v` is determined by `forcing(S)` at every step, not by an
//! initial state.
//!
//! Phase 4 reframed the mechanism as a **plate kinematic drift**:
//! `v_total = v_solver + v_drift`, where `v_drift` is the per-plate
//! velocity field (constructed once via [`build`]) added back to the
//! solver output after every Stokes solve. The drift then propagates
//! through advection of `S̃` / age field / etc., so the user-prescribed
//! plate motion is observable in the dynamics.
//!
//! This breaks strict momentum conservation but provides a controllable
//! forcing mechanism orthogonal to mantle convection. See
//! `docs/solver-scaling-step11-patch.md` (§4.12) for the full physical
//! note. The validity régime is `|v_drift| ≤ 1.0` so the drift stays a
//! moderate perturbation around the solver's solution.
//!
//! ## Configuration
//!
//! The mechanism is opt-in via [`PlateKinematicConfig`]:
//! [`PlateKinematicConfig::Zero`] (the default) short-circuits the
//! entire drift pipeline so Steps 0-10 regression scenarios stay
//! bit-identical to their pre-Step-11 baselines.
//! [`PlateKinematicConfig::PerPlate`] carries one `(vx, vy)` per plate
//! (indexed by `plate_id`) and a `boundary_smoothing_width` knob
//! (default 1.5 cells, same family as
//! [`super::init::InitMode::Uniform`]'s width).

pub mod field;
#[cfg(test)]
mod sanity_check;

pub use field::build;

use serde::{Deserialize, Serialize};

/// Step 11 — plate kinematic drift configuration. Stored on the
/// harness [`super::diagnostics::harness::BaselineConfig`].
///
/// Default [`PlateKinematicConfig::Zero`] preserves the Steps 0-10
/// zero-init contract — the harness must structurally short-circuit
/// the entire `plate_kinematic::field` build path when this variant
/// is set, so no allocations or computations are introduced for
/// existing regression scenarios.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PlateKinematicConfig {
    /// All plates start at rest. `vx[i] = vy[i] = 0` for every cell —
    /// the pre-Step-11 contract. The harness must short-circuit the
    /// `plate_kinematic::field::build` path entirely for this
    /// variant (no allocation, no per-cell loop) so the regression
    /// scenarios stay bit-identical to their Step 10 baselines.
    Zero,
    /// Per-plate velocity assignment with smoothstep blending across
    /// inter-plate boundaries.
    ///
    /// `velocities[p]` is the assigned `(vx, vy)` for plate id `p`,
    /// in the same nondimensional units as the dynamic velocity
    /// field (typical scale `O(1)`, range `[-1, 1]` per component
    /// — see issue D1).
    ///
    /// `boundary_smoothing_width` is the half-width (in cells) of
    /// the transition zone across an inter-plate boundary. Cells at
    /// `dist_to_boundary >= width` hold their plate's velocity
    /// exactly; cells closer than `width` blend toward the
    /// neighbouring plate's velocity via cubic smoothstep
    /// (`3t² − 2t³`). Default `1.5` cells — same family as
    /// [`super::init::InitMode::Uniform`]'s width knob.
    ///
    /// Validation contract (enforced at harness entry, Phase 3):
    /// `velocities.len()` must equal the Voronoï plate count.
    PerPlate { velocities: Vec<(f64, f64)>, boundary_smoothing_width: f64 },
}

impl Default for PlateKinematicConfig {
    fn default() -> Self {
        PlateKinematicConfig::Zero
    }
}

impl PlateKinematicConfig {
    /// `true` if the variant is structurally `Zero` — the harness
    /// uses this to take the short-circuit branch (no allocation, no
    /// per-cell init) at run start.
    pub fn is_zero(&self) -> bool {
        matches!(self, PlateKinematicConfig::Zero)
    }

    /// Default smoothing width used when
    /// [`PlateKinematicConfig::PerPlate`] is constructed without an
    /// explicit value (UI presets, panel reset, …). Matches the
    /// width family of [`super::init::InitMode::Uniform`].
    pub const DEFAULT_BOUNDARY_SMOOTHING_WIDTH: f64 = 1.5;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_zero() {
        assert!(matches!(PlateKinematicConfig::default(), PlateKinematicConfig::Zero));
        assert!(PlateKinematicConfig::default().is_zero());
    }

    #[test]
    fn per_plate_is_not_zero() {
        let cfg = PlateKinematicConfig::PerPlate {
            velocities: vec![(0.0, 0.0); 4],
            boundary_smoothing_width: 1.5,
        };
        assert!(
            !cfg.is_zero(),
            "PerPlate is structurally non-zero even with all-zero entries — \
             the variant itself triggers the build path so the regression \
             contract attaches to the Zero variant only"
        );
    }

    /// Round-trip through JSON to validate the panel / preset path
    /// (Step 8.6 v2 presets serialise the full BaselineConfig).
    #[test]
    fn json_roundtrip_zero() {
        let cfg = PlateKinematicConfig::Zero;
        let s = serde_json::to_string(&cfg).unwrap();
        let back: PlateKinematicConfig = serde_json::from_str(&s).unwrap();
        assert_eq!(cfg, back);
    }

    #[test]
    fn json_roundtrip_per_plate() {
        let cfg = PlateKinematicConfig::PerPlate {
            velocities: vec![(0.5, 0.0), (-0.5, 0.0), (0.0, 0.3)],
            boundary_smoothing_width: 2.0,
        };
        let s = serde_json::to_string(&cfg).unwrap();
        let back: PlateKinematicConfig = serde_json::from_str(&s).unwrap();
        assert_eq!(cfg, back);
    }
}
