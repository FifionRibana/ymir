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
    /// Produces visible convergent **and** divergent boundaries
    /// for the default 8-plate Voronoï layout at 64² (see
    /// [`crate::tectonics_v2::voronoi::VoronoiConfig::default`]).
    /// For other plate counts the vector is filled by cycling
    /// through the 8 hand-tuned vectors below — the resulting
    /// configuration retains the convergence/divergence signal
    /// at the cost of less geographic plausibility.
    ///
    /// See `docs/reports/c1_phase_1_1_advection/README.md` for
    /// the visual acceptance criteria this preset targets.
    pub fn preset_phase_1_1(num_plates: usize) -> Self {
        // Base set: 4 cardinal pairs producing two convergent
        // axes (E/W and N/S) and one diagonal divergent pair
        // (NE/SW). The two NW/SE diagonal plates round out the
        // 8-plate default.
        //
        // Magnitudes ~ 0.01 (non-dim length units per non-dim
        // time). At 64² with dx = 1/64 ≈ 0.0156, max|v| ≈ 0.0113
        // gives Δt_CFL = 0.5·dx/max|v| ≈ 0.69 non-dim/step; a
        // 300-step run advects ~200 non-dim time units which is
        // enough to move every plate ~2-3 cells worth — visible
        // signal without aliasing.
        let base: [(f64, f64); 8] = [
            (0.01, 0.00),    // plate 0 — east
            (-0.01, 0.00),   // plate 1 — west (converges with 0 on E/W axis)
            (0.00, 0.01),    // plate 2 — north
            (0.00, -0.01),   // plate 3 — south (converges with 2 on N/S axis)
            (0.008, 0.008),  // plate 4 — NE
            (-0.008, -0.008), // plate 5 — SW (diverges from 4)
            (0.005, -0.005), // plate 6 — SE
            (-0.005, 0.005), // plate 7 — NW
        ];

        let velocities = (0..num_plates).map(|p| base[p % base.len()]).collect();
        Self { velocities }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preset_default_8_plates_has_visible_convergence_pairs() {
        let k = PlateKinematics::preset_phase_1_1(8);
        assert_eq!(k.num_plates(), 8);
        // Plate 0 east-bound, plate 1 west-bound → convergent
        // (their relative velocity points into each other on
        // the E/W axis).
        assert!(k.velocities[0].0 > 0.0);
        assert!(k.velocities[1].0 < 0.0);
        // Plate 2/3 same on N/S axis.
        assert!(k.velocities[2].1 > 0.0);
        assert!(k.velocities[3].1 < 0.0);
    }

    #[test]
    fn preset_smaller_plate_count_cycles_base_pattern() {
        let k = PlateKinematics::preset_phase_1_1(4);
        assert_eq!(k.num_plates(), 4);
        assert_eq!(k.velocities[0], (0.01, 0.00));
        assert_eq!(k.velocities[3], (0.00, -0.01));
    }

    #[test]
    fn preset_larger_plate_count_cycles_base_pattern() {
        let k = PlateKinematics::preset_phase_1_1(12);
        assert_eq!(k.num_plates(), 12);
        assert_eq!(k.velocities[8], k.velocities[0]);
        assert_eq!(k.velocities[11], k.velocities[3]);
    }

    #[test]
    fn max_velocity_matches_largest_magnitude() {
        let k = PlateKinematics::preset_phase_1_1(8);
        let max = k.max_velocity();
        // Largest base entry is the diagonal pair (±0.008, ±0.008)
        // with magnitude √(2)·0.008 ≈ 0.01131, edging out the
        // cardinal pairs at magnitude 0.01. This is intentional —
        // see the `preset_phase_1_1` doc comment for the Δt_CFL
        // sizing rationale.
        let expected = 0.008_f64 * std::f64::consts::SQRT_2;
        assert!((max - expected).abs() < 1e-12, "max_velocity = {}", max);
    }

    #[test]
    fn max_velocity_is_zero_for_empty_kinematics() {
        let k = PlateKinematics { velocities: Vec::new() };
        assert_eq!(k.max_velocity(), 0.0);
    }
}
