//! Hard clamp `S̃ ≥ S_MIN` + flux tracking for mass balance.
//!
//! After advection + source/sink update, any cell whose updated
//! `S̃` would fall below `S_MIN` is raised to `S_MIN`. The
//! difference `S̃_post - S̃_pre_clamp` is an **artificial mass
//! flux** introduced by the clamp; tracking its integral is what
//! makes the Step 5 mass-balance residual (issue #89 D5) sound
//! when the clamp activates.
//!
//! `S_MIN = 0.05` corresponds to ~1.75 km of crustal thickness
//! (S* = 35 km), the thickness of very young mid-ocean ridge
//! material. Below that the solver's power-law rheology degenerates
//! non-physically; the clamp prevents unbounded thinning.

use super::super::field::Field2D;

/// Minimum admissible `S̃` value (issue #89 D8). Shared constant so
/// tests and report sections stay consistent with the harness.
pub const S_MIN: f64 = 0.05;

/// Per-step clamp statistics. All numbers are summed over the whole
/// field for the step in question.
#[derive(Clone, Copy, Debug, Default)]
pub struct ClampStats {
    /// Number of cells where the clamp fired (`S̃ < S_MIN` before).
    pub activations: usize,
    /// Total number of cells tested.
    pub cells: usize,
    /// Artificial mass flux injected by the clamp at this step:
    /// `Σ_cells (S̃_post - S̃_pre)`. Always ≥ 0 for a floor clamp.
    pub injected_flux: f64,
}

impl ClampStats {
    pub fn activation_fraction(&self) -> f64 {
        if self.cells == 0 {
            0.0
        } else {
            self.activations as f64 / self.cells as f64
        }
    }
}

/// Apply a hard floor `S̃ ≥ S_MIN` in place and return the
/// per-step clamp stats. Callers accumulate these into the run-
/// level [`crate::tectonics_v2::diagnostics::metrics::MassBudget`].
pub fn apply_clamp_with_tracking(s: &mut Field2D) -> ClampStats {
    let nx = s.nx();
    let ny = s.ny();
    let mut activations = 0usize;
    let mut injected = 0.0_f64;
    for cell in s.data_mut().iter_mut() {
        if *cell < S_MIN {
            injected += S_MIN - *cell;
            *cell = S_MIN;
            activations += 1;
        }
    }
    ClampStats {
        activations,
        cells: nx * ny,
        injected_flux: injected,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clamp_leaves_above_floor_cells_untouched() {
        let mut s = Field2D::new(3, 3);
        for v in s.data_mut().iter_mut() {
            *v = 1.0;
        }
        let before: Vec<f64> = s.data().to_vec();
        let stats = apply_clamp_with_tracking(&mut s);
        assert_eq!(stats.activations, 0);
        assert_eq!(stats.injected_flux, 0.0);
        for (a, b) in before.iter().zip(s.data().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn clamp_raises_below_floor_cells() {
        let mut s = Field2D::new(2, 2);
        s.set(0, 0, 0.01);
        s.set(1, 0, 0.05);
        s.set(0, 1, -0.2);
        s.set(1, 1, 1.0);
        let stats = apply_clamp_with_tracking(&mut s);
        // 0.01 and -0.2 below floor → 2 activations; 0.05 is at
        // the floor exactly and should not fire (strict `<`).
        assert_eq!(stats.activations, 2);
        assert_eq!(s.get(0, 0), S_MIN);
        assert_eq!(s.get(1, 0), 0.05);
        assert_eq!(s.get(0, 1), S_MIN);
        assert_eq!(s.get(1, 1), 1.0);
        // Flux = (0.05 - 0.01) + (0.05 - (-0.2)) = 0.04 + 0.25 = 0.29.
        assert!((stats.injected_flux - 0.29).abs() < 1e-12);
    }

    #[test]
    fn activation_fraction_is_count_over_cells() {
        let stats = ClampStats {
            activations: 3,
            cells: 12,
            injected_flux: 0.0,
        };
        assert!((stats.activation_fraction() - 0.25).abs() < 1e-12);
    }
}
