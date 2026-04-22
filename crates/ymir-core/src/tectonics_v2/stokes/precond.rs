//! Preconditioner applications with null-space projection wrapping.
//!
//! The thin-sheet momentum operator is SPD on the zero-mean velocity
//! subspace; a Jacobi (inverse-diagonal) preconditioner is adequate
//! for Step 0's constant-η regime. Wrapping with the 2-D velocity
//! projector before and after `M⁻¹` keeps CG search directions
//! orthogonal to the null space at every iteration.
//!
//! # Case (B) — diagonal supplied, not reconstructed here
//!
//! Step 4 (basal drag) codified the observation that this module is
//! case (B) of the Step 4 spec:
//!
//! - [`VelocityJacobi::from_diagonal`] takes `diag_vx` and `diag_vy`
//!   as **external slices**. No symbolic rewrite of the viscous
//!   stencil happens in this file; the module is a pure consumer.
//! - The analytical reconstruction lives in
//!   [`stokes::operator::momentum_diagonal`][md] — a symbolic rewrite
//!   of [`stokes::operator::apply_momentum`][am]'s stencil. Any new
//!   operator contribution (Step 4 basal drag's `Br · S̃²` diagonal;
//!   future Step 7/8 spike operators) must be added in both
//!   [`apply_momentum`][am] and [`momentum_diagonal`][md] with
//!   matching cell-to-face averaging conventions, or CG's
//!   preconditioner drifts silently from the assembled operator.
//!
//! Therefore basal drag does **not** modify `precond.rs`: it is
//! injected into the diagonal slice at construction time by the
//! caller (usually the solver harness). The
//! `tests/v2_precond_drag_diagonal` integration test pins the
//! consistency between the analytical diagonal and a unit-vector
//! probe of the assembled operator at 1e-14.
//!
//! [md]: crate::tectonics_v2::stokes::operator::momentum_diagonal
//! [am]: crate::tectonics_v2::stokes::operator::apply_momentum

use super::nullspace::project_velocity;

/// Precomputed reciprocal of the momentum-diagonal, one value per
/// velocity DOF. Invariant per sheet solve since η is held fixed
/// during the solve.
pub struct VelocityJacobi {
    inv_diag_vx: Vec<f64>,
    inv_diag_vy: Vec<f64>,
}

impl VelocityJacobi {
    pub fn from_diagonal(diag_vx: &[f64], diag_vy: &[f64], floor: f64) -> Self {
        let inv_diag_vx: Vec<f64> = diag_vx
            .iter()
            .map(|d| {
                let eff = d.abs().max(floor).copysign(d.signum());
                1.0 / eff
            })
            .collect();
        let inv_diag_vy: Vec<f64> = diag_vy
            .iter()
            .map(|d| {
                let eff = d.abs().max(floor).copysign(d.signum());
                1.0 / eff
            })
            .collect();
        Self { inv_diag_vx, inv_diag_vy }
    }

    /// Apply `z = M⁻¹ r` with null-space projection wrapping both ends.
    /// Slices are packed `[vx; vy]` of length `2·n_cells`.
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = self.inv_diag_vx.len();
        debug_assert_eq!(r.len(), 2 * n);
        debug_assert_eq!(z.len(), 2 * n);

        let (r_vx, r_vy) = r.split_at(n);
        let (z_vx, z_vy) = z.split_at_mut(n);
        let mut rx = r_vx.to_vec();
        let mut ry = r_vy.to_vec();
        project_velocity(&mut rx, &mut ry);
        for k in 0..n {
            z_vx[k] = self.inv_diag_vx[k] * rx[k];
            z_vy[k] = self.inv_diag_vy[k] * ry[k];
        }
        project_velocity(z_vx, z_vy);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn velocity_jacobi_produces_zero_mean_output() {
        let n = 64;
        let diag_vx = vec![6.0; n];
        let diag_vy = vec![6.0; n];
        let jac = VelocityJacobi::from_diagonal(&diag_vx, &diag_vy, 1e-20);
        let mut r = vec![0.0; 2 * n];
        for k in 0..n {
            r[k] = 3.0 + (k as f64).sin();
            r[n + k] = -1.0 + ((k as f64) * 0.7).cos();
        }
        let mut z = vec![0.0; 2 * n];
        jac.apply(&r, &mut z);
        let mean_zx: f64 = z[..n].iter().sum::<f64>() / n as f64;
        let mean_zy: f64 = z[n..].iter().sum::<f64>() / n as f64;
        assert!(mean_zx.abs() < 1e-12);
        assert!(mean_zy.abs() < 1e-12);
    }
}
