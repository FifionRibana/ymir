//! Preconditioner applications with null-space projection wrapping.
//!
//! Every `precond` closure passed to the CG solver projects out the
//! null-space components **before and after** the pointwise
//! `M⁻¹` step. Two array-wide means per application is O(N) and
//! negligible against the stencil cost.
//!
//! Block-diagonal structure (per Step 0 spec):
//! - velocity block: inverse of the assembled momentum diagonal
//!   (Jacobi). Wrapped with 2-D velocity null-space projection.
//! - pressure block: `diag(1/η)` viscosity-scaled mass matrix.
//!   Wrapped with 1-D pressure null-space projection.
//!
//! For constant η both diagonals are spatially uniform; the projection
//! is what prevents null-space contamination of Krylov search
//! directions.

use super::super::field::Field2D;
use super::nullspace::{project_pressure, project_velocity, subtract_mean};

/// Precomputed reciprocal of the momentum-diagonal, one value per
/// velocity DOF. Invariant per Stokes solve since η is held fixed
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

    /// Apply `z = M_v⁻¹ r` with null-space projection wrapping both ends.
    /// Slices are packed `[vx; vy]` of length `2·n_cells`.
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = self.inv_diag_vx.len();
        debug_assert_eq!(r.len(), 2 * n);
        debug_assert_eq!(z.len(), 2 * n);

        // Project the input so that the null-space is never multiplied
        // by diag⁻¹, which for large diagonal entries would otherwise
        // leave a noise floor at the level of the residual's null-space
        // projection (typically 10⁻¹² but nonzero).
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

/// Pressure-mass preconditioner `M_p = diag(η)`, i.e. `M_p⁻¹ r = r / η`.
/// For constant η this is identity scaled by `1/η`; the mean-projection
/// wrapping is what makes the preconditioner consistent with the
/// pressure null space.
pub struct PressureMass {
    inv_eta: Vec<f64>,
}

impl PressureMass {
    pub fn from_eta(eta: &Field2D, floor: f64) -> Self {
        let inv_eta: Vec<f64> = eta
            .data()
            .iter()
            .map(|&e| 1.0 / e.abs().max(floor))
            .collect();
        Self { inv_eta }
    }

    /// Apply `z = M_p⁻¹ r` with pressure null-space projection.
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = self.inv_eta.len();
        debug_assert_eq!(r.len(), n);
        debug_assert_eq!(z.len(), n);
        let mut rp = r.to_vec();
        project_pressure(&mut rp);
        for k in 0..n {
            z[k] = self.inv_eta[k] * rp[k];
        }
        project_pressure(z);
    }
}

/// Convenience: project pressure in place and discard its mean.
pub fn clean_pressure(p: &mut [f64]) {
    subtract_mean(p);
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
        // RHS with nonzero mean in both components.
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

    #[test]
    fn pressure_mass_produces_zero_mean_output() {
        let nx = 8;
        let ny = 8;
        let eta = Field2D::filled(nx, ny, 2.0);
        let pc = PressureMass::from_eta(&eta, 1e-20);
        let mut r = vec![0.0; nx * ny];
        for k in 0..r.len() {
            r[k] = 0.5 + (k as f64).sin();
        }
        let mut z = vec![0.0; nx * ny];
        pc.apply(&r, &mut z);
        let mean: f64 = z.iter().sum::<f64>() / z.len() as f64;
        assert!(mean.abs() < 1e-12);
    }
}
