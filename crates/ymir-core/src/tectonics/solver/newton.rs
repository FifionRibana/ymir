//! Quasi-Newton solver for the nonlinear Stokes problem.
//!
//! At each Newton iteration, the viscosity η is frozen at the current iterate
//! and the linearized Stokes operator A(η_frozen) is used as the Jacobian
//! approximation. This operator is SPD, so CG is used instead of BiCGSTAB.

use super::config::{NewtonConfig, PicardConfig, Preconditioner};
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::linear_solve::{apply_jacobi, solve_cg};
use super::picard::{compute_strain_rate, compute_viscosity};
use super::stokes::{apply_ssor, apply_stokes, compute_jacobi_precond, compute_rhs, StencilCoeffs};
use super::traction::TractionField;
use super::workspace::{pack_velocity, unpack_velocity, SolverWorkspace};

/// Result of a Newton solve.
pub struct NewtonResult {
    pub converged: bool,
    pub iterations: usize,
    pub final_residual: f64,
    pub total_linear_iterations: usize,
}

/// Compute the nonlinear residual F(v) = A(η(v))·v - b.
///
/// Side effects: updates `eta_out` and `strain_rate_out` from the velocity field.
fn compute_nonlinear_residual(
    v_packed: &[f64],
    b: &[f64],
    grid: &mut StaggeredGrid,
    picard_config: &PicardConfig,
    eta_out: &mut Field2D,
    strain_rate_out: &mut Field2D,
    residual: &mut [f64],
) {
    unpack_velocity(v_packed, grid);
    compute_strain_rate(grid, strain_rate_out);
    compute_viscosity(
        strain_rate_out,
        picard_config.power_law_n,
        picard_config.strain_rate_min,
        picard_config.eta_min,
        picard_config.eta_max,
        eta_out,
    );
    apply_stokes(v_packed, eta_out, grid, residual);
    for i in 0..residual.len() {
        residual[i] -= b[i];
    }
}

/// Solve for velocity using quasi-Newton with frozen viscosity.
///
/// At each Newton step, η is computed from the current v, then frozen.
/// The inner linear system A(η_frozen)·δv = -F(v) is solved by CG (SPD operator).
pub fn solve_velocity_newton(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    gravity_factor: f64,
    picard_config: &PicardConfig,
    newton_config: &NewtonConfig,
    ws: &mut SolverWorkspace,
) -> NewtonResult {
    let n = grid.n;
    let n2 = n * n;
    let n_dof = 2 * n2;
    let mut total_linear = 0usize;

    // Compute RHS (constant during Newton)
    compute_rhs(grid, plates, gravity_factor, &mut ws.rhs);

    // Project out null space
    let mean_vx: f64 = ws.rhs[..n2].iter().sum::<f64>() / n2 as f64;
    let mean_vy: f64 = ws.rhs[n2..n_dof].iter().sum::<f64>() / n2 as f64;
    for val in &mut ws.rhs[..n2] {
        *val -= mean_vx;
    }
    for val in &mut ws.rhs[n2..n_dof] {
        *val -= mean_vy;
    }

    // Pack current velocity as initial guess
    pack_velocity(grid, &mut ws.v_packed);

    for k in 0..newton_config.max_iterations {
        // 1. Compute F(vᵏ) → jfnk_f_v, also updates eta and strain_rate
        compute_nonlinear_residual(
            &ws.v_packed,
            &ws.rhs,
            grid,
            picard_config,
            &mut ws.eta,
            &mut ws.strain_rate,
            &mut ws.jfnk_f_v,
        );

        let f_norm: f64 = ws.jfnk_f_v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let b_norm: f64 = ws.rhs.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-14);

        if f_norm < newton_config.tolerance * b_norm {
            unpack_velocity(&ws.v_packed, grid);
            return NewtonResult {
                converged: true,
                iterations: k,
                final_residual: f_norm / b_norm,
                total_linear_iterations: total_linear,
            };
        }

        // 2. Set up -F(v) as RHS for the linear solve
        for i in 0..n_dof {
            ws.jfnk_neg_f[i] = -ws.jfnk_f_v[i];
        }

        // 3. Solve A(η_frozen)·δv = -F(vᵏ) via CG (SPD operator)
        ws.jfnk_delta_v.fill(0.0);

        // Build preconditioner from frozen η
        let eta_ref = &ws.eta;
        let grid_ref = &*grid;

        let linear_result = match newton_config.preconditioner {
            Preconditioner::Jacobi => {
                compute_jacobi_precond(eta_ref, grid_ref, &mut ws.jacobi_precond);
                let precond_ref = &ws.jacobi_precond;
                solve_cg(
                    &mut ws.jfnk_delta_v,
                    &ws.jfnk_neg_f,
                    |dv, out| apply_stokes(dv, eta_ref, grid_ref, out),
                    |r, z| apply_jacobi(precond_ref, r, z),
                    &mut ws.cg,
                    newton_config.cg_max_iter,
                    newton_config.cg_tolerance,
                )
            }
            Preconditioner::Ssor { omega } => {
                let stencil = StencilCoeffs::compute(eta_ref, grid_ref);
                solve_cg(
                    &mut ws.jfnk_delta_v,
                    &ws.jfnk_neg_f,
                    |dv, out| apply_stokes(dv, eta_ref, grid_ref, out),
                    |r, z| apply_ssor(r, &stencil, n, omega, z),
                    &mut ws.cg,
                    newton_config.cg_max_iter,
                    newton_config.cg_tolerance,
                )
            }
        };
        total_linear += linear_result.iterations;

        // 4. Update: vᵏ⁺¹ = vᵏ + δv
        for i in 0..n_dof {
            ws.v_packed[i] += ws.jfnk_delta_v[i];
        }

        // Project out null space from solution
        let mean_vx: f64 = ws.v_packed[..n2].iter().sum::<f64>() / n2 as f64;
        let mean_vy: f64 = ws.v_packed[n2..n_dof].iter().sum::<f64>() / n2 as f64;
        for val in &mut ws.v_packed[..n2] {
            *val -= mean_vx;
        }
        for val in &mut ws.v_packed[n2..n_dof] {
            *val -= mean_vy;
        }
    }

    unpack_velocity(&ws.v_packed, grid);
    NewtonResult {
        converged: false,
        iterations: newton_config.max_iterations,
        final_residual: f64::NAN,
        total_linear_iterations: total_linear,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::solver::grid::StaggeredGrid;
    use crate::tectonics::solver::traction::TractionField;

    #[test]
    fn newton_linear_viscosity_converges_in_one() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                grid.s.set(i, j, 1.0 + 0.3 * (2.0 * std::f64::consts::PI * x).sin());
            }
        }

        let plates = TractionField::uniform(n, 0.1, 0.0);
        let picard_config = PicardConfig {
            power_law_n: 1.0,
            strain_rate_min: 1e-6,
            ..PicardConfig::default()
        };
        let newton_config = NewtonConfig {
            max_iterations: 15,
            tolerance: 1e-8,
            ..NewtonConfig::default()
        };
        let mut ws = SolverWorkspace::new(n);

        let result =
            solve_velocity_newton(&mut grid, &plates, 1.0, &picard_config, &newton_config, &mut ws);
        assert!(result.converged, "Newton should converge for linear viscosity");
        assert!(
            result.iterations <= 2,
            "Linear problem should converge in <= 2 Newton iterations, got {}",
            result.iterations
        );
    }

    #[test]
    fn newton_power_law_converges() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let plates = TractionField::two_plates_convergent(n, 0.5);
        let picard_config = PicardConfig {
            power_law_n: 3.0,
            strain_rate_min: 1e-3,
            ..PicardConfig::default()
        };
        let newton_config = NewtonConfig {
            max_iterations: 30,
            tolerance: 1e-4,
            ..NewtonConfig::default()
        };
        let mut ws = SolverWorkspace::new(n);

        let result =
            solve_velocity_newton(&mut grid, &plates, 1.0, &picard_config, &newton_config, &mut ws);
        assert!(
            result.converged,
            "Newton should converge for power-law, got {} iterations",
            result.iterations
        );
        assert!(
            result.iterations <= 30,
            "Newton should converge in <= 30 iterations, got {}",
            result.iterations
        );
    }

    #[test]
    fn newton_and_picard_agree() {
        use crate::tectonics::solver::picard::solve_velocity_picard;

        let n = 16;
        let dx = 1.0 / n as f64;
        let n_dof = 2 * n * n;

        let picard_config = PicardConfig {
            max_iterations: 50,
            tolerance: 1e-6,
            relaxation: 0.7,
            cg_max_iter: 500,
            cg_tolerance: 1e-10,
            strain_rate_min: 1e-3,
            power_law_n: 3.0,
            ..PicardConfig::default()
        };

        // Run Picard
        let mut grid_p = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_p.s.set(i, j, 1.0);
            }
        }
        let plates = TractionField::two_plates_convergent(n, 0.5);
        let mut ws_p = SolverWorkspace::new(n);
        let picard_result =
            solve_velocity_picard(&mut grid_p, &plates, 1.0, &picard_config, &mut ws_p);
        assert!(picard_result.converged, "Picard should converge");
        let mut v_picard = vec![0.0; n_dof];
        pack_velocity(&grid_p, &mut v_picard);

        // Run Newton
        let mut grid_n = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_n.s.set(i, j, 1.0);
            }
        }
        let newton_config = NewtonConfig {
            max_iterations: 30,
            tolerance: 1e-6,
            ..NewtonConfig::default()
        };
        let mut ws_n = SolverWorkspace::new(n);
        let newton_result = solve_velocity_newton(
            &mut grid_n,
            &plates,
            1.0,
            &picard_config,
            &newton_config,
            &mut ws_n,
        );
        assert!(newton_result.converged, "Newton should converge");
        let mut v_newton = vec![0.0; n_dof];
        pack_velocity(&grid_n, &mut v_newton);

        // Compare velocity fields
        let mut diff_sq = 0.0;
        let mut norm_sq = 0.0;
        for i in 0..n_dof {
            diff_sq += (v_picard[i] - v_newton[i]).powi(2);
            norm_sq += v_picard[i].powi(2);
        }
        let rel_err = (diff_sq / norm_sq.max(1e-30)).sqrt();
        assert!(
            rel_err < 1e-3,
            "Picard and Newton should agree: rel_err = {rel_err}"
        );
    }
}
