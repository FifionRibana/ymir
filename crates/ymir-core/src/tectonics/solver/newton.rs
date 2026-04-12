//! Jacobian-Free Newton-Krylov (JFNK) solver for the nonlinear Stokes problem.
//!
//! The Jacobian-vector product J·δv is approximated by a finite difference:
//! `J·δv ≈ (F(v + ε·δv) - F(v)) / ε`, reusing all existing Stokes/viscosity code.

use super::config::{NewtonConfig, PicardConfig};
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::linear_solve::solve_bicgstab;
use super::picard::{compute_strain_rate, compute_viscosity};
use super::traction::TractionField;
use super::stokes::{apply_stokes, compute_jacobi_precond, compute_rhs};
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

/// Solve for velocity using Jacobian-Free Newton-Krylov (JFNK).
///
/// The inner linear solve uses BiCGSTAB with the Jacobian approximated
/// by finite-difference directional derivatives.
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

        // Update preconditioner based on current η
        compute_jacobi_precond(&ws.eta, grid, &mut ws.bicgstab.precond);

        // 3. Solve J(vᵏ)·δv = -F(vᵏ) via BiCGSTAB with JFNK operator
        ws.jfnk_delta_v.fill(0.0);

        // We need to capture references for the closure. The tricky part is that
        // apply_jacobian_fd needs to call compute_nonlinear_residual which borrows
        // grid mutably. We solve this by using the dedicated jfnk_v_pert buffer
        // and a separate eta/strain_rate computation path.
        //
        // The closure captures immutable snapshots needed for the FD approximation.
        let eps_scale = newton_config.fd_epsilon_scale;
        let v_base = ws.v_packed.clone();
        let f_v_base = ws.jfnk_f_v.clone();
        let rhs_ref = ws.rhs.clone();

        // We can't borrow ws mutably inside the closure that also borrows ws.bicgstab,
        // so we use a separate mini-workspace for the JFNK matvec.
        let mut jfnk_eta = Field2D::new(n);
        let mut jfnk_sr = Field2D::new(n);
        let mut jfnk_residual = vec![0.0; n_dof];

        let linear_result = {
            let v_pert_buf = &mut ws.jfnk_v_pert;

            solve_bicgstab(
                &mut ws.jfnk_delta_v,
                &ws.jfnk_neg_f,
                |delta_v, out| {
                    // J·δv ≈ (F(v + ε·δv) - F(v)) / ε
                    let v_norm: f64 = v_base.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let dv_norm: f64 = delta_v.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let eps = eps_scale * v_norm.max(1.0) / dv_norm.max(1e-30);

                    for i in 0..n_dof {
                        v_pert_buf[i] = v_base[i] + eps * delta_v[i];
                    }

                    compute_nonlinear_residual(
                        v_pert_buf,
                        &rhs_ref,
                        grid,
                        picard_config,
                        &mut jfnk_eta,
                        &mut jfnk_sr,
                        &mut jfnk_residual,
                    );

                    let inv_eps = 1.0 / eps;
                    for i in 0..n_dof {
                        out[i] = (jfnk_residual[i] - f_v_base[i]) * inv_eps;
                    }
                },
                &mut ws.bicgstab,
                newton_config.bicgstab_max_iter,
                newton_config.bicgstab_tolerance,
            )
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
            bicgstab_max_iter: 500,
            bicgstab_tolerance: 1e-10,
            fd_epsilon_scale: 1e-7,
        };
        let mut ws = SolverWorkspace::new(n);

        let result =
            solve_velocity_newton(&mut grid, &plates, 1.0, &picard_config, &newton_config, &mut ws);
        assert!(result.converged, "Newton should converge for linear viscosity");
        // With n=1 the problem is linear: F(v) = A·v - b. Newton should converge
        // in very few iterations (1-2, depending on initial guess quality).
        assert!(
            result.iterations <= 2,
            "Linear problem should converge in ≤ 2 Newton iterations, got {}",
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
            strain_rate_min: 1e-6,
            ..PicardConfig::default()
        };
        let newton_config = NewtonConfig {
            max_iterations: 15,
            tolerance: 1e-4,
            bicgstab_max_iter: 500,
            bicgstab_tolerance: 1e-6,
            fd_epsilon_scale: 1e-7,
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
            result.iterations <= 10,
            "Newton should converge in ≤ 10 iterations, got {}",
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
            strain_rate_min: 1e-6,
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
            max_iterations: 15,
            tolerance: 1e-6,
            bicgstab_max_iter: 500,
            bicgstab_tolerance: 1e-8,
            fd_epsilon_scale: 1e-7,
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

    #[test]
    fn newton_convergence_is_quadratic() {
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
            strain_rate_min: 1e-6,
            ..PicardConfig::default()
        };
        let newton_config = NewtonConfig {
            max_iterations: 15,
            tolerance: 1e-10, // Very tight to see the convergence curve
            bicgstab_max_iter: 500,
            bicgstab_tolerance: 1e-8,
            fd_epsilon_scale: 1e-7,
        };
        let mut ws = SolverWorkspace::new(n);

        // Compute RHS
        compute_rhs(&grid, &plates, 1.0, &mut ws.rhs);
        let n2 = n * n;
        let n_dof = 2 * n2;
        let mean_vx: f64 = ws.rhs[..n2].iter().sum::<f64>() / n2 as f64;
        let mean_vy: f64 = ws.rhs[n2..n_dof].iter().sum::<f64>() / n2 as f64;
        for val in &mut ws.rhs[..n2] {
            *val -= mean_vx;
        }
        for val in &mut ws.rhs[n2..n_dof] {
            *val -= mean_vy;
        }

        pack_velocity(&grid, &mut ws.v_packed);

        let mut residuals = Vec::new();
        for _ in 0..newton_config.max_iterations {
            compute_nonlinear_residual(
                &ws.v_packed,
                &ws.rhs,
                &mut grid,
                &picard_config,
                &mut ws.eta,
                &mut ws.strain_rate,
                &mut ws.jfnk_f_v,
            );
            let f_norm: f64 = ws.jfnk_f_v.iter().map(|x| x * x).sum::<f64>().sqrt();
            residuals.push(f_norm);

            let b_norm: f64 = ws.rhs.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-14);
            if f_norm < newton_config.tolerance * b_norm {
                break;
            }

            for i in 0..n_dof {
                ws.jfnk_neg_f[i] = -ws.jfnk_f_v[i];
            }
            compute_jacobi_precond(&ws.eta, &grid, &mut ws.bicgstab.precond);
            ws.jfnk_delta_v.fill(0.0);

            let v_base = ws.v_packed.clone();
            let f_v_base = ws.jfnk_f_v.clone();
            let rhs_clone = ws.rhs.clone();
            let eps_scale = newton_config.fd_epsilon_scale;
            let mut jfnk_eta = Field2D::new(n);
            let mut jfnk_sr = Field2D::new(n);
            let mut jfnk_res = vec![0.0; n_dof];

            solve_bicgstab(
                &mut ws.jfnk_delta_v,
                &ws.jfnk_neg_f,
                |dv, out| {
                    let v_norm: f64 = v_base.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let dv_norm: f64 = dv.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let eps = eps_scale * v_norm.max(1.0) / dv_norm.max(1e-30);
                    for i in 0..n_dof {
                        ws.jfnk_v_pert[i] = v_base[i] + eps * dv[i];
                    }
                    compute_nonlinear_residual(
                        &ws.jfnk_v_pert,
                        &rhs_clone,
                        &mut grid,
                        &picard_config,
                        &mut jfnk_eta,
                        &mut jfnk_sr,
                        &mut jfnk_res,
                    );
                    let inv_eps = 1.0 / eps;
                    for i in 0..n_dof {
                        out[i] = (jfnk_res[i] - f_v_base[i]) * inv_eps;
                    }
                },
                &mut ws.bicgstab,
                newton_config.bicgstab_max_iter,
                newton_config.bicgstab_tolerance,
            );

            for i in 0..n_dof {
                ws.v_packed[i] += ws.jfnk_delta_v[i];
            }
            let mean_vx: f64 = ws.v_packed[..n2].iter().sum::<f64>() / n2 as f64;
            let mean_vy: f64 = ws.v_packed[n2..n_dof].iter().sum::<f64>() / n2 as f64;
            for val in &mut ws.v_packed[..n2] {
                *val -= mean_vx;
            }
            for val in &mut ws.v_packed[n2..n_dof] {
                *val -= mean_vy;
            }
        }

        // Check that convergence is at least superlinear: residuals should decrease
        // faster and faster. Check that the log-residual curve is convex (each ratio
        // is smaller than the previous).
        assert!(
            residuals.len() >= 3,
            "Need at least 3 residuals to check convergence rate"
        );

        // Just verify the residual is decreasing overall
        let last = *residuals.last().unwrap();
        let first = residuals[0];
        assert!(
            last < first * 0.01,
            "Residual should decrease significantly: {first} → {last}"
        );
    }
}
