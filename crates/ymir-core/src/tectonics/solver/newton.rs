//! Jacobian-Free Newton-Krylov (JFNK) solver for the nonlinear Stokes problem.
//!
//! The Jacobian-vector product J·δv is approximated by finite difference:
//! `J·δv ≈ (F(v + ε·δv) - F(v)) / ε`, reusing all existing Stokes/viscosity code.
//! The JFNK operator is NOT SPD, so BiCGSTAB is used as the inner linear solver.
//! The preconditioner (Jacobi or SSOR) is built from A(η_frozen) — the Stokes
//! operator linearized at the current viscosity.

use tracing::{debug, warn};

use super::config::YieldingConfig;
use super::config::{NewtonConfig, PicardConfig, Preconditioner};
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::linear_solve::{apply_jacobi, solve_bicgstab};
use super::picard::{apply_eta_multiplier, apply_yielding, compute_strain_rate, compute_viscosity};
use super::stokes::{StencilCoeffs, apply_ssor, apply_stokes, compute_jacobi_precond, compute_rhs};
use super::traction::TractionField;
use super::workspace::{SolverWorkspace, pack_velocity, unpack_velocity};

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
#[allow(clippy::too_many_arguments)]
fn compute_nonlinear_residual(
    v_packed: &[f64],
    b: &[f64],
    grid: &mut StaggeredGrid,
    picard_config: &PicardConfig,
    eta_multiplier: &Field2D,
    plastic_strain: &Field2D,
    yielding: &YieldingConfig,
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
    apply_eta_multiplier(eta_multiplier, picard_config.eta_max, eta_out);
    apply_yielding(strain_rate_out, plastic_strain, yielding, eta_out);
    for val in eta_out.data_mut().iter_mut() {
        *val = val.clamp(picard_config.eta_min, picard_config.eta_max);
    }
    apply_stokes(v_packed, eta_out, grid, residual);
    for i in 0..residual.len() {
        residual[i] -= b[i];
    }
}

/// Solve for velocity using JFNK with inexact Newton and SSOR preconditioning.
///
/// At each Newton step:
/// 1. Compute nonlinear residual F(vᵏ) (updates η)
/// 2. Build preconditioner from A(η_frozen)
/// 3. Solve J·δv = -F(vᵏ) via BiCGSTAB with JFNK operator
/// 4. Update v ← v + δv
#[allow(clippy::too_many_arguments)]
pub fn solve_velocity_newton(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    gravity_factor: f64,
    rho_continental: f64,
    rho_mantle: f64,
    picard_config: &PicardConfig,
    yielding: &YieldingConfig,
    newton_config: &NewtonConfig,
    ws: &mut SolverWorkspace,
) -> NewtonResult {
    let nx = grid.nx();
    let ny = grid.ny();
    let n2 = nx * ny;
    let n_dof = 2 * n2;
    let mut total_linear = 0usize;

    // Compute RHS (constant during Newton)
    compute_rhs(grid, plates, gravity_factor, rho_continental, rho_mantle, &mut ws.rhs);

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

    // Clone eta_multiplier so we can pass grid mutably while reading the multiplier
    let eta_mult = grid.eta_multiplier.clone();
    let ps_snap = grid.plastic_strain.clone();

    let mut prev_f_norm = f64::MAX;

    for k in 0..newton_config.max_iterations {
        // 1. Compute F(vᵏ) → jfnk_f_v, also updates eta and strain_rate
        compute_nonlinear_residual(
            &ws.v_packed,
            &ws.rhs,
            grid,
            picard_config,
            &eta_mult,
            &ps_snap,
            yielding,
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

        // Inexact Newton: adapt inner tolerance to Newton progress
        let mut linear_tol = if newton_config.inexact {
            let ratio = f_norm / prev_f_norm.max(1e-30);
            // Eisenstat-Walker choice 2 (simplified)
            let adaptive = (0.9_f64).min(0.5 * ratio);
            adaptive.max(newton_config.tolerance * 0.1)
        } else {
            newton_config.cg_tolerance
        };
        linear_tol = linear_tol.min(0.1); // never allow the linear solve to be sloppier than 10%
        prev_f_norm = f_norm;

        // 3. Solve J(vᵏ)·δv = -F(vᵏ) via BiCGSTAB with JFNK operator
        ws.jfnk_delta_v.fill(0.0);

        // Snapshot current state for the finite-difference closure
        let eps_scale = newton_config.fd_epsilon_scale;
        let v_base = ws.v_packed.clone();
        let f_v_base = ws.jfnk_f_v.clone();
        let rhs_ref = ws.rhs.clone();

        // Separate mini-workspace for the JFNK matvec (avoids borrow conflicts
        // with ws.bicgstab which is also borrowed mutably by solve_bicgstab)
        let mut jfnk_eta = Field2D::new(nx, ny);
        let mut jfnk_sr = Field2D::new(nx, ny);
        let mut jfnk_residual = vec![0.0; n_dof];

        // Build preconditioner from frozen η (computed once per Newton step)
        let eta_ref = &ws.eta;
        let grid_ref = &*grid;

        let linear_result = match newton_config.preconditioner {
            Preconditioner::Jacobi => {
                compute_jacobi_precond(eta_ref, grid_ref, &mut ws.jacobi_precond);
                let precond_ref = &ws.jacobi_precond;
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
                            &eta_mult,
                            &ps_snap,
                            yielding,
                            &mut jfnk_eta,
                            &mut jfnk_sr,
                            &mut jfnk_residual,
                        );

                        let inv_eps = 1.0 / eps;
                        for i in 0..n_dof {
                            out[i] = (jfnk_residual[i] - f_v_base[i]) * inv_eps;
                        }
                    },
                    |r, z| apply_jacobi(precond_ref, r, z),
                    &mut ws.bicgstab,
                    newton_config.cg_max_iter,
                    linear_tol,
                )
            }
            Preconditioner::Ssor { omega } => {
                let stencil = StencilCoeffs::compute(eta_ref, grid_ref);
                let v_pert_buf = &mut ws.jfnk_v_pert;

                solve_bicgstab(
                    &mut ws.jfnk_delta_v,
                    &ws.jfnk_neg_f,
                    |delta_v, out| {
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
                            &eta_mult,
                            &ps_snap,
                            yielding,
                            &mut jfnk_eta,
                            &mut jfnk_sr,
                            &mut jfnk_residual,
                        );

                        let inv_eps = 1.0 / eps;
                        for i in 0..n_dof {
                            out[i] = (jfnk_residual[i] - f_v_base[i]) * inv_eps;
                        }
                    },
                    |r, z| apply_ssor(r, &stencil, nx, ny, omega, z),
                    &mut ws.bicgstab,
                    newton_config.cg_max_iter,
                    linear_tol,
                )
            }
        };
        total_linear += linear_result.iterations;

        debug!(
            newton_iter = k,
            f_norm,
            rel_residual = f_norm / b_norm,
            linear_iters = linear_result.iterations,
            linear_converged = linear_result.converged,
            linear_tol,
            "newton iteration"
        );

        if !linear_result.converged {
            warn!(
                newton_iter = k,
                linear_iters = linear_result.iterations,
                residual = linear_result.residual_norm,
                "linear solver did not converge within max iterations"
            );
        }

        // 4. Update: vᵏ⁺¹ = vᵏ + α·δv with backtracking line search.
        // Ensure the residual decreases by halving α up to 5 times.
        let v_old = ws.v_packed.clone();
        let mut alpha = 1.0_f64;
        let mut trial_eta = Field2D::new(nx, ny);
        let mut trial_sr = Field2D::new(nx, ny);
        let mut trial_residual = vec![0.0; n_dof];
        let max_backtracks = 5usize;
        let mut final_alpha = alpha;

        for _bt in 0..=max_backtracks {
            for i in 0..n_dof {
                ws.v_packed[i] = v_old[i] + alpha * ws.jfnk_delta_v[i];
            }

            compute_nonlinear_residual(
                &ws.v_packed,
                &ws.rhs,
                grid,
                picard_config,
                &eta_mult,
                &ps_snap,
                yielding,
                &mut trial_eta,
                &mut trial_sr,
                &mut trial_residual,
            );
            let f_trial: f64 = trial_residual.iter().map(|x| x * x).sum::<f64>().sqrt();

            final_alpha = alpha;
            if f_trial < f_norm || alpha <= 0.0625 {
                // Accept: residual decreased, or we've backtracked far enough.
                break;
            }
            alpha *= 0.5;
        }

        if final_alpha < 1.0 {
            debug!(newton_iter = k, alpha = final_alpha, "newton line search backtracked");
        }

        // 5. Project out null space from solution
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
        let mut grid = StaggeredGrid::new(n, n, dx);

        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                grid.s.set(i, j, 1.0 + 0.3 * (2.0 * std::f64::consts::PI * x).sin());
            }
        }

        let plates = TractionField::uniform(n, n, 0.1, 0.0);
        let picard_config =
            PicardConfig { power_law_n: 1.0, strain_rate_min: 1e-6, ..PicardConfig::default() };
        let newton_config =
            NewtonConfig { max_iterations: 15, tolerance: 1e-8, ..NewtonConfig::default() };
        let mut ws = SolverWorkspace::new(n, n);

        let result = solve_velocity_newton(
            &mut grid,
            &plates,
            1.0,
            0.0,
            0.0,
            &picard_config,
            &Default::default(),
            &newton_config,
            &mut ws,
        );
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
        let mut grid = StaggeredGrid::new(n, n, dx);

        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let plates = TractionField::two_plates_convergent(n, n, 0.5);
        let picard_config =
            PicardConfig { power_law_n: 3.0, strain_rate_min: 1e-3, ..PicardConfig::default() };
        let newton_config =
            NewtonConfig { max_iterations: 30, tolerance: 1e-4, ..NewtonConfig::default() };
        let mut ws = SolverWorkspace::new(n, n);

        let result = solve_velocity_newton(
            &mut grid,
            &plates,
            1.0,
            0.0,
            0.0,
            &picard_config,
            &Default::default(),
            &newton_config,
            &mut ws,
        );
        assert!(
            result.converged,
            "Newton should converge for power-law, got {} iterations",
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
        let mut grid_p = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_p.s.set(i, j, 1.0);
            }
        }
        let plates = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ws_p = SolverWorkspace::new(n, n);
        let picard_result = solve_velocity_picard(
            &mut grid_p,
            &plates,
            1.0,
            0.0,
            0.0,
            &picard_config,
            &Default::default(),
            &mut ws_p,
        );
        assert!(picard_result.converged, "Picard should converge");
        let mut v_picard = vec![0.0; n_dof];
        pack_velocity(&grid_p, &mut v_picard);

        // Run Newton
        let mut grid_n = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_n.s.set(i, j, 1.0);
            }
        }
        let newton_config =
            NewtonConfig { max_iterations: 30, tolerance: 1e-6, ..NewtonConfig::default() };
        let mut ws_n = SolverWorkspace::new(n, n);
        let newton_result = solve_velocity_newton(
            &mut grid_n,
            &plates,
            1.0,
            0.0,
            0.0,
            &picard_config,
            &Default::default(),
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
        assert!(rel_err < 1e-2, "Picard and Newton should agree: rel_err = {rel_err}");
    }

    #[test]
    fn inexact_newton_uses_fewer_total_linear_iterations() {
        let n = 16;
        let dx = 1.0 / n as f64;

        let picard_config =
            PicardConfig { power_law_n: 3.0, strain_rate_min: 1e-3, ..PicardConfig::default() };
        let plates = TractionField::two_plates_convergent(n, n, 0.5);

        // Exact Newton (tight inner tolerance)
        let mut grid_exact = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_exact.s.set(i, j, 1.0);
            }
        }
        let config_exact = NewtonConfig {
            max_iterations: 30,
            tolerance: 5e-2,
            inexact: false,
            ..NewtonConfig::default()
        };
        let mut ws_exact = SolverWorkspace::new(n, n);
        let r_exact = solve_velocity_newton(
            &mut grid_exact,
            &plates,
            1.0,
            0.0,
            0.0,
            &picard_config,
            &Default::default(),
            &config_exact,
            &mut ws_exact,
        );

        // Inexact Newton
        let mut grid_inexact = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_inexact.s.set(i, j, 1.0);
            }
        }
        let config_inexact = NewtonConfig {
            max_iterations: 30,
            tolerance: 5e-2,
            inexact: true,
            ..NewtonConfig::default()
        };
        let mut ws_inexact = SolverWorkspace::new(n, n);
        let r_inexact = solve_velocity_newton(
            &mut grid_inexact,
            &plates,
            1.0,
            0.0,
            0.0,
            &picard_config,
            &Default::default(),
            &config_inexact,
            &mut ws_inexact,
        );

        assert!(r_exact.converged, "Exact Newton should converge");
        assert!(r_inexact.converged, "Inexact Newton should converge");

        // Inexact should use fewer total linear iterations
        assert!(
            r_inexact.total_linear_iterations <= r_exact.total_linear_iterations,
            "Inexact should be cheaper or equal: {} vs {} linear iters",
            r_inexact.total_linear_iterations,
            r_exact.total_linear_iterations
        );
    }
}
