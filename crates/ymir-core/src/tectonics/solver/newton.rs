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

/// Outcome of a Newton solve, with semantic distinction between success and
/// failure modes. Downstream code (e.g. adaptive dt sub-stepping) uses this
/// to choose appropriate recovery strategies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonOutcome {
    /// Converged on the residual criterion: f_norm < tolerance * b_norm.
    ConvergedOnResidual,
    /// Converged on the state criterion: residual stagnated with descending
    /// trend and the velocity increment is below the state tolerance.
    /// The residual may still be slightly above the residual tolerance,
    /// but the physical state is effectively stable.
    ConvergedOnState,
    /// Failed by stagnation: the residual stopped decreasing without
    /// reaching tolerance, but did not oscillate or diverge.
    Stagnation,
    /// Failed by oscillation: detected via cosine of consecutive Newton
    /// steps below threshold for two iterations.
    Oscillation,
    /// Failed by divergence: f_norm increased significantly over recent
    /// iterations.
    Divergence,
    /// Failed by exhausting max_iterations without classification.
    MaxIterations,
}

/// Result of a Newton solve.
pub struct NewtonResult {
    pub outcome: NewtonOutcome,
    pub iterations: usize,
    pub final_residual: f64,
    pub total_linear_iterations: usize,
}

impl NewtonResult {
    /// Returns true for any successful convergence (residual or state).
    pub fn is_converged(&self) -> bool {
        matches!(
            self.outcome,
            NewtonOutcome::ConvergedOnResidual | NewtonOutcome::ConvergedOnState
        )
    }
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
    let mut residual_history: Vec<f64> = Vec::with_capacity(newton_config.max_iterations + 1);
    let mut prev_delta_v: Option<Vec<f64>> = None;
    let mut consecutive_anti_aligned = 0usize;
    let mut diverging = false;

    // b_norm is invariant over Newton iterations (rhs is constant), so
    // compute it once for the exhaustion-branch final_residual reporting.
    let b_norm_global: f64 = ws.rhs.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-14);

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
        residual_history.push(f_norm);

        // Divergence tracker: residual grew by more than 1.5× over the
        // last trend_window iterations. Latches true once tripped so
        // the exhaustion classifier can report Divergence.
        if residual_history.len() > newton_config.trend_window {
            let trend_idx = residual_history.len() - 1 - newton_config.trend_window;
            let r_old = residual_history[trend_idx];
            if r_old > 0.0 && f_norm > r_old * 1.5 {
                diverging = true;
            }
        }

        if f_norm < newton_config.tolerance * b_norm {
            unpack_velocity(&ws.v_packed, grid);
            let final_residual = f_norm / b_norm;
            debug!(
                outcome = ?NewtonOutcome::ConvergedOnResidual,
                iterations = k,
                final_residual,
                "newton solve completed"
            );
            return NewtonResult {
                outcome: NewtonOutcome::ConvergedOnResidual,
                iterations: k,
                final_residual,
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

        // Compute the effective Newton step actually applied (α·δv), used
        // below for the state-based convergence criterion and the
        // oscillation detector.
        let actual_step: Vec<f64> =
            ws.jfnk_delta_v.iter().map(|x| final_alpha * x).collect();

        // Oscillation detection: two consecutive Newton steps with a
        // strongly negative cosine signal back-and-forth motion. Skip the
        // check when either step is effectively zero (e.g. line search
        // gave up with α = 0) to avoid spurious NaN.
        if k >= newton_config.min_iterations_before_classification {
            if let Some(ref prev) = prev_delta_v {
                let dot: f64 =
                    actual_step.iter().zip(prev.iter()).map(|(a, b)| a * b).sum();
                let n_curr: f64 =
                    actual_step.iter().map(|x| x * x).sum::<f64>().sqrt();
                let n_prev: f64 = prev.iter().map(|x| x * x).sum::<f64>().sqrt();
                let denom = n_curr * n_prev;
                if denom > 1e-30 {
                    let cos_theta = dot / denom;
                    debug!(newton_iter = k, cos_theta, "newton step alignment");
                    if cos_theta < newton_config.oscillation_cosine_threshold {
                        consecutive_anti_aligned += 1;
                        if consecutive_anti_aligned >= 2 {
                            debug!(newton_iter = k, "oscillation detected, exiting");
                            unpack_velocity(&ws.v_packed, grid);
                            let final_residual = f_norm / b_norm;
                            debug!(
                                outcome = ?NewtonOutcome::Oscillation,
                                iterations = k + 1,
                                final_residual,
                                "newton solve completed"
                            );
                            return NewtonResult {
                                outcome: NewtonOutcome::Oscillation,
                                iterations: k + 1,
                                final_residual,
                                total_linear_iterations: total_linear,
                            };
                        }
                    } else {
                        consecutive_anti_aligned = 0;
                    }
                }
            }
        }
        prev_delta_v = Some(actual_step.clone());

        // State-based convergence has two independent acceptance paths:
        //   1. Physical state is frozen AND the residual trend is
        //      descending: Newton has found a true local minimum of |F|.
        //   2. Residual is near tolerance AND the recent history is flat:
        //      Newton cannot descend further through a non-smooth barrier
        //      in F(v) but the residual is effectively stable.
        if k >= newton_config.min_iterations_before_classification
            && residual_history.len() > newton_config.trend_window
        {
            let v_state_norm: f64 =
                v_old.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
            let step_norm: f64 =
                actual_step.iter().map(|x| x * x).sum::<f64>().sqrt();
            let relative_step = step_norm / v_state_norm;

            let trend_idx = residual_history.len() - 1 - newton_config.trend_window;
            let trend_descending = f_norm < residual_history[trend_idx];

            // Window covers (trend_window + 1) most recent residuals.
            let window = &residual_history[trend_idx..];
            let window_max = window.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let window_min = window.iter().cloned().fold(f64::INFINITY, f64::min);
            let window_spread = (window_max - window_min) / window_max.max(1e-30);

            debug!(
                newton_iter = k,
                relative_step,
                step_norm,
                v_state_norm,
                trend_descending,
                window_spread,
                f_norm_window_head = residual_history[trend_idx],
                "newton state criterion diagnostics"
            );

            let path_state_frozen =
                relative_step < newton_config.state_tolerance && trend_descending;
            let near_tolerance = f_norm
                < newton_config.tolerance
                    * b_norm
                    * newton_config.stagnation_residual_multiplier;
            let path_residual_stagnant = near_tolerance
                && window_spread < newton_config.stagnation_spread_threshold;

            if path_state_frozen || path_residual_stagnant {
                let reason = if path_state_frozen {
                    "state frozen with descending trend"
                } else {
                    "residual near tolerance with flat history"
                };
                debug!(
                    newton_iter = k,
                    relative_step,
                    f_norm,
                    window_spread,
                    reason,
                    "state-based convergence"
                );
                unpack_velocity(&ws.v_packed, grid);
                let final_residual = f_norm / b_norm;
                debug!(
                    outcome = ?NewtonOutcome::ConvergedOnState,
                    iterations = k + 1,
                    final_residual,
                    "newton solve completed"
                );
                return NewtonResult {
                    outcome: NewtonOutcome::ConvergedOnState,
                    iterations: k + 1,
                    final_residual,
                    total_linear_iterations: total_linear,
                };
            }
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

    // Classify the exhaustion outcome so downstream recovery logic can
    // react differently to Divergence vs Stagnation vs a generic
    // MaxIterations failure.
    let final_outcome = if diverging {
        NewtonOutcome::Divergence
    } else if residual_history.len() >= 3 {
        let recent: Vec<f64> = residual_history.iter().rev().take(3).copied().collect();
        let max_recent = recent.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_recent = recent.iter().cloned().fold(f64::INFINITY, f64::min);
        let spread = (max_recent - min_recent) / max_recent.max(1e-30);
        if spread < newton_config.stagnation_spread_threshold {
            NewtonOutcome::Stagnation
        } else {
            NewtonOutcome::MaxIterations
        }
    } else {
        NewtonOutcome::MaxIterations
    };

    let final_residual = residual_history
        .last()
        .copied()
        .map(|r| r / b_norm_global)
        .unwrap_or(f64::NAN);

    debug!(
        outcome = ?final_outcome,
        iterations = newton_config.max_iterations,
        final_residual,
        "newton solve completed"
    );

    NewtonResult {
        outcome: final_outcome,
        iterations: newton_config.max_iterations,
        final_residual,
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
        assert!(result.is_converged(), "Newton should converge for linear viscosity");
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
            result.is_converged(),
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
        assert!(newton_result.is_converged(), "Newton should converge");
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

        assert!(r_exact.is_converged(), "Exact Newton should converge");
        assert!(r_inexact.is_converged(), "Inexact Newton should converge");

        // Inexact should use fewer total linear iterations
        assert!(
            r_inexact.total_linear_iterations <= r_exact.total_linear_iterations,
            "Inexact should be cheaper or equal: {} vs {} linear iters",
            r_inexact.total_linear_iterations,
            r_exact.total_linear_iterations
        );
    }

    #[test]
    fn newton_state_convergence_on_stagnant_residual() {
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
        let newton_config = NewtonConfig {
            max_iterations: 25,
            tolerance: 1e-10,
            state_tolerance: 1e-3,
            trend_window: 3,
            ..NewtonConfig::default()
        };
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
        assert_eq!(
            result.outcome,
            NewtonOutcome::ConvergedOnState,
            "expected state-based convergence, got {:?} at iter {}",
            result.outcome,
            result.iterations
        );
        assert!(result.iterations < 25);
    }

    #[test]
    fn newton_accepts_state_when_residual_stagnates_near_tolerance() {
        // Configured so the residual-only criterion is unreachable and the
        // state_tolerance too tight for the moving-but-stable case: the
        // near-tolerance pathway is the only way out before max_iterations.
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
        let newton_config = NewtonConfig {
            max_iterations: 25,
            tolerance: 1e-6,
            state_tolerance: 1e-6,
            stagnation_residual_multiplier: 10.0,
            stagnation_spread_threshold: 0.15,
            trend_window: 3,
            min_iterations_before_classification: 3,
            ..NewtonConfig::default()
        };
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
            result.is_converged(),
            "expected convergence (residual or state), got {:?} at iter {}",
            result.outcome,
            result.iterations
        );
        assert!(result.iterations < 25, "should converge before exhausting iterations");
    }

    #[test]
    fn newton_outcome_is_converged_helper() {
        let r1 = NewtonResult {
            outcome: NewtonOutcome::ConvergedOnResidual,
            iterations: 5,
            final_residual: 0.001,
            total_linear_iterations: 50,
        };
        let r2 = NewtonResult {
            outcome: NewtonOutcome::ConvergedOnState,
            iterations: 8,
            final_residual: 0.06,
            total_linear_iterations: 80,
        };
        let r3 = NewtonResult {
            outcome: NewtonOutcome::Stagnation,
            iterations: 15,
            final_residual: 0.06,
            total_linear_iterations: 150,
        };
        let r4 = NewtonResult {
            outcome: NewtonOutcome::Oscillation,
            iterations: 7,
            final_residual: 0.08,
            total_linear_iterations: 70,
        };
        assert!(r1.is_converged());
        assert!(r2.is_converged());
        assert!(!r3.is_converged());
        assert!(!r4.is_converged());
    }
}
