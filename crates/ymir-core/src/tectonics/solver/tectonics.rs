//! Top-level orchestration of the thin viscous sheet simulation.

use super::advection::{compute_cfl_dt, compute_divergence_flux};
use super::config::{NonlinearSolver, TectonicsConfig};
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::newton::solve_velocity_newton;
use super::picard::solve_velocity_picard;
use super::traction::TractionField;
use super::workspace::{SolverWorkspace, StepStats};

/// Errors that can occur during a tectonic simulation run.
#[derive(Debug)]
pub enum SolverError {
    /// Nonlinear solver did not converge at the given timestep.
    NonlinearSolverDidNotConverge { step: usize },
    /// A NaN or Inf was detected in the solution.
    NumericalInstability { step: usize, field: &'static str },
    /// The simulation was cancelled via the progress callback.
    Cancelled,
}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolverError::NonlinearSolverDidNotConverge { step } => {
                write!(f, "Nonlinear solver did not converge at step {step}")
            }
            SolverError::NumericalInstability { step, field } => {
                write!(f, "NaN/Inf detected in {field} at step {step}")
            }
            SolverError::Cancelled => write!(f, "Simulation cancelled"),
        }
    }
}

impl std::error::Error for SolverError {}

/// Run a full tectonic simulation.
///
/// The `progress` callback is invoked after each timestep with
/// `(step, total, stats, s_field)`. The `s_field` is a reference to the
/// current crustal thickness field (safe to clone for snapshots).
/// Return `true` to continue, `false` to cancel.
pub fn run_tectonics<F>(
    config: &TectonicsConfig,
    plates: &TractionField,
    grid: &mut StaggeredGrid,
    workspace: &mut SolverWorkspace,
    mut progress: F,
) -> Result<(), SolverError>
where
    F: FnMut(usize, usize, &StepStats, &Field2D) -> bool,
{
    let n = grid.n;

    for step in 0..config.num_timesteps {
        // 1. Solve velocity — continuation only on first step (cold start)
        let need_continuation = config.continuation.enabled
            && config.picard.power_law_n > 1.0
            && step == 0;

        let t0 = std::time::Instant::now();
        let (converged, nl_iterations, linear_iterations) = if need_continuation {
            solve_with_continuation(grid, plates, config, workspace)
        } else {
            let result = solve_velocity_direct(grid, plates, config, workspace);
            // Fallback: if direct solve fails, try continuation as recovery
            if !result.0 && config.continuation.enabled && config.picard.power_law_n > 1.0 {
                solve_with_continuation(grid, plates, config, workspace)
            } else {
                result
            }
        };
        let solve_ms = t0.elapsed().as_millis();

        if !converged {
            return Err(SolverError::NonlinearSolverDidNotConverge { step });
        }

        // 2. CFL timestep (after velocity is known)
        let t1 = std::time::Instant::now();
        let dt_cfl = compute_cfl_dt(grid, config.cfl_factor);

        // 3. Adaptive timestep with retry on excessive clamping
        let s_backup: Vec<f64> = grid.s.data().to_vec();
        let mut dt = dt_cfl;
        let dt_min = dt * 0.01;
        let mut clamp_ratio = 0.0;

        for _retry in 0..5 {
            compute_divergence_flux(grid, &mut workspace.div_flux);

            for j in 0..n {
                for i in 0..n {
                    let s = grid.s.get(i, j) - dt * workspace.div_flux.get(i, j);
                    grid.s.set(i, j, s.clamp(config.s_min, config.s_max));
                }
            }

            let clamp_count = grid
                .s
                .data()
                .iter()
                .filter(|&&s| s <= config.s_min * 1.01 || s >= config.s_max * 0.99)
                .count();
            clamp_ratio = clamp_count as f64 / (n * n) as f64;

            if clamp_ratio < 0.05 || dt <= dt_min {
                break;
            }

            // Too many cells clamped — retry with smaller dt
            dt *= 0.5;
            grid.s.data_mut().copy_from_slice(&s_backup);
        }
        let advect_ms = t1.elapsed().as_millis();
        println!("Step {step}: solve={solve_ms}ms advect={advect_ms}ms nl={nl_iterations} lin={linear_iterations}");
        // 4. Update stats
        let mut max_v = 0.0_f64;
        let mut max_s = f64::NEG_INFINITY;
        let mut min_s = f64::INFINITY;
        for j in 0..n {
            for i in 0..n {
                let vx = grid.vx.get(i, j);
                let vy = grid.vy.get(i, j);
                max_v = max_v.max((vx * vx + vy * vy).sqrt());
                let s = grid.s.get(i, j);
                max_s = max_s.max(s);
                min_s = min_s.min(s);
            }
        }
        workspace.stats = StepStats {
            max_velocity: max_v,
            max_thickness: max_s,
            min_thickness: min_s,
            picard_iterations: nl_iterations,
            cg_iterations_last: linear_iterations,
            dt,
            clamp_ratio,
        };

        // 5. Callback — returns false to cancel
        if !progress(step, config.num_timesteps, &workspace.stats, &grid.s) {
            return Err(SolverError::Cancelled);
        }
    }

    Ok(())
}

/// Solve velocity directly (no continuation).
fn solve_velocity_direct(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    config: &TectonicsConfig,
    workspace: &mut SolverWorkspace,
) -> (bool, usize, usize) {
    match config.nonlinear_solver {
        NonlinearSolver::Picard => {
            let r = solve_velocity_picard(
                grid,
                plates,
                config.gravity_factor,
                &config.picard,
                workspace,
            );
            (r.converged, r.iterations, r.total_cg_iterations)
        }
        NonlinearSolver::Newton => {
            let r = solve_velocity_newton(
                grid,
                plates,
                config.gravity_factor,
                &config.picard,
                &config.newton,
                workspace,
            );
            (r.converged, r.iterations, r.total_linear_iterations)
        }
    }
}

/// Solve velocity using viscosity continuation: ramp n from 1 → target.
fn solve_with_continuation(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    config: &TectonicsConfig,
    workspace: &mut SolverWorkspace,
) -> (bool, usize, usize) {
    let target_eps = config.picard.strain_rate_min;
    let steps = &config.continuation.n_steps;
    let eps_start = config.continuation.eps_min_start.unwrap_or(target_eps);

    let mut total_nl = 0usize;
    let mut total_linear = 0usize;

    for (i, &n_exp) in steps.iter().enumerate() {
        // Interpolate ε_min from eps_start to target_eps
        let t = if steps.len() > 1 { i as f64 / (steps.len() - 1) as f64 } else { 1.0 };
        let eps_min = eps_start * (1.0 - t) + target_eps * t;

        // Adapt relaxation to nonlinearity level
        let relaxation = if n_exp <= 1.5 {
            0.9
        } else if n_exp <= 2.5 {
            0.6
        } else {
            0.4
        };

        let mut step_config = config.picard.clone();
        step_config.power_law_n = n_exp;
        step_config.strain_rate_min = eps_min;
        step_config.relaxation = relaxation;

        // Warm start: grid.vx/vy retain the solution from the previous step
        let (converged, iters, linear_iters) = match config.nonlinear_solver {
            NonlinearSolver::Picard => {
                let r = solve_velocity_picard(
                    grid,
                    plates,
                    config.gravity_factor,
                    &step_config,
                    workspace,
                );
                (r.converged, r.iterations, r.total_cg_iterations)
            }
            NonlinearSolver::Newton => {
                let r = solve_velocity_newton(
                    grid,
                    plates,
                    config.gravity_factor,
                    &step_config,
                    &config.newton,
                    workspace,
                );
                (r.converged, r.iterations, r.total_linear_iterations)
            }
        };
        total_nl += iters;
        total_linear += linear_iters;
        if !converged {
            return (false, total_nl, total_linear);
        }
    }

    (true, total_nl, total_linear)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::solver::config::{ContinuationConfig, NewtonConfig, PicardConfig};

    fn make_config(num_timesteps: usize) -> TectonicsConfig {
        TectonicsConfig {
            num_timesteps,
            gravity_factor: 1.0,
            cfl_factor: 0.1,
            s_min: 0.1,
            s_max: 2.5,
            nonlinear_solver: NonlinearSolver::Picard,
            picard: PicardConfig {
                max_iterations: 30,
                tolerance: 1e-3,
                relaxation: 1.0,
                cg_max_iter: 500,
                cg_tolerance: 1e-8,
                strain_rate_min: 1e-3,
                power_law_n: 1.0,
                eta_min: 1e-3,
                eta_max: 1e4,
            },
            newton: NewtonConfig::default(),
            continuation: ContinuationConfig { enabled: false, ..Default::default() },
        }
    }

    #[test]
    fn convergent_plates_thicken() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let plates = TractionField::two_plates_convergent(n, 1.0);
        let config = make_config(50);
        let mut ws = SolverWorkspace::new(n);

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        assert!(
            ws.stats.max_thickness > s_initial,
            "Convergent plates should thicken: max_s = {}",
            ws.stats.max_thickness
        );
    }

    #[test]
    fn divergent_plates_thin() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let plates = TractionField::two_plates_divergent(n, 1.0);
        let config = make_config(50);
        let mut ws = SolverWorkspace::new(n);

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        assert!(
            ws.stats.min_thickness < s_initial,
            "Divergent plates should thin: min_s = {}",
            ws.stats.min_thickness
        );
    }

    #[test]
    fn gpe_flattens_bump() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        // Background + broad Gaussian bump (σ²=0.02, well-resolved at 32²)
        let center = 0.5;
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                let r2 = (x - center).powi(2) + (y - center).powi(2);
                grid.s.set(i, j, 1.0 + 0.3 * (-r2 / 0.02).exp());
            }
        }

        let initial_var: f64 = {
            let mean = grid.s.data().iter().sum::<f64>() / (n * n) as f64;
            grid.s.data().iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n * n) as f64
        };

        let plates = TractionField::zero(n);
        let config = make_config(100);
        let mut ws = SolverWorkspace::new(n);

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        let final_var: f64 = {
            let mean = grid.s.data().iter().sum::<f64>() / (n * n) as f64;
            grid.s.data().iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n * n) as f64
        };

        assert!(
            final_var < initial_var,
            "GPE should flatten bump (reduce variance): initial_var = {initial_var}, final_var = {final_var}"
        );
    }

    #[test]
    fn no_nan_no_inf() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let plates = TractionField::two_plates_convergent(n, 0.5);
        let config = make_config(30);
        let mut ws = SolverWorkspace::new(n);

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok());

        for val in grid.s.data() {
            assert!(val.is_finite(), "S contains non-finite: {val}");
        }
        for val in grid.vx.data() {
            assert!(val.is_finite(), "vx contains non-finite: {val}");
        }
        for val in grid.vy.data() {
            assert!(val.is_finite(), "vy contains non-finite: {val}");
        }
    }

    #[test]
    fn continuation_enables_power_law_convergence() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let plates = TractionField::two_plates_convergent(n, 0.5);
        let config = TectonicsConfig {
            num_timesteps: 20,
            gravity_factor: 1.0,
            cfl_factor: 0.3,
            s_min: 0.1,
            s_max: 2.5,
            nonlinear_solver: NonlinearSolver::Picard,
            picard: PicardConfig {
                max_iterations: 60,
                power_law_n: 3.0,
                strain_rate_min: 1e-3,
                eta_min: 1e-3,
                eta_max: 1e4,
                relaxation: 0.5,
                ..PicardConfig::default()
            },
            newton: NewtonConfig::default(),
            continuation: ContinuationConfig::default(),
        };

        let mut ws = SolverWorkspace::new(n);
        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Continuation should enable convergence: {:?}", result.err());
        assert!(
            ws.stats.max_thickness > 1.0,
            "Convergent plates should thicken with power-law: max_s={}",
            ws.stats.max_thickness
        );
    }

    #[test]
    fn viscosity_clamp_works() {
        use crate::tectonics::solver::field::Field2D;
        use crate::tectonics::solver::picard::compute_viscosity;

        let n = 8;
        let mut strain = Field2D::new(n);
        let mut eta = Field2D::new(n);

        // Very low strain rate → very high viscosity without clamp
        strain.set(0, 0, 0.0);
        // Normal strain rate
        strain.set(1, 0, 0.1);

        compute_viscosity(&strain, 3.0, 1e-6, 1e-3, 1e3, &mut eta);

        assert!(eta.get(0, 0) <= 1e3 + 1e-10, "Should be clamped to eta_max: {}", eta.get(0, 0));
        assert!(eta.get(1, 0) >= 1e-3, "Normal cell below eta_min");
        assert!(eta.get(1, 0) <= 1e3, "Normal cell above eta_max");
    }
}
