//! Top-level orchestration of the thin viscous sheet simulation.

use super::advection::{compute_cfl_dt, compute_divergence_flux};
use super::config::{NonlinearSolver, TectonicsConfig};
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
        }
    }
}

impl std::error::Error for SolverError {}

/// Run a full tectonic simulation.
///
/// The `progress` callback is invoked after each timestep with (step, total, stats).
pub fn run_tectonics<F>(
    config: &TectonicsConfig,
    plates: &TractionField,
    grid: &mut StaggeredGrid,
    workspace: &mut SolverWorkspace,
    progress: F,
) -> Result<(), SolverError>
where
    F: Fn(usize, usize, &StepStats),
{
    let n = grid.n;

    for step in 0..config.num_timesteps {
        // 1. Solve velocity via selected nonlinear solver
        let (converged, nl_iterations, linear_iterations) = match config.nonlinear_solver {
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
        };

        if !converged {
            return Err(SolverError::NonlinearSolverDidNotConverge { step });
        }

        // 2. CFL timestep (after velocity is known)
        let dt = compute_cfl_dt(grid, config.cfl_factor);

        // 3. Compute divergence of flux
        compute_divergence_flux(grid, &mut workspace.div_flux);

        // 4. Euler explicit update: S_new = S - dt * div (no source terms yet)
        for j in 0..n {
            for i in 0..n {
                let s = grid.s.get(i, j) - dt * workspace.div_flux.get(i, j);
                grid.s.set(i, j, s);
            }
        }

        // 5. Clip S
        for j in 0..n {
            for i in 0..n {
                let s = grid.s.get(i, j).clamp(config.s_min, config.s_max);
                grid.s.set(i, j, s);
            }
        }

        // 6. Update stats
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
        };

        // 7. Callback
        progress(step, config.num_timesteps, &workspace.stats);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::solver::config::{NewtonConfig, PicardConfig};

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
                strain_rate_min: 1e-6,
                power_law_n: 1.0,
            },
            newton: NewtonConfig::default(),
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

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _| {});
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

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _| {});
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

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _| {});
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

        let result = run_tectonics(&config, &plates, &mut grid, &mut ws, |_, _, _| {});
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
}
