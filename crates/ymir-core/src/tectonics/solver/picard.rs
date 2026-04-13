//! Picard (fixed-point) iteration for nonlinear Stokes with power-law viscosity.

use rayon::prelude::*;
use tracing::debug;

use super::config::PicardConfig;
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::linear_solve::solve_cg;
use super::stokes::{apply_stokes, compute_jacobi_precond, compute_rhs};
use super::traction::TractionField;
use super::workspace::SolverWorkspace;

const PAR_THRESHOLD: usize = 64;

/// Result of a Picard solve.
pub struct PicardResult {
    pub converged: bool,
    pub iterations: usize,
    pub final_residual: f64,
    pub total_cg_iterations: usize,
}

/// Compute the second invariant of strain rate at cell centers.
#[allow(clippy::needless_range_loop)]
pub fn compute_strain_rate(grid: &StaggeredGrid, out: &mut Field2D) {
    let n = grid.n;
    let idx = &grid.idx;
    let inv_dx = 1.0 / grid.dx;

    let process_row = |j: usize, row: &mut [f64]| {
        let nj = idx.next(j);
        let pj = idx.prev(j);
        for i in 0..n {
            let ni = idx.next(i);
            let pi = idx.prev(i);

            let dvx_dx = (grid.vx.get(ni, j) - grid.vx.get(i, j)) * inv_dx;
            let dvy_dy = (grid.vy.get(i, nj) - grid.vy.get(i, j)) * inv_dx;

            let c00 = (grid.vx.get(i, j) - grid.vx.get(i, pj)) * inv_dx
                + (grid.vy.get(i, j) - grid.vy.get(pi, j)) * inv_dx;
            let c10 = (grid.vx.get(ni, j) - grid.vx.get(ni, pj)) * inv_dx
                + (grid.vy.get(ni, j) - grid.vy.get(i, j)) * inv_dx;
            let c01 = (grid.vx.get(i, nj) - grid.vx.get(i, j)) * inv_dx
                + (grid.vy.get(i, nj) - grid.vy.get(pi, nj)) * inv_dx;
            let c11 = (grid.vx.get(ni, nj) - grid.vx.get(ni, j)) * inv_dx
                + (grid.vy.get(ni, nj) - grid.vy.get(i, nj)) * inv_dx;

            let dvx_dy_plus_dvy_dx = 0.25 * (c00 + c10 + c01 + c11);

            row[i] = (0.5 * dvx_dx * dvx_dx
                + 0.5 * dvy_dy * dvy_dy
                + 0.25 * dvx_dy_plus_dvy_dx * dvx_dy_plus_dvy_dx)
                .sqrt();
        }
    };

    if n >= PAR_THRESHOLD {
        out.data_mut().par_chunks_mut(n).enumerate().for_each(|(j, row)| process_row(j, row));
    } else {
        for j in 0..n {
            let s = j * n;
            process_row(j, &mut out.data_mut()[s..s + n]);
        }
    }
}

/// Compute viscosity from strain rate: η = clamp((ε̇_II + ε_min)^(1/n - 1), η_min, η_max).
pub fn compute_viscosity(
    strain_rate: &Field2D,
    n_exp: f64,
    eps_min: f64,
    eta_min: f64,
    eta_max: f64,
    eta: &mut Field2D,
) {
    let exponent = 1.0 / n_exp - 1.0;
    let n = strain_rate.n();

    if n >= PAR_THRESHOLD {
        strain_rate.data().par_iter().zip(eta.data_mut().par_iter_mut()).for_each(
            |(&sr, eta_val)| {
                *eta_val = (sr + eps_min).powf(exponent).clamp(eta_min, eta_max);
            },
        );
    } else {
        for k in 0..n * n {
            let sr = strain_rate.data()[k];
            eta.data_mut()[k] = (sr + eps_min).powf(exponent).clamp(eta_min, eta_max);
        }
    }
}

/// Solve for velocity using Picard iteration.
pub fn solve_velocity_picard(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    gravity_factor: f64,
    rho_continental: f64,
    rho_mantle: f64,
    config: &PicardConfig,
    ws: &mut SolverWorkspace,
) -> PicardResult {
    let n = grid.n;
    let nn2 = 2 * n * n;
    let mut total_cg = 0usize;

    // Pack current velocity as initial guess
    super::workspace::pack_velocity(grid, &mut ws.v_packed);

    // Compute RHS (does not change during Picard iteration)
    compute_rhs(grid, plates, gravity_factor, rho_continental, rho_mantle, &mut ws.rhs);

    // Project out null space (constant mode) — periodic Stokes has a rank-2 null space
    let n2 = n * n;
    let mean_vx: f64 = ws.rhs[..n2].iter().sum::<f64>() / n2 as f64;
    let mean_vy: f64 = ws.rhs[n2..nn2].iter().sum::<f64>() / n2 as f64;
    for val in &mut ws.rhs[..n2] {
        *val -= mean_vx;
    }
    for val in &mut ws.rhs[n2..nn2] {
        *val -= mean_vy;
    }

    for k in 0..config.max_iterations {
        // Save previous velocity
        ws.v_prev.copy_from_slice(&ws.v_packed);

        // Unpack into grid for strain rate computation
        super::workspace::unpack_velocity(&ws.v_packed, grid);

        // Compute strain rate and viscosity
        compute_strain_rate(grid, &mut ws.strain_rate);
        compute_viscosity(
            &ws.strain_rate,
            config.power_law_n,
            config.strain_rate_min,
            config.eta_min,
            config.eta_max,
            &mut ws.eta,
        );

        // Update Jacobi preconditioner
        compute_jacobi_precond(&ws.eta, grid, &mut ws.jacobi_precond);

        // Solve A(η) v = b
        let eta_ref = &ws.eta;
        let grid_ref = &*grid;
        let precond_ref = &ws.jacobi_precond;
        let cg_result = solve_cg(
            &mut ws.v_packed,
            &ws.rhs,
            |v_in, v_out| apply_stokes(v_in, eta_ref, grid_ref, v_out),
            |r, z| super::linear_solve::apply_jacobi(precond_ref, r, z),
            &mut ws.cg,
            config.cg_max_iter,
            config.cg_tolerance,
        );
        total_cg += cg_result.iterations;

        // Project out null space from solution
        let mean_vx: f64 = ws.v_packed[..n2].iter().sum::<f64>() / n2 as f64;
        let mean_vy: f64 = ws.v_packed[n2..nn2].iter().sum::<f64>() / n2 as f64;
        for val in &mut ws.v_packed[..n2] {
            *val -= mean_vx;
        }
        for val in &mut ws.v_packed[n2..nn2] {
            *val -= mean_vy;
        }

        // Under-relaxation: v = α·v_new + (1-α)·v_prev
        let alpha = config.relaxation;
        for i in 0..nn2 {
            ws.v_packed[i] = alpha * ws.v_packed[i] + (1.0 - alpha) * ws.v_prev[i];
        }

        // Convergence check
        let mut diff_sq = 0.0;
        let mut v_sq = 0.0;
        for i in 0..nn2 {
            let d = ws.v_packed[i] - ws.v_prev[i];
            diff_sq += d * d;
            v_sq += ws.v_packed[i] * ws.v_packed[i];
        }
        let rel_change = diff_sq.sqrt() / v_sq.sqrt().max(1e-30);

        debug!(
            picard_iter = k,
            rel_change,
            cg_iters = cg_result.iterations,
            cg_converged = cg_result.converged,
            "picard iteration"
        );

        if rel_change < config.tolerance {
            super::workspace::unpack_velocity(&ws.v_packed, grid);
            return PicardResult {
                converged: true,
                iterations: k + 1,
                final_residual: rel_change,
                total_cg_iterations: total_cg,
            };
        }
    }

    super::workspace::unpack_velocity(&ws.v_packed, grid);
    PicardResult {
        converged: false,
        iterations: config.max_iterations,
        final_residual: f64::NAN,
        total_cg_iterations: total_cg,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_viscosity_one_iteration() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        // Set up a non-trivial S field
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                grid.s.set(i, j, 1.0 + 0.3 * (2.0 * std::f64::consts::PI * x).sin());
            }
        }

        let plates = TractionField::uniform(n, 0.1, 0.0);
        let config = PicardConfig {
            max_iterations: 50,
            tolerance: 1e-10,
            relaxation: 1.0,
            cg_max_iter: 500,
            cg_tolerance: 1e-10,
            strain_rate_min: 1e-6,
            power_law_n: 1.0,
            ..PicardConfig::default()
        };
        let mut ws = SolverWorkspace::new(n);

        let result = solve_velocity_picard(&mut grid, &plates, 1.0, 0.0, 0.0, &config, &mut ws);
        assert!(result.converged, "Should converge for linear viscosity");
        assert!(
            result.iterations <= 2,
            "Linear viscosity (n=1) should converge in ≤ 2 Picard iterations, got {}",
            result.iterations
        );
    }

    #[test]
    fn power_law_converges() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);

        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let plates = TractionField::two_plates_convergent(n, 0.5);
        let config = PicardConfig {
            max_iterations: 50,
            tolerance: 1e-4,
            relaxation: 0.7,
            cg_max_iter: 500,
            cg_tolerance: 1e-8,
            strain_rate_min: 1e-3,
            power_law_n: 3.0,
            ..PicardConfig::default()
        };
        let mut ws = SolverWorkspace::new(n);

        let result = solve_velocity_picard(&mut grid, &plates, 1.0, 0.0, 0.0, &config, &mut ws);
        assert!(
            result.converged,
            "Power-law Picard should converge, got {} iterations",
            result.iterations
        );
        assert!(
            result.iterations <= 50,
            "Should converge in ≤ 50 iterations, got {}",
            result.iterations
        );
    }
}
