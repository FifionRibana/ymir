//! Stokes subsystem: operator, null-space projection, preconditioner,
//! generic linear solver, and the pressure Schur-complement
//! coordinator.
//!
//! # Gauge fixing
//!
//! The fully periodic Stokes system has a 1-D pressure null space and
//! a 2-D velocity null space. Both are projected out **inside every
//! preconditioner application** (before and after `M⁻¹`), and also
//! once more after the full solve completes. See
//! [`nullspace`] and [`precond`].
//!
//! # Algorithm
//!
//! Pressure Schur complement with nested CG:
//!
//! ```text
//!   | A    -B^T | |v|   |f|
//!   | B      0  | |p| = |0|
//! ```
//!
//! Eliminating `v = A⁻¹(f + B^T p)` gives `S p = -B A⁻¹ f` with
//! `S = B A⁻¹ B^T` symmetric positive definite on the zero-mean
//! subspace. The outer solve runs CG on S; each `S·p` application
//! uses an inner CG for `A⁻¹·(B^T p)`.
//!
//! Inner tolerance is kept tighter than the outer tolerance so the
//! inexact-matvec perturbation stays below the outer convergence
//! criterion. A flexible outer iteration can replace this later if
//! inner costs dominate.

pub mod nullspace;
pub mod operator;
pub mod precond;
pub mod solver;

use std::cell::RefCell;

use super::field::Field2D;
use operator::{
    apply_divergence, apply_divergence_transpose, apply_momentum, momentum_diagonal, StokesGrid,
};
use precond::{PressureMass, VelocityJacobi};
use solver::{ConjugateGradient, LinearSolver, SolverStats};

pub use operator::StokesGrid as Grid;

/// Configuration for a Stokes solve.
#[derive(Clone, Copy, Debug)]
pub struct StokesConfig {
    pub outer_tol: f64,
    pub inner_tol: f64,
    pub outer_max_iter: usize,
    pub inner_max_iter: usize,
    /// Minimum absolute value used when inverting a diagonal entry in
    /// the Jacobi preconditioner.
    pub diag_floor: f64,
}

impl Default for StokesConfig {
    fn default() -> Self {
        Self {
            outer_tol: 1e-8,
            inner_tol: 1e-10,
            outer_max_iter: 200,
            inner_max_iter: 500,
            diag_floor: 1e-20,
        }
    }
}

/// Aggregate statistics from a Stokes solve.
#[derive(Clone, Copy, Debug, Default)]
pub struct StokesStats {
    pub outer_iterations: usize,
    pub outer_residual: f64,
    pub inner_iterations_total: usize,
    pub inner_solves: usize,
    pub inner_iterations_max: usize,
    pub mean_p_after: f64,
    pub mean_vx_after: f64,
    pub mean_vy_after: f64,
    pub converged: bool,
}

/// Counter accumulated across all nested inner solves in a single
/// outer Stokes solve.
#[derive(Clone, Copy, Debug, Default)]
struct InnerAcc {
    solves: usize,
    iter_total: usize,
    iter_max: usize,
}

/// Solve `A·(wx, wy) = (rhs_x, rhs_y)` for the packed velocity
/// vector, via preconditioned CG on the momentum block. The solution
/// is written in-place into `(wx, wy)` (which also serves as the
/// initial guess).
fn inner_velocity_solve(
    grid: &StokesGrid,
    eta: &Field2D,
    vjac: &VelocityJacobi,
    cg: &ConjugateGradient,
    rhs_x: &[f64],
    rhs_y: &[f64],
    wx: &mut [f64],
    wy: &mut [f64],
) -> SolverStats {
    let n = grid.n_cells();
    let mut b_pack = Vec::with_capacity(2 * n);
    b_pack.extend_from_slice(rhs_x);
    b_pack.extend_from_slice(rhs_y);
    let mut x_pack = vec![0.0; 2 * n];
    x_pack[..n].copy_from_slice(wx);
    x_pack[n..].copy_from_slice(wy);

    let mut tmp_ax = vec![0.0; n];
    let mut tmp_ay = vec![0.0; n];
    let mut matvec = |v: &[f64], out: &mut [f64]| {
        let (vx_in, vy_in) = v.split_at(n);
        let (out_x, out_y) = out.split_at_mut(n);
        apply_momentum(grid, eta, vx_in, vy_in, &mut tmp_ax, &mut tmp_ay);
        out_x.copy_from_slice(&tmp_ax);
        out_y.copy_from_slice(&tmp_ay);
    };
    let mut precond = |r: &[f64], z: &mut [f64]| {
        vjac.apply(r, z);
    };
    let stats = cg.solve(&mut matvec, &mut precond, &b_pack, &mut x_pack);
    wx.copy_from_slice(&x_pack[..n]);
    wy.copy_from_slice(&x_pack[n..]);
    stats
}

/// Apply the Schur complement `S p = B · A⁻¹ · B^T p` to a pressure
/// vector. Updates `inner_acc` with the inner-solve accounting.
fn schur_apply(
    grid: &StokesGrid,
    eta: &Field2D,
    vjac: &VelocityJacobi,
    cg: &ConjugateGradient,
    p_in: &[f64],
    out: &mut [f64],
    inner_acc: &RefCell<InnerAcc>,
) {
    let n = grid.n_cells();
    let mut btq_x = vec![0.0; n];
    let mut btq_y = vec![0.0; n];
    apply_divergence_transpose(grid, p_in, &mut btq_x, &mut btq_y);
    let mut w_vx = vec![0.0; n];
    let mut w_vy = vec![0.0; n];
    let s = inner_velocity_solve(grid, eta, vjac, cg, &btq_x, &btq_y, &mut w_vx, &mut w_vy);
    {
        let mut acc = inner_acc.borrow_mut();
        acc.solves += 1;
        acc.iter_total += s.iterations;
        acc.iter_max = acc.iter_max.max(s.iterations);
    }
    apply_divergence(grid, &w_vx, &w_vy, out);
    nullspace::subtract_mean(out);
}

/// Solve the periodic Stokes system on a MAC grid.
///
/// Inputs `fx`, `fy` are body-force components at velocity faces.
/// Outputs `vx`, `vy`, `p` are overwritten with the solution.
pub fn solve_stokes(
    grid: &StokesGrid,
    eta: &Field2D,
    fx: &[f64],
    fy: &[f64],
    vx: &mut [f64],
    vy: &mut [f64],
    p: &mut [f64],
    cfg: &StokesConfig,
) -> StokesStats {
    let n = grid.n_cells();
    assert_eq!(fx.len(), n);
    assert_eq!(fy.len(), n);
    assert_eq!(vx.len(), n);
    assert_eq!(vy.len(), n);
    assert_eq!(p.len(), n);

    // --- Preconditioners (constant during the solve; η is frozen) ---
    let mut diag_vx = vec![0.0; n];
    let mut diag_vy = vec![0.0; n];
    momentum_diagonal(grid, eta, &mut diag_vx, &mut diag_vy);
    let vjac = VelocityJacobi::from_diagonal(&diag_vx, &diag_vy, cfg.diag_floor);
    let pmass = PressureMass::from_eta(eta, cfg.diag_floor);

    let inner_cg = ConjugateGradient::new(cfg.inner_tol, cfg.inner_max_iter);
    let outer_cg = ConjugateGradient::new(cfg.outer_tol, cfg.outer_max_iter);

    let mut stats = StokesStats::default();
    let inner_acc = RefCell::new(InnerAcc::default());

    // --- Step 1: aux = A⁻¹ f ------------------------------------------
    let mut aux_vx = vec![0.0; n];
    let mut aux_vy = vec![0.0; n];
    let s1 = inner_velocity_solve(grid, eta, &vjac, &inner_cg, fx, fy, &mut aux_vx, &mut aux_vy);
    {
        let mut acc = inner_acc.borrow_mut();
        acc.solves += 1;
        acc.iter_total += s1.iterations;
        acc.iter_max = acc.iter_max.max(s1.iterations);
    }

    // --- Step 2: g = -B·aux --------------------------------------------
    let mut g = vec![0.0; n];
    apply_divergence(grid, &aux_vx, &aux_vy, &mut g);
    for gk in g.iter_mut() {
        *gk = -*gk;
    }
    nullspace::subtract_mean(&mut g);

    // --- Step 3: outer CG on S p = g -----------------------------------
    for pk in p.iter_mut() {
        *pk = 0.0;
    }

    let mut schur_mv = |q: &[f64], out: &mut [f64]| {
        schur_apply(grid, eta, &vjac, &inner_cg, q, out, &inner_acc);
    };
    let mut pmass_pc = |r: &[f64], z: &mut [f64]| {
        pmass.apply(r, z);
    };
    let outer_stats = outer_cg.solve(&mut schur_mv, &mut pmass_pc, &g, p);

    stats.outer_iterations = outer_stats.iterations;
    stats.outer_residual = outer_stats.final_residual;
    stats.converged = outer_stats.converged();

    // --- Step 4: v = A⁻¹ (f + B^T p) ----------------------------------
    let mut bt_x = vec![0.0; n];
    let mut bt_y = vec![0.0; n];
    let mut rhs_x = vec![0.0; n];
    let mut rhs_y = vec![0.0; n];
    apply_divergence_transpose(grid, p, &mut bt_x, &mut bt_y);
    for k in 0..n {
        rhs_x[k] = fx[k] + bt_x[k];
        rhs_y[k] = fy[k] + bt_y[k];
    }
    for vk in vx.iter_mut().chain(vy.iter_mut()) {
        *vk = 0.0;
    }
    let s4 = inner_velocity_solve(grid, eta, &vjac, &inner_cg, &rhs_x, &rhs_y, vx, vy);
    {
        let mut acc = inner_acc.borrow_mut();
        acc.solves += 1;
        acc.iter_total += s4.iterations;
        acc.iter_max = acc.iter_max.max(s4.iterations);
    }

    // --- Step 5: clean null-space of final iterates -------------------
    nullspace::project_pressure(p);
    nullspace::project_velocity(vx, vy);

    {
        let acc = inner_acc.borrow();
        stats.inner_solves = acc.solves;
        stats.inner_iterations_total = acc.iter_total;
        stats.inner_iterations_max = acc.iter_max;
    }
    stats.mean_p_after = nullspace::mean(p);
    stats.mean_vx_after = nullspace::mean(vx);
    stats.mean_vy_after = nullspace::mean(vy);

    stats
}
