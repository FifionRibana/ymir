//! Thin viscous sheet momentum solver (England & McKenzie 1982).
//!
//! Solves the depth-integrated horizontal momentum balance
//! ```text
//!   -∇·(2 η ε̇(v)) = f
//! ```
//! for the 2-D velocity `v = (v_x, v_y)` on a fully periodic MAC
//! grid. **No incompressibility constraint; no pressure unknown.**
//! `∇·v ≠ 0` is physically meaningful (thickening rate); pressure
//! does not enter as a Lagrange multiplier in this formulation.
//!
//! The operator `A v = -∇·(2 η ε̇(v))` is symmetric positive
//! definite on the zero-mean velocity subspace (modulo the 2-D
//! rigid-body translation null space). A single preconditioned
//! conjugate-gradient solve per time step suffices — no nested
//! iteration, no saddle point, no Schur complement.
//!
//! # Gauge fixing
//!
//! The periodic null space is handled by subtracting the mean of
//! each velocity component both inside every preconditioner
//! application (before and after `M⁻¹`) and once more after the
//! solve completes.
//!
//! # Solver trait
//!
//! CG is used behind a generic [`LinearSolver`][solver::LinearSolver]
//! trait so that Step 3's BiCGSTAB (needed once yielding makes the
//! system non-symmetric) can swap in as a drop-in replacement
//! without reshaping the caller side.

pub mod amg;
pub mod continuation;
pub mod nonlinear_solver;
pub mod nullspace;
pub mod operator;
pub mod picard;
pub mod precond;
pub mod snapshot;
pub mod solver;
pub mod sparse_assembly;

use super::field::Field2D;
use operator::{apply_momentum, momentum_diagonal, StokesGrid};
use precond::VelocityJacobi;
use solver::{ConjugateGradient, LinearSolver, SolverStats};

pub use operator::StokesGrid as Grid;

/// Configuration for a sheet solve.
#[derive(Clone, Copy, Debug)]
pub struct SheetConfig {
    pub tol: f64,
    pub max_iter: usize,
    /// Minimum absolute value used when inverting a diagonal entry in
    /// the Jacobi preconditioner. Protects against zero diagonals
    /// from degenerate η fields without annihilating near-singular
    /// information.
    pub diag_floor: f64,
}

impl Default for SheetConfig {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            max_iter: 1000,
            diag_floor: 1e-20,
        }
    }
}

/// Aggregate statistics from a sheet solve.
#[derive(Clone, Copy, Debug, Default)]
pub struct SheetStats {
    pub iterations: usize,
    pub final_residual: f64,
    pub initial_residual: f64,
    pub mean_vx_after: f64,
    pub mean_vy_after: f64,
    pub converged: bool,
}

/// Solve the thin-sheet momentum balance on a MAC grid.
///
/// Inputs `fx`, `fy` are body-force components at velocity faces.
/// Outputs `vx`, `vy` are overwritten with the solution; they also
/// serve as the initial guess.
///
/// `drag_diag` carries the optional Step-4 basal-drag contribution
/// `Br · S̃^exp` as a cell-centered field. `None` disables drag
/// entirely (zero-cost); `Some(&field)` is augmented onto both the
/// matvec and the preconditioner diagonal, with face averaging done
/// inside [`apply_momentum`] and [`momentum_diagonal`].
pub fn solve_sheet(
    grid: &StokesGrid,
    eta: &Field2D,
    drag_diag: Option<&Field2D>,
    fx: &[f64],
    fy: &[f64],
    vx: &mut [f64],
    vy: &mut [f64],
    cfg: &SheetConfig,
) -> SheetStats {
    let n = grid.n_cells();
    assert_eq!(fx.len(), n);
    assert_eq!(fy.len(), n);
    assert_eq!(vx.len(), n);
    assert_eq!(vy.len(), n);

    // --- Preconditioner (η and drag_diag are frozen during the solve) ---
    let mut diag_vx = vec![0.0; n];
    let mut diag_vy = vec![0.0; n];
    momentum_diagonal(grid, eta, drag_diag, &mut diag_vx, &mut diag_vy);
    let vjac = VelocityJacobi::from_diagonal(&diag_vx, &diag_vy, cfg.diag_floor);

    // --- Pack RHS and initial guess into [vx; vy] layout ---
    let mut b_pack = Vec::with_capacity(2 * n);
    b_pack.extend_from_slice(fx);
    b_pack.extend_from_slice(fy);
    // Gauge-fix the RHS: any null-space component is inconsistent with
    // the SPD reduced system.
    {
        let (bx, by) = b_pack.split_at_mut(n);
        nullspace::project_velocity(bx, by);
    }

    let mut x_pack = vec![0.0; 2 * n];
    x_pack[..n].copy_from_slice(vx);
    x_pack[n..].copy_from_slice(vy);

    // --- Matrix-vector and preconditioner closures ---
    let mut tmp_ax = vec![0.0; n];
    let mut tmp_ay = vec![0.0; n];
    let mut matvec = |v: &[f64], out: &mut [f64]| {
        let (vx_in, vy_in) = v.split_at(n);
        let (out_x, out_y) = out.split_at_mut(n);
        apply_momentum(grid, eta, drag_diag, vx_in, vy_in, &mut tmp_ax, &mut tmp_ay);
        out_x.copy_from_slice(&tmp_ax);
        out_y.copy_from_slice(&tmp_ay);
    };
    let mut precond = |r: &[f64], z: &mut [f64]| {
        vjac.apply(r, z);
    };

    // --- Solve ---
    let cg = ConjugateGradient::new(cfg.tol, cfg.max_iter);
    let cg_stats: SolverStats = cg.solve(&mut matvec, &mut precond, &b_pack, &mut x_pack);

    vx.copy_from_slice(&x_pack[..n]);
    vy.copy_from_slice(&x_pack[n..]);

    // --- Final null-space clean-up ---
    nullspace::project_velocity(vx, vy);

    SheetStats {
        iterations: cg_stats.iterations,
        final_residual: cg_stats.final_residual,
        initial_residual: cg_stats.initial_residual,
        mean_vx_after: nullspace::mean(vx),
        mean_vy_after: nullspace::mean(vy),
        converged: cg_stats.converged(),
    }
}
