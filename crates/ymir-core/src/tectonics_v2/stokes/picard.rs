//! Picard (fixed-point) nonlinear solver — the dual-track fallback to
//! Newton for transients where Newton diverges.
//!
//! At each outer step the Picard solver freezes the viscosity at the
//! current iterate's strain rate (`η_k = η(ε̇_II(v_k))`) and solves
//! the *linear* thin-sheet problem `-∇·(2 η_k ε̇(v_{k+1})) = rhs`. The
//! Jacobian used to derive the update is the Picard block of the
//! Newton Jacobian — the "Newton-extra" term is omitted. This makes
//! each outer step globally convergent (any SPD preconditioner on the
//! Picard block solves it) at the cost of loss of quadratic
//! convergence.
//!
//! The outer loop shares the same `NonlinearOutcome` reporting as
//! [`NewtonSolver`] so the test `picard_parity` can substitute one
//! for the other transparently.

use super::super::cratonic::CratonicState;
use super::super::field::Field2D;
use super::super::rheology::{self, StrainRate, ViscosityLaw};
use super::nonlinear_solver::{NonlinearOutcome, NonlinearSolver, NonlinearTrace};
use super::nullspace;
use super::operator::{apply_momentum, momentum_diagonal, StokesGrid};
use super::precond::VelocityJacobi;
use super::solver::{ConjugateGradient, LinearSolver, SolverStats};

#[derive(Clone, Copy, Debug)]
pub struct PicardConfig {
    pub rel_tol: f64,
    pub abs_tol: f64,
    pub max_outer_iters: u32,
    pub linear_tol: f64,
    pub linear_max_iter: usize,
    pub diag_floor: f64,
    /// Relaxation ω ∈ (0, 1]. `1.0` = pure Picard;
    /// lower values damp oscillations at the cost of speed.
    pub relaxation: f64,
}

impl Default for PicardConfig {
    fn default() -> Self {
        Self {
            rel_tol: 1.0e-6,
            abs_tol: 1.0e-10,
            max_outer_iters: 100,
            linear_tol: 1.0e-8,
            linear_max_iter: 2000,
            diag_floor: 1.0e-20,
            relaxation: 0.8,
        }
    }
}

pub struct PicardSolver {
    pub cfg: PicardConfig,
}

impl PicardSolver {
    pub fn new(cfg: PicardConfig) -> Self {
        Self { cfg }
    }
}

impl Default for PicardSolver {
    fn default() -> Self {
        Self { cfg: PicardConfig::default() }
    }
}

fn vec_norm(a: &[f64], b: &[f64]) -> f64 {
    let s = a.iter().map(|x| x * x).sum::<f64>() + b.iter().map(|x| x * x).sum::<f64>();
    s.sqrt()
}

fn compute_residual(
    grid: &StokesGrid,
    law: &ViscosityLaw,
    drag_diag: Option<&Field2D>,
    cratonic: Option<&CratonicState>,
    vx: &[f64],
    vy: &[f64],
    rhs_x: &[f64],
    rhs_y: &[f64],
    out_x: &mut [f64],
    out_y: &mut [f64],
) {
    let sr = StrainRate::compute(
        grid.nx, grid.ny, grid.dx, grid.dy,
        &grid.idx_x, &grid.idx_y, vx, vy,
    );
    let eta = rheology::build_eta_field(law, &sr.eps_ii_center, cratonic);
    apply_momentum(grid, &eta, drag_diag, vx, vy, out_x, out_y);
    for k in 0..out_x.len() {
        out_x[k] -= rhs_x[k];
        out_y[k] -= rhs_y[k];
    }
    nullspace::project_velocity(out_x, out_y);
}

impl NonlinearSolver for PicardSolver {
    fn solve(
        &self,
        grid: &StokesGrid,
        law: &ViscosityLaw,
        drag_diag: Option<&Field2D>,
        cratonic: Option<&CratonicState>,
        rhs_x: &[f64],
        rhs_y: &[f64],
        vx: &mut [f64],
        vy: &mut [f64],
        _linear_solver: &dyn LinearSolver,
    ) -> NonlinearOutcome {
        let n = grid.n_cells();
        let mut trace = NonlinearTrace::default();
        let mut r_x = vec![0.0; n];
        let mut r_y = vec![0.0; n];
        compute_residual(grid, law, drag_diag, cratonic, vx, vy, rhs_x, rhs_y, &mut r_x, &mut r_y);
        let r0_norm = vec_norm(&r_x, &r_y);
        trace.residuals.push(r0_norm);
        let abs_tol_eff = self.cfg.abs_tol.max(10.0 * self.cfg.linear_tol);
        if r0_norm <= abs_tol_eff {
            return NonlinearOutcome::Converged {
                outer_iters: 0,
                final_residual: r0_norm,
                linear_iters_total: 0,
                trace,
            };
        }

        let mut prev_resid = r0_norm;
        let mut linear_iters_total = 0u32;

        for k in 0..self.cfg.max_outer_iters {
            if prev_resid > 10.0 * r0_norm && r0_norm > 0.0 {
                return NonlinearOutcome::Diverged {
                    outer_iters: k,
                    last_residual: prev_resid,
                    trace,
                };
            }
            if trace.residuals.len() >= 4 {
                let m = trace.residuals.len();
                let deltas = [
                    (trace.residuals[m - 3] - trace.residuals[m - 2]) / trace.residuals[m - 3].max(1e-300),
                    (trace.residuals[m - 2] - trace.residuals[m - 1]) / trace.residuals[m - 2].max(1e-300),
                    (trace.residuals[m - 4] - trace.residuals[m - 3]) / trace.residuals[m - 4].max(1e-300),
                ];
                let conv_threshold = self.cfg.abs_tol.max(self.cfg.rel_tol * r0_norm);
                let far_from_convergence = prev_resid > 10.0 * conv_threshold;
                if far_from_convergence && deltas.iter().all(|d| d.abs() < 0.01) {
                    return NonlinearOutcome::Stalled { outer_iters: k, trace };
                }
            }

            // Freeze η at the current iterate.
            let sr = StrainRate::compute(
                grid.nx, grid.ny, grid.dx, grid.dy,
                &grid.idx_x, &grid.idx_y, vx, vy,
            );
            let eta = rheology::build_eta_field(law, &sr.eps_ii_center, cratonic);

            // Solve the Picard problem -∇·(2 η ε̇(v_{k+1})) = rhs.
            let mut diag_vx = vec![0.0; n];
            let mut diag_vy = vec![0.0; n];
            momentum_diagonal(grid, &eta, drag_diag, &mut diag_vx, &mut diag_vy);
            let vjac = VelocityJacobi::from_diagonal(&diag_vx, &diag_vy, self.cfg.diag_floor);

            let mut rhs_pack = Vec::with_capacity(2 * n);
            rhs_pack.extend_from_slice(rhs_x);
            rhs_pack.extend_from_slice(rhs_y);
            {
                let (bx, by) = rhs_pack.split_at_mut(n);
                nullspace::project_velocity(bx, by);
            }
            // Initial guess = current velocity (linear problem, warm start helps).
            let mut x_pack = vec![0.0; 2 * n];
            x_pack[..n].copy_from_slice(vx);
            x_pack[n..].copy_from_slice(vy);

            let cg = ConjugateGradient::new(self.cfg.linear_tol, self.cfg.linear_max_iter);
            let mut tmp_ax = vec![0.0; n];
            let mut tmp_ay = vec![0.0; n];
            let mut matvec = |v: &[f64], out: &mut [f64]| {
                let (vx_in, vy_in) = v.split_at(n);
                let (out_x, out_y) = out.split_at_mut(n);
                apply_momentum(grid, &eta, drag_diag, vx_in, vy_in, &mut tmp_ax, &mut tmp_ay);
                out_x.copy_from_slice(&tmp_ax);
                out_y.copy_from_slice(&tmp_ay);
            };
            let mut precond = |r: &[f64], z: &mut [f64]| vjac.apply(r, z);
            let stats: SolverStats = cg.solve(&mut matvec, &mut precond, &rhs_pack, &mut x_pack);
            trace.linear_iters.push(stats.iterations);
            linear_iters_total = linear_iters_total.saturating_add(stats.iterations as u32);

            // Relaxed update: v_{k+1} = v_k + ω (v_star - v_k).
            let (vs_x, vs_y) = x_pack.split_at(n);
            let omega = self.cfg.relaxation;
            for i in 0..n {
                vx[i] = vx[i] + omega * (vs_x[i] - vx[i]);
                vy[i] = vy[i] + omega * (vs_y[i] - vy[i]);
            }
            nullspace::project_velocity(vx, vy);
            trace.alphas.push(omega);

            // Residual at the new iterate.
            compute_residual(grid, law, drag_diag, cratonic, vx, vy, rhs_x, rhs_y, &mut r_x, &mut r_y);
            let r_norm = vec_norm(&r_x, &r_y);
            trace.residuals.push(r_norm);
            prev_resid = r_norm;
            if r_norm <= abs_tol_eff || r_norm <= self.cfg.rel_tol * r0_norm {
                return NonlinearOutcome::Converged {
                    outer_iters: k + 1,
                    final_residual: r_norm,
                    linear_iters_total,
                    trace,
                };
            }
        }
        NonlinearOutcome::CappedIters {
            max_iters_hit: self.cfg.max_outer_iters,
            last_residual: prev_resid,
            trace,
        }
    }
}
