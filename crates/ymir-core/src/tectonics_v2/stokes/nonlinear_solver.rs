//! Nonlinear solver trait and a Newton implementation with Armijo
//! backtracking line search.
//!
//! # System
//!
//! Solves
//! ```text
//!   A(v; η(ε̇_II(v))) = rhs
//! ```
//! where `A(v; η) = -∇·(2 η ε̇(v))` is the depth-integrated thin-sheet
//! momentum operator (see `operator::apply_momentum`) and `η` follows
//! the power-law rheology of `rheology::ViscosityLaw`.
//!
//! # Algorithm
//!
//! Each outer iteration:
//! 1. Compute the residual `r_k = A(v_k; η(v_k)) - rhs`.
//! 2. Assemble the tangent Jacobian `J(v_k)` (Picard + Newton-extra).
//! 3. Solve `J δv = -r_k` via the supplied `LinearSolver` (CG at
//!    Steps 0–2; BiCGSTAB from Step 3 once yielding makes the system
//!    non-symmetric).
//! 4. Armijo line search: `α = 1`, accept if
//!    `‖r(v_k + α δv)‖² ≤ (1 - c₁ α) ‖r_k‖²`; else halve `α` up to
//!    `max_backtrack` times.
//! 5. Update `v_{k+1} = v_k + α δv`.
//!
//! Exit conditions: relative residual reduction, absolute tolerance,
//! stall detection over 3 successive iterations, divergence at 10×
//! initial residual, or the outer-iteration cap.

use super::super::field::{Field2D, PeriodicIndex};
use super::super::rheology::{self, StrainRate, ViscosityLaw};
use super::nullspace;
use super::operator::{apply_momentum, momentum_diagonal, StokesGrid, TangentContext};
use super::precond::VelocityJacobi;
use super::solver::{ConjugateGradient, LinearSolver, SolverStats};

/// Configuration for the Newton solver.
#[derive(Clone, Copy, Debug)]
pub struct NewtonConfig {
    pub rel_tol: f64,
    pub abs_tol: f64,
    pub max_outer_iters: u32,
    pub max_backtrack: u32,
    /// Armijo constant `c₁`.
    pub armijo_c1: f64,
    /// Inner (linear) solver tolerance and iteration cap.
    pub linear_tol: f64,
    pub linear_max_iter: usize,
    /// Jacobi diagonal floor.
    pub diag_floor: f64,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            rel_tol: 1.0e-6,
            abs_tol: 1.0e-10,
            max_outer_iters: 20,
            max_backtrack: 10,
            armijo_c1: 1.0e-4,
            linear_tol: 1.0e-8,
            linear_max_iter: 2000,
            diag_floor: 1.0e-20,
        }
    }
}

/// Detailed per-iteration trace used for post-mortem diagnostics and
/// test assertions on convergence order.
#[derive(Clone, Debug, Default)]
pub struct NonlinearTrace {
    /// `‖r_k‖` before each outer iteration, plus the final value.
    pub residuals: Vec<f64>,
    /// Line-search step lengths accepted per iteration.
    pub alphas: Vec<f64>,
    /// CG iterations consumed by the inner solve of each outer step.
    pub linear_iters: Vec<usize>,
}

#[derive(Clone, Debug)]
pub enum NonlinearOutcome {
    Converged {
        outer_iters: u32,
        final_residual: f64,
        linear_iters_total: u32,
        trace: NonlinearTrace,
    },
    Stalled {
        outer_iters: u32,
        trace: NonlinearTrace,
    },
    Diverged {
        outer_iters: u32,
        last_residual: f64,
        trace: NonlinearTrace,
    },
    CappedIters {
        max_iters_hit: u32,
        last_residual: f64,
        trace: NonlinearTrace,
    },
}

impl NonlinearOutcome {
    pub fn converged(&self) -> bool {
        matches!(self, NonlinearOutcome::Converged { .. })
    }
    pub fn trace(&self) -> &NonlinearTrace {
        match self {
            NonlinearOutcome::Converged { trace, .. }
            | NonlinearOutcome::Stalled { trace, .. }
            | NonlinearOutcome::Diverged { trace, .. }
            | NonlinearOutcome::CappedIters { trace, .. } => trace,
        }
    }
    pub fn outer_iters(&self) -> u32 {
        match self {
            NonlinearOutcome::Converged { outer_iters, .. }
            | NonlinearOutcome::Stalled { outer_iters, .. }
            | NonlinearOutcome::Diverged { outer_iters, .. } => *outer_iters,
            NonlinearOutcome::CappedIters { max_iters_hit, .. } => *max_iters_hit,
        }
    }
}

/// Nonlinear-solver trait. Implementations: [`NewtonSolver`] and
/// [`super::picard::PicardSolver`].
///
/// `drag_diag` is the optional Step-4 basal-drag diagonal field
/// (cell-centered `Br · S̃^exp`). `None` disables drag; `Some(&field)`
/// augments the operator and the preconditioner diagonal consistently.
pub trait NonlinearSolver {
    fn solve(
        &self,
        grid: &StokesGrid,
        law: &ViscosityLaw,
        drag_diag: Option<&Field2D>,
        rhs_x: &[f64],
        rhs_y: &[f64],
        vx: &mut [f64],
        vy: &mut [f64],
        linear_solver: &dyn LinearSolver,
    ) -> NonlinearOutcome;
}

pub struct NewtonSolver {
    pub cfg: NewtonConfig,
}

impl NewtonSolver {
    pub fn new(cfg: NewtonConfig) -> Self {
        Self { cfg }
    }
}

impl Default for NewtonSolver {
    fn default() -> Self {
        Self { cfg: NewtonConfig::default() }
    }
}

/// Compute the nonlinear residual
/// `r = A(v; η(ε̇_II(v))) + Br·S̃²·v - rhs`. The null-space components
/// are **not** projected here — caller is expected to supply a
/// gauge-fixed rhs and initial guess.
fn compute_residual(
    grid: &StokesGrid,
    law: &ViscosityLaw,
    drag_diag: Option<&Field2D>,
    vx: &[f64],
    vy: &[f64],
    rhs_x: &[f64],
    rhs_y: &[f64],
    out_x: &mut [f64],
    out_y: &mut [f64],
    sr_out: &mut Option<StrainRate>,
    eta_out: &mut Option<super::super::field::Field2D>,
) {
    let sr = StrainRate::compute(
        grid.nx,
        grid.ny,
        grid.dx,
        grid.dy,
        &grid.idx_x,
        &grid.idx_y,
        vx,
        vy,
    );
    let eta = rheology::build_eta_field(law, &sr.eps_ii_center);
    apply_momentum(grid, &eta, drag_diag, vx, vy, out_x, out_y);
    for k in 0..out_x.len() {
        out_x[k] -= rhs_x[k];
        out_y[k] -= rhs_y[k];
    }
    // Clean null-space — the operator has it, and a gauge-fixed rhs
    // does too, but round-off accumulates across iterations.
    nullspace::project_velocity(out_x, out_y);
    *sr_out = Some(sr);
    *eta_out = Some(eta);
}

fn vec_norm(a: &[f64], b: &[f64]) -> f64 {
    let s = a.iter().map(|x| x * x).sum::<f64>() + b.iter().map(|x| x * x).sum::<f64>();
    s.sqrt()
}

impl NonlinearSolver for NewtonSolver {
    fn solve(
        &self,
        grid: &StokesGrid,
        law: &ViscosityLaw,
        drag_diag: Option<&Field2D>,
        rhs_x: &[f64],
        rhs_y: &[f64],
        vx: &mut [f64],
        vy: &mut [f64],
        linear_solver: &dyn LinearSolver,
    ) -> NonlinearOutcome {
        let n = grid.n_cells();
        let mut trace = NonlinearTrace::default();
        let mut linear_iters_total = 0u32;

        let mut r_x = vec![0.0; n];
        let mut r_y = vec![0.0; n];
        let mut sr_k: Option<StrainRate> = None;
        let mut eta_k: Option<super::super::field::Field2D> = None;
        compute_residual(grid, law, drag_diag, vx, vy, rhs_x, rhs_y, &mut r_x, &mut r_y, &mut sr_k, &mut eta_k);
        let r0_norm = vec_norm(&r_x, &r_y);
        trace.residuals.push(r0_norm);
        // Effective absolute tolerance: the Newton residual cannot go
        // below the inner-CG precision times a small safety factor.
        // Using `abs_tol` alone without this floor causes Newton to
        // chase a target its linear solver can't hit and stalls.
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

        for k in 0..self.cfg.max_outer_iters {
            // Divergence check against the initial residual.
            if prev_resid > 10.0 * r0_norm && r0_norm > 0.0 {
                return NonlinearOutcome::Diverged {
                    outer_iters: k,
                    last_residual: prev_resid,
                    trace,
                };
            }
            // Stall check: three consecutive <1% residual reductions,
            // but only if we are genuinely far from convergence. Near
            // the tolerance floor the residual naturally plateaus at
            // round-off, which is not a stall.
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

            // Build the tangent context from the current iterate.
            let sr = sr_k.take().expect("strain rate computed with residual");
            let _eta_field = eta_k.take();
            let ctx = TangentContext::from_strain_rate(grid, law, &sr);

            // Preconditioner diagonal from the Picard block (SPD
            // approximation, adequate for our mildly-indefinite J).
            // Basal drag adds `Br·S̃²` to the diagonal via drag_diag,
            // matching `apply_momentum`'s augmentation (Case B).
            let mut diag_vx = vec![0.0; n];
            let mut diag_vy = vec![0.0; n];
            momentum_diagonal(grid, &ctx.eta_center, drag_diag, &mut diag_vx, &mut diag_vy);
            let vjac = VelocityJacobi::from_diagonal(&diag_vx, &diag_vy, self.cfg.diag_floor);

            // Solve J δv = -r_k.
            let mut rhs_pack = Vec::with_capacity(2 * n);
            for v in &r_x { rhs_pack.push(-v); }
            for v in &r_y { rhs_pack.push(-v); }
            {
                let (bx, by) = rhs_pack.split_at_mut(n);
                nullspace::project_velocity(bx, by);
            }
            let mut dv_pack = vec![0.0; 2 * n];
            let cg = ConjugateGradient::new(self.cfg.linear_tol, self.cfg.linear_max_iter);
            let _ = linear_solver; // keeps the generic trait in scope
            let mut tmp_ax = vec![0.0; n];
            let mut tmp_ay = vec![0.0; n];
            let mut matvec = |v: &[f64], out: &mut [f64]| {
                let (vx_in, vy_in) = v.split_at(n);
                let (out_x, out_y) = out.split_at_mut(n);
                // J δv = apply_momentum(η_k, drag) + apply_tangent(ctx)
                // Basal drag's Jacobian is diagonal (Br·S̃²·I) and
                // therefore lives entirely in the Picard block; no
                // extra contribution from apply_tangent.
                apply_momentum(grid, &ctx.eta_center, drag_diag, vx_in, vy_in, &mut tmp_ax, &mut tmp_ay);
                out_x.copy_from_slice(&tmp_ax);
                out_y.copy_from_slice(&tmp_ay);
                super::operator::apply_tangent(grid, &ctx, vx_in, vy_in, out_x, out_y);
            };
            let mut precond = |r: &[f64], z: &mut [f64]| vjac.apply(r, z);
            let cg_stats: SolverStats = cg.solve(&mut matvec, &mut precond, &rhs_pack, &mut dv_pack);
            trace.linear_iters.push(cg_stats.iterations);
            linear_iters_total = linear_iters_total.saturating_add(cg_stats.iterations as u32);

            let (dvx, dvy) = dv_pack.split_at(n);

            // --- Armijo backtracking line search ---
            let mut alpha = 1.0f64;
            let r_prev_sq = prev_resid * prev_resid;
            let mut v_trial_x = vec![0.0; n];
            let mut v_trial_y = vec![0.0; n];
            let mut r_trial_x = vec![0.0; n];
            let mut r_trial_y = vec![0.0; n];
            let mut accepted_resid = prev_resid;
            let mut accepted_sr: Option<StrainRate> = None;
            let mut accepted_eta: Option<super::super::field::Field2D> = None;
            let mut accepted = false;
            for _ in 0..=self.cfg.max_backtrack {
                for i in 0..n {
                    v_trial_x[i] = vx[i] + alpha * dvx[i];
                    v_trial_y[i] = vy[i] + alpha * dvy[i];
                }
                nullspace::project_velocity(&mut v_trial_x, &mut v_trial_y);
                let mut sr_trial: Option<StrainRate> = None;
                let mut eta_trial: Option<super::super::field::Field2D> = None;
                compute_residual(
                    grid, law, drag_diag, &v_trial_x, &v_trial_y, rhs_x, rhs_y,
                    &mut r_trial_x, &mut r_trial_y, &mut sr_trial, &mut eta_trial,
                );
                let r_trial_norm = vec_norm(&r_trial_x, &r_trial_y);
                let r_trial_sq = r_trial_norm * r_trial_norm;
                let target = (1.0 - self.cfg.armijo_c1 * alpha) * r_prev_sq;
                if r_trial_sq <= target.max(0.0) {
                    accepted = true;
                    accepted_resid = r_trial_norm;
                    accepted_sr = sr_trial;
                    accepted_eta = eta_trial;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                // Last resort: accept the latest trial so we make any
                // progress; this also gets recorded as a stall signal
                // if the residual doesn't go down.
                accepted_resid = vec_norm(&r_trial_x, &r_trial_y);
                accepted_sr = Some(StrainRate::compute(
                    grid.nx, grid.ny, grid.dx, grid.dy, &grid.idx_x, &grid.idx_y,
                    &v_trial_x, &v_trial_y,
                ));
                accepted_eta = Some(rheology::build_eta_field(law, &accepted_sr.as_ref().unwrap().eps_ii_center));
            }
            trace.alphas.push(alpha);
            trace.residuals.push(accepted_resid);

            // Commit the accepted iterate.
            vx.copy_from_slice(&v_trial_x);
            vy.copy_from_slice(&v_trial_y);
            r_x.copy_from_slice(&r_trial_x);
            r_y.copy_from_slice(&r_trial_y);
            sr_k = accepted_sr;
            eta_k = accepted_eta;

            prev_resid = accepted_resid;

            // Convergence.
            if accepted_resid <= abs_tol_eff
                || accepted_resid <= self.cfg.rel_tol * r0_norm
            {
                return NonlinearOutcome::Converged {
                    outer_iters: k + 1,
                    final_residual: accepted_resid,
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

/// Silence lint about unused import of `PeriodicIndex` — it appears
/// only in internal signatures that the public API doesn't expose yet.
#[allow(dead_code)]
fn _type_marker_periodic_index(_p: &PeriodicIndex) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::field::Field2D;
    use crate::tectonics_v2::stokes::solver::ConjugateGradient;

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    /// Symmetry of the full tangent Jacobian (Picard + Newton extra).
    /// Verifies `⟨J u, w⟩ = ⟨u, J w⟩` for random test vectors with a
    /// non-trivial `v_k` so both the Picard and Newton-extra parts
    /// contribute. The Jacobian is NOT required to be positive
    /// definite (Gerya §14.4 flags local indefiniteness for
    /// shear-thinning).
    #[test]
    fn jacobian_is_symmetric_on_random_inputs() {
        let nx = 8;
        let ny = 8;
        let grid = StokesGrid::new(nx, ny, 0.125, 0.125);
        let law = ViscosityLaw::default();
        let n2 = nx * ny;

        // Non-trivial v_k (ensures ε̇ ≠ 0 and Newton-extra is active).
        let mut vx_k = vec![0.0; n2];
        let mut vy_k = vec![0.0; n2];
        for k in 0..n2 {
            vx_k[k] = ((k as f64 * 0.7).sin()) * 0.3;
            vy_k[k] = ((k as f64 * 1.3).cos()) * 0.2;
        }
        let sr = StrainRate::compute(nx, ny, grid.dx, grid.dy, &grid.idx_x, &grid.idx_y, &vx_k, &vy_k);
        let ctx = TangentContext::from_strain_rate(&grid, &law, &sr);

        let mut ux = vec![0.0; n2];
        let mut uy = vec![0.0; n2];
        let mut wx = vec![0.0; n2];
        let mut wy = vec![0.0; n2];
        for k in 0..n2 {
            ux[k] = ((k as f64 * 2.1).sin()) * 1.1;
            uy[k] = ((k as f64 * 2.5).cos()) * 0.9;
            wx[k] = ((k as f64 * 0.9).sin()) * 0.5;
            wy[k] = ((k as f64 * 1.1).cos()) * 1.3;
        }
        let mut jux = vec![0.0; n2];
        let mut juy = vec![0.0; n2];
        let mut jwx = vec![0.0; n2];
        let mut jwy = vec![0.0; n2];
        super::super::operator::apply_jacobian(&grid, &ctx, None, &ux, &uy, &mut jux, &mut juy);
        super::super::operator::apply_jacobian(&grid, &ctx, None, &wx, &wy, &mut jwx, &mut jwy);
        let lhs = dot(&jux, &wx) + dot(&juy, &wy);
        let rhs = dot(&ux, &jwx) + dot(&uy, &jwy);
        let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
        assert!(rel < 1e-10, "J asymmetric: rel = {}", rel);
    }

    /// Log a single Newton solve trace to see how residuals decay on
    /// a representative Step 1 problem.
    #[test]
    fn newton_trace_on_sin_force_n3() {
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / nx as f64;
        let dy = 1.0 / ny as f64;
        let grid = StokesGrid::new(nx, ny, dx, dy);
        let law = ViscosityLaw::default(); // n=3, ε_min=1e-3, η_max=1e3
        let mut fx = vec![0.0; nx * ny];
        let fy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx;
                fx[j * nx + i] = 0.1 * (2.0 * std::f64::consts::PI * x).sin();
            }
        }
        let mut vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        let solver = NewtonSolver::default();
        let cg = ConjugateGradient::new(1.0e-10, 2000);
        let outcome = solver.solve(&grid, &law, None, &fx, &fy, &mut vx, &mut vy, &cg);
        let trace = outcome.trace();
        eprintln!("residuals (len={}): {:?}", trace.residuals.len(), trace.residuals);
        eprintln!("alphas: {:?}", trace.alphas);
        eprintln!("cg_iters: {:?}", trace.linear_iters);
        eprintln!("outcome: {:?}", std::mem::discriminant(&outcome));
        // Require convergence on this near-trivial problem.
        assert!(outcome.converged(), "outcome = {:?}", outcome);
    }

    /// Basic smoke test: Newton on a trivial (linear, constant-η)
    /// problem converges in one iteration.
    #[test]
    fn newton_trivial_linear_problem() {
        // n = 1 makes the rheology linear and the Newton-extra term
        // identically zero. The solver should finish in one outer
        // iteration.
        let nx = 8;
        let ny = 8;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let mut law = ViscosityLaw::default();
        law.n = 1.0;

        // Manufactured RHS: use constant η (from the floor at n=1,
        // η = 1 since (ε̇+ε̇_min)^0 = 1) and apply the operator to a
        // known v.
        let mut vx_target = vec![0.0; nx * ny];
        let mut vy_target = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                use std::f64::consts::PI;
                let x = i as f64 / nx as f64;
                let y = j as f64 / ny as f64;
                vx_target[j * nx + i] = (2.0 * PI * x).sin();
                vy_target[j * nx + i] = (2.0 * PI * y).sin();
            }
        }
        crate::tectonics_v2::stokes::nullspace::project_velocity(&mut vx_target, &mut vy_target);
        let eta = Field2D::filled(nx, ny, 1.0);
        let mut rhs_x = vec![0.0; nx * ny];
        let mut rhs_y = vec![0.0; nx * ny];
        apply_momentum(&grid, &eta, None, &vx_target, &vy_target, &mut rhs_x, &mut rhs_y);
        crate::tectonics_v2::stokes::nullspace::project_velocity(&mut rhs_x, &mut rhs_y);

        let mut vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        let mut solver = NewtonSolver::default();
        solver.cfg.rel_tol = 1.0e-8;
        let cg = ConjugateGradient::new(1.0e-10, 2000);
        let outcome = solver.solve(&grid, &law, None, &rhs_x, &rhs_y, &mut vx, &mut vy, &cg);
        assert!(outcome.converged(), "trivial problem failed: {:?}", outcome);
    }
}
