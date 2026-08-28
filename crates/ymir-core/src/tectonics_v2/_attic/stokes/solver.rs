//! `LinearSolver` trait and a conjugate-gradient implementation.
//!
//! # Step 8.5b parallelisation
//!
//! All `O(N)` vector operations inside the CG inner loop (dot, norm,
//! axpy, and the `p = z + β p` update) are delegated to the
//! [`parallel_reduce`][super::parallel_reduce] helpers or to
//! `rayon::par_iter_mut`. The helpers use a fixed 16-chunk
//! sequential-reduce pattern, so CG is bit-identical across thread
//! counts on a given machine. This cost-dominant inner loop is the
//! main lever behind the ×4 Jacobi / ×2 AMG wallclock targets
//! documented in the Step 8.5b report.
//!
//! # Preconditioner dispatch — Step 8.5a Phase 4.3
//!
//! [`LinearSolverConfig`] selects between Jacobi-CG (default,
//! bit-parity with all pre-8.5a paths) and AMG-CG (Step 8.5a
//! opt-in; operates on the sparse Picard block via
//! `tectonics_v2::stokes::amg::AmgPreconditioner`). Enum dispatch
//! rather than a `dyn Preconditioner` trait object because (i)
//! the enum match compiles to a static jump and (ii) it avoids
//! permanently exposing an abstraction that only has two
//! variants.
//!
//! Reviewer contract (α.1): AmgCG is opt-in per regime. For
//! `step8`-like (η-contrast > 10⁴) callers should stay on
//! JacobiCG until Step 8.5a.2 delivers SA-AMG.
//!
//! Step 0 only ships CG; the trait is the integration point for
//! BiCGSTAB (Step 3) once plastic yielding makes the system
//! non-symmetric. Direct calls to this routine are a violation of the
//! trait discipline — every site that solves a linear system goes
//! through [`LinearSolver::solve`].
//!
//! The solver is matrix-free: callers supply closures for the
//! matrix–vector product and the preconditioner application. The
//! preconditioner closure is expected to internally wrap the
//! null-space projection (see [`super::precond`]), so CG search
//! directions remain orthogonal to the zero-frequency modes at every
//! iteration rather than only at the end.

use rayon::prelude::*;

use super::parallel_reduce::{par_axpy, par_dot, par_norm2};

/// Termination reason for a single linear solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverStatus {
    Converged,
    MaxIterations,
    Breakdown,
}

/// Per-solve statistics. All fields are advisory — the caller uses
/// them for diagnostics, not for retrying.
#[derive(Clone, Copy, Debug)]
pub struct SolverStats {
    pub iterations: usize,
    pub final_residual: f64,
    pub initial_residual: f64,
    pub status: SolverStatus,
}

impl SolverStats {
    pub fn converged(&self) -> bool {
        self.status == SolverStatus::Converged
    }
}

/// Abstraction over iterative linear solvers.
///
/// The closures take slices rather than a dedicated operator type to
/// keep the trait object-safe and to let callers assemble composite
/// operators without introducing extra trait-implementor types. This
/// is the integration point for BiCGSTAB at Step 3 (non-symmetric
/// system introduced by plastic yielding).
pub trait LinearSolver {
    fn solve(
        &self,
        matvec: &mut dyn FnMut(&[f64], &mut [f64]),
        precond: &mut dyn FnMut(&[f64], &mut [f64]),
        b: &[f64],
        x: &mut [f64],
    ) -> SolverStats;
}

/// Preconditioned conjugate gradient (Hestenes–Stiefel).
///
/// Assumes the operator is symmetric positive definite on the
/// orthogonal complement of the null space handled by the
/// preconditioner wrapper.
pub struct ConjugateGradient {
    pub tol: f64,
    pub max_iter: usize,
}

impl ConjugateGradient {
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self { tol, max_iter }
    }
}

impl LinearSolver for ConjugateGradient {
    fn solve(
        &self,
        matvec: &mut dyn FnMut(&[f64], &mut [f64]),
        precond: &mut dyn FnMut(&[f64], &mut [f64]),
        b: &[f64],
        x: &mut [f64],
    ) -> SolverStats {
        let n = b.len();
        assert_eq!(x.len(), n, "x and b must have the same length");

        let mut r = vec![0.0; n];
        let mut z = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut ap = vec![0.0; n];

        // r = b - A x — cell-local, deterministic across thread counts.
        matvec(x, &mut r);
        r.par_iter_mut().zip(b.par_iter()).for_each(|(ri, bi)| *ri = *bi - *ri);

        let b_norm = par_norm2(b).max(1.0);
        let r0_norm = par_norm2(&r);
        if r0_norm <= self.tol * b_norm {
            return SolverStats {
                iterations: 0,
                final_residual: r0_norm,
                initial_residual: r0_norm,
                status: SolverStatus::Converged,
            };
        }

        precond(&r, &mut z);
        p.copy_from_slice(&z);
        let mut rz = par_dot(&r, &z);

        // Step 12 follow-up — cooperative cancel check. The token is
        // bound by the v2 bridge thread before each `run_baseline*`
        // command (see `crates/ymir-viz/src/bridge/v2/thread.rs`); when
        // the UI sets the underlying `AtomicBool`, the next CG iter
        // whose `iter % CANCEL_CHECK_INTERVAL == 0` returns early with
        // `MaxIterations` status. Newton sees this as a non-converged
        // CG and (with its own cancel check at the iter top) returns
        // promptly; the harness's step loop then breaks at the
        // post-step callback. Total Stop-to-return latency drops from
        // one full step (≈ 5–25 s on 64² mantle-on) to a few
        // milliseconds (one CG iter window plus the Newton + step
        // unwind).
        const CANCEL_CHECK_INTERVAL: usize = 16;
        for iter in 1..=self.max_iter {
            if iter % CANCEL_CHECK_INTERVAL == 0 && crate::tectonics_v2::cancel::is_cancelled() {
                return SolverStats {
                    iterations: iter - 1,
                    final_residual: par_norm2(&r),
                    initial_residual: r0_norm,
                    status: SolverStatus::MaxIterations,
                };
            }
            matvec(&p, &mut ap);
            let pap = par_dot(&p, &ap);
            if !pap.is_finite() || pap <= 0.0 {
                return SolverStats {
                    iterations: iter - 1,
                    final_residual: par_norm2(&r),
                    initial_residual: r0_norm,
                    status: SolverStatus::Breakdown,
                };
            }
            let alpha = rz / pap;
            // x += α p   and   r -= α (A p)
            par_axpy(alpha, &p, x);
            par_axpy(-alpha, &ap, &mut r);
            let r_norm = par_norm2(&r);
            if r_norm <= self.tol * b_norm {
                return SolverStats {
                    iterations: iter,
                    final_residual: r_norm,
                    initial_residual: r0_norm,
                    status: SolverStatus::Converged,
                };
            }
            precond(&r, &mut z);
            let rz_new = par_dot(&r, &z);
            if !rz_new.is_finite() {
                return SolverStats {
                    iterations: iter,
                    final_residual: r_norm,
                    initial_residual: r0_norm,
                    status: SolverStatus::Breakdown,
                };
            }
            let beta = rz_new / rz;
            rz = rz_new;
            // p = z + β p   (cell-local, order-independent)
            p.par_iter_mut().zip(z.par_iter()).for_each(|(pi, zi)| *pi = *zi + beta * *pi);
        }

        SolverStats {
            iterations: self.max_iter,
            final_residual: par_norm2(&r),
            initial_residual: r0_norm,
            status: SolverStatus::MaxIterations,
        }
    }
}

/// Preconditioner selector for the CG inner solve.
///
/// Default is `JacobiCG`, preserving the pre-8.5a behaviour
/// byte-for-byte. Switching to `AmgCG(cfg)` selects the Option
/// B' Classical-RS V-cycle preconditioner on the Picard block;
/// the Newton tangent remains matrix-free in the CG matvec.
#[derive(Clone, Copy, Debug)]
pub enum LinearSolverConfig {
    JacobiCG,
    AmgCG(super::amg::AmgConfig),
}

impl Default for LinearSolverConfig {
    fn default() -> Self {
        Self::JacobiCG
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Diagonal SPD system: `A = diag(d)`, solve A x = b.
    /// Exercises the CG contract without coupling to the Stokes
    /// operator.
    #[test]
    fn cg_solves_diagonal_spd_system() {
        let d: Vec<f64> = (1..=16).map(|k| k as f64).collect();
        let b: Vec<f64> = (0..16).map(|k| 2.0 + (k as f64) * 0.1).collect();
        let mut x = vec![0.0; 16];

        let cg = ConjugateGradient::new(1e-12, 100);
        let d_mv = d.clone();
        let d_pc = d.clone();
        let stats = cg.solve(
            &mut |v: &[f64], out: &mut [f64]| {
                for k in 0..v.len() {
                    out[k] = d_mv[k] * v[k];
                }
            },
            &mut |r: &[f64], z: &mut [f64]| {
                for k in 0..r.len() {
                    z[k] = r[k] / d_pc[k];
                }
            },
            &b,
            &mut x,
        );
        assert!(stats.converged(), "stats = {:?}", stats);
        // CG with Jacobi preconditioner on a diagonal system converges in 1 iteration.
        assert!(stats.iterations <= 2);
        for k in 0..16 {
            let expected = b[k] / d[k];
            assert!((x[k] - expected).abs() < 1e-10);
        }
    }

    #[test]
    fn cg_terminates_when_initial_residual_is_small() {
        let cg = ConjugateGradient::new(1e-8, 100);
        let n = 8;
        let d = vec![1.0; n];
        let b = vec![0.0; n];
        let mut x = vec![0.0; n];
        let stats = cg.solve(
            &mut |v: &[f64], out: &mut [f64]| {
                for k in 0..n {
                    out[k] = d[k] * v[k];
                }
            },
            &mut |r: &[f64], z: &mut [f64]| {
                z.copy_from_slice(r);
            },
            &b,
            &mut x,
        );
        assert_eq!(stats.iterations, 0);
        assert!(stats.converged());
    }
}
