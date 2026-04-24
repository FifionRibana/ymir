//! `LinearSolver` trait and a conjugate-gradient implementation.
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

        // r = b - A x
        matvec(x, &mut r);
        for k in 0..n {
            r[k] = b[k] - r[k];
        }

        let b_norm = norm2(b).max(1.0);
        let r0_norm = norm2(&r);
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
        let mut rz = dot(&r, &z);

        for iter in 1..=self.max_iter {
            matvec(&p, &mut ap);
            let pap = dot(&p, &ap);
            if !pap.is_finite() || pap <= 0.0 {
                return SolverStats {
                    iterations: iter - 1,
                    final_residual: norm2(&r),
                    initial_residual: r0_norm,
                    status: SolverStatus::Breakdown,
                };
            }
            let alpha = rz / pap;
            for k in 0..n {
                x[k] += alpha * p[k];
                r[k] -= alpha * ap[k];
            }
            let r_norm = norm2(&r);
            if r_norm <= self.tol * b_norm {
                return SolverStats {
                    iterations: iter,
                    final_residual: r_norm,
                    initial_residual: r0_norm,
                    status: SolverStatus::Converged,
                };
            }
            precond(&r, &mut z);
            let rz_new = dot(&r, &z);
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
            for k in 0..n {
                p[k] = z[k] + beta * p[k];
            }
        }

        SolverStats {
            iterations: self.max_iter,
            final_residual: norm2(&r),
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

#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[inline]
pub fn norm2(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
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
