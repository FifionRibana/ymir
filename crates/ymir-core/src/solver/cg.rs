//! Preconditioned Conjugate Gradient solver.

use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::stokes::apply_stokes;

/// Workspace buffers for the CG solver.
pub struct CgWorkspace {
    pub r: Vec<f64>,
    pub p: Vec<f64>,
    pub ap: Vec<f64>,
    pub z: Vec<f64>,
    pub precond: Vec<f64>,
}

impl CgWorkspace {
    pub fn new(size: usize) -> Self {
        Self {
            r: vec![0.0; size],
            p: vec![0.0; size],
            ap: vec![0.0; size],
            z: vec![0.0; size],
            precond: vec![1.0; size],
        }
    }
}

/// Result of a CG solve.
pub struct CgResult {
    pub iterations: usize,
    pub residual_norm: f64,
    pub converged: bool,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Solve A·x = b using preconditioned Conjugate Gradient (Jacobi preconditioner).
///
/// `x` is the initial guess and will contain the solution on return.
/// `precond` in `ws` should be set to 1/diag(A) before calling.
#[allow(clippy::needless_range_loop)]
pub fn solve_cg(
    x: &mut [f64],
    b: &[f64],
    eta: &Field2D,
    grid: &StaggeredGrid,
    ws: &mut CgWorkspace,
    max_iter: usize,
    tolerance: f64,
) -> CgResult {
    let n = x.len();

    // r = b - A*x
    apply_stokes(x, eta, grid, &mut ws.ap);
    for i in 0..n {
        ws.r[i] = b[i] - ws.ap[i];
    }

    let b_norm = norm(b).max(1e-14);
    let tol = tolerance * b_norm;

    // z = M^{-1} r
    for i in 0..n {
        ws.z[i] = ws.precond[i] * ws.r[i];
    }

    // p = z
    ws.p.copy_from_slice(&ws.z);

    let mut rz = dot(&ws.r, &ws.z);

    for iter in 0..max_iter {
        let r_norm = norm(&ws.r);
        if r_norm < tol {
            return CgResult {
                iterations: iter,
                residual_norm: r_norm,
                converged: true,
            };
        }

        // ap = A*p
        apply_stokes(&ws.p, eta, grid, &mut ws.ap);

        let pap = dot(&ws.p, &ws.ap);
        if pap < 1e-30 {
            return CgResult {
                iterations: iter,
                residual_norm: norm(&ws.r),
                converged: false,
            };
        }

        let alpha = rz / pap;

        for i in 0..n {
            x[i] += alpha * ws.p[i];
            ws.r[i] -= alpha * ws.ap[i];
        }

        // z = M^{-1} r
        for i in 0..n {
            ws.z[i] = ws.precond[i] * ws.r[i];
        }

        let rz_new = dot(&ws.r, &ws.z);
        let beta = rz_new / rz.max(1e-30);
        rz = rz_new;

        for i in 0..n {
            ws.p[i] = ws.z[i] + beta * ws.p[i];
        }
    }

    CgResult {
        iterations: max_iter,
        residual_norm: norm(&ws.r),
        converged: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::stokes::compute_jacobi_precond;

    fn deterministic_rand(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*state >> 33) as f64 / (1u64 << 31) as f64 - 0.5
    }

    #[test]
    fn cg_converges() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        // Build a sinusoidal RHS
        let k = 2.0 * std::f64::consts::PI;
        let mut b = vec![0.0; nn2];
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                b[j * n + i] = (k * x).sin();
                b[n * n + j * n + i] = (k * y).sin(); // sin has zero mean on [0,1]
            }
        }

        // Project out null space (safety — sin should already have ~zero mean) (constant mode on periodic domain)
        let n2 = n * n;
        let mean_vx: f64 = b[..n2].iter().sum::<f64>() / n2 as f64;
        let mean_vy: f64 = b[n2..].iter().sum::<f64>() / n2 as f64;
        for val in &mut b[..n2] {
            *val -= mean_vx;
        }
        for val in &mut b[n2..] {
            *val -= mean_vy;
        }

        let mut x = vec![0.0; nn2];
        let mut ws = CgWorkspace::new(nn2);
        compute_jacobi_precond(&eta, &grid, &mut ws.precond);

        let result = solve_cg(&mut x, &b, &eta, &grid, &mut ws, 1000, 1e-8);
        assert!(result.converged, "CG did not converge: {:?} iters, residual={}", result.iterations, result.residual_norm);

        // Verify: ||Ax - b|| / ||b|| < 1e-8
        let mut ax = vec![0.0; nn2];
        apply_stokes(&x, &eta, &grid, &mut ax);
        let mut err_sq = 0.0;
        let mut b_sq = 0.0;
        for i in 0..nn2 {
            let e = ax[i] - b[i];
            err_sq += e * e;
            b_sq += b[i] * b[i];
        }
        let rel_err = (err_sq / b_sq).sqrt();
        assert!(rel_err < 1e-8, "CG solution inaccurate: rel_err={rel_err}");
    }

    #[test]
    fn jacobi_reduces_iterations() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let nn2 = 2 * n * n;

        // Variable η — harder to solve
        let mut eta = Field2D::new(n);
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                eta.set(i, j, 1.0 + 9.0 * x); // η varies 1..10
            }
        }

        let mut state = 55u64;
        let mut b: Vec<f64> = (0..nn2).map(|_| deterministic_rand(&mut state)).collect();

        // Project out null space
        let n2 = n * n;
        let mean_vx: f64 = b[..n2].iter().sum::<f64>() / n2 as f64;
        let mean_vy: f64 = b[n2..].iter().sum::<f64>() / n2 as f64;
        for val in &mut b[..n2] {
            *val -= mean_vx;
        }
        for val in &mut b[n2..] {
            *val -= mean_vy;
        }

        // Without preconditioner (identity)
        let mut x_no_prec = vec![0.0; nn2];
        let mut ws_no_prec = CgWorkspace::new(nn2);
        // precond stays at 1.0 (identity)
        let res_no_prec = solve_cg(&mut x_no_prec, &b, &eta, &grid, &mut ws_no_prec, 2000, 1e-6);

        // With Jacobi preconditioner
        let mut x_jac = vec![0.0; nn2];
        let mut ws_jac = CgWorkspace::new(nn2);
        compute_jacobi_precond(&eta, &grid, &mut ws_jac.precond);
        let res_jac = solve_cg(&mut x_jac, &b, &eta, &grid, &mut ws_jac, 2000, 1e-6);

        assert!(
            res_jac.converged,
            "Jacobi CG did not converge"
        );
        assert!(
            res_jac.iterations < res_no_prec.iterations,
            "Jacobi should reduce iterations: {} (Jacobi) vs {} (no precond)",
            res_jac.iterations,
            res_no_prec.iterations,
        );
    }
}
