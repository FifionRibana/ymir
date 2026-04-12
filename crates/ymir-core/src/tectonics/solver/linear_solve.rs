//! Linear solvers: Preconditioned CG and BiCGSTAB.
//!
//! Both solvers take a closure `apply_operator: impl Fn(&[f64], &mut [f64])`
//! for the matrix-vector product, making them reusable for Stokes (Picard)
//! and Jacobian-free Newton-Krylov (JFNK).

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

/// Workspace buffers for BiCGSTAB.
pub struct BiCgStabWorkspace {
    pub r: Vec<f64>,
    pub r_hat: Vec<f64>,
    pub p: Vec<f64>,
    pub v: Vec<f64>,
    pub s: Vec<f64>,
    pub t: Vec<f64>,
    pub p_hat: Vec<f64>,
    pub s_hat: Vec<f64>,
    pub precond: Vec<f64>,
}

impl BiCgStabWorkspace {
    pub fn new(size: usize) -> Self {
        Self {
            r: vec![0.0; size],
            r_hat: vec![0.0; size],
            p: vec![0.0; size],
            v: vec![0.0; size],
            s: vec![0.0; size],
            t: vec![0.0; size],
            p_hat: vec![0.0; size],
            s_hat: vec![0.0; size],
            precond: vec![1.0; size],
        }
    }
}

/// Result of a linear solve (CG or BiCGSTAB).
pub struct LinearSolveResult {
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

/// Solve A·x = b using preconditioned Conjugate Gradient.
///
/// `apply_operator` computes A·v → out. The Jacobi preconditioner
/// `ws.precond` (= 1/diag(A)) must be set before calling.
#[allow(clippy::needless_range_loop)]
pub fn solve_cg(
    x: &mut [f64],
    b: &[f64],
    mut apply_operator: impl FnMut(&[f64], &mut [f64]),
    ws: &mut CgWorkspace,
    max_iter: usize,
    tolerance: f64,
) -> LinearSolveResult {
    let n = x.len();

    // r = b - A*x
    apply_operator(x, &mut ws.ap);
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
            return LinearSolveResult { iterations: iter, residual_norm: r_norm, converged: true };
        }

        // ap = A*p
        apply_operator(&ws.p, &mut ws.ap);

        let pap = dot(&ws.p, &ws.ap);
        if pap < 1e-30 {
            return LinearSolveResult {
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

    LinearSolveResult { iterations: max_iter, residual_norm: norm(&ws.r), converged: false }
}

/// Solve A·x = b using preconditioned BiCGSTAB.
///
/// Works for non-symmetric operators (e.g. JFNK Jacobian approximation).
/// `apply_operator` computes A·v → out. The Jacobi preconditioner
/// `ws.precond` must be set before calling.
#[allow(clippy::needless_range_loop)]
pub fn solve_bicgstab(
    x: &mut [f64],
    b: &[f64],
    mut apply_operator: impl FnMut(&[f64], &mut [f64]),
    ws: &mut BiCgStabWorkspace,
    max_iter: usize,
    tolerance: f64,
) -> LinearSolveResult {
    let n = x.len();

    // r₀ = b - A·x₀
    apply_operator(x, &mut ws.t); // use t as temp
    for i in 0..n {
        ws.r[i] = b[i] - ws.t[i];
    }

    // r̂₀ = r₀
    ws.r_hat.copy_from_slice(&ws.r);

    let b_norm = norm(b).max(1e-14);
    let tol = tolerance * b_norm;

    let mut rho = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;

    ws.v.iter_mut().for_each(|x| *x = 0.0);
    ws.p.iter_mut().for_each(|x| *x = 0.0);

    for iter in 0..max_iter {
        let r_norm = norm(&ws.r);
        if r_norm < tol {
            return LinearSolveResult { iterations: iter, residual_norm: r_norm, converged: true };
        }

        let rho_new = dot(&ws.r_hat, &ws.r);
        if rho_new.abs() < 1e-30 {
            return LinearSolveResult { iterations: iter, residual_norm: r_norm, converged: false };
        }

        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;

        // p = r + β(p - ω·v)
        for i in 0..n {
            ws.p[i] = ws.r[i] + beta * (ws.p[i] - omega * ws.v[i]);
        }

        // p̂ = M⁻¹·p
        for i in 0..n {
            ws.p_hat[i] = ws.precond[i] * ws.p[i];
        }

        // v = A·p̂
        apply_operator(&ws.p_hat, &mut ws.v);

        let r_hat_v = dot(&ws.r_hat, &ws.v);
        if r_hat_v.abs() < 1e-30 {
            return LinearSolveResult { iterations: iter, residual_norm: r_norm, converged: false };
        }
        alpha = rho / r_hat_v;

        // s = r - α·v
        for i in 0..n {
            ws.s[i] = ws.r[i] - alpha * ws.v[i];
        }

        // Early convergence check on s
        let s_norm = norm(&ws.s);
        if s_norm < tol {
            for i in 0..n {
                x[i] += alpha * ws.p_hat[i];
            }
            return LinearSolveResult {
                iterations: iter + 1,
                residual_norm: s_norm,
                converged: true,
            };
        }

        // ŝ = M⁻¹·s
        for i in 0..n {
            ws.s_hat[i] = ws.precond[i] * ws.s[i];
        }

        // t = A·ŝ
        apply_operator(&ws.s_hat, &mut ws.t);

        // ω = <t, s> / <t, t>
        let tt = dot(&ws.t, &ws.t);
        if tt < 1e-30 {
            // t ≈ 0 means A·ŝ ≈ 0 — accept current step
            for i in 0..n {
                x[i] += alpha * ws.p_hat[i];
            }
            return LinearSolveResult {
                iterations: iter + 1,
                residual_norm: s_norm,
                converged: s_norm < tol,
            };
        }
        omega = dot(&ws.t, &ws.s) / tt;

        if omega.abs() < 1e-30 {
            return LinearSolveResult {
                iterations: iter + 1,
                residual_norm: norm(&ws.r),
                converged: false,
            };
        }

        // x = x + α·p̂ + ω·ŝ
        for i in 0..n {
            x[i] += alpha * ws.p_hat[i] + omega * ws.s_hat[i];
        }

        // r = s - ω·t
        for i in 0..n {
            ws.r[i] = ws.s[i] - omega * ws.t[i];
        }
    }

    LinearSolveResult { iterations: max_iter, residual_norm: norm(&ws.r), converged: false }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::solver::field::Field2D;
    use crate::tectonics::solver::grid::StaggeredGrid;
    use crate::tectonics::solver::stokes::{apply_stokes, compute_jacobi_precond};

    fn deterministic_rand(state: &mut u64) -> f64 {
        *state =
            state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (*state >> 33) as f64 / (1u64 << 31) as f64 - 0.5
    }

    fn project_null_space(b: &mut [f64], n: usize) {
        let n2 = n * n;
        let mean_vx: f64 = b[..n2].iter().sum::<f64>() / n2 as f64;
        let mean_vy: f64 = b[n2..].iter().sum::<f64>() / n2 as f64;
        for val in &mut b[..n2] {
            *val -= mean_vx;
        }
        for val in &mut b[n2..] {
            *val -= mean_vy;
        }
    }

    #[test]
    fn cg_converges() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        let k = 2.0 * std::f64::consts::PI;
        let mut b = vec![0.0; nn2];
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                b[j * n + i] = (k * x).sin();
                b[n * n + j * n + i] = (k * y).sin();
            }
        }
        project_null_space(&mut b, n);

        let mut x = vec![0.0; nn2];
        let mut ws = CgWorkspace::new(nn2);
        compute_jacobi_precond(&eta, &grid, &mut ws.precond);

        let result =
            solve_cg(&mut x, &b, |v, out| apply_stokes(v, &eta, &grid, out), &mut ws, 1000, 1e-8);
        assert!(
            result.converged,
            "CG did not converge: {} iters, residual={}",
            result.iterations, result.residual_norm
        );

        // Verify: ||Ax - b|| / ||b|| < 1e-8
        let mut ax = vec![0.0; nn2];
        apply_stokes(&x, &eta, &grid, &mut ax);
        let err: f64 = ax.iter().zip(&b).map(|(a, b)| (a - b).powi(2)).sum();
        let b_sq: f64 = b.iter().map(|v| v * v).sum();
        let rel_err = (err / b_sq).sqrt();
        assert!(rel_err < 1e-8, "CG solution inaccurate: rel_err={rel_err}");
    }

    #[test]
    fn jacobi_reduces_iterations() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let nn2 = 2 * n * n;

        let mut eta = Field2D::new(n);
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                eta.set(i, j, 1.0 + 9.0 * x);
            }
        }

        let mut state = 55u64;
        let mut b: Vec<f64> = (0..nn2).map(|_| deterministic_rand(&mut state)).collect();
        project_null_space(&mut b, n);

        // Without preconditioner
        let mut x_no_prec = vec![0.0; nn2];
        let mut ws_no_prec = CgWorkspace::new(nn2);
        let res_no_prec = solve_cg(
            &mut x_no_prec,
            &b,
            |v, out| apply_stokes(v, &eta, &grid, out),
            &mut ws_no_prec,
            2000,
            1e-6,
        );

        // With Jacobi
        let mut x_jac = vec![0.0; nn2];
        let mut ws_jac = CgWorkspace::new(nn2);
        compute_jacobi_precond(&eta, &grid, &mut ws_jac.precond);
        let res_jac = solve_cg(
            &mut x_jac,
            &b,
            |v, out| apply_stokes(v, &eta, &grid, out),
            &mut ws_jac,
            2000,
            1e-6,
        );

        assert!(res_jac.converged, "Jacobi CG did not converge");
        assert!(
            res_jac.iterations < res_no_prec.iterations,
            "Jacobi should reduce iterations: {} vs {}",
            res_jac.iterations,
            res_no_prec.iterations,
        );
    }

    #[test]
    fn bicgstab_converges_symmetric() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        let k = 2.0 * std::f64::consts::PI;
        let mut b = vec![0.0; nn2];
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                b[j * n + i] = (k * x).sin();
                b[n * n + j * n + i] = (k * y).sin();
            }
        }
        project_null_space(&mut b, n);

        let mut x = vec![0.0; nn2];
        let mut ws = BiCgStabWorkspace::new(nn2);
        compute_jacobi_precond(&eta, &grid, &mut ws.precond);

        let result = solve_bicgstab(
            &mut x,
            &b,
            |v, out| apply_stokes(v, &eta, &grid, out),
            &mut ws,
            1000,
            1e-8,
        );
        assert!(
            result.converged,
            "BiCGSTAB did not converge: {} iters, residual={}",
            result.iterations, result.residual_norm
        );

        // Verify: ||Ax - b|| / ||b|| < 1e-8
        let mut ax = vec![0.0; nn2];
        apply_stokes(&x, &eta, &grid, &mut ax);
        let err: f64 = ax.iter().zip(&b).map(|(a, b)| (a - b).powi(2)).sum();
        let b_sq: f64 = b.iter().map(|v| v * v).sum();
        let rel_err = (err / b_sq).sqrt();
        assert!(rel_err < 1e-8, "BiCGSTAB solution inaccurate: rel_err={rel_err}");
    }

    #[test]
    fn bicgstab_converges_nonsymmetric() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        let mut state = 42u64;
        let mut b: Vec<f64> = (0..nn2).map(|_| deterministic_rand(&mut state)).collect();
        project_null_space(&mut b, n);

        // Non-symmetric operator: Stokes + asymmetric perturbation
        let nonsym_op = |v_in: &[f64], v_out: &mut [f64]| {
            apply_stokes(v_in, &eta, &grid, v_out);
            for i in 0..nn2 {
                v_out[i] += 0.1 * v_in[(i + 7) % nn2];
            }
        };

        let mut x = vec![0.0; nn2];
        let mut ws = BiCgStabWorkspace::new(nn2);
        compute_jacobi_precond(&eta, &grid, &mut ws.precond);

        let result = solve_bicgstab(&mut x, &b, nonsym_op, &mut ws, 2000, 1e-6);
        assert!(
            result.converged,
            "BiCGSTAB should converge on non-symmetric system: {} iters, residual={}",
            result.iterations, result.residual_norm
        );
    }
}
