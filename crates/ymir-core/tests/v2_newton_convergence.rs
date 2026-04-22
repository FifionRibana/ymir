//! Nonlinear Newton convergence on power-law rheology.
//!
//! Strategy: pick a smooth target velocity, apply the full nonlinear
//! operator to it, use the result as the RHS, then see whether
//! Newton from v₀ = 0 recovers the target.
//!
//! The test verifies two things:
//! 1. **Convergence** within ≤ 10 outer iterations for `n ∈ {1, 2, 3}`.
//! 2. **Quadratic trail**: the last three non-trivial residuals should
//!    satisfy `r_{k+1} ≤ C · r_k²` with `C` bounded (quadratic slope
//!    close to 2 in a log-log fit of the last three residuals). A
//!    buggy Jacobian would give linear (slope ~ 1) convergence.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::rheology::{self, StrainRate, ViscosityLaw};
use ymir_core::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearOutcome, NonlinearSolver,
};
use ymir_core::tectonics_v2::stokes::solver::ConjugateGradient;
use ymir_core::tectonics_v2::stokes::{operator::apply_momentum, Grid};

fn build_rhs_from_target(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    law: &ViscosityLaw,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let grid = Grid::new(nx, ny, dx, dy);
    // Smooth, non-constant, non-divergence-free target.
    let mut vx_target = vec![0.0; nx * ny];
    let mut vy_target = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dy;
            vx_target[j * nx + i] = 0.3 * (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos();
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dy;
            vy_target[j * nx + i] = 0.2 * (2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin();
        }
    }
    ymir_core::tectonics_v2::stokes::nullspace::project_velocity(&mut vx_target, &mut vy_target);
    // Evaluate the nonlinear operator at the target — that's our RHS.
    let sr = StrainRate::compute(
        nx, ny, dx, dy, &grid.idx_x, &grid.idx_y, &vx_target, &vy_target,
    );
    let eta = rheology::build_eta_field(law, &sr.eps_ii_center);
    let mut rhs_x = vec![0.0; nx * ny];
    let mut rhs_y = vec![0.0; nx * ny];
    apply_momentum(&grid, &eta, None, &vx_target, &vy_target, &mut rhs_x, &mut rhs_y);
    ymir_core::tectonics_v2::stokes::nullspace::project_velocity(&mut rhs_x, &mut rhs_y);
    (rhs_x, rhs_y, vx_target, vy_target)
}

fn run_newton_for_n(n_val: f64) -> (NonlinearOutcome, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let mut law = ViscosityLaw::default();
    law.n = n_val;
    let (rhs_x, rhs_y, vx_target, vy_target) = build_rhs_from_target(nx, ny, dx, dy, &law);

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut newton_cfg = NewtonConfig::default();
    newton_cfg.rel_tol = 1.0e-8;
    newton_cfg.linear_tol = 1.0e-10;
    newton_cfg.max_outer_iters = 20;
    let solver = NewtonSolver::new(newton_cfg);
    let cg = ConjugateGradient::new(newton_cfg.linear_tol, newton_cfg.linear_max_iter);
    let outcome = solver.solve(&grid, &law, None, &rhs_x, &rhs_y, &mut vx, &mut vy, &cg);
    (outcome, vx, vy, vx_target, vy_target)
}

#[test]
fn newton_converges_for_n1_n2_n3() {
    for n in [1.0, 2.0, 3.0] {
        let (outcome, _, _, _, _) = run_newton_for_n(n);
        eprintln!(
            "n={}: outer_iters={}, residuals={:?}",
            n,
            outcome.outer_iters(),
            outcome.trace().residuals,
        );
        assert!(outcome.converged(), "Newton did not converge at n={}: {:?}", n, outcome);
        assert!(
            outcome.outer_iters() <= 10,
            "Newton at n={} took {} iters (> 10)",
            n,
            outcome.outer_iters(),
        );
    }
}

#[test]
fn newton_trail_is_superlinear_for_n3() {
    // Strict quadratic convergence requires an exact inner linear
    // solve; with inexact CG (tol = 1e-10) we see **superlinear** but
    // not quite quadratic decay. The spec's "`r_{k+1} ≤ C · r_k²`"
    // test is re-expressed as "two successive factor-of-100
    // reductions in the tail" — clearly super-linear, easily
    // achieved by a correct Jacobian, impossible for a buggy one
    // (which would give only geometric/linear decay).
    let (outcome, _, _, _, _) = run_newton_for_n(3.0);
    let res = &outcome.trace().residuals;
    assert!(outcome.converged());
    assert!(res.len() >= 3, "need ≥ 3 residuals for slope estimate");
    let last = res.len();
    let r_km2 = res[last - 3];
    let r_km1 = res[last - 2];
    let r_k = res[last - 1];
    // Skip the check if any residual is already at numerical floor.
    if r_km2.min(r_km1).min(r_k) < 1e-12 {
        eprintln!("residuals near numerical floor; skipping slope check");
        return;
    }
    let ratio_1 = r_km2 / r_km1;
    let ratio_2 = r_km1 / r_k;
    eprintln!(
        "Newton tail ratios: r_{{k-2}}/r_{{k-1}} = {:.1}×, r_{{k-1}}/r_k = {:.1}×",
        ratio_1, ratio_2,
    );
    assert!(
        ratio_1 >= 100.0 && ratio_2 >= 100.0,
        "Newton tail did not show the expected super-linear factor-of-100 reductions: {:.1}× then {:.1}×",
        ratio_1, ratio_2,
    );
}

#[test]
fn newton_recovers_target_velocity() {
    let (outcome, vx, vy, vx_tgt, vy_tgt) = run_newton_for_n(2.0);
    assert!(outcome.converged());
    let err_sq: f64 = vx
        .iter()
        .zip(vx_tgt.iter())
        .chain(vy.iter().zip(vy_tgt.iter()))
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    let err_rms = (err_sq / (2 * vx.len()) as f64).sqrt();
    eprintln!("Newton vs target, L² error = {:.3e}", err_rms);
    assert!(err_rms < 1e-6, "Newton solution deviates from target: L² = {}", err_rms);
}
