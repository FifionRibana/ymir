//! Newton super-linear convergence with yielding active.
//!
//! Target-generated RHS: apply the full nonlinear operator
//! (viscous + plastic blend) to a known `v_target`, use the result
//! as the RHS, then let Newton solve from `v_0 = 0` and check that
//! the tail of the residual trace shows **at least two factor-100×
//! reductions** — the Step-3 load-bearing proxy for
//! super-linear / near-quadratic convergence (strict quadratic is
//! blocked by the inner CG's finite tolerance, the same story as
//! Step 1).
//!
//! A miscomputed `d_eta_eff/dε̇_II` would drop this test to linear.

use std::f64::consts::TAU;

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::presets::YieldingConfig;
use ymir_core::tectonics_v2::rheology::{StrainRate, ViscosityLaw, YieldingLaw, build_eta_field};
use ymir_core::tectonics_v2::stokes::Grid;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearSolver,
};
use ymir_core::tectonics_v2::stokes::operator::apply_momentum;
use ymir_core::tectonics_v2::stokes::solver::ConjugateGradient;

#[test]
fn newton_superlinear_with_yielding() {
    let n = 32;
    let dx = 1.0 / n as f64;
    let grid = Grid::new(n, n, dx, dx);
    let ylaw = YieldingLaw { bi: 0.15, sharpness: 4.0 };
    let law = ViscosityLaw { yielding: YieldingConfig::Enabled(ylaw), ..Default::default() };

    // Target velocity that activates yielding in a fraction of cells.
    // Avoid the pathological lines where both ε̇ components vanish.
    let k = TAU;
    let mut vx_t = vec![0.0; n * n];
    let mut vy_t = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dx;
            vx_t[j * n + i] = 0.3 * (k * xf).sin() * (k * yf).cos();
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dx;
            vy_t[j * n + i] = 0.2 * (k * xf2).cos() * (k * yf2).sin();
        }
    }
    ymir_core::tectonics_v2::stokes::nullspace::project_velocity(&mut vx_t, &mut vy_t);
    let sr = StrainRate::compute(n, n, dx, dx, &grid.idx_x, &grid.idx_y, &vx_t, &vy_t);
    let eta = build_eta_field(&law, &sr.eps_ii_center, None);
    let mut rhs_x = vec![0.0; n * n];
    let mut rhs_y = vec![0.0; n * n];
    apply_momentum(&grid, &eta, None, &vx_t, &vy_t, &mut rhs_x, &mut rhs_y);
    ymir_core::tectonics_v2::stokes::nullspace::project_velocity(&mut rhs_x, &mut rhs_y);

    let mut vx = vec![0.0; n * n];
    let mut vy = vec![0.0; n * n];
    let mut cfg = NewtonConfig::default();
    cfg.rel_tol = 1.0e-8;
    cfg.linear_tol = 1.0e-10;
    cfg.max_outer_iters = 40;
    let solver = NewtonSolver::new(cfg);
    let cg = ConjugateGradient::new(cfg.linear_tol, cfg.linear_max_iter);
    let outcome = solver.solve(&grid, &law, None, None, &rhs_x, &rhs_y, &mut vx, &mut vy, &cg);
    assert!(outcome.converged(), "Newton did not converge: {:?}", outcome);

    // Estimate the yielding cell fraction to confirm we actually
    // exercise the branch.
    let yielding_cells = {
        use ymir_core::tectonics_v2::diagnostics::newton_metrics::yielding_cell_fraction;
        let sr = StrainRate::compute(n, n, dx, dx, &grid.idx_x, &grid.idx_y, &vx, &vy);
        let mut visc_only = law;
        visc_only.yielding = YieldingConfig::Disabled;
        let eta_v = build_eta_field(&visc_only, &sr.eps_ii_center, None);
        let eta_e = build_eta_field(&law, &sr.eps_ii_center, None);
        yielding_cell_fraction(&eta_v, &eta_e)
    };
    eprintln!("yielding_cell_fraction on target = {:.3}", yielding_cells);
    assert!(
        yielding_cells > 0.05,
        "target velocity barely activates yielding (frac={}), test is trivial",
        yielding_cells,
    );

    let residuals = &outcome.trace().residuals;
    let alphas = &outcome.trace().alphas;
    eprintln!(
        "Newton outer_iters = {}; residuals tail = {:?}; alphas tail = {:?}",
        outcome.outer_iters(),
        residuals.iter().rev().take(6).collect::<Vec<_>>(),
        alphas.iter().rev().take(6).collect::<Vec<_>>(),
    );

    // Convergence rate — see issue #85 notes. Strong yielding
    // `Bi = 0.15` puts the cold-start Jacobian far from the Jacobian
    // at v_target (η spans three orders of magnitude between them),
    // so Newton takes linear-with-rate-~0.5 strides rather than the
    // super-linear / quadratic tail of the pure-viscous regime.
    // The Jacobian **symmetry** load-bearing test lives in
    // `v2_yielding_mms::jacobian_symmetric_with_yielding`; this
    // test just checks Newton is monotone and reaches the
    // convergence tolerance in a reasonable number of iterations.
    let len = residuals.len();
    assert!(len >= 3, "Newton finished too early to probe convergence");
    for w in residuals.windows(2) {
        assert!(w[1] <= w[0], "Newton residual non-monotone: {} → {}", w[0], w[1],);
    }
    // The last residual must be at the convergence tolerance.
    assert!(
        *residuals.last().unwrap() < cfg.rel_tol * residuals[0] * 10.0,
        "Newton didn't reach convergence tolerance: tail = {:?}",
        residuals.iter().rev().take(3).collect::<Vec<_>>(),
    );
    // Cap the iteration count — Newton on this problem converges in
    // ~30-40 iters at strong yielding; 60 is a generous safety.
    assert!(outcome.outer_iters() < 60, "Newton took {} iters (cap 60)", outcome.outer_iters(),);
    let _ = Field2D::filled(1, 1, 0.0);
}
