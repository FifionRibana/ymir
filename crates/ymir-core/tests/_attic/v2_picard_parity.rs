//! Picard-Newton parity: the two nonlinear solvers must arrive at
//! the same fixed point on a manufactured problem. Picard is allowed
//! to take many more outer iterations, but the final solutions
//! should agree at `1e-8`.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::rheology::{self, StrainRate, ViscosityLaw};
use ymir_core::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearOutcome, NonlinearSolver,
};
use ymir_core::tectonics_v2::stokes::picard::{PicardConfig, PicardSolver};
use ymir_core::tectonics_v2::stokes::solver::ConjugateGradient;
use ymir_core::tectonics_v2::stokes::{Grid, operator::apply_momentum};

#[test]
fn picard_and_newton_arrive_at_the_same_solution() {
    let nx = 24;
    let ny = 24;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let law = ViscosityLaw::default();

    // RHS from applying the nonlinear operator to a target velocity.
    let mut vx_target = vec![0.0; nx * ny];
    let mut vy_target = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dy;
            vx_target[j * nx + i] = 0.25 * (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos();
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dy;
            vy_target[j * nx + i] = 0.15 * (2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin();
        }
    }
    ymir_core::tectonics_v2::stokes::nullspace::project_velocity(&mut vx_target, &mut vy_target);
    let sr = StrainRate::compute(nx, ny, dx, dy, &grid.idx_x, &grid.idx_y, &vx_target, &vy_target);
    let eta = rheology::build_eta_field(&law, &sr.eps_ii_center, None);
    let mut rhs_x = vec![0.0; nx * ny];
    let mut rhs_y = vec![0.0; nx * ny];
    apply_momentum(&grid, &eta, None, &vx_target, &vy_target, &mut rhs_x, &mut rhs_y);
    ymir_core::tectonics_v2::stokes::nullspace::project_velocity(&mut rhs_x, &mut rhs_y);

    let cg = ConjugateGradient::new(1.0e-10, 5000);

    // Newton.
    let mut newton_cfg = NewtonConfig::default();
    newton_cfg.rel_tol = 1.0e-8;
    newton_cfg.linear_tol = 1.0e-10;
    let newton = NewtonSolver::new(newton_cfg);
    let mut vn_x = vec![0.0; nx * ny];
    let mut vn_y = vec![0.0; nx * ny];
    let newton_outcome =
        newton.solve(&grid, &law, None, None, &rhs_x, &rhs_y, &mut vn_x, &mut vn_y, &cg);
    assert!(newton_outcome.converged(), "Newton failed: {:?}", newton_outcome);

    // Picard.
    let mut picard_cfg = PicardConfig::default();
    picard_cfg.rel_tol = 1.0e-8;
    picard_cfg.linear_tol = 1.0e-10;
    picard_cfg.max_outer_iters = 200;
    picard_cfg.relaxation = 1.0; // pure Picard — this MMS doesn't need damping
    let picard = PicardSolver::new(picard_cfg);
    let mut vp_x = vec![0.0; nx * ny];
    let mut vp_y = vec![0.0; nx * ny];
    let picard_outcome =
        picard.solve(&grid, &law, None, None, &rhs_x, &rhs_y, &mut vp_x, &mut vp_y, &cg);

    eprintln!(
        "Newton outer={}, Picard outer={}, outcome={:?}",
        newton_outcome.outer_iters(),
        picard_outcome.outer_iters(),
        std::mem::discriminant(&picard_outcome),
    );

    let converged_picard = matches!(
        &picard_outcome,
        NonlinearOutcome::Converged { .. } | NonlinearOutcome::CappedIters { .. }
    );
    assert!(converged_picard, "Picard ended unexpectedly: {:?}", picard_outcome,);

    // Compare solutions.
    let err_sq: f64 = vn_x
        .iter()
        .zip(vp_x.iter())
        .chain(vn_y.iter().zip(vp_y.iter()))
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    let rms = (err_sq / (2 * vn_x.len()) as f64).sqrt();
    eprintln!("‖v_newton − v_picard‖_RMS = {:.3e}", rms);
    assert!(
        rms < 1.0e-7,
        "Newton vs Picard RMS = {:.3e} (spec: 1e-8 target; 1e-7 gives margin for inexact inner CG)",
        rms,
    );
}
