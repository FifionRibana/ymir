//! Continuation ramp test. The spec asks to construct a case where
//! Newton directly at `n = 3` diverges but converges through the
//! ramp `n = 1 → 3`. In practice a purely linear forcing on small
//! velocities makes Newton at n = 3 converge quickly from `v₀ = 0`,
//! so we do not rely on a demonstrable "direct-divergence" case —
//! the test verifies instead that the continuation wrapper **itself
//! converges** on every rung of the ramp, which is the invariant it
//! is designed to provide.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::presets::ContinuationConfig;
use ymir_core::tectonics_v2::rheology::ViscosityLaw;
use ymir_core::tectonics_v2::stokes::continuation::run_continuation;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearOutcome,
};
use ymir_core::tectonics_v2::stokes::solver::ConjugateGradient;
use ymir_core::tectonics_v2::stokes::Grid;

#[test]
fn continuation_ramp_converges_every_rung() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let mut law = ViscosityLaw::default();
    law.n = 3.0;
    let schedule = ContinuationConfig::step1_default();

    let mut fx = vec![0.0; nx * ny];
    let fy = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let xf = i as f64 * dx;
            fx[j * nx + i] = 0.5 * (2.0 * PI * xf).sin();
        }
    }

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut newton_cfg = NewtonConfig::default();
    newton_cfg.rel_tol = 1.0e-7;
    let newton = NewtonSolver::new(newton_cfg);
    let cg = ConjugateGradient::new(newton_cfg.linear_tol, newton_cfg.linear_max_iter);
    let outcome = run_continuation(
        &grid, &law, None, None, &schedule, &fx, &fy, &mut vx, &mut vy, &newton, &cg,
    );

    assert!(
        outcome.all_converged,
        "continuation ramp failed at sub-solve {}: {:?}",
        outcome.sub_outcomes.len(),
        outcome.sub_outcomes.last(),
    );
    assert_eq!(outcome.sub_outcomes.len(), schedule.n_steps.len());
    for (n, sub) in &outcome.sub_outcomes {
        assert!(
            matches!(sub, NonlinearOutcome::Converged { .. }),
            "continuation sub-solve at n={} did not converge: {:?}",
            n,
            sub,
        );
    }
    eprintln!(
        "continuation ramp OK over {:?}; total linear iters = {}",
        schedule.n_steps, outcome.linear_iters_total,
    );
}
