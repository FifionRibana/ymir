//! MMS for the slab-mass ODE in isolation (no Stokes solve).
//!
//! The ODE `∂m̃/∂t̃ = Q̃ − m̃/τ̃` with constant `Q̃ = Q₀` and
//! initial `m(0) = 0` has the exact solution
//!
//! ```text
//!   m(t) = Q₀ · τ · (1 − e^{-t/τ})
//! ```
//!
//! and the half-life from `m = 0` to `m = ½ · Q₀·τ` is `τ · ln 2`.
//!
//! Acceptance criteria from the Step 7 spec:
//! - Convergence to `Q₀·τ` at `t = 5·τ` within 1%.
//! - Empirical half-life matches `τ · ln 2` within 5%.

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::slab::SlabState;

#[test]
fn converges_to_q_times_tau() {
    let nx = 8;
    let ny = 8;
    let tau = 0.5;
    let dt = 0.01;
    let q0 = 0.3;

    let mut state = SlabState::new_zero(nx, ny);
    let mut q = Field2D::new(nx, ny);
    for v in q.data_mut().iter_mut() {
        *v = q0;
    }

    for _ in 0..250 {
        state.step_ode(&q, dt, tau);
    }

    let asymptote = q0 * tau;
    for &m in state.m().data().iter() {
        let rel = (m - asymptote).abs() / asymptote;
        assert!(rel < 0.01, "m = {}, asymptote = {}, rel err = {}", m, asymptote, rel);
    }
}

#[test]
fn half_life_matches_tau_ln_two() {
    let nx = 4;
    let ny = 4;
    let tau = 0.5;
    let dt = 0.001;
    let q0 = 1.0;

    let mut state = SlabState::new_zero(nx, ny);
    let mut q = Field2D::new(nx, ny);
    for v in q.data_mut().iter_mut() {
        *v = q0;
    }

    // Analytic half-"approach" time: m(t) = Q₀·τ·(1 − e^{-t/τ}).
    // The instant m(t) = ½ · Q₀·τ is when e^{-t/τ} = ½,
    // i.e. t = τ · ln 2.
    let target = 0.5 * q0 * tau;
    let analytic_half = tau * std::f64::consts::LN_2;

    let mut steps = 0usize;
    while state.m().data()[0] < target && steps < 10_000 {
        state.step_ode(&q, dt, tau);
        steps += 1;
    }
    let observed = steps as f64 * dt;
    let rel = (observed - analytic_half).abs() / analytic_half;
    assert!(
        rel < 0.05,
        "observed half-life = {} (τ·ln2 = {}), rel err = {}",
        observed,
        analytic_half,
        rel,
    );
}

/// Pure decay (Q = 0) from m₀ = 1.0: the field should follow
/// `m(t) = e^{-t/τ}` cell-wise. We probe at three times against
/// the analytic curve.
#[test]
fn pure_decay_matches_exponential() {
    let nx = 4;
    let ny = 4;
    let tau = 0.5;
    let dt = 0.005;

    let mut state = SlabState::new_zero(nx, ny);
    for v in state.m_mut().data_mut().iter_mut() {
        *v = 1.0;
    }
    let q = Field2D::new(nx, ny); // zero

    let check_times: [f64; 3] = [0.25, 0.5, 1.0]; // in units of τ (t = 0.125, 0.25, 0.5).
    let mut step = 0usize;
    for &t_ratio in &check_times {
        let t: f64 = t_ratio * tau;
        let n_target: usize = (t / dt).round() as usize;
        while step < n_target {
            state.step_ode(&q, dt, tau);
            step += 1;
        }
        let analytic = (-t / tau).exp();
        let observed = state.m().data()[0];
        let rel = (observed - analytic).abs() / analytic;
        // Forward Euler has first-order error; with dt/τ = 0.01
        // we expect < 1% at t ≤ τ.
        assert!(rel < 0.01, "t={}: observed={}, analytic={}, rel={}", t, observed, analytic, rel,);
    }
}
