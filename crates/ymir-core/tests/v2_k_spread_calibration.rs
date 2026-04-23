//! Robustness of the `k_spread` bisection harness (Step 5).
//!
//! The full convergence-on-real-simulations case belongs in the
//! physics report (the iterations are printed there). These tests
//! pin the harness's three behaviours:
//!
//! 1. Nominal convergence — a monotone-increasing response lets the
//!    bisection land inside the target band within `max_iters`.
//! 2. Forced `MaxItersReached` — a narrow target band with
//!    `max_iters = 2` (insufficient to bisect the bracket down to
//!    the band width) returns `Err(MaxItersReached { .. })`
//!    deterministically.
//! 3. `OutOfRange` — a target above the bracket's maximum achievable
//!    response returns `Err(OutOfRange { .. })` without looping.
//!
//! All three use **synthetic responder closures** so the test
//! doesn't run real simulations; it verifies the bisection logic
//! alone.

use ymir_core::tectonics_v2::boundaries::{
    calibrate_k_spread, CalibrationError, KSpreadCalibration, K_SPREAD_BRACKET,
};

#[test]
fn nominal_convergence_on_linear_response() {
    let cfg = KSpreadCalibration::step5_default();
    let result = calibrate_k_spread(&cfg, |k| 0.4 * k).expect("should converge");
    assert!(result.final_s_oceanic_mean >= 0.18);
    assert!(result.final_s_oceanic_mean <= 0.22);
    assert!(!result.iterations.is_empty());
}

#[test]
fn max_iters_with_tight_band_terminates_as_max_iters_reached() {
    // Narrow target band + straddling response: bisection cannot
    // land inside the band within `max_iters = 2`.
    let cfg = KSpreadCalibration {
        bracket: K_SPREAD_BRACKET,
        target_range: (0.199, 0.201),
        convergence_tol: 0.001,
        simulation_time: 3.0,
        max_iters: 2,
    };
    let r = calibrate_k_spread(&cfg, |k| 0.2 * k + 0.03);
    assert!(
        matches!(r, Err(CalibrationError::MaxItersReached { .. })),
        "expected MaxItersReached, got {:?}",
        r,
    );
}

#[test]
fn out_of_range_target_above_bracket_max_response() {
    // Response `s(k) = 0.1 k` → at k = 1.0, s = 0.1 (well below
    // target 0.99–1.01). The bisection must not loop: after the
    // pre-bracket check it should return `OutOfRange`.
    let cfg = KSpreadCalibration {
        bracket: K_SPREAD_BRACKET,
        target_range: (0.99, 1.01),
        convergence_tol: 0.01,
        simulation_time: 3.0,
        max_iters: 20,
    };
    let r = calibrate_k_spread(&cfg, |k| 0.1 * k);
    assert!(
        matches!(r, Err(CalibrationError::OutOfRange { .. })),
        "expected OutOfRange, got {:?}",
        r,
    );
}
