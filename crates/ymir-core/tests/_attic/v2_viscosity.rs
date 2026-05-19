//! Integration tests for the `ViscosityLaw` and smooth-saturate parity.

use ymir_core::tectonics_v2::rheology::{smooth_saturate, ViscosityLaw};

fn approx_rel(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-12)
}

#[test]
fn smooth_saturate_parity_with_legacy() {
    // Legacy `tectonics/solver/smooth.rs::smooth_saturate` hardcodes
    // k = 4. Verify byte-equivalent output on a log-spaced grid.
    let x_max = 1.0e4;
    let xs: [f64; 9] = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5];
    for x in xs {
        let legacy = x / (1.0 + (x / x_max).powf(4.0)).powf(1.0 / 4.0);
        let here = smooth_saturate(x, x_max, 4.0);
        assert!((legacy - here).abs() < 1e-14);
    }
}

#[test]
fn eta_is_bounded_above_by_cap() {
    let law = ViscosityLaw::default();
    for x in [1e-6, 1e-3, 1.0, 1e6] {
        assert!(law.eta_effective(x) <= law.eta_max_cap);
    }
}

#[test]
fn derivative_matches_finite_difference_on_log_grid() {
    let law = ViscosityLaw::default();
    let h = 1.0e-7;
    for k in -5..6 {
        let x: f64 = 10.0_f64.powi(k);
        let analytic = law.d_eta_effective_d_eps_ii(x);
        let fd = (law.eta_effective(x + h) - law.eta_effective(x - h)) / (2.0 * h);
        assert!(
            approx_rel(analytic, fd, 1.0e-4),
            "d_eta at x={}: analytic={}, fd={}",
            x,
            analytic,
            fd,
        );
    }
}

#[test]
fn monotonicity_shear_thinning() {
    let law = ViscosityLaw::default();
    let xs: [f64; 7] = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0];
    for w in xs.windows(2) {
        let a = law.eta_effective(w[0]);
        let b = law.eta_effective(w[1]);
        assert!(a >= b, "η_eff not monotonically decreasing at ε̇={}", w[0]);
    }
}

#[test]
fn asymptotes_are_consistent() {
    // ε̇ → 0: η_newton bounded by ε̇_min^(1/n-1), saturated by the cap.
    let law = ViscosityLaw::default();
    let expected_floor = law.b_prefactor * (law.strain_rate_floor).powf(1.0 / law.n - 1.0);
    let cap = law.eta_max_cap;
    let expected_eff = expected_floor / (1.0 + (expected_floor / cap).powf(law.k_saturation))
        .powf(1.0 / law.k_saturation);
    let eta_at_zero = law.eta_effective(0.0);
    assert!(
        approx_rel(eta_at_zero, expected_eff, 1.0e-12),
        "η_eff(0) = {}, expected {}",
        eta_at_zero,
        expected_eff,
    );
    // ε̇ → ∞: η_eff → 0.
    assert!(law.eta_effective(1e12) < 1e-4);
}
