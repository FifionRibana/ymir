//! Smooth replacements for the non-differentiable operations used in the
//! viscosity pipeline. Newton with a finite-difference Jacobian oscillates
//! around corners of `clamp` and `min`; these helpers keep the asymptotic
//! behaviour intact while making the functions globally C¹.

/// Smooth saturation: `f(x) ≈ x` for `x ≪ x_max`, `f(x) → x_max` as
/// `x → ∞`. Strictly bounded by `x_max` for any positive input. C¹ globally.
///
/// With the sharpness `k = 4`:
/// * within 1% of `x` for `x ≤ 0.1 · x_max`
/// * `f(x_max) ≈ 0.84 · x_max`
/// * within 2% of `x_max` for `x ≥ 2 · x_max`
#[inline]
pub fn smooth_saturate(x: f64, x_max: f64) -> f64 {
    const K: f64 = 4.0;
    let ratio = x / x_max;
    x / (1.0 + ratio.powf(K)).powf(1.0 / K)
}

/// Smooth minimum of two positive values: approaches `min(a, b)` as
/// `p → ∞`; smaller `p` widens the transition band. `p = 4` gives a
/// transition width of about 30% around the equality point.
///
/// Guarded against zero or negative inputs via a small floor so
/// `0^(-p)` never produces NaN/Inf in the caller.
#[inline]
pub fn soft_min_harmonic(a: f64, b: f64, p: f64) -> f64 {
    const FLOOR: f64 = 1e-30;
    let inv_a = a.max(FLOOR).powf(-p);
    let inv_b = b.max(FLOOR).powf(-p);
    (inv_a + inv_b).powf(-1.0 / p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smooth_saturate_preserves_small_values() {
        let x_max = 1e4;
        for &x in &[1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0] {
            let saturated = smooth_saturate(x, x_max);
            let relative_error = (saturated - x).abs() / x;
            assert!(
                relative_error < 0.01,
                "smooth_saturate({x}, {x_max}) = {saturated} differs from {x} by more than 1%"
            );
        }
    }

    #[test]
    fn smooth_saturate_bounded_above() {
        // Mathematically strictly less than x_max; floating-point rounding
        // can land on x_max exactly once the ratio dominates, so we only
        // require <= x_max.
        let x_max = 1e4;
        for &x in &[1e3, 1e4, 1e5, 1e6, 1e10] {
            let saturated = smooth_saturate(x, x_max);
            assert!(
                saturated <= x_max,
                "smooth_saturate({x}, {x_max}) = {saturated} exceeds x_max"
            );
        }
    }

    #[test]
    fn smooth_saturate_monotonic() {
        let x_max = 1e4;
        let xs: Vec<f64> = (0..100).map(|i| 10.0_f64.powf(i as f64 / 10.0 - 5.0)).collect();
        for w in xs.windows(2) {
            let f1 = smooth_saturate(w[0], x_max);
            let f2 = smooth_saturate(w[1], x_max);
            assert!(f1 <= f2, "smooth_saturate not monotonic at x={}", w[0]);
        }
    }

    #[test]
    fn soft_min_harmonic_returns_approximately_min() {
        let p = 4.0;
        assert!((soft_min_harmonic(1.0, 100.0, p) - 1.0).abs() < 0.01);
        assert!((soft_min_harmonic(100.0, 1.0, p) - 1.0).abs() < 0.01);
        assert!((soft_min_harmonic(0.001, 1.0, p) - 0.001).abs() / 0.001 < 0.01);
    }

    #[test]
    fn soft_min_harmonic_smooth_at_equality() {
        let p = 4.0;
        let result = soft_min_harmonic(10.0, 10.0, p);
        assert!(
            result > 8.0 && result < 10.0,
            "soft_min_harmonic(10, 10, 4) = {result}; expected strictly between 8 and 10"
        );
    }

    #[test]
    fn soft_min_harmonic_handles_zero_input() {
        let p = 4.0;
        // With the floor guard, a zero input should not produce NaN/Inf.
        let result = soft_min_harmonic(0.0, 1.0, p);
        assert!(result.is_finite(), "soft_min_harmonic(0, 1, 4) = {result}");
    }
}
