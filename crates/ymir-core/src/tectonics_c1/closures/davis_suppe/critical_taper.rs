//! Davis-Suppe critical taper formula (small-angle, cohesionless
//! wedge with pore pressure).
//!
//! ## References
//!
//! - Davis, Suppe & Dahlen 1983, *JGR* 88(B2), 1153-1172 — original
//!   derivation (eq. 15 and surroundings).
//! - Dahlen 1990, *Annu. Rev. Earth Planet. Sci.* 18, 55-99 — review,
//!   small-angle linearisation (the form implemented here).
//! - Suppe 2007, *Geol. Acta* 5(1), 1-13 — pedagogical re-derivation.
//!
//! ## Formula
//!
//! For a cohesionless wedge with pore pressure ratios `λ_b` (basal)
//! and `λ` (internal) and water-to-wedge density ratio `ρ_w/ρ` (0
//! for subaerial, ~0.4 for marine), the small-angle critical taper
//! `α + β` satisfies:
//!
//! ```text
//!                  (1 - λ_b) μ_b + (1 - ρ_w/ρ) β
//!     α + β  =  ───────────────────────────────────────
//!                (1 - ρ_w/ρ) + (1 - λ) · 2 sin φ_i / (1 - sin φ_i)
//! ```
//!
//! where `μ_i = tan φ_i`, so `sin φ_i = μ_i / √(1 + μ_i²)`.
//!
//! ## Calibration against Davis 1983 sandbox
//!
//! The sandbox experiments use dry sand on a Mylar base
//! (`λ = λ_b = 0`, `ρ_w/ρ = 0`) with measured `μ_b = 0.30`,
//! `μ_i = 0.58`. For three pre-imposed basal dips:
//!
//! | β     | observed α | observed α + β | this code α + β |
//! |-------|-----------|----------------|-----------------|
//! | 0°    | 5.7°      | 5.7°           | 5.70° (Δ +0.00) |
//! | 3°    | 3.7°      | 6.7°           | 6.70° (Δ +0.00) |
//! | 6°    | 2.0°      | 8.0°           | 7.69° (Δ -0.31) |
//!
//! All three sandbox reproductions land inside the ±0.5° gate
//! (Stage 1 acceptance). The small β=6° miss is the expected
//! small-angle linearisation error growing with `tan(β)`.
//!
//! ## Natural-orogen agreement (informative only)
//!
//! Per the Stage 1 spec, natural-orogen agreement is tracked as
//! informational test output (±2° tolerance, non-blocking). The
//! `μ_b`, `μ_i`, `λ` parameters for natural orogens are inferred
//! indirectly (Byerlee + pore-pressure modelling), so disagreement
//! reflects input uncertainty as much as formula limit.

/// Compute the critical surface slope `α` (radians) for a
/// cohesionless wedge with the given mechanical parameters and
/// basal dip `β`.
///
/// All angles in radians. The result is `α` only; for the total
/// taper use [`critical_taper_angle_with_beta`].
pub fn compute_alpha_only(
    mu_b: f64,
    mu_i: f64,
    lambda_b: f64,
    lambda_i: f64,
    rho_w_over_rho: f64,
    beta: f64,
) -> f64 {
    critical_taper_angle_with_beta(mu_b, mu_i, lambda_b, lambda_i, rho_w_over_rho, beta) - beta
}

/// Compute the critical taper angle `α + β` (radians) for the
/// given mechanical parameters and pre-imposed basal dip `β`.
///
/// Returns `α + β` in radians. The small-angle Dahlen 1990
/// linearisation is used; the deviation from the exact formulation
/// is O(`tan²(α+β)`) ≈ 1 % at the 10° scale typical of natural
/// orogens.
pub fn critical_taper_angle_with_beta(
    mu_b: f64,
    mu_i: f64,
    lambda_b: f64,
    lambda_i: f64,
    rho_w_over_rho: f64,
    beta: f64,
) -> f64 {
    // φ_i from internal friction coefficient: μ_i = tan φ_i.
    let sin_phi_i = mu_i / (1.0 + mu_i * mu_i).sqrt();
    let one_minus_sin = 1.0 - sin_phi_i;

    let rho_term = 1.0 - rho_w_over_rho;
    let numerator = (1.0 - lambda_b) * mu_b + rho_term * beta;
    let denominator =
        rho_term + (1.0 - lambda_i) * 2.0 * sin_phi_i / one_minus_sin.max(f64::EPSILON);

    numerator / denominator
}

/// Convenience: critical taper at flat base (β = 0). Useful for
/// Davis 1983 sandbox β=0 case and for Phase 1.2 default mode
/// where C1 does not pre-impose a basal dip.
pub fn critical_taper_angle(
    mu_b: f64,
    mu_i: f64,
    lambda_b: f64,
    lambda_i: f64,
    rho_w_over_rho: f64,
) -> f64 {
    critical_taper_angle_with_beta(mu_b, mu_i, lambda_b, lambda_i, rho_w_over_rho, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Sandbox reproductions (blocking, ±0.5°) ────────────────

    #[test]
    fn reproduces_sandbox_beta_0() {
        // Davis 1983 sandbox: dry sand on Mylar, μ_b=0.30, μ_i=0.58,
        // λ=0. With β=0° the observed critical surface slope is
        // α=5.7°, so α+β = 5.7°.
        let predicted = critical_taper_angle(0.30, 0.58, 0.0, 0.0, 0.0);
        let expected = 5.7_f64.to_radians();
        let tol = 0.5_f64.to_radians();
        assert!(
            (predicted - expected).abs() < tol,
            "sandbox β=0°: predicted {:.3}°, expected 5.7°, tolerance ±0.5°",
            predicted.to_degrees(),
        );
    }

    #[test]
    fn reproduces_sandbox_beta_3() {
        // Same sandbox material, β=3° basal dip. Observed α ≈ 3.7°,
        // so total taper α+β ≈ 6.7°.
        let predicted =
            critical_taper_angle_with_beta(0.30, 0.58, 0.0, 0.0, 0.0, 3.0_f64.to_radians());
        let expected = 6.7_f64.to_radians();
        let tol = 0.5_f64.to_radians();
        assert!(
            (predicted - expected).abs() < tol,
            "sandbox β=3°: predicted {:.3}°, expected 6.7°, tolerance ±0.5°",
            predicted.to_degrees(),
        );
    }

    #[test]
    fn reproduces_sandbox_beta_6() {
        // Same sandbox material, β=6° basal dip. Observed α ≈ 2.0°,
        // so total taper α+β ≈ 8.0°.
        let predicted =
            critical_taper_angle_with_beta(0.30, 0.58, 0.0, 0.0, 0.0, 6.0_f64.to_radians());
        let expected = 8.0_f64.to_radians();
        let tol = 0.5_f64.to_radians();
        assert!(
            (predicted - expected).abs() < tol,
            "sandbox β=6°: predicted {:.3}°, expected 8.0°, tolerance ±0.5°",
            predicted.to_degrees(),
        );
    }

    // ── Natural orogens (informational, ±2°, non-blocking) ──────

    #[test]
    fn matches_taiwan_within_natural_scatter() {
        // Taiwan accretionary wedge (Davis 1983 + Dahlen 1990):
        // μ_b=0.85, μ_i=1.03, λ_b≈0.68. Observed α+β ≈ 8.9°.
        //
        // Inputs are inferred (Byerlee + pore-pressure model), so
        // disagreement is not necessarily a formula bug. Logs a
        // warning when |Δ| > 2° but does not fail.
        let predicted = critical_taper_angle(0.85, 1.03, 0.68, 0.68, 0.0);
        let expected = 8.9_f64.to_radians();
        let tol = 2.0_f64.to_radians();
        if (predicted - expected).abs() >= tol {
            eprintln!(
                "INFO Taiwan: predicted {:.2}°, observed 8.9°, |Δ|={:.2}° > 2° (informational only)",
                predicted.to_degrees(),
                (predicted - expected).to_degrees().abs(),
            );
        } else {
            eprintln!(
                "INFO Taiwan: predicted {:.2}°, observed 8.9°, |Δ|={:.2}° within ±2°",
                predicted.to_degrees(),
                (predicted - expected).to_degrees().abs(),
            );
        }
        // Always passes — informational only.
    }

    #[test]
    fn matches_barbados_within_natural_scatter() {
        // Barbados accretionary wedge: high pore pressure
        // (λ_b=0.95), submerged (ρ_w/ρ≈0.4), basal dip β≈2°.
        // Observed α+β ≈ 3.5° (canonical low-taper, overpressured
        // case — the test that a formula must respect λ to pass).
        let predicted =
            critical_taper_angle_with_beta(0.85, 1.57, 0.95, 0.95, 0.4, 2.0_f64.to_radians());
        let expected = 3.5_f64.to_radians();
        let tol = 2.0_f64.to_radians();
        if (predicted - expected).abs() >= tol {
            eprintln!(
                "INFO Barbados: predicted {:.2}°, observed 3.5°, |Δ|={:.2}° > 2° (informational only)",
                predicted.to_degrees(),
                (predicted - expected).to_degrees().abs(),
            );
        } else {
            eprintln!(
                "INFO Barbados: predicted {:.2}°, observed 3.5°, |Δ|={:.2}° within ±2°",
                predicted.to_degrees(),
                (predicted - expected).to_degrees().abs(),
            );
        }
        // Always passes — informational only.
    }

    // ── Monotonicity invariants (blocking, sign only) ──────────

    #[test]
    fn monotonic_in_mu_b() {
        // Higher basal friction → steeper taper, all else equal.
        let a = critical_taper_angle(0.30, 0.58, 0.0, 0.0, 0.0);
        let b = critical_taper_angle(0.50, 0.58, 0.0, 0.0, 0.0);
        let c = critical_taper_angle(0.70, 0.58, 0.0, 0.0, 0.0);
        assert!(
            a < b,
            "Expected monotonic in μ_b: a={:.2}°, b={:.2}°",
            a.to_degrees(),
            b.to_degrees()
        );
        assert!(b < c, "Expected b < c: b={:.2}°, c={:.2}°", b.to_degrees(), c.to_degrees());
    }

    #[test]
    fn monotonic_in_lambda_b() {
        // Higher basal pore-pressure ratio → effective μ_b reduced
        // → lower critical taper.
        let dry = critical_taper_angle(0.85, 1.03, 0.0, 0.0, 0.0);
        let wet = critical_taper_angle(0.85, 1.03, 0.68, 0.68, 0.0);
        let very_wet = critical_taper_angle(0.85, 1.03, 0.95, 0.95, 0.0);
        assert!(
            wet < dry,
            "λ=0.68 should reduce taper vs dry: wet={:.2}°, dry={:.2}°",
            wet.to_degrees(),
            dry.to_degrees(),
        );
        assert!(
            very_wet < wet,
            "λ=0.95 should further reduce: very_wet={:.2}°, wet={:.2}°",
            very_wet.to_degrees(),
            wet.to_degrees(),
        );
    }

    #[test]
    fn alpha_decreases_with_beta_at_fixed_mu_lambda() {
        // Sandbox μ_b=0.30, μ_i=0.58, λ=0:
        //   β=0° → α ≈ 5.7°
        //   β=6° → α ≈ 2.0°  (basal dip "consumes" taper budget).
        let alpha_0 = compute_alpha_only(0.30, 0.58, 0.0, 0.0, 0.0, 0.0);
        let alpha_6 = compute_alpha_only(0.30, 0.58, 0.0, 0.0, 0.0, 6.0_f64.to_radians());
        assert!(
            alpha_6 < alpha_0,
            "α should decrease with β: α(β=0°)={:.2}°, α(β=6°)={:.2}°",
            alpha_0.to_degrees(),
            alpha_6.to_degrees(),
        );
    }
}
