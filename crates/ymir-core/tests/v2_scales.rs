//! Integration-level coverage for `Scales` dim ↔ nondim round-trips.
//!
//! The unit tests in `tectonics_v2::scales` cover the same contract
//! with hard-coded numbers; this file exercises the full public
//! surface through the crate's external API (which is all a
//! downstream crate would see).

use ymir_core::tectonics_v2::scales::{Scales, GRAVITY, SECONDS_PER_MYR};

fn approx_rel(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
}

#[test]
fn default_scales_match_design_note() {
    let s = Scales::default();
    // η* = ρ*·g·τ*·S* per the scales derivation. Expected ≈ 1.07e24 Pa·s.
    let expected_eta = 3300.0 * GRAVITY * (30.0 * SECONDS_PER_MYR) * 35.0e3;
    assert!(approx_rel(s.viscosity, expected_eta, 1e-12));
    assert!(approx_rel(s.viscosity, 1.07e24, 5e-2));
}

#[test]
fn roundtrip_covers_every_quantity() {
    let s = Scales::default();
    for x in [1.0, 1.23e4, 9.87e7, 3.14e-3] {
        assert!(approx_rel(s.to_dim_length(s.to_nondim_length(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_thickness(s.to_nondim_thickness(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_time(s.to_nondim_time(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_velocity(s.to_nondim_velocity(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_viscosity(s.to_nondim_viscosity(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_stress(s.to_nondim_stress(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_pressure(s.to_nondim_pressure(x)), x, 1e-12));
        assert!(approx_rel(s.to_dim_body_force(s.to_nondim_body_force(x)), x, 1e-12));
    }
}

#[test]
fn alternate_primary_scales_produce_consistent_derived() {
    let s = Scales::from_primary(500.0e3, 50.0e3, 20.0 * SECONDS_PER_MYR, 3000.0);
    assert!(approx_rel(s.velocity, s.length / s.time, 1e-14));
    assert!(approx_rel(s.viscosity, s.density * GRAVITY * s.time * s.thickness, 1e-14));
    assert!(approx_rel(s.stress, s.viscosity * s.strain_rate, 1e-14));
    assert!(s.argand > 0.0);
}
