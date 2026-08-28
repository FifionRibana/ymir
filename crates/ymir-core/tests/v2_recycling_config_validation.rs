//! Step 6 — RecyclingConfig validation: fractions must sum to 1
//! within 1e-9 absolute, fractions must be ≥ 0, delay ≥ 1.

use ymir_core::tectonics_v2::recycling::{RecyclingConfig, RecyclingConfigError};

#[test]
fn default_validates() {
    RecyclingConfig::default().validate().unwrap();
}

#[test]
fn fractions_not_summing_to_one_rejected() {
    let cfg = RecyclingConfig {
        arc_fraction: 0.5,
        coll_v_fraction: 0.0,
        rift_v_fraction: 0.0,
        spread_fraction: 0.0,
        mantle_loss_fraction: 0.0,
        mantle_delay_steps: 10,
    };
    assert!(matches!(cfg.validate(), Err(RecyclingConfigError::FractionsDoNotSumToOne { .. })));
}

#[test]
fn negative_fraction_rejected() {
    let cfg = RecyclingConfig {
        arc_fraction: -0.1,
        coll_v_fraction: 0.1,
        rift_v_fraction: 0.0,
        spread_fraction: 1.0,
        mantle_loss_fraction: 0.0,
        mantle_delay_steps: 10,
    };
    assert!(matches!(cfg.validate(), Err(RecyclingConfigError::NegativeFraction { .. })));
}

#[test]
fn f64_rounding_absorbed_by_1em9_tolerance() {
    // Construct a config whose fractions sum to 1 but with f64
    // representation that may give 0.9999999999999998 or similar.
    let cfg = RecyclingConfig {
        arc_fraction: 0.1 + 0.1 + 0.1, // 0.30000000000000004 or similar
        coll_v_fraction: 0.1,
        rift_v_fraction: 0.0,
        spread_fraction: 0.60,
        mantle_loss_fraction: 0.0,
        mantle_delay_steps: 10,
    };
    cfg.validate().unwrap();
}
