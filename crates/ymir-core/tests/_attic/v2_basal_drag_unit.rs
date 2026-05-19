//! Step 4 — basal-drag module unit integration tests.
//!
//! These cover the `basal_drag.rs` module through the public API.
//! Pure-algebraic invariants live in the module's internal `tests`
//! submodule; this file exercises the external-callable surface.

use ymir_core::tectonics_v2::basal_drag::{
    build_drag_diagonal_field, drag_diagonal_at_cell, BasalDragConfig, BasalDragLaw,
};
use ymir_core::tectonics_v2::field::Field2D;

#[test]
fn drag_diagonal_at_cell_equals_br_times_s_squared_at_default_exponent() {
    // Default exponent is 2.0 (decision D1). Verify the algebraic
    // form to 1e-14 over a spread of Br and S values.
    let cases: [(f64, f64); 5] = [
        (0.01, 0.2),
        (0.05, 1.0),
        (0.10, 1.5),
        (0.20, 0.7),
        (0.30, 1.2),
    ];
    for (br, s) in cases {
        let expected = br * s * s;
        let got = drag_diagonal_at_cell(br, s, 2.0);
        assert!(
            (got - expected).abs() < 1e-14,
            "br={br}, s={s}: got {got}, expected {expected}",
        );
    }
}

#[test]
fn drag_diagonal_is_zero_at_zero_s_irrespective_of_exponent() {
    // Oceanic annihilation guarantee for any positive exponent.
    for exp in [1.0_f64, 1.5, 2.0, 3.0] {
        let got = drag_diagonal_at_cell(0.3, 0.0, exp);
        assert_eq!(got, 0.0, "expected 0 at s=0 for exponent {exp}, got {got}");
    }
}

#[test]
fn disabled_returns_none_field() {
    let s = Field2D::filled(8, 8, 1.0);
    assert!(
        build_drag_diagonal_field(&BasalDragConfig::Disabled, &s).is_none(),
        "Disabled must produce None — the short-circuit is load-bearing for zero-cost",
    );
}

#[test]
fn enabled_returns_populated_field_matching_algebra() {
    let mut s = Field2D::new(5, 7);
    for j in 0..7 {
        for i in 0..5 {
            s.set(i, j, 0.5 + 0.15 * (i + j) as f64);
        }
    }
    let law = BasalDragLaw { br: 0.08, s_exponent: 2.0 };
    let cfg = BasalDragConfig::Enabled(law);
    let field = build_drag_diagonal_field(&cfg, &s).expect("Enabled must produce Some");
    assert_eq!(field.nx(), 5);
    assert_eq!(field.ny(), 7);
    for j in 0..7 {
        for i in 0..5 {
            let sij = s.get(i, j);
            let expected = 0.08 * sij * sij;
            let got = field.get(i, j);
            assert!(
                (got - expected).abs() < 1e-14,
                "cell ({i},{j}) got {got}, expected {expected}",
            );
        }
    }
}

#[test]
fn parse_round_trips_known_tokens() {
    assert!(matches!(
        BasalDragConfig::parse("enabled").unwrap(),
        BasalDragConfig::Enabled(_),
    ));
    assert!(matches!(
        BasalDragConfig::parse("on").unwrap(),
        BasalDragConfig::Enabled(_),
    ));
    assert!(matches!(
        BasalDragConfig::parse("disabled").unwrap(),
        BasalDragConfig::Disabled,
    ));
    assert!(matches!(
        BasalDragConfig::parse("off").unwrap(),
        BasalDragConfig::Disabled,
    ));
}

#[test]
fn parse_rejects_unknown_tokens_with_message() {
    let err = BasalDragConfig::parse("yes").unwrap_err();
    assert!(err.contains("yes"), "error should echo the bad token: {err}");
    assert!(err.contains("disabled|enabled"), "error should list valid tokens: {err}");
}

#[test]
fn label_is_stable_across_enabled_variants() {
    // `label()` must only reflect the variant, never the parameters —
    // this keeps report wiring robust to Br values changing.
    let e1 = BasalDragConfig::Enabled(BasalDragLaw { br: 0.01, s_exponent: 2.0 });
    let e2 = BasalDragConfig::Enabled(BasalDragLaw { br: 0.30, s_exponent: 2.0 });
    assert_eq!(e1.label(), e2.label());
    assert_eq!(e1.label(), "enabled");
    assert_eq!(BasalDragConfig::Disabled.label(), "disabled");
}
