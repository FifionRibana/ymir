//! Issue #129 Phase 2 Track A Stage V — Stein-Stein 1992 oceanic
//! bathymetry closure validation tests.
//!
//! Four tests under default features, all on the public API surface
//! of [`ymir_core::tectonics_c1::closures::oceanic_bathymetry`]:
//!
//! 1. [`stein_stein_reproduces_5_age_points`] — **KEY quantitative
//!    anchor** against five published reference points from Stein
//!    & Stein 1992 (*Nature* 359:123-129). Tolerance ±50 m on the
//!    depth-in-meters formula. This is the load-bearing validation
//!    that distinguishes Phase 2 Track A from Phase 1.4 (Lague
//!    2014's "no universal K" framework) — Stein-Stein provides
//!    quantitative anchor points, so the closure can be tested
//!    against the paper instead of "visual review only".
//! 2. [`ridge_axis_shallower_than_asymptote`] — sanity ordering on
//!    `stein_stein_depth(0) < stein_stein_depth(1000)`. Direct
//!    consequence of the formula structure.
//! 3. [`disabled_no_op_on_altitude`] — `enabled = false` →
//!    altitude bit-identical to initial state. W4 closure-
//!    isolation discipline at the integration boundary.
//! 4. [`continental_cells_unmodified`] — plate-type filter
//!    correctness. Constructs a mixed continental + oceanic
//!    plate-type field, runs the apply function, asserts
//!    continental cells are byte-identical pre/post and at least
//!    one oceanic cell was modified.
//!
//! Pairs with Phase 2 Track A Stage A (`c1_phase_2_bathymetry_*`)
//! which adds acceptance tests with run-boundary observability of
//! Architecture C's transient altitude imprint.

use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::source_term::{
    apply_stein_stein_bathymetry, stein_stein_depth,
};
use ymir_core::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use ymir_core::tectonics_v2::field::Field2D;

/// Construct a mixed-plate-type fixture: top half continental,
/// bottom half oceanic; uniform initial altitude `0.5`; uniform
/// `age = 50.0` (`age_ma = 50 · 0.667 ≈ 33`, well into the S-S
/// old regime → depth ≈ 5000 m → altitude_nondim ≈ −1.0).
fn split_continental_oceanic(nx: usize, ny: usize) -> (GridF32, Field2D, PlateTypeField) {
    let altitude = GridF32::new(nx, ny, 0.5);
    let mut age = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            age.set(i, j, 50.0);
        }
    }
    let mut plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
    for j in 0..ny / 2 {
        for i in 0..nx {
            plate_type.set(i, j, PlateType::Continental);
        }
    }
    (altitude, age, plate_type)
}

/// **Test 1 — quantitative anchor (KEY).** Reproduce the 5
/// published Stein-Stein 1992 depth-age points to within ±50 m
/// tolerance. Stein & Stein 1992 §2 / Table 1 publishes ridge
/// depth `d_r = 2600 m`, young-regime `b = 365 m / √Ma`,
/// asymptotic depth `d_∞ = 5651 m`, exponential time constant
/// `α = 0.0278 Ma⁻¹`, continuity coefficient `C = 2473 m`, regime
/// crossover `t_c = 20 Ma`.
///
/// Reference depths computed from those parameters at five
/// canonical ages — load-bearing test for Phase 2 Track A's
/// paper-faithfulness claim.
#[test]
fn stein_stein_reproduces_5_age_points() {
    let params = SteinSteinParams::default();

    // 5 reference points spanning young √t regime, old exp
    // regime, and asymptotic saturation.
    let cases: [(f64, f64); 5] = [
        (0.0, 2600.0),   // Ridge axis. Young regime, t = 0.
        (10.0, 3754.0),  // Young √t regime: 2600 + 365·√10 ≈ 3754.23.
        (50.0, 5035.0),  // Old exp regime: 5651 - 2473·exp(-1.39) ≈ 5034.84.
        (100.0, 5498.0), // Old exp regime: 5651 - 2473·exp(-2.78) ≈ 5497.73.
        (150.0, 5613.0), // Old exp regime, near asymptote: 5651 - 2473·exp(-4.17) ≈ 5612.82.
    ];

    eprintln!("c1_phase_2 Stage V Test 1 — Stein-Stein 1992 5-point quantitative anchor:");
    eprintln!(
        "  Reference: Stein & Stein (1992), Nature 359:123-129, Table 1 (GDH1 plate model)"
    );
    eprintln!("  Tolerance: ±50.0 m");
    eprintln!();
    eprintln!("    {:>8} | {:>12} | {:>12} | {:>10}", "Age (Ma)", "S-S pub (m)", "Computed (m)", "Error (m)");
    eprintln!("    {:->8}-+-{:->12}-+-{:->12}-+-{:->10}", "", "", "", "");

    let mut max_error_m = 0.0_f64;
    for (age_ma, expected_depth_m) in cases.iter() {
        let computed = stein_stein_depth(*age_ma, &params);
        let error = (computed - expected_depth_m).abs();
        if error > max_error_m {
            max_error_m = error;
        }
        eprintln!(
            "    {:>8.1} | {:>12.1} | {:>12.3} | {:>10.3}",
            age_ma, expected_depth_m, computed, error
        );

        assert!(
            error < 50.0,
            "Stein-Stein 1992 reproduction at age {} Ma: computed depth {:.3} m differs from \
             published reference {:.1} m by {:.3} m (>50 m tolerance). Check `SteinSteinParams` \
             defaults vs paper Table 1.",
            age_ma,
            computed,
            expected_depth_m,
            error,
        );
    }

    eprintln!();
    eprintln!("  Max error across 5 anchor points: {max_error_m:.3} m");
    eprintln!("  Phase 2 Track A quantitative anchor: PASS (paper-faithful within ±50 m)");
}

/// **Test 2 — sanity ordering.** Ridge-axis (`t = 0`) depth is
/// shallower than late-asymptote (`t = 1000 Ma`) depth. Locks the
/// monotonicity of S-S subsidence — formula structure invariant
/// that must hold regardless of parameter tuning.
#[test]
fn ridge_axis_shallower_than_asymptote() {
    let params = SteinSteinParams::default();

    let ridge_depth = stein_stein_depth(0.0, &params);
    // 1000 Ma is far past saturation:
    //   exp(-0.0278 · 1000) ≈ 8 × 10⁻¹³ → the exp term contributes
    //   < 10⁻⁹ m. depth ≈ d_∞.
    let asymptote_depth = stein_stein_depth(1000.0, &params);

    eprintln!("c1_phase_2 Stage V Test 2 — ridge_axis_shallower_than_asymptote:");
    eprintln!("  ridge_depth (t = 0 Ma)        = {:.3} m (params.ridge_depth_m = {})", ridge_depth, params.ridge_depth_m);
    eprintln!(
        "  asymptote_depth (t = 1000 Ma) = {:.3} m (params.asymptotic_depth_m = {})",
        asymptote_depth, params.asymptotic_depth_m
    );

    assert!(
        ridge_depth < asymptote_depth,
        "S-S monotonicity broken: ridge depth {:.3} m ≥ asymptote depth {:.3} m. \
         The depth-age relation must be monotonically increasing per Stein-Stein 1992 §2.",
        ridge_depth,
        asymptote_depth,
    );
    assert!(
        (ridge_depth - params.ridge_depth_m).abs() < 1.0,
        "Ridge axis (t=0) depth should equal `ridge_depth_m = {}`; got {:.6} m \
         (residual {:.6} m > 1.0 m)",
        params.ridge_depth_m,
        ridge_depth,
        (ridge_depth - params.ridge_depth_m).abs(),
    );
    assert!(
        (asymptote_depth - params.asymptotic_depth_m).abs() < 1.0,
        "Late-asymptote (t=1000 Ma) depth should equal `asymptotic_depth_m = {}`; \
         got {:.6} m (residual {:.6} m > 1.0 m)",
        params.asymptotic_depth_m,
        asymptote_depth,
        (asymptote_depth - params.asymptotic_depth_m).abs(),
    );
}

/// **Test 3 — `enabled = false` integration no-op.** Mirrors the
/// W4 closure-isolation discipline applied to Phase 1.2 Davis-
/// Suppe, Phase 1.3 equilibrium-height, Phase 1.4 erosion. With
/// `enabled = false`, [`apply_stein_stein_bathymetry`] must early-
/// return; the `altitude` field is byte-identical to its
/// pre-call state regardless of `age`, `plate_type`, or any other
/// `params` field.
#[test]
fn disabled_no_op_on_altitude() {
    let (mut altitude, age, plate_type) = split_continental_oceanic(8, 8);
    let initial = altitude.data.clone();
    let params = SteinSteinParams {
        enabled: false,
        // Mix in non-default tunables to ensure the early-return
        // gate is the only thing protecting the no-op (not some
        // identity coincidence on defaults).
        ridge_depth_m: 9999.0,
        depth_scale_m: 1.0,
        ..SteinSteinParams::default()
    };

    apply_stein_stein_bathymetry(&mut altitude, &age, &plate_type, &params);

    for k in 0..initial.len() {
        assert_eq!(
            altitude.data[k], initial[k],
            "`enabled = false` must not touch any cell; mismatch at flat index {k}: \
             before = {}, after = {}",
            initial[k], altitude.data[k],
        );
    }

    eprintln!("c1_phase_2 Stage V Test 3 — disabled_no_op_on_altitude:");
    eprintln!(
        "  64 cells preserved bit-identical despite pathological params (ridge=9999, scale=1.0)"
    );
}

/// **Test 4 — `PlateType::Oceanic` filter correctness.** With a
/// mixed continental + oceanic plate-type field, apply the
/// closure with `enabled = true` and assert:
///
/// - Continental cells (top half) are byte-identical pre/post.
/// - At least one oceanic cell (bottom half) was modified.
///
/// This locks the plate-type discrimination that Stein-Stein's
/// "oceanic lithosphere subsides with age" semantics demand.
#[test]
fn continental_cells_unmodified() {
    let (mut altitude, age, plate_type) = split_continental_oceanic(8, 8);
    let initial = altitude.data.clone();
    let params = SteinSteinParams::default();

    apply_stein_stein_bathymetry(&mut altitude, &age, &plate_type, &params);

    // Continental top half (j < 4): byte-identical.
    let mut continental_count = 0_usize;
    let mut continental_changed = 0_usize;
    for j in 0..4 {
        for i in 0..8 {
            let idx = j * 8 + i;
            continental_count += 1;
            if altitude.data[idx] != initial[idx] {
                continental_changed += 1;
            }
        }
    }

    // Oceanic bottom half (j ≥ 4): at least one changed.
    let mut oceanic_count = 0_usize;
    let mut oceanic_changed = 0_usize;
    let mut sample_oceanic_altitude = f32::NAN;
    for j in 4..8 {
        for i in 0..8 {
            let idx = j * 8 + i;
            oceanic_count += 1;
            if (altitude.data[idx] - initial[idx]).abs() > 1e-9 {
                oceanic_changed += 1;
                if sample_oceanic_altitude.is_nan() {
                    sample_oceanic_altitude = altitude.data[idx];
                }
            }
        }
    }

    eprintln!("c1_phase_2 Stage V Test 4 — continental_cells_unmodified:");
    eprintln!(
        "  Continental cells: {continental_changed} / {continental_count} modified (expected 0)"
    );
    eprintln!(
        "  Oceanic cells:     {oceanic_changed} / {oceanic_count} modified (expected > 0)"
    );
    eprintln!(
        "  Sample oceanic altitude post-S-S: {sample_oceanic_altitude:.4} \
         (age = 50 → age_ma ≈ 33 → depth ≈ 5000 m → altitude ≈ −1.0)"
    );

    assert_eq!(
        continental_changed, 0,
        "S-S leaked into {continental_changed} continental cells; the `PlateType::Oceanic` \
         filter is wrong (expected 0 continental modifications, got {continental_changed})",
    );
    assert!(
        oceanic_changed > 0,
        "S-S apply did not touch any oceanic cell (0 / {oceanic_count} modified); the apply \
         function may have been silently no-op even with `enabled = true`",
    );
}
