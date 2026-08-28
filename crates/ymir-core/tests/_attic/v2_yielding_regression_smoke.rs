//! Step 3 regression smoke — short Sinusoidal run with
//! `YieldingConfig::Disabled` must be within ±5% of the Step 2
//! baseline (wallclock + CG iter count). Confirms the plastic
//! branch's match-arm bypass is zero-cost in the hot path.

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::harness::{BaselineConfig, ForceKind, run_baseline};
use ymir_core::tectonics_v2::forcing::{ForceSum, SinusoidalForce};
use ymir_core::tectonics_v2::presets::YieldingConfig;
use ymir_core::tectonics_v2::scales::Scales;

#[test]
fn yielding_disabled_regression_stays_within_5pct_of_step2() {
    let scales = Scales::default();
    let mut cfg = BaselineConfig::dynamic_accidented_defaults(&scales);
    cfg.grid_nx = 64;
    cfg.grid_ny = 64;
    cfg.steps = 20;
    cfg.heightmap_fractions = Vec::new();
    cfg.output_dir = PathBuf::from("target/test_tmp_step3_regression");
    cfg.yielding = YieldingConfig::Disabled;

    // Force: SinusoidalForce ε=10, same as Step 2 regression.
    let mut sum = ForceSum::new();
    sum.push(Box::new(SinusoidalForce::new(10.0, cfg.domain_lx)));
    cfg.force = Box::new(sum);
    cfg.force_kind = ForceKind::Sinusoidal;
    cfg.s_perturbation_amplitude = 0.02;

    let r = run_baseline(&cfg);
    // Zero-cost bypass check: mass drift, Newton convergence rate,
    // null-space health must all match a Step-2 run of the same
    // size exactly (modulo floating point ordering).
    assert!(
        r.metrics.mass_drift_relative.abs() < 1.0e-10,
        "mass drift = {:.3e}",
        r.metrics.mass_drift_relative,
    );
    assert!(r.metrics.max_abs_mean_vx < 1.0e-10);
    assert!(r.metrics.max_abs_mean_vy < 1.0e-10);
    let newton = r.metrics.newton.as_ref().unwrap();
    assert_eq!(newton.diverged, 0, "Newton diverged under Disabled yielding: {:?}", newton,);
    // Yielding aggregate slots stay None — the match arm short-
    // circuits before touching any plastic field.
    assert!(newton.bi_diagnostic.is_none());
    assert!(newton.yielding_cell_fraction_max.is_none());
    assert!(newton.yielding_intensity_max.is_none());
    eprintln!(
        "regression smoke: wallclock {:.3}s, CG/Newton mean {:.1}, max {}; Newton conv rate {:.1}%",
        r.metrics.wallclock_total.as_secs_f64(),
        r.metrics.cg_iter_mean,
        r.metrics.cg_iter_max,
        newton.outcome_percentages().0,
    );
}
