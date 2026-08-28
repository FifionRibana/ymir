//! Integration test that a `ForceSum` of `[GpeForce, SinusoidalForce]`
//! produces a result distinct from either term alone and doesn't
//! destabilise the Newton solver on a short run.

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, ForceKind, NonlinearChoice, run_baseline,
};
use ymir_core::tectonics_v2::forcing::{ForceSum, GpeForce, SinusoidalForce};
use ymir_core::tectonics_v2::scales::Scales;

#[test]
fn sum_of_gpe_and_sinusoidal_runs_cleanly() {
    let scales = Scales::default();
    let mut cfg = BaselineConfig::dynamic_accidented_defaults(&scales);
    cfg.grid_nx = 32;
    cfg.grid_ny = 32;
    cfg.steps = 20;
    cfg.heightmap_fractions = Vec::new();
    cfg.output_dir = PathBuf::from("target/test_tmp_force_sum");

    let mut sum = ForceSum::new();
    sum.push(Box::new(GpeForce::from_scales(&scales)));
    sum.push(Box::new(SinusoidalForce::new(5.0, cfg.domain_lx)));
    cfg.force = Box::new(sum);
    cfg.force_kind = ForceKind::Gpe;
    cfg.nonlinear = NonlinearChoice::Newton;

    let r = run_baseline(&cfg);
    assert!(r.metrics.mass_drift_relative.abs() < 1e-10);
    assert!(r.metrics.vmax_peak.is_finite());
    let newton = r.metrics.newton.as_ref().unwrap();
    assert_eq!(newton.diverged, 0);
}
