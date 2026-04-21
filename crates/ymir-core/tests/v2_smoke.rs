//! Coupled thin-sheet + advection smoke test at a small grid.
//!
//! The full 300-step × 64² / 128² run is driven by the
//! `step_baseline` binary. This test guards against silent breakage
//! by running a shorter loop at a smaller grid — enough to detect
//! NaN/blow-up and confirm S mass stays machine-conserved. Step 1
//! uses the Newton nonlinear solver with the default
//! `dynamic-accidented` preset.

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::harness::{run_baseline, BaselineConfig};

#[test]
fn coupled_loop_is_finite_and_mass_conserved() {
    let mut cfg = BaselineConfig::dynamic_accidented_defaults();
    cfg.grid_nx = 32;
    cfg.grid_ny = 32;
    cfg.steps = 50;
    cfg.heightmap_fractions = Vec::new();
    cfg.output_dir = PathBuf::from("target/test_tmp_smoke");
    let result = run_baseline(&cfg);

    assert!(
        result.metrics.mass_drift_relative.abs() < 1e-10,
        "mass drift = {:.3e}",
        result.metrics.mass_drift_relative,
    );
    assert!(result.metrics.wallclock_total.as_secs_f64().is_finite());
    assert!(result.metrics.vmax_peak.is_finite());
    assert!(result.metrics.mass_s_final.is_finite());

    assert!(result.metrics.max_abs_mean_vx < 1e-10);
    assert!(result.metrics.max_abs_mean_vy < 1e-10);

    // Step 1 additions: Newton must converge on every step.
    let newton = result.metrics.newton.as_ref().expect("Step 1 Newton aggregate present");
    assert_eq!(
        newton.diverged, 0,
        "Newton diverged on the smoke run: {:?}",
        newton,
    );
    assert!(newton.converged > 0, "no Newton solve recorded as converged");

    eprintln!(
        "smoke: wallclock = {:.3}s over {} steps; Newton outer mean/max = {:.1}/{}",
        result.metrics.wallclock_total.as_secs_f64(),
        cfg.steps,
        newton.outer_iters_mean(),
        newton.outer_iters_max(),
    );
}
