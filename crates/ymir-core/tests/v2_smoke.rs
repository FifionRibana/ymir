//! Coupled thin-sheet + advection smoke test at a small grid.
//!
//! The full 300-step × 64² / 128² run is driven by the
//! `step_baseline` binary. This test guards against silent breakage
//! by running a shorter loop at a smaller grid — enough to detect
//! NaN/blow-up and confirm S mass stays machine-conserved.

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::{run_baseline, BaselineConfig};
use ymir_core::tectonics_v2::stokes::SheetConfig;

#[test]
fn coupled_loop_is_finite_and_mass_conserved() {
    let cfg = BaselineConfig {
        seed: 42,
        grid_nx: 32,
        grid_ny: 32,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: 50,
        cfl_factor: 0.3,
        forcing_amplitude: 0.1,
        sheet: SheetConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from("target/test_tmp_smoke"),
    };
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

    eprintln!(
        "smoke: wallclock = {:.3}s over {} steps; CG iters mean/max = {:.1}/{}",
        result.metrics.wallclock_total.as_secs_f64(),
        cfg.steps,
        result.metrics.cg_iter_mean,
        result.metrics.cg_iter_max,
    );
}
