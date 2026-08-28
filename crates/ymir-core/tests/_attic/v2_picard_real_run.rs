//! Picard on a realistic scenario (preset `dynamic-accidented`,
//! short run). Not a numerical comparison with Newton — just a
//! guard against Picard quietly rotting on geometries that the MMS
//! parity test doesn't exercise.

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, NonlinearChoice, run_baseline,
};
use ymir_core::tectonics_v2::scales::Scales;

#[test]
fn picard_short_run_is_stable_on_dynamic_accidented() {
    let scales = Scales::default();
    let mut cfg = BaselineConfig::dynamic_accidented_defaults(&scales);
    cfg.grid_nx = 64;
    cfg.grid_ny = 64;
    cfg.steps = 15;
    cfg.nonlinear = NonlinearChoice::Picard;
    cfg.heightmap_fractions = Vec::new();
    cfg.output_dir = PathBuf::from("target/test_tmp_picard_real");

    let r = run_baseline(&cfg);
    assert!(
        r.metrics.mass_drift_relative.abs() < 1.0e-10,
        "mass drift = {:.3e}",
        r.metrics.mass_drift_relative,
    );
    assert!(r.metrics.vmax_peak.is_finite());
    assert!(r.metrics.mass_s_final.is_finite());

    let newton = r.metrics.newton.as_ref().expect("Newton aggregate present");
    // No Diverged outcomes tolerated.
    assert_eq!(newton.diverged, 0, "Picard diverged on a real run: {:?}", newton);

    eprintln!(
        "Picard real-run: wallclock {:.3}s; outer iters mean/max {:.1}/{}; outcomes {:?} conv/{} stall/{} div/{} cap",
        r.metrics.wallclock_total.as_secs_f64(),
        newton.outer_iters_mean(),
        newton.outer_iters_max(),
        newton.converged,
        newton.stalled,
        newton.diverged,
        newton.capped,
    );
}
