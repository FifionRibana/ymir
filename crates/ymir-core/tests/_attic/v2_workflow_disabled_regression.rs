//! Step 12 Phase 1 regression — `WorkflowConfig::Disabled` is
//! bit-identical to a direct `run_baseline` call (Step 11 standalone
//! contract).
//!
//! Acceptance criterion #15 of `step12_issue.md`:
//!
//! > Default state with `WorkflowConfig::Disabled` is equivalent to
//! > running tectonic only (no Phase A loop, no Phase B). Identical
//! > to running Step 11 directly.
//!
//! The test rationale is structural rather than statistical: under
//! `Disabled`, [`run_phase_a_cycle`] is implemented as
//! `run_baseline(cfg)` plus a wrap into `CycleOutput`. No additional
//! RNG consumption, no extra allocation, no Field2D mutation outside
//! what `run_baseline` does internally. The byte-equal contract is
//! therefore inherited from `run_baseline` determinism (see
//! `disabled_runs_are_bit_deterministic` of `v2_step8_regression_smoke`
//! for the precedent pattern).
//!
//! The test would catch any future contamination of the Disabled path
//! (e.g., an accidental allocation that happens to consume entropy
//! from the global RNG, a side-effect on `BaselineConfig`'s internal
//! state that flows into the next call).

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::harness::{run_baseline, BaselineConfig};
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::workflow::{
    run_phase_a_cycle, run_phase_a_loop, run_phase_b, WorkflowConfig,
};

fn build_test_config(scratch_subdir: &str) -> BaselineConfig {
    let scales = Scales::default();
    let mut cfg = BaselineConfig::dynamic_accidented_defaults(&scales);
    // Small grid + few steps so the test runs fast — this is a
    // structural regression check, not a physics validation, so
    // physical fidelity is irrelevant.
    cfg.grid_nx = 32;
    cfg.grid_ny = 32;
    cfg.steps = 20;
    cfg.total_time_nondim = 0.4;
    cfg.heightmap_fractions = Vec::new();
    cfg.output_dir = PathBuf::from(format!(
        "target/v2_workflow_disabled_regression/{}",
        scratch_subdir
    ));
    cfg
}

#[test]
fn workflow_disabled_run_phase_a_cycle_is_bit_identical_to_run_baseline() {
    let cfg_a = build_test_config("path_a");
    let cfg_b = build_test_config("path_b");

    let r_a = run_baseline(&cfg_a);
    let cycle_b = run_phase_a_cycle(&cfg_b, &WorkflowConfig::Disabled);

    let s_a = r_a.final_state.s_field.data();
    let s_b = cycle_b.baseline.final_state.s_field.data();
    assert_eq!(
        s_a, s_b,
        "s_field bytes must match between run_baseline and run_phase_a_cycle(Disabled)"
    );
    assert_eq!(r_a.final_state.vx, cycle_b.baseline.final_state.vx, "vx bytes mismatch");
    assert_eq!(r_a.final_state.vy, cycle_b.baseline.final_state.vy, "vy bytes mismatch");

    // Disabled never engages the low-res erosion path.
    assert_eq!(cycle_b.common.erosion_volume_removed, 0.0);

    // Cherry-picked metrics — full metrics struct is not Eq, so we
    // compare scalar invariants that suffice to flag any drift.
    assert_eq!(r_a.metrics.vmax_peak, cycle_b.baseline.metrics.vmax_peak);
    assert_eq!(r_a.metrics.cg_iter_mean, cycle_b.baseline.metrics.cg_iter_mean);
}

#[test]
fn workflow_disabled_run_phase_a_loop_returns_single_passthrough_cycle() {
    // Phase 4 changed run_phase_a_loop's signature to `&mut
    // BaselineConfig` because the Enabled branch mutates
    // `cfg.continuation` between cycles. The Disabled branch still
    // does not mutate cfg — the regression contract holds.
    let mut cfg = build_test_config("loop_single");
    let output = run_phase_a_loop(&mut cfg, &WorkflowConfig::Disabled);
    assert_eq!(
        output.cycles.len(),
        1,
        "Disabled loop must collapse to a single cycle"
    );
    assert_eq!(output.cycles[0].common.erosion_volume_removed, 0.0);
    // cfg.continuation is unchanged under Disabled — None remains None.
    assert!(
        cfg.continuation.is_none(),
        "Disabled loop must not mutate cfg.continuation"
    );
}

#[test]
fn workflow_disabled_run_phase_b_returns_none() {
    // Phase 5 added a `seed` parameter to run_phase_b for the
    // Enabled-branch FBM + erosion RNG. The Disabled branch ignores
    // it (returns None unconditionally), so the regression contract
    // is unchanged.
    let cfg = build_test_config("phase_b_skip");
    let r = run_baseline(&cfg);
    assert!(
        run_phase_b(&r.final_state.s_field, &WorkflowConfig::Disabled, cfg.seed).is_none(),
        "Disabled Phase B must be skipped (no allocation, no FBM seed consumption)"
    );
}
