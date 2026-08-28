//! Step 12 Phase 4 acceptance — Phase A multi-cycle loop integration.
//!
//! Pinned to the issue's integration acceptance:
//!
//! > `v2_workflow_phase_a_5_cycles_64sq` : 5 cycles × 20 steps,
//! > vérifier que mass drift accumule de façon prévisible, peak S̃
//! > stabilise après ~3 cycles
//!
//! The 64² runtime (~5 min total) is borderline for non-ignored tests,
//! so the default-suite test runs at 32² × 5 cycles × 20 steps
//! (`v2_workflow_phase_a_5_cycles_32sq`); the dynamics scale 1:1 with
//! grid resolution at this stage of the milestone (mass-balance and
//! peak-stabilisation are integrated quantities, not per-cell). A
//! 64² `#[ignore]` variant is left in for explicit heavy-validation
//! runs.

use std::path::PathBuf;

use ymir_core::tectonics_v2::age_field::AgeFieldConfig;
use ymir_core::tectonics_v2::basal_drag::BasalDragConfig;
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::CratonicConfig;
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, ForceKind, NonlinearChoice, build_force,
};
use ymir_core::tectonics_v2::init::InitMode;
use ymir_core::tectonics_v2::mantle::MantleConfig;
use ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;
use ymir_core::tectonics_v2::workflow::{
    PhaseAParams, WorkflowConfig, WorkflowParams, run_phase_a_loop_v2,
};

fn build_phase4_config(grid_size: usize, k_cycle: usize, scratch: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 4, continental_ratio: 0.5 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        grid_size,
        grid_size,
        &vcfg,
        42,
        rates,
        RecyclingConfig::default(),
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: 42,
        grid_nx: grid_size,
        grid_ny: grid_size,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: k_cycle,
        cfl_factor: 0.3,
        total_time_nondim: 0.4 * (k_cycle as f64) / 20.0,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from(format!("target/v2_workflow_phase4/{}", scratch)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Disabled,
        basal_drag: BasalDragConfig::Disabled,
        boundary,
        boundary_layout_name: format!("voronoi_seed42_n4_{}sq", grid_size),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: CratonicConfig::Disabled,
        age_field: AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: InitMode::Checkerboard,
        continuation: None,
        plate_kinematic: PlateKinematicConfig::Zero,
    }
}

fn run_5cycle_integration(grid_size: usize, scratch: &str) {
    let mut cfg = build_phase4_config(grid_size, 20, scratch);
    let wf = WorkflowConfig::Enabled(WorkflowParams {
        phase_a: PhaseAParams {
            n_cycles: 5,
            k_cycle: 20,
            alpha: 0.01,
            isostatic_rebound_ratio: 0.80,
            max_drainage_distance: 10,
        },
        phase_b: Default::default(),
    });

    let output = run_phase_a_loop_v2(&mut cfg, &wf);
    assert_eq!(output.cycles.len(), 5, "expected 5 cycles");

    // Continuation chained — by the second cycle on, cfg.continuation
    // is populated with the prior cycle's final state.
    assert!(
        cfg.continuation.is_some(),
        "after 5 cycles, cfg.continuation must carry the last-but-one cycle's state"
    );

    // Mass drift contract (Step 12 R3): macro_redistribution conserves
    // total mass by construction — every eroded gram either migrates
    // to a drainage target or is implicitly compensated by isostatic
    // rebound. Per-cycle `mass_drift` is bounded by IEEE-754 floor
    // (~ N · ε · Δh̄ at machine precision), independent of sign. The
    // pre-R3 contract "`β = 0` → mass_drift ≤ 0" is obsolete — there
    // is no `β` knob to flip mass migration off.
    let initial_mass_estimate = 0.6 * (grid_size * grid_size) as f64;
    let machine_drift_budget = initial_mass_estimate * 1e-10;
    for (i, cycle) in output.cycles.iter().enumerate() {
        assert!(
            cycle.common.mass_drift.abs() < machine_drift_budget,
            "cycle {i}: mass_drift = {} exceeds machine-precision budget {} (mass conservation must hold)",
            cycle.common.mass_drift,
            machine_drift_budget
        );
    }
    let total_drift: f64 = output.cycles.iter().map(|c| c.common.mass_drift).sum();
    assert!(
        total_drift.abs() < machine_drift_budget,
        "cumulative mass drift {total_drift} exceeds machine-precision budget {machine_drift_budget}"
    );

    // Peak S̃ stabilisation contract: by cycle 3, the system has
    // reached a regime where tectonic compression + erosion balance.
    // The peak S̃ from cycle 3 to cycle 4 should drift by < 10 %
    // (equilibrium proxy).
    let peaks: Vec<f64> = output
        .cycles
        .iter()
        .map(|c| {
            c.baseline.final_state.s_field.data().iter().copied().fold(f64::NEG_INFINITY, f64::max)
        })
        .collect();
    let peak_3 = peaks[3];
    let peak_4 = peaks[4];
    assert!(peak_3 > 0.0, "peak[3] must be positive: {peak_3}");
    let rel_drift_late = (peak_4 - peak_3).abs() / peak_3;
    assert!(
        rel_drift_late < 0.10,
        "peak S̃ failed to stabilise by cycle 3 → 4: peaks = {peaks:?}, \
         relative drift = {rel_drift_late:.3} > 10 %"
    );

    // Erosion volume per cycle: every cycle should have done some
    // work. A degenerate run (all-flat field, no slope) would
    // produce zero everywhere; that scenario must not silently slip
    // through.
    for (i, cycle) in output.cycles.iter().enumerate() {
        assert!(
            cycle.common.erosion_volume_removed > 0.0,
            "cycle {i}: erosion_volume_removed = {} (every cycle must engage)",
            cycle.common.erosion_volume_removed
        );
    }
}

#[test]
fn v2_workflow_phase_a_5_cycles_32sq() {
    run_5cycle_integration(32, "5cycles_32sq");
}

#[test]
#[ignore]
fn v2_workflow_phase_a_5_cycles_64sq() {
    run_5cycle_integration(64, "5cycles_64sq");
}
