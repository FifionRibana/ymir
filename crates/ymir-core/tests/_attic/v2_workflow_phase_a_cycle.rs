//! Step 12 Phase 3 acceptance — Phase A single-cycle orchestration.
//!
//! Two tests pinned to the issue's acceptance criteria:
//!
//! - **#1 / `v2_workflow_continuation_no_transient`** — D3 contract.
//!   Cycle 2 step 1 with continuation warm-start gives `peak|v|`
//!   within 10% of cycle 1's max. If a transient persists, it would
//!   indicate `ContinuationState` is missing a state field
//!   (m_subducted, slab buffer, …) that the warm-start needs.
//! - **#4 / `v2_workflow_cratonic_recompute`** — D4 retention rule.
//!   Two sub-tests:
//!   - *Loses craton on flipped plate*: a stringent
//!     `plate_area_min` (`0.99` within-plate fraction) flips every
//!     plate to oceanic-equivalent on the first cycle's recompute,
//!     erasing all cratons.
//!   - *Preserves stable craton*: default `plate_area_min = 0.10`
//!     under mild erosion keeps all plates retained, and the
//!     `craton_recomputation_change` metric reads `0.0` (the BFS
//!     distances are identical because no plate flipped).

use std::path::PathBuf;

use ymir_core::tectonics_v2::age_field::AgeFieldConfig;
use ymir_core::tectonics_v2::basal_drag::BasalDragConfig;
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::{CratonicConfig, CratonicConfigEnabled};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, BaselineConfig, ForceKind, NonlinearChoice,
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
    final_state_to_continuation, run_phase_a_cycle, PhaseAParams, WorkflowConfig, WorkflowParams,
};

fn build_minimal_cycle_config(steps: usize, scratch: &str) -> BaselineConfig {
    let scales = Scales::default();
    let mut cfg = BaselineConfig::dynamic_accidented_defaults(&scales);
    cfg.grid_nx = 32;
    cfg.grid_ny = 32;
    cfg.steps = steps;
    cfg.total_time_nondim = 0.4 * (steps as f64) / 20.0;
    cfg.heightmap_fractions = Vec::new();
    cfg.output_dir = PathBuf::from(format!("target/v2_workflow_phase3/{}", scratch));
    cfg
}

fn build_cratonic_cycle_config(scratch: &str, crcfg: CratonicConfigEnabled) -> BaselineConfig {
    let nx = 32;
    let ny = 32;
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    // 4 plates, 50% continental → expect ~2 continental plates with
    // observable cratons under default cratonic config.
    let vcfg = VoronoiConfig { num_plates: 4, continental_ratio: 0.5 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        nx,
        ny,
        &vcfg,
        42,
        rates,
        RecyclingConfig::default(),
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: 42,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: 3,
        cfl_factor: 0.3,
        total_time_nondim: 0.06,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from(format!("target/v2_workflow_phase3/{}", scratch)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Disabled,
        basal_drag: BasalDragConfig::Disabled,
        boundary,
        boundary_layout_name: "voronoi_seed42_n4".into(),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: CratonicConfig::Enabled(crcfg),
        age_field: AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: InitMode::Checkerboard,
        continuation: None,
        plate_kinematic: PlateKinematicConfig::Zero,
    }
}

#[test]
fn v2_workflow_continuation_no_transient() {
    // Minimal config: GPE + InitMode::Checkerboard, no cratonic, no
    // boundary, no slab, no mantle. Erosion runs on s_field but the
    // reclassify + craton recompute branches no-op (no plate_id /
    // plate_type populated under BoundaryConfig::Disabled).
    let wf = WorkflowConfig::Enabled(WorkflowParams {
        phase_a: PhaseAParams { alpha: 0.005, ..Default::default() },
        phase_b: Default::default(),
    });

    // Cycle 1 cold start, 5 steps.
    let cfg_1 = build_minimal_cycle_config(5, "cycle_1");
    let cycle_1 = run_phase_a_cycle(&cfg_1, &wf);
    let vmax_1 = cycle_1.baseline.metrics.vmax_peak;
    assert!(
        vmax_1 > 0.0,
        "cycle 1 must produce non-trivial dynamics: vmax_peak = {vmax_1}"
    );

    // Cycle 2: 1 step from continuation. With the warm-start contract
    // honoured, vmax_2 = peak|v| at step 1 of cycle 2 ≈ peak|v| in
    // cycle 1 (within 10%). Without warm-start (cold start), vmax_2
    // would either spike (transient) or dip (re-ramping) by far
    // more than 10%.
    let mut cfg_2 = build_minimal_cycle_config(1, "cycle_2");
    cfg_2.continuation = Some(final_state_to_continuation(&cycle_1.baseline.final_state));
    let cycle_2 = run_phase_a_cycle(&cfg_2, &wf);
    let vmax_2 = cycle_2.baseline.metrics.vmax_peak;

    let rel_diff = (vmax_2 - vmax_1).abs() / vmax_1;
    assert!(
        rel_diff < 0.10,
        "continuation transient detected: cycle 2 step 1 vmax = {vmax_2:.6e}, \
         cycle 1 vmax = {vmax_1:.6e}, relative diff = {rel_diff:.3} > 10%. \
         Investigate ContinuationState completeness — possibly missing \
         m_subducted, slab buffer, or other auxiliary state."
    );
}

#[test]
fn v2_workflow_cratonic_recompute_flips_eroded_plate() {
    // Step 12 Phase 3.5 — true D4 flip test. The two semantics are
    // now distinct parameters:
    //   - `plate_area_min = 0.10` (default, fraction-of-domain) →
    //     init's `build_cratonic_factor_field` retains every plate
    //     (each ~25% of the domain, well above 10%) and populates
    //     a non-trivial cratonic_factor on the continental plates.
    //   - `craton_retention_threshold = 0.95` (strict, within-plate)
    //     → the recompute requires 95% of each plate's cells to stay
    //     continental for the plate to keep its craton.
    //
    // Heavy erosion + 3 tectonic steps with k_sub=0.5 drop the
    // within-plate continental fraction below 0.95 for at least one
    // continental plate. That plate flips per the D4 retention rule
    // → BFS sources its cells → cratonic_factor = 0 on the flipped
    // plate. The Phase 3 commit message documented this scenario as
    // "untestable in single-cycle with overloaded parameter"; the
    // 3.5 disambiguation makes it a clean acceptance.
    let crcfg = CratonicConfigEnabled {
        plate_area_min: 0.10,
        craton_retention_threshold: 0.95,
        ..Default::default()
    };
    let cfg = build_cratonic_cycle_config("flip_eroded", crcfg);
    // alpha = 5.0 is aggressively non-physical — the test validates
    // the *mechanism* (D4 flip wiring), not realistic erosion
    // magnitudes. With realistic α (e.g. 0.01–0.05 from D8) a flip
    // emerges across multiple cycles; the multi-cycle natural-flip
    // scenario is Phase 4 integration territory.
    //
    // Step 12 R3 — `isostatic_rebound_ratio: 0.0` recovers the legacy
    // pre-R3 `β = 0` semantics for this wiring test (no rebound
    // compensation, full eroded-mass migration). The default 0.80
    // would attenuate the per-cell change by 5× and the
    // aggressive-but-bounded `α = 5.0` would no longer cross the
    // 0.95 retention threshold.
    let wf = WorkflowConfig::Enabled(WorkflowParams {
        phase_a: PhaseAParams {
            alpha: 5.0,
            isostatic_rebound_ratio: 0.0,
            ..Default::default()
        },
        phase_b: Default::default(),
    });

    let cycle = run_phase_a_cycle(&cfg, &wf);

    // The change metric must register > 0 — at least one plate
    // flipped. `measure_craton_change` counts cells whose factor
    // moved by more than 1e-9; a plate flip changes every cell of
    // that plate (the BFS source set rearranges), so the change
    // metric scales with the flipped plate's area-fraction of the
    // domain.
    let change = cycle
        .craton_recomputation_change
        .expect("change metric must be populated under CratonicConfig::Enabled");
    assert!(
        change > 0.0,
        "D4 retention rule must fire under strict craton_retention_threshold \
         + aggressive erosion: change = {change}. Investigate \
         continental_fraction per plate post-cycle if this surfaces."
    );
}

#[test]
fn v2_workflow_cratonic_recompute_preserves_stable_craton() {
    // Default `plate_area_min = 0.10` and mild erosion: continental
    // fraction in every retained plate stays ≥ 90 %, well above the
    // 10 % threshold. No plate flips → recomputed `per_plate_type`
    // matches the initial Voronoï → BFS distances identical →
    // cratonic_factor reproduced bit-for-bit. The change metric
    // reads `0.0`.
    let crcfg = CratonicConfigEnabled::default();
    let cfg = build_cratonic_cycle_config("cratonic_stable", crcfg);
    let wf = WorkflowConfig::Enabled(WorkflowParams {
        phase_a: PhaseAParams { alpha: 0.001, ..Default::default() },
        phase_b: Default::default(),
    });

    let cycle = run_phase_a_cycle(&cfg, &wf);
    let new_factor = cycle
        .baseline
        .final_state
        .cratonic_factor
        .as_ref()
        .expect("cratonic_factor must be populated");

    // Stable craton: at least one cell retains a high factor (deep
    // interior of a continental plate).
    let max_factor: f64 = new_factor.data().iter().copied().fold(0.0_f64, f64::max);
    assert!(
        max_factor > 0.5,
        "stable craton must retain interior factor > 0.5: got max = {max_factor}"
    );

    // BFS distances unchanged → recomputed factor identical.
    let change = cycle
        .craton_recomputation_change
        .expect("craton_recomputation_change must be populated");
    assert!(
        change < 1e-12,
        "mild erosion must not trigger any recompute change (no plate flip): \
         got {change}"
    );
}
