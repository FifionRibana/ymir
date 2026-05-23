//! Issue #125 Phase 1.3 H3 — integration tests for the C1 Phase A
//! path through the workflow wrapper.
//!
//! Lives under default features (no `v2_legacy` requirement) — the
//! C1 path is the mainline paradigm; the v2 path's parallel
//! regression tests live in `tests/_attic/v2_workflow_*`.
//!
//! Test inventory (per H3 user spec):
//!
//! - [`c1_phase_a_cycle_smoke_integration`] — public API smoke
//!   through the `tests/` directory. Catches missing pub re-
//!   exports that the in-module test in `phase_a_c1.rs` would
//!   silently bypass.
//! - [`c1_phase_a_decomposes_into_closures_then_post_tectonic`] —
//!   the C1 analogue of v2's bit-identical Disabled regression.
//!   Verifies that `run_phase_a_cycle_c1(input, Enabled(_))` is
//!   *exactly* `run_with_closures(state, …)` followed by
//!   `apply_post_tectonic(PostTectonicInput { … })` —
//!   "wrapper = building blocks, nothing more, nothing less".
//!   Bit-identical equality, no tolerance.
//! - [`c1_phase_a_with_disabled_closures_matches_phase_1_1`] —
//!   with all closures off and `wf = Disabled`, the wrapper
//!   reduces to plain Phase 1.1 advection; mass conservation
//!   `< 1e-10` must hold.
//! - [`c1_phase_a_with_cratonic_enabled_produces_factor`] — sanity
//!   check that the cratonic path is exercisable from C1, produces
//!   a `Field2D` of correct shape with values in `[0, 1]`.
//!
//! These are H3 smoke + regression tests, **not** Phase 1.3
//! acceptance tests — those land in Stage E3 (4 invariants on
//! equilibrium-height + Davis-Suppe interaction at 64²×300).

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams;
use ymir_core::tectonics_c1::closures::equilibrium_height::params::EquilibriumHeightParams;
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::cratonic::CratonicConfigEnabled;
use ymir_core::tectonics_v2::workflow::phase_a_common::{
    apply_post_tectonic, extract_per_plate_type, PostTectonicInput,
};
use ymir_core::tectonics_v2::workflow::{
    run_phase_a_cycle_c1, PhaseACycleInputC1, WorkflowConfig, WorkflowParams,
};

const GRID: usize = 32;
const SEED: u64 = 42;

fn make_time_loop_config(n_steps: usize) -> C1TimeLoopConfig {
    C1TimeLoopConfig {
        n_steps,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
    }
}

#[test]
fn c1_phase_a_cycle_smoke_integration() {
    // Smoke through the public API only. Lighter than the inline
    // Commit 4 test (32² × 50 steps vs 64² × 300) so it stays
    // sub-second under default `cargo test`.
    let mut state = init_c1_state_phase_1_1(GRID, SEED);
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let time_loop_config = make_time_loop_config(50);
    let iso_config = IsostasyConfig::default();
    let wf = WorkflowConfig::Enabled(WorkflowParams::default());

    let output = run_phase_a_cycle_c1(
        PhaseACycleInputC1 {
            state: &mut state,
            kinematics: &kinematics,
            closures: &closures,
            time_loop_config: &time_loop_config,
            iso_config: &iso_config,
            cratonic_config: None,
        },
        &wf,
    );

    assert!(
        state.s.data().iter().all(|v| v.is_finite()),
        "no NaN/Inf in S̃ after smoke run"
    );
    assert!(
        output.common.sea_level_normalized > 0.0,
        "Enabled post-pass must compute sea_level_normalized; got {}",
        output.common.sea_level_normalized
    );
}

#[test]
fn c1_phase_a_decomposes_into_closures_then_post_tectonic() {
    // The C1 analogue of v2's bit-identical Disabled regression:
    // `run_phase_a_cycle_c1(input, Enabled(_))` must equal
    // `run_with_closures` + `apply_post_tectonic` byte-for-byte.
    // No tolerance — exact equality on the S̃ buffer and on all
    // five `CycleOutputCommon` scalars.
    //
    // Surfaces order-sensitivity / hidden state in the wrapper if
    // the equality breaks. Per H3 W2: if numerical noise > 1e-15
    // appears, surface as architectural finding.

    let n_steps = 30;
    let kinematics_template = {
        let scratch = init_c1_state_phase_1_1(GRID, SEED);
        PlateKinematics::preset_phase_1_1(scratch.num_plates)
    };
    let closures = C1Closures::default();
    let time_loop_config = make_time_loop_config(n_steps);
    let iso_config = IsostasyConfig::default();
    let wf_params = WorkflowParams::default();
    let phase_a_params = wf_params.phase_a.clone();
    let wf = WorkflowConfig::Enabled(wf_params);

    // Path A — through the workflow wrapper.
    let mut state_a = init_c1_state_phase_1_1(GRID, SEED);
    let output_a = run_phase_a_cycle_c1(
        PhaseACycleInputC1 {
            state: &mut state_a,
            kinematics: &kinematics_template,
            closures: &closures,
            time_loop_config: &time_loop_config,
            iso_config: &iso_config,
            cratonic_config: None,
        },
        &wf,
    );

    // Path B — manual decomposition: run_with_closures then
    // apply_post_tectonic with the same paradigm-agnostic input
    // bundle the wrapper builds internally.
    let mut state_b = init_c1_state_phase_1_1(GRID, SEED);
    run_with_closures(
        &mut state_b,
        &kinematics_template,
        &time_loop_config,
        &closures,
        |_, _| {},
    );
    let original_per_plate_type = extract_per_plate_type(&state_b.plate_id, &state_b.plate_type);
    let output_b = apply_post_tectonic(PostTectonicInput {
        s_field: &mut state_b.s,
        plate_id: Some(&state_b.plate_id),
        plate_type: Some(&mut state_b.plate_type),
        previous_cratonic_factor: None,
        initial_per_plate_type: Some(&original_per_plate_type),
        params: &phase_a_params,
        iso_cfg: &iso_config,
        cratonic_cfg: None,
    });

    // Bit-identical state.
    assert_eq!(
        state_a.s.data(),
        state_b.s.data(),
        "S̃ field must be bit-identical between wrapper and manual decomposition"
    );

    // Bit-identical common scalars (no order-sensitivity allowed).
    assert_eq!(
        output_a.common.erosion_volume_removed, output_b.common.erosion_volume_removed,
        "erosion_volume_removed must be bit-identical"
    );
    assert_eq!(
        output_a.common.erosion_peak_delta_h, output_b.common.erosion_peak_delta_h,
        "erosion_peak_delta_h must be bit-identical"
    );
    assert_eq!(
        output_a.common.sea_level_normalized, output_b.common.sea_level_normalized,
        "sea_level_normalized must be bit-identical"
    );
    assert_eq!(
        output_a.common.mass_drift, output_b.common.mass_drift,
        "mass_drift must be bit-identical"
    );
    assert_eq!(
        output_a.common.craton_recomputation_change, output_b.common.craton_recomputation_change,
        "craton_recomputation_change must be bit-identical"
    );
}

#[test]
fn c1_phase_a_with_disabled_closures_matches_phase_1_1() {
    // C1 parallel of v2's `workflow_disabled_run_phase_a_cycle_is_
    // bit_identical_to_run_baseline` regression. Two layers:
    //
    //   (a) `wf = Disabled` → wrapper degenerates to a direct
    //       `run_with_closures` call (no post-pass). The two
    //       paths must produce bit-identical S̃ buffers.
    //   (b) With both closures disabled, the C1 tectonic step is
    //       advection-only, which is Phase 1.1 — mass conservation
    //       at machine-precision floor (1e-10 relative budget,
    //       same as Phase 1.1 invariant).
    let closures_off = C1Closures {
        davis_suppe: DavisSuppeParams {
            enabled: false,
            ..DavisSuppeParams::default()
        },
        equilibrium_height: EquilibriumHeightParams {
            enabled: false,
            ..EquilibriumHeightParams::default()
        },
    };
    let time_loop_config = make_time_loop_config(50);
    let iso_config = IsostasyConfig::default();

    let kinematics = {
        let scratch = init_c1_state_phase_1_1(GRID, SEED);
        PlateKinematics::preset_phase_1_1(scratch.num_plates)
    };

    // (a) Wrapper Disabled.
    let mut state_a = init_c1_state_phase_1_1(GRID, SEED);
    let initial_mass_a: f64 = state_a.s.data().iter().sum();
    let output_a = run_phase_a_cycle_c1(
        PhaseACycleInputC1 {
            state: &mut state_a,
            kinematics: &kinematics,
            closures: &closures_off,
            time_loop_config: &time_loop_config,
            iso_config: &iso_config,
            cratonic_config: None,
        },
        &WorkflowConfig::Disabled,
    );
    let final_mass_a: f64 = state_a.s.data().iter().sum();

    // Disabled output must be the default CycleOutputCommon.
    assert_eq!(output_a.common.sea_level_normalized, 0.0);
    assert_eq!(output_a.common.mass_drift, 0.0);
    assert_eq!(output_a.common.erosion_volume_removed, 0.0);
    assert_eq!(output_a.common.erosion_peak_delta_h, 0.0);
    assert!(output_a.common.craton_recomputation_change.is_none());
    assert!(output_a.new_cratonic_factor.is_none());

    // (b) Direct `run_with_closures` with the same closures_off.
    let mut state_b = init_c1_state_phase_1_1(GRID, SEED);
    let initial_mass_b: f64 = state_b.s.data().iter().sum();
    run_with_closures(
        &mut state_b,
        &kinematics,
        &time_loop_config,
        &closures_off,
        |_, _| {},
    );
    let final_mass_b: f64 = state_b.s.data().iter().sum();

    // Both paths must produce bit-identical S̃ — wrapper Disabled
    // is, literally, `run_with_closures` plus zero post-pass.
    assert_eq!(
        state_a.s.data(),
        state_b.s.data(),
        "Disabled wrapper must be bit-identical to direct run_with_closures"
    );

    // Phase 1.1 advection-only invariant: mass conserved at
    // machine-precision floor.
    let rel_drift_a = (final_mass_a - initial_mass_a).abs() / initial_mass_a;
    let rel_drift_b = (final_mass_b - initial_mass_b).abs() / initial_mass_b;
    assert!(
        rel_drift_a < 1e-10,
        "closures-off wrapper relative mass drift {rel_drift_a:.3e} exceeds 1e-10 — Phase 1.1 invariant broken"
    );
    assert!(
        rel_drift_b < 1e-10,
        "closures-off direct run relative mass drift {rel_drift_b:.3e} exceeds 1e-10 — Phase 1.1 invariant broken"
    );
}

#[test]
fn c1_phase_a_with_cratonic_enabled_produces_factor() {
    // Exercise the cratonic-recompute path: pass a real
    // `CratonicConfigEnabled` and verify the produced
    // `new_cratonic_factor` is `Some(Field2D)` with the right
    // shape and value range. Not a numerical / D4-rule test
    // (those are v2's `v2_workflow_cratonic_recompute_*` under
    // v2_legacy); just a "C1 can drive this codepath" sanity.
    let mut state = init_c1_state_phase_1_1(GRID, SEED);
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let time_loop_config = make_time_loop_config(30);
    let iso_config = IsostasyConfig::default();
    let cratonic_cfg = CratonicConfigEnabled::default();
    let wf = WorkflowConfig::Enabled(WorkflowParams::default());

    let output = run_phase_a_cycle_c1(
        PhaseACycleInputC1 {
            state: &mut state,
            kinematics: &kinematics,
            closures: &closures,
            time_loop_config: &time_loop_config,
            iso_config: &iso_config,
            cratonic_config: Some(&cratonic_cfg),
        },
        &wf,
    );

    let factor = output
        .new_cratonic_factor
        .expect("cratonic_config Enabled must produce a factor");
    assert_eq!(factor.nx(), GRID, "cratonic factor must match grid width");
    assert_eq!(factor.ny(), GRID, "cratonic factor must match grid height");
    for &v in factor.data() {
        assert!(
            (0.0..=1.0).contains(&v),
            "cratonic factor values must lie in [0, 1] (smoothstep range); got {v}"
        );
    }
}
