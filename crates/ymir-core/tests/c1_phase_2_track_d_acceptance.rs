//! Issue #132 Phase 2 Track D Stage A — formal acceptance tests
//! for the Track D boundary-evolution closure stack.
//!
//! Three tests at 64² × 300 steps with the Phase 2 R7 init +
//! Phase 1.1 kinematics, full Phase 2 + Track D stack:
//!
//! 1. [`boundary_events_fire_correctly_300_steps`] — Q-V.1
//!    Option B threshold: BOTH subduction AND accretion must fire
//!    at least once. Rifting splits NOT required (genuinely rare
//!    chewing-gum cut events, 0-3/seed per Stage V evidence).
//!    Architectural finding logged if either is zero (vs Stage V
//!    5/5 multi-seed evidence).
//! 2. [`ninth_bit_identical_preservation_phase_2_r7`] — regression
//!    guard reproducing the 9th bit-identical decomposition
//!    contract at Phase 2 R7 init + Track D disabled scope.
//!    Pattern from `c1_phase_a_decomposes_into_closures_then_post_tectonic`
//!    in `c1_phase_1_3_workflow.rs` (8th bit-identical contract);
//!    Stage E0 invariant preserved through Track D scope addition.
//! 3. [`acceptance_phase_2_gate_seed_diversity_300_steps`] —
//!    Phase 2 milestone gate proxy: 3 seeds × 300 steps, mean
//!    pairwise plate_id divergence > 30 %. Multi-seed variety
//!    signal forward toward §7.2 cross-track gate.

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::boundary_classification::classify_boundaries;
use ymir_core::tectonics_c1::closures::accretion::{
    apply_accretion_step, AccretionParams, ConvergenceTracker,
};
use ymir_core::tectonics_c1::closures::rifting::{
    apply_rifting_split, apply_rifting_thinning, DivergenceTracker, RiftingParams,
};
use ymir_core::tectonics_c1::closures::subduction::{
    apply_subduction_step, SubductionParams,
};
use ymir_core::tectonics_c1::init_r7::{init_c1_state_phase_2_r7, Phase2InitParams};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::cratonic::CratonicConfigEnabled;
use ymir_core::tectonics_v2::workflow::phase_a_common::{
    apply_post_tectonic, extract_per_plate_type, PostTectonicInput,
};
use ymir_core::tectonics_v2::workflow::{
    run_phase_a_cycle_c1, PhaseACycleInputC1, WorkflowConfig, WorkflowParams,
};

const GRID: usize = 64;
const N_STEPS: usize = 300;

fn make_config() -> C1TimeLoopConfig {
    C1TimeLoopConfig {
        n_steps: N_STEPS,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: IsostasyConfig::default(),
        drainage_max_distance: 30,
    }
}

fn phase_2_r7_init(seed: u64) -> (C1State, PlateKinematics) {
    let state = init_c1_state_phase_2_r7(GRID, seed, &Phase2InitParams::default());
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    (state, kinematics)
}

fn track_d_disabled_closures() -> C1Closures {
    C1Closures {
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
        ..C1Closures::default()
    }
}

// =========================================================================
// Test 1 — Q-V.1 Option B: both subduction AND accretion must fire
// =========================================================================

#[test]
fn boundary_events_fire_correctly_300_steps() {
    let seed = 42;
    let (mut state, mut kinematics) = phase_2_r7_init(seed);
    let closures = C1Closures::default();

    // Accumulate per-event-type counts via the Track-D-only
    // pipeline (parallel-state pattern from Stage V).
    let dt = 0.5 * (1.0 / GRID as f64) / kinematics.max_velocity().max(1e-12);

    let mut convergence_tracker = ConvergenceTracker::new();
    let mut divergence_tracker = DivergenceTracker::new();
    let mut total_subduction_cells = 0_usize;
    let mut total_arc_distributed = 0.0_f64;
    let mut total_reassignments = 0_usize;
    let mut total_accretion_merges = 0_usize;
    let mut total_rifting_thinned = 0_usize;
    let mut total_thinning_mass = 0.0_f64;
    let mut total_rifting_splits = 0_usize;
    let mut total_age_zeroed = 0_usize;

    for _step in 0..N_STEPS {
        let boundary = classify_boundaries(&state.plate_id, &kinematics);

        let sub = apply_subduction_step(
            &mut state.s,
            &mut state.plate_id,
            &mut state.plate_type,
            &boundary,
            &kinematics,
            &closures.subduction,
            dt,
        );
        total_subduction_cells += sub.cells_consumed;
        total_arc_distributed += sub.arc_mass_distributed;
        total_reassignments += sub.plate_ids_reassigned;

        let thin = apply_rifting_thinning(
            &mut state.s,
            &state.plate_type,
            &state.plate_id,
            &boundary,
            &kinematics,
            &closures.rifting,
            dt,
        );
        total_rifting_thinned += thin.cells_thinned;
        total_thinning_mass += thin.total_mass_removed;

        convergence_tracker.update(&state.plate_id, &kinematics);
        divergence_tracker.update(&state.plate_id, &kinematics);

        let acc = apply_accretion_step(
            &mut state.plate_id,
            &state.s,
            &mut kinematics,
            &convergence_tracker,
            &closures.accretion,
        );
        total_accretion_merges += acc.merges_count;

        let split = apply_rifting_split(
            &mut state.plate_id,
            &state.plate_type,
            &mut state.age,
            &state.s,
            &mut kinematics,
            &divergence_tracker,
            &closures.rifting,
        );
        total_rifting_splits += split.splits_count;
        total_age_zeroed += split.age_zeroed_cells;
    }

    eprintln!(
        "Stage A Test 1 — boundary_events_fire_correctly_300_steps (seed = {seed}, grid = {GRID}², steps = {N_STEPS}):"
    );
    eprintln!("  Subduction cells consumed   = {total_subduction_cells}");
    eprintln!("  Subduction arc distributed  = {:.4}", total_arc_distributed);
    eprintln!("  Subduction reassignments    = {total_reassignments}");
    eprintln!("  Accretion merges            = {total_accretion_merges}");
    eprintln!("  Rifting cells thinned       = {total_rifting_thinned}");
    eprintln!("  Rifting thinning mass       = {:.4}", total_thinning_mass);
    eprintln!("  Rifting splits              = {total_rifting_splits}");
    eprintln!("  Path 3.B age zeroed cells   = {total_age_zeroed}");
    eprintln!();
    eprintln!("Q-V.1 Option B threshold: subduction > 0 AND accretion > 0");

    if total_subduction_cells == 0 {
        eprintln!(
            "ARCHITECTURAL FINDING: zero subduction events at seed {seed} contradicts Stage V evidence (5/5 seeds with events). Investigate."
        );
    }
    if total_accretion_merges == 0 {
        eprintln!(
            "ARCHITECTURAL FINDING: zero accretion merges at seed {seed} contradicts Stage V evidence (6-10 merges per seed). Investigate."
        );
    }
    if total_rifting_splits == 0 {
        eprintln!(
            "  Note: zero rifting splits at seed {seed} (expected per Stage V evidence: 0-3/seed; rare chewing-gum cut)."
        );
    }

    assert!(
        total_subduction_cells > 0,
        "subduction must fire at least once over {N_STEPS} steps at seed {seed}"
    );
    assert!(
        total_accretion_merges > 0,
        "accretion must merge at least once over {N_STEPS} steps at seed {seed}"
    );
}

// =========================================================================
// Test 2 — 9th bit-identical decomposition preservation
//
// Reproduces the wrapper-equals-decomposition contract from
// `c1_phase_a_decomposes_into_closures_then_post_tectonic`
// (c1_phase_1_3_workflow.rs) at Phase 2 R7 init + Track D
// disabled scope. The contract: `run_phase_a_cycle_c1(input,
// Enabled(_))` is byte-for-byte equal to `run_with_closures(...)
// + apply_post_tectonic(...)` when Track D is off (= Track A/B
// regime). Stage E4's Track D wiring did NOT break the
// decomposition.
// =========================================================================

#[test]
fn ninth_bit_identical_preservation_phase_2_r7() {
    let seed = 42;
    let n_steps_short = 50; // Short run keeps the test sub-second.

    let closures = track_d_disabled_closures();
    // Issue #141: exercise the decomposition CONTRACT under the C1
    // production sea-level mode (P95-cap). The test is wrapper ==
    // decomposition (path-A == path-B), so the S̃ values are redefined
    // (P95-cap) but the CONTRACT must stay byte-exact; if it breaks,
    // P95-cap introduced a decomposition inconsistency (real bug).
    let iso_config = IsostasyConfig::c1_default();
    let time_loop_config = C1TimeLoopConfig {
        n_steps: n_steps_short,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    let workflow = WorkflowConfig::Enabled(WorkflowParams::default());

    // Path A — wrapper.
    let (mut state_a, mut kinematics_a) = phase_2_r7_init(seed);
    let _output_a = run_phase_a_cycle_c1(
        PhaseACycleInputC1 {
            state: &mut state_a,
            kinematics: &mut kinematics_a,
            closures: &closures,
            time_loop_config: &time_loop_config,
            iso_config: &iso_config,
            cratonic_config: None,
        },
        &workflow,
    );

    // Path B — manual decomposition (run_with_closures +
    // apply_post_tectonic).
    let (mut state_b, mut kinematics_b) = phase_2_r7_init(seed);
    run_with_closures(
        &mut state_b,
        &mut kinematics_b,
        &time_loop_config,
        &closures,
        |_, _| {},
    );
    let original_per_plate_type =
        extract_per_plate_type(&state_b.plate_id, &state_b.plate_type);
    let cratonic_cfg_b: Option<&CratonicConfigEnabled> = None;
    let _post_b = apply_post_tectonic(PostTectonicInput {
        s_field: &mut state_b.s,
        plate_id: Some(&state_b.plate_id),
        plate_type: Some(&mut state_b.plate_type),
        previous_cratonic_factor: None,
        initial_per_plate_type: Some(&original_per_plate_type),
        params: match &workflow {
            WorkflowConfig::Enabled(p) => &p.phase_a,
            WorkflowConfig::Disabled => panic!("workflow must be Enabled for this test"),
        },
        iso_cfg: &iso_config,
        cratonic_cfg: cratonic_cfg_b,
    });

    // Bit-identical S̃ comparison (no tolerance).
    let max_abs_delta = state_a
        .s
        .data()
        .iter()
        .zip(state_b.s.data().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);

    eprintln!(
        "Stage A Test 2 — ninth_bit_identical_preservation_phase_2_r7 (seed = {seed}, steps = {n_steps_short}):"
    );
    eprintln!("  max |S̃_wrapper − S̃_decomposition| = {max_abs_delta}");

    assert_eq!(
        state_a.s.data(),
        state_b.s.data(),
        "Phase 2 R7 + Track-D-disabled: wrapper vs decomposition must be byte-identical"
    );
    assert_eq!(state_a.plate_id.data(), state_b.plate_id.data());
}

// =========================================================================
// Test 3 — Phase 2 milestone gate proxy: multi-seed plate_id divergence
// =========================================================================

#[test]
fn acceptance_phase_2_gate_seed_diversity_300_steps() {
    let seeds: [u64; 3] = [42, 1337, 2026];
    let closures = C1Closures::default();
    let config = make_config();

    // Per-seed final states.
    let mut final_states: Vec<(u64, Vec<u16>)> = Vec::new();
    for &seed in seeds.iter() {
        let (mut state, mut kinematics) = phase_2_r7_init(seed);
        run_with_closures(
            &mut state,
            &mut kinematics,
            &config,
            &closures,
            |_, _| {},
        );
        final_states.push((seed, state.plate_id.data().to_vec()));
    }

    // Pairwise plate_id divergence: fraction of cells where the
    // two seeds' final plate_id differ.
    let total_cells = GRID * GRID;
    let mut divergences: Vec<(u64, u64, f64)> = Vec::new();
    for i in 0..final_states.len() {
        for j in (i + 1)..final_states.len() {
            let (seed_i, plate_id_i) = &final_states[i];
            let (seed_j, plate_id_j) = &final_states[j];
            let mismatch = plate_id_i
                .iter()
                .zip(plate_id_j.iter())
                .filter(|(a, b)| a != b)
                .count();
            let divergence = mismatch as f64 / total_cells as f64;
            divergences.push((*seed_i, *seed_j, divergence));
        }
    }

    let mean_divergence: f64 =
        divergences.iter().map(|&(_, _, d)| d).sum::<f64>() / divergences.len() as f64;

    eprintln!(
        "Stage A Test 3 — acceptance_phase_2_gate_seed_diversity_300_steps (grid = {GRID}², steps = {N_STEPS}):"
    );
    eprintln!("  Seeds: {seeds:?}");
    eprintln!("  Pairwise plate_id divergence:");
    for &(si, sj, d) in &divergences {
        eprintln!("    seed {si:>5} ↔ {sj:>5}: {:.1} %", d * 100.0);
    }
    eprintln!("  Mean pairwise divergence = {:.1} %", mean_divergence * 100.0);
    eprintln!("  Phase 2 milestone gate proxy threshold: > 30 %");

    if mean_divergence > 0.80 {
        eprintln!(
            "  Note: divergence > 80 % — threshold easily exceeded; consider raising if Phase 3+ wants stronger evidence."
        );
    } else if mean_divergence < 0.20 {
        eprintln!(
            "  ARCHITECTURAL FINDING: mean divergence {:.1} % below 20 %. Track D mutations may be too subtle to differentiate final plate_id across seeds. Investigate.",
            mean_divergence * 100.0
        );
    }

    assert!(
        mean_divergence > 0.30,
        "Phase 2 milestone gate proxy: mean pairwise plate_id divergence {:.1} % below 30 % threshold",
        mean_divergence * 100.0
    );
}
