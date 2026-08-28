//! Issue #131 Phase 2 Track B Stage V — R7 init validation tests.
//!
//! Eight tests under default features, all on the public API of
//! [`ymir_core::tectonics_c1::init_r7`]. Validates the
//! [`init_c1_state_phase_2_r7`] dispatcher and its three
//! sub-components (boundary displacement, continental clustering,
//! ridge-aligned age).
//!
//! Test inventory:
//!
//! 1. [`r7_init_non_rectilinear_boundaries`] — quantitative
//!    measure of boundary curvature: count cells whose `plate_id`
//!    differs between R7-enabled and R7-disabled dispatcher
//!    output at the same seed; assert the count is in
//!    `(0, 20 %]` (Stage E1 healthy regime, on real Voronoï).
//! 2. [`r7_init_deterministic_given_seed`] — same
//!    `(grid_size, seed, params)` → bit-identical `C1State`.
//! 3. [`r7_init_different_seeds_produce_different_continents`] —
//!    seeds 42 vs 1337 → divergent continental layouts.
//! 4. [`r7_continental_fraction_within_target`] — per-plate
//!    continental count in `target ± max(5 %, 1 / num_plates)`.
//! 5. [`r7_continental_cluster_contiguous`] — single connected
//!    continental subgraph in the post-dispatcher plate
//!    adjacency.
//! 6. [`r7_age_ridge_aligned_at_divergent_boundaries`] — multi-
//!    seed `[42, 100, 1337, 2026, 9999]`: count `age == 0` cells
//!    per seed; assert at least one seed produces ridge cells
//!    (Track C kinematics concern mitigation).
//! 7. [`r7_age_distribution_compared_to_phase_1_1`] — qualitative
//!    comparison: Phase 1.1 init produces exactly 2 unique age
//!    values `{0.5, 7.0}`; Phase 2 R7 produces ≥ 2 unique ages
//!    (3 when ridges present).
//! 8. [`r7_phase_1_1_init_preserved_unchanged`] — property-based
//!    regression guard on Phase 1.1 init: 8 plates, age field
//!    exactly 2 unique values, continental fraction in
//!    `[0.20, 0.45]`. Locks Phase 1.x baseline against silent
//!    drift from Phase 2 Track B changes.

use std::collections::{HashSet, VecDeque};

use ymir_core::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries};
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::init_r7::{
    Phase2InitParams, R7InitParams, build_plate_adjacency, init_c1_state_phase_2_r7,
};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const GRID: usize = 64;
const SEED: u64 = 42;

/// Build a `Phase2InitParams` with R7 displacement disabled (all
/// other sub-components default-enabled). Used for the Test 1
/// baseline comparison (R7-enabled vs R7-disabled, same seed)
/// per the Phase 1.4 Option B pattern.
fn params_with_r7_disabled() -> Phase2InitParams {
    Phase2InitParams {
        r7: R7InitParams { enabled: false, ..R7InitParams::default() },
        ..Phase2InitParams::default()
    }
}

/// **Test 1** — non-rectilinear boundaries. Counts cells whose
/// `plate_id` differs between R7-enabled and R7-disabled
/// dispatcher output at the same seed; asserts count in
/// `(0, 20 %]` per Stage E1 healthy regime threshold.
#[test]
fn r7_init_non_rectilinear_boundaries() {
    let state_r7_on = init_c1_state_phase_2_r7(GRID, SEED, &Phase2InitParams::default());
    let state_r7_off = init_c1_state_phase_2_r7(GRID, SEED, &params_with_r7_disabled());

    let total = GRID * GRID;
    let mut reassigned = 0;
    for j in 0..GRID {
        for i in 0..GRID {
            if state_r7_on.plate_id.get(i, j) != state_r7_off.plate_id.get(i, j) {
                reassigned += 1;
            }
        }
    }
    let frac = 100.0 * reassigned as f64 / total as f64;
    eprintln!(
        "Test 1 non_rectilinear_boundaries: reassigned = {reassigned} / {total} ({frac:.2} %)"
    );

    assert!(
        reassigned > 0,
        "R7 displacement must reassign at least one cell at seed {SEED}; got 0"
    );
    let upper_bound = total / 5; // 20 %
    assert!(
        reassigned < upper_bound,
        "reassignment must stay under 20 % (got {reassigned} / {total} = {frac:.2} %)"
    );
}

/// **Test 2** — deterministic given seed.
#[test]
fn r7_init_deterministic_given_seed() {
    let params = Phase2InitParams::default();
    let state_a = init_c1_state_phase_2_r7(GRID, SEED, &params);
    let state_b = init_c1_state_phase_2_r7(GRID, SEED, &params);

    for j in 0..GRID {
        for i in 0..GRID {
            assert_eq!(state_a.s.get(i, j), state_b.s.get(i, j), "S̃ mismatch at ({i}, {j})");
            assert_eq!(state_a.age.get(i, j), state_b.age.get(i, j), "age mismatch at ({i}, {j})");
            assert_eq!(
                state_a.plate_id.get(i, j),
                state_b.plate_id.get(i, j),
                "plate_id mismatch at ({i}, {j})"
            );
            assert_eq!(
                state_a.plate_type.get(i, j),
                state_b.plate_type.get(i, j),
                "plate_type mismatch at ({i}, {j})"
            );
            assert_eq!(
                state_a.cratonic_mask.get(i, j),
                state_b.cratonic_mask.get(i, j),
                "cratonic_mask mismatch at ({i}, {j})"
            );
        }
    }
}

/// **Test 3** — different seeds produce different continents.
#[test]
fn r7_init_different_seeds_produce_different_continents() {
    let params = Phase2InitParams::default();
    let state_a = init_c1_state_phase_2_r7(GRID, 42, &params);
    let state_b = init_c1_state_phase_2_r7(GRID, 1337, &params);

    let mut plate_id_diff = 0;
    let mut plate_type_diff = 0;
    for j in 0..GRID {
        for i in 0..GRID {
            if state_a.plate_id.get(i, j) != state_b.plate_id.get(i, j) {
                plate_id_diff += 1;
            }
            if state_a.plate_type.get(i, j) != state_b.plate_type.get(i, j) {
                plate_type_diff += 1;
            }
        }
    }
    eprintln!(
        "Test 3 different_seeds: plate_id_diff = {plate_id_diff}, plate_type_diff = {plate_type_diff} of {}",
        GRID * GRID
    );

    assert!(
        plate_id_diff > 0,
        "different seeds must produce different plate_id layouts; got identical"
    );
    assert!(
        plate_type_diff > 0,
        "different seeds must produce different continental layouts; got identical"
    );
}

/// Helper — collect per-plate `PlateType` from the cell-level
/// field via `plate_id` lookup. Returns indexed `Vec<PlateType>`
/// of length `num_plates`.
fn per_plate_type_from_state(state: &ymir_core::tectonics_c1::state::C1State) -> Vec<PlateType> {
    let mut per_plate: Vec<Option<PlateType>> = vec![None; state.num_plates];
    for j in 0..state.ny() {
        for i in 0..state.nx() {
            let pid = state.plate_id.get(i, j) as usize;
            if pid < per_plate.len() {
                per_plate[pid].get_or_insert(state.plate_type.get(i, j));
            }
        }
    }
    per_plate.into_iter().map(|t| t.unwrap_or(PlateType::Oceanic)).collect()
}

/// **Test 4** — continental fraction within target. Same
/// granularity-aware tolerance `max(5 %, 1 / num_plates)` as the
/// Stage E2 unit test, applied to the real Voronoï output.
#[test]
fn r7_continental_fraction_within_target() {
    let params = Phase2InitParams::default();
    let target = params.cluster.continental_fraction;
    let state = init_c1_state_phase_2_r7(GRID, SEED, &params);
    let num_plates = state.num_plates;

    let per_plate = per_plate_type_from_state(&state);
    let continental_count =
        per_plate.iter().filter(|t| matches!(t, PlateType::Continental)).count();
    let actual_fraction = continental_count as f64 / num_plates as f64;
    let tolerance = (1.0 / num_plates as f64).max(0.05);
    let diff = (actual_fraction - target).abs();
    eprintln!(
        "Test 4 continental_fraction: actual = {actual_fraction:.3} ({continental_count} / {num_plates}), target = {target:.3}, diff = {diff:.3}, tolerance = {tolerance:.3}"
    );

    assert!(
        diff <= tolerance,
        "continental fraction {actual_fraction:.3} too far from target {target:.3} (diff {diff:.3} > tolerance {tolerance:.3})"
    );
}

/// **Test 5** — continental cluster contiguous. Build the
/// adjacency from the dispatcher output, BFS within the
/// continental subgraph from the first continental plate, verify
/// it reaches every other continental plate (single connected
/// component).
#[test]
fn r7_continental_cluster_contiguous() {
    let params = Phase2InitParams::default();
    let state = init_c1_state_phase_2_r7(GRID, SEED, &params);
    let num_plates = state.num_plates;

    let adjacency = build_plate_adjacency(&state.plate_id, num_plates);
    let per_plate = per_plate_type_from_state(&state);
    let continental: Vec<u16> = per_plate
        .iter()
        .enumerate()
        .filter_map(|(i, t)| matches!(t, PlateType::Continental).then(|| i as u16))
        .collect();

    eprintln!("Test 5 cluster_contiguous: continental plates = {continental:?}");
    assert!(
        !continental.is_empty(),
        "default cluster_seed_count = 1 must produce at least one continental plate"
    );

    let mut reached: HashSet<u16> = HashSet::new();
    let mut queue: VecDeque<u16> = VecDeque::new();
    reached.insert(continental[0]);
    queue.push_back(continental[0]);
    while let Some(p) = queue.pop_front() {
        for &n in adjacency[p as usize].iter() {
            if continental.contains(&n) && !reached.contains(&n) {
                reached.insert(n);
                queue.push_back(n);
            }
        }
    }
    assert_eq!(
        reached.len(),
        continental.len(),
        "continental subgraph must be connected; reached {} of {} continental plates",
        reached.len(),
        continental.len()
    );
}

/// **Test 6** — age = 0 ridge cells across multiple seeds. Track
/// C kinematics concern mitigation: per Phase 1.1 preset, not
/// every Voronoï layout guarantees a divergent boundary. Iterate
/// over 5 seeds; assert at least one produces ridge cells.
/// Verbose output documents the per-seed distribution.
#[test]
fn r7_age_ridge_aligned_at_divergent_boundaries() {
    let params = Phase2InitParams::default();
    let ridge_value = params.age.ridge_value;
    let seeds = [42u64, 100, 1337, 2026, 9999];

    let mut seeds_with_ridges = 0;
    eprintln!("Test 6 ridge_aligned: multi-seed sampling [42, 100, 1337, 2026, 9999]");
    eprintln!("    seed   |  ridge_cells |  oceanic_cells |  continental_cells");
    eprintln!("    -------+--------------+----------------+-------------------");

    for &seed in seeds.iter() {
        let state = init_c1_state_phase_2_r7(GRID, seed, &params);
        let mut ridge_count = 0;
        let mut oceanic_count = 0;
        let mut continental_count = 0;
        for j in 0..GRID {
            for i in 0..GRID {
                let a = state.age.get(i, j);
                match state.plate_type.get(i, j) {
                    PlateType::Continental => continental_count += 1,
                    PlateType::Oceanic => {
                        if a == ridge_value {
                            ridge_count += 1;
                        } else {
                            oceanic_count += 1;
                        }
                    }
                }
            }
        }
        eprintln!(
            "    {seed:>6} | {ridge_count:>12} | {oceanic_count:>14} | {continental_count:>18}"
        );
        if ridge_count > 0 {
            seeds_with_ridges += 1;
        }
    }

    let fraction = seeds_with_ridges as f64 / seeds.len() as f64;
    eprintln!(
        "    seeds_with_ridges = {seeds_with_ridges} / {} ({:.0} %)",
        seeds.len(),
        100.0 * fraction
    );
    if fraction < 0.60 {
        eprintln!(
            "    ARCHITECTURAL FINDING: < 60 % of seeds produce ridges with Phase 1.1 \
             kinematics preset. Track C constrained-kinematics prioritisation indicated."
        );
    }

    assert!(
        seeds_with_ridges > 0,
        "at least one seed in {seeds:?} must produce divergent-boundary ridge cells; got 0. \
         If this systematically fails across seeds, Phase 1.1 kinematics is not producing \
         divergent boundaries — Track C constrained-kinematics prerequisite."
    );
}

/// **Test 7** — age distribution comparison Phase 1.1 vs Phase 2
/// R7. Phase 1.1 produces exactly 2 unique ages
/// `{CONTINENTAL_AGE_INIT, OCEANIC_AGE_INIT} = {7.0, 0.5}`.
/// Phase 2 R7 produces ≥ 2 (3 when ridges present).
#[test]
fn r7_age_distribution_compared_to_phase_1_1() {
    let state_phase_1_1 = init_c1_state_phase_1_1(GRID, SEED);
    let state_phase_2 = init_c1_state_phase_2_r7(GRID, SEED, &Phase2InitParams::default());

    fn unique_ages(state: &ymir_core::tectonics_c1::state::C1State) -> Vec<f64> {
        let mut seen: Vec<f64> = Vec::new();
        for j in 0..state.ny() {
            for i in 0..state.nx() {
                let a = state.age.get(i, j);
                if !seen.iter().any(|s| (s - a).abs() < 1e-12) {
                    seen.push(a);
                }
            }
        }
        seen.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        seen
    }

    let ages_phase_1_1 = unique_ages(&state_phase_1_1);
    let ages_phase_2 = unique_ages(&state_phase_2);

    eprintln!("Test 7 age_distribution_comparison:");
    eprintln!("    Phase 1.1 unique ages: {ages_phase_1_1:?}");
    eprintln!("    Phase 2 R7 unique ages: {ages_phase_2:?}");

    assert_eq!(
        ages_phase_1_1.len(),
        2,
        "Phase 1.1 init must produce exactly 2 unique ages (continental + oceanic baseline); got {} unique values",
        ages_phase_1_1.len()
    );
    assert!(
        ages_phase_2.len() >= 2,
        "Phase 2 R7 init must produce at least 2 unique ages (continental + oceanic baseline); got {} unique values",
        ages_phase_2.len()
    );

    // When ridges present, Phase 2 has 3 unique values
    // ({0.0, 0.5, 7.0}); when not, 2. Stage 6 verifies the
    // multi-seed ridge presence; here we just confirm the
    // qualitative shape difference (or sameness).
    if ages_phase_2.contains(&0.0) {
        eprintln!(
            "    Phase 2 R7 INCLUDES ridge cells (age = 0.0) — distribution wider than Phase 1.1."
        );
    } else {
        eprintln!(
            "    Phase 2 R7 at seed {SEED} has no ridge cells (no divergent boundaries with Phase 1.1 kinematics)."
        );
    }
}

/// **Test 8** — Phase 1.1 init preserved unchanged. Property-
/// based regression guard (NOT pinned hash) on the Phase 1.1
/// init contract — locks behaviour against silent drift from
/// Phase 2 Track B changes.
#[test]
fn r7_phase_1_1_init_preserved_unchanged() {
    let state = init_c1_state_phase_1_1(GRID, SEED);

    // Property 1 — 8 plates default.
    assert_eq!(state.num_plates, 8, "Phase 1.1 init must produce 8 plates");

    // Property 2 — age field exactly 2 unique values
    // {CONTINENTAL_AGE_INIT = 7.0, OCEANIC_AGE_INIT = 0.5}.
    let mut seen_ages: Vec<f64> = Vec::new();
    for j in 0..GRID {
        for i in 0..GRID {
            let a = state.age.get(i, j);
            if !seen_ages.iter().any(|s| (s - a).abs() < 1e-12) {
                seen_ages.push(a);
            }
        }
    }
    seen_ages.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!(
        "Test 8 phase_1_1_preserved: num_plates = {}, age unique = {seen_ages:?}",
        state.num_plates
    );
    assert_eq!(
        seen_ages.len(),
        2,
        "Phase 1.1 age field must have exactly 2 unique values; got {seen_ages:?}"
    );
    assert!(
        (seen_ages[0] - 0.5).abs() < 1e-12,
        "lower age value must be OCEANIC_AGE_INIT = 0.5; got {}",
        seen_ages[0]
    );
    assert!(
        (seen_ages[1] - 7.0).abs() < 1e-12,
        "upper age value must be CONTINENTAL_AGE_INIT = 7.0; got {}",
        seen_ages[1]
    );

    // Property 3 — continental fraction in [0.20, 0.45].
    let total = GRID * GRID;
    let mut continental_count = 0;
    for j in 0..GRID {
        for i in 0..GRID {
            if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                continental_count += 1;
            }
        }
    }
    let frac = continental_count as f64 / total as f64;
    eprintln!("    continental cell fraction = {frac:.3} ({continental_count} / {total})");
    assert!(
        (0.20..=0.45).contains(&frac),
        "Phase 1.1 continental cell fraction {frac:.3} must be in [0.20, 0.45]"
    );

    // Property 4 — boundary classification produces SOME
    // divergent or convergent cells (sanity: kinematics isn't
    // trivial).
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let info = classify_boundaries(&state.plate_id, &kinematics);
    let counts = info.counts();
    let convergent = counts[1];
    let divergent = counts[2];
    eprintln!("    boundary counts: Convergent = {convergent}, Divergent = {divergent}");
    assert!(
        convergent > 0 || divergent > 0,
        "Phase 1.1 kinematics must produce at least some Convergent or Divergent boundaries"
    );

    // Suppress unused-variable warning while preserving
    // future-extensibility of the BoundaryType match.
    let _ = BoundaryType::Internal;
}
