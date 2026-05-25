//! Issue #127 Phase 1.4 Stage I2 — isostasy E2E validation tests.
//!
//! Three tests under default features (the C1 path is mainline;
//! v2's `compute_isostasy` is itself default-features-on since
//! `tectonics::isostasy` predates the Phase 1.3 H2 refactor).
//!
//! Test inventory (per H1 architectural analysis — see commit
//! Stage I2 surface for the apply_post_tectonic structural
//! finding):
//!
//! - [`altitude_responds_to_s_changes_per_step`] — locks the
//!   property that the per-step erosion pipeline's altitude
//!   (stage 4a in [`run_with_closures`]) genuinely evolves
//!   across the run as `S̃` mutates under the three active
//!   closures. Samples multiple convergent-boundary cells to
//!   avoid pathological single-cell steady-state.
//! - [`compute_isostasy_deterministic_given_same_s`] — locks
//!   that `compute_isostasy` is a pure function: two calls on
//!   the same `S̃` field with the same config produce
//!   bit-identical `heightmap.data`. Used downstream by Stages
//!   E4 / D to safely re-compute altitude without caching.
//! - [`apply_post_tectonic_mutates_s_via_macro_redistribution`]
//!   — locks the architectural property that
//!   `apply_post_tectonic` Step 2 (`macro_redistribution::apply`)
//!   mutates `s_field` in place. **Per-step altitude (stage 4a,
//!   computed on post-erosion `S̃`) is therefore ≠ end-of-cycle
//!   altitude (computed on post-macro `S̃`) by construction.**
//!   This is documented in `time_loop.rs` § "End-of-cycle
//!   apply_post_tectonic consistency" as a deliberate design
//!   choice (clean per-step vs per-cycle concern separation
//!   over sharing computed altitude across the boundary). The
//!   test name self-documents the finding so future Phase 1.4
//!   reviewers / Phase 2+ developers see it as design, not bug.
//!
//! These are H1/I2 validation tests, **not** Phase 1.4
//! acceptance tests — those land in Stage E4 (4 invariants on
//! erosion behaviour at 64²×300).

use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::boundary_classification::{classify_boundaries, BoundaryType};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::workflow::phase_a_common::{
    apply_post_tectonic, extract_per_plate_type, PostTectonicInput,
};
use ymir_core::tectonics_v2::workflow::{PhaseAParams, WorkflowParams};

const GRID: usize = 64;
const SEED: u64 = 42;

/// Phase 1.4 closure stack — DS + EQ + erosion ON, S-S OFF. See
/// `c1_phase_1_4_erosion.rs::phase_1_4_closures` for the rationale
/// on holding the Phase 1.4 regime stable post-#129.
fn phase_1_4_closures() -> C1Closures {
    C1Closures {
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        ..C1Closures::default()
    }
}

#[test]
fn altitude_responds_to_s_changes_per_step() {
    // After 50 steps with all 3 Phase 1.4 closures enabled, the
    // altitude at convergent-boundary cells must differ from the
    // initial altitude by more than the test's tolerance. The
    // direction is **not** asserted: it could be up (Davis-Suppe
    // dominates), down (erosion dominates), or capped at h_eq
    // (equilibrium dominates) — all three are physically
    // legitimate outcomes of the joint dynamics. We only assert
    // that *something* moves.
    //
    // Multi-cell sampling (max_delta over ALL convergent cells)
    // protects against the edge case where a single sampled cell
    // happens to land at steady-state.
    let mut state = init_c1_state_phase_1_1(GRID, SEED);
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();

    // Snapshot altitude at cycle 0.
    let altitude_0 = compute_isostasy(&state.s, &iso_config);

    // Sample convergent boundary cells from the (static) plate
    // tessellation. classify_boundaries reads plate_id +
    // kinematics — both static across the run, so the sampled
    // set is invariant.
    let boundary = classify_boundaries(&state.plate_id, &kinematics);
    let mut convergent_cells: Vec<(usize, usize)> = Vec::new();
    for j in 0..state.ny() {
        for i in 0..state.nx() {
            if matches!(boundary.boundary_type.get(i, j), BoundaryType::Convergent) {
                convergent_cells.push((i, j));
            }
        }
    }
    assert!(
        !convergent_cells.is_empty(),
        "test premise: Phase 1.1 init must produce some convergent boundary cells"
    );

    // Run 50 steps with full Phase 1.4 closure stack.
    let closures = phase_1_4_closures();
    let config = C1TimeLoopConfig {
        n_steps: 50,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &kinematics, &config, &closures, |_, _| {});

    // Snapshot altitude at cycle 50.
    let altitude_50 = compute_isostasy(&state.s, &iso_config);

    // Compute max altitude delta over all convergent cells.
    let mut max_delta = 0.0_f32;
    for &(i, j) in &convergent_cells {
        let h0 = altitude_0.heightmap.get(i as i32, j as i32);
        let h50 = altitude_50.heightmap.get(i as i32, j as i32);
        let delta = (h50 - h0).abs();
        if delta > max_delta {
            max_delta = delta;
        }
    }
    eprintln!(
        "c1_phase_1_4 I2-T1: max altitude delta over {} convergent cells = {:.4} after 50 steps",
        convergent_cells.len(),
        max_delta
    );
    assert!(
        max_delta > 0.01,
        "altitude must evolve across 50 steps with closures active; max delta = {max_delta:.4} ≤ 0.01"
    );
}

#[test]
fn compute_isostasy_deterministic_given_same_s() {
    // `compute_isostasy` is a pure function — two calls on the
    // same `S̃` field with the same `IsostasyConfig` must produce
    // bit-identical `heightmap.data` (and all other output fields).
    //
    // Locks: no internal cache, no RNG, no global state. Used
    // downstream by Stages E4 / D / Phase 1.4 K calibration to
    // re-compute altitude at multiple points in the run without
    // worrying about staleness.
    let state = init_c1_state_phase_1_1(GRID, SEED);
    // Mutate `S̃` slightly so we're not testing on a trivially-
    // uniform field where Gaussian-blur identity is degenerate.
    let mut s = state.s.clone();
    s.set(10, 10, 1.5);
    s.set(20, 20, 0.05);
    let iso_config = IsostasyConfig::default();

    let isostasy_a = compute_isostasy(&s, &iso_config);
    let isostasy_b = compute_isostasy(&s, &iso_config);

    assert_eq!(
        isostasy_a.heightmap.data, isostasy_b.heightmap.data,
        "compute_isostasy must be deterministic — heightmap.data bit-identical"
    );
    assert_eq!(
        isostasy_a.heightmap.width, isostasy_b.heightmap.width,
        "heightmap dimensions must match"
    );
    assert_eq!(
        isostasy_a.heightmap.height, isostasy_b.heightmap.height,
        "heightmap dimensions must match"
    );
    assert_eq!(
        isostasy_a.sea_level_normalized, isostasy_b.sea_level_normalized,
        "sea_level_normalized must be bit-identical"
    );
    assert_eq!(
        isostasy_a.peak_altitude_m, isostasy_b.peak_altitude_m,
        "peak_altitude_m must be bit-identical"
    );
    assert_eq!(
        isostasy_a.max_depth_m, isostasy_b.max_depth_m,
        "max_depth_m must be bit-identical"
    );
    assert_eq!(
        isostasy_a.land_ratio, isostasy_b.land_ratio,
        "land_ratio must be bit-identical"
    );
}

#[test]
fn apply_post_tectonic_mutates_s_via_macro_redistribution() {
    // ARCHITECTURAL LOCK — DOCUMENTS THE DESIGN, NOT A BUG.
    //
    // `apply_post_tectonic` runs four steps:
    //   1. compute_sea_level_ref_s_space — reads s
    //   2. macro_redistribution::apply — **MUTATES s in place**
    //   3. reclassify_inplace — reads s, mutates plate_type
    //   4. (optional) recompute_cratonic_factor_for_cycle —
    //      reads s, produces new Field2D
    //
    // Step 2's mutation means: the per-step altitude (computed at
    // `run_with_closures` stage 4a, on post-erosion `S̃`) is
    // **structurally different** from end-of-cycle altitude
    // (computed externally on `S̃` post-everything, including
    // post-macro). This test asserts that divergence is
    // measurable — proving macro_redistribution actually
    // mutated `s`.
    //
    // If this test fails, EITHER:
    //   - macro_redistribution stopped mutating s (regression in
    //     `workflow::macro_redistribution::apply`), OR
    //   - macro_redistribution somehow happens to be a no-op on
    //     the Phase 1.1 init state (would invalidate the user's
    //     assumption in `time_loop.rs` § "End-of-cycle
    //     apply_post_tectonic consistency").
    //
    // Either way, the test surface forces a re-evaluation of the
    // documented design rather than letting the divergence go
    // silent.

    // Build a state with some Phase 1.4 dynamics — 30 steps of
    // run_with_closures so macro_redistribution has non-trivial
    // relief to operate on (the Phase 1.1 init alone is too
    // uniform to make the mutation reliably visible at the
    // 1e-6 tolerance).
    let mut state = init_c1_state_phase_1_1(GRID, SEED);
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    let closures = phase_1_4_closures();
    let config = C1TimeLoopConfig {
        n_steps: 30,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &kinematics, &config, &closures, |_, _| {});

    // Snapshot s + altitude *before* apply_post_tectonic.
    let s_pre: Vec<f64> = state.s.data().to_vec();
    let altitude_pre = compute_isostasy(&state.s, &iso_config);

    // Capture the pre-reclassification per-plate type for the D4
    // gate inside apply_post_tectonic.
    let original_per_plate_type = extract_per_plate_type(&state.plate_id, &state.plate_type);

    // Call apply_post_tectonic with the Enabled workflow params
    // — this is the path that runs macro_redistribution.
    let wf_params = WorkflowParams::default();
    let phase_a_params: PhaseAParams = wf_params.phase_a.clone();
    let _output = apply_post_tectonic(PostTectonicInput {
        s_field: &mut state.s,
        plate_id: Some(&state.plate_id),
        plate_type: Some(&mut state.plate_type),
        previous_cratonic_factor: None,
        initial_per_plate_type: Some(&original_per_plate_type),
        params: &phase_a_params,
        iso_cfg: &iso_config,
        cratonic_cfg: None,
    });

    // Snapshot s + altitude *after* apply_post_tectonic.
    let s_post: &[f64] = state.s.data();
    let altitude_post = compute_isostasy(&state.s, &iso_config);

    // Count mutated cells + compute max altitude delta — produces
    // the verbose evidence trail the user spec asked for.
    let mut mutated_cells = 0_usize;
    let mut max_s_delta = 0.0_f64;
    for (a, b) in s_pre.iter().zip(s_post.iter()) {
        let d = (a - b).abs();
        if d > 1e-12 {
            mutated_cells += 1;
        }
        if d > max_s_delta {
            max_s_delta = d;
        }
    }
    let mut max_altitude_delta = 0.0_f32;
    for k in 0..altitude_pre.heightmap.data.len() {
        let d = (altitude_pre.heightmap.data[k] - altitude_post.heightmap.data[k]).abs();
        if d > max_altitude_delta {
            max_altitude_delta = d;
        }
    }

    let total_cells = state.nx() * state.ny();
    eprintln!(
        "c1_phase_1_4 I2-T3 architectural lock — apply_post_tectonic mutates s via macro:"
    );
    eprintln!(
        "  S̃ cells mutated     = {mutated_cells} / {total_cells} ({:.1} %)",
        100.0 * mutated_cells as f64 / total_cells as f64
    );
    eprintln!("  S̃ max delta         = {max_s_delta:.4e}");
    eprintln!("  altitude max delta  = {max_altitude_delta:.4e}");
    eprintln!(
        "  → per-step altitude (stage 4a, on post-erosion S̃) differs from end-of-cycle"
    );
    eprintln!(
        "    altitude (post-everything S̃) by max {max_altitude_delta:.4e} due to"
    );
    eprintln!(
        "    macro_redistribution mutation. This is design, not bug — see"
    );
    eprintln!(
        "    `time_loop.rs` § \"End-of-cycle apply_post_tectonic consistency\"."
    );

    // Assertions — proves the architectural property.
    assert!(
        mutated_cells > 0,
        "apply_post_tectonic must mutate s via macro_redistribution; got 0 mutated cells \
         — macro_redistribution may have regressed to a no-op, or the test premise (relief \
         from 30 steps of closures) failed to generate non-trivial input"
    );
    assert!(
        max_s_delta > 1e-9,
        "S̃ mutation magnitude {max_s_delta:.4e} below 1e-9 — macro_redistribution effectively \
         absent"
    );
    assert!(
        max_altitude_delta > 1e-9,
        "altitude mutation magnitude {max_altitude_delta:.4e} below 1e-9 — Airy projection of \
         the S̃ mutation produced no visible altitude effect"
    );
}
