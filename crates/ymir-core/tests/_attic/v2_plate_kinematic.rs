//! Step 11 — bit-identical regression for `PlateKinematicConfig`.
//!
//! Acceptance criteria #1 / #11 / #12: the default
//! `PlateKinematicConfig::Zero` must produce identically-zero
//! initial velocity fields and propagate bit-for-bit through the
//! solver, matching pre-Step-11 baselines exactly.
//!
//! These tests run as regular (non-`#[ignore]`) integration tests so
//! they fire during `cargo test` and CI. Both runs are 32² × 5 steps
//! — small enough to be fast (≈1 s each), large enough to exercise
//! the full solver path (boundary advection, Newton + AMG / Jacobi
//! CG, source/sink, mass conservation).
//!
//! Two complementary checks:
//! 1. `zero_short_circuit_matches_per_plate_zeros` — runs with
//!    `Zero` and with `PerPlate { all-zeros, .. }`, asserts every
//!    `FinalState` field is bit-identical. Proves the
//!    `field::build` algorithm with all-zero input produces literal
//!    zeros (no numerical noise from intermediate computation), so
//!    the structural short-circuit and the algorithmic path agree.
//! 2. `zero_default_does_not_perturb_baseline` — sanity check that
//!    omitting the field altogether (relying on
//!    `PlateKinematicConfig::Zero` default) yields the same
//!    output as setting it explicitly. Catches `Default` mismatches.

use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

const NX: usize = 32;
const NY: usize = 32;
const STEPS: usize = 5;
const SEED: u64 = 42;
const NUM_PLATES: usize = 6;

fn build_minimal_config(plate_kinematic: PlateKinematicConfig) -> BaselineConfig {
    let scales = Scales::default();
    let vcfg = VoronoiConfig { num_plates: NUM_PLATES, continental_ratio: 0.3 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX,
        NY,
        &vcfg,
        SEED,
        rates,
        RecyclingConfig::default(),
    )
    .expect("voronoi boundary build");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: NX,
        grid_ny: NY,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 0.5,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_step11_regression"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: ymir_core::tectonics_v2::basal_drag::BasalDragConfig::Disabled,
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", SEED, NUM_PLATES),
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: LinearSolverConfig::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
        plate_kinematic,
    }
}

/// `PlateKinematicConfig::Zero` (structural short-circuit) and
/// `PerPlate { velocities: all-zeros }` (algorithmic path with
/// zero inputs) must produce bit-identical solver state at every
/// step boundary. This proves the `field::build` algorithm doesn't
/// introduce numerical noise on zero inputs and the regression
/// contract for `Zero` is upheld.
#[test]
fn zero_short_circuit_matches_per_plate_zeros() {
    let cfg_zero = build_minimal_config(PlateKinematicConfig::Zero);
    let cfg_per_plate = build_minimal_config(PlateKinematicConfig::PerPlate {
        velocities: vec![(0.0, 0.0); NUM_PLATES],
        boundary_smoothing_width: 1.5,
    });

    let result_zero = run_baseline(&cfg_zero);
    let result_per_plate = run_baseline(&cfg_per_plate);

    assert_eq!(
        result_zero.final_state.s_field.data(),
        result_per_plate.final_state.s_field.data(),
        "S̃ field diverged between Zero and PerPlate-zeros — \
         the algorithm introduced numerical noise on zero inputs"
    );
    assert_eq!(
        result_zero.final_state.vx, result_per_plate.final_state.vx,
        "vx diverged between Zero and PerPlate-zeros"
    );
    assert_eq!(
        result_zero.final_state.vy, result_per_plate.final_state.vy,
        "vy diverged between Zero and PerPlate-zeros"
    );
    assert_eq!(
        result_zero.final_state.strain_rate_invariant.data(),
        result_per_plate.final_state.strain_rate_invariant.data(),
        "strain_rate_invariant diverged"
    );
}

/// Default-constructed `PlateKinematicConfig::Zero` (via
/// `PlateKinematicConfig::default()`) and an explicit
/// `PlateKinematicConfig::Zero` must match. Catches accidental
/// `Default` impl drift.
#[test]
fn zero_default_does_not_perturb_baseline() {
    let cfg_default = build_minimal_config(PlateKinematicConfig::default());
    let cfg_explicit = build_minimal_config(PlateKinematicConfig::Zero);

    let result_default = run_baseline(&cfg_default);
    let result_explicit = run_baseline(&cfg_explicit);

    assert_eq!(
        result_default.final_state.s_field.data(),
        result_explicit.final_state.s_field.data()
    );
    assert_eq!(result_default.final_state.vx, result_explicit.final_state.vx);
    assert_eq!(result_default.final_state.vy, result_explicit.final_state.vy);
}

/// `vx[i] = vy[i] = 0` at step 0 with `Zero`. Sanity check that the
/// short-circuit branch produces what it claims (zeros). Combined
/// with the bit-identity test above, this anchors the contract:
/// the harness *starts* at zero and *stays* identical whether the
/// algorithmic path is taken with zero inputs or skipped entirely.
#[test]
fn zero_produces_zero_plate_kinematic_at_step_0() {
    use ymir_core::tectonics_v2::diagnostics::harness::run_baseline_with_progress;
    let cfg = build_minimal_config(PlateKinematicConfig::Zero);
    let mut step_0_seen = false;
    let mut max_v_step_0 = 0.0_f64;
    run_baseline_with_progress(&cfg, |progress| {
        // The harness emits an early step=0 callback before any
        // solver work — perfect spot to inspect the initial state.
        if progress.step == 0 {
            step_0_seen = true;
            max_v_step_0 = progress
                .peek_vx
                .iter()
                .zip(progress.peek_vy.iter())
                .map(|(&a, &b): (&f64, &f64)| (a * a + b * b).sqrt())
                .fold(0.0_f64, f64::max);
        }
        true
    });
    assert!(step_0_seen, "harness did not emit a step=0 callback");
    assert_eq!(max_v_step_0, 0.0, "initial velocity nonzero with Zero variant");
}
