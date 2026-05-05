//! Step 13 Phase 7 — solver-health acceptance #10.
//!
//! Acceptance #10: "CG iters ratio ≤ 1.1× existing modes baseline.
//! Small overhead acceptable for richer init, but no major
//! degradation."
//!
//! Method: build the same `BaselineConfig` (Step 8 shape — mantle
//! on, slab off, single_continent-like Voronoï layout) three times,
//! varying only `init_mode` between `Uniform` (baseline),
//! `RadialProfile { Smoothstep, defaults }`, and
//! `RadialProfileWithFBM { defaults }`. Run each via
//! [`run_baseline`] and read `metrics.cg_iter_mean`. Compare each
//! new mode's mean to the `Uniform` baseline.
//!
//! Heavy test — `#[ignore]` so it only runs when explicitly asked:
//!
//! ```text
//! cargo test --release -p ymir-core --test v2_step13_cg_ratio \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! Reduced 64² × 30 steps to keep wallclock < 60 s per run; the
//! ratio is a per-Newton-step average so the smaller step count
//! is statistically adequate.

use std::path::PathBuf;
use std::time::Instant;

use ymir_core::tectonics_v2::age_field::AgeFieldConfig;
use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::CratonicConfig;
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, BaselineResult, ForceKind, NonlinearChoice, build_force, run_baseline,
};
use ymir_core::tectonics_v2::init::{
    FBM_AMPLITUDE_DEFAULT, FBM_LACUNARITY_DEFAULT, FBM_OCTAVES_DEFAULT, FBM_PERSISTENCE_DEFAULT,
    FBM_SCALE_DEFAULT, FBM_SEED_DEFAULT, InitMode, ProfileShape,
};
use ymir_core::tectonics_v2::mantle::{
    COUPLING_DEFAULT, MF_DEFAULT, MantleConfig, NUM_MODES_DEFAULT,
};
use ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

const NX: usize = 64;
const NY: usize = 64;
const STEPS: usize = 30;
const SEED: u64 = 12;

fn build_cfg(init_mode: InitMode, label: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    // single_continent Voronoï shape — same as the
    // v2_step13_acceptance suite for layout consistency.
    let vcfg = VoronoiConfig {
        num_plates: 4,
        continental_ratio: 0.5,
    };
    let rates = BoundaryRates {
        k_sub: 0.5,
        k_arc: 0.0,
        k_spread: 0.0,
        k_coll_v: 0.0,
        k_rift_v: 0.0,
    };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX,
        NY,
        &vcfg,
        SEED,
        rates,
        RecyclingConfig::default(),
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: NX,
        grid_ny: NY,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 1.5,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: vec![],
        output_dir: PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../target/v2_step13_cg_ratio/{}", label)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.0,
        yielding: YieldingConfig::Enabled(YieldingLaw {
            bi: 0.15,
            ..Default::default()
        }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: 0.05,
            ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: format!("v2_step13_cg_ratio_{}", label),
        slab_pull: SlabPullConfig::Disabled,
        // Mantle on — exercises the dynamic regime (matches Step
        // 8/10 active baselines). Acceptance #10 wants headroom,
        // so we measure under load.
        mantle: MantleConfig::Enabled {
            mf: MF_DEFAULT,
            coupling: COUPLING_DEFAULT,
            num_modes: NUM_MODES_DEFAULT,
            seed: 7,
            evolution_rate: 0.0,
        },
        cratonic: CratonicConfig::Disabled,
        age_field: AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: LinearSolverConfig::default(),
        init_mode,
        continuation: None,
        plate_kinematic: PlateKinematicConfig::Zero,
    }
}

fn cg_iter_mean(label: &str, init_mode: InitMode) -> (f64, f64) {
    let cfg = build_cfg(init_mode, label);
    let t0 = Instant::now();
    let r: BaselineResult = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    let m = &r.metrics;
    let mean = m.cg_iter_mean;
    println!(
        "  {:<28} : cg_iter_mean = {:.2} (wallclock {:.2}s)",
        label, mean, dt
    );
    (mean, dt)
}

#[test]
#[ignore]
fn cg_ratio_acceptance() {
    println!("Step 13 Phase 7 — acceptance #10 (CG ratio ≤ 1.1× Uniform baseline)");
    println!(
        "  Config : 64² × {} steps, mantle on, slab off, yielding on, single_continent Voronoï (seed=12)",
        STEPS
    );
    println!();

    let (uniform_mean, _) = cg_iter_mean(
        "Uniform (baseline)",
        InitMode::Uniform { boundary_smoothing_width: 1.0 },
    );

    let (radial_mean, _) = cg_iter_mean(
        "RadialProfile",
        InitMode::RadialProfile {
            continental_value: 0.95,
            oceanic_value: 0.20,
            profile_shape: ProfileShape::Smoothstep,
        },
    );

    let (radial_fbm_mean, _) = cg_iter_mean(
        "RadialProfileWithFBM",
        InitMode::RadialProfileWithFBM {
            continental_value: 0.95,
            oceanic_value: 0.20,
            profile_shape: ProfileShape::Smoothstep,
            fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
            fbm_octaves: FBM_OCTAVES_DEFAULT,
            fbm_persistence: FBM_PERSISTENCE_DEFAULT,
            fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
            fbm_scale: FBM_SCALE_DEFAULT,
            fbm_seed: FBM_SEED_DEFAULT,
        },
    );

    let radial_ratio = radial_mean / uniform_mean.max(1e-12);
    let fbm_ratio = radial_fbm_mean / uniform_mean.max(1e-12);

    println!();
    println!(
        "  RadialProfile         / Uniform = {:.3}× (acceptance: ≤ 1.10)",
        radial_ratio
    );
    println!(
        "  RadialProfileWithFBM  / Uniform = {:.3}× (acceptance: ≤ 1.10)",
        fbm_ratio
    );

    const RATIO_LIMIT: f64 = 1.10;
    assert!(
        radial_ratio <= RATIO_LIMIT,
        "RadialProfile CG ratio {:.3}× exceeds 1.10× Uniform baseline",
        radial_ratio
    );
    assert!(
        fbm_ratio <= RATIO_LIMIT,
        "RadialProfileWithFBM CG ratio {:.3}× exceeds 1.10× Uniform baseline",
        fbm_ratio
    );
}
