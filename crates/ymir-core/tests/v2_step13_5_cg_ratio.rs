//! Step 13.5 Phase 6 — solver-health acceptance #9.
//!
//! Acceptance #9: "Solver runs with `apply_fbm_to_oceanic = true`
//! produce CG iters within ±10% of the disabled-flag baseline.
//! Acceptance `cg_ratio ∈ [0.90, 1.10]`."
//!
//! Method: build the same `BaselineConfig` (Step 13.5 shape —
//! mantle on, slab off, single_continent Voronoï layout) twice,
//! varying only the `apply_fbm_to_oceanic` flag inside the
//! `RadialProfileWithFBM` `init_mode`. Run each via
//! [`run_baseline`] and read `metrics.cg_iter_mean`. Compare the
//! enabled-flag mean to the disabled-flag baseline.
//!
//! Heavy test — `#[ignore]` so it only runs when explicitly asked:
//!
//! ```text
//! cargo test --release -p ymir-core --test v2_step13_5_cg_ratio \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! Reduced 64² × 30 steps (same shape as Step 13's CG ratio test)
//! to keep wallclock < 6 min per run; the ratio is a per-Newton-
//! step average so the smaller step count is statistically
//! adequate.

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
    FBM_AMPLITUDE_DEFAULT, FBM_AMPLITUDE_OCEANIC_DEFAULT, FBM_LACUNARITY_DEFAULT,
    FBM_OCTAVES_DEFAULT, FBM_PERSISTENCE_DEFAULT, FBM_SCALE_DEFAULT, FBM_SEED_DEFAULT, InitMode,
    ProfileShape,
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

fn radial_fbm_mode(apply_fbm_to_oceanic: bool) -> InitMode {
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
        apply_fbm_to_oceanic,
        fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
        fbm_scale_oceanic: None,
        fbm_seed_oceanic: None,
    }
}

fn build_cfg(init_mode: InitMode, label: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    // single_continent Voronoï shape — same as Step 13's CG ratio
    // suite for layout consistency.
    let vcfg = VoronoiConfig { num_plates: 4, continental_ratio: 0.5 };
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
            .join(format!("../../target/v2_step13_5_cg_ratio/{}", label)),
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
        boundary_layout_name: format!("v2_step13_5_cg_ratio_{}", label),
        slab_pull: SlabPullConfig::Disabled,
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
        "  {:<32} : cg_iter_mean = {:.2} (wallclock {:.2}s)",
        label, mean, dt
    );
    (mean, dt)
}

#[test]
#[ignore]
fn oceanic_fbm_cg_ratio_acceptance() {
    println!(
        "Step 13.5 Phase 6 — acceptance #9 (CG ratio ∈ [0.90, 1.10] vs disabled-flag baseline)"
    );
    println!(
        "  Config : 64² × {} steps, mantle on, slab off, yielding on, single_continent Voronoï (seed=12)",
        STEPS
    );
    println!(
        "  Init   : RadialProfileWithFBM with continental FBM enabled in both runs;"
    );
    println!(
        "           only `apply_fbm_to_oceanic` differs (false vs true)."
    );
    println!();

    let (baseline_mean, _) =
        cg_iter_mean("oceanic_disabled (Step 13)", radial_fbm_mode(false));

    let (oceanic_mean, _) =
        cg_iter_mean("oceanic_enabled (Step 13.5)", radial_fbm_mode(true));

    let ratio = oceanic_mean / baseline_mean.max(1e-12);

    println!();
    println!(
        "  oceanic_enabled / oceanic_disabled = {:.3}× (acceptance: ∈ [0.90, 1.10])",
        ratio
    );

    const RATIO_LIMIT_LOW: f64 = 0.90;
    const RATIO_LIMIT_HIGH: f64 = 1.10;
    assert!(
        ratio >= RATIO_LIMIT_LOW,
        "oceanic FBM CG ratio {:.3}× below 0.90× — unexpected speed-up; \
         diagnose before accepting (run was probably under-resolved or the \
         baseline was already very stiff)",
        ratio
    );
    assert!(
        ratio <= RATIO_LIMIT_HIGH,
        "oceanic FBM CG ratio {:.3}× above 1.10× — solver health degrades \
         when the oceanic FBM is enabled. Diagnose: maybe the bathymetric \
         heterogeneity stresses the Jacobi preconditioner past its working \
         range.",
        ratio
    );
}
