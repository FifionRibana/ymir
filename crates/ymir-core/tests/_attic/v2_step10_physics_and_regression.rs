//! Step 10 Phase 6 — physics baseline + Step 9 regression.
//! Marked `#[ignore]`; run explicitly:
//!
//! ```text
//! cargo test --release -p ymir-core \
//!     --test v2_step10_physics_and_regression \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! Two ignored tests:
//! - `step10_physics_baseline_64sq` — Step 8 shape (mantle on,
//!   slab off) at 64² × 100 steps with `AgeFieldConfig::Enabled`
//!   defaults. Reports the §4.11 metrics: per-region age means,
//!   event counts, the bound check (acceptance #1 / #11). The
//!   final `A` field is exported alongside `S̃` for the visual
//!   checkpoint embedded in `step10_physics_report.md`.
//! - `step10_regression_disabled_64sq` — same shape with
//!   `AgeFieldConfig::Disabled`, anchors the regression
//!   comparison (acceptance #12, #13, #14).

use std::path::PathBuf;
use std::time::Instant;

use ymir_core::tectonics_v2::age_field::{AgeFieldConfig, AgeFieldConfigEnabled};
use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::CratonicConfig;
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, BaselineResult, ForceKind, NonlinearChoice, build_force, run_baseline,
};
use ymir_core::tectonics_v2::mantle::{
    COUPLING_DEFAULT, MF_DEFAULT, MantleConfig, NUM_MODES_DEFAULT,
};
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
const STEPS: usize = 100;
const SEED: u64 = 42;

/// Step 8 shape (mantle on, slab off — same as Step 8 baseline +
/// Step 9 immunity demonstration regime). Step 10 baseline runs
/// the same shape because:
/// - `AgeFieldConfig` is passive; it does not interact with the
///   Stokes operator. Choosing the most-active dynamic regime
///   exposes the age field to the most ridge / arc / collision
///   events.
/// - Step 7 shape (no mantle) is too quiescent to generate
///   meaningful event counts; the age field would barely deviate
///   from `init_age + dt · steps`.
fn build_step10_config(age_field: AgeFieldConfig, label: &str) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").unwrap();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates = BoundaryRates {
        k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0,
    };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX, NY, &vcfg, SEED, rates, RecyclingConfig::default(),
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
        total_time_nondim: 6.0,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: vec![0.0, 1.0],
        output_dir: PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../docs/reports/step10_phase6_{}", label)),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw {
            bi: 0.15, ..Default::default()
        }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: 0.05, ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n8_step10_{}", SEED, label),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Enabled {
            mf: MF_DEFAULT,
            coupling: COUPLING_DEFAULT,
            num_modes: NUM_MODES_DEFAULT,
            seed: 7,
            evolution_rate: 0.0,
        },
        cratonic: CratonicConfig::Disabled,
        age_field,
        capture: None,
        linear_solver: LinearSolverConfig::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
            plate_kinematic: ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig::Zero,
    }
}

fn print_summary(label: &str, dt: f64, r: &BaselineResult) {
    let m = &r.metrics;
    let na = m.newton.as_ref().expect("newton aggregate");
    println!();
    println!("=== Step 10 {} (64x64 Step 8 shape) ===", label);
    println!("  wallclock                       : {:.2} s", dt);
    println!("  CG iters mean                   : {:.1}", m.cg_iter_mean);
    println!("  Newton outer mean               : {:.2}", na.outer_iters_mean());
    println!("  peak|v|                         : {:.3e}", m.vmax_peak);
    println!(
        "  yielding_cell_fraction_max      : {:.4}",
        na.yielding_cell_fraction_max.unwrap_or(0.0)
    );
    if let Some(ci) = na.continental_age_init_diagnostic {
        println!("  --- Step 10 age-field metrics ---");
        println!("  continental_age_init (config)   : {}", ci);
        println!(
            "  oceanic_age_init (config)       : {}",
            na.oceanic_age_init_diagnostic.unwrap_or(0.0)
        );
        println!(
            "  age_field [min, max, mean]      : [{:.4}, {:.4}, {:.4}]",
            na.age_field_min_final.unwrap_or(0.0),
            na.age_field_max_final.unwrap_or(0.0),
            na.age_field_mean_final.unwrap_or(0.0)
        );
        println!(
            "  age_at_continental_cells_mean   : {:.4}",
            na.age_at_continental_cells_mean_final.unwrap_or(0.0)
        );
        println!(
            "  age_at_oceanic_cells_mean       : {:.4}",
            na.age_at_oceanic_cells_mean_final.unwrap_or(0.0)
        );
        println!(
            "  ridge_resets total              : {}",
            na.age_ridge_resets_total.unwrap_or(0)
        );
        println!(
            "  arc_resets total                : {}",
            na.age_arc_resets_total.unwrap_or(0)
        );
        println!(
            "  collision_max_events total      : {}",
            na.age_collision_max_events_total.unwrap_or(0)
        );
        println!(
            "  collision_max_age mean          : {:.4}",
            na.age_collision_max_age_mean.unwrap_or(0.0)
        );
    }
}

#[test]
#[ignore]
fn step10_physics_baseline_64sq() {
    let cfg_age = AgeFieldConfig::Enabled(AgeFieldConfigEnabled::default());
    let cfg = build_step10_config(cfg_age, "baseline");
    let t0 = Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    print_summary("Enabled (defaults)", dt, &r);

    // Acceptance #1 / #11: A bounded.
    let na = r.metrics.newton.as_ref().expect("newton aggregate");
    let amin = na.age_field_min_final.expect("min populated");
    let amax = na.age_field_max_final.expect("max populated");
    let init_max = na
        .continental_age_init_diagnostic
        .unwrap_or(0.0)
        .max(na.oceanic_age_init_diagnostic.unwrap_or(0.0));
    let total_time = 6.0_f64; // total_time_nondim
    assert!(amin >= 0.0, "age_min = {} < 0", amin);
    assert!(
        amax <= init_max + total_time + 1e-6,
        "age_max = {} exceeds init_max + total_time = {}",
        amax,
        init_max + total_time
    );

    // Acceptance #10: mass conservation preserved.
    let mass_res = na.mass_conservation_residual.expect("mass residual populated");
    assert!(
        mass_res < 1e-6,
        "mass_conservation_residual = {} ≥ 1e-6",
        mass_res
    );

    // Acceptance #8 (soft): mean oceanic age generally smaller
    // than continental — ridge resets fire frequently on oceanic
    // boundary cells while continental cells are reset only on
    // arc/collision (rarer). This is informational only; we
    // print the values rather than asserting hard.
    let cont_mean = na.age_at_continental_cells_mean_final.unwrap_or(0.0);
    let oce_mean = na.age_at_oceanic_cells_mean_final.unwrap_or(0.0);
    println!(
        "  acceptance #8 informational: continental_mean = {:.3}, oceanic_mean = {:.3}",
        cont_mean, oce_mean
    );
}

#[test]
#[ignore]
fn step10_regression_disabled_64sq() {
    // Anchor for the regression report. With AgeFieldConfig::Disabled,
    // every Step 10 metric is None and the run reproduces the Step
    // 8-shape (mantle on) numerics from the Step 9 baseline.
    let cfg = build_step10_config(AgeFieldConfig::Disabled, "disabled_reference");
    let t0 = Instant::now();
    let r = run_baseline(&cfg);
    let dt = t0.elapsed().as_secs_f64();
    print_summary("Disabled (regression anchor)", dt, &r);

    let na = r.metrics.newton.as_ref().expect("newton aggregate");
    assert!(na.continental_age_init_diagnostic.is_none());
    assert!(na.age_field_min_final.is_none());
    assert!(na.age_field_max_final.is_none());
    assert!(na.age_ridge_resets_total.is_none());
    assert!(na.age_arc_resets_total.is_none());
    assert!(na.age_collision_max_events_total.is_none());
    assert!(na.age_collision_max_age_mean.is_none());
}

#[test]
#[ignore]
fn step10_disabled_runs_are_bit_deterministic() {
    // Acceptance #12: AgeFieldConfig::Disabled bit-deterministic.
    // Two identical runs must produce identical metrics.
    let r1 = run_baseline(&build_step10_config(AgeFieldConfig::Disabled, "det_a"));
    let r2 = run_baseline(&build_step10_config(AgeFieldConfig::Disabled, "det_b"));
    let na1 = r1.metrics.newton.as_ref().unwrap();
    let na2 = r2.metrics.newton.as_ref().unwrap();
    assert_eq!(
        na1.mass_conservation_residual.unwrap(),
        na2.mass_conservation_residual.unwrap(),
        "mass_conservation_residual not deterministic — Step 10 code path has a side-effect"
    );
    assert_eq!(r1.metrics.vmax_peak, r2.metrics.vmax_peak);
    assert_eq!(r1.metrics.cg_iter_mean, r2.metrics.cg_iter_mean);
}
