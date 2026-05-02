//! Step 8.5b 32² baseline measurement — Step 9 readiness check.
//!
//! 4 representative physics cases × 2 preconditioners × 5 runs at
//! 32² × 100 steps. Captures wallclock (mean ± std), CG iter
//! count, Newton outer iter count, Newton convergence rate, and
//! the Phase 5 extrapolation fallback rate. Output is printed to
//! stderr in a markdown-friendly format that the
//! `step8_5b_baseline_32sq.md` report ingests directly.
//!
//! `#[ignore]` because the run takes ~15-30 minutes total. Invoke
//! explicitly:
//!
//! ```bash
//! cargo test --release --test v2_baseline_32sq -- --ignored --nocapture
//! ```
//!
//! Default thread count (`available_parallelism()`) is intentional
//! per Step 9 readiness check spec — measure what a downstream
//! user would actually experience without tuning.

use std::time::Instant;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::mantle::{
    MantleConfig, COUPLING_DEFAULT, MF_DEFAULT, NUM_MODES_DEFAULT,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::amg::AmgConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

const NX: usize = 32;
const STEPS: usize = 100;
const RUNS: usize = 5;
const SEED: u64 = 42;

fn build_step3_config(linear_solver: LinearSolverConfig) -> BaselineConfig {
    let scales = Scales::default();
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: NX,
        grid_ny: NX,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_baseline_32sq_step3"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Disabled,
        boundary: BoundaryConfig::Disabled,
        boundary_layout_name: String::new(),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver,
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
    }
}

fn build_step6_config(linear_solver: LinearSolverConfig) -> BaselineConfig {
    let scales = Scales::default();
    let vcfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX,
        NX,
        &vcfg,
        SEED,
        rates,
        RecyclingConfig::default(),
    )
    .expect("step6 32² boundary build");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: NX,
        grid_ny: NX,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: STEPS,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_baseline_32sq_step6"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: 0.05,
            ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n8", SEED),
        slab_pull: SlabPullConfig::Disabled,
        mantle: MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver,
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
    }
}

fn build_step7_config(linear_solver: LinearSolverConfig) -> BaselineConfig {
    // Step 7 baseline shape minus slab-pull (Step 7 regression
    // convention exception — see Step 8 README note). Same geometry
    // as step6 but with a slab-pull-disabled marker; in practice
    // this is identical to step6 for the linear solver (slab pipe
    // does not run).
    let mut cfg = build_step6_config(linear_solver);
    cfg.boundary_layout_name = format!("voronoi_seed{}_n8_slab_off", SEED);
    cfg.output_dir = std::env::temp_dir().join("ymir_baseline_32sq_step7");
    cfg
}

fn build_step8_config(linear_solver: LinearSolverConfig) -> BaselineConfig {
    // Step 8 baseline accepted regime: mantle Enabled, slab held
    // Disabled per Step 8 co-calibration finding (slab + mantle
    // produced runaway G > 1; baseline ships slab-disabled
    // pending follow-up).
    let mut cfg = build_step6_config(linear_solver);
    cfg.mantle = MantleConfig::Enabled {
        mf: MF_DEFAULT,
        coupling: COUPLING_DEFAULT,
        num_modes: NUM_MODES_DEFAULT,
        seed: 7,
        evolution_rate: 0.0,
    };
    cfg.output_dir = std::env::temp_dir().join("ymir_baseline_32sq_step8");
    cfg
}

struct RunSummary {
    wc_mean: f64,
    wc_std: f64,
    cg_mean: f64,
    newton_mean: f64,
    fallback_rate: f64,
    convergence_pct: f64,
    peak_v: f64,
    mass_drift: f64,
    yf_max: f64,
    cg_capped_count: usize,
}

fn measure(label: &str, ls: LinearSolverConfig, mk: &dyn Fn(LinearSolverConfig) -> BaselineConfig) -> RunSummary {
    let mut wc = Vec::with_capacity(RUNS);
    let mut cg_mean = 0.0;
    let mut newton_mean = 0.0;
    let mut fallback_rate = 0.0;
    let mut converged_total = 0usize;
    let mut steps_total = 0usize;
    let mut peak_v = 0.0;
    let mut mass_drift = 0.0;
    let mut yf_max = 0.0;
    let mut cg_capped_count = 0usize;
    for run in 0..RUNS {
        let cfg = mk(ls.clone());
        let t0 = Instant::now();
        let r = run_baseline(&cfg);
        let dt = t0.elapsed().as_secs_f64();
        wc.push(dt);
        cg_mean = r.metrics.cg_iter_mean;
        let stats = r
            .metrics
            .extrapolation
            .as_ref()
            .expect("ExtrapolationStats present for steps > 0");
        newton_mean = stats.newton_outer_iters_mean();
        fallback_rate = stats.fallback_rate();
        let na = r
            .metrics
            .newton
            .as_ref()
            .expect("NewtonAggregate present");
        converged_total += na.converged;
        steps_total += na.converged + na.stalled + na.diverged + na.capped;
        peak_v = r.metrics.vmax_peak;
        mass_drift = r.metrics.mass_drift_relative;
        yf_max = r.metrics.yielding_cell_fraction.unwrap_or(0.0);
        cg_capped_count = r.metrics.cg_iter_max;
        eprintln!(
            "    [{}/{}] {} wc={:.2}s cg_mean={:.1} newton_mean={:.2} fallback={:.1}%",
            run + 1,
            RUNS,
            label,
            dt,
            cg_mean,
            newton_mean,
            fallback_rate * 100.0,
        );
    }
    let wc_mean = wc.iter().sum::<f64>() / wc.len() as f64;
    let wc_var = wc.iter().map(|x| (x - wc_mean).powi(2)).sum::<f64>() / wc.len() as f64;
    let wc_std = wc_var.sqrt();
    let convergence_pct = if steps_total > 0 {
        converged_total as f64 / steps_total as f64
    } else {
        0.0
    };
    RunSummary {
        wc_mean,
        wc_std,
        cg_mean,
        newton_mean,
        fallback_rate,
        convergence_pct,
        peak_v,
        mass_drift,
        yf_max,
        cg_capped_count,
    }
}

/// Focused re-run for step6 + step8 only — used when an OS-level
/// event (machine sleep) during the full sweep corrupted one
/// run's wallclock with a multi-hour outlier. Step3 and step7
/// data from the prior sweep are preserved when they are valid.
#[test]
#[ignore]
fn baseline_32sq_step6_step8_only() {
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    eprintln!("=== Step 8.5b 32² baseline RE-RUN (step6 + step8 only) ===");
    eprintln!("Hardware threads (available_parallelism): {}", threads);
    eprintln!("RAYON_NUM_THREADS env: {:?}", std::env::var("RAYON_NUM_THREADS").ok());
    eprintln!(
        "Configuration: NX={NX}, STEPS={STEPS}, RUNS={RUNS}, SEED={SEED}"
    );
    eprintln!();

    let cases: &[(&str, &dyn Fn(LinearSolverConfig) -> BaselineConfig)] = &[
        ("step6_voronoi", &build_step6_config),
        ("step8_mantle_on_slab_off", &build_step8_config),
    ];

    for (case, mk) in cases {
        eprintln!("--- {} ---", case);
        let jac = measure(
            &format!("{} Jacobi", case),
            LinearSolverConfig::JacobiCG,
            *mk,
        );
        let amg = measure(
            &format!("{} AMG", case),
            LinearSolverConfig::AmgCG(AmgConfig::default()),
            *mk,
        );
        eprintln!(
            "| {} | Jacobi | {:.2} ± {:.2} | {:.1} | {:.2} | {:.1} | {:.0} | {:.4e} | {:.4e} | {:.4e} | {} |",
            case,
            jac.wc_mean, jac.wc_std, jac.cg_mean, jac.newton_mean,
            jac.fallback_rate * 100.0, jac.convergence_pct * 100.0,
            jac.peak_v, jac.mass_drift, jac.yf_max, jac.cg_capped_count,
        );
        eprintln!(
            "| {} | AMG | {:.2} ± {:.2} | {:.1} | {:.2} | {:.1} | {:.0} | {:.4e} | {:.4e} | {:.4e} | {} |",
            case,
            amg.wc_mean, amg.wc_std, amg.cg_mean, amg.newton_mean,
            amg.fallback_rate * 100.0, amg.convergence_pct * 100.0,
            amg.peak_v, amg.mass_drift, amg.yf_max, amg.cg_capped_count,
        );
        let ratio = amg.wc_mean / jac.wc_mean.max(1e-12);
        eprintln!("    -> AMG/Jacobi wallclock ratio: {:.3}x", ratio);
        eprintln!();
    }
}

#[test]
#[ignore]
fn baseline_32sq_measurement() {
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    eprintln!("=== Step 8.5b 32² baseline measurement ===");
    eprintln!("Hardware threads (available_parallelism): {}", threads);
    eprintln!("RAYON_NUM_THREADS env: {:?}", std::env::var("RAYON_NUM_THREADS").ok());
    eprintln!(
        "Configuration: NX={NX}, STEPS={STEPS}, RUNS={RUNS}, SEED={SEED}"
    );
    eprintln!();

    let cases: &[(&str, &dyn Fn(LinearSolverConfig) -> BaselineConfig)] = &[
        ("step3_floor_yielding", &build_step3_config),
        ("step6_voronoi", &build_step6_config),
        ("step7_slab_off", &build_step7_config),
        ("step8_mantle_on_slab_off", &build_step8_config),
    ];

    eprintln!("| Case | Precond | Wallclock (s) | CG mean | Newton mean | Fallback % | Conv % | peak\\|v\\| | mass_drift | yf_max | CG max |");
    eprintln!("|---|---|---|---|---|---|---|---|---|---|---|");
    for (case, mk) in cases {
        eprintln!();
        eprintln!("--- {} ---", case);

        let jac = measure(
            &format!("{} Jacobi", case),
            LinearSolverConfig::JacobiCG,
            *mk,
        );
        let amg = measure(
            &format!("{} AMG", case),
            LinearSolverConfig::AmgCG(AmgConfig::default()),
            *mk,
        );

        eprintln!(
            "| {} | Jacobi | {:.2} ± {:.2} | {:.1} | {:.2} | {:.1} | {:.0} | {:.4e} | {:.4e} | {:.4e} | {} |",
            case,
            jac.wc_mean,
            jac.wc_std,
            jac.cg_mean,
            jac.newton_mean,
            jac.fallback_rate * 100.0,
            jac.convergence_pct * 100.0,
            jac.peak_v,
            jac.mass_drift,
            jac.yf_max,
            jac.cg_capped_count,
        );
        eprintln!(
            "| {} | AMG | {:.2} ± {:.2} | {:.1} | {:.2} | {:.1} | {:.0} | {:.4e} | {:.4e} | {:.4e} | {} |",
            case,
            amg.wc_mean,
            amg.wc_std,
            amg.cg_mean,
            amg.newton_mean,
            amg.fallback_rate * 100.0,
            amg.convergence_pct * 100.0,
            amg.peak_v,
            amg.mass_drift,
            amg.yf_max,
            amg.cg_capped_count,
        );
        let ratio = amg.wc_mean / jac.wc_mean.max(1e-12);
        eprintln!("    -> AMG/Jacobi wallclock ratio: {:.3}x", ratio);
    }
}
