//! Step 6 baseline driver.
//!
//! Runs three scenarios and emits three markdown reports:
//!
//! 1. **Physics** (`step6_physics_report.md`) — Voronoi tessellation
//!    (num_plates=8, continental_ratio=0.3, seed=42), dynamic
//!    boundary detection, Closed-mode recycling (default fractions).
//! 2. **Regression** (`step6_regression_report.md`) — Step 5 physics
//!    configuration mirrored: static `horizontal_oceanic_strip`,
//!    `RecyclingMode::Open` with Step 5 rates (k_spread=0.05).
//!    Compared to Step 5 physics at `[0.95, 1.05]`.
//! 3. **Voronoi sweep** (`step6_voronoi_sweep_report.md`) —
//!    `num_plates ∈ {4, 8, 12, 16}` with distinct seeds per point
//!    `{42, 43, 44, 45}` at 64² × 300 steps.
//!
//! ```bash
//! cargo run --release --bin step6_baseline -- \
//!     --grids 64,128 --steps 300 --output-dir docs/reports/
//! ```

use std::path::PathBuf;
use std::process::ExitCode;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{
    horizontal_oceanic_strip, BoundaryConfig, BoundaryRates,
};
use ymir_core::tectonics_v2::diagnostics::comparison::parse_step_report;
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::diagnostics::mms_bench;
use ymir_core::tectonics_v2::diagnostics::num_plates_sweep::run_num_plates_sweep;
use ymir_core::tectonics_v2::diagnostics::report::{
    default_previous_report_for, write_markdown_report, ReportInputs, ReportKind,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

struct Args {
    seed: u64,
    grids: Vec<(usize, usize)>,
    steps: usize,
    output_dir: PathBuf,
    preset: Preset,
    num_plates: usize,
    continental_ratio: f64,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        seed: 42,
        grids: vec![(64, 64), (128, 128)],
        steps: 300,
        output_dir: PathBuf::from("docs/reports/"),
        preset: Preset::by_name("dynamic-accidented")?,
        num_plates: 8,
        continental_ratio: 0.3,
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--seed" => { i += 1; a.seed = args[i].parse().map_err(|e| format!("bad --seed: {e}"))?; }
            "--grids" => {
                i += 1;
                let mut grids = Vec::new();
                for tok in args[i].split(',') {
                    let tok = tok.trim();
                    if tok.is_empty() { continue; }
                    if let Some((x, y)) = tok.split_once('x') {
                        let nx: usize = x.parse().map_err(|e| format!("{tok}: {e}"))?;
                        let ny: usize = y.parse().map_err(|e| format!("{tok}: {e}"))?;
                        grids.push((nx, ny));
                    } else {
                        let n: usize = tok.parse().map_err(|e| format!("{tok}: {e}"))?;
                        grids.push((n, n));
                    }
                }
                if !grids.is_empty() { a.grids = grids; }
            }
            "--steps" => { i += 1; a.steps = args[i].parse().map_err(|e| format!("bad --steps: {e}"))?; }
            "--output-dir" | "--output" => { i += 1; a.output_dir = PathBuf::from(&args[i]); }
            "--preset" => { i += 1; a.preset = Preset::by_name(&args[i])?; }
            "--num-plates" => { i += 1; a.num_plates = args[i].parse().map_err(|e| format!("bad --num-plates: {e}"))?; }
            "--continental-ratio" => { i += 1; a.continental_ratio = args[i].parse().map_err(|e| format!("bad --continental-ratio: {e}"))?; }
            "--help" | "-h" => {
                println!(
                    "Usage: step6_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset NAME] [--num-plates N] [--continental-ratio F] \
                     [--output-dir PATH]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    Ok(a)
}

fn run_voronoi_config(
    args: &Args,
    nx: usize,
    ny: usize,
) -> (
    ymir_core::tectonics_v2::diagnostics::metrics::Metrics,
    ymir_core::tectonics_v2::diagnostics::metrics::SolverConfigDump,
) {
    let scales = Scales::default();
    let vcfg = VoronoiConfig {
        num_plates: args.num_plates,
        continental_ratio: args.continental_ratio,
    };
    // k_sub retained as the only rate; k_spread/k_arc/k_coll_v/k_rift_v
    // are ignored in Closed mode, so we zero them to make the config
    // dump unambiguous.
    let rates = BoundaryRates {
        k_sub: 0.5,
        k_arc: 0.0,
        k_spread: 0.0,
        k_coll_v: 0.0,
        k_rift_v: 0.0,
    };
    let recycling_config = RecyclingConfig::default();
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        nx, ny, &vcfg, args.seed, rates, recycling_config,
    )
    .expect("recycling config valid");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: args.seed,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: args.steps,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: args.preset.clone(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: vec![0.0, 0.5, 1.0],
        output_dir: args.output_dir.join("step6_physics_heightmaps"),
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
        boundary_layout_name: format!("voronoi_seed{}_n{}", args.seed, args.num_plates),
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
    };
    println!("-- running Voronoi physics {}×{} for {} steps --", nx, ny, args.steps);
    let r = run_baseline(&cfg);
    println!(
        "  wallclock: {:.3}s; CG/Newton mean: {:.1}/{}; mass_conservation_residual: {:.3e}; vmax: {:.3e}",
        r.metrics.wallclock_total.as_secs_f64(),
        r.metrics.cg_iter_mean,
        r.metrics.cg_iter_max,
        r.metrics.newton.as_ref().and_then(|n| n.mass_conservation_residual).unwrap_or(0.0),
        r.metrics.vmax_peak,
    );
    (r.metrics, r.config_dump)
}

fn run_regression_config(
    args: &Args,
    nx: usize,
    ny: usize,
) -> (
    ymir_core::tectonics_v2::diagnostics::metrics::Metrics,
    ymir_core::tectonics_v2::diagnostics::metrics::SolverConfigDump,
) {
    let scales = Scales::default();
    let layout = horizontal_oceanic_strip(nx, ny);
    let rates = BoundaryRates::baseline_uncalibrated().with_k_spread(0.05);
    let boundary = layout.into_config(rates);
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: args.seed,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: args.steps,
        cfl_factor: 0.3,
        total_time_nondim: 6.0,
        preset: args.preset.clone(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: vec![0.0, 0.5, 1.0],
        output_dir: args.output_dir.join("step6_regression_heightmaps"),
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
        boundary_layout_name: "horizontal_oceanic_strip".into(),
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
    };
    println!("-- running Step 5-shape regression {}×{} for {} steps --", nx, ny, args.steps);
    let r = run_baseline(&cfg);
    println!(
        "  wallclock: {:.3}s; CG/Newton mean: {:.1}/{}; mass_balance_residual: {:.3e}",
        r.metrics.wallclock_total.as_secs_f64(),
        r.metrics.cg_iter_mean,
        r.metrics.cg_iter_max,
        r.metrics.newton.as_ref().and_then(|n| n.mass_balance_residual).unwrap_or(0.0),
    );
    (r.metrics, r.config_dump)
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => { eprintln!("argument error: {e}"); return ExitCode::from(2); }
    };
    let scales = Scales::default();
    println!("{}", scales.report());
    println!(
        "preset: {} | seed: {} | num_plates: {} | continental_ratio: {:.2} | grids: {:?} | steps: {}",
        args.preset.name, args.seed, args.num_plates, args.continental_ratio, args.grids, args.steps,
    );
    std::fs::create_dir_all(&args.output_dir).ok();

    println!("-- running MMS bench --");
    let mms = mms_bench::run_all();

    // -------- Physics --------
    println!("\n=== Step 6 physics baseline ===");
    let mut phys_metrics = Vec::new();
    let mut phys_configs = Vec::new();
    for (nx, ny) in &args.grids {
        let (m, c) = run_voronoi_config(&args, *nx, *ny);
        phys_metrics.push(m);
        phys_configs.push(c);
    }

    // -------- Regression --------
    println!("\n=== Step 6 regression (Step 5-shape + Open mode) ===");
    let mut reg_metrics = Vec::new();
    let mut reg_configs = Vec::new();
    for (nx, ny) in &args.grids {
        let (m, c) = run_regression_config(&args, *nx, *ny);
        reg_metrics.push(m);
        reg_configs.push(c);
    }

    // -------- Voronoi sweep --------
    println!("\n=== Voronoi num_plates × seed sweep ===");
    let sweep = run_num_plates_sweep(
        &[4, 8, 12, 16],
        &[42, 43, 44, 45],
        args.steps,
        &args.preset,
        0.2,
    );
    for p in &sweep.points {
        println!(
            "  num_plates={:>2} seed={}: plate_count={:?} cont_frac={:.3} s_oceanic={:.4} s_cont={:.4} residual={:.3e} wallclock={:.2}s",
            p.num_plates, p.seed,
            p.plate_count, p.continental_fraction.unwrap_or(f64::NAN),
            p.s_oceanic_mean.unwrap_or(f64::NAN),
            p.s_continental_interior_mean.unwrap_or(f64::NAN),
            p.mass_conservation_residual.unwrap_or(f64::NAN),
            p.wallclock_s,
        );
    }

    // -------- Reports --------
    let resolve_previous = |kind: ReportKind| {
        let path = default_previous_report_for(kind, &args.output_dir);
        if path.exists() {
            parse_step_report(&path).ok()
        } else {
            eprintln!("note: previous report not found at {:?} — comparison block omitted", path);
            None
        }
    };

    // Physics.
    let phys_justifications: Vec<String> = vec![String::new(); phys_metrics.len()];
    let phys_previous = resolve_previous(ReportKind::Step6Physics);
    let phys_output = args.output_dir.join("step6_physics_report.md");
    let phys_inputs = ReportInputs {
        kind: ReportKind::Step6Physics,
        seed: args.seed,
        scales: &scales,
        configs: &phys_configs,
        metrics: &phys_metrics,
        previous: phys_previous.as_ref(),
        suspect_justifications: &phys_justifications,
        mms: Some(&mms),
        ar_sweep: None,
        bi_sweep: None,
        br_sweep: None,
        regression_vmax_peak: None,
        k_sub_sweep: None,
        k_spread_calibration: None,
        boundary_layout_ascii: None,
        num_plates_sweep: None,
    };
    if let Err(e) = write_markdown_report(&phys_output, &phys_inputs) {
        eprintln!("failed to write physics report: {}", e);
        return ExitCode::from(1);
    }
    println!("physics report → {}", phys_output.display());

    // Regression.
    let reg_justifications: Vec<String> = vec![String::new(); reg_metrics.len()];
    let reg_previous = resolve_previous(ReportKind::Step6Regression);
    let reg_output = args.output_dir.join("step6_regression_report.md");
    let reg_inputs = ReportInputs {
        kind: ReportKind::Step6Regression,
        seed: args.seed,
        scales: &scales,
        configs: &reg_configs,
        metrics: &reg_metrics,
        previous: reg_previous.as_ref(),
        suspect_justifications: &reg_justifications,
        mms: Some(&mms),
        ar_sweep: None,
        bi_sweep: None,
        br_sweep: None,
        regression_vmax_peak: None,
        k_sub_sweep: None,
        k_spread_calibration: None,
        boundary_layout_ascii: None,
        num_plates_sweep: None,
    };
    if let Err(e) = write_markdown_report(&reg_output, &reg_inputs) {
        eprintln!("failed to write regression report: {}", e);
        return ExitCode::from(1);
    }
    println!("regression report → {}", reg_output.display());

    // Voronoi sweep report.
    let sweep_output = args.output_dir.join("step6_voronoi_sweep_report.md");
    let sweep_inputs = ReportInputs {
        kind: ReportKind::Step6VoronoiSweep,
        seed: args.seed,
        scales: &scales,
        configs: &[],
        metrics: &[],
        previous: None,
        suspect_justifications: &[],
        mms: None,
        ar_sweep: None,
        bi_sweep: None,
        br_sweep: None,
        regression_vmax_peak: None,
        k_sub_sweep: None,
        k_spread_calibration: None,
        boundary_layout_ascii: None,
        num_plates_sweep: Some(&sweep),
    };
    if let Err(e) = write_markdown_report(&sweep_output, &sweep_inputs) {
        eprintln!("failed to write sweep report: {}", e);
        return ExitCode::from(1);
    }
    println!("sweep report → {}", sweep_output.display());
    // TODO Step 6 follow-up: a dedicated sweep renderer would be cleaner.
    ExitCode::SUCCESS
}
