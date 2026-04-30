//! Binary: run the Step 5 scenarios (physics baseline, reference
//! variant, regression) + the k_spread calibration + the k_sub sweep
//! and emit three markdown reports.
//!
//! The reference variant is "Step 4 physics configuration with
//! yielding Enabled". It is produced on this branch because the
//! merged Step 4 physics ran with yielding `Disabled` for Br
//! isolation, which does not match the Step 5+ regression convention
//! ("activate all mechanisms through N-1"). The reference variant
//! serves as the comparison target for the Step 5 regression's
//! zero-cost-when-disabled invariant.
//!
//! ```bash
//! cargo run --release --bin step5_baseline -- \
//!     --seed 42 --grids 64,128 --steps 300 \
//!     --layout horizontal_oceanic_strip \
//!     --output-dir docs/reports/
//! ```

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{
    boundary_type_diversity, build_layout, calibrate_k_spread, BoundaryConfig, BoundaryLayout,
    BoundaryMechanismActive, BoundaryRates, KSpreadCalibration,
};
use ymir_core::tectonics_v2::diagnostics::comparison::{parse_step_report, StepReference};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::diagnostics::k_sub_sweep::{run_k_sub_sweep, KSubSweepResults};
use ymir_core::tectonics_v2::diagnostics::mms_bench;
use ymir_core::tectonics_v2::diagnostics::report::{
    default_previous_report_for, write_markdown_report, ReportInputs, ReportKind,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;

struct Args {
    seed: u64,
    grids: Vec<(usize, usize)>,
    steps: usize,
    output_dir: PathBuf,
    preset: Preset,
    layout_name: String,
    k_sub: f64,
    bi: f64,
    br: f64,
    /// Optional pre-computed `k_spread`. When supplied, the
    /// calibration step is skipped and the report's calibration
    /// section says so.
    k_spread_override: Option<f64>,
    /// Budget for the calibration simulation (in 64² steps per probe).
    /// Separate from `steps` to keep calibration cheap relative to
    /// the full baseline.
    calibration_steps: usize,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        seed: 42,
        grids: vec![(64, 64), (128, 128)],
        steps: 300,
        output_dir: PathBuf::from("docs/reports/"),
        preset: Preset::by_name("dynamic-accidented")?,
        layout_name: "horizontal_oceanic_strip".into(),
        k_sub: 0.5,
        bi: YieldingLaw::default().bi,
        br: BasalDragLaw::default().br,
        k_spread_override: None,
        // Match physics step count by default so the calibrated
        // k_spread lands the 300-step baseline in
        // `s_oceanic_mean ∈ [0.18, 0.22]` rather than the shorter-run
        // calibration target. Overridable via `--calibration-steps`.
        calibration_steps: 0, // 0 → "inherit from --steps", resolved after parse
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--seed" => { i += 1; a.seed = args[i].parse().map_err(|e| format!("bad --seed: {e}"))?; }
            "--grids" => {
                i += 1;
                let mut grids: Vec<(usize, usize)> = Vec::new();
                for tok in args[i].split(',') {
                    let tok = tok.trim();
                    if tok.is_empty() { continue; }
                    if let Some((x, y)) = tok.split_once('x') {
                        let nx: usize = x.parse().map_err(|e| format!("bad grid {tok}: {e}"))?;
                        let ny: usize = y.parse().map_err(|e| format!("bad grid {tok}: {e}"))?;
                        grids.push((nx, ny));
                    } else {
                        let n: usize = tok.parse().map_err(|e| format!("bad grid {tok}: {e}"))?;
                        grids.push((n, n));
                    }
                }
                if grids.is_empty() { return Err("--grids produced no grid sizes".into()); }
                a.grids = grids;
            }
            "--steps" => { i += 1; a.steps = args[i].parse().map_err(|e| format!("bad --steps: {e}"))?; }
            "--output-dir" | "--output" => { i += 1; a.output_dir = PathBuf::from(&args[i]); }
            "--preset" => { i += 1; a.preset = Preset::by_name(&args[i])?; }
            "--layout" => { i += 1; a.layout_name = args[i].clone(); }
            "--k-sub" => { i += 1; a.k_sub = args[i].parse().map_err(|e| format!("bad --k-sub: {e}"))?; }
            "--bi" => { i += 1; a.bi = args[i].parse().map_err(|e| format!("bad --bi: {e}"))?; }
            "--br" => { i += 1; a.br = args[i].parse().map_err(|e| format!("bad --br: {e}"))?; }
            "--k-spread" => { i += 1; a.k_spread_override = Some(args[i].parse().map_err(|e| format!("bad --k-spread: {e}"))?); }
            "--calibration-steps" => { i += 1; a.calibration_steps = args[i].parse().map_err(|e| format!("bad --calibration-steps: {e}"))?; }
            "--help" | "-h" => {
                println!(
                    "Usage: step5_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset NAME] [--layout NAME] [--k-sub F] [--bi F] [--br F] \
                     [--k-spread F] [--calibration-steps N] [--output-dir PATH]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    if a.calibration_steps == 0 {
        a.calibration_steps = a.steps;
    }
    Ok(a)
}

fn resolve_previous(
    kind: ReportKind,
    output_dir: &std::path::Path,
) -> Option<StepReference> {
    let path = default_previous_report_for(kind, output_dir);
    if !path.exists() {
        eprintln!("note: previous report not found at {:?} — comparison block omitted", path);
        return None;
    }
    match parse_step_report(&path) {
        Ok(r) => {
            println!("{:?} → parsed previous report {:?} ({} grids)", kind, path, r.grids.len());
            Some(r)
        }
        Err(e) => {
            eprintln!("warning: could not parse {:?}: {}", path, e);
            None
        }
    }
}

/// Ascii heatmap of a plate-type + flag layout. Used by the Step 5
/// physics report's "Layout visualization" section. Renders the
/// domain at sub-sampled resolution so the 64² grid fits on one
/// screen width.
fn layout_ascii(layout: &BoundaryLayout) -> String {
    let nx = layout.nx();
    let ny = layout.ny();
    // Keep every row so single-row features (rift, subduction) are
    // always visible. Horizontally sub-sample if the grid is wide.
    let stride_x = (nx / 64).max(1);
    let stride_y: usize = 1;
    let mut out = String::new();
    out.push_str("plate_types (.=Oceanic, #=Continental)       boundary_flags (.=None, r=Rift, s=Subd, S=OcSubd, C=ContColl)\n");
    for j in (0..ny).step_by(stride_y).rev() {
        for i in (0..nx).step_by(stride_x) {
            out.push(match layout.plate_types.get(i, j) {
                ymir_core::tectonics_v2::boundaries::PlateType::Oceanic => '.',
                ymir_core::tectonics_v2::boundaries::PlateType::Continental => '#',
            });
        }
        out.push_str("     ");
        for i in (0..nx).step_by(stride_x) {
            out.push(match layout.flags.get(i, j) {
                ymir_core::tectonics_v2::boundaries::BoundaryFlag::None => '.',
                ymir_core::tectonics_v2::boundaries::BoundaryFlag::Rift => 'r',
                ymir_core::tectonics_v2::boundaries::BoundaryFlag::Subduction => 's',
                ymir_core::tectonics_v2::boundaries::BoundaryFlag::OceanicSubduction => 'S',
                ymir_core::tectonics_v2::boundaries::BoundaryFlag::ContinentalCollision => 'C',
            });
        }
        out.push('\n');
    }
    out
}

/// Run one boundary-enabled 64²·calibration_steps simulation at the
/// supplied `k_spread` and return the final `s_oceanic_mean`. Used
/// as the calibration closure.
fn calibration_probe(
    args: &Args,
    k_spread: f64,
) -> f64 {
    let nx = 64;
    let ny = 64;
    let scales = Scales::default();
    let layout = build_layout(&args.layout_name, nx, ny).expect("valid layout");
    let layout_name = layout.name.to_string();
    let rates = BoundaryRates::baseline_uncalibrated()
        .with_k_spread(k_spread)
        .with_k_sub(args.k_sub);
    let boundary = layout.into_config(rates);
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    let cfg = BaselineConfig {
        seed: args.seed,
        grid_nx: nx,
        grid_ny: ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: args.calibration_steps,
        cfl_factor: 0.3,
        // Match the physics baseline's simulated time (6·τ*) so the
        // calibrated `k_spread` reproduces the physics run's
        // `s_oceanic_mean`. The spec-suggested `3·τ*` shortened
        // calibration but created a gap between the calibrator's
        // steady state and the 6·τ* physics run — the 64² baseline
        // saw `s_oceanic_mean = 0.2301` at 6·τ* while the 3·τ* probe
        // reported 0.2151 (in band). Match times to close the gap.
        total_time_nondim: 6.0,
        preset: args.preset.clone(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: Default::default(),
        picard_cfg: Default::default(),
        heightmap_fractions: Vec::new(),
        output_dir: PathBuf::from("target/step5_calibration_scratch"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding: YieldingConfig::Enabled(YieldingLaw { bi: args.bi, ..Default::default() }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: args.br,
            ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: layout_name,
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
    };
    let r = run_baseline(&cfg);
    r.metrics
        .newton
        .as_ref()
        .and_then(|n| n.s_oceanic_mean)
        .unwrap_or(f64::NAN)
}

struct ScenarioResult {
    configs: Vec<ymir_core::tectonics_v2::diagnostics::metrics::SolverConfigDump>,
    metrics: Vec<ymir_core::tectonics_v2::diagnostics::metrics::Metrics>,
    max_vmax: f64,
}

fn run_scenario_multi_grid(
    args: &Args,
    boundary_builder: impl Fn(usize, usize) -> (BoundaryConfig, String),
    yielding: YieldingConfig,
    basal_drag: BasalDragConfig,
    heightmap_subdir_name: &str,
    scales: &Scales,
) -> ScenarioResult {
    let mut configs = Vec::new();
    let mut metrics = Vec::new();
    let heightmap_subdir = args.output_dir.join(heightmap_subdir_name);
    for (nx, ny) in &args.grids {
        let (boundary, layout_name) = boundary_builder(*nx, *ny);
        let force = build_force(ForceKind::Gpe, scales, 10.0, 1.0);
        let cfg = BaselineConfig {
            seed: args.seed,
            grid_nx: *nx,
            grid_ny: *ny,
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
            output_dir: heightmap_subdir.clone(),
            force,
            force_kind: ForceKind::Gpe,
            sinusoidal_amplitude: 0.0,
            s_perturbation_amplitude: 0.2,
            yielding,
            basal_drag,
            boundary,
            boundary_layout_name: layout_name,
            slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
            mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
            cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
            age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
            capture: None,
            linear_solver: Default::default(),
            init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        };
        println!("-- running {}×{} for {} steps --", nx, ny, args.steps);
        let r = run_baseline(&cfg);
        println!(
            "  wallclock: {:.3}s; CG/Newton: {:.1}/{}; mass drift: {:.3e}; vmax: {:.3e}",
            r.metrics.wallclock_total.as_secs_f64(),
            r.metrics.cg_iter_mean,
            r.metrics.cg_iter_max,
            r.metrics.mass_drift_relative,
            r.metrics.vmax_peak,
        );
        if let Some(na) = &r.metrics.newton {
            if let (Some(so), Some(mbr)) = (na.s_oceanic_mean, na.mass_balance_residual) {
                println!(
                    "  s_oceanic_mean={:.4}, mass_balance_residual={:.3e}, clamp_frac_mean={:.3e}",
                    so,
                    mbr,
                    na.clamp_activation_fraction_mean.unwrap_or(0.0),
                );
            }
        }
        configs.push(r.config_dump);
        metrics.push(r.metrics);
    }
    let max_vmax = metrics.iter().map(|m| m.vmax_peak).fold(0.0_f64, f64::max);
    ScenarioResult { configs, metrics, max_vmax }
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => { eprintln!("argument error: {e}"); return ExitCode::from(2); }
    };
    let scales = Scales::default();
    println!("{}", scales.report());
    println!(
        "preset: {} | layout: {} | seed: {} | grids: {:?} | steps: {}",
        args.preset.name, args.layout_name, args.seed, args.grids, args.steps,
    );
    std::fs::create_dir_all(&args.output_dir).ok();

    // -------- MMS --------
    println!("-- running MMS bench --");
    let mms = mms_bench::run_all();

    // -------- k_spread calibration --------
    let calibration = match args.k_spread_override {
        Some(k) => {
            println!("--k-spread={} supplied → skipping calibration", k);
            ymir_core::tectonics_v2::boundaries::CalibrationResult {
                k_spread: k,
                iterations: vec![ymir_core::tectonics_v2::boundaries::calibration::CalibrationIter {
                    k_spread: k,
                    s_oceanic_mean: f64::NAN,
                }],
                final_s_oceanic_mean: f64::NAN,
            }
        }
        None => {
            println!("-- running k_spread calibration --");
            let cal_cfg = KSpreadCalibration::step5_default();
            let t0 = Instant::now();
            let calibration = match calibrate_k_spread(&cal_cfg, |k| {
                let s = calibration_probe(&args, k);
                println!("  probe k_spread={:.4} → s_oceanic_mean={:.4}", k, s);
                s
            }) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("calibration failed: {:?}", e);
                    return ExitCode::from(1);
                }
            };
            println!(
                "  calibration done in {:.1}s: k_spread={:.4}, s_oceanic_mean={:.4}",
                t0.elapsed().as_secs_f64(),
                calibration.k_spread,
                calibration.final_s_oceanic_mean,
            );
            calibration
        }
    };

    let k_spread_calibrated = calibration.k_spread;

    // -------- Physics baseline (Step 5) --------
    println!("\n=== Step 5 physics baseline ===");
    let physics_builder = |nx: usize, ny: usize| {
        let layout = build_layout(&args.layout_name, nx, ny).expect("layout");
        let name = layout.name.to_string();
        let rates = BoundaryRates::baseline_uncalibrated()
            .with_k_spread(k_spread_calibrated)
            .with_k_sub(args.k_sub);
        (layout.into_config(rates), name)
    };
    let physics = run_scenario_multi_grid(
        &args,
        physics_builder,
        YieldingConfig::Enabled(YieldingLaw { bi: args.bi, ..Default::default() }),
        BasalDragConfig::Enabled(BasalDragLaw { br: args.br, ..BasalDragLaw::default() }),
        "step5_physics_heightmaps",
        &scales,
    );

    // Annotate the physics metrics with `k_spread_calibrated` and
    // layout name so the report picks them up.
    let mut physics_metrics = physics.metrics;
    for m in physics_metrics.iter_mut() {
        if let Some(na) = m.newton.as_mut() {
            na.k_spread_calibrated = Some(k_spread_calibrated);
        }
    }

    // -------- k_sub sweep --------
    println!("\n=== k_sub sweep ===");
    let k_sub_values = [0.3_f64, 0.5, 0.7, 1.0];
    let k_sub_sweep: KSubSweepResults = run_k_sub_sweep(
        args.seed,
        args.steps,
        &args.preset,
        0.2,
        &k_sub_values,
        k_spread_calibrated,
        args.bi,
        args.br,
    );
    for p in &k_sub_sweep.points {
        println!(
            "  k_sub={:.2}: s_oceanic={:.8} peak|v|={:.3e} CG={:.1} wallclock={:.2}s",
            p.k_sub,
            p.s_oceanic_mean.unwrap_or(f64::NAN),
            p.peak_v,
            p.cg_iter_mean,
            p.wallclock_s,
        );
    }
    println!(
        "  monotonicity: s_oceanic strictly ↓ with k_sub? {}",
        k_sub_sweep.s_oceanic_mono_ok,
    );

    // -------- Reference variant (Step 4 physics-yielding-Enabled) --------
    println!("\n=== Step 5 reference variant (Step 4 physics with yielding Enabled) ===");
    let reference_builder = |_nx: usize, _ny: usize| {
        (BoundaryConfig::Disabled, String::new())
    };
    let reference_variant = run_scenario_multi_grid(
        &args,
        reference_builder,
        YieldingConfig::Enabled(YieldingLaw { bi: args.bi, ..Default::default() }),
        BasalDragConfig::Enabled(BasalDragLaw { br: args.br, ..BasalDragLaw::default() }),
        "step5_reference_variant_heightmaps",
        &scales,
    );

    // Emit the reference-variant report.
    let ref_output = args.output_dir.join("step5_reference_variant_report.md");
    let ref_justifications: Vec<String> = vec![String::new(); reference_variant.metrics.len()];
    let ref_previous = resolve_previous(ReportKind::Step5ReferenceVariant, &args.output_dir);
    let ref_inputs = ReportInputs {
        kind: ReportKind::Step5ReferenceVariant,
        seed: args.seed,
        scales: &scales,
        configs: &reference_variant.configs,
        metrics: &reference_variant.metrics,
        previous: ref_previous.as_ref(),
        suspect_justifications: &ref_justifications,
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
    if let Err(e) = write_markdown_report(&ref_output, &ref_inputs) {
        eprintln!("failed to write reference variant report: {}", e);
        return ExitCode::from(1);
    }
    println!("reference variant report → {}", ref_output.display());

    // -------- Regression (same as reference variant) --------
    // The regression run IS the reference variant; the variant's
    // own report is flavoured "physics" (for the header). The
    // regression report below is flavoured as a regression mirror,
    // compared against the variant. They consume the same metrics.
    let reg_output = args.output_dir.join("step5_regression_report.md");
    let reg_justifications: Vec<String> = vec![String::new(); reference_variant.metrics.len()];
    let reg_previous = resolve_previous(ReportKind::Step5Regression, &args.output_dir);
    let reg_inputs = ReportInputs {
        kind: ReportKind::Step5Regression,
        seed: args.seed,
        scales: &scales,
        configs: &reference_variant.configs,
        metrics: &reference_variant.metrics,
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

    // -------- Physics report (rendered last so it can reference sweeps) --------
    // Layout ascii for the physics report.
    let layout_viz = build_layout(&args.layout_name, 64, 64).ok().map(|l| layout_ascii(&l));
    let phys_output = args.output_dir.join("step5_physics_report.md");
    let phys_justifications: Vec<String> = vec![String::new(); physics_metrics.len()];
    let phys_previous = resolve_previous(ReportKind::Step5Physics, &args.output_dir);
    let phys_inputs = ReportInputs {
        kind: ReportKind::Step5Physics,
        seed: args.seed,
        scales: &scales,
        configs: &physics.configs,
        metrics: &physics_metrics,
        previous: phys_previous.as_ref(),
        suspect_justifications: &phys_justifications,
        mms: Some(&mms),
        ar_sweep: None,
        bi_sweep: None,
        br_sweep: None,
        regression_vmax_peak: Some(reference_variant.max_vmax),
        k_sub_sweep: Some(&k_sub_sweep),
        k_spread_calibration: Some(&calibration),
        boundary_layout_ascii: layout_viz,
        num_plates_sweep: None,
    };
    if let Err(e) = write_markdown_report(&phys_output, &phys_inputs) {
        eprintln!("failed to write physics report: {}", e);
        return ExitCode::from(1);
    }
    println!("physics report → {}", phys_output.display());

    // Minor diagnostic: echo the boundary diversity observed at 64²
    // physics so the operator sees the layout is effectively active.
    if let Some(na) = physics_metrics.first().and_then(|m| m.newton.as_ref()) {
        println!(
            "diversity={} (expected nonzero when the layout has ≥1 boundary-type cell)",
            na.boundary_type_diversity.unwrap_or(0),
        );
    }
    let _ = physics.max_vmax;
    let _ = boundary_type_diversity;
    let _ = BoundaryMechanismActive {
        sub: true, spread: true, coll_v: true, rift_v: true,
    };

    ExitCode::SUCCESS
}
