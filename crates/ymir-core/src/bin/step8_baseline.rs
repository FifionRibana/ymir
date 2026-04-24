//! Step 8 baseline driver.
//!
//! Runs three scenarios and emits three markdown reports:
//!
//! 1. **Physics** (`step8_physics_report.md`) — **Step 6 setup**
//!    (GPE + yielding + basal drag + Voronoi + dynamic detection
//!    + Closed recycling) plus `MantleConfig::Enabled (Mf=1.0,
//!    coupling=1.0, num_modes=6, seed=42, evolution_rate=0)`.
//!    **Slab-pull deliberately held Disabled** — see the
//!    "Step 8 regression convention exception" note in
//!    `tectonics_v2/README.md` and the "Slab+Mantle interaction
//!    instability finding" section of the generated physics
//!    report. Yielding STRICT > 0 last-chance acceptance is met
//!    in this setup.
//! 2. **Regression** (`step8_regression_report.md`) — identical
//!    setup with `MantleConfig::Disabled`, which reduces to
//!    **Step 6 physics**. Ratios compared vs
//!    `step6_physics_report.md`. Scalar parity expected by
//!    construction (no mantle, no slab = Step 6 bit-identity).
//! 3. **Mf sweep** (`step8_mf_sweep_report.md`) —
//!    `Mf ∈ {0.3, 0.6, 1.0, 1.5, 2.0}` at 64², **single seed**
//!    across all points (the Fourier pattern is fixed, only the
//!    amplitude varies — contrast with Step 6's `num_plates`
//!    sweep where topology varied). Slab-pull held Disabled
//!    here too.

use std::path::PathBuf;
use std::process::ExitCode;

use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::diagnostics::comparison::parse_step_report;
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, ForceKind, NonlinearChoice, build_force, run_baseline,
};
use ymir_core::tectonics_v2::diagnostics::mms_bench;
use ymir_core::tectonics_v2::diagnostics::report::{
    ReportInputs, ReportKind, default_previous_report_for, write_markdown_report,
};
use ymir_core::tectonics_v2::mantle::{
    COUPLING_DEFAULT, MF_DEFAULT, MantleConfig, NUM_MODES_DEFAULT,
};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

struct Args {
    seed: u64,
    grids: Vec<(usize, usize)>,
    steps: usize,
    output_dir: PathBuf,
    preset: Preset,
    num_plates: usize,
    continental_ratio: f64,
    mf: f64,
    coupling: f64,
    mantle_seed: u64,
    mantle_num_modes: usize,
    mf_sweep: bool,
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
        mf: MF_DEFAULT,
        coupling: COUPLING_DEFAULT,
        mantle_seed: 42,
        mantle_num_modes: NUM_MODES_DEFAULT,
        mf_sweep: true,
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--seed" => {
                i += 1;
                a.seed = args[i].parse().map_err(|e| format!("bad --seed: {e}"))?;
            }
            "--grids" => {
                i += 1;
                let mut grids = Vec::new();
                for tok in args[i].split(',') {
                    let tok = tok.trim();
                    if tok.is_empty() {
                        continue;
                    }
                    if let Some((x, y)) = tok.split_once('x') {
                        let nx: usize = x.parse().map_err(|e| format!("{tok}: {e}"))?;
                        let ny: usize = y.parse().map_err(|e| format!("{tok}: {e}"))?;
                        grids.push((nx, ny));
                    } else {
                        let n: usize = tok.parse().map_err(|e| format!("{tok}: {e}"))?;
                        grids.push((n, n));
                    }
                }
                if !grids.is_empty() {
                    a.grids = grids;
                }
            }
            "--steps" => {
                i += 1;
                a.steps = args[i].parse().map_err(|e| format!("bad --steps: {e}"))?;
            }
            "--output-dir" | "--output" => {
                i += 1;
                a.output_dir = PathBuf::from(&args[i]);
            }
            "--preset" => {
                i += 1;
                a.preset = Preset::by_name(&args[i])?;
            }
            "--num-plates" => {
                i += 1;
                a.num_plates = args[i].parse().map_err(|e| format!("bad --num-plates: {e}"))?;
            }
            "--continental-ratio" => {
                i += 1;
                a.continental_ratio =
                    args[i].parse().map_err(|e| format!("bad --continental-ratio: {e}"))?;
            }
            "--mf" => {
                i += 1;
                a.mf = args[i].parse().map_err(|e| format!("bad --mf: {e}"))?;
            }
            "--coupling" => {
                i += 1;
                a.coupling = args[i].parse().map_err(|e| format!("bad --coupling: {e}"))?;
            }
            "--mantle-seed" => {
                i += 1;
                a.mantle_seed = args[i].parse().map_err(|e| format!("bad --mantle-seed: {e}"))?;
            }
            "--mantle-num-modes" => {
                i += 1;
                a.mantle_num_modes =
                    args[i].parse().map_err(|e| format!("bad --mantle-num-modes: {e}"))?;
            }
            "--no-sweep" => {
                a.mf_sweep = false;
            }
            "--help" | "-h" => {
                println!(
                    "Usage: step8_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset NAME] [--num-plates N] [--continental-ratio F] \
                     [--mf F] [--coupling F] [--mantle-seed N] [--mantle-num-modes N] \
                     [--no-sweep] [--output-dir PATH]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    Ok(a)
}

/// Build a Step 8-flavour BaselineConfig: Step 7 physics
/// (Voronoï + Closed recycling; slab-pull Disabled per Step 8 exception)
/// with the supplied `MantleConfig`. The mantle toggle is the
/// only axis of variation between the physics and regression
/// scenarios.
fn build_baseline_config(
    args: &Args,
    nx: usize,
    ny: usize,
    mantle: MantleConfig,
    output_subdir: &str,
) -> Result<BaselineConfig, String> {
    let scales = Scales::default();
    let vcfg =
        VoronoiConfig { num_plates: args.num_plates, continental_ratio: args.continental_ratio };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let recycling_config = RecyclingConfig::default();
    let boundary =
        BoundaryConfig::enabled_voronoi_closed(nx, ny, &vcfg, args.seed, rates, recycling_config)
            .map_err(|e| format!("recycling config invalid: {:?}", e))?;
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    // Slab-pull is held Disabled for the Step 8 baseline, both
    // physics and regression, pending the slab+mantle co-
    // calibration issue raised during Step 8 development. See
    // `docs/reports/step8_physics_report.md §Slab+Mantle
    // interaction instability finding` and the regression
    // convention exception note in `tectonics_v2/README.md`.
    let slab = SlabPullConfig::Disabled;
    Ok(BaselineConfig {
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
        output_dir: args.output_dir.join(output_subdir),
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
        slab_pull: slab,
        mantle,
        capture: None,
        linear_solver: Default::default(),
    })
}

fn run_one(
    label: &str,
    cfg: BaselineConfig,
) -> (
    ymir_core::tectonics_v2::diagnostics::metrics::Metrics,
    ymir_core::tectonics_v2::diagnostics::metrics::SolverConfigDump,
) {
    let nx = cfg.grid_nx;
    let ny = cfg.grid_ny;
    let steps = cfg.steps;
    println!("-- running {} {}×{} for {} steps --", label, nx, ny, steps);
    let r = run_baseline(&cfg);
    let m = &r.metrics;
    let yielding = m.newton.as_ref().and_then(|n| n.yielding_cell_fraction_max).unwrap_or(0.0);
    let peak_v_mantle = m
        .newton
        .as_ref()
        .and_then(|n| n.peak_v_solved_mantle_run)
        .unwrap_or(0.0);
    let eps_ratio =
        m.newton.as_ref().and_then(|n| n.epsilon_ii_max_to_floor_ratio).unwrap_or(0.0);
    let div_max = m.newton.as_ref().and_then(|n| n.div_v_mantle_max).unwrap_or(0.0);
    println!(
        "  wallclock: {:.3}s; CG/Newton mean: {:.1}/{}; mass_conservation_residual: {:.3e}; peak|v|: {:.3e} (peak|v_solved|: {:.3e}); yielding: {:.3e}; ε̇_II/floor: {:.2e}; div_v_mantle_max: {:.2e}",
        m.wallclock_total.as_secs_f64(),
        m.cg_iter_mean,
        m.cg_iter_max,
        m.newton.as_ref().and_then(|n| n.mass_conservation_residual).unwrap_or(0.0),
        m.vmax_peak,
        peak_v_mantle,
        yielding,
        eps_ratio,
        div_max,
    );
    (r.metrics, r.config_dump)
}

#[derive(Clone, Copy)]
struct MfSweepRow {
    mf: f64,
    peak_v: f64,
    yielding_cell_fraction_max: f64,
    epsilon_ii_max_to_floor_ratio: f64,
    newton_converged_rate: f64,
    cg_iter_mean: f64,
    wallclock_s: f64,
    mass_conservation_residual: f64,
}

fn run_mf_sweep(args: &Args) -> Vec<MfSweepRow> {
    let mf_values = [0.3_f64, 0.6, 1.0, 1.5, 2.0];
    let mut rows = Vec::new();
    for &mf in &mf_values {
        let mantle = MantleConfig::Enabled {
            mf,
            coupling: args.coupling,
            num_modes: args.mantle_num_modes,
            seed: args.mantle_seed,
            evolution_rate: 0.0,
        };
        let cfg = match build_baseline_config(args, 64, 64, mantle, "step8_mf_sweep_heightmaps") {
            Ok(c) => c,
            Err(e) => {
                eprintln!("mf={mf}: build failed: {e}");
                continue;
            }
        };
        println!("-- Mf sweep: Mf = {} --", mf);
        let r = run_baseline(&cfg);
        let m = &r.metrics;
        let (conv_pct, _, _, _) = m
            .newton
            .as_ref()
            .map(|n| n.outcome_percentages())
            .unwrap_or((0.0, 0.0, 0.0, 0.0));
        rows.push(MfSweepRow {
            mf,
            peak_v: m
                .newton
                .as_ref()
                .and_then(|n| n.peak_v_solved_mantle_run)
                .unwrap_or(m.vmax_peak),
            yielding_cell_fraction_max: m
                .newton
                .as_ref()
                .and_then(|n| n.yielding_cell_fraction_max)
                .unwrap_or(0.0),
            epsilon_ii_max_to_floor_ratio: m
                .newton
                .as_ref()
                .and_then(|n| n.epsilon_ii_max_to_floor_ratio)
                .unwrap_or(0.0),
            newton_converged_rate: conv_pct,
            cg_iter_mean: m.cg_iter_mean,
            wallclock_s: m.wallclock_total.as_secs_f64(),
            mass_conservation_residual: m
                .newton
                .as_ref()
                .and_then(|n| n.mass_conservation_residual)
                .unwrap_or(0.0),
        });
    }
    rows
}

fn render_mf_sweep_markdown(rows: &[MfSweepRow], args: &Args) -> String {
    let mut s = String::new();
    s.push_str("# Step 8 — Mf sweep (peak|v_solved| scaling, yielding activation threshold)\n\n");
    s.push_str(&format!(
        "> Fixed: coupling = `{:.3}`, num_modes = `{}`, mantle_seed = `{}`, world seed = `{}`, \
         num_plates = `{}`, grid 64², steps = `{}`. **Single seed across all points** — the \
         Fourier pattern is fixed; only the amplitude `Mf` varies.\n\n",
        args.coupling,
        args.mantle_num_modes,
        args.mantle_seed,
        args.seed,
        args.num_plates,
        args.steps,
    ));
    s.push_str("| Mf | peak\\|v_solved\\| | yielding_cell_fraction_max | ε̇_II / ε̇_min | Newton converged (%) | CG iters (mean) | wallclock (s) | mass_conservation_residual |\n");
    s.push_str("|---|---|---|---|---|---|---|---|\n");
    for r in rows {
        s.push_str(&format!(
            "| `{:.2}` | `{:.3e}` | `{:.3e}` | `{:.3e}` | `{:.1}` | `{:.1}` | `{:.2}` | `{:.3e}` |\n",
            r.mf,
            r.peak_v,
            r.yielding_cell_fraction_max,
            r.epsilon_ii_max_to_floor_ratio,
            r.newton_converged_rate,
            r.cg_iter_mean,
            r.wallclock_s,
            r.mass_conservation_residual,
        ));
    }
    s.push('\n');
    // Monotonicity diagnostic.
    let mut monotonic = true;
    for w in rows.windows(2) {
        if w[1].peak_v + 1e-12 < w[0].peak_v {
            monotonic = false;
            break;
        }
    }
    if monotonic {
        s.push_str(
            "**Monotonicity: ✅ `peak|v_solved|` monotonically non-decreasing with `Mf`.** Expected from the linear coupling of mantle amplitude to forcing. Non-linear saturation (sub-linear growth) is acceptable at the top end — the full-field response includes viscous dissipation and (through Newton) the power-law rheology.\n\n"
        );
    } else {
        s.push_str("**Flag — non-monotonic `peak|v_solved|` with `Mf`.** The amplitude scaling is linear by construction (fixed pattern times `Mf`); a non-monotonic response signals either a solver convergence issue at a specific point, continuation-ramp interaction, or numerical noise dominating. Investigate the specific row before accepting the sweep.\n\n");
    }
    // Yielding activation threshold.
    let first_active = rows
        .iter()
        .position(|r| r.yielding_cell_fraction_max > 0.0)
        .map(|idx| rows[idx].mf);
    if let Some(mf_crit) = first_active {
        s.push_str(&format!(
            "**Yielding activation threshold (observed).** Yielding first fires at `Mf ≥ {:.2}` in this sweep. The critical `Mf` is a physical property measured, not prescribed; at smaller `Mf` the mantle bootstrap does not push `ε̇_II` above the regularisation floor.\n\n",
            mf_crit,
        ));
    } else {
        s.push_str("**Yielding activation threshold (observed).** Yielding did not fire at any `Mf ∈ {0.3, 0.6, 1.0, 1.5, 2.0}`. Combined with a similar result at Step 8 baseline (Mf=1.0), this would be a pathological signal — the mechanism hierarchy failed to initiate. Cross-reference the Step 8 physics report's yielding checkpoint section.\n\n");
    }
    s
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("argument error: {e}");
            return ExitCode::from(2);
        }
    };
    let scales = Scales::default();
    println!("{}", scales.report());
    println!(
        "preset: {} | seed: {} | num_plates: {} | grids: {:?} | steps: {} | Mf: {:.3} | coupling: {:.3} | num_modes: {} | mantle_seed: {}",
        args.preset.name,
        args.seed,
        args.num_plates,
        args.grids,
        args.steps,
        args.mf,
        args.coupling,
        args.mantle_num_modes,
        args.mantle_seed,
    );
    std::fs::create_dir_all(&args.output_dir).ok();

    println!("-- running MMS bench --");
    let mms = mms_bench::run_all();

    // -------- Physics --------
    println!("\n=== Step 8 physics baseline (mantle Enabled, slab-pull Disabled per Step 8 exception) ===");
    let mut phys_metrics = Vec::new();
    let mut phys_configs = Vec::new();
    for (nx, ny) in &args.grids {
        let mantle = MantleConfig::Enabled {
            mf: args.mf,
            coupling: args.coupling,
            num_modes: args.mantle_num_modes,
            seed: args.mantle_seed,
            evolution_rate: 0.0,
        };
        let cfg =
            match build_baseline_config(&args, *nx, *ny, mantle, "step8_physics_heightmaps") {
                Ok(c) => c,
                Err(e) => {
                    eprintln!("physics build failed: {e}");
                    return ExitCode::from(1);
                }
            };
        let (m, c) = run_one("Step 8 physics", cfg);
        phys_metrics.push(m);
        phys_configs.push(c);
    }

    // -------- Regression --------
    println!("\n=== Step 8 regression (mantle Disabled, mirror of Step 7 physics) ===");
    let mut reg_metrics = Vec::new();
    let mut reg_configs = Vec::new();
    for (nx, ny) in &args.grids {
        let cfg = match build_baseline_config(
            &args,
            *nx,
            *ny,
            MantleConfig::Disabled,
            "step8_regression_heightmaps",
        ) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("regression build failed: {e}");
                return ExitCode::from(1);
            }
        };
        let (m, c) = run_one("Step 8 regression", cfg);
        reg_metrics.push(m);
        reg_configs.push(c);
    }

    // -------- Mf sweep (optional) --------
    let mf_rows = if args.mf_sweep {
        println!("\n=== Step 8 Mf sweep (64² only, single mantle seed) ===");
        let rows = run_mf_sweep(&args);
        for r in &rows {
            println!(
                "  Mf={:.2}  peak|v_solved|={:.3e}  yielding={:.3e}  ε̇_II/floor={:.2e}  newton_conv={:.1}%  cg_mean={:.1}  wallclock={:.2}s",
                r.mf,
                r.peak_v,
                r.yielding_cell_fraction_max,
                r.epsilon_ii_max_to_floor_ratio,
                r.newton_converged_rate,
                r.cg_iter_mean,
                r.wallclock_s,
            );
        }
        Some(rows)
    } else {
        None
    };

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
    let phys_previous = resolve_previous(ReportKind::Step8Physics);
    let phys_output = args.output_dir.join("step8_physics_report.md");
    let phys_inputs = ReportInputs {
        kind: ReportKind::Step8Physics,
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
    let reg_previous = resolve_previous(ReportKind::Step8Regression);
    let reg_output = args.output_dir.join("step8_regression_report.md");
    let reg_inputs = ReportInputs {
        kind: ReportKind::Step8Regression,
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

    // Mf sweep report (minimal markdown — rendered directly,
    // not through the ReportInputs pipeline).
    if let Some(rows) = mf_rows {
        let sweep_md = render_mf_sweep_markdown(&rows, &args);
        let sweep_output = args.output_dir.join("step8_mf_sweep_report.md");
        if let Err(e) = std::fs::write(&sweep_output, sweep_md) {
            eprintln!("failed to write Mf sweep report: {}", e);
            return ExitCode::from(1);
        }
        println!("Mf sweep report → {}", sweep_output.display());
    }
    ExitCode::SUCCESS
}
