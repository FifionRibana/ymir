//! Step 7 baseline driver.
//!
//! Runs three scenarios and emits three markdown reports:
//!
//! 1. **Physics** (`step7_physics_report.md`) — Step 6 setup
//!    (Voronoi num_plates=8, Closed-mode recycling) with
//!    `SlabPullConfig::Enabled(Sp=1.5, τ_slab=0.5, k_slab_accum=1.0,
//!    ε=1e-6)`. `yielding_cell_fraction_max > 0` is a strict
//!    acceptance.
//! 2. **Regression** (`step7_regression_report.md`) — identical
//!    Step 6 physics setup with `SlabPullConfig::Disabled`.
//!    Zero-cost-when-disabled invariant; ratios vs Step 6 physics
//!    in `[0.95, 1.05]`.
//! 3. **Sp sweep** (`step7_sp_sweep_report.md`) —
//!    `Sp ∈ {0.5, 1.0, 1.5, 2.0, 3.0}` at 64² × `steps` with all
//!    other parameters fixed at baseline.
//!
//! ```bash
//! cargo run --release --bin step7_baseline -- \
//!     --grids 64,128 --steps 300 --output-dir docs/reports/
//! ```

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
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::{
    EPSILON_DEFAULT, K_SLAB_ACCUM_DEFAULT, SP_DEFAULT, SlabPullConfig, TAU_SLAB_DEFAULT,
};
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

struct Args {
    seed: u64,
    grids: Vec<(usize, usize)>,
    steps: usize,
    output_dir: PathBuf,
    preset: Preset,
    num_plates: usize,
    continental_ratio: f64,
    sp: f64,
    tau_slab: f64,
    k_slab_accum: f64,
    epsilon: f64,
    sp_sweep: bool,
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
        sp: SP_DEFAULT,
        tau_slab: TAU_SLAB_DEFAULT,
        k_slab_accum: K_SLAB_ACCUM_DEFAULT,
        epsilon: EPSILON_DEFAULT,
        sp_sweep: true,
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
            "--sp" => {
                i += 1;
                a.sp = args[i].parse().map_err(|e| format!("bad --sp: {e}"))?;
            }
            "--tau-slab" => {
                i += 1;
                a.tau_slab = args[i].parse().map_err(|e| format!("bad --tau-slab: {e}"))?;
            }
            "--k-slab-accum" => {
                i += 1;
                a.k_slab_accum = args[i].parse().map_err(|e| format!("bad --k-slab-accum: {e}"))?;
            }
            "--epsilon" => {
                i += 1;
                a.epsilon = args[i].parse().map_err(|e| format!("bad --epsilon: {e}"))?;
            }
            "--no-sweep" => {
                a.sp_sweep = false;
            }
            "--help" | "-h" => {
                println!(
                    "Usage: step7_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset NAME] [--num-plates N] [--continental-ratio F] \
                     [--sp F] [--tau-slab F] [--k-slab-accum F] [--epsilon F] \
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

/// Build a common Step 7-flavour BaselineConfig: Step 6 physics
/// (Voronoi + Closed recycling) with the supplied `SlabPullConfig`.
/// The slab toggle is the only axis of variation between the
/// physics and regression scenarios.
fn build_baseline_config(
    args: &Args,
    nx: usize,
    ny: usize,
    slab: SlabPullConfig,
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
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw { br: 0.05, ..BasalDragLaw::default() }),
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", args.seed, args.num_plates),
        slab_pull: slab,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: Default::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
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
    let peak_f_slab = m.newton.as_ref().and_then(|n| n.peak_f_slab_run).unwrap_or(0.0);
    println!(
        "  wallclock: {:.3}s; CG/Newton mean: {:.1}/{}; mass_conservation_residual: {:.3e}; vmax: {:.3e}; yielding: {:.3}; peak|f_slab|: {:.3e}",
        m.wallclock_total.as_secs_f64(),
        m.cg_iter_mean,
        m.cg_iter_max,
        m.newton.as_ref().and_then(|n| n.mass_conservation_residual).unwrap_or(0.0),
        m.vmax_peak,
        yielding,
        peak_f_slab,
    );
    (r.metrics, r.config_dump)
}

#[derive(Clone, Copy)]
struct SpSweepRow {
    sp: f64,
    peak_v: f64,
    m_subducted_max: f64,
    yielding_cell_fraction_max: f64,
    newton_converged_rate: f64,
    cg_iter_mean: f64,
    wallclock_s: f64,
}

fn run_sp_sweep(args: &Args) -> Vec<SpSweepRow> {
    let sp_values = [0.5_f64, 1.0, 1.5, 2.0, 3.0];
    let mut rows = Vec::new();
    for &sp in &sp_values {
        let slab = SlabPullConfig::Enabled {
            sp,
            tau_slab: args.tau_slab,
            k_slab_accum: args.k_slab_accum,
            epsilon: args.epsilon,
        };
        let cfg = match build_baseline_config(args, 64, 64, slab, "step7_sp_sweep_heightmaps") {
            Ok(c) => c,
            Err(e) => {
                eprintln!("sp={sp}: build failed: {e}");
                continue;
            }
        };
        println!("-- Sp sweep: Sp = {} --", sp);
        let r = run_baseline(&cfg);
        let m = &r.metrics;
        let (conv_pct, _, _, _) =
            m.newton.as_ref().map(|n| n.outcome_percentages()).unwrap_or((0.0, 0.0, 0.0, 0.0));
        rows.push(SpSweepRow {
            sp,
            peak_v: m.vmax_peak,
            m_subducted_max: m
                .newton
                .as_ref()
                .map(|n| n.slab_m_max_series.iter().cloned().fold(0.0_f64, f64::max))
                .unwrap_or(0.0),
            yielding_cell_fraction_max: m
                .newton
                .as_ref()
                .and_then(|n| n.yielding_cell_fraction_max)
                .unwrap_or(0.0),
            newton_converged_rate: conv_pct,
            cg_iter_mean: m.cg_iter_mean,
            wallclock_s: m.wallclock_total.as_secs_f64(),
        });
    }
    rows
}

fn render_sp_sweep_markdown(rows: &[SpSweepRow], args: &Args) -> String {
    let mut s = String::new();
    s.push_str("# Step 7 — Sp sweep (peak|v| monotonicity)\n\n");
    s.push_str(&format!(
        "> Fixed: τ_slab = `{:.3}`, k_slab_accum = `{:.3}`, ε = `{:.1e}`, seed = `{}`, num_plates = `{}`, grid 64², steps = `{}`.\n\n",
        args.tau_slab, args.k_slab_accum, args.epsilon, args.seed, args.num_plates, args.steps,
    ));
    s.push_str("| Sp | peak\\|v\\| | m_subducted_max | yielding_cell_fraction_max | Newton converged (%) | CG iters (mean) | wallclock (s) |\n");
    s.push_str("|---|---|---|---|---|---|---|\n");
    for r in rows {
        s.push_str(&format!(
            "| `{:.2}` | `{:.3e}` | `{:.3e}` | `{:.3e}` | `{:.1}` | `{:.1}` | `{:.2}` |\n",
            r.sp,
            r.peak_v,
            r.m_subducted_max,
            r.yielding_cell_fraction_max,
            r.newton_converged_rate,
            r.cg_iter_mean,
            r.wallclock_s,
        ));
    }
    s.push('\n');
    // Monotonicity + variation diagnostic.
    let mut monotonic = true;
    for w in rows.windows(2) {
        if w[1].peak_v + 1e-12 < w[0].peak_v {
            monotonic = false;
            break;
        }
    }
    let peak_v_min = rows.iter().map(|r| r.peak_v).fold(f64::INFINITY, f64::min);
    let peak_v_max = rows.iter().map(|r| r.peak_v).fold(0.0_f64, f64::max);
    let flat = peak_v_max - peak_v_min < 1.0e-10 * peak_v_max.max(1.0e-30);
    if !monotonic {
        s.push_str("**Flag — non-monotonic `peak|v|`.** The expected linear coupling `f_slab ∝ Sp · m̃` should yield monotonically non-decreasing `peak|v|` with Sp. A reversal suggests either (a) Sp-dependent convergence failures masking the signal, (b) interaction with the Newton continuation ramp, or (c) numerical noise dominating at low Sp. Investigate before submit.\n\n");
    } else if flat {
        s.push_str(
            "**Interpretation — flat across the Sp band (bootstrap failure regime).**\n\n\
             `peak|v|` is identical across `Sp ∈ [0.5, 3.0]` to f64 precision. This is the signature of the closed-loop gain `G = Sp · k_slab_accum · τ_slab / (η · L)` sitting `≪ 1` everywhere in the §4.8 target band, with the floor-dominated `η_newton ≈ 100`. The quiescent fixed point is linearly stable; the system remains at the Step 6 baseline regardless of `Sp`. `peak|f_slab|` does scale linearly with `Sp` (visible in the physics report), but the Stokes inversion `v = f · L²/η` damps it by `1/η ≈ 0.01`, so no measurable `peak|v|` response. Monotonicity is trivially satisfied (zero difference).\n\n\
             This is consistent with the amplifier-vs-initiator revision documented in `step7_physics_report.md §Yielding checkpoint`. A non-flat sweep is expected once Step 8 (mantle forcing) imposes `v_mantle` externally and breaks the floor-dominated regime; slab-pull will then amplify visibly.\n\n",
        );
    } else {
        s.push_str("**Interpretation.** `peak|v|` is monotonically non-decreasing with `Sp`, confirming the linear coupling of slab-pull strength (§4.8). A plateau at the top end would be the signature of τ_slab-limited saturation (m̃_steady ≈ Q·τ independent of Sp after some point).\n\n");
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
        "preset: {} | seed: {} | num_plates: {} | grids: {:?} | steps: {} | Sp: {:.3} | τ_slab: {:.3} | k_slab_accum: {:.3} | ε: {:.1e}",
        args.preset.name,
        args.seed,
        args.num_plates,
        args.grids,
        args.steps,
        args.sp,
        args.tau_slab,
        args.k_slab_accum,
        args.epsilon,
    );
    std::fs::create_dir_all(&args.output_dir).ok();

    println!("-- running MMS bench --");
    let mms = mms_bench::run_all();

    // -------- Physics --------
    println!("\n=== Step 7 physics baseline (slab-pull Enabled) ===");
    let mut phys_metrics = Vec::new();
    let mut phys_configs = Vec::new();
    for (nx, ny) in &args.grids {
        let slab = SlabPullConfig::Enabled {
            sp: args.sp,
            tau_slab: args.tau_slab,
            k_slab_accum: args.k_slab_accum,
            epsilon: args.epsilon,
        };
        let cfg = match build_baseline_config(&args, *nx, *ny, slab, "step7_physics_heightmaps") {
            Ok(c) => c,
            Err(e) => {
                eprintln!("physics build failed: {e}");
                return ExitCode::from(1);
            }
        };
        let (m, c) = run_one("Step 7 physics", cfg);
        phys_metrics.push(m);
        phys_configs.push(c);
    }

    // -------- Regression --------
    println!("\n=== Step 7 regression (slab-pull Disabled, mirror of Step 6 physics) ===");
    let mut reg_metrics = Vec::new();
    let mut reg_configs = Vec::new();
    for (nx, ny) in &args.grids {
        let cfg = match build_baseline_config(
            &args,
            *nx,
            *ny,
            SlabPullConfig::Disabled,
            "step7_regression_heightmaps",
        ) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("regression build failed: {e}");
                return ExitCode::from(1);
            }
        };
        let (m, c) = run_one("Step 7 regression", cfg);
        reg_metrics.push(m);
        reg_configs.push(c);
    }

    // -------- Sp sweep (optional) --------
    let sp_rows = if args.sp_sweep {
        println!("\n=== Step 7 Sp sweep (64² only) ===");
        let rows = run_sp_sweep(&args);
        for r in &rows {
            println!(
                "  Sp={:.2}  peak|v|={:.3e}  m_max={:.3e}  yielding={:.3}  newton_conv={:.1}%  cg_mean={:.1}  wallclock={:.2}s",
                r.sp,
                r.peak_v,
                r.m_subducted_max,
                r.yielding_cell_fraction_max,
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
    let phys_previous = resolve_previous(ReportKind::Step7Physics);
    let phys_output = args.output_dir.join("step7_physics_report.md");
    let phys_inputs = ReportInputs {
        kind: ReportKind::Step7Physics,
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
    let reg_previous = resolve_previous(ReportKind::Step7Regression);
    let reg_output = args.output_dir.join("step7_regression_report.md");
    let reg_inputs = ReportInputs {
        kind: ReportKind::Step7Regression,
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

    // Sp sweep report (minimal markdown — rendered directly, not
    // through the ReportInputs pipeline, since the sweep data
    // doesn't fit the per-grid-config structure).
    if let Some(rows) = sp_rows {
        let sweep_md = render_sp_sweep_markdown(&rows, &args);
        let sweep_output = args.output_dir.join("step7_sp_sweep_report.md");
        if let Err(e) = std::fs::write(&sweep_output, sweep_md) {
            eprintln!("failed to write Sp sweep report: {}", e);
            return ExitCode::from(1);
        }
        println!("Sp sweep report → {}", sweep_output.display());
    }
    ExitCode::SUCCESS
}
