//! Binary: run the Step 2 baseline scenarios and emit markdown
//! reports.
//!
//! ```bash
//! cargo run --release --bin step_baseline -- \
//!     --seed 42 --grids 64,128 --steps 300 \
//!     --preset dynamic-accidented --nonlinear-solver newton \
//!     --forcing both \
//!     --compare-to docs/reports/step1_report.md \
//!     --output-dir docs/reports/
//! ```
//!
//! With `--forcing both` (default), the binary emits
//! `step2_physics_report.md` (GpeForce) and
//! `step2_regression_report.md` (SinusoidalForce ε=10, mirror of
//! Step 1) side by side. Individual scenarios can be selected with
//! `--forcing gpe` or `--forcing sinusoidal`.

use std::path::PathBuf;
use std::process::ExitCode;

use ymir_core::tectonics_v2::diagnostics::ar_sweep::{self, ArSweepResults};
use ymir_core::tectonics_v2::diagnostics::bi_sweep::{self, BiSweepResults};
use ymir_core::tectonics_v2::diagnostics::br_sweep::{self, BrSweepResults};
use ymir_core::tectonics_v2::diagnostics::comparison::{parse_step_report, StepReference};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::diagnostics::mms_bench;
use ymir_core::tectonics_v2::diagnostics::report::{
    write_markdown_report, ReportInputs, ReportKind,
};
use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ForcingSelection {
    Gpe,
    Sinusoidal,
    Both,
}

impl ForcingSelection {
    fn parse(s: &str) -> Result<Self, String> {
        match s {
            "gpe" => Ok(Self::Gpe),
            "sinusoidal" => Ok(Self::Sinusoidal),
            "both" => Ok(Self::Both),
            other => Err(format!(
                "unknown --forcing value '{}'; expected gpe|sinusoidal|both",
                other,
            )),
        }
    }
}

struct Args {
    seed: u64,
    grids: Vec<(usize, usize)>,
    steps: usize,
    output_dir: PathBuf,
    preset: Preset,
    nonlinear: NonlinearChoice,
    /// Legacy shared-override for both scenarios. Still honoured for
    /// back-compat, but per-scenario flags below take precedence.
    compare_to: Option<PathBuf>,
    /// Physics-scenario comparison target. When `None`, the binary
    /// auto-detects via [`report::default_previous_report_for`], e.g.
    /// Step-4 physics → `step3_physics_report.md`. This is what
    /// prevents Step N from inheriting a leftover regression-mirror
    /// comparison that points at the wrong scenario's numbers.
    compare_to_physics: Option<PathBuf>,
    /// Regression-scenario comparison target. Same auto-detect
    /// fallback as `compare_to_physics`.
    compare_to_regression: Option<PathBuf>,
    forcing: ForcingSelection,
    sinusoidal_amplitude: f64,
    /// `None` → default per scenario (physics Enabled, regression
    /// Disabled). `Some(cfg)` → explicit override from CLI. The
    /// `--bi` flag is baked into the Enabled variant here so the
    /// scenario layer doesn't need to re-plumb it.
    yielding_override: Option<YieldingConfig>,
    /// `None` → default per scenario (Step-4 physics Enabled with
    /// `Br`, regression Disabled). `Some(cfg)` → explicit override.
    basal_drag_override: Option<BasalDragConfig>,
    /// Br value used when basal drag is Enabled (default 0.05). Read
    /// by `run_scenario` to inject into the per-scenario default.
    br: f64,
}

fn parse_args() -> Result<Args, String> {
    let mut seed: u64 = 42;
    let mut grids_str: String = "64,128".into();
    let mut steps: usize = 300;
    let mut output_dir: PathBuf = PathBuf::from("docs/reports/");
    let mut preset_name: String = "dynamic-accidented".into();
    let mut nonlinear_str: String = "newton".into();
    let mut compare_to: Option<PathBuf> = None;
    let mut compare_to_physics: Option<PathBuf> = None;
    let mut compare_to_regression: Option<PathBuf> = None;
    let mut forcing_str: String = "both".into();
    let mut sin_amp: f64 = 10.0;
    let mut yielding_str: Option<String> = None;
    let mut bi: f64 = YieldingLaw::default().bi;
    let mut basal_drag_str: Option<String> = None;
    let mut br: f64 = BasalDragLaw::default().br;

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--seed" => { i += 1; seed = args[i].parse().map_err(|e| format!("bad --seed: {e}"))?; }
            "--grids" => { i += 1; grids_str = args[i].clone(); }
            "--steps" => { i += 1; steps = args[i].parse().map_err(|e| format!("bad --steps: {e}"))?; }
            "--output-dir" | "--output" => { i += 1; output_dir = PathBuf::from(&args[i]); }
            "--preset" => { i += 1; preset_name = args[i].clone(); }
            "--nonlinear-solver" => { i += 1; nonlinear_str = args[i].clone(); }
            "--compare-to" => { i += 1; compare_to = Some(PathBuf::from(&args[i])); }
            "--compare-to-physics" => { i += 1; compare_to_physics = Some(PathBuf::from(&args[i])); }
            "--compare-to-regression" => { i += 1; compare_to_regression = Some(PathBuf::from(&args[i])); }
            "--forcing" => { i += 1; forcing_str = args[i].clone(); }
            "--sinusoidal-amplitude" => { i += 1; sin_amp = args[i].parse().map_err(|e| format!("bad --sinusoidal-amplitude: {e}"))?; }
            "--yielding-config" => { i += 1; yielding_str = Some(args[i].clone()); }
            "--bi" => { i += 1; bi = args[i].parse().map_err(|e| format!("bad --bi: {e}"))?; }
            "--basal-drag-config" => { i += 1; basal_drag_str = Some(args[i].clone()); }
            "--br" => { i += 1; br = args[i].parse().map_err(|e| format!("bad --br: {e}"))?; }
            "--help" | "-h" => {
                println!(
                    "Usage: step_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset NAME] [--nonlinear-solver newton|picard] \
                     [--forcing gpe|sinusoidal|both] [--sinusoidal-amplitude F] \
                     [--yielding-config enabled|disabled] [--bi F] \
                     [--basal-drag-config enabled|disabled] [--br F] \
                     [--compare-to PATH] \
                     [--compare-to-physics PATH] [--compare-to-regression PATH] \
                     [--output-dir PATH]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }

    let mut grids: Vec<(usize, usize)> = Vec::new();
    for tok in grids_str.split(',') {
        let tok = tok.trim();
        if tok.is_empty() { continue; }
        if let Some((a, b)) = tok.split_once('x') {
            let nx: usize = a.parse().map_err(|e| format!("bad grid {tok}: {e}"))?;
            let ny: usize = b.parse().map_err(|e| format!("bad grid {tok}: {e}"))?;
            grids.push((nx, ny));
        } else {
            let n: usize = tok.parse().map_err(|e| format!("bad grid {tok}: {e}"))?;
            grids.push((n, n));
        }
    }
    if grids.is_empty() {
        return Err("--grids produced no grid sizes".into());
    }

    let yielding_override = match yielding_str {
        Some(s) => {
            // Parse to a concrete variant, then re-inject the user's
            // `--bi` into the Enabled arm so the CLI's two flags
            // compose (bi only matters when Enabled).
            Some(match YieldingConfig::parse(&s)? {
                YieldingConfig::Disabled => YieldingConfig::Disabled,
                YieldingConfig::Enabled(_) => {
                    YieldingConfig::Enabled(YieldingLaw { bi, ..Default::default() })
                }
            })
        }
        None => None,
    };
    let basal_drag_override = match basal_drag_str {
        Some(s) => {
            // Mirror of yielding: parse the variant, then re-inject
            // --br into the Enabled arm.
            Some(match BasalDragConfig::parse(&s)? {
                BasalDragConfig::Disabled => BasalDragConfig::Disabled,
                BasalDragConfig::Enabled(_) => {
                    BasalDragConfig::Enabled(BasalDragLaw { br, ..BasalDragLaw::default() })
                }
            })
        }
        None => None,
    };
    Ok(Args {
        seed, grids, steps, output_dir,
        preset: Preset::by_name(&preset_name)?,
        nonlinear: NonlinearChoice::parse(&nonlinear_str)?,
        compare_to,
        compare_to_physics,
        compare_to_regression,
        forcing: ForcingSelection::parse(&forcing_str)?,
        sinusoidal_amplitude: sin_amp,
        yielding_override,
        basal_drag_override,
        br,
    })
}

fn run_scenario(
    args: &Args,
    kind: ForceKind,
    report_kind: ReportKind,
    report_file: &str,
    scales: &Scales,
    previous: Option<&StepReference>,
    mms: &mms_bench::MmsResults,
    ar_sweep_results: Option<&ArSweepResults>,
    bi_sweep_results: Option<&BiSweepResults>,
    br_sweep_results: Option<&BrSweepResults>,
    regression_vmax_peak: Option<f64>,
) -> Result<f64, String> {
    let heightmap_subdir = args.output_dir.join(match kind {
        ForceKind::Gpe => "step4_physics_heightmaps",
        ForceKind::Sinusoidal => "step4_regression_heightmaps",
    });

    let mut configs = Vec::new();
    let mut metrics = Vec::new();

    println!(
        "\n=== scenario: {} (writing {}) ===",
        kind.label(),
        report_file,
    );

    // Physics: amplified perturbation so GPE response is visible at
    // Ar = 0.1. Regression: 0.02 to preserve the Step-1 mirror.
    let s_amp = match kind {
        ForceKind::Gpe => 0.2,
        ForceKind::Sinusoidal => 0.02,
    };

    // Step 4: physics isolates the Br effect with yielding Disabled
    // (see issue #87 / prompt §Physique attendue). Regression stays
    // Disabled, as in Steps 2/3, to keep the mirror-of-previous-step
    // contract.
    //
    // The CLI's `--yielding-config` still overrides if provided —
    // it lets us re-run the Step-3 physics scenario for regression
    // comparison if needed.
    let yielding = args.yielding_override.unwrap_or(YieldingConfig::Disabled);

    // Basal drag default per-scenario: physics Enabled(Br),
    // regression Disabled. `--basal-drag-config` overrides.
    let basal_drag = args.basal_drag_override.unwrap_or_else(|| match kind {
        ForceKind::Gpe => BasalDragConfig::Enabled(BasalDragLaw {
            br: args.br,
            ..BasalDragLaw::default()
        }),
        ForceKind::Sinusoidal => BasalDragConfig::Disabled,
    });

    for (nx, ny) in &args.grids {
        let domain_lx = 1.0;
        let force = build_force(kind, scales, args.sinusoidal_amplitude, domain_lx);
        let base = BaselineConfig {
            seed: args.seed,
            grid_nx: *nx,
            grid_ny: *ny,
            domain_lx,
            domain_ly: 1.0,
            steps: args.steps,
            cfl_factor: 0.3,
            total_time_nondim: 6.0,
            preset: args.preset.clone(),
            nonlinear: args.nonlinear,
            newton_cfg: Default::default(),
            picard_cfg: Default::default(),
            heightmap_fractions: vec![0.0, 0.5, 1.0],
            output_dir: heightmap_subdir.clone(),
            force,
            force_kind: kind,
            sinusoidal_amplitude: args.sinusoidal_amplitude,
            s_perturbation_amplitude: s_amp,
            yielding,
            basal_drag,
            boundary: ymir_core::tectonics_v2::boundaries::BoundaryConfig::Disabled,
            boundary_layout_name: String::new(),
            slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
            mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
            cratonic: ymir_core::tectonics_v2::cratonic::CratonicConfig::Disabled,
            age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
            capture: None,
            linear_solver: Default::default(),
        };
        println!("-- running {}×{} for {} steps --", nx, ny, args.steps);
        let result = run_baseline(&base);
        println!(
            "  wallclock: {:.3}s; CG/Newton: {:.1}/{}; mass drift: {:.3e}; vmax: {:.3e}",
            result.metrics.wallclock_total.as_secs_f64(),
            result.metrics.cg_iter_mean,
            result.metrics.cg_iter_max,
            result.metrics.mass_drift_relative,
            result.metrics.vmax_peak,
        );
        if let Some(na) = &result.metrics.newton {
            let (c, s, d, cap) = na.outcome_percentages();
            println!(
                "  Newton outcomes — conv {:.1}%, stall {:.1}%, div {:.1}%, cap {:.1}%; η_max/η_min mean {:.2}, max {:.2}",
                c, s, d, cap,
                na.eta_contrast_mean(),
                na.eta_contrast_max(),
            );
        }
        configs.push(result.config_dump);
        metrics.push(result.metrics);
    }

    let output = args.output_dir.join(report_file);
    let justifications: Vec<String> = vec![String::new(); metrics.len()];
    let inputs = ReportInputs {
        kind: report_kind,
        seed: args.seed,
        scales,
        configs: &configs,
        metrics: &metrics,
        previous,
        suspect_justifications: &justifications,
        mms: Some(mms),
        ar_sweep: ar_sweep_results,
        bi_sweep: bi_sweep_results,
        br_sweep: br_sweep_results,
        regression_vmax_peak,
        k_sub_sweep: None,
        k_spread_calibration: None,
        boundary_layout_ascii: None,
        num_plates_sweep: None,
    };
    write_markdown_report(&output, &inputs)
        .map_err(|e| format!("failed to write report {:?}: {}", output, e))?;
    println!("report written to {}", output.display());
    // Max peak|v| across the scenario's grids — piped back into the
    // physics run as the regression reference for
    // `peak_v_damping_ratio`.
    let max_vmax = metrics.iter().map(|m| m.vmax_peak).fold(0.0_f64, f64::max);
    Ok(max_vmax)
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => { eprintln!("argument error: {e}"); return ExitCode::from(2); }
    };

    let scales = Scales::default();
    println!("{}", scales.report());
    println!(
        "preset: {} | nonlinear solver: {} | Ar (derived): {:.3} | grids: {:?} | steps: {}",
        args.preset.name,
        args.nonlinear.label(),
        scales.argand_number(),
        args.grids,
        args.steps,
    );

    // Resolve the "previous step" report to compare against, per
    // scenario. Priority:
    //   1. Explicit CLI override (`--compare-to-physics` /
    //      `--compare-to-regression`).
    //   2. Legacy shared `--compare-to` (applied to both scenarios).
    //   3. Auto-detect via `report::default_previous_report_for`,
    //      which for Step N {Physics,Regression} returns
    //      `{output}/step{N-1}_{physics,regression}_report.md`.
    //
    // The auto-detect step is what prevents future Step 5 from
    // inheriting a leftover regression-mirror comparison that
    // points at the wrong scenario (the bug Step 4 caught where
    // physics was being compared against Step-3 regression).
    let resolve_previous = |kind: ReportKind, override_path: Option<&PathBuf>| -> Option<StepReference> {
        let path = override_path
            .cloned()
            .or_else(|| args.compare_to.clone())
            .unwrap_or_else(|| ymir_core::tectonics_v2::diagnostics::report::default_previous_report_for(kind, &args.output_dir));
        if !path.exists() {
            eprintln!("note: previous-step report not found at {:?} — comparison block will be omitted", path);
            return None;
        }
        match parse_step_report(&path) {
            Ok(r) => {
                println!(
                    "{:?} → parsed previous report {:?} ({} grids)",
                    kind, path, r.grids.len(),
                );
                Some(r)
            }
            Err(e) => {
                eprintln!("warning: could not parse {:?}: {}", path, e);
                None
            }
        }
    };

    println!("-- running MMS bench for the reports --");
    let mms = mms_bench::run_all();
    println!(
        "  const-η slope: {:.3}; variable-η slope: {:.3}; Newton tail iters: {}",
        mms.const_eta.final_slope().unwrap_or(f64::NAN),
        mms.variable_eta.final_slope().unwrap_or(f64::NAN),
        mms.newton_tail.outer_iters,
    );

    let run_gpe = matches!(args.forcing, ForcingSelection::Gpe | ForcingSelection::Both);
    let run_sin = matches!(args.forcing, ForcingSelection::Sinusoidal | ForcingSelection::Both);

    // Step 4: GPE+basal-drag physics (yielding disabled to isolate
    // the Br effect), SinusoidalForce mirror for regression
    // (basal drag + yielding both disabled). Ar and Bi sweeps are
    // retained but not run by default — Step 4 runs the Br sweep.
    let _ = ar_sweep::run_ar_sweep; // keep symbol reachable without running the Step 2 sweep.
    let _ = bi_sweep::run_bi_sweep; // keep symbol reachable without running the Step 3 sweep.
    let br_sweep_res: Option<BrSweepResults> = if run_gpe {
        println!("-- running Br sweep (64²·{} steps × 5 points) --", args.steps);
        let br_values = [0.01_f64, 0.05, 0.10, 0.20, 0.30];
        let res = br_sweep::run_br_sweep(
            args.seed, args.steps, &args.preset, 0.2, &br_values,
        );
        for p in &res.points {
            println!(
                "  Br={:.3}: peak|v|={:.3e} CG={:.1} Newton_iter={:.1} conv={:.0}% wallclock {:.3}s",
                p.br, p.peak_v, p.cg_iter_mean, p.newton_iter_mean,
                p.newton_converged_pct, p.wallclock_s,
            );
        }
        Some(res)
    } else {
        None
    };
    let ar_sweep_res: Option<ArSweepResults> = None;
    let bi_sweep_res: Option<BiSweepResults> = None;

    // Run regression BEFORE physics so we can pipe its vmax_peak
    // into the physics report as the reference for
    // `peak_v_damping_ratio`. Each scenario resolves its own
    // previous-step report (auto-detect or CLI override) — the
    // physics run is compared against step{N-1}_physics_report.md,
    // the regression run against step{N-1}_regression_report.md.
    let regression_vmax: Option<f64> = if run_sin {
        let previous_regression = resolve_previous(
            ReportKind::Step4Regression,
            args.compare_to_regression.as_ref(),
        );
        match run_scenario(
            &args, ForceKind::Sinusoidal, ReportKind::Step4Regression,
            "step4_regression_report.md", &scales, previous_regression.as_ref(), &mms,
            None,
            None,
            None,
            None,
        ) {
            Ok(vmax) => Some(vmax),
            Err(e) => {
                eprintln!("regression run failed: {}", e);
                return ExitCode::from(1);
            }
        }
    } else {
        None
    };
    if run_gpe {
        let previous_physics = resolve_previous(
            ReportKind::Step4Physics,
            args.compare_to_physics.as_ref(),
        );
        if let Err(e) = run_scenario(
            &args, ForceKind::Gpe, ReportKind::Step4Physics,
            "step4_physics_report.md", &scales, previous_physics.as_ref(), &mms,
            ar_sweep_res.as_ref(),
            bi_sweep_res.as_ref(),
            br_sweep_res.as_ref(),
            regression_vmax,
        ) {
            eprintln!("physics run failed: {}", e);
            return ExitCode::from(1);
        }
    }

    ExitCode::SUCCESS
}
