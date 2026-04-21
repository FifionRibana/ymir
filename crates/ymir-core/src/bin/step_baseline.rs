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

use ymir_core::tectonics_v2::diagnostics::comparison::{parse_step_report, StepReference};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::diagnostics::mms_bench;
use ymir_core::tectonics_v2::diagnostics::report::{
    write_markdown_report, ReportInputs, ReportKind,
};
use ymir_core::tectonics_v2::presets::Preset;
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
    compare_to: Option<PathBuf>,
    forcing: ForcingSelection,
    sinusoidal_amplitude: f64,
}

fn parse_args() -> Result<Args, String> {
    let mut seed: u64 = 42;
    let mut grids_str: String = "64,128".into();
    let mut steps: usize = 300;
    let mut output_dir: PathBuf = PathBuf::from("docs/reports/");
    let mut preset_name: String = "dynamic-accidented".into();
    let mut nonlinear_str: String = "newton".into();
    let mut compare_to: Option<PathBuf> = None;
    let mut forcing_str: String = "both".into();
    let mut sin_amp: f64 = 10.0;

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
            "--forcing" => { i += 1; forcing_str = args[i].clone(); }
            "--sinusoidal-amplitude" => { i += 1; sin_amp = args[i].parse().map_err(|e| format!("bad --sinusoidal-amplitude: {e}"))?; }
            "--help" | "-h" => {
                println!(
                    "Usage: step_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset NAME] [--nonlinear-solver newton|picard] \
                     [--forcing gpe|sinusoidal|both] [--sinusoidal-amplitude F] \
                     [--compare-to PATH] [--output-dir PATH]"
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

    Ok(Args {
        seed, grids, steps, output_dir,
        preset: Preset::by_name(&preset_name)?,
        nonlinear: NonlinearChoice::parse(&nonlinear_str)?,
        compare_to,
        forcing: ForcingSelection::parse(&forcing_str)?,
        sinusoidal_amplitude: sin_amp,
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
) -> Result<(), String> {
    let heightmap_subdir = args.output_dir.join(match kind {
        ForceKind::Gpe => "step2_physics_heightmaps",
        ForceKind::Sinusoidal => "step2_regression_heightmaps",
    });

    let mut configs = Vec::new();
    let mut metrics = Vec::new();

    println!(
        "\n=== scenario: {} (writing {}) ===",
        kind.label(),
        report_file,
    );

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
    };
    write_markdown_report(&output, &inputs)
        .map_err(|e| format!("failed to write report {:?}: {}", output, e))?;
    println!("report written to {}", output.display());
    Ok(())
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

    let previous: Option<StepReference> = args.compare_to.as_ref().and_then(|p| {
        match parse_step_report(p) {
            Ok(r) => {
                println!("parsed previous report: {:?} ({} grids)", p, r.grids.len());
                Some(r)
            }
            Err(e) => { eprintln!("warning: could not parse {:?}: {}", p, e); None }
        }
    });

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

    if run_gpe {
        if let Err(e) = run_scenario(
            &args, ForceKind::Gpe, ReportKind::Step2Physics,
            "step2_physics_report.md", &scales, previous.as_ref(), &mms,
        ) {
            eprintln!("physics run failed: {}", e);
            return ExitCode::from(1);
        }
    }
    if run_sin {
        if let Err(e) = run_scenario(
            &args, ForceKind::Sinusoidal, ReportKind::Step2Regression,
            "step2_regression_report.md", &scales, previous.as_ref(), &mms,
        ) {
            eprintln!("regression run failed: {}", e);
            return ExitCode::from(1);
        }
    }

    ExitCode::SUCCESS
}
