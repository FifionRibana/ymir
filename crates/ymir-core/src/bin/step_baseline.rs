//! Binary: run the Step 1 baseline scenario and emit a markdown
//! report. Also supports re-running the Step 0 anchor when called
//! with a preset that selects the constant-η solver path; by default
//! Step 1 runs the power-law rheology + Newton solver.
//!
//! ```bash
//! cargo run --release --bin step_baseline -- \
//!     --seed 42 --grids 64,128 --steps 300 \
//!     --preset dynamic-accidented --nonlinear-solver newton \
//!     --output docs/reports/step1_report.md \
//!     --compare-to docs/reports/step0_report.md
//! ```

use std::path::PathBuf;
use std::process::ExitCode;

use ymir_core::tectonics_v2::diagnostics::comparison::{parse_step_report, StepReference};
use ymir_core::tectonics_v2::diagnostics::harness::{
    run_baseline, BaselineConfig, NonlinearChoice,
};
use ymir_core::tectonics_v2::diagnostics::report::{write_markdown_report, ReportInputs};
use ymir_core::tectonics_v2::presets::Preset;
use ymir_core::tectonics_v2::scales::Scales;

fn parse_args() -> Result<
    (u64, Vec<(usize, usize)>, usize, PathBuf, Preset, NonlinearChoice, Option<PathBuf>),
    String,
> {
    let mut seed: u64 = 42;
    let mut grids_str: String = "64,128".into();
    let mut steps: usize = 300;
    let mut output: PathBuf = PathBuf::from("docs/reports/step1_report.md");
    let mut preset_name: String = "dynamic-accidented".into();
    let mut nonlinear_str: String = "newton".into();
    let mut compare_to: Option<PathBuf> = None;

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--seed" => {
                i += 1;
                seed = args.get(i).ok_or("--seed requires a value")?
                    .parse().map_err(|e| format!("bad --seed: {e}"))?;
            }
            "--grids" => {
                i += 1;
                grids_str = args.get(i).ok_or("--grids requires a value")?.clone();
            }
            "--steps" => {
                i += 1;
                steps = args.get(i).ok_or("--steps requires a value")?
                    .parse().map_err(|e| format!("bad --steps: {e}"))?;
            }
            "--output" => {
                i += 1;
                output = PathBuf::from(args.get(i).ok_or("--output requires a value")?);
            }
            "--preset" => {
                i += 1;
                preset_name = args.get(i).ok_or("--preset requires a value")?.clone();
            }
            "--nonlinear-solver" => {
                i += 1;
                nonlinear_str = args.get(i).ok_or("--nonlinear-solver requires a value")?.clone();
            }
            "--compare-to" => {
                i += 1;
                compare_to = Some(PathBuf::from(
                    args.get(i).ok_or("--compare-to requires a value")?,
                ));
            }
            "--help" | "-h" => {
                println!(
                    "Usage: step_baseline [--seed N] [--grids N1,N2,...] [--steps N] \
                     [--preset {{dynamic-accidented|stable-shield|soft-planet}}] \
                     [--nonlinear-solver {{newton|picard}}] \
                     [--compare-to PATH] [--output PATH]"
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

    let preset = Preset::by_name(&preset_name)?;
    let nonlinear = NonlinearChoice::parse(&nonlinear_str)?;

    Ok((seed, grids, steps, output, preset, nonlinear, compare_to))
}

fn main() -> ExitCode {
    let (seed, grids, steps, output, preset, nonlinear, compare_to) = match parse_args() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("argument error: {e}");
            return ExitCode::from(2);
        }
    };

    let scales = Scales::default();
    println!("{}", scales.report());
    println!(
        "preset: {} | nonlinear solver: {} | grids: {:?} | steps: {}",
        preset.name, nonlinear.label(), grids, steps,
    );

    let previous: Option<StepReference> = compare_to.as_ref().and_then(|p| {
        match parse_step_report(p) {
            Ok(r) => {
                println!("parsed previous report: {:?} ({} grids)", p, r.grids.len());
                Some(r)
            }
            Err(e) => {
                eprintln!("warning: could not parse {:?}: {}", p, e);
                None
            }
        }
    });

    let heightmap_dir = output
        .parent()
        .unwrap_or_else(|| std::path::Path::new("."))
        .join("step1_heightmaps");

    let mut configs = Vec::new();
    let mut metrics = Vec::new();

    for (nx, ny) in &grids {
        println!("-- running baseline at {}×{} for {} steps --", nx, ny, steps);
        let mut base = BaselineConfig::dynamic_accidented_defaults();
        base.seed = seed;
        base.grid_nx = *nx;
        base.grid_ny = *ny;
        base.steps = steps;
        base.preset = preset.clone();
        base.nonlinear = nonlinear;
        base.output_dir = heightmap_dir.clone();

        let result = run_baseline(&base);
        println!(
            "  wallclock: {:.3}s; CG iters/Newton mean/max: {:.1}/{}; mass drift: {:.3e}",
            result.metrics.wallclock_total.as_secs_f64(),
            result.metrics.cg_iter_mean,
            result.metrics.cg_iter_max,
            result.metrics.mass_drift_relative,
        );
        if let Some(na) = &result.metrics.newton {
            let (c, s, d, cap) = na.outcome_percentages();
            println!(
                "  Newton outcomes — conv {:.1}%, stall {:.1}%, div {:.1}%, cap {:.1}%; η_max/η_min mean {:.1}, max {:.1}",
                c, s, d, cap,
                na.eta_contrast_mean(),
                na.eta_contrast_max(),
            );
        }
        configs.push(result.config_dump);
        metrics.push(result.metrics);
    }

    let justifications: Vec<String> = vec![String::new(); metrics.len()];
    let inputs = ReportInputs {
        seed,
        scales: &scales,
        configs: &configs,
        metrics: &metrics,
        previous: previous.as_ref(),
        suspect_justifications: &justifications,
    };
    if let Err(e) = write_markdown_report(&output, &inputs) {
        eprintln!("failed to write report: {}", e);
        return ExitCode::from(1);
    }
    println!("report written to {}", output.display());
    ExitCode::SUCCESS
}
