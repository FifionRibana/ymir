//! Binary: run the Step 0 baseline scenario and emit a markdown report.
//!
//! Usage:
//!
//! ```bash
//! cargo run --release --bin step_baseline -- \
//!     --seed 42 --grids 64,128 --steps 300 \
//!     --output docs/reports/step0_report.md
//! ```
//!
//! The binary accepts arbitrary seeds so the future T3 stochastic
//! harness can reuse it unchanged.

use std::path::PathBuf;
use std::process::ExitCode;

use ymir_core::tectonics_v2::diagnostics::{
    run_baseline, write_markdown_report, BaselineConfig,
};
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::stokes::StokesConfig;

fn parse_args() -> Result<(u64, Vec<(usize, usize)>, usize, PathBuf), String> {
    let mut seed: u64 = 42;
    let mut grids_str: String = "64,128".into();
    let mut steps: usize = 300;
    let mut output: PathBuf = PathBuf::from("docs/reports/step0_report.md");

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--seed" => {
                i += 1;
                seed = args
                    .get(i)
                    .ok_or("--seed requires a value")?
                    .parse()
                    .map_err(|e| format!("bad --seed: {e}"))?;
            }
            "--grids" => {
                i += 1;
                grids_str = args.get(i).ok_or("--grids requires a value")?.clone();
            }
            "--steps" => {
                i += 1;
                steps = args
                    .get(i)
                    .ok_or("--steps requires a value")?
                    .parse()
                    .map_err(|e| format!("bad --steps: {e}"))?;
            }
            "--output" => {
                i += 1;
                output = PathBuf::from(args.get(i).ok_or("--output requires a value")?);
            }
            "--help" | "-h" => {
                println!(
                    "Usage: step_baseline [--seed N] [--grids N1,N2,...] [--steps N] [--output PATH]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }

    // Grid sizes parsed as "N" for NxN, or "NxM" for rectangular.
    let mut grids: Vec<(usize, usize)> = Vec::new();
    for tok in grids_str.split(',') {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
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
    Ok((seed, grids, steps, output))
}

fn main() -> ExitCode {
    let (seed, grids, steps, output) = match parse_args() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("argument error: {e}");
            return ExitCode::from(2);
        }
    };

    let scales = Scales::default();
    println!("{}", scales.report());

    let mut configs = Vec::new();
    let mut metrics = Vec::new();

    let stokes_cfg = StokesConfig::default();

    // Heightmaps directory sits next to the report.
    let heightmap_dir = output
        .parent()
        .unwrap_or_else(|| std::path::Path::new("."))
        .join("step0_heightmaps");

    for (nx, ny) in &grids {
        println!("-- running baseline at {}×{} for {} steps --", nx, ny, steps);
        let cfg = BaselineConfig {
            seed,
            grid_nx: *nx,
            grid_ny: *ny,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps,
            cfl_factor: 0.3,
            forcing_amplitude: 0.1,
            stokes: stokes_cfg,
            heightmap_fractions: vec![0.0, 0.5, 1.0],
            output_dir: heightmap_dir.clone(),
        };
        let result = run_baseline(&cfg);
        println!(
            "  wallclock: {:.3}s; outer iters mean/max: {:.1}/{}; mass drift: {:.3e}",
            result.metrics.wallclock_total.as_secs_f64(),
            result.metrics.outer_iter_mean,
            result.metrics.outer_iter_max,
            result.metrics.mass_drift_relative,
        );
        configs.push(result.config_dump);
        metrics.push(result.metrics);
    }

    if let Err(e) = write_markdown_report(&output, seed, &scales, &configs, &metrics) {
        eprintln!("failed to write report to {:?}: {}", output, e);
        return ExitCode::from(1);
    }
    println!("report written to {}", output.display());
    ExitCode::SUCCESS
}
