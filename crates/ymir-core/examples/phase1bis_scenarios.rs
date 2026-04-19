//! Phase 1-bis scenario runner for issue #75.
//!
//! Runs the tectonic solver on a 64² seed=42 grid for 300 macro steps under
//! one of four configurations (A/B/C/D) and writes tracing diagnostic events
//! to a per-scenario log file. Prints the wallclock elapsed time on stdout as
//! the last line so a driver script can capture it trivially.
//!
//! The four scenarios isolate which extension drives solver cost:
//!
//! - **A**: bare thin-sheet — all extensions off (cratonic, yielding,
//!   slab-pull, mantle flow, dynamic boundaries, recycling, basal friction).
//! - **B**: everything on (current main default, Newton solver).
//! - **C**: default except `boundaries.slab_pull_enabled = false`.
//! - **D**: default except `mantle.enabled = false`.
//!
//! Usage:
//! ```text
//! cargo run --release --example phase1bis_scenarios -- <A|B|C|D> <rep_number> [log_dir]
//! ```
//!
//! The log directory defaults to `logs/`. One file per (scenario, rep) is
//! produced with name `phase1bis_<scenario>_<rep>.log`.

use std::env;
use std::fs::OpenOptions;
use std::path::PathBuf;
use std::sync::Mutex;
use std::time::Instant;

use tracing_subscriber::EnvFilter;

use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::plates::{PlateConfig, generate_plates};
use ymir_core::tectonics::solver::config::{NonlinearSolver, TectonicsConfig};
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::grid::StaggeredGrid;
use ymir_core::tectonics::solver::tectonics::{DynamicPlateContext, run_tectonics};
use ymir_core::tectonics::solver::workspace::SolverWorkspace;

#[derive(Clone, Copy, Debug)]
enum Scenario {
    A,
    B,
    C,
    D,
}

impl Scenario {
    fn from_arg(s: &str) -> Option<Self> {
        match s.to_ascii_uppercase().as_str() {
            "A" => Some(Self::A),
            "B" => Some(Self::B),
            "C" => Some(Self::C),
            "D" => Some(Self::D),
            _ => None,
        }
    }

    fn tag(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::B => "B",
            Self::C => "C",
            Self::D => "D",
        }
    }
}

fn build_config(scenario: Scenario, num_steps: usize) -> TectonicsConfig {
    // Start from the full-featured default, then strip per scenario.
    let mut cfg = TectonicsConfig::default();

    // Harmonize with the regime used in the Phase 1 reference logs
    // (Newton solver, adaptive dt on, dt_target = 2.0).
    cfg.num_timesteps = num_steps;
    cfg.nonlinear_solver = NonlinearSolver::Newton;

    match scenario {
        Scenario::A => {
            // Bare thin-sheet: every extension off. `boundaries.enabled =
            // false` skips the whole source-rate pipeline, which also
            // makes `slab_pull_enabled` a no-op.
            cfg.boundaries.enabled = false;
            cfg.boundaries.slab_pull_enabled = false;
            cfg.cratonic.enabled = false;
            cfg.cratonic.max_factor = 1.0;
            cfg.yielding.enabled = false;
            cfg.mantle.enabled = false;
            cfg.dynamic_boundaries = false;
            cfg.recycling.enabled = false;
            cfg.basal_friction = 0.0;
        }
        Scenario::B => {
            // Everything on — nothing to change.
        }
        Scenario::C => {
            cfg.boundaries.slab_pull_enabled = false;
        }
        Scenario::D => {
            cfg.mantle.enabled = false;
        }
    }

    cfg
}

fn install_subscriber(log_path: &PathBuf) {
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(log_path)
        .expect("failed to open log file");

    // Enable DEBUG on the diagnostic targets, INFO on everything else.
    // Using an EnvFilter lets the caller override via RUST_LOG.
    let default_filter = [
        "info",
        "ymir_core=info",
        "rhs_breakdown=debug",
        "eta_breakdown=debug",
        "residual_spatial=debug",
        "phase_timings=info",
    ]
    .join(",");
    let filter =
        EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(default_filter));

    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(Mutex::new(file))
        .with_ansi(false)
        .with_target(true)
        .init();
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <A|B|C|D> <rep_number> [log_dir]",
            args.first().map(String::as_str).unwrap_or("phase1bis_scenarios")
        );
        std::process::exit(2);
    }

    let scenario = Scenario::from_arg(&args[1]).unwrap_or_else(|| {
        eprintln!("unknown scenario `{}` (expected A, B, C, or D)", args[1]);
        std::process::exit(2);
    });
    let rep: usize = args[2].parse().unwrap_or_else(|_| {
        eprintln!("rep_number must be a non-negative integer, got `{}`", args[2]);
        std::process::exit(2);
    });

    let log_dir = args.get(3).cloned().unwrap_or_else(|| "logs".to_string());
    // Optional 4th arg: number of macro steps (Phase 1-bis uses 300 for
    // the full sweep; Phase 2 validation uses 120 for the reduced-scope
    // rerun described in issue #75).
    let num_steps: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(300);
    std::fs::create_dir_all(&log_dir).expect("failed to create log dir");
    let log_path: PathBuf =
        PathBuf::from(&log_dir).join(format!("phase1bis_{}_{:02}.log", scenario.tag(), rep));

    install_subscriber(&log_path);

    // Seed = 42, grid = 64 × 64, defaults otherwise — this is the
    // reference regime from the Phase 1 report.
    let plate_config = PlateConfig { grid_width: 64, grid_height: 64, ..Default::default() };
    let seed = WorldSeed::new(42);
    let init = generate_plates(&plate_config, &seed);

    let nx = init.grid_width;
    let ny = init.grid_height;
    let dx = 1.0 / nx as f64;
    let mut grid = StaggeredGrid::new(nx, ny, dx);
    for j in 0..ny {
        for i in 0..nx {
            grid.s.set(i, j, init.thickness.data[j * nx + i] as f64);
        }
    }

    let traction = init.to_traction_field();
    let num_plates = init.plates.len();
    let mut plate_ctx = DynamicPlateContext {
        ids: init.plate_ids.clone(),
        plates: init.plates.clone(),
        traction,
        next_id: num_plates,
        disp_x: Field2D::new(nx, ny),
        disp_y: Field2D::new(nx, ny),
    };

    let config = build_config(scenario, num_steps);
    grid.basal_friction = config.basal_friction;

    let mut ws = SolverWorkspace::new(nx, ny);

    eprintln!(
        "[phase1bis] scenario={} rep={:02} grid={}x{} steps={} solver={:?} log={}",
        scenario.tag(),
        rep,
        nx,
        ny,
        config.num_timesteps,
        config.nonlinear_solver,
        log_path.display()
    );

    let start = Instant::now();
    let result = run_tectonics(&config, &mut plate_ctx, &mut grid, &mut ws, |_, _, _, _| true);
    let elapsed = start.elapsed();

    match result {
        Ok(()) => {
            // Single machine-parseable summary line. The driver script
            // greps for `SUMMARY` to aggregate wallclock stats.
            println!(
                "SUMMARY scenario={} rep={:02} elapsed_ms={} steps={} outcome=ok",
                scenario.tag(),
                rep,
                elapsed.as_millis(),
                config.num_timesteps
            );
        }
        Err(e) => {
            println!(
                "SUMMARY scenario={} rep={:02} elapsed_ms={} steps={} outcome=err reason={}",
                scenario.tag(),
                rep,
                elapsed.as_millis(),
                config.num_timesteps,
                e
            );
        }
    }
}
