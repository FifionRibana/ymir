//! Baseline runner used by Step 0's `step_baseline` binary and the
//! smoke-test integration test.
//!
//! The harness sets up a clean, reproducible scenario (constant
//! η̃ = 1, mildly perturbed S̃, sinusoidal body-force placeholder)
//! and runs the coupled Stokes + advection loop for the configured
//! number of macro steps. Metrics are aggregated step-by-step.

use std::path::{Path, PathBuf};
use std::time::Instant;

use super::metrics::{
    condition_number_estimate, IterationHistogram, Metrics, SolverConfigDump,
};
use crate::grid::GridF32;
use crate::tectonics_v2::advection::{cfl_dt, integrated_mass, step_upwind};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::forcing::{sample_to_faces, SinusoidalForce};
use crate::tectonics_v2::stokes::{solve_stokes, Grid, StokesConfig};

/// Configuration for a baseline run.
#[derive(Clone, Debug)]
pub struct BaselineConfig {
    pub seed: u64,
    pub grid_nx: usize,
    pub grid_ny: usize,
    pub domain_lx: f64,
    pub domain_ly: f64,
    pub steps: usize,
    pub cfl_factor: f64,
    pub forcing_amplitude: f64,
    pub stokes: StokesConfig,
    /// Capture heightmaps at these fractional times (0.0 ≤ f ≤ 1.0).
    pub heightmap_fractions: Vec<f64>,
    /// Output directory for heightmaps. The harness appends filenames.
    pub output_dir: PathBuf,
}

impl Default for BaselineConfig {
    fn default() -> Self {
        Self {
            seed: 42,
            grid_nx: 64,
            grid_ny: 64,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps: 300,
            cfl_factor: 0.3,
            forcing_amplitude: 0.1,
            stokes: StokesConfig::default(),
            heightmap_fractions: vec![0.0, 0.5, 1.0],
            output_dir: PathBuf::from("docs/reports/step0_heightmaps"),
        }
    }
}

/// Return of a single baseline run.
#[derive(Clone, Debug)]
pub struct BaselineResult {
    pub metrics: Metrics,
    pub config_dump: SolverConfigDump,
}

/// Initialise S̃ with a low-amplitude spatial perturbation so the
/// advection is observable.
fn init_thickness(nx: usize, ny: usize, seed: u64) -> Field2D {
    use std::f64::consts::PI;
    // Deterministic, seed-aware but simple: a fixed spatial pattern
    // scaled by a small amplitude. The seed shifts the phase so two
    // seeds produce distinguishable initial conditions without
    // requiring an RNG dependency in Step 0.
    let phase = ((seed.wrapping_mul(2654435761u64)) as f64) / (u64::MAX as f64) * 2.0 * PI;
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let bump = 1.0 + 0.02 * ((2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos());
            s.set(i, j, bump);
        }
    }
    s
}

fn save_heightmap(s: &Field2D, path: &Path) -> Result<(), String> {
    let nx = s.nx();
    let ny = s.ny();
    let data: Vec<f32> = s.data().iter().map(|&v| v as f32).collect();
    let grid = GridF32::from_vec(nx, ny, data);
    grid.save_png_u16(path)
}

/// Run a single baseline configuration and return aggregated metrics.
pub fn run_baseline(cfg: &BaselineConfig) -> BaselineResult {
    let nx = cfg.grid_nx;
    let ny = cfg.grid_ny;
    let dx = cfg.domain_lx / nx as f64;
    let dy = cfg.domain_ly / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);

    let mut s = init_thickness(nx, ny, cfg.seed);
    let mut s_next = Field2D::new(nx, ny);
    let eta = Field2D::filled(nx, ny, 1.0);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut p = vec![0.0; nx * ny];
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    let force = SinusoidalForce::new(cfg.forcing_amplitude, cfg.domain_ly);

    let mass_initial = integrated_mass(&s);

    // Heightmap capture schedule.
    std::fs::create_dir_all(&cfg.output_dir).ok();
    let capture_steps: Vec<usize> = cfg
        .heightmap_fractions
        .iter()
        .map(|f| (f.clamp(0.0, 1.0) * cfg.steps as f64).round() as usize)
        .collect();
    let mut heightmap_paths: Vec<String> = Vec::new();

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut outer_iters: Vec<usize> = Vec::with_capacity(cfg.steps);
    let mut inner_iters: Vec<usize> = Vec::with_capacity(cfg.steps);
    let mut inner_max = 0usize;
    let mut max_abs_mean_p = 0.0f64;
    let mut max_abs_mean_vx = 0.0f64;
    let mut max_abs_mean_vy = 0.0f64;
    let mut vmax_peak = 0.0f64;

    let start = Instant::now();

    // Capture at t=0.
    if capture_steps.contains(&0) {
        let path = cfg.output_dir.join(format!("s_{}x{}_t0000.png", nx, ny));
        if save_heightmap(&s, &path).is_ok() {
            heightmap_paths.push(path.display().to_string().replace('\\', "/"));
        }
    }

    for step in 0..cfg.steps {
        // --- Stokes solve --------------------------------------------
        sample_to_faces(&force, nx, ny, dx, dy, &s, &mut fx, &mut fy);
        let stats = solve_stokes(
            &grid,
            &eta,
            fx.data(),
            fy.data(),
            &mut vx,
            &mut vy,
            &mut p,
            &cfg.stokes,
        );
        outer_iters.push(stats.outer_iterations);
        if stats.inner_solves > 0 {
            // Per-solve mean rounded to the nearest integer for histogramming.
            let mean_inner = stats.inner_iterations_total as f64 / stats.inner_solves as f64;
            inner_iters.push(mean_inner.round() as usize);
        }
        inner_max = inner_max.max(stats.inner_iterations_max);
        max_abs_mean_p = max_abs_mean_p.max(stats.mean_p_after.abs());
        max_abs_mean_vx = max_abs_mean_vx.max(stats.mean_vx_after.abs());
        max_abs_mean_vy = max_abs_mean_vy.max(stats.mean_vy_after.abs());

        let vmax_step = vx
            .iter()
            .chain(vy.iter())
            .fold(0.0_f64, |acc, &v| acc.max(v.abs()));
        vmax_peak = vmax_peak.max(vmax_step);

        // --- Advection step ----------------------------------------
        let dt = cfl_dt(dx, dy, &vx, &vy, cfg.cfl_factor);
        // Guard against a stationary-solution infinity: in that case skip advection.
        if dt.is_finite() {
            step_upwind(
                nx, ny, dx, dy, dt, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next,
            );
            std::mem::swap(&mut s, &mut s_next);
        }

        // Capture heightmaps at the scheduled step indices (step+1 because we've just completed a step).
        let completed = step + 1;
        if capture_steps.contains(&completed) {
            let path = cfg
                .output_dir
                .join(format!("s_{}x{}_t{:04}.png", nx, ny, completed));
            if save_heightmap(&s, &path).is_ok() {
                heightmap_paths.push(path.display().to_string().replace('\\', "/"));
            }
        }
    }

    let wallclock = start.elapsed();
    let mass_final = integrated_mass(&s);
    let drift = (mass_final - mass_initial) / mass_initial.abs().max(1.0);

    // κ estimate from mean outer iteration count (a single scalar is
    // more robust than any individual iter which can be an outlier).
    let iter_mean = if outer_iters.is_empty() {
        0.0
    } else {
        outer_iters.iter().sum::<usize>() as f64 / outer_iters.len() as f64
    };
    let kappa = condition_number_estimate(iter_mean.round() as usize, cfg.stokes.outer_tol);

    let outer_max = outer_iters.iter().copied().max().unwrap_or(0);
    let inner_mean = if inner_iters.is_empty() {
        0.0
    } else {
        inner_iters.iter().sum::<usize>() as f64 / inner_iters.len() as f64
    };

    let metrics = Metrics {
        grid_nx: nx,
        grid_ny: ny,
        steps: cfg.steps,
        wallclock_total: wallclock,
        wallclock_per_step_mean: if cfg.steps > 0 {
            wallclock / cfg.steps as u32
        } else {
            wallclock
        },
        kappa_estimate: kappa,
        eta_contrast: 1.0,
        outer_iter_mean: iter_mean,
        outer_iter_max: outer_max,
        outer_iter_histogram: IterationHistogram::from_samples(&outer_iters),
        inner_iter_mean: inner_mean,
        inner_iter_max: inner_max,
        mass_s_initial: mass_initial,
        mass_s_final: mass_final,
        mass_drift_relative: drift,
        max_abs_mean_p,
        max_abs_mean_vx,
        max_abs_mean_vy,
        vmax_peak,
        heightmap_paths,
        s_eq: None,
        boundary_type_diversity: None,
        yielding_cell_fraction: None,
        cratonic_stability: None,
        newton_outcome_distribution: None,
        age_field_stats: None,
    };

    let config_dump = SolverConfigDump {
        discretization: "MAC staggered (v face / P η S cell-centre)".into(),
        harmonic_averaging: "harmonic 4-point for η at corners".into(),
        preconditioner: "block-diag Jacobi (v) + diag(1/η) mass (P), null-space wrapped".into(),
        gauge_fixing: "mean(P), mean(vx), mean(vy) projected before & after every M^-1 and once post-solve".into(),
        outer_tol: cfg.stokes.outer_tol,
        inner_tol: cfg.stokes.inner_tol,
        outer_max_iter: cfg.stokes.outer_max_iter,
        inner_max_iter: cfg.stokes.inner_max_iter,
        cfl_factor: cfg.cfl_factor,
        grid_spacing_nondim: dx,
        body_force: format!("SinusoidalForce(ε={})", cfg.forcing_amplitude),
        seed: cfg.seed,
    };

    BaselineResult { metrics, config_dump }
}
