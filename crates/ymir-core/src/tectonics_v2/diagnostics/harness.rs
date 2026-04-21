//! Baseline runner for Step 1.
//!
//! Runs the coupled nonlinear-Stokes + advection loop with a
//! configurable [`Preset`], a startup continuation ramp, and a
//! choice of nonlinear solver (Newton or Picard). Aggregates the
//! same metrics as Step 0 plus Newton-specific statistics.

use std::path::{Path, PathBuf};
use std::time::Instant;

use super::metrics::{
    condition_number_estimate, IterationHistogram, Metrics, SolverConfigDump,
};
use super::newton_metrics::{cap_activation_fraction, eta_contrast, NewtonAggregate};
use crate::grid::GridF32;
use crate::tectonics_v2::advection::{cfl_dt, integrated_mass, step_upwind};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::forcing::{sample_to_faces, SinusoidalForce};
use crate::tectonics_v2::presets::Preset;
use crate::tectonics_v2::rheology::{self, StrainRate, ViscosityLaw};
use crate::tectonics_v2::stokes::continuation::run_continuation;
use crate::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearOutcome, NonlinearSolver,
};
use crate::tectonics_v2::stokes::picard::{PicardConfig, PicardSolver};
use crate::tectonics_v2::stokes::solver::ConjugateGradient;
use crate::tectonics_v2::stokes::Grid;

/// Which nonlinear solver to drive the main loop with.
#[derive(Clone, Copy, Debug)]
pub enum NonlinearChoice {
    Newton,
    Picard,
}

impl NonlinearChoice {
    pub fn label(&self) -> &'static str {
        match self {
            NonlinearChoice::Newton => "newton",
            NonlinearChoice::Picard => "picard",
        }
    }
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "newton" => Ok(NonlinearChoice::Newton),
            "picard" => Ok(NonlinearChoice::Picard),
            other => Err(format!("unknown nonlinear solver '{}'", other)),
        }
    }
}

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
    pub preset: Preset,
    pub nonlinear: NonlinearChoice,
    pub newton_cfg: NewtonConfig,
    pub picard_cfg: PicardConfig,
    pub heightmap_fractions: Vec<f64>,
    pub output_dir: PathBuf,
}

impl BaselineConfig {
    pub fn dynamic_accidented_defaults() -> Self {
        Self {
            seed: 42,
            grid_nx: 64,
            grid_ny: 64,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps: 300,
            cfl_factor: 0.3,
            forcing_amplitude: 10.0,
            preset: Preset::dynamic_accidented(),
            nonlinear: NonlinearChoice::Newton,
            newton_cfg: NewtonConfig::default(),
            picard_cfg: PicardConfig::default(),
            heightmap_fractions: vec![0.0, 0.5, 1.0],
            output_dir: PathBuf::from("docs/reports/step1_heightmaps"),
        }
    }
}

#[derive(Clone, Debug)]
pub struct BaselineResult {
    pub metrics: Metrics,
    pub config_dump: SolverConfigDump,
}

fn init_thickness(nx: usize, ny: usize, seed: u64) -> Field2D {
    use std::f64::consts::PI;
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

fn solve_nonlinear(
    grid: &Grid,
    law: &ViscosityLaw,
    rhs_x: &[f64],
    rhs_y: &[f64],
    vx: &mut [f64],
    vy: &mut [f64],
    choice: NonlinearChoice,
    newton_cfg: NewtonConfig,
    picard_cfg: PicardConfig,
    cg: &ConjugateGradient,
) -> NonlinearOutcome {
    match choice {
        NonlinearChoice::Newton => {
            let solver = NewtonSolver::new(newton_cfg);
            solver.solve(grid, law, rhs_x, rhs_y, vx, vy, cg)
        }
        NonlinearChoice::Picard => {
            let solver = PicardSolver::new(picard_cfg);
            solver.solve(grid, law, rhs_x, rhs_y, vx, vy, cg)
        }
    }
}

/// Drive a single baseline run.
pub fn run_baseline(cfg: &BaselineConfig) -> BaselineResult {
    let nx = cfg.grid_nx;
    let ny = cfg.grid_ny;
    let dx = cfg.domain_lx / nx as f64;
    let dy = cfg.domain_ly / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);

    let mut s = init_thickness(nx, ny, cfg.seed);
    let mut s_next = Field2D::new(nx, ny);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    let force = SinusoidalForce::new(cfg.forcing_amplitude, cfg.domain_lx);

    let mass_initial = integrated_mass(&s);

    std::fs::create_dir_all(&cfg.output_dir).ok();
    let capture_steps: Vec<usize> = cfg
        .heightmap_fractions
        .iter()
        .map(|f| (f.clamp(0.0, 1.0) * cfg.steps as f64).round() as usize)
        .collect();
    let mut heightmap_paths: Vec<String> = Vec::new();

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut newton_agg = NewtonAggregate::default();
    let mut max_abs_mean_vx = 0.0f64;
    let mut max_abs_mean_vy = 0.0f64;
    let mut vmax_peak = 0.0f64;

    // Final rheology (after continuation).
    let law_final = ViscosityLaw {
        n: cfg.preset.rheology.n,
        b_prefactor: cfg.preset.rheology.b_prefactor,
        strain_rate_floor: cfg.preset.rheology.strain_rate_floor,
        eta_max_cap: cfg.preset.rheology.eta_max_cap,
        k_saturation: cfg.preset.rheology.k_saturation,
    };

    let cg = ConjugateGradient::new(cfg.newton_cfg.linear_tol, cfg.newton_cfg.linear_max_iter);

    let start = Instant::now();

    if capture_steps.contains(&0) {
        let path = cfg.output_dir.join(format!("s_{}x{}_t0000.png", nx, ny));
        if save_heightmap(&s, &path).is_ok() {
            heightmap_paths.push(path.display().to_string().replace('\\', "/"));
        }
    }

    // --- Startup continuation (t = 0 only) ---
    sample_to_faces(&force, nx, ny, dx, dy, &s, &mut fx, &mut fy);
    let newton_solver = NewtonSolver::new(cfg.newton_cfg);
    let cont = run_continuation(
        &grid,
        &law_final,
        &cfg.preset.continuation,
        fx.data(),
        fy.data(),
        &mut vx,
        &mut vy,
        &newton_solver,
        &cg,
    );
    newton_agg.continuation_all_converged = Some(cont.all_converged);
    newton_agg.continuation_iters_used = cont.sub_outcomes.len() as u32;
    // Record each sub-solve outcome and its CG iters.
    for (_n, oc) in &cont.sub_outcomes {
        match oc {
            NonlinearOutcome::Converged { outer_iters, trace, .. } => {
                newton_agg.converged += 1;
                newton_agg.outer_iters.push(*outer_iters);
                for &cg_iters in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(cg_iters);
                }
            }
            NonlinearOutcome::Stalled { outer_iters, trace } => {
                newton_agg.stalled += 1;
                newton_agg.outer_iters.push(*outer_iters);
                for &cg_iters in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(cg_iters);
                }
            }
            NonlinearOutcome::Diverged { outer_iters, trace, .. } => {
                newton_agg.diverged += 1;
                newton_agg.outer_iters.push(*outer_iters);
                for &cg_iters in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(cg_iters);
                }
            }
            NonlinearOutcome::CappedIters { max_iters_hit, trace, .. } => {
                newton_agg.capped += 1;
                newton_agg.outer_iters.push(*max_iters_hit);
                for &cg_iters in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(cg_iters);
                }
            }
        }
    }

    // Track cap-activation during the ramp separately.
    let sr_after_ramp = StrainRate::compute(
        nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy,
    );
    let eta_after_ramp = rheology::build_eta_field(&law_final, &sr_after_ramp.eps_ii_center);
    newton_agg.cap_fraction_ramp_max = newton_agg
        .cap_fraction_ramp_max
        .max(cap_activation_fraction(&eta_after_ramp, law_final.eta_max_cap));

    // --- Steady-state loop ---
    for step in 0..cfg.steps {
        sample_to_faces(&force, nx, ny, dx, dy, &s, &mut fx, &mut fy);
        let outcome = solve_nonlinear(
            &grid,
            &law_final,
            fx.data(),
            fy.data(),
            &mut vx,
            &mut vy,
            cfg.nonlinear,
            cfg.newton_cfg,
            cfg.picard_cfg,
            &cg,
        );

        match &outcome {
            NonlinearOutcome::Converged { outer_iters, trace, .. } => {
                newton_agg.converged += 1;
                newton_agg.outer_iters.push(*outer_iters);
                for &c in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(c);
                }
            }
            NonlinearOutcome::Stalled { outer_iters, trace } => {
                newton_agg.stalled += 1;
                newton_agg.outer_iters.push(*outer_iters);
                for &c in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(c);
                }
            }
            NonlinearOutcome::Diverged { outer_iters, trace, .. } => {
                newton_agg.diverged += 1;
                newton_agg.outer_iters.push(*outer_iters);
                for &c in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(c);
                }
            }
            NonlinearOutcome::CappedIters { max_iters_hit, trace, .. } => {
                newton_agg.capped += 1;
                newton_agg.outer_iters.push(*max_iters_hit);
                for &c in &trace.linear_iters {
                    newton_agg.cg_iters_per_newton_step.push(c);
                }
            }
        }

        // Cap-activation fraction and η contrast at steady state.
        let sr = StrainRate::compute(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy);
        let eta_cc = rheology::build_eta_field(&law_final, &sr.eps_ii_center);
        newton_agg.cap_fraction_steady_max = newton_agg
            .cap_fraction_steady_max
            .max(cap_activation_fraction(&eta_cc, law_final.eta_max_cap));
        newton_agg.eta_contrast_samples.push(eta_contrast(&eta_cc));

        // Null-space health.
        let m_vx: f64 = vx.iter().sum::<f64>() / vx.len() as f64;
        let m_vy: f64 = vy.iter().sum::<f64>() / vy.len() as f64;
        max_abs_mean_vx = max_abs_mean_vx.max(m_vx.abs());
        max_abs_mean_vy = max_abs_mean_vy.max(m_vy.abs());

        let vmax_step = vx
            .iter()
            .chain(vy.iter())
            .fold(0.0_f64, |acc, &v| acc.max(v.abs()));
        vmax_peak = vmax_peak.max(vmax_step);

        // Advection.
        let dt = cfl_dt(dx, dy, &vx, &vy, cfg.cfl_factor);
        if dt.is_finite() {
            step_upwind(nx, ny, dx, dy, dt, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
            std::mem::swap(&mut s, &mut s_next);
        }

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

    let cg_iter_mean = newton_agg.cg_iters_per_newton_mean();
    let cg_iter_max = newton_agg.cg_iters_per_newton_max();
    let kappa = condition_number_estimate(
        cg_iter_mean.round() as usize,
        cfg.newton_cfg.linear_tol,
    );

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
        eta_contrast: newton_agg.eta_contrast_mean(),
        cg_iter_mean,
        cg_iter_max,
        cg_iter_histogram: IterationHistogram::from_samples(
            &newton_agg.cg_iters_per_newton_step,
        ),
        mass_s_initial: mass_initial,
        mass_s_final: mass_final,
        mass_drift_relative: drift,
        max_abs_mean_vx,
        max_abs_mean_vy,
        vmax_peak,
        heightmap_paths,
        newton: Some(newton_agg),
        s_eq: None,
        boundary_type_diversity: None,
        yielding_cell_fraction: None,
        cratonic_stability: None,
        newton_outcome_distribution: None,
        age_field_stats: None,
    };

    let continuation_str = format!("{:?}", cfg.preset.continuation.n_steps);
    let config_dump = SolverConfigDump {
        formulation: "thin viscous sheet (elliptic, no pressure) with power-law rheology"
            .into(),
        discretization: "MAC staggered (v face / η S cell-centre / ε̇_xy corner)".into(),
        eta_averaging: "arithmetic 4-point at corners (see operator.rs)".into(),
        preconditioner: "velocity Jacobi (Picard-block diagonal), null-space wrapped".into(),
        gauge_fixing: "mean(vx), mean(vy) projected before & after every M⁻¹ + post-solve"
            .into(),
        cg_tol: cfg.newton_cfg.linear_tol,
        cg_max_iter: cfg.newton_cfg.linear_max_iter,
        cfl_factor: cfg.cfl_factor,
        grid_spacing_nondim: dx,
        body_force: format!("SinusoidalForce(ε={})", cfg.forcing_amplitude),
        seed: cfg.seed,
        preset_name: cfg.preset.name.clone(),
        nonlinear_solver: cfg.nonlinear.label().into(),
        rheology_n: law_final.n,
        strain_rate_floor: law_final.strain_rate_floor,
        eta_max_cap: law_final.eta_max_cap,
        continuation_schedule: continuation_str,
        newton_rel_tol: cfg.newton_cfg.rel_tol,
        newton_max_outer_iters: cfg.newton_cfg.max_outer_iters,
    };

    BaselineResult { metrics, config_dump }
}
