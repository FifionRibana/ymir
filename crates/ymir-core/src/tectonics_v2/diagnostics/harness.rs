//! Baseline runner.
//!
//! A run is parametrised by a [`Preset`], a [`NonlinearChoice`] and a
//! **body force** supplied as a `Box<dyn BodyForce>`. Step 2 handles
//! both the GPE physics run and the Sinusoidal regression run through
//! the same harness — the only thing that differs between them is
//! which force is passed in.

use std::path::{Path, PathBuf};
use std::time::Instant;

use super::heightmap::{save_heightmap, HeightmapMetadata};
use super::metrics::{
    condition_number_estimate, IterationHistogram, Metrics, SolverConfigDump,
};
use super::newton_metrics::{
    cap_activation_fraction, eta_contrast, yielding_cell_fraction, yielding_intensity,
    NewtonAggregate,
};
use crate::tectonics_v2::advection::{cfl_dt, integrated_mass, step_upwind};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::forcing::{BodyForce, ForceSum, GpeForce, SimulationState, SinusoidalForce, VectorField};
use crate::tectonics_v2::presets::{Preset, YieldingConfig};
use crate::tectonics_v2::rheology::{self, StrainRate, ViscosityLaw};
use crate::tectonics_v2::scales::Scales;
use crate::tectonics_v2::stokes::continuation::run_continuation;
use crate::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearOutcome, NonlinearSolver,
};
use crate::tectonics_v2::stokes::picard::{PicardConfig, PicardSolver};
use crate::tectonics_v2::stokes::solver::ConjugateGradient;
use crate::tectonics_v2::stokes::Grid;

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

/// Which forcing scenario to run. The binary picks these and
/// constructs the corresponding `Box<dyn BodyForce>` to pass into
/// [`BaselineConfig::force`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ForceKind {
    Gpe,
    Sinusoidal,
}

impl ForceKind {
    pub fn label(&self) -> &'static str {
        match self {
            ForceKind::Gpe => "gpe",
            ForceKind::Sinusoidal => "sinusoidal",
        }
    }
}

/// Build a `Box<dyn BodyForce>` for the chosen scenario.
///
/// This centralises the default amplitudes/Ar choices so the binary
/// and the tests agree on what "gpe" or "sinusoidal" mean for a
/// baseline run.
pub fn build_force(kind: ForceKind, scales: &Scales, sin_amplitude: f64, domain_lx: f64) -> Box<dyn BodyForce> {
    match kind {
        ForceKind::Gpe => {
            let mut sum = ForceSum::new();
            sum.push(Box::new(GpeForce::from_scales(scales)));
            Box::new(sum)
        }
        ForceKind::Sinusoidal => {
            let mut sum = ForceSum::new();
            sum.push(Box::new(SinusoidalForce::new(sin_amplitude, domain_lx)));
            Box::new(sum)
        }
    }
}

pub struct BaselineConfig {
    pub seed: u64,
    pub grid_nx: usize,
    pub grid_ny: usize,
    pub domain_lx: f64,
    pub domain_ly: f64,
    pub steps: usize,
    pub cfl_factor: f64,
    /// Total simulated nondim time. The harness uses
    /// `dt_target = total_time / steps` as the intended timestep and
    /// clamps the CFL-derived dt to it when CFL would otherwise give
    /// a value orders of magnitude larger. Default: `6.0` (= 6·τ̃*).
    pub total_time_nondim: f64,
    pub preset: Preset,
    pub nonlinear: NonlinearChoice,
    pub newton_cfg: NewtonConfig,
    pub picard_cfg: PicardConfig,
    pub heightmap_fractions: Vec<f64>,
    pub output_dir: PathBuf,
    /// The body-force term driving the momentum RHS. Wrap composite
    /// forces in a `ForceSum`; the harness reads term names through
    /// the trait so the report can dump whatever the caller supplied.
    pub force: Box<dyn BodyForce>,
    pub force_kind: ForceKind,
    /// Amplitude for the Sinusoidal regression scenario (carried
    /// through to the config dump so the report stays
    /// self-descriptive).
    pub sinusoidal_amplitude: f64,
    /// Initial `S̃` perturbation amplitude around the `1.0` mean.
    /// Step 0/1 used `0.02`; Step 2 bumps the physics scenario to
    /// `0.2` so GPE drives an observable response on the 300-step
    /// window even at the honest thin-sheet `Ar = S*/L* = 0.1`. The
    /// regression scenario keeps `0.02` to preserve the mirror-of-
    /// Step-1 contract.
    pub s_perturbation_amplitude: f64,
    /// Plastic-yielding configuration (Step 3). Default is
    /// [`YieldingConfig::Disabled`] so Step 0/1/2 tests keep their
    /// pre-Step-3 behaviour without edits. The physics baseline and
    /// the Bi sweep pass `YieldingConfig::Enabled(..)` explicitly.
    pub yielding: YieldingConfig,
}

impl BaselineConfig {
    pub fn dynamic_accidented_defaults(scales: &Scales) -> Self {
        let force_kind = ForceKind::Gpe;
        let domain_lx = 1.0;
        let sin_amplitude = 10.0;
        Self {
            seed: 42,
            grid_nx: 64,
            grid_ny: 64,
            domain_lx,
            domain_ly: 1.0,
            steps: 300,
            cfl_factor: 0.3,
            total_time_nondim: 6.0,
            preset: Preset::dynamic_accidented(),
            nonlinear: NonlinearChoice::Newton,
            newton_cfg: NewtonConfig::default(),
            picard_cfg: PicardConfig::default(),
            heightmap_fractions: vec![0.0, 0.5, 1.0],
            output_dir: PathBuf::from("docs/reports/step2_heightmaps"),
            force: build_force(force_kind, scales, sin_amplitude, domain_lx),
            force_kind,
            sinusoidal_amplitude: sin_amplitude,
            s_perturbation_amplitude: 0.2,
            yielding: YieldingConfig::Disabled,
        }
    }
}

#[derive(Clone, Debug)]
pub struct BaselineResult {
    pub metrics: Metrics,
    pub config_dump: SolverConfigDump,
}

fn init_thickness(nx: usize, ny: usize, seed: u64, amplitude: f64) -> Field2D {
    use std::f64::consts::PI;
    let phase = ((seed.wrapping_mul(2654435761u64)) as f64) / (u64::MAX as f64) * 2.0 * PI;
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let bump = 1.0 + amplitude * ((2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos());
            s.set(i, j, bump);
        }
    }
    s
}

fn save_snapshot(s: &Field2D, path: &Path) -> Option<HeightmapMetadata> {
    save_heightmap(s, path).ok()
}

fn variance(field: &Field2D) -> f64 {
    let n = field.data().len() as f64;
    if n <= 0.0 {
        return 0.0;
    }
    let mean: f64 = field.data().iter().sum::<f64>() / n;
    field.data().iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n
}

fn max_abs_grad_s(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    s: &Field2D,
) -> f64 {
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let mut m = 0.0_f64;
    for j in 0..ny {
        for i in 0..nx {
            let ip = idx_x.next(i);
            let jp = idx_y.next(j);
            let gx = (s.get(ip, j) - s.get(i, j)) * inv_dx;
            let gy = (s.get(i, jp) - s.get(i, j)) * inv_dy;
            let mag = (gx * gx + gy * gy).sqrt();
            if mag > m {
                m = mag;
            }
        }
    }
    m
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

    let mut s = init_thickness(nx, ny, cfg.seed, cfg.s_perturbation_amplitude);
    let mut s_next = Field2D::new(nx, ny);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mass_initial = integrated_mass(&s);

    std::fs::create_dir_all(&cfg.output_dir).ok();
    let capture_steps: Vec<usize> = cfg
        .heightmap_fractions
        .iter()
        .map(|f| (f.clamp(0.0, 1.0) * cfg.steps as f64).round() as usize)
        .collect();
    let mut heightmap_paths: Vec<String> = Vec::new();
    let mut heightmap_metas: Vec<HeightmapMetadata> = Vec::new();

    let mut newton_agg = NewtonAggregate::default();
    let mut max_abs_mean_vx = 0.0f64;
    let mut max_abs_mean_vy = 0.0f64;
    let mut vmax_peak = 0.0f64;

    let law_final = ViscosityLaw {
        n: cfg.preset.rheology.n,
        b_prefactor: cfg.preset.rheology.b_prefactor,
        strain_rate_floor: cfg.preset.rheology.strain_rate_floor,
        eta_max_cap: cfg.preset.rheology.eta_max_cap,
        k_saturation: cfg.preset.rheology.k_saturation,
        yielding: cfg.yielding,
    };

    let cg = ConjugateGradient::new(cfg.newton_cfg.linear_tol, cfg.newton_cfg.linear_max_iter);

    // Time series collected at every macro step.
    let mut variance_series: Vec<f64> = Vec::with_capacity(cfg.steps + 1);
    let mut max_grad_s_series: Vec<f64> = Vec::with_capacity(cfg.steps + 1);
    variance_series.push(variance(&s));
    max_grad_s_series.push(max_abs_grad_s(nx, ny, dx, dy, &idx_x, &idx_y, &s));

    let start = Instant::now();

    if capture_steps.contains(&0) {
        let path = cfg.output_dir.join(format!("s_{}x{}_t0000.png", nx, ny));
        if let Some(md) = save_snapshot(&s, &path) {
            heightmap_paths.push(path.display().to_string().replace('\\', "/"));
            heightmap_metas.push(md);
        }
    }

    // --- Accumulate force from the configured BodyForce. ---
    // Helper closure to avoid repeating the zero-and-accumulate pattern.
    let sample_force = |fx: &mut Field2D, fy: &mut Field2D, s: &Field2D| {
        let mut out = VectorField { fx, fy };
        out.zero();
        let state = SimulationState {
            nx, ny, dx, dy,
            idx_x: &idx_x, idx_y: &idx_y,
            s,
        };
        cfg.force.accumulate(&state, &mut out);
    };

    // --- Startup continuation (t = 0 only) ---
    sample_force(&mut fx, &mut fy, &s);
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
    for (_n, oc) in &cont.sub_outcomes {
        record_outcome(oc, &mut newton_agg);
    }

    let sr_after_ramp = StrainRate::compute(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy);
    let eta_after_ramp = rheology::build_eta_field(&law_final, &sr_after_ramp.eps_ii_center);
    newton_agg.cap_fraction_ramp_max = newton_agg
        .cap_fraction_ramp_max
        .max(cap_activation_fraction(&eta_after_ramp, law_final.eta_max_cap));

    // --- Steady-state loop ---
    for step in 0..cfg.steps {
        sample_force(&mut fx, &mut fy, &s);
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
        record_outcome(&outcome, &mut newton_agg);

        let sr = StrainRate::compute(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy);
        let eta_cc = rheology::build_eta_field(&law_final, &sr.eps_ii_center);
        newton_agg.cap_fraction_steady_max = newton_agg
            .cap_fraction_steady_max
            .max(cap_activation_fraction(&eta_cc, law_final.eta_max_cap));
        newton_agg.eta_contrast_samples.push(eta_contrast(&eta_cc));

        // Step 3 — yielding diagnostics. Structural by-pass for
        // Disabled: no extra field build, no blend, counters stay
        // `None`. For Enabled we build a per-timestep `η_visc`
        // companion field to evaluate the "yielding is dominant"
        // and "yielding intensity" aggregates.
        if let crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) = cfg.yielding {
            // Pure-viscous reference at the same ε̇: matches
            // `law_final.eta_effective` with yielding flipped off.
            let mut eta_visc_only = law_final;
            eta_visc_only.yielding =
                crate::tectonics_v2::presets::YieldingConfig::Disabled;
            let eta_visc_cc = rheology::build_eta_field(&eta_visc_only, &sr.eps_ii_center);
            let frac = yielding_cell_fraction(&eta_visc_cc, &eta_cc);
            let intensity = yielding_intensity(&eta_visc_cc, &eta_cc);
            newton_agg.bi_diagnostic = Some(ylaw.bi);
            newton_agg.yielding_cell_fraction_max = Some(
                newton_agg.yielding_cell_fraction_max.unwrap_or(0.0).max(frac),
            );
            newton_agg.yielding_intensity_max = Some(
                newton_agg.yielding_intensity_max.unwrap_or(0.0).max(intensity),
            );

            // Floor-domination diagnostic: domain-level ε̇_II
            // aggregates at the **final** timestep (overwritten
            // each step; only the last one survives).
            let eps_data = sr.eps_ii_center.data();
            let mut sum = 0.0_f64;
            let mut max = 0.0_f64;
            let mut below = 0usize;
            let floor_band = 10.0 * law_final.strain_rate_floor;
            for &v in eps_data {
                sum += v;
                if v > max { max = v; }
                if v < floor_band { below += 1; }
            }
            let n = eps_data.len() as f64;
            newton_agg.eps_ii_mean_final = Some(sum / n);
            newton_agg.eps_ii_max_final = Some(max);
            newton_agg.eps_ii_floor_dominated_fraction_final =
                Some(below as f64 / n);
        }

        let m_vx: f64 = vx.iter().sum::<f64>() / vx.len() as f64;
        let m_vy: f64 = vy.iter().sum::<f64>() / vy.len() as f64;
        max_abs_mean_vx = max_abs_mean_vx.max(m_vx.abs());
        max_abs_mean_vy = max_abs_mean_vy.max(m_vy.abs());

        let vmax_step = vx
            .iter()
            .chain(vy.iter())
            .fold(0.0_f64, |acc, &v| acc.max(v.abs()));
        vmax_peak = vmax_peak.max(vmax_step);

        // Physical step target from total_time / steps. Clamp the
        // CFL-derived dt to that target so the placeholder-induced
        // tiny-v regimes don't walk through the run in a handful of
        // oversized steps (advection discretisation stability aside,
        // the physics cares about the ratio of dt to the dissipation
        // time τ ~ L²·η/(Ar·S²), not about CFL alone).
        let dt_cfl = cfl_dt(dx, dy, &vx, &vy, cfg.cfl_factor);
        let dt_target = if cfg.steps > 0 {
            cfg.total_time_nondim / cfg.steps as f64
        } else {
            cfg.total_time_nondim
        };
        let dt = dt_target.min(dt_cfl).max(0.0);
        if dt.is_finite() && dt > 0.0 {
            step_upwind(nx, ny, dx, dy, dt, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
            std::mem::swap(&mut s, &mut s_next);
        }

        variance_series.push(variance(&s));
        max_grad_s_series.push(max_abs_grad_s(nx, ny, dx, dy, &idx_x, &idx_y, &s));

        let completed = step + 1;
        if capture_steps.contains(&completed) {
            let path = cfg
                .output_dir
                .join(format!("s_{}x{}_t{:04}.png", nx, ny, completed));
            if let Some(md) = save_snapshot(&s, &path) {
                heightmap_paths.push(path.display().to_string().replace('\\', "/"));
                heightmap_metas.push(md);
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
        heightmap_metas,
        variance_series,
        max_grad_s_series,
        newton: Some(newton_agg),
        s_eq: None,
        boundary_type_diversity: None,
        yielding_cell_fraction: None,
        cratonic_stability: None,
        newton_outcome_distribution: None,
        age_field_stats: None,
    };

    let continuation_str = format!("{:?}", cfg.preset.continuation.n_steps);
    let force_name = cfg.force.name();
    let force_detail = match cfg.force_kind {
        ForceKind::Gpe => format!(
            "GpeForce (Ar = {:.3} from scales)",
            Scales::default().argand_number(),
        ),
        ForceKind::Sinusoidal => format!(
            "SinusoidalForce(ε = {}, Lx = {})",
            cfg.sinusoidal_amplitude, cfg.domain_lx,
        ),
    };
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
        body_force: format!("{} [{}]: {}", force_name, cfg.force_kind.label(), force_detail),
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

fn record_outcome(oc: &NonlinearOutcome, na: &mut NewtonAggregate) {
    match oc {
        NonlinearOutcome::Converged { outer_iters, trace, .. } => {
            na.converged += 1;
            na.outer_iters.push(*outer_iters);
            for &c in &trace.linear_iters {
                na.cg_iters_per_newton_step.push(c);
            }
        }
        NonlinearOutcome::Stalled { outer_iters, trace } => {
            na.stalled += 1;
            na.outer_iters.push(*outer_iters);
            for &c in &trace.linear_iters {
                na.cg_iters_per_newton_step.push(c);
            }
        }
        NonlinearOutcome::Diverged { outer_iters, trace, .. } => {
            na.diverged += 1;
            na.outer_iters.push(*outer_iters);
            for &c in &trace.linear_iters {
                na.cg_iters_per_newton_step.push(c);
            }
        }
        NonlinearOutcome::CappedIters { max_iters_hit, trace, .. } => {
            na.capped += 1;
            na.outer_iters.push(*max_iters_hit);
            for &c in &trace.linear_iters {
                na.cg_iters_per_newton_step.push(c);
            }
        }
    }
}
