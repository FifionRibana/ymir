//! Baseline runner.
//!
//! A run is parametrised by a [`Preset`], a [`NonlinearChoice`] and a
//! **body force** supplied as a `Box<dyn BodyForce>`. Step 2 handles
//! both the GPE physics run and the Sinusoidal regression run through
//! the same harness — the only thing that differs between them is
//! which force is passed in.

use std::path::{Path, PathBuf};
use std::time::Instant;

use super::heightmap::{HeightmapMetadata, save_heightmap};
use super::metrics::{IterationHistogram, Metrics, SolverConfigDump, condition_number_estimate};
use super::newton_metrics::{
    NewtonAggregate, cap_activation_fraction, eta_contrast, yielding_cell_fraction,
    yielding_intensity,
};
use crate::tectonics_v2::advection::{cfl_dt, integrated_mass, step_upwind};
use crate::tectonics_v2::basal_drag::{BasalDragConfig, build_drag_diagonal_field};
use crate::tectonics_v2::boundaries::{
    BoundaryConfig, BoundaryMechanismActive, apply_clamp_with_tracking, compute_source_sink_terms,
    div_v_cell, interface_mask, s_continental_collision_mean, s_continental_interior, s_oceanic,
};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::forcing::{MantleForce, SlabPullForce};
use crate::tectonics_v2::forcing::{
    BodyForce, ForceSum, GpeForce, SimulationState, SinusoidalForce, VectorField,
};
use crate::tectonics_v2::mantle::{
    build_mantle_diagonal_field, build_mantle_pattern, generate_stream_function, MantleConfig,
    MantlePattern, StreamFunctionConfig,
};
use crate::tectonics_v2::presets::{Preset, YieldingConfig};
use crate::tectonics_v2::rheology::{self, StrainRate, ViscosityLaw};
use crate::tectonics_v2::scales::Scales;
use crate::tectonics_v2::slab::{
    AccumulationConfig, ConvergenceDirectionConfig, SlabPullConfig, SlabState,
    compute_convergence_direction, compute_q_sub_conv,
};
use crate::tectonics_v2::stokes::Grid;
use crate::tectonics_v2::stokes::continuation::run_continuation;
use crate::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearOutcome, NonlinearSolver,
};
use crate::tectonics_v2::stokes::picard::{PicardConfig, PicardSolver};
use crate::tectonics_v2::stokes::solver::ConjugateGradient;

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
pub fn build_force(
    kind: ForceKind,
    scales: &Scales,
    sin_amplitude: f64,
    domain_lx: f64,
) -> Box<dyn BodyForce> {
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
    /// Basal-drag configuration (Step 4). Default is
    /// [`BasalDragConfig::Disabled`] so Step 0/1/2/3 harness callers
    /// remain structurally unchanged. Step 4 physics passes
    /// `BasalDragConfig::Enabled(BasalDragLaw { br, .. })`; the Br
    /// sweep varies the `br` field.
    pub basal_drag: BasalDragConfig,
    /// Boundary source/sink configuration (Step 5). Default is
    /// [`BoundaryConfig::Disabled`] — no Q evaluation, no clamp, no
    /// tracking. `Enabled` carries the plate-type field, the
    /// boundary-flag field, and the rate coefficients; the time loop
    /// then runs `Advect(S, v) + dt·Q(S, v) → clamp → S̃_next`.
    pub boundary: BoundaryConfig,
    /// Human-readable layout name for the report's config dump.
    /// Set when the boundary config is built from a named
    /// [`crate::tectonics_v2::boundaries::BoundaryLayout`];
    /// otherwise empty.
    pub boundary_layout_name: String,
    /// Slab-pull configuration (Step 7). Default is
    /// [`SlabPullConfig::Disabled`]: the full slab pipeline
    /// (`Q_sub_conv`, ODE, `n̂`, `SlabPullForce`, `m` advection)
    /// is structurally bypassed and Step 0-6 harness callers stay
    /// unchanged. The Step 7 regression verifies this invariant.
    pub slab_pull: SlabPullConfig,
    /// Mantle forcing configuration (Step 8). Default is
    /// [`MantleConfig::Disabled`]: the full mantle pipeline
    /// (stream-function pattern, `MantleForce` RHS term, and the
    /// `coupling · S̃` diagonal augmentation) is structurally
    /// bypassed and Step 0-7 harness callers stay unchanged. The
    /// Step 8 regression verifies this invariant in scalar
    /// parity with Step 7 physics (exact bit-identity by
    /// construction — no mantle contribution means the operator
    /// reproduces Step 7 exactly).
    pub mantle: MantleConfig,
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
            basal_drag: BasalDragConfig::Disabled,
            boundary: BoundaryConfig::Disabled,
            boundary_layout_name: String::new(),
            slab_pull: SlabPullConfig::Disabled,
            mantle: MantleConfig::Disabled,
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

/// Step 5 — initialize `S̃` with per-plate-type means (oceanic = 0.2,
/// continental = 1.0), each modulated by the same deterministic
/// sinusoidal perturbation pattern as [`init_thickness`]. Called
/// only when the run has a prescribed plate-type field via
/// `BoundaryConfig::Enabled`; Step 0-4 harness callers continue to
/// use the plate-type-agnostic `init_thickness`.
///
/// Physical justification: oceanic crust is ~7 km thick (S* scales
/// to 0.2 dimensionless), continental ~35 km (1.0). Starting every
/// cell at 1.0 forces the calibration loop to drain oceanic cells
/// from 1.0 to 0.2 within the 3·τ* simulation window — which it
/// cannot do at Step 5 baseline velocities (`peak|v| ≈ 1e-5`,
/// `|Δv_conv| ≈ 1e-5`, so `Q_sub ≈ 5e-6` per step, drain rate far
/// below the evolution window). Type-aware init sets the initial
/// state to the physical reference thickness, and the calibration
/// then adjusts `k_spread` around this.
fn init_thickness_plate_aware(
    nx: usize,
    ny: usize,
    seed: u64,
    amplitude: f64,
    plate_types: &crate::tectonics_v2::boundaries::PlateTypeField,
) -> Field2D {
    use crate::tectonics_v2::boundaries::PlateType;
    use std::f64::consts::PI;
    let phase = ((seed.wrapping_mul(2654435761u64)) as f64) / (u64::MAX as f64) * 2.0 * PI;
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let bump_scale = (2.0 * PI * x + phase).sin() * (2.0 * PI * y).cos();
            let mean = match plate_types.get(i, j) {
                PlateType::Oceanic => 0.2,
                PlateType::Continental => 1.0,
            };
            s.set(i, j, mean + amplitude * mean * bump_scale);
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

/// Return `(max|∇S̃| over interface cells, max|∇S̃| over the whole
/// domain)`. Both use the same forward-difference stencil as
/// [`max_abs_grad_s`]. Interface cells are those flagged by the
/// companion `interface_mask` field — oceanic cells adjacent to at
/// least one continental, or vice versa.
fn grad_s_interface_vs_global(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    s: &Field2D,
    interface: &[bool],
) -> (f64, f64) {
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let mut max_interface = 0.0_f64;
    let mut max_global = 0.0_f64;
    for j in 0..ny {
        let jp = idx_y.next(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let gx = (s.get(ip, j) - s.get(i, j)) * inv_dx;
            let gy = (s.get(i, jp) - s.get(i, j)) * inv_dy;
            let mag = (gx * gx + gy * gy).sqrt();
            if mag > max_global {
                max_global = mag;
            }
            if interface[j * nx + i] && mag > max_interface {
                max_interface = mag;
            }
        }
    }
    (max_interface, max_global)
}

/// Return `(peak|f_GPE| over interface cells, peak|f_GPE| over the
/// whole domain)` at the final timestep. Uses the same staggered
/// discretisation as `GpeForce::accumulate`. The interface peak is
/// the quantity monitored against issue #78 — a step-change jump
/// between steps signals the spike; a progressive rise is the
/// natural increase in S̃ heterogeneity as later mechanisms land.
fn f_gpe_interface_vs_global(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    s: &Field2D,
    interface: &[bool],
    ar: f64,
) -> (f64, f64) {
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let mut peak_interface = 0.0_f64;
    let mut peak_global = 0.0_f64;
    for j in 0..ny {
        for i in 0..nx {
            let im = idx_x.prev(i);
            let jm = idx_y.prev(j);
            let s_right = s.get(i, j);
            let s_left = s.get(im, j);
            let s_top = s.get(i, j);
            let s_bot = s.get(i, jm);
            let fx = -ar * 0.5 * (s_right + s_left) * (s_right - s_left) * inv_dx;
            let fy = -ar * 0.5 * (s_top + s_bot) * (s_top - s_bot) * inv_dy;
            let mag = (fx * fx + fy * fy).sqrt();
            if mag > peak_global {
                peak_global = mag;
            }
            // A face cell is "on an interface" if either side is
            // interface-classed — use the cell itself as the proxy
            // (its left/bottom faces inherit the classification).
            if interface[j * nx + i] && mag > peak_interface {
                peak_interface = mag;
            }
        }
    }
    (peak_interface, peak_global)
}

fn solve_nonlinear(
    grid: &Grid,
    law: &ViscosityLaw,
    drag_diag: Option<&Field2D>,
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
            solver.solve(grid, law, drag_diag, rhs_x, rhs_y, vx, vy, cg)
        }
        NonlinearChoice::Picard => {
            let solver = PicardSolver::new(picard_cfg);
            solver.solve(grid, law, drag_diag, rhs_x, rhs_y, vx, vy, cg)
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

    // Step 5: when boundary is Enabled, initialize S̃ per-plate-type
    // so oceanic cells start near their physical reference (0.2) and
    // continental cells start near 1.0. Without this, the calibration
    // loop cannot drain oceanic cells from 1.0 down to 0.2 within
    // the 3·τ* window (peak|v| ≈ 1e-5 at Step 5 baseline makes
    // Q_sub ≈ 5e-6 per step). Steps 0-4 callers and the Step 5
    // reference variant / regression (both `BoundaryConfig::Disabled`)
    // keep the plate-type-agnostic init.
    let mut s = match &cfg.boundary {
        BoundaryConfig::Enabled { geometry, .. } => init_thickness_plate_aware(
            nx,
            ny,
            cfg.seed,
            cfg.s_perturbation_amplitude,
            &geometry.plate_type,
        ),
        BoundaryConfig::Disabled => init_thickness(nx, ny, cfg.seed, cfg.s_perturbation_amplitude),
    };
    let mut s_next = Field2D::new(nx, ny);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mass_initial = integrated_mass(&s);

    // Step 5/6 — scratch fields + accumulators for the boundary
    // source/sink pipeline. Only allocated when `BoundaryConfig::Enabled`;
    // the `None` arm drops to a structural by-pass and the per-step
    // cost collapses to exactly Step 4's (modulo the scalar match).
    //
    // Step 6 additions:
    //   - `current_flag`: run-local mutable copy of the boundary-flag
    //     field. For static geometries (Step 5), equals `geometry.boundary_flag`
    //     forever. For Voronoi geometries (Step 6), updated by
    //     `detect_boundaries` each step before the source/sink loop.
    let boundary_enabled = matches!(cfg.boundary, BoundaryConfig::Enabled { .. });
    let mut current_flag: Option<crate::tectonics_v2::boundaries::BoundaryFlagField> =
        match &cfg.boundary {
            BoundaryConfig::Enabled { geometry, .. } => Some(geometry.boundary_flag.clone()),
            BoundaryConfig::Disabled => None,
        };
    let mut div_v_scratch = if boundary_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut q_field = if boundary_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut q_sub_scratch = if boundary_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let cell_area = dx * dy;
    let mut q_integral = 0.0_f64;
    let mut clamp_flux_integral = 0.0_f64;
    let mut clamp_activation_sum = 0.0_f64;
    let mut clamp_activation_max = 0.0_f64;
    let mut clamp_samples: usize = 0;

    // Step 6 Closed-mode state (only populated when the recycling
    // mode is Closed; Open leaves these at their defaults and the
    // Open arm in the time loop doesn't read them).
    let mut closed_state_accumulators =
        crate::tectonics_v2::recycling::ImmediateAccumulators::default();
    let mut closed_state_buffer: Option<crate::tectonics_v2::recycling::DelayedRecycler> = None;
    let mut closed_state_mantle_loss_integral = 0.0_f64;
    let mut closed_state_m_sub_total = 0.0_f64;
    let mut closed_state_arc_distributed = 0.0_f64;
    let mut closed_state_coll_v_distributed = 0.0_f64;
    let mut closed_state_rift_v_distributed = 0.0_f64;
    let mut closed_state_spread_distributed = 0.0_f64;
    // Flag counts captured at step 1 (after the first
    // detect_boundaries call) — informative for the Step 6
    // clarification "did detection fire on step 1?".
    let mut boundary_flag_counts_step1: Option<(usize, usize, usize, usize, usize)> = None;
    // Pre-compute the recycling config and build the buffer if
    // Closed mode. Validation panics at this point are programmer
    // errors (the builder `enabled_voronoi_closed` already validates
    // at construction; direct `BoundaryConfig::Enabled { .. }`
    // construction bypasses the builder and is considered unsafe).
    if let BoundaryConfig::Enabled { recycling_mode, .. } = &cfg.boundary {
        if let crate::tectonics_v2::boundaries::RecyclingModeInit::Closed(rcfg) = recycling_mode {
            rcfg.validate().expect(
                "RecyclingConfig invalid — use BoundaryConfig::enabled_voronoi_closed \
                 or validate() before run_baseline",
            );
            closed_state_buffer =
                Some(crate::tectonics_v2::recycling::DelayedRecycler::new(rcfg.mantle_delay_steps));
        }
    }
    // Tracking for #78 trajectory and other Step 6 samplers.
    let mut trajectory_samples: Vec<(usize, f64, f64, f64, f64, f64)> = Vec::new();
    let trajectory_steps: [usize; 5] = [1, 10, 50, 150, 300];
    let mut buffer_fill_samples: Vec<(usize, f64)> = Vec::new();
    let mut boundary_flag_transition_rate_series: Vec<f64> = Vec::new();
    let mut prev_boundary_flag: Option<crate::tectonics_v2::boundaries::BoundaryFlagField> = None;
    let mut clamp_activation_during_spinup_max: f64 = 0.0;
    let mut recycling_buffer_fill_sum: f64 = 0.0;
    let mut recycling_buffer_fill_max: f64 = 0.0;
    let mut recycling_buffer_fill_samples: usize = 0;
    let mut immediate_pending_max_observed: f64 = 0.0;
    // Step 6 perf: pre-compute per-run flags so the inner loop
    // avoids repeated `match &cfg.boundary` + `matches!` dispatches.
    // `is_dynamic` decides whether `detect_boundaries` + flag-clone
    // run each step; `is_closed` decides whether the recycling
    // pipeline tracks its buffer/pending stats each step;
    // `mantle_delay_for_spinup` short-circuits the spinup branch
    // when in Open mode (unused there).
    let (is_dynamic, is_closed, mantle_delay_for_spinup) = match &cfg.boundary {
        BoundaryConfig::Enabled { geometry, recycling_mode, .. } => {
            let dynamic = geometry.is_dynamic();
            match recycling_mode {
                crate::tectonics_v2::boundaries::RecyclingModeInit::Open => (dynamic, false, 0),
                crate::tectonics_v2::boundaries::RecyclingModeInit::Closed(rcfg) => {
                    (dynamic, true, rcfg.mantle_delay_steps)
                }
            }
        }
        BoundaryConfig::Disabled => (false, false, 0),
    };

    // Step 7 perf: pre-compute `is_slab_enabled` so the per-step
    // branches skip the match+dispatch when slab-pull is Disabled.
    // The scratch buffers are allocated once; `slab_state` is
    // `Some` only when enabled (zero allocation otherwise).
    let is_slab_enabled = matches!(cfg.slab_pull, SlabPullConfig::Enabled { .. });
    let mut slab_state: Option<SlabState> =
        if is_slab_enabled { Some(SlabState::new_zero(nx, ny)) } else { None };
    let mut slab_div_scratch: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_conv_scratch: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_q_sub_conv: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_n_x: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_n_y: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_f_x: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_f_y: Option<Field2D> =
        if is_slab_enabled { Some(Field2D::new(nx, ny)) } else { None };
    // Step 7 metrics accumulators (only meaningful when enabled).
    let mut slab_m_mean_series: Vec<f64> =
        Vec::with_capacity(if is_slab_enabled { cfg.steps } else { 0 });
    let mut slab_m_max_series: Vec<f64> =
        Vec::with_capacity(if is_slab_enabled { cfg.steps } else { 0 });
    let mut peak_f_slab_run: f64 = 0.0;
    let mut peak_f_gpe_run: f64 = 0.0;
    let mut f_slab_to_f_gpe_ratio_sum: f64 = 0.0;
    let mut f_slab_to_f_gpe_ratio_samples: usize = 0;

    // Step 8 — mantle forcing state, pre-allocated once.
    // The pattern is static at Step 8 (evolution_rate = 0, D6);
    // generate it once at run start. The diagonal field is
    // rebuilt each step because S̃ evolves with advection.
    let is_mantle_enabled = matches!(cfg.mantle, MantleConfig::Enabled { .. });
    let mantle_pattern: Option<MantlePattern> = if is_mantle_enabled {
        if let MantleConfig::Enabled { num_modes, seed, .. } = cfg.mantle {
            let psi = generate_stream_function(
                nx, ny,
                &StreamFunctionConfig { num_modes, seed },
            );
            Some(build_mantle_pattern(&psi, dx, dy, &idx_x, &idx_y))
        } else {
            None
        }
    } else {
        None
    };
    let mut slab_f_buf_mantle_x: Option<Field2D> =
        if is_mantle_enabled { Some(Field2D::new(nx, ny)) } else { None };
    let mut slab_f_buf_mantle_y: Option<Field2D> =
        if is_mantle_enabled { Some(Field2D::new(nx, ny)) } else { None };
    // Step 8 metric accumulators.
    let peak_v_mantle_pattern: f64 = mantle_pattern
        .as_ref()
        .map(|p| cfg.mantle.mf_or_zero() * p.peak_magnitude())
        .unwrap_or(0.0);
    let mut peak_v_solved_run: f64 = 0.0;
    let mut peak_f_mantle_run: f64 = 0.0;
    let mut alignment_sum: f64 = 0.0;
    let mut alignment_samples: usize = 0;
    let mut f_mantle_to_f_gpe_ratio_sum: f64 = 0.0;
    let mut f_mantle_to_f_slab_ratio_sum: f64 = 0.0;
    let mut f_mantle_ratio_samples: usize = 0;
    let mut eps_ii_max_to_floor_ratio_run: f64 = 0.0;
    let mut div_v_mantle_max_run: f64 = 0.0;
    // Sanity: check pattern div once at init — should be
    // essentially zero by construction. We track the max over
    // the run as well, though it is constant.
    if let Some(p) = mantle_pattern.as_ref() {
        div_v_mantle_max_run = crate::tectonics_v2::mantle::pattern::pattern_div_max(
            p, dx, dy, &idx_x, &idx_y,
        );
    }

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
        let state = SimulationState { nx, ny, dx, dy, idx_x: &idx_x, idx_y: &idx_y, s };
        cfg.force.accumulate(&state, &mut out);
    };

    // --- Startup continuation (t = 0 only) ---
    sample_force(&mut fx, &mut fy, &s);
    // Build the basal-drag diagonal field from the initial S̃. It is
    // constant throughout the continuation ramp (S̃ does not evolve at
    // t = 0).
    let drag_diag_init = build_drag_diagonal_field(&cfg.basal_drag, &s);
    let newton_solver = NewtonSolver::new(cfg.newton_cfg);
    let cont = run_continuation(
        &grid,
        &law_final,
        drag_diag_init.as_ref(),
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
    // Accumulators for the two Step-4 per-cell diagonal ratios. We
    // average per-step means (not raw per-cell ratios across steps),
    // which matches how the physics report describes the quantities
    // — "mean over cells, averaged across the run".
    let mut basal_drag_energy_ratio_sum = 0.0_f64;
    let mut drag_vs_visc_diagonal_ratio_sum = 0.0_f64;
    let mut drag_ratio_sample_count: usize = 0;

    for step in 0..cfg.steps {
        sample_force(&mut fx, &mut fy, &s);
        // Step 7 — slab pipeline (pre-solve).
        //
        // Ordering (prompt §"Ordering de la boucle temporelle"):
        //   2. Q_sub_conv from current v.
        //   3. ODE: m(t+Δt) = m(t) + Δt·[Q_sub_conv − m/τ].
        //   4. n̂_convergence from current v.
        //   5. Assemble SlabPullForce(m, n̂) on top of GpeForce's
        //      contribution already in (fx, fy).
        //
        // The `current v` at this point is the warm-start from the
        // previous iter's solve (or the continuation result on
        // step 0). `Q_sub_conv` needs `plate_type` + `boundary_flag`
        // from the active `CrustGeometry`; slab-pull without a
        // `BoundaryConfig::Enabled` is not meaningful so we gate on
        // that too.
        //
        // `dt_target = total_time / steps` is used for the ODE
        // step — it is the intended step length and matches the
        // actual `dt` whenever CFL is not binding (the expected
        // regime at Step 7 target parameters). The defensive cap
        // `min(Δt_CFL, 0.1·τ_slab)` mentioned in the spec is
        // enforced here by the `dt_target.min(0.1·τ_slab)`
        // clamp — it never binds with the Step 7 baseline
        // parameters (dt_target ≈ 0.02 < 0.05 = 0.1·0.5) but
        // guards against parameter changes in a sweep.
        let dt_target = if cfg.steps > 0 {
            cfg.total_time_nondim / cfg.steps as f64
        } else {
            cfg.total_time_nondim
        };
        if let (
            SlabPullConfig::Enabled { sp, tau_slab, k_slab_accum, epsilon },
            BoundaryConfig::Enabled { geometry, .. },
        ) = (&cfg.slab_pull, &cfg.boundary)
        {
            let slab = slab_state.as_mut().expect("enabled → slab_state");
            let div_scratch = slab_div_scratch.as_mut().expect("enabled → div_scratch");
            let conv_scratch = slab_conv_scratch.as_mut().expect("enabled → conv_scratch");
            let q_sub_conv_buf = slab_q_sub_conv.as_mut().expect("enabled → q_sub_conv");
            let n_x_buf = slab_n_x.as_mut().expect("enabled → n_x");
            let n_y_buf = slab_n_y.as_mut().expect("enabled → n_y");
            let f_slab_x_buf = slab_f_x.as_mut().expect("enabled → f_slab_x");
            let f_slab_y_buf = slab_f_y.as_mut().expect("enabled → f_slab_y");

            let plate_type = &geometry.plate_type;
            let boundary_flag =
                current_flag.as_ref().expect("BoundaryConfig::Enabled → current_flag");

            // 2. Q_sub_conv = k · max(0, -div v) on oceanic subducting cells.
            compute_q_sub_conv(
                nx,
                ny,
                dx,
                dy,
                &idx_x,
                &idx_y,
                &vx,
                &vy,
                plate_type,
                boundary_flag,
                div_scratch,
                conv_scratch,
                q_sub_conv_buf,
                &AccumulationConfig { k_slab_accum: *k_slab_accum },
            );

            // 3. ODE forward Euler. Stability-capped Δt.
            let dt_ode = dt_target.min(0.1 * *tau_slab);
            slab.step_ode(q_sub_conv_buf, dt_ode, *tau_slab);

            // Metric: series of m-field stats.
            let m_data = slab.m().data();
            let m_mean: f64 = m_data.iter().sum::<f64>() / m_data.len() as f64;
            let m_max: f64 = m_data.iter().cloned().fold(0.0_f64, f64::max);
            slab_m_mean_series.push(m_mean);
            slab_m_max_series.push(m_max);

            // 4. n̂_convergence from current v.
            compute_convergence_direction(
                nx,
                ny,
                dx,
                dy,
                &idx_x,
                &idx_y,
                &vx,
                &vy,
                div_scratch,
                n_x_buf,
                n_y_buf,
                &ConvergenceDirectionConfig { epsilon: *epsilon },
            );

            // 5. Assemble SlabPullForce into a dedicated buffer so
            //    we can measure `peak|f_slab|` and the ratio to
            //    `peak|f_GPE|` without aliasing. Then add into the
            //    main (fx, fy) RHS carrying the GPE contribution.
            let peak_f_gpe_step = {
                let mut peak = 0.0_f64;
                for (a, b) in fx.data().iter().zip(fy.data().iter()) {
                    peak = peak.max((a * a + b * b).sqrt());
                }
                peak
            };
            for v in f_slab_x_buf.data_mut().iter_mut() {
                *v = 0.0;
            }
            for v in f_slab_y_buf.data_mut().iter_mut() {
                *v = 0.0;
            }
            let slab_state_sim =
                SimulationState { nx, ny, dx, dy, idx_x: &idx_x, idx_y: &idx_y, s: &s };
            SlabPullForce::new(*sp, slab.m(), n_x_buf, n_y_buf).accumulate(
                &slab_state_sim,
                &mut VectorField { fx: f_slab_x_buf, fy: f_slab_y_buf },
            );
            let peak_f_slab_step = {
                let mut peak = 0.0_f64;
                for (a, b) in f_slab_x_buf.data().iter().zip(f_slab_y_buf.data().iter()) {
                    peak = peak.max((a * a + b * b).sqrt());
                }
                peak
            };
            peak_f_slab_run = peak_f_slab_run.max(peak_f_slab_step);
            peak_f_gpe_run = peak_f_gpe_run.max(peak_f_gpe_step);
            if peak_f_gpe_step > 0.0 {
                f_slab_to_f_gpe_ratio_sum += peak_f_slab_step / peak_f_gpe_step;
                f_slab_to_f_gpe_ratio_samples += 1;
            }
            // Add f_slab into the main RHS. `BodyForce` contract is
            // additive, but the harness chose to run them into
            // different buffers for diagnostic isolation.
            let fx_slice = fx.data_mut();
            let fy_slice = fy.data_mut();
            let fsx = f_slab_x_buf.data();
            let fsy = f_slab_y_buf.data();
            for k in 0..nx * ny {
                fx_slice[k] += fsx[k];
                fy_slice[k] += fsy[k];
            }
        }
        // Step 8 — mantle forcing pre-solve contribution.
        //
        // Formulation: f_mantle = coupling · S̃ · (Mf · v_pattern − v_solved).
        // Split into:
        //   • constant RHS part (handled by `MantleForce` body force,
        //     added into (fx, fy) alongside GPE + slab contributions)
        //   • `coupling · S̃` diagonal augmentation (summed with
        //     drag_diag into `total_diag`), which produces exact
        //     self-consistency at every Newton outer iteration by
        //     linearly coupling `v_solved` back into the operator.
        //
        // This mirrors the Step 4 basal-drag approach: the inner CG
        // sees only `A(v;η) + total_diag · I` — an SPD augmentation
        // — with no mantle-specific dispatch. `v_solved` converges
        // to the self-consistent balance via Newton's own outer loop.
        if is_mantle_enabled {
            if let MantleConfig::Enabled { mf, coupling, .. } = cfg.mantle {
                let pattern = mantle_pattern.as_ref().expect("enabled → pattern");
                let f_mantle_x = slab_f_buf_mantle_x.as_mut().expect("enabled → buf");
                let f_mantle_y = slab_f_buf_mantle_y.as_mut().expect("enabled → buf");
                for v in f_mantle_x.data_mut().iter_mut() { *v = 0.0; }
                for v in f_mantle_y.data_mut().iter_mut() { *v = 0.0; }
                let mantle_state =
                    SimulationState { nx, ny, dx, dy, idx_x: &idx_x, idx_y: &idx_y, s: &s };
                MantleForce::new(mf, coupling, pattern, &s).accumulate(
                    &mantle_state,
                    &mut VectorField { fx: f_mantle_x, fy: f_mantle_y },
                );
                let peak_f_mantle_step = {
                    let mut peak = 0.0_f64;
                    for (a, b) in f_mantle_x.data().iter().zip(f_mantle_y.data().iter()) {
                        let mag = (a * a + b * b).sqrt();
                        if mag > peak { peak = mag; }
                    }
                    peak
                };
                peak_f_mantle_run = peak_f_mantle_run.max(peak_f_mantle_step);
                // Ratios vs GPE / slab are captured here using the
                // peaks the slab branch already computed above. If
                // the slab branch did not run (slab disabled), we
                // fall back to the plain `peak|f_gpe|` on (fx, fy)
                // pre-mantle (these carry GPE only).
                let peak_f_gpe_fallback = {
                    let mut peak = 0.0_f64;
                    for (a, b) in fx.data().iter().zip(fy.data().iter()) {
                        let mag = (a * a + b * b).sqrt();
                        if mag > peak { peak = mag; }
                    }
                    peak
                };
                if peak_f_gpe_fallback > 0.0 {
                    f_mantle_to_f_gpe_ratio_sum += peak_f_mantle_step / peak_f_gpe_fallback;
                    peak_f_gpe_run = peak_f_gpe_run.max(peak_f_gpe_fallback);
                }
                // When slab was enabled, peak_f_slab_run now holds
                // the slab peak from this step's accumulation. If
                // both slab and mantle fire in the same step, we
                // compute the mantle/slab ratio from those peaks.
                if is_slab_enabled && peak_f_slab_run > 0.0 {
                    // Use the most recent slab sample — the loop
                    // already updated `peak_f_slab_run` above; we
                    // approximate the current slab peak as that
                    // running max (close enough for a telemetry
                    // ratio, and conservative in the sense of
                    // giving the smallest ratio).
                    f_mantle_to_f_slab_ratio_sum += peak_f_mantle_step / peak_f_slab_run;
                }
                f_mantle_ratio_samples += 1;
                // Add f_mantle into the main RHS.
                let fx_slice = fx.data_mut();
                let fy_slice = fy.data_mut();
                let fmx = f_mantle_x.data();
                let fmy = f_mantle_y.data();
                for k in 0..nx * ny {
                    fx_slice[k] += fmx[k];
                    fy_slice[k] += fmy[k];
                }
            }
        }
        // Rebuild the drag-diagonal field from the freshly-advected
        // S̃ (or from the initial S̃ on step 0). Disabled → None, no
        // allocation. Step 8: sum with mantle-diagonal contribution
        // so the Stokes operator sees `A + (drag + mantle) · I`.
        let drag_diag_step = build_drag_diagonal_field(&cfg.basal_drag, &s);
        let mantle_diag_step = build_mantle_diagonal_field(&cfg.mantle, &s);
        let total_diag_step: Option<Field2D> = match (drag_diag_step.as_ref(), mantle_diag_step.as_ref()) {
            (None, None) => None,
            (Some(d), None) => Some(d.clone()),
            (None, Some(m)) => Some(m.clone()),
            (Some(d), Some(m)) => {
                let mut t = Field2D::new(d.nx(), d.ny());
                let dd = d.data();
                let md = m.data();
                let td = t.data_mut();
                for k in 0..dd.len() {
                    td[k] = dd[k] + md[k];
                }
                Some(t)
            }
        };
        let outcome = solve_nonlinear(
            &grid,
            &law_final,
            total_diag_step.as_ref(),
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

        // Step 8 metric: peak|v_solved| per step, running max.
        if is_mantle_enabled {
            let peak_v_step = vx.iter().zip(vy.iter())
                .fold(0.0_f64, |acc, (a, b)| {
                    let m = (a * a + b * b).sqrt();
                    if m > acc { m } else { acc }
                });
            peak_v_solved_run = peak_v_solved_run.max(peak_v_step);
            // Alignment: <v_solved, Mf·v_pattern> / |Mf·v_pattern|².
            // Direction-aware magnitude alignment, matches the
            // contract probed in v2_mantle_relaxation.
            if let (Some(pat), MantleConfig::Enabled { mf, .. }) = (mantle_pattern.as_ref(), cfg.mantle) {
                let mut dot = 0.0_f64;
                let mut norm_m_sq = 0.0_f64;
                for k in 0..nx * ny {
                    let vmx = mf * pat.v_mantle_x.data()[k];
                    let vmy = mf * pat.v_mantle_y.data()[k];
                    dot += vx[k] * vmx + vy[k] * vmy;
                    norm_m_sq += vmx * vmx + vmy * vmy;
                }
                if norm_m_sq > 0.0 {
                    alignment_sum += dot / norm_m_sq;
                    alignment_samples += 1;
                }
            }
        }

        let sr = StrainRate::compute(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy);
        let eta_cc = rheology::build_eta_field(&law_final, &sr.eps_ii_center);
        newton_agg.cap_fraction_steady_max = newton_agg
            .cap_fraction_steady_max
            .max(cap_activation_fraction(&eta_cc, law_final.eta_max_cap));
        newton_agg.eta_contrast_samples.push(eta_contrast(&eta_cc));

        // Step 8 metric — the primary diagnostic that the system
        // has escaped floor-domination. `max(ε̇_II) / ε̇_min ≥ 1`
        // is the analytic condition for yielding dominance to be
        // achievable at Bi = 0.15 (the floor-dominated analysis
        // in the Step 3 report). We track the max over the run.
        if is_mantle_enabled {
            let eps_max = sr.eps_ii_center.data().iter().cloned().fold(0.0_f64, f64::max);
            let ratio = eps_max / law_final.strain_rate_floor.max(1.0e-300);
            if ratio > eps_ii_max_to_floor_ratio_run {
                eps_ii_max_to_floor_ratio_run = ratio;
            }
        }

        // Step 4 diagnostics — basal-drag ratios vs the viscous
        // diagonal η/Δx². Computed at cell centres (before
        // face-averaging) which is enough for an order-of-magnitude
        // diagnostic; the `build_drag_diagonal_field` output already
        // is Br · S̃^exp at cell centres.
        if let (BasalDragConfig::Enabled(law), Some(drag)) =
            (&cfg.basal_drag, drag_diag_step.as_ref())
        {
            let inv_dx2 = 1.0 / (dx * dx);
            let n_cells = eta_cc.data().len().max(1) as f64;
            let mut step_energy_ratio_sum = 0.0_f64;
            let mut step_ratio_sum = 0.0_f64;
            for j in 0..ny {
                for i in 0..nx {
                    let visc = eta_cc.get(i, j) * inv_dx2;
                    let drag_v = drag.get(i, j);
                    let denom = drag_v + visc;
                    let energy_ratio = if denom > 0.0 { drag_v / denom } else { 0.0 };
                    let ratio = if visc > 0.0 { drag_v / visc } else { 0.0 };
                    step_energy_ratio_sum += energy_ratio;
                    step_ratio_sum += ratio;
                }
            }
            basal_drag_energy_ratio_sum += step_energy_ratio_sum / n_cells;
            drag_vs_visc_diagonal_ratio_sum += step_ratio_sum / n_cells;
            drag_ratio_sample_count += 1;
            // Record the Br value once (constant over the run).
            newton_agg.br_diagnostic = Some(law.br);
        }

        // Step 3 — yielding diagnostics. Structural by-pass for
        // Disabled: no extra field build, no blend, counters stay
        // `None`. For Enabled we build a per-timestep `η_visc`
        // companion field to evaluate the "yielding is dominant"
        // and "yielding intensity" aggregates.
        if let crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) = cfg.yielding {
            // Pure-viscous reference at the same ε̇: matches
            // `law_final.eta_effective` with yielding flipped off.
            let mut eta_visc_only = law_final;
            eta_visc_only.yielding = crate::tectonics_v2::presets::YieldingConfig::Disabled;
            let eta_visc_cc = rheology::build_eta_field(&eta_visc_only, &sr.eps_ii_center);
            let frac = yielding_cell_fraction(&eta_visc_cc, &eta_cc);
            let intensity = yielding_intensity(&eta_visc_cc, &eta_cc);
            newton_agg.bi_diagnostic = Some(ylaw.bi);
            newton_agg.yielding_cell_fraction_max =
                Some(newton_agg.yielding_cell_fraction_max.unwrap_or(0.0).max(frac));
            newton_agg.yielding_intensity_max =
                Some(newton_agg.yielding_intensity_max.unwrap_or(0.0).max(intensity));

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
                if v > max {
                    max = v;
                }
                if v < floor_band {
                    below += 1;
                }
            }
            let n = eps_data.len() as f64;
            newton_agg.eps_ii_mean_final = Some(sum / n);
            newton_agg.eps_ii_max_final = Some(max);
            newton_agg.eps_ii_floor_dominated_fraction_final = Some(below as f64 / n);
        }

        let m_vx: f64 = vx.iter().sum::<f64>() / vx.len() as f64;
        let m_vy: f64 = vy.iter().sum::<f64>() / vy.len() as f64;
        max_abs_mean_vx = max_abs_mean_vx.max(m_vx.abs());
        max_abs_mean_vy = max_abs_mean_vy.max(m_vy.abs());

        let vmax_step = vx.iter().chain(vy.iter()).fold(0.0_f64, |acc, &v| acc.max(v.abs()));
        vmax_peak = vmax_peak.max(vmax_step);

        // Physical step target from total_time / steps. Clamp the
        // CFL-derived dt to that target so the placeholder-induced
        // tiny-v regimes don't walk through the run in a handful of
        // oversized steps (advection discretisation stability aside,
        // the physics cares about the ratio of dt to the dissipation
        // time τ ~ L²·η/(Ar·S²), not about CFL alone).
        let dt_cfl = cfl_dt(dx, dy, &vx, &vy, cfg.cfl_factor);
        // `dt_target` is computed at the top of the loop body (Step 7
        // pre-solve needs it for the ODE step). Reuse here.
        let dt = dt_target.min(dt_cfl).max(0.0);
        if dt.is_finite() && dt > 0.0 {
            step_upwind(nx, ny, dx, dy, dt, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
            // Boundary source/sinks — issue #89 D9 (Lie splitting):
            // `Q` is evaluated on `S̃(t)` and `ṽ(t)` (not post-
            // advection state) and accumulated into the advected
            // `s_next` as `s_next += dt·Q`. Clamp afterwards and
            // track the artificial flux for the mass balance.
            //
            // Step 6 refactor: boundary state is carried by
            // `CrustGeometry` (plate_type + plate_id + initial
            // boundary_flag); dynamic geometries update the run-local
            // `current_flag` before each source/sink call. Static
            // geometries (Step 5 layouts) leave `current_flag` at its
            // initial value forever — behaviour identical to Step 5.
            if let BoundaryConfig::Enabled { geometry, rates, recycling_mode } = &cfg.boundary {
                let div_v = div_v_scratch.as_mut().expect("boundary enabled → scratch");
                let q = q_field.as_mut().expect("boundary enabled → q");
                let q_sub = q_sub_scratch.as_mut().expect("boundary enabled → q_sub");
                let cur_flag = current_flag.as_mut().expect("boundary enabled → current_flag");

                // Update dynamic geometries' boundary_flag from the
                // current velocity field.
                if geometry.is_dynamic() {
                    crate::tectonics_v2::boundary_detection::detect_boundaries(
                        nx,
                        ny,
                        dx,
                        dy,
                        &idx_x,
                        &idx_y,
                        &vx,
                        &vy,
                        &geometry.plate_type,
                        &geometry.plate_id,
                        &geometry.detection_config,
                        cur_flag,
                    );
                }

                div_v_cell(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy, div_v);

                match recycling_mode {
                    crate::tectonics_v2::boundaries::RecyclingModeInit::Open => {
                        // Step 5 path: per-cell rate-based source/sinks.
                        // Arithmetic unchanged to preserve bit-parity
                        // with the committed Step 5 physics baseline.
                        compute_source_sink_terms(
                            &geometry.plate_type,
                            cur_flag,
                            rates,
                            div_v,
                            &idx_x,
                            &idx_y,
                            q_sub,
                            q,
                        );
                    }
                    crate::tectonics_v2::boundaries::RecyclingModeInit::Closed(rcfg) => {
                        // Step 6 closed-mode recycling. The Q field
                        // is built as the sum of:
                        //   Q_sub (per-cell rate, unchanged)
                        //   Q_arc + Q_coll_v + Q_rift_v (immediate
                        //     distribution of accumulator budgets
                        //     on eligible cells, rollover on absence)
                        //   Q_spread (delayed buffer output,
                        //     uniformly distributed on oceanic rift)
                        // Mantle loss never reaches the grid — it is
                        // a silent fraction of M_sub subtracted from
                        // the global budget (tracked separately for
                        // the mass conservation check).
                        use crate::tectonics_v2::boundaries::closed_mode::{
                            compute_q_sub_only, distribute_delayed, distribute_immediate,
                            integrate_sub_mass,
                        };
                        // Compute Q_sub only; write into `q` (Q
                        // field) directly. `q_sub` scratch is reused
                        // for clarity / mass accounting.
                        compute_q_sub_only(&geometry.plate_type, cur_flag, rates, div_v, q_sub);
                        // Copy Q_sub into q (Q starts as Q_sub).
                        for (qv, &sv) in q.data_mut().iter_mut().zip(q_sub.data().iter()) {
                            *qv = sv;
                        }
                        let m_sub_step = integrate_sub_mass(q_sub, dt, cell_area);
                        // Update accumulators from this step's budget.
                        closed_state_accumulators.arc_pending += rcfg.arc_fraction * m_sub_step;
                        closed_state_accumulators.coll_v_pending +=
                            rcfg.coll_v_fraction * m_sub_step;
                        closed_state_accumulators.rift_v_pending +=
                            rcfg.rift_v_fraction * m_sub_step;
                        // Distribute immediate (adds to Q). Track
                        // actually-distributed mass per class by
                        // snapshotting the pending values before/
                        // after the call — when distribution fires,
                        // the pending drops to 0 and the difference
                        // is the mass emitted this step. On rollover
                        // steps (no eligible cells) the diff is 0
                        // and pending carries over.
                        let pre_arc = closed_state_accumulators.arc_pending;
                        let pre_coll = closed_state_accumulators.coll_v_pending;
                        let pre_rift = closed_state_accumulators.rift_v_pending;
                        distribute_immediate(
                            &geometry.plate_type,
                            cur_flag,
                            &mut closed_state_accumulators,
                            &idx_x,
                            &idx_y,
                            dt,
                            cell_area,
                            q,
                        );
                        closed_state_arc_distributed +=
                            pre_arc - closed_state_accumulators.arc_pending;
                        closed_state_coll_v_distributed +=
                            pre_coll - closed_state_accumulators.coll_v_pending;
                        closed_state_rift_v_distributed +=
                            pre_rift - closed_state_accumulators.rift_v_pending;
                        // Delayed: deposit this step's spread
                        // fraction, advance-or-rollover, distribute
                        // the emerging mass.
                        let buffer =
                            closed_state_buffer.as_mut().expect("Closed mode → buffer initialised");
                        buffer.deposit(rcfg.spread_fraction * m_sub_step);
                        let emerged = distribute_delayed(
                            &geometry.plate_type,
                            cur_flag,
                            buffer,
                            dt,
                            cell_area,
                            q,
                        );
                        closed_state_spread_distributed += emerged;
                        // Mantle loss accumulator.
                        closed_state_mantle_loss_integral += rcfg.mantle_loss_fraction * m_sub_step;
                        closed_state_m_sub_total += m_sub_step;
                    }
                }
                // Integrate Q*dt into s_next and accumulate the
                // physical flux ∫Q dt dA.
                let mut q_sum_step = 0.0_f64;
                for (s_cell, &q_val) in s_next.data_mut().iter_mut().zip(q.data().iter()) {
                    *s_cell += dt * q_val;
                    q_sum_step += q_val;
                }
                q_integral += dt * q_sum_step * cell_area;
                // Apply the hard floor + track the injected flux.
                let clamp_stats = apply_clamp_with_tracking(&mut s_next);
                clamp_flux_integral += clamp_stats.injected_flux * cell_area;
                let frac = clamp_stats.activation_fraction();
                clamp_activation_sum += frac;
                if frac > clamp_activation_max {
                    clamp_activation_max = frac;
                }
                clamp_samples += 1;

                // Step 6 — per-step trackers. Gated on pre-computed
                // `is_dynamic` / `is_closed` to avoid re-matching
                // `cfg.boundary` each step. For Step 5-shape
                // regression (static + Open) both are false and
                // this block short-circuits to zero extra work
                // beyond the Step 5 path.
                if is_dynamic {
                    if let Some(prev) = prev_boundary_flag.as_ref() {
                        let mut diff = 0usize;
                        let n = cur_flag.data().len();
                        for (a, b) in prev.data().iter().zip(cur_flag.data().iter()) {
                            if a != b {
                                diff += 1;
                            }
                        }
                        boundary_flag_transition_rate_series.push(diff as f64 / n as f64);
                    }
                    prev_boundary_flag = Some(cur_flag.clone());
                    // Capture flag-type counts at the FIRST step
                    // (after detect_boundaries fired once). Used
                    // by the physics report to prove detection is
                    // actually running each step.
                    if step == 0 && boundary_flag_counts_step1.is_none() {
                        use crate::tectonics_v2::boundaries::BoundaryFlag;
                        let (mut n0, mut ns, mut nos, mut nr, mut nc) = (0, 0, 0, 0, 0);
                        for &f in cur_flag.data() {
                            match f {
                                BoundaryFlag::None => n0 += 1,
                                BoundaryFlag::Subduction => ns += 1,
                                BoundaryFlag::OceanicSubduction => nos += 1,
                                BoundaryFlag::Rift => nr += 1,
                                BoundaryFlag::ContinentalCollision => nc += 1,
                            }
                        }
                        boundary_flag_counts_step1 = Some((n0, ns, nos, nr, nc));
                    }
                }
                if is_closed {
                    if step < mantle_delay_for_spinup && frac > clamp_activation_during_spinup_max {
                        clamp_activation_during_spinup_max = frac;
                    }
                    if let Some(ref buf) = closed_state_buffer {
                        let f = buf.fill();
                        recycling_buffer_fill_sum += f;
                        if f > recycling_buffer_fill_max {
                            recycling_buffer_fill_max = f;
                        }
                        recycling_buffer_fill_samples += 1;
                    }
                    let pending_max = closed_state_accumulators.max_pending();
                    if pending_max > immediate_pending_max_observed {
                        immediate_pending_max_observed = pending_max;
                    }
                }
            }
            std::mem::swap(&mut s, &mut s_next);

            // Step 7 — advect `m_subducted` with the solved velocity.
            // Runs AFTER the Stokes solve and AFTER the S advection
            // swap, using the same `v` and `dt` as S (consistency
            // requirement of the Lie splitting; see the prompt's
            // "Advection utilise la velocity post-gauge" note).
            if is_slab_enabled {
                let slab = slab_state.as_mut().expect("enabled → slab_state");
                slab.advect(dx, dy, dt, &idx_x, &idx_y, &vx, &vy);
            }

            // Step 6 — #78 trajectory samples at {1, 10, 50, 150, 300}.
            // Fires at most 5 times per run; the `contains` check
            // alone is cheap but the inner interface_mask build is
            // O(n_cells) — so we gate on both `trajectory_steps` and
            // `boundary_enabled` (cached) to avoid building the mask
            // on Disabled runs.
            let completed = step + 1;
            if boundary_enabled && trajectory_steps.contains(&completed) {
                if let BoundaryConfig::Enabled { geometry, .. } = &cfg.boundary {
                    let inter = interface_mask(&geometry.plate_type, &idx_x, &idx_y);
                    let (grad_i, grad_g) =
                        grad_s_interface_vs_global(nx, ny, dx, dy, &idx_x, &idx_y, &s, &inter);
                    let (fg_i, fg_g) = f_gpe_interface_vs_global(
                        nx,
                        ny,
                        dx,
                        dy,
                        &idx_x,
                        &idx_y,
                        &s,
                        &inter,
                        Scales::default().argand_number(),
                    );
                    let buf_fill = closed_state_buffer.as_ref().map(|b| b.fill()).unwrap_or(0.0);
                    trajectory_samples.push((completed, grad_i, grad_g, fg_i, fg_g, buf_fill));
                    buffer_fill_samples.push((completed, buf_fill));
                }
            }
        }

        variance_series.push(variance(&s));
        max_grad_s_series.push(max_abs_grad_s(nx, ny, dx, dy, &idx_x, &idx_y, &s));

        let completed = step + 1;
        if capture_steps.contains(&completed) {
            let path = cfg.output_dir.join(format!("s_{}x{}_t{:04}.png", nx, ny, completed));
            if let Some(md) = save_snapshot(&s, &path) {
                heightmap_paths.push(path.display().to_string().replace('\\', "/"));
                heightmap_metas.push(md);
            }
        }
    }

    let wallclock = start.elapsed();
    let mass_final = integrated_mass(&s);
    let drift = (mass_final - mass_initial) / mass_initial.abs().max(1.0);

    // Commit the Step-4 drag ratio averages (only populated under
    // `BasalDragConfig::Enabled`).
    if drag_ratio_sample_count > 0 {
        let n = drag_ratio_sample_count as f64;
        newton_agg.basal_drag_energy_ratio = Some(basal_drag_energy_ratio_sum / n);
        newton_agg.drag_vs_visc_diagonal_ratio = Some(drag_vs_visc_diagonal_ratio_sum / n);
    }

    // Step 5 — boundary source/sink aggregates and final-state
    // stats. Entire block is skipped under `BoundaryConfig::Disabled`
    // (the Newton-agg Option slots stay `None`, which the report
    // writer treats as "not applicable").
    if let BoundaryConfig::Enabled { geometry, rates, .. } = &cfg.boundary {
        newton_agg.boundary_layout_name = None; // set by CLI via a separate hook if needed
        let pt: &crate::tectonics_v2::boundaries::PlateTypeField = &geometry.plate_type;
        // Post-run stats use the FINAL boundary_flag state. For
        // static geometries (Step 5) that equals the initial flag;
        // for Voronoi geometries (Step 6) it reflects the last
        // `detect_boundaries` call at step N.
        let fl_final = current_flag.as_ref().expect("boundary enabled → current_flag set at init");
        let fl = fl_final;
        let s_ocean = s_oceanic(&s, pt);
        let s_cont = s_continental_interior(&s, pt, fl);
        newton_agg.s_oceanic_mean = if s_ocean.count > 0 { Some(s_ocean.mean) } else { None };
        newton_agg.s_oceanic_std = if s_ocean.count > 0 { Some(s_ocean.std) } else { None };
        newton_agg.s_continental_interior_mean =
            if s_cont.count > 0 { Some(s_cont.mean) } else { None };
        newton_agg.s_continental_interior_std =
            if s_cont.count > 0 { Some(s_cont.std) } else { None };
        newton_agg.s_continental_collision_mean = s_continental_collision_mean(&s, pt, fl);
        // Diversity computation is mode-aware: Open mode reads
        // rate coefficients, Closed mode reads the recycling config
        // fractions (rates are typically zeroed for non-sub channels
        // in Closed mode since they're replaced by distributive
        // fractions — querying them would hide real activity).
        let mech_active = match &cfg.boundary {
            BoundaryConfig::Enabled { recycling_mode, .. } => match recycling_mode {
                crate::tectonics_v2::boundaries::RecyclingModeInit::Open => {
                    BoundaryMechanismActive::from(rates)
                }
                crate::tectonics_v2::boundaries::RecyclingModeInit::Closed(rcfg) => {
                    BoundaryMechanismActive::from_closed_mode(rates, rcfg)
                }
            },
            BoundaryConfig::Disabled => BoundaryMechanismActive::from(rates),
        };
        newton_agg.boundary_type_diversity =
            Some(crate::tectonics_v2::boundaries::boundary_type_diversity(pt, fl, mech_active));
        if clamp_samples > 0 {
            newton_agg.clamp_activation_fraction_mean =
                Some(clamp_activation_sum / clamp_samples as f64);
            newton_agg.clamp_activation_fraction_max = Some(clamp_activation_max);
        } else {
            newton_agg.clamp_activation_fraction_mean = Some(0.0);
            newton_agg.clamp_activation_fraction_max = Some(0.0);
        }
        newton_agg.q_integral = Some(q_integral);
        newton_agg.clamp_flux_integral = Some(clamp_flux_integral);
        // Mass-balance residual per issue #89 D5:
        //   |Δmass_obs − ∫Q − ∫clamp_flux| / max(|∫Q|+|∫clamp|, 1)
        // The `max(_, 1)` (on absolute mass scale) keeps the
        // denominator away from zero in the balanced layout. The
        // `reference_scale` is the total "expected traffic"
        // through the balance, which matches the spec's intent.
        let delta_mass_obs = (mass_final - mass_initial) * cell_area;
        let numerator = (delta_mass_obs - q_integral - clamp_flux_integral).abs();
        let denom_raw = q_integral.abs() + clamp_flux_integral.abs();
        // Keep the denominator away from zero: floor at the
        // machine-noise scale of `mass_initial * cell_area` so that
        // an exactly balanced layout does not divide by an
        // arbitrarily small "expected flux". The floor only kicks in
        // when the physical flux is genuinely zero.
        let denom_floor = 1.0e-12 * (mass_initial * cell_area).abs().max(1.0);
        let denom = denom_raw.max(denom_floor);
        newton_agg.mass_balance_residual = Some(numerator / denom);

        // Step 5 — #78 monitoring: max|∇S| and peak|f_GPE| on
        // oceanic/continental interface cells vs global reference.
        let inter = interface_mask(pt, &idx_x, &idx_y);
        let (grad_interface, grad_global) =
            grad_s_interface_vs_global(nx, ny, dx, dy, &idx_x, &idx_y, &s, &inter);
        newton_agg.max_grad_s_interface_final = Some(grad_interface);
        newton_agg.max_grad_s_global_final = Some(grad_global);
        let (fgpe_interface, fgpe_global) = f_gpe_interface_vs_global(
            nx,
            ny,
            dx,
            dy,
            &idx_x,
            &idx_y,
            &s,
            &inter,
            Scales::default().argand_number(),
        );
        newton_agg.peak_f_gpe_interface_final = Some(fgpe_interface);
        newton_agg.peak_f_gpe_global_final = Some(fgpe_global);

        // Step 6 — additional aggregates.
        newton_agg.plate_count = Some(geometry.distinct_plate_count() as u32);
        let n_cells = pt.data().len() as f64;
        let n_cont = geometry.continental_cell_count() as f64;
        newton_agg.plate_type_distribution = Some((1.0 - n_cont / n_cells, n_cont / n_cells));
        if !boundary_flag_transition_rate_series.is_empty() {
            let sum: f64 = boundary_flag_transition_rate_series.iter().sum();
            let cnt = boundary_flag_transition_rate_series.len() as f64;
            let max_val =
                boundary_flag_transition_rate_series.iter().copied().fold(0.0_f64, f64::max);
            newton_agg.boundary_flag_transition_rate_mean = Some(sum / cnt);
            newton_agg.boundary_flag_transition_rate_max = Some(max_val);
        }
        // Closed-mode aggregates: populated only when closed_state_buffer was built.
        if let Some(ref buf) = closed_state_buffer {
            if recycling_buffer_fill_samples > 0 {
                newton_agg.recycling_buffer_fill_mean =
                    Some(recycling_buffer_fill_sum / recycling_buffer_fill_samples as f64);
                newton_agg.recycling_buffer_fill_max = Some(recycling_buffer_fill_max);
            }
            newton_agg.recycling_buffer_fill_final = Some(buf.fill());
            newton_agg.immediate_pending_max = Some(immediate_pending_max_observed);
            newton_agg.immediate_pending_final = Some(closed_state_accumulators.sum());
            newton_agg.mantle_loss_integral = Some(closed_state_mantle_loss_integral);
            newton_agg.m_sub_total = Some(closed_state_m_sub_total);
            newton_agg.clamp_activation_during_spinup_max =
                Some(clamp_activation_during_spinup_max);
            // Absolute mass-conservation residual (Step 6):
            //   Δmass_obs + mantle_loss + buffer_fill + pending − clamp_flux = 0
            // We keep the raw integrated values in cell-area units
            // (q_integral and clamp_flux_integral already carry
            // `· cell_area`). mass_obs also in cell-area units.
            let delta_mass_obs = (mass_final - mass_initial) * cell_area;
            let residual_signed = delta_mass_obs
                + closed_state_mantle_loss_integral
                + buf.fill()
                + closed_state_accumulators.sum()
                - clamp_flux_integral;
            let denom = (mass_initial * cell_area).abs().max(1.0);
            newton_agg.mass_conservation_residual = Some(residual_signed.abs() / denom);
        }
        if !trajectory_samples.is_empty() {
            newton_agg.issue_78_trajectory = trajectory_samples.clone();
        }
        let _ = buffer_fill_samples; // superseded by issue_78_trajectory's last column

        // Step 6 final per-flag-type count. Populated whenever
        // boundary is Enabled; disambiguates `boundary_type_diversity`
        // by showing which flag types were actually present at end.
        {
            use crate::tectonics_v2::boundaries::BoundaryFlag;
            let mut none = 0usize;
            let mut sub = 0usize;
            let mut osub = 0usize;
            let mut rift = 0usize;
            let mut coll = 0usize;
            for &f in fl.data() {
                match f {
                    BoundaryFlag::None => none += 1,
                    BoundaryFlag::Subduction => sub += 1,
                    BoundaryFlag::OceanicSubduction => osub += 1,
                    BoundaryFlag::Rift => rift += 1,
                    BoundaryFlag::ContinentalCollision => coll += 1,
                }
            }
            newton_agg.boundary_flag_counts_final = Some((none, sub, osub, rift, coll));
            newton_agg.boundary_flag_counts_step1 = boundary_flag_counts_step1;
        }
        if is_closed {
            newton_agg.arc_distributed_integral = Some(closed_state_arc_distributed);
            newton_agg.coll_v_distributed_integral = Some(closed_state_coll_v_distributed);
            newton_agg.rift_v_distributed_integral = Some(closed_state_rift_v_distributed);
            newton_agg.spread_distributed_integral = Some(closed_state_spread_distributed);
        }
    }

    // Step 7 — slab-pull diagnostics. Populated only when the
    // slab pipeline actually fired; `Disabled` leaves every field
    // at `None` / empty so the report can structurally skip the
    // Step 7 section.
    if let SlabPullConfig::Enabled { sp, tau_slab, k_slab_accum, .. } = cfg.slab_pull {
        newton_agg.sp_diagnostic = Some(sp);
        newton_agg.tau_slab_diagnostic = Some(tau_slab);
        newton_agg.k_slab_accum_diagnostic = Some(k_slab_accum);
        newton_agg.slab_m_mean_series = slab_m_mean_series;
        newton_agg.slab_m_max_series = slab_m_max_series;
        newton_agg.peak_f_slab_run = Some(peak_f_slab_run);
        newton_agg.peak_f_gpe_run = Some(peak_f_gpe_run);
        if f_slab_to_f_gpe_ratio_samples > 0 {
            newton_agg.f_slab_to_f_gpe_ratio_mean =
                Some(f_slab_to_f_gpe_ratio_sum / f_slab_to_f_gpe_ratio_samples as f64);
        }
    }

    // Step 8 — mantle forcing diagnostics. Populated only when
    // the mantle pipeline actually fired.
    if let MantleConfig::Enabled { mf, coupling, num_modes, seed, .. } = cfg.mantle {
        newton_agg.mf_diagnostic = Some(mf);
        newton_agg.coupling_diagnostic = Some(coupling);
        newton_agg.mantle_num_modes = Some(num_modes);
        newton_agg.mantle_seed = Some(seed);
        newton_agg.peak_v_mantle_pattern = Some(peak_v_mantle_pattern);
        newton_agg.peak_v_solved_mantle_run = Some(peak_v_solved_run);
        if alignment_samples > 0 {
            newton_agg.v_solved_to_v_mantle_alignment =
                Some(alignment_sum / alignment_samples as f64);
        }
        newton_agg.peak_f_mantle_run = Some(peak_f_mantle_run);
        if f_mantle_ratio_samples > 0 {
            newton_agg.f_mantle_to_f_gpe_ratio_mean =
                Some(f_mantle_to_f_gpe_ratio_sum / f_mantle_ratio_samples as f64);
            if is_slab_enabled {
                newton_agg.f_mantle_to_f_slab_ratio_mean =
                    Some(f_mantle_to_f_slab_ratio_sum / f_mantle_ratio_samples as f64);
            }
        }
        newton_agg.epsilon_ii_max_to_floor_ratio = Some(eps_ii_max_to_floor_ratio_run);
        newton_agg.div_v_mantle_max = Some(div_v_mantle_max_run);
        // peak_f_gpe_run is already populated above when slab ran
        // this step too. If slab was disabled, populate it here as
        // a fallback so the Step 8 section can render the ratios.
        if !is_slab_enabled {
            newton_agg.peak_f_gpe_run = Some(peak_f_gpe_run);
        }
    }

    let cg_iter_mean = newton_agg.cg_iters_per_newton_mean();
    let cg_iter_max = newton_agg.cg_iters_per_newton_max();
    let kappa = condition_number_estimate(cg_iter_mean.round() as usize, cfg.newton_cfg.linear_tol);

    // Step 5 — bubble the run-level yielding fraction into the top-
    // level Metrics slot, and encode the boundary-type diversity as
    // a dormant-metric payload (the `BoundaryTypeCounts` variant
    // pre-dated the Step-5 design and does not carry the four mech
    // classes; we encode the diversity as a count under `subduction`
    // when appropriate and also expose the scalar via `newton_agg`).
    let yielding_cf_top_level = newton_agg.yielding_cell_fraction_max;
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
        cg_iter_histogram: IterationHistogram::from_samples(&newton_agg.cg_iters_per_newton_step),
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
        yielding_cell_fraction: yielding_cf_top_level,
        cratonic_stability: None,
        newton_outcome_distribution: None,
        age_field_stats: None,
    };

    let continuation_str = format!("{:?}", cfg.preset.continuation.n_steps);
    let force_name = cfg.force.name();
    let force_detail = match cfg.force_kind {
        ForceKind::Gpe => {
            format!("GpeForce (Ar = {:.3} from scales)", Scales::default().argand_number(),)
        }
        ForceKind::Sinusoidal => {
            format!("SinusoidalForce(ε = {}, Lx = {})", cfg.sinusoidal_amplitude, cfg.domain_lx,)
        }
    };
    let config_dump = SolverConfigDump {
        formulation: "thin viscous sheet (elliptic, no pressure) with power-law rheology".into(),
        discretization: "MAC staggered (v face / η S cell-centre / ε̇_xy corner)".into(),
        eta_averaging: "arithmetic 4-point at corners (see operator.rs)".into(),
        preconditioner: "velocity Jacobi (Picard-block diagonal), null-space wrapped".into(),
        gauge_fixing: "mean(vx), mean(vy) projected before & after every M⁻¹ + post-solve".into(),
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
        basal_drag_config: cfg.basal_drag.describe(),
        boundary_config: cfg.boundary.describe(),
        boundary_layout_name: cfg.boundary_layout_name.clone(),
        slab_pull_config: cfg.slab_pull.describe(),
        mantle_config: cfg.mantle.describe(),
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
