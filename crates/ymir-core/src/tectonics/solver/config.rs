//! Configuration for the thin viscous sheet tectonic solver.

/// Which nonlinear solver to use for the velocity solve.
#[derive(Clone, Copy, Debug, Default)]
pub enum NonlinearSolver {
    #[default]
    Picard,
    Newton,
}

/// Configuration for Picard (fixed-point) iteration of the nonlinear Stokes solve.
#[derive(Clone)]
pub struct PicardConfig {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub relaxation: f64,
    pub cg_max_iter: usize,
    pub cg_tolerance: f64,
    pub strain_rate_min: f64,
    pub power_law_n: f64,
    /// Minimum viscosity (prevents zero-viscosity zones).
    pub eta_min: f64,
    /// Maximum viscosity (prevents infinitely rigid zones).
    pub eta_max: f64,
}

impl Default for PicardConfig {
    fn default() -> Self {
        Self {
            max_iterations: 30,
            tolerance: 1e-4,
            relaxation: 0.7,
            cg_max_iter: 500,
            cg_tolerance: 1e-8,
            strain_rate_min: 1e-3,
            power_law_n: 3.0,
            eta_min: 1e-3,
            eta_max: 1e4,
        }
    }
}

/// Preconditioner type for linear solves.
#[derive(Clone, Copy, Debug, Default)]
pub enum Preconditioner {
    #[default]
    Jacobi,
    /// Symmetric Successive Over-Relaxation with parameter omega in (0, 2).
    Ssor { omega: f64 },
}

/// Configuration for JFNK (Jacobian-Free Newton-Krylov) iteration.
#[derive(Clone)]
pub struct NewtonConfig {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub cg_max_iter: usize,
    pub cg_tolerance: f64,
    /// Scale factor for the finite-difference step in the JFNK operator.
    pub fd_epsilon_scale: f64,
    pub preconditioner: Preconditioner,
    /// Use inexact Newton (Eisenstat-Walker adaptive inner tolerance).
    pub inexact: bool,
    /// Tolerance on the relative velocity increment for the state-based
    /// convergence criterion. Convergence is accepted on state if
    /// |Δv| / |v| < state_tolerance AND the residual is on a downward
    /// trend.
    pub state_tolerance: f64,
    /// Number of recent iterations to consider for the residual trend
    /// analysis. The trend is descending if the residual at iteration k
    /// is less than the residual at iteration k - trend_window.
    pub trend_window: usize,
    /// Cosine threshold below which two consecutive Newton steps are
    /// considered anti-aligned (oscillation indicator). Range (-1, 0).
    /// Two consecutive iterations below this threshold trigger the
    /// Oscillation outcome.
    pub oscillation_cosine_threshold: f64,
    /// Minimum number of Newton iterations before the state-based
    /// criterion or oscillation detection can fire. Prevents premature
    /// classification on the first few iterations where signals are noisy.
    pub min_iterations_before_classification: usize,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            tolerance: 5e-2,
            cg_max_iter: 500,
            cg_tolerance: 1e-6,
            fd_epsilon_scale: 1e-7,
            preconditioner: Preconditioner::default(),
            inexact: true,
            state_tolerance: 1e-4,
            trend_window: 3,
            oscillation_cosine_threshold: -0.5,
            min_iterations_before_classification: 3,
        }
    }
}

/// Configuration for viscosity continuation (ramp-up from linear to power-law).
#[derive(Clone, Debug)]
pub struct ContinuationConfig {
    /// Whether to use continuation. If false, solve directly with power_law_n.
    pub enabled: bool,
    /// Sequence of n exponents to solve through.
    pub n_steps: Vec<f64>,
    /// Ramp ε_min from this high value (easy) down to PicardConfig::strain_rate_min.
    /// If None, use PicardConfig::strain_rate_min for all steps.
    pub eps_min_start: Option<f64>,
}

impl Default for ContinuationConfig {
    fn default() -> Self {
        Self { enabled: true, n_steps: vec![1.0, 1.5, 2.0, 2.5, 3.0], eps_min_start: Some(1e-2) }
    }
}

/// Top-level configuration for the tectonic simulation.
#[derive(Clone)]
pub struct TectonicsConfig {
    pub num_timesteps: usize,
    pub gravity_factor: f64,
    pub cfl_factor: f64,
    pub s_min: f64,
    pub s_max: f64,
    pub nonlinear_solver: NonlinearSolver,
    pub picard: PicardConfig,
    pub newton: NewtonConfig,
    pub continuation: ContinuationConfig,
    pub boundaries: super::super::boundaries::BoundaryConfig,
    /// Whether plate boundaries are dynamically recomputed each timestep.
    /// When true, plate seeds are advected with the velocity field and the
    /// Voronoï partition is rebuilt, allowing plates to shrink and disappear.
    pub dynamic_boundaries: bool,
    /// Cratonic rigidity: spatial viscosity variation within continental plates.
    pub cratonic: CratonicConfig,
    /// Plastic yielding with strain weakening.
    pub yielding: YieldingConfig,
    /// Basal friction coefficient (mantle drag). Resists horizontal
    /// crustal motion proportionally to velocity × thickness.
    /// 0.0 = no friction. Typical range: 0.1-10.0. Default: 1.0.
    pub basal_friction: f64,
    /// Mantle convection flow configuration (continuous plate driving).
    pub mantle: super::super::mantle::MantleConfig,
    /// Conservative mass recycling configuration.
    pub recycling: super::super::recycling::RecyclingConfig,
}

/// Configuration for cratonic rigidity (spatial viscosity variation).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CratonicConfig {
    /// Enable spatial viscosity variation. Default: true.
    pub enabled: bool,
    /// Maximum viscosity multiplier at the center of continental plates.
    /// 1.0 = no effect. Range: 1.0-100.0, default: 20.0
    pub max_factor: f64,
    /// Controls how fast rigidity decays from center to edge.
    /// 1.0 = linear, 2.0 = quadratic, 3.0 = cubic.
    /// Range: 0.5-4.0, default: 2.0
    pub decay_power: f64,
}

impl Default for CratonicConfig {
    fn default() -> Self {
        Self { enabled: true, max_factor: 5.0, decay_power: 1.3 }
    }
}

/// Configuration for plastic yielding (Drucker-Prager-like).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct YieldingConfig {
    /// Enable plastic yielding. Default: true.
    pub enabled: bool,
    /// Base yield stress (solver units). When τ = 2ηε̇ exceeds this,
    /// viscosity is reduced so τ = τ_yield. Range: 1.0-500.0, default: 50.0
    pub yield_stress: f64,
    /// Enable strain weakening. Default: true.
    pub weakening_enabled: bool,
    /// Max fractional reduction of yield stress from accumulated plastic strain.
    /// 0.0 = none, 0.8 = yield stress can drop to 20%. Range: 0.0-0.9, default: 0.5
    pub weakening_fraction: f64,
    /// Reference plastic strain for full weakening. Range: 0.1-10.0, default: 1.0
    pub weakening_strain_ref: f64,
    /// Healing rate: plastic strain decreases over time (thermal recovery).
    /// 0.0 = no healing. Rate per timestep (scaled by dt). Default: 0.0
    pub healing_rate: f64,
}

impl Default for YieldingConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            yield_stress: 50.0,
            weakening_enabled: true,
            weakening_fraction: 0.5,
            weakening_strain_ref: 1.0,
            healing_rate: 0.0,
        }
    }
}

impl Default for TectonicsConfig {
    fn default() -> Self {
        Self {
            num_timesteps: 200,
            gravity_factor: 1.0,
            cfl_factor: 0.5,
            s_min: 0.1,
            s_max: 2.5,
            nonlinear_solver: NonlinearSolver::default(),
            picard: PicardConfig::default(),
            newton: NewtonConfig::default(),
            continuation: ContinuationConfig::default(),
            boundaries: Default::default(),
            dynamic_boundaries: true,
            cratonic: CratonicConfig::default(),
            yielding: YieldingConfig::default(),
            basal_friction: 0.05,
            mantle: Default::default(),
            recycling: Default::default(),
        }
    }
}
