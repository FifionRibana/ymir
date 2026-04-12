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
}

impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            tolerance: 1e-4,
            cg_max_iter: 500,
            cg_tolerance: 1e-6,
            fd_epsilon_scale: 1e-7,
            preconditioner: Preconditioner::default(),
            inexact: true,
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
        Self {
            enabled: true,
            n_steps: vec![1.0, 1.5, 2.0, 2.5, 3.0],
            eps_min_start: Some(1e-2),
        }
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
        }
    }
}
