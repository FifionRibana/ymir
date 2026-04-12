//! Configuration for the thin viscous sheet tectonic solver.

/// Which nonlinear solver to use for the velocity solve.
#[derive(Clone, Copy, Debug, Default)]
pub enum NonlinearSolver {
    #[default]
    Picard,
    Newton,
}

/// Configuration for Picard (fixed-point) iteration of the nonlinear Stokes solve.
pub struct PicardConfig {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub relaxation: f64,
    pub cg_max_iter: usize,
    pub cg_tolerance: f64,
    pub strain_rate_min: f64,
    pub power_law_n: f64,
}

impl Default for PicardConfig {
    fn default() -> Self {
        Self {
            max_iterations: 30,
            tolerance: 1e-4,
            relaxation: 0.7,
            cg_max_iter: 500,
            cg_tolerance: 1e-8,
            strain_rate_min: 1e-6,
            power_law_n: 3.0,
        }
    }
}

/// Configuration for Newton-Krylov (JFNK) iteration.
pub struct NewtonConfig {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub bicgstab_max_iter: usize,
    pub bicgstab_tolerance: f64,
    pub fd_epsilon_scale: f64,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            tolerance: 1e-4,
            bicgstab_max_iter: 500,
            bicgstab_tolerance: 1e-6,
            fd_epsilon_scale: 1e-7,
        }
    }
}

/// Top-level configuration for the tectonic simulation.
pub struct TectonicsConfig {
    pub num_timesteps: usize,
    pub gravity_factor: f64,
    pub cfl_factor: f64,
    pub s_min: f64,
    pub s_max: f64,
    pub nonlinear_solver: NonlinearSolver,
    pub picard: PicardConfig,
    pub newton: NewtonConfig,
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
        }
    }
}
