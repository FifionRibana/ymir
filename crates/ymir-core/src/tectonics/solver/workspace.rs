//! Pre-allocated workspace buffers for the solver.

use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::linear_solve::{BiCgStabWorkspace, CgWorkspace};

/// Statistics collected after each timestep.
#[derive(Clone, Debug)]
pub struct StepStats {
    pub max_velocity: f64,
    pub max_thickness: f64,
    pub min_thickness: f64,
    pub picard_iterations: usize,
    pub cg_iterations_last: usize,
    pub dt: f64,
    pub clamp_ratio: f64,
    /// True when the CFL retry loop exhausted all attempts without
    /// satisfying `clamp_ratio < 0.05`. The step is still accepted with
    /// the last attempted dt (smaller than the initial dt_cfl by up to
    /// `2^(MAX_RETRIES - 1)`) — the flag surfaces the degradation so
    /// the UI and callers can react.
    pub cfl_retry_exhausted: bool,
}

impl Default for StepStats {
    fn default() -> Self {
        Self {
            max_velocity: 0.0,
            max_thickness: 0.0,
            min_thickness: 0.0,
            picard_iterations: 0,
            cg_iterations_last: 0,
            dt: 0.0,
            clamp_ratio: 0.0,
            cfl_retry_exhausted: false,
        }
    }
}

/// All temporary buffers needed by the solver, pre-allocated to avoid per-step allocation.
pub struct SolverWorkspace {
    nx: usize,
    ny: usize,
    pub div_flux: Field2D,
    pub eta: Field2D,
    pub strain_rate: Field2D,
    pub v_packed: Vec<f64>,
    pub v_prev: Vec<f64>,
    pub rhs: Vec<f64>,
    pub jacobi_precond: Vec<f64>,
    pub cg: CgWorkspace,
    pub bicgstab: BiCgStabWorkspace,
    pub jfnk_v_pert: Vec<f64>,
    pub jfnk_f_v: Vec<f64>,
    pub jfnk_neg_f: Vec<f64>,
    pub jfnk_delta_v: Vec<f64>,
    pub source_rate: Field2D,
    pub boundary_field: Option<crate::tectonics::boundaries::BoundaryField>,
    pub stats: StepStats,
}

impl SolverWorkspace {
    pub fn new(nx: usize, ny: usize) -> Self {
        let nn2 = 2 * nx * ny;
        Self {
            nx,
            ny,
            div_flux: Field2D::new(nx, ny),
            eta: Field2D::new(nx, ny),
            strain_rate: Field2D::new(nx, ny),
            v_packed: vec![0.0; nn2],
            v_prev: vec![0.0; nn2],
            rhs: vec![0.0; nn2],
            jacobi_precond: vec![1.0; nn2],
            cg: CgWorkspace::new(nn2),
            bicgstab: BiCgStabWorkspace::new(nn2),
            jfnk_v_pert: vec![0.0; nn2],
            jfnk_f_v: vec![0.0; nn2],
            jfnk_neg_f: vec![0.0; nn2],
            jfnk_delta_v: vec![0.0; nn2],
            source_rate: Field2D::new(nx, ny),
            boundary_field: None,
            stats: StepStats::default(),
        }
    }

    pub fn resize_if_needed(&mut self, nx: usize, ny: usize) {
        if self.nx != nx || self.ny != ny {
            *self = Self::new(nx, ny);
        }
    }

    #[inline]
    pub fn nx(&self) -> usize {
        self.nx
    }

    #[inline]
    pub fn ny(&self) -> usize {
        self.ny
    }
}

/// Pack grid velocity fields (vx, vy) into a single flat vector.
/// Layout: first nx*ny elements = vx, next nx*ny = vy.
pub fn pack_velocity(grid: &StaggeredGrid, buf: &mut [f64]) {
    let n2 = grid.nx() * grid.ny();
    buf[..n2].copy_from_slice(grid.vx.data());
    buf[n2..2 * n2].copy_from_slice(grid.vy.data());
}

/// Unpack a flat vector into grid velocity fields.
pub fn unpack_velocity(buf: &[f64], grid: &mut StaggeredGrid) {
    let n2 = grid.nx() * grid.ny();
    grid.vx.data_mut().copy_from_slice(&buf[..n2]);
    grid.vy.data_mut().copy_from_slice(&buf[n2..2 * n2]);
}
