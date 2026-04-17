//! Staggered (MAC) grid for the thin viscous sheet solver.

use super::field::{Field2D, PeriodicIndex};
use crate::grid::GridF32;

/// Marker-And-Cell staggered grid with periodic boundary conditions.
///
/// - `s[i,j]` lives at cell center `((i+0.5)*dx, (j+0.5)*dx)`
/// - `vx[i,j]` lives at left vertical face `(i*dx, (j+0.5)*dx)`
/// - `vy[i,j]` lives at bottom horizontal face `((i+0.5)*dx, j*dx)`
///
/// With periodic BCs all fields are nx-by-ny. Cell spacing `dx` is
/// isotropic; rectangular domains are obtained by varying the cell
/// counts, not the spacing.
pub struct StaggeredGrid {
    nx: usize,
    ny: usize,
    pub dx: f64,
    pub s: Field2D,
    pub vx: Field2D,
    pub vy: Field2D,
    /// Crustal density at each cell center (kg/m³).
    /// Continental ≈ 2750, Oceanic ≈ 3000.
    pub rho: Field2D,
    /// Spatial viscosity multiplier (cratonic rigidity).
    /// 1.0 = normal, >1.0 = more rigid (continental interior).
    pub eta_multiplier: Field2D,
    /// Accumulated plastic strain at each cell.
    /// Grows when a cell yields; used for strain weakening.
    pub plastic_strain: Field2D,
    idx_x: PeriodicIndex,
    idx_y: PeriodicIndex,
    /// Basal friction coefficient (mantle drag). Set once before solving.
    pub basal_friction: f64,
}

impl StaggeredGrid {
    pub fn new(nx: usize, ny: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            s: Field2D::new(nx, ny),
            vx: Field2D::new(nx, ny),
            vy: Field2D::new(nx, ny),
            rho: Field2D::new(nx, ny),
            eta_multiplier: Field2D::filled(nx, ny, 1.0),
            plastic_strain: Field2D::new(nx, ny),
            idx_x: PeriodicIndex::new(nx),
            idx_y: PeriodicIndex::new(ny),
            basal_friction: 0.0,
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

    #[inline]
    pub fn idx_x(&self) -> &PeriodicIndex {
        &self.idx_x
    }

    #[inline]
    pub fn idx_y(&self) -> &PeriodicIndex {
        &self.idx_y
    }

    /// Extract the crustal thickness field as a GridF32 (f64 → f32).
    pub fn thickness_to_grid_f32(&self) -> GridF32 {
        let data: Vec<f32> = self.s.data().iter().map(|&v| v as f32).collect();
        GridF32::from_vec(self.nx, self.ny, data)
    }
}
