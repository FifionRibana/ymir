//! Staggered (MAC) grid for the thin viscous sheet solver.

use super::field::{Field2D, PeriodicIndex};
use crate::grid::GridF32;

/// Marker-And-Cell staggered grid with periodic boundary conditions.
///
/// - `s[i,j]` lives at cell center `((i+0.5)*dx, (j+0.5)*dx)`
/// - `vx[i,j]` lives at left vertical face `(i*dx, (j+0.5)*dx)`
/// - `vy[i,j]` lives at bottom horizontal face `((i+0.5)*dx, j*dx)`
///
/// With periodic BCs all three fields are N×N.
pub struct StaggeredGrid {
    pub n: usize,
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
    pub idx: PeriodicIndex,
}

impl StaggeredGrid {
    pub fn new(n: usize, dx: f64) -> Self {
        Self {
            n,
            dx,
            s: Field2D::new(n),
            vx: Field2D::new(n),
            vy: Field2D::new(n),
            rho: Field2D::new(n),
            eta_multiplier: Field2D::filled(n, 1.0),
            idx: PeriodicIndex::new(n),
        }
    }

    /// Extract the crustal thickness field as a GridF32 (f64 → f32).
    pub fn thickness_to_grid_f32(&self) -> GridF32 {
        let data: Vec<f32> = self.s.data().iter().map(|&v| v as f32).collect();
        GridF32::from_vec(self.n, self.n, data)
    }
}
