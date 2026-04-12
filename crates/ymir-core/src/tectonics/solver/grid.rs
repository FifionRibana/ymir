//! Staggered (MAC) grid for the thin viscous sheet solver.

use super::field::{Field2D, PeriodicIndex};

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
            idx: PeriodicIndex::new(n),
        }
    }
}
