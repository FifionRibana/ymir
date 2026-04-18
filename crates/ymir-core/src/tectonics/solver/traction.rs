//! Plate traction fields for the thin viscous sheet solver.

use super::field::Field2D;

/// Traction field imposed by plate motions.
#[derive(Clone)]
pub struct TractionField {
    pub tx: Field2D,
    pub ty: Field2D,
}

impl TractionField {
    /// No traction — for testing GPE-only dynamics.
    pub fn zero(nx: usize, ny: usize) -> Self {
        Self { tx: Field2D::new(nx, ny), ty: Field2D::new(nx, ny) }
    }

    /// Uniform traction everywhere.
    pub fn uniform(nx: usize, ny: usize, tx: f64, ty: f64) -> Self {
        Self {
            tx: Field2D::filled(nx, ny, tx),
            ty: Field2D::filled(nx, ny, ty),
        }
    }

    /// Two plates converging: left half pushes right (+speed), right half pushes left (-speed).
    pub fn two_plates_convergent(nx: usize, ny: usize, speed: f64) -> Self {
        let mut tx = Field2D::new(nx, ny);
        let ty = Field2D::new(nx, ny);
        let mid = nx / 2;
        for j in 0..ny {
            for i in 0..nx {
                if i < mid {
                    tx.set(i, j, speed);
                } else {
                    tx.set(i, j, -speed);
                }
            }
        }
        Self { tx, ty }
    }

    /// Two plates diverging: left half pushes left (-speed), right half pushes right (+speed).
    pub fn two_plates_divergent(nx: usize, ny: usize, speed: f64) -> Self {
        let mut tx = Field2D::new(nx, ny);
        let ty = Field2D::new(nx, ny);
        let mid = nx / 2;
        for j in 0..ny {
            for i in 0..nx {
                if i < mid {
                    tx.set(i, j, -speed);
                } else {
                    tx.set(i, j, speed);
                }
            }
        }
        Self { tx, ty }
    }
}
