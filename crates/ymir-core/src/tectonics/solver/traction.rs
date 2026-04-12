//! Plate traction fields for the thin viscous sheet solver.

use super::field::Field2D;

/// Traction field imposed by plate motions.
pub struct PlateField {
    pub tx: Field2D,
    pub ty: Field2D,
}

impl PlateField {
    /// No traction — for testing GPE-only dynamics.
    pub fn zero(n: usize) -> Self {
        Self {
            tx: Field2D::new(n),
            ty: Field2D::new(n),
        }
    }

    /// Uniform traction everywhere.
    pub fn uniform(n: usize, tx: f64, ty: f64) -> Self {
        Self {
            tx: Field2D::filled(n, tx),
            ty: Field2D::filled(n, ty),
        }
    }

    /// Two plates converging: left half pushes right (+speed), right half pushes left (-speed).
    pub fn two_plates_convergent(n: usize, speed: f64) -> Self {
        let mut tx = Field2D::new(n);
        let ty = Field2D::new(n);
        let mid = n / 2;
        for j in 0..n {
            for i in 0..n {
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
    pub fn two_plates_divergent(n: usize, speed: f64) -> Self {
        let mut tx = Field2D::new(n);
        let ty = Field2D::new(n);
        let mid = n / 2;
        for j in 0..n {
            for i in 0..n {
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
