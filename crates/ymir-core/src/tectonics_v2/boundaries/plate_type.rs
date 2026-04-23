//! Plate-type enum for Step 5.
//!
//! Distinguishes oceanic from continental lithosphere. At Step 5
//! the field is **prescribed statically** (set once at simulation
//! start from a layout generator and not updated) — dynamic
//! reclassification is Step 6.
//!
//! The enum has no `None` variant: every cell is either oceanic or
//! continental. "Non-boundary" cells are distinguished by the
//! companion [`super::boundary_flag::BoundaryFlag`] field, whose
//! `None` variant carries that meaning.

use super::super::field::Field2D;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlateType {
    Oceanic,
    Continental,
}

/// Cell-centred plate-type field.
///
/// Stored as `Vec<PlateType>` in row-major order, shape `nx × ny`
/// matching [`Field2D`]. No `Field2D<T>` generic exists in the
/// current tree (the legacy `Field2D` is `f64`-specialised), so
/// plate/flag fields are their own small type.
#[derive(Clone, Debug)]
pub struct PlateTypeField {
    nx: usize,
    ny: usize,
    data: Vec<PlateType>,
}

impl PlateTypeField {
    pub fn filled(nx: usize, ny: usize, t: PlateType) -> Self {
        Self {
            nx,
            ny,
            data: vec![t; nx * ny],
        }
    }

    pub fn nx(&self) -> usize { self.nx }
    pub fn ny(&self) -> usize { self.ny }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> PlateType {
        self.data[j * self.nx + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, t: PlateType) {
        self.data[j * self.nx + i] = t;
    }

    pub fn data(&self) -> &[PlateType] { &self.data }

    /// Render as an f64 heightmap (0.0 = Oceanic, 1.0 = Continental)
    /// for PNG layout visualisation in the report.
    pub fn to_heightmap(&self) -> Field2D {
        let mut f = Field2D::new(self.nx, self.ny);
        for j in 0..self.ny {
            for i in 0..self.nx {
                let v = match self.get(i, j) {
                    PlateType::Oceanic => 0.0,
                    PlateType::Continental => 1.0,
                };
                f.set(i, j, v);
            }
        }
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filled_has_uniform_type() {
        let f = PlateTypeField::filled(4, 3, PlateType::Oceanic);
        for j in 0..3 {
            for i in 0..4 {
                assert_eq!(f.get(i, j), PlateType::Oceanic);
            }
        }
    }

    #[test]
    fn set_updates_single_cell() {
        let mut f = PlateTypeField::filled(4, 3, PlateType::Continental);
        f.set(2, 1, PlateType::Oceanic);
        for j in 0..3 {
            for i in 0..4 {
                let expected = if (i, j) == (2, 1) {
                    PlateType::Oceanic
                } else {
                    PlateType::Continental
                };
                assert_eq!(f.get(i, j), expected);
            }
        }
    }

    #[test]
    fn heightmap_encodes_plate_type() {
        let mut f = PlateTypeField::filled(2, 2, PlateType::Oceanic);
        f.set(1, 1, PlateType::Continental);
        let hm = f.to_heightmap();
        assert_eq!(hm.get(0, 0), 0.0);
        assert_eq!(hm.get(1, 0), 0.0);
        assert_eq!(hm.get(0, 1), 0.0);
        assert_eq!(hm.get(1, 1), 1.0);
    }
}
