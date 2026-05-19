//! Per-cell state for the C1 lightweight tectonics model.
//!
//! Unlike v2, C1 does not solve a momentum equation. Velocity is
//! prescribed per-plate (see [`super::kinematics`]), and the state
//! evolves by advection plus future source terms.
//!
//! ## Field type rationale
//!
//! The v2 `Field2D` is `f64`-specialised (no `Field2D<T>` generic
//! exists in the tree, see comment in
//! `crates/ymir-core/src/tectonics_v2/boundaries/plate_type.rs:23-26`).
//! C1 follows the same convention:
//!
//! - `S̃` and `age` use [`Field2D`] directly.
//! - `plate_id` uses the v2 [`PlateIdField`] (internally `Vec<u16>`).
//! - `plate_type` uses the v2 [`PlateTypeField`].
//! - `cratonic_mask` is a small local [`BoolField`] since neither
//!   `Field2D` nor any v2 type covers `bool` and C1's mask is
//!   binary by design (§4.4: the v2 smoothstep amplification
//!   factor retires).

use crate::tectonics_v2::boundaries::plate_type::PlateTypeField;
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::PlateIdField;

/// Binary cratonic mask, row-major shape `nx × ny` matching
/// [`Field2D`]. Equivalent semantics to the v2 cratonic factor
/// thresholded at the binary cut-off, with the smoothstep
/// amplification factor explicitly retired (§4.4 design doc).
#[derive(Clone, Debug)]
pub struct BoolField {
    nx: usize,
    ny: usize,
    data: Vec<bool>,
}

impl BoolField {
    pub fn filled(nx: usize, ny: usize, value: bool) -> Self {
        Self { nx, ny, data: vec![value; nx * ny] }
    }

    pub fn from_vec(nx: usize, ny: usize, data: Vec<bool>) -> Self {
        assert_eq!(
            data.len(),
            nx * ny,
            "BoolField::from_vec: data length {} != nx*ny = {}*{} = {}",
            data.len(),
            nx,
            ny,
            nx * ny
        );
        Self { nx, ny, data }
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
    pub fn get(&self, i: usize, j: usize) -> bool {
        self.data[j * self.nx + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: bool) {
        self.data[j * self.nx + i] = value;
    }

    pub fn data(&self) -> &[bool] {
        &self.data
    }
}

/// Full C1 cell state.
///
/// Carries the dynamical fields (`s`, `age`) and the static-ish
/// classification fields (`plate_id`, `plate_type`, `cratonic_mask`).
/// Phase 1.1 evolves `s` and `age` by advection only.
pub struct C1State {
    /// Crust thickness, non-dimensional, per cell. Advected.
    pub s: Field2D,
    /// Geological age field, non-dimensional, per cell. Advected.
    pub age: Field2D,
    /// Voronoï plate index per cell. Static in Phase 1.1
    /// (boundary evolution lands in Phase 2).
    pub plate_id: PlateIdField,
    /// Continental / oceanic classification per cell. Static in
    /// Phase 1.1.
    pub plate_type: PlateTypeField,
    /// Binary cratonic mask per cell. Static in Phase 1.1; cells
    /// inside the mask transport rigidly with their plate (the
    /// per-cell velocity is the plate's velocity unchanged, same
    /// behaviour as non-cratonic cells in this stage).
    pub cratonic_mask: BoolField,
    /// Number of distinct plates in `plate_id`. Cached at init for
    /// fast access from the time loop and kinematics builder.
    pub num_plates: usize,
}

impl C1State {
    pub fn nx(&self) -> usize {
        self.s.nx()
    }

    pub fn ny(&self) -> usize {
        self.s.ny()
    }
}
