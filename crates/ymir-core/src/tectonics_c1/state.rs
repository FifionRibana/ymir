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

use crate::tectonics_c1::stats::C1StepStats;
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
/// Carries the dynamical fields (`s`, `age`, plus Phase 2 Track D
/// `plate_id` / `plate_type` / `kinematics` mutation surfaces) and
/// the slow-evolving `cratonic_mask`.
///
/// **Mutation status** (post-Track D, Issue #132):
///
/// - `s`, `age` — advected every step (Phase 1.1 baseline).
/// - `plate_id`, `plate_type` — mutate per-step under Track D:
///   subduction floor-trigger reassignment, accretion merges,
///   rifting splits. Static when Track D is disabled (Phase 1.x /
///   Track A/B regression mode).
/// - `cratonic_mask` — built at init by BFS over the
///   Continental-seeded plates; NOT recomputed per step. A
///   one-cycle lag exists when Track D mutates `plate_id`
///   mid-cycle (Issue #132 Stage S Q-S.2 Option (c) accepted
///   trade-off: cratons are ~100 Ma features, lag is ~0.67 Ma).
/// - `num_plates` — count at init time. With Track D enabled the
///   actual count of distinct plate ids in `plate_id` can drop
///   (merges) or grow (rifting splits) during a run; consumers
///   needing the live count should re-scan `plate_id.data()`.
/// - `last_step_stats` — refreshed every step by
///   `run_with_closures` just before the `on_step` callback fires
///   (Viz-D0 Option B, Issue #137).
pub struct C1State {
    /// Crust thickness, non-dimensional, per cell. Advected.
    pub s: Field2D,
    /// Geological age field, non-dimensional, per cell. Advected.
    /// Path 3.A ridge-aligned init (Track B Issue #131); Path 3.B
    /// event-driven `= 0` on rift-spawned cells (Track D Issue #132).
    pub age: Field2D,
    /// Voronoï plate index per cell. Mutates per-step when any
    /// Track D closure is enabled (subduction floor-trigger
    /// reassignment, accretion merge, rifting split). Static
    /// otherwise (Phase 1.x / Track A/B regression mode).
    pub plate_id: PlateIdField,
    /// Continental / oceanic classification per cell. Mutates
    /// per-step under Track D subduction reassignment (Oceanic →
    /// Continental on floor trigger). Accretion + rifting splits
    /// preserve per-cell type. Static under Track D disabled.
    pub plate_type: PlateTypeField,
    /// Binary cratonic mask per cell. Built at init by BFS over
    /// Continental-seeded plates; NOT recomputed per step. Cells
    /// inside the mask transport rigidly with their plate. One-
    /// cycle lag under Track D mutation is the accepted trade-off
    /// (Issue #132 Q-S.2 Option (c)).
    pub cratonic_mask: BoolField,
    /// Number of distinct plates **at init time**. With Track D
    /// enabled, accretion merges + rifting splits can change the
    /// live count; consumers needing live cardinality should
    /// re-scan `plate_id.data()`. Cached at init for fast access
    /// from the time loop and kinematics builder.
    pub num_plates: usize,
    /// Per-step diagnostic stats from the four Track D `apply_*_step`
    /// returns, captured by `run_with_closures` just before
    /// `on_step`. Default = all-zero when Track D closures are
    /// disabled or no event fired (Issue #137 Viz-D0 Option B).
    /// Lives outside the 9th bit-identical decomposition contract.
    pub last_step_stats: C1StepStats,
}

impl C1State {
    pub fn nx(&self) -> usize {
        self.s.nx()
    }

    pub fn ny(&self) -> usize {
        self.s.ny()
    }
}
