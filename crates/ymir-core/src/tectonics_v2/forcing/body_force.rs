//! `BodyForce` trait and the context structs it operates on.
//!
//! Step-0/1 carried a pointwise trait (`fx(x, y, s)`, `fy(x, y, s)`)
//! consumed by a `sample_to_faces` helper. Step 2 introduces
//! composite forces (`ForceSum`) and the `GpeForce` term, both of
//! which need access to the full thickness field at once — a
//! pointwise interface would force them to reassemble lookup state
//! per call.
//!
//! The new trait takes a whole-grid [`SimulationState`] and
//! accumulates into a mutable [`VectorField`]. Accumulation (not
//! assignment) is the defining contract: each term adds its
//! contribution to whatever the caller placed in `out`. Callers are
//! responsible for zeroing `out` before the first term.

use super::super::field::{Field2D, PeriodicIndex};

/// Read-only view of the solver state that forcing terms see.
///
/// Step 2 only needs grid geometry, periodic indexing and the
/// current thickness field. Future steps will extend this (velocity,
/// temperature, plate markers, …) as new terms come online;
/// `BodyForce` implementations MUST ignore fields they don't use so
/// this grows backwards-compatibly.
pub struct SimulationState<'a> {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub idx_x: &'a PeriodicIndex,
    pub idx_y: &'a PeriodicIndex,
    /// Current crustal thickness at cell centres.
    pub s: &'a Field2D,
}

/// Mutable view of the two face-centred force components.
///
/// Each `Field2D` is sized `nx × ny` with the MAC convention:
/// `fx[i,j]` lives at `(i·dx, (j+0.5)·dy)` (left vertical face of
/// cell `(i, j)`) and `fy[i,j]` at `((i+0.5)·dx, j·dy)` (bottom
/// horizontal face of cell `(i, j)`).
pub struct VectorField<'a> {
    pub fx: &'a mut Field2D,
    pub fy: &'a mut Field2D,
}

impl<'a> VectorField<'a> {
    /// Zero both components in place. Harness calls this before
    /// the first `accumulate` of each timestep.
    pub fn zero(&mut self) {
        for v in self.fx.data_mut().iter_mut() {
            *v = 0.0;
        }
        for v in self.fy.data_mut().iter_mut() {
            *v = 0.0;
        }
    }
}

/// A per-term body-force contribution to the momentum equation RHS.
///
/// # Contract
///
/// Implementations **add** their contribution to `out` — they must
/// not assume `out` starts at zero, and must not overwrite other
/// terms' contributions. This is what makes [`super::ForceSum`]
/// correct: a sum of terms is a sequential call of each term's
/// `accumulate` on the shared output.
pub trait BodyForce: Send + Sync {
    /// Add this term's contribution to `out`.
    fn accumulate(&self, state: &SimulationState, out: &mut VectorField);
    /// Short name used by diagnostic reports and logs.
    fn name(&self) -> &'static str;
}
