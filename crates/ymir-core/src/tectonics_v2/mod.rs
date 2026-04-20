//! Tectonics v2 — incremental rebuild of the tectonic solver from a
//! nondimensional core.
//!
//! Step 0 ships a linear, constant-viscosity Stokes solver coupled to a
//! passive advected thickness field on a fully periodic toroidal domain.
//! Subsequent steps add power-law rheology, GPE spreading, yielding,
//! basal drag, boundary sources, slab pull, mantle forcing, cratonic
//! immunity, and a geological age field, in that order.
//!
//! The module is strictly isolated from the legacy `tectonics/` module;
//! the only whitelisted imports are `Field2D` and `PeriodicIndex`, which
//! this crate re-exports from the local [`field`] submodule.

pub mod advection;
pub mod diagnostics;
pub mod field;
pub mod forcing;
pub mod scales;
pub mod stokes;
