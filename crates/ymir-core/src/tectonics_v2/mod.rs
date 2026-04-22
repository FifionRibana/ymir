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
//! the whitelisted imports are `Field2D`, `PeriodicIndex`, and
//! `soft_min_harmonic` (the latter re-exported by the [`rheology`]
//! submodule for the plastic-branch blend introduced at Step 3).

pub mod advection;
pub mod basal_drag;
pub mod diagnostics;
pub mod field;
pub mod forcing;
pub mod presets;
pub mod rheology;
pub mod scales;
pub mod stokes;

pub use rheology::soft_min_harmonic;
