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
//! the whitelisted imports are `Field2D`, `PeriodicIndex`,
//! `soft_min_harmonic` (re-exported by the [`rheology`] submodule for
//! the plastic-branch blend introduced at Step 3), and `RecyclingBuffer`
//! (re-exported by [`recycling`] for Step 6's delayed mantle
//! recycling; the wrapper [`recycling::DelayedRecycler`] adds rollover
//! and fill tracking on top).

pub mod advection;
pub mod age_field;
pub mod boundaries;
pub mod boundary_detection;
pub mod cancel;
pub mod cratonic;
pub mod diagnostics;
pub mod field;
pub mod forcing;
pub mod init;
pub mod plate_kinematic;
pub mod recycling;
pub mod scales;
pub mod voronoi;
pub mod workflow;

// Issue #117 — Stokes-coupled subtree retired to `_attic/` (gated by
// Cargo feature `v2_legacy`). Re-exports below restore the old module
// paths under feature so callers (harness, bridge build_config, R4/R5b/
// R6/R7.A tests) compile bit-identically once the feature is on. The
// default build sees no retired modules at all.
#[cfg(feature = "v2_legacy")]
pub mod _attic;

#[cfg(feature = "v2_legacy")]
pub use _attic::{basal_drag, mantle, presets, rheology, slab, stokes};

#[cfg(feature = "v2_legacy")]
pub use _attic::rheology::soft_min_harmonic;
