//! Issue #117 — retired body-force implementations.
//!
//! `SlabPullForce` and `MantleForce` were the two `BodyForce` impls that
//! coupled to Stokes-driven mechanisms (slab-pull from `_attic/slab/`,
//! mantle from `_attic/mantle/`). They moved here alongside their state
//! modules. The preserved force implementations (`GpeForce`,
//! `SinusoidalForce`, `ForceSum`) and the `BodyForce` trait stay in
//! `tectonics_v2/forcing/` since they have no Stokes coupling.

#![cfg(feature = "v2_legacy")]

pub mod mantle_force;
pub mod slab_pull;

pub use mantle_force::MantleForce;
pub use slab_pull::SlabPullForce;
