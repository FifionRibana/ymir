//! Body-force hierarchy.
//!
//! Step 0/1 lived in a single `forcing.rs`; Step 2 split it into a
//! directory so `GpeForce` and `ForceSum` can land without crowding
//! the one pointwise placeholder that used to live here.
//!
//! The public API is a trait `BodyForce` and a handful of
//! implementations. Composite forces use `ForceSum`.

pub mod body_force;
pub mod force_sum;
pub mod gpe;
pub mod sinusoidal;

pub use body_force::{BodyForce, SimulationState, VectorField};
pub use force_sum::ForceSum;
pub use gpe::GpeForce;
pub use sinusoidal::{SinusoidalForce, ZeroForce};

// Issue #117 — `MantleForce` and `SlabPullForce` retired to
// `_attic/forcing/` (gated by `v2_legacy`). Re-exports here restore
// the old `tectonics_v2::forcing::{mantle_force, slab_pull, …}` paths
// under feature so callers compile bit-identically.
#[cfg(feature = "v2_legacy")]
pub use super::_attic::forcing::{MantleForce, SlabPullForce, mantle_force, slab_pull};
