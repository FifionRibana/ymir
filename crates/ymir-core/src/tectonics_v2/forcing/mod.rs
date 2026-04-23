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
pub mod slab_pull;

pub use body_force::{BodyForce, SimulationState, VectorField};
pub use force_sum::ForceSum;
pub use gpe::GpeForce;
pub use sinusoidal::{SinusoidalForce, ZeroForce};
pub use slab_pull::SlabPullForce;
