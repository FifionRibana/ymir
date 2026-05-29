//! Rifting closure + split mechanism — divergent-boundary
//! continental thinning + "chewing-gum cut" plate split (C1
//! Phase 2 Track D, Issue #132).
//!
//! ## Physics in one paragraph
//!
//! At a divergent continental boundary, sustained extension
//! thins the lithosphere (McKenzie 1978 stretching-factor
//! framework). Once the thinning crosses a critical threshold
//! (`β = 1.4`, ~30 % crustal thinning) AND the extension has
//! been sustained over geological timescales (~50 Ma), the
//! lithosphere fails: a new plate boundary opens, the
//! continental margin separates from its parent, and a nascent
//! rift basin forms with the perpendicular-separation kinematics
//! that drove the rift.
//!
//! C1 implements this with two coupled algorithms:
//!
//! 1. **Thinning closure**
//!    ([`source_term::apply_rifting_thinning`]) — per-step
//!    negative `S̃` source on continental cells classified as
//!    `BoundaryType::Divergent`. Mirror of the Davis-Suppe
//!    orogenic source on convergent boundaries.
//! 2. **Split event** ([`split::apply_rifting_split`]) — fires
//!    when BOTH conditions hold:
//!    (a) `DivergenceTracker.get(a, b) >= split_time_threshold`
//!    (sustained extension); and
//!    (b) `min(S̃) at rift strip < split_thickness_threshold`
//!    (sub-threshold thinness). The "chewing-gum cut" framing
//!    (Q3.2 hybrid two-condition gate).
//!
//! ## Track D Q-decision history applied here
//!
//! - **Q3.2 — hybrid two-condition split** (vs time-only,
//!   vs thickness-only): both conditions required. Either alone
//!   is insufficient. Intuition: "stretched (time) + thinned
//!   (mass)" — Atlantic rifting needed both sustained motion AND
//!   accumulated thinning to mature into ocean basin opening.
//! - **Q3.4 — perpendicular velocity offset** (vs inherited
//!   velocity, vs random perturbation): new plate's velocity =
//!   parent's velocity + perpendicular `(perp_x, perp_y) ×
//!   split_velocity_offset` where the perpendicular is the
//!   right-hand-rule rotation of `v_rel = v_a − v_b`. Gives the
//!   new plate a deterministic separation direction that
//!   continues the rift's opening pattern.
//!
//! ## Track D mutation pattern (third closure)
//!
//! Rifting is the **third C1 closure to mutate `plate_id`**
//! (after subduction's floor-trigger reassignment and
//! accretion's merge). The split adds a NEW plate id, extending
//! `kinematics.velocities` by one — the first closure to grow
//! the plate count. The new pid is allocated as
//! `kinematics.velocities.len()`, so accretion's no-compaction
//! gaps (Stage E2) are preserved without collisions.
//!
//! Path 3.B event-driven `age = 0` extends Track B Path 3.A
//! (init-only) to the rifting event — every cell reassigned to
//! the new plate has its `age` reset, simulating the fresh rift
//! floor.
//!
//! ## Module layout
//!
//! - [`params`] — [`params::RiftingParams`] tunables
//!   (`thinning_rate = 1.0`, `split_time_threshold = 75`,
//!   `split_thickness_threshold = 0.7`, `split_velocity_offset
//!   = 0.005`, `plate_id_cap = 256`).
//! - [`source_term`] — [`source_term::apply_rifting_thinning`]
//!   closure portion + [`source_term::RiftingThinningStats`]
//!   diagnostics + 3 unit tests.
//! - [`split`] — [`split::DivergenceTracker`] symmetric mirror
//!   of [`crate::tectonics_c1::closures::accretion::ConvergenceTracker`],
//!   [`split::apply_rifting_split`] event + [`split::RiftingSplitStats`]
//!   + 5 unit tests (4 split-specific + 1 unified disabled).
//!
//! ## References
//!
//! - McKenzie, D. (1978). Some remarks on the development of
//!   sedimentary basins. *Earth and Planetary Science Letters*
//!   40(1), 25-32. doi:10.1016/0012-821X(78)90071-7 — the
//!   stretching-factor `β = 1.4` reference for the 30 %
//!   continental-thinning threshold used as
//!   `split_thickness_threshold = 0.7`.
//! - Buck, W. R. (1991). Modes of continental lithospheric
//!   extension. *Journal of Geophysical Research* 96(B12),
//!   20161-20178. doi:10.1029/91JB01485 — passive-margin
//!   morphology reference; framework that motivates Q3.4
//!   perpendicular velocity offset.

pub mod params;
pub mod source_term;
pub mod split;

pub use params::RiftingParams;
pub use source_term::{apply_rifting_thinning, RiftingThinningStats};
pub use split::{apply_rifting_split, DivergenceTracker, RiftingSplitStats};
