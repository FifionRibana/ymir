//! Accretion mechanism — sustained-convergence plate-id merge
//! (C1 Phase 2 Track D, Issue #132).
//!
//! ## Physics in one paragraph
//!
//! When two plates converge for a geologically meaningful duration
//! (Indian-Asian collision matured over ~50 Ma), the boundary
//! between them ceases to function as an active subduction or
//! orogenic margin and the two lithospheric blocks behave as a
//! single welded plate (Coney-Jones-Monger 1980 terrane accretion
//! framework). C1 implements this with a single event:
//!
//! 1. Per-pair sustained-convergence counter incremented each
//!    step the plate pair shows a net-convergent boundary verdict
//!    (see [`merge::ConvergenceTracker::update`]).
//! 2. When the counter `≥ merge_time_threshold`, the lower-index
//!    plate absorbs the higher-index plate's cells, and the
//!    surviving plate's velocity is the mass-weighted average of
//!    the two pre-merge velocities.
//!
//! Accretion **does not add a thickening source term** to `S̃` —
//! the Davis-Suppe closure (Phase 1.2, §5.1) is already
//! responsible for orogenic morphology at convergent boundaries
//! during the pre-merge accumulation phase. The merge itself only
//! resolves the boundary topology.
//!
//! ## Track D Q-decision history applied here
//!
//! - **Q2.4 — mass-weighted average velocity** (vs simple
//!   arithmetic average, vs momentum-conserving): the winner
//!   plate's new velocity is `v_new = (m_a v_a + m_b v_b) /
//!   (m_a + m_b)` where `m_p = Σ S̃ over plate-p cells`.
//!   Rationale: the smaller plate's pre-merge velocity should
//!   not impose itself on the larger plate's bulk motion;
//!   mass-weighting preserves angular momentum to first order.
//! - **Q3 (revised) — `merge_time_threshold = 50` steps**: the
//!   original design exploration set Q3 = 20. Stage E1 W7
//!   analytical revised to 50 for spurious-merge suppression on
//!   Phase 1.1 random-cycled kinematics. See
//!   [`params::AccretionParams`] docstring for the full revision
//!   rationale.
//!
//! ## Track D mutation pattern
//!
//! Accretion is the **second C1 closure to mutate `plate_id`**
//! (after subduction's floor-trigger reassignment), and the
//! **first to mutate `kinematics.velocities`**. The
//! `PlateKinematics.velocities` field is already `pub`, so no
//! API extension was needed for the kinematics mutation (W4 of
//! the Stage E2 spec — concern dismissed).
//!
//! `plate_type` is preserved per cell (a continental cell stays
//! continental after merge, an oceanic cell stays oceanic). The
//! merge is a `plate_id` re-assignment only.
//!
//! ## No plate-id compaction in E2
//!
//! Per the Stage E2 spec W3, this stage **leaves gaps** in the
//! plate-id space — if plates `0` and `1` merge, the surviving
//! ids are `{0, 2, 3, ...}` with `1` as an unused slot in
//! `kinematics.velocities`. Rifting (Stage E3) allocates new ids
//! starting at the highest currently-used index + 1, so it never
//! collides with a gap. If a future scenario needs the
//! lowest-available index, compaction can land in Stage E4 (the
//! integration stage) without rewriting this module.
//!
//! ## Module layout
//!
//! - [`params`] — [`params::AccretionParams`] tunables
//!   (`enabled = true`, `merge_time_threshold = 50`,
//!   `velocity_merge_method = MassWeightedAverage`).
//! - [`merge`] — [`merge::ConvergenceTracker`] sustained-
//!   convergence accumulator, [`merge::AccretionStats`] +
//!   [`merge::apply_accretion_step`] event + 6 unit tests.
//!
//! ## References
//!
//! - Coney, P. J., Jones, D. L. & Monger, J. W. H. (1980).
//!   Cordilleran suspect terranes. *Nature* 288, 329-333.
//!   doi:10.1038/288329a0

pub mod merge;
pub mod params;

pub use merge::{AccretionStats, ConvergenceTracker, apply_accretion_step};
pub use params::{AccretionParams, VelocityMergeMethod};
