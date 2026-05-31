//! Per-step diagnostic stats aggregated from the Track D closures
//! (Issue #137 Viz-D0 Option B).
//!
//! `C1StepStats` bundles the four per-step return values from
//! `apply_subduction_step`, `apply_accretion_step`,
//! `apply_rifting_thinning`, and `apply_rifting_split` into a
//! single field on [`crate::tectonics_c1::state::C1State`]
//! (`C1State.last_step_stats`). The time loop updates it just
//! before invoking the `on_step` callback, so any consumer
//! reading the post-step state sees the matching stats.
//!
//! ## Why a field, not a callback parameter
//!
//! Adding a fifth parameter to `on_step` would break every
//! existing caller of `run_with_closures` (12 test files + the
//! workflow). The field-on-state pattern keeps the callback
//! signature stable while exposing the diagnostic data to
//! consumers that need it (the Viz-0 bridge reads
//! `state.last_step_stats` in its `C1Snapshot::from_state`
//! adapter).
//!
//! ## 9th bit-identical decomposition contract
//!
//! `C1StepStats` lives OUTSIDE the bit-identical decomposition
//! contract — `c1_phase_a_decomposes_into_closures_then_post_tectonic`
//! compares `state.s` byte-for-byte, not the full `C1State`.
//! Adding a diagnostic field on the state struct does NOT
//! affect the `s` buffer in either decomposition path.
//!
//! ## Default = all-zero
//!
//! When Track D closures are disabled (Phase 1.x / Track A/B
//! regression mode), each `apply_*_step` returns its
//! `Default::default()` stats early. The captured
//! `state.last_step_stats` is therefore all-zero — diagnostic
//! consumers reading the field on a Track-D-disabled run see no
//! events, which is correct.

use crate::tectonics_c1::closures::accretion::AccretionStats;
use crate::tectonics_c1::closures::rifting::{RiftingSplitStats, RiftingThinningStats};
use crate::tectonics_c1::closures::subduction::SubductionStats;

/// Bundle of per-step Track D diagnostics. Captured by
/// `run_with_closures` into `C1State.last_step_stats` just
/// before the per-step `on_step` callback fires.
///
/// `Clone + Debug + Default` only — NOT `Copy` because
/// [`RiftingSplitStats.new_plate_ids_created: Vec<u16>`] cannot
/// be Copy. The field is mutated in place each step; no Copy
/// semantics required.
#[derive(Clone, Debug, Default)]
pub struct C1StepStats {
    pub subduction: SubductionStats,
    pub accretion: AccretionStats,
    pub rifting_thinning: RiftingThinningStats,
    pub rifting_split: RiftingSplitStats,
}
