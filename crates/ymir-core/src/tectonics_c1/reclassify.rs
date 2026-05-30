//! C1 facade for the workflow-Phase-A reclassification primitive
//! (Issue #137 Viz-0 Stage A bug fix).
//!
//! ## What this wrapper exposes
//!
//! [`c1_reclassify_plate_type`] is a thin public facade over the
//! internal
//! [`crate::tectonics_v2::workflow::phase_a_common::reclassify_inplace`]
//! (which is `pub(crate)`). The facade keeps viz-layer code (in
//! `ymir-viz/src/bridge/c1/snapshot.rs`) from reaching into the
//! workflow internals directly — viz calls
//! `tectonics_c1::c1_reclassify_plate_type` instead.
//!
//! ## Why a per-snapshot reclassify exists at all
//!
//! [`crate::tectonics_c1::time_loop::run_with_closures`] advects
//! `S̃` and `age` per-step but does NOT call `reclassify_inplace`
//! — that step lives end-of-cycle in the workflow Phase A wrapper
//! (`apply_post_tectonic`). The C1 viz worker (Issue #137) uses
//! `run_with_closures` directly (not the workflow wrapper), so
//! `state.plate_type` is effectively static through a run (modulo
//! Track D's ~20 subduction reassignments / 300 steps).
//!
//! Stage A user-feedback (Issue #137) surfaced that this produces a
//! visibly frozen coastline in the live altitude view even when
//! `S̃` is advecting. To restore the migrating coast on the
//! display, the C1 viz `C1Snapshot::from_state` runs reclassify
//! per snapshot on a **temp copy** of `plate_type` — the
//! simulation's actual `state.plate_type` is unchanged.
//!
//! ## Honest trade-off (snapshot-only, NOT influencing the sim)
//!
//! - **FIX**: the displayed coast (altitude view) migrates with
//!   `S̃` advection — the visual bug is resolved.
//! - **DOES NOT FIX**: the simulation is NOT reclassified. Track D
//!   closures (subduction's plate_type filter, etc.) continue to
//!   see the pre-reclassify `state.plate_type` (init-time + Track
//!   D ~20 reassignments).
//! - **Consequence**: a `run_with_closures` + viz-snapshot-reclassify
//!   run and a full `run_phase_a_cycle_c1(Enabled)` run would
//!   DIVERGE over time — Track D sees different plate_type inputs
//!   between the two paths. The viz shows the qualitative coast
//!   evolution, NOT the exact full-Phase-A trajectory.
//!
//! Acceptable for Viz-0 intuition-qualitative use case. The
//! per-step simulation-influencing reclassify is filed as
//! Viz-0-bis #6 (re-validate Track D event counts under that
//! regime).

use crate::tectonics_v2::boundaries::plate_type::PlateTypeField;
use crate::tectonics_v2::field::Field2D;

/// Snapshot-only reclassification of `plate_type` per the
/// workflow Phase A sea-level threshold (`s_cell > sea_level_ref
/// → Continental`, else `→ Oceanic`).
///
/// Delegates to the internal workflow function
/// `tectonics_v2::workflow::phase_a_common::reclassify_inplace`.
/// See the module docstring for the snapshot-only trade-off.
///
/// Caller is responsible for computing `sea_level_ref` via
/// `tectonics_v2::workflow::phase_a_common::compute_sea_level_ref_s_space`
/// (or equivalent S̃-space sea-level derivation).
pub fn c1_reclassify_plate_type(
    plate_type: &mut PlateTypeField,
    s: &Field2D,
    sea_level_ref: f64,
) {
    crate::tectonics_v2::workflow::phase_a_common::reclassify_inplace(
        plate_type,
        s,
        sea_level_ref,
    );
}
