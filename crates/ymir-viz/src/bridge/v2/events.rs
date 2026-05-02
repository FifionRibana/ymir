//! Step 8.6 — events streamed from the v2 solver thread back to the
//! Bevy main thread.
//!
//! `V2FinalState` is a thread-safe, owned snapshot of every raster
//! field a UI / screenshot harness might want to display. It mirrors
//! `ymir_core::tectonics_v2::diagnostics::harness::FinalState` but
//! carries the field payloads as flat `Vec<f64>` row-major buffers
//! (rather than `Field2D`) so consumers don't need to depend on the
//! v2-specific grid type.

use std::time::Duration;

use ymir_core::tectonics_v2::diagnostics::harness::{FinalState, StepProgress};
use ymir_core::tectonics_v2::diagnostics::metrics::Metrics;

use super::spec::V2RunSpec;

#[allow(clippy::large_enum_variant)]
pub enum V2Event {
    /// Emitted as soon as the worker thread accepts a run command.
    /// Carries the resolved spec so the UI can confirm what is
    /// running (especially handy when a preset re-derives some
    /// fields).
    Started { spec: V2RunSpec },
    /// Step 8.6 follow-up — emitted once per completed step. Carries
    /// the step counter, total step count, and a thread-safe owned
    /// snapshot of every raster the UI might render mid-run. The bridge
    /// fires this event from the harness step callback, so it lands at
    /// the same cadence as the simulation steps (typically one event
    /// per ~5-30 s on a 64² mantle-on regime).
    Progress {
        step: usize,
        total: usize,
        peek_state: V2FinalState,
    },
    /// Emitted at end of run with the full final-state snapshot and
    /// the harness-computed metrics. `elapsed` is the wallclock of
    /// the `run_baseline` call alone (not the queue-wait time).
    Completed {
        spec: V2RunSpec,
        final_state: V2FinalState,
        metrics: Box<Metrics>,
        elapsed: Duration,
    },
    /// Emitted on a panic-equivalent failure (currently any
    /// downstream `Err` from the harness, though `run_baseline` itself
    /// does not return a `Result` at Step 8.6 Phase 1; reserved for
    /// future cancellation-aware refactor).
    Failed { error: String },
}

/// Thread-safe owned snapshot of every raster field at end of run.
///
/// The buffers are row-major, shape `nx × ny`, indexed as
/// `data[j * nx + i]`. Optional fields are `None` when the
/// corresponding mechanism was disabled in the run config.
#[derive(Clone, Debug)]
pub struct V2FinalState {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub s_field: Vec<f64>,
    pub vx: Vec<f64>,
    pub vy: Vec<f64>,
    pub strain_rate_invariant: Vec<f64>,
    pub age_field: Option<Vec<f64>>,
    pub cratonic_factor: Option<Vec<f64>>,
    /// `plate_id` cast to `u16`. Useful for plate-boundary overlays.
    pub plate_id: Option<Vec<u16>>,
    /// Cell-centred plate-type tags as small `u8`: `0 = Oceanic`,
    /// `1 = Continental`. Avoids leaking the v2 enum into the bridge
    /// surface so non-rust consumers can read it as a plain raster.
    pub plate_type: Option<Vec<u8>>,
    /// Cell-centred boundary-flag tags as small `u8`:
    /// `0 = None`, `1 = Subduction`, `2 = OceanicSubduction`,
    /// `3 = Rift`, `4 = ContinentalCollision`. Same rationale as
    /// above.
    pub boundary_flag: Option<Vec<u8>>,
}

impl V2FinalState {
    /// Step 8.6 follow-up — build a thread-safe peek-state from the
    /// harness's per-step `StepProgress`. Mirrors `from_harness`'s
    /// shape. The strain-rate invariant is taken from
    /// `peek_strain_ii` (the harness already computes ε̇_II per step
    /// for the rheology pipeline, so the cost is one extra Field2D
    /// clone — negligible).
    pub fn from_step_progress(p: &StepProgress<'_>) -> Self {
        let nx = p.peek_s.nx();
        let ny = p.peek_s.ny();
        V2FinalState {
            nx,
            ny,
            // dx / dy are not threaded into StepProgress; consumers
            // that need world-space metrics fall back on the
            // `Completed` event's final_state. For mid-run viz the
            // raster grid alone is enough.
            dx: 0.0,
            dy: 0.0,
            s_field: p.peek_s.data().to_vec(),
            vx: p.peek_vx.to_vec(),
            vy: p.peek_vy.to_vec(),
            strain_rate_invariant: p.peek_strain_ii.data().to_vec(),
            age_field: p.peek_age.map(|f| f.data().to_vec()),
            cratonic_factor: p.peek_cratonic_factor.map(|f| f.data().to_vec()),
            plate_id: None,
            plate_type: None,
            boundary_flag: None,
        }
    }

    /// Convert from the harness-side `FinalState` (which carries
    /// `Field2D` and v2 enums) to the bridge-side payload (flat
    /// `Vec`s + small ints).
    pub fn from_harness(state: &FinalState) -> Self {
        use ymir_core::tectonics_v2::boundaries::{BoundaryFlag, PlateType};

        let plate_id = state
            .plate_id
            .as_ref()
            .map(|f| f.data().to_vec());
        let plate_type = state.plate_type.as_ref().map(|f| {
            f.data()
                .iter()
                .map(|t| match t {
                    PlateType::Oceanic => 0u8,
                    PlateType::Continental => 1u8,
                })
                .collect::<Vec<_>>()
        });
        let boundary_flag = state.boundary_flag.as_ref().map(|f| {
            f.data()
                .iter()
                .map(|b| match b {
                    BoundaryFlag::None => 0u8,
                    BoundaryFlag::Subduction => 1u8,
                    BoundaryFlag::OceanicSubduction => 2u8,
                    BoundaryFlag::Rift => 3u8,
                    BoundaryFlag::ContinentalCollision => 4u8,
                })
                .collect::<Vec<_>>()
        });

        V2FinalState {
            nx: state.nx,
            ny: state.ny,
            dx: state.dx,
            dy: state.dy,
            s_field: state.s_field.data().to_vec(),
            vx: state.vx.clone(),
            vy: state.vy.clone(),
            strain_rate_invariant: state.strain_rate_invariant.data().to_vec(),
            age_field: state.age_field.as_ref().map(|f| f.data().to_vec()),
            cratonic_factor: state.cratonic_factor.as_ref().map(|f| f.data().to_vec()),
            plate_id,
            plate_type,
            boundary_flag,
        }
    }
}
