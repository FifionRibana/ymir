//! Worker-thread event channel for `bridge::c1`.
//!
//! ## Event sequence per `RunBaseline`
//!
//! ```text
//!   Started { spec }
//!   StepCompleted { snapshot }     ← n_steps times
//!   Completed { spec, final_snapshot, elapsed }
//! ```
//!
//! Or, on panic / error (currently unreachable per Issue #137
//! Q-E1.2 — `run_with_closures` has no fallible paths):
//!
//! ```text
//!   Started { spec }
//!   Failed { error }
//! ```
//!
//! ## Backpressure (Q1, bounded events channel)
//!
//! The plugin spawns the events channel with `bounded::<C1Event>(2)`
//! (vs v2's `bounded(256)` per Issue #137 Stage S design). When the
//! UI thread lags (slow render system), the `tx.send` call inside
//! the per-step callback blocks until the channel drains — effectively
//! pausing the simulation. This is the intended "tight backpressure
//! = pause semantics" behaviour (W1 of Stage E2).

use std::time::Duration;

use ymir_core::tectonics_v2::workflow::PhaseAParams;

use super::snapshot::C1Snapshot;
use super::spec::C1RunSpec;

/// Which pipeline a run is exercising (Issue #139 Stage E2). Carried
/// on `Started` so the poll system can resolve the step-counter total
/// and the UI can show the per-cycle counter in workflow mode.
#[derive(Clone, Debug)]
pub enum C1RunKind {
    /// Gallery path — standalone `run_with_closures`, total tectonic
    /// steps = `spec.n_steps` (Issue #137 contract).
    Gallery,
    /// Workflow path — calibrated Phase A loop, total tectonic steps
    /// = `phase_a.n_cycles × phase_a.k_cycle` (NOT `n_steps`).
    Workflow { phase_a: PhaseAParams },
}

impl C1RunKind {
    /// Total **tectonic** steps for the run — the `/N` denominator in
    /// the UI step counter. Workflow's per-cycle `apply_post_tectonic`
    /// snapshots are extra and do NOT advance this count.
    pub fn total_tectonic_steps(&self, spec: &C1RunSpec) -> usize {
        match self {
            C1RunKind::Gallery => spec.n_steps,
            C1RunKind::Workflow { phase_a } => phase_a.n_cycles * phase_a.k_cycle,
        }
    }
}

#[derive(Clone, Debug)]
pub enum C1Event {
    /// Worker has accepted a run command and is about to start.
    /// `kind` distinguishes the gallery vs workflow pipeline and
    /// carries the workflow cadence for the UI counter. Caller can
    /// reset UI state on receipt.
    Started { spec: C1RunSpec, kind: C1RunKind },
    /// One C1 step has just completed; `snapshot` carries the
    /// post-step state + Track D stats (Issue #137 Viz-D0
    /// Option B).
    StepCompleted { snapshot: C1Snapshot },
    /// The run has completed. `final_snapshot` is the same as
    /// the last `StepCompleted` snapshot at `step = n_steps - 1`,
    /// re-emitted for convenience (callers can latch on
    /// `Completed` for final-state UI without tracking the last
    /// `StepCompleted`).
    Completed {
        spec: C1RunSpec,
        final_snapshot: C1Snapshot,
        elapsed: Duration,
    },
    /// Unreachable at MVP (Q-E1.2). Reserved for future panic
    /// catch or NaN detection paths.
    Failed { error: String },
}
