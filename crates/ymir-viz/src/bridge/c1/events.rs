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

use super::snapshot::C1Snapshot;
use super::spec::C1RunSpec;

#[derive(Clone, Debug)]
pub enum C1Event {
    /// Worker has accepted a `RunBaseline` command and is about
    /// to start the C1 time loop. Caller can reset UI state on
    /// receipt.
    Started { spec: C1RunSpec },
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
