//! Worker-thread command channel for `bridge::c1`.
//!
//! ## MVP cancel semantics (Issue #137 Q-E1.3 Option C)
//!
//! `Cancel` sets the shared `AtomicBool` cancel flag but does NOT
//! interrupt an in-flight `RunBaseline`. The C1 time loop
//! (`run_with_closures`) has no cancel hook — adding one would
//! require reopening `ymir-core` (Stage S anti-pattern: NE PAS
//! rouvrir core pour cancel token).
//!
//! C1 runs at the Viz-0 default (`grid_size = 64, n_steps = 300,
//! Track D enabled`) measure ~250 ms (Issue #132 Stage A). A
//! user-initiated `Cancel` will land between commands — the
//! current run completes naturally, then the cancel takes effect
//! on the NEXT `RunBaseline` (worker checks the flag and skips
//! the run). This is acceptable for MVP given the sub-second run
//! duration.
//!
//! Viz-0-bis candidate: add a thread-local cancel token to
//! `ymir-core::tectonics_c1::cancel` (mirroring
//! `tectonics_v2::cancel`), threaded through `run_with_closures`'s
//! per-step boundary. ~1 day effort.

use super::spec::C1RunSpec;

#[derive(Clone, Debug)]
pub enum C1Command {
    /// Launch a baseline run with `spec`. Worker emits
    /// `Started → StepCompleted × n_steps → Completed` (or
    /// `Failed` on panic, currently unreachable per Q-E1.2).
    RunBaseline { spec: C1RunSpec },
    /// Set the shared cancel flag. Does NOT interrupt the
    /// current run — see module docstring.
    Cancel,
}
