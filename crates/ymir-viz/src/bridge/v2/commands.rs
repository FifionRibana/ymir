//! Step 8.6 — commands sent from the Bevy main thread to the v2
//! solver thread.

use super::spec::V2RunSpec;

#[allow(clippy::large_enum_variant)]
pub enum V2Command {
    /// Run a full baseline simulation per the supplied spec. The
    /// thread blocks for the duration of the run and emits a final
    /// `V2Event::Completed` (or `V2Event::Failed`) when done.
    RunBaseline { spec: V2RunSpec },
    /// Cooperative cancellation request. Phase 1 is a no-op (the
    /// harness `run_baseline` does not currently expose a
    /// cancellation hook); Phase 5 will refactor `run_baseline` to
    /// accept a step callback and honour this signal.
    Cancel,
    /// Terminate the worker thread cleanly. Sent by the Bevy
    /// plugin's `Drop` impl so the test harness doesn't leak threads
    /// across fixtures.
    Shutdown,
}
