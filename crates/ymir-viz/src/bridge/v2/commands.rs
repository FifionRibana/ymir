//! Step 8.6 — commands sent from the Bevy main thread to the v2
//! solver thread.

use super::events::V2FinalState;
use super::spec::V2RunSpec;

#[allow(clippy::large_enum_variant)]
pub enum V2Command {
    /// Run a full baseline simulation per the supplied spec. The
    /// thread blocks for the duration of the run and emits a final
    /// `V2Event::Completed` (or `V2Event::Failed`) when done.
    RunBaseline { spec: V2RunSpec },
    /// Step 8.6 follow-up — continue from a prior run's final state.
    /// `spec` is the user-edited spec (may differ from the source
    /// run's spec, but voronoi-relevant fields — seed, num_plates,
    /// continental_ratio, grid dims — should match for the
    /// continuation to make physical sense). `from_state` carries
    /// the S̃ / vx / vy / age / cratonic_factor rasters; the bridge
    /// thread wraps them into a `harness::ContinuationState` and
    /// the harness uses them as the run's start state instead of
    /// computing init.
    ContinueRun { spec: V2RunSpec, from_state: V2FinalState },
    /// Pre-Phase-8h-follow-up cancellation command. Now that
    /// `V2SolverBridge::request_cancel` flips the shared
    /// `AtomicBool` directly, this variant is mainly kept for
    /// channel-level shutdown / tests; sending it has the same
    /// effect as setting the flag.
    Cancel,
    /// Terminate the worker thread cleanly. Sent by the Bevy
    /// plugin's `Drop` impl so the test harness doesn't leak threads
    /// across fixtures.
    Shutdown,
    /// Step 12 Phase 7 — run the Phase A multi-cycle loop end-to-end.
    /// `spec.workflow` must be `V2WorkflowSpec::On { .. }`; otherwise
    /// the thread emits `V2Event::Failed` (the orchestrator's
    /// `Disabled` branch is short-circuited at single-cycle, this
    /// command is the multi-cycle entry point). Streams
    /// `V2Event::Progress` events at the inner-cycle step granularity
    /// (one per harness step), `V2Event::WorkflowCycleCompleted`
    /// between cycles, and a final `V2Event::WorkflowPhaseACompleted`.
    /// Cancellation lands at the next harness step boundary.
    RunWorkflowPhaseA { spec: V2RunSpec },
    /// Step 12 Phase 7 — continue a Phase A multi-cycle loop from a
    /// prior run's final state. Cycle 1 receives `from_state` as the
    /// `ContinuationState` warm-start; subsequent cycles use the
    /// orchestrator's natural continuation chain. Other than that,
    /// behaviour matches `RunWorkflowPhaseA`.
    ContinueWorkflowPhaseA {
        spec: V2RunSpec,
        from_state: V2FinalState,
    },
    /// Step 12 Phase 7 — run Phase B HD finalization on a Phase A
    /// output (single shot, ~30-90 s at 2048² × 5M droplets). Emits
    /// a single `V2Event::WorkflowPhaseBCompleted` with the HD
    /// heightmap + sediment + grand-scale metrics. The bridge
    /// thread reuses `from_state.s_field` directly (the legacy
    /// `ContinuationState` shape carries every other field too,
    /// which Phase B does not consume).
    RunWorkflowPhaseB {
        spec: V2RunSpec,
        from_state: V2FinalState,
    },
}
