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

use ymir_core::tectonics_v2::workflow::PhaseAParams;

use super::hd::HdParams;
use super::spec::C1RunSpec;

#[derive(Clone, Debug)]
pub enum C1Command {
    /// Launch a **gallery-path** baseline run with `spec`: standalone
    /// `run_with_closures`, NO `apply_post_tectonic`. Worker emits
    /// `Started → StepCompleted × n_steps → Completed` (or `Failed`
    /// on panic, currently unreachable per Q-E1.2). Reproduces the
    /// Track A/B/D visual gallery (Issue #137 gallery contract).
    RunBaseline { spec: C1RunSpec },
    /// Launch a **workflow-path** run (Issue #139 Stage E2): the
    /// calibrated Phase A loop — `n_cycles` cycles of
    /// `run_with_closures(k_cycle)` each followed by
    /// `apply_post_tectonic` (sea-level → macro-redistribution →
    /// reclassify). `phase_a` defaults to `PhaseAParams::default()`
    /// (n_cycles 5, k_cycle 20, alpha 0.01) — the calibrated cadence.
    /// Total tectonic steps = `n_cycles × k_cycle`, NEVER `n_steps`
    /// (the A1-c over-erosion guard). The gallery `RunBaseline` path
    /// is left untouched (W4).
    RunWorkflow {
        spec: C1RunSpec,
        phase_a: PhaseAParams,
    },
    /// Launch the HD production chain (UI rewrite step b/5): tectonics →
    /// upscale → erosion → bathymetry → drainage → climate → biomes, via
    /// the cached `ymir-core` production functions, on the worker thread.
    /// Emits `HdStarted → (HdPhaseStarted → HdPhaseDone) × 4 → HdCompleted`
    /// (or `HdFailed`). Cancellable between phases.
    RunHd { spec: C1RunSpec, params: HdParams },
    /// Set the shared cancel flag. Does NOT interrupt the
    /// current run — see module docstring.
    Cancel,
}
