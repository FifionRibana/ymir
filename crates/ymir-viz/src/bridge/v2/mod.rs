//! Step 8.6 — `tectonics_v2`-bound bridge.
//!
//! Mirrors `bridge::legacy` / the existing `bridge` module pattern: a
//! background thread owns the v2 simulation harness and processes
//! `V2Command`s; results stream back as `V2Event`s via a crossbeam
//! channel. The bridge knows nothing about Bevy; the Bevy plugin
//! wrapper lives in `bridge::v2::plugin` (Phase 2).
//!
//! Phase 1 scope: spawn the thread, run a single baseline to
//! completion, return the final state. No live progress streaming
//! (the harness `run_baseline` is non-cancellable end-to-end at
//! Step 8.6 Phase 1 — refactoring it to take a step callback is
//! Phase 5 work, deferred until the visualization actually needs
//! intermediate snapshots).

pub mod build_config;
pub mod commands;
pub mod events;
pub mod plugin;
pub mod presets;
pub mod spec;
pub mod thread;

pub use commands::V2Command;
pub use events::{V2Event, V2FinalState};
pub use plugin::{V2BridgePlugin, V2RunState, V2SolverBridge};
pub use spec::{
    V2AgeFieldSpec, V2CratonicSpec, V2ForceKind, V2InitModeSpec, V2LinearSolverSpec, V2MantleSpec,
    V2RunSpec,
};
pub use thread::spawn_v2_thread;
