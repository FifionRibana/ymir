//! `tectonics_v2`-bound bridge — the only solver bridge after the
//! Step 8.6 Phase 8h sunset. A background thread owns the v2
//! simulation harness and processes `V2Command`s; results stream
//! back as `V2Event`s via a crossbeam channel. The bridge knows
//! nothing about Bevy; the Bevy plugin wrapper lives in
//! `bridge::v2::plugin`.
//!
//! `pub use` re-exports below give integration tests under
//! `crates/ymir-viz/tests/` a flat import path
//! (`ymir_viz::bridge_v2::*`). The bin does not reach all of them,
//! hence the `allow(unused_imports)` — these are public API
//! surface, not dead code.

#![allow(unused_imports)]

pub mod build_config;
pub mod commands;
pub mod events;
pub mod plugin;
pub mod presets;
pub mod snapshot;
pub mod spec;
pub mod thread;

pub use commands::V2Command;
pub use events::{V2Event, V2FinalState};
pub use plugin::{V2BridgePlugin, V2RunState, V2SolverBridge};
pub use snapshot::{V2RunSnapshot, V2ScalarMetrics, SNAPSHOT_FORMAT_VERSION};
pub use spec::{
    V2AgeFieldSpec, V2CratonicSpec, V2ForceKind, V2InitModeSpec, V2LinearSolverSpec, V2MantleSpec,
    V2RunSpec,
};
pub use thread::spawn_v2_thread;
