//! Bevy-side bridge for the C1 lightweight dynamic tectonics
//! solver (Issue #137 Viz-0 C1 integration).
//!
//! Mirrors `bridge::v2` shape but with smaller scope:
//! single-baseline `RunBaseline`, MVP cancel-between-runs only
//! (Q-E1.3 Option C), per-step `StepCompleted { snapshot }`
//! streaming, raw-fields snapshot for view-switch during pause.
//!
//! The Bevy plugin registration lives in `plugin.rs` (added
//! Stage E5). Stage E2 ships only the worker-thread surface
//! (this module's children).

pub mod commands;
pub mod events;
pub mod plugin;
pub mod snapshot;
pub mod spec;
pub mod thread;

pub use commands::C1Command;
pub use events::C1Event;
pub use plugin::{C1BridgePlugin, C1RunState, C1SolverBridge};
pub use snapshot::C1Snapshot;
pub use spec::C1RunSpec;
pub use thread::spawn_c1_thread;
