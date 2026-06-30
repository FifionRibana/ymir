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
pub mod hd;
pub mod plugin;
pub mod snapshot;
pub mod spec;
pub mod thread;

pub use events::C1RunKind;
pub use hd::{CacheRegime, HdParams, HdPhase, HdResult, HdState};
pub use plugin::{C1BridgePlugin, C1CumulativeStats, C1RunState, C1SolverBridge};
pub use spec::C1RunSpec;
