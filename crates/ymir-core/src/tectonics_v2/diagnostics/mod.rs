//! Diagnostics framework for the solver-reconstruction milestone.
//!
//! Step 0 ships the MVP: a [`Metrics`] struct with slots for every
//! metric the full milestone will eventually populate, a markdown
//! writer, and a baseline [`harness`] that runs the coupled Stokes +
//! advection loop and emits `step0_report.md`.
//!
//! Metrics that correspond to physics not yet active (yield fraction,
//! `S̃_eq`, boundary-type diversity, Newton outcomes, cratonic
//! stability, age field statistics) are `Option<...>` and left
//! `None` — the framework slot exists so that subsequent steps can
//! populate it without touching the report layout.

pub mod comparison;
pub mod harness;
pub mod metrics;
pub mod mms_bench;
pub mod newton_metrics;
pub mod report;

pub use harness::{BaselineConfig, BaselineResult, run_baseline};
pub use metrics::{IterationHistogram, Metrics, SolverConfigDump};
pub use newton_metrics::NewtonAggregate;
pub use report::write_markdown_report;
