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
pub mod heightmap;
pub mod metrics;
pub mod newton_metrics;

// Issue #117 — `report` aggregates renders from every sweep module
// (ar/bi/br/k_sub/num_plates) plus `mms_bench` so it transitively
// retires with them. The bins consuming `write_markdown_report`
// (step_baseline, step5-8_baseline) already gate via Cargo
// `required-features = ["v2_legacy"]`.
#[cfg(feature = "v2_legacy")]
pub mod report;

// Issue #117 — these diagnostic drivers are Stokes-coupled (they import
// from `_attic/{stokes, mantle, slab, rheology, basal_drag}`) and gate
// behind `v2_legacy`. The framework types they expose (`BaselineConfig`,
// `BaselineResult`, `run_baseline`) follow the same gate so callers
// (workflow/phase_a.rs, bridge/v2/*, tests) compile coherently.
#[cfg(feature = "v2_legacy")]
pub mod ar_sweep;
#[cfg(feature = "v2_legacy")]
pub mod bi_sweep;
#[cfg(feature = "v2_legacy")]
pub mod br_sweep;
#[cfg(feature = "v2_legacy")]
pub mod harness;
#[cfg(feature = "v2_legacy")]
pub mod k_sub_sweep;
#[cfg(feature = "v2_legacy")]
pub mod mms_bench;
#[cfg(feature = "v2_legacy")]
pub mod num_plates_sweep;

#[cfg(feature = "v2_legacy")]
pub use harness::{BaselineConfig, BaselineResult, run_baseline};

pub use metrics::{IterationHistogram, Metrics, SolverConfigDump};
pub use newton_metrics::NewtonAggregate;

#[cfg(feature = "v2_legacy")]
pub use report::write_markdown_report;
