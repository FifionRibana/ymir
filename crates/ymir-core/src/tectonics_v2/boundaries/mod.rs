//! Boundary sources/sinks — Step 5.
//!
//! Introduces the **first spatially heterogeneous cell typing** of
//! the milestone. Every cell carries a [`PlateType`] (Oceanic /
//! Continental) and a [`BoundaryFlag`] (None, Subduction, Rift,
//! ContinentalCollision, OceanicSubduction). Both fields are
//! **prescribed statically** at Step 5: see
//! [`layouts::build_layout`] for the synthetic generators that
//! produce them. Dynamic detection from the velocity field is
//! Step 6.
//!
//! The machinery provides:
//! - Five source/sink terms (`Q_sub`, `Q_arc`, `Q_spread`,
//!   `Q_coll-v`, `Q_rift-v`) via [`source_sink::compute_source_sink_terms`].
//! - A hard floor `S̃ ≥ S_MIN` with artificial-flux tracking for the
//!   mass-balance residual (see [`clamp::apply_clamp_with_tracking`]
//!   and issue #89 D5/D8).
//! - A bisection calibration for `k_spread` targeting
//!   `s_oceanic_mean ∈ [0.18, 0.22]` at steady state (issue #89 D4),
//!   in [`calibration::calibrate_k_spread`].
//!
//! When `BoundaryConfig::Disabled` is passed to the harness, the
//! whole pipeline is **structurally bypassed**: no Q evaluation, no
//! clamp, no tracking. This is the zero-cost-when-disabled invariant
//! that the Step 5 regression run verifies against the Step 4
//! physics baseline with yielding Enabled.

pub mod boundary_flag;
pub mod calibration;
pub mod clamp;
pub mod closed_mode;
pub mod crust_geometry;
pub mod layouts;
pub mod plate_type;
pub mod source_sink;
pub mod stats;

pub use boundary_flag::{
    BoundaryConfig, BoundaryFlag, BoundaryFlagField, BoundaryFlagFieldShared, BoundaryRates,
    PlateTypeFieldShared, RecyclingModeInit,
};
pub use calibration::{
    CalibrationError, CalibrationResult, K_SPREAD_BRACKET, KSpreadCalibration, calibrate_k_spread,
};
pub use clamp::{ClampStats, S_MIN, apply_clamp_with_tracking};
pub use closed_mode::{
    compute_q_sub_only, count_immediate_eligibilities, count_spread_eligibility,
    distribute_delayed, distribute_immediate, integrate_sub_mass,
};
pub use crust_geometry::{CrustGeometry, CrustGeometryShared, GeometryKind};
pub use layouts::{
    BoundaryLayout, balanced_sub_spread, build_layout, continental_collision_band,
    horizontal_oceanic_strip, vertical_rift_line,
};
pub use plate_type::{PlateType, PlateTypeField};
pub use source_sink::{compute_source_sink_terms, convergent_component, div_v_cell};
pub use stats::{
    BoundaryMechanismActive, MeanStd, boundary_type_diversity, interface_mask,
    s_continental_collision_mean, s_continental_interior, s_oceanic,
};
