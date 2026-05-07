//! Step 12 — interleaved tectonic-erosion workflow orchestrator.
//!
//! Pulls together three pieces that already exist in core (the v2
//! tectonic harness, legacy `compute_isostasy`, the rain-drop
//! `run_erosion`, the bicubic + FBM `upscale_with_fbm`) and adds the
//! one missing piece (a parametric low-res erosion algorithm) into a
//! two-phase pipeline:
//!
//! - **Phase A** (low-res loop) — `tectonic k_cycle steps → isostasy
//!   (for sea_level_normalized) → low-res erosion → reclassify →
//!   recompute craton`, repeated `N_cycles` times. Continuation
//!   between cycles is cheap by virtue of Step 8.6's
//!   [`crate::tectonics_v2::diagnostics::harness::ContinuationState`]
//!   feeding the next cycle's velocity warm start (no transient
//!   pay-back, see D3 of `step12_issue.md`).
//! - **Phase B** (HD finalization) — `Field2D → GridF32 → upscale +
//!   FBM → run_erosion`, once at the end of Phase A. Produces the
//!   2048² heightmap consumed by the Living Landz workflow.
//!
//! Phase 1 (this commit) ships only the type scaffolding plus a
//! [`WorkflowConfig::Disabled`] passthrough. The Disabled variant is
//! the default, structurally short-circuits every workflow branch,
//! and is the **bit-identical contract** for Steps 0–13.5 regression:
//! calling [`phase_a::run_phase_a_cycle`] with `Disabled` is exactly
//! [`crate::tectonics_v2::diagnostics::harness::run_baseline`] under
//! the hood, no extra allocation, no extra RNG consumption, no path
//! deviation — see `tests/v2_workflow_disabled_regression.rs`.
//!
//! The `Enabled(_)` variant carries the Phase A loop parameters
//! (`N_cycles`, `k_cycle`, `α`, `β`) and Phase B HD parameters
//! (target grid size, FBM config, erosion config, grand-scale
//! tolerance). The orchestration logic for `Enabled` lands in
//! Phases 2–5; Phase 1 stubs it as `unimplemented!`.

pub mod low_res_erosion;
pub mod phase_a;
pub mod phase_b;

use serde::{Deserialize, Serialize};

use crate::erosion::hydraulic::ErosionConfig;
use crate::terrain::upscale::FbmUpscaleConfig;

pub use low_res_erosion::ErosionStats;
pub use phase_a::{run_phase_a_cycle, run_phase_a_loop};
pub use phase_b::run_phase_b;

/// Top-level workflow on/off switch.
///
/// `Disabled` is the default and short-circuits the entire orchestrator
/// (Phase A loop + Phase B HD); the existing Step 11 / Step 13.5 entry
/// points are reached via [`crate::tectonics_v2::diagnostics::harness::run_baseline`]
/// directly. `Enabled(params)` activates the multi-cycle loop and the
/// HD finalization with the supplied parameters.
///
/// Default: `Disabled` — preserves Step 11 standalone behaviour
/// bit-for-bit. The Step 12 acceptance #15 regression test verifies
/// this by calling [`run_phase_a_cycle`] with `Disabled` and comparing
/// final-state byte-by-byte against a direct `run_baseline` call.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub enum WorkflowConfig {
    #[default]
    Disabled,
    Enabled(WorkflowParams),
}

/// Phase A + Phase B parameter bundle, only meaningful when wrapped in
/// [`WorkflowConfig::Enabled`]. Defaults are conservative starting
/// points (D8 of `step12_issue.md`), not calibrated values — Phase 8
/// reports document the empirically-tuned defaults future users should
/// prefer.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WorkflowParams {
    pub phase_a: PhaseAParams,
    pub phase_b: PhaseBParams,
}

impl Default for WorkflowParams {
    fn default() -> Self {
        Self {
            phase_a: PhaseAParams::default(),
            phase_b: PhaseBParams::default(),
        }
    }
}

/// Phase A loop parameters (D2 + D8).
///
/// Defaults: `N_cycles = 5`, `k_cycle = 20` → 100 effective tectonic
/// steps + 5 erosion passes — a starting point for exploration, not a
/// calibrated configuration.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhaseAParams {
    /// Number of cycles in the Phase A loop. Default: 5 (D8).
    pub n_cycles: usize,
    /// Tectonic steps per cycle. Default: 20 (D8).
    pub k_cycle: usize,
    /// Erosion rate `α` per cycle (cycle-time, not absolute). Default:
    /// 0.01. Range typically `[0.001, 0.05]`.
    pub alpha: f64,
    /// Sediment redistribution coefficient `β`. `0.0` = pure erosion
    /// (mass leaves the grid); `1.0` = full deposition downslope.
    /// Default: 0.0 (D8).
    pub beta: f64,
}

impl Default for PhaseAParams {
    fn default() -> Self {
        Self { n_cycles: 5, k_cycle: 20, alpha: 0.01, beta: 0.0 }
    }
}

/// Phase B HD finalization parameters (D5 + D8).
///
/// `hd_grid_size` overrides `fbm.target_size` for the upscale stage
/// (typically 2048 for Living Landz). The `grand_scale_tolerance` is
/// the D5 acceptance threshold — Phase 5 measures the actual deviation
/// without modifying `run_erosion`; if it lands above the tolerance,
/// remontée before any post-correction.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhaseBParams {
    /// Target HD grid resolution (longer axis). Default: 2048.
    pub hd_grid_size: usize,
    /// FBM upscale configuration. `target_size` is overridden by
    /// `hd_grid_size` at run time, so the field here may be left
    /// at default.
    #[serde(default)]
    pub fbm: FbmUpscaleConfig,
    /// HD rain-drop erosion configuration. `sea_level` is overridden
    /// at run time from the upstream isostasy normalisation.
    #[serde(default)]
    pub erosion: ErosionConfig,
    /// D5 grand-scale preservation tolerance:
    /// `‖S̃_HD_after - upscale(S̃_low_res)‖_∞ < tol`. Default: 0.10.
    pub grand_scale_tolerance: f64,
}

impl Default for PhaseBParams {
    fn default() -> Self {
        Self {
            hd_grid_size: 2048,
            fbm: FbmUpscaleConfig::default(),
            erosion: ErosionConfig::default(),
            grand_scale_tolerance: 0.10,
        }
    }
}

/// Output of a single Phase A cycle.
///
/// `baseline` carries the full
/// [`crate::tectonics_v2::diagnostics::harness::BaselineResult`] from
/// the cycle's tectonic run (final state, metrics, config dump);
/// `erosion_volume_removed` is the integrated `Δh` summed over
/// continental cells during the low-res erosion pass (`0.0` when
/// `WorkflowConfig::Disabled`).
#[derive(Debug)]
pub struct CycleOutput {
    pub baseline: crate::tectonics_v2::diagnostics::harness::BaselineResult,
    pub erosion_volume_removed: f64,
}

/// Output of the full Phase A loop (one [`CycleOutput`] per cycle).
///
/// When `WorkflowConfig::Disabled`, the vec contains exactly one
/// entry — the direct `run_baseline` passthrough — so the regression
/// contract holds: `output.cycles[0].baseline ≡ run_baseline(cfg)`.
#[derive(Debug)]
pub struct PhaseAOutput {
    pub cycles: Vec<CycleOutput>,
}

/// Output of Phase B HD finalization.
///
/// `grand_scale_deviation` is the D5 metric `‖HD_after -
/// upscale(low_res)‖_∞`, recorded for both the acceptance test
/// (Phase 5) and the physics report (Phase 8).
#[derive(Debug)]
pub struct PhaseBOutput {
    pub heightmap: crate::grid::GridF32,
    pub sediment: crate::grid::GridF32,
    pub slope: crate::grid::GridF32,
    pub grand_scale_deviation: f64,
}
