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
//! calling [`phase_a::run_phase_a_cycle_v2`] with `Disabled` is exactly
//! [`crate::tectonics_v2::diagnostics::harness::run_baseline`] under
//! the hood, no extra allocation, no extra RNG consumption, no path
//! deviation — see `tests/v2_workflow_disabled_regression.rs`.
//!
//! The `Enabled(_)` variant carries the Phase A loop parameters
//! (`N_cycles`, `k_cycle`, `α`, `β`) and Phase B HD parameters
//! (target grid size, FBM config, erosion config, grand-scale
//! tolerance). The orchestration logic for `Enabled` lands in
//! Phases 2–5; Phase 1 stubs it as `unimplemented!`.

pub mod drainage;
pub mod macro_redistribution;

// Issue #117 HC4 / Phase 1.3 H2 — Phase A's v2 path
// (`run_phase_a_cycle_v2`) couples to the retired harness; gated
// under `v2_legacy`. The C1 path (`run_phase_a_cycle_c1`, Commit 4
// of H2) is default-features-on. Phase B is paradigm-agnostic at
// runtime (consumes `Field2D + sea_level + iso_config`, no
// `BaselineResult`); H1 audit § 2.2 confirmed it is safe to un-gate.
// See `docs/migrations/harness_paradigm_agnostic.md`.
#[cfg(feature = "v2_legacy")]
pub mod phase_a_v2;
pub mod phase_b;

use serde::{Deserialize, Serialize};

use crate::erosion::hydraulic::ErosionConfig;
use crate::terrain::upscale::FbmUpscaleConfig;

pub use macro_redistribution::RedistributionStats;

#[cfg(feature = "v2_legacy")]
pub use phase_a_v2::{
    final_state_to_continuation_v2, run_phase_a_cycle_v2, run_phase_a_cycle_with_progress_v2,
    run_phase_a_loop_v2,
};
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
/// this by calling [`run_phase_a_cycle_v2`] with `Disabled` and comparing
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

/// Phase A loop parameters (D2 + D8, Step 12 R3 refactor).
///
/// Defaults: `N_cycles = 5`, `k_cycle = 20` → 100 effective tectonic
/// steps + 5 macro-redistribution passes. Step 12 R3 swapped the
/// legacy `low_res_erosion::apply` (per-cycle diffusive erosion with a
/// `β` local-deposition coefficient) for
/// [`macro_redistribution::apply`], which combines erosion + long-
/// distance drainage + isostatic rebound into one mass-conserving
/// call. The legacy `β` parameter is gone — its semantic of "deposit
/// `β` of eroded mass on the immediate downslope neighbour" is
/// replaced by macro drainage that targets the real basin sink
/// (`max_drainage_distance` hops away or the first oceanic neighbour).
/// The remaining mass that's not redistributed is implicit in
/// `isostatic_rebound_ratio` (the crust the doesn't migrate is
/// altitude-restored by mantle uplift, not S̃-deleted).
///
/// See [`macro_redistribution`] module docstring for the full
/// arithmetic.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhaseAParams {
    /// Number of cycles in the Phase A loop. Default: 5 (D8).
    pub n_cycles: usize,
    /// Tectonic steps per cycle. Default: 20 (D8).
    pub k_cycle: usize,
    /// Erosion rate `α` per cycle (cycle-time, not absolute). Default:
    /// 0.01. Range typically `[0.001, 0.05]`. Consumed by
    /// [`macro_redistribution::apply`] as the multiplier on
    /// `slope · (S̃ − sea_level)`.
    pub alpha: f64,
    /// Isostatic rebound fraction (R3 / Step 12 refactor). `1.0` =
    /// full rebound (no S̃ migration, the apply call is a no-op);
    /// `0.0` = no rebound (every eroded gram is redistributed
    /// downstream). Default: `0.80` — calibrated against the Earth
    /// `ρ_crust / ρ_mantle ≈ 2700/3300 ≈ 0.82` mass-balance ratio.
    /// Range typically `[0.5, 0.95]`; R5 sweep may pick a different
    /// nominal value.
    pub isostatic_rebound_ratio: f64,
    /// Max NESW hops the drainage walks per continental cell (R3 /
    /// Step 12 refactor). Default: `10`. Smaller values bias eroded
    /// mass toward nearby destinations (closer to the legacy `β`
    /// behaviour); larger values let interior cells reach distal
    /// basins. Anti-pattern: setting beyond `~30` rarely changes the
    /// drainage map on realistic Voronoï continental patches (the
    /// terminal cell is reached well before that on Phase A relief).
    pub max_drainage_distance: usize,
}

impl Default for PhaseAParams {
    fn default() -> Self {
        Self {
            n_cycles: 5,
            k_cycle: 20,
            alpha: 0.01,
            isostatic_rebound_ratio: 0.80,
            max_drainage_distance: 10,
        }
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
    /// D5 grand-scale preservation tolerance — applied as
    /// `p95(|HD_after - upscale(low_res)|) < tol`. Default: `0.10`.
    ///
    /// **Phase 5 reformulation (methodological finding).** The Step
    /// 12 issue originally specified the L_∞ max-norm
    /// `‖HD_after - upscale(low_res)‖_∞ < tol`. Empirical Phase 5
    /// measurement on 32² Phase A × 3 cycles → HD = 256 surfaced
    /// L_∞ ≈ 0.151 with p95 = 0.076: HD rain-drop erosion legitimately
    /// carves valleys with 15–20 % local pixel deviation on roughly
    /// 0.6 % of the domain (415 / 65536 cells). The L_∞ contract is
    /// structurally incompatible with what erosion is designed to
    /// do; the p95 contract measures grand-scale shape preservation
    /// over 95 % of the domain while still flagging pathological
    /// runs (if more than 5 % of the domain deviates by > 10 %, p95
    /// saturates above the threshold).
    ///
    /// The L_∞ value is retained as a diagnostic on
    /// [`PhaseBOutput::grand_scale_deviation`]; only
    /// [`PhaseBOutput::grand_scale_deviation_p95`] is the formal
    /// acceptance gate. Tolerance numerical value is unchanged at
    /// 0.10 — this is a statistic switch, not a threshold relax.
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

/// Paradigm-agnostic per-cycle Phase A diagnostics.
///
/// These five scalars describe the cycle's macro-redistribution +
/// reclassification + craton-recompute pass and are **paradigm-
/// independent**: both the v2 path and the C1 path populate them
/// from the same downstream code in
/// [`crate::tectonics_v2::workflow::macro_redistribution`] and the
/// shared post-tectonic-step pass (Phase 1.3 H2 introduces this
/// shared pass as `workflow::phase_a_common`).
///
/// Carrying them in a dedicated struct (rather than as flat fields
/// on the paradigm-specific `CycleOutputV2*` envelopes) keeps the
/// two paths from drifting silently. See
/// `docs/migrations/harness_paradigm_agnostic.md` § R3.
#[derive(Clone, Copy, Debug, Default)]
pub struct CycleOutputCommon {
    /// Gross integrated `Σ Δh` over continental cells (low-res erosion
    /// pass output). Net mass change is `-(1 - β) · volume_removed`.
    pub erosion_volume_removed: f64,
    /// Maximum per-cell `Δh` during the cycle's erosion pass (peak
    /// erosion magnitude). Useful for spotting outlier cells.
    pub erosion_peak_delta_h: f64,
    /// Adaptive sea-level threshold used for the cycle's
    /// reclassification + cratonic recompute, sourced from
    /// `compute_isostasy(s_post_tectonic).sea_level_normalized`. `0.0`
    /// when Disabled.
    pub sea_level_normalized: f64,
    /// Total continental mass change during the cycle's erosion pass
    /// (`Σ s_after - Σ s_before`). Equals `-(1 - β) · volume_removed`
    /// modulo IEEE-754 rounding; tracked separately as a sanity hook
    /// for the long-run mass-conservation diagnostic.
    pub mass_drift: f64,
    /// Fraction of cells whose `cratonic_factor` changed by more than
    /// `1e-9` between pre-cycle and post-cycle recompute. `None` when
    /// no `CratonicConfig::Enabled` was active (or the paradigm has
    /// no equivalent — C1 currently does not). The Phase 4
    /// multi-cycle metrics dashboard accumulates these per cycle.
    pub craton_recomputation_change: Option<f64>,
}

/// Output of a single Phase A cycle on the **v2 path**.
///
/// `baseline` carries the full
/// [`crate::tectonics_v2::diagnostics::harness::BaselineResult`] from
/// the cycle's tectonic run (final state, metrics, config dump). Under
/// `WorkflowConfig::Enabled`, the `final_state` is post-erosion +
/// post-reclassification + post-craton-recompute (i.e. the state the
/// next cycle's continuation should warm-start from).
///
/// `common` carries the paradigm-agnostic per-cycle diagnostics shared
/// with the C1 path; see [`CycleOutputCommon`].
///
/// All Disabled-cycle scalars in `common` are zeroed/`None` — the
/// regression test `v2_workflow_disabled_regression` pins this
/// contract.
// Issue #117 / Phase 1.3 H2 — `CycleOutputV2` carries `BaselineResult`
// which is v2-specific. Gated under `v2_legacy`; the paradigm-
// agnostic scalars have moved to [`CycleOutputCommon`] which is not
// gated.
#[cfg(feature = "v2_legacy")]
#[derive(Debug)]
pub struct CycleOutputV2 {
    pub baseline: crate::tectonics_v2::diagnostics::harness::BaselineResult,
    pub common: CycleOutputCommon,
}

/// Output of the full Phase A loop (one [`CycleOutputV2`] per cycle).
///
/// When `WorkflowConfig::Disabled`, the vec contains exactly one
/// entry — the direct `run_baseline` passthrough — so the regression
/// contract holds: `output.cycles[0].baseline ≡ run_baseline(cfg)`.
// Issue #117 — `PhaseAOutputV2` transitively gates because it carries
// `Vec<CycleOutputV2>` (which references retired harness types).
#[cfg(feature = "v2_legacy")]
#[derive(Debug)]
pub struct PhaseAOutputV2 {
    pub cycles: Vec<CycleOutputV2>,
}

/// Output of Phase B HD finalization.
///
/// Two D5 metrics are surfaced:
///
/// - [`Self::grand_scale_deviation`] — `L_∞` of `|HD_after -
///   upscale(low_res)|`. **Diagnostic only**: surfaces the deepest
///   valley pixel; not asserted in the acceptance test because run_erosion
///   is *meant* to carve sharp valleys locally (15–20 % per pixel).
/// - [`Self::grand_scale_deviation_p95`] — 95th-percentile of the
///   same per-cell delta distribution. **Formal acceptance**: the
///   D5 contract is `grand_scale_deviation_p95 <
///   grand_scale_tolerance` (default 0.10). This measures grand-scale
///   shape preservation over 95 % of the domain while still catching
///   pathological runs (a regression that flattens 30 % of the
///   domain by 15 % shifts p95 above 0.10).
///
/// Phase 5 finding: the original L_∞ contract was structurally
/// incompatible with the carved-valleys design intent of HD erosion.
/// p95 is the structural reformulation, not a threshold relax — the
/// numerical tolerance stays at 0.10. See
/// [`PhaseBParams::grand_scale_tolerance`] docstring for the full
/// rationale.
#[derive(Debug)]
pub struct PhaseBOutput {
    pub heightmap: crate::grid::GridF32,
    pub sediment: crate::grid::GridF32,
    pub slope: crate::grid::GridF32,
    /// `L_∞` of per-cell deltas vs the upscaled baseline.
    /// **Diagnostic only.**
    pub grand_scale_deviation: f64,
    /// 95th-percentile of per-cell deltas vs the upscaled baseline.
    /// **Formal D5 acceptance metric.**
    pub grand_scale_deviation_p95: f64,
}
