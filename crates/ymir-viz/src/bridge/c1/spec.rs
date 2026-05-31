//! `C1RunSpec` — the user-facing configuration for a C1 baseline
//! run dispatched through the `bridge::c1` worker thread.
//!
//! Minimal surface per Issue #137 Q-S.4 (6 fields):
//!
//! - `grid_size` — square grid `N × N`. Default 64.
//! - `seed` — RNG seed for `init_c1_state_phase_2_r7`. Default 42.
//! - `n_steps` — number of forward-Euler steps. Default 300 (forwarded
//!   to `C1TimeLoopConfig.n_steps` by the worker; single source of
//!   truth lives here).
//! - `init_params` — `Phase2InitParams` for the R7 boundary
//!   displacement + clustering + ridge-aligned age (Track B
//!   foundations). Default = `Phase2InitParams::default()`.
//! - `closures` — `C1Closures` with all 7 closures including Track D
//!   trio enabled by default.
//! - `drainage_max_distance` — BFS depth cap for the per-step
//!   drainage classifier (Phase 1.4 erosion path). Default 30.
//!
//! Implicit (NOT configurable in Viz-0):
//!
//! - `iso_config = IsostasyConfig::default()` (the gallery-anchored
//!   defaults; Architecture C altitude derivation matches Track A/B/D
//!   PNGs exactly).
//! - `dx = dy = 1.0 / grid_size` (derived).
//! - Kinematics preset = Phase 1.1 (only preset wired in Viz-0).
//!
//! Future Viz-0-bis candidates: JSON preset loader (mirror v2 `presets.rs`),
//! Track C kinematics presets, custom `IsostasyConfig`, run-locked
//! closure overrides via UI sliders.

use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
use ymir_core::tectonics_c1::closures::rifting::RiftingParams;
use ymir_core::tectonics_c1::closures::subduction::SubductionParams;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::time_loop::C1Closures;

#[derive(Clone, Debug)]
pub struct C1RunSpec {
    pub grid_size: usize,
    pub seed: u64,
    /// Total simulation steps across all cycles. Rounded DOWN to a
    /// multiple of `steps_per_cycle` by the worker (`n_cycles =
    /// n_steps / steps_per_cycle`, integer division). The UI
    /// exposes the rounded actual step count.
    pub n_steps: usize,
    /// Number of forward-Euler steps per Phase A cycle (A1-c worker
    /// design, Issue #137 Stage A revision). At the end of each
    /// cycle the worker runs `apply_post_tectonic` (sea-level +
    /// macro-redistribution + reclassify) before starting the next
    /// cycle. Smaller `steps_per_cycle` = more frequent coast
    /// reclassification (smoother visual coast migration); larger
    /// = closer to standalone closures (rare coast updates).
    /// Default 50 matches the gallery convention's implicit cycle
    /// scale.
    pub steps_per_cycle: usize,
    pub init_params: Phase2InitParams,
    pub closures: C1Closures,
    pub drainage_max_distance: usize,
}

impl Default for C1RunSpec {
    fn default() -> Self {
        Self {
            grid_size: 64,
            seed: 42,
            n_steps: 300,
            // Default 300 / 50 = 6 cycles. Matches the Phase 1.x
            // workflow tests' cycle scale.
            steps_per_cycle: 50,
            init_params: Phase2InitParams::default(),
            closures: C1Closures::default(),
            drainage_max_distance: 30,
        }
    }
}

impl C1RunSpec {
    /// Effective cycle count (`n_steps / steps_per_cycle`, rounded
    /// down). The worker runs exactly this many cycles; any
    /// remainder of `n_steps` is dropped (no partial cycles).
    pub fn n_cycles(&self) -> usize {
        if self.steps_per_cycle == 0 {
            0
        } else {
            self.n_steps / self.steps_per_cycle
        }
    }

    /// Actual step count after rounding (`n_cycles * steps_per_cycle`).
    /// UI should display this as the "effective" step count.
    pub fn effective_n_steps(&self) -> usize {
        self.n_cycles() * self.steps_per_cycle
    }
}

impl C1RunSpec {
    /// Convenience: a Phase-1.x-style spec with Track D disabled
    /// (subduction / accretion / rifting all `enabled: false`). Useful
    /// for UI testing without dynamic plate mutation.
    #[allow(dead_code)]
    pub fn track_d_disabled() -> Self {
        let mut spec = Self::default();
        spec.closures.subduction = SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        };
        spec.closures.accretion = AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        };
        spec.closures.rifting = RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        };
        spec
    }
}
