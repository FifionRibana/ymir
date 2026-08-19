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
//! closure overrides via UI sliders, opt-in "workflow mode" that
//! invokes `apply_post_tectonic` per cycle (currently OUT OF SCOPE
//! per Issue #137 Stage A revert — Viz-0 mirrors the standalone
//! `run_with_closures` path that produces the Track D visual gallery
//! PNGs).

use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
use ymir_core::tectonics_c1::closures::rifting::RiftingParams;
use ymir_core::tectonics_c1::closures::subduction::SubductionParams;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::time_loop::C1Closures;

#[derive(Clone, Debug)]
pub struct C1RunSpec {
    pub grid_size: usize,
    pub seed: u64,
    pub n_steps: usize,
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
            init_params: Phase2InitParams::default(),
            closures: C1Closures::default(),
            drainage_max_distance: 30,
        }
    }
}

impl C1RunSpec {
    /// M1 #190 production island config — the validated border-clean,
    /// ocean-surrounded continent. `num_plates = 16` + `seed_cluster_count = 3`
    /// (ISOLATION, not fragmentation) with `seed = 9` yield a single non-wrapping
    /// landmass that fits the 338 km export window at 8192² (~41 m/cell) under
    /// `target_land_fraction = 0.08` (set in
    /// [`ymir_core::terrain::upscale::FbmUpscaleConfig::c1_hd_production`]).
    ///
    /// Distinct from [`Self::default`], which stays the 8-plate scientific
    /// baseline the gallery / Stage-A acceptance tests pin. The HD "Générer" path
    /// uses THIS so a default run produces an island, not a torus-spanning blob.
    pub fn island_production() -> Self {
        let mut init_params = Phase2InitParams::default();
        init_params.num_plates = 16;
        init_params.cluster.seed_cluster_count = 3;
        Self { seed: 9, init_params, ..Self::default() }
    }

    /// Convenience: a Phase-1.x-style spec with Track D disabled
    /// (subduction / accretion / rifting all `enabled: false`). Useful
    /// for UI testing without dynamic plate mutation.
    #[allow(dead_code)]
    pub fn track_d_disabled() -> Self {
        let mut spec = Self::default();
        spec.closures.subduction =
            SubductionParams { enabled: false, ..SubductionParams::default() };
        spec.closures.accretion = AccretionParams { enabled: false, ..AccretionParams::default() };
        spec.closures.rifting = RiftingParams { enabled: false, ..RiftingParams::default() };
        spec
    }
}
