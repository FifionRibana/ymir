//! Step 8.6 — light-weight, cloneable run specification consumed by
//! the v2 bridge thread.
//!
//! The v2 harness `BaselineConfig` carries a `Box<dyn BodyForce>` which
//! is not `Clone`, so we cannot send it across a channel. Instead, the
//! bridge accepts this spec — a flat, serialisable, `Clone`-able
//! description of every v2 knob the UI exposes — and rebuilds a full
//! `BaselineConfig` inside the worker thread (`build_config::build`).
//!
//! `Serialize` / `Deserialize` enables Phase 4 preset files
//! (`presets/v2/*.json`) so users can share named tectonic
//! configurations.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

/// Mantle-forcing spec (Step 8). `Off` collapses to
/// `MantleConfig::Disabled`; `On` expands into the full
/// `MantleConfig::Enabled { .. }` variant.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum V2MantleSpec {
    Off,
    On {
        mf: f64,
        coupling: f64,
        num_modes: usize,
        seed: u64,
        evolution_rate: f64,
    },
}

impl Default for V2MantleSpec {
    fn default() -> Self {
        // `MF_DEFAULT = 1.0`, `COUPLING_DEFAULT = 1.0`,
        // `NUM_MODES_DEFAULT = 6` per `tectonics_v2::mantle`. Mantle
        // seed `7` matches the Step 8 baseline preset.
        V2MantleSpec::On {
            mf: 1.0,
            coupling: 1.0,
            num_modes: 6,
            seed: 7,
            evolution_rate: 0.0,
        }
    }
}

/// Cratonic-immunity spec (Step 9). `Off` collapses to
/// `CratonicConfig::Disabled`; `On` carries the three §9 knobs.
/// Defaults follow the Step 9 Phase 7b validated values.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum V2CratonicSpec {
    Off,
    On {
        cr: f64,
        k_viscous: f64,
        b_factor: f64,
    },
}

impl Default for V2CratonicSpec {
    fn default() -> Self {
        V2CratonicSpec::On { cr: 0.3, k_viscous: 5.0, b_factor: 8.0 }
    }
}

/// Age-field spec (Step 10). `Off` collapses to
/// `AgeFieldConfig::Disabled`. Default initial ages from §4.11.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum V2AgeFieldSpec {
    Off,
    On {
        continental_age_init: f64,
        oceanic_age_init: f64,
    },
}

impl Default for V2AgeFieldSpec {
    fn default() -> Self {
        V2AgeFieldSpec::On {
            continental_age_init: 7.0,
            oceanic_age_init: 0.5,
        }
    }
}

/// Linear-solver dispatch (Step 8.5a Phase 4.3). `Jacobi` is the
/// safe Step 9 default; `Amg` opts into Option B' (Picard-block
/// V-cycle, matrix-free Newton tangent).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum V2LinearSolverSpec {
    #[default]
    Jacobi,
    Amg,
}

/// Body-force scenario. `Gpe` is the physics default; `Sinusoidal`
/// is kept for the regression preset.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum V2ForceKind {
    Gpe,
    Sinusoidal { amplitude: f64 },
}

impl Default for V2ForceKind {
    fn default() -> Self {
        V2ForceKind::Gpe
    }
}

fn default_output_dir() -> PathBuf {
    std::env::temp_dir().join("ymir_v2_run")
}

fn default_preset_label() -> String {
    "unnamed".to_string()
}

fn default_capture_endpoints() -> bool {
    false
}

/// Full v2 run specification — every knob the UI exposes plus the
/// scratch directory for PNG snapshots (`heightmap_fractions` is
/// derived from `capture_at_end`: `[]` for "no PNGs" or `[0.0, 1.0]`
/// for first/last frame). The bridge thread converts this to a
/// fully-populated `BaselineConfig` via [`crate::bridge::v2::build_config::build`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct V2RunSpec {
    pub seed: u64,
    pub grid_nx: usize,
    pub grid_ny: usize,
    pub steps: usize,
    pub num_plates: usize,
    pub continental_ratio: f64,
    pub bi: f64,
    pub br: f64,
    pub mantle: V2MantleSpec,
    pub slab_enabled: bool,
    pub cratonic: V2CratonicSpec,
    pub age_field: V2AgeFieldSpec,
    pub linear_solver: V2LinearSolverSpec,
    pub force: V2ForceKind,
    pub s_perturbation_amplitude: f64,
    pub total_time_nondim: f64,
    pub cfl_factor: f64,
    /// Capture first + last PNG snapshots (S̃ + age) under
    /// `output_dir`. Set to `false` to skip disk I/O entirely
    /// (e.g. tests, interactive UI runs that consume the
    /// in-memory `FinalState` directly). Defaulted to `false` in
    /// preset files so a freshly-loaded preset stays interactive.
    #[serde(default = "default_capture_endpoints")]
    pub capture_endpoints: bool,
    /// Output directory for harness-side PNG dumps. Created on
    /// demand. Ignored when `capture_endpoints == false`.
    #[serde(default = "default_output_dir")]
    pub output_dir: PathBuf,
    /// Human-readable preset name (e.g. "active_medley"). Carried
    /// through to the boundary-layout-name field for downstream
    /// reports.
    #[serde(default = "default_preset_label")]
    pub preset_label: String,
}

impl V2RunSpec {
    /// "Active medley" preset — Step 8 shape with cratonic on and
    /// age-field on, the §4.11 validated regime. 64² × 100 steps is
    /// the canonical baseline; tests pass smaller `(grid, steps)`
    /// values to keep wallclock low.
    pub fn active_medley_defaults() -> Self {
        Self {
            seed: 42,
            grid_nx: 64,
            grid_ny: 64,
            steps: 100,
            num_plates: 8,
            continental_ratio: 0.3,
            bi: 0.15,
            br: 0.05,
            mantle: V2MantleSpec::default(),
            slab_enabled: false,
            cratonic: V2CratonicSpec::default(),
            age_field: V2AgeFieldSpec::default(),
            linear_solver: V2LinearSolverSpec::Jacobi,
            force: V2ForceKind::Gpe,
            s_perturbation_amplitude: 0.2,
            total_time_nondim: 6.0,
            cfl_factor: 0.3,
            capture_endpoints: false,
            output_dir: default_output_dir(),
            preset_label: "active_medley".to_string(),
        }
    }
}
