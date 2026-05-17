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
/// `CratonicConfig::Disabled`; `On` carries the §9 knobs.
/// Defaults follow the Step 9 Phase 7b validated values.
///
/// Phase 8d added `smoothing_width` and `plate_area_min` to expose
/// the geometry knobs that previously only lived in
/// `CratonicConfigEnabled` defaults. Both gain
/// `#[serde(default = …)]` so preset JSON files written before
/// Phase 8d (which only contain `cr`, `k_viscous`, `b_factor`)
/// continue to load with the core defaults.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum V2CratonicSpec {
    Off,
    On {
        cr: f64,
        k_viscous: f64,
        b_factor: f64,
        #[serde(default = "default_cratonic_smoothing_width")]
        smoothing_width: f64,
        #[serde(default = "default_cratonic_plate_area_min")]
        plate_area_min: f64,
    },
}

fn default_cratonic_smoothing_width() -> f64 {
    // `CratonicConfigEnabled::SMOOTHING_WIDTH_DEFAULT` (= 0.05).
    // Hardcoded here to keep the spec layer free of core constants.
    0.05
}

fn default_cratonic_plate_area_min() -> f64 {
    // `CratonicConfigEnabled::PLATE_AREA_MIN_DEFAULT` (= 0.10).
    0.10
}

impl Default for V2CratonicSpec {
    fn default() -> Self {
        V2CratonicSpec::On {
            cr: 0.3,
            k_viscous: 5.0,
            b_factor: 8.0,
            smoothing_width: default_cratonic_smoothing_width(),
            plate_area_min: default_cratonic_plate_area_min(),
        }
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

/// Step 13 — UI-side mirror of
/// [`ymir_core::tectonics_v2::init::ProfileShape`]. Same
/// `serde(tag = "kind", rename_all = "snake_case")` shape as the
/// core enum so v2 preset JSON round-trips identity.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum V2ProfileShape {
    /// Cubic smoothstep `3t² − 2t³`. Default.
    Smoothstep,
    /// Linear ramp `t`.
    Linear,
    /// Power profile `t^exponent`. UI clamps to `[0.3, 3.0]`.
    Pow { exponent: f64 },
}

impl Default for V2ProfileShape {
    fn default() -> Self {
        V2ProfileShape::Smoothstep
    }
}

impl V2ProfileShape {
    pub fn into_core(self) -> ymir_core::tectonics_v2::init::ProfileShape {
        use ymir_core::tectonics_v2::init::ProfileShape;
        match self {
            V2ProfileShape::Smoothstep => ProfileShape::Smoothstep,
            V2ProfileShape::Linear => ProfileShape::Linear,
            V2ProfileShape::Pow { exponent } => ProfileShape::Pow { exponent },
        }
    }

    pub fn ui_label(&self) -> &'static str {
        match self {
            V2ProfileShape::Smoothstep => "Smoothstep (cubic)",
            V2ProfileShape::Linear => "Linear",
            V2ProfileShape::Pow { .. } => "Pow",
        }
    }

    pub fn variant_index(&self) -> u8 {
        match self {
            V2ProfileShape::Smoothstep => 0,
            V2ProfileShape::Linear => 1,
            V2ProfileShape::Pow { .. } => 2,
        }
    }
}

/// Step 8.6 Phase 8a/8d + Step 13 — S̃ initialisation mode for v2
/// runs. Serialisable mirror of
/// [`ymir_core::tectonics_v2::init::InitMode`]. `Uniform` is the
/// default and matches TDD §4.2's prescription (flat per-plate-type,
/// smoothstep blending across boundaries). `Checkerboard` reproduces
/// the legacy sinusoidal-perturbation pattern bit-for-bit (regression
/// baseline). Step 13 adds `RadialProfile` and `RadialProfileWithFBM`
/// for continental-margin gradient + intra-plate FBM heterogeneity
/// (issue D1/D2).
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum V2InitModeSpec {
    Checkerboard,
    Uniform { boundary_smoothing_width: f64 },
    Gaussian { sigma_continental: f64, sigma_oceanic: f64 },
    Convolution { sigma: f64 },
    /// Step 13 Phase 2 — radial profile per continental plate.
    RadialProfile {
        continental_value: f64,
        oceanic_value: f64,
        profile_shape: V2ProfileShape,
    },
    /// Step 13 Phase 3 — radial profile + isotropic FBM noise on
    /// continental cells. Explicit `rename` overrides serde's
    /// `snake_case` default (which would expand "FBM" into
    /// `f_b_m`) so the on-disk JSON tag matches the core enum's
    /// `radial_profile_with_fbm`.
    ///
    /// Step 13.5 — extended with optional FBM on **oceanic** cells.
    /// Four new fields mirror the core enum's; all carry
    /// `#[serde(default)]` so legacy v2 preset JSON written before
    /// Step 13.5 (which lacks the oceanic keys) loads with the
    /// disabled-default behaviour and remains bit-identical to
    /// Step 13.
    #[serde(rename = "radial_profile_with_fbm")]
    RadialProfileWithFBM {
        continental_value: f64,
        oceanic_value: f64,
        profile_shape: V2ProfileShape,
        fbm_amplitude: f64,
        fbm_octaves: u8,
        fbm_persistence: f64,
        fbm_lacunarity: f64,
        fbm_scale: f64,
        fbm_seed: u64,
        #[serde(default)]
        apply_fbm_to_oceanic: bool,
        #[serde(default = "default_fbm_amplitude_oceanic")]
        fbm_amplitude_oceanic: f64,
        #[serde(default)]
        fbm_scale_oceanic: Option<f64>,
        #[serde(default)]
        fbm_seed_oceanic: Option<u64>,
    },
    /// Step 12 R7.A.1 — orogenic linear ridge per continental plate.
    /// Mirrors [`ymir_core::tectonics_v2::init::InitMode::Orogenic`].
    Orogenic {
        peak_value: f64,
        base_continental_value: f64,
        oceanic_value: f64,
        half_length_ratio: f64,
        width_sigma_ratio: f64,
        orientation: V2OrogenicOrientation,
    },
}

/// V2-side mirror of
/// [`ymir_core::tectonics_v2::init::OrogenicOrientation`]. Serialises
/// to the same `kind`-tagged JSON shape as the core enum so preset
/// round-trips identity.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum V2OrogenicOrientation {
    /// PCA principal axis of the plate's cells (periodic-aware).
    /// Default; falls back to `Fixed { angle_rad: 0.0 }` for plates
    /// with too few cells or rank-1 covariance.
    PlateMainAxisPca,
    /// Constant orientation for every plate.
    Fixed { angle_rad: f64 },
}

impl Default for V2OrogenicOrientation {
    fn default() -> Self {
        V2OrogenicOrientation::PlateMainAxisPca
    }
}

impl V2OrogenicOrientation {
    pub fn into_core(self) -> ymir_core::tectonics_v2::init::OrogenicOrientation {
        use ymir_core::tectonics_v2::init::OrogenicOrientation;
        match self {
            V2OrogenicOrientation::PlateMainAxisPca => OrogenicOrientation::PlateMainAxisPca,
            V2OrogenicOrientation::Fixed { angle_rad } => OrogenicOrientation::Fixed { angle_rad },
        }
    }
}

/// `#[serde(default)]` helper for `V2InitModeSpec::RadialProfile
/// WithFBM::fbm_amplitude_oceanic` — bare `Default::default()` on
/// `f64` is `0.0`, which would silently disable the oceanic FBM
/// perturbation when the user later flips the flag on. The
/// constant from `ymir_core::tectonics_v2::init` keeps the value
/// in one place.
fn default_fbm_amplitude_oceanic() -> f64 {
    ymir_core::tectonics_v2::init::FBM_AMPLITUDE_OCEANIC_DEFAULT
}

impl Default for V2InitModeSpec {
    fn default() -> Self {
        V2InitModeSpec::Uniform { boundary_smoothing_width: 1.0 }
    }
}

impl V2InitModeSpec {
    /// Translate to the core enum the harness understands. Pure
    /// data-shape conversion — no validation (the core
    /// implementation panics if a non-Checkerboard mode is paired
    /// with `BoundaryConfig::Disabled`, but the v2 bridge always
    /// runs with `BoundaryConfig::Enabled` so that path is
    /// unreachable here).
    pub fn into_core(self) -> ymir_core::tectonics_v2::init::InitMode {
        use ymir_core::tectonics_v2::init::InitMode;
        match self {
            V2InitModeSpec::Checkerboard => InitMode::Checkerboard,
            V2InitModeSpec::Uniform { boundary_smoothing_width } => InitMode::Uniform {
                boundary_smoothing_width,
            },
            V2InitModeSpec::Gaussian { sigma_continental, sigma_oceanic } => {
                InitMode::Gaussian { sigma_continental, sigma_oceanic }
            }
            V2InitModeSpec::Convolution { sigma } => InitMode::Convolution { sigma },
            V2InitModeSpec::RadialProfile {
                continental_value,
                oceanic_value,
                profile_shape,
            } => InitMode::RadialProfile {
                continental_value,
                oceanic_value,
                profile_shape: profile_shape.into_core(),
            },
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value,
                oceanic_value,
                profile_shape,
                fbm_amplitude,
                fbm_octaves,
                fbm_persistence,
                fbm_lacunarity,
                fbm_scale,
                fbm_seed,
                // Step 13.5 — oceanic FBM fields, mirrored from
                // the core enum. Phase 2 adds the spec-level
                // mirror; Phase 3 will wire these to the
                // parameter-panel UI controls.
                apply_fbm_to_oceanic,
                fbm_amplitude_oceanic,
                fbm_scale_oceanic,
                fbm_seed_oceanic,
            } => InitMode::RadialProfileWithFBM {
                continental_value,
                oceanic_value,
                profile_shape: profile_shape.into_core(),
                fbm_amplitude,
                fbm_octaves,
                fbm_persistence,
                fbm_lacunarity,
                fbm_scale,
                fbm_seed,
                apply_fbm_to_oceanic,
                fbm_amplitude_oceanic,
                fbm_scale_oceanic,
                fbm_seed_oceanic,
            },
            V2InitModeSpec::Orogenic {
                peak_value,
                base_continental_value,
                oceanic_value,
                half_length_ratio,
                width_sigma_ratio,
                orientation,
            } => InitMode::Orogenic {
                peak_value,
                base_continental_value,
                oceanic_value,
                half_length_ratio,
                width_sigma_ratio,
                orientation: orientation.into_core(),
            },
        }
    }

    /// Short label for the UI dropdown (also used in run-status
    /// strings). Distinct from the variant name so the user-facing
    /// label can drift independently of the on-disk schema.
    pub fn ui_label(&self) -> &'static str {
        match self {
            V2InitModeSpec::Checkerboard => "Checkerboard (legacy sinusoidal)",
            V2InitModeSpec::Uniform { .. } => "Uniform (TDD §4.2 default)",
            V2InitModeSpec::Gaussian { .. } => "Gaussian (peak at centroid)",
            V2InitModeSpec::Convolution { .. } => "Convolution (Gaussian blur)",
            V2InitModeSpec::RadialProfile { .. } => "RadialProfile (Step 13: gradient margins)",
            V2InitModeSpec::RadialProfileWithFBM { .. } => {
                "RadialProfileWithFBM (Step 13: gradient + FBM heterogeneity)"
            }
            V2InitModeSpec::Orogenic { .. } => "Orogenic (Step 12 R7.A.1: linear ridge)",
        }
    }

    /// Discriminant for the UI radio — ignores the inner numeric
    /// payload. Used to detect "user clicked a different mode" so
    /// the UI can swap the parameter widget block.
    pub fn variant_index(&self) -> u8 {
        match self {
            V2InitModeSpec::Checkerboard => 0,
            V2InitModeSpec::Uniform { .. } => 1,
            V2InitModeSpec::Gaussian { .. } => 2,
            V2InitModeSpec::Convolution { .. } => 3,
            V2InitModeSpec::RadialProfile { .. } => 4,
            V2InitModeSpec::RadialProfileWithFBM { .. } => 5,
            V2InitModeSpec::Orogenic { .. } => 6,
        }
    }

    /// Defaults for `RadialProfile` when picked from the dropdown
    /// the first time. Mirrors the core constants
    /// (`CONTINENTAL_VALUE_DEFAULT = 0.95`,
    /// `OCEANIC_VALUE_DEFAULT = 0.20`, `Smoothstep`).
    pub fn radial_profile_default() -> Self {
        V2InitModeSpec::RadialProfile {
            continental_value:
                ymir_core::tectonics_v2::init::CONTINENTAL_VALUE_DEFAULT,
            oceanic_value: ymir_core::tectonics_v2::init::OCEANIC_VALUE_DEFAULT,
            profile_shape: V2ProfileShape::Smoothstep,
        }
    }

    /// Defaults for `RadialProfileWithFBM` when picked from the
    /// dropdown the first time. Mirrors the core `FBM_*_DEFAULT`
    /// constants (amplitude 0.20, scale 0.10 since Step 13 Phase 6;
    /// oceanic FBM disabled by default since Step 13.5).
    pub fn radial_profile_fbm_default() -> Self {
        use ymir_core::tectonics_v2::init::{
            CONTINENTAL_VALUE_DEFAULT, FBM_AMPLITUDE_DEFAULT, FBM_AMPLITUDE_OCEANIC_DEFAULT,
            FBM_LACUNARITY_DEFAULT, FBM_OCTAVES_DEFAULT, FBM_PERSISTENCE_DEFAULT,
            FBM_SCALE_DEFAULT, FBM_SEED_DEFAULT, OCEANIC_VALUE_DEFAULT,
        };
        V2InitModeSpec::RadialProfileWithFBM {
            continental_value: CONTINENTAL_VALUE_DEFAULT,
            oceanic_value: OCEANIC_VALUE_DEFAULT,
            profile_shape: V2ProfileShape::Smoothstep,
            fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
            fbm_octaves: FBM_OCTAVES_DEFAULT,
            fbm_persistence: FBM_PERSISTENCE_DEFAULT,
            fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
            fbm_scale: FBM_SCALE_DEFAULT,
            fbm_seed: FBM_SEED_DEFAULT,
            // Step 13.5 — oceanic FBM disabled by default; the
            // amplitude/scale/seed values are written so a user
            // toggling the flag from the panel sees sensible
            // initial values rather than zeros.
            apply_fbm_to_oceanic: false,
            fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
            fbm_scale_oceanic: None,
            fbm_seed_oceanic: None,
        }
    }

    /// Defaults for `Orogenic` when picked from the dropdown the
    /// first time. Mirrors the core `OROGENIC_*_DEFAULT` constants
    /// (peak=1.20, base=0.85, oceanic=0.20, half_length_ratio=0.40,
    /// width_sigma_ratio=0.08, orientation=PlateMainAxisPca).
    pub fn orogenic_default() -> Self {
        use ymir_core::tectonics_v2::init::{
            OROGENIC_BASE_VALUE_DEFAULT, OROGENIC_HALF_LENGTH_RATIO_DEFAULT,
            OROGENIC_OCEANIC_VALUE_DEFAULT, OROGENIC_PEAK_VALUE_DEFAULT,
            OROGENIC_WIDTH_SIGMA_RATIO_DEFAULT,
        };
        V2InitModeSpec::Orogenic {
            peak_value: OROGENIC_PEAK_VALUE_DEFAULT,
            base_continental_value: OROGENIC_BASE_VALUE_DEFAULT,
            oceanic_value: OROGENIC_OCEANIC_VALUE_DEFAULT,
            half_length_ratio: OROGENIC_HALF_LENGTH_RATIO_DEFAULT,
            width_sigma_ratio: OROGENIC_WIDTH_SIGMA_RATIO_DEFAULT,
            orientation: V2OrogenicOrientation::PlateMainAxisPca,
        }
    }
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

/// Step 12 — Phase A loop parameters mirroring
/// [`ymir_core::tectonics_v2::workflow::PhaseAParams`]. Step 12 R3
/// refactor: dropped `β` (legacy local-deposition coefficient) in
/// favour of `isostatic_rebound_ratio` + `max_drainage_distance` —
/// the new macro-redistribution mechanism is mass-conserving by
/// construction (no `β`-driven mass loss path) and routes eroded
/// sediment along true drainage basins, not just immediate downslope
/// neighbours. See
/// [`ymir_core::tectonics_v2::workflow::macro_redistribution`] for
/// the algorithm.
///
/// Defaults: `α = 0.01`, `rebound = 0.80` (Earth `ρ_crust / ρ_mantle`),
/// `max_drainage_distance = 10` cells, `5 cycles × 20 steps`. R5
/// calibration sweep will refine these.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct V2PhaseAParams {
    pub n_cycles: usize,
    pub k_cycle: usize,
    pub alpha: f64,
    pub isostatic_rebound_ratio: f64,
    pub max_drainage_distance: usize,
}

impl Default for V2PhaseAParams {
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

/// Step 12 — Phase B HD finalization parameters mirroring the
/// user-relevant subset of
/// [`ymir_core::tectonics_v2::workflow::PhaseBParams`]. Only the
/// knobs the panel exposes are roundtripped; the remaining
/// `FbmUpscaleConfig` / `ErosionConfig` defaults are pinned to
/// `core::*::default()` at translation time.
///
/// Default `hd_grid_size = 2048`, `num_droplets = 5_000_000`
/// matches the issue's primary HD target. The `grand_scale_tolerance
/// = 0.10` is the Phase 5 reformulated p95 acceptance threshold.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct V2PhaseBParams {
    pub hd_grid_size: usize,
    pub num_droplets: usize,
    pub erosion_rate: f32,
    pub deposition_rate: f32,
    pub fbm_amplitude_base: f64,
    pub grand_scale_tolerance: f64,
}

impl Default for V2PhaseBParams {
    fn default() -> Self {
        // Matches `ymir_core::erosion::hydraulic::ErosionConfig::default()`
        // for `erosion_rate` / `deposition_rate` / `num_droplets` and
        // `ymir_core::terrain::upscale::FbmUpscaleConfig::default()`
        // for `fbm_amplitude_base`. The 0.10 tolerance is the Phase 5
        // p95 default.
        Self {
            hd_grid_size: 2048,
            num_droplets: 5_000_000,
            erosion_rate: 0.4,
            deposition_rate: 0.35,
            fbm_amplitude_base: 0.08,
            grand_scale_tolerance: 0.10,
        }
    }
}

/// Step 12 — workflow on/off spec mirroring
/// [`ymir_core::tectonics_v2::workflow::WorkflowConfig`]. `Off` is
/// the default for backward compatibility (legacy preset JSON files
/// without a `workflow` field deserialise as `Off` via
/// `#[serde(default)]` on `V2RunSpec.workflow`). `On` carries the
/// Phase A and Phase B parameter bundles.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum V2WorkflowSpec {
    #[default]
    Off,
    On {
        #[serde(default)]
        phase_a: V2PhaseAParams,
        #[serde(default)]
        phase_b: V2PhaseBParams,
    },
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
    /// Step 8.6 Phase 8a/8d — S̃ initialisation mode. Default
    /// `Uniform { boundary_smoothing_width: 1.0 }` aligns with
    /// TDD §4.2 (flat per-plate-type, smoothstep blending). The
    /// `#[serde(default)]` keeps existing preset JSON files
    /// (which predate this field) loading without modification.
    #[serde(default)]
    pub init_mode: V2InitModeSpec,
    /// Step 11 — plate kinematic drift. `Zero` (the default) is a
    /// no-op and bit-identical to pre-Step-11. `PerPlate` carries
    /// per-plate `(vx, vy)` velocities in `[-1, 1]` and a
    /// `boundary_smoothing_width` (cells, typical `1.5`–`6.0`)
    /// blended via smoothstep across inter-plate boundaries. See
    /// `docs/solver-scaling-step11-patch.md` §4.12 for the
    /// deformation/transport split that gives this field its
    /// semantics. `#[serde(default)]` keeps preset JSON files
    /// written before Step 11 loading unchanged (they default to
    /// `Zero`).
    #[serde(default)]
    pub plate_kinematic: V2PlateKinematicSpec,
    /// Step 12 — interleaved tectonic-erosion workflow.
    /// `V2WorkflowSpec::Off` (the default) is a structural no-op:
    /// the bridge runs single `RunBaseline` calls as before.
    /// `V2WorkflowSpec::On { phase_a, phase_b }` enables the multi-
    /// cycle Phase A loop and the HD Phase B finalization, dispatched
    /// by the new `V2Command::RunWorkflowPhaseA` /
    /// `RunWorkflowPhaseB` commands. `#[serde(default)]` keeps preset
    /// JSON files written before Step 12 loading unchanged.
    #[serde(default)]
    pub workflow: V2WorkflowSpec,
}

/// Step 11 — UI-side mirror of
/// [`ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig`].
/// The two enums share `serde(tag = "kind")` shape so a v2 preset
/// JSON round-trips between the panel state and the harness config
/// without bespoke conversion. `Zero` is the default and the
/// pre-Step-11 path; `PerPlate` exposes the per-plate velocity
/// vector + smoothing width to the UI sliders.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum V2PlateKinematicSpec {
    /// All plates at rest — no drift contribution. Bit-identical to
    /// pre-Step-11 paths.
    Zero,
    /// Per-plate velocity assignment with smoothstep blending across
    /// inter-plate boundaries. `velocities[p] = (vx, vy)` for plate
    /// id `p`. `boundary_smoothing_width` is in cells; the panel
    /// clamps it to `[0.5, 5.0]` (≥ `1.5` is the issue D2 default).
    PerPlate {
        velocities: Vec<(f64, f64)>,
        boundary_smoothing_width: f64,
    },
}

impl Default for V2PlateKinematicSpec {
    fn default() -> Self {
        V2PlateKinematicSpec::Zero
    }
}

impl V2PlateKinematicSpec {
    /// Default smoothing width when the user enables `PerPlate`
    /// from the panel. Matches the
    /// `PlateKinematicConfig::DEFAULT_BOUNDARY_SMOOTHING_WIDTH`
    /// constant on the core side so panel→harness round-trips are
    /// identity-default.
    pub const DEFAULT_BOUNDARY_SMOOTHING_WIDTH: f64 = 1.5;

    /// Build a fresh `PerPlate` with `n` plates' velocities all
    /// zeroed — the UI's "enable" path: the toggle becomes the
    /// non-`Zero` variant, but the user has not yet dragged any
    /// slider so the dynamics are equivalent to `Zero` until they
    /// do.
    pub fn per_plate_zero(num_plates: usize) -> Self {
        V2PlateKinematicSpec::PerPlate {
            velocities: vec![(0.0, 0.0); num_plates],
            boundary_smoothing_width: Self::DEFAULT_BOUNDARY_SMOOTHING_WIDTH,
        }
    }

    pub fn is_zero(&self) -> bool {
        matches!(self, V2PlateKinematicSpec::Zero)
    }

    /// Resize the `velocities` vector when the user changes
    /// `num_plates` in the panel. Preserves existing per-plate
    /// values; new slots default to `(0, 0)`. No-op for `Zero`.
    pub fn resize_to(&mut self, num_plates: usize) {
        if let V2PlateKinematicSpec::PerPlate { velocities, .. } = self {
            velocities.resize(num_plates, (0.0, 0.0));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Step 12 Phase 7a — V2WorkflowSpec round-trips through JSON.
    /// Off (default) and On with explicit phase A + phase B params.
    /// Catches schema drift between the panel state and the
    /// harness `WorkflowConfig` the bridge translates it into.
    #[test]
    fn workflow_spec_roundtrips_through_json() {
        let cases = [
            V2WorkflowSpec::Off,
            V2WorkflowSpec::On {
                phase_a: V2PhaseAParams::default(),
                phase_b: V2PhaseBParams::default(),
            },
            V2WorkflowSpec::On {
                phase_a: V2PhaseAParams {
                    n_cycles: 15,
                    k_cycle: 30,
                    alpha: 0.05,
                    isostatic_rebound_ratio: 0.70,
                    max_drainage_distance: 15,
                },
                phase_b: V2PhaseBParams {
                    hd_grid_size: 1024,
                    num_droplets: 1_000_000,
                    erosion_rate: 0.5,
                    deposition_rate: 0.3,
                    fbm_amplitude_base: 0.10,
                    grand_scale_tolerance: 0.10,
                },
            },
        ];
        for original in cases {
            let json = serde_json::to_string(&original).expect("serialize");
            let back: V2WorkflowSpec = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(original, back, "workflow spec roundtrip failed: {json}");
        }
    }

    /// Step 12 Phase 7a — defaults pinned: V2PhaseAParams matches the
    /// D8 issue defaults; V2PhaseBParams matches the issue-prescribed
    /// HD targets. The numerical values are documented to keep
    /// preset / report consumers and the panel's slider ranges in
    /// sync.
    ///
    /// Step 12 R3 refactor: dropped `β`, added
    /// `isostatic_rebound_ratio = 0.80` (Earth `ρ_crust / ρ_mantle`)
    /// and `max_drainage_distance = 10` cells.
    #[test]
    fn workflow_defaults_match_issue_d8() {
        let pa = V2PhaseAParams::default();
        assert_eq!(pa.n_cycles, 5);
        assert_eq!(pa.k_cycle, 20);
        assert_eq!(pa.alpha, 0.01);
        assert_eq!(pa.isostatic_rebound_ratio, 0.80);
        assert_eq!(pa.max_drainage_distance, 10);

        let pb = V2PhaseBParams::default();
        assert_eq!(pb.hd_grid_size, 2048);
        assert_eq!(pb.num_droplets, 5_000_000);
        assert_eq!(pb.grand_scale_tolerance, 0.10);

        assert_eq!(V2WorkflowSpec::default(), V2WorkflowSpec::Off);
    }

    /// Step 12 Phase 7a — legacy preset JSON without a `workflow`
    /// field must still load (defaults to Off via `#[serde(default)]`
    /// on `V2RunSpec.workflow`). This is the bit-identical contract
    /// for pre-Step-12 presets.
    #[test]
    fn legacy_preset_json_without_workflow_loads_with_off() {
        let json = r#"{
            "seed": 42,
            "grid_nx": 64,
            "grid_ny": 64,
            "steps": 100,
            "num_plates": 8,
            "continental_ratio": 0.3,
            "bi": 0.15,
            "br": 0.05,
            "mantle": { "kind": "off" },
            "slab_enabled": false,
            "cratonic": { "kind": "on", "cr": 0.3, "k_viscous": 5.0, "b_factor": 8.0 },
            "age_field": { "kind": "off" },
            "linear_solver": "jacobi",
            "force": { "kind": "gpe" },
            "s_perturbation_amplitude": 0.2,
            "total_time_nondim": 6.0,
            "cfl_factor": 0.3
        }"#;
        let recovered: V2RunSpec = serde_json::from_str(json).expect("legacy preset must load");
        assert_eq!(recovered.workflow, V2WorkflowSpec::Off);
    }

    /// Step 11 — round-trip every `V2PlateKinematicSpec` variant
    /// through JSON. Catches schema drift between the panel state
    /// and the harness config the bridge translates it into.
    #[test]
    fn plate_kinematic_spec_roundtrips_through_json() {
        let cases = [
            V2PlateKinematicSpec::Zero,
            V2PlateKinematicSpec::PerPlate {
                velocities: vec![(0.5, 0.0), (-0.5, 0.0), (0.0, 0.3)],
                boundary_smoothing_width: 1.5,
            },
            V2PlateKinematicSpec::PerPlate {
                velocities: vec![(0.0, 0.0); 8],
                boundary_smoothing_width: 4.5,
            },
        ];
        for original in cases {
            let json = serde_json::to_string(&original).expect("serialize");
            let back: V2PlateKinematicSpec = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(original, back, "roundtrip failed: {}", json);
        }
    }

    /// Step 11 — backward-compat: a preset JSON written before
    /// Step 11 (no `plate_kinematic` field) must still deserialise
    /// and default to `Zero` so the harness keeps the bit-identical
    /// regression contract.
    #[test]
    fn old_preset_without_plate_kinematic_defaults_to_zero() {
        let old_json = r#"{
            "seed": 42,
            "grid_nx": 64,
            "grid_ny": 64,
            "steps": 100,
            "num_plates": 8,
            "continental_ratio": 0.3,
            "bi": 0.15,
            "br": 0.05,
            "mantle": { "kind": "off" },
            "slab_enabled": false,
            "cratonic": { "kind": "off" },
            "age_field": { "kind": "off" },
            "linear_solver": "jacobi",
            "force": { "kind": "gpe" },
            "s_perturbation_amplitude": 0.2,
            "total_time_nondim": 6.0,
            "cfl_factor": 0.3
        }"#;
        let spec: V2RunSpec =
            serde_json::from_str(old_json).expect("old preset must deserialize");
        assert_eq!(spec.plate_kinematic, V2PlateKinematicSpec::Zero);
    }

    /// Phase 8d + Step 13 Phase 5 — round-trip every
    /// `V2InitModeSpec` variant through JSON. Catches schema drift
    /// between the Rust enum and any downstream JSON the user
    /// hand-edits.
    #[test]
    fn init_mode_spec_roundtrips_through_json() {
        let cases = [
            V2InitModeSpec::Checkerboard,
            V2InitModeSpec::Uniform { boundary_smoothing_width: 1.5 },
            V2InitModeSpec::Gaussian { sigma_continental: 6.0, sigma_oceanic: 4.5 },
            V2InitModeSpec::Convolution { sigma: 2.25 },
            V2InitModeSpec::RadialProfile {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Smoothstep,
            },
            V2InitModeSpec::RadialProfile {
                continental_value: 0.90,
                oceanic_value: 0.18,
                profile_shape: V2ProfileShape::Pow { exponent: 2.0 },
            },
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Smoothstep,
                fbm_amplitude: 0.10,
                fbm_octaves: 4,
                fbm_persistence: 0.5,
                fbm_lacunarity: 2.0,
                fbm_scale: 0.25,
                fbm_seed: 0x0FBA_5EED,
                // Step 13.5 — disabled-default oceanic FBM
                // exercises the most common path on roundtrip.
                apply_fbm_to_oceanic: false,
                fbm_amplitude_oceanic: 0.10,
                fbm_scale_oceanic: None,
                fbm_seed_oceanic: None,
            },
            // Step 13.5 — second case with oceanic FBM enabled
            // and explicit `Some(...)` for the optional fields,
            // so the roundtrip exercises both branches of the
            // serde encoding.
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Smoothstep,
                fbm_amplitude: 0.20,
                fbm_octaves: 4,
                fbm_persistence: 0.5,
                fbm_lacunarity: 2.0,
                fbm_scale: 0.10,
                fbm_seed: 0x0FBA_5EED,
                apply_fbm_to_oceanic: true,
                fbm_amplitude_oceanic: 0.12,
                fbm_scale_oceanic: Some(0.08),
                fbm_seed_oceanic: Some(0xC0FFEE_5EE_D),
            },
            // Step 12 R7.A.1 — orogenic mode, PCA orientation.
            V2InitModeSpec::Orogenic {
                peak_value: 1.20,
                base_continental_value: 0.85,
                oceanic_value: 0.20,
                half_length_ratio: 0.40,
                width_sigma_ratio: 0.08,
                orientation: V2OrogenicOrientation::PlateMainAxisPca,
            },
            // Step 12 R7.A.1 — orogenic mode, Fixed orientation
            // (exercises the second variant of OrogenicOrientation).
            V2InitModeSpec::Orogenic {
                peak_value: 1.20,
                base_continental_value: 0.85,
                oceanic_value: 0.20,
                half_length_ratio: 0.40,
                width_sigma_ratio: 0.08,
                orientation: V2OrogenicOrientation::Fixed { angle_rad: 0.7853981633974483 },
            },
        ];
        for original in cases {
            let json = serde_json::to_string(&original).expect("serialize");
            let back: V2InitModeSpec = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(original, back, "roundtrip failed: {}", json);
        }
    }

    /// Step 13 Phase 5 — explicit roundtrip test for the new
    /// `RadialProfile` / `RadialProfileWithFBM` variants. Covers
    /// every `ProfileShape` (Smoothstep, Linear, Pow) and asserts
    /// the on-disk JSON shape matches the documented schema (so a
    /// user hand-editing a preset gets predictable structure).
    #[test]
    fn v2_panel_radial_modes_serde_roundtrip() {
        // RadialProfile × every ProfileShape.
        let radial_cases = [
            V2InitModeSpec::RadialProfile {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Smoothstep,
            },
            V2InitModeSpec::RadialProfile {
                continental_value: 0.90,
                oceanic_value: 0.15,
                profile_shape: V2ProfileShape::Linear,
            },
            V2InitModeSpec::RadialProfile {
                continental_value: 0.85,
                oceanic_value: 0.25,
                profile_shape: V2ProfileShape::Pow { exponent: 0.5 },
            },
            V2InitModeSpec::RadialProfile {
                continental_value: 1.0,
                oceanic_value: 0.10,
                profile_shape: V2ProfileShape::Pow { exponent: 3.0 },
            },
        ];
        for c in radial_cases {
            let json = serde_json::to_string(&c).expect("serialize");
            let back: V2InitModeSpec = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(c, back, "RadialProfile roundtrip failed: {}", json);
        }

        // RadialProfileWithFBM with non-default FBM parameters to
        // catch any field-order mismatch between serialize and
        // deserialize.
        let fbm_case = V2InitModeSpec::RadialProfileWithFBM {
            continental_value: 0.92,
            oceanic_value: 0.18,
            profile_shape: V2ProfileShape::Pow { exponent: 1.7 },
            fbm_amplitude: 0.18,
            fbm_octaves: 6,
            fbm_persistence: 0.65,
            fbm_lacunarity: 2.3,
            fbm_scale: 0.18,
            fbm_seed: 0xCAFE_F00D,
            // Step 13.5 — disabled-default oceanic FBM in this
            // probe; the dedicated `v2_panel_radial_fbm_with_
            // oceanic_roundtrip` exercises the enabled path.
            apply_fbm_to_oceanic: false,
            fbm_amplitude_oceanic: 0.10,
            fbm_scale_oceanic: None,
            fbm_seed_oceanic: None,
        };
        let json = serde_json::to_string(&fbm_case).expect("serialize");
        let back: V2InitModeSpec = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(fbm_case, back, "RadialProfileWithFBM roundtrip failed: {}", json);

        // Schema probe: deserialize a hand-written preset fragment
        // matching the documented field shape (catches accidental
        // tag-renaming or field-ordering drift). This fragment
        // omits the Step 13.5 oceanic FBM fields, so the parsed
        // value carries the disabled defaults (`apply_fbm_to_
        // oceanic = false`, etc.) — an implicit legacy-preset
        // test; the dedicated
        // `v2_panel_radial_fbm_legacy_preset_load` covers the
        // contract explicitly.
        let hand_written = r#"{
            "kind": "radial_profile_with_fbm",
            "continental_value": 0.95,
            "oceanic_value": 0.20,
            "profile_shape": { "kind": "pow", "exponent": 2.0 },
            "fbm_amplitude": 0.10,
            "fbm_octaves": 4,
            "fbm_persistence": 0.5,
            "fbm_lacunarity": 2.0,
            "fbm_scale": 0.25,
            "fbm_seed": 264339693
        }"#;
        let parsed: V2InitModeSpec =
            serde_json::from_str(hand_written).expect("hand-written preset must parse");
        assert_eq!(
            parsed,
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Pow { exponent: 2.0 },
                fbm_amplitude: 0.10,
                fbm_octaves: 4,
                fbm_persistence: 0.5,
                fbm_lacunarity: 2.0,
                fbm_scale: 0.25,
                fbm_seed: 264339693,
                // Step 13.5 — these are the `#[serde(default)]`
                // values applied when the JSON omits the oceanic
                // FBM keys.
                apply_fbm_to_oceanic: false,
                fbm_amplitude_oceanic:
                    ymir_core::tectonics_v2::init::FBM_AMPLITUDE_OCEANIC_DEFAULT,
                fbm_scale_oceanic: None,
                fbm_seed_oceanic: None,
            }
        );
    }

    /// Step 13.5 Phase 2 — legacy v2 preset JSON written before the
    /// oceanic FBM extension must still deserialise without error,
    /// with all four new fields populated by their
    /// `#[serde(default)]` values: `apply_fbm_to_oceanic = false`,
    /// `fbm_amplitude_oceanic = FBM_AMPLITUDE_OCEANIC_DEFAULT`,
    /// `fbm_scale_oceanic = None`, `fbm_seed_oceanic = None`.
    /// `into_core()` converts the spec to an `InitMode` whose
    /// disabled flag short-circuits the oceanic FBM block — the
    /// run is bit-identical to its Step 13 form.
    #[test]
    fn v2_panel_radial_fbm_legacy_preset_load() {
        let legacy_json = r#"{
            "kind": "radial_profile_with_fbm",
            "continental_value": 0.95,
            "oceanic_value": 0.20,
            "profile_shape": { "kind": "smoothstep" },
            "fbm_amplitude": 0.20,
            "fbm_octaves": 4,
            "fbm_persistence": 0.5,
            "fbm_lacunarity": 2.0,
            "fbm_scale": 0.10,
            "fbm_seed": 264339693
        }"#;
        let parsed: V2InitModeSpec = serde_json::from_str(legacy_json)
            .expect("Step 13 preset must still load with Step 13.5's spec");
        match parsed {
            V2InitModeSpec::RadialProfileWithFBM {
                apply_fbm_to_oceanic,
                fbm_amplitude_oceanic,
                fbm_scale_oceanic,
                fbm_seed_oceanic,
                ..
            } => {
                assert!(
                    !apply_fbm_to_oceanic,
                    "legacy preset must default `apply_fbm_to_oceanic` to false"
                );
                assert_eq!(
                    fbm_amplitude_oceanic,
                    ymir_core::tectonics_v2::init::FBM_AMPLITUDE_OCEANIC_DEFAULT,
                    "legacy preset must default `fbm_amplitude_oceanic` to the placeholder constant"
                );
                assert_eq!(
                    fbm_scale_oceanic, None,
                    "legacy preset must default `fbm_scale_oceanic` to None"
                );
                assert_eq!(
                    fbm_seed_oceanic, None,
                    "legacy preset must default `fbm_seed_oceanic` to None"
                );
            }
            other => panic!("expected RadialProfileWithFBM, got {:?}", other),
        }

        // The parsed spec must convert to a core `InitMode` whose
        // disabled flag short-circuits the oceanic FBM path: the
        // resulting run is bit-identical to its Step 13 form.
        let core = parsed.into_core();
        match core {
            ymir_core::tectonics_v2::init::InitMode::RadialProfileWithFBM {
                apply_fbm_to_oceanic,
                ..
            } => assert!(
                !apply_fbm_to_oceanic,
                "into_core() must thread the disabled flag through the conversion"
            ),
            other => panic!("expected InitMode::RadialProfileWithFBM, got {:?}", other),
        }
    }

    /// Step 13.5 Phase 2 — JSON with the four new oceanic FBM
    /// fields explicitly populated must roundtrip
    /// (serialize → deserialize → compare) byte-for-byte. Locks
    /// the on-disk schema for the oceanic FBM extension.
    #[test]
    fn v2_panel_radial_fbm_with_oceanic_roundtrip() {
        let cases = [
            // Disabled flag with explicit oceanic params filled —
            // the params are written but unused by the run; this
            // case exercises the serde encoding of the fields.
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Smoothstep,
                fbm_amplitude: 0.20,
                fbm_octaves: 4,
                fbm_persistence: 0.5,
                fbm_lacunarity: 2.0,
                fbm_scale: 0.10,
                fbm_seed: 0x0FBA_5EED,
                apply_fbm_to_oceanic: false,
                fbm_amplitude_oceanic: 0.15,
                fbm_scale_oceanic: Some(0.12),
                fbm_seed_oceanic: Some(0xC0FFEE),
            },
            // Enabled flag with non-default oceanic params, all
            // four fields written.
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value: 0.92,
                oceanic_value: 0.18,
                profile_shape: V2ProfileShape::Pow { exponent: 2.0 },
                fbm_amplitude: 0.20,
                fbm_octaves: 4,
                fbm_persistence: 0.5,
                fbm_lacunarity: 2.0,
                fbm_scale: 0.10,
                fbm_seed: 0xCAFE_F00D,
                apply_fbm_to_oceanic: true,
                fbm_amplitude_oceanic: 0.08,
                fbm_scale_oceanic: Some(0.05),
                fbm_seed_oceanic: Some(0xDEADBEEF),
            },
            // Enabled flag with `None` for the optional fields,
            // exercising the default-derivation path of the spec.
            V2InitModeSpec::RadialProfileWithFBM {
                continental_value: 0.95,
                oceanic_value: 0.20,
                profile_shape: V2ProfileShape::Smoothstep,
                fbm_amplitude: 0.20,
                fbm_octaves: 4,
                fbm_persistence: 0.5,
                fbm_lacunarity: 2.0,
                fbm_scale: 0.10,
                fbm_seed: 0x0FBA_5EED,
                apply_fbm_to_oceanic: true,
                fbm_amplitude_oceanic: 0.10,
                fbm_scale_oceanic: None,
                fbm_seed_oceanic: None,
            },
        ];
        for original in cases {
            let json = serde_json::to_string(&original).expect("serialize");
            let back: V2InitModeSpec = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(
                original, back,
                "Step 13.5 RadialProfileWithFBM roundtrip with oceanic FBM failed: {}",
                json
            );

            // Cross-check the serialised JSON contains all four
            // oceanic FBM keys when the flag is on (lock the
            // documented schema; if the user toggles the flag
            // and edits a preset by hand they need the keys
            // visible).
            if let V2InitModeSpec::RadialProfileWithFBM {
                apply_fbm_to_oceanic: true,
                ..
            } = original
            {
                assert!(json.contains("\"apply_fbm_to_oceanic\":true"),
                    "serialised JSON must contain apply_fbm_to_oceanic=true: {}", json);
                assert!(json.contains("\"fbm_amplitude_oceanic\""),
                    "serialised JSON must contain fbm_amplitude_oceanic: {}", json);
            }
        }
    }

    /// Phase 8d — full-run-spec round-trip with the new `init_mode`
    /// + extended cratonic fields, plus a backward-compat probe:
    /// older preset JSON without `init_mode` must still deserialize
    /// (defaults to Uniform via `#[serde(default)]`).
    #[test]
    fn run_spec_roundtrips_and_old_json_loads_with_uniform_default() {
        let original = V2RunSpec::active_medley_defaults();
        let json = serde_json::to_string(&original).expect("serialize");
        let back: V2RunSpec = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back.init_mode, V2InitModeSpec::default());

        // Backward-compat: a preset JSON written before Phase 8d
        // (no `init_mode`, no `smoothing_width`/`plate_area_min` on
        // cratonic.On) must still deserialise.
        let old_json = r#"{
            "seed": 42,
            "grid_nx": 64,
            "grid_ny": 64,
            "steps": 100,
            "num_plates": 8,
            "continental_ratio": 0.3,
            "bi": 0.15,
            "br": 0.05,
            "mantle": { "kind": "off" },
            "slab_enabled": false,
            "cratonic": { "kind": "on", "cr": 0.3, "k_viscous": 5.0, "b_factor": 8.0 },
            "age_field": { "kind": "off" },
            "linear_solver": "jacobi",
            "force": { "kind": "gpe" },
            "s_perturbation_amplitude": 0.2,
            "total_time_nondim": 6.0,
            "cfl_factor": 0.3
        }"#;
        let recovered: V2RunSpec =
            serde_json::from_str(old_json).expect("legacy preset must still load");
        assert_eq!(recovered.init_mode, V2InitModeSpec::default());
        match recovered.cratonic {
            V2CratonicSpec::On { smoothing_width, plate_area_min, .. } => {
                assert!((smoothing_width - 0.05).abs() < 1e-12);
                assert!((plate_area_min - 0.10).abs() < 1e-12);
            }
            _ => panic!("expected cratonic On"),
        }
    }
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
            init_mode: V2InitModeSpec::default(),
            plate_kinematic: V2PlateKinematicSpec::default(),
            workflow: V2WorkflowSpec::default(),
        }
    }
}
