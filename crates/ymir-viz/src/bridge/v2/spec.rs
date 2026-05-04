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

/// Step 8.6 Phase 8a/8d — S̃ initialisation mode for v2 runs.
/// Serialisable mirror of [`ymir_core::tectonics_v2::init::InitMode`].
/// `Uniform` is the default and matches TDD §4.2's prescription
/// (flat per-plate-type, smoothstep blending across boundaries).
/// `Checkerboard` reproduces the legacy sinusoidal-perturbation
/// pattern bit-for-bit (regression baseline).
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum V2InitModeSpec {
    Checkerboard,
    Uniform { boundary_smoothing_width: f64 },
    Gaussian { sigma_continental: f64, sigma_oceanic: f64 },
    Convolution { sigma: f64 },
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

    /// Phase 8d — round-trip every `V2InitModeSpec` variant through
    /// JSON. Catches schema drift between the Rust enum and any
    /// downstream JSON the user hand-edits.
    #[test]
    fn init_mode_spec_roundtrips_through_json() {
        let cases = [
            V2InitModeSpec::Checkerboard,
            V2InitModeSpec::Uniform { boundary_smoothing_width: 1.5 },
            V2InitModeSpec::Gaussian { sigma_continental: 6.0, sigma_oceanic: 4.5 },
            V2InitModeSpec::Convolution { sigma: 2.25 },
        ];
        for original in cases {
            let json = serde_json::to_string(&original).expect("serialize");
            let back: V2InitModeSpec = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(original, back, "roundtrip failed: {}", json);
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
        }
    }
}
