use bevy::prelude::*;
use ymir_core::erosion::hydraulic::{ErosionConfig, ErosionStats};
use ymir_core::grid::GridF32;
use ymir_core::lakes::detection::{LakeConfig, LakeResult};
use ymir_core::tectonics::boundaries::BoundaryConfig;
use ymir_core::tectonics::plates::{PlateConfig, PlateInitResult};
use ymir_core::tectonics::solver::config::{NonlinearSolver, Preconditioner};
use ymir_core::terrain::flow::{FlowConfig, FlowResult, RiverConfig, RiverNetwork};

// ── View ─────────────────────────────────────────────────────────────────

#[derive(States, PartialEq, Eq, Clone, Copy, Default, Debug, Hash)]
pub enum ViewMode {
    #[default]
    Altitude,
    Slope,
    Tectonics,
    Biomes,
    Climate,
    Geology,
}

impl ViewMode {
    pub const ALL: &[ViewMode] = &[
        ViewMode::Altitude,
        ViewMode::Slope,
        ViewMode::Tectonics,
        ViewMode::Biomes,
        ViewMode::Climate,
        ViewMode::Geology,
    ];

    pub fn label(self) -> &'static str {
        match self {
            ViewMode::Altitude => "Altitude",
            ViewMode::Slope => "Slope",
            ViewMode::Tectonics => "Tectonics",
            ViewMode::Biomes => "Biomes",
            ViewMode::Climate => "Climate",
            ViewMode::Geology => "Geology",
        }
    }

    pub fn is_enabled(self) -> bool {
        matches!(self, ViewMode::Altitude | ViewMode::Slope | ViewMode::Tectonics)
    }
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct OverlayFlags {
    pub rivers: bool,
    pub lakes: bool,
    pub hillshade: bool,
    pub grid: bool,
    pub plates: bool,
}

impl Default for OverlayFlags {
    fn default() -> Self {
        Self { rivers: false, lakes: false, hillshade: true, grid: false, plates: false }
    }
}

#[derive(Resource, Clone, Debug, Default)]
pub struct ViewState {
    pub overlays: OverlayFlags,
}

// ── Pipeline ─────────────────────────────────────────────────────────────

#[derive(States, PartialEq, Eq, Clone, Copy, Debug, Default, Hash)]
pub enum PipelinePhase {
    Tectonics,
    Isostasy,
    UpscaleFbm,
    #[default]
    Erosion,
    Climate,
    Hydrology,
}

impl PipelinePhase {
    pub const ALL: &[PipelinePhase] = &[
        PipelinePhase::Tectonics,
        PipelinePhase::Isostasy,
        PipelinePhase::UpscaleFbm,
        PipelinePhase::Erosion,
        PipelinePhase::Hydrology,
        PipelinePhase::Climate,
    ];

    pub fn label(self) -> &'static str {
        match self {
            PipelinePhase::Tectonics => "Tectonics",
            PipelinePhase::Isostasy => "Isostasy",
            PipelinePhase::UpscaleFbm => "Upscale + FBM",
            PipelinePhase::Erosion => "Erosion",
            PipelinePhase::Hydrology => "Hydrology",
            PipelinePhase::Climate => "Climate",
        }
    }

    pub fn short_label(self) -> &'static str {
        match self {
            PipelinePhase::Tectonics => "TEC",
            PipelinePhase::Isostasy => "ISO",
            PipelinePhase::UpscaleFbm => "FBM",
            PipelinePhase::Erosion => "ERO",
            PipelinePhase::Hydrology => "HYD",
            PipelinePhase::Climate => "CLI",
        }
    }
}

#[derive(PartialEq, Eq, Clone, Copy, Debug, Default)]
#[allow(dead_code)]
pub enum PhaseStatus {
    #[default]
    Pending,
    Running(u8),
    Completed,
}

impl PhaseStatus {
    pub fn icon(self) -> &'static str {
        match self {
            PhaseStatus::Pending => "⚫",
            PhaseStatus::Running(_) => "🟡",
            PhaseStatus::Completed => "🟢",
        }
    }
}

#[derive(Resource, Clone, Debug)]
pub struct PipelineState {
    pub phases: Vec<(PipelinePhase, PhaseStatus, String)>,
}

impl Default for PipelineState {
    fn default() -> Self {
        Self {
            phases: PipelinePhase::ALL
                .iter()
                .map(|&p| (p, PhaseStatus::Pending, String::new()))
                .collect(),
        }
    }
}

// ── Parameters ───────────────────────────────────────────────────────────

#[derive(Resource, Clone, Debug)]
pub struct ErosionParams {
    pub erosion_rate: f32,
    pub deposition_rate: f32,
    pub inertia: f32,
    pub gravity: f32,
    pub evaporation_rate: f32,
    pub max_lifetime: u32,
    pub droplets_millions: f32,
    pub coastal_deposition: u32,
    pub min_slope: f32,
}

impl Default for ErosionParams {
    fn default() -> Self {
        Self {
            erosion_rate: 0.4,
            deposition_rate: 0.35,
            inertia: 0.08,
            gravity: 6.0,
            evaporation_rate: 0.015,
            max_lifetime: 150,
            droplets_millions: 5.0,
            coastal_deposition: 12,
            min_slope: 0.001,
        }
    }
}

/// Configuration for the thin viscous sheet solver (UI-facing).
#[derive(Resource, Clone, Debug)]
pub struct SolverConfig {
    pub num_timesteps: usize,
    pub gravity_factor: f64,
    pub cfl_factor: f64,
    pub power_law_n: f64,
    pub picard_relaxation: f64,
    pub nonlinear_solver: NonlinearSolver,
    pub continuation_enabled: bool,
    pub strain_rate_min: f64,
    pub eta_max: f64,
    pub preconditioner: Preconditioner,
    pub inexact_newton: bool,
    pub boundaries: BoundaryConfig,
    pub dynamic_boundaries: bool,
    pub cratonic: ymir_core::tectonics::solver::config::CratonicConfig,
    pub yielding: ymir_core::tectonics::solver::config::YieldingConfig,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            num_timesteps: 300,
            gravity_factor: 1.0,
            cfl_factor: 0.5,
            power_law_n: 3.0,
            picard_relaxation: 0.7,
            nonlinear_solver: NonlinearSolver::Newton,
            continuation_enabled: true,
            strain_rate_min: 1e-3,
            eta_max: 1e4,
            preconditioner: Preconditioner::Ssor { omega: 1.2 },
            inexact_newton: true,
            boundaries: BoundaryConfig::default(),
            dynamic_boundaries: true,
            cratonic: Default::default(),
            yielding: Default::default(),
        }
    }
}

#[derive(Resource, Clone, Debug)]
pub struct ClimateParams {
    pub base_temperature: f32,
    pub wind_direction_deg: f32,
    pub orographic_factor: f32,
    pub moisture_decay_km: f32,
}

impl Default for ClimateParams {
    fn default() -> Self {
        Self {
            base_temperature: 12.0,
            wind_direction_deg: 240.0,
            orographic_factor: 3.0,
            moisture_decay_km: 400.0,
        }
    }
}

#[derive(Resource, Clone, Debug)]
pub struct GenerationParamsUi {
    pub seed: u64,
    pub meters_per_pixel: f32,
    pub continent_size_km: f32,
    pub max_elevation_m: f32,
}

impl Default for GenerationParamsUi {
    fn default() -> Self {
        Self { seed: 42, meters_per_pixel: 40.0, continent_size_km: 300.0, max_elevation_m: 4000.0 }
    }
}

// ── Terrain data ─────────────────────────────────────────────────────────

#[derive(Resource)]
#[allow(dead_code)]
pub struct TerrainData {
    pub heightmap: GridF32,
    pub max_elevation: f32,
    pub dirty: bool,
}

#[derive(Resource, Clone, Debug, Default)]
pub struct TerrainStats {
    pub grid_width: usize,
    pub grid_height: usize,
    pub meters_per_pixel: f32,
    pub peak_altitude: f32,
    pub min_altitude: f32,
    pub land_ratio: f32,
    pub river_segments: usize,
    pub lake_count: usize,
}

// ── Tectonic state ────────────────────────────────────────────────────────

/// Live state for the tectonic plate visualization.
/// Regenerated on demand via the "Generate plates" UI button.
#[derive(Resource)]
pub struct TectonicState {
    pub init: PlateInitResult,
    pub config: PlateConfig,
    pub seed: u64,
    /// True when the GPU texture needs to be rebuilt from `init`.
    pub dirty: bool,
    /// Incremented on every regeneration — lets other systems detect stale caches
    /// without relying on Bevy's single-frame `is_changed()` window.
    pub generation: u64,
}

// ── Isostasy ─────────────────────────────────────────────────────────

#[derive(Resource, Clone, Debug)]
pub struct IsostasyParams {
    pub sea_level_fraction: f32,
    pub max_elevation_m: f32,
    pub max_depth_m: f32,
    pub altitude_smoothing_sigma: f32,
}

impl Default for IsostasyParams {
    fn default() -> Self {
        Self {
            sea_level_fraction: 0.4,
            max_elevation_m: 4000.0,
            max_depth_m: 500.0,
            altitude_smoothing_sigma: 2.0,
        }
    }
}

/// Cache for the isostasy result — recomputed when the source changes
/// or the sea level slider moves.
#[derive(Resource, Default)]
pub struct IsostasyCache {
    /// The computed land ratio, peak altitude, depth.
    pub land_ratio: f32,
    pub peak_altitude_m: f32,
    pub max_depth_m: f32,
    pub sea_level_normalized: f32,
    /// Sea level that was used for the last computation.
    pub computed_sea_level: f32,
    /// Whether valid data is available.
    pub valid: bool,
    /// The computed isostasy heightmap (stored for rendering and FBM input).
    pub heightmap: Option<GridF32>,
}

// ── UI Actions ──────────────────────────────────────────────────────────

/// Flags set by UI draw functions, consumed by Bevy systems.
#[derive(Resource, Default)]
pub struct UiActions {
    /// Set to true when the user clicks "Export".
    pub export_requested: bool,
    /// Set to a directory path when the user clicks "Load".
    pub load_requested: Option<std::path::PathBuf>,
    /// Feedback message displayed temporarily after export/load.
    pub last_message: Option<(String, std::time::Instant, bool)>,
    /// Cached list of export directories (invalidated after export).
    pub cached_dirs: Option<Vec<std::path::PathBuf>>,
    /// Set to true when the user clicks "Step" (single timestep).
    pub step_requested: bool,
    /// Set to true when the user clicks "Center Map".
    pub center_requested: bool,
    pub center_offset_changed: bool,
    /// Top-bar run button dispatches — consumed by parameter_panel / bridge systems.
    pub run_solver_requested: bool,
    pub run_fbm_requested: bool,
    pub run_erosion_requested: bool,
    pub run_hydrology_requested: bool,
}

// ── Centering ───────────────────────────────────────────────────────────

#[derive(Resource, Default)]
pub struct CenteringState {
    /// The original (unshifted) field, saved on first center operation.
    pub original_field: Option<ymir_core::tectonics::solver::field::Field2D>,
    pub original_plate_ids: Option<Vec<usize>>,
    pub original_plates: Option<Vec<ymir_core::tectonics::plates::Plate>>,
    pub original_grid_size: usize,
    /// Auto-centering shift (from circular mean).
    pub auto_shift: (i32, i32),
    /// Manual offset on top of auto shift.
    pub offset_x: i32,
    pub offset_y: i32,
}

// ── Dynamic plate boundaries ─────────────────────────────────────────────

/// Live state of plate boundaries during a dynamic simulation.
/// Updated from solver snapshots when `dynamic_boundaries` is enabled.
#[derive(Resource, Default)]
pub struct DynamicPlateIds {
    /// Current plate_ids (updated by solver snapshots).
    pub ids: Option<Vec<usize>>,
    /// Current plate seed positions and active flags.
    pub plates: Option<Vec<ymir_core::tectonics::plates::Plate>>,
    /// Number of active plates remaining.
    pub active_count: usize,
    /// Grid size (for indexing into ids).
    pub grid_size: usize,
}

// ── FBM Upscaling ───────────────────────────────────────────────────────

/// Current state of the FBM upscale computation.
#[derive(Default, Debug, Clone)]
pub enum FbmState {
    #[default]
    Idle,
    Running,
    Completed {
        elapsed: std::time::Duration,
    },
}

/// Cache for the FBM upscale result.
#[derive(Resource, Default)]
pub struct UpscaleCache {
    /// The upscaled heightmap (high resolution).
    pub heightmap: Option<GridF32>,
    /// The slope field at target resolution.
    pub slope: Option<GridF32>,
    /// Current computation state.
    pub state: FbmState,
    /// Pending command data — set by the UI, consumed by the dispatch system.
    pub pending_config: Option<ymir_core::terrain::upscale::FbmUpscaleConfig>,
    pub pending_seed: Option<ymir_core::seed::WorldSeed>,
    pub pending_sea_level: Option<f32>,
}

#[derive(Resource, Clone, Debug)]
pub struct FbmParams {
    pub target_size: usize,
    pub octaves: usize,
    pub lacunarity: f64,
    pub persistence: f64,
    pub amplitude_base: f64,
    pub amplitude_slope_factor: f64,
    pub max_anisotropy: f64,
    pub submarine_damping: f64,
    pub domain_warp_enabled: bool,
    pub domain_warp_strength: f64,
    pub domain_warp_frequency: f64,
    pub domain_warp_octaves: usize,
}

impl Default for FbmParams {
    fn default() -> Self {
        Self {
            target_size: 1024,
            octaves: 7,
            lacunarity: 2.0,
            persistence: 0.5,
            amplitude_base: 0.08,
            amplitude_slope_factor: 3.0,
            max_anisotropy: 3.0,
            submarine_damping: 0.3,
            domain_warp_enabled: false,
            domain_warp_strength: 0.4,
            domain_warp_frequency: 0.5,
            domain_warp_octaves: 3,
        }
    }
}

// ── Erosion ─────────────────────────────────────────────────────────────

#[derive(Default, Debug, Clone)]
pub enum ErosionState {
    #[default]
    Idle,
    Running {
        completed: usize,
        total: usize,
    },
    Completed {
        elapsed: std::time::Duration,
    },
}

#[derive(Resource, Default)]
pub struct ErosionCache {
    pub heightmap: Option<GridF32>,
    pub sediment: Option<GridF32>,
    pub stats: Option<ErosionStats>,
    pub state: ErosionState,
    pub pending_config: Option<ErosionConfig>,
    pub pending_seed: Option<ymir_core::seed::WorldSeed>,
}

// ── Flow / Rivers ───────────────────────────────────────────────────────

#[derive(Default, Debug, Clone)]
pub enum FlowState {
    #[default]
    Idle,
    Running,
    Completed {
        elapsed: std::time::Duration,
    },
}

#[derive(Resource, Default)]
pub struct FlowCache {
    pub result: Option<FlowResult>,
    pub state: FlowState,
    pub pending_config: Option<FlowConfig>,
    pub rivers: Option<RiverNetwork>,
    pub river_config: RiverConfig,
    pub rivers_dirty: bool,
}

// ── Run timer ────────────────────────────────────────────────────────────

/// Tracks the start time of the current pipeline operation for the status bar timer.
#[derive(Resource, Default)]
pub struct RunTimer {
    pub started_at: Option<std::time::Instant>,
    /// Elapsed time of the last completed operation (shown until next run).
    pub last_elapsed: Option<std::time::Duration>,
}

// ── Toasts ───────────────────────────────────────────────────────────────

#[derive(Resource, Default)]
pub struct Toasts {
    pub messages: Vec<(String, std::time::Instant, bool)>,
}

impl Toasts {
    pub fn add(&mut self, msg: String, success: bool) {
        self.messages.push((msg, std::time::Instant::now(), success));
    }

    pub fn prune(&mut self) {
        let now = std::time::Instant::now();
        self.messages.retain(|(_, t, _)| now.duration_since(*t).as_secs() < 5);
    }
}

// ── Lakes ────────────────────────────────────────────────────────────────

#[derive(Resource, Default)]
pub struct LakeCache {
    pub result: Option<LakeResult>,
    pub config: LakeConfig,
    pub dirty: bool,
}

// ── Cursor ───────────────────────────────────────────────────────────────

#[derive(Resource, Default)]
pub struct CursorWorldPos {
    pub pos: Option<Vec2>,
}
