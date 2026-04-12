use bevy::prelude::*;
use ymir_core::grid::GridF32;
use ymir_core::tectonics::plates::{PlateConfig, PlateInitResult};
use ymir_core::tectonics::solver::config::{NonlinearSolver, Preconditioner};

// ── View ─────────────────────────────────────────────────────────────────

#[derive(PartialEq, Eq, Clone, Copy, Default, Debug)]
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
    pub hillshade: bool,
    pub grid: bool,
}

impl Default for OverlayFlags {
    fn default() -> Self {
        Self { rivers: false, hillshade: true, grid: false }
    }
}

#[derive(Resource, Clone, Debug)]
pub struct ViewState {
    pub mode: ViewMode,
    pub overlays: OverlayFlags,
    pub selected_phase: PipelinePhase,
}

impl Default for ViewState {
    fn default() -> Self {
        Self {
            mode: ViewMode::Altitude,
            overlays: OverlayFlags::default(),
            selected_phase: PipelinePhase::Erosion,
        }
    }
}

// ── Pipeline ─────────────────────────────────────────────────────────────

#[derive(PartialEq, Eq, Clone, Copy, Debug, Default)]
pub enum PipelinePhase {
    Tectonics,
    Isostasy,
    UpscaleFbm,
    #[default]
    Erosion,
    Climate,
    Lakes,
}

impl PipelinePhase {
    pub const ALL: &[PipelinePhase] = &[
        PipelinePhase::Tectonics,
        PipelinePhase::Isostasy,
        PipelinePhase::UpscaleFbm,
        PipelinePhase::Erosion,
        PipelinePhase::Climate,
        PipelinePhase::Lakes,
    ];

    pub fn label(self) -> &'static str {
        match self {
            PipelinePhase::Tectonics => "Tectonics",
            PipelinePhase::Isostasy => "Isostasy",
            PipelinePhase::UpscaleFbm => "Upscale + FBM",
            PipelinePhase::Erosion => "Erosion",
            PipelinePhase::Climate => "Climate",
            PipelinePhase::Lakes => "Lakes",
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
            min_slope: 0.01,
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
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            num_timesteps: 300,
            gravity_factor: 1.0,
            cfl_factor: 0.5,
            power_law_n: 3.0,
            picard_relaxation: 0.7,
            nonlinear_solver: NonlinearSolver::Picard,
            continuation_enabled: true,
            strain_rate_min: 1e-3,
            eta_max: 1e4,
            preconditioner: Preconditioner::default(),
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

// ── Cursor ───────────────────────────────────────────────────────────────

#[derive(Resource, Default)]
pub struct CursorWorldPos {
    pub pos: Option<Vec2>,
}
