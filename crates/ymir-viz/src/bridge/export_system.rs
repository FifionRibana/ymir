//! Bevy systems for pipeline export and load, triggered by UiActions flags.

use std::path::Path;

use bevy::prelude::*;

use ymir_core::export::PipelineExport;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
use ymir_core::tectonics::plates::generate_plates;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::terrain::upscale::FbmUpscaleConfig;

use crate::state::{
    FbmParams, FbmState, GenerationParamsUi, IsostasyCache, IsostasyParams, TectonicState,
    UiActions, UpscaleCache,
};
use crate::visualization::render::TerrainDisplay;

pub fn handle_export(
    mut ui_actions: ResMut<UiActions>,
    tectonic_state: Option<Res<TectonicState>>,
    isostasy_params: Res<IsostasyParams>,
    terrain_display: Res<TerrainDisplay>,
    upscale_cache: Res<UpscaleCache>,
    fbm_params: Res<FbmParams>,
) {
    if !ui_actions.export_requested {
        return;
    }
    ui_actions.export_requested = false;

    let Some(tecto) = tectonic_state else {
        ui_actions.last_message =
            Some(("No tectonic state to export".into(), std::time::Instant::now(), false));
        return;
    };

    let Some(ref s_field) = terrain_display.s_field else {
        ui_actions.last_message =
            Some(("No solver result to export".into(), std::time::Instant::now(), false));
        return;
    };

    let output_root = Path::new("output");
    let grid_size = tecto.init.grid_size;
    let seed = tecto.seed;

    let mut export = PipelineExport::new(output_root, seed, grid_size, &tecto.config);

    // Save thickness (Field2D → GridF32)
    let n = s_field.n();
    let data: Vec<f32> = s_field.data().iter().map(|&v| v as f32).collect();
    let thickness = GridF32::from_vec(n, n, data);
    if let Err(e) = export.save_thickness(&thickness) {
        ui_actions.last_message =
            Some((format!("Export failed: {e}"), std::time::Instant::now(), false));
        return;
    }

    // Compute and save altitude
    let iso_config = IsostasyConfig {
        sea_level_fraction: isostasy_params.sea_level_fraction,
        max_elevation_m: isostasy_params.max_elevation_m,
        max_depth_m: isostasy_params.max_depth_m,
        altitude_smoothing_sigma: isostasy_params.altitude_smoothing_sigma,
        ..Default::default()
    };
    let iso_result = compute_isostasy(s_field, &iso_config);
    if let Err(e) = export.save_altitude(&iso_result, &iso_config) {
        ui_actions.last_message =
            Some((format!("Export failed: {e}"), std::time::Instant::now(), false));
        return;
    }

    // Save upscaled heightmap if available
    if let Some(ref heightmap) = upscale_cache.heightmap {
        let fbm_config = FbmUpscaleConfig {
            target_size: fbm_params.target_size,
            octaves: fbm_params.octaves,
            lacunarity: fbm_params.lacunarity,
            persistence: fbm_params.persistence,
            amplitude_base: fbm_params.amplitude_base,
            amplitude_slope_factor: fbm_params.amplitude_slope_factor,
            max_anisotropy: fbm_params.max_anisotropy,
            submarine_damping: fbm_params.submarine_damping,
            base_frequency: 1.0,
        };
        if let Err(e) = export.save_upscaled(heightmap, &fbm_config) {
            ui_actions.last_message =
                Some((format!("Export failed: {e}"), std::time::Instant::now(), false));
            return;
        }
    }

    let dir_name =
        export.dir.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    ui_actions.last_message =
        Some((format!("Exported to {dir_name}"), std::time::Instant::now(), true));
    ui_actions.cached_dirs = None;

    info!("Exported to {}", export.dir.display());
}

#[allow(clippy::too_many_arguments)]
pub fn handle_load(
    mut ui_actions: ResMut<UiActions>,
    mut terrain_display: ResMut<TerrainDisplay>,
    mut isostasy_cache: ResMut<IsostasyCache>,
    mut upscale_cache: ResMut<UpscaleCache>,
    mut gen_params: ResMut<GenerationParamsUi>,
    mut isostasy_params: ResMut<IsostasyParams>,
    mut fbm_params: ResMut<FbmParams>,
    tectonic_state: Option<ResMut<TectonicState>>,
    mut commands: Commands,
) {
    let Some(dir) = ui_actions.load_requested.take() else {
        return;
    };

    let export = match PipelineExport::load(&dir) {
        Ok(e) => e,
        Err(e) => {
            ui_actions.last_message =
                Some((format!("Load failed: {e}"), std::time::Instant::now(), false));
            return;
        }
    };

    let meta = &export.metadata;

    // ── Restore generation params ────────────────────────────────────────
    gen_params.seed = meta.seed;
    gen_params.meters_per_pixel = meta.meters_per_pixel;

    // ── Restore isostasy params ──────────────────────────────────────────
    if let Some(ref iso) = meta.isostasy {
        isostasy_params.sea_level_fraction = iso.config.sea_level_fraction;
        isostasy_params.max_elevation_m = iso.config.max_elevation_m;
        isostasy_params.max_depth_m = iso.config.max_depth_m;
        isostasy_params.altitude_smoothing_sigma = iso.config.altitude_smoothing_sigma;
    }

    // ── Restore FBM params ───────────────────────────────────────────────
    if let Some(ref upscale) = meta.upscale {
        fbm_params.target_size = upscale.target_size;
        fbm_params.octaves = upscale.octaves;
        fbm_params.lacunarity = upscale.lacunarity;
        fbm_params.persistence = upscale.persistence;
        fbm_params.amplitude_base = upscale.amplitude_base;
        fbm_params.amplitude_slope_factor = upscale.amplitude_slope_factor;
        fbm_params.max_anisotropy = upscale.max_anisotropy;
        fbm_params.submarine_damping = upscale.submarine_damping;
    }

    // ── Load thickness → TerrainDisplay ──────────────────────────────────
    match export.load_thickness() {
        Ok(thickness_grid) => {
            let n = thickness_grid.width;
            let mut field = Field2D::new(n);
            for j in 0..n {
                for i in 0..n {
                    field.set(i, j, thickness_grid.data[j * n + i] as f64);
                }
            }
            terrain_display.update_field(field);
            info!("Loaded thickness {}x{} from {}", n, n, dir.display());
        }
        Err(e) => warn!("No thickness to load: {e}"),
    }

    // Invalidate isostasy cache so it recomputes from the loaded field
    isostasy_cache.valid = false;

    // ── Load upscaled heightmap ──────────────────────────────────────────
    match export.load_upscaled() {
        Ok(heightmap) => {
            info!(
                "Loaded upscaled {}x{} from {}",
                heightmap.width,
                heightmap.height,
                dir.display()
            );
            upscale_cache.heightmap = Some(heightmap);
            upscale_cache.slope = None;
            upscale_cache.state = FbmState::Completed { elapsed: std::time::Duration::ZERO };
        }
        Err(_) => {
            upscale_cache.heightmap = None;
            upscale_cache.slope = None;
            upscale_cache.state = FbmState::Idle;
        }
    }

    // ── Restore TectonicState (always, to sync seed and plate config) ───
    let seed = meta.seed;
    let config = meta.plates.clone();
    let init = generate_plates(&config, &WorldSeed::new(seed));
    if let Some(mut state) = tectonic_state {
        state.init = init;
        state.config = config;
        state.seed = seed;
        state.dirty = true;
        state.generation = state.generation.wrapping_add(1);
    } else {
        commands
            .insert_resource(TectonicState { init, config, seed, dirty: true, generation: 1 });
    }

    let dir_name = dir.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    ui_actions.last_message = Some((format!("Loaded {dir_name}"), std::time::Instant::now(), true));
    ui_actions.cached_dirs = None;
}
