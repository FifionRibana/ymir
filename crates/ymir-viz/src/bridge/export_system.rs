//! Bevy systems for pipeline export and load, triggered by UiActions flags.

use std::path::Path;

use bevy::prelude::*;

use ymir_core::export::PipelineExport;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
use ymir_core::tectonics::plates::generate_plates;
use ymir_core::tectonics::solver::field::Field2D;

use crate::state::{IsostasyCache, IsostasyParams, TectonicState, UiActions};
use crate::visualization::render::TerrainDisplay;

pub fn handle_export(
    mut ui_actions: ResMut<UiActions>,
    tectonic_state: Option<Res<TectonicState>>,
    isostasy_params: Res<IsostasyParams>,
    terrain_display: Res<TerrainDisplay>,
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

    let dir_name =
        export.dir.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    ui_actions.last_message =
        Some((format!("Exported to {dir_name}"), std::time::Instant::now(), true));
    // Invalidate cached dir list so Load section sees the new export
    ui_actions.cached_dirs = None;

    info!("Exported to {}", export.dir.display());
}

pub fn handle_load(
    mut ui_actions: ResMut<UiActions>,
    mut terrain_display: ResMut<TerrainDisplay>,
    mut isostasy_cache: ResMut<IsostasyCache>,
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

    // Load thickness → TerrainDisplay
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

    // Restore TectonicState from metadata if not present
    if tectonic_state.is_none() {
        let seed = export.metadata.seed;
        let config = export.metadata.plates.clone();
        let init = generate_plates(&config, &WorldSeed::new(seed));
        commands.insert_resource(TectonicState { init, config, seed, dirty: true, generation: 1 });
    }

    let dir_name = dir.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    ui_actions.last_message = Some((format!("Loaded {dir_name}"), std::time::Instant::now(), true));
    ui_actions.cached_dirs = None;
}
