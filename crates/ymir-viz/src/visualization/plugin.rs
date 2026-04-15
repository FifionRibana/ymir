//! Visualization plugin: terrain sprite + texture update.

use bevy::prelude::*;

use super::erosion::update_erosion_texture;
use super::isostasy::{recompute_isostasy_cache, render_isostasy_texture};
use super::render::{TerrainDisplay, setup_solver_terrain_sprite, update_terrain_texture};
use super::rivers::render_river_overlay;
use super::upscale::update_upscale_texture;
use crate::state::{PipelinePhase, ViewMode};

pub struct SolverVisualizationPlugin;

impl Plugin for SolverVisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<TerrainDisplay>();
        app.add_systems(Startup, setup_solver_terrain_sprite);

        // Terrain texture (thickness view): always needed
        app.add_systems(Update, update_terrain_texture);

        // Isostasy cache recomputation: runs in ALL phases
        app.add_systems(Update, recompute_isostasy_cache);

        // Isostasy texture rendering: only when Altitude mode AND
        // Tectonics or Isostasy phase is selected
        app.add_systems(
            Update,
            render_isostasy_texture.run_if(in_state(ViewMode::Altitude)).run_if(
                |phase: Res<State<PipelinePhase>>| {
                    matches!(phase.get(), PipelinePhase::Tectonics | PipelinePhase::Isostasy)
                },
            ),
        );

        // Upscale texture rendering: only when UpscaleFbm phase + Altitude mode
        app.add_systems(
            Update,
            update_upscale_texture
                .run_if(in_state(PipelinePhase::UpscaleFbm))
                .run_if(in_state(ViewMode::Altitude)),
        );

        // Erosion texture rendering: only when Erosion phase + Altitude mode
        app.add_systems(
            Update,
            update_erosion_texture
                .run_if(in_state(PipelinePhase::Erosion))
                .run_if(in_state(ViewMode::Altitude)),
        );

        // River overlay: runs after terrain textures, blends onto existing image
        app.add_systems(
            Update,
            render_river_overlay
                .run_if(in_state(ViewMode::Altitude))
                .after(render_isostasy_texture)
                .after(update_upscale_texture)
                .after(update_erosion_texture),
        );
    }
}
