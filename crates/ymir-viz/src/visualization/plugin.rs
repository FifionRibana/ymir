//! Visualization plugin: terrain sprite + texture update.

use bevy::prelude::*;

use super::render::{setup_solver_terrain_sprite, update_terrain_texture, TerrainDisplay};

pub struct SolverVisualizationPlugin;

impl Plugin for SolverVisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<TerrainDisplay>();
        app.add_systems(Startup, setup_solver_terrain_sprite);
        app.add_systems(Update, update_terrain_texture);
    }
}
