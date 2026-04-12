use bevy::prelude::*;

mod bridge;
mod camera;
mod cursor_inspector;
mod state;
mod tectonic_view;
mod terrain_view;
mod ui;
mod visualization;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Ymir — Continent Generator".to_string(),
                resolution: (1280, 720).into(),
                ..default()
            }),
            ..default()
        }))
        .init_resource::<state::ViewState>()
        .init_resource::<state::PipelineState>()
        .init_resource::<state::ErosionParams>()
        .init_resource::<state::SolverConfig>()
        .init_resource::<state::ClimateParams>()
        .init_resource::<state::GenerationParamsUi>()
        .init_resource::<state::CursorWorldPos>()
        .add_plugins((
            camera::CameraPlugin,
            terrain_view::TerrainViewPlugin,
            tectonic_view::TectonicViewPlugin,
            ui::UiPlugin,
            cursor_inspector::CursorInspectorPlugin,
            bridge::TectonicsBridgePlugin,
            visualization::SolverVisualizationPlugin,
        ))
        .run();
}
