use bevy::prelude::*;

mod camera;
mod cursor_inspector;
mod state;
mod terrain_view;
mod ui;

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
        .init_resource::<state::TectonicsParams>()
        .init_resource::<state::ClimateParams>()
        .init_resource::<state::GenerationParamsUi>()
        .init_resource::<state::CursorWorldPos>()
        .add_plugins((
            camera::CameraPlugin,
            terrain_view::TerrainViewPlugin,
            ui::UiPlugin,
            cursor_inspector::CursorInspectorPlugin,
        ))
        .run();
}
