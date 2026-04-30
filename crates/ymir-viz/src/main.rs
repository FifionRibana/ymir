//! Ymir interactive visualisation binary — `tectonics_v2` only after
//! Step 8.6 Phase 8h sunset. The legacy `tectonics::` bridge and
//! every plugin/panel that drove its pipeline phases have been
//! removed; this binary now wires the v2 bridge plugin, the v2
//! visualization plugin, and the v2 UI plugin into a vanilla Bevy
//! app.

use bevy::prelude::*;

mod bridge;
mod camera;
mod ui;
mod visualization;

fn main() {
    let mut app = App::new();
    app.add_plugins(
        DefaultPlugins
            .set(WindowPlugin {
                primary_window: Some(Window {
                    title: "Ymir — Continent Generator (v2)".to_string(),
                    resolution: (1280, 720).into(),
                    ..default()
                }),
                ..default()
            })
            .set(bevy::log::LogPlugin {
                filter: "warn,ymir_core::tectonics_v2=info,ymir_viz=info".to_string(),
                level: bevy::log::Level::DEBUG,
                custom_layer: |_app| {
                    let file = std::fs::OpenOptions::new()
                        .create(true)
                        .write(true)
                        .truncate(true)
                        .open("ymir.log")
                        .expect("Failed to open ymir.log");

                    Some(Box::new(
                        tracing_subscriber::fmt::layer()
                            .with_writer(std::sync::Mutex::new(file))
                            .with_ansi(false),
                    ))
                },
                ..default()
            }),
    )
    .add_plugins((
        camera::CameraPlugin,
        bridge::v2::V2BridgePlugin,
        visualization::V2VisualizationPlugin,
        ui::UiPlugin,
    ));

    app.run();
}
