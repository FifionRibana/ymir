//! Ymir interactive visualisation binary — single-engine (C1) after the v2
//! legacy sunset. The Stokes-coupled `tectonics_v2` bridge, its visualization
//! sprite, its UI panels, and the per-phase legacy pipeline (isostasy / FBM /
//! erosion / hydrology + the climate/biome stubs) have all been removed. This
//! binary wires the C1 bridge + C1 visualization + the egui UI into a vanilla
//! Bevy app.

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
                    title: "Ymir — Continent Generator".to_string(),
                    resolution: (1280, 720).into(),
                    ..default()
                }),
                ..default()
            })
            .set(bevy::log::LogPlugin {
                filter: "warn,ymir_core=info,ymir_viz=info".to_string(),
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
        ui::UiPlugin,
        bridge::c1::C1BridgePlugin,
        visualization::C1VisualizationPlugin,
    ));

    app.run();
}
