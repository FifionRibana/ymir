use bevy::prelude::*;

mod bridge;
mod camera;
mod cursor_inspector;
mod state;
mod tectonic_view;
mod terrain_view;
mod ui;
mod visualization;

/// Selects which bridge plugin (legacy `tectonics::` vs new
/// `tectonics_v2::`) is wired into the app at startup. Step 8.6
/// ships both side-by-side; the legacy path is removed in Phase 8
/// (sunset commit).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BridgeMode {
    Legacy,
    V2,
}

impl BridgeMode {
    /// Resolved from `YMIR_BRIDGE` env var. Default is `V2` per the
    /// Step 8.6 issue (D1 sunset path).
    fn from_env() -> Self {
        match std::env::var("YMIR_BRIDGE").ok().as_deref() {
            Some("legacy") => BridgeMode::Legacy,
            Some("v2") | None => BridgeMode::V2,
            Some(other) => {
                eprintln!(
                    "[ymir-viz] unknown YMIR_BRIDGE='{}', falling back to v2 (valid: legacy|v2)",
                    other
                );
                BridgeMode::V2
            }
        }
    }
}

fn main() {
    let bridge_mode = BridgeMode::from_env();
    eprintln!("[ymir-viz] bridge mode: {:?}", bridge_mode);

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
                filter: "warn,ymir_core::tectonics=debug".to_string(),
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
    .init_resource::<state::ViewState>()
    .init_state::<state::ViewMode>()
    .init_state::<state::PipelinePhase>()
    .init_resource::<state::PipelineState>()
    .init_resource::<state::ErosionParams>()
    .init_resource::<state::SolverConfig>()
    .init_resource::<state::ClimateParams>()
    .init_resource::<state::GenerationParamsUi>()
    .init_resource::<state::IsostasyParams>()
    .init_resource::<state::IsostasyCache>()
    .init_resource::<state::FbmParams>()
    .init_resource::<state::UpscaleCache>()
    .init_resource::<state::ErosionCache>()
    .init_resource::<state::FlowCache>()
    .init_resource::<state::LakeCache>()
    .init_resource::<state::CenteringState>()
    .init_resource::<state::CursorWorldPos>()
    .init_resource::<state::Toasts>()
    .init_resource::<state::RunTimer>()
    .init_resource::<state::DynamicPlateIds>()
    .init_resource::<state::GridUiState>()
    .add_plugins((
        camera::CameraPlugin,
        terrain_view::TerrainViewPlugin,
        tectonic_view::TectonicViewPlugin,
        ui::UiPlugin,
        cursor_inspector::CursorInspectorPlugin,
        visualization::SolverVisualizationPlugin,
    ));

    // Step 8.6 D1 — bridge plugin dispatch. Both bridges coexist
    // through Phase 7; the sunset commit (Phase 8) removes the
    // legacy arm. Phase 5 adds the v2 visualization plugin so the
    // selected raster field actually paints to a sprite when in
    // v2 mode.
    match bridge_mode {
        BridgeMode::Legacy => {
            app.add_plugins(bridge::TectonicsBridgePlugin);
        }
        BridgeMode::V2 => {
            app.add_plugins(bridge::v2::V2BridgePlugin);
            app.add_plugins(visualization::V2VisualizationPlugin);
        }
    }

    app.run();
}
