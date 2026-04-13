pub mod parameter_panel;
pub mod pipeline_panel;
pub mod statistics_panel;

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPlugin, EguiPrimaryContextPass, egui};

use crate::state::{GenerationParamsUi, ViewMode, ViewState};

pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin::default())
            .add_systems(
                EguiPrimaryContextPass,
                (configure_egui_style, ui_top_bar, ui_right_panel, ui_bottom_bar).chain(),
            );
    }
}

fn configure_egui_style(mut contexts: EguiContexts, mut done: Local<bool>) {
    if *done {
        return;
    }
    let Ok(ctx) = contexts.ctx_mut() else { return };
    *done = true;

    let mut style = (*ctx.style()).clone();
    let copper = egui::Color32::from_rgb(0xB8, 0x73, 0x33);
    let bronze = egui::Color32::from_rgb(0xA0, 0x72, 0x4A);

    style.visuals.widgets.active.bg_fill = copper;
    style.visuals.widgets.hovered.bg_fill = bronze;
    style.visuals.selection.bg_fill = copper.linear_multiply(0.4);
    style.visuals.selection.stroke = egui::Stroke::new(1.0, copper);
    style.visuals.widgets.active.fg_stroke = egui::Stroke::new(1.0, copper);

    ctx.set_style(style);
}

fn ui_top_bar(
    mut contexts: EguiContexts,
    mut view_state: ResMut<ViewState>,
    gen_params: Res<GenerationParamsUi>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            for &mode in ViewMode::ALL {
                let active = view_state.mode == mode;
                let enabled = mode.is_enabled();
                let response =
                    ui.add_enabled(enabled, egui::Button::new(mode.label()).selected(active));
                if enabled && response.clicked() {
                    view_state.mode = mode;
                }
            }

            ui.separator();

            let mut hillshade = view_state.overlays.hillshade;
            if ui.checkbox(&mut hillshade, "Hillshade").changed() {
                view_state.overlays.hillshade = hillshade;
            }

            ui.add_enabled(false, egui::Checkbox::new(&mut false, "Rivers"));

            let mut grid = view_state.overlays.grid;
            if ui.checkbox(&mut grid, "Grid").changed() {
                view_state.overlays.grid = grid;
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.monospace(format!("Seed: {}", gen_params.seed));
            });
        });
    });
}

#[allow(clippy::too_many_arguments)]
fn ui_right_panel(
    mut contexts: EguiContexts,
    mut view_state: ResMut<ViewState>,
    pipeline_state: Res<crate::state::PipelineState>,
    mut erosion: ResMut<crate::state::ErosionParams>,
    mut climate: ResMut<crate::state::ClimateParams>,
    mut gen_params: ResMut<GenerationParamsUi>,
    stats: Res<crate::state::TerrainStats>,
    tectonic_state: Option<ResMut<crate::state::TectonicState>>,
    mut solver_config: ResMut<crate::state::SolverConfig>,
    mut bridge: ResMut<crate::bridge::SolverBridge>,
    mut isostasy_params: ResMut<crate::state::IsostasyParams>,
    isostasy_cache: Res<crate::state::IsostasyCache>,
    mut ui_actions: ResMut<crate::state::UiActions>,
    terrain_display: Res<crate::visualization::render::TerrainDisplay>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::SidePanel::right("right_panel")
        .exact_width(260.0)
        .show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                let has_terrain = terrain_display.s_field.is_some();
                pipeline_panel::draw(
                    ui,
                    &mut view_state,
                    &pipeline_state,
                    &mut ui_actions,
                    has_terrain,
                );
                ui.separator();
                parameter_panel::draw(
                    ui,
                    &view_state,
                    &mut erosion,
                    &mut climate,
                    &mut gen_params,
                    tectonic_state,
                    &mut solver_config,
                    &mut bridge,
                    &mut isostasy_params,
                    &isostasy_cache,
                );
                ui.separator();
                statistics_panel::draw(ui, &stats);
            });
        });
}

fn ui_bottom_bar(mut contexts: EguiContexts, time: Res<Time>) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            ui.label(
                egui::RichText::new("LMB: Pan  Scroll: Zoom").small().color(egui::Color32::GRAY),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let fps = 1.0 / time.delta_secs().max(0.001);
                ui.monospace(format!("{fps:.0} FPS"));
            });
        });
    });
}
