pub mod parameter_panel;
pub mod pipeline_panel;
pub mod statistics_panel;

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPlugin, EguiPrimaryContextPass, egui};

use crate::state::{GenerationParamsUi, PipelinePhase, ViewMode, ViewState};

pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin::default()).add_systems(
            EguiPrimaryContextPass,
            (configure_egui_style, ui_top_bar, ui_right_panel, ui_bottom_bar),
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
    view_mode: Res<State<ViewMode>>,
    mut next_view_mode: ResMut<NextState<ViewMode>>,
    gen_params: Res<GenerationParamsUi>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            for &mode in ViewMode::ALL {
                let active = *view_mode.get() == mode;
                let enabled = mode.is_enabled();
                let response =
                    ui.add_enabled(enabled, egui::Button::new(mode.label()).selected(active));
                if enabled && response.clicked() {
                    next_view_mode.set(mode);
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

            let mut plates = view_state.overlays.plates;
            if ui.checkbox(&mut plates, "Plates").changed() {
                view_state.overlays.plates = plates;
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.monospace(format!("Seed: {}", gen_params.seed));
            });
        });
    });
}

#[derive(bevy::ecs::system::SystemParam)]
struct UiRightPanelParams<'w> {
    erosion: ResMut<'w, crate::state::ErosionParams>,
    climate: ResMut<'w, crate::state::ClimateParams>,
    gen_params: ResMut<'w, GenerationParamsUi>,
    stats: Res<'w, crate::state::TerrainStats>,
    tectonic_state: Option<ResMut<'w, crate::state::TectonicState>>,
    solver_config: ResMut<'w, crate::state::SolverConfig>,
    bridge: ResMut<'w, crate::bridge::SolverBridge>,
    isostasy_params: ResMut<'w, crate::state::IsostasyParams>,
    isostasy_cache: Res<'w, crate::state::IsostasyCache>,
    fbm_params: ResMut<'w, crate::state::FbmParams>,
    upscale_cache: ResMut<'w, crate::state::UpscaleCache>,
    ui_actions: ResMut<'w, crate::state::UiActions>,
    terrain_display: Res<'w, crate::visualization::render::TerrainDisplay>,
}

fn ui_right_panel(
    mut contexts: EguiContexts,
    current_phase: Res<State<PipelinePhase>>,
    mut next_phase: ResMut<NextState<PipelinePhase>>,
    pipeline_state: Res<crate::state::PipelineState>,
    mut params: UiRightPanelParams,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::SidePanel::right("right_panel").exact_width(260.0).show(ctx, |ui| {
        egui::ScrollArea::vertical().show(ui, |ui| {
            let has_terrain = params.terrain_display.s_field.is_some();
            pipeline_panel::draw(
                ui,
                current_phase.get(),
                &mut next_phase,
                &pipeline_state,
                &mut params.ui_actions,
                has_terrain,
            );
            ui.separator();
            parameter_panel::draw(
                ui,
                current_phase.get(),
                &mut params.erosion,
                &mut params.climate,
                &mut params.gen_params,
                params.tectonic_state,
                &mut params.solver_config,
                &mut params.bridge,
                &mut params.isostasy_params,
                &params.isostasy_cache,
                &mut params.fbm_params,
                &mut params.upscale_cache,
            );
            ui.separator();
            statistics_panel::draw(ui, &params.stats);
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
