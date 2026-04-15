pub mod left_toolbar;
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
            (configure_egui_style, ui_top_bar, ui_right_panel, ui_bottom_bar, ui_toasts),
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
    current_phase: Res<State<PipelinePhase>>,
    bridge: Res<crate::bridge::SolverBridge>,
    upscale_cache: Res<crate::state::UpscaleCache>,
    erosion_cache: Res<crate::state::ErosionCache>,
    flow_cache: Res<crate::state::FlowCache>,
    isostasy_cache: Res<crate::state::IsostasyCache>,
    mut ui_actions: ResMut<crate::state::UiActions>,
    mut run_timer: ResMut<crate::state::RunTimer>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            // ── Left side: View modes + Overlays ──
            for &mode in ViewMode::ALL {
                let active = *view_mode.get() == mode;
                let enabled = mode.is_enabled();
                let response =
                    ui.add_enabled(enabled, egui::Button::new(mode.label()).selected(active));
                if enabled && response.clicked() {
                    next_view_mode.set(mode);
                }
            }

            ui.menu_button("\u{1f53d} Layers", |ui| {
                ui.checkbox(&mut view_state.overlays.hillshade, "Hillshade");
                ui.checkbox(&mut view_state.overlays.rivers, "Rivers");
                ui.checkbox(&mut view_state.overlays.lakes, "Lakes");
                ui.checkbox(&mut view_state.overlays.grid, "Grid");
                ui.checkbox(&mut view_state.overlays.plates, "Plates");
            });

            // ── Right side: Seed + Run / Step / Cancel + Progress ──
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let phase = *current_phase.get();
                let (run_label, run_enabled, is_running) = run_button_state(
                    phase,
                    &bridge,
                    &isostasy_cache,
                    &upscale_cache,
                    &erosion_cache,
                );

                if ui.add_enabled(run_enabled, egui::Button::new(run_label)).clicked() {
                    dispatch_run(phase, &mut ui_actions, &mut run_timer);
                }

                if phase == PipelinePhase::Tectonics {
                    if ui.add_enabled(!is_running, egui::Button::new("⏸ Step")).clicked() {
                        ui_actions.step_requested = true;
                    }
                }

                if is_running {
                    if ui.button("⏹ Cancel").clicked() {
                        bridge.cancel_flag.store(true, std::sync::atomic::Ordering::Relaxed);
                    }
                }

                ui.separator();
                ui.monospace(format!("Seed: {}", gen_params.seed));
            });
        });
    });
}

fn run_button_state(
    phase: PipelinePhase,
    bridge: &crate::bridge::SolverBridge,
    isostasy_cache: &crate::state::IsostasyCache,
    upscale_cache: &crate::state::UpscaleCache,
    erosion_cache: &crate::state::ErosionCache,
) -> (&'static str, bool, bool) {
    use crate::bridge::plugin::SolverState;
    use crate::state::{ErosionState, FbmState};

    let solver_running = matches!(bridge.state, SolverState::Running { .. });
    let erosion_running = matches!(erosion_cache.state, ErosionState::Running { .. });
    let fbm_running = matches!(upscale_cache.state, FbmState::Running);
    let any_running = solver_running || erosion_running || fbm_running;

    match phase {
        PipelinePhase::Tectonics => ("▶ Run Tectonics", !any_running, solver_running),
        PipelinePhase::Isostasy => ("⊕ Center Map", !any_running, false),
        PipelinePhase::UpscaleFbm => {
            ("▶ Run FBM", isostasy_cache.valid && !any_running, fbm_running)
        }
        PipelinePhase::Erosion => {
            ("▶ Run Erosion", upscale_cache.heightmap.is_some() && !any_running, erosion_running)
        }
        PipelinePhase::Hydrology => (
            "▶ Compute Hydrology",
            (erosion_cache.heightmap.is_some() || upscale_cache.heightmap.is_some())
                && !any_running,
            false,
        ),
        PipelinePhase::Climate => ("▶ Run Climate", false, false),
    }
}

fn dispatch_run(
    phase: PipelinePhase,
    ui_actions: &mut crate::state::UiActions,
    run_timer: &mut crate::state::RunTimer,
) {
    match phase {
        PipelinePhase::Tectonics => ui_actions.run_solver_requested = true,
        PipelinePhase::Isostasy => {
            ui_actions.center_requested = true;
            return; // Centering is instant — no timer
        }
        PipelinePhase::UpscaleFbm => ui_actions.run_fbm_requested = true,
        PipelinePhase::Erosion => ui_actions.run_erosion_requested = true,
        PipelinePhase::Hydrology => ui_actions.run_hydrology_requested = true,
        PipelinePhase::Climate => return,
    }
    run_timer.started_at = Some(std::time::Instant::now());
    run_timer.last_elapsed = None;
}

fn current_progress(
    phase: PipelinePhase,
    bridge: &crate::bridge::SolverBridge,
    upscale_cache: &crate::state::UpscaleCache,
    erosion_cache: &crate::state::ErosionCache,
    flow_cache: &crate::state::FlowCache,
) -> Option<(f32, String)> {
    use crate::bridge::plugin::SolverState;
    use crate::state::{ErosionState, FbmState, FlowState};

    match phase {
        PipelinePhase::Tectonics => {
            if let SolverState::Running { step, total_steps, .. } = &bridge.state {
                let frac = *step as f32 / (*total_steps).max(1) as f32;
                Some((frac, format!("{}/{}", step, total_steps)))
            } else {
                None
            }
        }
        PipelinePhase::UpscaleFbm => {
            if matches!(upscale_cache.state, FbmState::Running) {
                Some((0.0, "Running...".into()))
            } else {
                None
            }
        }
        PipelinePhase::Erosion => {
            if let ErosionState::Running { completed, total } = &erosion_cache.state {
                let frac = *completed as f32 / (*total).max(1) as f32;
                Some((frac, format!("{:.1}M/{:.1}M", *completed as f64 / 1e6, *total as f64 / 1e6)))
            } else {
                None
            }
        }
        PipelinePhase::Hydrology => {
            if matches!(flow_cache.state, FlowState::Running) {
                Some((0.0, "Computing flow...".into()))
            } else {
                None
            }
        }
        _ => None,
    }
}

fn completed_elapsed(
    phase: PipelinePhase,
    bridge: &crate::bridge::SolverBridge,
    upscale_cache: &crate::state::UpscaleCache,
    erosion_cache: &crate::state::ErosionCache,
    flow_cache: &crate::state::FlowCache,
) -> Option<std::time::Duration> {
    use crate::bridge::plugin::SolverState;
    use crate::state::{ErosionState, FbmState, FlowState};

    match phase {
        PipelinePhase::Tectonics => {
            if let SolverState::Completed { elapsed } = &bridge.state {
                Some(*elapsed)
            } else {
                None
            }
        }
        PipelinePhase::UpscaleFbm => {
            if let FbmState::Completed { elapsed } = &upscale_cache.state {
                Some(*elapsed)
            } else {
                None
            }
        }
        PipelinePhase::Erosion => {
            if let ErosionState::Completed { elapsed } = &erosion_cache.state {
                Some(*elapsed)
            } else {
                None
            }
        }
        PipelinePhase::Hydrology => {
            if let FlowState::Completed { elapsed } = &flow_cache.state {
                Some(*elapsed)
            } else {
                None
            }
        }
        _ => None,
    }
}

#[derive(bevy::ecs::system::SystemParam)]
struct UiRightPanelParams<'w> {
    erosion: ResMut<'w, crate::state::ErosionParams>,
    climate: ResMut<'w, crate::state::ClimateParams>,
    gen_params: ResMut<'w, GenerationParamsUi>,
    tectonic_state: Option<ResMut<'w, crate::state::TectonicState>>,
    solver_config: ResMut<'w, crate::state::SolverConfig>,
    bridge: ResMut<'w, crate::bridge::SolverBridge>,
    isostasy_params: ResMut<'w, crate::state::IsostasyParams>,
    isostasy_cache: Res<'w, crate::state::IsostasyCache>,
    fbm_params: ResMut<'w, crate::state::FbmParams>,
    upscale_cache: ResMut<'w, crate::state::UpscaleCache>,
    erosion_cache: ResMut<'w, crate::state::ErosionCache>,
    flow_cache: ResMut<'w, crate::state::FlowCache>,
    lake_cache: ResMut<'w, crate::state::LakeCache>,
    centering: ResMut<'w, crate::state::CenteringState>,
    ui_actions: ResMut<'w, crate::state::UiActions>,
    terrain_display: Res<'w, crate::visualization::render::TerrainDisplay>,
    run_timer: ResMut<'w, crate::state::RunTimer>,
}

fn ui_right_panel(
    mut contexts: EguiContexts,
    current_phase: Res<State<PipelinePhase>>,
    mut next_phase: ResMut<NextState<PipelinePhase>>,
    pipeline_state: Res<crate::state::PipelineState>,
    mut params: UiRightPanelParams,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    egui::SidePanel::left("pipeline_toolbar").exact_width(55.0).resizable(false).show(ctx, |ui| {
        left_toolbar::draw(ui, current_phase.get(), &mut next_phase);
    });

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
                &mut params.erosion_cache,
                &mut params.flow_cache,
                &mut params.lake_cache,
                &mut params.ui_actions,
                &mut params.centering,
                &mut params.run_timer,
            );
        });
    });
}

#[allow(clippy::too_many_arguments)]
fn ui_bottom_bar(
    mut contexts: EguiContexts,
    time: Res<Time>,
    isostasy_cache: Res<crate::state::IsostasyCache>,
    upscale_cache: Res<crate::state::UpscaleCache>,
    erosion_cache: Res<crate::state::ErosionCache>,
    flow_cache: Res<crate::state::FlowCache>,
    bridge: Res<crate::bridge::SolverBridge>,
    current_phase: Res<State<PipelinePhase>>,
    mut run_timer: ResMut<crate::state::RunTimer>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    let phase = *current_phase.get();

    let progress = current_progress(phase, &bridge, &upscale_cache, &erosion_cache, &flow_cache);
    let completed_elapsed =
        completed_elapsed(phase, &bridge, &upscale_cache, &erosion_cache, &flow_cache);

    // Detect completion: if timer is running but nothing is running anymore, stop it
    if run_timer.started_at.is_some() && progress.is_none() {
        run_timer.started_at = None;
    }

    egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            // Left: grid info
            let grid_info = match phase {
                PipelinePhase::Tectonics | PipelinePhase::Isostasy => isostasy_cache
                    .heightmap
                    .as_ref()
                    .map(|hm| (hm.width, hm.height, isostasy_cache.land_ratio)),
                PipelinePhase::UpscaleFbm => upscale_cache
                    .heightmap
                    .as_ref()
                    .map(|hm| (hm.width, hm.height, isostasy_cache.land_ratio)),
                _ => {
                    let hm = erosion_cache.heightmap.as_ref().or(upscale_cache.heightmap.as_ref());
                    hm.map(|hm| (hm.width, hm.height, isostasy_cache.land_ratio))
                }
            };

            if let Some((w, h, land)) = grid_info {
                let src = isostasy_cache.heightmap.as_ref().map(|hm| hm.width).unwrap_or(0);
                if src > 0 && w != src {
                    ui.monospace(format!("{src}\u{b2}\u{2192}{w}\u{b2}"));
                } else {
                    ui.monospace(format!("{w}\u{00d7}{h}"));
                }
                ui.separator();
                ui.monospace(format!("{:.0}% land", land * 100.0));
            }

            // Center: progress bar + timer
            if let Some((frac, label)) = progress {
                ui.separator();
                ui.add(egui::ProgressBar::new(frac).text(label).desired_width(180.0));
                if let Some(started) = run_timer.started_at {
                    let secs = started.elapsed().as_secs_f64();
                    ui.monospace(format!("{secs:.1}s"));
                    ctx.request_repaint();
                }
            } else if let Some(elapsed) = completed_elapsed {
                ui.separator();
                ui.monospace(format!("\u{2713} {:.1}s", elapsed.as_secs_f64()));
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let fps = 1.0 / time.delta_secs().max(0.001);
                ui.monospace(format!("{fps:.0} FPS"));
            });
        });
    });
}

fn ui_toasts(
    mut contexts: EguiContexts,
    mut toasts: ResMut<crate::state::Toasts>,
    mut ui_actions: ResMut<crate::state::UiActions>,
) {
    // Bridge: move last_message into toasts
    if let Some((msg, _when, success)) = ui_actions.last_message.take() {
        toasts.add(msg, success);
    }

    let Ok(ctx) = contexts.ctx_mut() else { return };
    toasts.prune();

    if toasts.messages.is_empty() {
        return;
    }

    egui::Area::new(egui::Id::new("toasts"))
        .anchor(egui::Align2::RIGHT_TOP, [-(260.0 + 10.0), 40.0])
        .show(ctx, |ui| {
            for (msg, created, success) in &toasts.messages {
                let age = created.elapsed().as_secs_f32();
                let alpha = if age > 4.0 { ((5.0 - age) * 255.0) as u8 } else { 255 };
                let bg = if *success {
                    egui::Color32::from_rgba_unmultiplied(20, 80, 20, alpha)
                } else {
                    egui::Color32::from_rgba_unmultiplied(120, 30, 30, alpha)
                };
                egui::Frame::new().fill(bg).corner_radius(4.0).inner_margin(8.0).show(ui, |ui| {
                    ui.label(
                        egui::RichText::new(msg)
                            .color(egui::Color32::from_white_alpha(alpha))
                            .small(),
                    );
                });
                ui.add_space(2.0);
            }
        });
}
