//! UI module — v2-only after Step 8.6 Phase 8h sunset.
//!
//! Pre-sunset legacy panels (`parameter_panel`, `pipeline_panel`,
//! `statistics_panel`, `left_toolbar`) drove the legacy pipeline
//! phases (TEC / ISO / FBM / ERO / HYD / CLI). Post-sunset the
//! binary surfaces only the v2 panels: a top bar (status + step
//! counter + wallclock), a left panel (real-time metrics dashboard),
//! and a right panel (parameter editor including init mode + run /
//! cancel / capture / export / import).

pub mod metrics_dashboard;
pub mod parameter_panel_v2;
pub mod pipeline_toolbar;
pub mod workflow_panel;

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPlugin, EguiPrimaryContextPass, egui};

use parameter_panel_v2::V2EditableSpec;
use workflow_panel::{WorkflowCycleHistory, WorkflowExportState};

pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin::default());
        app.init_resource::<V2EditableSpec>();
        app.init_resource::<WorkflowCycleHistory>();
        app.init_resource::<WorkflowExportState>();
        app.add_systems(
            EguiPrimaryContextPass,
            (
                configure_egui_style,
                pipeline_toolbar::ui_phase_toolbar,
                ui_v2_top_bar
                    .run_if(resource_exists::<crate::bridge::v2::V2SolverBridge>),
                ui_v2_right_panel
                    .run_if(resource_exists::<crate::bridge::v2::V2SolverBridge>),
                ui_v2_left_panel
                    .run_if(resource_exists::<crate::bridge::v2::V2SolverBridge>),
            ),
        );
    }
}

/// Step 8.6 v2 top bar — status badge + (during a run) step counter,
/// elapsed wall time, progress bar.
fn ui_v2_top_bar(mut contexts: EguiContexts, bridge: Res<crate::bridge::v2::V2SolverBridge>) {
    use crate::bridge::v2::V2RunState;
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::TopBottomPanel::top("v2_top_bar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new("Ymir v2").strong());
            ui.separator();
            let badge = match &bridge.state {
                V2RunState::Idle => ("Idle", egui::Color32::GRAY),
                V2RunState::Running { .. } => ("Running", egui::Color32::YELLOW),
                V2RunState::Completed { .. } => ("Completed", egui::Color32::GREEN),
                V2RunState::Imported { .. } => ("Imported", egui::Color32::LIGHT_BLUE),
                V2RunState::Failed { .. } => ("Failed", egui::Color32::RED),
                V2RunState::WorkflowPhaseACompleted { .. } => {
                    ("Phase A done", egui::Color32::LIGHT_GREEN)
                }
                V2RunState::WorkflowPhaseBCompleted { .. } => {
                    ("Phase B done", egui::Color32::from_rgb(0x80, 0xD0, 0xFF))
                }
            };
            ui.colored_label(badge.1, badge.0);

            if let V2RunState::Running { step, total, started_at, .. } = &bridge.state {
                ui.separator();
                let frac = if *total > 0 {
                    (*step as f32 / *total as f32).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let label = format!("{}/{}", step, total);
                ui.add(egui::ProgressBar::new(frac).text(label).desired_width(220.0));
                if let Some(start) = started_at {
                    let secs = start.elapsed().as_secs_f64();
                    ui.monospace(format!("{secs:.1}s"));
                    ctx.request_repaint();
                }
            } else if let V2RunState::Completed { elapsed, metrics, .. } = &bridge.state {
                ui.separator();
                ui.monospace(format!(
                    "\u{2713} {:.1}s · CG mean {:.0} · peak|v| {:.2e}",
                    elapsed.as_secs_f64(),
                    metrics.cg_iter_mean,
                    metrics.vmax_peak
                ));
            } else if let V2RunState::Imported {
                elapsed,
                scalar_metrics,
                exported_at,
                ..
            } = &bridge.state
            {
                ui.separator();
                ui.monospace(format!(
                    "\u{1f4c2} {} · {:.1}s · CG mean {:.0} · peak|v| {:.2e}",
                    exported_at,
                    elapsed.as_secs_f64(),
                    scalar_metrics.cg_iter_mean,
                    scalar_metrics.vmax_peak
                ));
            }
        });
    });
}

/// Step 8.6 v2 right panel — wraps the v2 parameter editor + the
/// active phase's collapsible config section.
#[allow(clippy::too_many_arguments)]
fn ui_v2_right_panel(
    mut contexts: EguiContexts,
    mut spec_state: ResMut<V2EditableSpec>,
    mut bridge: ResMut<crate::bridge::v2::V2SolverBridge>,
    mut viz: ResMut<crate::visualization::v2_viz::V2VizState>,
    active: Res<crate::pipeline::ActivePhase>,
    mut isostasy_params: ResMut<crate::phases::isostasy::IsostasyParams>,
    isostasy_cache: Res<crate::phases::isostasy::IsostasyCache>,
    mut fbm_params: ResMut<crate::phases::upscale_fbm::FbmParams>,
    fbm_cache: Res<crate::phases::upscale_fbm::FbmCache>,
    mut erosion_params: ResMut<crate::phases::erosion::ErosionParams>,
    erosion_cache: Res<crate::phases::erosion::ErosionCache>,
    mut hydrology_params: ResMut<crate::phases::hydrology::HydrologyParams>,
    hydrology_cache: Res<crate::phases::hydrology::HydrologyCache>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::SidePanel::right("v2_right_panel").exact_width(300.0).show(ctx, |ui| {
        egui::ScrollArea::vertical().show(ui, |ui| {
            parameter_panel_v2::draw(
                ui,
                &mut spec_state,
                &mut bridge,
                &mut viz,
                *active,
                &mut isostasy_params,
                &isostasy_cache,
                &mut fbm_params,
                &fbm_cache,
                &mut erosion_params,
                &erosion_cache,
                &mut hydrology_params,
                &hydrology_cache,
            );
        });
        // While any threaded phase worker is running, force egui /
        // Bevy to redraw every frame so the streamed in-progress
        // snapshots (erosion preview heightmaps, etc.) actually
        // hit the screen instead of waiting on the next user
        // input event. Mirrors the v2-Running repaint request in
        // `ui_v2_left_panel`.
        if matches!(
            erosion_cache.state,
            crate::phases::erosion::ErosionState::Running { .. }
        ) {
            ctx.request_repaint();
        }
    });
}

/// Step 8.6 Phase 8c — v2 left panel: real-time nondimensional metrics
/// dashboard. Live during a run (peek-state derived metrics + progress
/// + ETA), final summary post-run (Metrics struct), preview summary
/// (config recap) when the bridge is Idle and a `V2Preview` is in
/// flight.
fn ui_v2_left_panel(
    mut contexts: EguiContexts,
    bridge: Res<crate::bridge::v2::V2SolverBridge>,
    viz: Res<crate::visualization::v2_viz::V2VizState>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::SidePanel::left("v2_left_panel").default_width(280.0).show(ctx, |ui| {
        egui::ScrollArea::vertical().show(ui, |ui| {
            metrics_dashboard::draw(ui, &bridge, &viz);
        });
        if matches!(bridge.state, crate::bridge::v2::V2RunState::Running { .. }) {
            ctx.request_repaint();
        }
    });
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
