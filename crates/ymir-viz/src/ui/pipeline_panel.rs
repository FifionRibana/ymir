use bevy_egui::egui;

use crate::state::{PhaseStatus, PipelineState, ViewState};

pub fn draw(ui: &mut egui::Ui, view_state: &mut ViewState, pipeline: &PipelineState) {
    egui::CollapsingHeader::new("Pipeline").default_open(true).show(ui, |ui| {
        for (phase, status, timing) in &pipeline.phases {
            let selected = view_state.selected_phase == *phase;
            let text = format!("{} {}", status.icon(), phase.label());
            let response = ui.selectable_label(selected, &text);
            if response.clicked() {
                view_state.selected_phase = *phase;
            }
            if !timing.is_empty() {
                response.on_hover_text(timing);
            }
        }

        // Progress bar for running phase
        for (_, status, _) in &pipeline.phases {
            if let PhaseStatus::Running(pct) = status {
                ui.add(egui::ProgressBar::new(*pct as f32 / 100.0).show_percentage().animate(true));
                break;
            }
        }

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            if ui.button("Run all").clicked() {
                bevy::log::info!("Run all clicked");
            }
            if ui.button("Stop").clicked() {
                bevy::log::info!("Stop clicked");
            }
            if ui.button("Export").clicked() {
                bevy::log::info!("Export clicked");
            }
        });
    });
}
