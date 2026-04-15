use bevy::prelude::*;
use bevy_egui::egui;

use crate::state::PipelinePhase;

pub fn draw(
    ui: &mut egui::Ui,
    current_phase: &PipelinePhase,
    next_phase: &mut NextState<PipelinePhase>,
) {
    ui.vertical_centered(|ui| {
        ui.add_space(4.0);
        for &phase in PipelinePhase::ALL {
            let selected = *current_phase == phase;
            let btn = ui.add(
                egui::Button::new(egui::RichText::new(phase.short_label()).size(11.0).strong())
                    .selected(selected)
                    .min_size(egui::vec2(45.0, 32.0)),
            );
            if btn.clicked() {
                next_phase.set(phase);
            }
            btn.on_hover_text(phase.label());
        }
    });
}
