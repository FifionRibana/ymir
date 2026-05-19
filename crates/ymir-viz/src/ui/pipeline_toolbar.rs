//! Thin left-side phase navigation toolbar.
//!
//! Pre-Phase-8h sunset this lived in `ui/left_toolbar.rs`. The
//! restored version is intentionally narrow (55 px) so it sits
//! flush against the metrics dashboard left panel without crowding
//! the centre sprite. Buttons are 3-letter tags from
//! `PipelinePhase::short_label`; hovering surfaces the long label.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::pipeline::{ActivePhase, PipelinePhase};

pub fn ui_phase_toolbar(
    mut contexts: EguiContexts,
    mut active: ResMut<ActivePhase>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::SidePanel::left("phase_toolbar")
        .exact_width(55.0)
        .resizable(false)
        .show(ctx, |ui| {
            ui.add_space(6.0);
            for &phase in PipelinePhase::ALL {
                let is_active = active.0 == phase;
                let mut button = egui::Button::new(
                    egui::RichText::new(phase.short_label())
                        .monospace()
                        .strong(),
                )
                .min_size(egui::vec2(40.0, 28.0));
                if !phase.is_implemented() {
                    button = button.fill(egui::Color32::from_rgb(0x40, 0x40, 0x44));
                } else if is_active {
                    button = button.fill(egui::Color32::from_rgb(0xB8, 0x73, 0x33));
                }
                let response = ui.add(button).on_hover_text(phase.label());
                if response.clicked() {
                    active.0 = phase;
                }
                ui.add_space(2.0);
            }
        });
}
