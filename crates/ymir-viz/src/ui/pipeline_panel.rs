use std::fs;
use std::path::{Path, PathBuf};

use bevy_egui::egui;

use crate::state::{PhaseStatus, PipelineState, UiActions, ViewState};

pub fn draw(
    ui: &mut egui::Ui,
    view_state: &mut ViewState,
    pipeline: &PipelineState,
    ui_actions: &mut UiActions,
    has_terrain: bool,
) {
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
            if ui.button("Step").clicked() {
                ui_actions.step_requested = true;
            }
            if ui.button("Stop").clicked() {
                bevy::log::info!("Stop clicked");
            }
            if ui.add_enabled(has_terrain, egui::Button::new("Export")).clicked() {
                ui_actions.export_requested = true;
            }
        });

        // Feedback message
        if let Some((ref msg, when, success)) = ui_actions.last_message {
            if when.elapsed().as_secs() < 4 {
                let color = if success { egui::Color32::GREEN } else { egui::Color32::RED };
                ui.colored_label(color, msg);
            } else {
                ui_actions.last_message = None;
            }
        }

        // Load section
        ui.add_space(4.0);
        egui::CollapsingHeader::new("Load saved").show(ui, |ui| {
            let dirs =
                ui_actions.cached_dirs.get_or_insert_with(|| list_export_dirs(Path::new("output")));

            if dirs.is_empty() {
                ui.label("No saved exports in output/");
            } else {
                for dir in dirs.clone() {
                    let name = dir
                        .file_name()
                        .map(|n| n.to_string_lossy().to_string())
                        .unwrap_or_default();
                    ui.horizontal(|ui| {
                        ui.label(&name);
                        if ui.small_button("Load").clicked() {
                            ui_actions.load_requested = Some(dir);
                        }
                    });
                }
            }

            if ui.small_button("Refresh").clicked() {
                ui_actions.cached_dirs = None;
            }
        });
    });
}

fn list_export_dirs(output_root: &Path) -> Vec<PathBuf> {
    if !output_root.exists() {
        return vec![];
    }
    let mut dirs = vec![];
    if let Ok(entries) = fs::read_dir(output_root) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() && path.join("metadata.json").exists() {
                dirs.push(path);
            }
        }
    }
    dirs.sort();
    dirs
}
