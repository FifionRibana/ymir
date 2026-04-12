use bevy_egui::egui;

use crate::state::{
    ClimateParams, ErosionParams, GenerationParamsUi, PipelinePhase, TectonicsParams, ViewState,
};

pub fn draw(
    ui: &mut egui::Ui,
    view_state: &ViewState,
    erosion: &mut ErosionParams,
    tectonics: &mut TectonicsParams,
    climate: &mut ClimateParams,
    generation: &mut GenerationParamsUi,
) {
    egui::CollapsingHeader::new("Parameters").default_open(true).show(ui, |ui| {
        match view_state.selected_phase {
            PipelinePhase::Erosion => draw_erosion(ui, erosion),
            PipelinePhase::Tectonics => draw_tectonics(ui, tectonics),
            PipelinePhase::Climate => draw_climate(ui, climate),
            _ => {
                ui.label("No parameters for this phase");
            }
        }

        ui.add_space(8.0);
        ui.separator();
        ui.strong("Generation");
        slider_row_u64(ui, "Seed", &mut generation.seed);
        slider_row(ui, "m/pixel", &mut generation.meters_per_pixel, 20.0..=100.0);
        slider_row(ui, "Size (km)", &mut generation.continent_size_km, 50.0..=500.0);
        slider_row(ui, "Max elev (m)", &mut generation.max_elevation_m, 1000.0..=6000.0);
    });
}

fn draw_erosion(ui: &mut egui::Ui, p: &mut ErosionParams) {
    slider_row(ui, "Erosion rate", &mut p.erosion_rate, 0.0..=1.0);
    slider_row(ui, "Deposition rate", &mut p.deposition_rate, 0.0..=1.0);
    slider_row(ui, "Inertia", &mut p.inertia, 0.0..=0.5);
    slider_row(ui, "Gravity", &mut p.gravity, 1.0..=15.0);
    slider_row(ui, "Evaporation", &mut p.evaporation_rate, 0.001..=0.05);
    slider_row_u32(ui, "Max lifetime", &mut p.max_lifetime, 50..=300);
    slider_row(ui, "Droplets (M)", &mut p.droplets_millions, 0.5..=20.0);
    slider_row_u32(ui, "Coastal dep.", &mut p.coastal_deposition, 0..=30);
    slider_row(ui, "Min slope", &mut p.min_slope, 0.001..=0.1);
}

fn draw_tectonics(ui: &mut egui::Ui, p: &mut TectonicsParams) {
    slider_row(ui, "Viscosity", &mut p.viscosity, 0.1..=10.0);
    slider_row(ui, "Gravity factor", &mut p.gravity_factor, 0.1..=5.0);
    slider_row_u32(ui, "Timesteps", &mut p.num_timesteps, 50..=1000);
    slider_row_u32(ui, "Power law exp", &mut p.power_law_exponent, 1..=5);
    slider_row_u32(ui, "Plate count", &mut p.plate_count, 3..=15);
    slider_row(ui, "Continental ratio", &mut p.continental_ratio, 0.1..=0.6);
}

fn draw_climate(ui: &mut egui::Ui, p: &mut ClimateParams) {
    slider_row(ui, "Base temp (C)", &mut p.base_temperature, -10.0..=30.0);
    slider_row(ui, "Wind dir (deg)", &mut p.wind_direction_deg, 0.0..=360.0);
    slider_row(ui, "Orographic factor", &mut p.orographic_factor, 1.0..=5.0);
    slider_row(ui, "Moisture decay (km)", &mut p.moisture_decay_km, 100.0..=800.0);
}

fn slider_row(ui: &mut egui::Ui, label: &str, val: &mut f32, range: std::ops::RangeInclusive<f32>) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.monospace(format!("{val:.3}"));
        });
    });
    ui.add(egui::Slider::new(val, range).show_value(false));
}

fn slider_row_u32(
    ui: &mut egui::Ui,
    label: &str,
    val: &mut u32,
    range: std::ops::RangeInclusive<u32>,
) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.monospace(format!("{val}"));
        });
    });
    ui.add(egui::Slider::new(val, range).show_value(false));
}

fn slider_row_u64(ui: &mut egui::Ui, label: &str, val: &mut u64) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add(egui::DragValue::new(val).speed(1.0));
        });
    });
}
