use bevy::prelude::ResMut;
use bevy_egui::egui;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::plates::{PlateConfig, generate_plates};

use crate::state::{
    ClimateParams, ErosionParams, GenerationParamsUi, PipelinePhase, TectonicState,
    TectonicsParams, ViewState,
};

pub fn draw(
    ui: &mut egui::Ui,
    view_state: &ViewState,
    erosion: &mut ErosionParams,
    tectonics: &mut TectonicsParams,
    climate: &mut ClimateParams,
    generation: &mut GenerationParamsUi,
    tectonic_state: Option<ResMut<TectonicState>>,
) {
    egui::CollapsingHeader::new("Parameters").default_open(true).show(ui, |ui| {
        match view_state.selected_phase {
            PipelinePhase::Erosion => draw_erosion(ui, erosion),
            PipelinePhase::Tectonics => draw_tectonics(ui, tectonic_state),
            PipelinePhase::Climate => draw_climate(ui, climate),
            _ => {
                ui.label("No parameters for this phase");
            }
        }

        // Suppress unused-variable warning — TectonicsParams is kept for the
        // solver phase (future M1 work) and is intentionally not shown here.
        let _ = tectonics;

        ui.add_space(8.0);
        ui.separator();
        ui.strong("Generation");
        slider_row_u64(ui, "Seed", &mut generation.seed);
        slider_row(ui, "m/pixel", &mut generation.meters_per_pixel, 20.0..=100.0);
        slider_row(ui, "Size (km)", &mut generation.continent_size_km, 50.0..=500.0);
        slider_row(ui, "Max elev (m)", &mut generation.max_elevation_m, 1000.0..=6000.0);
    });
}

// ── Tectonics panel ───────────────────────────────────────────────────────

fn draw_tectonics(ui: &mut egui::Ui, tectonic_state: Option<ResMut<TectonicState>>) {
    let Some(mut state) = tectonic_state else {
        ui.label("Tectonic state not ready");
        return;
    };

    // Presets
    ui.strong("Presets");
    ui.horizontal_wrapped(|ui| {
        if ui.button("Default").clicked() {
            apply_preset(&mut state, PlateConfig::default());
        }
        if ui.button("Continent").clicked() {
            apply_preset(&mut state, PlateConfig::preset_single_continent());
        }
        if ui.button("Collision").clicked() {
            apply_preset(&mut state, PlateConfig::preset_collision());
        }
        if ui.button("Archipelago").clicked() {
            apply_preset(&mut state, PlateConfig::preset_archipelago());
        }
        if ui.button("Rift").clicked() {
            apply_preset(&mut state, PlateConfig::preset_rift());
        }
    });

    ui.add_space(4.0);
    ui.separator();

    // Sliders — modify config in place, no immediate regen
    slider_row_usize(ui, "Plates", &mut state.config.num_plates, 3..=15);
    slider_row(ui, "Cont. ratio", &mut state.config.continental_ratio, 0.1..=0.6);
    slider_row(ui, "Vel. min", &mut state.config.velocity_min, 0.1..=3.0);
    slider_row(ui, "Vel. max", &mut state.config.velocity_max, 0.5..=5.0);
    slider_row(ui, "Smoothing σ", &mut state.config.boundary_smoothing_sigma, 0.0..=5.0);

    // Grid size as discrete steps
    let grid_sizes = [64usize, 128, 256, 512];
    ui.horizontal(|ui| {
        ui.label("Grid size");
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            for &sz in &grid_sizes {
                if ui
                    .add(egui::Button::new(format!("{sz}")).selected(state.config.grid_size == sz))
                    .clicked()
                {
                    state.config.grid_size = sz;
                }
            }
        });
    });

    ui.add_space(4.0);

    // Seed + action buttons
    ui.horizontal(|ui| {
        ui.label("Seed");
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add(egui::DragValue::new(&mut state.seed).speed(1.0));
        });
    });

    ui.add_space(4.0);

    ui.horizontal(|ui| {
        if ui.button("⟳  Generate").clicked() {
            regenerate(&mut state);
        }
        if ui.button("🎲  Randomize").clicked() {
            state.seed = state.seed.wrapping_add(1);
            regenerate(&mut state);
        }
    });

    ui.add_space(6.0);
    ui.separator();

    // Info summary
    let init = &state.init;
    let cont = init.plates.iter().filter(|p| {
        p.plate_type == ymir_core::tectonics::plates::PlateType::Continental
    }).count();
    let oce = init.plates.len() - cont;
    let avg_vel: f32 = init
        .plates
        .iter()
        .map(|p| (p.velocity.0 * p.velocity.0 + p.velocity.1 * p.velocity.1).sqrt())
        .sum::<f32>()
        / init.plates.len().max(1) as f32;
    let t_min = init.thickness.min();
    let t_max = init.thickness.max();
    let t_mean = init.thickness.mean();

    egui::Grid::new("tectonic_info").num_columns(2).show(ui, |ui| {
        ui.label("Continental");
        ui.monospace(format!("{cont}"));
        ui.end_row();
        ui.label("Oceanic");
        ui.monospace(format!("{oce}"));
        ui.end_row();
        ui.label("Avg velocity");
        ui.monospace(format!("{avg_vel:.2}"));
        ui.end_row();
        ui.label("Thickness");
        ui.monospace(format!("{t_min:.2} / {t_mean:.2} / {t_max:.2}"));
        ui.end_row();
    });
}

fn apply_preset(state: &mut TectonicState, config: PlateConfig) {
    state.config = config;
    regenerate(state);
}

fn regenerate(state: &mut TectonicState) {
    let seed = WorldSeed::new(state.seed);
    state.init = generate_plates(&state.config, &seed);
    state.dirty = true;
    state.generation = state.generation.wrapping_add(1);
}

// ── Erosion & climate panels ──────────────────────────────────────────────

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

fn draw_climate(ui: &mut egui::Ui, p: &mut ClimateParams) {
    slider_row(ui, "Base temp (C)", &mut p.base_temperature, -10.0..=30.0);
    slider_row(ui, "Wind dir (deg)", &mut p.wind_direction_deg, 0.0..=360.0);
    slider_row(ui, "Orographic factor", &mut p.orographic_factor, 1.0..=5.0);
    slider_row(ui, "Moisture decay (km)", &mut p.moisture_decay_km, 100.0..=800.0);
}

// ── Slider helpers ────────────────────────────────────────────────────────

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

fn slider_row_usize(
    ui: &mut egui::Ui,
    label: &str,
    val: &mut usize,
    range: std::ops::RangeInclusive<usize>,
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
