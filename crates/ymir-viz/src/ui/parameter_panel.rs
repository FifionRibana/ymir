use std::sync::atomic::Ordering;

use bevy::prelude::ResMut;
use bevy_egui::egui;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::plates::{PlateConfig, generate_plates};
use ymir_core::tectonics::solver::config::{
    ContinuationConfig, NewtonConfig, NonlinearSolver, PicardConfig, Preconditioner,
    TectonicsConfig,
};
use ymir_core::tectonics::solver::tectonics::DynamicPlateContext;

use crate::bridge::commands::SolverCommand;
use crate::bridge::plugin::{SolverBridge, SolverState};
use crate::state::{
    ClimateParams, ErosionParams, GenerationParamsUi, IsostasyCache, IsostasyParams, PipelinePhase,
    SolverConfig, TectonicState, ViewState,
};

#[allow(clippy::too_many_arguments)]
pub fn draw(
    ui: &mut egui::Ui,
    view_state: &ViewState,
    erosion: &mut ErosionParams,
    climate: &mut ClimateParams,
    generation: &mut GenerationParamsUi,
    tectonic_state: Option<ResMut<TectonicState>>,
    solver_config: &mut SolverConfig,
    bridge: &mut SolverBridge,
    isostasy_params: &mut IsostasyParams,
    isostasy_cache: &IsostasyCache,
) {
    egui::CollapsingHeader::new("Parameters").default_open(true).show(ui, |ui| {
        match view_state.selected_phase {
            PipelinePhase::Erosion => draw_erosion(ui, erosion),
            PipelinePhase::Tectonics => {
                draw_tectonics(
                    ui,
                    tectonic_state,
                    solver_config,
                    bridge,
                    isostasy_params,
                    isostasy_cache,
                );
            }
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

// ── Tectonics panel ───────────────────────────────────────────────────────

fn draw_tectonics(
    ui: &mut egui::Ui,
    tectonic_state: Option<ResMut<TectonicState>>,
    solver_config: &mut SolverConfig,
    bridge: &mut SolverBridge,
    isostasy_params: &mut IsostasyParams,
    isostasy_cache: &IsostasyCache,
) {
    let Some(mut state) = tectonic_state else {
        ui.label("Tectonic state not ready");
        return;
    };

    // ── Plate generation ──
    ui.strong("Plate generation");
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

    slider_row_usize(ui, "Plates", &mut state.config.num_plates, 3..=15);
    slider_row(ui, "Cont. ratio", &mut state.config.continental_ratio, 0.1..=0.6);
    slider_row(ui, "Vel. min", &mut state.config.velocity_min, 0.1..=3.0);
    slider_row(ui, "Vel. max", &mut state.config.velocity_max, 0.5..=5.0);
    slider_row(ui, "Smoothing σ", &mut state.config.boundary_smoothing_sigma, 0.0..=5.0);

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

    // ── Info ──
    let init = &state.init;
    let cont = init
        .plates
        .iter()
        .filter(|p| p.plate_type == ymir_core::tectonics::plates::PlateType::Continental)
        .count();
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

    ui.add_space(6.0);
    ui.separator();

    // ── Solver ──
    ui.strong("Solver");

    slider_row_usize(ui, "Timesteps", &mut solver_config.num_timesteps, 50..=1000);
    slider_row_f64(ui, "Gravity", &mut solver_config.gravity_factor, 0.1..=5.0);
    slider_row_f64(ui, "CFL", &mut solver_config.cfl_factor, 0.05..=0.9);
    slider_row_f64(ui, "Relaxation", &mut solver_config.picard_relaxation, 0.3..=1.0);

    ui.horizontal(|ui| {
        ui.label("Viscosity");
        let mut idx = if solver_config.power_law_n < 2.0 { 0 } else { 1 };
        egui::ComboBox::from_id_salt("power_law")
            .selected_text(if idx == 0 { "Linear (n=1)" } else { "Power-law (n=3)" })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut idx, 0, "Linear (n=1)");
                ui.selectable_value(&mut idx, 1, "Power-law (n=3)");
            });
        solver_config.power_law_n = if idx == 0 { 1.0 } else { 3.0 };
    });

    ui.horizontal(|ui| {
        ui.label("Solver");
        let mut idx = match solver_config.nonlinear_solver {
            NonlinearSolver::Picard => 0,
            NonlinearSolver::Newton => 1,
        };
        egui::ComboBox::from_id_salt("nl_solver")
            .selected_text(if idx == 0 { "Picard" } else { "Newton (JFNK)" })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut idx, 0, "Picard");
                ui.selectable_value(&mut idx, 1, "Newton (JFNK)");
            });
        solver_config.nonlinear_solver =
            if idx == 0 { NonlinearSolver::Picard } else { NonlinearSolver::Newton };
    });

    ui.separator();
    ui.strong("Convergence");

    ui.checkbox(&mut solver_config.continuation_enabled, "Viscosity continuation");
    ui.add(
        egui::Slider::new(&mut solver_config.strain_rate_min, 1e-6..=1e-1)
            .text("ε_min")
            .logarithmic(true),
    );
    ui.add(
        egui::Slider::new(&mut solver_config.eta_max, 1e2..=1e6).text("η_max").logarithmic(true),
    );

    ui.horizontal(|ui| {
        ui.label("Precond.");
        let mut is_ssor = matches!(solver_config.preconditioner, Preconditioner::Ssor { .. });
        egui::ComboBox::from_id_salt("precond")
            .selected_text(if is_ssor { "SSOR" } else { "Jacobi" })
            .show_ui(ui, |ui| {
                if ui.selectable_value(&mut is_ssor, false, "Jacobi").changed() && !is_ssor {
                    solver_config.preconditioner = Preconditioner::Jacobi;
                }
                if ui.selectable_value(&mut is_ssor, true, "SSOR").changed() && is_ssor {
                    solver_config.preconditioner = Preconditioner::Ssor { omega: 1.2 };
                }
            });
    });
    if let Preconditioner::Ssor { ref mut omega } = solver_config.preconditioner {
        ui.add(egui::Slider::new(omega, 0.5..=1.9).text("omega"));
    }
    ui.checkbox(&mut solver_config.inexact_newton, "Inexact Newton");

    ui.add_space(4.0);
    ui.separator();
    ui.strong("Boundary processes");

    ui.checkbox(&mut solver_config.boundaries.enabled, "Enabled");
    ui.checkbox(&mut solver_config.dynamic_boundaries, "Dynamic boundaries");
    if solver_config.boundaries.enabled {
        slider_row_f64(
            ui,
            "Subduction rate",
            &mut solver_config.boundaries.subduction_rate,
            0.0..=2.0,
        );
        slider_row_f64(
            ui,
            "Volcanic arc rate",
            &mut solver_config.boundaries.volcanic_arc_rate,
            0.0..=1.0,
        );
        slider_row_f64(
            ui,
            "Spreading rate",
            &mut solver_config.boundaries.spreading_rate,
            0.0..=2.0,
        );
        slider_row_f64(
            ui,
            "Rift threshold",
            &mut solver_config.boundaries.rift_thickness_threshold,
            0.1..=0.8,
        );
        slider_row_f64(
            ui,
            "Source smooth. σ",
            &mut solver_config.boundaries.source_smoothing_sigma,
            0.0..=5.0,
        );
        slider_row_f64(
            ui,
            "Ocean ref. S",
            &mut solver_config.boundaries.oceanic_reference_thickness,
            0.1..=0.5,
        );
        slider_row_f64(
            ui,
            "Ocean restore",
            &mut solver_config.boundaries.oceanic_restore_rate,
            0.0..=1.0,
        );

        ui.separator();
        ui.strong("Density & slab pull");

        slider_row_f64(
            ui,
            "ρ continental",
            &mut solver_config.boundaries.rho_continental,
            2500.0..=3000.0,
        );
        slider_row_f64(ui, "ρ oceanic", &mut solver_config.boundaries.rho_oceanic, 2800.0..=3200.0);
        slider_row_f64(ui, "ρ mantle", &mut solver_config.boundaries.rho_mantle, 3100.0..=3500.0);

        ui.checkbox(&mut solver_config.boundaries.slab_pull_enabled, "Slab pull");
        if solver_config.boundaries.slab_pull_enabled {
            slider_row_f64(
                ui,
                "Pull factor",
                &mut solver_config.boundaries.slab_pull_factor,
                0.001..=0.5,
            );
        }
    }

    ui.add_space(4.0);

    let is_running = matches!(bridge.state, SolverState::Running { .. });

    ui.horizontal(|ui| {
        let run_btn = ui.add_enabled(!is_running, egui::Button::new("▶ Run solver"));
        if run_btn.clicked() {
            launch_solver(&state, solver_config, bridge);
        }
        if ui.add_enabled(is_running, egui::Button::new("⏹ Cancel")).clicked() {
            bridge.cancel_flag.store(true, Ordering::Relaxed);
        }
    });

    ui.add_space(4.0);

    // ── Solver status ──
    draw_solver_status(ui, &bridge.state);

    ui.add_space(6.0);
    ui.separator();

    // ── Isostasy ──
    ui.strong("Isostasy");
    slider_row(ui, "Sea level", &mut isostasy_params.sea_level_fraction, 0.0..=1.0);
    slider_row(ui, "Max elev (m)", &mut isostasy_params.max_elevation_m, 1000.0..=6000.0);
    slider_row(ui, "Max depth (m)", &mut isostasy_params.max_depth_m, 100.0..=1000.0);
    slider_row(ui, "Smoothing σ", &mut isostasy_params.altitude_smoothing_sigma, 0.0..=5.0);

    if isostasy_cache.valid {
        ui.add_space(4.0);
        egui::Grid::new("isostasy_info").num_columns(2).show(ui, |ui| {
            ui.label("Land ratio");
            ui.monospace(format!("{:.1}%", isostasy_cache.land_ratio * 100.0));
            ui.end_row();
            ui.label("Peak");
            ui.monospace(format!("{:.0} m", isostasy_cache.peak_altitude_m));
            ui.end_row();
            ui.label("Deepest");
            ui.monospace(format!("-{:.0} m", isostasy_cache.max_depth_m));
            ui.end_row();
        });
    }
}

fn draw_solver_status(ui: &mut egui::Ui, state: &SolverState) {
    match state {
        SolverState::Idle => {
            ui.small("⚪ Ready");
        }
        SolverState::Running { step, total_steps, stats } => {
            let frac = *step as f32 / (*total_steps).max(1) as f32;
            ui.small(format!("🟡 Running… step {}/{}", step, total_steps));
            ui.add(egui::ProgressBar::new(frac).show_percentage());
            if let Some(s) = stats {
                ui.small(format!(
                    "max_v={:.4}  S=[{:.3}, {:.3}]  picard={}  dt={:.2e}",
                    s.max_velocity, s.min_thickness, s.max_thickness, s.picard_iterations, s.dt
                ));
            }
        }
        SolverState::Completed { elapsed } => {
            ui.small(format!("🟢 Done in {:.1}s", elapsed.as_secs_f64()));
        }
        SolverState::Failed { error } => {
            ui.colored_label(egui::Color32::RED, format!("❌ {error}"));
        }
    }
}

fn launch_solver(
    tectonic_state: &TectonicState,
    solver_config: &SolverConfig,
    bridge: &mut SolverBridge,
) {
    let init = &tectonic_state.init;
    let traction = init.to_traction_field();
    let grid_size = init.grid_size;
    let dx = 1.0 / grid_size as f64;

    let initial_s = {
        let mut field = ymir_core::tectonics::solver::field::Field2D::new(grid_size);
        for j in 0..grid_size {
            for i in 0..grid_size {
                field.set(i, j, init.thickness.data[j * grid_size + i] as f64);
            }
        }
        field
    };

    let config = TectonicsConfig {
        num_timesteps: solver_config.num_timesteps,
        gravity_factor: solver_config.gravity_factor,
        cfl_factor: solver_config.cfl_factor,
        s_min: 0.1,
        s_max: 2.5,
        nonlinear_solver: solver_config.nonlinear_solver,
        picard: PicardConfig {
            power_law_n: solver_config.power_law_n,
            relaxation: solver_config.picard_relaxation,
            strain_rate_min: solver_config.strain_rate_min,
            eta_max: solver_config.eta_max,
            ..PicardConfig::default()
        },
        newton: NewtonConfig {
            preconditioner: solver_config.preconditioner,
            inexact: solver_config.inexact_newton,
            ..NewtonConfig::default()
        },
        continuation: ContinuationConfig {
            enabled: solver_config.continuation_enabled,
            ..ContinuationConfig::default()
        },
        boundaries: solver_config.boundaries.clone(),
        dynamic_boundaries: solver_config.dynamic_boundaries,
    };

    let plate_ctx =
        DynamicPlateContext { ids: init.plate_ids.clone(), plates: init.plates.clone(), traction };

    let _ = bridge.commands_tx.send(SolverCommand::RunTectonics {
        config,
        plate_ctx,
        initial_s,
        grid_size,
        dx,
    });

    bridge.state =
        SolverState::Running { step: 0, total_steps: solver_config.num_timesteps, stats: None };
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

fn slider_row_f64(
    ui: &mut egui::Ui,
    label: &str,
    val: &mut f64,
    range: std::ops::RangeInclusive<f64>,
) {
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
