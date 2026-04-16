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
    ClimateParams, ErosionCache, ErosionParams, ErosionState, FbmParams, FbmState, FlowCache,
    FlowState, GenerationParamsUi, IsostasyCache, IsostasyParams, LakeCache, PipelinePhase,
    SolverConfig, TectonicState, UpscaleCache,
};

#[allow(clippy::too_many_arguments)]
pub fn draw(
    ui: &mut egui::Ui,
    current_phase: &PipelinePhase,
    erosion: &mut ErosionParams,
    climate: &mut ClimateParams,
    generation: &mut GenerationParamsUi,
    tectonic_state: Option<ResMut<TectonicState>>,
    solver_config: &mut SolverConfig,
    bridge: &mut SolverBridge,
    isostasy_params: &mut IsostasyParams,
    isostasy_cache: &IsostasyCache,
    fbm_params: &mut FbmParams,
    upscale_cache: &mut UpscaleCache,
    erosion_cache: &mut ErosionCache,
    flow_cache: &mut FlowCache,
    lake_cache: &mut LakeCache,
    ui_actions: &mut crate::state::UiActions,
    centering: &mut crate::state::CenteringState,
    run_timer: &mut crate::state::RunTimer,
) {
    // Consume run_solver_requested before the match moves tectonic_state
    if ui_actions.run_solver_requested {
        ui_actions.run_solver_requested = false;
        if let Some(ref state) = tectonic_state {
            launch_solver(state, solver_config, bridge);
            run_timer.started_at = Some(std::time::Instant::now());
        }
    }

    egui::CollapsingHeader::new("Parameters").default_open(true).show(ui, |ui| {
        match *current_phase {
            PipelinePhase::Erosion => {
                draw_erosion(ui, erosion, isostasy_cache, upscale_cache, erosion_cache, generation);
            }
            PipelinePhase::Tectonics => {
                draw_tectonics(ui, tectonic_state, solver_config, bridge);
            }
            PipelinePhase::Isostasy => {
                draw_isostasy(ui, isostasy_params, isostasy_cache, ui_actions, centering);
            }
            PipelinePhase::Climate => draw_climate(ui, climate),
            PipelinePhase::UpscaleFbm => {
                draw_fbm(ui, fbm_params, isostasy_cache, upscale_cache, generation);
            }
            PipelinePhase::Hydrology => {
                draw_hydrology(
                    ui,
                    isostasy_cache,
                    upscale_cache,
                    erosion_cache,
                    flow_cache,
                    lake_cache,
                    generation,
                );
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

    // ── Consume top-bar run actions ──
    consume_top_bar_actions(
        ui_actions,
        solver_config,
        bridge,
        isostasy_cache,
        fbm_params,
        upscale_cache,
        erosion,
        erosion_cache,
        flow_cache,
        generation,
        run_timer,
    );
}

#[allow(clippy::too_many_arguments)]
fn consume_top_bar_actions(
    ui_actions: &mut crate::state::UiActions,
    _solver_config: &mut SolverConfig,
    _bridge: &mut SolverBridge,
    isostasy_cache: &IsostasyCache,
    fbm_params: &mut FbmParams,
    upscale_cache: &mut UpscaleCache,
    erosion_params: &mut ErosionParams,
    erosion_cache: &mut ErosionCache,
    flow_cache: &mut FlowCache,
    generation: &mut GenerationParamsUi,
    run_timer: &mut crate::state::RunTimer,
) {
    let mut started = false;
    if ui_actions.run_fbm_requested {
        ui_actions.run_fbm_requested = false;
        launch_fbm(fbm_params, isostasy_cache, upscale_cache, generation);
        started = true;
    }
    if ui_actions.run_erosion_requested {
        ui_actions.run_erosion_requested = false;
        launch_erosion(erosion_params, isostasy_cache, erosion_cache, generation);
        started = true;
    }
    if ui_actions.run_hydrology_requested {
        ui_actions.run_hydrology_requested = false;
        flow_cache.state = FlowState::Running;
        flow_cache.pending_config = Some(ymir_core::terrain::flow::FlowConfig {
            sea_level: isostasy_cache.sea_level_normalized,
        });
        started = true;
    }
    if started {
        run_timer.started_at = Some(std::time::Instant::now());
    }
}

// ── Tectonics panel ───────────────────────────────────────────────────────

fn draw_tectonics(
    ui: &mut egui::Ui,
    tectonic_state: Option<ResMut<TectonicState>>,
    solver_config: &mut SolverConfig,
    bridge: &mut SolverBridge,
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
    slider_row_f64(ui, "Basal friction", &mut solver_config.basal_friction, 0.0..=1.0);

    ui.checkbox(&mut solver_config.recycling.enabled, "Conservative recycling");
    if solver_config.recycling.enabled {
        slider_row_f64(ui, "Arc fraction", &mut solver_config.recycling.arc_fraction, 0.0..=0.5);
        slider_row_f64(ui, "Loss fraction", &mut solver_config.recycling.loss_fraction, 0.0..=0.3);
        ui.horizontal(|ui| {
            ui.label("Mantle delay");
            let mut delay = solver_config.recycling.mantle_delay as f32;
            if ui.add(egui::Slider::new(&mut delay, 1.0..=100.0)).changed() {
                solver_config.recycling.mantle_delay = delay as usize;
            }
        });
    }

    ui.checkbox(&mut solver_config.mantle.enabled, "Mantle flow");
    if solver_config.mantle.enabled {
        slider_row_f64(ui, "Mantle amplitude", &mut solver_config.mantle.amplitude, 0.0..=3.0);
        slider_row_f64(ui, "Mantle coupling", &mut solver_config.mantle.coupling, 0.0..=5.0);
    }

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
            "Ocean restore. Thr",
            &mut solver_config.boundaries.oceanic_restore_threshold,
            0.2..=0.6,
        );
        slider_row_f64(
            ui,
            "Ocean restore",
            &mut solver_config.boundaries.oceanic_restore_rate,
            0.0..=0.5,
        );

        ui.separator();
        ui.strong("Continental restore");
        slider_row_f64(
            ui,
            "Cont. min S",
            &mut solver_config.boundaries.continental_min_thickness,
            0.3..=0.8,
        );
        slider_row_f64(
            ui,
            "Cont. restore Thr",
            &mut solver_config.boundaries.continental_restore_threshold,
            0.05..=0.3,
        );
        slider_row_f64(
            ui,
            "Cont. restore",
            &mut solver_config.boundaries.continental_restore_rate,
            0.0..=0.1,
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
    ui.separator();
    ui.strong("Cratonic rigidity");
    ui.checkbox(&mut solver_config.cratonic.enabled, "Enabled");
    if solver_config.cratonic.enabled {
        slider_row_f64(ui, "Max factor", &mut solver_config.cratonic.max_factor, 1.0..=20.0);
        slider_row_f64(ui, "Decay power", &mut solver_config.cratonic.decay_power, 0.5..=4.0);
    }

    ui.add_space(4.0);
    ui.separator();
    ui.strong("Plastic yielding");
    ui.checkbox(&mut solver_config.yielding.enabled, "Enabled");
    if solver_config.yielding.enabled {
        slider_row_f64(ui, "Yield stress", &mut solver_config.yielding.yield_stress, 1.0..=500.0);
        ui.checkbox(&mut solver_config.yielding.weakening_enabled, "Strain weakening");
        if solver_config.yielding.weakening_enabled {
            slider_row_f64(
                ui,
                "Weakening frac.",
                &mut solver_config.yielding.weakening_fraction,
                0.0..=0.9,
            );
            slider_row_f64(
                ui,
                "Weakening ε_ref",
                &mut solver_config.yielding.weakening_strain_ref,
                0.1..=10.0,
            );
        }
    }

    ui.add_space(4.0);

    // ── Solver status ──
    draw_solver_status(ui, &bridge.state);
}

// ── Isostasy panel ───────────────────────────────────────────────────────

fn draw_isostasy(
    ui: &mut egui::Ui,
    isostasy_params: &mut IsostasyParams,
    isostasy_cache: &IsostasyCache,
    ui_actions: &mut crate::state::UiActions,
    centering: &mut crate::state::CenteringState,
) {
    ui.strong("Isostasy");

    if ui.button("⊕ Center Map").clicked() {
        ui_actions.center_requested = true;
    }

    if centering.original_field.is_some() {
        ui.horizontal(|ui| {
            ui.label("Offset X");
            if ui
                .add(egui::DragValue::new(&mut centering.offset_x).speed(1).range(-64..=64))
                .changed()
            {
                ui_actions.center_offset_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Offset Y");
            if ui
                .add(egui::DragValue::new(&mut centering.offset_y).speed(1).range(-64..=64))
                .changed()
            {
                ui_actions.center_offset_changed = true;
            }
        });
    }

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
        SolverState::Idle => {}
        SolverState::Running { stats, .. } => {
            // Diagnostic stats only (progress bar is in the top bar)
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
        cratonic: solver_config.cratonic.clone(),
        yielding: solver_config.yielding.clone(),
        basal_friction: solver_config.basal_friction,
        mantle: solver_config.mantle.clone(),
        recycling: solver_config.recycling.clone(),
    };

    let next_id = init.plates.len();
    let plate_ctx = DynamicPlateContext {
        ids: init.plate_ids.clone(),
        plates: init.plates.clone(),
        traction,
        next_id,
        disp_x: ymir_core::tectonics::solver::field::Field2D::new(grid_size),
        disp_y: ymir_core::tectonics::solver::field::Field2D::new(grid_size),
    };

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

#[allow(clippy::too_many_arguments)]
fn draw_erosion(
    ui: &mut egui::Ui,
    p: &mut ErosionParams,
    _isostasy_cache: &IsostasyCache,
    _upscale_cache: &UpscaleCache,
    erosion_cache: &mut ErosionCache,
    _generation: &GenerationParamsUi,
) {
    ui.strong("Hydraulic erosion");
    slider_row(ui, "Erosion rate", &mut p.erosion_rate, 0.0..=1.0);
    slider_row(ui, "Deposition rate", &mut p.deposition_rate, 0.0..=1.0);
    slider_row(ui, "Inertia", &mut p.inertia, 0.0..=0.5);
    slider_row(ui, "Gravity", &mut p.gravity, 1.0..=20.0);
    slider_row(ui, "Evaporation", &mut p.evaporation_rate, 0.001..=0.05);
    slider_row_u32(ui, "Max lifetime", &mut p.max_lifetime, 50..=300);
    slider_row(ui, "Droplets (M)", &mut p.droplets_millions, 0.5..=20.0);
    slider_row_u32(ui, "Coastal dep.", &mut p.coastal_deposition, 0..=30);
    slider_row(ui, "Min slope", &mut p.min_slope, 0.001..=0.1);

    ui.add_space(6.0);

    // Completion stats only (progress is shown in the top bar)
    if let ErosionState::Completed { elapsed } = &erosion_cache.state {
        ui.small(format!("🟢 Done in {:.1}s", elapsed.as_secs_f64()));
        if let Some(ref stats) = erosion_cache.stats {
            ui.small(format!(
                "Eroded: {:.1}  Deposited: {:.1}  Avg life: {:.0}",
                stats.total_eroded, stats.total_deposited, stats.avg_lifetime
            ));
        }
    }
}

// ── Hydrology panel ──────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn draw_hydrology(
    ui: &mut egui::Ui,
    _isostasy_cache: &IsostasyCache,
    _upscale_cache: &UpscaleCache,
    _erosion_cache: &mut ErosionCache,
    flow_cache: &mut FlowCache,
    lake_cache: &mut LakeCache,
    _generation: &GenerationParamsUi,
) {
    ui.strong("Flow & Rivers");

    if let FlowState::Completed { elapsed } = &flow_cache.state {
        ui.small(format!("🟢 Done in {:.2}s", elapsed.as_secs_f64()));
        if let Some(ref result) = flow_cache.result {
            ui.small(format!("{} basins", result.num_basins));
        }
        if let Some(ref rivers) = flow_cache.rivers {
            ui.small(format!("{} river segments", rivers.segments.len()));
        }
    }

    if flow_cache.result.is_some() {
        let rc = &mut flow_cache.river_config;
        let mut changed = false;

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.label("Stream thr.");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.monospace(format!("{:.0}", rc.stream_threshold));
            });
        });
        if ui
            .add(
                egui::Slider::new(&mut rc.stream_threshold, 100.0..=5000.0)
                    .show_value(false)
                    .logarithmic(true),
            )
            .changed()
        {
            changed = true;
        }

        ui.horizontal(|ui| {
            ui.label("River thr.");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.monospace(format!("{:.0}", rc.river_threshold));
            });
        });
        if ui
            .add(
                egui::Slider::new(&mut rc.river_threshold, 500.0..=20000.0)
                    .show_value(false)
                    .logarithmic(true),
            )
            .changed()
        {
            changed = true;
        }

        ui.horizontal(|ui| {
            ui.label("Major thr.");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.monospace(format!("{:.0}", rc.major_river_threshold));
            });
        });
        if ui
            .add(
                egui::Slider::new(&mut rc.major_river_threshold, 2000.0..=100000.0)
                    .show_value(false)
                    .logarithmic(true),
            )
            .changed()
        {
            changed = true;
        }

        if changed {
            flow_cache.rivers_dirty = true;
        }
    }

    // ── Lakes ──
    if flow_cache.result.is_some() {
        ui.add_space(4.0);
        ui.separator();
        ui.strong("Lakes");

        let lc = &mut lake_cache.config;
        let mut lake_changed = false;

        ui.horizontal(|ui| {
            ui.label("Min depth");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.monospace(format!("{:.4}", lc.min_depth));
            });
        });
        if ui
            .add(
                egui::Slider::new(&mut lc.min_depth, 0.0001..=0.01)
                    .show_value(false)
                    .logarithmic(true),
            )
            .changed()
        {
            lake_changed = true;
        }

        slider_row_usize(ui, "Min area", &mut lc.min_area, 5..=200);
        if ui.ctx().input(|i| i.pointer.any_released()) {
            lake_changed = true;
        }

        if lake_changed {
            lake_cache.dirty = true;
        }

        if let Some(ref result) = lake_cache.result {
            let total_area: usize = result.lakes.iter().map(|l| l.area).sum();
            let n = result.width * result.height;
            ui.small(format!(
                "{} lakes, total area {:.1}%",
                result.lakes.len(),
                total_area as f64 / n as f64 * 100.0
            ));
        }
    }
}

fn launch_erosion(
    params: &ErosionParams,
    isostasy_cache: &IsostasyCache,
    erosion_cache: &mut ErosionCache,
    generation: &GenerationParamsUi,
) {
    use ymir_core::erosion::hydraulic::ErosionConfig;

    // Clear previous result so the view falls back to the upscale heightmap
    // until the first erosion snapshot arrives.
    erosion_cache.heightmap = None;
    erosion_cache.sediment = None;
    erosion_cache.stats = None;
    erosion_cache.state = ErosionState::Running { completed: 0, total: 0 };

    let config = ErosionConfig {
        num_droplets: (params.droplets_millions * 1_000_000.0) as usize,
        deposition_rate: params.deposition_rate,
        erosion_rate: params.erosion_rate,
        inertia: params.inertia,
        gravity: params.gravity,
        evaporation_rate: params.evaporation_rate,
        max_lifetime: params.max_lifetime as usize,
        min_slope: params.min_slope,
        erosion_radius: 3,
        coastal_deposition_range: params.coastal_deposition as usize,
        sea_level: isostasy_cache.sea_level_normalized,
        batch_size: 50_000,
        reference_size: 256,
    };

    let seed = WorldSeed::new(generation.seed);

    erosion_cache.pending_config = Some(config);
    erosion_cache.pending_seed = Some(seed);
}

#[allow(clippy::too_many_arguments)]
fn draw_fbm(
    ui: &mut egui::Ui,
    p: &mut FbmParams,
    _isostasy_cache: &IsostasyCache,
    upscale_cache: &mut UpscaleCache,
    _generation: &GenerationParamsUi,
) {
    ui.strong("FBM Upscaling");

    let sizes = [256usize, 512, 1024, 2048];
    ui.horizontal(|ui| {
        ui.label("Target size");
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            for &sz in &sizes {
                if ui
                    .add(egui::Button::new(format!("{sz}")).selected(p.target_size == sz))
                    .clicked()
                {
                    p.target_size = sz;
                }
            }
        });
    });

    slider_row_usize(ui, "Octaves", &mut p.octaves, 1..=10);
    slider_row_f64(ui, "Lacunarity", &mut p.lacunarity, 1.5..=3.0);
    slider_row_f64(ui, "Persistence", &mut p.persistence, 0.3..=0.7);
    slider_row_f64(ui, "Amplitude", &mut p.amplitude_base, 0.01..=0.3);
    slider_row_f64(ui, "Slope factor", &mut p.amplitude_slope_factor, 0.0..=10.0);
    slider_row_f64(ui, "Anisotropy", &mut p.max_anisotropy, 1.0..=5.0);
    slider_row_f64(ui, "Sub. damping", &mut p.submarine_damping, 0.0..=1.0);

    ui.add_space(4.0);
    ui.checkbox(&mut p.domain_warp_enabled, "Domain warp");
    if p.domain_warp_enabled {
        slider_row_f64(ui, "Warp strength", &mut p.domain_warp_strength, 0.0..=1.0);
        slider_row_f64(ui, "Warp frequency", &mut p.domain_warp_frequency, 0.1..=2.0);
        slider_row_usize(ui, "Warp octaves", &mut p.domain_warp_octaves, 1..=5);
    }

    ui.add_space(6.0);

    // Completion info only (progress is in the top bar)
    if let FbmState::Completed { elapsed } = &upscale_cache.state {
        ui.small(format!("🟢 Done in {:.2}s", elapsed.as_secs_f64()));
        if let Some(ref hm) = upscale_cache.heightmap {
            ui.small(format!("Output: {}×{}", hm.width, hm.height));
        }
    }
}

fn launch_fbm(
    params: &FbmParams,
    isostasy_cache: &IsostasyCache,
    upscale_cache: &mut UpscaleCache,
    generation: &GenerationParamsUi,
) {
    // Build the coarse heightmap from the isostasy result.
    // The isostasy system stores the result in the texture, but we need the
    // actual GridF32. We re-run isostasy to get it. However, the isostasy
    // cache only stores stats, not the heightmap. We need to get it from the
    // terrain display's s_field via the bridge. Instead, we'll send the command
    // and let the thread handle it.
    //
    // Actually, we need the isostasy heightmap. The simplest path: the solver
    // thread receives the s_field (thickness), re-runs isostasy, then upscales.
    // But that couples the commands. Instead, let's grab the altitude from the
    // existing isostasy computation — but it's not stored as a resource.
    //
    // For now: we'll recompute isostasy in the thread from the thickness field.
    // This is cheap (~ms) compared to the FBM upscale (~seconds).
    //
    // HOWEVER, we don't have access to the s_field here. Let's use UiActions
    // pattern instead, or pass it differently.
    //
    // The cleanest approach: store the isostasy heightmap in IsostasyCache.
    // But that's a bigger refactor. For now, we'll send the coarse heightmap
    // via a different path.
    //
    // SIMPLEST: we already have the isostasy system that recomputes every frame.
    // Let's store the heightmap in UpscaleCache as the "source" when isostasy
    // updates, OR we add the heightmap to IsostasyCache.
    //
    // Let me just add the heightmap to IsostasyCache — it's the right thing.

    // For now this is wired as a flag; the actual command is sent from a system
    // that has access to the terrain display.
    upscale_cache.state = FbmState::Running;

    let config = ymir_core::terrain::upscale::FbmUpscaleConfig {
        target_size: params.target_size,
        octaves: params.octaves,
        lacunarity: params.lacunarity,
        persistence: params.persistence,
        amplitude_base: params.amplitude_base,
        amplitude_slope_factor: params.amplitude_slope_factor,
        max_anisotropy: params.max_anisotropy,
        submarine_damping: params.submarine_damping,
        base_frequency: 1.0,
        domain_warp_strength: if params.domain_warp_enabled {
            params.domain_warp_strength
        } else {
            0.0
        },
        domain_warp_frequency: params.domain_warp_frequency,
        domain_warp_octaves: params.domain_warp_octaves,
    };

    let seed = WorldSeed::new(generation.seed);

    // Store pending command data in the cache for the system to pick up
    upscale_cache.pending_config = Some(config);
    upscale_cache.pending_seed = Some(seed);
    upscale_cache.pending_sea_level = Some(isostasy_cache.sea_level_normalized);
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
