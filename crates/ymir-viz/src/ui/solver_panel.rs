//! Solver control panel: parameters, run button, progress bar.

use bevy_egui::egui;

use crate::bridge::plugin::{SolverBridge, SolverState};
use crate::bridge::commands::SolverCommand;
use crate::visualization::render::TerrainDisplay;

use ymir_core::tectonics::solver::config::{NonlinearSolver, PicardConfig, TectonicsConfig, NewtonConfig};
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::traction::TractionField;

use std::sync::atomic::Ordering;

/// Persistent UI state for the solver panel.
pub struct SolverUiState {
    pub grid_size: usize,
    pub num_timesteps: usize,
    pub gravity_factor: f64,
    pub cfl_factor: f64,
    pub power_law_n_index: usize, // 0 = linear (1.0), 1 = realistic (3.0)
    pub picard_tolerance: f64,
    pub picard_relaxation: f64,
    pub plate_preset: usize, // 0 = convergent, 1 = divergent
    pub plate_speed: f64,
    pub nonlinear_solver_index: usize, // 0 = Picard, 1 = Newton
}

impl Default for SolverUiState {
    fn default() -> Self {
        Self {
            grid_size: 128,
            num_timesteps: 300,
            gravity_factor: 1.0,
            cfl_factor: 0.5,
            power_law_n_index: 1,
            picard_tolerance: 1e-4,
            picard_relaxation: 0.7,
            plate_preset: 0,
            plate_speed: 0.02,
            nonlinear_solver_index: 0,
        }
    }
}

pub fn draw(
    ui: &mut egui::Ui,
    state: &mut SolverUiState,
    bridge: &mut SolverBridge,
    terrain_display: &mut TerrainDisplay,
) {
    // --- Solver state ---
    ui.heading("Solver");
    ui.separator();

    match &bridge.state {
        SolverState::Idle => {
            ui.label("⚪ Ready");
        }
        SolverState::Running {
            step,
            total_steps,
            stats,
        } => {
            let frac = *step as f32 / (*total_steps).max(1) as f32;
            ui.label(format!("🟡 Running… step {}/{}", step, total_steps));
            ui.add(egui::ProgressBar::new(frac).show_percentage());
            if let Some(s) = stats {
                ui.small(format!(
                    "max_v={:.4}  S=[{:.3}, {:.3}]  picard={}",
                    s.max_velocity, s.min_thickness, s.max_thickness, s.picard_iterations
                ));
            }
        }
        SolverState::Completed { elapsed } => {
            ui.label(format!("🟢 Done in {:.1}s", elapsed.as_secs_f64()));
        }
        SolverState::Failed { error } => {
            ui.colored_label(egui::Color32::RED, format!("❌ {error}"));
        }
    }

    ui.separator();

    let is_running = matches!(bridge.state, SolverState::Running { .. });

    // --- Grid config ---
    ui.collapsing("Grid", |ui| {
        ui.horizontal(|ui| {
            ui.label("Size:");
            for &sz in &[64usize, 128, 256] {
                if ui
                    .selectable_label(state.grid_size == sz, format!("{sz}"))
                    .clicked()
                {
                    state.grid_size = sz;
                }
            }
        });
    });

    // --- Tectonic parameters ---
    ui.collapsing("Parameters", |ui| {
        ui.add(
            egui::Slider::new(&mut state.num_timesteps, 50..=1000).text("timesteps"),
        );
        ui.add(
            egui::Slider::new(&mut state.gravity_factor, 0.1..=5.0).text("gravity"),
        );
        ui.add(
            egui::Slider::new(&mut state.cfl_factor, 0.05..=0.9).text("CFL"),
        );
        ui.add(
            egui::Slider::new(&mut state.picard_relaxation, 0.3..=1.0).text("relaxation"),
        );

        ui.horizontal(|ui| {
            ui.label("Viscosity:");
            egui::ComboBox::from_id_salt("power_law")
                .selected_text(if state.power_law_n_index == 0 {
                    "Linear (n=1)"
                } else {
                    "Power-law (n=3)"
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut state.power_law_n_index, 0, "Linear (n=1)");
                    ui.selectable_value(&mut state.power_law_n_index, 1, "Power-law (n=3)");
                });
        });

        ui.horizontal(|ui| {
            ui.label("Solver:");
            egui::ComboBox::from_id_salt("nl_solver")
                .selected_text(if state.nonlinear_solver_index == 0 {
                    "Picard"
                } else {
                    "Newton (JFNK)"
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut state.nonlinear_solver_index, 0, "Picard");
                    ui.selectable_value(&mut state.nonlinear_solver_index, 1, "Newton (JFNK)");
                });
        });
    });

    // --- Plate config ---
    ui.collapsing("Plates", |ui| {
        ui.horizontal(|ui| {
            ui.label("Preset:");
            egui::ComboBox::from_id_salt("plate_preset")
                .selected_text(match state.plate_preset {
                    0 => "Convergent",
                    1 => "Divergent",
                    _ => "Unknown",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut state.plate_preset, 0, "Convergent");
                    ui.selectable_value(&mut state.plate_preset, 1, "Divergent");
                });
        });
        ui.add(
            egui::Slider::new(&mut state.plate_speed, 0.001..=0.1)
                .text("speed")
                .logarithmic(true),
        );
    });

    ui.separator();

    // --- Controls ---
    ui.horizontal(|ui| {
        if ui
            .add_enabled(!is_running, egui::Button::new("▶ Run"))
            .clicked()
        {
            launch_solver(state, bridge, terrain_display);
        }
        if ui
            .add_enabled(is_running, egui::Button::new("⏹ Cancel"))
            .clicked()
        {
            bridge.cancel_flag.store(true, Ordering::Relaxed);
        }
        if ui
            .add_enabled(!is_running, egui::Button::new("↺ Reset"))
            .clicked()
        {
            bridge.state = SolverState::Idle;
            // Reset terrain to uniform S=1
            let n = state.grid_size;
            let field = Field2D::filled(n, 1.0);
            terrain_display.update_field(field);
        }
    });
}

fn launch_solver(
    state: &SolverUiState,
    bridge: &mut SolverBridge,
    terrain_display: &mut TerrainDisplay,
) {
    let n = state.grid_size;
    let dx = 1.0 / n as f64;

    let power_law_n = if state.power_law_n_index == 0 {
        1.0
    } else {
        3.0
    };

    let nonlinear_solver = if state.nonlinear_solver_index == 0 {
        NonlinearSolver::Picard
    } else {
        NonlinearSolver::Newton
    };

    let config = TectonicsConfig {
        num_timesteps: state.num_timesteps,
        gravity_factor: state.gravity_factor,
        cfl_factor: state.cfl_factor,
        s_min: 0.1,
        s_max: 2.5,
        nonlinear_solver,
        picard: PicardConfig {
            max_iterations: 50,
            tolerance: state.picard_tolerance,
            relaxation: state.picard_relaxation,
            cg_max_iter: 500,
            cg_tolerance: 1e-8,
            strain_rate_min: 1e-6,
            power_law_n,
        },
        newton: NewtonConfig::default(),
    };

    let plates = match state.plate_preset {
        0 => TractionField::two_plates_convergent(n, state.plate_speed),
        1 => TractionField::two_plates_divergent(n, state.plate_speed),
        _ => TractionField::zero(n),
    };

    // Use current field if available, otherwise uniform S=1
    let initial_s = terrain_display
        .s_field
        .clone()
        .unwrap_or_else(|| Field2D::filled(n, 1.0));

    let _ = bridge.commands_tx.send(SolverCommand::RunTectonics {
        config,
        plates,
        initial_s,
        grid_size: n,
        dx,
    });

    bridge.state = SolverState::Running {
        step: 0,
        total_steps: state.num_timesteps,
        stats: None,
    };
}
