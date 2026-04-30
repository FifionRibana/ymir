//! Step 8.6 Phase 3 — parameter panel for the v2 bridge.
//!
//! Mirrors the legacy `parameter_panel::draw` shape but exposes v2
//! nondimensional knobs (Bi, Br, Mf, Cr, K, B_factor) and the
//! Step 9/10 toggles (cratonic, age field, slab, mantle, linear
//! solver). Drives `V2SolverBridge` directly: Run / Cancel buttons
//! submit / signal the worker thread.
//!
//! After Phase 8 sunset (legacy removal) this becomes the *only*
//! parameter panel in the viz; until then it coexists with
//! `parameter_panel.rs` (legacy) under bridge-mode dispatch.

use bevy::prelude::*;
use bevy_egui::egui;

use crate::bridge::v2::{
    presets, V2AgeFieldSpec, V2CratonicSpec, V2ForceKind, V2InitModeSpec, V2LinearSolverSpec,
    V2MantleSpec, V2RunSpec, V2RunState, V2SolverBridge,
};
use crate::visualization::v2_viz::{V2Field, V2VizState};

/// UI-side mutable copy of the run spec the user is editing. Submitted
/// to the bridge on Run click.
#[derive(Resource, Clone, Debug)]
pub struct V2EditableSpec(pub V2RunSpec);

impl Default for V2EditableSpec {
    fn default() -> Self {
        V2EditableSpec(V2RunSpec::active_medley_defaults())
    }
}

pub fn draw(
    ui: &mut egui::Ui,
    spec_state: &mut V2EditableSpec,
    bridge: &mut V2SolverBridge,
    viz: &mut V2VizState,
) {
    let spec = &mut spec_state.0;

    ui.heading("Tectonics v2");
    ui.add_space(4.0);

    // ── Run status badge ────────────────────────────────────────────
    let status_label = match &bridge.state {
        V2RunState::Idle => "Idle".to_string(),
        V2RunState::Running { spec: s, step, total, started_at, .. } => {
            let elapsed = started_at
                .as_ref()
                .map(|t| t.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            format!(
                "Running {}/{} ({}² × {} cfg) — {:.1}s",
                step, total, s.grid_nx, s.steps, elapsed
            )
        }
        V2RunState::Completed { elapsed, metrics, .. } => format!(
            "Completed in {:.1}s — peak|v|={:.2e}, CG mean={:.0}",
            elapsed.as_secs_f64(),
            metrics.vmax_peak,
            metrics.cg_iter_mean
        ),
        V2RunState::Imported { exported_at, scalar_metrics, .. } => format!(
            "Imported (exported {}) — peak|v|={:.2e}, CG mean={:.0}",
            exported_at, scalar_metrics.vmax_peak, scalar_metrics.cg_iter_mean
        ),
        V2RunState::Failed { error } => format!("Failed: {}", error),
    };
    ui.colored_label(
        match &bridge.state {
            V2RunState::Failed { .. } => egui::Color32::RED,
            V2RunState::Running { .. } => egui::Color32::YELLOW,
            V2RunState::Completed { .. } => egui::Color32::GREEN,
            V2RunState::Imported { .. } => egui::Color32::LIGHT_BLUE,
            V2RunState::Idle => egui::Color32::GRAY,
        },
        status_label,
    );

    ui.add_space(8.0);
    ui.separator();

    // ── Preset dropdown (Phase 4) ──────────────────────────────────
    ui.horizontal(|ui| {
        ui.label("Preset:");
        egui::ComboBox::from_id_salt("v2_preset_dropdown")
            .selected_text(spec.preset_label.clone())
            .show_ui(ui, |ui| {
                for name in presets::list() {
                    if ui
                        .selectable_label(spec.preset_label == name, name)
                        .clicked()
                    {
                        match presets::load(name) {
                            Ok(loaded) => {
                                // Preserve the user's current
                                // output_dir choice — preset files
                                // always default to a fresh temp
                                // path which would shadow the
                                // user's intent.
                                let prev_dir = spec.output_dir.clone();
                                let prev_capture = spec.capture_endpoints;
                                *spec = loaded;
                                spec.output_dir = prev_dir;
                                spec.capture_endpoints = prev_capture;
                            }
                            Err(e) => {
                                eprintln!("[ymir-viz] preset load failed: {}", e);
                            }
                        }
                    }
                }
            });
    });
    ui.add_space(4.0);

    // ── Grid + seed + plates ────────────────────────────────────────
    egui::CollapsingHeader::new("Grid & seed")
        .default_open(true)
        .show(ui, |ui| {
            egui::ComboBox::from_label("Resolution")
                .selected_text(format!("{}²", spec.grid_nx))
                .show_ui(ui, |ui| {
                    for &res in &[32usize, 64, 128] {
                        if ui.selectable_label(spec.grid_nx == res, format!("{}²", res)).clicked() {
                            spec.grid_nx = res;
                            spec.grid_ny = res;
                        }
                    }
                });
            ui.add(egui::DragValue::new(&mut spec.seed).prefix("seed = "));
            ui.add(
                egui::Slider::new(&mut spec.num_plates, 3..=15)
                    .text("Voronoï plates"),
            );
            ui.add(
                egui::Slider::new(&mut spec.continental_ratio, 0.1..=0.7)
                    .text("Continental ratio"),
            );
            ui.add(
                egui::Slider::new(&mut spec.steps, 10..=500)
                    .text("Steps"),
            );
        });

    // ── Yielding + drag ────────────────────────────────────────────
    egui::CollapsingHeader::new("Yielding & drag")
        .default_open(true)
        .show(ui, |ui| {
            ui.add(
                egui::Slider::new(&mut spec.bi, 0.0..=0.5)
                    .text("Bi (yield)")
                    .step_by(0.005),
            );
            ui.add(
                egui::Slider::new(&mut spec.br, 0.0..=0.5)
                    .text("Br (basal drag)")
                    .step_by(0.005),
            );
        });

    // ── Initialisation (Phase 8a/8d) ───────────────────────────────
    egui::CollapsingHeader::new("Initialisation (S̃)")
        .default_open(true)
        .show(ui, |ui| {
            init_mode_widget(ui, &mut spec.init_mode);
        });

    // ── Mantle ─────────────────────────────────────────────────────
    egui::CollapsingHeader::new("Mantle")
        .default_open(true)
        .show(ui, |ui| {
            let mut mantle_on = matches!(spec.mantle, V2MantleSpec::On { .. });
            if ui.checkbox(&mut mantle_on, "Enabled").changed() {
                spec.mantle = if mantle_on {
                    V2MantleSpec::default()
                } else {
                    V2MantleSpec::Off
                };
            }
            if let V2MantleSpec::On {
                mf,
                num_modes,
                seed: m_seed,
                evolution_rate,
                ..
            } = &mut spec.mantle
            {
                ui.add(egui::Slider::new(mf, 0.0..=2.0).text("Mf").step_by(0.05));
                ui.add(egui::Slider::new(num_modes, 1..=12).text("modes"));
                ui.add(egui::DragValue::new(m_seed).prefix("mantle seed = "));
                // Phase 8d — exposed for calibration of mode evolution
                // through a run. `0.0` keeps the mantle pattern static
                // (Step 8 baseline default); >0 advects the modes at
                // that rate per nondim time unit.
                ui.add(
                    egui::Slider::new(evolution_rate, 0.0..=2.0)
                        .text("evolution_rate")
                        .step_by(0.05),
                );
            }
        });

    // ── Cratonic immunity ──────────────────────────────────────────
    egui::CollapsingHeader::new("Cratonic immunity")
        .default_open(true)
        .show(ui, |ui| {
            let mut crat_on = matches!(spec.cratonic, V2CratonicSpec::On { .. });
            if ui.checkbox(&mut crat_on, "Enabled").changed() {
                spec.cratonic = if crat_on {
                    V2CratonicSpec::default()
                } else {
                    V2CratonicSpec::Off
                };
            }
            if let V2CratonicSpec::On {
                cr,
                k_viscous,
                b_factor,
                smoothing_width,
                plate_area_min,
            } = &mut spec.cratonic
            {
                ui.add(egui::Slider::new(cr, 0.0..=1.0).text("Cr"));
                ui.add(egui::Slider::new(k_viscous, 1.0..=20.0).text("K (viscous)"));
                ui.add(egui::Slider::new(b_factor, 1.0..=20.0).text("B_factor (Bi mult)"));
                // Phase 8d — geometry knobs previously hidden behind
                // `CratonicConfigEnabled::default()`.
                ui.add(
                    egui::Slider::new(smoothing_width, 0.02..=0.20)
                        .text("smoothing_width")
                        .step_by(0.005),
                );
                ui.add(
                    egui::Slider::new(plate_area_min, 0.05..=0.20)
                        .text("plate_area_min")
                        .step_by(0.005),
                );
            }
        });

    // ── Age field ──────────────────────────────────────────────────
    egui::CollapsingHeader::new("Age field")
        .default_open(false)
        .show(ui, |ui| {
            let mut age_on = matches!(spec.age_field, V2AgeFieldSpec::On { .. });
            if ui.checkbox(&mut age_on, "Enabled").changed() {
                spec.age_field = if age_on {
                    V2AgeFieldSpec::default()
                } else {
                    V2AgeFieldSpec::Off
                };
            }
            if let V2AgeFieldSpec::On { continental_age_init, oceanic_age_init } =
                &mut spec.age_field
            {
                ui.add(
                    egui::Slider::new(continental_age_init, 0.0..=20.0)
                        .text("continental init"),
                );
                ui.add(
                    egui::Slider::new(oceanic_age_init, 0.0..=5.0)
                        .text("oceanic init"),
                );
            }
        });

    // ── Slab + linear solver + force ───────────────────────────────
    egui::CollapsingHeader::new("Solver options")
        .default_open(false)
        .show(ui, |ui| {
            ui.checkbox(
                &mut spec.slab_enabled,
                "Slab pull (forward-compat — currently a no-op pending §4.8 co-calibration)",
            );
            egui::ComboBox::from_label("Linear solver")
                .selected_text(match spec.linear_solver {
                    V2LinearSolverSpec::Jacobi => "JacobiCG",
                    V2LinearSolverSpec::Amg => "AmgCG",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut spec.linear_solver,
                        V2LinearSolverSpec::Jacobi,
                        "JacobiCG",
                    );
                    ui.selectable_value(
                        &mut spec.linear_solver,
                        V2LinearSolverSpec::Amg,
                        "AmgCG (Phase 5+)",
                    );
                });
            let mut sin_amp = match spec.force {
                V2ForceKind::Gpe => 0.0,
                V2ForceKind::Sinusoidal { amplitude } => amplitude,
            };
            let mut force_is_sin = matches!(spec.force, V2ForceKind::Sinusoidal { .. });
            if ui.checkbox(&mut force_is_sin, "Sinusoidal forcing (regression)").changed() {
                spec.force = if force_is_sin {
                    V2ForceKind::Sinusoidal { amplitude: 10.0 }
                } else {
                    V2ForceKind::Gpe
                };
            }
            if force_is_sin {
                if ui
                    .add(egui::Slider::new(&mut sin_amp, 0.0..=20.0).text("ε"))
                    .changed()
                {
                    spec.force = V2ForceKind::Sinusoidal { amplitude: sin_amp };
                }
            }
        });

    ui.add_space(8.0);
    ui.separator();

    // ── Field selector (D5 dropdown) + legend ──────────────────────
    ui.heading("Display");
    egui::ComboBox::from_id_salt("v2_field_dropdown")
        .selected_text(viz.field.label())
        .show_ui(ui, |ui| {
            for &f in V2Field::ALL {
                if ui.selectable_label(viz.field == f, f.label()).clicked() {
                    viz.field = f;
                }
            }
        });
    ui.label(egui::RichText::new(viz.field.legend_caption()).small());
    draw_legend_bar(ui, viz.field);

    // Phase 8b — Voronoï + velocity overlays. Default off so Phase 7
    // screenshots remain unchanged unless the user opts in. Toggling
    // either flag invalidates `viz.last_signature` (the
    // `overlay_bits` component), forcing a re-render on the next
    // frame.
    ui.add_space(4.0);
    ui.checkbox(&mut viz.show_voronoi_boundaries, "Show Voronoï boundaries");
    ui.checkbox(&mut viz.show_velocity_vectors, "Show velocity vectors");

    ui.add_space(8.0);

    // ── Run / Cancel ────────────────────────────────────────────────
    let is_running = matches!(bridge.state, V2RunState::Running { .. });
    ui.horizontal(|ui| {
        if ui
            .add_enabled(!is_running, egui::Button::new("\u{25b6} Run"))
            .clicked()
        {
            // Stamp a per-run output dir under the OS temp so each
            // click gets its own subdir if capture is enabled (Phase 6).
            let stamp = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0);
            spec.output_dir = std::env::temp_dir().join(format!("ymir_v2_run_{}", stamp));
            if let Err(e) = bridge.submit_run(spec.clone()) {
                eprintln!("[ymir-viz] failed to submit v2 run: {}", e);
            }
        }
        if ui
            .add_enabled(is_running, egui::Button::new("\u{23f9} Cancel"))
            .clicked()
        {
            bridge.request_cancel();
        }
        let can_capture = matches!(
            bridge.state,
            V2RunState::Completed { .. } | V2RunState::Imported { .. }
        );
        if ui
            .add_enabled(can_capture, egui::Button::new("\u{1f4f7} Capture"))
            .clicked()
        {
            viz.capture_requested = true;
        }
    });

    // Phase 6 — last screenshot status line
    if let Some(result) = &viz.last_capture {
        match result {
            Ok(path) => {
                ui.label(
                    egui::RichText::new(format!("Saved: {}", path.display()))
                        .small()
                        .color(egui::Color32::LIGHT_GREEN),
                );
            }
            Err(err) => {
                ui.label(
                    egui::RichText::new(format!("Capture failed: {}", err))
                        .small()
                        .color(egui::Color32::LIGHT_RED),
                );
            }
        }
    }

    // ── Phase 8e — Export / Import ─────────────────────────────────
    ui.add_space(8.0);
    ui.separator();
    ui.label(egui::RichText::new("Snapshot (Phase 8e)").strong());
    ui.horizontal(|ui| {
        let can_export = matches!(bridge.state, V2RunState::Completed { .. });
        if ui
            .add_enabled(can_export, egui::Button::new("\u{1f4be} Export run"))
            .on_hover_text("Save the current Completed run to JSON")
            .clicked()
        {
            viz.export_requested = true;
        }
    });
    if let Some(result) = &viz.last_export {
        match result {
            Ok(path) => ui.label(
                egui::RichText::new(format!("Exported to: {}", path.display()))
                    .small()
                    .color(egui::Color32::LIGHT_GREEN),
            ),
            Err(err) => ui.label(
                egui::RichText::new(format!("Export failed: {}", err))
                    .small()
                    .color(egui::Color32::LIGHT_RED),
            ),
        };
    }

    ui.add_space(4.0);
    ui.label("Import path:");
    ui.add(
        egui::TextEdit::singleline(&mut viz.import_path_buffer)
            .hint_text("/path/to/snapshot.json")
            .desired_width(f32::INFINITY),
    );
    if ui.button("\u{1f4c2} Import run").clicked() {
        let trimmed = viz.import_path_buffer.trim();
        if trimmed.is_empty() {
            viz.last_import = Some(Err("path is empty".to_string()));
        } else {
            viz.import_requested_path = Some(std::path::PathBuf::from(trimmed));
        }
    }
    if let Some(result) = &viz.last_import {
        match result {
            Ok(path) => ui.label(
                egui::RichText::new(format!("Imported from: {}", path.display()))
                    .small()
                    .color(egui::Color32::LIGHT_GREEN),
            ),
            Err(err) => ui.label(
                egui::RichText::new(format!("Import failed: {}", err))
                    .small()
                    .color(egui::Color32::LIGHT_RED),
            ),
        };
    }
}

/// Phase 8d — `InitMode` editor widget. Dropdown picks the variant;
/// the sub-block below adapts to the selection and exposes per-mode
/// numeric parameters with sensible slider ranges.
///
/// Switching variants resets the inner numeric payload to that
/// variant's defaults so the user does not see stale values from the
/// previous selection.
fn init_mode_widget(ui: &mut egui::Ui, mode: &mut V2InitModeSpec) {
    let current_idx = mode.variant_index();
    egui::ComboBox::from_id_salt("v2_init_mode_dropdown")
        .selected_text(mode.ui_label())
        .show_ui(ui, |ui| {
            if ui
                .selectable_label(current_idx == 0, "Checkerboard (legacy sinusoidal)")
                .clicked()
                && current_idx != 0
            {
                *mode = V2InitModeSpec::Checkerboard;
            }
            if ui
                .selectable_label(current_idx == 1, "Uniform (TDD §4.2 default)")
                .clicked()
                && current_idx != 1
            {
                *mode = V2InitModeSpec::Uniform { boundary_smoothing_width: 1.0 };
            }
            if ui
                .selectable_label(current_idx == 2, "Gaussian (peak at centroid)")
                .clicked()
                && current_idx != 2
            {
                *mode = V2InitModeSpec::Gaussian {
                    sigma_continental: 5.0,
                    sigma_oceanic: 5.0,
                };
            }
            if ui
                .selectable_label(current_idx == 3, "Convolution (Gaussian blur)")
                .clicked()
                && current_idx != 3
            {
                *mode = V2InitModeSpec::Convolution { sigma: 2.0 };
            }
        });

    match mode {
        V2InitModeSpec::Checkerboard => {
            ui.label(
                egui::RichText::new(
                    "Pre-Phase-8a sinusoidal perturbation. Required for \
                     Steps 0–10 numerical regression.",
                )
                .small()
                .weak(),
            );
        }
        V2InitModeSpec::Uniform { boundary_smoothing_width } => {
            ui.add(
                egui::Slider::new(boundary_smoothing_width, 0.0..=8.0)
                    .text("boundary smoothing width (cells)")
                    .step_by(0.25),
            );
            ui.label(
                egui::RichText::new(
                    "Flat per-plate-type. Smoothstep blending across \
                     inter-plate edges over `width` cells.",
                )
                .small()
                .weak(),
            );
        }
        V2InitModeSpec::Gaussian { sigma_continental, sigma_oceanic } => {
            ui.add(
                egui::Slider::new(sigma_continental, 0.5..=20.0)
                    .text("σ continental (cells)")
                    .step_by(0.25),
            );
            ui.add(
                egui::Slider::new(sigma_oceanic, 0.5..=20.0)
                    .text("σ oceanic (cells)")
                    .step_by(0.25),
            );
            ui.label(
                egui::RichText::new(
                    "Per-plate Gaussian, peaked at each Voronoï seed and \
                     decaying with periodic minimum-image distance.",
                )
                .small()
                .weak(),
            );
        }
        V2InitModeSpec::Convolution { sigma } => {
            ui.add(
                egui::Slider::new(sigma, 0.5..=8.0)
                    .text("σ (cells)")
                    .step_by(0.1),
            );
            ui.label(
                egui::RichText::new(
                    "Periodic Gaussian blur of the binary plate-type \
                     mask. Output stays inside [oceanic, continental].",
                )
                .small()
                .weak(),
            );
        }
    }
}

/// Phase 5 D5 — colorbar for the currently displayed field. The bar
/// is allocated as a 200×16 rectangle and painted with 32 sample
/// stops drawn from the matching colormap. This is decorative; the
/// raster sprite itself is the authoritative legend.
fn draw_legend_bar(ui: &mut egui::Ui, field: V2Field) {
    use crate::visualization::colormap::{
        age_colormap, cratonic_grayscale, hypsometric_colormap, log_hot,
    };
    let (rect, _) = ui.allocate_exact_size(egui::vec2(200.0, 16.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    const SAMPLES: usize = 32;
    let cell_w = rect.width() / SAMPLES as f32;
    for k in 0..SAMPLES {
        let t = k as f64 / (SAMPLES - 1) as f64;
        let rgba = match field {
            V2Field::SThickness => hypsometric_colormap(t),
            V2Field::Age => age_colormap(t),
            V2Field::Cratonic => cratonic_grayscale(t),
            V2Field::StrainRate | V2Field::VelocityMagnitude => log_hot(t),
        };
        let x0 = rect.left() + k as f32 * cell_w;
        let cell_rect =
            egui::Rect::from_min_size(egui::pos2(x0, rect.top()), egui::vec2(cell_w + 1.0, rect.height()));
        painter.rect_filled(
            cell_rect,
            0.0,
            egui::Color32::from_rgba_unmultiplied(rgba[0], rgba[1], rgba[2], rgba[3]),
        );
    }
}
