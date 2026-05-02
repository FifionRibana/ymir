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
use crate::phases;
use crate::pipeline::{ActivePhase, PipelinePhase};
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

#[allow(clippy::too_many_arguments)]
pub fn draw(
    ui: &mut egui::Ui,
    spec_state: &mut V2EditableSpec,
    bridge: &mut V2SolverBridge,
    viz: &mut V2VizState,
    active: ActivePhase,
    isostasy_params: &mut phases::isostasy::IsostasyParams,
    isostasy_cache: &phases::isostasy::IsostasyCache,
    fbm_params: &mut phases::upscale_fbm::FbmParams,
    fbm_cache: &phases::upscale_fbm::FbmCache,
    erosion_params: &mut phases::erosion::ErosionParams,
    erosion_cache: &phases::erosion::ErosionCache,
    hydrology_params: &mut phases::hydrology::HydrologyParams,
    hydrology_cache: &phases::hydrology::HydrologyCache,
) {
    let spec = &mut spec_state.0;

    ui.heading("Tectonics v2");
    ui.add_space(4.0);

    // ── Run status badge ────────────────────────────────────────────
    let is_preview = matches!(&bridge.state, V2RunState::Idle) && viz.preview.is_some();
    let status_label = match &bridge.state {
        V2RunState::Idle => {
            if is_preview {
                format!(
                    "Preview ({}) — {}² × {} steps · seed {}",
                    spec.preset_label, spec.grid_nx, spec.steps, spec.seed
                )
            } else {
                "Idle".to_string()
            }
        }
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
            V2RunState::Idle if is_preview => egui::Color32::from_rgb(0xC0, 0xA0, 0x60),
            V2RunState::Idle => egui::Color32::GRAY,
        },
        status_label,
    );

    ui.add_space(8.0);
    ui.separator();

    // ── Save / Load (common across phases) ─────────────────────────
    // Sits at the top of the panel — visible regardless of which
    // pipeline phase is active. Replaces the legacy text-input
    // import path with a list of `output/seed<S>_<R>/` directories
    // that contain a `snapshot.json`.
    draw_save_load_section(ui, viz, &bridge.state);

    ui.add_space(8.0);
    ui.separator();

    // ── V2 raster field selector (global — visible on every phase) ──
    // The V2Field dropdown drives the V2 raster painter
    // (`update_v2_texture`). It is meaningful on the Tectonics view
    // only, but the dropdown sits here so the user can pick a field
    // ahead of switching back to Tectonics — same shape as the
    // legacy left-toolbar field picker that lived above the phase
    // sub-views.
    ui.label(egui::RichText::new("V2 raster field").strong());
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

    ui.add_space(8.0);
    ui.separator();

    // ── V2 solver controls ─────────────────────────────────────────
    // Only show on the Tectonics phase. The other phases consume
    // the v2 final state but don't drive its parameters; hiding
    // these keeps the right panel focused on the active phase.
    let show_v2_controls = matches!(active.0, PipelinePhase::Tectonics);
    if show_v2_controls {
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
        // Step 8.6 follow-up — Continue button. Reuses the prior
        // run / imported snapshot's `final_state` as the start
        // state for a new run (S̃ + vx/vy + age + cratonic_factor),
        // letting the user "add 100 steps to the 100 already
        // simulated" or extend an imported run. The user-edited
        // `spec` provides the step count, total_time, and any
        // tweaked physics knobs; voronoi-relevant fields (seed,
        // num_plates, continental_ratio, grid dims) are
        // overridden from the source to keep the tessellation
        // consistent — otherwise the prior plate_id / plate_type
        // would no longer match the regenerated voronoi.
        let continue_source = match &bridge.state {
            V2RunState::Completed { spec: src, final_state, .. } => {
                Some((src.clone(), final_state.as_ref().clone()))
            }
            V2RunState::Imported { spec: src, final_state, .. } => {
                Some((src.clone(), final_state.as_ref().clone()))
            }
            _ => None,
        };
        let can_continue = continue_source.is_some();
        if ui
            .add_enabled(can_continue, egui::Button::new("\u{21bb} Continue"))
            .on_hover_text(
                "Continue from the prior run's final state \
                 (voronoi config locked from source; user spec \
                 provides additional steps + physics knobs).",
            )
            .clicked()
        {
            if let Some((source_spec, from_state)) = continue_source {
                let mut next_spec = spec.clone();
                next_spec.seed = source_spec.seed;
                next_spec.grid_nx = source_spec.grid_nx;
                next_spec.grid_ny = source_spec.grid_ny;
                next_spec.num_plates = source_spec.num_plates;
                next_spec.continental_ratio = source_spec.continental_ratio;
                next_spec.init_mode = source_spec.init_mode;
                let stamp = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map(|d| d.as_secs())
                    .unwrap_or(0);
                next_spec.output_dir =
                    std::env::temp_dir().join(format!("ymir_v2_continue_{}", stamp));
                if let Err(e) = bridge.submit_continue(next_spec, from_state) {
                    eprintln!("[ymir-viz] failed to submit continue run: {}", e);
                }
            }
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
    } // end `if show_v2_controls`

    // ── Active phase parameters ────────────────────────────────────
    // Each non-Tectonics phase exposes its own collapsible config
    // section here. Visible only when that phase is active so the
    // panel doesn't grow unbounded for the user.
    let v2_finished = matches!(
        bridge.state,
        V2RunState::Completed { .. } | V2RunState::Imported { .. }
    );
    match active.0 {
        PipelinePhase::Tectonics => {}
        PipelinePhase::Isostasy => {
            ui.add_space(8.0);
            ui.separator();
            phases::isostasy::draw_section(ui, isostasy_params, isostasy_cache, v2_finished);
        }
        PipelinePhase::UpscaleFbm => {
            ui.add_space(8.0);
            ui.separator();
            let iso_ready = isostasy_cache.result.is_some();
            phases::upscale_fbm::draw_section(ui, fbm_params, fbm_cache, iso_ready);
        }
        PipelinePhase::Erosion => {
            ui.add_space(8.0);
            ui.separator();
            let fbm_ready = fbm_cache.result.is_some();
            phases::erosion::draw_section(ui, erosion_params, erosion_cache, fbm_ready);
        }
        PipelinePhase::Hydrology => {
            ui.add_space(8.0);
            ui.separator();
            let erosion_ready = erosion_cache.result.is_some();
            phases::hydrology::draw_section(
                ui,
                hydrology_params,
                hydrology_cache,
                erosion_ready,
            );
        }
        PipelinePhase::Climate => {
            ui.add_space(8.0);
            ui.separator();
            phases::climate::draw_section(ui);
        }
        PipelinePhase::Biome => {
            ui.add_space(8.0);
            ui.separator();
            phases::biome::draw_section(ui);
        }
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
            V2Field::StrainRate | V2Field::VelocityMagnitude | V2Field::Slope => log_hot(t),
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

/// Step 8.6 follow-up — common Save / Load section. Sits at the top of
/// the right panel (above the v2 solver controls) and is visible on
/// every pipeline phase. The "Available runs" list scans
/// `output/` for sub-directories that contain a `snapshot.json` and
/// gives each one a Load button. The list is cached on
/// [`V2VizState::cached_run_dirs`] and invalidated by the Refresh
/// button or by a successful export (so a freshly written run shows up
/// without the user clicking Refresh).
fn draw_save_load_section(
    ui: &mut egui::Ui,
    viz: &mut V2VizState,
    bridge_state: &V2RunState,
) {
    ui.label(egui::RichText::new("Save / Load").strong());

    let can_export = matches!(
        bridge_state,
        V2RunState::Completed { .. } | V2RunState::Imported { .. }
    );
    let is_running = matches!(bridge_state, V2RunState::Running { .. });

    ui.horizontal(|ui| {
        if ui
            .add_enabled(can_export, egui::Button::new("\u{1f4be} Save current run"))
            .on_hover_text(
                "Save thickness + altitude + upscale + eroded + flow + \
                 lakes (whichever ran) plus snapshot.json to \
                 `output/seed<seed>_<resolution>/`.",
            )
            .clicked()
        {
            viz.export_requested = true;
        }
        if ui
            .button("\u{1f504} Refresh")
            .on_hover_text("Rescan output/ for available runs.")
            .clicked()
        {
            viz.cached_run_dirs = None;
        }
    });
    if let Some(result) = &viz.last_export {
        match result {
            Ok(path) => {
                ui.label(
                    egui::RichText::new(format!("Saved → {}", path.display()))
                        .small()
                        .color(egui::Color32::LIGHT_GREEN),
                );
            }
            Err(err) => {
                ui.label(
                    egui::RichText::new(format!("Save failed: {}", err))
                        .small()
                        .color(egui::Color32::LIGHT_RED),
                );
            }
        }
    }

    ui.add_space(4.0);
    ui.label("Available runs:");
    let dirs = viz
        .cached_run_dirs
        .get_or_insert_with(|| list_export_dirs(std::path::Path::new("output")));
    if dirs.is_empty() {
        ui.label(
            egui::RichText::new("(no saved runs in output/)")
                .small()
                .weak(),
        );
    } else {
        let dirs_clone = dirs.clone();
        for dir in dirs_clone {
            let name = dir
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new(&name).small().monospace());
                if ui
                    .add_enabled(!is_running, egui::Button::new("Load").small())
                    .on_hover_text("Replace the current bridge state with this run's snapshot.")
                    .clicked()
                {
                    viz.import_requested_path = Some(dir.join("snapshot.json"));
                }
            });
        }
    }
    if let Some(result) = &viz.last_import {
        match result {
            Ok(path) => {
                ui.label(
                    egui::RichText::new(format!("Loaded {}", path.display()))
                        .small()
                        .color(egui::Color32::LIGHT_GREEN),
                );
            }
            Err(err) => {
                ui.label(
                    egui::RichText::new(format!("Load failed: {}", err))
                        .small()
                        .color(egui::Color32::LIGHT_RED),
                );
            }
        }
    }
}

/// Walk `output/` and collect every immediate sub-directory that
/// contains a `snapshot.json` — i.e. the v2 round-trip artefact written
/// by `handle_v2_export`. Sorted by path so the UI ordering is stable
/// across rescans.
fn list_export_dirs(root: &std::path::Path) -> Vec<std::path::PathBuf> {
    if !root.exists() {
        return Vec::new();
    }
    let Ok(entries) = std::fs::read_dir(root) else {
        return Vec::new();
    };
    let mut dirs: Vec<std::path::PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_dir() && p.join("snapshot.json").exists())
        .collect();
    dirs.sort();
    dirs
}
