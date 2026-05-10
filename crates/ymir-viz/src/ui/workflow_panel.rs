//! Step 12 Phase 7b — workflow panel for the v2 bridge.
//!
//! Companion to [`crate::ui::parameter_panel_v2`] surfacing the Step 12
//! interleaved tectonic-erosion workflow knobs and run controls. Layout
//!
//! - Top-level toggle `V2WorkflowSpec::{Off, On}` — when `Off`, the
//!   panel collapses to a hint and the rest is hidden so the user can
//!   confirm the bridge is in single-baseline mode at a glance.
//! - Phase A section — `n_cycles`, `k_cycle`, `α`, `β` sliders with
//!   D8-issue defaults; `Run` / `Stop` / `Continue` buttons mapped to
//!   [`crate::bridge::v2::commands::V2Command::RunWorkflowPhaseA`] and
//!   friends.
//! - Phase B section — `hd_grid_size` (256/512/1024/2048),
//!   `num_droplets`, `fbm_amplitude_base`, `grand_scale_tolerance`
//!   sliders; `Run Phase B` button enabled only after Phase A has
//!   completed; `Export HD heightmap` button writes a 16-bit Luma PNG.
//! - Cycle-history table — per-cycle metrics streamed through
//!   [`crate::bridge::v2::events::V2Event::WorkflowCycleCompleted`].
//!   FIFO-capped at [`MAX_HISTORY`] so long runs do not grow the
//!   resource unboundedly.
//!
//! Phases 7b.1 through 7b.5 land the layout (sliders + sections),
//! the click handlers (`Run` / `Stop` / `Continue` / `Run Phase B`),
//! the cycle-history population (consumed from
//! `V2Event::WorkflowCycleCompleted`), the PNG export of the Phase B
//! HD heightmap, and the binary-side mount of the panel into the
//! existing v2 right-side panel respectively.

use std::path::PathBuf;

use bevy::prelude::*;
use bevy_egui::egui;

use crate::bridge::v2::{
    V2PhaseAParams, V2PhaseBParams, V2RunSpec, V2RunState, V2SolverBridge, V2WorkflowSpec,
};
use crate::ui::parameter_panel_v2::V2EditableSpec;

/// Per-cycle metrics history, populated by the Bevy event-poll system
/// from `V2Event::WorkflowCycleCompleted`. FIFO-capped at
/// [`MAX_HISTORY`] entries to bound memory regardless of run length.
#[derive(Resource, Default, Debug)]
pub struct WorkflowCycleHistory {
    pub cycles: Vec<CycleMetricsSnapshot>,
}

/// Snapshot of the per-cycle metrics carried by
/// `V2Event::WorkflowCycleCompleted`. Owned, `Clone`able so the table
/// can be rendered without holding the resource borrow.
#[derive(Clone, Debug)]
pub struct CycleMetricsSnapshot {
    pub cycle_idx: usize,
    pub n_cycles: usize,
    pub erosion_volume_removed: f64,
    pub sea_level_normalized: f64,
    pub mass_drift: f64,
    pub craton_recomputation_change: Option<f64>,
}

/// Maximum number of cycle snapshots retained. A 30-cycle run already
/// exceeds the issue D8 N_cycles cap; bounding the buffer prevents
/// memory growth on user-driven extended runs.
pub const MAX_HISTORY: usize = 30;

/// Last-result carrier for the Phase B HD heightmap PNG export action.
/// Held as a Bevy resource so the panel can render the success / failure
/// status across frames without re-doing the I/O each tick.
#[derive(Resource, Default, Debug)]
pub struct WorkflowExportState {
    pub last_export: Option<Result<PathBuf, String>>,
}

impl WorkflowCycleHistory {
    /// Append a snapshot, evicting from the front when the buffer
    /// overflows [`MAX_HISTORY`].
    pub fn push(&mut self, snap: CycleMetricsSnapshot) {
        self.cycles.push(snap);
        if self.cycles.len() > MAX_HISTORY {
            let drop_count = self.cycles.len() - MAX_HISTORY;
            self.cycles.drain(..drop_count);
        }
    }

    pub fn clear(&mut self) {
        self.cycles.clear();
    }
}

/// Render the workflow panel into the supplied `egui::Ui`.
pub fn draw(
    ui: &mut egui::Ui,
    spec_state: &mut V2EditableSpec,
    bridge: &mut V2SolverBridge,
    history: &mut WorkflowCycleHistory,
    export_state: &mut WorkflowExportState,
) {
    ui.heading("Workflow (Step 12)");
    ui.add_space(4.0);

    // Off/On toggle — drives `spec_state.0.workflow` directly.
    let mut workflow_on = matches!(spec_state.0.workflow, V2WorkflowSpec::On { .. });
    if ui
        .checkbox(&mut workflow_on, "Enable interleaved tectonic-erosion workflow")
        .on_hover_text(
            "Enables Phase A (multi-cycle harness + low-res erosion) and \
             Phase B (HD finalization). Off (default) keeps the bridge in \
             single-baseline mode — no behaviour change vs Step 11.",
        )
        .changed()
    {
        spec_state.0.workflow = if workflow_on {
            V2WorkflowSpec::On {
                phase_a: V2PhaseAParams::default(),
                phase_b: V2PhaseBParams::default(),
            }
        } else {
            V2WorkflowSpec::Off
        };
    }

    if !workflow_on {
        ui.add_space(4.0);
        ui.label(
            egui::RichText::new(
                "Workflow disabled — single-baseline runs only. Toggle on to \
                 access Phase A (multi-cycle interleaved erosion) and Phase B \
                 (HD finalization).",
            )
            .small()
            .weak(),
        );
        return;
    }

    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(
            "Phase A current scope: counter-isostasy + sea_level adaptation. \
             Border curvature emerges from Phase B HD pipeline. Step 12.X \
             follow-up planned for direct Phase A curvature mechanism.",
        )
        .small()
        .italics()
        .weak(),
    );
    ui.add_space(4.0);

    let (status_text, status_color) = workflow_status_label(&bridge.state);
    ui.colored_label(status_color, status_text);

    ui.add_space(8.0);
    ui.separator();

    // Phase A section. The mutable destructure of `spec_state.0.workflow`
    // is scoped to the slider call so the borrow ends before
    // `phase_a_buttons` reads `spec_state.0` immutably.
    egui::CollapsingHeader::new("Phase A — multi-cycle interleaved")
        .default_open(true)
        .show(ui, |ui| {
            if let V2WorkflowSpec::On { phase_a, .. } = &mut spec_state.0.workflow {
                phase_a_sliders(ui, phase_a);
            }
            ui.add_space(4.0);
            phase_a_buttons(ui, &spec_state.0, bridge);
        });

    ui.add_space(4.0);

    egui::CollapsingHeader::new("Phase B — HD finalization")
        .default_open(true)
        .show(ui, |ui| {
            if let V2WorkflowSpec::On { phase_b, .. } = &mut spec_state.0.workflow {
                phase_b_sliders(ui, phase_b);
            }
            ui.add_space(4.0);
            phase_b_buttons(ui, &spec_state.0, bridge, export_state);
        });

    ui.add_space(4.0);

    egui::CollapsingHeader::new("Cycle history (per-cycle metrics)")
        .default_open(false)
        .show(ui, |ui| {
            cycle_history_table(ui, history);
        });
}

/// Status badge string + colour. Reads `V2RunState` and discriminates
/// the workflow-relevant variants. `Running` maps to "Phase A running"
/// only when the spec carries `V2WorkflowSpec::On`; a baseline run
/// (Off) shows the generic "Running" badge instead.
fn workflow_status_label(state: &V2RunState) -> (String, egui::Color32) {
    match state {
        V2RunState::Running { spec, step, total, started_at, .. }
            if matches!(spec.workflow, V2WorkflowSpec::On { .. }) =>
        {
            let elapsed = started_at
                .as_ref()
                .map(|t| t.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            (
                format!("Phase A running — cycle {}/{} · {:.1}s", step, total, elapsed),
                egui::Color32::YELLOW,
            )
        }
        V2RunState::WorkflowPhaseACompleted { cycles_run, elapsed, .. } => (
            format!(
                "Phase A done — {} cycles in {:.1}s",
                cycles_run,
                elapsed.as_secs_f64()
            ),
            egui::Color32::LIGHT_GREEN,
        ),
        V2RunState::WorkflowPhaseBCompleted {
            hd_nx,
            hd_ny,
            grand_scale_deviation_p95,
            elapsed,
            ..
        } => (
            format!(
                "Phase B done — {hd_nx}\u{00d7}{hd_ny} HD in {:.1}s · p95 = {:.4}",
                elapsed.as_secs_f64(),
                grand_scale_deviation_p95
            ),
            egui::Color32::from_rgb(0x80, 0xD0, 0xFF),
        ),
        V2RunState::Failed { error } => {
            (format!("Failed: {}", error), egui::Color32::RED)
        }
        _ => ("Workflow idle — no Phase A run yet".to_string(), egui::Color32::GRAY),
    }
}

fn phase_a_sliders(ui: &mut egui::Ui, phase_a: &mut V2PhaseAParams) {
    ui.add(egui::Slider::new(&mut phase_a.n_cycles, 1..=30).text("N_cycles"))
        .on_hover_text(
            "Number of Phase A cycles. Default 5 (D8 conservative). Aggressive \
             demos try \u{2265} 15 to surface visible counter-isostasy effects \
             — see Phase 8 report.",
        );
    ui.add(
        egui::Slider::new(&mut phase_a.k_cycle, 5..=100).text("k_cycle (steps per cycle)"),
    )
    .on_hover_text("Tectonic harness steps run between successive erosion passes.");
    ui.add(
        egui::Slider::new(&mut phase_a.alpha, 0.001..=0.10)
            .text("\u{03b1} (erosion rate)")
            .step_by(0.001)
            .logarithmic(true),
    )
    .on_hover_text(
        "Erosion intensity per pass. Default 0.01 (D8 conservative). \
         Aggressive: \u{03b1} \u{2265} 0.05 with N_cycles \u{2265} 15 \
         produces visible terrain change.",
    );
    ui.add(
        egui::Slider::new(&mut phase_a.beta, 0.0..=1.0)
            .text("\u{03b2} (downhill redistribution)")
            .step_by(0.05),
    )
    .on_hover_text(
        "0 = pure diffusive (mass loss); 1 = mass-conserving downhill \
         redistribution to NESW priority neighbour.",
    );
}

fn phase_a_buttons(ui: &mut egui::Ui, spec: &V2RunSpec, bridge: &mut V2SolverBridge) {
    let is_running = matches!(bridge.state, V2RunState::Running { .. });
    let phase_a_done = matches!(bridge.state, V2RunState::WorkflowPhaseACompleted { .. });

    let can_run = !is_running;
    let can_stop = is_running;
    let can_continue = phase_a_done;

    ui.horizontal(|ui| {
        if ui
            .add_enabled(can_run, egui::Button::new("\u{25b6} Run Phase A"))
            .on_hover_text(
                "Submit a fresh Phase A multi-cycle run with the current spec. \
                 Streams Progress + WorkflowCycleCompleted events; the dashboard \
                 displays per-cycle metrics on the left panel.",
            )
            .clicked()
        {
            if let Err(e) = bridge.submit_workflow_phase_a(spec.clone()) {
                eprintln!("[ymir-viz] failed to submit workflow Phase A: {}", e);
            }
        }
        if ui
            .add_enabled(can_stop, egui::Button::new("\u{23f9} Stop"))
            .on_hover_text(
                "Signal the running Phase A loop to abort at the next cycle \
                 boundary (latency: one harness step \u{2248} 5\u{2013}25 s on \
                 64\u{00b2} mantle-on).",
            )
            .clicked()
        {
            bridge.request_cancel();
        }
        if ui
            .add_enabled(can_continue, egui::Button::new("\u{21bb} Continue"))
            .on_hover_text(
                "Continue the Phase A loop from the prior run's final state. \
                 Voronoi-relevant fields (seed, num_plates, grid dims) are \
                 locked from source; user spec contributes additional cycles \
                 + tweaked erosion knobs.",
            )
            .clicked()
        {
            // Clone the source spec + final state out of bridge.state in
            // a self-contained match so the immutable borrow is dropped
            // before the `submit_continue_workflow_phase_a` call below.
            let continue_source = match &bridge.state {
                V2RunState::WorkflowPhaseACompleted { spec: src, final_state, .. } => {
                    Some((src.clone(), final_state.as_ref().clone()))
                }
                _ => None,
            };
            if let Some((source_spec, from_state)) = continue_source {
                let mut next_spec = spec.clone();
                // Voronoi-relevant fields are locked from source so the
                // continuation's plate tessellation matches the prior
                // run's plate_id / plate_type rasters carried by
                // from_state. User-tweaked physics knobs are honoured.
                next_spec.seed = source_spec.seed;
                next_spec.grid_nx = source_spec.grid_nx;
                next_spec.grid_ny = source_spec.grid_ny;
                next_spec.num_plates = source_spec.num_plates;
                next_spec.continental_ratio = source_spec.continental_ratio;
                next_spec.init_mode = source_spec.init_mode;
                if let Err(e) =
                    bridge.submit_continue_workflow_phase_a(next_spec, from_state)
                {
                    eprintln!(
                        "[ymir-viz] failed to submit workflow continue Phase A: {}",
                        e
                    );
                }
            }
        }
    });
}

fn phase_b_sliders(ui: &mut egui::Ui, phase_b: &mut V2PhaseBParams) {
    egui::ComboBox::from_label("HD grid size")
        .selected_text(format!("{}\u{00b2}", phase_b.hd_grid_size))
        .show_ui(ui, |ui| {
            for &res in &[256usize, 512, 1024, 2048] {
                ui.selectable_value(&mut phase_b.hd_grid_size, res, format!("{}\u{00b2}", res));
            }
        });
    ui.add(
        egui::Slider::new(&mut phase_b.num_droplets, 100_000..=10_000_000)
            .text("droplets")
            .logarithmic(true),
    )
    .on_hover_text(
        "Number of rain-drop iterations during HD erosion. Default 5\u{00d7}10\u{2076} \
         matches the Phase B reference target for 2048\u{00b2}.",
    );
    ui.add(
        egui::Slider::new(&mut phase_b.fbm_amplitude_base, 0.01..=0.30)
            .text("FBM amplitude base")
            .step_by(0.01),
    )
    .on_hover_text(
        "Base amplitude of the anisotropic FBM noise injected during the \
         upscale. Default 0.08.",
    );
    ui.add(
        egui::Slider::new(&mut phase_b.grand_scale_tolerance, 0.05..=0.30)
            .text("D5 tolerance (p95)")
            .step_by(0.01),
    )
    .on_hover_text(
        "Acceptance threshold on p95(|HD - LR|). Default 0.10. The L_\u{221e} \
         diagnostic is reported alongside but not gated (Phase 5 reformulation).",
    );
}

fn phase_b_buttons(
    ui: &mut egui::Ui,
    spec: &V2RunSpec,
    bridge: &mut V2SolverBridge,
    export_state: &mut WorkflowExportState,
) {
    let phase_a_done = matches!(bridge.state, V2RunState::WorkflowPhaseACompleted { .. });
    let phase_b_done = matches!(bridge.state, V2RunState::WorkflowPhaseBCompleted { .. });
    let is_running = matches!(bridge.state, V2RunState::Running { .. });

    ui.horizontal(|ui| {
        if ui
            .add_enabled(
                phase_a_done && !is_running,
                egui::Button::new("\u{25b6} Run Phase B (HD)"),
            )
            .on_hover_text(
                "Run HD finalization on the Phase A final state. Single-shot, \
                 \u{2248} 30\u{2013}90 s at 2048\u{00b2} \u{00d7} 5\u{00d7}10\u{2076} droplets.",
            )
            .clicked()
        {
            // Same self-contained match pattern as the Phase A continue
            // button: clone out of `bridge.state` first, drop the borrow,
            // then call `submit_workflow_phase_b`.
            let phase_b_source = match &bridge.state {
                V2RunState::WorkflowPhaseACompleted { spec: src, final_state, .. } => {
                    Some((src.clone(), final_state.as_ref().clone()))
                }
                _ => None,
            };
            if let Some((source_spec, from_state)) = phase_b_source {
                let mut next_spec = spec.clone();
                // Lock voronoi-relevant fields from source so the HD
                // upscale / FBM / erosion act on a heightmap whose
                // plate-derived structure matches the Phase A run.
                next_spec.seed = source_spec.seed;
                next_spec.grid_nx = source_spec.grid_nx;
                next_spec.grid_ny = source_spec.grid_ny;
                next_spec.num_plates = source_spec.num_plates;
                if let Err(e) = bridge.submit_workflow_phase_b(next_spec, from_state) {
                    eprintln!("[ymir-viz] failed to submit workflow Phase B: {}", e);
                }
            }
        }
        if ui
            .add_enabled(
                phase_b_done,
                egui::Button::new("\u{1f4be} Export HD heightmap"),
            )
            .on_hover_text(
                "Save the HD heightmap as a 16-bit grayscale PNG under the \
                 OS temp dir. File path is reported in the panel after success.",
            )
            .clicked()
        {
            // Pull a snapshot out of bridge.state in a self-contained
            // match so the read borrow drops before we touch
            // `export_state` (which is a different resource, so this is
            // technically fine even with overlapping borrows; the
            // pattern is consistent with the other buttons for
            // readability).
            let snapshot = match &bridge.state {
                V2RunState::WorkflowPhaseBCompleted {
                    hd_nx,
                    hd_ny,
                    hd_heightmap,
                    ..
                } => Some((*hd_nx, *hd_ny, hd_heightmap.clone())),
                _ => None,
            };
            if let Some((nx, ny, hm)) = snapshot {
                export_state.last_export = Some(export_hd_heightmap_png(nx, ny, &hm));
            }
        }
    });
    if let Some(result) = &export_state.last_export {
        match result {
            Ok(path) => {
                ui.label(
                    egui::RichText::new(format!("Saved \u{2192} {}", path.display()))
                        .small()
                        .color(egui::Color32::LIGHT_GREEN),
                );
            }
            Err(err) => {
                ui.label(
                    egui::RichText::new(format!("Export failed: {}", err))
                        .small()
                        .color(egui::Color32::LIGHT_RED),
                );
            }
        }
    }
}

/// Convert the Phase B HD heightmap (row-major `f32` in `[0, 1]`) to a
/// 16-bit grayscale PNG and write it under the OS temp directory with a
/// timestamped filename. Returns the absolute path on success, the
/// underlying `image` error message on failure.
///
/// 16-bit Luma is the format Living Landz expects for terrain
/// heightmaps; the doubled bit-depth vs 8-bit Luma matters for the
/// gentle gradients run_erosion produces (8-bit shows visible step
/// banding on >2k² grids).
pub fn export_hd_heightmap_png(
    nx: usize,
    ny: usize,
    heightmap: &[f32],
) -> Result<PathBuf, String> {
    if heightmap.len() != nx * ny {
        return Err(format!(
            "heightmap length ({}) does not match nx \u{00d7} ny ({} \u{00d7} {} = {})",
            heightmap.len(),
            nx,
            ny,
            nx * ny
        ));
    }
    let pixels: Vec<u16> = heightmap
        .iter()
        .map(|&v| (v.clamp(0.0, 1.0) * 65535.0).round() as u16)
        .collect();
    let img = image::ImageBuffer::<image::Luma<u16>, Vec<u16>>::from_raw(
        nx as u32, ny as u32, pixels,
    )
    .ok_or_else(|| "ImageBuffer::from_raw rejected the buffer (size mismatch)".to_string())?;
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let path = std::env::temp_dir().join(format!("ymir_v2_phase_b_{}.png", stamp));
    img.save(&path).map_err(|e| e.to_string())?;
    Ok(path)
}

fn cycle_history_table(ui: &mut egui::Ui, history: &mut WorkflowCycleHistory) {
    if history.cycles.is_empty() {
        ui.label(
            egui::RichText::new("(no cycles yet — run Phase A to populate)")
                .small()
                .weak(),
        );
        return;
    }
    ui.horizontal(|ui| {
        ui.label(
            egui::RichText::new(format!(
                "{} cycles retained (FIFO cap {})",
                history.cycles.len(),
                MAX_HISTORY
            ))
            .small()
            .weak(),
        );
        if ui
            .small_button("Clear")
            .on_hover_text(
                "Empty the cycle-metrics history. New cycles from a subsequent \
                 Phase A run will populate it again.",
            )
            .clicked()
        {
            history.clear();
        }
    });
    ui.add_space(2.0);
    egui::Grid::new("workflow_cycle_history")
        .num_columns(5)
        .spacing([8.0, 2.0])
        .striped(true)
        .show(ui, |ui| {
            ui.label(egui::RichText::new("Cycle").strong());
            ui.label(egui::RichText::new("Erosion vol").strong());
            ui.label(egui::RichText::new("Sea level").strong());
            ui.label(egui::RichText::new("Mass drift").strong());
            ui.label(egui::RichText::new("\u{0394} craton").strong());
            ui.end_row();

            for snap in &history.cycles {
                ui.monospace(format!("{} / {}", snap.cycle_idx + 1, snap.n_cycles));
                ui.monospace(format!("{:.3e}", snap.erosion_volume_removed));
                ui.monospace(format!("{:.3}", snap.sea_level_normalized));
                ui.monospace(format!("{:.3e}", snap.mass_drift));
                match snap.craton_recomputation_change {
                    Some(c) => ui.monospace(format!("{:.3}", c)),
                    None => ui.monospace("\u{2014}"),
                };
                ui.end_row();
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// FIFO eviction kicks in once the buffer exceeds `MAX_HISTORY`. We
    /// push `MAX_HISTORY + 5` snapshots and verify that (a) the buffer
    /// length is exactly `MAX_HISTORY`, (b) the oldest 5 cycles were
    /// evicted from the front (so the surviving entries start at
    /// `cycle_idx == 5`).
    #[test]
    fn cycle_history_evicts_fifo_past_max() {
        let mut h = WorkflowCycleHistory::default();
        for i in 0..(MAX_HISTORY + 5) {
            h.push(CycleMetricsSnapshot {
                cycle_idx: i,
                n_cycles: MAX_HISTORY + 5,
                erosion_volume_removed: 0.0,
                sea_level_normalized: 0.0,
                mass_drift: 0.0,
                craton_recomputation_change: None,
            });
        }
        assert_eq!(h.cycles.len(), MAX_HISTORY);
        assert_eq!(h.cycles.first().map(|s| s.cycle_idx), Some(5));
        assert_eq!(h.cycles.last().map(|s| s.cycle_idx), Some(MAX_HISTORY + 4));
    }

    /// Phase 7b.4 — round-trip a tiny synthetic heightmap through the
    /// PNG export path: known input values map to deterministic 16-bit
    /// pixels, the file lands under temp_dir, and re-reading the
    /// pixels recovers the original ramp. Catches future regressions
    /// in the clamp / scaling formula.
    #[test]
    fn export_hd_heightmap_png_roundtrips_pixel_values() {
        // 2x2 heightmap: corner cells span the full [0, 1] range so
        // each clamp branch is exercised. The off-by-one is intentional
        // (`-0.1` and `1.1`) to validate the saturating clamps.
        let nx = 2;
        let ny = 2;
        let hm: Vec<f32> = vec![-0.1_f32, 0.5, 1.1, 0.0];
        let path = super::export_hd_heightmap_png(nx, ny, &hm)
            .expect("export must succeed for a 2x2 buffer");
        assert!(path.exists(), "exported PNG should exist on disk");

        // Re-read the saved PNG and verify the pixel values match the
        // documented mapping (clamp + 65535 scale).
        let img = image::open(&path).expect("re-open exported PNG").to_luma16();
        // -0.1 -> clamp 0 -> 0
        assert_eq!(img.get_pixel(0, 0).0[0], 0);
        // 0.5 -> 32768 (rounded from 32767.5)
        assert_eq!(img.get_pixel(1, 0).0[0], 32768);
        // 1.1 -> clamp 1 -> 65535
        assert_eq!(img.get_pixel(0, 1).0[0], 65535);
        // 0.0 -> 0
        assert_eq!(img.get_pixel(1, 1).0[0], 0);

        // Cleanup: remove the temp file so repeated runs don't pile up.
        let _ = std::fs::remove_file(&path);
    }

    /// Length-mismatch returns an `Err` rather than panicking. Defends
    /// against a future caller hand-mangling the `nx` / `ny` /
    /// `heightmap` triple (the trio is also locked together by the
    /// `WorkflowPhaseBCompleted` event payload, but the export
    /// function is publicly reachable from integration tests).
    #[test]
    fn export_hd_heightmap_png_rejects_length_mismatch() {
        let result = super::export_hd_heightmap_png(2, 2, &[0.0, 0.5, 1.0]);
        assert!(result.is_err(), "expected Err, got {:?}", result);
    }

    #[test]
    fn cycle_history_clear_resets() {
        let mut h = WorkflowCycleHistory::default();
        h.push(CycleMetricsSnapshot {
            cycle_idx: 0,
            n_cycles: 1,
            erosion_volume_removed: 1.0,
            sea_level_normalized: 0.5,
            mass_drift: 0.0,
            craton_recomputation_change: Some(0.1),
        });
        assert_eq!(h.cycles.len(), 1);
        h.clear();
        assert!(h.cycles.is_empty());
    }
}
