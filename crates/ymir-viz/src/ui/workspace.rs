//! HD workspace — the functional egui UI (UI rewrite step d1: works
//! end-to-end, not yet polished).
//!
//! Three zones (immediate-mode egui): a left controls panel (seed /
//! resolution / latitude / GENERATE), a central map (the selected HD layer
//! as an egui texture, kept square), and a right inspection panel (the
//! hovered cell's quantities via [`inspect_cell`]).
//!
//! It CONSUMES step b (the worker's `HdEvents` / `HdResult`) and step c
//! (`inspect_cell` / `RiverCellMap`) — it does not re-implement the worker
//! or the accessors. Non-blocking: GENERATE submits to the worker and the
//! UI polls `bridge.hd` each frame (an indeterminate per-phase waiter — the
//! HD phases are opaque blocks; the animated frieze is step d2).
//!
//! The map texture is rebuilt ONLY when the layer or the HD result changes
//! (never per frame). Layers use the canonical data palettes (`colormap`,
//! `Biome::color`, the drainage tiers). Polish — the styled pipeline frieze,
//! Standard/Expert density, fine copper layout, zoom/pan — is step d2.

use std::sync::Arc;

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};

use crate::bridge::c1::{
    inspect_cell, C1RunSpec, C1SolverBridge, CellInspection, HdParams, HdResult, HdState,
    RiverCellMap,
};
use crate::visualization::colormap::hypsometric_bipolar;
use ymir_core::climate::precipitation::precip_mm_per_year;
use ymir_core::tectonics_c1::drainage::{LakeType, Navigability};

pub struct HdWorkspacePlugin;

impl Plugin for HdWorkspacePlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<WorkspaceState>();
        app.add_systems(EguiPrimaryContextPass, draw_workspace);
    }
}

/// The HD layers the map can display, each with its canonical palette.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
enum HdLayer {
    #[default]
    Relief,
    Drainage,
    Precipitation,
    Temperature,
    Biomes,
}

impl HdLayer {
    const ALL: [HdLayer; 5] = [
        HdLayer::Relief,
        HdLayer::Drainage,
        HdLayer::Precipitation,
        HdLayer::Temperature,
        HdLayer::Biomes,
    ];
    fn label(self) -> &'static str {
        match self {
            HdLayer::Relief => "Relief",
            HdLayer::Drainage => "Drainage",
            HdLayer::Precipitation => "Précipitations",
            HdLayer::Temperature => "Température",
            HdLayer::Biomes => "Biomes",
        }
    }
}

#[derive(Resource)]
struct WorkspaceState {
    // Controls.
    seed: u64,
    resolution: usize,
    latitude: f32,
    layer: HdLayer,
    // Current HD product + derived state.
    current: Option<Arc<HdResult>>,
    river_map: Option<RiverCellMap>,
    // Cached map texture + what it was built from (layer + result identity).
    texture: Option<egui::TextureHandle>,
    tex_layer: Option<HdLayer>,
    // Latest inspected cell (hover).
    hover: Option<CellInspection>,
}

impl Default for WorkspaceState {
    fn default() -> Self {
        Self {
            seed: 42,
            resolution: 2048,
            latitude: 45.0,
            layer: HdLayer::Relief,
            current: None,
            river_map: None,
            texture: None,
            tex_layer: None,
            hover: None,
        }
    }
}

fn draw_workspace(
    mut contexts: EguiContexts,
    bridge: Res<C1SolverBridge>,
    mut ws: ResMut<WorkspaceState>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    // 1. Latch a freshly-completed HD result (identity via Arc pointer).
    if let HdState::Completed { result, .. } = &bridge.hd {
        let is_new = match &ws.current {
            Some(cur) => !Arc::ptr_eq(cur, result),
            None => true,
        };
        if is_new {
            ws.current = Some(result.clone());
            ws.river_map = Some(RiverCellMap::from_drainage(&result.drainage));
            ws.texture = None; // force rebuild
            ws.tex_layer = None;
            ws.hover = None;
        }
    }

    // 2. Left panel — controls + progress.
    let hd_running = matches!(bridge.hd, HdState::Running { .. });
    egui::SidePanel::left("controls").exact_width(240.0).show(ctx, |ui| {
        ui.heading("Ymir");
        ui.separator();

        ui.add_enabled_ui(!hd_running, |ui| {
            ui.horizontal(|ui| {
                ui.label("Seed");
                ui.add(egui::DragValue::new(&mut ws.seed));
                if ui.button("🎲").on_hover_text("Seed aléatoire").clicked() {
                    // Deterministic LCG step — no rng dependency needed.
                    ws.seed = ws
                        .seed
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                }
            });
            ui.horizontal(|ui| {
                ui.label("Résolution");
                for &r in &[512usize, 1024, 2048] {
                    ui.selectable_value(&mut ws.resolution, r, format!("{r}"));
                }
            });
            ui.add(egui::Slider::new(&mut ws.latitude, 0.0..=90.0).text("Latitude °"));
        });

        ui.separator();
        ui.horizontal(|ui| {
            if ui.add_enabled(!hd_running, egui::Button::new("⛰ Générer")).clicked() {
                let spec = C1RunSpec { seed: ws.seed, ..C1RunSpec::default() };
                let params = HdParams { target_size: ws.resolution, latitude_deg: ws.latitude };
                let _ = bridge.submit_hd(spec, params);
            }
            if ui.add_enabled(hd_running, egui::Button::new("■ Annuler")).clicked() {
                bridge.request_cancel();
            }
        });

        ui.separator();
        draw_progress(ui, &bridge.hd);
    });

    // 3. Right panel — inspection.
    egui::SidePanel::right("inspection").exact_width(240.0).show(ctx, |ui| {
        ui.heading("Inspection");
        ui.separator();
        match &ws.hover {
            Some(c) => draw_inspection(ui, c),
            None => {
                ui.weak(if ws.current.is_some() {
                    "Survolez la carte."
                } else {
                    "Générez un continent."
                });
            }
        }
    });

    // 4. Central panel — layer selector + the map.
    egui::CentralPanel::default().show(ctx, |ui| {
        ui.horizontal(|ui| {
            for l in HdLayer::ALL {
                if ui.selectable_label(ws.layer == l, l.label()).clicked() {
                    ws.layer = l;
                }
            }
        });
        ui.separator();

        if ws.current.is_none() {
            ui.centered_and_justified(|ui| {
                ui.weak("Aucun continent — cliquez « Générer ».");
            });
            return;
        }

        // (Re)build the texture if the layer or the result changed.
        if ws.texture.is_none() || ws.tex_layer != Some(ws.layer) {
            let layer = ws.layer;
            // Borrow the pieces we need without holding `ws` mutably across
            // the (immutable) reads.
            let img = {
                let hd = ws.current.as_ref().unwrap();
                let rm = ws.river_map.as_ref().unwrap();
                layer_color_image(hd, layer, rm)
            };
            let handle = ctx.load_texture("hd_map", img, egui::TextureOptions::NEAREST);
            ws.texture = Some(handle);
            ws.tex_layer = Some(layer);
        }

        let handle = ws.texture.as_ref().unwrap().clone();
        let avail = ui.available_size();
        let side = avail.x.min(avail.y).max(1.0); // keep it SQUARE
        let resp = ui.add(
            egui::Image::new(egui::load::SizedTexture::new(handle.id(), egui::vec2(side, side)))
                .sense(egui::Sense::hover()),
        );

        // Hover → cell → inspect_cell (step c).
        ws.hover = None;
        if let (Some(pos), Some(hd), Some(rm)) =
            (resp.hover_pos(), ws.current.clone(), ws.river_map.as_ref())
        {
            let rel = pos - resp.rect.min;
            let fx = (rel.x / resp.rect.width()).clamp(0.0, 0.999);
            let fy = (rel.y / resp.rect.height()).clamp(0.0, 0.999);
            let x = (fx * hd.width as f32) as usize;
            let y = (fy * hd.height as f32) as usize;
            ws.hover = Some(inspect_cell(&hd, rm, x, y));
        }
    });
}

/// Per-phase indeterminate progress (the HD phases are opaque blocks).
fn draw_progress(ui: &mut egui::Ui, hd: &HdState) {
    match hd {
        HdState::Idle => {
            ui.weak("Prêt.");
        }
        HdState::Running { current, done, .. } => {
            for r in done {
                ui.weak(format!(
                    "✓ {} — {} ({:.1}s)",
                    r.phase.label(),
                    r.regime.label(),
                    r.elapsed.as_secs_f32()
                ));
            }
            if let Some(p) = current {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("{} …", p.label()));
                });
            }
        }
        HdState::Completed { done, total, .. } => {
            for r in done {
                ui.weak(format!(
                    "✓ {} — {} ({:.1}s)",
                    r.phase.label(),
                    r.regime.label(),
                    r.elapsed.as_secs_f32()
                ));
            }
            ui.colored_label(
                egui::Color32::from_rgb(0x6E, 0xBE, 0x6E),
                format!("✓ Terminé ({:.1}s)", total.as_secs_f32()),
            );
        }
        HdState::Failed { error } => {
            ui.colored_label(egui::Color32::from_rgb(0xE1, 0x78, 0x46), format!("Échec : {error}"));
        }
    }
}

fn draw_inspection(ui: &mut egui::Ui, c: &CellInspection) {
    egui::Grid::new("cell_inspect").num_columns(2).striped(true).show(ui, |ui| {
        ui.label("Cellule");
        ui.label(format!("({}, {})", c.x, c.y));
        ui.end_row();
        ui.label("Altitude");
        ui.label(format!("{:.0} m", c.altitude_m));
        ui.end_row();
        if let Some(d) = c.depth_m {
            ui.label("Océan");
            ui.label(format!("−{d:.0} m"));
            ui.end_row();
        }
        ui.label("Température");
        ui.label(format!("{:.1} °C", c.temperature_c));
        ui.end_row();
        ui.label("Précipitations");
        ui.label(format!("{:.0} mm/an", c.precip_mm));
        ui.end_row();
        ui.label("Ruissellement");
        ui.label(format!("{:.0} mm/an", c.runoff_mm));
        ui.end_row();
        ui.label("Biome");
        ui.label(format!("{:?}", c.biome));
        ui.end_row();
        if let Some(r) = &c.river {
            ui.label("Rivière");
            ui.label(format!("{:?} ({:.0} km²)", r.navigability, r.drainage_km2));
            ui.end_row();
        }
        if let Some(l) = &c.lake {
            ui.label("Lac");
            let kind = match l.lake_type {
                LakeType::Exorheic => "exorhéique",
                LakeType::Endorheic => "endorhéique",
            };
            ui.label(format!("{kind} · {:.0} km²", l.area_km2));
            ui.end_row();
        }
    });
}

// ── Layer → RGBA image (canonical palettes) ──────────────────────────────

fn layer_color_image(hd: &HdResult, layer: HdLayer, river_map: &RiverCellMap) -> egui::ColorImage {
    let (w, h) = (hd.width, hd.height);
    let mut rgba = vec![0u8; w * h * 4];
    for k in 0..w * h {
        let c = match layer {
            HdLayer::Relief => relief_color(hd.eroded.data[k]),
            HdLayer::Precipitation => precip_color(precip_mm_per_year(hd.precipitation.data[k])),
            HdLayer::Temperature => temp_color(hd.temperature.data[k]),
            HdLayer::Biomes => {
                let [r, g, b] = hd.biomes[k].color();
                [r, g, b]
            }
            HdLayer::Drainage => drainage_color(hd, river_map, k),
        };
        rgba[k * 4] = c[0];
        rgba[k * 4 + 1] = c[1];
        rgba[k * 4 + 2] = c[2];
        rgba[k * 4 + 3] = 255;
    }
    egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba)
}

fn relief_color(norm: f32) -> [u8; 3] {
    let [r, g, b, _] = hypsometric_bipolar(norm, 0.5);
    [r, g, b]
}

fn precip_color(mm: f32) -> [u8; 3] {
    if mm < 250.0 {
        [225, 200, 140]
    } else if mm < 500.0 {
        [200, 195, 110]
    } else if mm < 800.0 {
        [150, 180, 90]
    } else if mm < 1500.0 {
        [80, 150, 200]
    } else {
        [30, 90, 200]
    }
}

fn temp_color(t: f32) -> [u8; 3] {
    if t < -5.0 {
        [225, 235, 248]
    } else if t < 5.0 {
        [90, 140, 205]
    } else if t < 20.0 {
        [110, 190, 110]
    } else {
        [225, 120, 70]
    }
}

/// Drainage: dimmed relief base, then lakes, then rivers by navigability.
fn drainage_color(hd: &HdResult, river_map: &RiverCellMap, k: usize) -> [u8; 3] {
    let (w, _h) = (hd.width, hd.height);
    let (x, y) = (k % w, k / w);
    // Base: dimmed hypsometric relief so the water stands out.
    let [br, bg, bb] = relief_color(hd.eroded.data[k]);
    let dim = |v: u8| (v as f32 * 0.55) as u8;
    let mut col = [dim(br), dim(bg), dim(bb)];
    // Lakes.
    if hd.drainage.lake_map[k] != 0 {
        col = [30, 90, 180];
    }
    // Rivers on top, coloured by navigability tier.
    if let Some(info) = river_map.at(x, y) {
        col = match info.navigability {
            Navigability::Ship => [20, 70, 200],
            Navigability::Barge => [40, 110, 230],
            Navigability::SmallBoat => [90, 160, 240],
            Navigability::NonNavigable => [90, 120, 160],
        };
    }
    col
}
