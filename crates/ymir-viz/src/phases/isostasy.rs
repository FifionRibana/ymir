//! Isostasy phase — Airy buoyancy on the v2 final S̃.
//!
//! Reads the v2 bridge's `Completed` / `Imported` final state, wraps
//! its `s_field` into a `Field2D`, calls
//! `ymir_core::tectonics::isostasy::compute_isostasy`, stores the
//! result in [`IsostasyCache`]. The render system paints the
//! resulting `GridF32` heightmap with the hypsometric colormap onto
//! the shared v2 sprite when `ActivePhase == Isostasy`.

use bevy::prelude::*;
use bevy_egui::egui;
use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig, IsostasyResult};
use ymir_core::tectonics_v2::field::Field2D;

use super::erosion::ErosionCache;
use super::upscale_fbm::FbmCache;
use crate::bridge::v2::{V2RunState, V2SolverBridge};
use crate::pipeline::{ActivePhase, PipelinePhase};
use crate::visualization::v2_viz::{V2VizSprite, V2VizState};

#[derive(Resource, Clone, Debug)]
pub struct IsostasyParams {
    pub rho_crust: f32,
    pub rho_mantle: f32,
    pub max_elevation_m: f32,
    pub max_depth_m: f32,
    pub sea_level_fraction: f32,
    pub altitude_smoothing_sigma: f32,
    /// Set by the panel's "Run Isostasy" button; consumed by
    /// `handle_isostasy_compute` on the next tick.
    pub recompute_requested: bool,
}

impl Default for IsostasyParams {
    fn default() -> Self {
        let cfg = IsostasyConfig::default();
        Self {
            rho_crust: cfg.rho_crust,
            rho_mantle: cfg.rho_mantle,
            max_elevation_m: cfg.max_elevation_m,
            max_depth_m: cfg.max_depth_m,
            sea_level_fraction: cfg.sea_level_fraction,
            altitude_smoothing_sigma: cfg.altitude_smoothing_sigma,
            recompute_requested: false,
        }
    }
}

#[derive(Resource, Default)]
pub struct IsostasyCache {
    pub result: Option<IsostasyResult>,
    pub last_status: Option<Result<String, String>>,
    /// `(nx, ny, params_hash)` of the last successfully-rendered
    /// frame, so the render system avoids redundant repaints.
    last_signature: Option<(usize, usize, u64)>,
}

impl IsostasyCache {
    /// Drop the cached render signature so the next render-system
    /// frame repaints unconditionally. Called after loading state
    /// from disk so the freshly populated `result` actually reaches
    /// the screen even if the dims + params are unchanged from the
    /// last paint.
    pub fn mark_dirty(&mut self) {
        self.last_signature = None;
    }
}

pub fn draw_section(
    ui: &mut egui::Ui,
    params: &mut IsostasyParams,
    cache: &IsostasyCache,
    can_run: bool,
) {
    ui.add_space(4.0);
    ui.heading("Isostasy");
    ui.label(
        egui::RichText::new("Airy buoyancy: tectonic S̃ → altitude (GridF32).")
            .small()
            .weak(),
    );
    ui.add_space(4.0);

    ui.add(
        egui::Slider::new(&mut params.rho_crust, 2500.0..=3000.0)
            .text("ρ_crust (kg/m³)")
            .step_by(10.0),
    );
    ui.add(
        egui::Slider::new(&mut params.rho_mantle, 3100.0..=3500.0)
            .text("ρ_mantle (kg/m³)")
            .step_by(10.0),
    );
    ui.add(
        egui::Slider::new(&mut params.max_elevation_m, 1000.0..=10_000.0)
            .text("max elevation (m)")
            .step_by(100.0),
    );
    ui.add(
        egui::Slider::new(&mut params.max_depth_m, 100.0..=2000.0)
            .text("max depth (m)")
            .step_by(50.0),
    );
    ui.add(
        egui::Slider::new(&mut params.sea_level_fraction, 0.0..=1.0)
            .text("sea-level fraction")
            .step_by(0.01),
    );
    ui.add(
        egui::Slider::new(&mut params.altitude_smoothing_sigma, 0.0..=8.0)
            .text("smoothing σ (cells)")
            .step_by(0.25),
    );

    ui.add_space(4.0);
    ui.horizontal(|ui| {
        if ui
            .add_enabled(can_run, egui::Button::new("⚙ Run Isostasy"))
            .on_hover_text(
                "Compute Airy isostasy on the latest v2 final state. \
                 Requires a Completed or Imported v2 run.",
            )
            .clicked()
        {
            params.recompute_requested = true;
        }
    });

    if let Some(result) = &cache.result {
        ui.label(
            egui::RichText::new(format!(
                "✓ Land ratio {:.0}% · peak {:.0} m · depth {:.0} m",
                result.land_ratio * 100.0,
                result.peak_altitude_m,
                result.max_depth_m,
            ))
            .small()
            .color(egui::Color32::LIGHT_GREEN),
        );
    } else if let Some(Err(err)) = &cache.last_status {
        ui.label(
            egui::RichText::new(format!("Error: {}", err))
                .small()
                .color(egui::Color32::LIGHT_RED),
        );
    } else {
        ui.label(
            egui::RichText::new(
                "No isostasy data yet — Run a v2 tectonics simulation \
                 first, then click ⚙ Run Isostasy.",
            )
            .small()
            .weak(),
        );
    }
}

/// Compute system — fires when the panel's "Run Isostasy" button has
/// been clicked. Reads the bridge's final state, runs the compute,
/// stores the result.
pub fn handle_isostasy_compute(
    bridge: Res<V2SolverBridge>,
    mut params: ResMut<IsostasyParams>,
    mut cache: ResMut<IsostasyCache>,
) {
    if !params.recompute_requested {
        return;
    }
    params.recompute_requested = false;

    let final_ref = match &bridge.state {
        V2RunState::Completed { final_state, .. }
        | V2RunState::Imported { final_state, .. } => final_state.as_ref(),
        _ => {
            cache.last_status = Some(Err(
                "No completed v2 run — start the tectonics solver first.".to_string(),
            ));
            return;
        }
    };

    let s_field = Field2D::from_vec(
        final_ref.nx,
        final_ref.ny,
        final_ref.s_field.clone(),
    );
    let cfg = IsostasyConfig {
        rho_crust: params.rho_crust,
        rho_mantle: params.rho_mantle,
        rho_water: IsostasyConfig::default().rho_water,
        max_elevation_m: params.max_elevation_m,
        max_depth_m: params.max_depth_m,
        sea_level_fraction: params.sea_level_fraction,
        altitude_smoothing_sigma: params.altitude_smoothing_sigma,
    };
    let result = compute_isostasy(&s_field, &cfg);
    cache.last_status = Some(Ok(format!(
        "Computed at {}×{} cells",
        result.heightmap.width, result.heightmap.height
    )));
    cache.result = Some(result);
    // Force a repaint on the next frame.
    cache.last_signature = None;
}

/// Render system — paints the cached isostasy heightmap onto the
/// shared v2 sprite when `ActivePhase == Isostasy`. Falls back to the
/// upstream V2 S field when isostasy hasn't run yet, so the user
/// always sees something on switching into the view.
#[allow(clippy::too_many_arguments)]
pub fn render_isostasy_phase(
    active: Res<ActivePhase>,
    mut cache: ResMut<IsostasyCache>,
    params: Res<IsostasyParams>,
    bridge: Res<V2SolverBridge>,
    fbm_cache: Res<FbmCache>,
    erosion_cache: Res<ErosionCache>,
    viz: Res<V2VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<V2VizSprite>>,
) {
    if active.0 != PipelinePhase::Isostasy {
        return;
    }
    let Some((grid, source_tag)) = super::select_grid_for_phase(
        PipelinePhase::Isostasy,
        &bridge,
        &cache,
        &fbm_cache,
        &erosion_cache,
    ) else {
        return;
    };
    let nx = grid.width;
    let ny = grid.height;
    let signature = (nx, ny, params_hash(&params) ^ (source_tag as u64));
    if cache.last_signature == Some(signature) {
        return;
    }
    super::paint_grid_to_v2_sprite(&grid, &viz, &mut images, &mut sprite_q);
    cache.last_signature = Some(signature);
}

fn params_hash(params: &IsostasyParams) -> u64 {
    use std::hash::Hasher;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    hasher.write(&params.rho_crust.to_le_bytes());
    hasher.write(&params.rho_mantle.to_le_bytes());
    hasher.write(&params.max_elevation_m.to_le_bytes());
    hasher.write(&params.max_depth_m.to_le_bytes());
    hasher.write(&params.sea_level_fraction.to_le_bytes());
    hasher.write(&params.altitude_smoothing_sigma.to_le_bytes());
    hasher.finish()
}
