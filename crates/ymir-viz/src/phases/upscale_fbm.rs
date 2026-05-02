//! UpscaleFbm phase — bicubic interpolation + anisotropic FBM noise
//! on the isostasy output.
//!
//! Reads the [`super::isostasy::IsostasyCache`] heightmap, calls
//! `ymir_core::terrain::upscale::upscale_with_fbm`, stores the
//! resulting high-resolution `GridF32` (typically 512–2048 px) plus
//! the slope field for downstream consumers (Erosion). The render
//! system paints the upscaled heightmap with the hypsometric
//! colormap when `ActivePhase == UpscaleFbm`.

use bevy::prelude::*;
use bevy_egui::egui;
use ymir_core::seed::WorldSeed;
use ymir_core::terrain::upscale::{upscale_with_fbm, FbmUpscaleConfig, UpscaleResult};

use super::erosion::ErosionCache;
use super::isostasy::IsostasyCache;
use crate::bridge::v2::V2SolverBridge;
use crate::pipeline::{ActivePhase, PipelinePhase};
use crate::visualization::v2_viz::{V2VizSprite, V2VizState};

#[derive(Resource, Clone, Debug)]
pub struct FbmParams {
    pub target_size: usize,
    pub octaves: usize,
    pub amplitude_base: f64,
    pub amplitude_slope_factor: f64,
    pub max_anisotropy: f64,
    pub submarine_damping: f64,
    pub domain_warp_strength: f64,
    pub seed: u64,
    pub recompute_requested: bool,
}

impl Default for FbmParams {
    fn default() -> Self {
        let cfg = FbmUpscaleConfig::default();
        Self {
            target_size: cfg.target_size,
            octaves: cfg.octaves,
            amplitude_base: cfg.amplitude_base,
            amplitude_slope_factor: cfg.amplitude_slope_factor,
            max_anisotropy: cfg.max_anisotropy,
            submarine_damping: cfg.submarine_damping,
            domain_warp_strength: cfg.domain_warp_strength,
            seed: 42,
            recompute_requested: false,
        }
    }
}

#[derive(Resource, Default)]
pub struct FbmCache {
    pub result: Option<UpscaleResult>,
    pub last_status: Option<Result<String, String>>,
    last_signature: Option<u64>,
}

impl FbmCache {
    pub fn mark_dirty(&mut self) {
        self.last_signature = None;
    }
}

pub fn draw_section(
    ui: &mut egui::Ui,
    params: &mut FbmParams,
    cache: &FbmCache,
    can_run: bool,
) {
    ui.add_space(4.0);
    ui.heading("Upscale + FBM");
    ui.label(
        egui::RichText::new(
            "Bicubic interpolation + anisotropic FBM noise on the \
             isostasy heightmap.",
        )
        .small()
        .weak(),
    );
    ui.add_space(4.0);

    let presets: &[(&str, usize)] = &[
        ("256²", 256),
        ("512²", 512),
        ("1024²", 1024),
        ("2048²", 2048),
    ];
    egui::ComboBox::from_label("target resolution")
        .selected_text(format!("{}²", params.target_size))
        .show_ui(ui, |ui| {
            for &(label, size) in presets {
                ui.selectable_value(&mut params.target_size, size, label);
            }
        });
    ui.add(
        egui::Slider::new(&mut params.octaves, 1..=10)
            .text("octaves"),
    );
    ui.add(
        egui::Slider::new(&mut params.amplitude_base, 0.0..=0.3)
            .text("amplitude base")
            .step_by(0.005),
    );
    ui.add(
        egui::Slider::new(&mut params.amplitude_slope_factor, 0.0..=8.0)
            .text("amplitude × slope factor")
            .step_by(0.1),
    );
    ui.add(
        egui::Slider::new(&mut params.max_anisotropy, 1.0..=5.0)
            .text("max anisotropy")
            .step_by(0.1),
    );
    ui.add(
        egui::Slider::new(&mut params.submarine_damping, 0.0..=1.0)
            .text("submarine damping")
            .step_by(0.05),
    );
    ui.add(
        egui::Slider::new(&mut params.domain_warp_strength, 0.0..=2.0)
            .text("domain warp")
            .step_by(0.05),
    );
    ui.add(egui::DragValue::new(&mut params.seed).prefix("seed = "));

    ui.add_space(4.0);
    ui.horizontal(|ui| {
        if ui
            .add_enabled(can_run, egui::Button::new("⚙ Run Upscale + FBM"))
            .on_hover_text(
                "Bicubic interpolate the isostasy heightmap to the \
                 target resolution and overlay anisotropic FBM noise. \
                 Requires a populated Isostasy cache.",
            )
            .clicked()
        {
            params.recompute_requested = true;
        }
    });

    if let Some(result) = &cache.result {
        ui.label(
            egui::RichText::new(format!(
                "✓ {}×{} upscaled · slope field cached",
                result.heightmap.width, result.heightmap.height
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
                "No upscale data yet — compute Isostasy first, then \
                 click ⚙ Run Upscale + FBM.",
            )
            .small()
            .weak(),
        );
    }
}

pub fn handle_fbm_compute(
    iso_cache: Res<IsostasyCache>,
    mut params: ResMut<FbmParams>,
    mut cache: ResMut<FbmCache>,
) {
    if !params.recompute_requested {
        return;
    }
    params.recompute_requested = false;

    let Some(iso) = iso_cache.result.as_ref() else {
        cache.last_status = Some(Err(
            "No isostasy data — run the Isostasy phase first.".to_string(),
        ));
        return;
    };

    let cfg = FbmUpscaleConfig {
        target_size: params.target_size,
        octaves: params.octaves,
        amplitude_base: params.amplitude_base,
        amplitude_slope_factor: params.amplitude_slope_factor,
        max_anisotropy: params.max_anisotropy,
        submarine_damping: params.submarine_damping,
        domain_warp_strength: params.domain_warp_strength,
        ..FbmUpscaleConfig::default()
    };
    let world_seed = WorldSeed::new(params.seed);
    let result = upscale_with_fbm(&iso.heightmap, iso.sea_level_normalized, &world_seed, &cfg);
    cache.last_status = Some(Ok(format!(
        "Upscaled to {}×{}",
        result.heightmap.width, result.heightmap.height
    )));
    cache.result = Some(result);
    cache.last_signature = None;
}

#[allow(clippy::too_many_arguments)]
pub fn render_upscale_phase(
    active: Res<ActivePhase>,
    mut cache: ResMut<FbmCache>,
    iso_cache: Res<IsostasyCache>,
    erosion_cache: Res<ErosionCache>,
    bridge: Res<V2SolverBridge>,
    viz: Res<V2VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<V2VizSprite>>,
) {
    if active.0 != PipelinePhase::UpscaleFbm {
        return;
    }
    let Some((grid, source_tag)) = super::select_grid_for_phase(
        PipelinePhase::UpscaleFbm,
        &bridge,
        &iso_cache,
        &cache,
        &erosion_cache,
    ) else {
        return;
    };
    let nx = grid.width;
    let ny = grid.height;
    let sea_level = iso_cache
        .result
        .as_ref()
        .map(|r| r.sea_level_normalized as f64)
        .unwrap_or(0.4);
    let signature = signature_hash(nx, ny, source_tag, sea_level);
    if cache.last_signature == Some(signature) {
        return;
    }
    super::paint_grid_to_v2_sprite(&grid, &viz, &mut images, &mut sprite_q);
    cache.last_signature = Some(signature);
}

fn signature_hash(nx: usize, ny: usize, source_tag: u8, sea: f64) -> u64 {
    use std::hash::Hasher;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    hasher.write_usize(nx);
    hasher.write_usize(ny);
    hasher.write_u8(source_tag);
    hasher.write(&sea.to_le_bytes());
    hasher.finish()
}
