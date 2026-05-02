//! Hydrology phase — flow accumulation + pit-fill + river extraction
//! + lake detection on the eroded heightmap.
//!
//! Reads [`super::erosion::ErosionCache`], computes:
//!   - `compute_flow` → `FlowResult` (filled, direction,
//!     accumulation, basins).
//!   - `extract_rivers` → `RiverNetwork`.
//!   - `detect_lakes` → `LakeResult`.
//!
//! Render system paints the eroded heightmap with the hypsometric
//! colormap, then overlays the accumulation field as a blue gradient
//! (rivers proportional to log accumulation) and lake cells as solid
//! water blue. The legacy bridge had separate
//! `render_river_overlay_on_image` / `render_lake_overlay_on_image`
//! helpers; the equivalents are inlined here.

use bevy::image::ImageSampler;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use bevy_egui::egui;
use ymir_core::lakes::detection::{detect_lakes, LakeConfig, LakeResult};
use ymir_core::terrain::flow::{
    compute_flow, extract_rivers, FlowConfig, FlowResult, RiverConfig, RiverNetwork,
};

use super::erosion::ErosionCache;
use super::isostasy::IsostasyCache;
use crate::pipeline::{ActivePhase, PipelinePhase};
use crate::visualization::colormap::hypsometric_colormap;
use crate::visualization::v2_viz::{V2VizSprite, V2VizState};

#[derive(Resource, Clone, Debug)]
pub struct HydrologyParams {
    /// Override sea_level used by `FlowConfig`. Defaults to the
    /// Isostasy cache's `sea_level_normalized` if present.
    pub sea_level_override: Option<f32>,
    pub stream_threshold: f32,
    pub river_threshold: f32,
    pub major_river_threshold: f32,
    pub lake_min_depth: f32,
    pub lake_min_area: usize,
    pub recompute_requested: bool,
}

impl Default for HydrologyParams {
    fn default() -> Self {
        let r = RiverConfig::default();
        let l = LakeConfig::default();
        Self {
            sea_level_override: None,
            stream_threshold: r.stream_threshold,
            river_threshold: r.river_threshold,
            major_river_threshold: r.major_river_threshold,
            lake_min_depth: l.min_depth,
            lake_min_area: l.min_area,
            recompute_requested: false,
        }
    }
}

#[derive(Resource, Default)]
pub struct HydrologyCache {
    pub flow: Option<FlowResult>,
    pub rivers: Option<RiverNetwork>,
    pub lakes: Option<LakeResult>,
    pub last_status: Option<Result<String, String>>,
    last_signature: Option<u64>,
}

impl HydrologyCache {
    pub fn mark_dirty(&mut self) {
        self.last_signature = None;
    }
}

pub fn draw_section(
    ui: &mut egui::Ui,
    params: &mut HydrologyParams,
    cache: &HydrologyCache,
    can_run: bool,
) {
    ui.add_space(4.0);
    ui.heading("Hydrology");
    ui.label(
        egui::RichText::new(
            "Flow accumulation + pit-fill + river extraction + lake \
             detection on the eroded heightmap.",
        )
        .small()
        .weak(),
    );
    ui.add_space(4.0);

    let mut sea_override_on = params.sea_level_override.is_some();
    let prev = sea_override_on;
    ui.checkbox(&mut sea_override_on, "Override sea_level");
    if sea_override_on != prev {
        params.sea_level_override = if sea_override_on { Some(0.1) } else { None };
    }
    if let Some(sea) = params.sea_level_override.as_mut() {
        ui.add(
            egui::Slider::new(sea, 0.0..=1.0)
                .text("sea_level (override)")
                .step_by(0.005),
        );
    }

    ui.add(
        egui::Slider::new(&mut params.stream_threshold, 100.0..=5000.0)
            .text("stream threshold")
            .logarithmic(true),
    );
    ui.add(
        egui::Slider::new(&mut params.river_threshold, 500.0..=20_000.0)
            .text("river threshold")
            .logarithmic(true),
    );
    ui.add(
        egui::Slider::new(&mut params.major_river_threshold, 2_000.0..=100_000.0)
            .text("major river threshold")
            .logarithmic(true),
    );
    ui.add(
        egui::Slider::new(&mut params.lake_min_depth, 0.0001..=0.05)
            .text("lake min depth")
            .logarithmic(true),
    );
    ui.add(
        egui::Slider::new(&mut params.lake_min_area, 1..=500)
            .text("lake min area (cells)"),
    );

    ui.add_space(4.0);
    ui.horizontal(|ui| {
        if ui
            .add_enabled(can_run, egui::Button::new("⚙ Run Hydrology"))
            .on_hover_text(
                "Compute flow + rivers + lakes on the eroded \
                 heightmap. Synchronous, ~1–5 s at 1024².",
            )
            .clicked()
        {
            params.recompute_requested = true;
        }
    });

    if let (Some(_), Some(rivers), Some(lakes)) = (&cache.flow, &cache.rivers, &cache.lakes) {
        ui.label(
            egui::RichText::new(format!(
                "✓ {} river segments · {} lakes",
                rivers.segments.len(),
                lakes.lakes.len(),
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
                "No hydrology data — Run Erosion first, then click \
                 ⚙ Run Hydrology.",
            )
            .small()
            .weak(),
        );
    }
}

pub fn handle_hydrology_compute(
    erosion_cache: Res<ErosionCache>,
    iso_cache: Res<IsostasyCache>,
    mut params: ResMut<HydrologyParams>,
    mut cache: ResMut<HydrologyCache>,
) {
    if !params.recompute_requested {
        return;
    }
    params.recompute_requested = false;

    let Some(eroded) = erosion_cache.result.as_ref() else {
        cache.last_status = Some(Err(
            "No eroded heightmap — run Erosion first.".to_string(),
        ));
        return;
    };

    let sea_level = params
        .sea_level_override
        .or_else(|| iso_cache.result.as_ref().map(|r| r.sea_level_normalized))
        .unwrap_or(0.1);

    let flow_cfg = FlowConfig { sea_level };
    let flow = compute_flow(&eroded.heightmap, &flow_cfg);

    let river_cfg = RiverConfig {
        stream_threshold: params.stream_threshold,
        river_threshold: params.river_threshold,
        major_river_threshold: params.major_river_threshold,
    };
    let rivers = extract_rivers(&flow, &river_cfg, eroded.heightmap.width, eroded.heightmap.height);

    let lake_cfg = LakeConfig {
        min_depth: params.lake_min_depth,
        min_area: params.lake_min_area,
    };
    let lakes = detect_lakes(
        &eroded.heightmap,
        &flow.filled,
        &flow.direction,
        &flow.basins,
        &lake_cfg,
    );

    cache.last_status = Some(Ok(format!(
        "{} river segments · {} lakes",
        rivers.segments.len(),
        lakes.lakes.len()
    )));
    cache.flow = Some(flow);
    cache.rivers = Some(rivers);
    cache.lakes = Some(lakes);
    cache.last_signature = None;
}

#[allow(clippy::too_many_arguments)]
pub fn render_hydrology_phase(
    active: Res<ActivePhase>,
    mut cache: ResMut<HydrologyCache>,
    erosion_cache: Res<ErosionCache>,
    iso_cache: Res<super::isostasy::IsostasyCache>,
    fbm_cache: Res<super::upscale_fbm::FbmCache>,
    bridge: Res<crate::bridge::v2::V2SolverBridge>,
    params: Res<HydrologyParams>,
    viz: Res<V2VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<V2VizSprite>>,
) {
    if active.0 != PipelinePhase::Hydrology {
        return;
    }
    let (Some(flow), Some(lakes), Some(eroded)) = (
        cache.flow.as_ref(),
        cache.lakes.as_ref(),
        erosion_cache.result.as_ref(),
    ) else {
        // Hydrology not yet computed (or its inputs missing) — fall
        // back to the latest upstream heightmap so the user still
        // sees terrain in the Hydrology view.
        let Some((grid, source_tag)) = super::select_grid_for_phase(
            PipelinePhase::Hydrology,
            &bridge,
            &iso_cache,
            &fbm_cache,
            &erosion_cache,
        ) else {
            return;
        };
        let nx = grid.width;
        let ny = grid.height;
        // Distinct signature space (high bit = 1) so the fallback
        // path's signature can never collide with the primary
        // overlay path below.
        let signature = signature_hash(nx, ny, params.stream_threshold, params.lake_min_depth)
            ^ (1u64 << 63)
            ^ (source_tag as u64);
        if cache.last_signature == Some(signature) {
            return;
        }
        super::paint_grid_to_v2_sprite(&grid, &viz, &mut images, &mut sprite_q);
        cache.last_signature = Some(signature);
        return;
    };
    let nx = eroded.heightmap.width;
    let ny = eroded.heightmap.height;
    let signature = signature_hash(nx, ny, params.stream_threshold, params.lake_min_depth);
    if cache.last_signature == Some(signature) {
        return;
    }

    let Some(handle) = viz.texture_handle.clone() else {
        return;
    };
    let Some(image) = images.get_mut(&handle) else {
        return;
    };

    if image.width() as usize != nx || image.height() as usize != ny {
        *image = Image::new(
            Extent3d {
                width: nx as u32,
                height: ny as u32,
                depth_or_array_layers: 1,
            },
            TextureDimension::D2,
            vec![0u8; nx * ny * 4],
            TextureFormat::Rgba8UnormSrgb,
            bevy::asset::RenderAssetUsages::MAIN_WORLD
                | bevy::asset::RenderAssetUsages::RENDER_WORLD,
        );
        image.sampler = ImageSampler::nearest();
        let new_size = sprite_size_for(nx, ny);
        for mut s in sprite_q.iter_mut() {
            s.custom_size = Some(new_size);
        }
    }

    // Base layer: eroded hypsometric.
    let mut rgba = vec![0u8; nx * ny * 4];
    for j in 0..ny {
        for i in 0..nx {
            let val = eroded.heightmap.data[j * nx + i] as f64;
            let pixel = hypsometric_colormap(val.clamp(0.0, 1.0));
            let img_row = ny - 1 - j;
            let idx = (img_row * nx + i) * 4;
            rgba[idx..idx + 4].copy_from_slice(&pixel);
        }
    }

    // River overlay — log-scaled flow accumulation. Cells above
    // `stream_threshold` get a blue gradient that scales with
    // log(accumulation / stream_threshold). Skip lake cells (drawn
    // in the lake overlay below).
    let acc = &flow.accumulation.data;
    let lake_map = &lakes.lake_map;
    let stream = params.stream_threshold.max(1.0);
    let major = params.major_river_threshold.max(stream * 2.0);
    let log_max = (major / stream).max(1.0).ln();
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            if lake_map[k] != 0 {
                continue;
            }
            let a = acc[k];
            if a < stream {
                continue;
            }
            let t = ((a / stream).max(1.0).ln() / log_max).clamp(0.0, 1.0);
            // shallow stream → light blue, major river → deep blue
            let r = (160.0 - 110.0 * t).round() as u8;
            let g = (200.0 - 90.0 * t).round() as u8;
            let b = (255.0 - 35.0 * t).round() as u8;
            let img_row = ny - 1 - j;
            let idx = (img_row * nx + i) * 4;
            rgba[idx] = r;
            rgba[idx + 1] = g;
            rgba[idx + 2] = b;
            rgba[idx + 3] = 255;
        }
    }

    // Lake overlay — solid water colour on cells in any detected lake.
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            if lake_map[k] == 0 {
                continue;
            }
            let img_row = ny - 1 - j;
            let idx = (img_row * nx + i) * 4;
            rgba[idx] = 30;
            rgba[idx + 1] = 90;
            rgba[idx + 2] = 175;
            rgba[idx + 3] = 255;
        }
    }

    if let Some(data) = image.data.as_mut() {
        if data.len() == rgba.len() {
            data.copy_from_slice(&rgba);
        } else {
            *data = rgba;
        }
    }
    cache.last_signature = Some(signature);
}

fn signature_hash(nx: usize, ny: usize, stream: f32, lake_min_depth: f32) -> u64 {
    use std::hash::Hasher;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    hasher.write_usize(nx);
    hasher.write_usize(ny);
    hasher.write(&stream.to_le_bytes());
    hasher.write(&lake_min_depth.to_le_bytes());
    hasher.finish()
}

const SPRITE_BASE_SIZE: f32 = 600.0;

fn sprite_size_for(grid_width: usize, grid_height: usize) -> Vec2 {
    if grid_width == 0 || grid_height == 0 {
        return Vec2::splat(SPRITE_BASE_SIZE);
    }
    let longer = grid_width.max(grid_height) as f32;
    Vec2::new(
        SPRITE_BASE_SIZE * grid_width as f32 / longer,
        SPRITE_BASE_SIZE * grid_height as f32 / longer,
    )
}
