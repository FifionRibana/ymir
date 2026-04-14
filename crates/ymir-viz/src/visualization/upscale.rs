//! Upscale visualization: renders the FBM-upscaled heightmap to the terrain
//! texture when the UpscaleFbm phase is active and a result is available.

use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

use super::colormap::hypsometric_colormap;
use super::render::TerrainDisplay;
use crate::state::{FbmState, IsostasyCache, PipelinePhase, UpscaleCache, ViewMode, ViewState};

/// Update the terrain texture with the upscaled heightmap when the
/// UpscaleFbm phase is selected and a result is available.
pub fn update_upscale_texture(
    upscale_cache: Res<UpscaleCache>,
    isostasy_cache: Res<IsostasyCache>,
    terrain_display: Res<TerrainDisplay>,
    view_state: Res<ViewState>,
    mut images: ResMut<Assets<Image>>,
) {
    // Only render when in Altitude mode and UpscaleFbm phase
    if view_state.mode != ViewMode::Altitude {
        return;
    }
    if view_state.selected_phase != PipelinePhase::UpscaleFbm {
        return;
    }
    if !matches!(upscale_cache.state, FbmState::Completed { .. }) {
        return;
    }
    // Only re-render when the cache actually changed
    if !upscale_cache.is_changed() {
        return;
    }

    let Some(ref heightmap) = upscale_cache.heightmap else {
        return;
    };
    let Some(handle) = &terrain_display.texture_handle else {
        return;
    };
    let Some(image) = images.get_mut(handle) else {
        return;
    };

    let n = heightmap.width;
    let sea = isostasy_cache.sea_level_normalized as f64;

    // Resize if needed
    if image.width() != n as u32 || image.height() != n as u32 {
        *image = Image::new(
            Extent3d { width: n as u32, height: n as u32, depth_or_array_layers: 1 },
            TextureDimension::D2,
            vec![0u8; n * n * 4],
            TextureFormat::Rgba8UnormSrgb,
            bevy::asset::RenderAssetUsages::MAIN_WORLD
                | bevy::asset::RenderAssetUsages::RENDER_WORLD,
        );
        image.sampler = bevy::image::ImageSampler::nearest();
    }

    for y in 0..n {
        for x in 0..n {
            let h = heightmap.data[y * n + x] as f64;
            let [r, g, b, a] = if h < sea {
                let depth = 1.0 - (h / sea).max(0.0);
                let r = (20.0 + (1.0 - depth) * 30.0) as u8;
                let g = (50.0 + (1.0 - depth) * 80.0) as u8;
                let b = (120.0 + (1.0 - depth) * 50.0) as u8;
                [r, g, b, 255u8]
            } else {
                let t = ((h - sea) / (1.0 - sea)).clamp(0.0, 1.0);
                hypsometric_colormap(0.4 + t * 0.6)
            };
            let _ = image.set_color_at(x as u32, (n - 1 - y) as u32, Color::srgba_u8(r, g, b, a));
        }
    }
}
