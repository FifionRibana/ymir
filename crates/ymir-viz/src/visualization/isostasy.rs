//! Isostasy recomputation system: updates the cache and terrain texture
//! when the sea level slider changes or a new solver result arrives.

use bevy::prelude::*;

use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};

use super::colormap::hypsometric_colormap;
use super::render::TerrainDisplay;
use crate::state::{IsostasyCache, IsostasyParams, ViewMode, ViewState};

/// Recompute isostasy when the source field or parameters change,
/// and update the terrain texture when in Altitude view mode.
pub fn update_isostasy(
    terrain_display: Res<TerrainDisplay>,
    isostasy_params: Res<IsostasyParams>,
    mut isostasy_cache: ResMut<IsostasyCache>,
    view_state: Res<ViewState>,
    mut images: ResMut<Assets<Image>>,
) {
    let Some(field) = &terrain_display.s_field else {
        return;
    };

    let needs_recompute =
        !isostasy_cache.valid || isostasy_params.is_changed() || terrain_display.is_changed();

    if !needs_recompute {
        return;
    }

    let config = IsostasyConfig {
        sea_level_fraction: isostasy_params.sea_level_fraction,
        max_elevation_m: isostasy_params.max_elevation_m,
        max_depth_m: isostasy_params.max_depth_m,
        altitude_smoothing_sigma: isostasy_params.altitude_smoothing_sigma,
        ..Default::default()
    };

    let result = compute_isostasy(field, &config);

    isostasy_cache.land_ratio = result.land_ratio;
    isostasy_cache.peak_altitude_m = result.peak_altitude_m;
    isostasy_cache.max_depth_m = result.max_depth_m;
    isostasy_cache.sea_level_normalized = result.sea_level_normalized;
    isostasy_cache.computed_sea_level = isostasy_params.sea_level_fraction;
    isostasy_cache.valid = true;

    // Update the texture when in Altitude view mode
    if view_state.mode != ViewMode::Altitude {
        return;
    }
    let Some(handle) = &terrain_display.texture_handle else {
        return;
    };
    let Some(image) = images.get_mut(handle) else {
        return;
    };

    let n = result.heightmap.width;
    let sea = result.sea_level_normalized as f64;

    // Resize if needed
    if image.width() != n as u32 || image.height() != n as u32 {
        use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
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
            let h = result.heightmap.data[y * n + x] as f64;
            let [r, g, b, a] = if h < sea {
                // Ocean: blue gradient
                let depth = 1.0 - (h / sea).max(0.0);
                let r = (20.0 + (1.0 - depth) * 30.0) as u8;
                let g = (50.0 + (1.0 - depth) * 80.0) as u8;
                let b = (120.0 + (1.0 - depth) * 50.0) as u8;
                [r, g, b, 255u8]
            } else {
                // Land: hypsometric ramp
                let t = ((h - sea) / (1.0 - sea)).clamp(0.0, 1.0);
                hypsometric_colormap(0.4 + t * 0.6)
            };
            // Y-flip: image row 0 is top, grid row 0 is bottom
            let _ = image.set_color_at(x as u32, (n - 1 - y) as u32, Color::srgba_u8(r, g, b, a));
        }
    }
}
