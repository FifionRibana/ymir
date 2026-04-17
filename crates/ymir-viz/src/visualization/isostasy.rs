//! Isostasy recomputation system: updates the cache and terrain texture
//! when the sea level slider changes or a new solver result arrives.

use bevy::prelude::*;

use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};

use super::colormap::{hypsometric_colormap, slope_color};
use super::render::{SolverTerrainSprite, TerrainDisplay, resize_terrain_image};
use crate::state::{IsostasyCache, IsostasyParams, ViewMode};

/// Recompute the isostasy cache when parameters change.
/// Runs in ALL phases — the cache is needed by FBM, erosion, etc.
pub fn recompute_isostasy_cache(
    terrain_display: Res<TerrainDisplay>,
    isostasy_params: Res<IsostasyParams>,
    mut isostasy_cache: ResMut<IsostasyCache>,
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
    isostasy_cache.heightmap = Some(result.heightmap);
}

/// Render the isostasy heightmap to the terrain texture.
/// Only runs when the active phase is Tectonics or Isostasy
/// (controlled by run_if in the scheduler).
pub fn render_isostasy_texture(
    isostasy_cache: Res<IsostasyCache>,
    terrain_display: Res<TerrainDisplay>,
    view_mode: Res<State<ViewMode>>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<SolverTerrainSprite>>,
) {
    if !isostasy_cache.is_changed() && !view_mode.is_changed() {
        return;
    }
    let Some(ref heightmap) = isostasy_cache.heightmap else {
        return;
    };
    let Some(handle) = &terrain_display.texture_handle else {
        return;
    };
    let Some(image) = images.get_mut(handle) else {
        return;
    };

    let nx = heightmap.width;
    let ny = heightmap.height;
    let sea = isostasy_cache.sea_level_normalized as f64;

    resize_terrain_image(image, nx, ny, &mut sprite_q);

    let is_slope = *view_mode.get() == ViewMode::Slope;

    for y in 0..ny {
        for x in 0..nx {
            let [r, g, b, a] = if is_slope {
                let (gx, gy) = heightmap.gradient_at(x, y);
                slope_color(gx, gy)
            } else {
                let h = heightmap.data[y * nx + x] as f64;
                if h < sea {
                    let depth = 1.0 - (h / sea).max(0.0);
                    let r = (20.0 + (1.0 - depth) * 30.0) as u8;
                    let g = (50.0 + (1.0 - depth) * 80.0) as u8;
                    let b = (120.0 + (1.0 - depth) * 50.0) as u8;
                    [r, g, b, 255u8]
                } else {
                    let t = ((h - sea) / (1.0 - sea)).clamp(0.0, 1.0);
                    hypsometric_colormap(0.4 + t * 0.6)
                }
            };
            let _ =
                image.set_color_at(x as u32, (ny - 1 - y) as u32, Color::srgba_u8(r, g, b, a));
        }
    }
}
