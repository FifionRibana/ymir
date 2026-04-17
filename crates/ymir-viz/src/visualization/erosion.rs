//! Erosion visualization: renders the eroded heightmap to the terrain texture
//! when the Erosion phase is active and a result is available.

use bevy::prelude::*;

use super::colormap::{hypsometric_colormap, slope_color};
use super::render::{SolverTerrainSprite, TerrainDisplay, resize_terrain_image};
use crate::state::{
    ErosionCache, ErosionState, FlowCache, IsostasyCache, LakeCache, UpscaleCache, ViewMode,
    ViewState,
};

/// Update the terrain texture with the eroded heightmap.
pub fn update_erosion_texture(
    erosion_cache: Res<ErosionCache>,
    upscale_cache: Res<UpscaleCache>,
    isostasy_cache: Res<IsostasyCache>,
    flow_cache: Res<FlowCache>,
    lake_cache: Res<LakeCache>,
    view_state: Res<ViewState>,
    view_mode: Res<State<ViewMode>>,
    terrain_display: Res<TerrainDisplay>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<SolverTerrainSprite>>,
) {
    if matches!(erosion_cache.state, ErosionState::Idle) {
        return;
    }
    if !erosion_cache.is_changed()
        && !flow_cache.is_changed()
        && !lake_cache.is_changed()
        && !view_state.is_changed()
        && !view_mode.is_changed()
    {
        return;
    }

    // Use erosion heightmap if available, otherwise fall back to upscale
    let heightmap = match erosion_cache.heightmap.as_ref() {
        Some(hm) => hm,
        None => match upscale_cache.heightmap.as_ref() {
            Some(hm) => hm,
            None => return,
        },
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

    // Lake overlay
    if view_state.overlays.lakes {
        if let Some(ref lr) = lake_cache.result {
            super::rivers::render_lake_overlay_on_image(lr, image);
        }
    }

    // River overlay (skips lake cells)
    if view_state.overlays.rivers {
        let lr = lake_cache.result.as_ref();
        super::rivers::render_river_overlay_on_image(&flow_cache, lr, image);
    }
}
