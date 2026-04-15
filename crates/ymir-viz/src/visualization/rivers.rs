//! River overlay: renders flow accumulation as blue pixels on the terrain texture.

use bevy::prelude::*;

use super::render::TerrainDisplay;
use crate::state::{FlowCache, FlowState, ViewState};

/// Render river overlay on the terrain texture.
/// Runs after terrain/erosion texture rendering, blends blue over existing pixels.
pub fn render_river_overlay(
    flow_cache: Res<FlowCache>,
    view_state: Res<ViewState>,
    terrain_display: Res<TerrainDisplay>,
    mut images: ResMut<Assets<Image>>,
) {
    if !view_state.overlays.rivers {
        return;
    }
    if !matches!(flow_cache.state, FlowState::Completed { .. }) {
        return;
    }
    let Some(ref result) = flow_cache.result else {
        return;
    };

    let Some(handle) = &terrain_display.texture_handle else {
        return;
    };
    let Some(image) = images.get_mut(handle) else {
        return;
    };

    let w = result.accumulation.width;
    let h = result.accumulation.height;

    if image.width() != w as u32 || image.height() != h as u32 {
        return;
    }

    let stream_thr = flow_cache.river_config.stream_threshold;
    let major_thr = flow_cache.river_config.major_river_threshold;
    let log_lo = stream_thr.ln();
    let log_hi = major_thr.ln();
    let log_range = (log_hi - log_lo).max(1e-6);

    let data = image.data.as_mut().unwrap();

    for y in 0..h {
        for x in 0..w {
            let flow = result.accumulation.data[y * w + x];
            if flow < stream_thr {
                continue;
            }

            let t = ((flow.ln() - log_lo) / log_range).clamp(0.0, 1.0);
            let alpha = 0.3 + t * 0.7;

            // Y-flip for Bevy
            let idx = ((h - 1 - y) * w + x) * 4;

            data[idx] = lerp_u8(data[idx], 30, alpha);
            data[idx + 1] = lerp_u8(data[idx + 1], 80, alpha);
            data[idx + 2] = lerp_u8(data[idx + 2], 200, alpha);
        }
    }
}

fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    (a as f32 + (b as f32 - a as f32) * t).clamp(0.0, 255.0) as u8
}
