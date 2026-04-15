//! River overlay: blends flow accumulation as blue pixels onto a terrain image.

use bevy::prelude::*;
use ymir_core::lakes::detection::LakeResult;

use crate::state::FlowCache;
use crate::state::FlowState;

/// Blend river overlay onto an already-rendered terrain image.
/// Lake cells are skipped to avoid drawing rivers inside lakes.
pub fn render_river_overlay_on_image(
    flow_cache: &FlowCache,
    lake_result: Option<&LakeResult>,
    image: &mut Image,
) {
    if !matches!(flow_cache.state, FlowState::Completed { .. }) {
        return;
    }
    let Some(ref result) = flow_cache.result else {
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

            // Skip lake cells
            if let Some(lr) = lake_result {
                if lr.lake_map[y * w + x] > 0 {
                    continue;
                }
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

/// Render lake cells as flat blue-green water on the terrain image.
pub fn render_lake_overlay_on_image(lake_result: &LakeResult, image: &mut Image) {
    let w = lake_result.width;
    let h = lake_result.height;

    if image.width() != w as u32 || image.height() != h as u32 {
        return;
    }

    let data = image.data.as_mut().unwrap();

    for y in 0..h {
        for x in 0..w {
            if lake_result.lake_map[y * w + x] == 0 {
                continue;
            }
            // Y-flip for Bevy
            let idx = ((h - 1 - y) * w + x) * 4;
            data[idx] = 40;
            data[idx + 1] = 90;
            data[idx + 2] = 140;
        }
    }
}

fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    (a as f32 + (b as f32 - a as f32) * t).clamp(0.0, 255.0) as u8
}
