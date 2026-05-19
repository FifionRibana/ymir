//! Pipeline phase modules.
//!
//! Each phase owns its own `Params` resource (sliders), `Cache`
//! resource (compute results), `draw_section(ui, ...)` egui helper,
//! and Bevy systems for compute + texture render. The right-side
//! parameter panel calls into the active phase's `draw_section`;
//! the central sprite is shared across all phases (V2VizState's
//! texture handle), with each phase's render system gated on
//! `ActivePhase`.

pub mod biome;
pub mod climate;
pub mod erosion;
pub mod hydrology;
pub mod isostasy;
pub mod upscale_fbm;

use bevy::image::ImageSampler;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use ymir_core::grid::GridF32;

use crate::bridge::v2::{V2RunState, V2SolverBridge};
use crate::pipeline::{ActivePhase, PipelinePhase};
use crate::visualization::colormap::hypsometric_colormap;
use crate::visualization::v2_viz::{V2VizSprite, V2VizState};

/// Step 8.6 follow-up — drop every cached "last rendered" signature
/// when the active phase changes. Without this, switching out of a
/// phase and back leaves the per-phase signature matching the
/// just-computed sprite contents — even though another phase has
/// painted over the texture in the meantime — so the render system
/// thinks "nothing to do" and the user is stuck with the wrong image.
pub fn invalidate_renders_on_phase_change(
    active: Res<ActivePhase>,
    mut viz: ResMut<V2VizState>,
    mut iso_cache: ResMut<isostasy::IsostasyCache>,
    mut fbm_cache: ResMut<upscale_fbm::FbmCache>,
    mut ero_cache: ResMut<erosion::ErosionCache>,
    mut hyd_cache: ResMut<hydrology::HydrologyCache>,
) {
    if !active.is_changed() {
        return;
    }
    viz.last_signature = None;
    iso_cache.mark_dirty();
    fbm_cache.mark_dirty();
    ero_cache.mark_dirty();
    hyd_cache.mark_dirty();
}

/// Step 8.6 follow-up — pick the GridF32 to paint for the given phase
/// view. Tries the phase's own cache first, then walks the upstream
/// chain so the user always sees *something* meaningful when they
/// switch into a phase view that has not been computed yet.
///
/// The companion `u8` is a "source tag" that callers fold into their
/// repaint signature so a fallback → primary upgrade still triggers a
/// re-render.
pub fn select_grid_for_phase(
    phase: PipelinePhase,
    bridge: &V2SolverBridge,
    iso: &isostasy::IsostasyCache,
    fbm: &upscale_fbm::FbmCache,
    ero: &erosion::ErosionCache,
) -> Option<(GridF32, u8)> {
    let s_field_grid = || s_field_to_grid(bridge);
    match phase {
        PipelinePhase::Isostasy => iso
            .result
            .as_ref()
            .map(|r| (r.heightmap.clone(), 0u8))
            .or_else(|| s_field_grid().map(|g| (g, 1u8))),
        PipelinePhase::UpscaleFbm => fbm
            .result
            .as_ref()
            .map(|r| (r.heightmap.clone(), 0u8))
            .or_else(|| iso.result.as_ref().map(|r| (r.heightmap.clone(), 1u8)))
            .or_else(|| s_field_grid().map(|g| (g, 2u8))),
        PipelinePhase::Erosion => ero
            .preview_heightmap
            .as_ref()
            .map(|h| (h.clone(), 0u8))
            .or_else(|| ero.result.as_ref().map(|r| (r.heightmap.clone(), 1u8)))
            .or_else(|| fbm.result.as_ref().map(|r| (r.heightmap.clone(), 2u8)))
            .or_else(|| iso.result.as_ref().map(|r| (r.heightmap.clone(), 3u8)))
            .or_else(|| s_field_grid().map(|g| (g, 4u8))),
        PipelinePhase::Hydrology => ero
            .result
            .as_ref()
            .map(|r| (r.heightmap.clone(), 0u8))
            .or_else(|| fbm.result.as_ref().map(|r| (r.heightmap.clone(), 1u8)))
            .or_else(|| iso.result.as_ref().map(|r| (r.heightmap.clone(), 2u8)))
            .or_else(|| s_field_grid().map(|g| (g, 3u8))),
        _ => None,
    }
}

fn s_field_to_grid(bridge: &V2SolverBridge) -> Option<GridF32> {
    match &bridge.state {
        V2RunState::Completed { final_state, .. }
        | V2RunState::Imported { final_state, .. } => {
            let data: Vec<f32> = final_state.s_field.iter().map(|&v| v as f32).collect();
            Some(GridF32::from_vec(final_state.nx, final_state.ny, data))
        }
        V2RunState::Running { peek_state: Some(p), .. } => {
            let data: Vec<f32> = p.s_field.iter().map(|&v| v as f32).collect();
            Some(GridF32::from_vec(p.nx, p.ny, data))
        }
        _ => None,
    }
}

/// Paint a `GridF32` to the v2 sprite using the hypsometric colormap.
/// Resizes the underlying image + sprite to the grid dimensions if
/// needed. Used by phase render systems for both their primary output
/// and their upstream-fallback path.
pub fn paint_grid_to_v2_sprite(
    grid: &GridF32,
    viz: &V2VizState,
    images: &mut Assets<Image>,
    sprite_q: &mut Query<&mut Sprite, With<V2VizSprite>>,
) {
    let nx = grid.width;
    let ny = grid.height;
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
    let mut rgba = vec![0u8; nx * ny * 4];
    for j in 0..ny {
        for i in 0..nx {
            let val = grid.data[j * nx + i] as f64;
            let pixel = hypsometric_colormap(val.clamp(0.0, 1.0));
            let img_row = ny - 1 - j;
            let idx = (img_row * nx + i) * 4;
            rgba[idx..idx + 4].copy_from_slice(&pixel);
        }
    }
    if let Some(data) = image.data.as_mut() {
        if data.len() == rgba.len() {
            data.copy_from_slice(&rgba);
        } else {
            *data = rgba;
        }
    }
}

const SPRITE_BASE_SIZE: f32 = 600.0;

pub(crate) fn sprite_size_for(grid_width: usize, grid_height: usize) -> Vec2 {
    if grid_width == 0 || grid_height == 0 {
        return Vec2::splat(SPRITE_BASE_SIZE);
    }
    let longer = grid_width.max(grid_height) as f32;
    Vec2::new(
        SPRITE_BASE_SIZE * grid_width as f32 / longer,
        SPRITE_BASE_SIZE * grid_height as f32 / longer,
    )
}
