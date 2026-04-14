//! Terrain sprite rendering: updates a Bevy Image from the solver's S field.

use bevy::image::ImageSampler;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

use ymir_core::tectonics::solver::field::Field2D;

use super::colormap::hypsometric_colormap;

/// Marker component for the solver terrain sprite.
#[derive(Component)]
pub struct SolverTerrainSprite;

/// Resource holding the current S field and display state.
#[derive(Resource)]
pub struct TerrainDisplay {
    pub s_field: Option<Field2D>,
    pub dirty: bool,
    pub s_range: (f64, f64),
    pub grid_size: usize,
    pub texture_handle: Option<Handle<Image>>,
}

impl Default for TerrainDisplay {
    fn default() -> Self {
        Self {
            s_field: None,
            dirty: false,
            s_range: (0.1, 2.5),
            grid_size: 128,
            texture_handle: None,
        }
    }
}

impl TerrainDisplay {
    pub fn update_field(&mut self, field: Field2D) {
        // Auto-detect range from data
        let mut s_min = f64::INFINITY;
        let mut s_max = f64::NEG_INFINITY;
        for val in field.data() {
            s_min = s_min.min(*val);
            s_max = s_max.max(*val);
        }
        if (s_max - s_min).abs() < 1e-10 {
            s_min -= 0.5;
            s_max += 0.5;
        }
        self.s_range = (s_min, s_max);
        self.grid_size = field.n();
        self.s_field = Some(field);
        self.dirty = true;
    }
}

pub fn setup_solver_terrain_sprite(
    mut commands: Commands,
    mut images: ResMut<Assets<Image>>,
    mut terrain_display: ResMut<TerrainDisplay>,
) {
    let n = terrain_display.grid_size;

    let mut image = Image::new_fill(
        Extent3d { width: n as u32, height: n as u32, depth_or_array_layers: 1 },
        TextureDimension::D2,
        &[40, 120, 160, 255],
        TextureFormat::Rgba8UnormSrgb,
        bevy::asset::RenderAssetUsages::MAIN_WORLD | bevy::asset::RenderAssetUsages::RENDER_WORLD,
    );
    image.sampler = ImageSampler::nearest();

    let handle = images.add(image);
    terrain_display.texture_handle = Some(handle.clone());

    commands.spawn((
        Sprite { image: handle, custom_size: Some(Vec2::new(600.0, 600.0)), ..default() },
        SolverTerrainSprite,
    ));
}

pub fn update_terrain_texture(
    mut terrain_display: ResMut<TerrainDisplay>,
    mut images: ResMut<Assets<Image>>,
) {
    if !terrain_display.dirty {
        return;
    }
    let Some(field) = &terrain_display.s_field else {
        return;
    };
    let Some(handle) = &terrain_display.texture_handle else {
        return;
    };
    let Some(image) = images.get_mut(handle) else {
        return;
    };

    let n = field.n();
    let (s_min, s_max) = terrain_display.s_range;
    let range = (s_max - s_min).max(1e-10);

    // Resize the image if grid size changed
    if image.width() != n as u32 || image.height() != n as u32 {
        *image = Image::new(
            Extent3d { width: n as u32, height: n as u32, depth_or_array_layers: 1 },
            TextureDimension::D2,
            vec![0u8; n * n * 4], // RGBA8 = 4 bytes per pixel
            TextureFormat::Rgba8UnormSrgb,
            bevy::asset::RenderAssetUsages::MAIN_WORLD
                | bevy::asset::RenderAssetUsages::RENDER_WORLD,
        );
        image.sampler = ImageSampler::nearest();
    }

    for y in 0..n {
        for x in 0..n {
            let s = field.get(x, y);
            let t = ((s - s_min) / range).clamp(0.0, 1.0);
            let [r, g, b, a] = hypsometric_colormap(t);
            // Y-flip: image row 0 is top, grid row 0 is bottom
            let _ = image.set_color_at(x as u32, (n - 1 - y) as u32, Color::srgba_u8(r, g, b, a));
        }
    }

    terrain_display.dirty = false;
}
