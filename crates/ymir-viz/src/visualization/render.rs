//! Terrain sprite rendering: updates a Bevy Image from the solver's S field.

use bevy::image::ImageSampler;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

use ymir_core::tectonics::solver::field::Field2D;

use super::colormap::hypsometric_colormap;

/// Marker component for the solver terrain sprite.
#[derive(Component)]
pub struct SolverTerrainSprite;

/// Base sprite size in world units. The longer grid axis renders at this
/// size; the shorter axis is scaled down proportionally to preserve the
/// grid's aspect ratio on screen.
pub const SPRITE_BASE_SIZE: f32 = 600.0;

/// Compute the sprite's on-screen size so the longer grid dimension fills
/// [`SPRITE_BASE_SIZE`] and the shorter dimension is proportional. Returns
/// `(SPRITE_BASE_SIZE, SPRITE_BASE_SIZE)` when either dimension is zero, to
/// avoid division-by-zero on uninitialised state.
pub fn sprite_size_for(grid_width: usize, grid_height: usize) -> Vec2 {
    if grid_width == 0 || grid_height == 0 {
        return Vec2::splat(SPRITE_BASE_SIZE);
    }
    let longer = grid_width.max(grid_height) as f32;
    Vec2::new(
        SPRITE_BASE_SIZE * grid_width as f32 / longer,
        SPRITE_BASE_SIZE * grid_height as f32 / longer,
    )
}

/// Reallocate the terrain image to `(grid_width, grid_height)` if its
/// current size doesn't match, and keep the sprite's `custom_size` in sync
/// with the new aspect ratio. Returns `true` if the image was actually
/// reallocated (callers may want to redraw the full texture).
pub fn resize_terrain_image(
    image: &mut Image,
    grid_width: usize,
    grid_height: usize,
    sprite_q: &mut Query<&mut Sprite, With<SolverTerrainSprite>>,
) -> bool {
    if image.width() == grid_width as u32 && image.height() == grid_height as u32 {
        return false;
    }
    *image = Image::new(
        Extent3d {
            width: grid_width as u32,
            height: grid_height as u32,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        vec![0u8; grid_width * grid_height * 4],
        TextureFormat::Rgba8UnormSrgb,
        bevy::asset::RenderAssetUsages::MAIN_WORLD | bevy::asset::RenderAssetUsages::RENDER_WORLD,
    );
    image.sampler = ImageSampler::nearest();

    let new_size = sprite_size_for(grid_width, grid_height);
    for mut sprite in sprite_q.iter_mut() {
        sprite.custom_size = Some(new_size);
    }
    true
}

/// Resource holding the current S field and display state.
#[derive(Resource)]
pub struct TerrainDisplay {
    pub s_field: Option<Field2D>,
    pub dirty: bool,
    pub s_range: (f64, f64),
    pub grid_width: usize,
    pub grid_height: usize,
    pub texture_handle: Option<Handle<Image>>,
}

impl Default for TerrainDisplay {
    fn default() -> Self {
        Self {
            s_field: None,
            dirty: false,
            s_range: (0.1, 2.5),
            grid_width: 128,
            grid_height: 128,
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
        self.grid_width = field.nx();
        self.grid_height = field.ny();
        self.s_field = Some(field);
        self.dirty = true;
    }
}

pub fn setup_solver_terrain_sprite(
    mut commands: Commands,
    mut images: ResMut<Assets<Image>>,
    mut terrain_display: ResMut<TerrainDisplay>,
) {
    let nx = terrain_display.grid_width;
    let ny = terrain_display.grid_height;

    let mut image = Image::new_fill(
        Extent3d { width: nx as u32, height: ny as u32, depth_or_array_layers: 1 },
        TextureDimension::D2,
        &[40, 120, 160, 255],
        TextureFormat::Rgba8UnormSrgb,
        bevy::asset::RenderAssetUsages::MAIN_WORLD | bevy::asset::RenderAssetUsages::RENDER_WORLD,
    );
    image.sampler = ImageSampler::nearest();

    let handle = images.add(image);
    terrain_display.texture_handle = Some(handle.clone());

    commands.spawn((
        Sprite {
            image: handle,
            custom_size: Some(sprite_size_for(nx, ny)),
            ..default()
        },
        SolverTerrainSprite,
    ));
}

pub fn update_terrain_texture(
    mut terrain_display: ResMut<TerrainDisplay>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<SolverTerrainSprite>>,
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

    let nx = field.nx();
    let ny = field.ny();
    let (s_min, s_max) = terrain_display.s_range;
    let range = (s_max - s_min).max(1e-10);

    resize_terrain_image(image, nx, ny, &mut sprite_q);

    for j in 0..ny {
        for i in 0..nx {
            let s = field.get(i, j);
            let t = ((s - s_min) / range).clamp(0.0, 1.0);
            let [r, g, b, a] = hypsometric_colormap(t);
            // Y-flip: image row 0 is top, grid row 0 is bottom.
            let _ =
                image.set_color_at(i as u32, (ny - 1 - j) as u32, Color::srgba_u8(r, g, b, a));
        }
    }

    terrain_display.dirty = false;
}
