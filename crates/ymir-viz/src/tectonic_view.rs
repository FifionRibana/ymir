//! Tectonic plate visualization — texture rendering, plate overlays, velocity arrows.
//!
//! Renders `TectonicState` into a Bevy sprite using a crustal-thickness colormap.
//! Plate boundaries are drawn as white pixels in the texture.
//! Velocity vectors are drawn as gizmo arrows over the grid.

use bevy::image::{ImageSampler, ImageSamplerDescriptor};
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::plates::{PlateConfig, PlateType, generate_plates};

use crate::state::{TectonicState, ViewMode, ViewState};

pub struct TectonicViewPlugin;

impl Plugin for TectonicViewPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_tectonic).add_systems(
            Update,
            (rebuild_tectonic_texture, toggle_tectonic_visibility, draw_velocity_arrows),
        );
    }
}

// ── Components & resources ────────────────────────────────────────────────

#[derive(Component)]
pub struct TectonicSprite;

#[derive(Resource)]
pub struct TectonicImageHandle(pub Handle<Image>);

// ── Thickness colormap ────────────────────────────────────────────────────

/// Map crustal thickness [0.1 … 2.0+] to an RGB color.
///
/// Blue  = thin oceanic crust, green = normal continental, red = active collision.
pub fn thickness_color(thickness: f32) -> [u8; 3] {
    const STOPS: &[(f32, u8, u8, u8)] = &[
        (0.10, 0x0d, 0x2b, 0x42), // deep blue — thin oceanic
        (0.25, 0x0d, 0x2b, 0x42),
        (0.40, 0x1a, 0x5a, 0x7e), // mid blue — ocean/transition
        (0.70, 0x4a, 0x6a, 0x3a), // olive green — thin continental
        (1.00, 0x7a, 0x8a, 0x4a), // yellow-green — normal continental
        (1.30, 0xb8, 0xa0, 0x40), // yellow — thickening (early collision)
        (1.60, 0xc8, 0x70, 0x30), // orange — active collision
        (2.00, 0xb0, 0x40, 0x20), // red — major collision
    ];

    if thickness <= STOPS[0].0 {
        let (_, r, g, b) = STOPS[0];
        return [r, g, b];
    }
    for i in 1..STOPS.len() {
        let (t0, r0, g0, b0) = STOPS[i - 1];
        let (t1, r1, g1, b1) = STOPS[i];
        if thickness <= t1 {
            let t = (thickness - t0) / (t1 - t0);
            return [lerp_u8(r0, r1, t), lerp_u8(g0, g1, t), lerp_u8(b0, b1, t)];
        }
    }
    [0x80, 0x10, 0x10] // dark red — beyond maximum
}

fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    (a as f32 + (b as f32 - a as f32) * t).clamp(0.0, 255.0) as u8
}

// ── Startup system ────────────────────────────────────────────────────────

fn setup_tectonic(mut commands: Commands, mut images: ResMut<Assets<Image>>) {
    let config = PlateConfig::default();
    let seed_val: u64 = 42;
    let init = generate_plates(&config, &WorldSeed::new(seed_val));

    let size = init.grid_size as u32;
    let image_handle = alloc_image(&mut images, size, size);

    commands.spawn((
        Sprite {
            image: image_handle.clone(),
            custom_size: Some(Vec2::new(size as f32, size as f32)),
            ..default()
        },
        TectonicSprite,
        Visibility::Hidden,
    ));

    commands.insert_resource(TectonicImageHandle(image_handle));
    commands.insert_resource(TectonicState { init, config, seed: seed_val, dirty: true, generation: 0 });
}

fn alloc_image(images: &mut Assets<Image>, width: u32, height: u32) -> Handle<Image> {
    let pixel_data = vec![0u8; (width * height * 4) as usize];
    let mut image = Image::new(
        Extent3d { width, height, depth_or_array_layers: 1 },
        TextureDimension::D2,
        pixel_data,
        TextureFormat::Rgba8UnormSrgb,
        bevy::asset::RenderAssetUsages::RENDER_WORLD | bevy::asset::RenderAssetUsages::MAIN_WORLD,
    );
    image.sampler = ImageSampler::Descriptor(ImageSamplerDescriptor::nearest());
    images.add(image)
}

// ── Rebuild texture ───────────────────────────────────────────────────────

fn rebuild_tectonic_texture(
    mut tectonic: ResMut<TectonicState>,
    image_handle: Option<ResMut<TectonicImageHandle>>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<TectonicSprite>>,
) {
    if !tectonic.dirty {
        return;
    }
    let Some(mut handle_res) = image_handle else { return };

    let size = tectonic.init.grid_size as u32;

    // Recreate the GPU image if grid dimensions changed
    let needs_resize = images
        .get(&handle_res.0)
        .map(|img| img.width() != size || img.height() != size)
        .unwrap_or(true);

    if needs_resize {
        let new_handle = alloc_image(&mut images, size, size);
        for mut sprite in sprite_q.iter_mut() {
            sprite.image = new_handle.clone();
            sprite.custom_size = Some(Vec2::new(size as f32, size as f32));
        }
        handle_res.0 = new_handle;
    }

    let Some(image) = images.get_mut(&handle_res.0) else { return };
    let data = image.data.as_mut().unwrap();

    let init = &tectonic.init;
    let s = init.grid_size;

    for y in 0..s {
        for x in 0..s {
            let idx = (y * s + x) * 4;
            let thickness = init.thickness.data[y * s + x];
            let plate_id = init.plate_ids[y * s + x];

            // Mark cell as boundary if any of its toroidal neighbors differ
            let right = init.plate_ids[y * s + (x + 1) % s];
            let down = init.plate_ids[((y + 1) % s) * s + x];
            let is_boundary = plate_id != right || plate_id != down;

            let [r, g, b] = if is_boundary {
                [220u8, 220u8, 220u8] // light grey boundary line
            } else {
                thickness_color(thickness)
            };

            data[idx] = r;
            data[idx + 1] = g;
            data[idx + 2] = b;
            data[idx + 3] = 255;
        }
    }

    tectonic.dirty = false;
}

// ── Visibility toggle ─────────────────────────────────────────────────────

fn toggle_tectonic_visibility(
    view_state: Res<ViewState>,
    mut tectonic_q: Query<&mut Visibility, With<TectonicSprite>>,
    mut terrain_q: Query<
        &mut Visibility,
        (With<crate::terrain_view::TerrainSprite>, Without<TectonicSprite>),
    >,
) {
    if !view_state.is_changed() {
        return;
    }
    let show_tectonic = view_state.mode == ViewMode::Tectonics;
    for mut vis in tectonic_q.iter_mut() {
        *vis = if show_tectonic { Visibility::Visible } else { Visibility::Hidden };
    }
    for mut vis in terrain_q.iter_mut() {
        *vis = if show_tectonic { Visibility::Hidden } else { Visibility::Visible };
    }
}

// ── Velocity arrow gizmos ─────────────────────────────────────────────────

fn draw_velocity_arrows(
    mut gizmos: Gizmos,
    tectonic: Option<Res<TectonicState>>,
    view_state: Res<ViewState>,
) {
    if view_state.mode != ViewMode::Tectonics {
        return;
    }
    let Some(tectonic) = tectonic else { return };
    let init = &tectonic.init;
    let half = init.grid_size as f32 / 2.0;

    for plate in &init.plates {
        // Convert grid coords to world coords (grid centre = world origin, Y flipped)
        let wx = plate.seed_x - half;
        let wy = -(plate.seed_y - half);
        let start = Vec2::new(wx, wy);

        // Scale velocity so arrows span ~8 grid cells
        let scale = 8.0_f32;
        let end = start + Vec2::new(plate.velocity.0 * scale, -plate.velocity.1 * scale);

        let color = match plate.plate_type {
            PlateType::Continental => Color::WHITE,
            PlateType::Oceanic => Color::srgb(0.0, 0.9, 0.9),
        };

        gizmos.arrow_2d(start, end, color);
    }
}
