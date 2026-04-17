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

use crate::state::{DynamicPlateIds, TectonicState, ViewMode, ViewState};
use crate::visualization::render::{SPRITE_BASE_SIZE, sprite_size_for};

pub struct TectonicViewPlugin;

impl Plugin for TectonicViewPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_tectonic).add_systems(
            Update,
            (
                rebuild_tectonic_texture,
                toggle_tectonic_visibility,
                draw_velocity_arrows,
                draw_plate_boundaries,
                draw_boundary_types,
            ),
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

    let nx = init.grid_width as u32;
    let ny = init.grid_height as u32;
    let image_handle = alloc_image(&mut images, nx, ny);

    commands.spawn((
        Sprite {
            image: image_handle.clone(),
            custom_size: Some(sprite_size_for(init.grid_width, init.grid_height)),
            ..default()
        },
        TectonicSprite,
        Visibility::Hidden,
    ));

    commands.insert_resource(TectonicImageHandle(image_handle));
    commands.insert_resource(TectonicState {
        init,
        config,
        seed: seed_val,
        dirty: true,
        generation: 0,
    });
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

    let nx = tectonic.init.grid_width;
    let ny = tectonic.init.grid_height;
    let tex_w = nx as u32;
    let tex_h = ny as u32;

    // Recreate the GPU image if grid dimensions changed
    let needs_resize = images
        .get(&handle_res.0)
        .map(|img| img.width() != tex_w || img.height() != tex_h)
        .unwrap_or(true);

    if needs_resize {
        let new_handle = alloc_image(&mut images, tex_w, tex_h);
        let new_sprite_size = sprite_size_for(nx, ny);
        for mut sprite in sprite_q.iter_mut() {
            sprite.image = new_handle.clone();
            sprite.custom_size = Some(new_sprite_size);
        }
        handle_res.0 = new_handle;
    }

    let Some(image) = images.get_mut(&handle_res.0) else { return };
    let data = image.data.as_mut().unwrap();

    let init = &tectonic.init;

    for y in 0..ny {
        for x in 0..nx {
            // Y-flip: texture row 0 is top of screen, grid row 0 is bottom.
            let idx = ((ny - 1 - y) * nx + x) * 4;
            let thickness = init.thickness.data[y * nx + x];
            let plate_id = init.plate_ids[y * nx + x];

            // Mark cell as boundary if any of its toroidal neighbors differ
            let right = init.plate_ids[y * nx + (x + 1) % nx];
            let down = init.plate_ids[((y + 1) % ny) * nx + x];
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
    view_mode: Res<State<ViewMode>>,
    mut tectonic_q: Query<&mut Visibility, With<TectonicSprite>>,
    mut terrain_q: Query<
        &mut Visibility,
        (With<crate::terrain_view::TerrainSprite>, Without<TectonicSprite>),
    >,
) {
    if !view_mode.is_changed() {
        return;
    }
    let show_tectonic = *view_mode.get() == ViewMode::Tectonics;
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
    view_mode: Res<State<ViewMode>>,
) {
    if *view_mode.get() != ViewMode::Tectonics {
        return;
    }
    let Some(tectonic) = tectonic else { return };
    let init = &tectonic.init;
    let nx = init.grid_width;
    let ny = init.grid_height;
    let sprite = sprite_size_for(nx, ny);
    // With aspect preserved, cell size is isotropic on screen.
    let cell_size = SPRITE_BASE_SIZE / nx.max(ny) as f32;
    let half_x = sprite.x / 2.0;
    let half_y = sprite.y / 2.0;

    for plate in &init.plates {
        let wx = plate.seed_x * cell_size - half_x;
        let wy = plate.seed_y * cell_size - half_y;
        let start = Vec2::new(wx, wy);

        let scale = cell_size * 8.0;
        let end = start + Vec2::new(plate.velocity.0 * scale, plate.velocity.1 * scale);

        let color = match plate.plate_type {
            PlateType::Continental => Color::WHITE,
            PlateType::Oceanic => Color::srgb(0.0, 0.9, 0.9),
        };

        gizmos.arrow_2d(start, end, color);
    }
}

// ── Plate boundary overlay (gizmos) ──────────────────────────────────────

fn draw_plate_boundaries(
    mut gizmos: Gizmos,
    plate_ids_res: Res<DynamicPlateIds>,
    view_state: Res<ViewState>,
    tectonic: Option<Res<TectonicState>>,
) {
    if !view_state.overlays.plates {
        return;
    }

    // Use dynamic plate_ids if available, otherwise fall back to static
    let (ids, nx, ny, plates_opt) = if let Some(ref ids) = plate_ids_res.ids {
        let nx = plate_ids_res.grid_width;
        let ny = plate_ids_res.grid_height;
        if nx == 0 || ny == 0 || ids.len() != nx * ny {
            return;
        }
        (ids.as_slice(), nx, ny, plate_ids_res.plates.as_deref())
    } else if let Some(ref tecto) = tectonic {
        let nx = tecto.init.grid_width;
        let ny = tecto.init.grid_height;
        (tecto.init.plate_ids.as_slice(), nx, ny, Some(tecto.init.plates.as_slice()))
    } else {
        return;
    };

    // Scale: grid coords → world coords. Sprite preserves aspect so cells
    // are isotropic on screen (cell_size same on both axes).
    let sprite = sprite_size_for(nx, ny);
    let cell_size = SPRITE_BASE_SIZE / nx.max(ny) as f32;
    let half_x = sprite.x / 2.0;
    let half_y = sprite.y / 2.0;
    let boundary_color = Color::srgba(1.0, 1.0, 1.0, 0.6);

    // Draw boundary lines between cells of different plates
    for j in 0..ny {
        for i in 0..nx {
            let my_id = ids[j * nx + i];

            // Check right neighbor
            let ni = (i + 1) % nx;
            if ids[j * nx + ni] != my_id {
                let x = (i as f32 + 1.0) * cell_size - half_x;
                let y1 = j as f32 * cell_size - half_y;
                let y2 = (j as f32 + 1.0) * cell_size - half_y;
                gizmos.line_2d(Vec2::new(x, y1), Vec2::new(x, y2), boundary_color);
            }

            // Check bottom neighbor
            let nj = (j + 1) % ny;
            if ids[nj * nx + i] != my_id {
                let y = (j as f32 + 1.0) * cell_size - half_y;
                let x1 = i as f32 * cell_size - half_x;
                let x2 = (i as f32 + 1.0) * cell_size - half_x;
                gizmos.line_2d(Vec2::new(x1, y), Vec2::new(x2, y), boundary_color);
            }
        }
    }

    // Draw seed markers
    let plates = plates_opt.or_else(|| tectonic.as_ref().map(|t| t.init.plates.as_slice()));
    if let Some(plates) = plates {
        for plate in plates {
            if !plate.active {
                continue;
            }

            let wx = plate.seed_x * cell_size - half_x;
            let wy = plate.seed_y * cell_size - half_y;
            let pos = Vec2::new(wx, wy);

            let color = match plate.plate_type {
                PlateType::Continental => Color::srgb(0.72, 0.45, 0.20),
                PlateType::Oceanic => Color::srgb(0.35, 0.55, 0.75),
            };

            gizmos.circle_2d(Isometry2d::from_translation(pos), cell_size * 1.5, color);

            // Velocity arrow scaled to world units
            let arrow_scale = cell_size * 5.0;
            let arrow_end =
                pos + Vec2::new(plate.velocity.0 * arrow_scale, plate.velocity.1 * arrow_scale);
            gizmos.arrow_2d(pos, arrow_end, color);
        }
    }
}

fn draw_boundary_types(
    mut gizmos: Gizmos,
    plate_ids_res: Res<DynamicPlateIds>,
    view_state: Res<ViewState>,
) {
    if !view_state.overlays.boundary_types {
        return;
    }

    let Some(ref bt) = plate_ids_res.boundary_types else {
        return;
    };
    let nx = plate_ids_res.grid_width;
    let ny = plate_ids_res.grid_height;
    if nx == 0 || ny == 0 || bt.len() != nx * ny {
        return;
    }

    use ymir_core::tectonics::boundaries::BoundaryType;

    let sprite = sprite_size_for(nx, ny);
    let cell_size = SPRITE_BASE_SIZE / nx.max(ny) as f32;
    let half_x = sprite.x / 2.0;
    let half_y = sprite.y / 2.0;
    let half_cell = cell_size * 0.5;

    for j in 0..ny {
        for i in 0..nx {
            let btype = bt[j * nx + i];
            let color = match btype {
                BoundaryType::None => continue,
                BoundaryType::Subduction => Color::srgba(1.0, 0.15, 0.15, 0.7),
                BoundaryType::OceanicSubduction => Color::srgba(0.0, 0.9, 0.9, 0.7),
                BoundaryType::ContinentalCollision => Color::srgba(0.9, 0.2, 0.7, 0.7),
                BoundaryType::Rift => Color::srgba(0.2, 0.5, 1.0, 0.7),
            };

            let x = i as f32 * cell_size - half_x + half_cell;
            let y = j as f32 * cell_size - half_y + half_cell;

            gizmos.rect_2d(
                Isometry2d::from_translation(Vec2::new(x, y)),
                Vec2::splat(cell_size * 0.8),
                color,
            );
        }
    }
}
