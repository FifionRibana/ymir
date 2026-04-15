use bevy::image::{ImageSampler, ImageSamplerDescriptor};
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use noise::{NoiseFn, Perlin};

use ymir_core::grid::GridF32;

use crate::state::{TerrainData, TerrainStats, ViewMode, ViewState};

pub struct TerrainViewPlugin;

impl Plugin for TerrainViewPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_terrain)
            .add_systems(Update, (rebuild_terrain_texture, update_terrain_stats));
    }
}

#[derive(Component)]
pub struct TerrainSprite;

#[derive(Resource)]
pub struct TerrainImageHandle(pub Handle<Image>);

// ── Hypsometric color ramp ───────────────────────────────────────────────

struct ColorStop {
    altitude: f32,
    r: u8,
    g: u8,
    b: u8,
}

const HYPSOMETRIC_RAMP: &[ColorStop] = &[
    ColorStop { altitude: -500.0, r: 0x0d, g: 0x2b, b: 0x42 },
    ColorStop { altitude: -200.0, r: 0x0d, g: 0x2b, b: 0x42 },
    ColorStop { altitude: 0.0, r: 0x1a, g: 0x4a, b: 0x6e },
    ColorStop { altitude: 1.0, r: 0x4a, g: 0x7a, b: 0x3a },
    ColorStop { altitude: 100.0, r: 0x4a, g: 0x7a, b: 0x3a },
    ColorStop { altitude: 300.0, r: 0x6b, g: 0x8f, b: 0x4a },
    ColorStop { altitude: 800.0, r: 0x8a, g: 0x7a, b: 0x5a },
    ColorStop { altitude: 1500.0, r: 0x9a, g: 0x8a, b: 0x70 },
    ColorStop { altitude: 2500.0, r: 0xb0, g: 0xa8, b: 0x90 },
    ColorStop { altitude: 4000.0, r: 0xd0, g: 0xc8, b: 0xc0 },
];

pub fn hypsometric_color(altitude_m: f32) -> [u8; 3] {
    if altitude_m <= HYPSOMETRIC_RAMP[0].altitude {
        let s = &HYPSOMETRIC_RAMP[0];
        return [s.r, s.g, s.b];
    }
    for i in 1..HYPSOMETRIC_RAMP.len() {
        let lo = &HYPSOMETRIC_RAMP[i - 1];
        let hi = &HYPSOMETRIC_RAMP[i];
        if altitude_m <= hi.altitude {
            let t = (altitude_m - lo.altitude) / (hi.altitude - lo.altitude);
            return [lerp_u8(lo.r, hi.r, t), lerp_u8(lo.g, hi.g, t), lerp_u8(lo.b, hi.b, t)];
        }
    }
    let s = HYPSOMETRIC_RAMP.last().unwrap();
    [s.r, s.g, s.b]
}

fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    let v = a as f32 + (b as f32 - a as f32) * t;
    v.clamp(0.0, 255.0) as u8
}

// ── Hillshade ────────────────────────────────────────────────────────────

fn hillshade_factor(gx: f32, gy: f32) -> f32 {
    // Light from NW: azimuth 315° → math angle 135°, elevation 45°
    let az = std::f32::consts::FRAC_PI_4 * 3.0;
    let alt = std::f32::consts::FRAC_PI_4;

    let light_x = az.cos() * alt.cos();
    let light_y = az.sin() * alt.cos();
    let light_z = alt.sin();

    let nx = -gx;
    let ny = -gy;
    let nz = 1.0;
    let len = (nx * nx + ny * ny + nz * nz).sqrt();

    let dot = (nx * light_x + ny * light_y + nz * light_z) / len;
    dot.clamp(0.3, 1.0)
}

// ── Slope heatmap ────────────────────────────────────────────────────────

fn slope_color(gx: f32, gy: f32) -> [u8; 3] {
    let [r, g, b, _] = crate::visualization::colormap::slope_color(gx, gy);
    [r, g, b]
}

// ── Test terrain generation ──────────────────────────────────────────────

fn generate_test_terrain() -> (GridF32, f32) {
    let size = 512;
    let max_elev = 4000.0_f32;
    let mut data = vec![0.0_f32; size * size];
    let perlin = Perlin::new(42);

    let cx = size as f32 / 2.0;
    let cy = size as f32 / 2.0;
    let max_dist = cx * 0.8;

    for y in 0..size {
        for x in 0..size {
            let dx = x as f32 - cx;
            let dy = y as f32 - cy;
            let dist = (dx * dx + dy * dy).sqrt();

            let cone = 1.0 - (dist / max_dist).clamp(0.0, 1.5);
            let nx = x as f64 / 80.0;
            let ny = y as f64 / 80.0;
            let n1 = perlin.get([nx, ny]) as f32 * 0.25;
            let n2 = perlin.get([nx * 2.5, ny * 2.5]) as f32 * 0.10;

            data[y * size + x] = (cone + n1 + n2) * max_elev;
        }
    }

    (GridF32::from_vec(size, size, data), max_elev)
}

// ── Systems ──────────────────────────────────────────────────────────────

fn setup_terrain(mut commands: Commands, mut images: ResMut<Assets<Image>>) {
    let (heightmap, max_elevation) = if let Ok(entries) = std::fs::read_dir("assets/maps") {
        let png = entries
            .filter_map(|e| e.ok())
            .find(|e| e.path().extension().is_some_and(|ext| ext == "png"));
        if let Some(entry) = png {
            match GridF32::load_png(&entry.path()) {
                Ok(grid) => {
                    info!("Loaded heightmap: {:?}", entry.path());
                    let max_elev = 4000.0;
                    let scaled = GridF32::from_vec(
                        grid.width,
                        grid.height,
                        grid.data.iter().map(|&v| v * max_elev).collect(),
                    );
                    (scaled, max_elev)
                }
                Err(e) => {
                    warn!("Failed to load heightmap: {e}. Using test terrain.");
                    generate_test_terrain()
                }
            }
        } else {
            generate_test_terrain()
        }
    } else {
        generate_test_terrain()
    };

    let w = heightmap.width as u32;
    let h = heightmap.height as u32;

    let stats = compute_stats(&heightmap, 40.0);
    commands.insert_resource(stats);

    let pixel_data = vec![0u8; (w * h * 4) as usize];
    let mut image = Image::new(
        Extent3d { width: w, height: h, depth_or_array_layers: 1 },
        TextureDimension::D2,
        pixel_data,
        TextureFormat::Rgba8UnormSrgb,
        bevy::asset::RenderAssetUsages::RENDER_WORLD | bevy::asset::RenderAssetUsages::MAIN_WORLD,
    );
    image.sampler = ImageSampler::Descriptor(ImageSamplerDescriptor::nearest());
    let image_handle = images.add(image);

    commands.spawn((
        Sprite {
            image: image_handle.clone(),
            custom_size: Some(Vec2::new(w as f32, h as f32)),
            ..default()
        },
        TerrainSprite,
    ));

    commands.insert_resource(TerrainImageHandle(image_handle));
    commands.insert_resource(TerrainData { heightmap, max_elevation, dirty: true });
}

fn rebuild_terrain_texture(
    mut terrain: ResMut<TerrainData>,
    view_state: Res<ViewState>,
    view_mode: Res<State<ViewMode>>,
    image_handle: Option<Res<TerrainImageHandle>>,
    mut images: ResMut<Assets<Image>>,
) {
    let Some(handle) = image_handle else { return };
    if !terrain.dirty && !view_state.is_changed() && !view_mode.is_changed() {
        return;
    }

    let Some(image) = images.get_mut(&handle.0) else { return };

    let hm = &terrain.heightmap;
    let w = hm.width;
    let h = hm.height;
    let data = image.data.as_mut().unwrap();
    let mode = *view_mode.get();
    let apply_hillshade = view_state.overlays.hillshade && mode == ViewMode::Altitude;

    for y in 0..h {
        for x in 0..w {
            let alt = hm.data[y * w + x];
            let idx = (y * w + x) * 4;

            let [r, g, b] = match mode {
                ViewMode::Altitude => hypsometric_color(alt),
                ViewMode::Slope => {
                    let (gx, gy) = hm.gradient_at(x, y);
                    slope_color(gx, gy)
                }
                _ => hypsometric_color(alt),
            };

            if apply_hillshade {
                let (gx, gy) = hm.gradient_at(x, y);
                let shade = hillshade_factor(gx, gy);
                data[idx] = (r as f32 * shade) as u8;
                data[idx + 1] = (g as f32 * shade) as u8;
                data[idx + 2] = (b as f32 * shade) as u8;
            } else {
                data[idx] = r;
                data[idx + 1] = g;
                data[idx + 2] = b;
            }
            data[idx + 3] = 255;
        }
    }

    terrain.dirty = false;
}

fn update_terrain_stats(
    terrain: Res<TerrainData>,
    gen_params: Res<crate::state::GenerationParamsUi>,
    mut stats: ResMut<TerrainStats>,
) {
    if !terrain.is_changed() {
        return;
    }
    *stats = compute_stats(&terrain.heightmap, gen_params.meters_per_pixel);
}

pub fn compute_stats(hm: &GridF32, meters_per_pixel: f32) -> TerrainStats {
    let land_count = hm.data.iter().filter(|&&v| v > 0.0).count();
    TerrainStats {
        grid_width: hm.width,
        grid_height: hm.height,
        meters_per_pixel,
        peak_altitude: hm.max(),
        min_altitude: hm.min(),
        land_ratio: land_count as f32 / hm.data.len() as f32,
        river_segments: 0,
        lake_count: 0,
    }
}
