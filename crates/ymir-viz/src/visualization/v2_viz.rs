//! Step 8.6 Phase 5 — v2 raster visualisation.
//!
//! Spawns a sprite that renders the currently-selected field from the
//! latest `V2SolverBridge::state == Completed` snapshot. Field
//! selection comes from `V2VizState::field` (driven by the parameter
//! panel dropdown). Colormap dispatch matches the issue D5 spec:
//! S̃ / age / cratonic on linear scales, ε̇_II / |v| on log scales.

use bevy::image::ImageSampler;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

use super::colormap::{
    age_colormap, cratonic_grayscale, hypsometric_colormap, log_hot, log_normalize,
};
use super::overlay;

use crate::bridge::v2::{V2FinalState, V2RunSpec, V2RunState, V2SolverBridge};

/// Currently displayed field (D5 dropdown). All variants are populated
/// by [`crate::bridge::v2::V2FinalState`] except where the corresponding
/// mechanism was disabled in the run config — `Age` and `Cratonic`
/// fall back to a flat zero canvas if the field is `None`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum V2Field {
    #[default]
    SThickness,
    Age,
    Cratonic,
    StrainRate,
    VelocityMagnitude,
}

impl V2Field {
    pub const ALL: &'static [V2Field] = &[
        V2Field::SThickness,
        V2Field::Age,
        V2Field::Cratonic,
        V2Field::StrainRate,
        V2Field::VelocityMagnitude,
    ];

    pub fn label(self) -> &'static str {
        match self {
            V2Field::SThickness => "S̃ (crustal thickness)",
            V2Field::Age => "Age field A",
            V2Field::Cratonic => "Cratonic factor",
            V2Field::StrainRate => "ε̇_II (log)",
            V2Field::VelocityMagnitude => "|v| (log)",
        }
    }

    /// Colorbar caption for the side legend (units / scale hint).
    pub fn legend_caption(self) -> &'static str {
        match self {
            V2Field::SThickness => "S̃: 0.2 (oceanic) → 1.5+ (collision)",
            V2Field::Age => "A: 0 (reset) → init_max + run_time",
            V2Field::Cratonic => "f: 0 (mobile) → 1 (cratonic core)",
            V2Field::StrainRate => "ε̇_II: 1e-3 → 1e1 (log)",
            V2Field::VelocityMagnitude => "|v|: 1e-5 → 10 (log)",
        }
    }
}

/// Marker for the v2 viz sprite. Distinct from [`super::render::SolverTerrainSprite`]
/// (legacy) so both can coexist during the parallel-bridges phase
/// without z-order fights.
#[derive(Component)]
pub struct V2VizSprite;

#[derive(Resource)]
pub struct V2VizState {
    pub field: V2Field,
    pub texture_handle: Option<Handle<Image>>,
    /// `(grid_nx, grid_ny, field_variant_index, phase_tag, step,
    /// overlay_bits)` for the last rendered frame. Triggers a
    /// re-render when any component differs (dropdown change,
    /// run-completion arrival, mid-run Progress event, overlay
    /// toggle change).
    ///
    /// `phase_tag` discriminates Running (= 1) vs Completed (= 2) so
    /// the final post-run repaint fires even if the final step
    /// counter equals the last Progress step counter.
    /// `overlay_bits` packs `show_voronoi_boundaries | (show_velocity_vectors << 1)`.
    pub last_signature: Option<(usize, usize, u8, u8, u32, u8)>,
    /// Phase 6 — set by the UI Capture button. Consumed (cleared)
    /// by `handle_v2_screenshot` after the PNG is written. The
    /// `bool` arm is "save the currently displayed field"; when we
    /// later add a "save all 5 fields" button it becomes a richer
    /// enum.
    pub capture_requested: bool,
    /// Phase 6 — last screenshot status surfaced to the UI as a
    /// toast-equivalent line below the capture button.
    pub last_capture: Option<Result<std::path::PathBuf, String>>,
    /// Phase 8b — toggle the Voronoï plate-boundary overlay (drawn
    /// in black over the field colour). Default `false` so the
    /// existing screenshots stay unchanged unless the user opts in.
    pub show_voronoi_boundaries: bool,
    /// Phase 8b — toggle the per-plate velocity-vector overlay
    /// (yellow arrow at each plate's centroid, length proportional
    /// to mean per-plate velocity). Default `false`.
    pub show_velocity_vectors: bool,
}

impl V2VizState {
    fn overlay_bits(&self) -> u8 {
        let mut b = 0u8;
        if self.show_voronoi_boundaries {
            b |= 0b01;
        }
        if self.show_velocity_vectors {
            b |= 0b10;
        }
        b
    }

    fn any_overlay(&self) -> bool {
        self.show_voronoi_boundaries || self.show_velocity_vectors
    }
}

impl Default for V2VizState {
    fn default() -> Self {
        Self {
            field: V2Field::default(),
            texture_handle: None,
            last_signature: None,
            capture_requested: false,
            last_capture: None,
            show_voronoi_boundaries: false,
            show_velocity_vectors: false,
        }
    }
}

const V2_SPRITE_BASE_SIZE: f32 = 600.0;

fn sprite_size_for(grid_width: usize, grid_height: usize) -> Vec2 {
    if grid_width == 0 || grid_height == 0 {
        return Vec2::splat(V2_SPRITE_BASE_SIZE);
    }
    let longer = grid_width.max(grid_height) as f32;
    Vec2::new(
        V2_SPRITE_BASE_SIZE * grid_width as f32 / longer,
        V2_SPRITE_BASE_SIZE * grid_height as f32 / longer,
    )
}

fn setup_v2_sprite(mut commands: Commands, mut images: ResMut<Assets<Image>>) {
    // Start with a 1×1 black placeholder; first run completion resizes
    // it to the actual grid dimensions.
    let mut image = Image::new_fill(
        Extent3d { width: 1, height: 1, depth_or_array_layers: 1 },
        TextureDimension::D2,
        &[10, 10, 10, 255],
        TextureFormat::Rgba8UnormSrgb,
        bevy::asset::RenderAssetUsages::MAIN_WORLD | bevy::asset::RenderAssetUsages::RENDER_WORLD,
    );
    image.sampler = ImageSampler::nearest();
    let handle = images.add(image);

    commands.insert_resource(V2VizState {
        field: V2Field::default(),
        texture_handle: Some(handle.clone()),
        last_signature: None,
        capture_requested: false,
        last_capture: None,
        show_voronoi_boundaries: false,
        show_velocity_vectors: false,
    });

    // z = 10 puts the v2 sprite above any leftover legacy sprite that
    // SolverVisualizationPlugin might have spawned at z=0. Sunset
    // (Phase 8) will remove the legacy sprite so this stays harmless.
    commands.spawn((
        Sprite { image: handle, custom_size: Some(sprite_size_for(1, 1)), ..default() },
        Transform::from_xyz(0.0, 0.0, 10.0),
        V2VizSprite,
    ));
}

#[allow(clippy::too_many_arguments)]
fn update_v2_texture(
    bridge: Res<V2SolverBridge>,
    mut viz: ResMut<V2VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<V2VizSprite>>,
) {
    // Step 8.6 follow-up — render either the most recent peek-state
    // (mid-run, Running with Some peek) or the final-state (Completed).
    // Idle / Failed / Running-without-peek paint nothing.
    let (state_ref, phase_tag, step_counter, spec_ref) = match &bridge.state {
        V2RunState::Completed { final_state, spec, .. } => {
            (final_state.as_ref(), 2u8, 0u32, Some(spec))
        }
        V2RunState::Running { peek_state: Some(peek), step, spec, .. } => {
            (peek.as_ref(), 1u8, *step as u32, Some(spec))
        }
        _ => return,
    };

    let nx = state_ref.nx;
    let ny = state_ref.ny;
    let field_idx = viz.field as u8;
    let overlay_bits = viz.overlay_bits();
    let signature = (nx, ny, field_idx, phase_tag, step_counter, overlay_bits);
    if viz.last_signature == Some(signature) {
        return;
    }

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

    // Phase 8b — produce field colours via the unified `field_to_rgba`
    // helper (same code path as PNG screenshots). Then overlay the
    // Voronoï boundaries / velocity arrows on the resulting RGBA8
    // buffer if the user has either toggle on. Final push: bulk-copy
    // into the texture's data buffer.
    let mut rgba = field_to_rgba_buf(state_ref, viz.field);

    if viz.any_overlay() {
        if let Some(spec) = spec_ref {
            let plate_id = regenerate_plate_id(spec, nx, ny);
            if viz.show_voronoi_boundaries {
                overlay::draw_voronoi_boundaries(&mut rgba, nx, ny, &plate_id);
            }
            if viz.show_velocity_vectors {
                overlay::draw_velocity_vectors(
                    &mut rgba,
                    nx,
                    ny,
                    &state_ref.vx,
                    &state_ref.vy,
                    &plate_id,
                    VELOCITY_ARROW_SCALE_CELLS,
                );
            }
        }
    }

    if let Some(data) = image.data.as_mut() {
        if data.len() == rgba.len() {
            data.copy_from_slice(&rgba);
        } else {
            *data = rgba;
        }
    }
    viz.last_signature = Some(signature);
}

/// Phase 8b — visible arrow scale in cells per unit non-dim velocity.
/// Picked empirically: at 64² with `peak|v̄_plate| ≈ 5` (active_medley
/// regime) this gives head separation ≈ 40 cells, ≈ 60% of the grid
/// — long enough to read direction, short enough to fit. Phase 8d will
/// expose this knob in the UI.
const VELOCITY_ARROW_SCALE_CELLS: f64 = 8.0;

fn regenerate_plate_id(spec: &V2RunSpec, nx: usize, ny: usize) -> Vec<u16> {
    use ymir_core::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};
    let cfg = VoronoiConfig {
        num_plates: spec.num_plates,
        continental_ratio: spec.continental_ratio,
    };
    let plates = generate_voronoi(nx, ny, &cfg, spec.seed);
    plates.plate_id.data().to_vec()
}

/// Internal sibling of [`field_to_rgba`] that returns just the buffer
/// (the `(nx, ny)` tuple is already known in the caller). Centralises
/// the field colour logic so the live render and the screenshot path
/// produce identical pixels prior to overlay application.
fn field_to_rgba_buf(state: &V2FinalState, field: V2Field) -> Vec<u8> {
    let (_, _, rgba) = field_to_rgba(state, field);
    rgba
}

// ── Step 8.6 Phase 6 — standalone PNG rendering ────────────────────────
//
// `field_to_rgba` is a pure function (no Bevy dependency) that paints
// the given field into a row-major RGBA8 buffer with the same Y-flip
// the on-screen sprite uses (image row 0 at top). This is what the
// Phase 6 screenshot system + the Phase 6 acceptance test consume.

/// Render `field` from `state` into a row-major RGBA8 buffer of size
/// `nx · ny · 4`. The Y-axis is flipped so image row 0 corresponds to
/// the top of the rendered sprite (grid row `ny - 1`). Returns
/// `(nx, ny, rgba_buffer)`.
pub fn field_to_rgba(state: &V2FinalState, field: V2Field) -> (usize, usize, Vec<u8>) {
    let nx = state.nx;
    let ny = state.ny;
    let mut rgba = vec![0u8; nx * ny * 4];

    let vmag_buf: Vec<f64>;
    let buf: &[f64] = match field {
        V2Field::SThickness => &state.s_field,
        V2Field::Age => match state.age_field.as_ref() {
            Some(b) => b,
            None => {
                fill_disabled(&mut rgba);
                return (nx, ny, rgba);
            }
        },
        V2Field::Cratonic => match state.cratonic_factor.as_ref() {
            Some(b) => b,
            None => {
                fill_disabled(&mut rgba);
                return (nx, ny, rgba);
            }
        },
        V2Field::StrainRate => &state.strain_rate_invariant,
        V2Field::VelocityMagnitude => {
            vmag_buf = (0..nx * ny)
                .map(|k| (state.vx[k].powi(2) + state.vy[k].powi(2)).sqrt())
                .collect();
            &vmag_buf
        }
    };

    let (vmin, vmax) = match field {
        V2Field::Cratonic => (0.0, 1.0),
        // S̃ uses the harness's clamp band [s_min, s_max] = [0.1, 2.5]
        // (mass conservation bound from the time loop) as the fixed
        // visualisation range. Auto-rescaling per frame would make
        // slow drift invisible — small advection in a quiescent
        // regime would re-normalise to the same colors each step. A
        // fixed range keeps oceanic-vs-continental contrast stable
        // across the whole run.
        V2Field::SThickness => (0.0, 2.5),
        _ => {
            let mut a = f64::INFINITY;
            let mut b = f64::NEG_INFINITY;
            for &v in buf {
                if v.is_finite() {
                    if v < a { a = v; }
                    if v > b { b = v; }
                }
            }
            if !a.is_finite() || !b.is_finite() { (0.0, 1.0) } else { (a, b) }
        }
    };
    let range = (vmax - vmin).max(1e-12);

    for j in 0..ny {
        for i in 0..nx {
            let v = buf[j * nx + i];
            let rgba_pixel = match field {
                V2Field::SThickness => hypsometric_colormap(((v - vmin) / range).clamp(0.0, 1.0)),
                V2Field::Age => age_colormap(((v - vmin) / range).clamp(0.0, 1.0)),
                V2Field::Cratonic => cratonic_grayscale(v),
                // ε̇_II: log scale `[1e-3, 1e2]`. Phase 7 diagnostic on
            // active_medley showed `peak ε̇_II ≈ 33` after 30 steps,
            // so bounding at 1e1 saturated the high tail. The
            // `1e-3` floor matches the strain-rate floor in the
            // rheology law (Step 3 onward).
            V2Field::StrainRate => log_hot(log_normalize(v, 1e-3, 1e2)),
                V2Field::VelocityMagnitude => log_hot(log_normalize(v, 1e-5, 1e1)),
            };
            // Y-flip: image row 0 maps to grid row (ny - 1 - j).
            let img_row = ny - 1 - j;
            let idx = (img_row * nx + i) * 4;
            rgba[idx..idx + 4].copy_from_slice(&rgba_pixel);
        }
    }

    (nx, ny, rgba)
}

fn fill_disabled(rgba: &mut [u8]) {
    for chunk in rgba.chunks_exact_mut(4) {
        chunk[0] = 40;
        chunk[1] = 40;
        chunk[2] = 40;
        chunk[3] = 255;
    }
}

/// Encode `field` from `state` as a PNG file at `path`. Creates the
/// parent directory on demand. Returns `Err` on I/O failure or
/// invalid path.
pub fn save_field_png(
    state: &V2FinalState,
    field: V2Field,
    path: &std::path::Path,
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let (nx, ny, rgba) = field_to_rgba(state, field);
    image::save_buffer(
        path,
        &rgba,
        nx as u32,
        ny as u32,
        image::ColorType::Rgba8,
    )
    .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
}

/// Build a screenshot filename per D6: `{preset}_{field}_{ts}.png`.
/// `ts` is unix-epoch seconds; combined with the preset / field name
/// it disambiguates rapid clicks.
pub fn screenshot_filename(preset_label: &str, field: V2Field) -> String {
    let ts = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let field_tag = match field {
        V2Field::SThickness => "s",
        V2Field::Age => "age",
        V2Field::Cratonic => "cratonic",
        V2Field::StrainRate => "strain",
        V2Field::VelocityMagnitude => "vmag",
    };
    let safe_preset: String = preset_label
        .chars()
        .map(|c| if c.is_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect();
    format!("{}_{}_{}.png", safe_preset, field_tag, ts)
}

/// Phase 6 — handle the "Capture" button. When `viz.capture_requested`
/// is set and the bridge has a `Completed` run, dump the currently
/// selected field to PNG under `<run output_dir>/screenshots/`. Falls
/// back to OS temp on a None / disabled `output_dir`.
fn handle_v2_screenshot(bridge: Res<V2SolverBridge>, mut viz: ResMut<V2VizState>) {
    if !viz.capture_requested {
        return;
    }
    viz.capture_requested = false;
    let V2RunState::Completed { spec, final_state, .. } = &bridge.state else {
        viz.last_capture = Some(Err("no completed run to capture".to_string()));
        return;
    };
    let dir = if spec.output_dir.as_os_str().is_empty() {
        std::env::temp_dir().join("ymir_v2_screenshots")
    } else {
        spec.output_dir.join("screenshots")
    };
    let filename = screenshot_filename(&spec.preset_label, viz.field);
    let path = dir.join(filename);
    match save_field_png(final_state, viz.field, &path) {
        Ok(()) => viz.last_capture = Some(Ok(path)),
        Err(e) => viz.last_capture = Some(Err(format!("{}", e))),
    }
}

pub struct V2VisualizationPlugin;

impl Plugin for V2VisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_v2_sprite);
        app.add_systems(Update, (update_v2_texture, handle_v2_screenshot));
    }
}
