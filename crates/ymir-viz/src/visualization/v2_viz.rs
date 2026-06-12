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

use crate::bridge::v2::{
    V2AgeFieldSpec, V2CratonicSpec, V2FinalState, V2RunSpec, V2RunState, V2SolverBridge,
};

/// Currently displayed field (D5 dropdown). All variants are populated
/// by [`crate::bridge::v2::V2FinalState`] except where the corresponding
/// mechanism was disabled in the run config — `Age` and `Cratonic`
/// fall back to a flat zero canvas if the field is `None`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum V2Field {
    #[default]
    SThickness,
    /// Step 12 R0 — Airy-isostasy altitude derived from S̃ per-frame
    /// (`compute_isostasy(&Field2D::from_vec(s_field), IsostasyConfig
    /// ::default()).heightmap`). The native normalised heightmap is
    /// remapped piecewise-linearly so the configured sea level lands
    /// at 0.5 in the hypsometric colormap (visually balances the
    /// ocean / land bands; the colormap's blue→green transition
    /// occupies [0.2, 0.4] and the green→brown→white land ramp
    /// [0.4, 1.0], so `0.5` sits a hair above the vegetation band —
    /// the design choice in the issue plan).
    Altitude,
    Age,
    Cratonic,
    StrainRate,
    VelocityMagnitude,
    /// Step 8.6 follow-up — Sobel gradient magnitude of the S̃ field
    /// (periodic boundaries). Restores the legacy "slope" view that
    /// was dropped during Phase 8h sunset; useful for spotting
    /// boundary deformation at a glance.
    Slope,
}

impl V2Field {
    pub const ALL: &'static [V2Field] = &[
        V2Field::SThickness,
        V2Field::Altitude,
        V2Field::Age,
        V2Field::Cratonic,
        V2Field::StrainRate,
        V2Field::VelocityMagnitude,
        V2Field::Slope,
    ];

    pub fn label(self) -> &'static str {
        match self {
            V2Field::SThickness => "S̃ (crustal thickness)",
            V2Field::Altitude => "Altitude (post-isostasy)",
            V2Field::Age => "Age field A",
            V2Field::Cratonic => "Cratonic factor",
            V2Field::StrainRate => "ε̇_II (log)",
            V2Field::VelocityMagnitude => "|v| (log)",
            V2Field::Slope => "Slope |∇S̃|",
        }
    }

    /// Colorbar caption for the side legend (units / scale hint).
    pub fn legend_caption(self) -> &'static str {
        match self {
            V2Field::SThickness => "S̃: 0.2 (oceanic) → 1.5+ (collision)",
            V2Field::Altitude => "Altitude: deep ocean → mountain peaks (sea level @ 0.5)",
            V2Field::Age => "A: 0 (reset) → init_max + run_time",
            V2Field::Cratonic => "f: 0 (mobile) → 1 (cratonic core)",
            V2Field::StrainRate => "ε̇_II: 1e-3 → 1e1 (log)",
            V2Field::VelocityMagnitude => "|v|: 1e-5 → 10 (log)",
            V2Field::Slope => "|∇S̃|: 0 → max (per-frame)",
        }
    }
}

/// Marker for the shared phase sprite. Originally introduced for
/// the v2 raster (Step 8.6 Phase 5); now reused by the isostasy /
/// upscale / erosion / hydrology phase render systems too — only
/// one sprite ever exists, the active phase decides what's
/// painted.
#[derive(Component)]
pub struct V2VizSprite;

#[derive(Resource)]
pub struct V2VizState {
    pub field: V2Field,
    pub texture_handle: Option<Handle<Image>>,
    /// `(grid_nx, grid_ny, field_variant_index, phase_tag, step,
    /// reserved)` for the last rendered frame. Triggers a re-render
    /// when any component differs (dropdown change, run-completion
    /// arrival, mid-run Progress event).
    ///
    /// `phase_tag` discriminates Running (= 1) vs Completed (= 2) so
    /// the final post-run repaint fires even if the final step
    /// counter equals the last Progress step counter. The trailing
    /// `u8` is reserved (was `overlay_bits` pre-Phase 8h post-sunset
    /// overlay refactor; now overlays are gizmos and don't touch the
    /// texture).
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
    /// Step 11 — multiplier on the arrow length (1.0 = the legacy
    /// `VELOCITY_ARROW_SCALE_CELLS` baseline tuned for `peak|v̄_plate|
    /// ≈ 5` in the active_medley regime). Drift-driven runs at
    /// `|drift| ≈ 0.5` produce arrows ~10× shorter at default scale,
    /// so the user can dial this up to `2.0`–`8.0` to make small
    /// drifts visible. Anti-pattern (per the issue): never auto-scale
    /// to the current `max|v|` — that would make arrows look the same
    /// across runs and break visual comparison.
    pub arrow_scale: f64,
    /// Phase 8e — set by the Export button. Consumed by the
    /// `handle_v2_export` system; surfaces the result via
    /// [`Self::last_export`].
    pub export_requested: bool,
    /// Phase 8e — last export status surfaced to the UI.
    pub last_export: Option<Result<std::path::PathBuf, String>>,
    /// Phase 8e — set by the Import button. Carries the source file
    /// path. Consumed by the `handle_v2_import` system on the next
    /// frame.
    pub import_requested_path: Option<std::path::PathBuf>,
    /// Phase 8e — last import status surfaced to the UI.
    pub last_import: Option<Result<std::path::PathBuf, String>>,
    /// Step 8.6 follow-up — lazily-cached list of `output/seed*_*` dirs
    /// containing a `snapshot.json`. `None` means "rescan on next
    /// draw"; cleared after a successful export so the freshly written
    /// run shows up immediately. Replaces the legacy manual TextEdit:
    /// the UI surfaces this as an "Available runs" list with per-row
    /// Load buttons.
    pub cached_run_dirs: Option<Vec<std::path::PathBuf>>,
    /// Step 8.6 follow-up — pre-Run preview of the editable spec.
    /// Recomputed on `V2EditableSpec` change while the bridge is in
    /// `Idle`; cleared / replaced once a run starts. Lets the user
    /// pick a configuration graphically before clicking Run.
    pub preview: Option<V2Preview>,
}

/// A snapshot of "what the run would start from" given the user's
/// current editable spec. Holds the same shape as a real run's
/// `final_state`, so the rendering pipeline (`update_v2_texture` +
/// `render_v2_overlay_gizmos`) treats it identically to a
/// `Completed` / `Imported` state.
#[derive(Clone)]
pub struct V2Preview {
    pub state: Box<V2FinalState>,
    pub spec: V2RunSpec,
    /// Cached hash of the spec used to build `state`; the
    /// `update_v2_preview` system compares this against the current
    /// `V2EditableSpec` hash to decide whether to recompute.
    pub spec_hash: u64,
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
            arrow_scale: 1.0,
            export_requested: false,
            last_export: None,
            import_requested_path: None,
            last_import: None,
            cached_run_dirs: None,
            preview: None,
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
        arrow_scale: 1.0,
        export_requested: false,
        last_export: None,
        import_requested_path: None,
        last_import: None,
        cached_run_dirs: None,
        preview: None,
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
    active: Res<crate::pipeline::ActivePhase>,
    mut viz: ResMut<V2VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<V2VizSprite>>,
) {
    // Only the Tectonics phase paints the v2 raster; other phases
    // own the sprite when they're active.
    if active.0 != crate::pipeline::PipelinePhase::Tectonics {
        return;
    }
    // Step 8.6 follow-up — render the most recent state available:
    //   Running peek_state (mid-run) → phase_tag = 1
    //   Completed final_state         → phase_tag = 2
    //   Imported final_state          → phase_tag = 3
    //   Idle + preview state          → phase_tag = 4 (pre-Run)
    // Idle without preview, Failed, Running-without-peek paint
    // nothing. The owned tuple sidesteps borrow-checker friction
    // between the source borrow and the trailing
    // `viz.last_signature` mutation; clone cost is one V2FinalState
    // (~36 KB at 64²) per re-render, well under the per-frame budget.
    let target = render_target(&bridge.state, viz.preview.as_ref());
    let Some((state, phase_tag, step_counter, _spec)) = target else { return };

    let nx = state.nx;
    let ny = state.ny;
    let field_idx = viz.field as u8;
    // Trailing u8 reserved (was Phase 8b `overlay_bits` before
    // overlays moved to gizmos at display resolution).
    let signature = (nx, ny, field_idx, phase_tag, step_counter, 0u8);
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

    // Phase 8h post-sunset — produce field colours via the unified
    // `field_to_rgba` helper (same code path as PNG screenshots) and
    // bulk-copy into the texture data buffer. Overlays are NOT burned
    // here anymore — they're drawn in a separate Bevy `Gizmos` pass at
    // display resolution (`render_v2_overlay_gizmos`). The export PNG
    // path keeps the buffer-burn (`overlay::draw_*` in Phase 8g test).
    let rgba = field_to_rgba_buf(&state, viz.field);

    if let Some(data) = image.data.as_mut() {
        if data.len() == rgba.len() {
            data.copy_from_slice(&rgba);
        } else {
            *data = rgba;
        }
    }
    viz.last_signature = Some(signature);
}

/// Resolve the `(state, phase_tag, step_counter, spec)` that the
/// rendering pipeline should paint this frame. Returns owned values
/// so the caller can mutate `V2VizState` afterwards without borrow
/// conflicts. Returns `None` for `Idle` without a preview, `Failed`,
/// or `Running` without a peek_state (precompute in flight).
fn render_target(
    bridge_state: &V2RunState,
    preview: Option<&V2Preview>,
) -> Option<(V2FinalState, u8, u32, V2RunSpec)> {
    match bridge_state {
        V2RunState::Completed { final_state, spec, .. } => {
            Some((final_state.as_ref().clone(), 2u8, 0u32, spec.clone()))
        }
        V2RunState::Imported { final_state, spec, .. } => {
            Some((final_state.as_ref().clone(), 3u8, 0u32, spec.clone()))
        }
        V2RunState::Running { peek_state: Some(peek), step, spec, .. } => {
            Some((peek.as_ref().clone(), 1u8, *step as u32, spec.clone()))
        }
        V2RunState::Idle => preview.map(|p| {
            // Use the spec_hash low 32 bits as the step_counter slot
            // so changes to the editable spec invalidate the
            // signature and force a re-render.
            (
                p.state.as_ref().clone(),
                4u8,
                p.spec_hash as u32,
                p.spec.clone(),
            )
        }),
        _ => None,
    }
}

/// Phase 8h post-sunset — overlays at display resolution.
///
/// Voronoï plate boundaries (black) and per-plate velocity arrows
/// (yellow) are rendered with Bevy `Gizmos` in world coords on top of
/// the v2 sprite. The sprite occupies a square or near-square region
/// centered at the world origin (size = [`sprite_size_for`]); the
/// overlay maps grid cells to that world rect so the visual density
/// of the overlay scales with the displayed sprite (≈ 600 px on the
/// long axis), not the underlying simulation grid.
///
/// Toggles default off (`V2VizState.show_*`); the system early-exits
/// when both are off so it has zero cost in the steady state.
fn render_v2_overlay_gizmos(
    bridge: Res<V2SolverBridge>,
    active: Res<crate::pipeline::ActivePhase>,
    viz: Res<V2VizState>,
    mut gizmos: Gizmos,
) {
    // Overlays follow the v2 sprite — only show on the Tectonics
    // phase (the other phases render different fields where the
    // Voronoï tessellation isn't a meaningful annotation).
    if active.0 != crate::pipeline::PipelinePhase::Tectonics {
        return;
    }
    if !viz.show_voronoi_boundaries && !viz.show_velocity_vectors {
        return;
    }

    let (state_ref, spec_ref) = match &bridge.state {
        V2RunState::Completed { final_state, spec, .. } => (final_state.as_ref(), spec),
        V2RunState::Imported { final_state, spec, .. } => (final_state.as_ref(), spec),
        V2RunState::Running { peek_state: Some(peek), spec, .. } => (peek.as_ref(), spec),
        V2RunState::Idle => match viz.preview.as_ref() {
            Some(p) => (p.state.as_ref(), &p.spec),
            None => return,
        },
        _ => return,
    };

    let nx = state_ref.nx;
    let ny = state_ref.ny;
    if nx == 0 || ny == 0 {
        return;
    }

    let sprite_size = sprite_size_for(nx, ny);
    let w = sprite_size.x;
    let h = sprite_size.y;
    let cell_w = w / nx as f32;
    let cell_h = h / ny as f32;
    // Grid cell (i, j) — j=0 at the bottom in world coords (matches the
    // sprite's Y-flip: image row 0 = grid row ny-1 at the top of world).
    let world_corner = |i: i32, j: i32| -> Vec2 {
        Vec2::new(-w * 0.5 + (i as f32) * cell_w, -h * 0.5 + (j as f32) * cell_h)
    };
    let world_centre_for_cell = |ci: f64, cj: f64| -> Vec2 {
        Vec2::new(
            -w * 0.5 + (ci as f32 + 0.5) * cell_w,
            -h * 0.5 + (cj as f32 + 0.5) * cell_h,
        )
    };

    let plate_id = regenerate_plate_id(spec_ref, nx, ny);

    if viz.show_voronoi_boundaries {
        let boundary_color = Color::srgba(0.05, 0.05, 0.05, 1.0);
        for j in 0..ny {
            for i in 0..nx {
                let id = plate_id[j * nx + i];
                let ip = (i + 1) % nx;
                let jp = (j + 1) % ny;
                if plate_id[j * nx + ip] != id {
                    let x = -w * 0.5 + (i + 1) as f32 * cell_w;
                    gizmos.line_2d(
                        Vec2::new(x, world_corner(0, j as i32).y),
                        Vec2::new(x, world_corner(0, j as i32 + 1).y),
                        boundary_color,
                    );
                }
                if plate_id[jp * nx + i] != id {
                    let y = -h * 0.5 + (j + 1) as f32 * cell_h;
                    gizmos.line_2d(
                        Vec2::new(world_corner(i as i32, 0).x, y),
                        Vec2::new(world_corner(i as i32 + 1, 0).x, y),
                        boundary_color,
                    );
                }
            }
        }
    }

    if viz.show_velocity_vectors {
        let arrow_color = Color::srgba(1.0, 0.94, 0.0, 1.0);
        // Periodic-aware centroid + mean velocity per plate (mirrors
        // `overlay::draw_velocity_vectors` but emits gizmos in world
        // coords instead of pixels in the texture buffer).
        use std::collections::HashMap;
        use std::f64::consts::PI;

        #[derive(Default)]
        struct Acc {
            cos_x: f64,
            sin_x: f64,
            cos_y: f64,
            sin_y: f64,
            vx: f64,
            vy: f64,
            count: usize,
        }
        let mut by_plate: HashMap<u16, Acc> = HashMap::new();
        let nx_f = nx as f64;
        let ny_f = ny as f64;
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                let entry = by_plate.entry(plate_id[idx]).or_default();
                let theta_x = 2.0 * PI * (i as f64 + 0.5) / nx_f;
                let theta_y = 2.0 * PI * (j as f64 + 0.5) / ny_f;
                entry.cos_x += theta_x.cos();
                entry.sin_x += theta_x.sin();
                entry.cos_y += theta_y.cos();
                entry.sin_y += theta_y.sin();
                entry.vx += state_ref.vx[idx];
                entry.vy += state_ref.vy[idx];
                entry.count += 1;
            }
        }

        for acc in by_plate.values() {
            if acc.count == 0 {
                continue;
            }
            let cx = (acc.sin_x.atan2(acc.cos_x) / (2.0 * PI)) * nx_f;
            let cy = (acc.sin_y.atan2(acc.cos_y) / (2.0 * PI)) * ny_f;
            let cx = cx.rem_euclid(nx_f);
            let cy = cy.rem_euclid(ny_f);

            let mvx = acc.vx / acc.count as f64;
            let mvy = acc.vy / acc.count as f64;
            // Step 11 — `viz.arrow_scale` is a user-facing
            // multiplier on the constant `VELOCITY_ARROW_SCALE_CELLS`
            // baseline. Scaling stays *fixed per-frame* (proportional
            // to the velocity magnitude, not normalised by max|v| —
            // see the issue's anti-pattern note: auto-norm makes
            // arrows look identical across runs and breaks visual
            // comparison).
            let scale = VELOCITY_ARROW_SCALE_CELLS * viz.arrow_scale;
            let head_dx = mvx * scale;
            let head_dy = mvy * scale;
            let head_len_cells = (head_dx * head_dx + head_dy * head_dy).sqrt();
            if head_len_cells < 1.0 {
                continue;
            }

            let p0 = world_centre_for_cell(cx, cy);
            let p1 = world_centre_for_cell(cx + head_dx, cy + head_dy);
            gizmos.line_2d(p0, p1, arrow_color);

            // Two arrowhead barbs at ±150° from the shaft. Length
            // capped at 4 cells in cell-space → translated to world.
            let theta = mvy.atan2(mvx);
            let barb_cells = (head_len_cells * 0.3).min(4.0);
            for &phi_off in &[5.0 * PI / 6.0, -5.0 * PI / 6.0] {
                let phi = theta + phi_off;
                let bx = cx + head_dx + barb_cells * phi.cos();
                let by = cy + head_dy + barb_cells * phi.sin();
                let pb = world_centre_for_cell(bx, by);
                gizmos.line_2d(p1, pb, arrow_color);
            }
        }
    }
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
    let slope_buf: Vec<f64>;
    let altitude_buf: Vec<f64>;
    let buf: &[f64] = match field {
        V2Field::SThickness => &state.s_field,
        V2Field::Altitude => {
            altitude_buf = compute_altitude_buf(state);
            &altitude_buf
        }
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
        V2Field::Slope => {
            slope_buf = sobel_magnitude_periodic(&state.s_field, nx, ny);
            &slope_buf
        }
    };

    let (vmin, vmax) = match field {
        V2Field::Cratonic => (0.0, 1.0),
        // Step 12 R0 — `compute_altitude_buf` already piecewise-remaps
        // the isostatic heightmap into `[0, 1]` with sea level at 0.5;
        // bounds are therefore fixed (no per-frame auto-rescale).
        V2Field::Altitude => (0.0, 1.0),
        V2Field::Slope => {
            // Auto-scale: slope range varies wildly between quiescent
            // (~1e-3) and active (~1e0) regimes. Auto per-frame keeps
            // contrast useful across both.
            let mut a = f64::INFINITY;
            let mut b = f64::NEG_INFINITY;
            for &v in buf {
                if v.is_finite() {
                    if v < a { a = v; }
                    if v > b { b = v; }
                }
            }
            if !a.is_finite() || !b.is_finite() {
                (0.0, 1.0)
            } else {
                (a, b)
            }
        }
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
                V2Field::Altitude => hypsometric_colormap(v.clamp(0.0, 1.0)),
                V2Field::Age => age_colormap(((v - vmin) / range).clamp(0.0, 1.0)),
                V2Field::Cratonic => cratonic_grayscale(v),
                // ε̇_II: log scale `[1e-3, 1e2]`. Phase 7 diagnostic on
            // active_medley showed `peak ε̇_II ≈ 33` after 30 steps,
            // so bounding at 1e1 saturated the high tail. The
            // `1e-3` floor matches the strain-rate floor in the
            // rheology law (Step 3 onward).
            V2Field::StrainRate => log_hot(log_normalize(v, 1e-3, 1e2)),
                V2Field::VelocityMagnitude => log_hot(log_normalize(v, 1e-5, 1e1)),
                V2Field::Slope => log_hot(((v - vmin) / range).clamp(0.0, 1.0)),
            };
            // Y-flip: image row 0 maps to grid row (ny - 1 - j).
            let img_row = ny - 1 - j;
            let idx = (img_row * nx + i) * 4;
            rgba[idx..idx + 4].copy_from_slice(&rgba_pixel);
        }
    }

    (nx, ny, rgba)
}

/// Step 12 R0 — build the per-frame altitude buffer from the v2 S̃
/// state via the legacy Airy-isostasy helper. Returns a `Vec<f64>` of
/// length `nx · ny` in `[0, 1]` with the configured sea level remapped
/// to `0.5` (visually balances ocean / land bands inside the
/// hypsometric colormap whose blue→green transition is `[0.2, 0.4]`).
///
/// Default `IsostasyConfig` (no smoothing override): `sea_level_fraction
/// = 0.4` ≈ 30 % land / 70 % ocean over the field's actual S̃ range.
/// Smoothing sigma stays at the legacy default `2.0` so the altitude
/// view mirrors what the Isostasy phase will compute next in the
/// pipeline — same numerical contract, no double-rendering surprise.
fn compute_altitude_buf(state: &V2FinalState) -> Vec<f64> {
    use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
    use ymir_core::tectonics::solver::field::Field2D;

    let field = Field2D::from_vec(state.nx, state.ny, state.s_field.clone());
    // Step 12 R7.A.1 — disable the post-isostasy Gaussian blur for the
    // *viz* path. `IsostasyConfig::default()` carries
    // `altitude_smoothing_sigma = 2.0`, which is the right call for
    // downstream pipelines (hydraulic erosion stability, final export)
    // but in the diagnostic viz it erases any S̃ structure narrower
    // than ~5 cells — including the orogenic ridge at 32² × 64²
    // (σ_ridge < 2 cells). Disabling the blur locally keeps the
    // viz altitude faithful to the S̃ field's discrete sampling.
    let cfg = IsostasyConfig {
        altitude_smoothing_sigma: 0.0,
        ..IsostasyConfig::default()
    };
    let iso = compute_isostasy(&field, &cfg);
    let sea_norm = iso.sea_level_normalized as f64;
    let sea_clamped = sea_norm.clamp(1e-6, 1.0 - 1e-6);
    iso.heightmap
        .data
        .iter()
        .map(|&h| {
            let h = h as f64;
            if h <= sea_clamped {
                0.5 * h / sea_clamped
            } else {
                0.5 + 0.5 * (h - sea_clamped) / (1.0 - sea_clamped)
            }
        })
        .collect()
}

fn fill_disabled(rgba: &mut [u8]) {
    for chunk in rgba.chunks_exact_mut(4) {
        chunk[0] = 40;
        chunk[1] = 40;
        chunk[2] = 40;
        chunk[3] = 255;
    }
}

/// Step 8.6 follow-up — periodic Sobel gradient magnitude. Used to
/// compute the `Slope` view of the V2 S̃ field. The Sobel kernel
/// (3×3, weights `[1, 2, 1]` × `[-1, 0, 1]`) is normalised by `1/8`
/// so the output has roughly the same units as a centred difference.
/// All neighbours wrap modulo `(nx, ny)` so the gradient is consistent
/// with the simulation's periodic boundaries.
fn sobel_magnitude_periodic(field: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny];
    if nx < 2 || ny < 2 {
        return out;
    }
    for j in 0..ny {
        let jp = (j + 1) % ny;
        let jm = (j + ny - 1) % ny;
        for i in 0..nx {
            let ip = (i + 1) % nx;
            let im = (i + nx - 1) % nx;
            let s_jm_im = field[jm * nx + im];
            let s_jm_ip = field[jm * nx + ip];
            let s_jm_i = field[jm * nx + i];
            let s_jp_im = field[jp * nx + im];
            let s_jp_ip = field[jp * nx + ip];
            let s_jp_i = field[jp * nx + i];
            let s_j_im = field[j * nx + im];
            let s_j_ip = field[j * nx + ip];
            let gx = (-s_jm_im + s_jm_ip)
                + 2.0 * (-s_j_im + s_j_ip)
                + (-s_jp_im + s_jp_ip);
            let gy = (-s_jm_im - 2.0 * s_jm_i - s_jm_ip)
                + (s_jp_im + 2.0 * s_jp_i + s_jp_ip);
            out[j * nx + i] = (gx * gx + gy * gy).sqrt() / 8.0;
        }
    }
    out
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
        V2Field::Altitude => "altitude",
        V2Field::Age => "age",
        V2Field::Cratonic => "cratonic",
        V2Field::StrainRate => "strain",
        V2Field::VelocityMagnitude => "vmag",
        V2Field::Slope => "slope",
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
    let (spec, final_state) = match &bridge.state {
        V2RunState::Completed { spec, final_state, .. } => (spec, final_state.as_ref()),
        V2RunState::Imported { spec, final_state, .. } => (spec, final_state.as_ref()),
        _ => {
            viz.last_capture = Some(Err("no completed/imported run to capture".to_string()));
            return;
        }
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

/// Step 8.6 follow-up — full pipeline export.
///
/// Writes the v2 final state PLUS every populated downstream phase
/// (Isostasy → Upscale → Erosion → Hydrology → Lakes) to a single
/// `output/seed<seed>_<resolution>/` directory matching the
/// pre-Phase-8h-sunset layout (`01_thickness.{raw,png}`,
/// `02_altitude.*`, `03_upscaled.*`, `04_eroded.*` + `04_sediment.*`,
/// `05_*` for flow, `06_*` for lakes, `metadata.json`). The v2
/// snapshot used by the Import button lands in the same directory
/// as `snapshot.json`.
///
/// Resolution string is `<grid>` for square grids, `<W>x<H>`
/// otherwise (matches `PipelineExport::new`'s convention).
#[allow(clippy::too_many_arguments)]
fn handle_v2_export(
    bridge: Res<V2SolverBridge>,
    mut viz: ResMut<V2VizState>,
    iso_cache: Res<crate::phases::isostasy::IsostasyCache>,
    iso_params: Res<crate::phases::isostasy::IsostasyParams>,
    fbm_cache: Res<crate::phases::upscale_fbm::FbmCache>,
    fbm_params: Res<crate::phases::upscale_fbm::FbmParams>,
    erosion_cache: Res<crate::phases::erosion::ErosionCache>,
    erosion_params: Res<crate::phases::erosion::ErosionParams>,
    hydrology_cache: Res<crate::phases::hydrology::HydrologyCache>,
    hydrology_params: Res<crate::phases::hydrology::HydrologyParams>,
) {
    if !viz.export_requested {
        return;
    }
    viz.export_requested = false;

    let (spec, final_state, metrics_box, elapsed) = match &bridge.state {
        V2RunState::Completed { spec, final_state, metrics, elapsed } => {
            (spec, final_state.as_ref(), Some(metrics), *elapsed)
        }
        V2RunState::Imported { spec, final_state, elapsed, .. } => {
            (spec, final_state.as_ref(), None, *elapsed)
        }
        _ => {
            viz.last_export = Some(Err(
                "Export requires a Completed or Imported run.".to_string(),
            ));
            return;
        }
    };

    use ymir_core::erosion::hydraulic::{ErosionConfig, ErosionResult};
    use ymir_core::export::PipelineExport;
    use ymir_core::grid::GridF32;
    use ymir_core::lakes::detection::LakeConfig;
    use ymir_core::tectonics::isostasy::IsostasyConfig;
    use ymir_core::tectonics::plates::PlateConfig;
    use ymir_core::terrain::flow::{FlowConfig, RiverConfig};
    use ymir_core::terrain::upscale::FbmUpscaleConfig;

    // Synthesise a `PlateConfig` from the v2 spec for metadata.
    // The v2 solver doesn't use `PlateConfig` — it generates plates
    // from `(seed, num_plates, continental_ratio)` via Voronoï —
    // but the export's metadata format expects this struct, so we
    // populate the fields that are meaningful and zero the
    // legacy-only ones.
    let plates = PlateConfig {
        num_plates: spec.num_plates,
        num_continental_plates: ((spec.num_plates as f32) * (spec.continental_ratio as f32))
            .round() as usize,
        continental_area_factor: 1.0,
        velocity_min: 0.0,
        velocity_max: 0.0,
        grid_width: spec.grid_nx,
        grid_height: spec.grid_ny,
        boundary_smoothing_sigma: 0.0,
    };

    let output_root = std::path::Path::new("output");
    let mut export = PipelineExport::new(
        output_root,
        spec.seed,
        spec.grid_nx,
        spec.grid_ny,
        &plates,
    );
    let dir = export.dir.clone();

    let mut errors: Vec<String> = Vec::new();
    let mut saved: Vec<&'static str> = Vec::new();

    // ── 01: thickness (always available — v2 final S̃) ──
    let thickness = GridF32::from_vec(
        final_state.nx,
        final_state.ny,
        final_state.s_field.iter().map(|&v| v as f32).collect(),
    );
    match export.save_thickness(&thickness) {
        Ok(()) => saved.push("thickness"),
        Err(e) => errors.push(format!("thickness: {}", e)),
    }

    // ── 02: altitude (Isostasy) ──
    if let Some(iso) = iso_cache.result.as_ref() {
        let cfg = IsostasyConfig {
            rho_crust: iso_params.rho_crust,
            rho_mantle: iso_params.rho_mantle,
            rho_water: IsostasyConfig::default().rho_water,
            max_elevation_m: iso_params.max_elevation_m,
            max_depth_m: iso_params.max_depth_m,
            sea_level_fraction: iso_params.sea_level_fraction,
            altitude_smoothing_sigma: iso_params.altitude_smoothing_sigma,
            // v2 path: keep MinMaxFraction (Issue #141 byte-compat).
            ..IsostasyConfig::default()
        };
        match export.save_altitude(iso, &cfg) {
            Ok(()) => saved.push("altitude"),
            Err(e) => errors.push(format!("altitude: {}", e)),
        }
    }

    // ── 03: upscaled (UpscaleFbm) ──
    if let Some(fbm) = fbm_cache.result.as_ref() {
        // Drive `target_size` from the *actual* heightmap dims, not
        // the current panel slider — the user may have changed the
        // target_size dropdown after running FBM but before clicking
        // Save, in which case the slider value would mis-describe
        // the file on disk and the load path would compute wrong
        // expected dimensions (`Size mismatch` on re-import).
        let actual_target = fbm.heightmap.width.max(fbm.heightmap.height);
        let cfg = FbmUpscaleConfig {
            target_size: actual_target,
            octaves: fbm_params.octaves,
            amplitude_base: fbm_params.amplitude_base,
            amplitude_slope_factor: fbm_params.amplitude_slope_factor,
            max_anisotropy: fbm_params.max_anisotropy,
            submarine_damping: fbm_params.submarine_damping,
            domain_warp_strength: fbm_params.domain_warp_strength,
            ..FbmUpscaleConfig::default()
        };
        match export.save_upscaled(&fbm.heightmap, &cfg) {
            Ok(()) => saved.push("upscaled"),
            Err(e) => errors.push(format!("upscaled: {}", e)),
        }
    }

    // ── 04: eroded + sediment (Erosion) ──
    if let Some(eroded) = erosion_cache.result.as_ref() {
        let sea_level = iso_cache
            .result
            .as_ref()
            .map(|r| r.sea_level_normalized)
            .unwrap_or(0.4);
        let cfg = ErosionConfig {
            num_droplets: erosion_params.num_droplets,
            deposition_rate: erosion_params.deposition_rate,
            erosion_rate: erosion_params.erosion_rate,
            inertia: erosion_params.inertia,
            evaporation_rate: erosion_params.evaporation_rate,
            max_lifetime: erosion_params.max_lifetime,
            erosion_radius: erosion_params.erosion_radius,
            coastal_deposition_range: erosion_params.coastal_deposition_range,
            sea_level,
            ..ErosionConfig::default()
        };
        // `save_eroded` borrows `&ErosionResult` — clone the cache
        // result into a fresh struct so the caller keeps its
        // borrowing rights.
        let result_view = ErosionResult {
            heightmap: eroded.heightmap.clone(),
            sediment: eroded.sediment.clone(),
            stats: eroded.stats.clone(),
        };
        match export.save_eroded(&result_view, &cfg) {
            Ok(()) => saved.push("eroded+sediment"),
            Err(e) => errors.push(format!("eroded: {}", e)),
        }
    }

    // ── 05: flow + rivers (Hydrology) ──
    if let (Some(flow), Some(rivers)) =
        (hydrology_cache.flow.as_ref(), hydrology_cache.rivers.as_ref())
    {
        let sea_level = hydrology_params
            .sea_level_override
            .or_else(|| iso_cache.result.as_ref().map(|r| r.sea_level_normalized))
            .unwrap_or(0.1);
        let flow_cfg = FlowConfig { sea_level };
        let river_cfg = RiverConfig {
            stream_threshold: hydrology_params.stream_threshold,
            river_threshold: hydrology_params.river_threshold,
            major_river_threshold: hydrology_params.major_river_threshold,
        };
        match export.save_flow(flow, &flow_cfg, &river_cfg, Some(rivers)) {
            Ok(()) => saved.push("flow+rivers"),
            Err(e) => errors.push(format!("flow: {}", e)),
        }
    }

    // ── 06: lakes (Hydrology) ──
    if let Some(lakes) = hydrology_cache.lakes.as_ref() {
        let lake_cfg = LakeConfig {
            min_depth: hydrology_params.lake_min_depth,
            min_area: hydrology_params.lake_min_area,
        };
        match export.save_lakes(lakes, &lake_cfg) {
            Ok(()) => saved.push("lakes"),
            Err(e) => errors.push(format!("lakes: {}", e)),
        }
    }

    // ── snapshot.json (V2 round-trip artefact) ──
    if let Some(metrics) = metrics_box {
        let snapshot = crate::bridge::v2::V2RunSnapshot::new(
            spec.clone(),
            final_state.clone(),
            metrics,
            elapsed,
        );
        match snapshot.save(&dir.join("snapshot.json")) {
            Ok(()) => saved.push("snapshot.json"),
            Err(e) => errors.push(format!("snapshot.json: {}", e)),
        }
    }

    if errors.is_empty() {
        viz.last_export = Some(Ok(dir.clone()));
        // Force the "Available runs" list in the UI to rescan on the
        // next frame so the freshly written directory appears.
        viz.cached_run_dirs = None;
        info!(
            "[export] wrote {} artefacts to {} → {}",
            saved.len(),
            dir.display(),
            saved.join(", ")
        );
    } else {
        let combined = errors.join(" · ");
        viz.last_export = Some(Err(combined.clone()));
        warn!(
            "[export] partial failure ({} ok, {} errors) at {} — {}",
            saved.len(),
            errors.len(),
            dir.display(),
            combined
        );
    }
}

/// Step 8.6 follow-up — pre-Run preview.
///
/// When the bridge is `Idle` the user is editing a spec; this
/// system runs each frame, hashes the editable spec, and recomputes
/// the preview when the hash differs from the cached one. The
/// preview is the same shape as a real run's `final_state` (Voronoï
/// + per-cell init S̃ + initial age + initial cratonic factor) so the
/// existing rendering pipeline (`update_v2_texture` +
/// `render_v2_overlay_gizmos`) draws it identically. Velocity / ε̇_II
/// stay zero — there's no solver run yet.
///
/// Skipped while the bridge is Running / Completed / Imported /
/// Failed (the live state takes priority). The preview is left in
/// place across non-Idle transitions so the user sees it again when
/// the bridge returns to Idle (currently it doesn't, but a future
/// "reset" command could).
fn update_v2_preview(
    spec_state: Res<crate::ui::parameter_panel_v2::V2EditableSpec>,
    bridge: Res<V2SolverBridge>,
    mut viz: ResMut<V2VizState>,
) {
    if !matches!(bridge.state, V2RunState::Idle) {
        return;
    }
    let h = spec_hash(&spec_state.0);
    if let Some(p) = viz.preview.as_ref() {
        if p.spec_hash == h {
            return;
        }
    }
    let spec_clone = spec_state.0.clone();
    let state = compute_preview_state(&spec_clone);
    viz.preview = Some(V2Preview {
        state: Box::new(state),
        spec: spec_clone,
        spec_hash: h,
    });
}

/// Hash the editable spec via JSON serialisation. `serde_json`
/// produces a deterministic byte sequence for the same spec, so the
/// resulting hash is stable across frames; cheap enough to compute
/// every frame at this scale (the spec serialises to a few KB).
fn spec_hash(spec: &V2RunSpec) -> u64 {
    use std::hash::Hasher;
    let json = serde_json::to_vec(spec).expect("V2RunSpec serializes");
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    hasher.write(&json);
    hasher.finish()
}

/// Build the "what would the run start from" snapshot for the given
/// spec. Reuses the same core helpers the harness calls at run start
/// (`generate_voronoi`, `init::init_s_field`,
/// `AgeFieldState::from_initial_thickness`,
/// `build_cratonic_factor_field`); velocity / ε̇_II / boundary_flag
/// are left at zero / `None` since they only become meaningful after
/// the first solve.
fn compute_preview_state(spec: &V2RunSpec) -> V2FinalState {
    use ymir_core::tectonics_v2::age_field::{AgeFieldConfigEnabled, AgeFieldState};
    use ymir_core::tectonics_v2::boundaries::PlateType;
    use ymir_core::tectonics_v2::cratonic::factor::build_cratonic_factor_field;
    use ymir_core::tectonics_v2::cratonic::CratonicConfigEnabled;
    use ymir_core::tectonics_v2::init::{init_s_field, InitContext, PlateInitData};
    use ymir_core::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

    let nx = spec.grid_nx;
    let ny = spec.grid_ny;
    let cfg = VoronoiConfig {
        num_plates: spec.num_plates,
        continental_ratio: spec.continental_ratio,
    };
    let plates = generate_voronoi(nx, ny, &cfg, spec.seed);

    let plate_data = PlateInitData {
        plate_id: &plates.plate_id,
        plate_type: &plates.plate_type,
        seed_coords: Some(&plates.seed_coords),
    };
    let init_ctx = InitContext {
        nx,
        ny,
        seed: spec.seed,
        amplitude: spec.s_perturbation_amplitude,
        plate_data: Some(plate_data),
    };
    let s = init_s_field(spec.init_mode.into_core(), &init_ctx);

    let age_field = match &spec.age_field {
        V2AgeFieldSpec::On { continental_age_init, oceanic_age_init } => {
            let cfg = AgeFieldConfigEnabled {
                continental_age_init: *continental_age_init,
                oceanic_age_init: *oceanic_age_init,
            };
            let state = AgeFieldState::from_initial_thickness(&s, &cfg);
            Some(state.current.data().to_vec())
        }
        V2AgeFieldSpec::Off => None,
    };

    let cratonic_factor = match &spec.cratonic {
        V2CratonicSpec::On {
            cr,
            k_viscous,
            b_factor,
            smoothing_width,
            plate_area_min,
        } => {
            let cfg = CratonicConfigEnabled {
                cr: *cr,
                k_viscous: *k_viscous,
                b_factor: *b_factor,
                smoothing_width: *smoothing_width,
                plate_area_min: *plate_area_min,
                ..CratonicConfigEnabled::default()
            };
            let f = build_cratonic_factor_field(&plates, &cfg);
            Some(f.data().to_vec())
        }
        V2CratonicSpec::Off => None,
    };

    let plate_type = plates
        .plate_type
        .data()
        .iter()
        .map(|t| match t {
            PlateType::Oceanic => 0u8,
            PlateType::Continental => 1u8,
        })
        .collect();

    // Step 11 — when the user has configured PerPlate drift, populate
    // the preview's `vx, vy` with the drift field so the velocity-
    // arrow overlay shows the prescribed plate motion *before* the
    // user clicks Run. Without this the preview's velocity field
    // would stay at zeros and the user could not visually verify
    // their per-plate slider settings until launching a run.
    //
    // For `Zero` the buffers stay at zero — same as pre-Step-11
    // preview semantics.
    let (vx, vy) = match &spec.plate_kinematic {
        crate::bridge::v2::V2PlateKinematicSpec::Zero => {
            (vec![0.0; nx * ny], vec![0.0; nx * ny])
        }
        crate::bridge::v2::V2PlateKinematicSpec::PerPlate {
            velocities,
            boundary_smoothing_width,
        } => ymir_core::tectonics_v2::plate_kinematic::field::build(
            nx,
            ny,
            &plates.plate_id,
            velocities,
            *boundary_smoothing_width,
        ),
    };

    V2FinalState {
        nx,
        ny,
        dx: 1.0 / nx as f64,
        dy: 1.0 / ny as f64,
        s_field: s.data().to_vec(),
        vx,
        vy,
        strain_rate_invariant: vec![0.0; nx * ny],
        age_field,
        cratonic_factor,
        plate_id: Some(plates.plate_id.data().to_vec()),
        plate_type: Some(plate_type),
        boundary_flag: None,
    }
}

/// Phase 8e — load a JSON snapshot from disk and replace the bridge
/// state with `V2RunState::Imported`. Step 8.6 follow-up: also walks
/// the snapshot's parent dir (`output/seed<S>_<R>/`) and re-populates
/// every phase cache that has artefacts on disk, so the user can
/// switch into Isostasy / UpscaleFbm / Erosion / Hydrology views and
/// see the loaded results — not just the tectonics output.
#[allow(clippy::too_many_arguments)]
fn handle_v2_import(
    mut bridge: ResMut<V2SolverBridge>,
    mut viz: ResMut<V2VizState>,
    mut iso_cache: ResMut<crate::phases::isostasy::IsostasyCache>,
    mut fbm_cache: ResMut<crate::phases::upscale_fbm::FbmCache>,
    mut erosion_cache: ResMut<crate::phases::erosion::ErosionCache>,
    mut hydrology_cache: ResMut<crate::phases::hydrology::HydrologyCache>,
) {
    let Some(path) = viz.import_requested_path.take() else {
        return;
    };
    match crate::bridge::v2::V2RunSnapshot::load(&path) {
        Ok(snap) => {
            bridge.state = V2RunState::Imported {
                spec: snap.spec,
                elapsed: std::time::Duration::from_secs_f64(snap.elapsed_seconds.max(0.0)),
                exported_at: snap.exported_at,
                scalar_metrics: snap.scalar_metrics,
                final_state: Box::new(snap.final_state),
            };
            // Pre-populate phase caches from any sibling artefacts in
            // the snapshot's parent dir. Best-effort: a missing /
            // corrupt phase artefact does not abort the import — the
            // tectonics state is already loaded above.
            if let Some(parent) = path.parent() {
                load_phase_caches_from_dir(
                    parent,
                    &mut iso_cache,
                    &mut fbm_cache,
                    &mut erosion_cache,
                    &mut hydrology_cache,
                );
            }
            viz.last_import = Some(Ok(path));
        }
        Err(e) => {
            viz.last_import = Some(Err(format!("{}", e)));
        }
    }
}

/// Walk a `seed<S>_<R>/` export directory and populate every phase
/// cache whose artefacts we find on disk. Each phase is independent —
/// a partial export (e.g. tectonics + isostasy only) loads what it
/// can and leaves the rest in their current state.
fn load_phase_caches_from_dir(
    dir: &std::path::Path,
    iso_cache: &mut crate::phases::isostasy::IsostasyCache,
    fbm_cache: &mut crate::phases::upscale_fbm::FbmCache,
    erosion_cache: &mut crate::phases::erosion::ErosionCache,
    hydrology_cache: &mut crate::phases::hydrology::HydrologyCache,
) {
    use ymir_core::erosion::hydraulic::ErosionResult;
    use ymir_core::export::PipelineExport;
    use ymir_core::grid::GridF32;
    use ymir_core::tectonics::isostasy::IsostasyResult;
    use ymir_core::terrain::upscale::UpscaleResult;

    let mut export = match PipelineExport::load(dir) {
        Ok(e) => e,
        Err(e) => {
            warn!("[import] no metadata.json in {}: {}", dir.display(), e);
            return;
        }
    };

    if let Some(iso_meta) = export.metadata.isostasy.clone() {
        match export.load_altitude() {
            Ok(heightmap) => {
                iso_cache.result = Some(IsostasyResult {
                    heightmap,
                    sea_level_normalized: iso_meta.sea_level_normalized,
                    peak_altitude_m: iso_meta.peak_altitude_m,
                    max_depth_m: iso_meta.max_depth_m,
                    land_ratio: iso_meta.land_ratio,
                });
                iso_cache.last_status = Some(Ok("loaded from disk".into()));
                iso_cache.mark_dirty();
                info!("[import] isostasy result restored");
            }
            Err(e) => warn!("[import] failed to load 02_altitude: {}", e),
        }
    }

    if export.metadata.upscale.is_some() {
        let upscaled = export.load_upscaled().or_else(|primary_err| {
            // Older exports could mis-record `target_size` in
            // metadata.json (e.g. user nudged the slider after running
            // FBM but before saving). Recover by inferring the actual
            // grid dims from the raw file size.
            let raw_path = dir.join("03_upscaled.raw");
            load_square_f32_grid_autosize(&raw_path).map(|h| {
                warn!(
                    "[import] 03_upscaled metadata mismatch ({}); recovered {}×{} from raw size",
                    primary_err, h.width, h.height
                );
                // Patch the in-memory metadata so the dependent
                // load_eroded / load_sediment paths use the correct
                // dims for the rest of this import.
                if let Some(u) = export.metadata.upscale.as_mut() {
                    u.target_size = h.width.max(h.height);
                }
                h
            })
        });
        match upscaled {
            Ok(heightmap) => {
                let nx = heightmap.width;
                let ny = heightmap.height;
                // Slope is not persisted; downstream callers (Erosion)
                // only consume `heightmap` so a zero placeholder is
                // safe. Re-running UpscaleFbm regenerates the real
                // slope from scratch.
                let slope = GridF32::from_vec(nx, ny, vec![0.0; nx * ny]);
                fbm_cache.result = Some(UpscaleResult { heightmap, slope, sediment: None });
                fbm_cache.last_status = Some(Ok("loaded from disk".into()));
                fbm_cache.mark_dirty();
                info!("[import] upscale_fbm result restored");
            }
            Err(e) => warn!("[import] failed to load 03_upscaled: {}", e),
        }
    }

    if let Some(ero_meta) = export.metadata.erosion.clone() {
        let eroded = export.load_eroded().or_else(|primary_err| {
            let raw_path = dir.join("04_eroded.raw");
            load_square_f32_grid_autosize(&raw_path).map(|h| {
                warn!(
                    "[import] 04_eroded metadata mismatch ({}); recovered {}×{}",
                    primary_err, h.width, h.height
                );
                h
            })
        });
        let sediment = export.load_sediment().or_else(|primary_err| {
            let raw_path = dir.join("04_sediment.raw");
            load_square_f32_grid_autosize(&raw_path).map(|h| {
                warn!(
                    "[import] 04_sediment metadata mismatch ({}); recovered {}×{}",
                    primary_err, h.width, h.height
                );
                h
            })
        });
        match (eroded, sediment) {
            (Ok(heightmap), Ok(sediment)) => {
                erosion_cache.result = Some(ErosionResult {
                    heightmap,
                    sediment,
                    stats: ero_meta.stats.clone(),
                });
                erosion_cache.state = crate::phases::erosion::ErosionState::Completed;
                erosion_cache.preview_heightmap = None;
                erosion_cache.progress = None;
                erosion_cache.last_status = Some(Ok("loaded from disk".into()));
                erosion_cache.mark_dirty();
                info!("[import] erosion result restored");
            }
            (Err(e), _) | (_, Err(e)) => {
                warn!("[import] failed to load 04_eroded / 04_sediment: {}", e)
            }
        }
    }

    if export.metadata.flow.is_some() {
        match export.load_flow() {
            Ok((flow_result, _river_cfg)) => {
                hydrology_cache.flow = Some(flow_result);
                // Rivers JSON is optional in older exports.
                match export.load_rivers() {
                    Ok(rivers) => hydrology_cache.rivers = Some(rivers),
                    Err(e) => warn!("[import] no 05_rivers.json: {}", e),
                }
                info!("[import] hydrology flow restored");
            }
            Err(e) => warn!("[import] failed to load flow artefacts: {}", e),
        }
    }

    if export.metadata.lakes.is_some() {
        match export.load_lakes() {
            Ok(lakes) => {
                hydrology_cache.lakes = Some(lakes);
                info!("[import] hydrology lakes restored");
            }
            Err(e) => warn!("[import] failed to load 06_lakes: {}", e),
        }
    }

    if hydrology_cache.flow.is_some() && hydrology_cache.lakes.is_some() {
        let segs = hydrology_cache
            .rivers
            .as_ref()
            .map(|r| r.segments.len())
            .unwrap_or(0);
        let lk = hydrology_cache
            .lakes
            .as_ref()
            .map(|l| l.lakes.len())
            .unwrap_or(0);
        hydrology_cache.last_status = Some(Ok(format!(
            "loaded from disk · {} river segments · {} lakes",
            segs, lk
        )));
        hydrology_cache.mark_dirty();
    }
}

/// Read an `f32`-raw grid file with no dimension priors — infer a
/// square `(side, side)` shape from the byte length. Used as a
/// last-resort recovery path when metadata-derived dimensions don't
/// match the on-disk file (older exports with the `target_size`
/// post-FBM-run mis-recorded). Errors out for non-square byte counts
/// since the export format only writes square v2 grids today.
fn load_square_f32_grid_autosize(
    path: &std::path::Path,
) -> Result<ymir_core::grid::GridF32, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("read: {}", e))?;
    if bytes.len() % 4 != 0 {
        return Err(format!("file size {} not a multiple of 4", bytes.len()));
    }
    let n = bytes.len() / 4;
    if n == 0 {
        return Err("empty file".into());
    }
    let side = (n as f64).sqrt().round() as usize;
    if side == 0 || side * side != n {
        return Err(format!(
            "{} f32 elements is not a perfect square (autosize fallback only handles square grids)",
            n
        ));
    }
    let data: Vec<f32> = bytes
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();
    Ok(ymir_core::grid::GridF32::from_vec(side, side, data))
}

pub struct V2VisualizationPlugin;

impl Plugin for V2VisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_v2_sprite);
        app.add_systems(
            Update,
            (
                update_v2_preview,
                update_v2_texture,
                render_v2_overlay_gizmos,
                handle_v2_screenshot,
                handle_v2_export,
                handle_v2_import,
            ),
        );
    }
}
