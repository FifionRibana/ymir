// Step d1: the coarse-gallery presentation is RETAINED but no longer
// registered (the HD workspace supersedes it; a live-coarse layer is a d2
// decision). Allow dead code so the unregistered systems don't warn.
#![allow(dead_code)]

//! C1 Bevy plugin — sprite setup, per-frame render system, egui
//! control panel, cell inspector.
//!
//! The C1 sprite (`C1VizSprite`) is the sole map render after the v2 sunset
//! (the old engine switcher + the parallel `V2VizSprite` were removed).
//!
//! ## Closure toggles run-locked (Q-E4.2)
//!
//! egui checkboxes are `.enabled(matches!(state, Idle |
//! Completed | Failed))` — disabled during `Running`. UI hint
//! visible to the user when they hover during a run.
//!
//! ## Init params dirty flow (Q-E4.3)
//!
//! No separate Init button — mirrors v2 pattern (params → Run
//! submits a fresh `RunBaseline` command, worker re-inits from
//! spec). Params editable when not Running.
//!
//! ## Velocity overlay caveat (Q-E4.4)
//!
//! Velocity arrows use init-time per-plate velocities (Stage E2
//! pre-run clone). Tooltip on the toggle: "init-time only;
//! mid-run mutations from accretion / rifting splits NOT
//! reflected (Viz-0-bis follow-up)".

use bevy::asset::RenderAssetUsages;
use bevy::image::ImageSampler;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use bevy_egui::input::EguiWantsInput;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

use ymir_core::grid::GridF32;
use ymir_core::tectonics_v2::workflow::PhaseAParams;

use crate::bridge::c1::{
    C1CumulativeStats, C1RunKind, C1RunSpec, C1RunState, C1SolverBridge, HdParams, HdState,
};
use crate::camera::CursorWorldPos;
use crate::visualization::c1_viz::{C1Field, derive_altitude_field, field_to_rgba};
use crate::visualization::overlay::{draw_velocity_vectors, draw_voronoi_boundaries};

const C1_SPRITE_BASE_SIZE: f32 = 600.0;
const C1_SPRITE_Z: f32 = 12.0;

/// Which C1 pipeline the panel will submit on Run (Issue #139 Stage
/// E3). Default `Gallery` (the Issue #137 contract). Run-locked: set
/// at Run submit, displayed but not editable during a run.
#[derive(Resource, Clone, Copy, Default, PartialEq, Eq, Debug)]
pub enum ActivePipeline {
    #[default]
    Gallery,
    Workflow,
}

/// Hypsometric lens (Issue #139 Stage E3). Maps the **non-dim**
/// altitude (the verification value) to display **meters** — a
/// cosmetic curve that NEVER touches the model or the non-dim value
/// (W3 global). Land side is a free tunable lens; ocean side is
/// anchored to the Stein-Stein `depth_scale_m = 5000` so meters
/// invert the bathymetry back to physical depth.
#[derive(Resource, Clone, Copy, Debug)]
pub struct HypsometricScale {
    /// Land elevation scale (m) — tunable cosmetic lens.
    pub land_scale_m: f32,
    /// Land gamma (curve shape) — tunable.
    pub land_gamma: f32,
    /// Ocean depth scale (m) — anchored to Stein-Stein
    /// `depth_scale_m`; shown, not tuned.
    pub ocean_scale_m: f32,
    /// Ocean gamma — anchored at 1.0 (linear inverse of the
    /// non-dim bathymetry).
    pub ocean_gamma: f32,
}

impl Default for HypsometricScale {
    fn default() -> Self {
        Self { land_scale_m: 8000.0, land_gamma: 1.0, ocean_scale_m: 5000.0, ocean_gamma: 1.0 }
    }
}

/// Map a non-dim altitude to display meters via the hypsometric lens.
/// Land (`alt ≥ 0`) and ocean (`alt < 0`) use independent
/// scale/gamma. **Cosmetic only** — the non-dim altitude remains the
/// verification value (Issue #139 W3 global).
fn hypsometric_meters(alt_nondim: f32, scale: &HypsometricScale) -> f32 {
    if alt_nondim >= 0.0 {
        scale.land_scale_m * alt_nondim.powf(scale.land_gamma)
    } else {
        -(scale.ocean_scale_m * (-alt_nondim).powf(scale.ocean_gamma))
    }
}

#[derive(Component)]
pub struct C1VizSprite;

/// Per-cell hover readout (Issue #139 Stage E1). `altitude_nondim`
/// is the verification value (W3 global: non-dim = truth, meters =
/// cosmetic). Meters are added by the Stage E3 hypsometric lens and
/// rendered SECOND, below the non-dim row.
#[derive(Clone, Copy, Debug)]
pub struct HoverReadout {
    pub i: usize,
    pub j: usize,
    pub s: f64,
    pub age: f64,
    pub altitude_nondim: f32,
}

/// Run-locked UI state — editable spec the user is composing,
/// plus current field selection and overlay toggles.
#[derive(Resource)]
pub struct C1VizState {
    pub texture_handle: Option<Handle<Image>>,
    pub field: C1Field,
    /// User-editable spec; submitted on Run button.
    pub pending_spec: C1RunSpec,
    /// User-editable Phase A cadence for workflow mode (Issue #139
    /// Stage E3). Defaults to `PhaseAParams::default()` (the
    /// calibrated 5×20). Only the gallery vs workflow dispatch reads
    /// it; ignored in gallery mode.
    pub pending_phase_a: PhaseAParams,
    /// Toggle Voronoi boundary overlay (cardinality-agnostic;
    /// Track D dynamic plate_id handled natively by `draw_voronoi_boundaries`).
    pub show_voronoi_boundaries: bool,
    /// Toggle velocity arrows (init-time only; Q-E4.4).
    pub show_velocity_vectors: bool,
    /// Arrow length factor for the velocity overlay. C1 Phase 1.1
    /// per-plate velocities have magnitude ~0.01 (non-dim) — far
    /// smaller than v2's ~5 cells/step. To get visible arrows
    /// (≥ MIN_ARROW_CELLS = 1 cell long), `arrow_scale` must be
    /// ~100 minimum; default 500 gives ~5-cell arrows. Slider
    /// range `[50, 3000]` logarithmic in the panel.
    pub arrow_scale: f64,
    /// `(nx, ny, step, field_idx, overlay_bits)` for change detection.
    pub last_signature: Option<(usize, usize, usize, u8, u8)>,
    /// Lazily-derived Architecture C altitude field for the hover
    /// inspector, keyed by `snapshot.step`. Recomputed only when a
    /// newer snapshot arrives (W2: invalidate on snapshot change);
    /// cursor movement within the same step is a cache hit. Makes
    /// per-cell altitude available in ALL views (W4), not just the
    /// Altitude view.
    pub altitude_cache: Option<(usize, GridF32)>,
    /// Current per-cell hover readout, or `None` when the cursor is
    /// off the map / over an egui panel / a different engine is
    /// active. Written by `update_c1_hover`, read by
    /// `c1_hover_panel`.
    pub hover: Option<HoverReadout>,
}

impl Default for C1VizState {
    fn default() -> Self {
        Self {
            texture_handle: None,
            field: C1Field::default(),
            pending_spec: C1RunSpec::island_production(), // M1 #190 production island default
            // Issue #141: cap=0.92 is COUPLED with n_cycles≈12 — the
            // P95-cap equilibrium (bounded limit cycle ~30% emergent)
            // is reached by the worst-case band-entry cycle 9 + margin.
            // n_cycles=5 would cut mid-overshoot. PhaseAParams CORE
            // default stays 5 (shared with v2); this is the C1 viz
            // workflow default only.
            pending_phase_a: PhaseAParams { n_cycles: 12, ..PhaseAParams::default() },
            show_voronoi_boundaries: false,
            show_velocity_vectors: false,
            // C1-tuned default: 0.01 × 500 = 5 cells, comfortably
            // above MIN_ARROW_CELLS = 1. See arrow_scale docstring.
            arrow_scale: 500.0,
            last_signature: None,
            altitude_cache: None,
            hover: None,
        }
    }
}

/// Invert the C1 sprite transform to map a world-space cursor
/// position to a cell `(i, j)`. Returns `None` when the cursor is
/// outside the sprite bounds.
///
/// The sprite is centred at the world origin (`Transform::from_xyz(
/// 0, 0, C1_SPRITE_Z)`) with `custom_size = sprite_size(nx, ny)` and
/// a nearest sampler over the `nx × ny` texel image, so
/// `sprite_local == world`. Bevy sprite texel row 0 renders at the
/// top and world `+Y` points up, hence `v = (half.y − world.y)` for
/// the vertical axis (NOT `world.y + half.y`). Mirrors the
/// `sprite_size` formula used by `update_c1_texture` so the mapping
/// stays correct across grid resizes (W1).
fn world_to_cell(world: Vec2, nx: usize, ny: usize) -> Option<(usize, usize)> {
    if nx == 0 || ny == 0 {
        return None;
    }
    let size = sprite_size(nx, ny);
    let half = size / 2.0;
    let u = (world.x + half.x) / size.x;
    let v = (half.y - world.y) / size.y;
    if !(0.0..1.0).contains(&u) || !(0.0..1.0).contains(&v) {
        return None;
    }
    let i = ((u * nx as f32) as usize).min(nx - 1);
    let j = ((v * ny as f32) as usize).min(ny - 1);
    Some((i, j))
}

fn sprite_size(nx: usize, ny: usize) -> Vec2 {
    if nx == 0 || ny == 0 {
        return Vec2::splat(C1_SPRITE_BASE_SIZE);
    }
    let longer = nx.max(ny) as f32;
    Vec2::new(C1_SPRITE_BASE_SIZE * nx as f32 / longer, C1_SPRITE_BASE_SIZE * ny as f32 / longer)
}

fn setup_c1_sprite(mut commands: Commands, mut images: ResMut<Assets<Image>>) {
    let mut image = Image::new_fill(
        Extent3d { width: 1, height: 1, depth_or_array_layers: 1 },
        TextureDimension::D2,
        &[10, 10, 10, 255],
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    );
    image.sampler = ImageSampler::nearest();
    let handle = images.add(image);

    commands.insert_resource(C1VizState {
        texture_handle: Some(handle.clone()),
        ..C1VizState::default()
    });

    commands.spawn((
        Sprite { image: handle, custom_size: Some(sprite_size(1, 1)), ..default() },
        Transform::from_xyz(0.0, 0.0, C1_SPRITE_Z),
        Visibility::Inherited,
        C1VizSprite,
    ));
}

/// Per-frame render: reads cached snapshot via
/// `bridge.state.latest_snapshot()`, runs `field_to_rgba`, applies
/// optional overlays into the RGBA buffer, writes to the Bevy Image.
fn update_c1_texture(
    bridge: Res<C1SolverBridge>,
    mut viz: ResMut<C1VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<C1VizSprite>>,
) {
    let Some(snapshot) = bridge.state.latest_snapshot() else {
        return;
    };
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    let field_idx = viz.field as u8;
    let overlay_bits =
        (viz.show_voronoi_boundaries as u8) | ((viz.show_velocity_vectors as u8) << 1);
    let sig = (nx, ny, snapshot.step, field_idx, overlay_bits);
    if viz.last_signature == Some(sig) {
        return;
    }

    let Some(handle) = viz.texture_handle.clone() else {
        return;
    };
    let Some(image) = images.get_mut(&handle) else {
        return;
    };

    // Resize Image + sprite if grid changed.
    if image.width() as usize != nx || image.height() as usize != ny {
        *image = Image::new(
            Extent3d { width: nx as u32, height: ny as u32, depth_or_array_layers: 1 },
            TextureDimension::D2,
            vec![0u8; nx * ny * 4],
            TextureFormat::Rgba8UnormSrgb,
            RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
        );
        image.sampler = ImageSampler::nearest();
        let new_size = sprite_size(nx, ny);
        for mut s in sprite_q.iter_mut() {
            s.custom_size = Some(new_size);
        }
    }

    let mut rgba = field_to_rgba(snapshot, viz.field);

    if viz.show_voronoi_boundaries {
        draw_voronoi_boundaries(&mut rgba, nx, ny, &snapshot.plate_id);
    }
    if viz.show_velocity_vectors {
        let (vx, vy) = snapshot.expand_per_cell_velocities();
        draw_velocity_vectors(&mut rgba, nx, ny, &vx, &vy, &snapshot.plate_id, viz.arrow_scale);
    }

    if let Some(data) = image.data.as_mut() {
        data.copy_from_slice(&rgba);
    }
    viz.last_signature = Some(sig);
}

/// Hover-to-inspect (Issue #139 Stage E1). Maps the world-space
/// cursor to a cell, refreshes the lazily-derived altitude cache
/// when the snapshot step changes, and writes the per-cell readout
/// into `C1VizState.hover`. Suppressed (`hover = None`) when:
/// - the pointer is over an egui panel (not the map),
/// - the cursor is off the sprite, or
/// - there is no snapshot yet.
fn update_c1_hover(
    egui_input: Res<EguiWantsInput>,
    cursor: Res<CursorWorldPos>,
    bridge: Res<C1SolverBridge>,
    mut viz: ResMut<C1VizState>,
) {
    if egui_input.is_pointer_over_area() {
        viz.hover = None;
        return;
    }
    let Some(world) = cursor.pos else {
        viz.hover = None;
        return;
    };
    let Some(snapshot) = bridge.state.latest_snapshot() else {
        viz.hover = None;
        return;
    };
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    let Some((i, j)) = world_to_cell(world, nx, ny) else {
        viz.hover = None;
        return;
    };

    // Refresh the altitude cache if the snapshot step changed (W2).
    let step = snapshot.step;
    let stale = match &viz.altitude_cache {
        Some((cached_step, _)) => *cached_step != step,
        None => true,
    };
    if stale {
        let altitude = derive_altitude_field(snapshot);
        viz.altitude_cache = Some((step, altitude));
    }

    let altitude_nondim =
        viz.altitude_cache.as_ref().map(|(_, grid)| grid.get(i as i32, j as i32)).unwrap_or(0.0);
    let s = snapshot.s[j * nx + i];
    let age = snapshot.age[j * nx + i];

    viz.hover = Some(HoverReadout { i, j, s, age, altitude_nondim });
}

/// Fixed corner panel rendering the hover readout (Issue #139 Stage
/// E1, anchored bottom-right). Non-dim altitude is shown FIRST as
/// the verification value (W3 global). Meters (Stage E3 hypsometric
/// lens) will be added below, labelled "(after hypsometric curve)".
fn c1_hover_panel(mut contexts: EguiContexts, viz: Res<C1VizState>, hypso: Res<HypsometricScale>) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    egui::Window::new("C1 Cell Inspector")
        .anchor(egui::Align2::RIGHT_BOTTOM, [-16.0, -16.0])
        .resizable(false)
        .collapsible(false)
        .show(ctx, |ui| match viz.hover {
            Some(h) => {
                egui::Grid::new("c1_hover_grid").show(ui, |ui| {
                    ui.label("cell (i, j)");
                    ui.label(format!("({}, {})", h.i, h.j));
                    ui.end_row();
                    ui.label("S̃");
                    ui.label(format!("{:.4}", h.s));
                    ui.end_row();
                    ui.label("age");
                    ui.label(format!("{:.4}", h.age));
                    ui.end_row();
                    // Non-dim altitude FIRST — the verification value.
                    ui.label("altitude (non-dim)");
                    ui.label(format!("{:+.4}", h.altitude_nondim));
                    ui.end_row();
                    // Meters SECOND — cosmetic hypsometric lens.
                    ui.label("altitude (m)");
                    let meters = hypsometric_meters(h.altitude_nondim, &hypso);
                    ui.label(format!("{meters:+.0}"));
                    ui.end_row();
                });
                ui.weak("non-dim = verification; m = after hypsometric curve");
            }
            None => {
                ui.weak("Hover over the map…");
            }
        });
}

/// egui control panel — init params, closure toggles, Run/Cancel,
/// field switcher, overlay toggles, live stats.
fn c1_control_panel(
    mut contexts: EguiContexts,
    bridge: Res<C1SolverBridge>,
    mut viz: ResMut<C1VizState>,
    mut active_pipeline: ResMut<ActivePipeline>,
    mut hypso: ResMut<HypsometricScale>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };

    let is_running = matches!(bridge.state, C1RunState::Running { .. });

    egui::Window::new("C1 Engine Controls").default_pos([20.0, 80.0]).default_width(320.0).show(
        ctx,
        |ui| {
            // Pipeline toggle (run-locked, Issue #139 Stage E3).
            ui.horizontal(|ui| {
                ui.label("Pipeline:");
                ui.add_enabled_ui(!is_running, |ui| {
                    let mut pipe = *active_pipeline;
                    if ui
                        .selectable_value(&mut pipe, ActivePipeline::Gallery, "Gallery")
                        .on_hover_text(
                            "Standalone run_with_closures (Issue #137 contract). \
                             Static plate_type / coast; total = n_steps.",
                        )
                        .clicked()
                    {
                        *active_pipeline = pipe;
                    }
                    if ui
                        .selectable_value(&mut pipe, ActivePipeline::Workflow, "Workflow")
                        .on_hover_text(
                            "Calibrated Phase A loop (per-cycle post-tectonic: \
                             sea-level + macro-redistribution + reclassify). \
                             Coast migrates per cycle; total = n_cycles × k_cycle.",
                        )
                        .clicked()
                    {
                        *active_pipeline = pipe;
                    }
                });
            });
            if is_running {
                ui.weak("(pipeline run-locked during a run)");
            }
            ui.separator();

            // Init parameters (run-locked).
            ui.collapsing("Init parameters", |ui| {
                ui.add_enabled_ui(!is_running, |ui| {
                    egui::Grid::new("c1_init").show(ui, |ui| {
                        ui.label("grid_size");
                        ui.add(
                            egui::DragValue::new(&mut viz.pending_spec.grid_size)
                                .range(16..=256)
                                .speed(2),
                        );
                        ui.end_row();

                        ui.label("seed");
                        ui.add(egui::DragValue::new(&mut viz.pending_spec.seed));
                        ui.end_row();

                        ui.label("n_steps");
                        ui.add_enabled(
                            *active_pipeline == ActivePipeline::Gallery,
                            egui::DragValue::new(&mut viz.pending_spec.n_steps)
                                .range(10..=3000)
                                .speed(10),
                        )
                        .on_hover_text("Gallery only. Workflow uses n_cycles × k_cycle.");
                        ui.end_row();
                    });
                });
            });

            // Workflow cadence (workflow mode only, run-locked).
            if *active_pipeline == ActivePipeline::Workflow {
                ui.collapsing("Workflow cadence", |ui| {
                    ui.add_enabled_ui(!is_running, |ui| {
                        egui::Grid::new("c1_workflow_cadence").show(ui, |ui| {
                            ui.label("n_cycles");
                            ui.add(
                                egui::Slider::new(&mut viz.pending_phase_a.n_cycles, 1..=20)
                                    .text("(default 5)"),
                            );
                            ui.end_row();

                            ui.label("k_cycle");
                            ui.add(
                                egui::Slider::new(&mut viz.pending_phase_a.k_cycle, 1..=100)
                                    .text("(default 20)"),
                            );
                            ui.end_row();
                        });
                        let total = viz.pending_phase_a.n_cycles * viz.pending_phase_a.k_cycle;
                        ui.weak(format!(
                            "total tectonic steps = {} × {} = {}",
                            viz.pending_phase_a.n_cycles, viz.pending_phase_a.k_cycle, total
                        ));
                        // Calibration warning — the A1-c lesson.
                        let off_calibration =
                            viz.pending_phase_a.k_cycle != 20 || viz.pending_phase_a.n_cycles != 5;
                        if off_calibration {
                            ui.colored_label(
                                egui::Color32::from_rgb(220, 160, 40),
                                format!(
                                    "⚠ alpha = {:.3} is calibrated for k_cycle = 20, \
                                     n_cycles = 5. Raising the cadence without lowering \
                                     alpha over-erodes the continent (A1-c lesson).",
                                    viz.pending_phase_a.alpha
                                ),
                            );
                        }
                    });
                });
            }

            // Closure toggles (run-locked, Q-E4.2).
            ui.collapsing("Closures", |ui| {
                ui.add_enabled_ui(!is_running, |ui| {
                    ui.checkbox(
                        &mut viz.pending_spec.closures.davis_suppe.enabled,
                        "Davis-Suppe orogeny",
                    );
                    ui.checkbox(
                        &mut viz.pending_spec.closures.equilibrium_height.enabled,
                        "Equilibrium height",
                    );
                    ui.checkbox(
                        &mut viz.pending_spec.closures.erosion.enabled,
                        "Stream-power erosion",
                    );
                    ui.checkbox(
                        &mut viz.pending_spec.closures.oceanic_bathymetry.enabled,
                        "Oceanic bathymetry (Track A)",
                    );
                    ui.separator();
                    ui.label("Track D (Issue #132):");
                    ui.checkbox(&mut viz.pending_spec.closures.subduction.enabled, "Subduction");
                    ui.checkbox(&mut viz.pending_spec.closures.accretion.enabled, "Accretion");
                    ui.checkbox(&mut viz.pending_spec.closures.rifting.enabled, "Rifting");
                });
                if is_running {
                    ui.weak("(disabled during run; closures are run-locked)");
                }
            });

            // Run controls.
            ui.separator();
            ui.horizontal(|ui| {
                if ui.add_enabled(!is_running, egui::Button::new("▶ Run")).clicked() {
                    // Dispatch by pipeline (set at submit; run-locked).
                    match *active_pipeline {
                        ActivePipeline::Gallery => {
                            let _ = bridge.submit_run(viz.pending_spec.clone());
                        }
                        ActivePipeline::Workflow => {
                            let _ = bridge.submit_workflow(
                                viz.pending_spec.clone(),
                                viz.pending_phase_a.clone(),
                            );
                        }
                    }
                }
                if ui.add_enabled(is_running, egui::Button::new("■ Cancel")).clicked() {
                    bridge.request_cancel();
                }
            });

            // HD production (step b/5). Coarse validation surface: a button
            // that drives the full cached chain on the worker + a per-phase
            // status with the HIT/MISS regime. (The polished layer rendering
            // is the UI rewrite step d.)
            ui.separator();
            let hd_running = matches!(bridge.hd, HdState::Running { .. });
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(!hd_running, egui::Button::new("⛰ Generate HD"))
                    .on_hover_text(
                        "Run the full HD chain (eroded → climate → drainage → biomes) \
                         on the worker. First build is slow (erosion at 2048²); cached after.",
                    )
                    .clicked()
                {
                    let _ = bridge.submit_hd(viz.pending_spec.clone(), HdParams::default());
                }
                if hd_running {
                    ui.spinner();
                }
            });
            match &bridge.hd {
                HdState::Idle => {}
                HdState::Running { current, done, .. } => {
                    for r in done {
                        ui.weak(format!(
                            "✓ {} — {} ({:.1}s)",
                            r.phase.label(),
                            r.regime.label(),
                            r.elapsed.as_secs_f32()
                        ));
                    }
                    if let Some(p) = current {
                        ui.label(format!("⏳ {} …", p.label()));
                    }
                }
                HdState::Completed { result, total, done } => {
                    for r in done {
                        ui.weak(format!(
                            "✓ {} — {} ({:.1}s)",
                            r.phase.label(),
                            r.regime.label(),
                            r.elapsed.as_secs_f32()
                        ));
                    }
                    ui.colored_label(
                        egui::Color32::from_rgb(0x6E, 0xBE, 0x6E),
                        format!(
                            "✓ HD prêt : {}×{} · {} rivières · {} lacs · biomes ✓ ({:.1}s)",
                            result.width,
                            result.height,
                            result.drainage.rivers.segments.len(),
                            result.drainage.lakes.len(),
                            total.as_secs_f32(),
                        ),
                    );
                }
                HdState::Failed { error } => {
                    ui.colored_label(
                        egui::Color32::from_rgb(0xE1, 0x78, 0x46),
                        format!("HD : {error}"),
                    );
                }
            }

            // Field switcher.
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("Field:");
                egui::ComboBox::from_id_salt("c1_field").selected_text(viz.field.label()).show_ui(
                    ui,
                    |ui| {
                        for &f in C1Field::ALL {
                            ui.selectable_value(&mut viz.field, f, f.label());
                        }
                    },
                );
            });
            ui.weak(viz.field.legend_caption());

            // Overlays.
            ui.separator();
            ui.collapsing("Overlays", |ui| {
                ui.checkbox(&mut viz.show_voronoi_boundaries, "Voronoi boundaries");
                ui.horizontal(|ui| {
                    let resp = ui.checkbox(&mut viz.show_velocity_vectors, "Velocity arrows");
                    resp.on_hover_text(
                        "Init-time per-plate velocities only. Mid-run \
                         mutations from accretion / rifting splits NOT \
                         reflected (Viz-0-bis Option B follow-up).",
                    );
                });
                if viz.show_velocity_vectors {
                    ui.add(
                        egui::Slider::new(&mut viz.arrow_scale, 50.0..=3000.0)
                            .logarithmic(true)
                            .text("arrow scale (cells / unit velocity)"),
                    );
                }
            });

            // Hypsometric lens (Issue #139 Stage E3). Drives ONLY the
            // hover meters readout — never the non-dim value or the
            // model (W3 global). Editable any time (cosmetic).
            ui.separator();
            ui.collapsing("Hypsometric (hover meters)", |ui| {
                egui::Grid::new("c1_hypso").show(ui, |ui| {
                    ui.label("land scale (m)");
                    ui.add(
                        egui::Slider::new(&mut hypso.land_scale_m, 1000.0..=12000.0)
                            .text("tunable"),
                    );
                    ui.end_row();

                    ui.label("land gamma");
                    ui.add(egui::Slider::new(&mut hypso.land_gamma, 0.3..=3.0).text("tunable"));
                    ui.end_row();

                    ui.label("ocean scale (m)");
                    ui.label(format!("{:.0}", hypso.ocean_scale_m));
                    ui.end_row();

                    ui.label("ocean gamma");
                    ui.label(format!("{:.1}", hypso.ocean_gamma));
                    ui.end_row();
                });
                ui.weak("ocean anchored to Stein-Stein depth_scale_m; non-dim untouched");
            });

            // Live stats — Stage A bug fix: cumulative totals
            // accumulated across the StepCompleted stream, NOT
            // per-step (per-step reads 0 for rare events like
            // accretion merges and rifting splits even when the
            // run has accumulated several). For high-frequency
            // events (subduction, thinning) the display also
            // shows the last-step value alongside the cumulative
            // total to keep the panel "alive" during the run.
            ui.separator();
            ui.collapsing("Live stats", |ui| {
                match &bridge.state {
                    C1RunState::Idle => {
                        ui.label("Idle. Press ▶ Run to start.");
                    }
                    C1RunState::Running {
                        kind, step, total, latest_snapshot, cumulative, ..
                    } => {
                        // Workflow shows a per-cycle counter; gallery
                        // shows the flat step counter.
                        match kind {
                            C1RunKind::Workflow { phase_a } => {
                                let k = phase_a.k_cycle.max(1);
                                // step is the tectonic step (1-based
                                // after +1); cycle index = ceil(step/k).
                                let cycle = (*step).div_ceil(k).min(phase_a.n_cycles);
                                ui.label(format!(
                                    "Step {}/{}  ·  cycle {}/{}",
                                    step, total, cycle, phase_a.n_cycles
                                ));
                            }
                            C1RunKind::Gallery => {
                                ui.label(format!("Step {}/{}", step + 1, total));
                            }
                        }
                        if let Some(snap) = latest_snapshot {
                            render_live_stats_grid(
                                ui,
                                snap.live_plate_count,
                                snap.num_plates,
                                Some(&snap.stats),
                                cumulative,
                                false,
                            );
                        }
                    }
                    C1RunState::Completed { elapsed, final_snapshot, cumulative, .. } => {
                        ui.label(format!(
                            "Completed in {:.2?} ({} steps)",
                            elapsed,
                            final_snapshot.step + 1
                        ));
                        render_live_stats_grid(
                            ui,
                            final_snapshot.live_plate_count,
                            final_snapshot.num_plates,
                            None,
                            cumulative,
                            true,
                        );
                    }
                    C1RunState::Failed { error } => {
                        ui.colored_label(egui::Color32::RED, error);
                    }
                }
            });
        },
    );
}

/// Render the Track D live-stats grid. `per_step` is the current
/// step's stats (only present during `Running`); `cumulative` is
/// the run-cumulative totals accumulated by `poll_c1_events`. When
/// `is_completed` is true, the label says "final" instead of
/// "current"; otherwise the layout is identical.
fn render_live_stats_grid(
    ui: &mut egui::Ui,
    live_plates: usize,
    init_plates: usize,
    per_step: Option<&ymir_core::tectonics_c1::stats::C1StepStats>,
    cumulative: &C1CumulativeStats,
    is_completed: bool,
) {
    let live_label = if is_completed { "final live plates" } else { "live plates" };
    egui::Grid::new("c1_live_stats").show(ui, |ui| {
        ui.label(live_label);
        ui.label(format!("{live_plates} (init {init_plates})"));
        ui.end_row();

        // High-frequency events: cumulative + per-step.
        let sub_per = per_step.map(|s| s.subduction.cells_consumed).unwrap_or(0);
        ui.label("subduction cells");
        if per_step.is_some() {
            ui.label(format!("{}/step ({} total)", sub_per, cumulative.subduction_cells));
        } else {
            ui.label(format!("{} total", cumulative.subduction_cells));
        }
        ui.end_row();

        let thin_per = per_step.map(|s| s.rifting_thinning.cells_thinned).unwrap_or(0);
        ui.label("rifting cells thinned");
        if per_step.is_some() {
            ui.label(format!("{}/step ({} total)", thin_per, cumulative.thinning_cells));
        } else {
            ui.label(format!("{} total", cumulative.thinning_cells));
        }
        ui.end_row();

        // Rare events: cumulative only (per-step is almost always 0).
        ui.label("accretion merges");
        ui.label(format!("{} total", cumulative.accretion_merges));
        ui.end_row();

        ui.label("rifting splits");
        if cumulative.new_plate_ids.is_empty() {
            ui.label(format!("{} total", cumulative.rifting_splits));
        } else {
            ui.label(format!(
                "{} total (new pids: {:?})",
                cumulative.rifting_splits, cumulative.new_plate_ids
            ));
        }
        ui.end_row();
    });
}

pub struct C1VisualizationPlugin;

impl Plugin for C1VisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<ActivePipeline>()
            .init_resource::<HypsometricScale>()
            .add_systems(Startup, setup_c1_sprite)
            .add_systems(Update, (update_c1_texture, update_c1_hover))
            .add_systems(EguiPrimaryContextPass, (c1_control_panel, c1_hover_panel));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Default 64² grid → square sprite 600×600, half = 300.
    const N: usize = 64;

    #[test]
    fn world_to_cell_centre_is_mid_grid() {
        // World origin = sprite centre → centre cell.
        assert_eq!(world_to_cell(Vec2::new(0.0, 0.0), N, N), Some((32, 32)));
    }

    #[test]
    fn world_to_cell_y_axis_is_top_down() {
        // W1 Y-flip guard: an OFF-CENTRE, non-symmetric point. World
        // +Y is up; texel row 0 is the TOP. A point in the upper
        // half (world.y > 0) must map to a SMALL j (near the top),
        // not a large one. size = 600, half = 300.
        //   world.y = +150 → v = (300 - 150)/600 = 0.25 → j = 16.
        //   world.x = -150 → u = (-150 + 300)/600 = 0.25 → i = 16.
        // If the vertical term were (world.y + half)/size (the wrong
        // sign), v would be 0.75 → j = 48, failing this assert.
        assert_eq!(world_to_cell(Vec2::new(-150.0, 150.0), N, N), Some((16, 16)));
        // Symmetric lower-right point → large i and j.
        assert_eq!(world_to_cell(Vec2::new(150.0, -150.0), N, N), Some((48, 48)));
    }

    #[test]
    fn world_to_cell_outside_is_none() {
        // Beyond the +x edge (half = 300).
        assert_eq!(world_to_cell(Vec2::new(400.0, 0.0), N, N), None);
        // Beyond the -y edge.
        assert_eq!(world_to_cell(Vec2::new(0.0, -400.0), N, N), None);
    }

    #[test]
    fn hypsometric_meters_land_ocean_split() {
        let s = HypsometricScale::default(); // land 8000/γ1, ocean 5000/γ1
        // Land: linear at γ=1 → land_scale × alt.
        assert!((hypsometric_meters(1.0, &s) - 8000.0).abs() < 1e-3);
        assert!((hypsometric_meters(0.5, &s) - 4000.0).abs() < 1e-3);
        // Sea level maps to 0.
        assert!(hypsometric_meters(0.0, &s).abs() < 1e-3);
        // Ocean: negative, anchored to depth_scale_m = 5000.
        assert!((hypsometric_meters(-0.5, &s) + 2500.0).abs() < 1e-3);
    }

    #[test]
    fn hypsometric_meters_land_gamma_curves() {
        // γ = 2 compresses low elevations: alt 0.5 → 8000·0.25 = 2000.
        let s = HypsometricScale { land_gamma: 2.0, ..HypsometricScale::default() };
        assert!((hypsometric_meters(0.5, &s) - 2000.0).abs() < 1e-3);
        // Ocean side is independent of land gamma.
        assert!((hypsometric_meters(-0.5, &s) + 2500.0).abs() < 1e-3);
    }

    #[test]
    fn world_to_cell_respects_non_square_resize() {
        // W1 resize correctness: nx=4, ny=8 → longer=8,
        // size = (600·4/8, 600·8/8) = (300, 600), half = (150, 300).
        // World origin → centre cell (nx/2, ny/2) = (2, 4).
        assert_eq!(world_to_cell(Vec2::new(0.0, 0.0), 4, 8), Some((2, 4)));
        // A point at x = +75 (u = 0.75 → i = 3), y = +150
        // (v = (300-150)/600 = 0.25 → j = 2).
        assert_eq!(world_to_cell(Vec2::new(75.0, 150.0), 4, 8), Some((3, 2)));
    }
}
