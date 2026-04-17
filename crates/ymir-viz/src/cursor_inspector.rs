use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

use crate::camera::MainCamera;
use crate::state::{
    CursorWorldPos, ErosionCache, FlowCache, GenerationParamsUi, IsostasyCache, IsostasyParams,
    PipelinePhase, TectonicState, UpscaleCache, ViewMode,
};
use crate::tectonic_view::thickness_color;
use crate::terrain_view::hypsometric_color;

pub struct CursorInspectorPlugin;

impl Plugin for CursorInspectorPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(
            EguiPrimaryContextPass,
            (cursor_inspector_overlay, minimap_overlay, scale_bar_overlay),
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn cursor_inspector_overlay(
    mut contexts: EguiContexts,
    cursor_pos: Res<CursorWorldPos>,
    tectonic: Option<Res<TectonicState>>,
    view_mode: Res<State<ViewMode>>,
    current_phase: Res<State<PipelinePhase>>,
    isostasy_cache: Res<IsostasyCache>,
    isostasy_params: Res<IsostasyParams>,
    upscale_cache: Res<UpscaleCache>,
    erosion_cache: Res<ErosionCache>,
    flow_cache: Res<FlowCache>,
    gen_params: Res<GenerationParamsUi>,
    dynamic_plates: Res<crate::state::DynamicPlateIds>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    egui::Area::new(egui::Id::new("cursor_inspector"))
        .anchor(egui::Align2::LEFT_BOTTOM, [65.0, -32.0])
        .show(ctx, |ui| {
            egui::Frame::new()
                .fill(egui::Color32::from_black_alpha(180))
                .corner_radius(4.0)
                .inner_margin(6.0)
                .show(ui, |ui| {
                    if *view_mode.get() == ViewMode::Tectonics {
                        draw_tectonic_info(ui, &cursor_pos, tectonic.as_deref(), &dynamic_plates);
                    } else {
                        draw_pipeline_terrain_info(
                            ui,
                            &cursor_pos,
                            current_phase.get(),
                            &isostasy_cache,
                            &isostasy_params,
                            &upscale_cache,
                            &erosion_cache,
                            &flow_cache,
                            &gen_params,
                        );
                    }
                });
        });
}

#[allow(clippy::too_many_arguments)]
fn draw_pipeline_terrain_info(
    ui: &mut bevy_egui::egui::Ui,
    cursor_pos: &CursorWorldPos,
    phase: &PipelinePhase,
    isostasy_cache: &IsostasyCache,
    isostasy_params: &IsostasyParams,
    upscale_cache: &UpscaleCache,
    erosion_cache: &ErosionCache,
    flow_cache: &FlowCache,
    gen_params: &GenerationParamsUi,
) {
    // Select active heightmap based on pipeline phase
    let heightmap = match *phase {
        PipelinePhase::Tectonics | PipelinePhase::Isostasy => isostasy_cache.heightmap.as_ref(),
        PipelinePhase::UpscaleFbm => upscale_cache.heightmap.as_ref(),
        PipelinePhase::Erosion | PipelinePhase::Hydrology => {
            erosion_cache.heightmap.as_ref().or(upscale_cache.heightmap.as_ref())
        }
        _ => None,
    };

    let Some(hm) = heightmap else {
        ui.monospace("x: ---  y: ---  alt: ---  slope: ---");
        return;
    };

    let Some(world) = cursor_pos.pos else {
        ui.monospace("x: ---  y: ---  alt: ---  slope: ---");
        return;
    };

    let gx = world.x + hm.width as f32 / 2.0;
    let gy = -world.y + hm.height as f32 / 2.0;
    let ix = gx as i32;
    let iy = gy as i32;

    if ix < 0 || iy < 0 || (ix as usize) >= hm.width || (iy as usize) >= hm.height {
        ui.monospace("x: ---  y: ---  alt: ---  slope: ---");
        return;
    }

    let x = ix as usize;
    let y = iy as usize;
    let h = hm.data[y * hm.width + x];
    let sea = isostasy_cache.sea_level_normalized;
    let max_elev = isostasy_params.max_elevation_m;
    let max_depth = isostasy_params.max_depth_m;

    let alt_m =
        if h >= sea { (h - sea) / (1.0 - sea) * max_elev } else { -((sea - h) / sea * max_depth) };

    let (grad_x, grad_y) = hm.gradient_at(x, y);
    let mpp = gen_params.meters_per_pixel;
    let slope_deg = ((grad_x * grad_x + grad_y * grad_y).sqrt() * mpp).atan().to_degrees();

    let mut text = format!("x:{x}  y:{y}  alt:{alt_m:.0}m  slope:{slope_deg:.1}°");

    // Flow accumulation if available
    if let Some(ref result) = flow_cache.result {
        let fw = result.accumulation.width;
        let fh = result.accumulation.height;
        if x < fw && y < fh {
            let flow = result.accumulation.data[y * fw + x] as u32;
            let basin = result.basins[y * fw + x];
            text.push_str(&format!("  flow:{flow}  basin:{basin}"));
        }
    }

    ui.monospace(text);
}

fn draw_tectonic_info(
    ui: &mut bevy_egui::egui::Ui,
    cursor_pos: &CursorWorldPos,
    tectonic: Option<&TectonicState>,
    dynamic_plates: &crate::state::DynamicPlateIds,
) {
    use ymir_core::tectonics::boundaries::BoundaryType;
    use ymir_core::tectonics::plates::PlateType;

    let Some(tectonic) = tectonic else {
        ui.monospace("x: ---  plate: ---  thickness: ---");
        return;
    };
    let init = &tectonic.init;
    let nx = init.grid_width;
    let ny = init.grid_height;

    if let Some(world) = cursor_pos.pos {
        let gx = world.x + nx as f32 / 2.0;
        let gy = -world.y + ny as f32 / 2.0;
        let ix = gx as i32;
        let iy = gy as i32;

        if ix >= 0 && iy >= 0 && (ix as usize) < nx && (iy as usize) < ny {
            let x = ix as usize;
            let y = iy as usize;
            let plate_id = init.plate_ids[y * nx + x];
            let plate = &init.plates[plate_id];
            let thickness = init.thickness.data[y * nx + x];
            let ptype = match plate.plate_type {
                PlateType::Continental => "Cont",
                PlateType::Oceanic => "Oce",
            };
            let speed =
                (plate.velocity.0 * plate.velocity.0 + plate.velocity.1 * plate.velocity.1).sqrt();
            ui.monospace(format!(
                "x: {}  y: {}  plate: {} ({})  t: {:.2}  v: {:.2}",
                x, y, plate_id, ptype, thickness, speed
            ));

            // Show boundary type if available
            if let Some(ref bt) = dynamic_plates.boundary_types {
                let bt_nx = dynamic_plates.grid_width;
                let bt_ny = dynamic_plates.grid_height;
                if bt_nx == nx
                    && bt_ny == ny
                    && let Some(&btype) = bt.get(y * nx + x)
                {
                    let bt_name = match btype {
                        BoundaryType::None => "Interior",
                        BoundaryType::Subduction => "Subduction",
                        BoundaryType::OceanicSubduction => "Ocean. Subduction",
                        BoundaryType::ContinentalCollision => "Collision",
                        BoundaryType::Rift => "Rift",
                    };
                    ui.monospace(format!("boundary: {}", bt_name));
                }
            }
        } else {
            ui.monospace("x: ---  plate: ---  thickness: ---");
        }
    } else {
        ui.monospace("x: ---  plate: ---  thickness: ---");
    }
}

// ── Minimap ──────────────────────────────────────────────────────────────

#[derive(Resource)]
struct MinimapTexture {
    texture: egui::TextureHandle,
    grid_width: usize,
    grid_height: usize,
    /// The view mode the texture was built for — triggers a rebuild on mode change.
    mode: ViewMode,
    /// Generation counter copied from TectonicState — triggers a rebuild when plates
    /// are regenerated while already in Tectonics mode.
    tectonic_generation: u64,
}

#[allow(clippy::too_many_arguments)]
fn minimap_overlay(
    mut contexts: EguiContexts,
    tectonic: Option<Res<TectonicState>>,
    view_mode: Res<State<ViewMode>>,
    current_phase: Res<State<PipelinePhase>>,
    isostasy_cache: Res<IsostasyCache>,
    upscale_cache: Res<UpscaleCache>,
    erosion_cache: Res<ErosionCache>,
    flow_cache: Res<FlowCache>,
    minimap: Option<Res<MinimapTexture>>,
    mut commands: Commands,
    camera_q: Query<(&Transform, &Projection), With<MainCamera>>,
    windows: Query<&Window>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    let mode = *view_mode.get();
    let phase = *current_phase.get();

    // Select the active heightmap based on pipeline phase
    let active_heightmap = if mode == ViewMode::Tectonics {
        None // handled separately via tectonic state
    } else {
        match phase {
            PipelinePhase::Tectonics | PipelinePhase::Isostasy => isostasy_cache.heightmap.as_ref(),
            PipelinePhase::UpscaleFbm => upscale_cache.heightmap.as_ref(),
            PipelinePhase::Erosion | PipelinePhase::Hydrology => {
                erosion_cache.heightmap.as_ref().or(upscale_cache.heightmap.as_ref())
            }
            _ => None,
        }
    };

    // Determine the grid dimensions to draw the viewport rectangle correctly
    let (grid_w, grid_h) = if mode == ViewMode::Tectonics {
        let Some(ref tec) = tectonic else { return };
        (tec.init.grid_width, tec.init.grid_height)
    } else {
        let Some(ref hm) = active_heightmap else { return };
        (hm.width, hm.height)
    };

    // Rebuild when: first run, mode changed, grid size changed, or underlying data changed.
    let tectonic_gen = tectonic.as_ref().map(|t| t.generation).unwrap_or(0);
    let needs_rebuild = minimap.is_none()
        || minimap.as_ref().is_some_and(|m| {
            m.mode != mode
                || m.grid_width != grid_w
                || m.grid_height != grid_h
                || (mode == ViewMode::Tectonics && m.tectonic_generation != tectonic_gen)
        })
        || isostasy_cache.is_changed()
        || upscale_cache.is_changed()
        || erosion_cache.is_changed()
        || flow_cache.is_changed();

    let texture_ref;
    let minimap_texture: &MinimapTexture;

    if needs_rebuild {
        let mini_size = 128usize;
        let mut pixels = vec![egui::Color32::BLACK; mini_size * mini_size];

        if mode == ViewMode::Tectonics {
            let Some(ref tec) = tectonic else { return };
            let init = &tec.init;
            let nx = init.grid_width;
            let ny = init.grid_height;
            for my in 0..mini_size {
                for mx in 0..mini_size {
                    let sx = (mx as f32 / mini_size as f32 * (nx - 1) as f32) as usize;
                    let sy = (my as f32 / mini_size as f32 * (ny - 1) as f32) as usize;
                    let thickness = init.thickness.data[sy * nx + sx];
                    let [r, g, b] = thickness_color(thickness);
                    pixels[my * mini_size + mx] = egui::Color32::from_rgb(r, g, b);
                }
            }
        } else if let Some(ref hm) = active_heightmap {
            for my in 0..mini_size {
                for mx in 0..mini_size {
                    let sx = mx as f32 / mini_size as f32 * (hm.width - 1) as f32;
                    let sy = my as f32 / mini_size as f32 * (hm.height - 1) as f32;
                    let alt = hm.sample_bilinear(sx, sy);
                    let [r, g, b] = hypsometric_color(alt);
                    pixels[my * mini_size + mx] = egui::Color32::from_rgb(r, g, b);
                }
            }
        }

        let image = egui::ColorImage {
            size: [mini_size, mini_size],
            pixels,
            source_size: egui::Vec2::new(mini_size as f32, mini_size as f32),
        };
        let texture = ctx.load_texture("minimap", image, egui::TextureOptions::LINEAR);
        texture_ref = MinimapTexture {
            texture,
            grid_width: grid_w,
            grid_height: grid_h,
            mode,
            tectonic_generation: tectonic_gen,
        };
        commands.insert_resource(MinimapTexture {
            texture: texture_ref.texture.clone(),
            grid_width: grid_w,
            grid_height: grid_h,
            mode,
            tectonic_generation: tectonic_gen,
        });
        minimap_texture = &texture_ref;
    } else {
        let Some(ref m) = minimap else { return };
        minimap_texture = m;
    }

    let mini_size = 128.0;

    egui::Area::new(egui::Id::new("minimap"))
        .anchor(egui::Align2::RIGHT_BOTTOM, [-(260.0 + 10.0), -32.0])
        .show(ctx, |ui| {
            egui::Frame::new()
                .fill(egui::Color32::from_black_alpha(200))
                .corner_radius(4.0)
                .inner_margin(2.0)
                .show(ui, |ui| {
                    let (response, painter) =
                        ui.allocate_painter(egui::vec2(mini_size, mini_size), egui::Sense::hover());
                    let rect = response.rect;

                    // Draw minimap image
                    painter.image(
                        minimap_texture.texture.id(),
                        rect,
                        egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                        egui::Color32::WHITE,
                    );

                    // Draw camera viewport rectangle
                    if let Ok((cam_transform, projection)) = camera_q.single()
                        && let Ok(window) = windows.single()
                    {
                        let Projection::Orthographic(ortho) = projection else {
                            return;
                        };
                        let scale = ortho.scale;
                        let win_w = window.width() * scale;
                        let win_h = window.height() * scale;

                        let cam_x = cam_transform.translation.x;
                        let cam_y = cam_transform.translation.y;

                        let gw = grid_w as f32;
                        let gh = grid_h as f32;

                        // Camera world coords → minimap coords
                        let nx = (cam_x + gw / 2.0) / gw;
                        let ny = (-cam_y + gh / 2.0) / gh;
                        let nw = win_w / gw;
                        let nh = win_h / gh;

                        let vr = egui::Rect::from_center_size(
                            egui::pos2(rect.min.x + nx * mini_size, rect.min.y + ny * mini_size),
                            egui::vec2(nw * mini_size, nh * mini_size),
                        );

                        painter.rect_stroke(
                            vr,
                            0.0,
                            egui::Stroke::new(1.0, egui::Color32::from_white_alpha(180)),
                            egui::StrokeKind::Outside,
                        );
                    }
                });
        });
}

// ── Scale bar ────────────────────────────────────────────────────────────

fn scale_bar_overlay(
    mut contexts: EguiContexts,
    camera_q: Query<&Projection, With<MainCamera>>,
    gen_params: Res<GenerationParamsUi>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    let Ok(projection) = camera_q.single() else { return };
    let Projection::Orthographic(ortho) = projection else { return };

    let mpp = gen_params.meters_per_pixel;
    // World units per screen pixel at current zoom
    let world_per_screen_px = ortho.scale;
    let meters_per_screen_px = world_per_screen_px * mpp;

    // Target: ~100 screen pixels for the bar
    let bar_pixels = 100.0;
    let raw_meters = meters_per_screen_px * bar_pixels;

    // Round to a nice number
    let (nice_meters, label) = nice_distance(raw_meters);
    let actual_bar_px = nice_meters / meters_per_screen_px;

    egui::Area::new(egui::Id::new("scale_bar")).anchor(egui::Align2::LEFT_TOP, [65.0, 32.0]).show(
        ctx,
        |ui| {
            egui::Frame::new()
                .fill(egui::Color32::from_black_alpha(180))
                .corner_radius(4.0)
                .inner_margin(6.0)
                .show(ui, |ui| {
                    ui.monospace(&label);
                    let (response, painter) =
                        ui.allocate_painter(egui::vec2(actual_bar_px, 6.0), egui::Sense::hover());
                    let r = response.rect;
                    let color = egui::Color32::from_white_alpha(200);
                    // Horizontal line
                    painter.line_segment(
                        [egui::pos2(r.min.x, r.center().y), egui::pos2(r.max.x, r.center().y)],
                        egui::Stroke::new(2.0, color),
                    );
                    // End ticks
                    painter.line_segment(
                        [egui::pos2(r.min.x, r.min.y), egui::pos2(r.min.x, r.max.y)],
                        egui::Stroke::new(1.0, color),
                    );
                    painter.line_segment(
                        [egui::pos2(r.max.x, r.min.y), egui::pos2(r.max.x, r.max.y)],
                        egui::Stroke::new(1.0, color),
                    );
                });
        },
    );
}

fn nice_distance(meters: f32) -> (f32, String) {
    let nice_values = [
        (50.0, "50 m"),
        (100.0, "100 m"),
        (200.0, "200 m"),
        (500.0, "500 m"),
        (1000.0, "1 km"),
        (2000.0, "2 km"),
        (5000.0, "5 km"),
        (10000.0, "10 km"),
        (20000.0, "20 km"),
        (50000.0, "50 km"),
        (100000.0, "100 km"),
    ];

    let mut best = nice_values[0];
    for &nv in &nice_values {
        if nv.0 <= meters * 1.5 {
            best = nv;
        }
    }
    (best.0, best.1.to_string())
}
