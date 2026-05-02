//! Step 8.6 Phase 8g — visual revalidation post-corrections.
//!
//! Re-runs the Phase 7 4-preset gallery (`active_medley`,
//! `convergence`, `divergence`, `subduction`) with every Step 8.6
//! Phase 8 correction live:
//!
//! - Phase 8a `InitMode::Uniform` (default) — flat per-plate-type
//!   instead of the legacy sinusoidal perturbation.
//! - Phase 8b Voronoï boundary + per-plate velocity arrows overlaid
//!   on every captured frame.
//! - Phase 8c metrics (no UI surface here, but the run goes through
//!   the dashboard-visible code path).
//! - Phase 8d exposed knobs are taken from the preset JSON.
//!
//! Per-step PNGs land in
//! `docs/reports/step8_6_phase8g_visuals/<preset>/` plus a
//! patchwork composite (5 fields × N capture steps) at the preset
//! root. Reviewer judgement (D7 of the issue) decides whether
//! Phase 8h sunset is authorised.
//!
//! Marked `#[ignore]`. Defaults: 32² × 100 steps × capture every 10.
//! Override via env vars:
//!
//! ```text
//! YMIR_PHASE8G_GRID=64 YMIR_PHASE8G_STEPS=200 YMIR_PHASE8G_INTERVAL=20 \
//!     cargo test --release -p ymir-viz --test v2_phase8g_visuals \
//!         --jobs 1 -- --ignored --nocapture
//! ```

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use image::{imageops, GenericImage, ImageBuffer, Rgba, RgbaImage};
use ymir_core::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};
use ymir_viz::bridge_v2::{
    presets, spawn_v2_thread, V2Command, V2Event, V2FinalState, V2RunSpec,
};
use ymir_viz::visualization::overlay::{draw_velocity_vectors, draw_voronoi_boundaries};
use ymir_viz::visualization::v2_viz::{field_to_rgba, V2Field};

const PRESETS: &[&str] = &["active_medley", "convergence", "divergence", "subduction"];
const FIELD_TAGS: &[(V2Field, &str)] = &[
    (V2Field::SThickness, "s"),
    (V2Field::Age, "age"),
    (V2Field::Cratonic, "cratonic"),
    (V2Field::StrainRate, "strain"),
    (V2Field::VelocityMagnitude, "vmag"),
];
const VELOCITY_ARROW_SCALE: f64 = 8.0;

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn save_field_with_overlays(
    state: &V2FinalState,
    field: V2Field,
    plate_id: &[u16],
    path: &Path,
) -> std::io::Result<()> {
    let (nx, ny, mut rgba) = field_to_rgba(state, field);
    draw_voronoi_boundaries(&mut rgba, nx, ny, plate_id);
    draw_velocity_vectors(
        &mut rgba,
        nx,
        ny,
        &state.vx,
        &state.vy,
        plate_id,
        VELOCITY_ARROW_SCALE,
    );
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    image::save_buffer(path, &rgba, nx as u32, ny as u32, image::ColorType::Rgba8)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
}

fn run_preset_capture(
    preset_name: &str,
    spec: V2RunSpec,
    interval: usize,
    out_root: &Path,
) -> Vec<u32> {
    println!(
        "[phase8g] running '{}' — {}² × {} steps, capture every {}",
        preset_name, spec.grid_nx, spec.steps, interval
    );

    let plate_id = {
        let cfg = VoronoiConfig {
            num_plates: spec.num_plates,
            continental_ratio: spec.continental_ratio,
        };
        generate_voronoi(spec.grid_nx, spec.grid_ny, &cfg, spec.seed)
            .plate_id
            .data()
            .to_vec()
    };

    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(256);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    let t0 = Instant::now();
    cmd_tx
        .send(V2Command::RunBaseline { spec: spec.clone() })
        .expect("send command");

    // Generous deadline: 25 s/step at 64² mantle-on × 200 steps × 1.5
    // safety margin = 7500 s. At 32² × 100 steps × 7 s × 1.5 ≈ 1050 s.
    let secs_per_step = match spec.grid_nx {
        32 => 7.0,
        64 => 25.0,
        n => 25.0 * (n as f64 / 64.0).powi(2),
    };
    let deadline = Instant::now()
        + Duration::from_secs_f64((spec.steps as f64 * secs_per_step * 1.5).max(120.0));

    let preset_dir = out_root.join(preset_name);
    let _ = std::fs::remove_dir_all(&preset_dir);
    std::fs::create_dir_all(&preset_dir).expect("create preset dir");

    let mut captured_steps: Vec<u32> = Vec::new();
    let mut completed = false;
    while Instant::now() < deadline {
        match evt_rx.recv_timeout(Duration::from_millis(500)) {
            Ok(V2Event::Started { .. }) => {}
            Ok(V2Event::Progress { step, total, peek_state }) => {
                let is_capture = step == 0
                    || step == total
                    || (interval > 0 && step % interval == 0);
                if is_capture {
                    for &(field, tag) in FIELD_TAGS {
                        let path = preset_dir.join(format!("step_{:04}_{}.png", step, tag));
                        save_field_with_overlays(&peek_state, field, &plate_id, &path)
                            .unwrap_or_else(|e| {
                                panic!("save {} failed: {}", path.display(), e)
                            });
                    }
                    captured_steps.push(step as u32);
                    println!(
                        "[phase8g]   '{}' captured step {}/{} ({:.1}s)",
                        preset_name,
                        step,
                        total,
                        t0.elapsed().as_secs_f64()
                    );
                }
            }
            Ok(V2Event::Completed { final_state, .. }) => {
                let final_dir = preset_dir.join("final");
                std::fs::create_dir_all(&final_dir).expect("create final dir");
                for &(field, tag) in FIELD_TAGS {
                    let path = final_dir.join(format!("{}.png", tag));
                    save_field_with_overlays(&final_state, field, &plate_id, &path)
                        .unwrap_or_else(|e| panic!("save final {} failed: {}", tag, e));
                }
                completed = true;
                break;
            }
            Ok(V2Event::Failed { error }) => panic!("'{}' failed: {}", preset_name, error),
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("channel disconnected"),
        }
    }

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");
    assert!(completed, "'{}' did not complete within deadline", preset_name);
    assert!(!captured_steps.is_empty(), "'{}' captured zero frames", preset_name);

    println!(
        "[phase8g] '{}' done — {} captured steps in {:.1}s",
        preset_name,
        captured_steps.len(),
        t0.elapsed().as_secs_f64()
    );

    captured_steps
}

// ── Patchwork composer (mirrors v2_phase7_patchwork) ──────────────

const SCALE: u32 = 4;
const GUTTER: u32 = 2;
const FIELDS_FOR_PATCHWORK: &[&str] = &["s", "age", "cratonic", "strain", "vmag"];

fn parse_step_index(filename: &str, field_tag: &str) -> Option<u32> {
    let suffix = format!("_{}.png", field_tag);
    let stripped = filename.strip_suffix(&suffix)?;
    let num_str = stripped.strip_prefix("step_")?;
    num_str.parse().ok()
}

fn collect_field_frames(preset_dir: &Path, field_tag: &str) -> Vec<(u32, PathBuf)> {
    let mut hits = Vec::new();
    let Ok(entries) = std::fs::read_dir(preset_dir) else {
        return hits;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let Some(name) = path.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        if let Some(step) = parse_step_index(name, field_tag) {
            hits.push((step, path));
        }
    }
    hits.sort_by_key(|(step, _)| *step);
    hits
}

fn grid_dims(n: u32) -> (u32, u32) {
    if n == 0 {
        return (1, 1);
    }
    let cols = (n as f64).sqrt().ceil() as u32;
    let rows = (n + cols - 1) / cols;
    (cols, rows)
}

fn build_field_patchwork(preset_dir: &Path, field_tag: &str, cell_src: u32) -> Option<RgbaImage> {
    let frames = collect_field_frames(preset_dir, field_tag);
    if frames.is_empty() {
        return None;
    }
    let cell = cell_src * SCALE;
    let n = frames.len() as u32;
    let (cols, rows) = grid_dims(n);
    let total_w = cell * cols + GUTTER * (cols + 1);
    let total_h = cell * rows + GUTTER * (rows + 1);
    let mut canvas: RgbaImage =
        ImageBuffer::from_fn(total_w, total_h, |_, _| Rgba([0, 0, 0, 255]));
    for (idx, (_step, path)) in frames.iter().enumerate() {
        let img = image::open(path)
            .unwrap_or_else(|e| panic!("open {} failed: {}", path.display(), e))
            .to_rgba8();
        let cell_img = imageops::resize(&img, cell, cell, imageops::FilterType::Nearest);
        let col = (idx as u32) % cols;
        let row = (idx as u32) / cols;
        let x = GUTTER + (cell + GUTTER) * col;
        let y = GUTTER + (cell + GUTTER) * row;
        canvas.copy_from(&cell_img, x, y).expect("copy_from");
    }
    Some(canvas)
}

fn stack_vertically(parts: &[RgbaImage]) -> RgbaImage {
    let max_w = parts.iter().map(|p| p.width()).max().unwrap_or(0);
    let total_h: u32 =
        parts.iter().map(|p| p.height()).sum::<u32>() + GUTTER * (parts.len() as u32 + 1);
    let mut canvas: RgbaImage =
        ImageBuffer::from_fn(max_w, total_h, |_, _| Rgba([0, 0, 0, 255]));
    let mut y = GUTTER;
    for p in parts {
        canvas.copy_from(p, 0, y).expect("copy_from");
        y += p.height() + GUTTER;
    }
    canvas
}

fn build_patchwork_for_preset(preset_dir: &Path, grid: u32) {
    let cell_src = grid;
    let mut parts = Vec::new();
    for &field in FIELDS_FOR_PATCHWORK {
        let Some(canvas) = build_field_patchwork(preset_dir, field, cell_src) else {
            println!("[phase8g] {}: no '{}' frames", preset_dir.display(), field);
            continue;
        };
        let out = preset_dir.join(format!("_{}_patchwork.png", field));
        canvas
            .save(&out)
            .unwrap_or_else(|e| panic!("save {} failed: {}", out.display(), e));
        println!("[phase8g] wrote {}", out.display());
        parts.push(canvas);
    }
    if parts.len() >= 2 {
        let combined = stack_vertically(&parts);
        let out = preset_dir.join("_all.png");
        combined
            .save(&out)
            .unwrap_or_else(|e| panic!("save {} failed: {}", out.display(), e));
        println!("[phase8g] wrote {}", out.display());
    }
}

#[test]
#[ignore]
fn v2_phase8g_visuals() {
    let grid = env_usize("YMIR_PHASE8G_GRID", 32);
    let steps = env_usize("YMIR_PHASE8G_STEPS", 100);
    let interval = env_usize("YMIR_PHASE8G_INTERVAL", 10).max(1);

    let out_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step8_6_phase8g_visuals");
    std::fs::create_dir_all(&out_root).expect("create root");

    println!(
        "[phase8g] grid={}² steps={} interval={} presets={:?} → {}",
        grid,
        steps,
        interval,
        PRESETS,
        out_root.display()
    );

    for &preset_name in PRESETS {
        let mut spec = presets::load(preset_name).expect("preset load");
        spec.grid_nx = grid;
        spec.grid_ny = grid;
        spec.steps = steps;
        // Keep dt comparable to the canonical 64² × 100 setup
        // (dt ≈ 0.06). At larger step counts the simulated time
        // grows proportionally so the per-step physics doesn't
        // change.
        spec.total_time_nondim = (steps as f64) * 0.06;
        spec.capture_endpoints = false;

        run_preset_capture(preset_name, spec, interval, &out_root);
        build_patchwork_for_preset(&out_root.join(preset_name), grid as u32);
    }

    println!("[phase8g] all presets done → {}", out_root.display());
}
