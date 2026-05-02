//! Step 8.6 Phase 7 follow-up — patchwork composer.
//!
//! Pure post-processing. Reads the per-step PNGs produced by
//! `v2_phase7_step_diagnostic` under
//! `docs/reports/step8_6_phase7_gallery/diagnostic/<preset>/` and
//! composes:
//!
//! 1. **Per-field patchwork** (one PNG per field per preset) — every
//!    `step_NNNN_<field>.png` for that field tiled in an
//!    approximately-square grid, sorted by step number ascending,
//!    left-to-right, top-to-bottom. Output: `_<field>_patchwork.png`
//!    inside the preset's diagnostic directory.
//! 2. **All-fields combined** — the 5 per-field patchworks stacked
//!    vertically into a single image. Output: `_all.png` inside the
//!    preset's diagnostic directory.
//!
//! Each source frame is 32×32; cells are upsampled 4× via
//! nearest-neighbour to 128×128, separated by a 2-px black gutter.
//!
//! Usage:
//!
//! ```text
//! # All presets that have diagnostic frames
//! cargo test --release -p ymir-viz --test v2_phase7_patchwork \
//!     --jobs 1 -- --ignored --nocapture
//!
//! # Restrict to a subset
//! YMIR_PATCHWORK_PRESETS=convergence,divergence \
//!   cargo test --release -p ymir-viz --test v2_phase7_patchwork \
//!     --jobs 1 -- --ignored --nocapture
//! ```
//!
//! The patchwork PNGs land in the same gitignored `diagnostic/`
//! tree, so they are regenerated on demand and not committed.

use std::path::{Path, PathBuf};

use image::{imageops, GenericImage, ImageBuffer, Rgba, RgbaImage};

const SCALE: u32 = 4;
const CELL_SRC: u32 = 32;
const CELL: u32 = CELL_SRC * SCALE; // 128
const GUTTER: u32 = 2;
const FIELDS: &[&str] = &["s", "age", "cratonic", "strain", "vmag"];

fn upscale(img: &RgbaImage) -> RgbaImage {
    imageops::resize(img, CELL, CELL, imageops::FilterType::Nearest)
}

fn parse_step_index(filename: &str, field_tag: &str) -> Option<u32> {
    // Expected pattern: `step_NNNN_<field>.png`
    let suffix = format!("_{}.png", field_tag);
    let stripped = filename.strip_suffix(&suffix)?;
    let num_str = stripped.strip_prefix("step_")?;
    num_str.parse().ok()
}

/// Collect `(step_index, path)` for every `step_*_<field>.png` under
/// `preset_dir`, sorted by step ascending.
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

/// Choose `(cols, rows)` such that `cols * rows >= n` and the grid is
/// as square as possible (cols = ceil(sqrt(n))).
fn grid_dims(n: u32) -> (u32, u32) {
    if n == 0 {
        return (1, 1);
    }
    let cols = (n as f64).sqrt().ceil() as u32;
    let rows = (n + cols - 1) / cols;
    (cols, rows)
}

/// Compose a single-field patchwork. Returns `None` if no frames.
fn build_field_patchwork(preset_dir: &Path, field_tag: &str) -> Option<RgbaImage> {
    let frames = collect_field_frames(preset_dir, field_tag);
    if frames.is_empty() {
        return None;
    }
    let n = frames.len() as u32;
    let (cols, rows) = grid_dims(n);
    let total_w = CELL * cols + GUTTER * (cols + 1);
    let total_h = CELL * rows + GUTTER * (rows + 1);
    let mut canvas: RgbaImage =
        ImageBuffer::from_fn(total_w, total_h, |_, _| Rgba([0, 0, 0, 255]));
    for (idx, (_step, path)) in frames.iter().enumerate() {
        let img = image::open(path)
            .unwrap_or_else(|e| panic!("open {} failed: {}", path.display(), e))
            .to_rgba8();
        let cell_img = upscale(&img);
        let col = (idx as u32) % cols;
        let row = (idx as u32) / cols;
        let x = GUTTER + (CELL + GUTTER) * col;
        let y = GUTTER + (CELL + GUTTER) * row;
        canvas
            .copy_from(&cell_img, x, y)
            .expect("copy_from cell into canvas");
    }
    Some(canvas)
}

/// Stack a list of patchworks vertically into a single image, padding
/// each row to the max width with black.
fn stack_vertically(parts: &[RgbaImage]) -> RgbaImage {
    let max_w = parts.iter().map(|p| p.width()).max().unwrap_or(0);
    let total_h: u32 = parts.iter().map(|p| p.height()).sum::<u32>()
        + GUTTER * (parts.len() as u32 + 1);
    let mut canvas: RgbaImage =
        ImageBuffer::from_fn(max_w, total_h, |_, _| Rgba([0, 0, 0, 255]));
    let mut y = GUTTER;
    for p in parts {
        canvas
            .copy_from(p, 0, y)
            .expect("copy_from row into stacked canvas");
        y += p.height() + GUTTER;
    }
    canvas
}

#[test]
#[ignore]
fn v2_phase7_patchwork() {
    let diag_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step8_6_phase7_gallery/diagnostic");

    let presets: Vec<String> = if let Ok(env) = std::env::var("YMIR_PATCHWORK_PRESETS") {
        env.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect()
    } else {
        // Default: every preset subdir under diagnostic/ that
        // contains at least one step PNG.
        let mut found = Vec::new();
        if let Ok(entries) = std::fs::read_dir(&diag_root) {
            for entry in entries.flatten() {
                if entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
                    if let Some(name) = entry.file_name().to_str() {
                        let preset_dir = diag_root.join(name);
                        let any_step = std::fs::read_dir(&preset_dir)
                            .ok()
                            .into_iter()
                            .flatten()
                            .flatten()
                            .any(|e| {
                                e.file_name()
                                    .to_str()
                                    .map(|n| n.starts_with("step_"))
                                    .unwrap_or(false)
                            });
                        if any_step {
                            found.push(name.to_string());
                        }
                    }
                }
            }
        }
        found.sort();
        found
    };

    if presets.is_empty() {
        panic!(
            "no diagnostic frames found under {} — run \
             v2_phase7_step_diagnostic first (per preset)",
            diag_root.display()
        );
    }

    for preset in &presets {
        let preset_dir = diag_root.join(preset);
        let mut field_patches = Vec::new();
        for &field in FIELDS {
            let Some(canvas) = build_field_patchwork(&preset_dir, field) else {
                println!("[patchwork] {}/{}: no frames, skipped", preset, field);
                continue;
            };
            let out = preset_dir.join(format!("_{}_patchwork.png", field));
            canvas
                .save(&out)
                .unwrap_or_else(|e| panic!("save {} failed: {}", out.display(), e));
            println!(
                "[patchwork] {}/{}: {}×{} px → {}",
                preset,
                field,
                canvas.width(),
                canvas.height(),
                out.display()
            );
            field_patches.push(canvas);
        }
        if field_patches.len() >= 2 {
            let combined = stack_vertically(&field_patches);
            let out = preset_dir.join("_all.png");
            combined
                .save(&out)
                .unwrap_or_else(|e| panic!("save {} failed: {}", out.display(), e));
            println!(
                "[patchwork] {}/_all: {}×{} px → {}",
                preset,
                combined.width(),
                combined.height(),
                out.display()
            );
        }
    }
}
