//! Heightmap PNG saver with **dynamic** 16-bit remap and a
//! separate colour-bar strip.
//!
//! Step 0/1 used `GridF32::save_png_u16` which clamps `[0, 1]`. That
//! is fine for fields that truly live in that range; Step 2 adds
//! GPE spreading, which drives `S̃` across a much wider range, and
//! the clamp then hides the signal. The saver here rescales each
//! field to its own observed `[min, max]` and emits the bounds both
//! as the return value (for the markdown report) and as a tiny
//! gradient strip next to the map.
//!
//! The colour-bar is a plain vertical gradient — no tick labels, no
//! font rendering. The report markdown provides the numeric bounds
//! alongside the image so the visual can be read.

use std::path::{Path, PathBuf};

use crate::grid::GridF32;
use crate::tectonics_v2::field::Field2D;

/// Result of saving a heightmap: image path, value-range and
/// colour-bar path. The harness funnels these into the report
/// markdown as a small metadata block per snapshot.
#[derive(Clone, Debug)]
pub struct HeightmapMetadata {
    pub png_path: PathBuf,
    pub colorbar_path: PathBuf,
    pub min: f64,
    pub max: f64,
    pub mean: f64,
}

fn field_stats(field: &Field2D) -> (f64, f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0_f64;
    let n = field.data().len();
    for &v in field.data() {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        sum += v;
    }
    let mean = if n > 0 { sum / n as f64 } else { 0.0 };
    (min, max, mean)
}

/// Save `field` to `png_path` with dynamic `[min, max]` remap, and
/// write a small colour-bar strip next to it. The colour-bar path
/// is derived from the image path by inserting `_colorbar` before
/// the extension.
pub fn save_heightmap(
    field: &Field2D,
    png_path: &Path,
) -> Result<HeightmapMetadata, String> {
    let nx = field.nx();
    let ny = field.ny();
    let (min, max, mean) = field_stats(field);
    let span = (max - min).max(1e-30);

    // Main field remap to u16.
    let data: Vec<f32> = field
        .data()
        .iter()
        .map(|&v| {
            // Uniform-field fallback: if span is tiny, pin to mid-grey.
            if max - min < 1e-10 {
                0.5
            } else {
                (((v - min) / span).clamp(0.0, 1.0)) as f32
            }
        })
        .collect();
    let grid = GridF32::from_vec(nx, ny, data);
    if let Some(parent) = png_path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    grid.save_png_u16(png_path)?;

    // Colour-bar strip: 16 cols × 256 rows, linear gradient.
    let cb_w = 16;
    let cb_h = 256;
    let mut cb = Vec::with_capacity(cb_w * cb_h);
    for j in 0..cb_h {
        let t = 1.0 - (j as f32) / (cb_h as f32 - 1.0);
        for _ in 0..cb_w {
            cb.push(t);
        }
    }
    let cb_grid = GridF32::from_vec(cb_w, cb_h, cb);
    let colorbar_path = {
        let stem = png_path.file_stem().unwrap_or_default().to_string_lossy().into_owned();
        let ext = png_path.extension().unwrap_or_default().to_string_lossy().into_owned();
        let sib = png_path.parent().unwrap_or_else(|| Path::new("."));
        sib.join(format!("{}_colorbar.{}", stem, ext))
    };
    cb_grid.save_png_u16(&colorbar_path)?;

    Ok(HeightmapMetadata {
        png_path: png_path.to_path_buf(),
        colorbar_path,
        min,
        max,
        mean,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_field_does_not_divide_by_zero() {
        let f = Field2D::filled(8, 8, 3.14);
        let tmp = std::env::temp_dir().join("v2_heightmap_uniform.png");
        let md = save_heightmap(&f, &tmp).unwrap();
        assert_eq!(md.min, md.max);
        assert!((md.mean - 3.14).abs() < 1e-10);
    }

    #[test]
    fn nonuniform_field_reports_true_bounds() {
        let mut f = Field2D::new(4, 4);
        for j in 0..4 {
            for i in 0..4 {
                f.set(i, j, 0.2 + 1.8 * i as f64 / 3.0);
            }
        }
        let tmp = std::env::temp_dir().join("v2_heightmap_range.png");
        let md = save_heightmap(&f, &tmp).unwrap();
        assert!((md.min - 0.2).abs() < 1e-12);
        assert!((md.max - 2.0).abs() < 1e-12);
    }
}
