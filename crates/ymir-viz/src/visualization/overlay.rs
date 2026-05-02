//! Step 8.6 Phase 8b — overlays drawn on top of the v2 field texture.
//!
//! Two routines:
//!
//! - [`draw_voronoi_boundaries`] paints inter-plate boundaries in
//!   contrasting black using the `plate_id` raster.
//! - [`draw_velocity_vectors`] draws one arrow per plate, anchored at
//!   the plate's periodic-aware centroid, in bright yellow. Arrow
//!   length is proportional to mean per-plate velocity scaled by
//!   `arrow_scale_cells`; tiny velocities are skipped so quiescent
//!   regimes don't clutter the view.
//!
//! Both routines mutate an `&mut [u8]` RGBA8 buffer in place,
//! row-major, with the same Y-flip the v2 sprite uses (image row 0 =
//! top of the rendered sprite, mapped to grid row `ny - 1`). They are
//! pure (no Bevy / no Image) so the diagnostic and screenshot paths
//! can opt in just as easily as the live-render path.

const BOUNDARY_RGBA: [u8; 4] = [0, 0, 0, 255];
const VELOCITY_RGBA: [u8; 4] = [255, 240, 0, 255];
/// Arrows shorter than this many cells worth of magnitude are
/// suppressed — they would render as a single pixel and add noise.
const MIN_ARROW_CELLS: f64 = 1.0;

/// Paint pixels that lie on a plate-id discontinuity in
/// [`BOUNDARY_RGBA`]. Discontinuity = any 4-periodic neighbour
/// (forward `(i+1, j)` or `(i, j+1)`) carries a different plate id.
/// Forward-only sweep avoids drawing each boundary twice.
///
/// `rgba` is `nx · ny · 4` bytes, row-major; `plate_id` is `nx · ny`
/// `u16` row-major in **grid-space** (no Y-flip). The function applies
/// the Y-flip when it indexes into the image buffer.
pub fn draw_voronoi_boundaries(rgba: &mut [u8], nx: usize, ny: usize, plate_id: &[u16]) {
    debug_assert_eq!(rgba.len(), nx * ny * 4);
    debug_assert_eq!(plate_id.len(), nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let id = plate_id[j * nx + i];
            let ip = (i + 1) % nx;
            let jp = (j + 1) % ny;
            let id_right = plate_id[j * nx + ip];
            let id_down = plate_id[jp * nx + i];
            if id != id_right || id != id_down {
                paint_pixel(rgba, nx, ny, i, j, &BOUNDARY_RGBA);
            }
        }
    }
}

/// Draw a per-plate velocity arrow at each plate's periodic-aware
/// centroid. Arrow tail = centroid, head = centroid + `(v̄ₓ, v̄ᵧ) ·
/// arrow_scale_cells` in cell units. `arrow_scale_cells` is a tuning
/// knob rather than a physical constant: at 64² with `peak|v| ≈ 5`
/// (active_medley regime) a value of `8.0` produces head separation
/// ≈ 40 cells (pleasant scale on a 600px sprite).
pub fn draw_velocity_vectors(
    rgba: &mut [u8],
    nx: usize,
    ny: usize,
    vx: &[f64],
    vy: &[f64],
    plate_id: &[u16],
    arrow_scale_cells: f64,
) {
    debug_assert_eq!(rgba.len(), nx * ny * 4);
    debug_assert_eq!(plate_id.len(), nx * ny);
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);

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
            // Periodic centroid via circular mean — accumulate
            // (cos θ, sin θ) where θ = 2π · cell_centre / domain.
            let theta_x = 2.0 * PI * (i as f64 + 0.5) / nx_f;
            let theta_y = 2.0 * PI * (j as f64 + 0.5) / ny_f;
            entry.cos_x += theta_x.cos();
            entry.sin_x += theta_x.sin();
            entry.cos_y += theta_y.cos();
            entry.sin_y += theta_y.sin();
            entry.vx += vx[idx];
            entry.vy += vy[idx];
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
        let head_dx = mvx * arrow_scale_cells;
        let head_dy = mvy * arrow_scale_cells;
        let head_len = (head_dx * head_dx + head_dy * head_dy).sqrt();
        if head_len < MIN_ARROW_CELLS {
            continue;
        }

        let x1 = cx + head_dx;
        let y1 = cy + head_dy;
        draw_line_periodic(rgba, nx, ny, cx, cy, x1, y1, &VELOCITY_RGBA);

        // Two arrowhead barbs at ±150° from the shaft direction,
        // length min(head_len * 0.3, 4 cells). Ensures heads are
        // visible without dominating short arrows.
        let theta = mvy.atan2(mvx);
        let barb_len = (head_len * 0.3).min(4.0);
        for &phi_off in &[5.0 * PI / 6.0, -5.0 * PI / 6.0] {
            let phi = theta + phi_off;
            let xb = x1 + barb_len * phi.cos();
            let yb = y1 + barb_len * phi.sin();
            draw_line_periodic(rgba, nx, ny, x1, y1, xb, yb, &VELOCITY_RGBA);
        }
    }
}

/// Bresenham-light line rasteriser with periodic wrap. Walks a
/// fractional point from `(x0, y0)` to `(x1, y1)` in 1-pixel steps
/// (max-axis driven) and paints the cell containing each sample.
fn draw_line_periodic(
    rgba: &mut [u8],
    nx: usize,
    ny: usize,
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    color: &[u8; 4],
) {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let steps = dx.abs().max(dy.abs()).ceil() as usize;
    let steps = steps.max(1);
    for k in 0..=steps {
        let t = k as f64 / steps as f64;
        let xf = x0 + t * dx;
        let yf = y0 + t * dy;
        let xi = (xf.floor() as i32).rem_euclid(nx as i32) as usize;
        let yj = (yf.floor() as i32).rem_euclid(ny as i32) as usize;
        paint_pixel(rgba, nx, ny, xi, yj, color);
    }
}

/// Apply Y-flip and write a single pixel at grid coords `(i, j)`.
#[inline]
fn paint_pixel(rgba: &mut [u8], nx: usize, ny: usize, i: usize, j: usize, color: &[u8; 4]) {
    let img_row = ny - 1 - j;
    let p = (img_row * nx + i) * 4;
    rgba[p..p + 4].copy_from_slice(color);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two-plate vertical strip: id=0 for i<nx/2, id=1 otherwise.
    /// `draw_voronoi_boundaries` should paint a vertical seam at the
    /// switch column and the wrap column.
    #[test]
    fn voronoi_boundaries_marks_only_inter_plate_edges() {
        let nx = 8;
        let ny = 8;
        let mut plate_id = vec![0u16; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                if i >= nx / 2 {
                    plate_id[j * nx + i] = 1;
                }
            }
        }
        let mut rgba = vec![32u8; nx * ny * 4]; // start opaque grey
        for chunk in rgba.chunks_exact_mut(4) {
            chunk[3] = 255;
        }
        draw_voronoi_boundaries(&mut rgba, nx, ny, &plate_id);

        let pixel = |i: usize, j: usize| -> [u8; 4] {
            let img_row = ny - 1 - j;
            let p = (img_row * nx + i) * 4;
            [rgba[p], rgba[p + 1], rgba[p + 2], rgba[p + 3]]
        };

        // Boundary cells = those with i = nx/2 - 1 (right neighbour
        // on the other plate) or i = nx - 1 (wrap-right neighbour
        // on the other plate). Forward-only sweep paints the
        // upstream cell of each pair.
        for j in 0..ny {
            assert_eq!(pixel(nx / 2 - 1, j), BOUNDARY_RGBA, "left seam at i=3,j={}", j);
            assert_eq!(pixel(nx - 1, j), BOUNDARY_RGBA, "wrap seam at i=7,j={}", j);
        }
        // Interior cells stay grey.
        for j in 0..ny {
            for i in [0, 1, 2, 4, 5, 6] {
                let p = pixel(i, j);
                assert_ne!(p, BOUNDARY_RGBA, "non-boundary cell ({}, {}) was painted", i, j);
            }
        }
    }

    /// Constant velocity per plate: every cell in plate 0 has
    /// `(v̄ₓ, v̄ᵧ) = (3, 0)`, plate 1 has `(0, 0)`. The arrow for
    /// plate 0 should run roughly horizontally from the centroid;
    /// plate 1 should be skipped (zero velocity below MIN_ARROW_CELLS).
    #[test]
    fn velocity_vectors_skip_quiescent_plates() {
        let nx = 16;
        let ny = 16;
        let mut plate_id = vec![0u16; nx * ny];
        let mut vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                if i >= nx / 2 {
                    plate_id[j * nx + i] = 1;
                } else {
                    vx[j * nx + i] = 3.0;
                }
            }
        }
        let mut rgba = vec![0u8; nx * ny * 4];
        for chunk in rgba.chunks_exact_mut(4) {
            chunk[3] = 255;
        }
        draw_velocity_vectors(&mut rgba, nx, ny, &vx, &vy, &plate_id, 1.5);

        let any_yellow = rgba
            .chunks_exact(4)
            .any(|c| c == VELOCITY_RGBA);
        assert!(any_yellow, "no arrow drawn for plate 0 with v̄ₓ=3");

        // Plate 1 quadrant (right half), excluding the column where
        // plate 0's arrow may extend through (since vy = 0 the
        // arrow stays on plate 0's row, but periodic wrap could
        // bring it through). We test that the *vertical centre* of
        // the right half is untouched.
        let img_row = ny - 1 - (ny / 2);
        let p_idx = (img_row * nx + 3 * nx / 4) * 4;
        let centre = &rgba[p_idx..p_idx + 4];
        assert_ne!(centre, VELOCITY_RGBA, "plate 1 centre painted unexpectedly");
    }

    /// Periodic centroid handling: a plate that wraps around the
    /// horizontal edge (cells i ∈ [0, 1] ∪ [14, 15]) has its
    /// circular-mean centroid near the seam (i ≈ 0 or 16 mod 16),
    /// not in the middle of the domain (i ≈ 8). Verify the arrow
    /// is drawn near the wrap, not at the centre.
    #[test]
    fn velocity_vectors_use_periodic_centroid_for_wrapping_plate() {
        let nx = 16;
        let ny = 4;
        // Plate 0: cells i ∈ {0, 1, 14, 15} for every j (seam-spanning).
        // Plate 1: everywhere else.
        let mut plate_id = vec![1u16; nx * ny];
        for j in 0..ny {
            for &i in &[0_usize, 1, 14, 15] {
                plate_id[j * nx + i] = 0;
            }
        }
        // Plate 0 has a strong +x velocity to ensure the arrow draws.
        let mut vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for &i in &[0_usize, 1, 14, 15] {
                vx[j * nx + i] = 5.0;
            }
        }
        let mut rgba = vec![0u8; nx * ny * 4];
        for chunk in rgba.chunks_exact_mut(4) {
            chunk[3] = 255;
        }
        draw_velocity_vectors(&mut rgba, nx, ny, &vx, &vy, &plate_id, 1.0);

        let pixel = |i: usize, j: usize| -> [u8; 4] {
            let img_row = ny - 1 - j;
            let p = (img_row * nx + i) * 4;
            [rgba[p], rgba[p + 1], rgba[p + 2], rgba[p + 3]]
        };

        // The centroid should be near i ∈ {0, 15} (wrap), so some
        // pixel in the wrap region must be yellow. The arithmetic
        // mean of cell indices {0.5, 1.5, 14.5, 15.5} is 8.0 (i.e.
        // anti-centroid) — drawing at i=8 would be wrong.
        let any_in_wrap = (0..ny).any(|j| {
            (0..3).any(|i| pixel(i, j) == VELOCITY_RGBA)
                || (13..16).any(|i| pixel(i, j) == VELOCITY_RGBA)
        });
        assert!(
            any_in_wrap,
            "expected arrow near i ∈ {{0..2, 13..15}}; centroid mis-computed"
        );
        let any_at_anti_centroid = (0..ny).any(|j| pixel(8, j) == VELOCITY_RGBA);
        assert!(
            !any_at_anti_centroid,
            "arrow drawn at the arithmetic mean (i=8) — should use circular mean"
        );
    }
}
