//! Marching-squares contour extraction over a [`GridF32`].
//!
//! Traces the iso-level `c` of a scalar field as a set of polylines. Used by
//! [`crate::export::vector`] to derive the `.ymir` coastline (sea-level isoline)
//! and cliff-edge (slope-threshold isoline) vector layers.
//!
//! # Coordinate space (matches the `.ymir` raster orientation invariant)
//!
//! Vertices are emitted in **CELL coordinates**: `x` is the column
//! (`0..width-1`, increasing EAST) and `y` is the row (`0..height-1`). This is
//! the *same* row-major indexing the container's rasters use — no re-orientation
//! happens anywhere in the exchange format, so a vector layer overlays its
//! raster exactly. The format DEFINES row 0 as the SOUTH edge (see
//! [`crate::export::container`]); this module honours that by never flipping.
//!
//! # Determinism
//!
//! Output geometry is a pure function of `(grid, iso)`:
//! - cells are scanned in a fixed order (row-major, `y` outer, `x` inner);
//! - edge crossings are computed from the two grid samples in a canonical
//!   endpoint order, so adjacent cells agree bit-for-bit on a shared edge;
//! - segments are joined in scan order (the hash map is used only for lookup,
//!   never iterated for ordering);
//! - the two ambiguous saddle cases (5, 10) are resolved by the cell-center
//!   average — a fixed rule, no randomness.

use std::collections::HashMap;

use crate::grid::GridF32;

/// A traced contour polyline in cell coordinates. A closed ring (island /
/// lakeshore) has its first and last vertex equal.
pub type Polyline = Vec<(f32, f32)>;

/// Exact key for a vertex, used to weld shared endpoints. Adjacent cells
/// compute a shared-edge crossing with identical inputs → identical `f32`
/// bits, so exact equality (not epsilon) is correct and deterministic.
#[inline]
fn key(p: (f32, f32)) -> (u32, u32) {
    (p.0.to_bits(), p.1.to_bits())
}

/// Trace the `iso` contour of `grid`, returning polylines in cell coordinates.
///
/// A cell corner counts as *inside* when its value is `>= iso`. Empty when the
/// grid is smaller than `2×2` (no cell to march).
pub fn marching_squares(grid: &GridF32, iso: f32) -> Vec<Polyline> {
    let w = grid.width;
    let h = grid.height;
    if w < 2 || h < 2 {
        return Vec::new();
    }

    // Crossing point on the edge between grid samples a=(ax,ay) and b=(bx,by).
    // Canonical (endpoint-order-independent when callers pass the lower-index
    // sample first) so neighbouring cells produce the SAME vertex on a shared
    // edge.
    let cross = |ax: usize, ay: usize, bx: usize, by: usize| -> (f32, f32) {
        let va = grid.data[ay * w + ax];
        let vb = grid.data[by * w + bx];
        let denom = vb - va;
        let t = if denom == 0.0 { 0.5 } else { (iso - va) / denom };
        (ax as f32 + t * (bx as f32 - ax as f32), ay as f32 + t * (by as f32 - ay as f32))
    };

    let mut segments: Vec<[(f32, f32); 2]> = Vec::new();

    for y in 0..h - 1 {
        for x in 0..w - 1 {
            let v0 = grid.data[y * w + x]; // TL (x,   y)
            let v1 = grid.data[y * w + x + 1]; // TR (x+1, y)
            let v2 = grid.data[(y + 1) * w + x + 1]; // BR (x+1, y+1)
            let v3 = grid.data[(y + 1) * w + x]; // BL (x,   y+1)

            let mut idx = 0u8;
            if v0 >= iso {
                idx |= 1;
            }
            if v1 >= iso {
                idx |= 2;
            }
            if v2 >= iso {
                idx |= 4;
            }
            if v3 >= iso {
                idx |= 8;
            }
            if idx == 0 || idx == 15 {
                continue; // fully inside or fully outside — no crossing
            }

            // Edge crossings (lazy — only the crossed edges are evaluated).
            // Callers pass the lower-index sample first (see `cross`).
            let top = || cross(x, y, x + 1, y); // v0–v1
            let right = || cross(x + 1, y, x + 1, y + 1); // v1–v2
            let bottom = || cross(x, y + 1, x + 1, y + 1); // v3–v2
            let left = || cross(x, y, x, y + 1); // v0–v3

            match idx {
                1 | 14 => segments.push([top(), left()]),
                2 | 13 => segments.push([top(), right()]),
                3 | 12 => segments.push([left(), right()]),
                4 | 11 => segments.push([right(), bottom()]),
                6 | 9 => segments.push([top(), bottom()]),
                7 | 8 => segments.push([left(), bottom()]),
                5 => {
                    // v0,v2 inside; v1,v3 outside — saddle.
                    if 0.25 * (v0 + v1 + v2 + v3) >= iso {
                        // center inside → v0..v2 connected: cut around v1 and v3.
                        segments.push([top(), right()]);
                        segments.push([left(), bottom()]);
                    } else {
                        segments.push([top(), left()]);
                        segments.push([right(), bottom()]);
                    }
                }
                10 => {
                    // v1,v3 inside; v0,v2 outside — saddle.
                    if 0.25 * (v0 + v1 + v2 + v3) >= iso {
                        // center inside → v1..v3 connected: cut around v0 and v2.
                        segments.push([top(), left()]);
                        segments.push([right(), bottom()]);
                    } else {
                        segments.push([top(), right()]);
                        segments.push([left(), bottom()]);
                    }
                }
                _ => unreachable!("marching-squares idx is 4-bit, 0/15 handled"),
            }
        }
    }

    join_segments(segments)
}

/// Weld segments sharing an endpoint into maximal polylines, deterministically.
fn join_segments(segments: Vec<[(f32, f32); 2]>) -> Vec<Polyline> {
    // point key -> segment indices touching it (insertion order = scan order).
    let mut adj: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for (i, s) in segments.iter().enumerate() {
        adj.entry(key(s[0])).or_default().push(i);
        adj.entry(key(s[1])).or_default().push(i);
    }

    let mut used = vec![false; segments.len()];
    let mut polylines = Vec::new();

    // Find an unused segment touching `p`; consume it; return its far endpoint.
    let next = |adj: &HashMap<(u32, u32), Vec<usize>>,
                used: &mut [bool],
                p: (f32, f32)|
     -> Option<(f32, f32)> {
        let k = key(p);
        for &i in adj.get(&k)? {
            if used[i] {
                continue;
            }
            let s = segments[i];
            if key(s[0]) == k {
                used[i] = true;
                return Some(s[1]);
            } else if key(s[1]) == k {
                used[i] = true;
                return Some(s[0]);
            }
        }
        None
    };

    for start in 0..segments.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut poly = vec![segments[start][0], segments[start][1]];
        // Extend forward from the tail…
        while let Some(p) = next(&adj, &mut used, *poly.last().unwrap()) {
            poly.push(p);
        }
        // …then backward from the head (open contours only; a closed ring is
        // already whole and its head is used up).
        while let Some(p) = next(&adj, &mut used, *poly.first().unwrap()) {
            poly.insert(0, p);
        }
        polylines.push(poly);
    }

    polylines
}

/// Per-cell terrain slope in DEGREES from a metric height field and the metric
/// cell size. `slope_deg = atan(|∇h| / cell_size_m).to_degrees()`, where `|∇h|`
/// is the central-difference gradient magnitude (metres per cell). Edge cells
/// use the grid's forward/backward difference (a known small border artifact).
pub fn slope_degree_field(height_m: &GridF32, cell_size_m: f32) -> GridF32 {
    let w = height_m.width;
    let h = height_m.height;
    let inv = if cell_size_m > 0.0 { 1.0 / cell_size_m } else { 0.0 };
    let mut out = vec![0.0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let (gx, gy) = height_m.gradient_at(x, y); // metres per cell
            let rise_run = (gx * gx + gy * gy).sqrt() * inv; // m / m
            out[y * w + x] = rise_run.atan().to_degrees();
        }
    }
    GridF32::from_vec(w, h, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A radial bump crossing the iso once yields a single CLOSED ring.
    #[test]
    fn radial_bump_is_a_closed_ring() {
        let (w, h) = (32usize, 32usize);
        let (cx, cy) = (15.5f32, 15.5f32);
        let mut data = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let d = ((x as f32 - cx).powi(2) + (y as f32 - cy).powi(2)).sqrt();
                // High in the middle, low at the edges; iso 0.5 is crossed once.
                data[y * w + x] = (1.0 - d / 10.0).clamp(0.0, 1.0);
            }
        }
        let grid = GridF32::from_vec(w, h, data);
        let rings = marching_squares(&grid, 0.5);
        assert_eq!(rings.len(), 1, "one bump → one contour");
        let ring = &rings[0];
        assert!(ring.len() > 4, "ring must have vertices");
        assert_eq!(ring.first(), ring.last(), "island contour must be a closed loop");
        // The land center lies within the ring's bounding box.
        let (mut xmin, mut xmax, mut ymin, mut ymax) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
        for &(x, y) in ring {
            xmin = xmin.min(x);
            xmax = xmax.max(x);
            ymin = ymin.min(y);
            ymax = ymax.max(y);
        }
        assert!(xmin < cx && cx < xmax && ymin < cy && cy < ymax, "ring encloses the land");
    }

    #[test]
    fn deterministic_across_runs() {
        let (w, h) = (20usize, 20usize);
        let mut data = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let d = ((x as f32 - 9.5).powi(2) + (y as f32 - 9.5).powi(2)).sqrt();
                data[y * w + x] = (1.0 - d / 7.0).clamp(0.0, 1.0);
            }
        }
        let grid = GridF32::from_vec(w, h, data);
        assert_eq!(marching_squares(&grid, 0.5), marching_squares(&grid, 0.5));
    }

    #[test]
    fn empty_when_no_crossing() {
        let grid = GridF32::new(8, 8, 1.0);
        assert!(marching_squares(&grid, 0.5).is_empty(), "all-inside → no contour");
    }
}
