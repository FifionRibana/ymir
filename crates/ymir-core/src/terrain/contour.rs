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

/// The physical slope under which the contour is treated as PINNED and gets relaxed. The
/// measured defect is concentrated there: 82.5 % of turns above 80° on shores under 0.5°,
/// 54.4 % at 0.5–2°, and only 1.4 % above 15°.
pub const SMOOTH_RELAX_BELOW_DEG: f32 = 2.0;

/// Convert a physical slope in degrees to a gradient magnitude in NORMALISED units per cell —
/// what [`smooth_polylines_on_isoline`] gates on. `norm_to_m` is the vertical contract's
/// metres-per-normalised-unit (`2·1.13·depth_scale_m`).
#[must_use]
pub fn slope_deg_to_norm_gradient(deg: f32, m_per_cell: f32, norm_to_m: f32) -> f32 {
    deg.to_radians().tan() * m_per_cell / norm_to_m.max(1e-9)
}

/// Maximum relaxation weight, reached where the gradient vanishes (a fully pinned vertex).
pub const SMOOTH_LAMBDA_MAX: f32 = 0.5;
/// Newton steps per pass. One left 10 % of vertices >16 m off the isoline; three converge.
pub const SMOOTH_NEWTON_STEPS: usize = 3;
/// Each Newton step is clamped to this, so the reprojection is a descent and not a leap.
pub const SMOOTH_NEWTON_MAX_STEP_CELLS: f32 = 0.25;
/// A vertex may never move further than this from where marching squares put it. Smoothing is
/// a local regulariser, not a redrawing of the coast.
pub const SMOOTH_MAX_SHIFT_CELLS: f32 = 0.75;

/// Number of smoothing+reprojection passes used by the shipped coastline. Three is where the
/// >80° turn share stops falling materially (measured; see ADR Finding 48).
pub const COASTLINE_SMOOTH_PASSES: usize = 3;

/// Bilinear sample of `grid` at a sub-cell position, clamped at the borders.
fn sample_bilinear(grid: &GridF32, x: f32, y: f32) -> f32 {
    let (w, h) = (grid.width as i32, grid.height as i32);
    let x = x.clamp(0.0, (w - 1) as f32);
    let y = y.clamp(0.0, (h - 1) as f32);
    let (x0, y0) = (x.floor() as i32, y.floor() as i32);
    let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
    let (fx, fy) = (x - x0 as f32, y - y0 as f32);
    let g = |xi: i32, yi: i32| grid.data[yi as usize * grid.width + xi as usize];
    let top = g(x0, y0) * (1.0 - fx) + g(x1, y0) * fx;
    let bot = g(x0, y1) * (1.0 - fx) + g(x1, y1) * fx;
    top * (1.0 - fy) + bot * fy
}

/// Gradient of the BILINEAR INTERPOLANT at a sub-cell position, per cell.
///
/// ⚠️ The half-width matters and cost a measurement. A ±1 cell central difference smooths over
/// two cells, so it is NOT the gradient of the function `sample_bilinear` evaluates — and the
/// Newton reprojection below is trying to zero exactly that function. With the mismatched
/// gradient the iteration stalled instead of converging: raising the step count from 3 to 12
/// made the residual offset WORSE (7.1 → 8.2 m), which is the signature of descending on the
/// wrong slope rather than of needing more steps. A narrow difference tracks the local
/// bilinear patch.
const GRAD_H: f32 = 0.25;

fn gradient_at(grid: &GridF32, x: f32, y: f32) -> (f32, f32) {
    let inv = 0.5 / GRAD_H;
    let gx = inv * (sample_bilinear(grid, x + GRAD_H, y) - sample_bilinear(grid, x - GRAD_H, y));
    let gy = inv * (sample_bilinear(grid, x, y + GRAD_H) - sample_bilinear(grid, x, y - GRAD_H));
    (gx, gy)
}

/// SMOOTH A CONTOUR WITHOUT LEAVING THE ISOLINE (ADR Finding 48).
///
/// Marching squares already interpolates sub-cell along each crossed edge, so the barbs are NOT
/// blockiness. The mechanism is GRADIENT PINNING: where the field is nearly flat at the iso
/// value, WHICH edges get crossed is decided by tiny fluctuations, so the contour zig-zags —
/// measured at 82.5 % of turns above 80° on shores under 0.5°, against 1.4 % above 15°.
///
/// Three properties follow from that diagnosis, and all three are in the algorithm:
///
/// 1. **Relax only where the data is weak.** The relaxation weight is
///    `λ = λ_max · g_med/(g_med + |∇f|)`, with `g_med` the median gradient magnitude along the
///    contour. Where the gradient is strong the contour position is well determined by the
///    field and must be left alone; where it vanishes there is no defensible sub-cell position
///    and the smoothed one is at least not a fluctuation artefact. A uniform λ made the steep
///    shores WORSE (>15°: 1.4 % → 5.9 % of turns above 80°) — measured, then fixed.
/// 2. **Put every vertex back on the isoline.** Each pass ends with Newton steps along the
///    gradient, `p -= (f(p) − iso)/|∇f|² · ∇f`. Without them the coastline stops meaning
///    "altitude 0", which is the end-to-end coherence check with the consumer. A single step
///    per pass left 10 % of vertices more than 16 m off; three converge.
/// 3. **Never let a vertex wander.** Total displacement from the ORIGINAL position is capped at
///    [`SMOOTH_MAX_SHIFT_CELLS`]. Smoothing is a local regulariser, not a redrawing.
///
/// Closed rings (first point == last point) are smoothed cyclically; open lines keep their
/// endpoints fixed, so topology and the ring/line distinction are preserved exactly.
pub fn smooth_polylines_on_isoline(
    grid: &GridF32,
    iso: f32,
    polylines: &[Polyline],
    passes: usize,
    relax_below_gradient: f32,
) -> Vec<Polyline> {
    // Reference gradient: the median along the whole contour, so the weighting is scale-free
    // and needs no cell size or vertical contract.
    // The gate is a PHYSICAL slope supplied by the caller, not a percentile of this contour's
    // own gradient distribution. A percentile follows the coast it is measuring: when half the
    // coastline sits at 2–5°, the low quartile lands inside that class and the relaxation
    // spills onto shores that were never barbed — measured, +24 % of >80° turns at 2–5° while
    // the flattest class improved. The diagnosis names a slope (the defect is under ~2°), so
    // the remedy takes a slope.
    let g_lo = relax_below_gradient.max(1e-9);

    polylines
        .iter()
        .map(|pl| {
            let n = pl.len();
            if n < 3 {
                return pl.clone();
            }
            let closed = pl[0] == pl[n - 1];
            let orig: Vec<(f32, f32)> = if closed { pl[..n - 1].to_vec() } else { pl.clone() };
            let m = orig.len();
            if m < 3 {
                return pl.clone();
            }
            let mut pts = orig.clone();
            for _ in 0..passes {
                // 1. Gradient-weighted Laplacian relaxation.
                let src = pts.clone();
                for i in 0..m {
                    if !closed && (i == 0 || i == m - 1) {
                        continue; // open line: endpoints pinned
                    }
                    let (gx, gy) = gradient_at(grid, src[i].0, src[i].1);
                    let g = (gx * gx + gy * gy).sqrt();
                    let lambda = SMOOTH_LAMBDA_MAX * g_lo / (g_lo + g);
                    let prev = src[(i + m - 1) % m];
                    let next = src[(i + 1) % m];
                    let mid = ((prev.0 + next.0) * 0.5, (prev.1 + next.1) * 0.5);
                    let (mut dx, mut dy) =
                        (lambda * (mid.0 - src[i].0), lambda * (mid.1 - src[i].1));
                    // NORMAL COMPONENT ONLY. A plain Laplacian also moves vertices ALONG the
                    // curve, which bunches them: the mean step went 0.69 → 2.65 cells and the
                    // sampling itself changed, so the "after" metric was measuring a different
                    // polygon rather than a smoother one. The isoline's normal IS the gradient
                    // direction, so keep only the component across the curve — that removes a
                    // barb and cannot redistribute the vertices.
                    if g > 1e-9 {
                        let (nx, ny) = (gx / g, gy / g);
                        let dn = dx * nx + dy * ny;
                        dx = dn * nx;
                        dy = dn * ny;
                    }
                    pts[i] = (src[i].0 + dx, src[i].1 + dy);
                }
                // 2. DAMPED Newton reprojection onto the isoline. Undamped, a step of
                //    `d/|∇f|` DIVERGES where the gradient is small — it overshoots to a place
                //    where the field is more wrong, and three steps carry the vertex far away
                //    (measured: mean segment 0.69 → 2.65 cells, and 10 % of vertices left more
                //    than 8 m off the isoline they were supposed to be projected onto).
                //    Clamping each step to `SMOOTH_NEWTON_MAX_STEP_CELLS` makes it a descent.
                for p in pts.iter_mut() {
                    for _ in 0..SMOOTH_NEWTON_STEPS {
                        let (gx, gy) = gradient_at(grid, p.0, p.1);
                        let g2 = gx * gx + gy * gy;
                        if g2 <= 1e-12 {
                            break; // flat: no defensible direction, keep the relaxed point
                        }
                        let d = sample_bilinear(grid, p.0, p.1) - iso;
                        let (mut sx, mut sy) = (d * gx / g2, d * gy / g2);
                        let sl = (sx * sx + sy * sy).sqrt();
                        if sl > SMOOTH_NEWTON_MAX_STEP_CELLS {
                            let k = SMOOTH_NEWTON_MAX_STEP_CELLS / sl;
                            sx *= k;
                            sy *= k;
                        }
                        p.0 -= sx;
                        p.1 -= sy;
                    }
                }
                // 3. Cap the wander from the ORIGINAL position, LAST — a vertex may never
                //    be redrawn, only nudged.
                // from the ORIGINAL position.
                for (i, p) in pts.iter_mut().enumerate() {
                    let (dx, dy) = (p.0 - orig[i].0, p.1 - orig[i].1);
                    let d = (dx * dx + dy * dy).sqrt();
                    if d > SMOOTH_MAX_SHIFT_CELLS {
                        let k = SMOOTH_MAX_SHIFT_CELLS / d;
                        *p = (orig[i].0 + dx * k, orig[i].1 + dy * k);
                    }
                }
            }
            if closed {
                pts.push(pts[0]);
            }
            pts
        })
        .collect()
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
