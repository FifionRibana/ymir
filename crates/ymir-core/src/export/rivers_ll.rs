//! `rivers_ll.json` — the rivers for Living Landz (ADR Findings 146–147).
//!
//! A NEW vector layer, alongside the unchanged `rivers.json`: the same drainage, regrouped into **main stems** (at each
//! confluence the upstream reach with the largest catchment continues the river; Hack's convention), each smoothed
//! into a polyline that stays in its valley floor, with the attributes a selection needs. No threshold is applied:
//! the selection is the consumer's.
//!
//! - **Smoothing**: Chaikin corner cutting on the D8 trace's cell centres, with fixed vertices: both endpoints, every
//!   confluence point and its receiver's trace neighbours, the tributary's second-to-last point, and every cell a
//!   spillway shares with a watercourse. A step between two fixed vertices stays the straight D8 step, so a tributary
//!   meets its receiver only at the confluence vertex (Finding 147-X). Each smoothed vertex stands at most
//!   `lateral_tol_cells` from the D8 trace (pulled back towards its nearest trace point by bisection otherwise). A
//!   spillway's runs over a watercourse's cells are replaced by that watercourse's own vertices (an exact overlap).
//!   Then a Ramer–Douglas–Peucker pass (fixed points kept) bounds the point count.
//! - **Per-vertex attributes** (Finding 147-V): catchment, accumulated discharge, slope, the construction's valley
//!   width W(A) and the bed width a·Q^b, read at the trace point each vertex comes from.
//! - **Coordinates**: continuous erosion-grid cells, the cell centre at `+0.5`, `y = 0` the SOUTH edge (the container's
//!   orientation invariant, `export/container.rs`).
//! - Living Landz owns the hex grid (Finding 147-V): the export carries no hex edge. [`HexGrid`] and
//!   [`river_hex_edges`] stay as the test reference (`docs/rivers_ll_hex_reference.json`).
//! - The format is documented in `docs/rivers_ll_format.md`; `format_version` is semver.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;
use serde::Serialize;

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, C1DrainageResult, LakeType, SegmentKind, own_end};
use crate::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

/// The format version this writer emits.
pub const RIVERS_LL_FORMAT_VERSION: &str = "0.3.0";

/// The bed width's coefficient `a` in `w = a·Q^b` (m per (m³/s)^b) — a DECISION (the author's), the default being
/// Ymir's own channel width (`rivers.json`'s `width_m`, drainage.rs `CHANNEL_WIDTH_A`). Written in the header.
pub const BED_WIDTH_A_DEFAULT: f32 = 5.0;
/// The bed width's exponent `b` — ANCHORED: downstream hydraulic geometry, Leopold & Maddock 1953 (USGS Professional
/// Paper 252).
pub const BED_WIDTH_B: f32 = 0.5;

/// The smoothing's and the annotations' parameters (declared in Findings 146–147).
#[derive(Clone, Copy, Debug)]
pub struct RiversLlParams {
    /// Chaikin iterations.
    pub chaikin_iters: u32,
    /// A smoothed vertex may stand at most this far from the D8 trace (cells).
    pub lateral_tol_cells: f32,
    /// Bisection steps when a vertex is pulled back.
    pub bisect_steps: u32,
    /// Ramer–Douglas–Peucker tolerance (cells) applied after the smoothing; 0 = none.
    pub rdp_eps_cells: f32,
    /// The bed width `a·Q^BED_WIDTH_B` (m).
    pub bed_width_a: f32,
    /// The construction's valley-floor width W(A) = `coef · A^exp` (m, A in km²).
    pub valley_width_coef_m: f32,
    pub valley_width_exp: f32,
    /// Two rivers are parallel when they run within this many cells of each other …
    pub parallel_k_cells: i32,
    /// … over more than this length (km), without sharing a confluence.
    pub parallel_min_km: f32,
    /// A spillway retraces when at least this share of its trace cells are watercourse cells.
    pub retrace_share: f32,
}

impl Default for RiversLlParams {
    fn default() -> Self {
        Self {
            chaikin_iters: 3,
            lateral_tol_cells: 0.5,
            bisect_steps: 8,
            rdp_eps_cells: 0.1,
            bed_width_a: BED_WIDTH_A_DEFAULT,
            // the témoin's W(10 km²) = 200 m, exponent 0.3 (`ValleyConstruction::f121`)
            valley_width_coef_m: 200.0 / 10f32.powf(0.3),
            valley_width_exp: 0.3,
            parallel_k_cells: 2,
            parallel_min_km: 2.0,
            retrace_share: 0.5,
        }
    }
}

/// How a river ends.
#[derive(Clone, Debug, Serialize, PartialEq)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum RiverEnd {
    Sea,
    Lake { lake_id: u32 },
    EndorheicLake { lake_id: u32 },
    /// Into a larger river (its id), at that river's vertex.
    Confluence { river_id: u32 },
    Terminal,
}

/// A tagged fall on the river (the gorge's format, ADR Finding 140). Empty while the gorge is paused.
#[derive(Clone, Debug, Serialize)]
pub struct FallLl {
    pub kind: String,
    pub x: f32,
    pub y: f32,
    pub height_m: f32,
}

/// ADR Finding 147-V — the attributes at each vertex of a river's polyline (columnar, parallel to `points`), read at
/// the D8 trace point the vertex comes from. Rounded to 4 significant digits.
#[derive(Clone, Debug, Default, Serialize)]
pub struct VertexAttrs {
    /// Drained area (km²): the flow accumulation.
    pub catchment_km2: Vec<f32>,
    /// Mean discharge (m³/s): the climate runoff accumulated along the drainage.
    pub discharge_m3s: Vec<f32>,
    /// The bed's slope along the river (m/m, positive downstream) over the neighbouring trace points.
    pub slope: Vec<f32>,
    /// The construction's valley-floor width W(A) (m).
    pub valley_width_m: Vec<f32>,
    /// The bed width a·Q^b (m), `a` and `b` in the header.
    pub bed_width_m: Vec<f32>,
}

/// One river (a main stem).
#[derive(Clone, Debug, Serialize)]
pub struct RiverLl {
    pub id: u32,
    pub kind: SegmentKind,
    /// The `rivers.json` segments it chains, upstream → downstream.
    pub segments: Vec<usize>,
    /// The smoothed polyline, upstream → downstream, in continuous erosion-grid cells.
    pub points: Vec<[f32; 2]>,
    /// Per-vertex attributes, parallel to `points`.
    pub vertex: VertexAttrs,
    pub length_km: f32,
    /// Geometric catchment at the mouth (km²).
    pub catchment_km2: f32,
    /// The largest Strahler order among its segments; `None` for a spillway (no hierarchy).
    pub strahler: Option<u8>,
    /// Mean discharge at the mouth (m³/s) and `rivers.json`'s channel width there (m).
    pub discharge_m3s: f32,
    pub width_m: f32,
    pub end: RiverEnd,
    /// Lakes touching its course (other than the one it ends in).
    pub lakes_crossed: Vec<u32>,
    pub falls: Vec<FallLl>,
    /// 1 = the longest river of its mouth basin.
    pub length_rank_in_basin: u32,
    /// The mouth basin's key: the last land cell (row-major index) on the D8 path from the river's end.
    pub mouth_basin: u32,
    /// ADR Finding 147-F — the higher-discharge river of this river's longest parallel pair (its own id when it is that
    /// river), `None` without a pair.
    pub parallel_of: Option<u32>,
    /// The lower-discharge river of at least one parallel pair.
    pub fusion_candidate: bool,
    /// A spillway whose trace mostly retraces watercourses (its overlapping runs are drawn on their vertices).
    pub retraces_spillway: bool,
}

#[derive(Serialize)]
struct WidthLaw {
    form: &'static str,
    a: f32,
    b: f32,
    a_status: &'static str,
    b_status: &'static str,
}

#[derive(Serialize)]
struct ValleyLaw {
    form: &'static str,
    coef_m: f32,
    exp: f32,
    status: &'static str,
}

#[derive(Serialize)]
struct RiversLlView<'a> {
    format_version: &'static str,
    coordinate_space: &'static str,
    km_per_cell: f32,
    discharge_source: &'static str,
    bed_width_law: WidthLaw,
    valley_width_law: ValleyLaw,
    rivers: &'a [RiverLl],
}

/// ADR Finding 146-H — Living Landz's hex grid (the author, 2026-10-05): flat-top hexes of `size_m` (centre to corner),
/// axial coordinates (q, r), hex (0, 0) centred on the map's bottom-left corner (`origin_m` = (0, 0) in metres from
/// that corner, `y` northward, as the container's `y = 0` = south). Since Finding 147 the grid lives in Living Landz;
/// this is the test reference.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct HexGrid {
    pub size_m: f32,
    pub orientation: &'static str,
    pub coordinates: &'static str,
    pub origin_m: [f32; 2],
}

impl HexGrid {
    /// The author's grid: 40 m flat-top, axial, (0, 0) at the bottom-left corner.
    pub fn living_landz() -> Self {
        Self { size_m: 40.0, orientation: "flat_top", coordinates: "axial", origin_m: [0.0, 0.0] }
    }

    /// The axial hex holding the point (metres from the map's bottom-left corner), by cube rounding.
    pub fn hex_of(&self, x_m: f32, y_m: f32) -> (i32, i32) {
        let (x, y) = (x_m - self.origin_m[0], y_m - self.origin_m[1]);
        let q = (2.0 / 3.0) * x / self.size_m;
        let r = (-x / 3.0 + 3f32.sqrt() / 3.0 * y) / self.size_m;
        let (cx, cz) = (q, r);
        let cy = -cx - cz;
        let (mut rx, ry, mut rz) = (cx.round(), cy.round(), cz.round());
        let (dx, dy, dz) = ((rx - cx).abs(), (ry - cy).abs(), (rz - cz).abs());
        if dx > dy && dx > dz {
            rx = -ry - rz;
        } else if dy <= dz {
            rz = -rx - ry;
        }
        (rx as i32, rz as i32)
    }
}

/// ADR Finding 146-H — every crossing of a hex edge by a river, columnar (one entry per crossing): the two hexes (axial,
/// upstream side first), the river, and at the crossing its catchment (km², the flow accumulation), discharge (m³/s,
/// the mouth's scaled by the catchment), local slope along the river (m/m) and the construction's floor width W(A) (m).
#[derive(Clone, Debug, Default, Serialize)]
pub struct HexEdges {
    pub q1: Vec<i32>,
    pub r1: Vec<i32>,
    pub q2: Vec<i32>,
    pub r2: Vec<i32>,
    pub river_id: Vec<u32>,
    pub catchment_km2: Vec<f32>,
    pub discharge_m3s: Vec<f32>,
    pub slope: Vec<f32>,
    pub width_m: Vec<f32>,
}

/// ADR Finding 146-H — the hex edges the rivers cross (the test reference since Finding 147; not exported). Each
/// river's polyline is walked every 0.1 hex size; a change of hex is a crossing. `width` is the construction's law
/// W(A) = `coef · A^exp` (the témoin: 200 m at 10 km², exponent 0.3).
#[allow(clippy::too_many_arguments)]
pub fn river_hex_edges(rivers: &[RiverLl], dr: &C1DrainageResult, field: &GridF32, ss: &SteinSteinParams, cell_km2: f32, hex: &HexGrid, width_coef_m: f32, width_exp: f32) -> HexEdges {
    let (w, h) = (field.width, field.height);
    let cell_m = cell_km2.sqrt() * 1000.0;
    let zm: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
    let acc = &dr.flow.accumulation.data;
    let mut e = HexEdges::default();
    let step = 0.1 * hex.size_m / cell_m; // in cells
    for r in rivers.iter() {
        let p = &r.points;
        if p.len() < 2 {
            continue;
        }
        // the polyline resampled every `step` cells
        let mut s: Vec<[f32; 2]> = Vec::new();
        for i in 0..p.len() - 1 {
            let (a, b) = (p[i], p[i + 1]);
            let l = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
            let m = (l / step).ceil().max(1.0) as usize;
            for j in 0..m {
                let tt = j as f32 / m as f32;
                s.push([a[0] + tt * (b[0] - a[0]), a[1] + tt * (b[1] - a[1])]);
            }
        }
        s.push(*p.last().expect("len >= 2"));
        let hx: Vec<(i32, i32)> = s.iter().map(|q| hex.hex_of(q[0] * cell_m, q[1] * cell_m)).collect();
        for i in 1..s.len() {
            let (a, b) = (hx[i - 1], hx[i]);
            if a == b {
                continue;
            }
            let c = [0.5 * (s[i - 1][0] + s[i][0]), 0.5 * (s[i - 1][1] + s[i][1])];
            let cell = |q: [f32; 2]| -> usize { (q[1].floor() as i64).rem_euclid(h as i64) as usize * w + (q[0].floor() as i64).rem_euclid(w as i64) as usize };
            let cc = cell(c);
            // the catchment: the largest accumulation in the crossing cell's 3×3, capped at the river's mouth
            let mut a_cells = 0f32;
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let q = ((cc / w) as i64 + dy).rem_euclid(h as i64) as usize * w + ((cc % w) as i64 + dx).rem_euclid(w as i64) as usize;
                    a_cells = a_cells.max(acc[q]);
                }
            }
            let area = (a_cells * cell_km2).min(r.catchment_km2).max(0.0);
            let q_m3s = if r.catchment_km2 > 0.0 { r.discharge_m3s * area / r.catchment_km2 } else { 0.0 };
            // the local slope along the river, over ±1 cell of arc (bilinear terrain)
            let back = s[i.saturating_sub((1.0 / step).ceil() as usize)];
            let fwd = s[(i + (1.0 / step).ceil() as usize).min(s.len() - 1)];
            let dist = ((fwd[0] - back[0]).powi(2) + (fwd[1] - back[1]).powi(2)).sqrt() * cell_m;
            let slope = if dist > 0.0 { (bilinear(&zm, w, h, back[0], back[1]) - bilinear(&zm, w, h, fwd[0], fwd[1])) / dist } else { 0.0 };
            e.q1.push(a.0);
            e.r1.push(a.1);
            e.q2.push(b.0);
            e.r2.push(b.1);
            e.river_id.push(r.id);
            e.catchment_km2.push(area);
            e.discharge_m3s.push(q_m3s);
            e.slope.push(slope);
            e.width_m.push(width_coef_m * area.max(0.0).powf(width_exp));
        }
    }
    e
}

/// The instruments of Findings 146-L / 147-T (the bench reads them; production ignores them). All on the FINAL
/// polyline unless said otherwise.
#[derive(Clone, Debug, Default)]
pub struct RiversLlStats {
    /// terrain − bed (m) over the samples: Chaikin before the constraint / the final polyline.
    pub above_bed_before: Vec<f32>,
    pub above_bed_after: Vec<f32>,
    /// The same on the raw D8 trace (the negative control).
    pub above_bed_raw: Vec<f32>,
    /// Lateral distance (cells) of each final sample from the D8 trace.
    pub lateral_cells: Vec<f32>,
    /// The same on Chaikin's output before the lateral constraint.
    pub lateral_before: Vec<f32>,
    /// Samples whose cell drains to another outlet than the river's own trace cell: final / raw.
    pub ridge_cross_smoothed: usize,
    pub ridge_cross_raw: usize,
    pub samples: usize,
    pub samples_raw: usize,
    /// Vertices pulled back by the constraint, of all smoothed vertices.
    pub pulled: usize,
    pub vertices: usize,
    /// The share of turns ≥ 45° on the raw trace and on the final polyline.
    pub stair_raw: f32,
    pub stair_smoothed: f32,
    /// Spillway trace points spliced onto a watercourse's vertices / spillway trace points in all.
    pub spliced_points: usize,
    pub spillway_points: usize,
    /// Reaches joining their receiver below its head (T3b), which end in a confluence instead of continuing it; of
    /// them, those larger than the receiver's reach above the junction.
    pub mid_joins: usize,
    pub mid_joins_larger: usize,
    /// The final polyline's ridge-crossing samples (river, point, whether the sample's cell is one of the river's own
    /// trace cells near it), for the diagnosis.
    pub ridge_samples: Vec<(u32, [f32; 2], bool)>,
    /// Parallel pairs found ((low id, high id), run length in cells).
    pub parallel_pairs: Vec<((u32, u32), usize)>,
}

/// Serialize the rivers as `rivers_ll.json` bytes; the header carries the width laws `params` used.
pub fn rivers_ll_json(rivers: &[RiverLl], km_per_cell: f32, params: &RiversLlParams) -> Vec<u8> {
    let view = RiversLlView {
        format_version: RIVERS_LL_FORMAT_VERSION,
        coordinate_space: "erosion_grid_cells_continuous",
        km_per_cell,
        discharge_source: "accumulated: the climate runoff (precipitation − evapotranspiration) accumulated down the D8 drainage, per trace point (C1DrainageResult::segment_discharge_profile_m3s)",
        bed_width_law: WidthLaw {
            form: "a * Q^b",
            a: params.bed_width_a,
            b: BED_WIDTH_B,
            a_status: "decision: the author's coefficient; the default is Ymir's own channel width (rivers.json width_m)",
            b_status: "anchored: downstream hydraulic geometry, Leopold & Maddock 1953 (USGS Professional Paper 252)",
        },
        valley_width_law: ValleyLaw {
            form: "coef_m * A_km2^exp",
            coef_m: params.valley_width_coef_m,
            exp: params.valley_width_exp,
            status: "the valley construction's floor width W(A) (ValleyConstruction, W(10 km²) = 200 m: a PROXY, Finding 120)",
        },
        rivers,
    };
    serde_json::to_vec(&view).expect("rivers_ll serialization is infallible")
}

/// Round to 4 significant digits (the per-vertex attributes; the serialized f32 then prints short).
fn sig4(x: f32) -> f32 {
    if x == 0.0 || !x.is_finite() {
        return if x.is_finite() { 0.0 } else { x };
    }
    let m = 10f64.powi(3 - (x.abs() as f64).log10().floor() as i32);
    ((x as f64 * m).round() / m) as f32
}

fn bilinear(z: &[f32], w: usize, h: usize, x: f32, y: f32) -> f32 {
    let (fx, fy) = (x - 0.5, y - 0.5);
    let (x0, y0) = (fx.floor(), fy.floor());
    let (tx, ty) = (fx - x0, fy - y0);
    let at = |xi: i64, yi: i64| -> f32 { z[yi.rem_euclid(h as i64) as usize * w + xi.rem_euclid(w as i64) as usize] };
    let (xi, yi) = (x0 as i64, y0 as i64);
    let a = at(xi, yi) * (1.0 - tx) + at(xi + 1, yi) * tx;
    let b = at(xi, yi + 1) * (1.0 - tx) + at(xi + 1, yi + 1) * tx;
    a * (1.0 - ty) + b * ty
}

/// Chaikin on an open polyline, `fixed[i]` vertices kept (and the endpoints). `origin` follows each point to its trace
/// index. A step between two fixed vertices stays straight (its new points are collinear with it).
fn chaikin(pts: &[[f32; 2]], fixed: &[bool], origin: &[usize]) -> (Vec<[f32; 2]>, Vec<bool>, Vec<usize>) {
    if pts.len() < 3 {
        return (pts.to_vec(), fixed.to_vec(), origin.to_vec());
    }
    let mut out = Vec::with_capacity(pts.len() * 2);
    let mut fx = Vec::with_capacity(pts.len() * 2);
    let mut org = Vec::with_capacity(pts.len() * 2);
    out.push(pts[0]);
    fx.push(true);
    org.push(origin[0]);
    for i in 0..pts.len() - 1 {
        let (a, b) = (pts[i], pts[i + 1]);
        let q = [0.75 * a[0] + 0.25 * b[0], 0.75 * a[1] + 0.25 * b[1]];
        let r = [0.25 * a[0] + 0.75 * b[0], 0.25 * a[1] + 0.75 * b[1]];
        if i > 0 && fixed[i] {
            out.push(a);
            fx.push(true);
            org.push(origin[i]);
        }
        out.push(q);
        fx.push(false);
        org.push(origin[i]);
        out.push(r);
        fx.push(false);
        org.push(origin[i + 1]);
    }
    out.push(pts[pts.len() - 1]);
    fx.push(true);
    org.push(origin[pts.len() - 1]);
    (out, fx, org)
}

/// Ramer–Douglas–Peucker keeping the fixed points; returns the kept indices.
fn rdp(pts: &[[f32; 2]], fixed: &[bool], eps: f32) -> Vec<usize> {
    if pts.len() <= 2 || eps <= 0.0 {
        return (0..pts.len()).collect();
    }
    let mut keep = vec![false; pts.len()];
    for (i, &f) in fixed.iter().enumerate() {
        keep[i] = f;
    }
    keep[0] = true;
    keep[pts.len() - 1] = true;
    let mut stack = Vec::new();
    // recurse between consecutive kept points
    let mut prev = 0usize;
    for i in 1..pts.len() {
        if keep[i] {
            stack.push((prev, i));
            prev = i;
        }
    }
    while let Some((s, e)) = stack.pop() {
        if e <= s + 1 {
            continue;
        }
        let (a, b) = (pts[s], pts[e]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let l2 = (dx * dx + dy * dy).max(1e-12);
        let mut best = (0f32, 0usize);
        for (i, p) in pts.iter().enumerate().take(e).skip(s + 1) {
            let t = (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0);
            let (px, py) = (a[0] + t * dx - p[0], a[1] + t * dy - p[1]);
            let d = (px * px + py * py).sqrt();
            if d > best.0 {
                best = (d, i);
            }
        }
        if best.0 > eps {
            keep[best.1] = true;
            stack.push((s, best.1));
            stack.push((best.1, e));
        }
    }
    (0..pts.len()).filter(|&i| keep[i]).collect()
}

fn poly_len(p: &[[f32; 2]]) -> f32 {
    p.windows(2).map(|s| ((s[1][0] - s[0][0]).powi(2) + (s[1][1] - s[0][1]).powi(2)).sqrt()).sum()
}

fn stair_share(p: &[[f32; 2]]) -> (usize, usize) {
    let mut turns = 0usize;
    let mut n = 0usize;
    for s in p.windows(3) {
        let (ax, ay) = (s[1][0] - s[0][0], s[1][1] - s[0][1]);
        let (bx, by) = (s[2][0] - s[1][0], s[2][1] - s[1][1]);
        let (la, lb) = ((ax * ax + ay * ay).sqrt(), (bx * bx + by * by).sqrt());
        if la < 1e-6 || lb < 1e-6 {
            continue;
        }
        n += 1;
        let c = ((ax * bx + ay * by) / (la * lb)).clamp(-1.0, 1.0);
        if c < 45f32.to_radians().cos() + 1e-4 {
            turns += 1;
        }
    }
    (turns, n)
}

/// The distance (cells) from `q` to the trace polyline `pts` near the trace index `org` (± 3 steps).
fn lateral(pts: &[[f32; 2]], q: [f32; 2], org: usize) -> f32 {
    lateral_span(pts, q, org, org)
}

/// The distance (cells) from `q` to the trace polyline `pts` between the trace indices `a ≤ b` (± 3 steps): a final
/// edge spans the origins of its two vertices, far apart after the simplification (Finding 147-T).
fn lateral_span(pts: &[[f32; 2]], q: [f32; 2], a: usize, b: usize) -> f32 {
    let n0 = pts.len();
    if n0 < 2 {
        return pts.first().map_or(0.0, |p| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt());
    }
    let lo = a.saturating_sub(3).min(n0 - 2);
    let hi = (b + 3).min(n0 - 1);
    let mut dmin = f32::INFINITY;
    for s in lo..hi {
        let (u, v) = (pts[s], pts[s + 1]);
        let (dx, dy) = (v[0] - u[0], v[1] - u[1]);
        let l2 = (dx * dx + dy * dy).max(1e-12);
        let tt = (((q[0] - u[0]) * dx + (q[1] - u[1]) * dy) / l2).clamp(0.0, 1.0);
        let d = ((u[0] + tt * dx - q[0]).powi(2) + (u[1] + tt * dy - q[1]).powi(2)).sqrt();
        dmin = dmin.min(d);
    }
    dmin
}

/// ADR Finding 146-P / 147-F — the parallel pairs: two rivers whose polylines (sampled every ≤ 1 cell) run within `k`
/// cells of each other (Chebyshev, on the samples' cells) over more than `l_cells` consecutive samples, the pairs in
/// `excluded` (a shared confluence, a spillway with the watercourse it retraces) left out. Returns each pair's longest
/// run (cells), sorted.
fn parallel_pairs(polys: &[&[[f32; 2]]], excluded: &HashSet<(u32, u32)>, k: i32, l_cells: usize) -> Vec<((u32, u32), usize)> {
    let samples: Vec<Vec<(i32, i32)>> = polys
        .par_iter()
        .map(|p| {
            let mut s = Vec::new();
            for e in 0..p.len().saturating_sub(1) {
                let (a, b) = (p[e], p[e + 1]);
                let l = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
                let m = l.ceil().max(1.0) as usize;
                for j in 0..m {
                    let t = j as f32 / m as f32;
                    s.push(((a[0] + t * (b[0] - a[0])).floor() as i32, (a[1] + t * (b[1] - a[1])).floor() as i32));
                }
            }
            if let Some(&l) = p.last() {
                s.push((l[0].floor() as i32, l[1].floor() as i32));
            }
            s
        })
        .collect();
    let mut grid: HashMap<(i32, i32), Vec<u32>> = HashMap::new();
    for (ri, s) in samples.iter().enumerate() {
        for &c in s {
            let v = grid.entry(c).or_default();
            if v.last() != Some(&(ri as u32)) {
                v.push(ri as u32);
            }
        }
    }
    let per_river: Vec<Vec<((u32, u32), usize)>> = samples
        .par_iter()
        .enumerate()
        .map(|(ri, s)| {
            let ri = ri as u32;
            let mut cur: HashMap<u32, usize> = HashMap::new();
            let mut best: HashMap<(u32, u32), usize> = HashMap::new();
            let mut near: HashSet<u32> = HashSet::new();
            for &(cx, cy) in s {
                near.clear();
                for dy in -k..=k {
                    for dx in -k..=k {
                        if let Some(l) = grid.get(&(cx + dx, cy + dy)) {
                            near.extend(l.iter().copied().filter(|&o| o != ri));
                        }
                    }
                }
                cur.retain(|o, _| near.contains(o));
                for &o in &near {
                    let run = cur.entry(o).or_insert(0);
                    *run += 1;
                    let key = (ri.min(o), ri.max(o));
                    if excluded.contains(&key) {
                        continue;
                    }
                    let b = best.entry(key).or_insert(0);
                    *b = (*b).max(*run);
                }
            }
            best.into_iter().filter(|e| e.1 > l_cells).collect()
        })
        .collect();
    let mut merged: HashMap<(u32, u32), usize> = HashMap::new();
    for v in per_river {
        for (key, len) in v {
            let b = merged.entry(key).or_insert(0);
            *b = (*b).max(len);
        }
    }
    let mut out: Vec<((u32, u32), usize)> = merged.into_iter().collect();
    out.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    out
}

/// A river's D8 trace: cell centres, bed (m), cells, and per point the accumulated discharge (m³/s) and the drained
/// area (km²).
struct Trace {
    pts: Vec<[f32; 2]>,
    bed: Vec<f32>,
    cells: Vec<usize>,
    q: Vec<f32>,
    area: Vec<f32>,
}

/// A spillway's run over a watercourse's trace: spillway trace indices `t0..=t1` are the watercourse `wr`'s trace
/// indices `k0..=k1` (`k1` may be below `k0`).
#[derive(Clone, Copy)]
struct Run {
    t0: usize,
    t1: usize,
    wr: usize,
    k0: usize,
    k1: usize,
}

/// A smoothed polyline with each vertex's trace origin and fixed flag.
struct Smoothed {
    pts: Vec<[f32; 2]>,
    org: Vec<usize>,
    fixed: Vec<bool>,
}

/// Build the main stems, smooth them in their valley and attach the attributes. `field` is the conditioned
/// (breached) field the drainage was computed on (normalised). With `with_stats`, the Findings 146-L / 147-T
/// instruments are filled (an O(n) outlet labelling more).
pub fn build_rivers_ll(dr: &C1DrainageResult, field: &GridF32, ss: &SteinSteinParams, cell_km2: f32, params: RiversLlParams, with_stats: bool) -> (Vec<RiverLl>, RiversLlStats) {
    let (w, h) = (field.width, field.height);
    let segs = &dr.rivers.segments;
    let ns = segs.len();
    let cell_m = cell_km2.sqrt() * 1000.0;
    let zm: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
    let acc = &dr.flow.accumulation.data;
    let recv = |c: usize| -> Option<usize> {
        let d = dr.flow.direction.get(c).copied().unwrap_or(DIR_NONE);
        if d == DIR_NONE {
            return None;
        }
        let x = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let y = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        Some(y * w + x)
    };
    let catch = |s: usize| dr.segment_catchment_cells.get(s).copied().unwrap_or(0.0);
    // the continuing upstream reach of each segment: the largest catchment (ties: the lowest index), among the reaches
    // that join it AT ITS HEAD. A reach joining it further down (the drainage's T3b confluence) cannot continue it: the
    // trace would jump back to the head and walk down again (Finding 147-X, amendment 2). It ends in a confluence on it.
    let cell_of = |p: (u32, u32)| p.1 as usize * w + p.0 as usize;
    let joins_head = |u: usize, d: usize| -> bool {
        match (segs[u].points.last(), segs[d].points.first()) {
            (Some(&ul), Some(&d0)) => ul == d0 || recv(cell_of(ul)) == Some(cell_of(d0)),
            _ => false,
        }
    };
    let (mut mid_joins, mut mid_joins_larger) = (0usize, 0usize);
    let mut cont_up = vec![usize::MAX; ns];
    for d in 0..ns {
        let mut best: Option<usize> = None;
        for &u in &segs[d].upstream {
            if u >= ns {
                continue;
            }
            if !joins_head(u, d) {
                mid_joins += 1;
                // larger than the receiver's own reach above the junction: Hack's convention is broken there
                let ul = segs[u].points.last().map(|&p| cell_of(p));
                let dc: Vec<usize> = segs[d].points.iter().map(|&p| cell_of(p)).collect();
                let j = ul.and_then(|c| dc.iter().position(|&x| x == c).or_else(|| recv(c).and_then(|r| dc.iter().position(|&x| x == r))));
                if let Some(j) = j
                    && j > 0
                    && catch(u) > acc[dc[j - 1]]
                {
                    mid_joins_larger += 1;
                }
                continue;
            }
            if best.is_none_or(|b| catch(u) > catch(b) || (catch(u) == catch(b) && u < b)) {
                best = Some(u);
            }
        }
        if let Some(b) = best {
            cont_up[d] = b;
        }
    }
    let mut seg_river = vec![u32::MAX; ns];
    let mut chains: Vec<Vec<usize>> = Vec::new();
    for s in 0..ns {
        if cont_up[s] != usize::MAX {
            continue; // not a head: it continues an upstream reach
        }
        let mut chain = vec![s];
        let mut c = s;
        while let Some(d) = segs[c].downstream {
            if d >= ns || cont_up[d] != c {
                break;
            }
            chain.push(d);
            c = d;
        }
        let rid = chains.len() as u32;
        for &x in &chain {
            seg_river[x] = rid;
        }
        chains.push(chain);
    }
    let nr = chains.len();
    let kind_of: Vec<SegmentKind> = chains.iter().map(|c| dr.segment_kind.get(*c.last().expect("non-empty")).copied().unwrap_or(SegmentKind::Watercourse)).collect();
    // the trace of each river (cell centres) with its bed, discharge and area, and its confluence end point
    let mut traces: Vec<Trace> = Vec::with_capacity(nr);
    for chain in &chains {
        let mut t = Trace { pts: Vec::new(), bed: Vec::new(), cells: Vec::new(), q: Vec::new(), area: Vec::new() };
        for &s in chain {
            let prof = dr.segment_profile_m.get(s);
            let qprof = dr.segment_discharge_profile_m3s.get(s);
            let q_seg = dr.segment_discharge_m3s.get(s).copied().unwrap_or(0.0);
            // the accumulation's scale to the segment's (possibly signified) catchment, at its own end (Finding 49)
            let own = segs[s].points.get(own_end(segs, s)).map_or(0.0, |&(x, y)| acc[y as usize * w + x as usize]);
            let k_s = if own > 0.0 { catch(s) / own } else { 1.0 };
            let a_cap = catch(s) * cell_km2;
            for (i, &(x, y)) in segs[s].points.iter().enumerate() {
                let c = y as usize * w + x as usize;
                let q = qprof.and_then(|p| p.get(i).copied()).unwrap_or(q_seg);
                let a = (acc[c] * k_s * cell_km2).min(a_cap).max(0.0);
                if t.cells.last() == Some(&c) {
                    // the junction cell shared with the previous reach: the downstream reach's values
                    *t.q.last_mut().expect("non-empty") = q;
                    *t.area.last_mut().expect("non-empty") = a;
                    continue;
                }
                t.cells.push(c);
                t.pts.push([x as f32 + 0.5, y as f32 + 0.5]);
                t.bed.push(prof.and_then(|p| p.get(i).copied()).unwrap_or(zm[c]));
                t.q.push(q);
                t.area.push(a);
            }
        }
        // a river ending in a confluence: its last point must be the receiving river's cell
        let last = *chain.last().expect("non-empty chain");
        if let Some(d) = segs[last].downstream
            && d < ns
            && let Some(&lc) = t.cells.last()
        {
            let dcells: Vec<usize> = segs[d].points.iter().map(|&(x, y)| y as usize * w + x as usize).collect();
            if !dcells.contains(&lc) {
                let r = recv(lc).filter(|r| dcells.contains(r)).or_else(|| dcells.first().copied());
                if let Some(r) = r {
                    let (q, a) = (*t.q.last().expect("non-empty"), *t.area.last().expect("non-empty"));
                    t.cells.push(r);
                    t.pts.push([(r % w) as f32 + 0.5, (r / w) as f32 + 0.5]);
                    t.bed.push(zm[r]);
                    t.q.push(q);
                    t.area.push(a);
                }
            }
        }
        traces.push(t);
    }
    // the fixed trace points (Finding 147-X): the endpoints; each confluence point on its receiver with the receiver's
    // trace neighbours; the tributary's second-to-last point
    let mut fixed: Vec<Vec<bool>> = traces.iter().map(|t| (0..t.pts.len()).map(|i| i == 0 || i + 1 == t.pts.len()).collect()).collect();
    for (ri, chain) in chains.iter().enumerate() {
        let last = *chain.last().expect("non-empty");
        let Some(d) = segs[last].downstream else { continue };
        if d >= ns || seg_river[d] == u32::MAX {
            continue;
        }
        let rj = seg_river[d] as usize;
        let Some(&lc) = traces[ri].cells.last() else { continue };
        let n0 = traces[ri].pts.len();
        if n0 >= 2 {
            fixed[ri][n0 - 2] = true;
        }
        if let Some(kk) = traces[rj].cells.iter().position(|&c| c == lc) {
            let nj = traces[rj].pts.len();
            for k2 in kk.saturating_sub(1)..=(kk + 1).min(nj - 1) {
                fixed[rj][k2] = true;
            }
        }
    }
    // the spillways' runs over watercourse traces (Finding 147-X): every shared cell is fixed on both rivers; a run of
    // ≥ 2 points is spliced onto the watercourse's vertices
    let spill_cells: HashSet<usize> = (0..nr).filter(|&ri| kind_of[ri] == SegmentKind::Spillway).flat_map(|ri| traces[ri].cells.iter().copied()).collect();
    let mut wc_at: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for ri in 0..nr {
        if kind_of[ri] == SegmentKind::Spillway {
            continue;
        }
        for (k2, c) in traces[ri].cells.iter().enumerate() {
            if spill_cells.contains(c) {
                wc_at.entry(*c).or_default().push((ri, k2));
            }
        }
    }
    let mut runs: Vec<Vec<Run>> = vec![Vec::new(); nr];
    let mut shared_pts = vec![0usize; nr];
    let mut retraced_wc: Vec<HashSet<usize>> = vec![HashSet::new(); nr];
    for ri in 0..nr {
        if kind_of[ri] != SegmentKind::Spillway {
            continue;
        }
        let cells = &traces[ri].cells;
        let n0 = cells.len();
        let mut t = 0usize;
        while t < n0 {
            let Some(list) = wc_at.get(&cells[t]) else {
                t += 1;
                continue;
            };
            // the longest consecutive match from t, over the candidates and both directions
            let mut best = (0usize, 0usize, 0usize, 1i64);
            for &(wr, k) in list {
                for dir in [1i64, -1] {
                    let wc = &traces[wr].cells;
                    let mut l = 1usize;
                    while t + l < n0 {
                        let kk = k as i64 + dir * l as i64;
                        if kk < 0 || kk as usize >= wc.len() || wc[kk as usize] != cells[t + l] {
                            break;
                        }
                        l += 1;
                    }
                    if l > best.0 || (l == best.0 && (wr, k) < (best.1, best.2)) {
                        best = (l, wr, k, dir);
                    }
                }
            }
            let (l, wr, k, dir) = best;
            let k_end = (k as i64 + dir * (l as i64 - 1)) as usize;
            fixed[ri][t] = true;
            fixed[ri][t + l - 1] = true;
            fixed[wr][k] = true;
            fixed[wr][k_end] = true;
            retraced_wc[ri].insert(wr);
            shared_pts[ri] += l;
            if l >= 2 {
                runs[ri].push(Run { t0: t, t1: t + l - 1, wr, k0: k, k1: k_end });
            }
            t += l;
        }
    }
    let mut stats = RiversLlStats { mid_joins, mid_joins_larger, ..Default::default() };
    // outlet labels for the ridge-crossing instrument (bench only)
    let outlet: Option<Vec<u32>> = with_stats.then(|| {
        let n = w * h;
        let mut lab = vec![u32::MAX; n];
        let mut path = Vec::new();
        for s in 0..n {
            if lab[s] != u32::MAX {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if lab[c] != u32::MAX {
                    v0 = lab[c];
                    break;
                }
                path.push(c);
                match recv(c) {
                    Some(r) if path.len() < 200_000 => c = r,
                    _ => {
                        v0 = c as u32;
                        break;
                    }
                }
            }
            for &q in &path {
                lab[q] = v0;
            }
        }
        lab
    });
    // the smoothing of one river (Chaikin + the lateral constraint), before any splice
    let smooth = |ri: usize, stats: &mut RiversLlStats| -> Smoothed {
        let t = &traces[ri];
        let n0 = t.pts.len();
        let (mut p, mut f, mut o) = (t.pts.clone(), fixed[ri].clone(), (0..n0).collect::<Vec<usize>>());
        for _ in 0..params.chaikin_iters {
            let r = chaikin(&p, &f, &o);
            p = r.0;
            f = r.1;
            o = r.2;
        }
        if with_stats {
            for i in 0..p.len() {
                stats.lateral_before.push(lateral(&t.pts, p[i], o[i]));
            }
        }
        // the nearest trace point of a smoothed point (its origin ± 2)
        let nearest = |q: [f32; 2], org: usize| -> usize {
            let lo = org.saturating_sub(2);
            let hi = (org + 2).min(n0 - 1);
            (lo..=hi)
                .min_by(|&a, &b| {
                    let da = (t.pts[a][0] - q[0]).powi(2) + (t.pts[a][1] - q[1]).powi(2);
                    let db = (t.pts[b][0] - q[0]).powi(2) + (t.pts[b][1] - q[1]).powi(2);
                    da.total_cmp(&db)
                })
                .unwrap_or(org.min(n0 - 1))
        };
        let tol = params.lateral_tol_cells;
        for i in 0..p.len() {
            stats.vertices += 1;
            if f[i] || lateral(&t.pts, p[i], o[i]) <= tol {
                continue;
            }
            stats.pulled += 1;
            let tp = t.pts[nearest(p[i], o[i])];
            let (mut lo, mut hi) = (0f32, 1f32); // fraction towards the trace point: 0 = the vertex, 1 = the trace point
            for _ in 0..params.bisect_steps {
                let mid = 0.5 * (lo + hi);
                let q = [p[i][0] + mid * (tp[0] - p[i][0]), p[i][1] + mid * (tp[1] - p[i][1])];
                if lateral(&t.pts, q, o[i]) > tol {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            p[i] = [p[i][0] + hi * (tp[0] - p[i][0]), p[i][1] + hi * (tp[1] - p[i][1])];
        }
        Smoothed { pts: p, org: o, fixed: f }
    };
    let simplify = |s: Smoothed| -> Smoothed {
        let keep = rdp(&s.pts, &s.fixed, params.rdp_eps_cells);
        Smoothed { pts: keep.iter().map(|&i| s.pts[i]).collect(), org: keep.iter().map(|&i| s.org[i]).collect(), fixed: keep.iter().map(|&i| s.fixed[i]).collect() }
    };
    let mut finals: Vec<Option<Smoothed>> = (0..nr).map(|_| None).collect();
    // the watercourses first …
    for ri in 0..nr {
        if kind_of[ri] == SegmentKind::Spillway {
            continue;
        }
        let s = smooth(ri, &mut stats);
        finals[ri] = Some(simplify(s));
    }
    // … then the spillways, their runs spliced onto the watercourses' final vertices
    for ri in 0..nr {
        if kind_of[ri] != SegmentKind::Spillway {
            continue;
        }
        stats.spillway_points += traces[ri].pts.len();
        let s = smooth(ri, &mut stats);
        let run_at: HashMap<usize, Run> = runs[ri].iter().map(|r| (r.t0, *r)).collect();
        let mut out = Smoothed { pts: Vec::new(), org: Vec::new(), fixed: Vec::new() };
        let mut skip_to: Option<usize> = None;
        for j in 0..s.pts.len() {
            if let Some(t1) = skip_to {
                if s.fixed[j] && s.org[j] == t1 {
                    skip_to = None;
                } else {
                    continue;
                }
            }
            out.pts.push(s.pts[j]);
            out.org.push(s.org[j]);
            out.fixed.push(s.fixed[j]);
            if !s.fixed[j] {
                continue;
            }
            let Some(run) = run_at.get(&s.org[j]) else { continue };
            let Some(wf) = finals[run.wr].as_ref() else { continue };
            let pa = wf.org.iter().zip(&wf.fixed).position(|(&o, &f)| f && o == run.k0);
            let pb = wf.org.iter().zip(&wf.fixed).position(|(&o, &f)| f && o == run.k1);
            let (Some(pa), Some(pb)) = (pa, pb) else { continue };
            let idx: Vec<usize> = if pb > pa { (pa + 1..pb).collect() } else { (pb + 1..pa).rev().collect() };
            for v in idx {
                out.pts.push(wf.pts[v]);
                // the spillway's own trace index at that watercourse point
                out.org.push((run.t0 + wf.org[v].abs_diff(run.k0)).min(run.t1));
                out.fixed.push(true);
            }
            stats.spliced_points += run.t1 - run.t0 + 1;
            skip_to = Some(run.t1);
        }
        finals[ri] = Some(simplify(out));
    }
    let finals: Vec<Smoothed> = finals.into_iter().map(|s| s.expect("every river smoothed")).collect();
    // the instruments on the final polylines
    if with_stats {
        let (mut st_raw, mut st_raw_n, mut st_sm, mut st_sm_n) = (0usize, 0usize, 0usize, 0usize);
        for ri in 0..nr {
            let t = &traces[ri];
            let n0 = t.pts.len();
            // the nearest trace point over the edge's origin span (± 2)
            let nearest = |q: [f32; 2], o0: usize, o1: usize| -> usize {
                let lo = o0.saturating_sub(2);
                let hi = (o1 + 2).min(n0 - 1);
                (lo..=hi)
                    .min_by(|&a, &b| {
                        let da = (t.pts[a][0] - q[0]).powi(2) + (t.pts[a][1] - q[1]).powi(2);
                        let db = (t.pts[b][0] - q[0]).powi(2) + (t.pts[b][1] - q[1]).powi(2);
                        da.total_cmp(&db)
                    })
                    .unwrap_or(o0.min(n0 - 1))
            };
            let sample = |poly: &[[f32; 2]], orgs: &[usize], above: &mut Vec<f32>, lat: Option<&mut Vec<f32>>| -> (usize, Vec<([f32; 2], bool)>) {
                let mut lat = lat;
                let (mut k, mut cross) = (0usize, Vec::new());
                for i in 0..poly.len().saturating_sub(1) {
                    let (a, b) = (poly[i], poly[i + 1]);
                    let l = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
                    let m = (l / 0.25).ceil().max(1.0) as usize;
                    for j in 0..m {
                        let tt = j as f32 / m as f32;
                        let q = [a[0] + tt * (b[0] - a[0]), a[1] + tt * (b[1] - a[1])];
                        let (o0, o1) = (orgs[i].min(orgs[i + 1]), orgs[i].max(orgs[i + 1]));
                        let nt = nearest(q, o0, o1);
                        above.push(bilinear(&zm, w, h, q[0], q[1]) - t.bed[nt]);
                        if let Some(ol) = outlet.as_ref() {
                            let qc = (q[1].floor() as i64).rem_euclid(h as i64) as usize * w + (q[0].floor() as i64).rem_euclid(w as i64) as usize;
                            if ol[qc] != ol[t.cells[nt]] {
                                let own = t.cells[o0.saturating_sub(3)..=(o1 + 3).min(n0 - 1)].contains(&qc);
                                cross.push((q, own));
                            }
                        }
                        if let Some(lv) = lat.as_deref_mut() {
                            lv.push(lateral_span(&t.pts, q, o0, o1));
                        }
                        k += 1;
                    }
                }
                (k, cross)
            };
            let raw_org: Vec<usize> = (0..n0).collect();
            let (k, c) = sample(&t.pts, &raw_org, &mut stats.above_bed_raw, None);
            stats.samples_raw += k;
            stats.ridge_cross_raw += c.len();
            let mut before = Vec::new();
            let (mut p, mut f, mut o) = (t.pts.clone(), fixed[ri].clone(), raw_org.clone());
            for _ in 0..params.chaikin_iters {
                let r = chaikin(&p, &f, &o);
                p = r.0;
                f = r.1;
                o = r.2;
            }
            sample(&p, &o, &mut before, None);
            stats.above_bed_before.extend(before);
            let fin = &finals[ri];
            let (k, c) = sample(&fin.pts, &fin.org, &mut stats.above_bed_after, Some(&mut stats.lateral_cells));
            stats.samples += k;
            stats.ridge_cross_smoothed += c.len();
            stats.ridge_samples.extend(c.into_iter().map(|(q, own)| (ri as u32, q, own)));
            let (a, b) = stair_share(&t.pts);
            st_raw += a;
            st_raw_n += b;
            let (a, b) = stair_share(&fin.pts);
            st_sm += a;
            st_sm_n += b;
        }
        stats.stair_raw = st_raw as f32 / st_raw_n.max(1) as f32;
        stats.stair_smoothed = st_sm as f32 / st_sm_n.max(1) as f32;
    }
    // attributes
    let cell_km = cell_km2.sqrt();
    let end_of = |ri: usize| -> RiverEnd {
        let chain = &chains[ri];
        let last = *chain.last().expect("non-empty");
        if let Some(d) = segs[last].downstream
            && d < ns
            && seg_river[d] != u32::MAX
        {
            return RiverEnd::Confluence { river_id: seg_river[d] };
        }
        let lc = *traces[ri].cells.last().expect("non-empty");
        let probe: Vec<usize> = std::iter::once(lc).chain(recv(lc)).collect();
        for &c in &probe {
            if field.data[c] <= C1_SEA_LEVEL_NORM && dr.lake_map.get(c).copied().unwrap_or(0) == 0 {
                return RiverEnd::Sea;
            }
            let lid = dr.lake_map.get(c).copied().unwrap_or(0);
            if lid != 0 {
                let endo = dr.lakes.iter().find(|l| l.base.id == lid).is_some_and(|l| l.lake_type == LakeType::Endorheic);
                return if endo { RiverEnd::EndorheicLake { lake_id: lid } } else { RiverEnd::Lake { lake_id: lid } };
            }
        }
        RiverEnd::Terminal
    };
    let mouth_of = |ri: usize| -> u32 {
        let mut c = *traces[ri].cells.last().expect("non-empty");
        let mut last_land = c;
        for _ in 0..200_000 {
            if field.data[c] > C1_SEA_LEVEL_NORM {
                last_land = c;
            } else {
                break;
            }
            match recv(c) {
                Some(r) => c = r,
                None => break,
            }
        }
        last_land as u32
    };
    let vertex_of = |ri: usize| -> VertexAttrs {
        let t = &traces[ri];
        let n0 = t.pts.len();
        let mut v = VertexAttrs::default();
        for &o in &finals[ri].org {
            let o = o.min(n0 - 1);
            let (a, b) = (o.saturating_sub(1), (o + 1).min(n0 - 1));
            let dist = ((t.pts[b][0] - t.pts[a][0]).powi(2) + (t.pts[b][1] - t.pts[a][1]).powi(2)).sqrt() * cell_m;
            let slope = if dist > 0.0 { (t.bed[a] - t.bed[b]) / dist } else { 0.0 };
            let (area, q) = (t.area[o], t.q[o].max(0.0));
            v.catchment_km2.push(sig4(area));
            v.discharge_m3s.push(sig4(q));
            v.slope.push(sig4(slope));
            v.valley_width_m.push(sig4(params.valley_width_coef_m * area.powf(params.valley_width_exp)));
            v.bed_width_m.push(sig4(params.bed_width_a * q.powf(BED_WIDTH_B)));
        }
        v
    };
    let mut rivers: Vec<RiverLl> = (0..nr)
        .map(|ri| {
            let chain = &chains[ri];
            let last = *chain.last().expect("non-empty");
            let kind = kind_of[ri];
            let end = end_of(ri);
            let end_lake = match end {
                RiverEnd::Lake { lake_id } | RiverEnd::EndorheicLake { lake_id } => lake_id,
                _ => 0,
            };
            let mut lakes: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
            for &c in &traces[ri].cells {
                for k in 0..9usize {
                    let q = if k == 8 { c } else { ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize };
                    let lid = dr.lake_map.get(q).copied().unwrap_or(0);
                    if lid != 0 && lid != end_lake {
                        lakes.insert(lid);
                    }
                }
            }
            let pts = finals[ri].pts.clone();
            RiverLl {
                id: ri as u32,
                kind,
                segments: chain.clone(),
                length_km: poly_len(&pts) * cell_km,
                points: pts,
                vertex: vertex_of(ri),
                catchment_km2: catch(last) * cell_km2,
                strahler: (kind == SegmentKind::Watercourse).then(|| chain.iter().map(|&s| segs[s].strahler_order).max().unwrap_or(0)),
                discharge_m3s: dr.segment_discharge_m3s.get(last).copied().unwrap_or(0.0),
                width_m: dr.segment_width_m.get(last).copied().unwrap_or(0.0),
                end,
                lakes_crossed: lakes.into_iter().collect(),
                falls: Vec::new(),
                length_rank_in_basin: 0,
                mouth_basin: mouth_of(ri),
                parallel_of: None,
                fusion_candidate: false,
                retraces_spillway: kind == SegmentKind::Spillway && shared_pts[ri] as f32 >= params.retrace_share * traces[ri].pts.len() as f32,
            }
        })
        .collect();
    // the length rank within the mouth basin
    let mut by_basin: HashMap<u32, Vec<usize>> = HashMap::new();
    for (i, r) in rivers.iter().enumerate() {
        by_basin.entry(r.mouth_basin).or_default().push(i);
    }
    for v in by_basin.values_mut() {
        v.sort_by(|&a, &b| rivers[b].length_km.total_cmp(&rivers[a].length_km).then(a.cmp(&b)));
        for (k, &i) in v.iter().enumerate() {
            rivers[i].length_rank_in_basin = k as u32 + 1;
        }
    }
    // ADR Finding 147-F — the parallel pairs, annotated (nothing removed)
    let mut excluded: HashSet<(u32, u32)> = HashSet::new();
    for r in &rivers {
        if let RiverEnd::Confluence { river_id } = r.end {
            excluded.insert((r.id.min(river_id), r.id.max(river_id)));
        }
    }
    for (ri, set) in retraced_wc.iter().enumerate() {
        for &wr in set {
            excluded.insert(((ri as u32).min(wr as u32), (ri as u32).max(wr as u32)));
        }
    }
    let polys: Vec<&[[f32; 2]]> = rivers.iter().map(|r| r.points.as_slice()).collect();
    let l_cells = (params.parallel_min_km / cell_km).ceil() as usize;
    let pairs = parallel_pairs(&polys, &excluded, params.parallel_k_cells, l_cells);
    let higher = |a: u32, b: u32| -> (u32, u32) {
        let (qa, qb) = (rivers[a as usize].discharge_m3s, rivers[b as usize].discharge_m3s);
        if qa > qb || (qa == qb && a < b) { (a, b) } else { (b, a) }
    };
    let mut longest: Vec<Option<(usize, u32)>> = vec![None; nr];
    let mut candidate = vec![false; nr];
    for &((a, b), len) in &pairs {
        let (_, lo) = higher(a, b);
        candidate[lo as usize] = true;
        for (me, other) in [(a, b), (b, a)] {
            let e = &mut longest[me as usize];
            if e.is_none_or(|(l, o)| len > l || (len == l && other < o)) {
                *e = Some((len, other));
            }
        }
    }
    let parallel_of: Vec<Option<u32>> = (0..nr).map(|ri| longest[ri].map(|(_, o)| higher(ri as u32, o).0)).collect();
    drop(polys);
    for (ri, r) in rivers.iter_mut().enumerate() {
        r.fusion_candidate = candidate[ri];
        r.parallel_of = parallel_of[ri];
    }
    if with_stats {
        stats.parallel_pairs = pairs;
    }
    (rivers, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Chaikin keeps the endpoints and every fixed vertex, and smooths a staircase.
    #[test]
    fn chaikin_keeps_fixed_points_and_smooths_a_staircase() {
        let pts: Vec<[f32; 2]> = (0..12).map(|i| [0.5 + (i / 2) as f32 + (i % 2) as f32, 0.5 + (i / 2) as f32]).collect();
        let mut fixed = vec![false; pts.len()];
        fixed[0] = true;
        fixed[6] = true;
        fixed[11] = true;
        let origin: Vec<usize> = (0..pts.len()).collect();
        let (mut p, mut f, mut o) = (pts.clone(), fixed.clone(), origin);
        for _ in 0..3 {
            let r = chaikin(&p, &f, &o);
            p = r.0;
            f = r.1;
            o = r.2;
        }
        assert_eq!(p[0], pts[0]);
        assert_eq!(*p.last().unwrap(), pts[11]);
        assert!(p.contains(&pts[6]), "the fixed confluence point stays a vertex");
        let (t0, n0) = stair_share(&pts);
        let (t1, n1) = stair_share(&p);
        assert!(t0 as f32 / n0 as f32 > 0.5, "negative control: the staircase turns at every step");
        assert!((t1 as f32 / n1 as f32) < 0.1, "smoothed: almost no turn ≥ 45° ({t1} of {n1})");
    }

    /// ADR Finding 147-X — a step between two fixed vertices stays the straight D8 step (every new point on it); the
    /// negative control: the same step with one end free leaves the line.
    #[test]
    fn chaikin_keeps_a_step_between_two_fixed_vertices_straight() {
        let pts: Vec<[f32; 2]> = vec![[0.5, 0.5], [1.5, 0.5], [2.5, 1.5], [3.5, 1.5], [4.5, 2.5]];
        let off_step = |fixed: &[bool]| -> f32 {
            let origin: Vec<usize> = (0..pts.len()).collect();
            let (mut p, mut f, mut o) = (pts.clone(), fixed.to_vec(), origin);
            for _ in 0..3 {
                let r = chaikin(&p, &f, &o);
                p = r.0;
                f = r.1;
                o = r.2;
            }
            // the points between vertices 1 and 2 (the diagonal step), their distance to that step's line
            let (a, b) = (pts[1], pts[2]);
            p.iter()
                .filter(|q| q[0] > a[0] + 1e-4 && q[0] < b[0] - 1e-4)
                .map(|q| ((b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0])).abs() / 2f32.sqrt())
                .fold(0f32, f32::max)
        };
        let both = off_step(&[true, true, true, false, true]);
        let one = off_step(&[true, false, true, false, true]);
        assert!(both < 1e-5, "two fixed ends: the step stays straight ({both})");
        assert!(one > 0.01, "negative control: one free end cuts the corner ({one})");
    }

    /// The flat-top axial hex of a point: the six neighbours' centres land in their hexes, and the origin's hex is (0, 0).
    #[test]
    fn hex_of_flat_top_axial() {
        let g = HexGrid::living_landz();
        assert_eq!(g.hex_of(0.0, 0.0), (0, 0));
        let s = g.size_m;
        let centre = |q: i32, r: i32| (s * 1.5 * q as f32, s * 3f32.sqrt() * (r as f32 + q as f32 / 2.0));
        for (q, r) in [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1), (3, -2), (-4, 7)] {
            let (x, y) = centre(q, r);
            assert_eq!(g.hex_of(x, y), (q, r), "the centre of ({q}, {r})");
            assert_eq!(g.hex_of(x + 0.4 * s, y - 0.3 * s), (q, r), "inside ({q}, {r})");
        }
    }

    /// RDP keeps the fixed points and the endpoints, and drops collinear ones.
    #[test]
    fn rdp_keeps_fixed_points() {
        let pts: Vec<[f32; 2]> = (0..10).map(|i| [i as f32, 0.0]).collect();
        let mut fixed = vec![false; 10];
        fixed[4] = true;
        assert_eq!(rdp(&pts, &fixed, 0.1), vec![0, 4, 9]);
    }

    /// ADR Finding 147-X, amendment 2 — a reach joining its receiver BELOW its head (the drainage's T3b confluence)
    /// ends in a confluence on it; only a reach joining at the head continues it. The negative control: the mid-joining
    /// reach is the larger one, so the old rule (the largest catchment) would have chained it and jumped back to the
    /// receiver's head.
    #[test]
    fn a_reach_joining_below_the_head_ends_in_a_confluence() {
        use crate::tectonics_c1::drainage::Navigability;
        use crate::terrain::flow::{FlowResult, RiverNetwork, RiverSegment};
        let seg = |points: Vec<(u32, u32)>, upstream: Vec<usize>, downstream: Option<usize>| RiverSegment {
            points,
            strahler_order: 1,
            avg_flow: 1.0,
            max_flow: 1.0,
            basin_id: 1,
            upstream,
            downstream,
        };
        // 0: the receiver, head at (1,1); 1: joins at its head; 2: joins at (3,1), mid-way, and is the larger
        let segs = vec![
            seg(vec![(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)], vec![1, 2], None),
            seg(vec![(1, 4), (1, 3), (1, 2), (1, 1)], vec![], Some(0)),
            seg(vec![(3, 5), (3, 4), (3, 3), (3, 2), (3, 1)], vec![], Some(0)),
        ];
        let n = segs.len();
        let dr = C1DrainageResult {
            flow: FlowResult {
                filled: GridF32::new(8, 8, 0.0),
                direction: vec![DIR_NONE; 64],
                accumulation: GridF32::new(8, 8, 1.0),
                basins: vec![0; 64],
                num_basins: 1,
            },
            segment_drainage_km2: vec![1.0; n],
            segment_navigability: vec![Navigability::NonNavigable; n],
            segment_discharge_m3s: vec![1.0; n],
            segment_width_m: vec![5.0; n],
            segment_profile_m: segs.iter().map(|s| vec![100.0; s.points.len()]).collect(),
            segment_discharge_profile_m3s: segs.iter().map(|s| vec![1.0; s.points.len()]).collect(),
            segment_catchment_cells: vec![20.0, 5.0, 10.0],
            segment_kind: vec![SegmentKind::Watercourse; n],
            segment_source_lake: vec![None; n],
            rivers: RiverNetwork { segments: segs },
            lakes: vec![],
            lake_map: vec![0; 64],
            width: 8,
            height: 8,
        };
        let field = GridF32::new(8, 8, 0.6);
        let (r, st) = build_rivers_ll(&dr, &field, &SteinSteinParams::default(), 1.0, RiversLlParams::default(), true);
        assert_eq!(r.len(), 2, "{:?}", r.iter().map(|x| &x.segments).collect::<Vec<_>>());
        let main = r.iter().find(|x| x.segments.contains(&0)).unwrap();
        assert_eq!(main.segments, vec![1, 0], "the head-joining reach continues the receiver");
        let trib = r.iter().find(|x| x.segments == vec![2]).unwrap();
        assert_eq!(trib.end, RiverEnd::Confluence { river_id: main.id }, "the mid-joining reach ends on it");
        assert_eq!(*trib.points.last().unwrap(), [3.5, 1.5], "at the receiver's junction cell");
        assert!(main.points.contains(&[3.5, 1.5]), "which is a vertex of the receiver");
        assert_eq!((st.mid_joins, st.mid_joins_larger), (1, 1));
        // the receiver's polyline never steps back westward (no jump back to its head)
        assert!(main.points.windows(2).all(|p| p[1][0] >= p[0][0] - 1e-6 || p[1][1] < p[0][1]), "{:?}", main.points);
    }

    /// The 4-significant-digit rounding of the per-vertex attributes.
    #[test]
    fn sig4_rounds_to_four_significant_digits() {
        assert_eq!(sig4(0.0), 0.0);
        assert!((sig4(123.456) - 123.5).abs() < 1e-4);
        assert!((sig4(0.00123456) - 0.001235).abs() < 1e-9);
        assert!((sig4(-98765.4) + 98770.0).abs() < 1e-1);
    }

    /// ADR Finding 147-F — two rivers side by side over a long run are a pair (the lower-discharge one a fusion
    /// candidate); excluded (a shared confluence) they are not; the negative control: far apart, no pair.
    #[test]
    fn parallel_pairs_find_side_by_side_rivers() {
        let a: Vec<[f32; 2]> = (0..60).map(|i| [i as f32 + 0.5, 10.5]).collect();
        let b: Vec<[f32; 2]> = (0..60).map(|i| [i as f32 + 0.5, 11.5]).collect();
        let c: Vec<[f32; 2]> = (0..60).map(|i| [i as f32 + 0.5, 40.5]).collect();
        let polys: Vec<&[[f32; 2]]> = vec![&a, &b, &c];
        let p = parallel_pairs(&polys, &HashSet::new(), 2, 40);
        assert_eq!(p.len(), 1, "{p:?}");
        assert_eq!(p[0].0, (0, 1));
        assert!(p[0].1 > 40);
        let ex: HashSet<(u32, u32)> = [(0, 1)].into_iter().collect();
        assert!(parallel_pairs(&polys, &ex, 2, 40).is_empty(), "an excluded pair is not reported");
    }
}
