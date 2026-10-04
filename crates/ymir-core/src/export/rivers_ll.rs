//! `rivers_ll.json` — the rivers for Living Landz (ADR Finding 146).
//!
//! A NEW vector layer, alongside the unchanged `rivers.json`: the same drainage, regrouped into **main stems** (at each
//! confluence the upstream reach with the largest catchment continues the river; Hack's convention), each smoothed
//! into a polyline that stays in its valley floor, with the attributes a selection needs. No threshold is applied:
//! the selection is the consumer's.
//!
//! - **Smoothing**: Chaikin corner cutting on the D8 trace's cell centres, both endpoints and every confluence point
//!   held fixed, so a tributary's last point is exactly a vertex of the river it joins. Each smoothed vertex must
//!   stand no more than `valley_tol_m` above the bed of its nearest trace point (bilinear terrain on the conditioned
//!   field); a vertex that fails is pulled back towards that trace point by bisection. Then a Ramer–Douglas–Peucker
//!   pass (fixed points kept) bounds the point count.
//! - **Coordinates**: continuous erosion-grid cells, the cell centre at `+0.5`, `y = 0` the SOUTH edge (the container's
//!   orientation invariant, `export/container.rs`).
//! - The format is documented in `docs/rivers_ll_format.md`; `format_version` is semver.

use serde::Serialize;

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, C1DrainageResult, LakeType, SegmentKind};
use crate::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

/// The format version this writer emits.
pub const RIVERS_LL_FORMAT_VERSION: &str = "0.2.0";

/// The smoothing's parameters (declared in Finding 146).
#[derive(Clone, Copy, Debug)]
pub struct RiversLlParams {
    /// Chaikin iterations.
    pub chaikin_iters: u32,
    /// A smoothed vertex may stand at most this far above its nearest trace point's bed (m).
    pub valley_tol_m: f32,
    /// Bisection steps when a vertex is pulled back.
    pub bisect_steps: u32,
    /// Ramer–Douglas–Peucker tolerance (cells) applied after the smoothing; 0 = none.
    pub rdp_eps_cells: f32,
}

impl Default for RiversLlParams {
    fn default() -> Self {
        Self { chaikin_iters: 3, valley_tol_m: 1.0, bisect_steps: 8, rdp_eps_cells: 0.1 }
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

/// One river (a main stem).
#[derive(Clone, Debug, Serialize)]
pub struct RiverLl {
    pub id: u32,
    pub kind: SegmentKind,
    /// The `rivers.json` segments it chains, upstream → downstream.
    pub segments: Vec<usize>,
    /// The smoothed polyline, upstream → downstream, in continuous erosion-grid cells.
    pub points: Vec<[f32; 2]>,
    pub length_km: f32,
    /// Geometric catchment at the mouth (km²).
    pub catchment_km2: f32,
    /// The largest Strahler order among its segments; `None` for a spillway (no hierarchy).
    pub strahler: Option<u8>,
    /// Mean discharge at the mouth (m³/s) and the bankfull width there (m).
    pub discharge_m3s: f32,
    pub width_m: f32,
    pub end: RiverEnd,
    /// Lakes touching its course (other than the one it ends in).
    pub lakes_crossed: Vec<u32>,
    pub falls: Vec<FallLl>,
    /// 1 = the longest river of its mouth basin.
    pub length_rank_in_basin: u32,
    /// The number of hex edges it crosses (0 when no hex grid is given), for a consumer's selection counts.
    pub hex_edges: u32,
    /// The mouth basin's key: the last land cell (row-major index) on the D8 path from the river's end.
    pub mouth_basin: u32,
}

#[derive(Serialize)]
struct RiversLlView<'a> {
    format_version: &'static str,
    coordinate_space: &'static str,
    km_per_cell: f32,
    rivers: &'a [RiverLl],
    #[serde(skip_serializing_if = "Option::is_none")]
    hex_grid: Option<&'a HexGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    hex_edges: Option<&'a HexEdges>,
}

/// ADR Finding 146-H — Living Landz's hex grid (the author, 2026-10-05): flat-top hexes of `size_m` (centre to corner),
/// axial coordinates (q, r), hex (0, 0) centred on the map's bottom-left corner (`origin_m` = (0, 0) in metres from
/// that corner, `y` northward, as the container's `y = 0` = south).
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

/// ADR Finding 146-H — the hex edges the rivers cross. Each river's polyline is walked every 0.1 hex size; a change of
/// hex is a crossing (a jump over a non-neighbour is subdivided). Sets each river's `hex_edges` count. `width` is the
/// construction's law W(A) = `coef · A^exp` (the témoin: 200 m at 10 km², exponent 0.3).
#[allow(clippy::too_many_arguments)]
pub fn river_hex_edges(rivers: &mut [RiverLl], dr: &C1DrainageResult, field: &GridF32, ss: &SteinSteinParams, cell_km2: f32, hex: &HexGrid, width_coef_m: f32, width_exp: f32) -> HexEdges {
    let (w, h) = (field.width, field.height);
    let cell_m = cell_km2.sqrt() * 1000.0;
    let zm: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
    let acc = &dr.flow.accumulation.data;
    let mut e = HexEdges::default();
    let step = 0.1 * hex.size_m / cell_m; // in cells
    for r in rivers.iter_mut() {
        let p = &r.points;
        if p.len() < 2 {
            continue;
        }
        // the polyline resampled every `step` cells, with each sample's arc length
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
        let mut count = 0u32;
        for i in 1..s.len() {
            let (a, b) = (hx[i - 1], hx[i]);
            if a == b {
                continue;
            }
            let c = [0.5 * (s[i - 1][0] + s[i][0]), 0.5 * (s[i - 1][1] + s[i][1])];
            let cell = |q: [f32; 2]| -> usize { (q[1].floor() as i64).rem_euclid(h as i64) as usize * w + (q[0].floor() as i64).rem_euclid(w as i64) as usize };
            let cc = cell(c);
            // the catchment: the largest accumulation in the crossing cell's 3×3 (the line sits within ~0.1 cell of the
            // D8 trace), capped at the river's mouth
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
            count += 1;
        }
        r.hex_edges = count;
    }
    e
}

/// The instruments of Finding 146-L (the bench reads them; production ignores them).
#[derive(Clone, Debug, Default)]
pub struct RiversLlStats {
    /// terrain − bed (m) over the samples, before / after the valley constraint.
    pub above_bed_before: Vec<f32>,
    pub above_bed_after: Vec<f32>,
    /// The same on the raw D8 trace (the negative control).
    pub above_bed_raw: Vec<f32>,
    /// Lateral distance (cells) of each smoothed sample from the D8 trace.
    pub lateral_cells: Vec<f32>,
    /// Samples whose cell drains to another outlet than the river's own trace cell: smoothed / raw.
    pub ridge_cross_smoothed: usize,
    pub ridge_cross_raw: usize,
    pub samples: usize,
    pub samples_raw: usize,
    /// Vertices pulled back by the constraint, of all smoothed vertices.
    pub pulled: usize,
    pub vertices: usize,
    /// The share of turns ≥ 45° on the raw trace and on the smoothed polyline.
    pub stair_raw: f32,
    pub stair_smoothed: f32,
}

/// Serialize the rivers as `rivers_ll.json` bytes, with the hex grid and its edges when given.
pub fn rivers_ll_json(rivers: &[RiverLl], km_per_cell: f32, hex: Option<(&HexGrid, &HexEdges)>) -> Vec<u8> {
    let view = RiversLlView {
        format_version: RIVERS_LL_FORMAT_VERSION,
        coordinate_space: "erosion_grid_cells_continuous",
        km_per_cell,
        rivers,
        hex_grid: hex.map(|x| x.0),
        hex_edges: hex.map(|x| x.1),
    };
    serde_json::to_vec(&view).expect("rivers_ll serialization is infallible")
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
/// index.
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

/// Ramer–Douglas–Peucker keeping the fixed points.
fn rdp(pts: &[[f32; 2]], fixed: &[bool], eps: f32) -> Vec<[f32; 2]> {
    if pts.len() <= 2 || eps <= 0.0 {
        return pts.to_vec();
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
    pts.iter().zip(keep.iter()).filter(|(_, k)| **k).map(|(p, _)| *p).collect()
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

/// Build the main stems, smooth them in their valley and attach the attributes. `field` is the conditioned
/// (breached) field the drainage was computed on (normalised). With `with_stats`, the Finding 146-L instruments are
/// filled (an O(n) outlet labelling more).
pub fn build_rivers_ll(dr: &C1DrainageResult, field: &GridF32, ss: &SteinSteinParams, cell_km2: f32, params: RiversLlParams, with_stats: bool) -> (Vec<RiverLl>, RiversLlStats) {
    let (w, h) = (field.width, field.height);
    let segs = &dr.rivers.segments;
    let ns = segs.len();
    let zm: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
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
    // the continuing upstream reach of each segment: the largest catchment (ties: the lowest index)
    let mut cont_up = vec![usize::MAX; ns];
    for d in 0..ns {
        let mut best: Option<usize> = None;
        for &u in &segs[d].upstream {
            if u < ns && best.is_none_or(|b| catch(u) > catch(b) || (catch(u) == catch(b) && u < b)) {
                best = Some(u);
            }
        }
        if let Some(b) = best {
            cont_up[d] = b;
        }
    }
    // river chains: from each segment that no downstream continues from it… i.e. heads of chains
    let mut is_cont = vec![false; ns];
    for d in 0..ns {
        if cont_up[d] != usize::MAX {
            is_cont[cont_up[d]] = true;
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
    // the trace of each river (cell centres) with its bed, and its confluence end point
    let mut traces: Vec<(Vec<[f32; 2]>, Vec<f32>, Vec<usize>)> = Vec::with_capacity(chains.len());
    for chain in &chains {
        let (mut pts, mut bed, mut cells): (Vec<[f32; 2]>, Vec<f32>, Vec<usize>) = (Vec::new(), Vec::new(), Vec::new());
        for &s in chain {
            let prof = dr.segment_profile_m.get(s);
            for (i, &(x, y)) in segs[s].points.iter().enumerate() {
                let c = y as usize * w + x as usize;
                if cells.last() == Some(&c) {
                    continue;
                }
                cells.push(c);
                pts.push([x as f32 + 0.5, y as f32 + 0.5]);
                bed.push(prof.and_then(|p| p.get(i).copied()).unwrap_or(zm[c]));
            }
        }
        // a river ending in a confluence: its last point must be the receiving river's cell
        let last = *chain.last().expect("non-empty chain");
        if let Some(d) = segs[last].downstream
            && d < ns
            && let Some(&lc) = cells.last()
        {
            let dcells: Vec<usize> = segs[d].points.iter().map(|&(x, y)| y as usize * w + x as usize).collect();
            if !dcells.contains(&lc) {
                let r = recv(lc).filter(|r| dcells.contains(r)).or_else(|| dcells.first().copied());
                if let Some(r) = r {
                    cells.push(r);
                    pts.push([(r % w) as f32 + 0.5, (r / w) as f32 + 0.5]);
                    bed.push(zm[r]);
                }
            }
        }
        traces.push((pts, bed, cells));
    }
    // the confluence points each river must keep as vertices: where a tributary ends on it
    let mut fixed_cells: Vec<std::collections::HashSet<usize>> = vec![std::collections::HashSet::new(); chains.len()];
    for (ri, chain) in chains.iter().enumerate() {
        let last = *chain.last().expect("non-empty");
        if let Some(d) = segs[last].downstream
            && d < ns
            && let Some(&lc) = traces[ri].2.last()
        {
            let rj = seg_river[d];
            if rj != u32::MAX {
                fixed_cells[rj as usize].insert(lc);
            }
        }
    }
    let tol = params.valley_tol_m;
    let mut stats = RiversLlStats::default();
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
    let (mut st_raw, mut st_raw_n, mut st_sm, mut st_sm_n) = (0usize, 0usize, 0usize, 0usize);
    let mut polylines: Vec<Vec<[f32; 2]>> = Vec::with_capacity(chains.len());
    for (ri, (pts, bed, cells)) in traces.iter().enumerate() {
        let n0 = pts.len();
        let fixed0: Vec<bool> = (0..n0).map(|i| i == 0 || i + 1 == n0 || fixed_cells[ri].contains(&cells[i])).collect();
        let origin0: Vec<usize> = (0..n0).collect();
        let (mut p, mut f, mut o) = (pts.clone(), fixed0.clone(), origin0);
        for _ in 0..params.chaikin_iters {
            let r = chaikin(&p, &f, &o);
            p = r.0;
            f = r.1;
            o = r.2;
        }
        // the nearest trace point of a smoothed point (its origin ± 2)
        let nearest = |q: [f32; 2], org: usize| -> usize {
            let lo = org.saturating_sub(2);
            let hi = (org + 2).min(n0 - 1);
            (lo..=hi).min_by(|&a, &b| {
                let da = (pts[a][0] - q[0]).powi(2) + (pts[a][1] - q[1]).powi(2);
                let db = (pts[b][0] - q[0]).powi(2) + (pts[b][1] - q[1]).powi(2);
                da.total_cmp(&db)
            }).unwrap_or(org.min(n0 - 1))
        };
        let above = |q: [f32; 2], org: usize| -> f32 {
            let t = nearest(q, org);
            bilinear(&zm, w, h, q[0], q[1]) - bed[t]
        };
        let sample = |poly: &[[f32; 2]], orgs: &[usize], out: &mut Vec<f32>, cross: &mut usize, lateral: Option<&mut Vec<f32>>| -> usize {
            let mut lat = lateral;
            let mut k = 0usize;
            for i in 0..poly.len().saturating_sub(1) {
                let (a, b) = (poly[i], poly[i + 1]);
                let l = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
                let m = (l / 0.25).ceil().max(1.0) as usize;
                for j in 0..m {
                    let t = j as f32 / m as f32;
                    let q = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])];
                    let org = if t < 0.5 { orgs[i] } else { orgs[i + 1] };
                    out.push(above(q, org));
                    if let Some(ol) = outlet.as_ref() {
                        let qc = (q[1].floor() as i64).rem_euclid(h as i64) as usize * w + (q[0].floor() as i64).rem_euclid(w as i64) as usize;
                        let tc = cells[nearest(q, org)];
                        if ol[qc] != ol[tc] {
                            *cross += 1;
                        }
                    }
                    if let Some(lv) = lat.as_deref_mut() {
                        // distance to the trace polyline near the origin
                        let lo = org.saturating_sub(3);
                        let hi = (org + 3).min(n0 - 1);
                        let mut dmin = f32::INFINITY;
                        for s in lo..hi {
                            let (u, v) = (pts[s], pts[s + 1]);
                            let (dx, dy) = (v[0] - u[0], v[1] - u[1]);
                            let l2 = (dx * dx + dy * dy).max(1e-12);
                            let tt = (((q[0] - u[0]) * dx + (q[1] - u[1]) * dy) / l2).clamp(0.0, 1.0);
                            let d = ((u[0] + tt * dx - q[0]).powi(2) + (u[1] + tt * dy - q[1]).powi(2)).sqrt();
                            dmin = dmin.min(d);
                        }
                        if dmin.is_finite() {
                            lv.push(dmin);
                        }
                    }
                    k += 1;
                }
            }
            k
        };
        if with_stats {
            let orgs_raw: Vec<usize> = (0..n0).collect();
            let mut c = 0usize;
            stats.samples_raw += sample(pts, &orgs_raw, &mut stats.above_bed_raw, &mut c, None);
            stats.ridge_cross_raw += c;
            let mut c2 = 0usize;
            let mut tmp = Vec::new();
            sample(&p, &o, &mut tmp, &mut c2, None);
            stats.above_bed_before.extend(tmp);
            let (a, b) = stair_share(pts);
            st_raw += a;
            st_raw_n += b;
        }
        // the valley constraint
        for i in 0..p.len() {
            stats.vertices += 1;
            if f[i] {
                continue;
            }
            if above(p[i], o[i]) > tol {
                stats.pulled += 1;
                let t = pts[nearest(p[i], o[i])];
                let (mut lo, mut hi) = (0f32, 1f32); // fraction towards the trace point: 0 = the vertex, 1 = the trace point
                for _ in 0..params.bisect_steps {
                    let mid = 0.5 * (lo + hi);
                    let q = [p[i][0] + mid * (t[0] - p[i][0]), p[i][1] + mid * (t[1] - p[i][1])];
                    if above(q, o[i]) > tol {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                p[i] = [p[i][0] + hi * (t[0] - p[i][0]), p[i][1] + hi * (t[1] - p[i][1])];
            }
        }
        if with_stats {
            let mut c3 = 0usize;
            stats.samples += sample(&p, &o, &mut stats.above_bed_after, &mut c3, Some(&mut stats.lateral_cells));
            stats.ridge_cross_smoothed += c3;
        }
        let simp = rdp(&p, &f, params.rdp_eps_cells);
        if with_stats {
            let (a, b) = stair_share(&simp);
            st_sm += a;
            st_sm_n += b;
        }
        polylines.push(simp);
    }
    stats.stair_raw = st_raw as f32 / st_raw_n.max(1) as f32;
    stats.stair_smoothed = st_sm as f32 / st_sm_n.max(1) as f32;
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
        let lc = *traces[ri].2.last().expect("non-empty");
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
        let mut c = *traces[ri].2.last().expect("non-empty");
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
    let mut rivers: Vec<RiverLl> = (0..chains.len())
        .map(|ri| {
            let chain = &chains[ri];
            let last = *chain.last().expect("non-empty");
            let kind = dr.segment_kind.get(last).copied().unwrap_or(SegmentKind::Watercourse);
            let end = end_of(ri);
            let end_lake = match end {
                RiverEnd::Lake { lake_id } | RiverEnd::EndorheicLake { lake_id } => lake_id,
                _ => 0,
            };
            let mut lakes: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
            for &c in &traces[ri].2 {
                for k in 0..9usize {
                    let q = if k == 8 { c } else { ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize };
                    let lid = dr.lake_map.get(q).copied().unwrap_or(0);
                    if lid != 0 && lid != end_lake {
                        lakes.insert(lid);
                    }
                }
            }
            RiverLl {
                id: ri as u32,
                kind,
                segments: chain.clone(),
                points: polylines[ri].clone(),
                length_km: poly_len(&polylines[ri]) * cell_km,
                catchment_km2: catch(last) * cell_km2,
                strahler: (kind == SegmentKind::Watercourse).then(|| chain.iter().map(|&s| segs[s].strahler_order).max().unwrap_or(0)),
                discharge_m3s: dr.segment_discharge_m3s.get(last).copied().unwrap_or(0.0),
                width_m: dr.segment_width_m.get(last).copied().unwrap_or(0.0),
                end,
                lakes_crossed: lakes.into_iter().collect(),
                falls: Vec::new(),
                length_rank_in_basin: 0,
                hex_edges: 0,
                mouth_basin: mouth_of(ri),
            }
        })
        .collect();
    // the length rank within the mouth basin
    let mut by_basin: std::collections::HashMap<u32, Vec<usize>> = std::collections::HashMap::new();
    for (i, r) in rivers.iter().enumerate() {
        by_basin.entry(r.mouth_basin).or_default().push(i);
    }
    for v in by_basin.values_mut() {
        v.sort_by(|&a, &b| rivers[b].length_km.total_cmp(&rivers[a].length_km).then(a.cmp(&b)));
        for (k, &i) in v.iter().enumerate() {
            rivers[i].length_rank_in_basin = k as u32 + 1;
        }
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
        let s = rdp(&pts, &fixed, 0.1);
        assert_eq!(s, vec![pts[0], pts[4], pts[9]]);
    }
}
