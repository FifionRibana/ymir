//! ADR 0001 Finding 121 — **valley construction: the SHAPE stage of the hybrid** (construct the
//! mature form along the trunks, then let a light incision supply the texture).
//!
//! **A bench / viz seam. `None` in production, byte-identical, nothing promoted.**
//!
//! ## Why this is not Option A again
//!
//! The C1 design document (`docs/c1_lightweight_dynamic_tectonics.md` §3.1) rejected "static with
//! invented rules — geometric templates modulated by procedural noise" for three reasons, copied
//! verbatim so the hybrid answers each:
//!
//! 1. *"cannot produce coherent chronological diversity (young vs old mountain ranges on the same
//!    continent)"*;
//! 2. *"parameter calibration drifts indefinitely with no physical anchor"*;
//! 3. *"'looks invented' is a real failure mode for users who can identify imposed geometry by
//!    eye"*.
//!
//! And ADR Finding 96 built the mature form as a FIELD (χ everywhere) and rejected it: 7.1× rougher
//! than the delivered field, cliffs across every divide, D8's axes stamped into the relief. This
//! module differs on each count: the valleys are laid ONLY along trunks of `A ≥ a_min_km2` (Finding
//! 120: 10 km²), the interfluves keep the tectonic + FBM field (`min(field, valley)`, so no divide
//! is ever built), and the cross-section is laid by the distance to a SMOOTHED skeleton, not to the
//! D8 cells.
//!
//! ## The laws, labelled
//!
//! | law | form | value | status |
//! |---|---|---|---|
//! | floor profile | `z = z_base + k·χ`, `χ = ∫ (A0/A)^0.5 dx`, A0 = 1 km², A clamped at `a_c_km2` | θ = 0.5 | ANCHORED — Harel et al. 2016 eqs. (4)–(7), PDF p. 15 |
//! | age `k` | one number, m/m at A0 | calibrated so the floor surface's paired p50 = 488 m (Finding 95's oracle; Finding 96's own calibration) | **PROXY** on a target |
//! | valley floor width | `W = a·A^0.3`, A in km², W in m | exponent 0.3 | exponent ANCHORED — Clubb et al. 2022 p. 437 |
//! | width coefficient `a` | | `W(10 km²) = 200 m` | **PROXY** (Clubb Table 2 p. 452 is unit-inconsistent, Finding 120) |
//! | walls | planar, 28° | | ANCHORED — Whipple & Tucker 1999 Table 1 p. 17,663 (colluvial slopes) |
//! | base | sea + 0.5 m | | ANCHORED — ADR Finding 83 |
//! | skeleton smoothing | moving average over ±`smooth_m` | 250 m | **PROXY** — against Finding 96's D8 stamping |
//! | trunk threshold | `A ≥ 10 km²` (4 194 cells at 8192²) | | DECISION — Finding 120 |
//!
//! Everything is in physical units: km², m, degrees. Nothing is expressed in cells or in Strahler
//! order (the author's amendment to Finding 120).

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, C1DrainageConfig, c1_drainage_windowed};
use crate::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, c1_metres_to_altitude_norm,
};
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, breach_monotone, compute_flow};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// **PROXY — the age `k` calibrated by `f121_hybrid`** (seed 1, 8192², 400 km): the floor surface
/// `base + k·χ` over the paired land of the pre-incision field has p50 = **488.1 m**, Finding 95's
/// oracle. Finding 96 calibrated its χ field on the same oracle and found 0.0719 with `A` clamped
/// at `A_c`; this skeleton gives 0.07186. A calibration on a TARGET, not a measurement of an age.
pub const F121_AGE_K: f32 = 0.07186;

/// The construction's parameters, all in physical units. See the module table for each label.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ValleyConstruction {
    /// Trunk threshold, km² of D8 drained area (in-domain, never the geo-ratio-signified area).
    /// DECISION (Finding 120).
    pub a_min_km2: f32,
    /// The age `k` of `z = z_base + k·χ` (m per m of χ, A0 = 1 km²). PROXY.
    pub age_k: f32,
    /// The drained area χ is clamped at (km²): Finding 96 B1′, the fluvial law is not applied
    /// below its own channel threshold. ANCHORED (Whipple & Tucker Table 1, via the anchoring
    /// report).
    pub a_c_km2: f32,
    /// `a` in `W = a·A^b` (m, A in km²). PROXY.
    pub width_coef_m: f32,
    /// `b` in `W = a·A^b`. ANCHORED (Clubb 2022 p. 437).
    pub width_exp: f32,
    /// Wall slope, degrees. ANCHORED (Whipple & Tucker 1999 p. 17,663).
    pub wall_deg: f32,
    /// Base level above the sea, m. ANCHORED (Finding 83).
    pub base_m: f32,
    /// Half-length of the skeleton smoothing, m. PROXY.
    pub smooth_m: f32,
    /// `None` = the BARE construction (no incision at all). `Some(f)` = ONE incision pass at
    /// `f × k_time` of the shipped stream power, AFTER the construction (the hybrid).
    pub light_k_time_fraction: Option<f32>,
    /// **ADR Finding 122-B — χ from the col of a closed depression.** `false` = Finding 121 (χ from
    /// the sea everywhere, byte-identical). `true`: a trunk whose D8 path enters a closed depression
    /// of the input field — an inland below-sea basin (`water_class` 2, Finding 85's spill level) or
    /// a lake the pre-incision drainage holds flat — takes that depression's SPILL LEVEL as its base
    /// and counts χ from the depression's col, so a constructed floor never lies below the water
    /// the depression will hold. Outside closed depressions nothing changes. Finding 122-A attributed
    /// the author's "brush" lake to exactly that: valleys built down to sea + 0.5 m inside a basin
    /// whose col stays at 459 m.
    #[serde(default)]
    pub basin_base: bool,
}

impl ValleyConstruction {
    /// Finding 121's values, with the age and the light-incision fraction supplied.
    pub fn f121(age_k: f32, light_k_time_fraction: Option<f32>) -> Self {
        Self {
            a_min_km2: 10.0,
            age_k,
            a_c_km2: 0.1,
            // W(10 km²) = 200 m ⇒ a = 200 / 10^0.3
            width_coef_m: 200.0 / 10f32.powf(0.3),
            width_exp: 0.3,
            wall_deg: 28.0,
            base_m: 0.5,
            smooth_m: 250.0,
            light_k_time_fraction,
            basin_base: false,
        }
    }

    /// Finding 121's values with Finding 122-B's `basin_base` on.
    pub fn f122(age_k: f32, light_k_time_fraction: Option<f32>) -> Self {
        Self { basin_base: true, ..Self::f121(age_k, light_k_time_fraction) }
    }

    /// Valley floor width (m) at drained area `a_km2`.
    pub fn width_m(&self, a_km2: f32) -> f32 {
        self.width_coef_m * a_km2.max(0.0).powf(self.width_exp)
    }
}

/// The skeleton a construction stands on: the D8 network of the BREACHED input field (Finding
/// 120's A0 instrument, reproduced exactly), χ along it, and the smoothed trunk polylines.
pub struct Skeleton {
    pub width: usize,
    pub height: usize,
    pub cell_m: f32,
    /// D8 drained area, km², per cell.
    pub area_km2: Vec<f32>,
    /// χ (m) of each land cell along its D8 path to its base; NaN on sea cells.
    pub chi_m: Vec<f32>,
    /// Altitude (m) of each land cell's base: `base_m` above the sea, or the altitude of a land
    /// terminal (a D8 cell with no receiver).
    pub base_alt_m: Vec<f32>,
    /// Trunk cells: land, `A ≥ a_min_km2`.
    pub trunk: Vec<bool>,
    /// Smoothed, densified trunk polylines: (x, y) in continuous cell units (cell centres at
    /// +0.5, UNWRAPPED along each line), floor altitude (m), half floor width (m).
    pub polylines: Vec<Vec<(f32, f32, f32, f32)>>,
}

impl Skeleton {
    /// Floor altitude (m) of a land cell at age `k`.
    pub fn floor_m(&self, k: usize, age_k: f32) -> f32 {
        self.base_alt_m[k] + age_k * self.chi_m[k]
    }
}

fn recv(dir: &[u8], k: usize, w: usize, h: usize) -> Option<usize> {
    let d = dir[k];
    if d == DIR_NONE {
        return None;
    }
    let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
    let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
    Some(ny * w + nx)
}

/// Build the skeleton of `field` (normalised altitude). `domain_km` is the span of the grid.
///
/// The flow is Finding 120's: `c1_drainage_windowed` → `breach_monotone` → the drainage chain's
/// own `compute_flow` on the breached field (same sea level, same flat perturbation).
pub fn skeleton(
    field: &GridF32,
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
    domain_km: f32,
) -> Skeleton {
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let cell_m = domain_km * 1000.0 / w as f32;
    let cell_km2 = (cell_m / 1000.0) * (cell_m / 1000.0);
    let dcfg = C1DrainageConfig::default();
    let bf = {
        let pre = c1_drainage_windowed(field, None, &dcfg, ss, domain_km);
        breach_monotone(field, &pre.flow.filled, &pre.lake_map, C1_SEA_LEVEL_NORM, w, h)
    };
    let flow = compute_flow(
        &bf,
        &FlowConfig {
            sea_level: C1_SEA_LEVEL_NORM,
            flat_perturbation: dcfg.flat_perturbation.clone(),
            dinf: dcfg.dinf,
        },
    );
    let land: Vec<bool> = bf.data.iter().map(|&v| v > C1_SEA_LEVEL_NORM).collect();
    let area_km2: Vec<f32> = flow.accumulation.data.iter().map(|&a| a * cell_km2).collect();
    let dir = &flow.direction;
    // ADR Finding 122-B -- the closed depressions and their spill levels: a priority flood seeded
    // on the OPEN OCEAN only, so an inland below-sea basin and a held lake both rise to their col.
    let spill: Option<Vec<f32>> = vc.basin_base.then(|| ocean_flood(&bf));
    let eps = 0.01 / c1_altitude_norm_to_metres(1.0, ss).max(1.0); // 1 cm, in norm units
    let dep = |k: usize| -> bool { spill.as_ref().is_some_and(|f| f[k] > bf.data[k] + eps) };
    let spill_m = |k: usize| -> f32 {
        c1_altitude_norm_to_metres(spill.as_ref().expect("basin_base")[k], ss)
    };

    // χ, integrated from each base upstream (memoised walk down the receivers)
    let mut chi = vec![f32::NAN; n];
    let mut base = vec![f32::NAN; n];
    let link_m = |a: usize, b: usize| -> f32 {
        let diag = (a % w != b % w) && (a / w != b / w);
        if diag { cell_m * std::f32::consts::SQRT_2 } else { cell_m }
    };
    let mut path: Vec<usize> = Vec::new();
    for s in 0..n {
        if !land[s] || !chi[s].is_nan() {
            continue;
        }
        path.clear();
        let mut c = s;
        let (b, mut x);
        loop {
            if !chi[c].is_nan() {
                b = base[c];
                x = chi[c];
                break;
            }
            // ADR Finding 122-B -- inside a closed depression: walk to its col; the base is the
            // spill level and χ is counted from the col.
            if land[c] && dep(c) {
                let lvl = spill_m(c);
                let mut terminal = false;
                loop {
                    path.push(c);
                    match recv(dir, c, w, h) {
                        None => {
                            terminal = true;
                            break;
                        }
                        Some(r) if !land[r] || !dep(r) => break, // the col (or the sea)
                        Some(r) => c = r,
                    }
                    if path.len() > n {
                        panic!("D8 cycle in a closed depression of the valley skeleton");
                    }
                }
                if terminal {
                    let t = path.pop().expect("pushed");
                    chi[t] = 0.0;
                    base[t] = lvl;
                }
                b = lvl;
                x = 0.0;
                break;
            }
            match recv(dir, c, w, h) {
                None => {
                    // a land terminal: its own altitude is the base
                    b = c1_altitude_norm_to_metres(bf.data[c], ss);
                    x = 0.0;
                    chi[c] = 0.0;
                    base[c] = b;
                    break;
                }
                Some(r) => {
                    path.push(c);
                    if !land[r] {
                        // ADR Finding 122-B -- an INLAND below-sea basin is not the sea: its base
                        // is its spill level. The open ocean keeps Finding 83's `base_m`.
                        b = if dep(r) { spill_m(r) } else { vc.base_m };
                        x = 0.0;
                        // the link c → r is integrated in the unwind below, with r = sea
                        break;
                    }
                    if path.len() > n {
                        panic!("D8 cycle in the valley skeleton");
                    }
                    c = r;
                }
            }
        }
        for i in (0..path.len()).rev() {
            let c = path[i];
            let r = recv(dir, c, w, h).expect("a path cell has a receiver");
            let a = area_km2[c].max(vc.a_c_km2);
            x += (1.0 / a).sqrt() * link_m(c, r);
            chi[c] = x;
            base[c] = b;
        }
    }

    // trunks and their polylines
    let trunk: Vec<bool> = (0..n).map(|k| land[k] && area_km2[k] >= vc.a_min_km2).collect();
    let mut donors = vec![0u8; n];
    for k in 0..n {
        if trunk[k]
            && let Some(r) = recv(dir, k, w, h)
            && trunk[r]
        {
            donors[r] = donors[r].saturating_add(1);
        }
    }
    let mut visited = vec![false; n];
    let mut raw_lines: Vec<Vec<usize>> = Vec::new();
    for head in 0..n {
        if !trunk[head] || donors[head] != 0 || visited[head] {
            continue;
        }
        let mut line = Vec::new();
        let mut c = head;
        loop {
            line.push(c);
            if visited[c] {
                break; // a junction: the line ends ON the parent's cell
            }
            visited[c] = true;
            match recv(dir, c, w, h) {
                Some(r) if trunk[r] => c = r,
                _ => break,
            }
        }
        if line.len() >= 2 {
            raw_lines.push(line);
        }
    }
    let half_k = (vc.smooth_m / cell_m).round().max(0.0) as usize;
    let polylines = raw_lines
        .iter()
        .map(|line| {
            // unwrap the torus along the line
            let mut pts: Vec<(f32, f32, f32, f32)> = Vec::with_capacity(line.len());
            let (mut px, mut py) = (0f32, 0f32);
            for (i, &c) in line.iter().enumerate() {
                let (mut x, mut y) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
                if i > 0 {
                    while x - px > w as f32 / 2.0 {
                        x -= w as f32;
                    }
                    while px - x > w as f32 / 2.0 {
                        x += w as f32;
                    }
                    while y - py > h as f32 / 2.0 {
                        y -= h as f32;
                    }
                    while py - y > h as f32 / 2.0 {
                        y += h as f32;
                    }
                }
                px = x;
                py = y;
                let zf = base[c] + vc.age_k * chi[c];
                pts.push((x, y, zf, 0.5 * vc.width_m(area_km2[c])));
            }
            // moving average of the POSITION only, endpoints pinned
            let sm = if half_k > 0 && pts.len() >= 2 * half_k + 3 {
                let mut o = pts.clone();
                for i in half_k..pts.len() - half_k {
                    let (mut sx, mut sy) = (0f32, 0f32);
                    for p in &pts[i - half_k..=i + half_k] {
                        sx += p.0;
                        sy += p.1;
                    }
                    let d = (2 * half_k + 1) as f32;
                    o[i].0 = sx / d;
                    o[i].1 = sy / d;
                }
                o
            } else {
                pts
            };
            // densify to ≤ 0.5 cell between samples
            let mut dense = Vec::with_capacity(sm.len() * 2);
            for i in 0..sm.len() {
                dense.push(sm[i]);
                if i + 1 < sm.len() {
                    let (a, b) = (sm[i], sm[i + 1]);
                    let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
                    let steps = (len / 0.5).ceil() as usize;
                    for s in 1..steps {
                        let t = s as f32 / steps as f32;
                        dense.push((
                            a.0 + t * (b.0 - a.0),
                            a.1 + t * (b.1 - a.1),
                            a.2 + t * (b.2 - a.2),
                            a.3 + t * (b.3 - a.3),
                        ));
                    }
                }
            }
            dense
        })
        .collect();
    Skeleton {
        width: w,
        height: h,
        cell_m,
        area_km2,
        chi_m: chi,
        base_alt_m: base,
        trunk,
        polylines,
    }
}

/// ADR Finding 122-B -- a priority flood seeded on the OPEN-OCEAN cells only (`water_class` 1):
/// every other cell rises to the lowest col on its way to the ocean, so an inland below-sea basin
/// and a closed land depression both read their spill level. Not the drainage chain's fill, which
/// takes every `h <= sea` cell as base level and never sees an inland basin (Finding 85).
fn ocean_flood(bf: &GridF32) -> Vec<f32> {
    use crate::lakes::connectivity::{WATER_CLASS_OCEAN, water_class};
    let (w, h) = (bf.width, bf.height);
    let n = w * h;
    let wc = water_class(bf, C1_SEA_LEVEL_NORM);
    let mut f = vec![f32::INFINITY; n];
    let mut heap = BinaryHeap::new();
    for k in 0..n {
        if wc[k] == WATER_CLASS_OCEAN {
            f[k] = bf.data[k];
            heap.push(Item(bf.data[k], k));
        }
    }
    while let Some(Item(v, c)) = heap.pop() {
        if v > f[c] {
            continue;
        }
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let nb = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                let nv = bf.data[nb].max(v);
                if nv < f[nb] {
                    f[nb] = nv;
                    heap.push(Item(nv, nb));
                }
            }
        }
    }
    // a grid with no open ocean (a test field): nothing is closed
    for (v, &b) in f.iter_mut().zip(bf.data.iter()) {
        if !v.is_finite() {
            *v = b;
        }
    }
    f
}

/// What the carve did, per cell.
pub struct CarveMasks {
    /// The cell was lowered by the construction.
    pub carved: Vec<bool>,
    /// The cell lies on a valley FLOOR (within `W/2` of the skeleton) and was lowered.
    pub floor: Vec<bool>,
}

#[derive(PartialEq)]
struct Item(f32, usize);
impl Eq for Item {}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Item {
    // min-heap on the key (distance), deterministic tie-break on the cell index
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.total_cmp(&self.0).then_with(|| o.1.cmp(&self.1))
    }
}

/// Lay the valleys: `out = min(field, V)`, with `V(x) = z_floor(p) + max(0, |x − p| − W(p)/2)·tan θ`
/// and `p` the skeleton sample NEAREST to `x`. **It can only lower a cell** (asserted in the tests).
///
/// ⚠️ **Nearest, not lowest.** Taking the minimum of the cones over every sample lets a DOWNSTREAM
/// sample, whose floor is lower, flatten the floor over half a valley width upstream of it — the
/// long profile would no longer be the χ law. So the cross-section is laid from the nearest sample
/// (vector propagation of the nearest source by Euclidean distance), and the minimum of two cones
/// is taken only between samples of DIFFERENT polylines, where it keeps the surface continuous
/// across the line where two valleys meet. Propagation stops where a wall reaches the terrain, so
/// a valley never cuts through a ridge into the next basin.
pub fn carve(
    field: &GridF32,
    sk: &Skeleton,
    vc: &ValleyConstruction,
    ss: &SteinSteinParams,
) -> (GridF32, CarveMasks) {
    let (w, h) = (sk.width, sk.height);
    assert_eq!((field.width, field.height), (w, h), "skeleton and field sizes differ");
    let n = w * h;
    let tan = vc.wall_deg.to_radians().tan();
    let zfield: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
    let mut src: Vec<(f32, f32, f32, f32)> = Vec::new();
    let mut line_of: Vec<u32> = Vec::new();
    for (li, l) in sk.polylines.iter().enumerate() {
        for &p in l {
            src.push(p);
            line_of.push(li as u32);
        }
    }
    // (distance m, cone value m, on the floor?)
    let geo = |s: usize, c: usize| -> (f32, f32, bool) {
        let (sx, sy, zf, hw) = src[s];
        let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
        let mut dx = (cx - sx).rem_euclid(w as f32);
        if dx > w as f32 / 2.0 {
            dx -= w as f32;
        }
        let mut dy = (cy - sy).rem_euclid(h as f32);
        if dy > h as f32 / 2.0 {
            dy -= h as f32;
        }
        let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
        (d, zf + (d - hw).max(0.0) * tan, d <= hw)
    };
    let mut bestd = vec![f32::INFINITY; n];
    let mut who = vec![u32::MAX; n];
    let mut heap = BinaryHeap::new();
    for (i, sp) in src.iter().enumerate() {
        let cx = (sp.0.floor() as i64).rem_euclid(w as i64) as usize;
        let cy = (sp.1.floor() as i64).rem_euclid(h as i64) as usize;
        let c = cy * w + cx;
        let (d, v, _) = geo(i, c);
        if d < bestd[c] && v < zfield[c] {
            bestd[c] = d;
            who[c] = i as u32;
            heap.push(Item(d, c));
        }
    }
    while let Some(Item(d, c)) = heap.pop() {
        if d > bestd[c] {
            continue;
        }
        let s = who[c] as usize;
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let nb = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                let (dd, vv, _) = geo(s, nb);
                if dd < bestd[nb] && vv < zfield[nb] {
                    bestd[nb] = dd;
                    who[nb] = s as u32;
                    heap.push(Item(dd, nb));
                }
            }
        }
    }
    let mut out = field.clone();
    let mut carved = vec![false; n];
    let mut floor = vec![false; n];
    for c in 0..n {
        if who[c] == u32::MAX {
            continue;
        }
        let me = who[c] as usize;
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        let (_, mut v, mut on_floor) = geo(me, c);
        // continuity across the line where two DIFFERENT valleys meet
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nb = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                let o = who[nb];
                if o != u32::MAX && line_of[o as usize] != line_of[me] {
                    let (_, vv, f) = geo(o as usize, c);
                    if vv < v {
                        v = vv;
                        on_floor = f;
                    }
                }
            }
        }
        if v < zfield[c] {
            let nv = c1_metres_to_altitude_norm(v, ss);
            if nv < out.data[c] {
                out.data[c] = nv;
                carved[c] = true;
                floor[c] = on_floor;
            }
        }
    }
    (out, CarveMasks { carved, floor })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 96×96 V-shaped basin draining north into a sea band: the flow converges on the centre
    /// column, which becomes the trunk.
    fn basin(ss: &SteinSteinParams) -> GridF32 {
        let (w, h) = (96usize, 96usize);
        let mut d = vec![0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let m = if y < 6 {
                    -50.0
                } else {
                    40.0 + 6.0 * (y as f32 - 6.0) + 9.0 * (x as f32 - 48.0).abs()
                };
                d[y * w + x] = c1_metres_to_altitude_norm(m, ss);
            }
        }
        GridF32 { width: w, height: h, data: d }
    }

    /// ADR Finding 122-B, rule 13 — with `basin_base`, no constructed floor upstream of an inland
    /// below-sea basin lies below that basin's spill level; without it (Finding 121), the valleys
    /// are built down toward sea + 0.5 m inside the basin's rim (the negative control).
    #[test]
    fn basin_base_keeps_the_floors_above_the_water_of_a_closed_basin() {
        let ss = SteinSteinParams::default();
        let mut f = basin(&ss);
        let w = f.width;
        // an inland basin on the trunk: a paraboloid dipping ~100 m below the sea, radius 8 cells
        for y in 0..f.height {
            for x in 0..w {
                let d2 = ((x as f32 - 48.0).powi(2) + (y as f32 - 50.0).powi(2)) / 64.0;
                if d2 < 1.0 {
                    let k = y * w + x;
                    let m = c1_altitude_norm_to_metres(f.data[k], &ss) - 400.0 * (1.0 - d2);
                    f.data[k] = c1_metres_to_altitude_norm(m, &ss);
                }
            }
        }
        let run = |basin: bool| {
            let mut vc = ValleyConstruction::f121(0.02, None);
            vc.a_min_km2 = 0.5;
            vc.smooth_m = 0.0;
            vc.basin_base = basin;
            let sk = skeleton(&f, &vc, &ss, 4.8);
            carve(&f, &sk, &vc, &ss)
        };
        // the lowest carved altitude upstream of the basin (beyond its rim)
        let low_up = |out: &GridF32, m: &CarveMasks| -> f32 {
            (60 * w..f.data.len())
                .filter(|&k| m.carved[k])
                .map(|k| c1_altitude_norm_to_metres(out.data[k], &ss))
                .fold(f32::INFINITY, f32::min)
        };
        let (o_off, m_off) = run(false);
        let (o_on, m_on) = run(true);
        let (lo_off, lo_on) = (low_up(&o_off, &m_off), low_up(&o_on, &m_on));
        // the basin's rim on the trunk: 40 + 6·(42 − 6) m ≈ 256 m, its lowest col
        assert!(
            lo_off < 150.0,
            "negative control: F121 must build below the basin's col, got {lo_off}"
        );
        assert!(lo_on > 240.0, "basin_base must keep the floors above the col, got {lo_on}");
        for k in 0..f.data.len() {
            assert!(o_on.data[k] <= f.data[k], "cell {k} RAISED with basin_base");
        }
    }

    /// ADR Finding 121, rule 13 — the construction only LOWERS, lays the floor the law states on
    /// the trunk, and leaves the field untouched beyond the walls.
    #[test]
    fn the_valley_construction_lowers_only_and_lays_its_floor() {
        let ss = SteinSteinParams::default();
        let f = basin(&ss);
        // 96 cells over 4.8 km: 50 m cells, 0.0025 km² per cell
        let domain_km = 4.8;
        let mut vc = ValleyConstruction::f121(0.02, None);
        vc.a_min_km2 = 0.5; // 200 cells: the centre column carries far more
        vc.smooth_m = 0.0; // a straight trunk: no smoothing needed, exact geometry
        let sk = skeleton(&f, &vc, &ss, domain_km);
        let n_trunk = sk.trunk.iter().filter(|&&b| b).count();
        assert!(n_trunk > 20, "the basin must carry a trunk, got {n_trunk} cells");
        let (out, masks) = carve(&f, &sk, &vc, &ss);
        // 1. never raises
        for k in 0..f.data.len() {
            assert!(out.data[k] <= f.data[k], "cell {k} was RAISED");
        }
        assert!(masks.carved.iter().any(|&b| b), "negative control: the construction must bite");
        // 2. on the trunk, where the law's floor is below the field, the cell IS the floor
        let mut checked = 0;
        for k in 0..f.data.len() {
            if sk.trunk[k] {
                let fl = sk.floor_m(k, vc.age_k);
                let z = c1_altitude_norm_to_metres(f.data[k], &ss);
                if fl < z - 1.0 {
                    let got = c1_altitude_norm_to_metres(out.data[k], &ss);
                    assert!(
                        (got - fl).abs() < 0.05,
                        "trunk cell {k}: {got} m against floor {fl} m"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 10, "the floor clause must be exercised, got {checked}");
        // 3. far from the trunk (the basin's edges), the field is untouched
        let w = f.width;
        for y in 10..90 {
            for x in [0usize, 1, 94, 95] {
                assert_eq!(out.data[y * w + x], f.data[y * w + x], "({x},{y}) was touched");
            }
        }
        // 4. the walls: no carved cell is steeper than the wall law allows toward the floor
        let tan = vc.wall_deg.to_radians().tan();
        for y in 20..80 {
            for x in 1..w - 1 {
                let k = y * w + x;
                if masks.carved[k] && !masks.floor[k] {
                    let (a, b) = (
                        c1_altitude_norm_to_metres(out.data[k - 1], &ss),
                        c1_altitude_norm_to_metres(out.data[k + 1], &ss),
                    );
                    let g = (a - b).abs() / (2.0 * sk.cell_m);
                    assert!(g <= tan * 1.05 + 1e-3, "({x},{y}) wall gradient {g} > tan 28°");
                }
            }
        }
    }
}
