//! ADR Finding 165-A/B -- **the continent's physiography, measured** (real metric: z in m, km cells):
//!
//! - [`classify`]: plain / plateau / hill / mountain. Mountain is Kapos et al. 2000 (UNEP-WCMC; the GMBA definition):
//!   z ≥ 2 500 m, or 1 500–2 500 m and slope ≥ 2°, or 1 000–1 500 m and (slope ≥ 5° or R7 ≥ 300 m), or 300–1 000 m and
//!   R7 ≥ 300 m, R7 the elevation range within 7 km. Then plateau (z ≥ [`PLATEAU_Z_M`]), hill (R7 ≥ [`HILL_R7_M`],
//!   below the plateau altitude), plain (F165-B2's one adjustment on Europe: hills are low). The sea is z ≤ 0; the relief and the slope read the sea surface as 0 m.
//! - [`macro_stats`]: the classes' fractions, the land deciles, the mountain components and their width over the
//!   continent's size (the author's « macro-filiform » test).
//! - [`shape`]: components above a threshold, their width (2 × the median distance transform), elongation and the
//!   share near a mask (A2).
//! - [`coast`]: the coastline's box-counting dimension, its length ratio, its aligned share and the 64² block index
//!   (B4).
//!
//! Declared in `docs/reports/relief_method/f165_continent/f165_declared.md` (A2, B2–B4).

use rayon::prelude::*;

/// B2: the hill threshold on R7 (m) and the plateau altitude (m), declared from Meybeck et al. 2001's class logic.
pub const HILL_R7_M: f32 = 150.0;
pub const PLATEAU_Z_M: f32 = 500.0;
/// B2: R7's radius (km), Kapos et al. 2000.
pub const RELIEF_RADIUS_KM: f32 = 7.0;

pub const SEA: u8 = 0;
pub const PLAIN: u8 = 1;
pub const PLATEAU: u8 = 2;
pub const HILL: u8 = 3;
pub const MOUNTAIN: u8 = 4;
pub const CLASS_NAMES: [&str; 5] = ["sea", "plain", "plateau", "hill", "mountain"];

/// The elevation range within `radius_km` (the sea read as 0 m), per cell.
pub fn local_relief(z: &[f32], w: usize, h: usize, cell_km: f32, radius_km: f32) -> Vec<f32> {
    let r = radius_km / cell_km;
    let ri = r.floor() as i64;
    let offs: Vec<(i64, i64)> = (-ri..=ri)
        .flat_map(|dy| (-ri..=ri).map(move |dx| (dx, dy)))
        .filter(|&(dx, dy)| ((dx * dx + dy * dy) as f32).sqrt() <= r)
        .collect();
    (0..w * h)
        .into_par_iter()
        .map(|k| {
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            let (mut lo, mut hi) = (f32::MAX, f32::MIN);
            for &(dx, dy) in &offs {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                let v = z[ny as usize * w + nx as usize].max(0.0);
                lo = lo.min(v);
                hi = hi.max(v);
            }
            hi - lo
        })
        .collect()
}

/// The slope (degrees, central differences, the sea read as 0 m).
pub fn slope_deg(z: &[f32], w: usize, h: usize, cell_km: f32) -> Vec<f32> {
    let c = cell_km * 1000.0;
    (0..w * h)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % w, k / w);
            let at = |x: usize, y: usize| z[y * w + x].max(0.0);
            let gx = (at((x + 1).min(w - 1), y) - at(x.saturating_sub(1), y)) / (2.0 * c);
            let gy = (at(x, (y + 1).min(h - 1)) - at(x, y.saturating_sub(1))) / (2.0 * c);
            (gx * gx + gy * gy).sqrt().atan().to_degrees()
        })
        .collect()
}

/// The macro classes (B2).
pub fn classify(z: &[f32], w: usize, h: usize, cell_km: f32) -> Vec<u8> {
    let r7 = local_relief(z, w, h, cell_km, RELIEF_RADIUS_KM);
    let s = slope_deg(z, w, h, cell_km);
    (0..w * h)
        .map(|k| {
            let (v, r, sl) = (z[k], r7[k], s[k]);
            if v <= 0.0 {
                return SEA;
            }
            let mountain = v >= 2500.0
                || (v >= 1500.0 && sl >= 2.0)
                || (v >= 1000.0 && (sl >= 5.0 || r >= 300.0))
                || (v >= 300.0 && r >= 300.0);
            // F165-B2, the one adjustment on Europe (before any world of ours): hills are LOW (Meybeck et al. 2001's
            // hills sit below its plateaus), so a non-mountain cell at or above the plateau altitude is a plateau
            if mountain {
                MOUNTAIN
            } else if v >= PLATEAU_Z_M {
                PLATEAU
            } else if r >= HILL_R7_M {
                HILL
            } else {
                PLAIN
            }
        })
        .collect()
}

/// The exact squared Euclidean distance (cells) from each `true` cell to the nearest `false` cell (or the outside),
/// Felzenszwalb & Huttenlocher 2012, separable.
pub fn edt(mask: &[bool], w: usize, h: usize) -> Vec<f32> {
    const INF: f64 = 1e18;
    fn pass(f: &[f64]) -> Vec<f64> {
        let n = f.len();
        let mut d = vec![0f64; n];
        let mut v = vec![0usize; n];
        let mut zz = vec![0f64; n + 1];
        let mut k = 0usize;
        v[0] = 0;
        zz[0] = f64::NEG_INFINITY;
        zz[1] = f64::INFINITY;
        let sep = |q: usize, p: usize| ((f[q] + (q * q) as f64) - (f[p] + (p * p) as f64)) / (2.0 * q as f64 - 2.0 * p as f64);
        for q in 1..n {
            let mut s = sep(q, v[k]);
            while s <= zz[k] {
                k -= 1;
                s = sep(q, v[k]);
            }
            k += 1;
            v[k] = q;
            zz[k] = s;
            zz[k + 1] = f64::INFINITY;
        }
        k = 0;
        for q in 0..n {
            while zz[k + 1] < q as f64 {
                k += 1;
            }
            let p = v[k];
            d[q] = (q as f64 - p as f64).powi(2) + f[p];
        }
        d
    }
    // pad by one cell of `false` so the outside counts as background
    let (pw, ph) = (w + 2, h + 2);
    let mut g = vec![0f64; pw * ph];
    for y in 0..h {
        for x in 0..w {
            if mask[y * w + x] {
                g[(y + 1) * pw + x + 1] = INF;
            }
        }
    }
    let cols: Vec<Vec<f64>> = (0..pw).into_par_iter().map(|x| pass(&(0..ph).map(|y| g[y * pw + x]).collect::<Vec<_>>())).collect();
    for x in 0..pw {
        for y in 0..ph {
            g[y * pw + x] = cols[x][y];
        }
    }
    let rows: Vec<Vec<f64>> = (0..ph).into_par_iter().map(|y| pass(&g[y * pw..(y + 1) * pw])).collect();
    let mut out = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            out[y * w + x] = rows[y + 1][x + 1].sqrt() as f32;
        }
    }
    out
}

/// 8-connected components of `mask` (labels from 1; 0 = none) and their sizes (index = label − 1).
pub fn components(mask: &[bool], w: usize, h: usize) -> (Vec<u32>, Vec<usize>) {
    let mut lab = vec![0u32; w * h];
    let mut sizes = Vec::new();
    let mut stack = Vec::new();
    for s in 0..w * h {
        if !mask[s] || lab[s] != 0 {
            continue;
        }
        let id = sizes.len() as u32 + 1;
        lab[s] = id;
        stack.push(s);
        let mut c = 0usize;
        while let Some(k) = stack.pop() {
            c += 1;
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            for dy in -1..=1i64 {
                for dx in -1..=1i64 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                        continue;
                    }
                    let j = ny as usize * w + nx as usize;
                    if mask[j] && lab[j] == 0 {
                        lab[j] = id;
                        stack.push(j);
                    }
                }
            }
        }
        sizes.push(c);
    }
    (lab, sizes)
}

/// Chebyshev distance (cells) from the `src` cells, capped at `cap` (`u32::MAX` beyond).
pub fn chebyshev(src: &[bool], w: usize, h: usize, cap: u32) -> Vec<u32> {
    let mut d = vec![u32::MAX; w * h];
    let mut front: Vec<usize> = (0..w * h).filter(|&k| src[k]).collect();
    for &k in &front {
        d[k] = 0;
    }
    let mut lvl = 0;
    while !front.is_empty() && lvl < cap {
        lvl += 1;
        let mut next = Vec::new();
        for &k in &front {
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            for dy in -1..=1i64 {
                for dx in -1..=1i64 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                        continue;
                    }
                    let j = ny as usize * w + nx as usize;
                    if d[j] == u32::MAX {
                        d[j] = lvl;
                        next.push(j);
                    }
                }
            }
        }
        front = next;
    }
    d
}

/// A2: the shape of a set of cells.
#[derive(Clone, Debug, Default)]
pub struct Shape {
    pub cells: usize,
    pub components: usize,
    /// 2 × the median distance transform over the cells (cells).
    pub width_cells: f32,
    /// √(λ₁/λ₂) of each component's coordinate covariance, area-weighted (components of ≥ 4 cells).
    pub elongation: f32,
    /// The share of the cells within `near_cells` (Chebyshev) of the `near` mask.
    pub near_share: f32,
}

/// A2's measures on `mask`.
pub fn shape(mask: &[bool], w: usize, h: usize, near: Option<&[bool]>, near_cells: u32) -> Shape {
    let cells = mask.iter().filter(|&&m| m).count();
    if cells == 0 {
        return Shape { width_cells: f32::NAN, elongation: f32::NAN, near_share: f32::NAN, ..Default::default() };
    }
    let d = edt(mask, w, h);
    let mut v: Vec<f32> = (0..w * h).filter(|&k| mask[k]).map(|k| d[k]).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let width_cells = 2.0 * v[v.len() / 2];
    let (lab, sizes) = components(mask, w, h);
    let mut acc = vec![[0f64; 5]; sizes.len()];
    for k in 0..w * h {
        if lab[k] != 0 {
            let (x, y) = ((k % w) as f64, (k / w) as f64);
            let a = &mut acc[lab[k] as usize - 1];
            a[0] += x;
            a[1] += y;
            a[2] += x * x;
            a[3] += y * y;
            a[4] += x * y;
        }
    }
    let (mut ew, mut es) = (0f64, 0f64);
    for (i, a) in acc.iter().enumerate() {
        let n = sizes[i] as f64;
        if sizes[i] < 4 {
            continue;
        }
        let (mx, my) = (a[0] / n, a[1] / n);
        let (sxx, syy, sxy) = (a[2] / n - mx * mx, a[3] / n - my * my, a[4] / n - mx * my);
        let tr = sxx + syy;
        let det = sxx * syy - sxy * sxy;
        let disc = ((tr * tr / 4.0) - det).max(0.0).sqrt();
        let (l1, l2) = (tr / 2.0 + disc, (tr / 2.0 - disc).max(1e-9));
        es += n * (l1 / l2).sqrt();
        ew += n;
    }
    let near_share = near.map_or(f32::NAN, |m| {
        let dd = chebyshev(m, w, h, near_cells + 1);
        (0..w * h).filter(|&k| mask[k] && dd[k] <= near_cells).count() as f32 / cells as f32
    });
    Shape {
        cells,
        components: sizes.len(),
        width_cells,
        elongation: if ew > 0.0 { (es / ew) as f32 } else { f32::NAN },
        near_share,
    }
}

/// B3: the macro quantities of a classified field.
#[derive(Clone, Debug, Default)]
pub struct MacroStats {
    /// plain, plateau, hill, mountain (share of the land).
    pub fractions: [f32; 4],
    pub deciles: Vec<f32>,
    pub land_km2: f64,
    /// Mountain components of ≥ 4 cells.
    pub mtn_components: usize,
    /// 2 × the median distance transform over the mountain cells (km).
    pub mtn_width_km: f32,
    /// `mtn_width_km` / √(land km²).
    pub mtn_ratio: f32,
}

/// B3 on a field (m) and its classes.
pub fn macro_stats(z: &[f32], classes: &[u8], w: usize, h: usize, cell_km: f32) -> MacroStats {
    let land = classes.iter().filter(|&&c| c != SEA).count();
    let mut out = MacroStats { land_km2: land as f64 * (cell_km * cell_km) as f64, ..Default::default() };
    for (i, c) in [PLAIN, PLATEAU, HILL, MOUNTAIN].iter().enumerate() {
        out.fractions[i] = classes.iter().filter(|&&x| x == *c).count() as f32 / land.max(1) as f32;
    }
    let mut alts: Vec<f32> = (0..w * h).filter(|&k| classes[k] != SEA).map(|k| z[k]).collect();
    alts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    out.deciles = (1..10).map(|d| if alts.is_empty() { f32::NAN } else { alts[(alts.len() - 1) * d / 10] }).collect();
    let m: Vec<bool> = classes.iter().map(|&c| c == MOUNTAIN).collect();
    let (_, sizes) = components(&m, w, h);
    out.mtn_components = sizes.iter().filter(|&&s| s >= 4).count();
    let s = shape(&m, w, h, None, 0);
    out.mtn_width_km = s.width_cells * cell_km;
    out.mtn_ratio = out.mtn_width_km / (out.land_km2.sqrt() as f32).max(1e-9);
    out
}

/// B4: the coastline's grain.
#[derive(Clone, Debug, Default)]
pub struct Coast {
    /// The box-counting dimension over 2–64 cells.
    pub dimension: f32,
    /// L(2 cells) / L(32 cells), L = N·ε.
    pub length_ratio: f32,
    /// The aligned share (±5° of 0/45/90/135°) of the coast's 8-step chords.
    pub aligned: f32,
    /// The share of land/sea edges on the `block` lattice over 1 / `block`.
    pub block_index: f32,
    pub edges: usize,
}

/// B4 on a land mask (`true` = land).
pub fn coast(land: &[bool], w: usize, h: usize, block: usize) -> Coast {
    let is_land = |x: i64, y: i64| x >= 0 && y >= 0 && x < w as i64 && y < h as i64 && land[y as usize * w + x as usize];
    // the coastline cells: land with a 4-neighbour sea (the outside is sea)
    let cl: Vec<bool> = (0..w * h)
        .map(|k| {
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            land[k] && (!is_land(x - 1, y) || !is_land(x + 1, y) || !is_land(x, y - 1) || !is_land(x, y + 1))
        })
        .collect();
    let mut pts = Vec::new();
    for s in [2usize, 4, 8, 16, 32, 64] {
        let (bw, bh) = (w.div_ceil(s), h.div_ceil(s));
        let mut hit = vec![false; bw * bh];
        for k in 0..w * h {
            if cl[k] {
                hit[(k / w / s) * bw + (k % w) / s] = true;
            }
        }
        pts.push((s as f64, hit.iter().filter(|&&b| b).count() as f64));
    }
    let lx: Vec<f64> = pts.iter().map(|p| (1.0 / p.0).ln()).collect();
    let ly: Vec<f64> = pts.iter().map(|p| p.1.max(1.0).ln()).collect();
    let (mx, my) = (lx.iter().sum::<f64>() / 6.0, ly.iter().sum::<f64>() / 6.0);
    let num: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let den: f64 = lx.iter().map(|a| (a - mx).powi(2)).sum();
    let length_ratio = (pts[0].1 * pts[0].0) / (pts[4].1 * pts[4].0).max(1.0);
    // the crack contour: directed land/sea edges, land on the left (counter-clockwise), followed into loops
    use std::collections::HashMap;
    let mut out_edges: HashMap<(i64, i64), Vec<(i64, i64)>> = HashMap::new();
    let (mut edges, mut on_lattice) = (0usize, 0usize);
    let b = block as i64;
    for y in 0..h as i64 {
        for x in 0..w as i64 {
            if !is_land(x, y) {
                continue;
            }
            let mut add = |a: (i64, i64), c: (i64, i64), lattice: bool| {
                out_edges.entry(a).or_default().push(c);
                edges += 1;
                if lattice {
                    on_lattice += 1;
                }
            };
            if !is_land(x, y - 1) {
                add((x, y), (x + 1, y), y % b == 0);
            }
            if !is_land(x + 1, y) {
                add((x + 1, y), (x + 1, y + 1), (x + 1) % b == 0);
            }
            if !is_land(x, y + 1) {
                add((x + 1, y + 1), (x, y + 1), (y + 1) % b == 0);
            }
            if !is_land(x - 1, y) {
                add((x, y + 1), (x, y), x % b == 0);
            }
        }
    }
    let (mut all, mut al) = (0f64, 0f64);
    let mut starts: Vec<(i64, i64)> = out_edges.keys().copied().collect();
    starts.sort();
    for s in starts {
        while let Some(first) = out_edges.get_mut(&s).and_then(|v| v.pop()) {
            let mut path = vec![s, first];
            let mut cur = first;
            let mut prev_dir = (first.0 - s.0, first.1 - s.1);
            while cur != s {
                let Some(v) = out_edges.get_mut(&cur) else { break };
                if v.is_empty() {
                    break;
                }
                // at a saddle, turn left first (keeps a land cell on the left)
                let left = (-prev_dir.1, prev_dir.0);
                let pick = v.iter().position(|&n| (n.0 - cur.0, n.1 - cur.1) == left).unwrap_or(0);
                let nxt = v.swap_remove(pick);
                prev_dir = (nxt.0 - cur.0, nxt.1 - cur.1);
                path.push(nxt);
                cur = nxt;
            }
            if path.len() < 33 {
                continue;
            }
            let mut i = 0;
            while i + 8 < path.len() {
                let (dx, dy) = ((path[i + 8].0 - path[i].0) as f32, (path[i + 8].1 - path[i].1) as f32);
                let len = (dx * dx + dy * dy).sqrt() as f64;
                if len > 0.0 {
                    let phi = dy.abs().atan2(dx.abs()).to_degrees() % 45.0;
                    all += len;
                    if phi.min(45.0 - phi) <= 5.0 {
                        al += len;
                    }
                }
                i += 8;
            }
        }
    }
    Coast {
        dimension: (num / den.max(1e-12)) as f32,
        length_ratio: length_ratio as f32,
        aligned: if all > 0.0 { (al / all) as f32 } else { f32::NAN },
        block_index: if edges > 0 { (on_lattice as f32 / edges as f32) * block as f32 } else { f32::NAN },
        edges,
    }
}

/// ADR Finding 167-C -- Hovius's spacing ratio on the mountain fronts (Hovius 1996, via Talling et al. 1997 eq. 1,
/// R = W/S), read per transverse basin (amended at its control): an outlet is a cell of the `mountain` class whose D8
/// receiver lies outside it; each mountain cell belongs to the outlet it drains to; for each outlet whose basin (inside
/// the class) covers ≥ `a_min_km2`, L = the largest distance from the outlet to a cell of its basin (the half-width W for
/// a transverse basin) and S = A / L (its width along the front, the outlet spacing); R = the mean of L / S = L² / A.
#[derive(Clone, Debug, Default)]
pub struct Hovius {
    pub ratio: f32,
    pub outlets: usize,
    /// The means of S and W (km).
    pub spacing_km: f32,
    pub half_width_km: f32,
}

/// [`Hovius`] on a field in metres and its mountain mask.
pub fn hovius(z: &crate::grid::GridF32, mountain: &[bool], cell_km: f32, a_min_km2: f32) -> Hovius {
    use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let n = z.width;
    let hy = super::hydro::hydrology(z, cell_km);
    let recv = |k: usize| -> Option<usize> {
        let d = hy.dir[k];
        if d == DIR_NONE {
            return None;
        }
        let (x, y) = ((k % n) as i64 + D8_DX[d as usize] as i64, (k / n) as i64 + D8_DY[d as usize] as i64);
        (x >= 0 && y >= 0 && x < n as i64 && y < n as i64).then(|| y as usize * n + x as usize)
    };
    // the outlet each mountain cell drains to (u32::MAX: none), memoised along the D8 paths
    const NONE: u32 = u32::MAX;
    const TODO: u32 = u32::MAX - 1;
    let mut lab = vec![TODO; n * n];
    for k in 0..n * n {
        if !mountain[k] {
            lab[k] = NONE;
        }
    }
    let mut path = Vec::new();
    for s in 0..n * n {
        if lab[s] != TODO {
            continue;
        }
        path.clear();
        let mut cur = s;
        let res = loop {
            if lab[cur] != TODO {
                break lab[cur];
            }
            path.push(cur);
            match recv(cur) {
                Some(r) if mountain[r] => cur = r,
                Some(_) => break cur as u32,
                None => break NONE,
            }
        };
        for &c in &path {
            lab[c] = res;
        }
    }
    let mut basins: std::collections::HashMap<u32, (usize, f32)> = std::collections::HashMap::new();
    for k in 0..n * n {
        let o = lab[k];
        if o == NONE || o == TODO {
            continue;
        }
        let (ox, oy) = ((o as usize % n) as f32, (o as usize / n) as f32);
        let d = (((k % n) as f32 - ox).powi(2) + ((k / n) as f32 - oy).powi(2)).sqrt() + 0.5;
        let e = basins.entry(o).or_insert((0, 0.0));
        e.0 += 1;
        e.1 = e.1.max(d);
    }
    let c2 = cell_km * cell_km;
    let rows: Vec<(f32, f32)> = basins
        .values()
        .filter(|(cells, _)| *cells as f32 * c2 >= a_min_km2)
        .map(|(cells, l)| (*cells as f32 * c2, l * cell_km))
        .collect();
    if rows.is_empty() {
        return Hovius { ratio: f32::NAN, ..Default::default() };
    }
    let m = rows.len() as f32;
    Hovius {
        ratio: rows.iter().map(|(a, l)| l * l / a).sum::<f32>() / m,
        outlets: rows.len(),
        spacing_km: rows.iter().map(|(a, l)| a / l).sum::<f32>() / m,
        half_width_km: rows.iter().map(|(_, l)| *l).sum::<f32>() / m,
    }
}

/// ADR Finding 168-A3 -- the periodic 8-connected components of `mask` on an n × n torus: per cell its component
/// (`u32::MAX` off the mask), per component (cells, wraps). A component wraps when a BFS that carries the unwrapped
/// coordinates reaches one of its cells at two different unwrapped positions: a cycle of non-zero winding number.
pub fn periodic_components(mask: &[bool], n: usize) -> (Vec<u32>, Vec<(usize, bool)>) {
    let mut id = vec![u32::MAX; n * n];
    let (mut ux, mut uy) = (vec![0i32; n * n], vec![0i32; n * n]);
    let mut comps = Vec::new();
    let mut stack = Vec::new();
    let ni = n as i32;
    for s in 0..n * n {
        if !mask[s] || id[s] != u32::MAX {
            continue;
        }
        let c = comps.len() as u32;
        let (mut cells, mut wraps) = (0usize, false);
        id[s] = c;
        ux[s] = (s % n) as i32;
        uy[s] = (s / n) as i32;
        stack.push(s);
        while let Some(k) = stack.pop() {
            cells += 1;
            for dy in -1..=1i32 {
                for dx in -1..=1i32 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let (x, y) = (ux[k] + dx, uy[k] + dy);
                    let m = y.rem_euclid(ni) as usize * n + x.rem_euclid(ni) as usize;
                    if !mask[m] {
                        continue;
                    }
                    if id[m] == u32::MAX {
                        id[m] = c;
                        ux[m] = x;
                        uy[m] = y;
                        stack.push(m);
                    } else if ux[m] != x || uy[m] != y {
                        wraps = true;
                    }
                }
            }
        }
        comps.push((cells, wraps));
    }
    (id, comps)
}

/// F168-A3: the periodic Chebyshev dilation of `mask` by `r` cells (separable circular windows of 2r + 1).
pub fn dilate_periodic(mask: &[bool], n: usize, r: usize) -> Vec<bool> {
    if 2 * r + 1 >= n {
        return vec![mask.iter().any(|&b| b); n * n];
    }
    let pass = |src: &[bool], stride: usize, step: usize| -> Vec<bool> {
        // stride: between lines; step: along a line
        let mut out = vec![false; n * n];
        for line in 0..n {
            let at = |i: usize| src[line * stride + (i % n) * step];
            let mut count = (0..=2 * r).filter(|&i| at(i + n - r)).count();
            for i in 0..n {
                out[line * stride + i * step] = count > 0;
                count -= at(i + n - r) as usize;
                count += at(i + r + 1) as usize;
            }
        }
        out
    };
    let rows = pass(mask, n, 1);
    pass(&rows, 1, n)
}

/// F168-A3: whether the largest land mass can be sailed around on the torus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circumnavigation {
    /// The largest land mass (periodic 8-connected), in cells.
    pub main_cells: usize,
    /// It wraps the torus (a non-zero winding number on either axis): not circumnavigable.
    pub wraps: bool,
    /// r*: the largest Chebyshev dilation of ALL the land after which the main mass's component still does not wrap;
    /// a sea loop at least 2r* + 1 cells wide surrounds it. `None` when it wraps or there is no land.
    pub sea_r: Option<usize>,
}

/// F168-A3, the circumnavigation test (declared before the measure): (1) the largest land mass does not wrap; (2) r*,
/// by bisection (the dilation is monotone, so the main component's wrapping is too).
pub fn circumnavigation(land: &[bool], n: usize) -> Circumnavigation {
    let (id, comps) = periodic_components(land, n);
    let Some(main) = (0..comps.len()).max_by(|a, b| comps[*a].0.cmp(&comps[*b].0).then(b.cmp(a))) else {
        return Circumnavigation { main_cells: 0, wraps: false, sea_r: None };
    };
    let (main_cells, wraps) = comps[main];
    if wraps {
        return Circumnavigation { main_cells, wraps, sea_r: None };
    }
    let rep = id.iter().position(|&c| c == main as u32).unwrap();
    let wraps_at = |r: usize| {
        let (id2, c2) = periodic_components(&dilate_periodic(land, n, r), n);
        c2[id2[rep] as usize].1
    };
    let (mut lo, mut hi) = (0usize, n / 2);
    if !wraps_at(hi) {
        return Circumnavigation { main_cells, wraps, sea_r: Some(hi) };
    }
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if wraps_at(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Circumnavigation { main_cells, wraps, sea_r: Some(lo) }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F168-A3: a centred square of side 20 on a 64² torus is circumnavigable with r* = 21 (side 62 still leaves a
    /// 2-cell gap; side 64 closes it); a full band and a diagonal line wrap; a ring around the seam does not.
    #[test]
    fn circumnavigation_reads_islands_bands_and_diagonals() {
        let n = 64;
        let sq: Vec<bool> = (0..n * n).map(|k| (22..42).contains(&(k % n)) && (22..42).contains(&(k / n))).collect();
        assert_eq!(circumnavigation(&sq, n), Circumnavigation { main_cells: 400, wraps: false, sea_r: Some(21) });
        let band: Vec<bool> = (0..n * n).map(|k| (10..20).contains(&(k / n))).collect();
        assert!(circumnavigation(&band, n).wraps);
        let diag: Vec<bool> = (0..n * n).map(|k| k % n == k / n).collect();
        assert!(circumnavigation(&diag, n).wraps);
        // a square straddling both seams (the torus corner) is one island, not a wrap
        let corner: Vec<bool> = (0..n * n).map(|k| ((k % n) + 5) % n < 10 && ((k / n) + 5) % n < 10).collect();
        let c = circumnavigation(&corner, n);
        assert!(!c.wraps && c.main_cells == 100 && c.sea_r == Some(26), "{c:?}");
        // the dilation is the periodic Chebyshev ball
        let one: Vec<bool> = (0..n * n).map(|k| k == 0).collect();
        let d = dilate_periodic(&one, n, 2);
        assert_eq!(d.iter().filter(|&&b| b).count(), 25);
        assert!(d[2 * n + 2] && d[(n - 2) * n + n - 2] && !d[3 * n]);
    }

    /// Hovius's ratio on a synthetic range: a ridge of half-width 24 cells with transverse V valleys every 12 cells,
    /// draining both flanks to the sea, reads R ≈ W/S = 2 (within ±35 %, the front's ends and the D8 staircase).
    #[test]
    fn hovius_reads_a_synthetic_range() {
        let (n, c) = (256usize, 0.2f32);
        let mut z = crate::grid::GridF32::new(n, n, -50.0);
        let mut m = vec![false; n * n];
        for y in 0..n {
            for x in 16..240 {
                let d = (y as f32 - 128.0).abs();
                if d > 24.0 {
                    continue;
                }
                let k = y * n + x;
                let v = ((x % 12) as f32 - 6.0).abs();
                z.data[k] = 10.0 + 40.0 * (24.0 - d) + 30.0 * v;
                m[k] = true;
            }
        }
        let h = hovius(&z, &m, c, 0.5);
        assert!(h.outlets > 20, "{h:?}");
        assert!((h.ratio - 2.0).abs() < 0.7, "{h:?}");
    }

    /// The classifier on synthetic terrain: a flat lowland reads plain, a high flat table plateau, a 2 600 m block
    /// mountain, a rough low terrain hill; the distance transform is exact on a rectangle; a lattice-aligned square
    /// island reads a coastline dimension near 1, all chords aligned and a block index far above 1, a disc does not.
    #[test]
    fn the_classes_the_distance_and_the_coast_read_as_built() {
        let (w, h, c) = (96usize, 96usize, 1.5625f32);
        let mut z = vec![-100f32; w * h];
        for y in 8..88 {
            for x in 8..88 {
                let k = y * w + x;
                z[k] = if x < 30 {
                    100.0
                } else if x < 50 {
                    800.0
                } else if x < 70 {
                    2600.0
                } else {
                    250.0 + 120.0 * ((x as f32 * 1.7).sin() * (y as f32 * 1.3).cos())
                };
            }
        }
        let cl = classify(&z, w, h, c);
        assert_eq!(cl[48 * w + 18], PLAIN);
        assert_eq!(cl[48 * w + 40], PLATEAU);
        assert_eq!(cl[48 * w + 60], MOUNTAIN);
        assert_eq!(cl[48 * w + 80], HILL);
        let mut m = vec![false; 40 * 30];
        for y in 5..25 {
            for x in 5..35 {
                m[y * 40 + x] = true;
            }
        }
        let d = edt(&m, 40, 30);
        assert!((d[15 * 40 + 20] - 10.0).abs() < 1e-4, "{}", d[15 * 40 + 20]);
        assert!((d[5 * 40 + 20] - 1.0).abs() < 1e-4);
        let n = 1024;
        let sq: Vec<bool> = (0..n * n).map(|k| (k % n) >= 128 && (k % n) < 896 && (k / n) >= 128 && (k / n) < 896).collect();
        let cs = coast(&sq, n, n, 32);
        assert!(cs.aligned > 0.95 && cs.block_index > 10.0 && (cs.dimension - 1.0).abs() < 0.15, "{cs:?}");
        let disc: Vec<bool> = (0..n * n).map(|k| (((k % n) as f32 - 511.3).powi(2) + ((k / n) as f32 - 512.6).powi(2)).sqrt() < 380.0).collect();
        let cd = coast(&disc, n, n, 32);
        assert!(cd.block_index < 3.0 && cd.aligned < 0.6, "{cd:?}");
    }
}
