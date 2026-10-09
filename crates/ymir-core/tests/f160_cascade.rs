//! ADR Finding 160 — the cascade, round 2: one physics level (P128, P256), then Schott 2024's amplification transposed
//! (S, S+ρ, S+U × faible / moyen / fort), the final retargeting, measured level by level against the production témoin
//! brought to the same resolutions. Declared in `docs/reports/relief_method/f160_amplification/f160_declared.md`.
//!
//! The témoin (8192², ~2.5 min) is reduced once to its block means and crops, cached under `target_f160/f160_cache`.
//!
//! Run: cargo test -p ymir-core --release --test f160_cascade -- --ignored --exact f160_cascade --nocapture

mod common;

use common::{CANONICAL_ORIGIN, Knobs, PSEED, aniso, build_world, pct, sorted};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, AmpLevel, Budget, Chain, Variant, ocean_mask, restrict};
use ymir_core::cascade::hydro::{Hydro, draw, hydrology};
use ymir_core::cascade::{CascadeConfig, physics_level, roll, shade_m_rgba};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const S_EQ: f32 = 0.024;
const DOMAIN_KM: f32 = 400.0;
const PEAK_TARGET_M: f32 = 2850.0;
const PEAK_OK: (f32, f32) = (0.85 * 2700.0, 1.15 * 3000.0);
const REF_LEVELS: [usize; 5] = [128, 256, 512, 1024, 2048];
/// F159's crops (512² origins): A range, B plateau (the amended rule), C basin + lake 4.
const CROPS: [(&str, usize, usize); 3] = [("A_chaine", 128, 320), ("B_plateau", 176, 192), ("C_bassin_lac4", 74, 160)];

fn cache_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target_f160/f160_cache")
}
fn out_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/reports/relief_method/f160_amplification/images")
}
fn save_grid(p: &Path, g: &GridF32) {
    let mut b = Vec::with_capacity(8 + g.data.len() * 4);
    b.extend_from_slice(&(g.width as u32).to_le_bytes());
    b.extend_from_slice(&(g.height as u32).to_le_bytes());
    for v in &g.data {
        b.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(p, b).unwrap();
}
fn load_grid(p: &Path) -> Option<GridF32> {
    let b = std::fs::read(p).ok()?;
    let w = u32::from_le_bytes(b[0..4].try_into().unwrap()) as usize;
    let h = u32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
    let data = b[8..].chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect();
    Some(GridF32 { width: w, height: h, data })
}
/// F159's D9: centred block means.
fn block_mean(g: &GridF32, n: usize) -> GridF32 {
    let f = g.width / n;
    let w = g.width;
    let mut out = GridF32::new(n, n, 0.0);
    out.data.par_iter_mut().enumerate().for_each(|(k, o)| {
        let (i, j) = (k % n, k / n);
        let mut s = 0f64;
        for dy in 0..f {
            let y = (j * f + w + dy - f / 2) % w;
            for dx in 0..f {
                let x = (i * f + w + dx - f / 2) % w;
                s += g.data[y * w + x] as f64;
            }
        }
        *o = (s / (f * f) as f64) as f32;
    });
    out
}
fn to_m(g: &GridF32, n2m: f32) -> GridF32 {
    GridF32 { width: g.width, height: g.height, data: g.data.iter().map(|v| (v - 0.5) * n2m).collect() }
}

// ───────────────────────────── the instruments (F159's, in metres, plus D7's) ─────────────────────────────

#[derive(Clone, Debug, Default)]
struct Reading {
    land: usize,
    walls28: f64,
    walls28_cells: usize,
    walls33: f64,
    crest_share: f64,
    crest_p50: f32,
    crest_p90: f32,
    facets: f64,
    facet_cells: usize,
    coherent: f32,
    r8: f32,
    r8_windows: usize,
    lambda: f32,
    relief_p50: f32,
    mean: f32,
    p90: f32,
    peak: f32,
    trunks: Vec<(usize, usize)>,
    bands: Vec<f32>,
    lakes: usize,
    lake_km2: f32,
}

fn land_of(z: &GridF32) -> Vec<bool> {
    let oc = ocean_mask(z);
    oc.iter().map(|o| !o).collect()
}

fn crest_mask(zm: &[f32], land: &[bool], n: usize) -> Vec<bool> {
    (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % n, k / n);
            if x == 0 || y == 0 || x + 1 == n || y + 1 == n || !land[k] {
                return false;
            }
            let z = zm[k];
            [(k - 1, k + 1), (k - n, k + n), (k - n - 1, k + n + 1), (k - n + 1, k + n - 1)].iter().any(|&(a, b)| z > zm[a] && z > zm[b])
        })
        .collect()
}
fn laplacian(zm: &[f32], n: usize, cell_km: f32) -> Vec<f32> {
    (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % n, k / n);
            if x == 0 || y == 0 || x + 1 == n || y + 1 == n {
                return 0.0;
            }
            (zm[k - 1] + zm[k + 1] + zm[k - n] + zm[k + n] - 4.0 * zm[k]) / (cell_km * cell_km)
        })
        .collect()
}
fn facets(crest: &[bool], land: &[bool], lap: &[f32], n: usize) -> (f64, usize) {
    let res: Vec<(bool, bool)> = (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % n, k / n);
            if !crest[k] || x < 3 || y < 3 || x + 3 >= n || y + 3 >= n {
                return (false, false);
            }
            let mut ring = Vec::with_capacity(16);
            for dy in -2i32..=2 {
                for dx in -2i32..=2 {
                    if dx.abs().max(dy.abs()) != 2 {
                        continue;
                    }
                    let j = ((y as i32 + dy) as usize) * n + (x as i32 + dx) as usize;
                    if land[j] && !crest[j] {
                        ring.push(lap[j].abs());
                    }
                }
            }
            if ring.is_empty() {
                return (true, false);
            }
            ring.sort_by(|a, b| a.partial_cmp(b).unwrap());
            (true, ring[ring.len() / 2] < 0.25 * lap[k].abs())
        })
        .collect();
    let crests = res.iter().filter(|r| r.0).count();
    let f = res.iter().filter(|r| r.1).count();
    (100.0 * f as f64 / crests.max(1) as f64, f)
}
fn facet_control(n: usize) -> (f64, f64) {
    let cell_km = DOMAIN_KM / n as f32;
    let (p, r) = (24usize, 11.0f32);
    let mut pyr = vec![0f32; n * n];
    let mut par = vec![0f32; n * n];
    for y in 0..n {
        for x in 0..n {
            let (fx, fy) = (x as f32, y as f32);
            let (u, v) = (fx * 0.866 + fy * 0.5, -fx * 0.5 + fy * 0.866);
            let (du, dv) = ((u.rem_euclid(p as f32)) - p as f32 / 2.0, (v.rem_euclid(p as f32)) - p as f32 / 2.0);
            pyr[y * n + x] = 100.0 + 800.0 * (1.0 - du.abs().max(dv.abs()) / r).max(0.0);
            par[y * n + x] = 100.0 + 800.0 * (1.0 - (du * du + dv * dv) / (r * r)).max(0.0);
        }
    }
    let land = vec![true; n * n];
    let f = |zm: &[f32]| {
        let crest = crest_mask(zm, &land, n);
        let lap = laplacian(zm, n, cell_km);
        facets(&crest, &land, &lap, n).0
    };
    (f(&pyr), f(&par))
}
fn valley_spacing(zm: &[f32], land: &[bool], n: usize) -> f32 {
    let lines: Vec<(usize, usize)> = (0..2 * n)
        .into_par_iter()
        .map(|li| {
            let idx = |i: usize| if li < n { li * n + i } else { i * n + (li - n) };
            let (mut len, mut minima) = (0usize, 0usize);
            let mut i = 0;
            while i < n {
                if !land[idx(i)] {
                    i += 1;
                    continue;
                }
                let s = i;
                while i < n && land[idx(i)] {
                    i += 1;
                }
                let run: Vec<f32> = (s..i).map(|t| zm[idx(t)]).collect();
                if run.len() < 16 {
                    continue;
                }
                len += run.len();
                for t in 1..run.len() - 1 {
                    if !(run[t] < run[t - 1] && run[t] <= run[t + 1]) {
                        continue;
                    }
                    let mut lm = run[t];
                    for u in (0..t).rev() {
                        if run[u] < run[t] {
                            break;
                        }
                        lm = lm.max(run[u]);
                    }
                    let mut rm = run[t];
                    for &v in &run[t + 1..] {
                        if v < run[t] {
                            break;
                        }
                        rm = rm.max(v);
                    }
                    if lm.min(rm) - run[t] >= 1.0 {
                        minima += 1;
                    }
                }
            }
            (len, minima)
        })
        .collect();
    let (l, m) = lines.iter().fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
    l as f32 / m.max(1) as f32
}
/// D7.1: the band-pass RMS, bands [2^k, 2^(k+1)] cells, B_r a periodic box blur of side r (the ocean read at 0 m).
fn bands(zm: &[f32], land: &[bool], n: usize) -> Vec<f32> {
    let blur = |src: &[f32], r: usize| -> Vec<f32> {
        if r <= 1 {
            return src.to_vec();
        }
        let h0 = r / 2;
        let pass = |s: &[f32], horiz: bool| -> Vec<f32> {
            (0..n * n)
                .into_par_iter()
                .map(|k| {
                    let (x, y) = (k % n, k / n);
                    let mut acc = 0f32;
                    for t in 0..r {
                        let o = (t + n - h0) % n;
                        acc += if horiz { s[y * n + (x + o) % n] } else { s[((y + o) % n) * n + x] };
                    }
                    acc / r as f32
                })
                .collect()
        };
        pass(&pass(src, true), false)
    };
    let z0: Vec<f32> = zm.iter().zip(land).map(|(v, l)| if *l { *v } else { 0.0 }).collect();
    let mut out = Vec::new();
    let mut r = 1usize;
    let mut prev = blur(&z0, 1);
    while 2 * r <= n / 8 {
        let next = blur(&z0, 2 * r);
        let (mut s, mut c) = (0f64, 0usize);
        for k in 0..n * n {
            if land[k] {
                let d = (prev[k] - next[k]) as f64;
                s += d * d;
                c += 1;
            }
        }
        out.push((s / c.max(1) as f64).sqrt() as f32);
        prev = next;
        r *= 2;
    }
    out
}

fn read(z: &GridF32) -> (Reading, Hydro) {
    let n = z.width;
    let cell_km = DOMAIN_KM / n as f32;
    let cell_m = cell_km * 1000.0;
    let land = land_of(z);
    let zm: Vec<f32> = z.data.iter().zip(&land).map(|(v, l)| if *l { *v } else { 0.0 }).collect();
    let mut r = Reading { land: land.iter().filter(|&&l| l).count(), ..Default::default() };
    let sl: Vec<(bool, bool, bool)> = (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % n, k / n);
            if x == 0 || y == 0 || x + 1 == n || y + 1 == n || !land[k] {
                return (false, false, false);
            }
            let gx = (zm[k + 1] - zm[k - 1]) / (2.0 * cell_m);
            let gy = (zm[k + n] - zm[k - n]) / (2.0 * cell_m);
            let s = (gx * gx + gy * gy).sqrt().atan().to_degrees();
            (true, (s - 28.0).abs() <= 0.5, (s - 33.0).abs() <= 0.5)
        })
        .collect();
    let inner = sl.iter().filter(|s| s.0).count().max(1);
    r.walls28_cells = sl.iter().filter(|s| s.1).count();
    r.walls28 = 100.0 * r.walls28_cells as f64 / inner as f64;
    r.walls33 = 100.0 * sl.iter().filter(|s| s.2).count() as f64 / inner as f64;
    let crest = crest_mask(&zm, &land, n);
    let lap = laplacian(&zm, n, cell_km);
    let curv = sorted(crest.iter().zip(&lap).filter(|(c, _)| **c).map(|(_, l)| -l).collect());
    r.crest_share = 100.0 * curv.len() as f64 / r.land.max(1) as f64;
    r.crest_p50 = pct(&curv, 0.5);
    r.crest_p90 = pct(&curv, 0.9);
    let (f, fc) = facets(&crest, &land, &lap, n);
    r.facets = f;
    r.facet_cells = fc;
    let gz = GridF32 { width: n, height: n, data: zm.clone() };
    let a = aniso(&gz, &land, 8);
    r.coherent = a.share_hi;
    r.r8 = a.r8;
    r.r8_windows = a.windows;
    r.lambda = valley_spacing(&zm, &land, n);
    let b = n / 32;
    let mut rel = Vec::new();
    for by in 0..32 {
        'blk: for bx in 0..32 {
            let (mut lo, mut hi) = (f32::MAX, f32::MIN);
            for y in by * b..(by + 1) * b {
                for x in bx * b..(bx + 1) * b {
                    let k = y * n + x;
                    if !land[k] {
                        continue 'blk;
                    }
                    lo = lo.min(zm[k]);
                    hi = hi.max(zm[k]);
                }
            }
            rel.push(hi - lo);
        }
    }
    r.relief_p50 = pct(&sorted(rel), 0.5);
    let alts = sorted(zm.iter().zip(&land).filter(|(_, l)| **l).map(|(v, _)| *v).collect());
    r.mean = alts.iter().sum::<f32>() / alts.len().max(1) as f32;
    r.p90 = pct(&alts, 0.9);
    r.peak = alts.last().copied().unwrap_or(0.0);
    let hy = hydrology(&gz, cell_km);
    for k in 0..n * n {
        if land[k] && hy.acc_km2[k] >= 100.0 {
            r.trunks.push((k % n, k / n));
        }
    }
    r.bands = bands(&zm, &land, n);
    r.lakes = hy.lakes.len();
    r.lake_km2 = hy.lakes.iter().map(|l| l.area_km2).sum();
    (r, hy)
}

fn trunk_dist(fine: &[(usize, usize)], coarse: &[(usize, usize)], ratio: f32, n_coarse: usize) -> (f32, f32, f32) {
    if fine.is_empty() || coarse.is_empty() {
        return (f32::NAN, f32::NAN, f32::NAN);
    }
    let nc = n_coarse as f32;
    let mut d: Vec<f32> = fine
        .par_iter()
        .map(|&(x, y)| {
            let (fx, fy) = (x as f32 / ratio, y as f32 / ratio);
            coarse
                .iter()
                .map(|&(cx, cy)| {
                    let mut dx = (fx - cx as f32).abs();
                    let mut dy = (fy - cy as f32).abs();
                    dx = dx.min(nc - dx);
                    dy = dy.min(nc - dy);
                    dx * dx + dy * dy
                })
                .fold(f32::MAX, f32::min)
                .sqrt()
        })
        .collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (pct(&d, 0.5), pct(&d, 0.9), *d.last().unwrap())
}

/// D6 / D7.2: level n+1 restricted − level n on the land of both: (bias m, mean |Δ| m).
fn drift(fine: &GridF32, coarse: &GridF32) -> (f32, f32) {
    let back = restrict(fine);
    let (lb, lc) = (land_of(&back), land_of(coarse));
    let (mut s, mut a, mut c) = (0f64, 0f64, 0usize);
    for k in 0..coarse.data.len() {
        if lb[k] && lc[k] {
            let d = (back.data[k] - coarse.data[k]) as f64;
            s += d;
            a += d.abs();
            c += 1;
        }
    }
    ((s / c.max(1) as f64) as f32, (a / c.max(1) as f64) as f32)
}

/// D7.6: the physics level's 5 largest depressions, followed at a level: (area km², spill m) of the lake containing
/// (or nearest within `f` cells of) the physics lowest point.
fn follow(hy: &Hydro, phys_lowest: &[(usize, usize)], phys_n: usize) -> Vec<Option<(f32, f32)>> {
    let n = hy.n;
    let f = n / phys_n;
    phys_lowest
        .iter()
        .map(|&(x, y)| {
            let (cx, cy) = (x * f + f / 2, y * f + f / 2);
            let r = f.max(2) as i32;
            let mut best: Option<(i32, u32)> = None;
            for dy in -r..=r {
                for dx in -r..=r {
                    let (xx, yy) = (cx as i32 + dx, cy as i32 + dy);
                    if xx < 0 || yy < 0 || xx >= n as i32 || yy >= n as i32 {
                        continue;
                    }
                    let id = hy.lake_id[yy as usize * n + xx as usize];
                    let d2 = dx * dx + dy * dy;
                    if id != 0 && best.is_none_or(|b| d2 < b.0) {
                        best = Some((d2, id));
                    }
                }
            }
            best.map(|(_, id)| {
                let l = &hy.lakes[id as usize - 1];
                (l.area_km2, l.spill_m)
            })
        })
        .collect()
}

// ───────────────────────────── images ─────────────────────────────

fn rgba_of(z: &GridF32, hy: Option<&Hydro>, ss: &SteinSteinParams) -> Vec<u8> {
    let n = z.width;
    let mut rgba = shade_m_rgba(z, DOMAIN_KM / n as f32 * 1000.0, ss.depth_scale_m as f32);
    if let Some(h) = hy {
        draw(&mut rgba, h, 2, true, true);
    }
    rgba
}
/// A crop (512² origin, 128 at 512²) of a north-up RGBA image of an n-level, resized (nearest) to `px`.
fn crop_tile(rgba: &[u8], n: usize, o512: (usize, usize), px: usize) -> image::RgbaImage {
    let f = n as f32 / 512.0;
    let side = (128.0 * f) as usize;
    let (ox, oy) = ((o512.0 as f32 * f) as usize, (o512.1 as f32 * f) as usize);
    // the image is north up: field row y is image row n − 1 − y
    let mut img = image::RgbaImage::new(side as u32, side as u32);
    for r in 0..side {
        let fy = oy + side - 1 - r;
        let ir = n - 1 - fy;
        for c in 0..side {
            let o = (ir * n + ox + c) * 4;
            img.put_pixel(c as u32, r as u32, image::Rgba([rgba[o], rgba[o + 1], rgba[o + 2], 255]));
        }
    }
    image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Nearest)
}
fn compose(rows: &[Vec<image::RgbaImage>], px: usize) -> image::RgbaImage {
    let cols = rows.iter().map(|r| r.len()).max().unwrap();
    let gap = 6;
    let mut out =
        image::RgbaImage::from_pixel((cols * (px + gap)) as u32, (rows.len() * (px + gap)) as u32, image::Rgba([24, 24, 24, 255]));
    for (r, row) in rows.iter().enumerate() {
        for (c, t) in row.iter().enumerate() {
            image::imageops::overlay(&mut out, t, (c * (px + gap)) as i64, (r * (px + gap)) as i64);
        }
    }
    out
}
/// The hypsometric curves (land altitude against the share of land above it), one colour per level, the target band.
fn hypsometry_png(curves: &[(String, Vec<f32>, [u8; 3])], path: &Path) {
    let (w, h) = (900u32, 520u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let zmax = 4000f32;
    let y_of = |m: f32| (h as f32 - 20.0 - (m / zmax).clamp(0.0, 1.0) * (h as f32 - 40.0)) as u32;
    for x in 40..w - 20 {
        for m in [PEAK_OK.0, PEAK_OK.1] {
            img.put_pixel(x, y_of(m), image::Rgba([230, 160, 160, 255]));
        }
        img.put_pixel(x, y_of(0.0), image::Rgba([0, 0, 0, 255]));
        for m in [1000.0, 2000.0, 3000.0] {
            img.put_pixel(x, y_of(m), image::Rgba([210, 210, 210, 255]));
        }
    }
    for (_, alts, c) in curves {
        if alts.is_empty() {
            continue;
        }
        for xi in 0..(w - 60) {
            let frac = xi as f32 / (w - 61) as f32;
            let idx = ((1.0 - frac) * (alts.len() - 1) as f32) as usize;
            let y = y_of(alts[idx]);
            for dy in 0..2 {
                img.put_pixel(40 + xi, (y + dy).min(h - 1), image::Rgba([c[0], c[1], c[2], 255]));
            }
        }
    }
    img.save(path).unwrap();
}

// ───────────────────────────── the run ─────────────────────────────

#[derive(Clone)]
struct LevelOut {
    n: usize,
    reading: Reading,
    result: GridF32,
    secs: [f64; 4],
    done: Budget,
    k: f32,
    s_ref: f32,
    depth: f32,
    carved_km3: f64,
    bound_km3: f64,
    deepest_m: f32,
    bound_cell_m: f32,
    smax_share: f64,
    lifted: (usize, f32),
    retarget_peak: f32,
    drift: Option<(f32, f32, f32)>,
    trunk_prev: Option<(f32, f32, f32)>,
    trunk_ref: (f32, f32, f32),
    follow: Vec<Option<(f32, f32)>>,
    retained_mb: f64,
}

struct ChainOut {
    name: String,
    levels: Vec<LevelOut>,
    final_retarget: Option<(GridF32, Reading, usize)>,
    stop: Option<(usize, String)>,
}

#[allow(clippy::too_many_arguments)]
fn run_chain(
    name: &str,
    phys_m: &GridF32,
    phys_reading: &Reading,
    cfg: AmpConfig,
    top: usize,
    refs: &BTreeMap<usize, (Reading, GridF32)>,
    controls: &BTreeMap<usize, (f64, f64)>,
    phys_lowest: &[(usize, usize)],
    keep_results: bool,
) -> ChainOut {
    let t0 = Instant::now();
    let mut ch = Chain::new(phys_m.clone(), cfg.clone());
    let mut outs: Vec<LevelOut> = Vec::new();
    let mut prev_reading = phys_reading.clone();
    let mut prev_result = phys_m.clone();
    let mut stop: Option<(usize, String)> = None;
    while ch.last_result().width < top {
        let lvl: AmpLevel = ch.next_level(&|| false).unwrap().clone();
        let n = lvl.n;
        let result = lvl.result();
        let (reading, hy) = read(&result);
        let c2 = (lvl.cell_km * lvl.cell_km) as f64;
        let land_cells = lvl.ocean.iter().filter(|o| !**o).count() as f64;
        let carved: f64 = lvl.sum_erosion.iter().zip(&lvl.ocean).filter(|(_, o)| !**o).map(|(v, _)| -*v as f64).sum();
        let bound_cell = lvl.done.erosion as f32 * lvl.k * cfg.s_max.powf(cfg.n_exp) * cfg.a_max.powf(cfg.m_exp);
        let deepest = lvl.sum_erosion.iter().fold(0f32, |a, &v| a.max(-v));
        let h0 = ch.physics_at(n);
        let rt = ymir_core::cascade::amplify::retarget(&result, &h0, &lvl.ocean, &cfg).0;
        let rt_peak = rt.data.iter().zip(&lvl.ocean).filter(|(_, o)| !**o).map(|(v, _)| *v).fold(f32::MIN, f32::max);
        let (bias, absd) = drift(&result, &prev_result);
        let rel = bias / prev_reading.mean.max(1.0);
        let tp = trunk_dist(&reading.trunks, &prev_reading.trunks, 2.0, n / 2);
        let tr = refs.get(&n).map_or((f32::NAN, f32::NAN, f32::NAN), |r| trunk_dist(&reading.trunks, &r.0.trunks, 1.0, n));
        let out = LevelOut {
            n,
            reading: reading.clone(),
            result: if keep_results { result.clone() } else { GridF32::new(1, 1, 0.0) },
            secs: lvl.secs,
            done: lvl.done,
            k: lvl.k,
            s_ref: lvl.s_ref,
            depth: lvl.depth_m,
            carved_km3: carved * c2 * 1e-3,
            bound_km3: bound_cell as f64 * land_cells * c2 * 1e-3,
            deepest_m: deepest,
            bound_cell_m: bound_cell,
            smax_share: 100.0 * lvl.smax_bind.0 as f64 / lvl.smax_bind.1.max(1) as f64,
            lifted: (lvl.lifted, lvl.lifted_mean_depth_m),
            retarget_peak: rt_peak,
            drift: Some((bias, absd, rel)),
            trunk_prev: Some(tp),
            trunk_ref: tr,
            follow: follow(&hy, phys_lowest, phys_m.width),
            retained_mb: lvl.retained_bytes() as f64 / 1e6,
        };
        // R
        if stop.is_none() {
            let mut why = Vec::new();
            if let Some(r) = refs.get(&n) {
                let ctrl = controls.get(&n).is_some_and(|c| c.0 >= 60.0 && c.1 <= 20.0);
                if reading.walls28 > r.0.walls28 && reading.walls28_cells > 20 {
                    why.push(format!("walls 28° {:.3} > {:.3} ({} cells)", reading.walls28, r.0.walls28, reading.walls28_cells));
                }
                if ctrl && reading.facets > r.0.facets && reading.facet_cells > 20 {
                    why.push(format!("facets {:.2} > {:.2} ({} crests)", reading.facets, r.0.facets, reading.facet_cells));
                }
            }
            if rel.abs() > 0.10 {
                why.push(format!("drift {:+.1} %", 100.0 * rel));
            }
            if rt_peak < PEAK_OK.0 || rt_peak > PEAK_OK.1 {
                why.push(format!("peak after retargeting {:.0} m", rt_peak));
            }
            if tp.1 > 1.5 {
                why.push(format!("trunk p90 {:.2} cells", tp.1));
            }
            if !why.is_empty() {
                stop = Some((n, why.join("; ")));
            }
        }
        eprintln!(
            "      {name} {n}² done ({:.1} s: E {:.1} · T {:.1} · D {:.1}) · peak {:.0} → retargeted {:.0} m · drift {:+.1} %",
            lvl.secs.iter().sum::<f64>(),
            lvl.secs[1],
            lvl.secs[2],
            lvl.secs[3],
            reading.peak,
            rt_peak,
            100.0 * rel
        );
        prev_reading = reading;
        prev_result = result;
        outs.push(out);
        // keep only the last level's internals
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let final_retarget = ch.retargeted().map(|(g, nc)| {
        let (r, _) = read(&g);
        (if keep_results { g } else { GridF32::new(1, 1, 0.0) }, r, nc)
    });
    eprintln!("   chain {name} to {top}² in {:.1} s", t0.elapsed().as_secs_f64());
    ChainOut { name: name.to_string(), levels: outs, final_retarget, stop }
}

fn print_chain(c: &ChainOut, refs: &BTreeMap<usize, (Reading, GridF32)>, controls: &BTreeMap<usize, (f64, f64)>) {
    eprintln!("\n   ═══ chain {} ═══", c.name);
    for l in &c.levels {
        let r = &l.reading;
        let rf = refs.get(&l.n).map(|x| &x.0);
        let ctrl = controls.get(&l.n).copied().unwrap_or((f64::NAN, f64::NAN));
        let g = |f: &dyn Fn(&Reading) -> String| rf.map_or("—".to_string(), f);
        eprintln!(
            "   ── {}² · E {} · T {} · D {} iterations · {:.1} s [upscale {:.2} · E {:.1} ({:.2} ms/it) · T {:.1} ({:.2} ms/it) · D {:.1} ({:.2} ms/it)] · retained {:.0} MB",
            l.n,
            l.done.erosion,
            l.done.talus,
            l.done.deposit,
            l.secs.iter().sum::<f64>(),
            l.secs[0],
            l.secs[1],
            1e3 * l.secs[1] / l.done.erosion.max(1) as f64,
            l.secs[2],
            1e3 * l.secs[2] / l.done.talus.max(1) as f64,
            l.secs[3],
            1e3 * l.secs[3] / l.done.deposit.max(1) as f64,
            l.retained_mb
        );
        eprintln!(
            "      k {:.3e} · s_ref (p90) {:.4} · d_L {:.0} m · carved {:.2} km³ of a bound {:.0} km³ ({:.2} %) · deepest {:.1} m of a per-cell bound {:.0} m · s ≥ s_max on {:.3} % of the erosion updates · D5 lifted {} cells (mean depth {:.0} m)",
            l.k, l.s_ref, l.depth, l.carved_km3, l.bound_km3, 100.0 * l.carved_km3 / l.bound_km3.max(1e-9), l.deepest_m, l.bound_cell_m, l.smax_share, l.lifted.0, l.lifted.1
        );
        let row = |name: &str, a: String, b: String| eprintln!("      {name:<36} cascade {a:>18} · témoin {b:>18}");
        row("peak / mean / p90 land (m)", format!("{:.0} / {:.0} / {:.0}", r.peak, r.mean, r.p90), g(&|x| format!("{:.0} / {:.0} / {:.0}", x.peak, x.mean, x.p90)));
        eprintln!("      peak after retargeting (evaluation): {:.0} m (target {} – {}, R band {:.0} – {:.0})", l.retarget_peak, 2700, 3000, PEAK_OK.0, PEAK_OK.1);
        if let Some((b, a, rel)) = l.drift {
            eprintln!("      drift (restricted − previous level): bias {:+.1} m ({:+.2} % of its mean land) · mean |Δ| {:.1} m", b, 100.0 * rel, a);
        }
        if let Some(t) = l.trunk_prev {
            eprintln!("      trunks ≥ 100 km² → previous level (cells of it): p50 {:.2} / p90 {:.2} / max {:.2}", t.0, t.1, t.2);
        }
        eprintln!("      trunks → témoin's at the same level (cells): p50 {:.2} / p90 {:.2} / max {:.2}", l.trunk_ref.0, l.trunk_ref.1, l.trunk_ref.2);
        row("walls 28° ± 0.5° (% / cells)", format!("{:.3} / {}", r.walls28, r.walls28_cells), g(&|x| format!("{:.3} / {}", x.walls28, x.walls28_cells)));
        row("band 33° ± 0.5° (%)", format!("{:.3}", r.walls33), g(&|x| format!("{:.3}", x.walls33)));
        row("crest cells (% land)", format!("{:.2}", r.crest_share), g(&|x| format!("{:.2}", x.crest_share)));
        row("crest −∇²z p50 / p90 (m/km²)", format!("{:.1} / {:.1}", r.crest_p50, r.crest_p90), g(&|x| format!("{:.1} / {:.1}", x.crest_p50, x.crest_p90)));
        row("crest facets (% / crests)", format!("{:.2} / {}", r.facets, r.facet_cells), g(&|x| format!("{:.2} / {}", x.facets, x.facet_cells)));
        eprintln!("      facet control: pyramids {:.1} % · paraboloids {:.1} %", ctrl.0, ctrl.1);
        row("coherent windows C > 0.7 (%)", format!("{:.2}", r.coherent), g(&|x| format!("{:.2}", x.coherent)));
        row("terrain R8 (windows)", format!("{:.4} ({})", r.r8, r.r8_windows), g(&|x| format!("{:.4} ({})", x.r8, x.r8_windows)));
        row("valley spacing λ (cells)", format!("{:.2}", r.lambda), g(&|x| format!("{:.2}", x.lambda)));
        row("relief p50, 12.5 km blocks (m)", format!("{:.0}", r.relief_p50), g(&|x| format!("{:.0}", x.relief_p50)));
        row("lakes (count / km²)", format!("{} / {:.0}", r.lakes, r.lake_km2), g(&|x| format!("{} / {:.0}", x.lakes, x.lake_km2)));
        let cell_km = DOMAIN_KM / l.n as f32;
        let bnames: Vec<String> = (0..r.bands.len()).map(|k| format!("{:.2}–{:.2} km", cell_km * (1 << k) as f32, cell_km * (2 << k) as f32)).collect();
        eprintln!(
            "      detail per band (RMS m), cascade / témoin: {}",
            (0..r.bands.len())
                .map(|k| format!("[{}] {:.1} / {}", bnames[k], r.bands[k], rf.map_or("—".to_string(), |x| format!("{:.1}", x.bands.get(k).copied().unwrap_or(f32::NAN)))))
                .collect::<Vec<_>>()
                .join(" · ")
        );
        eprintln!(
            "      the physics level's 5 largest depressions here (area km² / spill m): {}",
            l.follow.iter().map(|f| f.map_or("gone".to_string(), |(a, s)| format!("{a:.0} / {s:.0}"))).collect::<Vec<_>>().join(" · ")
        );
    }
    if let Some((_, r, nc)) = &c.final_retarget {
        eprintln!(
            "   final retargeting ({} constraint cells): peak {:.0} · mean {:.0} · p90 {:.0} m · walls 28° {:.3} % · facets {:.2} % · R8 {:.4}",
            nc, r.peak, r.mean, r.p90, r.walls28, r.facets, r.r8
        );
    }
    match &c.stop {
        Some((n, why)) => eprintln!("   R — {} FIRES at {}²: {}", c.name, n, why),
        None => eprintln!("   R — {} does not fire up to {}²", c.name, c.levels.last().map_or(0, |l| l.n)),
    }
}

#[test]
#[ignore]
fn f160_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    std::fs::create_dir_all(cache_dir()).unwrap();
    std::fs::create_dir_all(out_dir()).unwrap();
    eprintln!("\n==========  Finding 160 . the cascade, round 2: a physics level, then Schott 2024 transposed  ==========");
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let cpath = cache_dir().join("coarse64.bin");
    if !(cpath.exists() && REF_LEVELS.iter().all(|n| cache_dir().join(format!("ref{n}.bin")).exists())) {
        let t = Instant::now();
        let wd = build_world(temoin, None, PSEED, None);
        let coarse = c1_coarse_normalized_altitude(&wd.state, &IsostasyConfig::c1_default(), &ss, wd.cfg.target_land_fraction);
        save_grid(&cpath, &coarse);
        for n in REF_LEVELS {
            save_grid(&cache_dir().join(format!("ref{n}.bin")), &block_mean(&wd.heightmap, n));
        }
        drop(wd);
        eprintln!("   the témoin built and reduced ({:.0} s)", t.elapsed().as_secs_f64());
    }
    let coarse = roll(&load_grid(&cpath).unwrap(), off);
    // the témoin's readings per level (metres)
    let t = Instant::now();
    let mut refs: BTreeMap<usize, (Reading, GridF32)> = BTreeMap::new();
    for n in REF_LEVELS {
        let g = to_m(&load_grid(&cache_dir().join(format!("ref{n}.bin"))).unwrap(), n2m);
        let (r, _) = read(&g);
        refs.insert(n, (r, g));
    }
    let controls: BTreeMap<usize, (f64, f64)> = REF_LEVELS.iter().map(|&n| (n, facet_control(n))).collect();
    eprintln!("   the témoin read at {:?} ({:.0} s)", REF_LEVELS, t.elapsed().as_secs_f64());
    // ── D1: the physics levels ──
    let ccfg = CascadeConfig::declared(ss.depth_scale_m as f32);
    let mut phys: BTreeMap<usize, (GridF32, Reading, Vec<(usize, usize)>)> = BTreeMap::new();
    for (n, cap) in [(128usize, 400usize), (256, 300)] {
        let (rec, cal) = physics_level(&coarse, n, &ccfg, PEAK_TARGET_M, cap, &mut |_| {}, &|| false).unwrap();
        let pm = to_m(&rec.z, n2m);
        let (r, hy) = read(&pm);
        let mut big: Vec<&ymir_core::cascade::hydro::Lake> = hy.lakes.iter().collect();
        big.sort_by(|a, b| b.cells.cmp(&a.cells));
        let lowest: Vec<(usize, usize)> = big.iter().take(5).map(|l| (l.lowest % n, l.lowest / n)).collect();
        eprintln!(
            "\n   physics P{n}: U₀ trial {:.1e} → peak {:.0} m → U₀ = {:.3e} m/yr per km of h_iso · trial {} steps ({}) · run {} steps ({}, last ratio {:.4}) · {:.1} s · talus {:.2} % of land at the last step · D5 lifted {} cells (mean depth {:.0} m)",
            cal.u0_trial,
            cal.trial_peak_m,
            cal.u0,
            cal.steps_trial,
            if cal.equilibrium_trial { "eq" } else { "CAP" },
            rec.steps,
            if rec.equilibrium { "equilibrium" } else { "CAP" },
            rec.ratios.last().copied().unwrap_or(f32::NAN),
            cal.secs,
            100.0 * rec.talus_last_share,
            cal.lifted,
            cal.lifted_mean_depth_m
        );
        let rf = &refs[&n].0;
        eprintln!(
            "      peak / mean / p90 {:.0} / {:.0} / {:.0} m (témoin {:.0} / {:.0} / {:.0}) · relief p50 {:.0} (témoin {:.0}) · walls 28° {:.3} % · facets {:.2} % (témoin {:.2}) · R8 {:.4} (témoin {:.4}) · lakes {} / {:.0} km² · the 5 largest (km² / spill m): {}",
            r.peak, r.mean, r.p90, rf.peak, rf.mean, rf.p90, r.relief_p50, rf.relief_p50, r.walls28, r.facets, rf.facets, r.r8, rf.r8, r.lakes, r.lake_km2,
            big.iter().take(5).map(|l| format!("{:.0} / {:.0}", l.area_km2, l.spill_m)).collect::<Vec<_>>().join(" · ")
        );
        eprintln!(
            "      detail per band (RMS m), physics / témoin: {}",
            r.bands.iter().zip(&rf.bands).map(|(a, b)| format!("{a:.1} / {b:.1}")).collect::<Vec<_>>().join(" · ")
        );
        save_grid(&cache_dir().join(format!("phys{n}.bin")), &pm);
        phys.insert(n, (pm, r, lowest));
    }
    // ── D6: the chains ──
    let mut chains: Vec<ChainOut> = Vec::new();
    let mut runs: Vec<(String, usize, Variant, f32, usize)> = Vec::new();
    for v in [Variant::S, Variant::SRho, Variant::SU] {
        for (b, f) in [("faible", 0.5f32), ("moyen", 1.0), ("fort", 2.0)] {
            runs.push((format!("P128 {} {b}", v.label()), 128, v, f, 1024));
        }
    }
    for v in [Variant::S, Variant::SRho] {
        runs.push((format!("P256 {} moyen", v.label()), 256, v, 1.0, 1024));
    }
    for (name, pn, v, f, top) in &runs {
        let (pm, pr, low) = &phys[pn];
        let cfg = AmpConfig::declared(PSEED, *v, *f);
        let keep = *f == 1.0 && matches!(v, Variant::S | Variant::SRho);
        chains.push(run_chain(name, pm, pr, cfg, *top, &refs, &controls, low, keep));
    }
    for c in &chains {
        print_chain(c, &refs, &controls);
    }
    // the best (D-R): the highest level reached without R firing — the physics level counts as reached (corrected after
    // run 1, whose ranking ignored it and so fell back to the first chain), ties → the lowest R8 there
    let reach = |c: &ChainOut| -> (usize, f32) {
        let phys_n = c.levels.first().map_or(0, |l| l.n / 2);
        let top = match &c.stop {
            Some((n, _)) => c.levels.iter().filter(|l| l.n < *n).map(|l| l.n).max().unwrap_or(phys_n),
            None => c.levels.last().map_or(phys_n, |l| l.n),
        };
        let r8 = c.levels.iter().find(|l| l.n == top).map_or_else(|| phys[&phys_n].1.r8, |l| l.reading.r8);
        (top, r8)
    };
    let mut ranked: Vec<(usize, f32, String)> = chains.iter().map(|c| {
        let (t, r) = reach(c);
        (t, r, c.name.clone())
    }).collect();
    ranked.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.partial_cmp(&b.1).unwrap()));
    eprintln!("\n   the ranking (highest level reached without R, then the lowest R8 there):");
    for (t, r, nm) in &ranked {
        eprintln!("      {nm:<22} reaches {t}² (R8 there {r:.4})");
    }
    let best_name = ranked[0].2.clone();
    let best_phys = if best_name.starts_with("P256") { 256 } else { 128 };
    eprintln!("   the best: {best_name}; 2048² runs from P{best_phys}");
    // ── 2 048² for S, S+ρ, S+U moyen from the best physics level ──
    let mut big: Vec<ChainOut> = Vec::new();
    for v in [Variant::S, Variant::SRho, Variant::SU] {
        let (pm, pr, low) = &phys[&best_phys];
        let cfg = AmpConfig::declared(PSEED, v, 1.0);
        big.push(run_chain(&format!("P{best_phys} {} moyen → 2048", v.label()), pm, pr, cfg, 2048, &refs, &controls, low, true));
    }
    for c in &big {
        print_chain(c, &refs, &controls);
    }
    // ── S against S+ρ: the axis alignment ──
    eprintln!("\n   S against S+ρ (same physics level and budget): terrain R8 / coherent windows per level");
    for (a, b) in [(&chains[1], &chains[4]), (&big[0], &big[1])] {
        for (la, lb) in a.levels.iter().zip(&b.levels) {
            eprintln!(
                "      {}² · {} R8 {:.4} / coherent {:.2} % · {} R8 {:.4} / coherent {:.2} % · ΔR8 {:+.1} % · Δcoherent {:+.1} %",
                la.n,
                a.name,
                la.reading.r8,
                la.reading.coherent,
                b.name,
                lb.reading.r8,
                lb.reading.coherent,
                100.0 * (lb.reading.r8 - la.reading.r8) / la.reading.r8.max(1e-9),
                100.0 * (lb.reading.coherent - la.reading.coherent) / la.reading.coherent.max(1e-9)
            );
        }
    }
    // ── the cost, per process, extrapolated to 8 192² ──
    let best_big = &big[0];
    let last = best_big.levels.last().unwrap();
    let per = [
        last.secs[1] / last.done.erosion.max(1) as f64,
        last.secs[2] / last.done.talus.max(1) as f64,
        last.secs[3] / last.done.deposit.max(1) as f64,
    ];
    let measured: f64 = best_big.levels.iter().map(|l| l.secs.iter().sum::<f64>()).sum();
    let mut tot = [0f64; 3];
    for nn in [4096usize, 8192] {
        let b = Budget::moyen(nn);
        let sc = (nn as f64 / last.n as f64).powi(2);
        tot[0] += per[0] * sc * b.erosion as f64;
        tot[1] += per[1] * sc * b.talus as f64;
        tot[2] += per[2] * sc * b.deposit as f64;
    }
    let total = measured + tot.iter().sum::<f64>();
    eprintln!(
        "\n   cost: measured {} 128² → 2048² {:.0} s · per iteration at {}²: E {:.1} · T {:.1} · D {:.1} ms · ESTIMATE 4096² + 8192²: E {:.0} s · T {:.0} s · D {:.0} s → total {:.0} s = {:.1} min (alert above 1 h: {})",
        best_big.name,
        measured,
        last.n,
        1e3 * per[0],
        1e3 * per[1],
        1e3 * per[2],
        tot[0],
        tot[1],
        tot[2],
        total,
        total / 60.0,
        if total > 3600.0 { "YES" } else { "no" }
    );
    // ── I: the images ──
    let px = 512usize;
    let best = chains.iter().chain(big.iter()).find(|c| c.name.starts_with(&format!("P{best_phys} S moyen → 2048"))).unwrap();
    let srho = &big[1];
    let mut hyp: Vec<(String, Vec<f32>, [u8; 3])> = Vec::new();
    let colours = [[200, 60, 40], [220, 140, 30], [60, 150, 60], [40, 110, 200], [120, 60, 180], [30, 30, 30]];
    let (pm, _, _) = &phys[&best_phys];
    let land_alts = |g: &GridF32| sorted(g.data.iter().zip(land_of(g)).filter(|(_, l)| *l).map(|(v, _)| *v).collect());
    hyp.push((format!("physics {best_phys}²"), land_alts(pm), colours[0]));
    for (i, l) in best.levels.iter().enumerate() {
        hyp.push((format!("{}²", l.n), land_alts(&l.result), colours[(i + 1) % colours.len()]));
    }
    if let Some((g, _, _)) = &best.final_retarget {
        hyp.push(("retargeted".into(), land_alts(g), [0, 0, 0]));
    }
    hypsometry_png(&hyp, &out_dir().join("f160_hypsometry.png"));
    eprintln!(
        "   hypsometry (f160_hypsometry.png), colours: {}",
        hyp.iter().map(|(n, _, c)| format!("{n} rgb{:?}", c)).collect::<Vec<_>>().join(" · ")
    );
    // whole views per level (best and témoin), rivers (Strahler ≥ 2) and lakes drawn
    let mut views: Vec<(usize, GridF32)> = vec![(best_phys, pm.clone())];
    for l in &best.levels {
        views.push((l.n, l.result.clone()));
    }
    for (n, g) in &views {
        let (_, hy) = read(g);
        let img = image::RgbaImage::from_raw(*n as u32, *n as u32, rgba_of(g, Some(&hy), &ss)).unwrap();
        img.save(out_dir().join(format!("f160_best_{n}.png"))).unwrap();
        if let Some((_, rg)) = refs.get(n) {
            let (_, rh) = read(rg);
            image::RgbaImage::from_raw(*n as u32, *n as u32, rgba_of(rg, Some(&rh), &ss))
                .unwrap()
                .save(out_dir().join(format!("f160_temoin_{n}.png")))
                .unwrap();
        }
    }
    if let Some((g, _, _)) = &best.final_retarget {
        let (_, hy) = read(g);
        image::RgbaImage::from_raw(g.width as u32, g.width as u32, rgba_of(g, Some(&hy), &ss))
            .unwrap()
            .save(out_dir().join(format!("f160_best_{}_retargeted.png", g.width)))
            .unwrap();
    }
    // the crops: row 1 the best per level (physics → 2048, retargeted), row 2 the témoin per level
    for (name, ox, oy) in CROPS {
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for (n, g) in &views {
            let (_, hy) = read(g);
            r1.push(crop_tile(&rgba_of(g, Some(&hy), &ss), *n, (ox, oy), px));
            if let Some((_, rg)) = refs.get(n) {
                let (_, rh) = read(rg);
                r2.push(crop_tile(&rgba_of(rg, Some(&rh), &ss), *n, (ox, oy), px));
            }
        }
        if let Some((g, _, _)) = &best.final_retarget {
            let (_, hy) = read(g);
            r1.push(crop_tile(&rgba_of(g, Some(&hy), &ss), g.width, (ox, oy), px));
        }
        compose(&[r1, r2], px).save(out_dir().join(format!("f160_{name}.png"))).unwrap();
        eprintln!("   image f160_{name}.png · 512² origin ({ox}, {oy}) · row 1 the best (physics, 256², 512², 1024², 2048², retargeted) · row 2 the témoin (128²–2048²)");
    }
    // S against S+ρ on crop A, per level
    let (ox, oy) = (CROPS[0].1, CROPS[0].2);
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for (a, b) in best.levels.iter().zip(&srho.levels) {
        let (_, ha) = read(&a.result);
        let (_, hb) = read(&b.result);
        r1.push(crop_tile(&rgba_of(&a.result, Some(&ha), &ss), a.n, (ox, oy), px));
        r2.push(crop_tile(&rgba_of(&b.result, Some(&hb), &ss), b.n, (ox, oy), px));
    }
    compose(&[r1, r2], px).save(out_dir().join("f160_S_vs_Srho_A.png")).unwrap();
    eprintln!("   image f160_S_vs_Srho_A.png · row 1 S, row 2 S+ρ, per level");
    eprintln!("\n==========  end Finding 160 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 160, diagnostic A1 (NOT blind, written after the declared runs): the same chain (P128, S, moyen) with
/// k × 1 (the declaration), × 10 and × 100, to 1024², read with the same instruments. It does not change R's verdict
/// on the declared design; it asks whether the intensity is the missing lever.
///
/// Run (after `f160_cascade`, whose cache it reads): cargo test -p ymir-core --release --test f160_cascade --
/// --ignored --exact f160_diag --nocapture
#[test]
#[ignore]
fn f160_diag() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    eprintln!("\n==========  Finding 160 . diagnostic A1 (NOT blind): the intensity k × 1 / 10 / 100  ==========");
    let mut refs: BTreeMap<usize, (Reading, GridF32)> = BTreeMap::new();
    for n in [128usize, 256, 512, 1024] {
        let g = to_m(&load_grid(&cache_dir().join(format!("ref{n}.bin"))).expect("run f160_cascade first"), n2m);
        let (r, _) = read(&g);
        refs.insert(n, (r, g));
    }
    let controls: BTreeMap<usize, (f64, f64)> = [128usize, 256, 512, 1024].iter().map(|&n| (n, facet_control(n))).collect();
    let pm = load_grid(&cache_dir().join("phys128.bin")).expect("run f160_cascade first");
    let (pr, hy) = read(&pm);
    let mut big: Vec<&ymir_core::cascade::hydro::Lake> = hy.lakes.iter().collect();
    big.sort_by(|a, b| b.cells.cmp(&a.cells));
    let low: Vec<(usize, usize)> = big.iter().take(5).map(|l| (l.lowest % 128, l.lowest / 128)).collect();
    let px = 512usize;
    let (ox, oy) = (CROPS[0].1, CROPS[0].2);
    let mut rows = Vec::new();
    for boost in [1.0f32, 10.0, 100.0] {
        let mut cfg = AmpConfig::declared(PSEED, Variant::S, 1.0);
        cfg.k_boost = boost;
        let c = run_chain(&format!("P128 S moyen k×{boost}"), &pm, &pr, cfg, 1024, &refs, &controls, &low, true);
        print_chain(&c, &refs, &controls);
        let mut row = Vec::new();
        for l in &c.levels {
            let (_, h) = read(&l.result);
            row.push(crop_tile(&rgba_of(&l.result, Some(&h), &ss), l.n, (ox, oy), px));
        }
        let last = &c.levels.last().unwrap().result;
        let (_, h) = read(last);
        image::RgbaImage::from_raw(1024, 1024, rgba_of(last, Some(&h), &ss))
            .unwrap()
            .save(out_dir().join(format!("f160_diag_k{boost}_1024.png")))
            .unwrap();
        rows.push(row);
    }
    let mut trow = Vec::new();
    for n in [256usize, 512, 1024] {
        let g = &refs[&n].1;
        let (_, h) = read(g);
        trow.push(crop_tile(&rgba_of(g, Some(&h), &ss), n, (ox, oy), px));
    }
    rows.push(trow);
    compose(&rows, px).save(out_dir().join("f160_diag_A.png")).unwrap();
    eprintln!("   image f160_diag_A.png · crop A · rows k×1, k×10, k×100, the témoin · columns 256², 512², 1024²");
    eprintln!("\n==========  end Finding 160 diagnostic . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
