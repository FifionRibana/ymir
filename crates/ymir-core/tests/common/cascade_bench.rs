//! ADR Findings 161–162 -- the cascade benches' shared instruments (F161's `read`, in metres, with the D5 exclusion;
//! the trunks, the drift, the images). Copied from `f161_cascade.rs`, which keeps its own copy so its measurement stays
//! reproducible as committed.

#![allow(dead_code)]

use super::{aniso, pct, sorted};
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use ymir_core::cascade::amplify::{ocean_mask, restrict};
use ymir_core::cascade::hydro::{Hydro, draw, hydrology};
use ymir_core::cascade::measure::octave_rms;
use ymir_core::cascade::shade_m_rgba;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;

pub const DOMAIN_KM: f32 = 400.0;
pub const PEAK_OK: (f32, f32) = (0.85 * 2700.0, 1.15 * 3000.0);
pub const NOTICE: &str = "produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved";

pub fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}
pub fn cache_dir() -> PathBuf {
    root().join("target_f162/f162_cache")
}
pub fn out_dir() -> PathBuf {
    root().join("docs/reports/relief_method/f162_predictability/images")
}
pub fn save_grid(p: &Path, g: &GridF32) {
    let mut b = Vec::with_capacity(8 + g.data.len() * 4);
    b.extend_from_slice(&(g.width as u32).to_le_bytes());
    b.extend_from_slice(&(g.height as u32).to_le_bytes());
    for v in &g.data {
        b.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(p, b).unwrap();
}
pub fn load_grid(p: &Path) -> Option<GridF32> {
    let b = std::fs::read(p).ok()?;
    let w = u32::from_le_bytes(b[0..4].try_into().unwrap()) as usize;
    let h = u32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
    let data = b[8..].chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect();
    Some(GridF32 { width: w, height: h, data })
}
pub fn block_mean(g: &GridF32, n: usize) -> GridF32 {
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
pub fn to_m(g: &GridF32, n2m: f32) -> GridF32 {
    GridF32 { width: g.width, height: g.height, data: g.data.iter().map(|v| (v - 0.5) * n2m).collect() }
}

// ───────────────────────────── the instruments (F160's, in metres; F161-C3) ─────────────────────────────

#[derive(Clone, Debug, Default)]
pub struct Reading {
    pub n: usize,
    pub land: usize,
    pub textured: usize,
    pub slope_p50: f32,
    pub slope_p90: f32,
    pub slope28: f64,
    pub slope33: f64,
    pub walls28: f64,
    pub walls28_cells: usize,
    pub crest_share: f64,
    pub facets: f64,
    pub facet_cells: usize,
    pub coherent: f32,
    pub r8: f32,
    pub r8_windows: usize,
    pub lambda_km: f32,
    pub density: f32,
    pub relief_p50: f32,
    pub mean: f32,
    pub p90: f32,
    pub peak: f32,
    pub deciles: Vec<f32>,
    pub trunks: Vec<(usize, usize)>,
    pub octaves: Vec<f32>,
    pub lakes_d5: (usize, f32),
    pub lakes_other: (usize, f32),
}

pub fn crest_mask(zm: &[f32], land: &[bool], n: usize) -> Vec<bool> {
    (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % n, k / n);
            if x == 0 || y == 0 || x + 1 == n || y + 1 == n || !land[k] {
                return false;
            }
            let z = zm[k];
            [(k - 1, k + 1), (k - n, k + n), (k - n - 1, k + n + 1), (k - n + 1, k + n - 1)]
                .iter()
                .any(|&(a, b)| z > zm[a] && z > zm[b])
        })
        .collect()
}
pub fn laplacian(zm: &[f32], n: usize, cell_km: f32) -> Vec<f32> {
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
pub fn facets(crest: &[bool], land: &[bool], lap: &[f32], n: usize) -> (f64, usize) {
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
pub fn facet_control(n: usize, cell_km: f32) -> (f64, f64) {
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
pub fn valley_spacing(zm: &[f32], land: &[bool], n: usize) -> f32 {
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

/// F161-C3 on a field in metres (the sea ≤ 0 m connected to the border); `d5` excluded from the texture.
pub fn read(z: &GridF32, d5: Option<&[bool]>) -> (Reading, Hydro) {
    let n = z.width;
    let cell_km = DOMAIN_KM / n as f32 * if z.width == z.height { 1.0 } else { 1.0 };
    read_cell(z, d5, cell_km)
}
pub fn read_cell(z: &GridF32, d5: Option<&[bool]>, cell_km: f32) -> (Reading, Hydro) {
    let n = z.width;
    let cell_m = cell_km * 1000.0;
    let ocean = ocean_mask(z);
    let land: Vec<bool> = ocean.iter().map(|o| !o).collect();
    let tex: Vec<bool> = (0..n * n).map(|k| land[k] && !d5.is_some_and(|m| m[k])).collect();
    let zm: Vec<f32> = z.data.iter().zip(&land).map(|(v, l)| if *l { *v } else { 0.0 }).collect();
    let mut r = Reading {
        n,
        land: land.iter().filter(|&&l| l).count(),
        textured: tex.iter().filter(|&&l| l).count(),
        ..Default::default()
    };
    // slopes (central differences)
    let sl: Vec<Option<f32>> = (0..n * n)
        .into_par_iter()
        .map(|k| {
            let (x, y) = (k % n, k / n);
            if x == 0 || y == 0 || x + 1 == n || y + 1 == n || !tex[k] {
                return None;
            }
            let gx = (zm[k + 1] - zm[k - 1]) / (2.0 * cell_m);
            let gy = (zm[k + n] - zm[k - n]) / (2.0 * cell_m);
            Some((gx * gx + gy * gy).sqrt())
        })
        .collect();
    let s = sorted(sl.iter().flatten().copied().collect());
    r.slope_p50 = pct(&s, 0.5);
    r.slope_p90 = pct(&s, 0.9);
    let deg = |v: f32| v.atan().to_degrees();
    r.slope28 = 100.0 * s.iter().filter(|&&v| deg(v) > 28.0).count() as f64 / s.len().max(1) as f64;
    r.slope33 = 100.0 * s.iter().filter(|&&v| deg(v) > 33.0).count() as f64 / s.len().max(1) as f64;
    r.walls28_cells = s.iter().filter(|&&v| (deg(v) - 28.0).abs() <= 0.5).count();
    r.walls28 = 100.0 * r.walls28_cells as f64 / s.len().max(1) as f64;
    // crests and facets
    let crest = crest_mask(&zm, &tex, n);
    let lap = laplacian(&zm, n, cell_km);
    r.crest_share = 100.0 * crest.iter().filter(|&&c| c).count() as f64 / r.textured.max(1) as f64;
    let (f, fc) = facets(&crest, &tex, &lap, n);
    r.facets = f;
    r.facet_cells = fc;
    let gz = GridF32 { width: n, height: n, data: zm.clone() };
    let a = aniso(&gz, &tex, 8);
    r.coherent = a.share_hi;
    r.r8 = a.r8;
    r.r8_windows = a.windows;
    r.lambda_km = valley_spacing(&zm, &tex, n) * cell_km;
    // relief over 12.5 km blocks
    let b = ((12.5 / cell_km).round() as usize).max(2);
    let mut rel = Vec::new();
    for by in 0..n / b {
        'blk: for bx in 0..n / b {
            let (mut lo, mut hi) = (f32::MAX, f32::MIN);
            for y in by * b..(by + 1) * b {
                for x in bx * b..(bx + 1) * b {
                    let k = y * n + x;
                    if !tex[k] {
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
    r.deciles = (1..10).map(|d| pct(&alts, d as f64 / 10.0)).collect();
    // hydrology: trunks, drainage density, lakes (D5 apart)
    let hy = hydrology(&gz, cell_km);
    let mut chan = 0usize;
    for k in 0..n * n {
        if land[k] && hy.acc_km2[k] >= 100.0 {
            r.trunks.push((k % n, k / n));
        }
        if tex[k] && hy.strahler[k] > 0 {
            chan += 1;
        }
    }
    r.density = chan as f32 * cell_km / (r.textured.max(1) as f32 * cell_km * cell_km);
    let mut lake_d5 = vec![false; hy.lakes.len()];
    if let Some(m) = d5 {
        for k in 0..n * n {
            if hy.lake_id[k] != 0 && m[k] {
                lake_d5[hy.lake_id[k] as usize - 1] = true;
            }
        }
    }
    for (i, l) in hy.lakes.iter().enumerate() {
        let t = if lake_d5[i] { &mut r.lakes_d5 } else { &mut r.lakes_other };
        t.0 += 1;
        t.1 += l.area_km2;
    }
    r.octaves = octave_rms(&zm, n, r.textured);
    (r, hy)
}

pub fn trunk_dist(fine: &[(usize, usize)], coarse: &[(usize, usize)], ratio: f32, n_coarse: usize) -> (f32, f32, f32) {
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
pub fn drift(fine: &GridF32, coarse: &GridF32) -> (f32, f32) {
    let back = restrict(fine);
    let (ob, oc) = (ocean_mask(&back), ocean_mask(coarse));
    let (mut s, mut a, mut c) = (0f64, 0f64, 0usize);
    for k in 0..coarse.data.len() {
        if !ob[k] && !oc[k] {
            let d = (back.data[k] - coarse.data[k]) as f64;
            s += d;
            a += d.abs();
            c += 1;
        }
    }
    ((s / c.max(1) as f64) as f32, (a / c.max(1) as f64) as f32)
}

// ───────────────────────────── images ─────────────────────────────

pub fn rgba_of(z: &GridF32, hy: Option<&Hydro>, cell_m: f32, ss: &SteinSteinParams) -> Vec<u8> {
    let mut rgba = shade_m_rgba(z, cell_m, ss.depth_scale_m as f32);
    if let Some(h) = hy {
        draw(&mut rgba, h, 2, true, true);
    }
    rgba
}
/// A `side_frac` window (origin given as a fraction of the grid) of a north-up RGBA image, resized to `px`.
pub fn window(rgba: &[u8], n: usize, o: (f32, f32), frac: f32, px: usize) -> image::RgbaImage {
    let side = (n as f32 * frac) as usize;
    let (ox, oy) = ((o.0 * n as f32) as usize, (o.1 * n as f32) as usize);
    let mut img = image::RgbaImage::new(side as u32, side as u32);
    for r in 0..side {
        let fy = oy + side - 1 - r;
        let ir = n - 1 - fy;
        for c in 0..side {
            let off = (ir * n + ox + c) * 4;
            img.put_pixel(c as u32, r as u32, image::Rgba([rgba[off], rgba[off + 1], rgba[off + 2], 255]));
        }
    }
    image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Nearest)
}
pub fn whole(rgba: &[u8], n: usize, px: usize) -> image::RgbaImage {
    let img = image::RgbaImage::from_raw(n as u32, n as u32, rgba.to_vec()).unwrap();
    image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Triangle)
}
pub fn compose(rows: &[Vec<image::RgbaImage>]) -> image::RgbaImage {
    let gap = 6u32;
    let w: u32 = rows.iter().map(|r| r.iter().map(|t| t.width() + gap).sum::<u32>()).max().unwrap();
    let h: u32 = rows.iter().map(|r| r.iter().map(|t| t.height()).max().unwrap_or(0) + gap).sum();
    let mut out = image::RgbaImage::from_pixel(w, h, image::Rgba([24, 24, 24, 255]));
    let mut y = 0u32;
    for row in rows {
        let mut x = 0u32;
        for t in row {
            image::imageops::overlay(&mut out, t, x as i64, y as i64);
            x += t.width() + gap;
        }
        y += row.iter().map(|t| t.height()).max().unwrap_or(0) + gap;
    }
    out
}
/// A log-log line plot (wavelength km → RMS m) of named series.
pub fn spectra_png(series: &[(String, Vec<(f32, f32)>, [u8; 3])], path: &Path) {
    let (w, h) = (1000u32, 620u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let (x0, x1) = (0.2f32.ln(), 50f32.ln());
    let (y0, y1) = (0.5f32.ln(), 500f32.ln());
    let px = |lam: f32| (40.0 + (lam.ln() - x0) / (x1 - x0) * (w as f32 - 60.0)) as i64;
    let py = |v: f32| (h as f32 - 30.0 - (v.max(0.5).ln() - y0) / (y1 - y0) * (h as f32 - 50.0)) as i64;
    let mut put = |x: i64, y: i64, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as u32) < w && (y as u32) < h {
            img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    };
    for lam in [0.5f32, 1.0, 2.0, 5.0, 10.0, 20.0] {
        for y in 20..h as i64 - 30 {
            put(px(lam), y, [220, 220, 220]);
        }
    }
    for v in [1.0f32, 10.0, 100.0] {
        for x in 40..w as i64 - 20 {
            put(x, py(v), [220, 220, 220]);
        }
    }
    for (_, pts, c) in series {
        for win in pts.windows(2) {
            let (a, b) = (win[0], win[1]);
            let steps = 60;
            for s in 0..=steps {
                let t = s as f32 / steps as f32;
                let lam = (a.0.ln() * (1.0 - t) + b.0.ln() * t).exp();
                let v = (a.1.max(0.5).ln() * (1.0 - t) + b.1.max(0.5).ln() * t).exp();
                for d in -1..=1 {
                    put(px(lam), py(v) + d, *c);
                }
            }
        }
        for p in pts {
            for dx in -3..=3 {
                for dy in -3..=3 {
                    put(px(p.0) + dx, py(p.1) + dy, *c);
                }
            }
        }
    }
    img.save(path).unwrap();
}
pub fn hyps_png(curves: &[(String, Vec<f32>, [u8; 3])], path: &Path) {
    let (w, h) = (900u32, 520u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let y_of = |m: f32| (h as f32 - 20.0 - (m / 3500.0).clamp(0.0, 1.0) * (h as f32 - 40.0)) as u32;
    for x in 40..w - 20 {
        for m in [1000.0, 2000.0, 3000.0] {
            img.put_pixel(x, y_of(m), image::Rgba([215, 215, 215, 255]));
        }
        img.put_pixel(x, y_of(0.0), image::Rgba([0, 0, 0, 255]));
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

