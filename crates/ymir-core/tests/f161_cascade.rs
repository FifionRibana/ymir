//! ADR Finding 161 — the cascade, round 3: Corsica (Copernicus DEM GLO-30) as the texture reference; the physics level
//! at 256²; the N1 amplification (the implicit solver at n = 1, the area capped, the smooth retargeting after each
//! erosion block, k calibrated per level on Corsica's 2–4 cell octave), its variants (sans plafond, sans recalage,
//! talus au quart), measured level by level against Corsica at the same cell (and the témoin for information).
//! Declared in `docs/reports/relief_method/f161_corsica/f161_declared.md`.
//!
//! Needs `data/corsica/corse_<n>.bin` (`prep_corse.py`). The témoin is reduced once, cached under `target_f161`.
//!
//! Run: cargo test -p ymir-core --release --test f161_cascade -- --ignored --exact f161_cascade --nocapture

mod common;

use common::{CANONICAL_ORIGIN, Knobs, PSEED, aniso, build_world, pct, sorted};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, ocean_mask, restrict};
use ymir_core::cascade::hydro::{Hydro, draw, hydrology};
use ymir_core::cascade::measure::{octave_km, octave_rms};
use ymir_core::cascade::{CascadeConfig, physics_level, roll, shade_m_rgba};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const DOMAIN_KM: f32 = 400.0;
const PEAK_TARGET_M: f32 = 2850.0;
const PEAK_OK: (f32, f32) = (0.85 * 2700.0, 1.15 * 3000.0);
const LEVELS: [usize; 3] = [512, 1024, 2048];
const CROPS: [(&str, usize, usize); 3] = [("A_chaine", 128, 320), ("B_plateau", 176, 192), ("C_bassin_lac4", 74, 160)];
/// A 100 km Corsica window (256² origin of its 200 km grid): the Cinto–Rotondo spine and its valleys.
const CORSE_CROP: (usize, usize) = (72, 96);
const NOTICE: &str = "produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH \
2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn cache_dir() -> PathBuf {
    root().join("target_f161/f161_cache")
}
fn out_dir() -> PathBuf {
    root().join("docs/reports/relief_method/f161_corsica/images")
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

// ───────────────────────────── the instruments (F160's, in metres; F161-C3) ─────────────────────────────

#[derive(Clone, Debug, Default)]
struct Reading {
    n: usize,
    land: usize,
    textured: usize,
    slope_p50: f32,
    slope_p90: f32,
    slope28: f64,
    slope33: f64,
    walls28: f64,
    walls28_cells: usize,
    crest_share: f64,
    facets: f64,
    facet_cells: usize,
    coherent: f32,
    r8: f32,
    r8_windows: usize,
    lambda_km: f32,
    density: f32,
    relief_p50: f32,
    mean: f32,
    p90: f32,
    peak: f32,
    deciles: Vec<f32>,
    trunks: Vec<(usize, usize)>,
    octaves: Vec<f32>,
    lakes_d5: (usize, f32),
    lakes_other: (usize, f32),
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
            [(k - 1, k + 1), (k - n, k + n), (k - n - 1, k + n + 1), (k - n + 1, k + n - 1)]
                .iter()
                .any(|&(a, b)| z > zm[a] && z > zm[b])
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
fn facet_control(n: usize, cell_km: f32) -> (f64, f64) {
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

/// F161-C3 on a field in metres (the sea ≤ 0 m connected to the border); `d5` excluded from the texture.
fn read(z: &GridF32, d5: Option<&[bool]>) -> (Reading, Hydro) {
    let n = z.width;
    let cell_km = DOMAIN_KM / n as f32 * if z.width == z.height { 1.0 } else { 1.0 };
    read_cell(z, d5, cell_km)
}
fn read_cell(z: &GridF32, d5: Option<&[bool]>, cell_km: f32) -> (Reading, Hydro) {
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
fn drift(fine: &GridF32, coarse: &GridF32) -> (f32, f32) {
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

fn rgba_of(z: &GridF32, hy: Option<&Hydro>, cell_m: f32, ss: &SteinSteinParams) -> Vec<u8> {
    let mut rgba = shade_m_rgba(z, cell_m, ss.depth_scale_m as f32);
    if let Some(h) = hy {
        draw(&mut rgba, h, 2, true, true);
    }
    rgba
}
/// A `side_frac` window (origin given as a fraction of the grid) of a north-up RGBA image, resized to `px`.
fn window(rgba: &[u8], n: usize, o: (f32, f32), frac: f32, px: usize) -> image::RgbaImage {
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
fn whole(rgba: &[u8], n: usize, px: usize) -> image::RgbaImage {
    let img = image::RgbaImage::from_raw(n as u32, n as u32, rgba.to_vec()).unwrap();
    image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Triangle)
}
fn compose(rows: &[Vec<image::RgbaImage>]) -> image::RgbaImage {
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
fn spectra_png(series: &[(String, Vec<(f32, f32)>, [u8; 3])], path: &Path) {
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
fn hyps_png(curves: &[(String, Vec<f32>, [u8; 3])], path: &Path) {
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

// ───────────────────────────── the run ─────────────────────────────

fn row(name: &str, a: String, b: String, c: String) {
    eprintln!("      {name:<34} cascade {a:>16} · Corse {b:>16} · témoin {c:>16}");
}
fn print_level(tag: &str, r: &Reading, co: &Reading, te: Option<&Reading>, ctrl: (f64, f64)) {
    let t = |f: &dyn Fn(&Reading) -> String| te.map_or("—".to_string(), f);
    let cell_km = DOMAIN_KM / r.n as f32;
    eprintln!("   ── {tag} ({}², cell {:.0} m) · land {} (textured {}) · Corse land {}", r.n, cell_km * 1000.0, r.land, r.textured, co.land);
    row("slope p50 / p90 (m/m)", format!("{:.3} / {:.3}", r.slope_p50, r.slope_p90), format!("{:.3} / {:.3}", co.slope_p50, co.slope_p90), t(&|x| format!("{:.3} / {:.3}", x.slope_p50, x.slope_p90)));
    row("slopes > 28° / > 33° (%)", format!("{:.2} / {:.2}", r.slope28, r.slope33), format!("{:.2} / {:.2}", co.slope28, co.slope33), t(&|x| format!("{:.2} / {:.2}", x.slope28, x.slope33)));
    row("walls 28° ± 0.5° (% / cells)", format!("{:.3} / {}", r.walls28, r.walls28_cells), format!("{:.3} / {}", co.walls28, co.walls28_cells), t(&|x| format!("{:.3} / {}", x.walls28, x.walls28_cells)));
    row("crest facets (% / crests)", format!("{:.2} / {}", r.facets, r.facet_cells), format!("{:.2} / {}", co.facets, co.facet_cells), t(&|x| format!("{:.2} / {}", x.facets, x.facet_cells)));
    eprintln!("      facet control at this grid: pyramids {:.1} % · paraboloids {:.1} %", ctrl.0, ctrl.1);
    row("coherent windows C > 0.7 (%)", format!("{:.2}", r.coherent), format!("{:.2}", co.coherent), t(&|x| format!("{:.2}", x.coherent)));
    row("terrain R8 (windows)", format!("{:.4} ({})", r.r8, r.r8_windows), format!("{:.4} ({})", co.r8, co.r8_windows), t(&|x| format!("{:.4} ({})", x.r8, x.r8_windows)));
    row("valley spacing λ (km)", format!("{:.2}", r.lambda_km), format!("{:.2}", co.lambda_km), t(&|x| format!("{:.2}", x.lambda_km)));
    row("drainage density (km/km²)", format!("{:.3}", r.density), format!("{:.3}", co.density), t(&|x| format!("{:.3}", x.density)));
    row("relief p50, 12.5 km blocks (m)", format!("{:.0}", r.relief_p50), format!("{:.0}", co.relief_p50), t(&|x| format!("{:.0}", x.relief_p50)));
    row("peak / mean / p90 land (m)", format!("{:.0} / {:.0} / {:.0}", r.peak, r.mean, r.p90), format!("{:.0} / {:.0} / {:.0}", co.peak, co.mean, co.p90), t(&|x| format!("{:.0} / {:.0} / {:.0}", x.peak, x.mean, x.p90)));
    eprintln!(
        "      hypsometry deciles (m) cascade {} · Corse {}",
        r.deciles.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join("/"),
        co.deciles.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join("/")
    );
    eprintln!(
        "      lakes: D5 {} / {:.0} km² (apart) · other {} / {:.0} km² · Corse {} / {:.0} km²",
        r.lakes_d5.0, r.lakes_d5.1, r.lakes_other.0, r.lakes_other.1, co.lakes_other.0, co.lakes_other.1
    );
    let oct = |x: &Reading, j: usize| x.octaves.get(j).map_or("—".to_string(), |v| format!("{v:.1}"));
    eprintln!(
        "      spectrum by octave (RMS m), cascade / Corse / témoin: {}",
        (0..r.octaves.len().min(6))
            .map(|j| {
                let (a, b) = octave_km(j, cell_km);
                format!("[{a:.2}–{b:.2} km] {} / {} / {}", oct(r, j), oct(co, j), te.map_or("—".to_string(), |x| oct(x, j)))
            })
            .collect::<Vec<_>>()
            .join(" · ")
    );
}

/// R against Corsica (F161-R): the reasons it fires at this level, empty if none.
fn stop_reasons(r: &Reading, co: &Reading, ctrl_ok: bool, drift_rel: f32, trunk_p90: f32) -> Vec<String> {
    let mut why = Vec::new();
    if ctrl_ok && r.facets > 1.5 * co.facets && r.facet_cells > 20 {
        why.push(format!("facets {:.2} > 1.5 × {:.2}", r.facets, co.facets));
    }
    if r.walls28 > 1.5 * co.walls28 && r.walls28_cells > 20 {
        why.push(format!("walls 28° {:.3} > 1.5 × {:.3}", r.walls28, co.walls28));
    }
    let out = |a: f32, b: f32| a > 1.5 * b || a < b / 1.5;
    if out(r.slope_p50, co.slope_p50) {
        why.push(format!("slope p50 {:.3} vs {:.3}", r.slope_p50, co.slope_p50));
    }
    if out(r.slope_p90, co.slope_p90) {
        why.push(format!("slope p90 {:.3} vs {:.3}", r.slope_p90, co.slope_p90));
    }
    if out(r.lambda_km, co.lambda_km) {
        why.push(format!("λ {:.2} vs {:.2} km", r.lambda_km, co.lambda_km));
    }
    if drift_rel.abs() > 0.02 {
        why.push(format!("drift {:+.2} %", 100.0 * drift_rel));
    }
    if r.peak < PEAK_OK.0 || r.peak > PEAK_OK.1 {
        why.push(format!("peak {:.0} m", r.peak));
    }
    if trunk_p90 > 1.5 {
        why.push(format!("trunk p90 {:.2}", trunk_p90));
    }
    why
}

struct Lvl {
    n: usize,
    result: GridF32,
    d5: Vec<bool>,
    reading: Reading,
    secs: [f64; 4],
    iters: [usize; 3],
    k: f32,
    drift: (f32, f32, f32),
    trunk: (f32, f32, f32),
    recal_rms: f32,
}

#[allow(clippy::too_many_arguments)]
fn run_chain(
    name: &str,
    mut ch: Chain,
    phys_reading: &Reading,
    corse: &BTreeMap<usize, Reading>,
    temoin: &BTreeMap<usize, Reading>,
    controls: &BTreeMap<usize, (f64, f64)>,
    calibrate: bool,
) -> (Vec<Lvl>, Option<(usize, String)>, Chain, Vec<(usize, Vec<(f32, f32)>, f32)>) {
    let mut out: Vec<Lvl> = Vec::new();
    let mut stop = None;
    let mut cal_log = Vec::new();
    let mut prev_reading = phys_reading.clone();
    let mut prev_result = ch.physics_m.clone();
    for n in LEVELS {
        let co = &corse[&(n / 2)];
        if calibrate {
            let t = Instant::now();
            let target = co.octaves[0];
            let (k, trials) = ch.calibrate_next(target, 0.2, 8, &|| false);
            eprintln!(
                "      calibration {n}²: target {:.2} m (Corse octave 0) → k {:.3e} after {} trials ({:.0} s): {}",
                target,
                k,
                trials.len(),
                t.elapsed().as_secs_f64(),
                trials.iter().map(|(k, v)| format!("{k:.2e}→{v:.2}")).collect::<Vec<_>>().join(", ")
            );
            cal_log.push((n, trials, target));
        }
        let l = ch.next_level(&|| false).unwrap().clone();
        let result = l.result();
        let (r, _) = read(&result, Some(&l.d5));
        let (b, a) = drift(&result, &prev_result);
        let rel = b / prev_reading.mean.max(1.0);
        let tp = trunk_dist(&r.trunks, &prev_reading.trunks, 2.0, n / 2);
        let recal_rms = (l.sum_recalage.iter().map(|v| (*v as f64) * (*v as f64)).sum::<f64>() / (n * n) as f64).sqrt() as f32;
        let ctrl = controls[&n];
        let why = stop_reasons(&r, co, ctrl.0 >= 60.0 && ctrl.1 <= 20.0, rel, tp.1);
        eprintln!(
            "      {name} {n}² ({:.1} s: E {:.1} · T {:.1} · D {:.1}) · k {:.3e} · peak {:.0} · drift {:+.2} % · trunk p90 {:.2}{}",
            l.secs.iter().sum::<f64>(),
            l.secs[1],
            l.secs[2],
            l.secs[3],
            l.k,
            r.peak,
            100.0 * rel,
            tp.1,
            if why.is_empty() { String::new() } else { format!(" · R: {}", why.join("; ")) }
        );
        if stop.is_none() && !why.is_empty() {
            stop = Some((n, why.join("; ")));
        }
        let _ = temoin;
        out.push(Lvl {
            n,
            result: result.clone(),
            d5: l.d5.clone(),
            reading: r.clone(),
            secs: l.secs,
            iters: [l.done.erosion, l.done.talus, l.done.deposit],
            k: l.k,
            drift: (b, a, rel),
            trunk: tp,
            recal_rms,
        });
        prev_reading = r;
        prev_result = result;
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    (out, stop, ch, cal_log)
}

fn print_chain(name: &str, lv: &[Lvl], stop: &Option<(usize, String)>, corse: &BTreeMap<usize, Reading>, temoin: &BTreeMap<usize, Reading>, controls: &BTreeMap<usize, (f64, f64)>) {
    eprintln!("\n   ═══ chain {name} ═══");
    for l in lv {
        eprintln!(
            "   {}²: E {} · T {} · D {} iterations · {:.1} s [upscale {:.2} · E {:.1} ({:.0} ms/it) · T {:.1} ({:.2} ms/it) · D {:.1} ({:.2} ms/it)] · k {:.3e} · drift {:+.1} m ({:+.2} %) |Δ| {:.1} m · trunk p50/p90/max {:.2}/{:.2}/{:.2} · smooth retargeting RMS {:.1} m",
            l.n,
            l.iters[0],
            l.iters[1],
            l.iters[2],
            l.secs.iter().sum::<f64>(),
            l.secs[0],
            l.secs[1],
            1e3 * l.secs[1] / l.iters[0].max(1) as f64,
            l.secs[2],
            1e3 * l.secs[2] / l.iters[1].max(1) as f64,
            l.secs[3],
            1e3 * l.secs[3] / l.iters[2].max(1) as f64,
            l.k,
            l.drift.0,
            100.0 * l.drift.2,
            l.drift.1,
            l.trunk.0,
            l.trunk.1,
            l.trunk.2,
            l.recal_rms
        );
        print_level(name, &l.reading, &corse[&(l.n / 2)], temoin.get(&l.n), controls[&l.n]);
    }
    match stop {
        Some((n, why)) => eprintln!("   R — {name} FIRES at {n}²: {why}"),
        None => eprintln!("   R — {name} does not fire up to 2048²"),
    }
}

#[test]
#[ignore]
fn f161_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    std::fs::create_dir_all(cache_dir()).unwrap();
    std::fs::create_dir_all(out_dir()).unwrap();
    eprintln!("\n==========  Finding 161 . the cascade, round 3: Corsica as the reference, N1 calibrated on it  ==========");
    // ── C: Corsica ──
    let t = Instant::now();
    let mut corse: BTreeMap<usize, Reading> = BTreeMap::new();
    let mut corse_grids: BTreeMap<usize, GridF32> = BTreeMap::new();
    for nc in [128usize, 256, 512, 1024, 2048] {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).expect("data/corsica (prep_corse.py)");
        let (r, _) = read_cell(&g, None, 200.0 / nc as f32);
        corse.insert(nc, r);
        corse_grids.insert(nc, g);
    }
    eprintln!("   Corsica read at 1563 / 781 / 391 / 195 / 98 m ({:.0} s)", t.elapsed().as_secs_f64());
    for nc in [128usize, 256, 512, 1024, 2048] {
        let r = &corse[&nc];
        let cell_km = 200.0 / nc as f32;
        eprintln!(
            "   Corse {:.0} m: land {} · slope p50/p90 {:.3}/{:.3} · > 28° {:.2} % · > 33° {:.2} % · walls 28° {:.3} % ({}) · facets {:.2} % ({}) · R8 {:.4} ({}) · coherent {:.2} % · λ {:.2} km · density {:.3} km/km² · relief p50 {:.0} m · peak/mean/p90 {:.0}/{:.0}/{:.0} · lakes {} / {:.1} km²",
            cell_km * 1000.0, r.land, r.slope_p50, r.slope_p90, r.slope28, r.slope33, r.walls28, r.walls28_cells, r.facets, r.facet_cells, r.r8, r.r8_windows, r.coherent, r.lambda_km, r.density, r.relief_p50, r.peak, r.mean, r.p90, r.lakes_other.0, r.lakes_other.1
        );
        eprintln!(
            "      spectrum (RMS m): {}",
            r.octaves.iter().enumerate().map(|(j, v)| { let (a, b) = octave_km(j, cell_km); format!("[{a:.2}–{b:.2} km] {v:.1}") }).collect::<Vec<_>>().join(" · ")
        );
    }
    // ── the témoin, at the levels ──
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let cpath = cache_dir().join("coarse64.bin");
    let tlev = [256usize, 512, 1024, 2048];
    if !(cpath.exists() && tlev.iter().all(|n| cache_dir().join(format!("ref{n}.bin")).exists())) {
        let t = Instant::now();
        let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
        let wd = build_world(temoin, None, PSEED, None);
        let coarse = c1_coarse_normalized_altitude(&wd.state, &IsostasyConfig::c1_default(), &ss, wd.cfg.target_land_fraction);
        save_grid(&cpath, &coarse);
        for n in tlev {
            save_grid(&cache_dir().join(format!("ref{n}.bin")), &block_mean(&wd.heightmap, n));
        }
        drop(wd);
        eprintln!("   the témoin built and reduced ({:.0} s)", t.elapsed().as_secs_f64());
    }
    let coarse = roll(&load_grid(&cpath).unwrap(), off);
    let mut temoin: BTreeMap<usize, Reading> = BTreeMap::new();
    for n in tlev {
        let g = to_m(&load_grid(&cache_dir().join(format!("ref{n}.bin"))).unwrap(), n2m);
        temoin.insert(n, read(&g, None).0);
    }
    let controls: BTreeMap<usize, (f64, f64)> = [256usize, 512, 1024, 2048].iter().map(|&n| (n, facet_control(n, DOMAIN_KM / n as f32))).collect();
    // ── D1: the physics level at 256² ──
    let ccfg = CascadeConfig::declared(ss.depth_scale_m as f32);
    let (rec, cal) = physics_level(&coarse, 256, &ccfg, PEAK_TARGET_M, 300, &mut |_| {}, &|| false).unwrap();
    let pm = to_m(&rec.z, n2m);
    let (pr, _) = read(&pm, Some(&cal.lifted_mask));
    eprintln!(
        "\n   physics P256: U₀ {:.3e} (trial peak {:.0} m) · {} steps ({}) · D5 lifted {} cells (mean depth {:.0} m)",
        cal.u0,
        cal.trial_peak_m,
        rec.steps,
        if rec.equilibrium { "equilibrium" } else { "CAP" },
        cal.lifted,
        cal.lifted_mean_depth_m
    );
    print_level("physics 256²", &pr, &corse[&128], temoin.get(&256), controls[&256]);
    // ── N1, calibrated ──
    let mk = |v: Variant| Chain::new(pm.clone(), AmpConfig::declared(PSEED, v, 1.0)).with_d5(cal.lifted_mask.clone());
    let t = Instant::now();
    let (n1, n1_stop, n1_chain, cal_log) = run_chain("N1", mk(Variant::N1), &pr, &corse, &temoin, &controls, true);
    let n1_secs = t.elapsed().as_secs_f64();
    let ks: Vec<(usize, f32)> = LEVELS.iter().map(|&n| (n, n1_chain.cfg.n1_k_at(n))).collect();
    eprintln!("   N1 calibrated k: {} ({:.0} s including the calibration)", ks.iter().map(|(n, k)| format!("{n}² {k:.3e}")).collect::<Vec<_>>().join(" · "), n1_secs);
    let with_k = |v: Variant, tal: f32| {
        let mut c = mk(v);
        c.cfg.n1_k = ks.clone();
        c.cfg.talus_fine_scale = tal;
        c
    };
    let (nc, nc_stop, _, _) = run_chain("N1 sans plafond", with_k(Variant::N1NoCap, 1.0), &pr, &corse, &temoin, &controls, false);
    let (nr, nr_stop, _, _) = run_chain("N1 sans recalage", with_k(Variant::N1NoRecal, 1.0), &pr, &corse, &temoin, &controls, false);
    let (nt, nt_stop, _, _) = run_chain("N1 talus au quart", with_k(Variant::N1, 0.25), &pr, &corse, &temoin, &controls, false);
    for (name, lv, st) in [("N1", &n1, &n1_stop), ("N1 sans plafond", &nc, &nc_stop), ("N1 sans recalage", &nr, &nr_stop), ("N1 talus au quart", &nt, &nt_stop)] {
        print_chain(name, lv, st, &corse, &temoin, &controls);
    }
    // the variants against N1
    eprintln!("\n   the variants against N1, per level:");
    for i in 0..LEVELS.len() {
        let (a, c, r, q) = (&n1[i], &nc[i], &nr[i], &nt[i]);
        eprintln!(
            "      {}² · peak N1 {:.0} / sans recalage {:.0} · drift N1 {:+.2} % / sans recalage {:+.2} % · trunk p90 N1 {:.2} / sans plafond {:.2} · slope p50/p90 N1 {:.3}/{:.3} / talus ¼ {:.3}/{:.3} ({:+.1} % / {:+.1} %) · talus time {:.1} s / {:.1} s",
            a.n,
            a.reading.peak,
            r.reading.peak,
            100.0 * a.drift.2,
            100.0 * r.drift.2,
            a.trunk.1,
            c.trunk.1,
            a.reading.slope_p50,
            a.reading.slope_p90,
            q.reading.slope_p50,
            q.reading.slope_p90,
            100.0 * (q.reading.slope_p50 / a.reading.slope_p50 - 1.0),
            100.0 * (q.reading.slope_p90 / a.reading.slope_p90 - 1.0),
            a.secs[2],
            q.secs[2]
        );
    }
    // the cost (N1, measured; ESTIMATE per process to 8 192²: the per-iteration time at 2 048² × the cells × F160's budgets)
    let l = &n1[2];
    let per = [l.secs[1] / l.iters[0].max(1) as f64, l.secs[2] / l.iters[1].max(1) as f64, l.secs[3] / l.iters[2].max(1) as f64];
    let measured: f64 = n1.iter().map(|l| l.secs.iter().sum::<f64>()).sum();
    let mut est = [0f64; 3];
    for nn in [4096usize, 8192] {
        let sc = (nn as f64 / 2048.0).powi(2) * (nn as f64).ln() / 2048f64.ln();
        let b = n1_chain.cfg.budget(nn);
        est[0] += per[0] * sc * b.erosion as f64;
        est[1] += per[1] * (nn as f64 / 2048.0).powi(2) * b.talus as f64;
        est[2] += per[2] * (nn as f64 / 2048.0).powi(2) * b.deposit as f64;
    }
    let tot = measured + est.iter().sum::<f64>();
    eprintln!(
        "\n   cost: N1 512² → 2048² {:.0} s without the calibration ({:.0} s with it) · per iteration at 2048²: E {:.0} ms · T {:.1} ms · D {:.1} ms · ESTIMATE 4096² + 8192²: E {:.0} s · T {:.0} s · D {:.0} s → total {:.0} s = {:.1} min (alert above 1 h: {})",
        measured,
        n1_secs,
        1e3 * per[0],
        1e3 * per[1],
        1e3 * per[2],
        est[0],
        est[1],
        est[2],
        tot,
        tot / 60.0,
        if tot > 3600.0 { "YES" } else { "no" }
    );
    // ── I: the images ──
    let mut views: Vec<(usize, GridF32, Option<Vec<bool>>)> = vec![(256, pm.clone(), Some(cal.lifted_mask.clone()))];
    for l in &n1 {
        views.push((l.n, l.result.clone(), Some(l.d5.clone())));
    }
    for (n, g, _) in &views {
        let cell_m = DOMAIN_KM / *n as f32 * 1000.0;
        let (_, hy) = read(g, None);
        let a = whole(&rgba_of(g, Some(&hy), cell_m, &ss), *n, 1024);
        let cg = &corse_grids[&(n / 2)];
        let (_, chy) = read_cell(cg, None, 200.0 / (n / 2) as f32);
        let b = whole(&rgba_of(cg, Some(&chy), cell_m, &ss), n / 2, 512);
        compose(&[vec![a, b]]).save(out_dir().join(format!("f161_N1_vs_Corse_{n}.png"))).unwrap();
    }
    eprintln!("   images f161_N1_vs_Corse_<n>.png: N1 (400 km, 1024 px) beside Corsica (200 km, 512 px), same km per pixel · {NOTICE}");
    // crops: row 1 N1 per level, row 2 Corsica's spine window at the same cell, row 3 N1 sans recalage
    for (name, ox, oy) in CROPS {
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        let mut r3 = Vec::new();
        for (i, (n, g, _)) in views.iter().enumerate() {
            let cell_m = DOMAIN_KM / *n as f32 * 1000.0;
            let (_, hy) = read(g, None);
            r1.push(window(&rgba_of(g, Some(&hy), cell_m, &ss), *n, (ox as f32 / 512.0, oy as f32 / 512.0), 0.25, 512));
            let cg = &corse_grids[&(n / 2)];
            let (_, chy) = read_cell(cg, None, 200.0 / (n / 2) as f32);
            r2.push(window(&rgba_of(cg, Some(&chy), cell_m, &ss), n / 2, (CORSE_CROP.0 as f32 / 256.0, CORSE_CROP.1 as f32 / 256.0), 0.5, 512));
            if i > 0 {
                let gr = &nr[i - 1].result;
                let (_, hr) = read(gr, None);
                r3.push(window(&rgba_of(gr, Some(&hr), cell_m, &ss), *n, (ox as f32 / 512.0, oy as f32 / 512.0), 0.25, 512));
            }
        }
        compose(&[r1, r2, r3]).save(out_dir().join(format!("f161_{name}.png"))).unwrap();
    }
    eprintln!("   images f161_<crop>.png: row 1 N1 (physics 256², 512², 1024², 2048²), row 2 Corsica's spine window (100 km) at the same cell, row 3 N1 sans recalage (512²–2048²)");
    // the spectra
    let mut series: Vec<(String, Vec<(f32, f32)>, [u8; 3])> = Vec::new();
    let pts = |r: &Reading, cell_km: f32| -> Vec<(f32, f32)> {
        r.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, cell_km); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect()
    };
    let cols = [[200, 60, 40], [220, 140, 30], [60, 150, 60], [40, 110, 200]];
    for (i, (n, g, d5)) in views.iter().enumerate() {
        let (r, _) = read(g, d5.as_deref());
        series.push((format!("N1 {n}²"), pts(&r, DOMAIN_KM / *n as f32), cols[i % 4]));
    }
    for nc in [128usize, 256, 512, 1024] {
        series.push((format!("Corse {nc}"), pts(&corse[&nc], 200.0 / nc as f32), [0, 0, 0]));
    }
    for n in [512usize, 1024, 2048] {
        series.push((format!("témoin {n}²"), pts(&temoin[&n], DOMAIN_KM / n as f32), [150, 150, 150]));
    }
    spectra_png(&series, &out_dir().join("f161_spectra.png"));
    eprintln!("   image f161_spectra.png: RMS (m) against the octave's central wavelength (km), log-log; N1 in colour (red physics, orange 512², green 1024², blue 2048²), Corsica black, the témoin grey · {NOTICE}");
    let mut curves: Vec<(String, Vec<f32>, [u8; 3])> = Vec::new();
    for (i, (n, g, _)) in views.iter().enumerate() {
        let land = ocean_mask(g);
        curves.push((format!("N1 {n}²"), sorted(g.data.iter().zip(&land).filter(|(_, o)| !**o).map(|(v, _)| *v).collect()), cols[i % 4]));
    }
    let cg = &corse_grids[&1024];
    let cl = ocean_mask(cg);
    curves.push(("Corse 195 m".into(), sorted(cg.data.iter().zip(&cl).filter(|(_, o)| !**o).map(|(v, _)| *v).collect()), [0, 0, 0]));
    hyps_png(&curves, &out_dir().join("f161_hypsometry.png"));
    eprintln!("   image f161_hypsometry.png: land altitude against the share of land above it; N1 per level in colour, Corsica (195 m) black · {NOTICE}");
    let _ = &cal_log;
    eprintln!("\n==========  end Finding 161 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
