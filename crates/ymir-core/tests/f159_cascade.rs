//! ADR Finding 159 — the multi-scale cascade prototype (64² → 128² → 256² → 512²) on the témoin, measured level by
//! level against the production témoin brought to the same resolutions. Declared in
//! `docs/reports/relief_method/f159_cascade/f159_declared.md` (D1–D9, R, I); this bench prints every instrument of D8,
//! the stop rule R, the cost extrapolation, and writes the images of I.
//!
//! The témoin (8192², ~4 min) is reduced ONCE to its block means and crops, cached under `target_f159/f159_cache`
//! (delete it to rebuild).
//!
//! Run: cargo test -p ymir-core --release --test f159_cascade -- --ignored --exact f159_cascade --nocapture

mod common;

use common::{CANONICAL_ORIGIN, Knobs, PSEED, aniso, build_field_seed, build_world, pct, sorted};
use std::path::{Path, PathBuf};
use std::time::Instant;
use ymir_core::cascade::{Cascade, CascadeConfig, LevelRecord, SubStep, diff_rgba, roll, shade_rgba};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction, carve, skeleton};
use ymir_core::terrain::flow::{FlowConfig, compute_flow};

const SEA: f32 = 0.5;
const S_EQ: f32 = 0.024;
const DOMAIN_KM: f32 = 400.0;
const LEVELS: [usize; 3] = [128, 256, 512];
/// F157's lake 4 (data origin of its 615-cell crop at 8192²) — crop C's centre is the crop's centre.
const LAKE4: (usize, usize) = (1913 + 307, 3277 + 307);

fn cache_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target_f159/f159_cache")
}
fn out_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/reports/relief_method/f159_cascade/images")
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

/// D9 — centred block means: level pixel I ← fine pixels f·I − f/2 … f·I + f/2 − 1, wrapped.
fn block_mean(g: &GridF32, n: usize) -> GridF32 {
    let f = g.width / n;
    let w = g.width;
    let mut out = GridF32::new(n, n, 0.0);
    for j in 0..n {
        for i in 0..n {
            let mut s = 0f64;
            for dy in 0..f {
                let y = (j * f + w + dy - f / 2) % w;
                for dx in 0..f {
                    let x = (i * f + w + dx - f / 2) % w;
                    s += g.data[y * w + x] as f64;
                }
            }
            out.data[j * n + i] = (s / (f * f) as f64) as f32;
        }
    }
    out
}

/// D6 — the centred full-weighting restriction (1/4, 1/2, 1/4) of a 2N field to N.
fn restrict(g: &GridF32) -> GridF32 {
    let (w2, n) = (g.width, g.width / 2);
    let wt = [0.25f32, 0.5, 0.25];
    let mut out = GridF32::new(n, n, 0.0);
    for j in 0..n {
        for i in 0..n {
            let mut s = 0f32;
            for (a, wy) in wt.iter().enumerate() {
                let y = (2 * j + w2 + a - 1) % w2;
                for (b, wx) in wt.iter().enumerate() {
                    let x = (2 * i + w2 + b - 1) % w2;
                    s += wy * wx * g.data[y * w2 + x];
                }
            }
            out.data[j * n + i] = s;
        }
    }
    out
}

fn n2m(ss: &SteinSteinParams) -> f32 {
    2.0 * 1.13 * ss.depth_scale_m as f32
}

/// The per-level reading of every D8 instrument on one field.
#[derive(Clone, Debug, Default)]
struct Reading {
    land_cells: usize,
    walls28: f64,
    walls33: f64,
    crest_share: f64,
    crest_curv_p50: f32,
    crest_curv_p90: f32,
    facet_share: f64,
    coherent_share: f32,
    r8: f32,
    r8_windows: usize,
    lambda_cells: f32,
    relief_p50: f32,
    mean_alt: f32,
    alt_p90: f32,
    trunks: Vec<(usize, usize)>,
}

/// Slope (deg) by central differences, as F157-B7.
fn slope_deg(zm: &[f32], w: usize, k: usize, cell_m: f32) -> f32 {
    let gx = (zm[k + 1] - zm[k - 1]) / (2.0 * cell_m);
    let gy = (zm[k + w] - zm[k - w]) / (2.0 * cell_m);
    (gx * gx + gy * gy).sqrt().atan().to_degrees()
}

/// D8.2 — crest cells: land, higher than both neighbours along at least one of the four axes.
fn crest_mask(zm: &[f32], land: &[bool], w: usize, h: usize) -> Vec<bool> {
    let mut c = vec![false; w * h];
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if !land[k] {
                continue;
            }
            let z = zm[k];
            let pairs = [(k - 1, k + 1), (k - w, k + w), (k - w - 1, k + w + 1), (k - w + 1, k + w - 1)];
            c[k] = pairs.iter().any(|&(a, b)| z > zm[a] && z > zm[b]);
        }
    }
    c
}

/// The 5-point Laplacian (m / km²) on interior cells.
fn laplacian(zm: &[f32], w: usize, h: usize, cell_km: f32) -> Vec<f32> {
    let mut l = vec![0f32; w * h];
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            l[k] = (zm[k - 1] + zm[k + 1] + zm[k - w] + zm[k + w] - 4.0 * zm[k]) / (cell_km * cell_km);
        }
    }
    l
}

/// D8.3 — the facet share: crest cells whose ring at Chebyshev distance 2 (non-crest land cells) has a median |∇²z|
/// below 0.25 × the crest's own |∇²z|.
fn facet_share(crest: &[bool], land: &[bool], lap: &[f32], w: usize, h: usize) -> f64 {
    let (mut crests, mut facets) = (0usize, 0usize);
    let mut ring = Vec::with_capacity(16);
    for y in 3..h - 3 {
        for x in 3..w - 3 {
            let k = y * w + x;
            if !crest[k] {
                continue;
            }
            crests += 1;
            ring.clear();
            for dy in -2i32..=2 {
                for dx in -2i32..=2 {
                    if dx.abs().max(dy.abs()) != 2 {
                        continue;
                    }
                    let j = ((y as i32 + dy) as usize) * w + (x as i32 + dx) as usize;
                    if land[j] && !crest[j] {
                        ring.push(lap[j].abs());
                    }
                }
            }
            if ring.is_empty() {
                continue;
            }
            ring.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if ring[ring.len() / 2] < 0.25 * lap[k].abs() {
                facets += 1;
            }
        }
    }
    100.0 * facets as f64 / crests.max(1) as f64
}

/// D8.7 — λ (cells): row and column land runs ≥ 16 cells, local minima with a 1-D prominence ≥ 1 m.
fn valley_spacing(zm: &[f32], land: &[bool], w: usize, h: usize) -> f32 {
    let (mut len, mut minima) = (0usize, 0usize);
    let mut line = |idx: &dyn Fn(usize) -> usize, n: usize| {
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
    };
    for y in 0..h {
        line(&|x| y * w + x, w);
    }
    for x in 0..w {
        line(&|y| y * w + x, h);
    }
    len as f32 / minima.max(1) as f32
}

fn read(z: &GridF32, ss: &SteinSteinParams) -> Reading {
    let (w, h) = (z.width, z.height);
    let cell_km = DOMAIN_KM / w as f32;
    let cell_m = cell_km * 1000.0;
    let zm: Vec<f32> = z.data.iter().map(|&v| (v - SEA) * n2m(ss)).collect();
    let land: Vec<bool> = z.data.iter().map(|&v| v > SEA).collect();
    let mut r = Reading { land_cells: land.iter().filter(|&&l| l).count(), ..Default::default() };
    // walls
    let (mut n, mut w28, mut w33) = (0usize, 0usize, 0usize);
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if !land[k] {
                continue;
            }
            n += 1;
            let s = slope_deg(&zm, w, k, cell_m);
            if (s - 28.0).abs() <= 0.5 {
                w28 += 1;
            }
            if (s - 33.0).abs() <= 0.5 {
                w33 += 1;
            }
        }
    }
    r.walls28 = 100.0 * w28 as f64 / n.max(1) as f64;
    r.walls33 = 100.0 * w33 as f64 / n.max(1) as f64;
    // crests and facets
    let crest = crest_mask(&zm, &land, w, h);
    let lap = laplacian(&zm, w, h, cell_km);
    let curv = sorted(crest.iter().zip(&lap).filter(|(c, _)| **c).map(|(_, l)| -l).collect());
    r.crest_share = 100.0 * curv.len() as f64 / r.land_cells.max(1) as f64;
    r.crest_curv_p50 = pct(&curv, 0.5);
    r.crest_curv_p90 = pct(&curv, 0.9);
    r.facet_share = facet_share(&crest, &land, &lap, w, h);
    // coherence and R8 (windows of 8 cells)
    let a = aniso(z, &land, 8);
    r.coherent_share = a.share_hi;
    r.r8 = a.r8;
    r.r8_windows = a.windows;
    r.lambda_cells = valley_spacing(&zm, &land, w, h);
    // relief over 12.5 km all-land blocks, the land hypsometry
    let b = w / 32;
    let mut rel = Vec::new();
    for by in 0..32 {
        'blk: for bx in 0..32 {
            let (mut lo, mut hi) = (f32::MAX, f32::MIN);
            for y in by * b..(by + 1) * b {
                for x in bx * b..(bx + 1) * b {
                    let k = y * w + x;
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
    r.mean_alt = alts.iter().sum::<f32>() / alts.len().max(1) as f32;
    r.alt_p90 = pct(&alts, 0.9);
    // trunks (D8 area ≥ 100 km²)
    let flow = compute_flow(z, &FlowConfig { sea_level: SEA, ..Default::default() });
    let c2 = cell_km * cell_km;
    for k in 0..w * h {
        if land[k] && flow.accumulation.data[k] * c2 >= 100.0 {
            r.trunks.push((k % w, k / w));
        }
    }
    r
}

/// D8.6 — distances (in cells of the COARSER grid) from each trunk cell of `fine` to the nearest of `coarse`, where a
/// fine pixel p sits at coarse coordinate p / ratio; wrapped. (p50, p90, max).
fn trunk_dist(fine: &[(usize, usize)], coarse: &[(usize, usize)], ratio: f32, n_coarse: usize) -> (f32, f32, f32) {
    if fine.is_empty() || coarse.is_empty() {
        return (f32::NAN, f32::NAN, f32::NAN);
    }
    let nc = n_coarse as f32;
    let mut d: Vec<f32> = fine
        .iter()
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

/// D8.3's control (rule 13): pyramids (planar faces) must read ≥ 60 %, paraboloid hills ≤ 20 %, at the level's grid.
fn facet_control(n: usize) -> (f64, f64) {
    let cell_km = DOMAIN_KM / n as f32;
    let (p, r) = (24usize, 11.0f32);
    let mut pyr = vec![0f32; n * n];
    let mut par = vec![0f32; n * n];
    for y in 0..n {
        for x in 0..n {
            // an oblique lattice (rotated 30°) so the faces are not on the axes
            let (fx, fy) = (x as f32, y as f32);
            let (u, v) = (fx * 0.866 + fy * 0.5, -fx * 0.5 + fy * 0.866);
            let (du, dv) = ((u.rem_euclid(p as f32)) - p as f32 / 2.0, (v.rem_euclid(p as f32)) - p as f32 / 2.0);
            pyr[y * n + x] = 100.0 + 800.0 * (1.0 - du.abs().max(dv.abs()) / r).max(0.0);
            par[y * n + x] = 100.0 + 800.0 * (1.0 - (du * du + dv * dv) / (r * r)).max(0.0);
        }
    }
    let land = vec![true; n * n];
    let f = |zm: &[f32]| {
        let crest = crest_mask(zm, &land, n, n);
        let lap = laplacian(zm, n, n, cell_km);
        facet_share(&crest, &land, &lap, n, n)
    };
    (f(&pyr), f(&par))
}

/// The windows of I (on the témoin's 512² block mean): (A range, B plateau, C basin + lake 4), as 512² origins.
fn crops(ref512: &GridF32, ss: &SteinSteinParams) -> [(usize, usize); 3] {
    let n = 512;
    let zm: Vec<f32> = ref512.data.iter().map(|&v| (v - SEA) * n2m(ss)).collect();
    let land: Vec<bool> = ref512.data.iter().map(|&v| v > SEA).collect();
    let alts = sorted(zm.iter().zip(&land).filter(|(_, l)| **l).map(|(v, _)| *v).collect());
    let p75 = pct(&alts, 0.75);
    let (side, stride, b) = (128usize, 16usize, 16usize);
    let p50 = pct(&alts, 0.5);
    let (mut best_a, mut best_b) = ((0usize, 0usize, f32::MIN), (0usize, 0usize, f32::MAX));
    // AMENDMENT after run 1 (not blind, images only): the declared B rule found no window, and the first run fell back
    // to (0, 0) silently -- open sea. The fallback is ≥ 75 % land with a mean above the land p50.
    let mut best_b2 = (0usize, 0usize, f32::MAX);
    for oy in (0..=n - side).step_by(stride) {
        for ox in (0..=n - side).step_by(stride) {
            let (mut nl, mut sum) = (0usize, 0f64);
            for y in oy..oy + side {
                for x in ox..ox + side {
                    if land[y * n + x] {
                        nl += 1;
                        sum += zm[y * n + x] as f64;
                    }
                }
            }
            let mut rel = Vec::new();
            for by in (oy..oy + side).step_by(b) {
                'blk: for bx in (ox..ox + side).step_by(b) {
                    let (mut lo, mut hi) = (f32::MAX, f32::MIN);
                    for y in by..by + b {
                        for x in bx..bx + b {
                            if !land[y * n + x] {
                                continue 'blk;
                            }
                            lo = lo.min(zm[y * n + x]);
                            hi = hi.max(zm[y * n + x]);
                        }
                    }
                    rel.push(hi - lo);
                }
            }
            if rel.is_empty() {
                continue;
            }
            let rp = pct(&sorted(rel), 0.5);
            let frac = nl as f32 / (side * side) as f32;
            if frac >= 0.5 && rp > best_a.2 {
                best_a = (ox, oy, rp);
            }
            if frac >= 0.9 && (sum / nl as f64) as f32 > p75 && rp < best_b.2 {
                best_b = (ox, oy, rp);
            }
            if frac >= 0.75 && (sum / nl as f64) as f32 > p50 && rp < best_b2.2 {
                best_b2 = (ox, oy, rp);
            }
        }
    }
    let c = ((LAKE4.0 / 16).saturating_sub(side / 2), (LAKE4.1 / 16).saturating_sub(side / 2));
    assert!(best_a.2 > f32::MIN, "crop A: no window");
    let b = if best_b.2 < f32::MAX {
        eprintln!("   crop B: the declared rule (≥ 90 % land, mean > p75) chose its window");
        (best_b.0, best_b.1)
    } else {
        assert!(best_b2.2 < f32::MAX, "crop B: no window under the amended rule either");
        eprintln!("   crop B: the declared rule (≥ 90 % land, mean > p75) finds NO window; the amended rule (≥ 75 % land, mean > p50) is used");
        (best_b2.0, best_b2.1)
    };
    [(best_a.0, best_a.1), b, c]
}

/// A crop of a level field (window given at 512², scaled to the level), shaded and resized (nearest) to `px`.
fn tile(g: &GridF32, o512: (usize, usize), px: usize, ss: &SteinSteinParams) -> image::RgbaImage {
    let n = g.width;
    let f = n as f32 / 512.0;
    let (ox, oy, side) = ((o512.0 as f32 * f) as usize, (o512.1 as f32 * f) as usize, (128.0 * f) as usize);
    let mut c = GridF32::new(side, side, 0.0);
    for y in 0..side {
        for x in 0..side {
            c.data[y * side + x] = g.data[(oy + y) * n + ox + x];
        }
    }
    let rgba = shade_rgba(&c, DOMAIN_KM / n as f32 * 1000.0, ss.depth_scale_m as f32);
    let img = image::RgbaImage::from_raw(side as u32, side as u32, rgba).unwrap();
    image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Nearest)
}

/// A crop of a difference field (m, level grid), saturating at the p98 of |Δ| on the crop; (tile, saturation).
fn diff_tile(d: &[f32], n: usize, o512: (usize, usize), px: usize) -> (image::RgbaImage, f32) {
    let f = n as f32 / 512.0;
    let (ox, oy, side) = ((o512.0 as f32 * f) as usize, (o512.1 as f32 * f) as usize, (128.0 * f) as usize);
    let mut c = vec![0f32; side * side];
    for y in 0..side {
        for x in 0..side {
            c[y * side + x] = d[(oy + y) * n + ox + x];
        }
    }
    let sat = pct(&sorted(c.iter().map(|v| v.abs()).collect()), 0.98).max(1.0);
    let img = image::RgbaImage::from_raw(side as u32, side as u32, diff_rgba(&c, side, side, sat)).unwrap();
    (image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Nearest), sat)
}

fn compose(tiles: &[Vec<image::RgbaImage>], px: usize) -> image::RgbaImage {
    let cols = tiles.iter().map(|r| r.len()).max().unwrap();
    let gap = 6;
    let mut out = image::RgbaImage::from_pixel(
        (cols * (px + gap)) as u32,
        (tiles.len() * (px + gap)) as u32,
        image::Rgba([24, 24, 24, 255]),
    );
    for (r, row) in tiles.iter().enumerate() {
        for (c, t) in row.iter().enumerate() {
            image::imageops::overlay(&mut out, t, (c * (px + gap)) as i64, (r * (px + gap)) as i64);
        }
    }
    out
}

#[test]
#[ignore]
fn f159_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    std::fs::create_dir_all(cache_dir()).unwrap();
    std::fs::create_dir_all(out_dir()).unwrap();
    eprintln!("\n==========  Finding 159 . the multi-scale cascade 64² → 512² on the témoin  ==========");
    let temoin = Knobs {
        valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))),
        slope_floor_abs: Some(S_EQ),
        ..Knobs::passes(2)
    };
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    // ── the témoin reduced to the levels (cached), the coarse 64² field, the témoin's 8192² crop tiles ──
    let px = 512usize;
    let cpath = cache_dir().join("coarse64.bin");
    let have = LEVELS.iter().all(|n| cache_dir().join(format!("ref{n}.bin")).exists())
        && cpath.exists()
        && (0..3).all(|i| cache_dir().join(format!("temoin8192_crop{i}.png")).exists());
    if !have {
        let t = Instant::now();
        let wd = build_world(temoin, None, PSEED, None);
        assert_eq!(wd.heightmap.width, 8192);
        let coarse = c1_coarse_normalized_altitude(
            &wd.state,
            &IsostasyConfig::c1_default(),
            &ss,
            wd.cfg.target_land_fraction,
        );
        save_grid(&cpath, &coarse);
        for n in LEVELS {
            save_grid(&cache_dir().join(format!("ref{n}.bin")), &block_mean(&wd.heightmap, n));
        }
        let ref512 = load_grid(&cache_dir().join("ref512.bin")).unwrap();
        for (i, o) in crops(&ref512, &ss).iter().enumerate() {
            // the 8192² crop, shown at px by a block mean of 2048 / px (declared in I)
            let (ox, oy) = (o.0 * 16, o.1 * 16);
            let mut c = GridF32::new(2048, 2048, 0.0);
            for y in 0..2048 {
                for x in 0..2048 {
                    c.data[y * 2048 + x] = wd.heightmap.data[(oy + y) * 8192 + ox + x];
                }
            }
            let rgba = shade_rgba(&c, DOMAIN_KM / 8192.0 * 1000.0, ss.depth_scale_m as f32);
            let img = image::RgbaImage::from_raw(2048, 2048, rgba).unwrap();
            image::imageops::resize(&img, px as u32, px as u32, image::imageops::FilterType::Triangle)
                .save(cache_dir().join(format!("temoin8192_crop{i}.png")))
                .unwrap();
        }
        drop(wd);
        eprintln!("   the témoin built and reduced ({:.0} s)", t.elapsed().as_secs_f64());
    }
    let refs: Vec<GridF32> =
        LEVELS.iter().map(|n| load_grid(&cache_dir().join(format!("ref{n}.bin"))).unwrap()).collect();
    let coarse = roll(&load_grid(&cpath).unwrap(), off);
    // ── D8.3 for information: S1 and the construction at 512² (cached) ──
    let (s1p, conp) = (cache_dir().join("s1_512.bin"), cache_dir().join("construction_512.bin"));
    if !(s1p.exists() && conp.exists()) {
        let t = Instant::now();
        let vc = ValleyConstruction::new(F121_AGE_K, Some(0.1));
        let s1 = build_field_seed(
            Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) },
            PSEED,
        );
        save_grid(&s1p, &block_mean(&s1, 512));
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (built, _) = carve(&s1, &sk, &vc, &ss);
        drop((sk, s1));
        save_grid(&conp, &block_mean(&built, 512));
        eprintln!("   S1 and the construction built and reduced ({:.0} s)", t.elapsed().as_secs_f64());
    }
    // ── the cascade ──
    let cfg = CascadeConfig::declared(ss.depth_scale_m as f32);
    eprintln!(
        "\n   D: K {:.1e} (code) · dt {:.0e} yr · α {} · D_phys {:.3} m²/yr · eq < {} % for {} steps · caps {:?}",
        cfg.k,
        cfg.dt_yr,
        cfg.alpha,
        cfg.d_phys_m2_yr,
        cfg.eq_ratio * 100.0,
        cfg.eq_window,
        cfg.max_steps
    );
    for &n in &cfg.levels {
        eprintln!(
            "      level {n}² · cell {:.4} km · D {:.3} m²/yr · L_c {:.3} cells · weight {:.4}",
            cfg.cell_km(n),
            cfg.diffusivity(n),
            (cfg.diffusivity(n) / (cfg.k / 1000.0)).sqrt() / (cfg.cell_km(n) * 1000.0),
            cfg.diffusion_weight(n)
        );
    }
    let t_c = Instant::now();
    let mut cas = Cascade::new(coarse.clone(), cfg.clone());
    let mut last_print = 0usize;
    cas.run_all(
        &mut |p| {
            if p.step == 1 || p.step % 50 == 0 || p.step == p.max_steps {
                eprintln!(
                    "      {} {}² step {} / {} · ratio {:.4}",
                    if p.level.is_none() { "calibration" } else { "level" },
                    p.n,
                    p.step,
                    p.max_steps,
                    p.ratio
                );
                last_print = p.step;
            }
        },
        &|| false,
    );
    let cascade_secs = t_c.elapsed().as_secs_f64();
    let cal = cas.calibration.clone().unwrap();
    eprintln!(
        "\n   calibration (128²): U₀ trial {:.1e} → mean land {:.0} m against h_iso's {:.0} m → U₀ = {:.3e} m/yr per km \
         · {} steps ({}) · {:.1} s · talus moved {:.2} % of the land at the last step",
        cal.u0_trial,
        cal.trial_mean_m,
        cal.target_mean_m,
        cal.u0,
        cal.steps,
        if cal.equilibrium { "equilibrium" } else { "CAP" },
        cal.secs,
        100.0 * cal.talus_last_share
    );
    // ── M: per level ──
    let control: Vec<(f64, f64)> = LEVELS.iter().map(|&n| facet_control(n)).collect();
    let cr: Vec<Reading> = cas.levels.iter().map(|l| read(&l.z, &ss)).collect();
    let rr: Vec<Reading> = refs.iter().map(|g| read(g, &ss)).collect();
    let mut stop_at: Option<(usize, String)> = None;
    // the cost extrapolation (D8.9)
    let l512 = cas.levels.last().unwrap();
    let t_step = (l512.secs[1] + l512.secs[2] + l512.secs[3]) / l512.steps.max(1) as f64;
    let nl = |n: f64| n * n * (n * n).ln();
    let extra: f64 = [1024f64, 2048.0, 4096.0, 8192.0]
        .iter()
        .map(|&n| t_step * nl(n) / nl(512.0) * l512.steps as f64)
        .sum();
    let total_est = extra + cascade_secs;
    eprintln!(
        "\n   cost: 128² → 512² {:.1} s (calibration included) · t_step(512²) {:.3} s × {} steps · ESTIMATE to 8192² {:.0} s \
         = {:.1} min (15 min target, 1 h stop)",
        cascade_secs,
        t_step,
        l512.steps,
        total_est,
        total_est / 60.0
    );
    for (i, l) in cas.levels.iter().enumerate() {
        let (c, r) = (&cr[i], &rr[i]);
        let (pyr, par) = control[i];
        let ctrl_ok = pyr >= 60.0 && par <= 20.0;
        eprintln!(
            "\n   ── level {}² (cell {:.3} km) · {} steps ({}, last ratio {:.4}) · {:.1} s [upscale {:.2} · uplift {:.2} · \
             erosion {:.1} · diffusion {:.1}] · retained {:.1} MB · talus moved {:.2} % of the land at the last step",
            l.n,
            l.cell_km,
            l.steps,
            if l.equilibrium { "equilibrium" } else { "CAP" },
            l.ratios.last().copied().unwrap_or(f32::NAN),
            l.total_secs(),
            l.secs[0],
            l.secs[1],
            l.secs[2],
            l.secs[3],
            l.retained_bytes() as f64 / 1e6,
            100.0 * l.talus_last_share
        );
        let row = |name: &str, a: String, b: String| eprintln!("      {name:<34} cascade {a:>14} · témoin {b:>14}");
        row("land cells", c.land_cells.to_string(), r.land_cells.to_string());
        row("planar walls 28° ± 0.5° (%)", format!("{:.3}", c.walls28), format!("{:.3}", r.walls28));
        row("planar band 33° ± 0.5° (%)", format!("{:.3}", c.walls33), format!("{:.3}", r.walls33));
        row("crest cells (% land)", format!("{:.2}", c.crest_share), format!("{:.2}", r.crest_share));
        row(
            "crest −∇²z p50 / p90 (m/km²)",
            format!("{:.1} / {:.1}", c.crest_curv_p50, c.crest_curv_p90),
            format!("{:.1} / {:.1}", r.crest_curv_p50, r.crest_curv_p90),
        );
        row("crest facets (% crests)", format!("{:.2}", c.facet_share), format!("{:.2}", r.facet_share));
        eprintln!(
            "      facet control: pyramids {:.1} % (≥ 60) · paraboloids {:.1} % (≤ 20) → {}",
            pyr,
            par,
            if ctrl_ok { "PASSES" } else { "FAILS: the facet instrument is blind at this level" }
        );
        row("coherent windows C > 0.7 (%)", format!("{:.2}", c.coherent_share), format!("{:.2}", r.coherent_share));
        row(
            "terrain R8 (windows; floor 1/√N)",
            format!("{:.4} ({})", c.r8, c.r8_windows),
            format!("{:.4} ({})", r.r8, r.r8_windows),
        );
        eprintln!("      R8 noise floors: cascade {:.4} · témoin {:.4}", 1.0 / (c.r8_windows.max(1) as f32).sqrt(), 1.0 / (r.r8_windows.max(1) as f32).sqrt());
        let lc_cells = (cfg.diffusivity(l.n) / (cfg.k / 1000.0)).sqrt() / (l.cell_km * 1000.0);
        row("valley spacing λ (cells)", format!("{:.2}", c.lambda_cells), format!("{:.2}", r.lambda_cells));
        eprintln!(
            "      Perron's band 6.4–12.7 L_c = {:.2}–{:.2} cells (L_c {:.3} cells); transects read high by ≈ π/2",
            6.4 * lc_cells,
            12.7 * lc_cells,
            lc_cells
        );
        row("relief p50, 12.5 km blocks (m)", format!("{:.0}", c.relief_p50), format!("{:.0}", r.relief_p50));
        row("mean land / p90 (m)", format!("{:.0} / {:.0}", c.mean_alt, c.alt_p90), format!("{:.0} / {:.0}", r.mean_alt, r.alt_p90));
        row("trunk cells ≥ 100 km²", c.trunks.len().to_string(), r.trunks.len().to_string());
        let (a, b, m) = trunk_dist(&c.trunks, &r.trunks, 1.0, l.n);
        eprintln!(
            "      cascade trunks → témoin trunks (same level): p50 {:.2} / p90 {:.2} / max {:.2} cells ({:.2} / {:.2} / {:.2} km)",
            a, b, m, a * l.cell_km, b * l.cell_km, m * l.cell_km
        );
        let mut p90_prev = 0f32;
        if i > 0 {
            let prev = &cas.levels[i - 1];
            let (a, b, m) = trunk_dist(&c.trunks, &cr[i - 1].trunks, 2.0, prev.n);
            let (ra, rb, rm) = trunk_dist(&cr[i - 1].trunks, &c.trunks, 0.5, l.n);
            p90_prev = b;
            eprintln!(
                "      trunks {}² → {}² (cells of {}²): p50 {:.2} / p90 {:.2} / max {:.2} ({:.2} / {:.2} / {:.2} km) · reverse \
                 (cells of {}²) {:.2} / {:.2} / {:.2}",
                l.n, prev.n, prev.n, a, b, m, a * prev.cell_km, b * prev.cell_km, m * prev.cell_km, l.n, ra, rb, rm
            );
            let (rt_a, rt_b, rt_m) = trunk_dist(&rr[i].trunks, &rr[i - 1].trunks, 2.0, prev.n);
            eprintln!("      the témoin's own trunks {}² → {}²: p50 {:.2} / p90 {:.2} / max {:.2} cells", l.n, prev.n, rt_a, rt_b, rt_m);
            // D6: the coarse level kept? (full weighting)
            let back = restrict(&l.z);
            let (mut s, mut sa, mut nn) = (0f64, 0f64, 0usize);
            for k in 0..prev.n * prev.n {
                if prev.z.data[k] > SEA && back.data[k] > SEA {
                    let d = ((back.data[k] - prev.z.data[k]) * n2m(&ss)) as f64;
                    s += d;
                    sa += d.abs();
                    nn += 1;
                }
            }
            eprintln!(
                "      D6 the coarse level: {}² restricted − {}² on land: bias {:+.1} m · mean |Δ| {:.1} m ({} cells)",
                l.n, prev.n, s / nn.max(1) as f64, sa / nn.max(1) as f64, nn
            );
            let rel = (c.relief_p50 - cr[i - 1].relief_p50) / cr[i - 1].relief_p50.max(1.0) * 100.0;
            eprintln!("      relief p50 change {}² → {}²: {:+.1} %", prev.n, l.n, rel);
        }
        // R
        if stop_at.is_none() {
            let mut why = Vec::new();
            if c.walls28 > r.walls28 {
                why.push(format!("walls 28° {:.3} > {:.3}", c.walls28, r.walls28));
            }
            if c.walls33 > r.walls33 {
                why.push(format!("band 33° {:.3} > {:.3}", c.walls33, r.walls33));
            }
            if ctrl_ok && c.facet_share > r.facet_share {
                why.push(format!("facets {:.2} > {:.2}", c.facet_share, r.facet_share));
            }
            if i > 0 && p90_prev > 2.0 {
                why.push(format!("trunk p90 {:.2} > 2 cells", p90_prev));
            }
            if total_est > 3600.0 {
                why.push(format!("cost ESTIMATE {:.0} s > 1 h", total_est));
            }
            if !why.is_empty() {
                stop_at = Some((l.n, why.join("; ")));
            }
        }
    }
    let (s1, con) = (load_grid(&s1p).unwrap(), load_grid(&conp).unwrap());
    let (rs1, rcon) = (read(&s1, &ss), read(&con, &ss));
    eprintln!(
        "\n   D8.3 for information at 512² (block means): crest facets S1 {:.2} % · construction {:.2} % · témoin {:.2} % · \
         cascade {:.2} %; walls 28° S1 {:.3} · construction {:.3} %",
        rs1.facet_share, rcon.facet_share, rr[2].facet_share, cr[2].facet_share, rs1.walls28, rcon.walls28
    );
    match &stop_at {
        Some((n, why)) => eprintln!("\n   R — THE STOP RULE FIRES at {n}²: {why}"),
        None => eprintln!("\n   R — the stop rule does not fire at any level: the cascade reaches 512²"),
    }
    // ── I: the images ──
    let ref512 = &refs[2];
    let cs = crops(ref512, &ss);
    let names = ["A_chaine", "B_plateau", "C_bassin_lac4"];
    for (ci, o) in cs.iter().enumerate() {
        let tem8192 = image::open(cache_dir().join(format!("temoin8192_crop{ci}.png"))).unwrap().to_rgba8();
        let mut rows: Vec<Vec<image::RgbaImage>> = Vec::new();
        rows.push(cas.levels.iter().map(|l| tile(&l.z, *o, px, &ss)).chain(std::iter::once(tem8192)).collect());
        rows.push(refs.iter().map(|g| tile(g, *o, px, &ss)).collect());
        let last = cas.levels.last().unwrap();
        rows.push(SubStep::ALL.iter().map(|&s| tile(&last.substep_field(s), *o, px, &ss)).collect());
        // AMENDMENT after run 1 (not blind, images only): the cumulated sub-step fields are unreadable (the uplift alone
        // stacks kilometres), so a 4th row shows what each process DID at 512²: the level's total change (z − its
        // agrandissement), then Σ uplift, Σ erosion, Σ diffusion (m; red raised, blue lowered; saturation = the p98 of
        // |Δ| on the crop, printed)
        let n2 = n2m(&ss);
        let tot: Vec<f32> = last.z.data.iter().zip(&last.upscaled.data).map(|(a, b)| (a - b) * n2).collect();
        let parts: [Vec<f32>; 4] = [
            tot,
            last.sum_uplift.iter().map(|v| v * n2).collect(),
            last.sum_erosion.iter().map(|v| v * n2).collect(),
            last.sum_diffusion.iter().map(|v| v * n2).collect(),
        ];
        let mut sats = Vec::new();
        rows.push(
            parts
                .iter()
                .map(|d| {
                    let (t, s) = diff_tile(d, last.n, *o, px);
                    sats.push(s);
                    t
                })
                .collect(),
        );
        eprintln!("   image {} row 4 saturations (m): total {:.0} · uplift {:.0} · erosion {:.0} · diffusion {:.0}", names[ci], sats[0], sats[1], sats[2], sats[3]);
        compose(&rows, px).save(out_dir().join(format!("f159_{}.png", names[ci]))).unwrap();
        eprintln!(
            "   image {} · window 512² origin ({}, {}) = 8192² ({}, {}), 100 km, north up",
            names[ci],
            o.0,
            o.1,
            o.0 * 16,
            o.1 * 16
        );
    }
    // the whole map per level, and the difference with the previous level (Δ m, saturating at 300 m)
    for (i, l) in cas.levels.iter().enumerate() {
        let rgba = shade_rgba(&l.z, l.cell_km * 1000.0, ss.depth_scale_m as f32);
        image::RgbaImage::from_raw(l.n as u32, l.n as u32, rgba)
            .unwrap()
            .save(out_dir().join(format!("f159_level{}.png", l.n)))
            .unwrap();
        let prev_up = &l.upscaled;
        let d: Vec<f32> = l.z.data.iter().zip(&prev_up.data).map(|(a, b)| (a - b) * n2m(&ss)).collect();
        image::RgbaImage::from_raw(l.n as u32, l.n as u32, diff_rgba(&d, l.n, l.n, 300.0))
            .unwrap()
            .save(out_dir().join(format!("f159_level{}_diff.png", l.n)))
            .unwrap();
        let _ = i;
    }
    let _ = last_print;
    let _: &[LevelRecord] = &cas.levels;
    eprintln!("\n==========  end Finding 159 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
