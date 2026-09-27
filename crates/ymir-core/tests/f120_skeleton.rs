//! ADR 0001 Finding 120 — **the skeleton a valley construction would stand on: is it tectonic, and
//! which order deserves a valley?** Measurement only; no primitive is written here.
//!
//! ## A0 — do the trunks move between A1+B2 OFF and ON?
//!
//! The construction proposed by the author lays valleys along the EXISTING trunks. If those trunks
//! are the incision's own product, the construction inherits the process it replaces. So the two
//! worlds' trunks (Strahler ≥ 4, then ≥ 3) are laid over each other cell by cell: the share of
//! trunk cells common, the deviation (chessboard distance to the other world's nearest trunk
//! cell), and the mouths.
//!
//! ## A — where does the comb start?
//!
//! Per Strahler order, on the A1+B2 ON world: segment count, total length, **R8 network** at chord
//! 16 and 32 (Finding 111's instrument, isotropic floor 0.0402), **R8 terrain** in a ±16-cell band
//! around the order's segments (Finding 97's `aniso`, window 16), and **σ of the slope** of the
//! hillslope cells that drain FIRST into a segment of that order (D8 attribution).
//!
//! ## Rule 12, declared
//!
//! * **The network** is the exported one: `assemble_hd_drainage` on the breached field, head
//!   extension at `A_c`, the chain of Findings 110/111/119. It is read under BOTH `full_tree`
//!   settings: `false` is the viz/export layer, and Finding 110 recorded that it omits ~86 % of the
//!   tree, so its Strahler order ranks EXTRACTED segments; `true` is the hierarchy.
//! * **Only `Watercourse` segments carry an order**: Finding 110 recorded that `strahler_order` is
//!   meaningless on a spillway (the kind is authoritative) and that `clip_rivers_to_lakes` copies a
//!   parent's order onto every fragment.
//! * **The attribution flow** is the assembled drainage's own D8 direction on the breached field.
//!
//! Run: cargo test -p ymir-core --release --test f120_skeleton -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, pct, sorted};
use std::collections::VecDeque;
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, SegmentKind, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const CAP: u16 = 1000;
/// Finding 101's tile of maximum Δ(R8): the delivered field COMBED there (R8 0.1190 vs 0.0250).
const COMB_TILE: (usize, usize, usize) = (2048, 5120, 1024);

/// Circular statistics of axial directions (mod π): (R2, R8) — Finding 111's `axial`.
fn axial(thetas: &[f32]) -> (f32, f32) {
    if thetas.is_empty() {
        return (0.0, 0.0);
    }
    let (mut c2, mut s2, mut c8, mut s8) = (0f64, 0f64, 0f64, 0f64);
    for &t in thetas {
        let t = t as f64;
        c2 += (2.0 * t).cos();
        s2 += (2.0 * t).sin();
        c8 += (8.0 * t).cos();
        s8 += (8.0 * t).sin();
    }
    let n = thetas.len() as f64;
    (((c2 * c2 + s2 * s2).sqrt() / n) as f32, ((c8 * c8 + s8 * s8).sqrt() / n) as f32)
}

/// Chord directions of a polyline at scale `l` cells — Finding 111's `chords`.
fn chords(points: &[(u32, u32)], l: usize) -> Vec<f32> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i + l < points.len() {
        let (dx, dy) = (
            points[i + l].0 as f32 - points[i].0 as f32,
            points[i + l].1 as f32 - points[i].1 as f32,
        );
        if dx != 0.0 || dy != 0.0 {
            let mut t = dy.atan2(dx);
            if t < 0.0 {
                t += std::f32::consts::PI;
            }
            out.push(t);
        }
        i += l.max(1);
    }
    out
}

/// Multi-source BFS, 8-neighbour (chessboard) distance on the torus, capped.
fn dist_from(src: &[bool], w: usize, h: usize) -> Vec<u16> {
    let n = w * h;
    let mut d = vec![CAP; n];
    let mut q = VecDeque::new();
    for (k, &s) in src.iter().enumerate() {
        if s {
            d[k] = 0;
            q.push_back(k);
        }
    }
    while let Some(k) = q.pop_front() {
        if d[k] + 1 >= CAP {
            continue;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nk = (y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize;
                if d[nk] > d[k] + 1 {
                    d[nk] = d[k] + 1;
                    q.push_back(nk);
                }
            }
        }
    }
    d
}

/// The exported network of an eroded field, and the breached field it lives on.
fn network(f: &GridF32, full: bool, ss: &SteinSteinParams) -> (C1DrainageResult, GridF32) {
    let (w, h) = (f.width, f.height);
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = full;
    let pre = c1_drainage_windowed(f, None, &dcfg, ss, DOMAIN_KM);
    let bf = breach_monotone(f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let cl = c1_climate_placed(&bf, ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
    let dr =
        assemble_hd_drainage(&bf, &dc, Some(pre), &dcfg, ss, DOMAIN_KM, GEO_RATIO, None, false)
            .drainage;
    (dr, bf)
}

/// Per cell: the highest Watercourse order passing through it (0 = none).
fn order_map(dr: &C1DrainageResult, n: usize, w: usize) -> Vec<u8> {
    let mut o = vec![0u8; n];
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] != SegmentKind::Watercourse {
            continue;
        }
        for &(x, y) in &s.points {
            let k = y as usize * w + x as usize;
            o[k] = o[k].max(s.strahler_order);
        }
    }
    o
}

/// Mouth cells of Watercourse segments of order ≥ k that end next to the sea.
fn mouths(dr: &C1DrainageResult, bf: &GridF32, k: u8) -> Vec<usize> {
    let (w, h) = (bf.width, bf.height);
    let mut out = Vec::new();
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] != SegmentKind::Watercourse
            || s.strahler_order < k
            || s.downstream.is_some()
        {
            continue;
        }
        let Some(&(x, y)) = s.points.last() else { continue };
        let (x, y) = (x as i32, y as i32);
        let wet = (-1i32..=1).any(|dy| {
            (-1i32..=1).any(|dx| {
                bf.data[(y + dy).rem_euclid(h as i32) as usize * w
                    + (x + dx).rem_euclid(w as i32) as usize]
                    <= SEA
            })
        });
        if wet {
            out.push(y as usize * w + x as usize);
        }
    }
    out
}

fn p(v: &[f32], q: f64) -> f32 {
    if v.is_empty() { f32::NAN } else { pct(v, q) }
}

#[test]
#[ignore]
fn f120_skeleton() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 120 . the skeleton: A0 (tectonic?) and A (threshold)  =====");

    let off = build_field_seed(Knobs::passes(2), PSEED);
    let on = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (w, h) = (on.width, on.height);
    let n = w * h;

    for full in [false, true] {
        let (d_off, b_off) = network(&off, full, &ss);
        let (d_on, b_on) = network(&on, full, &ss);
        let (o_off, o_on) = (order_map(&d_off, n, w), order_map(&d_on, n, w));
        let max_o = |dr: &C1DrainageResult| {
            dr.rivers
                .segments
                .iter()
                .enumerate()
                .filter(|(i, _)| dr.segment_kind[*i] == SegmentKind::Watercourse)
                .map(|(_, s)| s.strahler_order)
                .max()
                .unwrap_or(0)
        };
        eprintln!(
            "\n   ════════ full_tree = {full} ({}) · max order OFF {} · ON {} ════════",
            if full { "the HIERARCHY" } else { "the viz/export layer" },
            max_o(&d_off),
            max_o(&d_on)
        );

        // ── A0 · the trunks laid over each other ──
        eprintln!(
            "   ── A0 · trunk cells OFF vs ON (chessboard cells; cap {CAP}) ──\n   {:>6} {:>9} {:>9} \
             {:>8} {:>8} {:>8} {:>8}   {:>16} {:>16}   {:>22}",
            "order≥",
            "OFF cells",
            "ON cells",
            "exact",
            "≤1",
            "≤2",
            "≤5",
            "OFF→ON p50/p90",
            "ON→OFF p50/p90",
            "mouths OFF/ON · ≤2 cells"
        );
        for k in [5u8, 4, 3] {
            let t_off: Vec<bool> = o_off.iter().map(|&o| o >= k).collect();
            let t_on: Vec<bool> = o_on.iter().map(|&o| o >= k).collect();
            let (c_off, c_on) =
                (t_off.iter().filter(|&&b| b).count(), t_on.iter().filter(|&&b| b).count());
            if c_off == 0 || c_on == 0 {
                eprintln!("   {k:>6} {c_off:>9} {c_on:>9}   (empty in one world)");
                continue;
            }
            let dist_on = dist_from(&t_on, w, h);
            let dist_off = dist_from(&t_off, w, h);
            let dv: Vec<f32> = (0..n).filter(|&c| t_off[c]).map(|c| dist_on[c] as f32).collect();
            let dv2: Vec<f32> = (0..n).filter(|&c| t_on[c]).map(|c| dist_off[c] as f32).collect();
            let share = |d: u16| {
                100.0 * dv.iter().filter(|&&v| v <= d as f32).count() as f32 / dv.len() as f32
            };
            let (s1, s2) = (sorted(dv.clone()), sorted(dv2.clone()));
            let (m_off, m_on) = (mouths(&d_off, &b_off, k), mouths(&d_on, &b_on, k));
            let mut msrc = vec![false; n];
            for &c in &m_on {
                msrc[c] = true;
            }
            let md = dist_from(&msrc, w, h);
            let same = m_off.iter().filter(|&&c| md[c] <= 2).count();
            let mdv = sorted(m_off.iter().map(|&c| md[c] as f32).collect());
            eprintln!(
                "   {k:>6} {c_off:>9} {c_on:>9} {:>7.1}% {:>7.1}% {:>7.1}% {:>7.1}%   {:>7.0} / {:<7.0} \
                 {:>7.0} / {:<7.0}   {:>4} / {:<4} · {same:>3} (moved p50 {:.0}, p90 {:.0})",
                share(0),
                share(1),
                share(2),
                share(5),
                p(&s1, 0.50),
                p(&s1, 0.90),
                p(&s2, 0.50),
                p(&s2, 0.90),
                m_off.len(),
                m_on.len(),
                p(&mdv, 0.50),
                p(&mdv, 0.90)
            );
        }

        // ── A · per order, ON world ──
        let dr = &d_on;
        let bf = &b_on;
        let om = &o_on;
        let land: Vec<bool> = bf.data.iter().map(|&v| v > SEA).collect();
        // slope, degrees, central differences
        let slope = |k: usize| -> f32 {
            let (x, y) = (k % w, k / w);
            let at = |xx: usize, yy: usize| bf.data[yy * w + xx];
            let (xm, xp) = ((x + w - 1) % w, (x + 1) % w);
            let (ym, yp) = ((y + h - 1) % h, (y + 1) % h);
            let gx = (at(xp, y) - at(xm, y)) * n2m / (2.0 * CELL_M);
            let gy = (at(x, yp) - at(x, ym)) * n2m / (2.0 * CELL_M);
            (gx * gx + gy * gy).sqrt().atan().to_degrees()
        };
        // D8 attribution: each hillslope cell → the order of the first river cell downstream
        let mut att = vec![255u8; n];
        for c in 0..n {
            if !land[c] {
                att[c] = 0;
            } else if om[c] > 0 {
                att[c] = om[c];
            }
        }
        let mut path = Vec::new();
        for s in 0..n {
            if att[s] != 255 {
                continue;
            }
            path.clear();
            let mut cur = s;
            let val;
            loop {
                if att[cur] != 255 {
                    val = if att[cur] == 254 { 0 } else { att[cur] };
                    break;
                }
                att[cur] = 254;
                path.push(cur);
                let d = dr.flow.direction[cur];
                if d == DIR_NONE {
                    val = 0;
                    break;
                }
                let nx = ((cur % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
                let ny = ((cur / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
                cur = ny * w + nx;
            }
            for &c in &path {
                att[c] = val;
            }
        }
        let mo = max_o(dr);
        eprintln!(
            "\n   ── A · per Strahler order, A1+B2 ON (R8 isotropic floor 0.0402, Finding 97) ──\n   \
             {:>5} {:>7} {:>10} {:>10} {:>10} {:>10} {:>10} {:>12} {:>12}",
            "order",
            "segs",
            "length km",
            "R8net c16",
            "R8net c32",
            "R2net c32",
            "R8 terrain",
            "hill σ(°)",
            "hill cells"
        );
        for o in 1..=mo {
            let mut th16 = Vec::new();
            let mut th32 = Vec::new();
            let (mut segs, mut len) = (0usize, 0f64);
            let mut band = vec![false; n];
            for (i, s) in dr.rivers.segments.iter().enumerate() {
                if dr.segment_kind[i] != SegmentKind::Watercourse || s.strahler_order != o {
                    continue;
                }
                segs += 1;
                th16.extend(chords(&s.points, 16));
                th32.extend(chords(&s.points, 32));
                for win in s.points.windows(2) {
                    let diag = win[0].0 != win[1].0 && win[0].1 != win[1].1;
                    len += if diag {
                        CELL_KM as f64 * std::f64::consts::SQRT_2
                    } else {
                        CELL_KM as f64
                    };
                }
                for &(x, y) in &s.points {
                    band[y as usize * w + x as usize] = true;
                }
            }
            let bd = dist_from(&band, w, h);
            let bmask: Vec<bool> = (0..n).map(|c| land[c] && bd[c] <= 16).collect();
            drop(bd);
            let a = aniso(bf, &bmask, 16);
            let (_, r8_16) = axial(&th16);
            let (r2_32, r8_32) = axial(&th32);
            // σ of the slope of the hillslope cells attributed to this order
            let (mut cnt, mut mean, mut m2) = (0f64, 0f64, 0f64);
            for c in 0..n {
                if land[c] && om[c] == 0 && att[c] == o {
                    let v = slope(c) as f64;
                    cnt += 1.0;
                    let dlt = v - mean;
                    mean += dlt / cnt;
                    m2 += dlt * (v - mean);
                }
            }
            let sd = if cnt > 1.0 { (m2 / (cnt - 1.0)).sqrt() } else { f64::NAN };
            eprintln!(
                "   {o:>5} {segs:>7} {len:>10.0} {r8_16:>10.4} {r8_32:>10.4} {r2_32:>10.4} {:>10.4} \
                 {:>6.2} (μ{:>5.1}) {:>12.0}",
                a.r8, sd, mean, cnt
            );
        }

        // ── controls: F87-A, and the comb tile ──
        if let Some(i) = (0..dr.rivers.segments.len())
            .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway)
            .max_by(|&a, &b| dr.segment_discharge_m3s[a].total_cmp(&dr.segment_discharge_m3s[b]))
        {
            let s = &dr.rivers.segments[i];
            let cells: Vec<usize> =
                s.points.iter().map(|&(x, y)| y as usize * w + x as usize).collect();
            let adj = cells.iter().map(|&c| om[c]).max().unwrap_or(0);
            // the Watercourse feeding into / out of it
            let up = s.upstream.iter().map(|&u| dr.rivers.segments[u].strahler_order).max();
            let down = s.downstream.map(|d| dr.rivers.segments[d].strahler_order);
            eprintln!(
                "\n   CONTROL · F87-A (highest-discharge spillway): **{:.2} m³/s**, {} cells · its own \
                 `strahler_order` {} (meaningless on a spillway) · Watercourse order ON ITS CELLS {} · \
                 upstream max {:?} · downstream {:?}",
                dr.segment_discharge_m3s[i],
                s.points.len(),
                s.strahler_order,
                adj,
                up,
                down
            );
        }
        let (tx, ty, ts) = COMB_TILE;
        let mut by = std::collections::BTreeMap::<u8, (usize, Vec<f32>)>::new();
        for (i, s) in dr.rivers.segments.iter().enumerate() {
            if dr.segment_kind[i] != SegmentKind::Watercourse {
                continue;
            }
            let inside: Vec<(u32, u32)> = s
                .points
                .iter()
                .copied()
                .filter(|&(x, y)| {
                    (x as usize) >= tx
                        && (x as usize) < tx + ts
                        && (y as usize) >= ty
                        && (y as usize) < ty + ts
                })
                .collect();
            if inside.is_empty() {
                continue;
            }
            let e = by.entry(s.strahler_order).or_insert((0, Vec::new()));
            e.0 += 1;
            e.1.extend(chords(&inside, 32));
        }
        eprintln!(
            "   CONTROL · comb tile ({tx}, {ty}) {ts}² — segments inside, by order: {}",
            by.iter()
                .map(|(o, (c, th))| format!(
                    "o{o} {c} seg (R8 c32 {:.3}, n {})",
                    axial(th).1,
                    th.len()
                ))
                .collect::<Vec<_>>()
                .join(" · ")
        );
        eprintln!("   ⏱ {:.1} s", t0.elapsed().as_secs_f64());
    }
    eprintln!(
        "\n==========  end Finding 120 skeleton . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
