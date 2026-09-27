//! ADR 0001 Finding 120 — **A0 and A in PHYSICAL units** (the author's amendment, 2026-09-25):
//! the skeleton and the valley threshold read by DRAINED AREA in km², never by Strahler order.
//!
//! *"L'ordre dépend du seuil en cellules, donc de la résolution"* (Findings 57–67). The first A0
//! (`f120_skeleton`) confirmed it on its own table: the A1+B2 ON world carries **2.7×** the
//! order ≥ 4 cells of OFF (6 850 against 2 512), because the closure raised the orders (Finding 110,
//! +80 % at order 4), so "order ≥ 4" names different RIVERS in the two worlds and a cell overlap
//! mixes membership with displacement. A drained area in km² names the same river in both.
//!
//! ## A0 — the trunks as `A ≥ A_t` on the D8 accumulation of the BREACHED field
//!
//! Three worlds on seed 1, 8192²: **PRE** (the tectonic + FBM field, no incision: the only skeleton a
//! construction replacing the heavy incision could stand on without running the pipeline twice —
//! declared extension of the round), **OFF** (delivered), **ON** (A1+B2). Each is breached as the
//! export is, and its own `compute_flow` accumulation (in cells, NOT multiplied by the geo ratio —
//! Finding 68: the ratio multiplies four SEGMENT arrays only) is thresholded in km² in-domain.
//! Pairs: OFF vs ON (the round's) and PRE vs ON (the construction's).
//!
//! ## A — per band of drained area, ON, the full exported tree
//!
//! River points of the exported network (`full_tree = true`: the only setting that carries the
//! heads below the extraction threshold) binned by the D8 area at the point. Per band: points,
//! length, runs, **R8 network** at chords of 16 and 32 cells (0.78 / 1.56 km), **R8 terrain** in a
//! ±16-cell band, and σ of the slope of the hillslope cells that drain FIRST into a river point of
//! that band. Strahler order reported beside, never used.
//!
//! Every threshold declares its cell count at 8192² (0.002 384 km² per cell); below 40 cells a
//! band is SUB-GRID (the seven-instance pattern) and is marked.
//!
//! Run: cargo test -p ymir-core --release --test f120_area -- --ignored --nocapture

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
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowResult, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const CAP: u16 = 1000;
const COMB_TILE: (usize, usize, usize) = (2048, 5120, 1024);
/// A0 trunk thresholds, km² in-domain.
const A_T: [f32; 4] = [1.0, 10.0, 100.0, 1000.0];
/// A bands, km² in-domain: [lo, hi).
const BANDS: [(f32, f32); 5] =
    [(0.1, 1.0), (1.0, 10.0), (10.0, 100.0), (100.0, 1000.0), (1000.0, f32::INFINITY)];

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

fn p(v: &[f32], q: f64) -> f32 {
    if v.is_empty() { f32::NAN } else { pct(v, q) }
}

fn dcfg(full: bool) -> C1DrainageConfig {
    let mut d = C1DrainageConfig::default();
    d.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    d.thresholds.full_tree = full;
    d
}

/// The breached field of an eroded one, and the D8 flow on it.
fn breached(f: &GridF32, ss: &SteinSteinParams) -> (GridF32, FlowResult) {
    let (w, h) = (f.width, f.height);
    let pre = c1_drainage_windowed(f, None, &dcfg(false), ss, DOMAIN_KM);
    let bf = breach_monotone(f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let flow = c1_drainage_windowed(&bf, None, &dcfg(false), ss, DOMAIN_KM).flow;
    (bf, flow)
}

/// Mouths of trunks ≥ `a_cells`: trunk cells whose D8 receiver is sea.
fn mouths(flow: &FlowResult, bf: &GridF32, a_cells: f32) -> Vec<usize> {
    let (w, h) = (bf.width, bf.height);
    let mut out = Vec::new();
    for k in 0..w * h {
        if bf.data[k] <= SEA || flow.accumulation.data[k] < a_cells {
            continue;
        }
        let d = flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        if bf.data[ny * w + nx] <= SEA {
            out.push(k);
        }
    }
    out
}

fn a0_pair(
    name: &str,
    (fa, ba): (&FlowResult, &GridF32),
    (fb, bb): (&FlowResult, &GridF32),
    cell_km2: f32,
) {
    let (w, h) = (ba.width, ba.height);
    let n = w * h;
    eprintln!(
        "\n   ── A0 · {name} · trunk = land cells with D8 area ≥ A_t (chessboard cells, cap {CAP}) ──\
         \n   {:>8} {:>7} {:>9} {:>9} {:>7} {:>7} {:>7} {:>7}   {:>14} {:>14}   {}",
        "A_t km²",
        "cells",
        "A cells",
        "B cells",
        "exact",
        "≤1",
        "≤2",
        "≤5",
        "A→B p50/p90",
        "B→A p50/p90",
        "mouths A/B · matched ≤2 cells"
    );
    for a_t in A_T {
        let ac = a_t / cell_km2;
        let ta: Vec<bool> =
            (0..n).map(|k| ba.data[k] > SEA && fa.accumulation.data[k] >= ac).collect();
        let tb: Vec<bool> =
            (0..n).map(|k| bb.data[k] > SEA && fb.accumulation.data[k] >= ac).collect();
        let (ca, cb) = (ta.iter().filter(|&&b| b).count(), tb.iter().filter(|&&b| b).count());
        let db = dist_from(&tb, w, h);
        let da = dist_from(&ta, w, h);
        let v1: Vec<f32> = (0..n).filter(|&k| ta[k]).map(|k| db[k] as f32).collect();
        let v2: Vec<f32> = (0..n).filter(|&k| tb[k]).map(|k| da[k] as f32).collect();
        let share = |d: u16| {
            100.0 * v1.iter().filter(|&&v| v <= d as f32).count() as f32 / v1.len().max(1) as f32
        };
        let (s1, s2) = (sorted(v1.clone()), sorted(v2.clone()));
        let (ma, mb) = (mouths(fa, ba, ac), mouths(fb, bb, ac));
        let mut src = vec![false; n];
        for &c in &mb {
            src[c] = true;
        }
        let md = dist_from(&src, w, h);
        let matched = ma.iter().filter(|&&c| md[c] <= 2).count();
        let mdv = sorted(ma.iter().map(|&c| md[c] as f32).collect());
        eprintln!(
            "   {a_t:>8.0} {:>7.0} {ca:>9} {cb:>9} {:>6.1}% {:>6.1}% {:>6.1}% {:>6.1}%   {:>6.0} / {:<6.0} \
             {:>6.0} / {:<6.0}   {:>4} / {:<4} · {matched:>4} (moved p50 {:.0}, p90 {:.0})",
            ac,
            share(0),
            share(1),
            share(2),
            share(5),
            p(&s1, 0.50),
            p(&s1, 0.90),
            p(&s2, 0.50),
            p(&s2, 0.90),
            ma.len(),
            mb.len(),
            p(&mdv, 0.50),
            p(&mdv, 0.90)
        );
    }
}

#[test]
#[ignore]
fn f120_area() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 120 . A0 and A by DRAINED AREA (km²)  ==========");
    eprintln!(
        "   cell = {cell_km2:.6} km² at 8192² · 40 cells = {:.3} km² (below: sub-grid)",
        40.0 * cell_km2
    );

    // ── A0 ──
    let on = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (w, h) = (on.width, on.height);
    let n = w * h;
    let (b_on, f_on) = breached(&on, &ss);
    {
        let off = build_field_seed(Knobs::passes(2), PSEED);
        let (b_off, f_off) = breached(&off, &ss);
        drop(off);
        a0_pair("OFF (A) vs ON (B)", (&f_off, &b_off), (&f_on, &b_on), cell_km2);
    }
    {
        let pre = build_field_seed(Knobs::no_incision(), PSEED);
        let (b_pre, f_pre) = breached(&pre, &ss);
        drop(pre);
        a0_pair(
            "PRE-INCISION (A) vs ON (B) — the construction's skeleton",
            (&f_pre, &b_pre),
            (&f_on, &b_on),
            cell_km2,
        );
    }
    eprintln!("   ⏱ A0 {:.1} s", t0.elapsed().as_secs_f64());
    drop(f_on);

    // ── A · the exported network, full tree, ON ──
    let dr: C1DrainageResult = {
        let cfg = dcfg(true);
        let pre = c1_drainage_windowed(&on, None, &cfg, &ss, DOMAIN_KM);
        let cl = c1_climate_placed(&b_on, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc =
            DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        assemble_hd_drainage(&b_on, &dc, Some(pre), &cfg, &ss, DOMAIN_KM, GEO_RATIO, None, false)
            .drainage
    };
    let bf = &b_on;
    let acc = &dr.flow.accumulation;
    let band_of = |k: usize| -> Option<usize> {
        let a = acc.data[k] * cell_km2;
        BANDS.iter().position(|&(lo, hi)| a >= lo && a < hi)
    };
    let land: Vec<bool> = bf.data.iter().map(|&v| v > SEA).collect();
    // per-cell band of the RIVER network (Watercourse points), 255 = not a river point
    let mut rb = vec![255u8; n];
    let mut ord_by_band: Vec<Vec<f32>> = vec![Vec::new(); BANDS.len()];
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] != SegmentKind::Watercourse {
            continue;
        }
        for &(x, y) in &s.points {
            let k = y as usize * w + x as usize;
            if let Some(b) = band_of(k) {
                rb[k] = b as u8;
                ord_by_band[b].push(s.strahler_order as f32);
            }
        }
    }
    // D8 attribution of hillslope cells to the band of the first river point downstream
    let mut att = vec![255u8; n];
    for c in 0..n {
        if !land[c] {
            att[c] = 254;
        } else if rb[c] != 255 {
            att[c] = rb[c];
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
                val = if att[cur] == 253 { 254 } else { att[cur] };
                break;
            }
            att[cur] = 253;
            path.push(cur);
            let d = dr.flow.direction[cur];
            if d == DIR_NONE {
                val = 254;
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
    let slope = |k: usize| -> f32 {
        let (x, y) = (k % w, k / w);
        let at = |xx: usize, yy: usize| bf.data[yy * w + xx];
        let (xm, xp) = ((x + w - 1) % w, (x + 1) % w);
        let (ym, yp) = ((y + h - 1) % h, (y + 1) % h);
        let gx = (at(xp, y) - at(xm, y)) * n2m / (2.0 * CELL_M);
        let gy = (at(x, yp) - at(x, ym)) * n2m / (2.0 * CELL_M);
        (gx * gx + gy * gy).sqrt().atan().to_degrees()
    };
    eprintln!(
        "\n   ── A · per band of drained area, A1+B2 ON, exported network (full tree) ──\n   \
         {:>14} {:>7} {:>9} {:>9} {:>7} {:>9} {:>9} {:>9} {:>10} {:>16} {:>11}  {}",
        "A km²",
        "cells@lo",
        "points",
        "length km",
        "runs",
        "R8 c16",
        "R8 c32",
        "R2 c32",
        "R8 terr",
        "hill σ(°) μ",
        "hill cells",
        "Strahler p10/p50/p90"
    );
    for (bi, &(lo, hi)) in BANDS.iter().enumerate() {
        let mut th16 = Vec::new();
        let mut th32 = Vec::new();
        let (mut pts, mut runs, mut len) = (0usize, 0usize, 0f64);
        let mut band = vec![false; n];
        for (i, s) in dr.rivers.segments.iter().enumerate() {
            if dr.segment_kind[i] != SegmentKind::Watercourse {
                continue;
            }
            let mut run: Vec<(u32, u32)> = Vec::new();
            let flush = |run: &mut Vec<(u32, u32)>,
                         th16: &mut Vec<f32>,
                         th32: &mut Vec<f32>,
                         runs: &mut usize,
                         len: &mut f64| {
                if run.len() >= 2 {
                    *runs += 1;
                    th16.extend(chords(run, 16));
                    th32.extend(chords(run, 32));
                    for wv in run.windows(2) {
                        let diag = wv[0].0 != wv[1].0 && wv[0].1 != wv[1].1;
                        *len += if diag {
                            CELL_KM as f64 * std::f64::consts::SQRT_2
                        } else {
                            CELL_KM as f64
                        };
                    }
                }
                run.clear();
            };
            for &(x, y) in &s.points {
                let k = y as usize * w + x as usize;
                if band_of(k) == Some(bi) {
                    run.push((x, y));
                    pts += 1;
                    band[k] = true;
                } else {
                    flush(&mut run, &mut th16, &mut th32, &mut runs, &mut len);
                }
            }
            flush(&mut run, &mut th16, &mut th32, &mut runs, &mut len);
        }
        let bd = dist_from(&band, w, h);
        let bmask: Vec<bool> = (0..n).map(|c| land[c] && bd[c] <= 16).collect();
        drop(bd);
        let a = aniso(bf, &bmask, 16);
        let (_, r8_16) = axial(&th16);
        let (r2_32, r8_32) = axial(&th32);
        let (mut cnt, mut mean, mut m2) = (0f64, 0f64, 0f64);
        for c in 0..n {
            if land[c] && rb[c] == 255 && att[c] == bi as u8 {
                let v = slope(c) as f64;
                cnt += 1.0;
                let dlt = v - mean;
                mean += dlt / cnt;
                m2 += dlt * (v - mean);
            }
        }
        let sd = if cnt > 1.0 { (m2 / (cnt - 1.0)).sqrt() } else { f64::NAN };
        let so = sorted(ord_by_band[bi].clone());
        let cells_lo = lo / cell_km2;
        eprintln!(
            "   {:>6}–{:<7} {:>7.0}{} {pts:>9} {len:>9.0} {runs:>7} {r8_16:>9.4} {r8_32:>9.4} \
             {r2_32:>9.4} {:>10.4} {:>7.2} {:>7.1} {cnt:>11.0}  {:.0} / {:.0} / {:.0}",
            lo,
            if hi.is_finite() { format!("{hi}") } else { "∞".into() },
            cells_lo,
            if cells_lo < 40.0 { " SUB" } else { "    " },
            a.r8,
            sd,
            mean,
            p(&so, 0.10),
            p(&so, 0.50),
            p(&so, 0.90)
        );
    }

    // ── controls ──
    if let Some(i) = (0..dr.rivers.segments.len())
        .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway)
        .max_by(|&a, &b| dr.segment_discharge_m3s[a].total_cmp(&dr.segment_discharge_m3s[b]))
    {
        let s = &dr.rivers.segments[i];
        let av: Vec<f32> = s
            .points
            .iter()
            .map(|&(x, y)| acc.data[y as usize * w + x as usize] * cell_km2)
            .collect();
        let sv = sorted(av);
        eprintln!(
            "\n   CONTROL · F87-A (highest-discharge spillway, {:.2} m³/s, {} cells): D8 area along it \
             min / p50 / max **{:.1} / {:.1} / {:.1} km²**",
            dr.segment_discharge_m3s[i],
            s.points.len(),
            p(&sv, 0.0),
            p(&sv, 0.50),
            p(&sv, 1.0)
        );
    }
    let (tx, ty, ts) = COMB_TILE;
    let mut by: Vec<(usize, Vec<f32>)> = vec![(0, Vec::new()); BANDS.len()];
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] != SegmentKind::Watercourse {
            continue;
        }
        let mut run: Vec<(u32, u32)> = Vec::new();
        let mut cur: Option<usize> = None;
        for &(x, y) in &s.points {
            let inside = (x as usize) >= tx
                && (x as usize) < tx + ts
                && (y as usize) >= ty
                && (y as usize) < ty + ts;
            let b = if inside { band_of(y as usize * w + x as usize) } else { None };
            if b != cur || !inside {
                if let Some(cb) = cur {
                    by[cb].0 += run.len();
                    by[cb].1.extend(chords(&run, 32));
                }
                run.clear();
            }
            cur = b;
            if b.is_some() {
                run.push((x, y));
            }
        }
        if let Some(cb) = cur {
            by[cb].0 += run.len();
            by[cb].1.extend(chords(&run, 32));
        }
    }
    eprintln!(
        "   CONTROL · comb tile ({tx}, {ty}) {ts}² — river points by band: {}",
        by.iter()
            .enumerate()
            .filter(|(_, (c, _))| *c > 0)
            .map(|(b, (c, th))| format!(
                "{}–{} km² {c} pts (R8 c32 {:.3}, n {})",
                BANDS[b].0,
                BANDS[b].1,
                axial(th).1,
                th.len()
            ))
            .collect::<Vec<_>>()
            .join(" · ")
    );
    eprintln!(
        "\n==========  end Finding 120 area . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
