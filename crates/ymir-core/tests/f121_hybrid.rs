//! ADR 0001 Finding 121 — **the hybrid at 8192²: the SHAPE constructed along the pre-incision
//! trunks, the TEXTURE left to a light incision.** No promotion.
//!
//! ## The three reasons Option A was rejected (C1 design doc §3.1, verbatim)
//!
//! 1. *"cannot produce coherent chronological diversity (young vs old mountain ranges on the same
//!    continent)"*;
//! 2. *"parameter calibration drifts indefinitely with no physical anchor"*;
//! 3. *"'looks invented' is a real failure mode for users who can identify imposed geometry by
//!    eye"*.
//!
//! This round answers (3) with the author's eye (block D) and the σ / R8 columns; (2) with the
//! labels — every law is ANCHORED or PROXY, and the one PROXY on a target (`k`) is calibrated on a
//! measured oracle, not tuned; (1) with block C3 (the age as a geometric parameter), which is
//! necessary for chronological diversity and NOT sufficient: a single `k` per continent says
//! nothing about young and old ranges side by side.
//!
//! ## Domain (rule 12)
//!
//! * **Skeleton** — the D8 network of the breached field the incision would receive (FBM + C-2
//!   craters), Finding 120's instrument, built INSIDE `upscale_from_c1` by the seam. Trunks
//!   `A ≥ 10 km²` in-domain (4 194 cells at 8192², declared).
//! * **The worlds are built through the production config** (`Knobs::valley` sets
//!   `FbmUpscaleConfig::valley_construction`, the field the viz sets), so a bench world and a viz
//!   world with the same parameters are the same world.
//! * **k is calibrated on the PRE build** (`Knobs::no_incision`), whose field also carries the C-2
//!   rim reconstruction and the bathymetry the in-pipeline construction does not see; the rims touch
//!   crater cells only and the bathymetry sub-sea cells only — declared, not measured.
//! * **The interfluve and wall masks** come from a carve of that same PRE build (the in-pipeline
//!   carve is not returned by the pipeline); declared.
//!
//! Run: cargo test -p ymir-core --release --test f121_hybrid -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, fill_field_m, land_u16, majority,
    over_dug_depression, pct, set_dump, sorted, take_bodies, to_mask,
};
use std::collections::{HashMap, HashSet, VecDeque};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::valley_construction::{ValleyConstruction, carve, skeleton};
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowResult, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const DELIVERED_P50_M: f32 = 424.2;
/// Finding 95's oracle, the target of `k` (PROXY). Finding 96 calibrated its χ field on it.
const ORACLE_P50_M: f32 = 488.1;
const A_TRUNK_KM2: f32 = 10.0;
const CAP: u16 = 1000;
const COMB_TILE: (usize, usize, usize) = (2048, 5120, 1024);

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

/// Finding 102's shoreline-development share (copied from `f109_toggle`).
fn dl_share(lake_map: &[u32], w: usize, h: usize, cell_km2: f32) -> f32 {
    let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
    for (k, &id) in lake_map.iter().enumerate().take(w * h) {
        if id != 0 {
            by.entry(id).or_default().push(k);
        }
    }
    let (mut n, mut hi) = (0usize, 0usize);
    for cells in by.values() {
        let a = cells.len() as f32 * cell_km2;
        if a < 1.0 {
            continue;
        }
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut per = 0usize;
        for &k in cells {
            let (x, y) = (k % w, k / w);
            for (dx, dy) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                let nx = (x as i32 + dx).rem_euclid(w as i32) as usize;
                let ny = (y as i32 + dy).rem_euclid(h as i32) as usize;
                if !inside.contains(&(ny * w + nx)) {
                    per += 1;
                }
            }
        }
        n += 1;
        let p = per as f32 * (400.0 / w as f32);
        if p / (2.0 * (std::f32::consts::PI * a).sqrt()) > 5.0 {
            hi += 1;
        }
    }
    100.0 * hi as f32 / n.max(1) as f32
}

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
    let mut d = vec![CAP; w * h];
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

/// Finding 120's A0 pair at one threshold: (exact %, ≤5 %, A→B p50, A→B p90, B→A p90, mouths A,
/// mouths matched ≤ 2 cells).
fn overlay(
    (fa, ba): (&FlowResult, &GridF32),
    (fb, bb): (&FlowResult, &GridF32),
    a_t: f32,
) -> (f32, f32, f32, f32, f32, usize, usize) {
    let (w, h) = (ba.width, ba.height);
    let n = w * h;
    let ac = a_t / (CELL_KM * CELL_KM);
    let ta: Vec<bool> = (0..n).map(|k| ba.data[k] > SEA && fa.accumulation.data[k] >= ac).collect();
    let tb: Vec<bool> = (0..n).map(|k| bb.data[k] > SEA && fb.accumulation.data[k] >= ac).collect();
    let db = dist_from(&tb, w, h);
    let da = dist_from(&ta, w, h);
    let v1: Vec<f32> = (0..n).filter(|&k| ta[k]).map(|k| db[k] as f32).collect();
    let v2: Vec<f32> = (0..n).filter(|&k| tb[k]).map(|k| da[k] as f32).collect();
    let sh =
        |d: f32| 100.0 * v1.iter().filter(|&&v| v <= d).count() as f32 / v1.len().max(1) as f32;
    let (s1, s2) = (sorted(v1.clone()), sorted(v2));
    let (ma, mb) = (mouths(fa, ba, ac), mouths(fb, bb, ac));
    let mut src = vec![false; n];
    for &c in &mb {
        src[c] = true;
    }
    let md = dist_from(&src, w, h);
    let matched = ma.iter().filter(|&&c| md[c] <= 2).count();
    (sh(0.0), sh(5.0), p(&s1, 0.5), p(&s1, 0.9), p(&s2, 0.9), ma.len(), matched)
}

fn breached_flow(f: &GridF32, ss: &SteinSteinParams) -> (GridF32, FlowResult) {
    let (w, h) = (f.width, f.height);
    let pre = c1_drainage_windowed(f, None, &dcfg(false), ss, DOMAIN_KM);
    let bf = breach_monotone(f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let flow = c1_drainage_windowed(&bf, None, &dcfg(false), ss, DOMAIN_KM).flow;
    (bf, flow)
}

struct Row {
    name: String,
    secs: f64,
    coast: i64,
    canyons_e: usize,
    canyons_b: usize,
    scanned: usize,
    lakes: usize,
    fam1: usize,
    exo: usize,
    lake_pct: f32,
    dl: f32,
    relief_p50: f32,
    r8_terr: f32,
    r8_net: f32,
    r8_net_trunk: f32,
    r8_comb_terr: f32,
    r8_comb_net: f32,
    sigma_p50: f32,
    sigma_p90: f32,
    inthalweg: f32,
    inthalweg_trunk: f32,
    trans_p50: f32,
    to_nothing: usize,
    unresolved: usize,
    inter: (f32, f32, f32),
    wall_p50: f32,
    skel: Option<(f32, f32, f32, f32, f32, usize, usize)>,
}

#[test]
#[ignore]
fn f121_hybrid() {
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 121 . the hybrid at 8192²  ==========");
    eprintln!(
        "   trunk threshold A ≥ {A_TRUNK_KM2} km² = {:.0} cells at 8192² (≥ 40: resolved)",
        A_TRUNK_KM2 / cell_km2
    );

    // ── PRE, the skeleton, and k ────────────────────────────────────────────
    let t = Instant::now();
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let pre_secs = t.elapsed().as_secs_f64();
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let land_pre: Vec<bool> = pre.data.iter().map(|&v| v > SEA).collect();
    let probe = ValleyConstruction::f121(0.0, None);
    let t = Instant::now();
    let sk = skeleton(&pre, &probe, &ss, DOMAIN_KM);
    let sk_secs = t.elapsed().as_secs_f64();
    // K1: the floor SURFACE `base + k·χ` over paired land has p50 = the oracle (Finding 96's
    // calibration, read literally from the round: "le relief p50 des fonds")
    let p50_at = |k: f32| -> f32 {
        let v: Vec<f32> = (0..n)
            .filter(|&c| land_pre[c] && !sk.chi_m[c].is_nan())
            .map(|c| sk.floor_m(c, k))
            .collect();
        pct(&sorted(v), 0.50)
    };
    let (mut lo, mut hi) = (0.0f32, 1.0f32);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if p50_at(mid) < ORACLE_P50_M {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let k_age = 0.5 * (lo + hi);
    let n_trunk = sk.trunk.iter().filter(|&&b| b).count();
    let ws: Vec<f32> =
        (0..n).filter(|&c| sk.trunk[c]).map(|c| probe.width_m(sk.area_km2[c])).collect();
    let ws = sorted(ws);
    eprintln!(
        "\n   ── k, the age (PROXY on Finding 95's oracle) ──\n   floor surface paired p50 = {:.1} m at \
         **k = {k_age:.5}** (Finding 96: 0.0646 all-land, 0.0719 clamped at A_c) · skeleton {:.1} s \
         · trunk cells {n_trunk} ({:.0} km) · {} polylines · law W p10/p50/p90 **{:.0} / {:.0} / {:.0} \
         m**",
        p50_at(k_age),
        sk_secs,
        n_trunk as f32 * CELL_KM,
        sk.polylines.len(),
        p(&ws, 0.10),
        p(&ws, 0.50),
        p(&ws, 0.90)
    );
    // the masks, from a carve of the PRE build at the retained k (declared)
    let vc1 = ValleyConstruction::f121(k_age, None);
    let t = Instant::now();
    let (pre_carved, masks) = carve(&pre, &sk, &vc1, &ss);
    let carve_secs = t.elapsed().as_secs_f64();
    let inter: Vec<bool> = (0..n).map(|c| land_pre[c] && !masks.carved[c]).collect();
    let wall: Vec<bool> = (0..n).map(|c| masks.carved[c] && !masks.floor[c]).collect();
    let carved_km2 = masks.carved.iter().filter(|&&b| b).count() as f32 * cell_km2;
    let floor_km2 = masks.floor.iter().filter(|&&b| b).count() as f32 * cell_km2;
    let land_km2 = land_pre.iter().filter(|&&b| b).count() as f32 * cell_km2;
    eprintln!(
        "   carve of PRE at k: {carve_secs:.1} s · carved **{carved_km2:.0} km²** ({:.1} % of land) · \
         floors {floor_km2:.0} km² · walls {:.0} km²",
        100.0 * carved_km2 / land_km2,
        carved_km2 - floor_km2
    );
    drop(pre_carved);
    let a_scale: Vec<f32> = sk.area_km2.clone();
    drop(sk);

    // ── references every world is scored against ────────────────────────────
    let auth = spurs(&land_u16(&pre, &ss), w);
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc =
            DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(false), &ss, &dc, DOMAIN_KM).0
    };
    let (b_pre, f_pre) = breached_flow(&pre, &ss);
    let slope_deg = |g: &GridF32, k: usize| -> f32 {
        let (x, y) = (k % w, k / w);
        let at = |xx: usize, yy: usize| g.data[yy * w + xx];
        let (xm, xp) = ((x + w - 1) % w, (x + 1) % w);
        let (ym, yp) = ((y + h - 1) % h, (y + 1) % h);
        let gx = (at(xp, y) - at(xm, y)) * n2m / (2.0 * CELL_M);
        let gy = (at(x, yp) - at(x, ym)) * n2m / (2.0 * CELL_M);
        (gx * gx + gy * gy).sqrt().atan().to_degrees()
    };

    // ── one world, every column ─────────────────────────────────────────────
    let read = |name: &str, g: &GridF32, secs: f64, control: bool| -> Row {
        let dd = c1_drainage_windowed(g, None, &dcfg(false), &ss, DOMAIN_KM);
        let bre = breach_monotone(g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate {
            precip_internal: &cl_e.precipitation,
            temperature: &cl_e.temperature,
        };
        let (fill_del, _) = fill_field_m(g, &dcfg(false), &ss, &dc_e, DOMAIN_KM);
        let (fill_bre, _) = fill_field_m(&bre, &dcfg(false), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let cr = common::f95_criteria(
            g,
            &bre,
            &pre,
            DELIVERED_P50_M,
            &ss,
            &dcfg(false),
            cell_km2,
            n2m,
            w,
            h,
        );
        set_dump(false);
        let (mut ce, mut cb) = (0usize, 0usize);
        for b in take_bodies() {
            let floor = *b
                .cells
                .iter()
                .min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c]))
                .expect("non-empty");
            ce += over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) as usize;
            cb += over_dug_depression(fill_bre[floor] - fill_pre[floor], b.rim) as usize;
        }
        drop((fill_del, fill_bre));
        let coast = spurs(&land_u16(&bre, &ss), w) as i64 - auth as i64;
        let land: Vec<bool> = (0..n).map(|k| g.data[k] > SEA).collect();
        let r8_terr = aniso(g, &land, 16).r8;
        let (tx, ty, ts) = COMB_TILE;
        let comb = {
            let mut crop = GridF32::new(ts, ts, 0.0);
            for y in 0..ts {
                for x in 0..ts {
                    crop.data[y * ts + x] = g.data[(ty + y) * w + tx + x];
                }
            }
            let cl: Vec<bool> = crop.data.iter().map(|&v| v > SEA).collect();
            aniso(&crop, &cl, 16).r8
        };
        // σ local: 3×3 heights in metres, land, stride 37 (Finding 96's population)
        let mut sig = Vec::new();
        for k in (0..n).step_by(37) {
            let (x, y) = (k % w, k / w);
            if !land[k] || x == 0 || y == 0 || x + 1 >= w || y + 1 >= h {
                continue;
            }
            let mut v = [0f32; 9];
            let mut i = 0;
            for dy in 0..3 {
                for dx in 0..3 {
                    v[i] = m(g, (y + dy - 1) * w + x + dx - 1);
                    i += 1;
                }
            }
            let mu = v.iter().sum::<f32>() / 9.0;
            sig.push((v.iter().map(|a| (a - mu) * (a - mu)).sum::<f32>() / 9.0).sqrt());
        }
        let sig = sorted(sig);
        // interfluves and walls, on the masks of the PRE carve
        let di: Vec<f32> =
            (0..n).filter(|&c| inter[c] && land[c]).map(|c| m(g, c) - m(&pre, c)).collect();
        let di = sorted(di);
        let wl: Vec<f32> = (0..n).filter(|&c| wall[c]).map(|c| slope_deg(g, c)).collect();
        let wl = sorted(wl);
        // the exported network
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate {
            precip_internal: &cl_b.precipitation,
            temperature: &cl_b.temperature,
        };
        let dr = assemble_hd_drainage(
            &bre,
            &dc_b,
            Some(dd),
            &dcfg(false),
            &ss,
            DOMAIN_KM,
            GEO_RATIO,
            None,
            false,
        )
        .drainage;
        let (mut th, mut th_tr, mut th_comb) = (Vec::new(), Vec::new(), Vec::new());
        let (mut pts, mut inth, mut pts_tr, mut inth_tr) = (0usize, 0usize, 0usize, 0usize);
        let mut offs = Vec::new();
        for (i, s) in dr.rivers.segments.iter().enumerate() {
            if dr.segment_kind[i] != SegmentKind::Watercourse {
                continue;
            }
            th.extend(chords(&s.points, 32));
            let big: Vec<(u32, u32)> = s
                .points
                .iter()
                .copied()
                .filter(|&(x, y)| {
                    dr.flow.accumulation.data[y as usize * w + x as usize] * cell_km2 >= A_TRUNK_KM2
                })
                .collect();
            th_tr.extend(chords(&big, 32));
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
            th_comb.extend(chords(&inside, 32));
            // Finding 15 / 113: the TRANSVERSE thalweg test on the breached field
            for j in 1..s.points.len().saturating_sub(1) {
                let (x, y) = (s.points[j].0 as i32, s.points[j].1 as i32);
                let (px, py) = (
                    s.points[j + 1].0 as i32 - s.points[j - 1].0 as i32,
                    s.points[j + 1].1 as i32 - s.points[j - 1].1 as i32,
                );
                if px == 0 && py == 0 {
                    continue;
                }
                let (qx, qy) = (-py.signum(), px.signum());
                let at = |dx: i32, dy: i32| -> f32 {
                    let nx = (x + dx).rem_euclid(w as i32) as usize;
                    let ny = (y + dy).rem_euclid(h as i32) as usize;
                    bre.data[ny * w + nx] * n2m
                };
                let (me, a, b) = (at(0, 0), at(qx, qy), at(-qx, -qy));
                let ok = me <= a && me <= b;
                pts += 1;
                inth += ok as usize;
                offs.push(me - a.min(b));
                let k = y as usize * w + x as usize;
                if dr.flow.accumulation.data[k] * cell_km2 >= A_TRUNK_KM2 {
                    pts_tr += 1;
                    inth_tr += ok as usize;
                }
            }
        }
        let offs = sorted(offs);
        let lakes: Vec<_> = dr.lakes.iter().filter(|l| l.area_km2 >= 1.0).collect();
        // E, once, on the control world: the highest-discharge spillway
        if control {
            if let Some(i) = (0..dr.rivers.segments.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway)
                .max_by(|&a, &b| {
                    dr.segment_discharge_m3s[a].total_cmp(&dr.segment_discharge_m3s[b])
                })
            {
                let s = &dr.rivers.segments[i];
                let d8: Vec<f32> = s
                    .points
                    .iter()
                    .map(|&(x, y)| {
                        dr.flow.accumulation.data[y as usize * w + x as usize] * cell_km2
                    })
                    .collect();
                let d8 = sorted(d8);
                eprintln!(
                    "\n   E · the highest-discharge spillway on {name}: **{:.2} m³/s**, {} cells · \
                     `segment_drainage_km2` **{:.0} km²** signified = **{:.1} km²** in-domain (÷ \
                     7.5²) · D8 area along its path p50 **{:.1} km²** · source lake {:?} ⇒ {}",
                    dr.segment_discharge_m3s[i],
                    s.points.len(),
                    dr.segment_drainage_km2[i],
                    dr.segment_drainage_km2[i] / (GEO_RATIO * GEO_RATIO),
                    p(&d8, 0.5),
                    dr.segment_source_lake[i],
                    "a below-sea basin's spillway, traced over its col OUTSIDE the D8 network \
                     (hd_assembly.rs: \"the raster reads ~0 along it\")"
                );
            }
        }
        Row {
            name: name.into(),
            secs,
            coast,
            canyons_e: ce,
            canyons_b: cb,
            scanned: cr.scanned,
            lakes: lakes.len(),
            fam1: lakes.iter().filter(|l| l.base.id < 1_000_001).count(),
            exo: lakes.iter().filter(|l| l.lake_type == LakeType::Exorheic).count(),
            lake_pct: cr.lake_pct,
            dl: dl_share(&dr.lake_map, w, h, cell_km2),
            relief_p50: cr.p50,
            r8_terr,
            r8_net: axial(&th).1,
            r8_net_trunk: axial(&th_tr).1,
            r8_comb_terr: comb,
            r8_comb_net: axial(&th_comb).1,
            sigma_p50: p(&sig, 0.50),
            sigma_p90: p(&sig, 0.90),
            inthalweg: 100.0 * inth as f32 / pts.max(1) as f32,
            inthalweg_trunk: 100.0 * inth_tr as f32 / pts_tr.max(1) as f32,
            trans_p50: p(&offs, 0.50),
            to_nothing: cr.to_nothing,
            unresolved: cr.unresolved,
            inter: (p(&di, 0.10), p(&di, 0.50), p(&di, 0.90)),
            wall_p50: p(&wl, 0.50),
            skel: None,
        }
    };

    // Printed AS EACH WORLD IS READ: the first run lost six worlds to a panic in the seventh.
    let print_row = |r: &Row| {
        eprintln!(
            "\n   ROW {:<22} | coast {:+} · Δ e/b {}/{} of {} · lakes≥1 {} (fam1 {}, exo {}) · D_L>5 {:.0} % \
             · lake {:.2} % · relief p50 {:.1} m · to_nothing {} · unresolved {} · build {:.0} s\n   \
             ROW {:<22} | R8 terr {:.4} · R8 net {:.4} (≥10 km² {:.4}) · comb tile terr {:.4} net {:.4} \
             · σ p50/p90 {:.2}/{:.2} m · in-thalweg {:.1} % (≥10 km² {:.1} %) · transverse p50 {:+.2} m \
             · interfluve Δ p10/p50/p90 {:.1}/{:.1}/{:.1} m · wall slope p50 {:.1}°",
            r.name,
            r.coast,
            r.canyons_e,
            r.canyons_b,
            r.scanned,
            r.lakes,
            r.fam1,
            r.exo,
            r.dl,
            r.lake_pct,
            r.relief_p50,
            r.to_nothing,
            r.unresolved,
            r.secs,
            r.name,
            r.r8_terr,
            r.r8_net,
            r.r8_net_trunk,
            r.r8_comb_terr,
            r.r8_comb_net,
            r.sigma_p50,
            r.sigma_p90,
            r.inthalweg,
            r.inthalweg_trunk,
            r.trans_p50,
            r.inter.0,
            r.inter.1,
            r.inter.2,
            r.wall_p50
        );
    };
    let try_read = |name: &str, g: &GridF32, secs: f64, control: bool| -> Option<Row> {
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            read(name, g, secs, control)
        })) {
            Ok(r) => {
                print_row(&r);
                Some(r)
            }
            Err(e) => {
                let msg = e
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_default();
                eprintln!(
                    "\n   ROW {name:<22} | ⛔ THE READ PANICKED — the world is reported, not skipped: {msg}"
                );
                None
            }
        }
    };
    // ADR Finding 122 -- `YMIR_F121_WORLDS=pre,on,c2a,c2basin` reads a subset (tags: pre, off, on,
    // c1, c2a = C2 /10, c2b = C2 /3, c3, c2basin = C2 /10 with Finding 122-B's `basin_base`). Unset,
    // the bench is Finding 121's, unchanged; `c2basin` runs only when asked for.
    let only: Option<Vec<String>> = std::env::var("YMIR_F121_WORLDS")
        .ok()
        .map(|v| v.split(',').map(|t| t.trim().to_string()).collect());
    let want = |tag: &str| only.as_ref().is_none_or(|v| v.iter().any(|t| t == tag));
    let asked = |tag: &str| only.as_ref().is_some_and(|v| v.iter().any(|t| t == tag));
    let mut rows: Vec<Row> = Vec::new();
    if want("pre") {
        rows.extend(try_read("PRE (no incision)", &pre, pre_secs, false));
    }
    for (tag, name, knobs) in [
        ("off", "OFF (delivered)", Knobs::passes(2)),
        ("on", "ON A1+B2 (control)", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
    ] {
        if !want(tag) {
            continue;
        }
        let t = Instant::now();
        let g = build_field_seed(knobs, PSEED);
        let secs = t.elapsed().as_secs_f64();
        let control = name.starts_with("ON");
        rows.extend(try_read(name, &g, secs, control));
    }

    // ── C1, and the skeleton control FIRST ──────────────────────────────────
    let build_c = |vc: ValleyConstruction, closure: bool| -> (GridF32, f64) {
        let t = Instant::now();
        let g = build_field_seed(
            Knobs {
                valley: Some(vc),
                slope_floor_abs: if closure { Some(S_EQ) } else { None },
                ..Knobs::passes(2)
            },
            PSEED,
        );
        (g, t.elapsed().as_secs_f64())
    };
    let pass = if !want("c1") {
        eprintln!(
            "\n   C1 not rebuilt (world selection) — its skeleton control PASSED at Finding 121 \
             (PRE→C1 p90 5 cells, mouths 155 of 184)"
        );
        true
    } else {
        let (c1, c1_secs) = build_c(vc1, false);
        let (b_c1, f_c1) = breached_flow(&c1, &ss);
        let sk10 = overlay((&f_pre, &b_pre), (&f_c1, &b_c1), A_TRUNK_KM2);
        let sk100 = overlay((&f_pre, &b_pre), (&f_c1, &b_c1), 100.0);
        drop((b_c1, f_c1));
        eprintln!(
            "\n   ── C1 · SKELETON CONTROL · PRE vs C1, trunks A ≥ A_t (chessboard cells) ──\n   \
         A ≥ 10 km²: exact **{:.1} %** · ≤ 5 **{:.1} %** · PRE→C1 p50/p90 **{:.0} / {:.0}** · C1→PRE p90 \
         **{:.0}** · mouths **{} of {}** matched\n   A ≥ 100 km²: exact **{:.1} %** · ≤ 5 **{:.1} %** · \
         PRE→C1 p50/p90 **{:.0} / {:.0}** · C1→PRE p90 **{:.0}** · mouths **{} of {}** matched\n   \
         (Finding 120, PRE vs ON at 10 km²: exact 32.0 %, ≤5 78.4 %, p90 10 / 34, mouths 97 of 184)",
            sk10.0,
            sk10.1,
            sk10.2,
            sk10.3,
            sk10.4,
            sk10.6,
            sk10.5,
            sk100.0,
            sk100.1,
            sk100.2,
            sk100.3,
            sk100.4,
            sk100.6,
            sk100.5
        );
        // declared BEFORE the run: C1 must reproduce PRE at least as well as ON does
        let pass = sk10.3 <= 10.0 && (sk10.6 as f32) >= 0.527 * sk10.5 as f32;
        eprintln!(
            "   ⇒ skeleton control (PRE→C1 p90 ≤ 10 cells AND mouths matched ≥ 52.7 %, ON's own values): \
         **{}**",
            if pass {
                "PASSES"
            } else {
                "**FAILS — STOP: the construction moved its own trunks**"
            }
        );
        if let Some(mut r) = try_read("C1 bare construction", &c1, c1_secs, false) {
            r.skel = Some(sk10);
            rows.push(r);
        }
        drop(c1);
        pass
    };

    if pass {
        // ── C2, two light settings ──
        for (tag, name, vc) in [
            ("c2a", "C2 hybrid k_time/10", ValleyConstruction::f121(k_age, Some(0.1))),
            ("c2b", "C2 hybrid k_time/3", ValleyConstruction::f121(k_age, Some(1.0 / 3.0))),
            ("c2basin", "C2/10 + basin base F122", ValleyConstruction::f122(k_age, Some(0.1))),
        ] {
            if (tag == "c2basin" && !asked(tag)) || (tag != "c2basin" && !want(tag)) {
                continue;
            }
            let (g, secs) = build_c(vc, true);
            let (bg, fg) = breached_flow(&g, &ss);
            let s = overlay((&f_pre, &b_pre), (&fg, &bg), A_TRUNK_KM2);
            drop((bg, fg));
            eprintln!(
                "   skeleton PRE vs {name}: A ≥ 10 km² ≤5 {:.1} % · p90 {:.0}/{:.0} · mouths {} of {}",
                s.1, s.3, s.4, s.6, s.5
            );
            if let Some(mut r) = try_read(name, &g, secs, false) {
                r.skel = Some(s);
                rows.push(r);
            }
        }
        // ── C3, three ages at k_time/10 (the lightest setting; k×1 is the row above) ──
        for (name, mult) in [("C3 age k×0.7 (/10)", 0.7f32), ("C3 age k×1.4 (/10)", 1.4)] {
            if !want("c3") {
                continue;
            }
            let vc = ValleyConstruction::f121(k_age * mult, Some(0.1));
            let (g, secs) = build_c(vc, true);
            rows.extend(try_read(name, &g, secs, false));
        }
    }

    // ── the tables ──────────────────────────────────────────────────────────
    eprintln!(
        "\n   ── the columns (8192², seed 1; Δ class = Finding 108 gate, eroded / breached) ──\n   \
         {:<24} {:>6} {:>7} {:>9} {:>6} {:>5} {:>4} {:>7} {:>6} {:>8} {:>7}",
        "world",
        "coast",
        "Δ e/b",
        "lakes≥1",
        "fam1",
        "exo",
        "D_L",
        "lake %",
        "p50 m",
        "to_noth",
        "build s"
    );
    for r in &rows {
        eprintln!(
            "   {:<24} {:>+6} {:>3}/{:<3} {:>9} {:>6} {:>5} {:>3.0}% {:>7.2} {:>6.0} {:>8} {:>7.0}",
            r.name,
            r.coast,
            r.canyons_e,
            r.canyons_b,
            r.lakes,
            r.fam1,
            r.exo,
            r.dl,
            r.lake_pct,
            r.relief_p50,
            r.to_nothing,
            r.secs
        );
    }
    eprintln!(
        "\n   {:<24} {:>8} {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>20} {:>8}",
        "world",
        "R8 terr",
        "R8 net",
        "R8 net≥10",
        "comb terr",
        "comb net",
        "σ p50 m",
        "σ p90 m",
        "inthalw %",
        "≥10 km² %",
        "interfluve Δ p10/50/90",
        "wall °"
    );
    for r in &rows {
        eprintln!(
            "   {:<24} {:>8.4} {:>8.4} {:>9.4} {:>9.4} {:>9.4} {:>9.2} {:>9.2} {:>9.1} {:>9.1} \
             {:>6.1}/{:>6.1}/{:>6.1} {:>8.1}",
            r.name,
            r.r8_terr,
            r.r8_net,
            r.r8_net_trunk,
            r.r8_comb_terr,
            r.r8_comb_net,
            r.sigma_p50,
            r.sigma_p90,
            r.inthalweg,
            r.inthalweg_trunk,
            r.inter.0,
            r.inter.1,
            r.inter.2,
            r.wall_p50
        );
    }
    eprintln!(
        "   (transverse offset p50 m: {}) · (unresolved: {}) · (bodies scanned: {})",
        rows.iter().map(|r| format!("{:+.2}", r.trans_p50)).collect::<Vec<_>>().join(" / "),
        rows.iter().map(|r| r.unresolved.to_string()).collect::<Vec<_>>().join(" / "),
        rows.iter().map(|r| r.scanned.to_string()).collect::<Vec<_>>().join(" / ")
    );
    for r in rows.iter().filter(|r| r.skel.is_some()) {
        let s = r.skel.unwrap();
        eprintln!(
            "   skeleton PRE vs {:<22} A ≥ 10 km²: exact {:.1} % · ≤5 {:.1} % · p50/p90 {:.0}/{:.0} · \
             back p90 {:.0} · mouths {} of {}",
            r.name, s.0, s.1, s.2, s.3, s.4, s.6, s.5
        );
    }
    let _ = a_scale;
    eprintln!("\n==========  end Finding 121 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 121 — **the 2048² smoke run: the primitive TURNS, and nothing else is read from it**
/// (the author's amendment: no verdict of texture, comb, entrenchment or coast at a resolution
/// other than the product's, Finding 89-C5). Through the production config, bare construction
/// against the pre-incision field: it must lower some land and raise nothing.
///
/// Run: cargo test -p ymir-core --release --test f121_hybrid -- --ignored smoke --nocapture
#[test]
#[ignore]
fn f121_smoke_2048() {
    let t = Instant::now();
    let pre = build_field_seed(Knobs { target: Some(2048), ..Knobs::no_incision() }, PSEED);
    let t_pre = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let c1 = build_field_seed(
        Knobs {
            target: Some(2048),
            valley: Some(ValleyConstruction::f121(0.0719, None)),
            ..Knobs::passes(2)
        },
        PSEED,
    );
    let t_c1 = t.elapsed().as_secs_f64();
    let raised = (0..pre.data.len()).filter(|&k| c1.data[k] > pre.data[k]).count();
    let lowered = (0..pre.data.len()).filter(|&k| c1.data[k] < pre.data[k]).count();
    eprintln!(
        "\n   smoke 2048² · PRE {t_pre:.1} s · C1 bare {t_c1:.1} s · lowered **{lowered}** cells · \
         raised **{raised}** (must be 0) · NO VERDICT at this resolution"
    );
    assert!(lowered > 0, "the construction must bite");
    assert_eq!(raised, 0, "the bare construction must never raise a cell");
}
