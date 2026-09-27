//! ADR 0001 Finding 124, part 4 — **break the plane: B2 (hierarchy) and B1 (profile), measured.**
//!
//! Finding 123's rule 18: the planar walls CREATE the comb's teeth (4.6 % have an A1+B2
//! counterpart; where one exists, the tooth is dead straight, 1.000, against 1.044). Two gated
//! remedies, on the definition (basin base) at the canonical framing:
//!
//! * **B2 — hierarchy**: tributaries built as valleys too, the same laws with a smaller W, down to
//!   `a_min_km2` = 1 km² and then A_c = 0.1 km² (42 cells at 8192²: resolved). The walls of the big
//!   trunks are then dissected by their tributaries' valleys.
//! * **B1 — profile**: concave foot, convex crest, and the tectonic + FBM detail the wall replaced
//!   (`WallProfile::f124`, every value PROXY).
//!
//! Columns: teeth (Finding 123-A's instrument on each world's OWN wall mask), the ex-teeth (C2/10's
//! teeth paired in each world, Finding 123's rule-18 pairing) — sinuosity and R8 at chord 8 —, R8
//! terrain, σ, the Δ class (Finding 108's gate), coast, the skeleton at 1 km² against PRE, cost.
//!
//! Run: cargo test -p ymir-core --release --test f124_walls -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, fill_field_m, land_u16, majority,
    over_dug_depression, pct, set_dump, sorted, take_bodies, to_mask,
};
use std::collections::{HashMap, VecDeque};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, SegmentKind, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::valley_construction::{
    F121_AGE_K, ValleyConstruction, WallProfile, carve, skeleton,
};
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowResult, RiverSegment, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const DELIVERED_P50_M: f32 = 424.2;
const CAP: u16 = 1000;

fn dcfg() -> C1DrainageConfig {
    let mut d = C1DrainageConfig::default();
    d.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    d.thresholds.full_tree = false;
    d
}

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

fn own(s: &RiverSegment) -> &[(u32, u32)] {
    if s.downstream.is_some() && s.points.len() >= 2 {
        &s.points[..s.points.len() - 1]
    } else {
        &s.points[..]
    }
}

fn sinuosity(pts: &[(u32, u32)]) -> f32 {
    if pts.len() < 4 {
        return f32::NAN;
    }
    let mut len = 0f32;
    for p in pts.windows(2) {
        let diag = p[0].0 != p[1].0 && p[0].1 != p[1].1;
        len += if diag { std::f32::consts::SQRT_2 } else { 1.0 };
    }
    let (a, b) = (pts[0], pts[pts.len() - 1]);
    let chord = ((a.0 as f32 - b.0 as f32).powi(2) + (a.1 as f32 - b.1 as f32).powi(2)).sqrt();
    if chord < 1.0 { f32::NAN } else { len / chord }
}

fn r8_chords(segs: &[&[(u32, u32)]], l: usize) -> f32 {
    let (mut c8, mut s8, mut n) = (0f64, 0f64, 0usize);
    for pts in segs {
        let mut i = 0usize;
        while i + l < pts.len() {
            let (dx, dy) =
                (pts[i + l].0 as f64 - pts[i].0 as f64, pts[i + l].1 as f64 - pts[i].1 as f64);
            if dx != 0.0 || dy != 0.0 {
                let t = dy.atan2(dx);
                c8 += (8.0 * t).cos();
                s8 += (8.0 * t).sin();
                n += 1;
            }
            i += l;
        }
    }
    if n == 0 { f32::NAN } else { ((c8 * c8 + s8 * s8).sqrt() / n as f64) as f32 }
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

fn breached_flow(f: &GridF32, ss: &SteinSteinParams) -> (GridF32, FlowResult) {
    let (w, h) = (f.width, f.height);
    let pre = c1_drainage_windowed(f, None, &dcfg(), ss, DOMAIN_KM);
    let bf = breach_monotone(f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let flow = c1_drainage_windowed(&bf, None, &dcfg(), ss, DOMAIN_KM).flow;
    (bf, flow)
}

/// Trunks A ≥ a_t km²: share of PRE trunk cells within 5 cells of the world's, p90, mouths matched.
fn skeleton_at(
    pre: (&FlowResult, &GridF32),
    world: (&FlowResult, &GridF32),
    a_t: f32,
) -> (f32, f32, usize, usize) {
    let (fa, ba) = pre;
    let (fb, bb) = world;
    let (w, h) = (ba.width, ba.height);
    let n = w * h;
    let ac = a_t / (CELL_KM * CELL_KM);
    let ta: Vec<bool> = (0..n).map(|k| ba.data[k] > SEA && fa.accumulation.data[k] >= ac).collect();
    let tb: Vec<bool> = (0..n).map(|k| bb.data[k] > SEA && fb.accumulation.data[k] >= ac).collect();
    let db = dist_from(&tb, w, h);
    let v: Vec<f32> = (0..n).filter(|&k| ta[k]).map(|k| db[k] as f32).collect();
    let within = 100.0 * v.iter().filter(|&&d| d <= 5.0).count() as f32 / v.len().max(1) as f32;
    let mouths = |f: &FlowResult, b: &GridF32| -> Vec<usize> {
        (0..n)
            .filter(|&k| {
                b.data[k] > SEA && f.accumulation.data[k] >= ac && f.direction[k] != DIR_NONE
            })
            .filter(|&k| {
                let d = f.direction[k] as usize;
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                b.data[ny * w + nx] <= SEA
            })
            .collect()
    };
    let (ma, mb) = (mouths(fa, ba), mouths(fb, bb));
    let mut src = vec![false; n];
    for &c in &mb {
        src[c] = true;
    }
    let md = dist_from(&src, w, h);
    let matched = ma.iter().filter(|&&c| md[c] <= 2).count();
    (within, pct(&sorted(v), 0.9), matched, ma.len())
}

#[test]
#[ignore]
fn f124_walls() {
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 124-4 . break the plane: B2 and B1  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let auth = spurs(&land_u16(&pre, &ss), w);
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc =
            DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let (b_pre, f_pre) = breached_flow(&pre, &ss);
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let worlds: [(&str, ValleyConstruction); 5] = [
        ("C2/10 (reference)", base),
        ("B2 tributaries ≥ 1 km²", ValleyConstruction { a_min_km2: 1.0, ..base }),
        ("B2 down to A_c 0.1 km²", ValleyConstruction { a_min_km2: 0.1, ..base }),
        ("B1 wall profile", ValleyConstruction { wall_profile: Some(WallProfile::f124()), ..base }),
        (
            "B2 (A_c) + B1",
            ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(WallProfile::f124()), ..base },
        ),
    ];
    let mut ref_teeth: Vec<(u32, u32, f32)> = Vec::new(); // own mouth x, y, own area — C2/10's teeth
    eprintln!(
        "\n   {:<24} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>7} {:>8} {:>12} {:>7}",
        "world",
        "teeth",
        "t len%",
        "ex paired",
        "ex sin50",
        "ex R8c8",
        "R8 terr",
        "σ p50",
        "Δ e/b",
        "coast",
        "skel≤5 1km²",
        "mouths 1km²",
        "build s"
    );
    for (label, vc) in worlds {
        let t = Instant::now();
        let g = build_field_seed(
            Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
            PSEED,
        );
        let secs = t.elapsed().as_secs_f64();
        // the world's OWN wall mask (carve of PRE with this world's parameters)
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let wall: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k]).collect();
        drop((sk, mk));
        // drainage, Δ class, coast
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate {
            precip_internal: &cl_e.precipitation,
            temperature: &cl_e.temperature,
        };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        let (fill_bre, _) = fill_field_m(&bre, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = common::f95_criteria(
            &g,
            &bre,
            &pre,
            DELIVERED_P50_M,
            &ss,
            &dcfg(),
            cell_km2,
            n2m,
            w,
            h,
        );
        set_dump(false);
        let (mut ce, mut cb) = (0usize, 0usize);
        for b in take_bodies() {
            let floor =
                *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            ce += over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) as usize;
            cb += over_dug_depression(fill_bre[floor] - fill_pre[floor], b.rim) as usize;
        }
        drop((fill_del, fill_bre));
        let coast = spurs(&land_u16(&bre, &ss), w) as i64 - auth as i64;
        let land: Vec<bool> = (0..n).map(|k| g.data[k] > SEA).collect();
        let r8t = aniso(&g, &land, 16).r8;
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
                    v[i] = m(&g, (y + dy - 1) * w + x + dx - 1);
                    i += 1;
                }
            }
            let mu = v.iter().sum::<f32>() / 9.0;
            sig.push((v.iter().map(|a| (a - mu) * (a - mu)).sum::<f32>() / 9.0).sqrt());
        }
        let sigma = pct(&sorted(sig), 0.5);
        // the skeleton at 1 km²
        let (bw, fw) = breached_flow(&g, &ss);
        let (within, _p90, matched, nm) = skeleton_at((&f_pre, &b_pre), (&fw, &bw), 1.0);
        drop((bw, fw));
        // the exported network: teeth, and the ex-teeth
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate {
            precip_internal: &cl_b.precipitation,
            temperature: &cl_b.temperature,
        };
        let dr = assemble_hd_drainage(
            &bre,
            &dc_b,
            Some(dd),
            &dcfg(),
            &ss,
            DOMAIN_KM,
            GEO_RATIO,
            None,
            false,
        )
        .drainage;
        let segs = &dr.rivers.segments;
        let is_wc =
            |i: usize| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2;
        let a_at =
            |(x, y): (u32, u32)| dr.flow.accumulation.data[y as usize * w + x as usize] * cell_km2;
        let seg_len = |s: &RiverSegment| -> f32 {
            s.points
                .windows(2)
                .map(|p| if p[0].0 != p[1].0 && p[0].1 != p[1].1 { 1.414 } else { 1.0 })
                .sum::<f32>()
        };
        let teeth: Vec<usize> = (0..segs.len())
            .filter(|&i| is_wc(i))
            .filter(|&i| {
                let p = &segs[i].points;
                p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32
                    >= 0.8 * p.len() as f32
            })
            .collect();
        let tot: f32 = (0..segs.len()).filter(|&i| is_wc(i)).map(|i| seg_len(&segs[i])).sum();
        let tlen: f32 = teeth.iter().map(|&i| seg_len(&segs[i])).sum();
        if ref_teeth.is_empty() {
            for &t in &teeth {
                let &(x, y) = own(&segs[t]).last().unwrap();
                ref_teeth.push((x, y, a_at((x, y))));
            }
        }
        // the ex-teeth: C2/10's teeth paired in this world (own mouth ≤ 3 cells, area ×/÷ 1.5)
        let mut by_mouth: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in (0..segs.len()).filter(|&i| is_wc(i)) {
            let &(x, y) = own(&segs[i]).last().unwrap();
            by_mouth.entry(y as usize * w + x as usize).or_default().push(i);
        }
        let mut paired: Vec<&[(u32, u32)]> = Vec::new();
        for &(mx, my, at) in &ref_teeth {
            let mut best: Option<(i64, usize)> = None;
            for dy in -3i64..=3 {
                for dx in -3i64..=3 {
                    let (x, y) = (mx as i64 + dx, my as i64 + dy);
                    if x < 0 || y < 0 || x >= w as i64 || y >= h as i64 {
                        continue;
                    }
                    if let Some(v) = by_mouth.get(&(y as usize * w + x as usize)) {
                        for &o in v {
                            let r = a_at(*own(&segs[o]).last().unwrap()) / at.max(1e-6);
                            let d = dx.abs().max(dy.abs());
                            if (1.0 / 1.5..=1.5).contains(&r) && best.is_none_or(|(bd, _)| d < bd) {
                                best = Some((d, o));
                            }
                        }
                    }
                }
            }
            if let Some((_, o)) = best {
                paired.push(own(&segs[o]));
            }
        }
        let sin = sorted(paired.iter().map(|p| sinuosity(p)).filter(|v| v.is_finite()).collect());
        eprintln!(
            "   {label:<24} {:>7} {:>6.1}% {:>8} {:>8.3} {:>8.4} {:>8.4} {:>8.2} {:>3}/{:<4} {:>+7} {:>7.1}% {:>5}/{:<6} {:>7.0}",
            teeth.len(),
            100.0 * tlen / tot.max(1e-6),
            paired.len(),
            if sin.is_empty() { f32::NAN } else { pct(&sin, 0.5) },
            r8_chords(&paired, 8),
            r8t,
            sigma,
            ce,
            cb,
            coast,
            within,
            matched,
            nm,
            secs
        );
    }
    eprintln!("\n==========  end 124-4 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
