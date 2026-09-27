//! ADR 0001 Finding 125 — **attribute before remedying: why B1 makes teeth, why B2 stripes, where
//! B2 → A_c digs.** Finding 124-P4 measured four remedies and each did something other than its
//! purpose. This bench writes NO remedy and modifies no primitive: it LOCATES each effect.
//!
//! Worlds (the definition, canonical framing, k_time/10 + A1+B2): the témoin C2/10 col; B1a (the
//! profile alone, `detail_gain` 0); B1b (the noise alone: a planar wall, `foot_m` 0, `crest_m` 0);
//! B1 (both, `WallProfile::f124`); B2 at 1 km²; B2 down to A_c; B2 (A_c) + B1 (run second, to locate
//! Finding 38's micro-pits before the others).
//!
//! Instruments (rule 18 on each line — the témoin is read with the SAME populations):
//! * R8 terrain DECOMPOSED by population: every fully-land 16×16 window's `e^{8iθ}` is split by the
//!   share of its cells in each class (interfluve / carved cells whose nearest trunk is ≥ 10, 3–10,
//!   1–3 or < 1 km²); the class contributions `Re(V_c · conj(V̂)) / N` add up to the total R8 exactly.
//! * the coast's spurs at the ERODED and at the BREACHED stage (Finding 90: spurs are born at the
//!   breach), each new spur (> 2 km from every témoin spur of the same stage) located against the
//!   world's COASTAL walls;
//! * the over-dug bodies (Finding 108's class) with their floor cell and the number of primitives
//!   (skeleton polylines) whose cone reaches below the field there;
//! * Finding 38's uncovered below-sea components (parsed from the invariant's own message), located;
//! * the 1 km² skeleton's R8 by band of drained area × the pre-incision terrain's local PLANARITY.
//!
//! Run: cargo test -p ymir-core --release --test f125_walls -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, build_field_seed, fill_field_m, land_u16, majority,
    over_dug_depression, pct, set_dump, sorted, take_bodies, to_mask,
};
use std::collections::VecDeque;
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
    CarveMasks, F121_AGE_K, Skeleton, ValleyConstruction, WallProfile, carve, skeleton,
};
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{RiverSegment, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const DELIVERED_P50_M: f32 = 424.2;
const CAP: u16 = 1000;
const NCLS: usize = 5;
const CLS_NAME: [&str; NCLS] = ["interfluve", "trunk ≥ 10", "sub 3–10", "sub 1–3", "sub < 1"];

fn dcfg() -> C1DrainageConfig {
    let mut d = C1DrainageConfig::default();
    d.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    d.thresholds.full_tree = false;
    d
}

fn own(s: &RiverSegment) -> &[(u32, u32)] {
    if s.downstream.is_some() && s.points.len() >= 2 { &s.points[..s.points.len() - 1] } else { &s.points[..] }
}

fn r8_chords(segs: &[&[(u32, u32)]], l: usize) -> f32 {
    let (mut c8, mut s8, mut n) = (0f64, 0f64, 0usize);
    for pts in segs {
        let mut i = 0usize;
        while i + l < pts.len() {
            let (dx, dy) = (pts[i + l].0 as f64 - pts[i].0 as f64, pts[i + l].1 as f64 - pts[i].1 as f64);
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
                let nk = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                if d[nk] > d[k] + 1 {
                    d[nk] = d[k] + 1;
                    q.push_back(nk);
                }
            }
        }
    }
    d
}

/// Spur mid-points (cell coordinates) of a land mask, Finding 75's detector.
fn spur_mids(land: &[bool], w: usize) -> Vec<(f32, f32)> {
    let (m, ww) = majority(land, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM)
        .0
        .iter()
        .map(|s| s.mid)
        .collect()
}

fn band(a_km2: f32) -> u8 {
    if a_km2 >= 10.0 {
        1
    } else if a_km2 >= 3.0 {
        2
    } else if a_km2 >= 1.0 {
        3
    } else {
        4
    }
}

/// The population class of every cell: sea 255; uncarved land 0; a carved cell takes the band of the
/// NEAREST trunk cell of its world's skeleton (chessboard BFS, a proxy for the carve's nearest sample).
fn classes(pre: &GridF32, sk: &Skeleton, m: &CarveMasks) -> Vec<u8> {
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let mut lab = vec![0u8; n];
    let mut seen = vec![false; n];
    let mut q = VecDeque::new();
    for k in 0..n {
        if sk.trunk[k] {
            lab[k] = band(sk.area_km2[k]);
            seen[k] = true;
            q.push_back(k);
        }
    }
    while let Some(k) = q.pop_front() {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nk = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                if !seen[nk] && m.carved[nk] {
                    seen[nk] = true;
                    lab[nk] = lab[k];
                    q.push_back(nk);
                }
            }
        }
    }
    (0..n)
        .map(|k| {
            if pre.data[k] <= SEA {
                255
            } else if m.carved[k] {
                if seen[k] { lab[k] } else { 1 }
            } else {
                0
            }
        })
        .collect()
}

/// R8 terrain decomposed by class: (total R8, per class (effective windows, class R8, contribution)).
fn r8_pop(f: &GridF32, cls: &[u8], win: usize) -> (f32, Vec<(f32, f32, f32)>) {
    let (w, h) = (f.width, f.height);
    let mut vr = [0f64; NCLS];
    let mut vi = [0f64; NCLS];
    let mut nc = [0f64; NCLS];
    let mut nwin = 0f64;
    for by in (0..h.saturating_sub(win)).step_by(win) {
        'win: for bx in (0..w.saturating_sub(win)).step_by(win) {
            let (mut jxx, mut jyy, mut jxy) = (0f64, 0f64, 0f64);
            let mut cnt = [0usize; NCLS];
            for y in by..by + win {
                for x in bx..bx + win {
                    let c = cls[y * w + x];
                    if c == 255 {
                        continue 'win;
                    }
                    cnt[c as usize] += 1;
                    let (gx, gy) = f.gradient_at(x, y);
                    jxx += (gx * gx) as f64;
                    jyy += (gy * gy) as f64;
                    jxy += (gx * gy) as f64;
                }
            }
            if 0.5 * (jxx + jyy) <= 1e-20 {
                continue;
            }
            let theta = 0.5 * (2.0 * jxy).atan2(jxx - jyy);
            let t = if theta < 0.0 { theta + std::f64::consts::PI } else { theta };
            let (zr, zi) = ((8.0 * t).cos(), (8.0 * t).sin());
            nwin += 1.0;
            for c in 0..NCLS {
                let fr = cnt[c] as f64 / (win * win) as f64;
                vr[c] += fr * zr;
                vi[c] += fr * zi;
                nc[c] += fr;
            }
        }
    }
    let (tr, ti) = (vr.iter().sum::<f64>(), vi.iter().sum::<f64>());
    let norm = (tr * tr + ti * ti).sqrt().max(1e-30);
    let total = (norm / nwin.max(1.0)) as f32;
    let per = (0..NCLS)
        .map(|c| {
            let r8c = ((vr[c] * vr[c] + vi[c] * vi[c]).sqrt() / nc[c].max(1e-9)) as f32;
            let contrib = ((vr[c] * tr + vi[c] * ti) / norm / nwin.max(1.0)) as f32;
            (nc[c] as f32, r8c, contrib)
        })
        .collect();
    (total, per)
}

/// Primitives (skeleton polylines) whose planar cone reaches below the field at cell `c`.
fn overlaps(sk: &Skeleton, zfield: &[f32], c: usize, tan: f32) -> usize {
    let (w, h) = (sk.width, sk.height);
    let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
    sk.polylines
        .iter()
        .filter(|l| {
            l.iter().any(|&(sx, sy, zf, hw)| {
                let mut dx = (cx - sx).rem_euclid(w as f32);
                if dx > w as f32 / 2.0 {
                    dx -= w as f32;
                }
                let mut dy = (cy - sy).rem_euclid(h as f32);
                if dy > h as f32 / 2.0 {
                    dy -= h as f32;
                }
                let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
                zf + (d - hw).max(0.0) * tan < zfield[c]
            })
        })
        .count()
}

/// RMS residual (m) of a least-squares plane over the (2r+1)² window centred on (x, y): low = planar.
fn planarity(zm: &[f32], w: usize, h: usize, x: usize, y: usize, r: i32) -> f32 {
    let (mut sz, mut sxz, mut syz, mut sxx, mut syy, mut n) = (0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
    for dy in -r..=r {
        for dx in -r..=r {
            let k = (y as i32 + dy).rem_euclid(h as i32) as usize * w + (x as i32 + dx).rem_euclid(w as i32) as usize;
            let z = zm[k] as f64;
            sz += z;
            sxz += dx as f64 * z;
            syz += dy as f64 * z;
            sxx += (dx * dx) as f64;
            syy += (dy * dy) as f64;
            n += 1.0;
        }
    }
    let (a, b, c) = (sxz / sxx, syz / syy, sz / n);
    let mut ss = 0f64;
    for dy in -r..=r {
        for dx in -r..=r {
            let k = (y as i32 + dy).rem_euclid(h as i32) as usize * w + (x as i32 + dx).rem_euclid(w as i32) as usize;
            let e = zm[k] as f64 - (a * dx as f64 + b * dy as f64 + c);
            ss += e * e;
        }
    }
    (ss / n).sqrt() as f32
}

/// The detail term B1 puts back: `field − box blur(field)` (m), radius r cells, on the torus.
fn detail(zm: &[f32], w: usize, h: usize, r: i64) -> Vec<f32> {
    let pass = |src: &[f32], horizontal: bool| -> Vec<f32> {
        let mut out = vec![0f32; w * h];
        let (len, lines) = if horizontal { (w, h) } else { (h, w) };
        let at = |line: usize, i: i64| -> usize {
            let i = i.rem_euclid(len as i64) as usize;
            if horizontal { line * w + i } else { i * w + line }
        };
        let win = (2 * r + 1) as f64;
        for line in 0..lines {
            let mut acc = 0f64;
            for i in -r..=r {
                acc += src[at(line, i)] as f64;
            }
            for i in 0..len as i64 {
                out[at(line, i)] = (acc / win) as f32;
                acc += src[at(line, i + r + 1)] as f64 - src[at(line, i - r)] as f64;
            }
        }
        out
    };
    let t = pass(&pass(zm, true), false);
    zm.iter().zip(&t).map(|(a, b)| a - b).collect()
}

struct World {
    name: &'static str,
    vc: ValleyConstruction,
}

#[test]
#[ignore]
fn f125_walls() {
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 125 . attribute the walls' four effects  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let zpre: Vec<f32> = (0..n).map(|k| m(&pre, k)).collect();
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let det = detail(&zpre, w, h, (500.0 / (CELL_KM * 1000.0)).round() as i64);
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let prof = WallProfile::f124();
    let worlds = [
        World { name: "C2/10 col (témoin)", vc: base },
        World { name: "B2 (A_c) + B1", vc: ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), ..base } },
        World {
            name: "B1a profile only",
            vc: ValleyConstruction { wall_profile: Some(WallProfile { detail_gain: 0.0, ..prof }), ..base },
        },
        World {
            name: "B1b noise only",
            vc: ValleyConstruction { wall_profile: Some(WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof }), ..base },
        },
        World { name: "B1 (both)", vc: ValleyConstruction { wall_profile: Some(prof), ..base } },
        World { name: "B2 ≥ 1 km²", vc: ValleyConstruction { a_min_km2: 1.0, ..base } },
        World { name: "B2 → A_c", vc: ValleyConstruction { a_min_km2: 0.1, ..base } },
    ];
    let mut g_t: Option<GridF32> = None;
    let mut spurs_t: Option<(Vec<(f32, f32)>, Vec<(f32, f32)>)> = None;
    let mut f38_cells: Vec<usize> = Vec::new();
    let mut summary: Vec<String> = Vec::new();
    for wd in worlds {
        let t = Instant::now();
        let g = build_field_seed(
            Knobs { valley: Some(wd.vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
            PSEED,
        );
        let secs = t.elapsed().as_secs_f64();
        let sk = skeleton(&pre, &wd.vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &wd.vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let tan = wd.vc.wall_deg.to_radians().tan();
        eprintln!("\n────────── {} · build {secs:.0} s ──────────", wd.name);
        // ── R8 terrain by population, the world and the témoin on the SAME classes ──
        let (r8w, popw) = r8_pop(&g, &cls, 16);
        let gt_ref = g_t.as_ref().unwrap_or(&g);
        let (r8t, popt) = r8_pop(gt_ref, &cls, 16);
        eprintln!("   R8 terrain: world {r8w:.4} · témoin (same classes) {r8t:.4} · Δ {:+.4}", r8w - r8t);
        for c in 0..NCLS {
            if popw[c].0 < 0.5 {
                continue;
            }
            eprintln!(
                "      {:<11} windows {:>8.1} · class R8 world {:.4} / témoin {:.4} · contribution world {:+.4} / témoin \
                 {:+.4} · Δ {:+.4}",
                CLS_NAME[c], popw[c].0, popw[c].1, popt[c].1, popw[c].2, popt[c].2, popw[c].2 - popt[c].2
            );
        }
        // ── the coast's spurs, eroded and breached, located against the COASTAL walls ──
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let se = spur_mids(&land_u16(&g, &ss), w);
        let sb = spur_mids(&land_u16(&bre, &ss), w);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let coastal_wall: Vec<bool> =
            (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&coastal_wall, w, h);
        let dcarved = dist_from(&mk.carved, w, h);
        let near = 2.0 / CELL_KM; // 2 km in cells
        let locate = |mine: &[(f32, f32)], refs: &[(f32, f32)]| -> String {
            let new: Vec<&(f32, f32)> = mine
                .iter()
                .filter(|a| refs.iter().all(|b| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt() > near))
                .collect();
            let cell = |p: &(f32, f32)| {
                (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize
            };
            let at_wall = new.iter().filter(|p| (dwall[cell(p)] as f32) <= near).count();
            let at_carve = new.iter().filter(|p| (dcarved[cell(p)] as f32) <= near).count();
            let first: Vec<String> =
                new.iter().take(6).map(|p| format!("({:.0},{:.0}) wall {} cells", p.0, p.1, dwall[cell(p)])).collect();
            format!(
                "{} spurs · NEW vs the témoin {} · of them within 2 km of a COASTAL wall {} · of any carved cell {} · e.g. {}",
                mine.len(),
                new.len(),
                at_wall,
                at_carve,
                first.join(", ")
            )
        };
        match &spurs_t {
            None => eprintln!("   spurs: eroded {} · breached {} (the témoin's reference)", se.len(), sb.len()),
            Some((te, tb)) => {
                eprintln!("   spurs ERODED: {}", locate(&se, te));
                eprintln!("   spurs BREACHED: {}", locate(&sb, tb));
            }
        }
        // ── Finding 108's class: the over-dug bodies, located, with their overlapping primitives ──
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = common::f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut ce = 0usize;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            let dfill = fill_del[floor] - fill_pre[floor];
            if over_dug_depression(dfill, b.rim) {
                ce += 1;
                eprintln!(
                    "   OVER-DUG body {} · {:.2} km² · floor ({},{}) · Δfill {dfill:.1} m · rim {:.1}° · class {} · primitives \
                     overlapping at the floor {} · carved {} / floor {}",
                    b.id,
                    b.km2,
                    floor % w,
                    floor / w,
                    b.rim,
                    CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea"),
                    overlaps(&sk, &zpre, floor, tan),
                    mk.carved[floor],
                    mk.floor[floor]
                );
            }
        }
        drop(fill_del);
        eprintln!("   over-dug bodies (eroded stage): {ce}");
        // ── the exported network: teeth and R8, or Finding 38's invariant, located ──
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage
        }));
        match res {
            Ok(dr) => {
                let segs = &dr.rivers.segments;
                let wall: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k]).collect();
                let wc: Vec<usize> = (0..segs.len())
                    .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2)
                    .collect();
                let teeth: Vec<usize> = wc
                    .iter()
                    .copied()
                    .filter(|&i| {
                        let p = &segs[i].points;
                        p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32 >= 0.8 * p.len() as f32
                    })
                    .collect();
                let mut by_cls = [0usize; NCLS];
                for &i in &teeth {
                    let &(x, y) = own(&segs[i]).last().unwrap();
                    let c = cls[y as usize * w + x as usize];
                    if (c as usize) < NCLS {
                        by_cls[c as usize] += 1;
                    }
                }
                let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
                eprintln!(
                    "   teeth {} (by class of the tooth's mouth: {}) · R8 network (chord 8) {:.4}",
                    teeth.len(),
                    (0..NCLS).filter(|&c| by_cls[c] > 0).map(|c| format!("{} {}", CLS_NAME[c], by_cls[c])).collect::<Vec<_>>().join(", "),
                    r8_chords(&pts, 8)
                );
                summary.push(format!("{:<20} teeth {:>5} · R8 terrain {r8w:.4} (Δ {:+.4}) · over-dug {ce} · build {secs:.0} s", wd.name, teeth.len(), r8w - r8t));
            }
            Err(e) => {
                let msg = e
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_default();
                eprintln!("   ⛔ the HD assembly PANICKED: {}", msg.lines().next().unwrap_or(""));
                if let Some(i) = msg.find("(floor cell, cells): [") {
                    let rest = &msg[i + 22..];
                    let body = &rest[..rest.find(']').unwrap_or(rest.len())];
                    for tup in body.split("),") {
                        let nums: Vec<usize> = tup
                            .split(|c: char| !c.is_ascii_digit())
                            .filter(|s| !s.is_empty())
                            .filter_map(|s| s.parse().ok())
                            .collect();
                        if nums.len() >= 2 {
                            f38_cells.push(nums[0]);
                        }
                    }
                }
                for &c in &f38_cells {
                    eprintln!(
                        "   F38 micro-pit floor ({},{}) · height {:.2} m (pre {:.1} m) · class {} · carved {} / floor {} · \
                         distance to the sea {} cells · the detail term there {:+.1} m · primitives overlapping {}",
                        c % w,
                        c / w,
                        m(&g, c),
                        zpre[c],
                        CLS_NAME.get(cls[c] as usize).copied().unwrap_or("sea"),
                        mk.carved[c],
                        mk.floor[c],
                        dsea[c],
                        det[c],
                        overlaps(&sk, &zpre, c, tan)
                    );
                }
                summary.push(format!("{:<20} ⛔ Finding 38 · R8 terrain {r8w:.4} · over-dug {ce} · build {secs:.0} s", wd.name));
            }
        }
        if !f38_cells.is_empty() && !wd.name.starts_with("B2 (A_c) + B1") {
            let hs: Vec<String> = f38_cells.iter().map(|&c| format!("({},{}) {:.2} m", c % w, c / w, m(&g, c))).collect();
            eprintln!("   heights at Finding 38's micro-pit floors in THIS world: {}", hs.join(" · "));
        }
        // ── the skeleton's R8 by band × planarity (B2 at 1 km² only) ──
        if (wd.vc.a_min_km2 - 1.0).abs() < 1e-6 {
            let mut per: [Vec<(f32, f64, f64)>; NCLS] = Default::default();
            for l in &sk.polylines {
                let mut i = 0usize;
                while i + 8 < l.len() {
                    let (a, b) = (l[i], l[i + 8]);
                    let (dx, dy) = ((b.0 - a.0) as f64, (b.1 - a.1) as f64);
                    if dx != 0.0 || dy != 0.0 {
                        let mid = l[i + 4];
                        let (x, y) = ((mid.0.floor() as i64).rem_euclid(w as i64) as usize, (mid.1.floor() as i64).rem_euclid(h as i64) as usize);
                        let c = band(sk.area_km2[y * w + x]) as usize;
                        let t = dy.atan2(dx);
                        per[c].push((planarity(&zpre, w, h, x, y, 10), (8.0 * t).cos(), (8.0 * t).sin()));
                    }
                    i += 8;
                }
            }
            eprintln!("   SKELETON (1 km²) R8 by band × planarity quartile of the pre-incision field (RMS plane residual over 21×21 cells):");
            for c in 1..NCLS {
                let v = &mut per[c];
                if v.len() < 40 {
                    continue;
                }
                v.sort_by(|a, b| a.0.total_cmp(&b.0));
                let q = v.len() / 4;
                let row: Vec<String> = (0..4)
                    .map(|k| {
                        let s = &v[k * q..if k == 3 { v.len() } else { (k + 1) * q }];
                        let (cr, ci) = (s.iter().map(|e| e.1).sum::<f64>(), s.iter().map(|e| e.2).sum::<f64>());
                        format!(
                            "Q{} (residual ≤ {:.1} m) R8 {:.3}",
                            k + 1,
                            s.last().map(|e| e.0).unwrap_or(0.0),
                            ((cr * cr + ci * ci).sqrt() / s.len() as f64) as f32
                        )
                    })
                    .collect();
                eprintln!("      {:<11} chords {:>6} · {}", CLS_NAME[c], v.len(), row.join(" · "));
            }
        }
        if g_t.is_none() {
            spurs_t = Some((se, sb));
            g_t = Some(g);
        }
    }
    eprintln!("\n   ── summary ──");
    for s in summary {
        eprintln!("   {s}");
    }
    let _ = pct(&sorted(vec![0.0]), 0.5);
    eprintln!("\n==========  end Finding 125 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 125, the prerequisite — the six states' eroded-cache digests under the CURRENT code
/// (no build): if they equal the viz's, the bench and the viz agree on the key and only the stored
/// reference is stale.
///
/// Run: cargo test -p ymir-core --release --test f125_walls -- --ignored f125_digests --nocapture
#[test]
#[ignore]
fn f125_digests() {
    let k = F121_AGE_K;
    let c2 = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    for (label, kn) in [
        ("livré", Knobs::passes(2)),
        ("A1+B2", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
        ("C1 nue", Knobs { valley: Some(ValleyConstruction::new(k, None)), ..Knobs::passes(2) }),
        ("C2 /10 col (défaut)", c2(ValleyConstruction::new(k, Some(0.1)))),
        ("C2 /3", c2(ValleyConstruction::new(k, Some(1.0 / 3.0)))),
        ("C2 /10 niveau mer", c2(ValleyConstruction::f121(k, Some(0.1)))),
    ] {
        eprintln!("   {label:<22} digest {}", common::bench_eroded_digest(kn, PSEED));
    }
    eprintln!("   Debug of the definition: {:?}", ValleyConstruction::new(k, Some(0.1)));
}
