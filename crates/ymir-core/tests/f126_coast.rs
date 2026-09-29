//! ADR 0001 Finding 126 — **the −1.00 m, the coast clause, the bent skeleton.**
//!
//! Part A runs first (the round's order: "une constante inconnue d'abord"): the two Finding 38 cells
//! of B2 (A_c) + B1 walked through the pipeline's OWN stages, through the benches' stage switches
//! (`upscale_from_c1`: FBM + C-2 craters → the construction → the light pass (k_time/10 + A1+B2) →
//! the droplet erosion → the C-2 rims → the bathymetry; the breach is the bench's, as in Finding
//! 125). Then the construction's terms at each cell, decomposed by switching the wall profile's
//! halves ON THE CONSTRUCTION'S OWN INPUT — and that input is checked to BE the pipeline's (the
//! direct carve must reproduce stage S2 bit for bit, or the decomposition is not read).
//!
//! Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f126_a --nocapture

mod common;

use common::{
    CELL_KM, Knobs, PSEED, SEA, build_field_seed, f95_criteria, fill_field_m, land_u16, majority,
    over_dug_depression, set_dump, take_bodies, to_mask,
};
use rayon::prelude::*;
use std::collections::VecDeque;
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::{WATER_CLASS_INLAND, WATER_CLASS_OCEAN, water_class};
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
const S_EQ: f32 = 0.024;
const GEO_RATIO: f32 = 7.5;
const DELIVERED_P50_M: f32 = 424.2;
const CAP: u16 = 1000;
/// Finding 125's two uncovered below-sea floor cells, (x, y).
const CELLS: [(usize, usize); 2] = [(4257, 2006), (4613, 2860)];
const NCLS: usize = 5;
const CLS_NAME: [&str; NCLS] = ["interfluve", "trunk ≥ 10", "sub 3–10", "sub 1–3", "sub < 1"];
/// Finding 126-C's positions on a wall, laid on the PLANAR geometry of the world's own skeleton
/// (the same bins for every world, so B1b's planar wall is read on B1's foot and crest).
const POS_NAME: [&str; 4] = ["floor", "foot", "straight", "crest"];

fn dcfg() -> C1DrainageConfig {
    let mut d = C1DrainageConfig::default();
    d.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    d.thresholds.full_tree = false;
    d
}

fn class_name(c: u8) -> &'static str {
    match c {
        WATER_CLASS_OCEAN => "OCEAN",
        WATER_CLASS_INLAND => "INLAND (enclosed)",
        _ => "land",
    }
}

/// Euclidean distance (cells) from `(x, y)` to the nearest cell with `wc == WATER_CLASS_OCEAN`,
/// by growing square rings until the ring lies beyond the best distance found. No torus (the
/// class itself is edge-seeded).
fn dist_to_ocean(wc: &[u8], w: usize, h: usize, x: usize, y: usize) -> f32 {
    let mut best = f32::INFINITY;
    let mut r = 0i64;
    while (r as f32) <= best && r < w.max(h) as i64 {
        for dy in -r..=r {
            for dx in -r..=r {
                if dx.abs() != r && dy.abs() != r {
                    continue;
                }
                let (nx, ny) = (x as i64 + dx, y as i64 + dy);
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                if wc[ny as usize * w + nx as usize] == WATER_CLASS_OCEAN {
                    best = best.min(((dx * dx + dy * dy) as f32).sqrt());
                }
            }
        }
        r += 1;
    }
    best
}

/// The land cells of `s1_land` that `g` puts at or below the sea, split by connection
/// (edge-connected OCEAN vs enclosed INLAND), with the enclosed components counted (8-connected,
/// over the NEW below-sea cells only) and the cells sitting at exactly −1.00 m.
fn below_census(
    s1_land: &[bool],
    g: &GridF32,
    masks: Option<&CarveMasks>,
    ss: &SteinSteinParams,
) -> String {
    let (w, h) = (g.width, g.height);
    let n = w * h;
    let wc = water_class(g, SEA);
    let new: Vec<bool> = (0..n).map(|k| s1_land[k] && g.data[k] <= SEA).collect();
    let total = new.iter().filter(|&&b| b).count();
    let ocean = (0..n).filter(|&k| new[k] && wc[k] == WATER_CLASS_OCEAN).count();
    let inland = (0..n).filter(|&k| new[k] && wc[k] == WATER_CLASS_INLAND).count();
    let at_m1 = (0..n)
        .filter(|&k| new[k] && (c1_altitude_norm_to_metres(g.data[k], ss) + 1.0).abs() < 0.005)
        .count();
    // enclosed components of NEW below-sea cells
    let mut seen = vec![false; n];
    let (mut comps, mut micro) = (0usize, 0usize);
    let mut q = VecDeque::new();
    for s in 0..n {
        if !new[s] || wc[s] != WATER_CLASS_INLAND || seen[s] {
            continue;
        }
        comps += 1;
        let mut size = 0usize;
        seen[s] = true;
        q.push_back(s);
        while let Some(k) = q.pop_front() {
            size += 1;
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if new[nk] && wc[nk] == WATER_CLASS_INLAND && !seen[nk] {
                        seen[nk] = true;
                        q.push_back(nk);
                    }
                }
            }
        }
        if size <= 4 {
            micro += 1;
        }
    }
    let split = masks.map_or(String::new(), |mk| {
        let walls = (0..n).filter(|&k| new[k] && mk.carved[k] && !mk.floor[k]).count();
        let floors = (0..n).filter(|&k| new[k] && mk.floor[k]).count();
        format!(" · on walls {walls} / on floors {floors} / uncarved {}", total - walls - floors)
    });
    format!(
        "land → ≤ sea: **{total}** cells ({:.2} km², domain){split} · edge-connected OCEAN {ocean} · \
         ENCLOSED {inland} in **{comps}** components ({micro} of ≤ 4 cells) · at exactly −1.00 m {at_m1}",
        total as f32 * CELL_KM * CELL_KM
    )
}

/// The skeleton sample nearest to cell `c` (torus distance, m): (line, index, distance).
fn nearest(sk: &Skeleton, c: usize) -> (usize, usize, f32) {
    let (w, h) = (sk.width, sk.height);
    let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
    let mut best = (usize::MAX, usize::MAX, f32::INFINITY);
    for (li, l) in sk.polylines.iter().enumerate() {
        for (i, &(sx, sy, _, _)) in l.iter().enumerate() {
            let mut dx = (cx - sx).rem_euclid(w as f32);
            if dx > w as f32 / 2.0 {
                dx -= w as f32;
            }
            let mut dy = (cy - sy).rem_euclid(h as f32);
            if dy > h as f32 / 2.0 {
                dy -= h as f32;
            }
            let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
            if d < best.2 {
                best = (li, i, d);
            }
        }
    }
    best
}

/// `field − box blur(field)` (m) at radius `r` cells on the torus — the construction's own detail
/// term (`valley_construction::box_blur_torus`, reproduced; checked through the recomposition).
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

#[test]
#[ignore]
fn f126_a_minus_one() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    eprintln!("\n==========  Finding 126-A . the −1.00 m, walked stage by stage  ==========");
    let prof = WallProfile::f124();
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let vc = ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), ..base };
    let world = Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let off = |k: Knobs| Knobs { erosion_off: true, bathymetry_off: true, ..k };

    // ── S1: the construction's input (FBM + C-2 craters; the C-2 rims also run, at the end) ──
    let t = Instant::now();
    let s1 = build_field_seed(Knobs { no_incision: true, ..off(Knobs::passes(2)) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let idx: Vec<usize> = CELLS.iter().map(|&(x, y)| y * w + x).collect();
    let s1_land: Vec<bool> = s1.data.iter().map(|&v| v > SEA).collect();
    let wc1 = water_class(&s1, SEA);
    eprintln!("   S1 built in {:.0} s", t.elapsed().as_secs_f64());

    // ── the pipeline's own stages for the world B2 (A_c) + B1 ──
    let mut rows: Vec<(String, Vec<(f32, u8)>, String)> = Vec::new();
    rows.push((
        "S1 input (FBM + craters)".into(),
        idx.iter().map(|&k| (m(&s1, k), wc1[k])).collect(),
        "—".into(),
    ));
    let stages: [(&str, Knobs); 4] = [
        ("S2 + construction", Knobs { valley: Some(vc), no_incision: true, ..off(Knobs::passes(2)) }),
        ("S3 + light pass (k/10, A1+B2)", off(world)),
        ("S4 + droplet erosion", Knobs { bathymetry_off: true, ..world }),
        ("S5 + bathymetry (the world)", world),
    ];
    let mut s2: Option<GridF32> = None;
    let mut s5: Option<GridF32> = None;
    for (label, k) in stages {
        let t = Instant::now();
        let g = build_field_seed(k, PSEED);
        let wc = water_class(&g, SEA);
        let census = below_census(&s1_land, &g, None, &ss);
        eprintln!("   {label} built in {:.0} s", t.elapsed().as_secs_f64());
        rows.push((label.into(), idx.iter().map(|&k| (m(&g, k), wc[k])).collect(), census));
        if label.starts_with("S2") {
            s2 = Some(g);
        } else if label.starts_with("S5") {
            s5 = Some(g);
        }
    }
    let s5 = s5.expect("S5");
    // ── S6: the bench's breach (Finding 125's), and whether a lake claims the cell ──
    let dd = c1_drainage_windowed(&s5, None, &dcfg(), &ss, DOMAIN_KM);
    let bre = breach_monotone(&s5, &dd.flow.filled, &dd.lake_map, SEA, w, h);
    let wc6 = water_class(&bre, SEA);
    rows.push((
        "S6 + breach (the HD assembly's input)".into(),
        idx.iter().map(|&k| (m(&bre, k), wc6[k])).collect(),
        below_census(&s1_land, &bre, None, &ss),
    ));
    eprintln!("\n   | stage | {} | {} | the world's census (cells land in S1) |",
        format!("({}, {})", CELLS[0].0, CELLS[0].1), format!("({}, {})", CELLS[1].0, CELLS[1].1));
    eprintln!("   |---|---|---|---|");
    for (label, v, census) in &rows {
        eprintln!(
            "   | {label} | {:.2} m · {} | {:.2} m · {} | {census} |",
            v[0].0,
            class_name(v[0].1),
            v[1].0,
            class_name(v[1].1)
        );
    }
    let wc5 = water_class(&s5, SEA);
    for (i, &k) in idx.iter().enumerate() {
        let (d1, d5) = (
            dist_to_ocean(&wc1, w, h, CELLS[i].0, CELLS[i].1),
            dist_to_ocean(&wc5, w, h, CELLS[i].0, CELLS[i].1),
        );
        eprintln!(
            "   cell ({}, {}): lake_map id at S5 {} · filled at S5 {:.2} m · distance to the OPEN OCEAN: S1 {d1:.0} cells \
             ({:.2} km domain) · S5 {d5:.0} cells ({:.2} km domain)",
            CELLS[i].0,
            CELLS[i].1,
            dd.lake_map[k],
            m(&dd.flow.filled, k),
            d1 * CELL_KM,
            d5 * CELL_KM,
        );
    }
    drop(wc5);
    drop(bre);
    drop(dd);
    drop(s5);

    // ── the direct carve on S1 must BE the pipeline's construction (S2), or nothing below is read ──
    let sk_ac = skeleton(&s1, &vc, &ss, DOMAIN_KM);
    let (c_full, mk_full) = carve(&s1, &sk_ac, &vc, &ss);
    let s2 = s2.expect("S2");
    let differ = (0..n).filter(|&k| c_full.data[k].to_bits() != s2.data[k].to_bits()).count();
    eprintln!(
        "\n   IDENTITY: the direct carve(skeleton(S1)) against the pipeline's S2 — {differ} cells differ of {n}{}",
        if differ == 0 { " ⇒ the decomposition below is read on the pipeline's own construction" } else { " ⇒ ⚠ NOT the pipeline's construction" }
    );
    for &k in &idx {
        eprintln!("      at ({}, {}): direct {:.2} m · S2 {:.2} m", k % w, k / w, m(&c_full, k), m(&s2, k));
    }
    drop(s2);

    // ── the construction's terms at each cell (the nearest sample, as `carve` lays it) ──
    let zfield: Vec<f32> = (0..n).map(|k| m(&s1, k)).collect();
    let det = detail(&zfield, w, h, (prof.detail_radius_m / sk_ac.cell_m).round().max(1.0) as i64);
    let tan = vc.wall_deg.to_radians().tan();
    let variants: [(&str, ValleyConstruction); 4] = [
        ("planar (B2 → A_c)", ValleyConstruction { wall_profile: None, ..vc }),
        ("profile only (B1a on A_c)", ValleyConstruction { wall_profile: Some(WallProfile { detail_gain: 0.0, ..prof }), ..vc }),
        ("noise only (B1b on A_c)", ValleyConstruction { wall_profile: Some(WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof }), ..vc }),
        ("both (B2 (A_c) + B1)", vc),
    ];
    let mut vals: Vec<Vec<f32>> = Vec::new();
    for (label, v) in variants {
        let (g, mk) = carve(&s1, &sk_ac, &v, &ss);
        vals.push(idx.iter().map(|&k| m(&g, k)).collect());
        eprintln!("   construction census · {label}: {}", below_census(&s1_land, &g, Some(&mk), &ss));
    }
    for (i, &k) in idx.iter().enumerate() {
        let (li, si, d) = nearest(&sk_ac, k);
        let (_, _, zf, hw) = sk_ac.polylines[li][si];
        let u = (d - hw).max(0.0);
        let v0 = zf + tan * u;
        let foot = prof.foot_m;
        let rise = if u < foot { tan * u * u / (2.0 * foot) } else { tan * (u - 0.5 * foot) };
        let v1 = zf + rise;
        let kc = prof.crest_m.max(1e-3);
        let crest = if v1 < zfield[k] { kc * (1.0 + (-(zfield[k] - v1) / kc).exp()).ln() } else { 0.0 };
        let v2 = v1 - crest;
        let v3 = v2 + prof.detail_gain * det[k];
        eprintln!(
            "\n   cell ({}, {}) · input {:.2} m · nearest sample: line {li} · floor zf {:.2} m · half-width {:.1} m · \
             distance {:.1} m · u = beyond the floor edge {:.1} m ({})",
            CELLS[i].0,
            CELLS[i].1,
            zfield[k],
            zf,
            hw,
            d,
            u,
            if d <= hw { "ON the floor" } else { "on the wall" }
        );
        eprintln!(
            "      planar zf + tan·u = {v0:.2} m · the profile's rise {rise:.2} m (planar {:.2}) → {v1:.2} m · crest −{crest:.2} m → \
             {v2:.2} m · detail {:+.2} m → **{v3:.2} m**",
            tan * u,
            det[k]
        );
        eprintln!(
            "      carve outputs: planar {:.2} m · profile only {:.2} m · noise only {:.2} m · both {:.2} m · \
             (recomposed planar {} the carve's, recomposed both {} the carve's)",
            vals[0][i],
            vals[1][i],
            vals[2][i],
            vals[3][i],
            if (vals[0][i] - v0.min(zfield[k])).abs() < 0.01 { "=" } else { "≠" },
            if (vals[3][i] - v3.min(zfield[k])).abs() < 0.01 { "=" } else { "≠" },
        );
        eprintln!("      carved {} / floor {}", mk_full.carved[k], mk_full.floor[k]);
    }

    // ── does the construction put land below the sea ELSEWHERE, in worlds where F38 does not fire? ──
    let sk_10 = skeleton(&s1, &base, &ss, DOMAIN_KM);
    for (label, v) in [
        ("témoin C2/10 col", base),
        ("B1a profile only", ValleyConstruction { wall_profile: Some(WallProfile { detail_gain: 0.0, ..prof }), ..base }),
        ("B1b noise only", ValleyConstruction { wall_profile: Some(WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof }), ..base }),
        ("B1 both", ValleyConstruction { wall_profile: Some(prof), ..base }),
    ] {
        let (g, mk) = carve(&s1, &sk_10, &v, &ss);
        eprintln!("   construction census · {label}: {}", below_census(&s1_land, &g, Some(&mk), &ss));
        let hs: Vec<String> = idx.iter().map(|&k| format!("{:.2} m", m(&g, k))).collect();
        eprintln!("      at the two cells: {}", hs.join(" · "));
    }
    eprintln!("\n==========  end Finding 126-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

// ════════════════════════════════ Parts B + C ════════════════════════════════
//
// B: the coast's base rate — spurs per km of coastline near a coastal wall and elsewhere, in every
// world — then the ONE clause (`wall_sea_floor_m`, gated) on B1a, B1 and B2 (A_c) + B1. C: where
// the teeth sit on their wall (floor / foot / straight / crest), with the wall's own area as the
// base rate. Finding 125's instruments throughout (the masks and classes of the world's skeleton on
// PRE, the same teeth, the same spur detector, the same over-dug class), so every row is comparable
// to Finding 125's.
//
// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f126_bc --nocapture

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

/// Finding 125's population classes (a carved cell takes its nearest trunk cell's band).
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

/// The skeleton's samples in buckets of `b` cells (torus-wrapped), for nearest-sample queries.
struct Buckets {
    b: usize,
    nbx: usize,
    nby: usize,
    w: usize,
    h: usize,
    heads: Vec<Vec<u32>>,
    pts: Vec<(f32, f32, f32, f32)>,
}

impl Buckets {
    fn new(sk: &Skeleton, b: usize) -> Self {
        let (w, h) = (sk.width, sk.height);
        let (nbx, nby) = (w.div_ceil(b), h.div_ceil(b));
        let mut heads = vec![Vec::new(); nbx * nby];
        let mut pts = Vec::new();
        for l in &sk.polylines {
            for &p in l {
                let (x, y) = (p.0.rem_euclid(w as f32), p.1.rem_euclid(h as f32));
                let (bx, by) = ((x as usize).min(w - 1) / b, (y as usize).min(h - 1) / b);
                heads[by * nbx + bx].push(pts.len() as u32);
                pts.push((x, y, p.2, p.3));
            }
        }
        Self { b, nbx, nby, w, h, heads, pts }
    }

    /// The sample nearest to the centre of cell `c` (torus): (sample, distance in cells).
    fn nearest(&self, c: usize) -> (usize, f32) {
        let (cx, cy) = ((c % self.w) as f32 + 0.5, (c / self.w) as f32 + 0.5);
        let (bx0, by0) = ((c % self.w) / self.b, (c / self.w) / self.b);
        let mut best = (usize::MAX, f32::INFINITY);
        let rmax = self.nbx.max(self.nby) as i64;
        let mut r = 0i64;
        while r <= rmax {
            if ((r - 1).max(0) as f32) * self.b as f32 > best.1 {
                break;
            }
            for dy in -r..=r {
                for dx in -r..=r {
                    if dx.abs() != r && dy.abs() != r {
                        continue;
                    }
                    let bx = (bx0 as i64 + dx).rem_euclid(self.nbx as i64) as usize;
                    let by = (by0 as i64 + dy).rem_euclid(self.nby as i64) as usize;
                    for &s in &self.heads[by * self.nbx + bx] {
                        let p = self.pts[s as usize];
                        let mut ddx = (cx - p.0).rem_euclid(self.w as f32);
                        if ddx > self.w as f32 / 2.0 {
                            ddx -= self.w as f32;
                        }
                        let mut ddy = (cy - p.1).rem_euclid(self.h as f32);
                        if ddy > self.h as f32 / 2.0 {
                            ddy -= self.h as f32;
                        }
                        let d = (ddx * ddx + ddy * ddy).sqrt();
                        if d < best.1 {
                            best = (s as usize, d);
                        }
                    }
                }
            }
            r += 1;
        }
        best
    }

    /// Finding 126-C's position of cell `c` on its wall, on the PLANAR geometry of the nearest
    /// sample: 0 floor (d ≤ W/2), 1 foot (u < `foot_m` beyond the floor edge), 3 crest (the planar
    /// wall within `crest_zone_m` of the terrain), 2 the straight segment between.
    fn position(&self, c: usize, zpre: &[f32], tan: f32, cell_m: f32, foot_m: f32, crest_zone_m: f32) -> usize {
        let (s, dc) = self.nearest(c);
        if s == usize::MAX {
            return 2;
        }
        let (_, _, zf, hw) = self.pts[s];
        let d = dc * cell_m;
        if d <= hw {
            return 0;
        }
        let u = d - hw;
        if u < foot_m {
            1
        } else if zpre[c] - (zf + tan * u) < crest_zone_m {
            3
        } else {
            2
        }
    }
}

fn shares(v: &[usize; 4]) -> String {
    let t = v.iter().sum::<usize>().max(1) as f32;
    (0..4).map(|i| format!("{} {:.1} %", POS_NAME[i], 100.0 * v[i] as f32 / t)).collect::<Vec<_>>().join(" / ")
}

/// Terrain R8 over fully-land (in PRE) 16×16 windows — Finding 125's `r8_pop` total.
fn r8_terrain(f: &GridF32, cls: &[u8], win: usize) -> f32 {
    let (w, h) = (f.width, f.height);
    let (mut vr, mut vi, mut nwin) = (0f64, 0f64, 0f64);
    for by in (0..h.saturating_sub(win)).step_by(win) {
        'win: for bx in (0..w.saturating_sub(win)).step_by(win) {
            let (mut jxx, mut jyy, mut jxy) = (0f64, 0f64, 0f64);
            for y in by..by + win {
                for x in bx..bx + win {
                    if cls[y * w + x] == 255 {
                        continue 'win;
                    }
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
            vr += (8.0 * t).cos();
            vi += (8.0 * t).sin();
            nwin += 1.0;
        }
    }
    ((vr * vr + vi * vi).sqrt() / nwin.max(1.0)) as f32
}

/// The coast of `g`: its spurs (eroded stage, Finding 75's detector on the u16 mask), each typed
/// PENINSULA or INLET (the cell halfway from the neck to the excursion's farthest sample is land or
/// sea) and placed NEAR (≤ `near` cells) or FAR from a coastal wall, and the coastline length in
/// each of the two zones (km, domain) — the base rate's denominators.
struct Coast {
    mids: Vec<(f32, f32)>,
    near: Vec<bool>,
    inlet: Vec<bool>,
    l_near_km: f32,
    l_far_km: f32,
}

fn coast(g: &GridF32, ss: &SteinSteinParams, dwall: &[u16], near: f32) -> Coast {
    let (w, h) = (g.width, g.height);
    let land = land_u16(g, ss);
    let (mj, ww) = majority(&land, w, 1);
    let polys = marching_squares(&to_mask(&mj, ww), 0.5);
    let (sp, _) = coast_spurs(&polys, CELL_KM, MIN_SPUR_KM, NECK_KM);
    let cell = |p: (f32, f32)| {
        (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize
    };
    let (mut l_near, mut l_far) = (0f32, 0f32);
    for pl in &polys {
        for s in pl.windows(2) {
            let len = ((s[1].0 - s[0].0).powi(2) + (s[1].1 - s[0].1).powi(2)).sqrt() * CELL_KM;
            let mid = (0.5 * (s[0].0 + s[1].0), 0.5 * (s[0].1 + s[1].1));
            if dwall[cell(mid)] as f32 <= near {
                l_near += len;
            } else {
                l_far += len;
            }
        }
    }
    let mut out = Coast { mids: Vec::new(), near: Vec::new(), inlet: Vec::new(), l_near_km: l_near, l_far_km: l_far };
    for s in &sp {
        let pl = &polys[s.poly];
        let idx: Vec<usize> = if s.j >= s.i { (s.i..=s.j).collect() } else { (s.i..pl.len()).chain(0..=s.j).collect() };
        let far = idx
            .iter()
            .map(|&i| pl[i])
            .max_by(|a, b| {
                let da = (a.0 - s.mid.0).powi(2) + (a.1 - s.mid.1).powi(2);
                let db = (b.0 - s.mid.0).powi(2) + (b.1 - s.mid.1).powi(2);
                da.total_cmp(&db)
            })
            .unwrap_or(s.mid);
        let probe = (0.5 * (s.mid.0 + far.0), 0.5 * (s.mid.1 + far.1));
        out.mids.push(s.mid);
        out.near.push(dwall[cell(s.mid)] as f32 <= near);
        out.inlet.push(!mj[cell(probe)]);
    }
    out
}

struct World {
    name: &'static str,
    vc: ValleyConstruction,
    /// the exported network (teeth, positions, F38) — Finding 125's HD assembly
    net: bool,
    /// Finding 108's over-dug class (the canyons)
    canyons: bool,
}

#[test]
#[ignore]
fn f126_bc() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 126-B/C . the coast's base rate, the clause, the teeth on their wall  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let zpre: Vec<f32> = (0..n).map(|k| c1_altitude_norm_to_metres(pre.data[k], &ss)).collect();
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let prof = WallProfile::f124();
    let (foot_m, crest_zone_m) = (prof.foot_m, 3.0 * prof.crest_m);
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let b1a = WallProfile { detail_gain: 0.0, ..prof };
    let b1b = WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof };
    let eps = Some(base.base_m);
    let worlds = [
        World { name: "témoin C2/10 col", vc: base, net: true, canyons: true },
        World { name: "B2 (A_c)", vc: ValleyConstruction { a_min_km2: 0.1, ..base }, net: false, canyons: false },
        World { name: "B1a profile", vc: ValleyConstruction { wall_profile: Some(b1a), ..base }, net: true, canyons: false },
        World { name: "B1b noise", vc: ValleyConstruction { wall_profile: Some(b1b), ..base }, net: true, canyons: false },
        World { name: "B1 both", vc: ValleyConstruction { wall_profile: Some(prof), ..base }, net: true, canyons: false },
        World {
            name: "B1a + B2 (A_c)",
            vc: ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(b1a), ..base },
            net: false,
            canyons: false,
        },
        World {
            name: "B2 (A_c) + B1",
            vc: ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), ..base },
            net: false,
            canyons: false,
        },
        World {
            name: "B1a + CLAUSE",
            vc: ValleyConstruction { wall_profile: Some(b1a), wall_sea_floor_m: eps, ..base },
            net: true,
            canyons: true,
        },
        World {
            name: "B1 + CLAUSE",
            vc: ValleyConstruction { wall_profile: Some(prof), wall_sea_floor_m: eps, ..base },
            net: true,
            canyons: true,
        },
        World {
            name: "B2 (A_c) + B1 + CLAUSE",
            vc: ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), wall_sea_floor_m: eps, ..base },
            net: true,
            canyons: true,
        },
    ];
    let near = 2.0 / CELL_KM; // 2 km (domain) in cells
    let mut g_t: Option<GridF32> = None;
    let mut spurs_t: Vec<(f32, f32)> = Vec::new();
    let mut summary: Vec<String> = Vec::new();
    for wd in worlds {
        let t = Instant::now();
        let g = build_field_seed(Knobs { valley: Some(wd.vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let secs = t.elapsed().as_secs_f64();
        eprintln!("\n────────── {} · build {secs:.0} s ──────────", wd.name);
        let sk = skeleton(&pre, &wd.vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &wd.vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let tan = wd.vc.wall_deg.to_radians().tan();
        // ── the coast and its base rate ──
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw));
        let co = coast(&g, &ss, &dwall, near);
        let nn = co.near.iter().filter(|&&b| b).count();
        let nf = co.mids.len() - nn;
        let is_new = |p: &(f32, f32)| spurs_t.iter().all(|b| ((p.0 - b.0).powi(2) + (p.1 - b.1).powi(2)).sqrt() > near);
        let new: Vec<usize> = if g_t.is_none() { Vec::new() } else { (0..co.mids.len()).filter(|&i| is_new(&co.mids[i])).collect() };
        let inl = co.inlet.iter().filter(|&&b| b).count();
        let new_inl = new.iter().filter(|&&i| co.inlet[i]).count();
        let new_near = new.iter().filter(|&&i| co.near[i]).count();
        eprintln!(
            "   COAST (eroded stage): {} spurs ({} inlets / {} peninsulas) · NEAR a coastal wall (≤ 2 km domain) {nn} over \
             {:.1} km → **{:.4} /km** · elsewhere {nf} over {:.1} km → **{:.4} /km** (km domain)",
            co.mids.len(),
            inl,
            co.mids.len() - inl,
            co.l_near_km,
            nn as f32 / co.l_near_km.max(1e-6),
            co.l_far_km,
            nf as f32 / co.l_far_km.max(1e-6)
        );
        if g_t.is_some() {
            eprintln!(
                "      NEW vs the témoin: {} · near a coastal wall {new_near} · inlets {new_inl} / peninsulas {}",
                new.len(),
                new.len() - new_inl
            );
        }
        // ── R8 terrain (Finding 125's windows), the world and the témoin on the SAME classes ──
        let r8w = r8_terrain(&g, &cls, 16);
        let r8t = r8_terrain(g_t.as_ref().unwrap_or(&g), &cls, 16);
        eprintln!("   R8 terrain {r8w:.4} · témoin on the same classes {r8t:.4}");
        let mut row = format!(
            "{:<24} spurs {:>4} (near {:.4}/km, far {:.4}/km, inlets {inl})",
            wd.name,
            co.mids.len(),
            nn as f32 / co.l_near_km.max(1e-6),
            nf as f32 / co.l_far_km.max(1e-6)
        );
        // ── the clause: what it moved at the construction stage, where ──
        if wd.vc.wall_sea_floor_m.is_some() {
            let (a, _) = carve(&pre, &sk, &wd.vc, &ss);
            let (b, _) = carve(&pre, &sk, &ValleyConstruction { wall_sea_floor_m: None, ..wd.vc }, &ss);
            let moved: Vec<usize> = (0..n).filter(|&k| a.data[k] != b.data[k]).collect();
            let wcb = water_class(&b, SEA);
            let ocean = moved.iter().filter(|&&k| wcb[k] == WATER_CLASS_OCEAN).count();
            let inland = moved.iter().filter(|&&k| wcb[k] == WATER_CLASS_INLAND).count();
            eprintln!(
                "   the CLAUSE moved {} cells at the construction stage (on PRE): below the sea without it {} (edge-connected \
                 OCEAN {ocean}, ENCLOSED {inland}), still land without it {}",
                moved.len(),
                ocean + inland,
                moved.len() - ocean - inland
            );
        }
        // ── the over-dug class (the canyons), located ──
        let drained = (wd.canyons || wd.net).then(|| {
            let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
            let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
            (dd, bre)
        });
        if wd.canyons {
            let bre0 = &drained.as_ref().expect("drained").1;
            let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
            let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
            let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
            set_dump(true);
            let _cr = f95_criteria(&g, bre0, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
            set_dump(false);
            let mut ce = 0usize;
            for b in take_bodies() {
                let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
                let dfill = fill_del[floor] - fill_pre[floor];
                if over_dug_depression(dfill, b.rim) {
                    ce += 1;
                    eprintln!(
                        "   OVER-DUG body {} · {:.2} km² (domain) · floor ({},{}) · Δfill {dfill:.1} m · rim {:.1}° · class {}",
                        b.id,
                        b.km2,
                        floor % w,
                        floor / w,
                        b.rim,
                        CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea")
                    );
                }
            }
            eprintln!("   over-dug bodies (eroded stage): {ce}");
            row += &format!(" · over-dug {ce}");
        }
        // ── the exported network: teeth, where they sit on their wall, or Finding 38 ──
        if wd.net {
            let (dd, bre) = drained.expect("drained");
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
                    let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
                    eprintln!("   teeth {} · R8 network (chord 8) {:.4}", teeth.len(), r8_chords(&pts, 8));
                    row += &format!(" · teeth {} · R8 net {:.4}", teeth.len(), r8_chords(&pts, 8));
                    // C — where on the wall: mouth, head and every point; the wall's own area as the base rate
                    let bk = Buckets::new(&sk, 32);
                    let pos = |c: usize| bk.position(c, &zpre, tan, sk.cell_m, foot_m, crest_zone_m);
                    let (mut mouth, mut head, mut all) = ([0usize; 4], [0usize; 4], [0usize; 4]);
                    for &i in &teeth {
                        let o = own(&segs[i]);
                        let c = |p: &(u32, u32)| p.1 as usize * w + p.0 as usize;
                        mouth[pos(c(o.last().expect("own")))] += 1;
                        head[pos(c(&o[0]))] += 1;
                        for p in o {
                            all[pos(c(p))] += 1;
                        }
                    }
                    let mut area = [0usize; 4];
                    for k in (0..n).step_by(13) {
                        if wall[k] {
                            area[pos(k)] += 1;
                        }
                    }
                    eprintln!("   C · teeth MOUTHS: {}", shares(&mouth));
                    eprintln!("   C · teeth HEADS:  {}", shares(&head));
                    eprintln!("   C · teeth POINTS: {}", shares(&all));
                    eprintln!("   C · the WALL's own area (every 13th wall cell, the base rate): {}", shares(&area));
                }
                Err(e) => {
                    let msg = e
                        .downcast_ref::<String>()
                        .cloned()
                        .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                        .unwrap_or_default();
                    eprintln!("   ⛔ the HD assembly PANICKED: {}", msg.lines().next().unwrap_or(""));
                    row += " · ⛔ Finding 38";
                }
            }
        }
        row += &format!(" · build {secs:.0} s");
        summary.push(row);
        if g_t.is_none() {
            spurs_t = co.mids.clone();
            g_t = Some(g);
        }
    }
    eprintln!("\n   ── summary ──");
    for s in summary {
        eprintln!("   {s}");
    }
    eprintln!("\n==========  end Finding 126-B/C . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

// ════════════════════════════════ Parts D + E ════════════════════════════════
//
// D is NOT Finding 113's snap. Finding 113 moved the EXPORTED river polyline toward a thalweg it
// already sat in (MFD ≡ D8 93 %, Finding 15) and smoothed that drawn line (R8 0.466 → 0.452 at ±8
// cells): the terrain never changed, so nothing the terrain's R8 reads could move. Here the
// SKELETON's geometry is changed BEFORE the construction digs, so the walls themselves move.
//
// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f126_d12 --nocapture
//      cargo test -p ymir-core --release --test f126_coast -- --ignored f126_e --nocapture
//      YMIR_F126_SMOOTH_W=<c> cargo test -p ymir-core --release --test f126_coast -- --ignored f126_d3 --nocapture

/// Per band (index = `band()`), the 8-th harmonic of the skeleton's chords: (count, Σcos 8θ, Σsin 8θ).
/// `Chord::Points(8)` is Finding 125's instrument (8 densified samples, ≈ 190 m); `Chord::Metres(m)`
/// cuts the line every `m` metres of arclength — a law in metres. The band is read at the chord's
/// middle sample.
#[derive(Clone, Copy)]
enum Chord {
    Points(usize),
    Metres(f32),
}

fn chord_r8(sk: &Skeleton, chord: Chord, keep: &dyn Fn(usize, usize) -> bool) -> [(usize, f64, f64); NCLS] {
    let (w, h) = (sk.width, sk.height);
    let mut acc = [(0usize, 0f64, 0f64); NCLS];
    let push = |a: (f32, f32), b: (f32, f32), mid: (f32, f32), acc: &mut [(usize, f64, f64); NCLS]| {
        let (dx, dy) = ((b.0 - a.0) as f64, (b.1 - a.1) as f64);
        if dx == 0.0 && dy == 0.0 {
            return;
        }
        let (x, y) = ((mid.0.floor() as i64).rem_euclid(w as i64) as usize, (mid.1.floor() as i64).rem_euclid(h as i64) as usize);
        let c = band(sk.area_km2[y * w + x]) as usize;
        let t = dy.atan2(dx);
        acc[c].0 += 1;
        acc[c].1 += (8.0 * t).cos();
        acc[c].2 += (8.0 * t).sin();
    };
    for (li, l) in sk.polylines.iter().enumerate() {
        match chord {
            Chord::Points(p) => {
                let mut i = 0usize;
                while i + p < l.len() {
                    if keep(li, i + p / 2) {
                        push((l[i].0, l[i].1), (l[i + p].0, l[i + p].1), (l[i + p / 2].0, l[i + p / 2].1), &mut acc);
                    }
                    i += p;
                }
            }
            Chord::Metres(m) => {
                let step = m / sk.cell_m;
                let (mut start, mut arc) = (0usize, 0f32);
                let mut marks: Vec<(usize, f32)> = vec![(0, 0.0)];
                for i in 1..l.len() {
                    arc += ((l[i].0 - l[i - 1].0).powi(2) + (l[i].1 - l[i - 1].1).powi(2)).sqrt();
                    marks.push((i, arc));
                    if arc - marks[start].1 >= step {
                        let half = marks[start].1 + 0.5 * (arc - marks[start].1);
                        let mid = marks[start..].iter().find(|&&(_, a)| a >= half).map_or(i, |&(j, _)| j);
                        if keep(li, mid) {
                            push((l[marks[start].0].0, l[marks[start].0].1), (l[i].0, l[i].1), (l[mid].0, l[mid].1), &mut acc);
                        }
                        start = marks.len() - 1;
                    }
                }
            }
        }
    }
    acc
}

fn r8_row(acc: &[(usize, f64, f64); NCLS]) -> String {
    (1..NCLS)
        .map(|c| {
            let (k, cr, ci) = acc[c];
            if k < 40 {
                format!("{} — ({k})", CLS_NAME[c])
            } else {
                format!("{} **{:.3}** ({k})", CLS_NAME[c], ((cr * cr + ci * ci).sqrt() / k as f64) as f32)
            }
        })
        .collect::<Vec<_>>()
        .join(" · ")
}

/// A provably ISOTROPIC relief on a DOME that drains to a surrounding sea, in metres: 256 plane waves
/// of uniform random direction and phase, wavelengths log-uniform in [4, 512] cells, amplitude
/// ∝ λ (every scale adds the same slope), scaled so the noise's RMS gradient is `rho` × the dome's
/// slope. Finding 97's isotropic synthetic is the same family (random-direction plane waves); the
/// spectrum here reaches the cell scale because the D8 tracé is decided there.
fn synth_dome(n: usize, cell_m: f32, rho: f32, seed: u64) -> Vec<f32> {
    let mut st = seed;
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let nw = 256usize;
    let r_dome = 0.45 * n as f32;
    let top = 1500.0f32;
    let slope = top / (r_dome * cell_m);
    let alpha = rho * slope * cell_m / (std::f32::consts::TAU * (nw as f32 / 2.0).sqrt());
    let waves: Vec<(f32, f32, f32, f32)> = (0..nw)
        .map(|_| {
            let ang = (rnd() * std::f64::consts::TAU) as f32;
            let lam = 4.0f32 * (128.0f32).powf(rnd() as f32);
            let k = std::f32::consts::TAU / lam;
            (k * ang.cos(), k * ang.sin(), (rnd() * std::f64::consts::TAU) as f32, alpha * lam)
        })
        .collect();
    let c = 0.5 * n as f32;
    let mut z = vec![0f32; n * n];
    z.par_iter_mut().enumerate().for_each(|(k, v)| {
        let (x, y) = ((k % n) as f32, (k / n) as f32);
        let r = ((x - c).powi(2) + (y - c).powi(2)).sqrt();
        if r >= r_dome {
            *v = -200.0;
            return;
        }
        let mut s = top * (1.0 - r / r_dome);
        for &(kx, ky, ph, a) in &waves {
            s += a * (kx * x + ky * y + ph).sin();
        }
        *v = s.max(1.0);
    });
    z
}

/// The roughness ratio the D8 tracé answers to: RMS |∇(z − trend)| / RMS |∇trend| over land, the
/// trend a box blur of ±1 km (20 cells). The same definition for a synthetic and for the real field,
/// so the D8 floor is read at the real terrain's own ratio.
fn roughness(zm: &[f32], land: &[bool], w: usize, h: usize) -> f32 {
    let det = detail(zm, w, h, 20);
    let (mut sd, mut st, mut k2) = (0f64, 0f64, 0usize);
    for y in (1..h - 1).step_by(7) {
        for x in (1..w - 1).step_by(7) {
            let k = y * w + x;
            if !(land[k] && land[k - 1] && land[k + 1] && land[k - w] && land[k + w]) {
                continue;
            }
            let g = |f: &dyn Fn(usize) -> f32| -> f64 {
                let (gx, gy) = (0.5 * (f(k + 1) - f(k - 1)), 0.5 * (f(k + w) - f(k - w)));
                (gx * gx + gy * gy) as f64
            };
            sd += g(&|i| det[i]);
            st += g(&|i| zm[i] - det[i]);
            k2 += 1;
        }
    }
    if k2 == 0 || st <= 0.0 { f32::NAN } else { (sd / st).sqrt() as f32 }
}

/// ADR Finding 126-D1 and D2 — the D8 floor of the skeleton, then the smoothing length.
#[test]
#[ignore]
fn f126_d12() {
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 126-D1/D2 . the skeleton's D8 floor, and the smoothing length  ==========");
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let b2 = ValleyConstruction { a_min_km2: 1.0, ..base };
    let variants: [(&str, ValleyConstruction); 5] = [
        ("raw D8 (no smoothing)", ValleyConstruction { smooth_m: 0.0, ..b2 }),
        ("L = 250 m (the existing, pinned)", b2),
        ("L = 0.5 W", ValleyConstruction { smooth_w: Some(0.5), ..b2 }),
        ("L = 1 W", ValleyConstruction { smooth_w: Some(1.0), ..b2 }),
        ("L = 2 W", ValleyConstruction { smooth_w: Some(2.0), ..b2 }),
    ];
    let all = |_: usize, _: usize| true;
    // ── D1: the floor, on isotropic synthetics at the SAME cell size (2048² over 100 km) ──
    let ns = 2048usize;
    let cell_m = CELL_KM * 1000.0;
    let dom_s = ns as f32 * CELL_KM;
    eprintln!("\n-- D1 . isotropic synthetics, {ns}² at {cell_m:.2} m (a {dom_s:.0} km domain), dome 1500 m over {:.1} km --", 0.45 * dom_s);
    for rho in [0.1f32, 0.3, 1.0, 3.0] {
        let zm = synth_dome(ns, cell_m, rho, 0x9E37_79B9_7F4A_7C15);
        let g = GridF32 { width: ns, height: ns, data: zm.iter().map(|&m| c1_metres_to_altitude_norm(m, &ss)).collect() };
        let land: Vec<bool> = g.data.iter().map(|&v| v > SEA).collect();
        let a = common::aniso(&g, &land, 16);
        eprintln!(
            "   ρ = {rho} (noise RMS gradient / dome slope) · MEASURED roughness ratio {:.3} · terrain R8 (w 16) **{:.4}**              (Finding 97's floor 0.0402) · R2 {:.4}",
            roughness(&zm, &land, ns, ns),
            a.r8,
            a.r2
        );
        for (label, vc) in &variants[..2] {
            let sk = skeleton(&g, vc, &ss, dom_s);
            eprintln!("      {label:<34} 8-pt: {}", r8_row(&chord_r8(&sk, Chord::Points(8), &all)));
            eprintln!("      {:<34} 200 m: {}", "", r8_row(&chord_r8(&sk, Chord::Metres(200.0), &all)));
        }
    }
    eprintln!("   (D1 in {:.0} s)", t0.elapsed().as_secs_f64());

    // ── D2: the real skeletons (PRE = Finding 125's instrument; S1 = the pipeline's input) ──
    let t = Instant::now();
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    eprintln!("\n-- D2 . the real skeleton at 1 km² (B2), R8 by band against L (L is the HALF-length, like `smooth_m`) · builds {:.0} s --", t.elapsed().as_secs_f64());
    for (fname, f) in [("PRE (Finding 125's instrument)", &pre), ("S1 (the pipeline's construction input)", &s1)] {
        let zm: Vec<f32> = f.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let land: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
        eprintln!("   on {fname} · MEASURED roughness ratio {:.3}:", roughness(&zm, &land, f.width, f.height));
        drop((zm, land));
        let raw = skeleton(f, &variants[0].1, &ss, DOMAIN_KM);
        // which raw samples the EXISTING rule smooths: a line needs ≥ 2·5 + 3 raw samples, and the
        // first / last 5 are pinned. Raw samples are the ones at cell centres in the unsmoothed line.
        let half = (250.0 / raw.cell_m).round() as usize;
        let covered: Vec<Vec<bool>> = raw
            .polylines
            .iter()
            .map(|l| {
                let is_raw: Vec<bool> = l.iter().map(|p| (p.0.fract() - 0.5).abs() < 1e-4 && (p.1.fract() - 0.5).abs() < 1e-4).collect();
                let nraw = is_raw.iter().filter(|&&b| b).count();
                let mut ri = 0usize;
                is_raw
                    .iter()
                    .map(|&r| {
                        let idx = ri;
                        if r {
                            ri += 1;
                        }
                        nraw >= 2 * half + 3 && idx >= half && idx + half < nraw
                    })
                    .collect()
            })
            .collect();
        let cov = |li: usize, i: usize| covered[li][i];
        let unc = |li: usize, i: usize| !covered[li][i];
        let r_cov = chord_r8(&raw, Chord::Points(8), &cov);
        let r_unc = chord_r8(&raw, Chord::Points(8), &unc);
        let share = |c: usize| 100.0 * r_unc[c].0 as f32 / (r_cov[c].0 + r_unc[c].0).max(1) as f32;
        eprintln!(
            "      the EXISTING rule leaves RAW: {:.1} % of the ≥ 10 chords, {:.1} % of 3–10, {:.1} % of 1–3 (short lines + pinned ends)",
            share(1),
            share(2),
            share(3)
        );
        eprintln!("      raw chords the existing rule WOULD smooth · 8-pt: {}", r8_row(&r_cov));
        eprintln!("      raw chords it leaves RAW                  · 8-pt: {}", r8_row(&r_unc));
        // the trunks' displacement: each raw trunk cell ≥ 10 km² to the nearest smoothed sample
        let trunk10: Vec<usize> = (0..raw.trunk.len()).filter(|&k| raw.trunk[k] && raw.area_km2[k] >= 10.0).step_by(7).collect();
        for (label, vc) in &variants {
            let sk = skeleton(f, vc, &ss, DOMAIN_KM);
            let bk = Buckets::new(&sk, 32);
            let mut disp: Vec<f32> = trunk10.iter().map(|&k| bk.nearest(k).1 * sk.cell_m).collect();
            disp.sort_by(f32::total_cmp);
            let q = |p: f32| disp[((disp.len() - 1) as f32 * p) as usize];
            eprintln!("      {label:<34} 8-pt: {}", r8_row(&chord_r8(&sk, Chord::Points(8), &all)));
            eprintln!("      {:<34} 200 m: {}", "", r8_row(&chord_r8(&sk, Chord::Metres(200.0), &all)));
            eprintln!(
                "      {:<34} trunks ≥ 10 km² moved: p50 {:.1} m · p90 {:.1} m · max {:.1} m (from the raw D8 cell centre)",
                "",
                q(0.5),
                q(0.9),
                q(1.0)
            );
        }
    }
    eprintln!("\n==========  end Finding 126-D1/D2 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// Primitives (skeleton polylines) whose planar cone reaches below the field at cell `c` — Finding
/// 125's instrument, verbatim.
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

/// ADR Finding 126-E — the cone count's reference: 200 random floor cells per population.
#[test]
#[ignore]
fn f126_e() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 126-E . the canyons' cone count, against its reference  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let zpre: Vec<f32> = (0..n).map(|k| c1_altitude_norm_to_metres(pre.data[k], &ss)).collect();
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let canyons: [(usize, usize); 4] = [(3892, 2310), (4588, 3240), (2281, 4563), (1835, 4562)];
    for (label, vc, pops) in [
        ("témoin C2/10 col (the asked reference)", base, vec![1u8]),
        ("B2 → A_c (the canyons' own world)", ValleyConstruction { a_min_km2: 0.1, ..base }, vec![1u8, 4]),
    ] {
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let tan = vc.wall_deg.to_radians().tan();
        eprintln!("\n   {label} · {} polylines", sk.polylines.len());
        let at: Vec<String> = canyons
            .iter()
            .map(|&(x, y)| {
                let c = y * w + x;
                format!("({x},{y}) {} [{}]", overlaps(&sk, &zpre, c, tan), CLS_NAME.get(cls[c] as usize).copied().unwrap_or("sea"))
            })
            .collect();
        eprintln!("      at the four canyon floors: {}", at.join(" · "));
        for pop in pops {
            let cand: Vec<usize> = (0..n).filter(|&k| mk.floor[k] && cls[k] == pop).collect();
            let mut st = 0x2545_F491_4F6C_DD1Du64;
            let mut pick: Vec<usize> = (0..200)
                .map(|_| {
                    st ^= st << 13;
                    st ^= st >> 7;
                    st ^= st << 17;
                    cand[(st % cand.len() as u64) as usize]
                })
                .collect();
            pick.sort_unstable();
            let mut v: Vec<usize> = pick.par_iter().map(|&c| overlaps(&sk, &zpre, c, tan)).collect();
            v.sort_unstable();
            let q = |p: f32| v[((v.len() - 1) as f32 * p) as usize];
            let rank = |x: usize| 100.0 * v.iter().filter(|&&a| a < x).count() as f32 / v.len() as f32;
            let ranks: Vec<String> = canyons
                .iter()
                .map(|&(x, y)| {
                    let k = overlaps(&sk, &zpre, y * w + x, tan);
                    format!("{k} → p{:.0}", rank(k))
                })
                .collect();
            eprintln!(
                "      200 random {} FLOOR cells (of {}): p10 {} · p50 {} · p90 {} · p99 {} · max {} · the canyons' ranks: {}",
                CLS_NAME[pop as usize],
                cand.len(),
                q(0.1),
                q(0.5),
                q(0.9),
                q(0.99),
                q(1.0),
                ranks.join(" · ")
            );
        }
    }
    eprintln!("\n==========  end Finding 126-E . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// R8 terrain decomposed by class — Finding 125's `r8_pop`, verbatim: (total, per class (windows,
/// class R8, contribution)).
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

/// Finding 124-P4's skeleton column: the share of PRE's 1 km² trunk cells within 5 cells of the
/// world's, and PRE's 1 km² mouths matched within 2 cells (both on breached D8 flows).
fn skeleton_pairing(pre: &GridF32, g: &GridF32, ss: &SteinSteinParams) -> (f32, usize, usize) {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let (w, h) = (g.width, g.height);
    let n = w * h;
    let flow_of = |f: &GridF32| {
        let d = c1_drainage_windowed(f, None, &dcfg(), ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, SEA, w, h);
        let fl = c1_drainage_windowed(&bf, None, &dcfg(), ss, DOMAIN_KM).flow;
        (bf, fl)
    };
    let (ba, fa) = flow_of(pre);
    let (bb, fb) = flow_of(g);
    let ac = 1.0 / (CELL_KM * CELL_KM);
    let ta: Vec<bool> = (0..n).map(|k| ba.data[k] > SEA && fa.accumulation.data[k] >= ac).collect();
    let tb: Vec<bool> = (0..n).map(|k| bb.data[k] > SEA && fb.accumulation.data[k] >= ac).collect();
    let db = dist_from(&tb, w, h);
    let v: Vec<u16> = (0..n).filter(|&k| ta[k]).map(|k| db[k]).collect();
    let within = 100.0 * v.iter().filter(|&&d| d <= 5).count() as f32 / v.len().max(1) as f32;
    let mouths = |f: &ymir_core::terrain::flow::FlowResult, b: &GridF32| -> Vec<usize> {
        (0..n)
            .filter(|&k| b.data[k] > SEA && f.accumulation.data[k] >= ac && f.direction[k] != DIR_NONE)
            .filter(|&k| {
                let d = f.direction[k] as usize;
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                b.data[ny * w + nx] <= SEA
            })
            .collect()
    };
    let (ma, mb) = (mouths(&fa, &ba), mouths(&fb, &bb));
    let mut src = vec![false; n];
    for &c in &mb {
        src[c] = true;
    }
    let md = dist_from(&src, w, h);
    (within, ma.iter().filter(|&&c| md[c] <= 2).count(), ma.len())
}

/// ADR Finding 126-D3 — B2 (1 km²) built on the skeleton smoothed at the best L of D2, against the
/// témoin and B2 as Findings 124/125 built it. `YMIR_F126_SMOOTH_W` = c (L = c·W).
#[test]
#[ignore]
fn f126_d3() {
    let c: f32 = std::env::var("YMIR_F126_SMOOTH_W").expect("YMIR_F126_SMOOTH_W").parse().expect("a number");
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 126-D3 . B2 on the skeleton smoothed at L = {c} W  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let spurs_of = |g: &GridF32| -> usize {
        let (m, ww) = majority(&land_u16(g, &ss), w, 1);
        coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
    };
    let auth = spurs_of(&pre);
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let b2 = ValleyConstruction { a_min_km2: 1.0, ..base };
    let worlds: [(&str, ValleyConstruction); 3] = [
        ("témoin C2/10 col", base),
        ("B2 ≥ 1 km² (L = 250 m, Findings 124/125)", b2),
        ("B2 ≥ 1 km², SMOOTHED (L = c·W)", ValleyConstruction { smooth_w: Some(c), ..b2 }),
    ];
    let mut g_t: Option<GridF32> = None;
    for (label, vc) in worlds {
        let t = Instant::now();
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let secs = t.elapsed().as_secs_f64();
        eprintln!("\n────────── {label} · build {secs:.0} s ──────────");
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let all = |_: usize, _: usize| true;
        eprintln!("   the skeleton (on PRE) · 8-pt: {}", r8_row(&chord_r8(&sk, Chord::Points(8), &all)));
        let (r8w, popw) = r8_pop(&g, &cls, 16);
        let (r8t, popt) = r8_pop(g_t.as_ref().unwrap_or(&g), &cls, 16);
        eprintln!("   R8 terrain: world {r8w:.4} · témoin (same classes) {r8t:.4} · Δ {:+.4}", r8w - r8t);
        for k in 0..NCLS {
            if popw[k].0 < 0.5 {
                continue;
            }
            eprintln!(
                "      {:<11} windows {:>8.1} · class R8 world {:.4} / témoin {:.4} · contribution world {:+.4} / témoin {:+.4} · Δ {:+.4}",
                CLS_NAME[k], popw[k].0, popw[k].1, popt[k].1, popw[k].2, popt[k].2, popw[k].2 - popt[k].2
            );
        }
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let coast = spurs_of(&bre) as i64 - auth as i64;
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut ce = 0usize;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            let dfill = fill_del[floor] - fill_pre[floor];
            if over_dug_depression(dfill, b.rim) {
                ce += 1;
                eprintln!(
                    "   OVER-DUG body {} · {:.2} km² (domain) · floor ({},{}) · Δfill {dfill:.1} m · rim {:.1}° · class {}",
                    b.id,
                    b.km2,
                    floor % w,
                    floor / w,
                    b.rim,
                    CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea")
                );
            }
        }
        drop(fill_del);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        let segs = &dr.rivers.segments;
        let wall: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k]).collect();
        let wc: Vec<usize> = (0..segs.len())
            .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2)
            .collect();
        let teeth = wc
            .iter()
            .filter(|&&i| {
                let p = &segs[i].points;
                p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32 >= 0.8 * p.len() as f32
            })
            .count();
        let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
        let r8n = r8_chords(&pts, 8);
        drop(dr);
        drop(bre);
        let (within, matched, nm) = skeleton_pairing(&pre, &g, &ss);
        eprintln!(
            "   teeth {teeth} · R8 network (chord 8) {r8n:.4} · over-dug {ce} · coast {coast:+} (spurs vs PRE) · skeleton ≤ 5 cells \
             {within:.1} % · mouths {matched}/{nm} · build {secs:.0} s"
        );
        if g_t.is_none() {
            g_t = Some(g);
        }
    }
    eprintln!("\n==========  end Finding 126-D3 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 126-A2 — which half of the PROFILE drowns a wall: the concave foot or the convex
/// crest (the author's "profil concave au pied sous le fond ?"), against the noise, by switching the
/// terms one at a time on the construction's input (S1), and where the drowned cells lie: the input
/// terrain's height there and their distance to the open ocean.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f126_a2 --nocapture
#[test]
#[ignore]
fn f126_a2() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 126-A2 . foot, crest or noise: which term drowns a wall  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let s1_land: Vec<bool> = s1.data.iter().map(|&v| v > SEA).collect();
    let wc1 = water_class(&s1, SEA);
    let ocean1: Vec<bool> = wc1.iter().map(|&c| c == WATER_CLASS_OCEAN).collect();
    let docean = dist_from(&ocean1, w, h);
    drop((wc1, ocean1));
    let idx: Vec<usize> = CELLS.iter().map(|&(x, y)| y * w + x).collect();
    let prof = WallProfile::f124();
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let terms: [(&str, Option<WallProfile>); 6] = [
        ("planar", None),
        ("foot only (crest 0, gain 0)", Some(WallProfile { crest_m: 0.0, detail_gain: 0.0, ..prof })),
        ("crest only (foot 0, gain 0)", Some(WallProfile { foot_m: 0.0, detail_gain: 0.0, ..prof })),
        ("noise only (foot 0, crest 0)", Some(WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof })),
        ("profile = foot + crest (gain 0)", Some(WallProfile { detail_gain: 0.0, ..prof })),
        ("all three (B1)", Some(prof)),
    ];
    for (skname, vc0) in [("C2/10 col skeleton (≥ 10 km²)", base), ("A_c skeleton (≥ 0.1 km², B2)", ValleyConstruction { a_min_km2: 0.1, ..base })] {
        let sk = skeleton(&s1, &vc0, &ss, DOMAIN_KM);
        eprintln!("\n   {skname}:");
        for (label, wp) in terms {
            let (g, mk) = carve(&s1, &sk, &ValleyConstruction { wall_profile: wp, ..vc0 }, &ss);
            let wc = water_class(&g, SEA);
            let drowned: Vec<usize> = (0..n).filter(|&k| s1_land[k] && g.data[k] <= SEA).collect();
            let walls = drowned.iter().filter(|&&k| mk.carved[k] && !mk.floor[k]).count();
            let ocean = drowned.iter().filter(|&&k| wc[k] == WATER_CLASS_OCEAN).count();
            let mut z: Vec<f32> = drowned.iter().map(|&k| c1_altitude_norm_to_metres(s1.data[k], &ss)).collect();
            let mut d: Vec<f32> = drowned.iter().map(|&k| docean[k] as f32 * CELL_KM).collect();
            z.sort_by(f32::total_cmp);
            d.sort_by(f32::total_cmp);
            let q = |v: &[f32], p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p) as usize] };
            let at: Vec<String> = idx.iter().map(|&k| format!("{:.2} m", c1_altitude_norm_to_metres(g.data[k], &ss))).collect();
            eprintln!(
                "      {label:<32} drowned {:>5} (walls {walls}, edge-connected OCEAN {ocean}, enclosed {}) · input terrain there \
                 p50 {:.1} / p90 {:.1} m · to the open ocean p50 {:.2} / p90 {:.2} km (domain) · at the two cells {}",
                drowned.len(),
                drowned.len() - ocean,
                q(&z, 0.5),
                q(&z, 0.9),
                q(&d, 0.5),
                q(&d, 0.9),
                at.join(" · ")
            );
        }
    }
    eprintln!("\n==========  end Finding 126-A2 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 126, rule 13 on this round's own code — the two gated fields left `None` (and the
/// smoothing moved into `smooth_positions`) leave the definition's world BIT-IDENTICAL: a fresh,
/// uncached build of "C2 /10 col (défaut)" must hash to the reference `bench_field_hashes.json` holds.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f126_identity --nocapture
#[test]
#[ignore]
fn f126_identity() {
    use ymir_core::tectonics_c1::bench_guard::field_hash;
    let vc = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let h = format!("{:016x}", field_hash(&g));
    eprintln!("   C2 /10 col (défaut), fresh build: field hash {h} (reference a8d2d538d692c2f0)");
    assert_eq!(h, "a8d2d538d692c2f0", "the gated-off fields moved the definition's world");
}

/// ADR Finding 126-D1, the instrument corrected — the roughness ratio with its ±1 km trend window
/// ENTIRELY on land (every cell ≥ 21 cells from the sea), so no coast enters the blur. The first
/// version let the sea's −200 m step into the trend near every coast and read 0.455 on a synthetic
/// built at ρ = 0.1: contaminated, and not read.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f126_rough --nocapture
#[test]
#[ignore]
fn f126_rough() {
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    let ss = SteinSteinParams::default();
    let rough = |zm: &[f32], land: &[bool], w: usize, h: usize| -> (f32, usize) {
        let sea: Vec<bool> = land.iter().map(|&b| !b).collect();
        let dsea = dist_from(&sea, w, h);
        let det = detail(zm, w, h, 20);
        let (mut sd, mut st, mut k2) = (0f64, 0f64, 0usize);
        for y in (1..h - 1).step_by(7) {
            for x in (1..w - 1).step_by(7) {
                let k = y * w + x;
                if dsea[k] < 22 {
                    continue;
                }
                let g = |f: &dyn Fn(usize) -> f32| -> f64 {
                    let (gx, gy) = (0.5 * (f(k + 1) - f(k - 1)), 0.5 * (f(k + w) - f(k - w)));
                    (gx * gx + gy * gy) as f64
                };
                sd += g(&|i| det[i]);
                st += g(&|i| zm[i] - det[i]);
                k2 += 1;
            }
        }
        (if k2 == 0 || st <= 0.0 { f32::NAN } else { (sd / st).sqrt() as f32 }, k2)
    };
    eprintln!("\n==========  Finding 126-D1 . the roughness ratio, the window on land only  ==========");
    let ns = 2048usize;
    let cell_m = CELL_KM * 1000.0;
    for rho in [0.1f32, 0.3, 1.0, 3.0] {
        let zm = synth_dome(ns, cell_m, rho, 0x9E37_79B9_7F4A_7C15);
        let g: Vec<f32> = zm.iter().map(|&m| c1_metres_to_altitude_norm(m, &ss)).collect();
        let land: Vec<bool> = g.iter().map(|&v| v > SEA).collect();
        let (r, k) = rough(&zm, &land, ns, ns);
        eprintln!("   synthetic ρ = {rho}: roughness ratio **{r:.3}** ({k} sample cells)");
    }
    for (label, k) in [
        ("PRE (no incision, Finding 125's instrument)", Knobs::no_incision()),
        ("S1 (the construction's input)", Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }),
    ] {
        let f = build_field_seed(k, PSEED);
        let zm: Vec<f32> = f.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let land: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
        let (r, k2) = rough(&zm, &land, f.width, f.height);
        eprintln!("   {label}: roughness ratio **{r:.3}** ({k2} sample cells)");
    }
}

// ════════════════════════════════ Finding 127 ════════════════════════════════
//
// The same file holds Finding 127's benches: they read Finding 126's instruments (the coast, the
// buckets, the classes, the synthetics), which a second file would have to copy.
//
// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f127_b0 --nocapture
//      cargo test -p ymir-core --release --test f126_coast -- --ignored f127_b12 --nocapture

/// Finding 127's chord R8 by band, the band read at the NEAREST TRUNK CELL (≤ 3 cells) of the chord's
/// middle sample — an off-grid sample would otherwise read a hillslope's area. Chords with no trunk
/// cell within 3 cells are not counted.
fn chord_r8_nt(sk: &Skeleton, chord: Chord) -> [(usize, f64, f64); NCLS] {
    chord_r8_polys(sk, &sk.polylines, chord)
}

/// [`chord_r8_nt`] over ANY polylines, the band read on `sk`'s trunk cells.
fn chord_r8_polys(sk: &Skeleton, polys: &[Vec<(f32, f32, f32, f32)>], chord: Chord) -> [(usize, f64, f64); NCLS] {
    let (w, h) = (sk.width, sk.height);
    let band_at = |mid: (f32, f32)| -> Option<usize> {
        let (mx, my) = (mid.0.rem_euclid(w as f32), mid.1.rem_euclid(h as f32));
        let (ix, iy) = (mx.floor() as i64, my.floor() as i64);
        let mut best: Option<(f32, usize)> = None;
        for dy in -3i64..=3 {
            for dx in -3i64..=3 {
                let c = (iy + dy).rem_euclid(h as i64) as usize * w + (ix + dx).rem_euclid(w as i64) as usize;
                if !sk.trunk[c] {
                    continue;
                }
                let d = ((dx as f32 + 0.5 - (mx - ix as f32)).powi(2) + (dy as f32 + 0.5 - (my - iy as f32)).powi(2)).sqrt();
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, c));
                }
            }
        }
        best.map(|(_, c)| band(sk.area_km2[c]) as usize)
    };
    let mut acc = [(0usize, 0f64, 0f64); NCLS];
    let push = |a: (f32, f32), b: (f32, f32), mid: (f32, f32), acc: &mut [(usize, f64, f64); NCLS]| {
        let (dx, dy) = ((b.0 - a.0) as f64, (b.1 - a.1) as f64);
        if dx == 0.0 && dy == 0.0 {
            return;
        }
        let Some(c) = band_at(mid) else { return };
        let t = dy.atan2(dx);
        acc[c].0 += 1;
        acc[c].1 += (8.0 * t).cos();
        acc[c].2 += (8.0 * t).sin();
    };
    for l in polys {
        match chord {
            Chord::Points(p) => {
                let mut i = 0usize;
                while i + p < l.len() {
                    push((l[i].0, l[i].1), (l[i + p].0, l[i + p].1), (l[i + p / 2].0, l[i + p / 2].1), &mut acc);
                    i += p;
                }
            }
            Chord::Metres(m) => {
                let step = m / sk.cell_m;
                let (mut i0, mut arc) = (0usize, 0f32);
                let mut arcs = vec![0f32];
                for i in 1..l.len() {
                    arc += ((l[i].0 - l[i - 1].0).powi(2) + (l[i].1 - l[i - 1].1).powi(2)).sqrt();
                    arcs.push(arc);
                    if arc - arcs[i0] >= step {
                        let half = arcs[i0] + 0.5 * (arc - arcs[i0]);
                        let mid = (i0..=i).find(|&j| arcs[j] >= half).unwrap_or(i);
                        push((l[i0].0, l[i0].1), (l[i].0, l[i].1), (l[mid].0, l[mid].1), &mut acc);
                        i0 = i;
                    }
                }
            }
        }
    }
    acc
}

fn r8_of(acc: &[(usize, f64, f64); NCLS], c: usize) -> f32 {
    let (k, cr, ci) = acc[c];
    if k < 40 { f32::NAN } else { ((cr * cr + ci * ci).sqrt() / k as f64) as f32 }
}

/// DIAGNOSTIC — the R8 of the traced PATHS themselves (all outcomes; and the kept ones), bands on the
/// trunk cells: does the descent print its own lattice?
fn paths_line(sk: &Skeleton) -> String {
    let Some(st) = &sk.trace_stats else { return String::new() };
    let conv = |keep: bool| -> Vec<Vec<(f32, f32, f32, f32)>> {
        st.paths.iter().filter(|(k, _)| !keep || *k).map(|(_, p)| p.iter().map(|&(x, y)| (x, y, 0.0, 0.0)).collect()).collect()
    };
    let row = |polys: &[Vec<(f32, f32, f32, f32)>]| {
        let a = chord_r8_polys(sk, polys, Chord::Points(8));
        format!("3–10 **{:.3}** ({}) · 1–3 **{:.3}** ({})", r8_of(&a, 2), a[2].0, r8_of(&a, 3), a[3].0)
    };
    format!("the traced PATHS' own R8: all {} · kept only {}", row(&conv(false)), row(&conv(true)))
}

fn trace_row(sk: &Skeleton) -> String {
    let Some(st) = &sk.trace_stats else { return "D8 (no tracé)".into() };
    let mut dv = st.dev_w.clone();
    dv.sort_by(f32::total_cmp);
    let q = |p: f32| if dv.is_empty() { f32::NAN } else { dv[((dv.len() - 1) as f32 * p) as usize] };
    format!(
        "retraced {} · SAME receiver {} ({:.1} %) · other {} · lost {} · flat fallback {:.2} % of steps · deviation from D8 p50 {:.2} W / p90 {:.2} W",
        st.segments,
        st.same_receiver,
        100.0 * st.same_receiver as f32 / st.segments.max(1) as f32,
        st.other_receiver,
        st.lost,
        100.0 * st.flat_steps as f32 / (st.steps + st.flat_steps).max(1) as f32,
        q(0.5),
        q(0.9)
    )
}

fn r8_line(sk: &Skeleton) -> String {
    let (a8, am) = (chord_r8_nt(sk, Chord::Points(8)), chord_r8_nt(sk, Chord::Metres(200.0)));
    (1..NCLS - 1)
        .map(|c| format!("{} **{:.3}** / {:.3} ({})", CLS_NAME[c], r8_of(&a8, c), r8_of(&am, c), a8[c].0))
        .collect::<Vec<_>>()
        .join(" · ")
}

/// ADR Finding 127-B0 — the off-grid tracé's own calibration, BLOCKING: on Finding 126's isotropic
/// synthetics, at the four roughnesses, the bicubic tracé must read a skeleton R8 ≤ 0.08 in the 1–3
/// and 3–10 km² bands (8-point chords), where D8 reads 0.24–0.99. The bilinear tracé is the negative
/// control (its gradient is piecewise constant per cell).
#[test]
#[ignore]
fn f127_b0() {
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    use ymir_core::tectonics_c1::valley_construction::{SkeletonTrace, TraceInterp};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 127-B0 . the off-grid tracé's calibration (blocking)  ==========");
    eprintln!("   R8 by band: 8-pt chord **bold** / 200 m chord (8-pt chord count); band at the nearest trunk cell");
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let b2 = ValleyConstruction { a_min_km2: 1.0, ..base };
    let ns = 2048usize;
    let cell_m = CELL_KM * 1000.0;
    let dom_s = ns as f32 * CELL_KM;
    let mut worst = 0f32;
    for rho in [0.1f32, 0.3, 1.0, 3.0] {
        let zm = synth_dome(ns, cell_m, rho, 0x9E37_79B9_7F4A_7C15);
        let g = GridF32 { width: ns, height: ns, data: zm.iter().map(|&m| c1_metres_to_altitude_norm(m, &ss)).collect() };
        eprintln!("\n   synthetic ρ = {rho}:");
        for (label, vc) in [
            ("D8 raw", ValleyConstruction { smooth_m: 0.0, ..b2 }),
            ("D8 + 250 m (the existing)", b2),
            ("TRACÉ bicubic", ValleyConstruction { skeleton_trace: Some(SkeletonTrace::f127(TraceInterp::Bicubic)), ..b2 }),
            ("TRACÉ bilinear (neg. control)", ValleyConstruction { skeleton_trace: Some(SkeletonTrace::f127(TraceInterp::Bilinear)), ..b2 }),
        ] {
            let t = Instant::now();
            let sk = skeleton(&g, &vc, &ss, dom_s);
            eprintln!("      {label:<30} {} · {:.1} s", r8_line(&sk), t.elapsed().as_secs_f64());
            if sk.trace_stats.is_some() {
                eprintln!("      {:<30} {}", "", trace_row(&sk));
                eprintln!("      {:<30} {}", "", paths_line(&sk));
            }
            if label.starts_with("TRACÉ bicubic") {
                let a8 = chord_r8_nt(&sk, Chord::Points(8));
                for c in [2usize, 3] {
                    let v = r8_of(&a8, c);
                    worst = worst.max(if v.is_nan() { f32::INFINITY } else { v });
                }
            }
        }
    }
    eprintln!(
        "\n   B0 VERDICT: the bicubic tracé's worst 1–3 / 3–10 km² R8 over the four roughnesses = {worst:.3} → {}",
        if worst <= 0.08 { "**PASS** (≤ 0.08)" } else { "**FAIL** (> 0.08): STOP, no B1–B3" }
    );
    eprintln!("\n==========  end Finding 127-B0 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 127-B1 and B2 — the tracé on the real skeleton (B2 at 1 km²): its topology (the share
/// of retraced segments that reach their D8 receiver first), its R8 against D1's floor, its deviation
/// from the D8 path in W, its flat fallback. On S1 (the construction's input) and on PRE (Findings
/// 125–126's instrument).
#[test]
#[ignore]
fn f127_b12() {
    use ymir_core::tectonics_c1::valley_construction::{SkeletonTrace, TraceInterp};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 127-B1/B2 . the tracé on the real skeleton  ==========");
    eprintln!("   R8 by band: 8-pt chord **bold** / 200 m chord (8-pt chord count); band at the nearest trunk cell");
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let b2 = ValleyConstruction { a_min_km2: 1.0, ..base };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    for (fname, f) in [("S1 (the construction's input)", &s1), ("PRE (Findings 125–126's instrument)", &pre)] {
        eprintln!("\n   on {fname}:");
        for (label, vc) in [
            ("D8 raw", ValleyConstruction { smooth_m: 0.0, ..b2 }),
            ("D8 + 250 m (the existing)", b2),
            ("TRACÉ bicubic", ValleyConstruction { skeleton_trace: Some(SkeletonTrace::f127(TraceInterp::Bicubic)), ..b2 }),
            ("TRACÉ bilinear (neg. control)", ValleyConstruction { skeleton_trace: Some(SkeletonTrace::f127(TraceInterp::Bilinear)), ..b2 }),
        ] {
            let t = Instant::now();
            let sk = skeleton(f, &vc, &ss, DOMAIN_KM);
            eprintln!("      {label:<30} {} · {:.0} s", r8_line(&sk), t.elapsed().as_secs_f64());
            if sk.trace_stats.is_some() {
                eprintln!("      {:<30} {}", "", trace_row(&sk));
            }
        }
    }
    eprintln!("\n==========  end Finding 127-B1/B2 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// Finding 124's σ: the p50 over land (every 37th cell) of the 3×3 standard deviation of heights (m).
fn sigma_p50(g: &GridF32, ss: &SteinSteinParams) -> f32 {
    let (w, h) = (g.width, g.height);
    let m = |k: usize| c1_altitude_norm_to_metres(g.data[k], ss);
    let mut sig = Vec::new();
    for k in (0..w * h).step_by(37) {
        let (x, y) = (k % w, k / w);
        if g.data[k] <= SEA || x == 0 || y == 0 || x + 1 >= w || y + 1 >= h {
            continue;
        }
        let mut v = [0f32; 9];
        let mut i = 0;
        for dy in 0..3 {
            for dx in 0..3 {
                v[i] = m((y + dy - 1) * w + x + dx - 1);
                i += 1;
            }
        }
        let mu = v.iter().sum::<f32>() / 9.0;
        sig.push((v.iter().map(|a| (a - mu) * (a - mu)).sum::<f32>() / 9.0).sqrt());
    }
    sig.sort_by(f32::total_cmp);
    sig[sig.len() / 2]
}

/// ADR Finding 127-A — the foot clause (`WallProfile::foot_quiet`, τ = 0.5), mur ↔ mer ON everywhere.
/// B1 and B2 (A_c) + B1 with and without it; the témoin; B1b with and without it (the control: no
/// foot, so the clause must be a bit-exact no-op). Finding 126's instruments.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f127_a --nocapture
#[test]
#[ignore]
fn f127_a() {
    use ymir_core::tectonics_c1::bench_guard::field_hash;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 127-A . bruit ↔ pied, the foot clause (mur ↔ mer ON everywhere)  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let zpre: Vec<f32> = (0..n).map(|k| c1_altitude_norm_to_metres(pre.data[k], &ss)).collect();
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let prof = WallProfile::f124();
    let quiet = WallProfile { foot_quiet: Some(0.5), ..prof };
    let (foot_m, crest_zone_m) = (prof.foot_m, 3.0 * prof.crest_m);
    let mm = Some(0.5f32);
    let base = ValleyConstruction { wall_sea_floor_m: mm, ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let knobs = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    // ── the two bit-identity controls first ──
    let t = Instant::now();
    let ht = format!("{:016x}", field_hash(&build_field_seed(knobs(base), PSEED)));
    eprintln!(
        "   CONTROL 1 · témoin + mur ↔ mer: field hash {ht} (the definition's reference a8d2d538d692c2f0) → {}",
        if ht == "a8d2d538d692c2f0" { "BIT-IDENTICAL (the clause never binds on a planar wall)" } else { "⚠ DIFFERS" }
    );
    let b1b = WallProfile { foot_m: 0.0, crest_m: 0.0, ..prof };
    let ha = format!("{:016x}", field_hash(&build_field_seed(knobs(ValleyConstruction { wall_profile: Some(b1b), ..base }), PSEED)));
    let hb = format!(
        "{:016x}",
        field_hash(&build_field_seed(knobs(ValleyConstruction { wall_profile: Some(WallProfile { foot_quiet: Some(0.5), ..b1b }), ..base }), PSEED))
    );
    eprintln!(
        "   CONTROL 2 · B1b (no foot) + mur ↔ mer, without / with the foot clause: {ha} / {hb} → {} · {:.0} s",
        if ha == hb { "BIT-IDENTICAL" } else { "⚠ DIFFERS" },
        t.elapsed().as_secs_f64()
    );
    let worlds: [(&str, ValleyConstruction); 5] = [
        ("témoin C2/10 col", base),
        ("B1 + mm", ValleyConstruction { wall_profile: Some(prof), ..base }),
        ("B1 + mm + FOOT CLAUSE", ValleyConstruction { wall_profile: Some(quiet), ..base }),
        ("B2 (A_c) + B1 + mm", ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), ..base }),
        ("B2 (A_c) + B1 + mm + FOOT CLAUSE", ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(quiet), ..base }),
    ];
    let near = 2.0 / CELL_KM;
    let mut g_t: Option<GridF32> = None;
    let mut summary = Vec::new();
    for (label, vc) in worlds {
        let t = Instant::now();
        let g = build_field_seed(knobs(vc), PSEED);
        let secs = t.elapsed().as_secs_f64();
        eprintln!("\n────────── {label} · build {secs:.0} s ──────────");
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let tan = vc.wall_deg.to_radians().tan();
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw));
        let co = coast(&g, &ss, &dwall, near);
        let nn = co.near.iter().filter(|&&b| b).count();
        let (rn, rf) = (nn as f32 / co.l_near_km.max(1e-6), (co.mids.len() - nn) as f32 / co.l_far_km.max(1e-6));
        let sig = sigma_p50(&g, &ss);
        let r8w = r8_terrain(&g, &cls, 16);
        let r8t = r8_terrain(g_t.as_ref().unwrap_or(&g), &cls, 16);
        eprintln!(
            "   COAST: {} spurs · near a coastal wall {nn} over {:.1} km → **{rn:.4} /km** · elsewhere **{rf:.4} /km** · σ p50 **{sig:.3} m** · \
             R8 terrain {r8w:.4} (témoin on the same classes {r8t:.4})",
            co.mids.len(),
            co.l_near_km
        );
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut ce = 0usize;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            let dfill = fill_del[floor] - fill_pre[floor];
            if over_dug_depression(dfill, b.rim) {
                ce += 1;
                eprintln!(
                    "   OVER-DUG body {} · {:.2} km² (domain) · floor ({},{}) · Δfill {dfill:.1} m · class {}",
                    b.id,
                    b.km2,
                    floor % w,
                    floor / w,
                    CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea")
                );
            }
        }
        drop(fill_del);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage
        }));
        let mut row = format!("{label:<34} spurs {:>4} (near {rn:.4}/km) · σ {sig:.3} m · R8 terr {r8w:.4} · over-dug {ce}", co.mids.len());
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
                let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
                let bk = Buckets::new(&sk, 32);
                let pos = |c: usize| bk.position(c, &zpre, tan, sk.cell_m, foot_m, crest_zone_m);
                let (mut mouth, mut head, mut all) = ([0usize; 4], [0usize; 4], [0usize; 4]);
                for &i in &teeth {
                    let o = own(&segs[i]);
                    let c = |p: &(u32, u32)| p.1 as usize * w + p.0 as usize;
                    mouth[pos(c(o.last().expect("own")))] += 1;
                    head[pos(c(&o[0]))] += 1;
                    for p in o {
                        all[pos(c(p))] += 1;
                    }
                }
                let mut area = [0usize; 4];
                for k in (0..n).step_by(13) {
                    if wall[k] {
                        area[pos(k)] += 1;
                    }
                }
                eprintln!("   teeth {} · R8 network (chord 8) {:.4} · Finding 38 holds", teeth.len(), r8_chords(&pts, 8));
                eprintln!("   C · teeth MOUTHS: {}", shares(&mouth));
                eprintln!("   C · teeth HEADS:  {}", shares(&head));
                eprintln!("   C · teeth POINTS: {}", shares(&all));
                eprintln!("   C · the WALL's own area (the base rate): {}", shares(&area));
                row += &format!(" · teeth {} · R8 net {:.4}", teeth.len(), r8_chords(&pts, 8));
            }
            Err(e) => {
                let msg = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
                eprintln!("   ⛔ the HD assembly PANICKED: {}", msg.lines().next().unwrap_or(""));
                row += " · ⛔ Finding 38";
            }
        }
        row += &format!(" · build {secs:.0} s");
        summary.push(row);
        if g_t.is_none() {
            g_t = Some(g);
        }
    }
    eprintln!("\n   ── summary ──");
    for s in summary {
        eprintln!("   {s}");
    }
    eprintln!("\n==========  end Finding 127-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 127-C — the four canyons of B2 → A_c at their TRUE col (Finding 119: the receiver of
/// `Lake::outlet` on the assembly's own flow, not the last water cell), read in the canyons' world
/// and at the same cells in the témoin and in B2 ≥ 1 km² (where no canyon forms). Planar worlds: the
/// mur ↔ mer clause never binds there (a planar wall is ≥ its floor ≥ sea + base_m), so these are
/// Finding 126's worlds bit for bit.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f127_c --nocapture
#[test]
#[ignore]
fn f127_c() {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 127-C . the canyons at their true col  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let base = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let recv = |dir: &[u8], c: usize| -> Option<usize> {
        let d = dir[c];
        if d == DIR_NONE {
            return None;
        }
        let nx = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        Some(ny * w + nx)
    };
    // (world, its eroded field, breached field, lake map, lake levels by id, masks + classes)
    struct Read {
        g: GridF32,
        bre: GridF32,
        lake_map: Vec<u32>,
        levels: std::collections::HashMap<u32, (f32, usize)>,
        carved: Vec<bool>,
        floor: Vec<bool>,
        cls: Vec<u8>,
        dir: Vec<u8>,
        fill_del: Vec<f32>,
    }
    let read = |vc: ValleyConstruction| -> Read {
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        let levels = dr
            .lakes
            .iter()
            .map(|l| (l.base.id, (l.level_m, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize)))
            .collect();
        Read { g, bre, lake_map: dr.lake_map, levels, carved: mk.carved, floor: mk.floor, cls, dir: dr.flow.direction, fill_del }
    };
    let prim = |r: &Read, c: usize| -> String {
        let what = if !r.carved[c] { "UNCARVED (the terrain: a rim)" } else if r.floor[c] { "a FLOOR" } else { "a WALL" };
        format!("{what} of {}", CLS_NAME.get(r.cls[c] as usize).copied().unwrap_or("sea"))
    };
    // ── the canyons' world first: the bodies, their lakes, their outlets, their true cols ──
    let t = Instant::now();
    let ac = read(ValleyConstruction { a_min_km2: 0.1, ..base });
    let dd_ac = c1_drainage_windowed(&ac.g, None, &dcfg(), &ss, DOMAIN_KM);
    let bre_ac = breach_monotone(&ac.g, &dd_ac.flow.filled, &dd_ac.lake_map, SEA, w, h);
    drop(dd_ac);
    set_dump(true);
    let _cr = f95_criteria(&ac.g, &bre_ac, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
    set_dump(false);
    drop(bre_ac);
    let mut canyons: Vec<(u32, f32, usize)> = Vec::new(); // (id, km², floor)
    for b in take_bodies() {
        let floor = *b.cells.iter().min_by(|&&a, &&c| ac.g.data[a].total_cmp(&ac.g.data[c])).expect("body");
        if over_dug_depression(ac.fill_del[floor] - fill_pre[floor], b.rim) {
            canyons.push((b.id, b.km2, floor));
        }
    }
    eprintln!("   B2 → A_c built and read in {:.0} s · {} over-dug bodies", t.elapsed().as_secs_f64(), canyons.len());
    let t = Instant::now();
    let tm = read(base);
    let b1 = read(ValleyConstruction { a_min_km2: 1.0, ..base });
    eprintln!("   témoin and B2 ≥ 1 km² built and read in {:.0} s", t.elapsed().as_secs_f64());
    for (id, km2, floor) in canyons {
        // the lake the assembly gave the body; `Lake::outlet` is its last WATER cell; the true col is
        // that cell's receiver on the assembly's own flow (Finding 119)
        let in_lake = |c: usize| ac.lake_map[c] == id;
        let Some(&(level, outlet)) = ac.levels.get(&id) else {
            eprintln!("\n   body {id}: not in the assembly's lakes (instrument check FAILED)");
            continue;
        };
        let Some(col) = recv(&ac.dir, outlet) else {
            eprintln!("\n   body {id}: its outlet has no receiver (instrument check FAILED)");
            continue;
        };
        let below = recv(&ac.dir, col);
        let s_col = below.map_or(f32::NAN, |r| {
            let diag = (col % w != r % w) && (col / w != r / w);
            (m(&ac.bre, col) - m(&ac.bre, r)) / (CELL_KM * 1000.0 * if diag { std::f32::consts::SQRT_2 } else { 1.0 })
        });
        eprintln!(
            "\n   CANYON body {id} · {km2:.2} km² (domain) · floor ({},{}) · lake level {level:.1} m · outlet (water) ({},{}) · TRUE COL ({},{}) · S at the col {s_col:.4} m/m · instrument: floor in the lake {} · outlet in the lake {} · col outside {}",
            floor % w,
            floor / w,
            outlet % w,
            outlet / w,
            col % w,
            col / w,
            in_lake(floor),
            in_lake(outlet),
            !in_lake(col)
        );
        eprintln!(
            "      PRE (pre-incision): floor {:.1} m · col {:.1} m · PRE's own fill at the floor {:.1} m (depth below PRE's sill) · at the col {:.1} m",
            m(&pre, floor),
            m(&pre, col),
            fill_pre[floor],
            fill_pre[col]
        );
        for (wname, r) in [("témoin C2/10 col", &tm), ("B2 ≥ 1 km²", &b1), ("B2 → A_c (the canyon)", &ac)] {
            let lid = r.lake_map[floor];
            eprintln!(
                "      {wname:<22} floor {:.1} m (breached {:.1}) · col {:.1} m (breached {:.1}) · fill above PRE at the floor Δ {:+.1} m · lake at the floor {} · the col is {} · the floor is {}",
                m(&r.g, floor),
                m(&r.bre, floor),
                m(&r.g, col),
                m(&r.bre, col),
                r.fill_del[floor] - fill_pre[floor],
                if lid == 0 { "none".to_string() } else { format!("id {lid}, level {:.1} m", r.levels.get(&lid).map_or(f32::NAN, |v| v.0)) },
                prim(r, col),
                prim(r, floor)
            );
        }
    }
    eprintln!("\n==========  end Finding 127-C . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

// ════════════════════════════════ Finding 128 ════════════════════════════════
//
// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f128_c0 --nocapture

/// ADR Finding 128-C0 — the implementation checked against the source BEFORE any calibration: Orlandini
/// et al. 2003's planar slope (relative drainage-area error εA, their Fig. 3a–c) and Paik 2008's plane
/// (cumulative lateral deviation, his Table 1), with the geometry declared in `reading_path_based.md`
/// §7: a 1:4 lateral:longitudinal plane (the TDD 14.04° off the cardinal axis), 30² and 60², sink cells
/// around it. Theoretical areas: the analytic strip of each cell's drainage lines, integrated on 8 × 8
/// sub-points per cell.
#[test]
#[ignore]
fn f128_c0() {
    use ymir_core::tectonics_c1::valley_construction::{accumulate_cells, ltd_directions};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    eprintln!("\n==========  Finding 128-C0 . D8-LTD against Orlandini 2003 and Paik 2008 (the planes)  ==========");
    let th = 0.25f32.atan();
    let t = (th.sin(), -th.cos()); // downslope: toward +x and north (−y)
    // DIAGNOSTIC (added after the first run failed; the geometry is unchanged): the same plane without the
    // 1 000 m offset, i.e. a smaller f32 rounding in the facet slopes. 1:4 makes EXACT ties in the
    // cumulative deviation (N N N NE, a tie at the second step), which rounding then breaks cell by cell.
    for (nin, offset) in [(30usize, 1000.0f32), (60, 1000.0), (30, 0.0), (60, 0.0)] {
        let n = nin + 2;
        let zm: Vec<f32> = (0..n * n)
            .map(|k| {
                let (x, y) = ((k % n) as f32, (k / n) as f32);
                offset + (y * th.cos() - x * th.sin()) * 10.0
            })
            .collect();
        let sink: Vec<bool> = (0..n * n).map(|k| {
            let (x, y) = (k % n, k / n);
            x == 0 || y == 0 || x == n - 1 || y == n - 1
        }).collect();
        let step = |c: usize, k: u8| ((c / n) as i32 + D8_DY[k as usize]) as usize * n + ((c % n) as i32 + D8_DX[k as usize]) as usize;
        let d8: Vec<u8> = (0..n * n)
            .map(|c| {
                if sink[c] {
                    return DIR_NONE;
                }
                let mut best = (DIR_NONE, 0f32);
                for k in 0..8u8 {
                    let m = step(c, k);
                    let s = (zm[c] - zm[m]) / if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                    if s > best.1 {
                        best = (k, s);
                    }
                }
                best.0
            })
            .collect();
        let (ltd, flats) = ltd_directions(&zm, &sink, &d8, n, n, 10.0);
        // theoretical areas (cells): every sub-point's exact drainage line, the cells it crosses
        let sub = 8usize;
        let mut at = vec![0f64; n * n];
        for c in (0..n * n).filter(|&c| !sink[c]) {
            for sy in 0..sub {
                for sx in 0..sub {
                    let (mut px, mut py) = (
                        (c % n) as f32 + (sx as f32 + 0.5) / sub as f32,
                        (c / n) as f32 + (sy as f32 + 0.5) / sub as f32,
                    );
                    let mut last = usize::MAX;
                    loop {
                        if px < 1.0 || py < 1.0 || px >= (n - 1) as f32 || py >= (n - 1) as f32 {
                            break;
                        }
                        let cc = py as usize * n + px as usize;
                        if cc != last {
                            at[cc] += 1.0 / (sub * sub) as f64;
                            last = cc;
                        }
                        px += 0.02 * t.0;
                        py += 0.02 * t.1;
                    }
                }
            }
        }
        let errs = |dir: &[u8]| -> (f64, f64, f64) {
            let acc = accumulate_cells(dir, &sink, n, n);
            let e: Vec<f64> = (0..n * n).filter(|&c| !sink[c]).map(|c| (acc[c] as f64 - at[c]) / at[c]).collect();
            let k = e.len() as f64;
            (e.iter().sum::<f64>() / k, e.iter().map(|v| v.abs()).sum::<f64>() / k, (e.iter().map(|v| v * v).sum::<f64>() / k).sqrt())
        };
        // Paik's cumulative lateral deviation: every downstream cell's distance to the start cell's exact line
        let lateral = |dir: &[u8]| -> f64 {
            let mut tot = 0f64;
            for c0 in (0..n * n).filter(|&c| !sink[c]) {
                let (x0, y0) = ((c0 % n) as f32 + 0.5, (c0 / n) as f32 + 0.5);
                let mut c = c0;
                while !sink[c] {
                    c = step(c, dir[c]);
                    let (x, y) = ((c % n) as f32 + 0.5, (c / n) as f32 + 0.5);
                    tot += ((x - x0) * t.1 - (y - y0) * t.0).abs() as f64;
                }
            }
            tot
        };
        let (e8, el) = (errs(&d8), errs(&ltd));
        let (l8, ll) = (lateral(&d8), lateral(&ltd));
        let (rmae, rrmse, rlat) = (el.1 / e8.1, el.2 / e8.2, ll / l8);
        let confl = |dir: &[u8]| {
            let mut indeg = vec![0u8; n * n];
            for c in (0..n * n).filter(|&c| !sink[c]) {
                let r = step(c, dir[c]);
                if !sink[r] {
                    indeg[r] = indeg[r].saturating_add(1);
                }
            }
            indeg.iter().filter(|&&d| d >= 2).count()
        };
        eprintln!(
            "\n   plane {nin}² (1:4, TDD 14.04° off north), z offset {offset} m · LTD flat fallbacks {flats} · confluences (cells with ≥ 2 donors): D8 {} · LTD {}",
            confl(&d8),
            confl(&ltd)
        );
        eprintln!("      D8    : ME {:+.3} · MAE {:.3} · RMSE {:.3} · Paik's cumulative lateral deviation {:.0} cells", e8.0, e8.1, e8.2, l8);
        eprintln!("      D8-LTD: ME {:+.3} · MAE {:.3} · RMSE {:.3} · Paik's cumulative lateral deviation {:.0} cells", el.0, el.1, el.2, ll);
        eprintln!(
            "      ratios LTD / D8: MAE **{rmae:.3}** (paper ≈ 0.44, pass [0.2, 0.9]) · RMSE **{rrmse:.3}** (≈ 0.13, [0.05, 0.4]) · lateral **{rlat:.3}** (Paik 0.04, [0.01, 0.15]) → {}",
            if (0.2..=0.9).contains(&rmae) && (0.05..=0.4).contains(&rrmse) && (0.01..=0.15).contains(&rlat) { "**PASS**" } else { "**FAIL**" }
        );
    }
    eprintln!("\n==========  end Finding 128-C0  ==========\n");
}

/// ADR Finding 128-A — the confluence clause (`ValleyConstruction::trunk_band`). mur ↔ mer and bruit ↔
/// pied are ON everywhere. Worlds: the témoin; B2 → A_c without and with the clause; B2 at 1 km² with
/// it; B2 (A_c) + B1 + the foot clause with it. Read: the canyons and their true cols (7 and 14 at the
/// cells Finding 127 found); Finding 121's trap (the trunk floors' slope and concavity against the
/// témoin, the source cones); canyon 13's dam (the témoin's drainage path walked under the clause);
/// canyon 3's instrument, repaired (the lowest exit of its lake on the assembly's flow).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f128_a --nocapture
#[test]
#[ignore]
fn f128_a() {
    use std::collections::HashMap;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 128-A . the confluence clause (mur ↔ mer and bruit ↔ pied ON)  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let prof = WallProfile { foot_quiet: Some(0.5), ..WallProfile::f124() };
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let at = |x: usize, y: usize| y * w + x;
    let fixed: [(&str, usize); 4] = [
        ("canyon 7's floor (4588,3240)", at(4588, 3240)),
        ("canyon 7's col (4551,3280)", at(4551, 3280)),
        ("canyon 14's floor (1835,4562)", at(1835, 4562)),
        ("canyon 14's col (1805,4580)", at(1805, 4580)),
    ];
    let recv = |dir: &[u8], c: usize| -> Option<usize> {
        let d = dir[c];
        if d == DIR_NONE {
            return None;
        }
        let nx = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        Some(ny * w + nx)
    };
    struct Read {
        g: GridF32,
        bre: GridF32,
        lake_map: Vec<u32>,
        levels: HashMap<u32, (f32, usize)>,
        dir: Vec<u8>,
        canyons: Vec<(u32, f32, usize)>,
    }
    let worlds: [(&str, ValleyConstruction); 5] = [
        ("témoin C2/10 col", base),
        ("B2 → A_c (no clause)", ac),
        ("B2 → A_c + CONFLUENCE CLAUSE", ValleyConstruction { trunk_band: true, ..ac }),
        ("B2 ≥ 1 km² + CONFLUENCE CLAUSE", ValleyConstruction { a_min_km2: 1.0, trunk_band: true, ..base }),
        (
            "B2 (A_c) + B1 + foot + CONFLUENCE CLAUSE",
            ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), trunk_band: true, ..base },
        ),
    ];
    let near = 2.0 / CELL_KM;
    let mut keep: HashMap<&str, Read> = HashMap::new();
    let mut g_t: Option<GridF32> = None;
    let mut summary = Vec::new();
    for (label, vc) in worlds {
        let t = Instant::now();
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let secs = t.elapsed().as_secs_f64();
        eprintln!("\n────────── {label} · build {secs:.0} s ──────────");
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw));
        let co = coast(&g, &ss, &dwall, near);
        let nn = co.near.iter().filter(|&&b| b).count();
        let rn = nn as f32 / co.l_near_km.max(1e-6);
        let sig = sigma_p50(&g, &ss);
        let r8w = r8_terrain(&g, &cls, 16);
        let r8t = r8_terrain(g_t.as_ref().unwrap_or(&g), &cls, 16);
        eprintln!(
            "   COAST {} spurs · near a coastal wall **{rn:.4} /km** · σ p50 **{sig:.3} m** · R8 terrain {r8w:.4} (témoin on the same classes {r8t:.4})",
            co.mids.len()
        );
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = Vec::new();
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            let dfill = fill_del[floor] - fill_pre[floor];
            if over_dug_depression(dfill, b.rim) {
                eprintln!(
                    "   OVER-DUG body {} · {:.2} km² (domain) · floor ({},{}) · Δfill {dfill:.1} m · class {}",
                    b.id,
                    b.km2,
                    floor % w,
                    floor / w,
                    CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea")
                );
                canyons.push((b.id, b.km2, floor));
            }
        }
        eprintln!("   over-dug bodies: {}", canyons.len());
        drop(fill_del);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage
        }));
        let mut row = format!(
            "{label:<42} over-dug {} · spurs near {rn:.4}/km · σ {sig:.3} m · R8 terr {r8w:.4}",
            canyons.len()
        );
        let heights: Vec<String> = fixed.iter().map(|(nm, k)| format!("{nm} {:.1} m", mm(&g, *k))).collect();
        eprintln!("   heights: {}", heights.join(" · "));
        match res {
            Ok(dr) => {
                let segs = &dr.rivers.segments;
                let wall: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k]).collect();
                let wc: Vec<usize> = (0..segs.len())
                    .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2)
                    .collect();
                let teeth = wc
                    .iter()
                    .filter(|&&i| {
                        let p = &segs[i].points;
                        p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32 >= 0.8 * p.len() as f32
                    })
                    .count();
                let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
                eprintln!("   teeth {teeth} · R8 network (chord 8) {:.4} · Finding 38 holds", r8_chords(&pts, 8));
                row += &format!(" · teeth {teeth} · R8 net {:.4}", r8_chords(&pts, 8));
                if label.starts_with("témoin") || label.starts_with("B2 → A_c") {
                    let levels = dr
                        .lakes
                        .iter()
                        .map(|l| (l.base.id, (l.level_m, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize)))
                        .collect();
                    keep.insert(label, Read { g: g.clone(), bre: bre.clone(), lake_map: dr.lake_map.clone(), levels, dir: dr.flow.direction.clone(), canyons: canyons.clone() });
                }
            }
            Err(e) => {
                let msg = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
                eprintln!("   ⛔ the HD assembly PANICKED: {}", msg.lines().next().unwrap_or(""));
                row += " · ⛔ Finding 38";
            }
        }
        row += &format!(" · build {secs:.0} s");
        summary.push(row);
        if g_t.is_none() {
            g_t = Some(g);
        }
    }
    // ── Finding 121's trap: the trunk floors under the clause against the témoin (same trunk cells) ──
    eprintln!("\n   ── Finding 121's trap ──");
    let sk_t = skeleton(&pre, &base, &ss, DOMAIN_KM);
    let (tm, a0, a1) = (&keep["témoin C2/10 col"], &keep["B2 → A_c (no clause)"], &keep["B2 → A_c + CONFLUENCE CLAUSE"]);
    let trunk_cells: Vec<usize> = (0..n).filter(|&k| sk_t.trunk[k] && sk_t.area_km2[k] >= 10.0 && pre.data[k] > SEA).collect();
    for (label, r) in [("B2 → A_c (no clause)", a0), ("B2 → A_c + CLAUSE", a1)] {
        let (mut dz, mut ratio) = (Vec::new(), Vec::new());
        let (mut lt, mut lw) = (Vec::new(), Vec::new());
        for &k in &trunk_cells {
            let Some(rc) = recv(&sk_t.direction, k) else { continue };
            if !sk_t.trunk[rc] {
                continue;
            }
            let dist = CELL_KM * 1000.0 * if (k % w != rc % w) && (k / w != rc / w) { std::f32::consts::SQRT_2 } else { 1.0 };
            let (st, sw) = ((mm(&tm.g, k) - mm(&tm.g, rc)) / dist, (mm(&r.g, k) - mm(&r.g, rc)) / dist);
            dz.push((mm(&r.g, k) - mm(&tm.g, k)).abs());
            if st > 1e-4 && sw > 1e-4 {
                ratio.push(sw / st);
                let la = sk_t.area_km2[k].ln();
                lt.push((la, st.ln()));
                lw.push((la, sw.ln()));
            }
        }
        dz.sort_by(f32::total_cmp);
        ratio.sort_by(f32::total_cmp);
        let q = |v: &[f32], p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p) as usize] };
        let theta = |v: &[(f32, f32)]| -> f32 {
            let k = v.len() as f32;
            let (mx, my) = (v.iter().map(|a| a.0).sum::<f32>() / k, v.iter().map(|a| a.1).sum::<f32>() / k);
            let (sxy, sxx) = (v.iter().map(|a| (a.0 - mx) * (a.1 - my)).sum::<f32>(), v.iter().map(|a| (a.0 - mx).powi(2)).sum::<f32>());
            -sxy / sxx
        };
        eprintln!(
            "   {label:<22} on {} trunk ≥ 10 km² cells: |Δz| vs the témoin p50 {:.2} / p90 {:.2} m · slope ratio (world / témoin) p50 {:.3} · concavity θ (S ∝ A^−θ) témoin {:.3} → world {:.3}",
            trunk_cells.len(),
            q(&dz, 0.5),
            q(&dz, 0.9),
            q(&ratio, 0.5),
            theta(&lt),
            theta(&lw)
        );
    }
    // the source cones: what the clause moves, and how much of it lies near a line's head
    let sk_ac = skeleton(&pre, &ac, &ss, DOMAIN_KM);
    let (c_free, _) = carve(&pre, &sk_ac, &ac, &ss);
    let (c_band, mkb) = carve(&pre, &sk_ac, &ValleyConstruction { trunk_band: true, ..ac }, &ss);
    let moved: Vec<usize> = (0..n).filter(|&k| c_free.data[k] != c_band.data[k]).collect();
    let mut head = vec![false; n];
    for l in &sk_ac.polylines {
        let (x, y) = ((l[0].0.floor() as i64).rem_euclid(w as i64) as usize, (l[0].1.floor() as i64).rem_euclid(h as i64) as usize);
        head[y * w + x] = true;
    }
    let dhead = dist_from(&head, w, h);
    let near_head = moved.iter().filter(|&&k| dhead[k] <= 3).count();
    let raised = moved.iter().filter(|&&k| c_band.data[k] > c_free.data[k]).count();
    eprintln!(
        "   the clause moves {} cells at the construction (on PRE, B2 → A_c): {} raised (a dam removed / a trunk floor laid higher than a smaller line's cone), {} lowered · {} of them within 3 cells of a line's head (the source cones) · all on a floor band: {}",
        moved.len(),
        raised,
        moved.len() - raised,
        near_head,
        moved.iter().all(|&k| mkb.floor[k])
    );
    drop((c_free, c_band, mkb, head, dhead));
    // ── canyon 13's dam: the témoin's drainage path from its floor, walked under the clause ──
    eprintln!("\n   ── canyon 13 (floor (2281,4563)) and canyon 3 ──");
    let f13 = at(2281, 4563);
    let in_a1 = a1.canyons.iter().any(|&(_, _, fl)| {
        let (dx, dy) = ((fl % w) as i64 - 2281, (fl / w) as i64 - 4563);
        dx * dx + dy * dy <= 400
    });
    eprintln!("   canyon 13 {} under the clause", if in_a1 { "PERSISTS" } else { "is GONE" });
    for (label, r) in [("B2 → A_c (no clause)", a0), ("B2 → A_c + CLAUSE", a1)] {
        let lid = r.lake_map[f13];
        let level = r.levels.get(&lid).map_or(f32::NAN, |v| v.0);
        let mut c = f13;
        let mut dam: Option<(usize, usize, f32)> = None;
        for step in 0..20_000 {
            let Some(nx) = recv(&tm.dir, c) else { break };
            c = nx;
            if tm.g.data[c] <= SEA {
                break;
            }
            if lid != 0 && r.lake_map[c] != lid && mm(&r.bre, c) > level + 0.01 {
                dam = Some((c, step, mm(&r.bre, c)));
                break;
            }
        }
        match dam {
            Some((c, step, z)) => eprintln!(
                "   {label}: lake {lid} at {level:.1} m · the témoin's path is DAMMED at ({},{}), {step} cells down, at {z:.1} m (témoin {:.1} m)",
                c % w,
                c / w,
                mm(&tm.bre, c)
            ),
            None => eprintln!("   {label}: lake {lid} (level {level:.1} m) · no dam on the témoin's path above the lake level"),
        }
    }
    // canyon 3's instrument, repaired: the lowest EXIT of its lake on the assembly's own flow
    for (label, r) in [("B2 → A_c (no clause)", a0), ("B2 → A_c + CLAUSE", a1)] {
        let f3 = at(3892, 2310);
        let lid = r.lake_map[f3];
        if lid == 0 {
            eprintln!("   {label}: canyon 3's floor carries no lake");
            continue;
        }
        let exits: Vec<(usize, usize)> = (0..n)
            .filter(|&k| r.lake_map[k] == lid)
            .filter_map(|k| recv(&r.dir, k).filter(|&rc| r.lake_map[rc] != lid).map(|rc| (k, rc)))
            .collect();
        let best = exits.iter().min_by(|a, b| r.bre.data[a.1].total_cmp(&r.bre.data[b.1]));
        let recorded = r.levels.get(&lid).map(|v| v.1);
        match best {
            Some(&(_, col)) => eprintln!(
                "   {label}: lake {lid} · {} exits on the assembly's flow · the lowest exit's receiver = the col ({},{}) at {:.1} m (témoin {:.1} m) · recorded Lake::outlet {} (in the lake: {})",
                exits.len(),
                col % w,
                col / w,
                mm(&r.bre, col),
                mm(&tm.bre, col),
                recorded.map_or("none".to_string(), |k| format!("({},{})", k % w, k / w)),
                recorded.is_some_and(|k| r.lake_map[k] == lid)
            ),
            None => eprintln!("   {label}: lake {lid} has no exit on the assembly's flow (a closed lake) · recorded outlet in the lake: {}", recorded.is_some_and(|k| r.lake_map[k] == lid)),
        }
    }
    eprintln!("\n   ── summary ──");
    for s in summary {
        eprintln!("   {s}");
    }
    eprintln!("\n==========  end Finding 128-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}


// ════════════════════════════════ Finding 129 ════════════════════════════════

/// ADR Finding 129-B1 — the reserved case: [O03]'s PARABOLIC VALLEY (Fig. 2d, Fig. 3j–l), which was not
/// used to find the tie bug. The geometry is declared in `f129_predictions.md` before the run: 30² interior
/// cells, the flow along +x to the east edge, z = (n − x) + (y − yc)²/30 (cell units, metres per cell),
/// sink cells around it. Theoretical areas come from the exact descent lines of 8 × 8 sub-points per cell.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f129_b1 --nocapture
#[test]
#[ignore]
fn f129_b1() {
    use ymir_core::tectonics_c1::valley_construction::{accumulate_cells, ltd_directions};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    eprintln!("\n==========  Finding 129-B1 . D8-LTD on the reserved case, Orlandini 2003's parabolic valley  ==========");
    let n = 32usize;
    let yc = n as f32 / 2.0;
    let zf = |x: f32, y: f32| (n as f32 - x) + (y - yc) * (y - yc) / 30.0;
    let zm: Vec<f32> = (0..n * n).map(|k| zf((k % n) as f32 + 0.5, (k / n) as f32 + 0.5) * 10.0).collect();
    let sink: Vec<bool> = (0..n * n).map(|k| {
        let (x, y) = (k % n, k / n);
        x == 0 || y == 0 || x == n - 1 || y == n - 1
    }).collect();
    let step = |c: usize, k: u8| ((c / n) as i32 + D8_DY[k as usize]) as usize * n + ((c % n) as i32 + D8_DX[k as usize]) as usize;
    let d8: Vec<u8> = (0..n * n)
        .map(|c| {
            if sink[c] {
                return DIR_NONE;
            }
            let mut best = (DIR_NONE, 0f32);
            for k in 0..8u8 {
                let s = (zm[c] - zm[step(c, k)]) / if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                if s > best.1 {
                    best = (k, s);
                }
            }
            best.0
        })
        .collect();
    let (ltd, flats) = ltd_directions(&zm, &sink, &d8, n, n, 10.0);
    // the exact descent direction: −∇z = (1, −2 (y − yc) / 30), in cell units
    let sub = 8usize;
    let mut at = vec![0f64; n * n];
    for c in (0..n * n).filter(|&c| !sink[c]) {
        for sy in 0..sub {
            for sx in 0..sub {
                let (mut px, mut py) = ((c % n) as f32 + (sx as f32 + 0.5) / sub as f32, (c / n) as f32 + (sy as f32 + 0.5) / sub as f32);
                let mut last = usize::MAX;
                loop {
                    if px < 1.0 || py < 1.0 || px >= (n - 1) as f32 || py >= (n - 1) as f32 {
                        break;
                    }
                    let cc = py as usize * n + px as usize;
                    if cc != last {
                        at[cc] += 1.0 / (sub * sub) as f64;
                        last = cc;
                    }
                    let (gx, gy) = (1.0f32, -2.0 * (py - yc) / 30.0);
                    let g = (gx * gx + gy * gy).sqrt();
                    px += 0.02 * gx / g;
                    py += 0.02 * gy / g;
                }
            }
        }
    }
    let errs = |dir: &[u8]| -> (f64, f64, f64) {
        let acc = accumulate_cells(dir, &sink, n, n);
        let e: Vec<f64> = (0..n * n).filter(|&c| !sink[c]).map(|c| (acc[c] as f64 - at[c]) / at[c]).collect();
        let k = e.len() as f64;
        (e.iter().sum::<f64>() / k, e.iter().map(|v| v.abs()).sum::<f64>() / k, (e.iter().map(|v| v * v).sum::<f64>() / k).sqrt())
    };
    let (e8, el) = (errs(&d8), errs(&ltd));
    let (rmae, rrmse) = (el.1 / e8.1, el.2 / e8.2);
    let pass = (0.19..=0.39).contains(&el.1) && (0.20..=0.60).contains(&el.2) && (0.3..=0.8).contains(&rmae) && (0.08..=0.40).contains(&rrmse);
    eprintln!("   parabolic valley 30² · LTD flat fallbacks {flats}");
    eprintln!("      D8    : ME {:+.3} · MAE {:.3} · RMSE {:.3}   ([O03] Fig. 3j–l at λ 0: ME ≈ 0.0, MAE ≈ 0.57, RMSE ≈ 2.1)", e8.0, e8.1, e8.2);
    eprintln!("      D8-LTD: ME {:+.3} · MAE **{:.3}** · RMSE **{:.3}**   ([O03] at λ 1: ME ≈ −0.26, MAE ≈ 0.29, RMSE ≈ 0.35)", el.0, el.1, el.2);
    eprintln!("      ratios LTD / D8: MAE **{rmae:.3}** (paper ≈ 0.51) · RMSE **{rrmse:.3}** (≈ 0.17) → **{}**", if pass { "PASS" } else { "FAIL" });
}

/// ADR Finding 129-B2 — the ties on the REAL pre-incision (S1 breached, the surface the skeleton reads):
/// the D8 ties (the two steepest of the eight neighbours' slopes equal, in f32 and after the u16
/// export's quantisation), and the D8-LTD choices the tie rule decides. Where: by local slope quartile.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f129_b2 --nocapture
#[test]
#[ignore]
fn f129_b2() {
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::ltd_directions_stats;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, FlowConfig, compute_flow};
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 129-B2 . the ties on the real pre-incision  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let dc = C1DrainageConfig::default();
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    drop(d);
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let zm: Vec<f32> = bf.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
    let sink: Vec<bool> = bf.data.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let land = sink.iter().filter(|&&s| !s).count();
    let cell_m = CELL_KM * 1000.0;
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    // D8 ties and the local slope (the steepest), on a height field
    let d8_ties = |z: &[f32]| -> (usize, usize, Vec<f32>) {
        let mut exact = 0usize;
        let mut rel = 0usize;
        let mut smax = vec![0f32; n];
        for c in (0..n).filter(|&c| !sink[c]) {
            let mut s: Vec<f32> = (0..8).map(|k| (z[c] - z[nb(c, k)]) / if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 }).collect();
            s.sort_by(|a, b| b.total_cmp(a));
            smax[c] = s[0] / cell_m;
            if s[0] > 0.0 {
                if s[0] == s[1] {
                    exact += 1;
                }
                if (s[0] - s[1]).abs() <= 1e-6 * s[0].abs() {
                    rel += 1;
                }
            }
        }
        (exact, rel, smax)
    };
    let (ex, rel, smax) = d8_ties(&zm);
    let hl = ymir_core::export::height::metric_height_u16(&bf, &ss);
    let zq: Vec<f32> = hl.codes.iter().map(|&c| hl.min_m + (c as f32 / 65535.0) * (hl.max_m - hl.min_m)).collect();
    let (exq, relq, _) = d8_ties(&zq);
    eprintln!(
        "   D8 on {land} land cells: exact ties (f32) {ex} ({:.4} %) · within 1e-6 relative {rel} ({:.4} %) · after the u16 quantisation ({:.3} m / code): exact {exq} ({:.3} %) · within 1e-6 {relq} ({:.3} %)",
        100.0 * ex as f32 / land as f32,
        100.0 * rel as f32 / land as f32,
        (hl.max_m - hl.min_m) / 65535.0,
        100.0 * exq as f32 / land as f32,
        100.0 * relq as f32 / land as f32
    );
    let (_, st) = ltd_directions_stats(&zm, &sink, &flow.direction, w, h, cell_m);
    eprintln!(
        "   D8-LTD: {} land cells · flat fallbacks {} ({:.3} %) · the tie rule decides {} choices ({:.3} %) · of them exact in f32 {}",
        land,
        st.flats,
        100.0 * st.flats as f32 / land as f32,
        st.ties,
        100.0 * st.ties as f32 / land as f32,
        st.exact_ties
    );
    // where: the local slope quartiles of land
    let mut sl: Vec<f32> = (0..n).filter(|&c| !sink[c]).map(|c| smax[c]).collect();
    sl.sort_by(f32::total_cmp);
    let q = |p: f32| sl[((sl.len() - 1) as f32 * p) as usize];
    let (q1, q2, q3) = (q(0.25), q(0.5), q(0.75));
    let mut bins = [0usize; 4];
    for &c in &st.tie_cells {
        let s = smax[c as usize];
        bins[if s <= q1 { 0 } else if s <= q2 { 1 } else if s <= q3 { 2 } else { 3 }] += 1;
    }
    let t = st.tie_cells.len().max(1) as f32;
    eprintln!(
        "   where the LTD ties sit, by local slope quartile of land (Q1 ≤ {q1:.4} · Q2 ≤ {q2:.4} · Q3 ≤ {q3:.4} m/m): {:.1} % / {:.1} % / {:.1} % / {:.1} %",
        100.0 * bins[0] as f32 / t,
        100.0 * bins[1] as f32 / t,
        100.0 * bins[2] as f32 / t,
        100.0 * bins[3] as f32 / t
    );
}

/// Slope-area θ: −slope of ln S on ln A over the skeleton's trunk ≥ `a_min` km² cells, S = the field's
/// drop to the cell's D8 receiver on that skeleton over the link length, cells with S > `smin`.
fn theta_sa(z: &[f32], sk: &Skeleton, a_min: f32, smin: f32) -> (f32, usize) {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let (w, h) = (sk.width, sk.height);
    let mut v: Vec<(f64, f64)> = Vec::new();
    for k in 0..w * h {
        if !sk.trunk[k] || sk.area_km2[k] < a_min {
            continue;
        }
        let d = sk.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let r = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize * w
            + ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        if !sk.trunk[r] {
            continue;
        }
        let dist = sk.cell_m * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
        let s = (z[k] - z[r]) / dist;
        if s > smin {
            v.push(((sk.area_km2[k] as f64).ln(), (s as f64).ln()));
        }
    }
    let k = v.len() as f64;
    let (mx, my) = (v.iter().map(|a| a.0).sum::<f64>() / k, v.iter().map(|a| a.1).sum::<f64>() / k);
    let sxy: f64 = v.iter().map(|a| (a.0 - mx) * (a.1 - my)).sum();
    let sxx: f64 = v.iter().map(|a| (a.0 - mx).powi(2)).sum();
    ((-sxy / sxx) as f32, v.len())
}

/// [`theta_sa`] on two fields over the SAME cells: those with S > `smin` in both (Finding 128's form).
fn theta_sa_paired(za: &[f32], zb: &[f32], sk: &Skeleton, a_min: f32, smin: f32) -> (f32, f32) {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let (w, h) = (sk.width, sk.height);
    let (mut va, mut vb): (Vec<(f64, f64)>, Vec<(f64, f64)>) = (Vec::new(), Vec::new());
    for k in 0..w * h {
        if !sk.trunk[k] || sk.area_km2[k] < a_min || sk.direction[k] == DIR_NONE {
            continue;
        }
        let d = sk.direction[k] as usize;
        let r = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize * w + ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
        if !sk.trunk[r] {
            continue;
        }
        let dist = sk.cell_m * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
        let (sa, sb) = ((za[k] - za[r]) / dist, (zb[k] - zb[r]) / dist);
        if sa > smin && sb > smin {
            let la = (sk.area_km2[k] as f64).ln();
            va.push((la, (sa as f64).ln()));
            vb.push((la, (sb as f64).ln()));
        }
    }
    let fit = |v: &[(f64, f64)]| -> f32 {
        let k = v.len() as f64;
        let (mx, my) = (v.iter().map(|a| a.0).sum::<f64>() / k, v.iter().map(|a| a.1).sum::<f64>() / k);
        let sxy: f64 = v.iter().map(|a| (a.0 - mx) * (a.1 - my)).sum();
        let sxx: f64 = v.iter().map(|a| (a.0 - mx).powi(2)).sum();
        (-sxy / sxx) as f32
    };
    (fit(&va), fit(&vb))
}

/// ADR Finding 129-A — θ by STAGE, which decides whether the confluence clause corrupts or restores the
/// construction's law (θ = 0.5 on its own D8 path), and what the clause does to the source cones.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f129_a --nocapture
#[test]
#[ignore]
fn f129_a() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 129-A . θ by stage, and the source cones  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let breach = |g: &GridF32| -> GridF32 {
        let d = c1_drainage_windowed(g, None, &dcfg(), &ss, DOMAIN_KM);
        breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, g.width, g.height)
    };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let sk = skeleton(&s1, &base, &ss, DOMAIN_KM); // trunk ≥ 10 km² cells are the same D8 cells at every a_min
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let sk_pre = skeleton(&pre, &base, &ss, DOMAIN_KM);
    let th = |z: &[f32]| theta_sa(z, &sk, 10.0, 1e-4);
    // ── the references first ──
    eprintln!("\n   ── references (the same instrument: slope-area over the S1 skeleton's trunk ≥ 10 km² cells, S > 1e-4) ──");
    eprintln!("      the construction's own law: θ = 0.500 (z = base + k·χ, χ = ∫ (A0/A)^0.5 dx)");
    eprintln!("      Harel et al. 2016: ⟨m/n⟩ = 0.51 ± 0.12 (PDF p. 23; 0.51 ± 0.14 at p. 36), by the χ-integral method");
    let (t, k) = th(&metres(&s1));
    eprintln!("      S1, the construction's input: θ **{t:.3}** ({k} cells)");
    let (t, k) = th(&metres(&pre));
    eprintln!("      PRE (no incision, bathymetry; production ships no droplet pass, upscale.rs:410): θ **{t:.3}** ({k} cells)");
    for (label, kn) in [("the delivered (livré)", Knobs::passes(2)), ("A1+B2", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) })] {
        let g = build_field_seed(kn, PSEED);
        let (t, k) = th(&metres(&g));
        eprintln!("      {label}: θ **{t:.3}** ({k} cells)");
    }
    // ── the stages ──
    eprintln!("\n   ── the stages (mur ↔ mer ON; planar walls, so bruit ↔ pied has nothing to act on) ──");
    let worlds: [(&str, ValleyConstruction); 3] = [
        ("témoin C2/10 col", base),
        ("B2 → A_c, no clause", ValleyConstruction { a_min_km2: 0.1, ..base }),
        ("B2 → A_c + CLAUSE", ValleyConstruction { a_min_km2: 0.1, trunk_band: true, ..base }),
    ];
    let mut s2s: Vec<GridF32> = Vec::new();
    let mut finals: Vec<GridF32> = Vec::new();
    for (label, vc) in worlds {
        let t = Instant::now();
        let s2 = build_field_seed(Knobs { valley: Some(vc), no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
        let (t1, k1) = th(&metres(&s2));
        let (t2, k2) = th(&metres(&breach(&s2)));
        let s3 = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
        let (t3, k3) = th(&metres(&s3));
        drop(s3);
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let (t4, k4) = th(&metres(&g));
        eprintln!(
            "   {label:<22} (i) as built **{t1:.3}** ({k1}) · (ii) + drainage/breach **{t2:.3}** ({k2}) · (iii-a) + the light pass **{t3:.3}** ({k3}) · (iii-b) the world **{t4:.3}** ({k4}) · {:.0} s",
            t.elapsed().as_secs_f64()
        );
        s2s.push(s2);
        finals.push(g);
    }
    // Finding 128's own paired form, for continuity: the PRE skeleton, cells with S > 1e-4 in BOTH worlds
    let (zt, za0, za1) = (metres(&finals[0]), metres(&finals[1]), metres(&finals[2]));
    let (p0t, p0w) = theta_sa_paired(&zt, &za0, &sk_pre, 10.0, 1e-4);
    let (p1t, p1w) = theta_sa_paired(&zt, &za1, &sk_pre, 10.0, 1e-4);
    eprintln!(
        "   Finding 128's paired form (PRE skeleton, S > 1e-4 in both): témoin {p0t:.3} → no clause {p0w:.3} · témoin {p1t:.3} → clause {p1w:.3} (F128 read 0.346 → 0.366 and 0.350 → 0.402)"
    );
    drop((zt, za0, za1));
    // ── the source cones: what the clause moves at the construction (the pipeline's own S2) ──
    eprintln!("\n   ── the source cones ──");
    let (z0, z1) = (metres(&s2s[1]), metres(&s2s[2]));
    let sk_ac = skeleton(&s1, &ValleyConstruction { a_min_km2: 0.1, ..base }, &ss, DOMAIN_KM);
    let (_, mk1) = carve(&s1, &sk_ac, &ValleyConstruction { a_min_km2: 0.1, trunk_band: true, ..base }, &ss);
    let mut head = vec![false; n];
    let mut heads: Vec<usize> = Vec::new();
    for (li, l) in sk_ac.polylines.iter().enumerate() {
        let (x, y) = ((l[0].0.floor() as i64).rem_euclid(w as i64) as usize, (l[0].1.floor() as i64).rem_euclid(h as i64) as usize);
        head[y * w + x] = true;
        heads.push(li);
    }
    let dhead = dist_from(&head, w, h);
    let moved: Vec<usize> = (0..n).filter(|&k| s2s[1].data[k] != s2s[2].data[k]).collect();
    let near: Vec<usize> = moved.iter().copied().filter(|&k| dhead[k] <= 3).collect();
    let mut dz: Vec<f32> = near.iter().map(|&k| z1[k] - z0[k]).collect();
    dz.sort_by(f32::total_cmp);
    let q = |v: &[f32], p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p) as usize] };
    let raised = dz.iter().filter(|&&d| d > 0.0).count();
    let on_floor = near.iter().filter(|&&k| mk1.floor[k]).count();
    let uncarved = near.iter().filter(|&&k| !mk1.carved[k]).count();
    eprintln!(
        "   the clause moves {} cells (the pipeline's S2, B2 → A_c) · {} within 3 cells of a line's head · of those: RAISED {raised} / lowered {} · Δz p10 {:+.1} / p50 {:+.1} / p90 {:+.1} m · on the larger line's FLOOR band {on_floor} · left UNCARVED {uncarved} · on a wall {}",
        moved.len(),
        near.len(),
        near.len() - raised,
        q(&dz, 0.1),
        q(&dz, 0.5),
        q(&dz, 0.9),
        near.len() - on_floor - uncarved
    );
    // the touched cones' profiles: the first 10 cells of each line whose head region moved
    let mut moved_mask = vec![false; n];
    for &k in &moved {
        moved_mask[k] = true;
    }
    let (mut touched, mut flatter, mut reversed) = (0usize, 0usize, 0usize);
    let (mut s_before, mut s_after) = (Vec::new(), Vec::new());
    for l in &sk_ac.polylines {
        let cells: Vec<usize> = l
            .iter()
            .map(|p| ((p.1.floor() as i64).rem_euclid(h as i64) as usize) * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize)
            .fold(Vec::new(), |mut acc, c| {
                if acc.last() != Some(&c) {
                    acc.push(c);
                }
                acc
            });
        if cells.len() < 11 || !cells[..4].iter().any(|&c| moved_mask[c]) {
            continue;
        }
        touched += 1;
        let (a, b) = (cells[0], cells[10]);
        let len = 10.0 * CELL_KM * 1000.0;
        let (sb, sa) = ((z0[a] - z0[b]) / len, (z1[a] - z1[b]) / len);
        s_before.push(sb);
        s_after.push(sa);
        if sa < 0.5 * sb {
            flatter += 1;
        }
        if sa < 0.0 {
            reversed += 1;
        }
    }
    s_before.sort_by(f32::total_cmp);
    s_after.sort_by(f32::total_cmp);
    eprintln!(
        "   touched source cones (a line whose first 4 cells moved): {touched} of {} lines · the head-to-10th-cell slope p50 {:.4} → {:.4} m/m · FLATTENED below half {flatter} · REVERSED (a pit at the head) {reversed}",
        sk_ac.polylines.len(),
        q(&s_before, 0.5),
        q(&s_after, 0.5)
    );
    eprintln!("\n==========  end Finding 129-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// A continuous steepest-descent trace (bicubic on `zm`), half a cell per step, from `(x, y)` for
/// `len` cells of arclength or until the sea / a flat / the step budget.
fn free_trace(zm: &[f32], w: usize, h: usize, cell_m: f32, sea_m: f32, x: f32, y: f32, len: f32) -> Vec<(f32, f32)> {
    use ymir_core::tectonics_c1::valley_construction::{TraceInterp, interp_grad};
    let mut p = (x, y);
    let mut out = vec![p];
    let mut arc = 0f32;
    let budget = (4.0 * len / 0.5) as usize + 8;
    for _ in 0..budget {
        if arc >= len {
            break;
        }
        let (z, gx, gy) = interp_grad(zm, w, h, cell_m, TraceInterp::Bicubic, p.0.rem_euclid(w as f32), p.1.rem_euclid(h as f32));
        let g = (gx * gx + gy * gy).sqrt();
        if z <= sea_m || g < 1e-3 {
            break;
        }
        p = (p.0 - 0.5 * gx / g, p.1 - 0.5 * gy / g);
        out.push(p);
        arc += 0.5;
    }
    out
}

/// R8 of a set of polylines' chords of `k` cells of arclength: (R8, chord count).
fn chords_r8_len(polys: &[Vec<(f32, f32)>], k: f32) -> (f32, usize) {
    let (mut c8, mut s8, mut m) = (0f64, 0f64, 0usize);
    for l in polys {
        let (mut i0, mut arc0, mut arc) = (0usize, 0f32, 0f32);
        for i in 1..l.len() {
            arc += ((l[i].0 - l[i - 1].0).powi(2) + (l[i].1 - l[i - 1].1).powi(2)).sqrt();
            if arc - arc0 >= k - 1e-3 {
                let t = ((l[i].1 - l[i0].1) as f64).atan2((l[i].0 - l[i0].0) as f64);
                c8 += (8.0 * t).cos();
                s8 += (8.0 * t).sin();
                m += 1;
                i0 = i;
                arc0 = arc;
            }
        }
    }
    (if m < 40 { f32::NAN } else { ((c8 * c8 + s8 * s8).sqrt() / m as f64) as f32 }, m)
}

/// The 1–3 km² prefixes of a skeleton's lines, as cell-centre chains (the raw D8 or LTD path cells).
fn branch_paths(sk: &Skeleton) -> Vec<Vec<(f32, f32)>> {
    let (w, h) = (sk.width, sk.height);
    sk.polylines
        .iter()
        .map(|l| {
            l.iter()
                .filter(|p| (p.0.fract() - 0.5).abs() < 1e-4 && (p.1.fract() - 0.5).abs() < 1e-4)
                .take_while(|p| {
                    let c = ((p.1.floor() as i64).rem_euclid(h as i64) as usize) * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize;
                    sk.area_km2[c] < 3.0
                })
                .map(|p| (p.0, p.1))
                .collect::<Vec<_>>()
        })
        .filter(|v| v.len() >= 2)
        .collect()
}

/// The concentration (Finding 112) and the cone (Paik 2008) of a pointer tree on a synthetic: one
/// receiver per land cell, the largest drained area (km²), the outlets (all and ≥ 1 km²), the confluences
/// (land cells with ≥ 2 donors), and the sector asymmetry around the dome's centre (coefficient of
/// variation of the drained area leaving through each of 16 azimuth sectors).
fn tree_row(dir: &[u8], area_km2: &[f32], sink: &[bool], w: usize, h: usize) -> String {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let n = w * h;
    let recv = |c: usize| -> Option<usize> {
        let d = dir[c];
        (d != DIR_NONE).then(|| ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize)
    };
    let land: Vec<usize> = (0..n).filter(|&c| !sink[c]).collect();
    let one = land.iter().filter(|&&c| recv(c).is_some()).count();
    let mut indeg = vec![0u8; n];
    let mut sector = [0f64; 16];
    let (cx, cy) = (w as f32 / 2.0, h as f32 / 2.0);
    let (mut outlets, mut big) = (0usize, 0usize);
    for &c in &land {
        if let Some(r) = recv(c) {
            if sink[r] {
                outlets += 1;
                if area_km2[c] >= 1.0 {
                    big += 1;
                }
                let az = ((c / w) as f32 + 0.5 - cy).atan2((c % w) as f32 + 0.5 - cx);
                let s = (((az + std::f32::consts::PI) / std::f32::consts::TAU * 16.0) as usize).min(15);
                sector[s] += area_km2[c] as f64;
            } else {
                indeg[r] = indeg[r].saturating_add(1);
            }
        }
    }
    let confl = land.iter().filter(|&&c| indeg[c] >= 2).count();
    let mean = sector.iter().sum::<f64>() / 16.0;
    let cv = (sector.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / 16.0).sqrt() / mean.max(1e-9);
    let amax = land.iter().map(|&c| area_km2[c]).fold(0f32, f32::max);
    format!(
        "one receiver on {one}/{} land cells · largest {amax:.1} km² · outlets {outlets} (≥ 1 km²: {big}) · confluences {confl} · sector asymmetry CV {cv:.3}",
        land.len()
    )
}

/// ADR Finding 129-C, the REFERENCES only (rule 15: the gate is written after them): on Finding 126's
/// isotropic synthetics and on the pure dome (ρ = 0, Paik's cone), the D8 tree (the high reference) and
/// the free continuous trace from the same 1–3 km² heads (the low reference), R8 by chord {1, 8, 16, 32}
/// cells, with the D8 tree's concentration and cone rows.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f129_c_refs --nocapture
#[test]
#[ignore]
fn f129_c_refs() {
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 129-C . the references (D8, the free trace), the pure dome  ==========");
    let ns = 2048usize;
    let cell_m = CELL_KM * 1000.0;
    let dom_s = ns as f32 * CELL_KM;
    let b2 = ValleyConstruction { a_min_km2: 1.0, smooth_m: 0.0, ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let sea_m = c1_altitude_norm_to_metres(SEA, &ss);
    for rho in [0.0f32, 0.1, 0.3, 1.0, 3.0] {
        let zm0 = synth_dome(ns, cell_m, rho, 0x9E37_79B9_7F4A_7C15);
        let g = GridF32 { width: ns, height: ns, data: zm0.iter().map(|&m| c1_metres_to_altitude_norm(m, &ss)).collect() };
        let d = c1_drainage_windowed(&g, None, &C1DrainageConfig::default(), &ss, dom_s);
        let bf = breach_monotone(&g, &d.flow.filled, &d.lake_map, SEA, ns, ns);
        let zm: Vec<f32> = bf.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
        let sink: Vec<bool> = bf.data.iter().map(|&v| v <= SEA).collect();
        let sk = skeleton(&g, &b2, &ss, dom_s);
        let paths = branch_paths(&sk);
        let d8r: Vec<String> = [1f32, 8.0, 16.0, 32.0].iter().map(|&k| { let (r, m) = chords_r8_len(&paths, k); format!("c{k:.0} **{r:.3}** ({m})") }).collect();
        let free: Vec<Vec<(f32, f32)>> = paths
            .iter()
            .map(|p| {
                let len: f32 = p.windows(2).map(|q| ((q[1].0 - q[0].0).powi(2) + (q[1].1 - q[0].1).powi(2)).sqrt()).sum();
                free_trace(&zm, ns, ns, cell_m, sea_m, p[0].0, p[0].1, len.max(8.0))
            })
            .collect();
        let frr: Vec<String> = [1f32, 8.0, 16.0, 32.0].iter().map(|&k| { let (r, m) = chords_r8_len(&free, k); format!("c{k:.0} **{r:.3}** ({m})") }).collect();
        eprintln!("\n   ρ = {rho}{} · {} branches (1–3 km²)", if rho == 0.0 { " (the PURE dome, Paik's cone)" } else { "" }, paths.len());
        eprintln!("      HIGH ref, the D8 tree   : {}", d8r.join(" · "));
        eprintln!("      LOW ref, the free trace : {}", frr.join(" · "));
        eprintln!("      the D8 tree: {}", tree_row(&sk.direction, &sk.area_km2, &sink, ns, ns));
    }
}

/// ADR Finding 129-A — the law's reference MEASURED (rule 15): the same slope-area instrument on the
/// construction's own law, `z = base + k·χ` on the S1 skeleton, then WHERE the as-built construction
/// (`carve` on S1, the pipeline's own call) leaves that law on the trunk ≥ 10 km² cells.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f129_a_law --nocapture
#[test]
#[ignore]
fn f129_a_law() {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 129-A . the law measured, and where the as-built construction leaves it  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let sk = skeleton(&s1, &base, &ss, DOMAIN_KM);
    let law: Vec<f32> = (0..n).map(|k| if sk.chi_m[k].is_finite() { sk.floor_m(k, base.age_k) } else { f32::NAN }).collect();
    let recv = |k: usize| -> Option<usize> {
        let d = sk.direction[k];
        (d != DIR_NONE).then(|| {
            ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize * w
                + ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize
        })
    };
    // the trunk links the instrument reads
    let links: Vec<(usize, usize)> = (0..n)
        .filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0)
        .filter_map(|k| recv(k).filter(|&r| sk.trunk[r]).map(|r| (k, r)))
        .collect();
    let base_jumps = links.iter().filter(|&&(k, r)| sk.base_alt_m[k] != sk.base_alt_m[r]).count();
    let (t, c) = theta_sa(&law, &sk, 10.0, 1e-4);
    let (tl, cl) = theta_sa(&law, &sk, 10.0, f32::MIN_POSITIVE);
    eprintln!(
        "   the instrument ON THE LAW (z = base + k·χ on the S1 skeleton, k = {:.5}): θ **{t:.3}** ({c} cells; S > 0: {tl:.3}, {cl}) · {} trunk links · links across two bases {base_jumps}",
        base.age_k,
        links.len()
    );
    let (t, c) = theta_sa(&zs1, &sk, 10.0, 1e-4);
    eprintln!("   the input S1 on the same cells: θ {t:.3} ({c})");
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let worlds: [(&str, ValleyConstruction); 3] = [
        ("témoin C2/10 col", base),
        ("B2 → A_c, no clause", ValleyConstruction { a_min_km2: 0.1, ..base }),
        ("B2 → A_c + CLAUSE", ValleyConstruction { a_min_km2: 0.1, trunk_band: true, ..base }),
    ];
    for (label, vc) in worlds {
        let t = Instant::now();
        let skv = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        // the law does not depend on a_min: the same D8 tree, areas, bases and χ
        let same = (0..n).all(|k| skv.chi_m[k].to_bits() == sk.chi_m[k].to_bits() && skv.base_alt_m[k].to_bits() == sk.base_alt_m[k].to_bits());
        let (built, mk) = carve(&s1, &skv, &vc, &ss);
        drop(skv);
        let zb = metres(&built);
        drop(built);
        let (tb, cb) = theta_sa(&zb, &sk, 10.0, 1e-4);
        // where: the trunk cells, as built against their own law
        let cells: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
        let mut dz: Vec<f32> = cells.iter().map(|&k| zb[k] - law[k]).collect();
        let uncarved = cells.iter().filter(|&&k| !mk.carved[k]).count();
        let unc_below = cells.iter().filter(|&&k| !mk.carved[k] && zs1[k] < law[k]).count();
        let above1 = dz.iter().filter(|&&d| d > 1.0).count();
        let below1 = dz.iter().filter(|&&d| d < -1.0).count();
        let (p10, p25, p50, p75, p90) = (q(&mut dz, 0.1), q(&mut dz, 0.25), q(&mut dz, 0.5), q(&mut dz, 0.75), q(&mut dz, 0.9));
        // θ on a partition of the SAME links (the instrument of `theta_sa`, restricted by a filter)
        let theta_on = |keep: &dyn Fn(usize, usize) -> bool| -> (f32, usize) {
            let v: Vec<(f64, f64)> = links
                .iter()
                .filter(|&&(k, r)| keep(k, r))
                .filter_map(|&(k, r)| {
                    let diag = (k % w != r % w) && (k / w != r / w);
                    let s = (zb[k] - zb[r]) / (sk.cell_m * if diag { std::f32::consts::SQRT_2 } else { 1.0 });
                    (s > 1e-4).then(|| ((sk.area_km2[k] as f64).ln(), (s as f64).ln()))
                })
                .collect();
            let m = v.len() as f64;
            let (mx, my) = (v.iter().map(|a| a.0).sum::<f64>() / m, v.iter().map(|a| a.1).sum::<f64>() / m);
            let sxy: f64 = v.iter().map(|a| (a.0 - mx) * (a.1 - my)).sum();
            let sxx: f64 = v.iter().map(|a| (a.0 - mx).powi(2)).sum();
            ((-sxy / sxx) as f32, v.len())
        };
        let (tall, call) = theta_on(&|_, _| true);
        let (tc, cc) = theta_on(&|k, r| mk.carved[k] && mk.carved[r]);
        let (tu, cu) = theta_on(&|k, r| !(mk.carved[k] && mk.carved[r]));
        let (ta, ca) = theta_on(&|k, r| (zb[k] - law[k]).abs() <= 1.0 && (zb[r] - law[r]).abs() <= 1.0);
        eprintln!(
            "   {label:<22} as built (carve on S1) θ **{tb:.3}** ({cb}) · the law independent of a_min: {same} · {:.0} s",
            t.elapsed().as_secs_f64()
        );
        eprintln!(
            "      {} trunk cells: as built − law p10 {p10:+.1} / p25 {p25:+.1} / p50 {p50:+.1} / p75 {p75:+.1} / p90 {p90:+.1} m · above the law by > 1 m {above1} ({:.1} %) · below by > 1 m {below1} ({:.1} %) · UNCARVED (min kept the field) {uncarved} ({:.1} %), of them the field below its law {unc_below}",
            cells.len(),
            100.0 * above1 as f32 / cells.len() as f32,
            100.0 * below1 as f32 / cells.len() as f32,
            100.0 * uncarved as f32 / cells.len() as f32
        );
        eprintln!(
            "      θ on all the links {tall:.3} ({call}; = theta_sa) · on the links lowered at both ends **{tc:.3}** ({cc}) · on the others **{tu:.3}** ({cu}) · on the links within ±1 m of the law at both ends **{ta:.3}** ({ca})"
        );
    }
    eprintln!("\n==========  end Finding 129-A (law) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 129-A — the source cones, BEFORE / AFTER separated: which heads the confluence clause
/// itself reverses or flattens (a head already reversed without the clause is not the clause's), and
/// the short lines Finding 129-A's census left out (< 11 cells). `carve` on S1, the pipeline's own
/// call (the pipeline's S2 adds only the active rims).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f129_a_cones --nocapture
#[test]
#[ignore]
fn f129_a_cones() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 129-A . the source cones, before / after separated  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let v0 = ValleyConstruction { a_min_km2: 0.1, ..base };
    let v1 = ValleyConstruction { a_min_km2: 0.1, trunk_band: true, ..base };
    let sk = skeleton(&s1, &v0, &ss, DOMAIN_KM);
    let (b0, _) = carve(&s1, &sk, &v0, &ss);
    let (b1, m1) = carve(&s1, &sk, &v1, &ss);
    let moved: Vec<bool> = (0..n).map(|k| b0.data[k] != b1.data[k]).collect();
    let (z0, z1) = (metres(&b0), metres(&b1));
    drop((b0, b1));
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let (mut short_all, mut short_touched, mut touched) = (0usize, 0usize, 0usize);
    let (mut rev_before, mut rev_still, mut rev_by, mut flat_by, mut steep_by, mut rev_all) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut head_dz, mut head_dz_rev) = (Vec::new(), Vec::new());
    let (mut head_on_floor, mut rev_by_head_lowered) = (0usize, 0usize);
    let mut short_head_dz = Vec::new();
    let (mut slope_b, mut slope_a) = (Vec::new(), Vec::new());
    let len_m = 10.0 * CELL_KM * 1000.0; // Finding 129-A's census: head → 10th cell over 10 cells
    for l in &sk.polylines {
        let cells: Vec<usize> = l
            .iter()
            .map(|p| ((p.1.floor() as i64).rem_euclid(h as i64) as usize) * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize)
            .fold(Vec::new(), |mut acc, c| {
                if acc.last() != Some(&c) {
                    acc.push(c);
                }
                acc
            });
        let hit = cells.iter().take(4).any(|&c| moved[c]);
        if cells.len() < 11 {
            short_all += 1;
            if hit {
                short_touched += 1;
                short_head_dz.push(z1[cells[0]] - z0[cells[0]]);
            }
            continue;
        }
        if !hit {
            continue;
        }
        touched += 1;
        let (a, b) = (cells[0], cells[10]);
        let (sb, sa) = (z0[a] - z0[b], z1[a] - z1[b]); // drops over the same 10 cells, m
        slope_b.push(sb / len_m);
        slope_a.push(sa / len_m);
        let dz = z1[a] - z0[a];
        head_dz.push(dz);
        if m1.floor[a] {
            head_on_floor += 1;
        }
        if sa < 0.0 {
            rev_all += 1;
        }
        if sb < 0.0 {
            rev_before += 1;
            if sa < 0.0 {
                rev_still += 1;
            }
        } else if sa < 0.0 {
            rev_by += 1;
            head_dz_rev.push(dz);
            if dz < 0.0 {
                rev_by_head_lowered += 1;
            }
        } else if sb > 0.0 && sa < 0.5 * sb {
            flat_by += 1;
        } else if sa > 1.5 * sb {
            steep_by += 1;
        }
    }
    eprintln!(
        "   {} lines on the A_c skeleton · ≥ 11 cells with a moved cell among the first 4: {touched} (Finding 129-A's census: 3 953 on S2) · < 11 cells: {short_all}, of them touched {short_touched}",
        sk.polylines.len()
    );
    eprintln!(
        "   the touched lines, head → 10th cell: reversed WITH the clause {rev_all} (Finding 129-A: 200 on S2) = already reversed without it {rev_before} (still reversed {rev_still}) + REVERSED BY the clause **{rev_by}** · FLATTENED below half BY the clause **{flat_by}** · steepened above × 1.5 {steep_by}"
    );
    eprintln!(
        "   the head cell's Δz (clause − none): p10 {:+.1} / p50 {:+.1} / p90 {:+.1} m · on a floor band under the clause {head_on_floor} ({:.1} %) · of the heads it reverses: Δz p50 {:+.1} m, lowered {rev_by_head_lowered}",
        q(&mut head_dz, 0.1),
        q(&mut head_dz, 0.5),
        q(&mut head_dz, 0.9),
        100.0 * head_on_floor as f32 / touched.max(1) as f32,
        q(&mut head_dz_rev, 0.5)
    );
    let mean = |v: &[f32]| v.iter().map(|&x| x as f64).sum::<f64>() / v.len().max(1) as f64;
    let (mb, ma) = (mean(&slope_b), mean(&slope_a));
    eprintln!(
        "   the touched lines' head → 10th-cell slope: MEAN {mb:.4} → {ma:.4} m/m ({:+.1} %) · p50 {:.4} → {:.4} · p25 {:.4} → {:.4}",
        100.0 * (ma / mb - 1.0),
        q(&mut slope_b, 0.5),
        q(&mut slope_a, 0.5),
        q(&mut slope_b, 0.25),
        q(&mut slope_a, 0.25)
    );
    eprintln!(
        "   the short touched lines' head Δz: p10 {:+.1} / p50 {:+.1} / p90 {:+.1} m",
        q(&mut short_head_dz, 0.1),
        q(&mut short_head_dz, 0.5),
        q(&mut short_head_dz, 0.9)
    );
    eprintln!("\n==========  end Finding 129-A (cones) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}


// ════════════════════════════════ Finding 130 ════════════════════════════════

/// [O03] Fig. 2h (D8, λ = 0) and Fig. 2l (D8-LTD, λ = 1): the DIAGONAL pointers of the parabolic valley,
/// read from the PDF's vector drawing (column, row from the top, dy; every diagonal goes east, dx = +1).
const O03_FIG2H_DIAG: &[(usize, usize, i32)] = &[(0,7,1), (0,8,1), (0,9,1), (0,10,1), (0,14,-1), (0,15,-1), (0,16,-1), (0,17,-1), (1,7,1), (1,8,1), (1,9,1), (1,10,1), (1,14,-1), (1,15,-1), (1,16,-1), (1,17,-1), (2,7,1), (2,8,1), (2,9,1), (2,10,1), (2,14,-1), (2,15,-1), (2,16,-1), (2,17,-1), (3,7,1), (3,8,1), (3,9,1), (3,10,1), (3,14,-1), (3,15,-1), (3,16,-1), (3,17,-1), (4,7,1), (4,8,1), (4,9,1), (4,10,1), (4,14,-1), (4,15,-1), (4,16,-1), (4,17,-1), (5,7,1), (5,8,1), (5,9,1), (5,10,1), (5,14,-1), (5,15,-1), (5,16,-1), (5,17,-1), (6,7,1), (6,8,1), (6,9,1), (6,10,1), (6,14,-1), (6,15,-1), (6,16,-1), (6,17,-1), (7,7,1), (7,8,1), (7,9,1), (7,10,1), (7,14,-1), (7,15,-1), (7,16,-1), (7,17,-1), (8,7,1), (8,8,1), (8,9,1), (8,10,1), (8,14,-1), (8,15,-1), (8,16,-1), (8,17,-1), (9,7,1), (9,8,1), (9,9,1), (9,10,1), (9,14,-1), (9,15,-1), (9,16,-1), (9,17,-1), (10,7,1), (10,8,1), (10,9,1), (10,10,1), (10,14,-1), (10,15,-1), (10,16,-1), (10,17,-1), (11,7,1), (11,8,1), (11,9,1), (11,10,1), (11,14,-1), (11,15,-1), (11,16,-1), (11,17,-1), (12,7,1), (12,8,1), (12,9,1), (12,10,1), (12,14,-1), (12,15,-1), (12,16,-1), (12,17,-1), (13,7,1), (13,8,1), (13,9,1), (13,10,1), (13,14,-1), (13,15,-1), (13,16,-1), (13,17,-1), (14,7,1), (14,8,1), (14,9,1), (14,10,1), (14,14,-1), (14,15,-1), (14,16,-1), (14,17,-1), (15,7,1), (15,8,1), (15,9,1), (15,10,1), (15,14,-1), (15,15,-1), (15,16,-1), (15,17,-1), (16,7,1), (16,8,1), (16,9,1), (16,10,1), (16,14,-1), (16,15,-1), (16,16,-1), (16,17,-1), (17,7,1), (17,8,1), (17,9,1), (17,10,1), (17,14,-1), (17,15,-1), (17,16,-1), (17,17,-1), (18,7,1), (18,8,1), (18,9,1), (18,10,1), (18,14,-1), (18,15,-1), (18,16,-1), (18,17,-1), (19,7,1), (19,8,1), (19,9,1), (19,10,1), (19,14,-1), (19,15,-1), (19,16,-1), (19,17,-1), (20,7,1), (20,8,1), (20,9,1), (20,10,1), (20,14,-1), (20,15,-1), (20,16,-1), (20,17,-1), (21,7,1), (21,8,1), (21,9,1), (21,10,1), (21,14,-1), (21,15,-1), (21,16,-1), (21,17,-1), (22,7,1), (22,8,1), (22,9,1), (22,10,1), (22,14,-1), (22,15,-1), (22,16,-1), (22,17,-1), (23,7,1), (23,8,1), (23,9,1), (23,10,1), (23,14,-1), (23,15,-1), (23,16,-1), (23,17,-1)];
const O03_FIG2L_DIAG: &[(usize, usize, i32)] = &[(0,2,1), (0,4,1), (0,6,1), (0,8,1), (0,9,1), (0,10,1), (0,14,-1), (0,15,-1), (0,16,-1), (0,18,-1), (0,20,-1), (0,22,-1), (1,2,1), (1,5,1), (1,8,1), (1,9,1), (1,10,1), (1,14,-1), (1,15,-1), (1,16,-1), (1,19,-1), (1,22,-1), (2,2,1), (2,5,1), (2,8,1), (2,9,1), (2,10,1), (2,14,-1), (2,15,-1), (2,16,-1), (2,19,-1), (2,22,-1), (3,2,1), (3,5,1), (3,8,1), (3,9,1), (3,10,1), (3,14,-1), (3,15,-1), (3,16,-1), (3,19,-1), (3,22,-1), (4,2,1), (4,5,1), (4,8,1), (4,9,1), (4,10,1), (4,11,1), (4,13,-1), (4,14,-1), (4,15,-1), (4,16,-1), (4,19,-1), (4,22,-1), (5,2,1), (5,5,1), (5,8,1), (5,9,1), (5,10,1), (5,14,-1), (5,15,-1), (5,16,-1), (5,19,-1), (5,22,-1), (6,2,1), (6,5,1), (6,8,1), (6,9,1), (6,10,1), (6,14,-1), (6,15,-1), (6,16,-1), (6,19,-1), (6,22,-1), (7,2,1), (7,5,1), (7,8,1), (7,9,1), (7,10,1), (7,14,-1), (7,15,-1), (7,16,-1), (7,19,-1), (7,22,-1), (8,2,1), (8,5,1), (8,8,1), (8,9,1), (8,10,1), (8,11,1), (8,13,-1), (8,14,-1), (8,15,-1), (8,16,-1), (8,19,-1), (8,22,-1), (9,2,1), (9,5,1), (9,8,1), (9,9,1), (9,10,1), (9,14,-1), (9,15,-1), (9,16,-1), (9,19,-1), (9,22,-1), (10,2,1), (10,5,1), (10,8,1), (10,9,1), (10,10,1), (10,14,-1), (10,15,-1), (10,16,-1), (10,19,-1), (10,22,-1), (11,2,1), (11,5,1), (11,8,1), (11,9,1), (11,10,1), (11,14,-1), (11,15,-1), (11,16,-1), (11,19,-1), (11,22,-1), (12,2,1), (12,5,1), (12,8,1), (12,9,1), (12,10,1), (12,11,1), (12,13,-1), (12,14,-1), (12,15,-1), (12,16,-1), (12,19,-1), (12,22,-1), (13,2,1), (13,5,1), (13,8,1), (13,9,1), (13,10,1), (13,14,-1), (13,15,-1), (13,16,-1), (13,19,-1), (13,22,-1), (14,2,1), (14,5,1), (14,8,1), (14,9,1), (14,10,1), (14,14,-1), (14,15,-1), (14,16,-1), (14,19,-1), (14,22,-1), (15,2,1), (15,5,1), (15,8,1), (15,9,1), (15,10,1), (15,14,-1), (15,15,-1), (15,16,-1), (15,19,-1), (15,22,-1), (16,2,1), (16,5,1), (16,8,1), (16,9,1), (16,10,1), (16,11,1), (16,13,-1), (16,14,-1), (16,15,-1), (16,16,-1), (16,19,-1), (16,22,-1), (17,2,1), (17,5,1), (17,8,1), (17,9,1), (17,10,1), (17,14,-1), (17,15,-1), (17,16,-1), (17,19,-1), (17,22,-1), (18,2,1), (18,5,1), (18,8,1), (18,9,1), (18,10,1), (18,14,-1), (18,15,-1), (18,16,-1), (18,19,-1), (18,22,-1), (19,2,1), (19,5,1), (19,8,1), (19,9,1), (19,10,1), (19,14,-1), (19,15,-1), (19,16,-1), (19,19,-1), (19,22,-1), (20,2,1), (20,5,1), (20,8,1), (20,9,1), (20,10,1), (20,11,1), (20,13,-1), (20,14,-1), (20,15,-1), (20,16,-1), (20,19,-1), (20,22,-1), (21,2,1), (21,5,1), (21,8,1), (21,9,1), (21,10,1), (21,14,-1), (21,15,-1), (21,16,-1), (21,19,-1), (21,22,-1), (22,2,1), (22,5,1), (22,8,1), (22,9,1), (22,10,1), (22,14,-1), (22,15,-1), (22,16,-1), (22,19,-1), (22,22,-1), (23,2,1), (23,5,1), (23,8,1), (23,9,1), (23,10,1), (23,14,-1), (23,15,-1), (23,16,-1), (23,19,-1), (23,22,-1)];

/// ADR Finding 130-B — B1 on [O03]'s parabolic valley as its figures define it (declared in
/// `f130_predictions.md` before the run): 25 × 25 cells, z = (25 − x) + q (y − 12.5)² in cell units,
/// q = 0.244 (the D8 interval's midpoint; 0.2205 and 0.2670 beside it), a CLOSED boundary with one outlet
/// east of (24, 12). Theoretical areas: the analytic flow-tube integral (C = (y − yc)·e^{2qx} constant on a
/// drainage line); the F129 8 × 8 tracing beside it. B0: our D8 must reproduce Fig. 2h's 625 pointers.
/// B1: D8's ME / MAE / RMSE within ±0.017 / ±0.034 / ±0.064 of −0.038 / 0.566 / 2.043, else STOP.
/// B2: LTD's within the same tolerance of −0.261 / 0.282 / 0.330; the ratio reported, not a gate.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_b --nocapture
#[test]
#[ignore]
fn f130_b() {
    use ymir_core::tectonics_c1::valley_construction::{accumulate_cells, ltd_directions_masked};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    eprintln!("\n==========  Finding 130-B . [O03]'s parabolic valley as its figures define it  ==========");
    const N: usize = 25;
    const P: usize = N + 2; // a ring of walls, one outlet
    let kdir = |dx: i32, dy: i32| (0..8u8).find(|&k| D8_DX[k as usize] == dx && D8_DY[k as usize] == dy).expect("a D8 step");
    let pad = |i: usize, j: usize| (j + 1) * P + (i + 1);
    let outlet = pad(N, 12); // east of (24, 12)
    let active: Vec<bool> = (0..P * P)
        .map(|k| {
            let (x, y) = (k % P, k / P);
            ((1..=N).contains(&x) && (1..=N).contains(&y)) || k == outlet
        })
        .collect();
    let sink: Vec<bool> = (0..P * P).map(|k| k == outlet).collect();
    let yc = 12.5f64;
    // the expected D8 of Fig. 2h
    let fig_d8 = |i: usize, j: usize| -> u8 {
        if let Some(&(_, _, dy)) = O03_FIG2H_DIAG.iter().find(|t| t.0 == i && t.1 == j) {
            return kdir(1, dy);
        }
        if i == N - 1 {
            return if j == 12 { kdir(1, 0) } else if j < 12 { kdir(0, 1) } else { kdir(0, -1) };
        }
        if (11..=13).contains(&j) {
            kdir(1, 0)
        } else if j < 12 {
            kdir(0, 1)
        } else {
            kdir(0, -1)
        }
    };
    // the analytic theoretical area of cell (i, j), in cells
    let area_t = |q: f64, i: usize, j: usize| -> f64 {
        let cc = |x: f64, y: f64| (y - yc) * (2.0 * q * x).exp();
        let corners = [
            cc(i as f64, j as f64),
            cc(i as f64 + 1.0, j as f64),
            cc(i as f64, j as f64 + 1.0),
            cc(i as f64 + 1.0, j as f64 + 1.0),
        ];
        let (cmin, cmax) = corners.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &c| (a.0.min(c), a.1.max(c)));
        let m = 20_000usize;
        let dc = (cmax - cmin) / m as f64;
        let mut a = 0.0;
        for s in 0..m {
            let c = cmin + (s as f64 + 0.5) * dc;
            // the band of the cell seen from the line's side of the axis
            let (lo, hi) = if c >= 0.0 { (j as f64 - yc, j as f64 + 1.0 - yc) } else { (yc - j as f64 - 1.0, yc - j as f64) };
            let ca = c.abs();
            if hi <= 0.0 {
                continue;
            }
            let xa = if ca == 0.0 { f64::NEG_INFINITY } else { (ca / hi).ln() / (2.0 * q) };
            let xb = if lo <= 0.0 { f64::INFINITY } else if ca == 0.0 { f64::NEG_INFINITY } else { (ca / lo).ln() / (2.0 * q) };
            let (u, v) = (xa.max(i as f64), xb.min(i as f64 + 1.0));
            if u > v {
                continue;
            }
            let x_lo = if ca == 0.0 { 0.0 } else { ((ca / (N as f64 / 2.0)).ln() / (2.0 * q)).max(0.0) };
            if x_lo >= v {
                continue;
            }
            a += ((-2.0 * q * x_lo).exp() - (-2.0 * q * v).exp()) / (2.0 * q) * dc;
        }
        a
    };
    // F129's form: 8 × 8 sub-points traced down the exact gradient, counted in every cell they cross
    let area_traced = |q: f64| -> Vec<f64> {
        let mut at = vec![0f64; N * N];
        let sub = 8usize;
        for c in 0..N * N {
            for sy in 0..sub {
                for sx in 0..sub {
                    let mut px = (c % N) as f64 + (sx as f64 + 0.5) / sub as f64;
                    let mut py = (c / N) as f64 + (sy as f64 + 0.5) / sub as f64;
                    let mut last = usize::MAX;
                    while px < N as f64 && (0.0..N as f64).contains(&py) {
                        let cc = py as usize * N + px as usize;
                        if cc != last {
                            at[cc] += 1.0 / (sub * sub) as f64;
                            last = cc;
                        }
                        let (gx, gy) = (1.0f64, -2.0 * q * (py - yc));
                        let g = (gx * gx + gy * gy).sqrt();
                        px += 0.01 * gx / g;
                        py += 0.01 * gy / g;
                    }
                }
            }
        }
        at
    };
    let stats = |a: &[f32], at: &[f64]| -> (f64, f64, f64) {
        let e: Vec<f64> = (0..N * N).map(|c| (a[pad(c % N, c / N)] as f64 - at[c]) / at[c]).collect();
        let k = e.len() as f64;
        (
            e.iter().sum::<f64>() / k,
            e.iter().map(|v| v.abs()).sum::<f64>() / k,
            (e.iter().map(|v| v * v).sum::<f64>() / k).sqrt(),
        )
    };
    let within = |v: (f64, f64, f64), r: (f64, f64, f64)| {
        (v.0 - r.0).abs() <= 0.017 && (v.1 - r.1).abs() <= 0.034 && (v.2 - r.2).abs() <= 0.064
    };
    for q in [0.244f64, 0.2205, 0.2670] {
        let z: Vec<f64> = (0..P * P)
            .map(|k| {
                let (x, y) = (k % P, k / P);
                if k == outlet {
                    return -1.0;
                }
                if !active[k] {
                    return 1e6;
                }
                let (cx, cy) = ((x - 1) as f64 + 0.5, (y - 1) as f64 + 0.5);
                (N as f64 - cx) + q * (cy - yc).powi(2)
            })
            .collect();
        // D8, on the domain only (the ring is a wall)
        let d8: Vec<u8> = (0..P * P)
            .map(|k| {
                if !active[k] || sink[k] {
                    return DIR_NONE;
                }
                let (x, y) = ((k % P) as i32, (k / P) as i32);
                let mut best = (DIR_NONE, 0f64);
                for kk in 0..8u8 {
                    let m = (y + D8_DY[kk as usize]) as usize * P + (x + D8_DX[kk as usize]) as usize;
                    if !active[m] {
                        continue;
                    }
                    let s = (z[k] - z[m]) / if kk % 2 == 1 { std::f64::consts::SQRT_2 } else { 1.0 };
                    if s > best.1 {
                        best = (kk, s);
                    }
                }
                best.0
            })
            .collect();
        let d8_match = (0..N * N).filter(|&c| d8[pad(c % N, c / N)] == fig_d8(c % N, c / N)).count();
        let (ltd, st) = ltd_directions_masked(&z, &sink, Some(&active), &d8, P, P, 1.0);
        let a8 = accumulate_cells(&d8, &sink, P, P);
        let al = accumulate_cells(&ltd, &sink, P, P);
        let at: Vec<f64> = (0..N * N).map(|c| area_t(q, c % N, c / N)).collect();
        let atr = area_traced(q);
        let (e8, el) = (stats(&a8, &at), stats(&al, &at));
        let (e8t, elt) = (stats(&a8, &atr), stats(&al, &atr));
        // LTD against Fig. 2l, diagonal by diagonal
        let ours: Vec<(usize, usize, i32)> = (0..N * N)
            .filter_map(|c| {
                let d = ltd[pad(c % N, c / N)];
                (d != DIR_NONE && D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0).then(|| (c % N, c / N, D8_DY[d as usize]))
            })
            .collect();
        let both = ours.iter().filter(|t| O03_FIG2L_DIAG.contains(t)).count();
        let pass1 = d8_match == N * N && within(e8, (-0.038, 0.566, 2.043));
        let pass2 = within(el, (-0.261, 0.282, 0.330));
        let tag = if q == 0.244 { " (DECLARED)" } else { " (sensitivity)" };
        eprintln!("\n   q = {q:.4}{tag}");
        eprintln!("      B0 · our D8 against Fig. 2h: {d8_match} / 625 pointers identical · LTD flat fallbacks {}", st.flats);
        eprintln!(
            "      D8    : ME {:+.3} · MAE **{:.3}** · RMSE **{:.3}**   ([O03] −0.038 / 0.566 / 2.043) · traced areas: {:+.3} / {:.3} / {:.3}",
            e8.0, e8.1, e8.2, e8t.0, e8t.1, e8t.2
        );
        eprintln!(
            "      D8-LTD: ME {:+.3} · MAE **{:.3}** · RMSE **{:.3}**   ([O03] −0.261 / 0.282 / 0.330) · traced areas: {:+.3} / {:.3} / {:.3}",
            el.0, el.1, el.2, elt.0, elt.1, elt.2
        );
        eprintln!(
            "      ratios LTD / D8 (reported, not a gate): MAE {:.3} (paper {:.3}) · RMSE {:.3} (paper {:.3})",
            el.1 / e8.1,
            0.282 / 0.566,
            el.2 / e8.2,
            0.330 / 2.043
        );
        eprintln!(
            "      LTD's diagonals against Fig. 2l: ours {} · the figure's {} · the same cell and direction {both}",
            ours.len(),
            O03_FIG2L_DIAG.len()
        );
        if q == 0.244 {
            let b2 = if !pass1 {
                "NOT JUDGED (B1 failed)"
            } else if pass2 {
                "PASS"
            } else {
                "FAIL"
            };
            eprintln!("      → B1 (D8) **{}** · B2 (LTD, absolute) **{b2}**", if pass1 { "PASS" } else { "FAIL" });
        }
    }
    eprintln!("\n==========  end Finding 130-B  ==========\n");
}

/// ADR Finding 130-B, diagnostic after the gate (no gate re-run): WHERE our D8 leaves Fig. 2h, and which
/// cells carry the D8 RMSE, at the declared q = 0.244.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_b_where --nocapture
#[test]
#[ignore]
fn f130_b_where() {
    use ymir_core::tectonics_c1::valley_construction::accumulate_cells;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    const N: usize = 25;
    const P: usize = N + 2;
    let q = 0.244f64;
    let yc = 12.5f64;
    let kdir = |dx: i32, dy: i32| (0..8u8).find(|&k| D8_DX[k as usize] == dx && D8_DY[k as usize] == dy).expect("a D8 step");
    let pad = |i: usize, j: usize| (j + 1) * P + (i + 1);
    let outlet = pad(N, 12);
    let active: Vec<bool> = (0..P * P)
        .map(|k| {
            let (x, y) = (k % P, k / P);
            ((1..=N).contains(&x) && (1..=N).contains(&y)) || k == outlet
        })
        .collect();
    let sink: Vec<bool> = (0..P * P).map(|k| k == outlet).collect();
    let fig_d8 = |i: usize, j: usize| -> u8 {
        if let Some(&(_, _, dy)) = O03_FIG2H_DIAG.iter().find(|t| t.0 == i && t.1 == j) {
            return kdir(1, dy);
        }
        if i == N - 1 {
            return if j == 12 { kdir(1, 0) } else if j < 12 { kdir(0, 1) } else { kdir(0, -1) };
        }
        if (11..=13).contains(&j) {
            kdir(1, 0)
        } else if j < 12 {
            kdir(0, 1)
        } else {
            kdir(0, -1)
        }
    };
    let z: Vec<f64> = (0..P * P)
        .map(|k| {
            let (x, y) = (k % P, k / P);
            if k == outlet {
                return -1.0;
            }
            if !active[k] {
                return 1e6;
            }
            (N as f64 - ((x - 1) as f64 + 0.5)) + q * (((y - 1) as f64 + 0.5) - yc).powi(2)
        })
        .collect();
    let d8: Vec<u8> = (0..P * P)
        .map(|k| {
            if !active[k] || sink[k] {
                return DIR_NONE;
            }
            let (x, y) = ((k % P) as i32, (k / P) as i32);
            let mut best = (DIR_NONE, 0f64);
            for kk in 0..8u8 {
                let m = (y + D8_DY[kk as usize]) as usize * P + (x + D8_DX[kk as usize]) as usize;
                if !active[m] {
                    continue;
                }
                let s = (z[k] - z[m]) / if kk % 2 == 1 { std::f64::consts::SQRT_2 } else { 1.0 };
                if s > best.1 {
                    best = (kk, s);
                }
            }
            best.0
        })
        .collect();
    eprintln!("\n==========  Finding 130-B . where our D8 leaves Fig. 2h (q = 0.244)  ==========");
    for c in 0..N * N {
        let (i, j) = (c % N, c / N);
        let (ours, fig) = (d8[pad(i, j)], fig_d8(i, j));
        if ours != fig {
            eprintln!(
                "   ({i},{j}): ours ({},{}) · the figure's ({},{})",
                D8_DX[ours as usize], D8_DY[ours as usize], D8_DX[fig as usize], D8_DY[fig as usize]
            );
        }
    }
    let a8 = accumulate_cells(&d8, &sink, P, P);
    let mut rows: Vec<(f64, usize, usize, f32)> = (0..N * N)
        .map(|c| {
            let (i, j) = (c % N, c / N);
            (0.0, i, j, a8[pad(i, j)])
        })
        .collect();
    // the traced theoretical area (the analytic and the traced agree to 0.01 in MAE / RMSE)
    let mut at = vec![0f64; N * N];
    let sub = 8usize;
    for c in 0..N * N {
        for sy in 0..sub {
            for sx in 0..sub {
                let mut px = (c % N) as f64 + (sx as f64 + 0.5) / sub as f64;
                let mut py = (c / N) as f64 + (sy as f64 + 0.5) / sub as f64;
                let mut last = usize::MAX;
                while px < N as f64 && (0.0..N as f64).contains(&py) {
                    let cc = py as usize * N + px as usize;
                    if cc != last {
                        at[cc] += 1.0 / (sub * sub) as f64;
                        last = cc;
                    }
                    let (gx, gy) = (1.0f64, -2.0 * q * (py - yc));
                    let g = (gx * gx + gy * gy).sqrt();
                    px += 0.01 * gx / g;
                    py += 0.01 * gy / g;
                }
            }
        }
    }
    for r in rows.iter_mut() {
        r.0 = (r.3 as f64 - at[r.2 * N + r.1]) / at[r.2 * N + r.1];
    }
    rows.sort_by(|a, b| b.0.abs().total_cmp(&a.0.abs()));
    let tot: f64 = rows.iter().map(|r| r.0 * r.0).sum();
    eprintln!("   the 12 largest |ε| (share of Σε²):");
    for r in rows.iter().take(12) {
        eprintln!("      ({},{}) A {:.0} · A_t {:.2} · ε {:+.2} · {:.1} %", r.1, r.2, r.3, at[r.2 * N + r.1], r.0, 100.0 * r.0 * r.0 / tot);
    }
    eprintln!("==========  end  ==========\n");
}
/// ADR Finding 130-B-DIAG (after B failed; NOT a gate): the outlet admits (24, 12) only, as Fig. 2h draws it. Otherwise f130_b: [O03]'s parabolic valley as its figures define it (declared in
/// `f130_predictions.md` before the run): 25 × 25 cells, z = (25 − x) + q (y − 12.5)² in cell units,
/// q = 0.244 (the D8 interval's midpoint; 0.2205 and 0.2670 beside it), a CLOSED boundary with one outlet
/// east of (24, 12). Theoretical areas: the analytic flow-tube integral (C = (y − yc)·e^{2qx} constant on a
/// drainage line); the F129 8 × 8 tracing beside it. B0: our D8 must reproduce Fig. 2h's 625 pointers.
/// B1: D8's ME / MAE / RMSE within ±0.017 / ±0.034 / ±0.064 of −0.038 / 0.566 / 2.043, else STOP.
/// B2: LTD's within the same tolerance of −0.261 / 0.282 / 0.330; the ratio reported, not a gate.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_b_diag --nocapture
#[test]
#[ignore]
fn f130_b_diag() {
    use ymir_core::tectonics_c1::valley_construction::{accumulate_cells, ltd_directions_masked};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    eprintln!("\n==========  Finding 130-B-DIAG . the single-inflow outlet, [O03]'s parabolic valley as its figures define it  ==========");
    const N: usize = 25;
    const P: usize = N + 2; // a ring of walls, one outlet
    let kdir = |dx: i32, dy: i32| (0..8u8).find(|&k| D8_DX[k as usize] == dx && D8_DY[k as usize] == dy).expect("a D8 step");
    let pad = |i: usize, j: usize| (j + 1) * P + (i + 1);
    let outlet = pad(N, 12); // east of (24, 12)
    let active: Vec<bool> = (0..P * P)
        .map(|k| {
            let (x, y) = (k % P, k / P);
            ((1..=N).contains(&x) && (1..=N).contains(&y)) || k == outlet
        })
        .collect();
    let sink: Vec<bool> = (0..P * P).map(|k| k == outlet).collect();
    let yc = 12.5f64;
    // the expected D8 of Fig. 2h
    let fig_d8 = |i: usize, j: usize| -> u8 {
        if let Some(&(_, _, dy)) = O03_FIG2H_DIAG.iter().find(|t| t.0 == i && t.1 == j) {
            return kdir(1, dy);
        }
        if i == N - 1 {
            return if j == 12 { kdir(1, 0) } else if j < 12 { kdir(0, 1) } else { kdir(0, -1) };
        }
        if (11..=13).contains(&j) {
            kdir(1, 0)
        } else if j < 12 {
            kdir(0, 1)
        } else {
            kdir(0, -1)
        }
    };
    // the analytic theoretical area of cell (i, j), in cells
    let area_t = |q: f64, i: usize, j: usize| -> f64 {
        let cc = |x: f64, y: f64| (y - yc) * (2.0 * q * x).exp();
        let corners = [
            cc(i as f64, j as f64),
            cc(i as f64 + 1.0, j as f64),
            cc(i as f64, j as f64 + 1.0),
            cc(i as f64 + 1.0, j as f64 + 1.0),
        ];
        let (cmin, cmax) = corners.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &c| (a.0.min(c), a.1.max(c)));
        let m = 20_000usize;
        let dc = (cmax - cmin) / m as f64;
        let mut a = 0.0;
        for s in 0..m {
            let c = cmin + (s as f64 + 0.5) * dc;
            // the band of the cell seen from the line's side of the axis
            let (lo, hi) = if c >= 0.0 { (j as f64 - yc, j as f64 + 1.0 - yc) } else { (yc - j as f64 - 1.0, yc - j as f64) };
            let ca = c.abs();
            if hi <= 0.0 {
                continue;
            }
            let xa = if ca == 0.0 { f64::NEG_INFINITY } else { (ca / hi).ln() / (2.0 * q) };
            let xb = if lo <= 0.0 { f64::INFINITY } else if ca == 0.0 { f64::NEG_INFINITY } else { (ca / lo).ln() / (2.0 * q) };
            let (u, v) = (xa.max(i as f64), xb.min(i as f64 + 1.0));
            if u > v {
                continue;
            }
            let x_lo = if ca == 0.0 { 0.0 } else { ((ca / (N as f64 / 2.0)).ln() / (2.0 * q)).max(0.0) };
            if x_lo >= v {
                continue;
            }
            a += ((-2.0 * q * x_lo).exp() - (-2.0 * q * v).exp()) / (2.0 * q) * dc;
        }
        a
    };
    // F129's form: 8 × 8 sub-points traced down the exact gradient, counted in every cell they cross
    let area_traced = |q: f64| -> Vec<f64> {
        let mut at = vec![0f64; N * N];
        let sub = 8usize;
        for c in 0..N * N {
            for sy in 0..sub {
                for sx in 0..sub {
                    let mut px = (c % N) as f64 + (sx as f64 + 0.5) / sub as f64;
                    let mut py = (c / N) as f64 + (sy as f64 + 0.5) / sub as f64;
                    let mut last = usize::MAX;
                    while px < N as f64 && (0.0..N as f64).contains(&py) {
                        let cc = py as usize * N + px as usize;
                        if cc != last {
                            at[cc] += 1.0 / (sub * sub) as f64;
                            last = cc;
                        }
                        let (gx, gy) = (1.0f64, -2.0 * q * (py - yc));
                        let g = (gx * gx + gy * gy).sqrt();
                        px += 0.01 * gx / g;
                        py += 0.01 * gy / g;
                    }
                }
            }
        }
        at
    };
    let stats = |a: &[f32], at: &[f64]| -> (f64, f64, f64) {
        let e: Vec<f64> = (0..N * N).map(|c| (a[pad(c % N, c / N)] as f64 - at[c]) / at[c]).collect();
        let k = e.len() as f64;
        (
            e.iter().sum::<f64>() / k,
            e.iter().map(|v| v.abs()).sum::<f64>() / k,
            (e.iter().map(|v| v * v).sum::<f64>() / k).sqrt(),
        )
    };
    let within = |v: (f64, f64, f64), r: (f64, f64, f64)| {
        (v.0 - r.0).abs() <= 0.017 && (v.1 - r.1).abs() <= 0.034 && (v.2 - r.2).abs() <= 0.064
    };
    for q in [0.244f64, 0.2205, 0.2670] {
        let z: Vec<f64> = (0..P * P)
            .map(|k| {
                let (x, y) = (k % P, k / P);
                if k == outlet {
                    // (24, 12) is at 0.5; the diagonal from (24, 11) breaks even at 0.5 + q - sqrt2 q
                    return 0.5 - 0.2 * q;
                }
                if !active[k] {
                    return 1e6;
                }
                let (cx, cy) = ((x - 1) as f64 + 0.5, (y - 1) as f64 + 0.5);
                (N as f64 - cx) + q * (cy - yc).powi(2)
            })
            .collect();
        // D8, on the domain only (the ring is a wall)
        let d8: Vec<u8> = (0..P * P)
            .map(|k| {
                if !active[k] || sink[k] {
                    return DIR_NONE;
                }
                let (x, y) = ((k % P) as i32, (k / P) as i32);
                let mut best = (DIR_NONE, 0f64);
                for kk in 0..8u8 {
                    let m = (y + D8_DY[kk as usize]) as usize * P + (x + D8_DX[kk as usize]) as usize;
                    if !active[m] {
                        continue;
                    }
                    let s = (z[k] - z[m]) / if kk % 2 == 1 { std::f64::consts::SQRT_2 } else { 1.0 };
                    if s > best.1 {
                        best = (kk, s);
                    }
                }
                best.0
            })
            .collect();
        let d8_match = (0..N * N).filter(|&c| d8[pad(c % N, c / N)] == fig_d8(c % N, c / N)).count();
        let (ltd, st) = ltd_directions_masked(&z, &sink, Some(&active), &d8, P, P, 1.0);
        let a8 = accumulate_cells(&d8, &sink, P, P);
        let al = accumulate_cells(&ltd, &sink, P, P);
        let at: Vec<f64> = (0..N * N).map(|c| area_t(q, c % N, c / N)).collect();
        let atr = area_traced(q);
        let (e8, el) = (stats(&a8, &at), stats(&al, &at));
        let (e8t, elt) = (stats(&a8, &atr), stats(&al, &atr));
        // LTD against Fig. 2l, diagonal by diagonal
        let ours: Vec<(usize, usize, i32)> = (0..N * N)
            .filter_map(|c| {
                let d = ltd[pad(c % N, c / N)];
                (d != DIR_NONE && D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0).then(|| (c % N, c / N, D8_DY[d as usize]))
            })
            .collect();
        let both = ours.iter().filter(|t| O03_FIG2L_DIAG.contains(t)).count();
        let pass1 = d8_match == N * N && within(e8, (-0.038, 0.566, 2.043));
        let pass2 = within(el, (-0.261, 0.282, 0.330));
        let tag = if q == 0.244 { " (DECLARED)" } else { " (sensitivity)" };
        eprintln!("\n   q = {q:.4}{tag}");
        eprintln!("      B0 · our D8 against Fig. 2h: {d8_match} / 625 pointers identical · LTD flat fallbacks {}", st.flats);
        eprintln!(
            "      D8    : ME {:+.3} · MAE **{:.3}** · RMSE **{:.3}**   ([O03] −0.038 / 0.566 / 2.043) · traced areas: {:+.3} / {:.3} / {:.3}",
            e8.0, e8.1, e8.2, e8t.0, e8t.1, e8t.2
        );
        eprintln!(
            "      D8-LTD: ME {:+.3} · MAE **{:.3}** · RMSE **{:.3}**   ([O03] −0.261 / 0.282 / 0.330) · traced areas: {:+.3} / {:.3} / {:.3}",
            el.0, el.1, el.2, elt.0, elt.1, elt.2
        );
        eprintln!(
            "      ratios LTD / D8 (reported, not a gate): MAE {:.3} (paper {:.3}) · RMSE {:.3} (paper {:.3})",
            el.1 / e8.1,
            0.282 / 0.566,
            el.2 / e8.2,
            0.330 / 2.043
        );
        eprintln!(
            "      LTD's diagonals against Fig. 2l: ours {} · the figure's {} · the same cell and direction {both}",
            ours.len(),
            O03_FIG2L_DIAG.len()
        );
        if q == 0.244 {
            let b2 = if !pass1 {
                "NOT JUDGED (B1 failed)"
            } else if pass2 {
                "PASS"
            } else {
                "FAIL"
            };
            eprintln!("      (diagnostic, not a gate) the D8 criterion would read **{}** · LTD **{b2}**", if pass1 { "PASS" } else { "FAIL" });
        }
    }
    eprintln!("\n==========  end Finding 130-B  ==========\n");
}

/// θ over the trunk links LOWERED AT BOTH ENDS (Finding 129's laid floor) with a percentile bootstrap over
/// the links (1 000 resamples, fixed seed): `(θ, lo 2.5 %, hi 97.5 %, links)`.
fn theta_laid_ci(zb: &[f32], carved: &[bool], sk: &Skeleton) -> (f32, f32, f32, usize) {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let (w, h) = (sk.width, sk.height);
    let mut v: Vec<(f64, f64)> = Vec::new();
    for k in 0..w * h {
        if !sk.trunk[k] || sk.area_km2[k] < 10.0 || sk.direction[k] == DIR_NONE {
            continue;
        }
        let d = sk.direction[k] as usize;
        let r = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize * w + ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
        if !sk.trunk[r] || !carved[k] || !carved[r] {
            continue;
        }
        let dist = sk.cell_m * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
        let s = (zb[k] - zb[r]) / dist;
        if s > 1e-4 {
            v.push(((sk.area_km2[k] as f64).ln(), (s as f64).ln()));
        }
    }
    let fit = |idx: &mut dyn Iterator<Item = usize>| -> f64 {
        let (mut n, mut sx, mut sy, mut sxx, mut sxy) = (0f64, 0f64, 0f64, 0f64, 0f64);
        for i in idx {
            let (x, y) = v[i];
            n += 1.0;
            sx += x;
            sy += y;
            sxx += x * x;
            sxy += x * y;
        }
        -((sxy - sx * sy / n) / (sxx - sx * sx / n))
    };
    let t = fit(&mut (0..v.len()));
    let mut rng = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    let m = v.len();
    let mut bs: Vec<f64> = (0..1000)
        .map(|_| {
            let picks: Vec<usize> = (0..m).map(|_| (next() % m as u64) as usize).collect();
            fit(&mut picks.into_iter())
        })
        .collect();
    bs.sort_by(f64::total_cmp);
    (t as f32, bs[25] as f32, bs[974] as f32, m)
}

/// ADR Finding 130-A0 / A1 / A2 (the law's part) — on `carve(S1)`, the pipeline's own call: WHERE the full
/// clause's raised trunk cells sit relative to the junction (A0), the references with a CI (A1: no clause, the
/// full clause), and the restricted clause on B2 → A_c and B2 + B1 + foot (A2's raised cells and θ laid).
/// Definitions declared in `f130_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_a01 --nocapture
#[test]
#[ignore]
fn f130_a01() {
    use ymir_core::tectonics_c1::valley_construction::carve_diag;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 130-A0/A1/A2 . the restricted clause, on the trunks (carve on S1)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let prof = WallProfile { foot_quiet: Some(0.5), ..WallProfile::f124() };
    let sk = skeleton(&s1, &base, &ss, DOMAIN_KM);
    let law: Vec<f32> = (0..n).map(|k| if sk.chi_m[k].is_finite() { sk.floor_m(k, base.age_k) } else { f32::NAN }).collect();
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let sk_ac = skeleton(&s1, &ac, &ss, DOMAIN_KM);
    let raised_of = |zb: &[f32]| -> Vec<bool> {
        let mut r = vec![false; n];
        for &k in &trunk10 {
            r[k] = zb[k] - law[k] > 1.0;
        }
        r
    };
    let row = |label: &str, zb: &[f32], carved: &[bool], raised: &[bool]| {
        let (t, lo, hi, m) = theta_laid_ci(zb, carved, &sk);
        let nr = trunk10.iter().filter(|&&k| raised[k]).count();
        eprintln!("   {label:<40} raised > 1 m above the law **{nr}** · θ laid **{t:.3}** [95 % {lo:.3}, {hi:.3}] ({m} links)");
        (nr, t, lo, hi)
    };
    // the witness
    let (bt, mt) = carve(&s1, &sk, &base, &ss);
    let zt = metres(&bt);
    drop(bt);
    row("témoin C2/10 col", &zt, &mt.carved, &raised_of(&zt));
    drop((zt, mt));
    for (tag, world) in [("B2 → A_c", ac), ("B2 (A_c) + B1 + foot", ValleyConstruction { wall_profile: Some(prof), ..ac })] {
        eprintln!("\n   ── {tag} ──");
        let (b0, m0) = carve(&s1, &sk_ac, &world, &ss);
        let z0 = metres(&b0);
        drop(b0);
        let r0 = raised_of(&z0);
        row("A1 · no clause (the floor)", &z0, &m0.carved, &r0);
        let (b1, m1, dg) = carve_diag(&s1, &sk_ac, &ValleyConstruction { trunk_band: true, ..world }, &ss);
        let z1 = metres(&b1);
        drop(b1);
        let r1 = raised_of(&z1);
        row("A1 · the FULL clause (Finding 128)", &z1, &m1.carved, &r1);
        let (b2, m2) = carve(&s1, &sk_ac, &ValleyConstruction { trunk_band: true, trunk_band_downstream: true, ..world }, &ss);
        let z2 = metres(&b2);
        drop(b2);
        let r2 = raised_of(&z2);
        row("A2 · the RESTRICTED clause", &z2, &m2.carved, &r2);
        if tag != "B2 → A_c" {
            continue;
        }
        // ── A0: where the full clause's own raised cells sit ──
        let pop: Vec<usize> = trunk10.iter().copied().filter(|&k| r1[k] && !r0[k]).collect();
        let banded = dg.banded.as_ref().expect("the clause is on");
        let junction_on = |small: u32, big: u32| -> Option<u32> {
            let mut l = small;
            for _ in 0..sk_ac.polylines.len() {
                match sk_ac.line_parent[l as usize] {
                    Some((p, j)) if p == big => return Some(j),
                    Some((p, _)) => l = p,
                    None => return None,
                }
            }
            None
        };
        let (mut no_owner, mut no_junction, mut up, mut down) = (0usize, 0usize, 0usize, 0usize);
        let (mut dist_w, mut dist_m, mut d_floor_j, mut d_floor_law) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let (mut within_1w, mut higher_than_j) = (0usize, 0usize);
        for &k in &pop {
            let s = banded[k];
            if s == u32::MAX {
                continue; // not laid by the clause: moved by a neighbour's change (counted below)
            }
            let who = dg.who[k];
            let big = dg.line_of[s as usize];
            let ps = dg.pos_of[s as usize] as usize;
            let line = &sk_ac.polylines[big as usize];
            d_floor_law.push(line[ps].2 - law[k]);
            if who == u32::MAX {
                no_owner += 1;
                continue;
            }
            let Some(j) = junction_on(dg.line_of[who as usize], big) else {
                no_junction += 1;
                continue;
            };
            let j = j as usize;
            if ps < j {
                up += 1;
            } else {
                down += 1;
            }
            let (a, b) = (ps.min(j), ps.max(j));
            let arc: f32 = (a..b).map(|i| ((line[i + 1].0 - line[i].0).powi(2) + (line[i + 1].1 - line[i].1).powi(2)).sqrt()).sum::<f32>() * sk.cell_m;
            let wj = 2.0 * line[j].3;
            dist_m.push(arc);
            dist_w.push(arc / wj);
            if arc <= wj {
                within_1w += 1;
            }
            let df = line[ps].2 - line[j].2;
            d_floor_j.push(df);
            if df > 0.0 {
                higher_than_j += 1;
            }
        }
        let not_laid = pop.iter().filter(|&&k| banded[k] == u32::MAX).count();
        let q = |v: &mut Vec<f32>, p: f32| -> f32 {
            if v.is_empty() {
                return f32::NAN;
            }
            v.sort_by(f32::total_cmp);
            v[((v.len() - 1) as f32 * p) as usize]
        };
        let with_j = up + down;
        eprintln!("\n   ── A0: the full clause's own raised cells (> 1 m above their law with it, ≤ 1 m without): {} ──", pop.len());
        eprintln!(
            "      laid by a larger line's sample {} · not laid by the clause (a neighbour moved them) {not_laid} · no owner {no_owner} · the owner's chain never meets the larger line (NO junction) **{no_junction}** · with a junction {with_j}",
            pop.len() - not_laid
        );
        eprintln!(
            "      of those with a junction: UPSTREAM **{up}** ({:.1} %) · at or downstream **{down}** ({:.1} %) · within 1 W of the junction {within_1w} ({:.1} %)",
            100.0 * up as f32 / with_j.max(1) as f32,
            100.0 * down as f32 / with_j.max(1) as f32,
            100.0 * within_1w as f32 / with_j.max(1) as f32
        );
        eprintln!(
            "      distance to the junction along the larger line: p10 / p50 / p90 {:.0} / {:.0} / {:.0} m = {:.2} / {:.2} / {:.2} W",
            q(&mut dist_m, 0.1),
            q(&mut dist_m, 0.5),
            q(&mut dist_m, 0.9),
            q(&mut dist_w, 0.1),
            q(&mut dist_w, 0.5),
            q(&mut dist_w, 0.9)
        );
        eprintln!(
            "      the laying sample's floor − the junction's floor: p10 / p50 / p90 {:+.1} / {:+.1} / {:+.1} m · higher than the junction {higher_than_j} ({:.1} %) · the laying floor − the cell's own law p10 / p50 / p90 {:+.1} / {:+.1} / {:+.1} m",
            q(&mut d_floor_j, 0.1),
            q(&mut d_floor_j, 0.5),
            q(&mut d_floor_j, 0.9),
            100.0 * higher_than_j as f32 / with_j.max(1) as f32,
            q(&mut d_floor_law, 0.1),
            q(&mut d_floor_law, 0.5),
            q(&mut d_floor_law, 0.9)
        );
    }
    eprintln!("\n==========  end Finding 130-A0/A1/A2 (law) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 130-P — the flats of the field the skeleton reads (S1 breached), and a resolution of them on
/// the skeleton's direction field ONLY (the hydrology's D8 is untouched). P0: where F129's 6.05 % sit, by
/// exclusive class, and whether `resolve_flats` (Garbrecht-Martz + the chain's `FlatPerturbation`) saw
/// them. P1: D8-LTD in f64 on `filled` + the Garbrecht-Martz gradient laid into the elevation, per flat,
/// below a tenth of the flat's smallest rise (declared in `f130_predictions.md` before the run).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_p --nocapture
#[test]
#[ignore]
fn f130_p() {
    use std::collections::VecDeque;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::ltd_directions_masked;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 130-P . the flats of the field the skeleton reads  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let cell_m = CELL_KM * 1000.0;
    let dc = C1DrainageConfig::default();
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    let lake_map = d.lake_map.clone();
    drop(d);
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let zm: Vec<f32> = bf.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
    let sink: Vec<bool> = bf.data.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let land = sink.iter().filter(|&&s| !s).count();
    // ── P0 ──
    let flat: Vec<bool> = (0..n).map(|c| !sink[c] && !(0..8).any(|k| zm[nb(c, k)] < zm[c])).collect();
    let nflat = flat.iter().filter(|&&f| f).count();
    let dsea = dist_from(&sink, w, h);
    let is_ocean: Vec<bool> = flow.filled.data.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
    let filled_moved = (0..n).filter(|&c| flow.filled.data[c] != bf.data[c]).count();
    let mut cls = [0usize; 5];
    let (mut in_needs, mut cardinal) = (0usize, 0usize);
    for c in (0..n).filter(|&c| flat[c]) {
        let k = if (0..8).any(|k| bf.data[nb(c, k)] < bf.data[c]) {
            0 // a rounding tie of the metres conversion
        } else if lake_map[c] != 0 {
            1
        } else if bf.data[c] != s1.data[c] {
            2
        } else if dsea[c] <= 20 {
            3
        } else {
            4
        };
        cls[k] += 1;
        if needs[c] {
            in_needs += 1;
        }
        if flow.direction[c] % 2 == 0 {
            cardinal += 1;
        }
    }
    let pc = |v: usize| 100.0 * v as f32 / nflat.max(1) as f32;
    eprintln!(
        "   {nflat} flat land cells ({:.3} % of {land}; Finding 129: 684 893) · `compute_flow` re-filled {filled_moved} cells of the breached field",
        100.0 * nflat as f32 / land as f32
    );
    eprintln!(
        "   exclusive classes: rounding ties of the metres conversion {} ({:.1} %) · pre-incision lake footprints {} ({:.1} %) · cells the breach changed {} ({:.1} %) · coastal (≤ 976 m of the sea) {} ({:.1} %) · other {} ({:.1} %)",
        cls[0], pc(cls[0]), cls[1], pc(cls[1]), cls[2], pc(cls[2]), cls[3], pc(cls[3]), cls[4], pc(cls[4])
    );
    eprintln!(
        "   in `resolve_flats`' flat set (Garbrecht-Martz, on compute_flow's filled): {in_needs} ({:.1} %) · their `compute_flow` pointer CARDINAL {cardinal} ({:.1} %)",
        pc(in_needs),
        pc(cardinal)
    );
    drop((zm, dsea));
    // ── P1 ──
    let f = &flow.filled.data;
    let mut comp = vec![u32::MAX; n];
    let mut stats: Vec<(f64, f64, f64)> = Vec::new(); // (g_min, g_max, gap)
    for c0 in 0..n {
        if !needs[c0] || comp[c0] != u32::MAX {
            continue;
        }
        let id = stats.len() as u32;
        let (mut gmin, mut gmax, mut gap) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY);
        let mut q = VecDeque::from([c0]);
        comp[c0] = id;
        while let Some(c) = q.pop_front() {
            gmin = gmin.min(g[c]);
            gmax = gmax.max(g[c]);
            for k in 0..8 {
                let m = nb(c, k);
                if !is_ocean[m] && f[m] > f[c] {
                    gap = gap.min(f[m] as f64 - f[c] as f64);
                }
                if needs[m] && f[m] == f[c] && comp[m] == u32::MAX {
                    comp[m] = id;
                    q.push_back(m);
                }
            }
        }
        stats.push((gmin, gmax, if gap.is_finite() { gap } else { 1e-7 }));
    }
    let z_ref: Vec<f64> = f.iter().map(|&v| v as f64).collect();
    let mut z_p1 = z_ref.clone();
    for c in 0..n {
        if comp[c] != u32::MAX {
            let (gmin, gmax, gap) = stats[comp[c] as usize];
            z_p1[c] += 0.1 * gap * (g[c] - gmin + 1.0) / (gmax - gmin + 2.0);
        }
    }
    eprintln!("   {} flats (connected equal-elevation sets of `resolve_flats`), {} cells", stats.len(), needs.iter().filter(|&&b| b).count());
    let (d_ref, st_ref) = ltd_directions_masked(&z_ref, &is_ocean, None, &flow.direction, w, h, cell_m);
    let (d_p1, st_p1) = ltd_directions_masked(&z_p1, &is_ocean, None, &flow.direction, w, h, cell_m);
    let land_f = is_ocean.iter().filter(|&&s| !s).count();
    let flat_ref: Vec<bool> = (0..n).map(|c| !is_ocean[c] && !(0..8).any(|k| is_ocean[nb(c, k)] || z_ref[nb(c, k)] < z_ref[c])).collect();
    let nfr = flat_ref.iter().filter(|&&b| b).count();
    let card = |dir: &[u8]| (0..n).filter(|&c| flat_ref[c] && dir[c] % 2 == 0).count();
    eprintln!(
        "   the reference (LTD f64 on filled, no increment): flat fallbacks **{}** ({:.3} % of {land_f} land) · flat cells by the elevation {nfr}",
        st_ref.flats,
        100.0 * st_ref.flats as f32 / land_f as f32
    );
    eprintln!(
        "   P1 (the Garbrecht-Martz gradient in the elevation): flat fallbacks **{}** ({:.4} % of land; C's condition < 0.1 %)",
        st_p1.flats,
        100.0 * st_p1.flats as f32 / land_f as f32
    );
    eprintln!(
        "   on the {nfr} former flats: CARDINAL pointers, reference {} ({:.1} %) → P1 {} ({:.1} %)",
        card(&d_ref),
        100.0 * card(&d_ref) as f32 / nfr.max(1) as f32,
        card(&d_p1),
        100.0 * card(&d_p1) as f32 / nfr.max(1) as f32
    );
    // the control on the non-flat cells: their steepest Tarboton facet, and their LTD pointer
    let facet = |z: &[f64], c: usize| -> Option<usize> {
        const FACETS: [(usize, usize); 8] = [(0, 1), (0, 7), (2, 1), (2, 3), (4, 3), (4, 5), (6, 5), (6, 7)];
        let e0 = z[c];
        let mut best: Option<(f32, usize)> = None;
        for (i, &(kc, kd)) in FACETS.iter().enumerate() {
            let (e1, e2) = (z[nb(c, kc)], z[nb(c, kd)]);
            let (s1, s2) = (((e0 - e1) as f32) / cell_m, ((e1 - e2) as f32) / cell_m);
            let r = s2.atan2(s1);
            let s = if r < 0.0 {
                s1
            } else if r > std::f32::consts::FRAC_PI_4 {
                ((e0 - e2) as f32) / (std::f32::consts::SQRT_2 * cell_m)
            } else {
                (s1 * s1 + s2 * s2).sqrt()
            };
            if s > 0.0 && best.is_none_or(|b| s > b.0) {
                best = Some((s, i));
            }
        }
        best.map(|b| b.1)
    };
    let nonflat: Vec<usize> = (0..n).filter(|&c| !is_ocean[c] && !flat_ref[c]).collect();
    let same_facet = nonflat.par_iter().filter(|&&c| facet(&z_ref, c) == facet(&z_p1, c)).count();
    let same_ptr = nonflat.iter().filter(|&&c| d_ref[c] == d_p1[c]).count();
    let adj = nonflat.iter().filter(|&&c| (0..8).any(|k| comp[nb(c, k)] != u32::MAX)).count();
    eprintln!(
        "   control on the {} non-flat land cells: steepest facet UNCHANGED {same_facet} ({:.4} %) · LTD pointer UNCHANGED {same_ptr} ({:.3} %) · of them adjacent to a flat {adj}",
        nonflat.len(),
        100.0 * same_facet as f64 / nonflat.len() as f64,
        100.0 * same_ptr as f64 / nonflat.len() as f64
    );
    eprintln!("\n==========  end Finding 130-P . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 130-A2 — the RESTRICTED confluence clause on the pipeline's worlds (light pass + A1+B2, as
/// Finding 128-A): the over-dug census, the cols 7 and 14, canyon 13's dam cell and the walk of the
/// témoin's path from 13's floor, the water bodies; beside them Finding 128-A's own columns. mur ↔ mer and
/// bruit ↔ pied ON. The law's part (raised cells, θ laid with its CI) is `f130_a01`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_a2w --nocapture
#[test]
#[ignore]
fn f130_a2w() {
    use std::collections::HashMap;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 130-A2 . the RESTRICTED clause, the worlds (mur ↔ mer and bruit ↔ pied ON)  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let prof = WallProfile { foot_quiet: Some(0.5), ..WallProfile::f124() };
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let at = |x: usize, y: usize| y * w + x;
    let fixed: [(&str, usize); 5] = [
        ("col 7 (4551,3280)", at(4551, 3280)),
        ("col 14 (1805,4580)", at(1805, 4580)),
        ("13's dam cell (2315,4614)", at(2315, 4614)),
        ("canyon 7's floor (4588,3240)", at(4588, 3240)),
        ("canyon 14's floor (1835,4562)", at(1835, 4562)),
    ];
    let recv = |dir: &[u8], c: usize| -> Option<usize> {
        let d = dir[c];
        if d == DIR_NONE {
            return None;
        }
        let nx = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        Some(ny * w + nx)
    };
    struct Read {
        g: GridF32,
        bre: GridF32,
        lake_map: Vec<u32>,
        levels: HashMap<u32, (f32, usize)>,
        dir: Vec<u8>,
        canyons: Vec<(u32, f32, usize)>,
    }
    let b1f = ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), ..base };
    let worlds: [(&str, ValleyConstruction); 6] = [
        ("témoin C2/10 col", base),
        ("B2 → A_c (no clause)", ac),
        ("B2 → A_c + FULL clause", ValleyConstruction { trunk_band: true, ..ac }),
        ("B2 → A_c + RESTRICTED clause", ValleyConstruction { trunk_band: true, trunk_band_downstream: true, ..ac }),
        ("B2 (A_c) + B1 + foot + FULL clause", ValleyConstruction { trunk_band: true, ..b1f }),
        ("B2 (A_c) + B1 + foot + RESTRICTED clause", ValleyConstruction { trunk_band: true, trunk_band_downstream: true, ..b1f }),
    ];
    let near = 2.0 / CELL_KM;
    let mut keep: HashMap<&str, Read> = HashMap::new();
    let mut g_t: Option<GridF32> = None;
    let mut summary = Vec::new();
    for (label, vc) in worlds {
        let t = Instant::now();
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let secs = t.elapsed().as_secs_f64();
        eprintln!("\n────────── {label} · build {secs:.0} s ──────────");
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw));
        let co = coast(&g, &ss, &dwall, near);
        let nn = co.near.iter().filter(|&&b| b).count();
        let rn = nn as f32 / co.l_near_km.max(1e-6);
        let sig = sigma_p50(&g, &ss);
        let r8w = r8_terrain(&g, &cls, 16);
        let r8t = r8_terrain(g_t.as_ref().unwrap_or(&g), &cls, 16);
        eprintln!(
            "   COAST {} spurs · near a coastal wall **{rn:.4} /km** · σ p50 **{sig:.3} m** · R8 terrain {r8w:.4} (témoin on the same classes {r8t:.4})",
            co.mids.len()
        );
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = Vec::new();
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            let dfill = fill_del[floor] - fill_pre[floor];
            if over_dug_depression(dfill, b.rim) {
                eprintln!(
                    "   OVER-DUG body {} · {:.2} km² (domain) · floor ({},{}) · Δfill {dfill:.1} m · class {}",
                    b.id,
                    b.km2,
                    floor % w,
                    floor / w,
                    CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea")
                );
                canyons.push((b.id, b.km2, floor));
            }
        }
        eprintln!("   over-dug bodies: {}", canyons.len());
        drop(fill_del);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage
        }));
        let mut row = format!(
            "{label:<42} over-dug {} {:?} · spurs near {rn:.4}/km · σ {sig:.3} m · R8 terr {r8w:.4}",
            canyons.len(),
            canyons.iter().map(|c| (c.2 % w, c.2 / w)).collect::<Vec<_>>()
        );
        let heights: Vec<String> = fixed.iter().map(|(nm, k)| format!("{nm} {:.1} m", mm(&g, *k))).collect();
        eprintln!("   heights: {}", heights.join(" · "));
        match res {
            Ok(dr) => {
                let segs = &dr.rivers.segments;
                let wall: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k]).collect();
                let wc: Vec<usize> = (0..segs.len())
                    .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2)
                    .collect();
                let teeth = wc
                    .iter()
                    .filter(|&&i| {
                        let p = &segs[i].points;
                        p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32 >= 0.8 * p.len() as f32
                    })
                    .count();
                let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
                eprintln!("   teeth {teeth} · R8 network (chord 8) {:.4} · water bodies (lakes) {} · Finding 38 holds", r8_chords(&pts, 8), dr.lakes.len());
                row += &format!(" · teeth {teeth} · R8 net {:.4} · lakes {}", r8_chords(&pts, 8), dr.lakes.len());
                {
                    let levels = dr
                        .lakes
                        .iter()
                        .map(|l| (l.base.id, (l.level_m, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize)))
                        .collect();
                    keep.insert(label, Read { g: g.clone(), bre: bre.clone(), lake_map: dr.lake_map.clone(), levels, dir: dr.flow.direction.clone(), canyons: canyons.clone() });
                }
            }
            Err(e) => {
                let msg = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
                eprintln!("   ⛔ the HD assembly PANICKED: {}", msg.lines().next().unwrap_or(""));
                row += " · ⛔ Finding 38";
            }
        }
        row += &format!(" · build {secs:.0} s");
        summary.push(row);
        if g_t.is_none() {
            g_t = Some(g);
        }
    }
    // ── canyon 13's dam: the témoin's drainage path from its floor, walked under the clause ──
    eprintln!("\n   ── canyon 13 (floor (2281,4563)): the témoin's path walked in each world ──");
    let f13 = at(2281, 4563);
    let tm = &keep["témoin C2/10 col"];
    let mut labels: Vec<&&str> = keep.keys().filter(|l| !l.starts_with("témoin")).collect();
    labels.sort();
    for label in labels {
        let r = &keep[*label];
        let persists = r.canyons.iter().any(|&(_, _, fl)| {
            let (dx, dy) = ((fl % w) as i64 - 2281, (fl / w) as i64 - 4563);
            dx * dx + dy * dy <= 400
        });
        eprintln!("   {label}: canyon 13 {}", if persists { "PERSISTS" } else { "is GONE" });
        let lid = r.lake_map[f13];
        let level = r.levels.get(&lid).map_or(f32::NAN, |v| v.0);
        let mut c = f13;
        let mut dam: Option<(usize, usize, f32)> = None;
        for step in 0..20_000 {
            let Some(nx) = recv(&tm.dir, c) else { break };
            c = nx;
            if tm.g.data[c] <= SEA {
                break;
            }
            if lid != 0 && r.lake_map[c] != lid && mm(&r.bre, c) > level + 0.01 {
                dam = Some((c, step, mm(&r.bre, c)));
                break;
            }
        }
        match dam {
            Some((c, step, z)) => eprintln!(
                "   {label}: lake {lid} at {level:.1} m · the témoin's path is DAMMED at ({},{}), {step} cells down, at {z:.1} m (témoin {:.1} m)",
                c % w,
                c / w,
                mm(&tm.bre, c)
            ),
            None => eprintln!("   {label}: lake {lid} (level {level:.1} m) · no dam on the témoin's path above the lake level"),
        }
    }
    eprintln!("\n   ── summary ──");
    for s in summary {
        eprintln!("   {s}");
    }
    eprintln!("\n==========  end Finding 130-A2 (worlds) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 130-E — the one-sided construction, MEASURED only: the witness's trunk ≥ 10 km² cells the
/// construction leaves UNCARVED with the field below its law (Finding 129's 32.1 %). After the drainage of the
/// as-built witness (stage ii): under a lake? In a low plain (≤ 100 m above the sea and steepest slope
/// ≤ 0.01)? Where along the trunk (the D8 distance to the base; coastal = within 2 km of the sea)? Every
/// number beside the same number for ALL the trunk ≥ 10 km² cells (the base rate). Nothing is built.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f130_e --nocapture
#[test]
#[ignore]
fn f130_e() {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 130-E . the trunk cells below their law (measured only)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let sk = skeleton(&s1, &base, &ss, DOMAIN_KM);
    let law: Vec<f32> = (0..n).map(|k| if sk.chi_m[k].is_finite() { sk.floor_m(k, base.age_k) } else { f32::NAN }).collect();
    let (bt, mt) = carve(&s1, &sk, &base, &ss);
    let zt = metres(&bt);
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    let pop: Vec<usize> = trunk10.iter().copied().filter(|&k| !mt.carved[k] && zs1[k] < law[k]).collect();
    // stage (ii): the drainage of the as-built witness
    let d = c1_drainage_windowed(&bt, None, &dcfg(), &ss, DOMAIN_KM);
    let lake = d.lake_map.clone();
    drop((d, bt));
    let sea_m = c1_altitude_norm_to_metres(SEA, &ss);
    let cell_m = CELL_KM * 1000.0;
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let steep = |c: usize| -> f32 {
        (0..8).map(|k| (zt[c] - zt[nb(c, k)]) / (cell_m * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 })).fold(0f32, f32::max)
    };
    let sea: Vec<bool> = (0..n).map(|k| s1.data[k] <= SEA).collect();
    let dsea = dist_from(&sea, w, h);
    let coast_cells = (2.0 / CELL_KM).round() as u16;
    // the D8 distance (m) of each land cell to its base, on the skeleton's own tree
    let mut dbase = vec![f32::NAN; n];
    let mut path = Vec::new();
    for s in 0..n {
        if sea[s] || !dbase[s].is_nan() {
            continue;
        }
        path.clear();
        let mut c = s;
        let mut acc;
        loop {
            if !dbase[c].is_nan() {
                acc = dbase[c];
                break;
            }
            path.push(c);
            let dd = sk.direction[c];
            if dd == DIR_NONE {
                acc = 0.0;
                path.pop();
                dbase[c] = 0.0;
                break;
            }
            let r = nb(c, dd as usize);
            if sea[r] {
                acc = 0.0;
                break;
            }
            c = r;
            if path.len() > n {
                panic!("cycle");
            }
        }
        for &p in path.iter().rev() {
            let dd = sk.direction[p] as usize;
            acc += cell_m * if dd % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
            dbase[p] = acc;
        }
    }
    let describe = |label: &str, cells: &[usize]| {
        let m = cells.len().max(1) as f32;
        let under = cells.iter().filter(|&&k| lake[k] != 0).count();
        let low = cells.iter().filter(|&&k| zt[k] - sea_m <= 100.0).count();
        let plain = cells.iter().filter(|&&k| zt[k] - sea_m <= 100.0 && steep(k) <= 0.01).count();
        let coastal = cells.iter().filter(|&&k| dsea[k] <= coast_cells).count();
        let mut db: Vec<f32> = cells.iter().map(|&k| dbase[k] / 1000.0).collect();
        db.sort_by(f32::total_cmp);
        let mut alt: Vec<f32> = cells.iter().map(|&k| zt[k] - sea_m).collect();
        alt.sort_by(f32::total_cmp);
        let mut ar: Vec<f32> = cells.iter().map(|&k| sk.area_km2[k]).collect();
        ar.sort_by(f32::total_cmp);
        let q = |v: &[f32], p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p) as usize] };
        eprintln!(
            "   {label:<44} {} cells · under a lake after drainage **{under}** ({:.1} %) · ≤ 100 m above the sea {low} ({:.1} %) · LOW PLAIN (and slope ≤ 0.01) **{plain}** ({:.1} %) · coastal (≤ 2 km) {coastal} ({:.1} %)",
            cells.len(),
            100.0 * under as f32 / m,
            100.0 * low as f32 / m,
            100.0 * plain as f32 / m,
            100.0 * coastal as f32 / m
        );
        eprintln!(
            "      D8 distance to the base p10 / p50 / p90 {:.1} / {:.1} / {:.1} km (domain) · altitude above the sea p10 / p50 / p90 {:.0} / {:.0} / {:.0} m · drained area p10 / p50 / p90 {:.0} / {:.0} / {:.0} km²",
            q(&db, 0.1), q(&db, 0.5), q(&db, 0.9), q(&alt, 0.1), q(&alt, 0.5), q(&alt, 0.9), q(&ar, 0.1), q(&ar, 0.5), q(&ar, 0.9)
        );
    };
    describe("ALL trunk ≥ 10 km² cells (the base rate)", &trunk10);
    describe("the UNCARVED cells below their law", &pop);
    let mut gap: Vec<f32> = pop.iter().map(|&k| law[k] - zs1[k]).collect();
    gap.sort_by(f32::total_cmp);
    let q = |v: &[f32], p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p) as usize] };
    eprintln!("   the law above the field on them: p10 / p50 / p90 {:.1} / {:.1} / {:.1} m", q(&gap, 0.1), q(&gap, 0.5), q(&gap, 0.9));
    eprintln!("\n==========  end Finding 130-E . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

// ════════════════════════════════ Finding 131 ════════════════════════════════

/// ADR Finding 131-B — D8 and D8-LTD on three ANALYTIC surfaces whose exact answer we compute (no number of
/// [O03] is used): a plane, a divergent cone (a = r/2 exact; Paik 2008's surface), a curved valley
/// z = a x² + b y (exact lines x = x₀·e^{2a(y−y₀)/b}, exact a in closed form). Declared in
/// `f131_predictions.md` before the run. `F131_B=d8` runs the references only (rule 15); `F131_B=all` adds LTD.
///
/// Run: F131_B=d8 cargo test -p ymir-core --release --test f126_coast -- --ignored f131_b --nocapture
#[test]
#[ignore]
fn f131_b() {
    use ymir_core::tectonics_c1::valley_construction::{accumulate_cells, ltd_directions_masked};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let with_ltd = std::env::var("F131_B").map(|v| v == "all").unwrap_or(false);
    eprintln!("\n==========  Finding 131-B . D8{} on three analytic surfaces  ==========", if with_ltd { " and D8-LTD" } else { " (the references)" });
    struct Surf {
        name: &'static str,
        p: usize,
        z: Vec<f64>,
        sink: Vec<bool>,
        active: Vec<bool>,
        // the cells the errors are read on, with their centre (x, y) in the surface's coordinates
        eval: Vec<(usize, f64, f64)>,
    }
    // the D8 steepest descent on the active cells
    let d8_of = |s: &Surf| -> Vec<u8> {
        let p = s.p;
        (0..p * p)
            .map(|k| {
                if !s.active[k] || s.sink[k] {
                    return DIR_NONE;
                }
                let (x, y) = ((k % p) as i32, (k / p) as i32);
                let mut best = (DIR_NONE, 0f64);
                for kk in 0..8u8 {
                    let (nx, ny) = (x + D8_DX[kk as usize], y + D8_DY[kk as usize]);
                    if nx < 0 || ny < 0 || nx >= p as i32 || ny >= p as i32 {
                        continue;
                    }
                    let m = ny as usize * p + nx as usize;
                    if !s.active[m] {
                        continue;
                    }
                    let d = (s.z[k] - s.z[m]) / if kk % 2 == 1 { std::f64::consts::SQRT_2 } else { 1.0 };
                    if d > best.1 {
                        best = (kk, d);
                    }
                }
                best.0
            })
            .collect()
    };
    let step = |p: usize, k: usize, d: u8| -> usize { ((k / p) as i32 + D8_DY[d as usize]) as usize * p + ((k % p) as i32 + D8_DX[d as usize]) as usize };
    let q = |v: &mut Vec<f64>, pp: f64| -> f64 {
        if v.is_empty() {
            return f64::NAN;
        }
        v.sort_by(f64::total_cmp);
        v[((v.len() - 1) as f64 * pp) as usize]
    };
    let angle = |a: (f64, f64), b: (f64, f64)| -> f64 { (a.0 * b.1 - a.1 * b.0).atan2(a.0 * b.0 + a.1 * b.1).abs().to_degrees() };
    // ── the plane: 60², z = −(y + x/4), open on its two downslope edges ──
    let plane = {
        let (n, p) = (60usize, 62usize);
        let (mut z, mut sink, mut active) = (vec![0f64; p * p], vec![false; p * p], vec![false; p * p]);
        let mut eval = Vec::new();
        for k in 0..p * p {
            let (ix, iy) = ((k % p) as i64 - 1, (k / p) as i64 - 1);
            let (x, y) = (ix as f64 + 0.5, iy as f64 + 0.5);
            z[k] = -(y + x / 4.0);
            let inside = (0..n as i64).contains(&ix) && (0..n as i64).contains(&iy);
            let out_edge = (ix == n as i64 && iy >= 0) || (iy == n as i64 && ix >= 0);
            active[k] = inside || out_edge;
            sink[k] = out_edge;
            if inside {
                eval.push((k, x, y));
            }
        }
        Surf { name: "plane 60² (14.04° off the axis)", p, z, sink, active, eval }
    };
    let pd = {
        let g = (0.25f64, 1.0f64);
        let nn = (g.0 * g.0 + g.1 * g.1).sqrt();
        (g.0 / nn, g.1 / nn)
    };
    // ── the cone: 101², apex at the centre of (50, 50), z = −r, r ≤ 48 ──
    let cone = {
        let p = 101usize;
        let (mut z, mut sink, active) = (vec![0f64; p * p], vec![false; p * p], vec![true; p * p]);
        let mut eval = Vec::new();
        for k in 0..p * p {
            let (x, y) = ((k % p) as f64 + 0.5 - 50.5, (k / p) as f64 + 0.5 - 50.5);
            let r = (x * x + y * y).sqrt();
            z[k] = -r;
            sink[k] = r > 48.0;
            if (2.0..=47.0).contains(&r) {
                eval.push((k, x, y));
            }
        }
        Surf { name: "cone r ≤ 48 (divergent)", p, z, sink, active, eval }
    };
    // ── the curve: z = a x² + b y, x in [−30, 30], y in [0, 60], outlet below y = 0 ──
    let (ca, cb) = (0.02f64, 1.0f64);
    let curve = {
        let (n, p) = (60usize, 62usize);
        let (mut z, mut sink, mut active) = (vec![0f64; p * p], vec![false; p * p], vec![false; p * p]);
        let mut eval = Vec::new();
        for k in 0..p * p {
            let (ix, iy) = ((k % p) as i64 - 1, (k / p) as i64 - 1);
            let (x, y) = (ix as f64 + 0.5 - 30.0, iy as f64 + 0.5);
            z[k] = ca * x * x + cb * y;
            let inside = (0..n as i64).contains(&ix) && (0..n as i64).contains(&iy);
            let outlet = iy == -1 && (0..n as i64).contains(&ix);
            active[k] = inside || outlet;
            sink[k] = outlet;
            if inside {
                eval.push((k, x, y));
            }
        }
        Surf { name: "curve z = 0.02 x² + y (60²)", p, z, sink, active, eval }
    };
    let curve_a = |x0: f64, y0: f64| -> f64 {
        let ymax = if x0 == 0.0 { 60.0 } else { (y0 + cb / (2.0 * ca) * (30.0 / x0.abs()).ln()).min(60.0) };
        let g = ((2.0 * ca * x0).powi(2) + cb * cb).sqrt();
        cb / (2.0 * ca) * ((2.0 * ca * (ymax - y0) / cb).exp() - 1.0) * g / cb
    };
    // the closed form against a narrow tube traced at high resolution (5 points)
    {
        let mut out = Vec::new();
        for &(x0, y0) in &[(0.3f64, 5.0f64), (5.0, 10.0), (-12.0, 30.0), (20.0, 2.0), (-28.0, 50.0)] {
            let g0 = ((2.0 * ca * x0).powi(2) + cb * cb).sqrt();
            let h = 1e-4;
            // two lines through p ± h·(the unit normal to the flow), traced upslope (dy > 0) to the domain's edge
            let nrm = (cb / g0, -2.0 * ca * x0 / g0);
            let trace = |sx: f64, sy: f64| -> Vec<(f64, f64)> {
                let (mut x, mut y) = (sx, sy);
                let mut pts = vec![(x, y)];
                let dy = 1e-3;
                while y < 60.0 && x.abs() < 30.0 {
                    x += 2.0 * ca * x / cb * dy;
                    y += dy;
                    pts.push((x, y));
                }
                pts
            };
            // two lines through (x0 ± h, y0), the SAME y, so each strip compares them at one y; the tube's
            // width across the flow at p is 2h · b / |grad z|
            let _ = nrm;
            let la = trace(x0 + h, y0);
            let lb = trace(x0 - h, y0);
            // the area between them, by horizontal strips (both parametrised by y)
            let m = la.len().min(lb.len());
            let mut area = 0.0;
            for i in 1..m {
                let (w0, w1) = ((la[i - 1].0 - lb[i - 1].0).abs(), (la[i].0 - lb[i].0).abs());
                area += 0.5 * (w0 + w1) * (la[i].1 - la[i - 1].1);
            }
            out.push(format!("({x0}, {y0}): closed {:.3} · tube {:.3}", curve_a(x0, y0), area / (2.0 * h * cb / g0)));
        }
        eprintln!("   the curve's closed-form a against a traced tube: {}", out.join(" · "));
    }
    let surfaces = [plane, cone, curve];
    for s in &surfaces {
        let p = s.p;
        let d8 = d8_of(s);
        let mut methods: Vec<(&str, Vec<u8>)> = vec![("D8", d8.clone())];
        if with_ltd {
            let (ltd, st) = ltd_directions_masked(&s.z, &s.sink, Some(&s.active), &d8, p, p, 1.0);
            eprintln!("   ({}: LTD flat fallbacks {})", s.name, st.flats);
            methods.push(("D8-LTD", ltd));
        }
        eprintln!("\n   ── {} ──", s.name);
        for (mname, dir) in &methods {
            let acc = accumulate_cells(dir, &s.sink, p, p);
            // the exact a at the downslope edge, and the exact descent at a point
            let (mut e, mut ang): (Vec<f64>, [Vec<f64>; 4]) = (Vec::new(), [Vec::new(), Vec::new(), Vec::new(), Vec::new()]);
            for &(k, x, y) in &s.eval {
                let a_exact = if s.name.starts_with("plane") {
                    let (ex, ey) = (x + 0.5 * pd.0, y + 0.5 * pd.1);
                    (ex / pd.0).min(ey / pd.1) // the upslope ray from (ex, ey) along −d to x = 0 or y = 0
                } else if s.name.starts_with("cone") {
                    ((x * x + y * y).sqrt() + 0.5) / 2.0
                } else {
                    let g = ((2.0 * ca * x).powi(2) + cb * cb).sqrt();
                    curve_a(x - 0.5 * 2.0 * ca * x / g, y - 0.5 * cb / g)
                };
                e.push((acc[k] as f64 - a_exact) / a_exact);
                // the chords
                for (ci, &kk) in [1usize, 8, 16, 32].iter().enumerate() {
                    let mut c = k;
                    let mut ok = true;
                    for _ in 0..kk {
                        let d = dir[c];
                        if d == DIR_NONE {
                            ok = false;
                            break;
                        }
                        c = step(p, c, d);
                        if s.sink[c] || !s.active[c] {
                            ok = false;
                            break;
                        }
                    }
                    if !ok {
                        continue;
                    }
                    let (gx, gy) = (((c % p) as f64 - (k % p) as f64), ((c / p) as f64 - (k / p) as f64));
                    let ex = if s.name.starts_with("plane") {
                        pd
                    } else if s.name.starts_with("cone") {
                        (x, y)
                    } else {
                        let yk = y + gy;
                        (x * (2.0 * ca * (yk - y) / cb).exp() - x, yk - y)
                    };
                    ang[ci].push(angle((gx, gy), ex));
                }
            }
            let m = e.len() as f64;
            let mae = e.iter().map(|v| v.abs()).sum::<f64>() / m;
            let rmse = (e.iter().map(|v| v * v).sum::<f64>() / m).sqrt();
            let chords: Vec<String> = [1, 8, 16, 32]
                .iter()
                .enumerate()
                .map(|(i, kk)| {
                    let mean = ang[i].iter().sum::<f64>() / ang[i].len().max(1) as f64;
                    let n_ = ang[i].len();
                    format!("{kk}: {mean:.2}° (p90 {:.1}°, {n_})", q(&mut ang[i], 0.9))
                })
                .collect();
            eprintln!("      {mname:<7} specific area: MAE **{mae:.3}** · RMSE **{rmse:.3}** ({} cells) · path angle by chord {}", e.len(), chords.join(" · "));
            if s.name.starts_with("cone") {
                let mut donors = vec![0u8; p * p];
                for k in 0..p * p {
                    if dir[k] != DIR_NONE {
                        donors[step(p, k, dir[k])] += 1;
                    }
                }
                let conf = (0..p * p).filter(|&k| !s.sink[k] && donors[k] >= 2).count();
                let mut sect = [0f64; 16];
                for k in 0..p * p {
                    let (x, y) = ((k % p) as f64 + 0.5 - 50.5, (k / p) as f64 + 0.5 - 50.5);
                    let r = (x * x + y * y).sqrt();
                    if r > 47.0 && r <= 48.0 {
                        let a = (y.atan2(x) + std::f64::consts::PI) / (2.0 * std::f64::consts::PI);
                        sect[((a * 16.0) as usize).min(15)] += acc[k] as f64;
                    }
                }
                let mean = sect.iter().sum::<f64>() / 16.0;
                let cv = (sect.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / 16.0).sqrt() / mean;
                // the same statistic on the exact a of the same ring cells: the ring's own discretisation
                let mut ex = [0f64; 16];
                for k in 0..p * p {
                    let (x, y) = ((k % p) as f64 + 0.5 - 50.5, (k / p) as f64 + 0.5 - 50.5);
                    let r = (x * x + y * y).sqrt();
                    if r > 47.0 && r <= 48.0 {
                        let a = (y.atan2(x) + std::f64::consts::PI) / (2.0 * std::f64::consts::PI);
                        ex[((a * 16.0) as usize).min(15)] += (r + 0.5) / 2.0;
                    }
                }
                let me = ex.iter().sum::<f64>() / 16.0;
                let cve = (ex.iter().map(|v| (v - me).powi(2)).sum::<f64>() / 16.0).sqrt() / me;
                eprintln!("               FALSE confluences (≥ 2 donors) **{conf}** · sector asymmetry (CV over 16 sectors of the outer ring) **{cv:.3}** (the exact a on the same ring: {cve:.3})");
            }
        }
    }
    eprintln!("\n==========  end Finding 131-B  ==========\n");
}

/// ADR Finding 131-P — do the lakes' flats touch what the construction lays? The pointers of the flat cells
/// of the pre-drainage lake footprints are replaced (identity: the hook's control; LTD's on the
/// Garbrecht-Martz elevation; an arbitrary acyclic choice), the skeleton rebuilt through
/// `skeleton_patched`, and the construction compared OUTSIDE the lake footprints. Declared in
/// `f131_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f131_p --nocapture
#[test]
#[ignore]
fn f131_p() {
    use std::collections::VecDeque;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{ltd_directions_masked, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 131-P . do the lakes' flats touch what is constructed?  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let cell_m = CELL_KM * 1000.0;
    let dc = C1DrainageConfig::default();
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    // the flats, their Garbrecht-Martz field and LTD's pointers on it (Finding 130-P1), computed once
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    let lake = d.lake_map.clone();
    drop(d);
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let f = flow.filled.data.clone();
    let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
    let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
    let nt = target.iter().filter(|&&b| b).count();
    let mut comp = vec![u32::MAX; n];
    let mut stats: Vec<(f64, f64, f64)> = Vec::new();
    for c0 in 0..n {
        if !needs[c0] || comp[c0] != u32::MAX {
            continue;
        }
        let id = stats.len() as u32;
        let (mut gmin, mut gmax, mut gap) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY);
        let mut q = VecDeque::from([c0]);
        comp[c0] = id;
        while let Some(c) = q.pop_front() {
            gmin = gmin.min(g[c]);
            gmax = gmax.max(g[c]);
            for k in 0..8 {
                let m = nb(c, k);
                if !is_ocean[m] && f[m] > f[c] {
                    gap = gap.min(f[m] as f64 - f[c] as f64);
                }
                if needs[m] && f[m] == f[c] && comp[m] == u32::MAX {
                    comp[m] = id;
                    q.push_back(m);
                }
            }
        }
        stats.push((gmin, gmax, if gap.is_finite() { gap } else { 1e-7 }));
    }
    let mut z = f.iter().map(|&v| v as f64).collect::<Vec<f64>>();
    for c in 0..n {
        if comp[c] != u32::MAX {
            let (gmin, gmax, gap) = stats[comp[c] as usize];
            z[c] += 0.1 * gap * (g[c] - gmin + 1.0) / (gmax - gmin + 2.0);
        }
    }
    let (ltd, _) = ltd_directions_masked(&z, &is_ocean, None, &flow.direction, w, h, cell_m);
    drop(z);
    // the arbitrary acyclic pointer: a seeded choice among the equal neighbours with a smaller flat_grad
    let arb: Vec<u8> = (0..n)
        .map(|c| {
            if !target[c] {
                return flow.direction[c];
            }
            let cand: Vec<u8> = (0..8u8)
                .filter(|&k| {
                    let m = nb(c, k as usize);
                    is_ocean[m] || f[m] < f[c] || (f[m] == f[c] && needs[m] && g[m] < g[c])
                })
                .collect();
            if cand.is_empty() {
                return flow.direction[c];
            }
            let hsh = (c as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(17) ^ 0xA5A5_5A5A;
            cand[(hsh % cand.len() as u64) as usize]
        })
        .collect();
    let ch_ltd = (0..n).filter(|&c| target[c] && ltd[c] != flow.direction[c]).count();
    let ch_arb = (0..n).filter(|&c| target[c] && arb[c] != flow.direction[c]).count();
    eprintln!("   {nt} flat lake cells patched · pointers changed there: LTD {ch_ltd} · arbitrary {ch_arb}");
    drop((flow, f, g, needs, comp));
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let hash_out = |g: &GridF32| -> u64 {
        let mut hh = 0xcbf2_9ce4_8422_2325u64;
        for k in 0..n {
            if lake[k] == 0 {
                hh ^= g.data[k].to_bits() as u64;
                hh = hh.wrapping_mul(0x100_0000_01b3);
            }
        }
        hh
    };
    for (label, vc) in [("témoin C2/10 col", base), ("B2 → A_c", ac)] {
        let t = Instant::now();
        let (b0, _) = carve(&s1, &skeleton(&s1, &vc, &ss, DOMAIN_KM), &vc, &ss);
        let variants: [(&str, Box<dyn Fn(&GridF32, &[u32], &mut Vec<u8>)>); 3] = [
            ("identity (the control)", Box::new(|_, _, _| {})),
            ("LTD on the lakes' flats", Box::new(|_, _, d: &mut Vec<u8>| {
                for c in 0..n {
                    if target[c] {
                        d[c] = ltd[c];
                    }
                }
            })),
            ("arbitrary on the lakes' flats", Box::new(|_, _, d: &mut Vec<u8>| {
                for c in 0..n {
                    if target[c] {
                        d[c] = arb[c];
                    }
                }
            })),
        ];
        let mut reference: Option<GridF32> = None;
        eprintln!("\n   ── {label} ── (no patch: outside-hash {:016x})", hash_out(&b0));
        for (vname, patch) in variants.iter() {
            let sk = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(patch.as_ref()));
            let (bp, _) = carve(&s1, &sk, &vc, &ss);
            let against = reference.as_ref().unwrap_or(&b0);
            let (mut out_diff, mut in_diff) = (0usize, 0usize);
            let mut near = Vec::new();
            for k in 0..n {
                if bp.data[k] != against.data[k] {
                    if lake[k] == 0 {
                        out_diff += 1;
                        near.push(k);
                    } else {
                        in_diff += 1;
                    }
                }
            }
            // where the outside changes sit: the distance to the nearest lake cell
            let lk: Vec<bool> = lake.iter().map(|&l| l != 0).collect();
            let dl = dist_from(&lk, w, h);
            let mut dd: Vec<f32> = near.iter().map(|&k| dl[k] as f32 * cell_m).collect();
            dd.sort_by(f32::total_cmp);
            let qq = |p: f32| if dd.is_empty() { f32::NAN } else { dd[((dd.len() - 1) as f32 * p) as usize] };
            eprintln!(
                "   {vname:<32} against {}: cells changed OUTSIDE the lake footprints **{out_diff}** · inside {in_diff} · outside-hash {:016x} · distance of the outside changes to a lake p50 / p90 / max {:.0} / {:.0} / {:.0} m",
                if reference.is_some() { "the identity" } else { "no patch" },
                hash_out(&bp),
                qq(0.5),
                qq(0.9),
                qq(1.0)
            );
            if reference.is_none() {
                reference = Some(bp);
            }
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n==========  end Finding 131-P . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 131-A — the CONCORDANT confluence on `carve(S1)` (the pipeline's own call): the raised cells and
/// θ laid with its CI against Finding 130's references, the replat it makes at the junctions, the tributaries'
/// long profiles and the source cones (Finding 121's trap). Declared in `f131_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f131_a_law --nocapture
#[test]
#[ignore]
fn f131_a_law() {
    use ymir_core::tectonics_c1::valley_construction::carve_diag;
    use ymir_core::terrain::flow::{D8_DX, D8_DY};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 131-A . the concordant confluence, on the trunks (carve on S1)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let cell_m = CELL_KM * 1000.0;
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let prof = WallProfile { foot_quiet: Some(0.5), ..WallProfile::f124() };
    let sk = skeleton(&s1, &base, &ss, DOMAIN_KM);
    let law: Vec<f32> = (0..n).map(|k| if sk.chi_m[k].is_finite() { sk.floor_m(k, base.age_k) } else { f32::NAN }).collect();
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let sk_ac = skeleton(&s1, &ac, &ss, DOMAIN_KM);
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let flat = |z: &[f32], c: usize| !(0..8).any(|k| z[nb(c, k)] < z[c]);
    let junction_on = |small: u32, big: u32| -> Option<u32> {
        let mut l = small;
        for _ in 0..sk_ac.polylines.len() {
            match sk_ac.line_parent[l as usize] {
                Some((p, j)) if p == big => return Some(j),
                Some((p, _)) => l = p,
                None => return None,
            }
        }
        None
    };
    let tor = |ax: f32, ay: f32, bx: f32, by: f32| -> f32 {
        let mut dx = (ax - bx).rem_euclid(w as f32);
        if dx > w as f32 / 2.0 {
            dx -= w as f32;
        }
        let mut dy = (ay - by).rem_euclid(h as f32);
        if dy > h as f32 / 2.0 {
            dy -= h as f32;
        }
        (dx * dx + dy * dy).sqrt()
    };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let mut verdict = Vec::new();
    for (tag, world, floor_ref) in [("B2 → A_c", ac, 86usize), ("B2 (A_c) + B1 + foot", ValleyConstruction { wall_profile: Some(prof), ..ac }, 66)] {
        eprintln!("\n   ── {tag} ──");
        let raised = |zb: &[f32]| trunk10.iter().filter(|&&k| zb[k] - law[k] > 1.0).count();
        let (b0, m0) = carve(&s1, &sk_ac, &world, &ss);
        let z0 = metres(&b0);
        drop(b0);
        let (t0_, lo0, hi0, l0) = theta_laid_ci(&z0, &m0.carved, &sk);
        eprintln!("   no clause (the reference)          raised > 1 m **{}** · θ laid **{t0_:.3}** [{lo0:.3}, {hi0:.3}] ({l0})", raised(&z0));
        let (b1, m1) = carve(&s1, &sk_ac, &ValleyConstruction { trunk_band: true, ..world }, &ss);
        let z1 = metres(&b1);
        drop(b1);
        let (t1, lo1, hi1, l1) = theta_laid_ci(&z1, &m1.carved, &sk);
        eprintln!("   the full clause (the reference)    raised > 1 m **{}** · θ laid **{t1:.3}** [{lo1:.3}, {hi1:.3}] ({l1})", raised(&z1));
        drop((z1, m1));
        let (bc, mc, dg) = carve_diag(&s1, &sk_ac, &ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..world }, &ss);
        let zc = metres(&bc);
        drop(bc);
        let (tc, loc, hic, lc) = theta_laid_ci(&zc, &mc.carved, &sk);
        let rc = raised(&zc);
        eprintln!("   the CONCORDANT clause              raised > 1 m **{rc}** · θ laid **{tc:.3}** [{loc:.3}, {hic:.3}] ({lc})");
        let moved_other = (0..n).filter(|&k| dg.banded.as_ref().unwrap()[k] == u32::MAX && zc[k] != z0[k]).count();
        // the clause's cells, by junction
        let banded = dg.banded.as_ref().expect("the clause is on");
        let mut by_j: std::collections::HashMap<(u32, u32), (Vec<usize>, u32)> = std::collections::HashMap::new();
        for k in 0..n {
            let s = banded[k];
            if s == u32::MAX {
                continue;
            }
            let big = dg.line_of[s as usize];
            let owner = dg.line_of[dg.who[k] as usize];
            let j = junction_on(owner, big).expect("a concordant cell has a junction");
            let e = by_j.entry((big, j)).or_insert((Vec::new(), owner));
            e.0.push(k);
        }
        let ncells: usize = by_j.values().map(|v| v.0.len()).sum();
        let (mut rep_c, mut rep_0) = (0usize, 0usize);
        let (mut len_w, mut flat_cells_c, mut flat_cells_0) = (Vec::new(), 0usize, 0usize);
        let (mut sl_c, mut sl_0) = (Vec::new(), Vec::new());
        for (&(big, j), (cells, owner)) in &by_j {
            let line = &sk_ac.polylines[big as usize];
            let (jx, jy, _, jhw) = line[j as usize];
            let wj = 2.0 * jhw;
            let fc: Vec<usize> = cells.iter().copied().filter(|&k| flat(&zc, k)).collect();
            let f0 = cells.iter().filter(|&&k| flat(&z0, k)).count();
            flat_cells_c += fc.len();
            flat_cells_0 += f0;
            if !fc.is_empty() {
                rep_c += 1;
                let far = fc.iter().map(|&k| tor((k % w) as f32 + 0.5, (k / w) as f32 + 0.5, jx, jy) * cell_m).fold(0f32, f32::max);
                len_w.push(far / wj);
            }
            if f0 > 0 {
                rep_0 += 1;
            }
            // the tributary's floor over its last 2 W before the junction
            let tl = &sk_ac.polylines[*owner as usize];
            let mut arc = 0f32;
            let mut i = tl.len() - 1;
            while i > 0 && arc < 2.0 * wj / cell_m {
                arc += ((tl[i].0 - tl[i - 1].0).powi(2) + (tl[i].1 - tl[i - 1].1).powi(2)).sqrt();
                i -= 1;
            }
            let cell = |p: (f32, f32, f32, f32)| (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize;
            let (a, b) = (cell(tl[i]), cell(tl[tl.len() - 1]));
            if arc > 0.0 {
                sl_c.push((zc[a] - zc[b]) / (arc * cell_m));
                sl_0.push((z0[a] - z0[b]) / (arc * cell_m));
            }
        }
        let nj = by_j.len().max(1);
        let flat_share = |v: &[f32]| v.iter().filter(|&&s| s <= 1e-4).count() as f32 / v.len().max(1) as f32 * 100.0;
        eprintln!(
            "   the clause's cells {ncells} at {} junctions · cells changed OUTSIDE them {moved_other} (must be 0)",
            by_j.len()
        );
        eprintln!(
            "   REPLAT: junctions with >= 1 flat clause cell **{rep_c}** ({:.1} %) against the same cells without the clause {rep_0} ({:.1} %) · flat clause cells {flat_cells_c} against {flat_cells_0} · replat length p50 / p90 / max {:.2} / {:.2} / {:.2} W",
            100.0 * rep_c as f32 / nj as f32,
            100.0 * rep_0 as f32 / nj as f32,
            q(&mut len_w.clone(), 0.5),
            q(&mut len_w.clone(), 0.9),
            q(&mut len_w, 1.0)
        );
        eprintln!(
            "   the tributaries' floor over their last 2 W: concordant p10 / p50 / p90 {:.4} / {:.4} / {:.4} (≤ 1e-4: {:.1} %) · without the clause {:.4} / {:.4} / {:.4} (≤ 1e-4: {:.1} %)",
            q(&mut sl_c.clone(), 0.1),
            q(&mut sl_c.clone(), 0.5),
            q(&mut sl_c.clone(), 0.9),
            flat_share(&sl_c),
            q(&mut sl_0.clone(), 0.1),
            q(&mut sl_0.clone(), 0.5),
            q(&mut sl_0.clone(), 0.9),
            flat_share(&sl_0)
        );
        // the source cones (Finding 129's census, concordant against none)
        let (mut sb, mut sa, mut touched, mut rev_by, mut flat_by) = (Vec::new(), Vec::new(), 0usize, 0usize, 0usize);
        for l in &sk_ac.polylines {
            let cells: Vec<usize> = l
                .iter()
                .map(|p| (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize)
                .fold(Vec::new(), |mut a, c| {
                    if a.last() != Some(&c) {
                        a.push(c);
                    }
                    a
                });
            if cells.len() < 11 || !cells.iter().take(4).any(|&c| zc[c] != z0[c]) {
                continue;
            }
            touched += 1;
            let (b_, a_) = ((z0[cells[0]] - z0[cells[10]]) / (10.0 * cell_m), (zc[cells[0]] - zc[cells[10]]) / (10.0 * cell_m));
            sb.push(b_);
            sa.push(a_);
            if b_ >= 0.0 && a_ < 0.0 {
                rev_by += 1;
            } else if b_ > 0.0 && a_ >= 0.0 && a_ < 0.5 * b_ {
                flat_by += 1;
            }
        }
        let mean = |v: &[f32]| v.iter().sum::<f32>() / v.len().max(1) as f32;
        eprintln!(
            "   source cones touched {touched} of {} lines · head → 10th-cell slope mean {:.4} → {:.4} · p50 {:.4} → {:.4} · reversed BY the clause {rev_by} · halved {flat_by}",
            sk_ac.polylines.len(),
            mean(&sb),
            mean(&sa),
            q(&mut sb.clone(), 0.5),
            q(&mut sa.clone(), 0.5)
        );
        let ok_raised = (rc as i64 - floor_ref as i64).abs() <= 20;
        let ok_theta = tc >= lo0 && tc <= hi0;
        let ok_replat = rep_c <= rep_0;
        verdict.push(format!(
            "{tag}: raised at the floor {} ({rc} vs {floor_ref} ± 20) · θ in the no-clause CI {} ({tc:.3} in [{lo0:.3}, {hi0:.3}]) · replat not up {} ({rep_c} vs {rep_0} junctions)",
            if ok_raised { "YES" } else { "NO" },
            if ok_theta { "YES" } else { "NO" },
            if ok_replat { "YES" } else { "NO" }
        ));
    }
    eprintln!("\n   ── the round's criteria on the construction (canyons and R8 around the junctions: f131_a_worlds) ──");
    for v in verdict {
        eprintln!("   {v}");
    }
    eprintln!("\n==========  end Finding 131-A (law) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 131-A — the CONCORDANT confluence clause on the pipeline's worlds, with R8 terrain and network within 2 W of its junctions; (light pass + A1+B2, as
/// Finding 128-A): the over-dug census, the cols 7 and 14, canyon 13's dam cell and the walk of the
/// témoin's path from 13's floor, the water bodies; beside them Finding 128-A's own columns. mur ↔ mer and
/// bruit ↔ pied ON. The law's part (raised cells, θ laid with its CI) is `f130_a01`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f131_a_worlds --nocapture
#[test]
#[ignore]
fn f131_a_worlds() {
    use std::collections::HashMap;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 131-A . the CONCORDANT clause, the worlds (mur ↔ mer and bruit ↔ pied ON)  ==========");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let prof = WallProfile { foot_quiet: Some(0.5), ..WallProfile::f124() };
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let at = |x: usize, y: usize| y * w + x;
    let fixed: [(&str, usize); 5] = [
        ("col 7 (4551,3280)", at(4551, 3280)),
        ("col 14 (1805,4580)", at(1805, 4580)),
        ("13's dam cell (2315,4614)", at(2315, 4614)),
        ("canyon 7's floor (4588,3240)", at(4588, 3240)),
        ("canyon 14's floor (1835,4562)", at(1835, 4562)),
    ];
    let recv = |dir: &[u8], c: usize| -> Option<usize> {
        let d = dir[c];
        if d == DIR_NONE {
            return None;
        }
        let nx = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        Some(ny * w + nx)
    };
    struct Read {
        g: GridF32,
        bre: GridF32,
        lake_map: Vec<u32>,
        levels: HashMap<u32, (f32, usize)>,
        dir: Vec<u8>,
        canyons: Vec<(u32, f32, usize)>,
    }
    let b1f = ValleyConstruction { a_min_km2: 0.1, wall_profile: Some(prof), ..base };
    let worlds: [(&str, ValleyConstruction); 7] = [
        ("témoin C2/10 col", base),
        ("B2 → A_c (no clause)", ac),
        ("B2 → A_c + FULL clause", ValleyConstruction { trunk_band: true, ..ac }),
        ("B2 → A_c + CONCORDANT clause", ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..ac }),
        ("B2 (A_c) + B1 + foot (no clause)", b1f),
        ("B2 (A_c) + B1 + foot + FULL clause", ValleyConstruction { trunk_band: true, ..b1f }),
        ("B2 (A_c) + B1 + foot + CONCORDANT clause", ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..b1f }),
    ];
    // ADR Finding 131-A -- the union of disks of radius 2 W around the concordant clause's junctions, per family
    let junction_mask = |fam: ValleyConstruction| -> Vec<bool> {
        use ymir_core::tectonics_c1::valley_construction::carve_diag;
        let vc = ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..fam };
        let skf = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, _, dg) = carve_diag(&pre, &skf, &vc, &ss);
        let banded = dg.banded.expect("the clause is on");
        let junction_on = |small: u32, big: u32| -> Option<u32> {
            let mut l = small;
            for _ in 0..skf.polylines.len() {
                match skf.line_parent[l as usize] {
                    Some((p, j)) if p == big => return Some(j),
                    Some((p, _)) => l = p,
                    None => return None,
                }
            }
            None
        };
        let mut js = std::collections::HashSet::new();
        for k in 0..n {
            let s = banded[k];
            if s != u32::MAX {
                let big = dg.line_of[s as usize];
                if let Some(j) = junction_on(dg.line_of[dg.who[k] as usize], big) {
                    js.insert((big, j));
                }
            }
        }
        let mut m = vec![false; n];
        for &(big, j) in &js {
            let (jx, jy, _, hw) = skf.polylines[big as usize][j as usize];
            let r = 2.0 * 2.0 * hw / (CELL_KM * 1000.0);
            for yy in (jy - r).floor() as i64..=(jy + r).ceil() as i64 {
                for xx in (jx - r).floor() as i64..=(jx + r).ceil() as i64 {
                    if ((xx as f32 + 0.5 - jx).powi(2) + (yy as f32 + 0.5 - jy).powi(2)).sqrt() <= r {
                        m[yy.rem_euclid(h as i64) as usize * w + xx.rem_euclid(w as i64) as usize] = true;
                    }
                }
            }
        }
        eprintln!("   the concordant clause's junctions: {} · the 2 W union covers {} cells", js.len(), m.iter().filter(|&&b| b).count());
        m
    };
    let mask_ac = junction_mask(ac);
    let mask_b1 = junction_mask(b1f);
    let near = 2.0 / CELL_KM;
    let mut keep: HashMap<&str, Read> = HashMap::new();
    let mut g_t: Option<GridF32> = None;
    let mut summary = Vec::new();
    for (label, vc) in worlds {
        let t = Instant::now();
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let secs = t.elapsed().as_secs_f64();
        eprintln!("\n────────── {label} · build {secs:.0} s ──────────");
        let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &sk, &vc, &ss);
        let cls = classes(&pre, &sk, &mk);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw));
        let co = coast(&g, &ss, &dwall, near);
        let nn = co.near.iter().filter(|&&b| b).count();
        let rn = nn as f32 / co.l_near_km.max(1e-6);
        let sig = sigma_p50(&g, &ss);
        let r8w = r8_terrain(&g, &cls, 16);
        let r8t = r8_terrain(g_t.as_ref().unwrap_or(&g), &cls, 16);
        eprintln!(
            "   COAST {} spurs · near a coastal wall **{rn:.4} /km** · σ p50 **{sig:.3} m** · R8 terrain {r8w:.4} (témoin on the same classes {r8t:.4})",
            co.mids.len()
        );
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = Vec::new();
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            let dfill = fill_del[floor] - fill_pre[floor];
            if over_dug_depression(dfill, b.rim) {
                eprintln!(
                    "   OVER-DUG body {} · {:.2} km² (domain) · floor ({},{}) · Δfill {dfill:.1} m · class {}",
                    b.id,
                    b.km2,
                    floor % w,
                    floor / w,
                    CLS_NAME.get(cls[floor] as usize).copied().unwrap_or("sea")
                );
                canyons.push((b.id, b.km2, floor));
            }
        }
        eprintln!("   over-dug bodies: {}", canyons.len());
        drop(fill_del);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage
        }));
        let mut row = format!(
            "{label:<42} over-dug {} {:?} · spurs near {rn:.4}/km · σ {sig:.3} m · R8 terr {r8w:.4}",
            canyons.len(),
            canyons.iter().map(|c| (c.2 % w, c.2 / w)).collect::<Vec<_>>()
        );
        let heights: Vec<String> = fixed.iter().map(|(nm, k)| format!("{nm} {:.1} m", mm(&g, *k))).collect();
        eprintln!("   heights: {}", heights.join(" · "));
        match res {
            Ok(dr) => {
                let segs = &dr.rivers.segments;
                let wall: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k]).collect();
                let wc: Vec<usize> = (0..segs.len())
                    .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2)
                    .collect();
                let teeth = wc
                    .iter()
                    .filter(|&&i| {
                        let p = &segs[i].points;
                        p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32 >= 0.8 * p.len() as f32
                    })
                    .count();
                let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
                eprintln!("   teeth {teeth} · R8 network (chord 8) {:.4} · water bodies (lakes) {} · Finding 38 holds", r8_chords(&pts, 8), dr.lakes.len());
                // ADR Finding 131-A -- R8 within 2 W of the concordant clause's junctions
                for (mname, mask) in [("B2 → A_c's junctions", &mask_ac), ("B2 + B1 + foot's junctions", &mask_b1)] {
                    if (label.starts_with("B2 → A_c") && !mname.starts_with("B2 → A_c")) || (label.starts_with("B2 (A_c)") && !mname.starts_with("B2 + B1")) {
                        continue;
                    }
                    let cls_m: Vec<u8> = (0..n).map(|k| if mask[k] && g.data[k] > SEA { 0 } else { 255 }).collect();
                    let r8t = r8_terrain(&g, &cls_m, 8);
                    let mut runs: Vec<Vec<(u32, u32)>> = Vec::new();
                    for &i in &wc {
                        let mut cur: Vec<(u32, u32)> = Vec::new();
                        for &(x, y) in own(&segs[i]) {
                            if mask[y as usize * w + x as usize] {
                                cur.push((x, y));
                            } else if !cur.is_empty() {
                                runs.push(std::mem::take(&mut cur));
                            }
                        }
                        if !cur.is_empty() {
                            runs.push(cur);
                        }
                    }
                    let rr: Vec<&[(u32, u32)]> = runs.iter().map(|r| r.as_slice()).collect();
                    let pts_in: usize = runs.iter().map(|r| r.len()).sum();
                    eprintln!("   within 2 W of {mname}: R8 terrain (8 × 8 windows) **{r8t:.4}** · R8 network (chord 8) **{:.4}** ({pts_in} watercourse points)", r8_chords(&rr, 8));
                    row += &format!(" · [{mname}] R8t {r8t:.4} R8n {:.4}", r8_chords(&rr, 8));
                }
                row += &format!(" · teeth {teeth} · R8 net {:.4} · lakes {}", r8_chords(&pts, 8), dr.lakes.len());
                {
                    let levels = dr
                        .lakes
                        .iter()
                        .map(|l| (l.base.id, (l.level_m, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize)))
                        .collect();
                    keep.insert(label, Read { g: g.clone(), bre: bre.clone(), lake_map: dr.lake_map.clone(), levels, dir: dr.flow.direction.clone(), canyons: canyons.clone() });
                }
            }
            Err(e) => {
                let msg = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
                eprintln!("   ⛔ the HD assembly PANICKED: {}", msg.lines().next().unwrap_or(""));
                row += " · ⛔ Finding 38";
            }
        }
        row += &format!(" · build {secs:.0} s");
        summary.push(row);
        if g_t.is_none() {
            g_t = Some(g);
        }
    }
    // ── canyon 13's dam: the témoin's drainage path from its floor, walked under the clause ──
    eprintln!("\n   ── canyon 13 (floor (2281,4563)): the témoin's path walked in each world ──");
    let f13 = at(2281, 4563);
    let tm = &keep["témoin C2/10 col"];
    let mut labels: Vec<&&str> = keep.keys().filter(|l| !l.starts_with("témoin")).collect();
    labels.sort();
    for label in labels {
        let r = &keep[*label];
        let persists = r.canyons.iter().any(|&(_, _, fl)| {
            let (dx, dy) = ((fl % w) as i64 - 2281, (fl / w) as i64 - 4563);
            dx * dx + dy * dy <= 400
        });
        eprintln!("   {label}: canyon 13 {}", if persists { "PERSISTS" } else { "is GONE" });
        let lid = r.lake_map[f13];
        let level = r.levels.get(&lid).map_or(f32::NAN, |v| v.0);
        let mut c = f13;
        let mut dam: Option<(usize, usize, f32)> = None;
        for step in 0..20_000 {
            let Some(nx) = recv(&tm.dir, c) else { break };
            c = nx;
            if tm.g.data[c] <= SEA {
                break;
            }
            if lid != 0 && r.lake_map[c] != lid && mm(&r.bre, c) > level + 0.01 {
                dam = Some((c, step, mm(&r.bre, c)));
                break;
            }
        }
        match dam {
            Some((c, step, z)) => eprintln!(
                "   {label}: lake {lid} at {level:.1} m · the témoin's path is DAMMED at ({},{}), {step} cells down, at {z:.1} m (témoin {:.1} m)",
                c % w,
                c / w,
                mm(&tm.bre, c)
            ),
            None => eprintln!("   {label}: lake {lid} (level {level:.1} m) · no dam on the témoin's path above the lake level"),
        }
    }
    eprintln!("\n   ── summary ──");
    for s in summary {
        eprintln!("   {s}");
    }
    eprintln!("\n==========  end Finding 131-A (worlds) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 131-P, where: the témoin's construction with LTD's pointers on the lakes' flats, against the
/// identity patch — each changed cell OUTSIDE the lake footprints classed upstream of a lake or not, and by
/// whether its skeleton χ or base changed. Declared in `f131_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f131_p_where --nocapture
#[test]
#[ignore]
fn f131_p_where() {
    use std::collections::VecDeque;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{ltd_directions_masked, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 131-P . where the lakes' flats act on the construction  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let cell_m = CELL_KM * 1000.0;
    let dc = C1DrainageConfig::default();
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    let lake = d.lake_map.clone();
    drop(d);
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let f = flow.filled.data.clone();
    let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
    let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
    let mut comp = vec![u32::MAX; n];
    let mut stats: Vec<(f64, f64, f64)> = Vec::new();
    for c0 in 0..n {
        if !needs[c0] || comp[c0] != u32::MAX {
            continue;
        }
        let id = stats.len() as u32;
        let (mut gmin, mut gmax, mut gap) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY);
        let mut q = VecDeque::from([c0]);
        comp[c0] = id;
        while let Some(c) = q.pop_front() {
            gmin = gmin.min(g[c]);
            gmax = gmax.max(g[c]);
            for k in 0..8 {
                let m = nb(c, k);
                if !is_ocean[m] && f[m] > f[c] {
                    gap = gap.min(f[m] as f64 - f[c] as f64);
                }
                if needs[m] && f[m] == f[c] && comp[m] == u32::MAX {
                    comp[m] = id;
                    q.push_back(m);
                }
            }
        }
        stats.push((gmin, gmax, if gap.is_finite() { gap } else { 1e-7 }));
    }
    let mut z = f.iter().map(|&v| v as f64).collect::<Vec<f64>>();
    for c in 0..n {
        if comp[c] != u32::MAX {
            let (gmin, gmax, gap) = stats[comp[c] as usize];
            z[c] += 0.1 * gap * (g[c] - gmin + 1.0) / (gmax - gmin + 2.0);
        }
    }
    let (ltd, _) = ltd_directions_masked(&z, &is_ocean, None, &flow.direction, w, h, cell_m);
    drop((z, flow, f, g, needs, comp));
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let id = |_: &GridF32, _: &[u32], _: &mut Vec<u8>| {};
    let lt = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = ltd[c];
            }
        }
    };
    let sk0 = skeleton_patched(&s1, &base, &ss, DOMAIN_KM, Some(&id));
    let sk1 = skeleton_patched(&s1, &base, &ss, DOMAIN_KM, Some(&lt));
    let (b0, _) = carve(&s1, &sk0, &base, &ss);
    let (b1, _) = carve(&s1, &sk1, &base, &ss);
    // upstream of a lake: the unpatched D8 path reaches a lake footprint cell
    let mut up = vec![0u8; n]; // 0 unknown, 1 upstream, 2 not
    let mut path = Vec::new();
    for s in 0..n {
        if up[s] != 0 {
            continue;
        }
        path.clear();
        let mut c = s;
        let verdict;
        loop {
            if up[c] != 0 {
                verdict = up[c];
                break;
            }
            if lake[c] != 0 {
                verdict = 1;
                break;
            }
            path.push(c);
            let dd = sk0.direction[c];
            if dd == DIR_NONE || s1.data[c] <= SEA {
                verdict = 2;
                break;
            }
            c = nb(c, dd as usize);
            if path.len() > n {
                panic!("cycle");
            }
        }
        for &p in &path {
            up[p] = verdict;
        }
    }
    let changed: Vec<usize> = (0..n).filter(|&k| lake[k] == 0 && b0.data[k] != b1.data[k]).collect();
    let upstream = changed.iter().filter(|&&k| up[k] == 1).count();
    let chi_ch = changed.iter().filter(|&&k| sk0.chi_m[k].to_bits() != sk1.chi_m[k].to_bits()).count();
    let base_ch = changed.iter().filter(|&&k| sk0.base_alt_m[k].to_bits() != sk1.base_alt_m[k].to_bits()).count();
    let up_chi = changed.iter().filter(|&&k| up[k] == 1 && sk0.chi_m[k].to_bits() != sk1.chi_m[k].to_bits()).count();
    let area_ch = changed.iter().filter(|&&k| sk0.area_km2[k].to_bits() != sk1.area_km2[k].to_bits()).count();
    let mut dz: Vec<f32> = changed.iter().map(|&k| (c1_altitude_norm_to_metres(b1.data[k], &ss) - c1_altitude_norm_to_metres(b0.data[k], &ss)).abs()).collect();
    dz.sort_by(f32::total_cmp);
    let qq = |p: f32| if dz.is_empty() { f32::NAN } else { dz[((dz.len() - 1) as f32 * p) as usize] };
    let m = changed.len().max(1) as f32;
    eprintln!(
        "   {} cells changed outside the lake footprints · UPSTREAM of a lake **{upstream}** ({:.1} %) · skeleton χ changed **{chi_ch}** ({:.1} %) · upstream AND χ changed {up_chi} ({:.1} %) · base changed {base_ch} ({:.1} %) · drained area changed {area_ch} ({:.1} %)",
        changed.len(),
        100.0 * upstream as f32 / m,
        100.0 * chi_ch as f32 / m,
        100.0 * up_chi as f32 / m,
        100.0 * base_ch as f32 / m,
        100.0 * area_ch as f32 / m
    );
    eprintln!("   |Δz| of the construction there: p50 / p90 / max {:.2} / {:.2} / {:.2} m", qq(0.5), qq(0.9), qq(1.0));
    eprintln!("\n==========  end Finding 131-P (where)  ==========\n");
}
