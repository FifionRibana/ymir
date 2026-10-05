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

// ════════════════════════════════ Finding 132 ════════════════════════════════

/// ADR Finding 132-P — the hidden parameter: the pointers of the lakes' flat cells, resolved two different
/// arbitrary (acyclic) ways and by LTD, and the construction compared outside the lakes (P1); the lakes whose flats
/// move χ against the delivered world's lakes (P2); the zone they move against Finding 131's junctions and cols
/// (P3). Declared in `f132_predictions.md` before the run. Writes the zone masks to the temp dir (`f132_zone_*.bin`).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f132_p --nocapture
#[test]
#[ignore]
fn f132_p() {
    use std::collections::{HashMap, VecDeque};
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, classify_lakes_water_balance};
    use ymir_core::tectonics_c1::valley_construction::{carve_diag, ltd_directions_masked, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 132-P . the lakes' flats, a hidden parameter  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let cell_m = CELL_KM * 1000.0;
    let cell_km2 = CELL_KM * CELL_KM;
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
    // LTD on the Garbrecht-Martz elevation (Finding 130-P1)
    let ltd = {
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
        ltd_directions_masked(&z, &is_ocean, None, &flow.direction, w, h, cell_m).0
    };
    let arbitrary = |seed: u64| -> Vec<u8> {
        (0..n)
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
                let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                cand[(hsh % cand.len() as u64) as usize]
            })
            .collect()
    };
    let (arb_a, arb_b) = (arbitrary(1), arbitrary(2));
    let nt = target.iter().filter(|&&b| b).count();
    let diff_ab = (0..n).filter(|&c| target[c] && arb_a[c] != arb_b[c]).count();
    eprintln!("   {nt} flat lake cells · pointers differing between the two arbitrary resolutions: {diff_ab}");
    drop((flow, f, g, needs));
    let lk: Vec<bool> = lake.iter().map(|&l| l != 0).collect();
    let dlake = dist_from(&lk, w, h);
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    for (label, tag, vc) in [("témoin C2/10 col", "temoin", base), ("B2 → A_c", "b2", ac)] {
        let t = Instant::now();
        eprintln!("\n   ── {label} ──");
        let id = |_: &GridF32, _: &[u32], _: &mut Vec<u8>| {};
        let sk_id = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&id));
        let (z_id, _) = carve(&s1, &sk_id, &vc, &ss);
        let (z_np, _) = carve(&s1, &skeleton(&s1, &vc, &ss, DOMAIN_KM), &vc, &ss);
        let ctrl = (0..n).filter(|&k| z_np.data[k] != z_id.data[k]).count();
        drop(z_np);
        eprintln!("   the control (identity patch against no patch): {ctrl} cells changed");
        let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
            for c in 0..n {
                if target[c] {
                    dd[c] = arb_a[c];
                }
            }
        };
        let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
            for c in 0..n {
                if target[c] {
                    dd[c] = arb_b[c];
                }
            }
        };
        let pl = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
            for c in 0..n {
                if target[c] {
                    dd[c] = ltd[c];
                }
            }
        };
        let (za, _) = carve(&s1, &skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pa)), &vc, &ss);
        let (zb, _) = carve(&s1, &skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pb)), &vc, &ss);
        let (zl, _) = carve(&s1, &skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pl)), &vc, &ss);
        // the lake each cell's unpatched D8 path enters first (0: none)
        let mut first = vec![u32::MAX; n];
        let mut path = Vec::new();
        for s in 0..n {
            if first[s] != u32::MAX {
                continue;
            }
            path.clear();
            let mut c = s;
            let v;
            loop {
                if first[c] != u32::MAX {
                    v = first[c];
                    break;
                }
                if lake[c] != 0 {
                    v = lake[c];
                    break;
                }
                path.push(c);
                let dd = sk_id.direction[c];
                if dd == DIR_NONE || s1.data[c] <= SEA {
                    v = 0;
                    break;
                }
                c = nb(c, dd as usize);
            }
            for &p in &path {
                first[p] = v;
            }
        }
        let describe = |name: &str, x: &GridF32, y: &GridF32| -> Vec<usize> {
            let ch: Vec<usize> = (0..n).filter(|&k| lake[k] == 0 && x.data[k] != y.data[k]).collect();
            let mut dz: Vec<f32> = ch.iter().map(|&k| (mm(x, k) - mm(y, k)).abs()).collect();
            let far = ch.iter().map(|&k| dlake[k] as f32 * cell_m).fold(0f32, f32::max);
            let (p50, p90, mx) = (q(&mut dz, 0.5), q(&mut dz, 0.9), q(&mut dz, 1.0));
            eprintln!(
                "   {name:<24} OUTSIDE the lakes: **{}** cells · |Δz| p50 {p50:.2} / **p90 {p90:.2}** / max {mx:.1} m · farthest from a lake {far:.0} m{}",
                ch.len(),
                if p90 < 1.0 { " · STOP RULE (p90 < 1 m)" } else { "" }
            );
            ch
        };
        let ch_ab = describe("arbitrary A vs B", &za, &zb);
        describe("LTD vs arbitrary A", &zl, &za);
        describe("LTD vs arbitrary B", &zl, &zb);
        // per lake (A vs B)
        let mut per: HashMap<u32, (usize, Vec<f32>, f32)> = HashMap::new();
        for &k in &ch_ab {
            let e = per.entry(first[k]).or_insert((0, Vec::new(), 0.0));
            e.0 += 1;
            e.1.push((mm(&za, k) - mm(&zb, k)).abs());
            e.2 = e.2.max(dlake[k] as f32 * cell_m);
        }
        // P2: the delivered world's lakes
        let gw = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let dd = c1_drainage_windowed(&gw, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&gw, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        let verdicts = classify_lakes_water_balance(&bre, &dr.flow, &dc_b, cell_km2, &dr.lakes, &dr.lake_map, None, w, h);
        let wc = water_class(&gw, SEA);
        let final_lake: HashMap<u32, usize> = dr.lakes.iter().enumerate().map(|(i, l)| (l.base.id, i)).collect();
        let mut rows: Vec<(u32, usize)> = per.iter().filter(|(k, _)| **k != 0).map(|(k, v)| (*k, v.0)).collect();
        rows.sort_by(|a, b| b.1.cmp(&a.1));
        let (mut present, mut absent) = (0usize, 0usize);
        eprintln!("   P2 — the lakes whose flats move the construction (A vs B), against this world's lakes ({} final lakes):", dr.lakes.len());
        for &(lid, cnt) in &rows {
            let flats: Vec<usize> = (0..n).filter(|&c| target[c] && lake[c] == lid).collect();
            let mut votes: HashMap<u32, usize> = HashMap::new();
            for &c in &flats {
                *votes.entry(dr.lake_map[c]).or_insert(0) += 1;
            }
            let (best, bv) = votes.iter().filter(|(k, _)| **k != 0).max_by_key(|(_, v)| **v).map(|(k, v)| (*k, *v)).unwrap_or((0, 0));
            let in_any = flats.iter().filter(|&&c| dr.lake_map[c] != 0).count();
            let wc2 = flats.iter().filter(|&&c| wc[c] == WATER_CLASS_INLAND).count();
            let e = per.get(&lid).expect("row");
            let mut dz = e.1.clone();
            let head = format!(
                "      pre-lake {lid:>8} · {} flat cells · {cnt} construction cells moved (p90 |Δz| {:.1} m, farthest {:.0} m) ·",
                flats.len(),
                q(&mut dz, 0.9),
                e.2
            );
            if bv as f32 >= 0.5 * flats.len() as f32 && best != 0 {
                present += 1;
                let l = &dr.lakes[final_lake[&best]];
                let a_eq = verdicts.iter().find(|v| v.id == best).map_or(f32::NAN, |v| v.a_eq_km2);
                eprintln!(
                    "{head} **PRESENT** as final lake {best} ({:.0} % of its flats): {:.1} km², depth {:.1} m, level {:.1} m, {:?}, a_eq {:.1} km²",
                    100.0 * bv as f32 / flats.len() as f32,
                    l.area_km2,
                    l.depth_m,
                    l.level_m,
                    l.lake_type,
                    a_eq
                );
            } else {
                absent += 1;
                eprintln!(
                    "{head} **ABSENT** (under any final lake {:.0} %, water class 2 {:.0} %)",
                    100.0 * in_any as f32 / flats.len().max(1) as f32,
                    100.0 * wc2 as f32 / flats.len().max(1) as f32
                );
            }
        }
        let none_attr = per.get(&0).map_or(0, |v| v.0);
        eprintln!(
            "   P2: {} lakes move the construction · PRESENT **{present}** · ABSENT **{absent}** · changed cells attributed to no lake {none_attr}{}",
            rows.len(),
            if present == 0 { " · STOP RULE (no lake present)" } else { "" }
        );
        // the zone (every cell whose z differs between A and B), written for f132_a
        let zone: Vec<u8> = (0..n).map(|k| (za.data[k] != zb.data[k]) as u8).collect();
        std::fs::write(std::env::temp_dir().join(format!("f132_zone_{tag}.bin")), &zone).expect("write the zone");
        let nz = zone.iter().filter(|&&b| b == 1).count();
        eprintln!("   zone P: {nz} cells ({:.1} % of the grid)", 100.0 * nz as f32 / n as f32);
        if tag == "b2" {
            let conc = ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..ac };
            let skc = skeleton(&s1, &conc, &ss, DOMAIN_KM);
            let (_, _, dg) = carve_diag(&s1, &skc, &conc, &ss);
            let banded = dg.banded.expect("the clause is on");
            let junction_on = |small: u32, big: u32| -> Option<u32> {
                let mut l = small;
                for _ in 0..skc.polylines.len() {
                    match skc.line_parent[l as usize] {
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
            let inz = js
                .iter()
                .filter(|&&(big, j)| {
                    let p = skc.polylines[big as usize][j as usize];
                    let c = (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize;
                    zone[c] == 1
                })
                .count();
            eprintln!(
                "   P3: the concordant clause's junctions IN the zone **{inz} of {}** ({:.1} %)",
                js.len(),
                100.0 * inz as f32 / js.len().max(1) as f32
            );
            for (nm, x, y) in [("col 7", 4551usize, 3280usize), ("col 14", 1805, 4580), ("canyon 14's floor", 1835, 4562), ("13's dam cell", 2315, 4614)] {
                let k = y * w + x;
                eprintln!(
                    "      {nm} ({x},{y}): {} · |Δz(A,B)| {:.2} m · its D8 path enters pre-lake {}",
                    if zone[k] == 1 { "IN the zone" } else { "outside" },
                    (mm(&za, k) - mm(&zb, k)).abs(),
                    first[k]
                );
            }
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n==========  end Finding 132-P . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 132-A — rule 16 on Finding 131's concordant clause, read-only: its five criteria split between
/// the junctions whose tributary's head lies in the band (ii) and the others (i); canyon 14's dam located; the
/// figures under the lakes' flats noise (zone P of `f132_p`). Declared in `f132_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f132_a --nocapture
#[test]
#[ignore]
fn f132_a() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::valley_construction::carve_diag;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 132-A . the concordant clause, by population (rule 16)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let cell_m = CELL_KM * 1000.0;
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let prof = WallProfile { foot_quiet: Some(0.5), ..WallProfile::f124() };
    let sk = skeleton(&s1, &base, &ss, DOMAIN_KM);
    let law: Vec<f32> = (0..n).map(|k| if sk.chi_m[k].is_finite() { sk.floor_m(k, base.age_k) } else { f32::NAN }).collect();
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let sk_ac = skeleton(&s1, &ac, &ss, DOMAIN_KM);
    let zone: Vec<u8> = std::fs::read(std::env::temp_dir().join("f132_zone_b2.bin")).unwrap_or_default();
    let has_zone = zone.len() == n;
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let cell_of = |p: (f32, f32, f32, f32)| (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize;
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
    // per family: the populations, and the masks for the pipeline worlds
    let mut fam_masks: Vec<(String, Vec<bool>, Vec<bool>, Vec<u32>)> = Vec::new(); // (family, mask (i), mask (ii), clause cell -> population 1/2)
    for (tag, world) in [("B2 → A_c", ac), ("B2 (A_c) + B1 + foot", ValleyConstruction { wall_profile: Some(prof), ..ac })] {
        eprintln!("\n   ── {tag} (construction) ──");
        let (b0, m0) = carve(&s1, &sk_ac, &world, &ss);
        let z0 = metres(&b0);
        let conc = ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..world };
        let (bc, mc, dg) = carve_diag(&s1, &sk_ac, &conc, &ss);
        let zc = metres(&bc);
        let banded = dg.banded.expect("the clause is on");
        // junction key and owner line of every clause cell
        let mut key_of: HashMap<(u32, u32), (Vec<usize>, u32)> = HashMap::new();
        for k in 0..n {
            let s = banded[k];
            if s == u32::MAX {
                continue;
            }
            let big = dg.line_of[s as usize];
            let owner = dg.line_of[dg.who[k] as usize];
            let j = junction_on(owner, big).expect("a clause cell has a junction");
            key_of.entry((big, j)).or_insert((Vec::new(), owner)).0.push(k);
        }
        // population (ii): the owner line has a moved cell among its first 4 cells
        let head_moved = |owner: u32| -> bool {
            let l = &sk_ac.polylines[owner as usize];
            let mut cells: Vec<usize> = Vec::new();
            for &p in l {
                let c = cell_of(p);
                if cells.last() != Some(&c) {
                    cells.push(c);
                }
                if cells.len() >= 4 {
                    break;
                }
            }
            cells.iter().any(|&c| zc[c] != z0[c])
        };
        let mut pop_of = vec![0u32; n];
        let (mut nj1, mut nj2) = (0usize, 0usize);
        let mut junc: Vec<((u32, u32), u32, Vec<usize>)> = Vec::new();
        for (key, (cells, owner)) in key_of {
            let p = if head_moved(owner) { 2 } else { 1 };
            if p == 1 {
                nj1 += 1;
            } else {
                nj2 += 1;
            }
            for &k in &cells {
                pop_of[k] = p;
            }
            junc.push((key, p, cells));
        }
        eprintln!("   junctions: (i) tributary head outside the band {nj1} · (ii) head in the band {nj2}");
        let flat = |z: &[f32], c: usize| !(0..8).any(|k| z[nb(c, k)] < z[c]);
        let (mut m1, mut m2) = (vec![false; n], vec![false; n]);
        for pop in [1u32, 2] {
            // the construction with the clause on this population only (exact: carve decides every cell alone)
            let zp: Vec<f32> = (0..n).map(|k| if pop_of[k] == pop { zc[k] } else { z0[k] }).collect();
            let carved: Vec<bool> = (0..n).map(|k| if pop_of[k] == pop { mc.carved[k] } else { m0.carved[k] }).collect();
            let raised = trunk10.iter().filter(|&&k| zp[k] - law[k] > 1.0).count();
            let (t, lo, hi, m) = theta_laid_ci(&zp, &carved, &sk);
            let (mut rep_c, mut rep_0, mut nj) = (0usize, 0usize, 0usize);
            let mut len_w = Vec::new();
            for (key, p, cells) in &junc {
                if *p != pop {
                    continue;
                }
                nj += 1;
                let fc: Vec<usize> = cells.iter().copied().filter(|&k| flat(&zp, k)).collect();
                if !fc.is_empty() {
                    rep_c += 1;
                    let (jx, jy, _, jhw) = sk_ac.polylines[key.0 as usize][key.1 as usize];
                    let far = fc
                        .iter()
                        .map(|&k| {
                            let (dx, dy) = ((k % w) as f32 + 0.5 - jx, (k / w) as f32 + 0.5 - jy);
                            (dx * dx + dy * dy).sqrt() * cell_m
                        })
                        .fold(0f32, f32::max);
                    len_w.push(far / (2.0 * jhw));
                }
                if cells.iter().any(|&k| flat(&z0, k)) {
                    rep_0 += 1;
                }
                // the 2 W disks of this population
                let (jx, jy, _, jhw) = sk_ac.polylines[key.0 as usize][key.1 as usize];
                let r = 4.0 * jhw / cell_m;
                let mm = if pop == 1 { &mut m1 } else { &mut m2 };
                for yy in (jy - r).floor() as i64..=(jy + r).ceil() as i64 {
                    for xx in (jx - r).floor() as i64..=(jx + r).ceil() as i64 {
                        if ((xx as f32 + 0.5 - jx).powi(2) + (yy as f32 + 0.5 - jy).powi(2)).sqrt() <= r {
                            mm[yy.rem_euclid(h as i64) as usize * w + xx.rem_euclid(w as i64) as usize] = true;
                        }
                    }
                }
            }
            let in_zone = if has_zone {
                junc.iter()
                    .filter(|(key, p, _)| *p == pop && zone[cell_of(sk_ac.polylines[key.0 as usize][key.1 as usize])] == 1)
                    .count()
            } else {
                0
            };
            eprintln!(
                "   population ({}) — the clause on it alone: raised > 1 m **{raised}** · θ laid **{t:.3}** [{lo:.3}, {hi:.3}] ({m}) · replat {rep_c} of {nj} junctions ({:.1} %) against {rep_0} without ({:.1} %) · replat length p50 {:.2} W · junctions in zone P {in_zone} ({:.1} %)",
                if pop == 1 { "i" } else { "ii" },
                100.0 * rep_c as f32 / nj.max(1) as f32,
                100.0 * rep_0 as f32 / nj.max(1) as f32,
                q(&mut len_w, 0.5),
                100.0 * in_zone as f32 / nj.max(1) as f32
            );
        }
        fam_masks.push((tag.to_string(), m1, m2, pop_of));
    }
    // ── the pipeline worlds: R8 around each population's junctions; canyon 14's dam ──
    eprintln!("\n   ── the pipeline worlds (light pass + A1+B2) ──");
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    // the témoin's drainage path, on its own assembly's flow (as Finding 128 walked canyon 13)
    let (tm_dir, tm_g) = {
        let g = build_field_seed(Knobs { valley: Some(base), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        (dr.flow.direction, g)
    };
    drop(pre);
    let b1f = ValleyConstruction { wall_profile: Some(prof), ..ac };
    for (fi, (tag, fam)) in [("B2 → A_c", ac), ("B2 (A_c) + B1 + foot", b1f)].into_iter().enumerate() {
        let (_, m1, m2, pop_of) = &fam_masks[fi];
        for (wname, vc) in [("no clause", fam), ("CONCORDANT", ValleyConstruction { trunk_band: true, trunk_band_concordant: true, ..fam })] {
            let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
            let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
            let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
            let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
            let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
            let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
            let segs = &dr.rivers.segments;
            let wc: Vec<usize> = (0..segs.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2)
                .collect();
            let mut line = format!("   {tag} · {wname:<11}");
            for (pname, mask) in [("(i)", m1), ("(ii)", m2)] {
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
                line += &format!(" · {pname} R8 terrain **{r8t:.4}** network **{:.4}**", r8_chords(&rr, 8));
            }
            eprintln!("{line}");
            // canyon 14: the témoin's path from its floor, walked in this world until a cell above the lake level
            if wname == "CONCORDANT" {
                // canyon 14's lake: the lake at its floor, or the nearest lake cell within 20 cells of it
                let mut f14 = 4562 * w + 1835;
                if dr.lake_map[f14] == 0 {
                    let mut best = (f32::INFINITY, f14);
                    for dy in -20i64..=20 {
                        for dx in -20i64..=20 {
                            let k = (4562 + dy) as usize * w + (1835 + dx) as usize;
                            let d = ((dx * dx + dy * dy) as f32).sqrt();
                            if dr.lake_map[k] != 0 && d < best.0 {
                                best = (d, k);
                            }
                        }
                    }
                    f14 = best.1;
                }
                let lid = dr.lake_map[f14];
                let level = dr.lakes.iter().find(|l| l.base.id == lid).map_or(f32::NAN, |l| l.level_m);
                let mut c = f14;
                let mut dam = None;
                for step in 0..20_000 {
                    let d = tm_dir[c];
                    if d == DIR_NONE {
                        break;
                    }
                    c = nb(c, d as usize);
                    if tm_g.data[c] <= SEA {
                        break;
                    }
                    if lid != 0 && dr.lake_map[c] != lid && c1_altitude_norm_to_metres(bre.data[c], &ss) > level + 0.01 {
                        dam = Some((c, step));
                        break;
                    }
                }
                match dam {
                    Some((c, step)) => eprintln!(
                        "      canyon 14: lake {lid} at {level:.1} m · the témoin's path is DAMMED at ({},{}), {step} cells down, at {:.1} m (témoin {:.1} m) · a clause cell: {} · in zone P: {}",
                        c % w,
                        c / w,
                        c1_altitude_norm_to_metres(bre.data[c], &ss),
                        c1_altitude_norm_to_metres(tm_g.data[c], &ss),
                        match pop_of[c] {
                            0 => "NO".to_string(),
                            p => format!("YES, population ({})", if p == 1 { "i" } else { "ii" }),
                        },
                        if has_zone { if zone[c] == 1 { "IN" } else { "out" } } else { "n/a" }
                    ),
                    None => eprintln!("      canyon 14: lake {lid} (level {level:.1} m) · no dam on the témoin's path above the lake level"),
                }
            }
        }
    }
    eprintln!("\n==========  end Finding 132-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// θ over the trunk ≥ 10 km² links (S > 1e-4), Finding 129's instrument, with a percentile bootstrap over the
/// links (1 000 resamples, fixed seed): `(θ, lo 2.5 %, hi 97.5 %, links)`.
fn theta_sa_ci(z: &[f32], sk: &Skeleton) -> (f32, f32, f32, usize) {
    let all = vec![true; z.len()];
    theta_laid_ci(z, &all, sk)
}

/// ADR Finding 132-P4 — the remedy test: with `lake_base` ON, the construction under two different arbitrary
/// resolutions of the lakes' flats (P1) must agree outside the lakes; any residual located. Témoin and B2 → A_c,
/// carve on S1. Declared in `f132_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f132_p4_remedy --nocapture
#[test]
#[ignore]
fn f132_p4_remedy() {
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 132-P4 . the remedy test: lake = base level, two arbitrary flats  ==========");
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
    let arbitrary = |seed: u64| -> Vec<u8> {
        (0..n)
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
                let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                cand[(hsh % cand.len() as u64) as usize]
            })
            .collect()
    };
    let (arb_a, arb_b) = (arbitrary(1), arbitrary(2));
    drop((flow, f, g, needs));
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    for (label, world) in [("témoin C2/10 col", base), ("B2 → A_c", ac)] {
        let t = Instant::now();
        eprintln!("\n   ── {label} ──");
        let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
            for c in 0..n {
                if target[c] {
                    dd[c] = arb_a[c];
                }
            }
        };
        let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
            for c in 0..n {
                if target[c] {
                    dd[c] = arb_b[c];
                }
            }
        };
        let id = |_: &GridF32, _: &[u32], _: &mut Vec<u8>| {};
        for (state, vc) in [("OFF", world), ("ON (lake = base level)", ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..world })] {
            let sk_id = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&id));
            let (z_id, _) = carve(&s1, &sk_id, &vc, &ss);
            let (z_np, _) = carve(&s1, &skeleton(&s1, &vc, &ss, DOMAIN_KM), &vc, &ss);
            let ctrl = (0..n).filter(|&k| z_np.data[k] != z_id.data[k]).count();
            drop(z_np);
            let (za, _) = carve(&s1, &skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pa)), &vc, &ss);
            let (zb, _) = carve(&s1, &skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pb)), &vc, &ss);
            // the first input lake each cell's (identity) D8 path enters
            let mut first = vec![u32::MAX; n];
            let mut path = Vec::new();
            for s in 0..n {
                if first[s] != u32::MAX {
                    continue;
                }
                path.clear();
                let mut c = s;
                let v;
                loop {
                    if first[c] != u32::MAX {
                        v = first[c];
                        break;
                    }
                    if lake[c] != 0 {
                        v = lake[c];
                        break;
                    }
                    path.push(c);
                    let dd = sk_id.direction[c];
                    if dd == DIR_NONE || s1.data[c] <= SEA {
                        v = 0;
                        break;
                    }
                    c = nb(c, dd as usize);
                }
                for &p in &path {
                    first[p] = v;
                }
            }
            let ch: Vec<usize> = (0..n).filter(|&k| lake[k] == 0 && za.data[k] != zb.data[k]).collect();
            let inside = (0..n).filter(|&k| lake[k] != 0 && za.data[k] != zb.data[k]).count();
            let up = ch.iter().filter(|&&k| first[k] != 0).count();
            let base_ch = ch.iter().count();
            let mut dz: Vec<f32> = ch.iter().map(|&k| (mm(&za, k) - mm(&zb, k)).abs()).collect();
            let (p50, p90, mx) = (q(&mut dz, 0.5), q(&mut dz, 0.9), q(&mut dz, 1.0));
            eprintln!(
                "   {state:<24} control {ctrl} · A vs B OUTSIDE the lakes **{base_ch}** cells (upstream of a lake {up}, attributed to no lake {}) · |Δz| p50 {p50:.2} / p90 {p90:.2} / max {mx:.1} m · inside the lakes {inside}",
                base_ch - up
            );
            if state.starts_with("ON") && !ch.is_empty() {
                // where the residual sits: distance to a lake, and whether the cell's skeleton base / χ / area moved
                let lk: Vec<bool> = lake.iter().map(|&l| l != 0).collect();
                let dl = dist_from(&lk, w, h);
                let ska = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pa));
                let skb = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pb));
                let chi_ch = ch.iter().filter(|&&k| ska.chi_m[k].to_bits() != skb.chi_m[k].to_bits()).count();
                let area_ch = ch.iter().filter(|&&k| ska.area_km2[k].to_bits() != skb.area_km2[k].to_bits()).count();
                let trunk_ch = ch.iter().filter(|&&k| ska.trunk[k] != skb.trunk[k]).count();
                let mut dd: Vec<f32> = ch.iter().map(|&k| dl[k] as f32 * cell_m).collect();
                eprintln!(
                    "      the residual: χ differs {chi_ch} · drained area differs {area_ch} · trunk membership differs {trunk_ch} · distance to a lake p50 / p90 / max {:.0} / {:.0} / {:.0} m",
                    q(&mut dd.clone(), 0.5),
                    q(&mut dd.clone(), 0.9),
                    q(&mut dd, 1.0)
                );
            }
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n==========  end Finding 132-P4 (remedy) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 132-P4 — rule 18's pair on the témoin C2/10 col, `lake_base` OFF and ON: canyons, coast, θ by stage
/// with a link bootstrap, R8 terrain and network, lakes, relief; the consistency counts (bases on lakes absent at
/// the end, final lakes without a base); rule 14 (the volume removed more / less, and where). Declared in
/// `f132_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f132_p4_worlds --nocapture
#[test]
#[ignore]
fn f132_p4_worlds() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::valley_construction::LakeBase;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 132-P4 . lake = base level, the témoin OFF / ON  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..off };
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM); // the θ instrument's skeleton (Finding 129)
    // the input lakes (the construction input's pre-drainage) and the lake each cell's path enters first
    let lake_in = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM).lake_map;
    let mut first = vec![u32::MAX; n];
    let mut path = Vec::new();
    for s in 0..n {
        if first[s] != u32::MAX {
            continue;
        }
        path.clear();
        let mut c = s;
        let v;
        loop {
            if first[c] != u32::MAX {
                v = first[c];
                break;
            }
            if lake_in[c] != 0 {
                v = lake_in[c];
                break;
            }
            path.push(c);
            let dd = sk.direction[c];
            if dd == DIR_NONE || s1.data[c] <= SEA {
                v = 0;
                break;
            }
            c = nb(c, dd as usize);
        }
        for &p in &path {
            first[p] = v;
        }
    }
    let near = 2.0 / CELL_KM;
    let mut worlds: Vec<(GridF32, Vec<u32>)> = Vec::new();
    for (label, vc) in [("témoin C2/10 col, lake_base OFF", off), ("témoin C2/10 col, lake_base ON", on)] {
        let t = Instant::now();
        eprintln!("\n────────── {label} ──────────");
        // θ by stage, with its CI
        let s2 = build_field_seed(Knobs { valley: Some(vc), no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
        let bre2 = {
            let d = c1_drainage_windowed(&s2, None, &dcfg(), &ss, DOMAIN_KM);
            breach_monotone(&s2, &d.flow.filled, &d.lake_map, SEA, w, h)
        };
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let mut th = Vec::new();
        for (stage, fld) in [("(i) as built", &s2), ("(ii) + breach", &bre2), ("(iii) the world", &g)] {
            let (tt, lo, hi, m) = theta_sa_ci(&metres(fld), &sk);
            th.push(format!("{stage} **{tt:.3}** [{lo:.3}, {hi:.3}] ({m})"));
        }
        eprintln!("   θ by stage: {}", th.join(" · "));
        drop((s2, bre2));
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mk) = carve(&pre, &skp, &vc, &ss);
        let cls = classes(&pre, &skp, &mk);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mk.carved[k] && !mk.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw));
        let co = coast(&g, &ss, &dwall, near);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let r8w = r8_terrain(&g, &cls, 16);
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
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons.push((floor % w, floor / w));
            }
        }
        drop(fill_del);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        let segs = &dr.rivers.segments;
        let wc: Vec<usize> = (0..segs.len()).filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2).collect();
        let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
        eprintln!(
            "   SUMMARY · over-dug (canyons) **{}** {:?} · spurs near a coastal wall {rn:.4} /km · R8 terrain {r8w:.4} · R8 network (chord 8) {:.4} · lakes **{}** ({:.0} km²) · Finding 38 holds",
            canyons.len(),
            canyons,
            r8_chords(&pts, 8),
            dr.lakes.len(),
            dr.lakes.iter().map(|l| l.area_km2).sum::<f32>()
        );
        // consistency (ON): bases on input lakes absent at the end; final lakes without a base
        if vc.lake_base.is_some() {
            let mut ids: Vec<u32> = lake_in.iter().copied().filter(|&l| l != 0).collect();
            ids.sort_unstable();
            ids.dedup();
            let mut absent = Vec::new();
            for &lid in &ids {
                let cells: Vec<usize> = (0..n).filter(|&c| lake_in[c] == lid).collect();
                let in_final = cells.iter().filter(|&&c| dr.lake_map[c] != 0).count();
                if (in_final as f32) < 0.5 * cells.len() as f32 {
                    absent.push((lid, cells.len(), 100.0 * in_final as f32 / cells.len() as f32));
                }
            }
            let mut unbased = Vec::new();
            for l in &dr.lakes {
                let cells: Vec<usize> = (0..n).filter(|&c| dr.lake_map[c] == l.base.id).collect();
                let covered = cells.iter().filter(|&&c| lake_in[c] != 0).count();
                if (covered as f32) < 0.5 * cells.len().max(1) as f32 {
                    unbased.push((l.base.id, l.area_km2, 100.0 * covered as f32 / cells.len().max(1) as f32));
                }
            }
            eprintln!(
                "   CONSISTENCY · {} input lakes carry a base · bases on a lake ABSENT at the end **{}** {:?} · final lakes WITHOUT a base **{}** {:?}",
                ids.len(),
                absent.len(),
                absent.iter().map(|a| (a.0, a.1, a.2.round() as i32)).collect::<Vec<_>>(),
                unbased.len(),
                unbased.iter().map(|u| (u.0, u.1.round() as i32, u.2.round() as i32)).collect::<Vec<_>>()
            );
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
        worlds.push((g, dr.lake_map.clone()));
    }
    // rule 14: the volume removed more / less, ON − OFF, and where
    let (goff, _) = &worlds[0];
    let (gon, _) = &worlds[1];
    let (zo, zn) = (metres(goff), metres(gon));
    let mut vol: HashMap<&str, (f64, f64, usize)> = HashMap::new();
    for k in 0..n {
        if s1.data[k] <= SEA || zo[k] == zn[k] {
            continue;
        }
        let place = if lake_in[k] != 0 {
            "inside an input lake"
        } else if first[k] != 0 {
            "upstream of an input lake"
        } else {
            "elsewhere"
        };
        let dv = (zn[k] - zo[k]) as f64 * (cell_km2 as f64) * 1e-3; // km³
        let e = vol.entry(place).or_insert((0.0, 0.0, 0));
        if dv > 0.0 {
            e.0 += dv; // less removed (ON higher)
        } else {
            e.1 += -dv; // more removed
        }
        e.2 += 1;
    }
    for place in ["upstream of an input lake", "inside an input lake", "elsewhere"] {
        let (less, more, c) = vol.get(place).copied().unwrap_or((0.0, 0.0, 0));
        eprintln!("   RULE 14 · {place:<28} {c} cells changed · removed LESS (ON higher) {less:.3} km³ · removed MORE {more:.3} km³");
    }
    eprintln!("\n==========  end Finding 132-P4 (worlds) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 132-P4, localisation (no correction): with `lake_base` ON, where the residual χ difference between
/// the two arbitrary resolutions of the lakes' flats is BORN. Each residual cell's D8 path (A's pointers) is walked
/// downstream to the first cell whose χ is equal in A and B; the cell just upstream of it is the origin, classed.
/// Témoin. Declared in `f132_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f132_p4_where --nocapture
#[test]
#[ignore]
fn f132_p4_where() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 132-P4 . where the residual is born  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
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
    let arbitrary = |seed: u64| -> Vec<u8> {
        (0..n)
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
                let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                cand[(hsh % cand.len() as u64) as usize]
            })
            .collect()
    };
    let (arb_a, arb_b) = (arbitrary(1), arbitrary(2));
    drop((flow, f, g, needs));
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), lake_base: Some(LakeBase::InputLakes), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_a[c];
            }
        }
    };
    let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_b[c];
            }
        }
    };
    let ska = skeleton_patched(&s1, &base, &ss, DOMAIN_KM, Some(&pa));
    let skb = skeleton_patched(&s1, &base, &ss, DOMAIN_KM, Some(&pb));
    let ptr_out = (0..n).filter(|&k| lake[k] == 0 && ska.direction[k] != skb.direction[k]).count();
    let ptr_in = (0..n).filter(|&k| lake[k] != 0 && ska.direction[k] != skb.direction[k]).count();
    let same = |k: usize| ska.chi_m[k].to_bits() == skb.chi_m[k].to_bits() && ska.base_alt_m[k].to_bits() == skb.base_alt_m[k].to_bits();
    let resid: Vec<usize> = (0..n).filter(|&k| lake[k] == 0 && ska.chi_m[k].is_finite() && !same(k)).collect();
    eprintln!("   pointers differing between A and B: outside the lakes {ptr_out} · inside {ptr_in} · residual cells (χ or base differs, outside the lakes) {}", resid.len());
    let mut origin_of: HashMap<usize, usize> = HashMap::new();
    let mut unresolved = 0usize;
    for &s in &resid {
        let mut c = s;
        let mut steps = 0;
        loop {
            let dd = ska.direction[c];
            if dd == DIR_NONE {
                unresolved += 1;
                break;
            }
            let r = nb(c, dd as usize);
            if same(r) || s1.data[r] <= SEA {
                *origin_of.entry(c).or_insert(0) += 1;
                break;
            }
            c = r;
            steps += 1;
            if steps > n {
                unresolved += 1;
                break;
            }
        }
    }
    let mut cls: HashMap<&str, (usize, usize)> = HashMap::new();
    for (&o, &cnt) in &origin_of {
        let r = nb(o, ska.direction[o] as usize);
        let k = if lake[o] != 0 {
            "origin IN an input lake"
        } else if ska.direction[o] != skb.direction[o] {
            "origin's pointer differs A / B"
        } else if lake[r] != 0 {
            "origin drains INTO an input lake (its lake entry)"
        } else if s1.data[r] <= SEA {
            "origin drains to the sea"
        } else {
            "other"
        };
        let e = cls.entry(k).or_insert((0, 0));
        e.0 += 1;
        e.1 += cnt;
    }
    let mut rows: Vec<_> = cls.into_iter().collect();
    rows.sort_by(|a, b| b.1 .1.cmp(&a.1 .1));
    eprintln!("   {} distinct origins · paths without an origin {unresolved}", origin_of.len());
    for (k, (o, c)) in rows {
        eprintln!("      {k:<52} origins {o} · residual cells behind them {c}");
    }
    // the largest origins, described
    let mut top: Vec<(usize, usize)> = origin_of.iter().map(|(&o, &c)| (o, c)).collect();
    top.sort_by(|a, b| b.1.cmp(&a.1));
    for &(o, cnt) in top.iter().take(8) {
        let r = nb(o, ska.direction[o] as usize);
        eprintln!(
            "      ({},{}) behind {cnt} · χ A/B {:.1} / {:.1} · base A/B {:.1} / {:.1} · lake {} · receiver ({},{}) lake {} χ {:.1} · pointer A/B {} / {}",
            o % w,
            o / w,
            ska.chi_m[o],
            skb.chi_m[o],
            ska.base_alt_m[o],
            skb.base_alt_m[o],
            lake[o],
            r % w,
            r / w,
            lake[r],
            ska.chi_m[r],
            ska.direction[o],
            skb.direction[o]
        );
    }
    eprintln!("\n==========  end Finding 132-P4 (where)  ==========\n");
}

// ════════════════════════════════ Finding 133 ════════════════════════════════

/// ADR Finding 133-R — reconciliation, measurements only: one lake per row (the input lakes, then the ON world's
/// unbased final lakes) with Finding 131/132's quantities, and R2's overlap of the F132 residual with the 19
/// unbased final lakes. Témoin and B2 → A_c. Declared in `f133_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f133_r --nocapture
#[test]
#[ignore]
fn f133_r() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 133-R . the reconciliation, one lake per row  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let dc = C1DrainageConfig::default();
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    let lake = d.lake_map.clone();
    drop(d);
    // basin_base's closed depressions on the breached input (the skeleton's own rule, 1 cm)
    let spill = ocean_flood(&bf);
    let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
    let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let f = flow.filled.data.clone();
    let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
    let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
    let arbitrary = |seed: u64| -> Vec<u8> {
        (0..n)
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
                let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                cand[(hsh % cand.len() as u64) as usize]
            })
            .collect()
    };
    let (arb_a, arb_b) = (arbitrary(1), arbitrary(2));
    drop((flow, f, g, needs));
    let mut ids: Vec<u32> = lake.iter().copied().filter(|&l| l != 0).collect();
    ids.sort_unstable();
    ids.dedup();
    let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_a[c];
            }
        }
    };
    let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_b[c];
            }
        }
    };
    let final_lakes = |vc: ValleyConstruction| -> (Vec<u32>, Vec<(u32, f32)>) {
        let gw = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let dd = c1_drainage_windowed(&gw, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&gw, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        let ls = dr.lakes.iter().map(|l| (l.base.id, l.area_km2)).collect();
        (dr.lake_map, ls)
    };
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    for (label, off) in [("témoin C2/10 col", base), ("B2 → A_c", ac)] {
        let t = Instant::now();
        eprintln!("\n   ── {label} ──");
        let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..off };
        // (a) OFF A vs B, attributed to the first input lake on the (OFF, identity) D8 path
        let id = |_: &GridF32, _: &[u32], _: &mut Vec<u8>| {};
        let sk0 = skeleton_patched(&s1, &off, &ss, DOMAIN_KM, Some(&id));
        let mut first = vec![u32::MAX; n];
        let mut path = Vec::new();
        for s in 0..n {
            if first[s] != u32::MAX {
                continue;
            }
            path.clear();
            let mut c = s;
            let v;
            loop {
                if first[c] != u32::MAX {
                    v = first[c];
                    break;
                }
                if lake[c] != 0 {
                    v = lake[c];
                    break;
                }
                path.push(c);
                let dd = sk0.direction[c];
                if dd == DIR_NONE || s1.data[c] <= SEA {
                    v = 0;
                    break;
                }
                c = nb(c, dd as usize);
            }
            for &p in &path {
                first[p] = v;
            }
        }
        drop(sk0);
        let (oa, _) = carve(&s1, &skeleton_patched(&s1, &off, &ss, DOMAIN_KM, Some(&pa)), &off, &ss);
        let (ob, _) = carve(&s1, &skeleton_patched(&s1, &off, &ss, DOMAIN_KM, Some(&pb)), &off, &ss);
        let mut moved: HashMap<u32, usize> = HashMap::new();
        for k in 0..n {
            if lake[k] == 0 && oa.data[k] != ob.data[k] {
                *moved.entry(first[k]).or_insert(0) += 1;
            }
        }
        drop((oa, ob));
        // (f) ON A vs B, the residual and its source
        let ska = skeleton_patched(&s1, &on, &ss, DOMAIN_KM, Some(&pa));
        let skb = skeleton_patched(&s1, &on, &ss, DOMAIN_KM, Some(&pb));
        let (na, _) = carve(&s1, &ska, &on, &ss);
        let (nb_, _) = carve(&s1, &skb, &on, &ss);
        let resid: Vec<usize> = (0..n).filter(|&k| lake[k] == 0 && na.data[k] != nb_.data[k]).collect();
        drop((na, nb_));
        let differs = |k: usize| ska.chi_m[k].to_bits() != skb.chi_m[k].to_bits() || ska.base_alt_m[k].to_bits() != skb.base_alt_m[k].to_bits();
        let mut origin: Vec<Option<usize>> = Vec::with_capacity(resid.len());
        for &s in &resid {
            let mut c = s;
            let mut seen = false;
            let mut last = None;
            for _ in 0..n {
                if differs(c) {
                    seen = true;
                    last = Some(c);
                } else if seen {
                    break;
                }
                let dd = ska.direction[c];
                if dd == DIR_NONE || s1.data[c] <= SEA {
                    break;
                }
                c = nb(c, dd as usize);
            }
            origin.push(if seen { last } else { None });
        }
        drop((ska, skb));
        // the worlds' final lakes
        let (lm_off, ls_off) = final_lakes(off);
        let (lm_on, ls_on) = final_lakes(on);
        // (b) present in the OFF world, (e) its OFF final lake is a below-sea basin
        let present_in = |lm: &[u32], lid: u32| -> (bool, u32) {
            let cells: Vec<usize> = (0..n).filter(|&c| target[c] && lake[c] == lid).collect();
            let cells = if cells.is_empty() { (0..n).filter(|&c| lake[c] == lid).collect() } else { cells };
            let mut votes: HashMap<u32, usize> = HashMap::new();
            for &c in &cells {
                *votes.entry(lm[c]).or_insert(0) += 1;
            }
            let (best, bv) = votes.iter().filter(|(k, _)| **k != 0).max_by_key(|(_, v)| **v).map(|(k, v)| (*k, *v)).unwrap_or((0, 0));
            (best != 0 && bv as f32 >= 0.5 * cells.len() as f32, best)
        };
        // the unbased final lakes of the ON world (no input lake covers >= 50 % of it)
        let mut unbased: Vec<u32> = Vec::new();
        for &(fid, _) in &ls_on {
            let cells: Vec<usize> = (0..n).filter(|&c| lm_on[c] == fid).collect();
            let covered = cells.iter().filter(|&&c| lake[c] != 0).count();
            if (covered as f32) < 0.5 * cells.len().max(1) as f32 {
                unbased.push(fid);
            }
        }
        // Finding 132's consistency rule: ALL of an input lake's cells, >= 50 % under any final lake
        let absent_on = ids
            .iter()
            .filter(|&&lid| {
                let cells: Vec<usize> = (0..n).filter(|&c| lake[c] == lid).collect();
                (cells.iter().filter(|&&c| lm_on[c] != 0).count() as f32) < 0.5 * cells.len() as f32
            })
            .count();
        // residual sources: the ON world's final lake containing the origin, and the input lake containing it
        let (mut src_final, mut src_input): (HashMap<u32, usize>, HashMap<u32, usize>) = (HashMap::new(), HashMap::new());
        let mut no_origin = 0usize;
        for o in &origin {
            match o {
                Some(c) => {
                    *src_final.entry(lm_on[*c]).or_insert(0) += 1;
                    *src_input.entry(lake[*c]).or_insert(0) += 1;
                }
                None => no_origin += 1,
            }
        }
        eprintln!("   input lake | OFF moves cells | PRESENT (OFF) | base | in a basin_base depression | its OFF final lake a below-sea basin | residual cells it sources");
        let (mut ta, mut tab, mut tdep, mut tsea) = (0usize, 0usize, 0usize, 0usize);
        for &lid in &ids {
            let cells: Vec<usize> = (0..n).filter(|&c| lake[c] == lid).collect();
            let dshare = cells.iter().filter(|&&c| dep[c]).count() as f32 / cells.len() as f32;
            let (pres, fid) = present_in(&lm_off, lid);
            let mv = moved.get(&lid).copied().unwrap_or(0);
            let below = fid >= 1_000_000;
            if mv > 0 {
                ta += 1;
                if pres {
                    tab += 1;
                }
            }
            if dshare >= 0.5 {
                tdep += 1;
            }
            if below {
                tsea += 1;
            }
            eprintln!(
                "   {lid:>6} ({:>6} cells) | {mv:>7} | {} (final {fid}) | yes | {:>5.1} % | {} | {}",
                cells.len(),
                if pres { "YES" } else { "no " },
                100.0 * dshare,
                if below { "YES" } else { "no " },
                src_input.get(&lid).copied().unwrap_or(0)
            );
        }
        eprintln!("   final lake (ON, unbased) | area km² | below-sea basin | residual cells it sources (origin inside it)");
        let mut r2 = 0usize;
        for &fid in &unbased {
            let a = ls_on.iter().find(|x| x.0 == fid).map_or(0.0, |x| x.1);
            let sc = src_final.get(&fid).copied().unwrap_or(0);
            r2 += sc;
            eprintln!("   {fid:>8} | {a:>7.1} | {} | {sc}", if fid >= 1_000_000 { "YES" } else { "no " });
        }
        let m = resid.len().max(1) as f32;
        eprintln!(
            "   TOTALS · input lakes {} (bases {}) · moving χ under OFF {ta} · of them PRESENT {tab} · bases on a lake absent at the end (ON) {absent_on} · unbased final lakes (ON) {} · in a basin_base depression {tdep} · in a below-sea basin (OFF) {tsea} · OFF final lakes {} / ON {}",
            ids.len(),
            ids.len(),
            unbased.len(),
            ls_off.len(),
            ls_on.len()
        );
        eprintln!(
            "   R2 · residual (ON F132, A vs B, outside the input lakes) **{}** cells · source = one of the {} unbased final lakes **{r2}** (**{:.1} %**){} · source elsewhere {} · laid by a moved sample (no differing cell on the path) {no_origin}",
            resid.len(),
            unbased.len(),
            100.0 * r2 as f32 / m,
            if (r2 as f32) < 0.5 * m { " · STOP RULE (< 50 %)" } else { "" },
            resid.len() - r2 - no_origin
        );
        let mut rest: Vec<(u32, usize)> = src_final.iter().filter(|(k, _)| !unbased.contains(k)).map(|(k, v)| (*k, *v)).collect();
        rest.sort_by(|a, b| b.1.cmp(&a.1));
        eprintln!("   the rest, by the final lake containing the origin (0 = no lake): {:?}", rest.iter().take(10).collect::<Vec<_>>());
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n==========  end Finding 133-R . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 133-F — the remedy test of the extended `lake_base` (basin lakes are present lakes): the two
/// arbitrary resolutions of the lakes' flats, compared OUTSIDE the present lakes (the input lakes and the
/// `basin_base` closed depressions); any residual located (its origin classed), not corrected. Témoin and
/// B2 → A_c, carve on S1. Declared in `f133_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f133_f_remedy --nocapture
#[test]
#[ignore]
fn f133_f_remedy() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 133-F . the remedy test: basin lakes are present lakes  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let dc = C1DrainageConfig::default();
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    let lake = d.lake_map.clone();
    drop(d);
    let spill = ocean_flood(&bf);
    let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
    let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
    let present: Vec<bool> = (0..n).map(|k| lake[k] != 0 || dep[k]).collect();
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let f = flow.filled.data.clone();
    let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
    let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
    let arbitrary = |seed: u64| -> Vec<u8> {
        (0..n)
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
                let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                cand[(hsh % cand.len() as u64) as usize]
            })
            .collect()
    };
    let (arb_a, arb_b) = (arbitrary(1), arbitrary(2));
    drop((flow, f, g, needs));
    let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_a[c];
            }
        }
    };
    let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_b[c];
            }
        }
    };
    let id = |_: &GridF32, _: &[u32], _: &mut Vec<u8>| {};
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let base = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let ac = ValleyConstruction { a_min_km2: 0.1, ..base };
    eprintln!("   present lakes: input lakes {} cells ∪ basin_base closed depressions {} cells", lake.iter().filter(|&&l| l != 0).count(), dep.iter().filter(|&&b| b).count());
    for (label, world) in [("témoin C2/10 col", base), ("B2 → A_c", ac)] {
        let t = Instant::now();
        eprintln!("\n   ── {label} ──");
        for (state, lb) in [("ON F132 (InputLakes)", LakeBase::InputLakes), ("ON EXTENDED (InputLakesAndBasins)", LakeBase::InputLakesAndBasins)] {
            let vc = ValleyConstruction { lake_base: Some(lb), ..world };
            let (z_np, _) = carve(&s1, &skeleton(&s1, &vc, &ss, DOMAIN_KM), &vc, &ss);
            let (z_id, _) = carve(&s1, &skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&id)), &vc, &ss);
            let ctrl = (0..n).filter(|&k| z_np.data[k] != z_id.data[k]).count();
            drop((z_np, z_id));
            let ska = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pa));
            let skb = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pb));
            let (za, _) = carve(&s1, &ska, &vc, &ss);
            let (zb, _) = carve(&s1, &skb, &vc, &ss);
            let ch: Vec<usize> = (0..n).filter(|&k| !present[k] && za.data[k] != zb.data[k]).collect();
            let inside = (0..n).filter(|&k| present[k] && za.data[k] != zb.data[k]).count();
            let mut dz: Vec<f32> = ch.iter().map(|&k| (mm(&za, k) - mm(&zb, k)).abs()).collect();
            let (p50, p90, mx) = (q(&mut dz, 0.5), q(&mut dz, 0.9), q(&mut dz, 1.0));
            eprintln!(
                "   {state:<34} control {ctrl} · A vs B OUTSIDE the present lakes **{}** cells · |Δz| p50 {p50:.2} / p90 {p90:.2} / max {mx:.1} m · inside {inside}",
                ch.len()
            );
            if lb == LakeBase::InputLakesAndBasins && !ch.is_empty() {
                // the residual, located: walk each cell's path to the first differing cell, then to where χ / base agree
                let differs = |k: usize| ska.chi_m[k].to_bits() != skb.chi_m[k].to_bits() || ska.base_alt_m[k].to_bits() != skb.base_alt_m[k].to_bits();
                let mut cls: HashMap<&str, usize> = HashMap::new();
                let mut origins: HashMap<usize, usize> = HashMap::new();
                for &s in &ch {
                    let mut c = s;
                    let (mut seen, mut last) = (false, None);
                    for _ in 0..n {
                        if differs(c) {
                            seen = true;
                            last = Some(c);
                        } else if seen {
                            break;
                        }
                        let dd = ska.direction[c];
                        if dd == DIR_NONE || s1.data[c] <= SEA {
                            break;
                        }
                        c = nb(c, dd as usize);
                    }
                    let k = match last {
                        None => "no differing cell on its path (laid by a moved sample)",
                        Some(o) => {
                            *origins.entry(o).or_insert(0) += 1;
                            if lake[o] != 0 {
                                "origin in an input lake"
                            } else if dep[o] {
                                "origin in a closed depression"
                            } else if ska.direction[o] != skb.direction[o] {
                                "origin's pointer differs"
                            } else {
                                "origin elsewhere"
                            }
                        }
                    };
                    *cls.entry(k).or_insert(0) += 1;
                }
                let mut rows: Vec<_> = cls.into_iter().collect();
                rows.sort_by(|a, b| b.1.cmp(&a.1));
                for (k, c) in rows {
                    eprintln!("      the residual: {k:<56} {c}");
                }
                let mut top: Vec<_> = origins.into_iter().collect();
                top.sort_by(|a, b| b.1.cmp(&a.1));
                for &(o, c) in top.iter().take(5) {
                    eprintln!(
                        "      origin ({},{}) behind {c} · lake {} · dep {} · χ A/B {:.1} / {:.1} · base A/B {:.1} / {:.1} · pointer A/B {} / {}",
                        o % w,
                        o / w,
                        lake[o],
                        dep[o],
                        ska.chi_m[o],
                        skb.chi_m[o],
                        ska.base_alt_m[o],
                        skb.base_alt_m[o],
                        ska.direction[o],
                        skb.direction[o]
                    );
                }
            }
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n==========  end Finding 133-F (remedy) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 133-T/L/E — the témoin C2/10 col OFF, ON F132 (InputLakes) and ON extended (InputLakesAndBasins):
/// θ by stage with its CI, θ on the carved links and on the links left under their law (T); the new lakes (L); the
/// removal withheld, with its denominators (E); the extended variant's consistency; the table's columns; the cost of
/// the construction stage and of one full second pass. Declared in `f133_predictions.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f133_worlds --nocapture
#[test]
#[ignore]
fn f133_worlds() {
    use std::collections::VecDeque;
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, classify_lakes_water_balance};
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let cell_m = CELL_KM * 1000.0;
    eprintln!("\n==========  Finding 133 . the témoin OFF / ON F132 / ON extended  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &dcfg(), &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    // the input lakes and the closed depressions (basin_base's rule) of the construction's input
    let (lake, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    // closed-depression components
    let mut dep_id = vec![0u32; n];
    let mut ndep = 0u32;
    for c0 in 0..n {
        if !dep[c0] || dep_id[c0] != 0 {
            continue;
        }
        ndep += 1;
        let mut q = VecDeque::from([c0]);
        dep_id[c0] = ndep;
        while let Some(c) = q.pop_front() {
            for k in 0..8 {
                let m = nb(c, k);
                if dep[m] && dep_id[m] == 0 {
                    dep_id[m] = ndep;
                    q.push_back(m);
                }
            }
        }
    }
    // upstream of a present lake (the OFF skeleton's D8 path reaches an input lake or a closed depression)
    let mut up = vec![0u8; n]; // 1 input lake, 2 closed depression, 3 neither
    let mut path = Vec::new();
    for s in 0..n {
        if up[s] != 0 {
            continue;
        }
        path.clear();
        let mut c = s;
        let v;
        loop {
            if up[c] != 0 {
                v = up[c];
                break;
            }
            if lake[c] != 0 {
                v = 1;
                break;
            }
            if dep[c] {
                v = 2;
                break;
            }
            path.push(c);
            let dd = sk.direction[c];
            if dd == DIR_NONE || s1.data[c] <= SEA {
                v = 3;
                break;
            }
            c = nb(c, dd as usize);
        }
        for &p in &path {
            up[p] = v;
        }
    }
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let near = 2.0 / CELL_KM;
    let mut worlds: Vec<(String, Vec<f32>, Vec<u32>)> = Vec::new(); // (label, z metres, final lake map)
    for (label, lb) in [("OFF", None), ("ON F132 (InputLakes)", Some(LakeBase::InputLakes)), ("ON EXTENDED (InputLakesAndBasins)", Some(LakeBase::InputLakesAndBasins))] {
        let vc = ValleyConstruction { lake_base: lb, ..off };
        eprintln!("\n────────── témoin C2/10 col · {label} ──────────");
        // the construction stage alone (and its cost)
        let tc = Instant::now();
        let skv = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (built, mk) = carve(&s1, &skv, &vc, &ss);
        let t_con = tc.elapsed().as_secs_f64();
        let zb = metres(&built);
        drop(built);
        // T
        let law: Vec<f32> = (0..n).map(|k| if skv.chi_m[k].is_finite() { skv.floor_m(k, vc.age_k) } else { f32::NAN }).collect();
        let (tl, lol, hil, ml) = theta_laid_ci(&zb, &mk.carved, &sk);
        let under: Vec<bool> = (0..n).map(|k| sk.trunk[k] && sk.area_km2[k] >= 10.0 && !mk.carved[k] && zs1[k] < law[k]).collect();
        let nunder = trunk10.iter().filter(|&&k| under[k]).count();
        let (tu, lou, hiu, mu) = theta_laid_ci(&zb, &under, &sk);
        eprintln!(
            "   T · θ on the links the construction LOWERS at both ends **{tl:.3}** [{lol:.3}, {hil:.3}] ({ml}) · trunk cells UNCARVED under their own law **{nunder}** ({:.1} %) · θ on those links {tu:.3} [{lou:.3}, {hiu:.3}] ({mu}) · construction stage {t_con:.0} s",
            100.0 * nunder as f32 / trunk10.len() as f32
        );
        drop((skv, mk, law, under));
        // θ by stage
        let s2 = build_field_seed(Knobs { valley: Some(vc), no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
        let s2_lakes = c1_drainage_windowed(&s2, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM).lake_map;
        let bre2 = {
            let d = c1_drainage_windowed(&s2, None, &dcfg(), &ss, DOMAIN_KM);
            breach_monotone(&s2, &d.flow.filled, &d.lake_map, SEA, w, h)
        };
        let tw = Instant::now();
        let g = build_field_seed(Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
        let t_build = tw.elapsed().as_secs_f64();
        let mut th = Vec::new();
        for (stage, fld) in [("(i) as built", &s2), ("(ii) + breach", &bre2), ("(iii) the world", &g)] {
            let (tt, lo, hi, m) = theta_sa_ci(&metres(fld), &sk);
            th.push(format!("{stage} **{tt:.3}** [{lo:.3}, {hi:.3}] ({m})"));
        }
        eprintln!("   θ by stage: {}", th.join(" · "));
        drop((s2, bre2));
        // the table's columns (Finding 128-A's instruments)
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let cls = classes(&pre, &skp, &mkp);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw, skp, mkp));
        let co = coast(&g, &ss, &dwall, near);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let r8w = r8_terrain(&g, &cls, 16);
        let ta = Instant::now();
        let dd = c1_drainage_windowed(&g, None, &dcfg(), &ss, DOMAIN_KM);
        let bre = breach_monotone(&g, &dd.flow.filled, &dd.lake_map, SEA, w, h);
        let t_pre_assembly = ta.elapsed().as_secs_f64();
        let cl_e = c1_climate_placed(&g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(&g, &dcfg(), &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(&g, &bre, &pre, DELIVERED_P50_M, &ss, &dcfg(), cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = Vec::new();
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons.push((floor % w, floor / w));
            }
        }
        drop(fill_del);
        let tb = Instant::now();
        let cl_b = c1_climate_placed(&bre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_b = DrainageClimate { precip_internal: &cl_b.precipitation, temperature: &cl_b.temperature };
        let dr = assemble_hd_drainage(&bre, &dc_b, Some(dd), &dcfg(), &ss, DOMAIN_KM, GEO_RATIO, None, false).drainage;
        let t_assembly = tb.elapsed().as_secs_f64() + t_pre_assembly;
        let segs = &dr.rivers.segments;
        let wc: Vec<usize> = (0..segs.len()).filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2).collect();
        let pts: Vec<&[(u32, u32)]> = wc.iter().map(|&i| own(&segs[i])).collect();
        let zg = metres(&g);
        // E: the removal relative to the construction's input
        let removal = |sel: &dyn Fn(usize) -> bool| -> f64 {
            (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((zs1[k] - zg[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum()
        };
        let r_all = removal(&|_| true);
        eprintln!(
            "   SUMMARY · canyons **{}** {:?} · spurs near a coastal wall {rn:.4} /km · R8 terrain {r8w:.4} · R8 network (chord 8) {:.4} · lakes **{}** ({:.0} km²) · removal (S1 − world, land) **{r_all:.1} km³** · Finding 38 holds",
            canyons.len(),
            canyons,
            r8_chords(&pts, 8),
            dr.lakes.len(),
            dr.lakes.iter().map(|l| l.area_km2).sum::<f32>()
        );
        eprintln!("   COST · construction stage (skeleton + carve) {t_con:.0} s · a full second pass: the world's build {t_build:.0} s + drainage, breach, climate, assembly {t_assembly:.0} s = **{:.0} s** (wall-clock)", t_build + t_assembly);
        // consistency (the ON variants): present set per variant
        if let Some(lbv) = lb {
            let pres = |k: usize| lake[k] != 0 || (lbv == LakeBase::InputLakesAndBasins && dep[k]);
            let mut absent = Vec::new();
            let mut ids: Vec<u32> = lake.iter().copied().filter(|&l| l != 0).collect();
            ids.sort_unstable();
            ids.dedup();
            for &lid in &ids {
                let cells: Vec<usize> = (0..n).filter(|&c| lake[c] == lid).collect();
                if (cells.iter().filter(|&&c| dr.lake_map[c] != 0).count() as f32) < 0.5 * cells.len() as f32 {
                    absent.push(format!("input {lid}"));
                }
            }
            if lbv == LakeBase::InputLakesAndBasins {
                for did in 1..=ndep {
                    let cells: Vec<usize> = (0..n).filter(|&c| dep_id[c] == did).collect();
                    if cells.len() >= 50 && (cells.iter().filter(|&&c| dr.lake_map[c] != 0).count() as f32) < 0.5 * cells.len() as f32 {
                        absent.push(format!("depression {did} ({} cells)", cells.len()));
                    }
                }
            }
            let mut unbased = Vec::new();
            for l in &dr.lakes {
                let cells: Vec<usize> = (0..n).filter(|&c| dr.lake_map[c] == l.base.id).collect();
                let cov = cells.iter().filter(|&&c| pres(c)).count();
                if (cov as f32) < 0.5 * cells.len().max(1) as f32 {
                    unbased.push((l.base.id, l.area_km2.round() as i32, (100.0 * cov as f32 / cells.len().max(1) as f32).round() as i32));
                }
            }
            eprintln!(
                "   CONSISTENCY · bases on a lake ABSENT at the end **{}** {:?} · final lakes WITHOUT a base **{}** {:?} · closed-depression components {ndep}",
                absent.len(),
                absent,
                unbased.len(),
                unbased
            );
        }
        // L: the new lakes against the OFF world
        if let Some((_, _, lm_off)) = worlds.first() {
            let verdicts = classify_lakes_water_balance(&bre, &dr.flow, &dc_b, cell_km2, &dr.lakes, &dr.lake_map, None, w, h);
            let mut rows = Vec::new();
            for l in &dr.lakes {
                let cells: Vec<usize> = (0..n).filter(|&c| dr.lake_map[c] == l.base.id).collect();
                if cells.is_empty() || (cells.iter().filter(|&&c| lm_off[c] != 0).count() as f32) >= 0.5 * cells.len() as f32 {
                    continue;
                }
                let (cx, cy) = (cells.iter().map(|&c| (c % w) as f64).sum::<f64>() / cells.len() as f64, cells.iter().map(|&c| (c / w) as f64).sum::<f64>() / cells.len() as f64);
                let is_input = cells.iter().filter(|&&c| lake[c] != 0).count() as f32 >= 0.5 * cells.len() as f32;
                let is_dep = cells.iter().filter(|&&c| dep[c]).count() as f32 >= 0.5 * cells.len() as f32;
                // upstream of a based lake: walk the OFF skeleton from the outlet, outside the lake itself
                let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
                let (mut c, mut steps, mut hit) = (o, 0usize, "none");
                for _ in 0..n {
                    if dr.lake_map[c] != l.base.id {
                        if lake[c] != 0 {
                            hit = "an input lake";
                            break;
                        }
                        if dep[c] {
                            hit = "a closed depression";
                            break;
                        }
                    }
                    let d8 = sk.direction[c];
                    if d8 == DIR_NONE || s1.data[c] <= SEA {
                        break;
                    }
                    c = nb(c, d8 as usize);
                    steps += 1;
                }
                let s2cov = cells.iter().filter(|&&c| s2_lakes[c] != 0).count();
                let stage = if s2cov as f32 >= 0.5 * cells.len() as f32 { "construction" } else { "light pass" };
                let a_eq = verdicts.iter().find(|v| v.id == l.base.id).map_or(f32::NAN, |v| v.a_eq_km2);
                rows.push(format!(
                    "      new lake {} at ({cx:.0},{cy:.0}) · {:.1} km² · depth {:.1} m · a_eq {a_eq:.1} · itself an input lake {} · a closed depression {} · downstream it meets {hit}{} · STAGE {stage} (at that stage {:.1} km²)",
                    l.base.id,
                    l.area_km2,
                    l.depth_m,
                    if is_input { "YES" } else { "no" },
                    if is_dep { "YES" } else { "no" },
                    if hit != "none" { format!(" at {:.0} m", steps as f32 * cell_m) } else { String::new() },
                    s2cov as f32 * cell_km2
                ));
            }
            eprintln!("   L · {} new lakes against OFF:", rows.len());
            for r in rows {
                eprintln!("{r}");
            }
        }
        worlds.push((label.to_string(), zg, dr.lake_map.clone()));
    }
    // E: the removal withheld, with its denominators
    let (z_off, z_132, z_ext) = (&worlds[0].1, &worlds[1].1, &worlds[2].1);
    let rem = |z: &[f32], sel: &dyn Fn(usize) -> bool| -> f64 {
        (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((zs1[k] - z[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum()
    };
    let r_off = rem(z_off, &|_| true);
    for (name, z, ext) in [("ON F132", z_132, false), ("ON EXTENDED", z_ext, true)] {
        let upstream = |k: usize| lake[k] != 0 || up[k] == 1 || (ext && (dep[k] || up[k] == 2));
        let r_off_up = rem(z_off, &upstream);
        let r_on = rem(z, &|_| true);
        let withheld = r_off - r_on;
        let parts: Vec<String> = [
            ("upstream of an input lake", &(|k: usize| up[k] == 1) as &dyn Fn(usize) -> bool),
            ("inside an input lake", &(|k: usize| lake[k] != 0) as &dyn Fn(usize) -> bool),
            ("in or upstream of a closed depression", &(|k: usize| dep[k] || up[k] == 2) as &dyn Fn(usize) -> bool),
            ("elsewhere", &(|k: usize| lake[k] == 0 && !dep[k] && up[k] == 3) as &dyn Fn(usize) -> bool),
        ]
        .iter()
        .map(|(p, sel)| format!("{p} {:.1} km³", rem(z_off, *sel) - rem(z, *sel)))
        .collect();
        eprintln!(
            "   E · {name}: removal withheld **{withheld:.1} km³** = **{:.1} %** of the témoin's total removal ({r_off:.1} km³) = {:.1} % of the témoin's removal in / upstream of its present lakes ({r_off_up:.1} km³) · where: {}",
            100.0 * withheld / r_off,
            100.0 * withheld / r_off_up.max(1e-9),
            parts.join(" · ")
        );
    }
    let _ = q;
    eprintln!("\n==========  end Finding 133 (worlds) . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 134-O — **the ocean moves.** The viz world (`new(k, Some(0.1))`, the guarded "C2 /10 col"), OFF vs
/// ON extended, stage by stage (rule 12): construction (+ rims), + light pass, + droplet erosion, + bathymetry.
/// The first stage where an OCEAN cell differs; the bathymetry's own mean depth on each world; the ratio
/// z_ON / z_OFF on the sea cells against `m_OFF / m_ON`. Declared in `f134_declared.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f134_ocean --nocapture
#[test]
#[ignore]
fn f134_ocean() {
    use ymir_core::tectonics_c1::valley_construction::LakeBase;
    use ymir_core::terrain::bathymetry::BathymetryProfile;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let depth_per_norm = 2.0 * 1.13 * ss.depth_scale_m as f32;
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 134-O . the ocean moves (viz world, OFF / ON extended)  ==========");
    eprintln!("   depth_per_norm {depth_per_norm:.3} m (the code's) · c1_altitude_norm_to_metres slope {n2m:.3} m");
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let knobs = |vc: ValleyConstruction, stage: usize| -> Knobs {
        let b = Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
        match stage {
            0 => Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..b },
            1 => Knobs { erosion_off: true, bathymetry_off: true, ..b },
            2 => Knobs { bathymetry_off: true, ..b },
            _ => b,
        }
    };
    let names = ["(a) construction (+ rims)", "(b) + light pass", "(c) + droplet erosion (pre-bathymetry)", "(d) + bathymetry = the eroded world"];
    let mean_depth = |g: &GridF32| -> (f64, u64) {
        let (mut s, mut c) = (0f64, 0u64);
        for &v in &g.data {
            if v <= SEA {
                s += ((SEA - v) * depth_per_norm) as f64;
                c += 1;
            }
        }
        (s / c as f64, c)
    };
    let mut pre: Option<(GridF32, GridF32)> = None;
    let mut first_ocean: Option<&str> = None;
    for (st, name) in names.iter().enumerate() {
        let t = Instant::now();
        let a = build_field_seed(knobs(off, st), PSEED);
        let b = build_field_seed(knobs(on, st), PSEED);
        let (w, h) = (a.width, a.height);
        let n = w * h;
        let wc = water_class(&a, SEA);
        let (mut d_all, mut d_ocean, mut d_inland, mut d_land, mut flip) = (0usize, 0usize, 0usize, 0usize, 0usize);
        for k in 0..n {
            let sa = a.data[k] <= SEA;
            let sb = b.data[k] <= SEA;
            if sa != sb {
                flip += 1;
            }
            if a.data[k].to_bits() == b.data[k].to_bits() {
                continue;
            }
            d_all += 1;
            if sa && sb {
                if wc[k] == WATER_CLASS_OCEAN {
                    d_ocean += 1;
                } else {
                    d_inland += 1;
                }
            } else {
                d_land += 1;
            }
        }
        let n_ocean = (0..n).filter(|&k| wc[k] == WATER_CLASS_OCEAN).count();
        if d_ocean > 0 && first_ocean.is_none() {
            first_ocean = Some(name);
        }
        eprintln!(
            "   {name:<42} cells differing **{d_all}** · OCEAN cells (both sea, edge-connected) **{d_ocean}** of {n_ocean} · INLAND sea cells {d_inland} · land / flipped {d_land} · sea↔land flips {flip} ({:.0} s)",
            t.elapsed().as_secs_f64()
        );
        if st == 2 {
            let (ma, ca) = mean_depth(&a);
            let (mb, cb) = mean_depth(&b);
            // where the mean moved: cells sea in both with a different depth, and the membership flips
            let (mut s_ocean, mut s_inland) = (0f64, 0f64);
            for k in 0..n {
                if a.data[k] <= SEA && b.data[k] <= SEA && a.data[k] != b.data[k] {
                    let dd = ((a.data[k] - b.data[k]) * depth_per_norm) as f64; // depth ON − depth OFF
                    if wc[k] == WATER_CLASS_OCEAN {
                        s_ocean += dd;
                    } else {
                        s_inland += dd;
                    }
                }
            }
            eprintln!(
                "   MEAN DEPTH (the bathymetry's own, over z ≤ sea at stage c): OFF **{ma:.4} m** over {ca} cells · ON **{mb:.4} m** over {cb} cells · m_OFF / m_ON **{:.7}** · Σ(depth ON − depth OFF) on OCEAN cells {:.3e} m, on INLAND sea cells {:.3e} m",
                ma / mb,
                s_ocean,
                s_inland
            );
            pre = Some((a, b));
        } else if st == 3 {
            let (pa, pb) = pre.take().expect("stage c");
            let p = BathymetryProfile::default();
            let floor_depth = SEA * depth_per_norm;
            let min_depth = p.shelf_min_depth_m.max(1.0);
            let (ma, _) = mean_depth(&pa);
            let (mb, _) = mean_depth(&pb);
            let pred = (ma / mb) as f32;
            let mut r: Vec<f32> = Vec::new();
            let (mut c_eq, mut c_eq_ok, mut c_ne, mut c_ne_ok, mut clamp, mut ocean_n) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
            for k in 0..n {
                if !(a.data[k] <= SEA && b.data[k] <= SEA) {
                    continue;
                }
                let (za, zb) = (c1_altitude_norm_to_metres(a.data[k], &ss), c1_altitude_norm_to_metres(b.data[k], &ss));
                if za == 0.0 {
                    continue;
                }
                if wc[k] == WATER_CLASS_OCEAN {
                    ocean_n += 1;
                }
                let rr = zb / za;
                r.push(rr);
                let (da, db) = ((SEA - a.data[k]) * depth_per_norm, (SEA - b.data[k]) * depth_per_norm);
                let clamped = (da - min_depth).abs() < 1e-3 || (db - min_depth).abs() < 1e-3 || (da - floor_depth).abs() < 1e-3 || (db - floor_depth).abs() < 1e-3;
                let ok = ((rr / pred) - 1.0).abs() <= 1e-5;
                if clamped {
                    clamp += 1;
                } else if pa.data[k].to_bits() == pb.data[k].to_bits() {
                    c_eq += 1;
                    c_eq_ok += ok as usize;
                } else {
                    c_ne += 1;
                    c_ne_ok += ok as usize;
                }
            }
            let m = r.len();
            let mut rs = r.clone();
            rs.sort_by(f32::total_cmp);
            let med = rs[m / 2];
            let within = r.iter().filter(|&&x| (x - med).abs() <= 1e-3).count();
            let q = |p: f64| rs[((m - 1) as f64 * p) as usize];
            eprintln!(
                "   RATIO z_ON / z_OFF on the {m} sea cells (both ≤ sea; {ocean_n} OCEAN): median **{med:.7}** · p0.1 {:.5} · p1 {:.5} · p99 {:.5} · p99.9 {:.5} · within 1e-3 of the median **{within} ({:.3} %)** · prediction m_OFF / m_ON {pred:.7}",
                q(0.001),
                q(0.01),
                q(0.99),
                q(0.999),
                100.0 * within as f64 / m as f64
            );
            eprintln!(
                "   ATTRIBUTION · pre-bathymetry depth c EQUAL: {c_eq} cells, of which **{c_eq_ok} ({:.3} %)** at the prediction within 1e-5 · c DIFFERENT: {c_ne} cells ({c_ne_ok} at the prediction) · clamped (shelf_min {min_depth:.1} m or floor {floor_depth:.1} m): {clamp}",
                100.0 * c_eq_ok as f64 / c_eq.max(1) as f64
            );
        }
    }
    eprintln!("   FIRST STAGE where an OCEAN cell differs: **{}**", first_ocean.unwrap_or("none"));
    eprintln!("\n==========  end Finding 134-O . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 134-D — **F133v's p90 of 1 177 m upstream of lake 2**, and lakes 11 and 1 beside it, on the viz
/// world reproduced by the bench's `run_hd` tail (`common::viz_hd_lakes_on`). Per lake: the field of Δz
/// (conditioned = with the pre-breach lakes' surfaces, and the eroded field without water), its sign, and its
/// split into the lake's ON footprint and the upstream cells (F133v's attribution), each into carved by the
/// construction vs kept at the terrain; the lake's depth, level and col beside it; a map. Declared in
/// `f134_declared.md` before the run. Maps to `F134_DIR` if set.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f134_d --nocapture
#[test]
#[ignore]
fn f134_d() {
    use common::{build_world, viz_hd_lakes_on};
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 134-D . the Δz upstream of the based lakes 2, 11, 1 (viz world)  ==========");
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let (eo, co, lm_off, lakes_off, dir, gw) = {
        let wd = build_world(kn(off), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(off), PSEED, 45.0, 40.0);
        eprintln!("   OFF · lakes {} · lake key {}", v.drainage.lakes.len(), v.digest);
        (metres(&wd.heightmap), metres(&v.conditioned), v.drainage.lake_map.clone(), v.drainage.lakes.clone(), v.drainage.flow.direction.clone(), wd.heightmap.width)
    };
    let (en, cn, lm_on, lakes_on) = {
        let wd = build_world(kn(on), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(on), PSEED, 45.0, 40.0);
        eprintln!("   ON  · lakes {} · lake key {}", v.drainage.lakes.len(), v.digest);
        (metres(&wd.heightmap), metres(&v.conditioned), v.drainage.lake_map.clone(), v.drainage.lakes.clone())
    };
    let _ = lakes_off;
    let n = eo.len();
    let (w, h) = (gw, n / gw);
    let changed_c = (0..n).filter(|&k| co[k] != cn[k]).count();
    eprintln!("   conditioned cells differing OFF vs ON: **{changed_c}** (F133v: 44 268 342)");
    // the construction input, the carve masks, the col
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (_, mk_off) = carve(&s1, &skeleton(&s1, &off, &ss, DOMAIN_KM), &off, &ss);
    let (_, mk_on) = carve(&s1, &skeleton(&s1, &on, &ss, DOMAIN_KM), &on, &ss);
    let spill: Vec<f32> = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        ocean_flood(&bf).iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect()
    };
    let zs1 = metres(&s1);
    drop(s1);
    // F133v's attribution: the OFF drainage's D8 path to the first ON final lake (no torus, as F133v)
    let mut first = vec![u32::MAX; n];
    let mut path = Vec::new();
    for s in 0..n {
        if first[s] != u32::MAX {
            continue;
        }
        path.clear();
        let mut c = s;
        let v;
        loop {
            if first[c] != u32::MAX {
                v = first[c];
                break;
            }
            if lm_on[c] != 0 {
                v = lm_on[c];
                break;
            }
            path.push(c);
            let d = dir[c];
            if d == DIR_NONE {
                v = 0;
                break;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
                v = 0;
                break;
            }
            c = y as usize * w + x as usize;
        }
        for &p in &path {
            first[p] = v;
        }
    }
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let dir_out = std::env::var("F134_DIR").ok();
    for (lid, centre) in [(2u32, (4019.0f64, 2839.0f64)), (11, (1429.0, 4287.0)), (1, (3491.0, 2500.0))] {
        let cells: Vec<usize> = (0..n).filter(|&k| lm_on[k] == lid).collect();
        let (cx, cy) = (
            cells.iter().map(|&k| (k % w) as f64).sum::<f64>() / cells.len().max(1) as f64,
            cells.iter().map(|&k| (k / w) as f64).sum::<f64>() / cells.len().max(1) as f64,
        );
        let l = lakes_on.iter().find(|l| l.base.id == lid).expect("ON lake");
        let col = cells.iter().map(|&k| spill[k]).fold(f32::NEG_INFINITY, f32::max);
        let mut zoff_in: Vec<f32> = cells.iter().map(|&k| co[k]).collect();
        let mut zon_in: Vec<f32> = cells.iter().map(|&k| cn[k]).collect();
        let mut zs1_in: Vec<f32> = cells.iter().map(|&k| zs1[k]).collect();
        let in_off_lake = cells.iter().filter(|&&k| lm_off[k] != 0).count();
        eprintln!(
            "\n   ── lake {lid} · centre ({cx:.0}, {cy:.0}) (F133v ({:.0}, {:.0})) · {:?} · {:.1} km² · ON level **{:.1} m** · depth **{:.1} m** · col (ocean_flood spill on S1, max over the footprint) **{col:.1} m** · footprint under an OFF lake {in_off_lake} of {} · conditioned z in the footprint p50 OFF {:.1} / ON {:.1} m · S1 p50 {:.1} m",
            centre.0,
            centre.1,
            l.lake_type,
            l.area_km2,
            l.level_m,
            l.depth_m,
            cells.len(),
            q(&mut zoff_in, 0.5),
            q(&mut zon_in, 0.5),
            q(&mut zs1_in, 0.5)
        );
        for (fname, zo, zn) in [("CONDITIONED (pre-breach lakes at their surface)", &co, &cn), ("ERODED (pre-breach, no water)", &eo, &en)] {
            eprintln!("      field {fname}:");
            let groups: [(&str, &dyn Fn(usize) -> bool); 6] = [
                ("footprint · carved", &|k| lm_on[k] == lid && (mk_off.carved[k] || mk_on.carved[k])),
                ("footprint · kept at the terrain", &|k| lm_on[k] == lid && !(mk_off.carved[k] || mk_on.carved[k])),
                ("upstream · carved (OFF and ON)", &|k| lm_on[k] == 0 && first[k] == lid && mk_off.carved[k] && mk_on.carved[k]),
                ("upstream · carved OFF only", &|k| lm_on[k] == 0 && first[k] == lid && mk_off.carved[k] && !mk_on.carved[k]),
                ("upstream · carved ON only", &|k| lm_on[k] == 0 && first[k] == lid && !mk_off.carved[k] && mk_on.carved[k]),
                ("upstream · kept at the terrain", &|k| lm_on[k] == 0 && first[k] == lid && !mk_off.carved[k] && !mk_on.carved[k]),
            ];
            for (gname, sel) in groups {
                let ks: Vec<usize> = (0..n).filter(|&k| sel(k)).collect();
                let dz: Vec<f32> = ks.iter().map(|&k| zn[k] - zo[k]).filter(|&d| d != 0.0).collect();
                let pos = dz.iter().filter(|&&d| d > 0.0).count();
                let mut a: Vec<f32> = dz.iter().map(|d| d.abs()).collect();
                let s: f64 = a.iter().map(|&d| d as f64 * cell_km2 as f64 * 1e-3).sum();
                let in_off = ks.iter().filter(|&&k| lm_off[k] != 0 && zn[k] != zo[k]).count();
                let (p50, p90) = (q(&mut a, 0.5), q(&mut a, 0.9));
                eprintln!(
                    "         {gname:<34} cells {:>7} · Δz ≠ 0 {:>7} · Δz > 0 **{:.1} %** · |Δz| p50 {p50:.1} / p90 **{p90:.1} m** · Σ|Δz|·A {s:.2} km³ · of the Δz ≠ 0, under an OFF lake {in_off}",
                    ks.len(),
                    dz.len(),
                    100.0 * pos as f32 / dz.len().max(1) as f32
                );
            }
            let mut up_all: Vec<f32> = (0..n).filter(|&k| lm_on[k] == 0 && first[k] == lid && zn[k] != zo[k]).map(|k| (zn[k] - zo[k]).abs()).collect();
            let m = up_all.len();
            eprintln!("         upstream, all (F133v's population): {m} cells with Δz ≠ 0 · p90 |Δz| **{:.1} m**", q(&mut up_all, 0.9));
        }
        // where the large upstream ones sit (conditioned)
        let big: Vec<usize> = (0..n).filter(|&k| lm_on[k] == 0 && first[k] == lid && (cn[k] - co[k]).abs() > 500.0).collect();
        if !big.is_empty() {
            let mut d_lake: Vec<f32> = Vec::new();
            for &k in &big {
                let (x, y) = ((k % w) as f64, (k / w) as f64);
                d_lake.push((((x - cx).powi(2) + (y - cy).powi(2)).sqrt() as f32) * CELL_KM);
            }
            let in_off = big.iter().filter(|&&k| lm_off[k] != 0).count();
            let sea_off = big.iter().filter(|&&k| co[k] <= 0.0).count();
            let sea_on = big.iter().filter(|&&k| cn[k] <= 0.0).count();
            let mut zo_b: Vec<f32> = big.iter().map(|&k| co[k]).collect();
            let mut zn_b: Vec<f32> = big.iter().map(|&k| cn[k]).collect();
            let mut zs_b: Vec<f32> = big.iter().map(|&k| zs1[k]).collect();
            eprintln!(
                "      |Δz| > 500 m upstream (conditioned): {} cells · under an OFF lake {in_off} · z ≤ 0 OFF {sea_off} / ON {sea_on} · z OFF p50 {:.1} m · z ON p50 {:.1} m · S1 p50 {:.1} m · distance to the lake's centre p50 {:.1} km",
                big.len(),
                q(&mut zo_b, 0.5),
                q(&mut zn_b, 0.5),
                q(&mut zs_b, 0.5),
                q(&mut d_lake, 0.5)
            );
        }
        // the map: Δz conditioned, ON outline white, OFF lakes black
        if let Some(d) = &dir_out {
            let half = 307i64;
            let side = (2 * half + 1) as u32;
            let mut img = image::RgbImage::new(side, side);
            for yy in 0..side as i64 {
                for xx in 0..side as i64 {
                    let (x, y) = (cx.round() as i64 - half + xx, cy.round() as i64 - half + yy);
                    if x < 0 || y < 0 || x >= w as i64 || y >= h as i64 {
                        continue;
                    }
                    let k = y as usize * w + x as usize;
                    let dz = (cn[k] - co[k]).clamp(-1000.0, 1000.0) / 1000.0;
                    let mut px = if dz >= 0.0 {
                        [255, (255.0 * (1.0 - dz)) as u8, (255.0 * (1.0 - dz)) as u8]
                    } else {
                        [(255.0 * (1.0 + dz)) as u8, (255.0 * (1.0 + dz)) as u8, 255]
                    };
                    let edge = |lm: &[u32]| -> bool {
                        lm[k] != 0 && [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)].iter().any(|&(ox, oy)| {
                            let (a, b) = (x + ox, y + oy);
                            a >= 0 && b >= 0 && a < w as i64 && b < h as i64 && lm[b as usize * w + a as usize] != lm[k]
                        })
                    };
                    if edge(&lm_off) {
                        px = [0, 0, 0];
                    }
                    if edge(&lm_on) {
                        px = if lm_on[k] == lid { [0, 160, 0] } else { [90, 90, 90] };
                    }
                    img.put_pixel(xx as u32, yy as u32, image::Rgb(px));
                }
            }
            let p = std::path::Path::new(d).join(format!("d_lake{lid}_dz.png"));
            img.save(&p).expect("png");
            eprintln!("      map: {} (red Δz > 0, blue < 0, ±1000 m; lake {lid} outline green, other ON lakes grey, OFF lakes black)", p.display());
        }
    }
    eprintln!("\n==========  end Finding 134-D . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 134-T / U — on the benches' témoin (F133's world): θ on the carved links by population (upstream of
/// a based lake, downstream, untouched) with its CI; the carved trunk cells more than 1 m off their law, located
/// against the based shores; the trunk cells newly under their law in ON, located (distance to a based shore,
/// law − terrain, the under-law run they sit in). Declared in `f134_declared.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f134_tu --nocapture
#[test]
#[ignore]
fn f134_tu() {
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 134-T / U . θ by population, the off-law cells, the new under-law cells (témoin)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let off = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let t10 = |k: usize| sk.trunk[k] && sk.area_km2[k] >= 10.0;
    let ntrunk = (0..n).filter(|&k| t10(k)).count();
    let (lake, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let recv = |c: usize| -> Option<usize> {
        let d = sk.direction[c];
        if d == DIR_NONE || s1.data[c] <= SEA { None } else { Some(nb(c, d as usize)) }
    };
    // populations per based set: 1 upstream, 2 downstream, 3 untouched, 4 inside (upstream wins over downstream)
    let pops = |ext: bool| -> (Vec<bool>, Vec<u8>) {
        let based: Vec<bool> = (0..n).map(|k| lake[k] != 0 || (ext && dep[k])).collect();
        let mut up = vec![0u8; n]; // 1 reaches the based set, 2 does not
        let mut path = Vec::new();
        for s in 0..n {
            if up[s] != 0 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v;
            loop {
                if up[c] != 0 {
                    v = up[c];
                    break;
                }
                if based[c] && c != s {
                    v = 1;
                    break;
                }
                path.push(c);
                match recv(c) {
                    Some(r) => c = r,
                    None => {
                        v = 2;
                        break;
                    }
                }
            }
            for &p in &path {
                up[p] = v;
            }
        }
        // downstream: Kahn over the skeleton's receivers
        let mut indeg = vec![0u8; n];
        for c in 0..n {
            if let Some(r) = recv(c) {
                indeg[r] = indeg[r].saturating_add(1);
            }
        }
        let mut down = vec![false; n];
        let mut stack: Vec<usize> = (0..n).filter(|&c| indeg[c] == 0).collect();
        while let Some(c) = stack.pop() {
            if let Some(r) = recv(c) {
                if down[c] || based[c] {
                    down[r] = true;
                }
                indeg[r] -= 1;
                if indeg[r] == 0 {
                    stack.push(r);
                }
            }
        }
        let pop: Vec<u8> = (0..n)
            .map(|k| if based[k] { 4 } else if up[k] == 1 { 1 } else if down[k] { 2 } else { 3 })
            .collect();
        (based, pop)
    };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    const PNAME: [&str; 5] = ["", "upstream of a based lake", "downstream of a based lake", "untouched basins", "inside"];
    let sets = [("F132 set (input lakes)", pops(false)), ("extended set (input lakes ∪ closed depressions)", pops(true))];
    let dists: Vec<Vec<u16>> = sets.iter().map(|(_, (b, _))| dist_from(b, w, h)).collect();
    for (si, (sname, (_, pop))) in sets.iter().enumerate() {
        let cnt: Vec<usize> = (1..=4u8).map(|p| (0..n).filter(|&k| t10(k) && pop[k] == p).count()).collect();
        eprintln!("   {sname}: trunk ≥ 10 km² cells by population: upstream {} · downstream {} · untouched {} · inside {} (of {ntrunk})", cnt[0], cnt[1], cnt[2], cnt[3]);
        let _ = si;
    }
    let mut under_off: Vec<bool> = Vec::new();
    for (label, lb) in [("OFF", None), ("ON F132 (InputLakes)", Some(LakeBase::InputLakes)), ("ON EXTENDED (InputLakesAndBasins)", Some(LakeBase::InputLakesAndBasins))] {
        let t = Instant::now();
        let vc = ValleyConstruction { lake_base: lb, ..off };
        let skv = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (built, mk) = carve(&s1, &skv, &vc, &ss);
        let zb = metres(&built);
        drop(built);
        let law: Vec<f32> = (0..n).map(|k| if skv.chi_m[k].is_finite() { skv.floor_m(k, vc.age_k) } else { f32::NAN }).collect();
        drop(skv);
        let (ta, la, ha, ma) = theta_laid_ci(&zb, &mk.carved, &sk);
        eprintln!("\n   ── {label} ── θ on the carved links, all **{ta:.3}** [{la:.3}, {ha:.3}] ({ma} links) ({:.0} s)", t.elapsed().as_secs_f64());
        // which based set reads this variant: F132 → set 0, extended → set 1, OFF → both
        let which: Vec<usize> = match lb {
            None => vec![0, 1],
            Some(LakeBase::InputLakes) => vec![0],
            Some(LakeBase::InputLakesAndBasins) => vec![1],
        };
        for &si in &which {
            let (sname, (_, pop)) = &sets[si];
            let mut row = Vec::new();
            for p in 1..=3u8 {
                let mask: Vec<bool> = (0..n).map(|k| mk.carved[k] && pop[k] == p).collect();
                let (tp, lp, hp, mp) = theta_laid_ci(&zb, &mask, &sk);
                row.push(format!("{} **{tp:.3}** [{lp:.3}, {hp:.3}] ({mp})", PNAME[p as usize]));
            }
            eprintln!("   T · θ by population ({sname}): {}", row.join(" · "));
            // the carved trunk cells more than 1 m off their law
            let dev: Vec<usize> = (0..n).filter(|&k| t10(k) && mk.carved[k] && law[k].is_finite() && (zb[k] - law[k]).abs() > 1.0).collect();
            let above = dev.iter().filter(|&&k| zb[k] > law[k]).count();
            let carved_t = (0..n).filter(|&k| t10(k) && mk.carved[k]).count();
            let mut dkm: Vec<f32> = dev.iter().map(|&k| dists[si][k] as f32 * CELL_KM).collect();
            let near2 = dkm.iter().filter(|&&d| d <= 2.0).count();
            let by: Vec<String> = (1..=4u8).map(|p| format!("{} {}", PNAME[p as usize], dev.iter().filter(|&&k| pop[k] == p).count())).collect();
            let mut off_m: Vec<f32> = dev.iter().map(|&k| (zb[k] - law[k]).abs()).collect();
            eprintln!(
                "   T · carved trunk cells > 1 m off their law: **{}** of {carved_t} carved trunk cells · above the law {above} · |z − law| p50 {:.1} / p90 {:.1} m · distance to a based shore p50 **{:.2} km** / p90 {:.2} km · within 2 km **{:.1} %** · by population: {}",
                dev.len(),
                q(&mut off_m, 0.5),
                q(&mut off_m, 0.9),
                q(&mut dkm, 0.5),
                q(&mut dkm, 0.9),
                100.0 * near2 as f32 / dev.len().max(1) as f32,
                by.join(" · ")
            );
        }
        // U
        let under: Vec<bool> = (0..n).map(|k| t10(k) && !mk.carved[k] && law[k].is_finite() && zs1[k] < law[k]).collect();
        let nu = under.iter().filter(|&&b| b).count();
        eprintln!("   U · trunk cells uncarved under their law: **{nu}** ({:.1} % of {ntrunk})", 100.0 * nu as f32 / ntrunk as f32);
        if lb.is_none() {
            under_off = under;
            continue;
        }
        let si = which[0];
        let (based, pop) = (&sets[si].1.0, &sets[si].1.1);
        let added: Vec<usize> = (0..n).filter(|&k| under[k] && !under_off[k]).collect();
        let gone = (0..n).filter(|&k| !under[k] && under_off[k]).count();
        let mut dkm: Vec<f32> = added.iter().map(|&k| dists[si][k] as f32 * CELL_KM).collect();
        let near2 = dkm.iter().filter(|&&d| d <= 2.0).count();
        let mut lt: Vec<f32> = added.iter().map(|&k| law[k] - zs1[k]).collect();
        let by: Vec<String> = (1..=4u8).map(|p| format!("{} {}", PNAME[p as usize], added.iter().filter(|&&k| pop[k] == p).count())).collect();
        // the under-law run each added cell sits in: walked downstream while under the law
        let mut rem = vec![u32::MAX; n]; // remaining run length (cells) downstream, memoised
        let mut end = vec![0u8; n]; // 1 based shore, 2 sea / no receiver, 3 leaves the law on the trunk
        let mut stackp = Vec::new();
        for &s in &added {
            stackp.clear();
            let mut c = s;
            let mut base_len = 0u32;
            let base_end: u8;
            loop {
                if rem[c] != u32::MAX {
                    base_len = rem[c];
                    base_end = end[c];
                    break;
                }
                stackp.push(c);
                match recv(c) {
                    None => {
                        base_end = 2;
                        break;
                    }
                    Some(r) => {
                        if based[r] {
                            base_end = 1;
                            break;
                        }
                        if !under[r] {
                            base_end = 3;
                            break;
                        }
                        c = r;
                    }
                }
            }
            for (i, &p) in stackp.iter().rev().enumerate() {
                rem[p] = base_len + i as u32 + 1;
                end[p] = base_end;
            }
        }
        let mut run_km: Vec<f32> = added.iter().map(|&k| rem[k] as f32 * CELL_KM).collect();
        let ends: Vec<String> = [(1u8, "at a based shore"), (2, "at the sea / no receiver"), (3, "on the trunk (leaves the law)")]
            .iter()
            .map(|&(e, nme)| format!("{nme} {:.1} %", 100.0 * added.iter().filter(|&&k| end[k] == e).count() as f32 / added.len().max(1) as f32))
            .collect();
        eprintln!(
            "   U · NEW under-law cells vs OFF: **{}** (+{:.1} points; {gone} left the set) · distance to a based shore p50 **{:.2} km** / p90 {:.2} km · within 2 km **{:.1} %** · law − terrain p50 **{:.1} m** / p90 {:.1} m · by population: {} · their run downstream p50 {:.2} km / p90 {:.2} km, ending {}",
            added.len(),
            100.0 * added.len() as f32 / ntrunk as f32,
            q(&mut dkm, 0.5),
            q(&mut dkm, 0.9),
            100.0 * near2 as f32 / added.len().max(1) as f32,
            q(&mut lt, 0.5),
            q(&mut lt, 0.9),
            by.join(" · "),
            q(&mut run_km, 0.5),
            q(&mut run_km, 0.9),
            ends.join(" · ")
        );
    }
    eprintln!("\n==========  end Finding 134-T / U . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 134-Δ — F133-F's residual (the témoin, extended lake base, A/B flat-pointer patches): its |Δz|
/// p50 / p90 / max, overall and per F133-F class, and the sample `carve_diag` laid each cell from in A and in B
/// (moved / same position other floor / none on one side). Declared in `f134_declared.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored f134_delta --nocapture
#[test]
#[ignore]
fn f134_delta() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, carve_diag, ocean_flood, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 134-Δ . the residual's laying sample (témoin, extended)  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let dc = C1DrainageConfig::default();
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
    let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
    let lake = d.lake_map.clone();
    drop(d);
    let spill = ocean_flood(&bf);
    let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
    let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
    let present: Vec<bool> = (0..n).map(|k| lake[k] != 0 || dep[k]).collect();
    let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
    let f = flow.filled.data.clone();
    let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
    let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
    let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
    let arbitrary = |seed: u64| -> Vec<u8> {
        (0..n)
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
                let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                cand[(hsh % cand.len() as u64) as usize]
            })
            .collect()
    };
    let (arb_a, arb_b) = (arbitrary(1), arbitrary(2));
    drop((flow, f, g, needs));
    let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_a[c];
            }
        }
    };
    let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_b[c];
            }
        }
    };
    let mm = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let vc = ValleyConstruction {
        lake_base: Some(LakeBase::InputLakesAndBasins),
        wall_sea_floor_m: Some(0.5),
        ..ValleyConstruction::new(F121_AGE_K, Some(0.1))
    };
    let ska = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pa));
    let skb = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pb));
    let (za, _, da) = carve_diag(&s1, &ska, &vc, &ss);
    let (zb, _, db) = carve_diag(&s1, &skb, &vc, &ss);
    let ch: Vec<usize> = (0..n).filter(|&k| !present[k] && za.data[k] != zb.data[k]).collect();
    let mut all: Vec<f32> = ch.iter().map(|&k| (mm(&za, k) - mm(&zb, k)).abs()).collect();
    eprintln!(
        "   the residual: **{}** cells outside the present lakes (F133-F: 12 509) · |Δz| p50 **{:.2}** / p90 **{:.2}** / max **{:.1} m**",
        ch.len(),
        q(&mut all, 0.5),
        q(&mut all, 0.9),
        q(&mut all, 1.0)
    );
    let flat = |sk: &Skeleton| -> Vec<(f32, f32, f32, f32)> { sk.polylines.iter().flatten().copied().collect() };
    let (sa, sb) = (flat(&ska), flat(&skb));
    eprintln!("   polylines A {} ({} samples) · B {} ({} samples)", ska.polylines.len(), sa.len(), skb.polylines.len(), sb.len());
    // F133-F's class (the origin walk) and the laying sample's class
    let differs = |k: usize| ska.chi_m[k].to_bits() != skb.chi_m[k].to_bits() || ska.base_alt_m[k].to_bits() != skb.base_alt_m[k].to_bits();
    let mut rows: HashMap<(String, &str), Vec<f32>> = HashMap::new();
    for &s in &ch {
        let mut c = s;
        let (mut seen, mut last) = (false, None);
        for _ in 0..n {
            if differs(c) {
                seen = true;
                last = Some(c);
            } else if seen {
                break;
            }
            let dd = ska.direction[c];
            if dd == DIR_NONE || s1.data[c] <= SEA {
                break;
            }
            c = nb(c, dd as usize);
        }
        let f133 = if last.is_none() { "F133-F: no differing cell on its path" } else { "F133-F: an origin on its path" };
        let (wa, wb) = (da.who[s], db.who[s]);
        let smp = if wa == u32::MAX || wb == u32::MAX {
            "no sample on one side".to_string()
        } else {
            let (pa_, pb_) = (sa[wa as usize], sb[wb as usize]);
            let pos = pa_.0.to_bits() == pb_.0.to_bits() && pa_.1.to_bits() == pb_.1.to_bits();
            let zf = pa_.2.to_bits() == pb_.2.to_bits();
            let hw = pa_.3.to_bits() == pb_.3.to_bits();
            match (pos, zf, hw) {
                (false, _, _) => "the sample MOVED (another position)".to_string(),
                (true, false, _) => "same position, another floor".to_string(),
                (true, true, false) => "same position and floor, another width".to_string(),
                (true, true, true) => "the same sample (position, floor, width)".to_string(),
            }
        };
        rows.entry((smp, f133)).or_default().push((mm(&za, s) - mm(&zb, s)).abs());
    }
    let mut keys: Vec<_> = rows.keys().cloned().collect();
    keys.sort();
    for k in keys {
        let v = rows.get_mut(&k).unwrap();
        let c = v.len();
        eprintln!(
            "   {:<44} · {:<38} **{c}** cells · |Δz| p50 {:.2} / p90 {:.2} / max {:.1} m",
            k.0,
            k.1,
            q(v, 0.5),
            q(v, 0.9),
            q(v, 1.0)
        );
    }
    eprintln!("\n==========  end Finding 134-Δ . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// θ's fit and bootstrap of `theta_laid_ci`, on an explicit list of (ln A, ln S) pairs (Finding 135-M: the list
/// is filtered by a lake mask before the fit). Same estimator, same resampling seed.
fn theta_fit_ci(v: &[(f64, f64)]) -> (f32, f32, f32, usize) {
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
    let m = v.len();
    if m < 3 {
        return (f32::NAN, f32::NAN, f32::NAN, m);
    }
    let t = fit(&mut (0..m));
    let mut rng = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    let mut bs: Vec<f64> = (0..1000)
        .map(|_| {
            let picks: Vec<usize> = (0..m).map(|_| (next() % m as u64) as usize).collect();
            fit(&mut picks.into_iter())
        })
        .collect();
    bs.sort_by(f64::total_cmp);
    (t as f32, bs[25] as f32, bs[974] as f32, m)
}

/// The OFF skeleton's trunk links ≥ 10 km² in `theta_laid_ci`'s order: (cell, receiver, ln A, step length m).
fn trunk_links(sk: &Skeleton) -> Vec<(usize, usize, f64, f32)> {
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let (w, h) = (sk.width, sk.height);
    let mut out = Vec::new();
    for k in 0..w * h {
        if !sk.trunk[k] || sk.area_km2[k] < 10.0 || sk.direction[k] == DIR_NONE {
            continue;
        }
        let d = sk.direction[k] as usize;
        let r = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize * w + ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
        if !sk.trunk[r] {
            continue;
        }
        let dist = sk.cell_m * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
        out.push((k, r, (sk.area_km2[k] as f64).ln(), dist));
    }
    out
}

/// θ on the links whose index passes `keep`, from the end altitudes `zk`, `zr` (m), as `theta_laid_ci` does.
fn theta_links(links: &[(usize, usize, f64, f32)], zk: &[f32], zr: &[f32], keep: &dyn Fn(usize) -> bool) -> (f32, f32, f32, usize) {
    let mut v = Vec::new();
    for (i, &(_, _, x, dist)) in links.iter().enumerate() {
        if !keep(i) {
            continue;
        }
        let s = (zk[i] - zr[i]) / dist;
        if s > 1e-4 {
            v.push((x, (s as f64).ln()));
        }
    }
    theta_fit_ci(&v)
}

/// ADR Finding 135-W1 — every configuration difference between the benches' témoin and the viz state
/// "C2 /10 col (défaut)": the upscale config (recursive JSON diff), the eroded digests, the drainage config.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f135_w1 --nocapture
#[test]
#[ignore]
fn f135_w1() {
    use common::{bench_eroded_digest, bench_upscale_cfg, viz_dcfg};
    use serde_json::Value;
    eprintln!("\n==========  Finding 135-W1 . the témoin against the viz state, key by key  ==========");
    fn diff(path: &str, a: &Value, b: &Value, out: &mut Vec<String>) {
        match (a, b) {
            (Value::Object(x), Value::Object(y)) => {
                let mut keys: Vec<&String> = x.keys().chain(y.keys()).collect();
                keys.sort();
                keys.dedup();
                for k in keys {
                    diff(&format!("{path}.{k}"), x.get(k).unwrap_or(&Value::Null), y.get(k).unwrap_or(&Value::Null), out);
                }
            }
            _ if a == b => {}
            _ => out.push(format!("{path}: témoin {a} · viz {b}")),
        }
    }
    let k = F121_AGE_K;
    let viz = ValleyConstruction::new(k, Some(0.1));
    let tem = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..viz };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let a = serde_json::to_value(bench_upscale_cfg(kn(tem))).expect("json");
    let b = serde_json::to_value(bench_upscale_cfg(kn(viz))).expect("json");
    let mut d = Vec::new();
    diff("upscale", &a, &b, &mut d);
    eprintln!("   UPSCALE CONFIG: **{} difference(s)**", d.len());
    for l in &d {
        eprintln!("      {l}");
    }
    let (dt, dv) = (bench_eroded_digest(kn(tem), PSEED), bench_eroded_digest(kn(viz), PSEED));
    eprintln!("   eroded digest: témoin {dt} · viz state {dv} (the guard's C2 /10 col: cf1d4539c7afe588)");
    assert_eq!(dv, "cf1d4539c7afe588", "the viz state's digest moved");
    let mut dd = Vec::new();
    diff("drainage", &serde_json::to_value(dcfg()).expect("json"), &serde_json::to_value(viz_dcfg()).expect("json"), &mut dd);
    eprintln!("   DRAINAGE CONFIG (the benches' dcfg() vs run_hd's): **{} difference(s)**", dd.len());
    for l in &dd {
        eprintln!("      {l}");
    }
    eprintln!("\n==========  end Finding 135-W1  ==========\n");
}

/// ADR Finding 135-W3 / M — the realigned témoin (the viz state + run_hd's tail), OFF / F132 / extended: F133's
/// table with the guard on every world, then the trunk metrics without the lakes (θ (i), (ii), (iii) with CI, θ on
/// the carved links, the under-law share). Declared in `f135_declared.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f135_w3m --nocapture
#[test]
#[ignore]
fn f135_w3m() {
    use common::{bench_eroded_digest, build_world, viz_dcfg, viz_hd_lakes_on};
    use ymir_core::tectonics_c1::bench_guard::{check_lakes, entries, field_hash, lake_fingerprint};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, FlowConfig, breach_monotone_protected, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let near = 2.0 / CELL_KM;
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 135-W3 / M . the realigned témoin (viz state + run_hd tail)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &vd, &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let links = trunk_links(&sk);
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    eprintln!("   the OFF skeleton: {} trunk links ≥ 10 km² · {} trunk cells", links.len(), trunk10.len());
    // present lakes on S1 and the flat-resolution patches (F133-F)
    let dc0 = C1DrainageConfig::default();
    let (lake, dep, target, arb_a, arb_b) = {
        let d = c1_drainage_windowed(&s1, None, &dc0, &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let lake = d.lake_map.clone();
        drop(d);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        drop(spill);
        let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc0.flat_perturbation.clone(), dinf: dc0.dinf });
        let f = flow.filled.data.clone();
        let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
        let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc0.flat_perturbation.as_ref(), w, h);
        let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
        let arbitrary = |seed: u64| -> Vec<u8> {
            (0..n)
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
                    let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                    cand[(hsh % cand.len() as u64) as usize]
                })
                .collect()
        };
        let (a, b) = (arbitrary(1), arbitrary(2));
        (lake, dep, target, a, b)
    };
    let present: Vec<bool> = (0..n).map(|k| lake[k] != 0 || dep[k]).collect();
    let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_a[c];
            }
        }
    };
    let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_b[c];
            }
        }
    };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let fields = entries();
    let at_links = |z: &[f32]| -> (Vec<f32>, Vec<f32>) { (links.iter().map(|l| z[l.0]).collect(), links.iter().map(|l| z[l.1]).collect()) };
    struct Var {
        label: &'static str,
        st: [(Vec<f32>, Vec<f32>); 3],
        zb: (Vec<f32>, Vec<f32>),
        carved_l: Vec<bool>,
        own_l: Vec<bool>,
        under_t: Vec<bool>,
        own_t: Vec<bool>,
    }
    let mut union = vec![false; n];
    let mut vars: Vec<Var> = Vec::new();
    let variants = [
        ("OFF", off),
        ("ON F132 (InputLakes)", ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..off }),
        ("ON EXTENDED (InputLakesAndBasins)", ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off }),
    ];
    for (label, vc) in variants {
        let t = Instant::now();
        eprintln!("\n────────── realigned témoin · {label} ──────────");
        let kw = Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
        let digest = bench_eroded_digest(kw, PSEED);
        let wd = build_world(kw, None, PSEED, None);
        let fh = format!("{:016x}", field_hash(&wd.heightmap));
        let fg = match fields.iter().find(|e| e.digest == digest) {
            Some(e) if e.field_hash == fh => format!("Match ({})", e.label),
            Some(e) => format!("**MISMATCH** (bench {})", e.field_hash),
            None => "NoReference".to_string(),
        };
        let v = viz_hd_lakes_on(&wd, kw, PSEED, 45.0, 40.0);
        let lg = check_lakes(&v.digest, &v.drainage.lake_map, &v.drainage.lakes);
        eprintln!(
            "   GUARD · field {digest} = {fh} → {fg} · lakes key {} · {} → {lg:?}",
            v.digest,
            lake_fingerprint(&v.drainage.lake_map, &v.drainage.lakes)
        );
        let g = &wd.heightmap;
        // the table
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(g, &vd, &ss, &dc_e, DOMAIN_KM);
        drop(cl_e);
        set_dump(true);
        let cr = f95_criteria(g, &v.conditioned, &pre, DELIVERED_P50_M, &ss, &vd, cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = Vec::new();
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons.push((floor % w, floor / w));
            }
        }
        drop(fill_del);
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let cls = classes(&pre, &skp, &mkp);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        drop((sea, dsea, cw, skp, mkp));
        let co = coast(g, &ss, &dwall, near);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let r8w = r8_terrain(g, &cls, 16);
        drop((co, dwall, cls));
        let segs = &v.drainage.rivers.segments;
        let wci: Vec<usize> = (0..segs.len()).filter(|&i| v.drainage.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2).collect();
        let pts: Vec<&[(u32, u32)]> = wci.iter().map(|&i| own(&segs[i])).collect();
        eprintln!(
            "   TABLE · canyons **{}** {:?} · spurs near a coastal wall **{rn:.4} /km** · R8 terrain **{r8w:.4}** · R8 network (chord 8) {:.4} · lakes **{}** ({:.1} km²) · relief p50 **{:.1} m**",
            canyons.len(),
            canyons,
            r8_chords(&pts, 8),
            v.drainage.lakes.len(),
            v.drainage.lakes.iter().map(|l| l.area_km2).sum::<f32>(),
            cr.p50
        );
        let lm = &v.drainage.lake_map;
        for k in 0..n {
            if lm[k] != 0 {
                union[k] = true;
            }
        }
        let own_l: Vec<bool> = links.iter().map(|l| lm[l.0] != 0 || lm[l.1] != 0).collect();
        let own_t: Vec<bool> = trunk10.iter().map(|&k| lm[k] != 0).collect();
        let st3 = at_links(&metres(g));
        drop(v);
        drop(wd);
        // (i) as built, (ii) its run_hd breach
        let w2 = build_world(Knobs { valley: Some(vc), no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, None, PSEED, None);
        let st1 = at_links(&metres(&w2.heightmap));
        let st2 = {
            let d2 = c1_drainage_windowed(&w2.heightmap, None, &vd, &ss, DOMAIN_KM);
            let prot = w2.volc.enabled.then(|| crater_protect_mask(&w2.craters, w, h));
            let b2 = breach_monotone_protected(&w2.heightmap, &d2.flow.filled, &d2.lake_map, 0.5, w, h, prot.as_deref());
            at_links(&metres(&b2))
        };
        drop(w2);
        // the construction on S1
        let skv = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (built, mk) = carve(&s1, &skv, &vc, &ss);
        let zb = at_links(&metres(&built));
        drop(built);
        let carved_l: Vec<bool> = links.iter().map(|l| mk.carved[l.0] && mk.carved[l.1]).collect();
        let under_t: Vec<bool> = trunk10
            .iter()
            .map(|&k| {
                let law = if skv.chi_m[k].is_finite() { skv.floor_m(k, vc.age_k) } else { f32::NAN };
                !mk.carved[k] && law.is_finite() && zs1[k] < law
            })
            .collect();
        drop((skv, mk));
        // Δz between flat resolutions
        let ska = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pa));
        let (za, _) = carve(&s1, &ska, &vc, &ss);
        drop(ska);
        let skb = skeleton_patched(&s1, &vc, &ss, DOMAIN_KM, Some(&pb));
        let (zbb, _) = carve(&s1, &skb, &vc, &ss);
        drop(skb);
        let mut dz: Vec<f32> = (0..n)
            .filter(|&k| !present[k] && za.data[k] != zbb.data[k])
            .map(|k| (c1_altitude_norm_to_metres(za.data[k], &ss) - c1_altitude_norm_to_metres(zbb.data[k], &ss)).abs())
            .collect();
        let out_input = (0..n).filter(|&k| lake[k] == 0 && za.data[k] != zbb.data[k]).count();
        drop((za, zbb));
        let nd = dz.len();
        eprintln!(
            "   TABLE · Δz between flat resolutions OUTSIDE the present lakes **{nd}** cells (outside the input lakes alone {out_input}) · |Δz| p50 {:.2} / p90 {:.2} / max {:.1} m ({:.0} s)",
            q(&mut dz, 0.5),
            q(&mut dz, 0.9),
            q(&mut dz, 1.0),
            t.elapsed().as_secs_f64()
        );
        vars.push(Var { label, st: [st1, st2, st3], zb, carved_l, own_l, under_t, own_t });
    }
    // M — the trunk metrics, with and without the lakes
    let ul: Vec<bool> = links.iter().map(|l| union[l.0] || union[l.1]).collect();
    let ut: Vec<bool> = trunk10.iter().map(|&k| union[k]).collect();
    eprintln!(
        "\n   M · the lake mask: the UNION of the three worlds' final lakes (run_hd tail) · {} cells · {} of {} links touch it · {} of {} trunk cells in it",
        union.iter().filter(|&&b| b).count(),
        ul.iter().filter(|&&b| b).count(),
        links.len(),
        ut.iter().filter(|&&b| b).count(),
        trunk10.len()
    );
    let fmt = |r: (f32, f32, f32, usize)| format!("**{:.3}** [{:.3}, {:.3}] ({})", r.0, r.1, r.2, r.3);
    let sn = ["(i) as built", "(ii) + run_hd breach", "(iii) the world"];
    for v in &vars {
        eprintln!("   ── {} ──", v.label);
        for (si, (zk, zr)) in v.st.iter().enumerate() {
            let all = theta_links(&links, zk, zr, &|_| true);
            let nol = theta_links(&links, zk, zr, &|i| !ul[i]);
            let own = if si == 2 { format!(" · without its OWN lakes {}", fmt(theta_links(&links, zk, zr, &|i| !v.own_l[i]))) } else { String::new() };
            eprintln!("      θ {:<22} all {} · WITHOUT the lakes (union) {}{own}", sn[si], fmt(all), fmt(nol));
        }
        let ca = theta_links(&links, &v.zb.0, &v.zb.1, &|i| v.carved_l[i]);
        let cn = theta_links(&links, &v.zb.0, &v.zb.1, &|i| v.carved_l[i] && !ul[i]);
        eprintln!("      θ on the carved links   all {} · WITHOUT the lakes (union) {}", fmt(ca), fmt(cn));
        let nall = v.under_t.iter().filter(|&&b| b).count();
        let tn = (0..trunk10.len()).filter(|&i| !ut[i]).count();
        let nno = (0..trunk10.len()).filter(|&i| v.under_t[i] && !ut[i]).count();
        let tno = (0..trunk10.len()).filter(|&i| !v.own_t[i]).count();
        let nown = (0..trunk10.len()).filter(|&i| v.under_t[i] && !v.own_t[i]).count();
        eprintln!(
            "      under-law share: all **{:.1} %** ({nall} / {}) · WITHOUT the lakes (union) **{:.1} %** ({nno} / {tn}) · without its own lakes {:.1} % ({nown} / {tno})",
            100.0 * nall as f32 / trunk10.len() as f32,
            trunk10.len(),
            100.0 * nno as f32 / tn.max(1) as f32,
            100.0 * nown as f32 / tno.max(1) as f32
        );
    }
    eprintln!("\n==========  end Finding 135-W3 / M . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// A long profile as an SVG line chart: series of (km, m) with a colour and a dash, a vertical col marker,
/// axes with ticks and labels. Finding 135-K2 (no plotting library in the toolchain).
#[allow(clippy::type_complexity)]
fn profile_svg(title: &str, series: &[(&str, &str, &str, Vec<(f32, f32)>)], col_km: Option<f32>, note: &str) -> String {
    let (wd, ht, ml, mr, mt, mb) = (1100.0f32, 560.0f32, 80.0f32, 30.0f32, 50.0f32, 90.0f32);
    let pts = series.iter().flat_map(|s| s.3.iter()).filter(|p| p.1.is_finite());
    let (mut x0, mut x1, mut y0, mut y1) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
    for p in pts {
        x0 = x0.min(p.0);
        x1 = x1.max(p.0);
        y0 = y0.min(p.1);
        y1 = y1.max(p.1);
    }
    let ystep = if y1 - y0 > 2000.0 { 500.0 } else if y1 - y0 > 800.0 { 200.0 } else { 100.0 };
    y0 = (y0 / ystep).floor() * ystep;
    y1 = (y1 / ystep).ceil() * ystep;
    let px = |x: f32| ml + (x - x0) / (x1 - x0).max(1e-3) * (wd - ml - mr);
    let py = |y: f32| mt + (y1 - y) / (y1 - y0).max(1e-3) * (ht - mt - mb);
    let mut s = format!("<svg xmlns='http://www.w3.org/2000/svg' width='{wd}' height='{ht}' font-family='sans-serif' font-size='12'>\n<rect width='100%' height='100%' fill='white'/>\n<text x='{ml}' y='24' font-size='15' font-weight='bold'>{title}</text>\n");
    let xstep = if x1 - x0 > 60.0 { 10.0 } else { 5.0 };
    let mut x = (x0 / xstep).ceil() * xstep;
    while x <= x1 {
        s += &format!("<line x1='{0}' y1='{1}' x2='{0}' y2='{2}' stroke='#ddd'/><text x='{0}' y='{3}' text-anchor='middle'>{x:.0}</text>\n", px(x), mt, ht - mb, ht - mb + 16.0);
        x += xstep;
    }
    let mut y = y0;
    while y <= y1 + 0.1 {
        s += &format!("<line x1='{0}' y1='{1}' x2='{2}' y2='{1}' stroke='#ddd'/><text x='{3}' y='{4}' text-anchor='end'>{y:.0}</text>\n", ml, py(y), wd - mr, ml - 6.0, py(y) + 4.0);
        y += ystep;
    }
    s += &format!("<text x='{}' y='{}' text-anchor='middle'>distance along the path (km, upstream → downstream)</text>\n", (ml + wd - mr) / 2.0, ht - mb + 36.0);
    s += &format!("<text x='18' y='{}' transform='rotate(-90 18 {})' text-anchor='middle'>altitude (m)</text>\n", (mt + ht - mb) / 2.0, (mt + ht - mb) / 2.0);
    if let Some(c) = col_km {
        s += &format!("<line x1='{0}' y1='{1}' x2='{0}' y2='{2}' stroke='black' stroke-width='1.5' stroke-dasharray='2,3'/><text x='{3}' y='{4}'>col</text>\n", px(c), mt, ht - mb, px(c) + 4.0, mt + 12.0);
    }
    for (li, (name, colr, dash, p)) in series.iter().enumerate() {
        let mut d = String::new();
        let mut pen = false;
        for &(xx, yy) in p {
            if !yy.is_finite() {
                pen = false;
                continue;
            }
            d += &format!("{}{:.1},{:.1} ", if pen { "L" } else { "M" }, px(xx), py(yy));
            pen = true;
        }
        s += &format!("<path d='{d}' fill='none' stroke='{colr}' stroke-width='1.6' stroke-dasharray='{dash}'/>\n");
        let ly = ht - 30.0 + 0.0 * li as f32;
        let lx = ml + li as f32 * 240.0;
        s += &format!("<line x1='{lx}' y1='{ly}' x2='{}' y2='{ly}' stroke='{colr}' stroke-width='2' stroke-dasharray='{dash}'/><text x='{}' y='{}'>{name}</text>\n", lx + 28.0, lx + 34.0, ly + 4.0);
    }
    s += &format!("<text x='{ml}' y='{}' font-size='11' fill='#555'>{note}</text>\n", ht - 8.0);
    s + "</svg>\n"
}

/// ADR Finding 135-K — the cuts through the cols (realigned témoin): every input lake and closed depression, its
/// col and the built floor there in OFF and extended (K1); the long profiles of lakes 1, 2, 11 (K2, SVG); lake 2's
/// Δz by distance upstream of its col (K3); which extended cuts remain (K4). Declared in `f135_declared.md`.
/// SVGs to `F135_DIR` if set.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f135_k --nocapture
#[test]
#[ignore]
fn f135_k() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 135-K . the cuts through the cols (realigned témoin)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nb = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let (lake, dep, zbf) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep, metres(&bf))
    };
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let (b_off, mk_off) = carve(&s1, &sk, &off, &ss);
    let zoff_b = metres(&b_off);
    drop(b_off);
    let (zext_b, mk_ext) = {
        let skx = skeleton(&s1, &ext, &ss, DOMAIN_KM);
        let (b, m) = carve(&s1, &skx, &ext, &ss);
        (metres(&b), m)
    };
    // the items: input lakes and closed-depression components
    let mut items: Vec<(String, Vec<usize>, bool)> = Vec::new(); // (name, cells, on the breached field)
    {
        let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
        for k in 0..n {
            if lake[k] != 0 {
                by.entry(lake[k]).or_default().push(k);
            }
        }
        let mut ids: Vec<u32> = by.keys().copied().collect();
        ids.sort_unstable();
        for id in ids {
            items.push((format!("input lake {id}"), by.remove(&id).unwrap(), false));
        }
        let mut seen = vec![false; n];
        let mut cid = 0;
        for s in 0..n {
            if !dep[s] || seen[s] {
                continue;
            }
            cid += 1;
            let mut comp = vec![s];
            seen[s] = true;
            let mut i = 0;
            while i < comp.len() {
                let c = comp[i];
                i += 1;
                for k in 0..8 {
                    let m = nb(c, k);
                    if dep[m] && !seen[m] {
                        seen[m] = true;
                        comp.push(m);
                    }
                }
            }
            comp.sort_unstable();
            items.push((format!("closed depression {cid}"), comp, true));
        }
    }
    // the outer ring, its lowest cell (the col) and the cuts
    let ring_of = |cells: &[usize]| -> Vec<usize> {
        let mut inside = HashMap::with_capacity(cells.len());
        for &c in cells {
            inside.insert(c, ());
        }
        let mut r: Vec<usize> = Vec::new();
        for &c in cells {
            for k in 0..8 {
                let m = nb(c, k);
                if !inside.contains_key(&m) {
                    r.push(m);
                }
            }
        }
        r.sort_unstable();
        r.dedup();
        r
    };
    struct Cut {
        name: String,
        cells: Vec<usize>,
        col: usize,
        level: f32,
        cut_off: f32,
        cut_ext: f32,
        colcut_off: f32,
        colcut_ext: f32,
    }
    let mut cuts: Vec<Cut> = Vec::new();
    for (name, cells, on_bf) in items {
        let ring = ring_of(&cells);
        let zr = |c: usize| if on_bf { zbf[c] } else { zs1[c] };
        let col = *ring.iter().min_by(|&&a, &&b| zr(a).total_cmp(&zr(b))).expect("ring");
        let level = zr(col);
        let mo = ring.iter().map(|&c| zoff_b[c]).fold(f32::INFINITY, f32::min);
        let me = ring.iter().map(|&c| zext_b[c]).fold(f32::INFINITY, f32::min);
        cuts.push(Cut { name, cells, col, level, cut_off: level - mo, cut_ext: level - me, colcut_off: level - zoff_b[col], colcut_ext: level - zext_b[col] });
    }
    let count = |f: &dyn Fn(&Cut) -> f32, t: f32| cuts.iter().filter(|c| f(c) > t).count();
    eprintln!("   K1 · {} items ({} input lakes, {} closed depressions)", cuts.len(), cuts.iter().filter(|c| c.name.starts_with("input")).count(), cuts.iter().filter(|c| c.name.starts_with("closed")).count());
    for (st, f) in [("OFF", &(|c: &Cut| c.cut_off) as &dyn Fn(&Cut) -> f32), ("EXTENDED", &(|c: &Cut| c.cut_ext) as &dyn Fn(&Cut) -> f32)] {
        eprintln!("   K1 · {st}: rims cut (level − lowest built ring cell) > 10 m **{}** · > 100 m **{}** · > 500 m **{}**", count(f, 10.0), count(f, 100.0), count(f, 500.0));
    }
    for (st, f) in [("OFF", &(|c: &Cut| c.colcut_off) as &dyn Fn(&Cut) -> f32), ("EXTENDED", &(|c: &Cut| c.colcut_ext) as &dyn Fn(&Cut) -> f32)] {
        eprintln!("   K1 · {st}: at the col cell itself > 10 m {} · > 100 m {} · > 500 m {}", count(f, 10.0), count(f, 100.0), count(f, 500.0));
    }
    let mut order: Vec<usize> = (0..cuts.len()).collect();
    order.sort_by(|&a, &b| cuts[b].cut_off.total_cmp(&cuts[a].cut_off));
    eprintln!("   K1 · every item with a cut > 10 m in either state (sorted by the OFF cut):");
    for &i in &order {
        let c = &cuts[i];
        if c.cut_off <= 10.0 && c.cut_ext <= 10.0 {
            continue;
        }
        eprintln!(
            "      {:<24} {:>8} cells ({:>7.2} km²) · col ({}, {}) at **{:.1} m** · built floor at the col OFF {:.1} / EXT {:.1} m · rim cut OFF **{:.1} m** / EXT **{:.1} m** · at the col OFF {:.1} / EXT {:.1} m",
            c.name,
            c.cells.len(),
            c.cells.len() as f32 * CELL_KM * CELL_KM,
            c.col % w,
            c.col / w,
            c.level,
            zoff_b[c.col],
            zext_b[c.col],
            c.cut_off,
            c.cut_ext,
            c.colcut_off,
            c.colcut_ext
        );
    }
    // the worlds (conditioned fields, final lakes, the OFF drainage)
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let (co, dir_off) = {
        let wd = build_world(kn(off), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(off), PSEED, 45.0, 40.0);
        (metres(&v.conditioned), v.drainage.flow.direction.clone())
    };
    let (cn, lm_ext, lakes_ext) = {
        let wd = build_world(kn(ext), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(ext), PSEED, 45.0, 40.0);
        (metres(&v.conditioned), v.drainage.lake_map.clone(), v.drainage.lakes.clone())
    };
    // K4 — the extended cuts that remain
    let present = |k: usize| lake[k] != 0 || dep[k];
    let mut unbased: Vec<(u32, Vec<usize>)> = Vec::new();
    for l in &lakes_ext {
        let cells: Vec<usize> = (0..n).filter(|&k| lm_ext[k] == l.base.id).collect();
        let cov = cells.iter().filter(|&&k| present(k)).count();
        if (cov as f32) < 0.5 * cells.len().max(1) as f32 {
            unbased.push((l.base.id, cells));
        }
    }
    let absent: Vec<bool> = cuts.iter().map(|c| (c.cells.iter().filter(|&&k| lm_ext[k] != 0).count() as f32) < 0.5 * c.cells.len() as f32).collect();
    eprintln!("   K4 · bases on a lake ABSENT at the end (items < 50 % under the extended final lakes): **{}** · unbased final lakes: **{}** {:?}", absent.iter().filter(|&&b| b).count(), unbased.len(), unbased.iter().map(|u| u.0).collect::<Vec<_>>());
    for (i, c) in cuts.iter().enumerate() {
        if c.cut_ext <= 10.0 {
            continue;
        }
        let ub = unbased.iter().find(|u| u.1.iter().any(|k| c.cells.binary_search(k).is_ok())).map(|u| u.0);
        let lab = if absent[i] {
            "a base on an ABSENT lake".to_string()
        } else if let Some(id) = ub {
            format!("overlaps the unbased final lake {id}")
        } else {
            "OTHER (a present, based lake)".to_string()
        };
        eprintln!("      EXT cut {:<24} rim cut **{:.1} m** · col ({}, {}) · {lab}", c.name, c.cut_ext, c.col % w, c.col / w);
    }
    for (id, cells) in &unbased {
        let mut sorted = cells.clone();
        sorted.sort_unstable();
        let ring = ring_of(&sorted);
        let col = *ring.iter().min_by(|&&a, &&b| zs1[a].total_cmp(&zs1[b])).expect("ring");
        let me = ring.iter().map(|&c| zext_b[c]).fold(f32::INFINITY, f32::min);
        let mo = ring.iter().map(|&c| zoff_b[c]).fold(f32::INFINITY, f32::min);
        eprintln!("      unbased final lake {id}: {} cells · S1 col ({}, {}) at {:.1} m · rim cut OFF {:.1} / EXT **{:.1} m**", cells.len(), col % w, col / w, zs1[col], zs1[col] - mo, zs1[col] - me);
    }
    // K2 — the long profiles, K3 — lake 2's Δz by distance upstream of its col
    let dir_out = std::env::var("F135_DIR").ok();
    for (lid, centre) in [(2u32, (4019.0f64, 2839.0f64)), (11, (1429.0, 4287.0)), (1, (3491.0, 2500.0))] {
        let cells: Vec<usize> = (0..n).filter(|&k| lm_ext[k] == lid).collect();
        let (cx, cy) = (cells.iter().map(|&k| (k % w) as f64).sum::<f64>() / cells.len() as f64, cells.iter().map(|&k| (k / w) as f64).sum::<f64>() / cells.len() as f64);
        let l = lakes_ext.iter().find(|l| l.base.id == lid).expect("lake");
        let mut sorted = cells.clone();
        sorted.sort_unstable();
        let ring = ring_of(&sorted);
        let col = *ring.iter().min_by(|&&a, &&b| zs1[a].total_cmp(&zs1[b])).expect("ring");
        // the exit: the largest-area skeleton cell of the footprint ∪ ring
        let exit = *sorted.iter().chain(ring.iter()).max_by(|&&a, &&b| sk.area_km2[a].total_cmp(&sk.area_km2[b])).expect("exit");
        let step = |d: usize| -> f32 { CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 } };
        let mut down = vec![exit];
        let mut c = exit;
        let mut len = 0.0;
        while len < 60.0 {
            let d = sk.direction[c];
            if d == DIR_NONE || s1.data[c] <= SEA {
                break;
            }
            len += step(d as usize);
            c = nb(c, d as usize);
            down.push(c);
        }
        let mut up = Vec::new();
        let mut c = exit;
        let mut len = 0.0;
        while len < 40.0 {
            let mut best: Option<usize> = None;
            for k in 0..8 {
                let m = nb(c, k);
                let dm = sk.direction[m];
                if dm != DIR_NONE && nb(m, dm as usize) == c && best.is_none_or(|b| sk.area_km2[m] > sk.area_km2[b]) {
                    best = Some(m);
                }
            }
            match best {
                Some(m) if sk.area_km2[m] >= 1.0 => {
                    len += CELL_KM * if (m % w != c % w) && (m / w != c / w) { std::f32::consts::SQRT_2 } else { 1.0 };
                    up.push(m);
                    c = m;
                }
                _ => break,
            }
        }
        up.reverse();
        let path: Vec<usize> = up.into_iter().chain(down).collect();
        let mut xs = vec![0f32; path.len()];
        for i in 1..path.len() {
            let (a, b) = (path[i - 1], path[i]);
            xs[i] = xs[i - 1] + CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
        }
        let col_on_path = path.iter().position(|&p| p == col);
        let near_col = path.iter().enumerate().min_by(|a, b| {
            let d = |p: usize| ((p % w) as f32 - (col % w) as f32).hypot((p / w) as f32 - (col / w) as f32);
            d(*a.1).total_cmp(&d(*b.1))
        });
        let (ci, cdist) = near_col.map(|(i, &p)| (i, ((p % w) as f32 - (col % w) as f32).hypot((p / w) as f32 - (col / w) as f32) * CELL_KM)).unwrap();
        let ser = |z: &[f32]| -> Vec<(f32, f32)> { path.iter().enumerate().map(|(i, &p)| (xs[i], z[p])).collect() };
        let lake_line: Vec<(f32, f32)> = path.iter().enumerate().map(|(i, &p)| (xs[i], if lm_ext[p] == lid { l.level_m } else { f32::NAN })).collect();
        let max_cut = path.iter().filter(|&&p| (cn[p] - co[p]) > 0.0).map(|&p| cn[p] - co[p]).fold(0f32, f32::max);
        eprintln!(
            "\n   K2 · lake {lid} · centre ({cx:.0}, {cy:.0}) (F134 ({:.0}, {:.0})) · level {:.1} m · depth {:.1} m · S1 col ({}, {}) at **{:.1} m** · the path: {} cells, {:.1} km, exit ({}, {}) · col {} the path (nearest path point {cdist:.2} km away, at {:.1} km) · OFF conditioned at that point **{:.1} m** vs EXT **{:.1} m** · max (EXT − OFF) along the path {max_cut:.1} m",
            centre.0, centre.1, l.level_m, l.depth_m, col % w, col / w, zs1[col], path.len(), xs[xs.len() - 1], exit % w, exit / w,
            if col_on_path.is_some() { "ON" } else { "off" }, xs[ci], co[path[ci]], cn[path[ci]]
        );
        if let Some(d) = &dir_out {
            let svg = profile_svg(
                &format!("Lake {lid} — long profile, OFF vs ON extended (realigned témoin, viz state C2 /10 col)"),
                &[
                    ("S1 (construction input)", "#999999", "5,4", ser(&zs1)),
                    ("OFF (conditioned)", "#d62728", "", ser(&co)),
                    ("ON extended (conditioned)", "#1f77b4", "", ser(&cn)),
                    ("ON lake level", "#17becf", "1,3", lake_line),
                ],
                Some(xs[ci]),
                &format!("Path on the OFF skeleton: max-area donors upstream, D8 receivers downstream. Col = lowest S1 cell of the lake's outer ring, at {:.1} m (marked at the nearest path point, {cdist:.2} km off). Lake {lid}: level {:.1} m, depth {:.1} m.", zs1[col], l.level_m, l.depth_m),
            );
            let p = std::path::Path::new(d).join(format!("k2_profile_lake{lid}.svg"));
            std::fs::write(&p, svg).expect("svg");
            eprintln!("      profile: {}", p.display());
        }
        if lid == 2 {
            // K3 — distance upstream of the col along the OFF drainage (D8 length to within 3 cells of the col)
            let target: Vec<bool> = {
                let mut t = vec![false; n];
                let (x0, y0) = ((col % w) as i64, (col / w) as i64);
                for dy in -3i64..=3 {
                    for dx in -3i64..=3 {
                        let (x, y) = (x0 + dx, y0 + dy);
                        if x >= 0 && y >= 0 && x < w as i64 && y < h as i64 {
                            t[y as usize * w + x as usize] = true;
                        }
                    }
                }
                t
            };
            let mut dist = vec![f32::NAN; n];
            let mut done = vec![false; n];
            let mut path = Vec::new();
            for s in 0..n {
                if done[s] {
                    continue;
                }
                path.clear();
                let mut c = s;
                let base;
                loop {
                    if done[c] {
                        base = dist[c];
                        break;
                    }
                    if target[c] {
                        done[c] = true;
                        dist[c] = 0.0;
                        base = 0.0;
                        break;
                    }
                    path.push(c);
                    let d = dir_off[c];
                    if d == DIR_NONE || path.len() > 20000 {
                        base = f32::NAN;
                        break;
                    }
                    let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                    if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
                        base = f32::NAN;
                        break;
                    }
                    c = y as usize * w + x as usize;
                }
                let mut acc = base;
                for &p in path.iter().rev() {
                    let d = dir_off[p];
                    if acc.is_finite() && d != DIR_NONE {
                        acc += CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                    } else {
                        acc = f32::NAN;
                    }
                    dist[p] = acc;
                    done[p] = true;
                }
            }
            eprintln!("   K3 · lake 2 · Δz = EXT − OFF (conditioned) by OFF-drainage distance upstream of the col (cells draining through the col in OFF; footprint excluded):");
            let mut first_bins = Vec::new();
            for b in 0..15 {
                let (lo, hi) = (2.0 * b as f32, 2.0 * (b + 1) as f32);
                let sel: Vec<usize> = (0..n).filter(|&k| dist[k] >= lo && dist[k] < hi && lm_ext[k] != lid).collect();
                let mut dz: Vec<f32> = sel.iter().map(|&k| cn[k] - co[k]).collect();
                let mut dzc: Vec<f32> = sel.iter().filter(|&&k| mk_off.carved[k] && mk_ext.carved[k]).map(|&k| cn[k] - co[k]).collect();
                let nc = dzc.len();
                let (p50, p50c) = (q_(&mut dz, 0.5), q_(&mut dzc, 0.5));
                first_bins.push(p50c);
                eprintln!(
                    "      {lo:>4.0}–{hi:<3.0} km · {:>7} cells · p50 Δz **{p50:>7.1} m** (p90 {:>7.1}) · carved in both {nc:>6} · their p50 Δz **{p50c:>7.1} m**",
                    sel.len(),
                    q_(&mut dz, 0.9)
                );
            }
            let _ = first_bins;
        }
    }
    eprintln!("\n==========  end Finding 135-K . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

fn q_(v: &mut Vec<f32>, p: f32) -> f32 {
    if v.is_empty() {
        return f32::NAN;
    }
    v.sort_by(f32::total_cmp);
    v[((v.len() - 1) as f32 * p) as usize]
}

/// ADR Finding 135-C — what `carve_diag` lays each constructed cell from: its `who` sample's cone (the
/// nearest-sample propagation) or another line's cone (the cross-line minimum). The cone is rebuilt in the bench
/// and must reproduce `out` bit for bit. Populations: the realigned témoin's constructed cells, the residual of path
/// dependence (extended, A/B), Finding 127's four dams (B2 → A_c realigned). Declared in `f135_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f135_c --nocapture
#[test]
#[ignore]
fn f135_c() {
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    use ymir_core::tectonics_c1::valley_construction::{CarveDiag, LakeBase, carve_diag, ocean_flood, skeleton_patched};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, FlowConfig, compute_flow, flat_resolution};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 135-C . carve_diag: the nearest sample vs the cross-line minimum  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zfield: Vec<f32> = s1.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    // 0 = not carved, 1 = laid by its `who` sample (nearest), 2 = by another line's cone (cross-line minimum)
    // returns (classes, the laying sample, mismatches of the rebuild)
    let classify = |sk: &Skeleton, vc: &ValleyConstruction, out: &GridF32, carved: &[bool], dg: &CarveDiag| -> (Vec<u8>, Vec<u32>, usize) {
        assert!(vc.wall_profile.is_none() && vc.wall_sea_floor_m.is_none() && !vc.trunk_band, "the rebuild covers the plain cone only");
        let src: Vec<(f32, f32, f32, f32)> = sk.polylines.iter().flatten().copied().collect();
        let tan = vc.wall_deg.to_radians().tan();
        let geo = |s: usize, c: usize| -> f32 {
            let (sx, sy, zf, hw) = src[s];
            let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
            let mut dx = (cx - sx).rem_euclid(w as f32);
            if dx > w as f32 / 2.0 {
                dx -= w as f32;
            }
            let mut dy = (cy - sy).rem_euclid(h as f32);
            if dy > h as f32 / 2.0 {
                dy -= h as f32;
            }
            let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
            let u = (d - hw).max(0.0);
            zf + tan * (u - 0.5 * 0.0)
        };
        let mut cls = vec![0u8; n];
        let mut by = vec![u32::MAX; n];
        let mut bad = 0usize;
        for c in 0..n {
            if !carved[c] {
                continue;
            }
            let me = dg.who[c] as usize;
            let mut v = geo(me, c);
            let mut from = me;
            let (x, y) = ((c % w) as i32, (c / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nbc = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                    let o = dg.who[nbc];
                    if o != u32::MAX && dg.line_of[o as usize] != dg.line_of[me] {
                        let vv = geo(o as usize, c);
                        if vv < v {
                            v = vv;
                            from = o as usize;
                        }
                    }
                }
            }
            if c1_metres_to_altitude_norm(v, &ss).to_bits() != out.data[c].to_bits() || v >= zfield[c] {
                bad += 1;
            }
            cls[c] = if from == me { 1 } else { 2 };
            by[c] = from as u32;
        }
        (cls, by, bad)
    };
    // 1 · the realigned témoin's constructed cells
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let (out, mk, dg) = carve_diag(&s1, &sk, &off, &ss);
    let (cls, _, bad) = classify(&sk, &off, &out, &mk.carved, &dg);
    let nc = mk.carved.iter().filter(|&&b| b).count();
    let n1 = cls.iter().filter(|&&c| c == 1).count();
    let n2 = cls.iter().filter(|&&c| c == 2).count();
    let fl1 = (0..n).filter(|&c| cls[c] == 1 && mk.floor[c]).count();
    let fl2 = (0..n).filter(|&c| cls[c] == 2 && mk.floor[c]).count();
    let t1 = (0..n).filter(|&c| cls[c] == 1 && sk.trunk[c]).count();
    let t2 = (0..n).filter(|&c| cls[c] == 2 && sk.trunk[c]).count();
    eprintln!(
        "   1 · realigned témoin (OFF): constructed cells **{nc}** · the rebuild reproduces `out` bit for bit on all but **{bad}** · laid by the NEAREST sample **{n1} ({:.2} %)** · by the CROSS-LINE minimum **{n2} ({:.2} %)** · on the floor: nearest {fl1} / cross {fl2} · on a skeleton trunk cell: nearest {t1} / cross {t2}",
        100.0 * n1 as f64 / nc as f64,
        100.0 * n2 as f64 / nc as f64
    );
    drop((out, mk, dg, cls, sk));
    // 2 · the residual of path dependence (extended, F133-F's A/B patches)
    let dc = C1DrainageConfig::default();
    let nbf = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let (present, target, arb_a, arb_b) = {
        let d = c1_drainage_windowed(&s1, None, &dc, &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let lake = d.lake_map.clone();
        drop(d);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let present: Vec<bool> = (0..n).map(|k| lake[k] != 0 || (bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps)).collect();
        let flow = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dc.flat_perturbation.clone(), dinf: dc.dinf });
        let f = flow.filled.data.clone();
        let is_ocean: Vec<bool> = f.iter().map(|&v| v <= C1_SEA_LEVEL_NORM).collect();
        let (needs, g) = flat_resolution(&flow.filled, &is_ocean, dc.flat_perturbation.as_ref(), w, h);
        let target: Vec<bool> = (0..n).map(|c| needs[c] && lake[c] != 0).collect();
        let arbitrary = |seed: u64| -> Vec<u8> {
            (0..n)
                .map(|c| {
                    if !target[c] {
                        return flow.direction[c];
                    }
                    let cand: Vec<u8> = (0..8u8)
                        .filter(|&k| {
                            let m = nbf(c, k as usize);
                            is_ocean[m] || f[m] < f[c] || (f[m] == f[c] && needs[m] && g[m] < g[c])
                        })
                        .collect();
                    if cand.is_empty() {
                        return flow.direction[c];
                    }
                    let hsh = ((c as u64) ^ seed.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
                    cand[(hsh % cand.len() as u64) as usize]
                })
                .collect()
        };
        let (a, b) = (arbitrary(1), arbitrary(2));
        (present, target, a, b)
    };
    let pa = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_a[c];
            }
        }
    };
    let pb = |_: &GridF32, _: &[u32], dd: &mut Vec<u8>| {
        for c in 0..n {
            if target[c] {
                dd[c] = arb_b[c];
            }
        }
    };
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let ska = skeleton_patched(&s1, &ext, &ss, DOMAIN_KM, Some(&pa));
    let (za, ma, da) = carve_diag(&s1, &ska, &ext, &ss);
    let (ca, _, bad_a) = classify(&ska, &ext, &za, &ma.carved, &da);
    drop((ska, da, ma));
    let skb = skeleton_patched(&s1, &ext, &ss, DOMAIN_KM, Some(&pb));
    let (zb, mb, db) = carve_diag(&s1, &skb, &ext, &ss);
    let (cb, _, bad_b) = classify(&skb, &ext, &zb, &mb.carved, &db);
    drop((skb, db, mb));
    let ch: Vec<usize> = (0..n).filter(|&k| !present[k] && za.data[k] != zb.data[k]).collect();
    let mut tab: HashMap<(u8, u8), usize> = HashMap::new();
    for &c in &ch {
        *tab.entry((ca[c], cb[c])).or_insert(0) += 1;
    }
    let name = |c: u8| ["not carved", "nearest", "cross-line"][c as usize];
    eprintln!("   2 · the residual of path dependence (realigned, extended): **{}** cells outside the present lakes · rebuild mismatches A {bad_a} / B {bad_b}", ch.len());
    let mut rows: Vec<_> = tab.into_iter().collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1));
    for ((a, b), c) in rows {
        eprintln!("      A {:<11} · B {:<11} **{c}** ({:.1} %)", name(a), name(b), 100.0 * c as f64 / ch.len() as f64);
    }
    drop((za, zb, ca, cb));
    // 3 · Finding 127's dams (B2 → A_c, realigned)
    let ac = ValleyConstruction { a_min_km2: 0.1, ..off };
    let skc = skeleton(&s1, &ac, &ss, DOMAIN_KM);
    let (oc, mc, dgc) = carve_diag(&s1, &skc, &ac, &ss);
    let (cc, byc, badc) = classify(&skc, &ac, &oc, &mc.carved, &dgc);
    let src: Vec<(f32, f32, f32, f32)> = skc.polylines.iter().flatten().copied().collect();
    eprintln!("   3 · Finding 127's dams, B2 → A_c realigned (rebuild mismatches {badc}): each true col ± 1 cell");
    for (body, (x, y)) in [(3, (3852usize, 2337usize)), (7, (4551, 3280)), (13, (2253, 4495)), (14, (1805, 4580))] {
        let mut row = Vec::new();
        for dy in -1i64..=1 {
            for dx in -1i64..=1 {
                let c = (y as i64 + dy) as usize * w + (x as i64 + dx) as usize;
                if cc[c] == 0 {
                    row.push("—".to_string());
                    continue;
                }
                let s = byc[c] as usize;
                let (sx, sy, zf, _) = src[s];
                let sc = (sy.floor() as i64).rem_euclid(h as i64) as usize * w + (sx.floor() as i64).rem_euclid(w as i64) as usize;
                row.push(format!("{} (line {} at {:.2} km², floor {:.0} m)", name(cc[c]), dgc.line_of[s], skc.area_km2[sc], zf));
            }
        }
        let c0 = y * w + x;
        eprintln!(
            "      body {body:>2} col ({x}, {y}) · S1 {:.1} m · built {:.1} m · skeleton area at the col {:.2} km² · the 3 × 3: {}",
            zfield[c0],
            c1_altitude_norm_to_metres(oc.data[c0], &ss),
            skc.area_km2[c0],
            row.join(" | ")
        );
    }
    let _ = mc;
    eprintln!("\n==========  end Finding 135-C . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 135-W2 — does `wall_sea_floor_m = 0.5` change the témoin's FIELD? The F131–F134 témoin and the viz
/// state, built and hashed; the cells that differ, and where (coast distance). Also the construction alone
/// (`carve(S1)`), where the clause acts.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f135_w2 --nocapture
#[test]
#[ignore]
fn f135_w2() {
    use ymir_core::tectonics_c1::bench_guard::field_hash;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 135-W2 . the témoin's wall_sea_floor_m, field by field  ==========");
    let viz = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let tem = ValleyConstruction { wall_sea_floor_m: Some(0.5), ..viz };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let a = build_field_seed(kn(tem), PSEED);
    let b = build_field_seed(kn(viz), PSEED);
    let n = a.data.len();
    let d: Vec<usize> = (0..n).filter(|&k| a.data[k].to_bits() != b.data[k].to_bits()).collect();
    let mut dz: Vec<f32> = d.iter().map(|&k| (c1_altitude_norm_to_metres(a.data[k], &ss) - c1_altitude_norm_to_metres(b.data[k], &ss)).abs()).collect();
    dz.sort_by(f32::total_cmp);
    eprintln!(
        "   the world: témoin {:016x} · viz state {:016x} · cells differing **{}** · |Δz| max {:.3} m",
        field_hash(&a),
        field_hash(&b),
        d.len(),
        dz.last().copied().unwrap_or(0.0)
    );
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let sk = skeleton(&s1, &viz, &ss, DOMAIN_KM);
    let (ca, _) = carve(&s1, &sk, &tem, &ss);
    let (cb, _) = carve(&s1, &sk, &viz, &ss);
    let dc = (0..n).filter(|&k| ca.data[k].to_bits() != cb.data[k].to_bits()).count();
    let below = (0..n).filter(|&k| cb.data[k] <= SEA && s1.data[k] > SEA).count();
    eprintln!("   the construction alone, carve(S1): cells differing **{dc}** · land cells of S1 the viz state's construction lays at or below the sea {below}");
    eprintln!("\n==========  end Finding 135-W2 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 135-K3 (redone) / M2 — K3's first instrument targeted the S1 col of the final lake 2, which OFF's
/// drainage does not pass (27 cells reached): the target is now the final lake 2's FOOTPRINT (distance upstream of the
/// lake along OFF's drainage). M2: where θ's fall survives the lakes' removal — on the construction `carve(S1)`, the
/// carved links and the under-law share without (a) the union of the final lakes, (b) the BASED set (input lakes ∪
/// closed depressions), (c) both. Added after the first results; declared in the report as such.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f135_km2 --nocapture
#[test]
#[ignore]
fn f135_km2() {
    use common::{build_world, viz_hd_lakes_on};
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 135-K3 (redone) / M2  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let f132 = ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..off };
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let mut union: Vec<bool> = Vec::new();
    let mut co = Vec::new();
    let mut cn = Vec::new();
    let mut dir_off = Vec::new();
    let mut lm_ext = Vec::new();
    for (i, vc) in [off, f132, ext].into_iter().enumerate() {
        let wd = build_world(kn(vc), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0);
        if union.is_empty() {
            union = vec![false; v.drainage.lake_map.len()];
        }
        for (k, &l) in v.drainage.lake_map.iter().enumerate() {
            if l != 0 {
                union[k] = true;
            }
        }
        match i {
            0 => {
                co = metres(&v.conditioned);
                dir_off = v.drainage.flow.direction.clone();
            }
            2 => {
                cn = metres(&v.conditioned);
                lm_ext = v.drainage.lake_map.clone();
            }
            _ => {}
        }
    }
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let (lake, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let based: Vec<bool> = (0..n).map(|k| lake[k] != 0 || dep[k]).collect();
    let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let links = trunk_links(&sk);
    let trunk10: Vec<usize> = (0..n).filter(|&k| sk.trunk[k] && sk.area_km2[k] >= 10.0).collect();
    let lu: Vec<bool> = links.iter().map(|l| union[l.0] || union[l.1]).collect();
    let lb: Vec<bool> = links.iter().map(|l| based[l.0] || based[l.1]).collect();
    eprintln!(
        "   masks · union of the final lakes {} cells · based set {} cells · links touching: union {} · based {} · both {} (of {})",
        union.iter().filter(|&&b| b).count(),
        based.iter().filter(|&&b| b).count(),
        lu.iter().filter(|&&b| b).count(),
        lb.iter().filter(|&&b| b).count(),
        (0..links.len()).filter(|&i| lu[i] && lb[i]).count(),
        links.len()
    );
    let fmt = |r: (f32, f32, f32, usize)| format!("**{:.3}** [{:.3}, {:.3}] ({})", r.0, r.1, r.2, r.3);
    let mut carved_both = (Vec::new(), Vec::new());
    for (label, vc) in [("OFF", off), ("ON F132", f132), ("ON EXTENDED", ext)] {
        let skv = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (built, mk) = carve(&s1, &skv, &vc, &ss);
        let zb = metres(&built);
        drop(built);
        let zk: Vec<f32> = links.iter().map(|l| zb[l.0]).collect();
        let zr: Vec<f32> = links.iter().map(|l| zb[l.1]).collect();
        let cl: Vec<bool> = links.iter().map(|l| mk.carved[l.0] && mk.carved[l.1]).collect();
        eprintln!("   ── {label} (construction, carve(S1)) ──");
        for (nm, keep) in [
            ("all", &(|_: usize| true) as &dyn Fn(usize) -> bool),
            ("without the final lakes (union)", &(|i: usize| !lu[i]) as &dyn Fn(usize) -> bool),
            ("without the BASED set", &(|i: usize| !lb[i]) as &dyn Fn(usize) -> bool),
            ("without both", &(|i: usize| !lu[i] && !lb[i]) as &dyn Fn(usize) -> bool),
        ] {
            let ta = theta_links(&links, &zk, &zr, keep);
            let tc = theta_links(&links, &zk, &zr, &|i| keep(i) && cl[i]);
            eprintln!("      θ {nm:<32} all trunk links {} · carved links {}", fmt(ta), fmt(tc));
        }
        let under: Vec<bool> = trunk10
            .iter()
            .map(|&k| {
                let law = if skv.chi_m[k].is_finite() { skv.floor_m(k, vc.age_k) } else { f32::NAN };
                !mk.carved[k] && law.is_finite() && zs1[k] < law
            })
            .collect();
        let share = |sel: &dyn Fn(usize) -> bool| {
            let t = (0..trunk10.len()).filter(|&i| sel(trunk10[i])).count();
            let u = (0..trunk10.len()).filter(|&i| sel(trunk10[i]) && under[i]).count();
            format!("{:.1} % ({u} / {t})", 100.0 * u as f32 / t.max(1) as f32)
        };
        eprintln!(
            "      under-law share: all {} · without union {} · without BASED {} · without both {}",
            share(&|_| true),
            share(&|k| !union[k]),
            share(&|k| !based[k]),
            share(&|k| !union[k] && !based[k])
        );
        if label == "OFF" {
            carved_both.0 = mk.carved.clone();
        } else if label == "ON EXTENDED" {
            carved_both.1 = mk.carved.clone();
        }
    }
    // K3 — lake 2: distance upstream of its FOOTPRINT along OFF's drainage
    let target: Vec<bool> = (0..n).map(|k| lm_ext[k] == 2).collect();
    let mut dist = vec![f32::NAN; n];
    let mut done = vec![false; n];
    let mut path = Vec::new();
    for s in 0..n {
        if done[s] {
            continue;
        }
        path.clear();
        let mut c = s;
        let base;
        loop {
            if done[c] {
                base = dist[c];
                break;
            }
            if target[c] {
                done[c] = true;
                dist[c] = 0.0;
                base = 0.0;
                break;
            }
            path.push(c);
            let d = dir_off[c];
            if d == DIR_NONE || path.len() > 20000 {
                base = f32::NAN;
                break;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
                base = f32::NAN;
                break;
            }
            c = y as usize * w + x as usize;
        }
        let mut acc = base;
        for &p in path.iter().rev() {
            let d = dir_off[p];
            acc = if acc.is_finite() && d != DIR_NONE { acc + CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 } } else { f32::NAN };
            dist[p] = acc;
            done[p] = true;
        }
    }
    eprintln!("   K3 · lake 2 · Δz = EXT − OFF (conditioned) by OFF-drainage distance upstream of its FOOTPRINT (the footprint excluded):");
    for b in 0..15 {
        let (lo, hi) = (2.0 * b as f32, 2.0 * (b + 1) as f32);
        let sel: Vec<usize> = (0..n).filter(|&k| dist[k] > lo.max(1e-6) - 1e-6 && dist[k] > 0.0 && dist[k] >= lo && dist[k] < hi).collect();
        let mut dz: Vec<f32> = sel.iter().map(|&k| cn[k] - co[k]).collect();
        let mut dzc: Vec<f32> = sel.iter().filter(|&&k| carved_both.0[k] && carved_both.1[k]).map(|&k| cn[k] - co[k]).collect();
        let nc = dzc.len();
        eprintln!(
            "      {lo:>4.0}–{hi:<3.0} km · {:>7} cells · p50 Δz **{:>7.1} m** (p10 {:>7.1}, p90 {:>7.1}) · carved in both {nc:>6} · their p50 Δz **{:>7.1} m**",
            sel.len(),
            q_(&mut dz.clone(), 0.5),
            q_(&mut dz.clone(), 0.1),
            q_(&mut dz, 0.9),
            q_(&mut dzc, 0.5)
        );
    }
    eprintln!("\n==========  end Finding 135-K3 / M2 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 136 — the steps below the lakes' outlets (K5 profiles on the ON water's path, K6 inventory, K7 the
/// notch and the lowering), the bases on absent lakes (L), the contexts available without geology (K8). The
/// realigned témoin (viz state + run_hd tail). Declared in `f136_declared.md` before the run. SVGs to `F136_DIR`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f136_k --nocapture
#[test]
#[ignore]
fn f136_k() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::HashMap;
    use ymir_core::seed::WorldSeed;
    use ymir_core::tectonics_c1::bench_guard::{check_lakes, field_hash};
    use ymir_core::tectonics_c1::closures::lithology::stamp_volcanic_k;
    use ymir_core::tectonics_c1::closures::volcanism::place_edifices;
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, carve_diag, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 136 . the outlet steps, the absent bases, the contexts (realigned témoin)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let f132 = ValleyConstruction { lake_base: Some(LakeBase::InputLakes), ..off };
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    // the worlds
    let (co, lm_off) = {
        let wd = build_world(kn(off), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(off), PSEED, 45.0, 40.0);
        eprintln!(
            "   OFF · field {:016x} · lakes {:?}",
            field_hash(&wd.heightmap),
            check_lakes(&v.digest, &v.drainage.lake_map, &v.drainage.lakes)
        );
        (metres(&v.conditioned), v.drainage.lake_map.clone())
    };
    let lm_f132 = {
        let wd = build_world(kn(f132), None, PSEED, None);
        viz_hd_lakes_on(&wd, kn(f132), PSEED, 45.0, 40.0).drainage.lake_map
    };
    let wd_on = build_world(kn(ext), None, PSEED, None);
    let v_on = viz_hd_lakes_on(&wd_on, kn(ext), PSEED, 45.0, 40.0);
    let (w, h) = (wd_on.heightmap.width, wd_on.heightmap.height);
    let n = w * h;
    let eo_n = metres(&wd_on.heightmap);
    let pre_lm_on = c1_drainage_windowed(&wd_on.heightmap, None, &viz_dcfg(), &ss, DOMAIN_KM).lake_map;
    let cn = metres(&v_on.conditioned);
    let lm_on = v_on.drainage.lake_map.clone();
    let lakes_on = v_on.drainage.lakes.clone();
    let dir_on = v_on.drainage.flow.direction.clone();
    let acc_on = v_on.drainage.flow.accumulation.data.clone();
    eprintln!("   ON extended · field {:016x} · {} final lakes", field_hash(&wd_on.heightmap), lakes_on.len());
    // the construction input, the input bodies, the ON construction
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let zs1 = metres(&s1);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let recv = |c: usize, dir: &[u8]| -> Option<usize> {
        let d = dir[c];
        if d == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
    let (lake, dep, zbf) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep, metres(&bf))
    };
    let sk_off = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let skn = skeleton(&s1, &ext, &ss, DOMAIN_KM);
    let (built, mk_on, dg_on) = carve_diag(&s1, &skn, &ext, &ss);
    let zbn = metres(&built);
    let fill_b: Vec<f32> = {
        let f = compute_flow(&built, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
        (0..n).map(|k| (f.filled.data[k] - built.data[k]) * n2m).collect()
    };
    drop(built);
    let src: Vec<(f32, f32, f32, f32)> = skn.polylines.iter().flatten().copied().collect();
    let tan = ext.wall_deg.to_radians().tan();
    let geo = |s: usize, c: usize| -> f32 {
        let (sx, sy, zf, hw) = src[s];
        let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
        let mut dx = (cx - sx).rem_euclid(w as f32);
        if dx > w as f32 / 2.0 {
            dx -= w as f32;
        }
        let mut dy = (cy - sy).rem_euclid(h as f32);
        if dy > h as f32 / 2.0 {
            dy -= h as f32;
        }
        let d = (dx * dx + dy * dy).sqrt() * skn.cell_m;
        zf + tan * ((d - hw).max(0.0) - 0.0)
    };
    // who laid a constructed cell: (nearest?, the laying sample, its line's area at the sample)
    let laid = |c: usize| -> Option<(bool, usize, f32)> {
        if !mk_on.carved[c] {
            return None;
        }
        let me = dg_on.who[c] as usize;
        let (mut v, mut from) = (geo(me, c), me);
        let (x, y) = ((c % w) as i32, (c / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nb = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                let o = dg_on.who[nb];
                if o != u32::MAX && dg_on.line_of[o as usize] != dg_on.line_of[me] {
                    let vv = geo(o as usize, c);
                    if vv < v {
                        v = vv;
                        from = o as usize;
                    }
                }
            }
        }
        let _ = c1_metres_to_altitude_norm(v, &ss);
        let (sx, sy, _, _) = src[from];
        let sc = (sy.floor() as i64).rem_euclid(h as i64) as usize * w + (sx.floor() as i64).rem_euclid(w as i64) as usize;
        Some((from == me, from, skn.area_km2[sc]))
    };
    // input bodies
    let mut bodies: Vec<(String, Vec<usize>, bool)> = Vec::new();
    {
        let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
        for k in 0..n {
            if lake[k] != 0 {
                by.entry(lake[k]).or_default().push(k);
            }
        }
        let mut ids: Vec<u32> = by.keys().copied().collect();
        ids.sort_unstable();
        for id in ids {
            bodies.push((format!("input lake {id}"), by.remove(&id).unwrap(), false));
        }
        let mut seen = vec![false; n];
        let mut cid = 0;
        for s in 0..n {
            if !dep[s] || seen[s] {
                continue;
            }
            cid += 1;
            let mut comp = vec![s];
            seen[s] = true;
            let mut i = 0;
            while i < comp.len() {
                let c = comp[i];
                i += 1;
                for k in 0..8 {
                    let m = nbt(c, k);
                    if dep[m] && !seen[m] {
                        seen[m] = true;
                        comp.push(m);
                    }
                }
            }
            comp.sort_unstable();
            bodies.push((format!("closed depression {cid}"), comp, true));
        }
    }
    let ring_of = |cells: &[usize]| -> Vec<usize> {
        let mut r: Vec<usize> = Vec::new();
        for &c in cells {
            for k in 0..8 {
                let m = nbt(c, k);
                if cells.binary_search(&m).is_err() {
                    r.push(m);
                }
            }
        }
        r.sort_unstable();
        r.dedup();
        r
    };
    let mut body_of = vec![u32::MAX; n];
    let mut lin = Vec::new();
    let mut notch = Vec::new();
    for (bi, (_, cells, on_bf)) in bodies.iter().enumerate() {
        for &c in cells {
            body_of[c] = bi as u32;
        }
        let ring = ring_of(cells);
        let zr = |c: usize| if *on_bf { zbf[c] } else { zs1[c] };
        let col = *ring.iter().min_by(|&&a, &&b| zr(a).total_cmp(&zr(b))).expect("ring");
        lin.push(zr(col));
        let (mc, mz) = ring.iter().map(|&c| (c, zbn[c])).min_by(|a, b| a.1.total_cmp(&b.1)).expect("ring");
        notch.push((zr(col) - mz, mc));
    }
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    // ── K8 inputs first (K6's steps are classified with them)
    let dsea = dist_from(&cn.iter().map(|&z| z <= 0.0).collect::<Vec<_>>(), w, h);
    let cw: Vec<bool> = (0..n).map(|k| mk_on.carved[k] && !mk_on.floor[k] && cn[k] > 0.0 && dsea[k] <= 3).collect();
    let dwall = dist_from(&cw, w, h);
    drop(cw);
    let volc_mask: Vec<bool> = {
        let edi = place_edifices(&wd_on.state, &wd_on.kin, &WorldSeed::new(PSEED), wd_on.volc.domain_km, &wd_on.volc);
        let mut k = vec![1.0f32; n];
        let kpc = wd_on.cfg.sample_size as f32 * wd_on.volc.domain_km / w as f32;
        stamp_volcanic_k(&mut k, &edi, wd_on.cfg.sample_origin, wd_on.cfg.sample_size, kpc, w, h, &wd_on.cfg.lithology);
        k.iter().map(|&v| v != 1.0).collect()
    };
    // hanging junctions
    let area = |c: usize| acc_on[c] * cell_km2;
    let mut donors: HashMap<usize, Vec<usize>> = HashMap::new();
    for c in 0..n {
        if area(c) >= 10.0 {
            if let Some(r) = recv(c, &dir_on) {
                donors.entry(r).or_default().push(c);
            }
        }
    }
    let mut junctions: Vec<(usize, usize, f32, &'static str)> = Vec::new(); // (junction, tributary, H, kind)
    let mut jkeys: Vec<usize> = donors.keys().copied().collect();
    jkeys.sort_unstable();
    for r in jkeys {
        let ds = &donors[&r];
        if ds.len() < 2 {
            continue;
        }
        let tr = *ds.iter().max_by(|&&a, &&b| area(a).total_cmp(&area(b))).unwrap();
        for &t in ds {
            if t == tr {
                continue;
            }
            let hgt = cn[t] - cn[r];
            let kind = match laid(r) {
                None => "uncarved",
                Some((_, _, la)) if la <= area(t) => "dam-type",
                Some(_) => "clean break",
            };
            junctions.push((r, t, hgt, kind));
        }
    }
    let mut hang = vec![false; n];
    for &(r, _, hg, _) in &junctions {
        if hg > 50.0 {
            hang[r] = true;
        }
    }
    let dhang = dist_from(&hang, w, h);
    // coastal-wall mouths
    let mut mouths: Vec<(usize, f32)> = Vec::new();
    for (i, s) in v_on.drainage.rivers.segments.iter().enumerate() {
        if v_on.drainage.segment_kind[i] != SegmentKind::Watercourse || s.points.is_empty() {
            continue;
        }
        let (x, y) = *s.points.last().unwrap();
        let p = y as usize * w + x as usize;
        if let Some(r) = recv(p, &dir_on) {
            if cn[r] <= 0.0 && cn[p] > 1.0 && dwall[p] <= 2 {
                mouths.push((p, cn[p]));
            }
        }
    }
    let mut mouth_m = vec![false; n];
    for &(p, _) in &mouths {
        mouth_m[p] = true;
    }
    let dmouth = dist_from(&mouth_m, w, h);
    // ── K5 / K6 / K7 per present lake
    let dir_out = std::env::var("F136_DIR").ok();
    let present: Vec<&ymir_core::tectonics_c1::drainage::C1Lake> =
        lakes_on.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).collect();
    eprintln!("\n   K6 · present lakes (ON final, ≥ 1 km², not crater): {}", present.len());
    let mut steps: Vec<(u32, f32, f32)> = Vec::new(); // (lake, H_on, H_off)
    let mut step_ctx: Vec<String> = Vec::new();
    for l in &present {
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        let Some(col) = recv(o, &dir_on) else {
            eprintln!("      lake {} · outlet ({}, {}) has no receiver", l.base.id, l.base.outlet.0, l.base.outlet.1);
            continue;
        };
        let mut path = vec![col];
        let mut d = vec![0f32];
        let mut why = "60 km";
        let mut c = col;
        while *d.last().unwrap() < 60.0 {
            if cn[c] <= 0.0 {
                why = "the sea";
                break;
            }
            match recv(c, &dir_on) {
                None => {
                    why = "no receiver / edge";
                    break;
                }
                Some(r) => {
                    if lm_on[r] != 0 && lm_on[r] != l.base.id {
                        why = "another lake";
                        path.push(r);
                        d.push(d.last().unwrap() + step_len(c, r));
                        break;
                    }
                    d.push(d.last().unwrap() + step_len(c, r));
                    path.push(r);
                    c = r;
                }
            }
        }
        let zon: Vec<f32> = path.iter().map(|&p| cn[p]).collect();
        let zof: Vec<f32> = path.iter().map(|&p| co[p]).collect();
        let slope_over = |z: &[f32], i: usize, km: f32| -> Option<f32> {
            let j = (i..path.len()).find(|&j| d[j] >= d[i] + km)?;
            Some((z[i] - z[j]) / ((d[j] - d[i]) * 1000.0))
        };
        let step = |z: &[f32], top: f32| -> (f32, f32, f32, Option<(usize, usize)>) {
            let mut zone = None;
            if let Some(s) = (0..path.len()).find(|&i| d[i] <= 5.0 && slope_over(z, i, 1.0).is_some_and(|v| v >= 0.05)) {
                let e = (s + 1..path.len()).find(|&i| slope_over(z, i, 1.0).is_none_or(|v| v < 0.02)).unwrap_or(path.len() - 1);
                zone = Some((s, e));
            }
            let s2 = (0..path.len()).filter(|&i| d[i] <= 20.0).filter_map(|i| slope_over(z, i, 2.0)).fold(0f32, f32::max);
            match zone {
                Some((_, e)) => (top - z[e], d[e], s2, zone),
                None => (0.0, 0.0, s2, None),
            }
        };
        let (h_on, l_on, s2_on, zone_on) = step(&zon, l.level_m);
        let (h_off, l_off, s2_off, _) = step(&zof, zof[0]);
        // K7
        let mut ov: HashMap<u32, usize> = HashMap::new();
        for k in 0..n {
            if lm_on[k] == l.base.id && body_of[k] != u32::MAX {
                *ov.entry(body_of[k]).or_insert(0) += 1;
            }
        }
        let body = ov.into_iter().max_by_key(|e| e.1).map(|e| e.0 as usize);
        let k7 = match body {
            None => "no input body".to_string(),
            Some(b) => {
                let li = lin[b];
                let stage = if zs1[col] < li - 1.0 {
                    "S1 (the input body's own bed)".to_string()
                } else if zbn[col] < li - 1.0 {
                    match laid(col) {
                        Some((true, _, la)) => format!("the construction, NEAREST sample (line {la:.2} km²)"),
                        Some((false, _, la)) => format!("the construction, CROSS-LINE minimum (line {la:.2} km²)"),
                        None => "the construction (uncarved?)".to_string(),
                    }
                } else if eo_n[col] < li - 1.0 {
                    "the light pass / erosion".to_string()
                } else if cn[col] < li - 1.0 {
                    "the breach".to_string()
                } else {
                    "none (the sill is at L_in)".to_string()
                };
                let (nd, nc) = notch[b];
                let dist_nc = ((nc % w) as f32 - (col % w) as f32).hypot((nc / w) as f32 - (col / w) as f32) * CELL_KM;
                format!(
                    "{} · L_in {li:.1} m → level {:.1} m · **lowering {:.1} m** · true col ({}, {}) z S1 {:.1} / built {:.1} / eroded {:.1} / conditioned {:.1} → **{stage}** · notch (ring cut) **{nd:.1} m**, {dist_nc:.2} km from the true col",
                    bodies[b].0,
                    l.level_m,
                    li - l.level_m,
                    col % w,
                    col / w,
                    zs1[col],
                    zbn[col],
                    eo_n[col],
                    cn[col]
                )
            }
        };
        // contexts of the step
        let mut ctx = Vec::new();
        if let Some((s, e)) = zone_on {
            let near = |dd: &[u16]| path[s..=e].iter().any(|&p| (dd[p] as f32) * CELL_KM <= 1.0);
            if near(&dhang) {
                ctx.push("hanging junction");
            }
            if near(&dmouth) {
                ctx.push("coastal-wall mouth");
            }
            if path[s..=e].iter().any(|&p| volc_mask[p]) {
                ctx.push("volcanic");
            }
        }
        if h_on > 50.0 {
            steps.push((l.base.id, h_on, h_off));
            step_ctx.push(if ctx.is_empty() { "lake outlet only".to_string() } else { format!("lake outlet + {}", ctx.join(" + ")) });
        }
        eprintln!(
            "      lake {:>7} {:>7.1} km² level {:>7.1} m · path {:>5.1} km ({why}) · ON step **H {h_on:>6.1} m** over {l_on:>5.2} km, max S(2 km) {s2_on:.3} · OFF on the same cells H {h_off:>6.1} m over {l_off:.2} km, max S(2 km) {s2_off:.3} · ctx {}",
            l.base.id,
            l.area_km2,
            l.level_m,
            d.last().unwrap(),
            if ctx.is_empty() { "-".to_string() } else { ctx.join("+") }
        );
        eprintln!("         K7 · {k7}");
        // K5 profiles
        if let (Some(dout), true) = (&dir_out, [1u32, 2, 11].contains(&l.base.id)) {
            // K2's OFF-skeleton path (Finding 135), downstream of its point nearest the S1 col
            let cells: Vec<usize> = (0..n).filter(|&k| lm_on[k] == l.base.id).collect();
            let ring = ring_of(&cells);
            let colk2 = *ring.iter().min_by(|&&a, &&b| zs1[a].total_cmp(&zs1[b])).unwrap();
            let exit = *cells.iter().chain(ring.iter()).max_by(|&&a, &&b| sk_off.area_km2[a].total_cmp(&sk_off.area_km2[b])).unwrap();
            let mut k2 = vec![exit];
            let mut c = exit;
            let mut len = 0.0;
            while len < 60.0 {
                let dd = sk_off.direction[c];
                if dd == DIR_NONE || s1.data[c] <= SEA {
                    break;
                }
                let r = nbt(c, dd as usize);
                len += step_len(c, r);
                k2.push(r);
                c = r;
            }
            let mut xk = vec![0f32; k2.len()];
            for i in 1..k2.len() {
                xk[i] = xk[i - 1] + step_len(k2[i - 1], k2[i]);
            }
            let ci = (0..k2.len()).min_by(|&a, &b| {
                let dd = |p: usize| ((p % w) as f32 - (colk2 % w) as f32).hypot((p / w) as f32 - (colk2 / w) as f32);
                dd(k2[a]).total_cmp(&dd(k2[b]))
            }).unwrap();
            let ser = |z: &[f32], pth: &[usize], x: &[f32], x0: f32| -> Vec<(f32, f32)> { pth.iter().enumerate().map(|(i, &p)| (x[i] - x0, z[p])).collect() };
            let svg = profile_svg(
                &format!("Lake {} — downstream of the outlet: the ON water's own path (solid) vs Finding 135's OFF-skeleton path (dashed)", l.base.id),
                &[
                    ("ON, on the ON path", "#1f77b4", "", ser(&cn, &path, &d, 0.0)),
                    ("OFF, on the ON path", "#d62728", "", ser(&co, &path, &d, 0.0)),
                    ("S1, on the ON path", "#999999", "2,3", ser(&zs1, &path, &d, 0.0)),
                    ("ON, on K2's path", "#1f77b4", "8,5", ser(&cn, &k2[ci..], &xk[ci..], xk[ci])),
                    ("OFF, on K2's path", "#d62728", "8,5", ser(&co, &k2[ci..], &xk[ci..], xk[ci])),
                ],
                Some(0.0),
                &format!(
                    "0 = the true col (ON path; ON step H {h_on:.0} m over {l_on:.2} km) / K2's point nearest the S1 col (K2 path). Lake {} level {:.1} m. Path end: {why}.",
                    l.base.id, l.level_m
                ),
            );
            let p = std::path::Path::new(dout).join(format!("k5_profile_lake{}.svg", l.base.id));
            std::fs::write(&p, svg).expect("svg");
            eprintln!("         K5 profile: {}", p.display());
        }
    }
    let cnt = |t: f32| steps.iter().filter(|s| s.1 > t).count();
    let cnt_off = |t: f32| steps.iter().filter(|s| s.2 > t).count();
    eprintln!(
        "   K6 · ON steps > 50 / > 200 / > 500 m: **{} / {} / {}** · (OFF on the same cells, among those: {} / {} / {}) · the canyons instrument sees 0 (W3: canyons 0 in ON extended)",
        cnt(50.0),
        cnt(200.0),
        cnt(500.0),
        cnt_off(50.0),
        cnt_off(200.0),
        cnt_off(500.0)
    );
    let mut cx: HashMap<String, usize> = HashMap::new();
    for c in &step_ctx {
        *cx.entry(c.clone()).or_insert(0) += 1;
    }
    eprintln!("   K8 · K6's steps (> 50 m) by context: {cx:?}");
    // ── L
    let absent: Vec<usize> = (0..bodies.len()).filter(|&b| (bodies[b].1.iter().filter(|&&k| lm_on[k] != 0).count() as f32) < 0.5 * bodies[b].1.len() as f32).collect();
    eprintln!("\n   L · {} input bodies absent at the end:", absent.len());
    let cell_m = CELL_KM * 1000.0;
    for &b in &absent {
        let cells = &bodies[b].1;
        let m = cells.len() as f32;
        let sh = |f: &dyn Fn(usize) -> bool| 100.0 * cells.iter().filter(|&&k| f(k)).count() as f32 / m;
        let (s_b, s_e, s_f) = (sh(&|k| fill_b[k] > 0.5), sh(&|k| pre_lm_on[k] != 0), sh(&|k| lm_on[k] != 0));
        let stage = if s_b < 50.0 {
            "the CONSTRUCTION (no bowl left in carve(S1))"
        } else if s_e < 50.0 {
            "the light pass / erosion (no lake in the eroded world)"
        } else {
            "the BALANCE / cleanup (a lake in the eroded world, none at the end)"
        };
        let mut sl: Vec<f32> = cells
            .iter()
            .map(|&k| {
                let (x, y) = (k % w, k / w);
                let (xm, xp, ym, yp) = (x.saturating_sub(1), (x + 1).min(w - 1), y.saturating_sub(1), (y + 1).min(h - 1));
                let gx = (cn[y * w + xp] - cn[y * w + xm]) / (cell_m * (xp - xm).max(1) as f32);
                let gy = (cn[yp * w + x] - cn[ym * w + x]) / (cell_m * (yp - ym).max(1) as f32);
                (gx * gx + gy * gy).sqrt()
            })
            .collect();
        let flat = sl.iter().filter(|&&s| s < 0.001).count() as f32 / m * 100.0;
        let mut zz: Vec<f32> = cells.iter().map(|&k| cn[k]).collect();
        eprintln!(
            "      {:<24} {:>7} cells ({:.2} km²) · L_in {:.1} m · under water: S1 100 % / construction {s_b:.0} % / eroded {s_e:.0} % / final {s_f:.0} % → **{stage}** · final floor: slope p50 **{:.4}** · under 0.1 % **{flat:.1} %** · median z − L_in **{:.1} m** · ring cut ON {:.1} m",
            bodies[b].0,
            cells.len(),
            m * cell_km2,
            lin[b],
            q(&mut sl, 0.5),
            q(&mut zz, 0.5) - lin[b],
            notch[b].0
        );
    }
    // their share of F135-M's links
    let links = trunk_links(&sk_off);
    let based = |k: usize| lake[k] != 0 || dep[k];
    let uni = |k: usize| lm_off[k] != 0 || lm_f132[k] != 0 || lm_on[k] != 0;
    let mut is_abs = vec![false; n];
    for &b in &absent {
        for &k in &bodies[b].1 {
            is_abs[k] = true;
        }
    }
    let set: Vec<&(usize, usize, f64, f32)> = links.iter().filter(|l| (based(l.0) || based(l.1)) && !(uni(l.0) || uni(l.1))).collect();
    let in_abs = set.iter().filter(|l| is_abs[l.0] || is_abs[l.1]).count();
    let in_lake = set.iter().filter(|l| !(is_abs[l.0] || is_abs[l.1]) && (lake[l.0] != 0 || lake[l.1] != 0)).count();
    let in_dep = set.len() - in_abs - in_lake;
    eprintln!(
        "   L · F135-M's links in the based set but no final lake: **{}** (F135: 5 175) · touching the {} absent bodies **{in_abs} ({:.1} %)** · another input lake's cells (shrunk margins) {in_lake} ({:.1} %) · a closed depression only {in_dep} ({:.1} %)",
        set.len(),
        absent.len(),
        100.0 * in_abs as f32 / set.len().max(1) as f32,
        100.0 * in_lake as f32 / set.len().max(1) as f32,
        100.0 * in_dep as f32 / set.len().max(1) as f32
    );
    // ── K8 counts
    let jc = |t: f32| junctions.iter().filter(|j| j.2 > t).count();
    let jk = |t: f32, k: &str| junctions.iter().filter(|j| j.2 > t && j.3 == k).count();
    eprintln!(
        "\n   K8 · junctions with a tributary ≥ 10 km²: {} · hanging > 10 / > 50 / > 200 m: **{} / {} / {}** · of those > 10 m: dam-type {} · clean break {} · uncarved {}",
        junctions.len(),
        jc(10.0),
        jc(50.0),
        jc(200.0),
        jk(10.0, "dam-type"),
        jk(10.0, "clean break"),
        jk(10.0, "uncarved")
    );
    let mut top: Vec<&(usize, usize, f32, &str)> = junctions.iter().filter(|j| j.2 > 50.0).collect();
    top.sort_by(|a, b| b.2.total_cmp(&a.2));
    for j in top.iter().take(12) {
        eprintln!("      hanging {:.1} m at ({}, {}) · tributary {:.1} km² · trunk {:.1} km² · {}", j.2, j.0 % w, j.0 / w, area(j.1), area(j.0), j.3);
    }
    let mut mh: Vec<f32> = mouths.iter().map(|m| m.1).collect();
    eprintln!("   K8 · mouths above the sea (> 1 m) within 2 cells of a coastal wall: **{}** · height p50 {:.1} / max {:.1} m", mouths.len(), q(&mut mh.clone(), 0.5), q(&mut mh, 1.0));
    let slope2 = |c: usize| -> Option<f32> {
        let mut e = c;
        let mut len = 0.0;
        while len < 2.0 {
            let r = recv(e, &dir_on)?;
            len += step_len(e, r);
            e = r;
        }
        Some((cn[c] - cn[e]) / (len * 1000.0))
    };
    let (mut si, mut so) = (Vec::new(), Vec::new());
    for c in 0..n {
        if area(c) >= 10.0 && cn[c] > 0.0 {
            if let Some(s) = slope2(c) {
                if volc_mask[c] { si.push(s) } else { so.push(s) }
            }
        }
    }
    eprintln!(
        "   K8 · volcanic mask {} cells · trunk cells (≥ 10 km²) inside **{}** / outside {} · slope over 2 km inside p50 {:.4} / p90 {:.4} / max {:.3} · outside p50 {:.4} / p90 {:.4} / max {:.3}",
        volc_mask.iter().filter(|&&b| b).count(),
        si.len(),
        so.len(),
        q(&mut si.clone(), 0.5),
        q(&mut si.clone(), 0.9),
        q(&mut si, 1.0),
        q(&mut so.clone(), 0.5),
        q(&mut so.clone(), 0.9),
        q(&mut so, 1.0)
    );
    eprintln!("\n==========  end Finding 136 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 136-K6b — AMENDED instrument, declared after K5's profiles showed K6's step stopping at short
/// benches (< 1 km): the drop from the lake's level to the first graded reach at least 2 km long (slope over 2 km
/// < 2 %), searched within the first 10 km of the ON path; and the drop over the first 5 km. Same paths as K6.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f136_k6b --nocapture
#[test]
#[ignore]
fn f136_k6b() {
    use common::{build_world, viz_hd_lakes_on};
    use ymir_core::tectonics_c1::drainage::LakeType;
    use ymir_core::tectonics_c1::valley_construction::LakeBase;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let ss = SteinSteinParams::default();
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn(ext), None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn(ext), PSEED, 45.0, 40.0);
    let (w, h) = (wd.heightmap.width, wd.heightmap.height);
    let cn: Vec<f32> = v.conditioned.data.iter().map(|&x| c1_altitude_norm_to_metres(x, &ss)).collect();
    let (lm, dir) = (&v.drainage.lake_map, &v.drainage.flow.direction);
    let recv = |c: usize| -> Option<usize> {
        let d = dir[c];
        if d == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
    eprintln!("\n==========  Finding 136-K6b . the drop to the first graded reach ≥ 2 km (amended)  ==========");
    let mut hs = Vec::new();
    for l in v.drainage.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        let Some(col) = recv(o) else { continue };
        let (mut path, mut d) = (vec![col], vec![0f32]);
        let mut c = col;
        while *d.last().unwrap() < 30.0 && cn[c] > 0.0 {
            let Some(r) = recv(c) else { break };
            d.push(d.last().unwrap() + step_len(c, r));
            path.push(r);
            if lm[r] != 0 && lm[r] != l.base.id {
                break;
            }
            c = r;
        }
        let z: Vec<f32> = path.iter().map(|&p| cn[p]).collect();
        let s2 = |i: usize| -> Option<f32> {
            let j = (i..path.len()).find(|&j| d[j] >= d[i] + 2.0)?;
            Some((z[i] - z[j]) / ((d[j] - d[i]) * 1000.0))
        };
        let g = (0..path.len()).find(|&i| d[i] <= 10.0 && s2(i).is_some_and(|s| s < 0.02));
        let (hg, dg) = match g {
            Some(i) => (l.level_m - z[i], d[i]),
            None => (f32::NAN, f32::NAN),
        };
        let i5 = (0..path.len()).rev().find(|&i| d[i] <= 5.0).unwrap_or(0);
        eprintln!(
            "   lake {:>4} {:>6.1} km² level {:>7.1} m · drop to the first graded reach ≥ 2 km **{hg:>6.1} m**, reached at {dg:>5.2} km · drop over the first 5 km {:>6.1} m",
            l.base.id,
            l.area_km2,
            l.level_m,
            l.level_m - z[i5]
        );
        if hg.is_finite() {
            hs.push(hg);
        }
    }
    let c = |t: f32| hs.iter().filter(|&&x| x > t).count();
    eprintln!("   K6b · drops > 50 / > 200 / > 500 m: **{} / {} / {}** of {} lakes with a graded reach within 10 km", c(50.0), c(200.0), c(500.0), hs.len());
    eprintln!("\n==========  end Finding 136-K6b  ==========\n");
}

/// F136-K6b's amended step instrument on one world: for every final lake ≥ 1 km² that is not a crater lake and has a
/// D8 receiver at `Lake::outlet`, the drop from its level to the first graded reach ≥ 2 km (slope over 2 km < 2 %)
/// within 10 km of the true col, the distance to it, and the drop over the first 5 km. Returns (lake id, km², level,
/// drop, distance km, drop over 5 km) and the number of lakes without a D8 receiver.
fn steps_k6b(cn: &[f32], d: &ymir_core::tectonics_c1::drainage::C1DrainageResult, w: usize, h: usize) -> (Vec<(u32, f32, f32, f32, f32, f32)>, usize) {
    use ymir_core::tectonics_c1::drainage::LakeType;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let dir = &d.flow.direction;
    let recv = |c: usize| -> Option<usize> {
        let k = dir[c];
        if k == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
    let (mut out, mut none) = (Vec::new(), 0usize);
    for l in d.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        let Some(col) = recv(o) else {
            none += 1;
            continue;
        };
        let (mut path, mut dd) = (vec![col], vec![0f32]);
        let mut c = col;
        while *dd.last().unwrap() < 30.0 && cn[c] > 0.0 {
            let Some(r) = recv(c) else { break };
            dd.push(dd.last().unwrap() + step_len(c, r));
            path.push(r);
            if d.lake_map[r] != 0 && d.lake_map[r] != l.base.id {
                break;
            }
            c = r;
        }
        let (dg, xg, d5) = profile_k6b(&path.iter().map(|&p| cn[p]).collect::<Vec<_>>(), &dd, l.level_m);
        out.push((l.base.id, l.area_km2, l.level_m, dg, xg, d5));
    }
    (out, none)
}

/// K6b on an explicit profile (z along the path, cumulative distance in km): (drop to the first graded reach ≥ 2 km
/// within 10 km, the distance to it, the drop over the first 5 km), from `top`.
fn profile_k6b(z: &[f32], d: &[f32], top: f32) -> (f32, f32, f32) {
    let s2 = |i: usize| -> Option<f32> {
        let j = (i..z.len()).find(|&j| d[j] >= d[i] + 2.0)?;
        Some((z[i] - z[j]) / ((d[j] - d[i]) * 1000.0))
    };
    let g = (0..z.len()).find(|&i| d[i] <= 10.0 && s2(i).is_some_and(|s| s < 0.02));
    let i5 = (0..z.len()).rev().find(|&i| d[i] <= 5.0).unwrap_or(0);
    match g {
        Some(i) => (top - z[i], d[i], top - z[i5]),
        None => (f32::NAN, f32::NAN, top - z[i5]),
    }
}

/// Q3's classes on the mean slope of a descent: wall > 45°, steep 10–45°, gorge < 10°.
fn slope_class(drop: f32, km: f32) -> (f32, &'static str) {
    if !(drop.is_finite() && km > 0.0) || drop <= 0.0 {
        return (f32::NAN, "—");
    }
    let deg = (drop / (km * 1000.0)).atan().to_degrees();
    (deg, if deg > 45.0 { "WALL" } else if deg >= 10.0 { "steep" } else { "gorge" })
}

/// ADR Finding 137-Q — Q1 (the steps below the lakes in the delivered worlds, F136-K6b's instrument, against ON
/// extended), Q2 (C-3's `production_k_field`: components, histogram, correlation length, correlations, map), Q4 (the
/// below-sea basins' spillways, the same instrument). Declared in `f137_declared.md` before the run. Map to
/// `F137_DIR`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f137_q --nocapture
#[test]
#[ignore]
fn f137_q() {
    use common::{build_world, viz_hd_lakes_on};
    use ymir_core::seed::WorldSeed;
    use ymir_core::tectonics_c1::closures::fracture::build_hd_density_k;
    use ymir_core::tectonics_c1::closures::lithology::{build_coarse_k, stamp_volcanic_k, upscale_k_to_hd};
    use ymir_core::tectonics_c1::closures::volcanism::place_edifices;
    use ymir_core::tectonics_c1::production_upscale::production_k_field;
    use ymir_core::tectonics_c1::valley_construction::LakeBase;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 137-Q . the delivered worlds' steps, C-3, the below-sea spillways  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let worlds: [(&str, Knobs); 3] = [
        ("livré (delivered stream power)", Knobs::passes(2)),
        ("A1+B2 (delivered + the age closure, F114's world)", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
        ("ON extended (construction)", Knobs { valley: Some(ext), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
    ];
    let mut keep = None;
    for (label, kn) in worlds {
        let t = Instant::now();
        let wd = build_world(kn, None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
        let (w, h) = (wd.heightmap.width, wd.heightmap.height);
        let cn = metres(&v.conditioned);
        let (st, none) = steps_k6b(&cn, &v.drainage, w, h);
        let c = |t: f32| st.iter().filter(|s| s.3 > t).count();
        eprintln!(
            "\n   Q1 · {label}: {} final lakes · ≥ 1 km² non-crater with a D8 outlet {} (without {none}) · steps (drop to the first graded reach ≥ 2 km) > 50 / > 200 / > 500 m: **{} / {} / {}** ({:.0} s)",
            v.drainage.lakes.len(),
            st.len(),
            c(50.0),
            c(200.0),
            c(500.0),
            t.elapsed().as_secs_f64()
        );
        for s in &st {
            let (deg, cls) = slope_class(s.3, s.4);
            eprintln!(
                "      lake {:>7} {:>7.1} km² level {:>7.1} m · drop {:>6.1} m at {:>5.2} km · mean slope {deg:>5.1}° {cls:<5} · drop over 5 km {:>6.1} m",
                s.0, s.1, s.2, s.3, s.4, s.5
            );
        }
        if label.starts_with("ON") {
            keep = Some((wd, v, cn));
        }
    }
    let (wd, v, cn) = keep.expect("ON");
    let (w, h) = (wd.heightmap.width, wd.heightmap.height);
    let n = w * h;
    // ── Q4 — the below-sea basins' spillways
    eprintln!("\n   Q4 · the below-sea-basin lakes (ON extended, ≥ 1 km², no D8 outlet) and their spillways:");
    let dir = &v.drainage.flow.direction;
    let recv = |c: usize| -> Option<usize> {
        let k = dir[c];
        if k == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
    let segs = &v.drainage.rivers.segments;
    let mut q4 = Vec::new();
    for l in v.drainage.lakes.iter().filter(|l| l.area_km2 >= 1.0 && l.base.id >= 1_000_000 && l.base.id < 2_000_000) {
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        if recv(o).is_some() {
            continue;
        }
        let sp: Vec<usize> = (0..segs.len())
            .filter(|&i| v.drainage.segment_kind[i] == SegmentKind::Spillway && v.drainage.segment_source_lake.get(i).copied().flatten() == Some(l.base.id))
            .collect();
        if sp.is_empty() {
            eprintln!("      lake {} · {:.1} km² · level {:.1} m · NO spillway segment", l.base.id, l.area_km2, l.level_m);
            continue;
        }
        // the spillway's points, then its downstream chain, then the D8 until the sea / 30 km
        let mut path: Vec<usize> = Vec::new();
        let mut si = Some(sp[0]);
        let mut guard = 0;
        while let Some(i) = si {
            for &(x, y) in &segs[i].points {
                let p = y as usize * w + x as usize;
                if path.last() != Some(&p) {
                    path.push(p);
                }
            }
            si = segs[i].downstream;
            guard += 1;
            if guard > 10_000 {
                break;
            }
        }
        let mut c = *path.last().unwrap();
        let mut len: f32 = path.windows(2).map(|p| step_len(p[0], p[1])).sum();
        while len < 30.0 && cn[c] > 0.0 {
            let Some(r) = recv(c) else { break };
            len += step_len(c, r);
            path.push(r);
            c = r;
        }
        let mut dd = vec![0f32; path.len()];
        for i in 1..path.len() {
            dd[i] = dd[i - 1] + step_len(path[i - 1], path[i]);
        }
        let z: Vec<f32> = path.iter().map(|&p| cn[p]).collect();
        let (dg, xg, d5) = profile_k6b(&z, &dd, l.level_m);
        let (deg, cls) = slope_class(dg, xg);
        let zend = *z.last().unwrap();
        eprintln!(
            "      lake {} · {:>6.1} km² · {:?} · level {:>7.1} m · spillway start z {:.1} m · path {:.1} km, end z {zend:.1} m · drop to the first graded reach ≥ 2 km **{dg:.1} m** at {xg:.2} km · {deg:.1}° {cls} · drop over 5 km {d5:.1} m · total drop to the path end {:.1} m",
            l.base.id,
            l.area_km2,
            l.lake_type,
            l.level_m,
            z[0],
            dd.last().unwrap(),
            l.level_m - zend
        );
        q4.push(dg);
    }
    let c4 = |t: f32| q4.iter().filter(|&&x| x > t).count();
    eprintln!("   Q4 · spillways measured {} · drops > 50 / > 200 / > 500 m: **{} / {} / {}**", q4.len(), c4(50.0), c4(200.0), c4(500.0));
    // ── Q2 — C-3's production_k_field
    let edi = place_edifices(&wd.state, &wd.kin, &WorldSeed::new(PSEED), wd.volc.domain_km, &wd.volc);
    let kf = production_k_field(&wd.state, Some(&wd.kin), &wd.cfg, &edi, &wd.volc, w, h).expect("C-3 is on in the viz state");
    let lith: Vec<f32> = upscale_k_to_hd(&build_coarse_k(&wd.state, &wd.cfg.lithology), w, h, wd.cfg.sample_origin, wd.cfg.sample_size);
    let volc: Vec<bool> = {
        let mut k = vec![1.0f32; n];
        let kpc = wd.cfg.sample_size as f32 * wd.volc.domain_km / w as f32;
        stamp_volcanic_k(&mut k, &edi, wd.cfg.sample_origin, wd.cfg.sample_size, kpc, w, h, &wd.cfg.lithology);
        k.iter().map(|&x| x != 1.0).collect()
    };
    let frac = build_hd_density_k(&wd.state, &wd.kin, &wd.cfg.fracture, None, w, h, wd.cfg.sample_origin, wd.cfg.sample_size);
    let land: Vec<bool> = cn.iter().map(|&z| z > 0.0).collect();
    let nl = land.iter().filter(|&&b| b).count();
    eprintln!("\n   Q2 · C-3 production_k_field on the viz state ({} cells, {nl} land):", n);
    let stats = |name: &str, f: &dyn Fn(usize) -> f32| {
        let mut v: Vec<f32> = (0..n).filter(|&k| land[k]).map(f).collect();
        v.sort_by(f32::total_cmp);
        let q = |p: f64| v[((v.len() - 1) as f64 * p) as usize];
        eprintln!(
            "      {name:<34} land min {:.3} · p1 {:.3} · p10 {:.3} · p50 {:.3} · p90 {:.3} · p99 {:.3} · max {:.3}",
            v[0],
            q(0.01),
            q(0.1),
            q(0.5),
            q(0.9),
            q(0.99),
            v[v.len() - 1]
        );
    };
    stats("K (the product)", &|k| kf[k]);
    stats("lithology classes (bilinear)", &|k| lith[k]);
    stats("fracture density factor (C-3b)", &|k| frac[k]);
    let hist: Vec<(f32, f32, usize)> = [(0.0, 1.0), (1.0, 1.0001), (1.0001, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 6.0), (6.0, 9.99), (9.99, 10.01), (10.01, 100.0)]
        .iter()
        .map(|&(a, b)| (a, b, (0..n).filter(|&k| land[k] && kf[k] >= a && kf[k] < b).count()))
        .collect();
    eprintln!("      histogram of K on land: {}", hist.iter().map(|(a, b, c)| format!("[{a}, {b}) {:.2} %", 100.0 * *c as f64 / nl as f64)).collect::<Vec<_>>().join(" · "));
    eprintln!(
        "      volcanic mask {:.2} % of land · lithology rift-soft (> 1.5) {:.2} % · fracture factor ≠ 1 {:.2} %",
        100.0 * (0..n).filter(|&k| land[k] && volc[k]).count() as f64 / nl as f64,
        100.0 * (0..n).filter(|&k| land[k] && lith[k] > 1.5).count() as f64 / nl as f64,
        100.0 * (0..n).filter(|&k| land[k] && (frac[k] - 1.0).abs() > 1e-4).count() as f64 / nl as f64
    );
    // correlation length: the autocorrelation along x and y, over land pairs
    let acf = |f: &dyn Fn(usize) -> f32, lag: usize, horiz: bool| -> f64 {
        let (mut sx, mut sy, mut sxx, mut syy, mut sxy, mut m) = (0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
        for y in (0..h).step_by(8) {
            for x in (0..w).step_by(8) {
                let (x2, y2) = if horiz { (x + lag, y) } else { (x, y + lag) };
                if x2 >= w || y2 >= h {
                    continue;
                }
                let (a, b) = (y * w + x, y2 * w + x2);
                if !(land[a] && land[b]) {
                    continue;
                }
                let (u, v) = (f(a) as f64, f(b) as f64);
                sx += u;
                sy += v;
                sxx += u * u;
                syy += v * v;
                sxy += u * v;
                m += 1.0;
            }
        }
        let cov = sxy / m - (sx / m) * (sy / m);
        cov / ((sxx / m - (sx / m).powi(2)).sqrt() * (syy / m - (sy / m).powi(2)).sqrt()).max(1e-12)
    };
    let lags = [1usize, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
    for (name, f) in [
        ("K (the product)", &(|k: usize| kf[k]) as &dyn Fn(usize) -> f32),
        ("log K", &(|k: usize| kf[k].ln()) as &dyn Fn(usize) -> f32),
        ("lithology classes", &(|k: usize| lith[k]) as &dyn Fn(usize) -> f32),
        ("fracture factor", &(|k: usize| frac[k]) as &dyn Fn(usize) -> f32),
        ("altitude (conditioned)", &(|k: usize| cn[k]) as &dyn Fn(usize) -> f32),
    ] {
        let row: Vec<(usize, f64, f64)> = lags.iter().map(|&l| (l, acf(f, l, true), acf(f, l, false))).collect();
        let e = (1.0f64).exp().recip();
        let lx = row.iter().find(|r| r.1 < e).map_or(f32::NAN, |r| r.0 as f32 * CELL_KM);
        let ly = row.iter().find(|r| r.2 < e).map_or(f32::NAN, |r| r.0 as f32 * CELL_KM);
        eprintln!(
            "      ACF {name:<24} first lag below 1/e: x **{lx:.1} km** · y **{ly:.1} km** · ACF(x) at 1/4/16/64/256/1024 cells: {}",
            [0usize, 2, 4, 6, 8, 10].iter().map(|&i| format!("{:.2}", row[i].1)).collect::<Vec<_>>().join(" / ")
        );
    }
    // correlations on land
    let precip = c1_climate_placed(&v.conditioned, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM).precipitation;
    let cell_m = CELL_KM * 1000.0;
    let slope = |k: usize| -> f32 {
        let (x, y) = (k % w, k / w);
        let (xm, xp, ym, yp) = (x.saturating_sub(1), (x + 1).min(w - 1), y.saturating_sub(1), (y + 1).min(h - 1));
        let gx = (cn[y * w + xp] - cn[y * w + xm]) / (cell_m * (xp - xm).max(1) as f32);
        let gy = (cn[yp * w + x] - cn[ym * w + x]) / (cell_m * (yp - ym).max(1) as f32);
        (gx * gx + gy * gy).sqrt()
    };
    let crat: Vec<f32> = {
        // the cratonic mask, nearest-sampled at the SAME (sx, sy) mapping as the altitude
        let (cw, ch) = (wd.state.plate_id.nx(), wd.state.plate_id.ny());
        let cm = &wd.state.cratonic_mask;
        (0..n)
            .map(|k| {
                let (x, y) = (k % w, k / w);
                let sx = wd.cfg.sample_origin[0] * cw as f64 + x as f64 * wd.cfg.sample_size * cw as f64 / w as f64;
                let sy = wd.cfg.sample_origin[1] * ch as f64 + y as f64 * wd.cfg.sample_size * ch as f64 / h as f64;
                let (i, j) = ((sx.floor() as i64).rem_euclid(cw as i64) as usize, (sy.floor() as i64).rem_euclid(ch as i64) as usize);
                if cm.get(i, j) { 1.0 } else { 0.0 }
            })
            .collect()
    };
    let pearson = |a: &dyn Fn(usize) -> f32, b: &dyn Fn(usize) -> f32| -> f64 {
        let (mut sa, mut sb, mut saa, mut sbb, mut sab, mut m) = (0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
        for k in (0..n).step_by(7) {
            if !land[k] {
                continue;
            }
            let (u, v) = (a(k) as f64, b(k) as f64);
            sa += u;
            sb += v;
            saa += u * u;
            sbb += v * v;
            sab += u * v;
            m += 1.0;
        }
        (sab / m - sa / m * sb / m) / ((saa / m - (sa / m).powi(2)).sqrt() * (sbb / m - (sb / m).powi(2)).sqrt()).max(1e-12)
    };
    let lk = |k: usize| kf[k].ln();
    eprintln!(
        "      Pearson r of log K (land) with: altitude {:.3} · slope {:.3} · precipitation {:.3} · volcanic mask {:.3} · cratonic mask {:.3} · lithology rift-soft {:.3} · fracture factor {:.3}",
        pearson(&lk, &|k| cn[k]),
        pearson(&lk, &slope),
        pearson(&lk, &|k| precip.data[k]),
        pearson(&lk, &|k| volc[k] as u8 as f32),
        pearson(&lk, &|k| crat[k]),
        pearson(&lk, &|k| (lith[k] > 1.5) as u8 as f32),
        pearson(&lk, &|k| frac[k])
    );
    // the map (1/8 resolution, north up, log K)
    if let Ok(dout) = std::env::var("F137_DIR") {
        let s = 8usize;
        let (iw, ih) = (w / s, h / s);
        let mut img = image::RgbImage::new(iw as u32, ih as u32);
        let lmax = kf.iter().cloned().fold(1.0f32, f32::max).ln().max(1e-3);
        for yy in 0..ih {
            for xx in 0..iw {
                let k = (h - 1 - yy * s) * w + xx * s; // north up
                let px = if !land[k] {
                    [30, 40, 70]
                } else {
                    let t = (kf[k].ln() / lmax).clamp(0.0, 1.0);
                    let base = [(60.0 + 195.0 * t) as u8, (160.0 * (1.0 - t) + 40.0) as u8, (60.0 * (1.0 - t)) as u8];
                    if volc[k] { [255, 255, 255] } else { base }
                };
                img.put_pixel(xx as u32, yy as u32, image::Rgb(px));
            }
        }
        let p = std::path::Path::new(&dout).join("q2_k_field_map.png");
        img.save(&p).expect("png");
        eprintln!("      map: {} (north up, 1/8 resolution; land coloured by log K from green = 1 to red = max {:.2}; volcanic discs white; sea dark blue)", p.display(), lmax.exp());
    }
    eprintln!("\n==========  end Finding 137-Q . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// The index of F136-K6b's graded cell (the first cell within 10 km whose slope over the next 2 km is < 2 %).
fn graded_index(z: &[f32], d: &[f32]) -> Option<usize> {
    let s2 = |i: usize| -> Option<f32> {
        let j = (i..z.len()).find(|&j| d[j] >= d[i] + 2.0)?;
        Some((z[i] - z[j]) / ((d[j] - d[i]) * 1000.0))
    };
    (0..z.len()).find(|&i| d[i] <= 10.0 && s2(i).is_some_and(|s| s < 0.02))
}

/// The max slope over 100 m (degrees) on the cells `0..=end` of a path.
fn max_slope_100m(z: &[f32], d: &[f32], end: usize) -> f32 {
    let mut m = 0f32;
    for i in 0..=end.min(z.len() - 1) {
        if let Some(j) = (i..z.len()).find(|&j| d[j] >= d[i] + 0.1) {
            m = m.max(((z[i] - z[j]) / ((d[j] - d[i]) * 1000.0)).atan().to_degrees());
        }
    }
    m
}

/// ADR Finding 138-T — the tables for the author's decisions (read-only, no production code): T1 the gorge by the
/// steepness factor m, T2 the short slopes of today's steps, T3 the lakes by the retreat parameter r, T4 the falls at
/// m = 10. ON extended (the témoin) and, for T2, the delivered world. Declared in `f138_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f138_t --nocapture
#[test]
#[ignore]
fn f138_t() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{HashMap, VecDeque};
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let k_law = F121_AGE_K;
    eprintln!("\n==========  Finding 138-T . the tables for the decisions (ON extended, témoin)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    // the delivered world, for T2
    let (zl, dl) = {
        let wd = build_world(Knobs::passes(2), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, Knobs::passes(2), PSEED, 45.0, 40.0);
        (metres(&v.conditioned), v.drainage)
    };
    let kn = Knobs { valley: Some(ext), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let (w, h) = (wd.heightmap.width, wd.heightmap.height);
    let n = w * h;
    let cn = metres(&v.conditioned);
    let dr = &v.drainage;
    let level: HashMap<u32, f32> = dr.lakes.iter().map(|l| (l.base.id, l.level_m)).collect();
    let recv = |dir: &[u8], c: usize| -> Option<usize> {
        let k = dir[c];
        if k == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
    // ── T2 first (today's steps), on both worlds
    eprintln!("\n   T2 · the max slope over 100 m on today's steps (F136-K6b drop > 50 m):");
    for (name, z, d) in [("ON extended", &cn, dr), ("livré", &zl, &dl)] {
        let mut rows = Vec::new();
        for l in d.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
            let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
            let Some(col) = recv(&d.flow.direction, o) else { continue };
            let (mut path, mut dd) = (vec![col], vec![0f32]);
            let mut c = col;
            while *dd.last().unwrap() < 30.0 && z[c] > 0.0 {
                let Some(r) = recv(&d.flow.direction, c) else { break };
                dd.push(dd.last().unwrap() + step_len(c, r));
                path.push(r);
                if d.lake_map[r] != 0 && d.lake_map[r] != l.base.id {
                    break;
                }
                c = r;
            }
            let zz: Vec<f32> = path.iter().map(|&p| z[p]).collect();
            let Some(g) = graded_index(&zz, &dd) else { continue };
            let drop = l.level_m - zz[g];
            if drop <= 50.0 {
                continue;
            }
            let ms = max_slope_100m(&zz, &dd, g);
            rows.push((l.base.id, drop, ms));
        }
        let c = |t: f32, big: bool| rows.iter().filter(|r| r.2 > t && (!big || r.1 > 200.0)).count();
        eprintln!(
            "      {name}: {} steps > 50 m · max slope over 100 m > 28° / > 40° / > 60°: **{} / {} / {}** · among the {} > 200 m: **{} / {} / {}**",
            rows.len(),
            c(28.0, false),
            c(40.0, false),
            c(60.0, false),
            rows.iter().filter(|r| r.1 > 200.0).count(),
            c(28.0, true),
            c(40.0, true),
            c(60.0, true)
        );
        for r in &rows {
            eprintln!("         lake {:>4} · drop {:>6.1} m · max slope over 100 m **{:>5.1}°**", r.0, r.1, r.2);
        }
    }
    // ── the paths and the next bases (ON extended)
    struct P {
        id: u32,
        km2: f32,
        lvl: f32,
        path: Vec<usize>,
        d: Vec<f32>,
        base: f32,
        kind: String,
    }
    let mut ps: Vec<P> = Vec::new();
    let extend = |path: &mut Vec<usize>, d: &mut Vec<f32>, me: u32| -> (f32, String) {
        let mut c = *path.last().unwrap();
        loop {
            if cn[c] <= 0.0 {
                return (0.0, "the sea".into());
            }
            if *d.last().unwrap() >= 30.0 {
                return (cn[c], "30 km (no base)".into());
            }
            let Some(r) = recv(&dr.flow.direction, c) else { return (cn[c], "no receiver (no base)".into()) };
            d.push(d.last().unwrap() + step_len(c, r));
            path.push(r);
            if dr.lake_map[r] != 0 && dr.lake_map[r] != me {
                let id = dr.lake_map[r];
                return (level.get(&id).copied().unwrap_or(cn[r]), format!("lake {id}"));
            }
            c = r;
        }
    };
    let segs = &dr.rivers.segments;
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        if let Some(col) = recv(&dr.flow.direction, o) {
            let (mut path, mut d) = (vec![col], vec![0f32]);
            let (base, kind) = extend(&mut path, &mut d, l.base.id);
            ps.push(P { id: l.base.id, km2: l.area_km2, lvl: l.level_m, path, d, base, kind });
        } else if l.base.id >= 1_000_000 && l.base.id < 2_000_000 {
            let sp: Vec<usize> = (0..segs.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway && dr.segment_source_lake.get(i).copied().flatten() == Some(l.base.id))
                .collect();
            if sp.is_empty() {
                continue;
            }
            let mut path: Vec<usize> = Vec::new();
            let mut si = Some(sp[0]);
            let mut guard = 0;
            while let Some(i) = si {
                for &(x, y) in &segs[i].points {
                    let p = y as usize * w + x as usize;
                    if path.last() != Some(&p) {
                        path.push(p);
                    }
                }
                si = segs[i].downstream;
                guard += 1;
                if guard > 10_000 {
                    break;
                }
            }
            let mut d = vec![0f32; path.len()];
            for i in 1..path.len() {
                d[i] = d[i - 1] + step_len(path[i - 1], path[i]);
            }
            // a spillway that already enters another lake or the sea
            let hit = (1..path.len()).find(|&i| (dr.lake_map[path[i]] != 0 && dr.lake_map[path[i]] != l.base.id) || cn[path[i]] <= 0.0);
            let (base, kind) = match hit {
                Some(i) => {
                    path.truncate(i + 1);
                    d.truncate(i + 1);
                    let q = path[i];
                    if cn[q] <= 0.0 { (0.0, "the sea".to_string()) } else { (level.get(&dr.lake_map[q]).copied().unwrap_or(cn[q]), format!("lake {}", dr.lake_map[q])) }
                }
                None => extend(&mut path, &mut d, l.base.id),
            };
            ps.push(P { id: l.base.id, km2: l.area_km2, lvl: l.level_m, path, d, base, kind: format!("spillway → {kind}") });
        }
    }
    // the construction's law below
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let skn = skeleton(&s1, &ext, &ss, DOMAIN_KM);
    let acc = &dr.flow.accumulation.data;
    // ── T1 and T4
    let tan28 = 28f32.to_radians().tan();
    let phi_of = |id: u32| -> f32 {
        // splitmix64 of (seed, lake id): FNV-1a's high bits barely moved with the id (first run: φ ≈ 0.45 everywhere)
        let mut z = PSEED ^ (id as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        let u = (z >> 11) as f64 / (1u64 << 53) as f64;
        if u < 1.0 / 3.0 { 0.0 } else { (0.1 + 0.4 * (u - 1.0 / 3.0) / (2.0 / 3.0)) as f32 }
    };
    eprintln!("\n   T1 · the gorge by m (L = today's level; S_g = min(m·S_loi(A), tan 28°)):");
    let mut t1 = [(0usize, 0usize); 3];
    let mut t4: Vec<(u32, f32, f32, f32, f32, bool, bool)> = Vec::new(); // id, phi, head, shortage, H_f, head tagged, shortage tagged
    let mut gorge_viol = 0usize;
    for p in &ps {
        // AMENDED (after the first run): A = max(the accumulation on the path's first km, the lake's inflow), and the law
        // below = B + k·χ along the path FROM THE NEXT BASE (the skeleton's floor_m read the lake's own base on the bed cells)
        let inflow: f32 = (0..n)
            .filter(|&c| dr.lake_map[c] != p.id && recv(&dr.flow.direction, c).is_some_and(|r| dr.lake_map[r] == p.id))
            .map(|c| acc[c] * cell_km2)
            .sum::<f32>()
            + p.km2;
        let a = p.path.iter().zip(p.d.iter()).filter(|e| *e.1 <= 1.0).map(|e| acc[*e.0] * cell_km2).fold(0f32, f32::max).max(inflow).max(0.1);
        let sl = k_law * a.powf(-0.5);
        let dtot = p.lvl - p.base;
        let x = *p.d.last().unwrap();
        let mut chi = vec![0f32; p.path.len()];
        for i in (0..p.path.len().saturating_sub(1)).rev() {
            let ai = (acc[p.path[i]] * cell_km2).max(a).max(0.1);
            chi[i] = chi[i + 1] + (1.0 / ai).sqrt() * (p.d[i + 1] - p.d[i]) * 1000.0;
        }
        let law = |i: usize| -> f32 { p.base + k_law * chi[i] };
        let mut cols = Vec::new();
        for (mi, m) in [3f32, 10.0, 30.0].into_iter().enumerate() {
            let sg = (m * sl).min(tan28);
            let fit = (0..p.path.len()).find(|&i| p.lvl - sg * p.d[i] * 1000.0 <= law(i));
            match fit {
                Some(i) => {
                    t1[mi].0 += 1;
                    cols.push(format!("m={m}: {:.2}° fits in {:.2} km", sg.atan().to_degrees(), p.d[i]));
                }
                None => {
                    t1[mi].1 += 1;
                    let short = (p.lvl - sg * x * 1000.0 - p.base).max(0.0);
                    cols.push(format!("m={m}: {:.2}° NO ROOM (needs {:.1} km) → fall {short:.0} m", sg.atan().to_degrees(), dtot / sg / 1000.0));
                }
            }
            if m == 10.0 {
                let phi = phi_of(p.id);
                let (dg, short) = match fit {
                    Some(i) => (p.lvl - law(i), 0.0),
                    None => (sg * x * 1000.0, (p.lvl - sg * x * 1000.0 - p.base).max(0.0)),
                };
                let head = phi * dg;
                let hf = (1.5 * sg * 100.0).max(10.0);
                // the designed gorge over any 2 path cells
                for i in 0..p.path.len().saturating_sub(2) {
                    if sg * (p.d[i + 2] - p.d[i]) * 1000.0 >= hf {
                        gorge_viol += 1;
                    }
                }
                t4.push((p.id, phi, head, short, hf, head > hf, short > hf));
            }
        }
        eprintln!(
            "      lake {:>7} {:>6.1} km² · L {:>7.1} m · base {:>7.1} m ({}) · D {:>6.1} m over X {:>5.2} km · A {:>8.1} km² · S_loi {:.4} ({:.2}°) · {}",
            p.id,
            p.km2,
            p.lvl,
            p.base,
            p.kind,
            dtot,
            x,
            a,
            sl,
            sl.atan().to_degrees(),
            cols.join(" · ")
        );
    }
    eprintln!(
        "   T1 · descents: {} · fit / no room: m = 3 {} / {} · m = 10 **{} / {}** · m = 30 {} / {}",
        ps.len(),
        t1[0].0,
        t1[0].1,
        t1[1].0,
        t1[1].1,
        t1[2].0,
        t1[2].1
    );
    let mut heads: Vec<f32> = t4.iter().filter(|t| t.2 > 0.0).map(|t| t.2).collect();
    let mut shorts: Vec<f32> = t4.iter().filter(|t| t.3 > 0.0).map(|t| t.3).collect();
    heads.sort_by(f32::total_cmp);
    shorts.sort_by(f32::total_cmp);
    let med = |v: &[f32]| if v.is_empty() { f32::NAN } else { v[(v.len() - 1) / 2] };
    eprintln!(
        "\n   T4 · m = 10, φ (splitmix64 of seed, lake): head falls > 0 **{}** (tagged > H_f {}) median {:.1} m, max {:.1} m · shortage falls **{}** (tagged {}) median {:.1} m, max {:.1} m · gorge cells over H_f (2-cell drops) **{gorge_viol}**",
        heads.len(),
        t4.iter().filter(|t| t.5).count(),
        med(&heads),
        heads.last().copied().unwrap_or(f32::NAN),
        shorts.len(),
        t4.iter().filter(|t| t.6).count(),
        med(&shorts),
        shorts.last().copied().unwrap_or(f32::NAN)
    );
    for t in &t4 {
        eprintln!("      lake {:>7} · φ {:.2} · head fall {:>6.1} m{} · shortage fall {:>6.1} m{} · H_f {:.1} m", t.0, t.1, t.2, if t.5 { " TAGGED" } else { "" }, t.3, if t.6 { " TAGGED" } else { "" }, t.4);
    }
    // ── T3 — the lakes by r
    let (lake_in, dep, zbf, inland) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        let wc = water_class(&s1, SEA);
        let inland: Vec<bool> = (0..n).map(|k| s1.data[k] <= SEA && wc[k] == WATER_CLASS_INLAND).collect();
        (d.lake_map, dep, metres(&bf), inland)
    };
    let zs1 = metres(&s1);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    // bodies: label every cell (input lakes by id; depressions and inland basins by component)
    let mut body = vec![0u32; n]; // > 0: an input lake id; ≥ 10 000 000: a component (on the breached field)
    for k in 0..n {
        if lake_in[k] != 0 {
            body[k] = lake_in[k];
        }
    }
    let mut next = 10_000_000u32;
    for s in 0..n {
        if body[s] != 0 || !(dep[s] || inland[s]) {
            continue;
        }
        next += 1;
        body[s] = next;
        let mut q = VecDeque::from([s]);
        while let Some(c) = q.pop_front() {
            for k in 0..8 {
                let m = nbt(c, k);
                if body[m] == 0 && (dep[m] || inland[m]) {
                    body[m] = next;
                    q.push_back(m);
                }
            }
        }
    }
    let mut cells_of: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..n {
        if body[k] != 0 {
            cells_of.entry(body[k]).or_default().push(k);
        }
    }
    let rs = [0f32, 0.5, 1.0, 1.5, 2.0];
    let mut tot = [0f64; 5];
    let mut cnt = [0usize; 5];
    let mut tot_d8 = [0f64; 5];
    eprintln!("\n   T3 · the lakes by r (L(r) between L_in, L_bed = min(L_in, today), L_floor; area grown from the body's lowest cell, cells < L(r)):");
    for p in &ps {
        let mut ov: HashMap<u32, usize> = HashMap::new();
        for k in 0..n {
            if dr.lake_map[k] == p.id && body[k] != 0 {
                *ov.entry(body[k]).or_insert(0) += 1;
            }
        }
        let Some((b, _)) = ov.into_iter().max_by_key(|e| e.1) else {
            eprintln!("      lake {} · no input body", p.id);
            continue;
        };
        let cells = &cells_of[&b];
        let on_bf = b >= 10_000_000;
        let zf = |c: usize| if on_bf { zbf[c] } else { zs1[c] };
        let mut sorted = cells.clone();
        sorted.sort_unstable();
        let mut ring = Vec::new();
        for &c in &sorted {
            for k in 0..8 {
                let m = nbt(c, k);
                if sorted.binary_search(&m).is_err() {
                    ring.push(m);
                }
            }
        }
        let lin = ring.iter().map(|&c| zf(c)).fold(f32::INFINITY, f32::min);
        let lbed = lin.min(p.lvl);
        let lowc = *cells.iter().min_by(|&&a, &&c| zf(a).total_cmp(&zf(c))).unwrap();
        let lfl = zf(lowc);
        let area_at = |lv: f32| -> f64 {
            if lv <= lfl {
                return 0.0;
            }
            let mut seen: HashMap<usize, ()> = HashMap::new();
            let mut q = VecDeque::from([lowc]);
            seen.insert(lowc, ());
            while let Some(c) = q.pop_front() {
                if seen.len() > 2_000_000 {
                    break;
                }
                for k in 0..8 {
                    let m = nbt(c, k);
                    if !seen.contains_key(&m) && zf(m) < lv {
                        seen.insert(m, ());
                        q.push_back(m);
                    }
                }
            }
            seen.len() as f64 * cell_km2 as f64
        };
        let mut row = Vec::new();
        for (ri, &r) in rs.iter().enumerate() {
            let lv = if r <= 1.0 { lin - r * (lin - lbed) } else { lbed - (r - 1.0) * (lbed - lfl) };
            let a = area_at(lv);
            tot[ri] += a;
            if p.id < 1_000_000 {
                tot_d8[ri] += a;
            }
            if a >= 1.0 {
                cnt[ri] += 1;
            }
            row.push(format!("r {r}: {lv:.1} m {a:.1} km²"));
        }
        eprintln!(
            "      lake {:>7} (body {}{}) · L_in {lin:.1} · L_bed {lbed:.1} · L_floor {lfl:.1} · today {:.1} km² · {}",
            p.id,
            if on_bf { "component " } else { "input lake " },
            if on_bf { b - 10_000_000 } else { b },
            p.km2,
            row.join(" · ")
        );
    }
    eprintln!(
        "   T3 · lakes ≥ 1 km² / total area by r: {} · ratio r=0 / r=1 **{:.2}**",
        rs.iter().enumerate().map(|(i, r)| format!("r {r}: {} lakes, {:.1} km²", cnt[i], tot[i])).collect::<Vec<_>>().join(" · "),
        tot[0] / tot[2].max(1e-9)
    );
    eprintln!(
        "   T3 · the 14 D8 lakes only: total area by r: {} · ratio r=0 / r=1 **{:.2}**",
        rs.iter().enumerate().map(|(i, r)| format!("r {r}: {:.1} km²", tot_d8[i])).collect::<Vec<_>>().join(" · "),
        tot_d8[0] / tot_d8[2].max(1e-9)
    );
    eprintln!("   T3 · the selector age each r would propose (DECISION options): (a) ×0.7 → r 0, ×1 → r 1, ×1.4 → r 2 · (b) ×0.7 → r 0.5, ×1 → r 1, ×1.4 → r 1.5");
    eprintln!("\n==========  end Finding 138-T . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 139 — the columns F138 did not produce (C1 T3's calibration, C2 the minimal steepness m_fit, C3 the
/// rule "steepen", C4 T4's variant), the retreat per lake (R) and the below-sea lakes' merged spill (C5). F138's
/// instruments with their amendments; the témoin, ON extended. Declared in `f139_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f139_t --nocapture
#[test]
#[ignore]
fn f139_t() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{HashMap, VecDeque};
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let k_law = F121_AGE_K;
    let tan28 = 28f32.to_radians().tan();
    eprintln!("\n==========  Finding 139 . the columns F138 did not produce, the retreat per lake, the merged spill  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let ext = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let kn = Knobs { valley: Some(ext), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let (w, h) = (wd.heightmap.width, wd.heightmap.height);
    let n = w * h;
    let cn = metres(&v.conditioned);
    let dr = &v.drainage;
    let acc = &dr.flow.accumulation.data;
    let level: HashMap<u32, f32> = dr.lakes.iter().map(|l| (l.base.id, l.level_m)).collect();
    let mut final_cells: HashMap<u32, usize> = HashMap::new();
    for &id in dr.lake_map.iter() {
        if id != 0 {
            *final_cells.entry(id).or_insert(0) += 1;
        }
    }
    let recv = |c: usize| -> Option<usize> {
        let k = dr.flow.direction[c];
        if k == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
    // ── the paths (F138)
    struct P {
        id: u32,
        km2: f32,
        lvl: f32,
        path: Vec<usize>,
        d: Vec<f32>,
        base: f32,
        d8: bool,
    }
    let extend = |path: &mut Vec<usize>, d: &mut Vec<f32>, me: u32| -> f32 {
        let mut c = *path.last().unwrap();
        loop {
            if cn[c] <= 0.0 {
                return 0.0;
            }
            if *d.last().unwrap() >= 30.0 {
                return cn[c];
            }
            let Some(r) = recv(c) else { return cn[c] };
            d.push(d.last().unwrap() + step_len(c, r));
            path.push(r);
            if dr.lake_map[r] != 0 && dr.lake_map[r] != me {
                return level.get(&dr.lake_map[r]).copied().unwrap_or(cn[r]);
            }
            c = r;
        }
    };
    let segs = &dr.rivers.segments;
    let mut ps: Vec<P> = Vec::new();
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        if let Some(col) = recv(o) {
            let (mut path, mut d) = (vec![col], vec![0f32]);
            let base = extend(&mut path, &mut d, l.base.id);
            ps.push(P { id: l.base.id, km2: l.area_km2, lvl: l.level_m, path, d, base, d8: true });
        } else if l.base.id >= 1_000_000 && l.base.id < 2_000_000 {
            let sp: Vec<usize> = (0..segs.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway && dr.segment_source_lake.get(i).copied().flatten() == Some(l.base.id))
                .collect();
            if sp.is_empty() {
                continue;
            }
            let mut path: Vec<usize> = Vec::new();
            let mut si = Some(sp[0]);
            let mut guard = 0;
            while let Some(i) = si {
                for &(x, y) in &segs[i].points {
                    let p = y as usize * w + x as usize;
                    if path.last() != Some(&p) {
                        path.push(p);
                    }
                }
                si = segs[i].downstream;
                guard += 1;
                if guard > 10_000 {
                    break;
                }
            }
            let mut d = vec![0f32; path.len()];
            for i in 1..path.len() {
                d[i] = d[i - 1] + step_len(path[i - 1], path[i]);
            }
            let hit = (1..path.len()).find(|&i| (dr.lake_map[path[i]] != 0 && dr.lake_map[path[i]] != l.base.id) || cn[path[i]] <= 0.0);
            let base = match hit {
                Some(i) => {
                    path.truncate(i + 1);
                    d.truncate(i + 1);
                    let q = path[i];
                    if cn[q] <= 0.0 { 0.0 } else { level.get(&dr.lake_map[q]).copied().unwrap_or(cn[q]) }
                }
                None => extend(&mut path, &mut d, l.base.id),
            };
            ps.push(P { id: l.base.id, km2: l.area_km2, lvl: l.level_m, path, d, base, d8: false });
        }
    }
    let splitmix = |id: u32| -> f32 {
        let mut z = PSEED ^ (id as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        let u = (z >> 11) as f64 / (1u64 << 53) as f64;
        if u < 1.0 / 3.0 { 0.0 } else { (0.1 + 0.4 * (u - 1.0 / 3.0) / (2.0 / 3.0)) as f32 }
    };
    // ── per descent: A, S_loi, the law, C2, C3, C4
    let mut area_of: HashMap<u32, f32> = HashMap::new();
    let mut mfits: Vec<f32> = Vec::new();
    let (mut fixed_falls, mut steep_falls) = (Vec::new(), Vec::new());
    let mut c4 = [(0usize, 0usize); 2];
    let mut c4_falls: [Vec<f32>; 2] = [Vec::new(), Vec::new()];
    eprintln!("\n   C2 / C3 / C4 · per descent (L = today's level, B = the next base, law = B + k·χ along the path):");
    for p in &ps {
        let inflow: f32 = (0..n).filter(|&c| dr.lake_map[c] != p.id && recv(c).is_some_and(|r| dr.lake_map[r] == p.id)).map(|c| acc[c] * cell_km2).sum::<f32>() + p.km2;
        let a = p.path.iter().zip(p.d.iter()).filter(|e| *e.1 <= 1.0).map(|e| acc[*e.0] * cell_km2).fold(0f32, f32::max).max(inflow).max(0.1);
        area_of.insert(p.id, a);
        let sl = k_law * a.powf(-0.5);
        let x = *p.d.last().unwrap();
        let mut chi = vec![0f32; p.path.len()];
        for i in (0..p.path.len().saturating_sub(1)).rev() {
            let ai = (acc[p.path[i]] * cell_km2).max(a).max(0.1);
            chi[i] = chi[i + 1] + (1.0 / ai).sqrt() * (p.d[i + 1] - p.d[i]) * 1000.0;
        }
        let law = |i: usize| -> f32 { p.base + k_law * chi[i] };
        // the required slope from a start level `top`
        let s_req = |top: f32| -> f32 {
            if top <= law(0) {
                return 0.0;
            }
            (1..p.path.len()).map(|i| (top - law(i)) / (p.d[i] * 1000.0)).fold(f32::INFINITY, f32::min).max(0.0)
        };
        let fits_at = |top: f32, s: f32| -> Option<usize> { (0..p.path.len()).find(|&i| top - s * p.d[i] * 1000.0 <= law(i)) };
        let sr = s_req(p.lvl);
        let mfit = sr / sl;
        mfits.push(mfit);
        let beyond = if sr > tan28 { (p.lvl - tan28 * x * 1000.0 - p.base).max(0.0) } else { 0.0 };
        // C3: fixed slope (m = 10) against steepen
        let s10 = (10.0 * sl).min(tan28);
        let fixed_fall = if fits_at(p.lvl, s10).is_none() { (p.lvl - s10 * x * 1000.0 - p.base).max(0.0) } else { 0.0 };
        if fixed_fall > 0.0 {
            fixed_falls.push(fixed_fall);
        }
        let s_st = s10.max(sr).min(tan28);
        let steep_fall = if fits_at(p.lvl, s_st).is_none() { (p.lvl - s_st * x * 1000.0 - p.base).max(0.0) } else { 0.0 };
        if steep_fall > 0.0 {
            steep_falls.push(steep_fall);
        }
        // C4: the head fall first, then the gorge from L'
        let phi = splitmix(p.id);
        let dg0 = match fits_at(p.lvl, s10) {
            Some(i) => p.lvl - law(i),
            None => s10 * x * 1000.0,
        };
        let top = p.lvl - phi * dg0;
        let mut c4row = Vec::new();
        for (ri, s) in [s10, s10.max(s_req(top)).min(tan28)].into_iter().enumerate() {
            match fits_at(top, s) {
                Some(i) => {
                    c4[ri].0 += 1;
                    c4row.push(format!("{} fits in {:.2} km", ["fixed", "steepen"][ri], p.d[i]));
                }
                None => {
                    c4[ri].1 += 1;
                    let f = (top - s * x * 1000.0 - p.base).max(0.0);
                    c4_falls[ri].push(f);
                    c4row.push(format!("{} NO ROOM → fall {f:.0} m", ["fixed", "steepen"][ri]));
                }
            }
        }
        eprintln!(
            "      lake {:>7} · L {:>7.1} → B {:>7.1} m over {:>5.2} km · A {:>7.1} km² · S_loi {:.4} · S_req {:.4} ({:.2}°) → **m_fit {mfit:>6.1}**{} · m = 10 fixed: fall {fixed_fall:.0} m · steepen ({:.2}°): fall {steep_fall:.0} m · C4 φ {phi:.2}, head {:.0} m: {}",
            p.id,
            p.lvl,
            p.base,
            x,
            a,
            sl,
            sr,
            sr.atan().to_degrees(),
            if beyond > 0.0 { format!(" BEYOND THE CAP, remainder {beyond:.0} m") } else { String::new() },
            s_st.atan().to_degrees(),
            phi * dg0,
            c4row.join(" · ")
        );
    }
    let mut ms = mfits.clone();
    ms.sort_by(f32::total_cmp);
    let pct = |v: &[f32], q: f64| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f64 * q) as usize] };
    eprintln!(
        "   C2 · m_fit over {} descents: min **{:.2}** · p50 **{:.2}** · max **{:.1}** · beyond the 28° cap: {}",
        ms.len(),
        ms[0],
        pct(&ms, 0.5),
        ms[ms.len() - 1],
        mfits.iter().zip(ps.iter()).filter(|(m, p)| **m * k_law * area_of[&p.id].powf(-0.5) > tan28).count()
    );
    let desc = |v: &[f32]| -> String {
        let mut s = v.to_vec();
        s.sort_by(f32::total_cmp);
        if s.is_empty() { "none".to_string() } else { format!("{} falls, median {:.0} m, max {:.0} m", s.len(), pct(&s, 0.5), s[s.len() - 1]) }
    };
    eprintln!("   C3 · m = 10, shortage falls · fixed slope (F138): {} · STEEPEN: **{}**", desc(&fixed_falls), desc(&steep_falls));
    eprintln!(
        "   C4 · m = 10, head fall first · fixed slope: fit {} / no room {} ({}) · steepen: fit {} / no room {} ({})",
        c4[0].0,
        c4[0].1,
        desc(&c4_falls[0]),
        c4[1].0,
        c4[1].1,
        desc(&c4_falls[1])
    );
    // ── bodies on the input (F138-T3)
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let zs1 = metres(&s1);
    let (lake_in, dep, zbf, inland, spill) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        let wc = water_class(&s1, SEA);
        let inland: Vec<bool> = (0..n).map(|k| s1.data[k] <= SEA && wc[k] == WATER_CLASS_INLAND).collect();
        (d.lake_map, dep, metres(&bf), inland, metres(&GridF32 { width: w, height: h, data: spill }))
    };
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let mut body = vec![0u32; n];
    for k in 0..n {
        if lake_in[k] != 0 {
            body[k] = lake_in[k];
        }
    }
    let mut next = 10_000_000u32;
    for s in 0..n {
        if body[s] != 0 || !(dep[s] || inland[s]) {
            continue;
        }
        next += 1;
        body[s] = next;
        let mut q = VecDeque::from([s]);
        while let Some(c) = q.pop_front() {
            for k in 0..8 {
                let m = nbt(c, k);
                if body[m] == 0 && (dep[m] || inland[m]) {
                    body[m] = next;
                    q.push_back(m);
                }
            }
        }
    }
    let mut cells_of: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..n {
        if body[k] != 0 {
            cells_of.entry(body[k]).or_default().push(k);
        }
    }
    let area_from = |low: usize, lv: f32, on_bf: bool| -> f64 {
        let zf = |c: usize| if on_bf { zbf[c] } else { zs1[c] };
        if lv <= zf(low) {
            return 0.0;
        }
        let mut seen: HashMap<usize, ()> = HashMap::new();
        let mut q = VecDeque::from([low]);
        seen.insert(low, ());
        while let Some(c) = q.pop_front() {
            if seen.len() > 2_000_000 {
                break;
            }
            for k in 0..8 {
                let m = nbt(c, k);
                if !seen.contains_key(&m) && zf(m) < lv {
                    seen.insert(m, ());
                    q.push_back(m);
                }
            }
        }
        seen.len() as f64 * cell_km2 as f64
    };
    struct B {
        id: u32,
        lin: f32,
        lbed: f32,
        lfl: f32,
        low: usize,
        on_bf: bool,
    }
    let mut bs: Vec<B> = Vec::new();
    eprintln!("\n   C1 · T3's area at r = 1 against the final area (run_hd lake map):");
    let mut off20 = 0usize;
    for p in &ps {
        let mut ov: HashMap<u32, usize> = HashMap::new();
        for k in 0..n {
            if dr.lake_map[k] == p.id && body[k] != 0 {
                *ov.entry(body[k]).or_insert(0) += 1;
            }
        }
        let Some((b, _)) = ov.into_iter().max_by_key(|e| e.1) else {
            eprintln!("      lake {:>7} · no input body", p.id);
            continue;
        };
        let cells = &cells_of[&b];
        let on_bf = b >= 10_000_000;
        let zf = |c: usize| if on_bf { zbf[c] } else { zs1[c] };
        let mut sorted = cells.clone();
        sorted.sort_unstable();
        let mut lin = f32::INFINITY;
        for &c in &sorted {
            for k in 0..8 {
                let m = nbt(c, k);
                if sorted.binary_search(&m).is_err() {
                    lin = lin.min(zf(m));
                }
            }
        }
        let lbed = lin.min(p.lvl);
        let low = *cells.iter().min_by(|&&a, &&c| zf(a).total_cmp(&zf(c))).unwrap();
        let a1 = area_from(low, lbed, on_bf);
        let fin = final_cells.get(&p.id).copied().unwrap_or(0) as f64 * cell_km2 as f64;
        let ratio = a1 / fin.max(1e-9);
        if (ratio - 1.0).abs() > 0.2 {
            off20 += 1;
        }
        eprintln!("      lake {:>7} ({}) · T3 at r = 1 {a1:>7.1} km² · final (run_hd mask) {fin:>7.1} km² · ratio **{ratio:.3}**{}", p.id, if p.d8 { "D8" } else { "spillway" }, if (ratio - 1.0).abs() > 0.2 { " · > 20 % OFF" } else { "" });
        bs.push(B { id: p.id, lin, lbed, lfl: zf(low), low, on_bf });
    }
    eprintln!("   C1 · lakes off by > 20 %: **{off20} of {}**{}", bs.len(), if off20 * 3 > bs.len() { " — MORE THAN A THIRD (said first)" } else { "" });
    // ── R — retreat per lake (14 D8 lakes)
    let d8: Vec<&B> = bs.iter().filter(|b| ps.iter().any(|p| p.id == b.id && p.d8)).collect();
    let mut a_sorted: Vec<f32> = d8.iter().map(|b| area_of[&b.id]).collect();
    a_sorted.sort_by(f32::total_cmp);
    let a_ref = a_sorted[(a_sorted.len() - 1) / 2];
    eprintln!("\n   R · retreat per lake, the {} D8 lakes · A_ref (median outlet area) {a_ref:.1} km² · map (a): r 0 → ×0.7, 1 → ×1, 2 → ×1.4 (piecewise linear)", d8.len());
    let age_of = |r: f32| if r <= 1.0 { 0.7 + 0.3 * r } else { 1.0 + 0.4 * (r - 1.0) };
    let lvl_at = |b: &B, r: f32| if r <= 1.0 { b.lin - r * (b.lin - b.lbed) } else { b.lbed - (r - 1.0) * (b.lbed - b.lfl) };
    for p_exp in [0f32, 0.25, 0.5] {
        let mut row = Vec::new();
        for rw in [0f32, 0.5, 1.0, 1.5, 2.0] {
            let (mut cnt, mut tot) = (0usize, 0f64);
            for b in &d8 {
                let rl = (rw * (area_of[&b.id] / a_ref).powf(p_exp)).min(2.0);
                let a = area_from(b.low, lvl_at(b, rl), b.on_bf);
                tot += a;
                if a >= 1.0 {
                    cnt += 1;
                }
            }
            row.push(format!("r_world {rw}: {cnt} lakes, {tot:.1} km²"));
        }
        let mut ages: Vec<f32> = Vec::new();
        let mut never = 0usize;
        let mut per = Vec::new();
        for b in &d8 {
            let re = 2.0 / (area_of[&b.id] / a_ref).powf(p_exp);
            if re <= 2.0 + 1e-6 {
                ages.push(age_of(re));
                per.push(format!("{}: ×{:.2}", b.id, age_of(re)));
            } else {
                never += 1;
                per.push(format!("{}: never (r_world {re:.2})", b.id));
            }
        }
        ages.sort_by(f32::total_cmp);
        let spread = if ages.is_empty() { 0.0 } else { (ages[ages.len() - 1] - ages[0]) / 0.7 };
        eprintln!("      p = {p_exp}: {}", row.join(" · "));
        eprintln!(
            "         emptying age: {} lakes within ×0.7–×1.4 (from ×{:.2} to ×{:.2}, spread **{:.0} %** of the selector's range) · never within the range: {never} · {}",
            ages.len(),
            ages.first().copied().unwrap_or(f32::NAN),
            ages.last().copied().unwrap_or(f32::NAN),
            100.0 * spread,
            per.join(" · ")
        );
    }
    eprintln!("      p's provenance: DECISION; the link to the celerity (F115, c = K·A^m) is ANCHORED in form, not in value");
    // ── C5 — the below-sea lakes' merged spill
    eprintln!("\n   C5 · the spillway lakes' merged spill (median of ocean_flood's spill on the breached S1 over the final footprint):");
    let mut ok1 = 0usize;
    let mut nsp = 0usize;
    for p in ps.iter().filter(|p| !p.d8) {
        nsp += 1;
        let mut v: Vec<f32> = (0..n).filter(|&k| dr.lake_map[k] == p.id).map(|k| spill[k]).collect();
        v.sort_by(f32::total_cmp);
        let lm = v[(v.len() - 1) / 2];
        let ok = (lm - p.lvl).abs() <= 1.0;
        if ok {
            ok1 += 1;
        }
        let bt = bs.iter().find(|b| b.id == p.id);
        eprintln!(
            "      lake {:>7} · today {:>7.1} m · merged spill on the input **{lm:>7.1} m** (Δ {:+.1} m{}) · the F138 body's ring {} · L(1) with the merged spill {:.1} m",
            p.id,
            p.lvl,
            lm - p.lvl,
            if ok { ", within 1 m" } else { "" },
            bt.map_or("—".to_string(), |b| format!("{:.1} m", b.lin)),
            lm.min(p.lvl)
        );
    }
    eprintln!("   C5 · within 1 m: **{ok1} of {nsp}**");
    eprintln!("\n==========  end Finding 139 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 140-G — the gorge's retreat, built gated, measured against its gates on the témoin: GORGE at ×0.7 /
/// ×0.85 / ×1 / ×1.2 / ×1.4 (p = 0.5) and at ×1 (p = 0, 0.25); OFF and ON at ×1 for the controls, rule 14 and
/// rule 18. Declared in `f140_declared.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f140_g --nocapture
#[test]
#[ignore]
fn f140_g() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::HashMap;
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::valley_construction::{GorgeFallKind, GorgeRetreat, LakeBase, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let tan28 = 28f32.to_radians().tan();
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 140-G . the gorge's retreat against its gates (témoin)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let (lm_in, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &vd, &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let sk_off = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let links = trunk_links(&sk_off);
    let lb: Vec<bool> = links.iter().map(|l| lm_in[l.0] != 0 || dep[l.0] || lm_in[l.1] != 0 || dep[l.1]).collect();
    drop(sk_off);
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let map_a = |age: f32| if age <= 1.0 { (age - 0.7) / 0.3 } else { 1.0 + (age - 1.0) / 0.4 };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    // the body set of the reference GORGE world (×1, p = 0.5), for the ON / OFF negative controls
    struct RefB {
        cells: Vec<usize>,
        ring: Vec<usize>,
        l_in: f32,
        exit: usize,
    }
    let mut ref_bodies: Vec<RefB> = Vec::new();
    let mut area_series: Vec<(f32, f64, f64)> = Vec::new(); // (age, all final lakes, input-lake bodies' lakes)
    let mut r18: Vec<String> = Vec::new();
    let mut t_on_con = 0f64;
    // the extras (canyons, coast, θ, removal) for one world
    let extras = |label: &str, vc: ValleyConstruction, wd_h: &GridF32, cond: &GridF32, carved: &[bool], bodies_cells: &[bool]| -> String {
        let g = wd_h;
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(g, &vd, &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(g, cond, &pre, DELIVERED_P50_M, &ss, &vd, cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = 0;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons += 1;
            }
        }
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        let co = coast(g, &ss, &dwall, 2.0 / CELL_KM);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let eo = metres(g);
        let zk: Vec<f32> = links.iter().map(|l| eo[l.0]).collect();
        let _ = zk;
        let rem = |sel: &dyn Fn(usize) -> bool| -> f64 {
            (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((zs1[k] - eo[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum()
        };
        let (rt, ri) = (rem(&|_| true), rem(&|k| bodies_cells[k]));
        let _ = carved;
        format!("canyons **{canyons}** · spurs near a coastal wall {rn:.4} /km · removal (S1 − world, land) **{rt:.1} km³** (inside the lakes' bodies {ri:.1}, elsewhere {:.1}) [{label}]", rt - ri)
    };
    // θ on the carved links with the construction's mask (F135-M2)
    let theta_carved = |zb: &[f32], carved: &[bool]| -> String {
        let zk: Vec<f32> = links.iter().map(|l| zb[l.0]).collect();
        let zr: Vec<f32> = links.iter().map(|l| zb[l.1]).collect();
        let cl: Vec<bool> = links.iter().map(|l| carved[l.0] && carved[l.1]).collect();
        let t = theta_links(&links, &zk, &zr, &|i| cl[i] && !lb[i]);
        format!("θ on the carved links, construction mask **{:.3}** [{:.3}, {:.3}] ({})", t.0, t.1, t.2, t.3)
    };
    // ── ON at ×1 first (its construction time, the controls come after GORGE ×1 gives the bodies)
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let (zb_on, carved_on) = {
        let t = Instant::now();
        let sk = skeleton(&s1, &on, &ss, DOMAIN_KM);
        let (b, mk) = carve(&s1, &sk, &on, &ss);
        t_on_con = t.elapsed().as_secs_f64();
        (metres(&b), mk.carved)
    };
    let runs: Vec<(f32, f32)> = vec![(1.0, 0.5), (0.7, 0.5), (0.85, 0.5), (1.2, 0.5), (1.4, 0.5), (1.0, 0.0), (1.0, 0.25)];
    for (age, p) in runs {
        let r = map_a(age);
        let vc = ValleyConstruction {
            lake_base: Some(LakeBase::InputLakesAndBasins),
            gorge_retreat: Some(GorgeRetreat::v3(r, p)),
            ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1))
        };
        let label = format!("GORGE ×{age} (r_world {r}, p {p})");
        eprintln!("\n────────── {label} ──────────");
        let t = Instant::now();
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (built, mk) = carve(&s1, &sk, &vc, &ss);
        let t_con = t.elapsed().as_secs_f64();
        let zb = metres(&built);
        drop(built);
        let bo = sk.gorge_body_of.clone().expect("the gate is on");
        let bodies = sk.gorge_bodies.clone();
        let rf = sk.rim_floor_m.clone().expect("the gate is on");
        let dir_sk = sk.direction.clone();
        let falls_sk = sk.gorge_falls.clone();
        drop(sk);
        let wd = build_world(kn(vc), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0);
        let cn = metres(&v.conditioned);
        let eo = metres(&wd.heightmap);
        let lm = &v.drainage.lake_map;
        let lakes = &v.drainage.lakes;
        let lvl_of: HashMap<u32, (f32, f32)> = lakes.iter().map(|l| (l.base.id, (l.level_m, l.area_km2))).collect();
        eprintln!(
            "   construction {t_con:.1} s (ON at ×1: {t_on_con:.1} s) · {} bodies · {} tagged falls in the world ({} in the skeleton) · {} final lakes",
            bodies.len(),
            wd.gorge_falls.len(),
            falls_sk.len(),
            lakes.len()
        );
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); bodies.len()];
        for k in 0..n {
            if bo[k] != u32::MAX {
                cells[bo[k] as usize].push(k);
            }
        }
        let ring_of = |bi: usize| -> Vec<usize> {
            let mut r: Vec<usize> = Vec::new();
            for &c in &cells[bi] {
                for k in 0..8 {
                    let m = nbt(c, k);
                    if bo[m] == u32::MAX && s1.data[m] > SEA {
                        r.push(m);
                    }
                }
            }
            r.sort_unstable();
            r.dedup();
            r
        };
        let is_lake_body: Vec<bool> = cells.iter().map(|cs| 2 * cs.iter().filter(|&&c| lm_in[c] != 0).count() >= cs.len()).collect();
        // upstream of each body on the skeleton (for G-drained)
        let mut up_body = vec![u32::MAX - 1; n]; // MAX-1 unknown, MAX none
        let mut path = Vec::new();
        for s in 0..n {
            if up_body[s] != u32::MAX - 1 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if up_body[c] != u32::MAX - 1 {
                    v0 = up_body[c];
                    break;
                }
                if bo[c] != u32::MAX && c != s {
                    v0 = bo[c];
                    break;
                }
                path.push(c);
                let d = dir_sk[c];
                if d == DIR_NONE || s1.data[c] <= SEA || path.len() > 50_000 {
                    v0 = u32::MAX;
                    break;
                }
                c = nbt(c, d as usize);
            }
            for &q in &path {
                up_body[q] = v0;
            }
        }
        let head_at: HashMap<u32, f32> = falls_sk.iter().filter(|f| f.kind == GorgeFallKind::WaterfallHead).map(|f| (f.exit_y * w as u32 + f.exit_x, f.height_m)).collect();
        let (mut nl, mut ok_lv, mut drained, mut drained_ok) = (0usize, 0usize, 0usize, 0usize);
        let mut bad_lv: Vec<String> = Vec::new();
        let (mut rim_n, mut rim_bad_b, mut rim_bad_w) = (0usize, 0usize, 0usize);
        let (mut ring_bound, mut ring_viol) = (0usize, 0usize);
        let mut ring_worst = 0f32;
        let mut ring_where: Vec<String> = Vec::new();
        let (mut s100_n, mut s100_bad) = (0usize, 0usize);
        let mut s100_list: Vec<f32> = Vec::new();
        let (mut drn_up, mut drn_bad) = (0usize, 0usize);
        let mut gorge_viol = 0usize;
        let mut lake_area_bodies: HashMap<u32, f32> = HashMap::new();
        for (bi, b) in bodies.iter().enumerate() {
            if b.col == u32::MAX {
                continue;
            }
            nl += 1;
            // G-levels
            let mut ov: HashMap<u32, usize> = HashMap::new();
            for &c in &cells[bi] {
                if lm[c] != 0 {
                    *ov.entry(lm[c]).or_insert(0) += 1;
                }
            }
            let best = ov.iter().max_by_key(|e| *e.1).map(|(&id, &c)| (id, c));
            if b.r_lake >= 1.999 {
                drained += 1;
                if best.is_none_or(|(_, c)| 2 * c < cells[bi].len()) {
                    drained_ok += 1;
                }
            } else if let Some((id, _)) = best {
                let (lv, ar) = lvl_of[&id];
                if (lv - b.level).abs() <= 1.0 {
                    ok_lv += 1;
                } else {
                    bad_lv.push(format!("body@({},{}) {:.1} km² r {:.2} L(r) {:.1} final {lv:.1} (A {:.0} km²)", b.low % w as u32, b.low / w as u32, cells[bi].len() as f32 * cell_km2, b.r_lake, b.level, b.a_out_km2));
                }
                if is_lake_body[bi] {
                    lake_area_bodies.insert(id, ar);
                }
            } else {
                bad_lv.push(format!("body@({},{}) {:.1} km² r {:.2} L(r) {:.1}: NO final lake", b.low % w as u32, b.low / w as u32, cells[bi].len() as f32 * cell_km2, b.r_lake, b.level));
            }
            // G-rim
            let ring = ring_of(bi);
            if b.r_lake <= 1.0 {
                rim_n += 1;
                let mb = ring.iter().map(|&c| zb[c]).fold(f32::INFINITY, f32::min);
                let mw = ring.iter().map(|&c| eo[c]).fold(f32::INFINITY, f32::min);
                if b.level - mb > 10.0 {
                    rim_bad_b += 1;
                }
                if b.level - mw > 10.0 {
                    rim_bad_w += 1;
                }
            }
            // G-ring
            for &c in &ring {
                if rf[c].is_finite() && (zb[c] - rf[c]).abs() < 0.05 && zs1[c] > rf[c] + 0.05 {
                    ring_bound += 1;
                    if c == b.col as usize {
                        continue;
                    }
                    for k in 0..8 {
                        let m = nbt(c, k);
                        if bo[m] != u32::MAX || rf[m].is_finite() {
                            continue;
                        }
                        let dist = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                        let sl = (zb[c] - zb[m]) / dist;
                        if sl > tan28 {
                            ring_viol += 1;
                            if sl > ring_worst {
                                ring_worst = sl;
                            }
                            if ring_where.len() < 6 {
                                ring_where.push(format!("({},{}) {:.0}°", c % w, c / w, sl.atan().to_degrees()));
                            }
                            break;
                        }
                    }
                }
            }
            // G-slope100 on the final world, from the final lake's true col
            if let Some((id, _)) = best
                && b.r_lake < 1.999
                && let Some(l) = lakes.iter().find(|l| l.base.id == id)
            {
                let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
                let rcv = |c: usize| -> Option<usize> {
                    let d = v.drainage.flow.direction[c];
                    if d == DIR_NONE {
                        return None;
                    }
                    let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                    if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
                };
                if let Some(col) = rcv(o) {
                    let (mut pth, mut dd) = (vec![col], vec![0f32]);
                    let mut c = col;
                    while *dd.last().unwrap() < 10.5 && cn[c] > 0.0 {
                        let Some(r) = rcv(c) else { break };
                        if lm[r] != 0 && lm[r] != id {
                            break;
                        }
                        dd.push(dd.last().unwrap() + CELL_KM * if (c % w != r % w) && (c / w != r / w) { std::f32::consts::SQRT_2 } else { 1.0 });
                        pth.push(r);
                        c = r;
                    }
                    let skip = if head_at.contains_key(&b.exit) { 2 } else { 0 };
                    let mut ms = 0f32;
                    for i in skip..pth.len() {
                        if dd[i] > 10.0 {
                            break;
                        }
                        if let Some(j) = (i..pth.len()).find(|&j| dd[j] >= dd[i] + 0.1) {
                            ms = ms.max((cn[pth[i]] - cn[pth[j]]) / ((dd[j] - dd[i]) * 1000.0));
                        }
                    }
                    s100_n += 1;
                    s100_list.push(ms.atan().to_degrees());
                    if ms > tan28 {
                        s100_bad += 1;
                    }
                }
            }
            // G-drained
            if b.r_lake >= 1.999 {
                for k in 0..n {
                    if up_body[k] == bi as u32 {
                        drn_up += 1;
                        if zb[k] < b.l_floor - 1.0 {
                            drn_bad += 1;
                        }
                    }
                }
            }
            // G-tag: the gorge's cells beyond the lip, on the construction
            let hf = (1.5 * b.slope * 100.0).max(10.0);
            let glen = if b.slope > 0.0 { (b.d_g / b.slope / 1000.0).min(20.0) } else { 0.0 };
            let mut gp = vec![b.col as usize];
            let mut sd = 0f32;
            let mut c = b.col as usize;
            while sd < glen {
                let d = dir_sk[c];
                if d == DIR_NONE {
                    break;
                }
                let r2 = nbt(c, d as usize);
                if s1.data[r2] <= SEA || bo[r2] != u32::MAX {
                    break;
                }
                sd += CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                gp.push(r2);
                c = r2;
            }
            for i in 1..gp.len().saturating_sub(2) {
                if zb[gp[i]] - zb[gp[i + 2]] > hf {
                    gorge_viol += 1;
                }
            }
            if age == 1.0 && p == 0.5 {
                ref_bodies.push(RefB { cells: cells[bi].clone(), ring, l_in: b.l_in, exit: b.exit as usize });
            }
        }
        let all_area: f64 = lakes.iter().filter(|l| !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).map(|l| l.area_km2 as f64).sum();
        let lb_area: f64 = lake_area_bodies.values().map(|&a| a as f64).sum();
        if p == 0.5 {
            area_series.push((age, all_area, lb_area));
        }
        eprintln!("   G-levels · bodies with an outflow {nl} · at L(r) ± 1 m **{ok_lv}** of {} not drained · drained (r_lake = 2) {drained}, absent at the end **{drained_ok}**", nl - drained);
        for l in bad_lv.iter().take(12) {
            eprintln!("      off: {l}");
        }
        eprintln!("   G-rim (r_lake ≤ 1) · {rim_n} bodies · cuts > 10 m below L(r) on the construction **{rim_bad_b}** · on the eroded world **{rim_bad_w}**");
        eprintln!("   G-ring · ring cells bound by the clamp **{ring_bound}** · with a slope > 28° to an outside neighbour **{ring_viol}** (worst {:.0}°) {:?}", ring_worst.atan().to_degrees(), ring_where);
        eprintln!(
            "   G-slope100 · {s100_n} outlet paths · max slope over 100 m > 28° **{s100_bad}** · p50 {:.1}° · max {:.1}°",
            q(&mut s100_list.clone(), 0.5),
            q(&mut s100_list, 1.0)
        );
        eprintln!("   G-drained · upstream cells of drained lakes {drn_up} · below L_floor − 1 m on the construction **{drn_bad}**");
        let mut fh: Vec<f32> = wd.gorge_falls.iter().filter(|f| f.kind == GorgeFallKind::WaterfallHead).map(|f| f.height_m).collect();
        let ns = wd.gorge_falls.iter().filter(|f| f.kind == GorgeFallKind::WaterfallShortage).count();
        eprintln!(
            "   G-tag · tagged falls: head **{}** (median {:.1} m, max {:.1} m) · shortage **{ns}** · gorge 2-cell drops over H_f beyond the lip **{gorge_viol}**",
            fh.len(),
            q(&mut fh.clone(), 0.5),
            q(&mut fh, 1.0)
        );
        eprintln!("   G-area · final non-crater lakes {all_area:.1} km² · the input-lake bodies' lakes {lb_area:.1} km²");
        let mut bcells = vec![false; n];
        for k in 0..n {
            bcells[k] = bo[k] != u32::MAX;
        }
        if age == 1.0 && p == 0.5 {
            let ex = extras(&label, vc, &wd.heightmap, &v.conditioned, &mk.carved, &bcells);
            let th = theta_carved(&zb, &mk.carved);
            eprintln!("   ×1 extras · {ex} · {th}");
            r18.push(format!("GORGE ×1 p 0.5: {} lakes {all_area:.1} km² · {ex} · {th}", lakes.len()));
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    // ── the controls on ON and OFF at ×1 (with the reference bodies)
    let mut ref_cells = vec![false; n];
    for b in &ref_bodies {
        for &c in &b.cells {
            ref_cells[c] = true;
        }
    }
    for (label, vc, zb_c) in [("ON ×1", on, Some((&zb_on, &carved_on))), ("OFF ×1", off, None)] {
        eprintln!("\n────────── {label} (controls) ──────────");
        let wd = build_world(kn(vc), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0);
        let cn = metres(&v.conditioned);
        let lm = &v.drainage.lake_map;
        if let Some((zb, carved)) = zb_c {
            let cuts = ref_bodies.iter().filter(|b| b.l_in - b.ring.iter().map(|&c| zb[c]).fold(f32::INFINITY, f32::min) > 10.0).count();
            eprintln!("   G-rim negative control · cuts > 10 m below L_in on the construction **{cuts}** of {} reference bodies", ref_bodies.len());
            // G-slope100 negative control
            let (mut nn, mut bad) = (0usize, 0usize);
            let mut lst = Vec::new();
            for l in v.drainage.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
                let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
                let rcv = |c: usize| -> Option<usize> {
                    let d = v.drainage.flow.direction[c];
                    if d == DIR_NONE {
                        return None;
                    }
                    let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                    if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
                };
                let Some(col) = rcv(o) else { continue };
                let (mut pth, mut dd) = (vec![col], vec![0f32]);
                let mut c = col;
                while *dd.last().unwrap() < 10.5 && cn[c] > 0.0 {
                    let Some(r) = rcv(c) else { break };
                    if lm[r] != 0 && lm[r] != l.base.id {
                        break;
                    }
                    dd.push(dd.last().unwrap() + CELL_KM * if (c % w != r % w) && (c / w != r / w) { std::f32::consts::SQRT_2 } else { 1.0 });
                    pth.push(r);
                    c = r;
                }
                let mut ms = 0f32;
                for i in 0..pth.len() {
                    if dd[i] > 10.0 {
                        break;
                    }
                    if let Some(j) = (i..pth.len()).find(|&j| dd[j] >= dd[i] + 0.1) {
                        ms = ms.max((cn[pth[i]] - cn[pth[j]]) / ((dd[j] - dd[i]) * 1000.0));
                    }
                }
                nn += 1;
                lst.push(ms.atan().to_degrees());
                if ms > tan28 {
                    bad += 1;
                }
            }
            eprintln!("   G-slope100 negative control · {nn} outlet paths · > 28° **{bad}** · max {:.1}°", q(&mut lst, 1.0));
            let ex = extras(label, vc, &wd.heightmap, &v.conditioned, carved, &ref_cells);
            let th = theta_carved(zb, carved);
            eprintln!("   ×1 extras · {ex} · {th}");
            r18.push(format!("{label}: {} lakes {:.1} km² · {ex} · {th}", v.drainage.lakes.len(), v.drainage.lakes.iter().filter(|l| !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).map(|l| l.area_km2).sum::<f32>()));
        } else {
            let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
            let (b, mk) = carve(&s1, &sk, &vc, &ss);
            let zb = metres(&b);
            let ex = extras(label, vc, &wd.heightmap, &v.conditioned, &mk.carved, &ref_cells);
            let th = theta_carved(&zb, &mk.carved);
            eprintln!("   ×1 extras · {ex} · {th}");
            r18.push(format!("{label}: {} lakes {:.1} km² · {ex} · {th}", v.drainage.lakes.len(), v.drainage.lakes.iter().filter(|l| !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).map(|l| l.area_km2).sum::<f32>()));
        }
        let _ = ref_bodies.iter().map(|b| b.exit).count();
    }
    area_series.sort_by(|a, b| a.0.total_cmp(&b.0));
    eprintln!("\n   G-area (p = 0.5) · by age: {}", area_series.iter().map(|a| format!("×{}: all {:.1} km², input-lake bodies {:.1} km²", a.0, a.1, a.2)).collect::<Vec<_>>().join(" · "));
    let mono = area_series.windows(2).all(|p| p[1].1 <= p[0].1 + 1e-6);
    let mono_b = area_series.windows(2).all(|p| p[1].2 <= p[0].2 + 1e-6);
    eprintln!("   G-area · non-increasing: all lakes **{mono}** · input-lake bodies **{mono_b}** (F139-R: 1 308.5 / 1 202.7 / 521.3 / 103.0 / 19.0)");
    eprintln!("\n   RULE 18 · {}", r18.join("\n   RULE 18 · "));
    eprintln!("\n==========  end Finding 140-G . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 140 — the diagnosis of f140_g's failure: GORGE ×1 p 0.5 made the HD assembly fire Finding 38's invariant
/// ("enclosed below-sea components carry no water body", floor cells (4379, 3551) and (2435, 3586)). Which bodies lay a
/// level below the sea, and are the new below-sea land cells theirs?
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f140_diag --nocapture
#[test]
#[ignore]
fn f140_diag() {
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    let ss = SteinSteinParams::default();
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, _h) = (s1.width, s1.height);
    let n = s1.data.len();
    let vc = ValleyConstruction {
        lake_base: Some(LakeBase::InputLakesAndBasins),
        gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)),
        ..ValleyConstruction::new(F121_AGE_K, Some(0.1))
    };
    let on = ValleyConstruction { gorge_retreat: None, ..vc };
    let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
    let (b, _) = carve(&s1, &sk, &vc, &ss);
    let (b0, _) = carve(&s1, &skeleton(&s1, &on, &ss, DOMAIN_KM), &on, &ss);
    let bo = sk.gorge_body_of.as_ref().unwrap();
    eprintln!("\n=== F140 diag · {} bodies", sk.gorge_bodies.len());
    for (i, g) in sk.gorge_bodies.iter().enumerate() {
        if g.level < 0.5 || g.l_floor < 0.5 {
            eprintln!(
                "   body {i} low ({},{}) {} cells · L_in {:.1} · L_bed {:.1} · L_floor {:.1} · r {:.2} · level {:.1} · A {:.0}",
                g.low as usize % w, g.low as usize / w, g.cells, g.l_in, g.l_bed, g.l_floor, g.r_lake, g.level, g.a_out_km2
            );
        }
    }
    let newsea: Vec<usize> = (0..n).filter(|&k| s1.data[k] > SEA && b.data[k] <= SEA && b0.data[k] > SEA).collect();
    eprintln!("   land cells the GORGE lays at or below the sea (and ON does not): {}", newsea.len());
    for (x, y) in [(4379usize, 3551usize), (2435, 3586)] {
        let k = y * w + x;
        let near: Vec<usize> = newsea.iter().copied().filter(|&c| ((c % w) as i64 - x as i64).abs() <= 30 && ((c / w) as i64 - y as i64).abs() <= 30).collect();
        let bodies: std::collections::BTreeSet<u32> = near.iter().map(|&c| bo[c]).collect();
        eprintln!(
            "   floor ({x},{y}): S1 {:.1} m · GORGE {:.1} m · ON {:.1} m · body {} · new below-sea cells within 30 cells {} (bodies {:?})",
            c1_altitude_norm_to_metres(s1.data[k], &ss),
            c1_altitude_norm_to_metres(b.data[k], &ss),
            c1_altitude_norm_to_metres(b0.data[k], &ss),
            bo[k],
            near.len(),
            bodies
        );
    }
}

/// ADR Finding 140 — the failure, stage by stage (rule 12): the two floor cells of Finding 38's uncovered components
/// ((4379, 3551) and (2435, 3586)) and their 30-cell neighbourhoods, through the construction (+ rims), the light pass,
/// the droplets, the bathymetry and the protected breach, GORGE ×1 p 0.5 against ON.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f140_diag2 --nocapture
#[test]
#[ignore]
fn f140_diag2() {
    use common::{build_world, viz_dcfg};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    use ymir_core::terrain::flow::breach_monotone_protected;
    let ss = SteinSteinParams::default();
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let gz = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)), ..on };
    let cells = [(4379usize, 3551usize), (2435, 3586)];
    for (label, vc) in [("ON", on), ("GORGE", gz)] {
        eprintln!("\n=== {label}");
        let b = Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
        for (st, kn) in [
            ("(a) construction + rims", Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..b }),
            ("(b) + light pass", Knobs { erosion_off: true, bathymetry_off: true, ..b }),
            ("(c) + droplets", Knobs { bathymetry_off: true, ..b }),
            ("(d) + bathymetry", b),
        ] {
            let wd = build_world(kn, None, PSEED, None);
            let g = &wd.heightmap;
            let w = g.width;
            let mut line = Vec::new();
            for &(x, y) in &cells {
                let k = y * w + x;
                let below = (0..61 * 61)
                    .filter(|&i| {
                        let (xx, yy) = (x + i % 61 - 30, y + i / 61 - 30);
                        g.data[yy * w + xx] <= SEA
                    })
                    .count();
                line.push(format!("({x},{y}) {:.1} m, ≤ sea within 30 cells {below}", c1_altitude_norm_to_metres(g.data[k], &ss)));
            }
            eprintln!("   {st:<26} {}", line.join(" · "));
            if st.starts_with("(d)") {
                let (gw, gh) = (g.width, g.height);
                let d = c1_drainage_windowed(g, None, &viz_dcfg(), &ss, DOMAIN_KM);
                let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, gw, gh));
                let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, gw, gh, prot.as_deref());
                let mut line = Vec::new();
                for &(x, y) in &cells {
                    let k = y * gw + x;
                    line.push(format!("({x},{y}) {:.1} m (pre-breach lake id {})", c1_altitude_norm_to_metres(c.data[k], &ss), d.lake_map[k]));
                }
                eprintln!("   {:<26} {}", "(e) protected breach", line.join(" · "));
            }
        }
    }
}

/// ADR Finding 140 — the failure's mechanism: body 4 (the floor cell (4379, 3551) lies in it). Its levels, how much of
/// its footprint the pre-breach drainage still holds as a lake (ON against GORGE), and the cells the protected breach
/// lowers below the sea there.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f140_diag3 --nocapture
#[test]
#[ignore]
fn f140_diag3() {
    use common::{build_world, viz_dcfg};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    use ymir_core::terrain::flow::breach_monotone_protected;
    let ss = SteinSteinParams::default();
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let w = s1.width;
    let n = s1.data.len();
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let gz = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)), ..on };
    let sk = skeleton(&s1, &gz, &ss, DOMAIN_KM);
    let bo = sk.gorge_body_of.clone().unwrap();
    let bi = bo[3551 * w + 4379];
    let b = sk.gorge_bodies[bi as usize].clone();
    let foot: Vec<usize> = (0..n).filter(|&k| bo[k] == bi).collect();
    eprintln!(
        "\n=== body {bi}: {} cells ({:.1} km²) · L_in {:.1} · L_bed {:.1} · L_floor {:.1} · A {:.0} km² · r_lake {:.2} · level {:.1} · slope {:.3} · head {:.1} m · col ({},{})",
        b.cells, b.cells as f32 * CELL_KM * CELL_KM, b.l_in, b.l_bed, b.l_floor, b.a_out_km2, b.r_lake, b.level, b.slope, b.head_fall_m, b.col as usize % w, b.col as usize / w
    );
    let below_level = foot.iter().filter(|&&k| c1_altitude_norm_to_metres(s1.data[k], &ss) < b.level).count();
    eprintln!("   S1 footprint cells below L(r): {below_level} ({:.1} km²)", below_level as f32 * CELL_KM * CELL_KM);
    drop(sk);
    for (label, vc) in [("ON", on), ("GORGE", gz)] {
        let kn = Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
        let wd = build_world(kn, None, PSEED, None);
        let g = &wd.heightmap;
        let (gw, gh) = (g.width, g.height);
        let d = c1_drainage_windowed(g, None, &viz_dcfg(), &ss, DOMAIN_KM);
        let held = foot.iter().filter(|&&k| d.lake_map[k] != 0).count();
        let ids: std::collections::BTreeSet<u32> = foot.iter().map(|&k| d.lake_map[k]).filter(|&x| x != 0).collect();
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, gw, gh));
        let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, gw, gh, prot.as_deref());
        let newsea: Vec<usize> = (0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect();
        let near = newsea.iter().filter(|&&k| ((k % w) as i64 - 4379).abs() <= 200 && ((k / w) as i64 - 3551).abs() <= 200).count();
        let mut lowest = (f32::INFINITY, 0usize);
        for &k in &foot {
            let z = c1_altitude_norm_to_metres(g.data[k], &ss);
            if z < lowest.0 && d.lake_map[k] == 0 {
                lowest = (z, k);
            }
        }
        eprintln!(
            "   {label}: pre-breach lakes over the footprint {held} cells ({:.1} km², ids {:?}) · land cells the breach lowers to ≤ sea: {} (within 200 cells of the floor cell: {near}) · lowest UNHELD footprint cell in the world {:.1} m at ({},{})",
            held as f32 * CELL_KM * CELL_KM,
            ids,
            newsea.len(),
            lowest.0,
            lowest.1 % w,
            lowest.1 / w
        );
    }
}

/// ADR Finding 141 — F140's crash attributed, nothing corrected. A by three definitions for the lakes with a body
/// (Part A); the land cells the protected breach lowers below the sea, through an instrumented copy of the breach
/// checked bit for bit against production (Part Br, ON and GORGE ×1); the cell 78 m under body 4's drained floor, stage
/// by stage (Part Z). Declared in `f141_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f141_attr --nocapture
#[test]
#[ignore]
fn f141_attr() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{BTreeMap, BinaryHeap, HashMap};
    use ymir_core::tectonics_c1::drainage::LakeType;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase, carve_diag};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

    // ── the instrumented copy of `breach_monotone_protected` (flow.rs), line for line, plus the records
    struct Pq(f32, usize);
    impl PartialEq for Pq {
        fn eq(&self, o: &Self) -> bool {
            self.0 == o.0 && self.1 == o.1
        }
    }
    impl Eq for Pq {}
    impl PartialOrd for Pq {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for Pq {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            o.0.partial_cmp(&self.0).unwrap_or(std::cmp::Ordering::Equal).then_with(|| o.1.cmp(&self.1))
        }
    }
    /// A ramp that lowered at least one land cell to ≤ sea.
    struct Ev {
        pit: usize,
        ci: usize,
        floor: f32,
        steps: u32,
        stop: u8, // 0 a base, 1 a cell already lower than the ramp, 2 the chain's end
        stop_cell: usize,
        crossed: u32,
        cross_step: u32,
    }
    struct Br {
        z: Vec<f32>,
        pit_of: Vec<u32>,
        step_of: Vec<u32>,
        backlink: Vec<usize>,
        evs: Vec<Ev>,
    }
    fn breach_instr(height: &GridF32, filled: &GridF32, lake_map: &[u32], sea_level: f32, w: usize, h: usize, protect: Option<&[bool]>) -> Br {
        use ymir_core::terrain::flow::{D8_DX, D8_DY};
        let n = w * h;
        let prot = |k: usize| protect.is_some_and(|p| p[k]);
        let mut z = height.data.clone();
        for k in 0..n {
            if lake_map[k] != 0 && !prot(k) {
                z[k] = filled.data[k];
            }
        }
        let is_base = |k: usize, z: &[f32]| z[k] <= sea_level || lake_map[k] != 0;
        const EPS: f32 = 1e-5;
        let mut visited = vec![false; n];
        let mut backlink = vec![usize::MAX; n];
        let mut pit_of = vec![u32::MAX; n];
        let mut step_of = vec![0u32; n];
        let mut evs: Vec<Ev> = Vec::new();
        let mut heap: BinaryHeap<Pq> = BinaryHeap::new();
        for k in 0..n {
            if prot(k) {
                visited[k] = true;
            }
        }
        for k in 0..n {
            if !visited[k] && is_base(k, &z) {
                visited[k] = true;
                heap.push(Pq(z[k], k));
            }
        }
        while let Some(Pq(_, ci)) = heap.pop() {
            let (cx, cy) = (ci % w, ci / w);
            for d in 0..8 {
                let nx = ((cx as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                let ny = ((cy as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                let nb = ny * w + nx;
                if visited[nb] {
                    continue;
                }
                visited[nb] = true;
                backlink[nb] = ci;
                if height.data[nb] < z[ci] {
                    let mut target = height.data[nb] - EPS;
                    let mut cur = ci;
                    let (mut steps, mut crossed, mut cross_step) = (0u32, 0u32, u32::MAX);
                    while cur != usize::MAX && !is_base(cur, &z) && z[cur] > target {
                        if target <= sea_level {
                            crossed += 1;
                            if cross_step == u32::MAX {
                                cross_step = steps;
                            }
                        }
                        z[cur] = target;
                        pit_of[cur] = nb as u32;
                        step_of[cur] = steps;
                        target -= EPS;
                        cur = backlink[cur];
                        steps += 1;
                    }
                    if crossed > 0 {
                        let stop = if cur == usize::MAX {
                            2
                        } else if is_base(cur, &z) {
                            0
                        } else {
                            1
                        };
                        evs.push(Ev { pit: nb, ci, floor: height.data[nb], steps, stop, stop_cell: cur, crossed, cross_step });
                    }
                }
                z[nb] = height.data[nb];
                heap.push(Pq(z[nb], nb));
            }
        }
        let mut visited2 = vec![false; n];
        let mut heap2: BinaryHeap<Pq> = BinaryHeap::new();
        for k in 0..n {
            if prot(k) {
                visited2[k] = true;
            }
        }
        for k in 0..n {
            if !visited2[k] && is_base(k, &z) {
                visited2[k] = true;
                heap2.push(Pq(z[k], k));
            }
        }
        while let Some(Pq(_, ci)) = heap2.pop() {
            let (cx, cy) = (ci % w, ci / w);
            for d in 0..8 {
                let nx = ((cx as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                let ny = ((cy as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                let nb = ny * w + nx;
                if visited2[nb] {
                    continue;
                }
                visited2[nb] = true;
                if z[nb] < z[ci] {
                    z[nb] = z[ci];
                }
                heap2.push(Pq(z[nb], nb));
            }
        }
        Br { z, pit_of, step_of, backlink, evs }
    }

    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let mm = |v: f32| c1_altitude_norm_to_metres(v, &ss);
    let a_ref = 418.7f32;
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let gz = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)), ..on };
    let kb = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    eprintln!("\n==========  Finding 141 . F140's crash attributed (nothing corrected)  ==========");

    // ── the construction's input and GORGE ×1's skeleton (the bodies, (i), (iii))
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let sk = skeleton(&s1, &gz, &ss, DOMAIN_KM);
    let bo = sk.gorge_body_of.clone().expect("the gate is on");
    let nbodies = sk.gorge_bodies.len();
    let srecv = |c: usize| -> Option<usize> {
        let d = sk.direction[c];
        if d == DIR_NONE {
            return None;
        }
        Some(((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize)
    };
    let nb8 = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    // (iii): the inflow into each body, and its largest inflow cell
    let mut inflow = vec![0f32; nbodies];
    let mut big_in: Vec<(f32, usize)> = vec![(0.0, usize::MAX); nbodies];
    for c in 0..n {
        if let Some(r) = srecv(c)
            && bo[r] != u32::MAX
            && bo[c] != bo[r]
        {
            let b = bo[r] as usize;
            inflow[b] += sk.area_km2[c];
            if sk.area_km2[c] > big_in[b].0 {
                big_in[b] = (sk.area_km2[c], c);
            }
        }
    }
    let a_iii: Vec<f32> = (0..nbodies).map(|b| inflow[b] + sk.gorge_bodies[b].cells as f32 * cell_km2).collect();
    let bi4 = bo[3551 * w + 4379];
    eprintln!("   GORGE ×1 skeleton: {nbodies} bodies · body 4 check: the floor cell (4379,3551) is in body {bi4}");

    // ── Part A — ON ×1's final lakes and F139's instrument (ii)
    {
        let wd = build_world(kb(on), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kb(on), PSEED, 45.0, 40.0);
        let dr = &v.drainage;
        let cn: Vec<f32> = v.conditioned.data.iter().map(|&x| mm(x)).collect();
        let acc = &dr.flow.accumulation.data;
        let wrecv = |c: usize| -> Option<usize> {
            let k = dr.flow.direction[c];
            if k == DIR_NONE {
                return None;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
        };
        let step_len = |a: usize, b: usize| CELL_KM * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
        // F139's path, only its first km is read
        let extend = |path: &mut Vec<usize>, d: &mut Vec<f32>, me: u32| {
            let mut c = *path.last().unwrap();
            loop {
                if cn[c] <= 0.0 || *d.last().unwrap() > 1.0 {
                    return;
                }
                let Some(r) = wrecv(c) else { return };
                d.push(d.last().unwrap() + step_len(c, r));
                path.push(r);
                if dr.lake_map[r] != 0 && dr.lake_map[r] != me {
                    return;
                }
                c = r;
            }
        };
        let segs = &dr.rivers.segments;
        let mut fcells: HashMap<u32, Vec<usize>> = HashMap::new();
        for c in 0..n {
            if dr.lake_map[c] != 0 {
                fcells.entry(dr.lake_map[c]).or_default().push(c);
            }
        }
        struct Row {
            id: u32,
            d8: bool,
            body: u32,
            a1: f32,
            a2: f32,
            a2_path: f32,
            a2_in: f32,
            a3: f32,
            outlet: usize,
        }
        let mut rows: Vec<Row> = Vec::new();
        let mut no_body = Vec::new();
        for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)) {
            let id = l.base.id;
            let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
            let (path, d, d8) = if let Some(col) = wrecv(o) {
                let (mut path, mut d) = (vec![col], vec![0f32]);
                extend(&mut path, &mut d, id);
                (path, d, true)
            } else if (1_000_000..2_000_000).contains(&id) {
                let sp: Vec<usize> = (0..segs.len())
                    .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway && dr.segment_source_lake.get(i).copied().flatten() == Some(id))
                    .collect();
                if sp.is_empty() {
                    continue;
                }
                let mut path: Vec<usize> = Vec::new();
                let mut si = Some(sp[0]);
                let mut guard = 0;
                while let Some(i) = si {
                    for &(x, y) in &segs[i].points {
                        let p = y as usize * w + x as usize;
                        if path.last() != Some(&p) {
                            path.push(p);
                        }
                    }
                    si = segs[i].downstream;
                    guard += 1;
                    if guard > 10_000 {
                        break;
                    }
                }
                let mut d = vec![0f32; path.len()];
                for i in 1..path.len() {
                    d[i] = d[i - 1] + step_len(path[i - 1], path[i]);
                }
                let hit = (1..path.len()).find(|&i| (dr.lake_map[path[i]] != 0 && dr.lake_map[path[i]] != id) || cn[path[i]] <= 0.0);
                match hit {
                    Some(i) => {
                        path.truncate(i + 1);
                        d.truncate(i + 1);
                    }
                    None => extend(&mut path, &mut d, id),
                }
                (path, d, false)
            } else {
                continue;
            };
            let a2_in: f32 = (0..n).filter(|&c| dr.lake_map[c] != id && wrecv(c).is_some_and(|r| dr.lake_map[r] == id)).map(|c| acc[c] * cell_km2).sum::<f32>() + l.area_km2;
            let a2_path = path.iter().zip(d.iter()).filter(|e| *e.1 <= 1.0).map(|e| acc[*e.0] * cell_km2).fold(0f32, f32::max);
            let a2 = a2_path.max(a2_in).max(0.1);
            let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
            for &c in fcells.get(&id).map_or(&[][..], |v| v.as_slice()) {
                if bo[c] != u32::MAX {
                    *ov.entry(bo[c]).or_insert(0) += 1;
                }
            }
            let Some((&b, _)) = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) else {
                no_body.push(id);
                continue;
            };
            rows.push(Row { id, d8, body: b, a1: sk.gorge_bodies[b as usize].a_out_km2, a2, a2_path, a2_in, a3: a_iii[b as usize], outlet: o });
        }
        let r_of = |a: f32, rw: f32| (rw * (a / a_ref).powf(0.5)).min(2.0);
        let age_of = |r: f32| if r <= 1.0 { 0.7 + 0.3 * r } else { 1.0 + 0.4 * (r - 1.0) };
        let empt = |a: f32| -> String {
            let re = 2.0 / (a / a_ref).max(1e-9).powf(0.5);
            if re <= 2.0 + 1e-6 { format!("×{:.2}", age_of(re)) } else { format!("never (r_world {re:.2})") }
        };
        let rs = |a: f32| format!("{:.2}/{:.2}/{:.2}/{:.2}", r_of(a, 0.0), r_of(a, 1.0), r_of(a, 1.5), r_of(a, 2.0));
        eprintln!(
            "\n   A · {} lakes with a body ({} D8, {} spillway) · without a construction body: {:?} · p = 0.5, A_ref {a_ref} km² · r_lake at ×0.7/×1/×1.2/×1.4",
            rows.len(),
            rows.iter().filter(|r| r.d8).count(),
            rows.iter().filter(|r| !r.d8).count(),
            no_body
        );
        let (mut near15, mut far2) = (0usize, 0usize);
        for r in &rows {
            let q32 = r.a3 / r.a2;
            let q12 = r.a1 / r.a2;
            if (1.0 / 1.5..=1.5).contains(&q32) {
                near15 += 1;
            }
            let div = !(0.5..=2.0).contains(&q12);
            if div {
                far2 += 1;
            }
            let gb = &sk.gorge_bodies[r.body as usize];
            eprintln!(
                "      lake {:>7} ({}) · body {:>3} ({:.1} km², construction r_lake {:.2}) · (i) {:>7.1} · (ii) {:>7.1} [path {:.1}, inflow {:.1}] · (iii) {:>7.1} km² · (iii)/(ii) {q32:.2} · (i)/(ii) {q12:.2}{}",
                r.id,
                if r.d8 { "D8" } else { "spillway" },
                r.body,
                gb.cells as f32 * cell_km2,
                gb.r_lake,
                r.a1,
                r.a2,
                r.a2_path,
                r.a2_in,
                r.a3,
                if div { " · **> ×2**" } else { "" }
            );
            eprintln!(
                "         r_lake (i) {} · (ii) {} · (iii) {} · emptying (i) {} · (ii) {} · (iii) {}",
                rs(r.a1),
                rs(r.a2),
                rs(r.a3),
                empt(r.a1),
                empt(r.a2),
                empt(r.a3)
            );
        }
        eprintln!("   A · (iii) within ×1.5 of (ii): **{near15} of {}** · (i) beyond ×2 of (ii): **{far2} of {}**", rows.len(), rows.len());
        for (lbl, f) in [("(i)", 0usize), ("(ii)", 1), ("(iii)", 2)] {
            let mut empty_at = [0usize; 4];
            for r in &rows {
                let a = [r.a1, r.a2, r.a3][f];
                for (j, rw) in [0f32, 1.0, 1.5, 2.0].into_iter().enumerate() {
                    if r_of(a, rw) >= 2.0 - 1e-6 {
                        empty_at[j] += 1;
                    }
                }
            }
            eprintln!("   A · {lbl}: lakes drained (r_lake = 2) at ×0.7/×1/×1.2/×1.4: {:?}", empty_at);
        }
        // the gap's stream: lake 4 (body bi4) and every lake where (i) and (ii) differ by more than ×2
        let touches = |c: usize| -> Option<u32> {
            if dr.lake_map[c] != 0 {
                return Some(dr.lake_map[c]);
            }
            (0..8).map(|k| dr.lake_map[nb8(c, k)]).find(|&x| x != 0)
        };
        eprintln!("\n   A · the gap's stream (the construction's col and its donors; the body's largest inflow; the final lake's outlet):");
        for r in rows.iter().filter(|r| r.body == bi4 || !(0.5..=2.0).contains(&(r.a1 / r.a2))) {
            let gb = &sk.gorge_bodies[r.body as usize];
            let (ex, col) = (gb.exit as usize, gb.col as usize);
            eprintln!(
                "      lake {} (body {}{}) · exit {} {:.1} km² · col {} {:.1} km² (in final lake: {:?}) · final outlet {} · col ↔ final outlet {:.1} km",
                r.id,
                r.body,
                if r.body == bi4 { ", = body 4" } else { "" },
                xy(ex),
                sk.area_km2[ex],
                xy(col),
                sk.area_km2[col],
                touches(col),
                xy(r.outlet),
                ((((col % w) as f32 - (r.outlet % w) as f32).powi(2) + ((col / w) as f32 - (r.outlet / w) as f32).powi(2)).sqrt()) * CELL_KM
            );
            let mut donors: Vec<usize> = (0..8).map(|k| nb8(col, k)).filter(|&d| srecv(d) == Some(col)).collect();
            donors.sort_by(|a, b| sk.area_km2[*b].total_cmp(&sk.area_km2[*a]));
            for d in donors.iter().take(4) {
                eprintln!(
                    "         col donor {} · {:.1} km² · {} · touches a final lake: {:?}",
                    xy(*d),
                    sk.area_km2[*d],
                    if bo[*d] == r.body { "IN the body (the exit's flow)".to_string() } else if bo[*d] == u32::MAX { "outside every body".to_string() } else { format!("in body {}", bo[*d]) },
                    touches(*d)
                );
            }
            let (ba, bc) = big_in[r.body as usize];
            if bc != usize::MAX {
                eprintln!(
                    "         largest inflow into the body {} · {:.1} km² ({:.0} % of (iii)) · touches the FINAL lake {}: {}",
                    xy(bc),
                    ba,
                    100.0 * ba / r.a3,
                    r.id,
                    touches(bc) == Some(r.id)
                );
            }
            // where the col's area goes in the final world: the col cell's ON accumulation
            eprintln!(
                "         ON world at the construction's col: accumulation {:.1} km² · lake id {} · the final lake's footprint {:.1} km² of the body's {:.1}",
                acc[col] * cell_km2,
                dr.lake_map[col],
                fcells.get(&r.id).map_or(&[][..], |v| v.as_slice()).iter().filter(|&&c| bo[c] == r.body).count() as f32 * cell_km2,
                gb.cells as f32 * cell_km2
            );
        }

        // ── Part Br, ON ×1
        br_report("ON ×1", &wd, Some(dr.lake_map.as_slice()), &bo, &ss, w, h);
    }

    // ── Part Br + Z — GORGE ×1, stage by stage
    let zc = 3674 * w + 4148;
    let gb4 = sk.gorge_bodies[bi4 as usize].clone();
    let foot4: Vec<usize> = (0..n).filter(|&k| bo[k] == bi4).collect();
    eprintln!(
        "\n   Z · cell {} · body {} (body 4 = {bi4}) · L_in {:.1} · L_bed {:.1} · L_floor {:.1} · L(r) {:.1} m · r_lake {:.2}",
        xy(zc),
        bo[zc],
        gb4.l_in,
        gb4.l_bed,
        gb4.l_floor,
        gb4.level,
        gb4.r_lake
    );
    let foot_min = |g: &[f32]| -> (f32, usize, usize) {
        let mut lo = (f32::INFINITY, 0usize);
        let mut under = 0usize;
        for &k in &foot4 {
            let z = mm(g[k]);
            if z < gb4.l_floor - 1.0 {
                under += 1;
            }
            if z < lo.0 {
                lo = (z, k);
            }
        }
        (lo.0, lo.1, under)
    };
    let zline = |lbl: &str, g: &[f32]| {
        let (lz, lk, under) = foot_min(g);
        let nmin = (0..8).map(|k| mm(g[nb8(zc, k)])).fold(f32::INFINITY, f32::min);
        eprintln!(
            "      {lbl:<34} z {:>7.1} m (neighbours' min {nmin:>7.1}) · body 4's lowest {lz:>7.1} m at {} · body-4 cells > 1 m under L_floor: {under}",
            mm(g[zc]),
            xy(lk)
        );
    };
    zline("S1 (the construction's input)", &s1.data);
    {
        let (cout, _, dg) = carve_diag(&s1, &sk, &gz, &ss);
        zline("carve_diag(S1), GORGE skeleton", &cout.data);
        let sk_on = skeleton(&s1, &on, &ss, DOMAIN_KM);
        let (cout_on, _, _) = carve_diag(&s1, &sk_on, &on, &ss);
        zline("carve_diag(S1), ON skeleton", &cout_on.data);
        drop(sk_on);
        eprintln!(
            "      the cell's skeleton: base {:.1} m · χ {:.1} · floor_m {:.1} · area {:.2} km² · rim_floor {:?}",
            sk.base_alt_m[zc],
            sk.chi_m[zc],
            sk.floor_m(zc, gz.age_k),
            sk.area_km2[zc],
            sk.rim_floor_m.as_ref().map(|r| r[zc])
        );
        let s = dg.who[zc];
        if s == u32::MAX {
            eprintln!("      who: none (no sample reached the cell)");
        } else {
            let (li, pi) = (dg.line_of[s as usize] as usize, dg.pos_of[s as usize] as usize);
            let line = &sk.polylines[li];
            let (sx, sy, zf, hw) = line[pi];
            let sc = (sy.floor() as i64).rem_euclid(h as i64) as usize * w + (sx.floor() as i64).rem_euclid(w as i64) as usize;
            let dist = ((((zc % w) as f32 + 0.5 - sx).powi(2) + ((zc / w) as f32 + 0.5 - sy).powi(2)).sqrt()) * CELL_KM * 1000.0;
            eprintln!(
                "      who: sample {s} · line {li} pos {pi} of {} · at ({sx:.1},{sy:.1}) cell {} · zf {zf:.1} m · half-width {hw:.0} m · {dist:.0} m from the cell · banded {:?}",
                line.len(),
                xy(sc),
                dg.banded.as_ref().map(|b| b[zc])
            );
            eprintln!(
                "         the sample's cell: base {:.1} · χ {:.1} · floor_m {:.1} · area {:.1} km² · body {}",
                sk.base_alt_m[sc],
                sk.chi_m[sc],
                sk.floor_m(sc, gz.age_k),
                sk.area_km2[sc],
                if bo[sc] == u32::MAX { "none".to_string() } else { bo[sc].to_string() }
            );
            let in4 = line.iter().filter(|p| bo[(p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize] == bi4).count();
            let zmin = line.iter().map(|p| p.2).fold(f32::INFINITY, f32::min);
            let (f0, fl) = (line[0], line[line.len() - 1]);
            eprintln!(
                "         its line: {} samples, {in4} on body-4 cells · from ({:.0},{:.0}) zf {:.1} to ({:.0},{:.0}) zf {:.1} · min zf {zmin:.1} · parent {:?}",
                line.len(),
                f0.0,
                f0.1,
                f0.2,
                fl.0,
                fl.1,
                fl.2,
                sk.line_parent.get(li)
            );
        }
    }
    let b = kb(gz);
    for (st, kn) in [
        ("(a) construction + rims", Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..b }),
        ("(b) + light pass", Knobs { erosion_off: true, bathymetry_off: true, ..b }),
        ("(c) + droplets", Knobs { bathymetry_off: true, ..b }),
        ("(d) + bathymetry", b),
    ] {
        let wd = build_world(kn, None, PSEED, None);
        zline(st, &wd.heightmap.data);
        if st.starts_with("(d)") {
            let z = br_report("GORGE ×1", &wd, None, &bo, &ss, w, h);
            zline("(e) protected breach", &z);
        }
    }
    eprintln!("\n==========  end Finding 141 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());

    /// Part Br for one world: the instrumented breach (checked bit for bit), its below-sea land cells, their ramps.
    /// `final_lakes`: the run_hd tail's lake map when it completes (ON), `None` when it panics (GORGE).
    fn br_report(label: &str, wd: &common::World, final_lakes: Option<&[u32]>, bo: &[u32], ss: &SteinSteinParams, w: usize, h: usize) -> Vec<f32> {
        use std::collections::BTreeMap;
        use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
        use ymir_core::terrain::flow::breach_monotone_protected;
        let mm = |v: f32| c1_altitude_norm_to_metres(v, ss);
        let xy = |c: usize| format!("({},{})", c % w, c / w);
        let g = &wd.heightmap;
        let n = w * h;
        let d = c1_drainage_windowed(g, None, &common::viz_dcfg(), ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let prod = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        let br = breach_instr(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        let diff = (0..n).filter(|&k| prod.data[k].to_bits() != br.z[k].to_bits()).count();
        eprintln!("\n   Br · {label} · the instrumented copy against production: {diff} cells differ{}", if diff == 0 { " (bit-identical: the records are readable)" } else { " — THE RECORDS ARE NOT READ" });
        if diff != 0 {
            return prod.data;
        }
        let pop: Vec<usize> = (0..n).filter(|&k| g.data[k] > SEA && br.z[k] <= SEA).collect();
        let wc = water_class(&GridF32 { width: w, height: h, data: br.z.clone() }, SEA);
        let inland = pop.iter().filter(|&&k| wc[k] == WATER_CLASS_INLAND).count();
        let carved = pop.iter().filter(|&&k| br.pit_of[k] != u32::MAX).count();
        let in_body = pop.iter().filter(|&&k| bo[k] != u32::MAX).count();
        eprintln!(
            "      land cells (eroded > sea) ≤ sea after the breach: **{}** · laid by a ramp: {carved} · inland after the breach: {inland} · ocean-connected: {} · on a construction body (GORGE skeleton's bodies): {in_body}",
            pop.len(),
            pop.len() - inland
        );
        let pre_cov = pop.iter().filter(|&&k| d.lake_map[k] != 0).count();
        eprintln!("      covered by a PRE-breach lake: {pre_cov}");
        if let Some(fl) = final_lakes {
            let mut ids: BTreeMap<u32, usize> = BTreeMap::new();
            for &k in &pop {
                if fl[k] != 0 {
                    *ids.entry(fl[k]).or_insert(0) += 1;
                }
            }
            let cov: usize = ids.values().sum();
            let inl_unc = pop.iter().filter(|&&k| fl[k] == 0 && wc[k] == WATER_CLASS_INLAND).count();
            eprintln!(
                "      covered by a FINAL lake: **{cov} of {}** ({:.1} %) · by lake {:?} · uncovered and inland: {inl_unc}",
                pop.len(),
                100.0 * cov as f64 / pop.len().max(1) as f64,
                ids
            );
        } else {
            let inl_unc = pop.iter().filter(|&&k| d.lake_map[k] == 0 && wc[k] == WATER_CLASS_INLAND).count();
            eprintln!("      the run_hd tail panics on this world (F140): no final lake map · uncovered by a pre-breach lake and inland: {inl_unc}");
        }
        // the ramps that crossed the sea
        let mut kinds: BTreeMap<String, (usize, u32)> = BTreeMap::new();
        let base_kind = |c: usize| -> String {
            if c == usize::MAX {
                return "the chain's end".to_string();
            }
            if d.lake_map[c] != 0 {
                return format!("pre-breach lake {}", d.lake_map[c]);
            }
            if g.data[c] <= SEA {
                if wc[c] == WATER_CLASS_INLAND { "an inland below-sea cell".to_string() } else { "the sea".to_string() }
            } else {
                "a land cell an earlier ramp lowered ≤ sea".to_string()
            }
        };
        let mut floors_ge_spill = 0usize;
        for e in &br.evs {
            let k = format!("stop {} at {}", ["base", "lower cell", "chain end"][e.stop as usize], base_kind(e.stop_cell));
            let en = kinds.entry(k).or_insert((0, 0));
            en.0 += 1;
            en.1 += e.crossed;
            if mm(d.flow.filled.data[e.pit]) > mm(e.floor) + 0.5 {
                floors_ge_spill += 1;
            }
        }
        eprintln!(
            "      ramps that crossed the sea: {} (crossed cells {}) · pits whose spill is > 0.5 m above their floor: {floors_ge_spill}",
            br.evs.len(),
            br.evs.iter().map(|e| e.crossed).sum::<u32>()
        );
        for (k, (c, x)) in &kinds {
            eprintln!("         {k}: {c} ramps, {x} cells");
        }
        let mut ev: Vec<&Ev> = br.evs.iter().collect();
        ev.sort_by(|a, b| b.crossed.cmp(&a.crossed).then(a.pit.cmp(&b.pit)));
        let mut steps_all: Vec<u32> = pop.iter().filter(|&&k| br.pit_of[k] != u32::MAX).map(|&k| br.step_of[k]).collect();
        steps_all.sort_unstable();
        if !steps_all.is_empty() {
            eprintln!(
                "      ramp steps from the pit at the below-sea cells: min {} · p50 {} · max {} (EPS = 1e-5 norm ≈ {:.3} m per step)",
                steps_all[0],
                steps_all[steps_all.len() / 2],
                steps_all[steps_all.len() - 1],
                mm(0.5 + 1e-5) - mm(0.5)
            );
        }
        eprintln!("      the top ramps by cells laid ≤ sea:");
        for e in ev.iter().take(12) {
            // the chain's terminal base from the stop cell
            let mut t = e.stop_cell;
            let mut guard = 0usize;
            while t != usize::MAX && br.backlink[t] != usize::MAX && guard < n {
                t = br.backlink[t];
                guard += 1;
            }
            let cov = final_lakes.map(|fl| pop.iter().filter(|&&k| br.pit_of[k] == e.pit as u32 && fl[k] != 0).count());
            eprintln!(
                "         pit {} floor {:.1} m · spill (pre-breach filled) {:.1} m · pre-breach lake at the pit {} · outlet {} {:.1} m · {} steps, the ramp crossed the sea at step {} · stops at {} ({}) · {} cells ≤ sea · body of the pit {} · final-lake cover {:?} · the flood chain's root {} ({})",
                xy(e.pit),
                mm(e.floor),
                mm(d.flow.filled.data[e.pit]),
                d.lake_map[e.pit],
                xy(e.ci),
                mm(g.data[e.ci]),
                e.steps,
                e.cross_step,
                if e.stop_cell == usize::MAX { "—".to_string() } else { xy(e.stop_cell) },
                base_kind(e.stop_cell),
                e.crossed,
                if bo[e.pit] == u32::MAX { "none".to_string() } else { bo[e.pit].to_string() },
                cov,
                if t == usize::MAX { "—".to_string() } else { xy(t) },
                if t == usize::MAX { "—".to_string() } else { base_kind(t) }
            );
        }
        br.z
    }
}

/// ADR Finding 141-A, amendment — body 4 (the floor cell (4379, 3551)) matches none of F139's 22 lakes. Which ON final
/// lake stands on it, why F139's set left it out, and its A by the three definitions. Declared in `f141_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f141_b4 --nocapture
#[test]
#[ignore]
fn f141_b4() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::BTreeMap;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let gz = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)), ..on };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let sk = skeleton(&s1, &gz, &ss, DOMAIN_KM);
    let bo = sk.gorge_body_of.clone().unwrap();
    let bi = bo[3551 * w + 4379];
    let b = sk.gorge_bodies[bi as usize].clone();
    let srecv = |c: usize| -> Option<usize> {
        let d = sk.direction[c];
        if d == DIR_NONE {
            return None;
        }
        Some(((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize)
    };
    let nb8 = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let (mut inflow, mut big) = (0f32, (0f32, usize::MAX));
    let mut ins: Vec<(f32, usize)> = Vec::new();
    for c in 0..n {
        if bo[c] != bi
            && let Some(r) = srecv(c)
            && bo[r] == bi
        {
            inflow += sk.area_km2[c];
            ins.push((sk.area_km2[c], c));
            if sk.area_km2[c] > big.0 {
                big = (sk.area_km2[c], c);
            }
        }
    }
    ins.sort_by(|a, b| b.0.total_cmp(&a.0));
    let a3 = inflow + b.cells as f32 * cell_km2;
    let (ex, col) = (b.exit as usize, b.col as usize);
    eprintln!(
        "\n=== body {bi} · {:.1} km² · L_in {:.1} · L_bed {:.1} · L_floor {:.1} · r_lake {:.2} · level {:.1} · (i) {:.1} km² · (iii) {a3:.1} km² · exit {} {:.1} km² · col {} {:.1} km²",
        b.cells as f32 * cell_km2, b.l_in, b.l_bed, b.l_floor, b.r_lake, b.level, b.a_out_km2, xy(ex), sk.area_km2[ex], xy(col), sk.area_km2[col]
    );
    for (a, c) in ins.iter().take(5) {
        eprintln!("   inflow {} · {:.1} km² · from body {}", xy(*c), a, if bo[*c] == u32::MAX { "none".to_string() } else { bo[*c].to_string() });
    }
    let mut donors: Vec<usize> = (0..8).map(|k| nb8(col, k)).filter(|&d| srecv(d) == Some(col)).collect();
    donors.sort_by(|a, b| sk.area_km2[*b].total_cmp(&sk.area_km2[*a]));
    for d in &donors {
        eprintln!("   col donor {} · {:.1} km² · body {}", xy(*d), sk.area_km2[*d], if bo[*d] == u32::MAX { "none".to_string() } else { bo[*d].to_string() });
    }
    // the exit's flow inside the body: the largest-area cells of the body, and which other body cells drain out
    let outs: Vec<usize> = (0..n).filter(|&c| bo[c] == bi && srecv(c).is_some_and(|r| bo[r] != bi)).collect();
    eprintln!("   body cells whose receiver leaves the body: {} (largest: {})", outs.len(), outs.iter().map(|&c| sk.area_km2[c]).fold(0f32, f32::max));
    let kn = Knobs { valley: Some(on), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let dr = &v.drainage;
    let wrecv = |c: usize| -> Option<usize> {
        let k = dr.flow.direction[c];
        if k == DIR_NONE {
            return None;
        }
        let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
    };
    let mut ids: BTreeMap<u32, usize> = BTreeMap::new();
    for c in 0..n {
        if bo[c] == bi && dr.lake_map[c] != 0 {
            *ids.entry(dr.lake_map[c]).or_insert(0) += 1;
        }
    }
    eprintln!("   ON final lakes on the body: {:?}", ids);
    let acc = &dr.flow.accumulation.data;
    for (&id, &cnt) in &ids {
        let Some(l) = dr.lakes.iter().find(|l| l.base.id == id) else {
            eprintln!("   lake {id}: in the map, not in the list");
            continue;
        };
        let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
        let sp = (0..dr.rivers.segments.len()).filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway && dr.segment_source_lake.get(i).copied().flatten() == Some(id)).count();
        let a2_in: f32 = (0..n).filter(|&c| dr.lake_map[c] != id && wrecv(c).is_some_and(|r| dr.lake_map[r] == id)).map(|c| acc[c] * cell_km2).sum::<f32>() + l.area_km2;
        let tot: usize = dr.lake_map.iter().filter(|&&x| x == id).count();
        eprintln!(
            "   lake {id} · {cnt} of its {tot} cells on the body · area {:.1} km² · type {:?} · level {:.1} m · outlet {} (receiver {:?}, the outlet's own lake id {}) · spillway segments {sp} · inflow + area {a2_in:.1} km² · accumulation at the outlet {:.1} km² · F139's set: {}",
            l.area_km2,
            l.lake_type,
            l.level_m,
            xy(o),
            wrecv(o).map(xy),
            dr.lake_map[o],
            acc[o] * cell_km2,
            if l.area_km2 < 1.0 { "no (< 1 km²)" } else if wrecv(o).is_some() { "D8 (should be in)" } else if (1_000_000..2_000_000).contains(&id) && sp > 0 { "spillway" } else { "NO: no D8 receiver at the outlet and not a spillway lake" }
        );
    }
}

/// ADR Finding 142-R — re-tabulate before building (F141's rule): T3 (level and area by r), R (lakes left and area by
/// r_world and p, the emptying ages) and T4 (head falls, steepen) for the 14 D8 lakes, with the definitions AS BUILT
/// (the bodies of `gorge_bodies`, A = the col's area, L_floor on S1), on the input alone. Also every body's class, for
/// the scope's criterion. Declared in `f142_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f142_r --nocapture
#[test]
#[ignore]
fn f142_r() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use ymir_core::tectonics_c1::drainage::LakeType;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let a_ref = 418.7f32;
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let gzp = |r: f32, p: f32| ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(r, p)), ..on };
    eprintln!("\n==========  Finding 142-R . the tables re-tabulated as built (14 D8 lakes, the input alone)  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let zs1: Vec<f32> = s1.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
    let sk = skeleton(&s1, &gzp(1.0, 0.5), &ss, DOMAIN_KM);
    let bo = sk.gorge_body_of.clone().unwrap();
    let bodies = sk.gorge_bodies.clone();
    drop(sk);
    // ── the final lakes (ON ×1), F139's D8 set, matched to the bodies by footprint
    let (d8, lake_map): (Vec<(u32, f32, f32)>, Vec<u32>) = {
        let kn = Knobs { valley: Some(on), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
        let wd = build_world(kn, None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
        let dr = &v.drainage;
        let wrecv = |c: usize| -> Option<usize> {
            let k = dr.flow.direction[c];
            if k == DIR_NONE {
                return None;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[k as usize], (c / w) as i32 + D8_DY[k as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
        };
        let mut d8 = Vec::new();
        let mut finals: HashMap<u32, (String, usize)> = HashMap::new();
        let mut cnt: HashMap<u32, usize> = HashMap::new();
        for &id in &dr.lake_map {
            if id != 0 {
                *cnt.entry(id).or_insert(0) += 1;
            }
        }
        for l in &dr.lakes {
            let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
            let crater = matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral);
            let kind = if crater {
                "crater"
            } else if l.area_km2 < 1.0 {
                "< 1 km²"
            } else if wrecv(o).is_some() {
                "D8"
            } else if (1_000_000..2_000_000).contains(&l.base.id) {
                "below-sea / spillway"
            } else {
                "other ≥ 1 km²"
            };
            if kind == "D8" {
                d8.push((l.base.id, l.level_m, l.area_km2));
            }
            finals.insert(l.base.id, (format!("{kind}, {:.1} km², level {:.1} m", l.area_km2, l.level_m), cnt.get(&l.base.id).copied().unwrap_or(0)));
        }
        // the overlap of each body with each final lake
        let mut ov: HashMap<(u32, u32), usize> = HashMap::new();
        for c in 0..n {
            if bo[c] != u32::MAX && dr.lake_map[c] != 0 {
                *ov.entry((bo[c], dr.lake_map[c])).or_insert(0) += 1;
            }
        }
        eprintln!("\n   every body (GORGE skeleton r 1, p 0.5) · its class and its final lakes by overlap (share of the body, share of the lake):");
        for (bi, b) in bodies.iter().enumerate() {
            let mut m: Vec<(u32, usize)> = ov.iter().filter(|e| e.0.0 == bi as u32).map(|e| (e.0.1, *e.1)).collect();
            m.sort_by(|a, b| b.1.cmp(&a.1));
            let ms: Vec<String> = m
                .iter()
                .take(3)
                .map(|(id, k)| {
                    let f = finals.get(id);
                    format!("lake {id} [{}] {:.0} % / {:.0} %", f.map_or("?".to_string(), |x| x.0.clone()), 100.0 * *k as f64 / b.cells as f64, 100.0 * *k as f64 / f.map_or(1, |x| x.1.max(1)) as f64)
                })
                .collect();
            eprintln!(
                "      body {bi:>2} · {} · {:>7.1} km² · ring touches a depression {} · the sea {} · L_in {:>7.1} · L_bed {:>7.1} · L_floor {:>7.1} · A {:>7.1} km² · col {} · {}",
                if b.input_lake { "input lake" } else { "depression" },
                b.cells as f32 * cell_km2,
                b.touches_depression,
                b.touches_below_sea,
                b.l_in,
                b.l_bed,
                b.l_floor,
                b.a_out_km2,
                if b.col == u32::MAX { "none".to_string() } else { xy(b.col as usize) },
                if ms.is_empty() { "no final lake".to_string() } else { ms.join(" · ") }
            );
        }
        (d8, dr.lake_map.clone())
    };
    // the D8 lakes' bodies: by the lake's largest overlap, kept if the body's share is > 50 %
    let mut rows: Vec<(u32, f32, f32, usize, f64, f64)> = Vec::new(); // (lake, today's level, area, body, body share, lake share)
    for &(id, lvl, km2) in &d8 {
        let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
        let mut tot = 0usize;
        for c in 0..n {
            if lake_map[c] == id {
                tot += 1;
                if bo[c] != u32::MAX {
                    *ov.entry(bo[c]).or_insert(0) += 1;
                }
            }
        }
        match ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) {
            Some((&b, &k)) => rows.push((id, lvl, km2, b as usize, k as f64 / bodies[b as usize].cells as f64, k as f64 / tot as f64)),
            None => eprintln!("   lake {id}: NO BODY"),
        }
    }
    let used: HashSet<usize> = rows.iter().map(|r| r.3).collect();
    eprintln!(
        "\n   G-bodies' input · {} D8 lakes · matched {} · distinct bodies {} · body share > 50 %: {}",
        d8.len(),
        rows.len(),
        used.len(),
        rows.iter().filter(|r| r.4 > 0.5).count()
    );
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let area_from = |low: usize, lv: f32| -> f64 {
        if lv <= zs1[low] {
            return 0.0;
        }
        let mut seen: HashSet<usize> = HashSet::from([low]);
        let mut q = VecDeque::from([low]);
        while let Some(c) = q.pop_front() {
            if seen.len() > 2_000_000 {
                break;
            }
            for k in 0..8 {
                let m = nbt(c, k);
                if zs1[m] < lv && seen.insert(m) {
                    q.push_back(m);
                }
            }
        }
        seen.len() as f64 * cell_km2 as f64
    };
    let lvl_at = |b: usize, r: f32| {
        let g = &bodies[b];
        if r <= 1.0 { g.l_in - r * (g.l_in - g.l_bed) } else { g.l_bed - (r - 1.0) * (g.l_bed - g.l_floor) }
    };
    // ── T3
    eprintln!("\n   T3 (as built: the body's L_in, L_bed, L_floor; the area under L(r) from the lowest cell on S1):");
    let rs = [0f32, 0.5, 1.0, 1.5, 2.0];
    let mut tot = [0f64; 5];
    for r in &rows {
        let g = &bodies[r.3];
        let mut cells = Vec::new();
        for (j, &rr) in rs.iter().enumerate() {
            let a = area_from(g.low as usize, lvl_at(r.3, rr));
            tot[j] += a;
            cells.push(format!("{:.1} m {a:.1} km²", lvl_at(r.3, rr)));
        }
        eprintln!(
            "      lake {:>3} (body {:>2}, body share {:.0} %, lake share {:.0} %) · today's level {:.1} m, area {:.1} km² · L_in {:.1} · L_bed {:.1} (Δ today {:+.1}) · L_floor {:.1} · r 0/0.5/1/1.5/2: {}",
            r.0,
            r.3,
            100.0 * r.4,
            100.0 * r.5,
            r.1,
            r.2,
            g.l_in,
            g.l_bed,
            g.l_bed - r.1,
            g.l_floor,
            cells.join(" · ")
        );
    }
    let f138 = [1308.5f64, 1223.6, 1145.7, 367.3, 0.0];
    eprintln!(
        "   T3 · total, the 14: {} · F138-T3: {}",
        tot.iter().map(|a| format!("{a:.1}")).collect::<Vec<_>>().join(" / "),
        f138.iter().map(|a| format!("{a:.1}")).collect::<Vec<_>>().join(" / ")
    );
    // ── R
    let age_of = |r: f32| if r <= 1.0 { 0.7 + 0.3 * r } else { 1.0 + 0.4 * (r - 1.0) };
    let mut a_sorted: Vec<f32> = rows.iter().map(|r| bodies[r.3].a_out_km2).collect();
    a_sorted.sort_by(f32::total_cmp);
    eprintln!(
        "\n   R · A = the col's area (as built) · the 14's median {:.1} km² (lower median, F139's rule) against the declared A_ref {a_ref}",
        a_sorted[(a_sorted.len() - 1) / 2]
    );
    let f139: [[(usize, f64); 5]; 3] = [
        [(14, 1308.5), (14, 1223.6), (14, 1145.7), (14, 367.3), (0, 0.0)],
        [(14, 1308.5), (14, 1214.3), (14, 845.5), (12, 135.4), (3, 5.2)],
        [(14, 1308.5), (14, 1202.7), (14, 521.3), (9, 103.0), (4, 19.0)],
    ];
    for (pi, p) in [0f32, 0.25, 0.5].into_iter().enumerate() {
        let mut row = Vec::new();
        for (j, &rw) in rs.iter().enumerate() {
            let (mut c, mut t) = (0usize, 0f64);
            for r in &rows {
                let rl = (rw * (bodies[r.3].a_out_km2 / a_ref).powf(p)).min(2.0);
                let a = area_from(bodies[r.3].low as usize, lvl_at(r.3, rl));
                t += a;
                if a >= 1.0 {
                    c += 1;
                }
            }
            row.push(format!("r_world {rw}: {c}, {t:.1} km² (F139 {}, {:.1})", f139[pi][j].0, f139[pi][j].1));
        }
        let mut per = Vec::new();
        let mut ages = Vec::new();
        for r in &rows {
            let re = 2.0 / (bodies[r.3].a_out_km2 / a_ref).powf(p);
            if re <= 2.0 + 1e-6 {
                ages.push(age_of(re));
                per.push(format!("{}: ×{:.2}", r.0, age_of(re)));
            } else {
                per.push(format!("{}: never ({re:.2})", r.0));
            }
        }
        ages.sort_by(f32::total_cmp);
        eprintln!("      p = {p}: {}", row.join(" · "));
        eprintln!(
            "         emptying: {} lakes, ×{:.2}–×{:.2} (spread {:.0} %) · {}",
            ages.len(),
            ages.first().copied().unwrap_or(f32::NAN),
            ages.last().copied().unwrap_or(f32::NAN),
            if ages.is_empty() { 0.0 } else { 100.0 * (ages[ages.len() - 1] - ages[0]) / 0.7 },
            per.join(" · ")
        );
    }
    // ── T4 (the construction's own falls, steepen), at r_world 1 with p 0 (r_lake 1: F139-C4's level) and p 0.5
    for p in [0f32, 0.5] {
        let sk = skeleton(&s1, &gzp(1.0, p), &ss, DOMAIN_KM);
        let mut heads = Vec::new();
        let mut tagged = 0usize;
        let mut line = Vec::new();
        for r in &rows {
            let g = &sk.gorge_bodies[r.3];
            let hf = (1.5 * g.slope * 100.0).max(10.0);
            let tag = g.head_fall_m > hf;
            if g.head_fall_m > 0.0 {
                heads.push(g.head_fall_m);
            }
            if tag {
                tagged += 1;
            }
            line.push(format!(
                "{}: φ {:.2} r {:.2} S {:.2}° D_g {:.0} head {:.1} m{} short {:.0}",
                r.0,
                g.phi,
                g.r_lake,
                g.slope.atan().to_degrees(),
                g.d_g,
                g.head_fall_m,
                if tag { " TAG" } else { "" },
                g.shortage_m
            ));
        }
        heads.sort_by(f32::total_cmp);
        eprintln!(
            "\n   T4 · r_world 1, p {p} (steepen, m 10) · head falls > 0: {} · tagged (> H_f): **{tagged}** · median {:.1} m · max {:.1} m · shortage falls: {}",
            heads.len(),
            if heads.is_empty() { f32::NAN } else { heads[heads.len() / 2] },
            heads.last().copied().unwrap_or(f32::NAN),
            rows.iter().filter(|r| sk.gorge_bodies[r.3].shortage_m > 0.5).count()
        );
        for l in line {
            eprintln!("      {l}");
        }
    }
    eprintln!("\n==========  end Finding 142-R . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 142-G — spec v4 (the D8 scope, the minimal plain) against its gates on the témoin: GORGE v4 at ×0.7 /
/// ×0.85 / ×1 / ×1.2 / ×1.4 (p = 0.5), ×1 (p = 0, 0.25) and ×1.4 (p = 0, the worst case); OFF and ON at ×1 for the
/// controls, rule 14 and rule 18. Each world runs isolated (a panic is a failed gate, the others go on). Declared in
/// `f142_declared.md` before the run.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f142_g --nocapture
#[test]
#[ignore]
fn f142_g() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::valley_construction::{GorgeFallKind, GorgeRetreat, LakeBase, gorge_plain, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone_protected};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let tan28 = 28f32.to_radians().tan();
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 142-G . spec v4 against its gates (témoin)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let (lm_in, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &vd, &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let sk_off = skeleton(&s1, &off, &ss, DOMAIN_KM);
    let links = trunk_links(&sk_off);
    let lb: Vec<bool> = links.iter().map(|l| lm_in[l.0] != 0 || dep[l.0] || lm_in[l.1] != 0 || dep[l.1]).collect();
    drop(sk_off);
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let map_a = |age: f32| if age <= 1.0 { (age - 0.7) / 0.3 } else { 1.0 + (age - 1.0) / 0.4 };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let msg = |e: Box<dyn std::any::Any + Send>| -> String {
        e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_else(|| "non-string panic".to_string())
    };
    // the breach's below-sea land cells of one world (the protected breach on the pre-breach drainage, as run_hd)
    let below_sea = |wd: &common::World| -> Vec<usize> {
        let g = &wd.heightmap;
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        (0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect()
    };
    // the extras (canyons, coast, θ, removal) for one world
    let extras = |label: &str, vc: ValleyConstruction, wd_h: &GridF32, cond: &GridF32, bodies_cells: &[bool]| -> String {
        let g = wd_h;
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(g, &vd, &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(g, cond, &pre, DELIVERED_P50_M, &ss, &vd, cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = 0;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons += 1;
            }
        }
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        let co = coast(g, &ss, &dwall, 2.0 / CELL_KM);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let eo = metres(g);
        let rem = |sel: &dyn Fn(usize) -> bool| -> f64 {
            (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((zs1[k] - eo[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum()
        };
        let dep_ = |sel: &dyn Fn(usize) -> bool| -> f64 {
            (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((eo[k] - zs1[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum()
        };
        let (rt, ri) = (rem(&|_| true), rem(&|k| bodies_cells[k]));
        let (dt, di) = (dep_(&|_| true), dep_(&|k| bodies_cells[k]));
        format!(
            "canyons **{canyons}** · spurs near a coastal wall {rn:.4} /km · removal (S1 − world, land) **{rt:.1} km³** (inside the kept bodies {ri:.1}, elsewhere {:.1}) · raise (world − S1, land) {dt:.2} km³ (inside the kept bodies {di:.2}) [{label}]",
            rt - ri
        )
    };
    let theta_carved = |zb: &[f32], carved: &[bool]| -> String {
        let zk: Vec<f32> = links.iter().map(|l| zb[l.0]).collect();
        let zr: Vec<f32> = links.iter().map(|l| zb[l.1]).collect();
        let cl: Vec<bool> = links.iter().map(|l| carved[l.0] && carved[l.1]).collect();
        let t = theta_links(&links, &zk, &zr, &|i| cl[i] && !lb[i]);
        format!("θ on the carved links, construction mask **{:.3}** [{:.3}, {:.3}] ({})", t.0, t.1, t.2, t.3)
    };
    // ── ON ×1: the 14 D8 lakes (G-bodies' reference) and the breach's below-sea set (G-sea's reference)
    let on_at = |age: f32| ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1)) };
    let on = on_at(1.0);
    eprintln!("\n────────── ON ×1 (the references and the controls) ──────────");
    let (d8, lake_on, sea_on, on_zb, on_carved, on_lakes_line, on_slope100) = {
        let wd = build_world(kn(on), None, PSEED, None);
        let sea_on = below_sea(&wd);
        let v = viz_hd_lakes_on(&wd, kn(on), PSEED, 45.0, 40.0);
        let dr = &v.drainage;
        let cn = metres(&v.conditioned);
        let lm = &dr.lake_map;
        let rcv = |c: usize| -> Option<usize> {
            let d = dr.flow.direction[c];
            if d == DIR_NONE {
                return None;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
        };
        let d8: Vec<u32> = dr
            .lakes
            .iter()
            .filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral) && rcv(l.base.outlet.1 as usize * w + l.base.outlet.0 as usize).is_some())
            .map(|l| l.base.id)
            .collect();
        // the G-slope100 negative control on ON
        let (mut nn, mut bad) = (0usize, 0usize);
        let mut lst = Vec::new();
        for l in dr.lakes.iter().filter(|l| d8.contains(&l.base.id)) {
            let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
            let Some(col) = rcv(o) else { continue };
            let (mut pth, mut dd) = (vec![col], vec![0f32]);
            let mut c = col;
            while *dd.last().unwrap() < 10.5 && cn[c] > 0.0 {
                let Some(r) = rcv(c) else { break };
                if lm[r] != 0 && lm[r] != l.base.id {
                    break;
                }
                dd.push(dd.last().unwrap() + CELL_KM * if (c % w != r % w) && (c / w != r / w) { std::f32::consts::SQRT_2 } else { 1.0 });
                pth.push(r);
                c = r;
            }
            let mut ms = 0f32;
            for i in 0..pth.len() {
                if dd[i] > 10.0 {
                    break;
                }
                if let Some(j) = (i..pth.len()).find(|&j| dd[j] >= dd[i] + 0.1) {
                    ms = ms.max((cn[pth[i]] - cn[pth[j]]) / ((dd[j] - dd[i]) * 1000.0));
                }
            }
            nn += 1;
            lst.push(ms.atan().to_degrees());
            if ms > tan28 {
                bad += 1;
            }
        }
        let sk = skeleton(&s1, &on, &ss, DOMAIN_KM);
        let (b, mk) = carve(&s1, &sk, &on, &ss);
        let line = format!(
            "{} lakes, non-crater {:.1} km²",
            dr.lakes.len(),
            dr.lakes.iter().filter(|l| !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).map(|l| l.area_km2).sum::<f32>()
        );
        (d8, lm.clone(), sea_on, metres(&b), mk.carved, line, format!("{nn} outlet paths (the 14 D8) · > 28° **{bad}** · max {:.1}°", q(&mut lst, 1.0)))
    };
    let sea_on_set: HashSet<usize> = sea_on.iter().copied().collect();
    eprintln!("   the 14 D8 lakes: {:?} · ON's below-sea land cells after the breach: {} · G-slope100 negative control: {on_slope100}", d8, sea_on.len());
    let mut d8_cells: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..n {
        if d8.contains(&lake_on[k]) {
            d8_cells.entry(lake_on[k]).or_default().push(k);
        }
    }
    // ── the GORGE v4 worlds
    let mut area_series: Vec<(f32, f64)> = Vec::new();
    let mut r18: Vec<String> = vec![format!("ON ×1: {on_lakes_line}")];
    let mut summary: Vec<String> = Vec::new();
    let runs: Vec<(f32, f32)> = vec![(1.0, 0.5), (0.7, 0.5), (0.85, 0.5), (1.2, 0.5), (1.4, 0.5), (1.0, 0.0), (1.0, 0.25), (1.4, 0.0)];
    for (age, p) in runs {
        let r = map_a(age);
        let label = format!("GORGE v4 ×{age} (r_world {r}, p {p})");
        eprintln!("\n────────── {label} ──────────");
        let t = Instant::now();
        let res = catch_unwind(AssertUnwindSafe(|| -> String {
            let vc = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(r, p)), ..on_at(age) };
            // the construction + the plain, and ON's construction at the same age (the cost, the gorge's removal)
            let tc = Instant::now();
            let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
            let (built0, mk) = carve(&s1, &sk, &vc, &ss);
            let mut built = built0.clone();
            let plain = gorge_plain(&mut built, &sk, &ss, 1.0);
            let t_con = tc.elapsed().as_secs_f64();
            let tn = Instant::now();
            let skn = skeleton(&s1, &on_at(age), &ss, DOMAIN_KM);
            let (bon, _) = carve(&s1, &skn, &on_at(age), &ss);
            let t_on = tn.elapsed().as_secs_f64();
            drop(skn);
            let zb0 = metres(&built0);
            let zb = metres(&built);
            let zon = metres(&bon);
            drop(bon);
            let bo = sk.gorge_body_of.clone().expect("the gate is on");
            let bodies = sk.gorge_bodies.clone();
            let rf = sk.rim_floor_m.clone().expect("the gate is on");
            let dir_sk = sk.direction.clone();
            let falls_sk = sk.gorge_falls.clone();
            let mut cells: Vec<Vec<usize>> = vec![Vec::new(); bodies.len()];
            for k in 0..n {
                if bo[k] != u32::MAX {
                    cells[bo[k] as usize].push(k);
                }
            }
            let mut bcells = vec![false; n];
            for k in 0..n {
                bcells[k] = bo[k] != u32::MAX;
            }
            eprintln!("   construction + plain {t_con:.1} s · ON's construction at ×{age} {t_on:.1} s · added **{:.1} s** · {} kept bodies · {} tagged falls in the skeleton", t_con - t_on, bodies.len(), falls_sk.len());
            // G-bodies
            let mut matched: BTreeMap<u32, (usize, f64)> = BTreeMap::new();
            let mut split = 0usize;
            for (&id, cs) in &d8_cells {
                let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
                for &c in cs {
                    if bo[c] != u32::MAX {
                        *ov.entry(bo[c]).or_insert(0) += 1;
                    }
                }
                if ov.len() > 1 {
                    split += 1;
                }
                if let Some((&b, &k)) = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) {
                    matched.insert(id, (b as usize, k as f64 / bodies[b as usize].cells as f64));
                }
            }
            let used: HashSet<usize> = matched.values().map(|v| v.0).collect();
            let ok11 = matched.len() == d8.len() && used.len() == d8.len() && matched.values().all(|v| v.1 > 0.5);
            let extra: Vec<String> = (0..bodies.len())
                .filter(|b| !used.contains(b))
                .map(|b| {
                    let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
                    for &c in &cells[b] {
                        if lake_on[c] != 0 {
                            *ov.entry(lake_on[c]).or_insert(0) += 1;
                        }
                    }
                    let best = ov.iter().max_by_key(|e| *e.1).map(|(id, k)| format!("ON lake {id} ({:.0} %)", 100.0 * *k as f64 / cells[b].len() as f64));
                    format!("body {b} ({:.1} km², r {:.2}) {}", cells[b].len() as f32 * cell_km2, bodies[b].r_lake, best.unwrap_or_else(|| "no ON lake".to_string()))
                })
                .collect();
            eprintln!(
                "   G-bodies · D8 lakes matched {} of {} · distinct bodies {} · body share > 50 % {} · lakes split over > 1 body {split} · **{}** · kept bodies matching no D8 lake: {} {:?}",
                matched.len(),
                d8.len(),
                used.len(),
                matched.values().filter(|v| v.1 > 0.5).count(),
                if ok11 && split == 0 { "1:1 HOLDS" } else { "FAILS" },
                extra.len(),
                extra
            );
            // G-deposit
            let dep_total: f64 = plain.iter().map(|p| p.volume_km3).sum();
            let rem_gorge: f64 = (0..n).filter(|&k| s1.data[k] > SEA).map(|k| ((zon[k] - zb0[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum();
            let rem_gorge_plain: f64 = (0..n).filter(|&k| s1.data[k] > SEA).map(|k| ((zon[k] - zb[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum();
            eprintln!(
                "   G-deposit · bodies with a plain (r_lake > 1): {} · deposited **{dep_total:.4} km³** · the gorge's removal against ON at ×{age} (construction, land) {rem_gorge:.3} km³ (after the plain {rem_gorge_plain:.3}) · deposit / removal **{:.1} %**",
                plain.len(),
                if rem_gorge > 0.0 { 100.0 * dep_total / rem_gorge } else { f64::NAN }
            );
            for pl in &plain {
                let id = matched.iter().find(|e| e.1.0 == pl.body as usize).map(|e| *e.0);
                eprintln!(
                    "      body {} (lake {:?}) r {:.2} · spill {:.1} m · raised {} cells ({:.2} km²) · deepest {:.1} m · {:.4} km³",
                    pl.body, id, pl.r_lake, pl.spill_m, pl.cells, pl.area_km2, pl.max_depth_m, pl.volume_km3
                );
            }
            // G-pits at a stage (AMENDED after the first run, declared): per body with r_lake > 1, a priority-flood fill
            // restricted to the footprint and seeded on its ring. The LAKE is the depression cells (fill > z + 0.1 m)
            // at the fill level of the body's lowest cell; a PIT is a depression cell at any other level (a closed
            // hollow outside the lake). The first run's measure took the lake as the cells at the ring's minimum and
            // counted a lake held by a sill inside the footprint as pits.
            let pits = |g: &GridF32| -> String {
                use std::cmp::Reverse;
                use std::collections::BinaryHeap;
                let mut tot = 0usize;
                let mut deep = 0f32;
                let mut per: Vec<String> = Vec::new();
                for (bi, b) in bodies.iter().enumerate() {
                    if b.r_lake <= 1.0 || cells[bi].is_empty() {
                        continue;
                    }
                    let me = bi as u32;
                    let mut fill: HashMap<usize, f32> = HashMap::with_capacity(cells[bi].len());
                    let mut heap: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
                    let mut ring_min = f32::INFINITY;
                    let mut seen: HashSet<usize> = HashSet::new();
                    for &c in &cells[bi] {
                        for k in 0..8 {
                            let m = nbt(c, k);
                            if bo[m] != me && seen.insert(m) {
                                ring_min = ring_min.min(g.data[m]);
                                heap.push(Reverse((g.data[m].to_bits(), m)));
                            }
                        }
                    }
                    while let Some(Reverse((lb, c))) = heap.pop() {
                        let lev = f32::from_bits(lb);
                        for k in 0..8 {
                            let m = nbt(c, k);
                            if bo[m] != me || fill.contains_key(&m) {
                                continue;
                            }
                            let l = lev.max(g.data[m]);
                            fill.insert(m, l);
                            heap.push(Reverse((l.to_bits(), m)));
                        }
                    }
                    let low = *cells[bi].iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("cells");
                    let lake = ((fill[&low] - g.data[low]) * n2m > 0.1).then(|| fill[&low]);
                    let (mut np, mut dp) = (0usize, 0f32);
                    let mut levels: HashSet<u32> = HashSet::new();
                    for &c in &cells[bi] {
                        let d = (fill[&c] - g.data[c]) * n2m;
                        if d > 0.1 && lake.is_none_or(|ll| ((fill[&c] - ll) * n2m).abs() > 0.01) {
                            np += 1;
                            dp = dp.max(d);
                            levels.insert(fill[&c].to_bits());
                        }
                    }
                    tot += np;
                    deep = deep.max(dp);
                    per.push(format!(
                        "b{bi} r {:.2} L(r) {:.1} · ring min {:.1} · lake {} · pits {np} ({} levels, {dp:.1} m)",
                        b.r_lake,
                        b.level,
                        c1_altitude_norm_to_metres(ring_min, &ss),
                        lake.map_or("none".to_string(), |ll| format!("{:.1} m", c1_altitude_norm_to_metres(ll, &ss))),
                        levels.len()
                    ));
                }
                format!("**{tot}** pit cells outside the lakes (deepest {deep:.1} m) · {}", per.join(" | "))
            };
            let any_plain = !plain.is_empty();
            let mut pit_lines = vec![format!("(bench) construction + plain: {}", pits(&built))];
            if any_plain {
                let b = kn(vc);
                for (st, k2) in [
                    ("(a) construction + plain + rims", Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..b }),
                    ("(b) + light pass", Knobs { erosion_off: true, bathymetry_off: true, ..b }),
                    ("(c) + droplets", Knobs { bathymetry_off: true, ..b }),
                ] {
                    let wst = build_world(k2, None, PSEED, None);
                    let mut s = format!("{st}: {}", pits(&wst.heightmap));
                    if st.starts_with("(a)") {
                        let dk: Vec<usize> = (0..n).filter(|&k| bcells[k] && wst.heightmap.data[k].to_bits() != built.data[k].to_bits()).collect();
                        let mx = dk.iter().map(|&k| ((wst.heightmap.data[k] - built.data[k]) * n2m).abs()).fold(0f32, f32::max);
                        let mut byb: BTreeMap<u32, usize> = BTreeMap::new();
                        for &k in &dk {
                            *byb.entry(bo[k]).or_insert(0) += 1;
                        }
                        s.push_str(&format!(" · the pipeline's (a) against the bench's construction + plain on the kept bodies: {} cells differ (max {mx:.2} m; by body {byb:?})", dk.len()));
                    }
                    pit_lines.push(s);
                }
            }
            let wd = build_world(kn(vc), None, PSEED, None);
            pit_lines.push(format!("(d) + bathymetry (the breach's input): {}", pits(&wd.heightmap)));
            for l in &pit_lines {
                eprintln!("   G-pits · {l}");
            }
            // G-sea
            let sea_w = below_sea(&wd);
            let sw: HashSet<usize> = sea_w.iter().copied().collect();
            let added = sw.difference(&sea_on_set).count();
            let missing = sea_on_set.difference(&sw).count();
            eprintln!(
                "   G-sea · land cells ≤ sea after the breach {} (ON {}) · added {added} · missing {missing} · **{}**",
                sea_w.len(),
                sea_on.len(),
                if added == 0 && missing == 0 { "EXACTLY ON's" } else { "DIFFERS" }
            );
            // the run_hd tail (F38's invariant), isolated
            let tail = catch_unwind(AssertUnwindSafe(|| viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0)));
            let v = match tail {
                Ok(v) => {
                    eprintln!("   F38 · the run_hd tail completes: the invariant does NOT fire");
                    v
                }
                Err(e) => {
                    let m = msg(e);
                    eprintln!("   F38 · **THE INVARIANT FIRES** in the run_hd tail (the HD assembly): {m}");
                    return format!("{label}: F38 FIRED · {}", m.chars().take(160).collect::<String>());
                }
            };
            let cn = metres(&v.conditioned);
            let eo = metres(&wd.heightmap);
            let lm = &v.drainage.lake_map;
            let lakes = &v.drainage.lakes;
            let lvl_of: HashMap<u32, (f32, f32)> = lakes.iter().map(|l| (l.base.id, (l.level_m, l.area_km2))).collect();
            let ring_of = |bi: usize| -> Vec<usize> {
                let mut r: Vec<usize> = Vec::new();
                for &c in &cells[bi] {
                    for k in 0..8 {
                        let m = nbt(c, k);
                        if bo[m] == u32::MAX && s1.data[m] > SEA {
                            r.push(m);
                        }
                    }
                }
                r.sort_unstable();
                r.dedup();
                r
            };
            let mut up_body = vec![u32::MAX - 1; n];
            let mut path = Vec::new();
            for s in 0..n {
                if up_body[s] != u32::MAX - 1 {
                    continue;
                }
                path.clear();
                let mut c = s;
                let v0;
                loop {
                    if up_body[c] != u32::MAX - 1 {
                        v0 = up_body[c];
                        break;
                    }
                    if bo[c] != u32::MAX && c != s {
                        v0 = bo[c];
                        break;
                    }
                    path.push(c);
                    let d = dir_sk[c];
                    if d == DIR_NONE || s1.data[c] <= SEA || path.len() > 50_000 {
                        v0 = u32::MAX;
                        break;
                    }
                    c = nbt(c, d as usize);
                }
                for &q2 in &path {
                    up_body[q2] = v0;
                }
            }
            let head_at: HashMap<u32, f32> = falls_sk.iter().filter(|f| f.kind == GorgeFallKind::WaterfallHead).map(|f| (f.exit_y * w as u32 + f.exit_x, f.height_m)).collect();
            let (mut nl, mut ok_lv, mut drained, mut drained_ok) = (0usize, 0usize, 0usize, 0usize);
            let mut bad_lv: Vec<String> = Vec::new();
            let (mut rim_n, mut rim_bad_b, mut rim_bad_w) = (0usize, 0usize, 0usize);
            let (mut ring_bound, mut ring_viol) = (0usize, 0usize);
            let mut ring_worst = 0f32;
            let mut ring_where: Vec<String> = Vec::new();
            let (mut s100_n, mut s100_bad) = (0usize, 0usize);
            let mut s100_list: Vec<f32> = Vec::new();
            let (mut drn_up, mut drn_bad) = (0usize, 0usize);
            let mut gorge_viol = 0usize;
            let mut area14 = 0f64;
            for (bi, b) in bodies.iter().enumerate() {
                if b.col == u32::MAX {
                    continue;
                }
                nl += 1;
                let mut ov: HashMap<u32, usize> = HashMap::new();
                for &c in &cells[bi] {
                    if lm[c] != 0 {
                        *ov.entry(lm[c]).or_insert(0) += 1;
                    }
                }
                let best = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))).map(|(&id, &c)| (id, c));
                let id14 = matched.iter().find(|e| e.1.0 == bi).map(|e| *e.0);
                if b.r_lake >= 1.999 {
                    drained += 1;
                    if best.is_none_or(|(_, c)| 2 * c < cells[bi].len()) {
                        drained_ok += 1;
                    } else if let Some((id, c)) = best {
                        bad_lv.push(format!("body {bi} (D8 lake {id14:?}) drained, yet final lake {id} covers {:.0} % at {:.1} m (L_floor {:.1})", 100.0 * c as f64 / cells[bi].len() as f64, lvl_of[&id].0, b.l_floor));
                    }
                } else if let Some((id, _)) = best {
                    let (lv, ar) = lvl_of[&id];
                    if (lv - b.level).abs() <= 1.0 {
                        ok_lv += 1;
                    } else {
                        bad_lv.push(format!("body {bi} (D8 lake {id14:?}) {:.1} km² r {:.2} L(r) {:.1} final {lv:.1} (A {:.0} km²)", cells[bi].len() as f32 * cell_km2, b.r_lake, b.level, b.a_out_km2));
                    }
                    if id14.is_some() {
                        area14 += ar as f64;
                    }
                } else {
                    bad_lv.push(format!("body {bi} (D8 lake {id14:?}) {:.1} km² r {:.2} L(r) {:.1}: NO final lake", cells[bi].len() as f32 * cell_km2, b.r_lake, b.level));
                }
                let ring = ring_of(bi);
                if b.r_lake <= 1.0 {
                    rim_n += 1;
                    let mb = ring.iter().map(|&c| zb[c]).fold(f32::INFINITY, f32::min);
                    let mw = ring.iter().map(|&c| eo[c]).fold(f32::INFINITY, f32::min);
                    if b.level - mb > 10.0 {
                        rim_bad_b += 1;
                    }
                    if b.level - mw > 10.0 {
                        rim_bad_w += 1;
                    }
                }
                for &c in &ring {
                    if rf[c].is_finite() && (zb[c] - rf[c]).abs() < 0.05 && zs1[c] > rf[c] + 0.05 {
                        ring_bound += 1;
                        if c == b.col as usize {
                            continue;
                        }
                        for k in 0..8 {
                            let m = nbt(c, k);
                            if bo[m] != u32::MAX || rf[m].is_finite() {
                                continue;
                            }
                            let dist = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                            let sl = (zb[c] - zb[m]) / dist;
                            if sl > tan28 {
                                ring_viol += 1;
                                if sl > ring_worst {
                                    ring_worst = sl;
                                }
                                if ring_where.len() < 6 {
                                    ring_where.push(format!("body {bi} ({},{}) {:.0}°", c % w, c / w, sl.atan().to_degrees()));
                                }
                                break;
                            }
                        }
                    }
                }
                if let Some((id, _)) = best
                    && b.r_lake < 1.999
                    && let Some(l) = lakes.iter().find(|l| l.base.id == id)
                {
                    let o = l.base.outlet.1 as usize * w + l.base.outlet.0 as usize;
                    let rcv = |c: usize| -> Option<usize> {
                        let d = v.drainage.flow.direction[c];
                        if d == DIR_NONE {
                            return None;
                        }
                        let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
                    };
                    if let Some(col) = rcv(o) {
                        let (mut pth, mut dd) = (vec![col], vec![0f32]);
                        let mut c = col;
                        while *dd.last().unwrap() < 10.5 && cn[c] > 0.0 {
                            let Some(r2) = rcv(c) else { break };
                            if lm[r2] != 0 && lm[r2] != id {
                                break;
                            }
                            dd.push(dd.last().unwrap() + CELL_KM * if (c % w != r2 % w) && (c / w != r2 / w) { std::f32::consts::SQRT_2 } else { 1.0 });
                            pth.push(r2);
                            c = r2;
                        }
                        let skip = if head_at.contains_key(&b.exit) { 2 } else { 0 };
                        let mut ms = 0f32;
                        for i in skip..pth.len() {
                            if dd[i] > 10.0 {
                                break;
                            }
                            if let Some(j) = (i..pth.len()).find(|&j| dd[j] >= dd[i] + 0.1) {
                                ms = ms.max((cn[pth[i]] - cn[pth[j]]) / ((dd[j] - dd[i]) * 1000.0));
                            }
                        }
                        s100_n += 1;
                        s100_list.push(ms.atan().to_degrees());
                        if ms > tan28 {
                            s100_bad += 1;
                        }
                    }
                }
                if b.r_lake >= 1.999 {
                    for k in 0..n {
                        if up_body[k] == bi as u32 {
                            drn_up += 1;
                            if zb[k] < b.l_floor - 1.0 {
                                drn_bad += 1;
                            }
                        }
                    }
                }
                let hf = (1.5 * b.slope * 100.0).max(10.0);
                let glen = if b.slope > 0.0 { (b.d_g / b.slope / 1000.0).min(20.0) } else { 0.0 };
                let mut gp = vec![b.col as usize];
                let mut sd = 0f32;
                let mut c = b.col as usize;
                while sd < glen {
                    let d = dir_sk[c];
                    if d == DIR_NONE {
                        break;
                    }
                    let r2 = nbt(c, d as usize);
                    if s1.data[r2] <= SEA || bo[r2] != u32::MAX {
                        break;
                    }
                    sd += CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                    gp.push(r2);
                    c = r2;
                }
                for i in 1..gp.len().saturating_sub(2) {
                    if zb[gp[i]] - zb[gp[i + 2]] > hf {
                        gorge_viol += 1;
                    }
                }
            }
            if p == 0.5 {
                area_series.push((age, area14));
            }
            eprintln!("   G-levels · kept bodies with an outflow {nl} · at L(r) ± 1 m **{ok_lv}** of {} not drained · drained (r_lake = 2) {drained}, absent at the end **{drained_ok}**", nl - drained);
            for l in bad_lv.iter().take(14) {
                eprintln!("      off: {l}");
            }
            eprintln!("   G-rim (r_lake ≤ 1) · {rim_n} bodies · cuts > 10 m below L(r) on the construction **{rim_bad_b}** · on the eroded world **{rim_bad_w}**");
            eprintln!("   G-ring · ring cells bound by the clamp **{ring_bound}** · with a slope > 28° to an outside neighbour **{ring_viol}** (worst {:.0}°) {:?}", ring_worst.atan().to_degrees(), ring_where);
            eprintln!(
                "   G-slope100 · {s100_n} outlet paths · max slope over 100 m > 28° **{s100_bad}** · p50 {:.1}° · max {:.1}°",
                q(&mut s100_list.clone(), 0.5),
                q(&mut s100_list, 1.0)
            );
            eprintln!("   G-drained · upstream cells of drained lakes {drn_up} · below L_floor − 1 m on the construction **{drn_bad}**");
            let mut fh: Vec<f32> = wd.gorge_falls.iter().filter(|f| f.kind == GorgeFallKind::WaterfallHead).map(|f| f.height_m).collect();
            let ns = wd.gorge_falls.iter().filter(|f| f.kind == GorgeFallKind::WaterfallShortage).count();
            eprintln!(
                "   G-tag · tagged falls: head **{}** (median {:.1} m, max {:.1} m) · shortage **{ns}** · gorge 2-cell drops over H_f beyond the lip **{gorge_viol}**",
                fh.len(),
                q(&mut fh.clone(), 0.5),
                q(&mut fh, 1.0)
            );
            let all_area: f64 = lakes.iter().filter(|l| !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).map(|l| l.area_km2 as f64).sum();
            eprintln!("   G-area · final non-crater lakes {all_area:.1} km² · the 14 D8 bodies' final lakes {area14:.1} km²");
            if (age == 1.0 && p == 0.5) || (age == 1.4 && p == 0.0) {
                let ex = extras(&label, vc, &wd.heightmap, &v.conditioned, &bcells);
                let th = theta_carved(&zb, &mk.carved);
                eprintln!("   extras · {ex} · {th}");
                r18.push(format!("{label}: {} lakes {all_area:.1} km² · deposit {dep_total:.4} km³ · {ex} · {th}", lakes.len()));
            }
            format!("{label}: completed · G-bodies {} · G-sea +{added}/−{missing} · deposit {dep_total:.4} km³", if ok11 && split == 0 { "1:1" } else { "FAILS" })
        }));
        let line = match res {
            Ok(s) => s,
            Err(e) => {
                let m = msg(e);
                eprintln!("   **THE WORLD PANICKED** (outside the run_hd tail): {m}");
                format!("{label}: PANICKED · {}", m.chars().take(160).collect::<String>())
            }
        };
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
        summary.push(line);
    }
    // ── the controls ON / OFF at ×1 (θ, extras)
    for (label, vc) in [("ON ×1", on), ("OFF ×1", off)] {
        eprintln!("\n────────── {label} (controls) ──────────");
        let wd = build_world(kn(vc), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0);
        let none = vec![false; n];
        let (zb, carved) = if label.starts_with("ON") {
            (on_zb.clone(), on_carved.clone())
        } else {
            let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
            let (b, mk) = carve(&s1, &sk, &vc, &ss);
            (metres(&b), mk.carved)
        };
        let ex = extras(label, vc, &wd.heightmap, &v.conditioned, &none);
        let th = theta_carved(&zb, &carved);
        eprintln!("   extras · {ex} · {th}");
        r18.push(format!(
            "{label}: {} lakes {:.1} km² · {ex} · {th}",
            v.drainage.lakes.len(),
            v.drainage.lakes.iter().filter(|l| !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral)).map(|l| l.area_km2).sum::<f32>()
        ));
    }
    area_series.sort_by(|a, b| a.0.total_cmp(&b.0));
    eprintln!("\n   G-area (p = 0.5) · the 14 D8 bodies' final lakes by age: {}", area_series.iter().map(|a| format!("×{}: {:.1} km²", a.0, a.1)).collect::<Vec<_>>().join(" · "));
    eprintln!("   G-area · non-increasing: **{}**", area_series.windows(2).all(|p| p[1].1 <= p[0].1 + 1e-6));
    eprintln!("\n   RULE 18 · {}", r18.join("\n   RULE 18 · "));
    eprintln!("\n   SUMMARY · {}", summary.join("\n   SUMMARY · "));
    eprintln!("\n==========  end Finding 142-G . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 142-G, attribution — θ on the carved links (the construction's mask) fell to 0.264 under GORGE v4 ×1,
/// p 0.5 (OFF 0.493). Which cells move it: the outlet paths the gorge RAISES (the invariant), or the catchments the
/// lakes' lowered base LOWERS (the rebase)? On the construction alone: ON ×1 against GORGE v4 ×1, with each class
/// masked out. Declared in `f142_declared.md` (amendment).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f142_theta --nocapture
#[test]
#[ignore]
fn f142_theta() {
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase, gorge_plain, ocean_flood};
    let ss = SteinSteinParams::default();
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let (lm_in, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let links = trunk_links(&skeleton(&s1, &off, &ss, DOMAIN_KM));
    let lb: Vec<bool> = links.iter().map(|l| lm_in[l.0] != 0 || dep[l.0] || lm_in[l.1] != 0 || dep[l.1]).collect();
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..off };
    let gz = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.5)), ..on };
    let sk_on = skeleton(&s1, &on, &ss, DOMAIN_KM);
    let (b_on, mk_on) = carve(&s1, &sk_on, &on, &ss);
    let sk_g = skeleton(&s1, &gz, &ss, DOMAIN_KM);
    let (mut b_g, mk_g) = carve(&s1, &sk_g, &gz, &ss);
    let _ = gorge_plain(&mut b_g, &sk_g, &ss, 1.0);
    let (z_on, z_g) = (metres(&b_on), metres(&b_g));
    // the classes, from the skeletons' bases: raised by the gorge, lowered by the rebase, unchanged
    let mut raised = vec![false; n];
    let mut lowered = vec![false; n];
    for k in 0..n {
        let (a, b) = (sk_on.base_alt_m[k], sk_g.base_alt_m[k]);
        if a.is_finite() && b.is_finite() {
            raised[k] = b > a + 0.5;
            lowered[k] = b < a - 0.5;
        }
    }
    eprintln!(
        "\n=== F142 θ attribution · cells whose base the gorge raises {} · lowers {} (of {} land)",
        raised.iter().filter(|&&x| x).count(),
        lowered.iter().filter(|&&x| x).count(),
        (0..n).filter(|&k| s1.data[k] > SEA).count()
    );
    let th = |zb: &[f32], carved: &[bool], extra: &dyn Fn(usize) -> bool| -> String {
        let zk: Vec<f32> = links.iter().map(|l| zb[l.0]).collect();
        let zr: Vec<f32> = links.iter().map(|l| zb[l.1]).collect();
        let t = theta_links(&links, &zk, &zr, &|i| carved[links[i].0] && carved[links[i].1] && !lb[i] && extra(i));
        format!("**{:.3}** [{:.3}, {:.3}] ({} links)", t.0, t.1, t.2, t.3)
    };
    let any = |i: usize, m: &[bool]| m[links[i].0] || m[links[i].1];
    eprintln!("   ON ×1, its own mask: {}", th(&z_on, &mk_on.carved, &|_| true));
    eprintln!("   GORGE v4 ×1, its own mask: {}", th(&z_g, &mk_g.carved, &|_| true));
    eprintln!("   GORGE v4 ×1, without the links the gorge RAISES: {}", th(&z_g, &mk_g.carved, &|i| !any(i, &raised)));
    eprintln!("   GORGE v4 ×1, without the links the rebase LOWERS: {}", th(&z_g, &mk_g.carved, &|i| !any(i, &lowered)));
    eprintln!("   GORGE v4 ×1, without both: {}", th(&z_g, &mk_g.carved, &|i| !any(i, &raised) && !any(i, &lowered)));
    eprintln!("   ON ×1 on the same links as the last line: {}", th(&z_on, &mk_on.carved, &|i| !any(i, &raised) && !any(i, &lowered)));
    eprintln!("   GORGE v4 ×1 on the raised links only: {} · on the lowered links only: {}", th(&z_g, &mk_g.carved, &|i| any(i, &raised)), th(&z_g, &mk_g.carved, &|i| any(i, &lowered)));
}

/// ADR Finding 142-G, attribution — G-sea failed at GORGE v4 ×1.4 (land cells ≤ sea after the breach beyond ON's 292).
/// Where are the added cells: in which kept body (or how far from one), ocean-connected or inland, covered by which
/// final lake, and how high they stood before the breach? Declared in `f142_declared.md` (amendment).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f142_sea --nocapture
#[test]
#[ignore]
fn f142_sea() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashSet};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    use ymir_core::terrain::flow::breach_monotone_protected;
    let ss = SteinSteinParams::default();
    let vd = viz_dcfg();
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let on_at = |age: f32| ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1)) };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let below = |wd: &common::World| -> (Vec<usize>, GridF32) {
        let g = &wd.heightmap;
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        ((0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect(), c)
    };
    let on_set: HashSet<usize> = {
        let wd = build_world(kn(on_at(1.0)), None, PSEED, None);
        below(&wd).0.into_iter().collect()
    };
    eprintln!("\n=== F142 G-sea attribution · ON ×1: {} cells", on_set.len());
    for (age, p) in [(1.4f32, 0.5f32), (1.4, 0.0)] {
        let r = if age <= 1.0 { (age - 0.7) / 0.3 } else { 1.0 + (age - 1.0) / 0.4 };
        let vc = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(r, p)), ..on_at(age) };
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let bo = sk.gorge_body_of.clone().unwrap();
        let bodies = sk.gorge_bodies.clone();
        drop(sk);
        let mut inb = vec![false; n];
        for k in 0..n {
            inb[k] = bo[k] != u32::MAX;
        }
        let dist = dist_from(&inb, w, h);
        let wd = build_world(kn(vc), None, PSEED, None);
        let (set, cond) = below(&wd);
        let added: Vec<usize> = set.iter().copied().filter(|k| !on_set.contains(k)).collect();
        let wc = water_class(&cond, SEA);
        let v = viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0);
        let lm = &v.drainage.lake_map;
        let mut by_body: BTreeMap<String, usize> = BTreeMap::new();
        let mut by_lake: BTreeMap<u32, usize> = BTreeMap::new();
        let (mut inland, mut covered) = (0usize, 0usize);
        let mut zpre: Vec<f32> = Vec::new();
        for &k in &added {
            let key = if bo[k] != u32::MAX {
                let b = &bodies[bo[k] as usize];
                format!("in body {} (r {:.2}, L_floor {:.1}, L(r) {:.1})", bo[k], b.r_lake, b.l_floor, b.level)
            } else {
                format!("outside, {}–{} cells from a body (cap {CAP})", (dist[k] as usize / 10) * 10, (dist[k] as usize / 10) * 10 + 9)
            };
            *by_body.entry(key).or_insert(0) += 1;
            if wc[k] == WATER_CLASS_INLAND {
                inland += 1;
            }
            if lm[k] != 0 {
                covered += 1;
                *by_lake.entry(lm[k]).or_insert(0) += 1;
            }
            zpre.push(c1_altitude_norm_to_metres(wd.heightmap.data[k], &ss));
        }
        zpre.sort_by(f32::total_cmp);
        let xs: Vec<usize> = added.iter().map(|k| k % w).collect();
        let ys: Vec<usize> = added.iter().map(|k| k / w).collect();
        eprintln!(
            "\n   GORGE v4 ×{age} p {p}: {} cells (added {}, missing {}) · inland {inland} · covered by a final lake {covered} {:?} · height before the breach p0/p50/p100 {:.1}/{:.1}/{:.1} m · bbox x {}–{} y {}–{}",
            set.len(),
            added.len(),
            on_set.iter().filter(|k| !set.contains(k)).count(),
            by_lake,
            zpre.first().copied().unwrap_or(f32::NAN),
            zpre.get(zpre.len() / 2).copied().unwrap_or(f32::NAN),
            zpre.last().copied().unwrap_or(f32::NAN),
            xs.iter().min().unwrap_or(&0),
            xs.iter().max().unwrap_or(&0),
            ys.iter().min().unwrap_or(&0),
            ys.iter().max().unwrap_or(&0)
        );
        for (k, c) in &by_body {
            eprintln!("      {k}: {c}");
        }
    }
}

/// ADR Finding 143-L — the exact scope (S), the localisations of F142's unlocated failures (G-ring, G-drained, G-tag,
/// G-sea at ×1.4 p 0.5), the rebased links' slope against their law (θ) and the head falls under caps (H). Declared in
/// `f143_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f143_l --nocapture
#[test]
#[ignore]
fn f143_l() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase, carve_diag, gorge_plain, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone_protected};

    // ── the instrumented copy of `breach_monotone_protected` (flow.rs), line for line, plus the records
    struct Pq(f32, usize);
    impl PartialEq for Pq {
        fn eq(&self, o: &Self) -> bool {
            self.0 == o.0 && self.1 == o.1
        }
    }
    impl Eq for Pq {}
    impl PartialOrd for Pq {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for Pq {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            o.0.partial_cmp(&self.0).unwrap_or(std::cmp::Ordering::Equal).then_with(|| o.1.cmp(&self.1))
        }
    }
    /// A ramp that lowered at least one land cell to ≤ sea.
    struct Ev {
        pit: usize,
        ci: usize,
        floor: f32,
        steps: u32,
        stop: u8, // 0 a base, 1 a cell already lower than the ramp, 2 the chain's end
        stop_cell: usize,
        crossed: u32,
        cross_step: u32,
    }
    struct Br {
        z: Vec<f32>,
        pit_of: Vec<u32>,
        step_of: Vec<u32>,
        backlink: Vec<usize>,
        evs: Vec<Ev>,
    }
    fn breach_instr(height: &GridF32, filled: &GridF32, lake_map: &[u32], sea_level: f32, w: usize, h: usize, protect: Option<&[bool]>) -> Br {
        use ymir_core::terrain::flow::{D8_DX, D8_DY};
        let n = w * h;
        let prot = |k: usize| protect.is_some_and(|p| p[k]);
        let mut z = height.data.clone();
        for k in 0..n {
            if lake_map[k] != 0 && !prot(k) {
                z[k] = filled.data[k];
            }
        }
        let is_base = |k: usize, z: &[f32]| z[k] <= sea_level || lake_map[k] != 0;
        const EPS: f32 = 1e-5;
        let mut visited = vec![false; n];
        let mut backlink = vec![usize::MAX; n];
        let mut pit_of = vec![u32::MAX; n];
        let mut step_of = vec![0u32; n];
        let mut evs: Vec<Ev> = Vec::new();
        let mut heap: BinaryHeap<Pq> = BinaryHeap::new();
        for k in 0..n {
            if prot(k) {
                visited[k] = true;
            }
        }
        for k in 0..n {
            if !visited[k] && is_base(k, &z) {
                visited[k] = true;
                heap.push(Pq(z[k], k));
            }
        }
        while let Some(Pq(_, ci)) = heap.pop() {
            let (cx, cy) = (ci % w, ci / w);
            for d in 0..8 {
                let nx = ((cx as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                let ny = ((cy as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                let nb = ny * w + nx;
                if visited[nb] {
                    continue;
                }
                visited[nb] = true;
                backlink[nb] = ci;
                if height.data[nb] < z[ci] {
                    let mut target = height.data[nb] - EPS;
                    let mut cur = ci;
                    let (mut steps, mut crossed, mut cross_step) = (0u32, 0u32, u32::MAX);
                    while cur != usize::MAX && !is_base(cur, &z) && z[cur] > target {
                        if target <= sea_level {
                            crossed += 1;
                            if cross_step == u32::MAX {
                                cross_step = steps;
                            }
                        }
                        z[cur] = target;
                        pit_of[cur] = nb as u32;
                        step_of[cur] = steps;
                        target -= EPS;
                        cur = backlink[cur];
                        steps += 1;
                    }
                    if crossed > 0 {
                        let stop = if cur == usize::MAX {
                            2
                        } else if is_base(cur, &z) {
                            0
                        } else {
                            1
                        };
                        evs.push(Ev { pit: nb, ci, floor: height.data[nb], steps, stop, stop_cell: cur, crossed, cross_step });
                    }
                }
                z[nb] = height.data[nb];
                heap.push(Pq(z[nb], nb));
            }
        }
        let mut visited2 = vec![false; n];
        let mut heap2: BinaryHeap<Pq> = BinaryHeap::new();
        for k in 0..n {
            if prot(k) {
                visited2[k] = true;
            }
        }
        for k in 0..n {
            if !visited2[k] && is_base(k, &z) {
                visited2[k] = true;
                heap2.push(Pq(z[k], k));
            }
        }
        while let Some(Pq(_, ci)) = heap2.pop() {
            let (cx, cy) = (ci % w, ci / w);
            for d in 0..8 {
                let nx = ((cx as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                let ny = ((cy as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                let nb = ny * w + nx;
                if visited2[nb] {
                    continue;
                }
                visited2[nb] = true;
                if z[nb] < z[ci] {
                    z[nb] = z[ci];
                }
                heap2.push(Pq(z[nb], nb));
            }
        }
        Br { z, pit_of, step_of, backlink, evs }
    }


    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let tan28 = 28f32.to_radians().tan();
    let vd = viz_dcfg();
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    eprintln!("\n==========  Finding 143-L . the exact scope, the localisations, the rebased links, the head falls  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let on_at = |age: f32| ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1)) };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let map_a = |age: f32| if age <= 1.0 { (age - 0.7) / 0.3 } else { 1.0 + (age - 1.0) / 0.4 };
    let on = on_at(1.0);
    // ── ON ×1: the 14 D8 lakes, ON's breach set, ON's construction
    let (d8, lake_on, sea_on, zd_on) = {
        let wd = build_world(kn(on), None, PSEED, None);
        let g = &wd.heightmap;
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        let sea_on: HashSet<usize> = (0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect();
        let v = viz_hd_lakes_on(&wd, kn(on), PSEED, 45.0, 40.0);
        let dr = &v.drainage;
        let rcv = |c: usize| -> Option<usize> {
            let dd = dr.flow.direction[c];
            if dd == DIR_NONE {
                return None;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[dd as usize], (c / w) as i32 + D8_DY[dd as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
        };
        let d8: Vec<u32> = dr
            .lakes
            .iter()
            .filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral) && rcv(l.base.outlet.1 as usize * w + l.base.outlet.0 as usize).is_some())
            .map(|l| l.base.id)
            .collect();
        (d8, dr.lake_map.clone(), sea_on, metres(g))
    };
    let sk_on = skeleton(&s1, &on, &ss, DOMAIN_KM);
    let (b_on, _) = carve(&s1, &sk_on, &on, &ss);
    let zb_on = metres(&b_on);
    drop(b_on);
    // ── S: the exact scope
    let proxy = skeleton(&s1, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.5)), ..on }, &ss, DOMAIN_KM);
    let bo_p = proxy.gorge_body_of.clone().unwrap();
    let mut pcells: Vec<Vec<usize>> = vec![Vec::new(); proxy.gorge_bodies.len()];
    for k in 0..n {
        if bo_p[k] != u32::MAX {
            pcells[bo_p[k] as usize].push(k);
        }
    }
    let mut keep: BTreeMap<u32, usize> = BTreeMap::new();
    for &id in &d8 {
        let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
        for k in 0..n {
            if lake_on[k] == id && bo_p[k] != u32::MAX {
                *ov.entry(bo_p[k]).or_insert(0) += 1;
            }
        }
        if let Some((&b, _)) = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) {
            keep.insert(id, b as usize);
        }
    }
    let mut list = [u32::MAX; 16];
    let mut lows: Vec<u32> = keep.values().map(|&b| proxy.gorge_bodies[b].low).collect();
    lows.sort_unstable();
    for (i, l) in lows.iter().enumerate() {
        list[i] = *l;
    }
    eprintln!("\n   S · the 14 D8 lakes {:?} → {} bodies · the exact list (lowest cells): {:?}", d8, keep.len(), lows);
    let kept: HashSet<usize> = keep.values().copied().collect();
    for (b, gb) in proxy.gorge_bodies.iter().enumerate() {
        if kept.contains(&b) {
            continue;
        }
        let ring: Vec<usize> = {
            let mut r: Vec<usize> = pcells[b].iter().flat_map(|&c| (0..8).map(move |k| (c, k))).map(|(c, k)| nbt(c, k)).filter(|&m| bo_p[m] != b as u32).collect();
            r.sort_unstable();
            r.dedup();
            r
        };
        let rmin = |z: &[f32]| ring.iter().map(|&c| z[c]).fold(f32::INFINITY, f32::min);
        let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
        for &c in &pcells[b] {
            if lake_on[c] != 0 {
                *ov.entry(lake_on[c]).or_insert(0) += 1;
            }
        }
        eprintln!(
            "      DROPPED body {b} at {} · {:.1} km² · input lake {} · ring touches a depression {} · the sea {} (why the proxy kept it: a lone input lake) · L_in {:.1} · A {:.1} km² · in ON: the ring's minimum on the construction {:.1} m ({:+.1} vs L_in), on the world (d) {:.1} m ({:+.1}) · ON's final lakes over it {:?}",
            xy(gb.low as usize),
            gb.cells as f32 * cell_km2,
            gb.input_lake,
            gb.touches_depression,
            gb.touches_below_sea,
            gb.l_in,
            gb.a_out_km2,
            rmin(&zb_on),
            rmin(&zb_on) - gb.l_in,
            rmin(&zd_on),
            rmin(&zd_on) - gb.l_in,
            ov
        );
    }
    let low_to_id: HashMap<u32, u32> = keep.iter().map(|(&id, &b)| (proxy.gorge_bodies[b].low, id)).collect();
    drop(proxy);
    let ex = |r: f32, p: f32| GorgeRetreat { scope_lows: Some(list), ..GorgeRetreat::v4(r, p) };
    // ── the exact GORGE ×1, p 0.5: the construction with carve_diag
    let vc1 = ValleyConstruction { gorge_retreat: Some(ex(1.0, 0.5)), ..on };
    let sk1 = skeleton(&s1, &vc1, &ss, DOMAIN_KM);
    let (c1o, mk1, dg1) = carve_diag(&s1, &sk1, &vc1, &ss);
    let mut b1 = c1o.clone();
    drop(c1o);
    let _ = gorge_plain(&mut b1, &sk1, &ss, 1.0);
    let zb1 = metres(&b1);
    let bo1 = sk1.gorge_body_of.clone().unwrap();
    let bodies1 = sk1.gorge_bodies.clone();
    let rf1 = sk1.rim_floor_m.clone().unwrap();
    eprintln!("   S · the exact scope keeps {} bodies (G-bodies: exactly 14 in 1:1 is checked in f143_c)", bodies1.len());
    let mut cells1: Vec<Vec<usize>> = vec![Vec::new(); bodies1.len()];
    for k in 0..n {
        if bo1[k] != u32::MAX {
            cells1[bo1[k] as usize].push(k);
        }
    }
    let line_at = |c: usize| -> Option<u32> {
        let s = dg1.banded.as_ref().map(|b| b[c]).filter(|&s| s != u32::MAX).unwrap_or(dg1.who[c]);
        (s != u32::MAX).then(|| dg1.line_of[s as usize])
    };
    // the (b) world: the light pass
    let eo_b = {
        let wb = build_world(Knobs { erosion_off: true, bathymetry_off: true, ..kn(vc1) }, None, PSEED, None);
        metres(&wb.heightmap)
    };
    // ── L-ring
    {
        let mut cat: BTreeMap<String, usize> = BTreeMap::new();
        let mut ex_: Vec<String> = Vec::new();
        let (mut nv_c, mut nv_b) = (0usize, 0usize);
        for (bi, b) in bodies1.iter().enumerate() {
            if b.col == u32::MAX {
                continue;
            }
            let lcol = line_at(b.col as usize);
            let mut ring: Vec<usize> = cells1[bi].iter().flat_map(|&c| (0..8).map(move |k| (c, k))).map(|(c, k)| nbt(c, k)).filter(|&m| bo1[m] == u32::MAX && s1.data[m] > SEA).collect();
            ring.sort_unstable();
            ring.dedup();
            for &c in &ring {
                if !(rf1[c].is_finite() && (zb1[c] - rf1[c]).abs() < 0.05 && zs1[c] > rf1[c] + 0.05) || c == b.col as usize {
                    continue;
                }
                let mut hit_c = false;
                let mut hit_b = false;
                for k in 0..8 {
                    let m = nbt(c, k);
                    if bo1[m] != u32::MAX || rf1[m].is_finite() {
                        continue;
                    }
                    let dist = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                    if !hit_c && (zb1[c] - zb1[m]) / dist > tan28 {
                        hit_c = true;
                        let cut = zs1[m] > zb1[m] + 0.5;
                        let lm_ = line_at(m);
                        let who = match (lm_, lcol) {
                            (Some(a), Some(bb)) if a == bb => "the OUTLET line",
                            (Some(_), _) => "ANOTHER line",
                            (None, _) => "no sample",
                        };
                        let key = format!("neighbour {} · laid by {who}", if cut { "cut by the construction" } else { "uncut (S1's own slope)" });
                        *cat.entry(key).or_insert(0) += 1;
                        if ex_.len() < 8 {
                            ex_.push(format!("body {bi} ring {} rf {:.1} ↔ {} S1 {:.1} → {:.1} ({who})", xy(c), rf1[c], xy(m), zs1[m], zb1[m]));
                        }
                    }
                    if !hit_b && (eo_b[c] - eo_b[m]) / dist > tan28 {
                        hit_b = true;
                    }
                }
                if hit_c {
                    nv_c += 1;
                }
                if hit_b {
                    nv_b += 1;
                }
            }
        }
        eprintln!("\n   L-ring (×1, p 0.5) · clamped ring cells with a > 28° drop to an outside neighbour: on the construction **{nv_c}** · the same cells on the light pass's world **{nv_b}**");
        for (k, v) in &cat {
            eprintln!("      {k}: {v}");
        }
        for e in &ex_ {
            eprintln!("      e.g. {e}");
        }
    }
    // ── L-tag
    {
        let dir = &sk1.direction;
        let mut cat: BTreeMap<String, usize> = BTreeMap::new();
        let mut ex_: Vec<String> = Vec::new();
        let mut per: BTreeMap<usize, usize> = BTreeMap::new();
        for (bi, b) in bodies1.iter().enumerate() {
            if b.col == u32::MAX {
                continue;
            }
            let lcol = line_at(b.col as usize);
            let hf = (1.5 * b.slope * 100.0).max(10.0);
            let glen = if b.slope > 0.0 { (b.d_g / b.slope / 1000.0).min(20.0) } else { 0.0 };
            let (mut gp, mut sd) = (vec![b.col as usize], vec![0f32]);
            let mut c = b.col as usize;
            while *sd.last().unwrap() < glen {
                let d = dir[c];
                if d == DIR_NONE {
                    break;
                }
                let r2 = nbt(c, d as usize);
                if s1.data[r2] <= SEA || bo1[r2] != u32::MAX {
                    break;
                }
                sd.push(sd.last().unwrap() + CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 });
                gp.push(r2);
                c = r2;
            }
            for i in 1..gp.len().saturating_sub(2) {
                let dz = zb1[gp[i]] - zb1[gp[i + 2]];
                if dz > hf {
                    *per.entry(bi).or_insert(0) += 1;
                    let fl = (mk1.floor[gp[i]], mk1.floor[gp[i + 1]], mk1.floor[gp[i + 2]]);
                    let lines: Vec<&str> = [gp[i], gp[i + 1], gp[i + 2]].iter().map(|&q| match (line_at(q), lcol) {
                        (Some(a), Some(bb)) if a == bb => "outlet",
                        (Some(_), _) => "other",
                        _ => "none",
                    }).collect();
                    let key = format!("floor flags {:?} · lines {:?}", fl, lines);
                    *cat.entry(key).or_insert(0) += 1;
                    if ex_.len() < 10 {
                        ex_.push(format!("body {bi} at {:.2} km from the col {} : {:.1} → {:.1} → {:.1} m (H_f {hf:.1})", sd[i], xy(gp[i]), zb1[gp[i]], zb1[gp[i + 1]], zb1[gp[i + 2]]));
                    }
                }
            }
        }
        eprintln!("\n   L-tag (×1, p 0.5, the construction) · 2-cell drops > H_f beyond the lip, by body {:?}", per);
        for (k, v) in &cat {
            eprintln!("      {k}: {v}");
        }
        for e in &ex_ {
            eprintln!("      e.g. {e}");
        }
    }
    // ── H: the head falls (×1, p 0.5) and the lip drop in C0
    {
        let caps: [(&str, Option<f32>); 4] = [("50 m", Some(50.0)), ("100 m", Some(100.0)), ("200 m", Some(200.0)), ("none", None)];
        eprintln!("\n   H (×1, p 0.5, the 14) · A, D_g, φ, the head fall, S, H_f · the lip drop in C0 (the ring's minimum on the light pass's world − L(r)):");
        let mut rows: Vec<(f32, f32, f32)> = Vec::new(); // (A, head, H_f)
        for (bi, b) in bodies1.iter().enumerate() {
            let mut ring: Vec<usize> = cells1[bi].iter().flat_map(|&c| (0..8).map(move |k| (c, k))).map(|(c, k)| nbt(c, k)).filter(|&m| bo1[m] == u32::MAX).collect();
            ring.sort_unstable();
            ring.dedup();
            let rb = ring.iter().map(|&c| eo_b[c]).fold(f32::INFINITY, f32::min);
            let hf = (1.5 * b.slope * 100.0).max(10.0);
            eprintln!(
                "      body {bi:>2} (lake {:?}, low {}) · A {:>7.1} km² · D_g {:>6.1} m · φ {:.2} · head **{:>6.1} m** · S {:.2}° · H_f {hf:.1} · r_lake {:.2} · lip drop in C0 **{:+.1} m**",
                low_to_id.get(&b.low),
                xy(b.low as usize),
                b.a_out_km2,
                b.d_g,
                b.phi,
                b.head_fall_m,
                b.slope.atan().to_degrees(),
                b.r_lake,
                rb - b.level
            );
            rows.push((b.a_out_km2, b.head_fall_m, hf));
        }
        let hcap_a = |a: f32| 100.0 * (a / 418.7).max(1e-6).powf(-0.5);
        let mut lines = Vec::new();
        for (name, cap) in caps.iter().map(|&(nm, c)| (nm.to_string(), c)).chain(std::iter::once(("H_cap(A) = 100 m·(A/418.7)^−0.5".to_string(), None::<f32>))) {
            let is_a = name.starts_with("H_cap");
            let mut tag: Vec<f32> = rows
                .iter()
                .map(|&(a, hd, hf)| {
                    let c = if is_a { hcap_a(a) } else { cap.unwrap_or(f32::INFINITY) };
                    (hd.min(c), hf)
                })
                .filter(|&(hd, hf)| hd > hf)
                .map(|(hd, _)| hd)
                .collect();
            tag.sort_by(f32::total_cmp);
            lines.push(format!(
                "cap {name}: tagged **{}** · median {:.1} m · max {:.1} m",
                tag.len(),
                if tag.is_empty() { f32::NAN } else { tag[tag.len() / 2] },
                tag.last().copied().unwrap_or(f32::NAN)
            ));
        }
        for l in lines {
            eprintln!("      {l}");
        }
        eprintln!(
            "      H_cap(A) per body: {}",
            rows.iter().map(|&(a, hd, _)| format!("A {a:.0} → cap {:.0} m (head {hd:.0})", hcap_a(a))).collect::<Vec<_>>().join(" · ")
        );
        eprintln!("      lakes with a head fall > 100 m today: {}", rows.iter().filter(|r| r.1 > 100.0).count());
    }
    // ── θ: the rebased links
    {
        let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
        let links = trunk_links(&skeleton(&s1, &off, &ss, DOMAIN_KM));
        let (lm_in, dep) = {
            let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
            let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
            let spill = ocean_flood(&bf);
            let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
            let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
            (d.lake_map, dep)
        };
        let inb: Vec<bool> = (0..n).map(|k| bo1[k] != u32::MAX).collect();
        let dist = dist_from(&inb, w, h);
        let lowered = |k: usize| sk_on.base_alt_m[k].is_finite() && sk1.base_alt_m[k].is_finite() && sk1.base_alt_m[k] < sk_on.base_alt_m[k] - 0.5;
        let bins = [(0.0f32, 1.0f32), (1.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 1e9)];
        let ratio = |z: &[f32], l: &(usize, usize, f64, f32)| -> Option<f32> {
            let (a, b) = (l.0, l.1);
            let len = CELL_KM * 1000.0 * if (a % w != b % w) && (a / w != b / w) { std::f32::consts::SQRT_2 } else { 1.0 };
            let s = (z[a] - z[b]) / len;
            let sl = vc1.age_k * sk1.area_km2[a].max(vc1.a_c_km2).powf(-0.5);
            (sl > 0.0).then(|| s / sl)
        };
        eprintln!("\n   θ · the rebased links (base > 0.5 m under ON's) against the unchanged: S / S_law by distance upstream of the nearest kept body");
        for (label, sel) in [("REBASED", true), ("unchanged", false)] {
            for (zlabel, z) in [("construction", &zb1), ("light pass (b)", &eo_b), ("ON construction", &zb_on)] {
                let mut row = Vec::new();
                for &(lo, hi) in &bins {
                    let mut v: Vec<f32> = links
                        .iter()
                        .enumerate()
                        .filter(|(i, l)| mk1.carved[l.0] && mk1.carved[l.1] && !(lm_in[l.0] != 0 || dep[l.0] || lm_in[l.1] != 0 || dep[l.1]) && (lowered(l.0) || lowered(l.1)) == sel && {
                            let dk = dist[l.0] as f32 * CELL_KM;
                            dk >= lo && dk < hi
                        } && *i < usize::MAX)
                        .filter_map(|(_, l)| ratio(z, l))
                        .collect();
                    v.sort_by(f32::total_cmp);
                    let out = v.iter().filter(|&&r| !(0.5..=2.0).contains(&r)).count();
                    row.push(format!(
                        "{lo:.0}–{} km: n {} · median **{:.2}** · outside [0.5, 2] {:.0} %",
                        if hi > 1e8 { "∞".to_string() } else { format!("{hi:.0}") },
                        v.len(),
                        if v.is_empty() { f32::NAN } else { v[v.len() / 2] },
                        if v.is_empty() { 0.0 } else { 100.0 * out as f64 / v.len() as f64 }
                    ));
                }
                eprintln!("      {label} on the {zlabel}: {}", row.join(" · "));
            }
        }
    }
    drop(sk_on);
    // ── L-drained (×1.2 p 0.5 and ×1.4 p 0, the construction)
    for (age, p) in [(1.2f32, 0.5f32), (1.4, 0.0)] {
        let vc = ValleyConstruction { gorge_retreat: Some(ex(map_a(age), p)), ..on_at(age) };
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (co, _, dg) = carve_diag(&s1, &sk, &vc, &ss);
        let zc = metres(&co);
        let mut bp = co.clone();
        drop(co);
        let _ = gorge_plain(&mut bp, &sk, &ss, 1.0);
        let zp = metres(&bp);
        drop(bp);
        let bo = sk.gorge_body_of.clone().unwrap();
        let bodies = sk.gorge_bodies.clone();
        let dir = sk.direction.clone();
        let mut up_body = vec![u32::MAX - 1; n];
        let mut path = Vec::new();
        for s in 0..n {
            if up_body[s] != u32::MAX - 1 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if up_body[c] != u32::MAX - 1 {
                    v0 = up_body[c];
                    break;
                }
                if bo[c] != u32::MAX && c != s {
                    v0 = bo[c];
                    break;
                }
                path.push(c);
                let d = dir[c];
                if d == DIR_NONE || s1.data[c] <= SEA || path.len() > 50_000 {
                    v0 = u32::MAX;
                    break;
                }
                c = nbt(c, d as usize);
            }
            for &q in &path {
                up_body[q] = v0;
            }
        }
        let mut cat: BTreeMap<String, usize> = BTreeMap::new();
        let mut ex_: Vec<String> = Vec::new();
        let (mut bad_c, mut bad_p) = (0usize, 0usize);
        for k in 0..n {
            let ub = up_body[k];
            if ub == u32::MAX || ub == u32::MAX - 1 || bodies[ub as usize].r_lake < 1.999 {
                continue;
            }
            let lf = bodies[ub as usize].l_floor;
            if zc[k] < lf - 1.0 {
                bad_c += 1;
            }
            if zp[k] < lf - 1.0 {
                bad_p += 1;
                let s = dg.who[k];
                let key = if s == u32::MAX {
                    "no sample (S1 itself below the floor)".to_string()
                } else {
                    let (sx, sy, zf, _) = sk.polylines[dg.line_of[s as usize] as usize][dg.pos_of[s as usize] as usize];
                    let sc = (sy.floor() as i64).rem_euclid(h as i64) as usize * w + (sx.floor() as i64).rem_euclid(w as i64) as usize;
                    let sb = bo[sc];
                    let ub2 = up_body[sc];
                    let where_ = if sb == ub {
                        "in the SAME drained body".to_string()
                    } else if sb != u32::MAX {
                        format!("in ANOTHER body (L(r) {:.1})", bodies[sb as usize].level)
                    } else if ub2 == ub {
                        "upstream of the same body".to_string()
                    } else {
                        "on a line draining ELSEWHERE".to_string()
                    };
                    if ex_.len() < 8 {
                        ex_.push(format!("cell {} z {:.1} < L_floor {lf:.1} of body {ub} · sample {} zf {zf:.1} · {where_}", xy(k), zp[k], xy(sc)));
                    }
                    format!("laid by a sample {where_} · zf {} the floor", if zf < lf - 1.0 { "UNDER" } else { "above" })
                };
                *cat.entry(key).or_insert(0) += 1;
            }
        }
        eprintln!("\n   L-drained (×{age}, p {p}) · upstream of drained bodies, under L_floor − 1 m: at `carve` **{bad_c}** · after the plain **{bad_p}**");
        for (kk, v) in &cat {
            eprintln!("      {kk}: {v}");
        }
        for e in &ex_ {
            eprintln!("      e.g. {e}");
        }
    }
    // ── L-sea (×1.4, p 0.5): the instrumented breach
    {
        let vc = ValleyConstruction { gorge_retreat: Some(ex(2.0, 0.5)), ..on_at(1.4) };
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (co, _) = carve(&s1, &sk, &vc, &ss);
        let mut bp = co;
        let _ = gorge_plain(&mut bp, &sk, &ss, 1.0);
        let zcon = metres(&bp);
        drop(bp);
        let bo = sk.gorge_body_of.clone().unwrap();
        drop(sk);
        let wd = build_world(kn(vc), None, PSEED, None);
        let g = &wd.heightmap;
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let prod = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        let br = breach_instr(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        let diff = (0..n).filter(|&k| prod.data[k].to_bits() != br.z[k].to_bits()).count();
        eprintln!("\n   L-sea (×1.4, p 0.5) · the instrumented breach against production: {diff} cells differ{}", if diff == 0 { " (bit-identical)" } else { " — NOT READ" });
        if diff == 0 {
            let set: Vec<usize> = (0..n).filter(|&k| g.data[k] > SEA && br.z[k] <= SEA).collect();
            let added: Vec<usize> = set.iter().copied().filter(|k| !sea_on.contains(k)).collect();
            let sea0: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
            let dsea = dist_from(&sea0, w, h);
            let wc = water_class(&GridF32 { width: w, height: h, data: br.z.clone() }, SEA);
            let mut by_pit: BTreeMap<u32, usize> = BTreeMap::new();
            for &k in &added {
                *by_pit.entry(br.pit_of[k]).or_insert(0) += 1;
            }
            let mut dd: Vec<f32> = added.iter().map(|&k| dsea[k] as f32 * CELL_KM).collect();
            dd.sort_by(f32::total_cmp);
            let mut zc: Vec<f32> = added.iter().map(|&k| zcon[k]).collect();
            zc.sort_by(f32::total_cmp);
            eprintln!(
                "      added {} (ocean-connected {}) · their height on the construction p0/p50/p100 {:.1}/{:.1}/{:.1} m · their distance from the coast of (d): p50 {:.2} km · max {:.2} km · pits {}",
                added.len(),
                added.iter().filter(|&&k| wc[k] == WATER_CLASS_OCEAN).count(),
                zc.first().copied().unwrap_or(f32::NAN),
                zc.get(zc.len() / 2).copied().unwrap_or(f32::NAN),
                zc.last().copied().unwrap_or(f32::NAN),
                dd.get(dd.len() / 2).copied().unwrap_or(f32::NAN),
                dd.last().copied().unwrap_or(f32::NAN),
                by_pit.len()
            );
            for (&pit, &cnt) in by_pit.iter() {
                if pit == u32::MAX {
                    eprintln!("      {cnt} cells laid by no ramp");
                    continue;
                }
                let pk = pit as usize;
                let e = br.evs.iter().find(|e| e.pit == pk);
                eprintln!(
                    "      pit {} floor {:.1} m · spill {:.1} m · body {} · the construction there {:.1} m · {} cells · ramp {} steps, crossed the sea at step {} · base {}",
                    xy(pk),
                    c1_altitude_norm_to_metres(g.data[pk], &ss),
                    c1_altitude_norm_to_metres(d.flow.filled.data[pk], &ss),
                    if bo[pk] == u32::MAX { "none".to_string() } else { bo[pk].to_string() },
                    zcon[pk],
                    cnt,
                    e.map_or(0, |e| e.steps),
                    e.map_or(0, |e| e.cross_step),
                    e.map_or("—".to_string(), |e| if e.stop_cell == usize::MAX { "chain end".to_string() } else if g.data[e.stop_cell] <= SEA { format!("the sea at {}", xy(e.stop_cell)) } else { format!("{} at {}", if d.lake_map[e.stop_cell] != 0 { "a lake" } else { "a cell lowered earlier" }, xy(e.stop_cell)) })
                );
            }
        }
    }
    eprintln!("\n==========  end Finding 143-L . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 143-C — the decisive controls (rule 18) on F142's failures, with the exact scope: at ×1 (p 0.5),
/// ×1.4 (p 0.5) and ×1.4 (p 0), C0 (v4 + the exact scope), C1 (φ = 0), C2 (the light pass does not modify the design)
/// and C3 (both). Every gate on the construction (C0, C1) and on the final world (all twelve). Each world runs isolated.
/// Declared in `f143_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f143_c --nocapture
#[test]
#[ignore]
fn f143_c() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::{C1_SEA_LEVEL_NORM, LakeType};
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase, gorge_plain, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone_protected};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let tan28 = 28f32.to_radians().tan();
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 143-C . the controls C0–C3 (exact scope)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let (lm_in, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &vd, &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let links = trunk_links(&skeleton(&s1, &off, &ss, DOMAIN_KM));
    let lb: Vec<bool> = links.iter().map(|l| lm_in[l.0] != 0 || dep[l.0] || lm_in[l.1] != 0 || dep[l.1]).collect();
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let map_a = |age: f32| if age <= 1.0 { (age - 0.7) / 0.3 } else { 1.0 + (age - 1.0) / 0.4 };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let msg = |e: Box<dyn std::any::Any + Send>| -> String {
        e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_else(|| "non-string panic".to_string())
    };
    let breach_of = |g: &GridF32, wd: Option<&common::World>| -> (Vec<usize>, GridF32, ymir_core::tectonics_c1::drainage::C1DrainageResult) {
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let prot = wd.and_then(|wd| wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h)));
        let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        ((0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect(), c, d)
    };
    let theta_on = |zb: &[f32], carved: &[bool]| -> String {
        let zk: Vec<f32> = links.iter().map(|l| zb[l.0]).collect();
        let zr: Vec<f32> = links.iter().map(|l| zb[l.1]).collect();
        let cl: Vec<bool> = links.iter().map(|l| carved[l.0] && carved[l.1]).collect();
        let t = theta_links(&links, &zk, &zr, &|i| cl[i] && !lb[i]);
        format!("**{:.3}** [{:.3}, {:.3}]", t.0, t.1, t.2)
    };
    let extras = |label: &str, vc: ValleyConstruction, g: &GridF32, cond: &GridF32, bodies_cells: &[bool]| -> String {
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(g, &vd, &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(g, cond, &pre, DELIVERED_P50_M, &ss, &vd, cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = 0;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons += 1;
            }
        }
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        let co = coast(g, &ss, &dwall, 2.0 / CELL_KM);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let eo = metres(g);
        let rem = |sel: &dyn Fn(usize) -> bool| -> f64 { (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((zs1[k] - eo[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum() };
        let rai = |sel: &dyn Fn(usize) -> bool| -> f64 { (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((eo[k] - zs1[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum() };
        let (rt, ri) = (rem(&|_| true), rem(&|k| bodies_cells[k]));
        format!(
            "canyons **{canyons}** · coast: spurs near a coastal wall {rn:.4} /km · removal {rt:.1} km³ (kept bodies {ri:.1}) · raise {:.2} km³ (kept bodies {:.2}) [{label}]",
            rai(&|_| true),
            rai(&|k| bodies_cells[k])
        )
    };
    // ── ON ×1: the 14 D8 lakes, the 292, and the exact list
    let on_at = |age: f32| ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1)) };
    let on = on_at(1.0);
    eprintln!("\n────────── ON ×1 (the references) ──────────");
    let (d8, lake_on, sea_on) = {
        let wd = build_world(kn(on), None, PSEED, None);
        let (sea, _, _) = breach_of(&wd.heightmap, Some(&wd));
        let v = viz_hd_lakes_on(&wd, kn(on), PSEED, 45.0, 40.0);
        let dr = &v.drainage;
        let rcv = |c: usize| -> Option<usize> {
            let dd = dr.flow.direction[c];
            if dd == DIR_NONE {
                return None;
            }
            let (x, y) = ((c % w) as i32 + D8_DX[dd as usize], (c / w) as i32 + D8_DY[dd as usize]);
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
        };
        let d8: Vec<u32> = dr
            .lakes
            .iter()
            .filter(|l| l.area_km2 >= 1.0 && !matches!(l.lake_type, LakeType::CraterAcidic | LakeType::CraterNeutral) && rcv(l.base.outlet.1 as usize * w + l.base.outlet.0 as usize).is_some())
            .map(|l| l.base.id)
            .collect();
        (d8, dr.lake_map.clone(), sea.into_iter().collect::<HashSet<usize>>())
    };
    let list = {
        let proxy = skeleton(&s1, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v4(1.0, 0.5)), ..on }, &ss, DOMAIN_KM);
        let bo = proxy.gorge_body_of.as_ref().unwrap();
        let mut lows: Vec<u32> = Vec::new();
        for &id in &d8 {
            let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
            for k in 0..n {
                if lake_on[k] == id && bo[k] != u32::MAX {
                    *ov.entry(bo[k]).or_insert(0) += 1;
                }
            }
            if let Some((&b, _)) = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) {
                lows.push(proxy.gorge_bodies[b as usize].low);
            }
        }
        lows.sort_unstable();
        let mut l = [u32::MAX; 16];
        for (i, x) in lows.iter().enumerate() {
            l[i] = *x;
        }
        l
    };
    let mut d8_cells: HashMap<u32, Vec<usize>> = HashMap::new();
    for k in 0..n {
        if d8.contains(&lake_on[k]) {
            d8_cells.entry(lake_on[k]).or_default().push(k);
        }
    }
    eprintln!("   the 14 D8 lakes {:?} · ON ×1's breach set {} cells · the exact list {:?}", d8, sea_on.len(), list.iter().filter(|&&x| x != u32::MAX).collect::<Vec<_>>());
    // ON's construction at each age (the G-sea control on the construction)
    let mut on_con_sea: HashMap<u32, HashSet<usize>> = HashMap::new();
    let mut summary: Vec<String> = Vec::new();
    let worlds: Vec<(f32, f32)> = vec![(1.0, 0.5), (1.4, 0.5), (1.4, 0.0)];
    for (age, p) in worlds {
        if !on_con_sea.contains_key(&age.to_bits()) {
            let sk = skeleton(&s1, &on_at(age), &ss, DOMAIN_KM);
            let (b, _) = carve(&s1, &sk, &on_at(age), &ss);
            let (set, _, _) = breach_of(&b, None);
            on_con_sea.insert(age.to_bits(), set.into_iter().collect());
        }
        for ctl in ["C0", "C1", "C2", "C3"] {
            let phi0 = ctl == "C1" || ctl == "C3";
            let frz = ctl == "C2" || ctl == "C3";
            let label = format!("{ctl} ×{age} (r_world {}, p {p})", map_a(age));
            eprintln!("\n────────── {label}{}{} ──────────", if phi0 { " · φ = 0" } else { "" }, if frz { " · the design frozen through the light pass" } else { "" });
            let t = Instant::now();
            let on_sea_c = on_con_sea[&age.to_bits()].clone();
            let res = catch_unwind(AssertUnwindSafe(|| -> String {
                let g4 = GorgeRetreat { scope_lows: Some(list), phi_zero: phi0, freeze_design: frz, ..GorgeRetreat::v4(map_a(age), p) };
                let vc = ValleyConstruction { gorge_retreat: Some(g4), ..on_at(age) };
                let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
                let (built0, mk) = carve(&s1, &sk, &vc, &ss);
                let mut built = built0;
                let plain = gorge_plain(&mut built, &sk, &ss, 1.0);
                let zb = metres(&built);
                let bo = sk.gorge_body_of.clone().expect("the gate is on");
                let bodies = sk.gorge_bodies.clone();
                let rf = sk.rim_floor_m.clone().expect("the gate is on");
                let dir_sk = sk.direction.clone();
                let mut cells: Vec<Vec<usize>> = vec![Vec::new(); bodies.len()];
                for k in 0..n {
                    if bo[k] != u32::MAX {
                        cells[bo[k] as usize].push(k);
                    }
                }
                let bcells: Vec<bool> = (0..n).map(|k| bo[k] != u32::MAX).collect();
                let ring_of = |bi: usize| -> Vec<usize> {
                    let mut r: Vec<usize> = cells[bi].iter().flat_map(|&c| (0..8).map(move |k| (c, k))).map(|(c, k)| nbt(c, k)).filter(|&m| bo[m] == u32::MAX && s1.data[m] > SEA).collect();
                    r.sort_unstable();
                    r.dedup();
                    r
                };
                let rings: Vec<Vec<usize>> = (0..bodies.len()).map(ring_of).collect();
                // G-bodies
                let mut matched: BTreeMap<u32, (usize, f64)> = BTreeMap::new();
                let mut split = 0usize;
                for (&id, cs) in &d8_cells {
                    let mut ov: BTreeMap<u32, usize> = BTreeMap::new();
                    for &c in cs {
                        if bo[c] != u32::MAX {
                            *ov.entry(bo[c]).or_insert(0) += 1;
                        }
                    }
                    if ov.len() > 1 {
                        split += 1;
                    }
                    if let Some((&b, &k)) = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))) {
                        matched.insert(id, (b as usize, k as f64 / bodies[b as usize].cells as f64));
                    }
                }
                let used: HashSet<usize> = matched.values().map(|v| v.0).collect();
                let gb_ok = bodies.len() == 14 && matched.len() == 14 && used.len() == 14 && split == 0 && matched.values().all(|v| v.1 > 0.5);
                let tagged_sk = sk.gorge_falls.len();
                drop(sk);
                // the upstream of each body (G-drained)
                let mut up_body = vec![u32::MAX - 1; n];
                let mut path = Vec::new();
                for s in 0..n {
                    if up_body[s] != u32::MAX - 1 {
                        continue;
                    }
                    path.clear();
                    let mut c = s;
                    let v0;
                    loop {
                        if up_body[c] != u32::MAX - 1 {
                            v0 = up_body[c];
                            break;
                        }
                        if bo[c] != u32::MAX && c != s {
                            v0 = bo[c];
                            break;
                        }
                        path.push(c);
                        let d = dir_sk[c];
                        if d == DIR_NONE || s1.data[c] <= SEA || path.len() > 50_000 {
                            v0 = u32::MAX;
                            break;
                        }
                        c = nbt(c, d as usize);
                    }
                    for &q2 in &path {
                        up_body[q2] = v0;
                    }
                }
                // the geometric gates on a field: G-ring, G-drained, G-tag's drops
                let geom = |z: &[f32]| -> (usize, usize, usize) {
                    let (mut ring_viol, mut drn_bad, mut gorge_viol) = (0usize, 0usize, 0usize);
                    for (bi, b) in bodies.iter().enumerate() {
                        if b.col == u32::MAX {
                            continue;
                        }
                        for &c in &rings[bi] {
                            if rf[c].is_finite() && (zb[c] - rf[c]).abs() < 0.05 && zs1[c] > rf[c] + 0.05 && c != b.col as usize {
                                for k in 0..8 {
                                    let m = nbt(c, k);
                                    if bo[m] != u32::MAX || rf[m].is_finite() {
                                        continue;
                                    }
                                    let dist = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                                    if (z[c] - z[m]) / dist > tan28 {
                                        ring_viol += 1;
                                        break;
                                    }
                                }
                            }
                        }
                        if b.r_lake >= 1.999 {
                            for k in 0..n {
                                if up_body[k] == bi as u32 && z[k] < b.l_floor - 1.0 {
                                    drn_bad += 1;
                                }
                            }
                        }
                        let hf = (1.5 * b.slope * 100.0).max(10.0);
                        let glen = if b.slope > 0.0 { (b.d_g / b.slope / 1000.0).min(20.0) } else { 0.0 };
                        let mut gp = vec![b.col as usize];
                        let mut sd = 0f32;
                        let mut c = b.col as usize;
                        while sd < glen {
                            let d = dir_sk[c];
                            if d == DIR_NONE {
                                break;
                            }
                            let r2 = nbt(c, d as usize);
                            if s1.data[r2] <= SEA || bo[r2] != u32::MAX {
                                break;
                            }
                            sd += CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                            gp.push(r2);
                            c = r2;
                        }
                        for i in 1..gp.len().saturating_sub(2) {
                            if z[gp[i]] - z[gp[i + 2]] > hf {
                                gorge_viol += 1;
                            }
                        }
                    }
                    (ring_viol, drn_bad, gorge_viol)
                };
                // G-pits on a field (F142's amended measure)
                let pits = |g: &GridF32| -> usize {
                    use std::cmp::Reverse;
                    use std::collections::BinaryHeap;
                    let mut tot = 0usize;
                    for (bi, b) in bodies.iter().enumerate() {
                        if b.r_lake <= 1.0 || cells[bi].is_empty() {
                            continue;
                        }
                        let me = bi as u32;
                        let mut fill: HashMap<usize, f32> = HashMap::with_capacity(cells[bi].len());
                        let mut heap: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
                        let mut seen: HashSet<usize> = HashSet::new();
                        for &c in &cells[bi] {
                            for k in 0..8 {
                                let m = nbt(c, k);
                                if bo[m] != me && seen.insert(m) {
                                    heap.push(Reverse((g.data[m].to_bits(), m)));
                                }
                            }
                        }
                        while let Some(Reverse((lb_, c))) = heap.pop() {
                            let lev = f32::from_bits(lb_);
                            for k in 0..8 {
                                let m = nbt(c, k);
                                if bo[m] != me || fill.contains_key(&m) {
                                    continue;
                                }
                                let l = lev.max(g.data[m]);
                                fill.insert(m, l);
                                heap.push(Reverse((l.to_bits(), m)));
                            }
                        }
                        let low = *cells[bi].iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("cells");
                        let lake = ((fill[&low] - g.data[low]) * n2m > 0.1).then(|| fill[&low]);
                        for &c in &cells[bi] {
                            let d = (fill[&c] - g.data[c]) * n2m;
                            if d > 0.1 && lake.is_none_or(|ll| ((fill[&c] - ll) * n2m).abs() > 0.01) {
                                tot += 1;
                            }
                        }
                    }
                    tot
                };
                // the lake gates on a drainage: G-levels, drained absent, G-slope100, G-area (the 14)
                let lake_gates = |lm: &[u32], lakes: &[(u32, f32, f32, usize)], dirs: &[u8], cnf: &[f32]| -> String {
                    let lvl_of: HashMap<u32, (f32, f32, usize)> = lakes.iter().map(|l| (l.0, (l.1, l.2, l.3))).collect();
                    let rcv = |c: usize| -> Option<usize> {
                        let d = dirs[c];
                        if d == DIR_NONE {
                            return None;
                        }
                        let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                        if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
                    };
                    let (mut ok, mut und, mut dr, mut dr_ok, mut s_n, mut s_bad) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
                    let mut area14 = 0f64;
                    let mut s_list: Vec<f32> = Vec::new();
                    for (bi, b) in bodies.iter().enumerate() {
                        let mut ov: HashMap<u32, usize> = HashMap::new();
                        for &c in &cells[bi] {
                            if lm[c] != 0 {
                                *ov.entry(lm[c]).or_insert(0) += 1;
                            }
                        }
                        let best = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))).map(|(&id, &c)| (id, c));
                        if b.r_lake >= 1.999 {
                            dr += 1;
                            if best.is_none_or(|(_, c)| 2 * c < cells[bi].len()) {
                                dr_ok += 1;
                            }
                            continue;
                        }
                        und += 1;
                        let Some((id, _)) = best else { continue };
                        let Some(&(lv, ar, o)) = lvl_of.get(&id) else { continue };
                        if (lv - b.level).abs() <= 1.0 {
                            ok += 1;
                        }
                        area14 += ar as f64;
                        if let Some(col) = rcv(o) {
                            let (mut pth, mut dd) = (vec![col], vec![0f32]);
                            let mut c = col;
                            while *dd.last().unwrap() < 10.5 && cnf[c] > 0.0 {
                                let Some(r2) = rcv(c) else { break };
                                if lm[r2] != 0 && lm[r2] != id {
                                    break;
                                }
                                dd.push(dd.last().unwrap() + CELL_KM * if (c % w != r2 % w) && (c / w != r2 / w) { std::f32::consts::SQRT_2 } else { 1.0 });
                                pth.push(r2);
                                c = r2;
                            }
                            let mut ms = 0f32;
                            for i in 0..pth.len() {
                                if dd[i] > 10.0 {
                                    break;
                                }
                                if let Some(j) = (i..pth.len()).find(|&j| dd[j] >= dd[i] + 0.1) {
                                    ms = ms.max((cnf[pth[i]] - cnf[pth[j]]) / ((dd[j] - dd[i]) * 1000.0));
                                }
                            }
                            s_n += 1;
                            s_list.push(ms.atan().to_degrees());
                            if ms > tan28 {
                                s_bad += 1;
                            }
                        }
                    }
                    format!(
                        "G-levels **{ok}/{und}** · drained absent **{dr_ok}/{dr}** · G-slope100 > 28° **{s_bad}/{s_n}** (max {:.1}°) · G-area (the 14) **{area14:.1} km²**",
                        q(&mut s_list, 1.0)
                    )
                };
                let rim = |z: &[f32]| -> usize { bodies.iter().enumerate().filter(|(bi, b)| b.r_lake <= 1.0 && b.col != u32::MAX && b.level - rings[*bi].iter().map(|&c| z[c]).fold(f32::INFINITY, f32::min) > 10.0).count() };
                let n_rim = bodies.iter().filter(|b| b.r_lake <= 1.0 && b.col != u32::MAX).count();
                let dep_total: f64 = plain.iter().map(|p| p.volume_km3).sum();
                eprintln!(
                    "   G-bodies {} ({} kept, {} matched, split {split}) · the plain {} bodies {dep_total:.4} km³ · tagged falls in the skeleton {tagged_sk}",
                    if gb_ok { "**exactly 14, 1:1**" } else { "**FAILS**" },
                    bodies.len(),
                    matched.len(),
                    plain.len()
                );
                // ── on the construction (C0, C1)
                let mut con_line = String::from("= the construction of the matching C0/C1 (not re-run)");
                if !frz {
                    let (set_c, cond_c, d_c) = breach_of(&built, None);
                    let lk: Vec<(u32, f32, f32, usize)> = d_c.lakes.iter().map(|l| (l.base.id, l.level_m, l.area_km2, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize)).collect();
                    let lg = lake_gates(&d_c.lake_map, &lk, &d_c.flow.direction, &metres(&cond_c));
                    let (rv, db, gv) = geom(&zb);
                    let sc: HashSet<usize> = set_c.into_iter().collect();
                    con_line = format!(
                        "G-pits **{}** · {lg} · G-rim **{}/{n_rim}** · G-ring **{rv}** · G-drained **{db}** · G-tag drops **{gv}** · G-sea vs ON's construction +{} / −{} · θ {}",
                        pits(&built),
                        rim(&zb),
                        sc.difference(&on_sea_c).count(),
                        on_sea_c.difference(&sc).count(),
                        theta_on(&zb, &mk.carved)
                    );
                }
                eprintln!("   ON THE CONSTRUCTION · {con_line}");
                // ── on the final world
                let wd = build_world(kn(vc), None, PSEED, None);
                let eo = metres(&wd.heightmap);
                let (set_w, _, _) = breach_of(&wd.heightmap, Some(&wd));
                let sw: HashSet<usize> = set_w.into_iter().collect();
                let (added, missing) = (sw.difference(&sea_on).count(), sea_on.difference(&sw).count());
                let pits_w = pits(&wd.heightmap);
                let (rv, db, gv) = geom(&eo);
                let tail = catch_unwind(AssertUnwindSafe(|| viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0)));
                let v = match tail {
                    Ok(v) => v,
                    Err(e) => {
                        let m = msg(e);
                        eprintln!("   F38 · **THE INVARIANT FIRES** in the run_hd tail: {m}");
                        return format!("{label}: F38 FIRED · G-pits {pits_w} · G-sea +{added}/−{missing}");
                    }
                };
                let lk: Vec<(u32, f32, f32, usize)> = v.drainage.lakes.iter().map(|l| (l.base.id, l.level_m, l.area_km2, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize)).collect();
                let lg = lake_gates(&v.drainage.lake_map, &lk, &v.drainage.flow.direction, &metres(&v.conditioned));
                let tagged_w = wd.gorge_falls.len();
                let ex = extras(&label, vc, &wd.heightmap, &v.conditioned, &bcells);
                let fin_line = format!(
                    "G-pits **{pits_w}** · {lg} · G-rim **{}/{n_rim}** · G-ring (eroded) **{rv}** · G-drained (eroded) **{db}** · G-tag drops (eroded) **{gv}** · tagged {tagged_w} · G-sea +{added} / −{missing} · F38 does not fire · θ (eroded) {} · {ex}",
                    rim(&eo),
                    theta_on(&eo, &mk.carved)
                );
                eprintln!("   ON THE FINAL WORLD · {fin_line}");
                format!("{label} · CONSTRUCTION {con_line} || FINAL {fin_line}")
            }));
            let line = match res {
                Ok(s) => s,
                Err(e) => {
                    let m = msg(e);
                    eprintln!("   **THE WORLD PANICKED**: {m}");
                    format!("{label}: PANICKED · {}", m.chars().take(160).collect::<String>())
                }
            };
            eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
            summary.push(line);
        }
    }
    // ── the controls ON / OFF ×1
    for (label, vc) in [("ON ×1", on), ("OFF ×1", off)] {
        eprintln!("\n────────── {label} (controls) ──────────");
        let wd = build_world(kn(vc), None, PSEED, None);
        let v = viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0);
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (b, mk) = carve(&s1, &sk, &vc, &ss);
        let none = vec![false; n];
        let ex = extras(label, vc, &wd.heightmap, &v.conditioned, &none);
        let line = format!("{label}: {} lakes · θ construction {} · θ eroded {} · {ex}", v.drainage.lakes.len(), theta_on(&metres(&b), &mk.carved), theta_on(&metres(&wd.heightmap), &mk.carved));
        eprintln!("   {line}");
        summary.push(line);
    }
    eprintln!("\n   SUMMARY · {}", summary.join("\n   SUMMARY · "));
    eprintln!("\n==========  end Finding 143-C . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 144 — spec v5: the scope computed on the input (S), the outlet-base bound (B), the capped head falls
/// (H), the catchment of the laying cones (Z), then phase 1 of the light pass's candidates P0–P3 on ×1 (p 0.5), ×1.4
/// (p 0.5) and ×1.4 (p 0), with the selection rule declared in `f144_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f144_p --nocapture
#[test]
#[ignore]
fn f144_p() {
    use common::{build_world, viz_dcfg, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase, carve_diag, gorge_design_distance_m, gorge_design_mask, gorge_plain, ocean_flood};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone_protected};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let tan28 = 28f32.to_radians().tan();
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 144 . spec v5: S, B, H, Z and phase 1 of P0–P3  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1 = metres(&s1);
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let (lm_in, dep) = {
        let d = c1_drainage_windowed(&s1, None, &C1DrainageConfig::default(), &ss, DOMAIN_KM);
        let bf = breach_monotone(&s1, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let spill = ocean_flood(&bf);
        let eps = 0.01 / c1_altitude_norm_to_metres(1.0, &ss).max(1.0);
        let dep: Vec<bool> = (0..n).map(|k| bf.data[k] > C1_SEA_LEVEL_NORM && spill[k] > bf.data[k] + eps).collect();
        (d.lake_map, dep)
    };
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &vd, &ss, &dc, DOMAIN_KM).0
    };
    let off = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let links = trunk_links(&skeleton(&s1, &off, &ss, DOMAIN_KM));
    let lb: Vec<bool> = links.iter().map(|l| lm_in[l.0] != 0 || dep[l.0] || lm_in[l.1] != 0 || dep[l.1]).collect();
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let map_a = |age: f32| if age <= 1.0 { (age - 0.7) / 0.3 } else { 1.0 + (age - 1.0) / 0.4 };
    let on_at = |age: f32| ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1)) };
    let on = on_at(1.0);
    let v5 = |age: f32, p: f32| ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v5(map_a(age), p)), ..on_at(age) };
    let q = |v: &mut Vec<f32>, p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(f32::total_cmp);
        v[((v.len() - 1) as f32 * p) as usize]
    };
    let msg = |e: Box<dyn std::any::Any + Send>| -> String {
        e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_else(|| "non-string panic".to_string())
    };
    let breach_set = |wd: &common::World| -> Vec<usize> {
        let g = &wd.heightmap;
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let c = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        (0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect()
    };
    let theta_fit = |z: &[f32], carved: &[bool], excl: &dyn Fn(usize) -> bool| -> (f64, f64, f64, usize) {
        let zk: Vec<f32> = links.iter().map(|l| z[l.0]).collect();
        let zr: Vec<f32> = links.iter().map(|l| z[l.1]).collect();
        let t = theta_links(&links, &zk, &zr, &|i| carved[links[i].0] && carved[links[i].1] && !lb[i] && !excl(i));
        (t.0 as f64, t.1 as f64, t.2 as f64, t.3 as usize)
    };
    let extras = |label: &str, vc: ValleyConstruction, g: &GridF32, cond: &GridF32, bodies_cells: &[bool]| -> (usize, String) {
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(g, &vd, &ss, &dc_e, DOMAIN_KM);
        set_dump(true);
        let _cr = f95_criteria(g, cond, &pre, DELIVERED_P50_M, &ss, &vd, cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = 0;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons += 1;
            }
        }
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        let co = coast(g, &ss, &dwall, 2.0 / CELL_KM);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        let eo = metres(g);
        let rem = |sel: &dyn Fn(usize) -> bool| -> f64 { (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((zs1[k] - eo[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum() };
        let rai = |sel: &dyn Fn(usize) -> bool| -> f64 { (0..n).filter(|&k| s1.data[k] > SEA && sel(k)).map(|k| ((eo[k] - zs1[k]).max(0.0) as f64) * cell_km2 as f64 * 1e-3).sum() };
        (
            canyons,
            format!(
                "canyons **{canyons}** · coast {rn:.4} /km · removal {:.1} km³ (kept bodies {:.1}) · raise {:.2} km³ (kept bodies {:.2}) [{label}]",
                rem(&|_| true),
                rem(&|k| bodies_cells[k]),
                rai(&|_| true),
                rai(&|k| bodies_cells[k])
            ),
        )
    };
    // ── references: ON ×1 (G-sea) and OFF ×1 (θ on the eroded world)
    let sea_on: HashSet<usize> = {
        let wd = build_world(kn(on), None, PSEED, None);
        breach_set(&wd).into_iter().collect()
    };
    let off_ci = {
        let wd = build_world(kn(off), None, PSEED, None);
        let sk = skeleton(&s1, &off, &ss, DOMAIN_KM);
        let (_, mk) = carve(&s1, &sk, &off, &ss);
        theta_fit(&metres(&wd.heightmap), &mk.carved, &|_| false)
    };
    eprintln!("   ON ×1's breach set {} cells · OFF ×1's θ on the eroded world {:.3} [{:.3}, {:.3}] ({} links)", sea_on.len(), off_ci.0, off_ci.1, off_ci.2, off_ci.3);
    // ── S
    {
        let list14: [u32; 14] = [20360735, 23179228, 25325689, 29042877, 30033743, 32780258, 33206583, 34376571, 34784695, 38292872, 41298524, 41535879, 47099127, 48588290];
        let sk5 = skeleton(&s1, &v5(1.0, 0.5), &ss, DOMAIN_KM);
        let sk3 = skeleton(&s1, &ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v3(1.0, 0.5)), ..on }, &ss, DOMAIN_KM);
        let kept: Vec<u32> = sk5.gorge_bodies.iter().map(|b| b.low).collect();
        eprintln!(
            "\n   S · v5 keeps **{}** bodies · of F143's 14: {} kept, missing {:?} · beyond the 14: {:?}",
            kept.len(),
            list14.iter().filter(|l| kept.contains(l)).count(),
            list14.iter().filter(|l| !kept.contains(l)).map(|&l| xy(l as usize)).collect::<Vec<_>>(),
            kept.iter().filter(|l| !list14.contains(l)).map(|&l| xy(l as usize)).collect::<Vec<_>>()
        );
        for b in sk3.gorge_bodies.iter().filter(|b| b.input_lake) {
            eprintln!(
                "      input lake at {} · {:.1} km² · in a depression {} · under the sea {} · ring touches a depression {} · v5 {}",
                xy(b.low as usize),
                b.cells as f32 * cell_km2,
                b.in_depression,
                b.in_sea,
                b.touches_depression,
                if kept.contains(&b.low) { "KEPT" } else { "dropped" }
            );
        }
    }
    // ── B (the constructions' bodies at every age)
    eprintln!("\n   B · the bodies the outlet's base raises (spec level → level, B):");
    for (age, p) in [(0.7f32, 0.5f32), (0.85, 0.5), (1.0, 0.5), (1.2, 0.5), (1.4, 0.5), (1.0, 0.0), (1.0, 0.25), (1.4, 0.0)] {
        let sk = skeleton(&s1, &v5(age, p), &ss, DOMAIN_KM);
        let b4 = sk.gorge_bodies.iter().find(|b| b.low == 29042877);
        let raised: Vec<String> = sk
            .gorge_bodies
            .iter()
            .filter(|b| b.bounded)
            .map(|b| format!("{} r {:.2} → {:.1} (B {:.1}, L_floor {:.1})", xy(b.low as usize), b.r_lake, b.level, b.b_out, b.l_floor))
            .collect();
        eprintln!(
            "      ×{age} p {p}: {} raised {:?} · lake 4: {}",
            raised.len(),
            raised,
            b4.map_or("not kept".to_string(), |b| format!("r {:.2} · L_floor {:.1} · B {:.1} · level {:.1} · bounded {}", b.r_lake, b.l_floor, b.b_out, b.level, b.bounded))
        );
    }
    // ── H and the G-tag instrument (×1, p 0.5, the construction)
    {
        let vc = v5(1.0, 0.5);
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (b0, _) = carve(&s1, &sk, &vc, &ss);
        let mut bp = b0;
        let _ = gorge_plain(&mut bp, &sk, &ss, 1.0);
        let zb = metres(&bp);
        let bo = sk.gorge_body_of.clone().unwrap();
        eprintln!("\n   H (×1, p 0.5) · A, D_g, φ, the head fall (capped), the cap, S, H_f:");
        let (mut tagged, mut heads) = (0usize, Vec::new());
        let (mut old_n, mut new_n) = (0usize, 0usize);
        for b in &sk.gorge_bodies {
            let hf = (1.5 * b.slope * 100.0).max(10.0);
            let tag = b.head_fall_m > hf;
            if tag {
                tagged += 1;
                heads.push(b.head_fall_m);
            }
            eprintln!(
                "      {} · A {:>7.1} · D_g {:>6.1} · φ {:.2} · head {:>6.1} m (cap {:.0}) · S {:.2}° · H_f {hf:.1}{}",
                xy(b.low as usize),
                b.a_out_km2,
                b.d_g,
                b.phi,
                b.head_fall_m,
                b.head_cap_m,
                b.slope.atan().to_degrees(),
                if tag { " TAG" } else { "" }
            );
            if b.col == u32::MAX {
                continue;
            }
            let glen = if b.slope > 0.0 { (b.d_g / b.slope / 1000.0).min(20.0) } else { 0.0 };
            let (mut gp, mut sd) = (vec![b.col as usize], vec![0f32]);
            let mut c = b.col as usize;
            while *sd.last().unwrap() < glen {
                let d = sk.direction[c];
                if d == DIR_NONE {
                    break;
                }
                let r2 = nbt(c, d as usize);
                if s1.data[r2] <= SEA || bo[r2] != u32::MAX {
                    break;
                }
                sd.push(sd.last().unwrap() + CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 });
                gp.push(r2);
                c = r2;
            }
            for i in 1..gp.len().saturating_sub(2) {
                let dz = zb[gp[i]] - zb[gp[i + 2]];
                if dz > hf {
                    old_n += 1;
                }
                if dz - b.slope * (sd[i + 2] - sd[i]) * 1000.0 > hf {
                    new_n += 1;
                }
            }
        }
        heads.sort_by(f32::total_cmp);
        eprintln!(
            "   H · tagged **{tagged}** of {} · median {:.1} m · max {:.1} m · G-tag on the construction: old instrument **{old_n}** · corrected **{new_n}** · (H_f in the code: `(1.5 * sg * 100.0).max(10.0)` with sg = slope_for(top), the real slope)",
            sk.gorge_bodies.len(),
            if heads.is_empty() { f32::NAN } else { heads[heads.len() / 2] },
            heads.last().copied().unwrap_or(f32::NAN)
        );
    }
    // ── Z (×1.4, p 0.5 and p 0): the laying cone's catchment
    for (age, p) in [(1.4f32, 0.5f32), (1.4, 0.0)] {
        let vc = v5(age, p);
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let (co, mk, dg) = carve_diag(&s1, &sk, &vc, &ss);
        let mut bp = co;
        let _ = gorge_plain(&mut bp, &sk, &ss, 1.0);
        let zp = metres(&bp);
        drop(bp);
        let bo = sk.gorge_body_of.clone().unwrap();
        let owner = sk.line_owner.clone().unwrap();
        let dir = &sk.direction;
        // the first line met down the D8 path (memoised)
        let mut first = vec![u32::MAX - 1; n];
        let mut path = Vec::new();
        for s in 0..n {
            if first[s] != u32::MAX - 1 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if first[c] != u32::MAX - 1 {
                    v0 = first[c];
                    break;
                }
                if owner[c] != u32::MAX {
                    v0 = owner[c];
                    break;
                }
                path.push(c);
                let d = dir[c];
                if d == DIR_NONE || path.len() > 50_000 {
                    v0 = u32::MAX;
                    break;
                }
                c = nbt(c, d as usize);
            }
            for &q2 in &path {
                first[q2] = v0;
            }
        }
        let chain_has = |start: u32, target: u32| -> bool {
            let mut l = start;
            for _ in 0..sk.polylines.len() + 1 {
                if l == u32::MAX {
                    return false;
                }
                if l == target {
                    return true;
                }
                match sk.line_parent.get(l as usize).copied().flatten() {
                    Some((pl, _)) => l = pl,
                    None => return false,
                }
            }
            false
        };
        let cone_line = |c: usize| -> Option<u32> {
            let s = dg.banded.as_ref().map(|b| b[c]).filter(|&s| s != u32::MAX).unwrap_or(dg.who[c]);
            (s != u32::MAX).then(|| dg.line_of[s as usize])
        };
        let outside = |c: usize| -> Option<bool> {
            let l = cone_line(c)?;
            let f = if owner[c] != u32::MAX { owner[c] } else { first[c] };
            Some(!chain_has(f, l))
        };
        // up_body for G-drained
        let mut up_body = vec![u32::MAX - 1; n];
        for s in 0..n {
            if up_body[s] != u32::MAX - 1 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if up_body[c] != u32::MAX - 1 {
                    v0 = up_body[c];
                    break;
                }
                if bo[c] != u32::MAX && c != s {
                    v0 = bo[c];
                    break;
                }
                path.push(c);
                let d = dir[c];
                if d == DIR_NONE || s1.data[c] <= SEA || path.len() > 50_000 {
                    v0 = u32::MAX;
                    break;
                }
                c = nbt(c, d as usize);
            }
            for &q2 in &path {
                up_body[q2] = v0;
            }
        }
        let drained: Vec<bool> = sk.gorge_bodies.iter().map(|b| b.r_lake >= 1.999).collect();
        let share = |sel: &dyn Fn(usize) -> bool| -> String {
            let (mut t, mut o, mut nos) = (0usize, 0usize, 0usize);
            for c in 0..n {
                if !sel(c) {
                    continue;
                }
                match outside(c) {
                    Some(true) => {
                        t += 1;
                        o += 1;
                    }
                    Some(false) => t += 1,
                    None => nos += 1,
                }
            }
            format!("{o} of {t} laid by a cone ({:.1} %) outside their catchment (+ {nos} without a laying sample)", if t > 0 { 100.0 * o as f64 / t as f64 } else { 0.0 })
        };
        let viol = |c: usize| -> bool {
            let ub = up_body[c];
            ub != u32::MAX && ub != u32::MAX - 1 && drained[ub as usize] && zp[c] < sk.gorge_bodies[ub as usize].l_floor - 1.0
        };
        eprintln!("\n   Z (×{age}, p {p}) · G-drained's violators: {}", share(&viol));
        eprintln!("      the drained bowls' footprints: {}", share(&|c| bo[c] != u32::MAX && drained[bo[c] as usize]));
        eprintln!("      the cells upstream of the drained bowls: {}", share(&|c| up_body[c] != u32::MAX && up_body[c] != u32::MAX - 1 && drained[up_body[c] as usize]));
        eprintln!("      every carved cell of the world: {}", share(&|c| mk.carved[c]));
    }
    // ── P, phase 1
    let mut held_by: BTreeMap<String, (usize, usize, bool)> = BTreeMap::new(); // (held, G-ring cells, G-pits in every world)
    let mut summary: Vec<String> = Vec::new();
    let mut d_t_m: f32 = 0.0;
    let mut p0_x1_eo: Option<Vec<f32>> = None;
    for (age, p) in [(1.0f32, 0.5f32), (1.4, 0.5), (1.4, 0.0)] {
        for cand in ["P0", "C2-for-d_t", "P1", "P2", "P3"] {
            if cand == "C2-for-d_t" && age != 1.0 {
                continue;
            }
            let mode: u8 = match cand {
                "P1" => 1,
                "P2" => 2,
                "P3" => 3,
                _ => 0,
            };
            let base = GorgeRetreat::v5(map_a(age), p);
            let g5 = GorgeRetreat { light_mode: mode, light_dt_m: if mode == 2 { d_t_m } else { 0.0 }, freeze_design: cand == "C2-for-d_t", ..base };
            let vc = ValleyConstruction { gorge_retreat: Some(g5), ..on_at(age) };
            let label = format!("{cand} ×{age} (r_world {}, p {p}){}", map_a(age), if mode == 2 { format!(" · d_t {d_t_m:.0} m") } else { String::new() });
            eprintln!("\n────────── {label} ──────────");
            let t = Instant::now();
            if cand == "C2-for-d_t" {
                // d_t: the width of C2's walls against P0, by distance band outside the design mask
                let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
                let mask = gorge_design_mask(&sk);
                let dist = gorge_design_distance_m(&mask, &sk);
                drop(sk);
                let wd = build_world(kn(vc), None, PSEED, None);
                let ec2 = metres(&wd.heightmap);
                let e0 = p0_x1_eo.take().expect("P0 ×1 first");
                let steep = |z: &[f32], c: usize| -> bool {
                    (0..8).any(|k| {
                        let m = nbt(c, k);
                        let dd = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                        (z[c] - z[m]) / dd > tan28
                    })
                };
                let mut band = vec![(0usize, 0usize, 0usize); 21];
                for c in 0..n {
                    let j = (dist[c] / (CELL_KM * 1000.0)).round() as usize;
                    if (1..=20).contains(&j) {
                        band[j].0 += 1;
                        band[j].1 += steep(&ec2, c) as usize;
                        band[j].2 += steep(&e0, c) as usize;
                    }
                }
                let mut chosen = 0usize;
                let mut rows = Vec::new();
                for j in 1..=20 {
                    let (tot, s2, s0) = band[j];
                    let ex = if tot > 0 { 100.0 * (s2 as f64 - s0 as f64) / tot as f64 } else { 0.0 };
                    rows.push(format!("{j}: {ex:+.2} pp"));
                    if chosen == 0 && ex < 0.5 {
                        chosen = j;
                    }
                }
                d_t_m = chosen.max(1) as f32 * CELL_KM * 1000.0;
                eprintln!("   d_t · C2's steep share minus P0's, by band (cells): {} · **d_t = {} cells = {d_t_m:.0} m**", rows.join(" · "), chosen.max(1));
                continue;
            }
            let res = catch_unwind(AssertUnwindSafe(|| -> (Vec<bool>, usize, String, Option<Vec<f32>>) {
                let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
                let (built0, mk) = carve(&s1, &sk, &vc, &ss);
                let mut built = built0;
                let plain = gorge_plain(&mut built, &sk, &ss, 1.0);
                let zb = metres(&built);
                let bo = sk.gorge_body_of.clone().expect("the gate is on");
                let bodies = sk.gorge_bodies.clone();
                let rf = sk.rim_floor_m.clone().expect("the gate is on");
                let dir_sk = sk.direction.clone();
                let mask = gorge_design_mask(&sk);
                let mut corridor = vec![false; n];
                for &c in &sk.gorge_path {
                    for dy in -2i32..=2 {
                        for dx in -2i32..=2 {
                            corridor[((c as usize / w) as i32 + dy).rem_euclid(h as i32) as usize * w + ((c as usize % w) as i32 + dx).rem_euclid(w as i32) as usize] = true;
                        }
                    }
                }
                drop(sk);
                let mut cells: Vec<Vec<usize>> = vec![Vec::new(); bodies.len()];
                for k in 0..n {
                    if bo[k] != u32::MAX {
                        cells[bo[k] as usize].push(k);
                    }
                }
                let bcells: Vec<bool> = (0..n).map(|k| bo[k] != u32::MAX).collect();
                let rings: Vec<Vec<usize>> = (0..bodies.len())
                    .map(|bi| {
                        let mut r: Vec<usize> = cells[bi].iter().flat_map(|&c| (0..8).map(move |k| (c, k))).map(|(c, k)| nbt(c, k)).filter(|&m| bo[m] == u32::MAX && s1.data[m] > SEA).collect();
                        r.sort_unstable();
                        r.dedup();
                        r
                    })
                    .collect();
                let mut up_body = vec![u32::MAX - 1; n];
                let mut path = Vec::new();
                for s in 0..n {
                    if up_body[s] != u32::MAX - 1 {
                        continue;
                    }
                    path.clear();
                    let mut c = s;
                    let v0;
                    loop {
                        if up_body[c] != u32::MAX - 1 {
                            v0 = up_body[c];
                            break;
                        }
                        if bo[c] != u32::MAX && c != s {
                            v0 = bo[c];
                            break;
                        }
                        path.push(c);
                        let d = dir_sk[c];
                        if d == DIR_NONE || s1.data[c] <= SEA || path.len() > 50_000 {
                            v0 = u32::MAX;
                            break;
                        }
                        c = nbt(c, d as usize);
                    }
                    for &q2 in &path {
                        up_body[q2] = v0;
                    }
                }
                // the geometric gates on a field
                let geom = |z: &[f32]| -> (usize, usize, usize, usize) {
                    let (mut ring_viol, mut drn_bad, mut tag_new, mut edge) = (0usize, 0usize, 0usize, 0usize);
                    for (bi, b) in bodies.iter().enumerate() {
                        if b.col == u32::MAX {
                            continue;
                        }
                        for &c in &rings[bi] {
                            if rf[c].is_finite() && (zb[c] - rf[c]).abs() < 0.05 && zs1[c] > rf[c] + 0.05 && c != b.col as usize {
                                for k in 0..8 {
                                    let m = nbt(c, k);
                                    if bo[m] != u32::MAX || rf[m].is_finite() {
                                        continue;
                                    }
                                    let dist = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                                    if (z[c] - z[m]) / dist > tan28 {
                                        ring_viol += 1;
                                        break;
                                    }
                                }
                            }
                        }
                        if b.r_lake >= 1.999 {
                            for k in 0..n {
                                if up_body[k] == bi as u32 && z[k] < b.l_floor - 1.0 {
                                    drn_bad += 1;
                                }
                            }
                        }
                        let hf = (1.5 * b.slope * 100.0).max(10.0);
                        let glen = if b.slope > 0.0 { (b.d_g / b.slope / 1000.0).min(20.0) } else { 0.0 };
                        let (mut gp, mut sd) = (vec![b.col as usize], vec![0f32]);
                        let mut c = b.col as usize;
                        while *sd.last().unwrap() < glen {
                            let d = dir_sk[c];
                            if d == DIR_NONE {
                                break;
                            }
                            let r2 = nbt(c, d as usize);
                            if s1.data[r2] <= SEA || bo[r2] != u32::MAX {
                                break;
                            }
                            sd.push(sd.last().unwrap() + CELL_KM * if d % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 });
                            gp.push(r2);
                            c = r2;
                        }
                        for i in 1..gp.len().saturating_sub(2) {
                            if (z[gp[i]] - z[gp[i + 2]]) - b.slope * (sd[i + 2] - sd[i]) * 1000.0 > hf {
                                tag_new += 1;
                            }
                        }
                    }
                    for c in 0..n {
                        if !mask[c] {
                            continue;
                        }
                        if (0..8).any(|k| {
                            let m = nbt(c, k);
                            let dd = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                            !mask[m] && (z[c] - z[m]) / dd > tan28
                        }) {
                            edge += 1;
                        }
                    }
                    (ring_viol, drn_bad, tag_new, edge)
                };
                let pits = |g: &GridF32| -> usize {
                    use std::cmp::Reverse;
                    use std::collections::BinaryHeap;
                    let mut tot = 0usize;
                    for (bi, b) in bodies.iter().enumerate() {
                        if b.r_lake <= 1.0 || cells[bi].is_empty() {
                            continue;
                        }
                        let me = bi as u32;
                        let mut fill: HashMap<usize, f32> = HashMap::with_capacity(cells[bi].len());
                        let mut heap: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
                        let mut seen: HashSet<usize> = HashSet::new();
                        for &c in &cells[bi] {
                            for k in 0..8 {
                                let m = nbt(c, k);
                                if bo[m] != me && seen.insert(m) {
                                    heap.push(Reverse((g.data[m].to_bits(), m)));
                                }
                            }
                        }
                        while let Some(Reverse((lb_, c))) = heap.pop() {
                            let lev = f32::from_bits(lb_);
                            for k in 0..8 {
                                let m = nbt(c, k);
                                if bo[m] != me || fill.contains_key(&m) {
                                    continue;
                                }
                                let l = lev.max(g.data[m]);
                                fill.insert(m, l);
                                heap.push(Reverse((l.to_bits(), m)));
                            }
                        }
                        let low = *cells[bi].iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("cells");
                        let lake = ((fill[&low] - g.data[low]) * n2m > 0.1).then(|| fill[&low]);
                        for &c in &cells[bi] {
                            let d = (fill[&c] - g.data[c]) * n2m;
                            if d > 0.1 && lake.is_none_or(|ll| ((fill[&c] - ll) * n2m).abs() > 0.01) {
                                tot += 1;
                            }
                        }
                    }
                    tot
                };
                let (_, _, tag_c, edge_c) = geom(&zb);
                let tb = Instant::now();
                let wd = build_world(kn(vc), None, PSEED, None);
                let t_build = tb.elapsed().as_secs_f64();
                let eo = metres(&wd.heightmap);
                let sw: HashSet<usize> = breach_set(&wd).into_iter().collect();
                let (added, missing) = (sw.difference(&sea_on).count(), sea_on.difference(&sw).count());
                let pits_w = pits(&wd.heightmap);
                let (rv, db, tg, edge) = geom(&eo);
                let tail = catch_unwind(AssertUnwindSafe(|| viz_hd_lakes_on(&wd, kn(vc), PSEED, 45.0, 40.0)));
                let v = match tail {
                    Ok(v) => v,
                    Err(e) => {
                        let m = msg(e);
                        eprintln!("   F38 · **THE INVARIANT FIRES**: {m}");
                        let held = vec![false; 14];
                        return (held, rv, format!("{label}: F38 FIRED · G-pits {pits_w} · G-sea +{added}/−{missing}"), None);
                    }
                };
                let lm = &v.drainage.lake_map;
                let lakes = &v.drainage.lakes;
                let cn = metres(&v.conditioned);
                let lvl_of: HashMap<u32, (f32, f32, usize)> = lakes.iter().map(|l| (l.base.id, (l.level_m, l.area_km2, l.base.outlet.1 as usize * w + l.base.outlet.0 as usize))).collect();
                let rcv = |c: usize| -> Option<usize> {
                    let d = v.drainage.flow.direction[c];
                    if d == DIR_NONE {
                        return None;
                    }
                    let (x, y) = ((c % w) as i32 + D8_DX[d as usize], (c / w) as i32 + D8_DY[d as usize]);
                    if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 { None } else { Some(y as usize * w + x as usize) }
                };
                let (mut ok_lv, mut und, mut dr, mut dr_ok, mut s_bad, mut s_n) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
                let (mut gb_ok, mut gb_n) = (0usize, 0usize);
                let mut used: HashSet<u32> = HashSet::new();
                let mut area = 0f64;
                let mut rim_bad = 0usize;
                let mut s_list: Vec<f32> = Vec::new();
                let mut off_lv: Vec<String> = Vec::new();
                for (bi, b) in bodies.iter().enumerate() {
                    let mut ov: HashMap<u32, usize> = HashMap::new();
                    for &c in &cells[bi] {
                        if lm[c] != 0 {
                            *ov.entry(lm[c]).or_insert(0) += 1;
                        }
                    }
                    let best = ov.iter().max_by_key(|e| (*e.1, std::cmp::Reverse(*e.0))).map(|(&id, &c)| (id, c));
                    gb_n += 1;
                    if b.r_lake >= 1.999 {
                        dr += 1;
                        if best.is_none_or(|(_, c)| 2 * c < cells[bi].len()) {
                            dr_ok += 1;
                            gb_ok += 1;
                        }
                        continue;
                    }
                    und += 1;
                    if b.r_lake <= 1.0 && b.col != u32::MAX && b.level - rings[bi].iter().map(|&c| eo[c]).fold(f32::INFINITY, f32::min) > 10.0 {
                        rim_bad += 1;
                    }
                    let Some((id, c)) = best else {
                        off_lv.push(format!("{} r {:.2} L {:.1}: no lake", xy(b.low as usize), b.r_lake, b.level));
                        continue;
                    };
                    if 2 * c > cells[bi].len() && used.insert(id) {
                        gb_ok += 1;
                    }
                    let Some(&(lv, ar, o)) = lvl_of.get(&id) else { continue };
                    area += ar as f64;
                    if (lv - b.level).abs() <= 1.0 {
                        ok_lv += 1;
                    } else {
                        off_lv.push(format!("{} r {:.2} L {:.1} final {lv:.1}", xy(b.low as usize), b.r_lake, b.level));
                    }
                    if let Some(col) = rcv(o) {
                        let (mut pth, mut dd) = (vec![col], vec![0f32]);
                        let mut c = col;
                        while *dd.last().unwrap() < 10.5 && cn[c] > 0.0 {
                            let Some(r2) = rcv(c) else { break };
                            if lm[r2] != 0 && lm[r2] != id {
                                break;
                            }
                            dd.push(dd.last().unwrap() + CELL_KM * if (c % w != r2 % w) && (c / w != r2 / w) { std::f32::consts::SQRT_2 } else { 1.0 });
                            pth.push(r2);
                            c = r2;
                        }
                        let mut ms = 0f32;
                        for i in 0..pth.len() {
                            if dd[i] > 10.0 {
                                break;
                            }
                            if let Some(j) = (i..pth.len()).find(|&j| dd[j] >= dd[i] + 0.1) {
                                ms = ms.max((cn[pth[i]] - cn[pth[j]]) / ((dd[j] - dd[i]) * 1000.0));
                            }
                        }
                        s_n += 1;
                        s_list.push(ms.atan().to_degrees());
                        if ms > tan28 {
                            s_bad += 1;
                        }
                    }
                }
                let th = theta_fit(&eo, &mk.carved, &|i| corridor[links[i].0] || corridor[links[i].1]);
                let th_in = th.0 >= off_ci.1 && th.0 <= off_ci.2;
                let (canyons, ex) = extras(&label, vc, &wd.heightmap, &v.conditioned, &bcells);
                let dep_total: f64 = plain.iter().map(|p| p.volume_km3).sum();
                let held = vec![
                    gb_ok == gb_n,
                    pits_w == 0,
                    ok_lv == und,
                    dr_ok == dr,
                    rim_bad == 0,
                    rv == 0,
                    s_bad == 0,
                    db == 0,
                    tg == 0,
                    added == 0 && missing == 0,
                    true,
                    canyons == 0,
                    th_in,
                    edge == 0,
                ];
                let names = ["G-bodies", "G-pits", "G-levels", "drained", "G-rim", "G-ring", "G-slope100", "G-drained", "G-tag", "G-sea", "F38", "canyons", "θ", "edge"];
                let failed: Vec<&str> = held.iter().zip(names.iter()).filter(|(h2, _)| !**h2).map(|(_, nm)| *nm).collect();
                let line = format!(
                    "**{} of 14 held** (failed {:?}) · {} kept bodies · G-bodies {gb_ok}/{gb_n} · G-pits {pits_w} · G-levels {ok_lv}/{und} · drained absent {dr_ok}/{dr} · G-rim {rim_bad} · G-ring {rv} · G-slope100 {s_bad}/{s_n} (max {:.1}°) · G-drained {db} · G-tag (corrected) {tg} (construction {tag_c}) · tagged {} · G-sea +{added}/−{missing} · F38 no · θ (corridors excluded) {:.3} [{:.3}, {:.3}] vs OFF eroded [{:.3}, {:.3}] {} · edge {edge} (construction {edge_c}) · G-area {area:.1} km² · plain {dep_total:.4} km³ · build {t_build:.0} s · {ex}",
                    held.iter().filter(|x| **x).count(),
                    failed,
                    bodies.len(),
                    q(&mut s_list, 1.0),
                    wd.gorge_falls.len(),
                    th.0,
                    th.1,
                    th.2,
                    off_ci.1,
                    off_ci.2,
                    if th_in { "IN" } else { "OUT" }
                );
                eprintln!("   {line}");
                for o in off_lv.iter().take(12) {
                    eprintln!("      off: {o}");
                }
                let keep_eo = (age == 1.0 && cand == "P0").then_some(eo);
                (held, rv, line, keep_eo)
            }));
            let (held, rv, line) = match res {
                Ok((hd, rv, l, keep)) => {
                    if keep.is_some() {
                        p0_x1_eo = keep;
                    }
                    (hd, rv, l)
                }
                Err(e) => {
                    let m = msg(e);
                    eprintln!("   **THE WORLD PANICKED**: {m}");
                    (vec![false; 14], 0, format!("{label}: PANICKED · {}", m.chars().take(160).collect::<String>()))
                }
            };
            eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
            let e = held_by.entry(cand.to_string()).or_insert((0, 0, true));
            e.0 += held.iter().filter(|x| **x).count();
            e.1 += rv;
            e.2 &= held[1];
            summary.push(format!("{label}: {line}"));
        }
    }
    eprintln!("\n   PHASE 1 · held of 42 (14 gates × 3 worlds) · G-ring cells · G-pits in every world:");
    for (c, (hd, rg, pits_all)) in &held_by {
        eprintln!("      {c}: **{hd}** held · G-ring {rg} · G-pits everywhere {pits_all}");
    }
    let cands: Vec<(&String, &(usize, usize, bool))> = held_by.iter().filter(|e| e.0 != "P0").collect();
    let any_pits = cands.iter().any(|e| e.1.2);
    let best = cands.iter().max_by(|a, b| a.1.0.cmp(&b.1.0).then(b.1.1.cmp(&a.1.1)));
    eprintln!(
        "   SELECTION (declared rule) · {}",
        if !any_pits { "NO candidate holds G-pits in all three worlds: REPORT AND STOP (no phase 2)".to_string() } else { format!("selected **{}**", best.map_or("—".to_string(), |e| e.0.clone())) }
    );
    eprintln!("\n   SUMMARY · {}", summary.join("\n   SUMMARY · "));
    eprintln!("\n==========  end Finding 144 phase 1 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 145 — the base's first defect: `carve_diag`'s cones laying cells outside their line's basin. M1 (Δz of
/// the foreign-laid cells against the own-basin carve), M2 (F127's dams), M3 (lowered divides, artificial captures),
/// then the candidate C1 (the own-basin carve) against the témoin at the construction stage. A test-only copy of
/// `carve_diag` for this configuration, checked bit for bit with its filter off. Declared in `f145_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f145_m --nocapture
#[test]
#[ignore]
fn f145_m() {
    use common::viz_dcfg;
    use std::collections::{BTreeMap, BinaryHeap};
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, carve_diag};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};

    // the heap item of carve_diag (min-heap on the distance, ties on the cell index), verbatim
    struct It(f32, usize);
    impl PartialEq for It {
        fn eq(&self, o: &Self) -> bool {
            self.cmp(o) == std::cmp::Ordering::Equal
        }
    }
    impl Eq for It {}
    impl PartialOrd for It {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for It {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            o.0.total_cmp(&self.0).then_with(|| o.1.cmp(&self.1))
        }
    }
    /// `carve_diag`'s path for a configuration with no confluence band, no wall profile, no wall-sea floor and no rim,
    /// with an optional filter: a sample may lay a cell only if its line's basin is the cell's basin.
    fn carve_own(field: &GridF32, sk: &Skeleton, vc: &ValleyConstruction, ss: &SteinSteinParams, filt: Option<(&[u32], &[u32])>) -> (GridF32, Vec<bool>) {
        assert!(!vc.trunk_band && vc.wall_profile.is_none() && vc.wall_sea_floor_m.is_none() && sk.rim_floor_m.is_none(), "the copy covers this configuration only");
        let (w, h) = (sk.width, sk.height);
        let n = w * h;
        let tan = vc.wall_deg.to_radians().tan();
        let zfield: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
        let mut src: Vec<(f32, f32, f32, f32)> = Vec::new();
        let mut line_of: Vec<u32> = Vec::new();
        for (li, l) in sk.polylines.iter().enumerate() {
            for &p in l.iter() {
                src.push(p);
                line_of.push(li as u32);
            }
        }
        let ok = |s: usize, c: usize| -> bool { filt.is_none_or(|(cell_b, line_b)| line_b[line_of[s] as usize] == cell_b[c]) };
        let foot = 0.0f32;
        let geo = |s: usize, c: usize| -> (f32, f32, bool, f32) {
            let (sx, sy, zf, hw) = src[s];
            let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
            let mut dx = (cx - sx).rem_euclid(w as f32);
            if dx > w as f32 / 2.0 {
                dx -= w as f32;
            }
            let mut dy = (cy - sy).rem_euclid(h as f32);
            if dy > h as f32 / 2.0 {
                dy -= h as f32;
            }
            let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
            let u = (d - hw).max(0.0);
            let rise = if foot > 0.0 && u < foot { tan * u * u / (2.0 * foot) } else { tan * (u - 0.5 * foot) };
            (d, zf + rise, d <= hw, u)
        };
        let mut bestd = vec![f32::INFINITY; n];
        let mut who = vec![u32::MAX; n];
        let mut heap = BinaryHeap::new();
        for (i, sp) in src.iter().enumerate() {
            let cx = (sp.0.floor() as i64).rem_euclid(w as i64) as usize;
            let cy = (sp.1.floor() as i64).rem_euclid(h as i64) as usize;
            let c = cy * w + cx;
            let (d, v, _, _) = geo(i, c);
            if d < bestd[c] && v < zfield[c] && ok(i, c) {
                bestd[c] = d;
                who[c] = i as u32;
                heap.push(It(d, c));
            }
        }
        while let Some(It(d, c)) = heap.pop() {
            if d > bestd[c] {
                continue;
            }
            let s = who[c] as usize;
            let (x, y) = ((c % w) as i32, (c / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let nb = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                    let (dd, vv, _, _) = geo(s, nb);
                    if dd < bestd[nb] && vv < zfield[nb] && ok(s, nb) {
                        bestd[nb] = dd;
                        who[nb] = s as u32;
                        heap.push(It(dd, nb));
                    }
                }
            }
        }
        let mut out = field.clone();
        let mut carved = vec![false; n];
        for c in 0..n {
            if who[c] == u32::MAX {
                continue;
            }
            let me = who[c] as usize;
            let (x, y) = ((c % w) as i32, (c / w) as i32);
            let (_, mut v, _, _) = geo(me, c);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nb = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                    let o = who[nb];
                    if o != u32::MAX && line_of[o as usize] != line_of[me] && ok(o as usize, c) {
                        let (_, vv, _, _) = geo(o as usize, c);
                        if vv < v {
                            v = vv;
                        }
                    }
                }
            }
            if v < zfield[c] {
                let nv = c1_metres_to_altitude_norm(v, ss);
                if nv < out.data[c] {
                    out.data[c] = nv;
                    carved[c] = true;
                }
            }
        }
        (out, carved)
    }

    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let tan28 = 28f32.to_radians().tan();
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 145 . the cones outside their basin (M1–M3, C1 at the construction)  ==========");
    let metres = |g: &GridF32| -> Vec<f32> { g.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect() };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    let fill_pre = {
        let cl = c1_climate_placed(&pre, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
        fill_field_m(&pre, &vd, &ss, &dc, DOMAIN_KM).0
    };
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let links = trunk_links(&skeleton(&s1, &temoin, &ss, DOMAIN_KM));
    // the route of every cell on a field: the skeleton's flow chain, labelled by the last land cell on its D8 path
    let routes = |f: &GridF32| -> Vec<u32> {
        let dcfg = C1DrainageConfig::default();
        let d = c1_drainage_windowed(f, None, &dcfg, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let fl = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dcfg.flat_perturbation.clone(), dinf: dcfg.dinf });
        let land: Vec<bool> = bf.data.iter().map(|&v| v > C1_SEA_LEVEL_NORM).collect();
        let mut lab = vec![u32::MAX - 1; n];
        let mut path = Vec::new();
        for s in 0..n {
            if lab[s] != u32::MAX - 1 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if !land[c] {
                    v0 = path.last().map_or(u32::MAX, |&p| p as u32);
                    break;
                }
                if lab[c] != u32::MAX - 1 {
                    v0 = lab[c];
                    break;
                }
                path.push(c);
                let dd = fl.direction[c];
                if dd == DIR_NONE || path.len() > 100_000 {
                    v0 = c as u32;
                    break;
                }
                c = nbt(c, dd as usize);
            }
            for &q in &path {
                lab[q] = v0;
            }
            if !land[s] {
                lab[s] = u32::MAX;
            }
        }
        lab
    };
    // the extras on a field (canyons, the relief via F95's dump, the coast)
    let extras = |label: &str, vc: ValleyConstruction, g: &GridF32| -> (usize, String) {
        let cl_e = c1_climate_placed(g, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dc_e = DrainageClimate { precip_internal: &cl_e.precipitation, temperature: &cl_e.temperature };
        let (fill_del, _) = fill_field_m(g, &vd, &ss, &dc_e, DOMAIN_KM);
        let cond = {
            let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
            breach_monotone(g, &d.flow.filled, &d.lake_map, 0.5, w, h)
        };
        eprintln!("   [F95 dump for {label}]");
        set_dump(true);
        let _cr = f95_criteria(g, &cond, &pre, DELIVERED_P50_M, &ss, &vd, cell_km2, n2m, w, h);
        set_dump(false);
        let mut canyons = 0;
        for b in take_bodies() {
            let floor = *b.cells.iter().min_by(|&&a, &&c| g.data[a].total_cmp(&g.data[c])).expect("body");
            if over_dug_depression(fill_del[floor] - fill_pre[floor], b.rim) {
                canyons += 1;
            }
        }
        let skp = skeleton(&pre, &vc, &ss, DOMAIN_KM);
        let (_, mkp) = carve(&pre, &skp, &vc, &ss);
        let sea: Vec<bool> = (0..n).map(|k| g.data[k] <= SEA).collect();
        let dsea = dist_from(&sea, w, h);
        let cw: Vec<bool> = (0..n).map(|k| mkp.carved[k] && !mkp.floor[k] && g.data[k] > SEA && dsea[k] <= 3).collect();
        let dwall = dist_from(&cw, w, h);
        let co = coast(g, &ss, &dwall, 2.0 / CELL_KM);
        let rn = co.near.iter().filter(|&&b| b).count() as f32 / co.l_near_km.max(1e-6);
        (canyons, format!("canyons **{canyons}** · coast: spurs near a coastal wall {rn:.4} /km"))
    };
    let r8_net = |g: &GridF32| -> f32 {
        let d = c1_drainage_windowed(g, None, &vd, &ss, DOMAIN_KM);
        let segs = &d.rivers.segments;
        let pts: Vec<&[(u32, u32)]> = (0..segs.len()).filter(|&i| d.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.len() >= 2).map(|i| segs[i].points.as_slice()).collect();
        r8_chords(&pts, 8)
    };
    let route_in = routes(&s1);
    let cols = [(3852usize, 2337usize), (4551, 3280), (2253, 4495), (1805, 4580)];
    for (wname, vc) in [
        ("the témoin (C2 /10 col)", temoin),
        ("ON extended", ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..temoin }),
        ("B2 → A_c (F127's canyons)", ValleyConstruction { a_min_km2: 0.1, ..temoin }),
    ] {
        eprintln!("\n────────── {wname} ──────────");
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        // the basin label: the last land cell on the skeleton's D8 path (= the route on the input)
        let basin = &route_in;
        // the line's basin: the majority over its sample cells
        let line_b: Vec<u32> = sk
            .polylines
            .iter()
            .map(|l| {
                let mut cnt: BTreeMap<u32, usize> = BTreeMap::new();
                for p in l {
                    let c = (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize;
                    *cnt.entry(basin[c]).or_insert(0) += 1;
                }
                cnt.into_iter().max_by_key(|e| (e.1, std::cmp::Reverse(e.0))).map_or(u32::MAX, |e| e.0)
            })
            .collect();
        let tp = Instant::now();
        let (g_prod, mk, dg) = carve_diag(&s1, &sk, &vc, &ss);
        let t_prod = tp.elapsed().as_secs_f64();
        let (g_copy, _) = carve_own(&s1, &sk, &vc, &ss, None);
        let diff = (0..n).filter(|&k| g_copy.data[k].to_bits() != g_prod.data[k].to_bits()).count();
        eprintln!("   the copy (filter off) against production carve_diag: {diff} cells differ{}", if diff == 0 { " (bit-identical: read)" } else { " — NOT READ" });
        drop(g_copy);
        if diff != 0 {
            continue;
        }
        let to = Instant::now();
        let (g_own, carved_own) = carve_own(&s1, &sk, &vc, &ss, Some((basin.as_slice(), line_b.as_slice())));
        let t_own = to.elapsed().as_secs_f64();
        let zp = metres(&g_prod);
        let zo = metres(&g_own);
        let carved_n = mk.carved.iter().filter(|&&x| x).count();
        let line_at = |c: usize| -> Option<u32> { (dg.who[c] != u32::MAX).then(|| dg.line_of[dg.who[c] as usize]) };
        let land: Vec<bool> = (0..n).map(|k| s1.data[k] > SEA).collect();
        let divide: Vec<bool> = (0..n).map(|c| land[c] && (0..8).any(|k| { let m = nbt(c, k); land[m] && basin[m] != basin[c] })).collect();
        let ddiv = dist_from(&divide, w, h);
        // M1
        let mut dz_all: Vec<f32> = Vec::new();
        let mut by_band: [(usize, usize, usize); 3] = [(0, 0, 0); 3]; // (foreign, deepened < −1, |Δz| ≤ 1)
        let (mut foreign, mut defect, mut legit) = (0usize, 0usize, 0usize);
        for c in 0..n {
            if !mk.carved[c] {
                continue;
            }
            let Some(l) = line_at(c) else { continue };
            if line_b[l as usize] == basin[c] {
                continue;
            }
            foreign += 1;
            let dz = zp[c] - zo[c];
            dz_all.push(dz);
            let b = if ddiv[c] <= 2 { 0 } else if ddiv[c] <= 10 { 1 } else { 2 };
            by_band[b].0 += 1;
            if dz < -1.0 {
                defect += 1;
                by_band[b].1 += 1;
            }
            if dz.abs() <= 1.0 {
                by_band[b].2 += 1;
                if ddiv[c] <= 2 {
                    legit += 1;
                }
            }
        }
        let mut s = dz_all.clone();
        s.sort_by(f32::total_cmp);
        let pc = |p: f64| if s.is_empty() { f32::NAN } else { s[((s.len() - 1) as f64 * p) as usize] };
        eprintln!(
            "   M1 · carved cells {carved_n} · laid by a FOREIGN-basin line **{foreign}** ({:.1} %) · Δz p5/p25/p50/p75/p95 {:.1}/{:.1}/{:.1}/{:.1}/{:.1} m",
            100.0 * foreign as f64 / carved_n.max(1) as f64,
            pc(0.05),
            pc(0.25),
            pc(0.5),
            pc(0.75),
            pc(0.95)
        );
        for (i, nm) in ["0–2 cells from a divide", "3–10", "> 10"].iter().enumerate() {
            let (f, dd, l) = by_band[i];
            eprintln!("      {nm}: {f} foreign · deepened (Δz < −1 m) {dd} · |Δz| ≤ 1 m {l}");
        }
        eprintln!(
            "   M1 · **THE REAL DEFECT (Δz < −1 m): {defect} cells = {:.2} % of the carved cells** · the legitimate (≤ 2 cells from a divide, |Δz| ≤ 1 m): {legit} = {:.1} % of the foreign-laid",
            100.0 * defect as f64 / carved_n.max(1) as f64,
            100.0 * legit as f64 / foreign.max(1) as f64
        );
        // M2: F127's cols (B2 → A_c)
        if wname.starts_with("B2") {
            for &(x, y) in &cols {
                let c = y * w + x;
                let l = line_at(c);
                eprintln!(
                    "   M2 · col {} · carved {} · laid by line {:?} ({}) · z_prod {:.1} m · z_own {:.1} m · Δz {:+.1} m · real defect {}",
                    xy(c),
                    mk.carved[c],
                    l,
                    match l {
                        Some(l) if line_b[l as usize] != basin[c] => "FOREIGN",
                        Some(_) => "own basin",
                        None => "no sample",
                    },
                    zp[c],
                    zo[c],
                    zp[c] - zo[c],
                    l.is_some_and(|l| line_b[l as usize] != basin[c]) && zp[c] - zo[c] < -1.0
                );
            }
        }
        // M3: lowered divides and artificial captures
        let mut low: Vec<f32> = (0..n).filter(|&c| divide[c] && zp[c] < zo[c] - 1.0).map(|c| zo[c] - zp[c]).collect();
        low.sort_by(f32::total_cmp);
        let rp = routes(&g_prod);
        let ro = routes(&g_own);
        let mut cap = 0usize;
        let mut pairs: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        let mut cap_own = 0usize;
        for c in 0..n {
            if !land[c] || route_in[c] == u32::MAX {
                continue;
            }
            let to_in = |r: u32| if r == u32::MAX || r == u32::MAX - 1 { u32::MAX } else { route_in[r as usize] };
            let (a, b, o) = (route_in[c], to_in(rp[c]), to_in(ro[c]));
            if b != a && o == a {
                cap += 1;
                *pairs.entry((a, b)).or_insert(0) += 1;
            }
            if o != a {
                cap_own += 1;
            }
        }
        let mut pv: Vec<((u32, u32), usize)> = pairs.into_iter().collect();
        pv.sort_by(|x, y| y.1.cmp(&x.1));
        eprintln!(
            "   M3 · divide cells lowered > 1 m by the cones the own-basin carve does not lay: **{}** · lowering p50 {:.1} m · p90 {:.1} m",
            low.len(),
            if low.is_empty() { f32::NAN } else { low[low.len() / 2] },
            if low.is_empty() { f32::NAN } else { low[(low.len() - 1) * 9 / 10] }
        );
        eprintln!(
            "   M3 · **ARTIFICIAL CAPTURES** (route changed under the témoin's carve, kept under the own-basin carve): **{cap} cells = {:.1} km²** in {} basin pairs · cells whose route changes even under the own-basin carve: {cap_own}",
            cap as f32 * cell_km2,
            pv.len()
        );
        for ((a, b), k) in pv.iter().take(6) {
            eprintln!("      {k} cells ({:.2} km²) from the basin of {} to that of {}", *k as f32 * cell_km2, xy(*a as usize), if *b == u32::MAX { "—".to_string() } else { xy(*b as usize) });
        }
        let stop = (defect as f64) < 0.01 * carved_n as f64 && cap == 0;
        eprintln!("   STOP RULE · real defect < 1 % and no capture: {}", if stop { "**TRIGGERED — C is not run on this world**" } else { "not triggered" });
        if stop {
            continue;
        }
        // C1 against the témoin's carve (rule 18), at the construction
        let walls = |z: &[f32]| -> usize {
            let mut k_ = 0usize;
            for c in 0..n {
                if !land[c] {
                    continue;
                }
                for k in [0usize, 1, 2, 3] {
                    let m = nbt(c, k);
                    if land[m] && basin[m] != basin[c] {
                        let dd = CELL_KM * 1000.0 * if k % 2 == 1 { std::f32::consts::SQRT_2 } else { 1.0 };
                        if ((z[c] - z[m]) / dd).abs() > tan28 {
                            k_ += 1;
                        }
                    }
                }
            }
            k_
        };
        let th = |z: &[f32], carved: &[bool]| -> String {
            let zk: Vec<f32> = links.iter().map(|l| z[l.0]).collect();
            let zr: Vec<f32> = links.iter().map(|l| z[l.1]).collect();
            let t = theta_links(&links, &zk, &zr, &|i| carved[links[i].0] && carved[links[i].1]);
            format!("{:.3} [{:.3}, {:.3}]", t.0, t.1, t.2)
        };
        let cls = classes(&pre, &sk, &mk);
        let (cy_p, ex_p) = extras(&format!("{wname}, the témoin's carve"), vc, &g_prod);
        let (cy_o, ex_o) = extras(&format!("{wname}, C1"), vc, &g_own);
        let _ = (cy_p, cy_o);
        eprintln!(
            "   C · the témoin's carve: walls across the divides **{}** · θ {} · R8 terrain {:.4} · R8 network {:.4} · {ex_p}",
            walls(&zp),
            th(&zp, &mk.carved),
            r8_terrain(&g_prod, &cls, 16),
            r8_net(&g_prod)
        );
        eprintln!(
            "   C · C1 (the own-basin carve): walls across the divides **{}** · θ {} · R8 terrain {:.4} · R8 network {:.4} · {ex_o} · captures under C1 (route changed against the input) {}",
            walls(&zo),
            th(&zo, &carved_own),
            r8_terrain(&g_own, &cls, 16),
            r8_net(&g_own),
            cap_own
        );
        if wname.starts_with("B2") {
            for &(x, y) in &cols {
                let c = y * w + x;
                eprintln!("   C · F127 col {} · the témoin's carve {:.1} m · C1 {:.1} m", xy(c), zp[c], zo[c]);
            }
        }
        eprintln!(
            "   C · COST · production carve {t_prod:.1} s · the bench's filtered carve {t_own:.1} s ({:+.1} s; the basin labels come free from compute_flow in production) · against run_hd 249.8 s",
            t_own - t_prod
        );
    }
    eprintln!("\n==========  end Finding 145 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 145-M3, amendment — the artificial captures refined: only between input basins of at least 10 km²,
/// against the own-basin carve's noise floor. Declared in `f145_declared.md` (amendment, non-blind).
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f145_cap --nocapture
#[test]
#[ignore]
fn f145_cap() {
    use std::collections::{BTreeMap, BinaryHeap, HashMap};
    use ymir_core::tectonics_c1::drainage::C1_SEA_LEVEL_NORM;
    use ymir_core::tectonics_c1::production_upscale::c1_metres_to_altitude_norm;
    use ymir_core::tectonics_c1::valley_construction::{LakeBase, carve_diag};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};
    // the heap item of carve_diag (min-heap on the distance, ties on the cell index), verbatim
    struct It(f32, usize);
    impl PartialEq for It {
        fn eq(&self, o: &Self) -> bool {
            self.cmp(o) == std::cmp::Ordering::Equal
        }
    }
    impl Eq for It {}
    impl PartialOrd for It {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for It {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            o.0.total_cmp(&self.0).then_with(|| o.1.cmp(&self.1))
        }
    }
    /// `carve_diag`'s path for a configuration with no confluence band, no wall profile, no wall-sea floor and no rim,
    /// with an optional filter: a sample may lay a cell only if its line's basin is the cell's basin.
    fn carve_own(field: &GridF32, sk: &Skeleton, vc: &ValleyConstruction, ss: &SteinSteinParams, filt: Option<(&[u32], &[u32])>) -> (GridF32, Vec<bool>) {
        assert!(!vc.trunk_band && vc.wall_profile.is_none() && vc.wall_sea_floor_m.is_none() && sk.rim_floor_m.is_none(), "the copy covers this configuration only");
        let (w, h) = (sk.width, sk.height);
        let n = w * h;
        let tan = vc.wall_deg.to_radians().tan();
        let zfield: Vec<f32> = field.data.iter().map(|&v| c1_altitude_norm_to_metres(v, ss)).collect();
        let mut src: Vec<(f32, f32, f32, f32)> = Vec::new();
        let mut line_of: Vec<u32> = Vec::new();
        for (li, l) in sk.polylines.iter().enumerate() {
            for &p in l.iter() {
                src.push(p);
                line_of.push(li as u32);
            }
        }
        let ok = |s: usize, c: usize| -> bool { filt.is_none_or(|(cell_b, line_b)| line_b[line_of[s] as usize] == cell_b[c]) };
        let foot = 0.0f32;
        let geo = |s: usize, c: usize| -> (f32, f32, bool, f32) {
            let (sx, sy, zf, hw) = src[s];
            let (cx, cy) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
            let mut dx = (cx - sx).rem_euclid(w as f32);
            if dx > w as f32 / 2.0 {
                dx -= w as f32;
            }
            let mut dy = (cy - sy).rem_euclid(h as f32);
            if dy > h as f32 / 2.0 {
                dy -= h as f32;
            }
            let d = (dx * dx + dy * dy).sqrt() * sk.cell_m;
            let u = (d - hw).max(0.0);
            let rise = if foot > 0.0 && u < foot { tan * u * u / (2.0 * foot) } else { tan * (u - 0.5 * foot) };
            (d, zf + rise, d <= hw, u)
        };
        let mut bestd = vec![f32::INFINITY; n];
        let mut who = vec![u32::MAX; n];
        let mut heap = BinaryHeap::new();
        for (i, sp) in src.iter().enumerate() {
            let cx = (sp.0.floor() as i64).rem_euclid(w as i64) as usize;
            let cy = (sp.1.floor() as i64).rem_euclid(h as i64) as usize;
            let c = cy * w + cx;
            let (d, v, _, _) = geo(i, c);
            if d < bestd[c] && v < zfield[c] && ok(i, c) {
                bestd[c] = d;
                who[c] = i as u32;
                heap.push(It(d, c));
            }
        }
        while let Some(It(d, c)) = heap.pop() {
            if d > bestd[c] {
                continue;
            }
            let s = who[c] as usize;
            let (x, y) = ((c % w) as i32, (c / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let nb = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                    let (dd, vv, _, _) = geo(s, nb);
                    if dd < bestd[nb] && vv < zfield[nb] && ok(s, nb) {
                        bestd[nb] = dd;
                        who[nb] = s as u32;
                        heap.push(It(dd, nb));
                    }
                }
            }
        }
        let mut out = field.clone();
        let mut carved = vec![false; n];
        for c in 0..n {
            if who[c] == u32::MAX {
                continue;
            }
            let me = who[c] as usize;
            let (x, y) = ((c % w) as i32, (c / w) as i32);
            let (_, mut v, _, _) = geo(me, c);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nb = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                    let o = who[nb];
                    if o != u32::MAX && line_of[o as usize] != line_of[me] && ok(o as usize, c) {
                        let (_, vv, _, _) = geo(o as usize, c);
                        if vv < v {
                            v = vv;
                        }
                    }
                }
            }
            if v < zfield[c] {
                let nv = c1_metres_to_altitude_norm(v, ss);
                if nv < out.data[c] {
                    out.data[c] = nv;
                    carved[c] = true;
                }
            }
        }
        (out, carved)
    }

    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 145-M3 amendment . the captures between basins of at least 10 km²  ==========");
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let xy = |c: usize| format!("({},{})", c % w, c / w);
    let nbt = |c: usize, k: usize| ((c / w) as i32 + D8_DY[k]).rem_euclid(h as i32) as usize * w + ((c % w) as i32 + D8_DX[k]).rem_euclid(w as i32) as usize;
    // the route of every cell on a field: the skeleton's flow chain, labelled by the last land cell on its D8 path
    let routes = |f: &GridF32| -> Vec<u32> {
        let dcfg = C1DrainageConfig::default();
        let d = c1_drainage_windowed(f, None, &dcfg, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, C1_SEA_LEVEL_NORM, w, h);
        let fl = compute_flow(&bf, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, flat_perturbation: dcfg.flat_perturbation.clone(), dinf: dcfg.dinf });
        let land: Vec<bool> = bf.data.iter().map(|&v| v > C1_SEA_LEVEL_NORM).collect();
        let mut lab = vec![u32::MAX - 1; n];
        let mut path = Vec::new();
        for s in 0..n {
            if lab[s] != u32::MAX - 1 {
                continue;
            }
            path.clear();
            let mut c = s;
            let v0;
            loop {
                if !land[c] {
                    v0 = path.last().map_or(u32::MAX, |&p| p as u32);
                    break;
                }
                if lab[c] != u32::MAX - 1 {
                    v0 = lab[c];
                    break;
                }
                path.push(c);
                let dd = fl.direction[c];
                if dd == DIR_NONE || path.len() > 100_000 {
                    v0 = c as u32;
                    break;
                }
                c = nbt(c, dd as usize);
            }
            for &q in &path {
                lab[q] = v0;
            }
            if !land[s] {
                lab[s] = u32::MAX;
            }
        }
        lab
    };
    let route_in = routes(&s1);
    let mut area: HashMap<u32, usize> = HashMap::new();
    for &r in &route_in {
        if r != u32::MAX {
            *area.entry(r).or_insert(0) += 1;
        }
    }
    let big = |r: u32| r != u32::MAX && area.get(&r).copied().unwrap_or(0) as f32 * cell_km2 >= 10.0;
    eprintln!("   input basins: {} · of them ≥ 10 km²: {} (holding {:.1} % of the land)", area.len(), area.keys().filter(|&&r| big(r)).count(),
        100.0 * area.iter().filter(|e| big(*e.0)).map(|e| *e.1).sum::<usize>() as f64 / area.values().sum::<usize>().max(1) as f64);
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    for (wname, vc) in [("the témoin (C2 /10 col)", temoin), ("ON extended", ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..temoin })] {
        eprintln!("\n────────── {wname} ──────────");
        let sk = skeleton(&s1, &vc, &ss, DOMAIN_KM);
        let line_b: Vec<u32> = sk
            .polylines
            .iter()
            .map(|l| {
                let mut cnt: BTreeMap<u32, usize> = BTreeMap::new();
                for p in l {
                    let c = (p.1.floor() as i64).rem_euclid(h as i64) as usize * w + (p.0.floor() as i64).rem_euclid(w as i64) as usize;
                    *cnt.entry(route_in[c]).or_insert(0) += 1;
                }
                cnt.into_iter().max_by_key(|e| (e.1, std::cmp::Reverse(e.0))).map_or(u32::MAX, |e| e.0)
            })
            .collect();
        let (g_prod, _, _) = carve_diag(&s1, &sk, &vc, &ss);
        let (g_own, _) = carve_own(&s1, &sk, &vc, &ss, Some((route_in.as_slice(), line_b.as_slice())));
        drop(sk);
        let rp = routes(&g_prod);
        let ro = routes(&g_own);
        let to_in = |r: u32| if r == u32::MAX || r == u32::MAX - 1 { u32::MAX } else { route_in[r as usize] };
        let (mut cap, mut noise, mut raw) = (0usize, 0usize, 0usize);
        let mut pairs: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        for c in 0..n {
            let a = route_in[c];
            if !big(a) {
                continue;
            }
            let (b, o) = (to_in(rp[c]), to_in(ro[c]));
            if b != a {
                raw += 1;
            }
            if b != a && big(b) && o == a {
                cap += 1;
                *pairs.entry((a, b)).or_insert(0) += 1;
            }
            if o != a && big(o) {
                noise += 1;
            }
        }
        let mut pv: Vec<((u32, u32), usize)> = pairs.into_iter().collect();
        pv.sort_by(|x, y| y.1.cmp(&x.1));
        eprintln!(
            "   cells of basins ≥ 10 km² routed elsewhere under the témoin's carve: {raw} · **captured into another basin ≥ 10 km² and kept by the own-basin carve: {cap} cells = {:.1} km²** in {} pairs · the noise floor (the own-basin carve moving cells between basins ≥ 10 km²): {noise} cells = {:.1} km²",
            cap as f32 * cell_km2,
            pv.len(),
            noise as f32 * cell_km2
        );
        for ((a, b), k) in pv.iter().take(10) {
            eprintln!(
                "      {:.2} km² from the basin of {} ({:.1} km²) to that of {} ({:.1} km²)",
                *k as f32 * cell_km2,
                xy(*a as usize),
                area[a] as f32 * cell_km2,
                xy(*b as usize),
                area[b] as f32 * cell_km2
            );
        }
    }
    eprintln!("\n==========  end Finding 145-M3 amendment . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 146 — the rivers for Living Landz on the témoin: `rivers_ll` built from the viz's HD drainage, the
/// smoothing's instruments (L), the parallel bundles (P), the selection attributes (A) and the cost. Declared in
/// `f146_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f146_rivers --nocapture
#[test]
#[ignore]
fn f146_rivers() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use ymir_core::export::rivers_ll::{RiverEnd, RiversLlParams, build_rivers_ll, rivers_ll_json};
    use ymir_core::export::hydro::rivers_json;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 146 . the rivers for Living Landz (témoin C2 /10 col)  ==========");
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let kn = Knobs { valley: Some(temoin), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let dr = &v.drainage;
    let (w, h) = (v.conditioned.width, v.conditioned.height);
    let n = w * h;
    let rj0 = rivers_json(dr, cell_km2);
    // production cost (no instruments), then the instrumented build
    let tp = Instant::now();
    let (rivers, _) = build_rivers_ll(dr, &v.conditioned, &ss, cell_km2, RiversLlParams::default(), false);
    let t_build = tp.elapsed().as_secs_f64();
    let tj = Instant::now();
    let bytes = rivers_ll_json(&rivers, CELL_KM, &RiversLlParams::default());
    let t_json = tj.elapsed().as_secs_f64();
    let (_, st) = build_rivers_ll(dr, &v.conditioned, &ss, cell_km2, RiversLlParams::default(), true);
    let rj1 = rivers_json(dr, cell_km2);
    eprintln!(
        "   {} segments → **{} rivers** · {} smoothed points · rivers_ll.json {:.1} MB · rivers.json {:.1} MB, byte-identical before/after the build: {}",
        dr.rivers.segments.len(),
        rivers.len(),
        rivers.iter().map(|r| r.points.len()).sum::<usize>(),
        bytes.len() as f64 / 1e6,
        rj0.len() as f64 / 1e6,
        rj0 == rj1
    );
    eprintln!("   COST · build {t_build:.2} s + serialize {t_json:.2} s per world (against run_hd 249.8 s)");
    // ── L
    let pct = |v: &[f32], p: f64| -> f32 {
        let mut s = v.to_vec();
        s.sort_by(f32::total_cmp);
        if s.is_empty() { f32::NAN } else { s[((s.len() - 1) as f64 * p) as usize] }
    };
    let out = |v: &[f32]| 100.0 * v.iter().filter(|&&x| x > 1.0).count() as f64 / v.len().max(1) as f64;
    eprintln!("\n   L · terrain − bed (m), samples every 0.25 cell: p50 / p99 / max · share of length above the 1 m tolerance");
    for (nm, vv) in [("raw D8 trace (negative control)", &st.above_bed_raw), ("smoothed, before the constraint", &st.above_bed_before), ("smoothed, after the constraint", &st.above_bed_after)] {
        eprintln!("      {nm}: {:.2} / {:.2} / {:.2} m · **{:.2} %**", pct(vv, 0.5), pct(vv, 0.99), pct(vv, 1.0), out(vv));
    }
    eprintln!(
        "   L · vertices pulled back by the constraint {} of {} ({:.1} %) · ridge crossings: raw {} of {} samples · smoothed **{}** of {} · lateral deviation from the trace p50 {:.3} · p99 {:.3} cell · stair index (turns ≥ 45°): raw {:.3} → smoothed {:.3}",
        st.pulled,
        st.vertices,
        100.0 * st.pulled as f64 / st.vertices.max(1) as f64,
        st.ridge_cross_raw,
        st.samples_raw,
        st.ridge_cross_smoothed,
        st.samples,
        pct(&st.lateral_cells, 0.5),
        pct(&st.lateral_cells, 0.99),
        st.stair_raw,
        st.stair_smoothed
    );
    // crossings between rivers (away from a shared confluence) and self-intersections, on a 4-cell hash
    let seg_inter = |a: [f32; 2], b: [f32; 2], c: [f32; 2], d: [f32; 2]| -> Option<[f32; 2]> {
        let r = [b[0] - a[0], b[1] - a[1]];
        let s = [d[0] - c[0], d[1] - c[1]];
        let den = r[0] * s[1] - r[1] * s[0];
        if den.abs() < 1e-9 {
            return None;
        }
        let t = ((c[0] - a[0]) * s[1] - (c[1] - a[1]) * s[0]) / den;
        let u = ((c[0] - a[0]) * r[1] - (c[1] - a[1]) * r[0]) / den;
        ((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u)).then(|| [a[0] + t * r[0], a[1] + t * r[1]])
    };
    let crossings = |polys: &[Vec<[f32; 2]>]| -> (usize, usize) {
        let mut hash: HashMap<(i32, i32), Vec<(usize, usize)>> = HashMap::new();
        for (ri, p) in polys.iter().enumerate() {
            for e in 0..p.len().saturating_sub(1) {
                let (a, b) = (p[e], p[e + 1]);
                let (x0, x1) = ((a[0].min(b[0]) / 4.0).floor() as i32, (a[0].max(b[0]) / 4.0).floor() as i32);
                let (y0, y1) = ((a[1].min(b[1]) / 4.0).floor() as i32, (a[1].max(b[1]) / 4.0).floor() as i32);
                for gx in x0..=x1 {
                    for gy in y0..=y1 {
                        hash.entry((gx, gy)).or_default().push((ri, e));
                    }
                }
            }
        }
        let ends: Vec<([f32; 2], [f32; 2])> = polys.iter().map(|p| (p[0], *p.last().unwrap())).collect();
        let near_end = |q: [f32; 2], ri: usize| -> bool {
            let (a, b) = ends[ri];
            ((q[0] - a[0]).powi(2) + (q[1] - a[1]).powi(2)).sqrt() < 0.05 || ((q[0] - b[0]).powi(2) + (q[1] - b[1]).powi(2)).sqrt() < 0.05
        };
        let mut seen: HashSet<(usize, usize, usize, usize)> = HashSet::new();
        let (mut cross, mut selfx) = (0usize, 0usize);
        for list in hash.values() {
            for i in 0..list.len() {
                for j in i + 1..list.len() {
                    let ((ra, ea), (rb, eb)) = (list[i], list[j]);
                    let key = if (ra, ea) < (rb, eb) { (ra, ea, rb, eb) } else { (rb, eb, ra, ea) };
                    if !seen.insert(key) {
                        continue;
                    }
                    if ra == rb && ea.abs_diff(eb) <= 1 {
                        continue;
                    }
                    let (pa, pb) = (&polys[ra], &polys[rb]);
                    if let Some(q) = seg_inter(pa[ea], pa[ea + 1], pb[eb], pb[eb + 1]) {
                        if ra == rb {
                            selfx += 1;
                        } else if !(near_end(q, ra) || near_end(q, rb)) {
                            cross += 1;
                        }
                    }
                }
            }
        }
        (cross, selfx)
    };
    let smoothed: Vec<Vec<[f32; 2]>> = rivers.iter().map(|r| r.points.clone()).collect();
    let raw: Vec<Vec<[f32; 2]>> = rivers
        .iter()
        .map(|r| {
            let mut p: Vec<[f32; 2]> = Vec::new();
            for &s in &r.segments {
                for &(x, y) in &dr.rivers.segments[s].points {
                    let q = [x as f32 + 0.5, y as f32 + 0.5];
                    if p.last() != Some(&q) {
                        p.push(q);
                    }
                }
            }
            p
        })
        .collect();
    let (cs, ss_) = crossings(&smoothed);
    let (cr, sr) = crossings(&raw);
    eprintln!("   L · crossings between two rivers away from a confluence: raw {cr} · smoothed **{cs}** · self-intersections: raw {sr} · smoothed **{ss_}**");
    // ── P: the parallel bundles
    let sk = skeleton(&build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED), &temoin, &ss, DOMAIN_KM);
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (_, mk) = carve(&s1, &sk, &temoin, &ss);
    drop(sk);
    let joined: HashSet<(u32, u32)> = rivers
        .iter()
        .filter_map(|r| match r.end {
            RiverEnd::Confluence { river_id } => Some((r.id.min(river_id), r.id.max(river_id))),
            _ => None,
        })
        .collect();
    // samples every cell along each smoothed river
    let samples: Vec<Vec<[f32; 2]>> = smoothed
        .iter()
        .map(|p| {
            let mut s = Vec::new();
            for e in 0..p.len().saturating_sub(1) {
                let (a, b) = (p[e], p[e + 1]);
                let l = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
                let m = l.ceil().max(1.0) as usize;
                for j in 0..m {
                    let t = j as f32 / m as f32;
                    s.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
                }
            }
            if let Some(&l) = p.last() {
                s.push(l);
            }
            s
        })
        .collect();
    let mut grid: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for (ri, s) in samples.iter().enumerate() {
        for q in s {
            grid.entry((q[0].floor() as i32, q[1].floor() as i32)).or_default().push(ri);
        }
    }
    let k = 2i32;
    let mut runs: HashMap<(u32, u32), (usize, Vec<usize>)> = HashMap::new(); // longest run (samples) and its cells
    for (ri, s) in samples.iter().enumerate() {
        let mut cur: HashMap<usize, (usize, Vec<usize>)> = HashMap::new();
        for q in s {
            let (cx, cy) = (q[0].floor() as i32, q[1].floor() as i32);
            let mut near: HashSet<usize> = HashSet::new();
            for dy in -k..=k {
                for dx in -k..=k {
                    if let Some(l) = grid.get(&(cx + dx, cy + dy)) {
                        for &o in l {
                            if o != ri {
                                near.insert(o);
                            }
                        }
                    }
                }
            }
            let cell = (cy.rem_euclid(h as i32) as usize) * w + cx.rem_euclid(w as i32) as usize;
            cur.retain(|o, _| near.contains(o));
            for &o in &near {
                let e = cur.entry(o).or_insert((0, Vec::new()));
                e.0 += 1;
                e.1.push(cell);
                let key = ((ri as u32).min(o as u32), (ri as u32).max(o as u32));
                if joined.contains(&key) {
                    continue;
                }
                let best = runs.entry(key).or_insert((0, Vec::new()));
                if e.0 > best.0 {
                    *best = (e.0, e.1.clone());
                }
            }
        }
    }
    let l_cells = (2.0 / CELL_KM).ceil() as usize;
    let mut pairs: Vec<((u32, u32), (usize, Vec<usize>))> = runs.into_iter().filter(|e| e.1.0 > l_cells).collect();
    pairs.sort_by(|a, b| b.1.0.cmp(&a.1.0).then(a.0.cmp(&b.0)));
    let near_volcano = |c: usize| -> bool {
        let (x, y) = ((c % w) as f32 + 0.5, (c / w) as f32 + 0.5);
        wd.craters.iter().any(|cr| ((x - cr.center_px.0).powi(2) + (y - cr.center_px.1).powi(2)).sqrt() < 4.0 * cr.radius_px)
    };
    let mut loc: BTreeMap<&str, usize> = BTreeMap::new();
    let mut by_area: BTreeMap<&str, (usize, f64)> = BTreeMap::new();
    let mut total_km = 0f64;
    for ((a, b), (len, cells)) in &pairs {
        let km = *len as f64 * CELL_KM as f64;
        total_km += km;
        for &c in cells {
            let cat = if near_volcano(c) {
                "volcanic edifice (≤ 4 crater radii)"
            } else if mk.carved[c] && !mk.floor[c] {
                "construction planar wall (F124)"
            } else if mk.floor[c] {
                "construction valley floor"
            } else {
                "elsewhere"
            };
            *loc.entry(cat).or_insert(0) += 1;
        }
        let small = rivers[*a as usize].catchment_km2.min(rivers[*b as usize].catchment_km2);
        let band = if small < 1.0 {
            "< 1 km²"
        } else if small < 10.0 {
            "1–10 km²"
        } else if small < 100.0 {
            "10–100 km²"
        } else {
            "≥ 100 km²"
        };
        let e = by_area.entry(band).or_insert((0, 0.0));
        e.0 += 1;
        e.1 += km;
    }
    let loc_tot: usize = loc.values().sum();
    eprintln!("\n   P · parallel pairs (within {k} cells over > 2 km, no shared confluence): **{}** · their runs {total_km:.1} km", pairs.len());
    eprintln!("      by the smaller river's catchment: {:?}", by_area.iter().map(|(b, (c, km))| format!("{b}: {c} pairs, {km:.1} km")).collect::<Vec<_>>());
    eprintln!("      where the runs lie: {:?}", loc.iter().map(|(c, k2)| format!("{c}: {:.1} %", 100.0 * *k2 as f64 / loc_tot.max(1) as f64)).collect::<Vec<_>>());
    for ((a, b), (len, cells)) in pairs.iter().take(8) {
        let (ra, rb) = (&rivers[*a as usize], &rivers[*b as usize]);
        eprintln!(
            "      e.g. rivers {a} ({:.1} km², {:.2} m³/s) and {b} ({:.1} km², {:.2} m³/s) · run {:.1} km at ({},{})",
            ra.catchment_km2,
            ra.discharge_m3s,
            rb.catchment_km2,
            rb.discharge_m3s,
            *len as f64 * CELL_KM as f64,
            cells[0] % w,
            cells[0] / w
        );
    }
    // the proposed rules (not applied)
    let lower: HashSet<u32> = pairs.iter().map(|((a, b), _)| if rivers[*a as usize].discharge_m3s >= rivers[*b as usize].discharge_m3s { *b } else { *a }).collect();
    let km_lower: f32 = lower.iter().map(|&r| rivers[r as usize].length_km).sum();
    eprintln!("   P · rule (a) merge each pair into its higher-discharge river: removes {} rivers, {km_lower:.0} km", lower.len());
    for thr in [1.0f32, 3.0, 10.0] {
        let gone: Vec<&_> = rivers.iter().filter(|r| r.catchment_km2 < thr).collect();
        let pairs_left = pairs.iter().filter(|((a, b), _)| rivers[*a as usize].catchment_km2 >= thr && rivers[*b as usize].catchment_km2 >= thr).count();
        eprintln!(
            "   P · rule (b) drop the rivers under {thr} km²: removes {} rivers, {:.0} km · parallel pairs left {pairs_left}",
            gone.len(),
            gone.iter().map(|r| r.length_km).sum::<f32>()
        );
    }
    // ── A: the selection attributes
    let mut lens: Vec<f32> = rivers.iter().map(|r| r.length_km).collect();
    lens.sort_by(f32::total_cmp);
    let mut areas: Vec<f32> = rivers.iter().map(|r| r.catchment_km2).collect();
    areas.sort_by(f32::total_cmp);
    let mut ends: BTreeMap<String, usize> = BTreeMap::new();
    for r in &rivers {
        let k2 = match r.end {
            RiverEnd::Sea => "sea",
            RiverEnd::Lake { .. } => "lake",
            RiverEnd::EndorheicLake { .. } => "endorheic lake",
            RiverEnd::Confluence { .. } => "confluence",
            RiverEnd::Terminal => "terminal",
        };
        *ends.entry(k2.to_string()).or_insert(0) += 1;
    }
    let mut strahler: BTreeMap<String, usize> = BTreeMap::new();
    for r in &rivers {
        *strahler.entry(r.strahler.map_or("spillway".to_string(), |s| s.to_string())).or_insert(0) += 1;
    }
    let p90_len = pct(&lens, 0.9);
    let p50_area = pct(&areas, 0.5);
    let long_small = rivers.iter().filter(|r| r.length_km >= p90_len && r.catchment_km2 < 10.0).count();
    eprintln!(
        "\n   A · length km p10/p50/p90/max {:.2}/{:.2}/{:.2}/{:.1} · catchment km² p10/p50/p90/max {:.2}/{:.2}/{:.1}/{:.0} · ends {ends:?} · Strahler {strahler:?} · rivers crossing a lake {} · long (≥ p90 length {p90_len:.1} km) and small (< 10 km²): **{long_small}** · median catchment {p50_area:.2} km²",
        pct(&lens, 0.1),
        pct(&lens, 0.5),
        p90_len,
        pct(&lens, 1.0),
        pct(&areas, 0.1),
        pct(&areas, 0.5),
        pct(&areas, 0.9),
        pct(&areas, 1.0),
        rivers.iter().filter(|r| !r.lakes_crossed.is_empty()).count()
    );
    for (thr_l, thr_a) in [(5.0f32, 0.0f32), (10.0, 0.0), (0.0, 10.0), (0.0, 100.0), (10.0, 10.0)] {
        let kept = rivers.iter().filter(|r| r.length_km >= thr_l && r.catchment_km2 >= thr_a).count();
        eprintln!("      e.g. length ≥ {thr_l} km and catchment ≥ {thr_a} km²: {kept} rivers");
    }
    let _ = n;
    eprintln!("\n==========  end Finding 146 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 146-H — the hex edges on Living Landz's grid (the author, 2026-10-05: 40 m, flat-top, axial, origin at
/// the bottom left), on the témoin: their count, the neighbour check, the attributes' ranges, the cost and the size.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f146_hex --nocapture
#[test]
#[ignore]
fn f146_hex() {
    use common::{build_world, viz_hd_lakes_on};
    use ymir_core::export::rivers_ll::{HexGrid, RiversLlParams, build_rivers_ll, river_hex_edges, rivers_ll_json};
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 146-H . the hex edges (40 m flat-top, axial, origin bottom-left)  ==========");
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let kn = Knobs { valley: Some(temoin), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let (rivers, _) = build_rivers_ll(&v.drainage, &v.conditioned, &ss, cell_km2, RiversLlParams::default(), false);
    let grid = HexGrid::living_landz();
    let te = Instant::now();
    let e = river_hex_edges(&rivers, &v.drainage, &v.conditioned, &ss, cell_km2, &grid, 200.0 / 10f32.powf(0.3), 0.3);
    let t_edges = te.elapsed().as_secs_f64();
    // since Finding 147 the edges are not exported: the size is the rivers' plus the edges' own JSON
    let bytes = [rivers_ll_json(&rivers, CELL_KM, &RiversLlParams::default()), serde_json::to_vec(&e).unwrap()].concat();
    let mut hex_count = vec![0u32; rivers.len()];
    for &r in &e.river_id {
        hex_count[r as usize] += 1;
    }
    let n = e.q1.len();
    let not_nb = (0..n)
        .filter(|&i| {
            let (dq, dr) = (e.q2[i] - e.q1[i], e.r2[i] - e.r1[i]);
            let ds = -dq - dr;
            (dq.abs() + dr.abs() + ds.abs()) / 2 != 1
        })
        .count();
    let pct = |v: &[f32], p: f64| -> f32 {
        let mut s = v.to_vec();
        s.sort_by(f32::total_cmp);
        if s.is_empty() { f32::NAN } else { s[((s.len() - 1) as f64 * p) as usize] }
    };
    let per: Vec<f32> = hex_count.iter().map(|&c| c as f32).collect();
    eprintln!(
        "   **{n} edge crossings** by {} rivers · not between neighbours: {not_nb} · per river p50 / p90 / max {:.0} / {:.0} / {:.0} · edges per km of river {:.1}",
        rivers.len(),
        pct(&per, 0.5),
        pct(&per, 0.9),
        pct(&per, 1.0),
        n as f32 / rivers.iter().map(|r| r.length_km).sum::<f32>().max(1e-6)
    );
    eprintln!(
        "   at the crossings · catchment km² p10/p50/p90 {:.1}/{:.1}/{:.0} · discharge m³/s p10/p50/p90 {:.2}/{:.2}/{:.1} · slope p10/p50/p90 {:.4}/{:.4}/{:.4} (negative {:.1} %) · width m p10/p50/p90 {:.0}/{:.0}/{:.0}",
        pct(&e.catchment_km2, 0.1),
        pct(&e.catchment_km2, 0.5),
        pct(&e.catchment_km2, 0.9),
        pct(&e.discharge_m3s, 0.1),
        pct(&e.discharge_m3s, 0.5),
        pct(&e.discharge_m3s, 0.9),
        pct(&e.slope, 0.1),
        pct(&e.slope, 0.5),
        pct(&e.slope, 0.9),
        100.0 * e.slope.iter().filter(|&&s| s < 0.0).count() as f64 / n.max(1) as f64,
        pct(&e.width_m, 0.1),
        pct(&e.width_m, 0.5),
        pct(&e.width_m, 0.9)
    );
    for (thr_l, thr_a) in [(0.0f32, 0.0f32), (5.0, 0.0), (10.0, 0.0), (0.0, 100.0), (10.0, 10.0)] {
        let kept: Vec<_> = rivers.iter().filter(|r| r.length_km >= thr_l && r.catchment_km2 >= thr_a).collect();
        eprintln!("      length ≥ {thr_l} km and catchment ≥ {thr_a} km²: {} rivers · {} hex edges", kept.len(), kept.iter().map(|r| hex_count[r.id as usize]).sum::<u32>());
    }
    eprintln!("   COST · hex edges {t_edges:.2} s per world · rivers_ll.json with the edges {:.1} MB", bytes.len() as f64 / 1e6);
    eprintln!("\n==========  end Finding 146-H . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 147-X, diagnosis before any fix — what the crossings between two smoothed rivers (F146: 662) are: a
/// spillway involved, two rivers sharing trace cells (a retrace), a D8 "X" (two diagonal steps across one 2×2 block), near
/// a river's end, or other. The témoin, F146's build.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f147_xdiag --nocapture
#[test]
#[ignore]
fn f147_xdiag() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use ymir_core::export::rivers_ll::{RiversLlParams, build_rivers_ll};
    use ymir_core::tectonics_c1::drainage::SegmentKind;
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let kn = Knobs { valley: Some(temoin), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let dr = &v.drainage;
    let (rivers, _) = build_rivers_ll(dr, &v.conditioned, &ss, cell_km2, RiversLlParams::default(), false);
    let w = v.conditioned.width;
    eprintln!("\n==========  Finding 147-X diagnosis . the crossings' classes  ==========");
    let cells_of: Vec<HashSet<(u32, u32)>> = rivers.iter().map(|r| r.segments.iter().flat_map(|&s| dr.rivers.segments[s].points.iter().copied()).collect()).collect();
    let seg_inter = |a: [f32; 2], b: [f32; 2], c: [f32; 2], d: [f32; 2]| -> Option<[f32; 2]> {
        let r = [b[0] - a[0], b[1] - a[1]];
        let s = [d[0] - c[0], d[1] - c[1]];
        let den = r[0] * s[1] - r[1] * s[0];
        if den.abs() < 1e-9 {
            return None;
        }
        let t = ((c[0] - a[0]) * s[1] - (c[1] - a[1]) * s[0]) / den;
        let u = ((c[0] - a[0]) * r[1] - (c[1] - a[1]) * r[0]) / den;
        ((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u)).then(|| [a[0] + t * r[0], a[1] + t * r[1]])
    };
    let polys: Vec<&Vec<[f32; 2]>> = rivers.iter().map(|r| &r.points).collect();
    let mut hash: HashMap<(i32, i32), Vec<(usize, usize)>> = HashMap::new();
    for (ri, p) in polys.iter().enumerate() {
        for e in 0..p.len().saturating_sub(1) {
            let (a, b) = (p[e], p[e + 1]);
            for gx in (a[0].min(b[0]) / 4.0).floor() as i32..=(a[0].max(b[0]) / 4.0).floor() as i32 {
                for gy in (a[1].min(b[1]) / 4.0).floor() as i32..=(a[1].max(b[1]) / 4.0).floor() as i32 {
                    hash.entry((gx, gy)).or_default().push((ri, e));
                }
            }
        }
    }
    let ends: Vec<([f32; 2], [f32; 2])> = polys.iter().map(|p| (p[0], *p.last().unwrap())).collect();
    let dist_end = |q: [f32; 2], ri: usize| -> f32 {
        let (a, b) = ends[ri];
        (((q[0] - a[0]).powi(2) + (q[1] - a[1]).powi(2)).sqrt()).min(((q[0] - b[0]).powi(2) + (q[1] - b[1]).powi(2)).sqrt())
    };
    // a D8 X: in the 2×2 block holding q, both diagonals are trace steps, one of each river
    let is_x = |q: [f32; 2], ra: usize, rb: usize| -> bool {
        let (bx, by) = ((q[0] - 0.5).floor() as i64, (q[1] - 0.5).floor() as i64);
        let c = |x: i64, y: i64| (x.max(0) as u32, y.max(0) as u32);
        let (p00, p11, p10, p01) = (c(bx, by), c(bx + 1, by + 1), c(bx + 1, by), c(bx, by + 1));
        let diag_a = cells_of[ra].contains(&p00) && cells_of[ra].contains(&p11);
        let anti_b = cells_of[rb].contains(&p10) && cells_of[rb].contains(&p01);
        let diag_b = cells_of[rb].contains(&p00) && cells_of[rb].contains(&p11);
        let anti_a = cells_of[ra].contains(&p10) && cells_of[ra].contains(&p01);
        (diag_a && anti_b) || (diag_b && anti_a)
    };
    let mut seen: HashSet<(usize, usize, usize, usize)> = HashSet::new();
    let mut cls: BTreeMap<&str, usize> = BTreeMap::new();
    let mut ex: Vec<String> = Vec::new();
    for list in hash.values() {
        for i in 0..list.len() {
            for j in i + 1..list.len() {
                let ((ra, ea), (rb, eb)) = (list[i], list[j]);
                if ra == rb {
                    continue;
                }
                let key = if (ra, ea) < (rb, eb) { (ra, ea, rb, eb) } else { (rb, eb, ra, ea) };
                if !seen.insert(key) {
                    continue;
                }
                let Some(q) = seg_inter(polys[ra][ea], polys[ra][ea + 1], polys[rb][eb], polys[rb][eb + 1]) else { continue };
                if dist_end(q, ra) < 0.05 || dist_end(q, rb) < 0.05 {
                    continue;
                }
                let spill = rivers[ra].kind == SegmentKind::Spillway || rivers[rb].kind == SegmentKind::Spillway;
                let shared = cells_of[ra].intersection(&cells_of[rb]).count();
                let c = if spill && shared > 0 {
                    "spillway sharing trace cells (a retrace)"
                } else if spill {
                    "spillway, no shared cell"
                } else if shared > 0 {
                    "two watercourses sharing trace cells"
                } else if is_x(q, ra, rb) {
                    "a D8 X (two diagonals across one 2×2 block)"
                } else if dist_end(q, ra).min(dist_end(q, rb)) < 3.0 {
                    "within 3 cells of a river's end"
                } else {
                    "other"
                };
                *cls.entry(c).or_insert(0) += 1;
                if ex.len() < 12 && (c == "other" || c.starts_with("a D8")) {
                    ex.push(format!("{c}: rivers {ra} / {rb} at ({:.1},{:.1})", q[0], q[1]));
                }
            }
        }
    }
    let tot: usize = cls.values().sum();
    eprintln!("   crossings away from an end: {tot}");
    for (k, c) in &cls {
        eprintln!("      {k}: {c}");
    }
    for e in &ex {
        eprintln!("      e.g. {e}");
    }
    let _ = w;
}

/// ADR Finding 147 — the rivers for Living Landz, round 2 (format 0.3.0) on the témoin: the crossings resolved in the
/// exported geometry (X), the per-vertex attributes and the bed-width table (V), the annotations (F), the valley-floor
/// instrument on the final polyline (T) and the cost. Declared in `f147_declared.md`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f147_rivers --nocapture
#[test]
#[ignore]
fn f147_rivers() {
    use common::{build_world, viz_hd_lakes_on};
    use std::collections::{BTreeMap, HashMap, HashSet};
    use ymir_core::export::rivers_ll::{RiverEnd, RiversLlParams, build_rivers_ll, rivers_ll_json};
    use ymir_core::export::hydro::rivers_json;
    use ymir_core::tectonics_c1::drainage::SegmentKind;
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 147 . the rivers for Living Landz, format 0.3.0 (témoin C2 /10 col)  ==========");
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let kn = Knobs { valley: Some(temoin), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let dr = &v.drainage;
    let rj0 = rivers_json(dr, cell_km2);
    let params = RiversLlParams::default();
    // production cost (no instruments), twice (the first warms the pool), then the instrumented build
    let mut t_build = 0f64;
    let mut rivers = Vec::new();
    for _ in 0..2 {
        let tp = Instant::now();
        rivers = build_rivers_ll(dr, &v.conditioned, &ss, cell_km2, params, false).0;
        t_build = tp.elapsed().as_secs_f64();
    }
    let tj = Instant::now();
    let bytes = rivers_ll_json(&rivers, CELL_KM, &params);
    let t_json = tj.elapsed().as_secs_f64();
    let (rivers_s, st) = build_rivers_ll(dr, &v.conditioned, &ss, cell_km2, params, true);
    let same = rivers.len() == rivers_s.len() && rivers.iter().zip(&rivers_s).all(|(a, b)| a.points == b.points && a.vertex.discharge_m3s == b.vertex.discharge_m3s);
    let rj1 = rivers_json(dr, cell_km2);
    let nv: usize = rivers.iter().map(|r| r.points.len()).sum();
    let aligned = rivers.iter().all(|r| {
        let n = r.points.len();
        r.vertex.catchment_km2.len() == n && r.vertex.discharge_m3s.len() == n && r.vertex.slope.len() == n && r.vertex.valley_width_m.len() == n && r.vertex.bed_width_m.len() == n
    });
    eprintln!(
        "   {} segments → **{} rivers** · {nv} vertices · per-vertex arrays aligned: {aligned} · instrumented build identical: {same} · rivers.json byte-identical: {}",
        dr.rivers.segments.len(),
        rivers.len(),
        rj0 == rj1
    );
    eprintln!("   V · **rivers_ll.json {:.2} MB** (format 0.3.0, no hex edge)", bytes.len() as f64 / 1e6);
    eprintln!(
        "   X · reaches joining their receiver below its head (T3b), now ending in a confluence: {} · of them larger than the receiver's reach above the junction (Hack's convention broken there): {}",
        st.mid_joins, st.mid_joins_larger
    );
    eprintln!("   COST · build {t_build:.2} s + serialize {t_json:.2} s per world (against run_hd 249.8 s)");
    let pct = |v: &[f32], p: f64| -> f32 {
        let mut s = v.to_vec();
        s.sort_by(f32::total_cmp);
        if s.is_empty() { f32::NAN } else { s[((s.len() - 1) as f64 * p) as usize] }
    };
    // ── T: the valley-floor instrument on the final polyline
    eprintln!(
        "\n   T · lateral deviation from the D8 trace, final polyline: p50 {:.3} · p99 **{:.3}** · max {:.3} cell ({} samples every 0.25 cell) · Chaikin before the constraint p99 {:.3} · max {:.3}",
        pct(&st.lateral_cells, 0.5),
        pct(&st.lateral_cells, 0.99),
        pct(&st.lateral_cells, 1.0),
        st.lateral_cells.len(),
        pct(&st.lateral_before, 0.99),
        pct(&st.lateral_before, 1.0)
    );
    eprintln!(
        "   T · ridge crossings: raw D8 trace {} of {} samples · final **{}** of {} · vertices pulled back by the 0.5-cell bound {} of {} · stair index (turns ≥ 45°): raw {:.3} → final {:.3}",
        st.ridge_cross_raw,
        st.samples_raw,
        st.ridge_cross_smoothed,
        st.samples,
        st.pulled,
        st.vertices,
        st.stair_raw,
        st.stair_smoothed
    );
    {
        // the ridge-crossing samples, diagnosed: on the river's own trace cell or not, by river kind, on an exact half
        // step of the D8 trace (a corner between two trace cells) or not
        let mut cls: BTreeMap<String, usize> = BTreeMap::new();
        for &(ri, q, own) in &st.ridge_samples {
            let r = &rivers_s[ri as usize];
            let corner = ((q[0] - q[0].round()).abs() < 1e-3) && ((q[1] - q[1].round()).abs() < 1e-3);
            let k = format!("{:?} · {} · {}", r.kind, if own { "its own trace cell" } else { "OFF its trace" }, if corner { "at a cell corner" } else { "not at a corner" });
            *cls.entry(k).or_insert(0) += 1;
        }
        eprintln!("   T · the {} final ridge-crossing samples by class: {cls:?}", st.ridge_samples.len());
        let dc: Vec<f32> = st.ridge_samples.iter().map(|&(_, q, _)| ((q[0] - q[0].round()).powi(2) + (q[1] - q[1].round()).powi(2)).sqrt()).collect();
        eprintln!(
            "   T · their distance to the nearest cell corner: max {:.4} cell · within 0.02 cell: {} of {}",
            dc.iter().copied().fold(0f32, f32::max),
            dc.iter().filter(|&&d| d <= 0.02).count(),
            dc.len()
        );
        for &(ri, q, own) in st.ridge_samples.iter().take(8) {
            let r = &rivers_s[ri as usize];
            eprintln!("      e.g. river {ri} ({:?}, end {:?}, {} lakes crossed) at ({:.3},{:.3}) own cell {own}", r.kind, r.end, r.lakes_crossed.len(), q[0], q[1]);
        }
    }
    let out1 = |v: &[f32]| 100.0 * v.iter().filter(|&&x| x > 1.0).count() as f64 / v.len().max(1) as f64;
    eprintln!(
        "   T · (F146's retired instrument, for the record) terrain − bed > 1 m: raw {:.2} % · Chaikin {:.2} % · final {:.2} %",
        out1(&st.above_bed_raw),
        out1(&st.above_bed_before),
        out1(&st.above_bed_after)
    );
    // ── X: the crossings on the final geometry
    let polys: Vec<&Vec<[f32; 2]>> = rivers.iter().map(|r| &r.points).collect();
    let seg_inter = |a: [f32; 2], b: [f32; 2], c: [f32; 2], d: [f32; 2]| -> Option<[f32; 2]> {
        let r = [b[0] - a[0], b[1] - a[1]];
        let s = [d[0] - c[0], d[1] - c[1]];
        let den = r[0] * s[1] - r[1] * s[0];
        if den.abs() < 1e-9 {
            return None;
        }
        let t = ((c[0] - a[0]) * s[1] - (c[1] - a[1]) * s[0]) / den;
        let u = ((c[0] - a[0]) * r[1] - (c[1] - a[1]) * r[0]) / den;
        ((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u)).then(|| [a[0] + t * r[0], a[1] + t * r[1]])
    };
    let mut hash: HashMap<(i32, i32), Vec<(usize, usize)>> = HashMap::new();
    for (ri, p) in polys.iter().enumerate() {
        for e in 0..p.len().saturating_sub(1) {
            let (a, b) = (p[e], p[e + 1]);
            for gx in (a[0].min(b[0]) / 4.0).floor() as i32..=(a[0].max(b[0]) / 4.0).floor() as i32 {
                for gy in (a[1].min(b[1]) / 4.0).floor() as i32..=(a[1].max(b[1]) / 4.0).floor() as i32 {
                    hash.entry((gx, gy)).or_default().push((ri, e));
                }
            }
        }
    }
    let d2 = |p: [f32; 2], q: [f32; 2]| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
    let near_end = |q: [f32; 2], ri: usize| -> bool { d2(q, polys[ri][0]) < 0.05 || d2(q, *polys[ri].last().unwrap()) < 0.05 };
    let cells_of: Vec<HashSet<(u32, u32)>> = rivers.iter().map(|r| r.segments.iter().flat_map(|&s| dr.rivers.segments[s].points.iter().copied()).collect()).collect();
    let mut seen: HashSet<(usize, usize, usize, usize)> = HashSet::new();
    let (mut f146_count, mut selfx) = (0usize, 0usize);
    let mut strict: Vec<(usize, usize, [f32; 2])> = Vec::new();
    let mut strict_seen: HashSet<(usize, usize, i64, i64)> = HashSet::new();
    for list in hash.values() {
        for i in 0..list.len() {
            for j in i + 1..list.len() {
                let ((ra, ea), (rb, eb)) = (list[i], list[j]);
                let key = if (ra, ea) < (rb, eb) { (ra, ea, rb, eb) } else { (rb, eb, ra, ea) };
                if !seen.insert(key) {
                    continue;
                }
                if ra == rb && ea.abs_diff(eb) <= 1 {
                    continue;
                }
                let (pa, pb) = (polys[ra], polys[rb]);
                let Some(q) = seg_inter(pa[ea], pa[ea + 1], pb[eb], pb[eb + 1]) else { continue };
                if ra == rb {
                    selfx += 1;
                    continue;
                }
                if near_end(q, ra) || near_end(q, rb) {
                    continue;
                }
                f146_count += 1;
                // strict: not a touch at a vertex both polylines share
                let at_a = d2(q, pa[ea]) < 1e-4 || d2(q, pa[ea + 1]) < 1e-4;
                let at_b = d2(q, pb[eb]) < 1e-4 || d2(q, pb[eb + 1]) < 1e-4;
                if at_a && at_b {
                    continue;
                }
                let (lo, hi) = (ra.min(rb), ra.max(rb));
                if strict_seen.insert((lo, hi, (q[0] * 1000.0).round() as i64, (q[1] * 1000.0).round() as i64)) {
                    strict.push((lo, hi, q));
                }
            }
        }
    }
    // crossings through a shared run of vertices (a polyline arriving on one side of the other and leaving on the other)
    let mut vmap: HashMap<(u32, u32), Vec<(usize, usize)>> = HashMap::new();
    for (ri, p) in polys.iter().enumerate() {
        for (vi, q) in p.iter().enumerate() {
            vmap.entry((q[0].to_bits(), q[1].to_bits())).or_default().push((ri, vi));
        }
    }
    let mut share_pairs: HashSet<(usize, usize)> = HashSet::new();
    for l in vmap.values() {
        for &(a, _) in l {
            for &(b, _) in l {
                if a < b {
                    share_pairs.insert((a, b));
                }
            }
        }
    }
    let side = |o: [f32; 2], d: [f32; 2], p: [f32; 2]| -> i32 {
        let c = d[0] * (p[1] - o[1]) - d[1] * (p[0] - o[0]);
        if c > 1e-6 { 1 } else if c < -1e-6 { -1 } else { 0 }
    };
    let (mut overlap_cross, mut overlap_touch, mut shared_runs) = (0usize, 0usize, 0usize);
    let mut sp: Vec<(usize, usize)> = share_pairs.into_iter().collect();
    sp.sort();
    for &(a, b) in &sp {
        let (pa, pb) = (polys[a], polys[b]);
        let idx_b: HashMap<(u32, u32), usize> = pb.iter().enumerate().map(|(i, q)| ((q[0].to_bits(), q[1].to_bits()), i)).collect();
        let kb = |q: [f32; 2]| idx_b.get(&(q[0].to_bits(), q[1].to_bits())).copied();
        let mut i = 0usize;
        while i < pa.len() {
            let Some(k0) = kb(pa[i]) else {
                i += 1;
                continue;
            };
            // the maximal run of consecutive shared vertices from i
            let (mut i1, mut k1) = (i, k0);
            while i1 + 1 < pa.len() {
                match kb(pa[i1 + 1]) {
                    Some(k) if k.abs_diff(k1) == 1 => {
                        i1 += 1;
                        k1 = k;
                    }
                    _ => break,
                }
            }
            shared_runs += 1;
            let interior_a = i > 0 && i1 + 1 < pa.len();
            let interior_b = k0.min(k1) > 0 && k0.max(k1) + 1 < pb.len();
            if interior_a && interior_b {
                let tan = |k: usize| [pb[k + 1][0] - pb[k - 1][0], pb[k + 1][1] - pb[k - 1][1]];
                let sb = side(pb[k0], tan(k0), pa[i - 1]);
                let sa = side(pb[k1], tan(k1), pa[i1 + 1]);
                if sb != 0 && sa != 0 && sb != sa {
                    overlap_cross += 1;
                    if strict.len() < 100_000 {
                        strict.push((a, b, pa[i]));
                    }
                } else {
                    overlap_touch += 1;
                }
            }
            i = i1 + 1;
        }
    }
    let n_proper = strict.len() - overlap_cross;
    eprintln!(
        "\n   X · F146's instrument (any intersection away from a river's end): {f146_count} · of which not a touch at a shared vertex: **{n_proper}** · shared vertex runs {shared_runs}: crossing through the run **{overlap_cross}**, touching {overlap_touch} · self-intersections **{selfx}**"
    );
    let mut cls: BTreeMap<String, usize> = BTreeMap::new();
    let mut ex: Vec<String> = Vec::new();
    for &(ra, rb, q) in &strict {
        let spill = rivers[ra].kind == SegmentKind::Spillway || rivers[rb].kind == SegmentKind::Spillway;
        let shared = cells_of[ra].intersection(&cells_of[rb]).count();
        let joined = matches!(rivers[ra].end, RiverEnd::Confluence { river_id } if river_id as usize == rb) || matches!(rivers[rb].end, RiverEnd::Confluence { river_id } if river_id as usize == ra);
        let c = format!(
            "{} · {} · {}",
            if spill { "a spillway involved" } else { "two watercourses" },
            if shared > 0 { "sharing trace cells" } else { "no shared cell" },
            if joined { "a confluence pair" } else { "not joined" }
        );
        *cls.entry(c.clone()).or_insert(0) += 1;
        if ex.len() < 12 {
            ex.push(format!("{c}: rivers {ra} ({:?}, {} pts) / {rb} ({:?}, {} pts) at ({:.2},{:.2})", rivers[ra].kind, polys[ra].len(), rivers[rb].kind, polys[rb].len(), q[0], q[1]));
        }
    }
    eprintln!("   X · the crossings left (proper + through a run), by class: {cls:?}");
    for e in &ex {
        eprintln!("      e.g. {e}");
    }
    // the local picture of the first crossings: both rivers' D8 trace cells and final vertices within 3 cells
    for &(ra, rb, q) in strict.iter().take(6) {
        eprintln!("      ── at ({:.2},{:.2}): rivers {ra} (end {:?}) / {rb} (end {:?})", q[0], q[1], rivers[ra].end, rivers[rb].end);
        for ri in [ra, rb] {
            let mut tr: Vec<(u32, u32)> = Vec::new();
            for &sg in &rivers[ri].segments {
                for &c in &dr.rivers.segments[sg].points {
                    if tr.last() != Some(&c) {
                        tr.push(c);
                    }
                }
            }
            let near_c: Vec<String> = tr
                .iter()
                .enumerate()
                .filter(|(_, c)| (c.0 as f32 + 0.5 - q[0]).abs() <= 3.0 && (c.1 as f32 + 0.5 - q[1]).abs() <= 3.0)
                .map(|(i, c)| format!("{i}:({},{})", c.0, c.1))
                .collect();
            let np = polys[ri].len();
            let near_v: Vec<String> = polys[ri]
                .iter()
                .enumerate()
                .filter(|(_, p)| (p[0] - q[0]).abs() <= 3.0 && (p[1] - q[1]).abs() <= 3.0)
                .map(|(i, p)| format!("{i}/{np}:({:.2},{:.2})", p[0], p[1]))
                .collect();
            eprintln!("         river {ri}: trace cells (of {}) {near_c:?}", tr.len());
            eprintln!("         river {ri}: vertices {near_v:?}");
        }
    }
    eprintln!(
        "   X · spillway trace points spliced onto a watercourse's vertices: {} of {} ({:.1} %)",
        st.spliced_points,
        st.spillway_points,
        100.0 * st.spliced_points as f64 / st.spillway_points.max(1) as f64
    );
    // ── V: the discharge, accumulated against the prorata
    let (mut n_v, mut n_gap, mut n_zero) = (0usize, 0usize, 0usize);
    let mut ratio: Vec<f32> = Vec::new();
    let mut by_kind: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for r in &rivers {
        for i in 0..r.points.len() {
            let qa = r.vertex.discharge_m3s[i];
            let qp = if r.catchment_km2 > 0.0 { r.discharge_m3s * r.vertex.catchment_km2[i] / r.catchment_km2 } else { 0.0 };
            n_v += 1;
            if qa <= 0.0 {
                n_zero += 1;
                continue;
            }
            ratio.push(qp / qa);
            let gap = (qp - qa).abs() / qa > 0.2;
            if gap {
                n_gap += 1;
            }
            let e = by_kind.entry(format!("{:?}", r.kind)).or_insert((0, 0));
            e.0 += 1;
            e.1 += gap as usize;
        }
    }
    eprintln!(
        "\n   V · discharge: accumulated (segment_discharge_profile_m3s) against the prorata (mouth Q × A/A_mouth): **{:.1} %** of the {n_v} vertices differ by > 20 % · prorata/accumulated p10/p50/p90 {:.2}/{:.2}/{:.2} · vertices with no accumulated discharge {n_zero} · by kind {:?}",
        100.0 * n_gap as f64 / (n_v - n_zero).max(1) as f64,
        pct(&ratio, 0.1),
        pct(&ratio, 0.5),
        pct(&ratio, 0.9),
        by_kind.iter().map(|(k, (n, g))| format!("{k}: {:.1} % of {n}", 100.0 * *g as f64 / (*n).max(1) as f64)).collect::<Vec<_>>()
    );
    let all = |f: &dyn Fn(&ymir_core::export::rivers_ll::VertexAttrs) -> &Vec<f32>| -> Vec<f32> { rivers.iter().flat_map(|r| f(&r.vertex).iter().copied()).collect() };
    let (ca, qa, sl, vw, bw) = (all(&|v| &v.catchment_km2), all(&|v| &v.discharge_m3s), all(&|v| &v.slope), all(&|v| &v.valley_width_m), all(&|v| &v.bed_width_m));
    eprintln!(
        "   V · per vertex p10/p50/p90/max · catchment km² {:.2}/{:.2}/{:.1}/{:.0} · discharge m³/s {:.3}/{:.3}/{:.2}/{:.0} · slope {:.4}/{:.4}/{:.4} (negative {:.1} %) · valley width m {:.0}/{:.0}/{:.0}/{:.0} · bed width m (a = 5) {:.1}/{:.1}/{:.1}/{:.0}",
        pct(&ca, 0.1),
        pct(&ca, 0.5),
        pct(&ca, 0.9),
        pct(&ca, 1.0),
        pct(&qa, 0.1),
        pct(&qa, 0.5),
        pct(&qa, 0.9),
        pct(&qa, 1.0),
        pct(&sl, 0.1),
        pct(&sl, 0.5),
        pct(&sl, 0.9),
        100.0 * sl.iter().filter(|&&s| s < 0.0).count() as f64 / sl.len().max(1) as f64,
        pct(&vw, 0.1),
        pct(&vw, 0.5),
        pct(&vw, 0.9),
        pct(&vw, 1.0),
        pct(&bw, 0.1),
        pct(&bw, 0.5),
        pct(&bw, 0.9),
        pct(&bw, 1.0)
    );
    // ── V: the bed-width table (the share of river km wider than one hex, both radius readings)
    let flat_corner = 40.0 * 3f32.sqrt(); // 69.28 m: flat to flat with a 40 m centre-to-corner radius
    let flat_edge = 80.0; // 80 m: flat to flat with a 40 m centre-to-edge radius
    for (label, keep) in [("all rivers", None), ("watercourses only (no spillway)", Some(SegmentKind::Watercourse))] {
        let mut tot = 0f64;
        let mut wide: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
        let coefs = [2.5f32, 3.5, 5.0, 7.0];
        for r in rivers.iter().filter(|r| keep.is_none_or(|k| r.kind == k)) {
            for e in 0..r.points.len().saturating_sub(1) {
                let l = d2(r.points[e], r.points[e + 1]) as f64 * CELL_KM as f64;
                let q = 0.5 * (r.vertex.discharge_m3s[e] + r.vertex.discharge_m3s[e + 1]);
                tot += l;
                for &a in &coefs {
                    let wb = a * q.max(0.0).sqrt();
                    let ent = wide.entry((a * 10.0) as u32).or_insert((0.0, 0.0));
                    if wb > flat_corner {
                        ent.0 += l;
                    }
                    if wb > flat_edge {
                        ent.1 += l;
                    }
                }
            }
        }
        eprintln!("   V · bed width w = a·Q^0.5, {label} ({tot:.0} km): share of km wider than 69.3 m (radius at the corner) / 80 m (radius at the edge)");
        for (a10, (k1, k2)) in &wide {
            let a = *a10 as f32 / 10.0;
            eprintln!(
                "      a = {a:.1}: **{:.2} %** / **{:.2} %** ({:.1} / {:.1} km) · needs Q > {:.0} / {:.0} m³/s",
                100.0 * k1 / tot.max(1e-9),
                100.0 * k2 / tot.max(1e-9),
                k1,
                k2,
                (flat_corner / a).powi(2),
                (flat_edge / a).powi(2)
            );
        }
    }
    // ── F: the annotations
    let n_pairs = st.parallel_pairs.len();
    let n_cand = rivers.iter().filter(|r| r.fusion_candidate).count();
    let n_par = rivers.iter().filter(|r| r.parallel_of.is_some()).count();
    let n_retr = rivers.iter().filter(|r| r.retraces_spillway).count();
    let n_spill = rivers.iter().filter(|r| r.kind == SegmentKind::Spillway).count();
    let km_cand: f32 = rivers.iter().filter(|r| r.fusion_candidate).map(|r| r.length_km).sum();
    let km_retr: f32 = rivers.iter().filter(|r| r.retraces_spillway).map(|r| r.length_km).sum();
    eprintln!(
        "\n   F · parallel pairs (≤ 2 cells over > 2 km, no shared confluence, no spillway with the watercourse it retraces): **{n_pairs}** · rivers with a parallel_of {n_par} · fusion candidates **{n_cand}** ({km_cand:.0} km) · retracing spillways **{n_retr}** of {n_spill} spillways ({km_retr:.0} km)"
    );
    let hidden = rivers.iter().filter(|r| r.fusion_candidate || r.retraces_spillway).count();
    eprintln!("   F · both filters on: {} rivers left of {} · {:.0} km of {:.0} km", rivers.len() - hidden, rivers.len(), rivers.iter().filter(|r| !(r.fusion_candidate || r.retraces_spillway)).map(|r| r.length_km).sum::<f32>(), rivers.iter().map(|r| r.length_km).sum::<f32>());
    eprintln!("\n==========  end Finding 147 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 147-R — the hex reference for Living Landz: 5 varied rivers of the témoin, their polylines and the hex
/// edges F146's computation expects on the assumed grid, written to `docs/rivers_ll_hex_reference.json`.
///
/// Run: cargo test -p ymir-core --release --test f126_coast -- --ignored --exact f147_hexref --nocapture
#[test]
#[ignore]
fn f147_hexref() {
    use common::{build_world, viz_hd_lakes_on};
    use ymir_core::export::rivers_ll::{HexGrid, RiverEnd, RiversLlParams, build_rivers_ll, river_hex_edges};
    use ymir_core::tectonics_c1::drainage::SegmentKind;
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let cell_m = CELL_KM * 1000.0;
    let temoin = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let kn = Knobs { valley: Some(temoin), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn, None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
    let (rivers, _) = build_rivers_ll(&v.drainage, &v.conditioned, &ss, cell_km2, RiversLlParams::default(), false);
    // the 5 picks (declared): the first river by id in each class, not already picked
    let used = std::cell::RefCell::new(Vec::<usize>::new());
    let pick = |name: &str, f: &dyn Fn(&ymir_core::export::rivers_ll::RiverLl) -> bool| -> (String, usize) {
        let i = rivers.iter().position(|r| f(r) && !used.borrow().contains(&(r.id as usize))).unwrap_or_else(|| panic!("no river for {name}"));
        used.borrow_mut().push(i);
        (name.to_string(), i)
    };
    let picks = [
        pick("a watercourse ending at the sea, 2–5 km", &|r| r.kind == SegmentKind::Watercourse && r.end == RiverEnd::Sea && (2.0..5.0).contains(&r.length_km)),
        pick("a tributary ending at a confluence, 1–3 km", &|r| matches!(r.end, RiverEnd::Confluence { .. }) && (1.0..3.0).contains(&r.length_km)),
        pick("a river ending in a lake, 1–5 km", &|r| matches!(r.end, RiverEnd::Lake { .. } | RiverEnd::EndorheicLake { .. }) && (1.0..5.0).contains(&r.length_km)),
        pick("a spillway, 0.5–5 km", &|r| r.kind == SegmentKind::Spillway && (0.5..5.0).contains(&r.length_km)),
        pick("a large river (catchment ≥ 50 km²) under 8 km", &|r| r.catchment_km2 >= 50.0 && r.length_km < 8.0),
    ];
    let grid = HexGrid::living_landz();
    let mut out_rivers = Vec::new();
    for (name, i) in &picks {
        let r = &rivers[*i];
        let e = river_hex_edges(std::slice::from_ref(r), &v.drainage, &v.conditioned, &ss, cell_km2, &grid, 200.0 / 10f32.powf(0.3), 0.3);
        let edges: Vec<serde_json::Value> = (0..e.q1.len()).map(|k| serde_json::json!([[e.q1[k], e.r1[k]], [e.q2[k], e.r2[k]]])).collect();
        let pts_m: Vec<[f32; 2]> = r.points.iter().map(|p| [p[0] * cell_m, p[1] * cell_m]).collect();
        let hex_first = grid.hex_of(pts_m[0][0], pts_m[0][1]);
        let hex_last = grid.hex_of(pts_m.last().unwrap()[0], pts_m.last().unwrap()[1]);
        eprintln!("   {name}: river {} · {:.2} km · {} vertices · {} hex edges", r.id, r.length_km, r.points.len(), edges.len());
        out_rivers.push(serde_json::json!({
            "case": name,
            "river_id": r.id,
            "kind": r.kind,
            "length_km": r.length_km,
            "points_cells": r.points,
            "points_m": pts_m,
            "hex_of_first_point": [hex_first.0, hex_first.1],
            "hex_of_last_point": [hex_last.0, hex_last.1],
            "expected_hex_edges": edges,
        }));
    }
    let doc = serde_json::json!({
        "what": "ADR Finding 147-R — a reference for Living Landz's own hex computation from the rivers_ll.json polylines: 5 rivers of the témoin world, their polylines and the hex edges Ymir's reference computation (Finding 146-H) finds.",
        "assumptions_TO_CONFIRM_by_the_author": {
            "size_m": 40.0,
            "size_reading": "the 40 m radius is read CENTRE TO CORNER (a hex is 80 m corner to corner, 69.28 m flat to flat). If Living Landz's 40 m is centre to EDGE (80 m flat to flat), the edges below do not apply.",
            "orientation": "flat_top",
            "coordinates": "axial (q, r); cube rounding; the six neighbours of (q, r) are (q±1, r), (q, r±1), (q+1, r−1), (q−1, r+1)",
            "origin": "hex (0, 0) is CENTRED on the map's bottom-left corner (x = 0, y = 0). If Living Landz puts that corner at a hex's corner or edge instead, every hex index shifts.",
            "y_axis": "y grows NORTHWARD from the bottom edge (the container's y = 0 = south)."
        },
        "conversion": {
            "metres_per_cell": cell_m,
            "x_m": "x_cells × metres_per_cell",
            "y_m": "y_cells × metres_per_cell",
            "axial_from_metres": "q = (2/3 · x_m) / size_m ; r = (−x_m/3 + √3/3 · y_m) / size_m ; then cube rounding"
        },
        "edge_rule": "each polyline is walked every 0.1 × size_m (4 m); each change of hex between two consecutive samples is one edge crossing [[q1, r1], [q2, r2]], upstream side first, in order along the river",
        "rivers": out_rivers,
    });
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/rivers_ll_hex_reference.json");
    std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    eprintln!("   written {}", path.display());
}
