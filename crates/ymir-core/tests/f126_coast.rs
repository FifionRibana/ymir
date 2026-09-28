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
