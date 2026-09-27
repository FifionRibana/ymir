//! ADR 0001 Finding 124, part 3 — **does the export tell the truth?** Partition, catchments, length,
//! and the rectangle, on the C2/10 world (the definition: basin base) at the CANONICAL framing.
//!
//! The microscope's objects are `aggregate_watercourses` of the viz, PORTED (with Finding 123's
//! max-A trunk and its 1 % tie): an object is every segment draining to one terminal, lakes
//! chained. Areas and lengths are printed SIGNIFIED (× 7.5², × 7.5, what the viz shows) AND drawn.
//!
//! Run: cargo test -p ymir-core --release --test f124_export -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::collections::HashMap;
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, LakeType, SegmentKind,
    c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::valley_construction::{
    F121_AGE_K, ValleyConstruction, carve, skeleton,
};
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;

fn dcfg() -> C1DrainageConfig {
    let mut d = C1DrainageConfig::default();
    d.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    d.thresholds.full_tree = false;
    d
}

fn inventory(g: &GridF32, ss: &SteinSteinParams) -> (C1DrainageResult, GridF32) {
    let (w, h) = (g.width, g.height);
    let d = c1_drainage_windowed(g, None, &dcfg(), ss, DOMAIN_KM);
    let bre = breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h);
    let cl = c1_climate_placed(&bre, ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dc = DrainageClimate { precip_internal: &cl.precipitation, temperature: &cl.temperature };
    let dr =
        assemble_hd_drainage(&bre, &dc, Some(d), &dcfg(), ss, DOMAIN_KM, GEO_RATIO, None, false)
            .drainage;
    (dr, bre)
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Sink {
    Sea,
    ExoLake,
    EndoLake,
    SubSeaSink,
    Unknown,
}

/// PORT of `workspace.rs::classify_sink`.
fn classify_sink(
    dr: &C1DrainageResult,
    wc: &[u8],
    w: usize,
    h: usize,
    x: u32,
    y: u32,
) -> (Sink, Option<u32>) {
    let endo: std::collections::HashSet<u32> =
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).map(|l| l.base.id).collect();
    let mut sea = false;
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let (nx, ny) = (x as i32 + dx, y as i32 + dy);
            if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                let k = ny as usize * w + nx as usize;
                let id = dr.lake_map[k];
                if id != 0 {
                    return (
                        if endo.contains(&id) { Sink::EndoLake } else { Sink::ExoLake },
                        Some(id),
                    );
                }
                if wc.get(k).copied() == Some(1) {
                    sea = true;
                }
            }
        }
    }
    let mk = y as usize * w + x as usize;
    if !sea && wc.get(mk).copied() == Some(2) {
        return (Sink::SubSeaSink, None);
    }
    (if sea { Sink::Sea } else { Sink::Unknown }, None)
}

struct Wc {
    trunk: Vec<usize>,
    segments: Vec<usize>,
    order: u8,
    catchment_km2: f32,
    length_km: f32,
    sink: Sink,
    kind: SegmentKind,
    /// exorheic lakes the trunk crosses (Finding 45's chaining) — Finding 93's control
    lake_hops: usize,
    /// per trunk step (source → mouth): reached from below by a lake / body hop (dissection)
    hopped: Vec<bool>,
}

/// PORT of `workspace.rs::aggregate_watercourses_opt(.., chain_lakes = true)`, minus the fields
/// this bench does not read. `trunk_rule` = the shipped longest-path climb, or (C2) the max-A climb.
fn aggregate(dr: &C1DrainageResult, bre: &GridF32, max_a_trunk: bool) -> Vec<Wc> {
    let (w, h) = (bre.width, bre.height);
    let wc = ymir_core::lakes::connectivity::water_class(bre, 0.5);
    let segs = &dr.rivers.segments;
    let n = segs.len();
    let q = |i: usize| dr.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
    let mut outlet_of: HashMap<u32, usize> = HashMap::new();
    for (i, sg) in segs.iter().enumerate() {
        let (sx, sy) = sg.points[0];
        if let (Sink::ExoLake, Some(id)) = classify_sink(dr, &wc, w, h, sx, sy) {
            let better = outlet_of.get(&id).map_or(true, |&j| q(i) > q(j));
            if better {
                outlet_of.insert(id, i);
            }
        }
    }
    let across_lake = |i: usize| -> Option<usize> {
        let (mx, my) = *segs[i].points.last()?;
        match classify_sink(dr, &wc, w, h, mx, my) {
            (Sink::ExoLake, Some(id)) => outlet_of.get(&id).copied().filter(|&o| o != i),
            _ => None,
        }
    };
    let mut root = vec![0usize; n];
    for i in 0..n {
        let mut j = i;
        let mut seen = 0usize;
        loop {
            seen += 1;
            if seen > n {
                break;
            }
            match segs[j].downstream {
                Some(k) if k < n => j = k,
                _ => match across_lake(j) {
                    Some(o) => j = o,
                    None => break,
                },
            }
        }
        root[i] = j;
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        groups.entry(root[i]).or_default().push(i);
    }
    let km_per_cell = DOMAIN_KM / w as f32;
    let mut pathlen = vec![0u32; n];
    {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| segs[i].strahler_order);
        for _ in 0..16 {
            let mut changed = false;
            for &i in &order {
                let up = segs[i]
                    .upstream
                    .iter()
                    .copied()
                    .filter(|&u| u < n)
                    .map(|u| pathlen[u])
                    .max()
                    .unwrap_or(0);
                let v = up + segs[i].points.len() as u32;
                if v != pathlen[i] {
                    pathlen[i] = v;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    let area = |i: usize| dr.segment_drainage_km2.get(i).copied().unwrap_or(0.0);
    // ADR Finding 123-C2 -- max drained area; two areas within 1 % (PROXY) are a tie -- an
    // inherited area (F42/F93) or two branches of one catchment differing by a few cells -- and
    // the tie goes to the shipped longest-path rule
    let by_area = |a: usize, b: usize| -> std::cmp::Ordering {
        let (x, y) = (area(a), area(b));
        if (x - y).abs() <= 0.01 * x.max(y) { pathlen[a].cmp(&pathlen[b]) } else { x.total_cmp(&y) }
    };
    let inflow_to_lake_of = |i: usize| -> Option<usize> {
        let (sx, sy) = segs[i].points[0];
        let id = match classify_sink(dr, &wc, w, h, sx, sy) {
            (Sink::ExoLake, Some(id)) => id,
            _ => return None,
        };
        segs.iter()
            .enumerate()
            .filter(|&(j, sg)| {
                j != i
                    && sg.downstream.is_none()
                    && sg.points.last().is_some_and(|&(mx, my)| {
                        matches!(classify_sink(dr, &wc, w, h, mx, my), (Sink::ExoLake, Some(k)) if k == id)
                    })
            })
            .max_by(|&(a, _), &(b, _)| {
                if max_a_trunk {
                    // ties are INHERITED areas (a clipped fragment carries its parent's, F42/F93):
                    // broken by the shipped longest-path rule
                    by_area(a, b)
                } else {
                    pathlen[a].cmp(&pathlen[b])
                }
            })
            .map(|(j, _)| j)
    };
    let mut group_list: Vec<(usize, Vec<usize>)> = groups.into_iter().collect();
    group_list.sort_by_key(|(r, _)| *r);
    let mut out = Vec::new();
    for (r, members) in group_list {
        let mut trunk = Vec::new();
        let mut hopped = Vec::new();
        let mut cur = r;
        let mut lake_hops = 0usize;
        let mut hop_next = false;
        for _ in 0..=n {
            trunk.push(cur);
            hopped.push(std::mem::take(&mut hop_next));
            let up = segs[cur].upstream.iter().copied().filter(|&u| u < n);
            let next = if max_a_trunk {
                up.max_by(|&a, &b| by_area(a, b))
            } else {
                up.max_by_key(|&u| pathlen[u])
            };
            match next {
                Some(nx) => cur = nx,
                None => match inflow_to_lake_of(cur) {
                    Some(u) => {
                        lake_hops += 1;
                        hop_next = true;
                        cur = u
                    }
                    None => break,
                },
            }
        }
        trunk.reverse();
        hopped.reverse();
        let pts: usize = trunk.iter().map(|&s| segs[s].points.len()).sum();
        let (mx, my) = *segs[r].points.last().unwrap();
        let (sink, _) = classify_sink(dr, &wc, w, h, mx, my);
        out.push(Wc {
            trunk,
            segments: members,
            order: segs[r].strahler_order,
            catchment_km2: area(r),
            length_km: pts as f32 * km_per_cell * GEO_RATIO,
            sink,
            kind: dr.segment_kind.get(r).copied().unwrap_or(SegmentKind::Watercourse),
            lake_hops,
            hopped,
        });
    }
    out
}

#[test]
#[ignore]
fn f124_export() {
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    let cell_km2 = CELL_KM * CELL_KM;
    let r2 = GEO_RATIO * GEO_RATIO;
    eprintln!("\n==========  Finding 124-3 . does the export tell the truth?  ==========");
    let vc = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    // the constructed floor mask (Findings 121-123's approximation), for the rectangle's census
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
    let (_, masks) = carve(&pre, &sk, &vc, &ss);
    drop((pre, sk));
    let g = build_field_seed(
        Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        PSEED,
    );
    let (dr, bre) = inventory(&g, &ss);
    drop(g);
    let segs = &dr.rivers.segments;
    let agg = aggregate(&dr, &bre, true);
    let land_km2 = bre.data.iter().filter(|&&v| v > SEA).count() as f32 * cell_km2;

    // ── (i) partition: which cells belong to more than one object ──
    let mut seg2obj = vec![usize::MAX; segs.len()];
    for (o, x) in agg.iter().enumerate() {
        for &s in &x.segments {
            seg2obj[s] = o;
        }
    }
    let mut owners: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        let pts: &[(u32, u32)] = if s.downstream.is_some() && s.points.len() >= 2 {
            &s.points[..s.points.len() - 1] // the junction point sits on the parent
        } else {
            &s.points[..]
        };
        for &(x, y) in pts {
            let v = owners.entry(y as usize * w + x as usize).or_default();
            if !v.contains(&seg2obj[i]) {
                v.push(seg2obj[i]);
            }
        }
    }
    let shared: Vec<(&usize, &Vec<usize>)> = owners.iter().filter(|(_, v)| v.len() > 1).collect();
    let mut pair: HashMap<(usize, usize), usize> = HashMap::new();
    for (_, v) in &shared {
        for a in 0..v.len() {
            for b in a + 1..v.len() {
                let k = (v[a].min(v[b]), v[a].max(v[b]));
                *pair.entry(k).or_default() += 1;
            }
        }
    }
    let mut pv: Vec<_> = pair.into_iter().collect();
    pv.sort_by(|a, b| b.1.cmp(&a.1));
    eprintln!(
        "\n   (i) river cells **{}** · in more than one object **{} ({:.2} %)** · object pairs sharing \
         cells {}",
        owners.len(),
        shared.len(),
        100.0 * shared.len() as f32 / owners.len().max(1) as f32,
        pv.len()
    );
    for ((a, b), c) in pv.iter().take(6) {
        let d = |o: usize| {
            let x = &agg[o];
            let (mx, my) = *segs[*x.trunk.last().unwrap()].points.last().unwrap();
            format!(
                "{:?} S{} {:.0} km² signified, mouth ({mx},{my})",
                x.kind, x.order, x.catchment_km2
            )
        };
        eprintln!("       {c:>6} cells · [{}] × [{}]", d(*a), d(*b));
    }

    // ── (ii) Σ catchments against the land ──
    let sum = |kind: SegmentKind| -> f32 {
        agg.iter().filter(|x| x.kind == kind).map(|x| x.catchment_km2 / r2).sum()
    };
    let (sw, ss_) = (sum(SegmentKind::Watercourse), sum(SegmentKind::Spillway));
    let big = agg.iter().max_by(|a, b| a.catchment_km2.total_cmp(&b.catchment_km2)).unwrap();
    eprintln!(
        "   (ii) land **{land_km2:.0} km²** drawn (= {:.0} signified) · Σ catchments: watercourses **{sw:.0}** \
         ({:.2}× land) · spillways **{ss_:.0}** ({:.2}× land) · the biggest object {:?} **{:.0} km² signified \
         = {:.0} km² drawn**",
        land_km2 * r2,
        sw / land_km2,
        ss_ / land_km2,
        big.kind,
        big.catchment_km2,
        big.catchment_km2 / r2
    );

    // ── (iv) length ──
    let diag = DOMAIN_KM * std::f32::consts::SQRT_2;
    let mut by_len: Vec<&Wc> = agg.iter().filter(|x| x.kind == SegmentKind::Watercourse).collect();
    by_len.sort_by(|a, b| b.length_km.total_cmp(&a.length_km));
    eprintln!(
        "   (iv) domain diagonal **{diag:.0} km drawn = {:.0} km signified** · the longest trunks:",
        diag * GEO_RATIO
    );
    for x in by_len.iter().take(5) {
        let pts: usize = x.trunk.iter().map(|&s| segs[s].points.len()).sum();
        let uniq: std::collections::HashSet<(u32, u32)> =
            x.trunk.iter().flat_map(|&s| segs[s].points.iter().copied()).collect();
        eprintln!(
            "       {:.0} km signified = **{:.0} km drawn** · {} segments, {} points ({} distinct cells) · \
             lake hops {} · {:.0} km² signified",
            x.length_km,
            x.length_km / GEO_RATIO,
            x.trunk.len(),
            pts,
            uniq.len(),
            x.lake_hops,
            x.catchment_km2
        );
    }

    // ── (v) the rectangles: blocks of river points ──
    let river: Vec<bool> = {
        let mut r = vec![false; n];
        for (i, s) in segs.iter().enumerate() {
            if dr.segment_kind[i] == SegmentKind::Watercourse {
                for &(x, y) in &s.points {
                    r[y as usize * w + x as usize] = true;
                }
            }
        }
        r
    };
    let dense: Vec<bool> = (0..n)
        .map(|k| {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            river[k]
                && (-1i32..=1).all(|dy| {
                    (-1i32..=1).all(|dx| {
                        river[(y + dy).rem_euclid(h as i32) as usize * w
                            + (x + dx).rem_euclid(w as i32) as usize]
                    })
                })
        })
        .collect();
    let wc = ymir_core::lakes::connectivity::water_class(&bre, SEA);
    let mut seen = vec![false; n];
    let mut blocks = Vec::new();
    for s0 in 0..n {
        if !dense[s0] || seen[s0] {
            continue;
        }
        let mut comp = vec![s0];
        let mut q = vec![s0];
        seen[s0] = true;
        while let Some(k) = q.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nk = (y + dy).rem_euclid(h as i32) as usize * w
                        + (x + dx).rem_euclid(w as i32) as usize;
                    if dense[nk] && !seen[nk] {
                        seen[nk] = true;
                        q.push(nk);
                        comp.push(nk);
                    }
                }
            }
        }
        if comp.len() >= 50 {
            blocks.push(comp);
        }
    }
    blocks.sort_by(|a, b| b.len().cmp(&a.len()));
    eprintln!(
        "   (v) blocks of river points (every 3×3 neighbour a river point), ≥ 50 cells: **{}**",
        blocks.len()
    );
    for b in blocks.iter().take(6) {
        let (mut x0, mut x1, mut y0, mut y1) = (w, 0, h, 0);
        for &k in b {
            x0 = x0.min(k % w);
            x1 = x1.max(k % w);
            y0 = y0.min(k / w);
            y1 = y1.max(k / w);
        }
        let bbox = (x1 - x0 + 1) * (y1 - y0 + 1);
        let frac = |f: &dyn Fn(usize) -> bool| {
            100.0 * b.iter().filter(|&&k| f(k)).count() as f32 / b.len() as f32
        };
        // flat: the breached height equals its four neighbours'
        let flat = |k: usize| {
            let (x, y) = (k % w, k / w);
            [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)].iter().all(|&(dx, dy)| {
                let nk = (y as i32 + dy).rem_euclid(h as i32) as usize * w
                    + (x as i32 + dx).rem_euclid(w as i32) as usize;
                bre.data[nk] == bre.data[k]
            })
        };
        let bset: std::collections::HashSet<usize> = b.iter().copied().collect();
        let nseg: std::collections::HashSet<usize> = segs
            .iter()
            .enumerate()
            .filter(|(i, s)| {
                dr.segment_kind[*i] == SegmentKind::Watercourse
                    && s.points.iter().any(|&(x, y)| bset.contains(&(y as usize * w + x as usize)))
            })
            .map(|(i, _)| i)
            .collect();
        let accs = sorted(b.iter().map(|&k| dr.flow.accumulation.data[k] * cell_km2).collect());
        // ADR Finding 124-3 -- the width the export gives these reaches, and the width read at
        // their LAST point (the confluence: what the export gave before the own-end correction)
        let wexp = sorted(nseg.iter().map(|&i| dr.segment_width_m[i]).collect());
        let wlast = sorted(
            nseg.iter()
                .map(|&i| 5.0 * dr.segment_discharge_profile_m3s[i].last().copied().unwrap_or(0.0).max(0.0).sqrt())
                .collect(),
        );
        eprintln!(
            "       {} cells · bbox ({x0},{y0})–({x1},{y1}) {}×{} · fill {:.0} % · lake {:.0} % · inland below-sea \
             {:.0} % · flat {:.0} % · constructed floor {:.0} % · carved {:.0} % · {} segments · D8 A p50 {:.2} km² · \
             width exported p50 {:.1} m p90 {:.1} m · at the last point (pre-F124) p50 {:.1} m p90 {:.1} m",
            b.len(),
            x1 - x0 + 1,
            y1 - y0 + 1,
            100.0 * b.len() as f32 / bbox as f32,
            frac(&|k| dr.lake_map[k] != 0),
            frac(&|k| wc[k] == 2),
            frac(&flat),
            frac(&|k| masks.floor[k]),
            frac(&|k| masks.carved[k]),
            nseg.len(),
            pct(&accs, 0.5),
            pct(&wexp, 0.5),
            pct(&wexp, 0.9),
            pct(&wlast, 0.5),
            pct(&wlast, 0.9)
        );
    }
    eprintln!("\n==========  end 124-3 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 124-3 — **the 1 093 km trunk, dissected**: its segments from the mouth up, their area,
/// kind and length, how many of its cells lie in the 2 459-cell block of river points, and the
/// trunk's straight-line reach.
///
/// Run: cargo test -p ymir-core --release --test f124_export -- --ignored f124_trunk --nocapture
#[test]
#[ignore]
fn f124_trunk() {
    let ss = SteinSteinParams::default();
    let vc = ValleyConstruction::new(F121_AGE_K, Some(0.1));
    let g = build_field_seed(
        Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        PSEED,
    );
    let (w, _h) = (g.width, g.height);
    let (dr, bre) = inventory(&g, &ss);
    let segs = &dr.rivers.segments;
    let agg = aggregate(&dr, &bre, true);
    for target in [1093.0f32, 1173.0] {
        let x = agg
            .iter()
            .filter(|x| x.kind == SegmentKind::Watercourse)
            .min_by(|a, b| (a.length_km - target).abs().total_cmp(&(b.length_km - target).abs()))
            .unwrap();
        let root = *x.trunk.last().unwrap();
        let (mx, my) = *segs[root].points.last().unwrap();
        let src = *segs[x.trunk[0]].points.first().unwrap();
        let chord = (((src.0 as f32 - mx as f32).powi(2) + (src.1 as f32 - my as f32).powi(2)).sqrt())
            * CELL_KM;
        let in_block = x
            .trunk
            .iter()
            .flat_map(|&s| segs[s].points.iter())
            .filter(|&&(px, py)| (4818..=4876).contains(&px) && (2979..=3026).contains(&py))
            .count();
        eprintln!(
            "\n   trunk {:.0} km signified ({:.0} km drawn) · mouth ({mx},{my}) · source ({},{}) · straight \
             line {chord:.1} km drawn · {} segments · {in_block} trunk points in the 2 459-cell block",
            x.length_km,
            x.length_km / GEO_RATIO,
            src.0,
            src.1,
            x.trunk.len()
        );
        eprintln!("   from the mouth up: segment · points · area signified · first → last point");
        for &s in x.trunk.iter().rev().take(14) {
            let p = &segs[s].points;
            eprintln!(
                "      {s:>6} · {:>4} pts · {:>8.0} km² · {:?} → {:?}",
                p.len(),
                dr.segment_drainage_km2[s],
                p.first().unwrap(),
                p.last().unwrap()
            );
        }
        let _ = w;
    }
}

/// Water bodies: every lake id and every connected inland below-sea component, UNIFIED when they
/// touch (a lake filling a below-sea basin up to its spill level is one body). `body[k]` = the
/// canonical body of cell k, or `u32::MAX` on land / open sea.
fn water_bodies(dr: &C1DrainageResult, wc: &[u8], w: usize, h: usize) -> Vec<u32> {
    let n = w * h;
    let max_lake = dr.lake_map.iter().copied().max().unwrap_or(0) as usize;
    // components of inland below-sea water (8-connectivity, torus)
    let mut comp = vec![u32::MAX; n];
    let mut nc = 0u32;
    for s in 0..n {
        if wc[s] != 2 || comp[s] != u32::MAX {
            continue;
        }
        comp[s] = nc;
        let mut st = vec![s];
        while let Some(k) = st.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nk = (y + dy).rem_euclid(h as i32) as usize * w + (x + dx).rem_euclid(w as i32) as usize;
                    if wc[nk] == 2 && comp[nk] == u32::MAX {
                        comp[nk] = nc;
                        st.push(nk);
                    }
                }
            }
        }
        nc += 1;
    }
    // node ids: lake L → L, component c → max_lake + 1 + c; union where a cell is both
    let nn = max_lake + 1 + nc as usize;
    let mut p: Vec<usize> = (0..nn).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    for k in 0..n {
        if dr.lake_map[k] != 0 && comp[k] != u32::MAX {
            let (a, b) = (find(&mut p, dr.lake_map[k] as usize), find(&mut p, max_lake + 1 + comp[k] as usize));
            if a != b {
                p[a] = b;
            }
        }
    }
    (0..n)
        .map(|k| {
            if dr.lake_map[k] != 0 {
                find(&mut p, dr.lake_map[k] as usize) as u32
            } else if comp[k] != u32::MAX {
                find(&mut p, max_lake + 1 + comp[k] as usize) as u32
            } else {
                u32::MAX
            }
        })
        .collect()
}

/// The body a point of a reach touches (its 3×3), if any.
fn body_near(body: &[u32], w: usize, h: usize, x: u32, y: u32) -> Option<u32> {
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let k = (y as i32 + dy).rem_euclid(h as i32) as usize * w + (x as i32 + dx).rem_euclid(w as i32) as usize;
            if body[k] != u32::MAX {
                return Some(body[k]);
            }
        }
    }
    None
}

/// ADR Finding 124-3 — **the corrected aggregation** (the microscope's objects), three gestures:
/// (a) the trunk climbs by the segment's OWN geometric area — the D8 accumulation at its own mouth
///     (one point back when that point is a junction on the parent or a mouth shared with another
///     reach, Finding 47), never the area it inherits at the junction: every tooth carried its
///     collector's area and the 1 093 km trunk zig-zagged up 66 wall columns;
/// (b) roots whose reaches share a cell are ONE object (union): the common downstream is one trunk;
/// (c) a reach dying in a water body that has a SPILLWAY is chained to it, as Finding 45 chains an
///     exorheic lake to its outlet (the lake and the below-sea core it floods are one body).
/// The object's displayed catchment is the sum of its DISJOINT roots' own runoff-equivalent areas
/// (the exported unit, Finding 47), a spillway keeping its basin's.
fn aggregate_v2(dr: &C1DrainageResult, bre: &GridF32) -> Vec<Wc> {
    let (w, h) = (bre.width, bre.height);
    let wc = ymir_core::lakes::connectivity::water_class(bre, 0.5);
    let body = water_bodies(dr, &wc, w, h);
    let segs = &dr.rivers.segments;
    let n = segs.len();
    let q = |i: usize| dr.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
    let is_spill = |i: usize| dr.segment_kind.get(i) == Some(&SegmentKind::Spillway);
    let mut terminus_count: HashMap<(u32, u32), usize> = HashMap::new();
    for sg in segs.iter().filter(|sg| sg.downstream.is_none()) {
        *terminus_count.entry(*sg.points.last().unwrap()).or_default() += 1;
    }
    let own_idx = |i: usize| -> usize {
        let p = &segs[i].points;
        let back = segs[i].downstream.is_some() || terminus_count.get(p.last().unwrap()).copied().unwrap_or(0) > 1;
        if back && p.len() >= 2 { p.len() - 2 } else { p.len() - 1 }
    };
    // (a) own geometric area (cells) for the climb; own runoff-equivalent area (signified km²) shown
    let own_cells: Vec<f32> = (0..n)
        .map(|i| {
            let (x, y) = segs[i].points[own_idx(i)];
            dr.flow.accumulation.data[y as usize * w + x as usize]
        })
        .collect();
    let q300 = ymir_core::tectonics_c1::drainage::runoff_km2_to_m3s(300.0);
    let own_shown: Vec<f32> = (0..n)
        .map(|i| {
            if is_spill(i) {
                dr.segment_drainage_km2.get(i).copied().unwrap_or(0.0)
            } else {
                dr.segment_discharge_profile_m3s
                    .get(i)
                    .and_then(|p| p.get(own_idx(i)))
                    .map(|&v| v / q300)
                    .unwrap_or_else(|| dr.segment_drainage_km2.get(i).copied().unwrap_or(0.0))
            }
        })
        .collect();
    // exorheic lake → its outlet reach (Finding 45); body → its spillway (c)
    let mut outlet_of: HashMap<u32, usize> = HashMap::new();
    let mut spill_of: HashMap<u32, usize> = HashMap::new();
    for (i, sg) in segs.iter().enumerate() {
        let (sx, sy) = sg.points[0];
        if is_spill(i) {
            // a spillway's body: where it starts, else the lake it names (as the viz)
            let b = body_near(&body, w, h, sx, sy).or_else(|| {
                let id = dr.segment_source_lake.get(i).copied().flatten()?;
                dr.lake_map.iter().position(|&v| v == id).map(|k| body[k])
            });
            if let Some(b) = b {
                if spill_of.get(&b).is_none_or(|&j| q(i) > q(j)) {
                    spill_of.insert(b, i);
                }
            }
        } else if let (Sink::ExoLake, Some(id)) = classify_sink(dr, &wc, w, h, sx, sy) {
            if outlet_of.get(&id).is_none_or(|&j| q(i) > q(j)) {
                outlet_of.insert(id, i);
            }
        }
    }
    let end_body = |i: usize| -> Option<u32> {
        let &(x, y) = segs[i].points.last()?;
        body_near(&body, w, h, x, y)
    };
    let across = |i: usize| -> Option<usize> {
        if let Some(&s) = end_body(i).and_then(|b| spill_of.get(&b)) {
            return Some(s).filter(|&o| o != i);
        }
        let &(mx, my) = segs[i].points.last()?;
        match classify_sink(dr, &wc, w, h, mx, my) {
            (Sink::ExoLake, Some(id)) => outlet_of.get(&id).copied().filter(|&o| o != i),
            _ => None,
        }
    };
    let mut root = vec![0usize; n];
    for i in 0..n {
        let (mut j, mut seen) = (i, 0usize);
        loop {
            seen += 1;
            if seen > n {
                break;
            }
            match segs[j].downstream {
                Some(k) if k < n => j = k,
                _ => match across(j) {
                    Some(o) => j = o,
                    None => break,
                },
            }
        }
        root[i] = j;
    }
    // (b) union of the roots of reaches sharing a cell (a junction point on its parent excluded)
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    let mut owner: HashMap<usize, usize> = HashMap::new();
    for i in 0..n {
        for &(x, y) in own_points(&segs[i]) {
            let k = y as usize * w + x as usize;
            match owner.get(&k) {
                Some(&o) => {
                    let (a, b) = (find(&mut parent, root[o]), find(&mut parent, root[i]));
                    if a != b {
                        parent[a] = b;
                    }
                }
                None => {
                    owner.insert(k, i);
                }
            }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let g = find(&mut parent, root[i]);
        groups.entry(g).or_default().push(i);
    }
    let mut pathlen = vec![0u32; n];
    {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| segs[i].strahler_order);
        for _ in 0..16 {
            let mut changed = false;
            for &i in &order {
                let up = segs[i].upstream.iter().copied().filter(|&u| u < n).map(|u| pathlen[u]).max().unwrap_or(0);
                let v = up + segs[i].points.len() as u32;
                if v != pathlen[i] {
                    pathlen[i] = v;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    // the shipped 1 % tie, now on OWN areas (a spillway ranks by its basin's cells)
    let climb_area = |i: usize| -> f32 {
        if is_spill(i) { own_shown[i] / (GEO_RATIO * GEO_RATIO * CELL_KM * CELL_KM) } else { own_cells[i] }
    };
    let by_area = |a: usize, b: usize| -> std::cmp::Ordering {
        let (x, y) = (climb_area(a), climb_area(b));
        if (x - y).abs() <= 0.01 * x.max(y) { pathlen[a].cmp(&pathlen[b]) } else { x.total_cmp(&y) }
    };
    // upstream of a reach's head: the largest reach CHAINED to it -- the exact inverse of `across`,
    // so the trunk never climbs into a reach the chaining gave to another object (the first run's
    // inverse was rebuilt from the head's body and lake separately, and an exorheic outlet of a
    // body that also has a spillway climbed into the spillway's inflows)
    let mut chained_to: HashMap<usize, Vec<usize>> = HashMap::new();
    for j in 0..n {
        if segs[j].downstream.is_none() {
            if let Some(o) = across(j) {
                chained_to.entry(o).or_default().push(j);
            }
        }
    }
    let inflow_of = |i: usize| -> Option<usize> {
        chained_to.get(&i)?.iter().copied().max_by(|&a, &c| by_area(a, c))
    };
    let km_per_cell = DOMAIN_KM / w as f32;
    let mut glist: Vec<(usize, Vec<usize>)> = groups.into_iter().collect();
    glist.sort_by_key(|(r, _)| *r);
    let mut out = Vec::new();
    for (_, members) in glist {
        let mut roots: Vec<usize> = members.iter().map(|&i| root[i]).collect();
        roots.sort();
        roots.dedup();
        let rep = *roots.iter().max_by(|&&a, &&b| by_area(a, b)).unwrap();
        let mut trunk = Vec::new();
        let mut hopped = Vec::new();
        let (mut cur, mut hops) = (rep, 0usize);
        let mut hop_next = false;
        let mut visited = std::collections::HashSet::new();
        while visited.insert(cur) {
            trunk.push(cur);
            hopped.push(std::mem::take(&mut hop_next));
            let next = segs[cur].upstream.iter().copied().filter(|&u| u < n).max_by(|&a, &b| by_area(a, b));
            match next {
                Some(nx) => cur = nx,
                None => match inflow_of(cur) {
                    Some(u) => {
                        hops += 1;
                        hop_next = true;
                        cur = u
                    }
                    None => break,
                },
            }
        }
        trunk.reverse();
        hopped.reverse();
        let path_cells = trunk_length_cells(segs, &trunk, w, h); // as the viz
        let (mx, my) = *segs[rep].points.last().unwrap();
        let (sink, _) = classify_sink(dr, &wc, w, h, mx, my);
        // a RIVER as soon as one member is a watercourse (as the viz)
        let river = members.iter().any(|&s| !is_spill(s));
        out.push(Wc {
            trunk,
            segments: members,
            order: segs[rep].strahler_order,
            catchment_km2: roots.iter().map(|&r| own_shown[r]).sum(),
            length_km: path_cells * km_per_cell * GEO_RATIO,
            sink,
            kind: if river { SegmentKind::Watercourse } else { SegmentKind::Spillway },
            lake_hops: hops,
            hopped,
        });
    }
    out
}

/// ADR Finding 124-3 — the LENGTH of a trunk, in cells: the path through its reaches' points in
/// order, each step Euclidean on the torus (1 or √2 on D8, a crossed water body by its chord), the
/// confluence cell a reach shares with its receiver counted ONCE. It was `points × cell`, which
/// counted every confluence twice (+7.6 % on a 61 km trunk) and a diagonal step as one cell.
fn trunk_length_cells(segs: &[ymir_core::terrain::flow::RiverSegment], trunk: &[usize], w: usize, h: usize) -> f32 {
    let (mut len, mut prev): (f32, Option<(u32, u32)>) = (0.0, None);
    for &s in trunk {
        for &p in &segs[s].points {
            if let Some(q) = prev {
                if q != p {
                    let dx = (p.0 as i64 - q.0 as i64).rem_euclid(w as i64);
                    let dy = (p.1 as i64 - q.1 as i64).rem_euclid(h as i64);
                    let (dx, dy) = (dx.min(w as i64 - dx) as f32, dy.min(h as i64 - dy) as f32);
                    len += (dx * dx + dy * dy).sqrt();
                }
            }
            prev = Some(p);
        }
    }
    len
}

/// A reach's points without the junction point it shares with its parent.
fn own_points(s: &ymir_core::terrain::flow::RiverSegment) -> &[(u32, u32)] {
    if s.downstream.is_some() && s.points.len() >= 2 { &s.points[..s.points.len() - 1] } else { &s.points[..] }
}

/// The GEOMETRIC partition of the land by the objects (instrument, not export): every land cell
/// follows D8 to the first river cell or water body it meets and takes that object; a water body
/// belongs to the object of the reach leaving it (outlet / spillway), else of the largest reach
/// dying in it. Cells reaching the open sea without meeting a river are the uncaptured coastal
/// fringe. Returns (per-object km² drawn, uncaptured km² drawn, land km² drawn).
fn geometric_partition(dr: &C1DrainageResult, bre: &GridF32, agg: &[Wc]) -> (Vec<f32>, f32, f32) {
    let (w, h) = (bre.width, bre.height);
    let n = w * h;
    let cell_km2 = CELL_KM * CELL_KM;
    let wc = ymir_core::lakes::connectivity::water_class(bre, 0.5);
    let body = water_bodies(dr, &wc, w, h);
    let segs = &dr.rivers.segments;
    const UNK: u32 = u32::MAX;
    const SEAO: u32 = u32::MAX - 1;
    let mut obj = vec![UNK; n];
    let mut seg2obj = vec![UNK; segs.len()];
    for (o, x) in agg.iter().enumerate() {
        for &s in &x.segments {
            seg2obj[s] = o as u32;
        }
    }
    // bodies: the object leaving it, else the largest reach dying in it
    let mut body_obj: HashMap<u32, (bool, f32, u32)> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        let (sx, sy) = s.points[0];
        let a = dr.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
        if let Some(b) = body_near(&body, w, h, sx, sy) {
            let e = body_obj.entry(b).or_insert((true, a, seg2obj[i]));
            if !e.0 || a > e.1 {
                *e = (true, a, seg2obj[i]);
            }
        }
        if s.downstream.is_none() {
            let &(mx, my) = s.points.last().unwrap();
            if let Some(b) = body_near(&body, w, h, mx, my) {
                let e = body_obj.entry(b).or_insert((false, a, seg2obj[i]));
                if !e.0 && a > e.1 {
                    *e = (false, a, seg2obj[i]);
                }
            }
        }
    }
    for i in 0..segs.len() {
        for &(x, y) in own_points(&segs[i]) {
            let k = y as usize * w + x as usize;
            if obj[k] == UNK {
                obj[k] = seg2obj[i];
            }
        }
    }
    for k in 0..n {
        if obj[k] == UNK {
            if body[k] != u32::MAX {
                obj[k] = body_obj.get(&body[k]).map_or(SEAO, |e| e.2);
            } else if wc[k] == 1 {
                obj[k] = SEAO;
            }
        }
    }
    let dir = &dr.flow.direction;
    let mut path = Vec::new();
    for s in 0..n {
        if obj[s] != UNK {
            continue;
        }
        path.clear();
        let mut cur = s;
        let mut steps = 0usize;
        let res = loop {
            if obj[cur] != UNK {
                break obj[cur];
            }
            path.push(cur);
            steps += 1;
            let d = dir[cur];
            if d == ymir_core::terrain::flow::DIR_NONE || steps > n {
                break SEAO;
            }
            let (x, y) = ((cur % w) as i32, (cur / w) as i32);
            let nx = (x + ymir_core::terrain::flow::D8_DX[d as usize]).rem_euclid(w as i32) as usize;
            let ny = (y + ymir_core::terrain::flow::D8_DY[d as usize]).rem_euclid(h as i32) as usize;
            cur = ny * w + nx;
        };
        for &k in &path {
            obj[k] = res;
        }
    }
    // COUNTED IN CELLS: 67 M f32 increments of 0.0024 km² lose the sum (the first run read
    // Σ + fringe = 1.149 of the land)
    let mut per = vec![0u64; agg.len()];
    let (mut unc, mut land) = (0u64, 0u64);
    for k in 0..n {
        if wc[k] == 1 {
            continue; // the open sea is not land
        }
        land += 1;
        match obj[k] {
            SEAO | UNK => unc += 1,
            o => per[o as usize] += 1,
        }
    }
    (per.iter().map(|&c| c as f32 * cell_km2).collect(), unc as f32 * cell_km2, land as f32 * cell_km2)
}

/// ADR Finding 124-3 — the corrected aggregation against the shipped one: shared cells, Σ displayed
/// catchments against the land, the geometric partition, the longest trunks (drawn length ≤ the
/// domain diagonal, ASSERTED), and Finding 93's control (no trunk loses its lake).
///
/// Run: cargo test -p ymir-core --release --test f124_export -- --ignored f124_objects --nocapture
#[test]
#[ignore]
fn f124_objects() {
    let ss = SteinSteinParams::default();
    let r2 = GEO_RATIO * GEO_RATIO;
    let t0 = Instant::now();
    let diag = DOMAIN_KM * std::f32::consts::SQRT_2;
    eprintln!("\n==========  Finding 124-3 . the corrected objects  ==========");
    for (name, kn) in [
        (
            "C2/10",
            Knobs {
                valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))),
                slope_floor_abs: Some(S_EQ),
                ..Knobs::passes(2)
            },
        ),
        ("delivered", Knobs::passes(2)),
    ] {
        let g = build_field_seed(kn, PSEED);
        let (dr, bre) = inventory(&g, &ss);
        drop(g);
        let w = bre.width;
        let segs = &dr.rivers.segments;
        let wcl = ymir_core::lakes::connectivity::water_class(&bre, 0.5);
        let land_drawn = wcl.iter().filter(|&&c| c != 1).count() as f32 * CELL_KM * CELL_KM;
        let longest_path = aggregate(&dr, &bre, false);
        let shipped = aggregate(&dr, &bre, true);
        let corrected = aggregate_v2(&dr, &bre);
        for (rule, agg) in [("longest-path (pre-F123)", &longest_path), ("shipped (F123 max-A)", &shipped), ("corrected (F124)", &corrected)] {
            let mut seg2obj = vec![usize::MAX; segs.len()];
            for (o, x) in agg.iter().enumerate() {
                for &s in &x.segments {
                    seg2obj[s] = o;
                }
            }
            let mut owners: HashMap<usize, Vec<usize>> = HashMap::new();
            for (i, s) in segs.iter().enumerate() {
                for &(x, y) in own_points(s) {
                    let v = owners.entry(y as usize * w + x as usize).or_default();
                    if !v.contains(&seg2obj[i]) {
                        v.push(seg2obj[i]);
                    }
                }
            }
            let shared = owners.values().filter(|v| v.len() > 1).count();
            let shown: f32 = agg.iter().map(|x| x.catchment_km2 / r2).sum();
            let (per, unc, land) = geometric_partition(&dr, &bre, agg);
            let geo_sum: f32 = per.iter().sum();
            let mut longest: Vec<usize> = (0..agg.len()).collect();
            longest.sort_by(|&a, &b| agg[b].length_km.total_cmp(&agg[a].length_km));
            // Finding 93's control: the shipped objects' trunks, keyed by their mouth reach
            let lost = if std::ptr::eq(agg, &corrected) {
                longest_path
                    .iter()
                    .filter(|x| {
                        let o = seg2obj[*x.trunk.last().unwrap()];
                        o != usize::MAX && agg[o].lake_hops < x.lake_hops
                    })
                    .count()
                    .to_string()
            } else {
                "-".into()
            };
            eprintln!(
                "\n   {name:<9} {rule:<24} objects {} · cells in > 1 object **{shared}** · Σ displayed catchments / land \
                 **{:.2}** · geometric partition Σ / land **{:.3}** (uncaptured fringe {:.3}) · lake hops Σ {} · \
                 trunks losing a lake vs pre-F123 **{lost}**",
                agg.len(),
                shown / land_drawn,
                geo_sum / land,
                unc / land,
                agg.iter().map(|x| x.lake_hops).sum::<usize>()
            );
            for &o in longest.iter().take(4) {
                let x = &agg[o];
                let s = segs[x.trunk[0]].points[0];
                let m = *segs[*x.trunk.last().unwrap()].points.last().unwrap();
                let chord = ((s.0 as f32 - m.0 as f32).powi(2) + (s.1 as f32 - m.1 as f32).powi(2)).sqrt() * CELL_KM;
                eprintln!(
                    "       {:>5.0} km sig = **{:>4.0} km drawn** · straight line {chord:>4.0} km · sinuosity {:>5.1} · \
                     {} trunk segments · {:?} · displayed {:>7.0} km² sig · geometric {:>7.0} km² sig · mouth ({},{})",
                    x.length_km,
                    x.length_km / GEO_RATIO,
                    x.length_km / GEO_RATIO / chord.max(CELL_KM),
                    x.trunk.len(),
                    x.kind,
                    x.catchment_km2,
                    per[o] * r2,
                    m.0,
                    m.1
                );
            }
            if std::ptr::eq(agg, &corrected) {
                for x in agg.iter() {
                    assert!(
                        x.length_km / GEO_RATIO <= diag,
                        "a trunk longer than the domain diagonal: {:.0} km drawn > {diag:.0} km",
                        x.length_km / GEO_RATIO
                    );
                }
            }
        }
    }
    eprintln!("\n==========  end 124-3 objects . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// ADR Finding 124-3 — dissect the corrected objects the first f124_objects run flagged: the long
/// trunks on tiny catchments (C2/10), the old 1 093 / 1 173 km trunks' mouths, and the delivered
/// trunks that lose a lake crossing against the pre-F123 rule.
///
/// Run: cargo test -p ymir-core --release --test f124_export -- --ignored f124_dissect --nocapture
#[test]
#[ignore]
fn f124_dissect() {
    let ss = SteinSteinParams::default();
    let r2 = GEO_RATIO * GEO_RATIO;
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 124-3 . dissection  ==========");
    for (name, kn, mouths) in [
        (
            "C2/10",
            Knobs {
                valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))),
                slope_floor_abs: Some(S_EQ),
                ..Knobs::passes(2)
            },
            vec![(4936u32, 3051u32), (3758, 4737)],
        ),
        ("delivered", Knobs::passes(2), vec![]),
    ] {
        let g = build_field_seed(kn, PSEED);
        let (dr, bre) = inventory(&g, &ss);
        drop(g);
        let (w, h) = (bre.width, bre.height);
        let segs = &dr.rivers.segments;
        let wcl = ymir_core::lakes::connectivity::water_class(&bre, 0.5);
        let body = water_bodies(&dr, &wcl, w, h);
        let lp = aggregate(&dr, &bre, false);
        let cor = aggregate_v2(&dr, &bre);
        let (per, unc, land) = geometric_partition(&dr, &bre, &cor);
        let shown: f32 = cor.iter().map(|x| x.catchment_km2 / r2).sum();
        let geo: f32 = per.iter().sum();
        eprintln!(
            "\n   {name}: geometric partition Σ / land {:.3} + uncaptured {:.3} = {:.3} · Σ displayed / Σ geometric \
             **{:.2}** (runoff-equivalent / area)",
            geo / land,
            unc / land,
            (geo + unc) / land,
            shown / geo
        );
        let mut s2o = vec![usize::MAX; segs.len()];
        for (o, x) in cor.iter().enumerate() {
            for &s in &x.segments {
                s2o[s] = o;
            }
        }
        let walk = |x: &Wc, o: usize| {
            let pts: usize = x.trunk.iter().map(|&s| segs[s].points.len()).sum();
            let uniq: std::collections::HashSet<(u32, u32)> =
                x.trunk.iter().flat_map(|&s| segs[s].points.iter().copied()).collect();
            let foreign = x.trunk.iter().filter(|&&s| s2o[s] != o).count();
            eprintln!(
                "      object #{o}: {} trunk segments, {pts} points, {} distinct cells, {foreign} trunk segments \
                 OUTSIDE the object, {} members, geometric {:.1} km² drawn, displayed {:.1} km² drawn",
                x.trunk.len(),
                uniq.len(),
                x.segments.len(),
                per[o],
                x.catchment_km2 / r2
            );
            for (t, &si) in x.trunk.iter().enumerate().take(14) {
                let sg = &segs[si];
                let own = if sg.downstream.is_some() && sg.points.len() >= 2 {
                    sg.points[sg.points.len() - 2]
                } else {
                    *sg.points.last().unwrap()
                };
                let a = dr.flow.accumulation.data[own.1 as usize * w + own.0 as usize] * cell_km2;
                let (fx, fy) = sg.points[0];
                let (lx, ly) = *sg.points.last().unwrap();
                eprintln!(
                    "        {t:>2} seg {si:>6} {:?} {}pts ({fx},{fy})→({lx},{ly}) down {:?} own A {a:.3} km² · body \
                     start {:?} end {:?} · object {} · {}",
                    dr.segment_kind[si],
                    sg.points.len(),
                    sg.downstream,
                    body_near(&body, w, h, fx, fy),
                    body_near(&body, w, h, lx, ly),
                    if s2o[si] == o { "same".to_string() } else { format!("#{}", s2o[si]) },
                    if x.hopped.get(t + 1).copied().unwrap_or(false) { "then HOP" } else { "" }
                );
            }
        };
        let mut by_len: Vec<usize> = (0..cor.len()).collect();
        by_len.sort_by(|&a, &b| cor[b].length_km.total_cmp(&cor[a].length_km));
        eprintln!("   -- the four longest corrected trunks");
        for &o in by_len.iter().take(4) {
            walk(&cor[o], o);
        }
        for (mx, my) in mouths {
            if let Some(i) =
                segs.iter().position(|sg| sg.downstream.is_none() && *sg.points.last().unwrap() == (mx, my))
            {
                let o = s2o[i];
                eprintln!(
                    "   -- the old trunk at mouth ({mx},{my}): now {:.0} km drawn",
                    cor[o].length_km / GEO_RATIO
                );
                walk(&cor[o], o);
            }
        }
        eprintln!("   -- trunks that lose a lake crossing against the pre-F123 longest path");
        for x in &lp {
            let o = s2o[*x.trunk.last().unwrap()];
            if o != usize::MAX && cor[o].lake_hops < x.lake_hops {
                eprintln!(
                    "      pre-F123: {} hops, {:.0} km drawn, {} trunk segs, mouth {:?} · corrected object: {} hops, \
                     {:.0} km drawn, mouth {:?}",
                    x.lake_hops,
                    x.length_km / GEO_RATIO,
                    x.trunk.len(),
                    segs[*x.trunk.last().unwrap()].points.last().unwrap(),
                    cor[o].lake_hops,
                    cor[o].length_km / GEO_RATIO,
                    segs[*cor[o].trunk.last().unwrap()].points.last().unwrap(),
                );
                walk(&cor[o], o);
            }
        }
    }
    eprintln!("\n==========  end 124-3 dissection  ==========\n");
}

/// ADR Finding 124-3 — the TOPOLOGY census of the exported network (rivers.json's `upstream` /
/// `downstream` and the per-reach values), the two defects the dissection exposed:
/// (T1) links the lake clip leaves pointing at a run the water never reaches: `j ∈ i.upstream`
///      while `j.downstream ≠ i`, or `j.downstream = i` while j's water enters i nowhere on i
///      (neither i's first point nor the D8 receiver of j's last point lies on i);
/// (T2) reaches that END ON THEIR RECEIVER'S CONFLUENCE and so read its union value: count, and the
///      inflation of the geometric area, the runoff area and the width against the own point.
///
/// Run: cargo test -p ymir-core --release --test f124_export -- --ignored f124_topology --nocapture
#[test]
#[ignore]
fn f124_topology() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 124-3 . topology census  ==========");
    for (name, kn) in [
        (
            "C2/10",
            Knobs {
                valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))),
                slope_floor_abs: Some(S_EQ),
                ..Knobs::passes(2)
            },
        ),
        ("delivered", Knobs::passes(2)),
    ] {
        let g = build_field_seed(kn, PSEED);
        let (dr, bre) = inventory(&g, &ss);
        drop(g);
        let (w, h) = (bre.width, bre.height);
        let segs = &dr.rivers.segments;
        let n = segs.len();
        let dir = &dr.flow.direction;
        let recv = |x: u32, y: u32| -> Option<(u32, u32)> {
            let d = dir[y as usize * w + x as usize];
            (d != ymir_core::terrain::flow::DIR_NONE).then(|| {
                (
                    (x as i32 + ymir_core::terrain::flow::D8_DX[d as usize]).rem_euclid(w as i32) as u32,
                    (y as i32 + ymir_core::terrain::flow::D8_DY[d as usize]).rem_euclid(h as i32) as u32,
                )
            })
        };
        // (T1)
        let (mut links, mut asym, mut misplaced, mut to_terminal) = (0usize, 0usize, 0usize, 0usize);
        for i in 0..n {
            for &j in &segs[i].upstream {
                if j >= n {
                    continue;
                }
                links += 1;
                if segs[j].downstream != Some(i) {
                    asym += 1;
                    if segs[j].downstream.is_none() {
                        to_terminal += 1;
                    }
                }
            }
        }
        let mut down_links = 0usize;
        for j in 0..n {
            let Some(i) = segs[j].downstream else { continue };
            if i >= n {
                continue;
            }
            down_links += 1;
            let &(lx, ly) = segs[j].points.last().unwrap();
            let enters = segs[i].points.contains(&(lx, ly))
                || recv(lx, ly).is_some_and(|r| segs[i].points.contains(&r));
            if !enters {
                misplaced += 1;
            }
        }
        eprintln!(
            "\n   {name}: {n} reaches · upstream links {links}: **{asym} asymmetric** ({to_terminal} naming a reach \
             whose downstream is None) · downstream links {down_links}: **{misplaced} whose water enters the \
             named reach nowhere**"
        );
        // (T3) a reach ending ON another reach's first cell but linked elsewhere (a stolen junction)
        let mut starts: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (i, sg) in segs.iter().enumerate() {
            starts.entry(sg.points[0]).or_default().push(i);
        }
        let (mut stolen, mut orphan_heads) = (0usize, std::collections::HashSet::new());
        for (j, sg) in segs.iter().enumerate() {
            let Some(r) = sg.points.last().and_then(|p| starts.get(p)) else { continue };
            for &i in r {
                if i != j && sg.downstream.is_some() && sg.downstream != Some(i) && segs[i].upstream.is_empty() {
                    stolen += 1;
                    orphan_heads.insert(i);
                }
            }
        }
        eprintln!(
            "   (T3) reaches ending on the FIRST cell of a reach that has no upstream, linked elsewhere: **{stolen}** \
             (orphaned heads {})",
            orphan_heads.len()
        );
        // every reach holding a cell in its INTERIOR (a cell can sit inside a watercourse AND a
        // spillway traced over it: the first census kept one holder per cell)
        let mut inner: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (i, sg) in segs.iter().enumerate() {
            if sg.points.len() >= 3 {
                for p in &sg.points[1..sg.points.len() - 1] {
                    inner.entry(*p).or_default().push(i);
                }
            }
        }
        let (mut t3b, mut t3b_spill, mut t3b_self_kind) = (0usize, 0usize, [0usize; 2]);
        for (j, sg) in segs.iter().enumerate() {
            let Some(hs) = sg.points.last().and_then(|p| inner.get(p)) else { continue };
            let hs: Vec<usize> = hs.iter().copied().filter(|&i| i != j).collect();
            if hs.is_empty() || sg.downstream.is_some_and(|d| hs.contains(&d)) {
                continue;
            }
            let wc_holders: Vec<usize> =
                hs.iter().copied().filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse).collect();
            if wc_holders.is_empty() {
                t3b_spill += 1; // only a spillway runs through that cell
            } else {
                t3b += 1;
                t3b_self_kind[usize::from(dr.segment_kind[j] == SegmentKind::Spillway)] += 1;
            }
        }
        eprintln!(
            "   (T3b) reaches ending on a cell INSIDE another reach, linked elsewhere: **{t3b}** held by a watercourse \
             ({} of them watercourses, {} spillways) · {t3b_spill} held by a spillway only",
            t3b_self_kind[0], t3b_self_kind[1]
        );
        // (T2)
        let q300 = ymir_core::tectonics_c1::drainage::runoff_km2_to_m3s(300.0);
        let (mut on_conf, mut small) = (0usize, 0usize);
        let (mut rg, mut rr, mut rw) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            let s = &segs[i];
            let Some(d) = s.downstream else { continue };
            // on the receiver's first cell OR one it runs through (T3b)
            if d >= n || s.points.len() < 2 || !segs[d].points.contains(s.points.last().unwrap()) {
                continue;
            }
            on_conf += 1;
            let own = s.points[s.points.len() - 2];
            let a_own = dr.flow.accumulation.data[own.1 as usize * w + own.0 as usize];
            let a_last = dr.flow.accumulation.data[s.points.last().unwrap().1 as usize * w + s.points.last().unwrap().0 as usize];
            rg.push(a_last / a_own.max(1.0));
            let prof = &dr.segment_discharge_profile_m3s[i];
            let q_own = prof[prof.len() - 2];
            let km2_own = q_own / q300;
            rr.push(dr.segment_drainage_km2[i] / km2_own.max(1e-6));
            let w_own = 5.0 * q_own.max(0.0).powf(0.5);
            rw.push(dr.segment_width_m[i] / w_own.max(1e-6));
            if a_own * cell_km2 < 1.0 {
                small += 1;
            }
        }
        let (rg, rr, rw) = (sorted(rg), sorted(rr), sorted(rw));
        let over = |v: &[f32], t: f32| 100.0 * v.iter().filter(|&&x| x > t).count() as f32 / v.len().max(1) as f32;
        eprintln!(
            "   reaches with a downstream: {down_links} · ENDING ON THE RECEIVER'S CONFLUENCE **{on_conf}** ({small} of \
             them with an own area < 1 km² drawn)\n   union / own — geometric area p50 {:.2} p90 {:.1} (> 2×: {:.0} %) · \
             runoff area p50 {:.2} p90 {:.1} (> 2×: {:.0} %) · width p50 {:.2} p90 {:.1} (> 2×: {:.0} %)",
            pct(&rg, 0.5),
            pct(&rg, 0.9),
            over(&rg, 2.0),
            pct(&rr, 0.5),
            pct(&rr, 0.9),
            over(&rr, 2.0),
            pct(&rw, 0.5),
            pct(&rw, 0.9),
            over(&rw, 2.0)
        );
    }
    eprintln!("\n==========  end 124-3 topology  ==========\n");
}

/// ADR Finding 124-3 — the corrected trunk HEADS that carry area (Finding 93's residual on the
/// delivered world): a head reach with no upstream and no water body at its start, yet with tens to
/// thousands of km² at its own end. Where does that area come from? For each such head: its start's
/// 5×5 (lake ids / water class), the D8 DONORS of its first point (cells whose receiver it is), their
/// accumulation, lake id, class and the reach holding them.
///
/// Run: cargo test -p ymir-core --release --test f124_export -- --ignored f124_heads --nocapture
#[test]
#[ignore]
fn f124_heads() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 124-3 . heads that carry area  ==========");
    let g = build_field_seed(Knobs::passes(2), PSEED);
    let (dr, bre) = inventory(&g, &ss);
    drop(g);
    let (w, h) = (bre.width, bre.height);
    let segs = &dr.rivers.segments;
    let wcl = ymir_core::lakes::connectivity::water_class(&bre, 0.5);
    let body = water_bodies(&dr, &wcl, w, h);
    let dir = &dr.flow.direction;
    let mut cell_seg: HashMap<usize, usize> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        for &(x, y) in &s.points {
            cell_seg.entry(y as usize * w + x as usize).or_insert(i);
        }
    }
    let cor = aggregate_v2(&dr, &bre);
    let lp = aggregate(&dr, &bre, false);
    let mut s2o = vec![usize::MAX; segs.len()];
    for (o, x) in cor.iter().enumerate() {
        for &s in &x.segments {
            s2o[s] = o;
        }
    }
    // the population: every corrected trunk head with no upstream and no body at its start
    let mut heads: Vec<(usize, f32)> = cor
        .iter()
        .map(|x| x.trunk[0])
        .filter(|&i| {
            let (sx, sy) = segs[i].points[0];
            segs[i].upstream.is_empty()
                && body_near(&body, w, h, sx, sy).is_none()
                && dr.segment_kind[i] == SegmentKind::Watercourse
        })
        .map(|i| {
            let (sx, sy) = segs[i].points[0];
            (i, dr.flow.accumulation.data[sy as usize * w + sx as usize] * cell_km2)
        })
        .collect();
    heads.sort_by(|a, b| b.1.total_cmp(&a.1));
    let big = heads.iter().filter(|h| h.1 > 10.0).count();
    eprintln!(
        "\n   corrected trunk heads with no upstream and no body at their start: {} · with > 10 km² drawn AT \
         THEIR FIRST POINT: **{big}**",
        heads.len()
    );
    let losers: std::collections::HashSet<usize> = lp
        .iter()
        .filter_map(|x| {
            let o = s2o[*x.trunk.last().unwrap()];
            (o != usize::MAX && cor[o].lake_hops < x.lake_hops).then(|| cor[o].trunk[0])
        })
        .collect();
    for &(i, a0) in heads.iter().filter(|h| h.1 > 10.0 || losers.contains(&h.0)).take(8) {
        let (sx, sy) = segs[i].points[0];
        eprintln!(
            "\n   head seg {i}{} · {} pts · start ({sx},{sy}) A {a0:.1} km² drawn · strahler {} · segment_kind {:?}",
            if losers.contains(&i) { " (a Finding 93 LOSER)" } else { "" },
            segs[i].points.len(),
            segs[i].strahler_order,
            dr.segment_kind[i]
        );
        let mut grid = String::new();
        for dy in -2i32..=2 {
            grid.push_str("        ");
            for dx in -2i32..=2 {
                let k = (sy as i32 + dy).rem_euclid(h as i32) as usize * w + (sx as i32 + dx).rem_euclid(w as i32) as usize;
                let id = dr.lake_map[k];
                let tag = if id != 0 {
                    format!("L{:<8}", id)
                } else {
                    format!("c{}{:<7}", wcl[k], if cell_seg.contains_key(&k) { "r" } else { "" })
                };
                grid.push_str(&tag);
            }
            grid.push('\n');
        }
        eprint!("{grid}");
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let (nx, ny) = ((sx as i32 + dx).rem_euclid(w as i32) as usize, (sy as i32 + dy).rem_euclid(h as i32) as usize);
                let k = ny * w + nx;
                let d = dir[k];
                if d == ymir_core::terrain::flow::DIR_NONE {
                    continue;
                }
                let rx = (nx as i32 + ymir_core::terrain::flow::D8_DX[d as usize]).rem_euclid(w as i32) as u32;
                let ry = (ny as i32 + ymir_core::terrain::flow::D8_DY[d as usize]).rem_euclid(h as i32) as u32;
                if (rx, ry) != (sx, sy) {
                    continue;
                }
                let sg = cell_seg.get(&k).copied();
                eprintln!(
                    "        donor ({nx},{ny}) A {:.1} km² · lake {} · class {} · reach {:?}{}",
                    dr.flow.accumulation.data[k] * cell_km2,
                    dr.lake_map[k],
                    wcl[k],
                    sg,
                    sg.map(|j| format!(
                        " ({:?}, {} pts, down {:?}, ends ({},{}))",
                        dr.segment_kind[j],
                        segs[j].points.len(),
                        segs[j].downstream,
                        segs[j].points.last().unwrap().0,
                        segs[j].points.last().unwrap().1
                    ))
                    .unwrap_or_default()
                );
            }
        }
    }
    eprintln!("\n==========  end 124-3 heads  ==========\n");
}

/// ADR Finding 124-3 (d) — **the rectangle, reach by reach.** The biggest dense block of river points
/// on C2/10 (bbox (4818,2979)–(4876,3026)) carries 66 reaches whose exported width reads 29.7 m p50
/// for a cell accumulation of 0.16 km² p50 — a width no 0.16 km² catchment can drive. Every reach
/// touching the bbox: points, first → last, own geometric area, own discharge, width, links; and the
/// D8 directions of the block's river cells.
///
/// Run: cargo test -p ymir-core --release --test f124_export -- --ignored f124_block --nocapture
#[test]
#[ignore]
fn f124_block() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let r2 = GEO_RATIO * GEO_RATIO;
    eprintln!("\n==========  Finding 124-3 (d) . the rectangle, reach by reach  ==========");
    let g = build_field_seed(
        Knobs {
            valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))),
            slope_floor_abs: Some(S_EQ),
            ..Knobs::passes(2)
        },
        PSEED,
    );
    let (dr, bre) = inventory(&g, &ss);
    drop(g);
    let w = bre.width;
    let segs = &dr.rivers.segments;
    let inside = |&(x, y): &(u32, u32)| (4818..=4876).contains(&x) && (2979..=3026).contains(&y);
    let mut rows: Vec<usize> = (0..segs.len())
        .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse && segs[i].points.iter().any(inside))
        .collect();
    let own_a = |i: usize| {
        let s = &segs[i];
        let p = if s.downstream.is_some() && s.points.len() >= 2 { s.points[s.points.len() - 2] } else { *s.points.last().unwrap() };
        dr.flow.accumulation.data[p.1 as usize * w + p.0 as usize] * cell_km2
    };
    rows.sort_by(|&a, &b| own_a(b).total_cmp(&own_a(a)));
    let q300 = ymir_core::tectonics_c1::drainage::runoff_km2_to_m3s(300.0);
    eprintln!("\n   {} reaches touch the bbox (own area in km² drawn; Q and width signified)", rows.len());
    for &i in rows.iter().take(40) {
        let s = &segs[i];
        let qa = dr.segment_discharge_m3s[i];
        // the runoff depth the discharge implies over the own geometric area (mm/yr, real units)
        let depth = qa / r2 / q300 * 300.0 / own_a(i).max(1e-6);
        let in_pts = s.points.iter().filter(|p| inside(p)).count();
        eprintln!(
            "      seg {i:>6} S{} {:>4} pts ({} in) ({},{})→({},{}) · own A {:>8.3} km² · Q {:>7.2} m³/s · width {:>5.1} m \
             · implied runoff {:>9.0} mm/yr · up {} · down {:?}",
            s.strahler_order,
            s.points.len(),
            in_pts,
            s.points[0].0,
            s.points[0].1,
            s.points.last().unwrap().0,
            s.points.last().unwrap().1,
            own_a(i),
            qa,
            dr.segment_width_m[i],
            depth,
            s.upstream.len(),
            s.downstream
        );
    }
    // the D8 directions of the bbox's river cells
    let mut hist = [0usize; 9];
    let mut riv = std::collections::HashSet::new();
    for &i in &rows {
        for p in segs[i].points.iter().filter(|p| inside(p)) {
            riv.insert(*p);
        }
    }
    for &(x, y) in &riv {
        let d = dr.flow.direction[y as usize * w + x as usize];
        hist[if d == ymir_core::terrain::flow::DIR_NONE { 8 } else { d as usize }] += 1;
    }
    eprintln!(
        "\n   D8 directions of the {} river cells in the bbox (0 N, 1 NE, 2 E, 3 SE, 4 S, 5 SW, 6 W, 7 NW, none): {:?}",
        riv.len(),
        hist
    );
    // the bbox's accumulation on its river cells, and the runoff depth of the land (the reference)
    let mut accs: Vec<f32> = riv.iter().map(|&(x, y)| dr.flow.accumulation.data[y as usize * w + x as usize] * cell_km2).collect();
    accs.sort_by(|a, b| a.total_cmp(b));
    eprintln!(
        "   accumulation on those cells (km² drawn): p10 {:.2} · p50 {:.2} · p90 {:.2} · max {:.2}",
        accs[accs.len() / 10],
        accs[accs.len() / 2],
        accs[accs.len() * 9 / 10],
        accs[accs.len() - 1]
    );
    eprintln!("\n==========  end 124-3 (d) block  ==========\n");
}
