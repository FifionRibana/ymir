//! ADR 0001 Finding 123 — **the teeth of the comb: are they flow lines on the constructed WALLS?**
//!
//! The author named, on C2/10: teeth perpendicular to the trunk and attached to it; the same comb
//! as spillways into a lake; round lake arms (Finding 122). Hypothesis: the teeth are D8 flow lines
//! on the constructed walls, planar by construction — Finding 111's comb on an ideal surface.
//!
//! ## Instrument (rule 12)
//!
//! * **Walls** — `carved && !floor` of the PRE carve at `F121_AGE_K` (Findings 121/122's mask
//!   approximation, declared).
//! * **Drawn segments** — the export's `Watercourse` segments (`full_tree = false`: the layer the
//!   viz draws). A TOOTH is a segment with ≥ 80 % of its cells on a wall; "on walls" is ≥ 50 %.
//! * **A** — the D8 drained area in-domain at the segment's cells (the geo ratio multiplies the
//!   exported arrays only, Finding 68).
//! * **Angle to the parent** — the acute angle between the tooth's last 8 cells and the parent's
//!   cells around the junction.
//! * **The microscope's watercourses** — `aggregate_watercourses` of the viz (`workspace.rs`),
//!   PORTED here line for line (root by downstream links crossing exorheic lakes; trunk = the
//!   longest upstream path in points at each confluence, continued up the longest inflow of an
//!   exorheic lake; `length_km = points × km/cell × ratio`).
//!
//! Run: cargo test -p ymir-core --release --test f123_teeth -- --ignored --nocapture

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
const TOOTH: f32 = 0.80;
const ON_WALL: f32 = 0.50;

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
        let mut cur = r;
        let mut lake_hops = 0usize;
        for _ in 0..=n {
            trunk.push(cur);
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
                        cur = u
                    }
                    None => break,
                },
            }
        }
        trunk.reverse();
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
        });
    }
    out
}

fn angle_deg(a: (f32, f32), b: (f32, f32)) -> f32 {
    let (na, nb) = ((a.0 * a.0 + a.1 * a.1).sqrt(), (b.0 * b.0 + b.1 * b.1).sqrt());
    if na == 0.0 || nb == 0.0 {
        return f32::NAN;
    }
    let c = ((a.0 * b.0 + a.1 * b.1) / (na * nb)).abs().min(1.0);
    c.acos().to_degrees()
}

#[test]
#[ignore]
fn f123_teeth() {
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 123-A . the teeth of the comb  ==========");
    let vc = ValleyConstruction::f121(F121_AGE_K, Some(0.1));
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
    let (_, masks) = carve(&pre, &sk, &vc, &ss);
    let wall: Vec<bool> = (0..n).map(|k| masks.carved[k] && !masks.floor[k]).collect();
    drop((pre, masks, sk));

    let c2 = build_field_seed(
        Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        PSEED,
    );
    let (dr, bre) = inventory(&c2, &ss);
    drop(c2);
    let segs = &dr.rivers.segments;
    let acc = |x: u32, y: u32| dr.flow.accumulation.data[y as usize * w + x as usize] * cell_km2;

    // ── per drawn segment ──
    struct S {
        wall: f32,
        amax: f32,
        len_km: f32,
        order: u8,
        angle: f32,
        trunk: bool,
    }
    let mut ss_: Vec<S> = Vec::new();
    for (i, s) in segs.iter().enumerate() {
        if dr.segment_kind[i] != SegmentKind::Watercourse || s.points.len() < 2 {
            ss_.push(S {
                wall: f32::NAN,
                amax: 0.0,
                len_km: 0.0,
                order: 0,
                angle: f32::NAN,
                trunk: false,
            });
            continue;
        }
        let on = s.points.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count();
        // the segment's OWN mouth: the last point is the junction ON the parent, whose area is the
        // parent's (the first run read the trunk's area there: teeth p50 29 km²)
        let own = if s.downstream.is_some() && s.points.len() >= 2 {
            &s.points[..s.points.len() - 1]
        } else {
            &s.points[..]
        };
        let amax = own.iter().map(|&(x, y)| acc(x, y)).fold(0f32, f32::max);
        let mut len = 0f64;
        for p in s.points.windows(2) {
            let diag = p[0].0 != p[1].0 && p[0].1 != p[1].1;
            len += if diag { CELL_KM as f64 * std::f64::consts::SQRT_2 } else { CELL_KM as f64 };
        }
        // angle to the parent at the junction
        let angle = s.downstream.map_or(f32::NAN, |d| {
            let pts = &s.points;
            let e = pts.len() - 1;
            let b = e.saturating_sub(8);
            let dir_t = (pts[e].0 as f32 - pts[b].0 as f32, pts[e].1 as f32 - pts[b].1 as f32);
            let par = &segs[d].points;
            let (jx, jy) = pts[e];
            let j = par
                .iter()
                .enumerate()
                .min_by_key(|(_, p)| {
                    (p.0 as i64 - jx as i64).abs() + (p.1 as i64 - jy as i64).abs()
                })
                .map_or(0, |(j, _)| j);
            let (a, c) = (j.saturating_sub(8), (j + 8).min(par.len() - 1));
            let dir_p = (par[c].0 as f32 - par[a].0 as f32, par[c].1 as f32 - par[a].1 as f32);
            angle_deg(dir_t, dir_p)
        });
        ss_.push(S {
            wall: on as f32 / s.points.len() as f32,
            amax,
            len_km: len as f32,
            order: s.strahler_order,
            angle,
            trunk: amax >= 10.0,
        });
    }
    let wcs: Vec<usize> =
        (0..segs.len()).filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse).collect();
    let teeth: Vec<usize> = wcs.iter().copied().filter(|&i| ss_[i].wall >= TOOTH).collect();
    let onwall: Vec<usize> = wcs.iter().copied().filter(|&i| ss_[i].wall >= ON_WALL).collect();
    let tot_len: f32 = wcs.iter().map(|&i| ss_[i].len_km).sum();
    let teeth_len: f32 = teeth.iter().map(|&i| ss_[i].len_km).sum();
    let ta = sorted(teeth.iter().map(|&i| ss_[i].amax).collect());
    let tang = sorted(teeth.iter().map(|&i| ss_[i].angle).filter(|a| a.is_finite()).collect());
    let tl = sorted(teeth.iter().map(|&i| ss_[i].len_km).collect());
    let mut tord: HashMap<u8, usize> = HashMap::new();
    for &i in &teeth {
        *tord.entry(ss_[i].order).or_default() += 1;
    }
    let mut to: Vec<_> = tord.into_iter().collect();
    to.sort();
    let trunks: Vec<usize> = wcs.iter().copied().filter(|&i| ss_[i].trunk).collect();
    let trw = sorted(trunks.iter().map(|&i| ss_[i].wall).collect());
    let p = |v: &[f32], q: f64| if v.is_empty() { f32::NAN } else { pct(v, q) };
    eprintln!(
        "\n   drawn Watercourse segments **{}** ({:.0} km in-domain) · on walls (≥ 50 %) **{}** · TEETH \
         (≥ 80 %) **{}** = **{:.1} %** of the on-wall ones · teeth length **{:.0} km = {:.1} %** of the \
         network\n   teeth: A max p10/p50/p90 **{:.3} / {:.3} / {:.3} km²** (A_c = 0.1) · length p50 **{:.2} \
         km** · angle to the parent p10/p50/p90 **{:.0} / {:.0} / {:.0}°** · Strahler {:?}\n   \
         CONTROL · trunk segments (A max ≥ 10 km², {}): wall share p50 **{:.1} %**, p90 {:.1} %",
        wcs.len(),
        tot_len,
        onwall.len(),
        teeth.len(),
        100.0 * teeth.len() as f32 / onwall.len().max(1) as f32,
        teeth_len,
        100.0 * teeth_len / tot_len.max(1e-6),
        p(&ta, 0.10),
        p(&ta, 0.50),
        p(&ta, 0.90),
        p(&tl, 0.50),
        p(&tang, 0.10),
        p(&tang, 0.50),
        p(&tang, 0.90),
        to,
        trunks.len(),
        100.0 * p(&trw, 0.50),
        100.0 * p(&trw, 0.90)
    );
    let non: Vec<usize> = wcs.iter().copied().filter(|&i| ss_[i].wall < ON_WALL).collect();
    let nang = sorted(non.iter().map(|&i| ss_[i].angle).filter(|a| a.is_finite()).collect());
    eprintln!(
        "   off-wall segments ({}): angle to the parent p50 **{:.0}°** · A max p50 **{:.3} km²**",
        non.len(),
        p(&nang, 0.50),
        p(&sorted(non.iter().map(|&i| ss_[i].amax).collect()), 0.50)
    );

    // ── the microscope's watercourses, and the "S2, 580 km², 654 km" entry ──
    let agg = aggregate(&dr, &bre, false);
    let agg_a = aggregate(&dr, &bre, true);
    let mut cand: Vec<(usize, f32)> = agg
        .iter()
        .enumerate()
        .filter(|(_, x)| x.kind == SegmentKind::Watercourse && x.order == 2)
        .map(|(i, x)| {
            let d =
                ((x.catchment_km2 - 580.0) / 580.0).abs() + ((x.length_km - 654.0) / 654.0).abs();
            (i, d)
        })
        .collect();
    cand.sort_by(|a, b| a.1.total_cmp(&b.1));
    eprintln!("\n   the microscope's S2 entries closest to (580 km², 654 km), signified:");
    for &(i, d) in cand.iter().take(5) {
        let x = &agg[i];
        let rep = x.trunk.len() - x.trunk.iter().collect::<std::collections::HashSet<_>>().len();
        let tooth_in_trunk = x.trunk.iter().filter(|&&s| ss_[s].wall >= TOOTH).count();
        let teeth_members = x.segments.iter().filter(|&&s| ss_[s].wall >= TOOTH).count();
        let trunk_no_teeth: f32 = x
            .trunk
            .iter()
            .filter(|&&s| ss_[s].wall < TOOTH)
            .map(|&s| segs[s].points.len() as f32 * CELL_KM * GEO_RATIO)
            .sum();
        // the same group under the max-A climb (C2 identity)
        let root = *x.trunk.last().unwrap();
        let alt = agg_a.iter().find(|y| y.trunk.last() == Some(&root));
        eprintln!(
            "   · {:.0} km² · trunk **{:.0} km** in {} segments ({} repeated, {} teeth) · members {} ({} \
             teeth) · sink {:?} · score {:.3}\n     trunk WITHOUT its teeth **{:.0} km** · max-A trunk \
             (C2 identity) **{}**",
            x.catchment_km2,
            x.length_km,
            x.trunk.len(),
            rep,
            tooth_in_trunk,
            x.segments.len(),
            teeth_members,
            x.sink,
            d,
            trunk_no_teeth,
            alt.map_or("—".into(), |y| format!(
                "{:.0} km in {} segments",
                y.length_km,
                y.trunk.len()
            ))
        );
    }
    // whole-list view of the identity change (F93's control: trunks through lakes must not regress)
    let (mut longer, mut shorter, mut same) = (0usize, 0usize, 0usize);
    let mut ratio = Vec::new();
    for x in agg.iter().filter(|x| x.kind == SegmentKind::Watercourse) {
        let root = *x.trunk.last().unwrap();
        if let Some(y) = agg_a.iter().find(|y| y.trunk.last() == Some(&root)) {
            match y.length_km.partial_cmp(&x.length_km) {
                Some(std::cmp::Ordering::Greater) => longer += 1,
                Some(std::cmp::Ordering::Less) => shorter += 1,
                _ => same += 1,
            }
            if x.length_km > 0.0 {
                ratio.push(y.length_km / x.length_km);
            }
        }
    }
    let ratio = sorted(ratio);
    eprintln!(
        "\n   identity (C2 read, NOT applied): max-A trunk vs shipped longest-path trunk over {} watercourses: \
         same **{same}** · shorter **{shorter}** · longer **{longer}** · length ratio p10/p50/p90 {:.2} / \
         {:.2} / {:.2}",
        same + shorter + longer,
        p(&ratio, 0.10),
        p(&ratio, 0.50),
        p(&ratio, 0.90)
    );

    // ── C1 read: the drawn layer at a threshold ──
    eprintln!("\n   ── C1 (read, not applied) · the drawn layer at a threshold ──");
    let s2_sea: Vec<&Wc> = agg
        .iter()
        .filter(|x| x.kind == SegmentKind::Watercourse && x.order == 2 && x.sink == Sink::Sea)
        .collect();
    for (label, keep) in [
        ("A ≥ 0.1 km²", Box::new(|i: usize| ss_[i].amax >= 0.1) as Box<dyn Fn(usize) -> bool>),
        ("A ≥ 0.5 km²", Box::new(|i: usize| ss_[i].amax >= 0.5)),
        ("A ≥ 1 km²", Box::new(|i: usize| ss_[i].amax >= 1.0)),
        ("Strahler ≥ 2", Box::new(|i: usize| ss_[i].order >= 2)),
    ] {
        let kept = wcs.iter().filter(|&&i| keep(i)).count();
        let kept_teeth = teeth.iter().filter(|&&i| keep(i)).count();
        let lost_s2 = s2_sea.iter().filter(|x| !x.segments.iter().any(|&s| keep(s))).count();
        eprintln!(
            "   {label:<14} segments kept **{kept} of {}** · teeth kept **{kept_teeth} of {}** · S2 \
             watercourses to the sea that vanish entirely **{lost_s2} of {}**",
            wcs.len(),
            teeth.len(),
            s2_sea.len()
        );
    }
    let _ = h;
    eprintln!("\n==========  end 123-A . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// The viz's auto-framing origin for seed 1 (the author's export manifest).
const VIZ_ORIGIN: [f64; 2] = [0.09375, 0.578_125];
/// Finding 109's recorded hash of seed 1's delivered field at the benches' origin (ADR L17300).
const F109_DELIVERED: u64 = 0x316c_a8f7_27ca_6f62;

/// ADR Finding 123 — **the identity guard's reference, built AT THE VIZ'S SETTINGS.**
///
/// 0. Control: the delivered field at the benches' origin must still hash to Finding 109's
///    `316ca8f727ca6f62` — every seam added since (a1_exempt, valley_construction, basin_base, the
///    config refactor) claims to be inert by default, and this is where that claim is checked.
/// 1. Is a change of framing a pure translation? The delivered field at the viz's origin against
///    the benches' field rolled by 768 cells.
/// 2. The five states at the viz's origin: digest + field hash → `data/bench_field_hashes.json`.
/// 3. The "S2, 580 km², 654 km" entry, looked for again in the C2 /10 world AT THE VIZ'S ORIGIN.
///
/// Run: cargo test -p ymir-core --release --test f123_teeth -- --ignored f123_guard --nocapture
#[test]
#[ignore]
fn f123_guard() {
    use common::{bench_eroded_digest, build_field_seed};
    use ymir_core::tectonics_c1::bench_guard::{GuardEntry, field_hash};
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 123 . the viz / bench identity guard  ==========");

    // 0 · the inertness control
    let del_b = build_field_seed(
        Knobs { origin: Some(common::FINDINGS_107_123_ORIGIN), ..Knobs::passes(2) },
        PSEED,
    );
    let hb = field_hash(&del_b);
    eprintln!(
        "   0 · delivered at the benches' origin: {hb:016x} vs Finding 109's {F109_DELIVERED:016x} ⇒ **{}**",
        if hb == F109_DELIVERED {
            "BIT-IDENTICAL — every seam since F109 is inert"
        } else {
            "**DIFFERENT**"
        }
    );
    assert_eq!(hb, F109_DELIVERED, "a seam added since Finding 109 is NOT inert by default");

    let k = F121_AGE_K;
    let with_closure = |vc: ValleyConstruction| Knobs {
        valley: Some(vc),
        slope_floor_abs: Some(S_EQ),
        ..Knobs::passes(2)
    };
    // ADR Finding 124 -- the viz's states: the construction's DEFINITION (`new`, basin base) and
    // Finding 121's historical sea base beside it
    let states: [(&str, Knobs); 6] = [
        ("livré", Knobs::passes(2)),
        ("A1+B2", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
        ("C1 nue", Knobs { valley: Some(ValleyConstruction::new(k, None)), ..Knobs::passes(2) }),
        ("C2 /10 col (défaut)", with_closure(ValleyConstruction::new(k, Some(0.1)))),
        ("C2 /3", with_closure(ValleyConstruction::new(k, Some(1.0 / 3.0)))),
        ("C2 /10 niveau mer", with_closure(ValleyConstruction::f121(k, Some(0.1)))),
    ];
    let mut out: Vec<GuardEntry> = Vec::new();
    // Both framings: the viz's auto-framing, and the benches' — the one EVERY Finding from 107 to
    // 123 was measured at. A framing is not a translation (block 1), so the labels say which.
    // ADR Finding 124 -- ONE framing, the canonical one (the viz's)
    let framings: [([f64; 2], &str); 1] = [(VIZ_ORIGIN, "cadrage canonique (viz)")];
    for (origin, fname) in framings {
        for (label, kn) in states {
            let kv = Knobs { origin: Some(origin), ..kn };
            let digest = bench_eroded_digest(kv, PSEED);
            let t = Instant::now();
            let g = build_field_seed(kv, PSEED);
            let hv = field_hash(&g);
            eprintln!(
                "   2 · {label:<11} {fname}: digest {digest} · field {hv:016x} ({:.0} s)",
                t.elapsed().as_secs_f64()
            );
            if label == "livré" && origin == VIZ_ORIGIN {
                // 1 · translation: viz cell x ↔ bench cell (x + 768) mod w
                let (w, h) = (g.width, g.height);
                let shift = (VIZ_ORIGIN[0] * w as f64).round() as usize;
                let (mut same, mut maxd) = (0usize, 0f32);
                for y in 0..h {
                    for x in 0..w {
                        let a = g.data[y * w + x];
                        let b = del_b.data[y * w + (x + shift) % w];
                        if a.to_bits() == b.to_bits() {
                            same += 1;
                        } else {
                            let d = (ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres(a, &ss)
                            - ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres(b, &ss))
                        .abs();
                            maxd = maxd.max(d);
                        }
                    }
                }
                eprintln!(
                    "   1 · framing = translation? viz-origin field vs the benches' field rolled by {shift} \
                 cells: **{same} of {} cells bit-identical ({:.3} %)**, max |Δ| **{maxd:.2} m**",
                    w * h,
                    100.0 * same as f64 / (w * h) as f64
                );
            }
            if label == "C2 /10 col (défaut)" && origin == VIZ_ORIGIN {
                // 3 · the 654 km entry, at the viz's origin
                let (dr, bre) = inventory(&g, &ss);
                let agg = aggregate(&dr, &bre, false);
                let mut cand: Vec<&Wc> = agg
                    .iter()
                    .filter(|x| x.kind == SegmentKind::Watercourse && x.order == 2)
                    .collect();
                cand.sort_by(|a, b| {
                    let s = |x: &Wc| {
                        ((x.catchment_km2 - 580.0) / 580.0).abs()
                            + ((x.length_km - 654.0) / 654.0).abs()
                    };
                    s(a).total_cmp(&s(b))
                });
                eprintln!(
                    "   3 · the S2 entries closest to (580 km², 654 km) in C2 /10 AT THE VIZ'S ORIGIN:"
                );
                for x in cand.iter().take(3) {
                    let (mx, my) =
                        *dr.rivers.segments[*x.trunk.last().unwrap()].points.last().unwrap();
                    eprintln!(
                        "       {:.0} km² · trunk {:.0} km in {} segments · sink {:?} · mouth (viz frame) ({mx}, {my})",
                        x.catchment_km2,
                        x.length_km,
                        x.trunk.len(),
                        x.sink
                    );
                }
            }
            out.push(GuardEntry {
                digest,
                label: format!("{label} · {fname}"),
                field_hash: format!("{hv:016x}"),
                seed: PSEED,
                target: g.width,
                origin,
            });
        }
    }
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/bench_field_hashes.json");
    std::fs::write(&path, serde_json::to_string_pretty(&out).expect("json") + "\n").expect("write");
    eprintln!("   wrote {} ({} entries)", path.display(), out.len());
    eprintln!(
        "\n==========  end Finding 123 guard . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}

/// Path length / chord of a polyline (cells): 1 = straight. A straight line drawn in D8 already
/// reads up to 1.08 (the 22.5° staircase), so only the comparison between populations means much.
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

/// Axial R8 of chord directions at scale `l`.
fn r8_chords(segs: &[&[(u32, u32)]], l: usize) -> (f32, usize) {
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
    if n == 0 { (f32::NAN, 0) } else { (((c8 * c8 + s8 * s8).sqrt() / n as f64) as f32, n) }
}

/// A segment's OWN polyline: without the junction point that sits on its parent.
fn own(s: &ymir_core::terrain::flow::RiverSegment) -> &[(u32, u32)] {
    if s.downstream.is_some() && s.points.len() >= 2 {
        &s.points[..s.points.len() - 1]
    } else {
        &s.points[..]
    }
}

/// ADR Finding 123 — **rule 18 on the teeth, and Finding 93's control for the trunk identity.**
///
/// 1. Each C2/10 tooth (≥ 80 % wall) paired with the A1+B2 segment whose OWN mouth lies within 3
///    cells and whose own area is within a factor 1.5: sinuosity and R8 (chord 8 cells, 0.39 km)
///    per population. Straight on both sides → the D8 lines of the tectonic slope (Finding 111):
///    B does not reach them. Sinuous in A1+B2, straight on the walls → the plane: B2 is the remedy.
/// 2. The trunk identity (max-A climb) against the shipped longest-path climb, on the delivered
///    world and on C2/10, with Finding 93's own columns: trunks of ONE reach, trunk length p50 /
///    p90 / max, and the exorheic-lake crossings (Finding 45's chaining) — which must not regress.
///
/// Run: cargo test -p ymir-core --release --test f123_teeth -- --ignored f123_rule18 --nocapture
#[test]
#[ignore]
fn f123_rule18() {
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!(
        "\n==========  Finding 123 . rule 18 on the teeth, and the trunk identity  =========="
    );
    let vc = ValleyConstruction::f121(F121_AGE_K, Some(0.1));
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let sk = skeleton(&pre, &vc, &ss, DOMAIN_KM);
    let (_, masks) = carve(&pre, &sk, &vc, &ss);
    let wall: Vec<bool> = (0..n).map(|k| masks.carved[k] && !masks.floor[k]).collect();
    drop((pre, masks, sk));

    let c2 = build_field_seed(
        Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
        PSEED,
    );
    let (dc2, bc2) = inventory(&c2, &ss);
    drop(c2);
    let on = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (don, _bon) = inventory(&on, &ss);
    drop(on);
    let a_at = |dr: &C1DrainageResult, (x, y): (u32, u32)| {
        dr.flow.accumulation.data[y as usize * w + x as usize] * cell_km2
    };

    // ── 1 · the pairing ──
    let is_wc = |dr: &C1DrainageResult, i: usize| dr.segment_kind[i] == SegmentKind::Watercourse;
    let teeth: Vec<usize> = (0..dc2.rivers.segments.len())
        .filter(|&i| is_wc(&dc2, i) && dc2.rivers.segments[i].points.len() >= 2)
        .filter(|&i| {
            let p = &dc2.rivers.segments[i].points;
            p.iter().filter(|&&(x, y)| wall[y as usize * w + x as usize]).count() as f32
                >= 0.8 * p.len() as f32
        })
        .collect();
    let mut by_mouth: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, s) in don.rivers.segments.iter().enumerate() {
        if is_wc(&don, i) {
            let &(x, y) = own(s).last().unwrap();
            by_mouth.entry(y as usize * w + x as usize).or_default().push(i);
        }
    }
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for &t in &teeth {
        let st = &dc2.rivers.segments[t];
        let &(mx, my) = own(st).last().unwrap();
        let at = a_at(&dc2, (mx, my)).max(1e-6);
        let mut best: Option<(i64, usize)> = None;
        for dy in -3i64..=3 {
            for dx in -3i64..=3 {
                let (x, y) = (mx as i64 + dx, my as i64 + dy);
                if x < 0 || y < 0 || x >= w as i64 || y >= h as i64 {
                    continue;
                }
                if let Some(v) = by_mouth.get(&(y as usize * w + x as usize)) {
                    for &o in v {
                        let so = &don.rivers.segments[o];
                        let ao = a_at(&don, *own(so).last().unwrap()).max(1e-6);
                        let r = ao / at;
                        let d = dx.abs().max(dy.abs());
                        if (1.0 / 1.5..=1.5).contains(&r) && best.is_none_or(|(bd, _)| d < bd) {
                            best = Some((d, o));
                        }
                    }
                }
            }
        }
        if let Some((_, o)) = best {
            pairs.push((t, o));
        }
    }
    let p = |v: Vec<f32>, q: f64| {
        let v = sorted(v.into_iter().filter(|x| x.is_finite()).collect());
        if v.is_empty() { f32::NAN } else { pct(&v, q) }
    };
    let sin_t: Vec<f32> =
        pairs.iter().map(|&(t, _)| sinuosity(own(&dc2.rivers.segments[t]))).collect();
    let sin_o: Vec<f32> =
        pairs.iter().map(|&(_, o)| sinuosity(own(&don.rivers.segments[o]))).collect();
    let tp: Vec<&[(u32, u32)]> = pairs.iter().map(|&(t, _)| own(&dc2.rivers.segments[t])).collect();
    let op: Vec<&[(u32, u32)]> = pairs.iter().map(|&(_, o)| own(&don.rivers.segments[o])).collect();
    let (r8t, nt) = r8_chords(&tp, 8);
    let (r8o, no) = r8_chords(&op, 8);
    // references: every A1+B2 segment, and C2/10's off-wall segments
    let all_on: Vec<&[(u32, u32)]> = (0..don.rivers.segments.len())
        .filter(|&i| is_wc(&don, i))
        .map(|i| own(&don.rivers.segments[i]))
        .collect();
    let (r8a, na) = r8_chords(&all_on, 8);
    let sin_all: Vec<f32> = all_on.iter().map(|pp| sinuosity(pp)).collect();
    eprintln!(
        "\n   ── 1 · rule 18: {} teeth · paired in A1+B2 (own mouth ≤ 3 cells, own area ×/÷ 1.5): **{} \
         ({:.1} %)** ──\n   {:<34} {:>12} {:>12} {:>14}\n   {:<34} {:>12.3} {:>12.3} {:>9.4} (n {})\n   \
         {:<34} {:>12.3} {:>12.3} {:>9.4} (n {})\n   {:<34} {:>12.3} {:>12.3} {:>9.4} (n {})",
        teeth.len(),
        pairs.len(),
        100.0 * pairs.len() as f32 / teeth.len().max(1) as f32,
        "population",
        "sinuosity p50",
        "p90",
        "R8 chord 8",
        "C2/10 teeth (paired)",
        p(sin_t.clone(), 0.5),
        p(sin_t, 0.9),
        r8t,
        nt,
        "their A1+B2 pairs",
        p(sin_o.clone(), 0.5),
        p(sin_o, 0.9),
        r8o,
        no,
        "every A1+B2 segment (reference)",
        p(sin_all.clone(), 0.5),
        p(sin_all, 0.9),
        r8a,
        na
    );

    // ── 2 · Finding 93's control for the trunk identity ──
    eprintln!(
        "\n   ── 2 · trunk identity: shipped longest-path climb vs max-A climb (Finding 93's columns) ──\n   \
         {:<16} {:<12} {:>9} {:>10} {:>10} {:>10} {:>12} {:>12}",
        "world", "rule", "entries", "one-reach", "len p50", "len p90", "len max km", "lake-hops Σ"
    );
    let off = build_field_seed(Knobs::passes(2), PSEED);
    let (doff, boff) = inventory(&off, &ss);
    drop(off);
    for (name, dr, bf) in [("delivered", &doff, &boff), ("C2/10", &dc2, &bc2)] {
        let mut hop_sets: Vec<Vec<usize>> = Vec::new();
        for (rule, max_a) in [("longest", false), ("max-A", true)] {
            let agg: Vec<Wc> = aggregate(dr, bf, max_a)
                .into_iter()
                .filter(|x| x.kind == SegmentKind::Watercourse)
                .collect();
            let one = agg.iter().filter(|x| x.trunk.len() == 1).count();
            let lens: Vec<f32> = agg.iter().map(|x| x.length_km).collect();
            let hops: usize = agg.iter().map(|x| x.lake_hops).sum();
            hop_sets.push(agg.iter().map(|x| x.lake_hops).collect());
            eprintln!(
                "   {name:<16} {rule:<12} {:>9} {:>9.1}% {:>10.1} {:>10.1} {:>12.1} {:>12}",
                agg.len(),
                100.0 * one as f32 / agg.len().max(1) as f32,
                p(lens.clone(), 0.5),
                p(lens.clone(), 0.9),
                lens.iter().copied().fold(0f32, f32::max),
                hops
            );
        }
        let lost =
            hop_sets[0].iter().zip(&hop_sets[1]).filter(|(a, b)| **a > 0 && **b == 0).count();
        eprintln!(
            "   {name:<16} entries whose trunk crossed a lake under the shipped rule and NO LONGER does: \
             **{lost}**"
        );
    }
    eprintln!(
        "\n==========  end Finding 123 rule 18 . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}

/// ADR Finding 123-C2 — **the three delivered trunks that stop crossing a lake under the max-A rule.**
/// Where the two climbs diverge, and the drained area of the branch each one took.
///
/// Run: cargo test -p ymir-core --release --test f123_teeth -- --ignored f123_identity_cases --nocapture
#[test]
#[ignore]
fn f123_identity_cases() {
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 123-C2 . the trunks that stop crossing a lake  ==========");
    let off = build_field_seed(Knobs::passes(2), PSEED);
    let (dr, bf) = inventory(&off, &ss);
    drop(off);
    let segs = &dr.rivers.segments;
    let a = |i: usize| dr.segment_drainage_km2.get(i).copied().unwrap_or(0.0);
    let old = aggregate(&dr, &bf, false);
    let new = aggregate(&dr, &bf, true);
    for (x, y) in old.iter().zip(&new) {
        if x.kind != SegmentKind::Watercourse || !(x.lake_hops > 0 && y.lake_hops == 0) {
            continue;
        }
        let root = *x.trunk.last().unwrap();
        let (mx, my) = *segs[root].points.last().unwrap();
        // first divergence, walking from the mouth up
        let (ro, rn): (Vec<usize>, Vec<usize>) =
            (x.trunk.iter().rev().copied().collect(), y.trunk.iter().rev().copied().collect());
        let d = ro.iter().zip(&rn).position(|(p, q)| p != q).unwrap_or(ro.len().min(rn.len()));
        let (bo, bn) = (ro.get(d).copied(), rn.get(d).copied());
        eprintln!(
            "   · mouth ({mx}, {my}) bench frame · {:.0} km² · shipped trunk **{:.0} km, {} lake hop(s)** · \
             max-A trunk **{:.0} km, 0 hops** · they diverge {} reaches above the mouth: shipped branch \
             {:?} ({:.0} km², {} pts) vs max-A branch {:?} ({:.0} km², {} pts)",
            x.catchment_km2,
            x.length_km,
            x.lake_hops,
            y.length_km,
            d,
            bo,
            bo.map_or(0.0, a),
            bo.map_or(0, |b| segs[b].points.len()),
            bn,
            bn.map_or(0.0, a),
            bn.map_or(0, |b| segs[b].points.len())
        );
    }
    eprintln!("\n==========  end 123-C2 cases  ==========\n");
}

/// ADR Finding 123-C2 — **Finding 93's control for the max-A trunk, ties broken by the longest
/// path**, on the delivered world and on C2/10 (the identity part of `f123_rule18` alone).
///
/// Run: cargo test -p ymir-core --release --test f123_teeth -- --ignored f123_identity_control --nocapture
#[test]
#[ignore]
fn f123_identity_control() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "
==========  Finding 123-C2 . the identity control (max-A, ties → longest)  =========="
    );
    let vc = ValleyConstruction::f121(F121_AGE_K, Some(0.1));
    for (name, kn) in [
        ("delivered", Knobs::passes(2)),
        ("C2/10", Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
    ] {
        let g = build_field_seed(kn, PSEED);
        let (dr, bf) = inventory(&g, &ss);
        drop(g);
        let mut hops: Vec<Vec<usize>> = Vec::new();
        for (rule, max_a) in [("longest", false), ("max-A", true)] {
            let agg: Vec<Wc> = aggregate(&dr, &bf, max_a)
                .into_iter()
                .filter(|x| x.kind == SegmentKind::Watercourse)
                .collect();
            let one = agg.iter().filter(|x| x.trunk.len() == 1).count();
            let lens = sorted(agg.iter().map(|x| x.length_km).collect());
            hops.push(agg.iter().map(|x| x.lake_hops).collect());
            eprintln!(
                "   {name:<10} {rule:<8} entries {} · one-reach {:.1} % · length p50/p90/max {:.1} / {:.1} /                  {:.1} km · lake hops Σ {}",
                agg.len(),
                100.0 * one as f32 / agg.len().max(1) as f32,
                pct(&lens, 0.5),
                pct(&lens, 0.9),
                lens.last().copied().unwrap_or(0.0),
                hops.last().unwrap().iter().sum::<usize>()
            );
        }
        let lost = hops[0].iter().zip(&hops[1]).filter(|(a, b)| **a > 0 && **b == 0).count();
        eprintln!("   {name:<10} trunks that stop crossing a lake: **{lost}**");
    }
    eprintln!(
        "
==========  end 123-C2 control  ==========
"
    );
}

/// ADR Finding 134, item 0.2 / 0.3 — **the guard covers the lakes.** For the six states of the viz
/// menu at the canonical framing, the HD tail is run as `run_hd` runs it ([`common::viz_hd_lakes_on`]:
/// the protected breach, the placed climate 45° / 40°, the H-1 infiltration field, the assembly,
/// the C-2 crater pass), the final lakes are fingerprinted under the final drainage's digest, and
/// `data/bench_lake_hashes.json` is written. Each world's eroded field is asserted equal to the
/// field guard's reference first (the lakes are guarded on the guarded world, or not at all).
///
/// Item 0.3: the lake-by-lake listing ([`lake_listing`]) of the default state (OFF) and of the
/// extended lake base (ON) is written to `F134_LAKES_DIR` (if set), to be diffed against the viz's.
///
/// Run: cargo test -p ymir-core --release --test f123_teeth -- --ignored f134_lake_guard --nocapture
#[test]
#[ignore]
fn f134_lake_guard() {
    use common::{bench_eroded_digest, build_world, viz_hd_lakes_on};
    use ymir_core::tectonics_c1::bench_guard::{
        LakeGuardEntry, entries, field_hash, lake_fingerprint, lake_listing,
    };
    use ymir_core::tectonics_c1::valley_construction::LakeBase;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 134 . the lake guard (six states, 45° / 40°)  ==========");
    let (lat, span) = (45.0f32, 40.0f32);
    let k = F121_AGE_K;
    let with_closure = |vc: ValleyConstruction| Knobs {
        valley: Some(vc),
        slope_floor_abs: Some(S_EQ),
        origin: Some(VIZ_ORIGIN),
        ..Knobs::passes(2)
    };
    let on_vc = ValleyConstruction {
        lake_base: Some(LakeBase::InputLakesAndBasins),
        ..ValleyConstruction::new(k, Some(0.1))
    };
    let at = |kn: Knobs| Knobs { origin: Some(VIZ_ORIGIN), ..kn };
    let states: [(&str, Knobs, bool); 7] = [
        ("livré", at(Knobs::passes(2)), true),
        ("A1+B2", at(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }), true),
        ("C1 nue", at(Knobs { valley: Some(ValleyConstruction::new(k, None)), ..Knobs::passes(2) }), true),
        ("C2 /10 col (défaut)", with_closure(ValleyConstruction::new(k, Some(0.1))), true),
        ("C2 /3", with_closure(ValleyConstruction::new(k, Some(1.0 / 3.0))), true),
        ("C2 /10 niveau mer", with_closure(ValleyConstruction::f121(k, Some(0.1))), true),
        // item 0.3 only: the extended lake base, not a menu state (no field reference: NOT written)
        ("C2 /10 col, lake_base ON", with_closure(on_vc), false),
    ];
    let fields = entries();
    let dir = std::env::var("F134_LAKES_DIR").ok();
    let mut out: Vec<LakeGuardEntry> = Vec::new();
    for (label, kn, guarded) in states {
        let t = Instant::now();
        let digest = bench_eroded_digest(kn, PSEED);
        let wd = build_world(kn, None, PSEED, None);
        let fh = format!("{:016x}", field_hash(&wd.heightmap));
        if guarded {
            let e = fields.iter().find(|e| e.digest == digest).expect("the field guard has this state");
            assert_eq!(e.field_hash, fh, "{label}: the world is not the field guard's");
        }
        let v = viz_hd_lakes_on(&wd, kn, PSEED, lat, span);
        let fp = lake_fingerprint(&v.drainage.lake_map, &v.drainage.lakes);
        let craters = v.drainage.lakes.iter().filter(|l| format!("{:?}", l.lake_type).starts_with("Crater")).count();
        // ADR Finding 149-K — the lakes the tagged-outlet rule keeps Exorheic: no segment source on their border, a
        // segment tagged with their id. Exactly the lakes whose type the fix changes (the only change is there).
        {
            let dd = &v.drainage;
            let (lw, lh) = (dd.width, dd.height);
            let borders = |sx: u32, sy: u32, id: u32| {
                (-1i32..=1).any(|dy| {
                    (-1i32..=1).any(|dx| {
                        let (nx, ny) = (sx as i32 + dx, sy as i32 + dy);
                        nx >= 0 && ny >= 0 && (nx as usize) < lw && (ny as usize) < lh && dd.lake_map[ny as usize * lw + nx as usize] == id
                    })
                })
            };
            let kept: Vec<String> = dd
                .lakes
                .iter()
                .filter(|l| l.lake_type == LakeType::Exorheic)
                .filter(|l| !dd.rivers.segments.iter().any(|sg| sg.points.first().is_some_and(|&(x, y)| borders(x, y, l.base.id))))
                .map(|l| {
                    let tags: Vec<usize> = (0..dd.rivers.segments.len()).filter(|&i| dd.segment_source_lake.get(i) == Some(&Some(l.base.id))).collect();
                    format!("{} ({:.1} km², tagged segments {:?})", l.base.id, l.area_km2, tags)
                })
                .collect();
            eprintln!(
                "   {label:<26} F149-K · Exorheic only through a tagged segment: {} {:?} · Unresolved left: {}",
                kept.len(),
                kept,
                dd.lakes.iter().filter(|l| l.lake_type == LakeType::Unresolved).count()
            );
        }
        eprintln!(
            "   {label:<26} field {digest} = {fh} · lakes key {} · **{} lakes** ({craters} crater) · {fp} ({:.0} s)",
            v.digest,
            v.drainage.lakes.len(),
            t.elapsed().as_secs_f64()
        );
        if let Some(d) = &dir {
            if label == "C2 /10 col (défaut)" || !guarded {
                let name = if guarded { "bench_lakes_OFF.txt" } else { "bench_lakes_ON.txt" };
                let mut txt = format!("# {label} · lakes key {} · {fp}\n", v.digest);
                for l in lake_listing(&v.drainage.lake_map, &v.drainage.lakes) {
                    txt.push_str(&l);
                    txt.push('\n');
                }
                std::fs::write(std::path::Path::new(d).join(name), txt).expect("write listing");
            }
        }
        if guarded {
            out.push(LakeGuardEntry {
                digest: v.digest,
                label: format!("{label} · cadrage canonique (viz) · 45° / 40°"),
                lakes_hash: fp,
                lakes: v.drainage.lakes.len(),
                latitude_deg: lat,
                span_deg: span,
            });
        }
    }
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/bench_lake_hashes.json");
    std::fs::write(&path, serde_json::to_string_pretty(&out).expect("json") + "\n").expect("write");
    eprintln!("   wrote {} ({} entries)", path.display(), out.len());
    eprintln!("\n==========  end Finding 134 lake guard . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
