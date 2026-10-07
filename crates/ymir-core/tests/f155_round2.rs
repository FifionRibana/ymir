//! ADR Finding 155 — lakes / gorges / falls resumed, round 2 of 4, Part Br: F154's amendment 2 (the breach's ramp
//! anchored at the pit's spill), confirmed blind on the guard's six states before any production change. Declared in
//! `docs/reports/lakes_gorges/f155_round2/f155_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f155_round2 -- --ignored --exact f155_br --nocapture

mod common;

use common::{Knobs, PSEED, SEA, build_world, viz_hd_lakes_on, viz_hd_lakes_with};
use std::collections::BinaryHeap;
use std::time::Instant;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::bench_guard::{lake_entries, lake_fingerprint};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};
use ymir_core::terrain::flow::{D8_DX, D8_DY};

const S_EQ: f32 = 0.024;
const VIZ_ORIGIN: [f64; 2] = [0.09375, 0.578_125];

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

/// F154's amendment-2 copy of `breach_monotone_protected`, untouched: `spill` anchors the ramp at
/// `max(height[nb], filled[nb])`.
#[allow(clippy::too_many_arguments)]
fn breach_spill(height: &GridF32, filled: &GridF32, lake_map: &[u32], sea_level: f32, w: usize, h: usize, protect: Option<&[bool]>, spill: bool) -> GridF32 {
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
                let mut target = if spill { height.data[nb].max(filled.data[nb]) - EPS } else { height.data[nb] - EPS };
                let mut cur = ci;
                while cur != usize::MAX && !is_base(cur, &z) && z[cur] > target {
                    z[cur] = target;
                    target -= EPS;
                    cur = backlink[cur];
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
    GridF32::from_vec(w, h, z)
}

/// Br — the blind confirmation on the guard's six states (the declared stop rule).
#[test]
#[ignore]
fn f155_br() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 155 Br . the spill anchor on the guard's six states (45° / 40°)  ==========");
    let k = F121_AGE_K;
    let with_closure = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), origin: Some(VIZ_ORIGIN), ..Knobs::passes(2) };
    let at = |kn: Knobs| Knobs { origin: Some(VIZ_ORIGIN), ..kn };
    let states: [(&str, Knobs); 6] = [
        ("livré", at(Knobs::passes(2))),
        ("A1+B2", at(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) })),
        ("C1 nue", at(Knobs { valley: Some(ValleyConstruction::new(k, None)), ..Knobs::passes(2) })),
        ("C2 /10 col (défaut, the témoin)", with_closure(ValleyConstruction::new(k, Some(0.1)))),
        ("C2 /3", with_closure(ValleyConstruction::new(k, Some(1.0 / 3.0)))),
        ("C2 /10 niveau mer", with_closure(ValleyConstruction::f121(k, Some(0.1)))),
    ];
    let guard = lake_entries();
    let mut fired: Vec<String> = Vec::new();
    for (label, kn) in states {
        let t = Instant::now();
        eprintln!("\n────────── {label} ──────────");
        let wd = build_world(kn, None, PSEED, None);
        let g = &wd.heightmap;
        let (w, h) = (g.width, g.height);
        let n = w * h;
        let vp = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
        let vs = viz_hd_lakes_with(&wd, kn, PSEED, 45.0, 40.0, &|a, b, c, d, e, f, p| breach_spill(a, b, c, d, e, f, p, true));
        let below = |c: &GridF32| -> Vec<bool> { (0..n).map(|k| g.data[k] > SEA && c.data[k] <= SEA).collect() };
        let (bp, bs) = (below(&vp.conditioned), below(&vs.conditioned));
        let new_below: Vec<usize> = (0..n).filter(|&k| bs[k] && !bp[k]).collect();
        let fp = lake_fingerprint(&vp.drainage.lake_map, &vp.drainage.lakes);
        let fs = lake_fingerprint(&vs.drainage.lake_map, &vs.drainage.lakes);
        let (mp, lp) = fp.split_once('/').expect("two halves");
        let (ms, ls) = fs.split_once('/').expect("two halves");
        let stored = guard.iter().find(|e| e.digest == vp.digest).map(|e| e.lakes_hash.clone());
        let (mut up, mut dn) = (0usize, 0usize);
        for k in 0..n {
            let dz = (vs.conditioned.data[k] - vp.conditioned.data[k]) * n2m;
            if dz > 0.01 {
                up += 1;
            } else if dz < -0.01 {
                dn += 1;
            }
        }
        let (sp, sn) = (vp.drainage.rivers.segments.len(), vs.drainage.rivers.segments.len());
        eprintln!(
            "   below-sea land cells: production **{}** · spill **{}** · NEW under the spill **{}**",
            bp.iter().filter(|&&b| b).count(),
            bs.iter().filter(|&&b| b).count(),
            new_below.len()
        );
        eprintln!(
            "   lakes: production {} · spill {} · footprint (mask hash) {mp} → {ms} **{}** · list hash {lp} → {ls} {} · production = the guard's stored hash: {}",
            vp.drainage.lakes.len(),
            vs.drainage.lakes.len(),
            if mp == ms { "SAME" } else { "DIFFERS" },
            if lp == ls { "SAME" } else { "DIFFERS" },
            stored.as_deref().map_or("no entry for this digest".to_string(), |s| (s == fp).to_string())
        );
        eprintln!("   conditioned cells > 1 cm: raised **{up}** · lowered **{dn}** · river segments {sp} → {sn} ({:+.2} %)", 100.0 * (sn as f64 / sp as f64 - 1.0));
        if lp != ls {
            // where the list differs: lake by lake
            let lv = |v: &common::VizLakes| -> Vec<(u32, String, u32, u32)> { v.drainage.lakes.iter().map(|l| (l.base.id, format!("{:?}", l.lake_type), l.level_m.to_bits(), l.area_km2.to_bits())).collect() };
            let (a, b) = (lv(&vp), lv(&vs));
            for (x, y) in a.iter().zip(b.iter()).filter(|(x, y)| x != y).take(12) {
                eprintln!(
                    "      list differs: {} {} {:.4} m {:.5} km² → {} {} {:.4} m {:.5} km²",
                    x.0,
                    x.1,
                    f32::from_bits(x.2),
                    f32::from_bits(x.3),
                    y.0,
                    y.1,
                    f32::from_bits(y.2),
                    f32::from_bits(y.3)
                );
            }
        }
        let fires = vp.drainage.lakes.len() != vs.drainage.lakes.len() || mp != ms || !new_below.is_empty();
        if fires {
            let mut why = Vec::new();
            if vp.drainage.lakes.len() != vs.drainage.lakes.len() || mp != ms {
                why.push("a lake changes (count or footprint)".to_string());
            }
            if !new_below.is_empty() {
                why.push(format!("{} new below-sea land cells, first {:?}", new_below.len(), new_below.iter().take(5).map(|&c| (c % w, c / w)).collect::<Vec<_>>()));
            }
            fired.push(format!("{label}: {}", why.join("; ")));
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n   STOP RULE (declared): {}", if fired.is_empty() { "does NOT fire on any state → production".to_string() } else { format!("FIRES: {}", fired.join(" | ")) });
    eprintln!("\n==========  end Finding 155 Br . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
