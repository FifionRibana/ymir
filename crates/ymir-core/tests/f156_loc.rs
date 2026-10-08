//! ADR Finding 156-Br, the stop rule's localisation (declared after `f156_br`, NOT blind): on the seeds that violate
//! the criterion, the lost cells that are not an old ramp's: their lakes, the lakes' levels, the cells' heights and
//! their distance to the nearest old-ramp cell.
//!
//! Run: F156_LOC_SEEDS=<seed,...> cargo test -p ymir-core --release --test f156_loc -- --ignored --exact f156_br_loc --nocapture

mod common;

use common::{Knobs, build_world, viz_dcfg, viz_hd_lakes_with};
use std::collections::{BinaryHeap, HashMap};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::c1_coarse_land_report;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
use ymir_core::tectonics_c1::drainage::c1_drainage_windowed;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone_protected};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const S_EQ: f32 = 0.024;
const SEEDS: [u64; 3] = [20_261_008_001, 20_261_008_002, 20_261_008_003];

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

/// `breach_monotone_protected` line for line (F154 amendment 2: `spill` anchors at `max(height[nb], filled[nb])`), with
/// the cells its ramp stage lowers.
#[allow(clippy::too_many_arguments)]
fn breach_rec(height: &GridF32, filled: &GridF32, lake_map: &[u32], sea_level: f32, w: usize, h: usize, protect: Option<&[bool]>, spill: bool) -> (GridF32, Vec<bool>) {
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
    let mut ramped = vec![false; n];
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
                    ramped[cur] = true;
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
    (GridF32::from_vec(w, h, z), ramped)
}

/// Br — the localisation of the lost cells outside the old ramps.
#[test]
#[ignore]
fn f156_br_loc() {
    let ss = SteinSteinParams::default();
    let m = |v: f32| c1_altitude_norm_to_metres(v, &ss);
    let vd = viz_dcfg();
    let run = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0 / 64.0, dy: 1.0 / 64.0, iso_config: IsostasyConfig::c1_default(), drainage_max_distance: 30 };
    let tlf = FbmUpscaleConfig::c1_hd_production(8192).target_land_fraction;
    let seeds: Vec<u64> = std::env::var("F156_LOC_SEEDS").unwrap_or_else(|_| SEEDS[0].to_string()).split(',').map(|s| s.trim().parse().expect("a seed")).collect();
    eprintln!("\n==========  Finding 156 Br . the lost cells outside the old ramps, localised  ==========");
    for seed in seeds {
        let land = c1_coarse_land_report(seed, 64, &Phase2InitParams::default(), &run, &C1Closures::default(), &ss, tlf);
        let (cx, cy) = land.topology.center_cell;
        let origin = [((cx as i64) - 32).rem_euclid(64) as f64 / 64.0, ((cy as i64) - 32).rem_euclid(64) as f64 / 64.0];
        let kn = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(S_EQ), origin: Some(origin), ..Knobs::passes(2) };
        eprintln!("\n────────── seed {seed} ──────────");
        let wd = build_world(kn, None, seed, None);
        let g = &wd.heightmap;
        let (w, h) = (g.width, g.height);
        let n = w * h;
        let d = c1_drainage_windowed(g, None, &vd, &ss, 400.0);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let (_, ramped) = breach_rec(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), false);
        let pre_lake: Vec<bool> = d.lake_map.iter().map(|&l| l != 0).collect();
        drop(d);
        let vp = viz_hd_lakes_with(&wd, kn, seed, 45.0, 40.0, &|a, b, c, d2, e, f, p| breach_monotone_protected(a, b, c, d2, e, f, p));
        let vs = viz_hd_lakes_with(&wd, kn, seed, 45.0, 40.0, &|a, b, c, d2, e, f, p| breach_rec(a, b, c, d2, e, f, p, true).0);
        let (lp, ls) = (&vp.drainage.lake_map, &vs.drainage.lake_map);
        let lake_of = |v: &common::VizLakes, id: u32| v.drainage.lakes.iter().find(|l| l.base.id == id).map(|l| (format!("{:?}", l.lake_type), l.level_m, l.area_km2));
        let mut cp: HashMap<u32, Vec<usize>> = HashMap::new();
        for (k, &l) in lp.iter().enumerate() {
            if l != 0 {
                cp.entry(l).or_default().push(k);
            }
        }
        // the distance (Chebyshev, cells) to the nearest old-ramp cell
        let mut dist = vec![u32::MAX; n];
        let mut q = std::collections::VecDeque::new();
        for k in 0..n {
            if ramped[k] {
                dist[k] = 0;
                q.push_back(k);
            }
        }
        while let Some(c) = q.pop_front() {
            for kk in 0..8 {
                let x = ((c % w) as i32 + D8_DX[kk]).rem_euclid(w as i32) as usize;
                let y = ((c / w) as i32 + D8_DY[kk]).rem_euclid(h as i32) as usize;
                let mm = y * w + x;
                if dist[mm] == u32::MAX {
                    dist[mm] = dist[c] + 1;
                    q.push_back(mm);
                }
            }
        }
        let mut ids: Vec<&u32> = cp.keys().collect();
        ids.sort();
        for &pid in ids {
            let cells = &cp[&pid];
            let mut ov: HashMap<u32, usize> = HashMap::new();
            for &c in cells {
                if ls[c] != 0 {
                    *ov.entry(ls[c]).or_insert(0) += 1;
                }
            }
            let sid = ov.into_iter().max_by_key(|e| (e.1, std::cmp::Reverse(e.0))).map(|e| e.0);
            let lost: Vec<usize> = cells.iter().copied().filter(|&c| sid.is_none_or(|s| ls[c] != s)).collect();
            let off: Vec<usize> = lost.iter().copied().filter(|&c| !ramped[c]).collect();
            if off.is_empty() {
                continue;
            }
            let mut dd: Vec<u32> = off.iter().map(|&c| dist[c]).collect();
            dd.sort_unstable();
            let mut zz: Vec<f32> = off.iter().map(|&c| m(g.data[c])).collect();
            zz.sort_by(f32::total_cmp);
            let mut cz: Vec<f32> = off.iter().map(|&c| m(vs.conditioned.data[c]) - m(vp.conditioned.data[c])).collect();
            cz.sort_by(f32::total_cmp);
            let (bx0, bx1, by0, by1) = (off.iter().map(|&c| c % w).min().unwrap(), off.iter().map(|&c| c % w).max().unwrap(), off.iter().map(|&c| c / w).min().unwrap(), off.iter().map(|&c| c / w).max().unwrap());
            eprintln!(
                "   production lake {pid} {:?} → spill lake {sid:?} {:?} · cells {} · lost {} · lost NOT an old ramp's **{}** · their distance to the nearest ramp cell (cells) p0 {} p50 {} p100 {} · their eroded height (m) p0 {:.2} p50 {:.2} p100 {:.2} · their conditioned Δz (spill − production, m) p0 {:+.3} p50 {:+.3} p100 {:+.3} · in a pre-breach lake {} · bbox x {bx0}–{bx1} y {by0}–{by1}",
                lake_of(&vp, pid),
                sid.and_then(|s| lake_of(&vs, s)),
                cells.len(),
                lost.len(),
                off.len(),
                dd[0],
                dd[dd.len() / 2],
                dd[dd.len() - 1],
                zz[0],
                zz[zz.len() / 2],
                zz[zz.len() - 1],
                cz[0],
                cz[cz.len() / 2],
                cz[cz.len() - 1],
                off.iter().filter(|&&c| pre_lake[c]).count()
            );
        }
    }
}
