//! ADR Finding 156-Br — F154's amendment 2 (the breach's ramp anchored at the pit's spill), confirmed blind on three NEW
//! seeds in the production state with the criterion written before the measurement. Declared in
//! `docs/reports/lakes_gorges/f156_round3/f156_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f156_round3 -- --ignored --exact f156_br --nocapture

mod common;

use common::{Knobs, SEA, build_world, viz_dcfg, viz_hd_lakes_with};
use std::collections::{BinaryHeap, HashMap};
use std::time::Instant;
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

/// Br — the blind confirmation on new seeds (the declared criterion and stop rule).
#[test]
#[ignore]
fn f156_br() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 156 Br . the spill anchor on three new seeds (C2 /10 col, 45° / 40°)  ==========");
    let run = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0 / 64.0, dy: 1.0 / 64.0, iso_config: IsostasyConfig::c1_default(), drainage_max_distance: 30 };
    let tlf = FbmUpscaleConfig::c1_hd_production(8192).target_land_fraction;
    let mut violations: Vec<String> = Vec::new();
    for seed in SEEDS {
        let t = Instant::now();
        // the viz's automatic framing (run_hd without a manual pan)
        let land = c1_coarse_land_report(seed, 64, &Phase2InitParams::default(), &run, &C1Closures::default(), &ss, tlf);
        let (cx, cy) = land.topology.center_cell;
        let g64 = 64i64;
        let origin = [((cx as i64) - g64 / 2).rem_euclid(g64) as f64 / 64.0, ((cy as i64) - g64 / 2).rem_euclid(g64) as f64 / 64.0];
        let kn = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(S_EQ), origin: Some(origin), ..Knobs::passes(2) };
        eprintln!("\n────────── seed {seed} · origin {origin:?} · landmasses {} ──────────", land.topology.num_landmasses);
        let wd = build_world(kn, None, seed, None);
        let g = &wd.heightmap;
        let (w, h) = (g.width, g.height);
        let n = w * h;
        // the old ramps' cells, and the instrument's check against production
        let d = c1_drainage_windowed(g, None, &vd, &ss, 400.0);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let prod = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
        let (zf, ramped) = breach_rec(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), false);
        let bit = zf.data.iter().zip(&prod.data).all(|(a, b)| a.to_bits() == b.to_bits());
        drop((zf, prod, d));
        // amendment after the first run (declared): each chain in `catch_unwind`, so a panic names its chain and the next
        // seed is still measured
        use std::panic::{AssertUnwindSafe, catch_unwind};
        let msg = |e: Box<dyn std::any::Any + Send>| -> String { e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default() };
        let vp = catch_unwind(AssertUnwindSafe(|| viz_hd_lakes_with(&wd, kn, seed, 45.0, 40.0, &|a, b, c, d2, e, f, p| breach_monotone_protected(a, b, c, d2, e, f, p))));
        let vs = catch_unwind(AssertUnwindSafe(|| viz_hd_lakes_with(&wd, kn, seed, 45.0, 40.0, &|a, b, c, d2, e, f, p| breach_rec(a, b, c, d2, e, f, p, true).0)));
        let (vp, vs) = match (vp, vs) {
            (Ok(a), Ok(b)) => (a, b),
            (p, s) => {
                let why = |r: &Result<common::VizLakes, Box<dyn std::any::Any + Send>>| if r.is_ok() { "ok".to_string() } else { "PANICKED".to_string() };
                eprintln!("   the run_hd tail: production {} · spill {}", why(&p), why(&s));
                if let Err(e) = p {
                    eprintln!("      production: {}", msg(e));
                }
                if let Err(e) = s {
                    eprintln!("      spill: {}", msg(e));
                }
                violations.push(format!("seed {seed}: a chain panicked"));
                eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
                continue;
            }
        };
        let (lp, ls) = (&vp.drainage.lake_map, &vs.drainage.lake_map);
        // the pairing by footprint
        let cells_of = |lm: &[u32]| -> HashMap<u32, Vec<usize>> {
            let mut m: HashMap<u32, Vec<usize>> = HashMap::new();
            for (k, &l) in lm.iter().enumerate() {
                if l != 0 {
                    m.entry(l).or_default().push(k);
                }
            }
            m
        };
        let best = |cells: &[usize], other: &[u32]| -> Option<u32> {
            let mut ov: HashMap<u32, usize> = HashMap::new();
            for &c in cells {
                if other[c] != 0 {
                    *ov.entry(other[c]).or_insert(0) += 1;
                }
            }
            ov.into_iter().max_by_key(|e| (e.1, std::cmp::Reverse(e.0))).map(|e| e.0)
        };
        let (cp, cs) = (cells_of(lp), cells_of(ls));
        // 2. gained: a spill lake's cell outside its paired production lake
        let mut gained = 0usize;
        let mut gained_by: Vec<String> = Vec::new();
        for (sid, cells) in &cs {
            let pid = best(cells, lp);
            let g2 = cells.iter().filter(|&&c| pid.is_none_or(|p| lp[c] != p)).count();
            if g2 > 0 {
                gained += g2;
                gained_by.push(format!("spill lake {sid} (paired {pid:?}): +{g2}"));
            }
        }
        // 3. lost: a production lake's cell outside its paired spill lake
        let mut lost: Vec<usize> = Vec::new();
        for cells in cp.values() {
            let sid = best(cells, ls);
            lost.extend(cells.iter().copied().filter(|&c| sid.is_none_or(|s| ls[c] != s)));
        }
        let lost_not_r = lost.iter().filter(|&&c| !ramped[c]).count();
        let r_lake = (0..n).filter(|&k| ramped[k] && lp[k] != 0).count();
        let mut lost_set = vec![false; n];
        for &c in &lost {
            lost_set[c] = true;
        }
        let r_lake_kept = (0..n).filter(|&k| ramped[k] && lp[k] != 0 && !lost_set[k]).count();
        // 4. below-sea land cells
        let below = |c: &GridF32| -> Vec<bool> { (0..n).map(|k| g.data[k] > SEA && c.data[k] <= SEA).collect() };
        let (bp, bs) = (below(&vp.conditioned), below(&vs.conditioned));
        let new_below = (0..n).filter(|&k| bs[k] && !bp[k]).count();
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
        let (np_, ns) = (vp.drainage.lakes.len(), vs.drainage.lakes.len());
        eprintln!("   the instrument's floor copy = production, bit for bit: {bit}");
        eprintln!("   1. lakes: production {np_} · spill {ns} → {}", if np_ == ns { "unchanged" } else { "**CHANGED**" });
        eprintln!("   2. cells gained by a lake: **{gained}** {:?}", gained_by.iter().take(8).collect::<Vec<_>>());
        eprintln!(
            "   3. cells lost by the lakes **{}** · of them NOT an old ramp's cell **{lost_not_r}** · the old ramps' cells under a production lake {r_lake}, of them still under a lake (not lost) {r_lake_kept} · {}",
            lost.len(),
            if lost_not_r == 0 && r_lake_kept == 0 { "EQUAL" } else { "GAP" }
        );
        eprintln!(
            "   4. below-sea land cells: production {} · spill {} · NEW **{new_below}**",
            bp.iter().filter(|&&b| b).count(),
            bs.iter().filter(|&&b| b).count()
        );
        eprintln!("   raised {up} · lowered {dn} · river segments {sp} → {sn} ({:+.2} %)", 100.0 * (sn as f64 / sp as f64 - 1.0));
        let mut v = Vec::new();
        if np_ != ns {
            v.push("the lake count changes");
        }
        if gained > 0 {
            v.push("a lake gains cells");
        }
        if lost_not_r > 0 {
            v.push("a lost cell is not an old ramp's");
        }
        if new_below > 0 {
            v.push("a below-sea land cell is created");
        }
        if !v.is_empty() {
            violations.push(format!("seed {seed}: {}", v.join(", ")));
        }
        eprintln!("   ({:.0} s)", t.elapsed().as_secs_f64());
    }
    eprintln!("\n   STOP RULE (declared): {}", if violations.is_empty() { "no violation on any seed → production".to_string() } else { format!("FIRES: {}", violations.join(" | ")) });
    eprintln!("\n==========  end Finding 156 Br . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
