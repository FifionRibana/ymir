//! ADR Finding 154 — lakes / gorges / falls resumed, round 1 of 4: S (the geology's substratum / surface split, its
//! gate on the témoin's rock grid), H (falls by rock, three rules tabulated), L (P4: the light pass not applied on the
//! kept bodies' catchments) and Br (a bench copy of the breach whose ramp aims at the overflow level). Declared in
//! `docs/reports/lakes_gorges/f154_resume/f154_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f154_resume -- --ignored --exact <name> --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, build_world, viz_hd_lakes_on};
use std::time::Instant;
use ymir_core::climate::precipitation::precip_mm_per_year;
use ymir_core::geology::rocks::ROCK_CLASSES;
use ymir_core::geology::{GeologyInputs, build_geology};
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const S_EQ: f32 = 0.024;

fn temoin() -> Knobs {
    Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }
}

fn c1_run() -> C1TimeLoopConfig {
    C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0 / 64.0, dy: 1.0 / 64.0, iso_config: IsostasyConfig::c1_default(), drainage_max_distance: 30 }
}

/// FNV-1a 64 over bytes.
fn fnv(b: &[u8]) -> u64 {
    let mut x: u64 = 0xcbf2_9ce4_8422_2325;
    for &v in b {
        x ^= v as u64;
        x = x.wrapping_mul(0x0000_0100_0000_01b3);
    }
    x
}

/// S — the témoin's final rock grid: its hash and class counts (the gate), and the build's timings.
#[test]
#[ignore]
fn f154_s() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 154 S . the témoin's rock grid  ==========");
    let wd = build_world(temoin(), None, PSEED, None);
    let v = viz_hd_lakes_on(&wd, temoin(), PSEED, 45.0, 40.0);
    let f = &v.conditioned;
    let z_m: Vec<f32> = f.data.iter().map(|&x| c1_altitude_norm_to_metres(x, &ss)).collect();
    let precip_mm: Vec<f32> = v.precipitation.data.iter().map(|&p| precip_mm_per_year(p)).collect();
    let run = c1_run();
    let tb = Instant::now();
    let product = build_geology(&GeologyInputs {
        seed: PSEED,
        grid: 64,
        init: &Phase2InitParams::default(),
        run: &run,
        closures: &C1Closures::default(),
        field: f,
        z_m: &z_m,
        cell_km: CELL_KM,
        drainage: &v.drainage,
        wetland: &v.wetland,
        temp_c: &v.temperature.data,
        precip_mm: &precip_mm,
        sample_origin: wd.cfg.sample_origin,
        sample_size: wd.cfg.sample_size,
        domain_km: 400.0,
        volcanism: Some(&wd.volc),
    });
    let t = product.timings;
    let mut counts = [0usize; 9];
    for &r in &product.rocks {
        counts[r as usize] += 1;
    }
    eprintln!(
        "   ROCK GRID · fnv1a64 **{:016x}** · {} cells · counts {:?}",
        fnv(&product.rocks),
        product.rocks.len(),
        ROCK_CLASSES.iter().map(|c| format!("{} {}", c.key, counts[c.id as usize])).collect::<Vec<_>>()
    );
    // the HEAD run (aa8ca60, F152's rocks code) printed dbd91290b6c6148c
    eprintln!("   S GATE · the rock grid equals HEAD's (F152's) dbd91290b6c6148c: **{}**", fnv(&product.rocks) == 0xdbd9_1290_b6c6_148c);
    eprintln!("   build_geology {:.2} s: history {:.2} · rocks {:.2} · context {:.2}", tb.elapsed().as_secs_f64(), t.history_s, t.rocks_s, t.context_s);
    eprintln!("\n==========  end Finding 154 S . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// Br, amendment AFTER the run (declared non-blind in `f154_declared.md`): ON's below-sea land cells under both anchors,
/// attributed. An instrumented copy of the breach records, for every land cell a ramp takes to ≤ sea, its ramp (the pit,
/// the outlet cell, the anchor), its step and its original height above the sea; cells the ramps never touch are
/// "other". Run on ON, production anchor (floor) and overflow anchor.
#[test]
#[ignore]
fn f154_br_attr() {
    use common::{SEA, build_world, viz_dcfg};
    use std::collections::{BTreeMap, BinaryHeap};
    use ymir_core::grid::GridF32;
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::c1_drainage_windowed;
    use ymir_core::tectonics_c1::valley_construction::LakeBase;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone_protected};
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
    /// The copy, recording the first ramp to take each cell to ≤ sea: (pit, outlet cell ci, anchor, step).
    #[allow(clippy::type_complexity)]
    fn breach_rec(height: &GridF32, filled: &GridF32, lake_map: &[u32], sea_level: f32, w: usize, h: usize, protect: Option<&[bool]>, overflow: bool) -> (Vec<f32>, Vec<Option<(usize, usize, f32, u32)>>) {
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
        let mut rec: Vec<Option<(usize, usize, f32, u32)>> = vec![None; n];
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
                    let anchor = if overflow { z[ci] } else { height.data[nb] };
                    let mut target = anchor - EPS;
                    let mut cur = ci;
                    let mut step = 0u32;
                    while cur != usize::MAX && !is_base(cur, &z) && z[cur] > target {
                        z[cur] = target;
                        if target <= sea_level && rec[cur].is_none() {
                            rec[cur] = Some((nb, ci, anchor, step));
                        }
                        target -= EPS;
                        step += 1;
                        cur = backlink[cur];
                    }
                }
                z[nb] = height.data[nb];
                heap.push(Pq(z[nb], nb));
            }
        }
        (z, rec)
    }
    let ss = SteinSteinParams::default();
    let m = |v: f32| c1_altitude_norm_to_metres(v, &ss);
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 154 Br, amendment: ON's below-sea land cells attributed  ==========");
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let k = Knobs { valley: Some(on), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(k, None, PSEED, None);
    let g = &wd.heightmap;
    let (w, h) = (g.width, g.height);
    let d = c1_drainage_windowed(g, None, &vd, &ss, 400.0);
    let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
    let prod = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
    for overflow in [false, true] {
        let (zr, rec) = breach_rec(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), overflow);
        let below: Vec<usize> = (0..w * h).filter(|&k| g.data[k] > SEA && zr[k] <= SEA).collect();
        if !overflow {
            let bp: Vec<usize> = (0..w * h).filter(|&k| g.data[k] > SEA && prod.data[k] <= SEA).collect();
            eprintln!("   the ramp stage's below-sea set equals production's (the mop-up fill raises no base): {}", bp == below);
        }
        let by_ramp = below.iter().filter(|&&k| rec[k].is_some()).count();
        let mut ramps: BTreeMap<(usize, usize), (f32, f32, f32, u32, usize)> = BTreeMap::new(); // (pit, ci) → (pit floor, anchor, ci's z, max step, cells)
        let mut above: Vec<f32> = Vec::new();
        for &c in &below {
            above.push(m(g.data[c]));
            if let Some((pit, ci, anchor, step)) = rec[c] {
                let e = ramps.entry((pit, ci)).or_insert((m(g.data[pit]), m(anchor), m(d.flow.filled.data[ci]), 0, 0));
                e.3 = e.3.max(step);
                e.4 += 1;
            }
        }
        above.sort_by(f32::total_cmp);
        let pc = |p: f32| if above.is_empty() { f32::NAN } else { above[((above.len() - 1) as f32 * p) as usize] };
        eprintln!(
            "\n   anchor {} · below-sea land cells **{}** · taken by a ramp {by_ramp} · their original height above the sea (m): p0 {:.3} · p50 {:.3} · p100 {:.3} · ramps {}",
            if overflow { "OVERFLOW (z[ci])" } else { "FLOOR (production)" },
            below.len(),
            pc(0.0),
            pc(0.5),
            pc(1.0),
            ramps.len()
        );
        let mut rv: Vec<_> = ramps.into_iter().collect();
        rv.sort_by_key(|e| std::cmp::Reverse(e.1.4));
        for ((pit, ci), (pf, an, cf, st, cells)) in rv.iter().take(12) {
            eprintln!(
                "      ramp from pit ({},{}) floor {pf:.2} m via ({},{}) [filled {cf:.2} m] · anchor {an:.2} m · {cells} cells ≤ sea · steps to the deepest {st}",
                pit % w,
                pit / w,
                ci % w,
                ci / w
            );
        }
    }
}

/// Br, second amendment AFTER the run (declared non-blind in `f154_declared.md`): the ramp anchored at the pit's SPILL on
/// the pre-breach flood (`filled[nb]`), since `f154_br_attr` showed that `z[ci]` is already lowered by the neighbouring
/// pits' ramps. On ON (below-sea cells, cells, lakes, rivers) and on the gorge P4 ×1.4 p .5 (G-sea).
#[test]
#[ignore]
fn f154_br_spill() {
    use common::{SEA, build_world, viz_dcfg, viz_hd_lakes_on, viz_hd_lakes_with};
    use std::collections::{BinaryHeap, HashMap, HashSet};
    use ymir_core::grid::GridF32;
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::c1_drainage_windowed;
    use ymir_core::tectonics_c1::valley_construction::{GorgeRetreat, LakeBase};
    use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone_protected};
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
    /// `breach_monotone_protected` line for line; `spill` anchors the ramp at `max(height[nb], filled[nb])`.
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
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let vd = viz_dcfg();
    eprintln!("\n==========  Finding 154 Br, amendment 2: the ramp anchored at the pit's spill (filled[nb])  ==========");
    let on_at = |age: f32| ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K * age, Some(0.1)) };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let wd = build_world(kn(on_at(1.0)), None, PSEED, None);
    let g = &wd.heightmap;
    let (w, h) = (g.width, g.height);
    let n = w * h;
    let below = |g: &GridF32, c: &GridF32| -> Vec<usize> { (0..n).filter(|&k| g.data[k] > SEA && c.data[k] <= SEA).collect() };
    let d = c1_drainage_windowed(g, None, &vd, &ss, 400.0);
    let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
    let prod = breach_monotone_protected(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref());
    let same = breach_spill(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), false);
    eprintln!("   the copy with the change OFF is bit-identical to production: {}", prod.data.iter().zip(&same.data).all(|(a, b)| a.to_bits() == b.to_bits()));
    drop(same);
    let sea_on: HashSet<usize> = below(g, &prod).into_iter().collect();
    drop(d);
    let vp = viz_hd_lakes_on(&wd, kn(on_at(1.0)), PSEED, 45.0, 40.0);
    let vs = viz_hd_lakes_with(&wd, kn(on_at(1.0)), PSEED, 45.0, 40.0, &|a, b, c, d2, e, f, p| breach_spill(a, b, c, d2, e, f, p, true));
    let bs = below(g, &vs.conditioned);
    eprintln!(
        "   ON · below-sea land cells: production **{}** (covered {}) · spill copy **{}** (covered {})",
        sea_on.len(),
        sea_on.iter().filter(|&&k| vp.drainage.lake_map[k] != 0).count(),
        bs.len(),
        bs.iter().filter(|&&k| vs.drainage.lake_map[k] != 0).count()
    );
    let (mut up, mut dn) = (0usize, 0usize);
    for k in 0..n {
        let dz = (vs.conditioned.data[k] - vp.conditioned.data[k]) * n2m;
        if dz > 0.01 {
            up += 1;
        } else if dz < -0.01 {
            dn += 1;
        }
    }
    eprintln!("   ON · conditioned cells changed > 1 cm: **{}** (raised {up}, lowered {dn})", up + dn);
    let cells_of = |lm: &[u32]| -> HashMap<u32, Vec<usize>> {
        let mut m: HashMap<u32, Vec<usize>> = HashMap::new();
        for (k, &l) in lm.iter().enumerate() {
            if l != 0 {
                m.entry(l).or_default().push(k);
            }
        }
        m
    };
    let best = |cells: &[usize], other: &[u32]| -> Option<(u32, usize)> {
        let mut ov: HashMap<u32, usize> = HashMap::new();
        for &c in cells {
            if other[c] != 0 {
                *ov.entry(other[c]).or_insert(0) += 1;
            }
        }
        ov.into_iter().max_by_key(|e| (e.1, std::cmp::Reverse(e.0)))
    };
    let (cp, cs) = (cells_of(&vp.drainage.lake_map), cells_of(&vs.drainage.lake_map));
    let lv = |v: &common::VizLakes| -> HashMap<u32, (f32, f32)> { v.drainage.lakes.iter().map(|l| (l.base.id, (l.level_m, l.area_km2))).collect() };
    let (lp, ls) = (lv(&vp), lv(&vs));
    let (mut gone, mut moved, mut same_n) = (0usize, 0usize, 0usize);
    let mut ids: Vec<&u32> = cp.keys().collect();
    ids.sort();
    for &id in &ids {
        let cells = &cp[id];
        let (l0, a0) = lp.get(id).copied().unwrap_or((f32::NAN, f32::NAN));
        match best(cells, &vs.drainage.lake_map).filter(|&(_, c)| 2 * c > cells.len()) {
            None => {
                gone += 1;
                eprintln!("      lake {id} ({a0:.2} km², {l0:.1} m): GONE");
            }
            Some((sid, _)) => {
                let (l1, a1) = ls.get(&sid).copied().unwrap_or((f32::NAN, f32::NAN));
                if (l1 - l0).abs() > 1.0 || (a1 / a0 - 1.0).abs() > 0.05 {
                    moved += 1;
                    eprintln!("      lake {id}: level {l0:.1} → {l1:.1} m · area {a0:.2} → {a1:.2} km²");
                } else {
                    same_n += 1;
                }
            }
        }
    }
    let new_n = cs.values().filter(|cells| best(cells, &vp.drainage.lake_map).is_none_or(|(_, c)| 2 * c <= cells.len())).count();
    eprintln!("   ON · lakes: production {} · copy {} · unchanged {same_n} · gone {gone} · moved (level > 1 m or area > 5 %) {moved} · new {new_n}", cp.len(), cs.len());
    let rset = |v: &common::VizLakes| -> HashSet<(u32, u32)> { v.drainage.rivers.segments.iter().flat_map(|s| s.points.iter().copied()).collect() };
    let (rp, rs) = (rset(&vp), rset(&vs));
    eprintln!(
        "   ON · rivers: segments {} → {} ({:+.2} %) · river cells in one only: {} (production only {}, copy only {})",
        vp.drainage.rivers.segments.len(),
        vs.drainage.rivers.segments.len(),
        100.0 * (vs.drainage.rivers.segments.len() as f64 / vp.drainage.rivers.segments.len() as f64 - 1.0),
        rp.symmetric_difference(&rs).count(),
        rp.difference(&rs).count(),
        rs.difference(&rp).count()
    );
    drop((vp, vs, wd));
    // the gorge P4 ×1.4 p .5: G-sea against ON's production set
    let g5 = GorgeRetreat { light_mode: 4, ..GorgeRetreat::v5(2.0, 0.5) };
    let vc = ValleyConstruction { gorge_retreat: Some(g5), ..on_at(1.4) };
    let wg = build_world(kn(vc), None, PSEED, None);
    let gg = &wg.heightmap;
    let d = c1_drainage_windowed(gg, None, &vd, &ss, 400.0);
    let prot = wg.volc.enabled.then(|| crater_protect_mask(&wg.craters, w, h));
    for spill in [false, true] {
        let c = breach_spill(gg, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), spill);
        let s: HashSet<usize> = below(gg, &c).into_iter().collect();
        eprintln!("   GORGE P4 ×1.4 p .5 · {} · below-sea land cells {} · G-sea +{}/−{}", if spill { "spill copy" } else { "production (the copy, change off)" }, s.len(), s.difference(&sea_on).count(), sea_on.difference(&s).count());
    }
}
