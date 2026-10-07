//! ADR Finding 155-F — the soft-lip rule on the témoin: v5 against v6 at ×1 p .5, the construction (every φ, the
//! tagged falls, the constructed field) and the production path (`build_world`'s field), plus the soft mask against
//! F154's substratum. Declared in `docs/reports/lakes_gorges/f155_round2/f155_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f155_soft -- --ignored --exact f155_f --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, build_world};
use std::time::Instant;
use ymir_core::geology::rocks::ROCK_CLASSES;
use ymir_core::geology::tectonic_sources;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, GorgeRetreat, LakeBase, ValleyConstruction, carve, skeleton, skeleton_with_soft_lips};

const S_EQ: f32 = 0.024;
const DOMAIN_KM: f32 = 400.0;

fn fnv(data: &[f32]) -> u64 {
    let mut x: u64 = 0xcbf2_9ce4_8422_2325;
    for v in data {
        for b in v.to_bits().to_le_bytes() {
            x ^= b as u64;
            x = x.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    x
}

#[test]
#[ignore]
fn f155_f() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 155 F . the soft-lip rule on the témoin (×1, p 0.5)  ==========");
    let on = ValleyConstruction { lake_base: Some(LakeBase::InputLakesAndBasins), ..ValleyConstruction::new(F121_AGE_K, Some(0.1)) };
    let v5 = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v5(1.0, 0.5)), ..on };
    let v6 = ValleyConstruction { gorge_retreat: Some(GorgeRetreat::v6(1.0, 0.5)), ..on };
    let kn = |vc: ValleyConstruction| Knobs { valley: Some(vc), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let (w, h) = (s1.width, s1.height);
    let n = w * h;
    let zs1: Vec<f32> = s1.data.iter().map(|&v| c1_altitude_norm_to_metres(v, &ss)).collect();
    // the production path first (its world gives the cfg and the volcanism)
    let tb = Instant::now();
    let w5 = build_world(kn(v5), None, PSEED, None);
    let t5 = tb.elapsed().as_secs_f64();
    let tb = Instant::now();
    let w6 = build_world(kn(v6), None, PSEED, None);
    let t6 = tb.elapsed().as_secs_f64();
    eprintln!(
        "   PRODUCTION PATH · build_world v5 {:016x} ({t5:.0} s) · v6 {:016x} ({t6:.0} s) · identical: **{}** · tagged falls {} / {}",
        fnv(&w5.heightmap.data),
        fnv(&w6.heightmap.data),
        w5.heightmap.data.iter().zip(&w6.heightmap.data).all(|(a, b)| a.to_bits() == b.to_bits()),
        w5.gorge_falls.len(),
        w6.gorge_falls.len()
    );
    // the soft cells from F154's substratum on S1
    let run = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0 / 64.0, dy: 1.0 / 64.0, iso_config: IsostasyConfig::c1_default(), drainage_max_distance: 30 };
    let ts = Instant::now();
    let src = tectonic_sources(PSEED, 64, &Phase2InitParams::default(), &run, &C1Closures::default(), Some(&w5.volc));
    let sub = src.substratum(&zs1, w, h, CELL_KM, w5.cfg.sample_origin, w5.cfg.sample_size);
    let t_soft = ts.elapsed().as_secs_f64();
    let soft: Vec<bool> = sub.iter().map(|&r| ROCK_CLASSES[r as usize].hardness == 0).collect();
    let land = (0..n).filter(|&k| s1.data[k] > SEA).count();
    let soft_land = (0..n).filter(|&k| s1.data[k] > SEA && soft[k]).count();
    eprintln!("   soft substratum cells (rift fill) on S1: {soft_land} land cells ({:.2} % of the land) · sources + substratum {t_soft:.2} s", 100.0 * soft_land as f64 / land as f64);
    // the construction, v5 against v6 with the soft cells
    let s5 = skeleton(&s1, &v5, &ss, DOMAIN_KM);
    let s6 = skeleton_with_soft_lips(&s1, &v6, &ss, DOMAIN_KM, Some(&soft));
    let mut same_phi = true;
    let mut soft_lips = 0;
    for (a, b) in s5.gorge_bodies.iter().zip(&s6.gorge_bodies) {
        let lip_soft = a.col != u32::MAX && soft[a.col as usize];
        soft_lips += lip_soft as usize;
        same_phi &= a.phi.to_bits() == b.phi.to_bits() && a.low == b.low;
        eprintln!(
            "      body ({},{}) · lip {} · φ v5 {:.3} · v6 {:.3}",
            a.low as usize % w,
            a.low as usize / w,
            if a.col == u32::MAX { "none".to_string() } else { ROCK_CLASSES[sub[a.col as usize] as usize].key.to_string() },
            a.phi,
            b.phi
        );
    }
    let (c5, _) = carve(&s1, &s5, &v5, &ss);
    let (c6, _) = carve(&s1, &s6, &v6, &ss);
    eprintln!(
        "   CONSTRUCTION · bodies {} / {} · soft lips **{soft_lips}** · every φ identical **{same_phi}** · falls {} / {} identical {} · carve v5 {:016x} · v6 {:016x} · identical **{}**",
        s5.gorge_bodies.len(),
        s6.gorge_bodies.len(),
        s5.gorge_falls.len(),
        s6.gorge_falls.len(),
        s5.gorge_falls == s6.gorge_falls,
        fnv(&c5.data),
        fnv(&c6.data),
        c5.data.iter().zip(&c6.data).all(|(a, b)| a.to_bits() == b.to_bits())
    );
    eprintln!("\n==========  end Finding 155 F . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// Br, the stop rule's localisation (declared in `f155_declared.md`'s amendment, after `f155_br`): on the states the rule
/// fired on, the lakes whose footprint or listing differ (`lake_listing`, line by line) and the land cells the spill
/// anchor still leaves below the sea, attributed to their ramps (the anchor's height above the sea, the ramp's length).
///
/// Run: cargo test -p ymir-core --release --test f155_soft -- --ignored --exact f155_br_loc --nocapture
#[test]
#[ignore]
fn f155_br_loc() {
    use common::{viz_hd_lakes_on, viz_hd_lakes_with};
    use std::collections::{BTreeMap, BinaryHeap, HashMap};
    use ymir_core::grid::GridF32;
    use ymir_core::tectonics_c1::bench_guard::lake_listing;
    use ymir_core::tectonics_c1::closures::volcanism::crater_protect_mask;
    use ymir_core::tectonics_c1::drainage::c1_drainage_windowed;
    use ymir_core::terrain::flow::{D8_DX, D8_DY};
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
    /// F154's amendment-2 copy (`spill`), recording for each cell a ramp takes to ≤ sea: (pit, anchor, step); with the
    /// mop-up fill, so its output is the copy's.
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    fn breach_rec(height: &GridF32, filled: &GridF32, lake_map: &[u32], sea_level: f32, w: usize, h: usize, protect: Option<&[bool]>, spill: bool) -> (GridF32, Vec<Option<(usize, f32, u32)>>) {
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
        let mut rec: Vec<Option<(usize, f32, u32)>> = vec![None; n];
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
                    let anchor = if spill { height.data[nb].max(filled.data[nb]) } else { height.data[nb] };
                    let mut target = anchor - EPS;
                    let mut cur = ci;
                    let mut step = 0u32;
                    while cur != usize::MAX && !is_base(cur, &z) && z[cur] > target {
                        z[cur] = target;
                        if target <= sea_level && rec[cur].is_none() {
                            rec[cur] = Some((nb, anchor, step));
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
        (GridF32::from_vec(w, h, z), rec)
    }
    let ss = SteinSteinParams::default();
    let m = |v: f32| c1_altitude_norm_to_metres(v, &ss);
    let viz_origin: [f64; 2] = [0.09375, 0.578_125];
    let states: Vec<(String, Knobs)> = std::env::var("F155_LOC_STATES")
        .unwrap_or_else(|_| "livré".to_string())
        .split(',')
        .map(|s| {
            let k = F121_AGE_K;
            let kn = match s.trim() {
                "livré" => Knobs::passes(2),
                "A1+B2" => Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
                "C1 nue" => Knobs { valley: Some(ValleyConstruction::new(k, None)), ..Knobs::passes(2) },
                "C2 /10 col" => Knobs { valley: Some(ValleyConstruction::new(k, Some(0.1))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
                "C2 /3" => Knobs { valley: Some(ValleyConstruction::new(k, Some(1.0 / 3.0))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
                "C2 /10 niveau mer" => Knobs { valley: Some(ValleyConstruction::f121(k, Some(0.1))), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) },
                o => panic!("unknown state {o}"),
            };
            (s.trim().to_string(), Knobs { origin: Some(viz_origin), ..kn })
        })
        .collect();
    eprintln!("\n==========  Finding 155 Br . the stop rule localised  ==========");
    for (label, kn) in states {
        eprintln!("\n────────── {label} ──────────");
        let wd = build_world(kn, None, PSEED, None);
        let g = &wd.heightmap;
        let (w, h) = (g.width, g.height);
        let n = w * h;
        let d = c1_drainage_windowed(g, None, &common::viz_dcfg(), &ss, DOMAIN_KM);
        let prot = wd.volc.enabled.then(|| crater_protect_mask(&wd.craters, w, h));
        let (zs, rec) = breach_rec(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), true);
        let (zf, recf) = breach_rec(g, &d.flow.filled, &d.lake_map, 0.5, w, h, prot.as_deref(), false);
        drop(d);
        // the below-sea land cells under each anchor, attributed
        for (name, z, r) in [("FLOOR (production)", &zf, &recf), ("SPILL (amendment 2)", &zs, &rec)] {
            let below: Vec<usize> = (0..n).filter(|&k| g.data[k] > SEA && z.data[k] <= SEA).collect();
            let by_ramp = below.iter().filter(|&&k| r[k].is_some()).count();
            let mut anchors: Vec<f32> = below.iter().filter_map(|&k| r[k].map(|e| m(e.1))).collect();
            anchors.sort_by(f32::total_cmp);
            let mut steps: Vec<u32> = below.iter().filter_map(|&k| r[k].map(|e| e.2)).collect();
            steps.sort_unstable();
            let pc = |v: &[f32], p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p) as usize] };
            let ps = |v: &[u32], p: f32| if v.is_empty() { 0 } else { v[((v.len() - 1) as f32 * p) as usize] };
            let mut tiles: BTreeMap<(usize, usize), usize> = BTreeMap::new();
            for &k in &below {
                *tiles.entry(((k % w) / 1024, (k / w) / 1024)).or_insert(0) += 1;
            }
            let mut tv: Vec<_> = tiles.into_iter().collect();
            tv.sort_by_key(|e| std::cmp::Reverse(e.1));
            eprintln!(
                "   {name} · below-sea land cells **{}** · by a ramp {by_ramp} · the ramp's anchor above the sea (m) p0 {:.2} · p50 {:.2} · p90 {:.2} · steps p50 {} · p90 {} · original height (m) p50 {:.2} · the 1024² tiles holding most {:?}",
                below.len(),
                pc(&anchors, 0.0),
                pc(&anchors, 0.5),
                pc(&anchors, 0.9),
                ps(&steps, 0.5),
                ps(&steps, 0.9),
                {
                    let mut o: Vec<f32> = below.iter().map(|&k| m(g.data[k])).collect();
                    o.sort_by(f32::total_cmp);
                    pc(&o, 0.5)
                },
                tv.iter().take(5).collect::<Vec<_>>()
            );
        }
        drop((zs, zf, rec, recf));
        // the lakes, line by line
        let vp = viz_hd_lakes_on(&wd, kn, PSEED, 45.0, 40.0);
        let vs = viz_hd_lakes_with(&wd, kn, PSEED, 45.0, 40.0, &|a, b, c, d2, e, f, p| breach_rec(a, b, c, d2, e, f, p, true).0);
        let (lp, ls) = (lake_listing(&vp.drainage.lake_map, &vp.drainage.lakes), lake_listing(&vs.drainage.lake_map, &vs.drainage.lakes));
        let diff: Vec<(usize, &String, &String)> = lp.iter().zip(ls.iter()).enumerate().filter(|(_, (a, b))| a != b).map(|(i, (a, b))| (i, a, b)).collect();
        eprintln!("   LAKES · listing lines {} / {} · lines that differ **{}**:", lp.len(), ls.len(), diff.len());
        for (i, a, b) in diff.iter().take(20) {
            eprintln!("      #{i}\n        production {a}\n        spill      {b}");
        }
        // per lake: the cells that change membership
        let mut moved: HashMap<u32, (usize, usize)> = HashMap::new(); // id → (cells lost, cells gained)
        for k in 0..n {
            let (a, b) = (vp.drainage.lake_map[k], vs.drainage.lake_map[k]);
            if a != b {
                if a != 0 {
                    moved.entry(a).or_insert((0, 0)).0 += 1;
                }
                if b != 0 {
                    moved.entry(b).or_insert((0, 0)).1 += 1;
                }
            }
        }
        let mut mv: Vec<_> = moved.into_iter().collect();
        mv.sort_by_key(|e| std::cmp::Reverse(e.1.0 + e.1.1));
        eprintln!("   LAKES · ids whose mask changes (lost / gained cells): {:?}", mv.iter().take(15).collect::<Vec<_>>());
    }
}
