//! ADR Finding 151 — geology v1, round 2 of 3: the fossil-belt prototype (F: past convergent / collision cells and the
//! sutures of the plate merges, recorded during the C1 loop through a read-only observer) and the measurements the
//! specification needs (S: the land shares of the rock classes' sources; the export sizes of 12 proxy chance grids).
//! Declared in `docs/reports/geology_v1/f151_spec/f151_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f151_geology -- --ignored --nocapture

mod common;

use common::{CANONICAL_ORIGIN, CELL_KM, DOMAIN_KM, Knobs, PSEED, SEA, build_field_seed, build_world, pct, sorted, viz_hd_lakes_on};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::lithology::{LithologyConfig, build_coarse_k, stamp_volcanic_k, upscale_k_to_hd};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::place_edifices;
use ymir_core::tectonics_c1::debug_labels::derive_tectonic_labels;
use ymir_core::tectonics_c1::drainage::LakeType;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures, run_with_closures_observed};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction, carve, skeleton};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const S_EQ: f32 = 0.024;

fn temoin_vc() -> ValleyConstruction {
    ValleyConstruction::new(F121_AGE_K, Some(0.1))
}

fn temoin() -> Knobs {
    Knobs { valley: Some(temoin_vc()), slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }
}

fn run_cfg() -> C1TimeLoopConfig {
    C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    }
}

/// The final coarse state's bits, for the read-only check.
fn state_bits(s: &C1State) -> Vec<u64> {
    let (nx, ny) = (s.nx(), s.ny());
    let mut v = Vec::with_capacity(nx * ny * 5);
    for j in 0..ny {
        for i in 0..nx {
            v.push(s.s.get(i, j).to_bits());
            v.push(s.age.get(i, j).to_bits());
            v.push(s.plate_id.get(i, j) as u64);
            v.push((s.plate_type.get(i, j) == PlateType::Continental) as u64);
            v.push(s.cratonic_mask.get(i, j) as u64);
        }
    }
    v
}

/// A 1024² north-up PNG (internal y = 0 is south), the sea greyed.
fn save_map(path: &std::path::Path, land: &dyn Fn(usize, usize) -> bool, col: &dyn Fn(usize, usize) -> [u8; 3]) {
    let n = 1024u32;
    let img = image::RgbImage::from_fn(n, n, |px, py| {
        let x = px as usize * 8 + 4;
        let y = (n - 1 - py) as usize * 8 + 4;
        let c = col(x, y);
        if land(x, y) { image::Rgb(c) } else { image::Rgb([(c[0] as u16 / 4 + 60) as u8, (c[1] as u16 / 4 + 62) as u8, (c[2] as u16 / 4 + 70) as u8]) }
    });
    img.save(path).expect("write the map");
}

/// Two-pass chamfer distance (cells) to the `true` cells of `m`.
fn chamfer(m: &[bool], w: usize, h: usize) -> Vec<f32> {
    let mut d: Vec<f32> = m.iter().map(|&b| if b { 0.0 } else { f32::INFINITY }).collect();
    let r2 = std::f32::consts::SQRT_2;
    for y in 0..h {
        for x in 0..w {
            let c = y * w + x;
            let mut v = d[c];
            if x > 0 { v = v.min(d[c - 1] + 1.0); }
            if y > 0 {
                v = v.min(d[c - w] + 1.0);
                if x > 0 { v = v.min(d[c - w - 1] + r2); }
                if x + 1 < w { v = v.min(d[c - w + 1] + r2); }
            }
            d[c] = v;
        }
    }
    for y in (0..h).rev() {
        for x in (0..w).rev() {
            let c = y * w + x;
            let mut v = d[c];
            if x + 1 < w { v = v.min(d[c + 1] + 1.0); }
            if y + 1 < h {
                v = v.min(d[c + w] + 1.0);
                if x + 1 < w { v = v.min(d[c + w + 1] + r2); }
                if x > 0 { v = v.min(d[c + w - 1] + r2); }
            }
            d[c] = v;
        }
    }
    d
}

/// The PNG (8-bit grey, deflate) size of a u8 grid.
fn png_bytes(g: &[u8], w: usize, h: usize) -> usize {
    use image::ImageEncoder;
    let mut out = Vec::new();
    image::codecs::png::PngEncoder::new(&mut out).write_image(g, w as u32, h as u32, image::ExtendedColorType::L8).expect("png");
    out.len()
}

/// Down-sample a u8 grid by `f` (the max of each f×f block: a chance zone never vanishes).
fn down(g: &[u8], w: usize, h: usize, f: usize) -> (Vec<u8>, usize, usize) {
    let (w2, h2) = (w / f, h / f);
    let mut o = vec![0u8; w2 * h2];
    for y in 0..h2 {
        for x in 0..w2 {
            let mut m = 0u8;
            for dy in 0..f {
                for dx in 0..f {
                    m = m.max(g[(y * f + dy) * w + x * f + dx]);
                }
            }
            o[y * w2 + x] = m;
        }
    }
    (o, w2, h2)
}

#[test]
#[ignore]
fn f151_geology() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/reports/geology_v1/f151_spec");
    std::fs::create_dir_all(&out).expect("the report folder");
    eprintln!("\n==========  Finding 151 . the fossil belts (F) and the spec's measurements (S)  ==========");
    // ── F: the plain run and the observed run, from the same init
    let tp = Instant::now();
    let mut s_plain = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut k_plain = PlateKinematics::preset_phase_1_1(s_plain.num_plates);
    run_with_closures(&mut s_plain, &mut k_plain, &run_cfg(), &C1Closures::default(), |_, _| {});
    let t_plain = tp.elapsed().as_secs_f64();
    let mut st = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let (nx, ny) = (st.nx(), st.ny());
    let nc = nx * ny;
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let mut conv = vec![0u32; nc];
    let mut coll = vec![0u32; nc];
    let mut upper = vec![0u32; nc];
    let mut div = vec![0u32; nc];
    let mut last_conv = vec![0u32; nc];
    let mut last_coll = vec![0u32; nc];
    let mut prev_id: Vec<u16> = (0..nc).map(|k| st.plate_id.get(k % nx, k / nx)).collect();
    let mut prev_type: Vec<bool> = (0..nc).map(|k| st.plate_type.get(k % nx, k / nx) == PlateType::Continental).collect();
    // (step, vanished plate, absorbing plate, suture cells, continental share)
    let mut sutures: Vec<(usize, u16, u16, Vec<usize>, f32)> = Vec::new();
    let mut splits = 0usize;
    let to = Instant::now();
    run_with_closures_observed(&mut st, &mut kin, &run_cfg(), &C1Closures::default(), |step, s, k| {
        let lab = derive_tectonic_labels(s, k);
        for c in 0..nc {
            let isconv = lab.subduction_upper[c] || lab.subduction_slab[c] || lab.collision[c];
            if isconv {
                conv[c] += 1;
                last_conv[c] = step as u32;
            }
            if lab.collision[c] {
                coll[c] += 1;
                last_coll[c] = step as u32;
            }
            if lab.subduction_upper[c] {
                upper[c] += 1;
            }
            if lab.divergent[c] {
                div[c] += 1;
            }
        }
        let now: Vec<u16> = (0..nc).map(|c| s.plate_id.get(c % nx, c / nx)).collect();
        let ids_prev: HashSet<u16> = prev_id.iter().copied().collect();
        let ids_now: HashSet<u16> = now.iter().copied().collect();
        splits += ids_now.difference(&ids_prev).count();
        for &a in ids_prev.difference(&ids_now) {
            // the absorbing plate: what the vanished plate's cells became (the majority)
            let mut votes: HashMap<u16, usize> = HashMap::new();
            for c in 0..nc {
                if prev_id[c] == a {
                    *votes.entry(now[c]).or_insert(0) += 1;
                }
            }
            let Some((&b, _)) = votes.iter().max_by_key(|e| (e.1, std::cmp::Reverse(*e.0))) else { continue };
            // the suture: the previous step's cells of a and b that touch each other (4-neighbours, periodic)
            let mut cells = Vec::new();
            for c in 0..nc {
                if prev_id[c] != a && prev_id[c] != b {
                    continue;
                }
                let other = if prev_id[c] == a { b } else { a };
                let (i, j) = ((c % nx) as i64, (c / nx) as i64);
                let touches = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)].iter().any(|(di, dj)| {
                    let q = ((j + dj).rem_euclid(ny as i64) as usize) * nx + (i + di).rem_euclid(nx as i64) as usize;
                    prev_id[q] == other
                });
                if touches {
                    cells.push(c);
                }
            }
            let cont = cells.iter().filter(|&&c| prev_type[c]).count() as f32 / cells.len().max(1) as f32;
            sutures.push((step, a, b, cells, cont));
        }
        prev_id = now;
        prev_type = (0..nc).map(|c| s.plate_type.get(c % nx, c / nx) == PlateType::Continental).collect();
    });
    let t_obs = to.elapsed().as_secs_f64();
    let identical = state_bits(&st) == state_bits(&s_plain);
    eprintln!("\n   F · read-only: the observed run's final state is bit-identical to the plain run's: **{identical}** · plain {t_plain:.2} s · observed and recorded {t_obs:.2} s (+{:.2} s)", t_obs - t_plain);
    let lab_end = derive_tectonic_labels(&st, &kin);
    let cont_end: Vec<bool> = (0..nc).map(|c| st.plate_type.get(c % nx, c / nx) == PlateType::Continental).collect();
    let craton: Vec<bool> = (0..nc).map(|c| st.cratonic_mask.get(c % nx, c / nx)).collect();
    let now_conv: Vec<bool> = (0..nc).map(|c| lab_end.subduction_upper[c] || lab_end.subduction_slab[c] || lab_end.collision[c]).collect();
    let cnt = |f: &dyn Fn(usize) -> bool| (0..nc).filter(|&c| f(c)).count();
    eprintln!(
        "   F · over 300 steps: ever convergent {} cells (continental at the end {}) · ever collision (C-C) {} · ever subduction upper plate {} · ever divergent {} · convergent at the end {} · collision at the end {}",
        cnt(&|c| conv[c] > 0),
        cnt(&|c| conv[c] > 0 && cont_end[c]),
        cnt(&|c| coll[c] > 0),
        cnt(&|c| upper[c] > 0),
        cnt(&|c| div[c] > 0),
        cnt(&|c| now_conv[c]),
        cnt(&|c| lab_end.collision[c])
    );
    for thr in [1u32, 5, 10, 20, 50] {
        eprintln!(
            "      collision ≥ {thr:>2} steps: {} cells ({} continental at the end, {} on the craton, {} convergent at the end) · any convergence ≥ {thr} steps: {} cells ({} continental)",
            cnt(&|c| coll[c] >= thr),
            cnt(&|c| coll[c] >= thr && cont_end[c]),
            cnt(&|c| coll[c] >= thr && craton[c]),
            cnt(&|c| coll[c] >= thr && now_conv[c]),
            cnt(&|c| conv[c] >= thr),
            cnt(&|c| conv[c] >= thr && cont_end[c])
        );
    }
    let last: Vec<f32> = (0..nc).filter(|&c| coll[c] > 0).map(|c| last_coll[c] as f32).collect();
    eprintln!("      the last collision step of the ever-collision cells: p10 / p50 / p90 {:.0} / {:.0} / {:.0} (of 300; ≈ 0.67 Ma per step, PROXY)", pct(&sorted(last.clone()), 0.1), pct(&sorted(last.clone()), 0.5), pct(&sorted(last), 0.9));
    eprintln!("   F · plate merges (sutures) {} · rift splits (new plate ids) {splits}", sutures.len());
    let mut suture_cells: HashSet<usize> = HashSet::new();
    for (step, a, b, cells, cont) in &sutures {
        eprintln!("      step {step:>3}: plate {a} into {b} · suture {} cells · continental {:.0} %", cells.len(), 100.0 * cont);
        suture_cells.extend(cells.iter().copied());
    }
    // the proposed fossil belt: collision ≥ 10 steps, or a suture cell on continental crust; not convergent at the end
    let fossil: Vec<bool> = (0..nc).map(|c| !now_conv[c] && (coll[c] >= 10 || (suture_cells.contains(&c) && cont_end[c]))).collect();
    eprintln!(
        "   F · proposed fossil belt (collision ≥ 10 steps, or a continental suture cell; not convergent at the end): **{} cells** ({} continental, {} on the craton) · suture cells {} ({} continental at the end)",
        cnt(&|c| fossil[c]),
        cnt(&|c| fossil[c] && cont_end[c]),
        cnt(&|c| fossil[c] && craton[c]),
        suture_cells.len(),
        suture_cells.iter().filter(|&&c| cont_end[c]).count()
    );
    // ── the témoin's HD world (the maps' land mask, S)
    let wd = build_world(temoin(), None, PSEED, None);
    let same_state = state_bits(&wd.state) == state_bits(&st);
    eprintln!("   F · the témoin's coarse state equals this run's: {same_state}");
    let f = &wd.heightmap;
    let (w, h) = (f.width, f.height);
    let n = w * h;
    let land = |x: usize, y: usize| f.data[y * w + x] > SEA;
    let (ox, oy) = (CANONICAL_ORIGIN[0] * nx as f64, CANONICAL_ORIGIN[1] * ny as f64);
    let coarse = |x: usize, y: usize| -> usize {
        let sx = ox + x as f64 * nx as f64 / w as f64;
        let sy = oy + y as f64 * ny as f64 / h as f64;
        ((sy.round() as i64).rem_euclid(ny as i64) as usize) * nx + (sx.round() as i64).rem_euclid(nx as i64) as usize
    };
    let mx = conv.iter().copied().max().unwrap_or(1).max(1) as f32;
    save_map(&out.join("map_f_ever_convergent_steps.png"), &land, &|x, y| {
        let c = coarse(x, y);
        if conv[c] == 0 { [110, 110, 110] } else { let t = (conv[c] as f32 / mx).sqrt(); [(80.0 + 170.0 * t) as u8, (60.0 + 40.0 * t) as u8, (160.0 - 120.0 * t) as u8] }
    });
    save_map(&out.join("map_f_ever_collision_and_sutures.png"), &land, &|x, y| {
        let c = coarse(x, y);
        if suture_cells.contains(&c) { [250, 230, 60] } else if coll[c] >= 10 { [230, 60, 200] } else if coll[c] > 0 { [150, 90, 150] } else if craton[c] { [205, 178, 100] } else { [110, 110, 110] }
    });
    save_map(&out.join("map_f_proposed_fossil_belt.png"), &land, &|x, y| {
        let c = coarse(x, y);
        if now_conv[c] { [245, 140, 30] } else if fossil[c] { [200, 60, 160] } else if craton[c] { [205, 178, 100] } else { [130, 130, 130] }
    });
    // ── S: the land shares of the rock classes' sources
    let v = viz_hd_lakes_on(&wd, temoin(), PSEED, 45.0, 40.0);
    let dr = &v.drainage;
    let lcfg = LithologyConfig { enabled: true, soft_multiplier: 10.0, volcanic_multiplier: 3.0, rift_age_threshold: 1.0 };
    let k_rift = upscale_k_to_hd(&build_coarse_k(&wd.state, &lcfg), w, h, CANONICAL_ORIGIN, 1.0);
    let edifices = place_edifices(&wd.state, &wd.kin, &WorldSeed::new(PSEED), DOMAIN_KM, &wd.volc);
    let mut k_volc = vec![1.0f32; n];
    stamp_volcanic_k(&mut k_volc, &edifices, CANONICAL_ORIGIN, 1.0, CELL_KM, w, h, &lcfg);
    let s1 = build_field_seed(Knobs { no_incision: true, erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let sk = skeleton(&s1, &temoin_vc(), &ss, DOMAIN_KM);
    let (_, mk) = carve(&s1, &sk, &temoin_vc(), &ss);
    drop(sk);
    drop(s1);
    let zm: Vec<f32> = v.conditioned.data.iter().map(|&x| c1_altitude_norm_to_metres(x, &ss)).collect();
    let cell_m = CELL_KM * 1000.0;
    let acc = &dr.flow.accumulation.data;
    let acc_1km2 = 1.0 / (CELL_KM * CELL_KM);
    let slope = |c: usize| -> f32 {
        let (x, y) = (c % w, c / w);
        let at = |xx: usize, yy: usize| zm[yy * w + xx];
        let gx = (at((x + 1).min(w - 1), y) - at(x.saturating_sub(1), y)) / (2.0 * cell_m);
        let gy = (at(x, (y + 1).min(h - 1)) - at(x, y.saturating_sub(1))) / (2.0 * cell_m);
        (gx * gx + gy * gy).sqrt()
    };
    let endo: HashSet<u32> = dr.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).map(|l| l.base.id).collect();
    let near_belt: Vec<bool> = {
        // the current belts: convergent at the end within 2 coarse cells
        let mut m = vec![false; nc];
        for c in 0..nc {
            if now_conv[c] {
                let (i, j) = ((c % nx) as i64, (c / nx) as i64);
                for dj in -2i64..=2 {
                    for di in -2i64..=2 {
                        m[((j + dj).rem_euclid(ny as i64) as usize) * nx + (i + di).rem_euclid(nx as i64) as usize] = true;
                    }
                }
            }
        }
        m
    };
    let classes = ["evaporites (endorheic lake)", "loose deposits", "volcanic", "rift-soft", "belt (fossil or current)", "craton", "default basement"];
    let class_of = |c: usize| -> usize {
        let cc = coarse(c % w, c / w);
        let lid = dr.lake_map[c];
        if lid != 0 && endo.contains(&lid) {
            0
        } else if mk.floor[c] || lid != 0 || (slope(c) < 0.02 && acc[c] >= acc_1km2) {
            1
        } else if k_volc[c] > 1.5 {
            2
        } else if k_rift[c] > 1.5 {
            3
        } else if fossil[cc] || near_belt[cc] {
            4
        } else if craton[cc] {
            5
        } else {
            6
        }
    };
    let mut counts = [0usize; 7];
    let mut cls = vec![255u8; n];
    let mut land_n = 0usize;
    for c in 0..n {
        if f.data[c] <= SEA {
            continue;
        }
        land_n += 1;
        let k = class_of(c);
        counts[k] += 1;
        cls[c] = k as u8;
    }
    eprintln!("\n   S · the land shares of the rock classes' sources (priority in this order; this measurement's thresholds, not the spec's) over {land_n} land cells:");
    for (k, name) in classes.iter().enumerate() {
        eprintln!("      {name}: {:.2} %", 100.0 * counts[k] as f64 / land_n as f64);
    }
    let pal = [[230, 220, 180], [200, 170, 120], [175, 70, 210], [40, 185, 175], [200, 60, 160], [205, 178, 100], [140, 140, 140]];
    save_map(&out.join("map_s_class_sources.png"), &land, &|x, y| {
        let k = cls[y * w + x];
        if k == 255 { [110, 110, 110] } else { pal[k as usize] }
    });
    // ── the export sizes: 12 PROXY chance grids (a field mask × 100 − 10 per km outside)
    let landv: Vec<bool> = f.data.iter().map(|&x| x > SEA).collect();
    let rivers: Vec<bool> = (0..n).map(|c| acc[c] >= 20.0 * acc_1km2).collect();
    let mut masks: BTreeMap<&str, Vec<bool>> = BTreeMap::new();
    masks.insert("01 iron (craton)", (0..n).map(|c| cls[c] == 5).collect());
    masks.insert("02 gold (belt + volcanic)", (0..n).map(|c| cls[c] == 4 || cls[c] == 2).collect());
    masks.insert("03 silver (belt)", (0..n).map(|c| cls[c] == 4).collect());
    masks.insert("04 copper (volcanic + rift)", (0..n).map(|c| cls[c] == 2 || cls[c] == 3).collect());
    masks.insert("05 tin (belt)", (0..n).map(|c| cls[c] == 4).collect());
    masks.insert("06 stone (craton + volcanic + default)", (0..n).map(|c| cls[c] == 5 || cls[c] == 2 || cls[c] == 6).collect());
    masks.insert("07 clay (loose deposits)", (0..n).map(|c| cls[c] == 1).collect());
    masks.insert("08 sand and gravel (rivers)", rivers.clone());
    masks.insert("09 salt (evaporites)", (0..n).map(|c| cls[c] == 0).collect());
    masks.insert("10 peat (loose, flat)", (0..n).map(|c| cls[c] == 1 && slope(c) < 0.01).collect());
    masks.insert("11 sulphur and obsidian (volcanic)", (0..n).map(|c| cls[c] == 2).collect());
    masks.insert("12 gems (craton + belt)", (0..n).map(|c| cls[c] == 5 || cls[c] == 4).collect());
    let per_km = 1.0 / CELL_KM;
    let (mut raw, mut p1, mut p2, mut p4) = (0usize, 0usize, 0usize, 0usize);
    let tz = Instant::now();
    for (name, m) in &masks {
        let d = chamfer(m, w, h);
        let g: Vec<u8> = (0..n).map(|c| if !landv[c] { 0 } else { (100.0 - 10.0 * d[c] / per_km).clamp(0.0, 100.0) as u8 }).collect();
        let nonzero = g.iter().filter(|&&x| x > 0).count();
        let b1 = png_bytes(&g, w, h);
        let (g2, w2, h2) = down(&g, w, h, 2);
        let b2 = png_bytes(&g2, w2, h2);
        let (g4, w4, h4) = down(&g, w, h, 4);
        let b4 = png_bytes(&g4, w4, h4);
        raw += n;
        p1 += b1;
        p2 += b2;
        p4 += b4;
        eprintln!("      {name}: non-zero {:.1} % of the grid · PNG full {:.2} MB · ½ {:.2} MB · ¼ {:.2} MB", 100.0 * nonzero as f64 / n as f64, b1 as f64 / 1e6, b2 as f64 / 1e6, b4 as f64 / 1e6);
    }
    eprintln!(
        "   S · the 12 PROXY chance grids: raw full {:.0} MB · ½ {:.0} MB · ¼ {:.0} MB · **PNG full {:.1} MB · ½ {:.1} MB · ¼ {:.1} MB** ({:.1} s for the 12 at full with distances)",
        raw as f64 / 1e6,
        raw as f64 / 4e6,
        raw as f64 / 16e6,
        p1 as f64 / 1e6,
        p2 as f64 / 1e6,
        p4 as f64 / 1e6,
        tz.elapsed().as_secs_f64()
    );
    let cls_png = png_bytes(&cls.iter().map(|&k| if k == 255 { 0 } else { k + 1 }).collect::<Vec<u8>>(), w, h);
    eprintln!("   S · the rock-class grid (u8, full resolution): raw {:.0} MB · PNG {:.2} MB", n as f64 / 1e6, cls_png as f64 / 1e6);
    eprintln!("\n==========  end Finding 151 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
