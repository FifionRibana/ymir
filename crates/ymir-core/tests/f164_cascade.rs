//! ADR Finding 164 — the cascade, round 6: bend the trunks inherited from the coarse grid with a smooth, invertible warp
//! of the upscaled field (W-phys once, W-tous at every level; two settings), arbitrate the physics level's MFD exponent
//! (p = 2, 4, 6, each with F163's roughness), and attribute the flat interfluves (by altitude, distance to the coast,
//! nearness to the D5 lakes). All chains R4+π to 2 048², judged against Corsica at 195 m.
//! Declared in `docs/reports/relief_method/f164_warp/f164_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f164_cascade -- --ignored --exact f164_cascade --nocapture

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, Knobs, PSEED, build_world, pct, sorted};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, ocean_mask};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, flat_map, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::{CascadeConfig, physics_level, roll};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const PEAK_TARGET_M: f32 = 2850.0;
const LEVELS: [usize; 3] = [512, 1024, 2048];
const CROPS: [(&str, usize, usize); 3] = [("A_chaine", 128, 320), ("B_plateau", 176, 192), ("C_bassin_lac4", 74, 160)];
const TIGHT: (usize, usize) = (176, 352);
const CORSE_SPINE: (f32, f32) = (72.0 / 256.0, 96.0 / 256.0);
const CORSE_FLANK: (f32, f32) = (100.0 / 256.0, 140.0 / 256.0);
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];
/// F164-D1: the two declared settings (A, L in previous-level cells).
const SETTINGS: [(&str, f32, f32); 2] = [("W1", 0.5, 4.0), ("W2", 0.5, 2.0)];

fn cache() -> PathBuf {
    root().join("target_f164/f164_cache")
}
fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f164_warp/images")
}
fn out15(a: f32, b: f32) -> bool {
    !(a.is_finite() && b.is_finite() && b != 0.0) || a > 1.5 * b || a < b / 1.5
}
fn lnr(a: f32, b: f32) -> f32 {
    if a > 0.0 && b > 0.0 && a.is_finite() && b.is_finite() { (a / b).ln().abs() } else { 0.0 }
}
fn fmt_pf(p: &Planform) -> String {
    format!(
        "A_dir {:.1} % (trunks {:.1} %) · S−1 {:.3} / {:.3} / {:.3} · F {:.2} %",
        100.0 * p.a_dir,
        100.0 * p.a_dir_trunks,
        p.excess(0),
        p.excess(1),
        p.excess(2),
        100.0 * p.flat
    )
}

/// F164-D2: one warp's controls.
#[derive(Clone, Debug)]
struct WarpCtl {
    lakes: (usize, usize),
    hyps: f32,
    land: (usize, usize),
    jmin: f32,
    jlow: f32,
}

fn land_deciles(g: &GridF32) -> (Vec<f32>, usize) {
    let oc = ocean_mask(g);
    let v = sorted(g.data.iter().zip(&oc).filter(|(_, o)| !**o).map(|(v, _)| *v).collect());
    ((1..10).map(|d| pct(&v, d as f64 / 10.0)).collect(), v.len())
}
fn warp_ctl(a: &GridF32, b: &GridF32, cell_km: f32, jmin: f32, jlow: f32) -> WarpCtl {
    let (qa, la) = land_deciles(a);
    let (qb, lb) = land_deciles(b);
    let hyps = qa.iter().zip(&qb).map(|(x, y)| (x - y).abs()).fold(0f32, f32::max) / qa[8].max(1.0);
    WarpCtl { lakes: (hydrology(a, cell_km).lakes.len(), hydrology(b, cell_km).lakes.len()), hyps, land: (la, lb), jmin, jlow }
}

#[allow(dead_code)]
struct Lvl {
    n: usize,
    result: GridF32,
    d5: Vec<bool>,
    r: Reading,
    p: Option<Predict>,
    pf: Option<Planform>,
    k: f32,
    drift: f32,
    trunk_p90: f32,
    pi_var: f64,
    warp: Option<WarpCtl>,
}

fn run_chain(name: &str, mut ch: Chain, phys_r: &Reading, corse: &BTreeMap<usize, Reading>) -> (Vec<Lvl>, f64, Vec<String>) {
    let t0 = Instant::now();
    let mut out: Vec<Lvl> = Vec::new();
    let mut prev_r = phys_r.clone();
    let mut prev = ch.physics_m.clone();
    let mut stops = Vec::new();
    for n in LEVELS {
        let target = corse[&(n / 2)].octaves[0];
        let (k, trials) = ch.calibrate_peak(target, &|| false);
        let l = ch.next_level(&|| false).unwrap().clone();
        let result = l.result();
        let cell_km = DOMAIN_KM / n as f32;
        let (r, _) = read(&result, Some(&l.d5));
        let p = (n >= 1024).then(|| predictability(&result, cell_km, Some(&l.d5)));
        let pf = (n >= 1024).then(|| planform(&result, cell_km, Some(&l.d5)));
        let warp = l.pre_warp.as_ref().zip(l.post_warp.as_ref()).map(|(a, b)| warp_ctl(a, b, cell_km, l.warp_jmin, l.warp_jlow));
        let (b, _) = drift(&result, &prev);
        let rel = b / prev_r.mean.max(1.0);
        let tp = trunk_dist(&r.trunks, &prev_r.trunks, 2.0, n / 2).1;
        let pi_var = l.pi.iter().filter(|v| **v != 0.0).map(|v| (*v as f64) * (*v as f64)).sum::<f64>()
            / l.pi.iter().filter(|v| **v != 0.0).count().max(1) as f64;
        let mut why = Vec::new();
        if rel.abs() > 0.02 {
            why.push(format!("drift {:+.2} %", 100.0 * rel));
        }
        if r.peak < PEAK_OK.0 || r.peak > PEAK_OK.1 {
            why.push(format!("peak {:.0} m", r.peak));
        }
        if n < 2048 && tp > 1.5 {
            why.push(format!("trunk p90 {tp:.2}"));
        }
        if !why.is_empty() {
            stops.push(format!("{n}²: {}", why.join("; ")));
        }
        eprintln!(
            "      {name} {n}² ({:.1} s; k {k:.3e} after {} trials) · peak {:.0} · drift {:+.2} % · trunk p90 {:.2} · octave 0 {:.2} m (Corse {:.2}){}{}",
            l.secs.iter().sum::<f64>(),
            trials.len(),
            r.peak,
            100.0 * rel,
            tp,
            r.octaves[0],
            target,
            warp.as_ref().map_or(String::new(), |w| format!(
                " · warp: lakes {} → {} · hypsometry Δ {:.3} % of q9 · land {} → {} ({:+.3} %) · min J {:.3} · J < 0.5 {:.3} %",
                w.lakes.0, w.lakes.1, 100.0 * w.hyps, w.land.0, w.land.1, 100.0 * (w.land.1 as f64 / w.land.0.max(1) as f64 - 1.0), w.jmin, 100.0 * w.jlow
            )),
            if why.is_empty() { String::new() } else { format!(" · per-level R: {}", why.join("; ")) }
        );
        out.push(Lvl { n, result: result.clone(), d5: l.d5.clone(), r: r.clone(), p, pf, k: l.k, drift: rel, trunk_p90: tp, pi_var, warp });
        prev_r = r;
        prev = result;
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    (out, t0.elapsed().as_secs_f64(), stops)
}

/// R at 2 048² (F164-R): the reasons, and Σ |ln ratio| over R's measures.
fn judge(lv: &[Lvl], co: &Reading, cp: &Predict, cpf: &Planform, ctrl_ok: bool) -> (Vec<String>, f32) {
    let l = &lv[2];
    let r = &l.r;
    let mut why = Vec::new();
    let mut s = 0f32;
    let mut m = |name: String, a: f32, b: f32, why: &mut Vec<String>| {
        s += lnr(a, b);
        if out15(a, b) {
            why.push(format!("{name} {a:.3} / {b:.3}"));
        }
    };
    if ctrl_ok && r.facets > 1.5 * co.facets && r.facet_cells > 20 {
        why.push(format!("facets {:.2} > 1.5 × {:.2}", r.facets, co.facets));
    }
    if r.walls28 > 1.5 * co.walls28 && r.walls28_cells > 20 {
        why.push(format!("walls 28° {:.3} > 1.5 × {:.3}", r.walls28, co.walls28));
    }
    m("slope p50".into(), r.slope_p50, co.slope_p50, &mut why);
    m("slope p90".into(), r.slope_p90, co.slope_p90, &mut why);
    m("λ km".into(), r.lambda_km, co.lambda_km, &mut why);
    let cell_km = DOMAIN_KM / l.n as f32;
    for j in 0..r.octaves.len().min(co.octaves.len()) {
        let (a, b) = octave_km(j, cell_km);
        if a >= 0.39 && b <= 12.6 {
            m(format!("octave {a:.2}–{b:.2} km"), r.octaves[j], co.octaves[j], &mut why);
        }
    }
    if let Some(p) = &l.p {
        for (i, ((name, x, _), (_, c, _))) in p.values().iter().zip(cp.values().iter()).enumerate() {
            if VALID_P[i] {
                m(name.to_string(), *x, *c, &mut why);
            }
        }
    }
    let band: f64 = (0..r.octaves.len())
        .filter(|&j| {
            let (a, b) = octave_km(j, cell_km);
            a >= 0.39 && b <= 12.6
        })
        .map(|j| (r.octaves[j] as f64).powi(2))
        .sum();
    let share = lv.iter().map(|x| x.pi_var).sum::<f64>() / band.max(1e-9);
    if share > 0.10 {
        why.push(format!("π share {:.1} %", 100.0 * share));
    }
    if let Some(pf) = &l.pf {
        m("A_dir".into(), pf.a_dir, cpf.a_dir, &mut why);
        for i in 0..3 {
            m(format!("S_{} − 1", WINDOWS_KM[i]), pf.excess(i), cpf.excess(i), &mut why);
        }
        m("F".into(), pf.flat, cpf.flat, &mut why);
    }
    for x in lv {
        if let Some(w) = &x.warp {
            if w.lakes.1 as f64 > 1.01 * w.lakes.0 as f64 {
                why.push(format!("warp {}²: depressions {} → {}", x.n, w.lakes.0, w.lakes.1));
            }
            if w.jmin <= 0.0 {
                why.push(format!("warp {}²: min J {:.3}", x.n, w.jmin));
            }
        }
    }
    (why, s)
}

/// F164-F: three numbers per band: its share of the land read, its share of the flat cells, its own flat rate.
fn bands(read: &[bool], flat: &[bool], key: &[f32], edges: &[f32]) -> Vec<(f32, f32, f32)> {
    let total = read.iter().filter(|&&r| r).count().max(1) as f32;
    let total_flat = read.iter().zip(flat).filter(|(r, f)| **r && **f).count().max(1) as f32;
    (0..edges.len() - 1)
        .map(|b| {
            let inb: Vec<usize> = (0..read.len()).filter(|&k| read[k] && key[k] >= edges[b] && key[k] < edges[b + 1]).collect();
            let f = inb.iter().filter(|&&k| flat[k]).count() as f32;
            (inb.len() as f32 / total, f / total_flat, f / inb.len().max(1) as f32)
        })
        .collect()
}
/// Chebyshev distance (cells) from the `src` cells, capped at `cap`.
fn cheb(src: &[bool], n: usize, cap: u32) -> Vec<u32> {
    let mut d = vec![u32::MAX; n * n];
    let mut front: Vec<usize> = (0..n * n).filter(|&k| src[k]).collect();
    for &k in &front {
        d[k] = 0;
    }
    let mut lvl = 0;
    while !front.is_empty() && lvl < cap {
        lvl += 1;
        let mut next = Vec::new();
        for &k in &front {
            let (x, y) = ((k % n) as i64, (k / n) as i64);
            for dy in -1..=1i64 {
                for dx in -1..=1i64 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx < 0 || ny < 0 || nx >= n as i64 || ny >= n as i64 {
                        continue;
                    }
                    let j = ny as usize * n + nx as usize;
                    if d[j] == u32::MAX {
                        d[j] = lvl;
                        next.push(j);
                    }
                }
            }
        }
        front = next;
    }
    d
}
fn attribute(tag: &str, g: &GridF32, cell_km: f32, d5: Option<&[bool]>) -> f32 {
    let n = g.width;
    let (read, flat) = flat_map(g, cell_km, d5);
    let alt_edges = [f32::MIN, 50.0, 100.0, 300.0, f32::MAX];
    let alt = bands(&read, &flat, &g.data, &alt_edges);
    let ocean = ocean_mask(g);
    let dc = cheb(&ocean, n, (40.0 / cell_km) as u32 + 2);
    let dkm: Vec<f32> = dc.iter().map(|&d| if d == u32::MAX { 1e9 } else { d as f32 * cell_km }).collect();
    let coast = bands(&read, &flat, &dkm, &[0.0, 2.0, 5.0, 10.0, 20.0, f32::MAX]);
    let f_all = read.iter().zip(&flat).filter(|(r, f)| **r && **f).count() as f32 / read.iter().filter(|&&r| r).count().max(1) as f32;
    eprintln!("   ── {tag}: F {:.2} % of the {} cells read", 100.0 * f_all, read.iter().filter(|&&r| r).count());
    for (i, name) in ["0–50 m", "50–100 m", "100–300 m", "> 300 m"].iter().enumerate() {
        eprintln!("      altitude {name:<10} land {:5.1} % · flat cells {:5.1} % · flat rate {:5.2} %", 100.0 * alt[i].0, 100.0 * alt[i].1, 100.0 * alt[i].2);
    }
    for (i, name) in ["0–2 km", "2–5 km", "5–10 km", "10–20 km", "> 20 km"].iter().enumerate() {
        eprintln!("      coast {name:<10} land {:5.1} % · flat cells {:5.1} % · flat rate {:5.2} %", 100.0 * coast[i].0, 100.0 * coast[i].1, 100.0 * coast[i].2);
    }
    if let Some(m) = d5 {
        if m.iter().any(|&b| b) {
            let dd = cheb(m, n, (2.0 / cell_km).ceil() as u32 + 1);
            let near: Vec<f32> = dd.iter().map(|&d| if d as f32 * cell_km <= 2.0 { 0.0 } else { 1.0 }).collect();
            let b = bands(&read, &flat, &near, &[0.0, 0.5, 2.0]);
            eprintln!(
                "      within 2 km of a D5 cell: land {:.1} % · flat cells {:.1} % · flat rate {:.2} % (the rest {:.2} %)",
                100.0 * b[0].0,
                100.0 * b[0].1,
                100.0 * b[0].2,
                100.0 * b[1].2
            );
        }
    }
    let low = alt[0].1 + alt[1].1;
    let mostly_low = low >= 0.6 && alt[0].2.max(alt[1].2) > 0.0 && ((alt[0].2 * alt[0].0 + alt[1].2 * alt[1].0) / (alt[0].0 + alt[1].0).max(1e-9)) >= 2.0 * alt[3].2;
    eprintln!(
        "      → {:.1} % of the flat cells below 100 m · the rate below 100 m {:.2} % against {:.2} % above 300 m → {}",
        100.0 * low,
        100.0 * (alt[0].2 * alt[0].0 + alt[1].2 * alt[1].0) / (alt[0].0 + alt[1].0).max(1e-9),
        100.0 * alt[3].2,
        if mostly_low { "MOSTLY LOW (the declared rule)" } else { "SPREAD OUT (the declared rule)" }
    );
    low
}

struct Run {
    name: String,
    p: f32,
    place: &'static str,
    setting: &'static str,
    lv: Vec<Lvl>,
    secs: f64,
    stops: Vec<String>,
    why: Vec<String>,
    sum_ln: f32,
}

#[test]
#[ignore]
fn f164_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    std::fs::create_dir_all(cache()).unwrap();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 164 . the cascade, round 6: the warp, the exponent p, the flat interfluves  ==========");
    // ── Corsica ──
    let mut corse: BTreeMap<usize, Reading> = BTreeMap::new();
    let mut corse_g: BTreeMap<usize, GridF32> = BTreeMap::new();
    for nc in [128usize, 256, 512, 1024] {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).expect("data/corsica (prep_corse.py)");
        let (r, _) = read_cell(&g, None, 200.0 / nc as f32);
        corse.insert(nc, r);
        corse_g.insert(nc, g);
    }
    let cg195 = &corse_g[&1024];
    let cpf195 = planform(cg195, 200.0 / 1024.0, None);
    let cp195 = predictability(cg195, 200.0 / 1024.0, None);
    let cpf391 = planform(&corse_g[&512], 200.0 / 512.0, None);
    eprintln!("   Corse 195 m: {}", fmt_pf(&cpf195));
    eprintln!("   Corse 391 m: {}", fmt_pf(&cpf391));
    // ── the témoin's coarse field ──
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let cpath = cache().join("coarse64.bin");
    if !cpath.exists() {
        let t = Instant::now();
        let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
        let wd = build_world(temoin, None, PSEED, None);
        let coarse = c1_coarse_normalized_altitude(&wd.state, &IsostasyConfig::c1_default(), &ss, wd.cfg.target_land_fraction);
        save_grid(&cpath, &coarse);
        drop(wd);
        eprintln!("   the témoin's coarse field built ({:.0} s)", t.elapsed().as_secs_f64());
    }
    let coarse = roll(&load_grid(&cpath).unwrap(), off);
    // ── D3: the physics levels, p = 2, 4, 6, each with F163's roughness ──
    let rough_a = 0.5 * corse[&128].octaves[0];
    let base = CascadeConfig { seed: PSEED, roughness_m: rough_a, ..CascadeConfig::declared(ss.depth_scale_m as f32) };
    let mut phys: BTreeMap<u32, (GridF32, Reading, Vec<bool>)> = BTreeMap::new();
    eprintln!("\n   D3 — the physics levels (roughness {rough_a:.1} m, U₀ on {PEAK_TARGET_M:.0} m):");
    for p in [2.0f32, 4.0, 6.0] {
        let c = CascadeConfig { mfd_exponent: Some(p), ..base.clone() };
        let (rec, cal) = physics_level(&coarse, 256, &c, PEAK_TARGET_M, 300, &mut |_| {}, &|| false).unwrap();
        let pm = to_m(&rec.z, n2m);
        let (r, _) = read(&pm, Some(&cal.lifted_mask));
        eprintln!("   p = {p}: {:.1} s · {} steps · peak {:.0} m · λ {:.2} km · facets {:.2} % · octaves {:.1} / {:.1} m", cal.secs, rec.steps, r.peak, r.lambda_km, r.facets, r.octaves[0], r.octaves[1]);
        phys.insert(p as u32, (pm, r, cal.lifted_mask.clone()));
    }
    let pi_table: Vec<(usize, f32)> = LEVELS.iter().map(|&n| (n, 0.15 * corse[&(n / 2)].octaves[0])).collect();
    let mk_cfg = |warp: Option<Warp>| AmpConfig {
        talus_fine_scale: 0.25,
        n1_recal_depth: 2,
        n1_pi_m: pi_table.clone(),
        warp,
        ..AmpConfig::declared(PSEED, Variant::N1, 1.0)
    };
    let co195 = &corse[&1024];
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    let mut runs: Vec<Run> = Vec::new();
    let go = |p: f32, place: &'static str, setting: &'static str, runs: &mut Vec<Run>| {
        let warp = SETTINGS.iter().find(|s| s.0 == setting).and_then(|s| match place {
            "W-phys" => Some(Warp { amp: s.1, corr: s.2, place: WarpPlace::Phys }),
            "W-tous" => Some(Warp { amp: s.1, corr: s.2, place: WarpPlace::All }),
            _ => None,
        });
        let name = if place == "none" { format!("p{p} · none") } else { format!("p{p} · {place} {setting}") };
        eprintln!("\n   chain {name}:");
        let (pm, pr, d5) = &phys[&(p as u32)];
        let ch = Chain::new(pm.clone(), mk_cfg(warp)).with_d5(d5.clone());
        let (lv, secs, stops) = run_chain(&name, ch, pr, &corse);
        let (why, sum_ln) = judge(&lv, co195, &cp195, &cpf195, ctrl_ok);
        eprintln!("      → {} · Σ|ln| {:.2} · {:.0} s", if why.is_empty() { "PASSES R".to_string() } else { format!("REJECTED ({}): {}", why.len(), why.join("; ")) }, sum_ln, secs);
        runs.push(Run { name, p, place, setting, lv, secs, stops, why, sum_ln });
    };
    // stage 1: the settings at p = 2
    go(2.0, "none", "", &mut runs);
    for (s, _, _) in SETTINGS {
        go(2.0, "W-phys", s, &mut runs);
        go(2.0, "W-tous", s, &mut runs);
    }
    let set_score = |s: &str, runs: &[Run]| -> (usize, f32) {
        let r: Vec<&Run> = runs.iter().filter(|r| r.setting == s).collect();
        (r.iter().map(|r| r.why.len()).sum(), r.iter().map(|r| r.lv[2].pf.as_ref().unwrap().a_dir).sum::<f32>() / r.len().max(1) as f32)
    };
    let (s1, s2) = (set_score("W1", &runs), set_score("W2", &runs));
    let carried: &'static str = if (s1.0, s1.1) <= (s2.0, s2.1) { "W1" } else { "W2" };
    eprintln!("\n   the settings at p = 2: W1 {} reasons, mean A_dir {:.1} % · W2 {} reasons, mean A_dir {:.1} % → carried: {carried}", s1.0, 100.0 * s1.1, s2.0, 100.0 * s2.1);
    // stage 2: p = 4 and 6
    for p in [4.0f32, 6.0] {
        go(p, "none", "", &mut runs);
        go(p, "W-phys", carried, &mut runs);
        go(p, "W-tous", carried, &mut runs);
    }
    // ── the table at 2048² ──
    eprintln!("\n   at 2048² against Corsica 195 m (facet control pyramids {:.1} % · paraboloids {:.1} %):", ctrl.0, ctrl.1);
    eprintln!("   Corse 195 m: slope p50/p90 {:.3}/{:.3} · λ {:.2} · facets {:.2} · octaves {} · {}", co195.slope_p50, co195.slope_p90, co195.lambda_km, co195.facets,
        (0..6).map(|j| format!("{:.1}", co195.octaves[j])).collect::<Vec<_>>().join("/"), fmt_pf(&cpf195));
    for r in &runs {
        let l = &r.lv[2];
        eprintln!(
            "   {:<18} slope {:.3}/{:.3} · λ {:.2} · facets {:.2} · walls {:.3} · octaves {} · {} · R2_g {:.3} · {:.0} s · {}",
            r.name, l.r.slope_p50, l.r.slope_p90, l.r.lambda_km, l.r.facets, l.r.walls28,
            (0..6).map(|j| format!("{:.1}", l.r.octaves[j])).collect::<Vec<_>>().join("/"),
            fmt_pf(l.pf.as_ref().unwrap()),
            l.p.as_ref().map_or(f32::NAN, |p| p.r2_gully),
            r.secs,
            if r.why.is_empty() { "PASSES".to_string() } else { format!("{} reasons", r.why.len()) }
        );
        let l1 = &r.lv[1];
        eprintln!("      1024² (Corse 391 m: {}): {}", fmt_pf(&cpf391), fmt_pf(l1.pf.as_ref().unwrap()));
        if !r.stops.is_empty() {
            eprintln!("      per-level R: {}", r.stops.join(" | "));
        }
    }
    let bi = (0..runs.len()).min_by(|&a, &b| runs[a].why.len().cmp(&runs[b].why.len()).then(runs[a].sum_ln.partial_cmp(&runs[b].sum_ln).unwrap())).unwrap();
    let best = &runs[bi];
    eprintln!("\n   → the best chain: {} ({} reasons, Σ|ln| {:.2})", best.name, best.why.len(), best.sum_ln);
    // ── F: the flat interfluves attributed ──
    eprintln!("\n   F — the flat interfluves attributed (195 m):");
    let low_best = attribute(&format!("{} 2048²", best.name), &best.lv[2].result, DOMAIN_KM / 2048.0, Some(&best.lv[2].d5));
    let low_corse = attribute("Corse 195 m", cg195, 200.0 / 1024.0, None);
    let _ = (low_best, low_corse);
    // ── I: the images ──
    let (bp, bplace) = (best.p, best.place);
    let pick = |p: f32, place: &str| -> Option<&Run> {
        runs.iter().find(|r| r.p == p && r.place == place && (place == "none" || r.setting == carried))
    };
    let set1: Vec<&Run> = ["none", "W-phys", "W-tous"].iter().filter_map(|pl| pick(bp, pl)).collect();
    let set2: Vec<&Run> = [2.0f32, 4.0, 6.0].iter().filter_map(|&p| pick(p, bplace)).collect();
    let n = 2048usize;
    let cell_m = DOMAIN_KM / n as f32 * 1000.0;
    let rgba = |r: &Run| -> Vec<u8> {
        let g = &r.lv[2].result;
        let (_, hy) = read(g, None);
        rgba_of(g, Some(&hy), cell_m, &ss)
    };
    let rg1: Vec<Vec<u8>> = set1.iter().map(|r| rgba(r)).collect();
    let rg2: Vec<Vec<u8>> = set2.iter().map(|r| rgba(r)).collect();
    let (_, chy) = read_cell(cg195, None, 200.0 / 1024.0);
    let crgba = rgba_of(cg195, Some(&chy), cell_m, &ss);
    for (name, ox, oy) in CROPS.iter().copied().chain(std::iter::once(("tight_25km", TIGHT.0, TIGHT.1))) {
        let frac = if name == "tight_25km" { 0.0625 } else { 0.25 };
        let corigin = if name == "tight_25km" { CORSE_FLANK } else { CORSE_SPINE };
        let mut rows = Vec::new();
        for rg in [&rg1, &rg2] {
            let mut row: Vec<image::RgbaImage> = rg.iter().map(|x| window(x, n, (ox as f32 / 512.0, oy as f32 / 512.0), frac, 448)).collect();
            row.push(window(&crgba, 1024, corigin, frac * 2.0, 448));
            rows.push(row);
        }
        compose(&rows).save(out().join(format!("f164_{name}.png"))).unwrap();
    }
    let mut row: Vec<image::RgbaImage> = rg1.iter().map(|x| whole(x, n, 640)).collect();
    row.push(whole(&crgba, 1024, 320));
    compose(&[row]).save(out().join("f164_whole_2048.png")).unwrap();
    eprintln!(
        "   images f164_<crop>.png (tight_25km, A_chaine, B_plateau, C_bassin_lac4): row 1 {} then Corsica; row 2 {} then Corsica; f164_whole_2048.png row 1 · rivers order ≥ 2, north up · {NOTICE}",
        set1.iter().map(|r| r.name.as_str()).collect::<Vec<_>>().join(", "),
        set2.iter().map(|r| r.name.as_str()).collect::<Vec<_>>().join(", ")
    );
    let mut series: Vec<(String, Vec<(f32, f32)>, [u8; 3])> = Vec::new();
    let cols = [[120, 120, 120], [200, 60, 40], [40, 110, 200]];
    let ck = DOMAIN_KM / 2048.0;
    for (i, r) in set1.iter().enumerate() {
        series.push((r.name.clone(), r.lv[2].r.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, ck); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect(), cols[i % 3]));
    }
    series.push(("Corse 195 m".into(), co195.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, 0.1953125); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect(), [0, 0, 0]));
    spectra_png(&series, &out().join("f164_spectra.png"));
    eprintln!("   image f164_spectra.png: {} grey / red / blue, Corsica black · {NOTICE}", set1.iter().map(|r| r.name.as_str()).collect::<Vec<_>>().join(", "));
    eprintln!("\n   cost: {} · total {:.0} s", runs.iter().map(|r| format!("{} {:.0} s", r.name, r.secs)).collect::<Vec<_>>().join(" · "), t0.elapsed().as_secs_f64());
    eprintln!("\n==========  end Finding 164 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
