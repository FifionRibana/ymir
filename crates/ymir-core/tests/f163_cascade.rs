//! ADR Finding 163 — the cascade, round 5: the physics level redone (a network less aligned on the grid, more
//! dissected) and what the author sees named (the directions, the sinuosity, the flat interfluves), each instrument
//! validated on its controls first; the physics-level variants measured at 256² against Corsica at 1 563 m, then the
//! best and P0 through F162's R4+π chain to 2 048², judged against Corsica at 195 m.
//! Declared in `docs/reports/relief_method/f163_physics_level/f163_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f163_cascade -- --ignored --exact f163_cascade --nocapture

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, Knobs, PSEED, build_world};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant};
use ymir_core::cascade::hydro::Hydro;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::planform::{
    Planform, WINDOWS_KM, control_plane, control_terraces, control_valleys, planform, sine_reference,
};
use ymir_core::cascade::predict::{Predict, control_fractal, predictability};
use ymir_core::cascade::{CascadeConfig, physics_level, roll};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const PEAK_TARGET_M: f32 = 2850.0;
const LEVELS: [usize; 3] = [512, 1024, 2048];
const CROPS: [(&str, usize, usize); 3] = [("A_chaine", 128, 320), ("B_plateau", 176, 192), ("C_bassin_lac4", 74, 160)];
/// F162's tight 25 km crop (512² origin, 32 cells there).
const TIGHT: (usize, usize) = (176, 352);
const CORSE_SPINE: (f32, f32) = (72.0 / 256.0, 96.0 / 256.0);
const CORSE_FLANK: (f32, f32) = (100.0 / 256.0, 140.0 / 256.0);
/// F162's validated predictability instruments (CV_λ, CV_L, σ_θ, R2_g; Horton and P_s set aside at F162's commit).
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];
const P_CELL: f32 = 1.5625;

fn cache() -> PathBuf {
    root().join("target_f163/f163_cache")
}
fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f163_physics_level/images")
}

fn out15(a: f32, b: f32) -> bool {
    !(a.is_finite() && b.is_finite() && b != 0.0) || a > 1.5 * b || a < b / 1.5
}
fn lnr(a: f32, b: f32) -> f32 {
    if a > 0.0 && b > 0.0 && a.is_finite() && b.is_finite() { (a / b).ln().abs() } else { 0.0 }
}

fn fmt_s(p: &Planform) -> String {
    (0..3)
        .map(|i| {
            if p.sinuosity[i].is_nan() {
                format!("S_{} —", WINDOWS_KM[i])
            } else {
                format!("S_{} {:.3} ({})", WINDOWS_KM[i], p.sinuosity[i], p.windows[i])
            }
        })
        .collect::<Vec<_>>()
        .join(" · ")
}
fn print_pf(tag: &str, p: &Planform) {
    eprintln!(
        "      {tag:<30} A_dir {:.1} % ({} pieces; trunks ≥ 100 km² {:.1} %, {}) · {} · F {:.2} %",
        100.0 * p.a_dir,
        p.pieces,
        100.0 * p.a_dir_trunks,
        p.trunk_pieces,
        fmt_s(p),
        100.0 * p.flat
    );
}

/// F163-I: an instrument's verdict at one cell.
#[derive(Clone, Copy, Debug)]
struct Valid {
    dir: bool,
    sinu: [bool; 3],
    flat: bool,
}

struct Phys {
    name: String,
    pm: GridF32,
    r: Reading,
    hy: Hydro,
    pf: Planform,
    secs: f64,
    steps: usize,
    u0: f32,
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
    secs: [f64; 4],
    drift: f32,
    trunk_p90: f32,
    pi_var: f64,
}

fn run_chain(
    name: &str,
    mut ch: Chain,
    phys_r: &Reading,
    corse: &BTreeMap<usize, Reading>,
    calibrate: bool,
) -> (Vec<Lvl>, f64, Vec<String>) {
    let t0 = Instant::now();
    let mut out: Vec<Lvl> = Vec::new();
    let mut prev_r = phys_r.clone();
    let mut prev = ch.physics_m.clone();
    let mut per_level_stop = Vec::new();
    for n in LEVELS {
        if calibrate {
            let t = Instant::now();
            let target = corse[&(n / 2)].octaves[0];
            let (k, trials) = ch.calibrate_peak(target, &|| false);
            eprintln!(
                "      {name} calibration {n}²: target {target:.2} m → k {k:.3e} ({} trials, {:.0} s)",
                trials.len(),
                t.elapsed().as_secs_f64()
            );
        }
        let l = ch.next_level(&|| false).unwrap().clone();
        let result = l.result();
        let (r, _) = read(&result, Some(&l.d5));
        let cell_km = DOMAIN_KM / n as f32;
        let p = (n >= 1024).then(|| predictability(&result, cell_km, Some(&l.d5)));
        let pf = (n >= 1024).then(|| planform(&result, cell_km, Some(&l.d5)));
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
            per_level_stop.push(format!("{n}²: {}", why.join("; ")));
        }
        eprintln!(
            "      {name} {n}² ({:.1} s) · k {:.3e} · peak {:.0} · drift {:+.2} % · trunk p90 {:.2} · octave 0 {:.2} m (Corse {:.2}){}",
            l.secs.iter().sum::<f64>(),
            l.k,
            r.peak,
            100.0 * rel,
            tp,
            r.octaves[0],
            corse[&(n / 2)].octaves[0],
            if why.is_empty() { String::new() } else { format!(" · per-level R: {}", why.join("; ")) }
        );
        out.push(Lvl {
            n,
            result: result.clone(),
            d5: l.d5.clone(),
            r: r.clone(),
            p,
            pf,
            k: l.k,
            secs: l.secs,
            drift: rel,
            trunk_p90: tp,
            pi_var,
        });
        prev_r = r;
        prev = result;
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    (out, t0.elapsed().as_secs_f64(), per_level_stop)
}

/// R at 2 048² (F163-R): F162's rule plus the three new instruments whose control passed.
#[allow(clippy::too_many_arguments)]
fn reject_2048(
    l: &Lvl,
    co: &Reading,
    cp: &Predict,
    cpf: &Planform,
    v: &Valid,
    ctrl_ok: bool,
    pi_share: f64,
) -> Vec<String> {
    let r = &l.r;
    let mut why = Vec::new();
    if ctrl_ok && r.facets > 1.5 * co.facets && r.facet_cells > 20 {
        why.push(format!("facets {:.2} > 1.5 × {:.2}", r.facets, co.facets));
    }
    if r.walls28 > 1.5 * co.walls28 && r.walls28_cells > 20 {
        why.push(format!("walls 28° {:.3} > 1.5 × {:.3}", r.walls28, co.walls28));
    }
    if out15(r.slope_p50, co.slope_p50) {
        why.push(format!("slope p50 {:.3} / {:.3}", r.slope_p50, co.slope_p50));
    }
    if out15(r.slope_p90, co.slope_p90) {
        why.push(format!("slope p90 {:.3} / {:.3}", r.slope_p90, co.slope_p90));
    }
    if out15(r.lambda_km, co.lambda_km) {
        why.push(format!("λ {:.2} / {:.2} km", r.lambda_km, co.lambda_km));
    }
    let cell_km = DOMAIN_KM / l.n as f32;
    for j in 0..r.octaves.len().min(co.octaves.len()) {
        let (a, b) = octave_km(j, cell_km);
        if a >= 0.39 && b <= 12.6 && out15(r.octaves[j], co.octaves[j]) {
            why.push(format!("octave {a:.2}–{b:.2} km {:.1} / {:.1} m", r.octaves[j], co.octaves[j]));
        }
    }
    if let Some(p) = &l.p {
        for (i, ((name, x, _), (_, c, _))) in p.values().iter().zip(cp.values().iter()).enumerate() {
            if VALID_P[i] && out15(*x, *c) {
                why.push(format!("{name} {x:.3} / {c:.3}"));
            }
        }
    }
    if pi_share > 0.10 {
        why.push(format!("π share {:.1} %", 100.0 * pi_share));
    }
    if let Some(pf) = &l.pf {
        if v.dir && out15(pf.a_dir, cpf.a_dir) {
            why.push(format!("A_dir {:.1} / {:.1} %", 100.0 * pf.a_dir, 100.0 * cpf.a_dir));
        }
        for i in 0..3 {
            if v.sinu[i] && out15(pf.excess(i), cpf.excess(i)) {
                why.push(format!("S_{} − 1 {:.3} / {:.3}", WINDOWS_KM[i], pf.excess(i), cpf.excess(i)));
            }
        }
        if v.flat && out15(pf.flat, cpf.flat) {
            why.push(format!("F {:.2} / {:.2} %", 100.0 * pf.flat, 100.0 * cpf.flat));
        }
    }
    why
}

#[test]
#[ignore]
fn f163_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    std::fs::create_dir_all(cache()).unwrap();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 163 . the cascade, round 5: the physics level redone; directions, sinuosity, flat interfluves  ==========");

    // ── C: Corsica ──
    let mut corse: BTreeMap<usize, Reading> = BTreeMap::new();
    let mut corse_g: BTreeMap<usize, GridF32> = BTreeMap::new();
    let mut corse_pf: BTreeMap<usize, Planform> = BTreeMap::new();
    let mut corse_p: BTreeMap<usize, Predict> = BTreeMap::new();
    for nc in [128usize, 256, 512, 1024] {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).expect("data/corsica (prep_corse.py)");
        let ck = 200.0 / nc as f32;
        let (r, _) = read_cell(&g, None, ck);
        if nc != 256 {
            let pf = planform(&g, ck, None);
            print_pf(&format!("Corse {:.0} m", ck * 1000.0), &pf);
            corse_pf.insert(nc, pf);
        }
        if nc >= 512 {
            corse_p.insert(nc, predictability(&g, ck, None));
        }
        corse.insert(nc, r);
        corse_g.insert(nc, g);
    }

    // ── I: the controls, at 256² (1 563 m) and 512² (195 m) ──
    eprintln!("\n   I — the controls:");
    let mut valid: BTreeMap<usize, Valid> = BTreeMap::new();
    for (cn, cc, corse_n) in [(256usize, P_CELL, 128usize), (512, 0.1953125f32, 1024)] {
        let plane = planform(&control_plane(cn, cc), cc, None);
        let frac = planform(&control_fractal(cn, cc, 7), cc, None);
        let terr = planform(&control_terraces(cn, cc), cc, None);
        let s0 = planform(&control_valleys(cn, cc, 0.0, None), cc, None);
        let s22 = planform(&control_valleys(cn, cc, 22.5, None), cc, None);
        let s30 = planform(&control_valleys(cn, cc, 30.0, None), cc, None);
        let sine = planform(&control_valleys(cn, cc, 0.0, Some((4.0, 32.0))), cc, None);
        let co = &corse_pf[&corse_n];
        eprintln!("   at {cn}² ({:.0} m), against Corsica at the same cell:", cc * 1000.0);
        print_pf("D8 plane at 10° (high)", &plane);
        print_pf("fractal (reported)", &frac);
        print_pf("terraced field (high F)", &terr);
        print_pf("straight valley 0°", &s0);
        print_pf("straight valley 22.5°", &s22);
        print_pf("straight valley 30°", &s30);
        print_pf("sine valley (a 4, P 32 cells)", &sine);
        eprintln!(
            "      sine reference (the continuous curve): {}",
            WINDOWS_KM.iter().map(|w| format!("S_{w} {:.3}", sine_reference(4.0, 32.0, w / cc))).collect::<Vec<_>>().join(" · ")
        );
        print_pf("Corsica", co);
        let dir = plane.a_dir >= 1.5 * co.a_dir;
        let mut sinu = [false; 3];
        for i in 0..3 {
            let ce = co.excess(i);
            sinu[i] = ce.is_finite() && s0.excess(i) <= ce / 1.5 && s22.excess(i) <= ce / 1.5;
        }
        let flat = terr.flat >= 1.5 * co.flat;
        eprintln!(
            "      → A_dir {} · S−1 {} · F {}",
            if dir { "PASSES" } else { "FAILS (blind)" },
            (0..3).map(|i| format!("{} km {}", WINDOWS_KM[i], if sinu[i] { "passes" } else if co.excess(i).is_nan() { "not read" } else { "FAILS" })).collect::<Vec<_>>().join(", "),
            if flat { "PASSES" } else { "FAILS (blind)" }
        );
        valid.insert(corse_n, Valid { dir, sinu, flat });
    }
    let v1563 = valid[&128];
    let v195 = valid[&1024];
    valid.insert(512, v195);

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

    // ── D: the physics-level variants at 256² ──
    let base = CascadeConfig { seed: PSEED, ..CascadeConfig::declared(ss.depth_scale_m as f32) };
    let co128 = &corse[&128];
    let cpf128 = &corse_pf[&128];
    let rough_a = 0.5 * co128.octaves[0];
    let fctl = facet_control(256, P_CELL);
    let fctl_ok = fctl.0 >= 60.0 && fctl.1 <= 20.0;
    eprintln!(
        "\n   D — the physics level at 256² (U₀ calibrated on {PEAK_TARGET_M:.0} m); Corsica 1563 m: λ {:.2} km · facets {:.2} % · octaves {:.1} / {:.1} m · roughness A = 0.5 × {:.1} = {:.1} m · facet control pyramids {:.1} % / paraboloids {:.1} %",
        co128.lambda_km, co128.facets, co128.octaves[0], co128.octaves[1], co128.octaves[0], rough_a, fctl.0, fctl.1
    );
    let mut lifted: Vec<bool> = Vec::new();
    let mut run_phys = |name: &str, cfg: &CascadeConfig| -> Phys {
        let (rec, cal) = physics_level(&coarse, 256, cfg, PEAK_TARGET_M, 300, &mut |_| {}, &|| false).unwrap();
        let pm = to_m(&rec.z, n2m);
        let (r, hy) = read(&pm, Some(&cal.lifted_mask));
        let pf = planform(&pm, P_CELL, Some(&cal.lifted_mask));
        lifted = cal.lifted_mask.clone();
        let p = Phys { name: name.to_string(), pm, r, hy, pf, secs: cal.secs, steps: rec.steps, u0: cal.u0 };
        eprintln!(
            "   ── {name:<22} ({:.1} s, {} steps, U₀ {:.3e}) peak {:.0} · λ {:.2} km · facets {:.2} % · octaves {:.1} / {:.1} m · slope p50/p90 {:.3}/{:.3} · density {:.3}",
            p.secs, p.steps, p.u0, p.r.peak, p.r.lambda_km, p.r.facets, p.r.octaves[0], p.r.octaves[1], p.r.slope_p50, p.r.slope_p90, p.r.density
        );
        print_pf(name, &p.pf);
        p
    };
    let measures = |p: &Phys| -> Vec<(&'static str, f32, f32, bool)> {
        vec![
            ("λ", p.r.lambda_km, co128.lambda_km, true),
            ("facets", p.r.facets as f32, co128.facets as f32, fctl_ok),
            ("oct 3.1–6.25", p.r.octaves[0], co128.octaves[0], true),
            ("oct 6.25–12.5", p.r.octaves[1], co128.octaves[1], true),
            ("A_dir", p.pf.a_dir, cpf128.a_dir, v1563.dir),
            ("S_10 − 1", p.pf.excess(2), cpf128.excess(2), v1563.sinu[2]),
            ("F", p.pf.flat, cpf128.flat, v1563.flat),
        ]
    };
    let score = |p: &Phys| -> (usize, f32) {
        let m = measures(p);
        (m.iter().filter(|x| x.3 && out15(x.1, x.2)).count(), m.iter().filter(|x| x.3).map(|x| lnr(x.1, x.2)).sum())
    };
    let mut phys: Vec<(Phys, CascadeConfig)> = Vec::new();
    let p0 = run_phys("P0", &base);
    phys.push((p0, base.clone()));
    let c_rr = CascadeConfig { receiver_tau: Some(0.5), ..base.clone() };
    phys.push((run_phys("P-mfd (receiver)", &c_rr), c_rr.clone()));
    let c_b = CascadeConfig { roughness_m: rough_a, ..base.clone() };
    phys.push((run_phys("P-bruit", &c_b), c_b.clone()));
    // P-dissection: the MFD exponent scan
    let mut scan: Vec<(Option<f32>, f32)> = vec![(Some(2.0), phys[0].0.r.lambda_km)];
    for p in [Some(1.0f32), Some(1.5), Some(3.0), Some(4.0), Some(6.0), None] {
        let c = CascadeConfig { mfd_exponent: p, ..base.clone() };
        let name = p.map_or("p-scan D8".to_string(), |v| format!("p-scan p = {v}"));
        let ph = run_phys(&name, &c);
        scan.push((p, ph.r.lambda_km));
        phys.push((ph, c));
    }
    let best_p = scan.iter().min_by(|a, b| lnr(a.1, co128.lambda_km).partial_cmp(&lnr(b.1, co128.lambda_km)).unwrap()).unwrap().0;
    let dis_reached = best_p != Some(2.0);
    let best_lam = scan.iter().find(|s| s.0 == best_p).unwrap().1;
    eprintln!(
        "   P-dissection: λ by p {} → kept p = {} (λ {:.2} km, {} of 9.08 km ×1.25){}",
        scan.iter().map(|s| format!("{}: {:.2}", s.0.map_or("D8".into(), |v| format!("{v}")), s.1)).collect::<Vec<_>>().join(" · "),
        best_p.map_or("D8".into(), |v| format!("{v}")),
        best_lam,
        if (best_lam / co128.lambda_km).ln().abs() <= 1.25f32.ln() { "within" } else { "NOT within" },
        if dis_reached { "" } else { " · NOT REACHED: P0's p is the nearest, P-dissection leaves the combinations" }
    );
    let mut combos: Vec<(String, CascadeConfig)> = vec![("P-mfd + P-bruit".into(), CascadeConfig { receiver_tau: Some(0.5), roughness_m: rough_a, ..base.clone() })];
    if dis_reached {
        combos.push(("P-mfd + P-dis".into(), CascadeConfig { receiver_tau: Some(0.5), mfd_exponent: best_p, ..base.clone() }));
        combos.push(("P-bruit + P-dis".into(), CascadeConfig { roughness_m: rough_a, mfd_exponent: best_p, ..base.clone() }));
        combos.push(("P-mfd + P-bruit + P-dis".into(), CascadeConfig { receiver_tau: Some(0.5), roughness_m: rough_a, mfd_exponent: best_p, ..base.clone() }));
    }
    for (name, c) in &combos {
        phys.push((run_phys(name, c), c.clone()));
    }
    let _ = &lifted;
    // the roughness's share
    let z0 = &phys[0].0.pm;
    for (p, c) in phys.iter().filter(|(_, c)| c.roughness_m > 0.0) {
        let band = (p.r.octaves[0] as f64).powi(2) + (p.r.octaves[1] as f64).powi(2);
        let a = (c.roughness_m as f64).powi(2) / band;
        let (mut sd, mut sz, mut mz, mut md, mut cnt) = (0f64, 0f64, 0f64, 0f64, 0usize);
        for k in 0..z0.data.len() {
            if p.pm.data[k] > 0.0 && z0.data[k] > 0.0 && !lifted[k] {
                md += (p.pm.data[k] - z0.data[k]) as f64;
                mz += p.pm.data[k] as f64;
                cnt += 1;
            }
        }
        md /= cnt.max(1) as f64;
        mz /= cnt.max(1) as f64;
        for k in 0..z0.data.len() {
            if p.pm.data[k] > 0.0 && z0.data[k] > 0.0 && !lifted[k] {
                sd += ((p.pm.data[k] - z0.data[k]) as f64 - md).powi(2);
                sz += (p.pm.data[k] as f64 - mz).powi(2);
            }
        }
        eprintln!("   roughness share, {}: (a) A² / final 3.1–12.5 km band {:.1} % · (b) var(z − z_P0) / var(z) {:.1} %", p.name, 100.0 * a, 100.0 * sd / sz.max(1e-9));
    }
    // the ranking at 256²
    eprintln!("\n   the ranking at 256² (measures outside ×1.5 of Corsica 1563 m, among the validated; Σ |ln ratio|):");
    let mut ranked: Vec<(usize, f32, usize)> = Vec::new();
    for (i, (p, _)) in phys.iter().enumerate() {
        let peak_ok = p.r.peak >= 2700.0 && p.r.peak <= 3000.0;
        let (o, s) = score(p);
        let m = measures(p);
        eprintln!(
            "   {:<26} {} out · Σ {:.2} · {}{}",
            p.name,
            o,
            s,
            m.iter().map(|x| format!("{} {:.3}/{:.3}{}{}", x.0, x.1, x.2, if x.3 { "" } else { " (blind)" }, if x.3 && out15(x.1, x.2) { " ✗" } else { "" })).collect::<Vec<_>>().join(" · "),
            if peak_ok { String::new() } else { format!(" · EXCLUDED: peak {:.0} m", p.r.peak) }
        );
        if peak_ok {
            ranked.push((o, s, i));
        }
    }
    ranked.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.partial_cmp(&b.1).unwrap()));
    let bi = ranked.iter().find(|r| r.2 != 0).map_or(0, |r| r.2);
    let p0_rank = ranked.iter().position(|r| r.2 == 0);
    eprintln!(
        "   → best variant: {} (P0 ranks {}){}",
        phys[bi].0.name,
        p0_rank.map_or("excluded".into(), |r| format!("{}", r + 1)),
        if ranked.first().is_some_and(|r| r.2 == 0) { " · NOTE: P0 itself ranks first; the best non-P0 variant runs as declared" } else { "" }
    );

    // ── images: the physics level alone, each variant ──
    let cell_m = P_CELL * 1000.0;
    let mut tiles: Vec<image::RgbaImage> = Vec::new();
    for (p, _) in &phys {
        tiles.push(whole(&rgba_of(&p.pm, Some(&p.hy), cell_m, &ss), 256, 384));
    }
    let cg = &corse_g[&128];
    let (_, chy) = read_cell(cg, None, P_CELL);
    tiles.push(whole(&rgba_of(cg, Some(&chy), cell_m, &ss), 128, 192));
    let rows: Vec<Vec<image::RgbaImage>> = tiles.chunks(4).map(|c| c.to_vec()).collect();
    compose(&rows).save(out().join("f163_physics_256.png")).unwrap();
    eprintln!(
        "   image f163_physics_256.png: {} then Corsica 1563 m (same km per pixel), rivers order ≥ 2 · {NOTICE}",
        phys.iter().map(|p| p.0.name.as_str()).collect::<Vec<_>>().join(", ")
    );

    // ── the chains ──
    let pi_table: Vec<(usize, f32)> = LEVELS.iter().map(|&n| (n, 0.15 * corse[&(n / 2)].octaves[0])).collect();
    let r4pi = || AmpConfig { talus_fine_scale: 0.25, n1_recal_depth: 2, n1_pi_m: pi_table.clone(), ..AmpConfig::declared(PSEED, Variant::N1, 1.0) };
    let mk = |pm: &GridF32, c: AmpConfig| Chain::new(pm.clone(), c).with_d5(lifted.clone());
    let mut f161 = AmpConfig::declared(PSEED, Variant::N1, 1.0);
    f161.n1_k = vec![(512, 1.0e3), (1024, 1.0e3), (2048, 177.8)];
    eprintln!("\n   chain N1 (F161 as committed) from P0:");
    let (n1_lv, n1_secs, n1_stop) = run_chain("N1", mk(&phys[0].0.pm, f161), &phys[0].0.r, &corse, false);
    eprintln!("\n   chain R4+π from P0:");
    let (p0_lv, p0_secs, p0_stop) = run_chain("R4+π·P0", mk(&phys[0].0.pm, r4pi()), &phys[0].0.r, &corse, true);
    let best_name = format!("R4+π·{}", phys[bi].0.name);
    eprintln!("\n   chain {best_name}:");
    let (b_lv, b_secs, b_stop) = run_chain(&best_name, mk(&phys[bi].0.pm, r4pi()), &phys[bi].0.r, &corse, true);
    let chains: Vec<(String, &Vec<Lvl>, f64, &Vec<String>)> =
        vec![("N1".into(), &n1_lv, n1_secs, &n1_stop), ("R4+π·P0".into(), &p0_lv, p0_secs, &p0_stop), (best_name.clone(), &b_lv, b_secs, &b_stop)];

    // ── the readings at 1024² and 2048² ──
    eprintln!("\n   the readings (1024² against Corsica 391 m, 2048² against Corsica 195 m):");
    for (name, lv, _, _) in &chains {
        for l in lv.iter().filter(|l| l.n >= 1024) {
            let cn = l.n / 2;
            let co = &corse[&cn];
            let r = &l.r;
            let cell_km = DOMAIN_KM / l.n as f32;
            eprintln!(
                "   ── {name} {}²: slope p50/p90 {:.3}/{:.3} ({:.3}/{:.3}) · facets {:.2} % ({:.2}) · walls 28° {:.3} % ({:.3}) · λ {:.2} km ({:.2}) · peak {:.0} · k {:.3e}",
                l.n, r.slope_p50, r.slope_p90, co.slope_p50, co.slope_p90, r.facets, co.facets, r.walls28, co.walls28, r.lambda_km, co.lambda_km, r.peak, l.k
            );
            eprintln!(
                "      spectrum (m), cascade / Corse: {}",
                (0..r.octaves.len().min(co.octaves.len()).min(6))
                    .map(|j| {
                        let (a, b) = octave_km(j, cell_km);
                        format!("[{a:.2}–{b:.2} km] {:.1} / {:.1}", r.octaves[j], co.octaves[j])
                    })
                    .collect::<Vec<_>>()
                    .join(" · ")
            );
            if let Some(p) = &l.p {
                let cp = &corse_p[&cn];
                eprintln!(
                    "      predictability: {}",
                    p.values().iter().zip(cp.values().iter()).enumerate().filter(|(i, _)| VALID_P[*i]).map(|(_, ((n, a, _), (_, c, _)))| format!("{n} {a:.3} ({c:.3})")).collect::<Vec<_>>().join(" · ")
                );
            }
            if let Some(pf) = &l.pf {
                print_pf(&format!("{name} {}²", l.n), pf);
                print_pf(&format!("Corse {:.0} m", cell_km * 1000.0), &corse_pf[&cn]);
            }
        }
    }
    // ── R at 2048² ──
    let co195 = &corse[&1024];
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    eprintln!("\n   R at 2048² (facet control: pyramids {:.1} % · paraboloids {:.1} %):", ctrl.0, ctrl.1);
    for (name, lv, secs, stops) in &chains {
        let l = &lv[2];
        let cell_km = DOMAIN_KM / l.n as f32;
        let band_var: f64 = (0..l.r.octaves.len())
            .filter(|&j| {
                let (a, b) = octave_km(j, cell_km);
                a >= 0.39 && b <= 12.6
            })
            .map(|j| (l.r.octaves[j] as f64).powi(2))
            .sum();
        let share = lv.iter().map(|x| x.pi_var).sum::<f64>() / band_var.max(1e-9);
        let why = reject_2048(l, co195, &corse_p[&1024], &corse_pf[&1024], &v195, ctrl_ok, share);
        eprintln!(
            "   {name:<24} ({:.0} s) π share {:.2} % · {}",
            secs,
            100.0 * share,
            if why.is_empty() { "PASSES R at 2048²".to_string() } else { format!("REJECTED: {}", why.join("; ")) }
        );
        if !stops.is_empty() {
            eprintln!("      per-level R: {}", stops.join(" | "));
        }
    }

    // ── images at 2048²: P0 and the best, rivers drawn, beside Corsica ──
    let shown: [(&str, &Vec<Lvl>); 2] = [("R4+π·P0", &p0_lv), (best_name.as_str(), &b_lv)];
    let n = 2048usize;
    let cell_m = DOMAIN_KM / n as f32 * 1000.0;
    let rgbas: Vec<Vec<u8>> = shown
        .iter()
        .map(|(_, lv)| {
            let g = &lv[2].result;
            let (_, hy) = read(g, None);
            rgba_of(g, Some(&hy), cell_m, &ss)
        })
        .collect();
    let cg = &corse_g[&1024];
    let (_, chy) = read_cell(cg, None, 200.0 / 1024.0);
    let crgba = rgba_of(cg, Some(&chy), cell_m, &ss);
    let mut row = Vec::new();
    for rg in &rgbas {
        row.push(whole(rg, n, 768));
    }
    row.push(whole(&crgba, 1024, 384));
    compose(&[row]).save(out().join("f163_whole_2048.png")).unwrap();
    for (name, ox, oy) in CROPS.iter().copied().chain(std::iter::once(("tight_25km", TIGHT.0, TIGHT.1))) {
        let frac = if name == "tight_25km" { 0.0625 } else { 0.25 };
        let corigin = if name == "tight_25km" { CORSE_FLANK } else { CORSE_SPINE };
        let mut row = Vec::new();
        for rg in &rgbas {
            row.push(window(rg, n, (ox as f32 / 512.0, oy as f32 / 512.0), frac, 512));
        }
        row.push(window(&crgba, 1024, corigin, frac * 2.0, 512));
        compose(&[row]).save(out().join(format!("f163_{name}.png"))).unwrap();
    }
    eprintln!(
        "   images f163_whole_2048.png and f163_<crop>.png (tight_25km, A_chaine, B_plateau, C_bassin_lac4): {} then Corsica 195 m (its spine, or its western flank for the 25 km crop), rivers order ≥ 2, north up · {NOTICE}",
        shown.iter().map(|s| s.0).collect::<Vec<_>>().join(", ")
    );
    let mut series: Vec<(String, Vec<(f32, f32)>, [u8; 3])> = Vec::new();
    let cols = [[120, 120, 120], [200, 60, 40], [40, 110, 200]];
    for (i, (name, lv, _, _)) in chains.iter().enumerate() {
        let l = &lv[2];
        let ck = DOMAIN_KM / 2048.0;
        series.push((
            name.clone(),
            l.r.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, ck); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect(),
            cols[i],
        ));
    }
    series.push((
        "Corse 195 m".into(),
        co195.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, 0.1953125); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect(),
        [0, 0, 0],
    ));
    spectra_png(&series, &out().join("f163_spectra.png"));
    eprintln!("   image f163_spectra.png: the 2048² octaves, N1 grey, R4+π·P0 red, {best_name} blue, Corsica (195 m) black · {NOTICE}");
    eprintln!("\n==========  end Finding 163 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// F163-I: the physics level alone (256²), the land cropped and enlarged ×3 (nearest), rivers order ≥ 2, beside Corsica
/// at 1 563 m at the same scale. Reruns the six variants that matter (each a few seconds) from the bench's cache.
///
/// Run (after `f163_cascade`): cargo test -p ymir-core --release --test f163_cascade -- --ignored --exact f163_physics_crops --nocapture
#[test]
#[ignore]
fn f163_physics_crops() {
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let coarse = roll(&load_grid(&cache().join("coarse64.bin")).expect("run f163_cascade first"), off);
    let cg = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let (cr, chy) = read_cell(&cg, None, P_CELL);
    let rough_a = 0.5 * cr.octaves[0];
    let base = CascadeConfig { seed: PSEED, ..CascadeConfig::declared(ss.depth_scale_m as f32) };
    let variants: Vec<(&str, CascadeConfig)> = vec![
        ("P0", base.clone()),
        ("P-mfd", CascadeConfig { receiver_tau: Some(0.5), ..base.clone() }),
        ("P-bruit", CascadeConfig { roughness_m: rough_a, ..base.clone() }),
        ("P-dis p6", CascadeConfig { mfd_exponent: Some(6.0), ..base.clone() }),
        ("P-mfd + P-dis", CascadeConfig { receiver_tau: Some(0.5), mfd_exponent: Some(6.0), ..base.clone() }),
        ("P-bruit + P-dis", CascadeConfig { roughness_m: rough_a, mfd_exponent: Some(6.0), ..base.clone() }),
    ];
    let cell_m = P_CELL * 1000.0;
    let mut tiles = Vec::new();
    for (name, c) in &variants {
        let (rec, cal) = physics_level(&coarse, 256, c, PEAK_TARGET_M, 300, &mut |_| {}, &|| false).unwrap();
        let pm = to_m(&rec.z, n2m);
        let (r, hy) = read(&pm, Some(&cal.lifted_mask));
        eprintln!("   {name}: peak {:.0} m · λ {:.2} km", r.peak, r.lambda_km);
        tiles.push(window(&rgba_of(&pm, Some(&hy), cell_m, &ss), 256, (24.0 / 256.0, 46.0 / 256.0), 168.0 / 256.0, 504));
    }
    let rows: Vec<Vec<image::RgbaImage>> = vec![tiles[..3].to_vec(), tiles[3..].to_vec()];
    compose(&rows).save(out().join("f163_physics_256_crops.png")).unwrap();
    let c = window(&rgba_of(&cg, Some(&chy), cell_m, &ss), 128, (0.0, 0.0), 1.0, 384);
    compose(&[vec![c]]).save(out().join("f163_physics_256_corse.png")).unwrap();
    eprintln!(
        "   images f163_physics_256_crops.png (row 1: {}; row 2: {}), the land 168 cells (262 km) ×3, and f163_physics_256_corse.png (Corsica 1563 m, 200 km, same scale) · {NOTICE}",
        variants[..3].iter().map(|v| v.0).collect::<Vec<_>>().join(", "),
        variants[3..].iter().map(|v| v.0).collect::<Vec<_>>().join(", ")
    );
}
