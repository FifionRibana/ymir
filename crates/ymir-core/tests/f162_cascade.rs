//! ADR Finding 162 — the cascade, round 4: the predictability measured against Corsica (six instruments, each
//! validated on a regular and a dispersed control), remedied by a noisy hardness (ρ) and an initial perturbation (π),
//! and the 1.5–6 km deficit filled by a 4×4 retargeting (R4); judged at 2048² against Corsica at 195 m.
//! Declared in `docs/reports/relief_method/f162_predictability/f162_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f162_cascade -- --ignored --exact f162_cascade --nocapture

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, Knobs, PSEED, build_world, sorted};
use std::collections::BTreeMap;
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, ocean_mask};
use ymir_core::cascade::measure::{octave_km, octave_rms};
use ymir_core::cascade::predict::{Predict, control_fractal, control_regular, predictability};
use ymir_core::cascade::{CascadeConfig, physics_level, roll};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const PEAK_TARGET_M: f32 = 2850.0;
const LEVELS: [usize; 3] = [512, 1024, 2048];
const CROPS: [(&str, usize, usize); 3] = [("A_chaine", 128, 320), ("B_plateau", 176, 192), ("C_bassin_lac4", 74, 160)];
/// The tight 25 km crop (512² origin, 32 cells there): the range's southern flank in crop A.
const TIGHT: (usize, usize) = (176, 352);
/// Corsica's windows (fractions of its 200 km grid): the spine (100 km) and a western flank (25 km).
const CORSE_SPINE: (f32, f32) = (72.0 / 256.0, 96.0 / 256.0);
const CORSE_FLANK: (f32, f32) = (100.0 / 256.0, 140.0 / 256.0);

fn print_predict(tag: &str, p: &Predict) {
    eprintln!(
        "      {tag:<26} CV_λ {:.3} ({} gaps) · CV_L {:.3} [o1 {:.3} o2 {:.3} o3 {:.3}] · R_b {:.2} CV {:.3} · R_l {:.2} CV {:.3} ({} basins) · θ {:.1}° σ {:.1}° ({} junctions) · P_s {:.3} dex · R2_g {:.3} ({} windows)",
        p.cv_lambda, p.n_gaps, p.cv_len, p.cv_len_by_order[0], p.cv_len_by_order[1], p.cv_len_by_order[2], p.rb_mean, p.cv_rb, p.rl_mean, p.cv_rl, p.n_basins, p.theta_mean, p.theta_std, p.n_junctions, p.peakiness, p.r2_gully, p.n_windows
    );
}

/// The ratio of a value to Corsica's, and whether it is outside ×1.5.
fn out15(a: f32, b: f32) -> bool {
    !(a.is_finite() && b.is_finite() && b != 0.0) || a > 1.5 * b || a < b / 1.5
}

#[allow(dead_code)]
struct Lvl {
    n: usize,
    result: GridF32,
    d5: Vec<bool>,
    r: Reading,
    p: Option<Predict>,
    k: f32,
    secs: [f64; 4],
    iters: [usize; 3],
    drift: f32,
    trunk_p90: f32,
    pi_var: f64,
}

#[allow(clippy::too_many_arguments)]
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
                "      {name} calibration {n}²: target {target:.2} m → k {k:.3e} ({} trials, {:.0} s): {}",
                trials.len(),
                t.elapsed().as_secs_f64(),
                trials.iter().map(|(k, v)| format!("{k:.2e}→{v:.2}")).collect::<Vec<_>>().join(", ")
            );
        }
        let l = ch.next_level(&|| false).unwrap().clone();
        let result = l.result();
        let (r, _) = read(&result, Some(&l.d5));
        let p = (n >= 1024).then(|| predictability(&result, DOMAIN_KM / n as f32, Some(&l.d5)));
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
            "      {name} {n}² ({:.1} s: E {:.1} · T {:.1} · D {:.1}) · k {:.3e} · peak {:.0} · drift {:+.2} % · trunk p90 {:.2} · octave 0 {:.2} m (Corse {:.2}){}",
            l.secs.iter().sum::<f64>(),
            l.secs[1],
            l.secs[2],
            l.secs[3],
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
            k: l.k,
            secs: l.secs,
            iters: [l.done.erosion, l.done.talus, l.done.deposit],
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

/// R at 2 048² (F162-R): the reasons a variant is rejected.
fn reject_2048(l: &Lvl, co: &Reading, cp: &Predict, valid: &[bool; 7], ctrl_ok: bool, pi_share: f64) -> Vec<String> {
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
        for (i, ((name, v, _), (_, c, _))) in p.values().iter().zip(cp.values().iter()).enumerate() {
            if valid[i] && out15(*v, *c) {
                why.push(format!("{name} {v:.3} / {c:.3}"));
            }
        }
    }
    if pi_share > 0.10 {
        why.push(format!("π share {:.1} %", 100.0 * pi_share));
    }
    why
}

#[test]
#[ignore]
fn f162_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let n2m = 2.0 * 1.13 * ss.depth_scale_m as f32;
    std::fs::create_dir_all(cache_dir()).unwrap();
    std::fs::create_dir_all(out_dir()).unwrap();
    eprintln!("\n==========  Finding 162 . the cascade, round 4: predictability against Corsica; ρ, π, R4  ==========");
    // ── P: the controls ──
    let cc = 0.1953125f32;
    let reg = predictability(&control_regular(512, cc), cc, None);
    let dis = predictability(&control_fractal(512, cc, 7), cc, None);
    eprintln!("\n   P — the controls (512², 195 m):");
    print_predict("regular (fishbone)", &reg);
    print_predict("dispersed (fractal)", &dis);
    let mut valid = [false; 7];
    for (i, ((name, rv, inv), (_, dv, _))) in reg.values().iter().zip(dis.values().iter()).enumerate() {
        // « more regular »: a smaller CV / σ / … except P_s and R2_g, where regular reads LARGER
        let ok = if *inv { *rv >= 1.5 * *dv } else { *rv * 1.5 <= *dv };
        valid[i] = ok && rv.is_finite() && dv.is_finite();
        eprintln!("      {name:<6} regular {rv:.3} · dispersed {dv:.3} → {}", if valid[i] { "PASSES (×1.5)" } else { "FAILS: blind, out of R" });
    }
    // ── C: Corsica ──
    let mut corse: BTreeMap<usize, Reading> = BTreeMap::new();
    let mut corse_p: BTreeMap<usize, Predict> = BTreeMap::new();
    let mut corse_g: BTreeMap<usize, GridF32> = BTreeMap::new();
    for nc in [256usize, 512, 1024, 4096] {
        let t = Instant::now();
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).expect("data/corsica (prep_corse.py)");
        let ck = 200.0 / nc as f32;
        let (r, _) = read_cell(&g, None, ck);
        let p = predictability(&g, ck, None);
        eprintln!(
            "   Corse {:.0} m ({:.0} s): slope p50/p90 {:.3}/{:.3} · facets {:.2} % · λ {:.2} km · R8 {:.4} · coherent {:.2} %",
            ck * 1000.0,
            t.elapsed().as_secs_f64(),
            r.slope_p50,
            r.slope_p90,
            r.facets,
            r.lambda_km,
            r.r8,
            r.coherent
        );
        print_predict(&format!("Corse {:.0} m", ck * 1000.0), &p);
        corse.insert(nc, r);
        corse_p.insert(nc, p);
        if nc != 4096 {
            corse_g.insert(nc, g);
        }
    }
    // ── the témoin ──
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let cpath = cache_dir().join("coarse64.bin");
    if !(cpath.exists() && [1024usize, 2048].iter().all(|n| cache_dir().join(format!("ref{n}.bin")).exists())) {
        let t = Instant::now();
        let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
        let wd = build_world(temoin, None, PSEED, None);
        let coarse = c1_coarse_normalized_altitude(&wd.state, &IsostasyConfig::c1_default(), &ss, wd.cfg.target_land_fraction);
        save_grid(&cpath, &coarse);
        for n in [1024usize, 2048] {
            save_grid(&cache_dir().join(format!("ref{n}.bin")), &block_mean(&wd.heightmap, n));
        }
        drop(wd);
        eprintln!("   the témoin built and reduced ({:.0} s)", t.elapsed().as_secs_f64());
    }
    for n in [1024usize, 2048] {
        let g = to_m(&load_grid(&cache_dir().join(format!("ref{n}.bin"))).unwrap(), n2m);
        print_predict(&format!("témoin {n}²"), &predictability(&g, DOMAIN_KM / n as f32, None));
    }
    let coarse = roll(&load_grid(&cpath).unwrap(), off);
    // ── the physics level ──
    let ccfg = CascadeConfig::declared(ss.depth_scale_m as f32);
    let (_rec, cal) = physics_level(&coarse, 256, &ccfg, PEAK_TARGET_M, 300, &mut |_| {}, &|| false).unwrap();
    let pm = to_m(&_rec.z, n2m);
    let (pr, _) = read(&pm, Some(&cal.lifted_mask));
    eprintln!("\n   physics P256: U₀ {:.3e} · {} steps · peak {:.0} m", cal.u0, _rec.steps, pr.peak);
    let base_cfg = |v: Variant| {
        let mut c = AmpConfig::declared(PSEED, v, 1.0);
        c.talus_fine_scale = 0.25;
        c
    };
    let mk = |c: AmpConfig| Chain::new(pm.clone(), c).with_d5(cal.lifted_mask.clone());
    let pi_table: Vec<(usize, f32)> = LEVELS.iter().map(|&n| (n, 0.15 * corse[&(n / 2)].octaves[0])).collect();
    eprintln!("   π amplitudes (0.15 × Corsica's octave 0): {}", pi_table.iter().map(|(n, a)| format!("{n}² {a:.2} m")).collect::<Vec<_>>().join(" · "));
    // ── F161's N1, as committed (k from F161, the talus at F160's budget) ──
    let mut f161 = AmpConfig::declared(PSEED, Variant::N1, 1.0);
    f161.n1_k = vec![(512, 1.0e3), (1024, 1.0e3), (2048, 177.8)];
    eprintln!("\n   F161's N1 (as committed):");
    let (f161_lv, _, _) = run_chain("F161 N1", mk(f161), &pr, &corse, false);
    // ── the variants ──
    let mut runs: Vec<(String, Vec<Lvl>, f64, Vec<String>)> = Vec::new();
    let variants: Vec<(&str, Box<dyn Fn() -> AmpConfig>)> = vec![
        ("N1b", Box::new(|| base_cfg(Variant::N1))),
        ("ρ", Box::new(|| AmpConfig { n1_rho: Some(0.4), ..base_cfg(Variant::N1) })),
        ("π", Box::new(|| AmpConfig { n1_pi_m: pi_table.clone(), ..base_cfg(Variant::N1) })),
        ("ρ+π", Box::new(|| AmpConfig { n1_rho: Some(0.4), n1_pi_m: pi_table.clone(), ..base_cfg(Variant::N1) })),
        ("R4", Box::new(|| AmpConfig { n1_recal_depth: 2, ..base_cfg(Variant::N1) })),
    ];
    for (name, mkcfg) in &variants {
        eprintln!("\n   chain {name}:");
        let (lv, secs, stops) = run_chain(name, mk(mkcfg()), &pr, &corse, true);
        runs.push((name.to_string(), lv, secs, stops));
    }
    // the best predictability variant among ρ, π, ρ+π (fewest validated instruments outside ×1.5, then Σ |log ratio|)
    let cp195 = &corse_p[&1024];
    let score = |lv: &[Lvl]| -> (usize, f32) {
        let p = lv[2].p.as_ref().unwrap();
        let mut out = 0;
        let mut s = 0f32;
        for (i, ((_, v, _), (_, c, _))) in p.values().iter().zip(cp195.values().iter()).enumerate() {
            if !valid[i] {
                continue;
            }
            if out15(*v, *c) {
                out += 1;
            }
            if v.is_finite() && c.is_finite() && *c > 0.0 && *v > 0.0 {
                s += (v / c).ln().abs();
            }
        }
        (out, s)
    };
    let mut cands: Vec<(usize, f32, String)> = runs.iter().filter(|r| ["ρ", "π", "ρ+π"].contains(&r.0.as_str())).map(|r| { let (o, s) = score(&r.1); (o, s, r.0.clone()) }).collect();
    cands.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.partial_cmp(&b.1).unwrap()));
    let base_score = score(&runs[0].1);
    eprintln!(
        "\n   the predictability ranking at 2048² (validated instruments outside ×1.5 of Corsica, Σ|ln ratio|): N1b {} / {:.2} · {}",
        base_score.0,
        base_score.1,
        cands.iter().map(|c| format!("{} {} / {:.2}", c.2, c.0, c.1)).collect::<Vec<_>>().join(" · ")
    );
    let best = cands[0].2.clone();
    if (cands[0].0, cands[0].1) >= base_score {
        eprintln!("   NOTE: no predictability variant beats N1b; R4 is combined with {best} as declared");
    }
    let best_cfg = |c: AmpConfig| -> AmpConfig {
        match best.as_str() {
            "ρ" => AmpConfig { n1_rho: Some(0.4), ..c },
            "π" => AmpConfig { n1_pi_m: pi_table.clone(), ..c },
            _ => AmpConfig { n1_rho: Some(0.4), n1_pi_m: pi_table.clone(), ..c },
        }
    };
    let combo = format!("R4+{best}");
    eprintln!("\n   chain {combo}:");
    let (lv, secs, stops) = run_chain(&combo, mk(best_cfg(AmpConfig { n1_recal_depth: 2, ..base_cfg(Variant::N1) })), &pr, &corse, true);
    runs.push((combo.clone(), lv, secs, stops));
    // ── the readings and R ──
    let co195 = &corse[&1024];
    let co391 = &corse[&512];
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    eprintln!("\n   the readings (at 1024² against Corsica 391 m, at 2048² against Corsica 195 m):");
    let mut all: Vec<(String, &Vec<Lvl>)> = vec![("F161 N1".to_string(), &f161_lv)];
    for r in &runs {
        all.push((r.0.clone(), &r.1));
    }
    for (name, lv) in &all {
        for l in lv.iter().filter(|l| l.n >= 1024) {
            let co = if l.n == 1024 { co391 } else { co195 };
            let r = &l.r;
            let cell_km = DOMAIN_KM / l.n as f32;
            eprintln!(
                "   ── {name} {}²: slope p50/p90 {:.3}/{:.3} ({:.3}/{:.3}) · facets {:.2} % ({:.2}) · walls 28° {:.3} % /{} ({:.3}) · λ {:.2} km ({:.2}) · R8 {:.4} ({:.4}) · coherent {:.2} % ({:.2}) · density {:.3} ({:.3}) · peak {:.0} · k {:.3e}",
                l.n, r.slope_p50, r.slope_p90, co.slope_p50, co.slope_p90, r.facets, co.facets, r.walls28, r.walls28_cells, co.walls28, r.lambda_km, co.lambda_km, r.r8, co.r8, r.coherent, co.coherent, r.density, co.density, r.peak, l.k
            );
            eprintln!(
                "      spectrum (m), cascade / Corse: {}",
                (0..r.octaves.len().min(co.octaves.len()).min(6))
                    .map(|j| { let (a, b) = octave_km(j, cell_km); format!("[{a:.2}–{b:.2} km] {:.1} / {:.1}", r.octaves[j], co.octaves[j]) })
                    .collect::<Vec<_>>()
                    .join(" · ")
            );
            if let Some(p) = &l.p {
                print_predict(&format!("{name} {}²", l.n), p);
            }
        }
    }
    eprintln!("\n   R at 2048² (facet control: pyramids {:.1} % · paraboloids {:.1} %):", ctrl.0, ctrl.1);
    let mut verdicts = Vec::new();
    for (name, lv) in &all {
        let l = &lv[2];
        // the π share: Σ_L A² over the final field's variance in 0.4–12.5 km
        let cell_km = DOMAIN_KM / l.n as f32;
        let band_var: f64 = (0..l.r.octaves.len())
            .filter(|&j| { let (a, b) = octave_km(j, cell_km); a >= 0.39 && b <= 12.6 })
            .map(|j| (l.r.octaves[j] as f64).powi(2))
            .sum();
        let pi_sum: f64 = lv.iter().map(|x| x.pi_var).sum();
        let share = pi_sum / band_var.max(1e-9);
        let why = reject_2048(l, co195, cp195, &valid, ctrl_ok, share);
        eprintln!(
            "   {name:<10} π share {:.2} % · {}",
            100.0 * share,
            if why.is_empty() { "PASSES R at 2048²".to_string() } else { format!("REJECTED: {}", why.join("; ")) }
        );
        verdicts.push((name.clone(), why.len()));
    }
    for r in &runs {
        if !r.3.is_empty() {
            eprintln!("   per-level R, {}: {}", r.0, r.3.join(" | "));
        }
    }
    // the R4 deficit: the 1.5–6 km octaves at 2048² against N1b and Corsica
    let oct_km = |lv: &[Lvl], lo: f32| -> f32 {
        let l = &lv[2];
        let cell_km = DOMAIN_KM / l.n as f32;
        (0..l.r.octaves.len()).find(|&j| (octave_km(j, cell_km).0 - lo).abs() < 0.05).map_or(f32::NAN, |j| l.r.octaves[j])
    };
    let cor_oct = |lo: f32| -> f32 { (0..co195.octaves.len()).find(|&j| (octave_km(j, 0.1953125).0 - lo).abs() < 0.05).map_or(f32::NAN, |j| co195.octaves[j]) };
    for lo in [1.5625f32, 3.125] {
        let (b, c) = (oct_km(&runs[0].1, lo), cor_oct(lo));
        eprintln!(
            "   the {:.2}–{:.2} km octave at 2048²: N1b {:.1} · R4 {:.1} · {} {:.1} · Corsica {:.1} m → R4 closes {:.0} % of N1b's deficit",
            lo,
            2.0 * lo,
            b,
            oct_km(&runs[4].1, lo),
            combo,
            oct_km(&runs[5].1, lo),
            c,
            100.0 * (oct_km(&runs[4].1, lo) - b) / (c - b)
        );
    }
    eprintln!("\n   cost per chain (calibration included): {}", runs.iter().map(|r| format!("{} {:.0} s", r.0, r.2)).collect::<Vec<_>>().join(" · "));
    let l = &runs[0].1[2];
    eprintln!(
        "   N1b at 2048² without the calibration: {:.1} s (E {:.1} · T {:.1} · D {:.1}; {} / {} / {} iterations)",
        l.secs.iter().sum::<f64>(),
        l.secs[1],
        l.secs[2],
        l.secs[3],
        l.iters[0],
        l.iters[1],
        l.iters[2]
    );
    // ── I: the images ──
    let pick = |name: &str| -> &Vec<Lvl> { &all.iter().find(|a| a.0 == name).unwrap().1 };
    let shown: Vec<String> = vec!["N1b".into(), best.clone(), "R4".into(), combo.clone()];
    for n in [1024usize, 2048] {
        let i = if n == 1024 { 1 } else { 2 };
        let cell_m = DOMAIN_KM / n as f32 * 1000.0;
        let mut row = Vec::new();
        for s in &shown {
            let g = &pick(s)[i].result;
            let (_, hy) = read(g, None);
            row.push(whole(&rgba_of(g, Some(&hy), cell_m, &ss), n, 768));
        }
        let cg = &corse_g[&(n / 2)];
        let (_, chy) = read_cell(cg, None, 200.0 / (n / 2) as f32);
        row.push(whole(&rgba_of(cg, Some(&chy), cell_m, &ss), n / 2, 384));
        compose(&[row]).save(out_dir().join(format!("f162_whole_{n}.png"))).unwrap();
    }
    eprintln!("   images f162_whole_<n>.png: {} then Corsica (same km per pixel) · {NOTICE}", shown.join(", "));
    for (name, ox, oy) in CROPS.iter().copied().chain(std::iter::once(("tight_25km", TIGHT.0, TIGHT.1))) {
        let frac = if name == "tight_25km" { 0.0625 } else { 0.25 };
        let cfrac = frac * 2.0;
        let corigin = if name == "tight_25km" { CORSE_FLANK } else { CORSE_SPINE };
        let mut rows = Vec::new();
        for n in [1024usize, 2048] {
            let i = if n == 1024 { 1 } else { 2 };
            let cell_m = DOMAIN_KM / n as f32 * 1000.0;
            let mut row = Vec::new();
            for s in &shown {
                let g = &pick(s)[i].result;
                let (_, hy) = read(g, None);
                row.push(window(&rgba_of(g, Some(&hy), cell_m, &ss), n, (ox as f32 / 512.0, oy as f32 / 512.0), frac, 384));
            }
            let cg = &corse_g[&(n / 2)];
            let (_, chy) = read_cell(cg, None, 200.0 / (n / 2) as f32);
            row.push(window(&rgba_of(cg, Some(&chy), cell_m, &ss), n / 2, corigin, cfrac, 384));
            rows.push(row);
        }
        compose(&rows).save(out_dir().join(format!("f162_{name}.png"))).unwrap();
    }
    eprintln!(
        "   images f162_<crop>.png: rows 1024² and 2048², columns {} then Corsica (spine, or its western flank for the 25 km crop) · {NOTICE}",
        shown.join(", ")
    );
    let mut series: Vec<(String, Vec<(f32, f32)>, [u8; 3])> = Vec::new();
    let cols = [[200, 60, 40], [220, 140, 30], [60, 150, 60], [40, 110, 200]];
    for (i, s) in shown.iter().enumerate() {
        let l = &pick(s)[2];
        let ck = DOMAIN_KM / 2048.0;
        series.push((s.clone(), l.r.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, ck); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect(), cols[i]));
    }
    series.push(("Corse 195 m".into(), co195.octaves.iter().enumerate().map(|(j, &v)| { let (a, b) = octave_km(j, 0.1953125); ((a * b).sqrt(), v) }).filter(|p| p.0 <= 30.0).collect(), [0, 0, 0]));
    spectra_png(&series, &out_dir().join("f162_spectra.png"));
    eprintln!("   image f162_spectra.png: the 2048² octaves, {} in red / orange / green / blue, Corsica (195 m) black · {NOTICE}", shown.join(", "));
    let _ = (&verdicts, sorted(vec![0.0f32]), ocean_mask(&pm).len(), octave_rms(&[0.0; 16], 4, 1));
    eprintln!("\n==========  end Finding 162 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
