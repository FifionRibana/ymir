//! ADR Finding 167 — the thickness → altitude mapping: physical isostasy (Airy, Whitehead & Clift 2009) and an absolute
//! sea level, in the cascade only; the frozen cascade to 2 048²; the regularity baseline (for F169), Hovius's ratio.
//! Declared in `docs/reports/relief_method/f167_isostasy/f167_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f167_cascade -- --ignored --exact f167_cascade --nocapture
//!      cargo test -p ymir-core --release --test f167_cascade -- --ignored --exact f167_inventory --nocapture

#![allow(dead_code)]

mod common;


use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, Knobs, PSEED, bench_upscale_cfg, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, ocean_mask};
use ymir_core::cascade::amplify::restrict;
use ymir_core::cascade::history::{
    Airy, C1Sources, C1Term, HistoryRun, Mapping, T_F166_YEARS, Timing, isostatic_datum, land_altitude_m, physics_history,
    run_c1_sources, run_c1_sources_mapped,
};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{CLASS_NAMES, Hovius, MOUNTAIN, MacroStats, classify, coast, hovius, macro_stats};
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, flat_map, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::{CascadeConfig, Rebound, physics_level, roll};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};
use ymir_core::tectonics_c1::time_loop::run_with_closures;
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const SEEDS: [u64; 4] = [PSEED, 42, 1, 9];
const CELL: f32 = 1.5625;
const N_PHYS: usize = 300;
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];

fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f167_isostasy/images")
}
fn n2m() -> f32 {
    2.0 * 1.13 * SteinSteinParams::default().depth_scale_m as f32
}
fn out15(a: f32, b: f32) -> bool {
    !(a.is_finite() && b.is_finite() && b != 0.0) || a > 1.5 * b || a < b / 1.5
}
fn up_u8(m: &[u8], n: usize, f: usize) -> Vec<u8> {
    let nn = n * f;
    (0..nn * nn).map(|k| m[(k / nn / f) * n + (k % nn) / f]).collect()
}
fn roll_mask(m: &[bool], n: usize, off: [usize; 2]) -> Vec<bool> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n]).collect()
}
fn roll_f(m: &[f64], n: usize, off: [usize; 2]) -> Vec<f32> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n] as f32).collect()
}
fn peak_m(z: &GridF32) -> f32 {
    let oc = ocean_mask(z);
    z.data.iter().zip(&oc).filter(|(_, o)| !**o).map(|(v, _)| *v).fold(f32::MIN, f32::max)
}
/// The process's total CPU time (all threads), seconds.
fn cpu_secs() -> f64 {
    let o = std::process::Command::new("powershell")
        .args(["-NoProfile", "-Command", &format!("(Get-Process -Id {}).TotalProcessorTime.TotalSeconds", std::process::id())])
        .output();
    o.ok().and_then(|o| String::from_utf8_lossy(&o.stdout).trim().replace(',', ".").parse().ok()).unwrap_or(f64::NAN)
}
struct Clock {
    wall: Instant,
    cpu: f64,
}
impl Clock {
    fn start() -> Self {
        Self { cpu: cpu_secs(), wall: Instant::now() }
    }
    fn stop(&self) -> (f64, f64) {
        let w = self.wall.elapsed().as_secs_f64();
        let c = cpu_secs() - self.cpu;
        (c, w)
    }
}
fn fmt_cost(c: (f64, f64)) -> String {
    format!("{:.1} s CPU / {:.1} s wall{}", c.0, c.1, if c.1 > 5.0 * c.0.max(0.1) { " SUSPENDED?" } else { "" })
}
fn class_png(classes: &[u8], w: usize, h: usize, scale: u32) -> image::RgbaImage {
    let col = |c: u8| match c {
        0 => [70u8, 110, 170],
        1 => [170, 210, 120],
        2 => [220, 190, 120],
        3 => [150, 160, 70],
        _ => [130, 80, 50],
    };
    let mut img = image::RgbaImage::new(w as u32, h as u32);
    for y in 0..h {
        for x in 0..w {
            let c = col(classes[y * w + x]);
            img.put_pixel(x as u32, (h - 1 - y) as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    }
    image::imageops::resize(&img, w as u32 * scale, h as u32 * scale, image::imageops::FilterType::Nearest)
}
fn print_macro(tag: &str, m: &MacroStats) {
    eprintln!(
        "   {tag:<24} plain {:5.1} % · plateau {:5.1} % · hill {:5.1} % · mountain {:5.1} % · land {:.0} km² · deciles {}",
        100.0 * m.fractions[0],
        100.0 * m.fractions[1],
        100.0 * m.fractions[2],
        100.0 * m.fractions[3],
        m.land_km2,
        m.deciles.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join("/")
    );
}

fn c1_cfg() -> C1TimeLoopConfig {
    C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    }
}
fn c1_run(seed: u64, exclude: &[C1Term]) -> (C1Sources, ymir_core::tectonics_c1::state::C1State, PlateKinematics) {
    let cfg = c1_cfg();
    let mut st = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let s = run_c1_sources(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &SteinSteinParams::default(), exclude);
    (s, st, kin)
}

struct Tex {
    r: Reading,
    p: Predict,
    pf: Planform,
}
fn tex(g: &GridF32, cell_km: f32, excl: &[bool]) -> Tex {
    let (r, _) = read_cell(g, Some(excl), cell_km);
    Tex { r, p: predictability(g, cell_km, Some(excl)), pf: planform(g, cell_km, Some(excl)) }
}
fn elements(a: &Tex, c: &Tex, cell_km: f32, ctrl_ok: bool) -> Vec<(String, f32, f32, bool)> {
    let mut v = Vec::new();
    v.push(("facets".into(), a.r.facets as f32, c.r.facets as f32, !(ctrl_ok && a.r.facets > 1.5 * c.r.facets && a.r.facet_cells > 20)));
    v.push(("walls 28°".into(), a.r.walls28 as f32, c.r.walls28 as f32, !(a.r.walls28 > 1.5 * c.r.walls28 && a.r.walls28_cells > 20)));
    for (n, x, y) in [("slope p50", a.r.slope_p50, c.r.slope_p50), ("slope p90", a.r.slope_p90, c.r.slope_p90), ("λ km", a.r.lambda_km, c.r.lambda_km)] {
        v.push((n.to_string(), x, y, !out15(x, y)));
    }
    for j in 0..a.r.octaves.len().min(c.r.octaves.len()) {
        let (lo, hi) = octave_km(j, cell_km);
        if lo >= 0.39 && hi <= 12.6 {
            v.push((format!("octave {lo:.2}–{hi:.2} km"), a.r.octaves[j], c.r.octaves[j], !out15(a.r.octaves[j], c.r.octaves[j])));
        }
    }
    for (i, ((n, x, _), (_, y, _))) in a.p.values().iter().zip(c.p.values().iter()).enumerate() {
        if VALID_P[i] {
            v.push((n.to_string(), *x, *y, !out15(*x, *y)));
        }
    }
    for i in 0..3 {
        v.push((format!("S_{} − 1", WINDOWS_KM[i]), a.pf.excess(i), c.pf.excess(i), !out15(a.pf.excess(i), c.pf.excess(i))));
    }
    v.push(("A_dir".into(), a.pf.a_dir, c.pf.a_dir, !out15(a.pf.a_dir, c.pf.a_dir)));
    v.push(("F".into(), a.pf.flat, c.pf.flat, !out15(a.pf.flat, c.pf.flat)));
    v
}

struct ChainOut {
    g: GridF32,
    d5: Vec<bool>,
    cost: (f64, f64),
    notes: Vec<String>,
    /// (level n, footprints of the new depressions at that level)
    fresh: Vec<(usize, Vec<Vec<usize>>)>,
}
fn new_depressions(a: &GridF32, b: &GridF32, cell_km: f32) -> (usize, usize, Vec<Vec<usize>>) {
    let (ha, hb) = (hydrology(a, cell_km), hydrology(b, cell_km));
    let mut fresh = vec![true; hb.lakes.len()];
    for k in 0..hb.lake_id.len() {
        if hb.lake_id[k] != 0 && ha.lake_id[k] != 0 {
            fresh[hb.lake_id[k] as usize - 1] = false;
        }
    }
    let mut prints: Vec<Vec<usize>> = vec![Vec::new(); hb.lakes.len()];
    for k in 0..hb.lake_id.len() {
        if hb.lake_id[k] != 0 && fresh[hb.lake_id[k] as usize - 1] {
            prints[hb.lake_id[k] as usize - 1].push(k);
        }
    }
    (ha.lakes.len(), hb.lakes.len(), prints.into_iter().filter(|p| !p.is_empty()).collect())
}
fn chain(pm: &GridF32, d5: &[bool], cfg: AmpConfig) -> ChainOut {
    let clk = Clock::start();
    let mut ch = Chain::new(pm.clone(), cfg).with_d5(d5.to_vec());
    let mut notes = Vec::new();
    let mut fresh = Vec::new();
    for _ in 0..3 {
        let l = ch.next_level(&|| false).unwrap().clone();
        let n = l.n;
        if let Some((a, b)) = l.pre_warp.as_ref().zip(l.post_warp.as_ref()) {
            let (na, nb, pr) = new_depressions(a, b, DOMAIN_KM / n as f32);
            notes.push(format!("{n}²: depressions {na} → {nb}, new {}", pr.len()));
            fresh.push((n, pr));
        }
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let last = ch.levels.last().unwrap();
    ChainOut { g: last.result(), d5: last.d5.clone(), cost: clk.stop(), notes, fresh }
}

struct Run {
    name: String,
    zm: GridF32,
    d5: Vec<bool>,
    classes: Vec<u8>,
    m: MacroStats,
    cost: (f64, f64),
    sum_up: Option<Vec<f32>>,
}
fn run_of(name: &str, zm: GridF32, d5: Vec<bool>, cost: (f64, f64), sum_up: Option<Vec<f32>>) -> Run {
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    Run { name: name.into(), zm, d5, classes, m, cost, sum_up }
}
fn hist_run(name: &str, snaps: &[GridF32], cfg: &CascadeConfig, t: f64) -> Run {
    let clk = Clock::start();
    let l = physics_history(snaps, 256, cfg, &HistoryRun { t_years: t, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History }, &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask(&zm);
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    run_of(name, zm, d5, cost, Some(l.rec.sum_uplift.iter().map(|v| v * k).collect()))
}


fn c1_run_m(seed: u64, exclude: &[C1Term], mapping: Mapping) -> (C1Sources, ymir_core::tectonics_c1::state::C1State) {
    let cfg = c1_cfg();
    let mut st = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let s = run_c1_sources_mapped(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &SteinSteinParams::default(), exclude, mapping);
    (s, st)
}
fn up_mask(m: &[bool], n: usize, f: usize) -> Vec<bool> {
    let nn = n * f;
    (0..nn * nn).map(|k| m[(k / nn / f) * n + (k % nn) / f]).collect()
}
fn land_of(z: &GridF32) -> Vec<bool> {
    ocean_mask(z).iter().map(|o| !o).collect()
}

#[test]
#[ignore]
fn f167_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let iso = IsostasyConfig::c1_default();
    let k = n2m();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 167 . the thickness → altitude mapping (Airy, absolute sea level), cascade only  ==========");
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
    let tlf = bench_upscale_cfg(temoin).target_land_fraction;
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let air = Airy::declared();
    eprintln!(
        "   the Airy mapping: S̃ = 0.2 → {:.2} km · 0.6 → {:.2} · 1.0 → {:.2} · 1.25 → {:.2} · 1.5 → {:.2} · 2.0 → {:.2} km (the normalised ceiling 5.65 km) · the sea at S̃ = {:.3} · a 5.5 km crust (L&M's ridge column) → {:.2} km",
        air.altitude_km(0.2), air.altitude_km(0.6), air.altitude_km(1.0), air.altitude_km(1.25), air.altitude_km(1.5), air.altitude_km(2.0),
        { let bl = (air.rho_m - air.rho_c) / air.rho_m * air.km_per_s; (bl - air.anchor_km) / bl },
        air.altitude_km(5.5 / 35.0)
    );
    // Europe, Corsica
    let e = load_grid(&root().join("data/europe/europe_1920.bin")).expect("data/europe");
    let me = macro_stats(&e.data, &classify(&e.data, e.width, e.height, CELL), e.width, e.height, CELL);
    print_macro("Europe", &me);
    let cg = load_grid(&root().join("data/corsica/corse_1024.bin")).unwrap();
    let c512 = load_grid(&root().join("data/corsica/corse_512.bin")).unwrap();
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let cclass = up_u8(&classify(&c128.data, 128, 128, CELL), 128, 8);
    let cmtn: Vec<bool> = cclass.iter().map(|&c| c == MOUNTAIN).collect();
    let ck = 200.0 / 1024.0;
    let c_whole = tex(&cg, ck, &vec![false; 1024 * 1024]);
    let c_mtn = tex(&cg, ck, &cmtn.iter().map(|m| !m).collect::<Vec<_>>());
    let c_coast = coast(&land_of(&cg), 1024, 1024, 32);
    let cp391 = predictability(&c512, 200.0 / 512.0, None);
    let hov_c: Vec<(f32, Hovius)> = [10.0f32, 25.0, 100.0].iter().map(|&a| (a, hovius(&cg, &cmtn, ck, a))).collect();
    eprintln!("   Corse 195 m: mountain-class slope p50/p90 {:.3}/{:.3} · coast {:?}", c_mtn.r.slope_p50, c_mtn.r.slope_p90, c_coast);
    for (a, h) in &hov_c {
        eprintln!("   Corse Hovius (A_min {a} km²): R {:.2} over {} basins · S {:.2} km · W {:.2} km", h.ratio, h.outlets, h.spacing_km, h.half_width_km);
    }
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let frozen = || AmpConfig {
        talus_fine_scale: 0.25,
        n1_recal_depth: 2,
        n1_pi_m: pi_table.clone(),
        warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }),
        n1_k: vec![(512, 5.623e2), (1024, 1.0e3), (2048, 3.162e2)],
        ..AmpConfig::declared(PSEED, Variant::N1, 1.0)
    };
    let base = |seed: u64| CascadeConfig { seed, roughness_m: rough_a, mfd_exponent: Some(6.0), ..CascadeConfig::declared(ss.depth_scale_m as f32) };
    let snaps = |c: &C1Sources| -> Vec<GridF32> { c.hist.normalized(tlf).iter().map(|g| roll(g, off)).collect() };
    let ero = [C1Term::Erosion];
    let ero_eh = [C1Term::Erosion, C1Term::EquilibriumHeight];
    let own = Airy { oceanic_s: None, ..air };
    let lo = Airy { anchor_km: air.anchor_km - 0.2, ..air };
    let hi = Airy { anchor_km: air.anchor_km + 0.2, ..air };
    let names = ["F164 (steady)", "F166 (C1 mapping)", "F167 (Airy)", "F167 − EH", "F167 own s", "F167 own s − EH", "anchor 0.25 km", "anchor 0.65 km"];
    let mut per_seed: Vec<(u64, Vec<Run>, ymir_core::tectonics_c1::state::C1State, C1Sources)> = Vec::new();
    let mut map_cost = Vec::new();
    for seed in SEEDS {
        let cfg = base(seed);
        let (c1m, _) = c1_run_m(seed, &ero, Mapping::C1);
        let clk = Clock::start();
        let (am, st) = c1_run_m(seed, &ero, Mapping::Airy(air));
        map_cost.push(clk.stop());
        let (aeh, _) = c1_run_m(seed, &ero_eh, Mapping::Airy(air));
        let (ao, _) = c1_run_m(seed, &ero, Mapping::Airy(own));
        let (aoeh, _) = c1_run_m(seed, &ero_eh, Mapping::Airy(own));
        let (alo, _) = c1_run_m(seed, &ero, Mapping::Airy(lo));
        let (ahi, _) = c1_run_m(seed, &ero, Mapping::Airy(hi));
        let (all, _) = c1_run_m(seed, &[], Mapping::C1);
        let clk = Clock::start();
        let (rec, cal) = physics_level(snaps(&all).last().unwrap(), 256, &cfg, 2850.0, 300, &mut |_| {}, &|| false).unwrap();
        let f164 = run_of(names[0], to_m(&rec.z, k), cal.lifted_mask.clone(), clk.stop(), None);
        let mut runs = vec![f164];
        for (i, c) in [&c1m, &am, &aeh, &ao, &aoeh, &alo, &ahi].iter().enumerate() {
            runs.push(hist_run(names[i + 1], &snaps(c), &cfg, T_F166_YEARS));
        }
        per_seed.push((seed, runs, st, am));
    }
    eprintln!("   the Airy record's cost per C1 run: {}", map_cost.iter().map(|c| fmt_cost(*c)).collect::<Vec<_>>().join(" · "));
    // ── 256² ──
    eprintln!("\n   C — 256² (Europe: plain 51.6 · plateau 4.7 · hill 18.9 · mountain 24.8):");
    for (seed, runs, st, _) in &per_seed {
        eprintln!("   ── seed {seed}:");
        let crat = up_mask(&roll_mask(st.cratonic_mask.data(), 64, off), 64, 4);
        let ocean_plate: Vec<bool> = up_mask(&roll_mask(&st.plate_type.data().iter().map(|t| *t == PlateType::Oceanic).collect::<Vec<_>>(), 64, off), 64, 4);
        let l164 = land_of(&runs[0].zm);
        let n164 = l164.iter().filter(|&&b| b).count() as f32;
        for r in runs {
            print_macro(&r.name, &r.m);
            let land = land_of(&r.zm);
            let lost = (0..65536).filter(|&i| l164[i] && !land[i]).count() as f32 / n164;
            let gained = (0..65536).filter(|&i| land[i] && !l164[i]).count() as f32 / n164;
            let nl = land.iter().filter(|&&b| b).count();
            let le0 = (0..65536).filter(|&i| land[i] && r.zm.data[i] <= 0.0).count() as f32 / nl.max(1) as f32;
            let cr = sorted((0..65536).filter(|&i| crat[i] && land[i]).map(|i| r.zm.data[i]).collect());
            let phantom = (0..65536).filter(|&i| ocean_plate[i] && r.zm.data[i] > 0.0).count();
            let hy = hydrology(&r.zm, CELL);
            let lake_cells = hy.lake_id.iter().filter(|&&l| l != 0).count();
            eprintln!(
                "      {:<20} peak {:.0} m · cratons' land p50 {:.0} m ({} cells) · land lost {:.1} % · gained {:.1} % vs F164 · interior ≤ 0 m {:.2} % · land on oceanic plate {} cells · lakes {} / {:.0} km² · {}",
                r.name,
                peak_m(&r.zm),
                pct(&cr, 0.5),
                cr.len(),
                100.0 * lost,
                100.0 * gained,
                100.0 * le0,
                phantom,
                hy.lakes.len(),
                lake_cells as f32 * CELL * CELL,
                fmt_cost(r.cost)
            );
        }
        // the decoupling tests
        for (a, b, tag) in [(2usize, 3usize, "main"), (4, 5, "own s")] {
            let (la, lb) = (land_of(&runs[a].zm), land_of(&runs[b].zm));
            let diff = (0..65536).filter(|&i| la[i] != lb[i]).count();
            let na = la.iter().filter(|&&x| x).count();
            eprintln!("      decoupling ({tag}): sources vs sources − EH differ on {diff} land cells = {:.3} % of the land → {}", 100.0 * diff as f32 / na.max(1) as f32, if (diff as f32) <= 0.001 * na as f32 { "IDENTICAL within 0.1 %" } else { "NOT within 0.1 %" });
        }
    }
    // ── the frozen chains: F164, F166, F167 ──
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    let ckm = DOMAIN_KM / 2048.0;
    let mut finals: Vec<(u64, Vec<ChainOut>)> = Vec::new();
    for (seed, runs, _, _) in &per_seed {
        let cs: Vec<ChainOut> = (0..3).map(|i| chain(&runs[i].zm, &runs[i].d5, frozen())).collect();
        finals.push((*seed, cs));
    }
    let mut successes = 0;
    let mut mtn_over = 0;
    let mut lost_over = 0;
    let mut reg_rows: Vec<String> = Vec::new();
    for ((seed, cs), (_, runs, st, am)) in finals.iter().zip(&per_seed) {
        eprintln!("\n   ── seed {seed} at 2048² (195 m) against Corsica:");
        let k166 = coast(&land_of(&cs[1].g), 2048, 2048, 32);
        for (ci, c) in cs.iter().enumerate() {
            let r = &runs[ci];
            let m8 = up_u8(&r.classes, 256, 8);
            let mtn8: Vec<bool> = m8.iter().map(|&x| x == MOUNTAIN).collect();
            let excl: Vec<bool> = (0..2048 * 2048).map(|i| c.d5[i] || !mtn8[i]).collect();
            let whole = tex(&c.g, ckm, &c.d5);
            let mtn = tex(&c.g, ckm, &excl);
            for (tag, x, cc) in [("whole", &whole, &c_whole), ("mountain", &mtn, &c_mtn)] {
                let el = elements(x, cc, ckm, ctrl_ok);
                let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
                eprintln!("      {:<18} {tag:<9} {} of {} within ×1.5 · out: {}", r.name, el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") });
            }
            eprintln!("      {:<18} slope guard: peak {:.0} m · mountain-class slope p50/p90 {:.3}/{:.3} (Corsica {:.3}/{:.3})", r.name, peak_m(&c.g), mtn.r.slope_p50, mtn.r.slope_p90, c_mtn.r.slope_p50, c_mtn.r.slope_p90);
            let (rd, fl) = flat_map(&c.g, ckm, Some(&c.d5));
            let per: Vec<String> = (1..5u8)
                .map(|cl| {
                    let cells = (0..2048 * 2048).filter(|&i| rd[i] && m8[i] == cl).count();
                    let f = (0..2048 * 2048).filter(|&i| rd[i] && fl[i] && m8[i] == cl).count();
                    format!("{} {:.2} %", CLASS_NAMES[cl as usize], 100.0 * f as f32 / cells.max(1) as f32)
                })
                .collect();
            eprintln!("      {:<18} F by class: {}", r.name, per.join(" · "));
            let kc = coast(&land_of(&c.g), 2048, 2048, 32);
            eprintln!("      {:<18} coast {kc:?}", r.name);
            // the regularity baseline
            let p2048 = &whole.p;
            let p1024 = {
                let g = restrict(&c.g);
                let d = (0..1024 * 1024).map(|i| c.d5[(i / 1024 * 2) * 2048 + (i % 1024) * 2]).collect::<Vec<_>>();
                predictability(&g, 2.0 * ckm, Some(&d))
            };
            let hv: Vec<String> = [10.0f32, 25.0, 100.0]
                .iter()
                .map(|&a| {
                    let h = hovius(&c.g, &mtn8, ckm, a);
                    format!("A ≥ {a}: R {:.2} ({} basins, S {:.1} km, W {:.1} km)", h.ratio, h.outlets, h.spacing_km, h.half_width_km)
                })
                .collect();
            reg_rows.push(format!(
                "   seed {seed} {:<18} 391 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | 195 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | Hovius {}",
                r.name, p1024.cv_lambda, p1024.cv_len, p1024.theta_mean, p1024.theta_std, p1024.r2_gully, p2048.cv_lambda, p2048.cv_len, p2048.theta_mean, p2048.theta_std, p2048.r2_gully, hv.join(" · ")
            ));
            if ci == 2 {
                let degraded = |a: f32, b: f32, cc: f32| (b - cc).abs() > (a - cc).abs() + 0.10 * a.abs();
                let deg: Vec<&str> = [
                    ("dimension", degraded(k166.dimension, kc.dimension, c_coast.dimension)),
                    ("length ratio", degraded(k166.length_ratio, kc.length_ratio, c_coast.length_ratio)),
                    ("aligned", degraded(k166.aligned, kc.aligned, c_coast.aligned)),
                    ("block index", degraded(k166.block_index, kc.block_index, c_coast.block_index)),
                ]
                .iter()
                .filter(|x| x.1)
                .map(|x| x.0)
                .collect();
                let l164 = land_of(&runs[0].zm);
                let lr = land_of(&r.zm);
                let lost = (0..65536).filter(|&i| l164[i] && !lr[i]).count() as f32 / l164.iter().filter(|&&b| b).count() as f32;
                let m = &r.m;
                let mtn_ok = (0.165..=0.372).contains(&m.fractions[3]);
                let low = m.fractions[0] + m.fractions[1];
                let low_ok = (0.375..=0.845).contains(&low);
                let succ = mtn_ok && low_ok && lost < 0.15 && deg.is_empty();
                successes += succ as usize;
                mtn_over += (m.fractions[3] > 0.372) as usize;
                lost_over += (lost > 0.15) as usize;
                eprintln!(
                    "      SUCCESS test (F167): mountain {:.1} % {mtn_ok} · plain + plateau {:.1} % {low_ok} · land lost {:.1} % {} · coast degraded vs F166 {} → {}",
                    100.0 * m.fractions[3],
                    100.0 * low,
                    100.0 * lost,
                    lost < 0.15,
                    if deg.is_empty() { "none".into() } else { deg.join(", ") },
                    if succ { "SUCCESS" } else { "no" }
                );
                // the mountain decomposition (256²): cratons, Davis-Suppe zones, arcs, the rest
                let term = |t: C1Term| -> Vec<bool> {
                    let ds = &am.terms.iter().find(|x| x.0 == t).unwrap().1;
                    up_mask(&roll_mask(&ds.iter().map(|&v| v > 0.05).collect::<Vec<_>>(), 64, off), 64, 4)
                };
                let (dsz, arc) = (term(C1Term::DavisSuppe), term(C1Term::Subduction));
                let crat = up_mask(&roll_mask(st.cratonic_mask.data(), 64, off), 64, 4);
                let mc: Vec<usize> = (0..65536).filter(|&i| r.classes[i] == MOUNTAIN).collect();
                let nm = mc.len().max(1) as f32;
                let c_cr = mc.iter().filter(|&&i| crat[i]).count() as f32 / nm;
                let c_ds = mc.iter().filter(|&&i| !crat[i] && dsz[i]).count() as f32 / nm;
                let c_arc = mc.iter().filter(|&&i| !crat[i] && !dsz[i] && arc[i]).count() as f32 / nm;
                eprintln!(
                    "      the mountain class (256²) decomposed: cratons {:.1} % · Davis-Suppe zones (Σ Δs > 0.05) {:.1} % · arcs {:.1} % · the rest {:.1} %",
                    100.0 * c_cr,
                    100.0 * c_ds,
                    100.0 * c_arc,
                    100.0 * (1.0 - c_cr - c_ds - c_arc)
                );
            }
        }
    }
    eprintln!("\n   the regularity baseline (Corsica 391 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | 195 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3}):", cp391.cv_lambda, cp391.cv_len, cp391.theta_mean, cp391.theta_std, cp391.r2_gully, c_whole.p.cv_lambda, c_whole.p.cv_len, c_whole.p.theta_mean, c_whole.p.theta_std, c_whole.p.r2_gully);
    for r in &reg_rows {
        eprintln!("{r}");
    }
    eprintln!(
        "\n   R: macro success on {successes} / 4 · mountains > 37.2 % on {mtn_over} / 4 · land lost > 15 % on {lost_over} / 4 → {}",
        if successes >= 3 { "SUCCESS: F168 follows" } else if mtn_over >= 2 { "the mountains still exceed 37.2 %: decomposed above, stop" } else if lost_over >= 2 { "the land lost exceeds 15 %: the anchor and its ±0.2 km reported, stop" } else { "no success, no stop condition met" }
    );
    // ── images ──
    let sh = |z: &GridF32, n: usize, px: usize| {
        let (_, hy) = read(z, None);
        whole(&rgba_of(z, Some(&hy), DOMAIN_KM / n as f32 * 1000.0, &ss), n, px)
    };
    let (mut r256, mut r2048, mut rcls) = (Vec::new(), Vec::new(), Vec::new());
    for ((_, cs), (_, runs, _, _)) in finals.iter().zip(&per_seed) {
        r256.push(vec![sh(&runs[1].zm, 256, 320), sh(&runs[2].zm, 256, 320)]);
        r2048.push(vec![sh(&cs[1].g, 2048, 448), sh(&cs[2].g, 2048, 448)]);
        rcls.push(vec![class_png(&runs[1].classes, 256, 256, 1), class_png(&runs[2].classes, 256, 256, 1)]);
    }
    compose(&r256).save(out().join("f167_seeds_256.png")).unwrap();
    compose(&r2048).save(out().join("f167_seeds_2048.png")).unwrap();
    compose(&rcls).save(out().join("f167_classes_seeds.png")).unwrap();
    // the mapping plotted: h(s), C1's land ramp at the témoin's final datum against Airy
    let am0 = &per_seed[0].3;
    let hs1 = isostatic_datum(&am0.s_final, &am0.craton, &iso);
    let (w, h) = (800u32, 500u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let px = |s: f64| (40.0 + s / 2.5 * (w as f64 - 60.0)) as i64;
    let py = |km: f64| (h as f64 - 30.0 - (km + 6.0) / 13.0 * (h as f64 - 50.0)) as i64;
    let mut put = |x: i64, y: i64, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as u32) < w && (y as u32) < h {
            img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    };
    for x in 40..(w as i64 - 20) {
        put(x, py(0.0), [0, 0, 0]);
        put(x, py(5.65), [200, 200, 200]);
    }
    for s in [0.5f64, 1.0, 1.5, 2.0] {
        for y in 20..(h as i64 - 30) {
            put(px(s), y, [220, 220, 220]);
        }
    }
    for i in 0..1000 {
        let s = 2.5 * i as f64 / 1000.0;
        for d in -1..=1 {
            put(px(s), py(air.altitude_km(s)) + d, [40, 110, 200]);
            let c1 = land_altitude_m(s, false, hs1, &iso) / 1000.0;
            if c1 >= 0.0 {
                put(px(s), py(c1.min(5.65)) + d, [200, 60, 40]);
            }
        }
    }
    img.save(out().join("f167_mapping.png")).unwrap();
    // coast crops: Corsica | F166 | F167
    let pick = |g: &GridF32, n: usize| -> (usize, usize) {
        let land = land_of(g);
        for oy in (0..n - 128).step_by(64) {
            for ox in (0..n - 128).step_by(64) {
                let c = (0..128 * 128).filter(|&i| land[(oy + i / 128) * n + ox + i % 128]).count() as f32 / 16384.0;
                if (0.4..=0.6).contains(&c) {
                    return (ox, oy);
                }
            }
        }
        (0, 0)
    };
    let (g166, g167) = (&finals[0].1[1].g, &finals[0].1[2].g);
    let (ox, oy) = pick(g167, 2048);
    let (cx, cy) = pick(&cg, 1024);
    let crop = |g: &GridF32, n: usize, x: usize, y: usize| {
        let (_, hy) = read(g, None);
        window(&rgba_of(g, Some(&hy), 195.3, &ss), n, (x as f32 / n as f32, y as f32 / n as f32), 128.0 / n as f32, 448)
    };
    compose(&[vec![crop(&cg, 1024, cx, cy), crop(g166, 2048, ox, oy), crop(g167, 2048, ox, oy)]]).save(out().join("f167_coast_25km.png")).unwrap();
    eprintln!(
        "   images: f167_seeds_256.png, f167_seeds_2048.png, f167_classes_seeds.png (rows = seeds {SEEDS:?}; columns F166 | F167), f167_mapping.png (h(S̃), 0–2.5 on x, −6 to +7 km on y; C1's land ramp at the témoin's final datum red, Airy blue; sea level black, the 5.65 km ceiling grey), f167_coast_25km.png (Corsica | F166 | F167; window ({ox}, {oy}) of 2048², Corsica ({cx}, {cy})) · north up · {NOTICE}"
    );
    eprintln!("\n==========  end Finding 167 . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F167-A1, the inventory: the distribution of s on oceanic and continental cells, at t₀ and at the end, per seed.
#[test]
#[ignore]
fn f167_inventory() {
    for seed in [PSEED, 42, 1, 9] {
        let cfg = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0 / 64.0, dy: 1.0 / 64.0, iso_config: IsostasyConfig::c1_default(), drainage_max_distance: 30 };
        let mut st = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let split = |st: &ymir_core::tectonics_c1::state::C1State| {
            let (mut o, mut c) = (Vec::new(), Vec::new());
            for (k, v) in st.s.data().iter().enumerate() {
                if st.plate_type.data()[k] == PlateType::Oceanic { o.push(*v as f32) } else { c.push(*v as f32) }
            }
            (sorted(o), sorted(c))
        };
        let (o0, c0) = split(&st);
        let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
        run_with_closures(&mut st, &mut kin, &cfg, &C1Closures::default(), |_, _| {});
        let (o1, c1) = split(&st);
        let q = |v: &[f32]| format!("p01 {:.3} p50 {:.3} p99 {:.3} max {:.3} · > 0.915: {} / {}", pct(v, 0.01), pct(v, 0.5), pct(v, 0.99), v.last().copied().unwrap_or(0.0), v.iter().filter(|&&x| x > 0.915).count(), v.len());
        eprintln!("   seed {seed}: t₀ oceanic s {} · continental s {}", q(&o0), q(&c0));
        eprintln!("   seed {seed}: end oceanic s {} · continental s {}", q(&o1), q(&c1));
    }
}
