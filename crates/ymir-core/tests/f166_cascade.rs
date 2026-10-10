//! ADR Finding 166 — the physics level driven by C1's tectonic sources only (the erosion done once, by the physics
//! level): A2 (the decomposition of F165's history by C1 term), the sources' drive with T calibrated once on the témoin,
//! the variant without the equilibrium-height sink, the frozen cascade to 2 048², the macro / texture / coast measures.
//! Declared in `docs/reports/relief_method/f166_sources/f166_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f166_cascade -- --ignored --exact f166_cascade --nocapture
//!      cargo test -p ymir-core --release --test f166_cascade -- --ignored --exact f166_terms_land --nocapture

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, Knobs, PSEED, bench_upscale_cfg, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, ocean_mask};
use ymir_core::cascade::history::{
    C1Sources, C1Term, HistoryRun, T_F165_YEARS, Timing, isostatic_datum, land_altitude_m, physics_history, run_c1_sources,
};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{CLASS_NAMES, MOUNTAIN, MacroStats, PLATEAU, chebyshev, classify, coast, components, macro_stats};
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, flat_map, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::{CascadeConfig, Rebound, diff_rgba, physics_level, roll};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};

const SEEDS: [u64; 4] = [PSEED, 42, 1, 9];
const CELL: f32 = 1.5625;
const N_PHYS: usize = 300;
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];

fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f166_sources/images")
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

#[test]
#[ignore]
fn f166_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let iso = IsostasyConfig::c1_default();
    let k = n2m();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 166 . the physics level driven by C1's tectonic sources only  ==========");
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
    let tlf = bench_upscale_cfg(temoin).target_land_fraction;
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    // Europe (F165, frozen)
    let e = load_grid(&root().join("data/europe/europe_1920.bin")).expect("data/europe");
    let ce = classify(&e.data, e.width, e.height, CELL);
    let me = macro_stats(&e.data, &ce, e.width, e.height, CELL);
    print_macro("Europe", &me);
    // Corsica
    let cg = load_grid(&root().join("data/corsica/corse_1024.bin")).unwrap();
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let cclass = up_u8(&classify(&c128.data, 128, 128, CELL), 128, 8);
    let ck = 200.0 / 1024.0;
    let c_whole = tex(&cg, ck, &vec![false; 1024 * 1024]);
    let c_mtn = tex(&cg, ck, &cclass.iter().map(|&c| c != MOUNTAIN).collect::<Vec<_>>());
    let c_coast = coast(&ocean_mask(&cg).iter().map(|o| !o).collect::<Vec<_>>(), 1024, 1024, 32);
    eprintln!("   Corse 195 m: mountain-class slope p50/p90 {:.3}/{:.3} · coast {:?}", c_mtn.r.slope_p50, c_mtn.r.slope_p90, c_coast);
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

    // ── A2: the decomposition of F165's history (the témoin), and A3 ──
    let clk = Clock::start();
    let (all, st, kin) = c1_run(PSEED, &[]);
    eprintln!("\n   A2 — F165's history decomposed by C1 term (the témoin; C1 with the observer: {})", fmt_cost(clk.stop()));
    let n = 64;
    let h_fin: Vec<f32> = all.hist.normalized(tlf).last().unwrap().data.iter().map(|v| (v - 0.5) * k).collect();
    let h_ini: Vec<f32> = all.hist.normalized(tlf)[0].data.iter().map(|v| (v - 0.5) * k).collect();
    let land: Vec<bool> = h_fin.iter().map(|&v| v > 0.0).collect();
    let bi = classify_boundaries(&st.plate_id, &kin);
    let src: Vec<bool> = (0..n * n).map(|i| !land[i] || bi.boundary_type.data()[i] == BoundaryType::Convergent || all.hist.suture[i]).collect();
    let dist = chebyshev(&src, n, n, 5);
    let interior: Vec<bool> = (0..n * n).map(|i| land[i] && dist[i] > 3).collect();
    let margin: Vec<bool> = (0..n * n).map(|i| land[i] && dist[i] <= 3).collect();
    let (hs0, hs1) = (isostatic_datum(&all.s0, &all.craton, &iso), isostatic_datum(&all.s_final, &all.craton, &iso));
    eprintln!("      the isostatic datum (raw h_sea): t₀ {hs0:.4} → final {hs1:.4} · land cells: interior {} · margin {}", interior.iter().filter(|&&b| b).count(), margin.iter().filter(|&&b| b).count());
    let mean_over = |v: &dyn Fn(usize) -> f64, m: &[bool]| -> f64 { (0..n * n).filter(|&i| m[i]).map(v).sum::<f64>() / m.iter().filter(|&&b| b).count().max(1) as f64 };
    let mut rows: Vec<(String, f64, f64, f64, f64)> = Vec::new();
    for (term, ds) in &all.terms {
        let dh = |i: usize| land_altitude_m(all.s_final[i], all.craton[i], hs1, &iso) - land_altitude_m(all.s_final[i] - ds[i], all.craton[i], hs1, &iso);
        rows.push((format!("{term:?}"), mean_over(&|i| ds[i], &interior), mean_over(&dh, &interior), mean_over(&|i| ds[i], &margin), mean_over(&dh, &margin)));
    }
    let datum = |i: usize| land_altitude_m(all.s_final[i], all.craton[i], hs1, &iso) - land_altitude_m(all.s_final[i], all.craton[i], hs0, &iso);
    rows.push(("datum drift".into(), 0.0, mean_over(&datum, &interior), 0.0, mean_over(&datum, &margin)));
    let obs = |i: usize| (h_fin[i] - h_ini[i]) as f64;
    for r in &rows {
        eprintln!("      {:<20} interior Σ Δs {:+.4} (Δh {:+.0} m) · margin Σ Δs {:+.4} (Δh {:+.0} m)", r.0, r.1, r.2, r.3, r.4);
    }
    eprintln!("      observed Δh_iso (final − t₀): interior {:+.0} m · margin {:+.0} m", mean_over(&obs, &interior), mean_over(&obs, &margin));
    let neg: f64 = rows.iter().map(|r| r.2.min(0.0)).sum();
    let surf_sink: f64 = rows.iter().filter(|r| r.0 == "Erosion" || r.0 == "EquilibriumHeight").map(|r| r.2.min(0.0)).sum();
    let q1 = surf_sink / neg.min(-1e-9);
    eprintln!("   → Q1: the surface term and the sinks make {:.0} % of the interior's negative Σ Δh ({:.0} of {:.0} m) → {}", 100.0 * q1, surf_sink, neg, if q1 >= 0.7 { "HOLDS" } else { "REFUTED → the « attribution » stop" });
    // A3
    eprintln!("\n   A3 — h_iso(t₀) on each seed's land:");
    let mut runs_c1: Vec<(u64, C1Sources, C1Sources, C1Sources)> = Vec::new();
    for seed in SEEDS {
        let (a, _, _) = c1_run(seed, &[]);
        let (s, _, _) = c1_run(seed, &[C1Term::Erosion]);
        let (s2, _, _) = c1_run(seed, &[C1Term::Erosion, C1Term::EquilibriumHeight]);
        let h0: Vec<f32> = a.hist.normalized(tlf)[0].data.iter().map(|v| (v - 0.5) * k).filter(|&v| v > 0.0).collect();
        let hs = sorted(h0);
        eprintln!(
            "   seed {seed}: land {} cells · mean {:.0} m · deciles {}",
            hs.len(),
            hs.iter().sum::<f32>() / hs.len().max(1) as f32,
            (1..10).map(|d| format!("{:.0}", pct(&hs, d as f64 / 10.0))).collect::<Vec<_>>().join("/")
        );
        runs_c1.push((seed, a, s, s2));
    }
    let snaps = |c: &C1Sources| -> Vec<GridF32> { c.hist.normalized(tlf).iter().map(|g| roll(g, off)).collect() };

    // ── T calibrated once on the témoin (sources) ──
    let s_tem = snaps(&runs_c1[0].2);
    let mut trials: Vec<(f64, f32)> = Vec::new();
    let clk = Clock::start();
    let mut peak_at = |t: f64| -> f32 {
        let r = hist_run("T trial", &s_tem, &base(PSEED), t);
        let p = peak_m(&r.zm);
        trials.push((t, p));
        eprintln!("      T trial {t:.3e} yr → peak {p:.0} m");
        p
    };
    eprintln!("\n   B3 — T, calibrated once on the témoin (sources):");
    let grid: Vec<f64> = [5.0f64, 6.0, 7.0, 8.0, 9.0].iter().map(|e| 10f64.powf(*e)).collect();
    let vals: Vec<f32> = grid.iter().map(|&t| peak_at(t)).collect();
    let ok = |v: f32| (2700.0..=3000.0).contains(&v);
    let mut t_cal = grid[vals.iter().enumerate().min_by(|a, b| (a.1 - 2850.0).abs().partial_cmp(&(b.1 - 2850.0).abs()).unwrap()).unwrap().0];
    if let Some(i) = (0..4).find(|&i| (vals[i] - 2850.0) * (vals[i + 1] - 2850.0) <= 0.0) {
        let (mut lo, mut hi, mut vlo) = (grid[i].log10(), grid[i + 1].log10(), vals[i]);
        if ok(vals[i]) {
            t_cal = grid[i];
        } else if ok(vals[i + 1]) {
            t_cal = grid[i + 1];
        } else {
            for _ in 0..7 {
                let mid = 0.5 * (lo + hi);
                let v = peak_at(10f64.powf(mid));
                t_cal = 10f64.powf(mid);
                if ok(v) {
                    break;
                }
                if (v - 2850.0) * (vlo - 2850.0) > 0.0 {
                    lo = mid;
                    vlo = v;
                } else {
                    hi = mid;
                }
            }
        }
    } else {
        eprintln!("      NOTE: no trial brackets 2 850 m; the nearest is kept");
    }
    let t_cost = clk.stop();
    eprintln!("   → T = {t_cal:.3e} yr ({} trials, {}) against C1's nominal 10–30 Ma (×{:.2}–×{:.2})", trials.len(), fmt_cost(t_cost), t_cal / 1e7, t_cal / 3e7);

    // ── the runs at 256² ──
    struct SeedRuns {
        seed: u64,
        r: Vec<Run>,
    }
    let mut all_runs: Vec<SeedRuns> = Vec::new();
    for (seed, a, s, s2) in &runs_c1 {
        let cfg = base(*seed);
        let clk = Clock::start();
        let (rec, cal) = physics_level(snaps(a).last().unwrap(), 256, &cfg, 2850.0, 300, &mut |_| {}, &|| false).unwrap();
        let f164 = run_of("F164 (steady)", to_m(&rec.z, k), cal.lifted_mask.clone(), clk.stop(), None);
        let f165 = hist_run("F165 (Δh_iso)", &snaps(a), &cfg, T_F165_YEARS);
        let src = hist_run("sources", &snaps(s), &cfg, t_cal);
        let src2 = hist_run("sources − EH", &snaps(s2), &cfg, t_cal);
        all_runs.push(SeedRuns { seed: *seed, r: vec![f164, f165, src, src2] });
    }
    eprintln!("\n   check: the témoin's F165 history recomputed: peak {:.0} m (F165: 2852), land {:.0} km² (F165: 20376)", peak_m(&all_runs[0].r[1].zm), all_runs[0].r[1].m.land_km2);
    eprintln!("\n   C — 256² (Europe: plain 51.6 · plateau 4.7 · hill 18.9 · mountain 24.8):");
    for sr in &all_runs {
        eprintln!("   ── seed {}:", sr.seed);
        let f164_land: Vec<bool> = ocean_mask(&sr.r[0].zm).iter().map(|o| !o).collect();
        let nf = f164_land.iter().filter(|&&b| b).count() as f32;
        for r in &sr.r {
            print_macro(&r.name, &r.m);
            let land: Vec<bool> = ocean_mask(&r.zm).iter().map(|o| !o).collect();
            let lost = (0..65536).filter(|&i| f164_land[i] && !land[i]).count() as f32 / nf;
            let gained = (0..65536).filter(|&i| land[i] && !f164_land[i]).count() as f32 / nf;
            let interior_le0 = (0..65536).filter(|&i| land[i] && r.zm.data[i] <= 0.0).count() as f32 / land.iter().filter(|&&b| b).count().max(1) as f32;
            let mtn: Vec<bool> = r.classes.iter().map(|&c| c == MOUNTAIN).collect();
            let (_, sizes) = components(&mtn, 256, 256);
            let big = sizes.iter().filter(|&&s| s as f32 * CELL * CELL > 100.0).count();
            let hy = hydrology(&r.zm, CELL);
            let lake_cells: Vec<bool> = hy.lake_id.iter().map(|&l| l != 0).collect();
            let la = lake_cells.iter().filter(|&&b| b).count();
            let on_plateau = (0..65536).filter(|&i| lake_cells[i] && r.classes[i] == PLATEAU).count();
            eprintln!(
                "      {:<22} peak {:.0} m · land lost {:.1} % · gained {:.1} % vs F164 · interior ≤ 0 m {:.2} % · mountain regions > 100 km² {} · lakes {} / {:.0} km² ({:.1} % on plateau cells) · {}",
                r.name,
                peak_m(&r.zm),
                100.0 * lost,
                100.0 * gained,
                100.0 * interior_le0,
                big,
                hy.lakes.len(),
                la as f32 * CELL * CELL,
                100.0 * on_plateau as f32 / la.max(1) as f32,
                fmt_cost(r.cost)
            );
            if let Some(u) = &r.sum_up {
                let pos: Vec<bool> = u.iter().map(|&v| v > 0.0).collect();
                let bidx = coast(&pos, 256, 256, 4);
                eprintln!("         the {{Σ U > 0}} boundary's block index (lattice 4 = the 64² blocks): {:.2}", bidx.block_index);
            }
        }
    }

    // ── the frozen chains ──
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    let ckm = DOMAIN_KM / 2048.0;
    let mut finals: Vec<(u64, Vec<ChainOut>)> = Vec::new();
    for sr in &all_runs {
        let cs: Vec<ChainOut> = sr.r.iter().map(|r| chain(&r.zm, &r.d5, frozen())).collect();
        finals.push((sr.seed, cs));
    }
    let mut successes = 0;
    let mut src_no_plateau = 0;
    for ((seed, cs), sr) in finals.iter().zip(&all_runs) {
        eprintln!("\n   ── seed {seed} at 2048² against Corsica 195 m:");
        let land_f164: Vec<bool> = ocean_mask(&cs[0].g).iter().map(|o| !o).collect();
        let k164 = coast(&land_f164, 2048, 2048, 32);
        for (ci, (c, r)) in cs.iter().zip(&sr.r).enumerate() {
            let m8 = up_u8(&r.classes, 256, 8);
            let excl: Vec<bool> = (0..2048 * 2048).map(|i| c.d5[i] || m8[i] != MOUNTAIN).collect();
            let whole = tex(&c.g, ckm, &c.d5);
            let mtn = tex(&c.g, ckm, &excl);
            for (tag, x, cc) in [("whole", &whole, &c_whole), ("mountain", &mtn, &c_mtn)] {
                let el = elements(x, cc, ckm, ctrl_ok);
                let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
                eprintln!("      {:<16} {tag:<9} {} of {} within ×1.5 · out: {}", r.name, el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") });
            }
            eprintln!("      {:<16} slope guard: peak {:.0} m · mountain-class slope p50/p90 {:.3}/{:.3} (Corsica {:.3}/{:.3})", r.name, peak_m(&c.g), mtn.r.slope_p50, mtn.r.slope_p90, c_mtn.r.slope_p50, c_mtn.r.slope_p90);
            let (rd, fl) = flat_map(&c.g, ckm, Some(&c.d5));
            let per: Vec<String> = (1..5u8)
                .map(|cl| {
                    let cells = (0..2048 * 2048).filter(|&i| rd[i] && m8[i] == cl).count();
                    let f = (0..2048 * 2048).filter(|&i| rd[i] && fl[i] && m8[i] == cl).count();
                    format!("{} {:.2} %", CLASS_NAMES[cl as usize], 100.0 * f as f32 / cells.max(1) as f32)
                })
                .collect();
            eprintln!("      {:<16} F by class: {}", r.name, per.join(" · "));
            let land: Vec<bool> = ocean_mask(&c.g).iter().map(|o| !o).collect();
            let kc = coast(&land, 2048, 2048, 32);
            eprintln!("      {:<16} coast {kc:?}", r.name);
            // the warp's depressions surviving to 2048²
            let hy = hydrology(&c.g, ckm);
            let surv: Vec<String> = c
                .fresh
                .iter()
                .map(|(nl, prints)| {
                    let f = 2048 / nl;
                    let alive = prints
                        .iter()
                        .filter(|p| p.iter().any(|&q| {
                            let (x, y) = (q % nl, q / nl);
                            (0..f).any(|dy| (0..f).any(|dx| hy.lake_id[(y * f + dy) * 2048 + x * f + dx] != 0))
                        }))
                        .count();
                    format!("{nl}²: {alive} of {} new still a lake", prints.len())
                })
                .collect();
            eprintln!("      {:<16} warp depressions: {} · {} · chain {}", r.name, c.notes.join(" | "), surv.join(" · "), fmt_cost(c.cost));
            if ci == 2 {
                // the success test (sources)
                let degraded = |a: f32, b: f32, cc: f32| (b - cc).abs() > (a - cc).abs() + 0.10 * a.abs();
                let deg: Vec<&str> = [
                    ("dimension", degraded(k164.dimension, kc.dimension, c_coast.dimension)),
                    ("length ratio", degraded(k164.length_ratio, kc.length_ratio, c_coast.length_ratio)),
                    ("aligned", degraded(k164.aligned, kc.aligned, c_coast.aligned)),
                    ("block index", degraded(k164.block_index, kc.block_index, c_coast.block_index)),
                ]
                .iter()
                .filter(|x| x.1)
                .map(|x| x.0)
                .collect();
                let f164_land: Vec<bool> = ocean_mask(&sr.r[0].zm).iter().map(|o| !o).collect();
                let l256: Vec<bool> = ocean_mask(&r.zm).iter().map(|o| !o).collect();
                let lost = (0..65536).filter(|&i| f164_land[i] && !l256[i]).count() as f32 / f164_land.iter().filter(|&&b| b).count() as f32;
                let m = &r.m;
                let mtn_ok = (0.165..=0.372).contains(&m.fractions[3]);
                let low = m.fractions[0] + m.fractions[1];
                let low_ok = (0.375..=0.845).contains(&low);
                let succ = mtn_ok && low_ok && lost < 0.15 && deg.is_empty();
                if succ {
                    successes += 1;
                    if m.fractions[1] < 0.01 {
                        src_no_plateau += 1;
                    }
                }
                eprintln!(
                    "      SUCCESS test (sources): mountain {:.1} % {} · plain + plateau {:.1} % {} · land lost {:.1} % {} · coast degraded {} → {}",
                    100.0 * m.fractions[3],
                    mtn_ok,
                    100.0 * low,
                    low_ok,
                    100.0 * lost,
                    lost < 0.15,
                    if deg.is_empty() { "none".into() } else { deg.join(", ") },
                    if succ { "SUCCESS" } else { "no" }
                );
            }
        }
    }
    eprintln!("\n   R: macro success on {successes} / 4 seeds ({src_no_plateau} of them with plateau < 1 %) · Q1 {}", if q1 >= 0.7 { "holds" } else { "refuted (the « attribution » stop)" });

    // ── seed 9's rim ──
    let s9 = &all_runs[3].r[2];
    let border: Vec<usize> = (0..65536).filter(|&i| { let (x, y) = (i % 256, i / 256); (x == 0 || y == 0 || x == 255 || y == 255) && s9.zm.data[i] > 0.0 }).collect();
    let alt_border = border.iter().map(|&i| s9.zm.data[i]).sum::<f32>() / border.len().max(1) as f32;
    let landc: Vec<usize> = (0..65536).filter(|&i| s9.zm.data[i] > 0.0).collect();
    let alt_land = landc.iter().map(|&i| s9.zm.data[i]).sum::<f32>() / landc.len().max(1) as f32;
    eprintln!(
        "\n   seed 9's rim (sources): {} land cells on the map border, mean {:.0} m against the land's {:.0} m. The cause: `erosion/stream_power.rs:651` sets a cell whose D8 receiver leaves the map as its own receiver (a fixed base node, never incised), and `terrain/flow.rs:256` routes D8 across the periodic wrap, so a border cell draining outward is such a node; under a positive U it only rises (the history does not hold the sea, so border land is not reset).",
        border.len(),
        alt_border,
        alt_land
    );

    // ── images ──
    let mut r256 = Vec::new();
    let mut r2048 = Vec::new();
    let mut rcls = Vec::new();
    let sh = |z: &GridF32, n: usize, px: usize| {
        let (_, hy) = read(z, None);
        whole(&rgba_of(z, Some(&hy), DOMAIN_KM / n as f32 * 1000.0, &ss), n, px)
    };
    for ((_, cs), sr) in finals.iter().zip(&all_runs) {
        r256.push((0..3).map(|i| sh(&sr.r[i].zm, 256, 320)).collect::<Vec<_>>());
        r2048.push((0..3).map(|i| sh(&cs[i].g, 2048, 448)).collect::<Vec<_>>());
        rcls.push((0..3).map(|i| class_png(&sr.r[i].classes, 256, 256, 1)).collect::<Vec<_>>());
    }
    compose(&r256).save(out().join("f166_seeds_256.png")).unwrap();
    compose(&r2048).save(out().join("f166_seeds_2048.png")).unwrap();
    compose(&rcls).save(out().join("f166_classes_seeds.png")).unwrap();
    // Σ uplift per term (the témoin, 64², Δh by the linear land ramp, rolled to the framing), ×6
    let mut tiles = Vec::new();
    for (term, ds) in &all.terms {
        let dh: Vec<f64> = (0..n * n).map(|i| land_altitude_m(all.s_final[i], all.craton[i], hs1, &iso) - land_altitude_m(all.s_final[i] - ds[i], all.craton[i], hs1, &iso)).collect();
        let d = roll_f(&dh, n, off);
        let mut a: Vec<f32> = d.iter().map(|v| v.abs()).collect();
        a.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let sat = a[((a.len() - 1) as f32 * 0.98) as usize].max(1.0);
        tiles.push(whole(&diff_rgba(&d, n, n, sat), n, 256));
        eprintln!("   term tile {term:?}: saturated at {sat:.0} m");
    }
    let dd: Vec<f64> = (0..n * n).map(datum).collect();
    let d = roll_f(&dd, n, off);
    let sat = d.iter().map(|v| v.abs()).fold(1.0f32, f32::max);
    tiles.push(whole(&diff_rgba(&d, n, n, sat), n, 256));
    let rows: Vec<Vec<image::RgbaImage>> = tiles.chunks(4).map(|c| c.to_vec()).collect();
    compose(&rows).save(out().join("f166_terms.png")).unwrap();
    // a plateau edge profile, if any plateau cell exists in the sources' runs
    let mut prof_done = false;
    for sr in &all_runs {
        let r = &sr.r[2];
        if let Some(i) = (0..65536).find(|&i| r.classes[i] == PLATEAU && i % 256 > 16 && i % 256 < 240) {
            let (x0, y0) = (i % 256 - 16, i / 256);
            let prof: Vec<f32> = (0..32).map(|dx| r.zm.data[y0 * 256 + x0 + dx]).collect();
            let (w, h) = (640u32, 300u32);
            let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
            let (lo, hi) = (prof.iter().copied().fold(f32::MAX, f32::min), prof.iter().copied().fold(f32::MIN, f32::max));
            for (j, v) in prof.iter().enumerate() {
                let x = 10 + j as u32 * 19;
                let y = h - 10 - ((v - lo) / (hi - lo).max(1.0) * (h - 20) as f32) as u32;
                let c = if r.classes[y0 * 256 + x0 + j] == PLATEAU { [220, 140, 30] } else { [60, 60, 60] };
                for dx in 0..14 {
                    for dy in 0..4 {
                        img.put_pixel((x + dx).min(w - 1), (y + dy).min(h - 1), image::Rgba([c[0], c[1], c[2], 255]));
                    }
                }
            }
            img.save(out().join("f166_plateau_profile.png")).unwrap();
            eprintln!("   image f166_plateau_profile.png: seed {}, row {y0}, x {x0}–{} (32 cells, 50 km), {:.0}–{:.0} m; plateau cells orange: {}", sr.seed, x0 + 31, lo, hi, prof.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join(" "));
            prof_done = true;
            break;
        }
    }
    if !prof_done {
        eprintln!("   no plateau cell in the sources' runs: no profile");
    }
    // the coast crops: the témoin's first 40–60 % land window, Corsica's
    let pick = |g: &GridF32, n: usize| -> (usize, usize) {
        let land: Vec<bool> = ocean_mask(g).iter().map(|o| !o).collect();
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
    let (f164g, srcg) = (&finals[0].1[0].g, &finals[0].1[2].g);
    let (ox, oy) = pick(srcg, 2048);
    let (cx, cy) = pick(&cg, 1024);
    let crop = |g: &GridF32, n: usize, x: usize, y: usize| {
        let (_, hy) = read(g, None);
        window(&rgba_of(g, Some(&hy), 195.3, &ss), n, (x as f32 / n as f32, y as f32 / n as f32), 128.0 / n as f32, 448)
    };
    compose(&[vec![crop(&cg, 1024, cx, cy), crop(f164g, 2048, ox, oy), crop(srcg, 2048, ox, oy)]]).save(out().join("f166_coast_25km.png")).unwrap();
    eprintln!(
        "   images: f166_seeds_256.png, f166_seeds_2048.png, f166_classes_seeds.png (rows = seeds {SEEDS:?}; columns F164 | F165 | sources), f166_terms.png (the témoin's Σ Δh per term: {}, then the datum drift; red +, blue −), f166_coast_25km.png (Corsica | F164 | sources; window ({ox}, {oy}) of 2048², Corsica ({cx}, {cy})) · north up · {NOTICE}",
        all.terms.iter().map(|t| format!("{:?}", t.0)).collect::<Vec<_>>().join(", ")
    );
    let _ = roll_mask;
    eprintln!("\n==========  end Finding 166 . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F166-C, the term image on the LAND only (the first run's tiles saturated on oceanic cells, where the linear land ramp
/// does not apply): the témoin's Σ Δh per C1 term over the final land, rolled to the framing, ×4, and their land means.
#[test]
#[ignore]
fn f166_terms_land() {
    let iso = IsostasyConfig::c1_default();
    let k = n2m();
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
    let tlf = bench_upscale_cfg(temoin).target_land_fraction;
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let (all, _, _) = c1_run(PSEED, &[]);
    let n = 64;
    let h_fin: Vec<f32> = all.hist.normalized(tlf).last().unwrap().data.iter().map(|v| (v - 0.5) * k).collect();
    let land: Vec<bool> = h_fin.iter().map(|&v| v > 0.0).collect();
    let (hs0, hs1) = (isostatic_datum(&all.s0, &all.craton, &iso), isostatic_datum(&all.s_final, &all.craton, &iso));
    let mut fields: Vec<(String, Vec<f64>)> = all
        .terms
        .iter()
        .map(|(term, ds)| {
            (format!("{term:?}"), (0..n * n).map(|i| if land[i] { land_altitude_m(all.s_final[i], all.craton[i], hs1, &iso) - land_altitude_m(all.s_final[i] - ds[i], all.craton[i], hs1, &iso) } else { 0.0 }).collect())
        })
        .collect();
    fields.push(("datum drift".into(), (0..n * n).map(|i| if land[i] { land_altitude_m(all.s_final[i], all.craton[i], hs1, &iso) - land_altitude_m(all.s_final[i], all.craton[i], hs0, &iso) } else { 0.0 }).collect()));
    let mut pool: Vec<f32> = fields.iter().flat_map(|f| (0..n * n).filter(|&i| land[i]).map(|i| f.1[i].abs() as f32).collect::<Vec<_>>()).collect();
    pool.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let sat = pool[((pool.len() - 1) as f32 * 0.98) as usize].max(1.0);
    let nl = land.iter().filter(|&&b| b).count() as f64;
    let mut tiles = Vec::new();
    for (name, f) in &fields {
        let d = roll_f(f, n, off);
        tiles.push(whole(&diff_rgba(&d, n, n, sat), n, 256));
        eprintln!("   {name:<18} land mean {:+.0} m · min {:+.0} · max {:+.0}", f.iter().sum::<f64>() / nl, f.iter().copied().fold(f64::MAX, f64::min), f.iter().copied().fold(f64::MIN, f64::max));
    }
    let rows: Vec<Vec<image::RgbaImage>> = tiles.chunks(4).map(|c| c.to_vec()).collect();
    compose(&rows).save(out().join("f166_terms.png")).unwrap();
    eprintln!("   image f166_terms.png: {} (land only, one shared saturation at the pooled p98, ±{sat:.0} m; red +, blue −), north up", fields.iter().map(|f| f.0.as_str()).collect::<Vec<_>>().join(", "));
}
