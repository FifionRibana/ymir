//! ADR Finding 169 -- a world that does not depend on its framing (the torus cascade, the offset at the end, the coast on
//! the torus), a continent that can be sailed around (the plates' selection), and the continent's age A.
//! Declared in `docs/reports/relief_method/f169_age/f169_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f169_cascade -- --ignored --exact f169_cascade --nocapture

#![allow(dead_code)]

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, PSEED, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::restrict;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, ocean_mask, ocean_mask_t};
use ymir_core::cascade::history::{
    Airy, C1Sources, C1Term, HistoryRun, Mapping, T_F166_YEARS, Timing, physics_history, run_c1_sources_mapped, step_durations,
};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{
    CLASS_NAMES, Circumnavigation, Coast, MOUNTAIN, MacroStats, PLAIN, PLATEAU, RELIEF_RADIUS_KM, circumnavigation, classify, coast,
    coast_torus, edt, hovius, local_relief, macro_stats,
};
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, flat_map, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::profile::{C1Profile, c1_run_config, select_continent};
use ymir_core::cascade::{CascadeConfig, Rebound, roll, upsample_to};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries, retarget_upper_plate_continental};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::land_topology::land_topology;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, run_with_closures_terms};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const SEEDS: [u64; 4] = [PSEED, 42, 1, 9];
const CELL: f32 = 1.5625;
const N_PHYS: usize = 300;
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];
const AGES: [(&str, Option<f64>); 4] = [("B0", None), ("young 20 Myr", Some(20.0e6)), ("intermediate 100 Myr", Some(100.0e6)), ("old 300 Myr", Some(300.0e6))];
const SHIFTS: [[usize; 2]; 3] = [[5, 0], [0, 23], [37, 50]];

fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f169_age/images")
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
fn up_mask(m: &[bool], n: usize, f: usize) -> Vec<bool> {
    let nn = n * f;
    (0..nn * nn).map(|k| m[(k / nn / f) * n + (k % nn) / f]).collect()
}
fn up_f(m: &[f32], n: usize, f: usize) -> Vec<f32> {
    let nn = n * f;
    (0..nn * nn).map(|k| m[(k / nn / f) * n + (k % nn) / f]).collect()
}
fn roll_mask(m: &[bool], n: usize, off: [usize; 2]) -> Vec<bool> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n]).collect()
}
fn roll_u8(m: &[u8], n: usize, off: [usize; 2]) -> Vec<u8> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n]).collect()
}
fn rolln(g: &GridF32, off64: [usize; 2]) -> GridF32 {
    let f = g.width / 64;
    roll(g, [off64[0] * f, off64[1] * f])
}
fn peak_m(z: &GridF32) -> f32 {
    z.data.iter().copied().fold(f32::MIN, f32::max)
}
/// The land on the torus (the ocean = the largest periodic sea component).
fn land_t(z: &GridF32) -> Vec<bool> {
    ocean_mask_t(z, Some([0, 0])).iter().map(|o| !o).collect()
}
fn land_w(z: &GridF32) -> Vec<bool> {
    ocean_mask(z).iter().map(|o| !o).collect()
}
fn cpu_secs() -> f64 {
    let pid = std::process::id();
    let o = std::process::Command::new("powershell").args(["-NoProfile", "-Command", &format!("(Get-Process -Id {pid}).TotalProcessorTime.TotalSeconds")]).output();
    o.ok().and_then(|o| String::from_utf8_lossy(&o.stdout).trim().replace(',', ".").parse().ok()).unwrap_or(f64::NAN)
}
struct Clock {
    wall: Instant,
    cpu: f64,
}
impl Clock {
    fn start() -> Self {
        Clock { wall: Instant::now(), cpu: cpu_secs() }
    }
    fn stop(&self) -> (f64, f64) {
        (self.wall.elapsed().as_secs_f64(), cpu_secs() - self.cpu)
    }
}
fn fmt_cost(c: (f64, f64)) -> String {
    format!("{:.1} s wall / {:.1} s CPU", c.0, c.1)
}
fn class_png(classes: &[u8], w: usize, h: usize, scale: u32) -> image::RgbaImage {
    let pal = [[40u8, 70, 120], [120, 170, 90], [200, 180, 110], [170, 120, 70], [235, 235, 235]];
    let mut img = image::RgbaImage::new(w as u32 * scale, h as u32 * scale);
    for y in 0..h as u32 * scale {
        for x in 0..w as u32 * scale {
            // north up: row 0 of the image is the grid's last row
            let gy = h - 1 - (y / scale) as usize;
            let c = pal[classes[gy * w + (x / scale) as usize] as usize];
            img.put_pixel(x, y, image::Rgba([c[0], c[1], c[2], 255]));
        }
    }
    img
}
fn print_macro(tag: &str, m: &MacroStats) {
    eprintln!(
        "   {tag:<26} plain {:5.1} % · plateau {:5.1} % · hill {:5.1} % · mountain {:5.1} % · land {:.0} km² · deciles {}",
        100.0 * m.fractions[0],
        100.0 * m.fractions[1],
        100.0 * m.fractions[2],
        100.0 * m.fractions[3],
        m.land_km2,
        m.deciles.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join("/")
    );
}
fn cheb_periodic(src: &[bool], n: usize) -> Vec<u32> {
    let mut d = vec![u32::MAX; n * n];
    let mut front: Vec<usize> = (0..n * n).filter(|&k| src[k]).collect();
    for &k in &front {
        d[k] = 0;
    }
    let mut lvl = 0;
    while !front.is_empty() {
        lvl += 1;
        let mut next = Vec::new();
        for &k in &front {
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
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
fn edt_periodic(land: &[bool], n: usize) -> Vec<f32> {
    let t = 3 * n;
    let tiled: Vec<bool> = (0..t * t).map(|k| land[((k / t) % n) * n + (k % t) % n]).collect();
    let d = edt(&tiled, t, t);
    (0..n * n).map(|k| d[(k / n + n) * t + k % n + n]).collect()
}
fn viz_offset(z64_norm: &GridF32) -> [usize; 2] {
    let (cx, cy) = land_topology(z64_norm, 0.5).center_cell;
    [(cx + 64 - 32) % 64, (cy + 64 - 32) % 64]
}
fn fmt_circ(c: &Circumnavigation, km: f32) -> String {
    match (c.wraps, c.sea_r) {
        (true, _) => format!("WRAPS (main {} cells)", c.main_cells),
        (false, Some(r)) => format!("circumnavigable (main {} cells) · r* {r} → sea loop ≥ {:.1} km", c.main_cells, (2 * r + 1) as f32 * km),
        (false, None) => format!("no land (main {} cells)", c.main_cells),
    }
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

/// The frozen chain (F164's setting) to 2 048².
fn frozen(pi_table: &[(usize, f32)], torus: Option<[usize; 2]>) -> AmpConfig {
    AmpConfig {
        talus_fine_scale: 0.25,
        n1_recal_depth: 2,
        n1_pi_m: pi_table.to_vec(),
        warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }),
        n1_k: vec![(512, 5.623e2), (1024, 1.0e3), (2048, 3.162e2)],
        torus,
        ..AmpConfig::declared(PSEED, Variant::N1, 1.0)
    }
}
fn chain(pm: &GridF32, d5: &[bool], cfg: AmpConfig) -> (GridF32, Vec<bool>, (f64, f64)) {
    let clk = Clock::start();
    let mut ch = Chain::new(pm.clone(), cfg).with_d5(d5.to_vec());
    for _ in 0..3 {
        ch.next_level(&|| false).unwrap();
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let last = ch.levels.last().unwrap();
    (last.result(), last.d5.clone(), clk.stop())
}
fn base_cfg(seed: u64, rough_a: f32, torus: Option<[usize; 2]>) -> CascadeConfig {
    CascadeConfig { seed, roughness_m: rough_a, mfd_exponent: Some(6.0), torus, ..CascadeConfig::declared(SteinSteinParams::default().depth_scale_m as f32) }
}

struct Run {
    zm: GridF32,
    d5: Vec<bool>,
    classes: Vec<u8>,
    m: MacroStats,
    cost: (f64, f64),
    rebound: Vec<f32>,
}
fn hist_run(snaps: &[GridF32], cfg: &CascadeConfig, age: Option<f64>) -> Run {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: None };
    let l = physics_history(snaps, 256, cfg, &run, &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, cfg.torus.map(|o| [o[0] * 4 % 256, o[1] * 4 % 256]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    Run { zm, d5, classes, m, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }
}

/// C1 under a profile, its Airy sources' record, in C1's own frame.
struct World {
    src: C1Sources,
    st: C1State,
    kin: PlateKinematics,
    snaps: Vec<GridF32>,
    cost: (f64, f64),
}
fn world(seed: u64, p: C1Profile) -> World {
    let cfg = c1_run_config();
    let mut st = p.init_state(64, seed);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let clk = Clock::start();
    let src = run_c1_sources_mapped(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &SteinSteinParams::default(), &[C1Term::Erosion], Mapping::Airy(Airy::declared()));
    let cost = clk.stop();
    let snaps = src.hist.normalized(None);
    World { src, st, kin, snaps, cost }
}

/// B1: the chronology of one C1 run -- per step and per term the Σ Δs per cell (f32), and the C-C collision cells at
/// each snapshot (every 10 steps).
struct Chrono {
    /// [term][step] -> per-cell Δs
    per_step: Vec<Vec<Vec<f32>>>,
    collisions: Vec<(usize, usize)>,
}
const CHRONO_TERMS: [C1Term; 4] = [C1Term::DavisSuppe, C1Term::Subduction, C1Term::RiftThinning, C1Term::Transport];
fn chrono(seed: u64, p: C1Profile) -> Chrono {
    let cfg = c1_run_config();
    let mut st = p.init_state(64, seed);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let cur = std::cell::RefCell::new(vec![vec![0f32; 4096]; 4]);
    let per_step = std::cell::RefCell::new(vec![Vec::new(); 4]);
    let mut collisions = Vec::new();
    let mut on_term = |t: C1Term, a: &[f64], b: &[f64]| {
        if let Some(i) = CHRONO_TERMS.iter().position(|x| *x == t) {
            let mut c = cur.borrow_mut();
            for k in 0..4096 {
                c[i][k] += (b[k] - a[k]) as f32;
            }
        }
    };
    run_with_closures_terms(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        |step, s, kn| {
            let mut c = cur.borrow_mut();
            let mut ps = per_step.borrow_mut();
            for i in 0..4 {
                ps[i].push(std::mem::replace(&mut c[i], vec![0f32; 4096]));
            }
            if (step + 1) % 10 == 0 {
                collisions.push((step + 1, collision_cells(s, kn).0));
            }
        },
        &mut on_term,
    );
    Chrono { per_step: per_step.into_inner(), collisions }
}
/// F165's rule: Convergent cells with a continental cell on both sides; and the convergent mask.
fn collision_cells(st: &C1State, kin: &PlateKinematics) -> (usize, Vec<bool>, usize) {
    let n = 64;
    let mut bi = classify_boundaries(&st.plate_id, kin);
    retarget_upper_plate_continental(&mut bi, &st.plate_id, &st.plate_type);
    let cont = |i: usize, j: usize| st.plate_type.get(i, j) == PlateType::Continental;
    let mut coll = 0;
    let mut conv = vec![false; n * n];
    for j in 0..n {
        for i in 0..n {
            if bi.boundary_type.get(i, j) == BoundaryType::Convergent {
                conv[j * n + i] = true;
                let me_ = st.plate_id.get(i, j);
                let nb = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
                if cont(i, j) && nb.iter().any(|&(x, y)| st.plate_id.get(x, y) != me_ && cont(x, y)) {
                    coll += 1;
                }
            }
        }
    }
    (coll, conv, bi.upper_plate_count())
}

/// A line plot of profiles (x km, z m).
fn profile_png(series: &[(Vec<f32>, [u8; 3])], km_per_px: f32, path: &std::path::Path) {
    let (w, h) = (900u32, 360u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let zmax = series.iter().flat_map(|s| s.0.iter()).copied().fold(500f32, f32::max);
    let zmin = series.iter().flat_map(|s| s.0.iter()).copied().fold(0f32, f32::min).max(-1000.0);
    let py = |z: f32| (h as f32 - 20.0 - (z - zmin) / (zmax - zmin) * (h as f32 - 40.0)) as i64;
    let mut put = |x: i64, y: i64, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as u32) < w && (y as u32) < h {
            img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    };
    for x in 20..w as i64 - 20 {
        put(x, py(0.0), [0, 0, 0]);
        put(x, py(1000.0), [210, 210, 210]);
        put(x, py(2000.0), [210, 210, 210]);
    }
    let _ = km_per_px;
    for (s, c) in series {
        for (i, &z) in s.iter().enumerate() {
            let x = 20 + (i as f32 / (s.len() - 1).max(1) as f32 * (w as f32 - 40.0)) as i64;
            for d in -1..=1 {
                put(x, py(z) + d, *c);
            }
        }
    }
    img.save(path).unwrap();
}

#[test]
#[ignore]
fn f169_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let k = n2m();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 169 . the torus cascade, the plates' selection, the continent's age  ==========");
    // references
    let e = load_grid(&root().join("data/europe/europe_1920.bin")).expect("data/europe");
    let me = macro_stats(&e.data, &classify(&e.data, e.width, e.height, CELL), e.width, e.height, CELL);
    print_macro("Europe", &me);
    let europe_alts = sorted(e.data.iter().copied().filter(|&v| v > 0.0).collect());
    let cg = load_grid(&root().join("data/corsica/corse_1024.bin")).unwrap();
    let c512 = load_grid(&root().join("data/corsica/corse_512.bin")).unwrap();
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let cclass = up_u8(&classify(&c128.data, 128, 128, CELL), 128, 8);
    let cmtn: Vec<bool> = cclass.iter().map(|&c| c == MOUNTAIN).collect();
    let ck = 200.0 / 1024.0;
    let c_whole = tex(&cg, ck, &vec![false; 1024 * 1024]);
    let c_mtn = tex(&cg, ck, &cmtn.iter().map(|m| !m).collect::<Vec<_>>());
    let c_coast = coast(&land_w(&cg), 1024, 1024, 32);
    let cp391 = predictability(&c512, 200.0 / 512.0, None);
    eprintln!("   Corse 195 m: mountain-class slope p50/p90 {:.3}/{:.3} · coast {:?}", c_mtn.r.slope_p50, c_mtn.r.slope_p90, c_coast);
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let ckm = DOMAIN_KM / 2048.0;
    let ctrl = facet_control(2048, ckm);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    // C2: the time step
    for (name, age) in AGES {
        if let Some(a) = age {
            let (d, r) = step_durations(T_F166_YEARS, a, N_PHYS);
            let cfg = base_cfg(0, rough_a, None);
            let cell_m = cfg.cell_km(256) * 1000.0;
            let w = cfg.diffusivity(256) * d[0] as f32 / (cell_m * cell_m);
            let subs = cfg.diffusion_substeps.max((w / 0.2).ceil() as usize);
            // the indicative « Courant » of the implicit incision at the oldest step, A = the whole land (upper bound)
            let f = cfg.k * d[0] as f32 * (100_000f32).sqrt() / cell_m;
            eprintln!("   C2 {name}: r {r:.6} · the oldest step d_1 = {:.3} Myr · the diffusion weight {w:.3} → {subs} sub-steps · the incision's K·A^m·dt/dx at A = 100 000 km² {f:.1} (implicit, indicative)", d[0] / 1e6);
        }
    }

    // ── A4: the selection ──
    eprintln!("\n   A4 — the continental plates' selection (land 50–62 %, circumnavigable, the largest r*):");
    let mut chosen: Vec<Option<C1Profile>> = Vec::new();
    for seed in SEEDS {
        let clk = Clock::start();
        let s = select_continent(seed);
        let c = clk.stop();
        for (i, cd) in s.candidates.iter().enumerate() {
            eprintln!(
                "      seed {seed} f {:.3} start {} (RNG rank {}): land {:.1} % · {} {}",
                cd.fraction,
                cd.start,
                cd.rng_rank,
                100.0 * cd.land,
                fmt_circ(&cd.circ, 6.25),
                if s.chosen == Some(i) { "← CHOSEN" } else if cd.qualifies() { "(qualifies)" } else { "" }
            );
        }
        eprintln!("      seed {seed}: {} · the selection's cost {}", if s.chosen.is_some() { "selected" } else { "NO CANDIDATE QUALIFIES" }, fmt_cost(c));
        chosen.push(s.chosen.map(|i| s.candidates[i].profile()));
    }
    let failed = chosen.iter().filter(|c| c.is_none()).count();
    if failed >= 2 {
        eprintln!("\n   R: the selection fails on {failed} / 4 seeds → STOP (next: more plates)");
        return;
    }

    // ── the worlds ──
    struct Done {
        name: &'static str,
        r: Run,
        g: GridF32,
        d5: Vec<bool>,
        cost: (f64, f64),
    }
    struct Seed {
        seed: u64,
        profile: C1Profile,
        w: World,
        ch: Chrono,
        off: [usize; 2],
        runs: Vec<Done>,
        f167_coast: Coast,
        inv: Vec<String>,
        inv_ok: bool,
    }
    let mut seeds: Vec<Seed> = Vec::new();
    for (si, seed) in SEEDS.iter().copied().enumerate() {
        let Some(profile) = chosen[si] else {
            eprintln!("   seed {seed}: no selected continent, skipped");
            continue;
        };
        let w = world(seed, profile);
        let ch = chrono(seed, profile);
        let off = viz_offset(w.snaps.last().unwrap());
        let mut runs = Vec::new();
        for (name, age) in AGES {
            let r = hist_run(&w.snaps, &base_cfg(seed, rough_a, Some([0, 0])), age);
            let (g, d5, cost) = chain(&r.zm, &r.d5, frozen(&pi_table, Some([0, 0])));
            eprintln!(
                "   seed {seed} {name}: C1 {} · physics {} · chain {} → {:.0} s CPU per world",
                fmt_cost(w.cost),
                fmt_cost(r.cost),
                fmt_cost(cost),
                w.cost.1 + r.cost.1 + cost.1
            );
            runs.push(Done { name, r, g, d5, cost });
        }
        // A1: the invariance at full scale on B0
        let mut inv = Vec::new();
        let mut inv_ok = true;
        for d in SHIFTS {
            let rolled: Vec<GridF32> = w.snaps.iter().map(|g| roll(g, d)).collect();
            let r = hist_run(&rolled, &base_cfg(seed, rough_a, Some(d)), None);
            let (g, _, _) = chain(&r.zm, &r.d5, frozen(&pi_table, Some(d)));
            let (a256, a2048) = (rolln(&runs[0].r.zm, d), rolln(&runs[0].g, d));
            let diff = |x: &GridF32, y: &GridF32| {
                let c = x.data.iter().zip(&y.data).filter(|(p, q)| p.to_bits() != q.to_bits()).count();
                let m = x.data.iter().zip(&y.data).map(|(p, q)| (p - q).abs()).fold(0f32, f32::max);
                (c, m)
            };
            let (c1, m1) = diff(&a256, &r.zm);
            let (c2, m2) = diff(&a2048, &g);
            inv_ok &= c1 == 0 && c2 == 0;
            inv.push(format!("shift {d:?}: 256² {c1} cells differ (max {m1} m) · 2048² {c2} cells differ (max {m2} m)"));
        }
        // the F167 reference coast, on the torus (production, the canonical roll, no torus, linear time)
        let coff = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
        let w0 = world(seed, C1Profile::Production);
        let sn0: Vec<GridF32> = w0.snaps.iter().map(|g| roll(g, coff)).collect();
        let r0 = hist_run(&sn0, &base_cfg(seed, rough_a, None), None);
        let (g0, _, _) = chain(&r0.zm, &r0.d5, frozen(&pi_table, None));
        let f167_coast = coast_torus(&land_w(&g0), 2048, 2048, 32);
        seeds.push(Seed { seed, profile, w, ch, off, runs, f167_coast, inv, inv_ok });
    }

    // ── A: infrastructure ──
    eprintln!("\n   A1 — the invariance test (bit for bit declared):");
    let mut inv_pass = 0;
    for s in &seeds {
        for l in &s.inv {
            eprintln!("      seed {} {l}", s.seed);
        }
        inv_pass += s.inv_ok as usize;
    }
    eprintln!("      → bit for bit on {inv_pass} / {} seeds", seeds.len());
    if inv_pass < seeds.len() {
        eprintln!("\n   R: the invariance test FAILS → STOP (no physical measure on a non-invariant cascade)");
        return;
    }

    // ── B1: the chronology ──
    eprintln!("\n   B1 — when C1 thickens (end land, 64²):");
    let mut thick_age: Vec<Vec<f32>> = Vec::new();
    for s in &seeds {
        let land64: Vec<bool> = s.w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
        let steps = s.ch.per_step[0].len();
        let tot = |ti: usize, a: usize, b: usize| -> f64 { (a..b).map(|st| (0..4096).filter(|&i| land64[i]).map(|i| s.ch.per_step[ti][st][i] as f64).sum::<f64>()).sum() };
        let sums: Vec<String> = CHRONO_TERMS
            .iter()
            .enumerate()
            .map(|(ti, t)| format!("{t:?} Σ {:.1} (thirds {:.1} / {:.1} / {:.1})", tot(ti, 0, steps), tot(ti, 0, 100), tot(ti, 100, 200), tot(ti, 200, steps)))
            .collect();
        let ds_all = tot(0, 0, steps).max(1e-12);
        let late = |f: f64| tot(0, steps - (steps as f64 * f).round() as usize, steps) / ds_all;
        eprintln!("      seed {}: {}", s.seed, sums.join(" · "));
        eprintln!(
            "      seed {}: the land's Σ Δs_DS received in the last 10 % {:.1} % · last 25 % {:.1} % · last 50 % {:.1} % (first half {:.1} %)",
            s.seed,
            100.0 * late(0.10),
            100.0 * late(0.25),
            100.0 * late(0.50),
            100.0 * (1.0 - late(0.50))
        );
        eprintln!("      seed {}: C-C collision cells at the snapshots (step, cells): {:?}", s.seed, s.ch.collisions);
        // the thickening age per cell (the step, 1-based, at which the cumulative Σ Δs_DS reaches half its total)
        let ages: Vec<f32> = (0..4096)
            .map(|i| {
                let total: f32 = (0..steps).map(|st| s.ch.per_step[0][st][i]).sum();
                if !land64[i] || total <= 0.05 {
                    return f32::NAN;
                }
                let mut c = 0f32;
                for st in 0..steps {
                    c += s.ch.per_step[0][st][i];
                    if c >= 0.5 * total {
                        return (st + 1) as f32;
                    }
                }
                steps as f32
            })
            .collect();
        let v = sorted(ages.iter().copied().filter(|a| a.is_finite()).collect());
        eprintln!(
            "      seed {}: the thickening age (step) over {} land cells with Σ Δs_DS > 0.05: p10 {:.0} · p50 {:.0} · p90 {:.0}",
            s.seed,
            v.len(),
            pct(&v, 0.1),
            pct(&v, 0.5),
            pct(&v, 0.9)
        );
        thick_age.push(ages);
    }

    // ── D: the measures ──
    let (mut successes, mut mtn_over, mut erased) = (0usize, 0usize, 0usize);
    let mut reg_rows = Vec::new();
    for (s, ages64) in seeds.iter().zip(&thick_age) {
        eprintln!("\n   ══ seed {} ({:?}), display offset {:?} ══", s.seed, s.profile, s.off);
        let st = &s.w.st;
        let crat = up_mask(st.cratonic_mask.data(), 64, 4);
        let ds = &s.w.src.terms.iter().find(|x| x.0 == C1Term::DavisSuppe).unwrap().1;
        let arc = &s.w.src.terms.iter().find(|x| x.0 == C1Term::Subduction).unwrap().1;
        let dsz = up_mask(&ds.iter().map(|&v| v > 0.05).collect::<Vec<_>>(), 64, 4);
        let arcm = up_mask(&arc.iter().map(|&v| v > 0.05).collect::<Vec<_>>(), 64, 4);
        let age256 = up_f(ages64, 64, 4);
        let (_, conv, upper) = collision_cells(st, &s.w.kin);
        let near: Vec<bool> = (0..4096).map(|i| conv[i] || s.w.src.hist.suture[i]).collect();
        let dnear = up_f(&cheb_periodic(&near, 64).iter().map(|&d| d as f32 * 6.25).collect::<Vec<_>>(), 64, 4);
        let land64: Vec<bool> = s.w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
        let c64 = circumnavigation(&land64, 64);
        eprintln!(
            "      C1: plates at the end {} · C-C collision cells {} · upper-plate cells {upper} · circumnavigation 64² {}",
            st.plate_id.data().iter().copied().collect::<std::collections::BTreeSet<u16>>().len(),
            collision_cells(st, &s.w.kin).0,
            fmt_circ(&c64, 6.25)
        );
        let t0m: Vec<f32> = upsample_to(&s.w.snaps[0], 256).data.iter().map(|v| (v - 0.5) * k).collect();
        let mut hyps = vec![("Europe".to_string(), europe_alts.clone(), [0u8, 0, 0])];
        for (ai, d) in s.runs.iter().enumerate() {
            let r = &d.r;
            eprintln!("   ── {}:", d.name);
            let land = land_t(&r.zm);
            let nl = land.iter().filter(|&&b| b).count().max(1);
            let land2048 = land_t(&d.g);
            let c2048 = circumnavigation(&land2048, 2048);
            // the rim: display window and C1's frame border
            let rim = |z: &GridF32, l: &[bool], n: usize, off: Option<[usize; 2]>| {
                let (zz, ll) = match off {
                    Some(o) => (rolln(z, o), roll_mask(l, n, [o[0] * n / 64, o[1] * n / 64])),
                    None => (z.clone(), l.to_vec()),
                };
                let b: Vec<usize> = (0..n * n).filter(|&k| ll[k] && (k % n == 0 || k % n == n - 1 || k / n == 0 || k / n == n - 1)).collect();
                let mb = b.iter().map(|&k| zz.data[k]).sum::<f32>() / b.len().max(1) as f32;
                let ml = (0..n * n).filter(|&k| ll[k]).map(|k| zz.data[k]).sum::<f32>() / ll.iter().filter(|&&x| x).count().max(1) as f32;
                (b.len(), mb / ml)
            };
            let (bd, rd) = rim(&d.g, &land2048, 2048, Some(s.off));
            let (bc, rc) = rim(&d.g, &land2048, 2048, None);
            eprintln!(
                "      land {:.1} % of the torus · circumnavigation 2048² {} · the rim 2048²: display window {bd} border land cells, ratio {rd:.2} · C1's frame border {bc} cells, ratio {rc:.2}",
                100.0 * nl as f32 / 65536.0,
                fmt_circ(&c2048, ckm)
            );
            print_macro(&format!("{} classes", d.name), &r.m);
            let peak = peak_m(&r.zm);
            let peak2048 = peak_m(&d.g);
            // the mountain class decomposed: cratons, Davis-Suppe zones by thickening third, arcs, the rest
            let mc: Vec<usize> = (0..65536).filter(|&i| r.classes[i] == MOUNTAIN).collect();
            let nm = mc.len().max(1) as f32;
            let third = |i: usize| -> Option<usize> { if age256[i].is_finite() { Some(((age256[i] - 1.0) / 100.0).floor().clamp(0.0, 2.0) as usize) } else { None } };
            let mut parts = [0f32; 6]; // cratons, DS 1st, 2nd, 3rd third, arcs, rest
            for &i in &mc {
                let p = if crat[i] {
                    0
                } else if dsz[i] {
                    1 + third(i).unwrap_or(2)
                } else if arcm[i] {
                    4
                } else {
                    5
                };
                parts[p] += 1.0;
            }
            eprintln!(
                "      peaks 256² {peak:.0} m · 2048² {peak2048:.0} m · the mountain class: cratons {:.1} % · DS thickened in the 1st third of the steps {:.1} % · 2nd {:.1} % · last {:.1} % · arcs {:.1} % · the rest {:.1} %",
                100.0 * parts[0] / nm,
                100.0 * parts[1] / nm,
                100.0 * parts[2] / nm,
                100.0 * parts[3] / nm,
                100.0 * parts[4] / nm,
                100.0 * parts[5] / nm
            );
            // plateaus and plains
            let dist = edt_periodic(&land, 256);
            let r7 = local_relief(&r.zm.data, 256, 256, CELL, RELIEF_RADIUS_KM);
            for (cl, nmc) in [(PLAIN, "plain"), (PLATEAU, "plateau")] {
                let cells: Vec<usize> = (0..65536).filter(|&i| r.classes[i] == cl).collect();
                if cells.is_empty() {
                    eprintln!("      {nmc}: none");
                    continue;
                }
                let mean = |f: &dyn Fn(usize) -> f32| cells.iter().map(|&i| f(i)).sum::<f32>() / cells.len() as f32;
                eprintln!(
                    "      {nmc}: {:.1} % of the land · distance to the coast {:.0} km · to a convergent cell / suture {:.0} km · altitude {:.0} m · R7 {:.0} m",
                    100.0 * cells.len() as f32 / nl as f32,
                    mean(&|i| dist[i] * CELL),
                    mean(&|i| dnear[i]),
                    mean(&|i| r.zm.data[i]),
                    mean(&|i| r7[i])
                );
            }
            // B2: « the rest » (B0)
            if ai == 0 {
                let rest: Vec<usize> = mc.iter().copied().filter(|&i| !crat[i] && !dsz[i] && !arcm[i]).collect();
                if !rest.is_empty() {
                    let mean = |f: &dyn Fn(usize) -> f32| rest.iter().map(|&i| f(i)).sum::<f32>() / rest.len() as f32;
                    let rb_land: f32 = (0..65536).filter(|&i| land[i]).map(|i| r.rebound[i]).sum();
                    let rb_rest: f32 = rest.iter().map(|&i| r.rebound[i]).sum();
                    eprintln!(
                        "      B2 « the rest » ({} cells = {:.1} % of the mountain class): altitude {:.0} m · R7 {:.0} m · distance to the coast {:.0} km · its share of Σ rebound {:.1} % (of {:.1} % of the land's cells) · mean rebound {:.0} m · its altitude at t₀ {:.0} m",
                        rest.len(),
                        100.0 * rest.len() as f32 / nm,
                        mean(&|i| r.zm.data[i]),
                        mean(&|i| r7[i]),
                        mean(&|i| dist[i] * CELL),
                        100.0 * rb_rest / rb_land.max(1e-6),
                        100.0 * rest.len() as f32 / nl as f32,
                        mean(&|i| r.rebound[i]),
                        mean(&|i| t0m[i])
                    );
                }
            }
            let hy = hydrology(&r.zm, CELL);
            eprintln!("      lakes (256²): {} / {:.0} km²", hy.lakes.len(), hy.lake_id.iter().filter(|&&l| l != 0).count() as f32 * CELL * CELL);
            // texture at 2048²
            let m8 = up_u8(&r.classes, 256, 8);
            let mtn8: Vec<bool> = m8.iter().map(|&x| x == MOUNTAIN).collect();
            let excl: Vec<bool> = (0..2048 * 2048).map(|i| d.d5[i] || !mtn8[i]).collect();
            let whole_t = tex(&d.g, ckm, &d.d5);
            let mtn_t = tex(&d.g, ckm, &excl);
            for (t, x, cc) in [("whole", &whole_t, &c_whole), ("mountain", &mtn_t, &c_mtn)] {
                let el = elements(x, cc, ckm, ctrl_ok);
                let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
                eprintln!("      2048² {t:<9} {} of {} within ×1.5 · out: {}", el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") });
            }
            eprintln!("      slope guard: mountain-class slope p50/p90 {:.3}/{:.3} (Corsica {:.3}/{:.3})", mtn_t.r.slope_p50, mtn_t.r.slope_p90, c_mtn.r.slope_p50, c_mtn.r.slope_p90);
            let (rd_, fl) = flat_map(&d.g, ckm, Some(&d.d5));
            let per: Vec<String> = (1..5u8)
                .map(|cl| {
                    let cells = (0..2048 * 2048).filter(|&i| rd_[i] && m8[i] == cl).count();
                    let f = (0..2048 * 2048).filter(|&i| rd_[i] && fl[i] && m8[i] == cl).count();
                    format!("{} {:.2} %", CLASS_NAMES[cl as usize], 100.0 * f as f32 / cells.max(1) as f32)
                })
                .collect();
            eprintln!("      F by class: {}", per.join(" · "));
            let kc = coast_torus(&land2048, 2048, 2048, 32);
            eprintln!("      coast (torus) {kc:?} · F167 on the torus {:?}", s.f167_coast);
            let mut alts = sorted((0..65536).filter(|&i| land[i]).map(|i| r.zm.data[i]).collect());
            alts.retain(|v| v.is_finite());
            hyps.push((d.name.to_string(), alts, [[120u8, 120, 120], [40, 150, 60], [40, 90, 200], [210, 50, 40]][ai]));
            if ai == 3 {
                // the regularity baseline (old)
                let p2048 = &whole_t.p;
                let p1024 = {
                    let g = restrict(&d.g);
                    let dd = (0..1024 * 1024).map(|i| d.d5[(i / 1024 * 2) * 2048 + (i % 1024) * 2]).collect::<Vec<_>>();
                    predictability(&g, 2.0 * ckm, Some(&dd))
                };
                let hv: Vec<String> = [10.0f32, 25.0, 100.0].iter().map(|&a| { let h = hovius(&d.g, &mtn8, ckm, a); format!("A ≥ {a}: R {:.2} ({} basins)", h.ratio, h.outlets) }).collect();
                reg_rows.push(format!(
                    "   seed {} old: 391 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | 195 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | Hovius {}",
                    s.seed, p1024.cv_lambda, p1024.cv_len, p1024.theta_mean, p1024.theta_std, p1024.r2_gully, p2048.cv_lambda, p2048.cv_len, p2048.theta_mean, p2048.theta_std, p2048.r2_gully, hv.join(" · ")
                ));
                let degraded = |a: f32, b: f32, cc: f32| (b - cc).abs() > (a - cc).abs() + 0.10 * a.abs();
                let fc = &s.f167_coast;
                let deg: Vec<&str> = [
                    ("dimension", degraded(fc.dimension, kc.dimension, c_coast.dimension)),
                    ("length ratio", degraded(fc.length_ratio, kc.length_ratio, c_coast.length_ratio)),
                    ("aligned", degraded(fc.aligned, kc.aligned, c_coast.aligned)),
                    ("block index", degraded(fc.block_index, kc.block_index, c_coast.block_index)),
                ]
                .iter()
                .filter(|x| x.1)
                .map(|x| x.0)
                .collect();
                let m = &r.m;
                let mtn_ok = (0.165..=0.372).contains(&m.fractions[3]);
                let low = m.fractions[0] + m.fractions[1];
                let low_ok = (0.375..=0.845).contains(&low);
                let circ = !c64.wraps && !c2048.wraps;
                let succ = mtn_ok && low_ok && circ && deg.is_empty();
                successes += succ as usize;
                mtn_over += (m.fractions[3] > 0.372) as usize;
                erased += (peak < 1500.0) as usize;
                eprintln!(
                    "      SUCCESS test (old): mountain {:.1} % {mtn_ok} · plain + plateau {:.1} % {low_ok} · circumnavigable {circ} · coast degraded vs F167 (torus) {} · median land {:.0} m (Europe 195) → {}",
                    100.0 * m.fractions[3],
                    100.0 * low,
                    if deg.is_empty() { "none".into() } else { deg.join(", ") },
                    m.deciles[4],
                    if succ { "SUCCESS" } else { "no" }
                );
            }
        }
        hyps_png(&hyps, &out().join(format!("f169_hyps_{}.png", s.seed)));
    }
    eprintln!(
        "\n   the regularity baseline (Corsica 391 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | 195 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3}):",
        cp391.cv_lambda, cp391.cv_len, cp391.theta_mean, cp391.theta_std, cp391.r2_gully, c_whole.p.cv_lambda, c_whole.p.cv_len, c_whole.p.theta_mean, c_whole.p.theta_std, c_whole.p.r2_gully
    );
    for r in &reg_rows {
        eprintln!("{r}");
    }
    let n = seeds.len();
    eprintln!(
        "\n   R: invariance bit for bit on {inv_pass} / {n} · selection on {} / 4 · macro success (old) on {successes} / {n} · mountains > 37.2 % (old) on {mtn_over} / {n} · peak < 1.5 km (old) on {erased} / {n} → {}",
        4 - failed,
        if successes >= 3 {
            "SUCCESS: F170 (the sea-level rise) follows"
        } else if erased >= 3 {
            "the old profile erases the high mountains: the chronology (B1) reported; next, the collisions' timing in C1"
        } else if mtn_over >= 1 && mtn_over >= n.saturating_sub(successes) {
            "the mountains still exceed 37.2 % on the old profile: decomposed by thickening age, stop"
        } else {
            "no success, no stop condition met as written"
        }
    );

    // ── images (north up, the offset at the end) ──
    let shade = |z: &GridF32, nn: usize, px: usize| {
        let (_, hy) = read(z, None);
        whole(&rgba_of(z, Some(&hy), DOMAIN_KM / nn as f32 * 1000.0, &ss), nn, px)
    };
    let (mut r256, mut r2048, mut rcls, mut rage) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (s, ages64) in seeds.iter().zip(&thick_age) {
        r256.push(s.runs.iter().map(|d| shade(&rolln(&d.r.zm, s.off), 256, 300)).collect::<Vec<_>>());
        r2048.push(s.runs.iter().map(|d| shade(&rolln(&d.g, s.off), 2048, 420)).collect::<Vec<_>>());
        rcls.push(s.runs.iter().map(|d| class_png(&roll_u8(&d.r.classes, 256, [s.off[0] * 4, s.off[1] * 4]), 256, 256, 1)).collect::<Vec<_>>());
        // the thickening-age map: blue = early, red = late, grey = land without Davis-Suppe, dark = sea
        let land64: Vec<bool> = s.w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
        let mut img = image::RgbaImage::new(256, 256);
        for y in 0..256u32 {
            for x in 0..256u32 {
                let gy = 63 - (y / 4) as usize;
                let (cx, cy) = (((x / 4) as usize + s.off[0]) % 64, (gy + s.off[1]) % 64);
                let i = cy * 64 + cx;
                let c = if !land64[i] {
                    [25, 35, 60]
                } else if ages64[i].is_finite() {
                    let t = (ages64[i] - 1.0) / 299.0;
                    [(40.0 + 200.0 * t) as u8, 60, (230.0 - 200.0 * t) as u8]
                } else {
                    [160, 160, 160]
                };
                img.put_pixel(x, y, image::Rgba([c[0], c[1], c[2], 255]));
            }
        }
        rage.push(vec![img]);
    }
    compose(&r256).save(out().join("f169_seeds_256.png")).unwrap();
    compose(&r2048).save(out().join("f169_seeds_2048.png")).unwrap();
    compose(&rcls).save(out().join("f169_classes_seeds.png")).unwrap();
    compose(&[rage.into_iter().flatten().collect()]).save(out().join("f169_thickening_age.png")).unwrap();
    // the cross-sections (old profile, the témoin or the first seed): the oldest and the most recent thickened belt
    if let (Some(s), Some(ages64)) = (seeds.first(), thick_age.first()) {
        let land64: Vec<bool> = s.w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
        let ds = &s.w.src.terms.iter().find(|x| x.0 == C1Term::DavisSuppe).unwrap().1;
        let cand: Vec<usize> = (0..4096).filter(|&i| land64[i] && ages64[i].is_finite() && ds[i] > 0.3).collect();
        let oldest = cand.iter().copied().min_by(|a, b| ages64[*a].partial_cmp(&ages64[*b]).unwrap());
        let newest = cand.iter().copied().max_by(|a, b| ages64[*a].partial_cmp(&ages64[*b]).unwrap());
        for (tag, cell) in [("old_belt", oldest), ("recent_belt", newest)] {
            let Some(c) = cell else { continue };
            let (cx, cy) = ((c % 64) * 32 + 16, (c / 64) * 32 + 16);
            let line = |g: &GridF32| -> Vec<f32> { (0..1024).map(|i| g.data[cy * 2048 + (cx + 2048 - 512 + i) % 2048]).collect() };
            profile_png(&[(line(&s.runs[0].g), [120, 120, 120]), (line(&s.runs[3].g), [210, 50, 40])], ckm, &out().join(format!("f169_section_{tag}.png")));
            eprintln!("   section {tag}: 64² cell ({}, {}) thickened at step {:.0} · W–E line of 200 km through it · B0 grey, old red", c % 64, c / 64, ages64[c]);
        }
    }
    eprintln!(
        "   images: f169_seeds_256.png, f169_seeds_2048.png, f169_classes_seeds.png (rows = the seeds; columns B0 | young | intermediate | old), f169_thickening_age.png (the seeds left to right; blue early → red late, grey = land without Davis-Suppe), f169_hyps_<seed>.png (Europe black, B0 grey, young green, intermediate blue, old red), f169_section_*.png · north up, the display offset applied at the end"
    );
    eprintln!("\n==========  end Finding 169 . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F169-A1–A3 at full scale, infrastructure only (no physical measure; run after the selection's stop): the invariance
/// test, the rim and the coast on the torus, on F168's L2 worlds (their land touches the window), torus on.
#[test]
#[ignore]
fn f169_infra() {
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    eprintln!("\n==========  Finding 169 . A1–A3 at full scale (F168 L2 worlds, torus on)  ==========");
    let mut pass = 0;
    for seed in SEEDS {
        let w = world(seed, C1Profile::L2);
        let off = viz_offset(w.snaps.last().unwrap());
        let r = hist_run(&w.snaps, &base_cfg(seed, rough_a, Some([0, 0])), None);
        let (g, _, cost) = chain(&r.zm, &r.d5, frozen(&pi_table, Some([0, 0])));
        let mut ok = true;
        for d in SHIFTS {
            let rolled: Vec<GridF32> = w.snaps.iter().map(|x| roll(x, d)).collect();
            let rs = hist_run(&rolled, &base_cfg(seed, rough_a, Some(d)), None);
            let (gs, _, _) = chain(&rs.zm, &rs.d5, frozen(&pi_table, Some(d)));
            let diff = |x: &GridF32, y: &GridF32| {
                let c = x.data.iter().zip(&y.data).filter(|(p, q)| p.to_bits() != q.to_bits()).count();
                let m = x.data.iter().zip(&y.data).map(|(p, q)| (p - q).abs()).fold(0f32, f32::max);
                (c, m)
            };
            let (c1, m1) = diff(&rolln(&r.zm, d), &rs.zm);
            let (c2, m2) = diff(&rolln(&g, d), &gs);
            ok &= c1 == 0 && c2 == 0;
            eprintln!("   seed {seed} shift {d:?}: 256² {c1} cells differ (max {m1} m) · 2048² {c2} cells differ (max {m2} m)");
        }
        pass += ok as usize;
        let land = land_t(&g);
        let rim = |o: Option<[usize; 2]>| {
            let (zz, ll) = match o {
                Some(o) => (rolln(&g, o), roll_mask(&land, 2048, [o[0] * 32, o[1] * 32])),
                None => (g.clone(), land.clone()),
            };
            let n = 2048;
            let b: Vec<usize> = (0..n * n).filter(|&k| ll[k] && (k % n == 0 || k % n == n - 1 || k / n == 0 || k / n == n - 1)).collect();
            let mb = b.iter().map(|&k| zz.data[k]).sum::<f32>() / b.len().max(1) as f32;
            let ml = (0..n * n).filter(|&k| ll[k]).map(|k| zz.data[k]).sum::<f32>() / ll.iter().filter(|&&x| x).count().max(1) as f32;
            (b.len(), mb, ml)
        };
        let (bd, mbd, ml) = rim(Some(off));
        let (bc, mbc, _) = rim(None);
        eprintln!(
            "   seed {seed}: the rim 2048² (F168 L2 at 2048²: 1.29–2.01×): display window {bd} border land cells, mean {mbd:.0} m vs land {ml:.0} m → ratio {:.2} · C1's frame border {bc} cells, mean {mbc:.0} m → ratio {:.2} · land {:.1} %",
            mbd / ml,
            mbc / ml,
            100.0 * land.iter().filter(|&&b| b).count() as f32 / (2048.0 * 2048.0)
        );
        eprintln!("   seed {seed}: coast on the torus {:?} · coast in the window (F168's measure) {:?} · peak {:.0} m · chain {}", coast_torus(&land, 2048, 2048, 32), coast(&land_w(&rolln(&g, off)), 2048, 2048, 32), peak_m(&g), fmt_cost(cost));
    }
    eprintln!("   → the invariance bit for bit on {pass} / 4 seeds");
    eprintln!("==========  end F169 infra  ==========\n");
}
