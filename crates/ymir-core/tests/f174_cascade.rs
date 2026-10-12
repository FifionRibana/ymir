//! ADR Finding 174 -- the coastline inherited from C1 (A1), real reference coasts by type (A2), the physics level's time step
//! (A3), the continent's perimeter at each C1 step (A4); the panel 8, 18, 19, 23, 29, K2, no sea. F173's helpers are kept.
//! Declared in `docs/reports/relief_method/f174_margin/f174_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f174_cascade -- --ignored --exact f174_measure --nocapture
#![allow(dead_code)]

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, PSEED, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::restrict;
use ymir_core::cascade::amplify::{AmpConfig, Chain, KHetero, Variant, Warp, WarpPlace, ocean_mask, ocean_mask_t};
use ymir_core::cascade::history::{
    Airy, C1Sources, C1Term, HistoryRun, Mapping, SeaCurve, SeaLow, T_F166_YEARS, Timing, physics_history, physics_history_rock,
    physics_history_sea, run_c1_sources_mapped,
    step_durations,
};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{
    CLASS_NAMES, Circumnavigation, Coast, MOUNTAIN, MacroStats, PLAIN, PLATEAU, RELIEF_RADIUS_KM, circumnavigation, classify, coast,
    coast_torus, edt, hovius, local_relief, macro_stats, periodic_components, slope_deg,
};
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, flat_map, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::profile::{C1Profile, Selection, c1_run_config, select_continent_default, select_continent_plates};
use ymir_core::cascade::{CascadeConfig, Rebound, Rock, roll, upsample_to};
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
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
/// F170-A1: the plate counts and their k brackets (⌈0.35 N⌉ … ⌊0.75 N⌋).
const PLATES: [(usize, std::ops::RangeInclusive<usize>); 2] = [(12, 5..=9), (16, 6..=12)];

/// F170-A4: the convergent + suture cells of a C1 end state (64²), × 6.25 km, per 10 000 km² of end land.
fn suture_density(w: &World) -> (usize, f32, f32) {
    let (_, conv, _) = collision_cells(&w.st, &w.kin);
    let near = (0..4096).filter(|&i| conv[i] || w.src.hist.suture[i]).count();
    let land = w.snaps.last().unwrap().data.iter().filter(|&&v| v > 0.5).count() as f32 * 6.25 * 6.25;
    (near, land, near as f32 * 6.25 / (land / 1.0e4))
}

fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f173_sea_curve/images")
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


// ════════════════════════════════════════ F171 ════════════════════════════════════════

/// F171-B2: the declared rock values.
const RHO_K: f32 = 10.0;
const COVER_M: f32 = 1500.0;

/// The cover map at 256² (m): d_c, or 0 on C1's cratonic shields (the init state's mask, nearest cell).
fn cover_map(seed: u64, p: C1Profile, d_c: f32) -> Vec<f32> {
    let st0 = p.init_state(64, seed);
    up_mask(st0.cratonic_mask.data(), 64, 4).iter().map(|&s| if s { 0.0 } else { d_c }).collect()
}
/// The physics level under a rock law.
fn hist_rock(snaps: &[GridF32], cfg: &CascadeConfig, age: Option<f64>, rock: Option<&Rock>) -> (Run, Vec<f32>) {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: None };
    let l = physics_history_rock(snaps, 256, cfg, &run, rock, &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, Some([0, 0]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    (Run { zm, d5, classes, m, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }, l.rec.e_cum_m)
}
fn chain_k(pm: &GridF32, d5: &[bool], cfg: AmpConfig, ln_k: Option<Vec<f32>>) -> (GridF32, Vec<bool>, (f64, f64)) {
    let clk = Clock::start();
    let mut ch = Chain::new(pm.clone(), cfg).with_d5(d5.to_vec());
    if let Some(l) = ln_k {
        ch = ch.with_ln_k(l);
    }
    for _ in 0..3 {
        ch.next_level(&|| false).unwrap();
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let last = ch.levels.last().unwrap();
    (last.result(), last.d5.clone(), clk.stop())
}

/// 0 = sea, 1 = cover, 2 = basement, 3 = shield.
fn rock_at_surface(land: &[bool], cover: &[f32], e_cum: &[f32]) -> Vec<u8> {
    (0..land.len())
        .map(|i| {
            if !land[i] {
                0
            } else if cover[i] <= 0.0 {
                3
            } else if !e_cum.is_empty() && e_cum[i] >= cover[i] {
                2
            } else {
                1
            }
        })
        .collect()
}

/// F171-D (declared): the lithological knickpoints on the 256² Strahler ≥ 3 network at the cover / basement contacts.
/// Returns (count, median height m, the knickpoints (cell, height)).
fn knickpoints(z: &GridF32, rock: &[u8]) -> (usize, f32, Vec<(usize, f32)>) {
    let n = 256;
    let hy = hydrology(z, CELL);
    let recv = |k: usize| -> Option<usize> {
        let d = hy.dir[k];
        if d == DIR_NONE {
            return None;
        }
        let (x, y) = ((k % n) as i32 + D8_DX[d as usize], (k / n) as i32 + D8_DY[d as usize]);
        Some(y.rem_euclid(n as i32) as usize * n + x.rem_euclid(n as i32) as usize)
    };
    // the main donor (the largest area) of each cell
    let mut main_donor = vec![usize::MAX; n * n];
    for k in 0..n * n {
        if let Some(r) = recv(k) {
            if main_donor[r] == usize::MAX || hy.acc_km2[k] > hy.acc_km2[main_donor[r]] {
                main_donor[r] = k;
            }
        }
    }
    let hard = |r: u8| r >= 2;
    let dist = |a: usize, b: usize| {
        let (dx, dy) = ((a % n) as i32 - (b % n) as i32, (a / n) as i32 - (b / n) as i32);
        let (dx, dy) = (dx.rem_euclid(n as i32).min((-dx).rem_euclid(n as i32)), dy.rem_euclid(n as i32).min((-dy).rem_euclid(n as i32)));
        ((dx * dx + dy * dy) as f32).sqrt() * CELL * 1000.0
    };
    // the mean slope and the drop over 3 cells from `start`, following `next`
    let reach = |start: usize, next: &dyn Fn(usize) -> Option<usize>| -> Option<(f32, f32)> {
        let (mut c, mut len, mut drop) = (start, 0f32, 0f32);
        for _ in 0..3 {
            let nx = next(c)?;
            len += dist(c, nx);
            drop += (z.data[c] - z.data[nx]).abs();
            c = nx;
        }
        Some((drop / len.max(1.0), drop))
    };
    let mut found = Vec::new();
    for k in 0..n * n {
        if hy.strahler[k] < 3 || rock[k] == 0 {
            continue;
        }
        let Some(r) = recv(k) else { continue };
        if rock[r] == 0 || hard(rock[k]) == hard(rock[r]) {
            continue;
        }
        // upstream of k (main donors) and downstream of r (receivers)
        let up = reach(k, &|c| { let d = main_donor[c]; (d != usize::MAX).then_some(d) });
        let dn = reach(r, &|c| recv(c));
        let (Some((su, hu)), Some((sd, hd))) = (up, dn) else { continue };
        let (s_hard, s_soft, h) = if hard(rock[r]) { (sd, su, hd) } else { (su, sd, hu) };
        if s_hard >= 2.0 * s_soft.max(1e-6) {
            found.push((k, h));
        }
    }
    let hs = sorted(found.iter().map(|f| f.1).collect());
    (found.len(), pct(&hs, 0.5), found)
}


// ════════════════════════════════════════ F172 ════════════════════════════════════════

const DELTA_M: f32 = 134.0;

/// The physics level under the rock law and (optionally) the lowstand.
fn hist_k2(snaps: &[GridF32], cfg: &CascadeConfig, age: Option<f64>, rock: &Rock, sea: Option<SeaLow>) -> (Run, Vec<f32>) {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: sea };
    let l = physics_history_rock(snaps, 256, cfg, &run, Some(rock), &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, Some([0, 0]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    (Run { zm, d5, classes, m, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }, l.rec.e_cum_m)
}
fn shifted(g: &GridF32, d: f32) -> GridF32 {
    GridF32 { width: g.width, height: g.height, data: g.data.iter().map(|v| v + d).collect() }
}
/// The frozen chain run on the field raised by Δ (its sea at the real −Δ), lowered by Δ at the end.
fn chain_low(zm: &GridF32, amp: AmpConfig, ln_k: Option<Vec<f32>>, delta: f32) -> (GridF32, Vec<bool>, (f64, f64)) {
    let zr = shifted(zm, delta);
    let oc = ocean_mask_t(&zr, Some([0, 0]));
    let d5: Vec<bool> = (0..zr.data.len()).map(|i| !oc[i] && zr.data[i] <= 0.0).collect();
    let (g, d5o, c) = chain_k(&zr, &d5, amp, ln_k);
    (shifted(&g, -delta), d5o, c)
}

/// Spearman's rank correlation (ties averaged).
fn spearman(a: &[f32], b: &[f32]) -> f32 {
    let rank = |v: &[f32]| -> Vec<f32> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&i, &j| v[i].partial_cmp(&v[j]).unwrap());
        let mut r = vec![0f32; v.len()];
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            for k in i..=j {
                r[idx[k]] = (i + j) as f32 / 2.0;
            }
            i = j + 1;
        }
        r
    };
    let (ra, rb) = (rank(a), rank(b));
    let n = a.len() as f32;
    let (ma, mb) = (ra.iter().sum::<f32>() / n, rb.iter().sum::<f32>() / n);
    let cov: f32 = ra.iter().zip(&rb).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let va: f32 = ra.iter().map(|x| (x - ma).powi(2)).sum();
    let vb: f32 = rb.iter().map(|y| (y - mb).powi(2)).sum();
    cov / (va * vb).sqrt().max(1e-12)
}

/// A1: the perimeter's boundary shares (convergent, divergent, transform, none) at 64², and Σ Δs_DS by the distance to the
/// coast (cells 0–7, ≥ 8).
fn perimeter(w: &World) -> ([f32; 4], Vec<(f32, usize)>) {
    let n = 64;
    let land: Vec<bool> = w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
    let mut bi = classify_boundaries(&w.st.plate_id, &w.kin);
    retarget_upper_plate_continental(&mut bi, &w.st.plate_id, &w.st.plate_type);
    let mut counts = [0usize; 4];
    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
            if !land[k] {
                continue;
            }
            let nb4 = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
            if !nb4.iter().any(|&(x, y)| !land[y * n + x]) {
                continue;
            }
            let mut best = 3usize;
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let b = bi.boundary_type.get((i + dx) % n, (j + dy) % n);
                    let c = match b {
                        BoundaryType::Convergent => 0,
                        BoundaryType::Divergent => 1,
                        BoundaryType::Internal => 3,
                        _ => 2,
                    };
                    best = best.min(c);
                }
            }
            counts[best] += 1;
        }
    }
    let tot = counts.iter().sum::<usize>().max(1) as f32;
    let ds = &w.src.terms.iter().find(|x| x.0 == C1Term::DavisSuppe).unwrap().1;
    let dist = edt_periodic(&land, 64);
    let mut bins = vec![(0f32, 0usize); 9];
    for k in 0..4096 {
        if land[k] {
            let b = ((dist[k] - 0.5).max(0.0).floor() as usize).min(8);
            bins[b].0 += ds[k] as f32;
            bins[b].1 += 1;
        }
    }
    ([counts[0] as f32 / tot, counts[1] as f32 / tot, counts[2] as f32 / tot, counts[3] as f32 / tot], bins.iter().map(|&(s, c)| (if c > 0 { s / c as f32 } else { f32::NAN }, c)).collect())
}

/// A2: (f10, f20, f_in, R, rampart) on a 256² field.
fn rampart(zm: &GridF32, classes: &[u8]) -> (f32, f32, f32, f32, bool) {
    let land = land_t(zm);
    let d = edt_periodic(&land, 256);
    let share = |f: &dyn Fn(usize) -> bool| -> f32 {
        let cells: Vec<usize> = (0..65536).filter(|&i| land[i] && f(i)).collect();
        cells.iter().filter(|&&i| classes[i] == MOUNTAIN).count() as f32 / cells.len().max(1) as f32
    };
    let (f10, f20, fin) = (share(&|i| d[i] * CELL <= 10.0), share(&|i| d[i] * CELL <= 20.0), share(&|i| d[i] * CELL > 20.0));
    let r = f20 / fin.max(0.05);
    (f10, f20, fin, r, f20 >= 0.40 && r >= 2.0)
}

/// C: the drowned valleys at 2 048² (the declared detector): (count, lengths p50 / p90 / max km, Spearman(length, area), the
/// five longest (length km, area km²)).
fn drowned_valleys(g: &GridF32, delta: f32) -> (usize, f32, f32, f32, f32, Vec<(f32, f32)>) {
    let n = 2048;
    let ocean = ocean_mask_t(g, Some([0, 0]));
    let dmask: Vec<bool> = (0..n * n).map(|k| ocean[k] && g.data[k] > -delta).collect();
    let low: Vec<bool> = (0..n * n).map(|k| ocean[k] && g.data[k] <= -delta).collect();
    let hy = hydrology(g, DOMAIN_KM / n as f32);
    let mut seen = vec![false; n * n];
    let mut out: Vec<(f32, f32)> = Vec::new();
    let ckm = DOMAIN_KM / n as f32;
    for s in 0..n * n {
        if !dmask[s] || seen[s] {
            continue;
        }
        // the component
        let mut comp = vec![s];
        seen[s] = true;
        let mut i = 0;
        while i < comp.len() {
            let k = comp[i];
            i += 1;
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if dmask[j] && !seen[j] {
                        seen[j] = true;
                        comp.push(j);
                    }
                }
            }
        }
        // the geodesic BFS from the cells adjacent to the lowstand ocean
        let mut dist: std::collections::HashMap<usize, u32> = std::collections::HashMap::new();
        let mut q = std::collections::VecDeque::new();
        let mut area = 0f32;
        for &k in &comp {
            let (x, y) = (k % n, k / n);
            let mut by_low = false;
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if low[j] {
                        by_low = true;
                    }
                    if !ocean[j] {
                        area = area.max(hy.acc_km2[j]);
                    }
                }
            }
            if by_low {
                dist.insert(k, 0);
                q.push_back(k);
            }
        }
        let inset: std::collections::HashSet<usize> = comp.iter().copied().collect();
        let mut maxd = 0u32;
        while let Some(k) = q.pop_front() {
            let dk = dist[&k];
            maxd = maxd.max(dk);
            let (x, y) = (k % n, k / n);
            for (dx, dy) in [(1usize, 0usize), (n - 1, 0), (0, 1), (0, n - 1)] {
                let j = ((y + dy) % n) * n + (x + dx) % n;
                if inset.contains(&j) && !dist.contains_key(&j) {
                    dist.insert(j, dk + 1);
                    q.push_back(j);
                }
            }
        }
        let len = maxd as f32 * ckm;
        if len >= 2.0 {
            out.push((len, area));
        }
    }
    let ls = sorted(out.iter().map(|o| o.0).collect());
    let sp = if out.len() >= 3 { spearman(&out.iter().map(|o| o.0).collect::<Vec<_>>(), &out.iter().map(|o| o.1).collect::<Vec<_>>()) } else { f32::NAN };
    out.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    (out.len(), pct(&ls, 0.5), pct(&ls, 0.9), ls.last().copied().unwrap_or(0.0), sp, out.into_iter().take(5).collect())
}

/// The land components ≥ 1 km² (periodic 8-connected) of a mask at 2 048².
fn islands(land: &[bool]) -> usize {
    let (_, comps) = periodic_components(land, 2048);
    let min_cells = (1.0 / (DOMAIN_KM / 2048.0).powi(2)).ceil() as usize;
    comps.iter().filter(|c| c.0 >= min_cells).count()
}

/// B4: the shelf -- for each coastal sea cell the distance to the nearest cell ≤ −Δ (km): (median, p90), and the mean depth
/// by the distance from the coast (0–20 km, 1 km bins).
fn shelf(g: &GridF32, delta: f32) -> (f32, f32, Vec<f32>) {
    let n = 2048;
    let ckm = DOMAIN_KM / n as f32;
    let ocean = ocean_mask_t(g, Some([0, 0]));
    let land: Vec<bool> = ocean.iter().map(|o| !o).collect();
    // distance to the cells ≤ −Δ (EDT of their complement), and to the land
    let notdeep: Vec<bool> = (0..n * n).map(|k| !(ocean[k] && g.data[k] <= -delta)).collect();
    let dd = edt(&notdeep, n, n);
    let dl = edt(&ocean, n, n);
    let mut w = Vec::new();
    for k in 0..n * n {
        if !ocean[k] {
            continue;
        }
        let (x, y) = (k % n, k / n);
        let nb = [y * n + (x + 1) % n, y * n + (x + n - 1) % n, ((y + 1) % n) * n + x, ((y + n - 1) % n) * n + x];
        if nb.iter().any(|&j| land[j]) {
            w.push(dd[k] * ckm);
        }
    }
    let w = sorted(w);
    let mut prof = vec![(0f64, 0usize); 21];
    for k in 0..n * n {
        if ocean[k] {
            let b = (dl[k] * ckm).floor() as usize;
            if b <= 20 {
                prof[b].0 += g.data[k] as f64;
                prof[b].1 += 1;
            }
        }
    }
    (pct(&w, 0.5), pct(&w, 0.9), prof.iter().map(|&(s, c)| if c > 0 { (s / c as f64) as f32 } else { f32::NAN }).collect())
}


// ════════════════════════════════════════ F173 ════════════════════════════════════════

const PANEL: [u64; 5] = [8, 18, 19, 23, 29];

fn load_curve() -> SeaCurve {
    let p = root().join("data/sealevel/spratt2016.txt");
    let text = std::fs::read_to_string(&p).expect("data/sealevel/spratt2016.txt (NOAA study 19982)");
    SeaCurve::from_noaa_spratt2016(&text).unwrap()
}
fn hist_curve(snaps: &[GridF32], cfg: &CascadeConfig, age: Option<f64>, rock: &Rock, curve: &SeaCurve) -> (Run, Vec<f32>) {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: None };
    let l = physics_history_sea(snaps, 256, cfg, &run, Some(rock), Some(curve), &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, Some([0, 0]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    eprintln!("      (curve: {} physics steps)", l.rec.steps);
    (Run { zm, d5, classes, m, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }, l.rec.e_cum_m)
}
/// C2: (W10, W_in, rampart_W) on a 256² field.
fn wall_index(zm: &GridF32) -> (f32, f32, bool) {
    let land = land_t(zm);
    let d = edt_periodic(&land, 256);
    let sl = slope_deg(&zm.data, 256, 256, CELL);
    let share = |f: &dyn Fn(usize) -> bool| -> f32 {
        let c: Vec<usize> = (0..65536).filter(|&i| land[i] && f(i)).collect();
        c.iter().filter(|&&i| sl[i] >= 8.0).count() as f32 / c.len().max(1) as f32
    };
    let (w10, win) = (share(&|i| d[i] * CELL <= 10.0), share(&|i| d[i] * CELL > 20.0));
    (w10, win, w10 >= 0.30 && w10 >= 2.0 * win)
}
/// C2's 2 048² reference: the share of the land within 10 km of the coast (the 256² EDT) whose slope exceeds 28°.
fn walls_coast_2048(zm256: &GridF32, g: &GridF32) -> f32 {
    let land = land_t(zm256);
    let d = edt_periodic(&land, 256);
    let lg = land_t(g);
    let sl = slope_deg(&g.data, 2048, 2048, DOMAIN_KM / 2048.0);
    let (mut n, mut w) = (0usize, 0usize);
    for k in 0..2048 * 2048 {
        let c = (k / 2048 / 8) * 256 + (k % 2048) / 8;
        if lg[k] && d[c] * CELL <= 10.0 {
            n += 1;
            if sl[k] > 28.0 {
                w += 1;
            }
        }
    }
    w as f32 / n.max(1) as f32
}
fn pearson(a: &[f32], b: &[f32]) -> f32 {
    let n = a.len() as f32;
    let (ma, mb) = (a.iter().sum::<f32>() / n, b.iter().sum::<f32>() / n);
    let cov: f32 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let va: f32 = a.iter().map(|x| (x - ma).powi(2)).sum();
    let vb: f32 = b.iter().map(|y| (y - mb).powi(2)).sum();
    cov / (va * vb).sqrt().max(1e-12)
}

/// C1: the coastal prism's history -- for today's land within 2 cells of the coast, Σ Δs_DS per step attributed by the
/// boundary within Chebyshev 2 at that step (C-O, C-C, none) and by the plate count at that step (≥ 10, 7–9, ≤ 6).
fn coastal_prism_history(seed: u64, p: C1Profile, end_land: &[bool]) -> ([f64; 3], [f64; 3]) {
    let n = 64;
    let dist = edt_periodic(end_land, 64);
    let coastal: Vec<bool> = (0..4096).map(|k| end_land[k] && dist[k] <= 2.0).collect();
    let cfg = c1_run_config();
    let mut st = p.init_state(64, seed);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let cur = std::cell::RefCell::new(vec![0f64; 4096]);
    let by_type = std::cell::RefCell::new([0f64; 3]);
    let by_plates = std::cell::RefCell::new([0f64; 3]);
    let mut on_term = |t: C1Term, a: &[f64], b: &[f64]| {
        if t == C1Term::DavisSuppe {
            let mut c = cur.borrow_mut();
            for k in 0..4096 {
                c[k] += b[k] - a[k];
            }
        }
    };
    run_with_closures_terms(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        |_, s, kn| {
            let ds = std::mem::replace(&mut *cur.borrow_mut(), vec![0f64; 4096]);
            let mut bi = classify_boundaries(&s.plate_id, kn);
            retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
            let plates = s.plate_id.data().iter().copied().collect::<std::collections::BTreeSet<u16>>().len();
            let pb = if plates >= 10 { 0 } else if plates >= 7 { 1 } else { 2 };
            let oceanic = |i: usize, j: usize| s.plate_type.get(i, j) == PlateType::Oceanic;
            let mut bt = by_type.borrow_mut();
            let mut bp = by_plates.borrow_mut();
            for k in 0..4096 {
                if !coastal[k] || ds[k] == 0.0 {
                    continue;
                }
                let (i, j) = (k % n, k / n);
                let mut ty = 2usize;
                for dy in 0..5usize {
                    for dx in 0..5usize {
                        let (x, y) = ((i + n + dx - 2) % n, (j + n + dy - 2) % n);
                        if bi.boundary_type.get(x, y) == BoundaryType::Convergent {
                            let nb = [((x + 1) % n, y), ((x + n - 1) % n, y), (x, (y + 1) % n), (x, (y + n - 1) % n)];
                            let co = oceanic(x, y) || nb.iter().any(|&(a, b)| oceanic(a, b));
                            ty = ty.min(if co { 0 } else { 1 });
                        }
                    }
                }
                bt[ty] += ds[k];
                bp[pb] += ds[k];
            }
        },
        &mut on_term,
    );
    (by_type.into_inner(), by_plates.into_inner())
}

// ════════════════════════════════════════ F174 ════════════════════════════════════════

const AGES174: [(&str, Option<f64>); 2] = [("B0", None), ("old 300 Myr", Some(300.0e6))];
const WINDOWS: [(&str, bool); 8] = [
    ("galicia_rias_baixas", true),
    ("dalmatia_split", true),
    ("bretagne_nord", true),
    ("corsica_west", true),
    ("landes_mimizan", false),
    ("vendee_sables", false),
    ("jutland_west", false),
    ("golfe_du_lion_agde", false),
];
const WCELL_KM: f32 = 50.0 / 256.0;

fn out174() -> PathBuf {
    root().join("docs/reports/relief_method/f174_margin/images")
}

/// The 4-connected sea component of the largest size in a window (`true` = sea), no wrap.
fn largest_sea_4(sea: &[bool], w: usize, h: usize) -> Vec<bool> {
    let mut lab = vec![u32::MAX; w * h];
    let mut sizes = Vec::new();
    let mut stack = Vec::new();
    for s in 0..w * h {
        if !sea[s] || lab[s] != u32::MAX {
            continue;
        }
        let id = sizes.len() as u32;
        lab[s] = id;
        stack.push(s);
        let mut c = 0usize;
        while let Some(k) = stack.pop() {
            c += 1;
            let (x, y) = (k % w, k / w);
            let mut push = |j: usize| {
                if sea[j] && lab[j] == u32::MAX {
                    lab[j] = id;
                    stack.push(j);
                }
            };
            if x > 0 {
                push(k - 1);
            }
            if x + 1 < w {
                push(k + 1);
            }
            if y > 0 {
                push(k - w);
            }
            if y + 1 < h {
                push(k + w);
            }
        }
        sizes.push(c);
    }
    let Some(best) = (0..sizes.len()).max_by_key(|&i| (sizes[i], std::cmp::Reverse(i))) else { return vec![false; w * h] };
    lab.iter().map(|&l| l == best as u32).collect()
}

/// A2: `physio::coast`'s four measures on a WINDOW: the window's border is never a coast (only the land/sea edges between two
/// in-window cells count, a coastline cell needs an in-window sea neighbour, the crack contours stop at the border).
fn coast_window(land: &[bool], w: usize, h: usize, block: usize) -> Coast {
    let (wi, hi) = (w as i64, h as i64);
    let inside = |x: i64, y: i64| x >= 0 && y >= 0 && x < wi && y < hi;
    let is_land = |x: i64, y: i64| land[y as usize * w + x as usize];
    let sea_nb = |x: i64, y: i64| inside(x, y) && !is_land(x, y);
    let cl: Vec<bool> = (0..w * h)
        .map(|k| {
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            land[k] && (sea_nb(x - 1, y) || sea_nb(x + 1, y) || sea_nb(x, y - 1) || sea_nb(x, y + 1))
        })
        .collect();
    let mut pts = Vec::new();
    for s in [2usize, 4, 8, 16, 32, 64] {
        let (bw, bh) = (w.div_ceil(s), h.div_ceil(s));
        let mut hit = vec![false; bw * bh];
        for k in 0..w * h {
            if cl[k] {
                hit[(k / w / s) * bw + (k % w) / s] = true;
            }
        }
        pts.push((s as f64, hit.iter().filter(|&&b| b).count() as f64));
    }
    let lx: Vec<f64> = pts.iter().map(|p| (1.0 / p.0).ln()).collect();
    let ly: Vec<f64> = pts.iter().map(|p| p.1.max(1.0).ln()).collect();
    let (mx, my) = (lx.iter().sum::<f64>() / 6.0, ly.iter().sum::<f64>() / 6.0);
    let num: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let den: f64 = lx.iter().map(|a| (a - mx).powi(2)).sum();
    let length_ratio = (pts[0].1 * pts[0].0) / (pts[4].1 * pts[4].0).max(1.0);
    use std::collections::HashMap;
    let mut out_edges: HashMap<(i64, i64), Vec<(i64, i64)>> = HashMap::new();
    let (mut edges, mut on_lattice) = (0usize, 0usize);
    let b = block as i64;
    for y in 0..hi {
        for x in 0..wi {
            if !is_land(x, y) {
                continue;
            }
            let mut add = |a: (i64, i64), c: (i64, i64), lattice: bool| {
                out_edges.entry(a).or_default().push(c);
                edges += 1;
                if lattice {
                    on_lattice += 1;
                }
            };
            if sea_nb(x, y - 1) {
                add((x, y), (x + 1, y), y % b == 0);
            }
            if sea_nb(x + 1, y) {
                add((x + 1, y), (x + 1, y + 1), (x + 1) % b == 0);
            }
            if sea_nb(x, y + 1) {
                add((x + 1, y + 1), (x, y + 1), (y + 1) % b == 0);
            }
            if sea_nb(x - 1, y) {
                add((x, y + 1), (x, y), x % b == 0);
            }
        }
    }
    let (mut all, mut al) = (0f64, 0f64);
    // the open contours (cut by the border) are followed from their heads (no incoming edge) first, then the loops
    let mut indeg: HashMap<(i64, i64), usize> = HashMap::new();
    for v in out_edges.values() {
        for &c in v {
            *indeg.entry(c).or_default() += 1;
        }
    }
    let mut starts: Vec<(i64, i64)> = out_edges.keys().copied().filter(|k| !indeg.contains_key(k)).collect();
    starts.sort();
    let mut rest: Vec<(i64, i64)> = out_edges.keys().copied().collect();
    rest.sort();
    starts.extend(rest);
    for s in starts {
        while let Some(first) = out_edges.get_mut(&s).and_then(|v| v.pop()) {
            let mut path = vec![s, first];
            let mut cur = first;
            let mut prev_dir = (first.0 - s.0, first.1 - s.1);
            while cur != s {
                let Some(v) = out_edges.get_mut(&cur) else { break };
                if v.is_empty() {
                    break;
                }
                let left = (-prev_dir.1, prev_dir.0);
                let pick = v.iter().position(|&n| (n.0 - cur.0, n.1 - cur.1) == left).unwrap_or(0);
                let nxt = v.swap_remove(pick);
                prev_dir = (nxt.0 - cur.0, nxt.1 - cur.1);
                path.push(nxt);
                cur = nxt;
            }
            if path.len() < 33 {
                continue;
            }
            let mut i = 0;
            while i + 8 < path.len() {
                let (dx, dy) = ((path[i + 8].0 - path[i].0) as f32, (path[i + 8].1 - path[i].1) as f32);
                let len = (dx * dx + dy * dy).sqrt() as f64;
                if len > 0.0 {
                    let phi = dy.abs().atan2(dx.abs()).to_degrees() % 45.0;
                    all += len;
                    if phi.min(45.0 - phi) <= 5.0 {
                        al += len;
                    }
                }
                i += 8;
            }
        }
    }
    Coast {
        dimension: (num / den.max(1e-12)) as f32,
        length_ratio: length_ratio as f32,
        aligned: if all > 0.0 { (al / all) as f32 } else { f32::NAN },
        block_index: if edges > 0 { (on_lattice as f32 / edges as f32) * block as f32 } else { f32::NAN },
        edges,
    }
}

/// The coastline cells (land with a 4-neighbour sea) of a window, no wrap, the border never a coast.
fn coast_cells_window(land: &[bool], w: usize, h: usize) -> Vec<bool> {
    (0..w * h)
        .map(|k| {
            let (x, y) = (k % w, k / w);
            land[k] && ((x > 0 && !land[k - 1]) || (x + 1 < w && !land[k + 1]) || (y > 0 && !land[k - w]) || (y + 1 < h && !land[k + w]))
        })
        .collect()
}
/// The coastline cells on the torus.
fn coast_cells_t(land: &[bool], n: usize) -> Vec<bool> {
    (0..n * n)
        .map(|k| {
            let (x, y) = (k % n, k / n);
            land[k] && (!land[y * n + (x + 1) % n] || !land[y * n + (x + n - 1) % n] || !land[((y + 1) % n) * n + x] || !land[((y + n - 1) % n) * n + x])
        })
        .collect()
}

/// The disk offsets within `r` cells.
fn disk(r: f32) -> Vec<(i64, i64)> {
    let ri = r.ceil() as i64;
    let mut v = Vec::new();
    for dy in -ri..=ri {
        for dx in -ri..=ri {
            if ((dx * dx + dy * dy) as f32) <= r * r {
                v.push((dx, dy));
            }
        }
    }
    v
}
/// A2's rule: for each listed coastline cell, the p90 of the land altitude within 5 km; `true` = rocky (> 200 m). `per` wraps
/// (our worlds' torus); otherwise only the in-window cells count.
fn rocky_cells(z: &[f32], land: &[bool], w: usize, h: usize, cells: &[usize], cell_km: f32, per: bool) -> Vec<bool> {
    let off = disk(5.0 / cell_km);
    let (wi, hi) = (w as i64, h as i64);
    cells
        .iter()
        .map(|&k| {
            let (x, y) = ((k % w) as i64, (k / w) as i64);
            let mut v: Vec<f32> = Vec::with_capacity(off.len());
            for &(dx, dy) in &off {
                let (mut a, mut b) = (x + dx, y + dy);
                if per {
                    a = a.rem_euclid(wi);
                    b = b.rem_euclid(hi);
                } else if a < 0 || b < 0 || a >= wi || b >= hi {
                    continue;
                }
                let j = b as usize * w + a as usize;
                if land[j] {
                    v.push(z[j]);
                }
            }
            if v.is_empty() {
                return false;
            }
            let i = ((v.len() - 1) as f32 * 0.9).round() as usize;
            let (_, p90, _) = v.select_nth_unstable_by(i, |a, b| a.partial_cmp(b).unwrap());
            *p90 > 200.0
        })
        .collect()
}

/// A real window: (altitude m, land = not the largest 4-connected sea component).
fn load_window(name: &str) -> (Vec<f32>, Vec<bool>) {
    let g = load_grid(&root().join(format!("data/coasts/{name}_256.bin"))).unwrap();
    assert_eq!((g.width, g.height), (256, 256));
    let sea: Vec<bool> = g.data.iter().map(|&v| v <= 0.0).collect();
    let ocean = largest_sea_4(&sea, 256, 256);
    (g.data.clone(), ocean.iter().map(|o| !o).collect())
}

#[derive(Clone, Debug, Default)]
struct Band {
    lo: [f32; 4],
    hi: [f32; 4],
}
fn meas(c: &Coast) -> [f32; 4] {
    [c.dimension, c.length_ratio, c.aligned, c.block_index]
}
fn in_band(b: &Band, m: &[f32; 4]) -> usize {
    (0..4).filter(|&i| m[i] >= b.lo[i] && m[i] <= b.hi[i]).count()
}
fn fm(m: &[f32; 4]) -> String {
    format!("{:.3} / {:.3} / {:.3} / {:.3}", m[0], m[1], m[2], m[3])
}

/// A2 on one of our worlds (2 048², m, sea at 0): the 64 windows of 256², eligible ≥ 128 coastline cells and a type by the
/// rule (≥ 80 %); per type the median of each measure, and the count of windows.
fn world_windows(g: &GridF32) -> ([Option<[f32; 4]>; 2], [usize; 3]) {
    let n = 2048;
    let land = land_t(g);
    let cl = coast_cells_t(&land, n);
    let mut per_type: [Vec<[f32; 4]>; 2] = [Vec::new(), Vec::new()];
    let mut counts = [0usize; 3];
    for wy in 0..8 {
        for wx in 0..8 {
            let idx: Vec<usize> = (0..256 * 256).map(|k| (wy * 256 + k / 256) * n + wx * 256 + k % 256).collect();
            let wl: Vec<bool> = idx.iter().map(|&k| land[k]).collect();
            let wcl = coast_cells_window(&wl, 256, 256);
            let cells: Vec<usize> = (0..256 * 256).filter(|&k| wcl[k]).map(|k| idx[k]).collect();
            if cells.len() < 128 {
                continue;
            }
            let _ = &cl;
            let rk = rocky_cells(&g.data, &land, n, n, &cells, WCELL_KM, true);
            let share = rk.iter().filter(|&&b| b).count() as f32 / rk.len() as f32;
            let ty = if share >= 0.8 { 0 } else if share <= 0.2 { 1 } else { 2 };
            counts[ty] += 1;
            if ty < 2 {
                per_type[ty].push(meas(&coast_window(&wl, 256, 256, 32)));
            }
        }
    }
    let med = |v: &Vec<[f32; 4]>| -> Option<[f32; 4]> {
        if v.is_empty() {
            return None;
        }
        let mut o = [0f32; 4];
        for i in 0..4 {
            let mut c: Vec<f32> = v.iter().map(|m| m[i]).filter(|x| x.is_finite()).collect();
            c.sort_by(|a, b| a.partial_cmp(b).unwrap());
            o[i] = if c.is_empty() { f32::NAN } else { c[c.len() / 2] };
        }
        Some(o)
    };
    ([med(&per_type[0]), med(&per_type[1])], counts)
}

/// A1: the final coastline's distance (cells of 256²) to the C1 contour, both on an n² torus: (≤ 1, ≤ 2, ≤ 5 shares, p50, p90).
fn inherited(final_land: &[bool], c1_land: &[bool], n: usize) -> (f32, f32, f32, f32, f32, Vec<f32>) {
    let fc = coast_cells_t(final_land, n);
    let cc = coast_cells_t(c1_land, n);
    let not_c1: Vec<bool> = cc.iter().map(|c| !c).collect();
    let d = edt_periodic(&not_c1, n);
    let f = (n / 256) as f32;
    let mut v: Vec<f32> = (0..n * n).filter(|&k| fc[k]).map(|k| d[k] / f).collect();
    let tot = v.len().max(1) as f32;
    let s = |t: f32| v.iter().filter(|&&x| x <= t).count() as f32 / tot;
    let (a, b, c) = (s(1.0), s(2.0), s(5.0));
    let dmap: Vec<f32> = (0..n * n).map(|k| if fc[k] { d[k] / f } else { -1.0 }).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p).round() as usize] };
    (a, b, c, q(0.5), q(0.9), dmap)
}

/// The land of a normalised field (the cascade's sea at 0.5) on the torus.
fn land_norm_t(z: &GridF32) -> Vec<bool> {
    land_t(&shifted(z, -0.5))
}

/// A1's image: north up, a 2 048² field downsampled ×2: land grey, the C1 contour black, the final coast coloured by its
/// distance (green ≤ 1, yellow ≤ 2, orange ≤ 5, red > 5 cells of 256²).
fn a1_png(final_land: &[bool], c1_land: &[bool], dmap: &[f32], n: usize, path: &std::path::Path) {
    let cc = coast_cells_t(c1_land, n);
    let m = n / 2;
    let mut img = image::RgbaImage::new(m as u32, m as u32);
    for y in 0..m {
        for x in 0..m {
            let mut px = [40u8, 70, 140, 255];
            let mut best: Option<f32> = None;
            let mut c1 = false;
            let mut land = 0;
            for dy in 0..2 {
                for dx in 0..2 {
                    let k = (2 * y + dy) * n + 2 * x + dx;
                    if final_land[k] {
                        land += 1;
                    }
                    if cc[k] {
                        c1 = true;
                    }
                    if dmap[k] >= 0.0 {
                        best = Some(best.map_or(dmap[k], |b: f32| b.max(dmap[k])));
                    }
                }
            }
            if land >= 2 {
                px = [200, 200, 190, 255];
            }
            if c1 {
                px = [0, 0, 0, 255];
            }
            if let Some(d) = best {
                px = if d <= 1.0 { [30, 170, 60, 255] } else if d <= 2.0 { [230, 210, 40, 255] } else if d <= 5.0 { [240, 130, 20, 255] } else { [220, 30, 30, 255] };
            }
            img.put_pixel(x as u32, (m - 1 - y) as u32, image::Rgba(px));
        }
    }
    img.save(path).unwrap();
}

/// The physics level under K2 with `steps` steps (A3).
fn hist_k2n(snaps: &[GridF32], cfg: &CascadeConfig, age: Option<f64>, rock: &Rock, steps: usize) -> (Run, Vec<f32>) {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: None };
    let l = physics_history_rock(snaps, 256, cfg, &run, Some(rock), &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, Some([0, 0]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    (Run { zm, d5, classes, m, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }, l.rec.e_cum_m)
}

/// A3: the share of `a`'s order-≥ 3 cells with no order-≥ 3 cell of `b` within its 3 × 3 (periodic).
fn displaced(a: &[u8], b: &[u8], n: usize) -> f32 {
    let mut tot = 0usize;
    let mut moved = 0usize;
    for k in 0..n * n {
        if a[k] < 3 {
            continue;
        }
        tot += 1;
        let (x, y) = (k % n, k / n);
        let mut hit = false;
        for dy in [n - 1, 0, 1] {
            for dx in [n - 1, 0, 1] {
                if b[((y + dy) % n) * n + (x + dx) % n] >= 3 {
                    hit = true;
                }
            }
        }
        if !hit {
            moved += 1;
        }
    }
    moved as f32 / tot.max(1) as f32
}

/// A4: per C1 step, the perimeter's shares (convergent, divergent, transform, none) and the same-plate share.
fn perimeter_series(seed: u64, p: C1Profile) -> Vec<([f32; 4], f32, usize)> {
    let n = 64;
    let cfg = c1_run_config();
    let mut st = p.init_state(64, seed);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let series = std::cell::RefCell::new(Vec::new());
    run_with_closures_terms(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        |_, s, kn| {
            let mut bi = classify_boundaries(&s.plate_id, kn);
            retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
            let cont = |i: usize, j: usize| s.plate_type.get(i, j) == PlateType::Continental;
            let mut counts = [0usize; 4];
            let mut same = 0usize;
            for j in 0..n {
                for i in 0..n {
                    if !cont(i, j) {
                        continue;
                    }
                    let nb4 = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
                    let oc: Vec<(usize, usize)> = nb4.iter().copied().filter(|&(x, y)| !cont(x, y)).collect();
                    if oc.is_empty() {
                        continue;
                    }
                    if oc.iter().any(|&(x, y)| s.plate_id.get(x, y) == s.plate_id.get(i, j)) {
                        same += 1;
                    }
                    let mut best = 3usize;
                    for dy in [n - 1, 0, 1] {
                        for dx in [n - 1, 0, 1] {
                            let c = match bi.boundary_type.get((i + dx) % n, (j + dy) % n) {
                                BoundaryType::Convergent => 0,
                                BoundaryType::Divergent => 1,
                                BoundaryType::Internal => 3,
                                _ => 2,
                            };
                            best = best.min(c);
                        }
                    }
                    counts[best] += 1;
                }
            }
            let tot = counts.iter().sum::<usize>().max(1) as f32;
            let plates = s.plate_id.data().iter().copied().collect::<std::collections::BTreeSet<u16>>().len();
            series.borrow_mut().push(([counts[0] as f32 / tot, counts[1] as f32 / tot, counts[2] as f32 / tot, counts[3] as f32 / tot], same as f32 / tot, plates));
        },
        &mut |_, _, _| {},
    );
    series.into_inner()
}

/// The windows' montage (north up), real coasts: the land shaded by altitude, the sea blue.
fn window_png(z: &[f32], land: &[bool], path: &std::path::Path) {
    let mx = z.iter().copied().fold(1f32, f32::max);
    let mut img = image::RgbaImage::new(256, 256);
    for k in 0..65536 {
        let (x, y) = (k % 256, k / 256);
        let v = (z[k].max(0.0) / mx).sqrt();
        let px = if land[k] { [(60.0 + 180.0 * v) as u8, (140.0 + 100.0 * v) as u8, (60.0 + 160.0 * v) as u8, 255] } else { [30, 60, 140, 255] };
        img.put_pixel(x as u32, 255 - y as u32, image::Rgba(px));
    }
    img.save(path).unwrap();
}

#[test]
#[ignore]
fn f174_measure() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out174()).unwrap();
    eprintln!("\n==========  Finding 174 . the coastline inherited from C1, real reference coasts, the time step, the active perimeter  ==========");
    // ── the window variant's control: on a closed island off the border it is physio::coast exactly ──
    let c256 = load_grid(&root().join("data/corsica/corse_256.bin")).unwrap();
    let l256: Vec<bool> = c256.data.iter().map(|&v| v > 0.0).collect();
    let (cw, cc) = (coast_window(&l256, 256, 256, 32), coast(&l256, 256, 256, 32));
    eprintln!("\n   control: corse_256 (closed island) window variant {} · physio::coast {}", fm(&meas(&cw)), fm(&meas(&cc)));
    assert_eq!(format!("{:?}", meas(&cw)), format!("{:?}", meas(&cc)), "the window variant on a closed island is physio::coast");
    // ── A2: the real windows ──
    eprintln!("\n── A2: the real coast windows (50 km, 195.3 m; dimension / length ratio / aligned / block index; the rule's rocky share) ──");
    let mut bands = [Band { lo: [f32::INFINITY; 4], hi: [f32::NEG_INFINITY; 4] }, Band { lo: [f32::INFINITY; 4], hi: [f32::NEG_INFINITY; 4] }];
    let mut rule_ok = 0;
    let mut tiles = Vec::new();
    for (name, rocky) in WINDOWS {
        let (z, land) = load_window(name);
        let c = coast_window(&land, 256, 256, 32);
        let m = meas(&c);
        let cl = coast_cells_window(&land, 256, 256);
        let cells: Vec<usize> = (0..65536).filter(|&k| cl[k]).collect();
        let rk = rocky_cells(&z, &land, 256, 256, &cells, WCELL_KM, false);
        let share = rk.iter().filter(|&&b| b).count() as f32 / rk.len().max(1) as f32;
        let rule = if share >= 0.8 { "rocky" } else if share <= 0.2 { "low" } else { "mixed" };
        let declared = if rocky { "rocky" } else { "low" };
        if rule == declared {
            rule_ok += 1;
        }
        let b = &mut bands[if rocky { 0 } else { 1 }];
        for i in 0..4 {
            b.lo[i] = b.lo[i].min(m[i]);
            b.hi[i] = b.hi[i].max(m[i]);
        }
        eprintln!(
            "   {name:<22} declared {declared:<5} · {} · coastline cells {} ({:.0} km) · land {:.0} % · the rule: rocky share {:.0} % → {rule}",
            fm(&m), cells.len(), cells.len() as f32 * WCELL_KM, 100.0 * land.iter().filter(|&&b| b).count() as f32 / 65536.0, 100.0 * share
        );
        let p = out174().join(format!("a2_{name}.png"));
        window_png(&z, &land, &p);
        tiles.push(image::open(&p).unwrap().to_rgba8());
    }
    let mut mont = image::RgbaImage::from_pixel(4 * 260, 2 * 260, image::Rgba([255, 255, 255, 255]));
    for (i, t) in tiles.iter().enumerate() {
        image::imageops::overlay(&mut mont, t, (2 + (i % 4) * 260) as i64, (2 + (i / 4) * 260) as i64);
    }
    mont.save(out174().join("a2_windows.png")).unwrap();
    for (i, t) in ["rocky", "low"].iter().enumerate() {
        let b = &bands[i];
        eprintln!(
            "   the {t} band: dimension {:.3}–{:.3} · length ratio {:.3}–{:.3} · aligned {:.3}–{:.3} · block index {:.3}–{:.3}",
            b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2], b.lo[3], b.hi[3]
        );
    }
    let p2_low = bands[1].lo[0] >= 1.02 && bands[1].hi[0] <= 1.08;
    let p2_rocky = bands[0].lo[0] >= 1.10 && bands[0].hi[0] <= 1.20;
    eprintln!("   the rule agrees with the declared type on {rule_ok} / 8 windows · P2: the low dimensions in 1.02–1.08 {p2_low} · the rocky in 1.10–1.20 {p2_rocky}");

    // ── the worlds: A1, A2 (ours), A3 ──
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let kh = KHetero::declared();
    let mut p1_ok = [0usize; 2];
    let mut p3_dz = 0usize;
    let mut p3_cls = true;
    let mut a4_all: Vec<(u64, Vec<([f32; 4], f32, usize)>)> = Vec::new();
    for seed in PANEL {
        let (np, tried) = select_continent_default(seed);
        let Some(np) = np else {
            eprintln!("   seed {seed}: no continent, skipped");
            continue;
        };
        let s = &tried.last().unwrap().1;
        let p = s.candidates[s.chosen.unwrap()].profile();
        eprintln!("\n   ══ seed {seed} ({p:?}, N {np}) ══");
        let w = world(seed, p);
        let cover = cover_map(seed, p, COVER_M);
        let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
        let rock = Rock { rho_k: RHO_K, cover_m: cover.clone(), ln_k: Some(ln_k.clone()) };
        let cfg = base_cfg(seed, rough_a, Some([0, 0]));
        let amp = || AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
        let last = w.snaps.last().unwrap();
        let c1_256 = land_norm_t(&upsample_to(last, 256));
        let c1_2048 = land_norm_t(&upsample_to(last, 2048));
        let mut dz_max = 0f32;
        for (ai, (name, age)) in AGES174.iter().enumerate() {
            let (r0, _) = hist_k2n(&w.snaps, &cfg, *age, &rock, N_PHYS);
            let (g0, _, c0) = chain_low(&r0.zm, amp(), Some(ln_k.clone()), 0.0);
            // A1
            let fl = land_t(&g0);
            let (a, b, c, q50, q90, dmap) = inherited(&fl, &c1_2048, 2048);
            let (pa, pb, pc, pq50, pq90, _) = inherited(&land_t(&r0.zm), &c1_256, 256);
            if a >= 0.9 {
                p1_ok[ai] += 1;
            }
            eprintln!(
                "   ── {name}: A1 at 2 048²: the coast within 1 / 2 / 5 cells of 256² of C1's contour {:.1} / {:.1} / {:.1} % (p50 {:.2}, p90 {:.2} cells) · the physics level alone (256²): {:.1} / {:.1} / {:.1} % (p50 {:.2}, p90 {:.2})",
                100.0 * a, 100.0 * b, 100.0 * c, q50, q90, 100.0 * pa, 100.0 * pb, 100.0 * pc, pq50, pq90
            );
            a1_png(&fl, &c1_2048, &dmap, 2048, &out174().join(format!("a1_seed{seed}_{}.png", if ai == 0 { "b0" } else { "old" })));
            // A2 on our world
            let (med, counts) = world_windows(&g0);
            let kt = coast_torus(&fl, 2048, 2048, 32);
            let fmo = |o: &Option<[f32; 4]>, b: &Band| match o {
                Some(m) => format!("{} ({} / 4 in the band)", fm(m), in_band(b, m)),
                None => "none".into(),
            };
            eprintln!(
                "      A2 on our world: windows rocky / low / mixed {} / {} / {} · rocky median {} · low median {} · the whole coast (coast_torus) {}",
                counts[0], counts[1], counts[2], fmo(&med[0], &bands[0]), fmo(&med[1], &bands[1]), fm(&meas(&kt))
            );
            // A3
            let (r1, _) = hist_k2n(&w.snaps, &cfg, *age, &rock, 2 * N_PHYS);
            let mut dz: Vec<f32> = r0.zm.data.iter().zip(&r1.zm.data).map(|(a, b)| (a - b).abs()).collect();
            let lz: Vec<f32> = (0..65536).filter(|&k| r0.zm.data[k] > 0.0 || r1.zm.data[k] > 0.0).map(|k| dz[k]).collect();
            dz.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mut lz = lz;
            lz.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let q = |v: &Vec<f32>, p: f32| v[((v.len() - 1) as f32 * p).round() as usize];
            let mx = *dz.last().unwrap();
            dz_max = dz_max.max(mx);
            let mv = (0..4).map(|i| (r1.m.fractions[i] - r0.m.fractions[i]).abs() * 100.0).fold(0f32, f32::max);
            let changed = (0..65536).filter(|&k| r0.classes[k] != r1.classes[k]).count() as f32 / 655.36;
            if mv >= 1.0 {
                p3_cls = false;
            }
            let (h0, h1) = (hydrology(&r0.zm, CELL), hydrology(&r1.zm, CELL));
            let dis = 0.5 * (displaced(&h0.strahler, &h1.strahler, 256) + displaced(&h1.strahler, &h0.strahler, 256));
            eprintln!(
                "      A3 (300 vs 600 steps): |Δz| all p50 / p99 / max {:.2} / {:.1} / {:.0} m · land p50 / p99 {:.2} / {:.1} m · the classes' largest move {:.2} points, cells changed {:.2} % · order-≥ 3 channels displaced {:.1} % · physics CPU {:.0} → {:.0} s",
                q(&dz, 0.5), q(&dz, 0.99), mx, q(&lz, 0.5), q(&lz, 0.99), mv, changed, 100.0 * dis, r0.cost.1, r1.cost.1
            );
            eprintln!("      cost (CPU per world): C1 {:.0} · physics {:.0} · chain {:.0} s", w.cost.1, r0.cost.1, c0.1);
        }
        if dz_max > 200.0 {
            p3_dz += 1;
        }
        // A4
        let ser = perimeter_series(seed, p);
        let mean = |r: std::ops::Range<usize>| -> ([f32; 4], f32) {
            let r = r.start..r.end.min(ser.len());
            let k = r.len().max(1) as f32;
            let mut o = [0f32; 4];
            let mut sp = 0f32;
            for i in r {
                for c in 0..4 {
                    o[c] += ser[i].0[c] / k;
                }
                sp += ser[i].1 / k;
            }
            (o, sp)
        };
        let (m50, s50) = mean(0..50);
        let (ma, sa) = mean(0..ser.len());
        eprintln!(
            "      A4 the perimeter (convergent / divergent / transform / none): first 50 steps {:.1} / {:.1} / {:.1} / {:.1} % · all {} steps {:.1} / {:.1} / {:.1} / {:.1} % · same-plate share {:.2} / {:.2} %",
            100.0 * m50[0], 100.0 * m50[1], 100.0 * m50[2], 100.0 * m50[3], ser.len(), 100.0 * ma[0], 100.0 * ma[1], 100.0 * ma[2], 100.0 * ma[3], 100.0 * s50, 100.0 * sa
        );
        let every: Vec<String> = (0..ser.len()).step_by(25).chain(std::iter::once(ser.len() - 1)).map(|i| format!("{}:{:.0}/{:.0}/{:.0}/{:.0}({})", i, 100.0 * ser[i].0[0], 100.0 * ser[i].0[1], 100.0 * ser[i].0[2], 100.0 * ser[i].0[3], ser[i].2)).collect();
        eprintln!("      A4 the series (step: C/D/T/none %, (plates)): {}", every.join(" "));
        a4_all.push((seed, ser));
    }
    let conv50: Vec<f32> = a4_all.iter().map(|(_, s)| s.iter().take(50).map(|x| x.0[0]).sum::<f32>() / 50f32.min(s.len() as f32)).collect();
    let mean50 = conv50.iter().sum::<f32>() / conv50.len().max(1) as f32;
    eprintln!(
        "\n   R: P1 (≥ 90 % within 1 cell on ≥ 4 / 5) B0 {} seeds → {} · old {} seeds → {} · P3: max |Δz| > 200 m on {p3_dz} seeds (≥ 3: {}), the classes < 1 point everywhere {p3_cls} · P4: the convergent share over the first 50 steps, the panel's mean {:.1} % (> 60 %: {}), per seed {:?}",
        p1_ok[0], p1_ok[0] >= 4, p1_ok[1], p1_ok[1] >= 4, p3_dz >= 3, 100.0 * mean50, mean50 > 0.6,
        conv50.iter().map(|v| format!("{:.1}", 100.0 * v)).collect::<Vec<_>>()
    );
    eprintln!("   total {:.0} s", t0.elapsed().as_secs_f64());
}
