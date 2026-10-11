//! ADR Finding 172 -- the sea-level rise (Δ 134 m over the history's last 2.58 Myr), the rampart diagnostics, and a screening
//! of 33 seeds, with K2, 16 plates (falling back to 20 / 24), on the torus.
//! Declared in `docs/reports/relief_method/f172_sea_level/f172_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f172_cascade -- --ignored --exact f172_rise --nocapture
//!      cargo test -p ymir-core --release --test f172_cascade -- --ignored --exact f172_screen --nocapture

#![allow(dead_code)]

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, PSEED, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::restrict;
use ymir_core::cascade::amplify::{AmpConfig, Chain, KHetero, Variant, Warp, WarpPlace, ocean_mask, ocean_mask_t};
use ymir_core::cascade::history::{
    Airy, C1Sources, C1Term, HistoryRun, Mapping, SeaLow, T_F166_YEARS, Timing, physics_history, physics_history_rock, run_c1_sources_mapped,
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
    root().join("docs/reports/relief_method/f172_sea_level/images")
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

#[test]
#[ignore]
fn f172_rise() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 172 . the sea-level rise (Δ 134 m over the last 2.58 Myr), the rampart  ==========");
    let cg = load_grid(&root().join("data/corsica/corse_1024.bin")).unwrap();
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let cclass = up_u8(&classify(&c128.data, 128, 128, CELL), 128, 8);
    let cmtn: Vec<bool> = cclass.iter().map(|&c| c == MOUNTAIN).collect();
    let ck = 200.0 / 1024.0;
    let c_whole = tex(&cg, ck, &vec![false; 1024 * 1024]);
    let c_mtn = tex(&cg, ck, &cmtn.iter().map(|m| !m).collect::<Vec<_>>());
    let c_coast = coast(&land_w(&cg), 1024, 1024, 32);
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let ckm = DOMAIN_KM / 2048.0;
    let ctrl = facet_control(2048, ckm);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    eprintln!("   Corsica's coast {c_coast:?}; the bands: dimension 1.06–1.18 · length ratio 1.15–1.60 · aligned 0.22–0.32 · block index 0.58–0.98");
    let band = |c: &Coast| -> usize {
        (1.06..=1.18).contains(&c.dimension) as usize + (1.15..=1.60).contains(&c.length_ratio) as usize + (0.22..=0.32).contains(&c.aligned) as usize + (0.58..=0.98).contains(&c.block_index) as usize
    };
    let kh = KHetero::declared();
    let mut succ_worlds: Vec<(u64, usize)> = Vec::new();
    let mut dim_gain_any = false;
    let mut macro_moves: Vec<f32> = Vec::new();
    let mut rows_img: Vec<Vec<image::RgbaImage>> = Vec::new();
    let mut crops: Vec<Vec<image::RgbaImage>> = Vec::new();
    let mut shelf_series: Vec<(Vec<f32>, [u8; 3])> = Vec::new();
    let mut a2_imgs: Vec<image::RgbaImage> = Vec::new();
    let pal = [[200u8, 50, 40], [40, 90, 200], [40, 150, 60], [140, 140, 140]];
    for (si, seed) in SEEDS.iter().copied().enumerate() {
        let (np, tried) = select_continent_default(seed);
        let Some(np) = np else {
            eprintln!("\n   seed {seed}: no continent qualifies even at 24 plates, skipped");
            continue;
        };
        let s = &tried.last().unwrap().1;
        let p = s.candidates[s.chosen.unwrap()].profile();
        eprintln!("\n   ══ seed {seed} ({p:?}, N {np}) ══");
        let w = world(seed, p);
        let off = viz_offset(w.snaps.last().unwrap());
        let cover = cover_map(seed, p, COVER_M);
        let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
        let rock = Rock { rho_k: RHO_K, cover_m: cover.clone(), ln_k: Some(ln_k.clone()) };
        let cfg = base_cfg(seed, rough_a, Some([0, 0]));
        let amp = || AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
        // A1: the perimeter
        let (shares, dsb) = perimeter(&w);
        eprintln!(
            "      A1 the perimeter (64²): convergent {:.1} % · divergent {:.1} % · transform {:.1} % · no boundary {:.1} % · Σ Δs_DS by the distance to the coast (cells 0…7, ≥ 8): {}",
            100.0 * shares[0],
            100.0 * shares[1],
            100.0 * shares[2],
            100.0 * shares[3],
            dsb.iter().map(|(m, c)| format!("{m:.2} [{c}]")).collect::<Vec<_>>().join(" · ")
        );
        let mut img_row = Vec::new();
        for (ai, (name, age)) in AGES.iter().enumerate() {
            // F171 K2 (before) and F172 (after the rise)
            let (r0, e0) = hist_k2(&w.snaps, &cfg, *age, &rock, None);
            let (g0, d50, c0) = chain_low(&r0.zm, amp(), Some(ln_k.clone()), 0.0);
            let (r1, _e1) = hist_k2(&w.snaps, &cfg, *age, &rock, Some(SeaLow::declared()));
            let (g1, d51, c1) = chain_low(&r1.zm, amp(), Some(ln_k.clone()), DELTA_M);
            let (k0, k1) = (coast_torus(&land_t(&g0), 2048, 2048, 32), coast_torus(&land_t(&g1), 2048, 2048, 32));
            let closer = (k1.dimension - c_coast.dimension).abs() < (k0.dimension - c_coast.dimension).abs()
                && (k1.length_ratio - c_coast.length_ratio).abs() < (k0.length_ratio - c_coast.length_ratio).abs()
                && band(&k1) >= band(&k0);
            dim_gain_any |= k1.dimension - k0.dimension >= 0.01;
            let mv = (0..4).map(|i| (r1.m.fractions[i] - r0.m.fractions[i]).abs() * 100.0).fold(0f32, f32::max);
            macro_moves.push(mv);
            if closer && mv < 3.0 {
                succ_worlds.push((seed, ai));
            }
            let (lt0, lt1) = (land_t(&g0), land_t(&g1));
            let n0 = lt0.iter().filter(|&&b| b).count() as f32;
            let lost = (0..2048 * 2048).filter(|&i| lt0[i] && !lt1[i]).count() as f32 / n0.max(1.0);
            let low_land: Vec<bool> = g1.data.iter().map(|&v| v > -DELTA_M).collect();
            let (isl_low, isl_after) = (islands(&low_land), islands(&lt1));
            let dv = drowned_valleys(&g1, DELTA_M);
            let (sw50, sw90, prof) = shelf(&g1, DELTA_M);
            let inland_low = (0..65536).filter(|&i| r1.zm.data[i] <= -DELTA_M && !ocean_mask_t(&r1.zm, Some([0, 0]))[i]).count();
            eprintln!("   ── seed {seed} {name}:");
            print_macro("F171 K2 (before)", &r0.m);
            print_macro("F172 (after the rise)", &r1.m);
            eprintln!(
                "      the classes move by ≤ {mv:.1} points · peaks {:.0} → {:.0} m · the land lost against F171 {:.2} % · islands ≥ 1 km² (lowstand land → after) {isl_low} → {isl_after} · inland 256² cells ≤ −Δ {inland_low}",
                peak_m(&g0),
                peak_m(&g1),
                100.0 * lost
            );
            eprintln!(
                "      the coast (torus): before {:.3} / {:.3} / {:.3} / {:.3} ({} in the bands) → after {:.3} / {:.3} / {:.3} / {:.3} ({} in the bands) · Δdimension {:+.3} · Δlength ratio {:+.3} → closer to Corsica: {closer}",
                k0.dimension, k0.length_ratio, k0.aligned, k0.block_index, band(&k0), k1.dimension, k1.length_ratio, k1.aligned, k1.block_index, band(&k1), k1.dimension - k0.dimension, k1.length_ratio - k0.length_ratio
            );
            eprintln!(
                "      drowned valleys (≥ 2 km): {} · length p50 / p90 / max {:.1} / {:.1} / {:.1} km · Spearman(length, river area) {:.2} · the 5 longest (km, km²): {} · the shelf (0 to −Δ) median / p90 {:.2} / {:.2} km",
                dv.0,
                dv.1,
                dv.2,
                dv.3,
                dv.4,
                dv.5.iter().map(|(l, a)| format!("{l:.1} km / {a:.0} km²")).collect::<Vec<_>>().join(", "),
                sw50,
                sw90
            );
            shelf_series.push((prof, pal[si % 4]));
            // A2 (before, F171's world) and A3 at B0 and old
            let (f10, f20, fin, rr, rp) = rampart(&r0.zm, &r0.classes);
            eprintln!("      A2 (F171 K2): the mountain share within 10 km {:.1} % · 20 km {:.1} % · interior {:.1} % · R {rr:.2} → rampart {rp}", 100.0 * f10, 100.0 * f20, 100.0 * fin);
            if ai == 0 || ai == 3 {
                let land = land_t(&r0.zm);
                let dist = edt_periodic(&land, 256);
                let rk = rock_at_surface(&land, &cover, &e0);
                let m8 = up_u8(&r0.classes, 256, 8);
                let sl = slope_deg(&g0.data, 2048, 2048, ckm);
                let mut parts = Vec::new();
                for (bn, f) in [("0–10 km", &(|i: usize| dist[i] * CELL <= 10.0) as &dyn Fn(usize) -> bool), ("10–20 km", &|i: usize| dist[i] * CELL > 10.0 && dist[i] * CELL <= 20.0), ("> 20 km", &|i: usize| dist[i] * CELL > 20.0), ("cover", &|i: usize| rk[i] == 1), ("basement", &|i: usize| rk[i] == 2), ("shield", &|i: usize| rk[i] == 3)] {
                    let sel: Vec<bool> = (0..2048 * 2048).map(|k| { let c = (k / 2048 / 8) * 256 + (k % 2048) / 8; m8[k] == MOUNTAIN && f(c) && !d50[k] }).collect();
                    let nsel = sel.iter().filter(|&&b| b).count();
                    let walls = (0..2048 * 2048).filter(|&k| sel[k] && sl[k] > 28.0).count();
                    let fac = if nsel > 2000 { let excl: Vec<bool> = sel.iter().map(|b| !b).collect(); read_cell(&g0, Some(&excl), ckm).0.facets } else { f64::NAN };
                    parts.push(format!("{bn}: {:.0} km² of mountain, walls 28° {:.2} %, facets {:.2}", nsel as f32 * ckm * ckm, 100.0 * walls as f32 / nsel.max(1) as f32, fac));
                }
                eprintln!("      A3 (F171 K2, the mountain class): {}", parts.join(" · "));
            }
            // texture and regularity (after)
            let m8 = up_u8(&r1.classes, 256, 8);
            let mtn8: Vec<bool> = m8.iter().map(|&x| x == MOUNTAIN).collect();
            let excl: Vec<bool> = (0..2048 * 2048).map(|i| d51[i] || !mtn8[i]).collect();
            let wt = tex(&g1, ckm, &d51);
            let mt = tex(&g1, ckm, &excl);
            for (t, x, cc) in [("whole", &wt, &c_whole), ("mountain", &mt, &c_mtn)] {
                let el = elements(x, cc, ckm, ctrl_ok);
                let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
                eprintln!("      after, 2048² {t}: {} of {} within ×1.5 · out: {}", el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") });
            }
            let p1024 = { let gg = restrict(&g1); let dd = (0..1024 * 1024).map(|i| d51[(i / 1024 * 2) * 2048 + (i % 1024) * 2]).collect::<Vec<_>>(); predictability(&gg, 2.0 * ckm, Some(&dd)) };
            let hy = hydrology(&r1.zm, CELL);
            eprintln!(
                "      after: slope guard {:.3}/{:.3} · regularity 391 m CV_λ {:.3} · CV_L {:.3} · σ_θ {:.1}° · R2_g {:.3} · lakes {} / {:.0} km² · cost before {:.0} / after {:.0} s CPU per world",
                mt.r.slope_p50, mt.r.slope_p90, p1024.cv_lambda, p1024.cv_len, p1024.theta_std, p1024.r2_gully,
                hy.lakes.len(), hy.lake_id.iter().filter(|&&l| l != 0).count() as f32 * CELL * CELL,
                w.cost.1 + r0.cost.1 + c0.1, w.cost.1 + r1.cost.1 + c1.1
            );
            // images: the whole maps (young and old), the crops (old), the A2 map (young)
            let shade = |z: &GridF32| { let (_, hy) = read(z, None); whole(&rgba_of(&rolln(z, off), Some(&hy), 195.3, &ss), 2048, 280) };
            if ai == 1 || ai == 3 {
                img_row.push(shade(&g0));
                img_row.push(shade(&g1));
            }
            if ai == 3 {
                // the window around the longest drowned valley's mouth region: the 128² window with the most drowned cells
                let ocean = ocean_mask_t(&g1, Some([0, 0]));
                let dm: Vec<bool> = (0..2048 * 2048).map(|k| ocean[k] && g1.data[k] > -DELTA_M).collect();
                let mut best = (0usize, 0usize, 0usize);
                for oy in (0..2048 - 128).step_by(64) {
                    for ox in (0..2048 - 128).step_by(64) {
                        let c = (0..128 * 128).filter(|&i| dm[(oy + i / 128) * 2048 + ox + i % 128]).count();
                        if c > best.2 { best = (ox, oy, c) }
                    }
                }
                let crop = |g: &GridF32, nn: usize, x: usize, y: usize| { let (_, hy) = read(g, None); window(&rgba_of(g, Some(&hy), 195.3, &ss), nn, (x as f32 / nn as f32, y as f32 / nn as f32), 128.0 / nn as f32, 300) };
                let cc = { let land = land_w(&cg); let mut b = (0usize, 0usize); 'o: for oy in (0..1024 - 128).step_by(64) { for ox in (0..1024 - 128).step_by(64) { let c = (0..128 * 128).filter(|&i| land[(oy + i / 128) * 1024 + ox + i % 128]).count() as f32 / 16384.0; if (0.4..=0.6).contains(&c) { b = (ox, oy); break 'o; } } } b };
                crops.push(vec![crop(&cg, 1024, cc.0, cc.1), crop(&g0, 2048, best.0, best.1), crop(&g1, 2048, best.0, best.1)]);
            }
            if ai == 1 {
                let land = land_t(&r0.zm);
                let dist = edt_periodic(&land, 256);
                let mut img = image::RgbaImage::new(256, 256);
                for y in 0..256u32 {
                    for x in 0..256u32 {
                        let gy = 255 - y as usize;
                        let k = ((gy + off[1] * 4) % 256) * 256 + (x as usize + off[0] * 4) % 256;
                        let c = if !land[k] { [25u8, 35, 60] } else if r0.classes[k] == MOUNTAIN { if dist[k] * CELL <= 20.0 { [235, 70, 50] } else { [235, 235, 235] } } else if dist[k] * CELL <= 20.0 { [150, 120, 90] } else { [110, 150, 90] };
                        img.put_pixel(x, y, image::Rgba([c[0], c[1], c[2], 255]));
                    }
                }
                a2_imgs.push(img);
            }
        }
        rows_img.push(img_row);
        // the sensitivity: Δ = 60 m, the témoin, old
        if seed == PSEED {
            let sl60 = SeaLow { delta_m: 60.0, ..SeaLow::declared() };
            let (r0, _) = hist_k2(&w.snaps, &cfg, Some(300.0e6), &rock, None);
            let (g0, _, _) = chain_low(&r0.zm, amp(), Some(ln_k.clone()), 0.0);
            let (r1, _) = hist_k2(&w.snaps, &cfg, Some(300.0e6), &rock, Some(sl60));
            let (g1, _, _) = chain_low(&r1.zm, amp(), Some(ln_k.clone()), 60.0);
            let (k0, k1) = (coast_torus(&land_t(&g0), 2048, 2048, 32), coast_torus(&land_t(&g1), 2048, 2048, 32));
            let dv = drowned_valleys(&g1, 60.0);
            eprintln!(
                "   sensitivity Δ = 60 m (the témoin, old): the coast {:.3} / {:.3} → {:.3} / {:.3} (dimension / length ratio) · drowned valleys {} (max {:.1} km) · the classes move ≤ {:.1} points",
                k0.dimension, k0.length_ratio, k1.dimension, k1.length_ratio, dv.0, dv.3,
                (0..4).map(|i| (r1.m.fractions[i] - r0.m.fractions[i]).abs() * 100.0).fold(0f32, f32::max)
            );
        }
    }
    // the success rule: ≥ 3 seeds with ≥ 2 ages closer and < 3 points
    let mut per_seed: std::collections::BTreeMap<u64, usize> = std::collections::BTreeMap::new();
    for (s, _) in &succ_worlds {
        *per_seed.entry(*s).or_default() += 1;
    }
    let seeds_ok = per_seed.values().filter(|&&c| c >= 2).count();
    eprintln!(
        "\n   R: the worlds closer to Corsica with the classes within 3 points: {:?} · seeds with ≥ 2 such ages {seeds_ok} · the classes' largest move {:.1} points · the dimension rose by ≥ 0.01 somewhere: {dim_gain_any} → {}",
        succ_worlds,
        macro_moves.iter().copied().fold(0f32, f32::max),
        if seeds_ok >= 3 { "SUCCESS: the rise acts as a mechanism; the screening follows" } else if !dim_gain_any { "NO EFFECT: stop before the screening" } else { "no success as written" }
    );
    compose(&rows_img).save(out().join("f172_before_after.png")).unwrap();
    compose(&crops).save(out().join("f172_coast_25km.png")).unwrap();
    compose(&[a2_imgs]).save(out().join("f172_rampart_map.png")).unwrap();
    profile_png(&shelf_series, 1.0, &out().join("f172_shelf_profile.png"));
    eprintln!("   images: f172_before_after.png (rows = the seeds; young before | after, old before | after), f172_coast_25km.png (Corsica | F171 | F172, old, the window with the most drowned cells), f172_rampart_map.png (young: red = mountains within 20 km of the coast, white = interior mountains, brown = coastal band, green = interior), f172_shelf_profile.png (the mean depth 0–20 km from the coast, one line per world, coloured by seed) · north up · {NOTICE}");
    eprintln!("\n==========  end Finding 172 (rise) . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F172-D: the screening of 33 seeds (physics 256² only), young and old, K2 and the sea-level rise.
#[test]
#[ignore]
fn f172_screen() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    std::fs::create_dir_all(out()).unwrap();
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let kh = KHetero::declared();
    let mut seeds: Vec<u64> = (1..=32).collect();
    seeds.push(PSEED);
    eprintln!("\n==========  Finding 172 . D — the screening of {} seeds (physics 256², young and old, K2, the rise)  ==========", seeds.len());
    let mut tsv = String::from("seed\tage\tN\tk\tfamily\tland_pct\tr_star\tmountain\tplain_plateau\thill\tmedian_m\tpeak_m\tR\trampart\tconv_perimeter\tmtn_regions_100km2\tcover\tbasement\tshield\tlakes\tlake_km2\tcost_s\n");
    struct Row {
        age: usize,
        mtn: f32,
        low: f32,
        med: f32,
        peak: f32,
        r: f32,
        rampart: bool,
        conv: f32,
        regions: usize,
        europe: bool,
    }
    let mut rows: Vec<(u64, Vec<Row>)> = Vec::new();
    let mut thumbs: Vec<Vec<image::RgbaImage>> = Vec::new();
    let (mut cost_sel, mut cost_all, mut none) = (0f64, 0f64, 0usize);
    for seed in seeds.iter().copied() {
        let clk = Clock::start();
        let (np, tried) = select_continent_default(seed);
        let cs = clk.stop().1;
        cost_sel += cs;
        cost_all += cs;
        let Some(np) = np else {
            none += 1;
            eprintln!("   seed {seed}: no continent qualifies even at 24 plates (selection {cs:.0} s CPU)");
            tsv.push_str(&format!("{seed}\t-\tnone\n"));
            continue;
        };
        let s = &tried.last().unwrap().1;
        let cand = s.candidates[s.chosen.unwrap()];
        let p = cand.profile();
        let w = world(seed, p);
        let (conv_share, _) = { let (sh, b) = perimeter(&w); (sh[0], b) };
        let cover = cover_map(seed, p, COVER_M);
        let rock = Rock { rho_k: RHO_K, cover_m: cover.clone(), ln_k: Some(kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]))) };
        let cfg = base_cfg(seed, rough_a, Some([0, 0]));
        let off = viz_offset(w.snaps.last().unwrap());
        let mut rr = Vec::new();
        let mut th = Vec::new();
        for (ai, age) in [(1usize, Some(20.0e6)), (3, Some(300.0e6))] {
            let c2 = Clock::start();
            let (r, e) = hist_k2(&w.snaps, &cfg, age, &rock, Some(SeaLow::declared()));
            let cw = c2.stop().1 + w.cost.1;
            cost_all += cw;
            let (_, _, _, rv, rp) = rampart(&r.zm, &r.classes);
            let land = land_t(&r.zm);
            let rk = rock_at_surface(&land, &cover, &e);
            let nl = land.iter().filter(|&&b| b).count().max(1) as f32;
            let share = |v: u8| rk.iter().filter(|&&x| x == v).count() as f32 / nl;
            let mtn_mask: Vec<bool> = r.classes.iter().map(|&c| c == MOUNTAIN).collect();
            let (_, comps) = periodic_components(&mtn_mask, 256);
            let regions = comps.iter().filter(|c| c.0 as f32 * CELL * CELL > 100.0).count();
            let hy = hydrology(&r.zm, CELL);
            let m = &r.m;
            let low = m.fractions[0] + m.fractions[1];
            let europe = (0.165..=0.372).contains(&m.fractions[3]) && (0.375..=0.845).contains(&low) && (98.0..=390.0).contains(&m.deciles[4]);
            let lake_km2 = hy.lake_id.iter().filter(|&&l| l != 0).count() as f32 * CELL * CELL;
            tsv.push_str(&format!(
                "{seed}\t{}\t{np}\t{}\t{}\t{:.1}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.0}\t{:.0}\t{:.2}\t{}\t{:.2}\t{regions}\t{:.2}\t{:.2}\t{:.2}\t{}\t{lake_km2:.0}\t{:.0}\n",
                if ai == 1 { "young" } else { "old" }, cand.k, if cand.compact { "compact" } else { "BFS" }, 100.0 * cand.land, cand.circ.sea_r.unwrap_or(0),
                100.0 * m.fractions[3], 100.0 * low, 100.0 * m.fractions[2], m.deciles[4], peak_m(&r.zm), rv, rp, conv_share, share(1), share(2), share(3), hy.lakes.len(), cw
            ));
            eprintln!(
                "   seed {seed} {} (N {np}, k {}, {} r* {}): land {:.1} % · mountain {:.1} % · plain + plateau {:.1} % · median {:.0} m · peak {:.0} m · R {rv:.2} rampart {rp} · convergent perimeter {:.1} % · mountain regions > 100 km² {regions} · Europe's bands {europe}",
                if ai == 1 { "young" } else { "old" }, cand.k, if cand.compact { "compact" } else { "BFS" }, cand.circ.sea_r.unwrap_or(0), 100.0 * cand.land,
                100.0 * m.fractions[3], 100.0 * low, m.deciles[4], peak_m(&r.zm), 100.0 * conv_share
            );
            rr.push(Row { age: ai, mtn: m.fractions[3], low, med: m.deciles[4], peak: peak_m(&r.zm), r: rv, rampart: rp, conv: conv_share, regions, europe });
            let (_, hyd) = read(&r.zm, None);
            th.push(whole(&rgba_of(&rolln(&r.zm, off), Some(&hyd), 1562.5, &ss), 256, 128));
            th.push(class_png(&roll_u8(&r.classes, 256, [off[0] * 4, off[1] * 4]), 256, 256, 1));
        }
        thumbs.push(th.into_iter().map(|i| image::imageops::resize(&i, 128, 128, image::imageops::FilterType::Nearest)).collect());
        rows.push((seed, rr));
    }
    std::fs::write(out().join("../catalogue.tsv"), &tsv).unwrap();
    compose(&thumbs).save(out().join("f172_catalogue.png")).unwrap();
    // the distributions
    for (ai, an) in [(1usize, "young"), (3, "old")] {
        let pick = |f: &dyn Fn(&Row) -> f32| -> Vec<f32> { sorted(rows.iter().flat_map(|(_, r)| r.iter().filter(|x| x.age == ai).map(f)).collect()) };
        let dec = |v: &[f32]| (1..10).map(|d| format!("{:.1}", pct(v, d as f64 / 10.0))).collect::<Vec<_>>().join("/");
        for (nm, f) in [("mountain %", &(|r: &Row| 100.0 * r.mtn) as &dyn Fn(&Row) -> f32), ("plain + plateau %", &|r: &Row| 100.0 * r.low), ("median m", &|r: &Row| r.med), ("peak m", &|r: &Row| r.peak), ("R", &|r: &Row| r.r), ("mountain regions", &|r: &Row| r.regions as f32), ("convergent perimeter %", &|r: &Row| 100.0 * r.conv)] {
            let v = pick(f);
            eprintln!("   D {an} {nm}: median {:.1} · deciles {}", pct(&v, 0.5), dec(&v));
        }
        let n = rows.len().max(1) as f32;
        let ramp = rows.iter().filter(|(_, r)| r.iter().any(|x| x.age == ai && x.rampart)).count() as f32 / n;
        let (cv, rv): (Vec<f32>, Vec<f32>) = rows.iter().flat_map(|(_, r)| r.iter().filter(|x| x.age == ai).map(|x| (x.conv, x.r))).unzip();
        eprintln!("   D {an}: the rampart on {:.0} % of the seeds · Spearman(convergent perimeter, R) {:.2}", 100.0 * ramp, spearman(&cv, &rv));
    }
    let n = rows.len().max(1) as f32;
    let eu_any = rows.iter().filter(|(_, r)| r.iter().any(|x| x.europe)).count() as f32 / n;
    let ramp_young = rows.iter().filter(|(_, r)| r.iter().any(|x| x.age == 1 && x.rampart)).count() as f32 / n;
    eprintln!(
        "\n   D: {} seeds screened ({none} with no qualified continent) · within Europe's three bands at ≥ 1 age: {:.0} % · the rampart at young on {:.0} %{} · cost {:.0} s CPU in all ({:.0} per seed), the selection {:.0} %",
        seeds.len(),
        100.0 * eu_any,
        100.0 * ramp_young,
        if ramp_young > 0.8 { " → a PROBABLE C1 BIAS (> 80 %)" } else { "" },
        cost_all,
        cost_all / seeds.len() as f64,
        100.0 * cost_sel / cost_all.max(1e-9)
    );
    eprintln!("   the catalogue: docs/reports/relief_method/f172_sea_level/catalogue.tsv and images/f172_catalogue.png (one row per seed in the table's order; young shaded | young classes | old shaded | old classes)");
    eprintln!("\n==========  end Finding 172 (screen) . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}
