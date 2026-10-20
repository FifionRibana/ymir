//! ADR Finding 186 -- the deposition by closures (Yuan's K_eff at 256², the cover transmitted to the amplifications) and the
//! chain's 5-minute budget to 8 192². F185's helpers are kept (its tests disabled). Declared in
//! `docs/reports/relief_method/f186_deposition_closure/f186_declared.md`; resumable (`.../f186_deposition_closure/results/`).
//!
//! Run: cargo test -p ymir-core --release --test f186_cascade -- --ignored --exact f186_budget --nocapture
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

// ════════════════════════════════════════ F174 B/C ════════════════════════════════════════

use ymir_core::cascade::margin::{Hetero, Margin, Ramp};
use ymir_core::cascade::profile::select_continent_default_margin;

/// The variants of Part B (the main margin first after Part A's « off »).
fn variants() -> Vec<(&'static str, Margin)> {
    let m = Margin::declared();
    vec![
        ("A (off)", Margin::off()),
        ("margin", m),
        ("w_rim x0.5", Margin { w_rim: m.w_rim * 0.5, ..m }),
        ("w_rim x2", Margin { w_rim: m.w_rim * 2.0, ..m }),
        ("w_n x2", Margin { ramp: Some(Ramp::declared().scaled(2.0)), ..m }),
        ("w_n /2", Margin { ramp: Some(Ramp::declared().scaled(0.5)), ..m }),
        ("B-hetero", Margin { hetero: Some(Hetero::declared()), ..m }),
    ]
}

/// C1 under a profile and a margin, its Airy sources' record.
fn world_m(seed: u64, p: C1Profile, m: &Margin) -> World {
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let clk = Clock::start();
    let src = run_c1_sources_mapped(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &SteinSteinParams::default(), &[C1Term::Erosion], Mapping::Airy(Airy::declared()));
    let cost = clk.stop();
    let snaps = src.hist.normalized(None);
    World { src, st, kin, snaps, cost }
}
fn cover_map_m(seed: u64, p: C1Profile, m: &Margin, d_c: f32) -> Vec<f32> {
    let st0 = p.init_state_margin(64, seed, m);
    up_mask(st0.cratonic_mask.data(), 64, 4).iter().map(|&s| if s { 0.0 } else { d_c }).collect()
}

/// B1's verification, per C1 step: (active share, rim-consumed share, plates); and the last step's perimeter map
/// (0 none, 1 continental, 2 passive perimeter, 3 active perimeter).
fn margin_series(seed: u64, p: C1Profile, m: &Margin) -> (Vec<(f32, f32, usize)>, Vec<u8>, (f32, f32)) {
    let n = 64;
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let series = std::cell::RefCell::new(Vec::new());
    let last = std::cell::RefCell::new(vec![0u8; n * n]);
    let classify_step = |s: &C1State, kn: &PlateKinematics| -> ((f32, f32, usize), Vec<u8>) {
        let mut bi = classify_boundaries(&s.plate_id, kn);
        retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
        let cont = |i: usize, j: usize| s.plate_type.get(i, j) == PlateType::Continental;
        let mut map = vec![0u8; n * n];
        let (mut per, mut act, mut cons) = (0usize, 0usize, 0usize);
        for j in 0..n {
            for i in 0..n {
                if !cont(i, j) {
                    continue;
                }
                map[j * n + i] = 1;
                let nb4 = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
                let oc: Vec<(usize, usize)> = nb4.iter().copied().filter(|&(x, y)| !cont(x, y)).collect();
                if oc.is_empty() {
                    continue;
                }
                per += 1;
                if oc.iter().any(|&(x, y)| s.plate_id.get(x, y) != s.plate_id.get(i, j)) {
                    cons += 1;
                }
                let mut a = false;
                for dy in [n - 1, 0, 1] {
                    for dx in [n - 1, 0, 1] {
                        if bi.boundary_type.get((i + dx) % n, (j + dy) % n) == BoundaryType::Convergent {
                            a = true;
                        }
                    }
                }
                if a {
                    act += 1;
                }
                map[j * n + i] = if a { 3 } else { 2 };
            }
        }
        let plates = s.plate_id.data().iter().copied().collect::<std::collections::BTreeSet<u16>>().len();
        ((act as f32 / per.max(1) as f32, cons as f32 / per.max(1) as f32, plates), map)
    };
    let t0 = classify_step(&st, &kin).0;
    run_with_closures_terms(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        |_, s, kn| {
            let (v, map) = classify_step(s, kn);
            series.borrow_mut().push(v);
            *last.borrow_mut() = map;
        },
        &mut |_, _, _| {},
    );
    (series.into_inner(), last.into_inner(), (t0.0, t0.1))
}

/// The shelf (2 048², m, sea at 0): per coastline cell the steepest-descent path through the ocean; the break is the first
/// point ≥ 1 km out whose next-1-km mean slope is ≥ 3 × the mean slope since the coast. (median width km, p90 width, median
/// break depth m, share with a shelf).
fn shelf_break(g: &GridF32) -> (f32, f32, f32, f32) {
    let n = g.width;
    let ck = DOMAIN_KM / n as f32;
    let land = land_t(g);
    let cl = coast_cells_t(&land, n);
    let z = &g.data;
    let idx = |x: i64, y: i64| (y.rem_euclid(n as i64) as usize) * n + x.rem_euclid(n as i64) as usize;
    let mut widths = Vec::new();
    let mut depths = Vec::new();
    let mut with = 0usize;
    let mut tot = 0usize;
    for k in 0..n * n {
        if !cl[k] {
            continue;
        }
        tot += 1;
        let (mut x, mut y) = ((k % n) as i64, (k / n) as i64);
        let mut s = vec![0f32];
        let mut zz = vec![z[k].max(0.0)];
        let mut cur = f32::INFINITY; // the first step goes to the lowest ocean neighbour
        loop {
            let mut best: Option<(i64, i64, f32, f32)> = None;
            for (dx, dy) in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)] {
                let j = idx(x + dx, y + dy);
                if land[j] {
                    continue;
                }
                let v = z[j];
                if v < cur && best.is_none_or(|b| v < b.2) {
                    best = Some((x + dx, y + dy, v, if dx != 0 && dy != 0 { std::f32::consts::SQRT_2 } else { 1.0 }));
                }
            }
            let Some((nx, ny, v, l)) = best else { break };
            x = nx;
            y = ny;
            cur = v;
            s.push(s.last().unwrap() + l * ck);
            zz.push(v);
            if *s.last().unwrap() > 60.0 || v < -3000.0 {
                break;
            }
        }
        let mut found: Option<(f32, f32)> = None;
        let mut jj = 0usize;
        for i in 1..s.len() {
            if s[i] < 1.0 {
                continue;
            }
            jj = jj.max(i);
            while jj < s.len() && s[jj] < s[i] + 1.0 {
                jj += 1;
            }
            if jj >= s.len() {
                break;
            }
            let before = (zz[0] - zz[i]) / s[i];
            let after = (zz[i] - zz[jj]) / (s[jj] - s[i]);
            if after > 0.0 && after >= 3.0 * before {
                found = Some((s[i], -zz[i]));
                break;
            }
        }
        match found {
            Some((w, d)) => {
                with += 1;
                widths.push(w);
                depths.push(d);
            }
            None => widths.push(0.0),
        }
    }
    widths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |v: &Vec<f32>, p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p).round() as usize] };
    (q(&widths, 0.5), q(&widths, 0.9), q(&depths, 0.5), with as f32 / tot.max(1) as f32)
}

/// A1 with the distance field kept: (≤ 1, ≤ 2, ≤ 5 shares, p50, p90, the coast-cell distance map, the full distance field in
/// cells of 256²).
fn inherited_full(final_land: &[bool], c1_land: &[bool], n: usize) -> (f32, f32, f32, f32, f32, Vec<f32>, Vec<f32>) {
    let fc = coast_cells_t(final_land, n);
    let cc = coast_cells_t(c1_land, n);
    let not_c1: Vec<bool> = cc.iter().map(|c| !c).collect();
    let f = (n / 256) as f32;
    let d: Vec<f32> = edt_periodic(&not_c1, n).iter().map(|v| v / f).collect();
    let mut v: Vec<f32> = (0..n * n).filter(|&k| fc[k]).map(|k| d[k]).collect();
    let tot = v.len().max(1) as f32;
    let s = |t: f32| v.iter().filter(|&&x| x <= t).count() as f32 / tot;
    let (a, b, c) = (s(1.0), s(2.0), s(5.0));
    let dmap: Vec<f32> = (0..n * n).map(|k| if fc[k] { d[k] } else { -1.0 }).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f32| if v.is_empty() { f32::NAN } else { v[((v.len() - 1) as f32 * p).round() as usize] };
    (a, b, c, q(0.5), q(0.9), dmap, d)
}

/// The dendritic sea arms: 8-connected components (periodic) of the ocean lying inside C1's land by ≥ 2 cells of 256², ≥ 1
/// km²; (count, penetration p50, max, km).
fn sea_arms(final_land: &[bool], c1_land: &[bool], d256: &[f32], n: usize) -> (usize, f32, f32) {
    let ck = DOMAIN_KM / n as f32;
    let cell256 = DOMAIN_KM / 256.0;
    let m: Vec<bool> = (0..n * n).map(|k| !final_land[k] && c1_land[k] && d256[k] >= 2.0).collect();
    let mut lab = vec![u32::MAX; n * n];
    let mut pens = Vec::new();
    let mut stack = Vec::new();
    let min_cells = (1.0 / (ck * ck)).ceil() as usize;
    for s in 0..n * n {
        if !m[s] || lab[s] != u32::MAX {
            continue;
        }
        lab[s] = 0;
        stack.push(s);
        let (mut cnt, mut pen) = (0usize, 0f32);
        while let Some(k) = stack.pop() {
            cnt += 1;
            pen = pen.max(d256[k] * cell256);
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if m[j] && lab[j] == u32::MAX {
                        lab[j] = 0;
                        stack.push(j);
                    }
                }
            }
        }
        if cnt >= min_cells {
            pens.push(pen);
        }
    }
    pens.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p50 = if pens.is_empty() { 0.0 } else { pens[pens.len() / 2] };
    (pens.len(), p50, pens.last().copied().unwrap_or(0.0))
}

/// The coast by the rock at the coastline: 256² windows of the 2 048² field, typed by the 256² surface rock (nearest ×8):
/// rocky (basement / shield ≥ 80 %), low (cover ≥ 80 %); per type the median measures and the window count; and the median
/// seaside slope (°) of the land coastline cells.
fn rock_windows(g: &GridF32, rock256: &[u8]) -> ([Option<[f32; 4]>; 2], [usize; 3], f32) {
    let n = 2048;
    let land = land_t(g);
    let cl = coast_cells_t(&land, n);
    let rk = |k: usize| rock256[((k / n) / 8) * 256 + (k % n) / 8];
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
            let typed: Vec<u8> = cells.iter().map(|&k| rk(k)).filter(|&r| r > 0).collect();
            if typed.is_empty() {
                counts[2] += 1;
                continue;
            }
            let hard = typed.iter().filter(|&&r| r >= 2).count() as f32 / typed.len() as f32;
            let ty = if hard >= 0.8 { 0 } else if hard <= 0.2 { 1 } else { 2 };
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
    let sl = slope_deg(&g.data, n, n, DOMAIN_KM / n as f32);
    let mut ss: Vec<f32> = (0..n * n).filter(|&k| cl[k]).map(|k| sl[k]).collect();
    ss.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let sm = if ss.is_empty() { f32::NAN } else { ss[ss.len() / 2] };
    let _ = &cl;
    ([med(&per_type[0]), med(&per_type[1])], counts, sm)
}

/// An altitude / bathymetry image (north up), 2 048² downsampled by `step`: the sea blue by depth (0 to −5 500 m), the land
/// green → brown → white (0 to 3 000 m).
fn alt_png(g: &GridF32, step: usize) -> image::RgbaImage {
    let n = g.width;
    let m = n / step;
    let land = land_t(g);
    let mut img = image::RgbaImage::new(m as u32, m as u32);
    for y in 0..m {
        for x in 0..m {
            let k = (y * step) * n + x * step;
            let z = g.data[k];
            let px = if !land[k] {
                let t = (-z / 5500.0).clamp(0.0, 1.0).sqrt();
                [(150.0 - 130.0 * t) as u8, (200.0 - 160.0 * t) as u8, (240.0 - 120.0 * t) as u8, 255]
            } else {
                let t = (z / 3000.0).clamp(0.0, 1.0);
                if t < 0.4 {
                    let u = t / 0.4;
                    [(70.0 + 110.0 * u) as u8, (140.0 + 20.0 * u) as u8, (60.0 + 30.0 * u) as u8, 255]
                } else {
                    let u = (t - 0.4) / 0.6;
                    [(180.0 + 70.0 * u) as u8, (160.0 + 90.0 * u) as u8, (90.0 + 160.0 * u) as u8, 255]
                }
            };
            img.put_pixel(x as u32, (m - 1 - y) as u32, image::Rgba(px));
        }
    }
    img
}
/// The 64² perimeter map (north up, ×6): continental grey, passive perimeter green, active red, ocean blue.
fn perim_png(map: &[u8]) -> image::RgbaImage {
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..4096 {
        let (x, y) = (k % 64, k / 64);
        let c = match map[k] {
            1 => [190, 190, 180, 255],
            2 => [30, 170, 60, 255],
            3 => [220, 30, 30, 255],
            _ => [40, 70, 140, 255],
        };
        for dy in 0..6 {
            for dx in 0..6 {
                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(c));
            }
        }
    }
    img
}
/// A 25 km crop (128² cells of 195 m, ×2, north up) of a field, around the 128-aligned block with the most coastline cells.
fn coast_crop(g: &GridF32) -> image::RgbaImage {
    let n = g.width;
    let land = land_t(g);
    let cl = coast_cells_t(&land, n);
    let b = n / 128;
    let mut best = (0usize, 0usize, 0usize);
    for by in 0..b {
        for bx in 0..b {
            let c = (0..128 * 128).filter(|&k| cl[(by * 128 + k / 128) * n + bx * 128 + k % 128]).count();
            if c > best.0 {
                best = (c, bx, by);
            }
        }
    }
    let sub = GridF32 { width: 128, height: 128, data: (0..128 * 128).map(|k| g.data[(best.2 * 128 + k / 128) * n + best.1 * 128 + k % 128]).collect() };
    image::imageops::resize(&alt_png_plain(&sub), 256, 256, image::imageops::FilterType::Nearest)
}
/// [`alt_png`] without the torus ocean (the sea = z ≤ 0), at full resolution.
fn alt_png_plain(g: &GridF32) -> image::RgbaImage {
    let n = g.width;
    let mut img = image::RgbaImage::new(n as u32, n as u32);
    for k in 0..n * n {
        let (x, y) = (k % n, k / n);
        let z = g.data[k];
        let px = if z <= 0.0 {
            let t = (-z / 5500.0).clamp(0.0, 1.0).sqrt();
            [(150.0 - 130.0 * t) as u8, (200.0 - 160.0 * t) as u8, (240.0 - 120.0 * t) as u8, 255]
        } else {
            let t = (z / 3000.0).clamp(0.0, 1.0);
            [(70.0 + 180.0 * t) as u8, (140.0 + 110.0 * t) as u8, (60.0 + 190.0 * t) as u8, 255]
        };
        img.put_pixel(x as u32, (n - 1 - y) as u32, image::Rgba(px));
    }
    img
}
/// A transect image: altitude (top, −5 500 … +2 500 m) and s (bottom, 0 … 2), each series a colour.
fn transect_png(z: &[(Vec<f32>, [u8; 3])], s: &[(Vec<f32>, [u8; 3])], path: &std::path::Path) {
    let (w, h) = (900u32, 560u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let mut put = |x: i64, y: i64, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as u32) < w && (y as u32) < h {
            img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    };
    let pz = |v: f32| (20.0 + (2500.0 - v.clamp(-5500.0, 2500.0)) / 8000.0 * 300.0) as i64;
    let ps = |v: f32| (360.0 + (2.0 - v.clamp(0.0, 2.0)) / 2.0 * 180.0) as i64;
    for x in 20..w as i64 - 20 {
        put(x, pz(0.0), [0, 0, 0]);
        put(x, pz(-2000.0), [210, 210, 210]);
        put(x, pz(-4000.0), [210, 210, 210]);
        put(x, ps(1.0), [0, 0, 0]);
        put(x, ps(0.2), [210, 210, 210]);
        put(x, ps(0.915), [160, 160, 230]);
    }
    for (series, f) in [(z, &pz as &dyn Fn(f32) -> i64), (s, &ps as &dyn Fn(f32) -> i64)] {
        for (v, c) in series {
            for (i, &q) in v.iter().enumerate() {
                let x = 20 + (i as f32 / (v.len() - 1).max(1) as f32 * (w as f32 - 40.0)) as i64;
                for d in -1..=1 {
                    put(x, f(q) + d, *c);
                }
            }
        }
    }
    img.save(path).unwrap();
}
/// A transect through a 64² perimeter cell `k` toward its ocean (the 4-direction with an oceanic neighbour), 60 km each side:
/// the 2 048² altitude and C1's final s (nearest).
fn transect(g: &GridF32, s64: &[f64], plate_type: &[PlateType], k: usize) -> (Vec<f32>, Vec<f32>) {
    let (i, j) = (k % 64, k / 64);
    let dirs = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)];
    let dir = dirs
        .iter()
        .copied()
        .find(|&(dx, dy)| plate_type[(((j as i64 + dy).rem_euclid(64)) as usize) * 64 + ((i as i64 + dx).rem_euclid(64)) as usize] == PlateType::Oceanic)
        .unwrap_or((1, 0));
    let n = g.width as i64;
    let (cx, cy) = ((i as i64 * 32 + 16) as f32, (j as i64 * 32 + 16) as f32);
    let ck = DOMAIN_KM / n as f32;
    let steps = (120.0 / ck) as i64;
    let mut zs = Vec::new();
    let mut ss = Vec::new();
    for t in 0..=steps {
        let u = (t as f32 - steps as f32 / 2.0) as f32;
        let (x, y) = ((cx + u * dir.0 as f32).round() as i64, (cy + u * dir.1 as f32).round() as i64);
        let (x, y) = (x.rem_euclid(n), y.rem_euclid(n));
        zs.push(g.data[(y * n + x) as usize]);
        ss.push(s64[((y / 32) * 64 + x / 32) as usize] as f32);
    }
    (zs, ss)
}





// ════════════════════════════════════════ F174 B′/C′ ════════════════════════════════════════

use ymir_core::cascade::margin::seaward_conflicts;

/// The land components (8-connected, periodic) of ≥ 1 km² other than the largest: (count, total km²).
fn islands_off_main(land: &[bool], n: usize) -> (usize, f32) {
    let ck = DOMAIN_KM / n as f32;
    let mut lab = vec![u32::MAX; n * n];
    let mut sizes = Vec::new();
    let mut stack = Vec::new();
    for s in 0..n * n {
        if !land[s] || lab[s] != u32::MAX {
            continue;
        }
        lab[s] = 0;
        stack.push(s);
        let mut c = 0usize;
        while let Some(k) = stack.pop() {
            c += 1;
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if land[j] && lab[j] == u32::MAX {
                        lab[j] = 0;
                        stack.push(j);
                    }
                }
            }
        }
        sizes.push(c);
    }
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    let min = (1.0 / (ck * ck)).ceil() as usize;
    let rest: Vec<usize> = sizes.into_iter().skip(1).filter(|&c| c >= min).collect();
    (rest.len(), rest.iter().sum::<usize>() as f32 * ck * ck)
}

/// The plates at t₀ (64², north up, ×6): each plate a colour; Part A's continental cells bright, the annexed band medium
/// (continental type) or dark (oceanic type), the other oceanic cells pale.
fn plates_png(a: &C1State, b: &C1State) -> image::RgbaImage {
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..4096 {
        let pid = b.plate_id.data()[k] as u32;
        let h = pid.wrapping_mul(2654435761);
        let base = [(60 + (h % 160)) as f32, (60 + ((h >> 8) % 160)) as f32, (60 + ((h >> 16) % 160)) as f32];
        let ca = a.plate_type.data()[k] == PlateType::Continental;
        let cb = b.plate_type.data()[k] == PlateType::Continental;
        let annexed = a.plate_id.data()[k] != b.plate_id.data()[k];
        let f = if ca { 1.0 } else if annexed && cb { 0.7 } else if annexed { 0.45 } else { 0.25 };
        let c = [(base[0] * f + 255.0 * (1.0 - f) * 0.2) as u8, (base[1] * f + 255.0 * (1.0 - f) * 0.2) as u8, (base[2] * f + 255.0 * (1.0 - f) * 0.6) as u8, 255];
        for d in 0..36 {
            img.put_pixel(((k % 64) * 6 + d % 6) as u32, ((63 - k / 64) * 6 + d / 6) as u32, image::Rgba(c));
        }
    }
    img
}

// ════════════════════════════════════════ F175 ════════════════════════════════════════

use ymir_core::cascade::history::{run_c1_sources_mapped_flags, run_c1_sources_mapped_oo};
use ymir_core::tectonics_c1::boundary_classification::clear_oo_seeds;
use ymir_core::tectonics_c1::closures::subduction::OoTracker;
use ymir_core::tectonics_c1::time_loop::{run_with_closures_terms_flags, run_with_closures_terms_oo};

fn out175() -> PathBuf {
    root().join("docs/reports/relief_method/f175_oo_subduction/images")
}

/// C1 under a profile, a margin and its O-O flag (`m.oo`), its Airy sources' record.
fn world_oo(seed: u64, p: C1Profile, m: &Margin) -> World {
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let clk = Clock::start();
    let src = run_c1_sources_mapped_ext(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &SteinSteinParams::default(), &[C1Term::Erosion], Mapping::Airy(Airy::declared()), m.ext(), m.arcs_emerge);
    let cost = clk.stop();
    let snaps = src.hist.normalized(None);
    World { src, st, kin, snaps, cost }
}

/// One initiation (a plate pair becoming subducting), located by a replica of C1's D_init tracker.
#[derive(Clone, Debug)]
struct Initiation {
    step: usize,
    near_continent: bool,
    /// The distance (cells of 64²) to the nearest cell consumed at the previous step (either rule); None when none.
    to_subduction: Option<f32>,
}

#[derive(Clone, Debug, Default)]
struct C1Series {
    active: Vec<f32>,
    map_t0: Vec<u8>,
    map_end: Vec<u8>,
    /// Σ over the run of the stats' O-O fields, and the C-O cells consumed.
    oo_consumed: usize,
    oo_mass: f64,
    oo_arc_cont: f64,
    oo_arc_ocean: f64,
    oo_arc_lost: f64,
    oo_reassigned: usize,
    oo_age: usize,
    oo_age_close: usize,
    oo_tie: usize,
    oo_seeds_cut: usize,
    oo_pairs_end: usize,
    co_consumed: usize,
    /// The DS seeds an O-O boundary would have had (the A7 cut on a copy), Σ over the run.
    oo_seeds_potential: usize,
    initiations: Vec<Initiation>,
    replica_pairs: usize,
    fingers_end: usize,
    /// Σ Δs per cell of the subduction term (arc deposits > 0) and of Davis-Suppe.
    sub_sum: Vec<f64>,
    ds_sum: Vec<f64>,
    andean: usize,
    andean_cells: usize,
}

fn perimeter_map(s: &C1State, kn: &PlateKinematics) -> (f32, Vec<u8>) {
    let n = 64;
    let mut bi = classify_boundaries(&s.plate_id, kn);
    retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
    let cont = |i: usize, j: usize| s.plate_type.get(i, j) == PlateType::Continental;
    let mut map = vec![0u8; n * n];
    let (mut per, mut act) = (0usize, 0usize);
    for j in 0..n {
        for i in 0..n {
            if !cont(i, j) {
                continue;
            }
            map[j * n + i] = 1;
            let nb4 = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
            if !nb4.iter().any(|&(x, y)| !cont(x, y)) {
                continue;
            }
            per += 1;
            let mut a = false;
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    if bi.boundary_type.get((i + dx) % n, (j + dy) % n) == BoundaryType::Convergent {
                        a = true;
                    }
                }
            }
            act += a as usize;
            map[j * n + i] = if a { 3 } else { 2 };
        }
    }
    (act as f32 / per.max(1) as f32, map)
}

/// F175: one C1 run with the margin and its O-O flag, observed per step.
fn c1_series(seed: u64, p: C1Profile, m: &Margin) -> C1Series {
    let n = 64;
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let dt = 0.5 * cfg.dx.min(cfg.dy) / kin.max_velocity().max(1e-12);
    let out = std::cell::RefCell::new(C1Series { sub_sum: vec![0.0; n * n], ds_sum: vec![0.0; n * n], ..Default::default() });
    {
        let (a0, m0) = perimeter_map(&st, &kin);
        let mut o = out.borrow_mut();
        o.active.push(a0);
        o.map_t0 = m0;
    }
    let consumed_prev = std::cell::RefCell::new(vec![false; n * n]);
    let replica = std::cell::RefCell::new(OoTracker::default());
    let mut on_term = |t: C1Term, a: &[f64], b: &[f64]| {
        let mut o = out.borrow_mut();
        if t == C1Term::Subduction {
            let mut c = consumed_prev.borrow_mut();
            for k in 0..n * n {
                o.sub_sum[k] += b[k] - a[k];
                c[k] = b[k] < a[k];
            }
        }
        if t == C1Term::DavisSuppe {
            for k in 0..n * n {
                o.ds_sum[k] += b[k] - a[k];
            }
        }
    };
    run_with_closures_terms_flags(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        m.oo,
        m.merge_at_collision,
        |step, s, kn| {
            let mut o = out.borrow_mut();
            let x = s.last_step_stats.subduction;
            o.oo_consumed += x.oo_cells_consumed;
            o.oo_mass += x.oo_mass_consumed;
            o.oo_arc_cont += x.oo_arc_continental;
            o.oo_arc_ocean += x.oo_arc_oceanic;
            o.oo_arc_lost += x.oo_arc_lost;
            o.oo_reassigned += x.oo_reassigned;
            o.oo_age += x.oo_age_decided;
            o.oo_age_close += x.oo_age_close;
            o.oo_tie += x.oo_tie_decided;
            o.oo_seeds_cut += x.oo_seeds_cleared;
            o.oo_pairs_end = x.oo_active_pairs;
            o.co_consumed += x.cells_consumed;
            let (a, mp) = perimeter_map(s, kn);
            o.active.push(a);
            o.map_end = mp;
            // the DS seeds an O-O boundary would have (a copy, cut)
            let mut bi = classify_boundaries(&s.plate_id, kn);
            retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
            let mut bi2 = classify_boundaries(&s.plate_id, kn);
            retarget_upper_plate_continental(&mut bi2, &s.plate_id, &s.plate_type);
            o.oo_seeds_potential += clear_oo_seeds(&mut bi2, &s.plate_id, &s.plate_type);
            // the replica of the D_init tracker (the loop's next step reads this state)
            if let Some(oo) = m.oo {
                let new = replica.borrow_mut().update(&s.plate_id, &s.plate_type, &bi, kn, dt, oo.d_init_domain(), step + 1);
                let cont: Vec<bool> = s.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
                let dcont = edt_periodic(&cont, 64);
                let prev = consumed_prev.borrow();
                for (pa, pb) in new {
                    // the pair's O-O contacts
                    let mut near = false;
                    let mut best: Option<f32> = None;
                    for k in 0..n * n {
                        let (i, j) = (k % n, k / n);
                        let pid = s.plate_id.get(i, j);
                        if (pid != pa && pid != pb) || bi.boundary_type.get(i, j) != BoundaryType::Convergent || s.plate_type.get(i, j) != PlateType::Oceanic {
                            continue;
                        }
                        let other = if pid == pa { pb } else { pa };
                        let nb4 = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
                        if !nb4.iter().any(|&(x, y)| s.plate_id.get(x, y) == other) {
                            continue;
                        }
                        if dcont[k] <= 2.0 {
                            near = true;
                        }
                        for q in 0..n * n {
                            if prev[q] {
                                let mut dx = ((q % n) as i64 - i as i64).rem_euclid(64);
                                let mut dy = ((q / n) as i64 - j as i64).rem_euclid(64);
                                dx = dx.min(64 - dx);
                                dy = dy.min(64 - dy);
                                let d = ((dx * dx + dy * dy) as f32).sqrt();
                                best = Some(best.map_or(d, |b: f32| b.min(d)));
                            }
                        }
                    }
                    o.initiations.push(Initiation { step: step + 1, near_continent: near, to_subduction: best });
                }
                o.replica_pairs = replica.borrow().active.len();
            }
        },
        &mut on_term,
    );
    let mut o = out.into_inner();
    // V5: one-cell fingers (an oceanic cell with ≥ 3 of its 4 neighbours on another plate)
    o.fingers_end = (0..n * n)
        .filter(|&k| {
            let (i, j) = (k % n, k / n);
            if st.plate_type.get(i, j) != PlateType::Oceanic {
                return false;
            }
            let pid = st.plate_id.get(i, j);
            [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)].iter().filter(|&&(x, y)| st.plate_id.get(x, y) != pid).count() >= 3
        })
        .count();
    // the Andean margins: 8-connected groups of ≥ 3 perimeter cells passive at t₀ and active at the end
    let newly: Vec<bool> = (0..n * n).map(|k| o.map_end[k] == 3 && o.map_t0[k] == 2).collect();
    let mut seen = vec![false; n * n];
    for s0 in 0..n * n {
        if !newly[s0] || seen[s0] {
            continue;
        }
        let mut stack = vec![s0];
        seen[s0] = true;
        let mut c = 0;
        while let Some(k) = stack.pop() {
            c += 1;
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if newly[j] && !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        if c >= 3 {
            o.andean += 1;
            o.andean_cells += c;
        }
    }
    o
}

/// W₁₀ by segment type (256²): the land within 10 km of the coast, typed by its nearest perimeter cell of C1's end map
/// (active / passive); the share with a slope ≥ 8°. (active, passive, cells active, cells passive).
fn w10_by_type(zm: &GridF32, map_end: &[u8]) -> (f32, f32, usize, usize) {
    let land = land_t(zm);
    let d = edt_periodic(&land, 256);
    let sl = slope_deg(&zm.data, 256, 256, CELL);
    let up = |v: u8| -> Vec<bool> { (0..65536).map(|k| map_end[((k / 256) / 4) * 64 + (k % 256) / 4] == v).collect() };
    let (act, pas) = (up(3), up(2));
    let da = edt_periodic(&act.iter().map(|b| !b).collect::<Vec<_>>(), 256);
    let dp = edt_periodic(&pas.iter().map(|b| !b).collect::<Vec<_>>(), 256);
    let (mut na, mut wa, mut np, mut wp) = (0usize, 0usize, 0usize, 0usize);
    for k in 0..65536 {
        if !land[k] || d[k] * CELL > 10.0 {
            continue;
        }
        let steep = sl[k] >= 8.0;
        if da[k] < dp[k] {
            na += 1;
            wa += steep as usize;
        } else {
            np += 1;
            wp += steep as usize;
        }
    }
    (wa as f32 / na.max(1) as f32, wp as f32 / np.max(1) as f32, na, np)
}

/// The margins at the end (64², north up, ×6): continental grey, passive green, active red, arc deposits (Σ subduction term
/// > 0.02) orange on the ocean, ocean blue.
fn margins_png(map_end: &[u8], sub_sum: &[f64]) -> image::RgbaImage {
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..4096 {
        let c = match map_end[k] {
            1 => [190, 190, 180, 255],
            2 => [30, 170, 60, 255],
            3 => [220, 30, 30, 255],
            _ if sub_sum[k] > 0.02 => [240, 150, 30, 255],
            _ => [40, 70, 140, 255],
        };
        for d in 0..36 {
            img.put_pixel(((k % 64) * 6 + d % 6) as u32, ((63 - k / 64) * 6 + d / 6) as u32, image::Rgba(c));
        }
    }
    img
}

/// The share of the islands' area (≥ 1 km², off the main land, 2 048²) lying on C1's arc deposits (Σ subduction term > 0, 64²).
fn islands_on_arcs(land: &[bool], sub_sum: &[f64]) -> f32 {
    let n = 2048;
    let mut lab = vec![u32::MAX; n * n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for s0 in 0..n * n {
        if !land[s0] || lab[s0] != u32::MAX {
            continue;
        }
        let id = comps.len() as u32;
        lab[s0] = id;
        let mut stack = vec![s0];
        let mut cells = Vec::new();
        while let Some(k) = stack.pop() {
            cells.push(k);
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if land[j] && lab[j] == u32::MAX {
                        lab[j] = id;
                        stack.push(j);
                    }
                }
            }
        }
        comps.push(cells);
    }
    comps.sort_by_key(|c| std::cmp::Reverse(c.len()));
    let min = (1.0 / (DOMAIN_KM / n as f32).powi(2)).ceil() as usize;
    let isl: Vec<&Vec<usize>> = comps.iter().skip(1).filter(|c| c.len() >= min).collect();
    let tot: usize = isl.iter().map(|c| c.len()).sum();
    let on: usize = isl.iter().map(|c| c.iter().filter(|&&k| sub_sum[((k / n) / 32) * 64 + (k % n) / 32] > 0.0).count()).sum();
    on as f32 / tot.max(1) as f32
}

// ── F175, the bench by tiers (the author's budget rule): each world's lines are written to `results/` as soon as ready, and a
// rerun reads them back instead of recomputing ──

/// F174 B′'s selected plates (`f174b2_bench_output.txt`): the fixed continent of every tier.
fn fixed_profile(seed: u64) -> C1Profile {
    let (k, start) = match seed {
        8 => (8, 2),
        18 => (8, 11),
        19 => (7, 1),
        23 => (8, 6),
        29 => (9, 11),
        _ => unreachable!("the panel"),
    };
    C1Profile::Selected { plates: 16, k, start, compact: true }
}
fn variants175() -> [(&'static str, Margin); 3] {
    [("Part A", Margin::off()), ("B′", Margin::seaward()), ("B′ + option 5", Margin::seaward_oo())]
}
fn results_dir() -> PathBuf {
    root().join("docs/reports/relief_method/f175_oo_subduction/results")
}
/// The resumable step: `name`'s lines from `results/<name>.txt`, else computed by `f`, written, and returned.
fn cached(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results_dir().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results_dir()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
/// The land components (8-connected, periodic) of a 64² or 256² mask other than the largest with ≥ `min_cells`: (count, cells).
fn islands_n(land: &[bool], n: usize, min_cells: usize) -> (usize, usize) {
    let mut lab = vec![false; n * n];
    let mut sizes = Vec::new();
    for s0 in 0..n * n {
        if !land[s0] || lab[s0] {
            continue;
        }
        lab[s0] = true;
        let mut stack = vec![s0];
        let mut c = 0;
        while let Some(k) = stack.pop() {
            c += 1;
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if land[j] && !lab[j] {
                        lab[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        sizes.push(c);
    }
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    let rest: Vec<usize> = sizes.into_iter().skip(1).filter(|&c| c >= min_cells).collect();
    (rest.len(), rest.iter().sum())
}

// ════════════════════════════════════════ F176 ════════════════════════════════════════

fn out176() -> PathBuf {
    root().join("docs/reports/relief_method/f176_accretion/images")
}
fn results176() -> PathBuf {
    root().join("docs/reports/relief_method/f176_accretion/results")
}
/// The resumable step (F176's results directory).
fn cached176(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results176().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results176()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
const SCREEN: [u64; 10] = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11];
fn variants176() -> [(&'static str, Margin); 3] {
    [("F175", Margin::seaward_oo()), ("F176", Margin::seaward_oo_collision()), ("F176 + B", Margin { arcs_emerge: true, ..Margin::seaward_oo_collision() })]
}
fn key176(tier: &str, seed: u64, name: &str) -> String {
    format!("{tier}_seed{seed}_{}", name.replace(' ', "_").replace('+', "p"))
}

/// The F176 per-step reading of one C1 run: the active share (3 × 3), the plates, and at the end the « within 3 cells » share,
/// the Andean groups per continent; plus F175's series (V1–V6) through `c1_series`.
#[derive(Clone, Debug, Default)]
struct Margins176 {
    plates: Vec<usize>,
    near3_end: f32,
    continents: usize,
    continents_andean: usize,
    andean_groups: usize,
    andean_km: f32,
    per_continent: Vec<(usize, usize, f32)>,
}
fn margins176(seed: u64, p: C1Profile, m: &Margin, ser: &C1Series) -> Margins176 {
    let n = 64;
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let plates = std::cell::RefCell::new(Vec::new());
    let last_consumed = std::cell::RefCell::new(vec![usize::MAX; n * n]);
    let step_now = std::cell::Cell::new(0usize);
    let mut on_term = |t: C1Term, a: &[f64], b: &[f64]| {
        if t == C1Term::Subduction {
            let mut lc = last_consumed.borrow_mut();
            for k in 0..n * n {
                if b[k] < a[k] {
                    lc[k] = step_now.get();
                }
            }
        }
    };
    run_with_closures_terms_flags(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        m.oo,
        m.merge_at_collision,
        |step, s, _| {
            step_now.set(step + 1);
            plates.borrow_mut().push(s.plate_id.data().iter().copied().collect::<std::collections::BTreeSet<u16>>().len());
        },
        &mut on_term,
    );
    let lc = last_consumed.into_inner();
    let recent: Vec<bool> = lc.iter().map(|&s| s != usize::MAX && s + 10 >= cfg.n_steps).collect();
    let d_recent = edt_periodic(&recent.iter().map(|r| !r).collect::<Vec<_>>(), 64);
    let per: Vec<usize> = (0..n * n).filter(|&k| ser.map_end[k] >= 2).collect();
    let near3 = per.iter().filter(|&&k| d_recent[k] <= 3.0).count() as f32 / per.len().max(1) as f32;
    // the continents: 4-connected components of continental cells (map values 1–3) with ≥ 20 cells
    let cont: Vec<bool> = ser.map_end.iter().map(|&v| v >= 1).collect();
    let mut comp = vec![usize::MAX; n * n];
    let mut sizes = Vec::new();
    for s0 in 0..n * n {
        if !cont[s0] || comp[s0] != usize::MAX {
            continue;
        }
        let id = sizes.len();
        comp[s0] = id;
        let mut stack = vec![s0];
        let mut c = 0;
        while let Some(k) = stack.pop() {
            c += 1;
            let (x, y) = (k % n, k / n);
            for j in [y * n + (x + 1) % n, y * n + (x + n - 1) % n, ((y + 1) % n) * n + x, ((y + n - 1) % n) * n + x] {
                if cont[j] && comp[j] == usize::MAX {
                    comp[j] = id;
                    stack.push(j);
                }
            }
        }
        sizes.push(c);
    }
    let big: Vec<usize> = (0..sizes.len()).filter(|&i| sizes[i] >= 20).collect();
    // the Andean groups (8-connected, ≥ 3 cells passive at t0 and active at the end), assigned to their continent
    let newly: Vec<bool> = (0..n * n).map(|k| ser.map_end[k] == 3 && ser.map_t0[k] == 2).collect();
    let mut seen = vec![false; n * n];
    let mut per_cont: std::collections::BTreeMap<usize, (usize, usize)> = std::collections::BTreeMap::new();
    for s0 in 0..n * n {
        if !newly[s0] || seen[s0] {
            continue;
        }
        let mut stack = vec![s0];
        seen[s0] = true;
        let mut cells = 0;
        while let Some(k) = stack.pop() {
            cells += 1;
            let (x, y) = (k % n, k / n);
            for dy in [n - 1, 0, 1] {
                for dx in [n - 1, 0, 1] {
                    let j = ((y + dy) % n) * n + (x + dx) % n;
                    if newly[j] && !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        if cells >= 3 && comp[s0] != usize::MAX {
            let e = per_cont.entry(comp[s0]).or_insert((0, 0));
            e.0 += 1;
            e.1 += cells;
        }
    }
    let per_continent: Vec<(usize, usize, f32)> = big.iter().map(|&c| {
        let (g, cells) = per_cont.get(&c).copied().unwrap_or((0, 0));
        (sizes[c], g, cells as f32 * 6.25)
    }).collect();
    Margins176 {
        plates: plates.into_inner(),
        near3_end: near3,
        continents: big.len(),
        continents_andean: per_continent.iter().filter(|x| x.1 > 0).count(),
        andean_groups: per_continent.iter().map(|x| x.1).sum(),
        andean_km: per_continent.iter().map(|x| x.2).sum(),
        per_continent,
    }
}

/// One C1 line of F176 (world_oo + c1_series + margins176), the margins image saved.
fn c1_line176(seed: u64, p: C1Profile, name: &str, m: &Margin, img_key: &str) -> String {
    let clk = Clock::start();
    let w = world_oo(seed, p, m);
    let ser = c1_series(seed, p, m);
    let mg = margins176(seed, p, m, &ser);
    let cost = clk.stop();
    let last = w.snaps.last().unwrap();
    let land64: Vec<bool> = last.data.iter().map(|&v| v > 0.5).collect();
    let circ = circumnavigation(&land64, 64);
    let (isl, isl_cells) = islands_n(&land64, 64, 1);
    let n_init = ser.initiations.len();
    let near = ser.initiations.iter().filter(|x| x.near_continent).count();
    let within = ser.initiations.iter().filter(|x| x.to_subduction.is_some_and(|d| d <= 32.0)).count();
    let cons_err = (ser.oo_arc_cont + ser.oo_arc_ocean + ser.oo_arc_lost - 0.5 * ser.oo_mass).abs();
    let every: Vec<String> = (0..ser.active.len()).step_by(50).map(|i| format!("{i}:{:.0}/{}", 100.0 * ser.active[i], if i == 0 { 0 } else { mg.plates[i - 1] })).collect();
    margins_png(&ser.map_end, &ser.sub_sum).save(out176().join(format!("margins_end_{img_key}.png"))).unwrap();
    format!(
        "   seed {seed} {name} ({p:?}): C1 end land {:.1} %, circumnavigable {} r* {:?}, islands (64², off the main land) {isl} ({isl_cells} cells) · active at t0 {:.1} % → end {:.1} % (within 3 cells of a recent consumption {:.1} %) · plates 16 → {} · continents (≥ 20 cells) {} of which with ≥ 1 Andean margin {} · Andean groups {} ({:.0} km) · per continent (cells, groups, km) {:?} · C-O cells consumed {} · O-O consumed {} (age-decided {}, |Δage| < 0.01 {}, tie {}), reassigned {}, pairs {} · arc mass continent / ocean / lost {:.3} / {:.3} / {:.3} (error {:.1e}) · DS O-O seeds cut {} · fingers {} · V1 {near} / {n_init} near a continent · V2 {within} / {n_init} within 200 km · C1 CPU {:.1} s\n      the active share / plates (step: % / n) {}\n",
        100.0 * land64.iter().filter(|&&b| b).count() as f32 / 4096.0, !circ.wraps, circ.sea_r, 100.0 * ser.active[0], 100.0 * ser.active.last().unwrap(), 100.0 * mg.near3_end,
        mg.plates.last().copied().unwrap_or(0), mg.continents, mg.continents_andean, mg.andean_groups, mg.andean_km, mg.per_continent,
        ser.co_consumed, ser.oo_consumed, ser.oo_age, ser.oo_age_close, ser.oo_tie, ser.oo_reassigned, ser.oo_pairs_end,
        ser.oo_arc_cont, ser.oo_arc_ocean, ser.oo_arc_lost, cons_err, ser.oo_seeds_cut, ser.fingers_end, cost.1, every.join(" ")
    )
}


// ════════════════════════════════════════ F177 ════════════════════════════════════════

use ymir_core::tectonics_c1::boundary_classification::{BoundaryInfo, oc_override_seed_mask};
use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;

fn out177() -> PathBuf {
    root().join("docs/reports/relief_method/f177_accretion_close/images")
}
fn results177() -> PathBuf {
    root().join("docs/reports/relief_method/f177_accretion_close/results")
}
/// The resumable step (F177's results directory).
fn cached177(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results177().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results177()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
/// The numbers of a `#KEY v v v` line of a results text.
fn data_line(s: &str, key: &str) -> Option<Vec<f64>> {
    s.lines().find(|l| l.starts_with(&format!("#{key} "))).map(|l| l.split_whitespace().skip(1).filter_map(|x| x.parse().ok()).collect())
}

const N64: usize = 64;
const CELL64_KM: f32 = 6.25;
fn nb4_64(k: usize) -> [usize; 4] {
    let (x, y) = (k % N64, k / N64);
    [y * N64 + (x + 1) % N64, y * N64 + (x + N64 - 1) % N64, ((y + 1) % N64) * N64 + x, ((y + N64 - 1) % N64) * N64 + x]
}
fn nb8_n(k: usize, n: usize) -> [usize; 8] {
    let (x, y) = (k % n, k / n);
    let mut o = [0usize; 8];
    let mut i = 0;
    for dy in [n - 1, 0, 1] {
        for dx in [n - 1, 0, 1] {
            if dx == 0 && dy == 0 {
                continue;
            }
            o[i] = ((y + dy) % n) * n + (x + dx) % n;
            i += 1;
        }
    }
    o
}
/// 8-connected components of a mask (periodic), each a list of cells.
fn comps8(mask: &[bool], n: usize) -> Vec<Vec<usize>> {
    let mut seen = vec![false; n * n];
    let mut out = Vec::new();
    for s0 in 0..n * n {
        if !mask[s0] || seen[s0] {
            continue;
        }
        seen[s0] = true;
        let mut stack = vec![s0];
        let mut cells = Vec::new();
        while let Some(k) = stack.pop() {
            cells.push(k);
            for j in nb8_n(k, n) {
                if mask[j] && !seen[j] {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        out.push(cells);
    }
    out
}

/// The state a C1 step starts from, read as the loop reads it (the boundary, the velocity flag before A7, Davis-Suppe's typed
/// wedge, the continent-ocean contact cells).
struct StepInfo {
    bi: BoundaryInfo,
    upper_vel: Vec<bool>,
    oc_wedge: Vec<bool>,
    co_cell: Vec<bool>,
    plate_id: Vec<u16>,
    plate_type: Vec<PlateType>,
    age: Vec<f64>,
}
fn step_info(s: &C1State, kn: &PlateKinematics, ds_max: f64, oo: bool) -> StepInfo {
    step_info_a(s, kn, ds_max, oo, None)
}
/// [`step_info`] with F178's Andean seeds (the reach, and the D_init tracker as the step's Davis-Suppe stage sees it).
fn step_info_a(s: &C1State, kn: &PlateKinematics, ds_max: f64, oo: bool, andean: Option<(usize, &OoTracker)>) -> StepInfo {
    let mut bi = classify_boundaries(&s.plate_id, kn);
    retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
    let upper_vel = bi.upper_plate_mask.data().to_vec();
    let uv = bi.upper_plate_mask.clone();
    if oo {
        clear_oo_seeds(&mut bi, &s.plate_id, &s.plate_type);
    }
    let seeds = andean.map(|(gap, tr)| andean_seeds(&s.plate_id, &s.plate_type, &s.age, &bi, &uv, kn, tr, gap));
    if let Some(sd) = &seeds {
        for (k, &on) in sd.data().iter().enumerate() {
            if on {
                bi.upper_plate_mask.set(k % N64, k / N64, true);
            }
        }
    }
    let mut oc_seed = oc_override_seed_mask(&bi, &s.plate_id, &s.plate_type);
    if let Some(sd) = &seeds {
        for (k, &on) in sd.data().iter().enumerate() {
            if on {
                oc_seed.set(k % N64, k / N64, true);
            }
        }
    }
    let (_, ocw) = wedge_distance_intra_plate_typed(&s.plate_id, &bi.upper_plate_mask, &oc_seed, ds_max);
    let plate_id: Vec<u16> = s.plate_id.data().to_vec();
    let plate_type: Vec<PlateType> = s.plate_type.data().to_vec();
    let co_cell = (0..N64 * N64)
        .map(|k| plate_type[k] == PlateType::Oceanic && nb4_64(k).iter().any(|&j| plate_id[j] != plate_id[k] && plate_type[j] == PlateType::Continental))
        .collect();
    StepInfo { bi, upper_vel, oc_wedge: ocw.data().to_vec(), co_cell, plate_id, plate_type, age: s.age.data().to_vec() }
}

/// F177's observer of one C1 run: the per-origin Σ Δs > 0 (Andean, collision, O-O arcs), Σ Δs of every term but the erosion,
/// the rift thinning, the merges' oceanic cells, the rift strips, the replayed D_init tracker, the band at t₀.
struct Obs {
    a_sum: Vec<f64>,
    c_sum: Vec<f64>,
    o_sum: Vec<f64>,
    src_sum: Vec<f64>,
    rift_sum: Vec<f64>,
    merged_ocean: Vec<bool>,
    strip: Vec<bool>,
    merges: usize,
    splits: usize,
    s0: Vec<f64>,
    band: Vec<bool>,
    tracker: OoTracker,
    st: C1State,
    kin: PlateKinematics,
    st0: C1State,
    kin0: PlateKinematics,
    ds_max: f64,
    /// F178: Davis-Suppe on O-C wedge cells only, and every positive subduction increment (the arcs).
    w_sum: Vec<f64>,
    arc_sum: Vec<f64>,
}
fn observe(seed: u64, p: C1Profile, m: &Margin) -> Obs {
    let nn = N64 * N64;
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let st0 = p.init_state_margin(64, seed, m);
    let st_off = p.init_state_margin(64, seed, &Margin::off());
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let kin0 = PlateKinematics::preset_phase_1_1(st.num_plates);
    let dt = 0.5 * cfg.dx.min(cfg.dy) / kin.max_velocity().max(1e-12);
    let ds_max = C1Closures::default().davis_suppe.scaled_to_grid(64).max_distance;
    let oo = m.oo.is_some();
    let band: Vec<bool> = (0..nn).map(|k| st.plate_id.data()[k] != st_off.plate_id.data()[k] || st.plate_type.data()[k] != st_off.plate_type.data()[k]).collect();
    let mut tracker = OoTracker::default();
    let info0 = step_info_a(&st, &kin, ds_max, oo, m.andean_gap.map(|g| (g, &tracker)));
    if let Some(o) = m.oo {
        tracker.update(&st.plate_id, &st.plate_type, &info0.bi, &kin, dt, o.d_init_domain(), 0);
    }
    let cur = std::cell::RefCell::new(info0);
    let tr = std::cell::RefCell::new(tracker);
    let acc = std::cell::RefCell::new((
        vec![0f64; nn],
        vec![0f64; nn],
        vec![0f64; nn],
        vec![0f64; nn],
        vec![0f64; nn],
        vec![false; nn],
        vec![false; nn],
        0usize,
        0usize,
        vec![0f64; nn],
        vec![0f64; nn],
    ));
    let s0 = st.s.data().to_vec();
    let mut on_term = |t: C1Term, a: &[f64], b: &[f64]| {
        let mut o = acc.borrow_mut();
        if t != C1Term::Erosion {
            for k in 0..nn {
                o.3[k] += b[k] - a[k];
            }
        }
        match t {
            C1Term::DavisSuppe => {
                let c = cur.borrow();
                for k in 0..nn {
                    let d = b[k] - a[k];
                    if d > 0.0 {
                        if c.oc_wedge[k] {
                            o.0[k] += d;
                            o.9[k] += d;
                        } else {
                            o.1[k] += d;
                        }
                    }
                }
            }
            C1Term::Subduction => {
                let c = cur.borrow();
                let consumed: Vec<usize> = (0..nn).filter(|&k| b[k] < a[k]).collect();
                for k in 0..nn {
                    let d = b[k] - a[k];
                    if d > 0.0 {
                        o.10[k] += d;
                    }
                    if d <= 0.0 || consumed.is_empty() {
                        continue;
                    }
                    let (x, y) = ((k % N64) as i64, (k / N64) as i64);
                    let mut best = (i64::MAX, false);
                    for &q in &consumed {
                        let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
                        dx = dx.min(64 - dx);
                        let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
                        dy = dy.min(64 - dy);
                        let d2 = dx * dx + dy * dy;
                        if d2 < best.0 || (d2 == best.0 && c.co_cell[q] && !best.1) {
                            best = (d2, c.co_cell[q]);
                        }
                    }
                    if best.1 {
                        o.0[k] += d;
                    } else {
                        o.2[k] += d;
                    }
                }
            }
            C1Term::RiftThinning => {
                for k in 0..nn {
                    o.4[k] += b[k] - a[k];
                }
            }
            _ => {}
        }
    };
    run_with_closures_terms_ext(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        m.ext(),
        |step, s, kn| {
            {
                let c = cur.borrow();
                let mut o = acc.borrow_mut();
                let now = s.plate_id.data();
                if s.last_step_stats.accretion.merges_count > 0 {
                    let present: std::collections::BTreeSet<u16> = now.iter().copied().collect();
                    let before: std::collections::BTreeSet<u16> = c.plate_id.iter().copied().collect();
                    for &gone in before.difference(&present) {
                        let mut cnt: std::collections::BTreeMap<u16, usize> = std::collections::BTreeMap::new();
                        for k in 0..nn {
                            if c.plate_id[k] == gone {
                                *cnt.entry(now[k]).or_insert(0) += 1;
                            }
                        }
                        if let Some((&w, _)) = cnt.iter().max_by_key(|(_, c)| **c) {
                            o.7 += 1;
                            for k in 0..nn {
                                if now[k] == w && s.plate_type.data()[k] == PlateType::Oceanic {
                                    o.5[k] = true;
                                }
                            }
                        }
                    }
                }
                for &np in &s.last_step_stats.rifting_split.new_plate_ids_created {
                    o.8 += 1;
                    for k in 0..nn {
                        if now[k] == np {
                            o.6[k] = true;
                        }
                    }
                }
            }
            let trb = tr.borrow();
            let ni = step_info_a(s, kn, ds_max, oo, m.andean_gap.map(|g| (g, &*trb)));
            drop(trb);
            if let Some(o) = m.oo {
                if step + 1 < cfg.n_steps {
                    tr.borrow_mut().update(&s.plate_id, &s.plate_type, &ni.bi, kn, dt, o.d_init_domain(), step + 1);
                }
            }
            *cur.borrow_mut() = ni;
        },
        &mut on_term,
    );
    let (a_sum, c_sum, o_sum, src_sum, rift_sum, merged_ocean, strip, merges, splits, w_sum, arc_sum) = acc.into_inner();
    Obs { a_sum, c_sum, o_sum, src_sum, rift_sum, merged_ocean, strip, merges, splits, s0, band, tracker: tr.into_inner(), st, kin, st0, kin0, ds_max, w_sum, arc_sum }
}

/// A1 at one state: per cell 0 not a coast, 1 active (walk), 2 passive (walk), 3 active (3 × 3, on a boundary), 4 passive (3 × 3),
/// 5 enclosed; the walk's depth (cells); the trench cells the active walks reached; the mixed ties.
struct A1 {
    cls: Vec<u8>,
    dist: Vec<u32>,
    trench: Vec<bool>,
    mixed: usize,
    /// F178: the majority rule at the walk's first depth (strictly more active than passive boundary cells; a tie → passive).
    major: Vec<u8>,
}
fn a1(info: &StepInfo, kn: &PlateKinematics, tracker: &OoTracker) -> A1 {
    let nn = N64 * N64;
    let pid = &info.plate_id;
    let pt = &info.plate_type;
    let bt = info.bi.boundary_type.data();
    let is_bnd = |k: usize| nb4_64(k).iter().any(|&j| pid[j] != pid[k]);
    let mut memo: Vec<Option<bool>> = vec![None; nn];
    let mut classify_b = |b: usize| -> bool {
        if let Some(v) = memo[b] {
            return v;
        }
        let v = (|| {
            if bt[b] != BoundaryType::Convergent {
                return false;
            }
            if nb4_64(b).iter().any(|&j| pid[j] != pid[b] && pt[j] == PlateType::Continental) {
                return false;
            }
            let (vxc, vyc) = kn.velocities[pid[b] as usize];
            let (x, y) = (b % N64, b / N64);
            let dirs = [((x + 1) % N64, y, 1.0, 0.0), ((x + N64 - 1) % N64, y, -1.0, 0.0), (x, (y + 1) % N64, 0.0, 1.0), (x, (y + N64 - 1) % N64, 0.0, -1.0)];
            let mut best: Option<(usize, f64)> = None;
            for (ni, nj, nxn, nyn) in dirs {
                let q = nj * N64 + ni;
                if pid[q] == pid[b] {
                    continue;
                }
                let (vxn, vyn) = kn.velocities[pid[q] as usize];
                let dot = (vxc - vxn) * nxn + (vyc - vyn) * nyn;
                if dot > 0.0 && tracker.is_active(pid[b], pid[q]) && best.is_none_or(|(_, d)| dot >= d) {
                    best = Some((q, dot));
                }
            }
            let Some((q, _)) = best else { return false };
            let (ac, an) = (info.age[b], info.age[q]);
            let sinks = if ac != an {
                ac > an
            } else {
                match (info.upper_vel[b], info.upper_vel[q]) {
                    (false, true) => true,
                    (true, false) => false,
                    _ => pid[b] < pid[q],
                }
            };
            !sinks
        })();
        memo[b] = Some(v);
        v
    };
    let mut cls = vec![0u8; nn];
    let mut major = vec![0u8; nn];
    let mut dist = vec![0u32; nn];
    let mut trench = vec![false; nn];
    let mut mixed = 0;
    for k in 0..nn {
        if pt[k] != PlateType::Continental || !nb4_64(k).iter().any(|&j| pt[j] == PlateType::Oceanic) {
            continue;
        }
        if is_bnd(k) {
            let (x, y) = (k % N64, k / N64);
            let mut conv = false;
            for dy in [N64 - 1, 0, 1] {
                for dx in [N64 - 1, 0, 1] {
                    if bt[((y + dy) % N64) * N64 + (x + dx) % N64] == BoundaryType::Convergent {
                        conv = true;
                    }
                }
            }
            cls[k] = if conv { 3 } else { 4 };
            major[k] = cls[k];
            continue;
        }
        let p = pid[k];
        let mut seen = vec![false; nn];
        seen[k] = true;
        let mut frontier: Vec<usize> = nb4_64(k).iter().copied().filter(|&j| pt[j] == PlateType::Oceanic && pid[j] == p).collect();
        for &j in &frontier {
            seen[j] = true;
        }
        let mut depth = 1u32;
        loop {
            let bnds: Vec<usize> = frontier.iter().copied().filter(|&q| is_bnd(q)).collect();
            if !bnds.is_empty() {
                let acts: Vec<bool> = bnds.iter().map(|&b| classify_b(b)).collect();
                let any = acts.iter().any(|&a| a);
                if any && acts.iter().any(|&a| !a) {
                    mixed += 1;
                }
                cls[k] = if any { 1 } else { 2 };
                major[k] = if acts.iter().filter(|&&a| a).count() * 2 > acts.len() { 1 } else { 2 };
                dist[k] = depth;
                for (b, a) in bnds.iter().zip(&acts) {
                    if *a {
                        trench[*b] = true;
                    }
                }
                break;
            }
            let mut next = Vec::new();
            for &q in &frontier {
                for j in nb4_64(q) {
                    if !seen[j] && pt[j] == PlateType::Oceanic && pid[j] == p {
                        seen[j] = true;
                        next.push(j);
                    }
                }
            }
            if next.is_empty() {
                cls[k] = 5;
                major[k] = 5;
                break;
            }
            frontier = next;
            depth += 1;
        }
    }
    A1 { cls, dist, trench, mixed, major }
}
fn is_active_cls(c: u8) -> bool {
    c == 1 || c == 3
}
fn is_passive_cls(c: u8) -> bool {
    c == 2 || c == 4
}

/// The continents (4-connected components of continental cells with ≥ 20 cells): the component id per cell (usize::MAX elsewhere)
/// and the qualifying ids.
fn continents64(pt: &[PlateType]) -> (Vec<usize>, Vec<usize>) {
    let nn = N64 * N64;
    let mut comp = vec![usize::MAX; nn];
    let mut sizes = Vec::new();
    for s0 in 0..nn {
        if pt[s0] != PlateType::Continental || comp[s0] != usize::MAX {
            continue;
        }
        let id = sizes.len();
        comp[s0] = id;
        let mut stack = vec![s0];
        let mut c = 0;
        while let Some(k) = stack.pop() {
            c += 1;
            for j in nb4_64(k) {
                if pt[j] == PlateType::Continental && comp[j] == usize::MAX {
                    comp[j] = id;
                    stack.push(j);
                }
            }
        }
        sizes.push(c);
    }
    let big = (0..sizes.len()).filter(|&i| sizes[i] >= 20).collect();
    (comp, big)
}
/// Per qualifying continent, the groups of ≥ 3 8-connected cells of `mask`: (continent id, groups, cells).
fn groups_per_continent(mask: &[bool], comp: &[usize], big: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut per: std::collections::BTreeMap<usize, (usize, usize)> = big.iter().map(|&c| (c, (0, 0))).collect();
    for g in comps8(mask, N64) {
        if g.len() >= 3 {
            if let Some(e) = per.get_mut(&comp[g[0]]) {
                e.0 += 1;
                e.1 += g.len();
            }
        }
    }
    per.into_iter().map(|(c, (g, n))| (c, g, n)).collect()
}

/// B4: the sea arms ≤ 2 cells wide (the sea cells a 3 × 3 opening removes), the one-cell-wide cells, the channels (8-connected).
struct Channels {
    narrow: Vec<bool>,
    one: Vec<bool>,
    chans: Vec<Vec<usize>>,
    ones: Vec<Vec<usize>>,
}
fn channels(land: &[bool]) -> Channels {
    let nn = N64 * N64;
    let sea: Vec<bool> = land.iter().map(|l| !l).collect();
    let kept: Vec<bool> = (0..nn).map(|k| sea[k] && nb8_n(k, N64).iter().all(|&j| sea[j])).collect();
    let opened: Vec<bool> = (0..nn).map(|k| kept[k] || nb8_n(k, N64).iter().any(|&j| kept[j])).collect();
    let narrow: Vec<bool> = (0..nn).map(|k| sea[k] && !opened[k]).collect();
    let one: Vec<bool> = (0..nn)
        .map(|k| {
            let [e, w, nth, s] = nb4_64(k);
            sea[k] && ((land[e] && land[w]) || (land[nth] && land[s]))
        })
        .collect();
    Channels { chans: comps8(&narrow, N64), ones: comps8(&one, N64), narrow, one }
}
/// The geodesic diameter (cells, 8-connected, two sweeps) of a component.
fn diameter(cells: &[usize], n: usize) -> usize {
    let set: std::collections::HashSet<usize> = cells.iter().copied().collect();
    let sweep = |s: usize| -> (usize, usize) {
        let mut d: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        d.insert(s, 0);
        let mut q = std::collections::VecDeque::from([s]);
        let mut far = (s, 0);
        while let Some(k) = q.pop_front() {
            let dk = d[&k];
            if dk > far.1 {
                far = (k, dk);
            }
            for j in nb8_n(k, n) {
                if set.contains(&j) && !d.contains_key(&j) {
                    d.insert(j, dk + 1);
                    q.push_back(j);
                }
            }
        }
        far
    };
    let (f, _) = sweep(cells[0]);
    sweep(f).1 + 1
}
/// B4's origin per cell: 1 rift, 2 ocean trapped by a merge, 3 the seaward band's remnant, 4 other.
fn origin_b4(o: &Obs) -> Vec<u8> {
    let pid = o.st.plate_id.data();
    let pt = o.st.plate_type.data();
    (0..N64 * N64)
        .map(|k| {
            if o.rift_sum[k] < 0.0 {
                1
            } else if pt[k] == PlateType::Oceanic && !nb4_64(k).iter().any(|&j| pid[j] != pid[k]) && o.merged_ocean[k] {
                2
            } else if o.band[k] {
                3
            } else {
                4
            }
        })
        .collect()
}
fn majority(cells: &[usize], org: &[u8]) -> u8 {
    let mut c = [0usize; 5];
    for &k in cells {
        c[org[k] as usize] += 1;
    }
    (1..5u8).max_by_key(|&i| (c[i as usize], std::cmp::Reverse(i))).unwrap()
}
const B4_NAMES: [&str; 5] = ["", "rift", "trapped by a merge", "band remnant", "other"];

/// B3: the origin per C1 cell (1 Andean, 2 collision, 3 O-O arcs, 0 other).
fn origin64(o: &Obs) -> Vec<u8> {
    (0..N64 * N64)
        .map(|k| {
            let v = [o.a_sum[k], o.c_sum[k], o.o_sum[k]];
            let (i, &mx) = v.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap().then(b.0.cmp(&a.0))).unwrap();
            if mx > 0.0 { i as u8 + 1 } else { 0 }
        })
        .collect()
}
const ORIGIN_NAMES: [&str; 4] = ["other", "Andean (C-O)", "collision (C-C)", "O-O arcs"];

/// The A1 map (64², ×6, north up): continental grey, active walk red, active 3 × 3 dark red, passive green (3 × 3 dark green),
/// enclosed yellow, the reached trench black, ocean blue.
fn a1_map_png(a: &A1, pt: &[PlateType]) -> image::RgbaImage {
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..N64 * N64 {
        let (x, y) = (k % N64, k / N64);
        let c = if a.trench[k] {
            [0, 0, 0, 255]
        } else {
            match a.cls[k] {
                1 => [230, 30, 30, 255],
                3 => [140, 10, 10, 255],
                2 => [40, 190, 70, 255],
                4 => [20, 110, 40, 255],
                5 => [240, 210, 40, 255],
                _ if pt[k] == PlateType::Continental => [190, 190, 180, 255],
                _ => [40, 70, 140, 255],
            }
        };
        for dy in 0..6 {
            for dx in 0..6 {
                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(c));
            }
        }
    }
    img
}
/// The channels map (64², ×6): land tan, sea blue, channel cells by class (rift magenta, trapped cyan, band yellow, other
/// white), one-cell-wide cells darker.
fn channels_png(land: &[bool], ch: &Channels, org: &[u8]) -> image::RgbaImage {
    let mut cls = vec![0u8; N64 * N64];
    for c in &ch.chans {
        let m = majority(c, org);
        for &k in c {
            cls[k] = m;
        }
    }
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..N64 * N64 {
        let (x, y) = (k % N64, k / N64);
        let mut c = if land[k] {
            [200, 180, 140]
        } else {
            match cls[k] {
                1 => [220, 40, 200],
                2 => [40, 220, 230],
                3 => [240, 220, 40],
                4 => [250, 250, 250],
                _ => [40, 70, 140],
            }
        };
        if ch.one[k] {
            c = [c[0] / 2 + 20, c[1] / 2 + 20, c[2] / 2 + 20];
        }
        for dy in 0..6 {
            for dx in 0..6 {
                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba([c[0], c[1], c[2], 255]));
            }
        }
    }
    img
}
/// The relief (256², grey shade, sea blue) with the mountain cells coloured by origin (Andean red, collision blue, O-O arcs
/// orange, other purple); `form` (if given) overrides: 1 cordillera red, 2 Andean not touching pink, 3 collision blue, 4 separated
/// collision navy, 5 arc orange, 6 other purple.
fn origin256_png(zm: &GridF32, classes: &[u8], org64: &[u8], form: Option<&[u8]>) -> image::RgbaImage {
    let n = 256;
    let mut img = image::RgbaImage::new(n as u32, n as u32);
    for k in 0..n * n {
        let (x, y) = (k % n, k / n);
        let z = zm.data[k];
        let mut px = if z <= 0.0 {
            let t = (-z / 5500.0).clamp(0.0, 1.0).sqrt();
            [(150.0 - 130.0 * t) as u8, (200.0 - 160.0 * t) as u8, (240.0 - 120.0 * t) as u8]
        } else {
            let t = (z / 3000.0).clamp(0.0, 1.0);
            let g = (110.0 + 130.0 * t) as u8;
            [g, g, g]
        };
        if classes[k] == MOUNTAIN && z > 0.0 {
            px = match form {
                Some(f) => match f[k] {
                    1 => [220, 30, 30],
                    2 => [240, 150, 170],
                    3 => [50, 90, 230],
                    4 => [20, 20, 110],
                    5 => [240, 150, 30],
                    _ => [150, 60, 170],
                },
                None => match org64[(y / 4) * 64 + x / 4] {
                    1 => [220, 30, 30],
                    2 => [50, 90, 230],
                    3 => [240, 150, 30],
                    _ => [150, 60, 170],
                },
            };
        }
        img.put_pixel(x as u32, (n - 1 - y) as u32, image::Rgba([px[0], px[1], px[2], 255]));
    }
    img
}
fn montage(tiles: &[image::RgbaImage], cols: usize, tw: u32, th: u32) -> image::RgbaImage {
    let rows = tiles.len().div_ceil(cols);
    let mut img = image::RgbaImage::from_pixel(cols as u32 * tw, rows as u32 * th, image::Rgba([255, 255, 255, 255]));
    for (i, t) in tiles.iter().enumerate() {
        image::imageops::overlay(&mut img, t, ((i % cols) as u32 * tw + 2) as i64, ((i / cols) as u32 * th + 2) as i64);
    }
    img
}

/// One C1 analysis (A1, B1, B2, B4, B5) of a world: the human lines, `#` data lines, and the A1 / channels images.
fn c1_177(seed: u64, p: C1Profile, m: &Margin, tag: &str) -> String {
    let clk = Clock::start();
    let o = observe(seed, p, m);
    let w = world_oo(seed, p, m);
    let cost = clk.stop();
    let nn = N64 * N64;
    let info_end = step_info(&o.st, &o.kin, o.ds_max, m.oo.is_some());
    let a = a1(&info_end, &o.kin, &o.tracker);
    let info_t0 = step_info(&o.st0, &o.kin0, o.ds_max, m.oo.is_some());
    let a0 = a1(&info_t0, &o.kin0, &OoTracker::default());
    let pt = o.st.plate_type.data();
    let (comp, big) = continents64(pt);
    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
    let act = (0..nn).filter(|&k| in_cont(k) && is_active_cls(a.cls[k])).count();
    let pas = (0..nn).filter(|&k| in_cont(k) && is_passive_cls(a.cls[k])).count();
    let enc = (0..nn).filter(|&k| in_cont(k) && a.cls[k] == 5).count();
    let act3 = (0..nn).filter(|&k| in_cont(k) && a.cls[k] == 3).count();
    let share = act as f32 / (act + pas).max(1) as f32;
    let act_mask: Vec<bool> = (0..nn).map(|k| in_cont(k) && is_active_cls(a.cls[k])).collect();
    let gpc = groups_per_continent(&act_mask, &comp, &big);
    let cont_active = gpc.iter().filter(|x| x.1 > 0).count();
    let newly: Vec<bool> = (0..nn).map(|k| act_mask[k] && is_passive_cls(a0.cls[k])).collect();
    let gpc_new = groups_per_continent(&newly, &comp, &big);
    let cont_new = gpc_new.iter().filter(|x| x.1 > 0).count();
    // the distance distribution (walk depth, active walk cells), cells
    let mut hist = vec![0usize; 12];
    let mut dists = Vec::new();
    for k in 0..nn {
        if in_cont(k) && a.cls[k] == 1 {
            hist[(a.dist[k] as usize).min(11)] += 1;
            dists.push(a.dist[k]);
        }
    }
    dists.sort_unstable();
    let med_walk = dists.get(dists.len() / 2).copied().unwrap_or(0);
    // the median over every active cell, the 3 × 3 ones at 1 cell
    let mut all_d: Vec<u32> = dists.clone();
    all_d.extend(std::iter::repeat_n(1, act3));
    all_d.sort_unstable();
    let med_all = all_d.get(all_d.len() / 2).copied().unwrap_or(0);
    // B2: the band (continental cells within Davis-Suppe's reach of a coast cell, assigned to the nearest by a 4-connected walk)
    let airy = Airy::declared();
    let dh = |k: usize, ds: f64| (airy.altitude_km(o.s0[k] + ds) - airy.altitude_km(o.s0[k])) * 1000.0;
    let mut lab = vec![0u8; nn];
    let mut dd = vec![u32::MAX; nn];
    let mut q = std::collections::VecDeque::new();
    for k in 0..nn {
        if in_cont(k) && (is_active_cls(a.cls[k]) || is_passive_cls(a.cls[k])) {
            lab[k] = if is_active_cls(a.cls[k]) { 1 } else { 2 };
            dd[k] = 0;
            q.push_back(k);
        }
    }
    while let Some(k) = q.pop_front() {
        if dd[k] as f64 + 1.0 > o.ds_max {
            continue;
        }
        for j in nb4_64(k) {
            if pt[j] == PlateType::Continental && dd[j] == u32::MAX {
                dd[j] = dd[k] + 1;
                lab[j] = lab[k];
                q.push_back(j);
            }
        }
    }
    let band_stat = |l: u8| -> (usize, f64, f64, f64) {
        let cells: Vec<usize> = (0..nn).filter(|&k| lab[k] == l).collect();
        let n = cells.len().max(1) as f64;
        let mean_dh = cells.iter().map(|&k| dh(k, o.src_sum[k])).sum::<f64>() / n;
        let mean_a = cells.iter().map(|&k| dh(k, o.a_sum[k])).sum::<f64>() / n;
        let pos: f64 = cells.iter().map(|&k| o.a_sum[k] + o.c_sum[k] + o.o_sum[k]).sum();
        let a_share = cells.iter().map(|&k| o.a_sum[k]).sum::<f64>() / pos.max(1e-12);
        (cells.len(), mean_dh, mean_a, a_share)
    };
    let (na, dha, aa, sha) = band_stat(1);
    let (np, dhp, ap, shp) = band_stat(2);
    // B4
    let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let land_f176: Vec<bool> = w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
    let land_diff = (0..nn).filter(|&k| land64[k] != land_f176[k]).count();
    let ch = channels(&land64);
    let org = origin_b4(&o);
    let mut by = [0usize; 5];
    let mut by_len = [0f32; 5];
    let mut ones_by = [0usize; 5];
    let mut chan_list = Vec::new();
    for c in &ch.chans {
        let mj = majority(c, &org);
        by[mj as usize] += 1;
        let len = diameter(c, N64) as f32 * CELL64_KM;
        by_len[mj as usize] += len;
        let (cx, cy) = (c.iter().map(|&k| k % N64).sum::<usize>() / c.len(), c.iter().map(|&k| k / N64).sum::<usize>() / c.len());
        chan_list.push(format!("({cx},{cy}) {} cells {len:.0} km {}", c.len(), B4_NAMES[mj as usize]));
    }
    for c in &ch.ones {
        ones_by[majority(c, &org) as usize] += 1;
    }
    // B5
    let mut sizes: std::collections::BTreeMap<u16, usize> = std::collections::BTreeMap::new();
    for &pidv in o.st.plate_id.data() {
        *sizes.entry(pidv).or_insert(0) += 1;
    }
    let mut sz: Vec<usize> = sizes.values().copied().collect();
    sz.sort_unstable_by(|a, b| b.cmp(a));
    let small = sz.iter().filter(|&&c| (c as f32) < 0.01 * nn as f32).count();
    let large = sz.iter().filter(|&&c| (c as f32) > 0.05 * nn as f32).count();
    let mid = sz.len() - small - large;
    a1_map_png(&a, pt).save(out177().join(format!("a1_{tag}.png"))).unwrap();
    channels_png(&land64, &ch, &org).save(out177().join(format!("channels_{tag}.png"))).unwrap();
    let mut s = format!(
        "   seed {seed} {tag} ({p:?}): A1 active {act} / passive {pas} (of which 3 × 3 active {act3}) / enclosed {enc} → share {:.1} % · mixed ties {} · continents {} with an active margin {cont_active} (newly active, F176's reading: {cont_new}) · groups per continent (id, groups, cells) {:?} · walk depth: median {med_walk} cells, all active (3 × 3 at 1) median {med_all} cells ({:.2} km world, {:.0} km terrestrial) · histogram (cells 0..11+) {:?}\n      B2: active band {na} cells, Δh_src {dha:.0} m, of which the C-O wedge {aa:.0} m ({:.0} % of the positive sources) · passive band {np} cells, Δh_src {dhp:.0} m (C-O {ap:.0} m, {:.0} %) · ratio active / passive {:.2}\n      B4: channels ≤ 2 cells {} (rift {}, trapped {}, band {}, other {}; length km {:.0} / {:.0} / {:.0} / {:.0}) · one-cell-wide groups {} (rift {}, trapped {}, band {}, other {}) · land 64² raw > 0 vs snaps > 0.5 differ on {land_diff} cells\n      B4 channels: {}\n      B5: plates at the end {} (< 1 % {small}, 1–5 % {mid}, > 5 % {large}) · sizes {:?} · merges {} · splits {}\n      C1 CPU {:.1} s\n",
        100.0 * share, a.mixed, big.len(), gpc, med_all as f32 * CELL64_KM, med_all as f32 * CELL64_KM * 7.5, hist,
        100.0 * sha, 100.0 * shp, dha / dhp,
        ch.chans.len(), by[1], by[2], by[3], by[4], by_len[1], by_len[2], by_len[3], by_len[4], ch.ones.len(), ones_by[1], ones_by[2], ones_by[3], ones_by[4],
        chan_list.join("; "), sz.len(), sz, o.merges, o.splits, cost.1
    );
    s += &format!("#A1 {} {act} {pas} {enc} {act3} {} {} {cont_active} {cont_new} {med_walk} {med_all}\n", share, a.mixed, big.len());
    s += &format!("#HIST {}\n", hist.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" "));
    s += &format!("#B2 {na} {dha} {aa} {sha} {np} {dhp} {ap} {shp}\n");
    s += &format!("#B4 {} {} {} {} {} {} {} {} {} {}\n", ch.chans.len(), by[1], by[2], by[3], by[4], ch.ones.len(), ones_by[1], ones_by[2], ones_by[3], ones_by[4]);
    s += &format!("#B5 {} {small} {mid} {large} {} {}\n", sz.len(), o.merges, o.splits);
    s
}

/// The selection under F176's rule, cached: (plates, k, start, compact) or None.
fn selection177(seed: u64) -> Option<C1Profile> {
    let s = cached177(&format!("sel_seed{seed}"), || {
        let m = Margin::seaward_oo_collision();
        let (np, tried) = select_continent_default_margin(seed, &m);
        match np {
            Some(np) => {
                let sel = &tried.last().unwrap().1;
                let cd = sel.candidates[sel.chosen.unwrap()];
                match cd.profile() {
                    C1Profile::Selected { plates, k, start, compact } => format!("#SEL {plates} {k} {start} {}\n   seed {seed}: selection N {np}, land {:.1} %, r* {:?}\n", compact as u8, 100.0 * cd.land, cd.circ.sea_r),
                    other => format!("#SELX\n   seed {seed}: unexpected profile {other:?}\n"),
                }
            }
            None => "#SELNONE\n".to_string(),
        }
    });
    let v = data_line(&s, "SEL")?;
    Some(C1Profile::Selected { plates: v[0] as usize, k: v[1] as usize, start: v[2] as u16, compact: v[3] != 0.0 })
}
fn panel_and_screen() -> Vec<u64> {
    PANEL.iter().chain(SCREEN.iter()).copied().collect()
}

/// F177 tier 1 -- C1: the 15 seeds under F176's selection (A1, B1, B2, B4, B5) and the panel at the fixed continent (B4's straits).
#[allow(dead_code)]
fn f177_tier1_c1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out177()).unwrap();
    eprintln!("\n==========  F177 tier 1 . C1: A1, the uplift, the channels, the plates  ==========");
    let seeds = panel_and_screen();
    // the selections (parallel, cached)
    std::thread::scope(|sc| {
        for &seed in &seeds {
            sc.spawn(move || {
                let _ = selection177(seed);
            });
        }
    });
    let m = Margin::seaward_oo_collision();
    std::thread::scope(|sc| {
        for &seed in &seeds {
            let m = m;
            sc.spawn(move || {
                if let Some(p) = selection177(seed) {
                    let tag = format!("sel_seed{seed}");
                    cached177(&format!("t1_{tag}"), || c1_177(seed, p, &m, &tag));
                }
            });
        }
        for seed in PANEL {
            let m = m;
            sc.spawn(move || {
                let tag = format!("fixed_seed{seed}");
                cached177(&format!("t1_{tag}"), || c1_177(seed, fixed_profile(seed), &m, &tag));
            });
        }
    });
    let mut shares = Vec::new();
    let (mut conts, mut conts_act, mut conts_new) = (0.0, 0.0, 0.0);
    let mut hist = vec![0f64; 12];
    let mut b4 = vec![0f64; 10];
    let mut b5 = vec![0f64; 6];
    let mut walk_meds = Vec::new();
    let mut all_meds = Vec::new();
    for &seed in &seeds {
        let sel = std::fs::read_to_string(results177().join(format!("sel_seed{seed}.txt"))).unwrap_or_default();
        eprint!("{}", sel.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        let Ok(s) = std::fs::read_to_string(results177().join(format!("t1_sel_seed{seed}.txt"))) else { continue };
        eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        let v = data_line(&s, "A1").unwrap();
        shares.push(v[0] * 100.0);
        conts += v[6];
        conts_act += v[7];
        conts_new += v[8];
        walk_meds.push(v[9]);
        all_meds.push(v[10]);
        for (h, x) in hist.iter_mut().zip(data_line(&s, "HIST").unwrap()) {
            *h += x;
        }
        for (h, x) in b4.iter_mut().zip(data_line(&s, "B4").unwrap()) {
            *h += x;
        }
        for (h, x) in b5.iter_mut().zip(data_line(&s, "B5").unwrap()) {
            *h += x;
        }
    }
    let med = |v: &mut Vec<f64>| -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if v.is_empty() { f64::NAN } else { v[v.len() / 2] }
    };
    let mut sh = shares.clone();
    eprintln!(
        "\n   B1 (15 seeds, the selection): the A1 share median {:.1} %, sorted {:?} · continents with an active margin {conts_act} / {conts} ({:.0} %), newly active (F176's reading) {conts_new} · the walk depth medians (cells) {:?}, all active {:?} · the walk depth histogram (cells 0..11+) {:?}\n   B4 (15 seeds): channels ≤ 2 cells {} (rift {}, trapped {}, band {}, other {}) · one-cell-wide groups {} (rift {}, trapped {}, band {}, other {})\n   B5 (15 seeds): plates {} (< 1 % {}, 1–5 % {}, > 5 % {}) · merges {} · splits {}",
        med(&mut sh), shares.iter().map(|x| (x * 10.0).round() / 10.0).collect::<Vec<_>>(), 100.0 * conts_act / conts.max(1.0), walk_meds, all_meds, hist,
        b4[0], b4[1], b4[2], b4[3], b4[4], b4[5], b4[6], b4[7], b4[8], b4[9], b5[0], b5[1], b5[2], b5[3], b5[4], b5[5]
    );
    for seed in PANEL {
        if let Ok(s) = std::fs::read_to_string(results177().join(format!("t1_fixed_seed{seed}.txt"))) {
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        }
    }
    // the A1 maps of the panel (the selection), the channels maps of the panel (the fixed continent)
    let tiles: Vec<image::RgbaImage> = PANEL.iter().filter_map(|s| image::open(out177().join(format!("a1_sel_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 388, 388).save(out177().join("a1_panel_selection.png")).unwrap();
    let tiles: Vec<image::RgbaImage> = PANEL.iter().filter_map(|s| image::open(out177().join(format!("channels_fixed_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 388, 388).save(out177().join("channels_panel_fixed.png")).unwrap();
    let tiles: Vec<image::RgbaImage> = panel_and_screen().iter().filter_map(|s| image::open(out177().join(format!("channels_sel_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 388, 388).save(out177().join("channels_15_selection.png")).unwrap();
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// The physics level at 256², young, for one world: (Run, the world).
fn physics_young(seed: u64, p: C1Profile, m: &Margin, rough_a: f32) -> (Run, World) {
    let kh = KHetero::declared();
    let w = world_oo(seed, p, m);
    let cover = cover_map_m(seed, p, m, COVER_M);
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k) };
    let cfg = base_cfg(seed, rough_a, Some([0, 0]));
    let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[1].1, &rock, N_PHYS);
    (r, w)
}
fn rough_a177() -> f32 {
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    0.5 * read_cell(&c128, None, CELL).0.octaves[0]
}

/// F177 tier 2 -- B3 on the panel at the fixed continent: the mountain cells at young by origin, F175 against F176.
#[allow(dead_code)]
fn f177_tier2_attribution() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out177()).unwrap();
    eprintln!("\n==========  F177 tier 2 . B3: the mountain gain by origin, the panel at the fixed continent, young  ==========");
    let rough_a = rough_a177();
    let km2 = CELL * CELL;
    for seed in PANEL {
        let p = fixed_profile(seed);
        for (name, m) in [("F175", Margin::seaward_oo()), ("F176", Margin::seaward_oo_collision())] {
            let key = format!("t2_seed{seed}_{name}");
            let s = cached177(&key, || {
                let o = observe(seed, p, &m);
                let org = origin64(&o);
                let (r, _) = physics_young(seed, p, &m, rough_a);
                let mut cnt = [0usize; 4];
                for k in 0..65536 {
                    if r.classes[k] == MOUNTAIN && r.zm.data[k] > 0.0 {
                        cnt[org[((k / 256) / 4) * 64 + (k % 256) / 4] as usize] += 1;
                    }
                }
                let tot: usize = cnt.iter().sum();
                origin256_png(&r.zm, &r.classes, &org, None).save(out177().join(format!("origin_seed{seed}_{name}.png"))).unwrap();
                format!(
                    "   seed {seed} {name} young: mountain {:.1} % of the land ({:.0} km²) · by origin: Andean {:.0} · collision {:.0} · O-O arcs {:.0} · other {:.0} km² · C1 cells by origin (A / C / O / other) {} / {} / {} / {}\n#B3 {} {} {} {} {}\n",
                    100.0 * r.m.fractions[3], tot as f32 * km2, cnt[1] as f32 * km2, cnt[2] as f32 * km2, cnt[3] as f32 * km2, cnt[0] as f32 * km2,
                    org.iter().filter(|&&v| v == 1).count(), org.iter().filter(|&&v| v == 2).count(), org.iter().filter(|&&v| v == 3).count(), org.iter().filter(|&&v| v == 0).count(),
                    tot, cnt[1], cnt[2], cnt[3], cnt[0]
                )
            });
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        }
    }
    let mut collision_major = 0;
    for seed in PANEL {
        let g = |name: &str| data_line(&std::fs::read_to_string(results177().join(format!("t2_seed{seed}_{name}.txt"))).unwrap(), "B3").unwrap();
        let (a, b) = (g("F175"), g("F176"));
        let d: Vec<f64> = (0..5).map(|i| (b[i] - a[i]) * km2 as f64).collect();
        let pos: f64 = d[1..].iter().filter(|x| **x > 0.0).sum();
        if d[0] > 0.0 && d[2] > 0.5 * d[0] {
            collision_major += 1;
        }
        eprintln!(
            "   seed {seed} the gain F176 − F175: total {:+.0} km² · Andean {:+.0} · collision {:+.0} · O-O arcs {:+.0} · other {:+.0} km² · collision's share of the net gain {:.0} % (of the positive gains {:.0} %)",
            d[0], d[1], d[2], d[3], d[4], 100.0 * d[2] / d[0], 100.0 * d[2].max(0.0) / pos.max(1e-9)
        );
    }
    eprintln!("   T3: collision gives > half of a positive gain on {collision_major} / 5 seeds");
    let mut tiles = Vec::new();
    for seed in PANEL {
        for name in ["F175", "F176"] {
            if let Ok(i) = image::open(out177().join(format!("origin_seed{seed}_{name}.png"))) {
                tiles.push(i.to_rgba8());
            }
        }
        if let Ok(i) = image::open(out177().join(format!("channels_fixed_seed{seed}.png"))) {
            tiles.push(image::imageops::resize(&i.to_rgba8(), 256, 256, image::imageops::FilterType::Nearest));
        }
    }
    montage(&tiles, 3, 260, 260).save(out177().join("attribution_panel_F175_F176_channels.png")).unwrap();
    eprintln!("   tier 2: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F177 tier 3 -- the formations on the 15 seeds under F176's selection, young.
#[allow(dead_code)]
fn f177_tier3_formations() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out177()).unwrap();
    eprintln!("\n==========  F177 tier 3 . the formations, 15 seeds, the selection, young  ==========");
    let rough_a = rough_a177();
    let km2 = CELL * CELL;
    let m = Margin::seaward_oo_collision();
    for seed in panel_and_screen() {
        let Some(p) = selection177(seed) else { continue };
        let s = cached177(&format!("t3_seed{seed}"), || {
            let o = observe(seed, p, &m);
            let org = origin64(&o);
            let info_end = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info_end, &o.kin, &o.tracker);
            let (r, w) = physics_young(seed, p, &m, rough_a);
            let n = 256;
            let org256 = |k: usize| org[((k / n) / 4) * 64 + (k % n) / 4];
            let mtn: Vec<bool> = (0..n * n).map(|k| r.classes[k] == MOUNTAIN && r.zm.data[k] > 0.0).collect();
            // the land components at 256² (8-connected); the largest
            let land: Vec<bool> = r.zm.data.iter().map(|&z| z > 0.0).collect();
            let lc = comps8(&land, n);
            let mut lab = vec![usize::MAX; n * n];
            for (i, c) in lc.iter().enumerate() {
                for &k in c {
                    lab[k] = i;
                }
            }
            let main = (0..lc.len()).max_by_key(|&i| lc[i].len()).unwrap_or(0);
            // B4's channels at the end of C1 (for the separated highlands' nearest channel)
            let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
            let ch = channels(&land64);
            let org4 = origin_b4(&o);
            let mut chan_cls = vec![0u8; N64 * N64];
            for c in &ch.chans {
                let mj = majority(c, &org4);
                for &k in c {
                    chan_cls[k] = mj;
                }
            }
            let mut form = vec![0u8; n * n];
            // the Andean groups and the cordilleras
            let and_mask: Vec<bool> = (0..n * n).map(|k| mtn[k] && org256(k) == 1).collect();
            let mut cords = Vec::new();
            for g in comps8(&and_mask, n) {
                let mut touched = std::collections::BTreeSet::new();
                for &k in &g {
                    let mut cand = vec![k];
                    cand.extend(nb8_n(k, n));
                    for c in cand {
                        let c64 = ((c / n) / 4) * 64 + (c % n) / 4;
                        if is_active_cls(a.cls[c64]) {
                            touched.insert(c64);
                        }
                    }
                }
                let f = if touched.is_empty() { 2 } else { 1 };
                for &k in &g {
                    form[k] = f;
                }
                if f == 1 {
                    cords.push((g.len(), touched.len()));
                }
            }
            // collision (separated or not), arcs, other
            let mut sep_classes = Vec::new();
            let col_mask: Vec<bool> = (0..n * n).map(|k| mtn[k] && org256(k) == 2).collect();
            for g in comps8(&col_mask, n) {
                let sep = lab[g[0]] != main;
                for &k in &g {
                    form[k] = if sep { 4 } else { 3 };
                }
                if sep {
                    // the nearest channel cell (64²) to the group
                    let mut best = (i64::MAX, 0u8);
                    for &k in &g {
                        let (x, y) = (((k % n) / 4) as i64, ((k / n) / 4) as i64);
                        for q in 0..N64 * N64 {
                            if chan_cls[q] == 0 {
                                continue;
                            }
                            let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
                            dx = dx.min(64 - dx);
                            let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
                            dy = dy.min(64 - dy);
                            if dx * dx + dy * dy < best.0 {
                                best = (dx * dx + dy * dy, chan_cls[q]);
                            }
                        }
                    }
                    sep_classes.push((g.len(), B4_NAMES[best.1 as usize], ((best.0 as f32).sqrt() * CELL64_KM) as i32));
                }
            }
            for k in 0..n * n {
                if mtn[k] && form[k] == 0 {
                    form[k] = if org256(k) == 3 { 5 } else { 6 };
                }
            }
            let area = |f: u8| form.iter().filter(|&&v| v == f).count() as f32 * km2;
            let (cord_km2, and_other, col, sep, arc, oth) = (area(1), area(2), area(3), area(4), area(5), area(6));
            let cord_len: f32 = cords.iter().map(|c| c.1 as f32 * CELL64_KM).sum();
            let only_collision = cord_km2 + and_other + arc + oth == 0.0 && col + sep > 0.0;
            origin256_png(&r.zm, &r.classes, &org, Some(&form)).save(out177().join(format!("formations_seed{seed}.png"))).unwrap();
            format!(
                "   seed {seed} ({p:?}) young: mountain {:.1} % of the land · Andean cordilleras {} (lengths along the coast km {:?}; total {cord_len:.0} km world, {:.0} km terrestrial; {cord_km2:.0} km²) · Andean not touching an active coast {and_other:.0} km² · collisional highland {col:.0} km² · separated {sep:.0} km² (groups: cells, nearest channel class, km {:?}) · arc mountain on the band {arc:.0} km² · other {oth:.0} km² · collision the only formation: {only_collision}\n#F {} {cord_len} {cord_km2} {and_other} {col} {sep} {arc} {oth} {}\n",
                100.0 * r.m.fractions[3], cords.len(), cords.iter().map(|c| c.1 as f32 * CELL64_KM).collect::<Vec<_>>(), cord_len * 7.5, sep_classes,
                cords.len(), only_collision as u8
            )
        });
        eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
    }
    let mut with_cord = 0;
    let mut only = 0;
    let mut n = 0;
    for seed in panel_and_screen() {
        if let Ok(s) = std::fs::read_to_string(results177().join(format!("t3_seed{seed}.txt"))) {
            let v = data_line(&s, "F").unwrap();
            n += 1;
            if v[0] > 0.0 {
                with_cord += 1;
            }
            if v[8] > 0.0 {
                only += 1;
            }
        }
    }
    eprintln!("   (v): seeds with ≥ 1 Andean cordillera {with_cord} / {n} (≥ 5 required) · seeds whose only mountain formation is the collisional / separated highland {only} (0 required)");
    let tiles: Vec<image::RgbaImage> = panel_and_screen().iter().filter_map(|s| image::open(out177().join(format!("formations_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 260, 260).save(out177().join("formations_15.png")).unwrap();
    eprintln!("   tier 3: {:.0} s", t0.elapsed().as_secs_f64());
}

/// A transect image without a clamp: the altitude's range from the data (top), s (bottom, 0 … max).
fn transect_png_free(z: &[(Vec<f32>, [u8; 3])], s: &[(Vec<f32>, [u8; 3])], path: &std::path::Path) -> (f32, f32) {
    let (w, h) = (900u32, 600u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
    let zmin = z.iter().flat_map(|(v, _)| v.iter().copied()).fold(f32::INFINITY, f32::min);
    let zmax = z.iter().flat_map(|(v, _)| v.iter().copied()).fold(f32::NEG_INFINITY, f32::max);
    let smax = s.iter().flat_map(|(v, _)| v.iter().copied()).fold(0f32, f32::max).max(1.0);
    let (lo, hi) = (zmin - 100.0, zmax + 100.0);
    let mut put = |x: i64, y: i64, c: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as u32) < w && (y as u32) < h {
            img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], 255]));
        }
    };
    let pz = |v: f32| (20.0 + (hi - v) / (hi - lo) * 340.0) as i64;
    let ps = |v: f32| (400.0 + (smax - v) / smax * 180.0) as i64;
    for x in 20..w as i64 - 20 {
        put(x, pz(0.0), [0, 0, 0]);
        let mut t = (lo / 1000.0).ceil() * 1000.0;
        while t < hi {
            put(x, pz(t), [215, 215, 215]);
            t += 1000.0;
        }
        put(x, ps(1.0), [0, 0, 0]);
        put(x, ps(0.915), [160, 160, 230]);
    }
    for (series, f) in [(z, &pz as &dyn Fn(f32) -> i64), (s, &ps as &dyn Fn(f32) -> i64)] {
        for (v, c) in series {
            for (i, &q) in v.iter().enumerate() {
                let x = 20 + (i as f32 / (v.len() - 1).max(1) as f32 * (w as f32 - 40.0)) as i64;
                for d in -1..=1 {
                    put(x, f(q) + d, *c);
                }
            }
        }
    }
    img.save(path).unwrap();
    (zmin, zmax)
}

/// F177 tier 4 -- B7: F176's tier-4 transect of seed 18 (young, 2 048², the fixed continent), redone without a clamp.
#[allow(dead_code)]
fn f177_tier4_transect() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out177()).unwrap();
    eprintln!("\n==========  F177 tier 4 . B7: the transect of seed 18 without a clamp  ==========");
    let seed = 18;
    let m = Margin::seaward_oo_collision();
    let p = fixed_profile(seed);
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let kh = KHetero::declared();
    let s = cached177("t4_seed18", || {
        let w = world_oo(seed, p, &m);
        let ser = c1_series(seed, p, &m);
        let cover = cover_map_m(seed, p, &m, COVER_M);
        let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
        let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k.clone()) };
        let cfg = base_cfg(seed, rough_a, Some([0, 0]));
        let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
        let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[1].1, &rock, N_PHYS);
        let (g, _, _) = chain_low(&r.zm, amp, Some(ln_k), 0.0);
        let act = (0..4096).find(|&k| ser.map_end[k] == 3 && ser.map_t0[k] == 2).or_else(|| (0..4096).find(|&k| ser.map_end[k] == 3));
        let pas = (0..4096).find(|&k| ser.map_end[k] == 2);
        let s64 = w.st.s.data().to_vec();
        let pt = w.st.plate_type.data().to_vec();
        let mut zs = Vec::new();
        let mut ss = Vec::new();
        let mut desc = Vec::new();
        for (name, cell, col) in [("passive", pas, [30u8, 140, 60]), ("Andean", act, [220, 30, 30])] {
            if let Some(k) = cell {
                let (a, b) = transect(&g, &s64, &pt, k);
                let mx = a.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                let at_max = a.iter().filter(|&&v| (v - mx).abs() < 1.0).count();
                let mut run = 0;
                let mut best = 0;
                for wv in a.windows(2) {
                    if (wv[1] - wv[0]).abs() < 0.01 {
                        run += 1;
                        best = best.max(run);
                    } else {
                        run = 0;
                    }
                }
                desc.push(format!("{name} (cell {k}): max {mx:.0} m, samples within 1 m of the max {at_max} of {}, the longest run of equal values (|Δ| < 0.01 m) {best} samples", a.len()));
                zs.push((a, col));
                ss.push((b, col));
            }
        }
        let (zmin, zmax) = transect_png_free(&zs, &ss, &out177().join("transect_seed18_noclip.png"));
        let peak = g.data.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        format!("   seed 18 F176 young (2 048²): the transect's range {zmin:.0} … {zmax:.0} m · {} · the world's 2 048² maximum {peak:.0} m\n", desc.join(" · "))
    });
    eprint!("{s}");
    eprintln!("   tier 4: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F177-B2's stop, its mechanism (C1 only, the 15 seeds under the selection): on the active band, the mean Δh by origin (Andean /
/// collision / O-O arcs / the rest), and the active coasts by kind: on a boundary (3 × 3), walk within the arc's reach
/// (`arc_distance` = 3 cells), walk beyond it.
#[allow(dead_code)]
fn f177_b2_mechanism() {
    let m = Margin::seaward_oo_collision();
    let arc_d = C1Closures::default().subduction.arc_distance as u32;
    let seeds = panel_and_screen();
    std::thread::scope(|sc| {
        for &seed in &seeds {
            sc.spawn(move || {
                let Some(p) = selection177(seed) else { return };
                cached177(&format!("b2mech_seed{seed}"), || {
                    let o = observe(seed, p, &m);
                    let info = step_info(&o.st, &o.kin, o.ds_max, true);
                    let a = a1(&info, &o.kin, &o.tracker);
                    let pt = o.st.plate_type.data();
                    let (comp, big) = continents64(pt);
                    let nn = N64 * N64;
                    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
                    let airy = Airy::declared();
                    let dh = |k: usize, ds: f64| (airy.altitude_km(o.s0[k] + ds) - airy.altitude_km(o.s0[k])) * 1000.0;
                    let mut lab = vec![0u8; nn];
                    let mut dd = vec![u32::MAX; nn];
                    let mut q = std::collections::VecDeque::new();
                    for k in 0..nn {
                        if in_cont(k) && (is_active_cls(a.cls[k]) || is_passive_cls(a.cls[k])) {
                            lab[k] = if is_active_cls(a.cls[k]) { 1 } else { 2 };
                            dd[k] = 0;
                            q.push_back(k);
                        }
                    }
                    while let Some(k) = q.pop_front() {
                        if dd[k] as f64 + 1.0 > o.ds_max {
                            continue;
                        }
                        for j in nb4_64(k) {
                            if pt[j] == PlateType::Continental && dd[j] == u32::MAX {
                                dd[j] = dd[k] + 1;
                                lab[j] = lab[k];
                                q.push_back(j);
                            }
                        }
                    }
                    let cells: Vec<usize> = (0..nn).filter(|&k| lab[k] == 1).collect();
                    let nc = cells.len().max(1) as f64;
                    let mean = |f: &dyn Fn(usize) -> f64| cells.iter().map(|&k| f(k)).sum::<f64>() / nc;
                    let (ma, mc, mo) = (mean(&|k| dh(k, o.a_sum[k])), mean(&|k| dh(k, o.c_sum[k])), mean(&|k| dh(k, o.o_sum[k])));
                    let tot = mean(&|k| dh(k, o.src_sum[k]));
                    let act: Vec<usize> = (0..nn).filter(|&k| in_cont(k) && is_active_cls(a.cls[k])).collect();
                    let on_b = act.iter().filter(|&&k| a.cls[k] == 3).count();
                    let within = act.iter().filter(|&&k| a.cls[k] == 1 && a.dist[k] <= arc_d).count();
                    let beyond = act.iter().filter(|&&k| a.cls[k] == 1 && a.dist[k] > arc_d).count();
                    // the coast cells of the active band that received any O-O arc, any C-O deposit or wedge
                    let coast_oo = act.iter().filter(|&&k| o.o_sum[k] > 0.0).count();
                    let coast_a = act.iter().filter(|&&k| o.a_sum[k] > 0.0).count();
                    format!(
                        "   seed {seed}: active band {} cells: Δh_src {tot:.0} m = Andean {ma:.0} + collision {mc:.0} + O-O arcs {mo:.0} + the rest {:.0} m · active coast cells {}: on a boundary (3 × 3) {on_b}, walk within the arc's reach (≤ {arc_d} cells) {within}, beyond it {beyond} · active coast cells with an O-O arc deposit {coast_oo}, with a C-O wedge or arc {coast_a}\n#B2M {tot} {ma} {mc} {mo} {} {on_b} {within} {beyond} {coast_oo} {coast_a}\n",
                        cells.len(), tot - ma - mc - mo, act.len(), act.len()
                    )
                });
            });
        }
    });
    for seed in panel_and_screen() {
        if let Ok(s) = std::fs::read_to_string(results177().join(format!("b2mech_seed{seed}.txt"))) {
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        }
    }
}

// ════════════════════════════════════════ F178 ════════════════════════════════════════

use ymir_core::cascade::history::run_c1_sources_mapped_ext;
use ymir_core::cascade::margin::{ANDEAN_GAP_CELLS, ARC_BAND_CELLS};
use ymir_core::tectonics_c1::andean_seed::{andean_seeds, overriding_trench_cells};
use ymir_core::tectonics_c1::closures::subduction::OoSubductionParams;
use ymir_core::tectonics_c1::state::BoolField;
use ymir_core::tectonics_c1::time_loop::run_with_closures_terms_ext;

fn out178() -> PathBuf {
    root().join("docs/reports/relief_method/f178_andean_reach/images")
}
fn results178() -> PathBuf {
    root().join("docs/reports/relief_method/f178_andean_reach/results")
}
fn cached178(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results178().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results178()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
/// F176's selection, as F177 cached it.
fn sel176(seed: u64) -> Option<C1Profile> {
    let s = std::fs::read_to_string(root().join(format!("docs/reports/relief_method/f177_accretion_close/results/sel_seed{seed}.txt"))).ok()?;
    let v = data_line(&s, "SEL")?;
    Some(C1Profile::Selected { plates: v[0] as usize, k: v[1] as usize, start: v[2] as u16, compact: v[3] != 0.0 })
}
fn variants178() -> Vec<(&'static str, Margin)> {
    let f176 = Margin::seaward_oo_collision();
    let f178 = Margin::seaward_oo_collision_andean();
    let a2 = Margin { oo: f178.oo.map(|o| OoSubductionParams { arc_band: Some(ARC_BAND_CELLS), ..o }), ..f178 };
    vec![
        ("F176", f176),
        ("F178", f178),
        ("F178_A2", a2),
        ("F178_gap3", Margin { andean_gap: Some(3), ..f178 }),
        ("F178_gap5", Margin { andean_gap: Some(5), ..f178 }),
    ]
}
fn to_boolfield(v: &[bool]) -> BoolField {
    let mut b = BoolField::filled(N64, N64, false);
    for (k, &on) in v.iter().enumerate() {
        if on {
            b.set(k % N64, k / N64, true);
        }
    }
    b
}
/// The active trench at the end: the overriding-side O-O cells of subducting pairs and the oceanic convergent cells touching
/// another plate's continent.
fn active_trench(o: &Obs, info: &StepInfo) -> Vec<bool> {
    let tr = overriding_trench_cells(&o.st.plate_id, &o.st.plate_type, &o.st.age, &info.bi, &to_boolfield(&info.upper_vel), &o.kin, &o.tracker);
    let bt = info.bi.boundary_type.data();
    (0..N64 * N64).map(|k| tr.data()[k] || (info.co_cell[k] && bt[k] == BoundaryType::Convergent)).collect()
}
const BANDS: [(u32, u32); 4] = [(0, 3), (3, 6), (6, 15), (15, 30)];
/// The uplift by origin by distance band: [coast type: active, passive][band][origin: C-O wedge, collision, arcs, other] mean Δh (m),
/// and the cell counts [type][band].
fn banded(o: &Obs, a: &A1, comp: &[usize], big: &[usize]) -> ([[[f64; 4]; 4]; 2], [[usize; 4]; 2]) {
    let nn = N64 * N64;
    let pt = o.st.plate_type.data();
    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
    let airy = Airy::declared();
    let dh = |k: usize, ds: f64| (airy.altitude_km(o.s0[k] + ds) - airy.altitude_km(o.s0[k])) * 1000.0;
    let mut lab = vec![0u8; nn];
    let mut dd = vec![u32::MAX; nn];
    let mut q = std::collections::VecDeque::new();
    for k in 0..nn {
        if in_cont(k) && (is_active_cls(a.cls[k]) || is_passive_cls(a.cls[k])) {
            lab[k] = if is_active_cls(a.cls[k]) { 1 } else { 2 };
            dd[k] = 0;
            q.push_back(k);
        }
    }
    while let Some(k) = q.pop_front() {
        for j in nb4_64(k) {
            if pt[j] == PlateType::Continental && dd[j] == u32::MAX && dd[k] + 1 < 30 {
                dd[j] = dd[k] + 1;
                lab[j] = lab[k];
                q.push_back(j);
            }
        }
    }
    let mut sum = [[[0f64; 4]; 4]; 2];
    let mut cnt = [[0usize; 4]; 2];
    for k in 0..nn {
        if lab[k] == 0 {
            continue;
        }
        let t = (lab[k] - 1) as usize;
        let Some(b) = BANDS.iter().position(|&(lo, hi)| dd[k] >= lo && dd[k] < hi) else { continue };
        let (w, c, ar) = (dh(k, o.w_sum[k]), dh(k, o.c_sum[k]), dh(k, o.arc_sum[k]));
        let tot = dh(k, o.src_sum[k]);
        for (i, v) in [w, c, ar, tot - w - c - ar].into_iter().enumerate() {
            sum[t][b][i] += v;
        }
        cnt[t][b] += 1;
    }
    for t in 0..2 {
        for b in 0..4 {
            for i in 0..4 {
                sum[t][b][i] /= cnt[t][b].max(1) as f64;
            }
        }
    }
    (sum, cnt)
}

/// One F178 C1 analysis: A1 (both tie rules), the continents, the banded uplift, the trench segments, the land and the
/// circumnavigation.
fn c1_178(seed: u64, p: C1Profile, m: &Margin, tag: &str) -> String {
    let clk = Clock::start();
    let o = observe(seed, p, m);
    let w = world_oo(seed, p, m);
    let cost = clk.stop();
    let nn = N64 * N64;
    let info = step_info(&o.st, &o.kin, o.ds_max, true);
    let a = a1(&info, &o.kin, &o.tracker);
    let pt = o.st.plate_type.data();
    let (comp, big) = continents64(pt);
    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
    let share = |cls: &[u8]| -> f64 {
        let act = (0..nn).filter(|&k| in_cont(k) && is_active_cls(cls[k])).count();
        let pas = (0..nn).filter(|&k| in_cont(k) && is_passive_cls(cls[k])).count();
        act as f64 / (act + pas).max(1) as f64
    };
    let (s_any, s_maj) = (share(&a.cls), share(&a.major));
    let act_mask: Vec<bool> = (0..nn).map(|k| in_cont(k) && is_active_cls(a.cls[k])).collect();
    let cont_active = groups_per_continent(&act_mask, &comp, &big).iter().filter(|x| x.1 > 0).count();
    let (b, bc) = banded(&o, &a, &comp, &big);
    let mean06 = |t: usize, i: usize| (b[t][0][i] * bc[t][0] as f64 + b[t][1][i] * bc[t][1] as f64) / (bc[t][0] + bc[t][1]).max(1) as f64;
    let tot06 = |t: usize| (0..4).map(|i| mean06(t, i)).sum::<f64>();
    let wedge_act06 = mean06(0, 0);
    let (up_a, up_p) = (tot06(0), tot06(1));
    // the trench segments
    let tr = active_trench(&o, &info);
    let mut lens: Vec<f32> = comps8(&tr, N64).iter().map(|c| diameter(c, N64) as f32 * CELL64_KM).collect();
    lens.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let q = |f: f32| if lens.is_empty() { 0.0 } else { lens[((lens.len() - 1) as f32 * f) as usize] };
    // the land and the circumnavigation (C1's end land, as the selection)
    let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let circ = circumnavigation(&land64, 64);
    let land = land64.iter().filter(|&&l| l).count() as f32 / nn as f32;
    let names = ["C-O wedge", "collision", "arcs", "other"];
    let mut tab = String::new();
    for (t, tn) in ["active", "passive"].iter().enumerate() {
        for (bi_, (lo, hi)) in BANDS.iter().enumerate() {
            tab += &format!(
                "      {tn} {lo}–{hi} cells ({} cells): {}\n",
                bc[t][bi_],
                (0..4).map(|i| format!("{} {:.0}", names[i], b[t][bi_][i])).collect::<Vec<_>>().join(" · ")
            );
        }
    }
    a1_map_png(&a, pt).save(out178().join(format!("a1_{tag}.png"))).unwrap();
    let mut s = format!(
        "   seed {seed} {tag}: A1 share any {:.1} % · majority {:.1} % · continents {} with an active margin {cont_active} · on 0–6 cells: active Δh {up_a:.0} m (C-O wedge {wedge_act06:.0} m), passive {up_p:.0} m, ratio {:.2} · trench segments {} (median {:.0} km, p90 {:.0} km, max {:.0} km world) · C1 end land {:.1} %, circumnavigable {} r* {:?} · merges {} splits {} · C1 CPU {:.1} s\n{tab}",
        100.0 * s_any, 100.0 * s_maj, big.len(), up_a / up_p, lens.len(), q(0.5), q(0.9), lens.last().copied().unwrap_or(0.0), 100.0 * land, !circ.wraps, circ.sea_r, o.merges, o.splits, cost.1
    );
    s += &format!("#A1 {s_any} {s_maj} {} {cont_active}\n", big.len());
    s += &format!("#UP {up_a} {up_p} {wedge_act06}\n");
    s += &format!("#TR {} {} {} {}\n", lens.len(), q(0.5), q(0.9), lens.last().copied().unwrap_or(0.0));
    s += &format!("#LAND {land} {} {}\n", !circ.wraps as u8, circ.sea_r.unwrap_or(0));
    s += &format!("#BAND {}\n", (0..2).flat_map(|t| (0..4).flat_map(move |bb| (0..4).map(move |i| (t, bb, i)))).map(|(t, bb, i)| b[t][bb][i].to_string()).collect::<Vec<_>>().join(" "));
    s
}

fn seeds15() -> Vec<u64> {
    PANEL.iter().chain(SCREEN.iter()).copied().collect()
}

/// F178 tier 1 -- C1 on the 15 seeds (F176's selection): F176, F178, F178 + A2, the reach at 3 and 5 cells.
#[allow(dead_code)]
fn f178_tier1_c1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out178()).unwrap();
    eprintln!("\n==========  F178 tier 1 . C1: A1, the banded uplift, the trenches  ==========");
    std::thread::scope(|sc| {
        for seed in seeds15() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                for (name, m) in variants178() {
                    let tag = format!("{name}_seed{seed}");
                    cached178(&format!("t1_{tag}"), || c1_178(seed, p, &m, &tag));
                }
            });
        }
    });
    let med = |mut v: Vec<f64>| -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if v.is_empty() { f64::NAN } else { v[v.len() / 2] }
    };
    for (name, _) in variants178() {
        let (mut any, mut maj, mut ratio, mut wedge) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let (mut conts, mut conts_a, mut wpos, mut r2, mut circ_lost, mut over70) = (0.0, 0.0, 0, 0, 0, 0);
        let mut tr_med = Vec::new();
        for seed in seeds15() {
            let Ok(s) = std::fs::read_to_string(results178().join(format!("t1_{name}_seed{seed}.txt"))) else { continue };
            if name == "F178" || name == "F176" {
                eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
            }
            let a = data_line(&s, "A1").unwrap();
            let u = data_line(&s, "UP").unwrap();
            let t = data_line(&s, "TR").unwrap();
            let l = data_line(&s, "LAND").unwrap();
            any.push(100.0 * a[0]);
            maj.push(100.0 * a[1]);
            conts += a[2];
            conts_a += a[3];
            if a[0] > 0.7 {
                over70 += 1;
            }
            ratio.push(u[0] / u[1]);
            wedge.push(u[2]);
            if u[2] > 0.0 {
                wpos += 1;
            }
            if u[0] >= 2.0 * u[1] {
                r2 += 1;
            }
            if l[1] == 0.0 {
                circ_lost += 1;
            }
            tr_med.push(t[1]);
        }
        eprintln!(
            "\n   {name}: A1 median any {:.1} % · majority {:.1} % (sorted majority {:?}) · continents with an active margin {conts_a} / {conts} · > 70 % active {over70} · C-O wedge on the active 0–6 cells positive on {wpos} / {} (median {:.0} m) · active ≥ 2 × passive on 0–6 cells {r2} · ratio median {:.2} · circumnavigation lost {circ_lost} · trench segment medians (km) {:?}\n",
            med(any.clone()), med(maj.clone()), { let mut m2 = maj.clone(); m2.sort_by(|a, b| a.partial_cmp(b).unwrap()); m2.iter().map(|x| (x * 10.0).round() / 10.0).collect::<Vec<_>>() }, any.len(), med(wedge.clone()), med(ratio.clone()), tr_med
        );
    }
    // T7: A2 against F178, the active coasts' uplift on 0–6 cells
    let mut ch = Vec::new();
    for seed in seeds15() {
        let g = |n: &str| std::fs::read_to_string(results178().join(format!("t1_{n}_seed{seed}.txt"))).ok().and_then(|s| data_line(&s, "UP"));
        if let (Some(a), Some(b)) = (g("F178"), g("F178_A2")) {
            ch.push((seed, 100.0 * (b[0] - a[0]) / a[0].abs().max(1e-9)));
        }
    }
    eprintln!("   T7: A2's change of the active coasts' 0–6 cell uplift (%) {:?}", ch.iter().map(|(s, c)| (*s, (c * 10.0).round() / 10.0)).collect::<Vec<_>>());
    let tiles: Vec<image::RgbaImage> = PANEL.iter().filter_map(|s| image::open(out178().join(format!("a1_F178_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 388, 388).save(out178().join("a1_panel_F178.png")).unwrap();
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// The formations of one world at 256² young (F177's complement definitions), its W₁₀ and macro: (the line, the formation map).
fn formations178(o: &Obs, a: &A1, r: &Run, w: &World) -> (String, Vec<u8>) {
    let km2 = CELL * CELL;
    let n = 256;
    let org = origin64(o);
    let org256 = |k: usize| org[((k / n) / 4) * 64 + (k % n) / 4];
    let mtn: Vec<bool> = (0..n * n).map(|k| r.classes[k] == MOUNTAIN && r.zm.data[k] > 0.0).collect();
    let land: Vec<bool> = r.zm.data.iter().map(|&z| z > 0.0).collect();
    let lc = comps8(&land, n);
    let mut lab = vec![usize::MAX; n * n];
    for (i, c) in lc.iter().enumerate() {
        for &k in c {
            lab[k] = i;
        }
    }
    let main = (0..lc.len()).max_by_key(|&i| lc[i].len()).unwrap_or(0);
    let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let ch = channels(&land64);
    let org4 = origin_b4(o);
    let mut chan_cls = vec![0u8; N64 * N64];
    for c in &ch.chans {
        let mj = majority(c, &org4);
        for &k in c {
            chan_cls[k] = mj;
        }
    }
    let mut form = vec![0u8; n * n];
    let and_mask: Vec<bool> = (0..n * n).map(|k| mtn[k] && org256(k) == 1).collect();
    let mut cords = Vec::new();
    for g in comps8(&and_mask, n) {
        let mut touched = std::collections::BTreeSet::new();
        for &k in &g {
            let mut cand = vec![k];
            cand.extend(nb8_n(k, n));
            for c in cand {
                let c64 = ((c / n) / 4) * 64 + (c % n) / 4;
                if is_active_cls(a.cls[c64]) {
                    touched.insert(c64);
                }
            }
        }
        let f = if touched.is_empty() { 2 } else { 1 };
        for &k in &g {
            form[k] = f;
        }
        if f == 1 {
            cords.push((g.len(), touched.len()));
        }
    }
    let mut sep_classes = Vec::new();
    let col_mask: Vec<bool> = (0..n * n).map(|k| mtn[k] && org256(k) == 2).collect();
    for g in comps8(&col_mask, n) {
        let sep = lab[g[0]] != main;
        for &k in &g {
            form[k] = if sep { 4 } else { 3 };
        }
        if sep {
            let mut best = (i64::MAX, 0u8);
            for &k in &g {
                let (x, y) = (((k % n) / 4) as i64, ((k / n) / 4) as i64);
                for q in 0..N64 * N64 {
                    if chan_cls[q] == 0 {
                        continue;
                    }
                    let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
                    dx = dx.min(64 - dx);
                    let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
                    dy = dy.min(64 - dy);
                    if dx * dx + dy * dy < best.0 {
                        best = (dx * dx + dy * dy, chan_cls[q]);
                    }
                }
            }
            sep_classes.push((g.len(), B4_NAMES[best.1 as usize], ((best.0 as f32).sqrt() * CELL64_KM) as i32));
        }
    }
    for k in 0..n * n {
        if mtn[k] && form[k] == 0 {
            form[k] = if org256(k) == 3 { 5 } else { 6 };
        }
    }
    let area = |f: u8| form.iter().filter(|&&v| v == f).count() as f32 * km2;
    let (cord_km2, and_other, col, sep, arc, oth) = (area(1), area(2), area(3), area(4), area(5), area(6));
    let cord_len: f32 = cords.iter().map(|c| c.1 as f32 * CELL64_KM).sum();
    let only_collision = cord_km2 + and_other + arc + oth == 0.0 && col + sep > 0.0;
    let map_end: Vec<u8> = a.cls.iter().map(|&c| if is_active_cls(c) { 3 } else if is_passive_cls(c) { 2 } else { 0 }).collect();
    let (wi, win, rw) = wall_index(&r.zm);
    let (wa, wp, na, np) = w10_by_type(&r.zm, &map_end);
    let mm = &r.m;
    let line = format!(
        "      formations: Andean cordilleras {} (lengths km {:?}; total {cord_len:.0} km world, {:.0} km terrestrial; {cord_km2:.0} km²) · Andean not touching an active coast {and_other:.0} km² · collisional highland {col:.0} km² · separated {sep:.0} km² (groups: cells, nearest channel, km {:?}) · arc mountain {arc:.0} km² · other {oth:.0} km² · collision the only formation: {only_collision}\n      W10 {wi:.3} (W_in {win:.3}, rampart {rw}) · active / passive {wa:.3} / {wp:.3} ({na} / {np} cells) · macro plain {:.1} · plateau {:.1} · hill {:.1} · mountain {:.1} % · land {:.0} km² · median {:.0} m\n#F {} {cord_len} {cord_km2} {and_other} {col} {sep} {arc} {oth} {}\n#W {wi} {wa} {wp} {}\n#M {} {} {} {} {} {}\n",
        cords.len(), cords.iter().map(|c| c.1 as f32 * CELL64_KM).collect::<Vec<_>>(), cord_len * 7.5, sep_classes,
        100.0 * mm.fractions[0], 100.0 * mm.fractions[1], 100.0 * mm.fractions[2], 100.0 * mm.fractions[3], mm.land_km2, mm.deciles.get(4).copied().unwrap_or(f32::NAN),
        cords.len(), only_collision as u8, rw as u8, mm.fractions[0], mm.fractions[1], mm.fractions[2], mm.fractions[3], mm.land_km2, mm.deciles.get(4).copied().unwrap_or(f32::NAN)
    );
    (line, form)
}

/// F178 tier 2 -- the physics level at 256², young, the 15 seeds: F178 and F176 (B3, the formations, W₁₀, the macro).
#[allow(dead_code)]
fn f178_tier2_physics() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out178()).unwrap();
    eprintln!("\n==========  F178 tier 2 . 256² young: B3, the formations, W10, the macro  ==========");
    let rough_a = rough_a177();
    for (name, m) in variants178().into_iter().take(2) {
        for seed in seeds15() {
            let Some(p) = sel176(seed) else { continue };
            let s = cached178(&format!("t2_{name}_seed{seed}"), || {
                let o = observe(seed, p, &m);
                let info = step_info(&o.st, &o.kin, o.ds_max, true);
                let a = a1(&info, &o.kin, &o.tracker);
                let (r, w) = physics_young(seed, p, &m, rough_a);
                let (line, form) = formations178(&o, &a, &r, &w);
                let org = origin64(&o);
                origin256_png(&r.zm, &r.classes, &org, None).save(out178().join(format!("origin_{name}_seed{seed}.png"))).unwrap();
                origin256_png(&r.zm, &r.classes, &org, Some(&form)).save(out178().join(format!("formations_{name}_seed{seed}.png"))).unwrap();
                format!("   seed {seed} {name} young ({p:?}): physics CPU {:.0} s\n{line}", r.cost.1)
            });
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        }
        let (mut with_cord, mut only, mut n) = (0, 0, 0);
        let mut w10_panel = Vec::new();
        let mut mtn = Vec::new();
        for seed in seeds15() {
            if let Ok(s) = std::fs::read_to_string(results178().join(format!("t2_{name}_seed{seed}.txt"))) {
                let f = data_line(&s, "F").unwrap();
                let wv = data_line(&s, "W").unwrap();
                let mv = data_line(&s, "M").unwrap();
                n += 1;
                with_cord += (f[0] > 0.0) as usize;
                only += (f[8] > 0.0) as usize;
                if PANEL.contains(&seed) {
                    w10_panel.push((seed, (wv[0] * 1000.0).round() / 1000.0));
                }
                mtn.push((seed, (mv[3] * 1000.0).round() / 10.0));
            }
        }
        eprintln!("\n   {name}: seeds with ≥ 1 Andean cordillera {with_cord} / {n} · collision the only formation on {only} · W10 young on the panel {w10_panel:?} · mountain % {mtn:?}\n");
        let tiles: Vec<image::RgbaImage> = seeds15().iter().filter_map(|s| image::open(out178().join(format!("formations_{name}_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
        montage(&tiles, 5, 260, 260).save(out178().join(format!("formations_15_{name}.png"))).unwrap();
        let tiles: Vec<image::RgbaImage> = seeds15().iter().filter_map(|s| image::open(out178().join(format!("origin_{name}_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
        montage(&tiles, 5, 260, 260).save(out178().join(format!("origin_15_{name}.png"))).unwrap();
    }
    eprintln!("   tier 2: {:.0} s", t0.elapsed().as_secs_f64());
}

/// One unclipped 2 048² transect pair (an active and a passive coast cell of the A1 map) for a world.
fn transect178(seed: u64, p: C1Profile, m: &Margin, active: Option<usize>, passive: Option<usize>, path: &std::path::Path) -> String {
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let kh = KHetero::declared();
    let w = world_oo(seed, p, m);
    let cover = cover_map_m(seed, p, m, COVER_M);
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k.clone()) };
    let cfg = base_cfg(seed, rough_a, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[1].1, &rock, N_PHYS);
    let (g, _, _) = chain_low(&r.zm, amp, Some(ln_k), 0.0);
    let s64 = w.st.s.data().to_vec();
    let pt = w.st.plate_type.data().to_vec();
    let (mut zs, mut ss, mut desc) = (Vec::new(), Vec::new(), Vec::new());
    for (nm, cell, col) in [("passive", passive, [30u8, 140, 60]), ("active", active, [220, 30, 30])] {
        if let Some(k) = cell {
            let (a, b) = transect(&g, &s64, &pt, k);
            let mx = a.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            let at_max = a.iter().filter(|&&v| (v - mx).abs() < 1.0).count();
            desc.push(format!("{nm} (cell {k}): max {mx:.0} m, samples within 1 m of the max {at_max} of {}", a.len()));
            zs.push((a, col));
            ss.push((b, col));
        }
    }
    let (zmin, zmax) = transect_png_free(&zs, &ss, path);
    format!("   the range {zmin:.0} … {zmax:.0} m · {}\n", desc.join(" · "))
}

/// F178 tier 3 -- the unclipped transects: F178 seed 18, the F178 seed with the longest new cordillera, F176's seed 18 (fixed).
#[allow(dead_code)]
fn f178_tier3_transects() {
    let t0 = Instant::now();
    eprintln!("\n==========  F178 tier 3 . the transects without a clamp (2 048², young)  ==========");
    let m178 = Margin::seaward_oo_collision_andean();
    // the active cell: the A1-active coast cell with the largest C-O wedge; the passive: the passive coast cell with the least
    let cells = |seed: u64, p: C1Profile, m: &Margin| -> (Option<usize>, Option<usize>) {
        let o = observe(seed, p, m);
        let info = step_info(&o.st, &o.kin, o.ds_max, true);
        let a = a1(&info, &o.kin, &o.tracker);
        let act = (0..4096).filter(|&k| is_active_cls(a.cls[k])).max_by(|&x, &y| o.w_sum[x].partial_cmp(&o.w_sum[y]).unwrap().then(y.cmp(&x)));
        let pas = (0..4096).filter(|&k| is_passive_cls(a.cls[k])).min_by(|&x, &y| o.src_sum[x].partial_cmp(&o.src_sum[y]).unwrap().then(x.cmp(&y)));
        (act, pas)
    };
    let p18 = sel176(18).unwrap();
    let s = cached178("t3_F178_seed18", || {
        let (a, pz) = cells(18, p18, &m178);
        format!("   seed 18 F178 (selection) young:\n{}", transect178(18, p18, &m178, a, pz, &out178().join("transect_F178_seed18.png")))
    });
    eprint!("{s}");
    // the seed with the longest cordillera under F178 among those without one under F176
    let cord = |name: &str, seed: u64| -> f64 { std::fs::read_to_string(results178().join(format!("t2_{name}_seed{seed}.txt"))).ok().and_then(|s| data_line(&s, "F")).map(|f| f[1]).unwrap_or(0.0) };
    let pick = seeds15().into_iter().filter(|&sd| sd != 18).max_by(|&x, &y| {
        let k = |sd: u64| (cord("F176", sd) == 0.0, cord("F178", sd));
        let (kx, ky) = (k(x), k(y));
        kx.0.cmp(&ky.0).then(kx.1.partial_cmp(&ky.1).unwrap())
    });
    if let Some(sd) = pick {
        let p = sel176(sd).unwrap();
        let s = cached178(&format!("t3_F178_seed{sd}"), || {
            let (a, pz) = cells(sd, p, &m178);
            format!("   seed {sd} F178 (selection, the longest new cordillera) young:\n{}", transect178(sd, p, &m178, a, pz, &out178().join(format!("transect_F178_seed{sd}.png"))))
        });
        eprint!("{s}");
    }
    // F176's image unclipped: seed 18 at the fixed continent, F176's cells (its 3 × 3 map)
    let m176 = Margin::seaward_oo_collision();
    let s = cached178("t3_F176_seed18_fixed", || {
        let ser = c1_series(18, fixed_profile(18), &m176);
        let act = (0..4096).find(|&k| ser.map_end[k] == 3 && ser.map_t0[k] == 2).or_else(|| (0..4096).find(|&k| ser.map_end[k] == 3));
        let pas = (0..4096).find(|&k| ser.map_end[k] == 2);
        format!("   seed 18 F176 (fixed continent, F176's transect cells) young:\n{}", transect178(18, fixed_profile(18), &m176, act, pas, &out178().join("transect_F176_seed18_fixed_noclip.png")))
    });
    eprint!("{s}");
    eprintln!("   tier 3: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F178's stop, its mechanism (C1 only; measures, no closure changed): on the active 0–3 cell band, the per-term Σ Δs over the run
/// (C1's own record, `C1Sources::terms`), and the land F178 adds over F176 (the type, the seaward band, the t₀ thickness).
#[allow(dead_code)]
fn f178_stop_mechanism() {
    let seeds = [23u64, 3, 4, 10, 18, 8];
    std::thread::scope(|sc| {
        for seed in seeds {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                cached178(&format!("stopmech_seed{seed}"), || {
                    let nn = N64 * N64;
                    let mut out = String::new();
                    let mut land_of = Vec::new();
                    for (name, m) in variants178().into_iter().take(2) {
                        let o = observe(seed, p, &m);
                        let w = world_oo(seed, p, &m);
                        let info = step_info(&o.st, &o.kin, o.ds_max, true);
                        let a = a1(&info, &o.kin, &o.tracker);
                        let pt = o.st.plate_type.data();
                        let (comp, big) = continents64(pt);
                        let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
                        // the active 0–3 band (nearest coast active, walk depth < 3)
                        let mut lab = vec![0u8; nn];
                        let mut dd = vec![u32::MAX; nn];
                        let mut q = std::collections::VecDeque::new();
                        for k in 0..nn {
                            if in_cont(k) && (is_active_cls(a.cls[k]) || is_passive_cls(a.cls[k])) {
                                lab[k] = if is_active_cls(a.cls[k]) { 1 } else { 2 };
                                dd[k] = 0;
                                q.push_back(k);
                            }
                        }
                        while let Some(k) = q.pop_front() {
                            for j in nb4_64(k) {
                                if pt[j] == PlateType::Continental && dd[j] == u32::MAX && dd[k] + 1 < 3 {
                                    dd[j] = dd[k] + 1;
                                    lab[j] = lab[k];
                                    q.push_back(j);
                                }
                            }
                        }
                        let band: Vec<usize> = (0..nn).filter(|&k| lab[k] == 1).collect();
                        let nb = band.len().max(1) as f64;
                        let term = |t: C1Term| -> f64 { w.src.terms.iter().find(|(x, _)| *x == t).map(|(_, v)| band.iter().map(|&k| v[k]).sum::<f64>() / nb).unwrap_or(f64::NAN) };
                        let s_end: f64 = band.iter().map(|&k| o.st.s.data()[k]).sum::<f64>() / nb;
                        let s0: f64 = band.iter().map(|&k| o.s0[k]).sum::<f64>() / nb;
                        out += &format!(
                            "   seed {seed} {name}: the active 0–3 band {} cells · mean Σ Δs (s units; 1 = 35 km of crust): Davis-Suppe {:+.2} · equilibrium height {:+.2} · erosion {:+.2} · transport {:+.2} · subduction {:+.2} · rift {:+.2} · s at t₀ {s0:.2} → C1's s at the end {s_end:.2} · the cascade's s_src (s minus the erosion) {:.2}\n",
                            band.len(), term(C1Term::DavisSuppe), term(C1Term::EquilibriumHeight), term(C1Term::Erosion), term(C1Term::Transport), term(C1Term::Subduction), term(C1Term::RiftThinning), s_end - term(C1Term::Erosion)
                        );
                        land_of.push((w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect::<Vec<bool>>(), o));
                    }
                    let (l176, _) = &land_of[0];
                    let (l178, o178) = &land_of[1];
                    let new: Vec<usize> = (0..nn).filter(|&k| l178[k] && !l176[k]).collect();
                    let lost = (0..nn).filter(|&k| !l178[k] && l176[k]).count();
                    let pt = o178.st.plate_type.data();
                    let n = new.len().max(1) as f32;
                    out += &format!(
                        "   seed {seed}: land F178 − F176: +{} cells, −{lost} · of the new land: continental type {:.0} %, in the seaward band at t₀ {:.0} %, mean s at t₀ {:.2} (sea level 0.915), with a C-O wedge increment {:.0} %\n",
                        new.len(),
                        100.0 * new.iter().filter(|&&k| pt[k] == PlateType::Continental).count() as f32 / n,
                        100.0 * new.iter().filter(|&&k| o178.band[k]).count() as f32 / n,
                        new.iter().map(|&k| o178.s0[k] as f32).sum::<f32>() / n,
                        100.0 * new.iter().filter(|&&k| o178.w_sum[k] > 0.0).count() as f32 / n
                    );
                    out
                });
            });
        }
    });
    for seed in seeds {
        if let Ok(s) = std::fs::read_to_string(results178().join(format!("stopmech_seed{seed}.txt"))) {
            eprint!("{s}");
        }
    }
}

// ════════════════════════════════════════ F179 ════════════════════════════════════════

use ymir_core::tectonics_c1::andean_seed::{andean_seeds_active, andean_seeds_arc, overriding_trench_partners};
use ymir_core::cascade::margin::D_ARC_CELLS;

fn out179() -> PathBuf {
    root().join("docs/reports/relief_method/f179_pumping/images")
}
fn results179() -> PathBuf {
    root().join("docs/reports/relief_method/f179_pumping/results")
}
fn cached179(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results179().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results179()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
fn profiles179() -> Vec<(&'static str, Margin)> {
    vec![
        ("F175", Margin::seaward_oo()),
        ("F176", Margin::seaward_oo_collision()),
        ("F178", Margin::seaward_oo_collision_andean()),
        ("F179", Margin::seaward_oo_collision_arc()),
    ]
}
/// Davis-Suppe's wedge width L (cells of 64²): its reach, `max_distance` (f179_declared, A2).
fn wedge_l() -> f64 {
    C1Closures::default().davis_suppe.scaled_to_grid(64).max_distance
}
const FAM_NAMES: [&str; 5] = ["none", "collision", "O-C boundary", "Andean", "other"];

/// The pumping record of one run.
struct Pump {
    /// Σ Davis-Suppe per family [1..=4] per cell.
    ds_fam: Vec<Vec<f64>>,
    steps: [usize; 5],
    exceed: [usize; 5],
    /// Σ max(0, inc − Δs_kin) per cell, and per family.
    excess: Vec<f64>,
    excess_fam: [f64; 5],
    /// The Andean seed instances over the run: their walk depth from the nearest trench (0..=10, 11 = beyond) and emerged at t₀.
    seed_d: [usize; 12],
    seed_emerged: usize,
    seed_total: usize,
    /// The cells seeded at least once.
    seeded: Vec<bool>,
    /// F180: the bound's record by family [1..=4]: seed-steps, bounded seed-steps, mass removed (s × cells).
    bnd_seeds: [usize; 5],
    bnd_bounded: [usize; 5],
    removed_fam: [f64; 5],
    /// F180 diagnostic: option 5's consumptions decided by age, of which |Δage| < 0.01, and by the tie rule (over the run).
    oo_age: usize,
    oo_age_close: usize,
    oo_tie: usize,
}

/// One step's start state for F179: the loop's seeds and wedge (as `step_info_a`), the families and the supply per wedge cell.
struct Step179 {
    info: StepInfo,
    fam: Vec<u8>,
    supply: Vec<f64>,
    andean: Vec<bool>,
    seed_depth: Vec<u32>,
}
fn step179(s: &C1State, kn: &PlateKinematics, ds_max: f64, m: &Margin, tracker: &OoTracker) -> Step179 {
    let nn = N64 * N64;
    let pid = s.plate_id.data();
    let pt = s.plate_type.data();
    let mut bi = classify_boundaries(&s.plate_id, kn);
    retarget_upper_plate_continental(&mut bi, &s.plate_id, &s.plate_type);
    let upper_vel = bi.upper_plate_mask.data().to_vec();
    let uv = bi.upper_plate_mask.clone();
    if m.oo.is_some() {
        clear_oo_seeds(&mut bi, &s.plate_id, &s.plate_type);
    }
    let oc_bnd = oc_override_seed_mask(&bi, &s.plate_id, &s.plate_type);
    let andean_bf = match (m.andean_active, m.andean_arc, m.andean_gap) {
        (Some((d_arc, gap)), _, _) => Some(andean_seeds_active(&s.plate_id, &s.plate_type, &s.age, &bi, &uv, kn, tracker, d_arc, gap).mask),
        (None, Some((d_arc, gap)), _) => Some(andean_seeds_arc(&s.plate_id, &s.plate_type, &s.age, &bi, &uv, kn, tracker, d_arc, gap)),
        (None, None, Some(gap)) => Some(andean_seeds(&s.plate_id, &s.plate_type, &s.age, &bi, &uv, kn, tracker, gap)),
        _ => None,
    };
    let andean: Vec<bool> = andean_bf.as_ref().map(|b| b.data().to_vec()).unwrap_or_else(|| vec![false; nn]);
    // the families of the seeds (before the Andean ones join the masks)
    let up0 = bi.upper_plate_mask.data().to_vec();
    let bt = bi.boundary_type.data().to_vec();
    let mut seed_fam = vec![0u8; nn];
    let mut seed_sup = vec![0f64; nn];
    let edge_supply = |k: usize| -> f64 {
        let (x, y) = (k % N64, k / N64);
        let (vxc, vyc) = kn.velocities[pid[k] as usize];
        let mut best = 0f64;
        for (q, nx_, ny_) in [(y * N64 + (x + 1) % N64, 1.0, 0.0), (y * N64 + (x + N64 - 1) % N64, -1.0, 0.0), (((y + 1) % N64) * N64 + x, 0.0, 1.0), (((y + N64 - 1) % N64) * N64 + x, 0.0, -1.0)] {
            if pid[q] == pid[k] {
                continue;
            }
            let (vxn, vyn) = kn.velocities[pid[q] as usize];
            let dot = (vxc - vxn) * nx_ + (vyc - vyn) * ny_;
            if dot > 0.0 {
                best = best.max(dot * s.s.data()[q]);
            }
        }
        best
    };
    for k in 0..nn {
        if !up0[k] {
            continue;
        }
        let cc = pt[k] == PlateType::Continental && bt[k] == BoundaryType::Convergent && nb4_64(k).iter().any(|&j| pid[j] != pid[k] && pt[j] == PlateType::Continental);
        seed_fam[k] = if oc_bnd.data()[k] { 2 } else if cc { 1 } else { 4 };
        seed_sup[k] = edge_supply(k);
    }
    // the Andean seeds' supply: the trench's (the partner's s × v_rel · n̂), carried by a walk through P's cells; their depth
    let mut seed_depth = vec![u32::MAX; nn];
    if andean.iter().any(|&a| a) {
        let tr = overriding_trench_cells(&s.plate_id, &s.plate_type, &s.age, &bi, &uv, kn, tracker);
        let mut sup = vec![0f64; nn];
        let mut q = std::collections::VecDeque::new();
        for k in 0..nn {
            if tr.data()[k] {
                seed_depth[k] = 0;
                sup[k] = edge_supply(k);
                q.push_back(k);
            }
        }
        while let Some(k) = q.pop_front() {
            for j in nb4_64(k) {
                if pid[j] == pid[k] && seed_depth[j] == u32::MAX && seed_depth[k] < 11 {
                    seed_depth[j] = seed_depth[k] + 1;
                    sup[j] = sup[k];
                    q.push_back(j);
                }
            }
        }
        for k in 0..nn {
            if andean[k] {
                seed_fam[k] = 3;
                seed_sup[k] = sup[k];
                bi.upper_plate_mask.set(k % N64, k / N64, true);
            }
        }
    }
    let mut oc_seed = oc_bnd.clone();
    for k in 0..nn {
        if andean[k] {
            oc_seed.set(k % N64, k / N64, true);
        }
    }
    let (_, ocw) = wedge_distance_intra_plate_typed(&s.plate_id, &bi.upper_plate_mask, &oc_seed, ds_max);
    // the nearest seed's family and supply (the wedge's metric: 8-connected, √2 diagonals, within the plate)
    let mut dist = vec![f64::INFINITY; nn];
    let mut fam = vec![0u8; nn];
    let mut supply = vec![0f64; nn];
    let mut heap = std::collections::BinaryHeap::new();
    for k in 0..nn {
        if bi.upper_plate_mask.data()[k] {
            dist[k] = 0.0;
            fam[k] = seed_fam[k].max(1);
            supply[k] = seed_sup[k];
            heap.push(std::cmp::Reverse(((0.0f64 * 1e6) as i64, k)));
        }
    }
    while let Some(std::cmp::Reverse((dk, k))) = heap.pop() {
        let d = dk as f64 / 1e6;
        if d > dist[k] + 1e-9 || d >= ds_max {
            continue;
        }
        let (x, y) = (k % N64, k / N64);
        for (dx, dy, c) in [(1i64, 0i64, 1.0), (-1, 0, 1.0), (0, 1, 1.0), (0, -1, 1.0), (1, 1, std::f64::consts::SQRT_2), (1, -1, std::f64::consts::SQRT_2), (-1, 1, std::f64::consts::SQRT_2), (-1, -1, std::f64::consts::SQRT_2)] {
            let j = (((y as i64 + dy).rem_euclid(64)) as usize) * N64 + ((x as i64 + dx).rem_euclid(64)) as usize;
            if pid[j] != pid[k] {
                continue;
            }
            let nd = d + c;
            if nd < dist[j] - 1e-9 {
                dist[j] = nd;
                fam[j] = fam[k];
                supply[j] = supply[k];
                heap.push(std::cmp::Reverse(((nd * 1e6) as i64, j)));
            }
        }
    }
    let plate_id: Vec<u16> = pid.to_vec();
    let plate_type: Vec<PlateType> = pt.to_vec();
    let co_cell = (0..nn).map(|k| plate_type[k] == PlateType::Oceanic && nb4_64(k).iter().any(|&j| plate_id[j] != plate_id[k] && plate_type[j] == PlateType::Continental)).collect();
    let info = StepInfo { bi, upper_vel, oc_wedge: ocw.data().to_vec(), co_cell, plate_id, plate_type, age: s.age.data().to_vec() };
    Step179 { info, fam, supply, andean, seed_depth }
}

/// F179's observer: F177/F178's sums (the `Obs`) plus the pumping record.
fn observe179(seed: u64, p: C1Profile, m: &Margin) -> (Obs, Pump) {
    let nn = N64 * N64;
    let cfg = c1_run_config();
    let mut st = p.init_state_margin(64, seed, m);
    let st0 = p.init_state_margin(64, seed, m);
    let st_off = p.init_state_margin(64, seed, &Margin::off());
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let kin0 = PlateKinematics::preset_phase_1_1(st.num_plates);
    let dt = 0.5 * cfg.dx.min(cfg.dy) / kin.max_velocity().max(1e-12);
    let ds_max = C1Closures::default().davis_suppe.scaled_to_grid(64).max_distance;
    let l = wedge_l();
    let band: Vec<bool> = (0..nn).map(|k| st.plate_id.data()[k] != st_off.plate_id.data()[k] || st.plate_type.data()[k] != st_off.plate_type.data()[k]).collect();
    let s0 = st.s.data().to_vec();
    let emerged0: Vec<bool> = (0..nn).map(|k| st.plate_type.data()[k] == PlateType::Continental && s0[k] >= 0.915).collect();
    let mut tracker = OoTracker::default();
    let first = step179(&st, &kin, ds_max, m, &tracker);
    if let Some(o) = m.oo {
        tracker.update(&st.plate_id, &st.plate_type, &first.info.bi, &kin, dt, o.d_init_domain(), 0);
    }
    let cur = std::cell::RefCell::new(first);
    let tr = std::cell::RefCell::new(tracker);
    let sums = std::cell::RefCell::new([(); 8].map(|_| vec![0f64; nn]));
    let pump = std::cell::RefCell::new(Pump {
        ds_fam: vec![vec![0f64; nn]; 5],
        steps: [0; 5],
        exceed: [0; 5],
        excess: vec![0f64; nn],
        excess_fam: [0.0; 5],
        seed_d: [0; 12],
        seed_emerged: 0,
        seed_total: 0,
        seeded: vec![false; nn],
        bnd_seeds: [0; 5],
        bnd_bounded: [0; 5],
        removed_fam: [0.0; 5],
        oo_age: 0,
        oo_age_close: 0,
        oo_tie: 0,
    });
    let record_seeds = |c: &Step179, pm: &mut Pump| {
        for k in 0..nn {
            if c.andean[k] {
                pm.seed_total += 1;
                pm.seed_d[(c.seed_depth[k] as usize).min(11)] += 1;
                pm.seed_emerged += emerged0[k] as usize;
                pm.seeded[k] = true;
            }
        }
    };
    record_seeds(&cur.borrow(), &mut pump.borrow_mut());
    // sums: 0 a (C-O wedge + C-O arcs), 1 c (non-O-C Davis-Suppe), 2 o (O-O arcs), 3 src (every term but the erosion),
    // 4 rift, 5 w (Davis-Suppe on oc_wedge), 6 arcs (every positive subduction increment)
    let mut on_term = |t: C1Term, a: &[f64], b: &[f64]| {
        let mut sm = sums.borrow_mut();
        if t != C1Term::Erosion {
            for k in 0..nn {
                sm[3][k] += b[k] - a[k];
            }
        }
        match t {
            C1Term::DavisSuppe => {
                let c = cur.borrow();
                let mut pm = pump.borrow_mut();
                for k in 0..nn {
                    let d = b[k] - a[k];
                    if d <= 0.0 {
                        continue;
                    }
                    if c.info.oc_wedge[k] {
                        sm[0][k] += d;
                        sm[5][k] += d;
                    } else {
                        sm[1][k] += d;
                    }
                    let f = if c.fam[k] == 0 { 4 } else { c.fam[k] as usize };
                    pm.ds_fam[f][k] += d;
                    pm.steps[f] += 1;
                    let kin_s = c.supply[k] * dt * 64.0 / l;
                    if d > kin_s {
                        pm.exceed[f] += 1;
                        pm.excess[k] += d - kin_s;
                        pm.excess_fam[f] += d - kin_s;
                    }
                }
            }
            C1Term::Subduction => {
                let c = cur.borrow();
                let consumed: Vec<usize> = (0..nn).filter(|&k| b[k] < a[k]).collect();
                for k in 0..nn {
                    let d = b[k] - a[k];
                    if d > 0.0 {
                        sm[6][k] += d;
                    }
                    if d <= 0.0 || consumed.is_empty() {
                        continue;
                    }
                    let (x, y) = ((k % N64) as i64, (k / N64) as i64);
                    let mut best = (i64::MAX, false);
                    for &q in &consumed {
                        let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
                        dx = dx.min(64 - dx);
                        let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
                        dy = dy.min(64 - dy);
                        let d2 = dx * dx + dy * dy;
                        if d2 < best.0 || (d2 == best.0 && c.info.co_cell[q] && !best.1) {
                            best = (d2, c.info.co_cell[q]);
                        }
                    }
                    if best.1 {
                        sm[0][k] += d;
                    } else {
                        sm[2][k] += d;
                    }
                }
            }
            C1Term::RiftThinning => {
                for k in 0..nn {
                    sm[4][k] += b[k] - a[k];
                }
            }
            _ => {}
        }
    };
    let merges = std::cell::Cell::new(0usize);
    let splits = std::cell::Cell::new(0usize);
    run_with_closures_terms_ext(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        m.ext(),
        |step, s, kn| {
            merges.set(merges.get() + s.last_step_stats.accretion.merges_count);
            splits.set(splits.get() + s.last_step_stats.rifting_split.splits_count);
            {
                // F180: the bound's record of this step, by the step's families
                let c = cur.borrow();
                let mut pmb = pump.borrow_mut();
                let sub = &s.last_step_stats.subduction;
                pmb.oo_age += sub.oo_age_decided;
                pmb.oo_age_close += sub.oo_age_close;
                pmb.oo_tie += sub.oo_tie_decided;
                let bs = &s.last_step_stats.ds_bound;
                for &(sd, bd) in &bs.seed_bounded {
                    let f = if c.fam[sd] == 0 { 4 } else { c.fam[sd] as usize };
                    pmb.bnd_seeds[f] += 1;
                    pmb.bnd_bounded[f] += bd as usize;
                }
                for (k, &r) in bs.removed.iter().enumerate() {
                    if r > 0.0 {
                        let f = if c.fam[k] == 0 { 4 } else { c.fam[k] as usize };
                        pmb.removed_fam[f] += r;
                    }
                }
            }
            let trb = tr.borrow();
            let ni = step179(s, kn, ds_max, m, &trb);
            drop(trb);
            if let Some(o) = m.oo {
                if step + 1 < cfg.n_steps {
                    tr.borrow_mut().update(&s.plate_id, &s.plate_type, &ni.info.bi, kn, dt, o.d_init_domain(), step + 1);
                }
            }
            if step + 1 < cfg.n_steps {
                record_seeds(&ni, &mut pump.borrow_mut());
            }
            *cur.borrow_mut() = ni;
        },
        &mut on_term,
    );
    let sm = sums.into_inner();
    let [a_sum, c_sum, o_sum, src_sum, rift_sum, w_sum, arc_sum, _] = sm;
    let o = Obs {
        a_sum,
        c_sum,
        o_sum,
        src_sum,
        rift_sum,
        merged_ocean: vec![false; nn],
        strip: vec![false; nn],
        merges: merges.get(),
        splits: splits.get(),
        s0,
        band,
        tracker: tr.into_inner(),
        st,
        kin,
        st0,
        kin0,
        ds_max,
        w_sum,
        arc_sum,
    };
    (o, pump.into_inner())
}

fn q_of(v: &mut Vec<f64>, f: f64) -> f64 {
    v.retain(|x| !x.is_nan());
    v.sort_by(|a, b| a.total_cmp(b));
    if v.is_empty() { f64::NAN } else { v[((v.len() - 1) as f64 * f) as usize] }
}

/// One F179 C1 analysis of a profile: A1–A3 (the pumping), B (the land, r*, the added land, the bands, A1, the seeds).
fn c1_179(seed: u64, p: C1Profile, name: &str, m: &Margin) -> String {
    let clk = Clock::start();
    let (o, pm) = observe179(seed, p, m);
    let w = world_oo(seed, p, m);
    let cost = clk.stop();
    let nn = N64 * N64;
    let info = step_info(&o.st, &o.kin, o.ds_max, true);
    let a = a1(&info, &o.kin, &o.tracker);
    let pt = o.st.plate_type.data();
    let (comp, big) = continents64(pt);
    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
    // the dominant family per cell
    let dom: Vec<u8> = (0..nn)
        .map(|k| {
            let (f, v) = (1..5).map(|f| (f, pm.ds_fam[f][k])).fold((0usize, 0f64), |acc, x| if x.1 > acc.1 { x } else { acc });
            if v > 0.0 { f as u8 } else { 0 }
        })
        .collect();
    let s_end = o.st.s.data();
    let mut s = format!("   seed {seed} {name}: C1 CPU {:.1} s · merges {} splits {}\n", cost.1, o.merges, o.splits);
    // A1: the pumping factor per family
    for f in 1..5u8 {
        let cells: Vec<usize> = (0..nn).filter(|&k| dom[k] == f).collect();
        let tot: f64 = cells.iter().map(|&k| (1..5).map(|g| pm.ds_fam[g][k]).sum::<f64>()).sum();
        let net_tot: f64 = cells.iter().map(|&k| s_end[k] - o.s0[k]).sum();
        let mut fac: Vec<f64> = cells.iter().filter(|&&k| s_end[k] - o.s0[k] > 1e-9).map(|&k| (1..5).map(|g| pm.ds_fam[g][k]).sum::<f64>() / (s_end[k] - o.s0[k])).collect();
        let nonpos = cells.len() - fac.len();
        let (md, p9, mx) = (q_of(&mut fac.clone(), 0.5), q_of(&mut fac.clone(), 0.9), q_of(&mut fac, 1.0));
        let exc = cells.iter().map(|&k| pm.excess[k]).sum::<f64>() / cells.len().max(1) as f64;
        s += &format!(
            "      A1/A2 {}: {} cells (net ≤ 0: {nonpos}) · pumping factor median {md:.2}, p90 {p9:.2}, max {mx:.1} · ΣDS / Σnet {:.2} · steps {} of which above the kinematic supply {:.0} % · excess per cell {exc:.3} s ({:.1} km of crust)\n",
            FAM_NAMES[f as usize], cells.len(), tot / net_tot, pm.steps[f as usize], 100.0 * pm.exceed[f as usize] as f64 / pm.steps[f as usize].max(1) as f64, exc * 35.0
        );
        s += &format!("#PF{f} {} {nonpos} {md} {p9} {mx} {} {} {} {exc}\n", cells.len(), tot / net_tot, pm.steps[f as usize], pm.exceed[f as usize]);
        // per-cell factors for pooling (cells with a positive net)
        s += &format!("#FAC{f} {}\n", cells.iter().filter(|&&k| s_end[k] - o.s0[k] > 1e-9).map(|&k| format!("{:.3}", (1..5).map(|g| pm.ds_fam[g][k]).sum::<f64>() / (s_end[k] - o.s0[k]))).collect::<Vec<_>>().join(" "));
    }
    // A3: s_src and its bounded estimate on the collision belts and on the active 0–3 band
    let s_src: Vec<f64> = (0..nn).map(|k| o.s0[k] + o.src_sum[k]).collect();
    let s_b: Vec<f64> = (0..nn).map(|k| s_src[k] - pm.excess[k]).collect();
    let mut lab = vec![0u8; nn];
    let mut dd = vec![u32::MAX; nn];
    let mut qd = std::collections::VecDeque::new();
    for k in 0..nn {
        if in_cont(k) && (is_active_cls(a.cls[k]) || is_passive_cls(a.cls[k])) {
            lab[k] = if is_active_cls(a.cls[k]) { 1 } else { 2 };
            dd[k] = 0;
            qd.push_back(k);
        }
    }
    while let Some(k) = qd.pop_front() {
        for j in nb4_64(k) {
            if pt[j] == PlateType::Continental && dd[j] == u32::MAX && dd[k] + 1 < 3 {
                dd[j] = dd[k] + 1;
                lab[j] = lab[k];
                qd.push_back(j);
            }
        }
    }
    let coll: Vec<usize> = (0..nn).filter(|&k| dom[k] == 1).collect();
    let act: Vec<usize> = (0..nn).filter(|&k| lab[k] == 1).collect();
    let st = |cells: &[usize], v: &[f64]| -> (f64, f64) {
        let mut x: Vec<f64> = cells.iter().map(|&k| v[k]).collect();
        (q_of(&mut x.clone(), 0.5), q_of(&mut x, 0.9))
    };
    let (c50, c90) = st(&coll, &s_src);
    let (cb50, cb90) = st(&coll, &s_b);
    let (a50, a90) = st(&act, &s_src);
    let (ab50, ab90) = st(&act, &s_b);
    s += &format!(
        "      A3: collision belts ({} cells) s_src median {c50:.2} (p90 {c90:.2}) → bounded {cb50:.2} ({cb90:.2}), −{:.0} % · active 0–3 band ({} cells) s_src {a50:.2} ({a90:.2}) → bounded {ab50:.2} ({ab90:.2}), −{:.0} %\n",
        coll.len(), 100.0 * (1.0 - cb50 / c50), act.len(), 100.0 * (1.0 - ab50 / a50)
    );
    s += &format!("#A3 {c50} {c90} {cb50} {cb90} {a50} {a90} {ab50} {ab90}\n");
    s += &format!("#A3C {}\n", coll.iter().map(|&k| format!("{:.3}:{:.3}", s_src[k], s_b[k])).collect::<Vec<_>>().join(" "));
    // B: the land, r*, A1, the bands, the seeds
    let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let circ = circumnavigation(&land64, 64);
    let land = land64.iter().filter(|&&l| l).count() as f64 / nn as f64;
    let share = |cls: &[u8]| -> f64 {
        let ac = (0..nn).filter(|&k| in_cont(k) && is_active_cls(cls[k])).count();
        let pa = (0..nn).filter(|&k| in_cont(k) && is_passive_cls(cls[k])).count();
        ac as f64 / (ac + pa).max(1) as f64
    };
    let (b, bc) = banded(&o, &a, &comp, &big);
    let mean06 = |t: usize, i: usize| (b[t][0][i] * bc[t][0] as f64 + b[t][1][i] * bc[t][1] as f64) / (bc[t][0] + bc[t][1]).max(1) as f64;
    let tot06 = |t: usize| (0..4).map(|i| mean06(t, i)).sum::<f64>();
    // the land added over F176 (F176's world)
    let w176 = world_oo(seed, p, &Margin::seaward_oo_collision());
    let l176: Vec<bool> = w176.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let added: Vec<usize> = (0..nn).filter(|&k| land64[k] && !l176[k]).collect();
    let na = added.len().max(1) as f64;
    let in_band = added.iter().filter(|&&k| o.band[k]).count() as f64 / na;
    let cont0 = added.iter().filter(|&&k| o.st0.plate_type.data()[k] == PlateType::Continental).count() as f64 / na;
    s += &format!(
        "      B: C1 end land {:.1} %, circumnavigable {} r* {:?} · land added over F176 {} cells ({:.0} % in the band at t₀, {:.0} % continental at t₀) · A1 any {:.1} %, majority {:.1} % · 0–6 cells: active {:.0} m (C-O wedge {:.0} m), passive {:.0} m · Andean seed instances {} ({:.0} % on cells emerged at t₀), cells seeded {} · walk depth from the trench (0..10, 11+) {:?}\n",
        100.0 * land, !circ.wraps, circ.sea_r, added.len(), 100.0 * in_band, 100.0 * cont0, 100.0 * share(&a.cls), 100.0 * share(&a.major), tot06(0), mean06(0, 0), tot06(1),
        pm.seed_total, 100.0 * pm.seed_emerged as f64 / pm.seed_total.max(1) as f64, pm.seeded.iter().filter(|&&b| b).count(), pm.seed_d
    );
    s += &format!(
        "#B {land} {} {} {} {in_band} {cont0} {} {} {} {} {} {} {} {}\n",
        !circ.wraps as u8, circ.sea_r.unwrap_or(0), added.len(), share(&a.cls), share(&a.major), tot06(0), mean06(0, 0), tot06(1), pm.seed_total, pm.seed_emerged, pm.seeded.iter().filter(|&&b| b).count()
    );
    s += &format!("#SD {}\n", pm.seed_d.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" "));
    // the images: the pumping factor by family (64², ×6), the seeds
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..nn {
        let (x, y) = (k % N64, k / N64);
        let net = s_end[k] - o.s0[k];
        let tot: f64 = (1..5).map(|g| pm.ds_fam[g][k]).sum();
        let base: [f64; 3] = match dom[k] {
            1 => [60.0, 110.0, 240.0],
            2 => [230.0, 40.0, 40.0],
            3 => [250.0, 160.0, 30.0],
            4 => [160.0, 70.0, 180.0],
            _ => if pt[k] == PlateType::Continental { [195.0, 195.0, 185.0] } else { [40.0, 70.0, 140.0] },
        };
        // brightness: the factor (log10, 1 → dim, ≥ 100 or net ≤ 0 → full)
        let t = if dom[k] == 0 { 1.0 } else if net <= 1e-9 { 1.0 } else { ((tot / net).max(1.0).log10() / 2.0).clamp(0.15, 1.0) };
        let px = if dom[k] == 0 { [base[0] as u8, base[1] as u8, base[2] as u8, 255] } else { [(base[0] * t) as u8, (base[1] * t) as u8, (base[2] * t) as u8, 255] };
        for dy in 0..6 {
            for dx in 0..6 {
                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(px));
            }
        }
    }
    img.save(out179().join(format!("pumping_{name}_seed{seed}.png"))).unwrap();
    let seeded: String = pm.seeded.iter().map(|&b| if b { '1' } else { '0' }).collect();
    s += &format!("#SEEDED {seeded}\n");
    s
}

fn seeds179() -> Vec<u64> {
    PANEL.iter().chain(SCREEN.iter()).copied().collect()
}

/// F179 -- C1, the 15 seeds of F176's selection × F175 / F176 / F178 / F179: the pumping (A1–A3) and the seed at the arc (B).
#[allow(dead_code)]
fn f179_c1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out179()).unwrap();
    eprintln!("\n==========  F179 . C1: the Davis-Suppe pumping and the Andean seed at the arc  ==========");
    std::thread::scope(|sc| {
        for seed in seeds179() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                for (name, m) in profiles179() {
                    cached179(&format!("c1_{name}_seed{seed}"), || c1_179(seed, p, name, &m));
                }
            });
        }
    });
    for (name, _) in profiles179() {
        eprintln!("\n   ── {name} ──");
        let mut pooled: Vec<Vec<f64>> = vec![Vec::new(); 5];
        let mut agg = [[0f64; 4]; 5];
        let mut a3 = Vec::new();
        let (mut lost, mut lands, mut anys, mut majs, mut added, mut emerged, mut tot_seeds) = (0, Vec::new(), Vec::new(), Vec::new(), Vec::new(), 0f64, 0f64);
        let mut coll_pairs: Vec<(f64, f64)> = Vec::new();
        let mut sd = [0f64; 12];
        for seed in seeds179() {
            let Ok(s) = std::fs::read_to_string(results179().join(format!("c1_{name}_seed{seed}.txt"))) else { continue };
            for f in 1..5 {
                if let Some(v) = data_line(&s, &format!("FAC{f}")) {
                    pooled[f].extend(v);
                }
                if let Some(v) = data_line(&s, &format!("PF{f}")) {
                    // #PF: cells, net ≤ 0, median, p90, max, ΣDS/Σnet, steps, steps above the supply, excess per cell
                    agg[f][0] += v[6];
                    agg[f][1] += v[7];
                    agg[f][2] += v[8] * v[0];
                    agg[f][3] += v[0];
                }
            }
            a3.push(data_line(&s, "A3").unwrap());
            if let Some(l) = s.lines().find(|l| l.starts_with("#A3C ")) {
                for t in l.split_whitespace().skip(1) {
                    if let Some((x, y)) = t.split_once(':') {
                        coll_pairs.push((x.parse().unwrap(), y.parse().unwrap()));
                    }
                }
            }
            let b = data_line(&s, "B").unwrap();
            lands.push(100.0 * b[0]);
            if b[1] == 0.0 {
                lost += 1;
            }
            added.push(b[3]);
            anys.push(100.0 * b[6]);
            majs.push(100.0 * b[7]);
            tot_seeds += b[11];
            emerged += b[12];
            for (i, v) in data_line(&s, "SD").unwrap().iter().enumerate() {
                sd[i] += v;
            }
        }
        for f in 1..5 {
            let mut v = pooled[f].clone();
            if v.is_empty() && agg[f][0] == 0.0 {
                continue;
            }
            eprintln!(
                "   {name} {}: pooled pumping factor (cells with a positive net) median {:.2}, p90 {:.2}, max {:.1} ({} cells) · steps above the kinematic supply {:.0} % · excess per family cell {:.3} s ({:.1} km)",
                FAM_NAMES[f], q_of(&mut v.clone(), 0.5), q_of(&mut v.clone(), 0.9), q_of(&mut v, 1.0), pooled[f].len(), 100.0 * agg[f][1] / agg[f][0].max(1.0), agg[f][2] / agg[f][3].max(1.0), 35.0 * agg[f][2] / agg[f][3].max(1.0)
            );
        }
        let mut cs: Vec<f64> = coll_pairs.iter().map(|p| p.0).collect();
        let mut cb: Vec<f64> = coll_pairs.iter().map(|p| p.1).collect();
        let (m0, m1) = (q_of(&mut cs, 0.5), q_of(&mut cb, 0.5));
        let mut act_s: Vec<f64> = a3.iter().map(|v| v[4]).collect();
        let mut act_b: Vec<f64> = a3.iter().map(|v| v[6]).collect();
        eprintln!(
            "   {name} A3: the collision belts' pooled s_src median {m0:.2} → bounded {m1:.2} (−{:.0} %) · the active 0–3 band, median of the seeds' medians {:.2} → {:.2}",
            100.0 * (1.0 - m1 / m0), q_of(&mut act_s, 0.5), q_of(&mut act_b, 0.5)
        );
        eprintln!(
            "   {name} B: land median {:.1} % (range {:.1}–{:.1}) · circumnavigation lost {lost} · land added over F176 median {:.0} cells · A1 any {:.1} %, majority {:.1} % · Andean seeds {:.0} instances, {:.0} % on cells emerged at t₀ · walk depth (0..10, 11+) {:?}",
            q_of(&mut lands.clone(), 0.5), q_of(&mut lands.clone(), 0.0), q_of(&mut lands, 1.0), q_of(&mut added, 0.5), q_of(&mut anys, 0.5), q_of(&mut majs, 0.5), tot_seeds, 100.0 * emerged / tot_seeds.max(1.0), sd
        );
    }
    // the images: the pumping maps of the panel (F176, F178, F179), the seeds F178 against F179
    let mut tiles = Vec::new();
    for name in ["F176", "F178", "F179"] {
        for seed in PANEL {
            if let Ok(i) = image::open(out179().join(format!("pumping_{name}_seed{seed}.png"))) {
                tiles.push(i.to_rgba8());
            }
        }
    }
    montage(&tiles, 5, 388, 388).save(out179().join("pumping_panel_F176_F178_F179.png")).unwrap();
    let mut tiles = Vec::new();
    for seed in PANEL {
        let rd = |n: &str| -> Option<(Vec<bool>, Vec<bool>)> {
            let s = std::fs::read_to_string(results179().join(format!("c1_{n}_seed{seed}.txt"))).ok()?;
            let l = s.lines().find(|l| l.starts_with("#SEEDED "))?;
            Some((l[8..].chars().map(|c| c == '1').collect(), Vec::new()))
        };
        let (Some((a, _)), Some((b, _))) = (rd("F178"), rd("F179")) else { continue };
        let p = sel176(seed).unwrap();
        let st0 = p.init_state_margin(64, seed, &Margin::seaward_oo_collision());
        let mut img = image::RgbaImage::new(384, 384);
        for k in 0..4096 {
            let (x, y) = (k % 64, k / 64);
            let emerged = st0.plate_type.data()[k] == PlateType::Continental && st0.s.data()[k] >= 0.915;
            let c = match (a[k], b[k]) {
                (true, true) => [170, 40, 200, 255],
                (true, false) => [230, 40, 40, 255],
                (false, true) => [40, 90, 240, 255],
                _ if emerged => [200, 185, 140, 255],
                _ if st0.plate_type.data()[k] == PlateType::Continental => [150, 170, 160, 255],
                _ => [40, 70, 140, 255],
            };
            for dy in 0..6 {
                for dx in 0..6 {
                    img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(c));
                }
            }
        }
        tiles.push(img);
    }
    montage(&tiles, 5, 388, 388).save(out179().join("seeds_F178_vs_F179_panel.png")).unwrap();
    eprintln!("   F179: {:.0} s", t0.elapsed().as_secs_f64());
}

// ════════════════════════════════════════ F180 ════════════════════════════════════════

fn out180() -> PathBuf {
    root().join("docs/reports/relief_method/f180_andean_mass/images")
}
fn results180() -> PathBuf {
    root().join("docs/reports/relief_method/f180_andean_mass/results")
}
fn cached180(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results180().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results180()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
/// The dominant Davis-Suppe family per cell over the run (0 none).
fn dominant(pm: &Pump) -> Vec<u8> {
    (0..N64 * N64)
        .map(|k| {
            let (f, v) = (1..5).map(|f| (f, pm.ds_fam[f][k])).fold((0usize, 0f64), |acc, x| if x.1 > acc.1 { x } else { acc });
            if v > 0.0 { f as u8 } else { 0 }
        })
        .collect()
}
/// F180 Part 0 -- the Andean coverage: the share of the continents' coast cells seeded at some step or within a periodic
/// Euclidean distance < 3 of a cell whose dominant family is Andean; with the coverage map (64², ×6).
fn coverage(o: &Obs, pm: &Pump, a: &A1) -> (f64, image::RgbaImage) {
    let nn = N64 * N64;
    let pt = o.st.plate_type.data();
    let (comp, big) = continents64(pt);
    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
    let dom = dominant(pm);
    let and: Vec<usize> = (0..nn).filter(|&k| dom[k] == 3).collect();
    let coast: Vec<usize> = (0..nn).filter(|&k| in_cont(k) && a.cls[k] >= 1 && a.cls[k] <= 5).collect();
    let near = |k: usize| -> bool {
        let (x, y) = ((k % N64) as i64, (k / N64) as i64);
        and.iter().any(|&q| {
            let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
            dx = dx.min(64 - dx);
            let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
            dy = dy.min(64 - dy);
            dx * dx + dy * dy < 9
        })
    };
    let covered: Vec<bool> = (0..nn).map(|k| coast.contains(&k) && (pm.seeded[k] || near(k))).collect();
    let cov = covered.iter().filter(|&&c| c).count() as f64 / coast.len().max(1) as f64;
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..nn {
        let (x, y) = (k % N64, k / N64);
        let c = if covered[k] {
            [220, 30, 30, 255]
        } else if coast.contains(&k) {
            if is_active_cls(a.major[k]) { [250, 220, 40, 255] } else { [40, 170, 70, 255] }
        } else if dom[k] == 3 {
            [245, 160, 40, 255]
        } else if pt[k] == PlateType::Continental {
            [195, 195, 185, 255]
        } else {
            [40, 70, 140, 255]
        };
        for dy in 0..6 {
            for dx in 0..6 {
                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(c));
            }
        }
    }
    (cov, img)
}
fn a1_shares(o: &Obs, a: &A1) -> (f64, f64) {
    let nn = N64 * N64;
    let (comp, big) = continents64(o.st.plate_type.data());
    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
    let share = |cls: &[u8]| -> f64 {
        let ac = (0..nn).filter(|&k| in_cont(k) && is_active_cls(cls[k])).count();
        let pa = (0..nn).filter(|&k| in_cont(k) && is_passive_cls(cls[k])).count();
        ac as f64 / (ac + pa).max(1) as f64
    };
    (share(&a.cls), share(&a.major))
}

/// F180 Part 0 -- the Andean coverage under F179 as committed, against the A1 shares (C1, the 15 seeds).
#[allow(dead_code)]
fn f180_part0_coverage() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out180()).unwrap();
    eprintln!("\n==========  F180 Part 0 . the Andean coverage under F179, against A1  ==========");
    let m = Margin::seaward_oo_collision_arc();
    std::thread::scope(|sc| {
        for seed in seeds179() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                cached180(&format!("p0_F179_seed{seed}"), || {
                    let (o, pm) = observe179(seed, p, &m);
                    let info = step_info(&o.st, &o.kin, o.ds_max, true);
                    let a = a1(&info, &o.kin, &o.tracker);
                    let (cov, img) = coverage(&o, &pm, &a);
                    let (any, maj) = a1_shares(&o, &a);
                    img.save(out180().join(format!("coverage_F179_seed{seed}.png"))).unwrap();
                    format!("   seed {seed} F179: the Andean coverage {:.1} % of the coast · A1 any {:.1} %, majority {:.1} % · coverage / majority {:.2}\n#COV {cov} {any} {maj}\n", 100.0 * cov, 100.0 * any, 100.0 * maj, cov / maj.max(1e-9))
                });
            });
        }
    });
    let (mut r2, mut over70, mut n) = (0, 0, 0);
    let mut ratios = Vec::new();
    for seed in seeds179() {
        if let Ok(s) = std::fs::read_to_string(results180().join(format!("p0_F179_seed{seed}.txt"))) {
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
            let v = data_line(&s, "COV").unwrap();
            n += 1;
            let r = v[0] / v[2].max(1e-9);
            ratios.push(r);
            r2 += (r >= 2.0) as usize;
            over70 += (v[0] > 0.7) as usize;
        }
    }
    eprintln!("\n   R0: coverage ≥ 2 × the majority A1 share on {r2} / {n} seeds · coverage above 70 % on {over70} · ratio median {:.2}", q_of(&mut ratios, 0.5));
    let tiles: Vec<image::RgbaImage> = PANEL.iter().filter_map(|s| image::open(out180().join(format!("coverage_F179_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 388, 388).save(out180().join("coverage_F179_panel.png")).unwrap();
    eprintln!("   Part 0: {:.0} s", t0.elapsed().as_secs_f64());
}

fn profiles180() -> Vec<(&'static str, Margin)> {
    let f179 = Margin::seaward_oo_collision_arc();
    vec![
        ("F176", Margin::seaward_oo_collision()),
        ("F179", f179),
        ("F179_A1", Margin { andean_active: Some((D_ARC_CELLS, ANDEAN_GAP_CELLS)), ..f179 }),
        ("F179_A2", Margin { ds_mass_bound: Some(1.0), ..f179 }),
        ("F180", Margin::seaward_oo_collision_f180()),
    ]
}
/// The pumping map (64², ×6): the dominant family coloured, its brightness the pumping factor (as F179).
fn pump_png(o: &Obs, pm: &Pump) -> image::RgbaImage {
    let nn = N64 * N64;
    let dom = dominant(pm);
    let pt = o.st.plate_type.data();
    let s_end = o.st.s.data();
    let mut img = image::RgbaImage::new(384, 384);
    for k in 0..nn {
        let (x, y) = (k % N64, k / N64);
        let net = s_end[k] - o.s0[k];
        let tot: f64 = (1..5).map(|g| pm.ds_fam[g][k]).sum();
        let base: [f64; 3] = match dom[k] {
            1 => [60.0, 110.0, 240.0],
            2 => [230.0, 40.0, 40.0],
            3 => [250.0, 160.0, 30.0],
            4 => [160.0, 70.0, 180.0],
            _ => if pt[k] == PlateType::Continental { [195.0, 195.0, 185.0] } else { [40.0, 70.0, 140.0] },
        };
        let t = if dom[k] == 0 || net <= 1e-9 { 1.0 } else { ((tot / net).max(1.0).log10() / 2.0).clamp(0.15, 1.0) };
        let px = [(base[0] * t) as u8, (base[1] * t) as u8, (base[2] * t) as u8, 255];
        for dy in 0..6 {
            for dx in 0..6 {
                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(px));
            }
        }
    }
    img
}

/// One F180 C1 analysis of a profile.
fn c1_180(seed: u64, p: C1Profile, name: &str, m: &Margin) -> String {
    let clk = Clock::start();
    let (o, pm) = observe179(seed, p, m);
    let w = world_oo(seed, p, m);
    let cost = clk.stop();
    let nn = N64 * N64;
    let info = step_info(&o.st, &o.kin, o.ds_max, true);
    let a = a1(&info, &o.kin, &o.tracker);
    let (cov, cov_img) = coverage(&o, &pm, &a);
    let (any, maj) = a1_shares(&o, &a);
    cov_img.save(out180().join(format!("coverage_{name}_seed{seed}.png"))).unwrap();
    pump_png(&o, &pm).save(out180().join(format!("pumping_{name}_seed{seed}.png"))).unwrap();
    let pt = o.st.plate_type.data();
    let (comp, big) = continents64(pt);
    let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let circ = circumnavigation(&land64, 64);
    let land = land64.iter().filter(|&&l| l).count() as f64 / nn as f64;
    let w176 = world_oo(seed, p, &Margin::seaward_oo_collision());
    let l176: Vec<bool> = w176.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
    let added: Vec<usize> = (0..nn).filter(|&k| land64[k] && !l176[k]).collect();
    let na = added.len().max(1) as f64;
    let in_band = added.iter().filter(|&&k| o.band[k]).count() as f64 / na;
    let cont0 = added.iter().filter(|&&k| o.st0.plate_type.data()[k] == PlateType::Continental).count() as f64 / na;
    let (b, bc) = banded(&o, &a, &comp, &big);
    let mean06 = |t: usize, i: usize| (b[t][0][i] * bc[t][0] as f64 + b[t][1][i] * bc[t][1] as f64) / (bc[t][0] + bc[t][1]).max(1) as f64;
    let tot06 = |t: usize| (0..4).map(|i| mean06(t, i)).sum::<f64>();
    // s_src by the dominant family
    let dom = dominant(&pm);
    let s_src: Vec<f64> = (0..nn).map(|k| o.s0[k] + o.src_sum[k]).collect();
    let mut sf = [f64::NAN; 5];
    let mut nf = [0usize; 5];
    for f in 1..5u8 {
        let mut v: Vec<f64> = (0..nn).filter(|&k| dom[k] == f).map(|k| s_src[k]).collect();
        nf[f as usize] = v.len();
        sf[f as usize] = q_of(&mut v, 0.5);
    }
    let mut s = format!(
        "   seed {seed} {name}: the Andean coverage {:.1} % · A1 any {:.1} %, majority {:.1} % · coverage / majority {:.2} · C1 end land {:.1} %, circumnavigable {} r* {:?} · land added over F176 {} cells ({:.0} % in the band, {:.0} % continental at t₀) · 0–6 cells: active {:.0} m (C-O wedge {:.0} m), passive {:.0} m · C1 CPU {:.1} s\n",
        100.0 * cov, 100.0 * any, 100.0 * maj, cov / maj.max(1e-9), 100.0 * land, !circ.wraps, circ.sea_r, added.len(), 100.0 * in_band, 100.0 * cont0, tot06(0), mean06(0, 0), tot06(1), cost.1
    );
    for f in 1..5 {
        s += &format!(
            "      {}: {} cells, s_src median {:.2} · seed-steps {} of which bounded {:.0} % · mass removed {:.2} (s × cells; {:.2} km of crust per family cell)\n",
            FAM_NAMES[f], nf[f], sf[f], pm.bnd_seeds[f], 100.0 * pm.bnd_bounded[f] as f64 / pm.bnd_seeds[f].max(1) as f64, pm.removed_fam[f], 35.0 * pm.removed_fam[f] / nf[f].max(1) as f64
        );
    }
    s += &format!("#COV {cov} {any} {maj}\n");
    s += &format!("#B {land} {} {} {} {in_band} {cont0}\n", !circ.wraps as u8, circ.sea_r.unwrap_or(0), added.len());
    s += &format!("#UP {} {} {}\n", tot06(0), mean06(0, 0), tot06(1));
    s += &format!("#SF {} {} {} {}\n", sf[1], sf[2], sf[3], sf[4]);
    s += &format!("#BND {}\n", (1..5).map(|f| format!("{} {} {}", pm.bnd_seeds[f], pm.bnd_bounded[f], pm.removed_fam[f])).collect::<Vec<_>>().join(" "));
    s
}

/// F180 tier 1 -- C1, the 15 seeds: F176, F179, F179 + A1, F179 + A2, F180.
#[allow(dead_code)]
fn f180_tier1_c1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out180()).unwrap();
    eprintln!("\n==========  F180 tier 1 . C1: the coverage, the land, the bands, s_src, the bound  ==========");
    std::thread::scope(|sc| {
        for seed in seeds179() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                for (name, m) in profiles180() {
                    cached180(&format!("t1_{name}_seed{seed}"), || c1_180(seed, p, name, &m));
                }
            });
        }
    });
    let rd = |name: &str, seed: u64| std::fs::read_to_string(results180().join(format!("t1_{name}_seed{seed}.txt"))).ok();
    for (name, _) in profiles180() {
        let (mut in_band_r, mut over70, mut lost, mut r0) = (0, 0, 0, 0);
        let (mut lands, mut ups, mut upps, mut covs, mut ratios, mut anys, mut majs) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut coll_ratio = Vec::new();
        let mut bnd = [[0f64; 3]; 5];
        for seed in seeds179() {
            let Some(s) = rd(name, seed) else { continue };
            if name == "F180" || name == "F179" {
                eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
            }
            let c = data_line(&s, "COV").unwrap();
            let b = data_line(&s, "B").unwrap();
            let u = data_line(&s, "UP").unwrap();
            let sf = data_line(&s, "SF").unwrap();
            let bn = data_line(&s, "BND").unwrap();
            let r = c[0] / c[2].max(1e-9);
            covs.push(100.0 * c[0]);
            ratios.push(r);
            anys.push(100.0 * c[1]);
            majs.push(100.0 * c[2]);
            in_band_r += (0.8..=1.5).contains(&r) as usize;
            r0 += (r >= 2.0) as usize;
            over70 += (c[0] > 0.7) as usize;
            lost += (b[1] == 0.0) as usize;
            lands.push(100.0 * b[0]);
            ups.push(u[0]);
            upps.push(u[2]);
            if let Some(s176) = rd("F176", seed) {
                let sf176 = data_line(&s176, "SF").unwrap();
                coll_ratio.push(sf[0] / sf176[0]);
            }
            for f in 1..5 {
                for i in 0..3 {
                    bnd[f][i] += bn[(f - 1) * 3 + i];
                }
            }
        }
        eprintln!(
            "\n   {name}: the Andean coverage median {:.1} % (> 70 % on {over70}) · coverage / majority median {:.2} ({r0} ≥ 2; {in_band_r} in 0.8–1.5) · A1 any {:.1} %, majority {:.1} % · land median {:.1} % ({:.1}–{:.1}) · circumnavigation lost {lost} · active 0–6 cells median {:.0} m, passive {:.0} m · collision s_src median / F176 {:.3} (median over seeds) · bounded seed-steps: {}\n",
            q_of(&mut covs, 0.5), q_of(&mut ratios, 0.5), q_of(&mut anys, 0.5), q_of(&mut majs, 0.5), q_of(&mut lands.clone(), 0.5), q_of(&mut lands.clone(), 0.0), q_of(&mut lands, 1.0), q_of(&mut ups, 0.5), q_of(&mut upps, 0.5), q_of(&mut coll_ratio, 0.5),
            (1..5).map(|f| format!("{} {:.0} % (removed {:.1})", FAM_NAMES[f], 100.0 * bnd[f][1] / bnd[f][0].max(1.0), bnd[f][2])).collect::<Vec<_>>().join(" · ")
        );
    }
    for (pre, tag) in [("coverage", "coverage_panel_F179_F180.png"), ("pumping", "pumping_panel_F179_F180.png")] {
        let mut tiles = Vec::new();
        for name in ["F179", "F180"] {
            for seed in PANEL {
                if let Ok(i) = image::open(out180().join(format!("{pre}_{name}_seed{seed}.png"))) {
                    tiles.push(i.to_rgba8());
                }
            }
        }
        montage(&tiles, 5, 388, 388).save(out180().join(tag)).unwrap();
    }
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F180 tier 2 -- 256², young, the 15 seeds: F176 and F180 (B3, the formations, W₁₀, the macro).
#[allow(dead_code)]
fn f180_tier2_physics() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out180()).unwrap();
    eprintln!("\n==========  F180 tier 2 . 256² young: B3, the formations, W10, the macro  ==========");
    let rough_a = rough_a177();
    for (name, m) in [("F176", Margin::seaward_oo_collision()), ("F180", Margin::seaward_oo_collision_f180())] {
        for seed in seeds179() {
            let Some(p) = sel176(seed) else { continue };
            let s = cached180(&format!("t2_{name}_seed{seed}"), || {
                let (o, _) = observe179(seed, p, &m);
                let info = step_info(&o.st, &o.kin, o.ds_max, true);
                let a = a1(&info, &o.kin, &o.tracker);
                let (r, w) = physics_young(seed, p, &m, rough_a);
                let (line, form) = formations178(&o, &a, &r, &w);
                let org = origin64(&o);
                // the relief at 256² with A1's majority-active coasts (red)
                let mut rel = alt_png_plain(&r.zm);
                for k in 0..65536 {
                    let (x, y) = (k % 256, k / 256);
                    if is_active_cls(a.major[(y / 4) * 64 + x / 4]) {
                        rel.put_pixel(x as u32, (255 - y) as u32, image::Rgba([230, 30, 30, 255]));
                    }
                }
                rel.save(out180().join(format!("relief_active_{name}_seed{seed}.png"))).unwrap();
                origin256_png(&r.zm, &r.classes, &org, None).save(out180().join(format!("origin_{name}_seed{seed}.png"))).unwrap();
                origin256_png(&r.zm, &r.classes, &org, Some(&form)).save(out180().join(format!("formations_{name}_seed{seed}.png"))).unwrap();
                format!("   seed {seed} {name} young: physics CPU {:.0} s\n{line}", r.cost.1)
            });
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        }
        let (mut with_cord, mut only, mut n, mut steep) = (0, 0, 0, 0);
        let mut w10_panel = Vec::new();
        let mut mtn = Vec::new();
        for seed in seeds179() {
            if let Ok(s) = std::fs::read_to_string(results180().join(format!("t2_{name}_seed{seed}.txt"))) {
                let f = data_line(&s, "F").unwrap();
                let wv = data_line(&s, "W").unwrap();
                let mv = data_line(&s, "M").unwrap();
                n += 1;
                with_cord += (f[0] > 0.0) as usize;
                only += (f[8] > 0.0) as usize;
                if PANEL.contains(&seed) {
                    w10_panel.push((seed, (wv[0] * 1000.0).round() / 1000.0, (wv[1] * 1000.0).round() / 1000.0, (wv[2] * 1000.0).round() / 1000.0));
                    steep += (wv[1] >= 2.0 * wv[2]) as usize;
                }
                mtn.push((seed, (mv[3] * 1000.0).round() / 10.0));
            }
        }
        eprintln!("\n   {name}: seeds with ≥ 1 Andean cordillera {with_cord} / {n} · collision the only formation on {only} · W10 young on the panel (seed, global, active, passive) {w10_panel:?} · active ≥ 2 × passive on {steep} panel seeds · mountain % {mtn:?}\n");
        for (pre, tag) in [("formations", "formations_15"), ("origin", "origin_15"), ("relief_active", "relief_active_15")] {
            let tiles: Vec<image::RgbaImage> = seeds179().iter().filter_map(|s| image::open(out180().join(format!("{pre}_{name}_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
            montage(&tiles, 5, 260, 260).save(out180().join(format!("{tag}_{name}.png"))).unwrap();
        }
    }
    eprintln!("   tier 2: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F180 tier 3 -- 2 048², young: the unclipped transects of seed 18 and of the seed with the longest new cordillera, and that
/// seed's whole map.
#[allow(dead_code)]
fn f180_tier3_cascade() {
    let t0 = Instant::now();
    eprintln!("\n==========  F180 tier 3 . 2 048² young: the transects without a clamp, a whole map  ==========");
    let m = Margin::seaward_oo_collision_f180();
    let cells = |seed: u64, p: C1Profile| -> (Option<usize>, Option<usize>) {
        let (o, _) = observe179(seed, p, &m);
        let info = step_info(&o.st, &o.kin, o.ds_max, true);
        let a = a1(&info, &o.kin, &o.tracker);
        let act = (0..4096).filter(|&k| is_active_cls(a.major[k])).max_by(|&x, &y| o.w_sum[x].total_cmp(&o.w_sum[y]).then(y.cmp(&x)));
        let pas = (0..4096).filter(|&k| is_passive_cls(a.cls[k])).min_by(|&x, &y| o.src_sum[x].total_cmp(&o.src_sum[y]).then(x.cmp(&y)));
        (act, pas)
    };
    let cord = |name: &str, seed: u64| -> f64 { std::fs::read_to_string(results180().join(format!("t2_{name}_seed{seed}.txt"))).ok().and_then(|s| data_line(&s, "F")).map(|f| f[1]).unwrap_or(0.0) };
    let pick = seeds179().into_iter().filter(|&sd| sd != 18).max_by(|&x, &y| {
        let k = |sd: u64| (cord("F176", sd) == 0.0, cord("F180", sd));
        let (kx, ky) = (k(x), k(y));
        kx.0.cmp(&ky.0).then(kx.1.total_cmp(&ky.1))
    });
    let mut todo = vec![(18u64, false)];
    if let Some(sd) = pick {
        todo.push((sd, true));
    }
    for (seed, with_map) in todo {
        let p = sel176(seed).unwrap();
        let s = cached180(&format!("t3_seed{seed}"), || {
            let (a, pz) = cells(seed, p);
            let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
            let corse_oct = |nc: usize| -> f32 {
                let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
                read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
            };
            let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
            let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
            let kh = KHetero::declared();
            let w = world_oo(seed, p, &m);
            let cover = cover_map_m(seed, p, &m, COVER_M);
            let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
            let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k.clone()) };
            let cfg = base_cfg(seed, rough_a, Some([0, 0]));
            let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
            let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[1].1, &rock, N_PHYS);
            let (g, _, _) = chain_low(&r.zm, amp, Some(ln_k), 0.0);
            if with_map {
                alt_png(&g, 4).save(out180().join(format!("map2048_F180_seed{seed}.png"))).unwrap();
            }
            let s64 = w.st.s.data().to_vec();
            let pt = w.st.plate_type.data().to_vec();
            let (mut zs, mut ss, mut desc) = (Vec::new(), Vec::new(), Vec::new());
            for (nm, cell, col) in [("passive", pz, [30u8, 140, 60]), ("active", a, [220, 30, 30])] {
                if let Some(k) = cell {
                    let (za, sb) = transect(&g, &s64, &pt, k);
                    let mx = za.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                    desc.push(format!("{nm} (cell {k}): max {mx:.0} m"));
                    zs.push((za, col));
                    ss.push((sb, col));
                }
            }
            let (zmin, zmax) = transect_png_free(&zs, &ss, &out180().join(format!("transect_F180_seed{seed}.png")));
            let peak = g.data.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            format!("   seed {seed} F180 young (2 048²){}: the transect's range {zmin:.0} … {zmax:.0} m · {} · the world's maximum {peak:.0} m\n", if with_map { " (the longest new cordillera; whole map saved)" } else { "" }, desc.join(" · "))
        });
        eprint!("{s}");
    }
    eprintln!("   tier 3: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F180's stop, its mechanism (C1, F180, the 15 seeds; measures only): on the coast cells the coverage counts, the active
/// (majority) against the passive ones -- the Andean Davis-Suppe they and their 3-cell neighbourhood received, and their
/// distance to the nearest cell ever seeded.
#[allow(dead_code)]
fn f180_stop_mechanism() {
    let m = Margin::seaward_oo_collision_f180();
    std::thread::scope(|sc| {
        for seed in seeds179() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                cached180(&format!("stopmech_seed{seed}"), || {
                    let nn = N64 * N64;
                    let (o, pm) = observe179(seed, p, &m);
                    let info = step_info(&o.st, &o.kin, o.ds_max, true);
                    let a = a1(&info, &o.kin, &o.tracker);
                    let pt = o.st.plate_type.data();
                    let (comp, big) = continents64(pt);
                    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
                    let dom = dominant(&pm);
                    let near = |k: usize, r2: i64| -> Vec<usize> {
                        let (x, y) = ((k % N64) as i64, (k / N64) as i64);
                        (0..nn).filter(|&q| {
                            let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
                            dx = dx.min(64 - dx);
                            let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
                            dy = dy.min(64 - dy);
                            dx * dx + dy * dy < r2
                        }).collect()
                    };
                    let seeded: Vec<usize> = (0..nn).filter(|&k| pm.seeded[k]).collect();
                    let dseed = |k: usize| -> f64 {
                        let (x, y) = ((k % N64) as i64, (k / N64) as i64);
                        seeded.iter().map(|&q| {
                            let mut dx = ((q % N64) as i64 - x).rem_euclid(64);
                            dx = dx.min(64 - dx);
                            let mut dy = ((q / N64) as i64 - y).rem_euclid(64);
                            dy = dy.min(64 - dy);
                            ((dx * dx + dy * dy) as f64).sqrt()
                        }).fold(f64::INFINITY, f64::min)
                    };
                    let mut out = String::new();
                    for (tag, sel) in [("active (majority)", 1u8), ("passive or mixed", 2u8)] {
                        let cells: Vec<usize> = (0..nn)
                            .filter(|&k| in_cont(k) && a.cls[k] >= 1 && a.cls[k] <= 5 && (is_active_cls(a.major[k]) == (sel == 1)))
                            .filter(|&k| pm.seeded[k] || near(k, 9).iter().any(|&q| dom[q] == 3))
                            .collect();
                        let n = cells.len().max(1) as f64;
                        let and_ds: f64 = cells.iter().map(|&k| near(k, 9).iter().map(|&q| pm.ds_fam[3][q]).sum::<f64>() / near(k, 9).len() as f64).sum::<f64>() / n;
                        let all_ds: f64 = cells.iter().map(|&k| near(k, 9).iter().map(|&q| (1..5).map(|f| pm.ds_fam[f][q]).sum::<f64>()).sum::<f64>() / near(k, 9).len() as f64).sum::<f64>() / n;
                        let mut ds: Vec<f64> = cells.iter().map(|&k| dseed(k)).collect();
                        out += &format!(
                            "   seed {seed} F180 covered {tag} coast cells {}: the Andean Σ DS in their 3-cell disc {:.3} s ({:.1} km of crust; every family {:.3} s) · the distance to the nearest seeded cell median {:.1} cells (p90 {:.1})\n",
                            cells.len(), and_ds, 35.0 * and_ds, all_ds, q_of(&mut ds.clone(), 0.5), q_of(&mut ds, 0.9)
                        );
                    }
                    out
                });
            });
        }
    });
    for seed in seeds179() {
        if let Ok(s) = std::fs::read_to_string(results180().join(format!("stopmech_seed{seed}.txt"))) {
            eprint!("{s}");
        }
    }
}

/// F180 diagnostic (C1, F180, the 15 seeds): the fragmentation of the majority-active coasts, option 5's polarity changes along
/// each subducting pair's boundary, the decisions on |Δage| < 0.01, and the net Andean uplift along the coasts (no threshold).
#[allow(dead_code)]
fn f180_diag_c1() {
    let m = Margin::seaward_oo_collision_f180();
    std::thread::scope(|sc| {
        for seed in seeds179() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                cached180(&format!("diag_seed{seed}"), || {
                    let nn = N64 * N64;
                    let (o, pm) = observe179(seed, p, &m);
                    let info = step_info(&o.st, &o.kin, o.ds_max, true);
                    let a = a1(&info, &o.kin, &o.tracker);
                    let pt = o.st.plate_type.data();
                    let pid = o.st.plate_id.data();
                    let (comp, big) = continents64(pt);
                    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
                    let coast: Vec<bool> = (0..nn).map(|k| in_cont(k) && a.cls[k] >= 1 && a.cls[k] <= 5).collect();
                    let act: Vec<bool> = (0..nn).map(|k| coast[k] && is_active_cls(a.major[k])).collect();
                    // the active segments (8-connected) and their gaps along the coast (8-connected walk through coast cells)
                    let segs = comps8(&act, N64);
                    let mut seg_of = vec![usize::MAX; nn];
                    for (i, sg) in segs.iter().enumerate() {
                        for &k in sg {
                            seg_of[k] = i;
                        }
                    }
                    let mut lens: Vec<f64> = segs.iter().map(|sg| sg.len() as f64).collect();
                    let mut gaps = Vec::new();
                    for (i, sg) in segs.iter().enumerate() {
                        let mut d = vec![usize::MAX; nn];
                        let mut q = std::collections::VecDeque::new();
                        for &k in sg {
                            d[k] = 0;
                            q.push_back(k);
                        }
                        let mut found = None;
                        while let Some(k) = q.pop_front() {
                            if seg_of[k] != usize::MAX && seg_of[k] != i {
                                found = Some(d[k] - 1);
                                break;
                            }
                            for j in nb8_n(k, N64) {
                                if coast[j] && d[j] == usize::MAX {
                                    d[j] = d[k] + 1;
                                    q.push_back(j);
                                }
                            }
                        }
                        if let Some(g) = found {
                            gaps.push(g as f64);
                        }
                    }
                    // option 5's polarity along each subducting pair's boundary
                    let partners = overriding_trench_partners(&o.st.plate_id, &o.st.plate_type, &o.st.age, &info.bi, &to_boolfield(&info.upper_vel), &o.kin, &o.tracker);
                    let mut pairs: std::collections::BTreeMap<(u16, u16), (Vec<bool>, Vec<bool>)> = std::collections::BTreeMap::new();
                    for k in 0..nn {
                        if let Some((q, _)) = partners[k] {
                            let (pa, pb) = (pid[k], pid[q]);
                            let key = (pa.min(pb), pa.max(pb));
                            let e = pairs.entry(key).or_insert((vec![false; nn], vec![false; nn]));
                            if pa == key.0 { e.0[k] = true } else { e.1[k] = true }
                        }
                    }
                    let mut rates = Vec::new();
                    let (mut ch_tot, mut len_tot) = (0usize, 0usize);
                    for (_, (oa, ob)) in &pairs {
                        let un: Vec<bool> = (0..nn).map(|k| oa[k] || ob[k]).collect();
                        let changes = comps8(oa, N64).len() + comps8(ob, N64).len() - comps8(&un, N64).len();
                        let len = un.iter().filter(|&&b| b).count();
                        ch_tot += changes;
                        len_tot += len;
                        rates.push(10.0 * changes as f64 / len.max(1) as f64);
                    }
                    let close = pm.oo_age_close as f64 / (pm.oo_age + pm.oo_tie).max(1) as f64;
                    // the net Andean uplift along the coasts: C1's net Δh × the Andean share of the Davis-Suppe in the 3-cell disc
                    let airy = Airy::declared();
                    let s_end = o.st.s.data();
                    let disc = |k: usize| -> Vec<usize> {
                        let (x, y) = ((k % N64) as i64, (k / N64) as i64);
                        let mut v = Vec::new();
                        for dy in -2i64..=2 {
                            for dx in -2i64..=2 {
                                if dx * dx + dy * dy < 9 {
                                    v.push((((y + dy).rem_euclid(64)) as usize) * N64 + ((x + dx).rem_euclid(64)) as usize);
                                }
                            }
                        }
                        v
                    };
                    let (mut na, mut np) = (Vec::new(), Vec::new());
                    for k in 0..nn {
                        if !coast[k] {
                            continue;
                        }
                        let dk = disc(k);
                        let and: f64 = dk.iter().map(|&q| pm.ds_fam[3][q]).sum();
                        let all: f64 = dk.iter().map(|&q| (1..5).map(|f| pm.ds_fam[f][q]).sum::<f64>()).sum();
                        let share = if all > 0.0 { and / all } else { 0.0 };
                        let net = (airy.altitude_km(s_end[k]) - airy.altitude_km(o.s0[k])) * 1000.0;
                        if act[k] { na.push(net * share) } else { np.push(net * share) }
                    }
                    let qs = |v: &mut Vec<f64>| format!("p10 {:.0} · p50 {:.0} · p90 {:.0} m ({} cells)", q_of(&mut v.clone(), 0.1), q_of(&mut v.clone(), 0.5), q_of(&mut v.clone(), 0.9), v.len());
                    format!(
                        "   seed {seed} F180: majority-active coast segments {} · lengths (cells) median {:.0}, p90 {:.0}, max {:.0} · gaps between neighbouring segments (coast cells) median {:.0}, p90 {:.0} · subducting pairs {} · polarity changes {ch_tot} over {len_tot} overriding trench cells ({:.2} per 10 cells; per pair median {:.2}) · option 5 decisions on |Δage| < 0.01 {:.0} % · net Andean uplift on the coast: active {} · passive {}\n#DG {} {} {} {} {} {} {}\n",
                        segs.len(), q_of(&mut lens.clone(), 0.5), q_of(&mut lens.clone(), 0.9), q_of(&mut lens.clone(), 1.0), q_of(&mut gaps.clone(), 0.5), q_of(&mut gaps.clone(), 0.9), pairs.len(),
                        10.0 * ch_tot as f64 / len_tot.max(1) as f64, q_of(&mut rates.clone(), 0.5), 100.0 * close, qs(&mut na), qs(&mut np),
                        q_of(&mut lens, 0.5), q_of(&mut gaps, 0.5), 10.0 * ch_tot as f64 / len_tot.max(1) as f64, close, q_of(&mut na.clone(), 0.5), q_of(&mut np.clone(), 0.5), q_of(&mut rates, 0.5)
                    )
                });
            });
        }
    });
    let mut v: Vec<Vec<f64>> = Vec::new();
    for seed in seeds179() {
        if let Ok(s) = std::fs::read_to_string(results180().join(format!("diag_seed{seed}.txt"))) {
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
            v.push(data_line(&s, "DG").unwrap());
        }
    }
    let col = |i: usize| -> Vec<f64> { v.iter().map(|r| r[i]).collect() };
    eprintln!(
        "\n   D4: the median active segment (median over seeds) {:.1} cells · polarity changes per 10 cells (median over seeds) {:.2} · |Δage| < 0.01 share median {:.0} % · net Andean uplift medians: active {:.0} m, passive {:.0} m",
        q_of(&mut col(0), 0.5), q_of(&mut col(2), 0.5), 100.0 * q_of(&mut col(3), 0.5), q_of(&mut col(4), 0.5), q_of(&mut col(5), 0.5)
    );
}

// ════════════════════════════════════════ F181 ════════════════════════════════════════

fn out181() -> PathBuf {
    root().join("docs/reports/relief_method/f181_adoption/images")
}
fn results181() -> PathBuf {
    root().join("docs/reports/relief_method/f181_adoption/results")
}
fn cached181(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results181().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results181()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}
fn physics_age(seed: u64, p: C1Profile, m: &Margin, rough_a: f32, age: usize) -> (Run, World) {
    let kh = KHetero::declared();
    let w = world_oo(seed, p, m);
    let cover = cover_map_m(seed, p, m, COVER_M);
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k) };
    let cfg = base_cfg(seed, rough_a, Some([0, 0]));
    let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[age].1, &rock, N_PHYS);
    (r, w)
}
/// The cordilleras' widths perpendicular to the coast: per cordillera group (`form == 1`), its area / its length along the coast
/// (the A1-active coast cells it touches × 6.25 km), world km.
fn cord_widths(form: &[u8], a: &A1) -> Vec<f64> {
    let n = 256;
    let mask: Vec<bool> = form.iter().map(|&f| f == 1).collect();
    comps8(&mask, n)
        .iter()
        .filter_map(|g| {
            let mut touched = std::collections::BTreeSet::new();
            for &k in g {
                let mut cand = vec![k];
                cand.extend(nb8_n(k, n));
                for c in cand {
                    let c64 = ((c / n) / 4) * 64 + (c % n) / 4;
                    if is_active_cls(a.cls[c64]) {
                        touched.insert(c64);
                    }
                }
            }
            (!touched.is_empty()).then(|| g.len() as f64 * (CELL * CELL) as f64 / (touched.len() as f64 * CELL64_KM as f64))
        })
        .collect()
}

/// F181 tier 1 -- 256², the 15 seeds × F176 / F180 × intermediate / old.
#[allow(dead_code)]
fn f181_tier1_ages() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out181()).unwrap();
    eprintln!("\n==========  F181 tier 1 . 256²: F176 and F180 at the intermediate and old ages  ==========");
    let rough_a = rough_a177();
    ages_runs(&[("F176", Margin::seaward_oo_collision()), ("F180", Margin::seaward_oo_collision_f180())], rough_a);
    ages_summary(&["F176", "F180"]);
    for name in ["F176", "F180"] {
        for ai in 2..4 {
            for pre in ["formations", "relief_active"] {
                let tiles: Vec<image::RgbaImage> = seeds179().iter().filter_map(|s| image::open(out181().join(format!("{pre}_{name}_age{ai}_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
                montage(&tiles, 5, 260, 260).save(out181().join(format!("{pre}_15_{name}_age{ai}.png"))).unwrap();
            }
        }
    }
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// The per-seed, per-age 256² runs of the profiles (resumable).
fn ages_runs(profiles: &[(&str, Margin)], rough_a: f32) {
    for &(name, m) in profiles {
        for seed in seeds179() {
            let Some(p) = sel176(seed) else { continue };
            if (2..4).all(|ai| results181().join(format!("t1_{name}_age{ai}_seed{seed}.txt")).exists()) {
                continue;
            }
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            for ai in 2..4 {
                cached181(&format!("t1_{name}_age{ai}_seed{seed}"), || {
                    let (r, w) = physics_age(seed, p, &m, rough_a, ai);
                    let (line, form) = formations178(&o, &a, &r, &w);
                    let mut widths = cord_widths(&form, &a);
                    let land: Vec<bool> = r.zm.data.iter().map(|&z| z > 0.0).collect();
                    let kc = coast_torus(&land, 256, 256, 4);
                    let circ = circumnavigation(&land, 256);
                    let org = origin64(&o);
                    origin256_png(&r.zm, &r.classes, &org, Some(&form)).save(out181().join(format!("formations_{name}_age{ai}_seed{seed}.png"))).unwrap();
                    let mut rel = alt_png_plain(&r.zm);
                    for k in 0..65536 {
                        let (x, y) = (k % 256, k / 256);
                        if is_active_cls(a.major[(y / 4) * 64 + x / 4]) {
                            rel.put_pixel(x as u32, (255 - y) as u32, image::Rgba([230, 30, 30, 255]));
                        }
                    }
                    rel.save(out181().join(format!("relief_active_{name}_age{ai}_seed{seed}.png"))).unwrap();
                    let wmed = q_of(&mut widths.clone(), 0.5);
                    format!(
                        "   seed {seed} {name} {}: physics CPU {:.0} s · cordillera widths (km world) median {wmed:.1} ({:.0} km terrestrial), range {:.1}–{:.1} · the coast's dimension {:.3} · the 256² land circumnavigable {} r* {:?}\n{line}#X {wmed} {} {} {}\n",
                        AGES[ai].0, r.cost.1, wmed * 7.5, q_of(&mut widths.clone(), 0.0), q_of(&mut widths, 1.0), kc.dimension, !circ.wraps, circ.sea_r, kc.dimension, !circ.wraps as u8, circ.sea_r.unwrap_or(0)
                    )
                });
            }
        }
    }
}

/// The summary lines of the profiles at the two ages.
fn ages_summary(names: &[&str]) {
    for &name in names {
        for ai in 2..4 {
            let (mut cord, mut only, mut n, mut lost, mut steep) = (0, 0, 0, 0, 0);
            let (mut w10p, mut mtn, mut wid, mut med, mut lnd, mut dim) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
            for seed in seeds179() {
                let Ok(s) = std::fs::read_to_string(results181().join(format!("t1_{name}_age{ai}_seed{seed}.txt"))) else { continue };
                let f = data_line(&s, "F").unwrap();
                let wv = data_line(&s, "W").unwrap();
                let mv = data_line(&s, "M").unwrap();
                let x = data_line(&s, "X").unwrap();
                n += 1;
                cord += (f[0] > 0.0) as usize;
                only += (f[8] > 0.0) as usize;
                lost += (x[2] == 0.0) as usize;
                if PANEL.contains(&seed) {
                    w10p.push((seed, (wv[0] * 1000.0).round() / 1000.0, (wv[1] * 1000.0).round() / 1000.0, (wv[2] * 1000.0).round() / 1000.0));
                    steep += (wv[1] >= 2.0 * wv[2]) as usize;
                }
                mtn.push(100.0 * mv[3]);
                med.push(mv[5]);
                lnd.push(mv[4]);
                dim.push(x[1]);
                if f[0] > 0.0 {
                    wid.push(x[0]);
                }
            }
            eprintln!(
                "\n   {name} {}: cordilleras on {cord} / {n} · collision the only formation on {only} · W10 panel (seed, global, active, passive) {w10p:?} · active ≥ 2 × passive on {steep} · mountain % median {:.1} (range {:.1}–{:.1}) · land median {:.0} km² · median altitude {:.0} m · coast dimension median {:.3} · cordillera width median of the seeds {:.1} km world ({:.0} terrestrial) · 256² circumnavigation lost {lost}\n",
                AGES[ai].0, q_of(&mut mtn.clone(), 0.5), q_of(&mut mtn.clone(), 0.0), q_of(&mut mtn, 1.0), q_of(&mut lnd, 0.5), q_of(&mut med, 0.5), q_of(&mut dim, 0.5), q_of(&mut wid.clone(), 0.5), 7.5 * q_of(&mut wid, 0.5)
            );
        }
    }
}

/// F181 Part C -- the « équerres »: the majority-active coast cells at sea (their C1 cell and its 8 neighbours sea), their groups,
/// their share, and the origin of their continental type (C1, F180, the 15 seeds).
#[allow(dead_code)]
fn f181_equerres() {
    let m = Margin::seaward_oo_collision_f180();
    std::thread::scope(|sc| {
        for seed in seeds179() {
            sc.spawn(move || {
                let Some(p) = sel176(seed) else { return };
                cached181(&format!("eq_seed{seed}"), || {
                    let nn = N64 * N64;
                    let (o, _) = observe179(seed, p, &m);
                    let w = world_oo(seed, p, &m);
                    let info = step_info(&o.st, &o.kin, o.ds_max, true);
                    let a = a1(&info, &o.kin, &o.tracker);
                    let pt = o.st.plate_type.data();
                    let pt0 = o.st0.plate_type.data();
                    let (comp, big) = continents64(pt);
                    let in_cont = |k: usize| comp[k] != usize::MAX && big.contains(&comp[k]);
                    let land64: Vec<bool> = w.src.hist.raw.last().unwrap().data.iter().map(|&v| v > 0.0).collect();
                    let act: Vec<bool> = (0..nn).map(|k| in_cont(k) && a.cls[k] >= 1 && a.cls[k] <= 5 && is_active_cls(a.major[k])).collect();
                    let off: Vec<bool> = (0..nn).map(|k| act[k] && !land64[k] && nb8_n(k, N64).iter().all(|&j| !land64[j])).collect();
                    let groups = comps8(&off, N64);
                    let n_act = act.iter().filter(|&&b| b).count();
                    let cells: Vec<usize> = (0..nn).filter(|&k| off[k]).collect();
                    let no = cells.len().max(1) as f64;
                    let cont0 = cells.iter().filter(|&&k| pt0[k] == PlateType::Continental).count() as f64 / no;
                    let retyped = cells.iter().filter(|&&k| pt0[k] == PlateType::Oceanic && pt[k] == PlateType::Continental).count() as f64 / no;
                    let band = cells.iter().filter(|&&k| o.band[k]).count() as f64 / no;
                    let mut send: Vec<f64> = cells.iter().map(|&k| o.st.s.data()[k]).collect();
                    // the same on every continental cell re-typed during the run, for scale
                    let all_retyped = (0..nn).filter(|&k| pt0[k] == PlateType::Oceanic && pt[k] == PlateType::Continental).count();
                    let mut img = image::RgbaImage::new(384, 384);
                    for k in 0..nn {
                        let (x, y) = (k % N64, k / N64);
                        let c = if off[k] {
                            [230, 30, 30, 255]
                        } else if act[k] {
                            [250, 200, 40, 255]
                        } else if pt0[k] == PlateType::Oceanic && pt[k] == PlateType::Continental {
                            [230, 120, 200, 255]
                        } else if land64[k] {
                            [200, 185, 140, 255]
                        } else if pt[k] == PlateType::Continental {
                            [150, 170, 160, 255]
                        } else {
                            [40, 70, 140, 255]
                        };
                        for dy in 0..6 {
                            for dx in 0..6 {
                                img.put_pixel((x * 6 + dx) as u32, ((63 - y) * 6 + dy) as u32, image::Rgba(c));
                            }
                        }
                    }
                    img.save(out181().join(format!("equerres_seed{seed}.png"))).unwrap();
                    format!(
                        "   seed {seed} F180: majority-active coast cells {n_act}, of which at sea {} ({:.1} %) in {} groups · of those: continental at t₀ {:.0} %, re-typed during the run {:.0} %, in the band at t₀ {:.0} %, s at the end median {:.3} (p90 {:.3}) · cells re-typed continental over the whole map {all_retyped}\n#EQ {n_act} {} {} {cont0} {retyped} {band}\n",
                        cells.len(), 100.0 * cells.len() as f64 / n_act.max(1) as f64, groups.len(), 100.0 * cont0, 100.0 * retyped, 100.0 * band, q_of(&mut send.clone(), 0.5), q_of(&mut send, 0.9),
                        cells.len(), groups.len()
                    )
                });
            });
        }
    });
    let (mut na, mut no, mut ng) = (0.0, 0.0, 0.0);
    let (mut rt, mut c0) = (0.0, 0.0);
    for seed in seeds179() {
        if let Ok(s) = std::fs::read_to_string(results181().join(format!("eq_seed{seed}.txt"))) {
            eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
            let v = data_line(&s, "EQ").unwrap();
            na += v[0];
            no += v[1];
            ng += v[2];
            rt += v[4] * v[1];
            c0 += v[3] * v[1];
        }
    }
    eprintln!("\n   the équerres: {ng} groups, {no} cells = {:.1} % of the majority-active coast cells · re-typed during the run {:.0} %, continental at t₀ {:.0} %", 100.0 * no / na.max(1.0), 100.0 * rt / no.max(1.0), 100.0 * c0 / no.max(1.0));
    let tiles: Vec<image::RgbaImage> = seeds179().iter().filter_map(|s| image::open(out181().join(format!("equerres_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
    montage(&tiles, 5, 388, 388).save(out181().join("equerres_15.png")).unwrap();
}

/// F181 tier 2 -- 2 048², old: seed 7's unclipped transect and the whole map of the seed with the longest cordillera at the old age.
#[allow(dead_code)]
fn f181_tier2_old() {
    let t0 = Instant::now();
    eprintln!("\n==========  F181 tier 2 . 2 048², old: the transect of seed 7, a whole map  ==========");
    let m = Margin::seaward_oo_collision_f180();
    let cord = |seed: u64| -> f64 { std::fs::read_to_string(results181().join(format!("t1_F180_age3_seed{seed}.txt"))).ok().and_then(|s| data_line(&s, "F")).map(|f| f[1]).unwrap_or(0.0) };
    let best = seeds179().into_iter().max_by(|&x, &y| cord(x).total_cmp(&cord(y))).unwrap();
    let mut todo = vec![(7u64, best == 7)];
    if best != 7 {
        todo.push((best, true));
    }
    for (seed, with_map) in todo {
        let p = sel176(seed).unwrap();
        let s = cached181(&format!("t2_seed{seed}"), || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let act = (0..4096).filter(|&k| is_active_cls(a.major[k])).max_by(|&x, &y| o.w_sum[x].total_cmp(&o.w_sum[y]).then(y.cmp(&x)));
            let pas = (0..4096).filter(|&k| is_passive_cls(a.cls[k])).min_by(|&x, &y| o.src_sum[x].total_cmp(&o.src_sum[y]).then(x.cmp(&y)));
            let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
            let corse_oct = |nc: usize| -> f32 {
                let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
                read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
            };
            let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
            let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
            let kh = KHetero::declared();
            let w = world_oo(seed, p, &m);
            let cover = cover_map_m(seed, p, &m, COVER_M);
            let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
            let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k.clone()) };
            let cfg = base_cfg(seed, rough_a, Some([0, 0]));
            let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
            let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[3].1, &rock, N_PHYS);
            let (g, _, _) = chain_low(&r.zm, amp, Some(ln_k), 0.0);
            if with_map {
                alt_png(&g, 4).save(out181().join(format!("map2048_F180_old_seed{seed}.png"))).unwrap();
            }
            let s64 = w.st.s.data().to_vec();
            let pt = w.st.plate_type.data().to_vec();
            let (mut zs, mut ss, mut desc) = (Vec::new(), Vec::new(), Vec::new());
            for (nm, cell, col) in [("passive", pas, [30u8, 140, 60]), ("active", act, [220, 30, 30])] {
                if let Some(k) = cell {
                    let (za, sb) = transect(&g, &s64, &pt, k);
                    let mx = za.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                    desc.push(format!("{nm} (cell {k}): max {mx:.0} m"));
                    zs.push((za, col));
                    ss.push((sb, col));
                }
            }
            let (zmin, zmax) = transect_png_free(&zs, &ss, &out181().join(format!("transect_F180_old_seed{seed}.png")));
            let peak = g.data.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            format!("   seed {seed} F180 old (2 048²){}: the transect's range {zmin:.0} … {zmax:.0} m · {} · the world's maximum {peak:.0} m\n", if with_map { " (the longest cordillera at the old age; whole map saved)" } else { "" }, desc.join(" · "))
        });
        eprint!("{s}");
    }
    eprintln!("   tier 2: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F181 -- the cordilleras' width at the young age (F180), for the comparison with the intermediate and old ages.
#[allow(dead_code)]
fn f181_young_widths() {
    let rough_a = rough_a177();
    let m = Margin::seaward_oo_collision_f180();
    let mut w = Vec::new();
    for seed in seeds179() {
        let Some(p) = sel176(seed) else { continue };
        let s = cached181(&format!("young_width_seed{seed}"), || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let (r, ww) = physics_age(seed, p, &m, rough_a, 1);
            let (_, form) = formations178(&o, &a, &r, &ww);
            let mut widths = cord_widths(&form, &a);
            format!("   seed {seed} F180 young: cordillera widths (km world) median {:.1} ({:.0} km terrestrial), range {:.1}–{:.1}\n#YW {}\n", q_of(&mut widths.clone(), 0.5), 7.5 * q_of(&mut widths.clone(), 0.5), q_of(&mut widths.clone(), 0.0), q_of(&mut widths.clone(), 1.0), q_of(&mut widths, 0.5))
        });
        eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
        w.push(data_line(&s, "YW").unwrap()[0]);
    }
    eprintln!("   young: cordillera width median of the seeds {:.1} km world ({:.0} terrestrial)", q_of(&mut w.clone(), 0.5), 7.5 * q_of(&mut w, 0.5));
}

/// F181 attribution -- « F176 + the bound alone » (no Andean seed) at the intermediate and old ages, and the split of the fall
/// between F176 and F180 into the bound's part (F176 + bound − F176) and the seed's part (F180 − (F176 + bound)).
#[allow(dead_code)]
fn f181_attribution() {
    let t0 = Instant::now();
    eprintln!("
==========  F181 attribution . F176 + the bound alone, intermediate and old  ==========");
    let rough_a = rough_a177();
    let m = Margin { ds_mass_bound: Some(1.0), ..Margin::seaward_oo_collision() };
    ages_runs(&[("F176B", m)], rough_a);
    ages_summary(&["F176", "F176B", "F180"]);
    for ai in 2..4 {
        let rd = |name: &str, seed: u64| std::fs::read_to_string(results181().join(format!("t1_{name}_age{ai}_seed{seed}.txt"))).ok();
        // per seed: the mountain share, the median altitude, the global W10
        let (mut tot, mut bnd, mut sd) = ([Vec::new(), Vec::new(), Vec::new()], [Vec::new(), Vec::new(), Vec::new()], [Vec::new(), Vec::new(), Vec::new()]);
        let mut cords = [0usize; 3];
        for seed in seeds179() {
            let (Some(a), Some(b), Some(c)) = (rd("F176", seed), rd("F176B", seed), rd("F180", seed)) else { continue };
            let q = |s: &str| -> [f64; 3] {
                let mv = data_line(s, "M").unwrap();
                let wv = data_line(s, "W").unwrap();
                [100.0 * mv[3], mv[5], wv[0]]
            };
            let (qa, qb, qc) = (q(&a), q(&b), q(&c));
            for i in 0..3 {
                tot[i].push(qc[i] - qa[i]);
                bnd[i].push(qb[i] - qa[i]);
                sd[i].push(qc[i] - qb[i]);
            }
            for (j, s) in [&a, &b, &c].iter().enumerate() {
                cords[j] += (data_line(s, "F").unwrap()[0] > 0.0) as usize;
            }
        }
        let names = ["mountain %", "median altitude m", "W10 global"];
        for i in 0..3 {
            let (t, b, s) = (q_of(&mut tot[i].clone(), 0.5), q_of(&mut bnd[i].clone(), 0.5), q_of(&mut sd[i].clone(), 0.5));
            let sumt: f64 = tot[i].iter().sum();
            let sumb: f64 = bnd[i].iter().sum();
            eprintln!(
                "   {} {}: F180 − F176 median {t:+.3} · the bound's part (F176 + bound − F176) median {b:+.3} · the seed's part (F180 − F176 + bound) median {s:+.3} · the bound's share of the summed change {:.0} %",
                AGES[ai].0, names[i], 100.0 * sumb / sumt
            );
        }
        eprintln!("   {} cordilleras: F176 {} · F176 + bound {} · F180 {} (of 15)", AGES[ai].0, cords[0], cords[1], cords[2]);
    }
    eprintln!("   attribution: {:.0} s", t0.elapsed().as_secs_f64());
}

// ==============================  F182 -- the sediment deposition at the physics level  ==============================

use std::cell::RefCell;
use ymir_core::cascade::deposition::{Deposition, FluvialCtx, FluvialHook, FluvialOut, Marine, SedRecord};
use ymir_core::cascade::history::physics_history_dep;
use ymir_core::terrain::flow::{D8_DIST, FlowConfig, compute_flow};

fn out182() -> PathBuf {
    root().join("docs/reports/relief_method/f182_deposition/images")
}
fn results182() -> PathBuf {
    root().join("docs/reports/relief_method/f182_deposition/results")
}
fn cached182(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results182().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results182()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}

/// One F182 world at 256²: the run, the deposition record, the field (norm).
struct Dep182 {
    r: Run,
    sed: Option<SedRecord>,
    /// km³: the incision's eroded volume over the level, the ocean (z < 0) at t₀ and at the end.
    vol: (f64, f64, f64),
}

fn hist_dep(snaps: &[GridF32], cfg: &CascadeConfig, age: Option<f64>, rock: &Rock, dep: Option<Deposition>, hook: Option<&FluvialHook>) -> Dep182 {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: age, sea_low: None };
    let l = physics_history_dep(snaps, 256, cfg, &run, Some(rock), None, dep, hook, &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, Some([0, 0]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    let ckm2 = (CELL * CELL) as f64;
    let eroded = -l.rec.sum_erosion.iter().map(|&v| v as f64).sum::<f64>() * k as f64 / 1000.0 * ckm2;
    let z0 = to_m(&l.z0, k);
    let ocean = |g: &GridF32| g.data.iter().filter(|&&z| z < 0.0).map(|&z| -z as f64).sum::<f64>() / 1000.0 * ckm2;
    let vol = (eroded, ocean(&z0), ocean(&zm));
    Dep182 { r: Run { zm, d5, classes, m, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }, sed: l.rec.sed.clone(), vol }
}

/// The rock and the config of a world (F181's `physics_age`), then the level with `dep`.
fn physics_dep(seed: u64, p: C1Profile, m: &Margin, w: &World, rough_a: f32, age: usize, dep: Option<Deposition>, hook: Option<&FluvialHook>) -> Dep182 {
    let kh = KHetero::declared();
    let cover = cover_map_m(seed, p, m, COVER_M);
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k) };
    let cfg = base_cfg(seed, rough_a, Some([0, 0]));
    hist_dep(&w.snaps, &cfg, AGES[age].1, &rock, dep, hook)
}

/// F182 A5 -- the lakes of a final field (m): 8-connected components (torus) of land cells with `filled > z`; per class (dust
/// 1–2 cells, 3–15, 16–63, large ≥ 64) the count and the area (km²); the mask.
fn lakes182(zm: &GridF32) -> ([usize; 4], [f64; 4], Vec<bool>) {
    let flow = compute_flow(zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let mask: Vec<bool> = (0..65536).map(|k| zm.data[k] > 0.0 && flow.filled.data[k] > zm.data[k]).collect();
    let (mut n, mut a) = ([0usize; 4], [0f64; 4]);
    for c in comps8(&mask, 256) {
        let cls = match c.len() {
            1..=2 => 0,
            3..=15 => 1,
            16..=63 => 2,
            _ => 3,
        };
        n[cls] += 1;
        a[cls] += c.len() as f64 * (CELL * CELL) as f64;
    }
    (n, a, mask)
}

/// F182 A5 -- the shelf (the declaration's 0.2): from each coastal cell of the main ocean, the sea floor's D8 steepest descent
/// (≤ 128 cells); the break at the start of the link with the largest increase of gradient; the width (km) and the break's
/// depth (m). Per coast class from A1 (the 64² parent and its 8 neighbours, the majority; 0 none, 1 active, 2 passive).
fn shelf182(zm: &GridF32, a: &A1) -> Vec<(u8, f64, f64, usize)> {
    let n = 256;
    let oc = ocean_mask_t(zm, Some([0, 0]));
    let mut out = Vec::new();
    for k in 0..n * n {
        if !oc[k] {
            continue;
        }
        let (x, y) = ((k % n) as i32, (k / n) as i32);
        let coastal = [(1, 0), (-1, 0), (0, 1), (0, -1)].iter().any(|&(dx, dy)| {
            let j = (y + dy).rem_euclid(n as i32) as usize * n + (x + dx).rem_euclid(n as i32) as usize;
            zm.data[j] > 0.0
        });
        if !coastal {
            continue;
        }
        let mut path = vec![k];
        let mut lens = Vec::new();
        while path.len() < 128 {
            let c = *path.last().unwrap();
            let (cx, cy) = ((c % n) as i32, (c / n) as i32);
            let mut best: Option<(usize, f32, f32)> = None;
            for d in 0..8 {
                let j = (cy + D8_DY[d]).rem_euclid(n as i32) as usize * n + (cx + D8_DX[d]).rem_euclid(n as i32) as usize;
                if !oc[j] {
                    continue;
                }
                let s = (zm.data[c] - zm.data[j]) / D8_DIST[d];
                if s > 0.0 && best.is_none_or(|b| s > b.1) {
                    best = Some((j, s, D8_DIST[d]));
                }
            }
            let Some((j, _, l)) = best else { break };
            path.push(j);
            lens.push(l as f64 * CELL as f64);
        }
        if path.len() < 3 {
            continue;
        }
        let g: Vec<f64> = (0..lens.len()).map(|i| (zm.data[path[i]] - zm.data[path[i + 1]]) as f64 / lens[i]).collect();
        let mut bi = 1;
        for i in 1..g.len() {
            if g[i] - g[i - 1] > g[bi] - g[bi - 1] {
                bi = i;
            }
        }
        let width: f64 = lens[..bi].iter().sum();
        let c64 = ((k / n) / 4) * 64 + (k % n) / 4;
        let (cx, cy) = ((c64 % 64) as i32, (c64 / 64) as i32);
        let (mut act, mut pas) = (0, 0);
        for dy in -1..=1 {
            for dx in -1..=1 {
                let j = (cy + dy).rem_euclid(64) as usize * 64 + (cx + dx).rem_euclid(64) as usize;
                act += is_active_cls(a.major[j]) as usize;
                pas += is_passive_cls(a.cls[j]) as usize;
            }
        }
        let cls = if act + pas == 0 { 0 } else if act > pas { 1 } else { 2 };
        out.push((cls, width, -(zm.data[path[bi]] as f64), k));
    }
    out
}

/// The deltas: components (8, torus) of land with the deposition that is sea without it and touches the main ocean.
fn deltas182(with: &GridF32, without: &GridF32) -> (usize, f64, f64) {
    let oc = ocean_mask_t(with, Some([0, 0]));
    let gain: Vec<bool> = (0..65536).map(|k| with.data[k] > 0.0 && without.data[k] <= 0.0).collect();
    let lost = (0..65536).filter(|&k| with.data[k] <= 0.0 && without.data[k] > 0.0).count();
    let mut cnt = 0;
    let mut area = 0.0;
    for c in comps8(&gain, 256) {
        if c.iter().any(|&k| nb8_n(k, 256).iter().any(|&j| oc[j])) {
            cnt += 1;
            area += c.len() as f64 * (CELL * CELL) as f64;
        }
    }
    (cnt, area, lost as f64 * (CELL * CELL) as f64)
}

/// A map of a positive field (m) on the land/sea relief: white → dark brown, log-scaled up to `top` m.
fn thick_png(zm: &GridF32, t: &[f32], top: f32) -> image::RgbaImage {
    let mut img = alt_png_plain(zm);
    for k in 0..65536 {
        if t[k] > 0.5 {
            let f = ((1.0 + t[k]).ln() / (1.0 + top).ln()).clamp(0.0, 1.0);
            let (x, y) = (k % 256, k / 256);
            img.put_pixel(x as u32, (255 - y) as u32, image::Rgba([(250.0 - 150.0 * f) as u8, (230.0 - 180.0 * f) as u8, (160.0 - 150.0 * f) as u8, 255]));
        }
    }
    img
}

/// The lakes in red on the relief.
fn lakes_png(zm: &GridF32, mask: &[bool]) -> image::RgbaImage {
    let mut img = alt_png_plain(zm);
    for k in 0..65536 {
        if mask[k] {
            let (x, y) = (k % 256, k / 256);
            img.put_pixel(x as u32, (255 - y) as u32, image::Rgba([220, 20, 20, 255]));
        }
    }
    img
}

/// The shelf breaks on the relief (yellow: passive, red: active, white: none).
fn shelf_png(zm: &GridF32, sh: &[(u8, f64, f64, usize)]) -> image::RgbaImage {
    let mut img = alt_png_plain(zm);
    for &(c, w, _, k) in sh {
        // the break cell is not kept; the coast cell is coloured by its width (dark = narrow)
        let f = (w / 60.0).clamp(0.0, 1.0);
        let px = match c {
            1 => [(120.0 + 135.0 * f) as u8, 20, 20, 255],
            2 => [(120.0 + 135.0 * f) as u8, (120.0 + 135.0 * f) as u8, 20, 255],
            _ => [(120.0 + 135.0 * f) as u8, (120.0 + 135.0 * f) as u8, (120.0 + 135.0 * f) as u8, 255],
        };
        let (x, y) = (k % 256, k / 256);
        img.put_pixel(x as u32, (255 - y) as u32, image::Rgba(px));
    }
    img
}

/// The measures of one world (with or without the deposition) as data lines.
fn measure182(tag: &str, d: &Dep182, o: &Obs, a: &A1, w: &World) -> (String, Vec<bool>) {
    let r = &d.r;
    let (ln, la, mask) = lakes182(&r.zm);
    let (line, _) = formations178(o, a, r, w);
    let land: Vec<bool> = r.zm.data.iter().map(|&z| z > 0.0).collect();
    let kc = coast_torus(&land, 256, 256, 4);
    let circ = circumnavigation(&land, 256);
    let sh = shelf182(&r.zm, a);
    let med = |c: u8| -> f64 {
        let mut v: Vec<f64> = sh.iter().filter(|s| c == 9 || s.0 == c).map(|s| s.1).collect();
        if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) }
    };
    let mut dep: Vec<f64> = sh.iter().filter(|s| s.0 == 2).map(|s| s.2).collect();
    let dmed = if dep.is_empty() { f64::NAN } else { q_of(&mut dep, 0.5) };
    let dec = &r.m.deciles;
    let mut s = format!(
        "      {tag}: lakes (dust 1–2 / 3–15 / 16–63 / large ≥ 64 cells) {:?}, area km² {:?} · below-sea inland {} cells · shelf median km: passive {:.1}, active {:.1}, all {:.1} ({} paths), the passive break's depth {dmed:.0} m · coast dimension {:.3}, length ratio {:.3} · circumnavigable {} · deciles {:?}\n{line}",
        ln, la.map(|v| v.round()), r.d5.iter().filter(|&&b| b).count(), med(2), med(1), med(9), sh.len(), kc.dimension, kc.length_ratio, !circ.wraps, dec.iter().map(|v| v.round()).collect::<Vec<_>>()
    );
    s += &format!("#L {} {} {} {} {} {} {} {}\n", ln[0], ln[1], ln[2], ln[3], la[0], la[1], la[2], la[3]);
    s += &format!("#S {} {} {} {} {}\n", med(2), med(1), med(9), dmed, sh.len());
    s += &format!("#C {} {} {} {}\n", kc.dimension, kc.length_ratio, !circ.wraps as u8, (land.iter().filter(|&&l| l).count()));
    s += &format!("#D {}\n", dec.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" "));
    if let Some(sr) = &d.sed {
        let b = sr.budget;
        let tot = b[2] + b[3] + b[4] + b[5] + b[6];
        let mut thick: Vec<f64> = sr.thickness_m.iter().zip(&land).filter(|(_, l)| **l).map(|(&t, _)| t as f64).collect();
        let tmed = q_of(&mut thick, 0.5);
        let tmax = sr.thickness_m.iter().copied().fold(0.0f32, f32::max);
        let it_max = sr.iterations.iter().copied().max().unwrap_or(0);
        s += &format!(
            "      the sediment: eroded {:.4e} + hillslope to the sea {:.4e} (norm × cells) → land {:.1} % · lakes {:.1} % · stranded {:.2} % · sea {:.1} % · exported {:.2} % · the mass control's worst step {:.2e} (f32 field {:.2e}) · thickness on land median {tmed:.1} m, max {tmax:.0} m · marine sub-steps {} in {:.1} s · iterations max {it_max}\n#B {} {} {} {} {} {} {} {} {} {} {}\n",
            b[0], b[1], 100.0 * b[2] / tot, 100.0 * b[3] / tot, 100.0 * b[4] / tot, 100.0 * b[5] / tot, 100.0 * b[6] / tot, sr.max_rel_residual, sr.max_field_residual, sr.marine_substeps, sr.secs_marine,
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], sr.max_rel_residual, sr.max_field_residual, sr.marine_substeps, sr.secs_marine
        );
    }
    (s, mask)
}

/// F182 -- the cost per world at the three ages (one panel seed), off and on: announced before the tiers.
#[allow(dead_code)]
fn f182_cost() {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let seed = PANEL[0];
    let p = sel176(seed).unwrap();
    let w = world_oo(seed, p, &m);
    eprintln!("\n==========  F182 the cost per world (seed {seed}, 256², F181's default)  ==========\n   C1: {}", fmt_cost(w.cost));
    for ai in 1..4 {
        let s = cached182(&format!("cost_age{ai}_seed{seed}"), || {
            let off = physics_dep(seed, p, &m, &w, rough_a, ai, None, None);
            let on = physics_dep(seed, p, &m, &w, rough_a, ai, Some(Deposition::declared()), None);
            let sr = on.sed.as_ref().unwrap();
            format!(
                "   {}: without {} · with {} (marine {} sub-steps, {:.1} s wall) · the mass control's worst step {:.2e}\n#K {} {} {} {} {} {}\n",
                AGES[ai].0, fmt_cost(off.r.cost), fmt_cost(on.r.cost), sr.marine_substeps, sr.secs_marine, sr.max_rel_residual,
                off.r.cost.0, off.r.cost.1, on.r.cost.0, on.r.cost.1, sr.marine_substeps, sr.secs_marine
            )
        });
        eprint!("{s}");
    }
}

/// F182 tier 1 -- 256², the 15 seeds, the intermediate age: F181's default without and with the deposition.
#[allow(dead_code)]
fn f182_tier1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out182()).unwrap();
    eprintln!("\n==========  F182 tier 1 . 256², intermediate: without / with the deposition  ==========");
    tier182(2, &seeds179(), "t1");
    summary182("t1", 2, &seeds179());
    for pre in ["lakes_off", "lakes_on", "thickness", "shelf_on"] {
        let tiles: Vec<image::RgbaImage> = seeds179().iter().filter_map(|s| image::open(out182().join(format!("{pre}_age2_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
        if !tiles.is_empty() {
            montage(&tiles, 5, 260, 260).save(out182().join(format!("{pre}_15_age2.png"))).unwrap();
        }
    }
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F182 tier 2 -- the young and old ages, the panel.
#[allow(dead_code)]
fn f182_tier2() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out182()).unwrap();
    eprintln!("\n==========  F182 tier 2 . 256², young and old, the panel: without / with the deposition  ==========");
    for ai in [1, 3] {
        tier182(ai, &PANEL, "t2");
        summary182("t2", ai, &PANEL);
    }
    eprintln!("   tier 2: {:.0} s", t0.elapsed().as_secs_f64());
}

/// The per-seed runs of a tier (resumable).
fn tier182(ai: usize, seeds: &[u64], tier: &str) {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    for &seed in seeds {
        let Some(p) = sel176(seed) else { continue };
        let name = format!("{tier}_age{ai}_seed{seed}");
        if results182().join(format!("{name}.txt")).exists() {
            eprint!("{}", std::fs::read_to_string(results182().join(format!("{name}.txt"))).unwrap());
            continue;
        }
        let s = cached182(&name, || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let w = world_oo(seed, p, &m);
            let off = physics_dep(seed, p, &m, &w, rough_a, ai, None, None);
            let on = physics_dep(seed, p, &m, &w, rough_a, ai, Some(Deposition::declared()), None);
            let (s_off, mk_off) = measure182("without", &off, &o, &a, &w);
            let (s_on, mk_on) = measure182("with", &on, &o, &a, &w);
            let (nd, ad, lost) = deltas182(&on.r.zm, &off.r.zm);
            lakes_png(&off.r.zm, &mk_off).save(out182().join(format!("lakes_off_age{ai}_seed{seed}.png"))).unwrap();
            lakes_png(&on.r.zm, &mk_on).save(out182().join(format!("lakes_on_age{ai}_seed{seed}.png"))).unwrap();
            let sr = on.sed.as_ref().unwrap();
            thick_png(&on.r.zm, &sr.thickness_m, 2000.0).save(out182().join(format!("thickness_age{ai}_seed{seed}.png"))).unwrap();
            shelf_png(&on.r.zm, &shelf182(&on.r.zm, &a)).save(out182().join(format!("shelf_on_age{ai}_seed{seed}.png"))).unwrap();
            shelf_png(&off.r.zm, &shelf182(&off.r.zm, &a)).save(out182().join(format!("shelf_off_age{ai}_seed{seed}.png"))).unwrap();
            format!(
                "   seed {seed} {}: physics without {} · with {}\n{s_off}{s_on}      deltas (land gained at the coast) {nd}, {ad:.0} km² · land lost to the sea {lost:.0} km²\n#G {nd} {ad} {lost} {} {}\n",
                AGES[ai].0, fmt_cost(off.r.cost), fmt_cost(on.r.cost), off.r.cost.1, on.r.cost.1
            )
        });
        eprint!("{s}");
    }
}

/// The data lines of a world's file: the first block is « without », the second « with ».
fn lines182(s: &str, key: &str) -> Vec<Vec<f64>> {
    s.lines().filter(|l| l.starts_with(&format!("#{key} "))).map(|l| l[key.len() + 2..].split_whitespace().map(|v| v.parse().unwrap_or(f64::NAN)).collect()).collect()
}

/// The tier's summary: the medians over the seeds, without → with, and the predictions' counts.
fn summary182(tier: &str, ai: usize, seeds: &[u64]) {
    let mut rows: Vec<(u64, [Vec<Vec<f64>>; 7], Vec<f64>, Vec<f64>)> = Vec::new();
    for &seed in seeds {
        let Ok(s) = std::fs::read_to_string(results182().join(format!("{tier}_age{ai}_seed{seed}.txt"))) else { continue };
        let ks = ["L", "S", "C", "D", "F", "W", "M"].map(|k| lines182(&s, k));
        let b = lines182(&s, "B").pop().unwrap_or_default();
        let g = lines182(&s, "G").pop().unwrap_or_default();
        rows.push((seed, ks, b, g));
    }
    let n = rows.len();
    let med = |f: &dyn Fn(&(u64, [Vec<Vec<f64>>; 7], Vec<f64>, Vec<f64>), usize) -> f64, i: usize| -> f64 {
        let mut v: Vec<f64> = rows.iter().map(|r| f(r, i)).filter(|x| x.is_finite()).collect();
        if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) }
    };
    let get = |key: usize, col: usize| move |r: &(u64, [Vec<Vec<f64>>; 7], Vec<f64>, Vec<f64>), i: usize| r.1[key].get(i).and_then(|v| v.get(col)).copied().unwrap_or(f64::NAN);
    let pair = |key: usize, col: usize| (med(&get(key, col), 0), med(&get(key, col), 1));
    let large = pair(0, 3);
    let dust = pair(0, 0);
    let la = (med(&|r, i| r.1[0].get(i).map_or(f64::NAN, |v| v[4] + v[5] + v[6] + v[7]), 0), med(&|r, i| r.1[0].get(i).map_or(f64::NAN, |v| v[4] + v[5] + v[6] + v[7]), 1));
    let large_drop: Vec<f64> = rows.iter().filter_map(|r| {
        let (a, b) = (r.1[0].first()?[3], r.1[0].get(1)?[3]);
        (a > 0.0).then(|| 100.0 * (b - a) / a)
    }).collect();
    let sp = pair(1, 0);
    let sa = pair(1, 1);
    let dim = pair(2, 0);
    let lr = pair(2, 1);
    let plain = pair(6, 0);
    let mtn = pair(6, 3);
    let medalt = pair(6, 5);
    let mut ddim = Vec::new();
    let mut dw10 = Vec::new();
    let (mut v_ok, mut circ_off, mut circ_on, mut in_band, mut in_harris) = (0, 0, 0, 0, 0);
    let mut dplain = Vec::new();
    let mut dmed = Vec::new();
    for r in &rows {
        let (c0, c1) = (&r.1[2], &r.1[2]);
        if let (Some(a), Some(b)) = (c0.first(), c1.get(1)) {
            ddim.push((b[0] - a[0]).abs());
            circ_off += (a[2] > 0.0) as usize;
            circ_on += (b[2] > 0.0) as usize;
        }
        if let (Some(a), Some(b)) = (r.1[5].first(), r.1[5].get(1)) {
            dw10.push((b[0] - a[0]).abs());
        }
        if let Some(f) = r.1[4].get(1) {
            v_ok += (f[8] == 0.0) as usize;
        }
        if let Some(s) = r.1[1].get(1) {
            in_band += (s[0] >= 4.65 && s[0] <= 16.4) as usize;
            in_harris += (s[0] >= 34.9 && s[0] <= 123.0) as usize;
        }
        if let (Some(a), Some(b)) = (r.1[6].first(), r.1[6].get(1)) {
            dplain.push(100.0 * (b[0] - a[0]));
            dmed.push(b[5] - a[5]);
        }
    }
    let shares: Vec<[f64; 4]> = rows.iter().filter(|r| r.2.len() >= 9).map(|r| {
        let b = &r.2;
        let t = b[2] + b[3] + b[4] + b[5] + b[6];
        [100.0 * (b[2] + b[4]) / t, 100.0 * b[3] / t, 100.0 * b[5] / t, b[7]]
    }).collect();
    let sh_med = |i: usize| -> f64 { let mut v: Vec<f64> = shares.iter().map(|s| s[i]).collect(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
    let worst = shares.iter().map(|s| s[3]).fold(0.0, f64::max);
    let deltas: Vec<f64> = rows.iter().filter_map(|r| r.3.get(1).copied()).collect();
    let nd: Vec<f64> = rows.iter().filter_map(|r| r.3.first().copied()).collect();
    let land_on: Vec<f64> = rows.iter().filter_map(|r| r.1[6].get(1).map(|v| v[4])).collect();
    let dpct: Vec<f64> = deltas.iter().zip(&land_on).map(|(d, l)| 100.0 * d / l).collect();
    let q = |v: &[f64], p: f64| -> f64 { let mut v = v.to_vec(); if v.is_empty() { f64::NAN } else { q_of(&mut v, p) } };
    eprintln!(
        "\n   F182 {tier} {} ({n} seeds), the medians without → with:\n      large lakes {:.0} → {:.0} (per-seed change median {:.0} %, range {:.0}…{:.0} %) · dust {:.0} → {:.0} · lake area {:.0} → {:.0} km²\n      shelf passive {:.1} → {:.1} km · active {:.1} → {:.1} km · in Harris ÷ 7.5 (4.65–16.4 km) {in_band} / {n} · in Harris unconverted (34.9–123.0 km) {in_harris} / {n}\n      plains {:.1} → {:.1} % (per-seed change median {:+.1} points) · mountain {:.1} → {:.1} % · median altitude {:.0} → {:.0} m (change median {:+.0} m)\n      coast dimension {:.3} → {:.3} (|Δ| max {:.3}) · length ratio {:.3} → {:.3} · deltas per seed median {:.0}, area median {:.2} % of the land\n      the sediment (median shares): land {:.1} % · lakes {:.1} % · sea {:.1} % · the mass control's worst step over the seeds {worst:.2e}\n      (v) holds on {v_ok} / {n} · W10 |Δ| max {:.3} · circumnavigable without {circ_off} / {n}, with {circ_on} / {n}\n",
        AGES[ai].0, large.0, large.1, q(&large_drop, 0.5), q(&large_drop, 0.0), q(&large_drop, 1.0), dust.0, dust.1, la.0, la.1,
        sp.0, sp.1, sa.0, sa.1, 100.0 * plain.0, 100.0 * plain.1, q(&dplain, 0.5), 100.0 * mtn.0, 100.0 * mtn.1, medalt.0, medalt.1, q(&dmed, 0.5),
        dim.0, dim.1, q(&ddim, 1.0), lr.0, lr.1, q(&nd, 0.5), q(&dpct, 0.5), sh_med(0), sh_med(1), sh_med(2), q(&dw10, 1.0)
    );
}

/// ADR Finding 183 (the author's decision) -- **abandoned as a reference**: its convergence is proven only for G / p̃ ≤ 1 (Yuan
/// et al. 2019 JGR p. 8) and G = 1.55; at F183 it never converged in 50 iterates (300 / 300 steps). The code is kept as a record.
/// F182 A3 -- Yuan's iterative coupling, the bench's reference (JGR eq. 19–21): each iterate the cascade's incision of the
/// field plus the previous iterate's deposit (eq. 20's form), the deposit from that iterate's erosion by the cascade's pass,
/// routed on the iterate's surface without its lake fill; stopped at eq. 21 (max |Δh| < 10⁻³ U Δt), at most 50 iterates;
/// mass-exact output z = before − e^K + D^K. Logs per step: (step, iterates, converged, the gap to the first iterate (the
/// splitting): Σ|D^K − D¹| / ΣD^K and the mean |z^K − z¹| on land, norm).
fn yuan_ref(c: &FluvialCtx, log: &RefCell<Vec<(usize, usize, bool, f64, f64)>>) -> FluvialOut {
    let tol = 1e-3 * c.u_dt;
    let n = c.before.data.len();
    let mut d_prev = vec![0f64; n];
    let mut lake_prev = vec![0f64; n];
    let mut h_prev: Option<GridF32> = None;
    let mut first: Option<(Vec<f64>, Vec<f64>)> = None;
    let mut it = 0;
    loop {
        it += 1;
        let mut surf = c.before.clone();
        for k in 0..n {
            if d_prev[k] != 0.0 {
                surf.data[k] = (surf.data[k] as f64 + d_prev[k]) as f32;
            }
        }
        let h = (c.incise)(&surf);
        let e: Vec<f64> = (0..n).map(|k| (surf.data[k] as f64 - h.data[k] as f64).max(0.0)).collect();
        let mut bare = h.clone();
        for k in 0..n {
            if lake_prev[k] != 0.0 {
                bare.data[k] = (bare.data[k] as f64 - lake_prev[k]) as f32;
            }
        }
        let p = (c.pass)(&bare, &e);
        let dmax = h_prev.as_ref().map_or(f64::INFINITY, |hp| (0..n).map(|k| (h.data[k] - hp.data[k]).abs() as f64).fold(0.0, f64::max));
        let dk: Vec<f64> = (0..n).map(|k| p.deposit(k)).collect();
        let zk: Vec<f64> = (0..n).map(|k| c.before.data[k] as f64 - e[k] + dk[k]).collect();
        if first.is_none() {
            first = Some((dk.clone(), zk.clone()));
        }
        let conv = dmax < tol;
        if conv || it >= 50 {
            let (d1, z1) = first.as_ref().unwrap();
            let sd: f64 = dk.iter().sum();
            let gap_d = if sd > 0.0 { dk.iter().zip(d1).map(|(a, b)| (a - b).abs()).sum::<f64>() / sd } else { 0.0 };
            let lnd: Vec<usize> = (0..n).filter(|&k| c.before.data[k] > 0.5).collect();
            let gap_z = lnd.iter().map(|&k| (zk[k] - z1[k]).abs()).sum::<f64>() / lnd.len().max(1) as f64;
            log.borrow_mut().push((c.step, it, conv, gap_d, gap_z));
            let mut z = c.before.clone();
            for k in 0..n {
                z.data[k] = zk[k] as f32;
            }
            return FluvialOut { z, e, pass: p, iterations: it };
        }
        lake_prev = p.lake.clone();
        d_prev = dk;
        h_prev = Some(h);
    }
}

/// F182 A3 -- the reference against the splitting, on 2 panel seeds at the intermediate age: the deposit maps (correlation, the
/// mass difference), the masses, the lake area, the median altitude; P7: the per-step gap against the step's duration.
#[allow(dead_code)]
fn f182_reference() {
    let t0 = Instant::now();
    eprintln!("\n==========  F182 A3 . the explicit splitting against Yuan's iterative coupling (bench only), 2 panel seeds, intermediate  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[2].1.unwrap(), N_PHYS);
    for seed in [PANEL[0], PANEL[1]] {
        let p = sel176(seed).unwrap();
        let s = cached182(&format!("ref_seed{seed}"), || {
            let w = world_oo(seed, p, &m);
            let split = physics_dep(seed, p, &m, &w, rough_a, 2, Some(Deposition::declared()), None);
            let log = RefCell::new(Vec::new());
            let hook = |c: &FluvialCtx| yuan_ref(c, &log);
            let refr = physics_dep(seed, p, &m, &w, rough_a, 2, Some(Deposition::declared()), Some(&hook));
            let (ss, rs) = (split.sed.as_ref().unwrap(), refr.sed.as_ref().unwrap());
            let landdep = |r: &SedRecord| -> Vec<f64> { (0..65536).map(|k| (r.fluvial_m[k] + r.lake_m[k]) as f64).collect() };
            let (a, b) = (landdep(ss), landdep(rs));
            let corr = pearson(&a.iter().map(|&v| v as f32).collect::<Vec<_>>(), &b.iter().map(|&v| v as f32).collect::<Vec<_>>());
            let (sa, sb) = (a.iter().sum::<f64>(), b.iter().sum::<f64>());
            let (la, lb) = (lakes182(&split.r.zm), lakes182(&refr.r.zm));
            let area = |l: &([usize; 4], [f64; 4], Vec<bool>)| l.1.iter().sum::<f64>();
            let lg = log.borrow();
            let its: Vec<f64> = lg.iter().map(|l| l.1 as f64).collect();
            let nconv = lg.iter().filter(|l| !l.2).count();
            // P7: the gap at the longest and the shortest steps, and its rank correlation with the duration
            let imax = (0..N_PHYS).max_by(|&x, &y| dur[x].total_cmp(&dur[y])).unwrap();
            let imin = (0..N_PHYS).min_by(|&x, &y| dur[x].total_cmp(&dur[y])).unwrap();
            let at = |i: usize| lg.iter().find(|l| l.0 == i).map(|l| (l.3, l.4 * n2m() as f64, l.1)).unwrap_or((f64::NAN, f64::NAN, 0));
            let dv: Vec<f64> = lg.iter().map(|l| dur[l.0.min(N_PHYS - 1)]).collect();
            let gv: Vec<f64> = lg.iter().map(|l| l.3).collect();
            let zv: Vec<f64> = lg.iter().map(|l| l.4).collect();
            let rho_d = spearman(&dv.iter().map(|&v| v as f32).collect::<Vec<_>>(), &gv.iter().map(|&v| v as f32).collect::<Vec<_>>());
            let rho_z = spearman(&dv.iter().map(|&v| v as f32).collect::<Vec<_>>(), &zv.iter().map(|&v| v as f32).collect::<Vec<_>>());
            let (gl, gs) = (at(imax), at(imin));
            let mut tl = String::new();
            for l in lg.iter() {
                tl += &format!("{} {} {} {} {} {}\n", l.0, dur[l.0.min(N_PHYS - 1)], l.1, l.2 as u8, l.3, l.4 * n2m() as f64);
            }
            std::fs::create_dir_all(results182()).unwrap();
            std::fs::write(results182().join(format!("ref_steps_seed{seed}.txt")), tl).unwrap();
            let bs = |r: &SedRecord| { let b = r.budget; let t = b[2] + b[3] + b[4] + b[5] + b[6]; (100.0 * (b[2] + b[4]) / t, 100.0 * b[3] / t, 100.0 * b[5] / t, t) };
            let (bsa, bsb) = (bs(ss), bs(rs));
            format!(
                "   seed {seed}: splitting {} · reference {} (iterates per step median {:.0}, max {:.0}; not converged in 50 on {nconv} / {} steps)\n      the land deposit maps (fluvial + lakes, m): r = {corr:.3} · mass {sa:.4e} → {sb:.4e} m·cells ({:+.1} %) · total deposit {:.4e} → {:.4e} norm × cells ({:+.1} %) · shares land / lakes / sea {:.1} / {:.1} / {:.1} % → {:.1} / {:.1} / {:.1} %\n      lake area {:.0} → {:.0} km² · median altitude {:.1} → {:.1} m · the reference's mass control {:.2e}\n      P7: the deposit gap Σ|D_Yuan − D_split| / ΣD_Yuan at the longest step ({:.2} Myr, step {imax}) {:.3} (mean |Δz| {:.2} m, {} iterates), at the shortest ({:.1} kyr, step {imin}) {:.3} (mean |Δz| {:.3} m, {} iterates); Spearman (duration, gap) {rho_d:.2}, (duration, |Δz|) {rho_z:.2}\n#R {corr} {sa} {sb} {} {} {} {} {} {} {} {}\n",
                fmt_cost(split.r.cost), fmt_cost(refr.r.cost), q_of(&mut its.clone(), 0.5), q_of(&mut its.clone(), 1.0), lg.len(),
                100.0 * (sb - sa) / sa, bsa.3, bsb.3, 100.0 * (bsb.3 - bsa.3) / bsa.3, bsa.0, bsa.1, bsa.2, bsb.0, bsb.1, bsb.2,
                area(&la), area(&lb), split.r.m.deciles[4], refr.r.m.deciles[4], rs.max_rel_residual,
                dur[imax] / 1e6, gl.0, gl.1, gl.2, dur[imin] / 1e3, gs.0, gs.1, gs.2,
                area(&la), area(&lb), split.r.m.deciles[4], refr.r.m.deciles[4], gl.0, gs.0, rho_d, rho_z
            )
        });
        eprint!("{s}");
    }
    eprintln!("   reference: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F182 -- the sensitivity (the panel, intermediate): G = 0 and 3.1 (Guerit's bounds), K_M = 300 and 500 (the Ogooué range).
#[allow(dead_code)]
fn f182_sensitivity() {
    let t0 = Instant::now();
    eprintln!("\n==========  F182 the sensitivity . the panel, intermediate: G 0 / 3.1, K_M 300 / 500  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let d = Deposition::declared();
    let variants = [
        ("G0", Deposition { g: 0.0, ..d }),
        ("G3.1", Deposition { g: 3.1, ..d }),
        ("K300", Deposition { marine: Some(Marine { k_m2_yr: 300.0, ..Marine::declared() }), ..d }),
        ("K500", Deposition { marine: Some(Marine { k_m2_yr: 500.0, ..Marine::declared() }), ..d }),
    ];
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let s = cached182(&format!("sens_seed{seed}"), || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let w = world_oo(seed, p, &m);
            let mut s = format!("   seed {seed}:\n");
            for (nm, dv) in variants {
                let r = physics_dep(seed, p, &m, &w, rough_a, 2, Some(dv), None);
                let (line, _) = measure182(nm, &r, &o, &a, &w);
                s += &line.replace("#L", &format!("#L{nm}")).replace("#S", &format!("#S{nm}")).replace("#C", &format!("#C{nm}")).replace("#M", &format!("#M{nm}")).replace("#B", &format!("#B{nm}"));
            }
            s
        });
        eprint!("{s}");
    }
    // the medians per variant
    for (nm, _) in variants {
        let (mut large, mut sp, mut plain, mut med, mut sea) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for seed in PANEL {
            let Ok(s) = std::fs::read_to_string(results182().join(format!("sens_seed{seed}.txt"))) else { continue };
            if let Some(l) = lines182(&s, &format!("L{nm}")).pop() { large.push(l[3]); }
            if let Some(v) = lines182(&s, &format!("S{nm}")).pop() { if v[0].is_finite() { sp.push(v[0]); } }
            if let Some(v) = lines182(&s, &format!("M{nm}")).pop() { plain.push(100.0 * v[0]); med.push(v[5]); }
            if let Some(b) = lines182(&s, &format!("B{nm}")).pop() { let t = b[2] + b[3] + b[4] + b[5] + b[6]; sea.push(100.0 * b[5] / t); }
        }
        let q = |v: &mut Vec<f64>| if v.is_empty() { f64::NAN } else { q_of(v, 0.5) };
        eprintln!("   {nm}: medians over the panel · large lakes {:.0} · shelf passive {:.1} km · plains {:.1} % · median altitude {:.0} m · sea share {:.1} %", q(&mut large), q(&mut sp), q(&mut plain), q(&mut med), q(&mut sea));
    }
    eprintln!("   sensitivity: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F182 tier 3 -- 2 048², one seed at the intermediate age, without and with the deposition: the whole maps, a crop on the largest
/// delta (or the widest land deposit), a transect across the shelf of a passive coast.
#[allow(dead_code)]
fn f182_tier3() {
    let t0 = Instant::now();
    eprintln!("\n==========  F182 tier 3 . 2 048², intermediate, one seed: without / with the deposition  ==========");
    let m = Margin::continent_default();
    // the seed: the panel seed of tier 1 with the largest delta area
    let delta = |sd: u64| -> f64 { std::fs::read_to_string(results182().join(format!("t1_age2_seed{sd}.txt"))).ok().and_then(|s| lines182(&s, "G").pop()).map(|g| g[1]).unwrap_or(-1.0) };
    let seed = PANEL.iter().copied().max_by(|&x, &y| delta(x).total_cmp(&delta(y)).then(y.cmp(&x))).unwrap();
    let p = sel176(seed).unwrap();
    let s = cached182(&format!("t3_seed{seed}"), || {
        let (o, _) = observe179(seed, p, &m);
        let info = step_info(&o.st, &o.kin, o.ds_max, true);
        let a = a1(&info, &o.kin, &o.tracker);
        let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
        let corse_oct = |nc: usize| -> f32 {
            let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
            read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
        };
        let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
        let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
        let kh = KHetero::declared();
        let w = world_oo(seed, p, &m);
        let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
        let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
        let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
        let on = physics_dep(seed, p, &m, &w, rough_a, 2, Some(Deposition::declared()), None);
        let (g0, _, _) = chain_low(&off.r.zm, amp.clone(), Some(ln_k.clone()), 0.0);
        let (g1, _, _) = chain_low(&on.r.zm, amp, Some(ln_k), 0.0);
        alt_png(&g0, 4).save(out182().join(format!("map2048_off_seed{seed}.png"))).unwrap();
        alt_png(&g1, 4).save(out182().join(format!("map2048_on_seed{seed}.png"))).unwrap();
        // the crop: centred on the 256² cell of the largest delta component, else of the thickest land deposit
        let gain: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && off.r.zm.data[k] <= 0.0).collect();
        let sr = on.sed.as_ref().unwrap();
        let centre = comps8(&gain, 256).into_iter().max_by_key(|c| c.len()).map(|c| c[c.len() / 2]).unwrap_or_else(|| {
            (0..65536).filter(|&k| on.r.zm.data[k] > 0.0).max_by(|&x, &y| sr.thickness_m[x].total_cmp(&sr.thickness_m[y])).unwrap()
        });
        let (cx, cy) = ((centre % 256) * 8, (centre / 256) * 8);
        for (nm, g) in [("off", &g0), ("on", &g1)] {
            let half = 192;
            let mut crop = GridF32::new(2 * half, 2 * half, 0.0);
            for y in 0..2 * half {
                for x in 0..2 * half {
                    let sx = (cx + 2048 - half + x) % 2048;
                    let sy = (cy + 2048 - half + y) % 2048;
                    crop.data[y * 2 * half + x] = g.data[sy * 2048 + sx];
                }
            }
            alt_png(&crop, 1).save(out182().join(format!("crop2048_{nm}_seed{seed}.png"))).unwrap();
        }
        // the transect: from the passive coast cell with the widest shelf (with the deposition), along its seaward path direction
        let sh = shelf182(&on.r.zm, &a);
        let pick = sh.iter().filter(|s| s.0 == 2).max_by(|x, y| x.1.total_cmp(&y.1));
        let mut desc = String::new();
        if let Some(&(_, wkm, dep, k)) = pick {
            // the direction: away from the land, the mean of the land 4-neighbours' opposite
            let (x, y) = ((k % 256) as i32, (k / 256) as i32);
            let (mut dx, mut dy) = (0i32, 0i32);
            for (ax, ay) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let j = (y + ay).rem_euclid(256) as usize * 256 + (x + ax).rem_euclid(256) as usize;
                if on.r.zm.data[j] > 0.0 {
                    dx -= ax;
                    dy -= ay;
                }
            }
            if dx == 0 && dy == 0 {
                dx = 1;
            }
            let len = 160usize; // 2048² cells (31 km world)
            let line = |g: &GridF32| -> Vec<f32> {
                (0..len * 2).map(|i| {
                    let t = i as i32 - len as i32 / 2;
                    let sx = ((x * 8 + 4) + t * dx).rem_euclid(2048) as usize;
                    let sy = ((y * 8 + 4) + t * dy).rem_euclid(2048) as usize;
                    g.data[sy * 2048 + sx]
                }).collect()
            };
            let (za, zb) = (line(&g0), line(&g1));
            let (zmin, zmax) = transect_png_free(&[(za, [120, 120, 120]), (zb, [200, 120, 30])], &[], &out182().join(format!("transect_shelf_seed{seed}.png")));
            desc = format!("the shelf transect at cell {k} (256² width {wkm:.1} km, break depth {dep:.0} m): range {zmin:.0} … {zmax:.0} m (grey without, orange with)");
        }
        let lnd = |g: &GridF32| g.data.iter().filter(|&&z| z > 0.0).count() as f64 * (0.1953125f64 * 0.1953125);
        let med = |g: &GridF32| { let mut v: Vec<f64> = g.data.iter().filter(|&&z| z > 0.0).map(|&z| z as f64).collect(); q_of(&mut v, 0.5) };
        format!("   seed {seed} (2 048², intermediate): land {:.0} → {:.0} km² · median altitude {:.0} → {:.0} m · crop centred on 256² cell {centre} · {desc}\n", lnd(&g0), lnd(&g1), med(&g0), med(&g1))
    });
    eprint!("{s}");
    eprintln!("   tier 3: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F182 -- the marine step per step (one panel seed, intermediate): the durations, the sub-steps, the sources.
#[allow(dead_code)]
fn f182_marine_steps() {
    let m = Margin::continent_default();
    let seed = PANEL[0];
    let p = sel176(seed).unwrap();
    let w = world_oo(seed, p, &m);
    let on = physics_dep(seed, p, &m, &w, rough_a177(), 2, Some(Deposition::declared()), None);
    let sr = on.sed.unwrap();
    let tot: f64 = sr.steps.iter().map(|s| s.0).sum();
    eprintln!("   {} marine steps, Σ dt {:.3} Myr, sub-steps {}", sr.steps.len(), tot / 1e6, sr.marine_substeps);
    for (i, s) in sr.steps.iter().enumerate() {
        if (24..32).contains(&i) || i % 50 == 0 {
            eprintln!("   step {i}: dt {:.0} yr · sub-steps {} · source {:.3e} · sea cells {} · eroded {:.3e} · sea sediment {:.3e}", s.0, s.1, s.2, s.3, s.4, s.5);
        }
    }
}

/// F182 -- the sea's filling (the panel, young and intermediate): per seed the eroded volume (without the deposition) against the
/// ocean's volume at t₀ and at the end, and with the deposition the step at which the sea falls under half, then vanishes.
#[allow(dead_code)]
fn f182_sea_fill() {
    let t0 = Instant::now();
    eprintln!("
==========  F182 the sea's filling . the panel, young and intermediate  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    for ai in [1, 2] {
        for seed in PANEL {
            let Some(p) = sel176(seed) else { continue };
            let s = cached182(&format!("fill_age{ai}_seed{seed}"), || {
                let w = world_oo(seed, p, &m);
                let off = physics_dep(seed, p, &m, &w, rough_a, ai, None, None);
                let on = physics_dep(seed, p, &m, &w, rough_a, ai, Some(Deposition::declared()), None);
                let sr = on.sed.as_ref().unwrap();
                let sea0 = sr.steps.first().map_or(0, |s| s.3);
                let half = sr.steps.iter().position(|s| 2 * s.3 < sea0);
                let gone = sr.steps.iter().position(|s| s.3 == 0);
                let (dur, _) = step_durations(T_F166_YEARS, AGES[ai].1.unwrap(), N_PHYS);
                let at = |i: Option<usize>| i.map_or("never".to_string(), |i| format!("step {i} ({:.1} Myr)", dur[..i].iter().sum::<f64>() / 1e6));
                let phi = Marine::declared().porosity as f64;
                let land_end = on.r.zm.data.iter().filter(|&&z| z > 0.0).count() as f64 / 655.36;
                format!(
                    "   seed {seed} {}: without the deposition, eroded {:.3e} km³ against the ocean {:.3e} km³ at t₀, {:.3e} km³ at the end (eroded / (1 − φ) / ocean₀ = {:.2}) · with it, the sea under half at {}, gone at {} · land at the end {land_end:.0} % (without {:.0} %)
#V {} {} {} {} {} {}
",
                    AGES[ai].0, off.vol.0, off.vol.1, off.vol.2, off.vol.0 / (1.0 - phi) / off.vol.1, at(half), at(gone),
                    off.r.zm.data.iter().filter(|&&z| z > 0.0).count() as f64 / 655.36,
                    off.vol.0, off.vol.1, off.vol.2, half.map_or(-1.0, |v| v as f64), gone.map_or(-1.0, |v| v as f64), land_end
                )
            });
            eprint!("{s}");
        }
    }
    eprintln!("   sea fill: {:.0} s", t0.elapsed().as_secs_f64());
}

// ==============================  F183 -- the sediments' accommodation  ==============================

fn out183() -> PathBuf {
    root().join("docs/reports/relief_method/f183_accommodation/images")
}
fn results183() -> PathBuf {
    root().join("docs/reports/relief_method/f183_accommodation/results")
}
fn cached183(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results183().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results183()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}

/// F183 tier 0 -- C1 alone: per oceanic plate C1 can recycle (no continental cell), the cumulative consumed share of its area over
/// the history; the recyclable share of the map. C1's run does not depend on the physics level's age (the age only maps its
/// steps to durations), so the three ages share these numbers.
#[allow(dead_code)]
fn f183_tier0() {
    let t0 = Instant::now();
    eprintln!("\n==========  F183 tier 0 . C1: the oceanic plates' consumption over the history (15 seeds)  ==========");
    let m = Margin::continent_default();
    let mut rows = Vec::new();
    for seed in seeds179() {
        let Some(p) = sel176(seed) else { continue };
        let s = cached183(&format!("t0_seed{seed}"), || {
            let w = world_oo(seed, p, &m);
            let c = &w.src.consumption;
            let n_int = c.plates.len();
            // the plates recyclable at the first interval, their initial area and their consumed cells summed over the intervals
            let first: Vec<(u16, usize)> = c.plates[0].iter().filter(|q| q.3).map(|q| (q.0, q.1)).collect();
            let mut per = Vec::new();
            let (mut a0, mut cons) = (0usize, 0usize);
            for &(pl, area) in &first {
                let tot: usize = c.plates.iter().map(|iv| iv.iter().find(|q| q.0 == pl && q.3).map_or(0, |q| q.2)).sum();
                // the product form, per interval: 1 − Π (1 − c_j / A_j)
                let keep: f64 = c.plates.iter().map(|iv| iv.iter().find(|q| q.0 == pl && q.3).map_or(1.0, |q| 1.0 - q.2 as f64 / q.1.max(1) as f64)).product();
                per.push(format!("{pl}: {area} cells, consumed {tot} ({:.1} %, product {:.1} %)", 100.0 * tot as f64 / area as f64, 100.0 * (1.0 - keep)));
                a0 += area;
                cons += tot;
            }
            let all: usize = c.plates[0].iter().map(|q| q.1).sum();
            let cont: usize = c.plates[0].iter().filter(|q| !q.3).map(|q| q.1).sum();
            let any_cons: usize = c.plates.iter().flatten().map(|q| q.2).sum();
            let rec_cons: usize = c.plates.iter().flatten().filter(|q| q.3).map(|q| q.2).sum();
            format!(
                "   seed {seed}: {n_int} intervals · recyclable plates {} covering {:.1} % of the map (the continent's plates {:.1} %) · their area-weighted consumed share over the history {:.1} % · consumed cells: all plates {any_cons}, recyclable {rec_cons}\n      per plate: {}\n#T {} {} {} {} {}\n",
                first.len(), 100.0 * a0 as f64 / all as f64, 100.0 * cont as f64 / all as f64, 100.0 * cons as f64 / a0.max(1) as f64, per.join(" · "),
                first.len(), a0 as f64 / all as f64, cons as f64 / a0.max(1) as f64, any_cons, rec_cons
            )
        });
        eprint!("{s}");
        if let Some(t) = lines182(&s, "T").pop() {
            rows.push((seed, t));
        }
    }
    let q = |i: usize, p: f64| { let mut v: Vec<f64> = rows.iter().map(|r| r.1[i]).collect(); if v.is_empty() { f64::NAN } else { q_of(&mut v, p) } };
    eprintln!(
        "\n   tier 0 ({} seeds): recyclable share of the map median {:.1} % (range {:.1}–{:.1} %) · consumed share of the recyclable plates over the history median {:.1} % (range {:.1}–{:.1} %)",
        rows.len(), 100.0 * q(1, 0.5), 100.0 * q(1, 0.0), 100.0 * q(1, 1.0), 100.0 * q(2, 0.5), 100.0 * q(2, 0.0), 100.0 * q(2, 1.0)
    );
    eprintln!("   tier 0: {:.0} s", t0.elapsed().as_secs_f64());
}

use ymir_core::cascade::deposition::{Compaction, Load};
use ymir_core::cascade::history::physics_history_acc;

/// The F183 variants: F182; + A1; + A2; + A3; + A1 + A2 + A3.
fn variants183() -> Vec<(&'static str, Deposition)> {
    let d = Deposition::declared();
    let (l, c) = (Some(Load::declared()), Some(Compaction::declared()));
    vec![
        ("F182", d),
        ("A1", Deposition { load: l, ..d }),
        ("A2", Deposition { compaction: c, ..d }),
        ("A3", Deposition { recycle: true, ..d }),
        ("A123", Deposition { load: l, compaction: c, recycle: true, ..d }),
    ]
}

/// One F183 world at 256² with C1's consumption for the recycling.
fn physics_acc(seed: u64, p: C1Profile, m: &Margin, w: &World, rough_a: f32, age: usize, dep: Deposition) -> Dep182 {
    physics_acc_h(seed, p, m, w, rough_a, age, dep, None)
}

/// The same with a fluvial hook (the reference).
fn physics_acc_h(seed: u64, p: C1Profile, m: &Margin, w: &World, rough_a: f32, age: usize, dep: Deposition, hook: Option<&FluvialHook>) -> Dep182 {
    let kh = KHetero::declared();
    let cover = cover_map_m(seed, p, m, COVER_M);
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k) };
    let cfg = base_cfg(seed, rough_a, Some([0, 0]));
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History, age_years: AGES[age].1, sea_low: None };
    let l = physics_history_acc(&w.snaps, 256, &cfg, &run, Some(&rock), None, Some(dep), hook, Some(&w.src.consumption), &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask_t(&zm, Some([0, 0]));
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let mm = macro_stats(&zm.data, &classes, 256, 256, CELL);
    let ckm2 = (CELL * CELL) as f64;
    let eroded = -l.rec.sum_erosion.iter().map(|&v| v as f64).sum::<f64>() * k as f64 / 1000.0 * ckm2;
    let z0 = to_m(&l.z0, k);
    let ocean = |g: &GridF32| g.data.iter().filter(|&&z| z < 0.0).map(|&z| -z as f64).sum::<f64>() / 1000.0 * ckm2;
    let vol = (eroded, ocean(&z0), ocean(&zm));
    Dep182 { r: Run { zm, d5, classes, m: mm, cost, rebound: l.rec.sum_rebound.iter().map(|v| v * k).collect() }, sed: l.rec.sed.clone(), vol }
}

/// F183 B -- the ocean's capacity at t₀ (km³ of solid) by local Airy columns filled to sea level (an upper bound on the load's
/// effect: no flexure): F182 (φ₀ fixed, no load), A2 (Athy, no load), A1 (φ₀ fixed, loaded), A1 + A2. A diagnostic in the bench
/// (the A1 + A2 column read from a table).
fn capacity183(z0: &GridF32) -> [f64; 4] {
    let c = Compaction::declared();
    let l = Load::declared();
    let (rg, rw, rm) = (l.rho_g as f64, l.rho_w as f64, Rebound::declared().rho_m as f64);
    let phi = c.phi0 as f64;
    let rb = (1.0 - phi) * rg + phi * rw;
    // A1 + A2: the column h with (h − d)(ρm − ρw) = S(h)(ρg − ρw); h ↦ h − S(h)(ρg − ρw)/(ρm − ρw) is increasing
    let dh = 1.0f64;
    let g: Vec<f64> = (0..60_001).map(|i| {
        let h = i as f64 * dh;
        h - c.solid(h) * (rg - rw) / (rm - rw)
    }).collect();
    let inv_g = |d: f64| -> f64 {
        let j = g.partition_point(|&v| v < d).clamp(1, g.len() - 1);
        let (g0, g1) = (g[j - 1], g[j]);
        (j - 1) as f64 * dh + (d - g0) / (g1 - g0) * dh
    };
    let ckm2 = (CELL * CELL) as f64;
    let mut cap = [0f64; 4];
    for &z in &z0.data {
        if z >= 0.0 {
            continue;
        }
        let d = -z as f64;
        cap[0] += (1.0 - phi) * d;
        cap[1] += c.solid(d);
        cap[2] += (1.0 - phi) * d * (rm - rw) / (rm - rb);
        cap[3] += c.solid(inv_g(d));
    }
    cap.map(|v| v / 1000.0 * ckm2)
}

/// The sea's area over the history (a share of the map), per variant (x: time, y: 0 … 0.7).
fn sea_chart(series: &[(Vec<(f64, f64)>, [u8; 3])], path: &std::path::Path) {
    let (w, h) = (420u32, 240u32);
    let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([252, 252, 252, 255]));
    let tmax = series.iter().flat_map(|s| s.0.iter().map(|p| p.0)).fold(1e-9, f64::max);
    for gy in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] {
        let y = (h as f64 - 10.0 - gy / 0.7 * (h as f64 - 20.0)) as u32;
        for x in 10..w - 10 {
            img.put_pixel(x, y, image::Rgba([220, 220, 220, 255]));
        }
    }
    for (pts, c) in series {
        for p in pts {
            let x = (10.0 + p.0 / tmax * (w as f64 - 20.0)) as u32;
            let y = (h as f64 - 10.0 - (p.1 / 0.7).clamp(0.0, 1.0) * (h as f64 - 20.0)) as u32;
            for dy in 0..2 {
                if y + dy < h && x < w {
                    img.put_pixel(x, y + dy, image::Rgba([c[0], c[1], c[2], 255]));
                }
            }
        }
    }
    img.save(path).unwrap();
}

/// The variants' colours: F182 grey, A1 blue, A2 green, A3 orange, A1 + A2 + A3 red.
const COLORS183: [[u8; 3]; 5] = [[120, 120, 120], [30, 90, 200], [30, 160, 70], [200, 120, 20], [200, 30, 30]];

/// F183 tier 1 -- 256², the panel at one age: the five variants (resumable per world and variant).
fn tier183(ai: usize) {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[ai].1.unwrap(), N_PHYS);
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let w = std::cell::OnceCell::new();
        let mut series = Vec::new();
        for (vi, (vn, dv)) in variants183().into_iter().enumerate() {
            let name = format!("t1_age{ai}_seed{seed}_{vn}");
            let s = cached183(&name, || {
                let w = w.get_or_init(|| world_oo(seed, p, &m));
                let r = physics_acc(seed, p, &m, w, rough_a, ai, dv);
                let sr = r.sed.as_ref().unwrap();
                let land: Vec<bool> = r.r.zm.data.iter().map(|&z| z > 0.0).collect();
                let circ = circumnavigation(&land, 256);
                let sea_end = 1.0 - land.iter().filter(|&&l| l).count() as f64 / 65536.0;
                let b = sr.budget;
                let input = b[0] + b[1];
                let tmax = sr.thickness_m.iter().copied().fold(0.0f32, f32::max);
                let kmax = sr.thickness_m.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).map_or(0, |x| x.0);
                let mut tcum = 0.0;
                let pts: Vec<String> = sr.steps.iter().enumerate().map(|(i, st)| {
                    tcum += dur[i.min(N_PHYS - 1)];
                    format!("{:.3}:{:.4}", tcum / 1e6, st.3 as f64 / 65536.0)
                }).collect();
                let gone = sr.steps.iter().position(|st| st.3 == 0);
                // the marine scheme's artifact: marine deposit standing above 100 m (the domain is fixed for a step)
                let mtot: f64 = sr.marine_m.iter().map(|&v| v.max(0.0) as f64).sum();
                let mhigh: f64 = (0..65536).filter(|&k| r.r.zm.data[k] > 100.0).map(|k| sr.marine_m[k].max(0.0) as f64).sum();
                let nhigh = (0..65536).filter(|&k| r.r.zm.data[k] > 100.0 && sr.marine_m[k] > 0.0).count();
                let mut th: Vec<f64> = sr.thickness_m.iter().map(|&t| t as f64).filter(|&t| t > 0.0).collect();
                let p99 = if th.is_empty() { 0.0 } else { q_of(&mut th, 0.99) };
                if vn == "A123" {
                    std::fs::create_dir_all(out183()).unwrap();
                    thick_png(&r.r.zm, &sr.thickness_m, 4000.0).save(out183().join(format!("thickness_A123_age{ai}_seed{seed}.png"))).unwrap();
                }
                format!(
                    "   seed {seed} {} {vn}: {} · the sea at the end {:.1} % of the map (gone at step {:?}) · circumnavigable {} · input {:.4e} (norm × cells) → land {:.1} % · lakes {:.1} % · sea {:.1} % · subducted {:.1} % · flux control {:.1e} · stock control {:.1e} · created (no A2) {:.2} % · the thickest cover {:.0} m at cell {kmax} · median altitude {:.0} m · compaction table error {:.1e} m\n#E {} {} {} {} {} {} {} {} {} {} {}\n#P {}\n",
                    AGES[ai].0, fmt_cost(r.r.cost), 100.0 * sea_end, gone, !circ.wraps, input,
                    100.0 * (b[2] + b[4]) / input, 100.0 * b[3] / input, 100.0 * b[5] / input, 100.0 * b[7] / input, sr.max_rel_residual, sr.max_stock_residual, 100.0 * b[10] / input,
                    tmax, r.r.m.deciles[4], sr.table_err_m,
                    sea_end, !circ.wraps as u8, input, b[7], sr.max_rel_residual, sr.max_stock_residual, tmax, gone.map_or(-1.0, |g| g as f64), r.r.m.deciles[4], r.r.m.fractions[0], b[10],
                    pts.join(" ")
                ) + &format!("      the marine deposit standing above 100 m: {:.1} % of the marine Σ on {nhigh} cells · the cover's 99th percentile {p99:.0} m
#H {} {nhigh} {p99}
", 100.0 * mhigh / mtot.max(1e-30), mhigh / mtot.max(1e-30))
            });
            eprint!("{}", s.lines().filter(|l| !l.starts_with("#P")).map(|l| format!("{l}\n")).collect::<String>());
            let pts: Vec<(f64, f64)> = s.lines().find(|l| l.starts_with("#P ")).map_or(Vec::new(), |l| {
                l[3..].split_whitespace().filter_map(|t| {
                    let (a, b) = t.split_once(':')?;
                    Some((a.parse().ok()?, b.parse().ok()?))
                }).collect()
            });
            series.push((pts, COLORS183[vi]));
        }
        std::fs::create_dir_all(out183()).unwrap();
        sea_chart(&series, &out183().join(format!("sea_age{ai}_seed{seed}.png")));
    }
}

/// The ratio eroded / capacity per world: the eroded volume without the deposition (F182's fill results), the capacity at t₀.
fn ratios183(ai: usize) {
    let m = Margin::continent_default();
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let s = cached183(&format!("cap_age{ai}_seed{seed}"), || {
            let w = world_oo(seed, p, &m);
            let f = std::fs::read_to_string(results182().join(format!("fill_age{ai}_seed{seed}.txt"))).unwrap_or_default();
            let v = lines182(&f, "V").pop().unwrap_or_default();
            let eroded = v.first().copied().unwrap_or(f64::NAN);
            let z0 = to_m(&upsample_to(&w.snaps[0], 256), n2m());
            let cap = capacity183(&z0);
            format!(
                "   seed {seed} {}: eroded without the deposition {:.3e} km³ · the ocean's solid capacity at t₀ (local Airy, km³): F182 {:.3e}, A2 {:.3e}, A1 {:.3e}, A1 + A2 {:.3e} · eroded / capacity {:.2}, {:.2}, {:.2}, {:.2}\n#Q {} {} {} {} {}\n",
                AGES[ai].0, eroded, cap[0], cap[1], cap[2], cap[3], eroded / cap[0], eroded / cap[1], eroded / cap[2], eroded / cap[3], eroded, cap[0], cap[1], cap[2], cap[3]
            )
        });
        eprint!("{s}");
    }
}

/// The tier's summary per variant: the sea kept, the circumnavigation, the subducted share, the controls, the thickest cover.
fn summary183(ai: usize) {
    for (vn, _) in variants183() {
        let (mut kept, mut circ, mut n) = (0, 0, 0);
        let (mut sea, mut subd, mut thick, mut fx, mut st) = (Vec::new(), Vec::new(), Vec::new(), 0f64, 0f64);
        for seed in PANEL {
            let Ok(s) = std::fs::read_to_string(results183().join(format!("t1_age{ai}_seed{seed}_{vn}.txt"))) else { continue };
            let Some(e) = lines182(&s, "E").pop() else { continue };
            n += 1;
            kept += (e[0] > 0.0) as usize;
            circ += (e[1] > 0.0) as usize;
            sea.push(100.0 * e[0]);
            subd.push(100.0 * e[3] / e[2]);
            thick.push(e[6]);
            fx = fx.max(e[4]);
            st = st.max(e[5]);
        }
        let q = |v: &mut Vec<f64>| if v.is_empty() { f64::NAN } else { q_of(v, 0.5) };
        let tmax = thick.iter().copied().fold(0.0, f64::max);
        eprintln!(
            "   {} {vn}: the sea kept on {kept} / {n} · circumnavigable {circ} / {n} · the sea at the end median {:.1} % · subducted / input median {:.1} % · the thickest cover median {:.0} m (max {tmax:.0}) · controls (worst) flux {fx:.1e}, stock {st:.1e}",
            AGES[ai].0, q(&mut sea), q(&mut subd), q(&mut thick)
        );
    }
}

/// F183 tier 1a -- young.
#[allow(dead_code)]
fn f183_tier1_young() {
    let t0 = Instant::now();
    eprintln!("\n==========  F183 tier 1 . 256², young, the panel: F182 / + A1 / + A2 / + A3 / + A1 + A2 + A3  ==========");
    ratios183(1);
    tier183(1);
    summary183(1);
    eprintln!("   tier 1 young: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F183 tier 1b -- intermediate.
#[allow(dead_code)]
fn f183_tier1_intermediate() {
    let t0 = Instant::now();
    eprintln!("\n==========  F183 tier 1 . 256², intermediate, the panel: F182 / + A1 / + A2 / + A3 / + A1 + A2 + A3  ==========");
    ratios183(2);
    tier183(2);
    summary183(2);
    eprintln!("   tier 1 intermediate: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F183 -- the thickest cover under A1 alone (seed 29, young): what built it.
#[allow(dead_code)]
fn f183_diag_thick() {
    let m = Margin::continent_default();
    let seed = 29u64;
    let p = sel176(seed).unwrap();
    let w = world_oo(seed, p, &m);
    let d = Deposition { load: Some(Load::declared()), ..Deposition::declared() };
    let r = physics_acc(seed, p, &m, &w, rough_a177(), 1, d);
    let sr = r.sed.as_ref().unwrap();
    let k = sr.thickness_m.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap().0;
    let (x, y) = (k % 256, k / 256);
    eprintln!("   cell {k} ({x}, {y}): thickness {:.0} m · solid {:.0} m · fluvial+stranded Σ {:.0} m · lake Σ {:.0} m · marine Σ {:.0} m · z {:.0} m", sr.thickness_m[k], sr.solid_m[k], sr.fluvial_m[k], sr.lake_m[k], sr.marine_m[k], r.r.zm.data[k]);
    for dy in -2i32..=2 {
        let row: Vec<String> = (-2i32..=2).map(|dx| {
            let j = ((y as i32 + dy).rem_euclid(256) * 256 + (x as i32 + dx).rem_euclid(256)) as usize;
            format!("{:>6.0}/{:>6.0}", r.r.zm.data[j], sr.thickness_m[j])
        }).collect();
        eprintln!("      {}", row.join("  "));
    }
    let mut big: Vec<(usize, f32)> = sr.thickness_m.iter().copied().enumerate().filter(|t| t.1 > 5000.0).collect();
    big.sort_by(|a, b| b.1.total_cmp(&a.1));
    eprintln!("   cells over 5 km: {} · the 10 thickest: {:?}", big.len(), big.iter().take(10).map(|t| (t.0, t.1.round(), sr.fluvial_m[t.0].round(), sr.lake_m[t.0].round(), sr.marine_m[t.0].round())).collect::<Vec<_>>());
}

/// F183 -- the baseline: the circumnavigation and the sea without the deposition (F181's default), the panel, young and
/// intermediate.
#[allow(dead_code)]
fn f183_baseline() {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    for ai in [1, 2] {
        for seed in PANEL {
            let Some(p) = sel176(seed) else { continue };
            let s = cached183(&format!("base_age{ai}_seed{seed}"), || {
                let w = world_oo(seed, p, &m);
                let r = physics_dep(seed, p, &m, &w, rough_a, ai, None, None);
                let land: Vec<bool> = r.r.zm.data.iter().map(|&z| z > 0.0).collect();
                let circ = circumnavigation(&land, 256);
                let sea = 1.0 - land.iter().filter(|&&l| l).count() as f64 / 65536.0;
                format!("   seed {seed} {} without the deposition: the sea {:.1} % · circumnavigable {} (r* {:?}) · median altitude {:.0} m\n#O {} {} {}\n", AGES[ai].0, 100.0 * sea, !circ.wraps, circ.sea_r, r.r.m.deciles[4], sea, !circ.wraps as u8, r.r.m.deciles[4])
            });
            eprint!("{s}");
        }
    }
}

/// F183: the combination A1 + A2 + A3.
fn a123() -> Deposition {
    Deposition { load: Some(Load::declared()), compaction: Some(Compaction::declared()), recycle: true, ..Deposition::declared() }
}

// ---- F182's tiers on the combination (tier 2 and 3) ----
/// F183 tier 2 (F182's tier 1 on A1 + A2 + A3) -- 256², the 15 seeds, the intermediate age: F181's default without and with the deposition.
#[allow(dead_code)]
fn f183_t2_tier1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out183()).unwrap();
    eprintln!("\n==========  F183 tier 2 (F182's tier 1, the deposition = A1 + A2 + A3) . 256², intermediate: without / with the deposition  ==========");
    tier183b(2, &seeds179(), "b1");
    summary183b("b1", 2, &seeds179());
    for pre in ["lakes_off", "lakes_on", "thickness", "shelf_on"] {
        let tiles: Vec<image::RgbaImage> = seeds179().iter().filter_map(|s| image::open(out183().join(format!("{pre}_age2_seed{s}.png"))).ok().map(|i| i.to_rgba8())).collect();
        if !tiles.is_empty() {
            montage(&tiles, 5, 260, 260).save(out183().join(format!("{pre}_15_age2.png"))).unwrap();
        }
    }
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F183 tier 2 (F182's tier 2 on A1 + A2 + A3) -- the young and old ages, the panel.
#[allow(dead_code)]
fn f183_t2_ages() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out183()).unwrap();
    eprintln!("\n==========  F183 tier 2 (F182's tier 2, A1 + A2 + A3) . 256², young and old, the panel: without / with the deposition  ==========");
    for ai in [1, 3] {
        tier183b(ai, &PANEL, "b2");
        summary183b("b2", ai, &PANEL);
    }
    eprintln!("   tier 2: {:.0} s", t0.elapsed().as_secs_f64());
}

/// The per-seed runs of a tier (resumable).
fn tier183b(ai: usize, seeds: &[u64], tier: &str) {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    for &seed in seeds {
        let Some(p) = sel176(seed) else { continue };
        let name = format!("{tier}_age{ai}_seed{seed}");
        if results183().join(format!("{name}.txt")).exists() {
            eprint!("{}", std::fs::read_to_string(results183().join(format!("{name}.txt"))).unwrap());
            continue;
        }
        let s = cached183(&name, || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let w = world_oo(seed, p, &m);
            let off = physics_dep(seed, p, &m, &w, rough_a, ai, None, None);
            let on = physics_acc(seed, p, &m, &w, rough_a, ai, a123());
            let (s_off, mk_off) = measure182("without", &off, &o, &a, &w);
            let (s_on, mk_on) = measure182("with", &on, &o, &a, &w);
            let (nd, ad, lost) = deltas182(&on.r.zm, &off.r.zm);
            lakes_png(&off.r.zm, &mk_off).save(out183().join(format!("lakes_off_age{ai}_seed{seed}.png"))).unwrap();
            lakes_png(&on.r.zm, &mk_on).save(out183().join(format!("lakes_on_age{ai}_seed{seed}.png"))).unwrap();
            let sr = on.sed.as_ref().unwrap();
            thick_png(&on.r.zm, &sr.thickness_m, 2000.0).save(out183().join(format!("thickness_age{ai}_seed{seed}.png"))).unwrap();
            shelf_png(&on.r.zm, &shelf182(&on.r.zm, &a)).save(out183().join(format!("shelf_on_age{ai}_seed{seed}.png"))).unwrap();
            shelf_png(&off.r.zm, &shelf182(&off.r.zm, &a)).save(out183().join(format!("shelf_off_age{ai}_seed{seed}.png"))).unwrap();
            format!(
                "   seed {seed} {}: physics without {} · with {}\n{s_off}{s_on}      deltas (land gained at the coast) {nd}, {ad:.0} km² · land lost to the sea {lost:.0} km²\n#G {nd} {ad} {lost} {} {}\n",
                AGES[ai].0, fmt_cost(off.r.cost), fmt_cost(on.r.cost), off.r.cost.1, on.r.cost.1
            )
        });
        eprint!("{s}");
    }
}

/// The tier's summary: the medians over the seeds, without → with, and the predictions' counts.
fn summary183b(tier: &str, ai: usize, seeds: &[u64]) {
    let mut rows: Vec<(u64, [Vec<Vec<f64>>; 7], Vec<f64>, Vec<f64>)> = Vec::new();
    for &seed in seeds {
        let Ok(s) = std::fs::read_to_string(results183().join(format!("{tier}_age{ai}_seed{seed}.txt"))) else { continue };
        let ks = ["L", "S", "C", "D", "F", "W", "M"].map(|k| lines182(&s, k));
        let b = lines182(&s, "B").pop().unwrap_or_default();
        let g = lines182(&s, "G").pop().unwrap_or_default();
        rows.push((seed, ks, b, g));
    }
    let n = rows.len();
    let med = |f: &dyn Fn(&(u64, [Vec<Vec<f64>>; 7], Vec<f64>, Vec<f64>), usize) -> f64, i: usize| -> f64 {
        let mut v: Vec<f64> = rows.iter().map(|r| f(r, i)).filter(|x| x.is_finite()).collect();
        if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) }
    };
    let get = |key: usize, col: usize| move |r: &(u64, [Vec<Vec<f64>>; 7], Vec<f64>, Vec<f64>), i: usize| r.1[key].get(i).and_then(|v| v.get(col)).copied().unwrap_or(f64::NAN);
    let pair = |key: usize, col: usize| (med(&get(key, col), 0), med(&get(key, col), 1));
    let large = pair(0, 3);
    let dust = pair(0, 0);
    let la = (med(&|r, i| r.1[0].get(i).map_or(f64::NAN, |v| v[4] + v[5] + v[6] + v[7]), 0), med(&|r, i| r.1[0].get(i).map_or(f64::NAN, |v| v[4] + v[5] + v[6] + v[7]), 1));
    let large_drop: Vec<f64> = rows.iter().filter_map(|r| {
        let (a, b) = (r.1[0].first()?[3], r.1[0].get(1)?[3]);
        (a > 0.0).then(|| 100.0 * (b - a) / a)
    }).collect();
    let sp = pair(1, 0);
    let sa = pair(1, 1);
    let dim = pair(2, 0);
    let lr = pair(2, 1);
    let plain = pair(6, 0);
    let mtn = pair(6, 3);
    let medalt = pair(6, 5);
    let mut ddim = Vec::new();
    let mut dw10 = Vec::new();
    let (mut v_ok, mut circ_off, mut circ_on, mut in_band, mut in_harris) = (0, 0, 0, 0, 0);
    let mut dplain = Vec::new();
    let mut dmed = Vec::new();
    for r in &rows {
        let (c0, c1) = (&r.1[2], &r.1[2]);
        if let (Some(a), Some(b)) = (c0.first(), c1.get(1)) {
            ddim.push((b[0] - a[0]).abs());
            circ_off += (a[2] > 0.0) as usize;
            circ_on += (b[2] > 0.0) as usize;
        }
        if let (Some(a), Some(b)) = (r.1[5].first(), r.1[5].get(1)) {
            dw10.push((b[0] - a[0]).abs());
        }
        if let Some(f) = r.1[4].get(1) {
            v_ok += (f[8] == 0.0) as usize;
        }
        if let Some(s) = r.1[1].get(1) {
            in_band += (s[0] >= 4.65 && s[0] <= 16.4) as usize;
            in_harris += (s[0] >= 34.9 && s[0] <= 123.0) as usize;
        }
        if let (Some(a), Some(b)) = (r.1[6].first(), r.1[6].get(1)) {
            dplain.push(100.0 * (b[0] - a[0]));
            dmed.push(b[5] - a[5]);
        }
    }
    let shares: Vec<[f64; 4]> = rows.iter().filter(|r| r.2.len() >= 9).map(|r| {
        let b = &r.2;
        let t = b[2] + b[3] + b[4] + b[5] + b[6];
        [100.0 * (b[2] + b[4]) / t, 100.0 * b[3] / t, 100.0 * b[5] / t, b[7]]
    }).collect();
    let sh_med = |i: usize| -> f64 { let mut v: Vec<f64> = shares.iter().map(|s| s[i]).collect(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
    let worst = shares.iter().map(|s| s[3]).fold(0.0, f64::max);
    let deltas: Vec<f64> = rows.iter().filter_map(|r| r.3.get(1).copied()).collect();
    let nd: Vec<f64> = rows.iter().filter_map(|r| r.3.first().copied()).collect();
    let land_on: Vec<f64> = rows.iter().filter_map(|r| r.1[6].get(1).map(|v| v[4])).collect();
    let dpct: Vec<f64> = deltas.iter().zip(&land_on).map(|(d, l)| 100.0 * d / l).collect();
    let q = |v: &[f64], p: f64| -> f64 { let mut v = v.to_vec(); if v.is_empty() { f64::NAN } else { q_of(&mut v, p) } };
    eprintln!(
        "\n   F182 {tier} {} ({n} seeds), the medians without → with:\n      large lakes {:.0} → {:.0} (per-seed change median {:.0} %, range {:.0}…{:.0} %) · dust {:.0} → {:.0} · lake area {:.0} → {:.0} km²\n      shelf passive {:.1} → {:.1} km · active {:.1} → {:.1} km · in Harris ÷ 7.5 (4.65–16.4 km) {in_band} / {n} · in Harris unconverted (34.9–123.0 km) {in_harris} / {n}\n      plains {:.1} → {:.1} % (per-seed change median {:+.1} points) · mountain {:.1} → {:.1} % · median altitude {:.0} → {:.0} m (change median {:+.0} m)\n      coast dimension {:.3} → {:.3} (|Δ| max {:.3}) · length ratio {:.3} → {:.3} · deltas per seed median {:.0}, area median {:.2} % of the land\n      the sediment (median shares): land {:.1} % · lakes {:.1} % · sea {:.1} % · the mass control's worst step over the seeds {worst:.2e}\n      (v) holds on {v_ok} / {n} · W10 |Δ| max {:.3} · circumnavigable without {circ_off} / {n}, with {circ_on} / {n}\n",
        AGES[ai].0, large.0, large.1, q(&large_drop, 0.5), q(&large_drop, 0.0), q(&large_drop, 1.0), dust.0, dust.1, la.0, la.1,
        sp.0, sp.1, sa.0, sa.1, 100.0 * plain.0, 100.0 * plain.1, q(&dplain, 0.5), 100.0 * mtn.0, 100.0 * mtn.1, medalt.0, medalt.1, q(&dmed, 0.5),
        dim.0, dim.1, q(&ddim, 1.0), lr.0, lr.1, q(&nd, 0.5), q(&dpct, 0.5), sh_med(0), sh_med(1), sh_med(2), q(&dw10, 1.0)
    );
}

/// ADR Finding 183 (the author's decision) -- **abandoned as a reference**: its convergence is proven only for G / p̃ ≤ 1 (Yuan
/// et al. 2019 JGR p. 8) and G = 1.55; at F183 it never converged in 50 iterates (300 / 300 steps). The code is kept as a record.
/// F182 A3 -- Yuan's iterative coupling, the bench's reference (JGR eq. 19–21): each iterate the cascade's incision of the
/// field plus the previous iterate's deposit (eq. 20's form), the deposit from that iterate's erosion by the cascade's pass,
/// routed on the iterate's surface without its lake fill; stopped at eq. 21 (max |Δh| < 10⁻³ U Δt), at most 50 iterates;
/// mass-exact output z = before − e^K + D^K. Logs per step: (step, iterates, converged, the gap to the first iterate (the
/// splitting): Σ|D^K − D¹| / ΣD^K and the mean |z^K − z¹| on land, norm).
fn yuan_ref183(c: &FluvialCtx, log: &RefCell<Vec<(usize, usize, bool, f64, f64)>>) -> FluvialOut {
    let tol = 1e-3 * c.u_dt;
    let n = c.before.data.len();
    let mut d_prev = vec![0f64; n];
    let mut lake_prev = vec![0f64; n];
    let mut h_prev: Option<GridF32> = None;
    let mut first: Option<(Vec<f64>, Vec<f64>)> = None;
    let mut it = 0;
    loop {
        it += 1;
        let mut surf = c.before.clone();
        for k in 0..n {
            if d_prev[k] != 0.0 {
                surf.data[k] = (surf.data[k] as f64 + d_prev[k] / (1.0 - c.p_dep)) as f32;
            }
        }
        let h = (c.incise)(&surf);
        let e: Vec<f64> = (0..n).map(|k| (surf.data[k] as f64 - h.data[k] as f64).max(0.0)).collect();
        let mut bare = h.clone();
        for k in 0..n {
            if lake_prev[k] != 0.0 {
                bare.data[k] = (bare.data[k] as f64 - lake_prev[k] / (1.0 - c.p_dep)) as f32;
            }
        }
        let p = (c.pass)(&bare, &e);
        let dmax = h_prev.as_ref().map_or(f64::INFINITY, |hp| (0..n).map(|k| (h.data[k] - hp.data[k]).abs() as f64).fold(0.0, f64::max));
        let dk: Vec<f64> = (0..n).map(|k| p.deposit(k)).collect();
        let zk: Vec<f64> = (0..n).map(|k| c.before.data[k] as f64 - e[k] + dk[k] / (1.0 - c.p_dep)).collect();
        if first.is_none() {
            first = Some((dk.clone(), zk.clone()));
        }
        let conv = dmax < tol;
        if conv || it >= 50 {
            let (d1, z1) = first.as_ref().unwrap();
            let sd: f64 = dk.iter().sum();
            let gap_d = if sd > 0.0 { dk.iter().zip(d1).map(|(a, b)| (a - b).abs()).sum::<f64>() / sd } else { 0.0 };
            let lnd: Vec<usize> = (0..n).filter(|&k| c.before.data[k] > 0.5).collect();
            let gap_z = lnd.iter().map(|&k| (zk[k] - z1[k]).abs()).sum::<f64>() / lnd.len().max(1) as f64;
            log.borrow_mut().push((c.step, it, conv, gap_d, gap_z));
            let mut z = c.before.clone();
            for k in 0..n {
                z.data[k] = zk[k] as f32;
            }
            return FluvialOut { z, e, pass: p, iterations: it };
        }
        lake_prev = p.lake.clone();
        d_prev = dk;
        h_prev = Some(h);
    }
}

/// F182 A3 -- the reference against the splitting, on 2 panel seeds at the intermediate age: the deposit maps (correlation, the
/// mass difference), the masses, the lake area, the median altitude; P7: the per-step gap against the step's duration.
#[allow(dead_code)]
fn f183_reference() {
    let t0 = Instant::now();
    eprintln!("\n==========  F183 (F182's A3 on A1 + A2 + A3) . the explicit splitting against Yuan's iterative coupling (bench only), 2 panel seeds, intermediate  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[2].1.unwrap(), N_PHYS);
    for seed in [PANEL[0], PANEL[1]] {
        let p = sel176(seed).unwrap();
        let s = cached183(&format!("ref183_seed{seed}"), || {
            let w = world_oo(seed, p, &m);
            let split = physics_acc(seed, p, &m, &w, rough_a, 2, a123());
            let log = RefCell::new(Vec::new());
            let hook = |c: &FluvialCtx| yuan_ref183(c, &log);
            let refr = physics_acc_h(seed, p, &m, &w, rough_a, 2, a123(), Some(&hook));
            let (ss, rs) = (split.sed.as_ref().unwrap(), refr.sed.as_ref().unwrap());
            let landdep = |r: &SedRecord| -> Vec<f64> { (0..65536).map(|k| (r.fluvial_m[k] + r.lake_m[k]) as f64).collect() };
            let (a, b) = (landdep(ss), landdep(rs));
            let corr = pearson(&a.iter().map(|&v| v as f32).collect::<Vec<_>>(), &b.iter().map(|&v| v as f32).collect::<Vec<_>>());
            let (sa, sb) = (a.iter().sum::<f64>(), b.iter().sum::<f64>());
            let (la, lb) = (lakes182(&split.r.zm), lakes182(&refr.r.zm));
            let area = |l: &([usize; 4], [f64; 4], Vec<bool>)| l.1.iter().sum::<f64>();
            let lg = log.borrow();
            let its: Vec<f64> = lg.iter().map(|l| l.1 as f64).collect();
            let nconv = lg.iter().filter(|l| !l.2).count();
            // P7: the gap at the longest and the shortest steps, and its rank correlation with the duration
            let imax = (0..N_PHYS).max_by(|&x, &y| dur[x].total_cmp(&dur[y])).unwrap();
            let imin = (0..N_PHYS).min_by(|&x, &y| dur[x].total_cmp(&dur[y])).unwrap();
            let at = |i: usize| lg.iter().find(|l| l.0 == i).map(|l| (l.3, l.4 * n2m() as f64, l.1)).unwrap_or((f64::NAN, f64::NAN, 0));
            let dv: Vec<f64> = lg.iter().map(|l| dur[l.0.min(N_PHYS - 1)]).collect();
            let gv: Vec<f64> = lg.iter().map(|l| l.3).collect();
            let zv: Vec<f64> = lg.iter().map(|l| l.4).collect();
            let rho_d = spearman(&dv.iter().map(|&v| v as f32).collect::<Vec<_>>(), &gv.iter().map(|&v| v as f32).collect::<Vec<_>>());
            let rho_z = spearman(&dv.iter().map(|&v| v as f32).collect::<Vec<_>>(), &zv.iter().map(|&v| v as f32).collect::<Vec<_>>());
            let (gl, gs) = (at(imax), at(imin));
            let mut tl = String::new();
            for l in lg.iter() {
                tl += &format!("{} {} {} {} {} {}\n", l.0, dur[l.0.min(N_PHYS - 1)], l.1, l.2 as u8, l.3, l.4 * n2m() as f64);
            }
            std::fs::create_dir_all(results183()).unwrap();
            std::fs::write(results183().join(format!("ref183_steps_seed{seed}.txt")), tl).unwrap();
            let bs = |r: &SedRecord| { let b = r.budget; let t = b[2] + b[3] + b[4] + b[5] + b[6]; (100.0 * (b[2] + b[4]) / t, 100.0 * b[3] / t, 100.0 * b[5] / t, t) };
            let (bsa, bsb) = (bs(ss), bs(rs));
            format!(
                "   seed {seed}: splitting {} · reference {} (iterates per step median {:.0}, max {:.0}; not converged in 50 on {nconv} / {} steps)\n      the land deposit maps (fluvial + lakes, m): r = {corr:.3} · mass {sa:.4e} → {sb:.4e} m·cells ({:+.1} %) · total deposit {:.4e} → {:.4e} norm × cells ({:+.1} %) · shares land / lakes / sea {:.1} / {:.1} / {:.1} % → {:.1} / {:.1} / {:.1} %\n      lake area {:.0} → {:.0} km² · median altitude {:.1} → {:.1} m · the reference's mass control {:.2e}\n      P7: the deposit gap Σ|D_Yuan − D_split| / ΣD_Yuan at the longest step ({:.2} Myr, step {imax}) {:.3} (mean |Δz| {:.2} m, {} iterates), at the shortest ({:.1} kyr, step {imin}) {:.3} (mean |Δz| {:.3} m, {} iterates); Spearman (duration, gap) {rho_d:.2}, (duration, |Δz|) {rho_z:.2}\n#R {corr} {sa} {sb} {} {} {} {} {} {} {} {}\n",
                fmt_cost(split.r.cost), fmt_cost(refr.r.cost), q_of(&mut its.clone(), 0.5), q_of(&mut its.clone(), 1.0), lg.len(),
                100.0 * (sb - sa) / sa, bsa.3, bsb.3, 100.0 * (bsb.3 - bsa.3) / bsa.3, bsa.0, bsa.1, bsa.2, bsb.0, bsb.1, bsb.2,
                area(&la), area(&lb), split.r.m.deciles[4], refr.r.m.deciles[4], rs.max_rel_residual,
                dur[imax] / 1e6, gl.0, gl.1, gl.2, dur[imin] / 1e3, gs.0, gs.1, gs.2,
                area(&la), area(&lb), split.r.m.deciles[4], refr.r.m.deciles[4], gl.0, gs.0, rho_d, rho_z
            )
        });
        eprint!("{s}");
    }
    eprintln!("   reference: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F182 -- the sensitivity (the panel, intermediate): G = 0 and 3.1 (Guerit's bounds), K_M = 300 and 500 (the Ogooué range).
#[allow(dead_code)]
fn f183_sensitivity() {
    let t0 = Instant::now();
    eprintln!("\n==========  F183 the sensitivity (A1 + A2 + A3) . the panel, intermediate: G 0 / 3.1, K_M 300 / 500  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let d = a123();
    let variants = [
        ("G0", Deposition { g: 0.0, ..d }),
        ("G3.1", Deposition { g: 3.1, ..d }),
        ("K300", Deposition { marine: Some(Marine { k_m2_yr: 300.0, ..Marine::declared() }), ..d }),
        ("K500", Deposition { marine: Some(Marine { k_m2_yr: 500.0, ..Marine::declared() }), ..d }),
    ];
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let s = cached183(&format!("sens183_seed{seed}"), || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let w = world_oo(seed, p, &m);
            let mut s = format!("   seed {seed}:\n");
            for (nm, dv) in variants {
                let r = physics_acc(seed, p, &m, &w, rough_a, 2, dv);
                let (line, _) = measure182(nm, &r, &o, &a, &w);
                s += &line.replace("#L", &format!("#L{nm}")).replace("#S", &format!("#S{nm}")).replace("#C", &format!("#C{nm}")).replace("#M", &format!("#M{nm}")).replace("#B", &format!("#B{nm}"));
            }
            s
        });
        eprint!("{s}");
    }
    // the medians per variant
    for (nm, _) in variants {
        let (mut large, mut sp, mut plain, mut med, mut sea) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for seed in PANEL {
            let Ok(s) = std::fs::read_to_string(results183().join(format!("sens183_seed{seed}.txt"))) else { continue };
            if let Some(l) = lines182(&s, &format!("L{nm}")).pop() { large.push(l[3]); }
            if let Some(v) = lines182(&s, &format!("S{nm}")).pop() { if v[0].is_finite() { sp.push(v[0]); } }
            if let Some(v) = lines182(&s, &format!("M{nm}")).pop() { plain.push(100.0 * v[0]); med.push(v[5]); }
            if let Some(b) = lines182(&s, &format!("B{nm}")).pop() { let t = b[2] + b[3] + b[4] + b[5] + b[6]; sea.push(100.0 * b[5] / t); }
        }
        let q = |v: &mut Vec<f64>| if v.is_empty() { f64::NAN } else { q_of(v, 0.5) };
        eprintln!("   {nm}: medians over the panel · large lakes {:.0} · shelf passive {:.1} km · plains {:.1} % · median altitude {:.0} m · sea share {:.1} %", q(&mut large), q(&mut sp), q(&mut plain), q(&mut med), q(&mut sea));
    }
    eprintln!("   sensitivity: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F182 tier 3 -- 2 048², one seed at the intermediate age, without and with the deposition: the whole maps, a crop on the largest
/// delta (or the widest land deposit), a transect across the shelf of a passive coast.
#[allow(dead_code)]
fn f183_tier3() {
    let t0 = Instant::now();
    eprintln!("\n==========  F183 tier 3 (A1 + A2 + A3) . 2 048², intermediate, one seed: without / with the deposition  ==========");
    let m = Margin::continent_default();
    // the seed: the panel seed of tier 1 with the largest delta area
    let delta = |sd: u64| -> f64 { std::fs::read_to_string(results183().join(format!("b1_age2_seed{sd}.txt"))).ok().and_then(|s| lines182(&s, "G").pop()).map(|g| g[1]).unwrap_or(-1.0) };
    let seed = PANEL.iter().copied().max_by(|&x, &y| delta(x).total_cmp(&delta(y)).then(y.cmp(&x))).unwrap();
    let p = sel176(seed).unwrap();
    let s = cached183(&format!("t3_183_seed{seed}"), || {
        let (o, _) = observe179(seed, p, &m);
        let info = step_info(&o.st, &o.kin, o.ds_max, true);
        let a = a1(&info, &o.kin, &o.tracker);
        let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
        let corse_oct = |nc: usize| -> f32 {
            let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
            read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
        };
        let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
        let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
        let kh = KHetero::declared();
        let w = world_oo(seed, p, &m);
        let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
        let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
        let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
        let on = physics_acc(seed, p, &m, &w, rough_a, 2, a123());
        let (g0, _, _) = chain_low(&off.r.zm, amp.clone(), Some(ln_k.clone()), 0.0);
        let (g1, _, _) = chain_low(&on.r.zm, amp, Some(ln_k), 0.0);
        alt_png(&g0, 4).save(out183().join(format!("map2048_off_seed{seed}.png"))).unwrap();
        alt_png(&g1, 4).save(out183().join(format!("map2048_on_seed{seed}.png"))).unwrap();
        // the crop: centred on the 256² cell of the largest delta component, else of the thickest land deposit
        let gain: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && off.r.zm.data[k] <= 0.0).collect();
        let sr = on.sed.as_ref().unwrap();
        let centre = comps8(&gain, 256).into_iter().max_by_key(|c| c.len()).map(|c| c[c.len() / 2]).unwrap_or_else(|| {
            (0..65536).filter(|&k| on.r.zm.data[k] > 0.0).max_by(|&x, &y| sr.thickness_m[x].total_cmp(&sr.thickness_m[y])).unwrap()
        });
        let (cx, cy) = ((centre % 256) * 8, (centre / 256) * 8);
        for (nm, g) in [("off", &g0), ("on", &g1)] {
            let half = 192;
            let mut crop = GridF32::new(2 * half, 2 * half, 0.0);
            for y in 0..2 * half {
                for x in 0..2 * half {
                    let sx = (cx + 2048 - half + x) % 2048;
                    let sy = (cy + 2048 - half + y) % 2048;
                    crop.data[y * 2 * half + x] = g.data[sy * 2048 + sx];
                }
            }
            alt_png(&crop, 1).save(out183().join(format!("crop2048_{nm}_seed{seed}.png"))).unwrap();
        }
        // the transect: from the passive coast cell with the widest shelf (with the deposition), along its seaward path direction
        let sh = shelf182(&on.r.zm, &a);
        let pick = sh.iter().filter(|s| s.0 == 2).max_by(|x, y| x.1.total_cmp(&y.1));
        let mut desc = String::new();
        if let Some(&(_, wkm, dep, k)) = pick {
            // the direction: away from the land, the mean of the land 4-neighbours' opposite
            let (x, y) = ((k % 256) as i32, (k / 256) as i32);
            let (mut dx, mut dy) = (0i32, 0i32);
            for (ax, ay) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let j = (y + ay).rem_euclid(256) as usize * 256 + (x + ax).rem_euclid(256) as usize;
                if on.r.zm.data[j] > 0.0 {
                    dx -= ax;
                    dy -= ay;
                }
            }
            if dx == 0 && dy == 0 {
                dx = 1;
            }
            let len = 160usize; // 2048² cells (31 km world)
            let line = |g: &GridF32| -> Vec<f32> {
                (0..len * 2).map(|i| {
                    let t = i as i32 - len as i32 / 2;
                    let sx = ((x * 8 + 4) + t * dx).rem_euclid(2048) as usize;
                    let sy = ((y * 8 + 4) + t * dy).rem_euclid(2048) as usize;
                    g.data[sy * 2048 + sx]
                }).collect()
            };
            let (za, zb) = (line(&g0), line(&g1));
            let (zmin, zmax) = transect_png_free(&[(za, [120, 120, 120]), (zb, [200, 120, 30])], &[], &out183().join(format!("transect_shelf_seed{seed}.png")));
            desc = format!("the shelf transect at cell {k} (256² width {wkm:.1} km, break depth {dep:.0} m): range {zmin:.0} … {zmax:.0} m (grey without, orange with)");
        }
        let lnd = |g: &GridF32| g.data.iter().filter(|&&z| z > 0.0).count() as f64 * (0.1953125f64 * 0.1953125);
        let med = |g: &GridF32| { let mut v: Vec<f64> = g.data.iter().filter(|&&z| z > 0.0).map(|&z| z as f64).collect(); q_of(&mut v, 0.5) };
        format!("   seed {seed} (2 048², intermediate): land {:.0} → {:.0} km² · median altitude {:.0} → {:.0} m · crop centred on 256² cell {centre} · {desc}\n", lnd(&g0), lnd(&g1), med(&g0), med(&g1))
    });
    eprint!("{s}");
    eprintln!("   tier 3: {:.0} s", t0.elapsed().as_secs_f64());
}

/// F183 -- the losing seeds of tier 2 (A1 + A2 + A3, intermediate): the marine deposit standing above sea level, the sea over
/// the history, the land's change.
#[allow(dead_code)]
fn f183_diag_losers() {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[2].1.unwrap(), N_PHYS);
    for seed in [7u64, 1, 18] {
        let Some(p) = sel176(seed) else { continue };
        let s = cached183(&format!("diag_losers_seed{seed}"), || {
            let w = world_oo(seed, p, &m);
            let r = physics_acc(seed, p, &m, &w, rough_a, 2, a123());
            let sr = r.sed.as_ref().unwrap();
            let mtot: f64 = sr.marine_m.iter().map(|&v| v.max(0.0) as f64).sum();
            let above = |h: f32| -> (f64, usize) {
                let v: f64 = (0..65536).filter(|&k| r.r.zm.data[k] > h).map(|k| sr.marine_m[k].max(0.0) as f64).sum();
                (v / mtot.max(1e-30), (0..65536).filter(|&k| r.r.zm.data[k] > h && sr.marine_m[k] > 0.0).count())
            };
            let (a0, n0) = above(0.0);
            let (a100, n100) = above(100.0);
            let (a500, n500) = above(500.0);
            let mut tcum = 0.0;
            let sea: Vec<String> = sr.steps.iter().enumerate().filter_map(|(i, st)| {
                tcum += dur[i.min(N_PHYS - 1)];
                (i % 30 == 0 || i == sr.steps.len() - 1).then(|| format!("{:.0} Myr {:.1} %", tcum / 1e6, 100.0 * st.3 as f64 / 65536.0))
            }).collect();
            // the land cells with marine sediment (were sea, built up by the marine step) and their median altitude
            let mut hl: Vec<f64> = (0..65536).filter(|&k| r.r.zm.data[k] > 0.0 && sr.marine_m[k] > 0.0).map(|k| r.r.zm.data[k] as f64).collect();
            let med = if hl.is_empty() { 0.0 } else { q_of(&mut hl, 0.5) };
            format!(
                "   seed {seed}: the marine Σ standing above 0 m {:.1} % ({n0} cells), above 100 m {:.1} % ({n100}), above 500 m {:.1} % ({n500}) · land with marine sediment: median altitude {med:.0} m · the sea over the history: {}\n",
                100.0 * a0, 100.0 * a100, 100.0 * a500, sea.join(" · ")
            )
        });
        eprint!("{s}");
    }
}

// ==============================  F184 -- the live marine domain and the fluvial sub-steps  ==============================

fn out184() -> PathBuf {
    root().join("docs/reports/relief_method/f184_marine_fluvial/images")
}
fn results184() -> PathBuf {
    root().join("docs/reports/relief_method/f184_marine_fluvial/results")
}
fn cached184(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results184().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results184()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}

/// The F184 variants: F183's combination (k = 1, the domain fixed); + A1 (the live domain); + A1 + A2 at k = 4 and 16.
fn variants184() -> Vec<(&'static str, Deposition)> {
    let d = a123();
    vec![
        ("F183", d),
        ("A1", Deposition { marine_live: true, ..d }),
        ("A1k4", Deposition { marine_live: true, fluvial_substeps: 4, ..d }),
        ("A1k16", Deposition { marine_live: true, fluvial_substeps: 16, ..d }),
    ]
}

/// F184 -- the lakes born of a deposit: lakes (`filled > z` on land, 8-connected) whose spill cell (the receiver, outside the lake,
/// of a lake cell) carries a cumulated fluvial deposit > 1 m.
fn lakes_born(zm: &GridF32, fluvial_m: &[f32]) -> (usize, usize) {
    let flow = compute_flow(zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let mask: Vec<bool> = (0..65536).map(|k| zm.data[k] > 0.0 && flow.filled.data[k] > zm.data[k]).collect();
    let comps = comps8(&mask, 256);
    let mut born = 0;
    for c in &comps {
        let spill = c.iter().filter_map(|&k| {
            let d = flow.direction[k];
            if d == DIR_NONE {
                return None;
            }
            let (x, y) = ((k % 256) as i32 + D8_DX[d as usize], (k / 256) as i32 + D8_DY[d as usize]);
            let j = y.rem_euclid(256) as usize * 256 + x.rem_euclid(256) as usize;
            (!mask[j]).then_some(j)
        });
        if spill.into_iter().any(|j| fluvial_m[j] > 1.0) {
            born += 1;
        }
    }
    (born, comps.len())
}

/// One world of F184's tier 1: the measures and the convergence statistics as data lines.
fn world184(seed: u64, ai: usize, vn: &str, d: Option<Deposition>, w: &World, p: C1Profile, m: &Margin, rough_a: f32, dur: &[f64]) -> String {
    let r = match d {
        None => physics_dep(seed, p, m, w, rough_a, ai, None, None),
        Some(dv) => physics_acc(seed, p, m, w, rough_a, ai, dv),
    };
    measure184(seed, ai, vn, &r, dur)
}

/// F185: the measures of a world already run (F184's `world184` split in two).
fn measure184(seed: u64, ai: usize, vn: &str, r: &Dep182, dur: &[f64]) -> String {
    let land: Vec<bool> = r.r.zm.data.iter().map(|&z| z > 0.0).collect();
    let circ = circumnavigation(&land, 256);
    let sea_end = 1.0 - land.iter().filter(|&&l| l).count() as f64 / 65536.0;
    let (ln, la, _) = lakes182(&r.r.zm);
    let (rc, _) = read_cell(&r.r.zm, None, CELL);
    let pf = planform(&r.r.zm, CELL, None);
    let mm = &r.r.m;
    let mut s = format!(
        "   seed {seed} {} {vn}: {} · the sea at the end {:.1} % · circumnavigable {} · median altitude {:.0} m · plains {:.1} % · mountain {:.1} % · lakes large {} / dust {} (area {:.0} km²) · density {:.4} · λ {:.2} km · S10 {:.4}\n#W {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n",
        AGES[ai].0, fmt_cost(r.r.cost), 100.0 * sea_end, !circ.wraps, mm.deciles[4], 100.0 * mm.fractions[0], 100.0 * mm.fractions[3], ln[3], ln[0], la.iter().sum::<f64>(),
        rc.density, rc.lambda_km, pf.sinuosity[2],
        sea_end, !circ.wraps as u8, mm.deciles[4], mm.fractions[0], mm.fractions[1], mm.fractions[2], mm.fractions[3], ln[3], ln[0], la[3], la.iter().sum::<f64>(), rc.density, rc.lambda_km, pf.sinuosity[2], r.r.cost.1
    );
    if let Some(sr) = &r.sed {
        let b = sr.budget;
        let input = b[0] + b[1];
        let mtot: f64 = sr.marine_m.iter().map(|&v| v.max(0.0) as f64).sum();
        let above = |h: f32| (0..65536).filter(|&k| r.r.zm.data[k] > h).map(|k| sr.marine_m[k].max(0.0) as f64).sum::<f64>() / mtot.max(1e-30);
        let (a0, a100) = (above(0.0), above(100.0));
        let tmax = sr.thickness_m.iter().copied().fold(0.0f32, f32::max);
        let cap = sr.capped / b[2].max(1e-30);
        let (born, nl) = lakes_born(&r.r.zm, &sr.fluvial_m);
        let fl_m3 = b[2] * n2m() as f64 * (CELL * CELL) as f64 * 1e6 / 1e9;
        let mut tcum = 0.0;
        let pts: Vec<String> = sr.steps.iter().enumerate().map(|(i, st)| {
            tcum += dur[i.min(N_PHYS - 1)];
            format!("{:.3}:{:.4}", tcum / 1e6, st.3 as f64 / 65536.0)
        }).collect();
        s += &format!(
            "      the marine Σ above 0 m {:.2} %, above 100 m {:.2} % · stranded by the marine step {:.2} % of the sea's solid · the fluvial deposit {:.3e} km³, under the cap {:.1} % · lakes born of a deposit {born} / {nl} · the thickest cover {tmax:.0} m · land / lakes / sea {:.1} / {:.1} / {:.1} % · subducted {:.1} % · controls flux {:.1e} stock {:.1e}\n#D {a0} {a100} {} {fl_m3} {cap} {born} {nl} {tmax} {} {} {} {} {}\n#P {}\n",
            100.0 * a0, 100.0 * a100, 100.0 * sr.marine_stranded / b[5].max(1e-30), fl_m3, 100.0 * cap,
            100.0 * (b[2] + b[4]) / input, 100.0 * b[3] / input, 100.0 * b[5] / input, 100.0 * b[7] / input, sr.max_rel_residual, sr.max_stock_residual,
            sr.marine_stranded / b[5].max(1e-30), b[3] / input, b[5] / input, b[7] / input, sr.max_rel_residual, sr.max_stock_residual,
            pts.join(" ")
        );
        if vn.starts_with("A1") {
            std::fs::create_dir_all(out184()).unwrap();
            thick_png(&r.r.zm, &sr.fluvial_m, 2000.0).save(out184().join(format!("fluvial_{vn}_age{ai}_seed{seed}.png"))).unwrap();
            let (_, _, mask) = lakes182(&r.r.zm);
            lakes_png(&r.r.zm, &mask).save(out184().join(format!("lakes_{vn}_age{ai}_seed{seed}.png"))).unwrap();
        }
    }
    s
}

/// F184 tier 1 -- the panel at one age, the chosen variants (resumable per world and variant); "off" is the world without the
/// deposition.
fn tier184(ai: usize, names: &[&str]) {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[ai].1.unwrap(), N_PHYS);
    let vs = variants184();
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let w = std::cell::OnceCell::new();
        for &vn in names {
            let d = if vn == "off" { None } else { Some(vs.iter().find(|v| v.0 == vn).unwrap().1) };
            let s = cached184(&format!("t1_age{ai}_seed{seed}_{vn}"), || {
                let w = w.get_or_init(|| world_oo(seed, p, &m));
                world184(seed, ai, vn, d, w, p, &m, rough_a, &dur)
            });
            eprint!("{}", s.lines().filter(|l| !l.starts_with("#P")).map(|l| format!("{l}\n")).collect::<String>());
        }
    }
}

fn read184(ai: usize, seed: u64, vn: &str) -> Option<(Vec<f64>, Option<Vec<f64>>, Vec<(f64, f64)>)> {
    let s = std::fs::read_to_string(results184().join(format!("t1_age{ai}_seed{seed}_{vn}.txt"))).ok()?;
    let w = lines182(&s, "W").pop()?;
    let d = lines182(&s, "D").pop();
    let pts = s.lines().find(|l| l.starts_with("#P ")).map_or(Vec::new(), |l| {
        l[3..].split_whitespace().filter_map(|t| {
            let (a, b) = t.split_once(':')?;
            Some((a.parse().ok()?, b.parse().ok()?))
        }).collect()
    });
    Some((w, d, pts))
}

/// The tier's summary: per variant against « off », the convergence (A1 vs A1k4, A1k4 vs A1k16), the E1–E3 expectations, the sea
/// charts.
fn summary184(ai: usize) {
    let names = ["F183", "A1", "A1k4", "A1k16"];
    for vn in names {
        let (mut n, mut kept, mut lost, mut a100, mut strd, mut cap, mut tmax, mut fx) = (0, 0, 0, Vec::new(), Vec::new(), Vec::new(), 0f64, 0f64);
        let (mut dmed, mut dpl, mut dmt, mut dla, mut born) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let (mut e1, mut e2, mut e3) = (0, 0, 0);
        for seed in PANEL {
            let (Some(off), Some(on)) = (read184(ai, seed, "off"), read184(ai, seed, vn)) else { continue };
            let (wo, wn) = (&off.0, &on.0);
            n += 1;
            kept += (wn[0] > 0.0) as usize;
            lost += (wo[1] > 0.0 && wn[1] == 0.0) as usize;
            dmed.push(wn[2] - wo[2]);
            dpl.push(100.0 * (wn[3] - wo[3]));
            dmt.push(100.0 * (wn[6] - wo[6]));
            dla.push(if wo[9] > 0.0 { 100.0 * (wn[9] - wo[9]) / wo[9] } else { f64::NAN });
            e1 += (wn[3] >= wo[3]) as usize;
            e2 += (wn[6] - wo[6] <= 0.02) as usize;
            e3 += (wn[9] < wo[9]) as usize;
            if let Some(d) = &on.1 {
                a100.push(100.0 * d[1]);
                strd.push(100.0 * d[2]);
                cap.push(100.0 * d[4]);
                born.push(d[5]);
                tmax = tmax.max(d[7]);
                fx = fx.max(d[11]).max(d[12]);
            }
        }
        let q = |v: &[f64]| { let mut v: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
        let mx = |v: &[f64]| v.iter().copied().fold(f64::NAN, f64::max);
        eprintln!(
            "   {} {vn} ({n} seeds): the sea kept {kept} / {n} · the circumnavigation lost (against off) {lost} · the marine Σ above 100 m median {:.2} % (max {:.2} %) · stranded median {:.2} % · under the cap median {:.1} % · lakes born of a deposit median {:.0} · the thickest cover max {tmax:.0} m · Δ median altitude median {:+.0} m · Δ plains median {:+.1} points · Δ mountain median {:+.1} points · Δ the large lakes' area median {:+.0} % · E1 {e1} / {n} · E2 {e2} / {n} · E3 {e3} / {n} · controls worst {fx:.1e}",
            AGES[ai].0, q(&a100), mx(&a100), q(&strd), q(&cap), q(&born), q(&dmed), q(&dpl), q(&dmt), q(&dla)
        );
    }
    // the convergence (the project's criterion, F174): density, λ, S10 within 5 %, the classes' largest move < 1 point
    for (a, b) in [("A1", "A1k4"), ("A1k4", "A1k16")] {
        let mut ok = 0;
        let mut n = 0;
        let mut worst = (0f64, 0f64, 0f64, 0f64);
        let (mut dvol, mut dmed) = (Vec::new(), Vec::new());
        for seed in PANEL {
            let (Some(ra), Some(rb)) = (read184(ai, seed, a), read184(ai, seed, b)) else { continue };
            let (x, y) = (&ra.0, &rb.0);
            let rel = |u: f64, v: f64| 100.0 * (v - u) / u.abs().max(1e-12);
            let (dd, dl, ds) = (rel(x[11], y[11]), rel(x[12], y[12]), rel(x[13], y[13]));
            let mv = (3..7).map(|i| 100.0 * (y[i] - x[i]).abs()).fold(0.0, f64::max);
            let pass = dd.abs() < 5.0 && dl.abs() < 5.0 && ds.abs() < 5.0 && mv < 1.0;
            n += 1;
            ok += pass as usize;
            worst = (worst.0.max(dd.abs()), worst.1.max(dl.abs()), worst.2.max(ds.abs()), worst.3.max(mv));
            if let (Some(da), Some(db)) = (&ra.1, &rb.1) {
                dvol.push(rel(da[3], db[3]));
            }
            dmed.push(rel(x[2], y[2]));
            eprintln!("      seed {seed} {a} → {b}: density {dd:+.1} % · λ {dl:+.1} % · S10 {ds:+.1} % · the classes' largest move {mv:.2} points · converged {pass}");
        }
        let q = |v: &[f64]| { let mut v = v.to_vec(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
        eprintln!(
            "   {} convergence {a} → {b}: holds on {ok} / {n} · worst density {:.1} %, λ {:.1} %, S10 {:.1} %, classes {:.2} points · the fluvial volume's change median {:+.1} % · the median altitude's change median {:+.1} %",
            AGES[ai].0, worst.0, worst.1, worst.2, worst.3, q(&dvol), q(&dmed)
        );
    }
    // the sea over time per variant (grey F183, blue A1, green A1k4, red A1k16)
    let cols: [[u8; 3]; 4] = [[120, 120, 120], [30, 90, 200], [30, 160, 70], [200, 30, 30]];
    for seed in PANEL {
        let series: Vec<(Vec<(f64, f64)>, [u8; 3])> = names.iter().zip(cols).filter_map(|(vn, c)| read184(ai, seed, vn).map(|r| (r.2, c))).collect();
        if !series.is_empty() {
            std::fs::create_dir_all(out184()).unwrap();
            sea_chart(&series, &out184().join(format!("sea_age{ai}_seed{seed}.png")));
        }
    }
}

/// F184 tier 1a -- intermediate: off, F183, A1, A1k4.
#[allow(dead_code)]
fn f184_t1_int_a() {
    let t0 = Instant::now();
    eprintln!("\n==========  F184 tier 1 . the panel, intermediate: off / F183 / + A1 / + A1 + A2 (k = 4)  ==========");
    tier184(2, &["off", "F183", "A1", "A1k4"]);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

/// F184 tier 1b -- intermediate: A1k16, then the summary.
#[allow(dead_code)]
fn f184_t1_int_b() {
    let t0 = Instant::now();
    eprintln!("\n==========  F184 tier 1 . the panel, intermediate: + A1 + A2 (k = 16)  ==========");
    tier184(2, &["A1k16"]);
    summary184(2);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

/// F184 tier 1c -- old: off, F183, A1.
#[allow(dead_code)]
fn f184_t1_old_a() {
    let t0 = Instant::now();
    eprintln!("\n==========  F184 tier 1 . the panel, old: off / F183 / + A1  ==========");
    tier184(3, &["off", "F183", "A1"]);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

/// F184 tier 1d -- old: A1k4.
#[allow(dead_code)]
fn f184_t1_old_b() {
    let t0 = Instant::now();
    eprintln!("\n==========  F184 tier 1 . the panel, old: + A1 + A2 (k = 4)  ==========");
    tier184(3, &["A1k4"]);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

/// F184 tier 1e -- old: A1k16, then the summary.
#[allow(dead_code)]
fn f184_t1_old_c() {
    let t0 = Instant::now();
    eprintln!("\n==========  F184 tier 1 . the panel, old: + A1 + A2 (k = 16)  ==========");
    tier184(3, &["A1k16"]);
    summary184(3);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

/// F184 -- the cost of one world per variant (seed 8, intermediate and old), announced before the tiers.
#[allow(dead_code)]
fn f184_cost() {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let seed = PANEL[0];
    let p = sel176(seed).unwrap();
    let w = world_oo(seed, p, &m);
    for ai in [2, 3] {
        for (vn, d) in variants184() {
            if vn == "A1k16" {
                continue;
            }
            let r = physics_acc(seed, p, &m, &w, rough_a, ai, d);
            eprintln!("   {} {vn}: {}", AGES[ai].0, fmt_cost(r.r.cost));
        }
    }
}

// ==============================  F185 -- the fluvial sub-step bounded by a duration  ==============================

fn out185() -> PathBuf {
    root().join("docs/reports/relief_method/f185_fluvial_dt/images")
}
fn results185() -> PathBuf {
    root().join("docs/reports/relief_method/f185_fluvial_dt/results")
}
fn cached185(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results185().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results185()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}

/// The F185 variants: F184's live configuration with Δt_f = 80, 20, 5 kyr.
fn variants185() -> Vec<(&'static str, Deposition)> {
    let d = Deposition { marine_live: true, ..a123() };
    vec![
        ("dt80", Deposition { fluvial_dt_max: 80.0e3, ..d }),
        ("dt20", Deposition { fluvial_dt_max: 20.0e3, ..d }),
        ("dt5", Deposition { fluvial_dt_max: 5.0e3, ..d }),
    ]
}

/// E1 / E2 split by zone: among the cells that are land in both runs, the zone in net deposition (the cover ≥ 1 m at the end) and
/// the rest (net erosion); per zone the plains' and the mountain's shares with and without the deposition (points).
fn zones185(on: &Dep182, off_classes: &[u8], off_z: &GridF32) -> [f64; 6] {
    let sr = on.sed.as_ref().unwrap();
    let (mut nd, mut ne) = (0f64, 0f64);
    let (mut pd, mut pd0, mut md, mut md0, mut pe, mut pe0, mut me, mut me0) = (0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
    for k in 0..65536 {
        if on.r.zm.data[k] <= 0.0 || off_z.data[k] <= 0.0 {
            continue;
        }
        let (c, c0) = (on.r.classes[k], off_classes[k]);
        let dep = sr.thickness_m[k] >= 1.0;
        let (p, p0, m, m0) = ((c == PLAIN) as u8 as f64, (c0 == PLAIN) as u8 as f64, (c == MOUNTAIN) as u8 as f64, (c0 == MOUNTAIN) as u8 as f64);
        if dep {
            nd += 1.0;
            pd += p;
            pd0 += p0;
            md += m;
            md0 += m0;
        } else {
            ne += 1.0;
            pe += p;
            pe0 += p0;
            me += m;
            me0 += m0;
        }
    }
    let f = |a: f64, n: f64| if n > 0.0 { 100.0 * a / n } else { f64::NAN };
    [100.0 * nd / (nd + ne).max(1.0), f(pd - pd0, nd), f(md - md0, nd), f(pe - pe0, ne), f(me - me0, ne), nd + ne]
}

/// F185 tier 1 -- the panel at one age, the chosen variants (resumable per world and variant); "off" is the world without the
/// deposition.
fn tier185(ai: usize, names: &[&str]) {
    tier185_seeds(ai, names, &PANEL);
}

/// The same on chosen seeds (so that a long tier splits under 45 min).
fn tier185_seeds(ai: usize, names: &[&str], seeds: &[u64]) {
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[ai].1.unwrap(), N_PHYS);
    let vs = variants185();
    for &seed in seeds {
        let Some(p) = sel176(seed) else { continue };
        let w = std::cell::OnceCell::new();
        let off = std::cell::OnceCell::new();
        for &vn in names {
            let s = cached185(&format!("t1_age{ai}_seed{seed}_{vn}"), || {
                let w = w.get_or_init(|| world_oo(seed, p, &m));
                // the world without the deposition, run once per seed (the zones read it)
                let o = off.get_or_init(|| physics_dep(seed, p, &m, w, rough_a, ai, None, None));
                if vn == "off" {
                    return measure184(seed, ai, vn, o, &dur);
                }
                let d = vs.iter().find(|v| v.0 == vn).unwrap().1;
                let on = physics_acc(seed, p, &m, w, rough_a, ai, d);
                let mut s = measure184(seed, ai, vn, &on, &dur);
                let z = zones185(&on, &o.r.classes, &o.r.zm);
                s += &format!(
                    "      zones: net deposition {:.1} % of the land · there Δ plains {:+.1} / Δ mountain {:+.1} points · in net erosion Δ plains {:+.1} / Δ mountain {:+.1} points\n#Z {} {} {} {} {}\n",
                    z[0], z[1], z[2], z[3], z[4], z[0], z[1], z[2], z[3], z[4]
                );
                s
            });
            eprint!("{}", s.lines().filter(|l| !l.starts_with("#P")).map(|l| format!("{l}\n")).collect::<String>());
        }
    }
}

fn read185(ai: usize, seed: u64, vn: &str) -> Option<(Vec<f64>, Option<Vec<f64>>, Vec<(f64, f64)>, Option<Vec<f64>>)> {
    let s = std::fs::read_to_string(results185().join(format!("t1_age{ai}_seed{seed}_{vn}.txt"))).ok()?;
    let w = lines182(&s, "W").pop()?;
    let d = lines182(&s, "D").pop();
    let z = lines182(&s, "Z").pop();
    let pts = s.lines().find(|l| l.starts_with("#P ")).map_or(Vec::new(), |l| {
        l[3..].split_whitespace().filter_map(|t| {
            let (a, b) = t.split_once(':')?;
            Some((a.parse().ok()?, b.parse().ok()?))
        }).collect()
    });
    Some((w, d, pts, z))
}

/// The tier's summary: per variant against « off », the zones, the convergence (80 → 20, 20 → 5), the extrapolated limit.
fn summary185(ai: usize) {
    let names = ["dt80", "dt20", "dt5"];
    let mut gap_med: Vec<(&str, f64, f64, f64)> = Vec::new();
    for vn in names {
        let (mut n, mut kept, mut lost, mut fx) = (0, 0, 0, 0f64);
        let (mut dmed, mut dpl, mut dmt, mut dla, mut e3, mut cost) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), 0, Vec::new());
        let mut zs: Vec<Vec<f64>> = Vec::new();
        for seed in PANEL {
            let (Some(off), Some(on)) = (read185(ai, seed, "off"), read185(ai, seed, vn)) else { continue };
            let (wo, wn) = (&off.0, &on.0);
            n += 1;
            kept += (wn[0] > 0.0) as usize;
            lost += (wo[1] > 0.0 && wn[1] == 0.0) as usize;
            dmed.push(wn[2] - wo[2]);
            dpl.push(100.0 * (wn[3] - wo[3]));
            dmt.push(100.0 * (wn[6] - wo[6]));
            dla.push(if wo[9] > 0.0 { 100.0 * (wn[9] - wo[9]) / wo[9] } else { f64::NAN });
            e3 += (wn[9] < wo[9]) as usize;
            cost.push(wn[14]);
            if let Some(d) = &on.1 {
                fx = fx.max(d[11]).max(d[12]);
            }
            if let Some(z) = &on.3 {
                zs.push(z.clone());
            }
        }
        if n == 0 {
            continue;
        }
        let q = |v: &[f64]| { let mut v: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
        let zq = |i: usize| q(&zs.iter().map(|z| z[i]).collect::<Vec<_>>());
        eprintln!(
            "   {} {vn} ({n} seeds): the sea kept {kept} / {n} · the circumnavigation lost (against off) {lost} · Δ median altitude median {:+.0} m · Δ plains {:+.1} · Δ mountain {:+.1} points · Δ the large lakes' area {:+.0} % · E3 {e3} / {n} · zones: net deposition {:.1} % of the land, Δ plains / mountain there {:+.1} / {:+.1}, in net erosion {:+.1} / {:+.1} points · controls worst {fx:.1e} · CPU per world median {:.0} s",
            AGES[ai].0, q(&dmed), q(&dpl), q(&dmt), q(&dla), zq(0), zq(1), zq(2), zq(3), zq(4), q(&cost)
        );
        gap_med.push((vn, q(&dmed), q(&dpl), q(&dmt)));
    }
    // the convergence (F174's criterion), 80 → 20 and 20 → 5
    for (a, b) in [("dt80", "dt20"), ("dt20", "dt5")] {
        let (mut ok, mut n) = (0, 0);
        let mut worst = (0f64, 0f64, 0f64, 0f64);
        for seed in PANEL {
            let (Some(ra), Some(rb)) = (read185(ai, seed, a), read185(ai, seed, b)) else { continue };
            let (x, y) = (&ra.0, &rb.0);
            let rel = |u: f64, v: f64| 100.0 * (v - u) / u.abs().max(1e-12);
            let (dd, dl, ds) = (rel(x[11], y[11]), rel(x[12], y[12]), rel(x[13], y[13]));
            let mv = (3..7).map(|i| 100.0 * (y[i] - x[i]).abs()).fold(0.0, f64::max);
            let pass = dd.abs() < 5.0 && dl.abs() < 5.0 && ds.abs() < 5.0 && mv < 1.0;
            n += 1;
            ok += pass as usize;
            worst = (worst.0.max(dd.abs()), worst.1.max(dl.abs()), worst.2.max(ds.abs()), worst.3.max(mv));
            eprintln!("      seed {seed} {a} → {b}: density {dd:+.1} % · λ {dl:+.1} % · S10 {ds:+.1} % · the classes' largest move {mv:.2} points · median altitude {:+.0} m · converged {pass}", y[2] - x[2]);
        }
        if n > 0 {
            eprintln!("   {} convergence {a} → {b}: holds on {ok} / {n} · worst density {:.1} %, λ {:.1} %, S10 {:.1} %, classes {:.2} points", AGES[ai].0, worst.0, worst.1, worst.2, worst.3);
        }
    }
    // the limit, extrapolated geometrically from the gaps to « off » at 80, 20 (and 5) kyr
    if gap_med.len() >= 2 {
        for (i, name) in [(1usize, "Δ median altitude (m)"), (2, "Δ plains (points)"), (3, "Δ mountain (points)")] {
            let v: Vec<f64> = gap_med.iter().map(|g| [g.1, g.2, g.3][i - 1]).collect();
            let lim = if v.len() >= 3 && (v[0] - 2.0 * v[1] + v[2]).abs() > 1e-9 {
                // Aitken's Δ² on the three
                v[2] - (v[2] - v[1]).powi(2) / (v[2] - 2.0 * v[1] + v[0])
            } else {
                f64::NAN
            };
            eprintln!("   {} {name}: {:?} · the limit (Aitken's Δ² on 80 / 20 / 5) {lim:.1}", AGES[ai].0, v.iter().map(|x| (x * 10.0).round() / 10.0).collect::<Vec<_>>());
        }
    }
    let cols: [[u8; 3]; 3] = [[30, 90, 200], [30, 160, 70], [200, 30, 30]];
    for seed in PANEL {
        let series: Vec<(Vec<(f64, f64)>, [u8; 3])> = names.iter().zip(cols).filter_map(|(vn, c)| read185(ai, seed, vn).map(|r| (r.2, c))).collect();
        if !series.is_empty() {
            std::fs::create_dir_all(out185()).unwrap();
            sea_chart(&series, &out185().join(format!("sea_age{ai}_seed{seed}.png")));
        }
    }
}

/// F185 tier 1 -- one test per (age, variants) so that each stays under 45 min of wall time.
#[allow(dead_code)]
fn f185_t1_int_a() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 tier 1 . the panel, intermediate: off / 80 kyr / 20 kyr  ==========");
    tier185(2, &["off", "dt80", "dt20"]);
    summary185(2);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

#[allow(dead_code)]
fn f185_t1_old_a() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 tier 1 . the panel, old: off / 80 kyr  ==========");
    tier185(3, &["off", "dt80"]);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

#[allow(dead_code)]
fn f185_t1_old_b() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 tier 1 . the panel, old: 20 kyr (seeds 8, 18, 19)  ==========");
    tier185_seeds(3, &["dt20"], &[8, 18, 19]);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

#[allow(dead_code)]
fn f185_t1_old_c() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 tier 1 . the panel, old: 20 kyr (seeds 23, 29)  ==========");
    tier185_seeds(3, &["dt20"], &[23, 29]);
    summary185(3);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

#[allow(dead_code)]
fn f185_t1_int_5() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 tier 1 . the panel, intermediate: 5 kyr  ==========");
    tier185(2, &["dt5"]);
    summary185(2);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

#[allow(dead_code)]
fn f185_t1_old_5() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 tier 1 . the panel, old: 5 kyr  ==========");
    tier185(3, &["dt5"]);
    summary185(3);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

/// F185 step 2 (the author's reduced plan): 80 against 20 kyr at the intermediate age, seed 23 run here, seeds 8, 18, 19 already
/// computed; F174's criterion on the four.
#[allow(dead_code)]
fn f185_step2() {
    let t0 = Instant::now();
    eprintln!("\n==========  F185 step 2 . intermediate, seeds 8, 18, 19, 23: 80 against 20 kyr  ==========");
    tier185_seeds(2, &["off", "dt80", "dt20"], &[8, 18, 19, 23]);
    summary185(2);
    eprintln!("   {:.0} s", t0.elapsed().as_secs_f64());
}

// ==============================  F186 -- the deposition by closures and the chain's budget  ==============================


fn out186() -> PathBuf {
    root().join("docs/reports/relief_method/f186_deposition_closure/images")
}
fn results186() -> PathBuf {
    root().join("docs/reports/relief_method/f186_deposition_closure/results")
}

/// F186 0.1 -- the chain's budget: seed 23, the intermediate age, F181's default profile, no deposition. Every stage timed (wall
/// and CPU); the amplification levels timed per operator; 4 096² and 8 192² extrapolated per operator (cost ∝ cells × the
/// budget's iterations, from the 2 048² level). 0.2's split of the deposition elements at 256² (one run of F184's live combination).
#[test]
#[ignore]
fn f186_budget() {
    let t_all = Instant::now();
    let m = Margin::continent_default();
    let seed = 23u64;
    let mut out = String::new();
    // the selection (the real chain recomputes it; the benches read it from F177's files)
    let clk = Clock::start();
    let (np, tried) = select_continent_default_margin(seed, &m);
    let c_sel = clk.stop();
    let cands: usize = tried.iter().map(|t| t.1.candidates.len()).sum();
    let runs: usize = tried.iter().map(|t| t.1.runs).sum();
    out += &format!("   the selection (select_continent_default_margin, {} plate counts tried, {cands} candidates, {runs} C1 runs, chose {np:?}): {}\n", tried.len(), fmt_cost(c_sel));
    let p = sel176(seed).unwrap();
    // C1 on the chosen continent
    let clk = Clock::start();
    let w = world_oo(seed, p, &m);
    let c_c1 = clk.stop();
    out += &format!("   C1 (64², the sources, Airy): {}\n", fmt_cost(c_c1));
    // the physics level, no deposition
    let rough_a = rough_a177();
    let clk = Clock::start();
    let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
    let c_phys = clk.stop();
    out += &format!("   the physics level (256², 300 steps, no deposition): {}\n", fmt_cost(c_phys));
    // the amplification levels, per operator
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let _ = c128;
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let kh = KHetero::declared();
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    let oc = ocean_mask_t(&off.r.zm, Some([0, 0]));
    let d5: Vec<bool> = (0..off.r.zm.data.len()).map(|i| !oc[i] && off.r.zm.data[i] <= 0.0).collect();
    let mut ch = Chain::new(off.r.zm.clone(), amp.clone()).with_d5(d5).with_ln_k(ln_k);
    let mut per_level: Vec<(usize, [f64; 4], [f64; 4], [usize; 3])> = Vec::new();
    for _ in 0..3 {
        let cfg = ch.cfg.clone();
        let t = |f: &mut dyn FnMut()| -> (f64, f64) {
            let c = Clock::start();
            f();
            c.stop()
        };
        let mut lvl_n = 0;
        let (u_w, u_c) = t(&mut || {
            let l = ch.start_level();
            lvl_n = l.n;
        });
        let b = cfg.budget(lvl_n);
        let (e_w, e_c) = t(&mut || ch.levels.last_mut().unwrap().erode(b.erosion, &cfg, &|| false));
        let (t_w, t_c) = t(&mut || ch.levels.last_mut().unwrap().talus(b.talus, &cfg, &|| false));
        let (d_w, d_c) = t(&mut || ch.levels.last_mut().unwrap().deposit(b.deposit, &cfg, &|| false));
        per_level.push((lvl_n, [u_w, e_w, t_w, d_w], [u_c, e_c, t_c, d_c], [b.erosion, b.talus, b.deposit]));
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let mut amp_wall = 0.0;
    for (n, w4, c4, b3) in &per_level {
        amp_wall += w4.iter().sum::<f64>();
        out += &format!(
            "   amplification {n}²: upscale {:.1} s · erosion {:.1} s ({} iterations) · talus {:.1} s ({}) · deposit {:.1} s ({}) · total {:.1} s wall / {:.1} s CPU\n",
            w4[0], w4[1], b3[0], w4[2], b3[1], w4[3], b3[2], w4.iter().sum::<f64>(), c4.iter().sum::<f64>()
        );
    }
    // the extrapolation to 4 096² and 8 192²: per operator, the 2 048² wall per cell-iteration × cells × the budget's iterations
    let (n2, w2, _, b2) = per_level.last().unwrap().clone();
    let cells2 = (n2 * n2) as f64;
    let per = [w2[0] / cells2, w2[1] / (cells2 * b2[0] as f64), w2[2] / (cells2 * b2[1] as f64), w2[3] / (cells2 * b2[2] as f64)];
    let mut ext_total = 0.0;
    for n in [4096usize, 8192] {
        let b = amp.budget(n);
        let cells = (n * n) as f64;
        let e = [per[0] * cells, per[1] * cells * b.erosion as f64, per[2] * cells * b.talus as f64, per[3] * cells * b.deposit as f64];
        ext_total += e.iter().sum::<f64>();
        out += &format!(
            "   amplification {n}² (EXTRAPOLATED: the 2 048² cost per cell-iteration × cells × the budget, {} / {} / {} iterations; Budget::moyen's 4 096² and 8 192² rows are a declared PROXY): upscale {:.0} s · erosion {:.0} s · talus {:.0} s · deposit {:.0} s · total {:.0} s\n",
            b.erosion, b.talus, b.deposit, e[0], e[1], e[2], e[3], e.iter().sum::<f64>()
        );
    }
    let measured = c_sel.0 + c_c1.0 + c_phys.0 + amp_wall;
    out += &format!(
        "   the chain to 2 048² (measured): {:.0} s wall · to 8 192² (with the extrapolation): {:.0} s wall, against the 300 s target ({:+.0} s)\n",
        measured, measured + ext_total, measured + ext_total - 300.0
    );
    // 0.2: the deposition elements at 256² (one run of F184's live combination: the erosion stage, the marine step)
    let d = Deposition { marine_live: true, ..a123() };
    let clk = Clock::start();
    let on = physics_acc(seed, p, &m, &w, rough_a, 2, d);
    let c_on = clk.stop();
    let sr = on.sed.as_ref().unwrap();
    out += &format!(
        "   0.2 · the physics level with F184's live combination (k = 1): {} · of which the marine step {:.1} s wall ({} sub-steps) · against no deposition {}\n",
        fmt_cost(c_on), sr.secs_marine, sr.marine_substeps, fmt_cost(c_phys)
    );
    out += &format!("   the bench: {:.0} s\n", t_all.elapsed().as_secs_f64());
    std::fs::create_dir_all(results186()).unwrap();
    std::fs::write(results186().join("budget_seed23.txt"), &out).unwrap();
    eprint!("{out}");
}

fn cached186(name: &str, f: impl FnOnce() -> String) -> String {
    let p = results186().join(format!("{name}.txt"));
    if let Ok(s) = std::fs::read_to_string(&p) {
        return s;
    }
    let s = f();
    std::fs::create_dir_all(results186()).unwrap();
    std::fs::write(&p, &s).unwrap();
    s
}

/// F186's configuration A: F183's combination, the live marine domain, the fluvial deposition as Yuan's K_eff.
fn a186() -> Deposition {
    Deposition { marine_live: true, fluvial_keff: Some(Deposition::G_DECLARED), ..a123() }
}

/// The medians of the slope (°) and of the altitude (m) on the land cells of a mask.
fn slope_alt(zm: &GridF32, mask: &[bool]) -> (f64, f64, usize) {
    let s = ymir_core::cascade::physio::slope_deg(&zm.data, 256, 256, CELL);
    let mut sv: Vec<f64> = (0..65536).filter(|&k| mask[k] && zm.data[k] > 0.0).map(|k| s[k] as f64).collect();
    let mut av: Vec<f64> = (0..65536).filter(|&k| mask[k] && zm.data[k] > 0.0).map(|k| zm.data[k] as f64).collect();
    let n = sv.len();
    if n == 0 {
        return (f64::NAN, f64::NAN, 0);
    }
    (q_of(&mut sv, 0.5), q_of(&mut av, 0.5), n)
}

/// F186 tier 1 -- 256², the panel, intermediate: without the deposition and with A (resumable per seed).
#[test]
#[ignore]
fn f186_tier1() {
    tier186("t1", "A (K_eff, the lakes, the live sea, the load, the compaction, the recycling)", a186());
}

/// F186 tier 1 with the author's fallback (the live marine step on 2 × 2 blocks).
#[test]
#[ignore]
fn f186_tier1_blocks() {
    tier186("t1b", "A with the marine step on 2 × 2 blocks", Deposition { marine_blocks: true, ..a186() });
}

fn tier186(tag: &str, what: &str, dep: Deposition) {
    let t0 = Instant::now();
    eprintln!("\n==========  F186 tier 1 . 256², the panel, intermediate: without / with {what}  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let (dur, _) = step_durations(T_F166_YEARS, AGES[2].1.unwrap(), N_PHYS);
    let mut rows = Vec::new();
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let s = cached186(&format!("{tag}_seed{seed}"), || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let w = world_oo(seed, p, &m);
            let clk = Clock::start();
            let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
            let c_off = clk.stop();
            let clk = Clock::start();
            let on = physics_acc(seed, p, &m, &w, rough_a, 2, dep);
            let c_on = clk.stop();
            let sr = on.sed.as_ref().unwrap();
            let mut s = format!("   seed {seed}: the physics level without {} · with A {} (the marine step {:.1} s wall, {} sub-steps)\n", fmt_cost(c_off), fmt_cost(c_on), sr.secs_marine, sr.marine_substeps);
            s += &measure184(seed, 2, "off", &off, &dur);
            s += &measure184(seed, 2, "A", &on, &dur);
            let z = zones185(&on, &off.r.classes, &off.r.zm);
            // the slopes and the altitudes in the deposition zone (the cover ≥ 1 m at the end) and in the erosion zone, both runs
            let dep: Vec<bool> = (0..65536).map(|k| sr.thickness_m[k] >= 1.0).collect();
            let ero: Vec<bool> = dep.iter().map(|d| !d).collect();
            let (sd_off, ad_off, _) = slope_alt(&off.r.zm, &dep);
            let (sd_on, ad_on, nd) = slope_alt(&on.r.zm, &dep);
            let (se_off, ae_off, _) = slope_alt(&off.r.zm, &ero);
            let (se_on, ae_on, _) = slope_alt(&on.r.zm, &ero);
            let sh_off = shelf182(&off.r.zm, &a);
            let sh_on = shelf182(&on.r.zm, &a);
            let med = |v: &[(u8, f64, f64, usize)], c: u8| -> f64 {
                let mut x: Vec<f64> = v.iter().filter(|t| t.0 == c).map(|t| t.1).collect();
                if x.is_empty() { f64::NAN } else { q_of(&mut x, 0.5) }
            };
            let b = sr.budget;
            let input = b[0] + b[1];
            s += &format!(
                "      zones: deposition {:.1} % of the land (Δ plains {:+.1} / Δ mountain {:+.1} points), erosion (Δ plains {:+.1} / Δ mountain {:+.1})\n      the deposition zone ({nd} cells): median slope {sd_off:.2}° → {sd_on:.2}°, median altitude {ad_off:.0} → {ad_on:.0} m · the erosion zone: slope {se_off:.2}° → {se_on:.2}°, altitude {ae_off:.0} → {ae_on:.0} m\n      the shelf, passive / active median: {:.1} / {:.1} km → {:.1} / {:.1} km (Harris unconverted 34.9–123.0, ÷ 7.5 4.65–16.4)\n      the sediment: lakes {:.1} % · sea {:.1} % · stranded {:.2} % · subducted {:.1} % · flux control {:.1e} · stock control {:.1e}\n#T {} {} {} {} {sd_off} {sd_on} {ad_off} {ad_on} {se_off} {se_on} {ae_off} {ae_on} {} {} {} {} {} {}\n",
                z[0], z[1], z[2], z[3], z[4], med(&sh_off, 2), med(&sh_off, 1), med(&sh_on, 2), med(&sh_on, 1),
                100.0 * b[3] / input, 100.0 * b[5] / input, 100.0 * b[4] / input, 100.0 * b[7] / input, sr.max_rel_residual, sr.max_stock_residual,
                c_off.0, c_off.1, c_on.0, c_on.1, med(&sh_on, 2), med(&sh_on, 1), b[3] / input, b[5] / input, sr.max_rel_residual, sr.max_stock_residual
            );
            std::fs::create_dir_all(out186()).unwrap();
            alt_png_plain(&off.r.zm).save(out186().join(format!("relief_off_seed{seed}.png"))).unwrap();
            alt_png_plain(&on.r.zm).save(out186().join(format!("relief_{tag}_seed{seed}.png"))).unwrap();
            thick_png(&on.r.zm, &sr.thickness_m, 4000.0).save(out186().join(format!("thickness_{tag}_seed{seed}.png"))).unwrap();
            let (_, _, mo) = lakes182(&off.r.zm);
            let (_, _, mn) = lakes182(&on.r.zm);
            lakes_png(&off.r.zm, &mo).save(out186().join(format!("lakes_off_seed{seed}.png"))).unwrap();
            lakes_png(&on.r.zm, &mn).save(out186().join(format!("lakes_{tag}_seed{seed}.png"))).unwrap();
            s
        });
        eprint!("{}", s.lines().filter(|l| !l.starts_with("#P")).map(|l| format!("{l}\n")).collect::<String>());
        rows.push((seed, s));
    }
    // the summary
    let (mut kept, mut lost, mut e3, mut n) = (0, 0, 0, 0);
    let (mut won, mut cpu, mut dmed, mut dpl, mut dmt) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (_, s) in &rows {
        let w = lines182(s, "W");
        let t = lines182(s, "T").pop();
        if w.len() < 2 {
            continue;
        }
        let (wo, wn) = (&w[0], &w[1]);
        n += 1;
        kept += (wn[0] > 0.0) as usize;
        lost += (wo[1] > 0.0 && wn[1] == 0.0) as usize;
        e3 += (wn[9] < wo[9]) as usize;
        dmed.push(wn[2] - wo[2]);
        dpl.push(100.0 * (wn[3] - wo[3]));
        dmt.push(100.0 * (wn[6] - wo[6]));
        if let Some(t) = t {
            won.push(t[2]);
            cpu.push(t[3]);
        }
    }
    let q = |v: &[f64]| { let mut v = v.to_vec(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
    let mx = |v: &[f64]| v.iter().copied().fold(f64::NAN, f64::max);
    eprintln!(
        "\n   F186 tier 1 ({n} seeds): the physics level with A, wall median {:.1} s, max {:.1} s (the budget 60 s) · the sea kept {kept} / {n} · the circumnavigation lost (against off) {lost} · E3 {e3} / {n} · Δ median altitude median {:+.0} m · Δ plains {:+.1} · Δ mountain {:+.1} points",
        q(&won), mx(&won), q(&dmed), q(&dpl), q(&dmt)
    );
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

use ymir_core::cascade::amplify::Transmit;

/// The amplification chain to 2 048² on a field (m), timed per level; with `bury`, each level's output is buried (the basement
/// is the working field). Returns the 2 048² output, the per-level wall seconds and the burial's own seconds.
fn chain186(field: &GridF32, amp: &AmpConfig, ln_k: &[f32], bury: Option<&Transmit>) -> (GridF32, Vec<f64>, f64) {
    let oc = ocean_mask_t(field, Some([0, 0]));
    let d5: Vec<bool> = (0..field.data.len()).map(|i| !oc[i] && field.data[i] <= 0.0).collect();
    let mut ch = Chain::new(field.clone(), amp.clone()).with_d5(d5).with_ln_k(ln_k.to_vec());
    let mut walls = Vec::new();
    let mut bury_s = 0.0;
    let mut last = GridF32::new(1, 1, 0.0);
    for _ in 0..3 {
        let t = Instant::now();
        ch.next_level(&|| false).unwrap();
        walls.push(t.elapsed().as_secs_f64());
        let r = ch.levels.last().unwrap().result();
        last = match bury {
            Some(tr) => {
                let t = Instant::now();
                let b = tr.bury(&r);
                bury_s += t.elapsed().as_secs_f64();
                b
            }
            None => r,
        };
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    (last, walls, bury_s)
}

/// A crop (2 × half cells) of a 2 048² field centred on a 256² cell.
fn crop186(g: &GridF32, centre: usize, half: usize) -> GridF32 {
    let (cx, cy) = ((centre % 256) * 8, (centre / 256) * 8);
    let mut c = GridF32::new(2 * half, 2 * half, 0.0);
    for y in 0..2 * half {
        for x in 0..2 * half {
            c.data[y * 2 * half + x] = g.data[((cy + 2048 - half + y) % 2048) * 2048 + (cx + 2048 - half + x) % 2048];
        }
    }
    c
}

/// F186 tier 2 -- the chain to 2 048², seed 23, intermediate: without the deposition, and with A (the marine blocks) and the
/// transmission; the stages timed against the 5-minute target (4 096² and 8 192² extrapolated, the selection from Part 0); the
/// plates at 2 048²; the fine cover's volume against the coarse.
#[test]
#[ignore]
fn f186_tier2() {
    let t0 = Instant::now();
    eprintln!("\n==========  F186 tier 2 . the chain to 2 048², seed 23, intermediate: without / with A (blocks) and the transmission  ==========");
    let m = Margin::continent_default();
    let seed = 23u64;
    let p = sel176(seed).unwrap();
    let rough_a = rough_a177();
    let clk = Clock::start();
    let w = world_oo(seed, p, &m);
    let c_c1 = clk.stop();
    let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
    let clk = Clock::start();
    let on = physics_acc(seed, p, &m, &w, rough_a, 2, Deposition { marine_blocks: true, transmit: true, ..a186() });
    let c_phys = clk.stop();
    let sr = on.sed.as_ref().unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let kh = KHetero::declared();
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    // without the deposition
    let (g0, w0, _) = chain186(&off.r.zm, &amp, &ln_k, None);
    // with A and the transmission: the chain on the basement, the outputs buried
    let basement = GridF32 { width: 256, height: 256, data: (0..65536).map(|k| on.r.zm.data[k] - sr.thickness_m[k]).collect() };
    let flow = compute_flow(&on.r.zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let lake: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && flow.filled.data[k] > on.r.zm.data[k]).collect();
    let tr = Transmit { surface_m: on.r.zm.clone(), cover_m: sr.thickness_m.clone(), lake: lake.clone() };
    let (g1, w1, bury_s) = chain186(&basement, &amp, &ln_k, Some(&tr));
    // the cover's volume, fine (the burial's added height over the basement chain's output) against coarse
    let (gb, _, _) = (g1.clone(), 0, 0);
    let (gbase, _, _) = chain186(&basement, &amp, &ln_k, None);
    // the volumes compared like for like: all the cells, and the land only (the fine land: the buried output > 0 m)
    let fa = 0.1953125f64 * 0.1953125 / 1000.0;
    let ca = CELL as f64 * CELL as f64 / 1000.0;
    let fine_all: f64 = gb.data.iter().zip(&gbase.data).map(|(a, b)| (a - b).max(0.0) as f64).sum::<f64>() * fa;
    let fine_land: f64 = gb.data.iter().zip(&gbase.data).filter(|(a, _)| **a > 0.0).map(|(a, b)| (a - b).max(0.0) as f64).sum::<f64>() * fa;
    let coarse_all: f64 = sr.thickness_m.iter().map(|&h| h as f64).sum::<f64>() * ca;
    let coarse_land: f64 = sr.thickness_m.iter().zip(&on.r.zm.data).filter(|(_, z)| **z > 0.0).map(|(h, _)| *h as f64).sum::<f64>() * ca;
    let (fine_vol, coarse_vol) = (fine_land, coarse_land);
    // the macro at 2 048², without and with
    let mm = |g: &GridF32| {
        let cl = classify(&g.data, 2048, 2048, 0.1953125);
        macro_stats(&g.data, &cl, 2048, 2048, 0.1953125)
    };
    let (m0, m1) = (mm(&g0), mm(&g1));
    // the budget: Part 0's selection (on mains) + C1 + the physics level with A + the levels + the burial; 4 096² and 8 192²
    // extrapolated from the 2 048² level with the burial's per-cell cost
    let sel_mains = 108.8;
    let amp_w: f64 = w1.iter().sum();
    let per2 = w1[2] / (2048.0 * 2048.0);
    let bury_per = bury_s / ((512.0f64 * 512.0) + (1024.0 * 1024.0) + (2048.0 * 2048.0));
    let ext = |n: f64, scale: f64| per2 * n * n * scale + bury_per * n * n;
    // the 2 048² level's cost per cell with its budget; 4 096² and 8 192² with Part 0's per-operator extrapolation ratio (154 / 41.7, 615 / 41.7)
    let (e4, e8) = (w1[2] * 154.0 / 41.7 + bury_per * 4096.0 * 4096.0, w1[2] * 615.0 / 41.7 + bury_per * 8192.0 * 8192.0);
    let _ = ext;
    let total = sel_mains + c_c1.0 + c_phys.0 + amp_w + bury_s + e4 + e8;
    let mut out = format!(
        "   the stages: the selection {sel_mains:.0} s (Part 0, mains) · C1 {:.1} s · the physics level with A {:.1} s · the amplifications 512² / 1 024² / 2 048² {:.1} / {:.1} / {:.1} s (without the deposition {:.1} / {:.1} / {:.1} s) · the burial {bury_s:.2} s · 4 096² ~{e4:.0} s and 8 192² ~{e8:.0} s (extrapolated, Part 0's ratios) · **the chain to 8 192² ≈ {total:.0} s against 300 s**\n",
        c_c1.0, c_phys.0, w1[0], w1[1], w1[2], w0[0], w0[1], w0[2]
    );
    out += &format!(
        "   at 2 048²: plains {:.1} → {:.1} % · plateau {:.1} → {:.1} % · hill {:.1} → {:.1} % · mountain {:.1} → {:.1} % · median {:.0} → {:.0} m · the cover's volume on land, fine {fine_vol:.0} km³ against coarse {coarse_vol:.0} km³ ({:+.1} %); on every cell, fine {fine_all:.0} against coarse {coarse_all:.0} km³ ({:+.1} %)\n",
        100.0 * m0.fractions[0], 100.0 * m1.fractions[0], 100.0 * m0.fractions[1], 100.0 * m1.fractions[1], 100.0 * m0.fractions[2], 100.0 * m1.fractions[2], 100.0 * m0.fractions[3], 100.0 * m1.fractions[3],
        m0.deciles[4], m1.deciles[4], 100.0 * (fine_vol / coarse_vol - 1.0), 100.0 * (fine_all / coarse_all - 1.0)
    );
    // the plates and the crops: a filled lake (a covered coarse lake), a delta (land gained at the coast), a plateau
    std::fs::create_dir_all(out186()).unwrap();
    alt_png(&g0, 4).save(out186().join("map2048_off_seed23.png")).unwrap();
    alt_png(&g1, 4).save(out186().join("map2048_A_seed23.png")).unwrap();
    let comp_centre = |mask: &[bool]| -> Option<usize> { comps8(mask, 256).into_iter().max_by_key(|c| c.len()).map(|c| c[c.len() / 2]) };
    let filled_lake: Vec<bool> = (0..65536).map(|k| lake[k] && sr.thickness_m[k] > 0.0).collect();
    let delta: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && off.r.zm.data[k] <= 0.0).collect();
    let cl256 = &on.r.classes;
    let plateau: Vec<bool> = (0..65536).map(|k| cl256[k] == PLATEAU).collect();
    for (nm, mask) in [("lake", &filled_lake), ("delta", &delta), ("plateau", &plateau)] {
        if let Some(c) = comp_centre(mask) {
            alt_png(&crop186(&g0, c, 160), 1).save(out186().join(format!("crop2048_{nm}_off_seed23.png"))).unwrap();
            alt_png(&crop186(&g1, c, 160), 1).save(out186().join(format!("crop2048_{nm}_A_seed23.png"))).unwrap();
            out += &format!("   crop « {nm} » centred on 256² cell {c}\n");
        } else {
            out += &format!("   crop « {nm} »: no such cell\n");
        }
    }
    out += &format!("   tier 2: {:.0} s\n", t0.elapsed().as_secs_f64());
    std::fs::create_dir_all(results186()).unwrap();
    std::fs::write(results186().join("t2_seed23.txt"), &out).unwrap();
    eprint!("{out}");
}

/// F186 step 2 -- tier 1 with the blocks and the author's prolongation (bilinear, capped by the room, the excess to the ring).
#[test]
#[ignore]
fn f186_tier1_blocks_bilinear() {
    tier186("t1d", "A with the marine step on 2 × 2 blocks, the bilinear prolongation", Deposition { marine_blocks: true, ..a186() });
}

/// F186 step 2 -- tier 1 with the blocks, the gains bilinear (the author's) and the losses per block in proportion to the sediment.
#[test]
#[ignore]
fn f186_tier1_blocks_final() {
    tier186("t1e", "A with the marine step on 2 × 2 blocks, the gains bilinear, the losses per block", Deposition { marine_blocks: true, ..a186() });
}

/// F186 step 2 -- tier 1 per cell with the active-cell sweep (bit for bit the full sweep): the physics level's time.
#[test]
#[ignore]
fn f186_tier1_active() {
    tier186("t1c", "A per cell, the marine step sweeping the active cells only", a186());
}

/// F186 step 4's retained configuration: A with the marine step on 2 × 2 blocks and the bilinear prolongation.
fn a186b() -> Deposition {
    Deposition { marine_blocks: true, ..a186() }
}

/// F186 step 4 -- the measures for the adoption, the panel at 256² (resumable per seed): at the intermediate age without the
/// deposition, with A (blocks), and with A without the load (E3); at the old age without and with A.
#[test]
#[ignore]
fn f186_step4() {
    let t0 = Instant::now();
    eprintln!("\n==========  F186 step 4 . 256², the panel: (v), W10, E3 without the load, the old age  ==========");
    let m = Margin::continent_default();
    let rough_a = rough_a177();
    let runs: [(&str, usize, Option<Deposition>); 5] = [
        ("off_100", 2, None),
        ("A_100", 2, Some(a186b())),
        ("A_noload_100", 2, Some(Deposition { load: None, ..a186b() })),
        ("off_300", 3, None),
        ("A_300", 3, Some(a186b())),
    ];
    for seed in PANEL {
        let Some(p) = sel176(seed) else { continue };
        let s = cached186(&format!("s4_seed{seed}"), || {
            let (o, _) = observe179(seed, p, &m);
            let info = step_info(&o.st, &o.kin, o.ds_max, true);
            let a = a1(&info, &o.kin, &o.tracker);
            let w = world_oo(seed, p, &m);
            let mut s = format!("   seed {seed}\n");
            for (nm, ai, dep) in runs {
                let clk = Clock::start();
                let r = match dep {
                    None => physics_dep(seed, p, &m, &w, rough_a, ai, None, None),
                    Some(d) => physics_acc(seed, p, &m, &w, rough_a, ai, d),
                };
                let c = clk.stop();
                let (line, mask) = measure182(nm, &r, &o, &a, &w);
                s += &format!("   [{nm}] {}\n{line}", fmt_cost(c));
                std::fs::create_dir_all(out186()).unwrap();
                if nm != "off_100" {
                    alt_png_plain(&r.r.zm).save(out186().join(format!("relief_s4_{nm}_seed{seed}.png"))).unwrap();
                    lakes_png(&r.r.zm, &mask).save(out186().join(format!("lakes_s4_{nm}_seed{seed}.png"))).unwrap();
                }
            }
            s
        });
        eprint!("{}", s.lines().filter(|l| !l.starts_with('#')).map(|l| format!("{l}\n")).collect::<String>());
    }
    summary186_step4();
    eprintln!("   step 4: {:.0} s", t0.elapsed().as_secs_f64());
}

/// The per-run data lines of a step-4 file, in the run order (off_100, A_100, A_noload_100, off_300, A_300).
fn s4_rows(s: &str, key: &str) -> Vec<Vec<f64>> {
    lines182(s, key)
}

fn summary186_step4() {
    let q = |v: &[f64]| { let mut v: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect(); if v.is_empty() { f64::NAN } else { q_of(&mut v, 0.5) } };
    let mx = |v: &[f64]| v.iter().copied().filter(|x| x.is_finite()).fold(f64::NAN, f64::max);
    let mut out = String::from("\n   F186 step 4, the panel's summary (medians over the seeds; the run order off_100, A_100, A_noload_100, off_300, A_300)\n");
    let (mut v_ok, mut n) = ([0usize; 5], 0usize);
    let mut w: [[Vec<f64>; 3]; 5] = Default::default();
    let mut dw: [Vec<f64>; 3] = Default::default();
    let (mut e3, mut e3_nl) = (0usize, 0usize);
    let (mut la_off, mut la_a, mut la_nl) = (Vec::new(), Vec::new(), Vec::new());
    let (mut sea, mut circ_off, mut circ_on, mut circ_lost) = ([Vec::new(), Vec::new()], 0usize, 0usize, 0usize);
    let mut old: [[Vec<f64>; 4]; 2] = Default::default();
    let mut per_seed = String::new();
    for seed in PANEL {
        let Ok(s) = std::fs::read_to_string(results186().join(format!("s4_seed{seed}.txt"))) else { continue };
        let (f, ww, l, c, mm) = (s4_rows(&s, "F"), s4_rows(&s, "W"), s4_rows(&s, "L"), s4_rows(&s, "C"), s4_rows(&s, "M"));
        if ww.len() < 5 || f.len() < 5 || l.len() < 5 || c.len() < 5 || mm.len() < 5 {
            continue;
        }
        n += 1;
        for i in 0..5 {
            v_ok[i] += (f[i][8] == 0.0) as usize;
            for j in 0..3 {
                w[i][j].push(ww[i][j]);
            }
        }
        for j in 0..3 {
            dw[j].push((ww[1][j] - ww[0][j]).abs());
        }
        let large = |r: &Vec<f64>| r[7];
        e3 += (large(&l[1]) < large(&l[0])) as usize;
        e3_nl += (large(&l[2]) < large(&l[0])) as usize;
        la_off.push(large(&l[0]));
        la_a.push(large(&l[1]));
        la_nl.push(large(&l[2]));
        for (i, ri) in [3usize, 4].iter().enumerate() {
            sea[i].push(100.0 * (1.0 - c[*ri][3] / 65536.0));
            // the macro: plain, mountain, median; the large lakes' count
            old[i][0].push(100.0 * mm[*ri][0]);
            old[i][1].push(100.0 * mm[*ri][3]);
            old[i][2].push(mm[*ri][5]);
            old[i][3].push(l[*ri][3]);
        }
        circ_off += (c[3][2] > 0.0) as usize;
        circ_on += (c[4][2] > 0.0) as usize;
        circ_lost += (c[3][2] > 0.0 && c[4][2] == 0.0) as usize;
        per_seed += &format!(
            "      seed {seed}: (v) {} / {} / {} · W10 off {:.3} ({:.3} / {:.3}) → A {:.3} ({:.3} / {:.3}) · the large lakes' area off {:.0} → A {:.0}, without the load {:.0} km² · old: sea {:.1} → {:.1} %, circumnavigable {} → {}, plains {:.1} → {:.1} %, mountain {:.1} → {:.1} %, median {:.0} → {:.0} m, large lakes {} → {}\n",
            (f[0][8] == 0.0) as u8, (f[1][8] == 0.0) as u8, (f[2][8] == 0.0) as u8,
            ww[0][0], ww[0][1], ww[0][2], ww[1][0], ww[1][1], ww[1][2], large(&l[0]), large(&l[1]), large(&l[2]),
            100.0 * (1.0 - c[3][3] / 65536.0), 100.0 * (1.0 - c[4][3] / 65536.0), c[3][2] > 0.0, c[4][2] > 0.0,
            100.0 * mm[3][0], 100.0 * mm[4][0], 100.0 * mm[3][3], 100.0 * mm[4][3], mm[3][5], mm[4][5], l[3][3], l[4][3]
        );
    }
    out += &per_seed;
    out += &format!(
        "   intermediate ({n} seeds): (v) holds off {} / A {} / A without the load {} · W10 global / active / passive: off {:.3} / {:.3} / {:.3} → A {:.3} / {:.3} / {:.3} (|Δ| max {:.3} / {:.3} / {:.3})\n",
        v_ok[0], v_ok[1], v_ok[2], q(&w[0][0]), q(&w[0][1]), q(&w[0][2]), q(&w[1][0]), q(&w[1][1]), q(&w[1][2]), mx(&dw[0]), mx(&dw[1]), mx(&dw[2])
    );
    out += &format!(
        "   E3 (the large lakes' area falls against off): A {e3} / {n}, A without the load {e3_nl} / {n} · medians off {:.0} · A {:.0} · without the load {:.0} km²\n",
        q(&la_off), q(&la_a), q(&la_nl)
    );
    out += &format!(
        "   old ({n} seeds): (v) off {} / A {} · the sea {:.1} → {:.1} % · circumnavigable off {circ_off} / A {circ_on} (lost against the same world without the deposition: {circ_lost}) · plains {:.1} → {:.1} % · mountain {:.1} → {:.1} % · median {:.0} → {:.0} m · large lakes {:.0} → {:.0}\n",
        v_ok[3], v_ok[4], q(&sea[0]), q(&sea[1]), q(&old[0][0]), q(&old[1][0]), q(&old[0][1]), q(&old[1][1]), q(&old[0][2]), q(&old[1][2]), q(&old[0][3]), q(&old[1][3])
    );
    std::fs::write(results186().join("s4_summary.txt"), &out).unwrap();
    eprint!("{out}");
}

/// F186 step 5 -- tier 2 again with step 2's marine (the blocks, the bilinear prolongation) and step 3's burial (the level L):
/// the chain to 2 048², seed 23, intermediate, without the deposition and with A and the transmission; the stages timed against the 5-minute target (4 096² and 8 192² extrapolated, the selection from Part 0); the
/// plates at 2 048²; the fine cover's volume against the coarse.
#[test]
#[ignore]
fn f186_step5() {
    let t0 = Instant::now();
    eprintln!("\n==========  F186 step 5 . the chain to 2 048², seed 23, intermediate: without / with A (blocks, bilinear) and the burial by the level L  ==========");
    let m = Margin::continent_default();
    let seed = 23u64;
    let p = sel176(seed).unwrap();
    let rough_a = rough_a177();
    let clk = Clock::start();
    let w = world_oo(seed, p, &m);
    let c_c1 = clk.stop();
    let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
    let clk = Clock::start();
    let on = physics_acc(seed, p, &m, &w, rough_a, 2, Deposition { transmit: true, ..a186b() });
    let c_phys = clk.stop();
    let sr = on.sed.as_ref().unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let kh = KHetero::declared();
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    // without the deposition
    let (g0, w0, _) = chain186(&off.r.zm, &amp, &ln_k, None);
    // with A and the transmission: the chain on the basement, the outputs buried
    let basement = GridF32 { width: 256, height: 256, data: (0..65536).map(|k| on.r.zm.data[k] - sr.thickness_m[k]).collect() };
    let flow = compute_flow(&on.r.zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let lake: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && flow.filled.data[k] > on.r.zm.data[k]).collect();
    let tr = Transmit { surface_m: on.r.zm.clone(), cover_m: sr.thickness_m.clone(), lake: lake.clone() };
    let (g1, w1, bury_s) = chain186(&basement, &amp, &ln_k, Some(&tr));
    // the cover's volume, fine (the burial's added height over the basement chain's output) against coarse
    let (gb, _, _) = (g1.clone(), 0, 0);
    let (gbase, _, _) = chain186(&basement, &amp, &ln_k, None);
    // the volumes compared like for like: all the cells, and the land only (the fine land: the buried output > 0 m)
    let fa = 0.1953125f64 * 0.1953125 / 1000.0;
    let ca = CELL as f64 * CELL as f64 / 1000.0;
    let fine_all: f64 = gb.data.iter().zip(&gbase.data).map(|(a, b)| (a - b).max(0.0) as f64).sum::<f64>() * fa;
    let fine_land: f64 = gb.data.iter().zip(&gbase.data).filter(|(a, _)| **a > 0.0).map(|(a, b)| (a - b).max(0.0) as f64).sum::<f64>() * fa;
    let coarse_all: f64 = sr.thickness_m.iter().map(|&h| h as f64).sum::<f64>() * ca;
    let coarse_land: f64 = sr.thickness_m.iter().zip(&on.r.zm.data).filter(|(_, z)| **z > 0.0).map(|(h, _)| *h as f64).sum::<f64>() * ca;
    let (fine_vol, coarse_vol) = (fine_land, coarse_land);
    // the macro at 2 048², without and with
    let mm = |g: &GridF32| {
        let cl = classify(&g.data, 2048, 2048, 0.1953125);
        macro_stats(&g.data, &cl, 2048, 2048, 0.1953125)
    };
    let (m0, m1) = (mm(&g0), mm(&g1));
    // the budget: Part 0's selection (on mains) + C1 + the physics level with A + the levels + the burial; 4 096² and 8 192²
    // extrapolated from the 2 048² level with the burial's per-cell cost
    let sel_mains = 108.8;
    let amp_w: f64 = w1.iter().sum();
    let per2 = w1[2] / (2048.0 * 2048.0);
    let bury_per = bury_s / ((512.0f64 * 512.0) + (1024.0 * 1024.0) + (2048.0 * 2048.0));
    let ext = |n: f64, scale: f64| per2 * n * n * scale + bury_per * n * n;
    // the 2 048² level's cost per cell with its budget; 4 096² and 8 192² with Part 0's per-operator extrapolation ratio (154 / 41.7, 615 / 41.7)
    let (e4, e8) = (w1[2] * 154.0 / 41.7 + bury_per * 4096.0 * 4096.0, w1[2] * 615.0 / 41.7 + bury_per * 8192.0 * 8192.0);
    let _ = ext;
    let total = sel_mains + c_c1.0 + c_phys.0 + amp_w + bury_s + e4 + e8;
    let mut out = format!(
        "   the stages: the selection {sel_mains:.0} s (Part 0, mains) · C1 {:.1} s · the physics level with A {:.1} s · the amplifications 512² / 1 024² / 2 048² {:.1} / {:.1} / {:.1} s (without the deposition {:.1} / {:.1} / {:.1} s) · the burial {bury_s:.2} s · 4 096² ~{e4:.0} s and 8 192² ~{e8:.0} s (extrapolated, Part 0's ratios) · **the chain to 8 192² ≈ {total:.0} s against 300 s**\n",
        c_c1.0, c_phys.0, w1[0], w1[1], w1[2], w0[0], w0[1], w0[2]
    );
    out += &format!(
        "   at 2 048²: plains {:.1} → {:.1} % · plateau {:.1} → {:.1} % · hill {:.1} → {:.1} % · mountain {:.1} → {:.1} % · median {:.0} → {:.0} m · the cover's volume on land, fine {fine_vol:.0} km³ against coarse {coarse_vol:.0} km³ ({:+.1} %); on every cell, fine {fine_all:.0} against coarse {coarse_all:.0} km³ ({:+.1} %)\n",
        100.0 * m0.fractions[0], 100.0 * m1.fractions[0], 100.0 * m0.fractions[1], 100.0 * m1.fractions[1], 100.0 * m0.fractions[2], 100.0 * m1.fractions[2], 100.0 * m0.fractions[3], 100.0 * m1.fractions[3],
        m0.deciles[4], m1.deciles[4], 100.0 * (fine_vol / coarse_vol - 1.0), 100.0 * (fine_all / coarse_all - 1.0)
    );
    // the plates and the crops: a filled lake (a covered coarse lake), a delta (land gained at the coast), a plateau
    std::fs::create_dir_all(out186()).unwrap();
    alt_png(&g0, 4).save(out186().join("s5_map2048_off_seed23.png")).unwrap();
    alt_png(&g1, 4).save(out186().join("s5_map2048_A_seed23.png")).unwrap();
    let comp_centre = |mask: &[bool]| -> Option<usize> { comps8(mask, 256).into_iter().max_by_key(|c| c.len()).map(|c| c[c.len() / 2]) };
    let filled_lake: Vec<bool> = (0..65536).map(|k| lake[k] && sr.thickness_m[k] > 0.0).collect();
    let delta: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && off.r.zm.data[k] <= 0.0).collect();
    let cl256 = &on.r.classes;
    let plateau: Vec<bool> = (0..65536).map(|k| cl256[k] == PLATEAU).collect();
    for (nm, mask) in [("lake", &filled_lake), ("delta", &delta), ("plateau", &plateau)] {
        if let Some(c) = comp_centre(mask) {
            alt_png(&crop186(&g0, c, 160), 1).save(out186().join(format!("s5_crop2048_{nm}_off_seed23.png"))).unwrap();
            alt_png(&crop186(&g1, c, 160), 1).save(out186().join(format!("s5_crop2048_{nm}_A_seed23.png"))).unwrap();
            if nm == "delta" {
                // the stripes' sources, separately: the physics level's surface bicubic at 2 048² (no amplification), and the
                // basement's chain without the burial
                let bic = upsample_to(&on.r.zm, 2048);
                alt_png(&crop186(&bic, c, 160), 1).save(out186().join("s5_crop2048_delta_bicubic_seed23.png")).unwrap();
                alt_png(&crop186(&gbase, c, 160), 1).save(out186().join("s5_crop2048_delta_basement_seed23.png")).unwrap();
                // the 2-cell oscillation at 256² around the delta, along x and along y (mean |z(i+1) − 2 z(i) + z(i−1)|, m) on the
                // sea cells within 20 cells: the surface with A, its sediment thickness, and the world without the deposition
                let (cx, cy) = ((c % 256) as i64, (c / 256) as i64);
                let osc = |g: &[f32], sea_of: &GridF32| -> (f64, f64) {
                    let at = |x: i64, y: i64| g[(y.rem_euclid(256) as usize) * 256 + x.rem_euclid(256) as usize] as f64;
                    let (mut ax, mut ay, mut nn) = (0f64, 0f64, 0usize);
                    for y in cy - 20..cy + 20 {
                        for x in cx - 20..cx + 20 {
                            if sea_of.data[(y.rem_euclid(256) as usize) * 256 + x.rem_euclid(256) as usize] > 0.0 {
                                continue;
                            }
                            ax += (at(x + 1, y) - 2.0 * at(x, y) + at(x - 1, y)).abs();
                            ay += (at(x, y + 1) - 2.0 * at(x, y) + at(x, y - 1)).abs();
                            nn += 1;
                        }
                    }
                    (ax / nn.max(1) as f64, ay / nn.max(1) as f64)
                };
                let (sx, sy) = osc(&on.r.zm.data, &on.r.zm);
                let (hx, hy) = osc(&sr.thickness_m, &on.r.zm);
                let (ox, oy) = osc(&off.r.zm.data, &off.r.zm);
                out += &format!("   the delta's 2-cell oscillation at 256² (sea cells within 20 cells; along x / along y, m): the surface with A {sx:.1} / {sy:.1} · its thickness {hx:.1} / {hy:.1} · the world without {ox:.1} / {oy:.1}\n");
            }
            out += &format!("   crop « {nm} » centred on 256² cell {c}\n");
        } else {
            out += &format!("   crop « {nm} »: no such cell\n");
        }
    }
    out += &format!("   step 5: {:.0} s\n", t0.elapsed().as_secs_f64());
    std::fs::create_dir_all(results186()).unwrap();
    std::fs::write(results186().join("s5_seed23.txt"), &out).unwrap();
    eprint!("{out}");
}

/// The deficit lines' index of a marine thickness: the deep-sea cells (z < −200 m) whose marine thickness is below 0.8 × the
/// median of their 8 neighbours (the neighbours' median > 50 m).
fn lines186(zm: &GridF32, t: &[f32]) -> usize {
    let mut n = 0;
    for k in 0..65536 {
        if zm.data[k] >= -200.0 {
            continue;
        }
        let (x, y) = ((k % 256) as i64, (k / 256) as i64);
        let mut v: Vec<f64> = Vec::with_capacity(8);
        for (dx, dy) in [(-1i64, -1i64), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)] {
            v.push(t[((y + dy).rem_euclid(256) as usize) * 256 + (x + dx).rem_euclid(256) as usize] as f64);
        }
        let m = q_of(&mut v, 0.5);
        if m > 50.0 && (t[k] as f64) < 0.8 * m {
            n += 1;
        }
    }
    n
}

/// F186 step 2's diagnosis -- the deficit lines of the blocks' prolongation, seed 23 intermediate: the blocks as they are, the
/// blocks without the recycling; the per-cell run's index from its own run for reference.
#[test]
#[ignore]
fn f186_lines() {
    let m = Margin::continent_default();
    let seed = 23u64;
    let p = sel176(seed).unwrap();
    let rough_a = rough_a177();
    let w = world_oo(seed, p, &m);
    let mut out = String::new();
    for (nm, d) in [
        ("cells_recycle_bilinear", a186()),
    ] {
        let r = physics_acc(seed, p, &m, &w, rough_a, 2, d);
        let sr = r.sed.as_ref().unwrap();
        let li = lines186(&r.r.zm, &sr.marine_m);
        out += &format!("   {nm}: deficit-line cells {li}\n");
        std::fs::create_dir_all(out186()).unwrap();
        thick_png(&r.r.zm, &sr.marine_m, 4000.0).save(out186().join(format!("lines_marine_{nm}_seed23.png"))).unwrap();
        eprint!("{out}");
    }
    std::fs::write(results186().join("lines_seed23_recycle_bilinear.txt"), &out).unwrap();
}

/// A grey PNG of a field crop (black lo, white hi).
fn grey_png(g: &GridF32, lo: f32, hi: f32) -> image::RgbaImage {
    let (w, h) = (g.width as u32, g.height as u32);
    let mut img = image::RgbaImage::new(w, h);
    for y in 0..h {
        for x in 0..w {
            let v = ((g.data[(y * w + x) as usize] - lo) / (hi - lo)).clamp(0.0, 1.0);
            let c = (255.0 * v) as u8;
            img.put_pixel(x, h - 1 - y, image::Rgba([c, c, c, 255]));
        }
    }
    img
}

/// The new high bumps of a world with the deposition against without (256², land, z_on − z_off > 150 m), by 8-connected
/// group: the group's size, its highest cell's altitude without / with, the cover H, its marine / lake / fluvial parts, the
/// basement's change (z_on − H − z_off: the incision's part, K_eff or the load), the lake flag.
fn bumps186(on: &Dep182, off: &Dep182, lake: &[bool]) -> String {
    let sr = on.sed.as_ref().unwrap();
    let m: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && on.r.zm.data[k] - off.r.zm.data[k] > 150.0).collect();
    let mut gs = comps8(&m, 256);
    gs.sort_by_key(|g| std::cmp::Reverse(g.len()));
    let mut s = format!("      {} groups with z_on − z_off > 150 m on land\n", gs.len());
    for g in gs.iter().take(12) {
        let k = *g.iter().max_by(|a, b| (on.r.zm.data[**a] - off.r.zm.data[**a]).total_cmp(&(on.r.zm.data[**b] - off.r.zm.data[**b]))).unwrap();
        let h = sr.thickness_m[k];
        let base = on.r.zm.data[k] - h - off.r.zm.data[k];
        let nl = g.iter().filter(|&&q| lake[q]).count();
        let nm = g.iter().filter(|&&q| sr.marine_m[q] > 1.0).count();
        s += &format!(
            "      group of {} cells, peak at 256² cell {k} (x {}, y {}): z {:.0} → {:.0} m (Δ {:+.0}) · cover H {h:.0} m (marine {:.0}, lake {:.0}, fluvial {:.0}) · basement Δ {base:+.0} m · cells in a lake {nl}, with a marine deposit > 1 m {nm}\n",
            g.len(), k % 256, k / 256, off.r.zm.data[k], on.r.zm.data[k], on.r.zm.data[k] - off.r.zm.data[k], sr.marine_m[k], sr.lake_m[k], sr.fluvial_m[k]
        );
    }
    s
}

/// F186 continuation, step 2 -- the diagnosis (no correction): the new white bumps with the deposition and the delta crop's
/// vertical edge, seed 23, intermediate, at 256² and on the chain to 2 048² (step 5's configuration, the blocks; the per-cell
/// marine at 256² for comparison).
#[test]
#[ignore]
fn f186_bumps() {
    let t0 = Instant::now();
    let m = Margin::continent_default();
    let seed = 23u64;
    let p = sel176(seed).unwrap();
    let rough_a = rough_a177();
    let w = world_oo(seed, p, &m);
    let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
    let mut out = String::from("\n==========  F186 continuation, step 2: the bumps and the delta's edge, seed 23  ==========\n");
    let mut kept = None;
    for (nm, d) in [("blocks", Deposition { transmit: true, ..a186b() }), ("cells", Deposition { transmit: true, ..a186() })] {
        let on = physics_acc(seed, p, &m, &w, rough_a, 2, d);
        let flow = compute_flow(&on.r.zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
        let lake: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && flow.filled.data[k] > on.r.zm.data[k]).collect();
        out += &format!("   {nm} (256²):\n{}", bumps186(&on, &off, &lake));
        if nm == "blocks" {
            kept = Some((on, lake));
        }
        eprint!("{out}");
    }
    // the chain (blocks, step 5's): the burial's added height at 2 048² and its parts, on the bumps and around the delta's edge
    let (on, lake) = kept.unwrap();
    let sr = on.sed.as_ref().unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let kh = KHetero::declared();
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    let basement = GridF32 { width: 256, height: 256, data: (0..65536).map(|k| on.r.zm.data[k] - sr.thickness_m[k]).collect() };
    let tr = Transmit { surface_m: on.r.zm.clone(), cover_m: sr.thickness_m.clone(), lake: lake.clone() };
    let (g1, _, _) = chain186(&basement, &amp, &ln_k, Some(&tr));
    let (gb, _, _) = chain186(&basement, &amp, &ln_k, None);
    let (g0, _, _) = chain186(&off.r.zm, &amp, &ln_k, None);
    // the bumps at 2 048²: per 256² peak cell of the 12 largest groups, the fine max without / basement chain / buried
    let mk: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && on.r.zm.data[k] - off.r.zm.data[k] > 150.0).collect();
    let mut gs = comps8(&mk, 256);
    gs.sort_by_key(|g| std::cmp::Reverse(g.len()));
    out += "   blocks at 2 048² (per group, the fine maximum over its cells): without the deposition / the basement chain / buried\n";
    for g in gs.iter().take(12) {
        let fmax = |f: &GridF32| -> f32 {
            let mut v = f32::MIN;
            for &k in g {
                let (cx, cy) = ((k % 256) * 8, (k / 256) * 8);
                for dy in 0..8 {
                    for dx in 0..8 {
                        v = v.max(f.data[(cy + dy) * 2048 + cx + dx]);
                    }
                }
            }
            v
        };
        out += &format!("      group of {} cells at 256² cell {}: {:.0} / {:.0} / {:.0} m\n", g.len(), g[0], fmax(&g0), fmax(&gb), fmax(&g1));
    }
    // the delta's crop (step 5's centre, 256² cell 8731): the burial's added height, the covered stencil, and the edge
    let c = 8731usize;
    let add = GridF32 { width: 2048, height: 2048, data: g1.data.iter().zip(&gb.data).map(|(a, b)| a - b).collect() };
    let cov = GridF32 {
        width: 2048,
        height: 2048,
        data: (0..2048 * 2048).map(|k| if sr.thickness_m[((k / 2048) / 8) * 256 + (k % 2048) / 8] > 0.0 { 1.0 } else { 0.0 }).collect(),
    };
    std::fs::create_dir_all(out186()).unwrap();
    grey_png(&crop186(&add, c, 160), 0.0, 300.0).save(out186().join("diag_delta_burial_added_seed23.png")).unwrap();
    grey_png(&crop186(&cov, c, 160), 0.0, 1.0).save(out186().join("diag_delta_cover_mask_seed23.png")).unwrap();
    alt_png(&crop186(&g1, c, 160), 1).save(out186().join("diag_delta_buried_seed23.png")).unwrap();
    alt_png(&crop186(&gb, c, 160), 1).save(out186().join("diag_delta_basement_seed23.png")).unwrap();
    // the edge: per crop column, the mean added height over the crop's land rows; the largest jump between neighbour columns
    let cr = crop186(&add, c, 160);
    let cz = crop186(&g1, c, 160);
    let colm: Vec<f64> = (0..320)
        .map(|x| {
            let v: Vec<f64> = (0..320).filter(|&y| cz.data[y * 320 + x] > 0.0).map(|y| cr.data[y * 320 + x] as f64).collect();
            if v.is_empty() { 0.0 } else { v.iter().sum::<f64>() / v.len() as f64 }
        })
        .collect();
    let (jx, jv) = (1..320).map(|x| (x, colm[x] - colm[x - 1])).max_by(|a, b| a.1.abs().total_cmp(&b.1.abs())).unwrap();
    let fine_x0 = ((c % 256) * 8 + 2048 - 160) % 2048;
    out += &format!(
        "   the delta crop: the burial's added height, column means on land, the largest jump {jv:+.0} m at crop column {jx} (fine x {}, 256² column {}, {} the 256² cell edge)\n",
        fine_x0 + jx, (fine_x0 + jx) / 8, if (fine_x0 + jx) % 8 == 0 { "on" } else { "off" }
    );
    out += &format!("   step 2's diagnosis: {:.0} s\n", t0.elapsed().as_secs_f64());
    std::fs::write(results186().join("diag_bumps_edge_seed23.txt"), &out).unwrap();
    eprint!("{out}");
}

/// F186 continuation, step 2 -- the new white bumps located at 2 048² (step 5's configuration, the blocks): the fine cells white
/// with the deposition (≥ 1 200 m, alt_png's white ramp from 1 200 m) and below 600 m without, grouped (8-connected, on the
/// 512² map grid as drawn: x right, y down); per group its map position, its peak with / without / the basement chain, and its
/// 256² parent's state (z with / without, the cover, the marine and lake parts, the basement's change against the worlds
/// without the load and without K_eff).
#[test]
#[ignore]
fn f186_white() {
    let t0 = Instant::now();
    let m = Margin::continent_default();
    let seed = 23u64;
    let p = sel176(seed).unwrap();
    let rough_a = rough_a177();
    let w = world_oo(seed, p, &m);
    let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
    let on = physics_acc(seed, p, &m, &w, rough_a, 2, Deposition { transmit: true, ..a186b() });
    let noload = physics_acc(seed, p, &m, &w, rough_a, 2, Deposition { transmit: true, load: None, ..a186b() });
    let nokeff = physics_acc(seed, p, &m, &w, rough_a, 2, Deposition { transmit: true, fluvial_keff: None, ..a186b() });
    let sr = on.sed.as_ref().unwrap();
    let flow = compute_flow(&on.r.zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let lake: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && flow.filled.data[k] > on.r.zm.data[k]).collect();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let kh = KHetero::declared();
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    let basement = GridF32 { width: 256, height: 256, data: (0..65536).map(|k| on.r.zm.data[k] - sr.thickness_m[k]).collect() };
    let tr = Transmit { surface_m: on.r.zm.clone(), cover_m: sr.thickness_m.clone(), lake: lake.clone() };
    let (g1, _, _) = chain186(&basement, &amp, &ln_k, Some(&tr));
    let (gb, _, _) = chain186(&basement, &amp, &ln_k, None);
    let (g0, _, _) = chain186(&off.r.zm, &amp, &ln_k, None);
    // the basement's interior cells below the sea at 256² (the chain lifts them, amplify.rs:903)
    let oc = ocean_mask_t(&basement, Some([0, 0]));
    let mut out = String::from("\n==========  F186 continuation, step 2: the new white bumps at 2 048², seed 23 (blocks)  ==========\n");
    let mut white = vec![false; 512 * 512];
    for y in 0..512 {
        for x in 0..512 {
            let k = (y * 4) * 2048 + x * 4;
            white[y * 512 + x] = g1.data[k] >= 1200.0 && g0.data[k] < 600.0;
        }
    }
    let mut gs = comps8(&white, 512);
    gs.sort_by_key(|g| std::cmp::Reverse(g.len()));
    out += &format!("   {} groups (map pixels, 512²)\n", gs.len());
    for g in gs.iter().take(15) {
        let pk = *g.iter().max_by(|a, b| g1.data[(**a / 512) * 4 * 2048 + (**a % 512) * 4].total_cmp(&g1.data[(**b / 512) * 4 * 2048 + (**b % 512) * 4])).unwrap();
        let (mx, my) = (pk % 512, pk / 512);
        let fk = (my * 4) * 2048 + mx * 4;
        let c = (my / 2) * 256 + mx / 2;
        let (xs, ys): (Vec<usize>, Vec<usize>) = g.iter().map(|&q| (q % 512, q / 512)).unzip();
        out += &format!(
            "      {} px, map x {}–{} y {}–{} (peak at x {mx} y {my}; 256² cell {c}): fine {:.0} m with, {:.0} without, the basement chain {:.0} · 256²: z {:.0} → {:.0} m, cover {:.0} m (marine {:.0}, lake {:.0}), in a lake {}, the basement {:.0} m below the sea and interior {} · the basement's change against without: {:+.0} m (without the load {:+.0}, without K_eff {:+.0})\n",
            g.len(), xs.iter().min().unwrap(), xs.iter().max().unwrap(), ys.iter().min().unwrap(), ys.iter().max().unwrap(),
            g1.data[fk], g0.data[fk], gb.data[fk],
            off.r.zm.data[c], on.r.zm.data[c], sr.thickness_m[c], sr.marine_m[c], sr.lake_m[c], lake[c],
            -basement.data[c], basement.data[c] <= 0.0 && !oc[c],
            basement.data[c] - off.r.zm.data[c],
            noload.r.zm.data[c] - noload.sed.as_ref().unwrap().thickness_m[c] - off.r.zm.data[c],
            nokeff.r.zm.data[c] - nokeff.sed.as_ref().unwrap().thickness_m[c] - off.r.zm.data[c]
        );
    }
    out += &format!("   {:.0} s\n", t0.elapsed().as_secs_f64());
    std::fs::write(results186().join("diag_white_seed23.txt"), &out).unwrap();
    eprint!("{out}");
}

/// F186 continuation, step 3 -- the chain to 2 048² with the per-cell marine step (the physical reference), the recycling
/// interpolated, and the burial by the level L:
/// the chain to 2 048², seed 23, intermediate, without the deposition and with A and the transmission; the stages timed against the 5-minute target (4 096² and 8 192² extrapolated, the selection from Part 0); the
/// plates at 2 048²; the fine cover's volume against the coarse.
#[test]
#[ignore]
fn f186_chain_cells() {
    let t0 = Instant::now();
    eprintln!("\n==========  F186 continuation 3 . the chain to 2 048², seed 23, intermediate: without / with A (the marine per cell) and the burial by the level L  ==========");
    let m = Margin::continent_default();
    let seed = 23u64;
    let p = sel176(seed).unwrap();
    let rough_a = rough_a177();
    let clk = Clock::start();
    let w = world_oo(seed, p, &m);
    let c_c1 = clk.stop();
    let off = physics_dep(seed, p, &m, &w, rough_a, 2, None, None);
    let clk = Clock::start();
    let on = physics_acc(seed, p, &m, &w, rough_a, 2, Deposition { transmit: true, ..a186() });
    let c_phys = clk.stop();
    let sr = on.sed.as_ref().unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let kh = KHetero::declared();
    let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
    let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
    // without the deposition
    let (g0, w0, _) = chain186(&off.r.zm, &amp, &ln_k, None);
    // with A and the transmission: the chain on the basement, the outputs buried
    let basement = GridF32 { width: 256, height: 256, data: (0..65536).map(|k| on.r.zm.data[k] - sr.thickness_m[k]).collect() };
    let flow = compute_flow(&on.r.zm, &FlowConfig { sea_level: 0.0, ..Default::default() });
    let lake: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && flow.filled.data[k] > on.r.zm.data[k]).collect();
    let tr = Transmit { surface_m: on.r.zm.clone(), cover_m: sr.thickness_m.clone(), lake: lake.clone() };
    let (g1, w1, bury_s) = chain186(&basement, &amp, &ln_k, Some(&tr));
    // the cover's volume, fine (the burial's added height over the basement chain's output) against coarse
    let (gb, _, _) = (g1.clone(), 0, 0);
    let (gbase, _, _) = chain186(&basement, &amp, &ln_k, None);
    // the volumes compared like for like: all the cells, and the land only (the fine land: the buried output > 0 m)
    let fa = 0.1953125f64 * 0.1953125 / 1000.0;
    let ca = CELL as f64 * CELL as f64 / 1000.0;
    let fine_all: f64 = gb.data.iter().zip(&gbase.data).map(|(a, b)| (a - b).max(0.0) as f64).sum::<f64>() * fa;
    let fine_land: f64 = gb.data.iter().zip(&gbase.data).filter(|(a, _)| **a > 0.0).map(|(a, b)| (a - b).max(0.0) as f64).sum::<f64>() * fa;
    let coarse_all: f64 = sr.thickness_m.iter().map(|&h| h as f64).sum::<f64>() * ca;
    let coarse_land: f64 = sr.thickness_m.iter().zip(&on.r.zm.data).filter(|(_, z)| **z > 0.0).map(|(h, _)| *h as f64).sum::<f64>() * ca;
    let (fine_vol, coarse_vol) = (fine_land, coarse_land);
    // the macro at 2 048², without and with
    let mm = |g: &GridF32| {
        let cl = classify(&g.data, 2048, 2048, 0.1953125);
        macro_stats(&g.data, &cl, 2048, 2048, 0.1953125)
    };
    let (m0, m1) = (mm(&g0), mm(&g1));
    // the budget: Part 0's selection (on mains) + C1 + the physics level with A + the levels + the burial; 4 096² and 8 192²
    // extrapolated from the 2 048² level with the burial's per-cell cost
    let sel_mains = 108.8;
    let amp_w: f64 = w1.iter().sum();
    let per2 = w1[2] / (2048.0 * 2048.0);
    let bury_per = bury_s / ((512.0f64 * 512.0) + (1024.0 * 1024.0) + (2048.0 * 2048.0));
    let ext = |n: f64, scale: f64| per2 * n * n * scale + bury_per * n * n;
    // the 2 048² level's cost per cell with its budget; 4 096² and 8 192² with Part 0's per-operator extrapolation ratio (154 / 41.7, 615 / 41.7)
    let (e4, e8) = (w1[2] * 154.0 / 41.7 + bury_per * 4096.0 * 4096.0, w1[2] * 615.0 / 41.7 + bury_per * 8192.0 * 8192.0);
    let _ = ext;
    let total = sel_mains + c_c1.0 + c_phys.0 + amp_w + bury_s + e4 + e8;
    let mut out = format!(
        "   the stages: the selection {sel_mains:.0} s (Part 0, mains) · C1 {:.1} s · the physics level with A {:.1} s · the amplifications 512² / 1 024² / 2 048² {:.1} / {:.1} / {:.1} s (without the deposition {:.1} / {:.1} / {:.1} s) · the burial {bury_s:.2} s · 4 096² ~{e4:.0} s and 8 192² ~{e8:.0} s (extrapolated, Part 0's ratios) · **the chain to 8 192² ≈ {total:.0} s against 300 s**\n",
        c_c1.0, c_phys.0, w1[0], w1[1], w1[2], w0[0], w0[1], w0[2]
    );
    out += &format!(
        "   at 2 048²: plains {:.1} → {:.1} % · plateau {:.1} → {:.1} % · hill {:.1} → {:.1} % · mountain {:.1} → {:.1} % · median {:.0} → {:.0} m · the cover's volume on land, fine {fine_vol:.0} km³ against coarse {coarse_vol:.0} km³ ({:+.1} %); on every cell, fine {fine_all:.0} against coarse {coarse_all:.0} km³ ({:+.1} %)\n",
        100.0 * m0.fractions[0], 100.0 * m1.fractions[0], 100.0 * m0.fractions[1], 100.0 * m1.fractions[1], 100.0 * m0.fractions[2], 100.0 * m1.fractions[2], 100.0 * m0.fractions[3], 100.0 * m1.fractions[3],
        m0.deciles[4], m1.deciles[4], 100.0 * (fine_vol / coarse_vol - 1.0), 100.0 * (fine_all / coarse_all - 1.0)
    );
    // the plates and the crops: a filled lake (a covered coarse lake), a delta (land gained at the coast), a plateau
    std::fs::create_dir_all(out186()).unwrap();
    alt_png(&g0, 4).save(out186().join("c3_map2048_off_seed23.png")).unwrap();
    alt_png(&g1, 4).save(out186().join("c3_map2048_A_seed23.png")).unwrap();
    let comp_centre = |mask: &[bool]| -> Option<usize> { comps8(mask, 256).into_iter().max_by_key(|c| c.len()).map(|c| c[c.len() / 2]) };
    let filled_lake: Vec<bool> = (0..65536).map(|k| lake[k] && sr.thickness_m[k] > 0.0).collect();
    let delta: Vec<bool> = (0..65536).map(|k| on.r.zm.data[k] > 0.0 && off.r.zm.data[k] <= 0.0).collect();
    let cl256 = &on.r.classes;
    let plateau: Vec<bool> = (0..65536).map(|k| cl256[k] == PLATEAU).collect();
    for (nm, mask) in [("lake", &filled_lake), ("delta", &delta), ("plateau", &plateau)] {
        if let Some(c) = comp_centre(mask) {
            alt_png(&crop186(&g0, c, 160), 1).save(out186().join(format!("c3_crop2048_{nm}_off_seed23.png"))).unwrap();
            alt_png(&crop186(&g1, c, 160), 1).save(out186().join(format!("c3_crop2048_{nm}_A_seed23.png"))).unwrap();
            if nm == "delta" {
                // the stripes' sources, separately: the physics level's surface bicubic at 2 048² (no amplification), and the
                // basement's chain without the burial
                let bic = upsample_to(&on.r.zm, 2048);
                alt_png(&crop186(&bic, c, 160), 1).save(out186().join("c3_crop2048_delta_bicubic_seed23.png")).unwrap();
                alt_png(&crop186(&gbase, c, 160), 1).save(out186().join("c3_crop2048_delta_basement_seed23.png")).unwrap();
                // the 2-cell oscillation at 256² around the delta, along x and along y (mean |z(i+1) − 2 z(i) + z(i−1)|, m) on the
                // sea cells within 20 cells: the surface with A, its sediment thickness, and the world without the deposition
                let (cx, cy) = ((c % 256) as i64, (c / 256) as i64);
                let osc = |g: &[f32], sea_of: &GridF32| -> (f64, f64) {
                    let at = |x: i64, y: i64| g[(y.rem_euclid(256) as usize) * 256 + x.rem_euclid(256) as usize] as f64;
                    let (mut ax, mut ay, mut nn) = (0f64, 0f64, 0usize);
                    for y in cy - 20..cy + 20 {
                        for x in cx - 20..cx + 20 {
                            if sea_of.data[(y.rem_euclid(256) as usize) * 256 + x.rem_euclid(256) as usize] > 0.0 {
                                continue;
                            }
                            ax += (at(x + 1, y) - 2.0 * at(x, y) + at(x - 1, y)).abs();
                            ay += (at(x, y + 1) - 2.0 * at(x, y) + at(x, y - 1)).abs();
                            nn += 1;
                        }
                    }
                    (ax / nn.max(1) as f64, ay / nn.max(1) as f64)
                };
                let (sx, sy) = osc(&on.r.zm.data, &on.r.zm);
                let (hx, hy) = osc(&sr.thickness_m, &on.r.zm);
                let (ox, oy) = osc(&off.r.zm.data, &off.r.zm);
                out += &format!("   the delta's 2-cell oscillation at 256² (sea cells within 20 cells; along x / along y, m): the surface with A {sx:.1} / {sy:.1} · its thickness {hx:.1} / {hy:.1} · the world without {ox:.1} / {oy:.1}\n");
            }
            out += &format!("   crop « {nm} » centred on 256² cell {c}\n");
        } else {
            out += &format!("   crop « {nm} »: no such cell\n");
        }
    }
    out += &format!("   continuation 3: {:.0} s\n", t0.elapsed().as_secs_f64());
    std::fs::create_dir_all(results186()).unwrap();
    std::fs::write(results186().join("c3_seed23.txt"), &out).unwrap();
    eprint!("{out}");
}
