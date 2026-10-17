//! ADR Finding 176 -- the accretion merges a pair only at collision (continent-continent edges), behind the continent profile's
//! flag, with option 5 and the seaward margin; the panel 8, 18, 19, 23, 29 plus 10 screened seeds; by tiers, resumable (the results
//! in `docs/reports/relief_method/f176_accretion/results/`). F175's helpers are kept. Declared in
//! `docs/reports/relief_method/f176_accretion/f176_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f176_cascade -- --ignored --exact f176_tier1_c1 --nocapture
//!      … f176_tier2_screen, f176_tier3_physics, f176_tier4_cascade
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
    let src = run_c1_sources_mapped_flags(&mut st, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &SteinSteinParams::default(), &[C1Term::Erosion], Mapping::Airy(Airy::declared()), m.oo, m.merge_at_collision, m.arcs_emerge);
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

/// Tier 1 -- C1 only, the panel at the fixed continent: F175 | F176 | F176 + B.
#[test]
#[ignore]
fn f176_tier1_c1() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out176()).unwrap();
    eprintln!("\n==========  F176 tier 1 . C1 only, the panel at the fixed continent (F174 B′'s plates)  ==========");
    for seed in PANEL {
        let p = fixed_profile(seed);
        for (name, m) in variants176() {
            let key = key176("t1", seed, name);
            let s = cached176(&key, || c1_line176(seed, p, name, &m, &key));
            eprint!("{s}");
        }
    }
    eprintln!("   tier 1: {:.0} s", t0.elapsed().as_secs_f64());
}

/// Tier 2 -- the screening: the 15 seeds, the selection redone once under F176's rule (parallel), C1 on the selected continent.
#[test]
#[ignore]
fn f176_tier2_screen() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out176()).unwrap();
    eprintln!("\n==========  F176 tier 2 . the screening: 15 seeds, the selection redone under F176  ==========");
    let seeds: Vec<u64> = PANEL.iter().chain(SCREEN.iter()).copied().collect();
    let m = Margin::seaward_oo_collision();
    std::thread::scope(|sc| {
        for &seed in &seeds {
            if results176().join(format!("t2_seed{seed}.txt")).exists() {
                continue;
            }
            sc.spawn(move || {
                cached176(&format!("t2_seed{seed}"), || {
                    let clk = Clock::start();
                    let (np, tried) = select_continent_default_margin(seed, &m);
                    let c = clk.stop();
                    match np {
                        Some(np) => {
                            let s = &tried.last().unwrap().1;
                            let cd = s.candidates[s.chosen.unwrap()];
                            let p = cd.profile();
                            let line = c1_line176(seed, p, "F176 (reselected)", &m, &format!("t2_seed{seed}"));
                            format!("   seed {seed}: selection N {np}, land {:.1} %, r* {:?} · {:.0} s wall\n{line}", 100.0 * cd.land, cd.circ.sea_r, c.0)
                        }
                        None => format!("   seed {seed}: selection NONE · {:.0} s wall\n", c.0),
                    }
                });
            });
        }
    });
    let (mut fails, mut circ_lost, mut over70) = (0, 0, 0);
    let (mut conts, mut conts_and) = (0usize, 0usize);
    let mut actives = Vec::new();
    for &seed in &seeds {
        let s = std::fs::read_to_string(results176().join(format!("t2_seed{seed}.txt"))).unwrap();
        eprint!("{s}");
        if s.contains("NONE") {
            fails += 1;
            continue;
        }
        // parse the numbers back from the line
        let grab = |key: &str| -> Option<f32> {
            let i = s.find(key)? + key.len();
            s[i..].split(|c: char| !(c.is_ascii_digit() || c == '.')).find(|x| !x.is_empty())?.parse().ok()
        };
        if s.contains("circumnavigable false") {
            circ_lost += 1;
        }
        let act = grab("% → end ").unwrap_or(f32::NAN);
        actives.push(act);
        if act > 70.0 {
            over70 += 1;
        }
        conts += grab("continents (≥ 20 cells) ").unwrap_or(0.0) as usize;
        conts_and += grab("of which with ≥ 1 Andean margin ").unwrap_or(0.0) as usize;
    }
    actives.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if actives.is_empty() { f32::NAN } else { actives[actives.len() / 2] };
    eprintln!(
        "\n   R: the selection fails on {fails} / 15 (stop at ≥ 3: {}) · circumnavigation lost on {circ_lost} (stop at ≥ 2: {}) · active > 70 % on {over70} (stop at ≥ 3: {}) · continents with ≥ 1 Andean margin {conts_and} / {conts} ({:.0} %; success ≥ 70 %) · the active share at the end: median {med:.1} % (success 10–50 %), sorted {:?} · tier 2: {:.0} s",
        fails >= 3, circ_lost >= 2, over70 >= 3, 100.0 * conts_and as f32 / conts.max(1) as f32, actives, t0.elapsed().as_secs_f64()
    );
}

/// Tier 3 -- the physics level at 256², the panel at the fixed continent: F175 | F176 young and old, F176 + B young.
#[test]
#[ignore]
fn f176_tier3_physics() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out176()).unwrap();
    eprintln!("\n==========  F176 tier 3 . the physics level at 256², the fixed continent  ==========");
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let kh = KHetero::declared();
    for seed in PANEL {
        let p = fixed_profile(seed);
        for (name, m) in variants176() {
            let key = key176("t3", seed, name);
            let ages: Vec<usize> = if name == "F176 + B" { vec![1] } else { vec![1, 3] };
            let s = cached176(&key, || {
                let w = world_oo(seed, p, &m);
                let ser = c1_series(seed, p, &m);
                let cover = cover_map_m(seed, p, &m, COVER_M);
                let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
                let rock = Rock { rho_k: RHO_K, cover_m: cover, ln_k: Some(ln_k) };
                let cfg = base_cfg(seed, rough_a, Some([0, 0]));
                let mut out = String::new();
                for &ai in &ages {
                    let (aname, age) = AGES[ai];
                    let (r, _) = hist_k2n(&w.snaps, &cfg, age, &rock, N_PHYS);
                    let land = land_t(&r.zm);
                    let (wi, win, rw) = wall_index(&r.zm);
                    let (wa, wp, na, np) = w10_by_type(&r.zm, &ser.map_end);
                    let (isl, isl_cells) = islands_n(&land, 256, 1);
                    let hy = hydrology(&r.zm, CELL);
                    let mm = &r.m;
                    if ai == 1 {
                        alt_png_plain(&r.zm).save(out176().join(format!("map256_young_seed{seed}_{}.png", name.replace(' ', "_").replace('+', "p")))).unwrap();
                    }
                    out += &format!(
                        "   seed {seed} {name} {aname}: W10 {wi:.3} (W_in {win:.3}, rampart {rw}) · active / passive {wa:.3} / {wp:.3} ({na} / {np} cells) · islands (256²) {isl} ({:.0} km²) · lakes {} ({:.0} km²) · macro plain {:.1} · plateau {:.1} · hill {:.1} · mountain {:.1} % · land {:.0} km² · median {:.0} m · physics CPU {:.0} s\n",
                        isl_cells as f32 * CELL * CELL, hy.lakes.len(), hy.lakes.iter().map(|l| l.area_km2).sum::<f32>(), 100.0 * mm.fractions[0], 100.0 * mm.fractions[1], 100.0 * mm.fractions[2], 100.0 * mm.fractions[3], mm.land_km2, mm.deciles.get(4).copied().unwrap_or(f32::NAN), r.cost.1
                    );
                }
                out
            });
            eprint!("{s}");
        }
    }
    // the maps F175 | F176 at young, one row per seed
    let mut rows: Vec<Vec<image::RgbaImage>> = Vec::new();
    for seed in PANEL {
        let mut row = Vec::new();
        for name in ["F175", "F176"] {
            if let Ok(i) = image::open(out176().join(format!("map256_young_seed{seed}_{name}.png"))) {
                row.push(i.to_rgba8());
            }
        }
        rows.push(row);
    }
    let (tw, th) = (260u32, 260u32);
    let mut img = image::RgbaImage::from_pixel(2 * tw, rows.len() as u32 * th, image::Rgba([255, 255, 255, 255]));
    for (ri, r) in rows.iter().enumerate() {
        for (ci, t) in r.iter().enumerate() {
            image::imageops::overlay(&mut img, t, (ci as u32 * tw + 2) as i64, (ri as u32 * th + 2) as i64);
        }
    }
    img.save(out176().join("maps256_young_F175_vs_F176.png")).unwrap();
    // the margins at the end, F176, the panel
    let mut mrow = Vec::new();
    for seed in PANEL {
        if let Ok(i) = image::open(out176().join(format!("margins_end_{}.png", key176("t1", seed, "F176")))) {
            mrow.push(i.to_rgba8());
        }
    }
    let mut mimg = image::RgbaImage::from_pixel(mrow.len() as u32 * 388, 388, image::Rgba([255, 255, 255, 255]));
    for (i, t) in mrow.iter().enumerate() {
        image::imageops::overlay(&mut mimg, t, (i * 388 + 2) as i64, 2);
    }
    mimg.save(out176().join("margins_end_F176_panel.png")).unwrap();
    eprintln!("   tier 3: {:.0} s", t0.elapsed().as_secs_f64());
}

/// Tier 4 -- the cascade to 2 048² on the panel world with the most Andean cells under F176 (young): the transect and the map.
#[test]
#[ignore]
fn f176_tier4_cascade() {
    let t0 = Instant::now();
    std::fs::create_dir_all(out176()).unwrap();
    eprintln!("\n==========  F176 tier 4 . the cascade to 2 048², one F176 world with an Andean margin, young  ==========");
    let m = Margin::seaward_oo_collision();
    let andean = |seed: u64| c1_series(seed, fixed_profile(seed), &m).andean_cells;
    let seed = PANEL.iter().copied().max_by_key(|&s| (andean(s), std::cmp::Reverse(s))).unwrap();
    let p = fixed_profile(seed);
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let kh = KHetero::declared();
    let s = cached176(&format!("t4_seed{seed}"), || {
        let w = world_oo(seed, p, &m);
        let ser = c1_series(seed, p, &m);
        let cover = cover_map_m(seed, p, &m, COVER_M);
        let ln_k = kh.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]));
        let rock = Rock { rho_k: RHO_K, cover_m: cover.clone(), ln_k: Some(ln_k.clone()) };
        let cfg = base_cfg(seed, rough_a, Some([0, 0]));
        let amp = AmpConfig { k_hetero: Some(kh), ..frozen(&pi_table, Some([0, 0])) };
        let (r, _) = hist_k2n(&w.snaps, &cfg, AGES[1].1, &rock, N_PHYS);
        let (g, _, cc) = chain_low(&r.zm, amp, Some(ln_k), 0.0);
        let fl = land_t(&g);
        let kt = coast_torus(&fl, 2048, 2048, 32);
        let (isl, isl_km2) = islands_off_main(&fl, 2048);
        alt_png(&g, 4).save(out176().join(format!("map2048_young_seed{seed}.png"))).unwrap();
        let act = (0..4096).find(|&k| ser.map_end[k] == 3 && ser.map_t0[k] == 2).or_else(|| (0..4096).find(|&k| ser.map_end[k] == 3));
        let pas = (0..4096).find(|&k| ser.map_end[k] == 2);
        let s64 = w.st.s.data().to_vec();
        let pt = w.st.plate_type.data().to_vec();
        let (mut zs, mut ss) = (Vec::new(), Vec::new());
        if let Some(k) = pas {
            let (a, b) = transect(&g, &s64, &pt, k);
            zs.push((a, [30, 140, 60]));
            ss.push((b, [30, 140, 60]));
        }
        if let Some(k) = act {
            let (a, b) = transect(&g, &s64, &pt, k);
            zs.push((a, [220, 30, 30]));
            ss.push((b, [220, 30, 30]));
        }
        transect_png(&zs, &ss, &out176().join(format!("transect_seed{seed}.png")));
        format!(
            "   seed {seed} F176 young (2 048²): the coast {} · islands ≥ 1 km² {isl} ({isl_km2:.0} km²) · transect: passive {} · Andean {} · chain CPU {:.0} s\n",
            fm(&meas(&kt)), pas.is_some(), act.is_some(), cc.1
        )
    });
    eprint!("{s}");
    eprintln!("   tier 4 (seed {seed}): {:.0} s", t0.elapsed().as_secs_f64());
}

/// F176-B diagnostic (C1 only): the arc load on oceanic cells under F176 -- its maximum, p99, and the cells above the 0.715 needed
/// to emerge under Airy (0.2 + arc ≥ 0.915).
#[test]
#[ignore]
fn f176_arc_load() {
    let m = Margin::seaward_oo_collision();
    for seed in PANEL {
        let p = fixed_profile(seed);
        let ser = c1_series(seed, p, &m);
        let mut v: Vec<f64> = (0..4096).filter(|&k| ser.map_end[k] == 0 && ser.sub_sum[k] > 0.0).map(|k| ser.sub_sum[k]).collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |p: f64| if v.is_empty() { 0.0 } else { v[((v.len() - 1) as f64 * p) as usize] };
        let line = format!(
            "   seed {seed} F176: oceanic cells with a net arc deposit {} · net arc per cell p50 / p99 / max {:.3} / {:.3} / {:.3} (s units; 0.715 to emerge) · above 0.715: {}\n",
            v.len(), q(0.5), q(0.99), v.last().copied().unwrap_or(0.0), v.iter().filter(|&&x| x >= 0.715).count()
        );
        let s = cached176(&format!("arcload_seed{seed}"), || line);
        eprint!("{s}");
    }
}
