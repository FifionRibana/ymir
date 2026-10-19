//! ADR Finding 171 -- a rock that depends on depth: a soft cover on a hard basement (K1), then the heterogeneous K (K2),
//! with 16 plates by default (falling back to 20 then 24), on the torus.
//! Declared in `docs/reports/relief_method/f171_rock/f171_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f171_cascade -- --ignored --exact f171_rock --nocapture

#![allow(dead_code)]

mod common;

// F181: the margin off explicitly (the continent profile's default margin changed at F181)
use ymir_core::cascade::margin::Margin as MarginOff181;
#[allow(unused_imports)]
use ymir_core::cascade::profile::{select_continent_default_margin as _f181_sel, select_continent_plates_margin as _f181_plates};
use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, PSEED, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::restrict;
use ymir_core::cascade::amplify::{AmpConfig, Chain, KHetero, Variant, Warp, WarpPlace, ocean_mask, ocean_mask_t};
use ymir_core::cascade::history::{
    Airy, C1Sources, C1Term, HistoryRun, Mapping, T_F166_YEARS, Timing, physics_history, physics_history_rock, run_c1_sources_mapped,
    step_durations,
};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{
    CLASS_NAMES, Circumnavigation, Coast, MOUNTAIN, MacroStats, PLAIN, PLATEAU, RELIEF_RADIUS_KM, circumnavigation, classify, coast,
    coast_torus, edt, hovius, local_relief, macro_stats,
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
    root().join("docs/reports/relief_method/f171_rock/images")
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

#[test]
#[ignore]
fn f171_rock() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let k = n2m();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 171 . a rock that depends on depth (K1), the heterogeneous K (K2)  ==========");
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
    let c_coast = coast(&land_w(&cg), 1024, 1024, 32);
    let cp391 = predictability(&c512, 200.0 / 512.0, None);
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let ckm = DOMAIN_KM / 2048.0;
    let ctrl = facet_control(2048, ckm);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    eprintln!(
        "   Corse 195 m: mountain slope p50/p90 {:.3}/{:.3} · coast {:?} · 391 m: CV_λ {:.3} · CV_L {:.3} · σ_θ {:.1}° · R2_g {:.3}",
        c_mtn.r.slope_p50, c_mtn.r.slope_p90, c_coast, cp391.cv_lambda, cp391.cv_len, cp391.theta_std, cp391.r2_gully
    );

    // ── A: the plates by default (16, falling back to 20 then 24) ──
    eprintln!("\n   A1 — the default selection (16 plates, then 20, then 24; 50–62 % land, r* ≥ 1):");
    let mut chosen: Vec<Option<(C1Profile, usize)>> = Vec::new();
    for seed in SEEDS {
        let clk = Clock::start();
        let (np, tried) = _f181_sel(seed, &MarginOff181::off());
        for (n, s) in &tried {
            let best = s.candidates.iter().filter(|c| (0.50..=0.62).contains(&c.land) && !c.circ.wraps).filter_map(|c| c.circ.sea_r).max();
            eprintln!(
                "      seed {seed} N {n}: {} candidates ({} distinct runs) · best r* in range {:?} · {}",
                s.candidates.len(),
                s.runs,
                best,
                match s.chosen {
                    Some(i) => { let c = &s.candidates[i]; format!("CHOSEN k {} {} start {} land {:.1} % r* {}", c.k, if c.compact { "compact" } else { "BFS" }, c.start, 100.0 * c.land, c.circ.sea_r.unwrap()) }
                    None => "none qualifies".into(),
                }
            );
        }
        eprintln!("      seed {seed}: retained N {np:?} · the selection's cost {}", fmt_cost(clk.stop()));
        chosen.push(np.map(|n| {
            let s = &tried.last().unwrap().1;
            (s.candidates[s.chosen.unwrap()].profile(), n)
        }));
    }
    let n_sel = chosen.iter().filter(|c| c.is_some()).count();
    let need = if n_sel == 4 { 3 } else { 2 };
    if n_sel < 3 {
        eprintln!("\n   R: fewer than 3 seeds selected → STOP");
        return;
    }

    // ── the runs ──
    struct Res {
        tag: String,
        zm: GridF32,
        classes: Vec<u8>,
        m: MacroStats,
        g: GridF32,
        rock: Vec<u8>,
        cost: f64,
        peak: f32,
        slope_ok: bool,
        mtn_slope: (f32, f32),
        coast: Coast,
        knick: (usize, f32, Vec<(usize, f32)>),
        reg: String,
        crop_tex: Option<GridF32>,
    }
    let kh = KHetero::declared();
    let measure = |tag: &str, seed: u64, w: &World, cover: &[f32], age: Option<f64>, rock: Option<Rock>, khet: Option<KHetero>| -> Res {
        let cfg = base_cfg(seed, rough_a, Some([0, 0]));
        let (r, e_cum) = hist_rock(&w.snaps, &cfg, age, rock.as_ref());
        let ln_k = khet.map(|kk| kk.band(256, CELL, 0, kk.resolved(CELL), PSEED, Some([0, 0])));
        let amp = AmpConfig { k_hetero: khet, ..frozen(&pi_table, Some([0, 0])) };
        let (g, d5, c) = chain_k(&r.zm, &r.d5, amp, ln_k);
        let land = land_t(&r.zm);
        let rk = rock_at_surface(&land, cover, &e_cum);
        let knick = knickpoints(&r.zm, &rk);
        let m8 = up_u8(&r.classes, 256, 8);
        let mtn8: Vec<bool> = m8.iter().map(|&x| x == MOUNTAIN).collect();
        let excl: Vec<bool> = (0..2048 * 2048).map(|i| d5[i] || !mtn8[i]).collect();
        let whole_t = tex(&g, ckm, &d5);
        let mtn_t = tex(&g, ckm, &excl);
        let mut lines = Vec::new();
        for (t, x, cc) in [("whole", &whole_t, &c_whole), ("mountain", &mtn_t, &c_mtn)] {
            let el = elements(x, cc, ckm, ctrl_ok);
            let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
            lines.push(format!("2048² {t}: {} of {} within ×1.5 · out: {}", el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") }));
        }
        let slope_ok = !out15(mtn_t.r.slope_p50, c_mtn.r.slope_p50) && !out15(mtn_t.r.slope_p90, c_mtn.r.slope_p90);
        let (rd_, fl) = flat_map(&g, ckm, Some(&d5));
        let per: Vec<String> = (1..5u8)
            .map(|cl| {
                let cells = (0..2048 * 2048).filter(|&i| rd_[i] && m8[i] == cl).count();
                let f = (0..2048 * 2048).filter(|&i| rd_[i] && fl[i] && m8[i] == cl).count();
                format!("{} {:.2} %", CLASS_NAMES[cl as usize], 100.0 * f as f32 / cells.max(1) as f32)
            })
            .collect();
        let p1024 = {
            let gg = restrict(&g);
            let dd = (0..1024 * 1024).map(|i| d5[(i / 1024 * 2) * 2048 + (i % 1024) * 2]).collect::<Vec<_>>();
            predictability(&gg, 2.0 * ckm, Some(&dd))
        };
        let hv: Vec<String> = [10.0f32, 25.0].iter().map(|&a| { let h = hovius(&g, &mtn8, ckm, a); format!("A ≥ {a}: R {:.2} ({} basins)", h.ratio, h.outlets) }).collect();
        let reg = format!(
            "391 m: CV_λ {:.3} · CV_L {:.3} · σ_θ {:.1}° · R2_g {:.3} | 195 m: σ_θ {:.1}° · R2_g {:.3} | Hovius {}",
            p1024.cv_lambda, p1024.cv_len, p1024.theta_std, p1024.r2_gully, whole_t.p.theta_std, whole_t.p.r2_gully, hv.join(" · ")
        );
        let nl = land.iter().filter(|&&b| b).count().max(1) as f32;
        let share = |v: u8| 100.0 * rk.iter().filter(|&&x| x == v).count() as f32 / nl;
        let pb: Vec<usize> = (0..65536).filter(|&i| rk[i] == 2 && r.classes[i] == PLAIN).collect();
        let hy = hydrology(&r.zm, CELL);
        let kc = coast_torus(&land_t(&g), 2048, 2048, 32);
        let cost = w.cost.1 + r.cost.1 + c.1;
        eprintln!("   ── seed {seed} {tag}:");
        print_macro(&format!("{tag} classes"), &r.m);
        eprintln!(
            "      peaks 256² {:.0} m · 2048² {:.0} m · the rock at the surface: cover {:.1} % · basement {:.1} % · shield {:.1} % · plains on the basement {} cells, altitude {:.0} m · knickpoints {} (median height {:.0} m) · lakes {} / {:.0} km² · E_cum p50 / p90 {:.0} / {:.0} m",
            peak_m(&r.zm),
            peak_m(&g),
            share(1),
            share(2),
            share(3),
            pb.len(),
            pb.iter().map(|&i| r.zm.data[i]).sum::<f32>() / pb.len().max(1) as f32,
            knick.0,
            knick.1,
            hy.lakes.len(),
            hy.lake_id.iter().filter(|&&l| l != 0).count() as f32 * CELL * CELL,
            { let v = sorted((0..65536).filter(|&i| land[i] && !e_cum.is_empty()).map(|i| e_cum[i]).collect()); pct(&v, 0.5) },
            { let v = sorted((0..65536).filter(|&i| land[i] && !e_cum.is_empty()).map(|i| e_cum[i]).collect()); pct(&v, 0.9) }
        );
        for l in &lines {
            eprintln!("      {l}");
        }
        eprintln!("      slope guard: mountain-class slope p50/p90 {:.3}/{:.3} (Corsica {:.3}/{:.3}) → {}", mtn_t.r.slope_p50, mtn_t.r.slope_p90, c_mtn.r.slope_p50, c_mtn.r.slope_p90, if slope_ok { "ok" } else { "FAILS" });
        eprintln!("      F by class: {}", per.join(" · "));
        eprintln!("      regularity: {reg}");
        eprintln!("      coast (torus) {kc:?} · cost {cost:.0} s CPU per world (physics {:.0} s, chain {:.0} s)", r.cost.1, c.1);
        Res {
            tag: tag.to_string(),
            peak: peak_m(&g),
            classes: r.classes.clone(),
            m: r.m.clone(),
            zm: r.zm,
            g,
            rock: rk,
            cost,
            slope_ok,
            mtn_slope: (mtn_t.r.slope_p50, mtn_t.r.slope_p90),
            coast: kc,
            knick,
            reg,
            crop_tex: None,
        }
    };

    struct SeedRes {
        seed: u64,
        off: [usize; 2],
        runs: Vec<Res>,
    }
    let mut all: Vec<SeedRes> = Vec::new();
    let (mut succ_k1, mut succ_k2, mut pene, mut offscale) = (0usize, 0usize, 0usize, 0usize);
    for (si, seed) in SEEDS.iter().copied().enumerate() {
        let Some((p, np)) = chosen[si] else {
            eprintln!("\n   seed {seed}: not selected even at 24 plates, skipped");
            continue;
        };
        eprintln!("\n   ══ seed {seed} ({p:?}, N {np}) ══");
        let w = world(seed, p);
        let off = viz_offset(w.snaps.last().unwrap());
        let land64: Vec<bool> = w.snaps.last().unwrap().data.iter().map(|&v| v > 0.5).collect();
        let c64 = circumnavigation(&land64, 64);
        eprintln!("      land 64² {:.1} % · {}", 100.0 * land64.iter().filter(|&&b| b).count() as f32 / 4096.0, fmt_circ(&c64, 6.25));
        let cover = cover_map(seed, p, COVER_M);
        let rock1 = Rock { rho_k: RHO_K, cover_m: cover.clone(), ln_k: None };
        let rock2 = |sigma: f32| Rock { ln_k: Some(KHetero { sigma, ..kh }.band(256, CELL, 0, kh.resolved(CELL), PSEED, Some([0, 0]))), ..rock1.clone() };
        let mut runs = Vec::new();
        let old = Some(300.0e6);
        runs.push(measure("K0 old", seed, &w, &cover, old, None, None));
        if seed == PSEED {
            for (name, age) in AGES.iter().take(3) {
                runs.push(measure(&format!("K0 {name}"), seed, &w, &cover, *age, None, None));
            }
        }
        for (name, age) in AGES {
            runs.push(measure(&format!("K1 {name}"), seed, &w, &cover, age, Some(rock1.clone()), None));
        }
        for (name, age) in AGES {
            runs.push(measure(&format!("K2 {name}"), seed, &w, &cover, age, Some(rock2(kh.sigma)), Some(kh)));
        }
        if seed == PSEED {
            // the sensitivity, the témoin, old profile
            for rho in [3.0f32, 30.0] {
                runs.push(measure(&format!("K1 old ρ_K {rho}"), seed, &w, &cover, old, Some(Rock { rho_k: rho, ..rock1.clone() }), None));
            }
            for dc in [750.0f32, 3000.0] {
                let cv = cover_map(seed, p, dc);
                runs.push(measure(&format!("K1 old d_c {dc:.0} m"), seed, &w, &cv, old, Some(Rock { cover_m: cv.clone(), ..rock1.clone() }), None));
            }
            for sg in [0.25f32, 1.0] {
                let k2 = KHetero { sigma: sg, ..kh };
                runs.push(measure(&format!("K2 old σ {sg}"), seed, &w, &cover, old, Some(rock2(sg)), Some(k2)));
            }
        }
        // the summaries
        let get = |t: &str| runs.iter().find(|r| r.tag == t).unwrap();
        for v in ["K1", "K2"] {
            eprintln!(
                "   seed {seed} {v}: the mountain share against A: {}",
                AGES.iter().map(|(n, _)| format!("{n} {:.1} %", 100.0 * get(&format!("{v} {n}")).m.fractions[3])).collect::<Vec<_>>().join(" · ")
            );
        }
        let k0 = get("K0 old");
        for v in ["K1", "K2"] {
            let r = get(&format!("{v} old 300 Myr"));
            let m = &r.m;
            let mtn_ok = (0.165..=0.372).contains(&m.fractions[3]);
            let low = m.fractions[0] + m.fractions[1];
            let low_ok = (0.375..=0.845).contains(&low);
            let med_ok = (98.0..=390.0).contains(&m.deciles[4]);
            let circ = !c64.wraps && c64.sea_r.is_some_and(|r| r >= 1);
            let degraded = |a: f32, b: f32, cc: f32| (b - cc).abs() > (a - cc).abs() + 0.10 * a.abs();
            let deg: Vec<&str> = [
                ("dimension", degraded(k0.coast.dimension, r.coast.dimension, c_coast.dimension)),
                ("length ratio", degraded(k0.coast.length_ratio, r.coast.length_ratio, c_coast.length_ratio)),
                ("aligned", degraded(k0.coast.aligned, r.coast.aligned, c_coast.aligned)),
                ("block index", degraded(k0.coast.block_index, r.coast.block_index, c_coast.block_index)),
            ]
            .iter()
            .filter(|x| x.1)
            .map(|x| x.0)
            .collect();
            let s = mtn_ok && low_ok && med_ok && circ && deg.is_empty();
            if v == "K1" {
                succ_k1 += s as usize;
                pene += (m.deciles[4] < 50.0) as usize;
            } else {
                succ_k2 += s as usize;
            }
            eprintln!(
                "   seed {seed} SUCCESS test ({v} old): mountain {:.1} % {mtn_ok} · plain + plateau {:.1} % {low_ok} · median {:.0} m {med_ok} · circumnavigable r* ≥ 1 {circ} · coast degraded vs K0 old {} → {}",
                100.0 * m.fractions[3],
                100.0 * low,
                m.deciles[4],
                if deg.is_empty() { "none".into() } else { deg.join(", ") },
                if s { "SUCCESS" } else { "no" }
            );
        }
        let off_b = ["K1 B0", "K1 young 20 Myr"].iter().any(|t| { let r = get(t); r.peak > 7000.0 || !r.slope_ok });
        offscale += off_b as usize;
        all.push(SeedRes { seed, off, runs });
    }
    let n = all.len();
    eprintln!(
        "\n   R: selected {n_sel} / 4 · macro success (old) K1 {succ_k1} / {n} · K2 {succ_k2} / {n} (need {need}) · peneplain (old K1 median < 50 m) {pene} / {n} · off-scale peaks or slope guard failing (K1 B0 / young) {offscale} / {n} → {}",
        if succ_k1.max(succ_k2) >= need {
            "SUCCESS: F172 (the sea-level rise) follows"
        } else if pene >= need {
            "a persistent peneplain: next candidate n > 1 or an incision threshold (Harel et al. 2016), nothing implemented"
        } else if offscale >= need {
            "peaks off the scale: reported, T not recalibrated"
        } else {
            "no success, no stop condition met as written"
        }
    );

    // ── images (north up, the offset at the end) ──
    let shade = |z: &GridF32, nn: usize, px: usize| {
        let (_, hy) = read(z, None);
        whole(&rgba_of(z, Some(&hy), DOMAIN_KM / nn as f32 * 1000.0, &ss), nn, px)
    };
    let cols = ["K0 old", "K1 old 300 Myr", "K2 old 300 Myr", "K1 B0", "K1 young 20 Myr", "K1 intermediate 100 Myr"];
    let (mut r256, mut rcls, mut r2048, mut rrock) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for s in &all {
        let get = |t: &str| s.runs.iter().find(|r| r.tag == t).unwrap();
        r256.push(cols.iter().map(|t| shade(&rolln(&get(t).zm, s.off), 256, 256)).collect::<Vec<_>>());
        rcls.push(cols.iter().map(|t| class_png(&roll_u8(&get(t).classes, 256, [s.off[0] * 4, s.off[1] * 4]), 256, 256, 1)).collect::<Vec<_>>());
        r2048.push(cols[..3].iter().map(|t| shade(&rolln(&get(t).g, s.off), 2048, 420)).collect::<Vec<_>>());
        let pal = [[25u8, 35, 60], [200, 175, 110], [150, 70, 60], [90, 90, 110]];
        rrock.push(
            ["K1 young 20 Myr", "K1 intermediate 100 Myr", "K1 old 300 Myr"]
                .iter()
                .map(|t| {
                    let rk = roll_u8(&get(t).rock, 256, [s.off[0] * 4, s.off[1] * 4]);
                    let mut img = image::RgbaImage::new(256, 256);
                    for y in 0..256u32 {
                        for x in 0..256u32 {
                            let c = pal[rk[(255 - y as usize) * 256 + x as usize] as usize];
                            img.put_pixel(x, y, image::Rgba([c[0], c[1], c[2], 255]));
                        }
                    }
                    img
                })
                .collect::<Vec<_>>(),
        );
    }
    compose(&r256).save(out().join("f171_seeds_256.png")).unwrap();
    compose(&rcls).save(out().join("f171_classes_seeds.png")).unwrap();
    compose(&r2048).save(out().join("f171_old_2048.png")).unwrap();
    compose(&rrock).save(out().join("f171_rock.png")).unwrap();
    // a river profile across a contact: the highest knickpoint of K1 intermediate (first seed with one)
    'prof: for s in &all {
        let r = s.runs.iter().find(|r| r.tag == "K1 intermediate 100 Myr").unwrap();
        let Some(&(kp, h)) = r.knick.2.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()) else { continue };
        let hy = hydrology(&r.zm, CELL);
        let nn = 256;
        let recv = |c: usize| -> Option<usize> {
            let d = hy.dir[c];
            if d == DIR_NONE { return None; }
            let (x, y) = ((c % nn) as i32 + D8_DX[d as usize], (c / nn) as i32 + D8_DY[d as usize]);
            Some(y.rem_euclid(nn as i32) as usize * nn + x.rem_euclid(nn as i32) as usize)
        };
        // walk up 40 cells by the largest donor, then down to the sea
        let mut start = kp;
        for _ in 0..40 {
            let donors: Vec<usize> = (0..nn * nn).filter(|&c| recv(c) == Some(start)).collect();
            match donors.into_iter().max_by(|a, b| hy.acc_km2[*a].partial_cmp(&hy.acc_km2[*b]).unwrap()) {
                Some(d) => start = d,
                None => break,
            }
        }
        let (mut c, mut zs, mut hard) = (start, Vec::new(), Vec::new());
        for _ in 0..400 {
            zs.push(r.zm.data[c]);
            hard.push(r.rock[c] >= 2);
            match recv(c) { Some(nx) if r.zm.data[c] > 0.0 => c = nx, _ => break }
        }
        let soft: Vec<f32> = zs.iter().zip(&hard).map(|(z, h)| if *h { f32::NAN } else { *z }).collect();
        let hardv: Vec<f32> = zs.iter().zip(&hard).map(|(z, h)| if *h { *z } else { f32::NAN }).collect();
        profile_png(&[(zs.clone(), [160, 160, 160]), (soft, [200, 160, 60]), (hardv, [150, 50, 40])], CELL, &out().join("f171_river_contact.png"));
        eprintln!("   river profile: seed {} K1 intermediate, through the knickpoint at cell ({}, {}) of height {h:.0} m; {} cells from 40 upstream to the sea (grey all, ochre cover, red basement)", s.seed, kp % nn, kp / nn, zs.len());
        break 'prof;
    }
    // 25 km crops: Corsica | K1 | K2, at B0 and old (the window of the largest mean slope on K1)
    let crop = |g: &GridF32, nn: usize, x: usize, y: usize| {
        let (_, hy) = read(g, None);
        window(&rgba_of(g, Some(&hy), 195.3, &ss), nn, (x as f32 / nn as f32, y as f32 / nn as f32), 128.0 / nn as f32, 300)
    };
    let pick_cg = {
        let land = land_w(&cg);
        let mut best = (0usize, 0usize, 0f32);
        for oy in (0..1024 - 128).step_by(64) {
            for ox in (0..1024 - 128).step_by(64) {
                let c = (0..128 * 128).filter(|&i| land[(oy + i / 128) * 1024 + ox + i % 128]).count() as f32 / 16384.0;
                if c > 0.95 {
                    let sl: f32 = (0..128 * 127).map(|i| (cg.data[(oy + i / 128) * 1024 + ox + i % 128 + 1024] - cg.data[(oy + i / 128) * 1024 + ox + i % 128]).abs()).sum();
                    if sl > best.2 { best = (ox, oy, sl) }
                }
            }
        }
        best
    };
    let mut crows = Vec::new();
    if let Some(s) = all.first() {
        for age in ["B0", "old 300 Myr"] {
            let g1 = &s.runs.iter().find(|r| r.tag == format!("K1 {age}")).unwrap().g;
            let g2 = &s.runs.iter().find(|r| r.tag == format!("K2 {age}")).unwrap().g;
            let mut best = (0usize, 0usize, f32::MIN);
            for oy in (0..2048 - 128).step_by(128) {
                for ox in (0..2048 - 128).step_by(128) {
                    let sl: f32 = (0..128 * 127).map(|i| (g1.data[(oy + i / 128) * 2048 + ox + i % 128 + 2048] - g1.data[(oy + i / 128) * 2048 + ox + i % 128]).abs()).sum();
                    if sl > best.2 { best = (ox, oy, sl) }
                }
            }
            crows.push(vec![crop(&cg, 1024, pick_cg.0, pick_cg.1), crop(g1, 2048, best.0, best.1), crop(g2, 2048, best.0, best.1)]);
        }
    }
    compose(&crows).save(out().join("f171_crops_25km.png")).unwrap();
    eprintln!(
        "   images: f171_seeds_256.png and f171_classes_seeds.png (rows = the selected seeds; columns {cols:?}), f171_old_2048.png (K0 | K1 | K2 old), f171_rock.png (K1 young | intermediate | old: ochre cover, red basement, grey shield, dark sea), f171_river_contact.png, f171_crops_25km.png (rows B0, old; Corsica | K1 | K2; the first seed) · north up, the display offset at the end · {NOTICE}"
    );
    eprintln!("\n==========  end Finding 171 . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F171-D, the river-profile image redone (the bench's plot was a sparse scatter): the témoin, K1 intermediate, the highest
/// knickpoint's river from 80 cells upstream to the sea, as a line coloured by the rock (ochre cover, red basement).
#[test]
#[ignore]
fn f171_river() {
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let (_, tried) = _f181_sel(PSEED, &MarginOff181::off());
    let s = &tried.last().unwrap().1;
    let p = s.candidates[s.chosen.unwrap()].profile();
    let w = world(PSEED, p);
    let cover = cover_map(PSEED, p, COVER_M);
    let rock = Rock { rho_k: RHO_K, cover_m: cover.clone(), ln_k: None };
    let (r, e_cum) = hist_rock(&w.snaps, &base_cfg(PSEED, rough_a, Some([0, 0])), Some(100.0e6), Some(&rock));
    let land = land_t(&r.zm);
    let rk = rock_at_surface(&land, &cover, &e_cum);
    let (_, _, kps) = knickpoints(&r.zm, &rk);
    let &(kp, h) = kps.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
    let hy = hydrology(&r.zm, CELL);
    let nn = 256;
    let recv = |c: usize| -> Option<usize> {
        let d = hy.dir[c];
        if d == DIR_NONE { return None; }
        let (x, y) = ((c % nn) as i32 + D8_DX[d as usize], (c / nn) as i32 + D8_DY[d as usize]);
        Some(y.rem_euclid(nn as i32) as usize * nn + x.rem_euclid(nn as i32) as usize)
    };
    let mut donors: Vec<Vec<usize>> = vec![Vec::new(); nn * nn];
    for c in 0..nn * nn {
        if let Some(r_) = recv(c) { donors[r_].push(c); }
    }
    let mut start = kp;
    for _ in 0..80 {
        match donors[start].iter().copied().max_by(|a, b| hy.acc_km2[*a].partial_cmp(&hy.acc_km2[*b]).unwrap()) {
            Some(d) => start = d,
            None => break,
        }
    }
    let (mut c, mut pts) = (start, Vec::new());
    for _ in 0..600 {
        pts.push((r.zm.data[c], rk[c] >= 2, c == kp));
        match recv(c) { Some(nx) if r.zm.data[c] > 0.0 => c = nx, _ => break }
    }
    let (w_, h_) = (900u32, 360u32);
    let mut img = image::RgbaImage::from_pixel(w_, h_, image::Rgba([250, 250, 250, 255]));
    let zmax = pts.iter().map(|p| p.0).fold(100f32, f32::max);
    let py = |z: f32| (h_ as f32 - 20.0 - z.max(0.0) / zmax * (h_ as f32 - 40.0)) as i64;
    let px = |i: usize| 20 + (i as f32 / (pts.len() - 1).max(1) as f32 * (w_ as f32 - 40.0)) as i64;
    let mut put = |x: i64, y: i64, cc: [u8; 3]| {
        if x >= 0 && y >= 0 && (x as u32) < w_ && (y as u32) < h_ { img.put_pixel(x as u32, y as u32, image::Rgba([cc[0], cc[1], cc[2], 255])); }
    };
    for x in 20..w_ as i64 - 20 { put(x, py(0.0), [0, 0, 0]); }
    for i in 0..pts.len().saturating_sub(1) {
        let col = if pts[i].1 { [170, 40, 30] } else { [200, 150, 40] };
        for s_ in 0..=20 {
            let t = s_ as f32 / 20.0;
            let (x, y) = (px(i) as f32 * (1.0 - t) + px(i + 1) as f32 * t, py(pts[i].0) as f32 * (1.0 - t) + py(pts[i + 1].0) as f32 * t);
            for d in -1..=1 { put(x as i64, y as i64 + d, col); }
        }
        if pts[i].2 { for d in -8..=8 { put(px(i), py(pts[i].0) + d, [30, 30, 200]); } }
    }
    img.save(out().join("f171_river_contact.png")).unwrap();
    eprintln!("   river profile: the témoin K1 intermediate, the knickpoint at ({}, {}) of height {h:.0} m (blue tick); {} cells, max {zmax:.0} m; ochre cover, red basement", kp % nn, kp / nn, pts.len());
}
