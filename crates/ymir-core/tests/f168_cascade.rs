//! ADR Finding 168 — C1 (b), part 1: a continent that fills the map and can be sailed around. The C1 « continent
//! profile » (the initial continental fraction, a design target tuned once on the témoin), the circumnavigation test,
//! the viz's offset per seed, Davis-Suppe's profile against the distance to a convergent boundary, the frozen cascade.
//! Declared in `docs/reports/relief_method/f168_continent_size/f168_declared.md`.
//!
//! Run: cargo test -p ymir-core --release --test f168_cascade -- --ignored --exact f168_trials --nocapture
//!      cargo test -p ymir-core --release --test f168_cascade -- --ignored --exact f168_cascade --nocapture

#![allow(dead_code)]

mod common;

use common::cascade_bench::*;
use common::{CANONICAL_ORIGIN, Knobs, PSEED, bench_upscale_cfg, pct, sorted};
use std::path::PathBuf;
use std::time::Instant;
use ymir_core::cascade::amplify::restrict;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, ocean_mask};
use ymir_core::cascade::history::{Airy, C1Sources, C1Term, HistoryRun, Mapping, T_F166_YEARS, Timing, physics_history, run_c1_sources_mapped};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{
    CLASS_NAMES, Circumnavigation, Hovius, MOUNTAIN, MacroStats, circumnavigation, classify, coast, dilate_periodic, edt, hovius,
    macro_stats, periodic_components,
};
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, flat_map, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::profile::C1Profile;
use ymir_core::cascade::{CascadeConfig, Rebound, roll};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries, retarget_upper_plate_continental};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::land_topology::land_topology;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const SEEDS: [u64; 4] = [PSEED, 42, 1, 9];
const CELL: f32 = 1.5625;
const N_PHYS: usize = 300;
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];
/// The trials on the témoin (declared): L0 = production's 0.29, then 3, 4, 5, 6 continental plates out of 8.
const TRIALS: [f64; 5] = [0.29, 0.375, 0.5, 0.625, 0.75];

fn out() -> PathBuf {
    root().join("docs/reports/relief_method/f168_continent_size/images")
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
fn roll_mask(m: &[bool], n: usize, off: [usize; 2]) -> Vec<bool> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n]).collect()
}
fn roll_f(m: &[f64], n: usize, off: [usize; 2]) -> Vec<f64> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n]).collect()
}
fn peak_m(z: &GridF32) -> f32 {
    z.data.iter().copied().fold(f32::MIN, f32::max)
}
fn land_of(z: &GridF32) -> Vec<bool> {
    ocean_mask(z).iter().map(|o| !o).collect()
}
fn above0(z: &GridF32) -> Vec<bool> {
    z.data.iter().map(|&v| v > 0.0).collect()
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
            let c = pal[classes[(y / scale) as usize * w + (x / scale) as usize] as usize];
            img.put_pixel(x, y, image::Rgba([c[0], c[1], c[2], 255]));
        }
    }
    img
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

/// The periodic Chebyshev distance (cells) from the `src` cells.
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
/// The periodic Euclidean distance (cells) of each land cell to the sea: the EDT of the 3 × 3 tiling, centre tile.
fn edt_periodic(land: &[bool], n: usize) -> Vec<f32> {
    let t = 3 * n;
    let tiled: Vec<bool> = (0..t * t).map(|k| land[((k / t) % n) * n + (k % t) % n]).collect();
    let d = edt(&tiled, t, t);
    (0..n * n).map(|k| d[(k / n + n) * t + k % n + n]).collect()
}
/// The viz's offset rule (`hd.rs:451-455`): the largest mass's circular-extent centre at the middle of the window.
fn viz_offset(z64_norm: &GridF32) -> [usize; 2] {
    let topo = land_topology(z64_norm, 0.5);
    let (cx, cy) = topo.center_cell;
    [(cx + 64 - 32) % 64, (cy + 64 - 32) % 64]
}

/// One world: a seed, a C1 profile, its sources-driven Airy record, and its framing.
struct World {
    seed: u64,
    name: &'static str,
    profile: C1Profile,
    src: C1Sources,
    st: C1State,
    kin: PlateKinematics,
    /// The continental cells at t₀ (64²), and the land at t₀ (64², the first snapshot above sea).
    cont0: usize,
    land0: Vec<bool>,
    /// The offset applied (the viz's rule on the world's own final 64² field, unless forced).
    off: [usize; 2],
    /// The snapshots, rolled by `off`.
    snaps: Vec<GridF32>,
    /// The unrolled final 64² field (normalised).
    last: GridF32,
    c1_cost: (f64, f64),
}
fn world(seed: u64, name: &'static str, profile: C1Profile, tlf: Option<f32>, force_off: Option<[usize; 2]>) -> World {
    let cfg = c1_cfg();
    let mut st = init_c1_state_phase_2_r7(64, seed, &profile.init_params());
    let cont0 = st.plate_type.data().iter().filter(|t| **t == PlateType::Continental).count();
    let mut kin = PlateKinematics::preset_phase_1_1(st.num_plates);
    let clk = Clock::start();
    let src = run_c1_sources_mapped(
        &mut st,
        &mut kin,
        &cfg,
        &C1Closures::default(),
        10,
        &cfg.iso_config,
        &SteinSteinParams::default(),
        &[C1Term::Erosion],
        Mapping::Airy(Airy::declared()),
    );
    let c1_cost = clk.stop();
    let un = src.hist.normalized(tlf);
    let last = un.last().unwrap().clone();
    let land0 = un[0].data.iter().map(|&v| v > 0.5).collect();
    let off = force_off.unwrap_or_else(|| viz_offset(&last));
    let snaps = un.iter().map(|g| roll(g, off)).collect();
    World { seed, name, profile, src, st, kin, cont0, land0, off, snaps, last, c1_cost }
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
}
fn chain(pm: &GridF32, d5: &[bool], cfg: AmpConfig) -> ChainOut {
    let clk = Clock::start();
    let mut ch = Chain::new(pm.clone(), cfg).with_d5(d5.to_vec());
    for _ in 0..3 {
        ch.next_level(&|| false).unwrap();
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let last = ch.levels.last().unwrap();
    ChainOut { g: last.result(), d5: last.d5.clone(), cost: clk.stop() }
}

struct Run {
    zm: GridF32,
    d5: Vec<bool>,
    classes: Vec<u8>,
    m: MacroStats,
    cost: (f64, f64),
}
fn hist_run(snaps: &[GridF32], cfg: &CascadeConfig) -> Run {
    let clk = Clock::start();
    let run = HistoryRun { t_years: T_F166_YEARS, steps: N_PHYS, rebound: Some(Rebound::declared()), timing: Timing::History };
    let l = physics_history(snaps, 256, cfg, &run, &mut |_| {}, &|| false).unwrap();
    let cost = clk.stop();
    let zm = to_m(&l.rec.z, n2m());
    let oc = ocean_mask(&zm);
    let d5 = (0..65536).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    Run { zm, d5, classes, m, cost }
}

fn temoin_tlf() -> Option<f32> {
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
    bench_upscale_cfg(temoin).target_land_fraction
}
fn base_cfg(seed: u64, rough_a: f32) -> CascadeConfig {
    CascadeConfig { seed, roughness_m: rough_a, mfd_exponent: Some(6.0), ..CascadeConfig::declared(SteinSteinParams::default().depth_scale_m as f32) }
}
fn rough_a() -> f32 {
    0.5 * read_cell(&load_grid(&root().join("data/corsica/corse_128.bin")).unwrap(), None, CELL).0.octaves[0]
}
fn fmt_circ(c: &Circumnavigation, km: f32) -> String {
    match (c.wraps, c.sea_r) {
        (true, _) => format!("WRAPS the torus (main {} cells) → not circumnavigable", c.main_cells),
        (false, Some(r)) => format!("circumnavigable (main {} cells) · r* {r} cells → a sea loop ≥ {} cells = {:.1} km wide", c.main_cells, 2 * r + 1, (2 * r + 1) as f32 * km),
        (false, None) => format!("no land (main {} cells)", c.main_cells),
    }
}

/// F168-B1: the trials on the témoin. The land share = the cells of the 256² history physics level above 0 m, over
/// the whole torus.
#[test]
#[ignore]
fn f168_trials() {
    let t0 = Instant::now();
    let tlf = temoin_tlf();
    let ra = rough_a();
    eprintln!("\n==========  Finding 168 . B1 — the initial continental fraction, trials on the témoin  ==========");
    let mut rows: Vec<(f64, f32)> = Vec::new();
    for f in TRIALS {
        let w = world(PSEED, "trial", if f == 0.29 { C1Profile::Production } else { C1Profile::Continent(f) }, tlf, None);
        let r = hist_run(&w.snaps, &base_cfg(PSEED, ra));
        let land = above0(&r.zm).iter().filter(|&&b| b).count() as f32 / 65536.0;
        let l64_0 = w.land0.iter().filter(|&&b| b).count() as f32 / 4096.0;
        let l64_1 = w.last.data.iter().filter(|&&v| v > 0.5).count() as f32 / 4096.0;
        eprintln!(
            "   f {f:.3}: {} continental plates of 8 (round(8 f)) · continental cells at t₀ {} ({:.1} %) · land 64² t₀ {:.1} % → end {:.1} % · land 256² (physics level) {:.1} % · peak {:.0} m · C1 {} · physics {}",
            ((8.0 * f).round() as usize).clamp(1, 8),
            w.cont0,
            100.0 * w.cont0 as f32 / 4096.0,
            100.0 * l64_0,
            100.0 * l64_1,
            100.0 * land,
            peak_m(&r.zm),
            fmt_cost(w.c1_cost),
            fmt_cost(r.cost)
        );
        rows.push((f, land));
    }
    let mono = rows.windows(2).all(|w| w[1].1 >= w[0].1);
    let near = |t: f32| rows.iter().copied().min_by(|a, b| (a.1 - t).abs().partial_cmp(&(b.1 - t).abs()).unwrap()).unwrap();
    let l1 = near(0.40);
    let l2 = rows.iter().copied().find(|r| (0.55..=0.60).contains(&r.1)).unwrap_or_else(|| near(0.575));
    let reachable = rows.iter().any(|r| (0.50..=0.65).contains(&r.1));
    eprintln!(
        "   monotone in f: {mono} · L1 = f {:.3} ({:.1} %) · L2 = f {:.3} ({:.1} %) · a trial within 50–65 %: {reachable} → {}",
        l1.0,
        100.0 * l1.1,
        l2.0,
        100.0 * l2.1,
        if mono && reachable { "the target is reachable" } else { "STOP: the land target is unreachable" }
    );
    eprintln!("==========  end F168 trials . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F168-C: the measures, 4 seeds × L0 / L1 / L2.
#[test]
#[ignore]
fn f168_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let tlf = temoin_tlf();
    let ra = rough_a();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 168 . C1 (b) part 1: the continent's size, the circumnavigation  ==========");
    let profiles: [(&'static str, C1Profile); 3] = [("L0", C1Profile::Production), ("L1", C1Profile::L1), ("L2", C1Profile::L2)];
    eprintln!("   the profiles: {profiles:?}");
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
    eprintln!("   Corse 195 m: mountain-class slope p50/p90 {:.3}/{:.3} · coast {:?}", c_mtn.r.slope_p50, c_mtn.r.slope_p90, c_coast);
    let corse_oct = |nc: usize| -> f32 {
        let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
        read_cell(&g, None, 200.0 / nc as f32).0.octaves[0]
    };
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct(n / 2))).collect();
    let frozen = || AmpConfig {
        talus_fine_scale: 0.25,
        n1_recal_depth: 2,
        n1_pi_m: pi_table.clone(),
        warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }),
        n1_k: vec![(512, 5.623e2), (1024, 1.0e3), (2048, 3.162e2)],
        ..AmpConfig::declared(PSEED, Variant::N1, 1.0)
    };
    // B3: the negative control — L0 under the canonical roll reproduces F167's témoin
    {
        let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
        let w = world(PSEED, "L0 canonical", C1Profile::Production, tlf, Some(off));
        let r = hist_run(&w.snaps, &base_cfg(PSEED, ra));
        eprintln!(
            "   B3 control: L0 under the canonical roll {off:?}: peak {:.0} m (F167 2 670 m) · land {:.0} km² (F167 31 091 km²) · mountain {:.1} % (F167 74.4 %)",
            peak_m(&r.zm),
            r.m.land_km2,
            100.0 * r.m.fractions[3]
        );
    }
    // the worlds
    struct Done {
        w: World,
        r: Run,
        c: ChainOut,
    }
    let mut all: Vec<Vec<Done>> = Vec::new();
    for seed in SEEDS {
        let mut row = Vec::new();
        for (name, p) in profiles {
            let w = world(seed, name, p, tlf, None);
            let r = hist_run(&w.snaps, &base_cfg(seed, ra));
            let c = chain(&r.zm, &r.d5, frozen());
            eprintln!(
                "   seed {seed} {name}: offset {:?} · C1 {} · physics {} · chain {} → {:.0} s CPU per world",
                w.off,
                fmt_cost(w.c1_cost),
                fmt_cost(r.cost),
                fmt_cost(c.cost),
                w.c1_cost.1 + r.cost.1 + c.cost.1
            );
            row.push(Done { w, r, c });
        }
        all.push(row);
    }
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    let ckm = DOMAIN_KM / 2048.0;
    let (mut successes, mut mtn_over, mut circ_lost) = (0usize, 0usize, 0usize);
    let mut ds_profiles: Vec<Vec<(String, Vec<(f32, usize)>)>> = Vec::new();
    let mut loops: Vec<Vec<Vec<bool>>> = Vec::new();
    let mut reg_rows: Vec<String> = Vec::new();
    for row in &all {
        let seed = row[0].w.seed;
        eprintln!("\n   ══ seed {seed} ══");
        let mut seed_ds = Vec::new();
        let mut seed_loops = Vec::new();
        let l0_coast = coast(&land_of(&row[0].c.g), 2048, 2048, 32);
        let mut mtn_vs_area: Vec<(f32, f32)> = Vec::new();
        for d in row {
            let (w, r, c) = (&d.w, &d.r, &d.c);
            let tag = w.name;
            eprintln!("   ── {tag} ({:?}), offset {:?}:", w.profile, w.off);
            // ── the land ──
            let land256 = above0(&r.zm);
            let nl = land256.iter().filter(|&&b| b).count();
            let (lab, comps) = periodic_components(&land256, 256);
            let big: Vec<usize> = comps.iter().map(|c| c.0).filter(|&s| s as f32 * CELL * CELL > 1000.0).collect();
            let _ = lab;
            eprintln!(
                "      land {:.1} % of the torus ({:.0} km²) · masses > 1 000 km²: {} ({}) · masses in all {}",
                100.0 * nl as f32 / 65536.0,
                nl as f32 * CELL * CELL,
                big.len(),
                big.iter().map(|s| format!("{:.0} km²", *s as f32 * CELL * CELL)).collect::<Vec<_>>().join(", "),
                comps.len()
            );
            // C1's land: t₀ → end (64²)
            let land64: Vec<bool> = w.last.data.iter().map(|&v| v > 0.5).collect();
            let tc = &w.src.hist.type_changed;
            let rift = &w.src.terms.iter().find(|x| x.0 == C1Term::RiftThinning).unwrap().1;
            let lost = (0..4096).filter(|&i| w.land0[i] && !land64[i]).count();
            let lost_rift = (0..4096).filter(|&i| w.land0[i] && !land64[i] && rift[i] < -0.05).count();
            let gained = (0..4096).filter(|&i| !w.land0[i] && land64[i]).count();
            let gained_tc = (0..4096).filter(|&i| !w.land0[i] && land64[i] && tc[i]).count();
            eprintln!(
                "      C1 (64²): continental cells at t₀ {} · land t₀ {} → end {} cells · lost {lost} (rift thinning < −0.05: {lost_rift}) · gained {gained} (type flipped: {gained_tc})",
                w.cont0,
                w.land0.iter().filter(|&&b| b).count(),
                land64.iter().filter(|&&b| b).count()
            );
            // ── the circumnavigation ──
            let c64 = circumnavigation(&land64, 64);
            let land2048 = above0(&c.g);
            let c2048 = circumnavigation(&land2048, 2048);
            let topo = land_topology(&w.last, 0.5);
            eprintln!("      circumnavigation 64² (end of C1): {}", fmt_circ(&c64, 6.25));
            eprintln!("      circumnavigation 2048²: {}", fmt_circ(&c2048, ckm));
            eprintln!("      land_topology: band x {} · band y {} · largest {:.0} km² · bbox {:.0} × {:.0} km", topo.wraps_x, topo.wraps_y, topo.largest_area_km2, topo.bbox_km.0, topo.bbox_km.1);
            // the sea loop (64², rolled): the component of the r*-dilated land holding the main mass
            let blob = match c64.sea_r {
                Some(rs) if !c64.wraps => {
                    let dil = dilate_periodic(&land64, 64, rs);
                    let (l2, _) = periodic_components(&dil, 64);
                    let (l1, cs1) = periodic_components(&land64, 64);
                    let main = (0..cs1.len()).max_by(|a, b| cs1[*a].0.cmp(&cs1[*b].0).then(b.cmp(a))).unwrap() as u32;
                    let rep = (0..4096).find(|&k| l1[k] == main).unwrap();
                    roll_mask(&(0..4096).map(|k| l2[k] == l2[rep]).collect::<Vec<_>>(), 64, w.off)
                }
                _ => vec![false; 4096],
            };
            seed_loops.push(blob);
            // the border (after the offset)
            let border = |n: usize, land: &[bool], z: &GridF32| {
                let cells: Vec<usize> = (0..n * n).filter(|&k| land[k] && (k % n == 0 || k % n == n - 1 || k / n == 0 || k / n == n - 1)).collect();
                let mean_b = cells.iter().map(|&k| z.data[k]).sum::<f32>() / cells.len().max(1) as f32;
                let nl = land.iter().filter(|&&b| b).count().max(1);
                let mean_l = (0..n * n).filter(|&k| land[k]).map(|k| z.data[k]).sum::<f32>() / nl as f32;
                (cells.len(), mean_b, mean_l)
            };
            let (b256, mb256, ml256) = border(256, &land256, &r.zm);
            let (b2048, mb2048, ml2048) = border(2048, &land2048, &c.g);
            eprintln!(
                "      border land after the offset: 256² {b256} cells (mean {mb256:.0} m vs land {ml256:.0} m) · 2048² {b2048} cells (mean {mb2048:.0} m vs land {ml2048:.0} m)"
            );
            // ── the interior ──
            let dist = edt_periodic(&land256, 256);
            let (i20, i40) = ((20.0 / CELL), (40.0 / CELL));
            eprintln!(
                "      interior: land > 20 km from the coast {:.1} % · > 40 km {:.1} % (max {:.0} km)",
                100.0 * dist.iter().filter(|&&v| v > i20).count() as f32 / nl.max(1) as f32,
                100.0 * dist.iter().filter(|&&v| v > i40).count() as f32 / nl.max(1) as f32,
                dist.iter().copied().fold(0f32, f32::max) * CELL
            );
            // ── macro ──
            print_macro(&format!("{tag} classes"), &r.m);
            mtn_vs_area.push((r.m.land_km2 as f32, r.m.fractions[3]));
            // ── C1's tectonics ──
            let st = &w.st;
            let n = 64;
            let plates: std::collections::BTreeSet<u16> = st.plate_id.data().iter().copied().collect();
            let mut bi = classify_boundaries(&st.plate_id, &w.kin);
            retarget_upper_plate_continental(&mut bi, &st.plate_id, &st.plate_type);
            let cont = |i: usize, j: usize| st.plate_type.get(i, j) == PlateType::Continental;
            let mut coll = 0usize;
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
            let ds = &w.src.terms.iter().find(|x| x.0 == C1Term::DavisSuppe).unwrap().1;
            let arc = &w.src.terms.iter().find(|x| x.0 == C1Term::Subduction).unwrap().1;
            let nl64 = land64.iter().filter(|&&b| b).count().max(1);
            let ds_share = (0..4096).filter(|&i| land64[i] && ds[i] > 0.05).count() as f32 / nl64 as f32;
            let crat = up_mask(&roll_mask(st.cratonic_mask.data(), 64, w.off), 64, 4);
            let crat_land: Vec<bool> = (0..65536).map(|i| crat[i] && land256[i]).collect();
            let cr = sorted((0..65536).filter(|&i| crat_land[i]).map(|i| r.zm.data[i]).collect());
            let cblock = coast(&crat_land, 256, 256, 4).block_index;
            eprintln!(
                "      C1 tectonics: plates at the end {} · C-C collision cells {coll} · upper-plate (active margin) cells {} · Davis-Suppe Σ Δs > 0.05 on {:.1} % of the land · cratons {} cells (64²), land p50 {:.0} m (256²), block index {:.2}",
                plates.len(),
                bi.upper_plate_count(),
                100.0 * ds_share,
                st.cratonic_mask.data().iter().filter(|&&c| c).count(),
                pct(&cr, 0.5),
                cblock
            );
            // the Davis-Suppe profile against the distance to a convergent cell or suture (64², unrolled)
            let near: Vec<bool> = (0..4096).map(|k| conv[k] || w.src.hist.suture[k]).collect();
            let dn = cheb_periodic(&near, 64);
            let mut bins = vec![(0f32, 0usize); 12];
            for i in 0..4096 {
                if land64[i] {
                    let b = (dn[i] as usize).min(11);
                    bins[b].0 += ds[i] as f32;
                    bins[b].1 += 1;
                }
            }
            let prof: Vec<(f32, usize)> = bins.iter().map(|&(s, c)| (if c > 0 { s / c as f32 } else { f32::NAN }, c)).collect();
            eprintln!(
                "      Σ Δs_DS by distance to a convergent cell / suture (cells of 6.25 km; mean [land cells]): {}",
                prof.iter().enumerate().map(|(i, (m, c))| format!("{}{i}: {m:.3} [{c}]", if i == 11 { ">" } else { "" })).collect::<Vec<_>>().join(" · ")
            );
            seed_ds.push((tag.to_string(), prof));
            // ── 2048²: peaks, slope guard, texture, F, coast ──
            let m8 = up_u8(&r.classes, 256, 8);
            let mtn8: Vec<bool> = m8.iter().map(|&x| x == MOUNTAIN).collect();
            let excl: Vec<bool> = (0..2048 * 2048).map(|i| c.d5[i] || !mtn8[i]).collect();
            let whole_t = tex(&c.g, ckm, &c.d5);
            let mtn_t = tex(&c.g, ckm, &excl);
            for (t, x, cc) in [("whole", &whole_t, &c_whole), ("mountain", &mtn_t, &c_mtn)] {
                let el = elements(x, cc, ckm, ctrl_ok);
                let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
                eprintln!("      2048² {t:<9} {} of {} within ×1.5 · out: {}", el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") });
            }
            eprintln!(
                "      slope guard: peak 256² {:.0} m · 2048² {:.0} m · mountain-class slope p50/p90 {:.3}/{:.3} (Corsica {:.3}/{:.3})",
                peak_m(&r.zm),
                peak_m(&c.g),
                mtn_t.r.slope_p50,
                mtn_t.r.slope_p90,
                c_mtn.r.slope_p50,
                c_mtn.r.slope_p90
            );
            let hy = hydrology(&r.zm, CELL);
            let lake_cells = hy.lake_id.iter().filter(|&&l| l != 0).count();
            eprintln!("      lakes (256²): {} / {:.0} km²", hy.lakes.len(), lake_cells as f32 * CELL * CELL);
            let (rd, fl) = flat_map(&c.g, ckm, Some(&c.d5));
            let per: Vec<String> = (1..5u8)
                .map(|cl| {
                    let cells = (0..2048 * 2048).filter(|&i| rd[i] && m8[i] == cl).count();
                    let f = (0..2048 * 2048).filter(|&i| rd[i] && fl[i] && m8[i] == cl).count();
                    format!("{} {:.2} %", CLASS_NAMES[cl as usize], 100.0 * f as f32 / cells.max(1) as f32)
                })
                .collect();
            eprintln!("      F by class: {}", per.join(" · "));
            let kc = coast(&land_of(&c.g), 2048, 2048, 32);
            eprintln!("      coast {kc:?}");
            if tag == "L2" {
                let p2048 = &whole_t.p;
                let p1024 = {
                    let g = restrict(&c.g);
                    let dd = (0..1024 * 1024).map(|i| c.d5[(i / 1024 * 2) * 2048 + (i % 1024) * 2]).collect::<Vec<_>>();
                    predictability(&g, 2.0 * ckm, Some(&dd))
                };
                let hv: Vec<String> = [10.0f32, 25.0, 100.0]
                    .iter()
                    .map(|&a| {
                        let h: Hovius = hovius(&c.g, &mtn8, ckm, a);
                        format!("A ≥ {a}: R {:.2} ({} basins)", h.ratio, h.outlets)
                    })
                    .collect();
                reg_rows.push(format!(
                    "   seed {seed} L2: 391 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | 195 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | Hovius {}",
                    p1024.cv_lambda, p1024.cv_len, p1024.theta_mean, p1024.theta_std, p1024.r2_gully, p2048.cv_lambda, p2048.cv_len, p2048.theta_mean, p2048.theta_std, p2048.r2_gully, hv.join(" · ")
                ));
                // R: the success test, against L0 (F167's configuration under the viz's offset)
                let degraded = |a: f32, b: f32, cc: f32| (b - cc).abs() > (a - cc).abs() + 0.10 * a.abs();
                let deg: Vec<&str> = [
                    ("dimension", degraded(l0_coast.dimension, kc.dimension, c_coast.dimension)),
                    ("length ratio", degraded(l0_coast.length_ratio, kc.length_ratio, c_coast.length_ratio)),
                    ("aligned", degraded(l0_coast.aligned, kc.aligned, c_coast.aligned)),
                    ("block index", degraded(l0_coast.block_index, kc.block_index, c_coast.block_index)),
                ]
                .iter()
                .filter(|x| x.1)
                .map(|x| x.0)
                .collect();
                let m = &r.m;
                let mtn_ok = (0.165..=0.372).contains(&m.fractions[3]);
                let low = m.fractions[0] + m.fractions[1];
                let low_ok = (0.375..=0.845).contains(&low);
                let circ = !c64.wraps && !c2048.wraps && c64.main_cells > 0;
                let succ = mtn_ok && low_ok && circ && deg.is_empty();
                successes += succ as usize;
                mtn_over += (m.fractions[3] > 0.372) as usize;
                circ_lost += (!circ) as usize;
                eprintln!(
                    "      SUCCESS test (L2): mountain {:.1} % {mtn_ok} · plain + plateau {:.1} % {low_ok} · circumnavigable (64² and 2048²) {circ} · coast degraded vs L0 {} → {}",
                    100.0 * m.fractions[3],
                    100.0 * low,
                    if deg.is_empty() { "none".into() } else { deg.join(", ") },
                    if succ { "SUCCESS" } else { "no" }
                );
                // the mountain class decomposed (256²): cratons, Davis-Suppe zones by distance to the margin, arcs, the rest
                let dsz = |lo: u32, hi: u32| -> Vec<bool> {
                    up_mask(&roll_mask(&(0..4096).map(|k| ds[k] > 0.05 && dn[k] >= lo && dn[k] <= hi).collect::<Vec<_>>(), 64, w.off), 64, 4)
                };
                let (ds_a, ds_b, ds_c) = (dsz(0, 5), dsz(6, 15), dsz(16, u32::MAX));
                let arcm = up_mask(&roll_mask(&arc.iter().map(|&v| v > 0.05).collect::<Vec<_>>(), 64, w.off), 64, 4);
                let mc: Vec<usize> = (0..65536).filter(|&i| r.classes[i] == MOUNTAIN).collect();
                let nm = mc.len().max(1) as f32;
                let sh = |f: &dyn Fn(usize) -> bool| mc.iter().filter(|&&i| f(i)).count() as f32 / nm;
                let c_cr = sh(&|i| crat[i]);
                let c_a = sh(&|i| !crat[i] && ds_a[i]);
                let c_b = sh(&|i| !crat[i] && ds_b[i]);
                let c_c = sh(&|i| !crat[i] && ds_c[i]);
                let c_arc = sh(&|i| !crat[i] && !ds_a[i] && !ds_b[i] && !ds_c[i] && arcm[i]);
                eprintln!(
                    "      the mountain class (256²) decomposed: cratons {:.1} % · Davis-Suppe zones ≤ 5 cells (31 km) from a convergent cell / suture {:.1} % · 6–15 cells {:.1} % · > 15 cells {:.1} % · arcs {:.1} % · the rest {:.1} %",
                    100.0 * c_cr,
                    100.0 * c_a,
                    100.0 * c_b,
                    100.0 * c_c,
                    100.0 * c_arc,
                    100.0 * (1.0 - c_cr - c_a - c_b - c_c - c_arc)
                );
            }
        }
        if mtn_vs_area.len() == 3 {
            let (a0, m0) = mtn_vs_area[0];
            let (a2, m2) = mtn_vs_area[2];
            eprintln!(
                "   seed {seed}: the mountain share against the land area: {} · slope L0 → L2 {:+.2} points per 10 000 km²",
                mtn_vs_area.iter().map(|(a, m)| format!("{a:.0} km² → {:.1} %", 100.0 * m)).collect::<Vec<_>>().join(" · "),
                100.0 * (m2 - m0) / ((a2 - a0) / 10000.0)
            );
        }
        ds_profiles.push(seed_ds);
        loops.push(seed_loops);
    }
    eprintln!(
        "\n   the regularity baseline (Corsica 391 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3} | 195 m: CV_λ {:.3} · CV_L {:.3} · θ {:.1}° σ {:.1}° · R2_g {:.3}):",
        cp391.cv_lambda, cp391.cv_len, cp391.theta_mean, cp391.theta_std, cp391.r2_gully, c_whole.p.cv_lambda, c_whole.p.cv_len, c_whole.p.theta_mean, c_whole.p.theta_std, c_whole.p.r2_gully
    );
    for r in &reg_rows {
        eprintln!("{r}");
    }
    eprintln!(
        "\n   R: macro success in L2 on {successes} / 4 · mountains > 37.2 % on {mtn_over} / 4 · circumnavigation lost on {circ_lost} / 4 → {}",
        if successes >= 3 {
            "SUCCESS: F169 (the sea-level rise) on L2"
        } else if mtn_over >= 2 {
            "the mountains still exceed 37.2 %: decomposed above, stop"
        } else if circ_lost >= 1 {
            "the circumnavigation is lost: report, propose, implement nothing"
        } else {
            "no success, no stop condition met"
        }
    );
    // ── images ──
    let shade = |z: &GridF32, n: usize, px: usize, blob: Option<&[bool]>| {
        let (_, hy) = read(z, None);
        let mut rgba = rgba_of(z, Some(&hy), DOMAIN_KM / n as f32 * 1000.0, &ss);
        if let Some(b) = blob {
            // the sea loop: the boundary of the r*-dilated main blob (64², drawn at n)
            let f = n / 64;
            for k in 0..n * n {
                let (x, y) = (k % n, k / n);
                let me_ = b[(y / f) * 64 + x / f];
                let edge = [((x + 1) % n, y), (x, (y + 1) % n)].iter().any(|&(a, c)| b[(c / f) * 64 + a / f] != me_);
                if edge {
                    rgba[4 * k..4 * k + 4].copy_from_slice(&[230, 40, 40, 255]);
                }
            }
        }
        whole(&rgba, n, px)
    };
    let (mut r256, mut r2048, mut rcls) = (Vec::new(), Vec::new(), Vec::new());
    for (row, lp) in all.iter().zip(&loops) {
        r256.push(row.iter().zip(lp).map(|(d, b)| shade(&d.r.zm, 256, 320, Some(b))).collect::<Vec<_>>());
        r2048.push(row.iter().map(|d| shade(&d.c.g, 2048, 448, None)).collect::<Vec<_>>());
        rcls.push(row.iter().map(|d| class_png(&d.r.classes, 256, 256, 1)).collect::<Vec<_>>());
    }
    compose(&r256).save(out().join("f168_seeds_256.png")).unwrap();
    compose(&r2048).save(out().join("f168_seeds_2048.png")).unwrap();
    compose(&rcls).save(out().join("f168_classes_seeds.png")).unwrap();
    // the Davis-Suppe profile: one panel per seed, L0 grey, L1 blue, L2 red
    let ymax = ds_profiles.iter().flatten().flat_map(|(_, p)| p.iter().map(|x| x.0)).filter(|v| v.is_finite()).fold(0.05f32, f32::max);
    let mut panels = Vec::new();
    for sp in &ds_profiles {
        let (w, h) = (420u32, 300u32);
        let mut img = image::RgbaImage::from_pixel(w, h, image::Rgba([250, 250, 250, 255]));
        let px = |b: f32| (30.0 + b / 11.0 * (w as f32 - 50.0)) as i64;
        let py = |v: f32| (h as f32 - 25.0 - v / ymax * (h as f32 - 45.0)) as i64;
        let mut put = |x: i64, y: i64, c: [u8; 3]| {
            if x >= 0 && y >= 0 && (x as u32) < w && (y as u32) < h {
                img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], 255]));
            }
        };
        for x in 30..w as i64 - 20 {
            put(x, py(0.0), [0, 0, 0]);
            put(x, py(0.05), [200, 200, 200]);
        }
        for (i, (_, p)) in sp.iter().enumerate() {
            let col = [[120u8, 120, 120], [40, 90, 200], [210, 50, 40]][i];
            for b in 0..11 {
                let (a, c) = (p[b].0, p[b + 1].0);
                if !(a.is_finite() && c.is_finite()) {
                    continue;
                }
                for s in 0..=40 {
                    let t = s as f32 / 40.0;
                    for dd in -1..=1 {
                        put(px(b as f32 + t), py(a + (c - a) * t) + dd, col);
                    }
                }
            }
        }
        panels.push(img);
    }
    compose(&[panels]).save(out().join("f168_ds_profile.png")).unwrap();
    // 25 km coast crops: Corsica | L0 | L1 | L2 per seed (the window picked on L2)
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
    let crop = |g: &GridF32, n: usize, x: usize, y: usize| {
        let (_, hy) = read(g, None);
        window(&rgba_of(g, Some(&hy), 195.3, &ss), n, (x as f32 / n as f32, y as f32 / n as f32), 128.0 / n as f32, 320)
    };
    let (cx, cy) = pick(&cg, 1024);
    let mut crows = Vec::new();
    for row in &all {
        // each world's own first half-land window (the coast differs per profile)
        let mut v = vec![crop(&cg, 1024, cx, cy)];
        v.extend(row.iter().map(|d| {
            let (ox, oy) = pick(&d.c.g, 2048);
            crop(&d.c.g, 2048, ox, oy)
        }));
        crows.push(v);
    }
    compose(&crows).save(out().join("f168_coast_25km.png")).unwrap();
    eprintln!(
        "   images: f168_seeds_256.png (the sea loop red), f168_seeds_2048.png, f168_classes_seeds.png (rows = seeds {SEEDS:?}; columns L0 | L1 | L2), f168_ds_profile.png (panels = seeds; x = distance 0–11 cells, y = mean Σ Δs_DS 0–{ymax:.2}; the 0.05 line grey; L0 grey, L1 blue, L2 red), f168_coast_25km.png (Corsica | L0 | L1 | L2 per seed, each world's first 40–60 % land window) · north up · {NOTICE}"
    );
    eprintln!("\n==========  end Finding 168 . {:.1} s wall  ==========\n", t0.elapsed().as_secs_f64());
}

/// F168, added AFTER the measure (declared in the report as an addition, to inform the proposal only): is the
/// circumnavigation already lost at t₀ (the first snapshot's land, and the continental-plate mask), or made by C1?
#[test]
#[ignore]
fn f168_t0() {
    let tlf = temoin_tlf();
    for seed in SEEDS {
        for (name, p) in [("L0", C1Profile::Production), ("L1", C1Profile::L1), ("L2", C1Profile::L2)] {
            let st0 = init_c1_state_phase_2_r7(64, seed, &p.init_params());
            let cont: Vec<bool> = st0.plate_type.data().iter().map(|t| *t == PlateType::Continental).collect();
            let w = world(seed, name, p, tlf, None);
            let end: Vec<bool> = w.last.data.iter().map(|&v| v > 0.5).collect();
            let (a, b, c) = (circumnavigation(&cont, 64), circumnavigation(&w.land0, 64), circumnavigation(&end, 64));
            eprintln!(
                "   seed {seed} {name}: continental plates t₀: {} | land t₀: {} | land end: {}",
                fmt_circ(&a, 6.25),
                fmt_circ(&b, 6.25),
                fmt_circ(&c, 6.25)
            );
        }
    }
}
