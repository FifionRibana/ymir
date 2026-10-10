//! ADR Finding 165 — the main bench: A2 (the diagnostic on four seeds), the T calibration, the historical physics
//! level and its controls, the frozen cascade to 2 048², the macro / texture / coast measures, and the images.

use super::common::cascade_bench::*;
use super::common::{CANONICAL_ORIGIN, Knobs, PSEED, bench_upscale_cfg, pct, sorted};
use super::{CELL, class_png, europe_grids, print_macro};
use std::collections::BTreeMap;
use std::time::Instant;
use ymir_core::cascade::amplify::{AmpConfig, Chain, Variant, Warp, WarpPlace, ocean_mask};
use ymir_core::cascade::history::{C1History, HistoryLevel, HistoryRun, Timing, physics_history, run_c1_recorded};
use ymir_core::cascade::hydro::hydrology;
use ymir_core::cascade::measure::octave_km;
use ymir_core::cascade::physio::{MOUNTAIN, MacroStats, Shape, classify, coast, macro_stats, shape};
use ymir_core::cascade::planform::{Planform, WINDOWS_KM, planform};
use ymir_core::cascade::predict::{Predict, predictability};
use ymir_core::cascade::{CascadeConfig, Rebound, diff_rgba, physics_level, roll, upsample_to};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::boundary_classification::{BoundaryType, classify_boundaries, retarget_upper_plate_continental};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_c1::valley_construction::{F121_AGE_K, ValleyConstruction};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const SEEDS: [u64; 4] = [PSEED, 42, 1, 9];
const N_PHYS: usize = 300;
const VALID_P: [bool; 7] = [true, true, false, false, true, false, true];

fn out() -> std::path::PathBuf {
    root().join("docs/reports/relief_method/f165_continent/images")
}
fn out15(a: f32, b: f32) -> bool {
    !(a.is_finite() && b.is_finite() && b != 0.0) || a > 1.5 * b || a < b / 1.5
}
fn n2m() -> f32 {
    2.0 * 1.13 * SteinSteinParams::default().depth_scale_m as f32
}
fn roll_mask(m: &[bool], n: usize, off: [usize; 2]) -> Vec<bool> {
    (0..n * n).map(|k| m[((k / n + off[1]) % n) * n + (k % n + off[0]) % n]).collect()
}
fn up_mask(m: &[bool], n: usize, f: usize) -> Vec<bool> {
    let nn = n * f;
    (0..nn * nn).map(|k| m[(k / nn / f) * n + (k % nn) / f]).collect()
}
fn up_u8(m: &[u8], n: usize, f: usize) -> Vec<u8> {
    let nn = n * f;
    (0..nn * nn).map(|k| m[(k / nn / f) * n + (k % nn) / f]).collect()
}
fn peak_m(z: &GridF32) -> f32 {
    let oc = ocean_mask(z);
    z.data.iter().zip(&oc).filter(|(_, o)| !**o).map(|(v, _)| *v).fold(f32::MIN, f32::max)
}

struct Seed {
    seed: u64,
    state: C1State,
    kin: PlateKinematics,
    hist: C1History,
    /// The normalised snapshots, rolled to the framing.
    snaps: Vec<GridF32>,
}

fn run_seed(seed: u64, tlf: Option<f32>) -> Seed {
    let ss = SteinSteinParams::default();
    let cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let hist = run_c1_recorded(&mut state, &mut kin, &cfg, &C1Closures::default(), 10, &cfg.iso_config, &ss);
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let snaps: Vec<GridF32> = hist.normalized(tlf).iter().map(|g| roll(g, off)).collect();
    // the last snapshot is the bench's coarse field
    let coarse = roll(&c1_coarse_normalized_altitude(&state, &cfg.iso_config, &ss, tlf), off);
    let err = coarse.data.iter().zip(&snaps.last().unwrap().data).map(|(a, b)| (a - b).abs()).fold(0f32, f32::max);
    assert!(err < 1e-5, "the final snapshot is the coarse field (max diff {err})");
    Seed { seed, state, kin, hist, snaps }
}

fn fmt_shape(s: &Shape) -> String {
    format!(
        "{} cells · {} components · width {:.2} cells ({:.1} km) · elongation {:.2} · near convergent/suture {:.0} %",
        s.cells,
        s.components,
        s.width_cells,
        s.width_cells * 6.25,
        s.elongation,
        100.0 * s.near_share
    )
}

struct Phys {
    name: String,
    zm: GridF32,
    d5: Vec<bool>,
    classes: Vec<u8>,
    m: MacroStats,
    secs: f64,
    rec_sums: Option<(Vec<f32>, Vec<f32>, Vec<f32>)>,
}

fn phys_of(name: &str, zm: GridF32, d5: Vec<bool>, secs: f64, sums: Option<(Vec<f32>, Vec<f32>, Vec<f32>)>) -> Phys {
    let classes = classify(&zm.data, 256, 256, CELL);
    let m = macro_stats(&zm.data, &classes, 256, 256, CELL);
    Phys { name: name.to_string(), zm, d5, classes, m, secs, rec_sums: sums }
}

fn history_phys(name: &str, s: &Seed, cfg: &CascadeConfig, run: &HistoryRun) -> (Phys, HistoryLevel) {
    let t = Instant::now();
    let l = physics_history(&s.snaps, 256, cfg, run, &mut |_| {}, &|| false).unwrap();
    let secs = t.elapsed().as_secs_f64();
    let k = n2m();
    let zm = to_m(&l.rec.z, k);
    let oc = ocean_mask(&zm);
    let d5: Vec<bool> = (0..256 * 256).map(|i| !oc[i] && zm.data[i] <= 0.0).collect();
    let sums = (
        l.rec.sum_uplift.iter().map(|v| v * k).collect(),
        l.rec.sum_erosion.iter().map(|v| v * k).collect(),
        l.rec.sum_rebound.iter().map(|v| v * k).collect(),
    );
    (phys_of(name, zm, d5, secs, Some(sums)), l)
}

/// The texture battery against Corsica (F163–F164's instruments); `None` excl = the whole world.
struct Tex {
    r: Reading,
    p: Predict,
    pf: Planform,
}
fn tex(g: &GridF32, cell_km: f32, excl: &[bool]) -> Tex {
    let (r, _) = read_cell(g, Some(excl), cell_km);
    let p = predictability(g, cell_km, Some(excl));
    let pf = planform(g, cell_km, Some(excl));
    Tex { r, p, pf }
}
/// The battery's elements and whether each is within ×1.5 of the reference.
fn elements(a: &Tex, c: &Tex, cell_km: f32, ctrl_ok: bool) -> Vec<(String, f32, f32, bool)> {
    let mut v = Vec::new();
    v.push(("facets".to_string(), a.r.facets as f32, c.r.facets as f32, !(ctrl_ok && a.r.facets > 1.5 * c.r.facets && a.r.facet_cells > 20)));
    v.push(("walls 28°".to_string(), a.r.walls28 as f32, c.r.walls28 as f32, !(a.r.walls28 > 1.5 * c.r.walls28 && a.r.walls28_cells > 20)));
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
    v.push(("A_dir".to_string(), a.pf.a_dir, c.pf.a_dir, !out15(a.pf.a_dir, c.pf.a_dir)));
    v.push(("F".to_string(), a.pf.flat, c.pf.flat, !out15(a.pf.flat, c.pf.flat)));
    v
}

/// New depressions after a warp: components of the after-field's lakes none of whose cells was a lake before.
fn new_depressions(a: &GridF32, b: &GridF32, cell_km: f32) -> (usize, usize, usize) {
    let (ha, hb) = (hydrology(a, cell_km), hydrology(b, cell_km));
    let mut fresh = vec![true; hb.lakes.len()];
    for k in 0..hb.lake_id.len() {
        if hb.lake_id[k] != 0 && ha.lake_id[k] != 0 {
            fresh[hb.lake_id[k] as usize - 1] = false;
        }
    }
    (ha.lakes.len(), hb.lakes.len(), fresh.iter().filter(|&&f| f).count())
}

struct ChainOut {
    g: GridF32,
    d5: Vec<bool>,
    secs: f64,
    notes: Vec<String>,
}

fn chain(pm: &GridF32, d5: &[bool], cfg: AmpConfig) -> ChainOut {
    let t = Instant::now();
    let mut ch = Chain::new(pm.clone(), cfg).with_d5(d5.to_vec());
    let mut prev = pm.clone();
    let mut notes = Vec::new();
    for _ in 0..3 {
        let l = ch.next_level(&|| false).unwrap().clone();
        let res = l.result();
        let n = l.n;
        let cell = DOMAIN_KM / n as f32;
        let (b, _) = drift(&res, &prev);
        let mean = {
            let oc = ocean_mask(&prev);
            let v: Vec<f32> = prev.data.iter().zip(&oc).filter(|(_, o)| !**o).map(|(v, _)| *v).collect();
            v.iter().sum::<f32>() / v.len().max(1) as f32
        };
        let dep = l.pre_warp.as_ref().zip(l.post_warp.as_ref()).map(|(a, b)| new_depressions(a, b, cell));
        notes.push(format!(
            "{n}²: peak {:.0} m · drift {:+.2} %{}",
            peak_m(&res),
            100.0 * b / mean.max(1.0),
            dep.map_or(String::new(), |(a, b, f)| format!(" · warp depressions {a} → {b}, NEW {f} ({:.1} % of before)", 100.0 * f as f32 / a.max(1) as f32))
        ));
        prev = res;
        let keep = ch.levels.len() - 1;
        ch.levels.drain(..keep);
    }
    let last = ch.levels.last().unwrap();
    ChainOut { g: last.result(), d5: last.d5.clone(), secs: t.elapsed().as_secs_f64(), notes }
}

#[test]
#[ignore]
fn f165_cascade() {
    let t0 = Instant::now();
    let ss = SteinSteinParams::default();
    let k = n2m();
    std::fs::create_dir_all(out()).unwrap();
    eprintln!("\n==========  Finding 165 . the continent's physiography: the physics level follows C1's history  ==========");
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
    let tlf = bench_upscale_cfg(temoin).target_land_fraction;
    eprintln!("   the benches' target_land_fraction: {tlf:?}");
    // ── Europe (B, frozen) ──
    let (e, _) = europe_grids();
    let ce = classify(&e.data, e.width, e.height, CELL);
    let me = macro_stats(&e.data, &ce, e.width, e.height, CELL);
    print_macro("Europe", &me);
    // ── Corsica (texture, coast) ──
    let cg = load_grid(&root().join("data/corsica/corse_1024.bin")).unwrap();
    let c128 = load_grid(&root().join("data/corsica/corse_128.bin")).unwrap();
    let cclass = up_u8(&classify(&c128.data, 128, 128, CELL), 128, 8);
    let ck = 200.0 / 1024.0;
    let c_whole = tex(&cg, ck, &vec![false; 1024 * 1024]);
    let c_mtn = tex(&cg, ck, &cclass.iter().map(|&c| c != MOUNTAIN).collect::<Vec<_>>());
    let c_coast = coast(&ocean_mask(&cg).iter().map(|o| !o).collect::<Vec<_>>(), 1024, 1024, 32);
    let c49 = load_grid(&root().join("data/corsica/corse_4096.bin")).unwrap();
    let c_coast49 = coast(&ocean_mask(&c49).iter().map(|o| !o).collect::<Vec<_>>(), 4096, 4096, 128);
    eprintln!("   Corse 195 m coast: {c_coast:?} · 49 m (lattice 128): {c_coast49:?}");
    let corse_oct: BTreeMap<usize, f32> = [256usize, 512, 1024]
        .iter()
        .map(|&nc| {
            let g = load_grid(&root().join(format!("data/corsica/corse_{nc}.bin"))).unwrap();
            (nc, read_cell(&g, None, 200.0 / nc as f32).0.octaves[0])
        })
        .collect();
    let rough_a = 0.5 * read_cell(&c128, None, CELL).0.octaves[0];
    let pi_table: Vec<(usize, f32)> = [512usize, 1024, 2048].iter().map(|&n| (n, 0.15 * corse_oct[&(n / 2)])).collect();
    let frozen = || AmpConfig {
        talus_fine_scale: 0.25,
        n1_recal_depth: 2,
        n1_pi_m: pi_table.clone(),
        warp: Some(Warp { amp: 0.5, corr: 4.0, place: WarpPlace::All }),
        n1_k: vec![(512, 5.623e2), (1024, 1.0e3), (2048, 3.162e2)],
        ..AmpConfig::declared(PSEED, Variant::N1, 1.0)
    };
    // ── the four seeds' C1 runs ──
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    let mut seeds = Vec::new();
    for s in SEEDS {
        let t = Instant::now();
        seeds.push(run_seed(s, tlf));
        eprintln!("   C1 seed {s} recorded ({:.1} s)", t.elapsed().as_secs_f64());
    }
    // ── A2 ──
    eprintln!("\n   A2 — the diagnostic (64², 6.25 km per cell; final state):");
    let mut a2_widths: Vec<(f32, f32, f32)> = Vec::new();
    for s in &seeds {
        let st = &s.state;
        let n = 64;
        let plates: std::collections::BTreeSet<u16> = st.plate_id.data().iter().copied().collect();
        let mut bi = classify_boundaries(&st.plate_id, &s.kin);
        retarget_upper_plate_continental(&mut bi, &st.plate_id, &st.plate_type);
        let cont = |i: usize, j: usize| st.plate_type.get(i, j) == PlateType::Continental;
        let (mut coll, mut rift) = (0usize, 0usize);
        let mut conv = vec![false; n * n];
        for j in 0..n {
            for i in 0..n {
                let b = bi.boundary_type.get(i, j);
                if b == BoundaryType::Convergent {
                    conv[j * n + i] = true;
                    let me = st.plate_id.get(i, j);
                    let nb = [((i + 1) % n, j), ((i + n - 1) % n, j), (i, (j + 1) % n), (i, (j + n - 1) % n)];
                    if cont(i, j) && nb.iter().any(|&(x, y)| st.plate_id.get(x, y) != me && cont(x, y)) {
                        coll += 1;
                    }
                }
                if b == BoundaryType::Divergent && cont(i, j) {
                    rift += 1;
                }
            }
        }
        let upper = bi.upper_plate_count();
        let craton = st.cratonic_mask.data().iter().filter(|&&c| c).count();
        let near_raw: Vec<bool> = (0..n * n).map(|k| conv[k] || s.hist.suture[k]).collect();
        let near = roll_mask(&near_raw, n, off);
        let changed = roll_mask(&s.hist.type_changed, n, off);
        let hm: Vec<f32> = s.snaps.last().unwrap().data.iter().map(|v| (v - 0.5) * k).collect();
        let land: Vec<bool> = hm.iter().map(|&v| v > 0.0).collect();
        let nl = land.iter().filter(|&&l| l).count();
        let alts = sorted((0..n * n).filter(|&i| land[i]).map(|i| hm[i]).collect());
        eprintln!(
            "   ── seed {}: plates {} (init {}) · collision {} · subduction upper plate {} · craton {} · rift {} · sutures (convergent at any step) {} · land {} cells, of which type changed {} ({:.1} %)",
            s.seed,
            plates.len(),
            st.num_plates,
            coll,
            upper,
            craton,
            rift,
            s.hist.suture.iter().filter(|&&b| b).count(),
            nl,
            (0..n * n).filter(|&i| land[i] && changed[i]).count(),
            100.0 * (0..n * n).filter(|&i| land[i] && changed[i]).count() as f32 / nl.max(1) as f32
        );
        eprintln!("      h_iso land deciles (m): {}", (1..10).map(|d| format!("{:.0}", pct(&alts, d as f64 / 10.0))).collect::<Vec<_>>().join("/"));
        let mut w1000 = f32::NAN;
        for thr in [500.0f32, 1000.0, 2000.0] {
            let m: Vec<bool> = hm.iter().map(|&v| v > thr).collect();
            let sh = shape(&m, n, n, Some(&near), 2);
            if thr == 1000.0 {
                w1000 = sh.width_cells;
                a2_widths.push((sh.width_cells, sh.near_share, f32::NAN));
            }
            eprintln!("      h_iso > {thr:.0} m: {}", fmt_shape(&sh));
        }
        let _ = w1000;
    }
    // ── the physics levels ──
    let base = |seed: u64| CascadeConfig { seed, roughness_m: rough_a, mfd_exponent: Some(6.0), ..CascadeConfig::declared(ss.depth_scale_m as f32) };
    let rb = Rebound::declared();
    eprintln!("\n   C — the physics levels (256²); rebound α {:.1} km (Te 25 km), ρc/ρm {:.3}", rb.alpha_km, rb.rho_c / rb.rho_m);
    // T calibrated once, on the témoin
    let s0 = &seeds[0];
    let mut trials: Vec<(f64, f32)> = Vec::new();
    let peak_at = |t: f64, trials: &mut Vec<(f64, f32)>| -> f32 {
        let (p, _) = history_phys("T trial", s0, &base(s0.seed), &HistoryRun { t_years: t, steps: N_PHYS, rebound: Some(rb), timing: Timing::History });
        let pk = peak_m(&p.zm);
        trials.push((t, pk));
        eprintln!("      T trial {:.3e} yr → peak {:.0} m", t, pk);
        pk
    };
    let grid: Vec<f64> = [5.0f64, 6.0, 7.0, 8.0, 9.0].iter().map(|e| 10f64.powf(*e)).collect();
    let vals: Vec<f32> = grid.iter().map(|&t| peak_at(t, &mut trials)).collect();
    let ok = |v: f32| (2700.0..=3000.0).contains(&v);
    let mut t_cal = grid[vals.iter().enumerate().min_by(|a, b| (a.1 - 2850.0).abs().partial_cmp(&(b.1 - 2850.0).abs()).unwrap()).unwrap().0];
    if let Some(i) = (0..4).find(|&i| (vals[i] - 2850.0) * (vals[i + 1] - 2850.0) <= 0.0) {
        let (mut lo, mut hi) = (grid[i].log10(), grid[i + 1].log10());
        let (mut vlo, mut vhi) = (vals[i], vals[i + 1]);
        if ok(vlo) {
            t_cal = grid[i];
        } else if ok(vhi) {
            t_cal = grid[i + 1];
        } else {
            for _ in 0..7 {
                let mid = 0.5 * (lo + hi);
                let v = peak_at(10f64.powf(mid), &mut trials);
                t_cal = 10f64.powf(mid);
                if ok(v) {
                    break;
                }
                if (v - 2850.0) * (vlo - 2850.0) > 0.0 {
                    lo = mid;
                    vlo = v;
                } else {
                    hi = mid;
                    vhi = v;
                }
            }
            let _ = vhi;
        }
    } else {
        eprintln!("      NOTE: no trial brackets 2 850 m; the nearest is kept");
    }
    let t_peak = trials.iter().find(|t| t.0 == t_cal).map_or(f32::NAN, |t| t.1);
    eprintln!("   → T = {:.3e} yr (peak {:.0} m; {} trials)", t_cal, t_peak, trials.len());
    // ── each seed: F164 steady, history, controls ──
    struct SeedOut {
        seed: u64,
        p: Vec<Phys>,
        hist_level: HistoryLevel,
    }
    let mut outs: Vec<SeedOut> = Vec::new();
    for s in &seeds {
        let cfg = base(s.seed);
        let t = Instant::now();
        let (rec, cal) = physics_level(s.snaps.last().unwrap(), 256, &cfg, 2850.0, 300, &mut |_| {}, &|| false).unwrap();
        let f164 = phys_of("F164 (steady)", to_m(&rec.z, k), cal.lifted_mask.clone(), t.elapsed().as_secs_f64(), None);
        let run = |rebound: Option<Rebound>, timing: Timing| HistoryRun { t_years: t_cal, steps: N_PHYS, rebound, timing };
        let (hist, hl) = history_phys("history", s, &cfg, &run(Some(rb), Timing::History));
        let (c2, _) = history_phys("control 2 (no rebound)", s, &cfg, &run(None, Timing::History));
        let (c3, _) = history_phys("control 3 (late 10 %)", s, &cfg, &run(Some(rb), Timing::Late10));
        let half = Rebound { alpha_km: 0.5 * rb.alpha_km, ..rb };
        let dbl = Rebound { alpha_km: 2.0 * rb.alpha_km, ..rb };
        let (sh, _) = history_phys("α × 0.5", s, &cfg, &run(Some(half), Timing::History));
        let (sd, _) = history_phys("α × 2", s, &cfg, &run(Some(dbl), Timing::History));
        eprintln!("\n   ── seed {} (256²):", s.seed);
        for p in [&f164, &hist, &c2, &c3, &sh, &sd] {
            let (r, _) = read(&p.zm, Some(&p.d5));
            let oc = ocean_mask(&p.zm);
            let interior_below = (0..256 * 256).filter(|&i| !oc[i] && p.zm.data[i] <= 0.0).count();
            print_macro(&p.name, &p.m);
            eprintln!(
                "      {:<22} peak {:.0} m · mean land {:.0} m · λ {:.2} km · density {:.3} km/km² · relief p50 {:.0} m · interior land ≤ 0 m {} cells ({:.2} % of land) · lakes {} / {:.0} km² · {:.1} s",
                p.name,
                r.peak,
                r.mean,
                r.lambda_km,
                r.density,
                r.relief_p50,
                interior_below,
                100.0 * interior_below as f32 / r.land.max(1) as f32,
                r.lakes_d5.0 + r.lakes_other.0,
                r.lakes_d5.1 + r.lakes_other.1,
                p.secs
            );
        }
        outs.push(SeedOut { seed: s.seed, p: vec![f164, hist, c2, c3, sh, sd], hist_level: hl });
    }
    // ── P2, P4, P6 at 256² ──
    eprintln!("\n   P2 / P4 / P6 (256²):");
    for (s, o) in seeds.iter().zip(&outs) {
        let up: Vec<Vec<f32>> = s.snaps.iter().map(|g| upsample_to(g, 256).data.iter().map(|v| (v - 0.5) * k).collect()).collect();
        let nsn = up.len();
        let (mut early, mut young) = (vec![false; 256 * 256], vec![false; 256 * 256]);
        for i in 0..256 * 256 {
            let tot = up[nsn - 1][i] - up[0][i];
            if tot < 500.0 {
                continue;
            }
            let first_half = up[nsn / 2][i] - up[0][i];
            let last_third = up[nsn - 1][i] - up[(2 * (nsn - 1)) / 3][i];
            early[i] = first_half > 0.5 * tot;
            young[i] = last_third > 0.5 * tot;
        }
        let maxin = |z: &GridF32, m: &[bool]| (0..256 * 256).filter(|&i| m[i]).map(|i| z.data[i]).fold(f32::NAN, f32::max);
        let (f164, hist) = (&o.p[0], &o.p[1]);
        let argmax = (0..256 * 256).max_by(|&a, &b| hist.zm.data[a].partial_cmp(&hist.zm.data[b]).unwrap()).unwrap();
        eprintln!(
            "   seed {}: P2 early belts {} cells (peak F164 {:.0} m → history {:.0} m, {:+.0} %) · young belts {} cells · the history's peak cell is {}",
            s.seed,
            early.iter().filter(|&&b| b).count(),
            maxin(&f164.zm, &early),
            maxin(&hist.zm, &early),
            100.0 * (maxin(&hist.zm, &early) / maxin(&f164.zm, &early) - 1.0),
            young.iter().filter(|&&b| b).count(),
            if young[argmax] { "YOUNG" } else if early[argmax] { "EARLY" } else { "neither" }
        );
        let crat = up_mask(&roll_mask(s.state.cratonic_mask.data(), 64, off), 64, 4);
        let p50 = |z: &GridF32| {
            let oc = ocean_mask(z);
            pct(&sorted((0..256 * 256).filter(|&i| crat[i] && !oc[i]).map(|i| z.data[i]).collect()), 0.5)
        };
        let (rh, _) = read(&hist.zm, Some(&hist.d5));
        let (r3, _) = read(&o.p[3].zm, Some(&o.p[3].d5));
        eprintln!(
            "      P4 cratons' land p50: history {:.0} m · no rebound {:.0} m ({:+.0} %) · P6 drainage density: history {:.3} · late 10 % {:.3} ({:.2}×)",
            p50(&hist.zm),
            p50(&o.p[2].zm),
            100.0 * (p50(&o.p[2].zm) / p50(&hist.zm) - 1.0),
            rh.density,
            r3.density,
            r3.density / rh.density
        );
    }
    // ── the frozen chains to 2 048² ──
    eprintln!("\n   D — the frozen cascade (p6 · W-tous W1, F164's k) to 2048²:");
    let ctrl = facet_control(2048, DOMAIN_KM / 2048.0);
    let ctrl_ok = ctrl.0 >= 60.0 && ctrl.1 <= 20.0;
    let mut finals: Vec<(u64, ChainOut, ChainOut, Vec<u8>, Vec<u8>)> = Vec::new();
    for o in &outs {
        let a = chain(&o.p[0].zm, &o.p[0].d5, frozen());
        let b = chain(&o.p[1].zm, &o.p[1].d5, frozen());
        eprintln!("   seed {} F164 ({:.0} s): {}", o.seed, a.secs, a.notes.join(" | "));
        eprintln!("   seed {} history ({:.0} s): {}", o.seed, b.secs, b.notes.join(" | "));
        finals.push((o.seed, a, b, o.p[0].classes.clone(), o.p[1].classes.clone()));
    }
    // ── texture, coast, success ──
    let ckm = DOMAIN_KM / 2048.0;
    let mut successes = 0usize;
    for (seed, a, b, ca, cb) in &finals {
        let o = outs.iter().find(|o| o.seed == *seed).unwrap();
        eprintln!("\n   ── seed {seed} at 2048² (195 m) against Corsica:");
        let mut p5 = true;
        for (name, run, cls) in [("F164", a, ca), ("history", b, cb)] {
            let whole = tex(&run.g, ckm, &run.d5);
            let m8 = up_u8(cls, 256, 8);
            let excl: Vec<bool> = (0..2048 * 2048).map(|i| run.d5[i] || m8[i] != MOUNTAIN).collect();
            let mtn = tex(&run.g, ckm, &excl);
            for (tag, x, c) in [("whole", &whole, &c_whole), ("mountain class", &mtn, &c_mtn)] {
                let el = elements(x, c, ckm, ctrl_ok);
                let fails: Vec<String> = el.iter().filter(|e| !e.3).map(|e| format!("{} {:.3}/{:.3}", e.0, e.1, e.2)).collect();
                eprintln!("      {name:<8} {tag:<15} {} of {} within ×1.5 · out: {}", el.len() - fails.len(), el.len(), if fails.is_empty() { "none".into() } else { fails.join("; ") });
                if name == "history" && tag == "mountain class" {
                    p5 = el.iter().filter(|e| e.0 != "A_dir" && e.0 != "F").all(|e| e.3);
                }
            }
            // F per class
            let (rd, fl) = ymir_core::cascade::planform::flat_map(&run.g, ckm, Some(&run.d5));
            let per: Vec<String> = (1..5u8)
                .map(|c| {
                    let cells = (0..2048 * 2048).filter(|&i| rd[i] && m8[i] == c).count();
                    let f = (0..2048 * 2048).filter(|&i| rd[i] && fl[i] && m8[i] == c).count();
                    format!("{} {:.2} %", ymir_core::cascade::physio::CLASS_NAMES[c as usize], 100.0 * f as f32 / cells.max(1) as f32)
                })
                .collect();
            eprintln!("      {name:<8} F by class: {}", per.join(" · "));
        }
        let land_a: Vec<bool> = ocean_mask(&a.g).iter().map(|o| !o).collect();
        let land_b: Vec<bool> = ocean_mask(&b.g).iter().map(|o| !o).collect();
        let (ka, kb) = (coast(&land_a, 2048, 2048, 32), coast(&land_b, 2048, 2048, 32));
        eprintln!("      coast F164 {ka:?}\n      coast history {kb:?}\n      coast Corsica {c_coast:?}");
        let degraded = |ma: f32, mb: f32, mc: f32| (mb - mc).abs() > (ma - mc).abs() + 0.10 * ma.abs();
        let deg: Vec<&str> = [
            ("dimension", degraded(ka.dimension, kb.dimension, c_coast.dimension)),
            ("length ratio", degraded(ka.length_ratio, kb.length_ratio, c_coast.length_ratio)),
            ("aligned", degraded(ka.aligned, kb.aligned, c_coast.aligned)),
            ("block index", degraded(ka.block_index, kb.block_index, c_coast.block_index)),
        ]
        .iter()
        .filter(|x| x.1)
        .map(|x| x.0)
        .collect();
        let na = land_a.iter().filter(|&&l| l).count() as f32;
        let gained = (0..2048 * 2048).filter(|&i| land_b[i] && !land_a[i]).count() as f32;
        let lost = (0..2048 * 2048).filter(|&i| land_a[i] && !land_b[i]).count() as f32;
        eprintln!("      coast degraded (> 10 %): {} · displacement at 2048²: gained {:.2} % · lost {:.2} % of F164's land", if deg.is_empty() { "none".into() } else { deg.join(", ") }, 100.0 * gained / na, 100.0 * lost / na);
        let (fa, fb) = (&o.p[0], &o.p[1]);
        let la: Vec<bool> = ocean_mask(&fa.zm).iter().map(|o| !o).collect();
        let lb: Vec<bool> = ocean_mask(&fb.zm).iter().map(|o| !o).collect();
        let n256 = la.iter().filter(|&&l| l).count() as f32;
        eprintln!(
            "      displacement at 256²: gained {:.2} % · lost {:.2} %",
            100.0 * (0..65536).filter(|&i| lb[i] && !la[i]).count() as f32 / n256,
            100.0 * (0..65536).filter(|&i| la[i] && !lb[i]).count() as f32 / n256
        );
        let m = &fb.m;
        let within = |a: f32, b: f32, f: f32| a > 0.0 && b > 0.0 && a <= f * b && a >= b / f;
        let fr = [
            within(m.fractions[3], me.fractions[3], 1.5),
            within(m.fractions[1], me.fractions[1], 1.5),
            within(m.fractions[0], me.fractions[0], 1.5),
        ];
        let ratio_ok = within(m.mtn_ratio, me.mtn_ratio, 2.0);
        let succ = fr.iter().all(|&b| b) && ratio_ok && p5 && deg.is_empty();
        if succ {
            successes += 1;
        }
        eprintln!(
            "      SUCCESS test: mountain {} · plateau {} · plain {} (×1.5 of Europe) · width ratio {} (×2) · P5 {} · coast {} → {}",
            fr[0], fr[1], fr[2], ratio_ok, p5, deg.is_empty(), if succ { "SUCCESS" } else { "no" }
        );
    }
    let plateau_low = outs.iter().filter(|o| o.p[1].m.fractions[1] < 0.5 * me.fractions[1]).count();
    let p1_holds = a2_widths.iter().filter(|w| w.0 <= 3.0 && w.1 >= 0.7).count();
    eprintln!(
        "\n   R: success on {successes} / 4 seeds · P1 holds on {p1_holds} / 4 seeds · the history's plateau fraction under 50 % of Europe's ({:.1} %) on {plateau_low} / 4 seeds",
        100.0 * 0.5 * me.fractions[1]
    );
    // ── images ──
    let cell_m = CELL * 1000.0;
    let mut r256 = Vec::new();
    let mut r2048 = Vec::new();
    let mut rcls = Vec::new();
    for ((seed, a, b, ca, cb), o) in finals.iter().zip(&outs) {
        let _ = seed;
        let sh = |z: &GridF32, n: usize, px: usize| {
            let (_, hy) = read(z, None);
            whole(&rgba_of(z, Some(&hy), DOMAIN_KM / n as f32 * 1000.0, &ss), n, px)
        };
        r256.push(vec![sh(&o.p[0].zm, 256, 384), sh(&o.p[1].zm, 256, 384)]);
        r2048.push(vec![sh(&a.g, 2048, 512), sh(&b.g, 2048, 512)]);
        rcls.push(vec![class_png(ca, 256, 256, 1), class_png(cb, 256, 256, 1)]);
        let _ = cell_m;
    }
    compose(&r256).save(out().join("f165_seeds_256.png")).unwrap();
    compose(&r2048).save(out().join("f165_seeds_2048.png")).unwrap();
    compose(&rcls).save(out().join("f165_classes_seeds.png")).unwrap();
    // three snapshots, Σ uplift / erosion / rebound (the témoin)
    let hl = &outs[0].hist_level;
    let snap_png = |g: &GridF32| {
        let zm = to_m(g, k);
        whole(&rgba_of(&zm, None, cell_m, &ss), 256, 384)
    };
    let ns = hl.snapshots.len();
    compose(&[vec![snap_png(&hl.snapshots[0]), snap_png(&hl.snapshots[ns / 2]), snap_png(&hl.snapshots[ns - 1])]]).save(out().join("f165_snapshots.png")).unwrap();
    if let Some((u, e, r)) = &outs[0].p[1].rec_sums {
        let tile = |d: &[f32]| {
            let mut a: Vec<f32> = d.iter().map(|v| v.abs()).collect();
            a.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let sat = a[((a.len() - 1) as f32 * 0.98) as usize].max(1.0);
            whole(&diff_rgba(d, 256, 256, sat), 256, 384)
        };
        compose(&[vec![tile(u), tile(e), tile(r)]]).save(out().join("f165_sums.png")).unwrap();
    }
    // the coast crops (25 km = 128 cells at 195 m): the first window of the témoin with 40–60 % land, and Corsica's
    let pick = |g: &GridF32, n: usize| -> (usize, usize) {
        let land: Vec<bool> = ocean_mask(g).iter().map(|o| !o).collect();
        for oy in (0..n - 128).step_by(64) {
            for ox in (0..n - 128).step_by(64) {
                let c = (0..128 * 128).filter(|&i| land[(oy + i / 128) * n + ox + i % 128]).count() as f32 / (128.0 * 128.0);
                if (0.4..=0.6).contains(&c) {
                    return (ox, oy);
                }
            }
        }
        (0, 0)
    };
    let (fa, fb) = (&finals[0].1.g, &finals[0].2.g);
    let (ox, oy) = pick(fa, 2048);
    let (cx, cy) = pick(&cg, 1024);
    let crop = |g: &GridF32, n: usize, x: usize, y: usize| {
        let (_, hy) = read(g, None);
        let rg = rgba_of(g, Some(&hy), 195.3, &ss);
        window(&rg, n, (x as f32 / n as f32, y as f32 / n as f32), 128.0 / n as f32, 448)
    };
    compose(&[vec![crop(&cg, 1024, cx, cy), crop(fa, 2048, ox, oy), crop(fb, 2048, ox, oy)]]).save(out().join("f165_coast_25km.png")).unwrap();
    eprintln!(
        "   images: f165_seeds_256.png and f165_seeds_2048.png (rows = seeds {SEEDS:?}; columns F164, history), f165_classes_seeds.png (same layout), f165_snapshots.png (témoin h_iso at t₀, mid, final), f165_sums.png (témoin Σ uplift, Σ erosion, Σ rebound; red +, blue −), f165_coast_25km.png (Corsica, F164, history; témoin window at ({ox}, {oy}) of 2048², Corsica at ({cx}, {cy}) of 1024²) · north up · {NOTICE}"
    );
    eprintln!("\n==========  end Finding 165 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}

/// F165-A2, the physics-level half (declared, missing from the first run): the shapes above 500 / 1 000 / 2 000 m on
/// the current physics level (F164's regime, 256², 1.5625 km cells), near the convergent cells or sutures within 8 cells
/// (2 coarse cells), beside h_iso's at 64².
#[test]
#[ignore]
fn f165_a2_physics() {
    let ss = SteinSteinParams::default();
    let k = n2m();
    let temoin = Knobs { valley: Some(ValleyConstruction::new(F121_AGE_K, Some(0.1))), slope_floor_abs: Some(0.024), ..Knobs::passes(2) };
    let tlf = bench_upscale_cfg(temoin).target_land_fraction;
    let rough_a = 0.5 * read_cell(&load_grid(&root().join("data/corsica/corse_128.bin")).unwrap(), None, super::CELL).0.octaves[0];
    let off = [(CANONICAL_ORIGIN[0] * 64.0).round() as usize, (CANONICAL_ORIGIN[1] * 64.0).round() as usize];
    eprintln!("\n   A2 — the physics level's shapes (F164's regime, 256²), beside h_iso (64²):");
    for seed in SEEDS {
        let s = run_seed(seed, tlf);
        let st = &s.state;
        let bi = classify_boundaries(&st.plate_id, &s.kin);
        let conv: Vec<bool> = bi.boundary_type.data().iter().map(|b| *b == BoundaryType::Convergent).collect();
        let near64 = roll_mask(&(0..64 * 64).map(|i| conv[i] || s.hist.suture[i]).collect::<Vec<_>>(), 64, off);
        let near256 = up_mask(&near64, 64, 4);
        let cfg = CascadeConfig { seed, roughness_m: rough_a, mfd_exponent: Some(6.0), ..CascadeConfig::declared(ss.depth_scale_m as f32) };
        let (rec, _) = physics_level(s.snaps.last().unwrap(), 256, &cfg, 2850.0, 300, &mut |_| {}, &|| false).unwrap();
        let zm = to_m(&rec.z, k);
        let hm: Vec<f32> = s.snaps.last().unwrap().data.iter().map(|v| (v - 0.5) * k).collect();
        for thr in [500.0f32, 1000.0, 2000.0] {
            let a = shape(&hm.iter().map(|&v| v > thr).collect::<Vec<_>>(), 64, 64, Some(&near64), 2);
            let b = shape(&zm.data.iter().map(|&v| v > thr).collect::<Vec<_>>(), 256, 256, Some(&near256), 8);
            eprintln!(
                "   seed {seed} > {thr:.0} m: h_iso width {:.1} km ({} comp., near {:.0} %) · physics width {:.1} km ({} comp., {} cells, near {:.0} %) · ratio h_iso / physics {:.2}",
                a.width_cells * 6.25,
                a.components,
                100.0 * a.near_share,
                b.width_cells * super::CELL,
                b.components,
                b.cells,
                100.0 * b.near_share,
                (a.width_cells * 6.25) / (b.width_cells * super::CELL)
            );
        }
    }
}
