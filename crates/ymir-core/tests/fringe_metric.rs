//! ADR Finding 52 — the turn-rate metric is a function of VERTEX SPACING, and a metric that
//! measures what the author sees.
//!
//! PART A settles a contradiction before anything is built on it. Two measurements of "the
//! coarse field" disagree: 43.5 % of turns above 80° tracing the RAW 64² grid (172 vertices),
//! against 0.35 % tracing the same field UPSCALED to 2048² (Finding 51, stage 1). If both are
//! right, the turn rate depends on how finely the contour is sampled and not on its shape — and
//! then "the coarse field is already at 43.5 %" cannot support "the fringes come from
//! tectonics", because the number would say the same thing about any coastline sampled coarsely.
//! The test: take ONE fixed coastline and resample it at decreasing vertex densities.
//!
//! PART B builds the metric the turn count cannot be: fringe LENGTH and fringe REGULARITY, in
//! kilometres, at an explicit physical scale.
//!
//! Run: cargo test -p ymir-core --release --test fringe_metric -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::contour::{Polyline, marching_squares};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;

fn coarse_state() -> (ymir_core::tectonics_c1::state::C1State, PlateKinematics, C1TimeLoopConfig) {
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    (state, kin, run_cfg)
}

/// The production HD field. `closures` = C-3 + C-3b on (the shipped product).
fn hd_field(target: usize, closures: bool) -> GridF32 {
    hd_stage(target, closures, 0.04, 3)
}

/// `incision`: 0 = none, 2 = relief-v2, 3 = relief-v3 (the shipped one).
fn hd_stage(target: usize, closures: bool, amplitude: f64, incision: u8) -> GridF32 {
    let ss = SteinSteinParams::default();
    let (state, kin, run_cfg) = coarse_state();
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: amplitude,
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: closures,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: closures,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    match incision {
        0 => cfg.stream_power = None,
        2 => {
            cfg.stream_power = Some(ymir_core::erosion::stream_power::StreamPowerConfig::relief_v2(
                (DOMAIN_KM / target as f32).powi(2),
                ss.depth_scale_m as f32,
            ))
        }
        _ => {}
    }
    let (up, _) = upscale_from_c1_with_progress(
        &state,
        &run_cfg.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        Some(&kin),
        &mut |_| {},
        &|| false,
    );
    up.heightmap
}

fn turn_rate(polys: &[Polyline]) -> (usize, usize, f32) {
    let (mut n, mut k) = (0usize, 0usize);
    for pl in polys {
        for w3 in pl.windows(3) {
            let (ax, ay) = (w3[1].0 - w3[0].0, w3[1].1 - w3[0].1);
            let (bx, by) = (w3[2].0 - w3[1].0, w3[2].1 - w3[1].1);
            let (la, lb) = ((ax * ax + ay * ay).sqrt(), (bx * bx + by * by).sqrt());
            if la <= 0.0 || lb <= 0.0 {
                continue;
            }
            let t = (((ax * bx + ay * by) / (la * lb)).clamp(-1.0, 1.0)).acos().to_degrees();
            n += 1;
            if t > 80.0 {
                k += 1;
            }
        }
    }
    (n, k, 100.0 * k as f32 / n.max(1) as f32)
}

/// Keep every `step`-th vertex — the SAME curve, sampled more coarsely.
fn decimate(polys: &[Polyline], step: usize) -> Vec<Polyline> {
    polys
        .iter()
        .map(|pl| {
            let mut out: Vec<(f32, f32)> = pl.iter().step_by(step).copied().collect();
            if pl.first() == pl.last() && out.len() > 2 && out.first() != out.last() {
                out.push(out[0]);
            }
            out
        })
        .filter(|pl| pl.len() >= 3)
        .collect()
}

// ─────────────────────────────── PART A ────────────────────────────────

#[test]
#[ignore]
fn turn_rate_is_a_function_of_vertex_spacing() {
    let ss = SteinSteinParams::default();
    let (state, _kin, run_cfg) = coarse_state();

    // The RAW 64² isostatic altitude — the field before any upscale.
    // The SAME normalisation production uses (`target_land_fraction: None` is the shipped
    // setting — `production_hd_config` never sets it, ADR Finding 43).
    let coarse = ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude(
        &state,
        &run_cfg.iso_config,
        &ss,
        None,
    );
    let cpolys = marching_squares(&coarse, SEA);
    let cverts: usize = cpolys.iter().map(|p| p.len()).sum();
    let (cn, ck, crate_) = turn_rate(&cpolys);
    eprintln!(
        "\n=====  PART A — is the turn rate a shape, or a sampling artefact?  =====\n\n\
         RAW 64² contour: {} rings, {cverts} vertices, {ck}/{cn} turns > 80° = {crate_:.1} %\n\
         (a 64² grid over {DOMAIN_KM} km is 6.25 km per cell, so consecutive vertices are ~6 km \
         apart)",
        cpolys.len()
    );

    // The SAME continent at 2048², production, decimated to progressively coarser spacings.
    let hd = hd_field(2048, true);
    let hpolys = marching_squares(&hd, SEA);
    let km_per_cell = DOMAIN_KM / 2048.0;
    eprintln!(
        "\nONE FIXED 2048² COASTLINE, resampled — the shape never changes, only the spacing:\n\
         {:>6} {:>10} {:>12} {:>12}",
        "step", "vertices", "spacing km", "turns > 80°"
    );
    for step in [1usize, 2, 4, 8, 16, 32, 64] {
        let d = decimate(&hpolys, step);
        let v: usize = d.iter().map(|p| p.len()).sum();
        let (_, _, r) = turn_rate(&d);
        eprintln!(
            "{:>6} {:>10} {:>12.2} {:>11.1}%",
            step,
            v,
            step as f32 * 0.7 * km_per_cell, // ~0.7 cell mean step × decimation
            r
        );
    }
    eprintln!(
        "\n⇒ If the rate climbs toward the coarse figure as the spacing approaches 6 km, then\n  \
         '43.5 % at 64²' and '0.35 % upscaled' are THE SAME COASTLINE measured two ways, and the\n  \
         number cannot attribute the fringes to any stage."
    );
}

// ─────────────────────────────── PART B ────────────────────────────────
// The metric the turn count cannot be. Three quantities, all in KILOMETRES or dimensionless,
// all at an EXPLICIT physical scale, chosen to match what the eye distinguishes:
//
//   1. SPUR LENGTH   — how far a fringe runs out from the shore body. Found as a "neck": two
//                      points close in SPACE but far apart ALONG the coastline. The arc between
//                      them is the spur. This is what "long fringes" means, and a turn count
//                      cannot see it — a spur of any length has the same two sharp turns.
//   2. REGULARITY    — the coefficient of variation of the spacing between consecutive spurs.
//                      A COMB is regular (low CV); an organic ria coast is irregular (high CV).
//   3. PARALLELISM   — the circular concentration R of the spur AXES. Parallel spurs read as
//                      manufactured; scattered axes read as natural.

struct Spurs {
    coast_km: f32,
    count: usize,
    med_len_km: f32,
    p90_len_km: f32,
    /// Share of the whole coastline that lies inside a spur.
    share_pct: f32,
    /// Coefficient of variation of the spacing between consecutive spur roots.
    spacing_cv: f32,
    /// Circular concentration of the spur axes (axial, doubled-angle). 1 = all parallel.
    ///
    /// ⚠️ GLOBAL, and therefore BLIND to the actual defect. The renders show the spurs arranged
    /// in BUNDLES, each bundle internally parallel but pointing a different way from the next —
    /// so the global resultant cancels to ~0.03 and reads "isotropic" while the picture shows a
    /// comb. Use `local_axis_r` to judge parallelism.
    axis_r: f32,
    /// MEDIAN over spurs of the axial concentration among each spur and its 8 nearest
    /// neighbours. This is what "parallel fringes" means: locally aligned, globally scattered.
    local_axis_r: f32,
}

/// A spur is an excursion that comes back on itself: an arc of at least `min_spur_km` whose
/// endpoints are within `neck_km` of each other. `neck_km` is the width of the isthmus that
/// makes it a spur rather than a bay.
/// ⚠️ `MAX_SPUR_KM` is the search WINDOW, and it is an instrument ceiling: a spur longer than
/// this cannot be found, so every length statistic saturates at it. The first version used
/// `min_spur_km * 8` = 8 km and the `max` column read 7.99–8.00 km at EVERY setting including
/// the reference — a constant, not a measurement. The `max` column is what exposed it, which is
/// the argument for reporting a tail rather than a median alone.
const MAX_SPUR_KM: f32 = 50.0;

fn spurs(polys: &[Polyline], km_per_cell: f32, min_spur_km: f32, neck_km: f32) -> Spurs {
    let neck_cells = neck_km / km_per_cell;
    let min_arc_cells = min_spur_km / km_per_cell;
    let mut lens: Vec<f32> = Vec::new();
    let mut axes: Vec<(f32, f32, f32)> = Vec::new(); // (mid x, mid y, axis angle)
    let mut roots: Vec<f32> = Vec::new(); // arc position of each spur root, for the spacing CV
    let (mut c2, mut s2, mut naxis) = (0.0f64, 0.0f64, 0usize);
    let mut coast_cells = 0.0f32;
    let mut in_spur_cells = 0.0f32;

    for pl in polys {
        // cumulative arc length along the ring
        let mut arc = vec![0.0f32; pl.len()];
        for i in 1..pl.len() {
            arc[i] = arc[i - 1]
                + ((pl[i].0 - pl[i - 1].0).powi(2) + (pl[i].1 - pl[i - 1].1).powi(2)).sqrt();
        }
        coast_cells += *arc.last().unwrap_or(&0.0);
        if pl.len() < 8 {
            continue;
        }
        // Walk forward; at each i find the FARTHEST j that closes a neck, then skip past it so
        // spurs do not overlap.
        let mut i = 0usize;
        while i < pl.len() {
            let mut best: Option<usize> = None;
            let mut j = i + 1;
            while j < pl.len() && (arc[j] - arc[i]) * km_per_cell < MAX_SPUR_KM {
                if arc[j] - arc[i] >= min_arc_cells {
                    let d = ((pl[j].0 - pl[i].0).powi(2) + (pl[j].1 - pl[i].1).powi(2)).sqrt();
                    if d <= neck_cells {
                        best = Some(j);
                    }
                }
                j += 1;
            }
            match best {
                Some(j) => {
                    let len = (arc[j] - arc[i]) * km_per_cell;
                    lens.push(len);
                    roots.push(arc[i] * km_per_cell);
                    in_spur_cells += arc[j] - arc[i];
                    // Axis: from the neck midpoint to the farthest point on the sub-path.
                    let mid = ((pl[i].0 + pl[j].0) * 0.5, (pl[i].1 + pl[j].1) * 0.5);
                    let mut far = (0.0f32, (0.0f32, 0.0f32));
                    for p in &pl[i..=j] {
                        let d = ((p.0 - mid.0).powi(2) + (p.1 - mid.1).powi(2)).sqrt();
                        if d > far.0 {
                            far = (d, *p);
                        }
                    }
                    let th = ((far.1.1 - mid.1) as f64).atan2((far.1.0 - mid.0) as f64);
                    c2 += (2.0 * th).cos();
                    s2 += (2.0 * th).sin();
                    naxis += 1;
                    axes.push((mid.0, mid.1, th as f32));
                    i = j; // no overlap
                }
                None => i += 1,
            }
        }
    }
    lens.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q = |f: f32| -> f32 {
        if lens.is_empty() {
            0.0
        } else {
            lens[(((lens.len() - 1) as f64 * f as f64).round() as usize).min(lens.len() - 1)]
        }
    };
    // Spacing CV between consecutive spur roots (per ring order, which is arc order).
    let mut gaps: Vec<f32> = Vec::new();
    for w in roots.windows(2) {
        if w[1] > w[0] {
            gaps.push(w[1] - w[0]);
        }
    }
    let (mean, cv) = if gaps.is_empty() {
        (0.0, 0.0)
    } else {
        let m = gaps.iter().sum::<f32>() / gaps.len() as f32;
        let var = gaps.iter().map(|g| (g - m) * (g - m)).sum::<f32>() / gaps.len() as f32;
        (m, if m > 0.0 { var.sqrt() / m } else { 0.0 })
    };
    let _ = mean;
    // LOCAL parallelism: for each spur, the axial resultant over itself and its 8 nearest
    // neighbours by midpoint. Bundles that are internally parallel score high here even when
    // the global resultant is ~0.
    const K: usize = 8;
    let mut local: Vec<f32> = Vec::new();
    for (i, a) in axes.iter().enumerate() {
        let mut d: Vec<(f32, usize)> = axes
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(j, b)| (((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt(), j))
            .collect();
        if d.len() < K {
            continue;
        }
        d.select_nth_unstable_by(K - 1, |x, y| {
            x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal)
        });
        let (mut cc, mut sc) = (2.0 * a.2, 0.0f32);
        cc = cc.cos();
        sc = (2.0 * a.2).sin();
        let (mut cs, mut ss_) = (cc, sc);
        for (_, j) in d.iter().take(K) {
            cs += (2.0 * axes[*j].2).cos();
            ss_ += (2.0 * axes[*j].2).sin();
        }
        let n = (K + 1) as f32;
        local.push(((cs / n).powi(2) + (ss_ / n).powi(2)).sqrt());
    }
    local.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let local_axis_r = if local.is_empty() { 0.0 } else { local[local.len() / 2] };

    Spurs {
        coast_km: coast_cells * km_per_cell,
        count: lens.len(),
        med_len_km: q(0.5),
        p90_len_km: q(0.9),
        share_pct: 100.0 * in_spur_cells / coast_cells.max(1e-6),
        spacing_cv: cv,
        axis_r: if naxis == 0 {
            0.0
        } else {
            (((c2 / naxis as f64).powi(2) + (s2 / naxis as f64).powi(2)).sqrt()) as f32
        },
        local_axis_r,
    }
}

fn report(label: &str, s: &Spurs) {
    eprintln!(
        "{label:<34} {:>9.0} {:>8} {:>10.2} {:>10.2} {:>9.1}% {:>10.2} {:>9.3}",
        s.coast_km, s.count, s.med_len_km, s.p90_len_km, s.share_pct, s.spacing_cv, s.axis_r
    );
}

fn header() {
    eprintln!(
        "{:<34} {:>9} {:>8} {:>10} {:>10} {:>10} {:>10} {:>9}",
        "config", "coast km", "spurs", "med km", "p90 km", "in-spur", "spacing CV", "axis R"
    );
}

#[test]
#[ignore]
fn fringe_metric_ranks_the_configs() {
    let km2048 = DOMAIN_KM / 2048.0;
    // Scale choice, stated: a spur must run at least 1 km out-and-back and close through a neck
    // no wider than 0.6 km. Those are the dimensions the author points at in the export.
    let (min_spur_km, neck_km) = (1.0f32, 0.6f32);
    eprintln!(
        "\n=====  PART B — FRINGE LENGTH AND REGULARITY (2048², spur ≥ {min_spur_km} km, neck ≤ {neck_km} km)  ====="
    );
    eprintln!(
        "\nVALIDATION TARGET: the author judges the CURRENT export worse than the pre-chantier\n\
         state. Pre-chantier = the same relief with C-3 and C-3b OFF. If the metric cannot put\n\
         those two in that order it is not ready to judge a remedy.\n"
    );
    header();
    let off = marching_squares(&hd_field(2048, false), SEA);
    let on = marching_squares(&hd_field(2048, true), SEA);
    let s_off = spurs(&off, km2048, min_spur_km, neck_km);
    let s_on = spurs(&on, km2048, min_spur_km, neck_km);
    report("pre-chantier (C-3, C-3b OFF)", &s_off);
    report("current export (both ON)", &s_on);
    eprintln!(
        "\n{:<34} {:>9.0} {:>8} {:>10.2} {:>10.2} {:>9.1} {:>10.2} {:>9.3}   <- delta",
        "",
        s_on.coast_km - s_off.coast_km,
        s_on.count as i64 - s_off.count as i64,
        s_on.med_len_km - s_off.med_len_km,
        s_on.p90_len_km - s_off.p90_len_km,
        s_on.share_pct - s_off.share_pct,
        s_on.spacing_cv - s_off.spacing_cv,
        s_on.axis_r - s_off.axis_r
    );
    eprintln!(
        "\nVERDICT ON THE METRIC: it ranks the two in the author's order if the spur COUNT and\n\
         the IN-SPUR SHARE rise. A falling spacing CV would say the fringes also became more\n\
         REGULAR (comb-like); a rising axis R would say they became more PARALLEL."
    );
}

/// PART C — where do the SPURS enter? Same stage bisection, now on a metric that measures the
/// symptom. This is the "defect or property" question: if the upscaled coarse field already has
/// long spurs, the upscale manufactures them from an angular 64² contour and the remedy belongs
/// in the interpolation. If it does not, the coarse contour's angularity is a PROPERTY of a 64²
/// grid that the upscale smooths away, and the spurs come from further down the chain.
#[test]
#[ignore]
fn spur_stage_bisection() {
    let km = DOMAIN_KM / 2048.0;
    let (min_spur_km, neck_km) = (1.0f32, 0.6f32);
    eprintln!(
        "
=====  PART C — WHERE DO THE SPURS ENTER? (2048², spur >= {min_spur_km} km)  ====="
    );
    header();
    for (label, amp, inc, clos) in [
        ("1 coarse upscaled (no FBM)", 0.0f64, 0u8, false),
        ("2 + FBM", 0.04, 0, false),
        ("3 + relief-v2", 0.04, 2, false),
        ("4 + relief-v3", 0.04, 3, false),
        ("5 final (+ C-3 + C-3b)", 0.04, 3, true),
    ] {
        let polys = marching_squares(&hd_stage(2048, clos, amp, inc), SEA);
        report(label, &spurs(&polys, km, min_spur_km, neck_km));
    }
    eprintln!(
        "
⇒ The stage where the SPUR COUNT and the IN-SPUR SHARE jump is where the fringes are
           created. Compare with Finding 51, which located the TURN-rate jump at the same place —
           if the two agree, the earlier attribution survives its metric being replaced."
    );
}

// ─────────────────────────────── PART D ────────────────────────────────
// THE WARP SWEEP, and the anchoring the metric was missing.
//
// ANCHORING — what an organic coast looks like, in reference values that are DERIVED rather
// than invented. The author's requirement was "indentations of VARIED lengths and regularity
// near zero"; both have exact mathematical references:
//
//   • SPACING REGULARITY. For a POISSON process — events placed completely at random along the
//     shore — the coefficient of variation of the gaps is EXACTLY 1.0. So CV ≈ 1 means "as
//     irregular as chance"; CV ≪ 1 is a comb (evenly spaced teeth); CV > 1 is clustered.
//     TARGET: CV ≈ 1, and certainly not below ~0.6.
//   • LENGTH VARIETY. For an exponential length distribution — the maximum-entropy choice for a
//     positive quantity with a given mean, i.e. "as varied as it can be without further
//     structure" — p90/median = ln(10)/ln(2) = 3.32. A ratio near 1 means every spur is the
//     SAME length, which is the signature of a manufactured comb.
//     TARGET: p90/median ≈ 3.3; a ratio under ~2 is suspiciously uniform.
//   • PARALLELISM. Axial concentration R = 0 is isotropic. TARGET: R < 0.1.
//
// WHICH METRIC IS BEING OPTIMISED, stated so the recommendation cannot be misread:
//   OPTIMISED → the LENGTH VARIETY ratio p90/median, towards 3.3, and the spacing CV towards 1.
//   REPORTED ONLY → coastline length and spur count. Reducing the warp shortens the coast
//   MECHANICALLY, and a shorter coast is NOT a better one — a too-smooth trace is as wrong as a
//   combed one, and the author has already rejected one.

const EXP_P90_OVER_MEDIAN: f32 = 3.322; // ln(10)/ln(2)
const POISSON_CV: f32 = 1.0;

fn hd_warp(target: usize, warp: f64) -> GridF32 {
    let ss = SteinSteinParams::default();
    let (state, kin, run_cfg) = coarse_state();
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: 0.04,
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: true,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: true,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    cfg.coast_warp_strength = warp;
    let (up, _) = upscale_from_c1_with_progress(
        &state,
        &run_cfg.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        Some(&kin),
        &mut |_| {},
        &|| false,
    );
    up.heightmap
}

fn sweep_row(label: &str, s: &Spurs) {
    let variety = if s.med_len_km > 0.0 { s.p90_len_km / s.med_len_km } else { 0.0 };
    let flag = |ok: bool| if ok { " " } else { "!" };
    eprintln!(
        "{label:<26} {:>9.0} {:>7} {:>8.2} {:>8.2} {:>8.2}{} {:>8.2}{} {:>7.3}{}",
        s.coast_km,
        s.count,
        s.med_len_km,
        s.p90_len_km,
        variety,
        flag(variety >= 2.0),
        s.spacing_cv,
        flag(s.spacing_cv >= 0.6),
        s.axis_r,
        flag(s.axis_r < 0.1),
    );
}

fn sweep_header() {
    eprintln!(
        "{:<26} {:>9} {:>7} {:>8} {:>8} {:>9} {:>9} {:>8}",
        "coast_warp_strength", "coast km", "spurs", "med km", "p90 km", "VARIETY", "CV", "axis R"
    );
    eprintln!(
        "{:<26} {:>9} {:>7} {:>8} {:>8} {:>9} {:>9} {:>8}",
        "  ORGANIC TARGET", "-", "-", "-", "-", ">= 3.3", "~ 1.0", "< 0.1"
    );
}

#[test]
#[ignore]
fn warp_strength_sweep() {
    let (min_spur_km, neck_km) = (1.0f32, 0.6f32);
    eprintln!(
        "\n=====  PART D — coast_warp_strength SWEEP  =====\n\n\
         OPTIMISED: length VARIETY (p90/median, exponential reference {EXP_P90_OVER_MEDIAN:.2}) \
         and spacing CV\n            (Poisson reference {POISSON_CV:.1}). REPORTED ONLY: coastline \
         length and spur count —\n            a shorter coast is not a better one.\n\
         '!' marks a value outside the organic range."
    );
    for target in [2048usize, 8192] {
        let km = DOMAIN_KM / target as f32;
        eprintln!("\n--- {target}² ({:.0} m/cell) ---", km * 1000.0);
        sweep_header();
        for warp in [1.5f64, 1.0, 0.5, 0.25, 0.0] {
            let polys = marching_squares(&hd_warp(target, warp), SEA);
            let s = spurs(&polys, km, min_spur_km, neck_km);
            let tag = if (warp - 1.5).abs() < 1e-9 {
                format!("{warp:.2} (SHIPPED)")
            } else {
                format!("{warp:.2}")
            };
            sweep_row(&tag, &s);
        }
    }
}

/// PART E — the CONTOUR RELAXATION, re-measured with the metric that sees the symptom.
///
/// Finding 48 evaluated it with the TURN COUNT and got −0.4 %, a number with no meaning: the
/// turn count cannot see spur length, which is the thing the relaxation would have to shorten.
/// The relaxation only post-processes the polyline, so all four settings come from ONE terrain
/// build — the sweep is nearly free, and it closes the question properly rather than on an
/// instrument now known to be unfit.
#[test]
#[ignore]
fn relaxation_passes_sweep() {
    use ymir_core::terrain::contour::{
        SMOOTH_RELAX_BELOW_DEG, slope_deg_to_norm_gradient, smooth_polylines_on_isoline,
    };
    let (min_spur_km, neck_km) = (1.0f32, 0.6f32);
    eprintln!(
        "
=====  PART E — contour relaxation passes, on the spur metric  =====

         Same OPTIMISED / REPORTED split as Part D. The relaxation would have to raise the
         length VARIETY towards 3.3 to be worth anything; Finding 48 could not see that."
    );
    for target in [2048usize, 8192] {
        let km = DOMAIN_KM / target as f32;
        let ss = SteinSteinParams::default();
        let field = hd_field(target, true);
        let raw = marching_squares(&field, SEA);
        let gate = slope_deg_to_norm_gradient(
            SMOOTH_RELAX_BELOW_DEG,
            km * 1000.0,
            2.0 * 1.13 * ss.depth_scale_m as f32,
        );
        eprintln!(
            "
--- {target}² ({:.0} m/cell) ---",
            km * 1000.0
        );
        sweep_header();
        for passes in [0usize, 1, 2, 4] {
            let polys = if passes == 0 {
                raw.clone()
            } else {
                smooth_polylines_on_isoline(&field, SEA, &raw, passes, gate)
            };
            let s = spurs(&polys, km, min_spur_km, neck_km);
            let tag = if passes == 0 {
                "0 (SHIPPED, relaxation off)".to_string()
            } else {
                format!("{passes} pass(es)")
            };
            sweep_row(&tag, &s);
        }
    }
}

// ─────────────────────────────── PART F ────────────────────────────────
// THE WARP SWEEP WITH FRINGE LENGTH AS THE PRIMARY METRIC, and a PNG per setting.
//
// PRIMARY (optimised): FRINGE LENGTH — p50, p90 and MAX. It is the LONG spurs that are visible,
//   and a median can improve while the tail stays, so the tail is reported.
// REPORTED ONLY (never optimised): regularity CV, axis R, coastline length, spur count. A
//   shorter coast is not a better one and a more regular one is not either — those numbers are
//   here to catch a remedy that trades the symptom for a different defect.
//
// AND A PNG PER SETTING. Three wrong attributions and two remedies evaluated against metrics
// that did not measure the symptom; a picture per setting would have caught each sooner. The
// crop is chosen ONCE, on the shipped config, as the window with the most spur vertices, and
// then held FIXED for every panel — otherwise the panels are not comparable.

const CROP: usize = 640;

/// Land/sea render with the contour drawn on top: sea 0, land 0.65, coastline 1.0.
fn render_crop(
    field: &GridF32,
    polys: &[Polyline],
    ox: usize,
    oy: usize,
    path: &std::path::Path,
) -> Result<(), String> {
    let mut img = GridF32::new(CROP, CROP, 0.0);
    for y in 0..CROP {
        for x in 0..CROP {
            let (sx, sy) = (ox + x, oy + y);
            if sx < field.width && sy < field.height && field.data[sy * field.width + sx] > SEA {
                img.set(x, y, 0.65);
            }
        }
    }
    for pl in polys {
        for &(px, py) in pl.iter() {
            let (ix, iy) = (px.round() as isize - ox as isize, py.round() as isize - oy as isize);
            if ix >= 0 && iy >= 0 && (ix as usize) < CROP && (iy as usize) < CROP {
                img.set(ix as usize, iy as usize, 1.0);
            }
        }
    }
    img.save_png_u8(path)
}

/// The window holding the most spur vertices — chosen once, then held fixed.
fn densest_window(polys: &[Polyline], w: usize, h: usize, km_per_cell: f32) -> (usize, usize) {
    let (min_spur_km, neck_km) = (1.0f32, 0.6f32);
    let _ = spurs(polys, km_per_cell, min_spur_km, neck_km); // same parameters as the metric
    let step = CROP / 2;
    let (mut best, mut bxy) = (0usize, (0usize, 0usize));
    let mut oy = 0;
    while oy + CROP <= h {
        let mut ox = 0;
        while ox + CROP <= w {
            let mut n = 0usize;
            for pl in polys {
                for &(px, py) in pl.iter() {
                    let (x, y) = (px as usize, py as usize);
                    if x >= ox && x < ox + CROP && y >= oy && y < oy + CROP {
                        n += 1;
                    }
                }
            }
            if n > best {
                best = n;
                bxy = (ox, oy);
            }
            ox += step;
        }
        oy += step;
    }
    bxy
}

fn len_row(label: &str, s: &Spurs, km_per_cell: f32, max_km: f32) {
    eprintln!(
        "{label:<24} {:>8.2} {:>8.2} {:>8.2} | {:>7.1} {:>7.1} {:>7.1} | {:>7} {:>8.0} {:>6.2} {:>6.3}",
        s.med_len_km,
        s.p90_len_km,
        max_km,
        s.med_len_km / km_per_cell,
        s.p90_len_km / km_per_cell,
        max_km / km_per_cell,
        s.count,
        s.coast_km,
        s.spacing_cv,
        s.local_axis_r
    );
}

fn len_header() {
    eprintln!(
        "{:<24} {:>27} | {:>25} | {:>32}",
        "", "PRIMARY: fringe length km", "the same, in CELLS", "REPORTED ONLY (not optimised)"
    );
    eprintln!(
        "{:<24} {:>8} {:>8} {:>8} | {:>7} {:>7} {:>7} | {:>7} {:>8} {:>6} {:>6}",
        "coast_warp_strength",
        "p50",
        "p90",
        "max",
        "p50",
        "p90",
        "max",
        "spurs",
        "coast km",
        "CV",
        "LOCAL R"
    );
}

/// Longest spur, needed for the `max` column the median can hide.
fn max_spur_km(polys: &[Polyline], km_per_cell: f32, min_spur_km: f32, neck_km: f32) -> f32 {
    let neck = neck_km / km_per_cell;
    let min_arc = min_spur_km / km_per_cell;
    let mut best = 0.0f32;
    for pl in polys {
        let mut arc = vec![0.0f32; pl.len()];
        for i in 1..pl.len() {
            arc[i] = arc[i - 1]
                + ((pl[i].0 - pl[i - 1].0).powi(2) + (pl[i].1 - pl[i - 1].1).powi(2)).sqrt();
        }
        let mut i = 0usize;
        while i < pl.len() {
            let mut hit = None;
            let mut j = i + 1;
            while j < pl.len() && (arc[j] - arc[i]) * km_per_cell < MAX_SPUR_KM {
                if arc[j] - arc[i] >= min_arc {
                    let d = ((pl[j].0 - pl[i].0).powi(2) + (pl[j].1 - pl[i].1).powi(2)).sqrt();
                    if d <= neck {
                        hit = Some(j);
                    }
                }
                j += 1;
            }
            match hit {
                Some(j) => {
                    best = best.max((arc[j] - arc[i]) * km_per_cell);
                    i = j;
                }
                None => i += 1,
            }
        }
    }
    best
}

#[test]
#[ignore]
fn warp_sweep_fringe_length_with_renders() {
    let (min_spur_km, neck_km) = (1.0f32, 0.6f32);
    // A test's cwd is the CRATE root, not the workspace root — the first run wrote to
    // `crates/ymir-core/exports/` where nobody would look for it.
    let out = std::path::Path::new("../../exports/coastal_fringes");
    std::fs::create_dir_all(out).expect("output dir");
    eprintln!(
        "\n=====  PART F — warp sweep, FRINGE LENGTH primary, with renders  =====\n\n\
         PRIMARY (optimised): fringe length p50 / p90 / MAX. The tail is reported because a\n\
         median can improve while the long spurs the author sees remain.\n\
         REPORTED ONLY: spur count, coastline length, CV, axis R — a shorter or more regular\n\
         coast is NOT a better one; these are here to catch a remedy that trades one defect for\n\
         another.\n\
         PNGs: exports/coastal_fringes/ — sea black, land grey, coastline white; ONE fixed crop."
    );

    for target in [2048usize, 8192] {
        let km = DOMAIN_KM / target as f32;
        eprintln!("\n--- {target}² ({:.0} m/cell) ---", km * 1000.0);
        len_header();

        // The crop, chosen once on the shipped config and then held fixed.
        let shipped = hd_warp(target, 1.5);
        let shipped_polys = marching_squares(&shipped, SEA);
        let (ox, oy) = densest_window(&shipped_polys, target, target, km);
        eprintln!(
            "{:<24} crop {CROP}×{CROP} at ({ox},{oy}) — the densest spur window, FIXED for every panel",
            ""
        );

        for warp in [1.5f64, 1.0, 0.5, 0.25, 0.0] {
            let field =
                if (warp - 1.5).abs() < 1e-9 { shipped.clone() } else { hd_warp(target, warp) };
            let polys = marching_squares(&field, SEA);
            let s = spurs(&polys, km, min_spur_km, neck_km);
            let mx = max_spur_km(&polys, km, min_spur_km, neck_km);
            let tag = if (warp - 1.5).abs() < 1e-9 {
                "1.50 (SHIPPED)".to_string()
            } else {
                format!("{warp:.2}")
            };
            len_row(&tag, &s, km, mx);
            if target == 8192 {
                let p = out.join(format!("warp_{:.2}_8192.png", warp));
                render_crop(&field, &polys, ox, oy, &p).expect("render");
            }
        }

        // THE REFERENCE PANEL: the coarse field upscaled, no FBM, no incision — what a coast
        // with no manufactured spurs looks like at this resolution.
        let coarse = hd_stage(target, false, 0.0, 0);
        let cpolys = marching_squares(&coarse, SEA);
        let cs = spurs(&cpolys, km, min_spur_km, neck_km);
        let cmx = max_spur_km(&cpolys, km, min_spur_km, neck_km);
        len_row("REFERENCE coarse", &cs, km, cmx);
        if target == 8192 {
            render_crop(&coarse, &cpolys, ox, oy, &out.join("reference_coarse_8192.png"))
                .expect("render");
        }
    }
    eprintln!(
        "\n⇒ The question: does any warp bring the fringe LENGTH near the coarse reference\n  \
         WITHOUT erasing the indentation? Read the p90 and max columns against the REFERENCE row,\n  \
         then look at the panels."
    );
}
