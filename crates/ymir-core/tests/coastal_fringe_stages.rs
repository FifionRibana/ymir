//! ADR Finding 51 — STAGE BISECTION on the fringe metric. Where does the 82 % baseline enter?
//!
//! The same instrument that attributed 90 682 depressions to the FBM in one table: measure the
//! metric after EACH stage of the build and read where it jumps. Attribution before remedy.
//!
//!   1. coarse post-isostasy   `amplitude_base = 0`, no incision — the bilinear C-1 altitude
//!   2. + FBM upscale          production amplitude, still no incision
//!   3. + relief-v2 incision   the pre-v3 closure set (nonlinear hillslope + lateral widening)
//!   4. + relief-v3 incision   the shipped incision, no K closures
//!   5. final                  + C-3 lithology + C-3b fracture (the shipped product)
//!
//! Two interpretation traps, both watched here:
//!   - the coastline is extracted AT SEA LEVEL, so a stage can raise the fringe rate WITHOUT
//!     touching the coast, merely by FLATTENING terrain near zero. The near-zero slope
//!     distribution is reported beside the rate; if they move together the mechanism is the
//!     "near-flat terrain at sea level" verdict again.
//!   - if the COARSE field already fringes, the subject stops being erosion at all and moves to
//!     the tectonic stage.
//!
//! Run: cargo test -p ymir-core --release --test coastal_fringe_stages -- --ignored --nocapture

use ymir_core::erosion::stream_power::StreamPowerConfig;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
/// Band around sea level in metres within which "near-zero" slopes are measured — the terrain
/// the coastline is actually traced through.
const NEAR_SEA_BAND_M: f32 = 20.0;

#[derive(Clone, Copy, PartialEq)]
enum Stage {
    Coarse,
    Fbm,
    ReliefV2,
    ReliefV3,
    Final,
}

impl Stage {
    fn label(self) -> &'static str {
        match self {
            Stage::Coarse => "1 coarse (no FBM, no incision)",
            Stage::Fbm => "2 + FBM upscale",
            Stage::ReliefV2 => "3 + relief-v2 incision",
            Stage::ReliefV3 => "4 + relief-v3 incision",
            Stage::Final => "5 final (+ C-3 + C-3b)",
        }
    }
}

fn build(target: usize, st: Stage) -> GridF32 {
    let ss = SteinSteinParams::default();
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
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let on = st == Stage::Final;
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: if st == Stage::Coarse { 0.0 } else { 0.04 },
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: on,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: on,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    let cell_km2 = (DOMAIN_KM / target as f32).powi(2);
    match st {
        Stage::Coarse | Stage::Fbm => cfg.stream_power = None,
        Stage::ReliefV2 => {
            cfg.stream_power = Some(StreamPowerConfig::relief_v2(cell_km2, ss.depth_scale_m as f32))
        }
        Stage::ReliefV3 | Stage::Final => {} // production_hd_config already sets relief-v3
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

fn sample(g: &GridF32, x: f32, y: f32) -> f32 {
    let (w, h) = (g.width as i32, g.height as i32);
    let x = x.clamp(0.0, (w - 1) as f32);
    let y = y.clamp(0.0, (h - 1) as f32);
    let (x0, y0) = (x.floor() as i32, y.floor() as i32);
    let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
    let (fx, fy) = (x - x0 as f32, y - y0 as f32);
    let v = |a: i32, b: i32| g.data[b as usize * g.width + a as usize];
    let t = v(x0, y0) * (1.0 - fx) + v(x1, y0) * fx;
    let b = v(x0, y1) * (1.0 - fx) + v(x1, y1) * fx;
    t * (1.0 - fy) + b * fy
}

fn slope_deg(g: &GridF32, x: f32, y: f32, m_per_cell: f32, ss: &SteinSteinParams) -> f32 {
    let to_m = |v: f32| c1_altitude_norm_to_metres(v, ss);
    let gx = 0.5 * (to_m(sample(g, x + 1.0, y)) - to_m(sample(g, x - 1.0, y)));
    let gy = 0.5 * (to_m(sample(g, x, y + 1.0)) - to_m(sample(g, x, y - 1.0)));
    ((gx * gx + gy * gy).sqrt() / m_per_cell).atan().to_degrees()
}

const EDGES: [f32; 4] = [0.5, 2.0, 5.0, 15.0];
const NAMES: [&str; 5] = ["<0.5", "0.5-2", "2-5", "5-15", ">15"];
fn class_of(d: f32) -> usize {
    EDGES.iter().position(|&e| d < e).unwrap_or(4)
}

struct Row {
    rings: usize,
    tiny: usize,
    turns: usize,
    fringes: usize,
    len_km: f32,
    per_class: [(usize, usize); 5],
    /// Near-sea-level land cells and the share of them under 0.5° / 2°.
    near_n: usize,
    near_flat05: f32,
    near_flat2: f32,
    near_median_slope: f32,
}

fn measure(g: &GridF32, m_per_cell: f32, ss: &SteinSteinParams) -> Row {
    let polys = marching_squares(g, SEA);
    let (mut turns, mut fringes, mut tiny) = (0usize, 0usize, 0usize);
    let mut len = 0.0f32;
    let mut per_class = [(0usize, 0usize); 5];
    for pl in &polys {
        if pl.len() < 20 {
            tiny += 1;
        }
        for w2 in pl.windows(2) {
            len += ((w2[1].0 - w2[0].0).powi(2) + (w2[1].1 - w2[0].1).powi(2)).sqrt();
        }
        for w3 in pl.windows(3) {
            let (ax, ay) = (w3[1].0 - w3[0].0, w3[1].1 - w3[0].1);
            let (bx, by) = (w3[2].0 - w3[1].0, w3[2].1 - w3[1].1);
            let (la, lb) = ((ax * ax + ay * ay).sqrt(), (bx * bx + by * by).sqrt());
            if la <= 0.0 || lb <= 0.0 {
                continue;
            }
            let t = (((ax * bx + ay * by) / (la * lb)).clamp(-1.0, 1.0)).acos().to_degrees();
            turns += 1;
            let k = class_of(slope_deg(g, w3[1].0, w3[1].1, m_per_cell, ss));
            per_class[k].0 += 1;
            if t > 80.0 {
                fringes += 1;
                per_class[k].1 += 1;
            }
        }
    }
    // Near-zero slope distribution: LAND cells within NEAR_SEA_BAND_M of sea level. This is the
    // terrain the contour is traced through, and a stage can raise the fringe rate purely by
    // flattening it.
    let (w, h) = (g.width, g.height);
    let mut slopes: Vec<f32> = Vec::new();
    let step = if w > 4096 { 4 } else { 2 }; // subsample: the distribution, not every cell
    for y in (1..h - 1).step_by(step) {
        for x in (1..w - 1).step_by(step) {
            let v = g.data[y * w + x];
            if v <= SEA {
                continue;
            }
            if c1_altitude_norm_to_metres(v, ss) > NEAR_SEA_BAND_M {
                continue;
            }
            slopes.push(slope_deg(g, x as f32, y as f32, m_per_cell, ss));
        }
    }
    slopes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = slopes.len();
    let f = |t: f32| 100.0 * slopes.iter().filter(|&&s| s < t).count() as f32 / n.max(1) as f32;
    Row {
        rings: polys.len(),
        tiny,
        turns,
        fringes,
        len_km: len * m_per_cell / 1000.0,
        per_class,
        near_n: n,
        near_flat05: f(0.5),
        near_flat2: f(2.0),
        near_median_slope: if n == 0 { 0.0 } else { slopes[n / 2] },
    }
}

fn run(target: usize) {
    let ss = SteinSteinParams::default();
    let m_per_cell = DOMAIN_KM * 1000.0 / target as f32;
    eprintln!(
        "\n==========  FRINGE STAGE BISECTION — {target}^2, {m_per_cell:.0} m/cell  =========="
    );
    eprintln!(
        "\n{:<32} {:>7} {:>7} {:>9} {:>8} {:>9}",
        "stage", "rings", "tiny", "fringes", "rate %", "coast km"
    );
    let stages = [Stage::Coarse, Stage::Fbm, Stage::ReliefV2, Stage::ReliefV3, Stage::Final];
    let rows: Vec<Row> =
        stages.iter().map(|&s| measure(&build(target, s), m_per_cell, &ss)).collect();
    let mut prev: Option<&Row> = None;
    for (i, r) in rows.iter().enumerate() {
        let rate = 100.0 * r.fringes as f32 / r.turns.max(1) as f32;
        eprintln!(
            "{:<32} {:>7} {:>7} {:>9} {:>7.2}% {:>9.0}",
            stages[i].label(),
            r.rings,
            r.tiny,
            r.fringes,
            rate,
            r.len_km
        );
        if let Some(p) = prev {
            let prate = 100.0 * p.fringes as f32 / p.turns.max(1) as f32;
            eprintln!(
                "{:<32} {:>+7} {:>+7} {:>+9} {:>+7.2}  {:>+9.0}   <- JUMP",
                "",
                r.rings as i64 - p.rings as i64,
                r.tiny as i64 - p.tiny as i64,
                r.fringes as i64 - p.fringes as i64,
                rate - prate,
                r.len_km - p.len_km
            );
        }
        prev = Some(r);
    }

    eprintln!("\nPER SLOPE CLASS — turns > 80 deg (the fringes live on the LOW-slope shores)");
    eprint!("{:<32}", "stage");
    for n in NAMES {
        eprint!("{n:>10}");
    }
    eprintln!();
    for (i, r) in rows.iter().enumerate() {
        eprint!("{:<32}", stages[i].label());
        for c in 0..5 {
            let (n, k) = r.per_class[c];
            if n == 0 {
                eprint!("{:>10}", "-");
            } else {
                eprint!("{:>9.1}%", 100.0 * k as f32 / n as f32);
            }
        }
        eprintln!();
    }

    eprintln!(
        "\nNEAR-ZERO SLOPE DISTRIBUTION — land within {NEAR_SEA_BAND_M:.0} m of sea level.\n\
         If this moves WITH the fringe rate, the stage is flattening the coastal terrain rather\n\
         than crenellating the coast, and it is the same mechanism as the 'near-flat at sea\n\
         level' verdict."
    );
    eprintln!(
        "{:<32} {:>10} {:>11} {:>11} {:>13}",
        "stage", "cells", "< 0.5 deg", "< 2 deg", "median slope"
    );
    for (i, r) in rows.iter().enumerate() {
        eprintln!(
            "{:<32} {:>10} {:>10.1}% {:>10.1}% {:>12.2}°",
            stages[i].label(),
            r.near_n,
            r.near_flat05,
            r.near_flat2,
            r.near_median_slope
        );
    }
}

#[test]
#[ignore]
fn fringe_stages_2048() {
    run(2048);
}

#[test]
#[ignore]
fn fringe_stages_8192() {
    run(8192);
}
