//! ADR Finding 56 — the invariant suite for the slope-dependent channel head.
//!
//! `A_c` governs the hillslope/channel split EVERYWHERE, so the law changes the network far
//! beyond the coast. A coastal gain paid for by a degraded network would be a bad trade and the
//! coastline metrics cannot see it. This runs the WHOLE production hydro chain — conditioning,
//! climate, drainage, lakes — with the law OFF and ON, at both resolutions, and reports what
//! moved.
//!
//! Climate: the arid-hot test bed (25°, span 10) — the standing config for below-sea basins.
//! Note the TERRAIN is climate-independent (`upscale_from_c1_with_progress` takes no climate),
//! so the coastline results elsewhere are the same at any latitude; climate matters HERE,
//! because the lakes and the discharge depend on it.
//!
//! Run: cargo test -p ymir-core --release --test channel_head_law_invariants -- --ignored --nocapture

use std::collections::{HashMap, HashSet};

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::{CHANNEL_HEAD_S_MIN, CHANNEL_HEAD_S_REF, ChannelHeadLaw};
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, LakeType, Navigability,
    c1_drainage_windowed, exorheic_lakes_missing_outlet,
};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::breach_monotone;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const LAT: f32 = 25.0;
const SPAN: f32 = 10.0;
const SEA: f32 = 0.5;

fn terrain(target: usize, law: bool) -> GridF32 {
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
    if law {
        if let Some(sp) = cfg.stream_power.as_mut() {
            sp.a_c_slope_law =
                Some(ChannelHeadLaw { s_ref: CHANNEL_HEAD_S_REF, s_min: CHANNEL_HEAD_S_MIN });
        }
    }
    upscale_from_c1_with_progress(
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
    )
    .0
    .heightmap
}

/// The production hydro chain: pre-breach drainage → conditioning → climate → final drainage.
fn hydro(raw: &GridF32, ss: &SteinSteinParams) -> (GridF32, C1DrainageResult) {
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre = c1_drainage_windowed(raw, None, &dcfg, ss, DOMAIN_KM);
    let field = breach_monotone(raw, &pre.flow.filled, &pre.lake_map, SEA, raw.width, raw.height);
    let climate = c1_climate_placed(&field, ss, LAT, SPAN, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let dr = c1_drainage_windowed(&field, Some(&dclim), &dcfg, ss, DOMAIN_KM);
    (field, dr)
}

fn pct(v: &mut Vec<f32>, f: f32) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        return 0.0;
    }
    v[(((v.len() - 1) as f64 * f as f64).round() as usize).min(v.len() - 1)]
}

fn report(
    label: &str,
    target: usize,
    field: &GridF32,
    dr: &C1DrainageResult,
    ss: &SteinSteinParams,
) {
    let (w, h) = (field.width, field.height);
    let km = DOMAIN_KM / target as f32;
    let m_per_cell = km * 1000.0;
    let cell_km2 = km * km;
    let to_m = |v: f32| c1_altitude_norm_to_metres(v, ss);
    eprintln!("\n──────── {label} ────────");

    // ── LAKE INVARIANTS, over the WHOLE population ──────────────────────────────────
    let mut v_footprint = 0usize; // a lake cell above its own level
    let mut v_depth = 0usize; // depth != level − floor
    let mut cells_by_lake: HashMap<u32, usize> = HashMap::new();
    for (k, &id) in dr.lake_map.iter().enumerate() {
        if id != 0 {
            *cells_by_lake.entry(id).or_insert(0) += 1;
        }
        let _ = k;
    }
    for lk in &dr.lakes {
        let lvl = lk.level_m;
        let mut floor = f32::INFINITY;
        for (k, &id) in dr.lake_map.iter().enumerate() {
            if id == lk.base.id {
                let a = to_m(field.data[k]);
                floor = floor.min(a);
                if a > lvl + 0.5 {
                    v_footprint += 1;
                }
            }
        }
        if floor.is_finite() && ((lvl - floor) - lk.depth_m).abs() > 1.0 {
            v_depth += 1;
        }
    }
    let missing_outlet = exorheic_lakes_missing_outlet(dr);
    // no mouth without membership: a segment that starts at a lake shore must name a live id
    let live: HashSet<u32> = dr.lakes.iter().map(|l| l.base.id).collect();
    let mut orphan_mouths = 0usize;
    for s in &dr.rivers.segments {
        let (sx, sy) = s.points[0];
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (sx as i32 + dx, sy as i32 + dy);
                if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                    let id = dr.lake_map[ny as usize * w + nx as usize];
                    if id != 0 && !live.contains(&id) {
                        orphan_mouths += 1;
                    }
                }
            }
        }
    }
    // monotone long profiles, flat lake crossings tolerated
    let mut non_monotone = 0usize;
    for prof in &dr.segment_profile_m {
        if prof.windows(2).any(|p| p[1] > p[0] + 1.0) {
            non_monotone += 1;
        }
    }
    let (exo, endo) = (
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Exorheic).count(),
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).count(),
    );
    eprintln!(
        "LAKES  {} total ({exo} exo, {endo} endo) | water {:.0} km²\n  \
         INVARIANTS  footprint-above-level {v_footprint} | depth≠level−floor {v_depth} | \
         exorheic-without-outlet {} | orphan mouths {orphan_mouths} | non-monotone profiles \
         {non_monotone} of {}",
        dr.lakes.len(),
        cells_by_lake.values().sum::<usize>() as f32 * cell_km2,
        missing_outlet.len(),
        dr.rivers.segments.len()
    );

    // ── SHAPE METRICS ───────────────────────────────────────────────────────────────
    let win = (1000.0 / m_per_cell).round().max(1.0) as usize; // PHYSICAL 1 km window
    let (mut fr, mut slopes) = (Vec::new(), Vec::new());
    let step = if w > 4096 { 4 } else { 2 };
    for y in (win..h - win).step_by(step) {
        for x in (win..w - win).step_by(step) {
            let k = y * w + x;
            if field.data[k] <= SEA {
                continue;
            }
            let gx = 0.5 * (to_m(field.data[k + 1]) - to_m(field.data[k - 1]));
            let gy = 0.5 * (to_m(field.data[k + w]) - to_m(field.data[k - w]));
            slopes.push(((gx * gx + gy * gy).sqrt() / m_per_cell).atan().to_degrees());
            // floor / local ridge over the 1 km window
            let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
            for j in (y - win)..=(y + win) {
                for i in (x - win)..=(x + win) {
                    let a = to_m(field.data[j * w + i]);
                    lo = lo.min(a);
                    hi = hi.max(a);
                }
            }
            if hi > lo + 1.0 {
                fr.push((to_m(field.data[k]) - lo) / (hi - lo));
            }
        }
    }
    let n_s = slopes.len().max(1);
    let gt30 = 100.0 * slopes.iter().filter(|&&s| s > 30.0).count() as f32 / n_s as f32;
    let gt45 = 100.0 * slopes.iter().filter(|&&s| s > 45.0).count() as f32 / n_s as f32;
    let land_cells = field.data.iter().filter(|&&v| v > SEA).count();
    let chan_cells: usize = dr.rivers.segments.iter().map(|s| s.points.len()).sum();
    eprintln!(
        "SHAPE  floor/local-ridge (1 km window) p25 {:.3} | p50 {:.3} | p75 {:.3}\n  \
         slope >30° {gt30:.2} % | >45° {gt45:.2} % | median {:.2}°\n  \
         drainage density {:.3} km/km² ({} channel cells over {:.0} km² of land)",
        pct(&mut fr.clone(), 0.25),
        pct(&mut fr.clone(), 0.50),
        pct(&mut fr.clone(), 0.75),
        pct(&mut slopes.clone(), 0.50),
        chan_cells as f32 * km / (land_cells as f32 * cell_km2).max(1e-6),
        chan_cells,
        land_cells as f32 * cell_km2
    );

    // ── NETWORK: Strahler histogram, W/D per order, confluences, navigability ───────
    let mut by_order: HashMap<u8, Vec<usize>> = HashMap::new();
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        by_order.entry(s.strahler_order).or_default().push(i);
    }
    let confluences = dr.rivers.segments.iter().filter(|s| s.upstream.len() >= 2).count();
    let mut orders: Vec<u8> = by_order.keys().copied().collect();
    orders.sort_unstable();
    eprint!("NETWORK  confluences {confluences} | Strahler");
    for o in &orders {
        eprint!("  S{o}:{}", by_order[o].len());
    }
    eprintln!();
    eprint!("  W/D per order (width_m / depth-to-local-ridge over 1 km):");
    for o in &orders {
        let mut wd: Vec<f32> = Vec::new();
        for &i in &by_order[o] {
            let wm = dr.segment_width_m[i];
            let (px, py) = dr.rivers.segments[i].points[dr.rivers.segments[i].points.len() / 2];
            let (x, y) = (px as usize, py as usize);
            if x < win || y < win || x + win >= w || y + win >= h {
                continue;
            }
            let mut hi = f32::NEG_INFINITY;
            for j in (y - win)..=(y + win) {
                for i2 in (x - win)..=(x + win) {
                    hi = hi.max(to_m(field.data[j * w + i2]));
                }
            }
            let d = hi - to_m(field.data[y * w + x]);
            if d > 1.0 && wm > 0.0 {
                wd.push(wm / d);
            }
        }
        if !wd.is_empty() {
            eprint!("  S{o}:{:.1}", pct(&mut wd, 0.5));
        }
    }
    eprintln!();
    let mut nav: HashMap<&str, usize> = HashMap::new();
    for n in &dr.segment_navigability {
        *nav.entry(match n {
            Navigability::NonNavigable => "non-nav",
            Navigability::SmallBoat => "boat",
            Navigability::Barge => "barge",
            Navigability::Ship => "ship",
        })
        .or_insert(0) += 1;
    }
    eprintln!(
        "  NAVIGABILITY  non-nav {} | boat {} | barge {} | ship {}",
        nav.get("non-nav").copied().unwrap_or(0),
        nav.get("boat").copied().unwrap_or(0),
        nav.get("barge").copied().unwrap_or(0),
        nav.get("ship").copied().unwrap_or(0)
    );

    // ── CLOSED DEPRESSIONS — C-1's conditioning target ─────────────────────────────
    let pits =
        dr.flow.filled.data.iter().zip(field.data.iter()).filter(|(f, r)| **f > **r + 1e-6).count();
    eprintln!("DEPRESSIONS  {pits} filled cells ({:.0} km²)", pits as f32 * cell_km2);
}

#[test]
#[ignore]
fn channel_head_law_invariants() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  CHANNEL-HEAD LAW — INVARIANT SUITE (arid-hot {LAT}°/span {SPAN})  =========="
    );
    for target in [2048usize, 8192] {
        eprintln!("\n╔══════ {target}² ══════╗");
        for law in [false, true] {
            let raw = terrain(target, law);
            let (field, dr) = hydro(&raw, &ss);
            report(if law { "LAW ON" } else { "SHIPPED (law off)" }, target, &field, &dr, &ss);
        }
    }
}
