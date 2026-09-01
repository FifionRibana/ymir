//! H-1c — APPLYING the water balance: endorheic basins settle at their evaporative
//! equilibrium instead of stagnating at their col. This measures the plain area gained BY
//! THE BALANCE ALONE, so that whatever H-2 (sill incision) adds later is attributable to
//! incision and to nothing else — the same discipline as keeping C-3's hard basement at ×1.
//!
//! Reported per climate, both resolutions:
//!   - lakes and water area BEFORE (all at sill) vs AFTER (endorheic at equilibrium);
//!   - the FLOOR EXPOSED (km²) and — the decisive question — ITS SLOPE DISTRIBUTION. An
//!     endorheic floor SHOULD be flat, but if erosion gullied it before the basin filled it
//!     will not be, and the "buildable plain" argument evaporates. Measured, not assumed.
//!   - the CONTIGUOUS PLAIN metric (200 m physical bridging, ×1.03 convergence — see
//!     `plain_metric.rs`) before and after.
//!
//! Run: cargo test -p ymir-core --test h1c_endorheic_shrink h1c -- --ignored --nocapture

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, apply_lake_water_balance, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::{FlowConfig, compute_flow};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;
const SPAN_DEG: f32 = 10.0;
const MAX_PLAIN_DEG: f32 = 5.0;
/// The validated grid-stable bridging distance (metres) — see ADR "contiguous plain".
const BRIDGE_M: f32 = 200.0;
const HEX_KM2: f32 = 2.598_076_2;

fn slope_deg_at(field: &GridF32, i: usize, j: usize, cell_m: f32) -> f32 {
    let (gx, gy) = field.gradient_at(i, j);
    ((gx * gx + gy * gy).sqrt() * FULL_RANGE_M / cell_m).atan().to_degrees()
}

fn chamfer_dt(mask: &[bool], w: usize, h: usize) -> Vec<f32> {
    const BIG: f32 = 1e9;
    const D2: f32 = 1.414_213_6;
    let mut d: Vec<f32> = mask.iter().map(|&m| if m { 0.0 } else { BIG }).collect();
    for j in 0..h {
        for i in 0..w {
            let k = j * w + i;
            let mut v = d[k];
            if j > 0 {
                v = v.min(d[k - w] + 1.0);
                if i > 0 {
                    v = v.min(d[k - w - 1] + D2);
                }
                if i + 1 < w {
                    v = v.min(d[k - w + 1] + D2);
                }
            }
            if i > 0 {
                v = v.min(d[k - 1] + 1.0);
            }
            d[k] = v;
        }
    }
    for j in (0..h).rev() {
        for i in (0..w).rev() {
            let k = j * w + i;
            let mut v = d[k];
            if j + 1 < h {
                v = v.min(d[k + w] + 1.0);
                if i + 1 < w {
                    v = v.min(d[k + w + 1] + D2);
                }
                if i > 0 {
                    v = v.min(d[k + w - 1] + D2);
                }
            }
            if i + 1 < w {
                v = v.min(d[k + 1] + 1.0);
            }
            d[k] = v;
        }
    }
    d
}

/// Largest contiguous plain (km², hex) using the validated definition: flat + dry land,
/// morphological closing at `BRIDGE_M`, then connected components.
fn plain(field: &GridF32, lake_map: &[u32], km_per_cell: f32) -> (f32, u64, f32) {
    let (w, h) = (field.width, field.height);
    let cell_m = km_per_cell * 1000.0;
    let cell_km2 = km_per_cell * km_per_cell;
    let mut ok = vec![false; w * h];
    for j in 1..h - 1 {
        for i in 1..w - 1 {
            let k = j * w + i;
            if field.data[k] <= 0.5 || lake_map[k] != 0 {
                continue;
            }
            if slope_deg_at(field, i, j, cell_m) < MAX_PLAIN_DEG {
                ok[k] = true;
            }
        }
    }
    let r = BRIDGE_M / cell_m;
    let closed = {
        let dt = chamfer_dt(&ok, w, h);
        let dil: Vec<bool> = dt.iter().map(|&v| v <= r).collect();
        let inv: Vec<bool> = dil.iter().map(|&m| !m).collect();
        let dt2 = chamfer_dt(&inv, w, h);
        dt2.iter().map(|&v| v > r).collect::<Vec<bool>>()
    };
    let n = w * h;
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn find(p: &mut Vec<u32>, mut x: u32) -> u32 {
        while p[x as usize] != x {
            let g = p[p[x as usize] as usize];
            p[x as usize] = g;
            x = g;
        }
        x
    }
    let mut uni = |p: &mut Vec<u32>, a: u32, b: u32| {
        let (ra, rb) = (find(p, a), find(p, b));
        if ra != rb {
            p[rb as usize] = ra;
        }
    };
    for j in 0..h {
        for i in 0..w {
            let k = j * w + i;
            if !closed[k] {
                continue;
            }
            let rr = j * w + (i + 1) % w;
            if closed[rr] {
                uni(&mut parent, k as u32, rr as u32);
            }
            let dd = ((j + 1) % h) * w + i;
            if closed[dd] {
                uni(&mut parent, k as u32, dd as u32);
            }
        }
    }
    let mut sizes = std::collections::HashMap::new();
    let mut total = 0u64;
    for k in 0..n {
        if closed[k] {
            total += 1;
            let rt = find(&mut parent, k as u32);
            *sizes.entry(rt).or_insert(0u64) += 1;
        }
    }
    let largest = sizes.values().copied().max().unwrap_or(0);
    (
        largest as f32 * cell_km2,
        ((largest as f32 * cell_km2) / HEX_KM2) as u64,
        total as f32 * cell_km2,
    )
}

fn run(target: usize) {
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
    let km_per_cell = DOMAIN_KM / target as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let cell_m = km_per_cell * 1000.0;
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
    cfg.sample_origin = [0.09375, 0.578125];
    cfg.sample_size = 1.0;
    cfg.erosion = None;
    cfg.lithology = LithologyConfig {
        enabled: true,
        soft_multiplier: 10.0,
        volcanic_multiplier: 3.0,
        rift_age_threshold: 1.0,
    };
    cfg.fracture = FractureConfig {
        enabled: true,
        amplitude: 6.0,
        decay_km: 25.0,
        domain_km: DOMAIN_KM,
        ..Default::default()
    };
    cfg.stream_power = {
        let mut sp = ymir_core::erosion::stream_power::StreamPowerConfig::relief_v3(
            cell_km2,
            ss.depth_scale_m as f32,
        );
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;
        Some(sp)
    };
    let (up, _c) = upscale_from_c1_with_progress(
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
    let raw = up.heightmap;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field = ymir_core::terrain::flow::breach_monotone(
        &raw,
        &pre.flow.filled,
        &pre.lake_map,
        0.5,
        raw.width,
        raw.height,
    );
    let (w, h) = (field.width, field.height);
    let flow =
        compute_flow(&field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });

    let (p0_lg, p0_hex, p0_tot) = plain(&field, &pre.lake_map, km_per_cell);
    let w0: f32 = pre.lakes.iter().map(|l| l.area_km2).sum();
    eprintln!("\n================  H-1c ENDORHEIC EQUILIBRIUM  {target}²  ================");
    eprintln!(
        "BEFORE (all lakes at their sill): {} lakes = {:.0} km² water | plain largest {:.0} km² ({} hex), total {:.0} km²",
        pre.lakes.len(),
        w0,
        p0_lg,
        p0_hex,
        p0_tot
    );
    eprintln!(
        "\n{:>11} | {:>13} | {:>10} | {:>12} | {:>26} | {:>22}",
        "climate",
        "lakes after",
        "water km²",
        "floor km²",
        "exposed-floor slope (°)",
        "plain largest / total"
    );
    for (label, lat) in
        [("tropical", 10.0f32), ("arid-hot", 25.0), ("humid", 45.0), ("arid-cold", 65.0)]
    {
        let climate =
            c1_climate_placed(&field, &ss, lat, SPAN_DEG, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let mut lm = pre.lake_map.clone();
        let after = apply_lake_water_balance(
            &field, &flow, &dclim, cell_km2, &ss, &pre.lakes, &mut lm, None, w, h,
        );
        let w1: f32 = after.iter().map(|l| l.area_km2).sum();
        // Newly exposed floor: was lake, now dry.
        let mut slopes = Vec::new();
        let mut exposed = 0u64;
        for j in 1..h - 1 {
            for i in 1..w - 1 {
                let k = j * w + i;
                if pre.lake_map[k] != 0 && lm[k] == 0 {
                    exposed += 1;
                    slopes.push(slope_deg_at(&field, i, j, cell_m));
                }
            }
        }
        slopes.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |p: f32| {
            if slopes.is_empty() {
                0.0
            } else {
                slopes[((slopes.len() as f32 * p) as usize).min(slopes.len() - 1)]
            }
        };
        let under5 = 100.0 * slopes.iter().filter(|&&s| s < MAX_PLAIN_DEG).count() as f32
            / slopes.len().max(1) as f32;
        let (p1_lg, p1_hex, p1_tot) = plain(&field, &lm, km_per_cell);
        let endo = after.iter().filter(|l| l.lake_type == LakeType::Endorheic).count();
        eprintln!(
            "{label:>11} | {:>5} ({endo:>2} endo) | {w1:>10.0} | {:>12.0} | p50 {:>4.1} p90 {:>4.1} <5° {:>3.0}% | {:>8.0} / {:>7.0} ({} hex)",
            after.len(),
            exposed as f32 * cell_km2,
            q(0.5),
            q(0.9),
            under5,
            p1_lg,
            p1_tot,
            p1_hex
        );
    }
    eprintln!(
        "(plain BEFORE was {:.0} km² largest / {:.0} total. The gain here is due to the BALANCE ALONE —\n H-2's sill incision has not run, so whatever it adds later is attributable to incision only.)",
        p0_lg, p0_tot
    );
}

#[test]
#[ignore]
fn h1c() {
    run(2048);
}

#[test]
#[ignore]
fn h1c_8192() {
    run(8192);
}
