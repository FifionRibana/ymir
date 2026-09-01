//! GRID-STABLE contiguous plain metric (pre-H-2). The naive metric — connected components
//! of land under 5° — is NOT a property of the continent but of the grid step: measured
//! 1074 km² (largest piece) at 2048² against 216 km² at 8192², for a comparable total flat
//! area (9891 vs 7493 km²). Fine relief resolves and fragments what the coarse grid
//! smoothed over. H-2 will be judged on the plain area it gains, so the measure must
//! converge BEFORE H-2 starts.
//!
//! Two candidate definitions, implemented and COMPARED (not argued):
//!   a) TOLERANT TO MICRO-OBSTACLES — a morphological CLOSING over a stated PHYSICAL
//!      distance (metres, converted to cells at each resolution) bridges passable cuts
//!      (one builds either side of a ten-metre gully) before the component pass;
//!   b) AT THE GAME'S RESOLUTION — the flat mask is resampled to the hex pitch (the
//!      consumer's grain) before the component pass.
//!
//! THE TEST: agreement between 2048² and 8192². Whichever converges wins; if neither does,
//! that is the finding.
//!
//! Run: cargo test -p ymir-core --test plain_metric plain_metric -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;
const MAX_PLAIN_DEG: f32 = 5.0;
/// Living Landz hex: 1 km EDGE → pitch (centre-to-centre) √3 km, area (3√3/2) km².
const HEX_PITCH_KM: f32 = 1.732_050_8;
const HEX_KM2: f32 = 2.598_076_2;

/// The raw flat-land predicate: land, dry (not a lake), slope < `max_deg`.
fn flat_mask(field: &GridF32, lake_map: &[u32], km_per_cell: f32, max_deg: f32) -> Vec<bool> {
    let (w, h) = (field.width, field.height);
    let cell_m = km_per_cell * 1000.0;
    let mut ok = vec![false; w * h];
    for j in 1..h - 1 {
        for i in 1..w - 1 {
            let k = j * w + i;
            if field.data[k] <= 0.5 || lake_map[k] != 0 {
                continue;
            }
            let (gx, gy) = field.gradient_at(i, j);
            let slope_m = (gx * gx + gy * gy).sqrt() * FULL_RANGE_M / cell_m;
            if slope_m.atan().to_degrees() < max_deg {
                ok[k] = true;
            }
        }
    }
    ok
}

/// Two-pass chamfer distance (in CELLS) to the nearest `true` cell of `mask`.
fn chamfer_dt(mask: &[bool], w: usize, h: usize) -> Vec<f32> {
    const BIG: f32 = 1e9;
    const D1: f32 = 1.0;
    const D2: f32 = 1.414_213_6;
    let mut d: Vec<f32> = mask.iter().map(|&m| if m { 0.0 } else { BIG }).collect();
    // forward
    for j in 0..h {
        for i in 0..w {
            let k = j * w + i;
            let mut v = d[k];
            if j > 0 {
                v = v.min(d[k - w] + D1);
                if i > 0 {
                    v = v.min(d[k - w - 1] + D2);
                }
                if i + 1 < w {
                    v = v.min(d[k - w + 1] + D2);
                }
            }
            if i > 0 {
                v = v.min(d[k - 1] + D1);
            }
            d[k] = v;
        }
    }
    // backward
    for j in (0..h).rev() {
        for i in (0..w).rev() {
            let k = j * w + i;
            let mut v = d[k];
            if j + 1 < h {
                v = v.min(d[k + w] + D1);
                if i + 1 < w {
                    v = v.min(d[k + w + 1] + D2);
                }
                if i > 0 {
                    v = v.min(d[k + w - 1] + D2);
                }
            }
            if i + 1 < w {
                v = v.min(d[k + 1] + D1);
            }
            d[k] = v;
        }
    }
    d
}

/// Morphological CLOSING of `mask` by a radius of `r_cells` (dilate then erode) — bridges
/// cuts narrower than ~2·r without growing the outer boundary.
fn closing(mask: &[bool], w: usize, h: usize, r_cells: f32) -> Vec<bool> {
    if r_cells <= 0.0 {
        return mask.to_vec();
    }
    let dt = chamfer_dt(mask, w, h);
    let dilated: Vec<bool> = dt.iter().map(|&v| v <= r_cells).collect();
    let inv: Vec<bool> = dilated.iter().map(|&m| !m).collect();
    let dt2 = chamfer_dt(&inv, w, h);
    dt2.iter().map(|&v| v > r_cells).collect()
}

/// Connected components (periodic 4-connectivity, the `land_topology` union-find). Returns
/// (largest_cells, component_sizes_desc, total_cells).
fn components(mask: &[bool], w: usize, h: usize) -> (u64, Vec<u64>, u64) {
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
    let mut union = |p: &mut Vec<u32>, a: u32, b: u32| {
        let (ra, rb) = (find(p, a), find(p, b));
        if ra != rb {
            p[rb as usize] = ra;
        }
    };
    for j in 0..h {
        for i in 0..w {
            let k = j * w + i;
            if !mask[k] {
                continue;
            }
            let r = j * w + (i + 1) % w;
            if mask[r] {
                union(&mut parent, k as u32, r as u32);
            }
            let d = ((j + 1) % h) * w + i;
            if mask[d] {
                union(&mut parent, k as u32, d as u32);
            }
        }
    }
    let mut sizes = std::collections::HashMap::new();
    let mut total = 0u64;
    for k in 0..n {
        if mask[k] {
            total += 1;
            let r = find(&mut parent, k as u32);
            *sizes.entry(r).or_insert(0u64) += 1;
        }
    }
    let mut v: Vec<u64> = sizes.into_values().collect();
    v.sort_unstable_by(|a, b| b.cmp(a));
    (v.first().copied().unwrap_or(0), v, total)
}

/// Definition (b): resample the flat mask to the HEX PITCH (the consumer's grain) — a hex
/// cell counts as plain when the MAJORITY of the fine cells it covers are flat — then run
/// the component pass on that coarse grid. Both resolutions land on the same pitch.
fn hex_resample(mask: &[bool], w: usize, h: usize, km_per_cell: f32) -> (Vec<bool>, usize, usize) {
    let step = (HEX_PITCH_KM / km_per_cell).round().max(1.0) as usize;
    let (cw, ch) = (w / step, h / step);
    let mut out = vec![false; cw * ch];
    for cj in 0..ch {
        for ci in 0..cw {
            let (mut hit, mut tot) = (0u32, 0u32);
            for dj in 0..step {
                for di in 0..step {
                    let k = (cj * step + dj) * w + (ci * step + di);
                    tot += 1;
                    if mask[k] {
                        hit += 1;
                    }
                }
            }
            out[cj * cw + ci] = hit * 2 > tot;
        }
    }
    (out, cw, ch)
}

fn build_terrain(target: usize) -> (GridF32, Vec<u32>) {
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
    (field, pre.lake_map)
}

#[test]
#[ignore]
fn plain_metric() {
    eprintln!("\n=== CONTIGUOUS PLAIN — grid stability of two definitions ===");
    eprintln!(
        "(land <{MAX_PLAIN_DEG}°, dry; hex = 1 km edge → pitch {HEX_PITCH_KM:.2} km, {HEX_KM2:.2} km²)"
    );
    // Bridging distances in METRES (physical, converted per resolution).
    let bridges_m = [0.0f32, 100.0, 200.0, 400.0, 800.0];
    let mut naive = [0.0f32; 2];
    let mut a_res: [Vec<f32>; 2] = [vec![], vec![]];
    let mut b_res = [0.0f32; 2];
    let mut totals = [0.0f32; 2];

    for (ri, &target) in [2048usize, 8192].iter().enumerate() {
        let km_per_cell = DOMAIN_KM / target as f32;
        let cell_km2 = km_per_cell * km_per_cell;
        let (field, lake_map) = build_terrain(target);
        let mask = flat_mask(&field, &lake_map, km_per_cell, MAX_PLAIN_DEG);
        let (w, h) = (field.width, field.height);
        let (largest, sizes, total) = components(&mask, w, h);
        naive[ri] = largest as f32 * cell_km2;
        totals[ri] = total as f32 * cell_km2;
        eprintln!(
            "\n--- {target}²  (cell {:.0} m) ---\nNAIVE: largest {:.0} km² ({} hex) | total flat {:.0} km² | components {} | top5 {:?} km²",
            km_per_cell * 1000.0,
            naive[ri],
            (naive[ri] / HEX_KM2) as u64,
            totals[ri],
            sizes.len(),
            sizes.iter().take(5).map(|&s| (s as f32 * cell_km2) as u32).collect::<Vec<_>>()
        );
        // (a) closing at physical bridging distances
        for &bm in &bridges_m {
            let r_cells = bm / (km_per_cell * 1000.0);
            let m2 = closing(&mask, w, h, r_cells);
            let (lg, s, tt) = components(&m2, w, h);
            let km2 = lg as f32 * cell_km2;
            a_res[ri].push(km2);
            eprintln!(
                "  (a) bridge {bm:>5.0} m (r={r_cells:.1} cells): largest {km2:>8.0} km² = {:>6} hex | total {:.0} km² | pieces {}",
                (km2 / HEX_KM2) as u64,
                tt as f32 * cell_km2,
                s.len()
            );
            // Size DISTRIBUTION for the candidate that converges (200 m).
            if (bm - 200.0).abs() < 1.0 {
                let big: Vec<u32> =
                    s.iter().take(10).map(|&c| (c as f32 * cell_km2) as u32).collect();
                let over = |t: f32| s.iter().filter(|&&c| c as f32 * cell_km2 >= t).count();
                eprintln!(
                    "      top10 km²: {big:?}\n      pieces ≥1000 km²: {} | ≥500: {} | ≥100: {} | ≥10: {}",
                    over(1000.0),
                    over(500.0),
                    over(100.0),
                    over(10.0)
                );
            }
        }
        // (b) hex-pitch resample
        let (hm, cw, ch) = hex_resample(&mask, w, h, km_per_cell);
        let (lg, _s, tt) = components(&hm, cw, ch);
        let hex_cell_km2 = (HEX_PITCH_KM / km_per_cell).round()
            * km_per_cell
            * ((HEX_PITCH_KM / km_per_cell).round() * km_per_cell);
        b_res[ri] = lg as f32 * hex_cell_km2;
        eprintln!(
            "  (b) hex-pitch resample ({cw}×{ch}): largest {:.0} km² = {} hex | total {:.0} km²",
            b_res[ri],
            (b_res[ri] / HEX_KM2) as u64,
            tt as f32 * hex_cell_km2
        );
    }

    // ── CONVERGENCE VERDICT
    let ratio = |a: f32, b: f32| if a.min(b) > 0.0 { a.max(b) / a.min(b) } else { f32::INFINITY };
    eprintln!("\n=== CONVERGENCE (2048² vs 8192²; 1.00 = identical, lower is better) ===");
    eprintln!(
        "  NAIVE                : {:>8.0} vs {:>8.0} km²  → ×{:.2}",
        naive[0],
        naive[1],
        ratio(naive[0], naive[1])
    );
    for (i, &bm) in bridges_m.iter().enumerate() {
        eprintln!(
            "  (a) bridge {bm:>5.0} m    : {:>8.0} vs {:>8.0} km²  → ×{:.2}",
            a_res[0][i],
            a_res[1][i],
            ratio(a_res[0][i], a_res[1][i])
        );
    }
    eprintln!(
        "  (b) hex pitch        : {:>8.0} vs {:>8.0} km²  → ×{:.2}",
        b_res[0],
        b_res[1],
        ratio(b_res[0], b_res[1])
    );
    eprintln!(
        "  total flat area      : {:>8.0} vs {:>8.0} km²  → ×{:.2}",
        totals[0],
        totals[1],
        ratio(totals[0], totals[1])
    );
    eprintln!(
        "\n(A usable metric gives ~the SAME value at both resolutions. Report the smallest bridging distance that\n converges; if neither definition does, that IS the finding.)"
    );
}
