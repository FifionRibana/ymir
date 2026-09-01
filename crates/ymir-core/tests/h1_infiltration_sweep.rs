//! H-1 infiltration measurement (closures roadmap H-1). Runs the WHOLE production chain
//! (production terrain via `upscale_from_c1_with_progress` + relief-v3 breach + the real
//! climate + the real drainage) and reports, per CLIMATE and with infiltration OFF vs ON:
//!   - the lake population and the EXORHEIC/ENDORHEIC split — the direct measure of
//!     whether infiltration does real work;
//!   - total water area, and the level change on the ten largest lakes;
//!   - the infiltrated fraction ACTUALLY APPLIED per permeability class, against the
//!     published ranges (BFI 0.4–0.8) — THE SWEEP CRITERION;
//!   - CONTIGUOUS PLAIN AREA: land under 5° in ONE piece, excluding water, km² + hex count;
//!   - closed depressions (C-1's conditioning must survive).
//!
//! THE SWEEP CRITERION, stated before measuring: success is NOT "lakes decrease". It is
//! that the applied fraction stays INSIDE the published range AND the exorheic/endorheic
//! split moves. If the range must be left to obtain an effect, the finding is that
//! infiltration is NOT the dominant lever and sill incision (H-2) carries the dial.
//!
//! Run: cargo test -p ymir-core --test h1_infiltration_sweep h1_sweep -- --ignored --nocapture

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::infiltration::{
    InfiltrationConfig, K_FRACTURED, K_MATRIX_HARD, K_MATRIX_RIFT, K_MATRIX_VOLCANIC,
    build_hd_infiltration,
};
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, c1_drainage_windowed, c1_drainage_windowed_infil,
    classify_lakes_water_balance,
};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;
/// Climatic span (deg) kept NARROW so each latitude sits inside its own belt — the arid
/// test bed the campaign uses (lat 25 + span 10). Wider spans blur the four climates.
const SPAN_DEG: f32 = 10.0;
/// Hex assumption for the plain metric: a 1 km-EDGE hex → area = (3√3/2)·1² km². Reported
/// alongside km² so the author can rescale to the real Living Landz hex.
const HEX_KM2: f32 = 2.598_076_2;

fn count_closed_depressions(field: &GridF32) -> usize {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let thr = 0.1f32 / FULL_RANGE_M;
    let mut seen = vec![false; n];
    let (mut count, mut stack) = (0usize, Vec::new());
    for s in 0..n {
        if seen[s] || flow.filled.data[s] - field.data[s] <= thr {
            continue;
        }
        count += 1;
        seen[s] = true;
        stack.push(s);
        while let Some(k) = stack.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if !seen[nk] && flow.filled.data[nk] - field.data[nk] > thr {
                            seen[nk] = true;
                            stack.push(nk);
                        }
                    }
                }
            }
        }
    }
    count
}

/// Union-find (the `land_topology` pattern) over cells that are LAND, DRY (not lake) and
/// under `max_deg` of slope → the largest CONTIGUOUS PLAIN. Returns (largest_km2,
/// largest_hexes, total_flat_km2). Periodic 4-connectivity, like `land_topology`.
fn contiguous_plain(
    field: &GridF32,
    lake_map: &[u32],
    cell_km2: f32,
    km_per_cell: f32,
    max_deg: f32,
) -> (f32, u64, f32) {
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let cell_m = km_per_cell * 1000.0;
    let mut ok = vec![false; n];
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
    // Union-find over the flat-land mask.
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
            if !ok[k] {
                continue;
            }
            let r = j * w + (i + 1) % w;
            if ok[r] {
                union(&mut parent, k as u32, r as u32);
            }
            let d = ((j + 1) % h) * w + i;
            if ok[d] {
                union(&mut parent, k as u32, d as u32);
            }
        }
    }
    let mut sizes = std::collections::HashMap::new();
    let mut total = 0u64;
    for k in 0..n {
        if ok[k] {
            total += 1;
            let r = find(&mut parent, k as u32);
            *sizes.entry(r).or_insert(0u64) += 1;
        }
    }
    let largest = sizes.values().copied().max().unwrap_or(0);
    (
        largest as f32 * cell_km2,
        (largest as f32 * cell_km2 / HEX_KM2) as u64,
        total as f32 * cell_km2,
    )
}

/// Classify the SHIPPED (pre-breach, geometry-only) lakes with the climate water balance —
/// the H-1 missing link — and report the split + the area that flips endorheic.
/// `infil = None` isolates the reclassification alone; `Some` adds infiltration on top.
fn reclassify(
    field: &GridF32,
    flow: &ymir_core::terrain::flow::FlowResult,
    lakes: &[ymir_core::tectonics_c1::drainage::C1Lake],
    lake_map: &[u32],
    climate: &ymir_core::climate::ClimateResult,
    cell_km2: f32,
    infil: Option<&[f32]>,
) -> (usize, usize, f32, f32) {
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let v = classify_lakes_water_balance(
        field,
        flow,
        &dclim,
        cell_km2,
        lakes,
        lake_map,
        infil,
        field.width,
        field.height,
    );
    let endo = v.iter().filter(|x| x.lake_type == LakeType::Endorheic).count();
    let exo = v.iter().filter(|x| x.lake_type == LakeType::Exorheic).count();
    // area that becomes endorheic, and the area those basins would shrink to (REPORTED,
    // not applied — adopting the equilibrium level is H-2 territory).
    let endo_km2: f32 =
        v.iter().filter(|x| x.lake_type == LakeType::Endorheic).map(|x| x.a_sill_km2).sum();
    let endo_eq_km2: f32 =
        v.iter().filter(|x| x.lake_type == LakeType::Endorheic).map(|x| x.a_eq_km2).sum();
    (exo, endo, endo_km2, endo_eq_km2)
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

    // Production terrain, with the closures the author runs (C-3 lithology + C-3b fracture),
    // so the permeability field reads the same state production erodes.
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let litho = LithologyConfig {
        enabled: true,
        soft_multiplier: 10.0,
        volcanic_multiplier: 3.0,
        rift_age_threshold: 1.0,
    };
    let frac = FractureConfig {
        enabled: true,
        amplitude: 6.0,
        decay_km: 25.0,
        domain_km: DOMAIN_KM,
        ..Default::default()
    };
    let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
    cfg.sample_origin = [0.09375, 0.578125];
    cfg.sample_size = 1.0;
    cfg.erosion = None;
    cfg.lithology = litho.clone();
    cfg.fracture = frac.clone();
    cfg.stream_power = {
        let mut sp = ymir_core::erosion::stream_power::StreamPowerConfig::relief_v3(
            cell_km2,
            ss.depth_scale_m as f32,
        );
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;
        Some(sp)
    };
    let (up, _craters) = upscale_from_c1_with_progress(
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
    // Production breach conditioning (relief-v3) — the terrain the drainage actually sees.
    let raw = up.heightmap;
    let dcfg = {
        let mut d = C1DrainageConfig::default();
        d.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
        d.thresholds.full_tree = false;
        d
    };
    let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field = ymir_core::terrain::flow::breach_monotone(
        &raw,
        &pre.flow.filled,
        &pre.lake_map,
        0.5,
        raw.width,
        raw.height,
    );

    // The H-1 field, built exactly as production builds it.
    let icfg = InfiltrationConfig { enabled: true, ..Default::default() };
    let infil = build_hd_infiltration(
        &state,
        &kin,
        &litho,
        &frac,
        &icfg,
        &edifices,
        target,
        target,
        cfg.sample_origin,
        cfg.sample_size,
        km_per_cell,
    );

    eprintln!("\n================  H-1 INFILTRATION  {target}²  ================");
    // ── THE SWEEP CRITERION: the fraction actually applied, per permeability class,
    //    against the published range (BFI 0.4–0.8).
    let f_of = |k: f32| icfg.fraction_for_k(k);
    eprintln!(
        "applied fraction by class (f_cap {:.2}, k_ref {:.1e} m/day):\n  intact crystalline K={:.0e} → f={:.3} | rift K={:.0e} → f={:.3} | volcaniclastic K={:.0e} → f={:.3} | fully fractured K={:.0e} → f={:.3}",
        icfg.f_cap,
        icfg.k_ref_m_per_day,
        K_MATRIX_HARD,
        f_of(K_MATRIX_HARD),
        K_MATRIX_RIFT,
        f_of(K_MATRIX_RIFT),
        K_MATRIX_VOLCANIC,
        f_of(K_MATRIX_VOLCANIC),
        K_FRACTURED,
        f_of(K_FRACTURED),
    );
    let mut fs: Vec<f32> = infil.clone();
    fs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f32| fs[((fs.len() as f32 * p) as usize).min(fs.len() - 1)];
    let mean: f32 = infil.iter().sum::<f32>() / infil.len() as f32;
    eprintln!(
        "  field-wide f: mean {mean:.3} | p10 {:.3} | p50 {:.3} | p90 {:.3} | p99 {:.3}   [published BFI 0.4–0.8]",
        q(0.1),
        q(0.5),
        q(0.9),
        q(0.99)
    );

    // The SHIPPED surface-lake population: the pre-breach, CLIMATE-FREE lakes production
    // carries into the final result (all exorheic by pure geometry — no water balance).
    let geo_exo = pre.lakes.iter().filter(|l| l.lake_type == LakeType::Exorheic).count();
    let geo_endo = pre.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).count();
    let geo_km2: f32 = pre.lakes.iter().map(|l| l.area_km2).sum();
    let (plain_km2, plain_hex, flat_total) =
        contiguous_plain(&field, &pre.lake_map, cell_km2, km_per_cell, 5.0);
    eprintln!(
        "\nSHIPPED surface lakes (pre-breach, geometry only, NO water balance): {} lakes = {:.0} km² — exo {geo_exo} / endo {geo_endo}",
        pre.lakes.len(),
        geo_km2
    );
    eprintln!(
        "closed depressions (C-1): {} | contiguous plain <5° (dry land): largest {:.0} km² = {} hex | all flat {:.0} km²",
        count_closed_depressions(&field),
        plain_km2,
        plain_hex,
        flat_total
    );
    eprintln!(
        "\n{:>11} | {:>22} | {:>22}",
        "climate", "reclassified (no infil)", "reclassified + H-1 infil"
    );
    eprintln!("{:>11} | {:>22} | {:>22}", "", "exo/endo  endo km²", "exo/endo  endo km²");
    for (label, lat) in
        [("tropical", 10.0f32), ("arid-hot", 25.0), ("humid", 45.0), ("arid-cold", 65.0)]
    {
        let climate =
            c1_climate_placed(&field, &ss, lat, SPAN_DEG, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let fin = c1_drainage_windowed_infil(&field, Some(&dclim), &dcfg, &ss, DOMAIN_KM, None);
        let (e0, n0, a0, q0) =
            reclassify(&field, &fin.flow, &pre.lakes, &pre.lake_map, &climate, cell_km2, None);
        let (e1, n1, a1, q1) = reclassify(
            &field,
            &fin.flow,
            &pre.lakes,
            &pre.lake_map,
            &climate,
            cell_km2,
            Some(infil.as_slice()),
        );
        eprintln!(
            "{label:>11} | {e0:>4}/{n0:<4} {a0:>9.0} | {e1:>4}/{n1:<4} {a1:>9.0}   (endorheic basins would shrink to {q0:.0} → {q1:.0} km² — reported, NOT applied)"
        );
    }
    eprintln!(
        "(CRITERION: the applied fraction must stay inside the published BFI range AND the exo/endo split\n must move. If an effect needs leaving the range, infiltration is NOT the dominant lever — H-2 carries it.)"
    );
}

#[test]
#[ignore]
fn h1_sweep() {
    run(2048);
}

#[test]
#[ignore]
fn h1_sweep_8192() {
    run(8192);
}
