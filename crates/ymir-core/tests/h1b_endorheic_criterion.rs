//! H-1b — VERIFY the endorheic criterion before H-2 is dimensioned on it. 40 endorheic out
//! of 67 in a HUMID climate is physically suspect: an endorheic basin means evaporation
//! absorbs the ENTIRE catchment inflow (Dead Sea / Great Salt Lake regime), and France,
//! Scotland and Scandinavia have essentially none.
//!
//! TWO THINGS ARE CHECKED, on the lakes the balance declares endorheic in humid:
//!
//! 1. **Is the AREA test equivalent to the physical LEVEL test?** `a_eq ≥ a_sill` compares
//!    the evaporative equilibrium AREA with the area AT THE SILL. Since the hypsometry is
//!    monotone (area grows with level), that SHOULD be equivalent to "the level reaches the
//!    sill" — this reports both so the equivalence is measured, not assumed.
//!
//! 2. **Is the EVAPORATION term right?** The code uses `a_eq = inflow / PE` — the FULL
//!    potential evaporation. But a lake also RECEIVES precipitation on its own surface, so
//!    the net loss per unit area is `PE − P`, not `PE`. Using full PE overstates the loss,
//!    shrinks `a_eq`, and pushes basins toward endorheic — massively in a humid climate
//!    where `PE − P → 0`. That is exactly the `net_evap` of Finding 39, absent from the
//!    surface path. This reports the verdict under BOTH formulations.
//!
//! Consequence either way: this does not remove H-2, it SIZES it. A stricter-than-physics
//! criterion means MORE exorheic lakes, hence more sills to incise.
//!
//! Run: cargo test -p ymir-core --test h1b_endorheic_criterion -- --ignored --nocapture

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, c1_drainage_windowed, potential_evaporation_mm,
    runoff_accumulation,
};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::{FlowConfig, compute_flow};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const SPAN_DEG: f32 = 10.0;
const FULL_RANGE_M: f32 = 11302.0;
/// Normalised altitude → metres (sea at 0.5, the C1 vertical contract).
fn to_m(norm: f32) -> f32 {
    (norm - 0.5) * FULL_RANGE_M
}

fn build(target: usize) -> (GridF32, Vec<u32>, Vec<ymir_core::tectonics_c1::drainage::C1Lake>) {
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
    (field, pre.lake_map, pre.lakes)
}

fn audit(target: usize, lat: f32, label: &str) {
    let ss = SteinSteinParams::default();
    let km_per_cell = DOMAIN_KM / target as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let (field, lake_map, lakes) = build(target);
    let (w, h) = (field.width, field.height);
    let climate =
        c1_climate_placed(&field, &ss, lat, SPAN_DEG, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let flow =
        compute_flow(&field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let runoff = runoff_accumulation(&field, &flow, &dclim, cell_km2, None, None, w, h);

    eprintln!(
        "\n================  ENDORHEIC CRITERION AUDIT — {label} {lat}° @ {target}²  ================"
    );
    eprintln!(
        "{:>5} | {:>9} | {:>7} {:>7} {:>8} | {:>9} {:>9} | {:>9} | {:>8} {:>8} | verdict",
        "lake",
        "inflow",
        "PE",
        "P",
        "PE−P",
        "a_eq(PE)",
        "a_eq(net)",
        "a_sill",
        "lvl_eq",
        "lvl_sill"
    );
    let (mut endo_pe, mut endo_net, mut n) = (0usize, 0usize, 0usize);
    let (mut area_pe, mut area_net) = (0.0f32, 0.0f32);
    let mut shown = 0usize;
    for lk in &lakes {
        let mut cells: Vec<(usize, f32)> =
            (0..w * h).filter(|&k| lake_map[k] == lk.base.id).map(|k| (k, field.data[k])).collect();
        if cells.is_empty() {
            continue;
        }
        cells.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        n += 1;
        let inflow = cells.iter().map(|&(k, _)| runoff[k]).fold(0.0f32, f32::max);
        let floor_k = cells[0].0;
        let pe = potential_evaporation_mm(climate.temperature.data[floor_k]).max(1.0);
        // Precipitation falling ON the lake surface (mean over its cells) — the term the
        // current criterion omits.
        let p: f32 = cells
            .iter()
            .map(|&(k, _)| precip_mm_per_year(climate.precipitation.data[k]))
            .sum::<f32>()
            / cells.len() as f32;
        let net = (pe - p).max(1.0); // net loss per unit lake area
        let a_eq_pe = inflow / pe;
        let a_eq_net = inflow / net;
        let a_sill = lk.base.area as f32 * cell_km2;
        let is_endo_pe = a_eq_pe < a_sill;
        let is_endo_net = a_eq_net < a_sill;
        if is_endo_pe {
            endo_pe += 1;
            area_pe += a_sill;
        }
        if is_endo_net {
            endo_net += 1;
            area_net += a_sill;
        }
        // Equilibrium LEVEL from the hypsometry vs the SILL level — the physical test.
        let n_eq = ((a_eq_pe / cell_km2).floor() as usize).clamp(1, cells.len());
        let lvl_eq = to_m(cells[n_eq - 1].1);
        let lvl_sill = to_m(lk.base.surface_elevation);
        if is_endo_pe && shown < 12 {
            shown += 1;
            eprintln!(
                "{:>5} | {inflow:>9.0} | {pe:>7.0} {p:>7.0} {net:>8.0} | {a_eq_pe:>9.2} {a_eq_net:>9.2} | {a_sill:>9.2} | {lvl_eq:>8.0} {lvl_sill:>8.0} | {}",
                lk.base.id,
                if is_endo_net { "endo both" } else { "ENDO(PE) → EXO(net)" }
            );
        }
    }
    eprintln!(
        "\n  lakes {n} | ENDORHEIC with FULL PE (current): {endo_pe} ({area_pe:.0} km²)  →  with NET (PE−P): {endo_net} ({area_net:.0} km²)"
    );
    eprintln!(
        "  ⇒ {} lake(s) flip back to EXORHEIC once the precipitation falling ON the lake is counted.",
        endo_pe.saturating_sub(endo_net)
    );
}

#[test]
#[ignore]
fn endorheic_criterion_humid() {
    audit(2048, 45.0, "humid");
}

#[test]
#[ignore]
fn endorheic_criterion_all_climates() {
    for (label, lat) in
        [("tropical", 10.0f32), ("arid-hot", 25.0), ("humid", 45.0), ("arid-cold", 65.0)]
    {
        audit(2048, lat, label);
    }
}
