//! C-2 crater lakes — the FOUR-CLIMATE table at 8192² on the PRODUCTION terrain
//! (`upscale_from_c1_with_progress`, the real generation function, not a raster
//! reconstruction). Per climate: active craters, craters holding a lake, their
//! area/depth, and the inflow + a_eq behind each. States plainly whether the
//! proportion is plausible (a majority of active craters holding is a RED FLAG).
//!
//! Run: cargo test -p ymir-core --test c2_crater_8192 --release -- --ignored --nocapture

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{
    VolcanismConfig, detect_crater_lakes, place_edifices,
};
use ymir_core::tectonics_c1::drainage::{
    DrainageClimate, LakeType, potential_evaporation_mm, runoff_accumulation,
};
use ymir_core::terrain::flow::{FlowConfig, FlowResult, compute_flow};

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const TARGET: usize = 8192;

fn climate_row(
    label: &str,
    lat: f32,
    span: f32,
    eroded: &GridF32,
    flow: &FlowResult,
    craters: &[ymir_core::tectonics_c1::closures::volcanism::CraterRecord],
    ss: &SteinSteinParams,
    cell_km2: f32,
) {
    let (w, h) = (eroded.width, eroded.height);
    let climate = c1_climate_placed(eroded, ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let runoff = runoff_accumulation(eroded, flow, &dclim, cell_km2, None, w, h);

    // detect on a fresh (zero) lake_map so we measure ONLY the crater lakes.
    let mut lakes = Vec::new();
    let mut lake_map = vec![0u32; w * h];
    let (held, dry) = detect_crater_lakes(
        eroded,
        &flow.filled,
        &runoff,
        &climate.temperature,
        craters,
        cell_km2,
        ss,
        &mut lakes,
        &mut lake_map,
    );
    let active = craters.iter().filter(|c| c.active).count();

    // per-active-crater inflow + a_eq (max over the footprint, as the balance uses).
    let mut per = Vec::new();
    for c in craters.iter().filter(|c| c.active) {
        let (cx, cy, rc) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
        let (i0, i1) = ((cx - rc) as usize, ((cx + rc) as usize).min(w - 1));
        let (j0, j1) = ((cy - rc) as usize, ((cy + rc) as usize).min(h - 1));
        let (mut inflow, mut floor_k, mut floor_v) = (0.0f32, 0usize, f32::MAX);
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d > rc {
                    continue;
                }
                let k = j * w + i;
                inflow = inflow.max(runoff[k]);
                if eroded.data[k] < floor_v {
                    floor_v = eroded.data[k];
                    floor_k = k;
                }
            }
        }
        let a_eq = inflow / potential_evaporation_mm(climate.temperature.data[floor_k]).max(1.0);
        per.push((2.0 * rc * (DOMAIN_KM / TARGET as f32), inflow, a_eq, a_eq / cell_km2));
    }

    eprintln!(
        "\n--- {label} (lat {lat}°) --- active {active} | holding {held} | dry {dry}  {}",
        if held * 2 > active { "*** MAJORITY of active — RED FLAG ***" } else { "(minority)" }
    );
    let mut la: Vec<&_> = lakes.iter().collect();
    la.sort_by(|a, b| b.area_km2.partial_cmp(&a.area_km2).unwrap());
    for l in &la {
        eprintln!("    lake {:.2} km² · {:.0} m · {:?}", l.area_km2, l.depth_m, l.lake_type);
    }
    for (dia, inflow, a_eq, a_eq_cells) in &per {
        eprintln!(
            "    active crater Ø{dia:.2}km | inflow {inflow:.1} | a_eq {a_eq:.3} km² = {a_eq_cells:.0} cells"
        );
    }
}

#[test]
#[ignore]
fn c2_crater_climates() {
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::FbmUpscaleConfig;

    let ss = SteinSteinParams::default();
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let cell_km2 = (DOMAIN_KM / TARGET as f32).powi(2);

    // Production terrain (relief-v3 + reconstruction), built once (climate-free).
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = FbmUpscaleConfig::c1_hd_production(TARGET);
    cfg.amplitude_base = 0.04;
    cfg.flow_conditioning = 0.1;
    cfg.sample_origin = [0.09375, 0.578125];
    cfg.sample_size = 1.0;
    cfg.erosion = None;
    cfg.stream_power = {
        let mut sp = ymir_core::erosion::stream_power::StreamPowerConfig::relief_v3(
            cell_km2,
            ss.depth_scale_m as f32,
        );
        sp.mfd_exponent = Some(2.0);
        sp.iterations = 2;
        Some(sp)
    };
    let (up, craters) = upscale_from_c1_with_progress(
        &state,
        &run.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        None,
        &mut |_| {},
        &|| false,
    );
    // Match production: relief-v3 breach with ACTIVE crater bowls PROTECTED, so the
    // measurement reflects the shipped path (the breach otherwise re-breaches them).
    let raw = up.heightmap;
    let (rw, rh) = (raw.width, raw.height);
    let prebreach = ymir_core::tectonics_c1::drainage::c1_drainage_windowed(
        &raw,
        None,
        &ymir_core::tectonics_c1::drainage::C1DrainageConfig::default(),
        &ss,
        DOMAIN_KM,
    );
    let mut protect = vec![false; rw * rh];
    for c in craters.iter().filter(|c| c.active) {
        let (cx, cy, r) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
        let (i0, i1) = ((cx - r).floor().max(0.0) as usize, ((cx + r).ceil() as usize).min(rw - 1));
        let (j0, j1) = ((cy - r).floor().max(0.0) as usize, ((cy + r).ceil() as usize).min(rh - 1));
        for j in j0..=j1 {
            for i in i0..=i1 {
                if ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt() <= r {
                    protect[j * rw + i] = true;
                }
            }
        }
    }
    let eroded = ymir_core::terrain::flow::breach_monotone_protected(
        &raw,
        &prebreach.flow.filled,
        &prebreach.lake_map,
        0.5,
        rw,
        rh,
        Some(&protect),
    );
    let flow =
        compute_flow(&eroded, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });

    eprintln!(
        "\n=== C-2 CRATER LAKES per climate (production {TARGET}²) ===\n\
         cell {:.0} m; {} craters ({} active). Thresholds: general lake_min_area 5 km² = {:.0} cells; \
         crater floor CRATER_LAKE_MIN_CELLS 4 = {:.3} km²",
        DOMAIN_KM / TARGET as f32 * 1000.0,
        craters.len(),
        craters.iter().filter(|c| c.active).count(),
        5.0 / cell_km2,
        4.0 * cell_km2,
    );
    for (label, lat, span) in [
        ("arid-hot", 25.0f32, 10.0f32),
        ("humid", 45.0, 40.0),
        ("tropical", 10.0, 20.0),
        ("arid-cold(export cfg)", 65.0, 10.0),
    ] {
        climate_row(label, lat, span, &eroded, &flow, &craters, &ss, cell_km2);
    }
    let _ = LakeType::CraterNeutral;
}
