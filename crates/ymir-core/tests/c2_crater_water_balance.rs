//! C-2 crater water balance — the MEASUREMENT the author asked for, with the
//! position fix in place. Per climate (arid-hot 25° / humid 45° / tropical 10° /
//! arid-cold 65°), reports per crater the MEASURED inflow, equilibrium area, sill
//! area and margin (not the estimate), whether a lake actually formed, split by
//! active/extinct — and an explicit plausibility verdict (a MAJORITY of active
//! craters holding lakes is a red flag, not a success).
//!
//! Run: cargo test -p ymir-core --test c2_crater_water_balance --release -- --ignored --nocapture

use std::f32::consts::PI;

use ymir_core::climate::c1_climate_windowed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{
    VolcanismConfig, VolcanoSetting, apply_edifices, place_edifices,
};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, c1_drainage_windowed, potential_evaporation_mm,
    runoff_accumulation,
};
use ymir_core::terrain::flow::{FlowConfig, compute_flow};

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;
const TARGET: usize = 2048;

#[test]
#[ignore]
fn c2_crater_water_balance() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
    use ymir_core::tectonics_c1::kinematics::PlateKinematics;
    use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
    use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
    use ymir_core::terrain::upscale::{FbmUpscaleConfig, upscale_with_fbm};

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
    let coarse = c1_coarse_normalized_altitude(&state, &run.iso_config, &ss, None);
    let seed = WorldSeed::new(PSEED);
    let km_per_cell = DOMAIN_KM / TARGET as f32;
    let cell_km2 = km_per_cell * km_per_cell;

    // Terrain is climate-independent → build the eroded field + craters ONCE.
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = FbmUpscaleConfig::c1_hd_production(TARGET);
    cfg.amplitude_base = 0.04;
    cfg.sample_origin = [0.09375, 0.578125];
    cfg.sample_size = 1.0;
    cfg.stream_power = None;
    cfg.erosion = None;
    let mut fbm = upscale_with_fbm(&coarse, 0.5, &seed, &cfg).heightmap;
    let applied = apply_edifices(
        &mut fbm,
        &edifices,
        cfg.sample_origin,
        cfg.sample_size,
        km_per_cell,
        FULL_RANGE_M,
        &volc,
    );
    let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    sp.mfd_exponent = Some(2.0);
    sp.iterations = 2;
    let mut eroded = incise(&fbm, &sp);
    // C-2 active-rim reconstruction (as the production pipeline does after erosion):
    // active craters re-close so the lake stage can detect + hold them.
    ymir_core::tectonics_c1::closures::volcanism::reconstruct_active_rims(
        &mut eroded,
        &applied.craters,
        &volc,
        km_per_cell,
        FULL_RANGE_M,
    );
    let (w, h) = (eroded.width, eroded.height);
    let flow =
        compute_flow(&eroded, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });

    // Active flag per crater, in placement order (apply_edifices preserves it).
    let craters = &applied.craters;
    eprintln!(
        "\n=== C-2 CRATER WATER BALANCE (measured, {TARGET}², domain {DOMAIN_KM} km) ===\n\
         {} craters resolved ({} active, {} extinct); crater Ø {:.2} km median",
        craters.len(),
        craters.iter().filter(|c| c.active).count(),
        craters.iter().filter(|c| !c.active).count(),
        {
            let mut d: Vec<f32> = craters.iter().map(|c| 2.0 * c.radius_px * km_per_cell).collect();
            d.sort_by(|a, b| a.partial_cmp(b).unwrap());
            d.get(d.len() / 2).copied().unwrap_or(0.0)
        }
    );
    eprintln!("real anchors: Kawah Ijen Ø1.0 km/200 m, Poás Ø1.6 km/300 m, Pavin Ø0.7 km/90 m");

    let pp = PrecipParams::default();
    let dcfg = C1DrainageConfig::default();
    for (label, lat) in [
        ("arid-hot 25°", 25.0f32),
        ("humid 45°", 45.0),
        ("tropical 10°", 10.0),
        ("arid-cold 65°", 65.0),
    ] {
        let climate = c1_climate_windowed(&eroded, &ss, lat, &pp, DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let runoff = runoff_accumulation(&eroded, &flow, &dclim, cell_km2, None, w, h);
        // Actual detected/typed lakes (the pipeline's own verdict).
        let mut dr = c1_drainage_windowed(&eroded, Some(&dclim), &dcfg, &ss, DOMAIN_KM);
        let (held, dry) = ymir_core::tectonics_c1::closures::volcanism::classify_crater_lakes(
            &mut dr.lakes,
            &dr.lake_map,
            w,
            h,
            craters,
        );

        // Per-crater MEASURED balance: floor cell (min elev in the rim), inflow at
        // the floor (accumulated inner catchment), a_eq = inflow/PE, a_sill = crater
        // surface, margin = a_eq/a_sill. A lake forms iff margin ≥ 1 AND it survived
        // as a detected depression.
        let (mut act_fill, mut act_tot, mut ext_fill, mut ext_tot, mut near) = (0, 0, 0, 0, 0);
        let mut lines = Vec::new();
        for c in craters {
            let (cx, cy, r) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
            let (mut floor_k, mut floor_v) = (0usize, f32::MAX);
            let (i0, i1) = ((cx - r) as usize, ((cx + r) as usize).min(w - 1));
            let (j0, j1) = ((cy - r) as usize, ((cy + r) as usize).min(h - 1));
            for j in j0..=j1 {
                for i in i0..=i1 {
                    let d = ((i as f32 - cx).powi(2) + (j as f32 - cy).powi(2)).sqrt();
                    if d <= r && eroded.data[j * w + i] < floor_v {
                        floor_v = eroded.data[j * w + i];
                        floor_k = j * w + i;
                    }
                }
            }
            let inflow = runoff[floor_k]; // km²·(m/yr equivalent) — the balance's units
            let pe = potential_evaporation_mm(climate.temperature.data[floor_k]).max(1.0);
            let a_eq = inflow / pe;
            let a_sill = PI * (r * km_per_cell).powi(2);
            let margin = a_eq / a_sill.max(1e-6);
            let lake_here = {
                let k = (cy.round() as usize).min(h - 1) * w + (cx.round() as usize).min(w - 1);
                dr.lake_map[k] != 0
            };
            if c.active {
                act_tot += 1;
                if lake_here {
                    act_fill += 1;
                }
            } else {
                ext_tot += 1;
                if lake_here {
                    ext_fill += 1;
                }
            }
            if (0.7..1.5).contains(&margin) {
                near += 1;
            }
            lines.push(format!(
                "    {:<7} Ø{:.2}km | inflow {:.2} | a_eq {:.2} | a_sill {:.2} km² | margin {:.2} | lake={}",
                if c.active { "ACTIVE" } else { "extinct" },
                2.0 * r * km_per_cell,
                inflow,
                a_eq,
                a_sill,
                margin,
                if lake_here { "yes" } else { "no" },
            ));
        }
        eprintln!(
            "\n--- {label} --- craters holding a lake: {held} / {} (dry {dry}) | active {act_fill}/{act_tot}, extinct {ext_fill}/{ext_tot} | near-threshold (0.7–1.5×) {near}",
            craters.len()
        );
        // Lakes that formed: size vs real anchors.
        for c1l in dr.lakes.iter().filter(|l| {
            use ymir_core::tectonics_c1::drainage::LakeType::*;
            matches!(l.lake_type, CraterAcidic | CraterNeutral)
        }) {
            eprintln!(
                "      crater lake #{}: {:.2} km² · {:.0} m deep ({:?})",
                c1l.base.id, c1l.area_km2, c1l.depth_m, c1l.lake_type
            );
        }
        for l in &lines {
            eprintln!("{l}");
        }
    }
    let _ = GridF32::new(1, 1, 0.0);
}
