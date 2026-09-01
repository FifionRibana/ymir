//! Elucidate the identical-count anomaly: the H-1 benches inherited
//! `c1_hd_production`'s `amplitude_base = 0.16` while production overrides it to 0.04, yet
//! after adding the override the pre-breach lake footprint came out IDENTICAL TO THE CELL
//! (922 329). Either the parameter did not take effect, or it does not reach the terrain.
//! This builds the SAME terrain at both amplitudes and compares them directly.
//!
//! Run: cargo test -p ymir-core --test amplitude_anomaly -- --ignored --nocapture

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

fn build(target: usize, amplitude: f64) -> (GridF32, usize, usize) {
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
    cfg.amplitude_base = amplitude;
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
    let lake_cells = pre.lake_map.iter().filter(|&&v| v != 0).count();
    (raw, lake_cells, pre.lakes.len())
}

#[test]
#[ignore]
fn amplitude_reaches_the_terrain() {
    let target = 2048usize;
    eprintln!("\n=== AMPLITUDE ANOMALY — does amplitude_base reach the terrain? ===");
    let (a16, lc16, nl16) = build(target, 0.16);
    let (a04, lc04, nl04) = build(target, 0.04);

    let n = a16.data.len();
    let diff = a16.data.iter().zip(a04.data.iter()).filter(|(x, y)| (*x - *y).abs() > 1e-9).count();
    let maxd =
        a16.data.iter().zip(a04.data.iter()).map(|(x, y)| (x - y).abs()).fold(0.0f32, f32::max);
    let stat = |g: &GridF32| {
        let m: f32 = g.data.iter().sum::<f32>() / g.data.len() as f32;
        let v: f32 = g.data.iter().map(|x| (x - m) * (x - m)).sum::<f32>() / g.data.len() as f32;
        (m, v.sqrt())
    };
    let (m16, s16) = stat(&a16);
    let (m04, s04) = stat(&a04);
    eprintln!("heightmap 0.16: mean {m16:.6} sd {s16:.6} | 0.04: mean {m04:.6} sd {s04:.6}");
    eprintln!(
        "cells differing: {diff} / {n} ({:.2} %) | max |Δ| (norm) {maxd:.6} = {:.1} m",
        100.0 * diff as f32 / n as f32,
        maxd * 11302.0
    );
    eprintln!("pre-breach lake CELLS: 0.16 → {lc16} | 0.04 → {lc04}");
    eprintln!("pre-breach lake COUNT: 0.16 → {nl16} | 0.04 → {nl04}");
    if diff == 0 {
        eprintln!(
            "\n⇒ THE TERRAIN IS IDENTICAL: `amplitude_base` does NOT reach the terrain on this path.\n  The anomaly is real and the parameter is inert here — find where it is dropped."
        );
    } else {
        eprintln!(
            "\n⇒ The terrain DOES change with amplitude_base. The identical footprint counts in the\n  earlier run therefore came from a stale binary, not from an inert parameter."
        );
    }
}
