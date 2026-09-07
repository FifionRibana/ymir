//! ADR Finding 56 step 1 — MEASURE the slope at which the current `A_c` is correct.
//!
//! `A_c(S) = C/S²` must be calibrated on the regime where the present threshold works: the
//! HILLSLOPES. `C = A_c_ref · S_ref²`, so the apron inherits the right ratio without anyone
//! choosing a factor — calibrating on the apron instead would mean recalibrating two regimes
//! and would reintroduce the arbitrary multiplier the uniform sweep just refuted.
//!
//! `S_ref` is therefore the slope AT THE CHANNEL HEADS the model currently produces: cells whose
//! accumulation sits just at `A_c = 0.1 km²`. Measured, not rounded.
//!
//! Run: cargo test -p ymir-core --release --test channel_head_slope -- --ignored --nocapture

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
use ymir_core::terrain::flow::{FlowConfig, compute_flow};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;

fn production(target: usize) -> GridF32 {
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
    let cfg = production_hd_config(&ProductionHdOpts {
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

#[test]
#[ignore]
fn measure_channel_head_slope() {
    let ss = SteinSteinParams::default();
    let a_c = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    eprintln!(
        "\n=====  S_ref — the slope at the CURRENT channel heads (A_c = {a_c} km²)  =====\n\
         Calibrating C = A_c · S_ref² on the HILLSLOPE regime, where the present threshold is\n\
         correct, so the apron inherits the ratio the law imposes rather than a chosen factor."
    );
    for target in [2048usize, 8192] {
        let f = production(target);
        let (w, h) = (f.width, f.height);
        let km = DOMAIN_KM / target as f32;
        let m_per_cell = km * 1000.0;
        let cell_km2 = km * km;
        let head_cells = a_c / cell_km2;
        let flow = compute_flow(&f, &FlowConfig::default());
        let to_m = |v: f32| c1_altitude_norm_to_metres(v, &ss);

        // Channel HEADS: accumulation within ±25 % of A_c, on land. Their slope is the
        // calibration point — it is where the model currently decides a channel begins.
        let mut heads: Vec<f32> = Vec::new();
        // …and the whole LAND slope distribution, for the S_min discussion.
        let mut land: Vec<f32> = Vec::new();
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let k = y * w + x;
                if f.data[k] <= 0.5 {
                    continue;
                }
                let gx = 0.5 * (to_m(f.data[k + 1]) - to_m(f.data[k - 1]));
                let gy = 0.5 * (to_m(f.data[k + w]) - to_m(f.data[k - w]));
                let s = (gx * gx + gy * gy).sqrt() / m_per_cell; // dimensionless gradient
                land.push(s);
                let a = flow.accumulation.data[k];
                if a >= head_cells * 0.75 && a <= head_cells * 1.25 {
                    heads.push(s);
                }
            }
        }
        let q = |v: &mut Vec<f32>, p: f32| {
            v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            v[(((v.len() - 1) as f64 * p as f64).round() as usize).min(v.len() - 1)]
        };
        let n = heads.len();
        let (h25, h50, h75) = (q(&mut heads, 0.25), q(&mut heads, 0.50), q(&mut heads, 0.75));
        eprintln!(
            "\n--- {target}² ---\n  channel-head cells: {n}\n  \
             S at the heads: p25 {:.4} ({:.2}°) | p50 {:.4} ({:.2}°) | p75 {:.4} ({:.2}°)\n  \
             ⇒ C = A_c · S_p50² = {a_c} × {:.4}² = {:.3e} km²",
            h25,
            h25.atan().to_degrees(),
            h50,
            h50.atan().to_degrees(),
            h75,
            h75.atan().to_degrees(),
            h50,
            a_c * h50 * h50
        );
        let (l01, l05, l10, l50) =
            (q(&mut land, 0.01), q(&mut land, 0.05), q(&mut land, 0.10), q(&mut land, 0.50));
        eprintln!(
            "  LAND slope distribution: p1 {:.4} ({:.2}°) | p5 {:.4} ({:.2}°) | p10 {:.4} ({:.2}°) \
             | p50 {:.4} ({:.2}°)",
            l01,
            l01.atan().to_degrees(),
            l05,
            l05.atan().to_degrees(),
            l10,
            l10.atan().to_degrees(),
            l50,
            l50.atan().to_degrees()
        );
        // What each candidate S_min would cost: the share of land below it, and the A_c it caps at.
        eprintln!("  S_min candidates — share of land below, and the A_c the cap implies:");
        for smin_deg in [0.5f32, 1.0, 1.5, 2.0] {
            let smin = smin_deg.to_radians().tan();
            let below = land.iter().filter(|&&s| s < smin).count();
            eprintln!(
                "    {smin_deg:>4.1}° (S = {smin:.4}) → {:>5.1} % of land | A_c capped at {:.2} km²",
                100.0 * below as f32 / land.len() as f32,
                a_c * (h50 / smin).powi(2)
            );
        }
    }
}
