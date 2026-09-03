//! DIAGNOSTIC — Finding 42 follow-up. The exported terrain's HYPSOMETRY is not
//! resolution-stable: mean land altitude 287 m at 2048² against 693 m at 8192², p50 ×2.68,
//! p10 ×5.0, for the same seed, the same 400 km domain and near-identical emerged fraction
//! (14.95 % / 16.43 %). Everything altitude-dependent reads a different world per grid.
//!
//! THE BISECTION: the terrain is built in three stages, and this bench reports the full
//! hypsometry after EACH, at BOTH resolutions, through `production_hd_config` (so the
//! composition, the cap and the relief-v3 parameters are exactly the shipped ones):
//!
//!   1. COARSE ONLY  — `amplitude_base = 0`, no incision: the bilinearly upscaled C-1
//!                     isostatic altitude. The FBM contributes nothing (the amplitude
//!                     composition is `min(base, cap)` and base is 0).
//!   2. + FBM        — production amplitude, still no incision. The difference from (1) IS
//!                     the FBM's added relief, i.e. the OBSERVABLE EFFECT of the C-1 relief
//!                     budget cap. If that difference has the same distribution at both
//!                     resolutions, the cap is invariant — measured, not reconstructed.
//!   3. PRODUCTION   — the shipped config, incision included.
//!
//! Reading: inflation already present at (1) indicts the coarse upscale; appearing at (2)
//! indicts the FBM amplitude / the C-1 cap; appearing at (3) indicts stream-power / relief_v3.
//!
//! Run: cargo test -p ymir-core --release --test hypsometry_bisection -- --ignored --nocapture

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
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;

/// Which stage of the build to stop at.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Stage {
    CoarseOnly,
    PlusFbm,
    Production,
}

impl Stage {
    fn label(self) -> &'static str {
        match self {
            Stage::CoarseOnly => "1. COARSE ONLY (no FBM, no incision)",
            Stage::PlusFbm => "2. + FBM         (no incision)",
            Stage::Production => "3. PRODUCTION    (FBM + relief-v3 incision)",
        }
    }
}

struct Hyps {
    n_land: usize,
    emerged_pct: f64,
    mean_m: f64,
    p10: f32,
    p50: f32,
    p90: f32,
    p99: f32,
    max_m: f32,
    mean_norm_above: f64,
}

/// Hypsometry over the EMERGED cells only (the population every altitude-dependent
/// consumer reads). `mean_norm_above` is the raw normalised mean minus the 0.5 sea level —
/// an audit of the linear metric conversion.
fn hypsometry(field: &ymir_core::grid::GridF32, ss: &SteinSteinParams) -> Hyps {
    let mut alts: Vec<f32> = Vec::new();
    let mut sum_norm = 0.0f64;
    for &v in &field.data {
        if v > 0.5 {
            alts.push(c1_altitude_norm_to_metres(v, ss));
            sum_norm += (v - 0.5) as f64;
        }
    }
    alts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = alts.len().max(1);
    let q = |f: f32| alts[((n as f32 - 1.0) * f).max(0.0) as usize];
    Hyps {
        n_land: alts.len(),
        emerged_pct: 100.0 * alts.len() as f64 / field.data.len() as f64,
        mean_m: alts.iter().map(|&v| v as f64).sum::<f64>() / n as f64,
        p10: q(0.10),
        p50: q(0.50),
        p90: q(0.90),
        p99: q(0.99),
        max_m: *alts.last().unwrap_or(&0.0),
        mean_norm_above: sum_norm / n as f64,
    }
}

fn build(target: usize, stage: Stage) -> ymir_core::grid::GridF32 {
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

    // THE shipped config, then exactly one stage-isolating change.
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
    match stage {
        Stage::CoarseOnly => {
            cfg.amplitude_base = 0.0;
            cfg.stream_power = None;
        }
        Stage::PlusFbm => cfg.stream_power = None,
        Stage::Production => {}
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

fn row(h: &Hyps) -> String {
    format!(
        "{:>10} {:>7.2} {:>9.0} {:>7.0} {:>7.0} {:>7.0} {:>7.0} {:>7.0} {:>11.5}",
        h.n_land, h.emerged_pct, h.mean_m, h.p10, h.p50, h.p90, h.p99, h.max_m, h.mean_norm_above
    )
}

#[test]
#[ignore]
fn hypsometry_raw_vs_eroded() {
    let ss = SteinSteinParams::default();
    let stages = [Stage::CoarseOnly, Stage::PlusFbm, Stage::Production];
    let mut table: Vec<(Stage, usize, Hyps)> = Vec::new();
    // The FBM's added relief per cell — the OBSERVABLE effect of the C-1 cap.
    let mut fbm_delta: Vec<(usize, Vec<f32>)> = Vec::new();

    for target in [2048usize, 8192] {
        let mut coarse_field = None;
        for st in stages {
            let f = build(target, st);
            if st == Stage::CoarseOnly {
                coarse_field = Some(f.clone());
            }
            if st == Stage::PlusFbm {
                let c = coarse_field.as_ref().unwrap();
                let mut d: Vec<f32> = (0..f.data.len())
                    .map(|k| {
                        c1_altitude_norm_to_metres(f.data[k], &ss)
                            - c1_altitude_norm_to_metres(c.data[k], &ss)
                    })
                    .collect();
                d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                fbm_delta.push((target, d));
            }
            table.push((st, target, hypsometry(&f, &ss)));
        }
    }

    eprintln!(
        "\n==========  HYPSOMETRY BISECTION — seed {PSEED}, {DOMAIN_KM} km  ==========\n\
         over the EMERGED cells (altitude > sea level), metres unless stated\n"
    );
    eprintln!(
        "{:<44} {:>10} {:>7} {:>9} {:>7} {:>7} {:>7} {:>7} {:>7} {:>11}",
        "stage @ resolution",
        "land",
        "emerg%",
        "mean",
        "p10",
        "p50",
        "p90",
        "p99",
        "max",
        "norm>0.5"
    );
    for st in stages {
        for target in [2048usize, 8192] {
            let h = &table.iter().find(|(s, t, _)| *s == st && *t == target).unwrap().2;
            eprintln!("{:<44} {}", format!("{} @ {}²", st.label(), target), row(h));
        }
        // The ratio is the whole point: does THIS stage carry the inflation?
        let a = &table.iter().find(|(s, t, _)| *s == st && *t == 2048).unwrap().2;
        let b = &table.iter().find(|(s, t, _)| *s == st && *t == 8192).unwrap().2;
        eprintln!(
            "{:<44} {:>10} {:>7.2} {:>9.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>11.2}",
            "   ↳ 8192²/2048² RATIO",
            "-",
            b.emerged_pct / a.emerged_pct.max(1e-9),
            b.mean_m / a.mean_m.max(1e-9),
            b.p10 / a.p10.max(1e-6),
            b.p50 / a.p50.max(1e-6),
            b.p90 / a.p90.max(1e-6),
            b.p99 / a.p99.max(1e-6),
            b.max_m / a.max_m.max(1e-6),
            b.mean_norm_above / a.mean_norm_above.max(1e-9)
        );
        eprintln!();
    }

    eprintln!("--- THE FBM'S ADDED RELIEF (stage 2 − stage 1), over ALL cells, metres ---");
    eprintln!(
        "   this IS the observable effect of the C-1 relief budget cap \
         `β·slope_mag/(nscale·S)`.\n   Its inputs are COARSE-ONLY by construction \
         (`slope_map = compute_terrain_analysis(coarse)` on the 64²\n   grid, `nscale` from \
         `src_max`), so an identical distribution here CLEARS the cap.\n"
    );
    eprintln!(
        "{:>8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "n", "mean", "p1", "p50", "p99", "min", "max"
    );
    for (target, d) in &fbm_delta {
        let n = d.len();
        let q = |f: f32| d[((n as f32 - 1.0) * f) as usize];
        eprintln!(
            "{:>8} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10.2}",
            target,
            d.iter().map(|&v| v as f64).sum::<f64>() / n as f64,
            q(0.01),
            q(0.50),
            q(0.99),
            d[0],
            d[n - 1]
        );
    }
    // Diagnostic only: report, never gate.
}

/// Stage-3 SUB-BISECTION. The first bisection put the whole inflation in the incision
/// (stages 1 and 2 are invariant to 1.00 on every percentile). This isolates WHICH term of
/// `relief_v3` carries it: each variant changes exactly ONE thing from the shipped config,
/// and the diagnostic is the 8192Â²/2048Â² RATIO â€” the variant whose ratio collapses toward
/// 1.00 is the culprit.
#[derive(Clone, Copy)]
enum Variant {
    Shipped,
    MfdOff,
    NoDiffusion,
    NoTalus,
    Iters1,
    Iters4,
    Iters8,
    NoLateral,
    /// THE decisive one: remove the fluvial/hillslope REGIME SPLIT so every cell incises.
    /// `min_area_cells = A_c/cell_km2` is physically constant (0.1 km²) but the FRACTION of
    /// land cells that clears it is not: 55.4 % at 2048² against 10.1 % at 8192². If the
    /// ratio collapses here, the split is the carrier.
    NoRegimeSplit,
}

impl Variant {
    fn label(self) -> &'static str {
        match self {
            Variant::Shipped => "shipped relief-v3",
            Variant::MfdOff => "MFD off (D8 area)",
            Variant::NoDiffusion => "diffusion = 0",
            Variant::NoTalus => "talus off",
            Variant::Iters1 => "iterations = 1",
            Variant::Iters4 => "iterations = 4",
            Variant::Iters8 => "iterations = 8",
            Variant::NoLateral => "lateral_erosion = 0",
            Variant::NoRegimeSplit => "no regime split (A_c=0)",
        }
    }
    fn apply(self, sp: &mut ymir_core::erosion::stream_power::StreamPowerConfig) {
        match self {
            Variant::Shipped => {}
            Variant::MfdOff => sp.mfd_exponent = None,
            Variant::NoDiffusion => sp.diffusion = 0.0,
            Variant::NoTalus => sp.talus_slope = 0.0,
            Variant::Iters1 => sp.iterations = 1,
            Variant::Iters4 => sp.iterations = 4,
            Variant::Iters8 => sp.iterations = 8,
            Variant::NoLateral => sp.lateral_erosion = 0.0,
            Variant::NoRegimeSplit => sp.min_area_cells = 0.0,
        }
    }
}

fn build_variant(target: usize, v: Variant) -> ymir_core::grid::GridF32 {
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
    if let Some(sp) = cfg.stream_power.as_mut() {
        v.apply(sp);
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

#[test]
#[ignore]
fn incision_term_sub_bisection() {
    let ss = SteinSteinParams::default();
    // Pre-incision reference (stage 2), so "relief retained" is measurable per variant.
    let ref_2048 = hypsometry(&build(2048, Stage::PlusFbm), &ss).mean_m;
    let ref_8192 = hypsometry(&build(8192, Stage::PlusFbm), &ss).mean_m;
    eprintln!(
        "\n==========  STAGE-3 SUB-BISECTION â€” which incision term inflates?  ==========\n\
         pre-incision mean land altitude: {ref_2048:.0} m @ 2048Â² | {ref_8192:.0} m @ 8192Â² \
         (invariant)\n\n\
         'retained' = post-incision mean / pre-incision mean. The RATIO column is the \
         diagnostic:\n hunting a variant whose 8192Â²/2048Â² ratio collapses toward 1.00.\n"
    );
    eprintln!(
        "{:<22} {:>9} {:>9} {:>8} {:>9} {:>9} {:>8} {:>8}",
        "variant", "mean2048", "mean8192", "RATIO", "ret.2048", "ret.8192", "p50 2k", "p50 8k"
    );
    for v in [
        Variant::Shipped,
        Variant::MfdOff,
        Variant::NoDiffusion,
        Variant::NoTalus,
        Variant::NoLateral,
        Variant::NoRegimeSplit,
        Variant::Iters1,
        Variant::Iters4,
        Variant::Iters8,
    ] {
        let a = hypsometry(&build_variant(2048, v), &ss);
        let b = hypsometry(&build_variant(8192, v), &ss);
        eprintln!(
            "{:<22} {:>9.0} {:>9.0} {:>8.2} {:>8.1}% {:>8.1}% {:>8.0} {:>8.0}",
            v.label(),
            a.mean_m,
            b.mean_m,
            b.mean_m / a.mean_m.max(1e-9),
            100.0 * a.mean_m / ref_2048,
            100.0 * b.mean_m / ref_8192,
            a.p50,
            b.p50
        );
    }
    // Diagnostic only: report, never gate.
}

/// The decisive variant on its own (the full sub-bisection takes ~10 min; this is the pair
/// that settles the mechanism). `min_area_cells = 0` makes EVERY cell incise, removing the
/// fluvial/hillslope regime split â€” the one term whose *reach over the land* is
/// resolution-dependent even though its threshold is physically constant.
#[test]
#[ignore]
fn regime_split_is_the_carrier() {
    let ss = SteinSteinParams::default();
    let ref_2048 = hypsometry(&build(2048, Stage::PlusFbm), &ss).mean_m;
    let ref_8192 = hypsometry(&build(8192, Stage::PlusFbm), &ss).mean_m;
    eprintln!(
        "\n=====  REGIME SPLIT â€” the decisive variant  =====\n\
         pre-incision {ref_2048:.0} m @ 2048Â² | {ref_8192:.0} m @ 8192Â²\n"
    );
    eprintln!(
        "{:<26} {:>9} {:>9} {:>8} {:>9} {:>9}",
        "variant", "mean2048", "mean8192", "RATIO", "ret.2048", "ret.8192"
    );
    for v in [Variant::Shipped, Variant::NoRegimeSplit] {
        let a = hypsometry(&build_variant(2048, v), &ss);
        let b = hypsometry(&build_variant(8192, v), &ss);
        eprintln!(
            "{:<26} {:>9.0} {:>9.0} {:>8.2} {:>8.1}% {:>8.1}%",
            v.label(),
            a.mean_m,
            b.mean_m,
            b.mean_m / a.mean_m.max(1e-9),
            100.0 * a.mean_m / ref_2048,
            100.0 * b.mean_m / ref_8192
        );
    }
}

/// MECHANISM-2 SPLIT. With the regime partition removed (`A_c = 0`, mechanism 1 out of the
/// way) a residual Ã—1.85 remains â€” 37 % of the excess. Two candidates were visible in the
/// first sub-bisection and this separates them:
///
///  - **MFD dispersal**: the partition is applied PER CELL, so over the same physical path it
///    composes 4Ã— more often at 8192Â², diluting `A` and hence `f = KÂ·dtÂ·A^m/dist_m`.
///  - **Iteration saturation**: a fixed 2 sweeps, whose reach in the landscape is not a
///    physical quantity.
///
/// Each variant sits ON TOP of `A_c = 0`, so whatever moves the ratio toward 1.00 owns the
/// residual. `A_c = 0` remains a DIAGNOSTIC instrument and is never a candidate config (it
/// planes 2048Â² from 282 m to 172 m).
#[derive(Clone, Copy)]
enum Mech2 {
    NoSplit,
    NoSplitMfdOff,
    NoSplitIters1,
    NoSplitMfdOffIters1,
}

impl Mech2 {
    fn label(self) -> &'static str {
        match self {
            Mech2::NoSplit => "A_c=0 (residual baseline)",
            Mech2::NoSplitMfdOff => "A_c=0 + MFD off",
            Mech2::NoSplitIters1 => "A_c=0 + iterations=1",
            Mech2::NoSplitMfdOffIters1 => "A_c=0 + MFD off + iters=1",
        }
    }
    fn apply(self, sp: &mut ymir_core::erosion::stream_power::StreamPowerConfig) {
        sp.min_area_cells = 0.0;
        match self {
            Mech2::NoSplit => {}
            Mech2::NoSplitMfdOff => sp.mfd_exponent = None,
            Mech2::NoSplitIters1 => sp.iterations = 1,
            Mech2::NoSplitMfdOffIters1 => {
                sp.mfd_exponent = None;
                sp.iterations = 1;
            }
        }
    }
}

fn build_mech2(target: usize, v: Mech2) -> ymir_core::grid::GridF32 {
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
    if let Some(sp) = cfg.stream_power.as_mut() {
        v.apply(sp);
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

#[test]
#[ignore]
fn mechanism2_split() {
    let ss = SteinSteinParams::default();
    let ref_2048 = hypsometry(&build(2048, Stage::PlusFbm), &ss).mean_m;
    let ref_8192 = hypsometry(&build(8192, Stage::PlusFbm), &ss).mean_m;
    eprintln!(
        "\n=====  MECHANISM-2 SPLIT (MFD dilution vs iteration saturation)  =====\n\
         pre-incision {ref_2048:.0} m @ 2048Â² | {ref_8192:.0} m @ 8192Â²\n\
         every variant sits on A_c = 0, so mechanism 1 is out of the way. Whatever pulls the\n\
         RATIO toward 1.00 owns the residual (baseline Ã—1.85). 'excess' = mean8192 âˆ’ mean2048.\n"
    );
    eprintln!(
        "{:<30} {:>9} {:>9} {:>8} {:>9} {:>9} {:>9}",
        "variant", "mean2048", "mean8192", "RATIO", "excess m", "ret.2048", "ret.8192"
    );
    let mut base_excess = 0.0f64;
    for (i, v) in
        [Mech2::NoSplit, Mech2::NoSplitMfdOff, Mech2::NoSplitIters1, Mech2::NoSplitMfdOffIters1]
            .into_iter()
            .enumerate()
    {
        let a = hypsometry(&build_mech2(2048, v), &ss);
        let b = hypsometry(&build_mech2(8192, v), &ss);
        let excess = b.mean_m - a.mean_m;
        if i == 0 {
            base_excess = excess;
        }
        eprintln!(
            "{:<30} {:>9.0} {:>9.0} {:>8.2} {:>9.0} {:>8.1}% {:>8.1}%",
            v.label(),
            a.mean_m,
            b.mean_m,
            b.mean_m / a.mean_m.max(1e-9),
            excess,
            100.0 * a.mean_m / ref_2048,
            100.0 * b.mean_m / ref_8192
        );
        if i > 0 {
            eprintln!(
                "{:<30} â†’ removes {:.0}% of the residual excess ({:.0} m of {:.0} m)",
                "",
                100.0 * (base_excess - excess) / base_excess.max(1e-9),
                base_excess - excess,
                base_excess
            );
        }
    }
    // Diagnostic only: report, never gate.
}
