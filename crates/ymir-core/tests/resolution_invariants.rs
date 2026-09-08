//! ADR block F — guards for the two quantities block A/C **attributed**, and nothing else.
//!
//! 1. **The pre-incision field is resolution-invariant** in mean altitude, altitude quantiles
//!    and emerged fraction. Measured at 865.4 m / 865.5 m and 16.88 % / 16.88 % over a 4×
//!    resolution step (block A). This is the fact that killed "the upscale is the root cause",
//!    so if it ever stops holding, that conclusion must be revisited rather than inherited.
//!
//! 2. **A physically dimensioned AREA constant must not be sub-cell** at the coarsest
//!    resolution in its validity domain. `A_c = 0.1 km²` is 2.62 cells at 2048² — it passes
//!    every dimensional audit and is still not a threshold, because a channel head cannot be
//!    resolved by 2.6 cells. This is the "physically dimensioned but sub-cell" pattern (ADR
//!    method catalogue), distinct from "cells instead of metres" and invisible to it.
//!
//! Run: `cargo test -p ymir-core --release --test resolution_invariants`

use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
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
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;

/// The coarsest grid in the validity domain. 2048² is the project's calibration resolution
/// (`HILLSLOPE_REF_CELL_M` is defined as `400 km / 2048`), so it is the grid at which every
/// area constant has to remain resolvable.
const COARSEST_TARGET: usize = 2048;

/// The floor, in cells, below which an area constant is not a threshold on the coarsest grid.
///
/// **DECLARED, not derived.** 10 cells is the smallest footprint on which a channel head can be
/// distinguished from its own discretisation: a D8 neighbourhood is 9 cells, so anything at or
/// under that is a single-neighbourhood object and the "threshold" is really the stencil. It is
/// a stated convention, and `A_c = 0.1 km²` fails it at 2.62 — which is the point.
const MIN_AREA_CELLS_FLOOR: f32 = 10.0;

fn no_incision_terrain(target: usize) -> GridF32 {
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
    cfg.stream_power = None; // pre-incision: the quantity being guarded
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

/// `(mean m, p50 m, p90 m, land fraction)` — mean in **f64**: an f32 running sum over the land
/// cells of an 8192² grid saturates and read 13 m low (Finding 56b).
fn stats(f: &GridF32, ss: &SteinSteinParams) -> (f64, f32, f32, f64) {
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, ss);
    let mut a: Vec<f32> = f.data.iter().filter(|&&x| x > SEA).map(|&x| to_m(x)).collect();
    a.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let n = a.len().max(1);
    let pick = |q: f64| a[(((n - 1) as f64 * q).round() as usize).min(n - 1)];
    (
        a.iter().map(|&x| x as f64).sum::<f64>() / n as f64,
        pick(0.50),
        pick(0.90),
        a.len() as f64 / f.data.len() as f64,
    )
}

/// F1 — the pre-incision field does not depend on the output grid.
///
/// Run over a 2× step (1024² → 2048²) so the guard is affordable in the default suite; block A's
/// diagnostic (`hypsometry_work_attribution`) carries the same check over the 4× step to 8192²,
/// where it measured 865.4 → 865.5 m. Guarding the cheap step here and the expensive one there
/// is a cost decision, and it is stated rather than hidden: a regression that appears ONLY
/// between 2048² and 8192² would be caught by the diagnostic, not by this test.
#[test]
fn the_pre_incision_field_is_resolution_invariant() {
    let ss = SteinSteinParams::default();
    let (lo, hi) = (no_incision_terrain(1024), no_incision_terrain(COARSEST_TARGET));
    let (m_lo, p50_lo, p90_lo, land_lo) = stats(&lo, &ss);
    let (m_hi, p50_hi, p90_hi, land_hi) = stats(&hi, &ss);
    eprintln!(
        "[F1] pre-incision  1024²: mean {m_lo:.2} m p50 {p50_lo:.0} p90 {p90_lo:.0} land \
         {:.4} %\n[F1] pre-incision  2048²: mean {m_hi:.2} m p50 {p50_hi:.0} p90 {p90_hi:.0} \
         land {:.4} %",
        100.0 * land_lo,
        100.0 * land_hi
    );
    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(1e-9);
    // 1 % is generous against the 0.01 % measured over the 4× step, and is there to tolerate
    // the bicubic derivative converging, not to tolerate a mechanism appearing.
    assert!(
        rel(m_lo, m_hi) < 0.01,
        "pre-incision MEAN altitude is no longer resolution-invariant: {m_lo:.2} vs {m_hi:.2} m \
         ({:.2} %). Block A's conclusion — that the upscale is NOT the root cause of the \
         erosion-work divergence — rests on this. Re-open it rather than relaxing the bound.",
        100.0 * rel(m_lo, m_hi)
    );
    assert!(
        rel(land_lo, land_hi) < 0.01,
        "pre-incision EMERGED FRACTION is no longer resolution-invariant: {:.4} vs {:.4} %",
        100.0 * land_lo,
        100.0 * land_hi
    );
    for (name, a, b) in [("p50", p50_lo, p50_hi), ("p90", p90_lo, p90_hi)] {
        assert!(
            rel(a as f64, b as f64) < 0.02,
            "pre-incision altitude {name} is no longer resolution-invariant: {a:.0} vs {b:.0} m"
        );
    }
}

/// NEGATIVE CONTROL for F1 (rule 1): the comparison must be able to FAIL. Without it, a `stats`
/// that returned a constant, or two builds that silently produced the same grid, would pass F1
/// for the wrong reason — which is precisely how the first A2 table read "invariant" on a
/// quantity that could not vary.
#[test]
fn f1_can_detect_a_field_that_is_not_invariant() {
    let ss = SteinSteinParams::default();
    let base = no_incision_terrain(1024);
    // A field that differs by a KNOWN amount: every land cell raised by ~113 m (0.01 norm).
    let mut raised = base.clone();
    for x in raised.data.iter_mut() {
        if *x > SEA {
            *x += 0.01;
        }
    }
    let (m0, _, _, l0) = stats(&base, &ss);
    let (m1, _, _, l1) = stats(&raised, &ss);
    let rel = (m0 - m1).abs() / m0;
    eprintln!("[F1 control] mean {m0:.2} → {m1:.2} m ({:.2} %)", 100.0 * rel);
    assert!(rel > 0.01, "the F1 comparison cannot see a 113 m shift — it would pass on any field");
    assert!((l0 - l1).abs() < 1e-12, "the control must move the ALTITUDE, not the land mask");
}

/// F2 — `A_c` is SUB-CELL on the calibration grid, and this test **pins that known state**.
///
/// The first version was an `#[ignore]`d assertion of the property the code *ought* to have.
/// That guards nothing: an ignored test does not run, so it protects neither against a drift
/// nor against a silent correction. Inverted here — it asserts the CURRENT, measured,
/// documented state, so it **passes today and fails the day anyone touches `A_c`**, which is
/// when someone needs to read the finding.
///
/// The state being pinned: 2.6214 cells at 2048², under the declared 10-cell floor, and the
/// hillslope length it implies is `√0.1 km² = 316 m = 1.62 cells`. That is not a convention
/// question — a hillslope 1.6 cells wide cannot exist on the grid, by arithmetic.
#[test]
fn a_c_is_known_to_be_sub_cell_on_the_calibration_grid() {
    let cell_km2 = (DOMAIN_KM / COARSEST_TARGET as f32).powi(2);
    let cell_m = DOMAIN_KM / COARSEST_TARGET as f32 * 1000.0;
    let cells = RELIEF_V1_A_C_KM2 / cell_km2;
    let hillslope_len_m = RELIEF_V1_A_C_KM2.sqrt() * 1000.0;
    let hillslope_cells = hillslope_len_m / cell_m;
    eprintln!(
        "[F2] A_c = {RELIEF_V1_A_C_KM2} km² = {cells:.4} cells at {COARSEST_TARGET}² \
         ({cell_km2:.7} km²/cell); implied hillslope length {hillslope_len_m:.1} m = \
         {hillslope_cells:.2} cells"
    );
    assert!(
        cells < MIN_AREA_CELLS_FLOOR,
        "A_c is NO LONGER sub-cell at {COARSEST_TARGET}² ({cells:.4} cells, floor \
         {MIN_AREA_CELLS_FLOOR}). If that is deliberate, ADR Findings 57-58 and this test must \
         be updated together — the 'physically dimensioned but sub-cell' attribution, the \
         `S_ref` circularity and the B1 result (holding A_c in CELLS collapses the work ratio \
         from 3.18 to 1.23) all rest on this number."
    );
    assert!(
        hillslope_cells < 2.0,
        "the hillslope length implied by A_c is now {hillslope_cells:.2} cells; the claim that \
         the calibration grid cannot carry a hillslope regime needs re-deriving"
    );
}

/// The arithmetic behind F2, asserted so the numbers quoted in the ADR cannot drift from the
/// code. This one PASSES — it states what is, not what should be.
#[test]
fn the_channel_head_cell_count_is_recorded_correctly() {
    for (target, expect) in [(2048usize, 2.6214f32), (8192, 41.9430)] {
        let cell_km2 = (DOMAIN_KM / target as f32).powi(2);
        let cells = RELIEF_V1_A_C_KM2 / cell_km2;
        assert!(
            (cells - expect).abs() < 1e-3,
            "A_c cell count at {target}²: {cells:.4}, recorded {expect:.4}"
        );
    }
    // The ratio is the square of the resolution step — the whole mechanism in one line.
    let lo = RELIEF_V1_A_C_KM2 / (DOMAIN_KM / 2048.0f32).powi(2);
    let hi = RELIEF_V1_A_C_KM2 / (DOMAIN_KM / 8192.0f32).powi(2);
    assert!((hi / lo - 16.0).abs() < 1e-3, "the cell-count ratio must be (8192/2048)² = 16");
}
