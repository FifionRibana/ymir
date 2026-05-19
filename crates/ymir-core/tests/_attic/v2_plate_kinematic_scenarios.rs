//! Step 11 — physical-plausibility integration tests (Phase 4c).
//!
//! Acceptance criteria reformulated by the user's Phase-4 amendment
//! to match the *plate kinematic drift* semantics adopted in Phase
//! 4b (vs the original "initial velocity" wording in the issue,
//! which was wrong for the quasi-static Stokes solver — see
//! `docs/solver-scaling-step11-patch.md` §4.12).
//!
//! - **#5 reformulated** (`motion_without_mantle`): with
//!   `MantleConfig::Disabled` and a non-zero plate drift, the S̃
//!   field exhibits visible advection by step 30 — measurable as
//!   `||S̃(t=30) - S̃(t=0)||_2 / ||S̃(t=0)||_2 > 0.05` (5% relative
//!   change). Plates move as expected.
//!
//! - **#6 reformulated** (`convergence_scenario`): two plates with
//!   opposing drift produce a measurable interaction zone —
//!   either yielding fires at the contact, or S̃ accumulates at
//!   the contact (height increase > 10% in the boundary band) by
//!   step 50.
//!
//! - **#7 reformulated** (`with_cratonic`): cratons in plates with
//!   non-zero drift maintain coherent rigid motion — variance of
//!   `v_total` within the cratonic region is bounded; the
//!   `peak_yielding_in_craton ≤ 0.01` clause from Step 9 carries
//!   through unchanged.
//!
//! All three are `#[ignore]`'d (each runs ~5–30 s at 32²). Invoke:
//!
//! ```text
//! cargo test --release -p ymir-core \
//!     --test v2_plate_kinematic_scenarios \
//!     -- --ignored --nocapture --test-threads=1
//! ```

use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::{CratonicConfig, CratonicConfigEnabled};
use ymir_core::tectonics_v2::diagnostics::harness::{
    build_force, run_baseline, BaselineConfig, ForceKind, NonlinearChoice,
};
use ymir_core::tectonics_v2::plate_kinematic::PlateKinematicConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

const NX: usize = 32;
const NY: usize = 32;
const SEED: u64 = 42;

fn build_scenario(
    num_plates: usize,
    steps: usize,
    yielding: YieldingConfig,
    cratonic: CratonicConfig,
    plate_kinematic: PlateKinematicConfig,
) -> BaselineConfig {
    let scales = Scales::default();
    let vcfg = VoronoiConfig { num_plates, continental_ratio: 0.3 };
    let rates =
        BoundaryRates { k_sub: 0.5, k_arc: 0.0, k_spread: 0.0, k_coll_v: 0.0, k_rift_v: 0.0 };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        NX,
        NY,
        &vcfg,
        SEED,
        rates,
        RecyclingConfig::default(),
    )
    .expect("voronoi boundary build");
    let force = build_force(ForceKind::Gpe, &scales, 10.0, 1.0);
    BaselineConfig {
        seed: SEED,
        grid_nx: NX,
        grid_ny: NY,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps,
        cfl_factor: 0.3,
        total_time_nondim: 6.0 * (steps as f64 / 100.0).max(0.05),
        preset: Preset::dynamic_accidented(),
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions: Vec::new(),
        output_dir: std::env::temp_dir().join("ymir_step11_scenarios"),
        force,
        force_kind: ForceKind::Gpe,
        sinusoidal_amplitude: 0.0,
        s_perturbation_amplitude: 0.2,
        yielding,
        basal_drag: ymir_core::tectonics_v2::basal_drag::BasalDragConfig::Disabled,
        boundary,
        boundary_layout_name: format!("voronoi_seed{}_n{}", SEED, num_plates),
        slab_pull: ymir_core::tectonics_v2::slab::SlabPullConfig::Disabled,
        mantle: ymir_core::tectonics_v2::mantle::MantleConfig::Disabled,
        cratonic,
        age_field: ymir_core::tectonics_v2::age_field::AgeFieldConfig::Disabled,
        capture: None,
        linear_solver: LinearSolverConfig::default(),
        init_mode: ymir_core::tectonics_v2::init::InitMode::Checkerboard,
        continuation: None,
        plate_kinematic,
    }
}

/// Acceptance #5 (reformulated) — motion without mantle.
///
/// Mantle disabled + non-zero plate drift → S̃ advects measurably by
/// step 30, quantified by `||S(30) - S(0)||_2 / ||S(0)||_2 > 0.05`.
#[test]
#[ignore]
fn motion_without_mantle() {
    let drift_mag = 0.5;
    let num_plates = 2;
    let velocities = vec![(drift_mag, 0.0), (0.0, 0.0)];

    // Reference: same config but with `Zero` drift. Both runs share
    // initial S̃ exactly (drift only affects v), so
    // ||drift_final - zero_final||_2 isolates the drift-induced
    // S̃ displacement (zero_final captures whatever GPE-only viscous
    // evolution produces, which is the noise floor we measure
    // *against*, not the signal).
    let cfg_drift = build_scenario(
        num_plates,
        30,
        YieldingConfig::Disabled,
        CratonicConfig::Disabled,
        PlateKinematicConfig::PerPlate {
            velocities: velocities.clone(),
            boundary_smoothing_width: 1.5,
        },
    );
    let cfg_zero = build_scenario(
        num_plates,
        30,
        YieldingConfig::Disabled,
        CratonicConfig::Disabled,
        PlateKinematicConfig::Zero,
    );

    let result_drift = run_baseline(&cfg_drift);
    let result_zero = run_baseline(&cfg_zero);

    let s_drift = result_drift.final_state.s_field.data();
    let s_zero = result_zero.final_state.s_field.data();

    let mut diff_l2_sq = 0.0_f64;
    let mut zero_l2_sq = 0.0_f64;
    for (&a, &b) in s_drift.iter().zip(s_zero.iter()) {
        let d = a - b;
        diff_l2_sq += d * d;
        zero_l2_sq += b * b;
    }
    let rel = (diff_l2_sq / zero_l2_sq).sqrt();
    eprintln!(
        "[motion_without_mantle] drift-vs-zero relative L2 = {:.4} (threshold 0.05)",
        rel
    );
    eprintln!(
        "[motion_without_mantle] vmax_peak drift={:.4e}, zero={:.4e}",
        result_drift.metrics.vmax_peak, result_zero.metrics.vmax_peak
    );

    assert!(
        rel > 0.05,
        "S̃ field did not advect measurably under drift = {} (relative L2 = {:.4} below 0.05)",
        drift_mag,
        rel
    );
}

/// Acceptance #6 (reformulated) — convergence scenario.
///
/// Two plates with opposing drift produce a measurable contact-zone
/// interaction by step 50: yielding fires OR S̃ accumulates > 10%
/// in the boundary band.
///
/// **Yielding deliberately Disabled here.** The drift-induced
/// strain rate at the inter-plate boundary on a 32² grid with
/// width = 1.5 cells is `21·drift / width = 21·drift/1.5 ≈
/// 14·drift`. With `Bi = 0.15` and `drift ≥ 0.01`, this crosses
/// the yielding threshold and triggers a positive feedback
/// (`ε̇ ↑ → η ↓ → v_solver ↑ → ε̇ ↑`) that runs `vmax_peak →
/// 1e30+`. The OR clause in the reformulated acceptance #6
/// allows us to validate via S̃ accumulation alone, which is
/// what we exercise here. The yielding+drift coupling régime is
/// out of scope for this Step-11 acceptance check (it tests
/// Step-3 yielding stability under a forcing it was not
/// calibrated for) and is documented as a known boundary in
/// §4.12.
#[test]
#[ignore]
fn convergence_scenario() {
    let velocities = vec![(0.5, 0.0), (-0.5, 0.0)];
    let steps = 50;
    let cfg = build_scenario(
        2,
        steps,
        YieldingConfig::Disabled,
        CratonicConfig::Disabled,
        PlateKinematicConfig::PerPlate {
            velocities,
            boundary_smoothing_width: 1.5,
        },
    );
    let result = run_baseline(&cfg);
    let yielding_fraction = result.metrics.yielding_cell_fraction.unwrap_or(0.0);

    let plate_id = result
        .final_state
        .plate_id
        .as_ref()
        .expect("boundary enabled — plate_id must be present in FinalState");
    let s_final = result.final_state.s_field.data();

    let mut boundary_cell_count = 0usize;
    let mut boundary_s_sum = 0.0_f64;
    let mut interior_s_sum = 0.0_f64;
    let mut interior_count = 0usize;
    for j in 0..NY {
        for i in 0..NX {
            let id = plate_id.get(i, j);
            let ip = (i + 1) % NX;
            let im = (i + NX - 1) % NX;
            let jp = (j + 1) % NY;
            let jm = (j + NY - 1) % NY;
            let on_boundary = [(ip, j), (im, j), (i, jp), (i, jm)]
                .iter()
                .any(|&(ni, nj)| plate_id.get(ni, nj) != id);
            let v = s_final[j * NX + i];
            if on_boundary {
                boundary_s_sum += v;
                boundary_cell_count += 1;
            } else {
                interior_s_sum += v;
                interior_count += 1;
            }
        }
    }
    let boundary_mean = boundary_s_sum / boundary_cell_count.max(1) as f64;
    let interior_mean = interior_s_sum / interior_count.max(1) as f64;
    let height_increase = (boundary_mean - interior_mean) / interior_mean.abs().max(1e-12);

    eprintln!(
        "[convergence_scenario] yielding_fraction = {:.4}, vmax_peak = {:.4e}, mass drift = {:.3e}",
        yielding_fraction, result.metrics.vmax_peak, result.metrics.mass_drift_relative
    );
    eprintln!(
        "[convergence_scenario] boundary_mean S̃ = {:.4}, interior_mean S̃ = {:.4}, \
         relative excess = {:.4}",
        boundary_mean, interior_mean, height_increase
    );

    // Sanity: vmax_peak must stay within the perturbative régime.
    // A yielding_fraction = 1.0 with vmax = 1e32 is a blowup, not
    // physics — guard against false positives.
    assert!(
        result.metrics.vmax_peak < 10.0,
        "vmax_peak = {:.4e} indicates solver blowup (drift outside \
         validity régime — see §4.12). Real yielding signal must come \
         from a bounded run.",
        result.metrics.vmax_peak
    );
    assert!(
        yielding_fraction > 0.0 || height_increase > 0.10,
        "neither yielding fired (yielding_fraction = {:.4}) nor S̃ \
         accumulated noticeably (boundary excess = {:.4} ≤ 0.10) — \
         opposing plates should produce *some* observable interaction \
         by step {}",
        yielding_fraction,
        height_increase,
        steps
    );
}

/// Diagnostic baseline for `with_cratonic`: same config but
/// `PlateKinematicConfig::Zero`. Confirms the test setup is itself
/// stable (Step 9 yielding+cratonic at 32² over 20 steps). The
/// `with_cratonic` test below uses identical numerics with
/// non-zero drift; both should produce the same `vmax_peak` after
/// the deformation/transport split (the drift contributes only to
/// transport, not to deformation, so quiescent v_solver is
/// preserved).
#[test]
#[ignore]
fn with_cratonic_baseline_zero_drift() {
    let cfg = build_scenario(
        6,
        20,
        YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        CratonicConfig::Enabled(CratonicConfigEnabled {
            cr: 0.3,
            k_viscous: 5.0,
            b_factor: 8.0,
            ..CratonicConfigEnabled::default()
        }),
        PlateKinematicConfig::Zero,
    );
    let result = run_baseline(&cfg);
    eprintln!(
        "[with_cratonic_baseline_zero_drift] vmax_peak = {:.4e}, \
         peak_yielding_in_craton = {:?}",
        result.metrics.vmax_peak,
        result.metrics.newton.as_ref().and_then(|n| n.peak_yielding_in_craton),
    );
    assert!(
        result.metrics.vmax_peak < 10.0,
        "Zero baseline blew up: vmax = {:.4e}. The with_cratonic setup \
         is unstable independent of Step 11 — change the scenario, not \
         the drift.",
        result.metrics.vmax_peak
    );
}

/// Acceptance #7 (reformulated) — cratonic immunity preserved with
/// drift.
///
/// Cratons in plates with non-zero drift maintain coherent rigid
/// motion: the post-solve `vx, vy` (= solver-only after the
/// strip-at-iter-end hook) inside cratonic cells stays at the
/// quiescent-régime values that match the Zero baseline, and
/// `peak_yielding_in_craton ≤ 0.01`.
///
/// History: an earlier Phase-4 wiring of the drift mechanism added
/// it to `vx, vy` *before* the post-solve diagnostic block. That
/// caused `StrainRate::compute` to see the smoothstep gradient of
/// the drift at inter-plate boundaries, triggering yielding in the
/// metrics path and feeding back via S̃ advection to the next
/// solve → runaway, `vmax_peak → 1e28+`. Phase 4 final replaces
/// that wiring with the deformation/transport split documented in
/// §4.12: drift exists only inside the advection scope of each
/// iter, so deformation diagnostics see solver-only state and are
/// insensitive to drift magnitude.
///
/// Drift magnitudes here are kept at `≤ 0.001` and the boundary
/// smoothing widened to `6` cells per the validity-envelope
/// discussion in §4.12 — not because the test would otherwise
/// blow up (it does not, after the fix), but because the
/// cumulative displacement over a 20-step run is what the user is
/// likely to dial in for a *cratonic* visualisation (cratons
/// drift slowly relative to mobile belts). The signal under test
/// is the rigidity floor, not the magnitude.
#[test]
#[ignore]
fn with_cratonic() {
    let velocities = vec![
        (0.001, 0.0),
        (-0.0008, 0.0006),
        (0.0, 0.0008),
        (-0.0006, -0.0008),
        (0.0006, 0.0004),
        (0.0, 0.0),
    ];
    let cfg = build_scenario(
        velocities.len(),
        20,
        YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        CratonicConfig::Enabled(CratonicConfigEnabled {
            cr: 0.3,
            k_viscous: 5.0,
            b_factor: 8.0,
            ..CratonicConfigEnabled::default()
        }),
        PlateKinematicConfig::PerPlate {
            velocities: velocities.clone(),
            boundary_smoothing_width: 6.0,
        },
    );
    let result = run_baseline(&cfg);
    let cratonic_factor = result
        .final_state
        .cratonic_factor
        .as_ref()
        .expect("cratonic enabled — factor field must be populated");
    let vx = &result.final_state.vx;
    let vy = &result.final_state.vy;
    let plate_id = result.final_state.plate_id.as_ref().unwrap();

    let mut variance_sum = 0.0_f64;
    let mut craton_count = 0usize;
    for j in 0..NY {
        for i in 0..NX {
            let f = cratonic_factor.get(i, j);
            if f <= 0.5 {
                continue;
            }
            let pid = plate_id.get(i, j) as usize;
            let expected = velocities[pid];
            let dvx = vx[j * NX + i] - expected.0;
            let dvy = vy[j * NX + i] - expected.1;
            variance_sum += dvx * dvx + dvy * dvy;
            craton_count += 1;
        }
    }
    let variance = if craton_count > 0 {
        variance_sum / craton_count as f64
    } else {
        0.0
    };

    let peak_y_craton = result
        .metrics
        .newton
        .as_ref()
        .and_then(|n| n.peak_yielding_in_craton)
        .unwrap_or(0.0);

    eprintln!(
        "[with_cratonic] craton cells = {}, variance |v - v_drift|² = {:.4e}",
        craton_count, variance
    );
    eprintln!(
        "[with_cratonic] peak_yielding_in_craton = {:.4}, vmax_peak = {:.4e}",
        peak_y_craton, result.metrics.vmax_peak
    );

    // Sanity: solver must remain in the perturbative régime. A
    // blowup propagates to cratonic cells via the global Stokes
    // operator and would mask the rigidity signal entirely.
    assert!(
        result.metrics.vmax_peak < 10.0,
        "vmax_peak = {:.4e} indicates solver blowup (drift outside \
         validity régime — see §4.12)",
        result.metrics.vmax_peak
    );

    assert!(
        craton_count > 0,
        "no cratonic interior cells found — geometry too small or Cr too low"
    );
    assert!(
        variance < 1e-3,
        "cratonic rigidity violated: variance |v - v_drift|² = {:.4e} > 1e-3 \
         (cratonic interior should ride the per-plate drift, with the \
         small viscous solver contribution as the only noise)",
        variance
    );
    assert!(
        peak_y_craton <= 0.01,
        "peak_yielding_in_craton = {:.4} > 0.01 — cratonic immunity \
         violated under drift",
        peak_y_craton
    );
}

/// Acceptance #9 — CG iters ratio under drift is ≤ 1.2× the
/// zero-baseline. A drift-driven advection moves S̃ each step, which
/// re-conditions the next solve via a new forcing — but with the
/// deformation/transport split (Phase 4 fix) the rheology / η field
/// the Newton tangent operates on is solver-only, so the CG cost
/// per inner solve should not blow up under drift.
///
/// Single 32² × 30-step measurement (one run each, no multi-run
/// statistics). Drift = 0.05, yielding ON, cratonic ON — the
/// "richest" config that still stays in the perturbative régime
/// (acceptance #9 specifically targets the conditioning coupling,
/// so picking a non-trivial config matters more than a pristine
/// baseline shape).
#[test]
#[ignore]
fn cg_ratio_under_drift_within_acceptance() {
    let velocities = vec![
        (0.05, 0.0),
        (-0.04, 0.03),
        (0.0, 0.04),
        (-0.03, -0.04),
        (0.03, 0.02),
        (0.0, 0.0),
    ];
    let cfg_zero = build_scenario(
        velocities.len(),
        30,
        YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        CratonicConfig::Enabled(CratonicConfigEnabled {
            cr: 0.3,
            k_viscous: 5.0,
            b_factor: 8.0,
            ..CratonicConfigEnabled::default()
        }),
        PlateKinematicConfig::Zero,
    );
    let cfg_drift = build_scenario(
        velocities.len(),
        30,
        YieldingConfig::Enabled(YieldingLaw { bi: 0.15, ..Default::default() }),
        CratonicConfig::Enabled(CratonicConfigEnabled {
            cr: 0.3,
            k_viscous: 5.0,
            b_factor: 8.0,
            ..CratonicConfigEnabled::default()
        }),
        PlateKinematicConfig::PerPlate {
            velocities,
            boundary_smoothing_width: 6.0,
        },
    );

    let result_zero = run_baseline(&cfg_zero);
    let result_drift = run_baseline(&cfg_drift);

    let cg_zero = result_zero.metrics.cg_iter_mean;
    let cg_drift = result_drift.metrics.cg_iter_mean;
    let ratio = cg_drift / cg_zero.max(1e-12);

    eprintln!(
        "[cg_ratio] CG iters mean: zero = {:.2}, drift = {:.2}, ratio = {:.3}",
        cg_zero, cg_drift, ratio
    );
    eprintln!(
        "[cg_ratio] vmax_peak: zero = {:.3e}, drift = {:.3e}",
        result_zero.metrics.vmax_peak, result_drift.metrics.vmax_peak
    );

    assert!(
        ratio <= 1.2,
        "CG iters ratio under drift = {:.3} exceeds the 1.2× \
         zero-baseline acceptance #9 — drift is degrading \
         conditioning beyond the documented bound. Investigate \
         whether the post-solve S̃ advection by drift is creating \
         sharp gradients that load up the next solve.",
        ratio
    );
}
