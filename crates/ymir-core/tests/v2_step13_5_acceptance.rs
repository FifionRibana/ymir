//! Step 13.5 Phase 4 — empirical calibration sweep for the oceanic
//! FBM extension. Mirror of Step 13 Phase 6's
//! `fbm_calibration_probe` (continental side); produces the table
//! that grounds the Phase 5 default value of
//! `fbm_amplitude_oceanic`.
//!
//! Heavy `#[ignore]` test — invoked explicitly:
//!
//! ```text
//! cargo test --release -p ymir-core --test v2_step13_5_acceptance \
//!     -- --ignored --nocapture
//! ```
//!
//! Voronoï layout: `single_continent` (seed=12, 4 plates, 50 %
//! continental) at 64², same as Step 13's calibration so the
//! oceanic and continental probes are directly comparable.
//!
//! Per `(fbm_amplitude_oceanic, fbm_scale_oceanic)` combination
//! we measure:
//!
//! - **σ_fbm_oceanic_isolated** — std-dev of
//!   `S̃_with_oceanic_FBM − S̃_without_oceanic_FBM` over oceanic
//!   cells. Direct measurement of the FBM contribution, isolated
//!   from the (uniform) oceanic baseline. Same metric reformulation
//!   that Step 13 Phase 6 introduced for acceptance #7 (vacuous-
//!   truth guard).
//! - **max(S̃_oceanic) / min(S̃_oceanic)** — confirms the
//!   `OCEANIC_CLAMP_MAX = 0.49` upper bound is honoured (D7) and
//!   the lower clamp at 0.0 is never reached at sane amplitudes.
//! - **clipping fraction** — share of oceanic cells whose value
//!   landed at the upper clamp. Above ~10 % clipping, the
//!   distribution is saturated and the (amp, scale) pair is no
//!   longer in the linear-σ regime.

use ymir_core::tectonics_v2::boundaries::PlateType;
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::init::{
    FBM_AMPLITUDE_DEFAULT, FBM_AMPLITUDE_OCEANIC_DEFAULT, FBM_LACUNARITY_DEFAULT,
    FBM_OCTAVES_DEFAULT, FBM_PERSISTENCE_DEFAULT, FBM_SCALE_DEFAULT, FBM_SEED_DEFAULT, InitContext,
    InitMode, OCEANIC_CLAMP_MAX, PlateInitData, ProfileShape, init_s_field,
};
use ymir_core::tectonics_v2::voronoi::{VoronoiConfig, VoronoiPlates, generate_voronoi};

const SC_SEED: u64 = 12;
const SC_NUM_PLATES: usize = 4;
const SC_CONTINENTAL_RATIO: f64 = 0.5;
const CONTINENTAL_VALUE: f64 = 0.95;
const OCEANIC_VALUE: f64 = 0.20;

fn build_voronoi(nx: usize, ny: usize) -> VoronoiPlates {
    generate_voronoi(
        nx,
        ny,
        &VoronoiConfig { num_plates: SC_NUM_PLATES, continental_ratio: SC_CONTINENTAL_RATIO },
        SC_SEED,
    )
}

fn make_ctx<'a>(plates: &'a VoronoiPlates, nx: usize, ny: usize) -> InitContext<'a> {
    InitContext {
        nx,
        ny,
        seed: SC_SEED,
        amplitude: 0.0,
        plate_data: Some(PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        }),
    }
}

/// Build a `RadialProfileWithFBM` with the Step 13 continental
/// defaults; only the oceanic block changes per probe call.
fn build_with_oceanic(amp_oceanic: f64, scale_oceanic: f64) -> InitMode {
    InitMode::RadialProfileWithFBM {
        continental_value: CONTINENTAL_VALUE,
        oceanic_value: OCEANIC_VALUE,
        profile_shape: ProfileShape::Smoothstep,
        fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
        fbm_octaves: FBM_OCTAVES_DEFAULT,
        fbm_persistence: FBM_PERSISTENCE_DEFAULT,
        fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
        fbm_scale: FBM_SCALE_DEFAULT,
        fbm_seed: FBM_SEED_DEFAULT,
        apply_fbm_to_oceanic: true,
        fbm_amplitude_oceanic: amp_oceanic,
        fbm_scale_oceanic: Some(scale_oceanic),
        fbm_seed_oceanic: None, // derive via XOR — exercise the default-derivation path
    }
}

/// Per-cell stats over oceanic cells: (count, σ_diff, max(s),
/// min(s), clipping_fraction). `σ_diff` is std-dev of
/// `s_enabled − s_baseline` where the baseline has oceanic FBM
/// disabled (uniform `oceanic_value`).
fn oceanic_stats(
    plates: &VoronoiPlates,
    s_baseline: &Field2D,
    s_enabled: &Field2D,
    nx: usize,
    ny: usize,
) -> (usize, f64, f64, f64, f64) {
    let mut diffs: Vec<f64> = Vec::new();
    let mut max_s = f64::NEG_INFINITY;
    let mut min_s = f64::INFINITY;
    let mut clipped = 0usize;
    for j in 0..ny {
        for i in 0..nx {
            if !matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                continue;
            }
            let v = s_enabled.get(i, j);
            let b = s_baseline.get(i, j);
            diffs.push(v - b);
            if v > max_s {
                max_s = v;
            }
            if v < min_s {
                min_s = v;
            }
            if (v - OCEANIC_CLAMP_MAX).abs() < 1e-12 {
                clipped += 1;
            }
        }
    }
    let n = diffs.len();
    let mean = diffs.iter().sum::<f64>() / n.max(1) as f64;
    let var = diffs.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / n.max(1) as f64;
    let sigma = var.sqrt();
    let clip_frac = clipped as f64 / n.max(1) as f64;
    (n, sigma, max_s, min_s, clip_frac)
}

#[test]
#[ignore]
fn fbm_oceanic_calibration_probe() {
    let nx = 64;
    let ny = 64;
    let plates = build_voronoi(nx, ny);
    let ctx = make_ctx(&plates, nx, ny);

    // Baseline build with oceanic FBM disabled — every oceanic
    // cell sits at OCEANIC_VALUE = 0.20 by construction.
    let baseline_mode = InitMode::RadialProfileWithFBM {
        continental_value: CONTINENTAL_VALUE,
        oceanic_value: OCEANIC_VALUE,
        profile_shape: ProfileShape::Smoothstep,
        fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
        fbm_octaves: FBM_OCTAVES_DEFAULT,
        fbm_persistence: FBM_PERSISTENCE_DEFAULT,
        fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
        fbm_scale: FBM_SCALE_DEFAULT,
        fbm_seed: FBM_SEED_DEFAULT,
        apply_fbm_to_oceanic: false,
        fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
        fbm_scale_oceanic: None,
        fbm_seed_oceanic: None,
    };
    let s_baseline = init_s_field(baseline_mode, &ctx);

    eprintln!(
        "Step 13.5 Phase 4 — fbm_oceanic_calibration_probe on \
         single_continent (64², seed=12, {} plates, {:.0}% continental)",
        SC_NUM_PLATES,
        SC_CONTINENTAL_RATIO * 100.0
    );
    eprintln!("Per (amplitude, scale) cell: σ_fbm_oceanic_isolated, max(S̃), clip%.");
    eprintln!("Target: σ ∈ [0.02, 0.08], max ≤ 0.49, clip% small.");
    eprintln!();
    eprintln!("amp \\ scale  |    0.05         0.10         0.15         0.20");

    for amp in [0.05_f64, 0.10, 0.15, 0.20, 0.25] {
        eprint!("amp={:.2}      |", amp);
        for scale in [0.05_f64, 0.10, 0.15, 0.20] {
            let s = init_s_field(build_with_oceanic(amp, scale), &ctx);
            let (_, sigma, max_s, _min_s, clip_frac) =
                oceanic_stats(&plates, &s_baseline, &s, nx, ny);
            eprint!("  σ={:.4}/{:.3}/{:>3.0}%", sigma, max_s, clip_frac * 100.0);
        }
        eprintln!();
    }
    eprintln!();
    eprintln!(
        "Reading: σ = std(S_with_oceanic − S_without_oceanic) over oceanic cells; \
         max = max(S̃_oceanic); clip% = oceanic cells saturated at {}.",
        OCEANIC_CLAMP_MAX
    );
}

/// Step 13.5 acceptance #7 — oceanic FBM contribution measurable.
/// With the Phase 5 calibrated defaults
/// (`FBM_AMPLITUDE_OCEANIC_DEFAULT = 0.15`,
/// `fbm_scale_oceanic = None` ⇒ reuse `fbm_scale = 0.10`,
/// `fbm_seed_oceanic = None` ⇒ XOR derive), the FBM-isolated
/// oceanic std-dev `σ_fbm_oceanic_isolated` lands in the issue's
/// target band `[0.02, 0.08]` on `single_continent` 64².
///
/// Reformulated as Step 13's acceptance #7 — direct measurement
/// of the FBM contribution (subtract the FBM-disabled baseline)
/// rather than total `σ(S̃_oceanic)`. With the Step 13.5 default
/// oceanic baseline being uniform at `oceanic_value`, the two
/// metrics agree numerically here, but the FBM-isolated form
/// is robust to any future change of the oceanic baseline (e.g.,
/// per-plate variation in Step 14+).
fn run_oceanic_fbm_amplitude_target(nx: usize, ny: usize, label: &str) {
    let plates = build_voronoi(nx, ny);
    let ctx = make_ctx(&plates, nx, ny);

    // FBM-disabled baseline (Step 13 oceanic uniform).
    let baseline_mode = InitMode::RadialProfileWithFBM {
        continental_value: CONTINENTAL_VALUE,
        oceanic_value: OCEANIC_VALUE,
        profile_shape: ProfileShape::Smoothstep,
        fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
        fbm_octaves: FBM_OCTAVES_DEFAULT,
        fbm_persistence: FBM_PERSISTENCE_DEFAULT,
        fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
        fbm_scale: FBM_SCALE_DEFAULT,
        fbm_seed: FBM_SEED_DEFAULT,
        apply_fbm_to_oceanic: false,
        fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
        fbm_scale_oceanic: None,
        fbm_seed_oceanic: None,
    };
    let s_baseline = init_s_field(baseline_mode, &ctx);

    // FBM-enabled with the Phase 5 calibrated defaults.
    let enabled_mode = InitMode::RadialProfileWithFBM {
        continental_value: CONTINENTAL_VALUE,
        oceanic_value: OCEANIC_VALUE,
        profile_shape: ProfileShape::Smoothstep,
        fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
        fbm_octaves: FBM_OCTAVES_DEFAULT,
        fbm_persistence: FBM_PERSISTENCE_DEFAULT,
        fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
        fbm_scale: FBM_SCALE_DEFAULT,
        fbm_seed: FBM_SEED_DEFAULT,
        apply_fbm_to_oceanic: true,
        fbm_amplitude_oceanic: FBM_AMPLITUDE_OCEANIC_DEFAULT,
        fbm_scale_oceanic: None,
        fbm_seed_oceanic: None,
    };
    let s_enabled = init_s_field(enabled_mode, &ctx);

    let (count, sigma, max_s, _min_s, clip_frac) =
        oceanic_stats(&plates, &s_baseline, &s_enabled, nx, ny);

    eprintln!(
        "Step 13.5 acceptance #7 — oceanic_fbm_amplitude_target on \
         single_continent {} (seed=12, {} plates, {:.0}% continental):",
        label,
        SC_NUM_PLATES,
        SC_CONTINENTAL_RATIO * 100.0
    );
    eprintln!("  oceanic cells           : {}", count);
    eprintln!("  σ_fbm_oceanic_isolated  : {:.4} (target [0.02, 0.08])", sigma);
    eprintln!("  max(S̃_oceanic)          : {:.4} (≤ {} required)", max_s, OCEANIC_CLAMP_MAX);
    eprintln!("  clip fraction           : {:.2}% (≤ 0.49 clamp)", clip_frac * 100.0);

    assert!(
        count > 0,
        "[{}] no oceanic cells found — single_continent setup invariant violated",
        label
    );
    assert!(
        sigma >= 0.02,
        "[{}] σ_fbm_oceanic_isolated = {:.4} below the 0.02 lower bound — \
         the FBM contribution is too small to be visible. \
         Default amplitude {} may need raising.",
        label,
        sigma,
        FBM_AMPLITUDE_OCEANIC_DEFAULT,
    );
    assert!(
        sigma <= 0.08,
        "[{}] σ_fbm_oceanic_isolated = {:.4} above the 0.08 upper bound — \
         the FBM contribution is approaching the threshold. \
         Default amplitude {} may need lowering.",
        label,
        sigma,
        FBM_AMPLITUDE_OCEANIC_DEFAULT,
    );
    assert!(
        max_s <= OCEANIC_CLAMP_MAX + 1e-15,
        "[{}] max(S̃_oceanic) = {:.6} > OCEANIC_CLAMP_MAX = {} — \
         threshold protection violated, volcanic islands would emerge \
         (out of scope for Step 13.5).",
        label,
        max_s,
        OCEANIC_CLAMP_MAX,
    );
    assert!(
        clip_frac < 0.05,
        "[{}] clip fraction = {:.2}% — sustained saturation at the {} clamp \
         indicates the amplitude default is too aggressive for this grid. \
         Mechanism is healthy (no threshold crossing) but the perturbation \
         is shape-distorted.",
        label,
        clip_frac * 100.0,
        OCEANIC_CLAMP_MAX,
    );
}

#[test]
fn oceanic_fbm_amplitude_target_64sq() {
    run_oceanic_fbm_amplitude_target(64, 64, "64²");
}

/// Phase 7 — same acceptance at 32² (milestone validation grid pair).
/// Single_continent at 32² has ~512 oceanic cells (vs ~928 at 64²),
/// well above the 150-cell threshold introduced in Step 13's
/// "small-sample noise floor" caveat. The σ measurement is
/// statistically reliable; the assertion stays strict (no
/// largest-plate-only relaxation here, because we measure aggregate
/// `σ` over all oceanic cells, not per-plate).
#[test]
fn oceanic_fbm_amplitude_target_32sq() {
    run_oceanic_fbm_amplitude_target(32, 32, "32²");
}
