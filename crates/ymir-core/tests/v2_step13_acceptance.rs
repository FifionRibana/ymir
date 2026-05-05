//! Step 13 Phase 6 — physical-plausibility acceptance tests for
//! the new init modes.
//!
//! Two checks, run on the [`single_continent`] preset Voronoï layout
//! (seed=12, 4 plates, 50 % continental, 64²) chosen for its few
//! large continental plates — small plates would not give the
//! gradient enough cells to be measurable.
//!
//! - **Acceptance #6 — margins gradient visible**:
//!   `RadialProfile` with `Smoothstep` and the spec defaults
//!   (`continental=0.95`, `oceanic=0.20`) must produce at least 2
//!   continental cells with `S̃ ∈ [0.5, 0.7]`. These are the
//!   "intermediate" cells that visually witness a smooth gradient
//!   between the boundary floor (`S̃ → 0.20`) and the interior peak
//!   (`S̃ → 0.95`).
//!
//! - **Acceptance #7 — intra-plate heterogeneity (Phase 6
//!   amendment)**: `RadialProfileWithFBM` with the spec defaults
//!   (amplitude=0.10, scale=0.10) must produce
//!   `σ_fbm_isolated ≥ 0.040` over the *interior cells* (`t = d /
//!   L_plate > 0.5`) of **every continental plate**.
//!   `σ_fbm_isolated = std_dev(S_with_fbm − S_radial_only)` —
//!   directly measures the FBM contribution rather than the
//!   total `σ(S̃)` (which is dominated by the radial gradient and
//!   passes vacuously). Calibration of `fbm_scale = 0.10`
//!   discussed in `radial_profile_fbm` module docstring and in
//!   `docs/solver-scaling-step13-patch.md` §4.13.
//!
//! No `#[ignore]` — these tests are cheap (init-only, 64²) and
//! belong in the regular acceptance suite.

use ymir_core::tectonics_v2::boundaries::PlateType;
use ymir_core::tectonics_v2::init::{
    init_s_field, FBM_AMPLITUDE_DEFAULT, FBM_LACUNARITY_DEFAULT, FBM_OCTAVES_DEFAULT,
    FBM_PERSISTENCE_DEFAULT, FBM_SCALE_DEFAULT, FBM_SEED_DEFAULT, InitContext, InitMode,
    PlateInitData, ProfileShape,
};
use ymir_core::tectonics_v2::voronoi::{
    VoronoiConfig, VoronoiPlates, compute_dist_to_inter_plate_boundary, generate_voronoi,
};

/// `single_continent.json` preset parameters — large continental
/// plates suitable for clean gradient measurement.
const NX: usize = 64;
const NY: usize = 64;
const SEED: u64 = 12;
const NUM_PLATES: usize = 4;
const CONTINENTAL_RATIO: f64 = 0.5;
const CONTINENTAL_VALUE: f64 = 0.95;
const OCEANIC_VALUE: f64 = 0.20;

fn build_voronoi() -> VoronoiPlates {
    generate_voronoi(
        NX,
        NY,
        &VoronoiConfig {
            num_plates: NUM_PLATES,
            continental_ratio: CONTINENTAL_RATIO,
        },
        SEED,
    )
}

fn make_ctx(plates: &VoronoiPlates) -> InitContext<'_> {
    InitContext {
        nx: NX,
        ny: NY,
        seed: SEED,
        amplitude: 0.0,
        plate_data: Some(PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        }),
    }
}

/// Per-plate L_plate = max BFS distance over cells of that plate.
/// Cells with non-finite distance count as `INFINITY` and are
/// skipped here (single-plate-on-torus degenerate case — saturates
/// to 1.0 in radial_profile, but isn't relevant for these tests).
fn per_plate_l(plates: &VoronoiPlates, dist: &ymir_core::tectonics_v2::field::Field2D) -> Vec<f64> {
    let mut max_d: Vec<f64> = vec![0.0; plates.num_plates];
    for j in 0..NY {
        for i in 0..NX {
            let pid = plates.plate_id.get(i, j) as usize;
            let d = dist.get(i, j);
            if d.is_finite() && d > max_d[pid] {
                max_d[pid] = d;
            }
        }
    }
    max_d
}

/// Acceptance #6 — margins gradient: with the default
/// `RadialProfile { Smoothstep, 0.95, 0.20 }` on `single_continent`,
/// at least 2 continental cells lie in the intermediate-value band
/// `[0.5, 0.7]`.
#[test]
fn radial_profile_margins_gradient_visible() {
    let plates = build_voronoi();
    let ctx = make_ctx(&plates);
    let s = init_s_field(
        InitMode::RadialProfile {
            continental_value: CONTINENTAL_VALUE,
            oceanic_value: OCEANIC_VALUE,
            profile_shape: ProfileShape::Smoothstep,
        },
        &ctx,
    );

    let (lo, hi) = (0.5_f64, 0.7_f64);
    let mut intermediate_count = 0usize;
    let mut continental_count = 0usize;
    let mut min_in_band = f64::INFINITY;
    let mut max_in_band = f64::NEG_INFINITY;
    for j in 0..NY {
        for i in 0..NX {
            if !matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                continue;
            }
            continental_count += 1;
            let v = s.get(i, j);
            if v >= lo && v <= hi {
                intermediate_count += 1;
                if v < min_in_band {
                    min_in_band = v;
                }
                if v > max_in_band {
                    max_in_band = v;
                }
            }
        }
    }

    eprintln!(
        "Acceptance #6 — RadialProfile{{Smoothstep, 0.95, 0.20}} on single_continent (64², seed=12, 4 plates, 50% continental):"
    );
    eprintln!(
        "  continental cells   : {}, intermediate [0.5, 0.7] : {} ({:.1} %)",
        continental_count,
        intermediate_count,
        100.0 * intermediate_count as f64 / continental_count.max(1) as f64
    );
    if intermediate_count > 0 {
        eprintln!(
            "  min/max value of intermediate cells : {:.4} / {:.4}",
            min_in_band, max_in_band
        );
    }

    assert!(
        intermediate_count >= 2,
        "expected ≥ 2 continental cells with S̃ ∈ [0.5, 0.7] for the gradient \
         zone to be visible; got {}",
        intermediate_count
    );
}

/// Acceptance #7 (Phase 6 amendment) — intra-plate heterogeneity:
/// with the spec defaults (`fbm_amplitude=0.10, fbm_scale=0.10`)
/// on `single_continent`, the **FBM-isolated** std-dev
/// `σ(S̃_FBM − S̃_radial)` over interior cells (`t > 0.5`) is
/// `≥ 0.040` for every continental plate.
///
/// Lower bound only — there is no upper-bound assertion. The Phase
/// 5 UI clamps `fbm_amplitude` to `[0.0, 0.40]` which already
/// caps how much heterogeneity FBM can introduce; the algorithm
/// also clamps S̃ to `[0, 1]`. The acceptance is therefore "FBM
/// must actually contribute ≥ 0.040 std-dev to the interior" —
/// a vacuous-truth guard against the original
/// `σ_total ∈ [0.04, 0.10]` formulation, which was satisfied by
/// the radial gradient alone (FBM contribution ≈ 0.018 with the
/// pre-amendment `fbm_scale = 0.25`).
///
/// The test also reports `σ_radial` and `σ_total` for context.
#[test]
fn radial_profile_with_fbm_intra_plate_heterogeneity() {
    let plates = build_voronoi();
    let ctx = make_ctx(&plates);

    let s_fbm = init_s_field(
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
        },
        &ctx,
    );
    let s_radial = init_s_field(
        InitMode::RadialProfile {
            continental_value: CONTINENTAL_VALUE,
            oceanic_value: OCEANIC_VALUE,
            profile_shape: ProfileShape::Smoothstep,
        },
        &ctx,
    );

    let bfs = compute_dist_to_inter_plate_boundary(NX, NY, &plates.plate_id);
    let l_plate = per_plate_l(&plates, &bfs.distance);

    // Per-plate metrics. The largest continental plate drives
    // acceptance #7; the others are reported for context.
    eprintln!(
        "Acceptance #7 — RadialProfileWithFBM{{amplitude={}, …}} on single_continent (64², seed=12):",
        FBM_AMPLITUDE_DEFAULT
    );
    eprintln!(
        "  Per-plate breakdown (interior = continental ∧ t > 0.5):"
    );
    eprintln!(
        "  pid type           cells_int  σ_radial   σ_total    σ_fbm_iso  L_plate"
    );

    let mut continental_metrics: Vec<(usize, usize, f64, f64, f64, f64)> = Vec::new();
    for pid in 0..plates.num_plates {
        let pt = plates.per_plate_type[pid];
        let lp = l_plate[pid];
        if lp <= 0.0 {
            continue;
        }
        let (cells_int, sigma_radial, sigma_total, sigma_fbm_iso) =
            interior_stats(pid, &plates, &bfs.distance, lp, &s_radial, &s_fbm);
        eprintln!(
            "  {:>3} {:<14} {:>9}  {:>8.4}   {:>8.4}   {:>8.4}   {:>5.2}",
            pid,
            format!("{:?}", pt),
            cells_int,
            sigma_radial,
            sigma_total,
            sigma_fbm_iso,
            lp
        );
        if matches!(pt, PlateType::Continental) && cells_int > 0 {
            continental_metrics.push((pid, cells_int, sigma_radial, sigma_total, sigma_fbm_iso, lp));
        }
    }
    assert!(
        !continental_metrics.is_empty(),
        "no continental plate with non-degenerate interior cells found — \
         the single_continent test setup needs revisiting"
    );

    eprintln!(
        "  Acceptance: σ_fbm_isolated ≥ 0.040 required for every continental plate."
    );

    const SIGMA_FBM_LOWER_BOUND: f64 = 0.040;
    let mut failures: Vec<String> = Vec::new();
    for &(pid, cells_int, sigma_radial, sigma_total, sigma_fbm_iso, lp) in &continental_metrics {
        if sigma_fbm_iso < SIGMA_FBM_LOWER_BOUND {
            failures.push(format!(
                "pid={} σ_fbm_iso={:.4} < 0.040 (cells_int={}, L_plate={:.2}, \
                 σ_radial={:.4}, σ_total={:.4})",
                pid, sigma_fbm_iso, cells_int, lp, sigma_radial, sigma_total
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "σ_fbm_isolated below the 0.040 lower bound on the following continental \
         plates — FBM is not introducing measurable intra-plate variability with \
         current defaults (fbm_amplitude=0.10, fbm_scale=0.10):\n  {}",
        failures.join("\n  ")
    );
}

/// Diagnostic probe — sweep `(fbm_scale, fbm_amplitude)` and
/// report σ_fbm_isolated on the largest continental plate so the
/// reviewer can pick a calibration empirically. `#[ignore]` so it
/// only runs on demand:
///
/// ```text
/// cargo test --release -p ymir-core --test v2_step13_acceptance \
///     fbm_calibration_probe -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn fbm_calibration_probe() {
    let plates = build_voronoi();
    let ctx = make_ctx(&plates);
    let bfs = compute_dist_to_inter_plate_boundary(NX, NY, &plates.plate_id);
    let l_plate = per_plate_l(&plates, &bfs.distance);

    let s_radial = init_s_field(
        InitMode::RadialProfile {
            continental_value: CONTINENTAL_VALUE,
            oceanic_value: OCEANIC_VALUE,
            profile_shape: ProfileShape::Smoothstep,
        },
        &ctx,
    );

    eprintln!(
        "FBM calibration probe — single_continent (64², seed=12, 4 plates, 50%% continental)"
    );
    eprintln!("σ_fbm_isolated on largest continental plate (interior t > 0.5):");
    eprintln!("  amp \\ scale   0.05      0.10      0.20      0.25");

    let largest_pid = (0..plates.num_plates)
        .filter(|&pid| matches!(plates.per_plate_type[pid], PlateType::Continental))
        .max_by_key(|&pid| (l_plate[pid] * 100.0) as i32)
        .expect("a continental plate");

    for &amp in &[0.10_f64, 0.15, 0.20, 0.25] {
        eprint!("  amp={:.2}   ", amp);
        for &scale in &[0.05_f64, 0.10, 0.20, 0.25] {
            let s_fbm = init_s_field(
                InitMode::RadialProfileWithFBM {
                    continental_value: CONTINENTAL_VALUE,
                    oceanic_value: OCEANIC_VALUE,
                    profile_shape: ProfileShape::Smoothstep,
                    fbm_amplitude: amp,
                    fbm_octaves: FBM_OCTAVES_DEFAULT,
                    fbm_persistence: FBM_PERSISTENCE_DEFAULT,
                    fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
                    fbm_scale: scale,
                    fbm_seed: FBM_SEED_DEFAULT,
                },
                &ctx,
            );
            let (_, _, _, sigma_fbm_iso) =
                interior_stats(largest_pid, &plates, &bfs.distance, l_plate[largest_pid], &s_radial, &s_fbm);
            eprint!("{:>8.4}  ", sigma_fbm_iso);
        }
        eprintln!();
    }
    eprintln!(
        "Largest continental plate: pid={}, L_plate={:.1} cells",
        largest_pid, l_plate[largest_pid]
    );
}

/// Compute (cells_int, σ_radial, σ_total, σ_fbm_iso) over the
/// interior cells (continental, plate=pid, t > 0.5) of plate `pid`.
fn interior_stats(
    pid: usize,
    plates: &VoronoiPlates,
    dist: &ymir_core::tectonics_v2::field::Field2D,
    l_plate: f64,
    s_radial: &ymir_core::tectonics_v2::field::Field2D,
    s_fbm: &ymir_core::tectonics_v2::field::Field2D,
) -> (usize, f64, f64, f64) {
    let mut interior_radial: Vec<f64> = Vec::new();
    let mut interior_total: Vec<f64> = Vec::new();
    let mut interior_fbm_iso: Vec<f64> = Vec::new();
    for j in 0..NY {
        for i in 0..NX {
            if plates.plate_id.get(i, j) as usize != pid {
                continue;
            }
            if !matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                continue;
            }
            let d = dist.get(i, j);
            if !d.is_finite() {
                continue;
            }
            let t = d / l_plate;
            if t <= 0.5 {
                continue;
            }
            let r = s_radial.get(i, j);
            let f = s_fbm.get(i, j);
            interior_radial.push(r);
            interior_total.push(f);
            interior_fbm_iso.push(f - r);
        }
    }
    let sigma = |xs: &[f64]| -> f64 {
        if xs.is_empty() {
            return 0.0;
        }
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        let var = xs.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / xs.len() as f64;
        var.sqrt()
    };
    (
        interior_radial.len(),
        sigma(&interior_radial),
        sigma(&interior_total),
        sigma(&interior_fbm_iso),
    )
}
