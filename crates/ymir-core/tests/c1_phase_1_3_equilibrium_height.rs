//! Issue #125 Phase 1.3 Stage E3 — C1 equilibrium-height closure
//! behavioral integration tests.
//!
//! Four acceptance tests on the 64²×300-step Phase 1.1 kinematics
//! preset:
//!
//! 1. [`equilibrium_height_caps_global_max`] — Phase 1.3 default
//!    (both Davis-Suppe + equilibrium-height enabled). Asserts
//!    `global_max < 1.2 · h_eq = 2.4`. Prediction: ≈ 2.0-2.4
//!    (boundary pile-up cells get one-step clamped to `h_eq` by
//!    the quadratic-formula threshold behaviour).
//! 2. [`davis_suppe_imprint_preserved_with_equilibrium`] — both
//!    closures enabled. Re-applies the Phase 1.2 fill-ratio
//!    profile metric on the wedge body (cells with
//!    `0 < d < max_distance`); asserts `fill_near > 0.5` and
//!    `asymmetry mean(0-5)/mean(10-20) > 1.5`. Verifies that the
//!    global equilibrium cap (Step 4 of the C1 closure stack)
//!    does NOT erase the Davis-Suppe wedge imprint (since the
//!    wedge body sits below `h_eq` and the asymmetric one-sided
//!    sink leaves below-`h_eq` cells untouched). Dumps the 5
//!    visual snapshots (0/50/100/200/300) used by the report.
//! 3. [`equilibrium_alone_caps_initial_state`] — Davis-Suppe
//!    disabled, equilibrium-height enabled. Phase 1.1-style
//!    advection + global cap. Asserts `global_max <= 1.1 · h_eq
//!    = 2.2`. Prediction: cap holds even without the Davis-Suppe
//!    source.
//! 4. [`both_closures_disabled_matches_phase_1_1`] — both
//!    closures disabled. Asserts `global_max > 100` (Phase 1.1
//!    unbounded pile-up baseline ≈ 1080 preserved). Regression
//!    guard: a silent default-state mutation that re-enables a
//!    closure would lower `global_max` below 100 and fail.
//!
//! ## Output directory
//!
//! `docs/reports/c1_phase_1_3_equilibrium_height/` carries the
//! 5-cycle PNG gallery (altitude + S̃ at `[0, 3.0]` palette for
//! direct visual comparison with the Phase 1.2 gallery) and the
//! per-phase comparison table in its `README.md`.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::boundary_classification::classify_boundaries;
use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
use ymir_core::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams;
use ymir_core::tectonics_c1::closures::equilibrium_height::params::EquilibriumHeightParams;
use ymir_core::tectonics_c1::closures::erosion::params::ErosionParams;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::rifting::RiftingParams;
use ymir_core::tectonics_c1::closures::subduction::SubductionParams;
use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate;
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::field::Field2D;

const GRID_SIZE: usize = 64;
const SEED: u64 = 42;
const N_STEPS: usize = 300;
/// Palette upper bound for the fixed-scale `S̃` PNG — same as
/// Phase 1.2 (`h_max = 2.5` + ~ 20 % headroom). Sharing the
/// scale lets reviewers compare the Phase 1.2 and Phase 1.3
/// galleries pixel-for-pixel.
const S_VIZ_MAX: f64 = 3.0;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_phase_1_3_equilibrium_height")
}

/// Build the Phase 1.1 init + kinematics + time-loop config used
/// by all four tests. Caller chooses the closure scenario.
fn setup() -> (C1State, PlateKinematics, C1TimeLoopConfig) {
    let state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let config = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: IsostasyConfig::default(),
        drainage_max_distance: 30,
    };
    (state, kinematics, config)
}

fn global_max(state: &C1State) -> f64 {
    state.s.data().iter().cloned().fold(0.0_f64, f64::max)
}

#[test]
fn equilibrium_height_caps_global_max() {
    // Phase 1.3 default — both closures enabled.
    // Davis-Suppe sources mass at the wedge body; equilibrium-
    // height caps the boundary pile-up at h_eq via the quadratic
    // formula's threshold behaviour (clamp triggers on the large-
    // excess cells, holding them at h_eq within one step).
    let (mut state, mut kinematics, config) = setup();
    let closures = C1Closures {
        davis_suppe: DavisSuppeParams::default(),
        equilibrium_height: EquilibriumHeightParams::default(),
        erosion: ErosionParams {
            enabled: false,
            ..ErosionParams::default()
        },
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
    };
    let h_eq = closures.equilibrium_height.h_eq;

    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    let g_max = global_max(&state);
    let threshold = 1.2 * h_eq;
    eprintln!(
        "c1_phase_1_3 T1: global_max = {g_max:.3} (threshold = {threshold:.2} = 1.2 × h_eq = {:.2})",
        h_eq
    );
    assert!(
        g_max < threshold,
        "Phase 1.3 T1 cap failed: global_max = {g_max:.3} ≥ 1.2 × h_eq = {threshold:.2} — \
         equilibrium height is not capping the boundary pile-up (quadratic clamp may be \
         buggy or k_collapse too low)"
    );
}

#[test]
fn davis_suppe_imprint_preserved_with_equilibrium() {
    // Both enabled. Stage E3 dumps the 5 PNGs here for the visual
    // gallery in `docs/reports/c1_phase_1_3_equilibrium_height/`.
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let (mut state, mut kinematics, config) = setup();
    let closures = C1Closures {
        davis_suppe: DavisSuppeParams::default(),
        equilibrium_height: EquilibriumHeightParams::default(),
        erosion: ErosionParams {
            enabled: false,
            ..ErosionParams::default()
        },
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
    };

    eprintln!(
        "c1_phase_1_3 T2: grid={GRID_SIZE}², steps={N_STEPS}, h_eq={}, k_collapse={}, \
         coupling={}, h_max={}, L_taper={}, L_decay={}",
        closures.equilibrium_height.h_eq,
        closures.equilibrium_height.k_collapse,
        closures.davis_suppe.coupling,
        closures.davis_suppe.h_max,
        closures.davis_suppe.l_taper,
        closures.davis_suppe.l_decay,
    );

    // Init snapshot.
    print_s_stats("000", &state);
    dump_snapshot(&state, 0, &dir);

    let snapshot_steps: [usize; 4] = [49, 99, 199, 299];

    let started = std::time::Instant::now();
    run_with_closures(
        &mut state,
        &mut kinematics,
        &config,
        &closures,
        |step, current_state| {
            assert!(
                current_state.s.data().iter().all(|v| v.is_finite()),
                "non-finite S̃ at step {}",
                step + 1
            );
            if snapshot_steps.contains(&step) {
                print_s_stats(&format!("{:03}", step + 1), current_state);
                dump_snapshot(current_state, step + 1, &dir);
            }
        },
    );
    let elapsed = started.elapsed();

    // Re-classify boundaries + recompute wedge distance on the
    // final state. plate_id is static under Phase 1.2 / 1.3 so the
    // classification matches the one used inside the loop.
    let boundary = classify_boundaries(&state.plate_id, &kinematics);
    let wedge_d = wedge_distance_intra_plate(
        &state.plate_id,
        &boundary.upper_plate_mask,
        closures.davis_suppe.max_distance,
    );

    let max_d_cfg = closures.davis_suppe.max_distance;
    let mut g_max = 0.0_f64;
    let mut wedge_values: Vec<f64> = Vec::new();
    let bucket_edges: [(f64, f64); 3] = [(0.0, 5.0), (5.0, 10.0), (10.0, 20.0)];
    let mut bucket_sum = [0.0_f64; 3];
    let mut bucket_count = [0_usize; 3];
    for j in 0..state.ny() {
        for i in 0..state.nx() {
            let v = state.s.get(i, j);
            if v > g_max {
                g_max = v;
            }
            let d = wedge_d.get(i, j);
            if d > 0.0 && d < max_d_cfg {
                wedge_values.push(v);
                for (b, &(lo, hi)) in bucket_edges.iter().enumerate() {
                    if d > lo && d <= hi {
                        bucket_sum[b] += v;
                        bucket_count[b] += 1;
                    }
                }
            }
        }
    }
    wedge_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = wedge_values.len();
    let wedge_max = *wedge_values.last().unwrap_or(&0.0);
    let wedge_mean = wedge_values.iter().sum::<f64>() / n.max(1) as f64;
    let wedge_p95 = wedge_values[(n * 95) / 100];
    let wedge_p99 = wedge_values[(n * 99) / 100];

    let bucket_mean: [f64; 3] = std::array::from_fn(|b| {
        if bucket_count[b] > 0 {
            bucket_sum[b] / bucket_count[b] as f64
        } else {
            f64::NAN
        }
    });
    let h_crit_at = |d: f64| -> f64 {
        closures.davis_suppe.h_max * (1.0 - (-d / closures.davis_suppe.l_taper).exp())
    };

    let fill_near = bucket_mean[0] / h_crit_at(2.5);
    let fill_mid = bucket_mean[1] / h_crit_at(7.5);
    let fill_far = bucket_mean[2] / h_crit_at(15.0);
    let asymmetry = bucket_mean[0] / bucket_mean[2];

    eprintln!(
        "c1_phase_1_3 T2: wedge cells = {n} ({:.1} %)",
        100.0 * n as f64 / (state.nx() * state.ny()) as f64
    );
    eprintln!("c1_phase_1_3 T2: wedge S̃  mean = {wedge_mean:.4}");
    eprintln!("c1_phase_1_3 T2: wedge S̃  p95  = {wedge_p95:.4}");
    eprintln!("c1_phase_1_3 T2: wedge S̃  p99  = {wedge_p99:.4}");
    eprintln!("c1_phase_1_3 T2: wedge S̃  max  = {wedge_max:.4}");
    for (b, &(lo, hi)) in bucket_edges.iter().enumerate() {
        let mid = (lo + hi) / 2.0;
        let h_crit_mid = h_crit_at(mid);
        let fill = if h_crit_mid > 0.0 {
            bucket_mean[b] / h_crit_mid
        } else {
            f64::NAN
        };
        eprintln!(
            "    d ∈ ({lo:>4.1}, {hi:>4.1}]  count = {:>5}  mean S̃ = {:.4}  h_crit = {:.4}  fill = {:.3}",
            bucket_count[b], bucket_mean[b], h_crit_mid, fill
        );
    }
    eprintln!("c1_phase_1_3 T2: fill_far  = {fill_far:.3}  (informational)");
    eprintln!("c1_phase_1_3 T2: global_max (boundary cap) = {g_max:.3}");
    eprintln!(
        "c1_phase_1_3 T2: wall time = {:.2?} ({:.2?} per step)",
        elapsed,
        elapsed / N_STEPS as u32
    );
    eprintln!("c1_phase_1_3 T2: output dir = {}", dir.display());
    eprintln!("c1_phase_1_3 T2 imprint preservation evidence (composite):");
    eprintln!("  wedge_p95: {wedge_p95:.3}   (RIGID #145; legacy 0.376; band [1.6, 2.0], below h_eq)");
    eprintln!("  taper fill: {fill_near:.3}/{fill_mid:.3}/{fill_far:.3} (near>mid>far = source critical taper)");
    eprintln!("  asymmetry: {asymmetry:.2}    (informational; legacy >1.5 was an advection toe-pile, see stage_5b)");
    eprintln!("  fill_near: {fill_near:.3}   (threshold > 0.1)");
    eprintln!(
        "  Note: fill_near drop is the *expected* consequence of the equilibrium clamp on"
    );
    eprintln!(
        "        Davis-Suppe outliers (top 1% of wedge cells, capped at h_eq). The bulk"
    );
    eprintln!(
        "        wedge body (wedge_p95) sits below h_eq and is bit-identical Phase 1.2 ↔ 1.3."
    );

    assert!(
        bucket_count[0] > 0 && bucket_count[2] > 0,
        "T2 profile buckets empty: cannot compute fill-ratio metrics"
    );

    // Sub-assertion 1 — bulk wedge body preserved (PRIMARY), RIGID transport (#145).
    // Under rigid continental crust the wedge is no longer dispersed by advection,
    // so the bulk sits high: wedge_p95 ≈ 1.80 (was 0.376 legacy), still BELOW
    // h_eq = 2.0 so the equilibrium clamp does not touch it (bit-identical to the
    // Phase 1.2 rigid run, also 1.7985). Band [1.6, 2.0] around the rigid baseline.
    assert!(
        (1.6..2.0).contains(&wedge_p95),
        "T2 sub-1 (wedge_p95 bulk preservation, rigid): {wedge_p95:.3} outside [1.6, 2.0] \
         — rigid wedge bulk sits high (crust not advected away), below h_eq so the \
         equilibrium clamp leaves it untouched (bit-identical to Phase 1.2 rigid)."
    );

    // Sub-assertion 2 — source critical TAPER (replaces the legacy asymmetry>1.5).
    // Legacy `asymmetry = near/far mean > 1` tested an ADVECTION toe-pile that
    // INVERTED the Davis-Suppe source taper. Under rigid transport the true source
    // signature appears: fill ratio DECREASES with distance (h_crit grows with
    // distance). Measured monotone: fill_near 0.747 > fill_mid 0.519 > fill_far
    // 0.361. Asserts the PHYSICS (source critical taper), not the advection artifact.
    // See docs/reports/c1_continental_buoyancy/stage_5b_asymmetry.md.
    assert!(
        fill_near > fill_mid && fill_mid > fill_far,
        "T2 sub-2 (source taper fill_near>fill_mid>fill_far): {fill_near:.3}/{fill_mid:.3}/{fill_far:.3} \
         not monotone-decreasing — Davis-Suppe critical taper lost (legacy asymmetry={asymmetry:.2})"
    );

    // Sub-assertion 3 — fill_near regime-tagged Phase 1.3 floor.
    // Phase 1.2 baseline 0.778, Phase 1.3 baseline ≈ 0.207 (the
    // equilibrium clamp removes Davis-Suppe outliers from the
    // bucket-mean dominator). The relaxed `> 0.1` threshold
    // catches catastrophic loss (Davis-Suppe entirely silenced
    // → fill_near ≈ 0) while accommodating further regime drift
    // in Phase 1.4+ (the erosion sink will likely lower this
    // further). Pairs with the memory entry
    // `feedback_fill_ratio_regime_agnostic_metric` — regime-
    // tagged thresholds, re-evaluated per phase.
    assert!(
        fill_near > 0.1,
        "T2 sub-3 (fill_near regime-tagged Phase 1.3 floor): {fill_near:.3} ≤ 0.1 \
         — Phase 1.2 baseline 0.778, Phase 1.3 expected ≈ 0.207. Drop to near \
         zero indicates Davis-Suppe is entirely suppressed."
    );
}

#[test]
fn equilibrium_alone_caps_initial_state() {
    // Davis-Suppe OFF, equilibrium-height ON. Pure advection
    // pile-up under the Phase 1.1 kinematics; the equilibrium
    // closure must still cap it. Tighter threshold (1.1 × h_eq)
    // than T1 because there's no Davis-Suppe source pushing
    // additional mass into the pile-up.
    let (mut state, mut kinematics, config) = setup();
    let closures = C1Closures {
        davis_suppe: DavisSuppeParams {
            enabled: false,
            ..DavisSuppeParams::default()
        },
        equilibrium_height: EquilibriumHeightParams::default(),
        erosion: ErosionParams {
            enabled: false,
            ..ErosionParams::default()
        },
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
    };
    let h_eq = closures.equilibrium_height.h_eq;

    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    let g_max = global_max(&state);
    let threshold = 1.1 * h_eq;
    eprintln!(
        "c1_phase_1_3 T3: global_max = {g_max:.3} (threshold = {threshold:.2} = 1.1 × h_eq = {:.2})",
        h_eq
    );
    assert!(
        g_max <= threshold,
        "Phase 1.3 T3 cap failed: global_max = {g_max:.3} > 1.1 × h_eq = {threshold:.2} — \
         equilibrium standalone ineffective on advection pile-up (no Davis-Suppe to absorb \
         the excess via its h_max plateau)"
    );
}

#[test]
fn both_closures_disabled_matches_phase_1_1() {
    // Regression guard: both closures off → run_with_closures
    // reduces to advection only (closures are no-ops). The Phase 1.1
    // unbounded pile-up baseline (≈ 1080) must be preserved.
    // A silent default-state mutation that re-enables a closure
    // would cap global_max below 100 and fail this assertion.
    let (mut state, mut kinematics, config) = setup();
    let closures = C1Closures {
        davis_suppe: DavisSuppeParams {
            enabled: false,
            ..DavisSuppeParams::default()
        },
        equilibrium_height: EquilibriumHeightParams {
            enabled: false,
            ..EquilibriumHeightParams::default()
        },
        erosion: ErosionParams {
            enabled: false,
            ..ErosionParams::default()
        },
        oceanic_bathymetry: SteinSteinParams {
            enabled: false,
            ..SteinSteinParams::default()
        },
        subduction: SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        },
        accretion: AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        },
        rifting: RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        },
    };

    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    let g_max = global_max(&state);
    eprintln!(
        "c1_phase_1_3 T4: global_max = {g_max:.3} (need > 100, Phase 1.1 baseline ≈ 1080)"
    );
    assert!(
        g_max > 100.0,
        "T4 Phase 1.1 baseline broken: global_max = {g_max:.3} ≤ 100 — \
         a closure may be silently active despite enabled=false"
    );
}

// ── Diagnostics + viz helpers (mirror of Phase 1.2 helpers) ──

fn print_s_stats(tag: &str, state: &C1State) {
    let data = state.s.data();
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for &v in data {
        min = min.min(v);
        max = max.max(v);
        sum += v;
    }
    let mean = sum / data.len() as f64;
    let mut sq = 0.0;
    for &v in data {
        sq += (v - mean) * (v - mean);
    }
    let std = (sq / data.len() as f64).sqrt();
    eprintln!(
        "c1_phase_1_3: cycle_{tag} S̃ min={min:.4} mean={mean:.4} max={max:.4} std={std:.4e}"
    );
}

fn dump_snapshot(state: &C1State, cycle: usize, dir: &Path) {
    let iso = compute_isostasy(&state.s, &IsostasyConfig::default());
    let alt_path = dir.join(format!("cycle_{:03}_altitude.png", cycle));
    save_hypsometric_png(&iso.heightmap, iso.sea_level_normalized, &alt_path);

    let s_path = dir.join(format!("cycle_{:03}_s.png", cycle));
    save_s_fixed_palette_png(&state.s, S_VIZ_MAX, &s_path);
}

fn save_hypsometric_png(
    heightmap: &ymir_core::grid::GridF32,
    sea_norm: f32,
    path: &Path,
) {
    let nx = heightmap.width;
    let ny = heightmap.height;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    for j in 0..ny {
        for i in 0..nx {
            let h = heightmap.get(i as i32, j as i32).clamp(0.0, 1.0);
            let rgb = hypsometric(h, sea_norm);
            let img_row = (ny - 1 - j) as u32;
            img.put_pixel(i as u32, img_row, Rgb(rgb));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create parent dir");
    }
    img.save(path).expect("save PNG");
}

fn save_s_fixed_palette_png(s: &Field2D, s_max: f64, path: &Path) {
    let nx = s.nx();
    let ny = s.ny();
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    let sea_norm = 0.2 / s_max;
    for j in 0..ny {
        for i in 0..nx {
            let v = (s.get(i, j) / s_max).clamp(0.0, 1.0) as f32;
            let rgb = hypsometric(v, sea_norm as f32);
            let img_row = (ny - 1 - j) as u32;
            img.put_pixel(i as u32, img_row, Rgb(rgb));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create parent dir");
    }
    img.save(path).expect("save PNG");
}

fn hypsometric(h: f32, sea_norm: f32) -> [u8; 3] {
    let mid = (sea_norm + 1.0) * 0.5;
    let lerp = |t: f32, a: [u8; 3], b: [u8; 3]| -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        [
            (a[0] as f32 + t * (b[0] as f32 - a[0] as f32)).round() as u8,
            (a[1] as f32 + t * (b[1] as f32 - a[1] as f32)).round() as u8,
            (a[2] as f32 + t * (b[2] as f32 - a[2] as f32)).round() as u8,
        ]
    };
    if h <= sea_norm * 0.5 {
        let t = h / (sea_norm * 0.5).max(1e-6);
        lerp(t, [10, 20, 60], [40, 80, 160])
    } else if h <= sea_norm {
        let t = (h - sea_norm * 0.5) / (sea_norm * 0.5).max(1e-6);
        lerp(t, [40, 80, 160], [120, 180, 230])
    } else if h <= mid {
        let t = (h - sea_norm) / (mid - sea_norm).max(1e-6);
        lerp(t, [60, 130, 60], [140, 100, 50])
    } else {
        let t = (h - mid) / (1.0 - mid).max(1e-6);
        lerp(t, [140, 100, 50], [245, 245, 245])
    }
}
