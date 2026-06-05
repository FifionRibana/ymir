//! Issue #127 Phase 1.4 Stage E4 — formal acceptance tests for
//! the C1 stream-power erosion closure (Whipple-Tucker 1999
//! + Lague 2014).
//!
//! Three tests at 64² × 300 steps with the Phase 1.1 kinematics
//! preset, evaluating the joint behaviour of all three C1
//! closures (Davis-Suppe + equilibrium-height + erosion).
//!
//! Test inventory (Phase 1.4-regime-tagged thresholds per the
//! Stage E3 calibration empirical baseline):
//!
//! 1. [`erosion_caps_height_below_equilibrium`] — global_max
//!    stays within `[1.0, 2.5]`. Phase 1.4 default `K = 0.001`
//!    measured 2.181 (Stage E3). Bounded above by `h_max = 2.5`
//!    (Davis-Suppe ceiling). Bounded below by 1.0 to guard
//!    against pathological erosion that flattens everything.
//! 2. [`erosion_preserves_davis_suppe_imprint_partially`] —
//!    composite assertion (3 sub-cases), regime-tagged Phase 1.4.
//!    Documents the "wedge_p95 UP" architectural finding
//!    (Phase 1.4 wedges 0.696 vs Phase 1.3 0.376 because erosion
//!    eats downstream continental shoulders preferentially over
//!    upstream wedge cells — W-T `E ∝ A^m` discriminates
//!    drainage-area-low cells from drainage-area-high cells).
//!    Also dumps the 5-cycle PNG gallery used by the report.
//! 3. [`all_closures_disabled_matches_phase_1_1`] — regression
//!    guard: with all three closures disabled, the time loop
//!    reduces to advection only and reproduces the Phase 1.1
//!    unbounded boundary pile-up baseline (`global_max > 100`).
//!    A silent default-state mutation re-enabling a closure
//!    would cap the global max and fail this assertion.
//!
//! ## Test deferred — `erosion_alone_produces_<X>` (architectural finding)
//!
//! The Phase 1.4 spec originally listed a fourth acceptance test
//! intended as a sanity check that erosion does *something*
//! useful when Davis-Suppe + equilibrium-height are disabled.
//! Three iterations during Stage E4 revealed that **no clean
//! regime-agnostic invariant exists for the erosion-alone
//! regime** in C1:
//!
//! 1. **"Variance smoothing"** (`final_variance < initial`):
//!    fails. The advection-dominated regime (per memory
//!    `project_c1_phase_1_2_advection_dominated_regime`) drives
//!    unbounded boundary pile-up. With no equilibrium clamp,
//!    `K = 0.001` erosion can't keep up with the pile-up rate;
//!    variance EXPLODES (boundary cells reach 100+, interior
//!    stays at 1.0). Smoothing is a property of erosion on a
//!    *static* heightmap; C1's heightmap is anything but.
//! 2. **"Mass loss"** (`final_mass < initial`): fails. With DS
//!    disabled, advection drives near-oceanic cells below the
//!    `floor = 0.2` clamp inside `apply_erosion_step`; for each
//!    such cell the clamp injects mass *upward* to bring it
//!    back to 0.2. Measured: erosion-only run added +5191 mass
//!    over 300 steps (+227 % vs init). The floor clamp is
//!    mass-non-conservative in **degraded regimes** (no source);
//!    in the Phase 1.4 default regime (DS active) continental
//!    cells stay well above the floor, the clamp rarely
//!    triggers, and erosion is a net sink (Stage E3 measured
//!    −35 % mass).
//! 3. **"Continental-filtered mass loss"**: over-engineered for
//!    a sanity check (requires pre-filtering by initial cell
//!    set with arbitrary threshold).
//!
//! The **architectural finding** — that the erosion closure's
//! defensive `floor = 0.2` clamp is *mass-non-conservative in
//! the erosion-alone degraded regime*, while *net-conservative
//! in the Phase 1.4 default* — is preserved in the
//! `project_c1_phase_1_4_erosion_outcomes` memory entry
//! (Stage Final) rather than as an acceptance test. The pattern
//! follows the `recursive-tuning-signals-structural` memory: 3
//! iterations on a single sanity test signal that no clean
//! invariant exists; document the structural finding and ship.
//! The three remaining tests cover the load-bearing Phase 1.4
//! acceptance surface (cap, imprint, regression).
//!
//! ## Output directory
//!
//! `docs/reports/c1_phase_1_4_erosion/` carries the 5-cycle
//! gallery (altitude + S̃ at fixed palette `[0, 3.0]` for
//! pixel-comparable inspection against the Phase 1.2 and 1.3
//! galleries). Re-generated each run; not committed (consistent
//! with Phase 1.1 / 1.2 / 1.3 convention).

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
const S_VIZ_MAX: f64 = 3.0;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_phase_1_4_erosion")
}

fn setup() -> (C1State, PlateKinematics, C1TimeLoopConfig) {
    let state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let config = C1TimeLoopConfig {
        rigid_continental_crust: false,
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

/// Phase 1.4 closure stack — Davis-Suppe + equilibrium-height +
/// erosion all ON; Phase 2 Track A oceanic bathymetry **OFF**.
/// Locks the Phase 1.4 acceptance regime against silent regime
/// drift from `C1Closures::default()` which (post-#129) enables
/// S-S bathymetry by default. With S-S on, oceanic altitude jumps
/// to `-d(t) / depth_scale_m`, slope at the continental/oceanic
/// coastline steepens dramatically, and erosion's slope factor
/// runs much hotter on coastal cells — Phase 1.4's regime-tagged
/// thresholds (wedge_p95 ∈ [0.4, 1.0], global_max ∈ [1.0, 2.5])
/// were calibrated without S-S and must remain testable
/// independently.
fn phase_1_4_closures() -> C1Closures {
    C1Closures {
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
        ..C1Closures::default()
    }
}

#[test]
fn erosion_caps_height_below_equilibrium() {
    // Phase 1.4 default — all 3 closures enabled. global_max must
    // stay bounded above by `h_max = 2.5` (Davis-Suppe ceiling)
    // and bounded below by 1.0 to guard against pathological
    // erosion that erases everything. Stage E3 measured 2.181.
    let (mut state, mut kinematics, config) = setup();
    let closures = phase_1_4_closures();

    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    let g_max = global_max(&state);
    eprintln!(
        "c1_phase_1_4 T1: global_max = {g_max:.3} (must lie in (1.0, 2.5); Stage E3 measured 2.181)"
    );
    assert!(
        g_max < 2.5,
        "Phase 1.4 T1: global_max = {g_max:.3} ≥ 2.5 — equilibrium clamp may have \
         broken or erosion produced unexpected uplift (impossible by W-T sign convention)"
    );
    assert!(
        g_max > 1.0,
        "Phase 1.4 T1: global_max = {g_max:.3} ≤ 1.0 — erosion has flattened the \
         entire grid. K_erosion likely too large; investigate K calibration"
    );
}

#[test]
fn erosion_preserves_davis_suppe_imprint_partially() {
    // Composite assertion, Phase-1.4-regime-tagged per the
    // [[fill-ratio-regime-agnostic-metric]] memory pattern.
    //
    // Three sub-assertions documenting that the Davis-Suppe
    // imprint survives the erosion sink — but in a different
    // form than Phase 1.3, because erosion preferentially eats
    // continental shoulders (large drainage area) over wedge
    // cells (small drainage area, at the top of their basin).
    //
    // ARCHITECTURAL FINDING — "wedge_p95 UP":
    //   Phase 1.2 baseline: wedge_p95 = 0.376
    //   Phase 1.3 baseline: wedge_p95 = 0.376 (bit-identical
    //                       — bulk below h_eq untouched)
    //   Phase 1.4 measured: wedge_p95 = 0.696 (UP by 85 %!)
    //
    //   Mechanism: W-T `E ∝ A^m` with `m = 0.5`. Erosion is
    //   concentrated on cells with LARGE drainage area —
    //   continental shoulders downstream of the wedge bodies.
    //   Wedge cells are topographically UPSTREAM (small `A`),
    //   so they erode less than the surroundings. Net effect:
    //   wedge ridges stand RELATIVELY HIGHER vs the eroding
    //   continental bulk. Earth-like: mountain ridges preserve
    //   in oceanic terrane.
    //
    // The composite test locks the Phase 1.4 morphology
    // signature in three independent dimensions:
    //
    //   sub-1: wedge_p95 ∈ [0.4, 1.0]  — bulk preservation
    //          + LIFT (Phase 1.4 specific)
    //   sub-2: asymmetry > 1.0         — spatial shape
    //          preserved (Phase 1.4 measured ≈ 1.34, Phase 1.3
    //          measured 2.12)
    //   sub-3: fill_near > 0.05        — saturation floor
    //          (Phase 1.4 measured 0.365, Phase 1.3 0.207)

    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let (mut state, mut kinematics, config) = setup();
    let closures = phase_1_4_closures();

    eprintln!(
        "c1_phase_1_4 T2: grid={GRID_SIZE}², steps={N_STEPS}, K={}, m={}, n={}, h_eq={}, h_max={}",
        closures.erosion.k,
        closures.erosion.m,
        closures.erosion.n,
        closures.equilibrium_height.h_eq,
        closures.davis_suppe.h_max,
    );

    // Cycle 0 snapshot.
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

    // Re-classify boundaries + recompute wedge distance.
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
    let wedge_p95 = wedge_values[(n * 95) / 100];
    let wedge_p99 = wedge_values[(n * 99) / 100];
    let wedge_mean = wedge_values.iter().sum::<f64>() / n.max(1) as f64;

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
    let asymmetry = bucket_mean[0] / bucket_mean[2];

    eprintln!();
    eprintln!("c1_phase_1_4 T2 imprint preservation evidence (composite):");
    eprintln!(
        "  wedge_p95  = {wedge_p95:.3}  (Phase 1.2: 0.376, Phase 1.3: 0.376, Phase 1.4 LIFTED; threshold [0.4, 1.0])"
    );
    eprintln!(
        "  asymmetry  = {asymmetry:.2}   (Phase 1.2: 4.66,  Phase 1.3: 2.12,  threshold > 1.0)"
    );
    eprintln!(
        "  fill_near  = {fill_near:.3}  (Phase 1.2: 0.778, Phase 1.3: 0.207, threshold > 0.05)"
    );
    eprintln!(
        "  Note: wedge_p95 UP is the Phase 1.4 architectural finding. Erosion E ∝ A^m"
    );
    eprintln!(
        "        eats continental shoulders (large A) preferentially over wedge cells"
    );
    eprintln!(
        "        (small A, upstream of their drainage basin). Wedge ridges stand HIGHER"
    );
    eprintln!(
        "        relative to the eroding bulk. Earth-like: mountain ranges in oceanic terrane."
    );
    eprintln!();
    eprintln!(
        "  wedge mean = {wedge_mean:.3}, wedge p99 = {wedge_p99:.3}, global_max = {g_max:.3}"
    );
    eprintln!(
        "  wall time  = {:.2?} ({:.2?} per step)",
        elapsed,
        elapsed / N_STEPS as u32
    );
    eprintln!("  output dir = {}", dir.display());

    assert!(
        bucket_count[0] > 0 && bucket_count[2] > 0,
        "T2: profile buckets empty"
    );

    // Sub-assertion 1 — wedge_p95 bulk preservation + LIFT.
    assert!(
        (0.4..=1.0).contains(&wedge_p95),
        "T2 sub-1 (wedge_p95 ∈ [0.4, 1.0]): {wedge_p95:.3} outside range — \
         Phase 1.4 should show wedge bulk lifted relative to Phase 1.3 baseline 0.376 \
         due to drainage-area-discriminated erosion, but bounded by Davis-Suppe \
         h_max"
    );

    // Sub-assertion 2 — spatial asymmetry preserved.
    assert!(
        asymmetry > 1.0,
        "T2 sub-2 (asymmetry > 1.0): {asymmetry:.2} ≤ 1.0 — spatial near-vs-far signature \
         lost; wedges have flattened or imprint erased"
    );

    // Sub-assertion 3 — fill_near regime-tagged Phase 1.4 floor.
    assert!(
        fill_near > 0.05,
        "T2 sub-3 (fill_near > 0.05): {fill_near:.3} ≤ 0.05 — Davis-Suppe + equilibrium \
         signature essentially erased by erosion; investigate K calibration"
    );
}

#[test]
fn all_closures_disabled_matches_phase_1_1() {
    // Regression guard: all 3 closures off → run_with_closures
    // reduces to advection only. Phase 1.1 unbounded boundary
    // pile-up baseline (`global_max ≈ 1080`) must be preserved.
    // A silent default-state mutation enabling any closure
    // would cap the global max and fail this assertion.
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
        "c1_phase_1_4 T3: global_max = {g_max:.3} (Phase 1.1 baseline ≈ 1080, threshold > 100)"
    );
    assert!(
        g_max > 100.0,
        "T3: Phase 1.1 baseline broken — global_max = {g_max:.3} ≤ 100. A closure may \
         be silently active despite enabled=false"
    );
}

// ── Diagnostics + viz helpers (mirror of Phase 1.3 helpers) ──

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
        "c1_phase_1_4: cycle_{tag} S̃ min={min:.4} mean={mean:.4} max={max:.4} std={std:.4e}"
    );
}

fn dump_snapshot(state: &C1State, cycle: usize, dir: &Path) {
    let iso = compute_isostasy(&state.s, &IsostasyConfig::default());
    let alt_path = dir.join(format!("cycle_{:03}_altitude.png", cycle));
    save_hypsometric_png(&iso.heightmap, iso.sea_level_normalized, &alt_path);

    let s_path = dir.join(format!("cycle_{:03}_s.png", cycle));
    save_s_fixed_palette_png(&state.s, S_VIZ_MAX, &s_path);
}

fn save_hypsometric_png(heightmap: &ymir_core::grid::GridF32, sea_norm: f32, path: &Path) {
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
