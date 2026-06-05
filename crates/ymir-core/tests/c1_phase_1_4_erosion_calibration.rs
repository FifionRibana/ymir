//! Issue #127 Phase 1.4 Stage E3 — erosion K calibration helper.
//!
//! Single `#[ignore]`'d test that exercises the full Phase 1.4
//! pipeline (all 3 closures enabled, default K = 0.001) at
//! 64² × 300 steps and produces:
//!
//! 1. A 5-cycle PNG gallery (cycle 0 / 50 / 100 / 200 / 300, 2
//!    palettes per cycle = 10 files) under
//!    `docs/reports/c1_phase_1_4_erosion/`.
//! 2. Per-cycle stats on `S̃` distribution.
//! 3. Final-cycle Phase-1.3-style metrics: `global_max`,
//!    `wedge_p95`, `wedge_p99`, fill-ratio bucket profile,
//!    `asymmetry mean(near) / mean(far)`.
//! 4. Erosion-specific metrics: total mass removed by erosion
//!    over the run (via mass-balance vs Phase 1.3 cycle 300
//'    baseline).
//!
//! Mark as `#[ignore]` because it generates artefacts and prints
//! verbose diagnostics — not part of the routine regression sweep.
//! Invocation:
//!
//! ```bash
//! cargo test --release -p ymir-core \
//!     --test c1_phase_1_4_erosion_calibration \
//!     -- --ignored --nocapture
//! ```
//!
//! The `c1_phase_1_4_erosion.rs` file (Stage E4 scope) will carry
//! the formal acceptance tests with assertions; this file is a
//! calibration tool retained for future K reviews.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::boundary_classification::classify_boundaries;
use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
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

#[test]
#[ignore]
fn erosion_calibration_visual_review() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let mut state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    // Phase 1.4 calibration tool — keep S-S OFF so the K calibration
    // record matches the Phase 1.4 PNG gallery committed in
    // `docs/reports/c1_phase_1_4_erosion/`. Phase 2 Track A gets its
    // own gallery (Stage D, Issue #129).
    let closures = C1Closures {
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
    };
    let config = C1TimeLoopConfig {
        rigid_continental_crust: false,
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };

    eprintln!(
        "c1_phase_1_4 E3 calibration: grid={GRID_SIZE}², steps={N_STEPS}, K={}, m={}, n={}, \
         floor={}, h_eq={}, h_max={}",
        closures.erosion.k,
        closures.erosion.m,
        closures.erosion.n,
        closures.erosion.floor,
        closures.equilibrium_height.h_eq,
        closures.davis_suppe.h_max,
    );

    let initial_mass: f64 = state.s.data().iter().sum();

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
            if snapshot_steps.contains(&step) {
                print_s_stats(&format!("{:03}", step + 1), current_state);
                dump_snapshot(current_state, step + 1, &dir);
            }
        },
    );
    let elapsed = started.elapsed();
    let final_mass: f64 = state.s.data().iter().sum();
    let mass_delta = final_mass - initial_mass;
    let mass_delta_rel = mass_delta / initial_mass;

    // Final-cycle Phase 1.3-style metrics.
    let boundary = classify_boundaries(&state.plate_id, &kinematics);
    let wedge_d = wedge_distance_intra_plate(
        &state.plate_id,
        &boundary.upper_plate_mask,
        closures.davis_suppe.max_distance,
    );
    let max_d_cfg = closures.davis_suppe.max_distance;

    let mut global_max = 0.0_f64;
    let mut wedge_values: Vec<f64> = Vec::new();
    let bucket_edges: [(f64, f64); 3] = [(0.0, 5.0), (5.0, 10.0), (10.0, 20.0)];
    let mut bucket_sum = [0.0_f64; 3];
    let mut bucket_count = [0_usize; 3];
    for j in 0..state.ny() {
        for i in 0..state.nx() {
            let v = state.s.get(i, j);
            if v > global_max {
                global_max = v;
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
    let fill_far = bucket_mean[2] / h_crit_at(15.0);
    let asymmetry = bucket_mean[0] / bucket_mean[2];

    eprintln!();
    eprintln!("c1_phase_1_4 E3 calibration — Phase 1.3 vs Phase 1.4 metric comparison:");
    eprintln!(
        "  initial mass   = {initial_mass:.2}     (sum S̃ over {} cells)",
        state.nx() * state.ny()
    );
    eprintln!(
        "  final mass     = {final_mass:.2}     (delta {mass_delta:+.2} = {:+.3} % vs init)",
        100.0 * mass_delta_rel
    );
    eprintln!(
        "  global_max     = {global_max:.3}    (Phase 1.3 baseline: 2.18; should be ≤ h_eq = 2.0 ± clamp tolerance)"
    );
    eprintln!(
        "  wedge cells    = {n} ({:.1} %)",
        100.0 * n as f64 / (state.nx() * state.ny()) as f64
    );
    eprintln!(
        "  wedge mean     = {wedge_mean:.4}   (Phase 1.3 baseline: 0.166)"
    );
    eprintln!(
        "  wedge p95      = {wedge_p95:.4}   (Phase 1.3 baseline: 0.376)"
    );
    eprintln!(
        "  wedge p99      = {wedge_p99:.4}   (Phase 1.3 baseline: 2.17)"
    );
    eprintln!(
        "  wedge max      = {wedge_max:.4}   (Phase 1.3 baseline: 2.18)"
    );
    eprintln!();
    eprintln!("  per-distance bucket fill profile (mean S̃ per bucket vs h_crit at mid):");
    for (b, &(lo, hi)) in bucket_edges.iter().enumerate() {
        let mid = (lo + hi) / 2.0;
        let h_crit_mid = h_crit_at(mid);
        let fill = if h_crit_mid > 0.0 {
            bucket_mean[b] / h_crit_mid
        } else {
            f64::NAN
        };
        eprintln!(
            "    d ∈ ({lo:>4.1}, {hi:>4.1}]  count={:>4}  mean S̃={:.4}  h_crit(mid)={:.4}  fill={:.3}",
            bucket_count[b], bucket_mean[b], h_crit_mid, fill
        );
    }
    eprintln!();
    eprintln!(
        "  fill_near (d∈0-5)   = {fill_near:.3}   (Phase 1.3 baseline: 0.207)"
    );
    eprintln!(
        "  fill_far  (d∈10-20) = {fill_far:.3}   (Phase 1.3 baseline: 0.046)"
    );
    eprintln!(
        "  asymmetry (near/far) = {asymmetry:.3}   (Phase 1.3 baseline: 2.12)"
    );
    eprintln!();
    eprintln!(
        "  wall time      = {:.2?} ({:.2?} per step) — Phase 1.3 baseline: 29 ms / 96 µs",
        elapsed,
        elapsed / N_STEPS as u32
    );
    eprintln!("  output dir     = {}", dir.display());
}

// ── Diagnostics + viz helpers (mirror Phase 1.3 helpers) ──

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
        "c1_phase_1_4 E3: cycle_{tag}  S̃ min={min:.4} mean={mean:.4} max={max:.4} std={std:.4e}"
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
