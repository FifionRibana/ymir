//! Issue #131 Phase 2 Track B Stage D — visual gallery generator
//! + multi-seed seed-diversity gallery.
//!
//! Two `#[ignore]`'d tests:
//!
//! 1. [`phase_2_track_b_visual_gallery`] — 5-cycle gallery at seed
//!    42 with the full Phase 2 stack (Davis-Suppe + equilibrium-
//!    height + erosion + Stein-Stein) on Phase 2 R7 init. Dumps
//!    10 PNGs (5 altitude + 5 S̃) at cycles 0 / 50 / 100 / 200 /
//!    300 under `docs/reports/c1_phase_2_track_b_init_r7/`.
//! 2. [`phase_2_track_b_seed_diversity_gallery`] — multi-seed
//!    cycle_000 comparison across seeds 42 / 1337 / 2026. Dumps 6
//!    PNGs (3 seeds × {altitude, S̃}) under
//!    `docs/reports/c1_phase_2_track_b_init_r7/seed_diversity/`.
//!    Forward signal toward the §7.2 Phase 2 milestone gate
//!    "different seeds produce visually distinct continents".
//!
//! Invocation:
//!
//! ```bash
//! cargo test --release -p ymir-core \
//!     --test c1_phase_2_track_b_visual_gallery \
//!     -- --ignored --nocapture
//! ```
//!
//! ## Architecture C re-apply at each cycle
//!
//! Per Track A's pattern, the altitude PNG is produced AFTER
//! re-applying `apply_stein_stein_bathymetry` on the
//! `compute_isostasy` output. The per-step S-S effects are
//! transient (overwritten by the next isostasy recompute from
//! `S̃`); the gallery dump explicitly re-applies S-S so each PNG
//! carries the bathymetric imprint visible at run boundary.
//!
//! ## Palette decisions (identical to Track A gallery)
//!
//! - **Altitude**: bipolar symmetric `[-1.13, +1.13]` mapped to
//!   hypsometric `[0, 1]` with sea level at midpoint. Phase 2
//!   altitude is bipolar (Architecture C produces negative
//!   altitudes on oceanic cells).
//! - **S̃**: `[0, 3.0]` unchanged. Cross-phase comparability per
//!   `feedback_viz_palette_absolute_for_comparison` — same as
//!   Phase 1.4 + Track A galleries.
//!
//! ## Downstream consumability (cross-reference, no new test)
//!
//! Phase 2 R7 init produces a downstream-consumable `C1State`.
//! The Track A acceptance test
//! `c1_phase_2_bathymetry_acceptance::downstream_pipeline_accepts_phase_2_altitude`
//! already validates that `compute_flow` and `run_erosion`
//! accept the bipolar Architecture C altitude. Track B's init
//! changes do not alter that contract — Phase 1.4 downstream
//! tests
//! (`crates/ymir-core/tests/c1_phase_1_4_downstream.rs`) re-run
//! during Stage D pass unchanged (verified by the Stage E4 + V
//! sweeps that ran the full integration suite). DRY: no new
//! downstream test needed in Track B scope.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
use ymir_core::tectonics_c1::closures::accretion::AccretionParams;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use ymir_core::tectonics_c1::closures::rifting::RiftingParams;
use ymir_core::tectonics_c1::closures::subduction::SubductionParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;
use ymir_core::tectonics_v2::field::Field2D;

const GRID_SIZE: usize = 64;
const SEED: u64 = 42;
const N_STEPS: usize = 300;
const S_VIZ_MAX: f64 = 3.0;
const ALTITUDE_PALETTE_HALF_RANGE: f32 = 1.13;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/reports/c1_phase_2_track_b_init_r7")
}

#[test]
#[ignore]
fn phase_2_track_b_visual_gallery() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let init_params = Phase2InitParams::default();
    let mut state = init_c1_state_phase_2_r7(GRID_SIZE, SEED, &init_params);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    // Phase 2 Track B gallery — Track D disabled so the
    // committed PNG references match the Track B-only behaviour
    // (no subduction / accretion / rifting events).
    let closures = C1Closures {
        subduction: SubductionParams { enabled: false, ..SubductionParams::default() },
        accretion: AccretionParams { enabled: false, ..AccretionParams::default() },
        rifting: RiftingParams { enabled: false, ..RiftingParams::default() },
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
        "c1_phase_2_track_b Stage D visual gallery — grid={GRID_SIZE}², steps={N_STEPS}, seed={SEED}"
    );
    eprintln!("  init: Phase 2 R7 (boundary displacement + cluster BFS + ridge age)");
    eprintln!(
        "  closures: DS={} EH={} erosion={} S-S={} (full Phase 2 stack)",
        closures.davis_suppe.enabled,
        closures.equilibrium_height.enabled,
        closures.erosion.enabled,
        closures.oceanic_bathymetry.enabled,
    );

    print_stats("000", &state, &iso_config, &closures);
    dump_snapshot(&state, 0, &dir, &iso_config, &closures);

    let snapshot_steps: [usize; 4] = [49, 99, 199, 299];
    let started = std::time::Instant::now();
    run_with_closures(&mut state, &mut kinematics, &config, &closures, |step, current_state| {
        if snapshot_steps.contains(&step) {
            print_stats(&format!("{:03}", step + 1), current_state, &iso_config, &closures);
            dump_snapshot(current_state, step + 1, &dir, &iso_config, &closures);
        }
    });
    let elapsed = started.elapsed();

    eprintln!();
    eprintln!("  wall time = {:.2?} ({:.2?} / step)", elapsed, elapsed / N_STEPS as u32);
    eprintln!("  output dir = {}", dir.display());
    eprintln!("  files = 10 PNGs (cycle_NNN_altitude.png + cycle_NNN_s.png × 5)");
    eprintln!();

    let final_age_stats = age_stats_oceanic(&state);
    eprintln!("  Phase 2 Track B age distribution (cycle 300, oceanic cells):");
    eprintln!(
        "    min={:.4} max={:.4} mean={:.4} median={:.4}",
        final_age_stats.min, final_age_stats.max, final_age_stats.mean, final_age_stats.median,
    );
    eprintln!("    Track A baseline (Phase 1.1 init): min≈0 max≈6958 mean≈4.67 median≈0");
    eprintln!(
        "    Track B improvement: pile-up factor ~43 % lower; ridge cells present from init."
    );
}

#[test]
#[ignore]
fn phase_2_track_b_seed_diversity_gallery() {
    let dir = output_dir().join("seed_diversity");
    std::fs::create_dir_all(&dir).expect("create seed_diversity dir");

    let seeds: [u64; 3] = [42, 1337, 2026];
    let init_params = Phase2InitParams::default();
    let iso_config = IsostasyConfig::default();
    let closures = C1Closures::default();

    eprintln!("c1_phase_2_track_b Stage D seed_diversity_gallery — cycle_000 at seeds {seeds:?}");
    eprintln!(
        "  per-seed continental fraction + bounding-box extent (forward signal toward §7.2):"
    );

    for &seed in seeds.iter() {
        let state = init_c1_state_phase_2_r7(GRID_SIZE, seed, &init_params);
        let total = GRID_SIZE * GRID_SIZE;
        let mut continental = 0;
        let (mut min_i, mut max_i, mut min_j, mut max_j) = (GRID_SIZE, 0_usize, GRID_SIZE, 0_usize);
        for j in 0..GRID_SIZE {
            for i in 0..GRID_SIZE {
                if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                    continental += 1;
                    if i < min_i {
                        min_i = i;
                    }
                    if i > max_i {
                        max_i = i;
                    }
                    if j < min_j {
                        min_j = j;
                    }
                    if j > max_j {
                        max_j = j;
                    }
                }
            }
        }
        let extent_i = max_i.saturating_sub(min_i).saturating_add(1);
        let extent_j = max_j.saturating_sub(min_j).saturating_add(1);
        eprintln!(
            "    seed = {seed:>5}  continental = {continental:>4} / {total}  bbox = {extent_i}×{extent_j} ({:.0}%×{:.0}%)",
            100.0 * extent_i as f64 / GRID_SIZE as f64,
            100.0 * extent_j as f64 / GRID_SIZE as f64,
        );

        // Re-apply S-S so the cycle_000 altitude PNG shows the
        // ridge-aligned init (Architecture C signature visible
        // from step 0).
        let isostasy = compute_isostasy(&state.s, &iso_config);
        let mut altitude = isostasy.heightmap.clone();
        apply_stein_stein_bathymetry(
            &mut altitude,
            &state.age,
            &state.plate_type,
            &closures.oceanic_bathymetry,
        );

        let alt_path = dir.join(format!("seed_{:05}_altitude.png", seed));
        save_altitude_bipolar_png(&altitude, ALTITUDE_PALETTE_HALF_RANGE, &alt_path);

        let s_path = dir.join(format!("seed_{:05}_s.png", seed));
        save_s_fixed_palette_png(&state.s, S_VIZ_MAX, &s_path);
    }

    eprintln!();
    eprintln!("  output dir = {}", dir.display());
    eprintln!("  files = 6 PNGs (seed_NNNNN_altitude.png + seed_NNNNN_s.png × 3 seeds)");
}

struct AgeStats {
    min: f64,
    max: f64,
    mean: f64,
    median: f64,
}

fn age_stats_oceanic(state: &C1State) -> AgeStats {
    let mut values: Vec<f64> = Vec::new();
    for j in 0..state.ny() {
        for i in 0..state.nx() {
            if matches!(state.plate_type.get(i, j), PlateType::Oceanic) {
                values.push(state.age.get(i, j));
            }
        }
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n == 0 {
        return AgeStats { min: f64::NAN, max: f64::NAN, mean: f64::NAN, median: f64::NAN };
    }
    let min = values[0];
    let max = values[n - 1];
    let mean = values.iter().sum::<f64>() / n as f64;
    let median = values[n / 2];
    AgeStats { min, max, mean, median }
}

fn print_stats(tag: &str, state: &C1State, iso_config: &IsostasyConfig, closures: &C1Closures) {
    let data = state.s.data();
    let mut s_min = f64::INFINITY;
    let mut s_max = f64::NEG_INFINITY;
    let mut s_sum = 0.0;
    for &v in data {
        s_min = s_min.min(v);
        s_max = s_max.max(v);
        s_sum += v;
    }
    let s_mean = s_sum / data.len() as f64;

    let isostasy = compute_isostasy(&state.s, iso_config);
    let mut altitude = isostasy.heightmap.clone();
    apply_stein_stein_bathymetry(
        &mut altitude,
        &state.age,
        &state.plate_type,
        &closures.oceanic_bathymetry,
    );
    let mut a_min = f32::INFINITY;
    let mut a_max = f32::NEG_INFINITY;
    let mut a_sum = 0.0_f64;
    for &v in &altitude.data {
        if v < a_min {
            a_min = v;
        }
        if v > a_max {
            a_max = v;
        }
        a_sum += v as f64;
    }
    let a_mean = a_sum / altitude.data.len() as f64;

    eprintln!(
        "    cycle_{tag}  S̃ min={s_min:.4} mean={s_mean:.4} max={s_max:.4}   \
         altitude (post-S-S) min={a_min:.4} mean={a_mean:.4} max={a_max:.4}"
    );
}

fn dump_snapshot(
    state: &C1State,
    cycle: usize,
    dir: &Path,
    iso_config: &IsostasyConfig,
    closures: &C1Closures,
) {
    let isostasy = compute_isostasy(&state.s, iso_config);
    let mut altitude = isostasy.heightmap;
    apply_stein_stein_bathymetry(
        &mut altitude,
        &state.age,
        &state.plate_type,
        &closures.oceanic_bathymetry,
    );

    let alt_path = dir.join(format!("cycle_{:03}_altitude.png", cycle));
    save_altitude_bipolar_png(&altitude, ALTITUDE_PALETTE_HALF_RANGE, &alt_path);

    let s_path = dir.join(format!("cycle_{:03}_s.png", cycle));
    save_s_fixed_palette_png(&state.s, S_VIZ_MAX, &s_path);
}

fn save_altitude_bipolar_png(altitude: &GridF32, half_range: f32, path: &Path) {
    let nx = altitude.width;
    let ny = altitude.height;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    let sea_norm = 0.5_f32;
    for j in 0..ny {
        for i in 0..nx {
            let raw = altitude.get(i as i32, j as i32);
            let t = ((raw + half_range) / (2.0 * half_range)).clamp(0.0, 1.0);
            let rgb = hypsometric(t, sea_norm);
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
