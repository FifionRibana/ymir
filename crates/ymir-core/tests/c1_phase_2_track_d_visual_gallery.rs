//! Issue #132 Phase 2 Track D Stage A — visual gallery generator.
//!
//! Two `#[ignore]`'d tests that exercise the full Phase 2 + Track D
//! stack at 64²×300 steps under Phase 2 R7 init and produce:
//!
//! 1. A 5-cycle PNG gallery at seed 42 (cycle 0 / 50 / 100 / 200 /
//!    300, 2 fields per cycle = 10 files) under
//!    `docs/reports/c1_phase_2_track_d_boundary_evolution/`.
//! 2. A multi-seed diversity gallery at cycle 300 across 3 seeds
//!    (6 files: altitude + S̃ × 3) under the `seed_diversity/`
//!    subdir. Forward signal toward §7.2 cross-track Phase 2
//!    milestone gate.
//!
//! Architecture C re-apply S-S at each cycle preserved from Track A
//! pattern. Palette continuity preserved from Track A/B (Q-V.3
//! Option A): altitude `[-1.13, +1.13]` symmetric, `S̃` `[0, 3.0]`.
//! Architectural finding documented inline if clip artifacts
//! visible at the palette edges (pattern Phase 1.4 floor-clamp
//! regime-specific). Invocation:
//!
//! ```bash
//! cargo test --release -p ymir-core \
//!     --test c1_phase_2_track_d_visual_gallery \
//!     -- --ignored --nocapture
//! ```
//!
//! ## Architecture C re-application reminder
//!
//! Stein-Stein adjustments to altitude are TRANSIENT within
//! `run_with_closures` — the next `compute_isostasy` call
//! regenerates altitude from S̃, overwriting the in-loop S-S
//! adjustment. To inspect the bathymetric imprint at each
//! snapshot, S-S is re-applied at the snapshot point (matching
//! the Track A Stage D pattern).
//!
//! ## PNG file naming
//!
//! - Main gallery: `cycle_NNN_altitude.png`, `cycle_NNN_s.png` for
//!   `NNN ∈ {000, 050, 100, 200, 300}`.
//! - Diversity gallery: `seed_NNNNN_altitude.png`,
//!   `seed_NNNNN_s.png` for `NNNNN ∈ {00042, 01337, 02026}` in
//!   the `seed_diversity/` subdir.
//!
//! Files NOT committed (Phase 1.x + Track A/B convention). PNG
//! references are regenerated each run.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use ymir_core::tectonics_c1::init_r7::{init_c1_state_phase_2_r7, Phase2InitParams};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::field::Field2D;

const GRID_SIZE: usize = 64;
const SEED: u64 = 42;
const N_STEPS: usize = 300;
const S_VIZ_MAX: f64 = 3.0;
const ALTITUDE_PALETTE_HALF_RANGE: f32 = 1.13;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_phase_2_track_d_boundary_evolution")
}

#[test]
#[ignore]
fn phase_2_track_d_visual_gallery() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let init_params = Phase2InitParams::default();
    let mut state = init_c1_state_phase_2_r7(GRID_SIZE, SEED, &init_params);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };

    eprintln!(
        "c1_phase_2_track_d Stage A visual gallery — grid={GRID_SIZE}², steps={N_STEPS}, seed={SEED}"
    );
    eprintln!(
        "  closures: DS={} EH={} erosion={} S-S={} subduction={} accretion={} rifting={}",
        closures.davis_suppe.enabled,
        closures.equilibrium_height.enabled,
        closures.erosion.enabled,
        closures.oceanic_bathymetry.enabled,
        closures.subduction.enabled,
        closures.accretion.enabled,
        closures.rifting.enabled,
    );
    eprintln!(
        "  palettes: altitude [-{ALTITUDE_PALETTE_HALF_RANGE}, +{ALTITUDE_PALETTE_HALF_RANGE}], S̃ [0, {S_VIZ_MAX}]"
    );

    print_stats("000", &state, &iso_config, &closures);
    dump_snapshot(&state, 0, &dir, &iso_config, &closures);

    let snapshot_steps: [usize; 4] = [49, 99, 199, 299];
    let started = std::time::Instant::now();
    run_with_closures(
        &mut state,
        &mut kinematics,
        &config,
        &closures,
        |step, current_state| {
            if snapshot_steps.contains(&step) {
                print_stats(
                    &format!("{:03}", step + 1),
                    current_state,
                    &iso_config,
                    &closures,
                );
                dump_snapshot(current_state, step + 1, &dir, &iso_config, &closures);
            }
        },
    );
    let elapsed = started.elapsed();
    let per_step_us = elapsed.as_secs_f64() * 1.0e6 / N_STEPS as f64;

    eprintln!();
    eprintln!(
        "  wall time      = {:.2?} ({:.2?} / step ≈ {per_step_us:.1} µs)",
        elapsed,
        elapsed / N_STEPS as u32
    );
    eprintln!("  output dir     = {}", dir.display());
    eprintln!(
        "  files          = 10 PNGs (cycle_NNN_altitude.png + cycle_NNN_s.png × 5)"
    );
    eprintln!();
    eprintln!("  Phase 3 optimisation forward signal:");
    if per_step_us > 800.0 {
        eprintln!(
            "    ARCHITECTURAL FINDING: per-step cost {per_step_us:.1} µs exceeds Stage E4 budget 800 µs."
        );
        eprintln!(
            "    Source likely Track D per-step boundary recompute (~200 µs)."
        );
        eprintln!(
            "    Phase 3+ optimisation: conditional skip when no Track D event fired previous step."
        );
    } else {
        eprintln!(
            "    Per-step cost {per_step_us:.1} µs within 800 µs Stage E4 budget. Phase 3 optimisation NOT prioritised."
        );
    }
    eprintln!();
    eprintln!("  Post-run plate count = {} (was {} at init)", kinematics.velocities.len(), state.num_plates);
}

#[test]
#[ignore]
fn phase_2_track_d_seed_diversity_gallery() {
    let dir = output_dir().join("seed_diversity");
    std::fs::create_dir_all(&dir).expect("create seed_diversity dir");

    let seeds: [u64; 3] = [42, 1337, 2026];
    let init_params = Phase2InitParams::default();
    let iso_config = IsostasyConfig::default();
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };

    eprintln!(
        "c1_phase_2_track_d Stage A seed_diversity_gallery — cycle_300 at seeds {seeds:?}"
    );
    eprintln!("  Per-seed final-state stats (forward signal toward §7.2 cross-track gate):");

    for &seed in seeds.iter() {
        let mut state = init_c1_state_phase_2_r7(GRID_SIZE, seed, &init_params);
        let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        run_with_closures(
            &mut state,
            &mut kinematics,
            &config,
            &closures,
            |_, _| {},
        );

        let total = GRID_SIZE * GRID_SIZE;
        let mut continental = 0;
        for j in 0..GRID_SIZE {
            for i in 0..GRID_SIZE {
                if matches!(
                    state.plate_type.get(i, j),
                    ymir_core::tectonics_v2::boundaries::plate_type::PlateType::Continental
                ) {
                    continental += 1;
                }
            }
        }
        let plates_remaining = {
            let mut seen = std::collections::HashSet::new();
            for &pid in state.plate_id.data() {
                seen.insert(pid);
            }
            seen.len()
        };
        let new_plates_count =
            kinematics.velocities.len().saturating_sub(state.num_plates);

        eprintln!(
            "    seed = {seed:>5}  continental = {continental:>4} / {total}  plates_remaining = {plates_remaining:>3}  new_plate_ids (rift) = {new_plates_count}"
        );

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
    eprintln!(
        "  files = 6 PNGs (seed_NNNNN_altitude.png + seed_NNNNN_s.png × 3 seeds)"
    );
}

fn print_stats(
    tag: &str,
    state: &C1State,
    iso_config: &IsostasyConfig,
    closures: &C1Closures,
) {
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

    let clip_low = a_min < -ALTITUDE_PALETTE_HALF_RANGE;
    let clip_high = a_max > ALTITUDE_PALETTE_HALF_RANGE;
    let clip_marker = match (clip_low, clip_high) {
        (false, false) => "",
        (true, false) => "  [clip-low]",
        (false, true) => "  [clip-high]",
        (true, true) => "  [clip both]",
    };

    eprintln!(
        "    cycle_{tag}  S̃ min={s_min:.4} mean={s_mean:.4} max={s_max:.4}   \
         altitude (post-S-S) min={a_min:.4} mean={a_mean:.4} max={a_max:.4}{clip_marker}"
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
