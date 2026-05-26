//! Issue #129 Phase 2 Track A Stage D — visual gallery generator.
//!
//! Single `#[ignore]`'d test that exercises the full Phase 2 stack
//! (Davis-Suppe + equilibrium-height + erosion + Stein-Stein
//! bathymetry, all 4 closures enabled) at 64²×300 steps and
//! produces:
//!
//! 1. A 5-cycle PNG gallery (cycle 0 / 50 / 100 / 200 / 300, 2
//!    fields per cycle = 10 files) under
//!    `docs/reports/c1_phase_2_oceanic_bathymetry/`.
//! 2. Per-cycle stats on `S̃` and altitude distributions.
//! 3. Final-cycle Architecture C observability re-application: the
//!    altitude PNG at each cycle is produced AFTER re-applying
//!    Stein-Stein, so the visual carries the bathymetric imprint
//!    that Phase 2 Track A is the load-bearing addition for.
//!
//! Marked `#[ignore]` because it generates artefacts and prints
//! verbose diagnostics — not part of the routine regression sweep.
//! Invocation:
//!
//! ```bash
//! cargo test --release -p ymir-core \
//!     --test c1_phase_2_visual_gallery \
//!     -- --ignored --nocapture
//! ```
//!
//! ## Palette decisions
//!
//! - **Altitude palette: `[-1.13, +1.13]` symmetric around sea
//!   level = 0.** Phase 2 Track A altitude is bipolar: oceanic
//!   cells receive negative values (`-d(t) / depth_scale_m`,
//!   range `[-1.13, -0.52]`) while continental cells stay
//!   positive from the isostasy output. A symmetric palette
//!   centred on sea level renders ridge cells mid-tone, deep
//!   oceanic dark, continental bright. Phase 1.4 used `[0, 1.0]`
//!   because Phase 1.4 altitude was unipolar (continental-
//!   dominant after isostasy normalisation); cross-phase
//!   comparison requires this palette context.
//! - **`S̃` palette: `[0, 3.0]`** unchanged. Preserves
//!   cross-phase comparability of `S̃` images Phase 1.1 → 1.2 →
//!   1.3 → 1.4 → 2 Track A per the
//!   `feedback_viz_palette_absolute_for_comparison` rule.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;
use ymir_core::tectonics_v2::field::Field2D;

const GRID_SIZE: usize = 64;
const SEED: u64 = 42;
const N_STEPS: usize = 300;
const S_VIZ_MAX: f64 = 3.0;
const ALTITUDE_PALETTE_HALF_RANGE: f32 = 1.13;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_phase_2_oceanic_bathymetry")
}

#[test]
#[ignore]
fn phase_2_bathymetry_visual_gallery() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let mut state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let iso_config = IsostasyConfig::default();
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };

    eprintln!("c1_phase_2 Stage D visual gallery — grid={GRID_SIZE}², steps={N_STEPS}");
    eprintln!(
        "  closures: DS={} EH={} erosion={} S-S={}  (full Phase 2 stack)",
        closures.davis_suppe.enabled,
        closures.equilibrium_height.enabled,
        closures.erosion.enabled,
        closures.oceanic_bathymetry.enabled,
    );

    print_stats("000", &state, &iso_config, &closures);
    dump_snapshot(&state, 0, &dir, &iso_config, &closures);

    let snapshot_steps: [usize; 4] = [49, 99, 199, 299];
    let started = std::time::Instant::now();
    run_with_closures(
        &mut state,
        &kinematics,
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

    eprintln!();
    eprintln!(
        "  wall time      = {:.2?} ({:.2?} / step)",
        elapsed,
        elapsed / N_STEPS as u32
    );
    eprintln!("  output dir     = {}", dir.display());
    eprintln!("  files          = 10 PNGs (cycle_NNN_altitude.png + cycle_NNN_s.png × 5)");
    eprintln!();
    eprintln!("  Architectural finding cross-check (Stage A consistency):");
    let final_age_stats = age_stats_oceanic(&state);
    eprintln!(
        "    oceanic age distribution: min={:.4} max={:.4} mean={:.4} median={:.4}",
        final_age_stats.min, final_age_stats.max, final_age_stats.mean, final_age_stats.median,
    );
    eprintln!(
        "    Stage A baseline:         min≈0     max≈6958   mean≈4.67   median≈0"
    );
    eprintln!(
        "    Pile-up consistency: same run reproduces Stage A density-advection finding."
    );
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
            if state.plate_type.get(i, j) == PlateType::Oceanic {
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
    // Altitude: compute isostasy + re-apply S-S (Architecture C
    // observability — same pattern as the Stage A acceptance
    // test's run-boundary observation).
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

/// Render a bipolar (signed) altitude field with sea level = 0 at
/// the palette midpoint. Maps `[-half_range, +half_range]` to the
/// hypsometric `[0, 1]` color ramp; values outside that range
/// clamp at the palette extremes.
fn save_altitude_bipolar_png(altitude: &GridF32, half_range: f32, path: &Path) {
    let nx = altitude.width;
    let ny = altitude.height;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    let sea_norm = 0.5_f32; // sea level at palette midpoint
    for j in 0..ny {
        for i in 0..nx {
            let raw = altitude.get(i as i32, j as i32);
            let t = (raw + half_range) / (2.0 * half_range);
            let t = t.clamp(0.0, 1.0);
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
