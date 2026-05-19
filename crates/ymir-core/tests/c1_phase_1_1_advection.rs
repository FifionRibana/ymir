//! Issue #120 — C1 Phase 1.1 advection-only sanity integration test.
//!
//! Runs the C1 prototype at 64² for 300 forward-Euler steps, dumps
//! altitude PNGs at cycles 0, 50, 100, 200, 300 through the existing
//! `tectonics::isostasy::compute_isostasy` heightmap converter, and
//! asserts mass conservation < 1e-6 drift.
//!
//! The visual outputs land in
//! `docs/reports/c1_phase_1_1_advection/` for human review. Phase 1.1
//! is a transport-correctness check, **not** a plausible-continent
//! check: the PNGs are expected to show convergence-zone thickening,
//! divergence-zone thinning, rigid cratonic-cell transport, and
//! nothing more (no orogenic structure, no isostatic balance, no
//! erosion — those land in Phase 1.2-1.4).
//!
//! ## PNG renderer
//!
//! Self-contained inline helper (~50 lines) using the `image` crate
//! that's already a dependency of ymir-core. The viz library
//! (`ymir-viz/src/visualization/v2_viz.rs`) ships full hypsometric
//! rendering but is gated under the `v2_legacy` feature, so we
//! bring a minimal renderer here. C1 Phase 1.4 (UI + production)
//! will reintroduce a paradigm-agnostic viz that this test can
//! migrate to.

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::init::init_c1_state_phase_1_1;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_advection_only, C1TimeLoopConfig};
use ymir_core::tectonics_v2::field::Field2D;

const GRID_SIZE: usize = 64;
const SEED: u64 = 42;
const N_STEPS: usize = 300;

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_phase_1_1_advection")
}

#[test]
fn c1_phase_1_1_advection_sanity() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");

    let mut state = init_c1_state_phase_1_1(GRID_SIZE, SEED);
    let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
    let config = C1TimeLoopConfig {
        n_steps: N_STEPS,
        dx: 1.0 / GRID_SIZE as f64,
        dy: 1.0 / GRID_SIZE as f64,
    };

    let initial_mass: f64 = state.s.data().iter().sum();
    eprintln!(
        "c1_phase_1_1: grid={}², steps={}, max|v|={:.4e}, dt={:.4e}",
        GRID_SIZE,
        N_STEPS,
        kinematics.max_velocity(),
        0.5 * config.dx / kinematics.max_velocity()
    );
    eprintln!("c1_phase_1_1: initial mass = {:.6}", initial_mass);

    // Snapshot the initial state before the loop runs (cycle 0).
    dump_snapshot(&state, 0, &dir);

    // Snapshot steps fire on the post-step-N state. Step indices
    // are 0-based inside the callback so e.g. step=49 means "after
    // 50 forward-Euler updates".
    let snapshot_steps: [usize; 4] = [49, 99, 199, 299];

    let started = std::time::Instant::now();
    print_s_stats("000", &state);
    run_advection_only(&mut state, &kinematics, &config, |step, current_state| {
        // NaN guard.
        assert!(
            current_state.s.data().iter().all(|v| v.is_finite()),
            "non-finite S̃ at step {}",
            step + 1
        );
        if snapshot_steps.contains(&step) {
            print_s_stats(&format!("{:03}", step + 1), current_state);
            dump_snapshot(current_state, step + 1, &dir);
        }
    });
    let elapsed = started.elapsed();

    let final_mass: f64 = state.s.data().iter().sum();
    let drift = (final_mass - initial_mass).abs() / initial_mass;
    eprintln!("c1_phase_1_1: final mass   = {:.6}", final_mass);
    eprintln!("c1_phase_1_1: mass drift    = {:.3e} (threshold 1e-6)", drift);
    eprintln!("c1_phase_1_1: wall time     = {:.2?} ({:.2?} per step)",
        elapsed,
        elapsed / N_STEPS as u32);
    eprintln!("c1_phase_1_1: output dir    = {}", dir.display());

    assert!(
        drift < 1e-6,
        "mass conservation drift {:.3e} exceeds 1e-6 threshold over {} steps",
        drift,
        N_STEPS
    );
}

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
        "c1_phase_1_1: cycle_{tag} S̃ min={min:.4} mean={mean:.4} max={max:.4} std={std:.4e}"
    );
}

/// Absolute palette bound for the fixed-scale `S̃` snapshot. Initial
/// `S̃` ranges in `[0.2, 1.0]` (oceanic / continental). Phase 1.1
/// transport accumulates mass at convergence zones with no closure
/// to bound it; cells above `S_VIZ_MAX` saturate at the top of the
/// palette, signalling "extreme pile-up that a closure will need to
/// absorb in Phase 1.2+".
const S_VIZ_MAX: f64 = 2.0;

/// Dump a snapshot pair for the given state:
///
/// - `cycle_NNN_altitude.png` — Airy-isostasy heightmap through
///   `compute_isostasy`, rendered with a 4-stop hypsometric palette.
///   Useful for "what the downstream Phase B pipeline sees" but uses
///   the v2 isostasy's **per-frame** rescale, so the visual is
///   palette-relative not absolute.
/// - `cycle_NNN_s.png` — direct `S̃` render with a **fixed** absolute
///   palette `[0, S_VIZ_MAX]` (palette discipline per the
///   `viz-palette-absolute-for-inter-run-comparison` memory entry).
///   This is the transport-correctness diagnostic: cells above the
///   initial continental level (1.0) glow brown / saturate white at
///   convergence zones; oceanic regions stay deep blue.
fn dump_snapshot(state: &C1State, cycle: usize, dir: &Path) {
    let iso = compute_isostasy(&state.s, &IsostasyConfig::default());
    let alt_path = dir.join(format!("cycle_{:03}_altitude.png", cycle));
    save_hypsometric_png(&iso.heightmap, iso.sea_level_normalized, &alt_path);

    let s_path = dir.join(format!("cycle_{:03}_s.png", cycle));
    save_s_fixed_palette_png(&state.s, S_VIZ_MAX, &s_path);
}

/// Direct `S̃` PNG with the same hypsometric palette family as
/// altitude, but normalised against a fixed `S_VIZ_MAX` rather than
/// a per-frame `[min, max]` rescale.
fn save_s_fixed_palette_png(s: &Field2D, s_max: f64, path: &Path) {
    let nx = s.nx();
    let ny = s.ny();
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    // Sea level marker at the oceanic init value (`OCEANIC_VALUE_DEFAULT = 0.2`)
    // — gives the palette a reference shore line. Continental cells
    // (S̃ around 1.0) and accumulated pile-ups (> 1.0) read as land.
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

/// Render a `GridF32` heightmap in `[0, 1]` to a PNG using a
/// 4-stop hypsometric palette: deep blue → light blue at sea
/// level, then green → brown → white above.
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
            // Y-flip: image row 0 maps to grid row (ny - 1 - j). The
            // underlying Field2D convention has j=0 at the bottom;
            // PNG row 0 sits at the top — flip to keep "north up".
            let img_row = (ny - 1 - j) as u32;
            img.put_pixel(i as u32, img_row, Rgb(rgb));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create parent dir");
    }
    img.save(path).expect("save PNG");
}

/// Simple 4-stop hypsometric palette.
///
/// - `[0, sea_norm/2]`         deep blue → mid blue
/// - `[sea_norm/2, sea_norm]`  mid blue → light blue (sea level)
/// - `[sea_norm, mid]`         green coast → brown mid-elevation
/// - `[mid, 1.0]`              brown → white peak
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

// Silence the unused-import warning on the Field2D pull when only
// the symbol is needed inside test helpers.
#[allow(dead_code)]
fn _field2d_unused_witness(f: &Field2D) -> usize {
    f.nx() * f.ny()
}
