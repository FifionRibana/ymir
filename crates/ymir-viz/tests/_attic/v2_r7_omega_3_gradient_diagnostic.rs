//! Step 12 R7 ω.3 D — gradient ∇S̃ distribution diagnostic.
//!
//! Builds the init S̃ field for the three R7.A.2.4 init modes
//! (Radial / Orogenic σ=0.10 / Composite) at the exact R7.A.2.4
//! configuration (64² × seed=42 × num_plates=8 × continental_ratio=0.3
//! × amplitude=0), computes the Sobel-periodic gradient magnitude
//! field per mode, and dumps:
//!
//! - `summary.md` — comparison table across the three modes
//! - `<mode>/stats.md` — full distribution stats per mode
//! - `<mode>/slope.png` — slope render with a **shared palette**
//!   `[0, global_max]` across the three modes (per
//!   `feedback_viz_palette_absolute_for_comparison`)
//!
//! Output root: `docs/reports/step12_r7_omega3_gradient_diagnostic/`
//!
//! Run with:
//! ```bash
//! cargo test --release -p ymir-viz \
//!   --test v2_r7_omega_3_gradient_diagnostic \
//!   -- --ignored --nocapture
//! ```
//!
//! No production code is touched; this is a pure diagnostic test.

use std::path::PathBuf;

use image::{ImageBuffer, Rgb};

use ymir_core::tectonics_v2::init::{init_s_field, InitContext, PlateInitData};
use ymir_core::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

use ymir_viz::bridge::v2::V2InitModeSpec;

const NX: usize = 64;
const NY: usize = 64;
const SEED: u64 = 42;
const NUM_PLATES: usize = 8;
const CONT_RATIO: f64 = 0.3;

fn out_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step12_r7_omega3_gradient_diagnostic")
}

fn build_init_s(mode_spec: V2InitModeSpec) -> Vec<f64> {
    let cfg = VoronoiConfig {
        num_plates: NUM_PLATES,
        continental_ratio: CONT_RATIO,
    };
    let plates = generate_voronoi(NX, NY, &cfg, SEED);
    let plate_data = PlateInitData {
        plate_id: &plates.plate_id,
        plate_type: &plates.plate_type,
        seed_coords: Some(&plates.seed_coords),
    };
    let ctx = InitContext {
        nx: NX,
        ny: NY,
        seed: SEED,
        amplitude: 0.0,
        plate_data: Some(plate_data),
    };
    init_s_field(mode_spec.into_core(), &ctx).data().to_vec()
}

/// Sobel-periodic per-cell gradient magnitude. Local copy of the
/// (private) `sobel_magnitude_periodic` in `v2_viz` so this test
/// keeps zero production touch.
fn sobel_periodic(field: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let mut out = vec![0.0_f64; nx * ny];
    for j in 0..ny {
        let jp = (j + 1) % ny;
        let jm = (j + ny - 1) % ny;
        for i in 0..nx {
            let ip = (i + 1) % nx;
            let im = (i + nx - 1) % nx;
            let s_jm_im = field[jm * nx + im];
            let s_jm_ip = field[jm * nx + ip];
            let s_jm_i = field[jm * nx + i];
            let s_jp_im = field[jp * nx + im];
            let s_jp_ip = field[jp * nx + ip];
            let s_jp_i = field[jp * nx + i];
            let s_j_im = field[j * nx + im];
            let s_j_ip = field[j * nx + ip];
            let gx = (-s_jm_im + s_jm_ip)
                + 2.0 * (-s_j_im + s_j_ip)
                + (-s_jp_im + s_jp_ip);
            let gy = (-s_jm_im - 2.0 * s_jm_i - s_jm_ip)
                + (s_jp_im + 2.0 * s_jp_i + s_jp_ip);
            out[j * nx + i] = (gx * gx + gy * gy).sqrt() / 8.0;
        }
    }
    out
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((p / 100.0) * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

#[derive(Clone, Copy, Debug)]
struct GradStats {
    min: f64,
    mean: f64,
    max: f64,
    p50: f64,
    p90: f64,
    p95: f64,
    p99: f64,
    frac_above_0_01: f64,
    frac_above_0_05: f64,
    frac_above_0_10: f64,
    frac_above_0_20: f64,
}

fn compute_stats(grad: &[f64]) -> GradStats {
    let n = grad.len() as f64;
    let mut sorted: Vec<f64> = grad.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let sum: f64 = grad.iter().sum();
    let count_above = |t: f64| grad.iter().filter(|&&v| v > t).count() as f64 / n;
    GradStats {
        min: sorted[0],
        mean: sum / n,
        max: sorted[sorted.len() - 1],
        p50: percentile(&sorted, 50.0),
        p90: percentile(&sorted, 90.0),
        p95: percentile(&sorted, 95.0),
        p99: percentile(&sorted, 99.0),
        frac_above_0_01: count_above(0.01),
        frac_above_0_05: count_above(0.05),
        frac_above_0_10: count_above(0.10),
        frac_above_0_20: count_above(0.20),
    }
}

/// Simple "hot" colormap (black → red → yellow → white) with linear
/// remap from `[vmin, vmax]` to `[0, 1]`. PNG row 0 maps to grid row
/// `ny - 1` (Y-flip, matches the convention in `v2_viz::field_to_rgba`).
fn save_slope_png(grad: &[f64], path: &PathBuf, vmin: f64, vmax: f64) {
    let nx = NX;
    let ny = NY;
    let range = (vmax - vmin).max(1e-12);
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    for j in 0..ny {
        for i in 0..nx {
            let v = grad[j * nx + i];
            let t = ((v - vmin) / range).clamp(0.0, 1.0);
            let r = (t * 3.0).min(1.0);
            let g = ((t - 1.0 / 3.0) * 3.0).clamp(0.0, 1.0);
            let b = ((t - 2.0 / 3.0) * 3.0).clamp(0.0, 1.0);
            let img_row = ny - 1 - j;
            img.put_pixel(
                i as u32,
                img_row as u32,
                Rgb([
                    (r * 255.0).round() as u8,
                    (g * 255.0).round() as u8,
                    (b * 255.0).round() as u8,
                ]),
            );
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create dir");
    }
    img.save(path).expect("write PNG");
}

fn write_stats_md(path: &PathBuf, label: &str, stats: &GradStats, vmax_shared: f64) {
    let body = format!(
        "# R7 ω.3 D — {label} init |∇S̃| stats\n\
\n\
64² × seed=42 × num_plates=8 × continental_ratio=0.3 × amplitude=0\n\
\n\
Sobel-periodic per-cell magnitude (no Δx normalization — values are\n\
in units of \"δS̃ across one cell\"). PNG palette shared across the\n\
three modes: `[0, {vmax_shared:.4e}]`.\n\
\n\
## Distribution\n\
\n\
| metric | value |\n\
|---|---|\n\
| min  | {:.4e} |\n\
| mean | {:.4e} |\n\
| max  | {:.4e} |\n\
| p50  | {:.4e} |\n\
| p90  | {:.4e} |\n\
| p95  | {:.4e} |\n\
| p99  | {:.4e} |\n\
\n\
## Fractions above threshold\n\
\n\
| threshold | fraction of cells |\n\
|---|---|\n\
| > 0.01 | {:.4} |\n\
| > 0.05 | {:.4} |\n\
| > 0.10 | {:.4} |\n\
| > 0.20 | {:.4} |\n",
        stats.min,
        stats.mean,
        stats.max,
        stats.p50,
        stats.p90,
        stats.p95,
        stats.p99,
        stats.frac_above_0_01,
        stats.frac_above_0_05,
        stats.frac_above_0_10,
        stats.frac_above_0_20,
    );
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("create dir");
    }
    std::fs::write(path, body).expect("write stats.md");
}

#[test]
#[ignore]
fn r7_omega_3_d_gradient_diagnostic_three_modes() {
    // R7.A.2.4 config: Run A = radial_profile_default;
    // Run B = orogenic_default with width_sigma_ratio overridden to 0.10;
    // Run C = composite_default (which already bakes in σ=0.10).
    let mut orogenic_sigma_10 = V2InitModeSpec::orogenic_default();
    if let V2InitModeSpec::Orogenic { width_sigma_ratio, .. } = &mut orogenic_sigma_10 {
        *width_sigma_ratio = 0.10;
    }

    let modes: Vec<(&str, V2InitModeSpec)> = vec![
        ("radial", V2InitModeSpec::radial_profile_default()),
        ("orogenic_sigma_10", orogenic_sigma_10),
        ("composite", V2InitModeSpec::composite_default()),
    ];

    // Pass 1 — build init S̃, compute slope, gather stats + global max
    let mut grads: Vec<(&str, Vec<f64>)> = Vec::with_capacity(modes.len());
    let mut all_stats: Vec<(&str, GradStats)> = Vec::with_capacity(modes.len());
    let mut global_max: f64 = 0.0;
    for (label, spec) in &modes {
        let s = build_init_s(*spec);
        let grad = sobel_periodic(&s, NX, NY);
        let stats = compute_stats(&grad);
        global_max = global_max.max(stats.max);
        eprintln!(
            "{label:>18}  max={:.4e}  p90={:.4e}  p99={:.4e}  frac>0.05={:.4}  frac>0.10={:.4}",
            stats.max, stats.p90, stats.p99, stats.frac_above_0_05, stats.frac_above_0_10,
        );
        grads.push((label, grad));
        all_stats.push((label, stats));
    }

    // Pass 2 — render PNGs with shared bounds + per-mode stats.md
    let root = out_root();
    for ((label, grad), (_, stats)) in grads.iter().zip(all_stats.iter()) {
        let dir = root.join(label);
        save_slope_png(grad, &dir.join("slope.png"), 0.0, global_max);
        write_stats_md(&dir.join("stats.md"), label, stats, global_max);
    }

    // Summary table across the three modes
    let mut summary = format!(
        "# R7 ω.3 D — gradient ∇S̃ distribution comparison\n\
\n\
64² × seed=42 × num_plates=8 × continental_ratio=0.3 × amplitude=0\n\
\n\
Sobel-periodic per-cell magnitude. Shared palette for all\n\
`slope.png` files: `[0, {global_max:.4e}]`.\n\
\n\
## Summary table\n\
\n\
| mode | min | mean | max | p90 | p99 | frac>0.05 | frac>0.10 | frac>0.20 |\n\
|---|---|---|---|---|---|---|---|---|\n",
    );
    for (label, stats) in &all_stats {
        summary.push_str(&format!(
            "| {} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3} | {:.3} | {:.3} |\n",
            label,
            stats.min,
            stats.mean,
            stats.max,
            stats.p90,
            stats.p99,
            stats.frac_above_0_05,
            stats.frac_above_0_10,
            stats.frac_above_0_20,
        ));
    }
    std::fs::create_dir_all(&root).expect("create root");
    std::fs::write(root.join("summary.md"), summary).expect("write summary.md");

    eprintln!("\nGlobal max |∇S̃| across 3 modes: {:.4e}", global_max);
    eprintln!("Output root: {}", root.display());
}
