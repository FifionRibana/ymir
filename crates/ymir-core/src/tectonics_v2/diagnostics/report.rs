//! Markdown writer for Step 0's baseline report.
//!
//! The report header explicitly disclaims any comparison — Step 0 is
//! the anchor, not a regression. Subsequent steps' PRs introduce
//! diff-against-previous logic.

use std::io::Write;
use std::path::Path;

use super::metrics::{Metrics, SolverConfigDump};
use crate::tectonics_v2::scales::Scales;

/// Write a full Step 0 baseline report.
pub fn write_markdown_report(
    path: &Path,
    seed: u64,
    scales: &Scales,
    configs: &[SolverConfigDump],
    metrics: &[Metrics],
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = std::fs::File::create(path)?;
    write!(f, "{}", build_markdown(seed, scales, configs, metrics))
}

pub fn build_markdown(
    seed: u64,
    scales: &Scales,
    configs: &[SolverConfigDump],
    metrics: &[Metrics],
) -> String {
    let mut out = String::new();
    out.push_str("# Step 0 — Nondim Stokes core + S advection (baseline)\n\n");
    out.push_str("> **Baseline reference for milestone \"Solver reconstruction\".**\n");
    out.push_str("> Do NOT compare against earlier reports — this is the first report.\n");
    out.push_str("> Subsequent steps' reports will diff against this one.\n\n");
    out.push_str(&format!("- Seed: `{}`\n", seed));
    out.push_str("- Entry-condition decisions archived in `tectonics_v2/README.md`.\n\n");

    out.push_str("## Physical scales\n\n");
    out.push_str("```\n");
    out.push_str(&scales.report());
    out.push_str("\n```\n\n");

    for (cfg, m) in configs.iter().zip(metrics.iter()) {
        out.push_str(&format!("## Grid {}×{}\n\n", m.grid_nx, m.grid_ny));
        out.push_str(&cfg.render_markdown());
        out.push_str("\n### Timing\n\n");
        out.push_str(&format!(
            "- wallclock total: `{:.3} s`\n- wallclock per step (mean): `{:.3} ms`\n- steps: `{}`\n\n",
            m.wallclock_total.as_secs_f64(),
            m.wallclock_per_step_mean.as_secs_f64() * 1.0e3,
            m.steps,
        ));

        out.push_str("### Solver health\n\n");
        if m.kappa_estimate.is_finite() {
            out.push_str(&format!(
                "- κ(A) estimate from outer CG iterations: `{:.2e}` (active metric)\n",
                m.kappa_estimate,
            ));
        } else {
            out.push_str(
                "- κ(A) estimate: N/A — outer CG converged in 0 iterations (the Kolmogorov-like placeholder forcing produces an exactly divergence-free velocity from A⁻¹f, so the Schur complement problem is trivially satisfied by p=0). The framework slot is exercised; real κ estimates come online at Step 2 when GPE spreading makes the Schur-complement nontrivial.\n",
            );
        }
        out.push_str(&format!(
            "- effective η_max/η_min over run: `{:.3}` (placeholder; trivially 1.0 at Step 0)\n",
            m.eta_contrast,
        ));
        out.push_str(&format!(
            "- outer CG iterations — mean: `{:.1}`, max: `{}`\n",
            m.outer_iter_mean, m.outer_iter_max,
        ));
        out.push_str("- outer CG iteration histogram (5 bins):\n\n");
        let hist = &m.outer_iter_histogram;
        out.push_str("  | bin ≤ | count |\n  |---|---|\n");
        for b in 0..5 {
            out.push_str(&format!("  | {} | {} |\n", hist.bin_edges[b], hist.counts[b]));
        }
        out.push('\n');
        out.push_str(&format!(
            "- inner CG iterations (per inner solve) — mean: `{:.1}`, max: `{}`\n\n",
            m.inner_iter_mean, m.inner_iter_max,
        ));

        out.push_str("### Mass conservation of S\n\n");
        out.push_str(&format!(
            "- initial mass: `{:.9e}`\n- final mass: `{:.9e}`\n- relative drift: `{:.3e}`\n\n",
            m.mass_s_initial, m.mass_s_final, m.mass_drift_relative,
        ));

        out.push_str("### Null-space health (post-solve means)\n\n");
        out.push_str(&format!(
            "- max |mean(P)| across solves: `{:.3e}`\n- max |mean(vx)|: `{:.3e}`\n- max |mean(vy)|: `{:.3e}`\n\n",
            m.max_abs_mean_p, m.max_abs_mean_vx, m.max_abs_mean_vy,
        ));

        out.push_str("### Velocity magnitude\n\n");
        out.push_str(&format!("- peak |v|: `{:.3e}`\n\n", m.vmax_peak));

        out.push_str("### Heightmaps of S\n\n");
        if m.heightmap_paths.is_empty() {
            out.push_str("- (none recorded)\n\n");
        } else {
            for path in &m.heightmap_paths {
                out.push_str(&format!("- `{}`\n", path));
            }
            out.push('\n');
        }

        out.push_str("### Dormant metrics (inactive at Step 0)\n\n");
        out.push_str("| metric | activated at |\n|---|---|\n");
        out.push_str("| S̃_eq (active-orogen mean thickness) | Step 5+ |\n");
        out.push_str("| boundary type diversity | Step 5 |\n");
        out.push_str("| yielding cell fraction | Step 3 |\n");
        out.push_str("| cratonic stability | Step 9 |\n");
        out.push_str("| Newton outcome distribution | Step 1 |\n");
        out.push_str("| age field stats | Step 10 |\n\n");
    }

    out.push_str("---\n");
    out.push_str("*Generated by `cargo run --release --bin step_baseline`.*\n");
    out
}

#[cfg(test)]
mod tests {
    use super::super::metrics::{IterationHistogram, Metrics};
    use super::*;
    use std::time::Duration;

    fn fake_config() -> SolverConfigDump {
        SolverConfigDump {
            discretization: "MAC staggered".into(),
            harmonic_averaging: "on".into(),
            preconditioner: "block-diag Jacobi + null-space proj".into(),
            gauge_fixing: "mean(P), mean(vx), mean(vy)".into(),
            outer_tol: 1e-8,
            inner_tol: 1e-10,
            outer_max_iter: 200,
            inner_max_iter: 500,
            cfl_factor: 0.3,
            grid_spacing_nondim: 1.0 / 64.0,
            body_force: "SinusoidalForce ε=0.1".into(),
            seed: 42,
        }
    }

    fn fake_metrics() -> Metrics {
        let mut m = Metrics::empty(64, 64, 300);
        m.wallclock_total = Duration::from_millis(1234);
        m.wallclock_per_step_mean = Duration::from_micros(4113);
        m.kappa_estimate = 1.2e4;
        m.outer_iter_mean = 12.0;
        m.outer_iter_max = 18;
        m.outer_iter_histogram = IterationHistogram::from_samples(&[10, 12, 14, 18, 12]);
        m.inner_iter_mean = 30.0;
        m.inner_iter_max = 55;
        m.mass_s_initial = 1024.0;
        m.mass_s_final = 1024.0;
        m.mass_drift_relative = 3.2e-15;
        m.max_abs_mean_p = 1e-13;
        m.max_abs_mean_vx = 5e-14;
        m.max_abs_mean_vy = 4e-14;
        m.vmax_peak = 0.034;
        m
    }

    #[test]
    fn markdown_contains_header_banner() {
        let s = build_markdown(42, &Scales::default(), &[fake_config()], &[fake_metrics()]);
        assert!(s.contains("Baseline reference"));
        assert!(s.contains("Do NOT compare"));
        assert!(s.contains("Solver configuration"));
        assert!(s.contains("κ(A) estimate"));
        assert!(s.contains("Dormant metrics"));
    }
}
