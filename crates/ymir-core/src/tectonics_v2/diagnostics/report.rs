//! Markdown writer for Step 1's baseline report.
//!
//! Compared to Step 0, the report now:
//! - identifies as Step 1 and references the Step 0 anchor;
//! - adds a Newton-specific section (outcome distribution, iters,
//!   cap-activation fractions during ramp and steady state,
//!   continuation outcome, effective `η_max/η_min`);
//! - emits a "Comparison vs Step 0" block when a reference report
//!   is supplied, with the CG-iter ratio classified in paliers
//!   (idéal / acceptable / suspect / fail) per the spec.

use std::io::Write;
use std::path::Path;

use super::comparison::{render_grid_comparison, StepReference};
use super::metrics::{Metrics, SolverConfigDump};
use super::mms_bench::MmsResults;
use crate::tectonics_v2::scales::Scales;

pub struct ReportInputs<'a> {
    pub seed: u64,
    pub scales: &'a Scales,
    pub configs: &'a [SolverConfigDump],
    pub metrics: &'a [Metrics],
    pub previous: Option<&'a StepReference>,
    /// One justification string per grid entry, for suspect-tier
    /// CG-iter ratios. Pass `""` to leave empty.
    pub suspect_justifications: &'a [String],
    /// Discretisation validation results (MMS slopes + Newton tail).
    pub mms: Option<&'a MmsResults>,
}

pub fn write_markdown_report(
    path: &Path,
    inputs: &ReportInputs,
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = std::fs::File::create(path)?;
    write!(f, "{}", build_markdown(inputs))
}

pub fn build_markdown(inputs: &ReportInputs) -> String {
    let mut out = String::new();
    out.push_str("# Step 1 — Power-law rheology + Newton solver (baseline)\n\n");
    out.push_str("> **Step 1 reference run for milestone \"Solver reconstruction\".**\n");
    out.push_str("> Compared against Step 0 (`docs/reports/step0_report.md`).\n");
    out.push_str("> Subsequent steps' reports will diff against this one.\n\n");
    out.push_str(&format!("- Seed: `{}`\n", inputs.seed));
    out.push_str("- Formulation: thin viscous sheet + power-law rheology (n > 1). Linear internal system per Newton iteration is symmetric (variational structure); CG suffices.\n\n");

    out.push_str("## Physical scales\n\n");
    out.push_str("```\n");
    out.push_str(&inputs.scales.report());
    out.push_str("\n```\n\n");

    if let Some(mms) = inputs.mms {
        out.push_str(&super::mms_bench::render_markdown(mms));
    }

    for (idx, (cfg, m)) in inputs.configs.iter().zip(inputs.metrics.iter()).enumerate() {
        out.push_str(&format!("## Grid {}×{}\n\n", m.grid_nx, m.grid_ny));
        out.push_str(&cfg.render_markdown());

        out.push_str("\n### Timing\n\n");
        out.push_str(&format!(
            "- wallclock total: `{:.3} s`\n- wallclock per step (mean): `{:.3} ms`\n- steps: `{}`\n\n",
            m.wallclock_total.as_secs_f64(),
            m.wallclock_per_step_mean.as_secs_f64() * 1.0e3,
            m.steps,
        ));

        out.push_str("### Linear-solver health (CG inside Newton)\n\n");
        if m.kappa_estimate.is_finite() {
            out.push_str(&format!(
                "- κ(A) estimate from CG iterations (per Newton step): `{:.2e}`\n",
                m.kappa_estimate,
            ));
        } else {
            out.push_str("- κ(A) estimate: N/A (CG converged on the initial guess)\n");
        }
        out.push_str(&format!(
            "- CG iterations per Newton step — mean: `{:.1}`, max: `{}`\n",
            m.cg_iter_mean, m.cg_iter_max,
        ));
        out.push_str("- CG iteration histogram (5 bins):\n\n");
        let hist = &m.cg_iter_histogram;
        out.push_str("  | bin ≤ | count |\n  |---|---|\n");
        for b in 0..5 {
            out.push_str(&format!("  | {} | {} |\n", hist.bin_edges[b], hist.counts[b]));
        }
        out.push('\n');

        if let Some(na) = &m.newton {
            out.push_str("### Newton (nonlinear) health\n\n");
            let (pc, ps, pd, pcap) = na.outcome_percentages();
            out.push_str(&format!(
                "- outcome distribution — Converged: `{:.1}%`, Stalled: `{:.1}%`, Diverged: `{:.1}%`, CappedIters: `{:.1}%`\n",
                pc, ps, pd, pcap,
            ));
            out.push_str(&format!(
                "- Newton outer iters per timestep — mean: `{:.1}`, max: `{}`\n",
                na.outer_iters_mean(),
                na.outer_iters_max(),
            ));
            out.push_str(&format!(
                "- effective η_max/η_min over run — mean: `{:.2}`, max: `{:.2}`\n",
                na.eta_contrast_mean(),
                na.eta_contrast_max(),
            ));
            out.push_str(&format!(
                "- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `{:.3}%`; steady state: `{:.3}%` (spec target < 1%)\n",
                100.0 * na.cap_fraction_ramp_max,
                100.0 * na.cap_fraction_steady_max,
            ));
            match na.continuation_all_converged {
                Some(true) => out.push_str(&format!(
                    "- continuation ramp: ✅ all {} sub-solves converged\n",
                    na.continuation_iters_used,
                )),
                Some(false) => out.push_str(&format!(
                    "- continuation ramp: ❌ failed after {} sub-solves\n",
                    na.continuation_iters_used,
                )),
                None => out.push_str("- continuation ramp: not run\n"),
            }
            out.push('\n');
        }

        out.push_str("### Mass conservation of S\n\n");
        out.push_str(&format!(
            "- initial mass: `{:.9e}`\n- final mass: `{:.9e}`\n- relative drift: `{:.3e}`\n\n",
            m.mass_s_initial, m.mass_s_final, m.mass_drift_relative,
        ));

        out.push_str("### Null-space health\n\n");
        out.push_str(&format!(
            "- max |mean(vx)| across solves: `{:.3e}`\n- max |mean(vy)|: `{:.3e}`\n\n",
            m.max_abs_mean_vx, m.max_abs_mean_vy,
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

        // --- Comparison against the previous step, if supplied ---
        if let Some(prev) = inputs.previous {
            if let Some(prev_grid) = prev
                .grids
                .iter()
                .find(|g| g.grid == (m.grid_nx, m.grid_ny))
            {
                out.push_str("### Comparison with Step 0\n\n");
                let justification = inputs
                    .suspect_justifications
                    .get(idx)
                    .map(|s| s.as_str())
                    .unwrap_or("");
                let justif = if justification.is_empty() { None } else { Some(justification) };
                out.push_str(&render_grid_comparison(
                    "Step 0",
                    prev_grid,
                    m.cg_iter_mean,
                    m.wallclock_total.as_secs_f64(),
                    m.mass_drift_relative,
                    m.max_abs_mean_vx,
                    m.max_abs_mean_vy,
                    justif,
                ));
            }
        }

        out.push_str("### Dormant metrics (inactive at Step 1)\n\n");
        out.push_str("| metric | activated at |\n|---|---|\n");
        out.push_str("| S̃_eq (active-orogen mean thickness) | Step 5+ |\n");
        out.push_str("| boundary type diversity | Step 5 |\n");
        out.push_str("| yielding cell fraction | Step 3 |\n");
        out.push_str("| cratonic stability | Step 9 |\n");
        out.push_str("| age field stats | Step 10 |\n\n");
    }

    out.push_str("---\n");
    out.push_str("*Generated by `cargo run --release --bin step_baseline`.*\n");
    out
}
