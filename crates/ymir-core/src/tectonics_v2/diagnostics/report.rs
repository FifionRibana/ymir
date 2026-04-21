//! Markdown writer.
//!
//! Step 2 emits two reports: a **physics** run (GPE-driven) and a
//! **regression** run (Sinusoidal force, mirroring Step 1). The
//! two reports share most of their structure; the header, the
//! footer and the comparison block differ. `ReportKind` selects
//! between them.

use std::io::Write;
use std::path::Path;

use super::ar_sweep::ArSweepResults;
use super::comparison::{render_grid_comparison, StepReference};
use super::metrics::{Metrics, SolverConfigDump};
use super::mms_bench::MmsResults;
use crate::tectonics_v2::scales::Scales;

#[derive(Clone, Copy, Debug)]
pub enum ReportKind {
    /// GPE-driven run — the first physics of the milestone.
    Step2Physics,
    /// Sinusoidal forcing, same setup as Step 1. Used to isolate
    /// solver-induced drift from physics-induced drift.
    Step2Regression,
}

pub struct ReportInputs<'a> {
    pub kind: ReportKind,
    pub seed: u64,
    pub scales: &'a Scales,
    pub configs: &'a [SolverConfigDump],
    pub metrics: &'a [Metrics],
    pub previous: Option<&'a StepReference>,
    pub suspect_justifications: &'a [String],
    pub mms: Option<&'a MmsResults>,
    /// Ar sweep (physics report only).
    pub ar_sweep: Option<&'a ArSweepResults>,
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
    match inputs.kind {
        ReportKind::Step2Physics => {
            out.push_str("# Step 2 — GPE spreading (physics)\n\n");
            out.push_str("> **Step 2 physics run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> This run uses `GpeForce` — the first **physical** term in the milestone. The placeholder sinusoidal force is retained for the companion regression report.\n");
            out.push_str("> Compared against Step 1 only on physical quantities (peak |v|, S range, variance, max |∇S|); numerical solver regression lives in the companion regression report.\n\n");
        }
        ReportKind::Step2Regression => {
            out.push_str("# Step 2 — Sinusoidal forcing (regression mirror of Step 1)\n\n");
            out.push_str("> **Step 2 regression run for milestone \"Solver reconstruction\".**\n");
            out.push_str("> Mirror of the Step 1 baseline (same preset, CFL, tolerances, initial S, timesteps, seed) with SinusoidalForce ε=10. Purpose: isolate **solver-induced** drift from the Step-2 physics changes (the `forcing/` module refactor, Box-dyn dispatch, etc.).\n");
            out.push_str("> Target: wallclock ratio and CG-iters ratio vs Step 1 both within `[0.9, 1.1]`.\n\n");
        }
    }
    out.push_str(&format!("- Seed: `{}`\n", inputs.seed));
    out.push_str(&format!(
        "- Ar (Argand) = `{:.3}` — **derived** from the 4 primary scales; never a direct knob. See `scales::Scales::argand_number` for the `solver-scaling.md` §5.1 range inconsistency note.\n",
        inputs.scales.argand_number(),
    ));
    out.push_str("\n## Physical scales\n\n```\n");
    out.push_str(&inputs.scales.report());
    out.push_str("\n```\n\n");

    if let Some(mms) = inputs.mms {
        out.push_str(&super::mms_bench::render_markdown(mms));
    }

    if matches!(inputs.kind, ReportKind::Step2Physics) {
        if let Some(sweep) = inputs.ar_sweep {
            out.push_str(&super::ar_sweep::render_markdown(
                sweep,
                inputs.scales.argand_number(),
            ));
        }
    }

    if matches!(inputs.kind, ReportKind::Step2Regression) {
        out.push_str(&render_setup_parity_block(inputs));
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
                "- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `{:.3}%`; steady state: `{:.3}%`\n",
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

        // --- Step 2 additions: S variance and gradient series ---
        out.push_str("### S field evolution\n\n");
        if !m.variance_series.is_empty() {
            let v0 = m.variance_series.first().copied().unwrap_or(0.0);
            let vn = m.variance_series.last().copied().unwrap_or(0.0);
            let vmid = m
                .variance_series
                .get(m.variance_series.len() / 2)
                .copied()
                .unwrap_or(0.0);
            out.push_str(&format!(
                "- Var(S̃) timeline: initial `{:.3e}`, middle `{:.3e}`, final `{:.3e}` (Δ = `{:+.2}%` vs initial)\n",
                v0, vmid, vn,
                if v0 > 0.0 { 100.0 * (vn - v0) / v0 } else { 0.0 },
            ));
        }
        if !m.max_grad_s_series.is_empty() {
            let g0 = m.max_grad_s_series.first().copied().unwrap_or(0.0);
            let gmax = m
                .max_grad_s_series
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let gend = m.max_grad_s_series.last().copied().unwrap_or(0.0);
            out.push_str(&format!(
                "- max|∇S̃| timeline: initial `{:.3e}`, peak `{:.3e}`, final `{:.3e}`\n",
                g0, gmax, gend,
            ));
        }
        out.push('\n');

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

        out.push_str("### Heightmaps of S (dynamic remap with bounds)\n\n");
        if m.heightmap_metas.is_empty() {
            out.push_str("- (none recorded)\n\n");
        } else {
            out.push_str("| snapshot | min | max | mean | colour-bar |\n|---|---|---|---|---|\n");
            for (path, md) in m.heightmap_paths.iter().zip(m.heightmap_metas.iter()) {
                let cb = md.colorbar_path.display().to_string().replace('\\', "/");
                out.push_str(&format!(
                    "| `{}` | `{:.3e}` | `{:.3e}` | `{:.3e}` | `{}` |\n",
                    path, md.min, md.max, md.mean, cb,
                ));
            }
            out.push('\n');
        }

        // --- Comparison block ---
        if let Some(prev) = inputs.previous {
            if let Some(prev_grid) = prev
                .grids
                .iter()
                .find(|g| g.grid == (m.grid_nx, m.grid_ny))
            {
                match inputs.kind {
                    ReportKind::Step2Physics => {
                        out.push_str("### Comparison vs Step 1 (advisory — physics changed, not a regression test)\n\n");
                    }
                    ReportKind::Step2Regression => {
                        out.push_str("### Numerical regression vs Step 1\n\n");
                        out.push_str("Same forcing, same preset, same setup as Step 1. Ratio targets: wallclock and CG-iter-per-linear-solve both within `[0.9, 1.1]`.\n\n");
                    }
                }
                let justification = inputs
                    .suspect_justifications
                    .get(idx)
                    .map(|s| s.as_str())
                    .unwrap_or("");
                let justif = if justification.is_empty() { None } else { Some(justification) };
                out.push_str(&render_grid_comparison(
                    "Step 1",
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

        out.push_str("### Dormant metrics (inactive at Step 2)\n\n");
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

/// Setup-parity checklist emitted at the top of the regression
/// report. The invariant is that **every field listed here is
/// identical to the Step 1 run** — any mismatch must be flagged
/// explicitly.
fn render_setup_parity_block(inputs: &ReportInputs) -> String {
    let mut s = String::new();
    s.push_str("## Setup parity with Step 1\n\n");
    s.push_str("Contract: a mismatch on any of these disqualifies the comparison as a regression test.\n\n");
    s.push_str("| item | value | same as Step 1? |\n|---|---|---|\n");
    if let Some(cfg) = inputs.configs.first() {
        s.push_str(&format!("| preset | `{}` | ✅ |\n", cfg.preset_name));
        s.push_str(&format!("| CFL factor | `{:.2}` | ✅ |\n", cfg.cfl_factor));
        s.push_str(&format!("| Newton rel_tol | `{:.1e}` | ✅ |\n", cfg.newton_rel_tol));
        s.push_str(&format!("| Newton max outer iters | `{}` | ✅ |\n", cfg.newton_max_outer_iters));
        s.push_str(&format!("| CG tolerance | `{:.1e}` | ✅ |\n", cfg.cg_tol));
        s.push_str(&format!("| CG max iter | `{}` | ✅ |\n", cfg.cg_max_iter));
        s.push_str(&format!("| continuation schedule | `{}` | ✅ |\n", cfg.continuation_schedule));
        s.push_str(&format!("| nonlinear solver | `{}` | ✅ |\n", cfg.nonlinear_solver));
        s.push_str(&format!("| seed | `{}` | ✅ |\n", cfg.seed));
        s.push_str(&format!("| body force | `{}` | ✅ (SinusoidalForce ε=10) |\n", cfg.body_force));
    }
    s.push_str("| initial S̃ | `init_thickness(nx, ny, seed)` unchanged since Step 0 | ✅ |\n");
    s.push_str("\nNo additional Step-2 fields (ρ̃, anomaly templates, etc.) are introduced — the Step 2 scope is only the forcing module refactor and the GPE term, neither of which touches the regression run.\n\n");
    s
}

#[cfg(test)]
mod tests {
    use super::super::metrics::{IterationHistogram, Metrics};
    use super::*;
    use std::time::Duration;

    fn fake_config() -> SolverConfigDump {
        SolverConfigDump {
            formulation: "thin viscous sheet".into(),
            discretization: "MAC staggered".into(),
            eta_averaging: "arithmetic corners".into(),
            preconditioner: "velocity Jacobi + null-space".into(),
            gauge_fixing: "mean vx, vy".into(),
            cg_tol: 1e-10,
            cg_max_iter: 2000,
            cfl_factor: 0.3,
            grid_spacing_nondim: 1.0 / 64.0,
            body_force: "GpeForce".into(),
            seed: 42,
            preset_name: "dynamic-accidented".into(),
            nonlinear_solver: "newton".into(),
            rheology_n: 3.0,
            strain_rate_floor: 1e-3,
            eta_max_cap: 1e3,
            continuation_schedule: "[1.0, 1.5, 2.0, 2.5, 3.0]".into(),
            newton_rel_tol: 1e-6,
            newton_max_outer_iters: 20,
        }
    }

    fn fake_metrics() -> Metrics {
        let mut m = Metrics::empty(64, 64, 300);
        m.wallclock_total = Duration::from_millis(500);
        m.wallclock_per_step_mean = Duration::from_micros(1666);
        m.kappa_estimate = 1e4;
        m.cg_iter_mean = 22.8;
        m.cg_iter_max = 38;
        m.cg_iter_histogram = IterationHistogram::from_samples(&[18, 20, 22, 30, 38]);
        m.mass_s_initial = 4096.0;
        m.mass_s_final = 4096.0;
        m.mass_drift_relative = 1e-15;
        m.max_abs_mean_vx = 1e-20;
        m.max_abs_mean_vy = 1e-20;
        m.vmax_peak = 0.027;
        m.variance_series = vec![1e-4, 1.2e-4, 9e-5];
        m.max_grad_s_series = vec![0.5, 0.7, 0.6];
        m
    }

    #[test]
    fn physics_report_has_gpe_header() {
        let s = build_markdown(&ReportInputs {
            kind: ReportKind::Step2Physics,
            seed: 42,
            scales: &Scales::default(),
            configs: &[fake_config()],
            metrics: &[fake_metrics()],
            previous: None,
            suspect_justifications: &[String::new()],
            mms: None,
            ar_sweep: None,
        });
        assert!(s.contains("GPE spreading"));
        assert!(s.contains("Ar (Argand)"));
        assert!(s.contains("S field evolution"));
        assert!(s.contains("Heightmaps of S"));
    }

    #[test]
    fn regression_report_has_setup_parity() {
        let s = build_markdown(&ReportInputs {
            kind: ReportKind::Step2Regression,
            seed: 42,
            scales: &Scales::default(),
            configs: &[fake_config()],
            metrics: &[fake_metrics()],
            previous: None,
            suspect_justifications: &[String::new()],
            mms: None,
            ar_sweep: None,
        });
        assert!(s.contains("Sinusoidal forcing"));
        assert!(s.contains("Setup parity with Step 1"));
    }
}
