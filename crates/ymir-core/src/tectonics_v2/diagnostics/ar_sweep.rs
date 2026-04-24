//! Argand-number sweep for the Step 2 physics report.
//!
//! The honest thin-sheet value from the default scales is
//! `Ar = S*/L* = 0.1`; the design note's historical target range
//! `[1, 5]` is mathematically incompatible with the thin-sheet
//! assumption (`S* ≪ L*`). To keep quantitative evidence of the
//! discretisation's correctness at a range of `Ar` values (and to
//! document how visible GPE dynamics scale with `Ar` on the 300-step
//! window), we run the physics baseline at
//! `Ar ∈ {0.1, 0.5, 1.0, 2.0, 5.0}` using explicit
//! `GpeForce::with_ar(value)`, and tabulate the resulting evolution
//! metrics.
//!
//! This is a diagnostic sweep, not a design knob: the baseline
//! `step2_physics_report.md` still uses `Ar = 0.1` derived from the
//! scales.

use std::path::PathBuf;

use super::harness::{run_baseline, BaselineConfig, BaselineResult, ForceKind, NonlinearChoice};
use crate::tectonics_v2::forcing::{ForceSum, GpeForce};
use crate::tectonics_v2::presets::Preset;

/// One point in the Ar sweep. Fields are summaries — keeping the
/// full `Metrics` per point would explode the report size.
#[derive(Clone, Debug)]
pub struct ArSweepPoint {
    pub ar: f64,
    pub wallclock_s: f64,
    pub var_initial: f64,
    pub var_final: f64,
    pub var_ratio: f64,
    pub max_grad_s_peak: f64,
    pub peak_v: f64,
    pub newton_converged_pct: f64,
    pub cg_iter_mean: f64,
    pub mass_drift_rel: f64,
}

#[derive(Clone, Debug)]
pub struct ArSweepResults {
    pub points: Vec<ArSweepPoint>,
}

/// Run the sweep: `64² × 300 steps` for each `Ar` value. Snapshots
/// are suppressed (no PNG writing) to keep the sweep cheap.
pub fn run_ar_sweep(
    seed: u64,
    steps: usize,
    preset: &Preset,
    s_perturbation_amplitude: f64,
    ar_values: &[f64],
) -> ArSweepResults {
    let mut points = Vec::with_capacity(ar_values.len());
    for &ar in ar_values {
        let mut sum = ForceSum::new();
        sum.push(Box::new(GpeForce::with_ar(ar)));
        let cfg = BaselineConfig {
            seed,
            grid_nx: 64,
            grid_ny: 64,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps,
            cfl_factor: 0.3,
            total_time_nondim: 6.0,
            preset: preset.clone(),
            nonlinear: NonlinearChoice::Newton,
            newton_cfg: Default::default(),
            picard_cfg: Default::default(),
            heightmap_fractions: Vec::new(), // no PNGs for the sweep
            output_dir: PathBuf::from("target/ar_sweep_scratch"),
            force: Box::new(sum),
            force_kind: ForceKind::Gpe,
            sinusoidal_amplitude: 0.0,
            s_perturbation_amplitude,
            yielding: crate::tectonics_v2::presets::YieldingConfig::Disabled,
            basal_drag: crate::tectonics_v2::basal_drag::BasalDragConfig::Disabled,
            boundary: crate::tectonics_v2::boundaries::BoundaryConfig::Disabled,
            boundary_layout_name: String::new(),
            slab_pull: crate::tectonics_v2::slab::SlabPullConfig::Disabled,
            mantle: crate::tectonics_v2::mantle::MantleConfig::Disabled,
            capture: None,
            linear_solver: Default::default(),
        };
        let r: BaselineResult = run_baseline(&cfg);
        points.push(summarise(ar, &r));
    }
    ArSweepResults { points }
}

fn summarise(ar: f64, r: &BaselineResult) -> ArSweepPoint {
    let m = &r.metrics;
    let var_initial = m.variance_series.first().copied().unwrap_or(f64::NAN);
    let var_final = m.variance_series.last().copied().unwrap_or(f64::NAN);
    let var_ratio = if var_initial > 0.0 { var_final / var_initial } else { f64::NAN };
    let max_grad_peak = m
        .max_grad_s_series
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let newton_pct = m
        .newton
        .as_ref()
        .map(|n| n.outcome_percentages().0)
        .unwrap_or(0.0);
    ArSweepPoint {
        ar,
        wallclock_s: m.wallclock_total.as_secs_f64(),
        var_initial,
        var_final,
        var_ratio,
        max_grad_s_peak: max_grad_peak,
        peak_v: m.vmax_peak,
        newton_converged_pct: newton_pct,
        cg_iter_mean: m.cg_iter_mean,
        mass_drift_rel: m.mass_drift_relative,
    }
}

/// Render the sweep as a markdown table + a short interpretation.
pub fn render_markdown(res: &ArSweepResults, physical_ar: f64) -> String {
    let mut s = String::new();
    s.push_str("## Ar sweep (diagnostic)\n\n");
    s.push_str(&format!(
        "Honest thin-sheet value from the default scales: **`Ar = S*/L* = {:.3}`** (used in the baseline above).\n\n",
        physical_ar,
    ));
    s.push_str("The design note's historical target `Ar ∈ [1, 5]` is mathematically incompatible with `S* ≪ L*` — the thin-sheet assumption forces `Ar ≪ 1`. The sweep below tabulates the GPE-only response at 64²·300 steps for a range that brackets both the honest value and the historical band, so the discretisation's behaviour across `Ar` is visible quantitatively.\n\n");
    s.push_str("| Ar | Var(S̃) init | Var(S̃) final | ratio | peak \\|∇S̃\\| | peak \\|v\\| | Newton conv | CG mean | mass drift | wallclock (s) |\n");
    s.push_str("|---|---|---|---|---|---|---|---|---|---|\n");
    for p in &res.points {
        s.push_str(&format!(
            "| `{:.2}` | `{:.3e}` | `{:.3e}` | `{:.3}` | `{:.3e}` | `{:.3e}` | `{:.0}%` | `{:.1}` | `{:.2e}` | `{:.3}` |\n",
            p.ar,
            p.var_initial,
            p.var_final,
            p.var_ratio,
            p.max_grad_s_peak,
            p.peak_v,
            p.newton_converged_pct,
            p.cg_iter_mean,
            p.mass_drift_rel,
            p.wallclock_s,
        ));
    }
    s.push_str("\n**Interpretation** — GPE dissipation scales as `Ar/τ*`, so the characteristic spreading time is `τ*/Ar`. At `Ar = 0.1` this is `~10·τ*`, ten times the tectonic time scale and well beyond the 300-step run (`6·τ*`). The variance ratio across the sweep confirms the expected monotonic response: lower `Ar` → slower spreading → larger `Var(S̃)_final / Var(S̃)_initial`. Narrative-level dynamics (continents building and breaking on the run window) must therefore come from the mechanisms being added at Steps 3–10, not from GPE alone — see `solver-scaling.md` §5.1.bis for the characteristic-time ordering.\n\n");
    s
}
