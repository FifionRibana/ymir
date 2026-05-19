//! `k_sub` sweep for the Step 5 physics report.
//!
//! Runs the baseline (64²·N steps, `GpeForce`, yielding Enabled,
//! basal drag Enabled, boundary Enabled with `horizontal_oceanic_strip`
//! layout) across `k_sub ∈ {0.3, 0.5, 0.7, 1.0}` and tabulates per-
//! point aggregates. The physics prediction is that `s_oceanic_mean`
//! is **strictly decreasing** with `k_sub` (more subduction → more
//! oceanic mass consumed per unit convergent motion); the sweep's
//! primary acceptance invariant is that strict monotonicity.
//!
//! All other parameters are held at baseline; `k_spread` is expected
//! to have been pre-calibrated for the baseline `k_sub = 0.5` and
//! the same value is reused across the sweep. As `k_sub` moves away
//! from 0.5, `s_oceanic_mean` drifts out of `[0.18, 0.22]` — that is
//! the expected response, not a failure.

use std::path::PathBuf;

use super::harness::{run_baseline, BaselineConfig, BaselineResult, ForceKind, NonlinearChoice};
use crate::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use crate::tectonics_v2::boundaries::{horizontal_oceanic_strip, BoundaryRates};
use crate::tectonics_v2::forcing::{ForceSum, GpeForce};
use crate::tectonics_v2::presets::{Preset, YieldingConfig};
use crate::tectonics_v2::rheology::YieldingLaw;

#[derive(Clone, Debug)]
pub struct KSubSweepPoint {
    pub k_sub: f64,
    pub wallclock_s: f64,
    pub newton_converged_pct: f64,
    pub cg_iter_mean: f64,
    pub newton_iter_mean: f64,
    pub peak_v: f64,
    pub s_oceanic_mean: Option<f64>,
    pub s_continental_interior_mean: Option<f64>,
    pub s_continental_collision_mean: Option<f64>,
    pub clamp_activation_fraction_mean: Option<f64>,
    pub mass_balance_residual: Option<f64>,
    pub yielding_cell_fraction: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct KSubSweepResults {
    pub points: Vec<KSubSweepPoint>,
    pub peak_v_mono_ok: bool,
    pub s_oceanic_mono_ok: bool,
    pub mass_residual_ok: bool,
}

pub fn run_k_sub_sweep(
    seed: u64,
    steps: usize,
    preset: &Preset,
    s_perturbation_amplitude: f64,
    k_sub_values: &[f64],
    k_spread_calibrated: f64,
    bi: f64,
    br: f64,
) -> KSubSweepResults {
    let nx = 64;
    let ny = 64;
    let mut points = Vec::with_capacity(k_sub_values.len());
    for &k_sub in k_sub_values {
        let mut sum = ForceSum::new();
        sum.push(Box::new(GpeForce::with_ar(0.1)));
        let layout = horizontal_oceanic_strip(nx, ny);
        let layout_name = layout.name;
        let rates = BoundaryRates::baseline_uncalibrated()
            .with_k_spread(k_spread_calibrated)
            .with_k_sub(k_sub);
        let boundary = layout.into_config(rates);
        let cfg = BaselineConfig {
            seed,
            grid_nx: nx,
            grid_ny: ny,
            domain_lx: 1.0,
            domain_ly: 1.0,
            steps,
            cfl_factor: 0.3,
            total_time_nondim: 6.0,
            preset: preset.clone(),
            nonlinear: NonlinearChoice::Newton,
            newton_cfg: Default::default(),
            picard_cfg: Default::default(),
            heightmap_fractions: Vec::new(),
            output_dir: PathBuf::from("target/k_sub_sweep_scratch"),
            force: Box::new(sum),
            force_kind: ForceKind::Gpe,
            sinusoidal_amplitude: 0.0,
            s_perturbation_amplitude,
            yielding: YieldingConfig::Enabled(YieldingLaw { bi, ..Default::default() }),
            basal_drag: BasalDragConfig::Enabled(BasalDragLaw { br, ..BasalDragLaw::default() }),
            boundary,
            boundary_layout_name: layout_name.to_string(),
            slab_pull: crate::tectonics_v2::slab::SlabPullConfig::Disabled,
            mantle: crate::tectonics_v2::mantle::MantleConfig::Disabled,
            cratonic: crate::tectonics_v2::cratonic::CratonicConfig::Disabled,
            age_field: crate::tectonics_v2::age_field::AgeFieldConfig::Disabled,
            capture: None,
            linear_solver: Default::default(),
            init_mode: crate::tectonics_v2::init::InitMode::Checkerboard,
            continuation: None,
            plate_kinematic: crate::tectonics_v2::plate_kinematic::PlateKinematicConfig::Zero,
        };
        let r: BaselineResult = run_baseline(&cfg);
        points.push(summarise(k_sub, &r));
    }
    // Strict monotonicity: s_oceanic_mean must decrease with k_sub.
    let mut s_ocean_mono = true;
    let mut peak_v_mono = true;
    for w in points.windows(2) {
        let a = w[0].s_oceanic_mean.unwrap_or(f64::NAN);
        let b = w[1].s_oceanic_mean.unwrap_or(f64::NAN);
        if !(a.is_finite() && b.is_finite() && b < a) {
            s_ocean_mono = false;
        }
        // peak|v| monotonicity is a weaker invariant (drag may
        // wobble); we check it but do not gate on it.
        if w[1].peak_v > w[0].peak_v * 1.10 {
            peak_v_mono = false;
        }
    }
    let mass_residual_ok = points
        .iter()
        .all(|p| p.mass_balance_residual.map(|r| r < 0.01).unwrap_or(false));
    KSubSweepResults {
        points,
        peak_v_mono_ok: peak_v_mono,
        s_oceanic_mono_ok: s_ocean_mono,
        mass_residual_ok,
    }
}

fn summarise(k_sub: f64, r: &BaselineResult) -> KSubSweepPoint {
    let m = &r.metrics;
    let newton = m.newton.as_ref();
    let (conv_pct, newton_iter_mean) = newton
        .map(|n| (n.outcome_percentages().0, n.outer_iters_mean()))
        .unwrap_or((0.0, 0.0));
    KSubSweepPoint {
        k_sub,
        wallclock_s: m.wallclock_total.as_secs_f64(),
        newton_converged_pct: conv_pct,
        cg_iter_mean: m.cg_iter_mean,
        newton_iter_mean,
        peak_v: m.vmax_peak,
        s_oceanic_mean: newton.and_then(|n| n.s_oceanic_mean),
        s_continental_interior_mean: newton.and_then(|n| n.s_continental_interior_mean),
        s_continental_collision_mean: newton.and_then(|n| n.s_continental_collision_mean),
        clamp_activation_fraction_mean: newton.and_then(|n| n.clamp_activation_fraction_mean),
        mass_balance_residual: newton.and_then(|n| n.mass_balance_residual),
        yielding_cell_fraction: newton.and_then(|n| n.yielding_cell_fraction_max),
    }
}

/// Render the sweep as a markdown table + monotonicity verdicts.
pub fn render_markdown(res: &KSubSweepResults) -> String {
    let mut s = String::new();
    s.push_str("## k_sub sweep (diagnostic)\n\n");
    s.push_str("Baseline `k_sub = 0.5` (preset `dynamic-accidented`, layout `horizontal_oceanic_strip`, `k_spread` pre-calibrated). The sweep covers `k_sub ∈ {0.3, 0.5, 0.7, 1.0}` at 64²·N steps with `GpeForce (Ar = 0.1)` + yielding Enabled + basal drag Enabled + boundary Enabled. Physical prediction: higher `k_sub` consumes more oceanic mass per unit convergent motion, so `s_oceanic_mean` strictly decreases with `k_sub`. That strict monotonicity is the acceptance invariant of the Step 5 sweep (issue #89).\n\n");
    s.push_str("| k_sub | s_oceanic_mean | s_cont_interior | s_cont_collision | peak \\|v\\| | CG iters | Newton iters | clamp frac mean | mass_balance_res | Newton conv | wallclock (s) |\n");
    s.push_str("|---|---|---|---|---|---|---|---|---|---|---|\n");
    for p in &res.points {
        // 6 decimals on s_oceanic_mean — monotonicity at the Step 5
        // baseline with GPE-only forcing lives at the 5th–6th digit
        // (peak|v| ≈ 1e-5 → Q_sub ≈ k_sub·1e-5 per step, drain tiny
        // but strict). Fewer decimals would hide the physical signal.
        let soc = p.s_oceanic_mean.map(|v| format!("`{:.6}`", v)).unwrap_or_else(|| "`—`".into());
        let sci = p
            .s_continental_interior_mean
            .map(|v| format!("`{:.4}`", v))
            .unwrap_or_else(|| "`—`".into());
        let scc = p
            .s_continental_collision_mean
            .map(|v| format!("`{:.4}`", v))
            .unwrap_or_else(|| "`—`".into());
        let cfm = p
            .clamp_activation_fraction_mean
            .map(|v| format!("`{:.3e}`", v))
            .unwrap_or_else(|| "`—`".into());
        let mbr = p
            .mass_balance_residual
            .map(|v| format!("`{:.3e}`", v))
            .unwrap_or_else(|| "`—`".into());
        s.push_str(&format!(
            "| `{:.2}` | {} | {} | {} | `{:.3e}` | `{:.1}` | `{:.1}` | {} | {} | `{:.0}%` | `{:.2}` |\n",
            p.k_sub,
            soc,
            sci,
            scc,
            p.peak_v,
            p.cg_iter_mean,
            p.newton_iter_mean,
            cfm,
            mbr,
            p.newton_converged_pct,
            p.wallclock_s,
        ));
    }

    s.push_str(&format!(
        "\n**Monotonicity of `s_oceanic_mean` vs k_sub** (strictly decreasing): {}\n",
        if res.s_oceanic_mono_ok {
            "✅ strictly decreasing across the 4 points, as required."
        } else {
            "❌ violation detected — `s_oceanic_mean` did not strictly decrease. Remonter (issue #89 acceptance): investigate the source/sink assembly or the layout's divergence profile before treating the Step 5 physics as validated."
        },
    ));
    s.push_str(&format!(
        "**mass_balance_residual < 1% across all points**: {}\n",
        if res.mass_residual_ok {
            "✅ residual bounded at every point; the flux accounting (Q + clamp) holds uniformly."
        } else {
            "❌ violation — at least one sweep point has `mass_balance_residual ≥ 1%`. The flux accounting under-captures either the physical `Q` integral or the artificial clamp contribution; re-check [`apply_clamp_with_tracking`] and the `q_integral` accumulator before merging."
        },
    ));
    s.push_str("\n**Interpretation** — lower `k_sub` leaves more oceanic mass in place; higher `k_sub` drains the strip, pushing `s_oceanic_mean` toward `S_MIN = 0.05`. The `s_continental_collision_mean` column is tracked telemetry (issue #89 acceptance notes): collision-row thickening can drift outside any reference band and that is the expected physics, not a failure.\n\n");
    s
}
