//! Metric structure for a single baseline run.
//!
//! Active Step-0 metrics are bare fields; metrics that come online at
//! later steps are `Option<...>` and stay `None` in Step 0's report.
//! This distinguishes "framework slot for a future metric" from
//! "forgot to compute this".

use std::time::Duration;

/// Histogram of iteration counts bucketed into 5 bins.
#[derive(Clone, Debug, Default)]
pub struct IterationHistogram {
    /// Upper bounds of each bin (inclusive).
    pub bin_edges: [usize; 5],
    /// Number of solves landing in each bin.
    pub counts: [usize; 5],
}

impl IterationHistogram {
    pub fn from_samples(samples: &[usize]) -> Self {
        if samples.is_empty() {
            return Self::default();
        }
        let min = *samples.iter().min().unwrap();
        let max = *samples.iter().max().unwrap();
        if max == min {
            let edges = [min, min, min, min, min];
            let mut counts = [0usize; 5];
            counts[0] = samples.len();
            return Self { bin_edges: edges, counts };
        }
        let span = max - min;
        let mut edges = [0usize; 5];
        for (k, e) in edges.iter_mut().enumerate() {
            *e = min + ((k + 1) * span) / 5;
        }
        edges[4] = max;
        let mut counts = [0usize; 5];
        for &s in samples {
            let mut placed = false;
            for (b, e) in edges.iter().enumerate() {
                if s <= *e {
                    counts[b] += 1;
                    placed = true;
                    break;
                }
            }
            if !placed {
                counts[4] += 1;
            }
        }
        Self { bin_edges: edges, counts }
    }
}

/// Full solver-configuration record captured alongside any wallclock
/// number. **Mandatory** per lesson 3 of the session handoff:
/// wallclock without solver config is not a comparable metric.
#[derive(Clone, Debug)]
pub struct SolverConfigDump {
    pub formulation: String,
    pub discretization: String,
    pub eta_averaging: String,
    pub preconditioner: String,
    pub gauge_fixing: String,
    pub cg_tol: f64,
    pub cg_max_iter: usize,
    pub cfl_factor: f64,
    pub grid_spacing_nondim: f64,
    pub body_force: String,
    pub seed: u64,
    // --- Step 1 additions ---
    pub preset_name: String,
    pub nonlinear_solver: String,
    pub rheology_n: f64,
    pub strain_rate_floor: f64,
    pub eta_max_cap: f64,
    pub continuation_schedule: String,
    pub newton_rel_tol: f64,
    pub newton_max_outer_iters: u32,
}

impl SolverConfigDump {
    pub fn render_markdown(&self) -> String {
        format!(
            "### Solver configuration\n\n\
             | field | value |\n\
             |---|---|\n\
             | formulation | {} |\n\
             | discretization | {} |\n\
             | η averaging to corners | {} |\n\
             | preconditioner | {} |\n\
             | gauge fixing | {} |\n\
             | preset | `{}` |\n\
             | nonlinear solver | `{}` |\n\
             | rheology `n` (after continuation) | {:.2} |\n\
             | rheology `ε̇_min` | {:.1e} |\n\
             | rheology `η_max` (soft cap) | {:.1e} |\n\
             | continuation schedule | `{}` |\n\
             | Newton rel tol | {:.1e} |\n\
             | Newton max outer iters | {} |\n\
             | CG tolerance | {:.1e} |\n\
             | CG max iter | {} |\n\
             | CFL factor | {:.2} |\n\
             | grid spacing (nondim) | {:.6} |\n\
             | body force | {} |\n\
             | seed | {} |\n",
            self.formulation,
            self.discretization,
            self.eta_averaging,
            self.preconditioner,
            self.gauge_fixing,
            self.preset_name,
            self.nonlinear_solver,
            self.rheology_n,
            self.strain_rate_floor,
            self.eta_max_cap,
            self.continuation_schedule,
            self.newton_rel_tol,
            self.newton_max_outer_iters,
            self.cg_tol,
            self.cg_max_iter,
            self.cfl_factor,
            self.grid_spacing_nondim,
            self.body_force,
            self.seed,
        )
    }
}

/// Per-grid metrics for one baseline run.
#[derive(Clone, Debug)]
pub struct Metrics {
    pub grid_nx: usize,
    pub grid_ny: usize,
    pub steps: usize,

    // ---- Timing ----
    pub wallclock_total: Duration,
    pub wallclock_per_step_mean: Duration,

    // ---- Solver health (Step 0 active) ----
    /// Condition-number estimate from CG iteration count: with a
    /// tolerance `tol` and `k` CG iterations to converge from a zero
    /// initial guess, one has (Saad, *Iterative Methods*, §6.7)
    /// `k ≈ ½ √κ · ln(2/tol)`, so `κ ≈ (2 k / ln(2/tol))²`. For
    /// `tol = 1e-8`, `ln(2/tol) ≈ 19.1`, hence the simplified
    /// `κ ≈ (k / 9.6)² ~ 0.01 k²` commonly used as a reporting proxy.
    /// Reported directly from the mean CG iteration count over the run.
    pub kappa_estimate: f64,

    /// Trivially 1.0 at Step 0 (constant η). Kept as framework slot
    /// for Step 1 onward.
    pub eta_contrast: f64,

    /// CG iterations per sheet solve: mean, max, histogram.
    pub cg_iter_mean: f64,
    pub cg_iter_max: usize,
    pub cg_iter_histogram: IterationHistogram,

    // ---- Mass conservation (Step 0 active) ----
    pub mass_s_initial: f64,
    pub mass_s_final: f64,
    pub mass_drift_relative: f64,

    // ---- Null-space health after every solve (Step 0 active) ----
    pub max_abs_mean_vx: f64,
    pub max_abs_mean_vy: f64,

    // ---- Velocity magnitude (advisory) ----
    pub vmax_peak: f64,

    // ---- Heightmap snapshots (Step 0 active) ----
    pub heightmap_paths: Vec<String>,
    /// Step-2 addition: per-snapshot metadata (min, max, mean,
    /// colorbar path). Parallel to `heightmap_paths`.
    pub heightmap_metas: Vec<super::heightmap::HeightmapMetadata>,

    // ---- Physical series (Step 2) ----
    /// Variance of S at every macro step (`steps + 1` entries,
    /// including the initial state before the first solve).
    pub variance_series: Vec<f64>,
    /// Max of `|∇S|` at every macro step (same layout).
    pub max_grad_s_series: Vec<f64>,

    // ---- Step 1 active: Newton aggregate. None at Step 0 (no nonlinear solve). ----
    pub newton: Option<super::newton_metrics::NewtonAggregate>,

    // ---- Dormant metrics ----
    pub s_eq: Option<f64>,
    pub boundary_type_diversity: Option<BoundaryTypeCounts>,
    pub yielding_cell_fraction: Option<f64>,
    pub cratonic_stability: Option<f64>,
    pub newton_outcome_distribution: Option<NewtonOutcomeCounts>,
    pub age_field_stats: Option<AgeFieldStats>,
}

#[derive(Clone, Debug, Default)]
pub struct BoundaryTypeCounts {
    pub subduction: usize,
    pub collision: usize,
    pub rift: usize,
}

#[derive(Clone, Debug, Default)]
pub struct NewtonOutcomeCounts {
    pub converged: usize,
    pub stagnation: usize,
    pub oscillation: usize,
    pub max_iter: usize,
}

#[derive(Clone, Debug, Default)]
pub struct AgeFieldStats {
    pub min: f64,
    pub mean: f64,
    pub max: f64,
}

impl Metrics {
    pub fn empty(grid_nx: usize, grid_ny: usize, steps: usize) -> Self {
        Self {
            grid_nx,
            grid_ny,
            steps,
            wallclock_total: Duration::ZERO,
            wallclock_per_step_mean: Duration::ZERO,
            kappa_estimate: 0.0,
            eta_contrast: 1.0,
            cg_iter_mean: 0.0,
            cg_iter_max: 0,
            cg_iter_histogram: IterationHistogram::default(),
            mass_s_initial: 0.0,
            mass_s_final: 0.0,
            mass_drift_relative: 0.0,
            max_abs_mean_vx: 0.0,
            max_abs_mean_vy: 0.0,
            vmax_peak: 0.0,
            heightmap_paths: Vec::new(),
            heightmap_metas: Vec::new(),
            variance_series: Vec::new(),
            max_grad_s_series: Vec::new(),
            newton: None,
            s_eq: None,
            boundary_type_diversity: None,
            yielding_cell_fraction: None,
            cratonic_stability: None,
            newton_outcome_distribution: None,
            age_field_stats: None,
        }
    }
}

/// `κ ≈ (2·iter / ln(2/tol))²` — CG converges geometrically at rate
/// `(√κ - 1)/(√κ + 1)`; matching against the residual reduction
/// `tol/2` gives the formula above. Returns `NaN` for degenerate
/// inputs rather than a misleading small κ.
pub fn condition_number_estimate(iterations: usize, tol: f64) -> f64 {
    if iterations == 0 || tol <= 0.0 || tol >= 1.0 {
        return f64::NAN;
    }
    let num = 2.0 * iterations as f64;
    let den = (2.0 / tol).ln();
    (num / den).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_of_single_value() {
        let h = IterationHistogram::from_samples(&[5, 5, 5, 5]);
        assert_eq!(h.counts[0], 4);
        for b in 1..5 {
            assert_eq!(h.counts[b], 0);
        }
    }

    #[test]
    fn histogram_distributes_across_bins() {
        let samples: Vec<usize> = (0..25).collect();
        let h = IterationHistogram::from_samples(&samples);
        let total: usize = h.counts.iter().sum();
        assert_eq!(total, 25);
    }

    #[test]
    fn kappa_estimate_is_finite_for_reasonable_inputs() {
        let k = condition_number_estimate(30, 1e-8);
        assert!(k > 0.0 && k.is_finite(), "k = {}", k);
    }
}
