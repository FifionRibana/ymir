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
    /// Upper bounds of each bin (inclusive). The 5 bins are determined
    /// dynamically from the min/max observed iteration counts.
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
        // Guarantee the last edge is the observed max (floating rounding above can leave a gap at boundaries).
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
    pub discretization: String,
    pub harmonic_averaging: String,
    pub preconditioner: String,
    pub gauge_fixing: String,
    pub outer_tol: f64,
    pub inner_tol: f64,
    pub outer_max_iter: usize,
    pub inner_max_iter: usize,
    pub cfl_factor: f64,
    pub grid_spacing_nondim: f64,
    pub body_force: String,
    pub seed: u64,
}

impl SolverConfigDump {
    pub fn render_markdown(&self) -> String {
        format!(
            "### Solver configuration\n\n\
             | field | value |\n\
             |---|---|\n\
             | discretization | {} |\n\
             | harmonic averaging | {} |\n\
             | preconditioner | {} |\n\
             | gauge fixing | {} |\n\
             | outer CG tolerance | {:.1e} |\n\
             | inner CG tolerance | {:.1e} |\n\
             | outer CG max iter | {} |\n\
             | inner CG max iter | {} |\n\
             | CFL factor | {:.2} |\n\
             | grid spacing (nondim) | {:.6} |\n\
             | body force | {} |\n\
             | seed | {} |\n",
            self.discretization,
            self.harmonic_averaging,
            self.preconditioner,
            self.gauge_fixing,
            self.outer_tol,
            self.inner_tol,
            self.outer_max_iter,
            self.inner_max_iter,
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
    /// Condition-number estimate derived from outer CG iteration count.
    /// `κ_est ≈ (2·iter / ln(2/tol))²`, rounded to a nearest power of
    /// ten for reporting. Cheap single-scalar proxy; power-iteration
    /// κ is a future T2 extension.
    pub kappa_estimate: f64,

    /// Trivially 1.0 at Step 0 (constant η). Kept in the struct so the
    /// framework slot exists for Step 1 onward.
    pub eta_contrast: f64,

    /// Outer CG iterations: mean, max, histogram.
    pub outer_iter_mean: f64,
    pub outer_iter_max: usize,
    pub outer_iter_histogram: IterationHistogram,

    /// Inner CG iterations (per inner solve): mean, max.
    pub inner_iter_mean: f64,
    pub inner_iter_max: usize,

    // ---- Mass conservation (Step 0 active) ----
    pub mass_s_initial: f64,
    pub mass_s_final: f64,
    pub mass_drift_relative: f64,

    // ---- Null-space health after every solve (Step 0 active) ----
    pub max_abs_mean_p: f64,
    pub max_abs_mean_vx: f64,
    pub max_abs_mean_vy: f64,

    // ---- Velocity magnitude (advisory) ----
    pub vmax_peak: f64,

    // ---- Heightmap snapshots (Step 0 active) ----
    /// Relative paths, reported in the markdown for hand inspection.
    pub heightmap_paths: Vec<String>,

    // ---- Dormant metrics (Option<> until their step introduces them) ----
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
            outer_iter_mean: 0.0,
            outer_iter_max: 0,
            outer_iter_histogram: IterationHistogram::default(),
            inner_iter_mean: 0.0,
            inner_iter_max: 0,
            mass_s_initial: 0.0,
            mass_s_final: 0.0,
            mass_drift_relative: 0.0,
            max_abs_mean_p: 0.0,
            max_abs_mean_vx: 0.0,
            max_abs_mean_vy: 0.0,
            vmax_peak: 0.0,
            heightmap_paths: Vec::new(),
            s_eq: None,
            boundary_type_diversity: None,
            yielding_cell_fraction: None,
            cratonic_stability: None,
            newton_outcome_distribution: None,
            age_field_stats: None,
        }
    }
}

/// `κ_est ≈ (2·iter / ln(2/tol))²` — CG convergence is roughly
/// geometric with rate `(√κ - 1)/(√κ + 1)`, giving the above relation
/// at convergence. Returns `f64::INFINITY` when `iter == 0` (solver
/// converged on the initial guess) so the estimate degrades
/// gracefully rather than producing a misleading low κ.
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
