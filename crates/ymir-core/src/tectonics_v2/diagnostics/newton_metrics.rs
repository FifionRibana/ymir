//! Step-1 Newton-specific aggregates.
//!
//! A single baseline run aggregates, per time step, the outcome of
//! the nonlinear solve and the number of CG iterations the inner
//! linear solve consumed. The totals are what ends up in the
//! markdown report.

#[derive(Clone, Debug, Default)]
pub struct NewtonAggregate {
    /// Outcome counts. Ordered so the report can render a stable row
    /// set even when one outcome doesn't appear.
    pub converged: usize,
    pub stalled: usize,
    pub diverged: usize,
    pub capped: usize,

    /// Newton outer iterations per timestep.
    pub outer_iters: Vec<u32>,
    /// CG iterations inside Newton outer iteration, aggregated over
    /// the run (one sample per Newton step).
    pub cg_iters_per_newton_step: Vec<usize>,

    /// Cap-activation counters — fraction of cells where
    /// `η_eff > 0.9·η_max`. Split into "ramp" (during the startup
    /// continuation) and "steady" (after), as the spec requires.
    pub cap_fraction_ramp_max: f64,
    pub cap_fraction_steady_max: f64,

    /// Effective η_max / η_min ratio, sampled once per timestep and
    /// reduced to mean/max across the run.
    pub eta_contrast_samples: Vec<f64>,

    /// Startup continuation outcome: `Some(true)` if every ramp
    /// sub-solve converged; `Some(false)` if any sub-solve failed;
    /// `None` if continuation wasn't run (no step in this run).
    pub continuation_all_converged: Option<bool>,
    pub continuation_iters_used: u32,

    /// Step 3 — plastic yielding metrics. `None` when yielding is
    /// `Disabled`. When `Enabled`, `bi_diagnostic` carries the run's
    /// Bi value and the two cell-fraction / intensity aggregates are
    /// max-over-timestep values (peaks picked over the run).
    pub bi_diagnostic: Option<f64>,
    /// Fraction of cells where `η_eff < 0.5 · η_visc`. Max over the run.
    pub yielding_cell_fraction_max: Option<f64>,
    /// Mean of `(η_visc / η_eff − 1)` over cells where
    /// `η_eff < 0.9 · η_visc`. Max over the run; zero if no cell
    /// meets the threshold during the run.
    pub yielding_intensity_max: Option<f64>,
    /// Domain-level `ε̇_II` aggregates at the **final timestep**,
    /// intended for the floor-dominated-regime diagnostic in the
    /// Step 3 physics report. `None` when yielding is Disabled.
    pub eps_ii_mean_final: Option<f64>,
    pub eps_ii_max_final: Option<f64>,
    /// Fraction of cells where `ε̇_II < 10·ε̇_min` at the final
    /// timestep — "how much of the domain sits in the floor-
    /// dominated band at the end of the run".
    pub eps_ii_floor_dominated_fraction_final: Option<f64>,
}

impl NewtonAggregate {
    /// Summary helper for the markdown writer.
    pub fn outer_iters_mean(&self) -> f64 {
        if self.outer_iters.is_empty() {
            0.0
        } else {
            self.outer_iters.iter().map(|&v| v as f64).sum::<f64>() / self.outer_iters.len() as f64
        }
    }
    pub fn outer_iters_max(&self) -> u32 {
        self.outer_iters.iter().copied().max().unwrap_or(0)
    }
    pub fn cg_iters_per_newton_mean(&self) -> f64 {
        if self.cg_iters_per_newton_step.is_empty() {
            0.0
        } else {
            self.cg_iters_per_newton_step.iter().sum::<usize>() as f64
                / self.cg_iters_per_newton_step.len() as f64
        }
    }
    pub fn cg_iters_per_newton_max(&self) -> usize {
        self.cg_iters_per_newton_step.iter().copied().max().unwrap_or(0)
    }
    pub fn eta_contrast_mean(&self) -> f64 {
        if self.eta_contrast_samples.is_empty() {
            1.0
        } else {
            self.eta_contrast_samples.iter().sum::<f64>()
                / self.eta_contrast_samples.len() as f64
        }
    }
    pub fn eta_contrast_max(&self) -> f64 {
        self.eta_contrast_samples
            .iter()
            .copied()
            .fold(1.0_f64, f64::max)
    }

    /// Percentage outcome distribution; rounds to one decimal place.
    pub fn outcome_percentages(&self) -> (f64, f64, f64, f64) {
        let total = (self.converged + self.stalled + self.diverged + self.capped) as f64;
        if total == 0.0 {
            return (0.0, 0.0, 0.0, 0.0);
        }
        (
            100.0 * self.converged as f64 / total,
            100.0 * self.stalled as f64 / total,
            100.0 * self.diverged as f64 / total,
            100.0 * self.capped as f64 / total,
        )
    }
}

/// Compute the fraction of cells in `eta_cc` where the soft cap is
/// close to active — `η_eff > 0.9 · η_max_cap`.
pub fn cap_activation_fraction(eta_cc: &crate::tectonics_v2::field::Field2D, eta_max_cap: f64) -> f64 {
    let n = eta_cc.data().len();
    if n == 0 {
        return 0.0;
    }
    let threshold = 0.9 * eta_max_cap;
    let count = eta_cc.data().iter().filter(|&&v| v > threshold).count();
    count as f64 / n as f64
}

/// Fraction of cells where `η_eff < 0.5 · η_visc`.
///
/// Step 3 primary yielding metric (issue #85). This definition
/// captures "yielding is the dominant branch" rather than
/// "yielding is present anywhere" — the legacy `η_p < η_v`
/// criterion saturated to ~1.0 as soon as the plastic branch was
/// defined, carrying no diagnostic signal. See the `#75`
/// discussion.
pub fn yielding_cell_fraction(
    eta_visc: &crate::tectonics_v2::field::Field2D,
    eta_eff: &crate::tectonics_v2::field::Field2D,
) -> f64 {
    let n = eta_visc.data().len();
    if n == 0 {
        return 0.0;
    }
    let mut count = 0usize;
    for (&ev, &ee) in eta_visc.data().iter().zip(eta_eff.data().iter()) {
        if ee < 0.5 * ev {
            count += 1;
        }
    }
    count as f64 / n as f64
}

/// Mean of `(η_visc / η_eff − 1)` over cells where
/// `η_eff < 0.9 · η_visc`. Zero if no cell meets the threshold.
/// Captures "how much the yielding softens the already-yielding
/// zones", orthogonal to `yielding_cell_fraction` which captures
/// "how widespread the yielding is".
pub fn yielding_intensity(
    eta_visc: &crate::tectonics_v2::field::Field2D,
    eta_eff: &crate::tectonics_v2::field::Field2D,
) -> f64 {
    let mut sum = 0.0_f64;
    let mut count = 0usize;
    for (&ev, &ee) in eta_visc.data().iter().zip(eta_eff.data().iter()) {
        if ee < 0.9 * ev && ee > 0.0 {
            sum += ev / ee - 1.0;
            count += 1;
        }
    }
    if count == 0 { 0.0 } else { sum / count as f64 }
}

/// Compute `η_max / η_min` from an η field.
pub fn eta_contrast(eta_cc: &crate::tectonics_v2::field::Field2D) -> f64 {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &v in eta_cc.data() {
        if v > 0.0 {
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
        }
    }
    if min.is_finite() && min > 0.0 {
        max / min
    } else {
        1.0
    }
}
