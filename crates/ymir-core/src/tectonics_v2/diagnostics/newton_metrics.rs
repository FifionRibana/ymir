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

    // ---- Step 4 — basal-drag diagnostics ----
    //
    // All three fields are `None` when `BasalDragConfig::Disabled`;
    // under `Enabled(law)`, `br_diagnostic` carries `law.br` and the
    // two ratios are means across the run (mean per step of per-cell
    // means). `peak_v_damping_ratio` is NOT computed here — it's a
    // cross-run quantity (physics vs regression) computed at report
    // rendering time.
    pub br_diagnostic: Option<f64>,
    /// `mean_cells(Br·S̃² / (Br·S̃² + η/Δx²))`, averaged across the
    /// run. Saturates at 1 when drag dominates; baseline Step-4 value
    /// is expected `≪ 1` (drag much smaller than viscous diagonal).
    pub basal_drag_energy_ratio: Option<f64>,
    /// `mean_cells(Br·S̃² / (η/Δx²))`, averaged across the run.
    /// Linear in Br; baseline Step-4 target band `[10⁻⁶, 10⁻⁴]` at
    /// 128² per the Step-4 spec algebra.
    pub drag_vs_visc_diagonal_ratio: Option<f64>,

    // ---- Step 5 — boundary source/sink diagnostics ----
    //
    // All fields are `None` under `BoundaryConfig::Disabled`; under
    // `Enabled` they are populated by the harness at the end of the
    // run. The mean/std fields are domain-wise stats on the final
    // `S̃`; the mass-balance residual and clamp activation fraction
    // are integrated over the whole run.
    pub boundary_layout_name: Option<&'static str>,
    pub s_oceanic_mean: Option<f64>,
    pub s_oceanic_std: Option<f64>,
    pub s_continental_interior_mean: Option<f64>,
    pub s_continental_interior_std: Option<f64>,
    pub s_continental_collision_mean: Option<f64>,
    /// Count of boundary-flag variants active (Q ≠ 0) during the run.
    /// Integer 0..=4 (None, plus Subduction|OceanicSubduction counted
    /// once, Rift, ContinentalCollision).
    pub boundary_type_diversity: Option<u32>,
    /// Mean over steps of `clamp_activation_fraction`.
    pub clamp_activation_fraction_mean: Option<f64>,
    /// Max over steps of `clamp_activation_fraction`.
    pub clamp_activation_fraction_max: Option<f64>,
    /// Relative residual per issue #89 D5:
    /// `|Δmass_observed − ∫Q − ∫clamp_flux| / max(|∫Q|+|∫clamp_flux|, 1)`.
    pub mass_balance_residual: Option<f64>,
    /// Integrated physical source/sink flux `Σ_steps dt·Σ_cells Q(cell,t)`.
    pub q_integral: Option<f64>,
    /// Integrated artificial clamp flux
    /// `Σ_steps dt·Σ_cells (S̃_post_clamp − S̃_pre_clamp)`.
    pub clamp_flux_integral: Option<f64>,
    /// `max |∇S̃|` on the interface cells (oceanic cells adjacent to
    /// continental, or vice versa) at the final timestep. Monitoring
    /// for issue #78.
    pub max_grad_s_interface_final: Option<f64>,
    /// `peak |f_GPE|` on the interface cells at the final timestep.
    /// Companion to `max_grad_s_interface_final`.
    pub peak_f_gpe_interface_final: Option<f64>,
    /// `max |∇S̃|` globally at the final timestep. Reference value
    /// for the interface number.
    pub max_grad_s_global_final: Option<f64>,
    /// `peak |f_GPE|` globally at the final timestep.
    pub peak_f_gpe_global_final: Option<f64>,
    /// Calibrated value of `k_spread`. `None` when calibration was
    /// not run (e.g., Step 5 regression, boundary Disabled, or when
    /// the CLI supplies an explicit `--k-spread`).
    pub k_spread_calibrated: Option<f64>,

    // ---- Step 6 — Voronoi + dynamic detection + Closed recycling ----
    //
    // All `None` outside Step 6 or when the relevant path is not
    // exercised. `plate_count` / `plate_type_distribution` are
    // populated whenever boundary is Enabled (for both static and
    // Voronoi geometries). Recycling / buffer fields are populated
    // only under `RecyclingModeInit::Closed`.
    pub plate_count: Option<u32>,
    pub plate_type_distribution: Option<(f64, f64)>, // (oceanic_frac, continental_frac)
    /// Time series of the fraction of cells whose `boundary_flag`
    /// changed vs the previous step. Populated under dynamic
    /// geometries (Voronoi) — for static geometries it would be
    /// identically zero, so we leave it `None` to signal "not
    /// applicable".
    pub boundary_flag_transition_rate_mean: Option<f64>,
    pub boundary_flag_transition_rate_max: Option<f64>,
    /// `recycling_buffer_fill` diagnostic: mean and max of the
    /// buffer's in-transit mass over the run.
    pub recycling_buffer_fill_mean: Option<f64>,
    pub recycling_buffer_fill_max: Option<f64>,
    /// Max observed `max(arc_pending, coll_v_pending, rift_v_pending)`
    /// over the run. Non-zero means some class had no eligible cell
    /// at some step and rolled over.
    pub immediate_pending_max: Option<f64>,
    /// Final immediate accumulator sum (a component of the
    /// mass-conservation residual).
    pub immediate_pending_final: Option<f64>,
    /// Final buffer fill (in-transit mass at end of run).
    pub recycling_buffer_fill_final: Option<f64>,
    /// Integrated mantle loss over the run
    /// (`Σ_steps mantle_loss_fraction · M_sub_step`). Zero when
    /// `mantle_loss_fraction = 0`.
    pub mantle_loss_integral: Option<f64>,
    /// Total subducted mass over the run — the denominator for the
    /// mantle_loss observed-fraction diagnostic.
    pub m_sub_total: Option<f64>,
    /// Closed-mode absolute mass-conservation residual per the
    /// Step 6 bilan:
    ///   `|Δmass_obs + mantle_loss_integral + buffer_fill_final +
    ///     pending_final - clamp_flux_integral| / initial_mass`
    /// Distinct from Step 5's `mass_balance_residual`; the Step 6
    /// version includes the buffer + pending terms. Target < 1e-6
    /// when `mantle_loss_fraction = 0`.
    pub mass_conservation_residual: Option<f64>,
    /// Max `clamp_activation_fraction` during the buffer spin-up
    /// window (first `mantle_delay_steps` steps). Zero if spin-up
    /// is safe.
    pub clamp_activation_during_spinup_max: Option<f64>,
    /// #78 trajectory samples: `(step, max|∇S̃|_interface,
    /// max|∇S̃|_global, peak|f_GPE|_interface, peak|f_GPE|_global,
    /// buffer_fill)` at steps {1, 10, 50, 150, 300}.
    pub issue_78_trajectory: Vec<(usize, f64, f64, f64, f64, f64)>,
    /// Per-flag-type cell count at the **final** step. Keys:
    /// `(none, subduction, oceanic_subduction, rift, continental_collision)`.
    /// Informative when `boundary_type_diversity` is suspicious —
    /// e.g. `diversity = 1` could mean either "only subduction
    /// detected" or "only rift detected"; this breakdown
    /// disambiguates.
    pub boundary_flag_counts_final: Option<(usize, usize, usize, usize, usize)>,
    /// Per-flag-type cell count at the **first** step (post-first
    /// `detect_boundaries` call). Compared against `_final` to show
    /// whether detection produced flags at step 1 and how they
    /// evolved. For static geometries, `_step1` and `_final` are
    /// identical (no dynamic update).
    pub boundary_flag_counts_step1: Option<(usize, usize, usize, usize, usize)>,
    /// Integrated arc return budget: `Σ_steps (arc_fraction · M_sub_step
    /// actually distributed this step)`. Equal to `arc_fraction ·
    /// M_sub_total` in the steady state where eligible cells always
    /// exist; deviates during rollover.
    pub arc_distributed_integral: Option<f64>,
    /// Integrated coll_v return budget. Same interpretation.
    pub coll_v_distributed_integral: Option<f64>,
    /// Integrated rift_v return budget. Same interpretation.
    pub rift_v_distributed_integral: Option<f64>,
    /// Integrated spread return emitted on rift oceanic cells.
    pub spread_distributed_integral: Option<f64>,
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
