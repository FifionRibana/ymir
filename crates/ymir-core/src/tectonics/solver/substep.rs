//! Support types for adaptive time-stepping (issue #52).
//!
//! A macro tectonic step targets a fixed geological duration `dt_target`,
//! consumed via one or more sub-steps of adaptive size. This module holds
//! the types that bridge the sub-step loop with the rest of the solver:
//!
//! * [`SubstepResult`] — per-sub-step outcome (Newton outcome + clamp ratio
//!   + iteration counts) with an `is_successful` helper.
//! * [`AccumulatedSubstepStats`] — running totals across all sub-steps of
//!   a macro step, used to build the final [`StepStats`].
//! * [`StateSnapshot`] — rollback buffer capturing every piece of mutable
//!   state a sub-step can touch, so a failed attempt can be undone and
//!   retried with a smaller dt.
//! * [`grow_dt`] / [`reduction_for`] — pure helpers that implement the
//!   per-outcome growth/reduction policy from [`AdaptiveDtConfig`].
//!
//! Nothing here drives the solver yet; wiring is done in a follow-up.
//!
//! [`StepStats`]: super::workspace::StepStats
//! [`AdaptiveDtConfig`]: super::config::AdaptiveDtConfig

use super::config::AdaptiveDtConfig;
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::newton::NewtonOutcome;
use super::tectonics::DynamicPlateContext;
use super::traction::TractionField;
use crate::tectonics::plates::Plate;

/// Outcome of a single sub-step: Newton outcome plus the clamp ratio that
/// results from the subsequent advection. A sub-step is considered a
/// success only if Newton converged AND advective clamping stayed below
/// the 5% threshold — otherwise the caller rolls back and retries with a
/// smaller dt.
#[derive(Debug, Clone, Copy)]
pub struct SubstepResult {
    pub outcome: NewtonOutcome,
    pub newton_iterations: usize,
    pub linear_iterations: usize,
    pub clamp_ratio: f64,
    pub max_velocity: f64,
}

impl SubstepResult {
    /// Returns true when this sub-step should be committed.
    ///
    /// Both conditions must hold:
    /// * Newton reached a converged outcome
    ///   ([`NewtonOutcome::ConvergedOnResidual`] or
    ///   [`NewtonOutcome::ConvergedOnState`]).
    /// * Advective clamping stayed below 5% of cells — otherwise the
    ///   sub-step was too aggressive for the CFL regime and must be retried.
    pub fn is_successful(&self) -> bool {
        matches!(
            self.outcome,
            NewtonOutcome::ConvergedOnResidual | NewtonOutcome::ConvergedOnState
        ) && self.clamp_ratio < 0.05
    }
}

/// Running totals across all sub-steps of a macro step. Cleared at the
/// start of each macro step and merged into [`StepStats`] at the end.
///
/// [`StepStats`]: super::workspace::StepStats
#[derive(Default, Debug, Clone)]
pub struct AccumulatedSubstepStats {
    pub dt_total: f64,
    pub newton_iters_total: usize,
    pub linear_iters_total: usize,
    pub max_velocity: f64,
    pub max_clamp_ratio: f64,
    pub substep_count: usize,
}

impl AccumulatedSubstepStats {
    pub fn merge(&mut self, result: &SubstepResult, dt_sub: f64) {
        self.dt_total += dt_sub;
        self.newton_iters_total += result.newton_iterations;
        self.linear_iters_total += result.linear_iterations;
        self.max_velocity = self.max_velocity.max(result.max_velocity);
        self.max_clamp_ratio = self.max_clamp_ratio.max(result.clamp_ratio);
        self.substep_count += 1;
    }
}

/// Per-outcome reduction factor applied to `dt_current` after a failed
/// sub-step. Mild for Stagnation, moderate for Oscillation, aggressive
/// for Divergence, and a shared default for the other failure modes.
///
/// The two successful outcomes (`ConvergedOnResidual`, `ConvergedOnState`)
/// fall through to `default_reduction`, but this helper is only called on
/// failed sub-steps so that branch never runs in practice.
pub fn reduction_for(result: &SubstepResult, cfg: &AdaptiveDtConfig) -> f64 {
    match result.outcome {
        NewtonOutcome::Stagnation => cfg.stagnation_reduction,
        NewtonOutcome::Oscillation => cfg.oscillation_reduction,
        NewtonOutcome::Divergence => cfg.divergence_reduction,
        _ => cfg.default_reduction,
    }
}

/// Post-success growth factor for `dt_current`. Easy wins (few Newton
/// iterations) grow `dt` aggressively; normal wins grow it modestly;
/// hard wins (near `normal_iters`) leave it unchanged so we don't try
/// to speed up when Newton is already struggling.
pub fn grow_dt(dt_current: f64, result: &SubstepResult, cfg: &AdaptiveDtConfig) -> f64 {
    if result.newton_iterations < cfg.easy_iters {
        dt_current * cfg.easy_growth
    } else if result.newton_iterations < cfg.normal_iters {
        dt_current * cfg.normal_growth
    } else {
        dt_current
    }
}

/// Snapshot of every piece of mutable solver state touched by a sub-step.
///
/// A failed sub-step restores this snapshot so the next attempt starts
/// from the exact pre-sub-step state, independent of whatever partial
/// work the failed attempt did (Newton velocity field, clamped advection,
/// plastic strain accumulation, boundary-source side effects, etc.).
///
/// The snapshot is taken inside the sub-step loop, not at macro-step
/// start, so successful sub-steps are preserved and only the current
/// attempt is rolled back.
pub struct StateSnapshot {
    pub s: Field2D,
    pub vx: Field2D,
    pub vy: Field2D,
    pub rho: Field2D,
    pub eta_multiplier: Field2D,
    pub plastic_strain: Field2D,
    pub basal_friction: f64,
    pub plate_ids: Vec<usize>,
    pub plates: Vec<Plate>,
    pub traction: TractionField,
    pub disp_x: Field2D,
    pub disp_y: Field2D,
    pub next_id: usize,
}

impl StateSnapshot {
    /// Capture the full mutable state that a sub-step can mutate.
    pub fn capture(grid: &StaggeredGrid, plate_ctx: &DynamicPlateContext) -> Self {
        Self {
            s: grid.s.clone(),
            vx: grid.vx.clone(),
            vy: grid.vy.clone(),
            rho: grid.rho.clone(),
            eta_multiplier: grid.eta_multiplier.clone(),
            plastic_strain: grid.plastic_strain.clone(),
            basal_friction: grid.basal_friction,
            plate_ids: plate_ctx.ids.clone(),
            plates: plate_ctx.plates.clone(),
            traction: plate_ctx.traction.clone(),
            disp_x: plate_ctx.disp_x.clone(),
            disp_y: plate_ctx.disp_y.clone(),
            next_id: plate_ctx.next_id,
        }
    }

    /// Overwrite `grid` and `plate_ctx` with the captured state.
    pub fn restore(&self, grid: &mut StaggeredGrid, plate_ctx: &mut DynamicPlateContext) {
        grid.s.data_mut().copy_from_slice(self.s.data());
        grid.vx.data_mut().copy_from_slice(self.vx.data());
        grid.vy.data_mut().copy_from_slice(self.vy.data());
        grid.rho.data_mut().copy_from_slice(self.rho.data());
        grid.eta_multiplier
            .data_mut()
            .copy_from_slice(self.eta_multiplier.data());
        grid.plastic_strain
            .data_mut()
            .copy_from_slice(self.plastic_strain.data());
        grid.basal_friction = self.basal_friction;
        plate_ctx.ids.clone_from(&self.plate_ids);
        plate_ctx.plates.clone_from(&self.plates);
        plate_ctx.traction = self.traction.clone();
        plate_ctx
            .disp_x
            .data_mut()
            .copy_from_slice(self.disp_x.data());
        plate_ctx
            .disp_y
            .data_mut()
            .copy_from_slice(self.disp_y.data());
        plate_ctx.next_id = self.next_id;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::plates::PlateType;

    fn fake_result(outcome: NewtonOutcome, iters: usize, clamp: f64) -> SubstepResult {
        SubstepResult {
            outcome,
            newton_iterations: iters,
            linear_iterations: iters * 10,
            clamp_ratio: clamp,
            max_velocity: 0.001,
        }
    }

    #[test]
    fn is_successful_requires_convergence_and_low_clamp() {
        assert!(fake_result(NewtonOutcome::ConvergedOnResidual, 4, 0.01).is_successful());
        assert!(fake_result(NewtonOutcome::ConvergedOnState, 10, 0.03).is_successful());
        // Newton converged but advection clamped too many cells
        assert!(!fake_result(NewtonOutcome::ConvergedOnResidual, 4, 0.12).is_successful());
        // Newton failed
        assert!(!fake_result(NewtonOutcome::Stagnation, 15, 0.01).is_successful());
        assert!(!fake_result(NewtonOutcome::Oscillation, 7, 0.01).is_successful());
        assert!(!fake_result(NewtonOutcome::Divergence, 9, 0.01).is_successful());
        assert!(!fake_result(NewtonOutcome::MaxIterations, 15, 0.01).is_successful());
    }

    #[test]
    fn reduction_for_matches_outcome_specific_factor() {
        let cfg = AdaptiveDtConfig::default();
        assert_eq!(
            reduction_for(&fake_result(NewtonOutcome::Stagnation, 15, 0.01), &cfg),
            cfg.stagnation_reduction
        );
        assert_eq!(
            reduction_for(&fake_result(NewtonOutcome::Oscillation, 7, 0.01), &cfg),
            cfg.oscillation_reduction
        );
        assert_eq!(
            reduction_for(&fake_result(NewtonOutcome::Divergence, 9, 0.01), &cfg),
            cfg.divergence_reduction
        );
        assert_eq!(
            reduction_for(&fake_result(NewtonOutcome::MaxIterations, 15, 0.01), &cfg),
            cfg.default_reduction
        );
    }

    #[test]
    fn all_reduction_factors_are_strict_contractions() {
        // Invariant the sub-step loop relies on: every reduction factor is
        // strictly below 1.0, so dt_current decreases on every failed
        // attempt and the loop cannot spin indefinitely above the floor.
        let cfg = AdaptiveDtConfig::default();
        for r in [
            cfg.stagnation_reduction,
            cfg.oscillation_reduction,
            cfg.divergence_reduction,
            cfg.default_reduction,
        ] {
            assert!(r > 0.0 && r < 1.0, "reduction factor {r} must be in (0, 1)");
        }
    }

    #[test]
    fn grow_dt_scales_with_difficulty() {
        let cfg = AdaptiveDtConfig::default();
        let easy = fake_result(NewtonOutcome::ConvergedOnResidual, 2, 0.01);
        let normal = fake_result(NewtonOutcome::ConvergedOnResidual, 5, 0.01);
        let hard = fake_result(NewtonOutcome::ConvergedOnState, 12, 0.01);

        let after_easy = grow_dt(1.0, &easy, &cfg);
        let after_normal = grow_dt(1.0, &normal, &cfg);
        let after_hard = grow_dt(1.0, &hard, &cfg);

        assert_eq!(after_easy, cfg.easy_growth);
        assert_eq!(after_normal, cfg.normal_growth);
        assert_eq!(after_hard, 1.0);
        assert!(after_easy > after_normal);
        assert!(after_normal > after_hard);
    }

    #[test]
    fn accumulated_stats_track_totals_and_maxes() {
        let mut acc = AccumulatedSubstepStats::default();
        acc.merge(
            &SubstepResult {
                outcome: NewtonOutcome::ConvergedOnResidual,
                newton_iterations: 4,
                linear_iterations: 30,
                clamp_ratio: 0.02,
                max_velocity: 0.01,
            },
            1.5,
        );
        acc.merge(
            &SubstepResult {
                outcome: NewtonOutcome::ConvergedOnState,
                newton_iterations: 7,
                linear_iterations: 45,
                clamp_ratio: 0.04,
                max_velocity: 0.008,
            },
            0.5,
        );
        assert_eq!(acc.substep_count, 2);
        assert!((acc.dt_total - 2.0).abs() < 1e-12);
        assert_eq!(acc.newton_iters_total, 11);
        assert_eq!(acc.linear_iters_total, 75);
        assert!((acc.max_velocity - 0.01).abs() < 1e-12);
        assert!((acc.max_clamp_ratio - 0.04).abs() < 1e-12);
    }

    fn make_minimal_ctx(nx: usize, ny: usize) -> DynamicPlateContext {
        DynamicPlateContext {
            ids: vec![0; nx * ny],
            plates: vec![Plate {
                id: 0,
                plate_type: PlateType::Continental,
                velocity: (0.1, -0.2),
                seed_x: 1.5,
                seed_y: 2.5,
                active: true,
                subducted_mass: 3.0,
                cell_count: nx * ny,
                mean_thickness: 1.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            }],
            traction: TractionField::uniform(nx, ny, 0.3, -0.1),
            next_id: 7,
            disp_x: Field2D::filled(nx, ny, 0.01),
            disp_y: Field2D::filled(nx, ny, -0.02),
        }
    }

    #[test]
    fn state_snapshot_roundtrip_is_lossless() {
        let nx = 4;
        let ny = 4;
        let dx = 0.25;
        let mut grid = StaggeredGrid::new(nx, ny, dx);
        grid.s.set(1, 2, 1.5);
        grid.vx.set(2, 1, 0.2);
        grid.vy.set(0, 3, -0.1);
        grid.rho.set(2, 2, 2800.0);
        grid.eta_multiplier.set(3, 3, 5.0);
        grid.plastic_strain.set(1, 1, 0.05);
        grid.basal_friction = 0.15;

        let mut plate_ctx = make_minimal_ctx(nx, ny);

        let snap = StateSnapshot::capture(&grid, &plate_ctx);

        // Scramble everything the snapshot should restore.
        for val in grid.s.data_mut() {
            *val = 999.0;
        }
        for val in grid.vx.data_mut() {
            *val = -999.0;
        }
        for val in grid.vy.data_mut() {
            *val = -999.0;
        }
        for val in grid.rho.data_mut() {
            *val = 0.0;
        }
        for val in grid.eta_multiplier.data_mut() {
            *val = 0.0;
        }
        for val in grid.plastic_strain.data_mut() {
            *val = 99.0;
        }
        grid.basal_friction = 12.0;
        plate_ctx.ids[0] = 42;
        plate_ctx.plates[0].active = false;
        plate_ctx.plates[0].subducted_mass = 999.0;
        plate_ctx.next_id = 999;
        for val in plate_ctx.disp_x.data_mut() {
            *val = 1.0;
        }
        for val in plate_ctx.disp_y.data_mut() {
            *val = 1.0;
        }

        snap.restore(&mut grid, &mut plate_ctx);

        assert_eq!(grid.s.get(1, 2), 1.5);
        assert_eq!(grid.vx.get(2, 1), 0.2);
        assert_eq!(grid.vy.get(0, 3), -0.1);
        assert_eq!(grid.rho.get(2, 2), 2800.0);
        assert_eq!(grid.eta_multiplier.get(3, 3), 5.0);
        assert_eq!(grid.plastic_strain.get(1, 1), 0.05);
        assert!((grid.basal_friction - 0.15).abs() < 1e-12);
        assert_eq!(plate_ctx.ids[0], 0);
        assert!(plate_ctx.plates[0].active);
        assert!((plate_ctx.plates[0].subducted_mass - 3.0).abs() < 1e-12);
        assert_eq!(plate_ctx.next_id, 7);
        assert!((plate_ctx.disp_x.get(0, 0) - 0.01).abs() < 1e-12);
        assert!((plate_ctx.disp_y.get(0, 0) - -0.02).abs() < 1e-12);
    }
}
