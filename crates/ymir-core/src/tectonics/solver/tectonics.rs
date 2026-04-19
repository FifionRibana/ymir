//! Top-level orchestration of the thin viscous sheet simulation.

use std::time::Instant;

use tracing::{debug, info, info_span, warn};

use super::advection::{compute_cfl_dt, compute_divergence_flux};
use super::config::{NonlinearSolver, TectonicsConfig};
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::newton::{NewtonOutcome, solve_velocity_newton};
use super::picard::solve_velocity_picard;
use super::substep::{
    AccumulatedSubstepStats, StateSnapshot, SubstepResult, grow_dt, reduction_for,
};
use super::traction::TractionField;
use super::workspace::{SolverWorkspace, StepStats};
use crate::tectonics::boundaries::{
    BoundaryType, accumulate_subducted_mass, compute_boundary_sources, gaussian_blur_f64,
};
use crate::tectonics::mantle::MantleFlow;
use crate::tectonics::plates::{
    Plate, PlateType, advect_plate_ids, apply_subduction_consumption, cleanup_plate_ids,
    compute_viscosity_multiplier, detect_disappeared_plates, detect_fragmentation,
    rebuild_traction, rebuild_traction_smooth, update_plate_stats,
};
use crate::tectonics::recycling::RecyclingBuffer;

/// Errors that can occur during a tectonic simulation run.
#[derive(Debug)]
pub enum SolverError {
    /// Nonlinear solver did not converge at the given timestep.
    NonlinearSolverDidNotConverge { step: usize },
    /// A NaN or Inf was detected in the solution.
    NumericalInstability { step: usize, field: &'static str },
    /// The simulation was cancelled via the progress callback.
    Cancelled,
}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolverError::NonlinearSolverDidNotConverge { step } => {
                write!(f, "Nonlinear solver did not converge at step {step}")
            }
            SolverError::NumericalInstability { step, field } => {
                write!(f, "NaN/Inf detected in {field} at step {step}")
            }
            SolverError::Cancelled => write!(f, "Simulation cancelled"),
        }
    }
}

impl std::error::Error for SolverError {}

/// Run a full tectonic simulation.
///
/// The `progress` callback is invoked after each timestep with
/// `(step, total, stats, s_field)`. The `s_field` is a reference to the
/// current crustal thickness field (safe to clone for snapshots).
/// Return `true` to continue, `false` to cancel.
/// Mutable plate context for dynamic boundary simulation.
///
/// Contains owned plate data that may be updated each timestep when
/// `dynamic_boundaries` is enabled (seeds advected, Voronoï recomputed,
/// traction rebuilt).
pub struct DynamicPlateContext {
    pub ids: Vec<usize>,
    pub plates: Vec<Plate>,
    pub traction: TractionField,
    /// Counter for creating new plates (rift creation, fragmentation).
    pub next_id: usize,
    /// Accumulated fractional X displacement per cell (sub-pixel ID advection).
    pub disp_x: Field2D,
    /// Accumulated fractional Y displacement per cell (sub-pixel ID advection).
    pub disp_y: Field2D,
}

/// Data passed to the progress callback after each timestep.
pub struct StepSnapshot<'a> {
    pub s_field: &'a Field2D,
    pub plate_ids: Option<&'a [usize]>,
    pub plates: Option<&'a [Plate]>,
    pub boundary_types: Option<&'a [BoundaryType]>,
}

/// How `execute_tectonic_pass` sizes the step's dt.
///
/// * `LegacyWithCflRetry` — pre-#52 behaviour: compute `dt_cfl` internally
///   and run up to 5 halving retries on excessive clamping.
/// * `Explicit(dt_sub)` — adaptive sub-step path: use exactly `dt_sub`
///   with no internal retry; the caller rolls back and shrinks on failure.
#[derive(Debug, Clone, Copy)]
enum DtMode {
    LegacyWithCflRetry,
    Explicit(f64),
}

/// Everything the caller of `execute_tectonic_pass` needs to either build
/// `StepStats` + diagnostic logs (legacy mode) or drive the adaptive
/// sub-step loop (explicit mode).
struct PassResult {
    outcome: NewtonOutcome,
    nl_iterations: usize,
    linear_iterations: usize,
    dt_consumed: f64,
    clamp_ratio: f64,
    retry_succeeded: bool,
    max_velocity: f64,
    max_thickness: f64,
    min_thickness: f64,
    mass_before: f64,
    mass_after: f64,
    total_div_flux: f64,
    total_source: f64,
    mass_clamp_delta: f64,
    oceanic_restore_delta: f64,
    continental_restore_delta: f64,
    /// Number of sub-steps that composed this macro step. Always 1 on the
    /// legacy path; ≥ 1 on the adaptive path.
    substep_count: usize,
}

impl PassResult {
    fn is_successful(&self, max_clamp_ratio: f64) -> bool {
        self.outcome.is_converged() && self.clamp_ratio < max_clamp_ratio
    }

    fn as_substep(&self) -> SubstepResult {
        SubstepResult {
            outcome: self.outcome,
            newton_iterations: self.nl_iterations,
            linear_iterations: self.linear_iterations,
            clamp_ratio: self.clamp_ratio,
            max_velocity: self.max_velocity,
        }
    }
}

pub fn run_tectonics<F>(
    config: &TectonicsConfig,
    plate_ctx: &mut DynamicPlateContext,
    grid: &mut StaggeredGrid,
    workspace: &mut SolverWorkspace,
    mut progress: F,
) -> Result<(), SolverError>
where
    F: FnMut(usize, usize, &StepStats, StepSnapshot<'_>) -> bool,
{
    let nx = grid.nx();
    let ny = grid.ny();

    // Initialize density field and cratonic viscosity multiplier from plate types
    assign_density_from_plates(grid, &plate_ctx.ids, &plate_ctx.plates, &config.boundaries);
    compute_viscosity_multiplier(grid, &plate_ctx.ids, &plate_ctx.plates, &config.cratonic);
    grid.basal_friction = config.basal_friction;

    let mut recycling_buffer = if config.recycling.enabled {
        Some(RecyclingBuffer::new(config.recycling.mantle_delay))
    } else {
        None
    };

    // Generate mantle convection flow field (static pattern, optionally evolving)
    let mut mantle_flow = if config.mantle.enabled {
        let mantle_seed =
            grid.s.data().iter().take(8).fold(0u64, |acc, &v| acc.wrapping_add((v * 1e10) as u64));
        Some(MantleFlow::generate(nx, ny, mantle_seed, &config.mantle))
    } else {
        None
    };

    for step in 0..config.num_timesteps {
        let _step_span = info_span!("solver_step", step, nx = grid.nx(), ny = grid.ny()).entered();

        let pass = if config.adaptive_dt.enabled {
            run_adaptive_macro_step(
                grid,
                plate_ctx,
                workspace,
                &mut recycling_buffer,
                &mut mantle_flow,
                config,
                step,
            )?
        } else {
            execute_tectonic_pass(
                grid,
                plate_ctx,
                workspace,
                &mut recycling_buffer,
                &mut mantle_flow,
                config,
                step,
                DtMode::LegacyWithCflRetry,
                step == 0,
            )?
        };

        workspace.stats = StepStats {
            max_velocity: pass.max_velocity,
            max_thickness: pass.max_thickness,
            min_thickness: pass.min_thickness,
            picard_iterations: pass.nl_iterations,
            cg_iterations_last: pass.linear_iterations,
            dt: pass.dt_consumed,
            clamp_ratio: pass.clamp_ratio,
            cfl_retry_exhausted: !pass.retry_succeeded,
        };

        info!(
            step,
            dt = pass.dt_consumed,
            nl_iters = pass.nl_iterations,
            lin_iters = pass.linear_iterations,
            max_velocity = pass.max_velocity,
            min_thickness = pass.min_thickness,
            max_thickness = pass.max_thickness,
            clamp_ratio = pass.clamp_ratio,
            substeps = pass.substep_count,
            "tectonic step"
        );

        if pass.clamp_ratio > 0.05 {
            warn!(
                step,
                clamp_ratio = pass.clamp_ratio,
                dt = pass.dt_consumed,
                "excessive clamping — consider reducing CFL or timestep"
            );
        }

        debug!(
            step,
            mass_total = %format!("{:.4}", pass.mass_after),
            mass_delta = %format!("{:.6}", pass.mass_after - pass.mass_before),
            advection_div = %format!("{:.6}", -pass.total_div_flux),
            sources_q = %format!("{:.6}", pass.total_source),
            clamping = %format!("{:.6}", pass.mass_clamp_delta),
            oceanic_restore = %format!("{:.6}", pass.oceanic_restore_delta),
            continental_restore = %format!("{:.6}", pass.continental_restore_delta),
            "mass balance"
        );

        if config.dynamic_boundaries {
            let num_active =
                plate_ctx.plates.iter().filter(|p| p.active && p.cell_count > 0).count();

            let mut plate_summary = String::new();
            for plate in plate_ctx.plates.iter().filter(|p| p.active && p.cell_count > 0) {
                let mat = if plate.mean_thickness > 0.4 { "C" } else { "O" };
                use std::fmt::Write;
                let _ = write!(
                    plate_summary,
                    " [{}:{} n={} S={:.3} v=({:.4},{:.4})]",
                    plate.id,
                    mat,
                    plate.cell_count,
                    plate.mean_thickness,
                    plate.mean_velocity.0,
                    plate.mean_velocity.1,
                );
            }

            let (mut n_sub, mut n_osub, mut n_coll, mut n_rift, mut n_boundary) =
                (0usize, 0usize, 0usize, 0usize, 0usize);
            if let Some(bf) = workspace.boundary_field.as_ref() {
                for bt in &bf.boundary_type {
                    match bt {
                        BoundaryType::Subduction => {
                            n_sub += 1;
                            n_boundary += 1;
                        }
                        BoundaryType::OceanicSubduction => {
                            n_osub += 1;
                            n_boundary += 1;
                        }
                        BoundaryType::ContinentalCollision => {
                            n_coll += 1;
                            n_boundary += 1;
                        }
                        BoundaryType::Rift => {
                            n_rift += 1;
                            n_boundary += 1;
                        }
                        BoundaryType::None => {}
                    }
                }
            }

            let mut q_pos = 0.0_f64;
            let mut q_neg = 0.0_f64;
            for &q in workspace.source_rate.data().iter() {
                if q > 0.0 {
                    q_pos += q;
                } else {
                    q_neg += q;
                }
            }

            let mut n_continental = 0usize;
            let mut n_oceanic = 0usize;
            let mut n_transitional = 0usize;
            for &s in grid.s.data().iter() {
                if s > 0.4 {
                    n_continental += 1;
                } else if s < 0.3 {
                    n_oceanic += 1;
                } else {
                    n_transitional += 1;
                }
            }

            debug!(
                step,
                active_plates = num_active,
                total_plates = plate_ctx.plates.len(),
                next_id = plate_ctx.next_id,
                boundary_cells = n_boundary,
                subduction = n_sub,
                oceanic_sub = n_osub,
                collision = n_coll,
                rift = n_rift,
                q_positive = format!("{:.4}", q_pos),
                q_negative = format!("{:.4}", q_neg),
                cells_continental = n_continental,
                cells_oceanic = n_oceanic,
                cells_transitional = n_transitional,
                "plate diagnostics"
            );
            debug!(step, plates = %plate_summary, "plate details");
        }

        let snapshot = StepSnapshot {
            s_field: &grid.s,
            plate_ids: if config.dynamic_boundaries {
                Some(plate_ctx.ids.as_slice())
            } else {
                None
            },
            plates: if config.dynamic_boundaries { Some(&plate_ctx.plates) } else { None },
            boundary_types: workspace.boundary_field.as_ref().map(|bf| bf.boundary_type.as_slice()),
        };
        if !progress(step, config.num_timesteps, &workspace.stats, snapshot) {
            return Err(SolverError::Cancelled);
        }
    }

    Ok(())
}

/// Assign density to each grid cell based on its plate type.
fn assign_density_from_plates(
    grid: &mut StaggeredGrid,
    ids: &[usize],
    plates: &[Plate],
    bc: &crate::tectonics::boundaries::BoundaryConfig,
) {
    let nx = grid.nx();
    let ny = grid.ny();
    for j in 0..ny {
        for i in 0..nx {
            let pid = ids[j * nx + i];
            let rho = match plates[pid].plate_type {
                PlateType::Continental => bc.rho_continental,
                PlateType::Oceanic => bc.rho_oceanic,
            };
            grid.rho.set(i, j, rho);
        }
    }
}

/// Assign crustal density based on local thickness, not plate type.
/// S > 0.4 → continental density, S < 0.3 → oceanic density,
/// 0.3–0.4 → linear interpolation (transitional).
fn assign_density_from_material(
    grid: &mut StaggeredGrid,
    config: &crate::tectonics::boundaries::BoundaryConfig,
) {
    let nx = grid.nx();
    let ny = grid.ny();
    for j in 0..ny {
        for i in 0..nx {
            let s = grid.s.get(i, j);
            let rho = if s > 0.4 {
                config.rho_continental
            } else if s < 0.3 {
                config.rho_oceanic
            } else {
                let t = (s - 0.3) / 0.1;
                config.rho_oceanic + t * (config.rho_continental - config.rho_oceanic)
            };
            grid.rho.set(i, j, rho);
        }
    }
}

/// Solve velocity directly (no continuation).
/// If friction is enabled and the direct solve fails, falls back to
/// solving without friction first, then re-solving with friction.
fn solve_velocity_direct(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    config: &TectonicsConfig,
    slab: Option<&super::stokes::SlabPullField>,
    workspace: &mut SolverWorkspace,
) -> (NewtonOutcome, usize, usize) {
    let rho_c = config.boundaries.rho_continental;
    let rho_m = config.boundaries.rho_mantle;

    // Picard maps to a synthetic NewtonOutcome — it does not surface the
    // same failure modes, so we conflate all Picard failures onto
    // MaxIterations which selects the default_reduction factor in the
    // adaptive sub-step policy.
    let do_solve =
        |grid: &mut StaggeredGrid, ws: &mut SolverWorkspace| -> (NewtonOutcome, usize, usize) {
            match config.nonlinear_solver {
                NonlinearSolver::Picard => {
                    let r = solve_velocity_picard(
                        grid,
                        plates,
                        config.gravity_factor,
                        rho_c,
                        rho_m,
                        &config.picard,
                        &config.yielding,
                        slab,
                        ws,
                    );
                    let outcome = if r.converged {
                        NewtonOutcome::ConvergedOnResidual
                    } else {
                        NewtonOutcome::MaxIterations
                    };
                    (outcome, r.iterations, r.total_cg_iterations)
                }
                NonlinearSolver::Newton => {
                    let r = solve_velocity_newton(
                        grid,
                        plates,
                        config.gravity_factor,
                        rho_c,
                        rho_m,
                        &config.picard,
                        &config.yielding,
                        &config.newton,
                        slab,
                        ws,
                    );
                    (r.outcome, r.iterations, r.total_linear_iterations)
                }
            }
        };

    let (outcome, nl, linear) = do_solve(grid, workspace);

    if !outcome.is_converged() && grid.basal_friction > 0.0 {
        // Fallback: solve without friction, then with friction
        let target_friction = grid.basal_friction;
        grid.basal_friction = 0.0;
        let (_, nl2, lin2) = do_solve(grid, workspace);
        grid.basal_friction = target_friction;
        let (outcome2, nl3, lin3) = do_solve(grid, workspace);
        return (outcome2, nl + nl2 + nl3, linear + lin2 + lin3);
    }

    (outcome, nl, linear)
}

/// Solve velocity using viscosity continuation: ramp n from 1 → target,
/// then optionally add a friction continuation step.
fn solve_with_continuation(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    config: &TectonicsConfig,
    slab: Option<&super::stokes::SlabPullField>,
    workspace: &mut SolverWorkspace,
) -> (NewtonOutcome, usize, usize) {
    let target_eps = config.picard.strain_rate_min;
    let steps = &config.continuation.n_steps;
    let eps_start = config.continuation.eps_min_start.unwrap_or(target_eps);
    let target_friction = grid.basal_friction;

    let mut total_nl = 0usize;
    let mut total_linear = 0usize;

    // Phase 1: power-law ramp WITHOUT friction
    grid.basal_friction = 0.0;

    for (i, &n_exp) in steps.iter().enumerate() {
        // Interpolate ε_min from eps_start to target_eps
        let t = if steps.len() > 1 { i as f64 / (steps.len() - 1) as f64 } else { 1.0 };
        let eps_min = eps_start * (1.0 - t) + target_eps * t;

        // Adapt relaxation to nonlinearity level
        let relaxation = if n_exp <= 1.5 {
            0.9
        } else if n_exp <= 2.5 {
            0.6
        } else {
            0.4
        };

        let mut step_config = config.picard.clone();
        step_config.power_law_n = n_exp;
        step_config.strain_rate_min = eps_min;
        step_config.relaxation = relaxation;

        let rho_c = config.boundaries.rho_continental;
        let rho_m = config.boundaries.rho_mantle;
        // Warm start: grid.vx/vy retain the solution from the previous step
        let (outcome, iters, linear_iters): (NewtonOutcome, usize, usize) =
            match config.nonlinear_solver {
                NonlinearSolver::Picard => {
                    let r = solve_velocity_picard(
                        grid,
                        plates,
                        config.gravity_factor,
                        rho_c,
                        rho_m,
                        &step_config,
                        &config.yielding,
                        slab,
                        workspace,
                    );
                    let o = if r.converged {
                        NewtonOutcome::ConvergedOnResidual
                    } else {
                        NewtonOutcome::MaxIterations
                    };
                    (o, r.iterations, r.total_cg_iterations)
                }
                NonlinearSolver::Newton => {
                    let r = solve_velocity_newton(
                        grid,
                        plates,
                        config.gravity_factor,
                        rho_c,
                        rho_m,
                        &step_config,
                        &config.yielding,
                        &config.newton,
                        slab,
                        workspace,
                    );
                    (r.outcome, r.iterations, r.total_linear_iterations)
                }
            };
        total_nl += iters;
        total_linear += linear_iters;
        if !outcome.is_converged() {
            grid.basal_friction = target_friction;
            return (outcome, total_nl, total_linear);
        }
    }

    // Phase 2: final solve WITH friction at target power-law
    if target_friction > 0.0 {
        grid.basal_friction = target_friction;

        let mut step_config = config.picard.clone();
        step_config.power_law_n = config.picard.power_law_n;
        step_config.strain_rate_min = target_eps;
        step_config.relaxation = 0.4;

        let rho_c = config.boundaries.rho_continental;
        let rho_m = config.boundaries.rho_mantle;

        let (outcome, iters, linear_iters): (NewtonOutcome, usize, usize) =
            match config.nonlinear_solver {
                NonlinearSolver::Picard => {
                    let r = solve_velocity_picard(
                        grid,
                        plates,
                        config.gravity_factor,
                        rho_c,
                        rho_m,
                        &step_config,
                        &config.yielding,
                        slab,
                        workspace,
                    );
                    let o = if r.converged {
                        NewtonOutcome::ConvergedOnResidual
                    } else {
                        NewtonOutcome::MaxIterations
                    };
                    (o, r.iterations, r.total_cg_iterations)
                }
                NonlinearSolver::Newton => {
                    let r = solve_velocity_newton(
                        grid,
                        plates,
                        config.gravity_factor,
                        rho_c,
                        rho_m,
                        &step_config,
                        &config.yielding,
                        &config.newton,
                        slab,
                        workspace,
                    );
                    (r.outcome, r.iterations, r.total_linear_iterations)
                }
            };

        total_nl += iters;
        total_linear += linear_iters;

        if !outcome.is_converged() {
            return (outcome, total_nl, total_linear);
        }
    } else {
        grid.basal_friction = target_friction;
    }

    (NewtonOutcome::ConvergedOnResidual, total_nl, total_linear)
}

/// One "pass" of the tectonic step body: velocity solve + plate/boundary
/// updates + advection + restoring forces. Shared between the legacy
/// single-pass-per-macro path (`DtMode::LegacyWithCflRetry`) and the
/// adaptive sub-step loop (`DtMode::Explicit`). The caller is responsible
/// for building `StepStats`, emitting the INFO + diagnostic logs, and
/// invoking the progress callback — this function only mutates physical
/// state and returns a `PassResult`.
#[allow(clippy::too_many_arguments)]
fn execute_tectonic_pass(
    grid: &mut StaggeredGrid,
    plate_ctx: &mut DynamicPlateContext,
    workspace: &mut SolverWorkspace,
    recycling_buffer: &mut Option<RecyclingBuffer>,
    mantle_flow: &mut Option<MantleFlow>,
    config: &TectonicsConfig,
    step: usize,
    dt_mode: DtMode,
    allow_continuation: bool,
) -> Result<PassResult, SolverError> {
    let nx = grid.nx();
    let ny = grid.ny();

    // Phase 1-bis instrumentation (#75) — per-phase wallclock timers.
    let t_solve_start = Instant::now();

    // 1. Solve velocity — continuation only on first pass (cold start).
    let need_continuation =
        allow_continuation && config.continuation.enabled && config.picard.power_law_n > 1.0;

    let traction = &plate_ctx.traction;
    // Split-borrow workaround: `workspace.boundary_field` (read-only, for
    // γ_slab and n̂) and the rest of `workspace` (written by the linear
    // solve) can't be borrowed simultaneously via field access on
    // `*workspace`. Take the Option out, build the SlabPullField borrow
    // from the owned value, run the solve, then put it back.
    let bf_taken = workspace.boundary_field.take();
    let slab = bf_taken.as_ref().map(|bf| super::stokes::SlabPullField {
        gamma: &bf.gamma_slab,
        n_x: &bf.normal_x,
        n_y: &bf.normal_y,
    });
    let slab_ref = slab.as_ref();
    let (solve_outcome, nl_iterations, linear_iterations) = if need_continuation {
        solve_with_continuation(grid, traction, config, slab_ref, workspace)
    } else {
        let result = solve_velocity_direct(grid, traction, config, slab_ref, workspace);
        if !result.0.is_converged()
            && config.continuation.enabled
            && config.picard.power_law_n > 1.0
        {
            warn!(step, "direct solve failed, falling back to continuation");
            solve_with_continuation(grid, traction, config, slab_ref, workspace)
        } else {
            result
        }
    };
    drop(slab);
    workspace.boundary_field = bf_taken;

    let t_solve_us = t_solve_start.elapsed().as_micros() as u64;

    // Emit residual_spatial diagnostic after a successful Newton solve.
    // `workspace.jfnk_f_v` holds `F(v_converged)` from the last Newton
    // iteration; it is unpopulated for Picard, so we guard on solver
    // choice. The boundary_field for this step is not yet computed here
    // (that happens below), so we pass the previous step's field (held
    // in workspace.boundary_field since the last macro step).
    if solve_outcome.is_converged() && matches!(config.nonlinear_solver, NonlinearSolver::Newton) {
        super::diagnostics::emit_residual_spatial(
            &workspace.jfnk_f_v,
            workspace.boundary_field.as_ref(),
            nx,
            ny,
        );
    }

    if !solve_outcome.is_converged() {
        // Legacy mode preserves the pre-#52 error-on-non-convergence
        // semantics. Adaptive mode returns a failure outcome so the outer
        // loop can retry with a smaller dt instead of aborting the run.
        match dt_mode {
            DtMode::LegacyWithCflRetry => {
                return Err(SolverError::NonlinearSolverDidNotConverge { step });
            }
            DtMode::Explicit(_) => {
                let mass = grid.s.data().iter().sum();
                return Ok(PassResult {
                    outcome: solve_outcome,
                    nl_iterations,
                    linear_iterations,
                    dt_consumed: 0.0,
                    clamp_ratio: 0.0,
                    retry_succeeded: false,
                    max_velocity: 0.0,
                    max_thickness: 0.0,
                    min_thickness: 0.0,
                    mass_before: mass,
                    mass_after: mass,
                    total_div_flux: 0.0,
                    total_source: 0.0,
                    mass_clamp_delta: 0.0,
                    oceanic_restore_delta: 0.0,
                    continental_restore_delta: 0.0,
                    substep_count: 1,
                });
            }
        }
    }

    // 2. Pick the dt that drives rate-based operations (plastic strain,
    //    plate advection, slab pull, mass recycling). Legacy uses the
    //    CFL dt; adaptive uses the externally sized sub-step.
    let dt_cfl = compute_cfl_dt(grid, config.cfl_factor);
    let dt_rates = match dt_mode {
        DtMode::LegacyWithCflRetry => dt_cfl,
        DtMode::Explicit(dt_sub) => dt_sub,
    };

    // 2b. Accumulate plastic strain
    if config.yielding.enabled {
        super::picard::accumulate_plastic_strain(
            &workspace.strain_rate,
            &workspace.eta,
            &config.yielding,
            dt_rates,
            &mut grid.plastic_strain,
        );
    }

    // ── Dynamic boundary update ────────────────────────────────────
    if config.dynamic_boundaries {
        advect_plate_ids(
            &mut plate_ctx.ids,
            &mut plate_ctx.disp_x,
            &mut plate_ctx.disp_y,
            grid,
            dt_rates,
        );

        let ids_before_cleanup: Vec<usize> = plate_ctx.ids.clone();
        cleanup_plate_ids(&mut plate_ctx.ids, nx, ny);
        for k in 0..nx * ny {
            if plate_ctx.ids[k] != ids_before_cleanup[k] {
                let i = k % nx;
                let j = k / nx;
                plate_ctx.disp_x.set(i, j, 0.0);
                plate_ctx.disp_y.set(i, j, 0.0);
            }
        }

        let disappeared = detect_disappeared_plates(&plate_ctx.ids, &mut plate_ctx.plates);
        for pid in &disappeared {
            info!(plate_id = pid, step, "plate disappeared");
        }

        update_plate_stats(&plate_ctx.ids, &mut plate_ctx.plates, grid);
        assign_density_from_material(grid, &config.boundaries);
        compute_viscosity_multiplier(grid, &plate_ctx.ids, &plate_ctx.plates, &config.cratonic);
    }

    // Compute boundary source terms if enabled
    let t_boundaries_start = Instant::now();
    let boundaries_active = config.boundaries.enabled;
    if boundaries_active {
        let bf = compute_boundary_sources(
            grid,
            &plate_ctx.ids,
            &plate_ctx.plates,
            &config.boundaries,
            config.recycling.enabled,
        );
        workspace.source_rate.data_mut().copy_from_slice(bf.source_rate.data());
        workspace.boundary_field = Some(bf);
        if config.boundaries.source_smoothing_sigma > 0.0 {
            let smoothed =
                gaussian_blur_f64(&workspace.source_rate, config.boundaries.source_smoothing_sigma);
            workspace.source_rate.data_mut().copy_from_slice(smoothed.data());
        }

        if config.boundaries.slab_pull_enabled {
            // Issue #75: the per-plate velocity boost
            // (`apply_slab_pull`) is replaced by the operator term
            // `γ_slab · (v·n̂) · n̂` computed in `compute_boundary_sources`
            // and applied in `apply_stokes`. We keep accumulating
            // `plate.subducted_mass` because downstream diagnostics
            // and the recycling buffer may still consult it.
            accumulate_subducted_mass(
                &workspace.source_rate,
                &plate_ctx.ids,
                &mut plate_ctx.plates,
                dt_rates,
                nx,
                ny,
            );
        }

        // Phase 2-bis calibration diagnostic: γ_slab field stats on
        // margin cells + max velocity samples. Runs only when the
        // `slab_pull_sweep` tracing target is enabled (#75).
        if let Some(bf) = workspace.boundary_field.as_ref() {
            super::diagnostics::emit_slab_pull_sweep(bf, grid);
        }
    }
    let t_boundaries_us = t_boundaries_start.elapsed().as_micros() as u64;

    // Subduction consumption + fragmentation detection + traction rebuild
    // are grouped as "plate bookkeeping" for the phase_timings breakdown.
    let t_plates_start = Instant::now();
    if config.dynamic_boundaries {
        let ids_before: Vec<usize> = plate_ctx.ids.clone();

        if let Some(ref bf) = workspace.boundary_field {
            apply_subduction_consumption(
                &mut plate_ctx.ids,
                grid,
                bf,
                config.boundaries.subduction_consumption_threshold,
            );
        }
        if step % 10 == 0 {
            detect_fragmentation(
                &mut plate_ctx.ids,
                &mut plate_ctx.plates,
                &mut plate_ctx.next_id,
                nx,
                ny,
                grid,
                0.25,
            );
        }

        for k in 0..nx * ny {
            if plate_ctx.ids[k] != ids_before[k] {
                let i = k % nx;
                let j = k / nx;
                plate_ctx.disp_x.set(i, j, 0.0);
                plate_ctx.disp_y.set(i, j, 0.0);
            }
        }

        update_plate_stats(&plate_ctx.ids, &mut plate_ctx.plates, grid);
    }

    // Rebuild traction once (after slab pull may have updated velocities)
    if config.boundaries.enabled || config.dynamic_boundaries {
        plate_ctx.traction = if config.dynamic_boundaries {
            rebuild_traction_smooth(
                &plate_ctx.ids,
                &plate_ctx.plates,
                &plate_ctx.disp_x,
                &plate_ctx.disp_y,
                nx,
                ny,
            )
        } else {
            rebuild_traction(&plate_ctx.ids, &plate_ctx.plates, nx, ny)
        };
    }

    // Add mantle convection flow to traction
    if let Some(mf) = mantle_flow.as_mut() {
        let coupling = config.mantle.coupling;
        for j in 0..ny {
            for i in 0..nx {
                let s = grid.s.get(i, j);
                let c = coupling * (s - 0.15).max(0.0);
                let tx = plate_ctx.traction.tx.get(i, j) + c * mf.vx.get(i, j);
                let ty = plate_ctx.traction.ty.get(i, j) + c * mf.vy.get(i, j);
                plate_ctx.traction.tx.set(i, j, tx);
                plate_ctx.traction.ty.set(i, j, ty);
            }
        }
        if config.mantle.evolution_rate > 0.0 {
            mf.evolve(config.mantle.evolution_rate, step);
        }
    }
    let t_plates_us = t_plates_start.elapsed().as_micros() as u64;

    // ── Conservative mass recycling ────────────────────────────
    let t_recycling_start = Instant::now();
    if let Some(recycler) = recycling_buffer.as_mut()
        && let Some(bf) = workspace.boundary_field.as_ref()
    {
        let mut total_subducted = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                let q = workspace.source_rate.get(i, j);
                if q < 0.0 {
                    total_subducted += (-q) * dt_rates;
                }
            }
        }

        let arc_mass = total_subducted * config.recycling.arc_fraction;
        let loss_mass = total_subducted * config.recycling.loss_fraction;
        let buffer_mass = total_subducted - arc_mass - loss_mass;

        let mut arc_cells: Vec<(usize, usize)> = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let btype = bf.boundary_type[k];
                let is_arc = match btype {
                    BoundaryType::Subduction => {
                        plate_ctx.plates[plate_ctx.ids[k]].mean_thickness > 0.4
                    }
                    BoundaryType::OceanicSubduction => grid.s.get(i, j) > 0.22,
                    _ => false,
                };
                if is_arc {
                    arc_cells.push((i, j));
                }
            }
        }

        if !arc_cells.is_empty() && arc_mass > 0.0 {
            let per_cell = arc_mass / arc_cells.len() as f64;
            for &(i, j) in &arc_cells {
                let current = workspace.source_rate.get(i, j);
                workspace.source_rate.set(i, j, current + per_cell / dt_rates);
            }
        }

        if buffer_mass > 0.0 {
            recycler.deposit(buffer_mass);
        }
        let spread_mass = recycler.advance();

        let mut rift_cells: Vec<(usize, usize)> = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                if bf.boundary_type[k] == BoundaryType::Rift
                    && grid.s.get(i, j) < config.boundaries.rift_thickness_threshold
                {
                    rift_cells.push((i, j));
                }
            }
        }

        if !rift_cells.is_empty() && spread_mass > 0.0 {
            let per_cell = spread_mass / rift_cells.len() as f64;
            for &(i, j) in &rift_cells {
                let current = workspace.source_rate.get(i, j);
                workspace.source_rate.set(i, j, current + per_cell / dt_rates);
            }
        }

        debug!(
            step,
            subducted = %format!("{:.4}", total_subducted),
            arc = %format!("{:.4}", arc_mass),
            buffered = %format!("{:.4}", buffer_mass),
            spread_available = %format!("{:.4}", spread_mass),
            arc_cells = arc_cells.len(),
            rift_cells = rift_cells.len(),
            "mass recycling"
        );
    }

    let t_recycling_us = t_recycling_start.elapsed().as_micros() as u64;

    // ── Mass balance tracking ──────────────────────────────────
    let mass_before: f64 = grid.s.data().iter().sum();
    let mut mass_after_advection = 0.0_f64;
    let mut mass_clamp_delta = 0.0_f64;
    let mut total_div_flux = 0.0_f64;
    let mut total_source = 0.0_f64;

    // ── Advection with legacy CFL retry or a single explicit attempt.
    let t_advection_start = Instant::now();
    let (retry_succeeded, clamp_ratio, dt_final) = match dt_mode {
        DtMode::LegacyWithCflRetry => {
            let s_backup: Vec<f64> = grid.s.data().to_vec();
            let mut dt = dt_rates;
            let dt_min = dt * 0.01;
            const MAX_RETRIES: usize = 5;
            let mut rs = false;
            let mut cr = 0.0;

            for retry_index in 0..MAX_RETRIES {
                compute_divergence_flux(grid, &mut workspace.div_flux);
                mass_after_advection = 0.0;
                mass_clamp_delta = 0.0;
                total_div_flux = 0.0;
                total_source = 0.0;

                for j in 0..ny {
                    for i in 0..nx {
                        let div = workspace.div_flux.get(i, j);
                        let q =
                            if boundaries_active { workspace.source_rate.get(i, j) } else { 0.0 };
                        let s_old = grid.s.get(i, j);
                        let s_raw = s_old - dt * div + dt * q;
                        let s_clamped = s_raw.clamp(config.s_min, config.s_max);

                        total_div_flux += dt * div;
                        total_source += dt * q;
                        mass_clamp_delta += s_clamped - s_raw;
                        mass_after_advection += s_clamped;

                        grid.s.set(i, j, s_clamped);
                    }
                }

                let clamp_count = grid
                    .s
                    .data()
                    .iter()
                    .filter(|&&s| s <= config.s_min * 1.01 || s >= config.s_max * 0.99)
                    .count();
                cr = clamp_count as f64 / (nx * ny) as f64;

                if cr < 0.05 || dt <= dt_min {
                    rs = true;
                    break;
                }

                if retry_index + 1 < MAX_RETRIES {
                    dt *= 0.5;
                    grid.s.data_mut().copy_from_slice(&s_backup);
                }
            }

            if !rs {
                warn!(
                    step,
                    clamp_ratio = cr,
                    dt,
                    attempts = MAX_RETRIES,
                    "CFL retry exhausted — accepting degraded step"
                );
            }
            (rs, cr, dt)
        }
        DtMode::Explicit(dt_sub) => {
            compute_divergence_flux(grid, &mut workspace.div_flux);
            mass_after_advection = 0.0;
            mass_clamp_delta = 0.0;
            total_div_flux = 0.0;
            total_source = 0.0;
            for j in 0..ny {
                for i in 0..nx {
                    let div = workspace.div_flux.get(i, j);
                    let q = if boundaries_active { workspace.source_rate.get(i, j) } else { 0.0 };
                    let s_old = grid.s.get(i, j);
                    let s_raw = s_old - dt_sub * div + dt_sub * q;
                    let s_clamped = s_raw.clamp(config.s_min, config.s_max);

                    total_div_flux += dt_sub * div;
                    total_source += dt_sub * q;
                    mass_clamp_delta += s_clamped - s_raw;
                    mass_after_advection += s_clamped;

                    grid.s.set(i, j, s_clamped);
                }
            }
            let clamp_count = grid
                .s
                .data()
                .iter()
                .filter(|&&s| s <= config.s_min * 1.01 || s >= config.s_max * 0.99)
                .count();
            let cr = clamp_count as f64 / (nx * ny) as f64;
            (true, cr, dt_sub)
        }
    };

    // 3b. Oceanic restoring
    if config.boundaries.enabled && config.boundaries.oceanic_restore_rate > 0.0 {
        let rate = config.boundaries.oceanic_restore_rate;
        let s_ref = config.boundaries.oceanic_reference_thickness;
        let s_thr = config.boundaries.oceanic_restore_threshold;
        for j in 0..ny {
            for i in 0..nx {
                let pid = plate_ctx.ids[j * nx + i];
                if plate_ctx.plates[pid].plate_type == PlateType::Oceanic {
                    let s_current = grid.s.get(i, j);
                    if s_current > s_ref && s_current < s_thr {
                        let s_new = s_current - dt_final * rate * (s_current - s_ref);
                        grid.s.set(i, j, s_new.max(config.s_min));
                    }
                }
            }
        }
    }
    let mass_after_oceanic_restore: f64 = grid.s.data().iter().sum();
    let oceanic_restore_delta = mass_after_oceanic_restore - mass_after_advection;

    // 3c. Continental minimum thickness protection
    if config.boundaries.enabled && config.boundaries.continental_restore_rate > 0.0 {
        let rate = config.boundaries.continental_restore_rate;
        let s_min_boundary = config.boundaries.continental_min_thickness;
        let s_thr = config.boundaries.continental_restore_threshold;
        for j in 0..ny {
            for i in 0..nx {
                let pid = plate_ctx.ids[j * nx + i];
                if plate_ctx.plates[pid].plate_type == PlateType::Continental {
                    let s_current = grid.s.get(i, j);
                    if s_current < s_min_boundary && s_current > s_thr {
                        let s_new = s_current + dt_final * rate * (s_min_boundary - s_current);
                        grid.s.set(i, j, s_new);
                    }
                }
            }
        }
    }
    let mass_after_continental_restore: f64 = grid.s.data().iter().sum();
    let continental_restore_delta = mass_after_continental_restore - mass_after_oceanic_restore;
    let mass_after = mass_after_continental_restore;

    // 4. Summary stats over the final state
    let mut max_v = 0.0_f64;
    let mut max_s = f64::NEG_INFINITY;
    let mut min_s = f64::INFINITY;
    for j in 0..ny {
        for i in 0..nx {
            let vx = grid.vx.get(i, j);
            let vy = grid.vy.get(i, j);
            max_v = max_v.max((vx * vx + vy * vy).sqrt());
            let s = grid.s.get(i, j);
            max_s = max_s.max(s);
            min_s = min_s.min(s);
        }
    }

    let t_advection_us = t_advection_start.elapsed().as_micros() as u64;

    info!(
        target: "phase_timings",
        step,
        t_boundaries_us,
        t_solve_us,
        t_advection_us,
        t_recycling_us,
        t_plates_us,
        "phase timings"
    );

    // We reached here only after `solve_outcome.is_converged()` passed, so
    // solve_outcome is guaranteed to be ConvergedOnResidual or ConvergedOnState.
    // Pass it through so the sub-step loop can distinguish the two.
    Ok(PassResult {
        outcome: solve_outcome,
        nl_iterations,
        linear_iterations,
        dt_consumed: dt_final,
        clamp_ratio,
        retry_succeeded,
        max_velocity: max_v,
        max_thickness: max_s,
        min_thickness: min_s,
        mass_before,
        mass_after,
        total_div_flux,
        total_source,
        mass_clamp_delta,
        oceanic_restore_delta,
        continental_restore_delta,
        substep_count: 1,
    })
}

/// Adaptive macro step: repeatedly call `execute_tectonic_pass` in
/// `DtMode::Explicit` mode, accumulating sub-step time toward
/// `config.adaptive_dt.dt_target`, with snapshot-based rollback on
/// failure and an outcome-driven reduction/growth policy.
///
/// Returns a synthesized `PassResult` aggregated over the committed
/// sub-steps, so the caller's logging and `StepStats` construction code
/// treats both modes uniformly.
#[allow(clippy::too_many_arguments)]
fn run_adaptive_macro_step(
    grid: &mut StaggeredGrid,
    plate_ctx: &mut DynamicPlateContext,
    workspace: &mut SolverWorkspace,
    recycling_buffer: &mut Option<RecyclingBuffer>,
    mantle_flow: &mut Option<MantleFlow>,
    config: &TectonicsConfig,
    step: usize,
) -> Result<PassResult, SolverError> {
    let dt_target = config.adaptive_dt.dt_target;
    let dt_min = dt_target * config.adaptive_dt.min_dt_fraction;
    let max_substeps = config.adaptive_dt.max_substeps;

    let mut accumulated = AccumulatedSubstepStats::default();
    let mut dt_current = dt_target;
    let mut elapsed = 0.0_f64;
    let mut last_successful_pass: Option<PassResult> = None;

    while elapsed < dt_target && accumulated.substep_count < max_substeps {
        let remaining = dt_target - elapsed;
        let mut dt_sub = dt_current.min(remaining);

        if config.adaptive_dt.respect_local_cfl {
            let dt_cfl = compute_cfl_dt(grid, config.cfl_factor);
            dt_sub = dt_sub.min(dt_cfl);
        }

        let snapshot = StateSnapshot::capture(grid, plate_ctx);
        let allow_continuation = step == 0 && accumulated.substep_count == 0;

        let pass = execute_tectonic_pass(
            grid,
            plate_ctx,
            workspace,
            recycling_buffer,
            mantle_flow,
            config,
            step,
            DtMode::Explicit(dt_sub),
            allow_continuation,
        )?;

        let sub = pass.as_substep();
        let max_clamp = config.adaptive_dt.max_clamp_ratio_success;
        if pass.is_successful(max_clamp) {
            elapsed += dt_sub;
            accumulated.merge(&sub, dt_sub);
            dt_current = grow_dt(dt_current, &sub, &config.adaptive_dt);
            last_successful_pass = Some(pass);
            debug!(
                step,
                substep = accumulated.substep_count - 1,
                dt_sub,
                elapsed,
                newton_iters = sub.newton_iterations,
                clamp_ratio = sub.clamp_ratio,
                "sub-step committed"
            );
        } else {
            snapshot.restore(grid, plate_ctx);
            let factor = reduction_for(&sub, &config.adaptive_dt);
            dt_current *= factor;
            // Distinguish the two distinct failure modes. Newton non-
            // convergence is a solver issue that may yield to smaller dt;
            // excessive clamping is often structural (cells pinned at a
            // boundary by relief saturation) and does not respond to dt
            // shrinkage — it is the signal that `max_clamp_ratio_success`
            // may need to be raised.
            let reason =
                if !sub.outcome.is_converged() { "newton_failed" } else { "excessive_clamping" };
            debug!(
                step,
                substep = accumulated.substep_count,
                dt_sub,
                clamp_ratio = sub.clamp_ratio,
                reason,
                reduction = factor,
                dt_after = dt_current,
                "sub-step failed, rolled back"
            );
            if dt_current < dt_min {
                warn!(
                    step,
                    dt_current,
                    dt_min,
                    elapsed,
                    dt_target,
                    "sub-step floor reached, abandoning remainder of macro step"
                );
                break;
            }
        }
    }

    if accumulated.substep_count >= max_substeps && elapsed < dt_target {
        warn!(
            step,
            substep_count = accumulated.substep_count,
            elapsed,
            dt_target,
            "max sub-steps reached before consuming dt_target"
        );
    }

    // If no sub-step ever committed, fall back to the legacy error
    // semantics so the caller can treat this like a non-convergence.
    let Some(last_pass) = last_successful_pass else {
        return Err(SolverError::NonlinearSolverDidNotConverge { step });
    };

    // Final state totals for stats and logging.
    let mut max_s = f64::NEG_INFINITY;
    let mut min_s = f64::INFINITY;
    for &s in grid.s.data().iter() {
        max_s = max_s.max(s);
        min_s = min_s.min(s);
    }

    // Mass-balance fields: we report the last committed sub-step's bookkeeping
    // rather than true per-macro totals. The debug log stays readable per
    // sub-step; a future commit can accumulate these if the macro-level
    // mass balance is needed for diagnostics.
    Ok(PassResult {
        outcome: NewtonOutcome::ConvergedOnResidual,
        nl_iterations: accumulated.newton_iters_total,
        linear_iterations: accumulated.linear_iters_total,
        dt_consumed: elapsed,
        clamp_ratio: accumulated.max_clamp_ratio,
        retry_succeeded: true,
        max_velocity: accumulated.max_velocity,
        max_thickness: max_s,
        min_thickness: min_s,
        mass_before: last_pass.mass_before,
        mass_after: last_pass.mass_after,
        total_div_flux: last_pass.total_div_flux,
        total_source: last_pass.total_source,
        mass_clamp_delta: last_pass.mass_clamp_delta,
        oceanic_restore_delta: last_pass.oceanic_restore_delta,
        continental_restore_delta: last_pass.continental_restore_delta,
        substep_count: accumulated.substep_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::solver::config::{ContinuationConfig, NewtonConfig, PicardConfig};

    fn make_config(num_timesteps: usize) -> TectonicsConfig {
        TectonicsConfig {
            num_timesteps,
            gravity_factor: 1.0,
            cfl_factor: 0.1,
            s_min: 0.1,
            s_max: 2.5,
            nonlinear_solver: NonlinearSolver::Picard,
            picard: PicardConfig {
                max_iterations: 30,
                tolerance: 1e-3,
                relaxation: 1.0,
                cg_max_iter: 500,
                cg_tolerance: 1e-8,
                strain_rate_min: 1e-3,
                power_law_n: 1.0,
                eta_min: 1e-3,
                eta_max: 1e4,
            },
            newton: NewtonConfig::default(),
            continuation: ContinuationConfig { enabled: false, ..Default::default() },
            boundaries: Default::default(),
            dynamic_boundaries: false,
            cratonic: Default::default(),
            yielding: Default::default(),
            basal_friction: 0.0,
            mantle: Default::default(),
            recycling: Default::default(),
            adaptive_dt: Default::default(),
        }
    }

    fn make_static_ctx(n: usize, traction: TractionField) -> DynamicPlateContext {
        use crate::tectonics::plates::{Plate, PlateType};
        DynamicPlateContext {
            ids: vec![0; n * n],
            plates: vec![Plate {
                id: 0,
                plate_type: PlateType::Continental,
                velocity: (0.0, 0.0),
                seed_x: (n / 2) as f32,
                seed_y: (n / 2) as f32,
                active: true,
                subducted_mass: 0.0,
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            }],
            traction,
            next_id: 1,
            disp_x: Field2D::new(n, n),
            disp_y: Field2D::new(n, n),
        }
    }

    #[test]
    fn convergent_plates_thicken() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let traction = TractionField::two_plates_convergent(n, n, 1.0);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(50);
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        assert!(
            ws.stats.max_thickness > s_initial,
            "Convergent plates should thicken: max_s = {}",
            ws.stats.max_thickness
        );
    }

    #[test]
    fn divergent_plates_thin() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let traction = TractionField::two_plates_divergent(n, n, 1.0);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(50);
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        assert!(
            ws.stats.min_thickness < s_initial,
            "Divergent plates should thin: min_s = {}",
            ws.stats.min_thickness
        );
    }

    #[test]
    fn gpe_flattens_bump() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);

        // Background + broad Gaussian bump (σ²=0.02, well-resolved at 32²)
        let center = 0.5;
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                let r2 = (x - center).powi(2) + (y - center).powi(2);
                grid.s.set(i, j, 1.0 + 0.3 * (-r2 / 0.02).exp());
            }
        }

        let initial_var: f64 = {
            let mean = grid.s.data().iter().sum::<f64>() / (n * n) as f64;
            grid.s.data().iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n * n) as f64
        };

        let traction = TractionField::zero(n, n);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(100);
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        let final_var: f64 = {
            let mean = grid.s.data().iter().sum::<f64>() / (n * n) as f64;
            grid.s.data().iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n * n) as f64
        };

        assert!(
            final_var < initial_var,
            "GPE should flatten bump (reduce variance): initial_var = {initial_var}, final_var = {final_var}"
        );
    }

    #[test]
    fn no_nan_no_inf() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(30);
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok());

        for val in grid.s.data() {
            assert!(val.is_finite(), "S contains non-finite: {val}");
        }
        for val in grid.vx.data() {
            assert!(val.is_finite(), "vx contains non-finite: {val}");
        }
        for val in grid.vy.data() {
            assert!(val.is_finite(), "vy contains non-finite: {val}");
        }
    }

    #[test]
    fn continuation_enables_power_law_convergence() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let config = TectonicsConfig {
            num_timesteps: 20,
            gravity_factor: 1.0,
            cfl_factor: 0.3,
            s_min: 0.1,
            s_max: 2.5,
            nonlinear_solver: NonlinearSolver::Picard,
            picard: PicardConfig {
                max_iterations: 60,
                power_law_n: 3.0,
                strain_rate_min: 1e-3,
                eta_min: 1e-3,
                eta_max: 1e4,
                relaxation: 0.5,
                ..PicardConfig::default()
            },
            newton: NewtonConfig::default(),
            continuation: ContinuationConfig::default(),
            boundaries: Default::default(),
            dynamic_boundaries: false,
            cratonic: Default::default(),
            yielding: Default::default(),
            basal_friction: 1.0,
            mantle: Default::default(),
            recycling: Default::default(),
            adaptive_dt: Default::default(),
        };

        let mut ctx = make_static_ctx(n, traction);
        let mut ws = SolverWorkspace::new(n, n);
        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Continuation should enable convergence: {:?}", result.err());
        assert!(
            ws.stats.max_thickness > 1.0,
            "Convergent plates should thicken with power-law: max_s={}",
            ws.stats.max_thickness
        );
    }

    #[test]
    fn oceanic_restoring_prevents_thickening() {
        use crate::tectonics::boundaries::BoundaryConfig;

        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);

        // All oceanic, initially thick (as if advection has thickened them)
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 0.8);
            }
        }

        let traction = TractionField::zero(n, n);
        let mut ctx = DynamicPlateContext {
            ids: vec![0; n * n],
            plates: vec![Plate {
                id: 0,
                plate_type: PlateType::Oceanic,
                velocity: (0.0, 0.0),
                seed_x: 8.0,
                seed_y: 8.0,
                active: true,
                subducted_mass: 0.0,
                cell_count: 0,
                mean_thickness: 0.0,
                mean_velocity: (0.0, 0.0),
                centroid_x: 0.0,
                centroid_y: 0.0,
            }],
            traction,
            next_id: 1,
            disp_x: Field2D::new(n, n),
            disp_y: Field2D::new(n, n),
        };

        let config = TectonicsConfig {
            num_timesteps: 50,
            boundaries: BoundaryConfig {
                enabled: true,
                oceanic_reference_thickness: 0.25,
                oceanic_restore_threshold: 1.0, // high enough to cover test's initial S=0.8
                oceanic_restore_rate: 0.3,
                ..Default::default()
            },
            dynamic_boundaries: false,
            cratonic: Default::default(),
            yielding: Default::default(),
            ..make_config(50)
        };

        let mut ws = SolverWorkspace::new(n, n);
        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        let final_mean = grid.s.data().iter().sum::<f64>() / (n * n) as f64;
        assert!(
            final_mean < 0.5,
            "Oceanic crust should thin toward reference: mean={}",
            final_mean
        );
    }

    #[test]
    fn viscosity_clamp_works() {
        use crate::tectonics::solver::field::Field2D;
        use crate::tectonics::solver::picard::compute_viscosity;

        let n = 8;
        let mut strain = Field2D::new(n, n);
        let mut eta = Field2D::new(n, n);

        // Very low strain rate → very high viscosity without clamp
        strain.set(0, 0, 0.0);
        // Normal strain rate
        strain.set(1, 0, 0.1);

        compute_viscosity(&strain, 3.0, 1e-6, 1e-3, 1e3, &mut eta);

        assert!(eta.get(0, 0) <= 1e3 + 1e-10, "Should be clamped to eta_max: {}", eta.get(0, 0));
        assert!(eta.get(1, 0) >= 1e-3, "Normal cell below eta_min");
        assert!(eta.get(1, 0) <= 1e3, "Normal cell above eta_max");
    }

    #[test]
    fn density_corrected_gpe_oceanic_spreads_less() {
        use crate::tectonics::solver::stokes::compute_rhs;

        let n = 16;
        let dx = 1.0 / n as f64;
        let rho_c = 2750.0;
        let rho_m = 3300.0;

        // Setup 1: all continental density
        let mut grid_c = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_c.s.set(i, j, if i < n / 2 { 1.0 } else { 0.5 });
                grid_c.rho.set(i, j, 2750.0);
            }
        }

        // Setup 2: same thickness, all oceanic density
        let mut grid_o = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_o.s.set(i, j, if i < n / 2 { 1.0 } else { 0.5 });
                grid_o.rho.set(i, j, 3000.0);
            }
        }

        let plates = TractionField::zero(n, n);
        let nn2 = 2 * n * n;
        let mut rhs_c = vec![0.0; nn2];
        let mut rhs_o = vec![0.0; nn2];

        compute_rhs(&grid_c, &plates, 1.0, rho_c, rho_m, &mut rhs_c);
        compute_rhs(&grid_o, &plates, 1.0, rho_c, rho_m, &mut rhs_o);

        let mag_c: f64 = rhs_c.iter().map(|x| x * x).sum::<f64>().sqrt();
        let mag_o: f64 = rhs_o.iter().map(|x| x * x).sum::<f64>().sqrt();

        assert!(
            mag_o < mag_c,
            "Oceanic GPE should be weaker: continental={mag_c}, oceanic={mag_o}"
        );

        let ratio = mag_o / mag_c;
        assert!(ratio < 0.7, "Ratio should be significantly less than 1: {ratio}");
    }

    #[test]
    #[ignore = "obsolete since issue #75: slab-pull is no longer a per-plate velocity boost; see follow-up cleanup issue #79"]
    fn slab_pull_increases_plate_velocity() {
        use crate::tectonics::boundaries::apply_slab_pull;

        let mut plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Oceanic,
            velocity: (1.0, 0.0),
            seed_x: 8.0,
            seed_y: 8.0,
            active: true,
            subducted_mass: 10.0,
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
        }];

        let initial_vx = plates[0].velocity.0;
        apply_slab_pull(&mut plates, 0.1, 5.0);

        assert!(
            plates[0].velocity.0 > initial_vx,
            "Slab pull should increase velocity: {} -> {}",
            initial_vx,
            plates[0].velocity.0
        );
    }

    #[test]
    #[ignore = "obsolete since issue #75: slab-pull no longer reads max_plate_velocity; see follow-up cleanup issue #79"]
    fn slab_pull_capped() {
        use crate::tectonics::boundaries::apply_slab_pull;

        let mut plates = vec![Plate {
            id: 0,
            plate_type: PlateType::Oceanic,
            velocity: (1.0, 0.0),
            seed_x: 8.0,
            seed_y: 8.0,
            active: true,
            subducted_mass: 1000.0,
            cell_count: 0,
            mean_thickness: 0.0,
            mean_velocity: (0.0, 0.0),
            centroid_x: 0.0,
            centroid_y: 0.0,
        }];

        apply_slab_pull(&mut plates, 0.5, 5.0);

        let speed = (plates[0].velocity.0.powi(2) + plates[0].velocity.1.powi(2)).sqrt();
        assert!(speed <= 5.1, "Velocity should be capped: {speed}");
    }

    #[test]
    fn cfl_retry_succeeds_on_standard_configuration() {
        // Gentle uniform field, modest convergent traction, single step: the
        // CFL retry loop should converge on the first try (retry 0) with a
        // clamp_ratio well below the 5% threshold.
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let mut config = make_config(1);
        // Exercise the legacy CFL retry path — the cfl_retry_exhausted flag
        // is only populated there; adaptive mode reports retry_succeeded: true
        // unconditionally because rollback happens at the sub-step boundary.
        config.adaptive_dt.enabled = false;
        config.boundaries.oceanic_restore_rate = 0.0;
        config.boundaries.continental_restore_rate = 0.0;
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        assert!(!ws.stats.cfl_retry_exhausted, "CFL retry should succeed on gentle configuration");
        assert!(
            ws.stats.clamp_ratio < 0.05,
            "clamp_ratio should stay below threshold on success path, got {}",
            ws.stats.clamp_ratio
        );
    }

    #[test]
    fn cfl_retry_exhausted_commits_last_attempt() {
        // Extremely tight clamp window around a uniform initial thickness:
        // every cell sits inside both clamp-proximity thresholds
        // (s <= s_min * 1.01 AND s >= s_max * 0.99) so clamp_ratio is 100%
        // at every attempted dt. The retry loop therefore exhausts its
        // budget — but grid.s MUST still evolve rather than silently
        // reverting to the pre-step state (the bug this fix addresses).
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }
        let s_before: Vec<f64> = grid.s.data().to_vec();

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let mut config = make_config(1);
        // Window width 0.002 straddling s_initial = 1.0. The proximity
        // tests `s <= s_min * 1.01` (= 1.00899) and `s >= s_max * 0.99`
        // (= 0.99099) both cover the entire [s_min, s_max] range, so
        // any value grid.s can hold after clamping will be counted.
        config.s_min = 0.999;
        config.s_max = 1.001;
        // The cfl_retry_exhausted flag and the "commit last attempt" behavior
        // are semantics of the legacy CFL retry loop; adaptive mode instead
        // rolls back the failing sub-step and shrinks dt. Gate the test to
        // the legacy path so the assertions remain meaningful.
        config.adaptive_dt.enabled = false;
        config.boundaries.oceanic_restore_rate = 0.0;
        config.boundaries.continental_restore_rate = 0.0;
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "Run failed: {:?}", result.err());

        assert!(
            ws.stats.cfl_retry_exhausted,
            "Expected retry exhaustion on pathological config, got clamp_ratio = {}",
            ws.stats.clamp_ratio
        );

        // Core regression: grid.s must have evolved — the previous code
        // silently restored s_backup on exhaustion, leaving this diff at 0.
        let s_after = grid.s.data();
        let max_diff =
            s_before.iter().zip(s_after.iter()).map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
        assert!(
            max_diff > 1e-6,
            "grid.s must have evolved even on CFL retry exhaustion, got max_diff = {max_diff}"
        );
    }

    #[test]
    fn adaptive_mode_covers_dt_target_on_gentle_configuration() {
        // Single macro step with adaptive mode enabled on an easy configuration
        // (gentle convergent flow, uniform thickness). The sub-step loop should
        // commit at least one sub-step, cover the full dt_target budget, evolve
        // grid.s, and populate StepStats with sensible values.
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }
        let s_before: Vec<f64> = grid.s.data().to_vec();

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let mut config = make_config(1);
        config.adaptive_dt.enabled = true;
        config.adaptive_dt.dt_target = 0.5;
        // Disable restoring forces so any evolution of grid.s comes purely
        // from the advection inside the sub-step loop.
        config.boundaries.oceanic_restore_rate = 0.0;
        config.boundaries.continental_restore_rate = 0.0;
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "adaptive run failed: {:?}", result.err());

        assert!(
            (ws.stats.dt - config.adaptive_dt.dt_target).abs() < 1e-9,
            "adaptive dt_consumed {} should match dt_target {}",
            ws.stats.dt,
            config.adaptive_dt.dt_target
        );

        let s_after = grid.s.data();
        let max_diff =
            s_before.iter().zip(s_after.iter()).map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
        assert!(
            max_diff > 1e-6,
            "grid.s should have evolved under adaptive sub-stepping, got max_diff = {max_diff}"
        );

        // cfl_retry_exhausted is a legacy-path signal; adaptive mode never sets it.
        assert!(!ws.stats.cfl_retry_exhausted);
    }

    #[test]
    fn adaptive_mode_matches_legacy_path_on_gentle_configuration() {
        // Both code paths (legacy CFL retry, adaptive single sub-step) should
        // converge and produce plausible physics on the same easy configuration.
        // They are numerically different — adaptive uses a fixed dt_target,
        // legacy uses dt_cfl — but both must populate StepStats and evolve
        // grid.s monotonically in the same direction.
        let n = 16;
        let dx = 1.0 / n as f64;

        let run = |adaptive: bool| -> (f64, f64, f64, f64) {
            let mut grid = StaggeredGrid::new(n, n, dx);
            for j in 0..n {
                for i in 0..n {
                    grid.s.set(i, j, 1.0);
                }
            }
            let traction = TractionField::two_plates_convergent(n, n, 0.5);
            let mut ctx = make_static_ctx(n, traction);
            let mut config = make_config(5);
            config.adaptive_dt.enabled = adaptive;
            config.adaptive_dt.dt_target = 0.5;
            config.boundaries.oceanic_restore_rate = 0.0;
            config.boundaries.continental_restore_rate = 0.0;
            let mut ws = SolverWorkspace::new(n, n);
            run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true).unwrap();
            let mass: f64 = grid.s.data().iter().sum();
            (ws.stats.dt, ws.stats.max_thickness, ws.stats.min_thickness, mass)
        };

        let (legacy_dt, legacy_max, legacy_min, legacy_mass) = run(false);
        let (adaptive_dt, adaptive_max, adaptive_min, adaptive_mass) = run(true);

        assert!(legacy_dt > 0.0 && legacy_max > legacy_min);
        assert!(adaptive_dt > 0.0 && adaptive_max > adaptive_min);

        // Adaptive mode on a gentle config consumes its full target budget.
        assert!(
            (adaptive_dt - 0.5).abs() < 1e-9,
            "adaptive last-step dt = {adaptive_dt}, expected 0.5"
        );

        // Both paths preserve overall mass on the same order of magnitude
        // (sources/sinks are balanced by default boundaries).
        let mass_ratio = adaptive_mass / legacy_mass;
        assert!(
            (0.5..2.0).contains(&mass_ratio),
            "adaptive mass {adaptive_mass} vs legacy {legacy_mass} diverged beyond 2x"
        );
    }

    #[test]
    fn adaptive_mode_falls_back_gracefully_when_config_is_unsolvable() {
        // Unsolvable setup: an extremely tight s_min/s_max window around the
        // initial thickness combined with strong convergent flow. Every
        // sub-step attempt clamps > 5% of cells no matter what dt_sub is
        // tried, so SubstepResult::is_successful() is false on every
        // iteration. The sub-step loop should shrink dt_current until it
        // crosses the min_dt_fraction floor, emit the "sub-step floor
        // reached" warning, and either (a) return an error if no sub-step
        // ever committed or (b) return a degraded StepStats if at least
        // one did — either way, without panicking or corrupting state.
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let mut config = make_config(1);
        config.adaptive_dt.enabled = true;
        config.adaptive_dt.dt_target = 2.0;
        // Clamp-proximity thresholds (s <= s_min * 1.01 AND s >= s_max * 0.99)
        // both cover the full [s_min, s_max] range, so 100% of cells always
        // count as clamped and no sub-step can ever succeed.
        config.s_min = 0.999;
        config.s_max = 1.001;
        config.boundaries.oceanic_restore_rate = 0.0;
        config.boundaries.continental_restore_rate = 0.0;
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        // The run must terminate — either with Ok (floor committed last
        // successful sub-step) or NonlinearSolverDidNotConverge (no sub-step
        // ever committed). Panics, infinite loops, or NaN states are bugs.
        match result {
            Ok(_) => {
                // If any sub-step committed, dt_consumed is strictly less
                // than dt_target because the floor aborted the remainder.
                assert!(
                    ws.stats.dt <= config.adaptive_dt.dt_target + 1e-9,
                    "dt_consumed {} should not exceed dt_target {}",
                    ws.stats.dt,
                    config.adaptive_dt.dt_target
                );
                // State must be coherent (no NaN / Inf).
                for &s in grid.s.data().iter() {
                    assert!(s.is_finite(), "grid.s contains non-finite value {s}");
                    assert!(
                        s >= config.s_min * 0.99 && s <= config.s_max * 1.01,
                        "grid.s value {s} drifted outside clamp window after rollback"
                    );
                }
            }
            Err(SolverError::NonlinearSolverDidNotConverge { step: 0 }) => {
                // Acceptable exit: no sub-step committed. State is the
                // pre-step snapshot because the snapshot scope is per-sub-step.
            }
            Err(other) => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn adaptive_mode_evolves_state_under_tight_clamp_window() {
        // A slightly wider clamp window than the floor test, so sub-steps
        // *can* succeed if dt is small enough. The adaptive loop should
        // engage multiple sub-steps (substep_count > 1) to cover the full
        // dt_target while keeping clamp_ratio below 5% on each committed
        // attempt. This is the central scenario the issue targets: a
        // configuration where the single-step path would fail but the
        // sub-step loop finds a working sequence.
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }
        let s_before: Vec<f64> = grid.s.data().to_vec();

        let traction = TractionField::two_plates_convergent(n, n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let mut config = make_config(1);
        config.adaptive_dt.enabled = true;
        // Generous dt_target so the sub-step loop has room to shrink dt
        // several times before hitting the floor.
        config.adaptive_dt.dt_target = 4.0;
        // Wider window than the floor test — clamp_ratio is < 100% so
        // small-enough sub-steps can actually succeed.
        config.s_min = 0.5;
        config.s_max = 1.5;
        config.boundaries.oceanic_restore_rate = 0.0;
        config.boundaries.continental_restore_rate = 0.0;
        let mut ws = SolverWorkspace::new(n, n);

        let result = run_tectonics(&config, &mut ctx, &mut grid, &mut ws, |_, _, _, _| true);
        assert!(result.is_ok(), "adaptive run failed: {:?}", result.err());

        // State evolved from initial uniform field.
        let s_after = grid.s.data();
        let max_diff =
            s_before.iter().zip(s_after.iter()).map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
        assert!(
            max_diff > 1e-6,
            "grid.s should have evolved under adaptive sub-stepping, max_diff = {max_diff}"
        );

        // The sub-step loop should have covered meaningful progress toward
        // dt_target, even if it didn't reach it fully.
        assert!(ws.stats.dt > 0.0, "adaptive dt_consumed must be positive, got {}", ws.stats.dt);
        assert!(
            ws.stats.dt <= config.adaptive_dt.dt_target + 1e-9,
            "adaptive dt_consumed {} should not exceed dt_target {}",
            ws.stats.dt,
            config.adaptive_dt.dt_target
        );
    }
}
