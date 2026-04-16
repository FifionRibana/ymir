//! Top-level orchestration of the thin viscous sheet simulation.

use tracing::{debug, info, info_span, warn};

use super::advection::{compute_cfl_dt, compute_divergence_flux};
use super::config::{NonlinearSolver, TectonicsConfig};
use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::newton::solve_velocity_newton;
use super::picard::solve_velocity_picard;
use super::traction::TractionField;
use super::workspace::{SolverWorkspace, StepStats};
use crate::tectonics::boundaries::{
    BoundaryType, accumulate_subducted_mass, apply_slab_pull, compute_boundary_sources,
    gaussian_blur_f64,
};
use crate::tectonics::plates::{
    Plate, PlateType, advect_plate_ids, apply_subduction_consumption, compute_viscosity_multiplier,
    detect_disappeared_plates, detect_fragmentation, rebuild_traction, update_plate_stats,
};

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
}

/// Data passed to the progress callback after each timestep.
pub struct StepSnapshot<'a> {
    pub s_field: &'a Field2D,
    pub plate_ids: Option<&'a [usize]>,
    pub plates: Option<&'a [Plate]>,
    pub boundary_types: Option<&'a [BoundaryType]>,
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
    let n = grid.n;

    // Initialize density field and cratonic viscosity multiplier from plate types
    assign_density_from_plates(grid, &plate_ctx.ids, &plate_ctx.plates, &config.boundaries);
    compute_viscosity_multiplier(grid, &plate_ctx.ids, &plate_ctx.plates, &config.cratonic);
    grid.basal_friction = config.basal_friction;

    for step in 0..config.num_timesteps {
        // 1. Solve velocity — continuation only on first step (cold start)
        let need_continuation =
            config.continuation.enabled && config.picard.power_law_n > 1.0 && step == 0;

        let _step_span = info_span!("solver_step", step, n = grid.n).entered();

        let traction = &plate_ctx.traction;
        let (converged, nl_iterations, linear_iterations) = if need_continuation {
            solve_with_continuation(grid, traction, config, workspace)
        } else {
            let result = solve_velocity_direct(grid, traction, config, workspace);
            // Fallback: if direct solve fails, try continuation as recovery
            if !result.0 && config.continuation.enabled && config.picard.power_law_n > 1.0 {
                warn!(step, "direct solve failed, falling back to continuation");
                solve_with_continuation(grid, traction, config, workspace)
            } else {
                result
            }
        };

        if !converged {
            return Err(SolverError::NonlinearSolverDidNotConverge { step });
        }

        // 2. CFL timestep (after velocity is known)
        let dt_cfl = compute_cfl_dt(grid, config.cfl_factor);

        // 2b. Accumulate plastic strain (after solve, with known dt)
        if config.yielding.enabled {
            super::picard::accumulate_plastic_strain(
                &workspace.strain_rate,
                &workspace.eta,
                &config.yielding,
                dt_cfl,
                &mut grid.plastic_strain,
            );
        }

        // ── Dynamic boundary update ────────────────────────────────────
        if config.dynamic_boundaries {
            // Advect plate IDs as material property (replaces seed advection + Voronoï)
            advect_plate_ids(&mut plate_ctx.ids, grid, dt_cfl);

            let disappeared = detect_disappeared_plates(&plate_ctx.ids, &mut plate_ctx.plates);
            for pid in &disappeared {
                info!(plate_id = pid, step, "plate disappeared");
            }

            // Recompute plate statistics from advected IDs
            update_plate_stats(&plate_ctx.ids, &mut plate_ctx.plates, grid);

            // Assign density from material thickness, not plate type
            assign_density_from_material(grid, &config.boundaries);
            compute_viscosity_multiplier(grid, &plate_ctx.ids, &plate_ctx.plates, &config.cratonic);
        }
        // ──────────────────────────────────────────────────────────────

        // 3. Adaptive timestep with retry on excessive clamping
        let s_backup: Vec<f64> = grid.s.data().to_vec();
        let mut dt = dt_cfl;
        let dt_min = dt * 0.01;
        let mut clamp_ratio = 0.0;

        // Compute boundary source terms if enabled
        let boundaries_active = config.boundaries.enabled;
        if boundaries_active {
            let bf = compute_boundary_sources(
                grid,
                &plate_ctx.ids,
                &plate_ctx.plates,
                &config.boundaries,
            );
            workspace.source_rate.data_mut().copy_from_slice(bf.source_rate.data());
            workspace.boundary_field = Some(bf);
            if config.boundaries.source_smoothing_sigma > 0.0 {
                let smoothed = gaussian_blur_f64(
                    &workspace.source_rate,
                    config.boundaries.source_smoothing_sigma,
                );
                workspace.source_rate.data_mut().copy_from_slice(smoothed.data());
            }

            // Slab pull: accumulate subducted mass and update plate velocities
            if config.boundaries.slab_pull_enabled {
                accumulate_subducted_mass(
                    &workspace.source_rate,
                    &plate_ctx.ids,
                    &mut plate_ctx.plates,
                    dt_cfl,
                    n,
                );
                apply_slab_pull(
                    &mut plate_ctx.plates,
                    config.boundaries.slab_pull_factor,
                    config.boundaries.max_plate_velocity,
                );
            }
        }

        // Subduction consumption + fragmentation detection
        if config.dynamic_boundaries {
            if let Some(ref bf) = workspace.boundary_field {
                apply_subduction_consumption(
                    &mut plate_ctx.ids,
                    grid,
                    bf,
                    config.boundaries.subduction_consumption_threshold,
                );
            }
            // Detect continental breakup: check if any plate has been split
            // into disconnected components by a rift zone (thin band of S < threshold).
            if step % 10 == 0 {
                detect_fragmentation(
                    &mut plate_ctx.ids,
                    &mut plate_ctx.plates,
                    &mut plate_ctx.next_id,
                    n,
                    grid,
                    0.25, // Only truly oceanic-thin rift zones break continental connectivity
                );
            }
            // Recompute stats after ID changes
            update_plate_stats(&plate_ctx.ids, &mut plate_ctx.plates, grid);
        }

        // Rebuild traction once (after slab pull may have updated velocities)
        if config.boundaries.enabled || config.dynamic_boundaries {
            plate_ctx.traction = rebuild_traction(&plate_ctx.ids, &plate_ctx.plates, n);
        }

        // ── Mass balance tracking ──────────────────────────────────
        let mass_before: f64 = grid.s.data().iter().sum();
        let mut mass_after_advection = 0.0_f64;
        let mut mass_clamp_delta = 0.0_f64;
        let mut total_div_flux = 0.0_f64;
        let mut total_source = 0.0_f64;

        for _retry in 0..5 {
            compute_divergence_flux(grid, &mut workspace.div_flux);

            mass_after_advection = 0.0;
            mass_clamp_delta = 0.0;
            total_div_flux = 0.0;
            total_source = 0.0;

            for j in 0..n {
                for i in 0..n {
                    let div = workspace.div_flux.get(i, j);
                    let q = if boundaries_active { workspace.source_rate.get(i, j) } else { 0.0 };
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
            clamp_ratio = clamp_count as f64 / (n * n) as f64;

            if clamp_ratio < 0.05 || dt <= dt_min {
                break;
            }

            // Too many cells clamped — retry with smaller dt
            dt *= 0.5;
            grid.s.data_mut().copy_from_slice(&s_backup);
        }

        // 3b. Oceanic restoring force — dense oceanic crust (ρ ≈ 3000 kg/m³)
        // that thickens beyond its equilibrium value becomes gravitationally
        // unstable relative to the underlying mantle, and sinks back down.
        // This is a proxy for thermogravitational instability (cold, dense
        // lithosphere subducting spontaneously) which the thin sheet cannot
        // model directly due to the absence of a vertical dimension.
        // The rate is scaled by dt for physical consistency across different
        // timestep sizes. Only applies to cells below the thickness threshold
        // to avoid eroding legitimate continental margins that happen to sit
        // on an oceanic plate_id.
        if config.boundaries.enabled && config.boundaries.oceanic_restore_rate > 0.0 {
            let rate = config.boundaries.oceanic_restore_rate;
            let s_ref = config.boundaries.oceanic_reference_thickness;
            let s_thr = config.boundaries.oceanic_restore_threshold;
            for j in 0..n {
                for i in 0..n {
                    let pid = plate_ctx.ids[j * n + i];
                    if plate_ctx.plates[pid].plate_type == PlateType::Oceanic {
                        let s_current = grid.s.get(i, j);
                        if s_current > s_ref && s_current < s_thr {
                            let s_new = s_current - dt * rate * (s_current - s_ref);
                            grid.s.set(i, j, s_new.max(config.s_min));
                        }
                    }
                }
            }
        }

        let mass_after_oceanic_restore: f64 = grid.s.data().iter().sum();
        let oceanic_restore_delta = mass_after_oceanic_restore - mass_after_advection;

        // 3c. Continental minimum thickness protection.
        // Prevents continental crust from thinning to oblivion by modeling
        // buoyancy-driven resistance. Cells below continental_restore_threshold
        // are considered fully rifted (new ocean floor) and not restored.
        if config.boundaries.enabled && config.boundaries.continental_restore_rate > 0.0 {
            let rate = config.boundaries.continental_restore_rate;
            let s_min = config.boundaries.continental_min_thickness;
            let s_thr = config.boundaries.continental_restore_threshold;
            for j in 0..n {
                for i in 0..n {
                    let pid = plate_ctx.ids[j * n + i];
                    if plate_ctx.plates[pid].plate_type == PlateType::Continental {
                        let s_current = grid.s.get(i, j);
                        if s_current < s_min && s_current > s_thr {
                            let s_new = s_current + dt * rate * (s_min - s_current);
                            grid.s.set(i, j, s_new);
                        }
                    }
                }
            }
        }

        let mass_after_continental_restore: f64 = grid.s.data().iter().sum();
        let continental_restore_delta = mass_after_continental_restore - mass_after_oceanic_restore;
        let mass_after = mass_after_continental_restore;

        // 4. Update stats
        let mut max_v = 0.0_f64;
        let mut max_s = f64::NEG_INFINITY;
        let mut min_s = f64::INFINITY;
        for j in 0..n {
            for i in 0..n {
                let vx = grid.vx.get(i, j);
                let vy = grid.vy.get(i, j);
                max_v = max_v.max((vx * vx + vy * vy).sqrt());
                let s = grid.s.get(i, j);
                max_s = max_s.max(s);
                min_s = min_s.min(s);
            }
        }
        workspace.stats = StepStats {
            max_velocity: max_v,
            max_thickness: max_s,
            min_thickness: min_s,
            picard_iterations: nl_iterations,
            cg_iterations_last: linear_iterations,
            dt,
            clamp_ratio,
        };

        info!(
            step,
            dt,
            nl_iters = nl_iterations,
            lin_iters = linear_iterations,
            max_velocity = max_v,
            min_thickness = min_s,
            max_thickness = max_s,
            clamp_ratio,
            "tectonic step"
        );

        if clamp_ratio > 0.05 {
            warn!(step, clamp_ratio, dt, "excessive clamping — consider reducing CFL or timestep");
        }

        // ── Mass balance diagnostic (always runs) ──────────────────
        debug!(
            step,
            mass_total = %format!("{:.4}", mass_after),
            mass_delta = %format!("{:.6}", mass_after - mass_before),
            advection_div = %format!("{:.6}", -total_div_flux),
            sources_q = %format!("{:.6}", total_source),
            clamping = %format!("{:.6}", mass_clamp_delta),
            oceanic_restore = %format!("{:.6}", oceanic_restore_delta),
            continental_restore = %format!("{:.6}", continental_restore_delta),
            "mass balance"
        );

        // ── Diagnostic: per-step plate & boundary summary ──────────────
        if config.dynamic_boundaries {
            let total_cells = n * n;
            let num_active =
                plate_ctx.plates.iter().filter(|p| p.active && p.cell_count > 0).count();

            // Per-plate summary
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

            // Boundary type counts
            let (mut n_sub, mut n_osub, mut n_coll, mut n_rift, mut n_boundary) =
                (0usize, 0usize, 0usize, 0usize, 0usize);
            if let Some(ref bf) = workspace.boundary_field {
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

            // Source term summary
            let mut q_pos = 0.0_f64;
            let mut q_neg = 0.0_f64;
            for &q in workspace.source_rate.data().iter() {
                if q > 0.0 {
                    q_pos += q;
                } else {
                    q_neg += q;
                }
            }

            // Thickness distribution
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
        // ──────────────────────────────────────────────────────────────

        // 5. Callback — returns false to cancel
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
    let n = grid.n;
    for j in 0..n {
        for i in 0..n {
            let pid = ids[j * n + i];
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
    let n = grid.n;
    for j in 0..n {
        for i in 0..n {
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
    workspace: &mut SolverWorkspace,
) -> (bool, usize, usize) {
    let rho_c = config.boundaries.rho_continental;
    let rho_m = config.boundaries.rho_mantle;

    let do_solve =
        |grid: &mut StaggeredGrid, ws: &mut SolverWorkspace| match config.nonlinear_solver {
            NonlinearSolver::Picard => {
                let r = solve_velocity_picard(
                    grid,
                    plates,
                    config.gravity_factor,
                    rho_c,
                    rho_m,
                    &config.picard,
                    &config.yielding,
                    ws,
                );
                (r.converged, r.iterations, r.total_cg_iterations)
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
                    ws,
                );
                (r.converged, r.iterations, r.total_linear_iterations)
            }
        };

    let (converged, nl, linear) = do_solve(grid, workspace);

    if !converged && grid.basal_friction > 0.0 {
        // Fallback: solve without friction, then with friction
        let target_friction = grid.basal_friction;
        grid.basal_friction = 0.0;
        let (_, nl2, lin2) = do_solve(grid, workspace);
        grid.basal_friction = target_friction;
        let (conv2, nl3, lin3) = do_solve(grid, workspace);
        return (conv2, nl + nl2 + nl3, linear + lin2 + lin3);
    }

    (converged, nl, linear)
}

/// Solve velocity using viscosity continuation: ramp n from 1 → target,
/// then optionally add a friction continuation step.
fn solve_with_continuation(
    grid: &mut StaggeredGrid,
    plates: &TractionField,
    config: &TectonicsConfig,
    workspace: &mut SolverWorkspace,
) -> (bool, usize, usize) {
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
        let (converged, iters, linear_iters) = match config.nonlinear_solver {
            NonlinearSolver::Picard => {
                let r = solve_velocity_picard(
                    grid,
                    plates,
                    config.gravity_factor,
                    rho_c,
                    rho_m,
                    &step_config,
                    &config.yielding,
                    workspace,
                );
                (r.converged, r.iterations, r.total_cg_iterations)
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
                    workspace,
                );
                (r.converged, r.iterations, r.total_linear_iterations)
            }
        };
        total_nl += iters;
        total_linear += linear_iters;
        if !converged {
            grid.basal_friction = target_friction;
            return (false, total_nl, total_linear);
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

        let (converged, iters, linear_iters) = match config.nonlinear_solver {
            NonlinearSolver::Picard => {
                let r = solve_velocity_picard(
                    grid,
                    plates,
                    config.gravity_factor,
                    rho_c,
                    rho_m,
                    &step_config,
                    &config.yielding,
                    workspace,
                );
                (r.converged, r.iterations, r.total_cg_iterations)
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
                    workspace,
                );
                (r.converged, r.iterations, r.total_linear_iterations)
            }
        };

        total_nl += iters;
        total_linear += linear_iters;

        if !converged {
            return (false, total_nl, total_linear);
        }
    } else {
        grid.basal_friction = target_friction;
    }

    (true, total_nl, total_linear)
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
        }
    }

    #[test]
    fn convergent_plates_thicken() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let traction = TractionField::two_plates_convergent(n, 1.0);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(50);
        let mut ws = SolverWorkspace::new(n);

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
        let mut grid = StaggeredGrid::new(n, dx);
        let s_initial = 1.0;
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, s_initial);
            }
        }

        let traction = TractionField::two_plates_divergent(n, 1.0);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(50);
        let mut ws = SolverWorkspace::new(n);

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
        let mut grid = StaggeredGrid::new(n, dx);

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

        let traction = TractionField::zero(n);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(100);
        let mut ws = SolverWorkspace::new(n);

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
        let mut grid = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let traction = TractionField::two_plates_convergent(n, 0.5);
        let mut ctx = make_static_ctx(n, traction);
        let config = make_config(30);
        let mut ws = SolverWorkspace::new(n);

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
        let mut grid = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let traction = TractionField::two_plates_convergent(n, 0.5);
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
        };

        let mut ctx = make_static_ctx(n, traction);
        let mut ws = SolverWorkspace::new(n);
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
        let mut grid = StaggeredGrid::new(n, dx);

        // All oceanic, initially thick (as if advection has thickened them)
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 0.8);
            }
        }

        let traction = TractionField::zero(n);
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

        let mut ws = SolverWorkspace::new(n);
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
        let mut strain = Field2D::new(n);
        let mut eta = Field2D::new(n);

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
        let mut grid_c = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_c.s.set(i, j, if i < n / 2 { 1.0 } else { 0.5 });
                grid_c.rho.set(i, j, 2750.0);
            }
        }

        // Setup 2: same thickness, all oceanic density
        let mut grid_o = StaggeredGrid::new(n, dx);
        for j in 0..n {
            for i in 0..n {
                grid_o.s.set(i, j, if i < n / 2 { 1.0 } else { 0.5 });
                grid_o.rho.set(i, j, 3000.0);
            }
        }

        let plates = TractionField::zero(n);
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
}
