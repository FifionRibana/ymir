//! Diagnostic statistics helpers for the Phase 1-bis instrumentation
//! of issue #75.
//!
//! Emits tracing `debug!` events on the dedicated targets
//! `rhs_breakdown`, `eta_breakdown`, `residual_spatial`, and
//! `phase_timings`. Expensive per-field percentile computations are
//! gated with `tracing::enabled!(target, Level::DEBUG)` so the solver
//! pays no cost when those targets are filtered out by the subscriber.
//!
//! **Default behaviour:** a standard info-level subscriber will see
//! only the `phase_timings` info events (cheap, always useful). The
//! debug-level breakdowns fire only when the scenario runner (or a
//! user-configured subscriber) enables the diagnostic targets.

use tracing::debug;

use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::traction::TractionField;
use crate::tectonics::boundaries::{BoundaryField, BoundaryType};

/// Phase 2-bis gamma_slab field stats (#75): min, max, mean of γ on
/// margin cells (where `boundary_type ∈ Subduction|OceanicSubduction`),
/// plus the plain-grid max velocity and max velocity restricted to
/// those same margin cells. Fires at debug level on target
/// `slab_pull_sweep` once per macro step; cheap enough to leave on.
pub fn emit_slab_pull_sweep(boundary_field: &BoundaryField, grid: &StaggeredGrid) {
    if !tracing::enabled!(target: "slab_pull_sweep", tracing::Level::DEBUG) {
        return;
    }
    let n2 = boundary_field.n * boundary_field.n;
    let mut g_min = f64::INFINITY;
    let mut g_max = 0.0_f64;
    let mut g_sum = 0.0_f64;
    let mut g_count = 0usize;
    let mut margin_v_max = 0.0_f64;
    let mut v_global_max = 0.0_f64;

    let nx = grid.nx();
    let ny = grid.ny();
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();

    for k in 0..n2.min(boundary_field.gamma_slab.data().len()) {
        let j = k / nx;
        let i = k % nx;
        let vxc = 0.5 * (grid.vx.get(i, j) + grid.vx.get(idx_x.next(i), j));
        let vyc = 0.5 * (grid.vy.get(i, j) + grid.vy.get(i, idx_y.next(j)));
        let v_mag = (vxc * vxc + vyc * vyc).sqrt();
        if v_mag > v_global_max {
            v_global_max = v_mag;
        }
        let is_margin = matches!(
            boundary_field.boundary_type[k],
            BoundaryType::Subduction | BoundaryType::OceanicSubduction
        );
        if is_margin {
            let g = boundary_field.gamma_slab.data()[k];
            if g > g_max {
                g_max = g;
            }
            if g < g_min {
                g_min = g;
            }
            g_sum += g;
            g_count += 1;
            if v_mag > margin_v_max {
                margin_v_max = v_mag;
            }
        }
    }

    let g_mean = if g_count > 0 { g_sum / g_count as f64 } else { 0.0 };
    let g_min_out = if g_count > 0 { g_min } else { 0.0 };

    debug!(
        target: "slab_pull_sweep",
        gamma_margin_min = g_min_out,
        gamma_margin_max = g_max,
        gamma_margin_mean = g_mean,
        margin_cells = g_count,
        margin_v_max,
        v_global_max,
        "slab-pull sweep"
    );
}

/// Floor for ratio denominators. Keeps spike_ratio finite when the
/// median of a field is identically zero (e.g. `T_plates` on a
/// scenario where boundaries are disabled).
const DENOM_FLOOR: f64 = 1e-20;

/// L2 norm, max |x|, and two percentiles (p50, p95) of a slice.
#[derive(Debug, Clone, Copy)]
struct ScalarDist {
    norm: f64,
    max_abs: f64,
    max_cell: usize,
    p50: f64,
    p95: f64,
}

fn scalar_dist(values: &[f64]) -> ScalarDist {
    let mut norm_sq = 0.0_f64;
    let mut max_abs = 0.0_f64;
    let mut max_cell = 0usize;
    for (i, &v) in values.iter().enumerate() {
        norm_sq += v * v;
        let a = v.abs();
        if a > max_abs {
            max_abs = a;
            max_cell = i;
        }
    }

    // Sort absolute values to get percentiles. O(N log N) is fine
    // for 64² or even 512² — this runs once per instrumentation call.
    let mut abs: Vec<f64> = values.iter().map(|v| v.abs()).collect();
    abs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let pick = |p: f64| -> f64 {
        if abs.is_empty() {
            return 0.0;
        }
        let idx = ((abs.len() - 1) as f64 * p).round() as usize;
        abs[idx.min(abs.len() - 1)]
    };

    ScalarDist { norm: norm_sq.sqrt(), max_abs, max_cell, p50: pick(0.50), p95: pick(0.95) }
}

/// RHS breakdown: recompute GPE and T_plates contributions separately
/// on the same stencil used by `compute_rhs`, then emit a debug event.
///
/// This must match the code path in `stokes::compute_rhs` verbatim
/// (same interpolation, same null-space behaviour) so the summed
/// metrics agree with the total. The recomputation is intentional —
/// it keeps `compute_rhs` itself untouched in the non-diagnostic case.
pub fn emit_rhs_breakdown(
    grid: &StaggeredGrid,
    plates: &TractionField,
    gravity_factor: f64,
    rho_continental: f64,
    rho_mantle: f64,
) {
    // Skip the expensive decomposition if the subscriber filters
    // this target out. Step/nx/ny come from the enclosing
    // `solver_step` tracing span and do not need to be passed in.
    if !tracing::enabled!(target: "rhs_breakdown", tracing::Level::DEBUG) {
        return;
    }

    let nx = grid.nx();
    let ny = grid.ny();
    let n2 = nx * ny;
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();
    let inv_dx = 1.0 / grid.dx;
    let use_density = rho_mantle > 0.0;

    // GPE potential at a cell (same as stokes::compute_rhs).
    let gpe = |i: usize, j: usize| -> f64 {
        let s = grid.s.get(i, j);
        if use_density {
            let rho = grid.rho.get(i, j);
            let buoyancy = rho * (1.0 - rho / rho_mantle);
            let ref_buoyancy = rho_continental * (1.0 - rho_continental / rho_mantle);
            (buoyancy / ref_buoyancy) * s * s
        } else {
            s * s
        }
    };

    let mut gpe_only = vec![0.0_f64; 2 * n2];
    let mut tp_only = vec![0.0_f64; 2 * n2];

    for j in 0..ny {
        let pj = idx_y.prev(j);
        for i in 0..nx {
            let pi = idx_x.prev(i);
            let k = j * nx + i;

            let dpdx = -gravity_factor * (gpe(i, j) - gpe(pi, j)) * inv_dx;
            let tx = 0.5 * (plates.tx.get(pi, j) + plates.tx.get(i, j));
            gpe_only[k] = dpdx;
            tp_only[k] = tx;

            let dpdy = -gravity_factor * (gpe(i, j) - gpe(i, pj)) * inv_dx;
            let ty = 0.5 * (plates.ty.get(i, pj) + plates.ty.get(i, j));
            gpe_only[n2 + k] = dpdy;
            tp_only[n2 + k] = ty;
        }
    }

    // Combined total (avoids a third pass).
    let total: Vec<f64> = gpe_only.iter().zip(tp_only.iter()).map(|(a, b)| a + b).collect();

    let gpe_stats = scalar_dist(&gpe_only);
    let tp_stats = scalar_dist(&tp_only);
    let tot_stats = scalar_dist(&total);

    let gpe_spike = gpe_stats.max_abs / gpe_stats.p50.max(DENOM_FLOOR);
    let tp_spike = tp_stats.max_abs / tp_stats.p50.max(DENOM_FLOOR);
    let tot_spike = tot_stats.max_abs / tot_stats.p50.max(DENOM_FLOOR);

    debug!(
        target: "rhs_breakdown",
        gpe_rhs_norm = gpe_stats.norm,
        gpe_rhs_max_abs = gpe_stats.max_abs,
        gpe_rhs_max_cell = gpe_stats.max_cell,
        gpe_rhs_p95 = gpe_stats.p95,
        gpe_rhs_p50 = gpe_stats.p50,
        gpe_rhs_spike_ratio = gpe_spike,
        tplates_rhs_norm = tp_stats.norm,
        tplates_rhs_max_abs = tp_stats.max_abs,
        tplates_rhs_max_cell = tp_stats.max_cell,
        tplates_rhs_p95 = tp_stats.p95,
        tplates_rhs_p50 = tp_stats.p50,
        tplates_rhs_spike_ratio = tp_spike,
        total_rhs_norm = tot_stats.norm,
        total_rhs_max_abs = tot_stats.max_abs,
        total_rhs_max_cell = tot_stats.max_cell,
        total_rhs_p95 = tot_stats.p95,
        total_rhs_p50 = tot_stats.p50,
        total_rhs_spike_ratio = tot_spike,
        "rhs breakdown"
    );
}

/// Viscosity breakdown: counts yielded and near-saturation cells by
/// recomputing `eta_plastic` per cell from the same inputs
/// `apply_yielding` uses.
///
/// Emit once per nonlinear outer iteration, after `apply_yielding`.
pub fn emit_eta_breakdown(
    newton_iter: usize,
    eta: &Field2D,
    strain_rate: &Field2D,
    plastic_strain: &Field2D,
    yielding: &super::config::YieldingConfig,
    eta_max: f64,
) {
    if !tracing::enabled!(target: "eta_breakdown", tracing::Level::DEBUG) {
        return;
    }

    let eta_data = eta.data();
    if eta_data.is_empty() {
        return;
    }

    let mut sorted: Vec<f64> = eta_data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = sorted.len();
    let pick = |p: f64| sorted[((n - 1) as f64 * p).round() as usize];

    let eta_min_actual = sorted[0];
    let eta_max_actual = sorted[n - 1];
    let eta_p05 = pick(0.05);
    let eta_p50 = pick(0.50);
    let eta_p95 = pick(0.95);
    let eta_ratio = eta_max_actual / eta_min_actual.max(DENOM_FLOOR);

    // Count yielded and saturated cells.
    let sr = strain_rate.data();
    let ps = plastic_strain.data();
    let yielding_active = yielding.enabled && yielding.yield_stress > 0.0;

    let saturation_threshold = 0.95 * eta_max;
    let mut yielded = 0usize;
    let mut saturated = 0usize;
    for k in 0..n {
        if eta_data[k] >= saturation_threshold {
            saturated += 1;
        }
        if yielding_active {
            let sr_k = sr[k];
            if sr_k < 1e-20 {
                continue;
            }
            let local_yield = if yielding.weakening_enabled {
                let w =
                    (ps[k] / yielding.weakening_strain_ref).min(1.0) * yielding.weakening_fraction;
                yielding.yield_stress * (1.0 - w)
            } else {
                yielding.yield_stress
            };
            let eta_plastic = local_yield / (2.0 * sr_k);
            if eta_data[k] < 1.01 * eta_plastic {
                yielded += 1;
            }
        }
    }

    debug!(
        target: "eta_breakdown",
        newton_iter,
        eta_min_actual,
        eta_max_actual,
        eta_ratio,
        eta_p05,
        eta_p50,
        eta_p95,
        yielding_cells_count = yielded,
        yielding_cells_fraction = (yielded as f64) / (n as f64),
        saturated_cells_count = saturated,
        saturated_cells_fraction = (saturated as f64) / (n as f64),
        "eta breakdown"
    );
}

/// Residual localization: fraction of L2 energy of `F(v_converged)`
/// concentrated in cells flagged as a non-None boundary type.
///
/// Residual is the flat `[vx; vy]` layout of length `2·nx·ny`; the
/// per-cell contribution is `|F_vx(i,j)|² + |F_vy(i,j)|²`.
pub fn emit_residual_spatial(
    residual: &[f64],
    boundary_field: Option<&BoundaryField>,
    nx: usize,
    ny: usize,
) {
    if !tracing::enabled!(target: "residual_spatial", tracing::Level::DEBUG) {
        return;
    }

    let n2 = nx * ny;
    if residual.len() != 2 * n2 {
        return;
    }

    let mut per_cell: Vec<f64> = Vec::with_capacity(n2);
    let mut total_energy = 0.0_f64;
    let mut max_per_cell = 0.0_f64;
    let mut max_cell = 0usize;

    for k in 0..n2 {
        let e = residual[k] * residual[k] + residual[n2 + k] * residual[n2 + k];
        per_cell.push(e);
        total_energy += e;
        if e > max_per_cell {
            max_per_cell = e;
            max_cell = k;
        }
    }

    // p99 over absolute |F| (sqrt of per-cell energy to match "residual
    // magnitude" intuition).
    let mut abs_mag: Vec<f64> = per_cell.iter().map(|e| e.sqrt()).collect();
    abs_mag.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p99 = if !abs_mag.is_empty() {
        abs_mag[((abs_mag.len() - 1) as f64 * 0.99).round() as usize]
    } else {
        0.0
    };

    let (boundary_energy, boundary_cell_count) = match boundary_field {
        Some(bf) => {
            let mut be = 0.0_f64;
            let mut count = 0usize;
            for k in 0..n2 {
                if bf.boundary_type[k] != BoundaryType::None {
                    be += per_cell[k];
                    count += 1;
                }
            }
            (be, count)
        }
        None => (0.0, 0),
    };

    let localization =
        if total_energy > DENOM_FLOOR { boundary_energy / total_energy } else { 0.0 };

    let boundary_fraction = (boundary_cell_count as f64) / (n2 as f64);

    debug!(
        target: "residual_spatial",
        residual_norm_total = total_energy.sqrt(),
        residual_max_cell = max_cell,
        residual_max_abs = max_per_cell.sqrt(),
        residual_p99 = p99,
        residual_localization = localization,
        boundary_cell_fraction = boundary_fraction,
        "residual spatial"
    );
}
