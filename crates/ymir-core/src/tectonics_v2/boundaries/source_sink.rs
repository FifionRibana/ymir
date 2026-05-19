//! Source/sink computation for Step 5.
//!
//! Computes the per-cell `Q̃(S̃, ṽ)` rate field from a
//! plate-type field, a boundary-flag field, and rate coefficients
//! (see issue #89 D2). The five terms are:
//!
//! - `Q_sub   = -k_sub · |Δṽ_conv|` on oceanic cells flagged as
//!   (Oceanic)Subduction.
//! - `Q_arc   = +k_arc · Σ_neighbours_subducting |Q_sub_neighbour|`
//!   on continental cells with at least one subducting neighbour.
//! - `Q_spread= +k_spread` constant on cells flagged `Rift` that are
//!   oceanic (mid-ocean ridge).
//! - `Q_coll-v= +k_coll_v` on continental `ContinentalCollision` cells.
//! - `Q_rift-v= +k_rift_v` on continental `Rift` cells.
//!
//! The divergence `div(v)_cell` is computed directly from MAC face
//! velocities with no interpolation (see [`div_v_cell`]), and
//! `|Δṽ_conv|_cell = max(0, -div(v)_cell)` extracts the convergent
//! component that drives subduction consumption. Divergent regions
//! (`max(0, +div)`) are tracked for diagnostic use; `Q_spread` is
//! constant at Step 5 and not modulated by `|Δṽ_div|` — the flag
//! field carries the "where spreading happens" information, and
//! modulation by local divergence is Step 6 work.

use super::boundary_flag::{BoundaryFlag, BoundaryFlagField, BoundaryRates};
use super::plate_type::{PlateType, PlateTypeField};
use super::super::field::{Field2D, PeriodicIndex};

/// Compute `div(v)_cell[i,j] = (vx[i+1,j]-vx[i,j])/dx + (vy[i,j+1]-vy[i,j])/dy`
/// on the staggered MAC grid, directly from face velocities.
///
/// `vx[i,j]` is the velocity at the **left** vertical face of cell
/// `(i,j)` — the same face indexing convention used throughout
/// `tectonics_v2` (see [`super::super::stokes::operator`]). So the
/// right face of cell `(i,j)` is `vx[i+1,j]`, wrapped periodically.
pub fn div_v_cell(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    vx: &[f64],
    vy: &[f64],
    out: &mut Field2D,
) {
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);
    debug_assert_eq!(out.nx(), nx);
    debug_assert_eq!(out.ny(), ny);
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let lin = |i: usize, j: usize| j * nx + i;
    for j in 0..ny {
        let jp = idx_y.next(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let d = (vx[lin(ip, j)] - vx[lin(i, j)]) * inv_dx
                + (vy[lin(i, jp)] - vy[lin(i, j)]) * inv_dy;
            out.set(i, j, d);
        }
    }
}

/// `|Δṽ_conv|_cell = max(0, -div(v)_cell)`. Non-allocating:
/// writes into `out`.
pub fn convergent_component(div_v: &Field2D, out: &mut Field2D) {
    debug_assert_eq!(div_v.nx(), out.nx());
    debug_assert_eq!(div_v.ny(), out.ny());
    for (s, d) in out.data_mut().iter_mut().zip(div_v.data().iter()) {
        *s = (-*d).max(0.0);
    }
}

/// Compute the combined source/sink rate field `Q̃(i,j)` from the
/// boundary state. Writes into `out`; cells not touched by any
/// mechanism receive 0.
///
/// Evaluation order matches issue #89 D2:
/// 1. First pass: `Q_sub` (needs conv component), `Q_spread`,
///    `Q_coll-v`, `Q_rift-v`.
/// 2. Second pass: `Q_arc` (needs `Q_sub` already computed on
///    neighbouring cells).
///
/// The two-pass structure is required: `Q_arc[i,j]` depends on
/// `|Q_sub[i±1,j]|` and `|Q_sub[i,j±1]|`, so `Q_sub` must exist on
/// the whole grid before the arc loop runs.
pub fn compute_source_sink_terms(
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
    rates: &BoundaryRates,
    div_v: &Field2D,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    q_sub_scratch: &mut Field2D,
    out: &mut Field2D,
) {
    let nx = plate_types.nx();
    let ny = plate_types.ny();
    debug_assert_eq!(flags.nx(), nx);
    debug_assert_eq!(flags.ny(), ny);
    debug_assert_eq!(div_v.nx(), nx);
    debug_assert_eq!(div_v.ny(), ny);
    debug_assert_eq!(out.nx(), nx);
    debug_assert_eq!(out.ny(), ny);
    debug_assert_eq!(q_sub_scratch.nx(), nx);
    debug_assert_eq!(q_sub_scratch.ny(), ny);

    // Pass 1 — per-cell terms that don't need neighbour Q values.
    for j in 0..ny {
        for i in 0..nx {
            let t = plate_types.get(i, j);
            let f = flags.get(i, j);
            let mut q = 0.0_f64;
            let mut q_sub_here = 0.0_f64;

            if f.is_subduction() && matches!(t, PlateType::Oceanic) {
                // |Δv_conv| = max(0, -div(v))
                let conv = (-div_v.get(i, j)).max(0.0);
                q_sub_here = -rates.k_sub * conv;
                q += q_sub_here;
            }
            match (t, f) {
                (PlateType::Oceanic, BoundaryFlag::Rift) => {
                    q += rates.k_spread;
                }
                (PlateType::Continental, BoundaryFlag::Rift) => {
                    q += rates.k_rift_v;
                }
                (PlateType::Continental, BoundaryFlag::ContinentalCollision) => {
                    q += rates.k_coll_v;
                }
                _ => {}
            }
            out.set(i, j, q);
            q_sub_scratch.set(i, j, q_sub_here);
        }
    }

    // Pass 2 — Q_arc on continental cells neighbouring a subducting
    // cell. Four-neighbour stencil with periodic wrap.
    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            if !matches!(plate_types.get(i, j), PlateType::Continental) {
                continue;
            }
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);
            let mut arc_sum = 0.0_f64;
            for (ni, nj) in [(ip, j), (im, j), (i, jp), (i, jm)] {
                if flags.get(ni, nj).is_subduction() {
                    arc_sum += q_sub_scratch.get(ni, nj).abs();
                }
            }
            if arc_sum > 0.0 {
                out.set(i, j, out.get(i, j) + rates.k_arc * arc_sum);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn div_v_matches_analytic_for_sinusoidal_velocity() {
        // vx(x,y) = sin(2πx), vy(x,y) = -sin(2πy).
        // On staggered grid: vx at face (i dx, (j+½)dy), vy at
        // face ((i+½)dx, j dy). With our convention vx[i,j] is the
        // left face of cell (i,j), at x = i·dx; vy[i,j] is the
        // bottom face, at y = j·dy.
        // Expected: div(v)(x,y) = 2π cos(2πx) - 2π cos(2πy),
        // sampled at cell centre ((i+½)dx, (j+½)dy).
        let nx = 64;
        let ny = 64;
        let dx = 1.0 / nx as f64;
        let dy = 1.0 / ny as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x_face = i as f64 * dx;
                let y_face = j as f64 * dy;
                vx[j * nx + i] = (2.0 * PI * x_face).sin();
                vy[j * nx + i] = -(2.0 * PI * y_face).sin();
            }
        }
        let mut div = Field2D::new(nx, ny);
        div_v_cell(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy, &mut div);
        let mut err = 0.0_f64;
        let mut count = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                let xc = (i as f64 + 0.5) * dx;
                let yc = (j as f64 + 0.5) * dy;
                let expected =
                    2.0 * PI * (2.0 * PI * xc).cos() - 2.0 * PI * (2.0 * PI * yc).cos();
                err += (div.get(i, j) - expected).powi(2);
                count += 1;
            }
        }
        let rms = (err / count as f64).sqrt();
        // Finite-difference error decays like O(dx²); at 64² with
        // this test function the RMS is ~0.07 (magnitude of deriv ~2π
        // ≈ 6.28). Bound of 0.2 is ample.
        assert!(rms < 0.2, "div rms error = {}", rms);
    }

    #[test]
    fn convergent_component_extracts_negative_div() {
        let nx = 3;
        let ny = 3;
        let mut div = Field2D::new(nx, ny);
        // div = [-1, 0, +2] in row-major, all rows identical.
        for j in 0..ny {
            for i in 0..nx {
                let v = match i {
                    0 => -1.0,
                    1 => 0.0,
                    _ => 2.0,
                };
                div.set(i, j, v);
            }
        }
        let mut conv = Field2D::new(nx, ny);
        convergent_component(&div, &mut conv);
        for j in 0..ny {
            assert_eq!(conv.get(0, j), 1.0); // -div = +1 → conv = 1
            assert_eq!(conv.get(1, j), 0.0);
            assert_eq!(conv.get(2, j), 0.0); // div = +2 → conv = 0
        }
    }

    #[test]
    fn q_sub_annihilates_on_divergent_subduction_cell() {
        // Subduction flag on a cell whose local divergence is
        // positive (pulling apart): Q_sub must be 0, not -k_sub·|div|.
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let plate_types = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 1, BoundaryFlag::OceanicSubduction);
        let mut div = Field2D::new(nx, ny);
        div.set(1, 1, 0.5); // positive: divergent
        let rates = BoundaryRates::baseline_uncalibrated();
        let mut q = Field2D::new(nx, ny);
        let mut q_sub = Field2D::new(nx, ny);
        compute_source_sink_terms(
            &plate_types,
            &flags,
            &rates,
            &div,
            &idx_x,
            &idx_y,
            &mut q_sub,
            &mut q,
        );
        // Oceanic cell, so no Q_coll-v / Q_rift-v. Oceanic subduction
        // flag but div > 0 → Q_sub = 0. Q_spread only on Rift. Result
        // is exactly 0.
        assert_eq!(q.get(1, 1), 0.0);
        assert_eq!(q_sub.get(1, 1), 0.0);
    }

    #[test]
    fn q_sub_matches_k_sub_times_conv_on_convergent_cell() {
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let plate_types = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 1, BoundaryFlag::OceanicSubduction);
        let mut div = Field2D::new(nx, ny);
        div.set(1, 1, -0.8); // convergent: conv = 0.8
        let rates = BoundaryRates::baseline_uncalibrated();
        let mut q = Field2D::new(nx, ny);
        let mut q_sub = Field2D::new(nx, ny);
        compute_source_sink_terms(
            &plate_types,
            &flags,
            &rates,
            &div,
            &idx_x,
            &idx_y,
            &mut q_sub,
            &mut q,
        );
        let expected = -rates.k_sub * 0.8;
        assert!((q.get(1, 1) - expected).abs() < 1e-14);
        assert!((q_sub.get(1, 1) - expected).abs() < 1e-14);
    }

    #[test]
    fn q_arc_activates_on_continental_neighbour_of_subduction() {
        // 4×4 oceanic domain with a continental cell at (2, 2)
        // neighbouring an oceanic-subduction cell at (1, 2).
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut plate_types = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        plate_types.set(2, 2, PlateType::Continental);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 2, BoundaryFlag::OceanicSubduction);
        let mut div = Field2D::new(nx, ny);
        div.set(1, 2, -1.0); // convergent
        let rates = BoundaryRates::baseline_uncalibrated();
        let mut q = Field2D::new(nx, ny);
        let mut q_sub = Field2D::new(nx, ny);
        compute_source_sink_terms(
            &plate_types,
            &flags,
            &rates,
            &div,
            &idx_x,
            &idx_y,
            &mut q_sub,
            &mut q,
        );
        // Continental cell (2,2) has neighbour (1,2) subduction with
        // |Q_sub| = k_sub. No other subduction neighbours.
        let expected = rates.k_arc * (rates.k_sub * 1.0);
        assert!((q.get(2, 2) - expected).abs() < 1e-14);
    }

    #[test]
    fn q_spread_is_constant_and_oceanic_only() {
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut plate_types = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        plate_types.set(2, 2, PlateType::Continental);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 1, BoundaryFlag::Rift); // oceanic rift
        flags.set(2, 2, BoundaryFlag::Rift); // continental rift
        let div = Field2D::new(nx, ny);
        let rates = BoundaryRates::baseline_uncalibrated();
        let mut q = Field2D::new(nx, ny);
        let mut q_sub = Field2D::new(nx, ny);
        compute_source_sink_terms(
            &plate_types,
            &flags,
            &rates,
            &div,
            &idx_x,
            &idx_y,
            &mut q_sub,
            &mut q,
        );
        // Oceanic rift → Q_spread = k_spread. Continental rift →
        // Q_rift-v = k_rift_v (smaller).
        assert!((q.get(1, 1) - rates.k_spread).abs() < 1e-14);
        assert!((q.get(2, 2) - rates.k_rift_v).abs() < 1e-14);
    }
}
