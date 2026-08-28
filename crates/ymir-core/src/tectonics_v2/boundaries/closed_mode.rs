//! Step 6 Closed-mode source/sink pipeline.
//!
//! The closed path:
//!
//! 1. [`compute_q_sub_only`]: fills `Q_sub` on oceanic subducting
//!    cells (`-k_sub · max(0, -div(v))`). Other cells zeroed.
//! 2. Caller integrates `M_sub_step = Σ |Q_sub| · Δt · cell_area`.
//! 3. Caller updates the [`super::super::recycling::ImmediateAccumulators`]
//!    by `config.X_fraction · M_sub_step` per class.
//! 4. [`distribute_immediate`]: for each class (arc, coll_v,
//!    rift_v), if eligible cells exist, distribute the full
//!    accumulator uniformly over eligible cells and zero the
//!    accumulator. Otherwise, keep the accumulator (rollover).
//! 5. Caller deposits `config.spread_fraction · M_sub_step` into the
//!    [`super::super::recycling::DelayedRecycler`] buffer.
//! 6. [`distribute_delayed`]: advance the buffer; if oceanic rift
//!    cells exist, distribute the emerging mass uniformly. Otherwise,
//!    roll it over via `advance_or_rollover(false)`.
//!
//! Eligibility rules (same as Step 5 Open mode's neighbour lookup):
//!
//! - `arc`: continental cells with ≥ 1 subducting neighbour.
//! - `coll_v`: continental cells with `BoundaryFlag::ContinentalCollision`.
//! - `rift_v`: continental cells with `BoundaryFlag::Rift`.
//! - `spread`: oceanic cells with `BoundaryFlag::Rift`.

use super::super::field::{Field2D, PeriodicIndex};
use super::super::recycling::{DelayedRecycler, ImmediateAccumulators};
use super::boundary_flag::{BoundaryFlag, BoundaryFlagField, BoundaryRates};
use super::plate_type::{PlateType, PlateTypeField};

/// Compute `Q_sub` on every oceanic subduction-flagged cell.
/// `out[i,j] = -k_sub · max(0, -div(v)[i,j])` when eligible, else 0.
/// No other source/sink terms written — caller composes the full Q.
pub fn compute_q_sub_only(
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
    rates: &BoundaryRates,
    div_v: &Field2D,
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

    for j in 0..ny {
        for i in 0..nx {
            let t = plate_types.get(i, j);
            let f = flags.get(i, j);
            let q = if f.is_subduction() && matches!(t, PlateType::Oceanic) {
                let conv = (-div_v.get(i, j)).max(0.0);
                -rates.k_sub * conv
            } else {
                0.0
            };
            out.set(i, j, q);
        }
    }
}

/// Total subducted mass this step: `Σ |Q_sub| · Δt · dA`.
#[inline]
pub fn integrate_sub_mass(q_sub: &Field2D, dt: f64, cell_area: f64) -> f64 {
    let mut sum = 0.0_f64;
    for &v in q_sub.data() {
        sum += v.abs();
    }
    sum * dt * cell_area
}

/// Count cells eligible for each immediate-distribution class.
/// Returns `(arc_count, coll_v_count, rift_v_count)`.
pub fn count_immediate_eligibilities(
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
) -> (usize, usize, usize) {
    let nx = plate_types.nx();
    let ny = plate_types.ny();
    let mut arc = 0usize;
    let mut coll = 0usize;
    let mut rift = 0usize;
    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);
            let t = plate_types.get(i, j);
            if !matches!(t, PlateType::Continental) {
                continue;
            }
            let f = flags.get(i, j);
            // Arc eligibility: continental cell with ≥ 1 subducting
            // neighbour. Excludes own flag — the continental cell
            // itself is `None` under Step 6's detection rule (see
            // `boundary_detection::classify_convergent`).
            let mut has_sub_neighbour = false;
            for (ni, nj) in [(ip, j), (im, j), (i, jp), (i, jm)] {
                if flags.get(ni, nj).is_subduction() {
                    has_sub_neighbour = true;
                    break;
                }
            }
            if has_sub_neighbour {
                arc += 1;
            }
            if matches!(f, BoundaryFlag::ContinentalCollision) {
                coll += 1;
            }
            if matches!(f, BoundaryFlag::Rift) {
                rift += 1;
            }
        }
    }
    (arc, coll, rift)
}

/// Count oceanic cells flagged `Rift` (eligible for spread
/// distribution).
pub fn count_spread_eligibility(plate_types: &PlateTypeField, flags: &BoundaryFlagField) -> usize {
    plate_types
        .data()
        .iter()
        .zip(flags.data().iter())
        .filter(|&(&t, &f)| matches!(t, PlateType::Oceanic) && matches!(f, BoundaryFlag::Rift))
        .count()
}

/// Distribute the arc/coll_v/rift_v accumulator budgets onto the
/// eligible cells for this step, writing per-cell `Q` values into
/// the shared `out` field (added to existing contents — the caller
/// must zero it before calling if needed, or compose it with
/// Q_sub).
///
/// Units: `out[i,j] += M_class_pending / (n_class_active · dt · dA)`
/// so that `Σ out[i,j]_class · dt · dA = M_class_pending`.
///
/// If no eligible cells exist for a class, the accumulator is
/// preserved (rollover); no Q contribution is emitted for that
/// class.
pub fn distribute_immediate(
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
    accumulators: &mut ImmediateAccumulators,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    dt: f64,
    cell_area: f64,
    out: &mut Field2D,
) {
    let (n_arc, n_coll, n_rift) = count_immediate_eligibilities(plate_types, flags, idx_x, idx_y);
    let inv_norm = 1.0 / (dt * cell_area);

    let q_arc_per_cell = if n_arc > 0 {
        let v = accumulators.arc_pending * inv_norm / n_arc as f64;
        accumulators.arc_pending = 0.0;
        v
    } else {
        0.0
    };
    let q_coll_per_cell = if n_coll > 0 {
        let v = accumulators.coll_v_pending * inv_norm / n_coll as f64;
        accumulators.coll_v_pending = 0.0;
        v
    } else {
        0.0
    };
    let q_rift_per_cell = if n_rift > 0 {
        let v = accumulators.rift_v_pending * inv_norm / n_rift as f64;
        accumulators.rift_v_pending = 0.0;
        v
    } else {
        0.0
    };

    let nx = plate_types.nx();
    let ny = plate_types.ny();
    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            if !matches!(plate_types.get(i, j), PlateType::Continental) {
                continue;
            }
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);
            let f = flags.get(i, j);
            let mut add = 0.0_f64;
            if q_arc_per_cell != 0.0 {
                let has_sub_neighbour = flags.get(ip, j).is_subduction()
                    || flags.get(im, j).is_subduction()
                    || flags.get(i, jp).is_subduction()
                    || flags.get(i, jm).is_subduction();
                if has_sub_neighbour {
                    add += q_arc_per_cell;
                }
            }
            if q_coll_per_cell != 0.0 && matches!(f, BoundaryFlag::ContinentalCollision) {
                add += q_coll_per_cell;
            }
            if q_rift_per_cell != 0.0 && matches!(f, BoundaryFlag::Rift) {
                add += q_rift_per_cell;
            }
            if add != 0.0 {
                out.set(i, j, out.get(i, j) + add);
            }
        }
    }
}

/// Distribute the mass emerging from the delayed buffer over
/// oceanic rift cells. Buffer advances with rollover semantics
/// ([`DelayedRecycler::advance_or_rollover`]).
///
/// Returns `M_emerging` (zero if rolled over this step).
pub fn distribute_delayed(
    plate_types: &PlateTypeField,
    flags: &BoundaryFlagField,
    buffer: &mut DelayedRecycler,
    dt: f64,
    cell_area: f64,
    out: &mut Field2D,
) -> f64 {
    let n_spread = count_spread_eligibility(plate_types, flags);
    let can_distribute = n_spread > 0;
    let m_emerging = buffer.advance_or_rollover(can_distribute);
    if !can_distribute || m_emerging == 0.0 {
        return m_emerging;
    }
    let q_per_cell = m_emerging / (dt * cell_area * n_spread as f64);
    let nx = plate_types.nx();
    let ny = plate_types.ny();
    for j in 0..ny {
        for i in 0..nx {
            if matches!(plate_types.get(i, j), PlateType::Oceanic)
                && matches!(flags.get(i, j), BoundaryFlag::Rift)
            {
                out.set(i, j, out.get(i, j) + q_per_cell);
            }
        }
    }
    m_emerging
}

#[cfg(test)]
mod tests {
    use super::super::super::recycling::RecyclingConfig;
    use super::*;

    #[test]
    fn q_sub_only_zero_on_non_subducting_cells() {
        let pt = PlateTypeField::filled(4, 4, PlateType::Oceanic);
        let flags = BoundaryFlagField::filled(4, 4, BoundaryFlag::None);
        let rates = BoundaryRates::baseline_uncalibrated();
        let mut div = Field2D::new(4, 4);
        for v in div.data_mut().iter_mut() {
            *v = -1.0;
        } // convergent everywhere
        let mut out = Field2D::new(4, 4);
        compute_q_sub_only(&pt, &flags, &rates, &div, &mut out);
        for &v in out.data() {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn q_sub_only_fires_on_oceanic_subducting_with_convergence() {
        let pt = PlateTypeField::filled(4, 4, PlateType::Oceanic);
        let mut flags = BoundaryFlagField::filled(4, 4, BoundaryFlag::None);
        flags.set(1, 1, BoundaryFlag::OceanicSubduction);
        let rates = BoundaryRates::baseline_uncalibrated();
        let mut div = Field2D::new(4, 4);
        div.set(1, 1, -0.5);
        let mut out = Field2D::new(4, 4);
        compute_q_sub_only(&pt, &flags, &rates, &div, &mut out);
        let expected = -rates.k_sub * 0.5;
        assert!((out.get(1, 1) - expected).abs() < 1e-14);
    }

    #[test]
    fn distribute_immediate_drains_accumulators_when_eligible() {
        // 4×4 domain, oceanic subducting at (1,1), continental at
        // (2,1) — continental has subducting neighbour to its west
        // → arc-eligible.
        let nx = 4;
        let ny = 4;
        let mut pt = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        pt.set(2, 1, PlateType::Continental);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 1, BoundaryFlag::OceanicSubduction);
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        let mut accs =
            ImmediateAccumulators { arc_pending: 0.1, coll_v_pending: 0.0, rift_v_pending: 0.0 };
        let mut out = Field2D::new(nx, ny);
        let dt = 0.01;
        let cell_area = 1.0 / (nx * ny) as f64;
        distribute_immediate(&pt, &flags, &mut accs, &idx_x, &idx_y, dt, cell_area, &mut out);
        // After distribution: arc_pending = 0, cell (2,1) received
        // Q_arc = 0.1 / (1 · dt · dA).
        assert_eq!(accs.arc_pending, 0.0);
        let expected = 0.1 / (dt * cell_area);
        assert!((out.get(2, 1) - expected).abs() < 1e-12);
    }

    #[test]
    fn distribute_immediate_rolls_over_when_no_eligible_cell() {
        // All oceanic, no continental → arc has no recipient; the
        // pending value must be preserved.
        let nx = 4;
        let ny = 4;
        let pt = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(0, 0, BoundaryFlag::OceanicSubduction);
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut accs =
            ImmediateAccumulators { arc_pending: 0.1, coll_v_pending: 0.05, rift_v_pending: 0.02 };
        let mut out = Field2D::new(nx, ny);
        distribute_immediate(&pt, &flags, &mut accs, &idx_x, &idx_y, 0.01, 0.01, &mut out);
        // No continental cells → arc, coll, rift all roll over.
        assert_eq!(accs.arc_pending, 0.1);
        assert_eq!(accs.coll_v_pending, 0.05);
        assert_eq!(accs.rift_v_pending, 0.02);
        for &v in out.data() {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn distribute_delayed_rolls_over_when_no_rift() {
        let nx = 4;
        let ny = 4;
        let pt = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        let mut buffer = DelayedRecycler::new(3);
        buffer.deposit(0.5);
        // Advance without rift cells — must rollover.
        let mut out = Field2D::new(nx, ny);
        let emerged = distribute_delayed(&pt, &flags, &mut buffer, 0.01, 0.01, &mut out);
        assert_eq!(emerged, 0.0);
        for &v in out.data() {
            assert_eq!(v, 0.0);
        }
        assert!((buffer.fill() - 0.5).abs() < 1e-14);
    }

    #[test]
    fn config_validate_default_passes() {
        RecyclingConfig::default().validate().unwrap();
    }
}
