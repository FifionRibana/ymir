//! Step 10 — event-driven resets of the age field at boundary
//! cells, per `solver-scaling.md` §4.11 / D3:
//!
//! - **Ridge** (`BoundaryFlag::Rift` on oceanic): `A := 0`. Fresh
//!   crust at a mid-ocean ridge inherits the simulation start age
//!   (zero by §4.11 / D3 — "new oceanic crust").
//! - **Continental rift** (`BoundaryFlag::Rift` on continental):
//!   §4.11 specifies the ridge case explicitly; for continental
//!   rifts the geological intuition is also "fresh material" so
//!   we extend the same `A := 0` reset. This is a conservative
//!   choice, easy to revisit if the §4.11 patch refines it.
//! - **Arc** (continental cell adjacent to a subducting cell;
//!   matches the `Q_arc` source-term computation in
//!   [`crate::tectonics_v2::boundaries::source_sink`]):
//!   `A := 0`. Volcanic resurfacing emplaces fresh material.
//! - **Continental collision** (`BoundaryFlag::ContinentalCollision`):
//!   `A := max(A_self, max_{n ∈ N(i,j)} A_n)`, restricted to
//!   continental neighbours. The resulting cell carries the age
//!   of the older protolith ("the scar is as old as the older
//!   protolith"). Per D4 we use `max`, not a mass-weighted
//!   average — the weighted alternative is deferred.
//! - **Subduction** cells: NO action. The cell is consumed by the
//!   `S̃` recycling pipeline; its `A` value is overwritten on the
//!   next step when the cell is replenished from neighbouring
//!   crust via advection. Explicit zero-out would double-count
//!   the consumption.
//!
//! `apply_age_events` is called **after** the advection step (which
//! has already applied the quiescent `+dt` source) and **after**
//! the `S̃` source/sink + clamp pipeline has updated `S̃`. The
//! quiescent growth applied by advection is intentionally
//! overwritten at boundary cells — the per-cell semantics is
//! "either the cell experienced a boundary event (overwrite) or
//! it did not (keep advected value with quiescent growth)".

use crate::tectonics_v2::boundaries::{BoundaryFlag, BoundaryFlagField, PlateType, PlateTypeField};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

use super::init::is_continental_thickness;

/// Per-step counts of event-reset firings, sampled by the harness
/// for the report. All optional under the structural by-pass.
#[derive(Clone, Copy, Debug, Default)]
pub struct AgeEventCounts {
    pub ridge_resets: u32,
    pub arc_resets: u32,
    pub collision_max_events: u32,
    /// Sum of the max ages produced by collision events at this
    /// step. Used by the harness to compute
    /// `collision_max_age_mean = sum / count`.
    pub collision_max_age_sum: f64,
}

impl AgeEventCounts {
    pub fn collision_max_age_mean(&self) -> f64 {
        if self.collision_max_events == 0 {
            0.0
        } else {
            self.collision_max_age_sum / self.collision_max_events as f64
        }
    }
}

/// Apply the §4.11 / D3 boundary-event resets to the age field.
///
/// Inputs:
/// - `flags`: cell-centred boundary-flag field (already detected
///   for the current step by the boundary-detection pipeline).
/// - `plate_type`: continental vs oceanic per cell, used to
///   distinguish `Q_arc`-eligible cells (continental adjacent to
///   subducting) from generic boundary cells.
/// - `s`: current `S̃` field, used as the dynamic-classification
///   fallback when `plate_type` is not Voronoï-tracked. A cell
///   with `S̃ > 0.5` is treated as continental for the arc
///   detection, mirroring the `is_continental_thickness`
///   threshold in `init.rs`.
/// - `a`: age field; mutated in place. `current` after advection.
///
/// Returns the per-step event counts for diagnostics.
pub fn apply_age_events(
    flags: &BoundaryFlagField,
    plate_type: &PlateTypeField,
    s: &Field2D,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    a: &mut Field2D,
) -> AgeEventCounts {
    let nx = a.nx();
    let ny = a.ny();
    debug_assert_eq!(flags.nx(), nx);
    debug_assert_eq!(flags.ny(), ny);
    debug_assert_eq!(plate_type.nx(), nx);
    debug_assert_eq!(plate_type.ny(), ny);
    debug_assert_eq!(s.nx(), nx);
    debug_assert_eq!(s.ny(), ny);

    let mut counts = AgeEventCounts::default();

    // Pass 1 — detect arc cells. An "arc" cell is continental AND
    // has at least one subducting neighbour. This mirrors the
    // Q_arc computation in `source_sink::compute_source_sink_terms`
    // (pass 2 there) so our reset semantics matches the volcanic-
    // resurfacing source term that put fresh mass on the cell.
    //
    // We materialise the arc mask in a small Vec so the second
    // pass can write through `a` mutably without an aliasing
    // borrow — `a` is read indirectly inside the arc detection
    // (we look at neighbours of subducting cells), so doing this
    // in a single pass would require reading-and-writing in the
    // same loop on a single buffer.
    let mut is_arc = vec![false; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            // Cell (i, j) is an arc candidate if it is itself
            // continental (per `plate_type` Voronoï tag OR per the
            // S̃ > 0.5 threshold — we use the union to be robust
            // to advection driving plate-type drift).
            let self_is_continental = matches!(plate_type.get(i, j), PlateType::Continental)
                || is_continental_thickness(s.get(i, j));
            if !self_is_continental {
                continue;
            }
            // And has at least one subducting neighbour (4-cell
            // periodic stencil — same as Q_arc).
            let neigh = [
                (idx_x.next(i), j),
                (idx_x.prev(i), j),
                (i, idx_y.next(j)),
                (i, idx_y.prev(j)),
            ];
            if neigh.iter().any(|&(ni, nj)| flags.get(ni, nj).is_subduction()) {
                is_arc[j * nx + i] = true;
            }
        }
    }

    // Pass 2 — apply resets. Order matters when a cell has two
    // applicable rules:
    //   1. Collision overwrites with the max-of-neighbours value
    //      (potentially large).
    //   2. Ridge / arc overwrite with 0.
    // Intersecting cases (a cell flagged both Rift and adjacent to
    // subduction) take the *most-recent* action, which we choose
    // to be ridge/arc → 0 (volcanic / mid-ocean signature wins
    // over distant subduction). We loop in a fixed cell order
    // and write each cell at most once.
    //
    // Special case: collision cells need the OLD age field for the
    // max, because if an earlier cell in the loop has already
    // been reset to 0, taking its post-reset value would lose the
    // protolith age. We snapshot `a` before the pass.
    let a_snapshot = a.clone();

    for j in 0..ny {
        for i in 0..nx {
            let flag = flags.get(i, j);

            // 1. Continental collision — A := max over self + 4
            //    continental neighbours, using the pre-pass
            //    snapshot.
            if matches!(flag, BoundaryFlag::ContinentalCollision) {
                let mut max_age = a_snapshot.get(i, j);
                let neigh = [
                    (idx_x.next(i), j),
                    (idx_x.prev(i), j),
                    (i, idx_y.next(j)),
                    (i, idx_y.prev(j)),
                ];
                for (ni, nj) in neigh {
                    let n_is_continental =
                        matches!(plate_type.get(ni, nj), PlateType::Continental)
                            || is_continental_thickness(s.get(ni, nj));
                    if n_is_continental {
                        let na = a_snapshot.get(ni, nj);
                        if na > max_age {
                            max_age = na;
                        }
                    }
                }
                a.set(i, j, max_age);
                counts.collision_max_events += 1;
                counts.collision_max_age_sum += max_age;
                continue; // do not apply ridge/arc to a collision cell
            }

            // 2. Ridge — Rift flag, A := 0.
            if matches!(flag, BoundaryFlag::Rift) {
                a.set(i, j, 0.0);
                counts.ridge_resets += 1;
                continue;
            }

            // 3. Arc — continental cell adjacent to subducting,
            //    A := 0.
            if is_arc[j * nx + i] {
                a.set(i, j, 0.0);
                counts.arc_resets += 1;
                continue;
            }

            // 4. Quiescent / Subduction / no-event: no action
            //    (advection's `+dt` quiescent growth already
            //    applied for non-subducting cells; subducting
            //    cells will be replenished by advection next
            //    step).
        }
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn periodic(nx: usize, ny: usize) -> (PeriodicIndex, PeriodicIndex) {
        (PeriodicIndex::new(nx), PeriodicIndex::new(ny))
    }

    #[test]
    fn ridge_cell_is_reset_to_zero() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y) = periodic(nx, ny);
        let mut a = Field2D::filled(nx, ny, 3.0);
        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 2, BoundaryFlag::Rift);
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let s = Field2D::filled(nx, ny, 0.2); // oceanic
        let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);
        assert_eq!(a.get(1, 2), 0.0);
        assert_eq!(counts.ridge_resets, 1);
        // Other cells unchanged.
        assert_eq!(a.get(0, 0), 3.0);
        assert_eq!(a.get(2, 1), 3.0);
    }

    #[test]
    fn collision_cell_takes_max_of_continental_neighbours() {
        let nx = 5;
        let ny = 5;
        let (idx_x, idx_y) = periodic(nx, ny);
        let mut a = Field2D::filled(nx, ny, 1.0);
        // Plant non-trivial ages on the collision cell's
        // neighbours: the right neighbour has a much older age
        // and should win the max.
        a.set(2, 2, 4.0); // self
        a.set(3, 2, 9.0); // east neighbour — the protolith
        a.set(1, 2, 3.0);
        a.set(2, 3, 5.0);
        a.set(2, 1, 2.0);

        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(2, 2, BoundaryFlag::ContinentalCollision);
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let s = Field2D::filled(nx, ny, 1.0); // continental everywhere

        let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);
        assert_eq!(a.get(2, 2), 9.0, "collision should pick max neighbour age");
        assert_eq!(counts.collision_max_events, 1);
        assert!((counts.collision_max_age_mean() - 9.0).abs() < 1e-12);
    }

    #[test]
    fn arc_cell_continental_adjacent_to_subducting_resets_to_zero() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y) = periodic(nx, ny);
        let mut a = Field2D::filled(nx, ny, 6.0);

        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(2, 1, BoundaryFlag::OceanicSubduction); // subducting

        // Cell (2, 2) is continental and adjacent to the
        // subducting cell (2, 1) — it is an arc cell.
        let mut plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        plate_type.set(2, 2, PlateType::Continental);
        let mut s = Field2D::filled(nx, ny, 0.2);
        s.set(2, 2, 1.0);

        let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);
        assert_eq!(a.get(2, 2), 0.0, "arc cell should reset");
        assert_eq!(counts.arc_resets, 1);
        // Subducting cell itself: no action (mass goes via S̃
        // recycling). a stays at 6.0 from the initial fill.
        assert_eq!(a.get(2, 1), 6.0);
    }

    #[test]
    fn quiescent_cells_are_unchanged() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y) = periodic(nx, ny);
        let mut a = Field2D::filled(nx, ny, 2.5);
        let flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let s = Field2D::filled(nx, ny, 1.0);
        let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);
        for v in a.data() {
            assert_eq!(*v, 2.5);
        }
        assert_eq!(counts.ridge_resets, 0);
        assert_eq!(counts.arc_resets, 0);
        assert_eq!(counts.collision_max_events, 0);
    }

    #[test]
    fn collision_uses_pre_pass_snapshot_not_partially_updated_field() {
        // Two adjacent collision cells where the loop order would
        // matter if we read the live `a`. With the snapshot
        // semantics, each cell sees the OLD ages of its
        // neighbours regardless of iteration order.
        let nx = 5;
        let ny = 1;
        let (idx_x, idx_y) = periodic(nx, ny);
        let mut a = Field2D::new(nx, ny);
        a.set(0, 0, 10.0);
        a.set(1, 0, 1.0); // collision-1
        a.set(2, 0, 1.0); // collision-2
        a.set(3, 0, 5.0);
        a.set(4, 0, 0.0);

        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 0, BoundaryFlag::ContinentalCollision);
        flags.set(2, 0, BoundaryFlag::ContinentalCollision);
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let s = Field2D::filled(nx, ny, 1.0);

        apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);

        // collision-1 max(self=1, east=1, west=10) = 10
        assert_eq!(a.get(1, 0), 10.0);
        // collision-2 max(self=1, east=5, west=1 — *not* 10
        // because we use the snapshot, not the just-updated value)
        assert_eq!(a.get(2, 0), 5.0);
    }

    #[test]
    fn ridge_and_arc_can_coexist_only_one_action_per_cell() {
        // Edge case: a Rift-flagged cell that is also continental
        // and adjacent to a subducting cell. Ridge takes
        // precedence (matches the spec interpretation: an
        // explicit Rift flag overrides the implicit arc
        // detection). The cell is reset to 0 either way.
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y) = periodic(nx, ny);
        let mut a = Field2D::filled(nx, ny, 4.0);

        let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        flags.set(1, 0, BoundaryFlag::OceanicSubduction);
        flags.set(1, 1, BoundaryFlag::Rift); // continental cell flagged
                                              // both Rift and arc-eligible
        let mut plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        plate_type.set(1, 1, PlateType::Continental);
        let mut s = Field2D::filled(nx, ny, 0.2);
        s.set(1, 1, 1.0);

        let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);
        assert_eq!(a.get(1, 1), 0.0);
        // Counted as ridge (the explicit Rift flag wins the
        // attribution).
        assert_eq!(counts.ridge_resets, 1);
        assert_eq!(counts.arc_resets, 0);
    }
}
