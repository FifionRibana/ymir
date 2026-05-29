//! Subduction event — oceanic mass consumption + arc volcanism +
//! floor-triggered `plate_id` reassignment.
//!
//! ## Per-cell algorithm
//!
//! For each cell `c` where `boundary_type[c] == Convergent` AND
//! `plate_type[c] == Oceanic`:
//!
//! 1. **Find the continental neighbour with the largest positive
//!    `v_rel · n̂`.** Iterates the 4-connected neighbours; for each
//!    neighbour `n` with `plate_id[n] != plate_id[c]` and
//!    `plate_type[n] == Continental`, compute the dot product of
//!    the per-plate relative velocity onto the outward normal `n̂`
//!    pointing from `c` toward `n`. The neighbour producing the
//!    largest positive value (most-convergent direction) is the
//!    "upper" continental side of this subduction edge. If no
//!    continental neighbour exists with positive convergence, the
//!    cell is skipped (oceanic-oceanic convergent boundary, not
//!    handled by this closure).
//!
//! 2. **Consume `Δs = consumption_rate × convergence × dt`,
//!    clamped at `s_before`** so `S̃` cannot go negative. The
//!    clamping is a defensive guard — at default parameters and
//!    Phase 1.1 kinematics, a single step's `Δs` is far below the
//!    initial oceanic baseline.
//!
//! 3. **Compute `arc_mass = Δs × arc_efficiency`.** The complement
//!    (`Δs × (1 − arc_efficiency)`) is the fraction lost to the
//!    deeper mantle, out of model.
//!
//! 4. **If `s[c] < plate_id_reassign_threshold` after step 2**,
//!    reassign `plate_id[c]` to the continental neighbour's plate
//!    id and `plate_type[c]` to `Continental`. The cell's
//!    (small) remaining `S̃` value stays in place — its mass is
//!    absorbed by the continental plate via the `plate_type`
//!    promotion, no S̃ reset needed for mass conservation.
//!
//! 5. **Distribute `arc_mass` via BFS** up to `arc_distance` cells
//!    from `c`. The BFS visits 4-connected neighbours layer by
//!    layer; cells reached at depth ≥ 1 with `plate_type ==
//!    Continental` collect arc mass at `per_cell = arc_mass /
//!    n_continental_cells_reached`. If the BFS finds **zero**
//!    continental cells (e.g., consuming cell isolated in an ocean
//!    or the continental neighbour was just promoted away in this
//!    step), the arc mass is **lost** — surfaced via
//!    `SubductionStats.arc_mass_distributed < total_mass_consumed
//!    × arc_efficiency`.
//!
//! Iteration order is row-major `(j, i)` — deterministic given the
//! same input. Mutations within a step propagate to later cells in
//! the same step (`plate_type` promotion may affect a later cell's
//! BFS continental count), which is a conscious choice to keep the
//! algorithm single-pass and simple.
//!
//! ## Architectural concerns surfaced Stage E1 (W7)
//!
//! 1. **Arc mass loss when BFS finds 0 continental cells.** Logged
//!    via the gap `arc_mass_distributed < total_mass_consumed ×
//!    arc_efficiency`. Mass conservation diagnostic in Stage E4
//!    accounts for this by tracking the two stats independently.
//!    Not a bug — graceful handling per the user spec.
//! 2. **Ambiguous boundaries excluded.** Triple junctions (cells
//!    with both Convergent and Divergent neighbours) are skipped
//!    by the `BoundaryType::Convergent` filter. They constitute a
//!    small fraction of cells (<5 % at typical Phase 1.1 init) and
//!    excluding them parallels the Davis-Suppe convention. Revisit
//!    if Stage A event-count diagnostic shows insufficient activity.
//! 3. **Single-pass mutation propagation.** A cell reassigned to
//!    Continental in row `j` can be a continental BFS target for a
//!    later row's subduction. Acceptable for this iteration; if
//!    Stage A shows determinism / ordering artefacts, switch to a
//!    collect-then-apply two-pass pattern.

use std::collections::{HashSet, VecDeque};

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::voronoi::PlateIdField;

use crate::tectonics_c1::boundary_classification::{BoundaryInfo, BoundaryType};
use crate::tectonics_c1::kinematics::PlateKinematics;

use super::params::SubductionParams;

/// Per-step subduction diagnostics. Returned by
/// `apply_subduction_step` for the Stage E4 mass-balance test.
///
/// **Mass conservation invariant** (test-only): the per-step
/// difference of `Σ s` over the grid equals
/// `total_mass_consumed − arc_mass_distributed`, modulo the
/// `arc_mass_lost = total_mass_consumed × arc_efficiency −
/// arc_mass_distributed` fraction that BFS couldn't place on a
/// continental cell.
#[derive(Clone, Copy, Debug, Default)]
pub struct SubductionStats {
    /// Number of distinct oceanic boundary cells that received a
    /// consumption term this step.
    pub cells_consumed: usize,
    /// Sum of `Δs` over all consuming cells. The amount actually
    /// removed from oceanic `S̃` (post-clamp).
    pub total_mass_consumed: f64,
    /// Sum of arc volcanism mass actually placed on continental
    /// cells via BFS distribution. Less than
    /// `total_mass_consumed × arc_efficiency` when BFS finds 0
    /// continental cells (arc-mass-lost case).
    pub arc_mass_distributed: f64,
    /// Number of cells reassigned from Oceanic to Continental
    /// (floor-triggered, `S̃ < plate_id_reassign_threshold` after
    /// consumption). First Track D mutation of `plate_id` /
    /// `plate_type`.
    pub plate_ids_reassigned: usize,
}

/// Apply one forward-Euler step of the subduction closure to the
/// state.
///
/// Mutates `s` (consumption + arc), `plate_id` (floor-triggered
/// reassignment), and `plate_type` (Oceanic → Continental on
/// reassignment). All other state passed read-only.
///
/// Returns [`SubductionStats`] for mass-balance accounting (used by
/// Stage E4's per-step diagnostic test).
///
/// Returns immediately with `SubductionStats::default()` when
/// `params.enabled == false` — bit-identical no-op (W4 closure-
/// isolation discipline).
pub fn apply_subduction_step(
    s: &mut Field2D,
    plate_id: &mut PlateIdField,
    plate_type: &mut PlateTypeField,
    boundary_info: &BoundaryInfo,
    kinematics: &PlateKinematics,
    params: &SubductionParams,
    dt: f64,
) -> SubductionStats {
    if !params.enabled {
        return SubductionStats::default();
    }

    let nx = s.nx();
    let ny = s.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    // 4-neighbour offsets and outward normals matching the
    // convention in `boundary_classification::classify_boundaries`
    // (normal points from central cell toward neighbour, so
    // `v_rel · n̂ > 0` ⇔ convergent).
    let neighbours: [(i32, i32, f64, f64); 4] = [
        (1, 0, 1.0, 0.0),
        (-1, 0, -1.0, 0.0),
        (0, 1, 0.0, 1.0),
        (0, -1, 0.0, -1.0),
    ];

    let mut stats = SubductionStats::default();

    for j in 0..ny {
        for i in 0..nx {
            // Filter 1: cell must be on a Convergent boundary.
            if !matches!(
                boundary_info.boundary_type.get(i, j),
                BoundaryType::Convergent
            ) {
                continue;
            }
            // Filter 2: cell must be Oceanic (subducting side).
            if plate_type.get(i, j) != PlateType::Oceanic {
                continue;
            }

            let pid_c = plate_id.get(i, j);
            let (vx_c, vy_c) = kinematics.velocities[pid_c as usize];

            // Pick the continental neighbour with the largest
            // positive convergence dot product.
            let mut best_continental_pid: Option<u16> = None;
            let mut best_convergence = 0.0_f64;

            for &(di, dj, nx_norm, ny_norm) in neighbours.iter() {
                let ni = if di > 0 {
                    idx_x.next(i)
                } else if di < 0 {
                    idx_x.prev(i)
                } else {
                    i
                };
                let nj = if dj > 0 {
                    idx_y.next(j)
                } else if dj < 0 {
                    idx_y.prev(j)
                } else {
                    j
                };
                let pid_n = plate_id.get(ni, nj);
                if pid_n == pid_c {
                    continue;
                }
                if plate_type.get(ni, nj) != PlateType::Continental {
                    continue;
                }
                let (vx_n, vy_n) = kinematics.velocities[pid_n as usize];
                let vrel_x = vx_c - vx_n;
                let vrel_y = vy_c - vy_n;
                let dot = vrel_x * nx_norm + vrel_y * ny_norm;
                if dot > best_convergence {
                    best_convergence = dot;
                    best_continental_pid = Some(pid_n);
                }
            }

            let Some(continental_pid) = best_continental_pid else {
                // No continental neighbour with positive convergence
                // — oceanic-oceanic convergent boundary, not handled
                // by this closure.
                continue;
            };

            // Consumption — clamped at s_before so S̃ can't go
            // negative even at pathological consumption_rate.
            let s_before = s.get(i, j);
            let delta_s_proposed = params.consumption_rate * best_convergence * dt;
            let delta_s = delta_s_proposed.min(s_before);
            let s_after = s_before - delta_s;
            s.set(i, j, s_after);

            stats.cells_consumed += 1;
            stats.total_mass_consumed += delta_s;

            let arc_mass = delta_s * params.arc_efficiency;

            // Floor-triggered plate_id reassignment.
            if s_after < params.plate_id_reassign_threshold {
                plate_id.set(i, j, continental_pid);
                plate_type.set(i, j, PlateType::Continental);
                stats.plate_ids_reassigned += 1;
            }

            // Distribute arc mass via BFS up to arc_distance cells.
            let distributed = distribute_arc_mass(
                s,
                plate_type,
                i,
                j,
                arc_mass,
                params.arc_distance,
                &idx_x,
                &idx_y,
            );
            stats.arc_mass_distributed += distributed;
        }
    }

    stats
}

/// BFS-based arc-volcanism distribution.
///
/// From the consuming oceanic cell `(origin_i, origin_j)`, walk
/// 4-connected neighbours up to `arc_distance` cells. Among the
/// cells visited at depth ≥ 1, those with `plate_type ==
/// Continental` collect `arc_mass / n_continental_reached` each.
///
/// Returns the total mass actually placed (`= arc_mass` when at
/// least one continental cell was reached, `0` otherwise — the
/// arc-mass-lost case).
fn distribute_arc_mass(
    s: &mut Field2D,
    plate_type: &PlateTypeField,
    origin_i: usize,
    origin_j: usize,
    arc_mass: f64,
    arc_distance: usize,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
) -> f64 {
    if arc_distance == 0 || arc_mass <= 0.0 {
        return 0.0;
    }

    let mut visited: HashSet<(usize, usize)> = HashSet::new();
    let mut queue: VecDeque<(usize, usize, usize)> = VecDeque::new();
    queue.push_back((origin_i, origin_j, 0));
    visited.insert((origin_i, origin_j));

    let mut continental_cells: Vec<(usize, usize)> = Vec::new();

    while let Some((i, j, depth)) = queue.pop_front() {
        if depth > 0 && plate_type.get(i, j) == PlateType::Continental {
            continental_cells.push((i, j));
        }
        if depth < arc_distance {
            let neighbours = [
                (idx_x.next(i), j),
                (idx_x.prev(i), j),
                (i, idx_y.next(j)),
                (i, idx_y.prev(j)),
            ];
            for (ni, nj) in neighbours {
                if visited.insert((ni, nj)) {
                    queue.push_back((ni, nj, depth + 1));
                }
            }
        }
    }

    if continental_cells.is_empty() {
        return 0.0;
    }

    let per_cell = arc_mass / continental_cells.len() as f64;
    for (ci, cj) in &continental_cells {
        let new_s = s.get(*ci, *cj) + per_cell;
        s.set(*ci, *cj, new_s);
    }
    arc_mass
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::tectonics_c1::boundary_classification::classify_boundaries;
    use crate::tectonics_c1::kinematics::PlateKinematics;

    /// Two-plate east-west fixture used by most of the unit tests.
    ///
    /// Plate 0 occupies `i < nx/2` with velocity `(+v, 0)`.
    /// Plate 1 occupies `i >= nx/2` with velocity `(-v, 0)`.
    /// `plate_type_left` and `plate_type_right` are the per-plate
    /// types assigned to all cells of that plate.
    ///
    /// The Phase 1.1 baseline S̃ values are applied per plate type:
    /// continental cells start at `1.0`, oceanic at `0.2`.
    ///
    /// Returns `(s, plate_id, plate_type, boundary_info, kinematics)`
    /// ready to feed `apply_subduction_step`.
    fn two_plate_east_west_fixture(
        nx: usize,
        ny: usize,
        plate_type_left: PlateType,
        plate_type_right: PlateType,
        v_left: f64,
        v_right: f64,
    ) -> (
        Field2D,
        PlateIdField,
        PlateTypeField,
        BoundaryInfo,
        PlateKinematics,
    ) {
        let mut s = Field2D::new(nx, ny);
        let mut plate_id = PlateIdField::new(nx, ny);
        let mut plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        for j in 0..ny {
            for i in 0..nx {
                let (pid, pt) = if i < nx / 2 {
                    (0_u16, plate_type_left)
                } else {
                    (1_u16, plate_type_right)
                };
                plate_id.set(i, j, pid);
                plate_type.set(i, j, pt);
                let s_init = match pt {
                    PlateType::Continental => 1.0,
                    PlateType::Oceanic => 0.2,
                };
                s.set(i, j, s_init);
            }
        }
        let kinematics = PlateKinematics {
            velocities: vec![(v_left, 0.0), (v_right, 0.0)],
        };
        let boundary_info = classify_boundaries(&plate_id, &kinematics);
        (s, plate_id, plate_type, boundary_info, kinematics)
    }

    #[test]
    fn subduction_consumes_oceanic_mass_at_convergent_oceanic_continental() {
        // Plate 0 = Continental, Plate 1 = Oceanic. They converge
        // at the i=3/i=4 boundary (8-column grid).
        let nx = 8;
        let ny = 4;
        let (mut s, mut plate_id, mut plate_type, boundary_info, kinematics) =
            two_plate_east_west_fixture(
                nx,
                ny,
                PlateType::Continental,
                PlateType::Oceanic,
                0.01,
                -0.01,
            );
        let params = SubductionParams::default();
        let dt = 0.69;

        let s_before_oceanic = s.get(4, 0);

        let stats = apply_subduction_step(
            &mut s,
            &mut plate_id,
            &mut plate_type,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        // The first column of plate 1 (i = 4) is on the convergent
        // boundary AND oceanic — it must have been consumed.
        let s_after_oceanic = s.get(4, 0);
        assert!(
            s_after_oceanic < s_before_oceanic,
            "oceanic boundary cell should be consumed: before {s_before_oceanic}, after {s_after_oceanic}"
        );

        // Stats sanity: at least one cell consumed, some mass moved.
        assert!(
            stats.cells_consumed >= ny,
            "expected at least {ny} cells consumed (one per row), got {}",
            stats.cells_consumed
        );
        assert!(
            stats.total_mass_consumed > 0.0,
            "total_mass_consumed must be positive when consumption fires"
        );
        assert!(
            stats.arc_mass_distributed > 0.0,
            "arc_mass_distributed must be positive when continental neighbours are reachable"
        );
    }

    #[test]
    fn subduction_arc_distribution_targets_continental_neighbors() {
        // Same fixture as test 1 — the cells of plate 0 (continental)
        // adjacent to the convergent boundary must receive arc mass.
        let nx = 8;
        let ny = 4;
        let (mut s, mut plate_id, mut plate_type, boundary_info, kinematics) =
            two_plate_east_west_fixture(
                nx,
                ny,
                PlateType::Continental,
                PlateType::Oceanic,
                0.01,
                -0.01,
            );
        let params = SubductionParams::default();
        let dt = 0.69;

        // Snapshot continental cells' initial S̃.
        let s_before: Vec<f64> = (0..nx / 2)
            .map(|i| s.get(i, 0))
            .collect();

        let _stats = apply_subduction_step(
            &mut s,
            &mut plate_id,
            &mut plate_type,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        // The cell immediately inside plate 0 (i = 3) should have
        // received arc mass — within BFS distance 1 of the
        // convergent oceanic cell at (4, 0).
        let s_after_adjacent = s.get(3, 0);
        let s_before_adjacent = s_before[3];
        assert!(
            s_after_adjacent > s_before_adjacent,
            "continental cell adjacent to subduction must receive arc mass: before {s_before_adjacent}, after {s_after_adjacent}"
        );
    }

    #[test]
    fn subduction_no_op_at_oceanic_oceanic_boundary() {
        // Both plates Oceanic — convergent boundary, but no
        // continental neighbour means subduction can't fire (no
        // arc destination, no upper plate).
        let nx = 8;
        let ny = 4;
        let (mut s, mut plate_id, mut plate_type, boundary_info, kinematics) =
            two_plate_east_west_fixture(
                nx,
                ny,
                PlateType::Oceanic,
                PlateType::Oceanic,
                0.01,
                -0.01,
            );
        let params = SubductionParams::default();
        let dt = 0.69;

        let s_before: Vec<f64> = s.data().to_vec();
        let plate_id_before: Vec<u16> = plate_id.data().to_vec();
        let plate_type_before: Vec<PlateType> = plate_type.data().to_vec();

        let stats = apply_subduction_step(
            &mut s,
            &mut plate_id,
            &mut plate_type,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        assert_eq!(stats.cells_consumed, 0);
        assert_eq!(stats.total_mass_consumed, 0.0);
        assert_eq!(stats.arc_mass_distributed, 0.0);
        assert_eq!(stats.plate_ids_reassigned, 0);
        assert_eq!(s.data(), s_before.as_slice());
        assert_eq!(plate_id.data(), plate_id_before.as_slice());
        for k in 0..plate_type.data().len() {
            assert_eq!(plate_type.data()[k], plate_type_before[k]);
        }
    }

    #[test]
    fn subduction_no_op_at_continental_continental_boundary() {
        // Both plates Continental — convergent boundary, but
        // subduction filter (plate_type == Oceanic) rejects every
        // cell.
        let nx = 8;
        let ny = 4;
        let (mut s, mut plate_id, mut plate_type, boundary_info, kinematics) =
            two_plate_east_west_fixture(
                nx,
                ny,
                PlateType::Continental,
                PlateType::Continental,
                0.01,
                -0.01,
            );
        let params = SubductionParams::default();
        let dt = 0.69;

        let s_before: Vec<f64> = s.data().to_vec();

        let stats = apply_subduction_step(
            &mut s,
            &mut plate_id,
            &mut plate_type,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        assert_eq!(stats.cells_consumed, 0);
        assert_eq!(stats.total_mass_consumed, 0.0);
        assert_eq!(stats.arc_mass_distributed, 0.0);
        assert_eq!(stats.plate_ids_reassigned, 0);
        assert_eq!(s.data(), s_before.as_slice());
    }

    #[test]
    fn subduction_disabled_no_op() {
        // params.enabled = false — bit-identical no-op on all
        // mutable state, regardless of fixture activity.
        let nx = 8;
        let ny = 4;
        let (mut s, mut plate_id, mut plate_type, boundary_info, kinematics) =
            two_plate_east_west_fixture(
                nx,
                ny,
                PlateType::Continental,
                PlateType::Oceanic,
                0.01,
                -0.01,
            );
        let params = SubductionParams {
            enabled: false,
            ..SubductionParams::default()
        };
        let dt = 0.69;

        let s_before: Vec<f64> = s.data().to_vec();
        let plate_id_before: Vec<u16> = plate_id.data().to_vec();
        let plate_type_before: Vec<PlateType> = plate_type.data().to_vec();

        let stats = apply_subduction_step(
            &mut s,
            &mut plate_id,
            &mut plate_type,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        assert_eq!(stats.cells_consumed, 0);
        assert_eq!(stats.total_mass_consumed, 0.0);
        assert_eq!(stats.arc_mass_distributed, 0.0);
        assert_eq!(stats.plate_ids_reassigned, 0);
        assert_eq!(s.data(), s_before.as_slice());
        assert_eq!(plate_id.data(), plate_id_before.as_slice());
        for k in 0..plate_type.data().len() {
            assert_eq!(plate_type.data()[k], plate_type_before[k]);
        }
    }

    #[test]
    fn subduction_plate_id_reassignment_below_floor() {
        // Pathologically high consumption rate drives the oceanic
        // boundary cell's S̃ from 0.2 down to below the floor (0.05)
        // in a single step. The cell must be reassigned to the
        // continental plate.
        //
        // At v_rel · n̂ = 0.02 and dt = 1.0, default
        // consumption_rate = 0.5 would give Δs = 0.01 — well above
        // the floor. To force reassignment, override
        // consumption_rate to push Δs above (s_before − floor) =
        // 0.15. With Δs = consumption_rate × 0.02 × 1.0 > 0.15,
        // consumption_rate > 7.5 suffices. Use 50.0 for clear
        // margin.
        let nx = 8;
        let ny = 4;
        let (mut s, mut plate_id, mut plate_type, boundary_info, kinematics) =
            two_plate_east_west_fixture(
                nx,
                ny,
                PlateType::Continental,
                PlateType::Oceanic,
                0.01,
                -0.01,
            );
        let params = SubductionParams {
            consumption_rate: 50.0,
            ..SubductionParams::default()
        };
        let dt = 1.0;

        let stats = apply_subduction_step(
            &mut s,
            &mut plate_id,
            &mut plate_type,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        assert!(
            stats.plate_ids_reassigned >= ny,
            "expected ≥ {ny} reassignments (one per row), got {}",
            stats.plate_ids_reassigned
        );

        // The previously-oceanic cell at (4, 0) is now part of
        // plate 0 (continental) and typed Continental.
        assert_eq!(
            plate_id.get(4, 0),
            0,
            "reassigned cell must adopt the adjacent continental plate's id"
        );
        assert_eq!(
            plate_type.get(4, 0),
            PlateType::Continental,
            "reassigned cell must become Continental"
        );
        // Its S̃ must be below the floor (post-clamp).
        assert!(
            s.get(4, 0) < params.plate_id_reassign_threshold,
            "reassigned cell S̃ ({}) must be below floor ({})",
            s.get(4, 0),
            params.plate_id_reassign_threshold
        );
    }
}
