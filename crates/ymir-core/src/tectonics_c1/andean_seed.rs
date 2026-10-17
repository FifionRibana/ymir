//! ADR Finding 178-A1 -- **the Andean seed**: Davis-Suppe's continent-ocean wedge starts on the overriding plate's continental
//! front, within the arc-trench distance of an active ocean-ocean trench.
//!
//! F177 (B2) found that an active margin whose trench sits offshore, inside the continent's own oceanic band, gets no relief:
//! the continent-ocean seed ([`super::boundary_classification::oc_override_seed_mask`]) needs the continent to touch the
//! subducting plate, and option 5's A7 cuts the ocean-ocean seeds. On Earth the arc stands where the slab reaches ~100 km
//! depth, 85–470 km behind the trench (Syracuse & Abers 2006, p. 10), and the Andean shortening is in the overriding plate.
//! Here a continental cell of the overriding plate P becomes a continent-ocean seed when it is on P's continental front facing
//! an active trench (the first continental cells of P a walk from the trench through P's cells reaches) within `gap` cells.
//! Davis-Suppe itself (its formula, reach, coupling, decay) is unchanged; only where the wedge starts changes. Declared in
//! `docs/reports/relief_method/f178_andean_reach/f178_declared.md`.

use std::collections::VecDeque;

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::voronoi::PlateIdField;

use super::boundary_classification::{BoundaryInfo, BoundaryType};
use super::closures::subduction::OoTracker;
use super::kinematics::PlateKinematics;
use super::state::BoolField;

/// The overriding side of the active ocean-ocean trenches: an oceanic convergent cell b of plate P, touching no continental
/// cell of another plate, whose partner (the ocean-ocean contact of a subducting pair with the largest v_rel · n̂, the last of
/// equals in E, W, N, S order — as `apply_oo_subduction_step`) is not out-aged by b: b does not sink by option 5's rule (the
/// older cell sinks; ties: the velocity flag read before A7, then the lower plate id).
#[allow(clippy::too_many_arguments)]
pub fn overriding_trench_cells(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    age: &Field2D,
    boundary: &BoundaryInfo,
    upper_vel: &BoolField,
    kinematics: &PlateKinematics,
    tracker: &OoTracker,
) -> BoolField {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let (idx_x, idx_y) = (PeriodicIndex::new(nx), PeriodicIndex::new(ny));
    let mut out = BoolField::filled(nx, ny, false);
    for j in 0..ny {
        for i in 0..nx {
            if !matches!(boundary.boundary_type.get(i, j), BoundaryType::Convergent) || plate_type.get(i, j) != PlateType::Oceanic {
                continue;
            }
            let pid_c = plate_id.get(i, j);
            let (vx_c, vy_c) = kinematics.velocities[pid_c as usize];
            let mut touches_continent = false;
            let mut best: Option<(usize, usize, f64)> = None;
            for (ni, nj, nx_n, ny_n) in [
                (idx_x.next(i), j, 1.0, 0.0),
                (idx_x.prev(i), j, -1.0, 0.0),
                (i, idx_y.next(j), 0.0, 1.0),
                (i, idx_y.prev(j), 0.0, -1.0),
            ] {
                let pid_n = plate_id.get(ni, nj);
                if pid_n == pid_c {
                    continue;
                }
                if plate_type.get(ni, nj) == PlateType::Continental {
                    touches_continent = true;
                    break;
                }
                let (vx_n, vy_n) = kinematics.velocities[pid_n as usize];
                let dot = (vx_c - vx_n) * nx_n + (vy_c - vy_n) * ny_n;
                if dot > 0.0 && tracker.is_active(pid_c, pid_n) && best.is_none_or(|b| dot >= b.2) {
                    best = Some((ni, nj, dot));
                }
            }
            if touches_continent {
                continue;
            }
            let Some((ni, nj, _)) = best else { continue };
            let pid_n = plate_id.get(ni, nj);
            let (ac, an) = (age.get(i, j), age.get(ni, nj));
            let sinks = if ac != an {
                ac > an
            } else {
                match (upper_vel.get(i, j), upper_vel.get(ni, nj)) {
                    (false, true) => true,
                    (true, false) => false,
                    _ => pid_c < pid_n,
                }
            };
            if !sinks {
                out.set(i, j, true);
            }
        }
    }
    out
}

/// ADR Finding 178-A1 -- the Andean seeds: P's continental cells first reached by a walk (4-connected, periodic) from P's
/// overriding trench cells through P's own cells, at a walk depth ≤ `gap` (the trench cell at 0). `gap = 0` seeds nothing.
#[allow(clippy::too_many_arguments)]
pub fn andean_seeds(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    age: &Field2D,
    boundary: &BoundaryInfo,
    upper_vel: &BoolField,
    kinematics: &PlateKinematics,
    tracker: &OoTracker,
    gap: usize,
) -> BoolField {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let (idx_x, idx_y) = (PeriodicIndex::new(nx), PeriodicIndex::new(ny));
    let mut seeds = BoolField::filled(nx, ny, false);
    if gap == 0 {
        return seeds;
    }
    let trench = overriding_trench_cells(plate_id, plate_type, age, boundary, upper_vel, kinematics, tracker);
    let mut depth = vec![usize::MAX; nx * ny];
    let mut q = VecDeque::new();
    for j in 0..ny {
        for i in 0..nx {
            if trench.get(i, j) {
                depth[j * nx + i] = 0;
                q.push_back((i, j));
            }
        }
    }
    while let Some((i, j)) = q.pop_front() {
        let d = depth[j * nx + i];
        let p = plate_id.get(i, j);
        for (ni, nj) in [(idx_x.next(i), j), (idx_x.prev(i), j), (i, idx_y.next(j)), (i, idx_y.prev(j))] {
            let k = nj * nx + ni;
            if plate_id.get(ni, nj) != p || depth[k] != usize::MAX {
                continue;
            }
            depth[k] = d + 1;
            if plate_type.get(ni, nj) == PlateType::Continental {
                // the front: a continental cell stops the walk
                if d + 1 <= gap {
                    seeds.set(ni, nj, true);
                }
            } else if d + 1 < gap {
                q.push_back((ni, nj));
            }
        }
    }
    seeds
}
