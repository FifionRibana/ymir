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

/// ADR Finding 179-B2 -- **the Andean seed at the arc**: for each active trench on which P overrides
/// ([`overriding_trench_cells`], grouped into 8-connected segments of one plate), a walk (4-connected, periodic) from the
/// segment through P's cells of any type; d_front = the depth of the first continental cell of P reached; the seeds are P's
/// continental cells at exactly depth d = max(`d_arc`, d_front), if d ≤ `gap`. No reading of the sea level. `gap = 0` seeds
/// nothing. Declared in `docs/reports/relief_method/f179_pumping/f179_declared.md` (d_arc = 4: 180 km, the lower bound of the
/// middle 50 % of arc-trench distances, Syracuse & Abers 2006 p. 10, ÷ 7.5 ÷ 6.25 = 3.84).
#[allow(clippy::too_many_arguments)]
pub fn andean_seeds_arc(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    age: &Field2D,
    boundary: &BoundaryInfo,
    upper_vel: &BoolField,
    kinematics: &PlateKinematics,
    tracker: &OoTracker,
    d_arc: usize,
    gap: usize,
) -> BoolField {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let (idx_x, idx_y) = (PeriodicIndex::new(nx), PeriodicIndex::new(ny));
    let mut seeds = BoolField::filled(nx, ny, false);
    if gap == 0 {
        return seeds;
    }
    let trench = overriding_trench_cells(plate_id, plate_type, age, boundary, upper_vel, kinematics, tracker);
    let n = nx * ny;
    let mut seen = vec![false; n];
    for s0 in 0..n {
        let (i0, j0) = (s0 % nx, s0 / nx);
        if !trench.get(i0, j0) || seen[s0] {
            continue;
        }
        // the segment: 8-connected trench cells of the same plate
        let p = plate_id.get(i0, j0);
        let mut seg = vec![s0];
        seen[s0] = true;
        let mut k = 0;
        while k < seg.len() {
            let (i, j) = (seg[k] % nx, seg[k] / nx);
            for (ni, nj) in [
                (idx_x.next(i), j),
                (idx_x.prev(i), j),
                (i, idx_y.next(j)),
                (i, idx_y.prev(j)),
                (idx_x.next(i), idx_y.next(j)),
                (idx_x.next(i), idx_y.prev(j)),
                (idx_x.prev(i), idx_y.next(j)),
                (idx_x.prev(i), idx_y.prev(j)),
            ] {
                let q = nj * nx + ni;
                if !seen[q] && trench.get(ni, nj) && plate_id.get(ni, nj) == p {
                    seen[q] = true;
                    seg.push(q);
                }
            }
            k += 1;
        }
        // the walk through P's cells from the segment, to depth `gap`
        let mut depth = vec![usize::MAX; n];
        let mut q = VecDeque::new();
        for &c in &seg {
            depth[c] = 0;
            q.push_back(c);
        }
        let mut d_front = usize::MAX;
        let mut by_depth: Vec<Vec<usize>> = vec![Vec::new(); gap + 1];
        while let Some(c) = q.pop_front() {
            let d = depth[c];
            let (i, j) = (c % nx, c / nx);
            if d > 0 && plate_type.get(i, j) == PlateType::Continental {
                d_front = d_front.min(d);
                by_depth[d].push(c);
            }
            if d == gap {
                continue;
            }
            for (ni, nj) in [(idx_x.next(i), j), (idx_x.prev(i), j), (i, idx_y.next(j)), (i, idx_y.prev(j))] {
                let k2 = nj * nx + ni;
                if depth[k2] == usize::MAX && plate_id.get(ni, nj) == p {
                    depth[k2] = d + 1;
                    q.push_back(k2);
                }
            }
        }
        if d_front == usize::MAX {
            continue;
        }
        let d = d_arc.max(d_front);
        if d <= gap {
            for &c in &by_depth[d] {
                seeds.set(c % nx, c / nx, true);
            }
        }
    }
    seeds
}

/// ADR Finding 180 -- the overriding trench cells with their sinking partner: per cell, `Some((partner cell, v_rel · n̂))` when
/// it is an overriding-side cell of an active ocean-ocean trench ([`overriding_trench_cells`]'s rule), else `None`.
#[allow(clippy::too_many_arguments)]
pub fn overriding_trench_partners(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    age: &Field2D,
    boundary: &BoundaryInfo,
    upper_vel: &BoolField,
    kinematics: &PlateKinematics,
    tracker: &OoTracker,
) -> Vec<Option<(usize, f64)>> {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let (idx_x, idx_y) = (PeriodicIndex::new(nx), PeriodicIndex::new(ny));
    let mut out = vec![None; nx * ny];
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
            let Some((ni, nj, dot)) = best else { continue };
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
                out[j * nx + i] = Some((nj * nx + ni, dot));
            }
        }
    }
    out
}

/// ADR Finding 180-A1 -- a set of Andean seeds with, for each, the trench cell that feeds it.
#[derive(Clone, Debug)]
pub struct AndeanSeedSet {
    pub mask: BoolField,
    /// (seed cell, trench cell), row-major indices.
    pub feeds: Vec<(usize, usize)>,
}

/// ADR Finding 180-A1 -- **the Andean seed behind the active coasts only**. For every coast cell of a plate P (continental, with
/// an oceanic 4-neighbour, no 4-neighbour of another plate) A1's walk goes through P's oceanic cells (breadth-first, E, W, N, S)
/// to its first depth with boundary cells; the coast is active by the majority rule when strictly more of them are overriding
/// trench cells than not. Its trench is the first overriding cell of that depth; d_front its depth. The seed is P's continental
/// cell on the walk's path at depth d = max(`d_arc`, d_front) from the trench, if d ≤ `gap`: the coast cell itself, or, when
/// d_front < d_arc, the cell reached by continuing the path's last step straight inland for d_arc − d_front cells (if
/// continental of P). Declared in `docs/reports/relief_method/f180_andean_mass/f180_declared.md`.
#[allow(clippy::too_many_arguments)]
pub fn andean_seeds_active(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    age: &Field2D,
    boundary: &BoundaryInfo,
    upper_vel: &BoolField,
    kinematics: &PlateKinematics,
    tracker: &OoTracker,
    d_arc: usize,
    gap: usize,
) -> AndeanSeedSet {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let n = nx * ny;
    let mut mask = BoolField::filled(nx, ny, false);
    let mut feeds = Vec::new();
    if gap == 0 {
        return AndeanSeedSet { mask, feeds };
    }
    let partners = overriding_trench_partners(plate_id, plate_type, age, boundary, upper_vel, kinematics, tracker);
    let pid = |k: usize| plate_id.get(k % nx, k / nx);
    let typ = |k: usize| plate_type.get(k % nx, k / nx);
    let nb4 = |k: usize| -> [usize; 4] {
        let (x, y) = (k % nx, k / nx);
        [y * nx + (x + 1) % nx, y * nx + (x + nx - 1) % nx, ((y + 1) % ny) * nx + x, ((y + ny - 1) % ny) * nx + x]
    };
    let is_bnd = |k: usize| nb4(k).iter().any(|&j| pid(j) != pid(k));
    let mut stamp = vec![usize::MAX; n];
    let mut pred = vec![usize::MAX; n];
    for k in 0..n {
        if typ(k) != PlateType::Continental || !nb4(k).iter().any(|&j| typ(j) == PlateType::Oceanic) || is_bnd(k) {
            continue;
        }
        let p = pid(k);
        stamp[k] = k;
        let mut frontier: Vec<usize> = Vec::new();
        for j in nb4(k) {
            if typ(j) == PlateType::Oceanic && pid(j) == p && stamp[j] != k {
                stamp[j] = k;
                pred[j] = k;
                frontier.push(j);
            }
        }
        let mut depth = 1usize;
        let mut found: Option<(usize, usize)> = None; // (trench, depth)
        while !frontier.is_empty() {
            let bnds: Vec<usize> = frontier.iter().copied().filter(|&q| is_bnd(q)).collect();
            if !bnds.is_empty() {
                let act: Vec<usize> = bnds.iter().copied().filter(|&b| partners[b].is_some()).collect();
                if act.len() * 2 > bnds.len() {
                    found = Some((act[0], depth));
                }
                break;
            }
            if depth >= gap {
                break;
            }
            let mut next = Vec::new();
            for &q in &frontier {
                for j in nb4(q) {
                    if stamp[j] != k && typ(j) == PlateType::Oceanic && pid(j) == p {
                        stamp[j] = k;
                        pred[j] = q;
                        next.push(j);
                    }
                }
            }
            frontier = next;
            depth += 1;
        }
        let Some((trench, d_front)) = found else { continue };
        let d = d_arc.max(d_front);
        if d > gap {
            continue;
        }
        let seed = if d_front >= d_arc {
            Some(k)
        } else {
            // the path's cell next to the coast (depth 1), then straight inland
            let mut c = trench;
            while pred[c] != k {
                c = pred[c];
            }
            let (cx, cy) = ((c % nx) as i64, (c / nx) as i64);
            let (kx, ky) = ((k % nx) as i64, (k / nx) as i64);
            let wrap = |v: i64, m: usize| -> i64 {
                let m = m as i64;
                let r = v.rem_euclid(m);
                if r > m / 2 { r - m } else { r }
            };
            let (dx, dy) = (wrap(kx - cx, nx), wrap(ky - cy, ny));
            let steps = (d_arc - d_front) as i64;
            let (tx, ty) = ((kx + dx * steps).rem_euclid(nx as i64) as usize, (ky + dy * steps).rem_euclid(ny as i64) as usize);
            let t = ty * nx + tx;
            (typ(t) == PlateType::Continental && pid(t) == p).then_some(t)
        };
        if let Some(sd) = seed {
            if !mask.get(sd % nx, sd / nx) {
                mask.set(sd % nx, sd / nx, true);
            }
            feeds.push((sd, trench));
        }
    }
    // one feed per seed: a seed reached from several coasts keeps the first (the lowest coast index)
    feeds.sort_by_key(|f| f.0);
    feeds.dedup_by_key(|f| f.0);
    AndeanSeedSet { mask, feeds }
}

/// ADR Finding 180-A2 -- the feeding trench of seeds that carry none (F179's arc seeds under the bound): the nearest overriding
/// trench cell of the same plate, by a breadth-first walk through the plate's cells (4-connected; ties: the lower trench index).
pub fn andean_feeds_nearest(plate_id: &PlateIdField, seeds: &BoolField, partners: &[Option<(usize, f64)>]) -> Vec<(usize, usize)> {
    let (nx, ny) = (plate_id.nx(), plate_id.ny());
    let n = nx * ny;
    let mut lab = vec![usize::MAX; n];
    let mut q = VecDeque::new();
    for k in 0..n {
        if partners[k].is_some() {
            lab[k] = k;
            q.push_back(k);
        }
    }
    while let Some(k) = q.pop_front() {
        let (x, y) = (k % nx, k / nx);
        for j in [y * nx + (x + 1) % nx, y * nx + (x + nx - 1) % nx, ((y + 1) % ny) * nx + x, ((y + ny - 1) % ny) * nx + x] {
            if lab[j] == usize::MAX && plate_id.get(j % nx, j / nx) == plate_id.get(x, y) {
                lab[j] = lab[k];
                q.push_back(j);
            }
        }
    }
    (0..n).filter(|&k| seeds.get(k % nx, k / nx) && lab[k] != usize::MAX).map(|k| (k, lab[k])).collect()
}
