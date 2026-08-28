//! Step 13 Phase 1 — shared inter-plate-boundary distance utility.
//!
//! Computes BFS distance from every cell to the nearest inter-plate
//! boundary on the periodic torus, with **Chebyshev 8-connectivity**,
//! plus the plate id of the across-boundary neighbour propagated from
//! the BFS seed. Used by:
//!
//! - [`super::super::init::InitMode::Uniform`] — S̃ initialisation
//!   blending across plate boundaries (Step 8.6 Phase 8a).
//! - [`super::super::plate_kinematic::field::build`] — drift velocity
//!   field smoothing across plate boundaries (Step 11).
//! - `init::radial_profile` — continental margin gradient (Step 13
//!   Phase 2, this issue).
//!
//! ## Why cratonic uses a different BFS
//!
//! [`super::super::cratonic::factor`] computes a *different* distance
//! field for a *different* physical purpose, with a *different*
//! algorithm. The Step 13 issue (D4) initially proposed unifying both
//! under one utility — that proved incorrect during Phase 1 review.
//! Both algorithms remain in their respective modules; only the
//! genuinely shared one (the three callers above) is extracted here.
//!
//! | Property        | This utility                             | `cratonic::factor`                                                    |
//! |-----------------|------------------------------------------|------------------------------------------------------------------------|
//! | Connectivity    | Chebyshev (8-neighbour, periodic)        | Manhattan (4-neighbour, periodic)                                      |
//! | Sources         | Cells whose 4-NESW neighbour is on a different plate (= boundary cells) | Cells of every non-retained plate (oceanic + small-continental — entire plate interiors) |
//! | Distance dtype  | `f64`, sentinel `INFINITY`               | `u32`, sentinel `u32::MAX`                                             |
//! | Physical meaning| Distance to nearest plate edge           | Distance to nearest oceanic-or-excluded crust                          |
//!
//! Forcing both into a single parameterised utility would either shift
//! cratonic's Step 9 numerical baseline (changing connectivity or
//! source predicate) or pollute the API with parameters that only one
//! caller would ever use. Split kept; documented here so the next
//! reviewer doesn't repeat the unification attempt.
//!
//! ## Algorithm
//!
//! 1. **Seed**: every cell whose 4-NESW periodic neighbour belongs to
//!    a different plate is a source at distance `0.0`. The propagated
//!    "across-boundary plate id" is the first such neighbour found in
//!    East-West-South-North order — deterministic.
//! 2. **Propagate**: Chebyshev 8-neighbour BFS on the periodic torus.
//!    Each hop adds `+1.0` to the distance. The target plate id
//!    propagates unchanged from the seed of each BFS chain.
//! 3. **Output**: a [`Field2D`] of distances and a `Vec<u16>` of
//!    target plate ids (row-major, `j * nx + i`).
//!
//! Cells that the BFS never reaches (the degenerate case where a
//! single plate spans the entire torus and there is no inter-plate
//! boundary) keep `distance = f64::INFINITY` and
//! `target_plate_id = u16::MAX`. Callers handle this naturally via
//! their `dist >= width` short-circuit.
//!
//! ## Why `target_plate_id` instead of a typed value
//!
//! The two existing callers propagate different value types
//! (`init::Uniform` propagates a per-plate `S̃` reference value;
//! `plate_kinematic::field::build` propagates a `(vx, vy)` pair).
//! Returning the across-boundary plate id keeps the utility free of
//! caller-specific value types — each caller indexes its own per-plate
//! lookup (`per_plate_value` for init, `velocities` for kinematic) on
//! `target_plate_id[idx]`. Bit-identical with the pre-refactor
//! per-cell propagation because the propagated plate id of a BFS chain
//! determines the propagated value uniquely (per-plate properties are
//! constant within a plate).
//!
//! ## Determinism
//!
//! Output is byte-identical for fixed inputs. BFS seeds in row-major
//! order; the FIFO queue and `nd < dist[nidx]` first-write-wins
//! tie-break make the propagation deterministic. Verified by
//! [`tests::deterministic_same_inputs`].

use std::collections::VecDeque;

use super::super::field::Field2D;
use super::PlateIdField;

/// Output of [`compute_dist_to_inter_plate_boundary`].
#[derive(Clone)]
pub struct InterPlateBoundaryDist {
    /// BFS distance from each cell to the nearest inter-plate
    /// boundary, in cell units. `0.0` at boundary cells, increasing
    /// by `1.0` per Chebyshev hop into plate interiors.
    /// `f64::INFINITY` for cells the BFS never reaches.
    pub distance: Field2D,

    /// Per-cell plate id of the across-boundary neighbour propagated
    /// from the BFS seed, row-major (`j * nx + i`). For a boundary
    /// source cell, this is the plate id of its first NESW neighbour
    /// belonging to a different plate (E, W, S, N order). For an
    /// interior cell, it is inherited unchanged from the seed of its
    /// BFS chain. `u16::MAX` for cells the BFS never reached.
    ///
    /// Callers index this into a per-plate lookup table to recover
    /// the "across-boundary value" without the utility needing to
    /// know what kind of value the caller is propagating.
    pub target_plate_id: Vec<u16>,
}

/// Compute the BFS distance to the nearest inter-plate boundary on the
/// periodic torus, using Chebyshev 8-connectivity, plus the plate id
/// of the across-boundary neighbour for each cell.
///
/// See module docstring for the algorithm rationale and the comparison
/// with [`super::super::cratonic::factor`] (which uses a different
/// BFS and is not refactored).
pub fn compute_dist_to_inter_plate_boundary(
    nx: usize,
    ny: usize,
    plate_id: &PlateIdField,
) -> InterPlateBoundaryDist {
    let n = nx * ny;
    let mut dist_buf = vec![f64::INFINITY; n];
    let mut target = vec![u16::MAX; n];
    let mut queue: VecDeque<(usize, usize)> = VecDeque::new();

    // Seed: every cell whose 4-NESW periodic neighbour belongs to a
    // different plate is a source at distance 0. The propagated
    // across-boundary plate id is the first such neighbour in
    // East-West-South-North order (deterministic). The NESW probe is
    // the same for all callers — the BFS itself is 8-connected.
    for j in 0..ny {
        for i in 0..nx {
            let id = plate_id.get(i, j);
            let ip = (i + 1) % nx;
            let im = (i + nx - 1) % nx;
            let jp = (j + 1) % ny;
            let jm = (j + ny - 1) % ny;
            let mut nb_pid: Option<u16> = None;
            for &(ni, nj) in &[(ip, j), (im, j), (i, jp), (i, jm)] {
                let other = plate_id.get(ni, nj);
                if other != id {
                    nb_pid = Some(other);
                    break;
                }
            }
            if let Some(pid) = nb_pid {
                let idx = j * nx + i;
                dist_buf[idx] = 0.0;
                target[idx] = pid;
                queue.push_back((i, j));
            }
        }
    }

    // Propagate: Chebyshev 8-neighbour BFS on the periodic torus.
    // Each hop adds 1.0; target plate id propagates unchanged from
    // the seed of each BFS chain (first-write-wins tie-break).
    while let Some((i, j)) = queue.pop_front() {
        let idx = j * nx + i;
        let d = dist_buf[idx];
        let tpid = target[idx];
        for dj in [-1_i32, 0, 1] {
            for di in [-1_i32, 0, 1] {
                if di == 0 && dj == 0 {
                    continue;
                }
                let ni = ((i as i32 + di).rem_euclid(nx as i32)) as usize;
                let nj = ((j as i32 + dj).rem_euclid(ny as i32)) as usize;
                let nidx = nj * nx + ni;
                let nd = d + 1.0;
                if nd < dist_buf[nidx] {
                    dist_buf[nidx] = nd;
                    target[nidx] = tpid;
                    queue.push_back((ni, nj));
                }
            }
        }
    }

    let mut distance = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            distance.set(i, j, dist_buf[j * nx + i]);
        }
    }

    InterPlateBoundaryDist { distance, target_plate_id: target }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

    fn make_plates(
        nx: usize,
        ny: usize,
        seed: u64,
        num_plates: usize,
    ) -> crate::tectonics_v2::voronoi::VoronoiPlates {
        let cfg = VoronoiConfig { num_plates, continental_ratio: 0.3 };
        generate_voronoi(nx, ny, &cfg, seed)
    }

    /// Boundary-source cells (4-NESW neighbour on a different plate)
    /// have distance exactly 0; their target_plate_id matches the
    /// plate id of their first NESW different-plate neighbour.
    #[test]
    fn boundary_cells_have_zero_distance() {
        let nx = 32;
        let ny = 32;
        let plates = make_plates(nx, ny, 42, 6);
        let bfs = compute_dist_to_inter_plate_boundary(nx, ny, &plates.plate_id);

        for j in 0..ny {
            for i in 0..nx {
                let id = plates.plate_id.get(i, j);
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                let mut expected_target: Option<u16> = None;
                for &(ni, nj) in &[(ip, j), (im, j), (i, jp), (i, jm)] {
                    let other = plates.plate_id.get(ni, nj);
                    if other != id {
                        expected_target = Some(other);
                        break;
                    }
                }
                let idx = j * nx + i;
                if let Some(tgt) = expected_target {
                    assert_eq!(
                        bfs.distance.get(i, j),
                        0.0,
                        "boundary cell ({},{}) expected dist 0, got {}",
                        i,
                        j,
                        bfs.distance.get(i, j)
                    );
                    assert_eq!(
                        bfs.target_plate_id[idx], tgt,
                        "boundary cell ({},{}) expected target plate {}, got {}",
                        i, j, tgt, bfs.target_plate_id[idx]
                    );
                }
            }
        }
    }

    /// Single-plate-on-torus degenerate case: no inter-plate boundary
    /// exists, BFS never seeds, all cells keep the sentinel
    /// (`INFINITY`, `u16::MAX`).
    #[test]
    fn single_plate_degenerate_case() {
        let nx = 16;
        let ny = 16;
        let plates = make_plates(nx, ny, 3, 1);
        let bfs = compute_dist_to_inter_plate_boundary(nx, ny, &plates.plate_id);
        for j in 0..ny {
            for i in 0..nx {
                assert_eq!(bfs.distance.get(i, j), f64::INFINITY);
                assert_eq!(bfs.target_plate_id[j * nx + i], u16::MAX);
            }
        }
    }

    /// Determinism: same `plate_id` field → byte-identical output.
    #[test]
    fn deterministic_same_inputs() {
        let nx = 32;
        let ny = 32;
        let plates_a = make_plates(nx, ny, 42, 5);
        let plates_b = make_plates(nx, ny, 42, 5);
        let a = compute_dist_to_inter_plate_boundary(nx, ny, &plates_a.plate_id);
        let b = compute_dist_to_inter_plate_boundary(nx, ny, &plates_b.plate_id);
        assert_eq!(a.distance.data(), b.distance.data());
        assert_eq!(a.target_plate_id, b.target_plate_id);
    }

    /// Distances grow by at most 1 per Chebyshev step — a basic BFS
    /// invariant. Equivalently, `|d[i, j] - d[ni, nj]| <= 1` for
    /// every 8-neighbour pair.
    #[test]
    fn distances_grow_at_most_one_per_step() {
        let nx = 32;
        let ny = 32;
        let plates = make_plates(nx, ny, 42, 4);
        let bfs = compute_dist_to_inter_plate_boundary(nx, ny, &plates.plate_id);
        for j in 0..ny {
            for i in 0..nx {
                let d = bfs.distance.get(i, j);
                if !d.is_finite() {
                    continue;
                }
                for dj in [-1_i32, 0, 1] {
                    for di in [-1_i32, 0, 1] {
                        if di == 0 && dj == 0 {
                            continue;
                        }
                        let ni = ((i as i32 + di).rem_euclid(nx as i32)) as usize;
                        let nj = ((j as i32 + dj).rem_euclid(ny as i32)) as usize;
                        let nd = bfs.distance.get(ni, nj);
                        if !nd.is_finite() {
                            continue;
                        }
                        assert!(
                            (d - nd).abs() <= 1.0 + 1e-12,
                            "Chebyshev neighbours ({},{}) and ({},{}) differ by {} > 1",
                            i,
                            j,
                            ni,
                            nj,
                            (d - nd).abs()
                        );
                    }
                }
            }
        }
    }
}
