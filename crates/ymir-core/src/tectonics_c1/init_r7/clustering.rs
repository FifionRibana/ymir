//! R7 continental clustering — BFS cluster-based plate-type
//! assignment. Phase 2 Track B sub-component 2 (Issue #131).
//!
//! See [`super`] for module-level rationale.
//!
//! ## Design
//!
//! [`assign_continental_clusters`] OVERRIDES the per-plate type
//! produced by `tectonics_v2::voronoi::generate_voronoi` (whose
//! independent Bernoulli draws scatter continental plates
//! uniformly across the torus). The BFS-expansion approach
//! produces **geographically clustered** continental plates per
//! §6.2 of the design doc, satisfying the §2.4 viewport-cadrable
//! requirement.
//!
//! Default `continental_fraction = 0.29` matches Earth's
//! continental-to-oceanic ratio and is **independent** of the
//! `VoronoiConfig::default().continental_ratio = 0.30` (which
//! remains unchanged so Phase 1.x init keeps its regression
//! baseline). Phase 2 Track B's BFS output supersedes the
//! Voronoï Bernoulli ratio when the Track B pipeline is
//! invoked.
//!
//! ## Periodic-boundary wrap risk
//!
//! C1's torus is fully periodic. The plate adjacency graph
//! derived by [`build_plate_adjacency`] uses `rem_euclid` wrapping
//! so plates on opposite sides of the domain that share a periodic
//! edge are correctly adjacent. **Consequence:** BFS from a single
//! seed could grow a continental cluster that wraps around the
//! torus and therefore does not render as a single connected blob
//! in a planar 512² viewport.
//!
//! Mitigation: Stage A acceptance test
//! `acceptance_track_b2_continent_cadrable` verifies wrap-free
//! behaviour on a real 512² init. Stage E2 unit tests assert the
//! "single connected component" property on the adjacency graph
//! (algorithm correctness) but cannot detect torus-wrap on
//! synthetic adjacency without a coordinate model. Two-layer
//! testing per W7 surface decision.

use std::collections::{HashSet, VecDeque};

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::tectonics_v2::boundaries::plate_type::PlateType;
use crate::tectonics_v2::voronoi::PlateIdField;

/// Continental clustering parameters.
///
/// Defaults target Earth-like continental fraction with a single
/// contiguous continental cluster (cadrable per §2.4 viewport
/// requirement). When `enabled = false`,
/// [`assign_continental_clusters`] is a no-op (W4 closure-
/// isolation discipline).
#[derive(Clone, Copy, Debug)]
pub struct ContinentalClusterParams {
    /// Master enable/disable.
    pub enabled: bool,
    /// Target continental fraction of `num_plates`. Default 0.29
    /// (Earth-like). The BFS expansion stops once the round-trip
    /// integer count `round(num_plates · continental_fraction)`
    /// is reached.
    pub continental_fraction: f64,
    /// Number of BFS seed plates. `1` produces a single
    /// contiguous cluster (cadrable). Larger values produce
    /// multiple disjoint continents.
    pub seed_cluster_count: usize,
    /// Deterministic seed for the RNG that picks BFS seed
    /// plate(s) from the adjacency graph.
    pub seed: u64,
}

impl Default for ContinentalClusterParams {
    fn default() -> Self {
        Self {
            enabled: true,
            continental_fraction: 0.29,
            seed_cluster_count: 1,
            seed: 0,
        }
    }
}

/// Build the per-plate adjacency graph from a `PlateIdField`.
///
/// Two plates are adjacent if any of their cells share a
/// 4-neighbour edge under periodic (`rem_euclid`) wraparound. The
/// output `Vec<Vec<u16>>` is indexed by plate id; each inner
/// `Vec<u16>` holds the sorted unique neighbour ids.
///
/// `Vec<Vec<u16>>` rather than `HashMap<u16, HashSet<u16>>`
/// because plate ids are dense `0..num_plates` (Voronoï convention)
/// — `Vec`-indexed access is faster and avoids hash overhead at
/// 8-plate typical scale.
///
/// O(`nx · ny`) time, single pass over cells.
pub fn build_plate_adjacency(
    plate_id_field: &PlateIdField,
    num_plates: usize,
) -> Vec<Vec<u16>> {
    let nx = plate_id_field.nx();
    let ny = plate_id_field.ny();
    let mut adjacency_sets: Vec<HashSet<u16>> =
        (0..num_plates).map(|_| HashSet::new()).collect();

    let offsets: [(i32, i32); 4] = [(1, 0), (-1, 0), (0, 1), (0, -1)];
    let nx_i = nx as i32;
    let ny_i = ny as i32;

    for j in 0..ny {
        for i in 0..nx {
            let pid_c = plate_id_field.get(i, j);
            for &(di, dj) in offsets.iter() {
                let ni = (i as i32 + di).rem_euclid(nx_i) as usize;
                let nj = (j as i32 + dj).rem_euclid(ny_i) as usize;
                let pid_n = plate_id_field.get(ni, nj);
                if pid_n != pid_c {
                    let bucket = pid_c as usize;
                    if bucket < num_plates {
                        adjacency_sets[bucket].insert(pid_n);
                    }
                }
            }
        }
    }

    adjacency_sets
        .into_iter()
        .map(|s| {
            let mut v: Vec<u16> = s.into_iter().collect();
            v.sort_unstable();
            v
        })
        .collect()
}

/// Assign plate types via BFS cluster-based growth.
///
/// Algorithm:
///
/// 1. Initialise all plates as `PlateType::Oceanic`.
/// 2. Sample `params.seed_cluster_count` BFS seed plates
///    deterministically from `params.seed` via `ChaCha8Rng`.
/// 3. BFS expand from the seeds simultaneously across the
///    adjacency graph. Mark each visited plate as
///    `PlateType::Continental`.
/// 4. Stop when the target count
///    `round(num_plates · params.continental_fraction)` is
///    reached, OR when BFS exhausts (disconnected adjacency
///    component — target unreachable from these seeds).
///
/// Inputs:
/// - `per_plate_type`: mutable slice indexed by plate id,
///   typically a copy of `VoronoiPlates::per_plate_type`
///   constructed by the Stage E4 dispatcher.
/// - `adjacency`: per-plate adjacency from
///   [`build_plate_adjacency`].
/// - `params`: clustering tunables.
///
/// `params.enabled = false` is a no-op (W4 closure-isolation).
pub fn assign_continental_clusters(
    per_plate_type: &mut [PlateType],
    adjacency: &[Vec<u16>],
    params: &ContinentalClusterParams,
) {
    if !params.enabled {
        return;
    }
    let num_plates = per_plate_type.len();
    if num_plates == 0 {
        return;
    }
    debug_assert_eq!(
        adjacency.len(),
        num_plates,
        "adjacency length ({}) must equal num_plates ({})",
        adjacency.len(),
        num_plates
    );

    // Target continental plate count. Clamped to `[1, num_plates]`
    // — a non-zero seed always plants at least one continental
    // plate; a fraction `>= 1.0` clamps to all plates continental.
    let target_continental =
        ((num_plates as f64 * params.continental_fraction).round() as usize)
            .clamp(1, num_plates);

    // Step 1 — initialise all plates as Oceanic.
    for t in per_plate_type.iter_mut() {
        *t = PlateType::Oceanic;
    }

    // Step 2 — pick BFS seed plates deterministically.
    let mut rng = ChaCha8Rng::seed_from_u64(params.seed);
    let seed_count = params.seed_cluster_count.clamp(1, num_plates);
    let mut available: Vec<u16> = (0..num_plates as u16).collect();
    let mut seeds: Vec<u16> = Vec::with_capacity(seed_count);
    for _ in 0..seed_count {
        if available.is_empty() {
            break;
        }
        let idx = (rng.random::<u64>() as usize) % available.len();
        seeds.push(available.swap_remove(idx));
    }

    // Step 3 — BFS from seeds, marking Continental as we go.
    let mut visited: Vec<bool> = vec![false; num_plates];
    let mut queue: VecDeque<u16> = VecDeque::new();
    let mut continental_count = 0_usize;
    for &s in seeds.iter() {
        let idx = s as usize;
        if !visited[idx] {
            visited[idx] = true;
            per_plate_type[idx] = PlateType::Continental;
            continental_count += 1;
            queue.push_back(s);
            if continental_count >= target_continental {
                break;
            }
        }
    }

    // Step 4 — BFS expansion.
    while continental_count < target_continental {
        let Some(p) = queue.pop_front() else {
            // Adjacency component exhausted. Target unreachable
            // from these seeds; documented behaviour — the
            // disconnected-Voronoï unit test exercises this.
            break;
        };
        for &n in adjacency[p as usize].iter() {
            let idx = n as usize;
            if !visited[idx] {
                visited[idx] = true;
                per_plate_type[idx] = PlateType::Continental;
                continental_count += 1;
                queue.push_back(n);
                if continental_count >= target_continental {
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fully-connected adjacency for `n` plates (complete graph).
    fn complete_adjacency(n: usize) -> Vec<Vec<u16>> {
        (0..n)
            .map(|i| (0..n).filter(|j| *j != i).map(|j| j as u16).collect())
            .collect()
    }

    /// Ring adjacency: 0-1-2-...-(n-1)-0.
    fn ring_adjacency(n: usize) -> Vec<Vec<u16>> {
        (0..n)
            .map(|i| {
                let prev = ((i + n - 1) % n) as u16;
                let next = ((i + 1) % n) as u16;
                vec![prev.min(next), prev.max(next)]
            })
            .collect()
    }

    /// Disconnected adjacency: 4 plates split into two components
    /// `{0, 1}` (connected) and `{2, 3}` (connected). No edge
    /// between components.
    fn disconnected_adjacency() -> Vec<Vec<u16>> {
        vec![vec![1], vec![0], vec![3], vec![2]]
    }

    /// Test 1 — deterministic given seed.
    #[test]
    fn clustering_deterministic_given_seed() {
        let adj = complete_adjacency(8);
        let mut a = vec![PlateType::Continental; 8]; // start non-default
        let mut b = vec![PlateType::Continental; 8];
        let params = ContinentalClusterParams::default();

        assign_continental_clusters(&mut a, &adj, &params);
        assign_continental_clusters(&mut b, &adj, &params);

        for i in 0..8 {
            assert_eq!(
                a[i], b[i],
                "same (adj, params) must produce bit-identical output at plate {i}"
            );
        }
    }

    /// Test 2 — continental fraction within target. Tolerance
    /// `max(5 %, 1 / num_plates)` accommodates discrete plate
    /// granularity (at `num_plates = 8`, one plate = 12.5 %).
    #[test]
    fn clustering_continental_fraction_within_target() {
        let num_plates = 8;
        let adj = complete_adjacency(num_plates);
        let mut types = vec![PlateType::Oceanic; num_plates];
        let params = ContinentalClusterParams {
            continental_fraction: 0.29,
            ..ContinentalClusterParams::default()
        };

        assign_continental_clusters(&mut types, &adj, &params);

        let continental_count = types
            .iter()
            .filter(|t| matches!(t, PlateType::Continental))
            .count();
        let actual_fraction = continental_count as f64 / num_plates as f64;
        let tolerance = (1.0 / num_plates as f64).max(0.05);
        let diff = (actual_fraction - params.continental_fraction).abs();
        eprintln!(
            "clustering fraction test: actual = {actual_fraction:.3} \
             ({continental_count} / {num_plates}), target = {:.3}, diff = {diff:.3}, tolerance = {tolerance:.3}",
            params.continental_fraction
        );
        assert!(
            diff <= tolerance,
            "continental fraction {actual_fraction:.3} too far from target {:.3} \
             (diff = {diff:.3} > tolerance {tolerance:.3})",
            params.continental_fraction
        );
    }

    /// Test 3 — single seed → single connected continental
    /// subgraph in the adjacency. Verified by BFS reachability
    /// check on the continental subgraph after assignment.
    #[test]
    fn clustering_single_contiguous_component() {
        let num_plates = 8;
        let adj = ring_adjacency(num_plates);
        let mut types = vec![PlateType::Oceanic; num_plates];
        let params = ContinentalClusterParams {
            continental_fraction: 0.5,
            seed_cluster_count: 1,
            seed: 42,
            ..ContinentalClusterParams::default()
        };

        assign_continental_clusters(&mut types, &adj, &params);

        let continental: Vec<u16> = types
            .iter()
            .enumerate()
            .filter_map(|(i, t)| matches!(t, PlateType::Continental).then(|| i as u16))
            .collect();
        assert!(
            !continental.is_empty(),
            "single-seed assignment must produce at least one continental plate"
        );

        // BFS reachability within the continental subgraph from
        // continental[0]: should reach every other continental plate.
        let mut reached: HashSet<u16> = HashSet::new();
        let mut queue: VecDeque<u16> = VecDeque::new();
        reached.insert(continental[0]);
        queue.push_back(continental[0]);
        while let Some(p) = queue.pop_front() {
            for &n in adj[p as usize].iter() {
                if continental.contains(&n) && !reached.contains(&n) {
                    reached.insert(n);
                    queue.push_back(n);
                }
            }
        }
        assert_eq!(
            reached.len(),
            continental.len(),
            "single-seed continental subgraph must be connected; \
             reached {} of {} continental plates",
            reached.len(),
            continental.len()
        );
    }

    /// Test 4 — `enabled = false` no-op.
    #[test]
    fn clustering_disabled_no_op() {
        let adj = complete_adjacency(8);
        // Pre-populate with a non-uniform pattern so a forgotten
        // branch can't pass on a default-Oceanic input.
        let initial = vec![
            PlateType::Continental,
            PlateType::Oceanic,
            PlateType::Continental,
            PlateType::Oceanic,
            PlateType::Oceanic,
            PlateType::Continental,
            PlateType::Oceanic,
            PlateType::Oceanic,
        ];
        let mut types = initial.clone();
        let params = ContinentalClusterParams {
            enabled: false,
            continental_fraction: 0.99, // pathological if it ran
            ..ContinentalClusterParams::default()
        };

        assign_continental_clusters(&mut types, &adj, &params);

        for i in 0..8 {
            assert_eq!(
                types[i], initial[i],
                "`enabled = false` must not touch plate {i}"
            );
        }
    }

    /// Test 5 — disconnected adjacency robustness. With a seed in
    /// component `{0, 1}` and `continental_fraction = 0.75`
    /// (target 3 of 4 plates), BFS can only reach `{0, 1}` (2
    /// plates). The function must complete without panic and
    /// produce a valid (if undersized) continental cluster.
    #[test]
    fn clustering_handles_disconnected_voronoi() {
        let adj = disconnected_adjacency();
        let mut types = vec![PlateType::Oceanic; 4];
        // Seed pick is deterministic via ChaCha8Rng(0). Component
        // {0, 1} or {2, 3} — either is fine for this test.
        let params = ContinentalClusterParams {
            continental_fraction: 0.75,
            seed_cluster_count: 1,
            seed: 0,
            ..ContinentalClusterParams::default()
        };

        assign_continental_clusters(&mut types, &adj, &params);

        let continental_count = types
            .iter()
            .filter(|t| matches!(t, PlateType::Continental))
            .count();
        // Disconnected component has 2 plates; target was 3.
        // Achievable count is 2.
        assert_eq!(
            continental_count, 2,
            "disconnected adjacency: BFS should reach exactly the seed component (2 plates); got {continental_count}"
        );
        // No panic, no infinite loop — the absence of those is
        // the principal assertion.
    }

    /// Smoke test — `build_plate_adjacency` on a real 2-plate
    /// `PlateIdField` produces a non-empty adjacency.
    #[test]
    fn build_adjacency_on_two_plate_field() {
        let nx = 16;
        let ny = 16;
        let mut p = PlateIdField::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                p.set(i, j, if i < nx / 2 { 0 } else { 1 });
            }
        }
        let adj = build_plate_adjacency(&p, 2);
        assert_eq!(adj.len(), 2);
        assert_eq!(adj[0], vec![1], "plate 0 should be adjacent to plate 1");
        assert_eq!(adj[1], vec![0], "plate 1 should be adjacent to plate 0");
    }
}
