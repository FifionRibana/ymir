//! Distance-to-convergent-boundary field for the C1 orogenic
//! closure (Issue #123 Stage 3).
//!
//! ## What this module does
//!
//! For each cell on the grid, computes the **octile distance**
//! (8-neighbour Dijkstra with cardinal cost 1 and diagonal cost
//! √2) to the nearest cell classified as
//! [`BoundaryType::Convergent`] by
//! [`super::boundary_classification::classify_boundaries`]. Cells
//! beyond `max_distance` from any convergent cell keep the
//! sentinel value `max_distance`.
//!
//! The result is consumed by the Davis-Suppe source term in
//! Stage 4: `h_critical(d) = h_max · (1 − exp(−d / L_taper))` and
//! `∂S̃/∂t ∝ exp(−d / L_decay)`.
//!
//! ## Distance metric
//!
//! Octile distance approximates Euclidean distance with a max
//! relative error of `(√2 − 1) ≈ 41 %` on pure-diagonal paths
//! and `~ 8 %` on the worst-case "knight's-move" path. For the
//! C1 closure profile that decays exponentially with distance,
//! this approximation is more than sufficient. True Euclidean
//! distance via Eikonal-equation marching is the upgrade path
//! if Phase 4 visual review surfaces staircase artefacts.
//!
//! ## Periodic boundaries
//!
//! Wraps through [`PeriodicIndex`] (reuses v2 infrastructure
//! verbatim per §4.8 design doc). A convergent boundary at the
//! east edge of the grid pulls cells on the west edge into its
//! influence zone, as the underlying physics requires.
//!
//! ## Complexity
//!
//! Dijkstra with a binary heap is `O(N log N)` where `N = nx · ny`.
//! At 64² ≈ 4096 cells the per-step cost is sub-millisecond, well
//! below the Phase 1.1 250 ms / 300-step budget. The
//! `max_distance` cap also prunes the search frontier — once a
//! cell exits the heap with `d ≥ max_distance` we stop relaxing
//! from it.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

use super::boundary_classification::{BoundaryInfo, BoundaryType};

/// √2 in `f64` (avoid the `f64::sqrt` cost in the inner loop).
const DIAG: f64 = std::f64::consts::SQRT_2;

/// Local newtype to give `f64` a total order for the priority
/// queue. Safe because all distances in this module are finite
/// non-negative — we never push NaN.
#[derive(Copy, Clone, PartialEq, PartialOrd)]
struct OrdF64(f64);

impl Eq for OrdF64 {}

impl Ord for OrdF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        // Distances are always finite non-negative; `partial_cmp`
        // never returns `None` here.
        self.0.partial_cmp(&other.0).expect("OrdF64 never holds NaN in distance_field")
    }
}

/// Compute the per-cell octile distance to the nearest
/// `Convergent` cell in `boundary`. Cells beyond `max_distance`
/// keep the sentinel value `max_distance`.
pub fn distance_to_convergent_boundary(boundary: &BoundaryInfo, max_distance: f64) -> Field2D {
    let nx = boundary.boundary_type.nx();
    let ny = boundary.boundary_type.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut dist = Field2D::filled(nx, ny, max_distance);

    // Seed the BFS with every Convergent cell at distance 0. Push
    // them all into the priority queue up front; the relaxation
    // loop below handles the rest.
    let mut heap: BinaryHeap<std::cmp::Reverse<(OrdF64, usize, usize)>> = BinaryHeap::new();
    for j in 0..ny {
        for i in 0..nx {
            if matches!(boundary.boundary_type.get(i, j), BoundaryType::Convergent) {
                dist.set(i, j, 0.0);
                heap.push(std::cmp::Reverse((OrdF64(0.0), i, j)));
            }
        }
    }

    // 8-neighbour offsets + edge costs.
    let neighbours: [(i32, i32, f64); 8] = [
        (1, 0, 1.0),
        (-1, 0, 1.0),
        (0, 1, 1.0),
        (0, -1, 1.0),
        (1, 1, DIAG),
        (1, -1, DIAG),
        (-1, 1, DIAG),
        (-1, -1, DIAG),
    ];

    while let Some(std::cmp::Reverse((OrdF64(d), i, j))) = heap.pop() {
        // Stale entry — a shorter path to (i, j) was found
        // after this entry was queued.
        if d > dist.get(i, j) {
            continue;
        }
        // Cells at or beyond the cap don't extend the frontier.
        if d >= max_distance {
            continue;
        }
        for &(di, dj, cost) in neighbours.iter() {
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
            let new_d = d + cost;
            if new_d < dist.get(ni, nj) && new_d < max_distance {
                dist.set(ni, nj, new_d);
                heap.push(std::cmp::Reverse((OrdF64(new_d), ni, nj)));
            }
        }
    }

    dist
}

/// Diagnostic helper: count cells whose distance is strictly less
/// than `threshold`. Useful for the Stage 3 cycle-0 verbose
/// output and for future Stage 4 source-term coverage stats.
pub fn count_cells_within(distance: &Field2D, threshold: f64) -> usize {
    distance.data().iter().filter(|&&d| d < threshold).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_c1::boundary_classification::{BoundaryTypeField, BoundaryType};
    use crate::tectonics_c1::state::BoolField;

    /// Build a [`BoundaryInfo`] where each cell takes its boundary
    /// type from a closure. Upper-plate mask is set false
    /// everywhere — distance-field is mask-agnostic.
    fn build_boundary<F: Fn(usize, usize) -> BoundaryType>(
        nx: usize,
        ny: usize,
        f: F,
    ) -> BoundaryInfo {
        let mut bt = BoundaryTypeField::filled(nx, ny, BoundaryType::Internal);
        for j in 0..ny {
            for i in 0..nx {
                bt.set(i, j, f(i, j));
            }
        }
        BoundaryInfo {
            boundary_type: bt,
            upper_plate_mask: BoolField::filled(nx, ny, false),
        }
    }

    #[test]
    fn distance_zero_at_convergent_cells() {
        // Single isolated Convergent cell at (4, 4) on an 8x8
        // grid. Distance to itself is 0.
        let boundary =
            build_boundary(8, 8, |i, j| if (i, j) == (4, 4) { BoundaryType::Convergent } else { BoundaryType::Internal });
        let dist = distance_to_convergent_boundary(&boundary, 10.0);
        assert_eq!(dist.get(4, 4), 0.0);
    }

    #[test]
    fn distance_one_at_immediate_neighbours_cardinal_and_sqrt2_diagonal() {
        // Isolated Convergent at the centre of a 16x16 grid (large
        // enough that periodic wrap doesn't compete for short
        // distances). Cardinal neighbours land at d=1; diagonal at
        // d=√2.
        let center = 8;
        let boundary = build_boundary(16, 16, |i, j| {
            if (i, j) == (center, center) {
                BoundaryType::Convergent
            } else {
                BoundaryType::Internal
            }
        });
        let dist = distance_to_convergent_boundary(&boundary, 10.0);
        // Cardinal neighbours
        assert_eq!(dist.get(center + 1, center), 1.0);
        assert_eq!(dist.get(center - 1, center), 1.0);
        assert_eq!(dist.get(center, center + 1), 1.0);
        assert_eq!(dist.get(center, center - 1), 1.0);
        // Diagonal neighbours
        assert!((dist.get(center + 1, center + 1) - DIAG).abs() < 1e-12);
        assert!((dist.get(center - 1, center + 1) - DIAG).abs() < 1e-12);
        assert!((dist.get(center + 1, center - 1) - DIAG).abs() < 1e-12);
        assert!((dist.get(center - 1, center - 1) - DIAG).abs() < 1e-12);
    }

    #[test]
    fn distance_periodic_wraparound() {
        // Convergent cell at (0, 4) on a 16x8 grid. The west
        // neighbour wraps around to (15, 4); the diagonal NW
        // wraps to (15, 5). Both should be 1.0 and √2 cells
        // away respectively.
        let boundary = build_boundary(16, 8, |i, j| {
            if (i, j) == (0, 4) { BoundaryType::Convergent } else { BoundaryType::Internal }
        });
        let dist = distance_to_convergent_boundary(&boundary, 10.0);
        assert_eq!(dist.get(15, 4), 1.0, "west wrap should be 1.0");
        assert!((dist.get(15, 5) - DIAG).abs() < 1e-12, "NW wrap should be √2");
        assert!((dist.get(15, 3) - DIAG).abs() < 1e-12, "SW wrap should be √2");
    }

    #[test]
    fn distance_capped_at_max() {
        // Single Convergent cell at the centre of a 32x32 grid
        // with max_distance = 3.0. Cells at distance > 3 must
        // retain the sentinel value 3.0.
        let nx = 32;
        let ny = 32;
        let center_i = 16;
        let center_j = 16;
        let boundary = build_boundary(nx, ny, |i, j| {
            if (i, j) == (center_i, center_j) {
                BoundaryType::Convergent
            } else {
                BoundaryType::Internal
            }
        });
        let max_d = 3.0;
        let dist = distance_to_convergent_boundary(&boundary, max_d);

        // Corner cell (0, 0) — chessboard distance to (16, 16)
        // is 16. Way past the cap → sentinel.
        assert_eq!(dist.get(0, 0), max_d);
        // Cell at (center + 4, center) — cardinal distance 4 > 3,
        // sentinel.
        assert_eq!(dist.get(center_i + 4, center_j), max_d);
        // Cell at (center + 2, center + 2) — diagonal distance
        // 2·√2 ≈ 2.83 < 3, should be reached.
        assert!(
            dist.get(center_i + 2, center_j + 2) < max_d,
            "cell at 2·√2 distance should land below cap; got {}",
            dist.get(center_i + 2, center_j + 2)
        );
    }

    /// Cycle-0 stats for the actual Phase 1.1 init state — pre-
    /// view of the orogenic-source-term coverage that Stage 4
    /// will exercise.
    #[test]
    fn phase_1_1_init_state_cycle_0_distance_stats() {
        use crate::tectonics_c1::boundary_classification::classify_boundaries;
        use crate::tectonics_c1::init::init_c1_state_phase_1_1;
        use crate::tectonics_c1::kinematics::PlateKinematics;

        let state = init_c1_state_phase_1_1(64, 42);
        let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let boundary = classify_boundaries(&state.plate_id, &kinematics);

        // Defaults aligned with Stage 4 plan: L_decay = 6.0,
        // max_distance = 5 · L_decay = 30.
        let l_decay = 6.0_f64;
        let max_d = 5.0 * l_decay;
        let dist = distance_to_convergent_boundary(&boundary, max_d);

        let data = dist.data();
        let total = data.len() as f64;
        let count_zero = data.iter().filter(|&&d| d == 0.0).count();
        let count_within_l_decay = count_cells_within(&dist, l_decay);
        let count_within_max = count_cells_within(&dist, max_d);

        let mut min = f64::INFINITY;
        let mut max = 0.0_f64;
        let mut sum = 0.0_f64;
        for &d in data {
            if d < min {
                min = d;
            }
            if d > max {
                max = d;
            }
            sum += d;
        }
        let mean = sum / total;

        eprintln!("Phase 1.1 cycle-0 distance-to-Convergent stats (64² seed 42, L_decay = {l_decay:.1}, max = {max_d:.1}):");
        eprintln!("  min                = {min:.4}");
        eprintln!("  mean               = {mean:.4}");
        eprintln!("  max                = {max:.4}");
        eprintln!(
            "  d = 0   (seeds)    = {} cells   ({:.1} %)  — Convergent count",
            count_zero,
            100.0 * count_zero as f64 / total
        );
        eprintln!(
            "  d < L_decay        = {} cells   ({:.1} %)  — primary influence zone",
            count_within_l_decay,
            100.0 * count_within_l_decay as f64 / total
        );
        eprintln!(
            "  d < max_distance   = {} cells   ({:.1} %)  — any source contribution",
            count_within_max,
            100.0 * count_within_max as f64 / total
        );

        // Stage 2 reported 289 Convergent cells; verify our seed
        // count matches.
        assert_eq!(count_zero, 289, "Convergent seed count should match Stage 2 stats");
        // Internal sanity: the primary influence zone (d <
        // L_decay) is bigger than the seed count.
        assert!(count_within_l_decay > count_zero);
        // The max distance after relaxation should not exceed the
        // configured cap.
        assert!(max <= max_d);
    }
}
