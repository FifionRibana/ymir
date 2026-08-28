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

/// **Intra-plate** Dijkstra constrained to stay on the seed's
/// plate during expansion. Returns per-cell octile distance to the
/// nearest [`super::state::BoolField`]-flagged cell **of the same
/// plate**, capped at `max_distance`.
///
/// ## Why
///
/// Stage 4 of Issue #123 surfaced that the Davis-Suppe wedge body
/// lives in the **interior of the upper plate**, not on the
/// convergent boundary itself. Seeding a generic distance-to-
/// boundary BFS from [`BoundaryType::Convergent`] cells gives `d=0`
/// for boundary cells, but those are the wrong cells to apply
/// `h_critical(d) = h_max · (1 − exp(−d/L_taper))` to (the formula
/// returns 0 at `d=0`, so the source term would *thin* the
/// boundary instead of thickening the upper-plate interior — anti-
/// geological).
///
/// This function builds a wedge-body-aware distance: seeds are
/// [`super::boundary_classification::BoundaryInfo::upper_plate_mask`]
/// cells (Convergent + faster-plate side), and expansion is
/// blocked at plate boundaries — so each plate's distance field
/// is independent. Stage 4 then applies the orogenic source term
/// where `d > 0 && d < max_distance` (interior cells reachable
/// from a same-plate upper-plate seed).
///
/// ## Plates without an upper-plate seed
///
/// The Phase 1.1 hand-tuned 8-plate preset has 4 cardinal-vs-
/// cardinal symmetric pairs (`|v|=0.01` each) where the faster-
/// plate heuristic returns no upper-plate flag. Those plates have
/// **no seeds in this function** and stay at `max_distance`
/// throughout — correct behaviour, and the resulting wedge body
/// is intentionally silent on those boundaries until the preset
/// or the kinematics heuristic is revised.
pub fn wedge_distance_intra_plate(
    plate_id: &crate::tectonics_v2::voronoi::PlateIdField,
    upper_plate_mask: &super::state::BoolField,
    max_distance: f64,
) -> Field2D {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut dist = Field2D::filled(nx, ny, max_distance);

    let mut heap: BinaryHeap<std::cmp::Reverse<(OrdF64, usize, usize)>> = BinaryHeap::new();
    for j in 0..ny {
        for i in 0..nx {
            if upper_plate_mask.get(i, j) {
                dist.set(i, j, 0.0);
                heap.push(std::cmp::Reverse((OrdF64(0.0), i, j)));
            }
        }
    }

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
        if d > dist.get(i, j) {
            continue;
        }
        if d >= max_distance {
            continue;
        }
        // Intra-plate constraint: only step to neighbours sharing
        // this cell's plate id. By induction, this restricts each
        // connected reachable region to a single plate (the seed's
        // plate).
        let plate_c = plate_id.get(i, j);
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
            if plate_id.get(ni, nj) != plate_c {
                continue;
            }
            let new_d = d + cost;
            if new_d < dist.get(ni, nj) && new_d < max_distance {
                dist.set(ni, nj, new_d);
                heap.push(std::cmp::Reverse((OrdF64(new_d), ni, nj)));
            }
        }
    }

    dist
}

/// #155 maillon 1b-i — typed variant of [`wedge_distance_intra_plate`]
/// that ALSO propagates each seed's convergence type to its reachable
/// wedge cells. Identical Dijkstra (same distances), plus a companion
/// `is_oc` mask: each reached cell inherits the O-C flag of the seed on
/// its shortest path (the nearest upper-plate seed). The Davis-Suppe
/// step then routes geometry by type — O-C cells get the margin-peaked
/// ridge profile (Andes), C-C / velocity-fallback cells keep the
/// rising-to-plateau dome (Tibet).
///
/// `oc_seed_mask` (from
/// [`super::boundary_classification::oc_override_seed_mask`]) must align
/// with the `true` cells of `upper_plate_mask`: it flags which of those
/// seeds are O-C continental overrides.
pub fn wedge_distance_intra_plate_typed(
    plate_id: &crate::tectonics_v2::voronoi::PlateIdField,
    upper_plate_mask: &super::state::BoolField,
    oc_seed_mask: &super::state::BoolField,
    max_distance: f64,
) -> (Field2D, super::state::BoolField) {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut dist = Field2D::filled(nx, ny, max_distance);
    let mut is_oc = super::state::BoolField::filled(nx, ny, false);

    let mut heap: BinaryHeap<std::cmp::Reverse<(OrdF64, usize, usize)>> = BinaryHeap::new();
    for j in 0..ny {
        for i in 0..nx {
            if upper_plate_mask.get(i, j) {
                dist.set(i, j, 0.0);
                is_oc.set(i, j, oc_seed_mask.get(i, j));
                heap.push(std::cmp::Reverse((OrdF64(0.0), i, j)));
            }
        }
    }

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
        if d > dist.get(i, j) {
            continue;
        }
        if d >= max_distance {
            continue;
        }
        let plate_c = plate_id.get(i, j);
        let oc_c = is_oc.get(i, j);
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
            if plate_id.get(ni, nj) != plate_c {
                continue;
            }
            let new_d = d + cost;
            if new_d < dist.get(ni, nj) && new_d < max_distance {
                dist.set(ni, nj, new_d);
                is_oc.set(ni, nj, oc_c);
                heap.push(std::cmp::Reverse((OrdF64(new_d), ni, nj)));
            }
        }
    }

    (dist, is_oc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_c1::boundary_classification::{BoundaryType, BoundaryTypeField};
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
        BoundaryInfo { boundary_type: bt, upper_plate_mask: BoolField::filled(nx, ny, false) }
    }

    #[test]
    fn distance_zero_at_convergent_cells() {
        // Single isolated Convergent cell at (4, 4) on an 8x8
        // grid. Distance to itself is 0.
        let boundary = build_boundary(8, 8, |i, j| {
            if (i, j) == (4, 4) { BoundaryType::Convergent } else { BoundaryType::Internal }
        });
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

    // ─────────────── wedge_distance_intra_plate tests ──────────

    /// Build a [`PlateIdField`] from a closure (test helper local
    /// to this `mod tests` block).
    fn build_plate_id_inline<F: Fn(usize, usize) -> u16>(
        nx: usize,
        ny: usize,
        f: F,
    ) -> crate::tectonics_v2::voronoi::PlateIdField {
        let mut p = crate::tectonics_v2::voronoi::PlateIdField::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                p.set(i, j, f(i, j));
            }
        }
        p
    }

    #[test]
    fn single_upper_plate_cell_expands_within_plate() {
        // 16x16 grid, entirely plate 0. Single upper_plate_mask
        // cell at the centre. The whole plate is reachable, so
        // distances grow as octile distance from the centre.
        let nx = 16;
        let ny = 16;
        let center = 8;
        let plate_id = build_plate_id_inline(nx, ny, |_, _| 0);
        let mut mask = BoolField::filled(nx, ny, false);
        mask.set(center, center, true);

        // max_distance set well above the worst-case octile
        // periodic distance to the centre (8·√2 ≈ 11.31) so the
        // single-plate / no-boundary case really does reach the
        // full grid.
        let max_d = 20.0;
        let dist = wedge_distance_intra_plate(&plate_id, &mask, max_d);

        // Seed itself.
        assert_eq!(dist.get(center, center), 0.0);
        // Cardinal neighbour
        assert_eq!(dist.get(center + 1, center), 1.0);
        // Diagonal neighbour
        assert!((dist.get(center + 1, center + 1) - DIAG).abs() < 1e-12);
        // No plate boundary anywhere → expansion reaches the
        // entire grid within cap.
        for j in 0..ny {
            for i in 0..nx {
                assert!(
                    dist.get(i, j) < max_d,
                    "single-plate grid: cell ({i},{j}) should be reached, got {}",
                    dist.get(i, j),
                );
            }
        }
    }

    #[test]
    fn plate_boundary_blocks_expansion() {
        // 16x16 split: i < 8 → plate 0, i >= 8 → plate 1. Single
        // upper_plate_mask cell at (4, 8) (plate 0 side, well
        // away from the boundary). Plate 1 cells must stay at
        // max_distance because the intra-plate constraint blocks
        // expansion across i=8.
        let nx = 16;
        let ny = 16;
        let plate_id = build_plate_id_inline(nx, ny, |i, _| if i < 8 { 0 } else { 1 });
        let mut mask = BoolField::filled(nx, ny, false);
        mask.set(4, 8, true);

        let max_d = 20.0;
        let dist = wedge_distance_intra_plate(&plate_id, &mask, max_d);

        // Plate 0 side: the seed and nearby cells are reached.
        assert_eq!(dist.get(4, 8), 0.0);
        assert!(dist.get(0, 8) < max_d, "plate 0 cell should be reachable");
        // Plate 1 side: every cell stays at sentinel (no seed
        // belongs to plate 1, and intra-plate constraint stops
        // any seed of plate 0 from crossing). Note: periodic wrap
        // around the East edge could in principle reach plate 1
        // through plate 0's east extent, but the wrap goes 0→15
        // staying on plate 1 → no plate-0 seed reaches plate 1.
        for j in 0..ny {
            for i in 8..nx {
                assert_eq!(
                    dist.get(i, j),
                    max_d,
                    "plate 1 cell ({i},{j}) must stay at sentinel; got {}",
                    dist.get(i, j),
                );
            }
        }
    }

    #[test]
    fn multi_plate_independent_expansion() {
        // 16x16 split into 2 plates. Both plates have one
        // upper_plate_mask cell each. Each seed expands within
        // its own plate independently; the two distance fields
        // are isolated by the intra-plate constraint.
        let nx = 16;
        let ny = 16;
        let plate_id = build_plate_id_inline(nx, ny, |i, _| if i < 8 { 0 } else { 1 });
        let mut mask = BoolField::filled(nx, ny, false);
        mask.set(2, 4, true); // plate 0 seed
        mask.set(12, 12, true); // plate 1 seed

        let max_d = 20.0;
        let dist = wedge_distance_intra_plate(&plate_id, &mask, max_d);

        assert_eq!(dist.get(2, 4), 0.0);
        assert_eq!(dist.get(12, 12), 0.0);

        // Cell near plate 0 seed (on plate 0): reached.
        assert!(dist.get(0, 0) < max_d, "plate 0 cell (0,0) should be reachable from plate 0 seed");
        // Cell near plate 1 seed (on plate 1): reached.
        assert!(
            dist.get(15, 15) < max_d,
            "plate 1 cell (15,15) should be reachable from plate 1 seed"
        );

        // Cross-plate sanity: a cell on plate 0 with the same
        // (i, j) as the plate 1 seed would be (12, 12) but
        // plate_id(12, 12) == 1, so that exact cell is plate 1.
        // Choose a different plate 0 cell visually close to the
        // plate 1 seed — e.g., (7, 12). Its distance to plate 0
        // seed (2, 4) is octile distance (5, 8) ≈ 5 + 3·√2 ≈ 9.24.
        // Its distance to plate 1 seed would have been small (~5
        // cardinal) but the intra-plate constraint blocks it.
        let d_07_12 = dist.get(7, 12);
        assert!(
            d_07_12 > 5.0,
            "plate 0 cell (7,12) should be reached via plate 0 seed only; got {}",
            d_07_12
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

        eprintln!(
            "Phase 1.1 cycle-0 distance-to-Convergent stats (64² seed 42, L_decay = {l_decay:.1}, max = {max_d:.1}):"
        );
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

    /// Cycle-0 stats for the **wedge** distance (Stage 3.1). The
    /// per-plate breakdown is the key remontée: it confirms that
    /// the 4 cardinal-vs-cardinal symmetric pairs are silent (no
    /// upper-plate seed → no expansion) and the 4 asymmetric
    /// pairs do carry seeds + interior expansion.
    #[test]
    fn phase_1_1_init_state_cycle_0_wedge_distance_stats() {
        use crate::tectonics_c1::boundary_classification::classify_boundaries;
        use crate::tectonics_c1::init::init_c1_state_phase_1_1;
        use crate::tectonics_c1::kinematics::PlateKinematics;

        let state = init_c1_state_phase_1_1(64, 42);
        let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let boundary = classify_boundaries(&state.plate_id, &kinematics);

        let l_decay = 6.0_f64;
        let max_d = 5.0 * l_decay;
        let dist = wedge_distance_intra_plate(&state.plate_id, &boundary.upper_plate_mask, max_d);
        let data = dist.data();
        let total = data.len() as f64;

        let seeds_total = boundary.upper_plate_count();
        let reached_total = count_cells_within(&dist, max_d);
        let primary_total = count_cells_within(&dist, l_decay);

        let mut min = f64::INFINITY;
        let mut max_reached = 0.0_f64;
        let mut sum = 0.0_f64;
        let mut reached_sample = 0_usize;
        for &d in data {
            if d < max_d {
                if d < min {
                    min = d;
                }
                if d > max_reached {
                    max_reached = d;
                }
                sum += d;
                reached_sample += 1;
            }
        }
        let mean = if reached_sample > 0 { sum / reached_sample as f64 } else { f64::NAN };

        eprintln!(
            "Phase 1.1 cycle-0 wedge distance stats (64² seed 42, intra-plate, L_decay = {l_decay:.1}, max = {max_d:.1}):"
        );
        eprintln!(
            "  seeds (upper_plate_mask)  = {seeds_total} cells ({:.2} %)",
            100.0 * seeds_total as f64 / total
        );
        eprintln!(
            "  reached (d < max)         = {reached_total} cells ({:.2} %)",
            100.0 * reached_total as f64 / total
        );
        eprintln!(
            "  primary (d < L_decay)     = {primary_total} cells ({:.2} %)",
            100.0 * primary_total as f64 / total
        );
        if reached_sample > 0 {
            eprintln!("  reached cells min/mean/max = {min:.3} / {mean:.3} / {max_reached:.3}");
        } else {
            eprintln!("  no cells reached — upper_plate_mask is empty?");
        }

        // Per-plate breakdown.
        let num_plates = state.num_plates;
        let mut per_plate_seeds = vec![0_usize; num_plates];
        let mut per_plate_reached = vec![0_usize; num_plates];
        let mut per_plate_mean_sum = vec![0.0_f64; num_plates];
        for j in 0..state.ny() {
            for i in 0..state.nx() {
                let pid = state.plate_id.get(i, j) as usize;
                if boundary.upper_plate_mask.get(i, j) {
                    per_plate_seeds[pid] += 1;
                }
                let d = dist.get(i, j);
                if d < max_d {
                    per_plate_reached[pid] += 1;
                    per_plate_mean_sum[pid] += d;
                }
            }
        }

        eprintln!("  per-plate breakdown:");
        let (vx_kin, vy_kin) = (
            kinematics.velocities.iter().map(|(x, _)| *x).collect::<Vec<f64>>(),
            kinematics.velocities.iter().map(|(_, y)| *y).collect::<Vec<f64>>(),
        );
        for p in 0..num_plates {
            let seeds = per_plate_seeds[p];
            let reached = per_plate_reached[p];
            let mean_d =
                if reached > 0 { per_plate_mean_sum[p] / reached as f64 } else { f64::NAN };
            let mag = (vx_kin[p] * vx_kin[p] + vy_kin[p] * vy_kin[p]).sqrt();
            eprintln!(
                "    plate {p}: v=({:>+.4}, {:>+.4}), |v|={:.4} — seeds={:>3}  reached={:>4}  mean_d={:.2}",
                vx_kin[p], vy_kin[p], mag, seeds, reached, mean_d
            );
        }

        // Sanity asserts.
        assert_eq!(
            seeds_total, 111,
            "wedge seed count should match Stage 2 upper_plate_mask count"
        );
        assert!(reached_total > seeds_total, "intra-plate expansion should reach interior cells");
        // At least one plate should be totally silent (no seeds),
        // confirming the symmetric-pair finding.
        let silent_plates = per_plate_seeds.iter().filter(|&&s| s == 0).count();
        assert!(
            silent_plates > 0,
            "expected at least one plate with no upper-plate seeds (symmetric kinematics)"
        );
    }
}
