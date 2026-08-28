//! Step 8.5b Phase 3 — algebraic greedy graph coloring for RBGS.
//!
//! RBGS needs a partition of the unknowns into **colors** such that
//! no two unknowns of the same color share a non-zero off-diagonal
//! entry. Once such a partition exists, all unknowns of the same
//! color can be updated in parallel without data races (no writer
//! thread for color `c` touches a cell that another writer thread
//! for color `c` depends on).
//!
//! # Algorithm (per reviewer Q2 validation)
//!
//! Classical greedy coloring, rows scanned in ascending index
//! order, colors assigned in ascending order:
//!
//! ```text
//! for i = 0 .. n:
//!     used ← { color[j] : j ≠ i, a[i, j] ≠ 0, color[j] defined }
//!     color[i] ← min c ≥ 0 with c ∉ used
//! ```
//!
//! Deterministic by construction — the same matrix always yields
//! the same coloring, regardless of parallelism elsewhere.
//!
//! # Cost
//!
//! `O(nnz)` per matrix. For a 9-point stencil on 4 096 cells this
//! is ~37 000 scalar operations — negligible compared to even one
//! CG matvec at the same level. Coloring is computed once per
//! hierarchy level at build time, never in the inner loop.
//!
//! # Expected number of colors
//!
//! On a structured grid the level-0 Picard block carries a 9-point
//! stencil; the greedy scan then returns a **4-coloring** (a
//! `(i mod 2, j mod 2)` partition would also work — the greedy
//! discovers it). On coarser levels, Galerkin coarsening expands
//! the stencil; the greedy often still returns 4 or 6 colors but
//! nothing guarantees this a priori. Hierarchies are instrumented
//! with [`max_colors_in_hierarchy`] so the report can flag any
//! level with > 4 colors (D2 watch point). > 4 colors does not
//! change correctness; only a mild amount of achievable
//! parallelism within a sweep.

use super::super::sparse_assembly::CsrMatrix;

/// Return the color partition of `a` as a vector of sorted row
/// indices per color. `result[c]` is the list of rows coloured `c`.
///
/// Rows with no off-diagonal non-zero entries (isolated vertices)
/// receive color 0.
pub fn greedy_coloring(a: &CsrMatrix) -> Vec<Vec<usize>> {
    let n = a.n_rows;
    debug_assert_eq!(a.n_cols, n, "greedy_coloring requires a square matrix");

    let mut color = vec![usize::MAX; n];

    // Scratch "used" set reused per row. Small `Vec<bool>` is faster
    // than a `BTreeSet` for the typical 4–8 colors encountered.
    let mut used = Vec::<bool>::new();

    for i in 0..n {
        used.clear();
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let j = a.col_idx[k];
            if j == i {
                continue;
            }
            let cj = color[j];
            if cj != usize::MAX {
                if cj >= used.len() {
                    used.resize(cj + 1, false);
                }
                used[cj] = true;
            }
        }
        let mut c = 0usize;
        while c < used.len() && used[c] {
            c += 1;
        }
        color[i] = c;
    }

    let n_colors = color.iter().copied().filter(|&c| c != usize::MAX).max().unwrap_or(0) + 1;
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_colors];
    for (i, &c) in color.iter().enumerate() {
        if c != usize::MAX {
            groups[c].push(i);
        }
    }
    groups
}

/// Maximum color count across every level of the hierarchy.
///
/// Reported to the Step 8.5b performance log per D2 watch point:
/// any level > 4 is worth noting in the report (classical
/// Ruge-Stüben matrices tend to admit 4-colorings on our
/// geometry; anomalies signal a denser-than-expected coarse
/// operator).
pub fn max_colors_in_hierarchy(levels_colors: &[Vec<Vec<usize>>]) -> usize {
    levels_colors.iter().map(|c| c.len()).max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::super::strong_connections::build_poisson_laplacian_csr;
    use super::*;

    fn validate_coloring(a: &CsrMatrix, colors: &[Vec<usize>]) {
        let n = a.n_rows;
        let mut membership = vec![usize::MAX; n];
        for (c, grp) in colors.iter().enumerate() {
            for &i in grp {
                assert_eq!(membership[i], usize::MAX, "row {i} assigned to two colors",);
                membership[i] = c;
            }
        }
        for i in 0..n {
            assert_ne!(membership[i], usize::MAX, "row {i} missing from partition");
        }
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k];
                if j == i {
                    continue;
                }
                assert_ne!(
                    membership[i], membership[j],
                    "row {i} and its non-zero neighbour {j} share colour {} \
                     — greedy invariant violated",
                    membership[i],
                );
            }
        }
    }

    #[test]
    fn poisson_5pt_colorable_with_2_colors() {
        // 5-point Laplacian: (i, j) connects to (i±1, j) and (i, j±1).
        // Classical checkerboard gives exactly 2 colors; greedy
        // naturally recovers it because rows 0, 2, 4, ... are all
        // mutually non-adjacent.
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let colors = greedy_coloring(&a);
        validate_coloring(&a, &colors);
        assert_eq!(colors.len(), 2, "5-point Poisson expected 2 colors, got {}", colors.len());
    }

    #[test]
    fn isolated_vertex_gets_color_zero() {
        let a = CsrMatrix {
            n_rows: 3,
            n_cols: 3,
            row_ptr: vec![0, 1, 2, 3],
            col_idx: vec![0, 1, 2],
            values: vec![1.0, 2.0, 3.0],
        };
        let colors = greedy_coloring(&a);
        validate_coloring(&a, &colors);
        assert_eq!(colors.len(), 1);
        assert_eq!(colors[0], vec![0, 1, 2]);
    }

    #[test]
    fn coloring_is_deterministic() {
        let a = build_poisson_laplacian_csr(6);
        let c1 = greedy_coloring(&a);
        for _ in 0..50 {
            let c2 = greedy_coloring(&a);
            assert_eq!(c1.len(), c2.len());
            for (a, b) in c1.iter().zip(c2.iter()) {
                assert_eq!(a, b);
            }
        }
    }

    #[test]
    fn coloring_handles_dense_row() {
        // Row 0 is fully connected to everything else — greedy will
        // put it alone in color 0 and place everything else in color 1.
        let n = 5;
        let mut row_ptr = vec![0];
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        // Row 0: connects to 0, 1, 2, 3, 4
        col_idx.extend_from_slice(&[0, 1, 2, 3, 4]);
        values.extend_from_slice(&[5.0, 1.0, 1.0, 1.0, 1.0]);
        row_ptr.push(col_idx.len());
        // Rows 1..5: diagonal only + connection to 0
        for i in 1..n {
            col_idx.extend_from_slice(&[0, i]);
            values.extend_from_slice(&[1.0, 2.0]);
            row_ptr.push(col_idx.len());
        }
        let a = CsrMatrix { n_rows: n, n_cols: n, row_ptr, col_idx, values };
        let colors = greedy_coloring(&a);
        validate_coloring(&a, &colors);
        assert_eq!(colors.len(), 2);
        assert_eq!(colors[0], vec![0]);
        assert_eq!(colors[1], vec![1, 2, 3, 4]);
    }

    #[test]
    fn groups_are_index_sorted() {
        // Greedy iterates i in ascending order and appends to the
        // group that wins. So each group[c] must already be sorted.
        let a = build_poisson_laplacian_csr(8);
        let colors = greedy_coloring(&a);
        for grp in &colors {
            assert!(grp.windows(2).all(|w| w[0] < w[1]), "group not sorted: {grp:?}");
        }
    }
}
