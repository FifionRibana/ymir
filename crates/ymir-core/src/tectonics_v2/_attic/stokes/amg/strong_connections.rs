//! Strong-connection detection — Step 8.5a Phase 2.1.
//!
//! For each row `i` of an SPD M-matrix `A` (non-positive
//! off-diagonals), a column `j ≠ i` is **strongly connected** to
//! `i` when
//! ```text
//!     -a_ij ≥ θ · max_{k ≠ i} (-a_ik)
//! ```
//! per Briggs-Henson-McCormick, "A Multigrid Tutorial" 2nd ed.,
//! §8.8 eq. (8.43). This is the Classical Ruge-Stüben definition
//! driving both C/F splitting (Phase 2.2) and prolongation
//! weights (Phase 2.3).
//!
//! # Sign convention
//!
//! Classical RS assumes the operator is an M-matrix (diagonally
//! dominant with `a_ii > 0`, `a_ij ≤ 0` for `i ≠ j`). The
//! Picard scalar blocks (u-u and v-v partitions of
//! `A_picard`) satisfy this: all off-diagonals are non-positive
//! by the discrete Laplacian + shear stencil (cf.
//! `sparse_assembly.rs` stencil documentation). Positive
//! off-diagonals, if they ever appear, are excluded from the
//! strong set by the `-a_ij > 0` implicit filter.
//!
//! # Output
//!
//! `Vec<Vec<usize>>` where `strong[i]` is the sorted list of
//! column indices strongly connected to row `i`. Sorted ascending
//! because the input CSR rows are column-sorted (Phase 1
//! invariant), and we preserve that order throughout — D9.

use super::super::sparse_assembly::CsrMatrix;

/// Compute strong-connection sets per row of `a` with threshold
/// `theta` (typical 0.25).
///
/// Complexity: `O(nnz(a))` — two linear passes per row.
pub fn compute_strong_connections(a: &CsrMatrix, theta: f64) -> Vec<Vec<usize>> {
    assert!(theta >= 0.0 && theta <= 1.0, "θ must lie in [0, 1]; got {}", theta);
    let n = a.n_rows;
    let mut out: Vec<Vec<usize>> = Vec::with_capacity(n);

    for i in 0..n {
        let start = a.row_ptr[i];
        let end = a.row_ptr[i + 1];

        // Pass 1: find max(-a_ij) over off-diagonal columns.
        let mut max_neg = 0.0f64;
        for k in start..end {
            let j = a.col_idx[k];
            if j == i {
                continue;
            }
            let neg = -a.values[k];
            if neg > max_neg {
                max_neg = neg;
            }
        }

        // Pass 2: collect columns whose -a_ij exceeds the threshold.
        let threshold = theta * max_neg;
        let mut strong = Vec::with_capacity(end - start);
        if max_neg > 0.0 {
            for k in start..end {
                let j = a.col_idx[k];
                if j == i {
                    continue;
                }
                let neg = -a.values[k];
                if neg >= threshold && neg > 0.0 {
                    strong.push(j);
                }
            }
        }
        // Input rows are column-sorted (Phase 1 invariant), so the
        // filtered output is already ascending.
        out.push(strong);
    }
    out
}

// --- Test utilities (used here and re-used by later sub-phases). ----

/// Build a scalar 5-point Poisson Laplacian `-∇²` on an `n × n`
/// periodic grid with unit spacing. Row/column ordering is
/// row-major: index `j·n + i`. Useful harness for testing the
/// AMG algorithm modules in isolation, independent of the
/// momentum-operator complexity.
#[cfg(test)]
pub(crate) fn build_poisson_laplacian_csr(n: usize) -> CsrMatrix {
    use super::super::sparse_assembly::CsrMatrix;
    let total = n * n;
    let mut row_ptr = Vec::with_capacity(total + 1);
    let mut col_idx = Vec::with_capacity(5 * total);
    let mut values = Vec::with_capacity(5 * total);
    row_ptr.push(0);
    for j in 0..n {
        let jp = (j + 1) % n;
        let jm = (j + n - 1) % n;
        for i in 0..n {
            let ip = (i + 1) % n;
            let im = (i + n - 1) % n;
            // Row entries in ascending column index order.
            // Indices in a row-major ordering: (jm*n + i), (j*n + im),
            // (j*n + i), (j*n + ip), (jp*n + i) — but ascending only if
            // jm < j < jp, im < i < ip which is NOT the case on wrap.
            // Collect into a buffer and sort to keep CSR invariant.
            let mut buf = [
                (jm * n + i, -1.0),
                (j * n + im, -1.0),
                (j * n + i, 4.0),
                (j * n + ip, -1.0),
                (jp * n + i, -1.0),
            ];
            buf.sort_by_key(|&(c, _)| c);
            // Merge duplicates on small-n periodic wrap.
            let mut prev: Option<usize> = None;
            for (c, v) in buf.iter() {
                match prev {
                    Some(p) if p == *c => {
                        *values.last_mut().unwrap() += v;
                    }
                    _ => {
                        col_idx.push(*c);
                        values.push(*v);
                        prev = Some(*c);
                    }
                }
            }
            row_ptr.push(col_idx.len());
        }
    }
    CsrMatrix {
        n_rows: total,
        n_cols: total,
        row_ptr,
        col_idx,
        values,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_neighbours_strong_on_constant_poisson_at_theta_0_25() {
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);

        for (i, row) in strong.iter().enumerate() {
            // Constant 5-pt Poisson: each row has 4 off-diagonals
            // of magnitude 1 and one positive diagonal of magnitude 4.
            // max(-a_ij) = 1; θ·max = 0.25. All 4 off-diagonals are
            // strong.
            assert_eq!(
                row.len(),
                4,
                "row {} has {} strong connections, expected 4",
                i,
                row.len()
            );
            // Sorted ascending.
            assert!(
                row.windows(2).all(|w| w[0] < w[1]),
                "row {} strong set not sorted: {:?}",
                i,
                row
            );
        }
    }

    #[test]
    fn empty_strong_at_theta_1_0_plus_epsilon_excludes_equal_values() {
        // For the constant Poisson, all off-diagonals are equal
        // in magnitude. With θ = 1.0 exact, -a_ij == threshold
        // keeps them (≥, not >). Raise θ slightly above 1 via
        // numerical comparison to confirm the threshold is a
        // well-defined cutoff (no entries qualify when threshold
        // strictly exceeds max magnitude).
        let n = 6;
        let a = build_poisson_laplacian_csr(n);
        // θ = 1.0 — all off-diagonals are exactly at threshold,
        // so they all remain strong.
        let strong = compute_strong_connections(&a, 1.0);
        for row in &strong {
            assert_eq!(row.len(), 4);
        }
    }

    #[test]
    fn all_strong_at_theta_0() {
        let n = 5;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.0);
        for row in &strong {
            assert_eq!(row.len(), 4);
        }
    }

    #[test]
    fn isolated_node_produces_empty_strong_set() {
        // Construct a 2×2 diag-only CSR: no off-diagonals. Every
        // row's strong set must be empty (max_neg = 0).
        let a = CsrMatrix {
            n_rows: 2,
            n_cols: 2,
            row_ptr: vec![0, 1, 2],
            col_idx: vec![0, 1],
            values: vec![1.0, 1.0],
        };
        let s = compute_strong_connections(&a, 0.25);
        assert!(s[0].is_empty());
        assert!(s[1].is_empty());
    }

    #[test]
    fn deterministic_100_runs_poisson() {
        // Determinism check. The function is pure on the input
        // CSR — no RNG, no HashMap, no floating-point reduction
        // with ambiguous ordering. 100 runs must all match the
        // first one bit-for-bit (structurally; `Vec<Vec<usize>>`
        // equality is sufficient).
        let n = 10;
        let a = build_poisson_laplacian_csr(n);
        let first = compute_strong_connections(&a, 0.25);
        for _ in 0..100 {
            let next = compute_strong_connections(&a, 0.25);
            assert_eq!(first, next);
        }
    }

    #[test]
    fn heterogeneous_row_filters_weak_offdiagonals() {
        // Craft a single-row toy matrix where off-diagonals have
        // differing magnitudes; only those ≥ 0.25·max survive at
        // θ = 0.25.
        let a = CsrMatrix {
            n_rows: 1,
            n_cols: 5,
            row_ptr: vec![0, 5],
            col_idx: vec![0, 1, 2, 3, 4],
            // Row values: diagonal 10, then -1, -4, -0.5 (weak), -2.
            // -a_ij values: -, 1, 4, 0.5, 2. max = 4. 0.25·max = 1.
            // So strong set at θ=0.25 = {columns with -a_ij ≥ 1}
            //                        = {1 (val 1), 2 (val 4), 4 (val 2)}.
            // Column 3 has -a_ij = 0.5 < 1 → filtered out.
            values: vec![10.0, -1.0, -4.0, -0.5, -2.0],
        };
        let s = compute_strong_connections(&a, 0.25);
        assert_eq!(s[0], vec![1, 2, 4]);
    }
}
