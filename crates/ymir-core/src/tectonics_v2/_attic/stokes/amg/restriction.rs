//! Restriction matrix construction — Step 8.5a Phase 2.3.
//!
//! For Galerkin coarsening, `R = Pᵀ`. This preserves the
//! variational property `A_coarse = R · A_fine · P` and keeps the
//! coarse operator SPD if the fine-grid one is, which is a
//! precondition for symmetric-Gauss-Seidel on each level.
//!
//! The transpose is done in two passes, O(nnz(P)):
//! 1. Count column occurrences of `p.col_idx` to build
//!    `r.row_ptr` by prefix sum.
//! 2. Fill `r.col_idx` / `r.values` by walking `p` and placing
//!    each `(i, k, v)` entry into the current write head of
//!    row `k` in `r`.
//!
//! Output invariant: rows of `r` are in strictly ascending
//! column-index order (= ascending i-index of the original fine
//! row), inherited from `p`'s row-major traversal.

use super::super::sparse_assembly::CsrMatrix;

pub fn transpose_to_restriction(p: &CsrMatrix) -> CsrMatrix {
    let n_fine = p.n_rows;
    let n_coarse = p.n_cols;
    let nnz = p.values.len();

    // --- Pass 1: column count ---
    let mut row_counts = vec![0usize; n_coarse];
    for &c in &p.col_idx {
        row_counts[c] += 1;
    }

    // Prefix-sum into row_ptr of R.
    let mut row_ptr = vec![0usize; n_coarse + 1];
    for k in 0..n_coarse {
        row_ptr[k + 1] = row_ptr[k] + row_counts[k];
    }
    debug_assert_eq!(row_ptr[n_coarse], nnz);

    // --- Pass 2: fill ---
    let mut col_idx = vec![0usize; nnz];
    let mut values = vec![0.0f64; nnz];
    let mut write_heads = row_ptr[..n_coarse].to_vec();

    for i in 0..n_fine {
        let start = p.row_ptr[i];
        let end = p.row_ptr[i + 1];
        for k in start..end {
            let col = p.col_idx[k];
            let val = p.values[k];
            let dest = write_heads[col];
            col_idx[dest] = i;
            values[dest] = val;
            write_heads[col] += 1;
        }
    }

    CsrMatrix {
        n_rows: n_coarse,
        n_cols: n_fine,
        row_ptr,
        col_idx,
        values,
    }
}

#[cfg(test)]
mod tests {
    use super::super::prolongation::build_prolongation;
    use super::super::splitting::classical_rs_splitting;
    use super::super::strong_connections::{
        build_poisson_laplacian_csr, compute_strong_connections,
    };
    use super::*;

    fn apply_csr_t(p: &CsrMatrix, x: &[f64], y: &mut [f64]) {
        // Reference: compute `y = Pᵀ · x` by walking P row-wise and
        // scattering instead of constructing the transpose.
        for v in y.iter_mut() {
            *v = 0.0;
        }
        for i in 0..p.n_rows {
            let start = p.row_ptr[i];
            let end = p.row_ptr[i + 1];
            for k in start..end {
                y[p.col_idx[k]] += p.values[k] * x[i];
            }
        }
    }

    #[test]
    fn restriction_matches_pt_product() {
        let n = 10;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let p = build_prolongation(&a, &strong, &cf);
        let r = transpose_to_restriction(&p);
        assert_eq!(r.n_rows, p.n_cols);
        assert_eq!(r.n_cols, p.n_rows);

        // Compare R·x with the reference Pᵀ·x for seeded random
        // vectors. Bitwise equality expected because both compute
        // the same sum-of-products in the same column-major order.
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let mut rng = ChaCha8Rng::seed_from_u64(7);
        let x: Vec<f64> = (0..r.n_cols).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let mut y_r = vec![0.0; r.n_rows];
        r.apply(&x, &mut y_r);
        let mut y_pt = vec![0.0; r.n_rows];
        apply_csr_t(&p, &x, &mut y_pt);
        for k in 0..y_r.len() {
            let rel = (y_r[k] - y_pt[k]).abs() / y_pt[k].abs().max(1e-300);
            assert!(
                rel < 1e-14,
                "R·x[{}] = {:.6e}, Pᵀ·x[{}] = {:.6e}",
                k,
                y_r[k],
                k,
                y_pt[k]
            );
        }
    }

    #[test]
    fn restriction_is_column_sorted() {
        let n = 6;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let p = build_prolongation(&a, &strong, &cf);
        let r = transpose_to_restriction(&p);

        for i in 0..r.n_rows {
            let start = r.row_ptr[i];
            let end = r.row_ptr[i + 1];
            if end - start <= 1 {
                continue;
            }
            for k in start + 1..end {
                assert!(
                    r.col_idx[k - 1] < r.col_idx[k],
                    "R row {} not column-sorted: {:?}",
                    i,
                    &r.col_idx[start..end]
                );
            }
        }
    }
}
