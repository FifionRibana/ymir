//! Dense LU (Doolittle with partial pivoting) for the coarsest
//! grid — Step 8.5a Phase 2.5.
//!
//! The coarsest level of an AMG hierarchy has ≤ `min_coarse_
//! unknowns` (default 50) entries; factorising it via a sparse
//! direct method would be overkill. The dense O(n³) Doolittle
//! factorisation runs in ~125 k operations for n = 50 and is
//! below the noise floor of a V-cycle's total cost.
//!
//! Reimplemented here instead of pulling in `nalgebra` — per
//! QA3 resolution, no new workspace dependency is added solely
//! for the coarse solve.
//!
//! # Algorithm
//!
//! `A = P · L · U` where L is unit lower triangular, U is upper
//! triangular, P is row permutation (partial pivoting).
//!
//! ```text
//! for k in 0..n:
//!     p = argmax_{i ≥ k} |A[i, k]|   (ties → lowest index)
//!     swap rows k and p (tracked via `perm`)
//!     for i in (k+1)..n:
//!         A[i, k] /= A[k, k]           (= L[i, k])
//!         for j in (k+1)..n:
//!             A[i, j] -= A[i, k] · A[k, j]
//! ```
//!
//! Back-substitution and forward-substitution are standard.
//!
//! # Determinism (D9)
//!
//! - Pivot ties by **lowest row index** (reviewer's vigilance
//!   point 2 echoed at the coarsest level).
//! - In-place dense factorisation is inherently deterministic;
//!   no PRNG, no HashMap, no unstable reduction.

use super::super::sparse_assembly::CsrMatrix;

/// Factorisation of a small dense SPD matrix via Doolittle LU
/// with partial pivoting.
#[derive(Clone, Debug)]
pub struct LuFactorisation {
    pub n: usize,
    /// Combined L (strict lower) + U (upper incl diag), row-major.
    pub lu: Vec<f64>,
    /// Row-permutation from partial pivoting; `perm[i]` is the
    /// original row index now sitting at row `i`.
    pub perm: Vec<usize>,
}

impl LuFactorisation {
    /// Factorise the CSR matrix `a` into in-place Doolittle LU.
    ///
    /// Panics if `a` is not square or is rank-deficient (a zero
    /// pivot column means the operator is singular even after
    /// permutation).
    pub fn factor(a: &CsrMatrix) -> Self {
        let n = a.n_rows;
        assert_eq!(a.n_cols, n, "coarse_solve requires a square matrix");

        // --- Densify ---
        let mut lu = vec![0.0f64; n * n];
        for i in 0..n {
            let start = a.row_ptr[i];
            let end = a.row_ptr[i + 1];
            for k in start..end {
                let j = a.col_idx[k];
                lu[i * n + j] = a.values[k];
            }
        }

        // --- Doolittle with partial pivoting ---
        let mut perm: Vec<usize> = (0..n).collect();
        for k in 0..n {
            // Pivot selection: largest |A[i, k]| for i ≥ k; ties
            // by lowest i.
            let mut best_row = k;
            let mut best_val = lu[k * n + k].abs();
            for i in (k + 1)..n {
                let v = lu[i * n + k].abs();
                if v > best_val {
                    best_val = v;
                    best_row = i;
                }
            }
            if best_val < 1e-300 {
                panic!(
                    "coarse-grid operator is singular at column {} (max |pivot| = {:.3e})",
                    k, best_val
                );
            }
            // Swap rows k and best_row in both lu and perm.
            if best_row != k {
                for j in 0..n {
                    lu.swap(k * n + j, best_row * n + j);
                }
                perm.swap(k, best_row);
            }

            // Elimination.
            let pivot = lu[k * n + k];
            for i in (k + 1)..n {
                let factor = lu[i * n + k] / pivot;
                lu[i * n + k] = factor; // store L[i, k]
                for j in (k + 1)..n {
                    let l_ik = factor;
                    let u_kj = lu[k * n + j];
                    lu[i * n + j] -= l_ik * u_kj;
                }
            }
        }

        Self { n, lu, perm }
    }

    /// Solve `A · x = b` using the stored factorisation.
    ///
    /// `b` is read-only; `x` is written. Both slices must have
    /// length `self.n`.
    pub fn solve(&self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        debug_assert_eq!(b.len(), n);
        debug_assert_eq!(x.len(), n);

        // Apply permutation: y ← P · b via perm.
        let mut y = vec![0.0f64; n];
        for i in 0..n {
            y[i] = b[self.perm[i]];
        }

        // Forward substitution: L · z = y. L has unit diagonal
        // (implicit, not stored).
        for i in 0..n {
            let mut acc = y[i];
            for j in 0..i {
                acc -= self.lu[i * n + j] * y[j];
            }
            y[i] = acc;
        }

        // Back substitution: U · x = z.
        for i in (0..n).rev() {
            let mut acc = y[i];
            for j in (i + 1)..n {
                acc -= self.lu[i * n + j] * x[j];
            }
            x[i] = acc / self.lu[i * n + i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a small dense CSR from a row-major `Vec<f64>`.
    fn csr_from_dense(n: usize, dense: &[f64]) -> CsrMatrix {
        let mut row_ptr = vec![0usize; n + 1];
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let v = dense[i * n + j];
                if v != 0.0 {
                    col_idx.push(j);
                    values.push(v);
                }
            }
            row_ptr[i + 1] = col_idx.len();
        }
        CsrMatrix { n_rows: n, n_cols: n, row_ptr, col_idx, values }
    }

    #[test]
    fn lu_solves_2x2_spd_exactly() {
        // A = [[4, 1], [1, 3]],  b = [1, 2]   → x = [1/11, 7/11]
        let a = csr_from_dense(2, &[4.0, 1.0, 1.0, 3.0]);
        let lu = LuFactorisation::factor(&a);
        let mut x = vec![0.0; 2];
        lu.solve(&[1.0, 2.0], &mut x);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-14);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-14);
    }

    #[test]
    fn lu_solves_diagonal_system() {
        let a = csr_from_dense(
            4,
            &[2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 7.0],
        );
        let lu = LuFactorisation::factor(&a);
        let b = [4.0, 9.0, 25.0, 49.0];
        let mut x = vec![0.0; 4];
        lu.solve(&b, &mut x);
        assert_eq!(x, vec![2.0, 3.0, 5.0, 7.0]);
    }

    #[test]
    fn lu_is_deterministic_across_factorisations() {
        // Non-trivial 4x4 SPD (shifted Hilbert-like). Factorise
        // 10 times, every LU must be bitwise identical.
        let mut dense = Vec::with_capacity(16);
        for i in 0..4 {
            for j in 0..4 {
                let v = 1.0 / (1.0 + (i + j) as f64) + if i == j { 2.0 } else { 0.0 };
                dense.push(v);
            }
        }
        let a = csr_from_dense(4, &dense);
        let first = LuFactorisation::factor(&a);
        for _ in 0..10 {
            let next = LuFactorisation::factor(&a);
            assert_eq!(first.n, next.n);
            assert_eq!(first.perm, next.perm);
            for k in 0..first.lu.len() {
                assert_eq!(first.lu[k].to_bits(), next.lu[k].to_bits());
            }
        }
    }

    #[test]
    fn lu_partial_pivoting_selects_lowest_index_on_tie() {
        // Matrix where rows 0 and 2 have the same magnitude in
        // column 0: |A[0,0]| = |A[2,0]| = 2. Partial pivoting
        // must leave row 0 at position 0 (lowest-index tiebreak),
        // not swap with row 2.
        let a = csr_from_dense(3, &[2.0, 1.0, 1.0, 0.0, 3.0, 2.0, 2.0, 0.0, 4.0]);
        let lu = LuFactorisation::factor(&a);
        // perm[0] should still be 0 — no tie-induced swap.
        assert_eq!(lu.perm[0], 0);
    }

    #[test]
    fn lu_solve_roundtrips_random_spd() {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let n = 12;
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        // Build random M = C·Cᵀ (symmetric positive-definite by
        // construction). Use small `C` magnitudes to keep condition
        // number tame.
        let mut c = vec![0.0f64; n * n];
        for v in c.iter_mut() {
            *v = rng.random::<f64>() * 0.5 - 0.25;
        }
        let mut dense = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0;
                for k in 0..n {
                    acc += c[i * n + k] * c[j * n + k];
                }
                dense[i * n + j] = acc + if i == j { 1.0 } else { 0.0 };
            }
        }
        let a = csr_from_dense(n, &dense);
        let lu = LuFactorisation::factor(&a);

        // Solve A · x = b; then verify A · x = b to 1e-12 relative.
        let b: Vec<f64> = (0..n).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let mut x = vec![0.0; n];
        lu.solve(&b, &mut x);
        // Reference: compute A · x directly from dense.
        let mut ax = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                ax[i] += dense[i * n + j] * x[j];
            }
        }
        for k in 0..n {
            let rel = (ax[k] - b[k]).abs() / b[k].abs().max(1e-10);
            assert!(rel < 1e-12, "A·x[{}] = {:.6e}, b[{}] = {:.6e}", k, ax[k], k, b[k]);
        }
    }
}
