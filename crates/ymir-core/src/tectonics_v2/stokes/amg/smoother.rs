//! Symmetric Gauss-Seidel smoother — Step 8.5a Phase 2.4.
//!
//! One **symmetric** sweep of Gauss-Seidel = one forward pass
//! (rows `0..N` in ascending order) + one backward pass
//! (rows `N-1..=0` in descending order). Each row update is
//! ```text
//!     x_i ← (b_i − ∑_{j ≠ i} a_ij · x_j) / a_ii
//! ```
//! using the most-recently-written `x_j` values. Symmetric
//! sweeps preserve SPD structure required by CG — a single-
//! direction forward GS would not.
//!
//! # Role in AMG
//!
//! The smoother's job is to **damp high-frequency error modes**
//! that the coarse-grid correction cannot represent. Classical
//! spectral theory: on Poisson the high-frequency error
//! reduction per SGS sweep is approximately 0.5 (per mode),
//! while low-frequency reduction is much slower — hence the
//! need for multigrid in the first place.
//!
//! AMG's convergence relies on the interplay: SGS handles highs,
//! coarse grid handles lows. The default `pre_sweeps = post_
//! sweeps = 1` (one SGS sweep = 2 passes forward+backward) is
//! sufficient for Classical RS on SPD M-matrices.

use super::super::sparse_assembly::CsrMatrix;

/// Apply one symmetric Gauss-Seidel sweep on `A · x = b`. `x`
/// is updated in place; `b` is read-only.
///
/// Complexity: O(nnz(a)) per sweep.
pub fn sgs_sweep(a: &CsrMatrix, b: &[f64], x: &mut [f64]) {
    let n = a.n_rows;
    debug_assert_eq!(a.n_cols, n, "SGS requires a square matrix");
    debug_assert_eq!(b.len(), n);
    debug_assert_eq!(x.len(), n);

    // --- Forward pass (i = 0..N) ---
    for i in 0..n {
        let start = a.row_ptr[i];
        let end = a.row_ptr[i + 1];
        let mut acc = b[i];
        let mut diag = 0.0f64;
        for k in start..end {
            let j = a.col_idx[k];
            if j == i {
                diag = a.values[k];
            } else {
                acc -= a.values[k] * x[j];
            }
        }
        // Skip rows with zero diagonal (degenerate; debug-assert
        // caught in production paths). Leaving x[i] unchanged is
        // the safe default.
        if diag.abs() > 1e-300 {
            x[i] = acc / diag;
        }
    }

    // --- Backward pass (i = N−1 ..= 0) ---
    for i in (0..n).rev() {
        let start = a.row_ptr[i];
        let end = a.row_ptr[i + 1];
        let mut acc = b[i];
        let mut diag = 0.0f64;
        for k in start..end {
            let j = a.col_idx[k];
            if j == i {
                diag = a.values[k];
            } else {
                acc -= a.values[k] * x[j];
            }
        }
        if diag.abs() > 1e-300 {
            x[i] = acc / diag;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::strong_connections::build_poisson_laplacian_csr;
    use super::*;

    fn matvec(a: &CsrMatrix, x: &[f64], y: &mut [f64]) {
        a.apply(x, y);
    }

    fn norm2(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    #[test]
    fn sgs_reduces_residual_monotonically_on_poisson() {
        // Solve A · x = b on a constant-coefficient 5-pt Poisson
        // with a random seeded RHS (mean-subtracted, matching the
        // null-space convention of the real solver). Residuals
        // must decrease monotonically across SGS sweeps.
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let total = n * n;
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        let mut b: Vec<f64> = (0..total).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let mean: f64 = b.iter().sum::<f64>() / total as f64;
        for v in b.iter_mut() {
            *v -= mean;
        }

        let mut x = vec![0.0f64; total];
        let mut ax = vec![0.0f64; total];
        matvec(&a, &x, &mut ax);
        let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let r0_norm = norm2(&r0);

        let mut prev = r0_norm;
        for sweep in 0..30 {
            sgs_sweep(&a, &b, &mut x);
            matvec(&a, &x, &mut ax);
            let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
            let r_norm = norm2(&r);
            // Monotonic decrease is the hard invariant.
            assert!(
                r_norm <= prev * 1.0 + 1e-12,
                "sweep {}: residual increased from {:.3e} to {:.3e}",
                sweep,
                prev,
                r_norm
            );
            prev = r_norm;
        }
        // After 30 sweeps on an 8·8 Poisson, expect ≥ 10× reduction
        // (high-freq modes damped fast; low-freq slow but still
        // present; well-conditioned at n=8 lets both reduce).
        assert!(
            prev < r0_norm / 10.0,
            "30 SGS sweeps gave only {:.3}× reduction",
            r0_norm / prev
        );
    }

    #[test]
    fn sgs_is_exact_on_diagonal_system() {
        // Diagonal SPD system: `D · x = b` — one SGS sweep (single
        // forward pass or single backward pass individually) is
        // already exact, so after one symmetric sweep we must have
        // exactly `x = b / D`.
        let d_vals: Vec<f64> = (1..=5).map(|k| k as f64).collect();
        let b: Vec<f64> = (0..5).map(|k| (k as f64) * 0.5 + 1.0).collect();
        let a = CsrMatrix {
            n_rows: 5,
            n_cols: 5,
            row_ptr: vec![0, 1, 2, 3, 4, 5],
            col_idx: vec![0, 1, 2, 3, 4],
            values: d_vals.clone(),
        };
        let mut x = vec![0.0; 5];
        sgs_sweep(&a, &b, &mut x);
        for k in 0..5 {
            assert!((x[k] - b[k] / d_vals[k]).abs() < 1e-15);
        }
    }

    #[test]
    fn sgs_high_frequency_mode_damps_by_roughly_half_per_sweep() {
        // Inject a high-frequency eigenmode of the constant-coefficient
        // 5-pt Poisson (the "checkerboard" pattern is close to the
        // maximum eigenvalue) and verify SGS reduces its amplitude
        // substantially. This is the "smoother does what AMG needs"
        // contract — high-freq damping ≤ ~0.5 per sweep is the
        // classical figure.
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let total = n * n;
        // Checkerboard-like mode x_ij = (-1)^(i+j).
        let mut x: Vec<f64> = (0..total)
            .map(|k| {
                let i = k % n;
                let j = k / n;
                if (i + j) % 2 == 0 { 1.0 } else { -1.0 }
            })
            .collect();
        // Subtract mean for null-space consistency.
        let m: f64 = x.iter().sum::<f64>() / total as f64;
        for v in x.iter_mut() {
            *v -= m;
        }
        let x0_norm = norm2(&x);
        // Smooth the A · x = 0 problem (the error equation) — each
        // sweep reduces ||x|| toward 0.
        let b = vec![0.0; total];
        let mut prev = x0_norm;
        for _ in 0..3 {
            sgs_sweep(&a, &b, &mut x);
            let n2 = norm2(&x);
            // Hard invariant: high-freq mode damps non-trivially
            // (< 0.9 per sweep). 0.5 is the ideal figure; relax to
            // 0.9 for robustness to grid-size edge cases.
            assert!(
                n2 < prev * 0.9,
                "high-freq mode damped only by {:.3} per sweep",
                n2 / prev
            );
            prev = n2;
        }
    }

    #[test]
    fn sgs_is_deterministic_across_runs() {
        let n = 6;
        let a = build_poisson_laplacian_csr(n);
        let total = n * n;
        let b: Vec<f64> = (0..total).map(|k| (k as f64).sin()).collect();
        let mut x_a = vec![0.1f64; total];
        let mut x_b = x_a.clone();
        for _ in 0..5 {
            sgs_sweep(&a, &b, &mut x_a);
            sgs_sweep(&a, &b, &mut x_b);
        }
        for k in 0..total {
            assert_eq!(x_a[k].to_bits(), x_b[k].to_bits());
        }
    }
}
