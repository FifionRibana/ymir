//! Prolongation matrix construction — Step 8.5a Phase 2.3.
//!
//! Builds `P : ℝ^{n_coarse} → ℝ^{n_fine}` so that
//! `v_fine ← P · v_coarse` interpolates coarse-grid corrections
//! back to the fine grid. The mapping is the Ruge-Stüben classical
//! formula (Briggs-Henson-McCormick §8.8.3, eq. 8.45–8.48), which
//! routes each F-point's interpolation through its strong
//! C-neighbours `C_i = {k ∈ S_i : cf[k] = Coarse}` and applies a
//! secondary correction through the strong F-F dependencies
//! `F_i^s = {j ∈ S_i : cf[j] = Fine}` that share at least one
//! C-neighbour with `i` (Pass 2 of the splitting guarantees this).
//!
//! # Formula
//!
//! For a Coarse point `i`, `P[i, c(i)] = 1` (direct injection,
//! where `c(i)` is `i`'s position in the consecutive coarse
//! numbering).
//!
//! For a Fine point `i`,
//! ```text
//!   P[i, c(k)] = -w_{ik} / d_i         for k ∈ C_i
//!
//!   w_{ik} = a_ik + ∑_{j ∈ F_i^s}  a_ij · a_jk / Σ_j
//!   Σ_j    = ∑_{l ∈ C_i ∩ S_j} a_jl         (F-F common-C sum)
//!   d_i    = a_ii + ∑_{m ∈ N_i^w} a_im      (diagonal + weak off-diagonals)
//! ```
//! where `N_i^w = {m ≠ i : a_im ≠ 0, m ∉ S_i ∪ {i}}` collects the
//! weak off-diagonal indices — their contribution is absorbed
//! into the effective diagonal (eq. 8.47), following the Schur
//! interpretation of "coarse-grid approximation of the local row
//! equation".
//!
//! # Determinism
//!
//! All per-row sums are accumulated in ascending column-index
//! order (inherits the CSR sort invariant from Phase 1). No PRNG,
//! no floating-point reduction with ambiguous order.

use super::super::sparse_assembly::CsrMatrix;
use super::splitting::CfType;

/// Build the prolongation CSR from the fine-grid operator
/// `a`, the strong-connection structure `strong`, and the C/F
/// labelling `cf`.
///
/// Output dimensions: `n_rows = cf.len()` (fine), `n_cols = |C|`
/// (coarse).
pub fn build_prolongation(a: &CsrMatrix, strong: &[Vec<usize>], cf: &[CfType]) -> CsrMatrix {
    let n = cf.len();
    assert_eq!(a.n_rows, n);
    assert_eq!(a.n_cols, n);
    assert_eq!(strong.len(), n);

    // --- Coarse numbering: c(i) = consecutive index among C-points
    //     in ascending i order. `coarse_of[i] = Some(c)` if C,
    //     else None. ---
    let mut coarse_of: Vec<Option<usize>> = vec![None; n];
    let mut n_coarse = 0;
    for i in 0..n {
        if cf[i] == CfType::Coarse {
            coarse_of[i] = Some(n_coarse);
            n_coarse += 1;
        }
    }

    // --- Build P row-by-row. Each row has at most |C_i| + 1
    //     entries in ascending column (= ascending c) order. ---
    let mut row_ptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut col_idx: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    row_ptr.push(0);

    // Helper: look up a_jk in row j of `a`, returning 0 if absent.
    let lookup = |j: usize, k: usize| -> f64 {
        let start = a.row_ptr[j];
        let end = a.row_ptr[j + 1];
        let slice = &a.col_idx[start..end];
        match slice.binary_search(&k) {
            Ok(offset) => a.values[start + offset],
            Err(_) => 0.0,
        }
    };

    for i in 0..n {
        match cf[i] {
            CfType::Coarse => {
                // Identity injection at its coarse slot.
                let c = coarse_of[i].expect("C-point must have coarse index");
                col_idx.push(c);
                values.push(1.0);
            }
            CfType::Fine => {
                // Classify i's strong neighbours into C_i and F_i^s.
                // Sorted ascending (inherits strong sort order).
                let mut c_nbrs: Vec<usize> = Vec::new();
                let mut f_nbrs: Vec<usize> = Vec::new();
                for &j in &strong[i] {
                    match cf[j] {
                        CfType::Coarse => c_nbrs.push(j),
                        CfType::Fine => f_nbrs.push(j),
                        CfType::Undecided => unreachable!(),
                    }
                }
                // Guarded by splitting Pass 2: F-points must have ≥ 1
                // strong C-neighbour. If not, the operator row has no
                // interpolation sources — fall back to identity-zero
                // (v_fine(i) = 0 contribution). This is a defensive
                // path; the splitting test suite ensures it's unused.
                if c_nbrs.is_empty() {
                    // Empty row in P — Pass 2 should prevent this.
                    // Debug-assert to catch regressions.
                    debug_assert!(
                        false,
                        "F-point {} lost its C-neighbours; splitting Pass 2 violated",
                        i
                    );
                    row_ptr.push(col_idx.len());
                    continue;
                }

                // d_i = a_ii + ∑_{m ∈ N_i^w} a_im
                // where N_i^w = nonzero row-i entries whose column
                // is neither i, a C_i member, nor an F_i^s member.
                let mut d_i = 0.0f64;
                let start = a.row_ptr[i];
                let end = a.row_ptr[i + 1];
                // Use BTreeSet-free path by leveraging sorted vectors.
                for k in start..end {
                    let col = a.col_idx[k];
                    if col == i {
                        d_i += a.values[k];
                    } else if !c_nbrs.binary_search(&col).is_ok()
                        && !f_nbrs.binary_search(&col).is_ok()
                    {
                        d_i += a.values[k];
                    }
                }

                // Guard against d_i = 0 (degenerate operator row).
                // In practice d_i > 0 for SPD M-matrices (diagonal
                // dominance). Fallback: treat d_i as 1 to avoid
                // divide-by-zero, flag with debug-assert.
                if d_i.abs() < 1e-300 {
                    debug_assert!(
                        false,
                        "effective diagonal d_i vanished at row {}; operator non-SPD?",
                        i
                    );
                }
                let inv_d = if d_i.abs() > 1e-300 { 1.0 / d_i } else { 0.0 };

                // For each C-neighbour k, compute w_ik.
                for &k in &c_nbrs {
                    let a_ik = lookup(i, k);
                    let mut w = a_ik;
                    // F-F common-C correction.
                    for &j in &f_nbrs {
                        let a_ij = lookup(i, j);
                        if a_ij == 0.0 {
                            continue;
                        }
                        let a_jk = lookup(j, k);
                        if a_jk == 0.0 {
                            continue;
                        }
                        // Σ_j = ∑_{l ∈ C_i ∩ S_j} a_jl, accumulated
                        // in the ascending order of l ∈ C_i ∩ S_j.
                        // Traverse c_nbrs in ascending order and
                        // check membership in strong[j] via binary
                        // search.
                        let s_j = &strong[j];
                        let mut sum_j = 0.0f64;
                        for &l in &c_nbrs {
                            if s_j.binary_search(&l).is_ok() {
                                sum_j += lookup(j, l);
                            }
                        }
                        if sum_j.abs() > 1e-300 {
                            w += a_ij * a_jk / sum_j;
                        }
                    }
                    let p_val = -w * inv_d;
                    col_idx.push(coarse_of[k].expect("C-neighbour must have coarse index"));
                    values.push(p_val);
                }
            }
            CfType::Undecided => unreachable!(),
        }
        row_ptr.push(col_idx.len());
    }

    CsrMatrix { n_rows: n, n_cols: n_coarse, row_ptr, col_idx, values }
}

#[cfg(test)]
mod tests {
    use super::super::splitting::{CfType, classical_rs_splitting};
    use super::super::strong_connections::{
        build_poisson_laplacian_csr, compute_strong_connections,
    };
    use super::*;

    #[test]
    fn prolongation_on_poisson_preserves_constants() {
        // Crucial invariant: interpolating the coarse-grid
        // constant vector `1` must reproduce the fine-grid
        // constant vector `1` everywhere (exact polynomial
        // reproduction for constants, per Ruge-Stüben design).
        let n = 12;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let p = build_prolongation(&a, &strong, &cf);

        let n_coarse = p.n_cols;
        let ones_coarse = vec![1.0f64; n_coarse];
        let mut fine = vec![0.0f64; p.n_rows];
        p.apply(&ones_coarse, &mut fine);

        for (i, v) in fine.iter().enumerate() {
            assert!(
                (*v - 1.0).abs() < 1e-12,
                "P·1_coarse fine entry {} = {:.6e}, expected 1",
                i,
                v
            );
        }
    }

    #[test]
    fn prolongation_is_identity_on_coarse_rows() {
        // For every C-point i, row i of P has exactly one entry:
        // value 1.0 at column c(i).
        let n = 10;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let p = build_prolongation(&a, &strong, &cf);

        for i in 0..p.n_rows {
            if cf[i] != CfType::Coarse {
                continue;
            }
            let start = p.row_ptr[i];
            let end = p.row_ptr[i + 1];
            assert_eq!(end - start, 1, "C-row {} has {} entries", i, end - start);
            assert_eq!(p.values[start], 1.0);
        }
    }

    #[test]
    fn prolongation_dimensions_match_coarse_count() {
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let n_coarse = cf.iter().filter(|&&c| c == CfType::Coarse).count();
        let p = build_prolongation(&a, &strong, &cf);
        assert_eq!(p.n_rows, a.n_rows);
        assert_eq!(p.n_cols, n_coarse);
    }

    #[test]
    fn prolongation_is_column_sorted() {
        let n = 6;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let p = build_prolongation(&a, &strong, &cf);
        for i in 0..p.n_rows {
            let start = p.row_ptr[i];
            let end = p.row_ptr[i + 1];
            if end - start <= 1 {
                continue;
            }
            for k in start + 1..end {
                assert!(
                    p.col_idx[k - 1] < p.col_idx[k],
                    "row {} column {:?} not sorted",
                    i,
                    &p.col_idx[start..end]
                );
            }
        }
    }

    #[test]
    fn prolongation_is_deterministic_across_runs() {
        let n = 10;
        let a = build_poisson_laplacian_csr(n);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let first = build_prolongation(&a, &strong, &cf);
        for _ in 0..10 {
            let next = build_prolongation(&a, &strong, &cf);
            assert_eq!(first.n_rows, next.n_rows);
            assert_eq!(first.n_cols, next.n_cols);
            assert_eq!(first.row_ptr, next.row_ptr);
            assert_eq!(first.col_idx, next.col_idx);
            for k in 0..first.values.len() {
                assert_eq!(first.values[k].to_bits(), next.values[k].to_bits());
            }
        }
    }
}
