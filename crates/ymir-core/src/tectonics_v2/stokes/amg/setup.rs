//! AMG hierarchy setup — Step 8.5a Phase 2.6.
//!
//! Orchestrates the sub-phase modules (strong_connections,
//! splitting, prolongation, restriction, coarse_solve) to build
//! a multi-level hierarchy from a fine-grid CSR `A_0`:
//!
//! ```text
//!   loop over k in 0..max_levels:
//!     if a[k].n_rows ≤ min_coarse_unknowns: break
//!     strong[k] = compute_strong_connections(a[k], θ)
//!     cf[k]     = classical_rs_splitting(strong[k])
//!     if no C-points created: break  (hierarchy stuck)
//!     p[k]      = build_prolongation(a[k], strong[k], cf[k])
//!     r[k]      = transpose_to_restriction(p[k])
//!     a[k+1]    = r[k] · a[k] · p[k]     (Galerkin product)
//!   lu         = LuFactorisation::factor(a[last])
//! ```
//!
//! # Determinism (D9)
//!
//! Every operation in the setup chain is deterministic: sparse
//! operations preserve column order via sort-and-merge, BTreeMap
//! accumulators enforce ascending iteration in the Galerkin
//! product, and the sub-phase modules carry their own D9 tests.
//!
//! # Reviewer's vigilance point — coarsest-ratio monitoring
//!
//! Instrumentation (coarse_solve wallclock / V-cycle wallclock
//! fraction) lands in Phase 2.7 benchmark integration. If the
//! ratio exceeds 30 % on any benchmark case, the hierarchy is
//! stopping too early — revise `max_levels` or
//! `min_coarse_unknowns`.

use std::collections::BTreeMap;

use super::super::sparse_assembly::CsrMatrix;
use super::coarse_solve::LuFactorisation;
use super::prolongation::build_prolongation;
use super::restriction::transpose_to_restriction;
use super::splitting::{classical_rs_splitting, CfType};
use super::strong_connections::compute_strong_connections;
use super::{AmgConfig, AmgHierarchy, AmgLevel};

/// Extract the diagonal `N×N` block of a `(2N)×(2N)` CSR matrix:
/// rows and columns in `[offset, offset + n)` — used by Option
/// B' to get the `u-u` and `v-v` scalar blocks from `A_picard`.
pub fn extract_diagonal_block(a: &CsrMatrix, offset: usize, n: usize) -> CsrMatrix {
    assert!(offset + n <= a.n_rows);
    assert!(offset + n <= a.n_cols);
    let mut row_ptr = Vec::with_capacity(n + 1);
    let mut col_idx = Vec::new();
    let mut values = Vec::new();
    row_ptr.push(0);
    for i in 0..n {
        let src_row = offset + i;
        let start = a.row_ptr[src_row];
        let end = a.row_ptr[src_row + 1];
        for k in start..end {
            let src_col = a.col_idx[k];
            // Keep only entries falling in the target block.
            if src_col >= offset && src_col < offset + n {
                col_idx.push(src_col - offset);
                values.push(a.values[k]);
            }
        }
        row_ptr.push(col_idx.len());
    }
    CsrMatrix { n_rows: n, n_cols: n, row_ptr, col_idx, values }
}

/// Galerkin coarsening: compute `R · A · P` as a new CSR.
///
/// Uses a `BTreeMap<usize, f64>` per output row so entries emerge
/// in ascending column order (D9, matches the Phase 1 canonical
/// CSR invariant).
pub fn galerkin_coarsen(r: &CsrMatrix, a: &CsrMatrix, p: &CsrMatrix) -> CsrMatrix {
    assert_eq!(r.n_cols, a.n_rows, "R columns must match A rows");
    assert_eq!(a.n_cols, p.n_rows, "A columns must match P rows");
    let n_coarse_rows = r.n_rows;
    let n_coarse_cols = p.n_cols;
    let mut row_ptr = Vec::with_capacity(n_coarse_rows + 1);
    let mut col_idx = Vec::new();
    let mut values = Vec::new();
    row_ptr.push(0);
    for i in 0..n_coarse_rows {
        let mut acc: BTreeMap<usize, f64> = BTreeMap::new();
        let r_start = r.row_ptr[i];
        let r_end = r.row_ptr[i + 1];
        for rk in r_start..r_end {
            let j = r.col_idx[rk];
            let r_val = r.values[rk];
            let a_start = a.row_ptr[j];
            let a_end = a.row_ptr[j + 1];
            for ak in a_start..a_end {
                let k = a.col_idx[ak];
                let a_val = a.values[ak];
                let ra = r_val * a_val;
                let p_start = p.row_ptr[k];
                let p_end = p.row_ptr[k + 1];
                for pk in p_start..p_end {
                    let l = p.col_idx[pk];
                    let p_val = p.values[pk];
                    *acc.entry(l).or_insert(0.0) += ra * p_val;
                }
            }
        }
        // BTreeMap iteration is ascending → canonical CSR row.
        for (l, v) in acc.iter() {
            col_idx.push(*l);
            values.push(*v);
        }
        row_ptr.push(col_idx.len());
    }
    CsrMatrix { n_rows: n_coarse_rows, n_cols: n_coarse_cols, row_ptr, col_idx, values }
}

/// Build the full multi-level hierarchy from a single scalar SPD
/// block `a_0`.
pub fn build_hierarchy(a_0: CsrMatrix, cfg: AmgConfig) -> AmgHierarchy {
    let mut levels: Vec<AmgLevel> = Vec::new();
    let mut current = a_0;

    for _level_idx in 0..cfg.max_levels {
        // Stop if the current grid is already coarse enough —
        // this level becomes the coarsest, LU-factored.
        if current.n_rows <= cfg.min_coarse_unknowns {
            let lu = LuFactorisation::factor(&current);
            levels.push(AmgLevel {
                a: current,
                p: None,
                r: None,
                coarse_lu: Some(lu),
            });
            return AmgHierarchy { levels };
        }

        // Coarsen one step.
        let strong = compute_strong_connections(&current, cfg.strong_connection_threshold);
        let cf = classical_rs_splitting(&strong);
        let n_coarse = cf.iter().filter(|&&c| c == CfType::Coarse).count();

        // Guard against degenerate coarsening (zero C-points or
        // "no reduction"). If the splitting produced a coarse
        // grid no smaller than fine, stop and factorise current.
        if n_coarse == 0 || n_coarse >= current.n_rows {
            let lu = LuFactorisation::factor(&current);
            levels.push(AmgLevel {
                a: current,
                p: None,
                r: None,
                coarse_lu: Some(lu),
            });
            return AmgHierarchy { levels };
        }

        let p = build_prolongation(&current, &strong, &cf);
        let r = transpose_to_restriction(&p);
        let a_next = galerkin_coarsen(&r, &current, &p);

        levels.push(AmgLevel {
            a: current,
            p: Some(p),
            r: Some(r),
            coarse_lu: None,
        });
        current = a_next;
    }

    // Hit max_levels cap — factorise the current level as coarsest.
    let lu = LuFactorisation::factor(&current);
    levels.push(AmgLevel {
        a: current,
        p: None,
        r: None,
        coarse_lu: Some(lu),
    });
    AmgHierarchy { levels }
}

#[cfg(test)]
mod tests {
    use super::super::strong_connections::build_poisson_laplacian_csr;
    use super::*;

    #[test]
    fn hierarchy_on_poisson_has_multiple_levels() {
        let a = build_poisson_laplacian_csr(16); // 256 unknowns
        let h = build_hierarchy(a, AmgConfig::default());
        // 256 → 128 (≈) → 64 → ≤ 50 → coarsest (LU-factored).
        assert!(h.levels.len() >= 2, "hierarchy has {} levels, expected ≥ 2", h.levels.len());
        // Coarsest level must carry an LU factorisation and no P/R.
        let last = h.levels.last().unwrap();
        assert!(last.coarse_lu.is_some());
        assert!(last.p.is_none());
        assert!(last.r.is_none());
        // All non-coarsest levels must carry P/R.
        for (i, lvl) in h.levels[..h.levels.len() - 1].iter().enumerate() {
            assert!(lvl.p.is_some(), "level {} missing P", i);
            assert!(lvl.r.is_some(), "level {} missing R", i);
            assert!(lvl.coarse_lu.is_none(), "level {} has spurious LU", i);
        }
    }

    #[test]
    fn coarse_matrix_size_decreases_monotonically() {
        let a = build_poisson_laplacian_csr(16);
        let h = build_hierarchy(a, AmgConfig::default());
        for w in h.levels.windows(2) {
            assert!(
                w[1].a.n_rows < w[0].a.n_rows,
                "level size increased: {} -> {}",
                w[0].a.n_rows,
                w[1].a.n_rows
            );
        }
    }

    #[test]
    fn galerkin_preserves_variational_property_on_constants() {
        // R · A · P applied to P·1_coarse should land at the same
        // inner product as A on P·1_coarse. Weaker check: if 1 is
        // in the null-space of the fine operator, 1 is also in
        // the null-space of the coarse operator (up to f64).
        let a = build_poisson_laplacian_csr(8);
        let h = build_hierarchy(a, AmgConfig::default());
        // The fine Poisson with periodic BC has the constant in
        // its null-space. The coarse operator inherits this.
        let a_coarse = &h.levels[1].a;
        let ones = vec![1.0; a_coarse.n_rows];
        let mut y = vec![0.0; a_coarse.n_rows];
        a_coarse.apply(&ones, &mut y);
        let max_y = y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(max_y < 1e-10, "coarse A·1 = max {:.3e}, expected ≈ 0", max_y);
    }

    #[test]
    fn extract_block_preserves_values_and_is_zero_off_block() {
        // Fake 4×4 with cross-coupling: extract the top-left 2×2.
        let a = CsrMatrix {
            n_rows: 4,
            n_cols: 4,
            row_ptr: vec![0, 4, 8, 12, 16],
            col_idx: vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
            values: vec![
                1.0, 2.0, 7.0, 8.0,
                3.0, 4.0, 9.0, 10.0,
                11.0, 12.0, 5.0, 6.0,
                13.0, 14.0, 7.0, 8.0,
            ],
        };
        let top_left = extract_diagonal_block(&a, 0, 2);
        assert_eq!(top_left.n_rows, 2);
        assert_eq!(top_left.n_cols, 2);
        assert_eq!(top_left.values, vec![1.0, 2.0, 3.0, 4.0]);

        let bot_right = extract_diagonal_block(&a, 2, 2);
        assert_eq!(bot_right.values, vec![5.0, 6.0, 7.0, 8.0]);
    }

    #[test]
    fn build_hierarchy_is_deterministic() {
        let a = build_poisson_laplacian_csr(10);
        let h1 = build_hierarchy(a.clone(), AmgConfig::default());
        let h2 = build_hierarchy(a, AmgConfig::default());
        assert_eq!(h1.levels.len(), h2.levels.len());
        for k in 0..h1.levels.len() {
            assert_eq!(h1.levels[k].a.row_ptr, h2.levels[k].a.row_ptr);
            assert_eq!(h1.levels[k].a.col_idx, h2.levels[k].a.col_idx);
            for i in 0..h1.levels[k].a.values.len() {
                assert_eq!(
                    h1.levels[k].a.values[i].to_bits(),
                    h2.levels[k].a.values[i].to_bits()
                );
            }
        }
    }
}
