//! Strong-connection detection — Step 8.5a Phase 2.1.
//!
//! For each row `i` of an SPD matrix `A`, a column `j ≠ i` is
//! **strongly connected** to `i` if
//! ```text
//!     -a_ij ≥ θ · max_{k ≠ i} (-a_ik)
//! ```
//! (Briggs-Henson-McCormick Ch. 8.8, eq. 8.43). This is the
//! Ruge-Stüben classical definition for SPD systems where the
//! off-diagonal entries are non-positive. Strong connections
//! drive both the C/F splitting (Phase 2.2) and the prolongation
//! weights (Phase 2.3).
//!
//! Per-row output is a `Vec<usize>` of column indices that are
//! strong for that row. Order is ascending (by the CSR column
//! layout invariant from Phase 1) — D9 determinism.
//!
//! # Phase 2.1 status — stub
//!
//! Real implementation lands in the next commit.

use super::super::sparse_assembly::CsrMatrix;

/// Compute strong-connection sets per row of `a` using threshold
/// `theta`. Returns a nested vector where index `i` is the sorted
/// list of column indices strongly connected to row `i`.
///
/// Phase 2.1 stub.
pub fn compute_strong_connections(_a: &CsrMatrix, _theta: f64) -> Vec<Vec<usize>> {
    panic!("compute_strong_connections — lands in Phase 2.1");
}
