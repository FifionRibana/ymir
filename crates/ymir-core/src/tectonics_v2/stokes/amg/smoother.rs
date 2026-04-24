//! Symmetric Gauss-Seidel smoother — Step 8.5a Phase 2.4.
//!
//! One "sweep" of Symmetric Gauss-Seidel is one forward traversal
//! (rows 0..N in ascending order) + one backward traversal (rows
//! N-1..0 in descending order), each updating
//! ```text
//!     x_i ← (1 / a_ii) · (b_i − ∑_{j ≠ i} a_ij · x_j)
//! ```
//! using the most-recently-updated components. Symmetric sweeps
//! preserve SPD structure required by CG; one-directional GS
//! alone would not. D8 default: 1 symmetric sweep pre + 1
//! symmetric sweep post (equivalent to 2 unsymmetric passes each
//! side in the AMG literature's alternate convention).
//!
//! Test target (Phase 2.4 gate): on a Poisson constant-coeff
//! problem, a single SGS sweep reduces the residual norm by
//! approximately 0.5 (the classical spectral-radius-based rate
//! estimate).
//!
//! # Phase 2.4 status — stub

use super::super::sparse_assembly::CsrMatrix;

/// Apply one symmetric Gauss-Seidel sweep: `x` is updated in
/// place, `b` is the right-hand side.
///
/// Phase 2.4 stub.
pub fn sgs_sweep(_a: &CsrMatrix, _b: &[f64], _x: &mut [f64]) {
    panic!("sgs_sweep — lands in Phase 2.4 (symmetric Gauss-Seidel)");
}
