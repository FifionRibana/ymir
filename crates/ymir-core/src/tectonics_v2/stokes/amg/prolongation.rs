//! Prolongation matrix construction — Step 8.5a Phase 2.3.
//!
//! Builds the interpolation matrix `P` that maps coarse-grid
//! corrections to fine-grid corrections:
//! ```text
//!     v_fine = P · v_coarse
//! ```
//! For a C-point `i`, `P[i, c(i)] = 1` where `c(i)` is its index
//! in the coarse-grid ordering. For an F-point, the Ruge-Stüben
//! classical formula assigns weights to each strongly-connected
//! C-point neighbour proportional to the off-diagonal magnitudes.
//! See Briggs-Henson-McCormick eq. (8.45)-(8.48).
//!
//! # Phase 2.3 status — stub

use super::super::sparse_assembly::CsrMatrix;
use super::splitting::CfType;

/// Build the prolongation CSR from the fine-grid operator, strong
/// connections, and the C/F labelling.
///
/// Phase 2.3 stub.
pub fn build_prolongation(
    _a: &CsrMatrix,
    _strong: &[Vec<usize>],
    _cf: &[CfType],
) -> CsrMatrix {
    panic!("build_prolongation — lands in Phase 2.3 (Ruge-Stüben classical formula)");
}
