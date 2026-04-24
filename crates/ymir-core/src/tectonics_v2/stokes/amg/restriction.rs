//! Restriction matrix construction — Step 8.5a Phase 2.3.
//!
//! For Galerkin coarsening, the restriction operator is the
//! transpose of the prolongation: `R = Pᵀ`. This preserves the
//! variational property `A_coarse = R · A_fine · P`, which keeps
//! the coarse operator SPD if the fine-grid one is.
//!
//! # Phase 2.3 status — stub

use super::super::sparse_assembly::CsrMatrix;

/// Build the restriction operator as the CSR transpose of `p`.
///
/// Phase 2.3 stub.
pub fn transpose_to_restriction(_p: &CsrMatrix) -> CsrMatrix {
    panic!("transpose_to_restriction — lands in Phase 2.3");
}
