//! Direct LU solve on the coarsest level — Step 8.5a Phase 2.5.
//!
//! Reimplements Doolittle LU with partial pivoting (~80 LOC) to
//! avoid pulling `nalgebra` into the workspace solely for this.
//! The coarsest matrix is always ≤ `min_coarse_unknowns` (default
//! 50), so O(n³) factorisation is ~125 k ops and the solve is
//! negligible per V-cycle.
//!
//! Determinism: Doolittle's natural implementation is O(n³)
//! deterministic by construction. Partial pivoting chooses the
//! largest-magnitude entry in the active column; ties by lowest
//! row index (per D9).
//!
//! # Phase 2.5 status — stub

use super::super::sparse_assembly::CsrMatrix;

/// Factorisation of a small dense SPD matrix via Doolittle LU.
#[derive(Debug)]
pub struct LuFactorisation {
    pub n: usize,
    /// Combined L (strict lower) + U (upper incl diag), row-major.
    pub lu: Vec<f64>,
    /// Row-permutation from partial pivoting.
    pub perm: Vec<usize>,
}

impl LuFactorisation {
    /// Factorise the coarse matrix. Phase 2.5 stub.
    pub fn factor(_a: &CsrMatrix) -> Self {
        panic!("LuFactorisation::factor — lands in Phase 2.5 (Doolittle LU)");
    }

    /// Solve `a · x = b` using the factorisation. Phase 2.5 stub.
    pub fn solve(&self, _b: &[f64], _x: &mut [f64]) {
        panic!("LuFactorisation::solve — lands in Phase 2.5");
    }
}
