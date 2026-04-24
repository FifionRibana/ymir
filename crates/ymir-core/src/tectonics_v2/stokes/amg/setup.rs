//! AMG hierarchy setup — Step 8.5a Phase 2.6 integration.
//!
//! Orchestrates the sub-phase modules (strong_connections,
//! splitting, prolongation, restriction, coarse_solve) to build
//! a multi-level hierarchy from an initial fine-grid CSR `A_0`:
//!
//! ```text
//!   for k in 0..max_levels:
//!     if a[k].n_rows <= min_coarse_unknowns: break
//!     strong[k]   = compute_strong_connections(a[k], θ)
//!     cf[k]       = classical_rs_splitting(strong[k])
//!     p[k]        = build_prolongation(a[k], strong[k], cf[k])
//!     r[k]        = transpose_to_restriction(p[k])
//!     a[k+1]      = r[k] · a[k] · p[k]     # Galerkin product
//!   levels.last().coarse_lu = LuFactorisation::factor(a.last())
//! ```
//!
//! Coarsest-grid ratio monitoring (reviewer's vigilance point 3):
//! if `coarse_solve_time / vcycle_time_total > 30 %` on the
//! benchmark suite, the hierarchy is stopping too early — revise
//! `max_levels` or `min_coarse_unknowns`. Instrumentation added
//! in Phase 2.7.
//!
//! # Phase 2.6 status — stub

use super::super::sparse_assembly::CsrMatrix;
use super::{AmgConfig, AmgHierarchy};

/// Build the multi-level hierarchy for a single SPD scalar block.
///
/// Phase 2.6 stub.
pub fn build_hierarchy(_a: CsrMatrix, _cfg: AmgConfig) -> AmgHierarchy {
    panic!("build_hierarchy — lands in Phase 2.6 (setup orchestration)");
}
