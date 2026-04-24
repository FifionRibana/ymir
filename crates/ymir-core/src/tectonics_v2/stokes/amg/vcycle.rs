//! V-cycle orchestration — Step 8.5a Phase 2.6.
//!
//! Textbook recursive V-cycle:
//!
//! ```text
//!   fn vcycle(level k):
//!     if k is the coarsest:
//!         x_k ← lu.solve(b_k)
//!         return
//!     pre_smooth(a_k, b_k, x_k, pre_sweeps)
//!     residual_k ← b_k - a_k · x_k
//!     b_{k+1}   ← r_k · residual_k         # restriction
//!     x_{k+1}   ← 0
//!     vcycle(level k+1)
//!     x_k       += p_k · x_{k+1}            # prolongation
//!     post_smooth(a_k, b_k, x_k, post_sweeps)
//! ```
//!
//! Convergence-target: Phase 2.6 gate — V-cycle on Poisson
//! constant η converges in ≤ 3 CG-wrapped iterations, single-
//! V-cycle residual reduction ≥ 10× (classical AMG figure).
//! Heterogeneous tests follow in Phase 2.7.
//!
//! # Phase 2.6 status — stub

use super::AmgHierarchy;

/// Apply one V-cycle to the hierarchy, updating `x` in-place
/// toward the solution of `A_0 · x = b`.
///
/// Phase 2.6 stub.
pub fn v_cycle(_h: &AmgHierarchy, _b: &[f64], _x: &mut [f64]) {
    panic!("v_cycle — lands in Phase 2.6 (V-cycle orchestration)");
}
