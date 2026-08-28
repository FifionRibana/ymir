//! C/F splitting via the Classical Ruge-Stüben two-pass algorithm
//! — Step 8.5a Phase 2.2.
//!
//! Partitions the fine-grid DOFs into **C-points** (kept on the
//! coarser level) and **F-points** (interpolated from C-points).
//! Reference: Briggs-Henson-McCormick, "A Multigrid Tutorial",
//! 2nd ed., §8.8.2 — the original Ruge-Stüben 1987 two-pass
//! algorithm (not CLJP / PMIS variants which target parallel
//! execution; 8.5a is single-threaded).
//!
//! # Pass 1 — Standard colouring (importance-driven)
//!
//! - λ_i = |S_iᵀ| = number of points that i strongly influences
//!   (= i's "importance" as a candidate C-point).
//! - While any point is `Undecided`:
//!   1. Pick `i*` with maximum λ. Ties broken by **lowest
//!      index** (reviewer's vigilance point 2). Mark as Coarse.
//!   2. For each undecided `j ∈ S_{i*}ᵀ` (points i* influences):
//!      mark as `Fine`; increment λ for each undecided `k ∈ S_j`
//!      (j's influencers become more important candidates).
//!   3. For each undecided `j ∈ S_{i*}` (points influencing i*):
//!      decrement λ (they've been "served" by this new C).
//!
//! # Pass 2 — F-F strong-connection enforcement
//!
//! For each F-point i, every strong F-F pair (i, j) must share a
//! common strong C-neighbour for the Ruge-Stüben prolongation to
//! have an interpolation formula. When no such C-neighbour exists,
//! promote i (the point encountered first in ascending-index
//! order) to C. This preserves the "every F-point has a full
//! interpolation stencil" invariant that prolongation needs.
//!
//! # Determinism (D9 — reviewer's top-priority invariant)
//!
//! - λ ordering: `BTreeSet<(-λ, i)>` — `iter().next()` returns
//!   max λ, lowest i at every step. No HashMap, no randomised
//!   tiebreak.
//! - Pass 2 promotion order: ascending `i`, break on first missing
//!   common C-neighbour.
//! - Verified by `deterministic_100_runs` test: 100 independent
//!   invocations with the same `strong` input produce the exact
//!   same `Vec<CfType>`.
//!
//! # Complexity
//!
//! Pass 1: O(|E| log N) via BTreeSet operations, where |E| = total
//! strong-connection count. Pass 2: O(|E| · max_row_width) from
//! the common-neighbour check. Setup is one-shot per Newton outer
//! iter (η frozen), so the cost is amortised per V-cycle.

use std::collections::BTreeSet;

/// Label for a fine-grid DOF after splitting.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CfType {
    Coarse,
    Fine,
    /// Used only during the splitting passes; no output entry
    /// carries `Undecided`.
    Undecided,
}

/// Run Classical Ruge-Stüben two-pass C/F splitting on the
/// supplied strong-connection structure (per-row sorted column
/// lists, as emitted by [`super::strong_connections::
/// compute_strong_connections`]).
///
/// Returns one [`CfType`] per row. All entries are `Coarse` or
/// `Fine`; `Undecided` never leaks into the output.
pub fn classical_rs_splitting(strong: &[Vec<usize>]) -> Vec<CfType> {
    let n = strong.len();
    let mut cf = vec![CfType::Undecided; n];

    // --- Transpose the strong-connection structure so that
    //     strong_t[j] = { i : j ∈ strong[i] } = "points that j
    //     strongly influences". Kept ascending-sorted for D9.
    let mut strong_t: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for &j in &strong[i] {
            strong_t[j].push(i);
        }
    }
    for v in strong_t.iter_mut() {
        // Each j receives at most one push per i (strong[i] has
        // no duplicates), so strong_t[j] already has no dupes.
        // Sort to make the traversal order deterministic.
        v.sort();
    }

    // --- λ_i initialisation and active-set BTreeSet ---
    // Active set orders by (-λ, i) so the first entry is
    // always max-λ with lowest i.
    let mut lambda: Vec<i64> = strong_t.iter().map(|v| v.len() as i64).collect();
    let mut active: BTreeSet<(i64, usize)> = BTreeSet::new();
    for i in 0..n {
        active.insert((-lambda[i], i));
    }

    // --- Pass 1 ---
    while let Some(&entry) = active.iter().next() {
        active.remove(&entry);
        let (neg_li, i_star) = entry;
        debug_assert_eq!(-neg_li, lambda[i_star]);
        debug_assert_eq!(cf[i_star], CfType::Undecided);

        cf[i_star] = CfType::Coarse;

        // (a) Points that i_star influences → mark F, bump their
        // own strong-influencers' λ by 1.
        let influenced: Vec<usize> = strong_t[i_star].clone();
        for &j in &influenced {
            if cf[j] != CfType::Undecided {
                continue;
            }
            // Remove j from active, mark Fine.
            active.remove(&(-lambda[j], j));
            cf[j] = CfType::Fine;

            // Bump λ for undecided k that strongly influence j.
            for &k in &strong[j] {
                if cf[k] != CfType::Undecided {
                    continue;
                }
                active.remove(&(-lambda[k], k));
                lambda[k] += 1;
                active.insert((-lambda[k], k));
            }
        }

        // (b) Points that strongly influence i_star → decrement
        // their λ (already-selected C-neighbour reduces their
        // "importance" as future C candidates).
        let influencers: Vec<usize> = strong[i_star].clone();
        for &j in &influencers {
            if cf[j] != CfType::Undecided {
                continue;
            }
            active.remove(&(-lambda[j], j));
            lambda[j] -= 1;
            active.insert((-lambda[j], j));
        }
    }

    // Sanity: Pass 1 must decide every point.
    debug_assert!(cf.iter().all(|c| *c != CfType::Undecided));

    // --- Pass 2: F-F strong connections need a common C-neighbour ---
    //
    // For each F-point i, for each F-point j in strong[i], verify
    // that some k in strong[i] ∩ strong[j] is Coarse. If no such
    // k exists, promote i to C (and break inner loop; i's later
    // F-F checks are moot since i is no longer F).
    for i in 0..n {
        if cf[i] != CfType::Fine {
            continue;
        }
        // Collect i's C-neighbours from strong[i] (as a BTreeSet
        // so the intersection check is O(log N) per j).
        let c_nbrs_i: BTreeSet<usize> =
            strong[i].iter().filter(|&&k| cf[k] == CfType::Coarse).copied().collect();

        for &j in &strong[i] {
            if cf[j] != CfType::Fine {
                continue;
            }
            let has_common = strong[j].iter().any(|k| c_nbrs_i.contains(k));
            if !has_common {
                // Promote i to C to break the orphan F-F link.
                cf[i] = CfType::Coarse;
                break;
            }
        }
    }

    cf
}

#[cfg(test)]
mod tests {
    use super::super::strong_connections::{
        build_poisson_laplacian_csr, compute_strong_connections,
    };
    use super::*;

    fn count_coarse_fine(cf: &[CfType]) -> (usize, usize) {
        let c = cf.iter().filter(|x| **x == CfType::Coarse).count();
        let f = cf.iter().filter(|x| **x == CfType::Fine).count();
        (c, f)
    }

    #[test]
    fn no_undecided_in_output_on_poisson() {
        let a = build_poisson_laplacian_csr(8);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        assert!(cf.iter().all(|c| *c != CfType::Undecided));
        let (c_count, f_count) = count_coarse_fine(&cf);
        assert_eq!(c_count + f_count, cf.len());
        assert!(c_count > 0, "splitting produced zero C-points");
        assert!(f_count > 0, "splitting produced zero F-points");
    }

    #[test]
    fn every_f_has_at_least_one_c_neighbour() {
        // After Pass 2 enforcement, every F-point must have a C
        // neighbour among its strong set (else prolongation would
        // have no interpolation sources).
        let a = build_poisson_laplacian_csr(12);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        for i in 0..cf.len() {
            if cf[i] != CfType::Fine {
                continue;
            }
            let has_c = strong[i].iter().any(|&j| cf[j] == CfType::Coarse);
            assert!(has_c, "F-point {} has no C-neighbour in its strong set {:?}", i, strong[i],);
        }
    }

    #[test]
    fn coarsening_ratio_is_reasonable_on_poisson() {
        // Classical RS on a regular Poisson tile produces a
        // coarsening ratio of roughly 1/2 (close to red-black-like
        // coarsening). We allow wide margins to keep the test
        // robust to small tie-breaking differences. The hard
        // invariant: both sets are non-empty.
        let a = build_poisson_laplacian_csr(16);
        let strong = compute_strong_connections(&a, 0.25);
        let cf = classical_rs_splitting(&strong);
        let (c, f) = count_coarse_fine(&cf);
        let ratio = c as f64 / (c + f) as f64;
        assert!(
            (0.25..=0.75).contains(&ratio),
            "coarsening ratio {:.3} out of [0.25, 0.75]",
            ratio
        );
    }

    /// Reviewer's top-priority invariant: D9 determinism.
    ///
    /// 100 invocations with identical input must produce byte-
    /// identical output. If `BTreeSet` ordering or λ tie-breaking
    /// ever leaks non-determinism, this test catches it.
    #[test]
    fn deterministic_100_runs_on_poisson() {
        let a = build_poisson_laplacian_csr(10);
        let strong = compute_strong_connections(&a, 0.25);
        let first = classical_rs_splitting(&strong);
        for run in 0..100 {
            let next = classical_rs_splitting(&strong);
            assert_eq!(first, next, "splitting non-deterministic at run {} — D9 violated", run);
        }
    }

    #[test]
    fn lowest_index_wins_on_lambda_tie() {
        // Construct a tiny fully-symmetric strong structure where
        // every node influences every other → all λ equal, 3 nodes.
        // The lowest-index (0) MUST be picked first as the initial
        // C-point.
        let strong = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let cf = classical_rs_splitting(&strong);
        assert_eq!(cf[0], CfType::Coarse, "lowest-index tiebreak violated");
    }

    #[test]
    fn isolated_points_are_coarse_by_default() {
        // No strong connections → λ_i = 0 for all i. The algorithm
        // picks points one by one with λ=0; each becomes C (no
        // influenced points to mark F). So all points should end
        // up C.
        let strong = vec![Vec::new(); 4];
        let cf = classical_rs_splitting(&strong);
        for (i, t) in cf.iter().enumerate() {
            assert_eq!(*t, CfType::Coarse, "isolated point {} marked {:?}", i, t);
        }
    }
}
