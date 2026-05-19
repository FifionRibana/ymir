//! Full Multigrid orchestration — Step 8.5a Phase 4.1.
//!
//! Brandt 1977's "Multi-Level Adaptive Solutions" recipe:
//!
//! 1. Restrict the RHS level-by-level to the coarsest grid.
//! 2. Solve directly on the coarsest grid.
//! 3. Prolongate the coarse solution as the initial guess for the
//!    next-finer level.
//! 4. Apply a V-cycle at that level to refine.
//! 5. Repeat (prolongate + V-cycle) up to the finest grid.
//!
//! Intuition: the V-cycle alone, starting from `x = 0` on the fine
//! grid, spends its first pass representing the low-frequency
//! solution that a cheaper coarse-grid solve already captures.
//! FMG hands the V-cycle a "good enough" initial guess, so the
//! V-cycle only has to refine the high-frequency residual — often
//! in a single pass.
//!
//! Gate (Phase 4.1): FMG achieves ≥ 2× additional iteration
//! reduction vs V-cycle-alone on Poisson; test
//! [`tests::fmg_beats_v_cycle_on_poisson_by_at_least_2x`].

use super::vcycle::v_cycle_level;
use super::{AmgConfig, AmgHierarchy};

/// Apply one FMG cycle: restrict `b` down, solve coarsest, then
/// prolongate + V-cycle on the way up, updating `x` in place.
pub fn fmg_cycle(h: &AmgHierarchy, cfg: &AmgConfig, b: &[f64], x: &mut [f64]) {
    let n_levels = h.levels.len();
    debug_assert!(n_levels >= 1, "hierarchy must have at least one level");

    // Start with a zeroed fine-grid state; FMG's entire point is
    // to feed the V-cycle a non-trivial initial guess.
    for v in x.iter_mut() {
        *v = 0.0;
    }

    if n_levels == 1 {
        // Degenerate case — only the coarse level exists. Solve
        // and return.
        let lu = h.levels[0]
            .coarse_lu
            .as_ref()
            .expect("single-level hierarchy must be factored");
        lu.solve(b, x);
        return;
    }

    // --- Phase 1: restrict b to every level ---
    let mut b_levels: Vec<Vec<f64>> = Vec::with_capacity(n_levels);
    b_levels.push(b.to_vec());
    for k in 0..(n_levels - 1) {
        let r = h.levels[k]
            .r
            .as_ref()
            .expect("non-coarsest level must have R");
        let mut b_next = vec![0.0f64; r.n_rows];
        r.apply(&b_levels[k], &mut b_next);
        b_levels.push(b_next);
    }

    // --- Phase 2: solve the coarsest grid directly ---
    let coarsest = n_levels - 1;
    let lu = h.levels[coarsest]
        .coarse_lu
        .as_ref()
        .expect("coarsest level must have LU");
    let mut x_current = vec![0.0f64; b_levels[coarsest].len()];
    lu.solve(&b_levels[coarsest], &mut x_current);

    // --- Phase 3: prolongate + V-cycle at each finer level ---
    for k in (0..coarsest).rev() {
        let p = h.levels[k]
            .p
            .as_ref()
            .expect("non-coarsest level must have P");
        let mut x_k = vec![0.0f64; p.n_rows];
        p.apply(&x_current, &mut x_k);
        v_cycle_level(h, cfg, k, &b_levels[k], &mut x_k);
        x_current = x_k;
    }

    // Copy the finest-level result into the caller's `x`.
    debug_assert_eq!(x.len(), x_current.len());
    x.copy_from_slice(&x_current);
}

#[cfg(test)]
mod tests {
    use super::super::setup::build_hierarchy;
    use super::super::strong_connections::build_poisson_laplacian_csr;
    use super::super::vcycle::v_cycle;
    use super::*;

    fn norm2(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    fn residual_norm(
        a: &super::super::super::sparse_assembly::CsrMatrix,
        b: &[f64],
        x: &[f64],
    ) -> f64 {
        let mut ax = vec![0.0; b.len()];
        a.apply(x, &mut ax);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        norm2(&r)
    }

    fn seeded_zero_mean_rhs(n: usize, seed: u64) -> Vec<f64> {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let total = n * n;
        let mut b: Vec<f64> = (0..total).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let m: f64 = b.iter().sum::<f64>() / total as f64;
        for v in b.iter_mut() {
            *v -= m;
        }
        b
    }

    #[test]
    fn fmg_converges_on_poisson() {
        // Invariant: one FMG application must achieve at least
        // ≥ 20× residual reduction (= 5% of initial). Measured
        // baseline on the current implementation is ~66× so the
        // margin is safe. This is the "FMG does useful work in
        // one shot" sanity check, not a tight performance gate.
        let n = 16;
        let a = build_poisson_laplacian_csr(n);
        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);
        let b = seeded_zero_mean_rhs(n, 42);
        let r0 = residual_norm(&a, &b, &vec![0.0; n * n]);
        let mut x = vec![0.0f64; n * n];
        fmg_cycle(&h, &cfg, &b, &mut x);
        let r1 = residual_norm(&a, &b, &x);
        let reduction = r1 / r0;
        eprintln!("[fmg] single-shot reduction ratio = {:.3e}", reduction);
        assert!(
            reduction < 0.05,
            "FMG single application gave only {:.3e} reduction ({:.3e} → {:.3e})",
            reduction,
            r0,
            r1
        );
    }

    #[test]
    fn fmg_beats_v_cycle_on_poisson_by_at_least_2x() {
        // Phase 4.1 gate: FMG one-shot residual reduction must be
        // ≥ 2× that of V-cycle one-shot from the same zero initial
        // guess.
        let n = 16;
        let a = build_poisson_laplacian_csr(n);
        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);
        let b = seeded_zero_mean_rhs(n, 7);
        let r0 = residual_norm(&a, &b, &vec![0.0; n * n]);

        let mut x_v = vec![0.0f64; n * n];
        v_cycle(&h, &cfg, &b, &mut x_v);
        let r_v = residual_norm(&a, &b, &x_v);

        let mut x_f = vec![0.0f64; n * n];
        fmg_cycle(&h, &cfg, &b, &mut x_f);
        let r_f = residual_norm(&a, &b, &x_f);

        let v_reduction = r_v / r0;
        let f_reduction = r_f / r0;
        let advantage = v_reduction / f_reduction;
        eprintln!(
            "[fmg vs v-cycle] r_0 = {:.3e}, r_v = {:.3e} (×{:.3e}), \
             r_fmg = {:.3e} (×{:.3e}), FMG advantage = {:.2}×",
            r0, r_v, v_reduction, r_f, f_reduction, advantage,
        );
        assert!(
            advantage >= 2.0,
            "FMG advantage over V-cycle = {:.2}× (expected ≥ 2×)",
            advantage
        );
    }

    #[test]
    fn fmg_leaves_zero_rhs_at_zero() {
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);
        let b = vec![0.0f64; n * n];
        let mut x = vec![0.0f64; n * n];
        fmg_cycle(&h, &cfg, &b, &mut x);
        for v in x.iter() {
            assert!(v.abs() < 1e-13, "FMG drifted from zero: {:.3e}", v);
        }
    }

    #[test]
    fn fmg_is_deterministic() {
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);
        let b = seeded_zero_mean_rhs(n, 999);
        let mut x1 = vec![0.0f64; n * n];
        let mut x2 = vec![0.0f64; n * n];
        fmg_cycle(&h, &cfg, &b, &mut x1);
        fmg_cycle(&h, &cfg, &b, &mut x2);
        for k in 0..(n * n) {
            assert_eq!(x1[k].to_bits(), x2[k].to_bits());
        }
    }
}
