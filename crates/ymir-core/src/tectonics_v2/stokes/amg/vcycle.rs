//! V-cycle orchestration — Step 8.5a Phase 2.6.
//!
//! Textbook recursive V-cycle on the hierarchy built by
//! [`super::setup::build_hierarchy`]:
//!
//! ```text
//!   fn vcycle(level k):
//!     if k is the coarsest:
//!         x_k ← lu.solve(b_k)
//!         return
//!     for _ in 0..pre_sweeps: sgs_sweep(a_k, b_k, x_k)
//!     residual_k ← b_k − a_k · x_k
//!     b_{k+1}    ← r_k · residual_k         # restriction
//!     x_{k+1}    ← 0
//!     vcycle(level k+1)
//!     x_k        += p_k · x_{k+1}           # prolongation
//!     for _ in 0..post_sweeps: sgs_sweep(a_k, b_k, x_k)
//! ```
//!
//! Recursion depth is bounded by `cfg.max_levels` (default 7),
//! so Rust's stack handles it without issue.

use super::smoother::rbgs_sweep;
use super::{AmgConfig, AmgHierarchy};

/// Apply one V-cycle starting at level 0: update `x` in-place
/// toward the solution of `levels[0].a · x = b`.
pub fn v_cycle(h: &AmgHierarchy, cfg: &AmgConfig, b: &[f64], x: &mut [f64]) {
    debug_assert!(!h.levels.is_empty(), "hierarchy is empty");
    v_cycle_level(h, cfg, 0, b, x);
}

/// Apply a V-cycle starting from level `k` of the hierarchy. Public
/// for use by FMG (`super::fmg`); external callers should prefer
/// the level-0 entry point [`v_cycle`].
pub fn v_cycle_level(
    h: &AmgHierarchy,
    cfg: &AmgConfig,
    k: usize,
    b: &[f64],
    x: &mut [f64],
) {
    let lvl = &h.levels[k];

    // Coarsest level: direct LU solve.
    if k == h.levels.len() - 1 {
        if let Some(lu) = lvl.coarse_lu.as_ref() {
            lu.solve(b, x);
        } else {
            debug_assert!(false, "coarsest level missing LU factorisation");
        }
        return;
    }

    let p = lvl.p.as_ref().expect("non-coarsest level must have P");
    let r = lvl.r.as_ref().expect("non-coarsest level must have R");

    // --- Pre-smooth ---
    for _ in 0..cfg.pre_smooth_sweeps {
        rbgs_sweep(&lvl.a, &lvl.colors, b, x);
    }

    // --- Residual r = b - A · x ---
    let n_fine = lvl.a.n_rows;
    let mut ax = vec![0.0f64; n_fine];
    lvl.a.apply(x, &mut ax);
    let mut residual = vec![0.0f64; n_fine];
    for i in 0..n_fine {
        residual[i] = b[i] - ax[i];
    }

    // --- Restrict: b_coarse = R · residual ---
    let n_coarse = r.n_rows;
    let mut b_coarse = vec![0.0f64; n_coarse];
    r.apply(&residual, &mut b_coarse);

    // --- Coarse solve (recurse) ---
    let mut x_coarse = vec![0.0f64; n_coarse];
    v_cycle_level(h, cfg, k + 1, &b_coarse, &mut x_coarse);

    // --- Prolongate and correct: x += P · x_coarse ---
    let mut correction = vec![0.0f64; n_fine];
    p.apply(&x_coarse, &mut correction);
    for i in 0..n_fine {
        x[i] += correction[i];
    }

    // --- Post-smooth ---
    for _ in 0..cfg.post_smooth_sweeps {
        rbgs_sweep(&lvl.a, &lvl.colors, b, x);
    }
}

#[cfg(test)]
mod tests {
    use super::super::setup::build_hierarchy;
    use super::super::strong_connections::build_poisson_laplacian_csr;
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

    #[test]
    fn vcycle_on_poisson_converges_fast() {
        // Random-seeded RHS with zero mean (null-space consistent);
        // run V-cycles until the residual drops by ≥ 10⁶ or we hit
        // a cycle budget. Classical AMG figure on Poisson: per-
        // V-cycle convergence factor ≈ 0.1 (= ρ_vcycle), so 6
        // V-cycles should suffice for 10⁶ reduction.
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let n = 16;
        let a = build_poisson_laplacian_csr(n);
        let total = n * n;
        let mut rng = ChaCha8Rng::seed_from_u64(7);
        let mut b: Vec<f64> = (0..total).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let mean: f64 = b.iter().sum::<f64>() / total as f64;
        for v in b.iter_mut() {
            *v -= mean;
        }

        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);

        let mut x = vec![0.0; total];
        let r0 = residual_norm(&a, &b, &x);
        let mut last = r0;
        let mut ratios: Vec<f64> = Vec::new();
        for k in 0..10 {
            v_cycle(&h, &cfg, &b, &mut x);
            let r_new = residual_norm(&a, &b, &x);
            ratios.push(r_new / last);
            if r_new < r0 * 1e-6 {
                eprintln!(
                    "v_cycle converged in {} cycles; per-cycle ratios {:?}",
                    k + 1,
                    ratios
                );
                return;
            }
            last = r_new;
        }
        panic!(
            "V-cycle failed to reduce residual by 10⁶ in 10 cycles; final {:.3e}/{:.3e}, ratios {:?}",
            last, r0, ratios
        );
    }

    #[test]
    fn vcycle_leaves_zero_rhs_at_zero() {
        // If b = 0 and x = 0, one V-cycle must leave x at 0 (no
        // spurious drift from the smoother or coarse-grid solve).
        let a = build_poisson_laplacian_csr(8);
        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);
        let n = 64;
        let b = vec![0.0f64; n];
        let mut x = vec![0.0f64; n];
        v_cycle(&h, &cfg, &b, &mut x);
        for v in x.iter() {
            assert!(v.abs() < 1e-13, "V-cycle drift from zero: {:.3e}", v);
        }
    }

    #[test]
    fn vcycle_is_deterministic() {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let cfg = AmgConfig::default();
        let h = build_hierarchy(a.clone(), cfg);
        let total = n * n;
        let mut rng = ChaCha8Rng::seed_from_u64(999);
        let b: Vec<f64> = (0..total).map(|_| rng.random::<f64>()).collect();
        let mut x1 = vec![0.0; total];
        let mut x2 = vec![0.0; total];
        for _ in 0..3 {
            v_cycle(&h, &cfg, &b, &mut x1);
            v_cycle(&h, &cfg, &b, &mut x2);
        }
        for k in 0..total {
            assert_eq!(x1[k].to_bits(), x2[k].to_bits());
        }
    }
}
