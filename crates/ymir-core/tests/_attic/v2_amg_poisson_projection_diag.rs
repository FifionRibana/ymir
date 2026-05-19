//! Phase 2.7 diagnostic — reviewer-requested investigation of the
//! `poisson_constant` gate miss (AMG = 5 iters vs gate ≤ 3).
//!
//! Three measurements to separate "structural AMG convergence rate"
//! from "null-space projection contributes noise":
//!
//! 1. `|mean(rhs)|` after construction — expected ~ ε_mach · ‖rhs‖.
//! 2. Idempotency: `|mean(subtract_mean(x))|` on an arbitrary x —
//!    must be ≤ ε_mach · ‖x‖. Drift beyond that = bug elsewhere.
//! 3. AMG iter count with projection vs without — if the "without"
//!    variant is drastically fewer, the projection is introducing
//!    systematic error; if similar, the 5-iter figure is the
//!    structural V-cycle convergence rate on a single-mode RHS.
//!
//! Outputs go to stderr via `eprintln!` — `cargo test -- --nocapture`
//! to see them.

use ymir_core::tectonics_v2::stokes::amg::setup::build_hierarchy;
use ymir_core::tectonics_v2::stokes::amg::vcycle::v_cycle;
use ymir_core::tectonics_v2::stokes::amg::AmgConfig;
use ymir_core::tectonics_v2::stokes::nullspace;
use ymir_core::tectonics_v2::stokes::solver::{ConjugateGradient, LinearSolver};
use ymir_core::tectonics_v2::stokes::sparse_assembly::CsrMatrix;

/// 5-pt periodic Laplacian on nx² grid — scalar-coefficient 1.
fn poisson_csr(n: usize) -> CsrMatrix {
    let total = n * n;
    let inv_dx2 = (n as f64) * (n as f64); // dx = 1/n → inv_dx² = n²
    let inv_dy2 = inv_dx2;
    let mut row_ptr = vec![0usize];
    let mut col_idx = Vec::new();
    let mut values = Vec::new();
    for j in 0..n {
        let jp = (j + 1) % n;
        let jm = (j + n - 1) % n;
        for i in 0..n {
            let ip = (i + 1) % n;
            let im = (i + n - 1) % n;
            let diag = 2.0 * inv_dx2 + 2.0 * inv_dy2;
            let mut buf = [
                (jm * n + i, -inv_dy2),
                (j * n + im, -inv_dx2),
                (j * n + i, diag),
                (j * n + ip, -inv_dx2),
                (jp * n + i, -inv_dy2),
            ];
            buf.sort_by_key(|&(c, _)| c);
            let mut prev: Option<usize> = None;
            for (c, v) in buf.iter() {
                match prev {
                    Some(p) if p == *c => *values.last_mut().unwrap() += v,
                    _ => {
                        col_idx.push(*c);
                        values.push(*v);
                        prev = Some(*c);
                    }
                }
            }
            row_ptr.push(col_idx.len());
        }
    }
    CsrMatrix { n_rows: total, n_cols: total, row_ptr, col_idx, values }
}

fn sin_sin_rhs(n: usize) -> Vec<f64> {
    let mut rhs = vec![0.0f64; n * n];
    for j in 0..n {
        for i in 0..n {
            let x = (i as f64 + 0.5) / n as f64;
            let y = (j as f64 + 0.5) / n as f64;
            rhs[j * n + i] = (2.0 * std::f64::consts::PI * x).sin()
                * (2.0 * std::f64::consts::PI * y).sin();
        }
    }
    // Match the bench: subtract mean once at construction.
    nullspace::subtract_mean(&mut rhs);
    rhs
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[test]
fn reviewer_requested_projection_diagnostic() {
    let n = 64;
    let total = n * n;
    let a = poisson_csr(n);
    let rhs = sin_sin_rhs(n);

    // --- Measurement 1: |mean(rhs)| after construction ---
    let mean_rhs: f64 = rhs.iter().sum::<f64>() / total as f64;
    let norm_rhs = norm2(&rhs);
    let eps_mach_bound = f64::EPSILON * norm_rhs;
    eprintln!(
        "[diag 1] |mean(rhs)| = {:.3e}  ‖rhs‖₂ = {:.3e}  ε_mach·‖rhs‖ = {:.3e}",
        mean_rhs.abs(),
        norm_rhs,
        eps_mach_bound
    );
    // sin(2πx)·sin(2πy) has exact zero mean on a uniform grid of
    // full periods — residual should be ≤ ε_mach · ‖rhs‖.
    assert!(
        mean_rhs.abs() <= 1e-14 * norm_rhs,
        "rhs mean too large: {:.3e}",
        mean_rhs
    );

    // --- Measurement 2: subtract_mean idempotency on a seeded
    // arbitrary vector ---
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    let mut rng = ChaCha8Rng::seed_from_u64(42);
    let mut x: Vec<f64> = (0..total).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
    let norm_x0 = norm2(&x);
    nullspace::subtract_mean(&mut x);
    let mean_after_once = x.iter().sum::<f64>() / total as f64;
    nullspace::subtract_mean(&mut x);
    let mean_after_twice = x.iter().sum::<f64>() / total as f64;
    eprintln!(
        "[diag 2] subtract_mean idempotency: mean after 1st = {:.3e}, after 2nd = {:.3e}  \
         (ε_mach·‖x‖ = {:.3e})",
        mean_after_once.abs(),
        mean_after_twice.abs(),
        f64::EPSILON * norm_x0
    );
    // The 2nd call should leave the mean at ε_mach · ‖x‖ relative
    // order. If it doesn't, there is a bug in `subtract_mean`.
    assert!(
        mean_after_twice.abs() <= 1e-14 * norm_x0,
        "subtract_mean non-idempotent: 2nd-call mean = {:.3e} (expected ε_mach-bounded)",
        mean_after_twice,
    );

    // --- Measurement 3: AMG iter count WITH vs WITHOUT projection ---
    let cfg = AmgConfig::default();
    let h = build_hierarchy(a.clone(), cfg);

    // Variant WITH projection (current bench precond).
    let iters_with = {
        let cg = ConjugateGradient::new(1e-8, 4000);
        let mut x = vec![0.0f64; total];
        let a_ref = a.clone();
        let h_ref = &h;
        let cfg_ref = &cfg;
        let mut matvec = |v: &[f64], out: &mut [f64]| a_ref.apply(v, out);
        let mut precond = |r: &[f64], z: &mut [f64]| {
            for v in z.iter_mut() {
                *v = 0.0;
            }
            let mut r_proj = r.to_vec();
            nullspace::subtract_mean(&mut r_proj);
            v_cycle(h_ref, cfg_ref, &r_proj, z);
            nullspace::subtract_mean(z);
        };
        let stats = cg.solve(&mut matvec, &mut precond, &rhs, &mut x);
        stats.iterations
    };

    // Variant WITHOUT projection — apply V-cycle directly, no
    // subtract_mean. If the hierarchy preserves the null-space
    // naturally on a zero-mean RHS, no projection is needed.
    let iters_without = {
        let cg = ConjugateGradient::new(1e-8, 4000);
        let mut x = vec![0.0f64; total];
        let a_ref = a.clone();
        let h_ref = &h;
        let cfg_ref = &cfg;
        let mut matvec = |v: &[f64], out: &mut [f64]| a_ref.apply(v, out);
        let mut precond = |r: &[f64], z: &mut [f64]| {
            for v in z.iter_mut() {
                *v = 0.0;
            }
            v_cycle(h_ref, cfg_ref, r, z);
        };
        let stats = cg.solve(&mut matvec, &mut precond, &rhs, &mut x);
        stats.iterations
    };

    eprintln!(
        "[diag 3] poisson_constant CG iter count: WITH projection = {}, \
         WITHOUT projection = {}",
        iters_with, iters_without,
    );
    eprintln!(
        "[diag 3] Delta = {}; positive → projection slows CG, negative → projection accelerates",
        iters_with as i64 - iters_without as i64,
    );

    // --- Measurement 3.5: norm of what projection removed per iter ---
    // Measure on the first residual: what mean does `r_0 = rhs - A·0 = rhs`
    // carry into the precond?
    let mean_r0 = rhs.iter().sum::<f64>() / total as f64;
    eprintln!(
        "[diag 3.5] mean(r_0) entering precond = {:.3e}  \
         ε_mach-normalised = {:.3e}",
        mean_r0.abs(),
        mean_r0.abs() / f64::EPSILON.max(1e-30),
    );

    // --- Measurement 3.6: V-cycle output mean on zero-mean input ---
    // Apply V-cycle to rhs (which is zero-mean); observe the output
    // mean before projection. If the V-cycle drifts the mean by
    // much more than ε_mach · ‖rhs‖, the projection has "real work"
    // to do — otherwise it's just noise removal.
    let mut z_nopproj = vec![0.0f64; total];
    v_cycle(&h, &cfg, &rhs, &mut z_nopproj);
    let mean_z = z_nopproj.iter().sum::<f64>() / total as f64;
    let norm_z = norm2(&z_nopproj);
    eprintln!(
        "[diag 3.6] V-cycle(rhs) output: ‖z‖₂ = {:.3e}, mean(z) = {:.3e}, \
         mean/‖z‖ = {:.3e}, ε_mach = {:.3e}",
        norm_z,
        mean_z.abs(),
        mean_z.abs() / norm_z.max(1e-30),
        f64::EPSILON,
    );
}
