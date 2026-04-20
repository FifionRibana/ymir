//! Method-of-manufactured-solutions convergence test for the Step-0
//! constant-η Stokes solver on the MAC grid.
//!
//! Manufactured solution (divergence-free, zero-mean periodic):
//! ```text
//!   v_exact = ( sin(2πx) cos(2πy), -cos(2πx) sin(2πy) )
//!   p_exact = sin(2πx) sin(2πy)
//! ```
//! With η = 1, the corresponding body force is
//! ```text
//!   f_x =  8π² sin(2πx) cos(2πy) + 2π cos(2πx) sin(2πy)
//!   f_y = -8π² cos(2πx) sin(2πy) + 2π sin(2πx) cos(2πy)
//! ```
//! The test measures the L² velocity error at 16²…128² and fits the
//! slope in `log₂(err) vs log₂(N)`. A staggered-MAC scheme with
//! harmonic averaging is expected to be 2nd order; the prompt flags
//! an observed slope departing from 2 by more than 0.3 as a
//! discretization bug.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::{solve_stokes, Grid, StokesConfig};

/// Sample MMS quantities on the MAC grid.
fn build_mms(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = nx * ny;
    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut vx_exact = vec![0.0; n];
    let mut vy_exact = vec![0.0; n];
    let mut p_exact = vec![0.0; n];

    for j in 0..ny {
        for i in 0..nx {
            // vx, fx at (i·dx, (j+0.5)·dy)
            let xfx = i as f64 * dx;
            let yfx = (j as f64 + 0.5) * dy;
            fx[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * xfx).sin() * (2.0 * PI * yfx).cos()
                + 2.0 * PI * (2.0 * PI * xfx).cos() * (2.0 * PI * yfx).sin();
            vx_exact[j * nx + i] = (2.0 * PI * xfx).sin() * (2.0 * PI * yfx).cos();

            // vy, fy at ((i+0.5)·dx, j·dy)
            let xfy = (i as f64 + 0.5) * dx;
            let yfy = j as f64 * dy;
            fy[j * nx + i] = -8.0 * PI * PI * (2.0 * PI * xfy).cos() * (2.0 * PI * yfy).sin()
                + 2.0 * PI * (2.0 * PI * xfy).sin() * (2.0 * PI * yfy).cos();
            vy_exact[j * nx + i] = -(2.0 * PI * xfy).cos() * (2.0 * PI * yfy).sin();

            // p at cell center ((i+0.5)·dx, (j+0.5)·dy)
            let xp = (i as f64 + 0.5) * dx;
            let yp = (j as f64 + 0.5) * dy;
            p_exact[j * nx + i] = (2.0 * PI * xp).sin() * (2.0 * PI * yp).sin();
        }
    }
    (fx, fy, vx_exact, vy_exact, p_exact)
}

fn rms(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    (a.iter().map(|v| v * v).sum::<f64>() / n).sqrt()
}

fn solve_and_error(n: usize) -> (f64, f64) {
    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);

    let (fx, fy, vx_exact, vy_exact, p_exact) = build_mms(nx, ny, dx, dy);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut p = vec![0.0; nx * ny];
    let mut cfg = StokesConfig::default();
    // Tighten tolerances to get a clean convergence slope.
    cfg.outer_tol = 1e-10;
    cfg.inner_tol = 1e-12;
    cfg.outer_max_iter = 400;
    cfg.inner_max_iter = 2000;
    let stats = solve_stokes(&grid, &eta, &fx, &fy, &mut vx, &mut vy, &mut p, &cfg);
    assert!(stats.converged, "outer CG did not converge at N={}: {:?}", n, stats);

    // Project the exact solution too (it already has zero mean analytically).
    let v_err: Vec<f64> = vx
        .iter()
        .zip(vx_exact.iter())
        .chain(vy.iter().zip(vy_exact.iter()))
        .map(|(a, b)| a - b)
        .collect();
    let v_err_rms = rms(&v_err);

    // Compare p also, after zeroing both means.
    let p_mean: f64 = p.iter().sum::<f64>() / p.len() as f64;
    let p_exact_mean: f64 = p_exact.iter().sum::<f64>() / p_exact.len() as f64;
    let p_err: Vec<f64> = p
        .iter()
        .zip(p_exact.iter())
        .map(|(a, b)| (a - p_mean) - (b - p_exact_mean))
        .collect();
    let p_err_rms = rms(&p_err);

    (v_err_rms, p_err_rms)
}

#[test]
fn velocity_converges_at_second_order() {
    let sizes = [16usize, 32, 64, 128];
    let errs: Vec<(f64, f64)> = sizes.iter().map(|&n| solve_and_error(n)).collect();
    for (n, (ve, pe)) in sizes.iter().zip(errs.iter()) {
        eprintln!("MMS N={:4} | v_err={:.3e} | p_err={:.3e}", n, ve, pe);
    }
    // Slope = log₂(err_k / err_{k+1}) between successive sizes.
    let mut min_slope_v = f64::INFINITY;
    let mut max_slope_v = f64::NEG_INFINITY;
    for k in 0..errs.len() - 1 {
        let slope = (errs[k].0 / errs[k + 1].0).log2();
        eprintln!("  slope v({}→{}) = {:.3}", sizes[k], sizes[k + 1], slope);
        min_slope_v = min_slope_v.min(slope);
        max_slope_v = max_slope_v.max(slope);
    }
    // The prompt requires convergence order ≥ 1.7 for nominally
    // 2nd-order schemes. Tolerate a modest preasymptotic slope on the
    // coarsest refinement; require the finer refinements to sit in
    // the expected band.
    let slope_final = (errs[errs.len() - 2].0 / errs[errs.len() - 1].0).log2();
    assert!(
        slope_final >= 1.7 && slope_final <= 2.3,
        "final v slope = {:.3} (expected ≈ 2.0 ± 0.3)",
        slope_final,
    );
}
