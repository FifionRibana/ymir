//! Method-of-manufactured-solutions convergence test for the Step-0
//! constant-η thin viscous sheet solver on the MAC grid.
//!
//! Manufactured solution — deliberately **NOT** divergence-free, so
//! both the normal-stress (`-η ∇²v_i`) and the grad-div
//! (`-η ∂_i(∇·v)`) parts of the operator are exercised:
//! ```text
//!   v_exact = ( sin(2π x),  sin(2π y) )
//!   ∇·v     = 2π ( cos(2π x) + cos(2π y) )
//! ```
//! With η = 1, `-∇·(2 η ε̇(v)) = -η ( ∇²v + ∇(∇·v) )`, component-wise:
//! ```text
//!   f_x = -η ∇² v_x - η ∂_x(∇·v) = 4π² sin(2π x) + 4π² sin(2π x) = 8π² sin(2π x)
//!   f_y = -η ∇² v_y - η ∂_y(∇·v) = 4π² sin(2π y) + 4π² sin(2π y) = 8π² sin(2π y)
//! ```
//! A staggered-MAC scheme with harmonic η-averaging is expected to
//! be 2nd order; a slope departing from 2 by more than 0.3 flags a
//! discretization bug — in particular, dropping the grad-div term
//! would show up here as a visible order loss.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::{Grid, SheetConfig, solve_sheet};

fn build_mms(nx: usize, ny: usize, dx: f64, dy: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = nx * ny;
    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut vx_exact = vec![0.0; n];
    let mut vy_exact = vec![0.0; n];

    for j in 0..ny {
        for i in 0..nx {
            // vx, fx at (i·dx, (j+0.5)·dy).
            let xfx = i as f64 * dx;
            fx[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * xfx).sin();
            vx_exact[j * nx + i] = (2.0 * PI * xfx).sin();

            // vy, fy at ((i+0.5)·dx, j·dy).
            let yfy = j as f64 * dy;
            fy[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * yfy).sin();
            vy_exact[j * nx + i] = (2.0 * PI * yfy).sin();
        }
    }
    (fx, fy, vx_exact, vy_exact)
}

fn rms(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    (a.iter().map(|v| v * v).sum::<f64>() / n).sqrt()
}

fn solve_and_error(n: usize) -> f64 {
    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);

    let (fx, fy, vx_exact, vy_exact) = build_mms(nx, ny, dx, dy);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut cfg = SheetConfig::default();
    cfg.tol = 1e-12;
    cfg.max_iter = 5000;
    let stats = solve_sheet(&grid, &eta, None, &fx, &fy, &mut vx, &mut vy, &cfg);
    assert!(stats.converged, "CG did not converge at N={}: {:?}", n, stats);

    let v_err: Vec<f64> = vx
        .iter()
        .zip(vx_exact.iter())
        .chain(vy.iter().zip(vy_exact.iter()))
        .map(|(a, b)| a - b)
        .collect();
    rms(&v_err)
}

#[test]
fn velocity_converges_at_second_order() {
    let sizes = [16usize, 32, 64, 128];
    let errs: Vec<f64> = sizes.iter().map(|&n| solve_and_error(n)).collect();
    for (n, e) in sizes.iter().zip(errs.iter()) {
        eprintln!("MMS N={:4} | v_err={:.3e}", n, e);
    }
    for k in 0..errs.len() - 1 {
        let slope = (errs[k] / errs[k + 1]).log2();
        eprintln!("  slope v({}→{}) = {:.3}", sizes[k], sizes[k + 1], slope);
    }
    let slope_final = (errs[errs.len() - 2] / errs[errs.len() - 1]).log2();
    assert!(
        slope_final >= 1.7 && slope_final <= 2.3,
        "final v slope = {:.3} (expected ≈ 2.0 ± 0.3)",
        slope_final,
    );
}
