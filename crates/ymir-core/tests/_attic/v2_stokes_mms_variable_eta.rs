//! MMS convergence test with prescribed spatially varying η.
//!
//! Manufactured solution (divergence-free-enough for the test — the
//! thin-sheet operator does not require it):
//! ```text
//!   v_exact = ( sin(2πx)·cos(2πy), -cos(2πx)·sin(2πy) )
//!   η(x,y)  = 1 + 0.5·sin(2πx)·cos(2πy)
//! ```
//! For this particular v_exact the shear strain rate ε̇_xy is
//! identically zero (since ∂_y vx = -∂_x vy); hence σ_xy ≡ 0 and the
//! momentum residual simplifies to `f_x = -∂_x h`, `f_y = +∂_y h`
//! with `h(x,y) = 4π · η(x,y) · cos(2πx)·cos(2πy)`. Expanding:
//! ```text
//!   h = 4π·cos(2πx)cos(2πy) + π·sin(4πx)·cos²(2πy)
//!   ∂_x h = -8π²·sin(2πx)cos(2πy) + 4π²·cos(4πx)·cos²(2πy)
//!   ∂_y h = -8π²·cos(2πx)sin(2πy) - 2π²·sin(4πx)·sin(4πy)
//! ```
//! so
//! ```text
//!   f_x =  8π²·sin(2πx)cos(2πy) - 4π²·cos(4πx)·cos²(2πy)
//!   f_y = -8π²·cos(2πx)sin(2πy) - 2π²·sin(4πx)·sin(4πy)
//! ```
//! This test only exercises the **linear** η-variable path of
//! `apply_momentum` (via `solve_sheet`); the nonlinear Newton path is
//! covered by `v2_newton_convergence`.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::{Grid, SheetConfig, solve_sheet};

fn build_mms(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
) -> (Field2D, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut eta = Field2D::new(nx, ny);
    let n = nx * ny;
    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut vx_exact = vec![0.0; n];
    let mut vy_exact = vec![0.0; n];

    for j in 0..ny {
        for i in 0..nx {
            // η at cell centre.
            let xc = (i as f64 + 0.5) * dx;
            let yc = (j as f64 + 0.5) * dy;
            let e = 1.0 + 0.5 * (2.0 * PI * xc).sin() * (2.0 * PI * yc).cos();
            eta.set(i, j, e);

            // vx, fx at (i dx, (j+0.5) dy).
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dy;
            vx_exact[j * nx + i] = (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos();
            fx[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos()
                - 4.0 * PI * PI * (4.0 * PI * xf).cos() * (2.0 * PI * yf).cos().powi(2);

            // vy, fy at ((i+0.5) dx, j dy).
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dy;
            vy_exact[j * nx + i] = -(2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin();
            fy[j * nx + i] = -8.0 * PI * PI * (2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin()
                - 2.0 * PI * PI * (4.0 * PI * xf2).sin() * (4.0 * PI * yf2).sin();
        }
    }
    (eta, fx, fy, vx_exact, vy_exact)
}

fn rms(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    (a.iter().map(|v| v * v).sum::<f64>() / n).sqrt()
}

fn error_for_grid(n: usize) -> f64 {
    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let (eta, fx, fy, vx_ex, vy_ex) = build_mms(nx, ny, dx, dy);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut cfg = SheetConfig::default();
    cfg.tol = 1.0e-12;
    cfg.max_iter = 10_000;
    let stats = solve_sheet(&grid, &eta, None, &fx, &fy, &mut vx, &mut vy, &cfg);
    assert!(stats.converged, "CG did not converge at N={}: {:?}", n, stats);

    let err: Vec<f64> = vx
        .iter()
        .zip(vx_ex.iter())
        .chain(vy.iter().zip(vy_ex.iter()))
        .map(|(a, b)| a - b)
        .collect();
    rms(&err)
}

#[test]
fn variable_eta_second_order_convergence() {
    let sizes = [32usize, 64, 128];
    let errs: Vec<f64> = sizes.iter().map(|&n| error_for_grid(n)).collect();
    for (n, e) in sizes.iter().zip(errs.iter()) {
        eprintln!("variable-η MMS N={:4} | v_err={:.3e}", n, e);
    }
    for k in 0..errs.len() - 1 {
        let slope = (errs[k] / errs[k + 1]).log2();
        eprintln!("  slope v({}→{}) = {:.3}", sizes[k], sizes[k + 1], slope);
    }
    let slope_final = (errs[errs.len() - 2] / errs[errs.len() - 1]).log2();
    assert!(slope_final >= 1.7, "final slope = {:.3} (expected ≥ 1.7)", slope_final,);
}
