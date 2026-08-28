//! Step 4 — method-of-manufactured-solutions convergence test for
//! the basal-drag-augmented momentum operator.
//!
//! # Manufactured solution
//!
//! ```text
//!   v(x, y)   = ( sin(2π x) · cos(2π y),  -cos(2π x) · sin(2π y) )
//!   S̃(x, y)   = 1 + 0.3 · cos(2π x) · cos(2π y)
//!   η(x, y)   = 1        (constant)
//!   Br        = 0.1
//! ```
//!
//! `S̃` is **even in x and even in y** so that `S̃² · v_exact` has
//! zero mean on the periodic unit square (the drag forcing then
//! contains no null-space component, and the gauge-projected RHS
//! consumed by `solve_sheet` agrees with the pointwise RHS built
//! here). A `sin(2πx)·cos(2πy)` variant of `S̃` *appears* to
//! preserve symmetry but its `S̃²` cross-term `0.6·sin·cos·sin·cos`
//! integrates to `0.15`, producing a non-zero RHS mean that the
//! solver removes — CG then converges to a mean-corrected v ≠
//! v_exact and slope tests fail. Keep `S̃` here even-even.
//!
//! `v` is deliberately divergence-free: `∂_x vx + ∂_y vy = 0`. With
//! `η = 1` constant, the viscous operator reduces to `-η∇²v` (the
//! grad-div term vanishes), so component-wise:
//!
//! ```text
//!   ∇² vx  = -8π² · sin(2πx)·cos(2πy) = -8π² · vx
//!   ∇² vy  = -8π² · (-cos(2πx)·sin(2πy)) = -8π² · vy
//!   -∇² v  = 8π² · v
//! ```
//!
//! The basal-drag term contributes `Br · S̃²(x, y) · v(x, y)` (a
//! positive diagonal coupling), so the continuous RHS is
//!
//! ```text
//!   f(x, y) = -∇·(2 η ε̇(v)) + Br · S̃²(x, y) · v(x, y)
//!           = (8π² + Br · S̃²(x, y)) · v(x, y)
//! ```
//!
//! evaluated pointwise at the MAC face locations. The face-averaging
//! convention inside `apply_momentum` introduces an O(dx²) error in
//! resolving the face drag from the cell-centered S̃ — consistent
//! with the 2nd-order discretisation of the stencil. Thus the overall
//! velocity error is O(dx²); the test asserts slope ≥ 1.7 at
//! N ∈ {32, 64, 128}.
//!
//! # Why this is load-bearing
//!
//! If the drag augmentation in `apply_momentum` (face-interpolation
//! convention) or `momentum_diagonal` (analytical diagonal
//! reconstruction, case B) drifts from the correct staggered
//! arithmetic, the assembled operator will not match the analytic
//! RHS — slope drops from ~2 to ~1 or the solve fails to converge.
//! This is the principal MMS garde-fou for the Step 4 operator
//! integration.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::basal_drag::{
    BasalDragConfig, BasalDragLaw, build_drag_diagonal_field,
};
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::{Grid, SheetConfig, solve_sheet};

const BR: f64 = 0.1;

/// `S̃(x, y) = 1 + 0.3 · cos(2π x) · cos(2π y)` — even-even to keep
/// the drag forcing null-space-free on the periodic domain.
fn s_tilde(x: f64, y: f64) -> f64 {
    1.0 + 0.3 * (2.0 * PI * x).cos() * (2.0 * PI * y).cos()
}

/// Build the MMS setup on an N×N grid and return (η, S̃ field,
/// RHS fx, RHS fy, exact vx, exact vy).
fn build_mms(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
) -> (Field2D, Field2D, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut eta = Field2D::filled(nx, ny, 1.0);
    // S̃ at cell centres — this is what `build_drag_diagonal_field`
    // consumes. The `apply_momentum` drag loop averages Br·S̃² from
    // cells to faces, which converges to the face-wise continuous
    // value at O(dx²).
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let xc = (i as f64 + 0.5) * dx;
            let yc = (j as f64 + 0.5) * dy;
            s.set(i, j, s_tilde(xc, yc));
        }
    }
    let _ = &mut eta;

    let n = nx * ny;
    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut vx_exact = vec![0.0; n];
    let mut vy_exact = vec![0.0; n];
    // Viscous contribution factor: -∇²v = 8π²·v for v sinusoidal.
    let visc_factor = 8.0 * PI * PI;

    for j in 0..ny {
        for i in 0..nx {
            // vx at (i·dx, (j+0.5)·dy).
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dy;
            let vx_e = (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos();
            vx_exact[j * nx + i] = vx_e;
            let drag_x = BR * s_tilde(xf, yf).powi(2);
            fx[j * nx + i] = (visc_factor + drag_x) * vx_e;

            // vy at ((i+0.5)·dx, j·dy).
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dy;
            let vy_e = -(2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin();
            vy_exact[j * nx + i] = vy_e;
            let drag_y = BR * s_tilde(xf2, yf2).powi(2);
            fy[j * nx + i] = (visc_factor + drag_y) * vy_e;
        }
    }
    (eta, s, fx, fy, vx_exact, vy_exact)
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
    let (eta, s, fx, fy, vx_ex, vy_ex) = build_mms(nx, ny, dx, dy);
    let drag_cfg = BasalDragConfig::Enabled(BasalDragLaw { br: BR, s_exponent: 2.0 });
    let drag_diag =
        build_drag_diagonal_field(&drag_cfg, &s).expect("Enabled drag config must produce a field");

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut cfg = SheetConfig::default();
    cfg.tol = 1.0e-12;
    cfg.max_iter = 10_000;
    let stats = solve_sheet(&grid, &eta, Some(&drag_diag), &fx, &fy, &mut vx, &mut vy, &cfg);
    assert!(stats.converged, "CG did not converge at N={n}: {stats:?}",);

    let err: Vec<f64> = vx
        .iter()
        .zip(vx_ex.iter())
        .chain(vy.iter().zip(vy_ex.iter()))
        .map(|(a, b)| a - b)
        .collect();
    rms(&err)
}

#[test]
fn basal_drag_mms_converges_at_second_order() {
    let sizes = [32usize, 64, 128];
    let errs: Vec<f64> = sizes.iter().map(|&n| error_for_grid(n)).collect();
    for (n, e) in sizes.iter().zip(errs.iter()) {
        eprintln!("basal-drag MMS N={n:4} | v_err={e:.3e}");
    }
    for k in 0..errs.len() - 1 {
        let slope = (errs[k] / errs[k + 1]).log2();
        eprintln!("  slope v({}→{}) = {:.3}", sizes[k], sizes[k + 1], slope);
    }
    let slope_final = (errs[errs.len() - 2] / errs[errs.len() - 1]).log2();
    // Spec acceptance: slope ≥ 1.7 at the finest pair.
    assert!(
        slope_final >= 1.7,
        "final basal-drag MMS slope = {slope_final:.3} (expected ≥ 1.7). \
         A slope drop flags an inconsistency in the drag augmentation \
         (apply_momentum vs momentum_diagonal, face-averaging, or \
         sign convention). Errors: {errs:?}, sizes: {sizes:?}.",
    );
}
