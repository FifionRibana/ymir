//! Null-space verification for the periodic thin-sheet solve.
//!
//! The solve is fed a body force with deliberately nonzero mean in
//! each component. The expected behaviour: the projector kills those
//! null-space components and the iterates are returned with machine
//! precision on `|mean(vx)|` and `|mean(vy)|`. The thin-sheet
//! formulation has **no pressure unknown**, hence no pressure-null-
//! space check.

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::{solve_sheet, Grid, SheetConfig};

#[test]
fn solve_returns_zero_mean_velocity() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);

    let mut fx = vec![0.0; nx * ny];
    let mut fy = vec![0.0; nx * ny];
    use std::f64::consts::PI;
    for j in 0..ny {
        for i in 0..nx {
            let x = i as f64 * dx;
            let y = (j as f64 + 0.5) * dy;
            fx[j * nx + i] = 0.3 + (2.0 * PI * y).sin();
            let x2 = (i as f64 + 0.5) * dx;
            fy[j * nx + i] = -0.5 + (2.0 * PI * x2).cos();
            // Silence unused-var lint on x when only x2 is consumed.
            let _ = x;
        }
    }

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let cfg = SheetConfig::default();
    let stats = solve_sheet(&grid, &eta, None, &fx, &fy, &mut vx, &mut vy, &cfg);

    assert!(stats.converged, "CG did not converge: {:?}", stats);

    let mean_abs = |v: &[f64]| {
        let m: f64 = v.iter().sum::<f64>() / v.len() as f64;
        m.abs()
    };
    let mvx = mean_abs(&vx);
    let mvy = mean_abs(&vy);
    assert!(mvx < 1e-10, "|mean(vx)| = {:.3e}", mvx);
    assert!(mvy < 1e-10, "|mean(vy)| = {:.3e}", mvy);
}

#[test]
fn zero_forcing_yields_zero_velocity() {
    let nx = 16;
    let ny = 16;
    let grid = Grid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
    let eta = Field2D::filled(nx, ny, 1.0);
    let fx = vec![0.0; nx * ny];
    let fy = vec![0.0; nx * ny];
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let cfg = SheetConfig::default();
    solve_sheet(&grid, &eta, None, &fx, &fy, &mut vx, &mut vy, &cfg);
    let peak = vx.iter().chain(vy.iter()).fold(0.0_f64, |a, &v| a.max(v.abs()));
    assert!(peak < 1e-12, "nonzero velocity from zero RHS: peak = {}", peak);
}
