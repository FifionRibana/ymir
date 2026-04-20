//! Null-space verification for the periodic Stokes solve.
//!
//! The solve is fed a body force with deliberately nonzero mean in
//! each component. The expected behaviour: the projector kills those
//! null-space components and the iterates are returned with machine
//! precision on `|mean(P)|`, `|mean(vx)|`, `|mean(vy)|`.

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::{solve_stokes, Grid, StokesConfig};

#[test]
fn solve_returns_zero_mean_pressure_and_velocity() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);

    // RHS with deliberately nonzero mean in every component.
    let mut fx = vec![0.0; nx * ny];
    let mut fy = vec![0.0; nx * ny];
    use std::f64::consts::PI;
    for j in 0..ny {
        for i in 0..nx {
            let x = i as f64 * dx;
            let y = (j as f64 + 0.5) * dy;
            fx[j * nx + i] = 0.3 + (2.0 * PI * y).sin(); // mean = 0.3
            let y2 = j as f64 * dy;
            fy[j * nx + i] = -0.5 + (2.0 * PI * x).cos(); // mean = -0.5 + 0 = -0.5
        }
    }

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut p = vec![0.0; nx * ny];
    let cfg = StokesConfig::default();
    let stats = solve_stokes(&grid, &eta, &fx, &fy, &mut vx, &mut vy, &mut p, &cfg);

    assert!(stats.converged, "outer CG did not converge: {:?}", stats);

    let mean_abs = |v: &[f64]| {
        let m: f64 = v.iter().sum::<f64>() / v.len() as f64;
        m.abs()
    };
    let mp = mean_abs(&p);
    let mvx = mean_abs(&vx);
    let mvy = mean_abs(&vy);
    assert!(mp < 1e-10, "|mean(P)| = {:.3e}", mp);
    assert!(mvx < 1e-10, "|mean(vx)| = {:.3e}", mvx);
    assert!(mvy < 1e-10, "|mean(vy)| = {:.3e}", mvy);
}

#[test]
fn zero_forcing_yields_zero_solution() {
    let nx = 16;
    let ny = 16;
    let grid = Grid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
    let eta = Field2D::filled(nx, ny, 1.0);
    let fx = vec![0.0; nx * ny];
    let fy = vec![0.0; nx * ny];
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut p = vec![0.0; nx * ny];
    let cfg = StokesConfig::default();
    solve_stokes(&grid, &eta, &fx, &fy, &mut vx, &mut vy, &mut p, &cfg);
    let peak = vx.iter().chain(vy.iter()).chain(p.iter()).fold(0.0_f64, |a, &v| a.max(v.abs()));
    assert!(peak < 1e-12, "nonzero solution from zero RHS: peak = {}", peak);
}
