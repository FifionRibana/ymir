//! Analytic tests for first-order upwind S advection.
//!
//! 1. Mass conservation: prescribed divergence-free shear flow, one
//!    flow period. Relative drift < 1e-12 required.
//! 2. Numerical diffusion: report L² error without gating. Logged for
//!    future reference when higher-order schemes come online.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::advection::{cfl_dt, integrated_mass, step_upwind};
use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};

/// Build an initial Gaussian bump centered at (0.5, 0.5) with width σ.
fn gaussian_bump(nx: usize, ny: usize, sigma: f64) -> Field2D {
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let r2 = (x - 0.5).powi(2) + (y - 0.5).powi(2);
            s.set(i, j, 1.0 + 0.5 * (-r2 / (2.0 * sigma * sigma)).exp());
        }
    }
    s
}

/// Divergence-free periodic shear flow: vx = sin(2π y), vy = 0.
/// (∂_x vx + ∂_y vy = 0.)
fn shear_velocity(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>) {
    let mut vx = vec![0.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    for j in 0..ny {
        // vx is at face (i·dx, (j+0.5)·dy). Sample at (j+0.5)/ny.
        for i in 0..nx {
            let y = (j as f64 + 0.5) / ny as f64;
            vx[j * nx + i] = (2.0 * PI * y).sin();
        }
    }
    (vx, vy)
}

#[test]
fn mass_is_conserved_under_shear_flow() {
    let nx = 64;
    let ny = 64;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut s = gaussian_bump(nx, ny, 0.1);
    let mass0 = integrated_mass(&s);
    let (vx, vy) = shear_velocity(nx, ny);

    // The bump is advected in x; a full period of the imposed flow is
    // Lx / max|vx| = 1 / 1 = 1.
    let t_end = 1.0;
    let dt = cfl_dt(dx, dy, &vx, &vy, 0.3);
    let n_steps = (t_end / dt).ceil() as usize;
    let dt_exact = t_end / n_steps as f64;

    let mut s_next = Field2D::new(nx, ny);
    for _ in 0..n_steps {
        step_upwind(nx, ny, dx, dy, dt_exact, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
        std::mem::swap(&mut s, &mut s_next);
    }

    let mass1 = integrated_mass(&s);
    let rel_drift = (mass1 - mass0).abs() / mass0.abs().max(1.0);
    assert!(rel_drift < 1e-12, "mass drift = {:.3e}", rel_drift);
}

#[test]
fn numerical_diffusion_reported_without_gating() {
    // After one period the exact solution in a pure-translation frame
    // returns to the initial condition. The L² error is a measure of
    // numerical diffusion from first-order upwind. Report, don't gate.
    let nx = 64;
    let ny = 64;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    // Use a pure x-translation flow so the reference solution after
    // one period is the initial field itself.
    let vy = vec![0.0; nx * ny];
    let vx = vec![1.0; nx * ny];
    let mut s = gaussian_bump(nx, ny, 0.1);
    let s0 = s.data().to_vec();

    let dt = cfl_dt(dx, dy, &vx, &vy, 0.5);
    let n_steps = (1.0 / dt).ceil() as usize;
    let dt_exact = 1.0 / n_steps as f64;

    let mut s_next = Field2D::new(nx, ny);
    for _ in 0..n_steps {
        step_upwind(nx, ny, dx, dy, dt_exact, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
        std::mem::swap(&mut s, &mut s_next);
    }
    let l2 = s.data().iter().zip(s0.iter()).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt()
        / (nx as f64 * ny as f64).sqrt();
    eprintln!("first-order upwind L² error after one period: {:.3e}", l2);
    // Sanity bound: shouldn't be crazy (<1, given amplitude < 0.5).
    assert!(l2 < 0.5, "suspicious L² error {}", l2);
}
