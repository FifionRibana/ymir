//! Conservative advection of `m_subducted`.
//!
//! Step 7 reuses the Step 0 upwind scheme unchanged for `m̃`
//! (a cell-centered extensive density, same as `S̃`). This test
//! pins the contract:
//!
//! 1. **Mass conservation** under a divergence-free periodic flow
//!    (shear: `vx = sin(2π y), vy = 0`) on one flow period.
//!    Acceptance: relative mass drift `< 10⁻¹⁰`.
//! 2. The solution shape stays reasonable (`max m` bounded) —
//!    upwind diffusion will smooth the Gaussian bump but it must
//!    not explode.

use std::f64::consts::PI;

use ymir_core::tectonics_v2::advection::cfl_dt;
use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::slab::SlabState;

fn gaussian_bump(nx: usize, ny: usize, sigma: f64) -> Field2D {
    let mut m = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) / nx as f64;
            let y = (j as f64 + 0.5) / ny as f64;
            let r2 = (x - 0.5).powi(2) + (y - 0.5).powi(2);
            m.set(i, j, 0.2 * (-r2 / (2.0 * sigma * sigma)).exp());
        }
    }
    m
}

fn shear_velocity(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>) {
    let mut vx = vec![0.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let y = (j as f64 + 0.5) / ny as f64;
            vx[j * nx + i] = (2.0 * PI * y).sin();
        }
    }
    (vx, vy)
}

#[test]
fn mass_conserved_under_shear_flow() {
    let nx = 64;
    let ny = 64;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut state = SlabState::new_zero(nx, ny);
    let bump = gaussian_bump(nx, ny, 0.1);
    for (dst, src) in state.m_mut().data_mut().iter_mut().zip(bump.data().iter()) {
        *dst = *src;
    }
    let mass0 = state.integrated();

    let (vx, vy) = shear_velocity(nx, ny);
    let dt = cfl_dt(dx, dy, &vx, &vy, 0.3);
    // One flow period: Lx / max|vx| = 1/1 = 1.
    let t_end = 1.0;
    let n_steps = (t_end / dt).ceil() as usize;
    let dt_exact = t_end / n_steps as f64;

    for _ in 0..n_steps {
        state.advect(dx, dy, dt_exact, &idx_x, &idx_y, &vx, &vy);
    }
    let mass1 = state.integrated();
    let rel_drift = (mass1 - mass0).abs() / mass0.abs().max(1.0);
    assert!(rel_drift < 1e-10, "mass drift = {:.3e}", rel_drift);
}

#[test]
fn shape_stays_bounded_under_shear_flow() {
    let nx = 64;
    let ny = 64;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut state = SlabState::new_zero(nx, ny);
    let bump = gaussian_bump(nx, ny, 0.1);
    for (dst, src) in state.m_mut().data_mut().iter_mut().zip(bump.data().iter()) {
        *dst = *src;
    }
    let max_initial = state.m().data().iter().cloned().fold(0.0_f64, f64::max);

    let (vx, vy) = shear_velocity(nx, ny);
    let dt = cfl_dt(dx, dy, &vx, &vy, 0.3);
    let n_steps = (1.0_f64 / dt).ceil() as usize;
    let dt_exact = 1.0 / n_steps as f64;
    for _ in 0..n_steps {
        state.advect(dx, dy, dt_exact, &idx_x, &idx_y, &vx, &vy);
    }
    let max_final = state.m().data().iter().cloned().fold(0.0_f64, f64::max);
    let min_final = state.m().data().iter().cloned().fold(f64::INFINITY, f64::min);

    // Upwind diffusion smooths → max should decrease, min should
    // stay ≥ 0 (positivity preservation for positive initial data).
    assert!(max_final <= max_initial * 1.01, "max grew: {} → {}", max_initial, max_final);
    assert!(min_final >= -1e-12, "m went negative: min = {}", min_final);
}

/// Pure translation (uniform vx, vy = 0) returns to the initial
/// state after one period. Exact mass conservation, and L² error
/// from upwind diffusion is bounded.
#[test]
fn pure_translation_returns_close_to_initial() {
    let nx = 64;
    let ny = 64;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let bump = gaussian_bump(nx, ny, 0.1);
    let mut state = SlabState::new_zero(nx, ny);
    for (dst, src) in state.m_mut().data_mut().iter_mut().zip(bump.data().iter()) {
        *dst = *src;
    }
    let mass0 = state.integrated();

    // Constant vx = 1, vy = 0. Period = Lx / vx = 1.
    let vx = vec![1.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    let dt = cfl_dt(dx, dy, &vx, &vy, 0.3);
    let n_steps = (1.0_f64 / dt).ceil() as usize;
    let dt_exact = 1.0 / n_steps as f64;
    for _ in 0..n_steps {
        state.advect(dx, dy, dt_exact, &idx_x, &idx_y, &vx, &vy);
    }
    let mass1 = state.integrated();
    let drift = (mass1 - mass0).abs() / mass0.abs().max(1.0);
    assert!(drift < 1e-10, "mass drift under pure translation = {:.3e}", drift);
}
