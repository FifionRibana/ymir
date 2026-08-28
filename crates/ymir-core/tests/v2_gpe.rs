//! Integration tests for `GpeForce` (the Step 2 body force).
//!
//! Unit-level properties (linearity in Ar, additive accumulation,
//! uniform-S → zero, integral-zero on periodic domain, smooth-field
//! analytic match) live in `src/tectonics_v2/forcing/gpe.rs`. This
//! file adds the whole-crate integration checks: smooth-field order
//! of accuracy under refinement, and sanity on a mildly perturbed
//! `S̃` that the baseline uses.

use std::f64::consts::TAU;

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::forcing::{BodyForce, GpeForce, SimulationState, VectorField};

fn build_state<'a>(
    nx: usize,
    ny: usize,
    dx: f64,
    idx_x: &'a PeriodicIndex,
    idx_y: &'a PeriodicIndex,
    s: &'a Field2D,
) -> SimulationState<'a> {
    SimulationState { nx, ny, dx, dy: dx, idx_x, idx_y, s }
}

#[test]
fn smooth_field_force_converges_at_second_order() {
    // Known closed-form gradient of ½·S² for S = 1 + α sin(kx) cos(ky):
    //   ∂_x(½ S²) = S · α k cos(kx) cos(ky)
    //   ∂_y(½ S²) = -S · α k sin(kx) sin(ky)
    // So f_x = -Ar·S·α k cos(kx)cos(ky) at (i dx, (j+½) dy).
    let ar = 2.0;
    let alpha = 0.1;
    let k = TAU;

    let sizes = [32usize, 64, 128];
    let errors: Vec<f64> = sizes
        .iter()
        .map(|&n| {
            let dx = 1.0 / n as f64;
            let idx_x = PeriodicIndex::new(n);
            let idx_y = PeriodicIndex::new(n);
            let mut s = Field2D::new(n, n);
            for j in 0..n {
                for i in 0..n {
                    let x = (i as f64 + 0.5) * dx;
                    let y = (j as f64 + 0.5) * dx;
                    s.set(i, j, 1.0 + alpha * (k * x).sin() * (k * y).cos());
                }
            }
            let mut fx = Field2D::new(n, n);
            let mut fy = Field2D::new(n, n);
            GpeForce::with_ar(ar).accumulate(
                &build_state(n, n, dx, &idx_x, &idx_y, &s),
                &mut VectorField { fx: &mut fx, fy: &mut fy },
            );
            // Analytic: evaluate at each face and compute L² error.
            let mut sq_err = 0.0_f64;
            let mut count = 0_usize;
            for j in 0..n {
                for i in 0..n {
                    let xf = i as f64 * dx;
                    let yf = (j as f64 + 0.5) * dx;
                    let s_at = 1.0 + alpha * (k * xf).sin() * (k * yf).cos();
                    let dsdx = alpha * k * (k * xf).cos() * (k * yf).cos();
                    let analytic = -ar * s_at * dsdx;
                    let numeric = fx.data()[j * n + i];
                    sq_err += (numeric - analytic).powi(2);
                    count += 1;
                }
            }
            (sq_err / count as f64).sqrt()
        })
        .collect();

    for (n, e) in sizes.iter().zip(errors.iter()) {
        eprintln!("GPE smooth-S RMS error at N={}: {:.3e}", n, e);
    }
    let slope = (errors[errors.len() - 2] / errors[errors.len() - 1]).log2();
    eprintln!("final slope (64→128): {:.3}", slope);
    assert!(slope >= 1.7, "GPE slope = {:.3} (expected ≥ 1.7)", slope);
}

#[test]
fn force_is_small_on_mildly_perturbed_baseline_s_field() {
    // The baseline init uses a 2% perturbation; the GPE force should
    // match the analytic order-of-magnitude (~Ar · α² · k) and not
    // carry surprises.
    let nx = 64;
    let ny = 64;
    let dx = 1.0 / nx as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dx;
            s.set(i, j, 1.0 + 0.02 * ((TAU * x).sin() * (TAU * y).cos()));
        }
    }
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    // Use Ar = 0.1 (derived from default scales).
    GpeForce::with_ar(0.1).accumulate(
        &build_state(nx, ny, dx, &idx_x, &idx_y, &s),
        &mut VectorField { fx: &mut fx, fy: &mut fy },
    );
    let peak = fx.data().iter().chain(fy.data().iter()).fold(0.0_f64, |a, &v| a.max(v.abs()));
    // Expected ~ Ar · (1 + α) · α · k ~ 0.1·1·0.02·2π ~ 1.3e-2. Sanity bound ≤ 0.1.
    assert!(peak < 0.1, "GPE peak on mild field = {}, too large", peak);
    assert!(peak > 1e-4, "GPE peak on mild field = {}, suspiciously small", peak);
}
