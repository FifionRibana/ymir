//! Validation of `div(v)_cell` and the conv/div split on the
//! staggered MAC grid (Step 5).
//!
//! - Analytic divergence: `v = (sin(2πx), -sin(2πy))` →
//!   `div(v) = 2π cos(2πx) - 2π cos(2πy)`.
//! - Expected order of convergence at refinement: 2 (bound set to
//!   ≥ 1.7, matching the rest of the milestone's MMS gates).
//! - Complementarity: `max(0, -div) + max(0, +div) = |div|` exactly.

use std::f64::consts::PI;
use ymir_core::tectonics_v2::boundaries::{convergent_component, div_v_cell};
use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};

fn run_divergence_at(nx: usize, ny: usize) -> f64 {
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let x_face = i as f64 * dx;
            let y_face = j as f64 * dy;
            vx[j * nx + i] = (2.0 * PI * x_face).sin();
            vy[j * nx + i] = -(2.0 * PI * y_face).sin();
        }
    }
    let mut div = Field2D::new(nx, ny);
    div_v_cell(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy, &mut div);
    // RMS error vs analytic at cell centres.
    let mut sum_sq = 0.0;
    let mut n = 0usize;
    for j in 0..ny {
        for i in 0..nx {
            let xc = (i as f64 + 0.5) * dx;
            let yc = (j as f64 + 0.5) * dy;
            let expected = 2.0 * PI * (2.0 * PI * xc).cos() - 2.0 * PI * (2.0 * PI * yc).cos();
            let got = div.get(i, j);
            sum_sq += (got - expected).powi(2);
            n += 1;
        }
    }
    (sum_sq / n as f64).sqrt()
}

#[test]
fn divergence_converges_at_order_at_least_1_7() {
    let e32 = run_divergence_at(32, 32);
    let e64 = run_divergence_at(64, 64);
    let e128 = run_divergence_at(128, 128);
    let slope_32_to_64 = (e32 / e64).log2();
    let slope_64_to_128 = (e64 / e128).log2();
    println!(
        "div MMS: e32={} e64={} e128={} slopes {} {}",
        e32, e64, e128, slope_32_to_64, slope_64_to_128
    );
    assert!(
        slope_64_to_128 >= 1.7,
        "final slope {} below 1.7 bound (e64={}, e128={})",
        slope_64_to_128,
        e64,
        e128,
    );
}

#[test]
fn convergent_and_divergent_components_sum_to_absolute_divergence() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            vx[j * nx + i] = (2.0 * PI * i as f64 * dx).sin();
            vy[j * nx + i] = -(2.0 * PI * j as f64 * dy).sin();
        }
    }
    let mut div = Field2D::new(nx, ny);
    div_v_cell(nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy, &mut div);
    let mut conv = Field2D::new(nx, ny);
    let mut divg = Field2D::new(nx, ny);
    convergent_component(&div, &mut conv);
    // `divergent_component` is symmetric: max(0, +div). Implement
    // inline from the primitive.
    for (out, &d) in divg.data_mut().iter_mut().zip(div.data().iter()) {
        *out = d.max(0.0);
    }
    for ((&c, &dv), &d) in conv.data().iter().zip(divg.data().iter()).zip(div.data().iter()) {
        let sum = c + dv;
        assert!((sum - d.abs()).abs() < 1e-14, "conv+div != |div|");
    }
}
