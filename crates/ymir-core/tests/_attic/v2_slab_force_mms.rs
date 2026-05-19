//! MMS for `SlabPullForce` assembly.
//!
//! Setup (Step 7 spec, Phase 2.2):
//! - `m_subducted(x, y) = sin(2π x)` at cell centre `(i+½, j+½)·dx`.
//! - `n̂_convergence = (1, 0)` uniform.
//! - Analytic: `f_slab(x, y) = Sp · sin(2π x) · (1, 0)`.
//!
//! The face discretisation of `SlabPullForce` places `fx` at the
//! left vertical face `(i·dx, (j+½)·dy)`. Face interpolation
//! `½(m[i-1,j] + m[i,j])` is a 2nd-order accurate midpoint
//! estimate of `m(i·dx, (j+½)·dy) = sin(2π · i·dx)`.
//!
//! Acceptance: convergence slope ≥ 1.7 at `N ∈ {32, 64, 128}`,
//! targeting the nominal 2.0 with headroom.

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::forcing::{BodyForce, SimulationState, SlabPullForce, VectorField};

fn rms_err_at_resolution(nx: usize) -> f64 {
    let ny = nx;
    let dx = 1.0 / nx as f64;
    let sp = 1.5;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    // S field — unused by SlabPullForce, but the SimulationState
    // struct still wants it.
    let s = Field2D::filled(nx, ny, 1.0);

    // m̃(x, y) = sin(2π x) at cell centre.
    let mut m = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x_cell = (i as f64 + 0.5) * dx;
            m.set(i, j, (2.0 * std::f64::consts::PI * x_cell).sin());
        }
    }

    // n̂ = (1, 0) uniform.
    let mut n_x = Field2D::new(nx, ny);
    for v in n_x.data_mut().iter_mut() {
        *v = 1.0;
    }
    let n_y = Field2D::new(nx, ny);

    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);

    let state = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
    SlabPullForce::new(sp, &m, &n_x, &n_y)
        .accumulate(&state, &mut VectorField { fx: &mut fx, fy: &mut fy });

    // Analytic fx at face (i·dx, (j+½)·dy): Sp · sin(2π · i·dx).
    // fy analytic is 0.
    let mut sum_sq_err = 0.0;
    let mut count = 0;
    for j in 0..ny {
        for i in 0..nx {
            let x_face = i as f64 * dx;
            let analytic_fx = sp * (2.0 * std::f64::consts::PI * x_face).sin();
            let numeric_fx = fx.data()[j * nx + i];
            sum_sq_err += (analytic_fx - numeric_fx).powi(2);
            sum_sq_err += fy.data()[j * nx + i].powi(2); // analytic fy = 0
            count += 2;
        }
    }
    (sum_sq_err / count as f64).sqrt()
}

#[test]
fn assembly_mms_second_order() {
    let resolutions = [32, 64, 128];
    let errs: Vec<f64> = resolutions.iter().map(|&n| rms_err_at_resolution(n)).collect();
    eprintln!("slab-force MMS errs (N ∈ {:?}): {:?}", resolutions, errs);

    for w in errs.windows(2) {
        let slope = (w[0] / w[1]).log2();
        assert!(
            slope >= 1.7,
            "MMS slope {:.3} < 1.7 between successive refinements (errs: {:?})",
            slope,
            errs,
        );
    }
}

/// Uniform `m` uniform `n̂ = (1, 0)` ⇒ face-averaged m ≡ m₀,
/// so fx = Sp · m₀ on every face, independent of the grid — the
/// exact arithmetic case (no discretisation error).
#[test]
fn uniform_fields_are_exact() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let sp = 2.0;
    let m0 = 0.37;

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let s = Field2D::filled(nx, ny, 1.0);

    let m = Field2D::filled(nx, ny, m0);
    let mut n_x = Field2D::new(nx, ny);
    for v in n_x.data_mut().iter_mut() {
        *v = 1.0;
    }
    let n_y = Field2D::new(nx, ny);

    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    let state = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
    SlabPullForce::new(sp, &m, &n_x, &n_y)
        .accumulate(&state, &mut VectorField { fx: &mut fx, fy: &mut fy });

    let expected = sp * m0;
    for &v in fx.data().iter() {
        assert!((v - expected).abs() < 1e-14, "fx = {}, expected {}", v, expected);
    }
    for &v in fy.data().iter() {
        assert_eq!(v, 0.0);
    }
}
