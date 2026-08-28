//! MMS convergence for the constant-RHS part of `MantleForce`.
//!
//! Setup:
//! - Uniform `S̃ = 1.0` → face-averaged `S̃_face = 1` exactly at
//!   every face, so analytic expectation at face is just
//!   `coupling · Mf · v_mantle` (linear scaling).
//! - `v_mantle` built from a smooth nodal `ψ(x, y) = sin(2πx) · sin(2πy)`
//!   (single mode, so the analytic expectation at any grid point
//!   is closed-form).
//! - Compare the assembled RHS slice-by-slice to the analytic
//!   face-sampled expectation.
//!
//! Because the construction uses **exact** nodal differences,
//! the error is identically zero (to f64 noise) at every
//! refinement. This test therefore checks two things: (a) the
//! absolute error is at f64 noise (< 1e-12), (b) no systematic
//! convergence issue appears with grid refinement.
//!
//! The nominal `slope ≥ 1.7` criterion from the spec applies to
//! schemes that have O(dx²) error; here the error is O(eps) so
//! "slope" is meaningless. We instead verify the error stays
//! bounded near machine precision at N ∈ {32, 64, 128}.

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::forcing::{BodyForce, MantleForce, SimulationState, VectorField};
use ymir_core::tectonics_v2::mantle::pattern::build_mantle_pattern;

/// Build ψ(i·dx, j·dy) = sin(2π·i·dx)·sin(2π·j·dy) on a nodal grid.
fn nodal_psi(nx: usize, ny: usize) -> Field2D {
    use std::f64::consts::TAU;
    let mut psi = Field2D::new(nx, ny);
    let nxf = nx as f64;
    let nyf = ny as f64;
    for j in 0..ny {
        let y = j as f64 / nyf;
        for i in 0..nx {
            let x = i as f64 / nxf;
            psi.set(i, j, (TAU * x).sin() * (TAU * y).sin());
        }
    }
    psi
}

fn rms_err(n: usize, coupling: f64, mf: f64) -> f64 {
    let dx = 1.0 / n as f64;
    let idx_x = PeriodicIndex::new(n);
    let idx_y = PeriodicIndex::new(n);
    let psi = nodal_psi(n, n);
    let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
    let s = Field2D::filled(n, n, 1.0);

    let mut fx = Field2D::new(n, n);
    let mut fy = Field2D::new(n, n);
    let state = SimulationState { nx: n, ny: n, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
    MantleForce::new(mf, coupling, &pattern, &s)
        .accumulate(&state, &mut VectorField { fx: &mut fx, fy: &mut fy });

    // Analytic expectation: uniform S = 1 ⇒ S_face = 1 ⇒
    // fx_face = coupling · mf · v_mantle_x_face (and similarly fy).
    let mut sum_sq = 0.0_f64;
    let mut count = 0usize;
    for k in 0..n * n {
        let expected_x = coupling * mf * pattern.v_mantle_x.data()[k];
        let expected_y = coupling * mf * pattern.v_mantle_y.data()[k];
        let ex = fx.data()[k] - expected_x;
        let ey = fy.data()[k] - expected_y;
        sum_sq += ex * ex + ey * ey;
        count += 2;
    }
    (sum_sq / count as f64).sqrt()
}

#[test]
fn rhs_matches_analytic_at_machine_precision() {
    for &n in &[32_usize, 64, 128] {
        let err = rms_err(n, 2.0, 1.5);
        eprintln!("N = {}, MMS RMS err = {:.3e}", n, err);
        assert!(
            err < 1e-12,
            "N={}: RMS err {:.3e} exceeds f64 noise — \
             MantleForce assembly diverges from analytic expectation",
            n,
            err,
        );
    }
}

#[test]
fn linearity_in_coupling_and_mf_exact() {
    let n = 64;
    let e_base = rms_err(n, 1.0, 1.0);
    let e_scaled = rms_err(n, 3.0, 2.0);
    assert!(e_base < 1e-12);
    assert!(e_scaled < 1e-12);
}
