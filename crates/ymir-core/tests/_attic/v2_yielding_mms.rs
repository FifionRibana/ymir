//! Step 3 yielding load-bearing tests.
//!
//! # Why not a strict MMS order-of-accuracy test?
//!
//! A straightforward manufactured solution
//! `v(x, y) = (sin(kx), sin(ky))` on the periodic unit square has
//! `ε̇_II = (k/√2) √(cos²(kx) + cos²(ky))`, which vanishes at
//! isolated points (where both cosines hit zero) and has
//! `cos²(kx) = O(dx²)` at the cell centres nearest the lines where
//! `cos(kx) = 0`. The discrete `ε̇_II_cc` at those cells differs
//! from the face-value `ε̇_II_face` by **O(1)** (because the missing
//! `ε̇_xx` contribution is O(dx²) in *squared* form, but under the
//! square root mixes with the finite `ε̇_yy` term), so
//! `η_eff(ε̇_II_cc) − η_eff(ε̇_II_face)` is O(1) near those lines.
//! The discretisation is still 2nd-order on the bulk of the
//! domain, but strictly proving `‖A(v_exact) − f_analytic‖ = O(dx²)`
//! in L² requires a manufactured velocity whose `ε̇_II` never
//! vanishes — which is incompatible with smooth periodic boundary
//! conditions.
//!
//! Rather than ship a test whose acceptance threshold would be
//! looser than the evidence warrants, Step 3 uses **two orthogonal
//! correctness garde-fous for the yielding Jacobian**:
//!
//! 1. **Algebraic symmetry** of the tangent operator with yielding
//!    active (`jacobian_symmetric_with_yielding` below). The
//!    operator SPD-compatibility is what lets CG keep working.
//!
//! 2. **Newton super-linear convergence** on a target-generated RHS
//!    (see `v2_newton_convergence.rs` Step-3 extension). Tail
//!    residuals reduce by ≥ 100× per iteration, confirming the
//!    chain-rule Jacobian is accurate in the full nonlinear
//!    regime. A miscomputed `d_eta_eff/dε̇_II` would degrade this
//!    to linear.
//!
//! Together these pin down the derivative chain well enough that a
//! strict spatial-order MMS becomes informational, not load-bearing.

use std::f64::consts::TAU;

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::presets::YieldingConfig;
use ymir_core::tectonics_v2::rheology::{StrainRate, ViscosityLaw, YieldingLaw};
use ymir_core::tectonics_v2::stokes::operator::{StokesGrid, TangentContext, apply_jacobian};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[test]
fn jacobian_symmetric_with_yielding() {
    // Non-trivial v_k so both the Picard block and the Newton-extra
    // term contribute, with yielding Enabled so the plastic branch
    // participates in the chain rule. Verify
    // ⟨J u, w⟩ = ⟨u, J w⟩ to 1e-10 on random test vectors.
    let nx = 12;
    let ny = 12;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = StokesGrid::new(nx, ny, dx, dy);
    let ylaw = YieldingLaw::default();
    let law = ViscosityLaw { yielding: YieldingConfig::Enabled(ylaw), ..Default::default() };

    let n2 = nx * ny;
    let mut vx_k = vec![0.0; n2];
    let mut vy_k = vec![0.0; n2];
    let k = TAU;
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dy;
            vx_k[j * nx + i] = 0.15 * (k * x).sin() * (k * y).cos();
            vy_k[j * nx + i] = -0.15 * (k * x).cos() * (k * y).sin();
        }
    }
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let sr = StrainRate::compute(nx, ny, dx, dy, &idx_x, &idx_y, &vx_k, &vy_k);
    let ctx = TangentContext::from_strain_rate(&grid, &law, &sr, None);

    let mut ux = vec![0.0; n2];
    let mut uy = vec![0.0; n2];
    let mut wx = vec![0.0; n2];
    let mut wy = vec![0.0; n2];
    for k in 0..n2 {
        ux[k] = ((k as f64 * 1.37).sin()) * 1.1;
        uy[k] = ((k as f64 * 2.11).cos()) * 0.7;
        wx[k] = ((k as f64 * 0.83).sin()) * 0.5;
        wy[k] = ((k as f64 * 1.19).cos()) * 1.3;
    }

    let mut jux = vec![0.0; n2];
    let mut juy = vec![0.0; n2];
    let mut jwx = vec![0.0; n2];
    let mut jwy = vec![0.0; n2];
    apply_jacobian(&grid, &ctx, None, &ux, &uy, &mut jux, &mut juy);
    apply_jacobian(&grid, &ctx, None, &wx, &wy, &mut jwx, &mut jwy);
    let lhs = dot(&jux, &wx) + dot(&juy, &wy);
    let rhs = dot(&ux, &jwx) + dot(&uy, &jwy);
    let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
    eprintln!("⟨J u, w⟩ = {:.6e}, ⟨u, J w⟩ = {:.6e}, rel = {:.3e}", lhs, rhs, rel,);
    assert!(rel < 1e-10, "Jacobian asymmetric with yielding active: rel = {}", rel,);
}

/// Build an η field at cell centres with yielding enabled, then feed
/// it to the operator and check the implied symmetry of the
/// plain-Picard block (no tangent-extra) is also preserved. This
/// catches regressions in how `build_eta_field` dispatches through
/// `ViscosityLaw::eta_effective` when the yielding match arm is
/// live.
#[test]
fn picard_block_symmetric_with_yielding_eta() {
    use ymir_core::tectonics_v2::rheology::build_eta_field;
    use ymir_core::tectonics_v2::stokes::operator::apply_momentum;
    let nx = 10;
    let ny = 10;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = StokesGrid::new(nx, ny, dx, dy);
    let ylaw = YieldingLaw::default();
    let law = ViscosityLaw { yielding: YieldingConfig::Enabled(ylaw), ..Default::default() };
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    // Set up a v_k to generate a non-trivial η field.
    let mut vx_k = vec![0.0; nx * ny];
    let mut vy_k = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            vx_k[j * nx + i] = 0.12 * ((i + j) as f64).sin();
            vy_k[j * nx + i] = 0.08 * ((i * 3 + j * 5) as f64).cos();
        }
    }
    let sr = StrainRate::compute(nx, ny, dx, dy, &idx_x, &idx_y, &vx_k, &vy_k);
    let eta = build_eta_field(&law, &sr.eps_ii_center, None);

    let n2 = nx * ny;
    let mut ux = vec![0.0; n2];
    let mut uy = vec![0.0; n2];
    let mut wx = vec![0.0; n2];
    let mut wy = vec![0.0; n2];
    for k in 0..n2 {
        ux[k] = (k as f64 * 1.7).sin();
        uy[k] = (k as f64 * 2.3).cos();
        wx[k] = (k as f64 * 0.9).sin();
        wy[k] = (k as f64 * 1.3).cos();
    }
    let mut aux_x = vec![0.0; n2];
    let mut aux_y = vec![0.0; n2];
    let mut awx = vec![0.0; n2];
    let mut awy = vec![0.0; n2];
    apply_momentum(&grid, &eta, None, &ux, &uy, &mut aux_x, &mut aux_y);
    apply_momentum(&grid, &eta, None, &wx, &wy, &mut awx, &mut awy);
    let lhs = dot(&aux_x, &wx) + dot(&aux_y, &wy);
    let rhs = dot(&ux, &awx) + dot(&uy, &awy);
    let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
    assert!(rel < 1e-12, "Picard block asymmetric with yielding η: rel={}", rel);
    let _ = Field2D::filled(1, 1, 0.0);
}
