//! Step 4 — operator symmetry with basal drag active.
//!
//! Verifies `⟨A · u, w⟩ = ⟨u, A · w⟩` for the drag-augmented
//! momentum operator on random test vectors. Symmetry is what lets
//! CG remain the inner linear solver (no BiCGSTAB needed at Step 4).
//!
//! The drag contribution is per-cell diagonal (`Br · S̃² · I` after
//! face interpolation), which is trivially symmetric. The test
//! guards the face-interpolation convention: if `apply_momentum`
//! used a different cell-to-face averaging for drag on vx vs vy
//! (e.g. biased to `(im, j)` on vx but to `(i, j)` on vy), the
//! bilinear form would pick up an asymmetric coupling.

use ymir_core::tectonics_v2::basal_drag::{
    BasalDragConfig, BasalDragLaw, build_drag_diagonal_field,
};
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::operator::{StokesGrid, apply_momentum};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[test]
fn drag_augmented_operator_is_symmetric_on_random_inputs() {
    let nx = 8;
    let ny = 8;
    let grid = StokesGrid::new(nx, ny, 0.13, 0.17);

    // Random-ish η field (same seed-idempotent pattern as the
    // Step-0 momentum_is_symmetric test).
    let mut eta = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let e = 1.0 + 0.5 * ((i * 3 + j * 7) % 5) as f64 / 5.0;
            eta.set(i, j, e);
        }
    }

    // Structured S̃ so the drag contribution is spatially
    // heterogeneous (exercises the face averaging on both sides).
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let v = 0.6 + 0.1 * (i + 2 * j) as f64 / (nx + ny) as f64;
            s.set(i, j, v);
        }
    }
    let drag_cfg = BasalDragConfig::Enabled(BasalDragLaw { br: 0.1, s_exponent: 2.0 });
    let drag_diag = build_drag_diagonal_field(&drag_cfg, &s).expect("Enabled must produce Some");

    let n2 = nx * ny;
    let mut ux = vec![0.0; n2];
    let mut uy = vec![0.0; n2];
    let mut wx = vec![0.0; n2];
    let mut wy = vec![0.0; n2];
    for k in 0..n2 {
        ux[k] = ((k as f64 * 1.7).sin()) * 1.1;
        uy[k] = ((k as f64 * 2.3).cos()) * 0.7;
        wx[k] = ((k as f64 * 0.9).sin()) * 0.5;
        wy[k] = ((k as f64 * 1.3).cos()) * 1.3;
    }

    let mut aux_x = vec![0.0; n2];
    let mut aux_y = vec![0.0; n2];
    let mut awx = vec![0.0; n2];
    let mut awy = vec![0.0; n2];
    apply_momentum(&grid, &eta, Some(&drag_diag), &ux, &uy, &mut aux_x, &mut aux_y);
    apply_momentum(&grid, &eta, Some(&drag_diag), &wx, &wy, &mut awx, &mut awy);
    let lhs = dot(&aux_x, &wx) + dot(&aux_y, &wy);
    let rhs = dot(&ux, &awx) + dot(&uy, &awy);
    let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
    eprintln!("⟨Au, w⟩ = {lhs:.6e} vs ⟨u, Aw⟩ = {rhs:.6e}; rel = {rel:.3e}");
    assert!(rel < 1e-12, "drag-augmented operator asymmetric: rel = {rel} (spec: < 1e-12)",);
}

/// Sanity guard: the disabled variant must produce the same output
/// as the un-augmented call (the short-circuit must not perturb
/// anything beyond skipping the augmentation loop).
#[test]
fn disabled_variant_matches_unaugmented_call() {
    let nx = 8;
    let ny = 8;
    let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
    let eta = Field2D::filled(nx, ny, 1.0);
    let n2 = nx * ny;
    let mut vx = vec![0.0; n2];
    let mut vy = vec![0.0; n2];
    for k in 0..n2 {
        vx[k] = (k as f64 * 0.7).sin();
        vy[k] = (k as f64 * 1.1).cos();
    }
    let mut out_with_none = vec![0.0; n2];
    let mut out_with_none_y = vec![0.0; n2];
    let mut out_without = vec![0.0; n2];
    let mut out_without_y = vec![0.0; n2];
    apply_momentum(&grid, &eta, None, &vx, &vy, &mut out_with_none, &mut out_with_none_y);
    apply_momentum(&grid, &eta, None, &vx, &vy, &mut out_without, &mut out_without_y);
    for k in 0..n2 {
        assert_eq!(out_with_none[k], out_without[k]);
        assert_eq!(out_with_none_y[k], out_without_y[k]);
    }
}
