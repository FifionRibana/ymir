//! Null-space preservation with `mean(f_slab) ≠ 0`.
//!
//! The Step 7 spec explicitly forbids subtracting the mean of
//! `f_slab` before the Stokes solve — the null-space projector
//! installed in the preconditioner operates on `v`, not `f`, and
//! removing the mean preemptively would mask physical information.
//!
//! This test builds an asymmetric slab configuration that leaves
//! `mean(f_slab) ≠ 0`, runs `solve_sheet`, and checks that the
//! solved velocity has null mean components to `< 1e-15`. A
//! successful pass confirms the documented design: the projector
//! handles the null space correctly even when the driving force
//! has a net DC component.

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::forcing::{BodyForce, SimulationState, SlabPullForce, VectorField};
use ymir_core::tectonics_v2::stokes::{Grid, SheetConfig, solve_sheet};

#[test]
fn solver_projects_null_space_even_when_slab_force_has_nonzero_mean() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let s = Field2D::filled(nx, ny, 1.0);

    // Asymmetric m field: strictly positive on the left half, zero
    // on the right half. Combined with a uniform n̂ = (1, 0), this
    // delivers a force field with clearly non-zero x-mean.
    let mut m = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx / 2 {
            m.set(i, j, 0.6);
        }
    }
    let mut n_x = Field2D::new(nx, ny);
    for v in n_x.data_mut().iter_mut() {
        *v = 1.0;
    }
    let mut n_y = Field2D::new(nx, ny);
    // Add a constant y bias to also test y null-space.
    for v in n_y.data_mut().iter_mut() {
        *v = 0.3;
    }

    // Assemble f_slab into fx, fy.
    let mut fx_field = Field2D::new(nx, ny);
    let mut fy_field = Field2D::new(nx, ny);
    let state = SimulationState { nx, ny, dx, dy, idx_x: &idx_x, idx_y: &idx_y, s: &s };
    SlabPullForce::new(1.5, &m, &n_x, &n_y)
        .accumulate(&state, &mut VectorField { fx: &mut fx_field, fy: &mut fy_field });

    // Sanity check: mean(f) is indeed non-zero before the solve —
    // otherwise we'd be testing a trivial case.
    let mean_fx: f64 = fx_field.data().iter().sum::<f64>() / (nx * ny) as f64;
    let mean_fy: f64 = fy_field.data().iter().sum::<f64>() / (nx * ny) as f64;
    assert!(
        mean_fx.abs() > 1e-6,
        "|mean(fx)| = {:.3e} — test setup does not deliver a non-zero-mean force",
        mean_fx.abs(),
    );
    assert!(
        mean_fy.abs() > 1e-6,
        "|mean(fy)| = {:.3e} — test setup does not deliver a non-zero-mean force",
        mean_fy.abs(),
    );
    eprintln!("mean(fx) = {:.3e}, mean(fy) = {:.3e}", mean_fx, mean_fy);

    // Pass raw slices — the solver's internal `nullspace::project_velocity`
    // is responsible for DC removal from the RHS before CG.
    let fx = fx_field.data().to_vec();
    let fy = fy_field.data().to_vec();
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
    eprintln!("|mean(vx)| = {:.3e}, |mean(vy)| = {:.3e}", mvx, mvy);
    assert!(mvx < 1e-15, "|mean(vx)| = {:.3e} (spec: < 1e-15)", mvx,);
    assert!(mvy < 1e-15, "|mean(vy)| = {:.3e} (spec: < 1e-15)", mvy,);
}
