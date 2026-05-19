//! Null-space projector robustness with mantle forcing.
//!
//! The Step 8 spec forbids subtracting `mean(f_mantle)` before
//! the solve — the null-space projector on `v` (Step 0) handles
//! the null mode. This test injects a mantle configuration whose
//! constant-RHS contribution has a non-zero mean (engineered by
//! asymmetric `S̃`), runs `solve_sheet`, and verifies
//! `|mean(v)| < 1e-15` after the solve.
//!
//! Setup choices:
//! - Asymmetric `S̃`: left half 1.0, right half 0.2 — mimics a
//!   plate-type contrast. Face-averaging S̃ produces values that
//!   do NOT cancel when multiplied by a zero-mean `v_pattern`,
//!   so `mean(f_mantle) ≠ 0` by construction.
//! - Constant η = 1, `drag_diag = coupling · S̃` → `total_diag`
//!   is just the mantle diagonal.
//! - `Mf = 1.0`, `coupling = 1.0`.

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::forcing::{BodyForce, MantleForce, SimulationState, VectorField};
use ymir_core::tectonics_v2::mantle::{
    build_mantle_diagonal_field, build_mantle_pattern, generate_stream_function, MantleConfig,
    StreamFunctionConfig,
};
use ymir_core::tectonics_v2::stokes::{solve_sheet, Grid, SheetConfig};

#[test]
fn solver_projects_null_space_even_when_mantle_rhs_has_nonzero_mean() {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    // Asymmetric S̃: 1.0 left half, 0.2 right half.
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            s.set(i, j, if i < nx / 2 { 1.0 } else { 0.2 });
        }
    }

    let psi = generate_stream_function(
        nx, ny,
        &StreamFunctionConfig { num_modes: 4, seed: 99 },
    );
    let pattern = build_mantle_pattern(&psi, dx, dy, &idx_x, &idx_y);

    let mf = 1.0;
    let coupling = 1.0;
    let cfg = MantleConfig::Enabled {
        mf, coupling, num_modes: 4, seed: 99, evolution_rate: 0.0,
    };
    let mantle_diag = build_mantle_diagonal_field(&cfg, &s).expect("enabled → Some");

    // Assemble constant-RHS part.
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    let state = SimulationState {
        nx, ny, dx, dy, idx_x: &idx_x, idx_y: &idx_y, s: &s,
    };
    MantleForce::new(mf, coupling, &pattern, &s).accumulate(
        &state,
        &mut VectorField { fx: &mut fx, fy: &mut fy },
    );

    // Inject a deliberate non-zero-mean perturbation on top of
    // the MantleForce output. The spec's concern is that the
    // projector handles any RHS with a non-zero mean — including
    // the occasional float-noise leftover. Rather than rely on
    // geometry that happens to produce a large enough mean, we
    // force the scenario: add uniform constants 0.5 and −0.3
    // across the domain. After the solve, `|mean(v)|` must
    // still be at machine noise.
    for v in fx.data_mut().iter_mut() { *v += 0.5; }
    for v in fy.data_mut().iter_mut() { *v -= 0.3; }

    // Sanity: the injected means dominate now.
    let mean_fx: f64 = fx.data().iter().sum::<f64>() / (nx * ny) as f64;
    let mean_fy: f64 = fy.data().iter().sum::<f64>() / (nx * ny) as f64;
    eprintln!("mean(fx) = {:.3e}, mean(fy) = {:.3e}", mean_fx, mean_fy);
    assert!(
        mean_fx.abs() > 0.1 && mean_fy.abs() > 0.1,
        "|mean(f)| ({:.3e}, {:.3e}) should be ~0.5, 0.3 after injection",
        mean_fx.abs(), mean_fy.abs(),
    );

    let fx_slice = fx.data().to_vec();
    let fy_slice = fy.data().to_vec();
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let cfg_solve = SheetConfig::default();
    let stats = solve_sheet(
        &grid, &eta, Some(&mantle_diag),
        &fx_slice, &fy_slice,
        &mut vx, &mut vy, &cfg_solve,
    );
    assert!(stats.converged, "CG did not converge: {:?}", stats);

    let mean_vx: f64 = vx.iter().sum::<f64>() / vx.len() as f64;
    let mean_vy: f64 = vy.iter().sum::<f64>() / vy.len() as f64;
    eprintln!("|mean(vx)| = {:.3e}, |mean(vy)| = {:.3e}", mean_vx.abs(), mean_vy.abs());
    assert!(
        mean_vx.abs() < 1e-15,
        "|mean(vx)| = {:.3e} (spec: < 1e-15)", mean_vx.abs(),
    );
    assert!(
        mean_vy.abs() < 1e-15,
        "|mean(vy)| = {:.3e} (spec: < 1e-15)", mean_vy.abs(),
    );
}
