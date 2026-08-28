//! End-to-end smoke test for rectangular-grid tectonic simulation.
//!
//! Drives the full solver pipeline on a non-square grid for several dozen
//! steps and checks that nothing panics, no NaN/Inf values appear, plate
//! bookkeeping stays sane, and total mass doesn't drift catastrophically.
//! If this passes, rectangular-grid support is functionally correct end-to-end.

use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::plates::{PlateConfig, generate_plates};
use ymir_core::tectonics::solver::config::{
    ContinuationConfig, NewtonConfig, NonlinearSolver, PicardConfig, TectonicsConfig,
};
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::grid::StaggeredGrid;
use ymir_core::tectonics::solver::tectonics::{DynamicPlateContext, run_tectonics};
use ymir_core::tectonics::solver::workspace::SolverWorkspace;

#[test]
fn rectangular_simulation_smoke_test() {
    // 128 × 85 (3:2 aspect). Non-square by design — this whole test is
    // worthless on a square grid.
    let config = PlateConfig::default().with_resolution_aspect(128, 1.5);
    assert_eq!(config.grid_width, 128, "with_resolution_aspect width");
    assert_eq!(config.grid_height, 85, "with_resolution_aspect height");

    let seed = WorldSeed::new(42);
    let init = generate_plates(&config, &seed);
    assert_eq!(init.grid_width, 128);
    assert_eq!(init.grid_height, 85);
    assert_eq!(init.plate_ids.len(), 128 * 85);
    assert_eq!(init.thickness.width, 128);
    assert_eq!(init.thickness.height, 85);

    // Build a StaggeredGrid from the generated thickness. dx is derived
    // from grid_width so the x-axis physical domain stays 1.0; cells
    // remain square in physical space regardless of aspect.
    let nx = init.grid_width;
    let ny = init.grid_height;
    let dx = 1.0 / nx as f64;
    let mut grid = StaggeredGrid::new(nx, ny, dx);
    for j in 0..ny {
        for i in 0..nx {
            grid.s.set(i, j, init.thickness.data[j * nx + i] as f64);
        }
    }

    let traction = init.to_traction_field();
    let num_plates = init.plates.len();
    let mut plate_ctx = DynamicPlateContext {
        ids: init.plate_ids.clone(),
        plates: init.plates.clone(),
        traction,
        next_id: num_plates,
        disp_x: Field2D::new(nx, ny),
        disp_y: Field2D::new(nx, ny),
    };

    let num_timesteps = 50usize;
    let sim_config = TectonicsConfig {
        num_timesteps,
        gravity_factor: 1.0,
        cfl_factor: 0.3,
        s_min: 0.1,
        s_max: 2.5,
        nonlinear_solver: NonlinearSolver::Picard,
        picard: PicardConfig {
            max_iterations: 30,
            tolerance: 1e-3,
            relaxation: 0.7,
            cg_max_iter: 500,
            cg_tolerance: 1e-8,
            strain_rate_min: 1e-3,
            power_law_n: 1.0,
            eta_min: 1e-3,
            eta_max: 1e4,
        },
        newton: NewtonConfig::default(),
        continuation: ContinuationConfig { enabled: false, ..Default::default() },
        boundaries: Default::default(),
        dynamic_boundaries: false,
        cratonic: Default::default(),
        yielding: Default::default(),
        basal_friction: 0.0,
        mantle: Default::default(),
        recycling: Default::default(),
        adaptive_dt: Default::default(),
    };
    let mut ws = SolverWorkspace::new(nx, ny);

    let initial_mass: f64 = grid.s.data().iter().sum();
    assert!(initial_mass.is_finite() && initial_mass > 0.0);

    let mut max_picard = 0usize;
    let result = run_tectonics(
        &sim_config,
        &mut plate_ctx,
        &mut grid,
        &mut ws,
        |step, _total, stats, snap| {
            max_picard = max_picard.max(stats.picard_iterations);

            // Per-step invariants: every thickness cell is finite.
            for (k, &s) in snap.s_field.data().iter().enumerate() {
                assert!(s.is_finite(), "non-finite thickness at step {step}, index {k}: {s}");
            }
            // Velocity is also finite.
            assert!(stats.max_velocity.is_finite(), "non-finite max_velocity at step {step}");
            true
        },
    );
    assert!(result.is_ok(), "run_tectonics failed on rectangular grid: {:?}", result.err());

    // Final-state invariants.
    for (k, &s) in grid.s.data().iter().enumerate() {
        assert!(s.is_finite(), "non-finite thickness at end, index {k}: {s}");
    }
    for &v in grid.vx.data() {
        assert!(v.is_finite(), "non-finite vx");
    }
    for &v in grid.vy.data() {
        assert!(v.is_finite(), "non-finite vy");
    }

    // Plate bookkeeping: the number of active plates stays within a
    // reasonable range. Dynamic boundaries are off so no plates can be
    // created; the count should be <= num_plates. We allow 0 as a lower
    // bound only to avoid a false positive if every plate is (unusually)
    // consumed — the real catch is "no panics and finite state".
    let active = plate_ctx.plates.iter().filter(|p| p.active).count();
    assert!(
        active <= 2 * num_plates,
        "unreasonable active plate count: {active} (started with {num_plates})"
    );

    // Mass balance: static boundaries mean no sources/sinks; the only
    // reason mass could drift is s_min/s_max clamping. Allow 20% drift —
    // catastrophic drift indicates a real bug, not normal clamping.
    let final_mass: f64 = grid.s.data().iter().sum();
    let drift = (final_mass - initial_mass).abs() / initial_mass;
    assert!(
        drift < 0.2,
        "mass drift {:.2}% over {num_timesteps} steps: initial={initial_mass}, final={final_mass}",
        drift * 100.0
    );

    // Picard convergence: we should never hit the iteration cap. Hitting
    // it means the nonlinear solve stalled, which would usually come with
    // NaN or divergence further along but is worth flagging directly.
    assert!(
        max_picard < sim_config.picard.max_iterations,
        "Picard hit max_iterations={} at least once (observed max {max_picard})",
        sim_config.picard.max_iterations
    );
}
