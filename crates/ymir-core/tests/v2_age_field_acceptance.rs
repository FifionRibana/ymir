//! Step 10 acceptance tests for the geological age field.
//!
//! Five named tests per the issue's "Definition of done":
//!
//! - `v2_age_field_initialization`: continental cells start at
//!   `continental_age_init`, oceanic cells at `oceanic_age_init`.
//! - `v2_age_field_bounds`: `age_min ≥ 0` and
//!   `age_max ≤ age_init_max + simulation_time` at every step
//!   (acceptance #1).
//! - `v2_age_field_advection_mms`: prescribe a smooth `A(x, 0)` +
//!   uniform `v`, verify analytic transport at the expected
//!   convergence order (acceptance #2, slope ≥ 0.95).
//! - `v2_age_field_quiescent_growth`: with `v = 0` (no events),
//!   `A` grows linearly at `dt` per step (acceptance #3).
//! - `v2_age_field_ridge_reset`: ridge cells drop to `0` after the
//!   reset event (acceptance #4).
//! - `v2_age_field_collision_max`: collision cells receive
//!   `max(neighbour ages)` (acceptance #5).
//!
//! These tests run as fast unit-style cases (no full harness
//! invocation) where possible, hitting the public API of the
//! `age_field` module directly. The harness-level acceptance is
//! verified by the Step 10 physics test
//! (`v2_step10_physics_and_regression`) at 64² × 100 steps.

use ymir_core::tectonics_v2::age_field::advection::{step_age_advect, QUIESCENT_GROWTH_RATE};
use ymir_core::tectonics_v2::age_field::events::apply_age_events;
use ymir_core::tectonics_v2::age_field::{AgeFieldConfigEnabled, AgeFieldState};
use ymir_core::tectonics_v2::boundaries::{BoundaryFlag, BoundaryFlagField, PlateType, PlateTypeField};
use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};

#[test]
fn v2_age_field_initialization() {
    // Build a small S̃ field with a 2x2 continental block in a sea
    // of oceanic cells. AgeFieldState::from_initial_thickness
    // should classify per the `S̃ > 0.5` threshold.
    let nx = 6;
    let ny = 6;
    let mut s = Field2D::filled(nx, ny, 0.2); // oceanic
    for j in 2..4 {
        for i in 2..4 {
            s.set(i, j, 1.0); // continental
        }
    }
    let cfg = AgeFieldConfigEnabled::default();
    let state = AgeFieldState::from_initial_thickness(&s, &cfg);
    for j in 0..ny {
        for i in 0..nx {
            let in_continental = (2..4).contains(&i) && (2..4).contains(&j);
            let expected = if in_continental {
                cfg.continental_age_init
            } else {
                cfg.oceanic_age_init
            };
            assert_eq!(
                state.current.get(i, j),
                expected,
                "cell ({},{}) age: got {}, expected {}",
                i, j, state.current.get(i, j), expected
            );
        }
    }
}

#[test]
fn v2_age_field_bounds() {
    // Acceptance #1: `age_min ≥ 0` and
    // `age_max ≤ age_init_max + simulation_time` at every step,
    // for a sequence of advection + event applications.
    //
    // We drive the system manually with a non-trivial velocity
    // (sinusoidal, divergence-non-zero) and a couple of boundary
    // events, then assert the bound holds after each step.
    let nx = 8;
    let ny = 8;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let cfg = AgeFieldConfigEnabled::default();
    let s_init = Field2D::filled(nx, ny, 1.0); // continental
    let mut state = AgeFieldState::from_initial_thickness(&s_init, &cfg);

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    for k in 0..nx * ny {
        vx[k] = ((k as f64 * 0.7).sin()) * 0.3;
        vy[k] = ((k as f64 * 1.3).cos()) * 0.2;
    }
    let dt = 0.05_f64;
    let n_steps = 20usize;
    let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
    // Plant a single ridge at (3, 3) and a collision at (5, 5).
    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    flags.set(3, 3, BoundaryFlag::Rift);
    flags.set(5, 5, BoundaryFlag::ContinentalCollision);

    let age_init_max = state.age_init_max();

    for step in 0..n_steps {
        step_age_advect(
            nx, ny, 1.0, 1.0, dt, &idx_x, &idx_y,
            &state.current, &vx, &vy, &mut state.next,
        );
        std::mem::swap(&mut state.current, &mut state.next);
        let _ = apply_age_events(&flags, &plate_type, &s_init, &idx_x, &idx_y, &mut state.current);

        let amin = state.current.data().iter().cloned().fold(f64::INFINITY, f64::min);
        let amax = state.current.data().iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        // Bound: age_max ≤ age_init_max + simulation_time so far.
        // We add a small ε to allow for numerical advection
        // contributions from upwind borrowing, which can transport
        // age forward but never above the global supremum.
        let elapsed = (step + 1) as f64 * dt;
        assert!(amin >= 0.0, "step {}: age_min = {} < 0", step, amin);
        assert!(
            amax <= age_init_max + elapsed + 1e-12,
            "step {}: age_max = {} exceeds age_init_max + elapsed = {}",
            step, amax, age_init_max + elapsed
        );
    }
}

#[test]
fn v2_age_field_advection_mms() {
    // Acceptance #2: smooth `A(x, 0) = sin(2π x)`, uniform vx = 1
    // (vy = 0), advect for `T` and compare the final state to the
    // analytic Lagrangian solution `A(x, T) = sin(2π (x - T)) + T`
    // (the `+ T` is the quiescent-growth contribution since
    // QUIESCENT_GROWTH_RATE = 1). Convergence order is first-
    // order upwind — slope ≥ 0.95 across grid refinement.
    fn run(nx: usize) -> f64 {
        let ny = 4;
        let dx = 1.0 / nx as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut a = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * dx;
                a.set(i, j, (2.0 * std::f64::consts::PI * x).sin());
            }
        }
        let vx = vec![1.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        // CFL-safe `dt`. Take total time = 0.1 (one tenth of a
        // wave period at v = 1 on the unit torus).
        let total_t = 0.1_f64;
        let dt_max = 0.4 * dx; // CFL = 0.4
        let n_steps = (total_t / dt_max).ceil() as usize;
        let dt = total_t / n_steps as f64;
        let mut a_next = Field2D::new(nx, ny);
        for _ in 0..n_steps {
            step_age_advect(nx, ny, dx, dx, dt, &idx_x, &idx_y, &a, &vx, &vy, &mut a_next);
            std::mem::swap(&mut a, &mut a_next);
        }
        // Analytic: A(x, T) = sin(2π(x - T)) + T·QUIESCENT_GROWTH_RATE
        // Compute RMS error.
        let mut sse = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * dx;
                let analytic =
                    (2.0 * std::f64::consts::PI * (x - total_t)).sin() + total_t * QUIESCENT_GROWTH_RATE;
                let err = a.get(i, j) - analytic;
                sse += err * err;
            }
        }
        (sse / (nx * ny) as f64).sqrt()
    }

    let nxs = [16, 32, 64, 128];
    let errs: Vec<f64> = nxs.iter().map(|&n| run(n)).collect();
    eprintln!("MMS errs by N: {:?}", errs);
    // First-order upwind: error halves when grid doubles.
    // Slope = log2(err_n / err_{n+1}). Acceptance #2: ≥ 0.95.
    for w in errs.windows(2) {
        let slope = (w[0] / w[1]).log2();
        eprintln!("  slope between successive Ns: {:.3}", slope);
        assert!(
            slope >= 0.95,
            "MMS convergence slope {} < 0.95 (first-order upwind expected)",
            slope
        );
    }
}

#[test]
fn v2_age_field_quiescent_growth() {
    // Acceptance #3: with v = 0 and no boundary events, A grows
    // by exactly `dt` per step starting from the initial value.
    let nx = 4;
    let ny = 4;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let cfg = AgeFieldConfigEnabled { continental_age_init: 5.0, oceanic_age_init: 0.0 };
    let s = Field2D::filled(nx, ny, 1.0); // continental
    let mut state = AgeFieldState::from_initial_thickness(&s, &cfg);

    let vx = vec![0.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    let dt = 0.05_f64;
    let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
    let flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);

    for step in 1..=10 {
        step_age_advect(
            nx, ny, 1.0, 1.0, dt, &idx_x, &idx_y,
            &state.current, &vx, &vy, &mut state.next,
        );
        std::mem::swap(&mut state.current, &mut state.next);
        let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut state.current);
        assert_eq!(counts.ridge_resets, 0);
        assert_eq!(counts.arc_resets, 0);
        assert_eq!(counts.collision_max_events, 0);

        let expected = 5.0 + (step as f64) * dt;
        for v in state.current.data() {
            assert!(
                (v - expected).abs() < 1e-12,
                "step {}: cell got {}, expected {}",
                step, v, expected
            );
        }
    }
}

#[test]
fn v2_age_field_ridge_reset() {
    // Acceptance #4: at known ridge cells, A drops to 0 after the
    // reset event. (Detailed coverage is in the events.rs unit
    // tests; this test exercises the harness-shape pipeline:
    // advect with v = 0, then apply ridge → expect A = 0 at the
    // ridge cell, A = unchanged elsewhere.)
    let nx = 5;
    let ny = 5;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let cfg = AgeFieldConfigEnabled { continental_age_init: 4.0, oceanic_age_init: 1.0 };
    let s = Field2D::filled(nx, ny, 0.2); // oceanic
    let mut state = AgeFieldState::from_initial_thickness(&s, &cfg);

    let vx = vec![0.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    let dt = 0.1_f64;
    let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    flags.set(2, 2, BoundaryFlag::Rift);

    step_age_advect(
        nx, ny, 1.0, 1.0, dt, &idx_x, &idx_y,
        &state.current, &vx, &vy, &mut state.next,
    );
    std::mem::swap(&mut state.current, &mut state.next);
    // Before events: every cell at oceanic_age_init + dt = 1.1.
    let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut state.current);
    assert_eq!(counts.ridge_resets, 1);
    assert_eq!(state.current.get(2, 2), 0.0, "ridge cell not reset");
    // Other cells: 1.1 (1.0 init + 0.1 quiescent growth).
    for j in 0..ny {
        for i in 0..nx {
            if (i, j) == (2, 2) {
                continue;
            }
            assert!(
                (state.current.get(i, j) - 1.1).abs() < 1e-12,
                "non-ridge cell ({},{}) age {} != 1.1",
                i, j, state.current.get(i, j)
            );
        }
    }
}

#[test]
fn v2_age_field_collision_max() {
    // Acceptance #5: at known collision cells, A receives the max
    // age over the contributing continental neighbours. We arrange
    // a 1D-like setup with a collision cell flanked by neighbours
    // of clearly different ages and verify the max wins.
    let nx = 5;
    let ny = 1;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut a = Field2D::new(nx, ny);
    a.set(0, 0, 2.0);
    a.set(1, 0, 8.0); // west neighbour — oldest
    a.set(2, 0, 1.0); // collision cell
    a.set(3, 0, 4.0);
    a.set(4, 0, 0.0);

    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    flags.set(2, 0, BoundaryFlag::ContinentalCollision);
    let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
    let s = Field2D::filled(nx, ny, 1.0); // continental everywhere

    let counts = apply_age_events(&flags, &plate_type, &s, &idx_x, &idx_y, &mut a);
    assert_eq!(counts.collision_max_events, 1);
    assert_eq!(a.get(2, 0), 8.0, "collision should pick the oldest neighbour");
    assert!((counts.collision_max_age_mean() - 8.0).abs() < 1e-12);
}
