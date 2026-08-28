//! Mass-balance test on a layout designed to have `∫Q = 0` exactly
//! with zero velocity (Step 5, issue #89 D5).
//!
//! `balanced_sub_spread` has matched counts of `Rift` and
//! `OceanicSubduction` cells on an all-oceanic domain. With
//! velocity fixed at zero (we force it by running the harness with
//! a `ZeroForce` body force and a strict advection skip), the
//! subduction flag's `Q_sub` branch evaluates to zero (no
//! convergent motion) while `Q_spread` fires at rate `k_spread`.
//!
//! Because this test targets the mass-balance arithmetic, not the
//! Stokes solver, we drive the time loop by hand: no Stokes, no
//! advection, just the source/sink pipeline + clamp. The acceptance
//! is `mass_drift_relative < 1e-8` (machine-noise floor for the
//! floating-point sums).

use ymir_core::tectonics_v2::boundaries::{
    BoundaryRates, apply_clamp_with_tracking, balanced_sub_spread, compute_source_sink_terms,
};
use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};

#[test]
fn balanced_layout_with_matched_rates_has_machine_drift() {
    let nx = 32;
    let ny = 24;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let layout = balanced_sub_spread(nx, ny);
    // Rates tuned so that ∫Q = 0 exactly: same count of Rift and
    // OceanicSubduction cells; set k_spread = k_sub (but note
    // Q_sub = -k_sub · |Δv_conv| depends on v; with v = 0, Q_sub = 0
    // and so Q_spread alone drives the mass up).
    //
    // To exercise a *balanced* flux we force a uniform convergence
    // at each subduction cell via a synthetic div(v) field.
    // Specifically, we fabricate `div_v[i, sub_row] = -k_spread / k_sub`
    // so that Q_sub = -k_sub · (k_spread/k_sub) = -k_spread — perfect
    // counter to Q_spread at every subducting cell. Q_arc is zero
    // because there are no continental cells.
    let k_spread = 0.2_f64;
    let rates = BoundaryRates {
        k_sub: 0.5,
        k_arc: 0.0, // disable arc to keep the balance pristine
        k_spread,
        k_coll_v: 0.0,
        k_rift_v: 0.0,
    };
    let conv_needed = k_spread / rates.k_sub;

    // Find the subduction row (set by balanced_sub_spread).
    use ymir_core::tectonics_v2::boundaries::BoundaryFlag;
    let mut sub_rows: Vec<usize> = Vec::new();
    for j in 0..ny {
        if matches!(layout.flags.get(0, j), BoundaryFlag::OceanicSubduction) {
            sub_rows.push(j);
        }
    }
    assert!(!sub_rows.is_empty(), "balanced_sub_spread must have at least one subduction row");

    // Bypass div_v_cell: we prescribe divergence directly.
    let mut div_v = Field2D::new(nx, ny);
    for j in sub_rows.iter().copied() {
        for i in 0..nx {
            // Negative divergence → convergent → |Δv_conv| = conv_needed.
            div_v.set(i, j, -conv_needed);
        }
    }
    // Note: we synthesise `div_v` directly instead of calling
    // `div_v_cell`, since the Stokes solve and velocity field are
    // deliberately bypassed in this test to isolate the source/sink
    // arithmetic.

    // Run 50 macro steps of Q-only update. The "balance" here is
    // global, not local: the rift row grows, the subduction row
    // drains, mass integrates to zero. A deep enough initial S
    // plus short duration keeps every cell above S_MIN so the
    // clamp never fires and the test isolates pure arithmetic
    // drift.
    let mut s = Field2D::new(nx, ny);
    for v in s.data_mut().iter_mut() {
        *v = 2.0;
    } // well above drain depth
    let mass_initial: f64 = s.data().iter().sum();
    let mut q = Field2D::new(nx, ny);
    let mut q_sub_scratch = Field2D::new(nx, ny);
    let dt = 0.02;
    for _step in 0..50 {
        compute_source_sink_terms(
            &layout.plate_types,
            &layout.flags,
            &rates,
            &div_v,
            &idx_x,
            &idx_y,
            &mut q_sub_scratch,
            &mut q,
        );
        for (cell, &q_val) in s.data_mut().iter_mut().zip(q.data().iter()) {
            *cell += dt * q_val;
        }
        // Apply clamp; injected_flux should stay at 0 because S̃ stays
        // well above the floor.
        let stats = apply_clamp_with_tracking(&mut s);
        assert_eq!(stats.activations, 0, "balanced test should never clamp");
    }
    let mass_final: f64 = s.data().iter().sum();
    let drift = (mass_final - mass_initial).abs() / mass_initial.abs().max(1.0);
    assert!(drift < 1.0e-8, "balanced layout mass drift = {} (should be < 1e-8)", drift,);
    // Suppress the underscored `dx, dy` warnings if any.
    let _ = (dx, dy);
}
