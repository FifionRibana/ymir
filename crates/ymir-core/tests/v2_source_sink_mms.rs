//! MMS on the Step-5 source/sink terms, with zero velocity.
//!
//! Isolates the arithmetic of `compute_source_sink_terms` from
//! advection and from the Stokes solve. On a uniform `S̃ = 1` field
//! with `ṽ = 0`:
//!
//! - `Q_sub` is zero everywhere (no convergent motion).
//! - `Q_spread, Q_coll-v, Q_rift-v` evaluate to their constant rates
//!   on flagged cells.
//! - `Q_arc` is zero because `Q_sub` is zero.
//!
//! After one macro step of `S̃_next = S̃ + Δt·Q`, each cell's change
//! must equal `Δt · Q` to machine precision, and cells flagged `None`
//! must be strictly unchanged.
//!
//! Layout: `continental_collision_band` (all-continental; one
//! collision row). That exercises `Q_coll-v` and keeps the test
//! pointwise-checkable: every cell on the collision row must change
//! by exactly `Δt · k_coll-v`; every other cell must stay at 1.0.

use ymir_core::tectonics_v2::boundaries::{
    BoundaryFlag, BoundaryRates, PlateType, compute_source_sink_terms, continental_collision_band,
};
use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};

#[test]
fn source_sink_on_static_continental_collision_band_exact_cell_delta() {
    let nx = 16;
    let ny = 12;
    let layout = continental_collision_band(nx, ny);
    let rates = BoundaryRates::baseline_uncalibrated();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    // Zero velocity → zero divergence → zero |Δv_conv| everywhere,
    // so Q_sub branch is inactive even if flags were present.
    let div_v = Field2D::new(nx, ny);
    let mut q = Field2D::new(nx, ny);
    let mut q_sub_scratch = Field2D::new(nx, ny);
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

    // One macro step of S += dt·Q, no advection (v = 0).
    let dt = 0.02_f64;
    let mut s = Field2D::new(nx, ny);
    for v in s.data_mut().iter_mut() {
        *v = 1.0;
    }
    let before: Vec<f64> = s.data().to_vec();
    for (cell, &q_val) in s.data_mut().iter_mut().zip(q.data().iter()) {
        *cell += dt * q_val;
    }

    let coll_j = ny / 2;
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            let expected = match (layout.plate_types.get(i, j), layout.flags.get(i, j)) {
                (PlateType::Continental, BoundaryFlag::ContinentalCollision) => {
                    1.0 + dt * rates.k_coll_v
                }
                _ => 1.0,
            };
            let got = s.data()[idx];
            assert!(
                (got - expected).abs() < 1e-14,
                "cell ({},{}): expected {}, got {}",
                i,
                j,
                expected,
                got,
            );
            // Sanity: non-collision cells are bit-identical to
            // before. This is what "strictly unchanged" means —
            // even a rounding-error drift would violate the Lie
            // splitting assumption that `None` cells do not carry Q.
            if j != coll_j {
                assert_eq!(got, before[idx], "None cell should be unchanged bit-exactly");
            }
        }
    }
}

#[test]
fn source_sink_all_none_layout_has_zero_q() {
    // A layout with `BoundaryFlag::None` everywhere and uniform
    // continental plates: every Q term is zero.
    let nx = 8;
    let ny = 6;
    use ymir_core::tectonics_v2::boundaries::{BoundaryFlagField, PlateTypeField};
    let plate_types = PlateTypeField::filled(nx, ny, PlateType::Continental);
    let flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    let rates = BoundaryRates::baseline_uncalibrated();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let div_v = Field2D::new(nx, ny);
    let mut q = Field2D::new(nx, ny);
    let mut q_sub_scratch = Field2D::new(nx, ny);
    compute_source_sink_terms(
        &plate_types,
        &flags,
        &rates,
        &div_v,
        &idx_x,
        &idx_y,
        &mut q_sub_scratch,
        &mut q,
    );
    for &v in q.data() {
        assert_eq!(v, 0.0);
    }
}
