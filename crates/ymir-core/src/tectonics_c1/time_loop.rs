//! C1 time loop — forward-Euler advection of `S̃` and `age`.
//!
//! ## Phase 1.1 contract
//!
//! Advection only, no closures. At each step:
//!
//! 1. Compute `Δt = 0.5 · min(dx, dy) / max(|v|, 1e-12)` (CFL,
//!    safety factor 0.5 from §4.6 of the design doc).
//! 2. Build per-cell velocity field `(vx, vy)` from per-plate
//!    kinematics — every cell of plate `p` gets `velocities[p]`
//!    (cratonic and non-cratonic cells alike at this stage).
//! 3. Advect `S̃` and `age` independently using v2's
//!    [`crate::tectonics_v2::advection::step_upwind`] scheme.
//! 4. Invoke the caller-supplied `on_step` callback for
//!    diagnostics (snapshot dumping, mass-drift checks, …).
//!
//! Closures join in Phase 1.2+ between steps 3 and 4.
//!
//! ## Reuse, do not reinvent
//!
//! The upwind scheme is `step_upwind` from
//! `crate::tectonics_v2::advection`. No reimplementation here
//! (W1 watchpoint of Issue #120).

use crate::tectonics_v2::advection::step_upwind;
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

use super::kinematics::PlateKinematics;
use super::state::C1State;

/// Tunables for [`run_advection_only`].
///
/// `dx` and `dy` are the cell-size in non-dimensional length units.
/// For the typical unit-domain `1×1` non-dim setup they both equal
/// `1.0 / grid_size`.
#[derive(Clone, Copy, Debug)]
pub struct C1TimeLoopConfig {
    pub n_steps: usize,
    pub dx: f64,
    pub dy: f64,
}

/// Run advection-only forward in time. The callback fires once
/// per step **after** the state has been updated, with the
/// 0-based step index already incremented (so first invocation
/// is `on_step(0, …)` representing "post-step-1 state").
pub fn run_advection_only<F>(
    state: &mut C1State,
    kinematics: &PlateKinematics,
    config: &C1TimeLoopConfig,
    mut on_step: F,
) where
    F: FnMut(usize, &C1State),
{
    let nx = state.nx();
    let ny = state.ny();
    let n_cells = nx * ny;

    // CFL Δt. `max_velocity` returns 0 only for empty kinematics;
    // the floor guards against div-by-zero. Per-plate magnitudes
    // are constant in Phase 1.1 so Δt is constant across the run.
    let max_v = kinematics.max_velocity().max(1e-12);
    let dx_min = config.dx.min(config.dy);
    let dt = 0.5 * dx_min / max_v;

    // Per-step velocity field is identical across steps in
    // Phase 1.1 (constant kinematics, static plate_id). Build
    // once and reuse — saves N_steps × N_cells worth of churn
    // (per W5: avoid obvious pessimisation).
    let mut vx = vec![0.0_f64; n_cells];
    let mut vy = vec![0.0_f64; n_cells];
    fill_velocity_field(&mut vx, &mut vy, state, kinematics);

    // Pre-allocated scratch buffers for the upwind sweep.
    let mut s_next = Field2D::new(nx, ny);
    let mut age_next = Field2D::new(nx, ny);
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    for step in 0..config.n_steps {
        // Advect S̃.
        step_upwind(
            nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.s, &vx, &vy, &mut s_next,
        );
        // Advect age.
        step_upwind(
            nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.age, &vx, &vy, &mut age_next,
        );

        // Swap rather than reassign — keeps the same backing
        // allocations alive for the next iteration.
        std::mem::swap(&mut state.s, &mut s_next);
        std::mem::swap(&mut state.age, &mut age_next);

        on_step(step, state);
    }
}

/// Populate the per-cell velocity slices from the plate-id map
/// and the per-plate kinematics.
fn fill_velocity_field(
    vx: &mut [f64],
    vy: &mut [f64],
    state: &C1State,
    kinematics: &PlateKinematics,
) {
    let nx = state.nx();
    let ny = state.ny();
    for j in 0..ny {
        for i in 0..nx {
            let plate = state.plate_id.get(i, j) as usize;
            let (vx_p, vy_p) = kinematics.velocities[plate];
            let k = j * nx + i;
            vx[k] = vx_p;
            vy[k] = vy_p;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
    use crate::tectonics_v2::voronoi::PlateIdField;

    use super::super::state::{BoolField, C1State};

    /// Synthetic single-plate state with uniform S̃ — advection
    /// under uniform velocity is a pure translation, total mass
    /// (sum of `S̃` over all cells) is conserved exactly under
    /// periodic boundaries.
    fn uniform_single_plate_state(nx: usize, ny: usize) -> C1State {
        let mut s = Field2D::new(nx, ny);
        let mut age = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                s.set(i, j, 1.0);
                age.set(i, j, 0.0);
            }
        }
        let plate_id = PlateIdField::new(nx, ny); // all zeros
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let cratonic_mask = BoolField::filled(nx, ny, false);
        C1State { s, age, plate_id, plate_type, cratonic_mask, num_plates: 1 }
    }

    #[test]
    fn uniform_field_advection_conserves_mass_exactly() {
        let nx = 16;
        let ny = 16;
        let mut state = uniform_single_plate_state(nx, ny);
        let initial_mass: f64 = state.s.data().iter().sum();

        let kinematics = PlateKinematics { velocities: vec![(0.01, 0.005)] };
        let config = C1TimeLoopConfig { n_steps: 100, dx: 1.0 / nx as f64, dy: 1.0 / ny as f64 };

        run_advection_only(&mut state, &kinematics, &config, |_, _| {});

        let final_mass: f64 = state.s.data().iter().sum();
        let drift = (final_mass - initial_mass).abs();
        assert!(
            drift < 1e-12,
            "uniform-field advection should conserve mass exactly under periodic BCs; drift = {:.3e}",
            drift
        );
    }

    #[test]
    fn callback_fires_once_per_step() {
        let nx = 8;
        let ny = 8;
        let mut state = uniform_single_plate_state(nx, ny);
        let kinematics = PlateKinematics { velocities: vec![(0.01, 0.0)] };
        let config = C1TimeLoopConfig { n_steps: 7, dx: 1.0 / nx as f64, dy: 1.0 / ny as f64 };

        let mut steps_seen = Vec::new();
        run_advection_only(&mut state, &kinematics, &config, |s, _| {
            steps_seen.push(s);
        });
        assert_eq!(steps_seen, vec![0, 1, 2, 3, 4, 5, 6]);
    }
}
