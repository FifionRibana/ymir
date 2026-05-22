//! C1 time loop — forward-Euler advection of `S̃` and `age`
//! plus per-step closure source terms (Phase 1.2+).
//!
//! ## Phase 1.1 contract — [`run_advection_only`]
//!
//! Advection only, no closures. Preserved verbatim for the W4
//! closure-disabled regression: a Phase 1.2 run with all
//! closures disabled is **not** equivalent to a `run_with_closures`
//! call with `params.enabled = false` (the latter still computes
//! boundary + wedge distance once at setup — small overhead, but
//! not bit-identical). The clean closure-OFF baseline is to call
//! [`run_advection_only`] directly.
//!
//! ## Phase 1.2 contract — [`run_with_closures`]
//!
//! Adds the Davis-Suppe orogenic source term after each advection
//! step. Per-step structure:
//!
//! 1. CFL Δt.
//! 2. Advect `S̃` and `age` (same as Phase 1.1).
//! 3. Apply Davis-Suppe source term on upper-plate interior cells
//!    via [`super::closures::davis_suppe::source_term::apply_davis_suppe_step`].
//! 4. Diagnostic callback.
//!
//! ## Static-classification optimisation
//!
//! In Phase 1.2 the `plate_id` field does **not** advect — only
//! `S̃` and `age` move under advection. As a consequence the
//! Stage 2 boundary classification and the Stage 3.1 intra-plate
//! wedge distance are **static throughout the run** and are
//! computed once outside the loop. Phase 2 (boundary evolution)
//! will lift them back inside the loop when the underlying
//! geometry actually changes per step.
//!
//! ## Reuse, do not reinvent
//!
//! The upwind scheme is `step_upwind` from
//! `crate::tectonics_v2::advection`. No reimplementation here
//! (W1 watchpoint of Issue #120).

use crate::tectonics_v2::advection::step_upwind;
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

use super::boundary_classification::classify_boundaries;
use super::closures::davis_suppe::source_term::{apply_davis_suppe_step, DavisSuppeParams};
use super::distance_field::wedge_distance_intra_plate;
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

/// Bundle of all C1 closures and their parameters.
///
/// Each closure is independently togglable via its own `enabled`
/// flag. The bundle grows alongside the C1 milestone:
///
/// - Phase 1.2 (Issue #123) — `davis_suppe` (this issue).
/// - Phase 1.3 — `equilibrium_height` (Molnar-Lyon-Caen). TBA.
/// - Phase 1.4 — `macro_erosion` (Whipple-Tucker) and isostasy
///   hook. TBA.
/// - Phase 2 — `parsons_sclater` (oceanic bathymetry). TBA.
#[derive(Clone, Copy, Debug)]
pub struct C1Closures {
    pub davis_suppe: DavisSuppeParams,
}

impl Default for C1Closures {
    fn default() -> Self {
        Self { davis_suppe: DavisSuppeParams::default() }
    }
}

/// Run the C1 forward-Euler time loop with per-step closure source
/// terms applied after each advection update.
///
/// Phase 1.2 wires only the Davis-Suppe orogenic closure; future
/// phases extend [`C1Closures`] and add new `apply_*_step` calls
/// in the per-step body of this function.
///
/// ## Static-classification optimisation (Phase 1.2)
///
/// Since `plate_id` is not advected, the boundary classification
/// and the intra-plate wedge-distance field are static throughout
/// the run. They are computed **once before** the loop and reused
/// every step. Phase 2 (boundary evolution) will move them back
/// inside the loop.
pub fn run_with_closures<F>(
    state: &mut C1State,
    kinematics: &PlateKinematics,
    config: &C1TimeLoopConfig,
    closures: &C1Closures,
    mut on_step: F,
) where
    F: FnMut(usize, &C1State),
{
    let nx = state.nx();
    let ny = state.ny();
    let n_cells = nx * ny;

    let max_v = kinematics.max_velocity().max(1e-12);
    let dx_min = config.dx.min(config.dy);
    let dt = 0.5 * dx_min / max_v;

    let mut vx = vec![0.0_f64; n_cells];
    let mut vy = vec![0.0_f64; n_cells];
    fill_velocity_field(&mut vx, &mut vy, state, kinematics);

    // Static-classification optimisation: plate_id is invariant
    // in Phase 1.2, so the boundary verdict and the intra-plate
    // wedge distance are computed once here, reused every step.
    let boundary = classify_boundaries(&state.plate_id, kinematics);
    let wedge_d = wedge_distance_intra_plate(
        &state.plate_id,
        &boundary.upper_plate_mask,
        closures.davis_suppe.max_distance,
    );

    let mut s_next = Field2D::new(nx, ny);
    let mut age_next = Field2D::new(nx, ny);
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    for step in 0..config.n_steps {
        // 1. Advection (Phase 1.1 unchanged).
        step_upwind(
            nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.s, &vx, &vy, &mut s_next,
        );
        step_upwind(
            nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.age, &vx, &vy,
            &mut age_next,
        );
        std::mem::swap(&mut state.s, &mut s_next);
        std::mem::swap(&mut state.age, &mut age_next);

        // 2. Davis-Suppe orogenic source term.
        apply_davis_suppe_step(
            &mut state.s,
            &state.plate_id,
            &boundary,
            &wedge_d,
            kinematics,
            &closures.davis_suppe,
            dt,
        );

        // 3. (Future closures land here.)

        on_step(step, state);
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
