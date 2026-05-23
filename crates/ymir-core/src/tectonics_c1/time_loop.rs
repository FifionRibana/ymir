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
//! ## Phase 1.2 / 1.3 contract — [`run_with_closures`]
//!
//! Adds per-step closure source / sink terms after each advection
//! step. Per-step structure:
//!
//! 1. CFL Δt.
//! 2. Advect `S̃` and `age` (same as Phase 1.1).
//! 3. Apply Davis-Suppe orogenic source term on upper-plate
//!    interior cells via
//!    [`super::closures::davis_suppe::source_term::apply_davis_suppe_step`]
//!    (Phase 1.2).
//! 4. Apply equilibrium-height gravitational sink globally via
//!    [`super::closures::equilibrium_height::source_term::apply_equilibrium_height_step`]
//!    (Phase 1.3). Strict ordering: AFTER Davis-Suppe — reversing
//!    would oscillate around `h_eq` instead of converging.
//! 5. Diagnostic callback.
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

use crate::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use crate::tectonics_v2::advection::step_upwind;
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::workflow::drainage::compute_drainage_targets;
use crate::tectonics_v2::workflow::phase_a_common::compute_sea_level_ref_s_space;

use super::boundary_classification::classify_boundaries;
use super::closures::davis_suppe::source_term::{apply_davis_suppe_step, DavisSuppeParams};
use super::closures::equilibrium_height::params::EquilibriumHeightParams;
use super::closures::equilibrium_height::source_term::apply_equilibrium_height_step;
use super::closures::erosion::params::ErosionParams;
use super::closures::erosion::source_term::{apply_erosion_step, compute_drainage_areas};
use super::distance_field::wedge_distance_intra_plate;
use super::kinematics::PlateKinematics;
use super::state::C1State;

/// Tunables for [`run_advection_only`] and [`run_with_closures`].
///
/// `dx` and `dy` are the cell-size in non-dimensional length units.
/// For the typical unit-domain `1×1` non-dim setup they both equal
/// `1.0 / grid_size`.
///
/// `iso_config` and `drainage_max_distance` are consumed by
/// [`run_with_closures`]'s per-step isostasy + drainage +
/// stream-power-erosion path (Phase 1.4). They are unused by
/// [`run_advection_only`] but must be present in the config for
/// the shared struct surface.
#[derive(Clone, Debug)]
pub struct C1TimeLoopConfig {
    pub n_steps: usize,
    pub dx: f64,
    pub dy: f64,
    /// Isostasy parameters consumed by the per-step erosion path
    /// (Phase 1.4). Used by [`compute_isostasy`] for the altitude
    /// heightmap and by
    /// [`compute_sea_level_ref_s_space`] for the
    /// drainage-classification sea-level threshold.
    pub iso_config: IsostasyConfig,
    /// Maximum drainage path length (cells) for
    /// [`compute_drainage_targets`]. Default `30` mirrors the
    /// Phase 1.2 + 1.3 default for wedge / drainage distances.
    pub drainage_max_distance: usize,
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
/// flag (W4 isolation discipline). The bundle grows alongside the
/// C1 milestone:
///
/// - Phase 1.2 (Issue #123) — `davis_suppe` (orogenic source).
/// - Phase 1.3 (Issue #125) — `equilibrium_height` (gravitational
///   collapse sink, Molnar-Lyon-Caen).
/// - Phase 1.4 (Issue #127) — `erosion` (stream-power incision
///   sink, Whipple-Tucker 1999 / Lague 2014).
/// - Phase 2 — `parsons_sclater` (oceanic bathymetry). TBA.
///
/// ## Default-behaviour caveat
///
/// `C1Closures::default()` enables **all three** closures
/// (Davis-Suppe + equilibrium-height + erosion). Tests written
/// for a prior-phase regime (where a closure-specific observable
/// is load-bearing — e.g. the Phase 1.2 unbounded boundary pile-up
/// `global_max ≈ 2297`, or the Phase 1.3 `wedge_p95 = 0.376`
/// preservation) must explicitly disable the later-phase
/// closures:
///
/// ```ignore
/// let closures = C1Closures {
///     davis_suppe: DavisSuppeParams::default(),
///     equilibrium_height: EquilibriumHeightParams {
///         enabled: false,
///         ..EquilibriumHeightParams::default()
///     },
///     erosion: ErosionParams {
///         enabled: false,
///         ..ErosionParams::default()
///     },
/// };
/// ```
#[derive(Clone, Copy, Debug)]
pub struct C1Closures {
    pub davis_suppe: DavisSuppeParams,
    pub equilibrium_height: EquilibriumHeightParams,
    pub erosion: ErosionParams,
}

impl Default for C1Closures {
    fn default() -> Self {
        Self {
            davis_suppe: DavisSuppeParams::default(),
            equilibrium_height: EquilibriumHeightParams::default(),
            erosion: ErosionParams::default(),
        }
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

        // 2. Davis-Suppe orogenic source term (Phase 1.2).
        apply_davis_suppe_step(
            &mut state.s,
            &state.plate_id,
            &boundary,
            &wedge_d,
            kinematics,
            &closures.davis_suppe,
            dt,
        );

        // 3. Equilibrium height sink (Phase 1.3).
        // Order critical: AFTER Davis-Suppe.
        // Reverse order causes oscillation:
        //   - Equilibrium caps at h_eq
        //   - Davis-Suppe re-injects mass above h_eq
        //   - Next step: excess from re-injection
        //   - Etc., oscillation around h_eq instead of stable equilibrium.
        apply_equilibrium_height_step(&mut state.s, &closures.equilibrium_height, dt);

        // 4. Stream-power erosion sink (Phase 1.4 — Issue #127).
        //
        // Three-stage pipeline (W1 strict order):
        //   (a) Isostasy — `S̃ → altitude` heightmap via Airy.
        //       Needed for the slope magnitude in W-T eq. (1).
        //   (b) Drainage — classify each cell's drainage target
        //       via `compute_drainage_targets` (operates on the
        //       cellular `S̃` with the Phase 3.5 S̃-space sea-
        //       level threshold), then convert per-cell targets
        //       into transitive drainage *areas* via
        //       `compute_drainage_areas`.
        //   (c) Erosion — apply W-T `E = K · A^m · S^n` with the
        //       safety floor at the oceanic baseline.
        //
        // Order critical: AFTER equilibrium-height. The erosion
        // closure reads altitude *after* the height cap has been
        // applied; reversing would produce visible incision on
        // boundary pile-up cells that the equilibrium clamp is
        // about to remove anyway — wasted computation, and the
        // cap-then-erode order matches the W-T citing of
        // Molnar-Lyon-Caen for `h_effective ≈ min(h_collapse,
        // h_erosion)`.
        //
        // Skipped entirely when `closures.erosion.enabled` is
        // `false` (W4 closure-isolation discipline). The
        // `apply_erosion_step` early-return guards the per-step
        // overhead so the Phase 1.2 / Phase 1.3 regression tests
        // remain bit-identical when erosion is disabled.
        if closures.erosion.enabled {
            let isostasy = compute_isostasy(&state.s, &config.iso_config);
            let sea_level_ref = compute_sea_level_ref_s_space(&state.s, &config.iso_config);
            let drainage_map = compute_drainage_targets(
                &state.s,
                sea_level_ref,
                config.drainage_max_distance,
            );
            let drainage_areas = compute_drainage_areas(&drainage_map);
            apply_erosion_step(
                &mut state.s,
                &isostasy.heightmap,
                &drainage_areas,
                &closures.erosion,
                dt,
                config.dx,
            );
        }

        // 5. (Future closures land here — Phase 2 oceanic bathymetry.)

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
        let config = C1TimeLoopConfig {
            n_steps: 100,
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            iso_config: IsostasyConfig::default(),
            drainage_max_distance: 30,
        };

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
        let config = C1TimeLoopConfig {
            n_steps: 7,
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            iso_config: IsostasyConfig::default(),
            drainage_max_distance: 30,
        };

        let mut steps_seen = Vec::new();
        run_advection_only(&mut state, &kinematics, &config, |s, _| {
            steps_seen.push(s);
        });
        assert_eq!(steps_seen, vec![0, 1, 2, 3, 4, 5, 6]);
    }
}
