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
//! ## Phase 1.2 / 1.3 / 1.4 contract — [`run_with_closures`]
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
//! 5. Apply stream-power erosion sink via the per-step isostasy
//!    + drainage-targets + drainage-areas + erosion pipeline
//!    (Phase 1.4 — see [`run_with_closures`] for the per-stage
//!    breakdown). Skipped entirely when `closures.erosion.enabled`
//!    is `false`.
//! 6. Diagnostic callback.
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
use crate::tectonics_v2::advection::{step_upwind, step_upwind_masked};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::workflow::drainage::compute_drainage_targets;
use crate::tectonics_v2::workflow::phase_a_common::compute_sea_level_ref_s_space;

use super::boundary_classification::classify_boundaries;
use super::closures::accretion::{apply_accretion_step, AccretionParams, ConvergenceTracker};
use super::closures::davis_suppe::source_term::{apply_davis_suppe_step, DavisSuppeParams};
use super::closures::equilibrium_height::params::EquilibriumHeightParams;
use super::closures::equilibrium_height::source_term::apply_equilibrium_height_step;
use super::closures::erosion::params::ErosionParams;
use super::closures::erosion::source_term::{apply_erosion_step, compute_drainage_areas};
use super::closures::oceanic_bathymetry::params::SteinSteinParams;
use super::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use super::closures::rifting::{
    apply_rifting_split, apply_rifting_thinning, DivergenceTracker, RiftingParams,
};
use super::closures::subduction::{apply_subduction_step, SubductionParams};
use super::distance_field::wedge_distance_intra_plate;
use super::kinematics::PlateKinematics;
use super::state::C1State;
use crate::tectonics_v2::boundaries::plate_type::PlateType;

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
    /// **TRANSITIONAL fix flag (Issue #145) — NOT a permanent mode.**
    ///
    /// `false` = the legacy advection (continental crust advected
    /// uniformly like a fluid). This is **physically incorrect** — it
    /// destroys continents (a convergent velocity field sweeps all
    /// crust into pile-ups and empties the interiors; cratons 1.00→0
    /// over a run, proven in the piste-4 investigation). It is kept
    /// ONLY transitionally, as the byte-identical A/B reference while
    /// the Phase 1.x/2 closures are re-validated on the corrected
    /// transport.
    ///
    /// `true` = the fix: continental crust (`plate_type==Continental`)
    /// is **rigid / non-subducting** — its advection velocity is
    /// zeroed (upstream v-masking in `fill_velocity_field`; still
    /// participates in face flux, so mass is conserved). Preserves
    /// ~92 % of cratonic area on the diagnostic.
    ///
    /// **EXIT PLAN:** once the full Phase 1.x/2 re-validation on rigid
    /// transport is complete, the default flips to `true` (or the flag
    /// is removed and rigidity becomes unconditional). `false` must
    /// NOT fossilise into a "legitimate mode" — it is the old broken
    /// behaviour, retained only for the re-validation transition. See
    /// `docs/reports/c1_continental_buoyancy/`.
    pub rigid_continental_crust: bool,
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
    let rigid_mask: Option<Vec<bool>> = if config.rigid_continental_crust {
        apply_continental_rigidity(&mut vx, &mut vy, state);
        Some(continental_rigid_mask(state))
    } else {
        None
    };

    // Pre-allocated scratch buffers for the upwind sweep.
    let mut s_next = Field2D::new(nx, ny);
    let mut age_next = Field2D::new(nx, ny);
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    for step in 0..config.n_steps {
        // Advect S̃ (no-flux rigid boundary when rigid; #145).
        match &rigid_mask {
            Some(m) => step_upwind_masked(
                nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.s, &vx, &vy, m, &mut s_next,
            ),
            None => step_upwind(
                nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.s, &vx, &vy, &mut s_next,
            ),
        }
        // Advect age (same no-flux rigid boundary; #145).
        match &rigid_mask {
            Some(m) => step_upwind_masked(
                nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.age, &vx, &vy, m,
                &mut age_next,
            ),
            None => step_upwind(
                nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.age, &vx, &vy,
                &mut age_next,
            ),
        }

        // Swap rather than reassign — keeps the same backing
        // allocations alive for the next iteration.
        std::mem::swap(&mut state.s, &mut s_next);
        std::mem::swap(&mut state.age, &mut age_next);

        on_step(step, state);
    }
}

/// Issue #145 — continental rigidity (the buoyancy fix). Zero the
/// advection velocity on `plate_type==Continental` cells so buoyant
/// continental crust does NOT advect (does not subduct / funnel into
/// convergence). Applied AFTER [`fill_velocity_field`], at the same
/// hook. Conservative: rigid cells still participate in `step_upwind`
/// face flux (mass cancels), so this is not a freeze and does not
/// inject/destroy mass. No-op unless `rigid_continental_crust`.
fn apply_continental_rigidity(vx: &mut [f64], vy: &mut [f64], state: &C1State) {
    for (k, t) in state.plate_type.data().iter().enumerate() {
        if matches!(t, PlateType::Continental) {
            vx[k] = 0.0;
            vy[k] = 0.0;
        }
    }
}

/// Issue #145 — rigid mask (`plate_type==Continental`) for the no-flux boundary
/// in [`step_upwind_masked`]. Rebuilt wherever the velocity field is rebuilt
/// (once in gallery, per-step on the Track-D path).
fn continental_rigid_mask(state: &C1State) -> Vec<bool> {
    state
        .plate_type
        .data()
        .iter()
        .map(|t| matches!(t, PlateType::Continental))
        .collect()
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
/// - Phase 2 Track A (Issue #129) — `oceanic_bathymetry`
///   (Stein-Stein 1992 post-isostasy depth assignment,
///   **Architecture C** — see [`SteinSteinParams`] and
///   [`super::closures::oceanic_bathymetry`] module docstring).
///
/// ## Default-behaviour caveat
///
/// `C1Closures::default()` enables **all seven** closures
/// (Davis-Suppe + equilibrium-height + erosion + oceanic
/// bathymetry + subduction + accretion + rifting). Tests written
/// for a prior-phase regime (where a closure-specific observable
/// is load-bearing — e.g. the Phase 1.2 unbounded boundary
/// pile-up `global_max ≈ 2297`, the Phase 1.3 `wedge_p95 = 0.376`
/// preservation, or the Track B 8th bit-identical decomposition
/// contract) must explicitly disable the later-phase closures.
/// Phase 1.x test fixtures use the `phase_1_X_closures()` helper
/// pattern (see [`crates/ymir-core/tests/c1_phase_1_4_erosion.rs`])
/// to disable forward-only closures with `enabled: false` overrides.
///
/// **Track D defaults** (subduction, accretion, rifting): each
/// `Default::default()` ships with `enabled: true` per its own
/// docstring, so `C1Closures::default()` includes them. Phase 1.x
/// + Track A + Track B tests must add the trio of `enabled: false`
/// overrides to preserve the regression invariants. Track D
/// integration tests + Phase 2 milestone gate runs leave the
/// defaults intact.
///
/// ```ignore
/// let closures = C1Closures {
///     subduction: SubductionParams {
///         enabled: false,
///         ..SubductionParams::default()
///     },
///     accretion: AccretionParams {
///         enabled: false,
///         ..AccretionParams::default()
///     },
///     rifting: RiftingParams {
///         enabled: false,
///         ..RiftingParams::default()
///     },
///     ..C1Closures::default()
/// };
/// ```
#[derive(Clone, Copy, Debug)]
pub struct C1Closures {
    pub davis_suppe: DavisSuppeParams,
    pub equilibrium_height: EquilibriumHeightParams,
    pub erosion: ErosionParams,
    pub oceanic_bathymetry: SteinSteinParams,
    /// Track D — subduction (Issue #132 Stage E1).
    pub subduction: SubductionParams,
    /// Track D — accretion (Issue #132 Stage E2).
    pub accretion: AccretionParams,
    /// Track D — rifting (Issue #132 Stage E3).
    pub rifting: RiftingParams,
}

impl Default for C1Closures {
    fn default() -> Self {
        Self {
            davis_suppe: DavisSuppeParams::default(),
            equilibrium_height: EquilibriumHeightParams::default(),
            erosion: ErosionParams::default(),
            oceanic_bathymetry: SteinSteinParams::default(),
            subduction: SubductionParams::default(),
            accretion: AccretionParams::default(),
            rifting: RiftingParams::default(),
        }
    }
}

/// Run the C1 forward-Euler time loop with per-step closure
/// source / sink terms applied after each advection update.
///
/// ## Per-step pipeline
///
/// Each iteration of the `0..config.n_steps` loop executes the
/// following stages, in strict order:
///
/// | # | Stage | Phase | Complexity per step |
/// |---|-------|-------|---------------------|
/// | 1 | Advection of `S̃` + `age` via `step_upwind` | 1.1 | `O(N)` |
/// | 2 | Davis-Suppe orogenic source — `apply_davis_suppe_step` | 1.2 | `O(N)` |
/// | 3 | Equilibrium-height sink — `apply_equilibrium_height_step` | 1.3 | `O(N)` |
/// | 4 | Stream-power erosion (4 sub-steps, gated by `closures.erosion.enabled`): | 1.4 | `O(N · max_d + N log N)` |
/// |   | 4a — `compute_isostasy` → altitude heightmap | | `O(N)` + Gaussian blur |
/// |   | 4b — `compute_sea_level_ref_s_space` → S̃-space threshold | | `O(N)` |
/// |   | 4c — `compute_drainage_targets` → `DrainageMap` | | `O(N · max_distance)` BFS |
/// |   | 4d — `compute_drainage_areas` → `Vec<u32>` transitive areas | | `O(N log N)` sort + iter |
/// |   | 4e — `apply_erosion_step` (W-T `K · A^m · S^n`) | | `O(N)` |
/// | 5 | Diagnostic `on_step(step, &state)` callback | — | caller |
///
/// `N = nx · ny` is the cell count. `max_distance` is bounded by
/// `config.drainage_max_distance` (default `30`). Stage 4 dominates
/// the per-step cost (~ 320 µs at 64² × Phase-1.1 kinematics vs
/// ~ 50 µs for stages 1-3 combined — see § Performance below).
///
/// ## Why per-step isostasy
///
/// Stages 4a → 4d are coupled: erosion (4e) needs **altitude**
/// (from 4a) AND **drainage areas** (from 4d), where drainage
/// targets (4c) classify oceanic vs continental cells using the
/// **current** `S̃` distribution via `sea_level_ref` (4b). Running
/// 4a once at start-of-run would compute drainage on stale
/// altitude — after a few hundred steps of Davis-Suppe source +
/// equilibrium clamp, the altitude field has shifted enough that
/// the drainage classification would be wrong. Per-step
/// recomputation is the simplest defensible choice.
///
/// Cost is bounded — § Performance below shows the full pipeline
/// at 64² takes ~ 110 ms for 300 steps. The Gaussian blur inside
/// `compute_isostasy` (default σ = 2.0) is the single largest
/// contributor (~ 150 µs/step). If profiling identifies it as a
/// bottleneck at 512², the blur can be moved to `apply_post_
/// tectonic` (end-of-cycle only) once the C1 cycle pattern
/// requires per-cycle altitude smoothing rather than per-step.
///
/// ## End-of-cycle `apply_post_tectonic` consistency
///
/// The C1 workflow wrapper
/// [`crate::tectonics_v2::workflow::phase_a_c1::run_phase_a_cycle_c1`]
/// runs this loop, then invokes
/// [`crate::tectonics_v2::workflow::phase_a_common::apply_post_tectonic`]
/// at the end of the cycle. The post-tectonic pass re-runs:
///
/// - **Sea-level**: same Phase 3.5 formula via the helper
///   `compute_sea_level_ref_s_space` extracted in Stage E0.
///   Per-step (4b) and end-of-cycle reuse the same code.
/// - **Macro-redistribution**: in-place mass redistribution under
///   drainage targets. Not run per step (only end-of-cycle).
/// - **Reclassification + cratonic recompute**: end-of-cycle only.
///
/// This means the per-step pipeline (stages 4a-4e) duplicates the
/// isostasy + sea-level computation that `apply_post_tectonic`
/// will redo. The redundancy is **mildly costly but
/// architecturally cleaner**: the time loop only knows about
/// per-step needs (erosion), the workflow wrapper only knows about
/// per-cycle needs (macro mass + reclass + craton). Sharing
/// computed altitude across the boundary would require a
/// `&mut FinalState`-like envelope that C1 deliberately avoids
/// (see `phase_a_c1` module docstring on the asymmetric API).
///
/// ## Performance
///
/// At 64² × 300 steps, Phase 1.1 kinematics, default closures
/// (DS + EH + erosion all enabled):
///
/// - Phase 1.3 baseline (DS + EH, no erosion): **29 ms / 300 steps
///   = 96 µs/step**.
/// - Phase 1.4 measurement (this commit chain, all 3 closures):
///   **~110 ms / 300 steps = ~367 µs/step**.
/// - ~ 3.8× slowdown vs Phase 1.3 baseline; well within the
///   user-spec acceptable range and well below the design-doc
///   §2.3 < 10 s / 512² target.
///
/// Per-step breakdown estimate at 64²:
///
/// | Stage | Cost | Comment |
/// |---|---|---|
/// | Advection (1) | ~ 50 µs | unchanged from Phase 1.1 |
/// | DS + EH (2, 3) | ~ 7 µs | small per-cell formulas |
/// | Isostasy (4a) | ~ 150 µs | Gaussian blur σ = 2 dominates |
/// | Sea-level (4b) | ~ 5 µs | single min/max pass |
/// | Drainage targets (4c) | ~ 100 µs | BFS bounded by `max_distance` |
/// | Drainage areas (4d) | ~ 30 µs | `O(N log N)` sort |
/// | Erosion (4e) | ~ 30 µs | linear scan with skip-if-flat |
/// | **Total** | **~ 372 µs** | matches measured 367 µs/step |
///
/// When `closures.erosion.enabled = false`, stages 4a-4e are
/// **skipped entirely** (single branch at the top of the
/// erosion block — not the closure's internal early-return).
/// In this configuration the loop reduces to Phase 1.3
/// behaviour bit-identically — see § Closure isolation below.
///
/// ## Closure isolation (W4 discipline)
///
/// Each closure has its own `enabled` flag. The time loop honours
/// these flags at two levels:
///
/// - **Davis-Suppe** (stage 2) and **equilibrium-height** (stage
///   3): the `apply_*_step` functions early-return on `!enabled`.
///   The per-step overhead is a single branch comparison.
/// - **Erosion** (stage 4): the entire `if closures.erosion.
///   enabled { … }` block (4a-4e) is skipped. The expensive
///   isostasy + drainage precomputation does NOT run when erosion
///   is off. This preserves bit-identical regression for
///   Phase 1.2 / 1.3 tests that disable erosion explicitly via
///   `erosion: ErosionParams { enabled: false, .. }`.
///
/// The bit-identical decomposition contract — pinned by
/// `c1_phase_a_decomposes_into_closures_then_post_tectonic` in
/// `crates/ymir-core/tests/c1_phase_1_3_workflow.rs` — depends on
/// this two-level isolation: with erosion disabled, the wrapper
/// and the manual `run_with_closures + apply_post_tectonic`
/// decomposition produce byte-identical `S̃` buffers.
///
/// ## Static-classification optimisation (Phase 1.2)
///
/// Since `plate_id` is not advected, the boundary classification
/// and the intra-plate wedge-distance field are static throughout
/// the run. They are computed **once before** the loop and reused
/// every step. Phase 2 (boundary evolution) will move them back
/// inside the loop when the underlying geometry actually changes
/// per step. The erosion stage 4 has no equivalent static
/// pre-computation today — `S̃` mutates every step, so altitude,
/// drainage targets, and drainage areas must be re-derived.
///
/// ## See also
///
/// - [`crate::tectonics_v2::workflow::phase_a_common::apply_post_tectonic`]
///   — the end-of-cycle post-tectonic pass that wraps this loop
///   in the C1 workflow (and the v2 workflow under `v2_legacy`).
/// - [`crate::tectonics_v2::workflow::phase_a_c1::run_phase_a_cycle_c1`]
///   — the C1 paradigm Phase A entry that orchestrates this loop
///   + `apply_post_tectonic`.
/// - `docs/c1_lightweight_dynamic_tectonics.md` §5.1 (MVP
///   closures), §7.1 (Phase 1 prototype plan), §11
///   (implicit physical scales).
pub fn run_with_closures<F>(
    state: &mut C1State,
    kinematics: &mut PlateKinematics,
    config: &C1TimeLoopConfig,
    closures: &C1Closures,
    mut on_step: F,
) where
    F: FnMut(usize, &C1State),
{
    let nx = state.nx();
    let ny = state.ny();
    let n_cells = nx * ny;

    // #147 mesh invariance — rescale the Davis-Suppe length scales
    // from their reference-grid (64²) calibration to the actual grid,
    // so the wedge has a FIXED PHYSICAL width independent of the mesh.
    // Used for BOTH the upstream wedge-distance `max_distance` and the
    // kernel, so they stay consistent. At nx==64 this is the identity.
    let davis_suppe = closures.davis_suppe.scaled_to_grid(nx);

    let max_v = kinematics.max_velocity().max(1e-12);
    let dx_min = config.dx.min(config.dy);
    let dt = 0.5 * dx_min / max_v;

    let mut vx = vec![0.0_f64; n_cells];
    let mut vy = vec![0.0_f64; n_cells];
    fill_velocity_field(&mut vx, &mut vy, state, kinematics);
    // #145 — rigid mask for the no-flux boundary (rebuilt per-step on Track D).
    let mut rigid_mask: Option<Vec<bool>> = if config.rigid_continental_crust {
        apply_continental_rigidity(&mut vx, &mut vy, state);
        Some(continental_rigid_mask(state))
    } else {
        None
    };

    // Track D enabled-flag drives per-step boundary / wedge_d
    // recompute (Stage E4 Q-E3.1). When any Track D closure is
    // enabled, `plate_id` may mutate per step (subduction floor-
    // trigger reassignment, accretion merge, rifting split), so
    // the static cache of `boundary` + `wedge_d` becomes stale.
    // Recompute every step costs ~ 200 µs at 64² — well inside
    // the Phase 2 W4 budget.
    let any_track_d_enabled = closures.subduction.enabled
        || closures.accretion.enabled
        || closures.rifting.enabled;

    // Outer cache for the Track-D-disabled path (Phase 1.x +
    // Track A/B regression — `plate_id` is invariant, the cache
    // stays valid for the full run).
    let outer_cache: Option<(_, _)> = if !any_track_d_enabled {
        let bi = classify_boundaries(&state.plate_id, kinematics);
        let wd = wedge_distance_intra_plate(
            &state.plate_id,
            &bi.upper_plate_mask,
            davis_suppe.max_distance,
        );
        Some((bi, wd))
    } else {
        None
    };

    // Trackers — internal to the run when Track D is enabled
    // (Q-E1.2 Option (c) confirmed). Dropped at function exit;
    // not part of the saveable `C1State`.
    let mut convergence_tracker = if closures.accretion.enabled {
        Some(ConvergenceTracker::new())
    } else {
        None
    };
    let mut divergence_tracker = if closures.rifting.enabled {
        Some(DivergenceTracker::new())
    } else {
        None
    };

    let mut s_next = Field2D::new(nx, ny);
    let mut age_next = Field2D::new(nx, ny);
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    for step in 0..config.n_steps {
        // Per-step boundary geometry. When Track D is enabled,
        // recompute from the current `plate_id` (which may have
        // been mutated by the previous step's Track D events).
        // When Track D is disabled, reuse the outer cache.
        let per_step_cache: Option<(_, _)> = if any_track_d_enabled {
            let bi = classify_boundaries(&state.plate_id, kinematics);
            let wd = wedge_distance_intra_plate(
                &state.plate_id,
                &bi.upper_plate_mask,
                davis_suppe.max_distance,
            );
            Some((bi, wd))
        } else {
            None
        };
        let (boundary, wedge_d): (&_, &_) = match (&per_step_cache, &outer_cache) {
            (Some((b, w)), _) => (b, w),
            (None, Some((b, w))) => (b, w),
            (None, None) => unreachable!(
                "boundary cache must be present: either per-step (Track D) or outer (cached)"
            ),
        };

        // Per-step velocity field. When Track D may have mutated
        // `kinematics` or `plate_id` in the previous step, rebuild.
        // Phase 1.x / Track A/B path keeps the initial fill (no
        // mutation possible).
        if any_track_d_enabled {
            fill_velocity_field(&mut vx, &mut vy, state, kinematics);
            // Track D mutates plate_type (subduction Oceanic→Continental,
            // accretion) → re-apply rigidity + rebuild the rigid mask so
            // newly-continental crust becomes rigid this step (point 3).
            if config.rigid_continental_crust {
                apply_continental_rigidity(&mut vx, &mut vy, state);
                rigid_mask = Some(continental_rigid_mask(state));
            }
        }

        // 1. Advection (no-flux rigid boundary when rigid; #145).
        match &rigid_mask {
            Some(m) => {
                step_upwind_masked(
                    nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.s, &vx, &vy, m,
                    &mut s_next,
                );
                step_upwind_masked(
                    nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.age, &vx, &vy, m,
                    &mut age_next,
                );
            }
            None => {
                step_upwind(
                    nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.s, &vx, &vy,
                    &mut s_next,
                );
                step_upwind(
                    nx, ny, config.dx, config.dy, dt, &idx_x, &idx_y, &state.age, &vx, &vy,
                    &mut age_next,
                );
            }
        }
        std::mem::swap(&mut state.s, &mut s_next);
        std::mem::swap(&mut state.age, &mut age_next);

        // 2. Davis-Suppe orogenic source term (Phase 1.2).
        apply_davis_suppe_step(
            &mut state.s,
            &state.plate_id,
            boundary,
            wedge_d,
            kinematics,
            &davis_suppe,
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
        // Stage 4 — altitude preparation (isostasy + S-S bathymetry)
        // is shared between the erosion path and the oceanic-
        // bathymetry path. Both closures consume altitude:
        //
        //  - Erosion (4c-4e) reads altitude for the slope-magnitude
        //    factor of W-T `E = K · A^m · S^n`.
        //  - Stein-Stein (4b) reads age + plate_type to OVERWRITE
        //    altitude on oceanic cells with `−d(t) / depth_scale_m`
        //    (Architecture C — see [`apply_stein_stein_bathymetry`]
        //    sign-convention paragraph).
        //
        // If both are enabled, S-S runs FIRST so the drainage /
        // erosion sub-pipeline sees S-S-modulated altitude — the
        // intended cross-closure interaction (older oceanic
        // plateaus drain toward ridges; ridges are local highs in
        // the oceanic basin).
        //
        // If only oceanic_bathymetry is enabled (erosion off), S-S
        // still runs each step but its effect is **transient** —
        // altitude is recomputed from S̃ via `compute_isostasy` at
        // the next call, so the closure modifies state only through
        // downstream consumers in the same stage. Acceptance tests
        // must re-apply S-S at the run boundary to observe its
        // imprint (see Stage A test 1 of `c1_phase_2_bathymetry_*`).
        if closures.erosion.enabled || closures.oceanic_bathymetry.enabled {
            let isostasy = compute_isostasy(&state.s, &config.iso_config);
            let mut altitude = isostasy.heightmap;

            // 4b. Stein-Stein 1992 oceanic bathymetry (Phase 2
            // Track A — Issue #129). Internal `enabled = false`
            // early-return so the call is free when the closure is
            // disabled.
            apply_stein_stein_bathymetry(
                &mut altitude,
                &state.age,
                &state.plate_type,
                &closures.oceanic_bathymetry,
            );

            // 4c-4e. Erosion sub-pipeline. Drainage classification
            // consumes `state.s` in S̃ space (sea-level threshold in
            // S̃ units, not altitude — so S-S adjustment does NOT
            // shift the drainage classification), but the erosion
            // formula's slope factor reads `altitude` directly and
            // therefore sees the S-S modulation.
            if closures.erosion.enabled {
                let sea_level_ref =
                    compute_sea_level_ref_s_space(&state.s, &config.iso_config);
                let drainage_map = compute_drainage_targets(
                    &state.s,
                    sea_level_ref,
                    config.drainage_max_distance,
                );
                let drainage_areas = compute_drainage_areas(&drainage_map);
                apply_erosion_step(
                    &mut state.s,
                    &altitude,
                    &drainage_areas,
                    &closures.erosion,
                    dt,
                    config.dx,
                );
            }
        }

        // 5. Track D boundary-evolution events (Issue #132).
        //
        // Order:
        //   5a. Subduction — consumes oceanic, distributes arc;
        //       may mutate `plate_id` + `plate_type` via floor-
        //       trigger reassignment.
        //   5b. Rifting thinning — negative `S̃` source on
        //       divergent continental cells.
        //   5c. Trackers update — re-scan post-subduction
        //       `plate_id` to identify currently-convergent and
        //       currently-divergent pairs.
        //   5d. Accretion merge — when convergence count ≥
        //       threshold, merge plate pair; mutates `plate_id`
        //       and `kinematics.velocities[winner]`.
        //   5e. Rifting split — when both time + thickness
        //       conditions hold, spawn new plate from rift strip;
        //       mutates `plate_id`, `kinematics.velocities`
        //       (push), `age` (Path 3.B = 0).
        //
        // `plate_type` is preserved per cell through accretion +
        // rifting events (continental cells stay continental;
        // oceanic stay oceanic). Subduction's floor-trigger
        // reassignment is the only path that re-types a cell, and
        // it does so in-place on the (i, j) being processed —
        // per-cell `plate_type` ↔ `plate_id` consistency is
        // preserved at all times. No "recompute plate_type from
        // plate_id" step is needed (the spec called it out as 5g
        // but it would be a no-op given how subduction handles
        // re-typing inline).
        // Capture Track D per-step diagnostics into `state.last_step_stats`
        // (Issue #137 Viz-D0 Option B). Default = all-zero when a closure
        // is disabled (the `apply_*_step` functions early-return their
        // `Default` stats on `!enabled`).
        let subduction_stats = if closures.subduction.enabled {
            apply_subduction_step(
                &mut state.s,
                &mut state.plate_id,
                &mut state.plate_type,
                boundary,
                kinematics,
                &closures.subduction,
                dt,
            )
        } else {
            Default::default()
        };

        let rifting_thinning_stats = if closures.rifting.enabled {
            apply_rifting_thinning(
                &mut state.s,
                &state.plate_type,
                &state.plate_id,
                boundary,
                kinematics,
                &closures.rifting,
                dt,
            )
        } else {
            Default::default()
        };

        if let Some(tracker) = convergence_tracker.as_mut() {
            tracker.update(&state.plate_id, kinematics);
        }
        if let Some(tracker) = divergence_tracker.as_mut() {
            tracker.update(&state.plate_id, kinematics);
        }

        let accretion_stats = if closures.accretion.enabled {
            apply_accretion_step(
                &mut state.plate_id,
                &state.s,
                kinematics,
                convergence_tracker
                    .as_ref()
                    .expect("convergence_tracker allocated when accretion enabled"),
                &closures.accretion,
            )
        } else {
            Default::default()
        };

        let rifting_split_stats = if closures.rifting.enabled {
            apply_rifting_split(
                &mut state.plate_id,
                &state.plate_type,
                &mut state.age,
                &state.s,
                kinematics,
                divergence_tracker
                    .as_ref()
                    .expect("divergence_tracker allocated when rifting enabled"),
                &closures.rifting,
            )
        } else {
            Default::default()
        };

        // Aggregate into `state.last_step_stats` BEFORE `on_step` so the
        // callback (and viz snapshot) sees the matching stats.
        state.last_step_stats = crate::tectonics_c1::stats::C1StepStats {
            subduction: subduction_stats,
            accretion: accretion_stats,
            rifting_thinning: rifting_thinning_stats,
            rifting_split: rifting_split_stats,
        };

        // Cratonic-factor staleness (Q-S.2 Option (c) confirmed).
        // The `cratonic_mask` BoolField on `C1State` was built by
        // BFS over the initial-Continental-seeded plates at init
        // time (`init.rs::build_phase_1_1_cratonic_mask`). It is
        // NOT recomputed when Track D mutates `plate_id`. Result:
        // a cell that was reassigned to a different continental
        // plate via subduction (or moved into a new rifted plate
        // via Path 3.B split) carries the cratonic flag of its
        // OLD plate identity until the end-of-cycle
        // `apply_post_tectonic` cratonic recompute runs (see
        // `workflow/phase_a_c1.rs::run_phase_a_cycle_c1`).
        //
        // The lag is at most one cycle (~0.67 Ma at Phase 1.1
        // `dt`). Acceptable per Stage S Q-S.2 trade-off:
        // cratonic mask is a slow-evolving feature (Earth
        // cratons preserved over ~100 Ma timescales), so a 0.67-
        // Ma lag is physically negligible. Alternative was a per-
        // step cratonic-mask rebuild — invasive and unprofiled.
        // Pattern matches Phase 1.4's "designed trade-off, not
        // bug" floor-clamp finding.

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
        C1State {
            s,
            age,
            plate_id,
            plate_type,
            cratonic_mask,
            num_plates: 1,
            last_step_stats: Default::default(),
        }
    }

    #[test]
    fn uniform_field_advection_conserves_mass_exactly() {
        let nx = 16;
        let ny = 16;
        let mut state = uniform_single_plate_state(nx, ny);
        let initial_mass: f64 = state.s.data().iter().sum();

        let kinematics = PlateKinematics { velocities: vec![(0.01, 0.005)] };
        let config = C1TimeLoopConfig {
        rigid_continental_crust: false,
            n_steps: 100,
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            iso_config: IsostasyConfig::default(),
            drainage_max_distance: 30,
        };

        // `run_advection_only` keeps the `&PlateKinematics`
        // signature — only `run_with_closures` was upgraded to
        // `&mut` for Track D.
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
        rigid_continental_crust: false,
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
