//! Davis-Suppe orogenic source term applied per time step.
//!
//! ## Formula
//!
//! For each interior cell `c` belonging to the wedge body of an
//! upper plate (`d > 0` from a same-plate upper-plate boundary
//! seed; `d < max_distance`; `boundary_type[c] != Convergent`):
//!
//! ```text
//!     ∂S̃/∂t  =  coupling · |v_plate(c)|
//!               · max(0, h_critical(d) − S̃(c))
//!               · exp(−d / L_decay)
//! ```
//!
//! where:
//!
//! - `h_critical(d) = h_max · (1 − exp(−d / L_taper))` — the
//!   classical Davis-Suppe wedge profile, growing from 0 at the
//!   boundary toward `h_max` deep in the wedge.
//! - `|v_plate(c)|` is the cell's plate-velocity magnitude. Phase
//!   1.2 uses this as the convergence-rate proxy; Phase 3
//!   Lallemand will refine to the true relative-velocity normal
//!   component at the closest boundary.
//! - The `exp(−d / L_decay)` factor decays the source away from
//!   the boundary so the wedge stays localised.
//!
//! ## Why `boundary_type != Convergent` is skipped
//!
//! Stage 3.1's intra-plate wedge distance gives `d = 0` to every
//! `BoundaryType::Convergent` cell on the upper-plate side. The
//! `h_critical(0) = 0` value would make the source term go
//! negative (`h_crit − h_continental ≈ −1`), thinning the
//! boundary instead of thickening the interior — anti-geological.
//! The skip is the architectural fix surfaced before Stage 4 and
//! is locked in by [`tests::source_term_skips_boundary_cells`].
//!
//! ## Mass conservation
//!
//! Source-only by design (W3 of Issue #123). Mass grows under
//! this closure; Phase 1.4's stream-power erosion will be the
//! corresponding sink. Phase 1.2 outputs should be read with
//! that in mind: the test asserts wedge formation, not mass
//! balance.

use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::PlateIdField;

use crate::tectonics_c1::boundary_classification::{BoundaryInfo, BoundaryType};
use crate::tectonics_c1::kinematics::PlateKinematics;

/// Davis-Suppe orogenic closure tunables.
///
/// Defaults selected for Phase 1.2 visual demonstration at 64²,
/// `dt ≈ 0.69` non-dim/step (Phase 1.1 timestep), 300-step run.
/// See `apply_davis_suppe_step` doc for the calibration argument.
#[derive(Clone, Copy, Debug)]
pub struct DavisSuppeParams {
    /// Master enable/disable. When `false`, [`apply_davis_suppe_step`]
    /// is a no-op (W4 watchpoint: closure-disabled regression must
    /// reproduce Phase 1.1 behaviour bit-identically).
    pub enabled: bool,
    /// Source-term coupling constant. Units `1 / (velocity ·
    /// thickness · time)` in the implementation; calibrated
    /// dimensionlessly here.
    pub coupling: f64,
    /// Plateau height for the wedge: `h_critical(∞) = h_max`.
    pub h_max: f64,
    /// Characteristic length over which `h_critical(d)` rises
    /// from 0 to `h_max`. In **cells at [`DS_REFERENCE_GRID`]** —
    /// the runtime rescales it ∝ `nx` via [`DavisSuppeParams::scaled_to_grid`]
    /// so the wedge has a FIXED PHYSICAL width independent of the
    /// mesh (Issue #147 mesh invariance). At `nx == 64` the rescale
    /// is the identity.
    pub l_taper: f64,
    /// Characteristic length over which the source-term decays
    /// away from the boundary. In **cells at [`DS_REFERENCE_GRID`]**
    /// (see `l_taper`; rescaled ∝ `nx` at runtime).
    pub l_decay: f64,
    /// Outside this distance the source term is zero. In **cells at
    /// [`DS_REFERENCE_GRID`]** (rescaled ∝ `nx` at runtime). Must
    /// match the `max_distance` argument used when computing the
    /// wedge distance field upstream — both come from the SAME
    /// [`DavisSuppeParams::scaled_to_grid`] result.
    pub max_distance: f64,
}

/// Reference resolution at which the Davis-Suppe length scales
/// (`l_taper`, `l_decay`, `max_distance`) were calibrated. The
/// runtime treats the stored values as *physical* widths expressed
/// in reference-grid cells: a length `l` cells at `DS_REFERENCE_GRID`
/// is the physical fraction `l / DS_REFERENCE_GRID` of the domain,
/// re-sampled to `l · nx / DS_REFERENCE_GRID` cells at the actual
/// grid `nx` (Issue #147). This makes the wedge width mesh-invariant
/// while preserving the calibrated 64² behaviour bit-for-bit.
pub const DS_REFERENCE_GRID: f64 = 64.0;

impl DavisSuppeParams {
    /// Rescale the cell-valued length parameters from
    /// [`DS_REFERENCE_GRID`] to the actual grid `nx`, yielding a
    /// FIXED PHYSICAL wedge width independent of the mesh (Issue
    /// #147). `factor = nx / DS_REFERENCE_GRID`; at `nx == 64` this
    /// is the identity (byte-identical to pre-#147). Non-length
    /// fields (`coupling`, `h_max`, `enabled`) are unchanged — they
    /// are intensive, not lengths.
    #[must_use]
    pub fn scaled_to_grid(&self, nx: usize) -> DavisSuppeParams {
        let factor = nx as f64 / DS_REFERENCE_GRID;
        DavisSuppeParams {
            l_taper: self.l_taper * factor,
            l_decay: self.l_decay * factor,
            max_distance: self.max_distance * factor,
            ..*self
        }
    }
}

impl Default for DavisSuppeParams {
    fn default() -> Self {
        Self {
            enabled: true,
            coupling: 2.0,
            h_max: 2.5,
            l_taper: 4.0,
            l_decay: 6.0,
            max_distance: 30.0, // = 5 × l_decay
        }
    }
}

/// Davis-Suppe critical taper profile: `h_critical(d) = h_max · (1
/// − exp(−d / L_taper))`. Exposed for diagnostics and unit tests.
#[inline]
pub fn h_critical(d: f64, params: &DavisSuppeParams) -> f64 {
    params.h_max * (1.0 - (-d / params.l_taper).exp())
}

/// Apply one forward-Euler step of the Davis-Suppe orogenic
/// closure to the `S̃` field.
///
/// Inputs (all per-cell aligned to the same `nx × ny` grid):
/// - `s`: mutable `S̃` field, updated in place
/// - `plate_id`: each cell's plate index (used for the `v_plate`
///   proxy)
/// - `boundary`: Stage 2 classification — `boundary_type` is read
///   to skip Convergent cells; `upper_plate_mask` is **not** read
///   here (the wedge body is identified by `wedge_distance > 0`,
///   filtered by the intra-plate Dijkstra upstream)
/// - `wedge_distance`: Stage 3.1 intra-plate distance field
/// - `kinematics`: per-plate velocity bundle
/// - `params`: closure tunables
/// - `dt`: time step in the same units as `coupling × v` produces
pub fn apply_davis_suppe_step(
    s: &mut Field2D,
    plate_id: &PlateIdField,
    boundary: &BoundaryInfo,
    wedge_distance: &Field2D,
    kinematics: &PlateKinematics,
    params: &DavisSuppeParams,
    dt: f64,
) {
    apply_davis_suppe_step_inner(
        s, plate_id, boundary, wedge_distance, kinematics, params, dt, None,
    );
}

/// #155 maillon 1b-i — geometry-routed Davis-Suppe step. Identical to
/// [`apply_davis_suppe_step`] except cells flagged in `oc_wedge` (their
/// nearest upper-plate seed is an O-C subduction; see
/// [`super::super::distance_field::wedge_distance_intra_plate_typed`])
/// use the **margin-peaked ridge** target `h_critical_oc(d) = h_max ·
/// exp(−d / l_taper)` — peaking near the margin and decaying inland
/// (Andes) — instead of the rising-to-plateau dome (Tibet) the C-C /
/// velocity-fallback cells keep. Reuses `l_taper` as the ridge length
/// (no new knob): the O-C profile is the literal inverse of the C-C one.
pub fn apply_davis_suppe_step_routed(
    s: &mut Field2D,
    plate_id: &PlateIdField,
    boundary: &BoundaryInfo,
    wedge_distance: &Field2D,
    kinematics: &PlateKinematics,
    params: &DavisSuppeParams,
    dt: f64,
    oc_wedge: &crate::tectonics_c1::state::BoolField,
) {
    apply_davis_suppe_step_inner(
        s, plate_id, boundary, wedge_distance, kinematics, params, dt, Some(oc_wedge),
    );
}

#[allow(clippy::too_many_arguments)]
fn apply_davis_suppe_step_inner(
    s: &mut Field2D,
    plate_id: &PlateIdField,
    boundary: &BoundaryInfo,
    wedge_distance: &Field2D,
    kinematics: &PlateKinematics,
    params: &DavisSuppeParams,
    dt: f64,
    oc_wedge: Option<&crate::tectonics_c1::state::BoolField>,
) {
    if !params.enabled {
        return;
    }
    let nx = s.nx();
    let ny = s.ny();
    let h_max = params.h_max;
    let l_taper = params.l_taper;
    let l_decay = params.l_decay;
    let coupling = params.coupling;
    let max_d = params.max_distance;
    for j in 0..ny {
        for i in 0..nx {
            // Architectural skip — boundary cells must NOT receive
            // the source term (W4 lock; see module docstring + the
            // `source_term_skips_boundary_cells` test).
            if matches!(boundary.boundary_type.get(i, j), BoundaryType::Convergent) {
                continue;
            }
            let d = wedge_distance.get(i, j);
            // Out-of-reach cells: either too far from any same-plate
            // upper-plate seed, or the plate has no upper seed at
            // all (silent plate per Stage 3.1 finding).
            if d >= max_d {
                continue;
            }
            // h_critical and current S̃. #155 1b-i: O-C cells (nearest
            // upper-plate seed is a subduction) use the margin-peaked
            // ridge target h_max·exp(−d/l_taper) — peaks near the margin,
            // decays inland (Andes); C-C / fallback keep the classical
            // rising-to-plateau h_max·(1−exp(−d/l_taper)) (Tibet).
            let is_oc = oc_wedge.is_some_and(|m| m.get(i, j));
            let h_crit = if is_oc {
                h_max * (-d / l_taper).exp()
            } else {
                h_max * (1.0 - (-d / l_taper).exp())
            };
            let h_now = s.get(i, j);
            // Saturation: no source when already at or above
            // critical taper — `max(0, h_crit - h_now)`.
            let driving = h_crit - h_now;
            if driving <= 0.0 {
                continue;
            }
            // Convergence-rate proxy: cell plate's velocity
            // magnitude. Phase 3 refines.
            let pid = plate_id.get(i, j) as usize;
            let (vx, vy) = kinematics.velocities[pid];
            let v_mag = (vx * vx + vy * vy).sqrt();
            // Decay envelope.
            let envelope = (-d / l_decay).exp();
            // Forward-Euler step.
            let ds_dt = coupling * v_mag * driving * envelope;
            s.set(i, j, h_now + ds_dt * dt);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::tectonics_c1::boundary_classification::{BoundaryInfo, BoundaryTypeField};
    use crate::tectonics_c1::state::BoolField;

    /// Build the minimal triplet needed by `apply_davis_suppe_step`
    /// for a single-cell hand-crafted test scenario.
    ///
    /// Layout: 4×4 grid, one plate (plate 0), one cell flagged
    /// Convergent at the eastern edge, the target test cell at
    /// `(target_i, target_j)` with the configured wedge distance.
    fn build_test_scenario(
        nx: usize,
        ny: usize,
        target_i: usize,
        target_j: usize,
        wedge_d: f64,
        initial_s: f64,
    ) -> (
        Field2D,
        PlateIdField,
        BoundaryInfo,
        Field2D,
        PlateKinematics,
    ) {
        let mut s = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                s.set(i, j, initial_s);
            }
        }
        let plate_id = PlateIdField::new(nx, ny); // all 0
        // Boundary: a single Convergent cell at (nx-1, target_j)
        // so the test target is a wedge-body cell.
        let mut bt = BoundaryTypeField::filled(nx, ny, BoundaryType::Internal);
        bt.set(nx - 1, target_j, BoundaryType::Convergent);
        let upper = BoolField::filled(nx, ny, false);
        let boundary = BoundaryInfo {
            boundary_type: bt,
            upper_plate_mask: upper,
        };
        // Wedge distance: max everywhere except the target cell
        // (= wedge_d) and the Convergent cell (= 0). The
        // `apply_davis_suppe_step` only reads the cells it visits
        // so other entries don't matter for these tests.
        let mut wd = Field2D::filled(nx, ny, 50.0);
        wd.set(target_i, target_j, wedge_d);
        wd.set(nx - 1, target_j, 0.0);
        let kin = PlateKinematics {
            velocities: vec![(0.01, 0.0)], // plate 0 at |v| = 0.01
        };
        (s, plate_id, boundary, wd, kin)
    }

    #[test]
    fn source_term_bounded_at_critical() {
        // Cell at d = 4 with S̃ already at h_critical(4). The
        // source term should produce zero change.
        let d = 4.0_f64;
        let params = DavisSuppeParams::default();
        let h_crit_target = h_critical(d, &params);
        let (mut s, plate_id, boundary, wd, kin) =
            build_test_scenario(8, 8, 3, 4, d, h_crit_target);
        let s_before = s.get(3, 4);
        apply_davis_suppe_step(&mut s, &plate_id, &boundary, &wd, &kin, &params, 1.0);
        let s_after = s.get(3, 4);
        assert!(
            (s_after - s_before).abs() < 1e-12,
            "S̃ at h_crit should not change; before={s_before}, after={s_after}"
        );
    }

    #[test]
    fn source_term_active_below_critical() {
        // Cell at d = 4 with S̃ = 0.2 (oceanic initial). Source
        // should be positive and approach h_critical(4).
        let d = 4.0_f64;
        let initial_s = 0.2;
        let params = DavisSuppeParams::default();
        let (mut s, plate_id, boundary, wd, kin) =
            build_test_scenario(8, 8, 3, 4, d, initial_s);
        apply_davis_suppe_step(&mut s, &plate_id, &boundary, &wd, &kin, &params, 1.0);
        let s_after = s.get(3, 4);
        assert!(
            s_after > initial_s,
            "S̃ below h_crit should grow; before={initial_s}, after={s_after}"
        );
        let h_crit_target = h_critical(d, &params);
        // Forward Euler should not overshoot for one step at this
        // coupling × dt = 2.0 × 1.0 = 2 (the relaxation factor is
        // coupling × v × decay × dt = 2 × 0.01 × exp(-4/6) ~ 0.010,
        // well below 1 → linear regime).
        assert!(
            s_after < h_crit_target,
            "single-step relaxation should not overshoot h_crit; after={s_after}, h_crit={h_crit_target}"
        );
    }

    #[test]
    fn source_term_disabled_when_flag_off() {
        // params.enabled = false → no change anywhere, regardless
        // of distance or initial S̃.
        let d = 4.0_f64;
        let params = DavisSuppeParams { enabled: false, ..DavisSuppeParams::default() };
        let initial_s = 0.2;
        let (mut s, plate_id, boundary, wd, kin) =
            build_test_scenario(8, 8, 3, 4, d, initial_s);
        let s_before = s.data().to_vec();
        apply_davis_suppe_step(&mut s, &plate_id, &boundary, &wd, &kin, &params, 1.0);
        let s_after = s.data();
        for k in 0..s_after.len() {
            assert_eq!(
                s_before[k], s_after[k],
                "disabled closure must not touch any cell; mismatch at flat index {k}"
            );
        }
    }

    #[test]
    fn source_term_skips_boundary_cells() {
        // Architectural lock from Stage 3.1: cells classified
        // Convergent must NOT receive the source term, even when
        // wedge_distance(cell) = 0. Verifies the skip and protects
        // the anti-thinning fix.
        let nx = 8;
        let ny = 8;
        let conv_i = nx - 1; // Convergent cell at (7, 4)
        let conv_j = 4;
        let initial_s = 1.0;
        // The Convergent cell starts at h = 1 with d = 0 — without
        // the skip, the source would compute h_crit(0) = 0 and
        // (h_crit - h) = -1 < 0, yielding ds_dt = 0 anyway by the
        // `driving <= 0` short-circuit. To make the test
        // unambiguous, set initial_s = 0 on the Convergent cell
        // so the source would have been positive without the skip.
        let mut s = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                s.set(i, j, initial_s);
            }
        }
        s.set(conv_i, conv_j, 0.0); // force a non-saturated state
        let plate_id = PlateIdField::new(nx, ny);
        let mut bt = BoundaryTypeField::filled(nx, ny, BoundaryType::Internal);
        bt.set(conv_i, conv_j, BoundaryType::Convergent);
        let upper = BoolField::filled(nx, ny, false);
        let boundary = BoundaryInfo {
            boundary_type: bt,
            upper_plate_mask: upper,
        };
        // Without the skip, a Convergent cell with d = 0 would
        // see: h_crit(0) = 0, driving = 0 - 0 = 0 → still 0.
        // So we force d = 1 on the Convergent cell to manufacture
        // a positive driving term. This is artificial — the
        // intra-plate Dijkstra would not put d > 0 on a seed
        // cell — but it exposes the skip's responsibility
        // independently.
        let mut wd = Field2D::filled(nx, ny, 50.0);
        wd.set(conv_i, conv_j, 1.0);
        let kin = PlateKinematics { velocities: vec![(0.01, 0.0)] };
        let params = DavisSuppeParams::default();
        apply_davis_suppe_step(&mut s, &plate_id, &boundary, &wd, &kin, &params, 1.0);
        let s_after = s.get(conv_i, conv_j);
        assert_eq!(
            s_after, 0.0,
            "Convergent boundary cell must not be modified by Davis-Suppe; got {s_after}"
        );
    }

    /// One-step dry-run on the real Phase 1.1 init state. Reports
    /// the source-term coverage metrics needed to calibrate the
    /// `coupling` value before the Stage 5 300-step run.
    #[test]
    fn phase_1_1_one_step_dry_run_calibration() {
        use crate::tectonics_c1::boundary_classification::classify_boundaries;
        use crate::tectonics_c1::distance_field::wedge_distance_intra_plate;
        use crate::tectonics_c1::init::init_c1_state_phase_1_1;

        let state = init_c1_state_phase_1_1(64, 42);
        let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let boundary = classify_boundaries(&state.plate_id, &kinematics);
        let params = DavisSuppeParams::default();
        let wedge_d = wedge_distance_intra_plate(
            &state.plate_id,
            &boundary.upper_plate_mask,
            params.max_distance,
        );

        // Snapshot S̃ before, run one step, diff.
        let dt = 0.69; // Phase 1.1 timestep
        let s_before = state.s.data().to_vec();
        let mut s = state.s.clone();
        apply_davis_suppe_step(
            &mut s,
            &state.plate_id,
            &boundary,
            &wedge_d,
            &kinematics,
            &params,
            dt,
        );
        let s_after = s.data();

        // Per-cell deltas.
        let mut active_count = 0_usize;
        let mut mass_added = 0.0_f64;
        let mut max_source = 0.0_f64;
        for k in 0..s_after.len() {
            let delta = s_after[k] - s_before[k];
            if delta > 0.0 {
                active_count += 1;
                mass_added += delta;
                if delta > max_source {
                    max_source = delta;
                }
            }
        }
        let mean_source = if active_count > 0 {
            mass_added / active_count as f64
        } else {
            0.0
        };
        let total_cells = s_after.len();
        let pct = |n: usize| 100.0 * n as f64 / total_cells as f64;

        eprintln!(
            "Phase 1.2 Stage 4 dry-run (1 step, dt={dt:.2}, coupling={:.2}, h_max={:.2}, L_taper={:.1}, L_decay={:.1}):",
            params.coupling, params.h_max, params.l_taper, params.l_decay
        );
        eprintln!(
            "  active cells   = {active_count} ({:.1} %)",
            pct(active_count)
        );
        eprintln!("  total mass added = {mass_added:.4}");
        eprintln!("  max source       = {max_source:.4}");
        eprintln!("  mean source      = {mean_source:.4}");

        // Project to 300 steps with naive linear extrapolation
        // (overestimate; the relaxation slows as h approaches h_crit).
        let projected_max_300 = max_source * 300.0;
        let projected_mean_300 = mean_source * 300.0;
        eprintln!(
            "  linear projection to 300 steps (upper bound — relaxation will slow):"
        );
        eprintln!(
            "    max single-cell  delta_max(300) ≈ {projected_max_300:.3}  (vs h_max = {:.2})",
            params.h_max
        );
        eprintln!(
            "    typical  delta_mean(300)        ≈ {projected_mean_300:.3}"
        );

        // Sanity asserts. At least some cells should be active
        // (the silent-plate finding tells us not 100% of cells, but
        // certainly > 0).
        assert!(active_count > 0, "expected some cells to receive source");
        // No cell should suddenly overshoot h_max in a single step
        // (forward-Euler stability check).
        assert!(
            max_source < params.h_max,
            "single-step delta {max_source} should not approach h_max {} — coupling × dt may be too high",
            params.h_max
        );
    }
}
