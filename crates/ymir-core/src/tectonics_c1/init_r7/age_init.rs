//! R7 ridge-aligned age = 0 initialisation. Phase 2 Track B
//! sub-component 3, **Path 3.A (init-only)** (Issue #131).
//!
//! See [`super`] for module-level rationale.
//!
//! ## Why this exists
//!
//! Phase 2 Track A
//! ([[c1-phase-2-track-a-outcomes]]) shipped Stein-Stein 1992
//! oceanic bathymetry under Architecture C and empirically
//! discovered that the C1 Phase 1.1 age field is advected as a
//! **density** (`∂_t·age + ∇·(age·v) = 0`) — same flux-form
//! upwind as `S̃`. After 300 steps from uniform oceanic
//! `age = 0.5`, the age distribution exhibits ~1000× pile-up at
//! convergent boundaries (max ≈ 6958, median ≈ 0). The S-S
//! closure is paper-faithful per individual age value (Stage V
//! anchor ±0.421 m), but the **input distribution** is
//! dominated by the pile-up artifact rather than a smooth
//! `√t → exp(-α·t)` profile.
//!
//! Path 3.A (this module) addresses the input distribution at
//! **init time only** without changing the advection PDE: detect
//! oceanic cells adjacent to divergent boundaries (where new
//! oceanic lithosphere physically forms at mid-ocean ridges),
//! set their age to 0. Non-boundary oceanic cells keep the
//! Phase 1.1 baseline. Per-step advection then carries the
//! `age = 0` ridge cells AWAY from the ridge (correct sign per
//! S-S 1992 — older oceanic cells should be deeper) — the
//! flux-form advection's density semantics still apply, but the
//! initial distribution starts geophysically sane.
//!
//! ## Trichotomy decision tree
//!
//! For each cell `c`:
//!
//! 1. `plate_type(c) == Continental` → `continental_baseline`
//!    (7.0). Continental cells keep the Phase 1.1 baseline
//!    regardless of boundary classification — Stein-Stein
//!    bathymetry is oceanic-only, and continental-rifting
//!    (McKenzie-Buck, design doc §5.2) is Phase 3 scope.
//! 2. `plate_type(c) == Oceanic` AND
//!    `boundary_type(c) == Divergent` → `ridge_value` (0.0).
//!    Oceanic cells on a divergent boundary sit at the
//!    "spreading ridge" interpretation — newly-formed oceanic
//!    lithosphere with zero age.
//! 3. Otherwise (oceanic, non-divergent boundary) →
//!    `oceanic_baseline` (0.5). Preserves Phase 1.1 baseline for
//!    interior oceanic cells.
//!
//! Precedence: continental > divergent > oceanic. Mixed cases
//! (e.g., a continental cell on a divergent boundary in a future
//! rifting scenario) fall through the precedence ladder
//! correctly without ambiguity.
//!
//! ## Strict vs loose ridge semantics
//!
//! Strict (Path 3.A this module): `age = 0` ONLY on cells with
//! `BoundaryType::Divergent` — typically 1–2 cells thick at
//! init. Per-step density-advection broadens the ridge zone
//! over time, replicating the expected geophysical picture
//! "ridge axis ages outward toward subduction zones".
//!
//! Loose (Path 3.B, NOT implemented here): `age = 0` on cells
//! within `k` BFS hops of a divergent boundary. Wider initial
//! ridge zone, lower per-step age-decay sensitivity. Deferred
//! pending Stage A acceptance test evidence.
//!
//! ## Fallback paths (NOT in scope this stage)
//!
//! - **Path 3.B — per-step ridge detection**. Run
//!   `classify_boundaries` each time step, set `age = 0` on
//!   newly-divergent cells. Couples age field with kinematics
//!   per step (more expensive but tracks dynamic boundary
//!   evolution; relevant once Track C is implemented).
//! - **Path 3.C — Lagrangian advection of age**. Replace the
//!   flux-form `∂_t·age + ∇·(age·v) = 0` with the
//!   non-conservative `∂_t·age + v·∇·age = 0` upwind. Different
//!   PDE: age moves with the velocity field as a per-cell
//!   attribute rather than a density. Eliminates pile-up at
//!   convergent boundaries entirely.
//!
//! Path 3.A ships in this PR. Escalation criterion for Stage A:
//!
//! - Stage A Spearman age-altitude correlation ≥ -0.4 (vs Track
//!   A baseline -0.476) → Path 3.A SUFFICIENT, ship.
//! - Stage A Spearman degraded → escalate to Path 3.B or 3.C
//!   per Phase 2 Track B-bis follow-up issue.

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::PlateIdField;

use super::super::boundary_classification::{BoundaryInfo, BoundaryType};

/// Parameters for ridge-aligned age initialisation.
///
/// Defaults preserve Phase 1.1 baselines `(continental = 7.0,
/// oceanic = 0.5)` and introduce `ridge_value = 0.0` for cells
/// classified as divergent-boundary by
/// [`crate::tectonics_c1::boundary_classification::classify_boundaries`].
///
/// Phase 1.1 baselines come from
/// `crate::tectonics_c1::init::{CONTINENTAL_AGE_INIT, OCEANIC_AGE_INIT}`
/// (7.0 and 0.5 non-dim respectively); preserved verbatim so
/// non-divergent oceanic cells get exactly the Phase 1.1 value.
#[derive(Clone, Copy, Debug)]
pub struct AgeInitParams {
    /// Age value for continental cells (any boundary
    /// classification). Default `7.0` per Phase 1.1
    /// `CONTINENTAL_AGE_INIT`.
    pub continental_baseline: f64,
    /// Age value for oceanic cells NOT on a divergent boundary.
    /// Default `0.5` per Phase 1.1 `OCEANIC_AGE_INIT`.
    pub oceanic_baseline: f64,
    /// Age value for oceanic cells classified as divergent-
    /// boundary by `classify_boundaries`. Default `0.0`
    /// (mid-ocean ridge axis, S-S 1992 `d = ridge_depth_m`).
    pub ridge_value: f64,
}

impl Default for AgeInitParams {
    fn default() -> Self {
        Self {
            continental_baseline: 7.0,
            oceanic_baseline: 0.5,
            ridge_value: 0.0,
        }
    }
}

/// Initialise the age field with ridge-aligned `age = 0` at
/// oceanic-divergent-boundary cells (Path 3.A).
///
/// Trichotomy decision tree applied per cell:
///
/// 1. `plate_type == Continental` → `params.continental_baseline`.
/// 2. `plate_type == Oceanic && boundary_type == Divergent` →
///    `params.ridge_value`.
/// 3. Otherwise → `params.oceanic_baseline`.
///
/// `boundary_info` is the caller-computed
/// [`BoundaryInfo`] from
/// [`crate::tectonics_c1::boundary_classification::classify_boundaries`].
/// Passing it instead of `kinematics` avoids re-running the
/// O(N · 4) classification when the caller already needs it for
/// other init work (Stage E4 dispatcher pattern).
///
/// Determinism: same `(plate_type, boundary_info, params)` →
/// bit-identical `age` field. No RNG; pure function over the
/// input fields.
///
/// `plate_id` is currently unused but kept on the signature for
/// future Path 3.B (per-step ridge detection) symmetry — same
/// signature shape eases the migration.
pub fn init_age_field_ridge_aligned(
    age: &mut Field2D,
    _plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    boundary_info: &BoundaryInfo,
    params: &AgeInitParams,
) {
    let nx = age.nx();
    let ny = age.ny();
    debug_assert_eq!(plate_type.nx(), nx, "plate_type nx mismatch");
    debug_assert_eq!(plate_type.ny(), ny, "plate_type ny mismatch");
    debug_assert_eq!(boundary_info.boundary_type.nx(), nx, "boundary nx mismatch");
    debug_assert_eq!(boundary_info.boundary_type.ny(), ny, "boundary ny mismatch");

    for j in 0..ny {
        for i in 0..nx {
            let value = match plate_type.get(i, j) {
                PlateType::Continental => params.continental_baseline,
                PlateType::Oceanic => {
                    match boundary_info.boundary_type.get(i, j) {
                        BoundaryType::Divergent => params.ridge_value,
                        _ => params.oceanic_baseline,
                    }
                }
            };
            age.set(i, j, value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_c1::boundary_classification::classify_boundaries;
    use crate::tectonics_c1::kinematics::PlateKinematics;

    /// Build a 2-plate horizontal-split `PlateIdField`:
    /// `plate_id = 0` for `i < nx/2`, `plate_id = 1` otherwise.
    fn two_plate_field(nx: usize, ny: usize) -> PlateIdField {
        let mut p = PlateIdField::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                p.set(i, j, if i < nx / 2 { 0 } else { 1 });
            }
        }
        p
    }

    /// Build a uniform `PlateTypeField` set to `PlateType::Oceanic`.
    fn all_oceanic(nx: usize, ny: usize) -> PlateTypeField {
        PlateTypeField::filled(nx, ny, PlateType::Oceanic)
    }

    /// Test 1 — age = 0 at divergent boundary only.
    ///
    /// 2 oceanic plates moving apart (plate 0 west, plate 1
    /// east). Columns nx/2 - 1 and nx/2 on the divergent
    /// boundary; remaining columns interior. Verify the strict
    /// 1-cell-thick semantics.
    #[test]
    fn age_init_ridge_aligned_at_divergent_boundary_only() {
        let nx = 16;
        let ny = 8;
        let plate_id = two_plate_field(nx, ny);
        let plate_type = all_oceanic(nx, ny);
        // Plate 0 west, plate 1 east → divergent across the
        // i = nx / 2 boundary.
        let kinematics = PlateKinematics {
            velocities: vec![(-0.02, 0.0), (0.01, 0.0)],
        };
        let boundary_info = classify_boundaries(&plate_id, &kinematics);
        let params = AgeInitParams::default();

        let mut age = Field2D::new(nx, ny);
        init_age_field_ridge_aligned(
            &mut age,
            &plate_id,
            &plate_type,
            &boundary_info,
            &params,
        );

        let left_edge = nx / 2 - 1;
        let right_edge = nx / 2;
        for j in 0..ny {
            assert_eq!(
                age.get(left_edge, j),
                params.ridge_value,
                "left-edge divergent oceanic cell ({left_edge}, {j}) must have ridge_value"
            );
            assert_eq!(
                age.get(right_edge, j),
                params.ridge_value,
                "right-edge divergent oceanic cell ({right_edge}, {j}) must have ridge_value"
            );
            // Interior cells (away from boundary AND wraparound)
            // get oceanic_baseline. Check col 3 (clearly interior
            // of plate 0) and col 12 (clearly interior of plate 1).
            assert_eq!(
                age.get(3, j),
                params.oceanic_baseline,
                "interior oceanic cell (3, {j}) must have oceanic_baseline"
            );
            assert_eq!(
                age.get(12, j),
                params.oceanic_baseline,
                "interior oceanic cell (12, {j}) must have oceanic_baseline"
            );
        }
    }

    /// Test 2 — no age = 0 at convergent boundary.
    ///
    /// Mirror of test 1 with plates converging. Boundary cells
    /// should be Convergent, not Divergent → no age = 0.
    #[test]
    fn age_init_no_age_zero_at_convergent_boundary() {
        let nx = 16;
        let ny = 8;
        let plate_id = two_plate_field(nx, ny);
        let plate_type = all_oceanic(nx, ny);
        // Plate 0 east, plate 1 west → CONVERGENT across
        // i = nx / 2 boundary.
        let kinematics = PlateKinematics {
            velocities: vec![(0.02, 0.0), (-0.01, 0.0)],
        };
        let boundary_info = classify_boundaries(&plate_id, &kinematics);
        let params = AgeInitParams::default();

        let mut age = Field2D::new(nx, ny);
        init_age_field_ridge_aligned(
            &mut age,
            &plate_id,
            &plate_type,
            &boundary_info,
            &params,
        );

        // No cell should have ridge_value = 0.0 in this scenario.
        // (Periodic wraparound at i = 0 / nx - 1 is Divergent in
        // mirror — that IS a divergent boundary geometrically. So
        // assert the DESIGNED convergent boundary at i = nx/2 has
        // NO ridge cells.)
        let left_edge = nx / 2 - 1;
        let right_edge = nx / 2;
        for j in 0..ny {
            assert_ne!(
                age.get(left_edge, j),
                params.ridge_value,
                "designed convergent boundary ({left_edge}, {j}) must NOT have ridge_value"
            );
            assert_ne!(
                age.get(right_edge, j),
                params.ridge_value,
                "designed convergent boundary ({right_edge}, {j}) must NOT have ridge_value"
            );
        }
    }

    /// Test 3 — baseline elsewhere.
    ///
    /// Continental plates with no divergent boundary in the
    /// fixture: every cell should get `continental_baseline`.
    /// Verifies the trichotomy precedence (continental > divergent
    /// > oceanic).
    #[test]
    fn age_init_baseline_elsewhere() {
        let nx = 16;
        let ny = 8;
        let plate_id = two_plate_field(nx, ny);
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        // Mixed kinematics — at least some divergent boundary
        // exists, but plate_type is Continental everywhere, so
        // continental_baseline must override per the trichotomy.
        let kinematics = PlateKinematics {
            velocities: vec![(-0.02, 0.0), (0.01, 0.0)],
        };
        let boundary_info = classify_boundaries(&plate_id, &kinematics);
        let params = AgeInitParams::default();

        let mut age = Field2D::new(nx, ny);
        init_age_field_ridge_aligned(
            &mut age,
            &plate_id,
            &plate_type,
            &boundary_info,
            &params,
        );

        for j in 0..ny {
            for i in 0..nx {
                assert_eq!(
                    age.get(i, j),
                    params.continental_baseline,
                    "continental cell ({i}, {j}) must get continental_baseline regardless of boundary classification"
                );
            }
        }
    }

    /// Test 4 — deterministic given identical input.
    #[test]
    fn age_init_deterministic() {
        let nx = 16;
        let ny = 8;
        let plate_id = two_plate_field(nx, ny);
        let plate_type = all_oceanic(nx, ny);
        let kinematics = PlateKinematics {
            velocities: vec![(-0.02, 0.0), (0.01, 0.0)],
        };
        let boundary_info = classify_boundaries(&plate_id, &kinematics);
        let params = AgeInitParams::default();

        let mut age_a = Field2D::new(nx, ny);
        let mut age_b = Field2D::new(nx, ny);
        init_age_field_ridge_aligned(
            &mut age_a,
            &plate_id,
            &plate_type,
            &boundary_info,
            &params,
        );
        init_age_field_ridge_aligned(
            &mut age_b,
            &plate_id,
            &plate_type,
            &boundary_info,
            &params,
        );

        for j in 0..ny {
            for i in 0..nx {
                assert_eq!(
                    age_a.get(i, j),
                    age_b.get(i, j),
                    "same input must produce bit-identical age field at ({i}, {j})"
                );
            }
        }
    }
}
