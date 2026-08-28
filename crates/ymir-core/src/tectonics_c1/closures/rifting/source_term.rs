//! Rifting thinning closure — negative `S̃` source on divergent
//! continental boundaries.
//!
//! ## Per-cell algorithm
//!
//! For each cell `c` where `boundary_type[c] == Divergent` AND
//! `plate_type[c] == Continental`:
//!
//! 1. Find the **most-divergent neighbour edge magnitude** —
//!    among the 4-connected neighbours with `plate_id != pid_c`,
//!    take the largest `|v_rel · n̂|` where `v_rel · n̂ < 0`
//!    (divergent sign per `boundary_classification`'s convention).
//! 2. Compute `Δs = thinning_rate × |divergence| × dt`,
//!    **clamped at `s_before`** so `S̃` cannot go negative.
//! 3. `s[c] -= Δs`.
//!
//! This is the **mirror of the Davis-Suppe orogenic source**
//! (Phase 1.2) but with the opposite sign and the opposite
//! boundary-type filter. Continental cells thin under sustained
//! divergent stretching, matching the McKenzie 1978 stretching-
//! factor formulation at the macro scale.
//!
//! Oceanic cells are NOT thinned — mid-ocean ridges create new
//! oceanic crust via a different mechanism (Track B Path 3.A
//! ridge-aligned `age = 0` init, plus the Track B continental
//! clustering that excludes ocean cells from the rift event).

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::voronoi::PlateIdField;

use crate::tectonics_c1::boundary_classification::{BoundaryInfo, BoundaryType};
use crate::tectonics_c1::kinematics::PlateKinematics;

use super::params::RiftingParams;

/// Per-step diagnostics for `apply_rifting_thinning`.
#[derive(Clone, Copy, Debug, Default)]
pub struct RiftingThinningStats {
    /// Number of continental-divergent cells that received a
    /// thinning increment this step.
    pub cells_thinned: usize,
    /// Total `Δs` applied (positive — represents the magnitude
    /// of mass REMOVED from continental S̃).
    pub total_mass_removed: f64,
}

/// Apply one forward-Euler step of the rifting thinning closure.
///
/// Mutates `s` only — `plate_id`, `plate_type`, and `kinematics`
/// are passed read-only. Returns
/// [`RiftingThinningStats`] for diagnostics + the Stage E4
/// mass-balance check.
///
/// Returns immediately with `RiftingThinningStats::default()` when
/// `params.enabled == false` — bit-identical no-op.
pub fn apply_rifting_thinning(
    s: &mut Field2D,
    plate_type: &PlateTypeField,
    plate_id: &PlateIdField,
    boundary_info: &BoundaryInfo,
    kinematics: &PlateKinematics,
    params: &RiftingParams,
    dt: f64,
) -> RiftingThinningStats {
    if !params.enabled {
        return RiftingThinningStats::default();
    }

    let nx = s.nx();
    let ny = s.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let neighbours: [(i32, i32, f64, f64); 4] =
        [(1, 0, 1.0, 0.0), (-1, 0, -1.0, 0.0), (0, 1, 0.0, 1.0), (0, -1, 0.0, -1.0)];

    let mut stats = RiftingThinningStats::default();

    for j in 0..ny {
        for i in 0..nx {
            if !matches!(boundary_info.boundary_type.get(i, j), BoundaryType::Divergent) {
                continue;
            }
            if plate_type.get(i, j) != PlateType::Continental {
                continue;
            }

            let pid_c = plate_id.get(i, j);
            let (vx_c, vy_c) = kinematics.velocities[pid_c as usize];

            // Find the most-divergent neighbour edge magnitude.
            // "Most divergent" = the neighbour with the largest
            // negative dot product (= largest positive |v_rel · n̂|
            // pointing AWAY from the boundary).
            let mut best_div_magnitude = 0.0_f64;
            for &(di, dj, nx_norm, ny_norm) in neighbours.iter() {
                let ni = if di > 0 {
                    idx_x.next(i)
                } else if di < 0 {
                    idx_x.prev(i)
                } else {
                    i
                };
                let nj = if dj > 0 {
                    idx_y.next(j)
                } else if dj < 0 {
                    idx_y.prev(j)
                } else {
                    j
                };
                let pid_n = plate_id.get(ni, nj);
                if pid_n == pid_c {
                    continue;
                }
                let (vx_n, vy_n) = kinematics.velocities[pid_n as usize];
                let vrel_x = vx_c - vx_n;
                let vrel_y = vy_c - vy_n;
                let dot = vrel_x * nx_norm + vrel_y * ny_norm;
                if dot < 0.0 {
                    let mag = -dot;
                    if mag > best_div_magnitude {
                        best_div_magnitude = mag;
                    }
                }
            }

            if best_div_magnitude <= 0.0 {
                // No divergent neighbour edge — verdict was
                // BoundaryType::Divergent because of a non-edge
                // case (rare).
                continue;
            }

            let delta_proposed = params.thinning_rate * best_div_magnitude * dt;
            let s_before = s.get(i, j);
            let delta = delta_proposed.min(s_before);
            s.set(i, j, s_before - delta);

            stats.cells_thinned += 1;
            stats.total_mass_removed += delta;
        }
    }

    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::tectonics_c1::boundary_classification::classify_boundaries;

    /// Three-plate east-west fixture mirroring the accretion tests'
    /// helper. Required because the 2-plate-on-torus pathology
    /// (Stage E2 finding) ties conv ↔ div counts and would make
    /// per-pair verdicts symmetric. `nx` must be divisible by 3.
    ///
    /// Returns `(s, plate_id, plate_type, kinematics, boundary_info)`
    /// — all the inputs `apply_rifting_thinning` needs.
    fn three_plate_fixture(
        nx: usize,
        ny: usize,
        plate_types: [PlateType; 3],
        velocities: [(f64, f64); 3],
    ) -> (Field2D, PlateIdField, PlateTypeField, PlateKinematics, BoundaryInfo) {
        assert_eq!(nx % 3, 0, "three_plate_fixture: nx must be divisible by 3");
        let mut s = Field2D::new(nx, ny);
        let mut plate_id = PlateIdField::new(nx, ny);
        let mut plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        let third = nx / 3;
        for j in 0..ny {
            for i in 0..nx {
                let p = if i < third {
                    0_usize
                } else if i < 2 * third {
                    1_usize
                } else {
                    2_usize
                };
                plate_id.set(i, j, p as u16);
                plate_type.set(i, j, plate_types[p]);
                let s_init = match plate_types[p] {
                    PlateType::Continental => 1.0,
                    PlateType::Oceanic => 0.2,
                };
                s.set(i, j, s_init);
            }
        }
        let kinematics = PlateKinematics { velocities: velocities.to_vec() };
        let boundary_info = classify_boundaries(&plate_id, &kinematics);
        (s, plate_id, plate_type, kinematics, boundary_info)
    }

    #[test]
    fn rifting_thinning_at_divergent_continental() {
        // Three continental plates. Pair (0, 1) at interior i = 2/3
        // diverges (plate 0 moves west, plate 1 moves east).
        // Continental cells at i = 2 and i = 3 should be thinned.
        let nx = 9;
        let ny = 4;
        let (mut s, plate_id, plate_type, kinematics, boundary_info) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        let params = RiftingParams::default();
        let dt = 0.69;

        let s_before_2 = s.get(2, 0);
        let s_before_3 = s.get(3, 0);

        let stats = apply_rifting_thinning(
            &mut s,
            &plate_type,
            &plate_id,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        assert!(
            stats.cells_thinned >= ny * 2,
            "expected ≥ {} cells thinned (2 boundary cols × {} rows), got {}",
            ny * 2,
            ny,
            stats.cells_thinned
        );
        assert!(stats.total_mass_removed > 0.0);
        assert!(
            s.get(2, 0) < s_before_2,
            "cell (2, 0) should be thinned: before {s_before_2}, after {}",
            s.get(2, 0)
        );
        assert!(
            s.get(3, 0) < s_before_3,
            "cell (3, 0) should be thinned: before {s_before_3}, after {}",
            s.get(3, 0)
        );
    }

    #[test]
    fn rifting_thinning_no_op_at_convergent() {
        // Pair (0, 1) at interior i = 2/3 CONVERGES (plate 0
        // east, plate 1 west). The convergent cells (i = 2 and
        // i = 3) must NOT be thinned.
        //
        // Note: other pairs in the three-plate fixture (the
        // (1, 2) interior and (0, 2) periodic wrap) may produce
        // divergent verdicts and trigger thinning elsewhere on
        // the grid — that's expected behaviour, not a regression.
        // The test asserts on the convergent boundary cells only.
        let nx = 9;
        let ny = 4;
        let (mut s, plate_id, plate_type, kinematics, boundary_info) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(0.01, 0.0), (-0.01, 0.0), (0.0, 0.0)],
        );
        let params = RiftingParams::default();
        let dt = 0.69;

        let s_2_before: Vec<f64> = (0..ny).map(|j| s.get(2, j)).collect();
        let s_3_before: Vec<f64> = (0..ny).map(|j| s.get(3, j)).collect();

        let _stats = apply_rifting_thinning(
            &mut s,
            &plate_type,
            &plate_id,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        for j in 0..ny {
            assert_eq!(
                s.get(2, j),
                s_2_before[j],
                "convergent cell (2, {j}) (plate 0 side) must not be thinned"
            );
            assert_eq!(
                s.get(3, j),
                s_3_before[j],
                "convergent cell (3, {j}) (plate 1 side) must not be thinned"
            );
        }
    }

    #[test]
    fn rifting_thinning_no_op_at_oceanic() {
        // Pair (0, 1) DIVERGES at i = 2/3 but plates 0 and 1 are
        // Oceanic — the Continental filter must reject these
        // cells. Plate 2 (Continental) sits at i = 6..8; its
        // boundary with plates 0 / 1 will produce divergent
        // verdicts on continental cells which DO get thinned —
        // that's correct behaviour, not a regression. Assertion
        // focuses on the oceanic divergent cells only.
        let nx = 9;
        let ny = 4;
        let (mut s, plate_id, plate_type, kinematics, boundary_info) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Oceanic, PlateType::Oceanic, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        let params = RiftingParams::default();
        let dt = 0.69;

        let s_2_before: Vec<f64> = (0..ny).map(|j| s.get(2, j)).collect();
        let s_3_before: Vec<f64> = (0..ny).map(|j| s.get(3, j)).collect();

        let _stats = apply_rifting_thinning(
            &mut s,
            &plate_type,
            &plate_id,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );

        for j in 0..ny {
            assert_eq!(
                s.get(2, j),
                s_2_before[j],
                "oceanic divergent cell (2, {j}) must not be thinned"
            );
            assert_eq!(
                s.get(3, j),
                s_3_before[j],
                "oceanic divergent cell (3, {j}) must not be thinned"
            );
        }
    }
}
