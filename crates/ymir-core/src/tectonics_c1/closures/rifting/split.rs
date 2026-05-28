//! Rifting split event — "chewing-gum cut" two-condition gate +
//! new plate_id allocation + Path 3.B event-driven `age = 0`.
//!
//! ## DivergenceTracker — symmetric mirror of ConvergenceTracker
//!
//! Per-pair counter keyed by canonical `(a, b)` with `a < b`,
//! incremented when the pair's boundary is net-divergent
//! (`div > conv`), reset to 0 otherwise. The only behavioural
//! difference from
//! [`crate::tectonics_c1::closures::accretion::ConvergenceTracker`]
//! is the sign of the verdict — everything else (canonical pair
//! ordering, `DOT_FLOOR`, reset-on-non-verdict semantics) is
//! identical. The Stage E2 architectural finding "2-plate-on-torus
//! pathology" applies here too: unit tests must use three_plate
//! fixtures.
//!
//! ## Chewing-gum cut — two-condition split gate (Q3.2)
//!
//! [`apply_rifting_split`] fires for a plate pair `(a, b)` when
//! BOTH conditions hold simultaneously:
//!
//! 1. **Time condition** — `divergence_tracker.get(a, b) >=
//!    split_time_threshold` (sustained extension).
//! 2. **Thickness condition** — the minimum `S̃` across the
//!    rift-strip cells of plate `a` (continental cells of `a`
//!    touching plate `b`) is `< split_thickness_threshold`
//!    (thinning below McKenzie β = 1.4 threshold).
//!
//! Either condition alone is insufficient. The intuition is
//! "stretched (time) + thinned (mass)" — both required for the
//! lithosphere to fail.
//!
//! ## Partitioning — pair.0 loses its rift strip
//!
//! For deterministic asymmetry, the **lower-index plate `a`
//! loses cells** to the newly-spawned plate. Specifically: the
//! continental cells of plate `a` directly adjacent (4-connected)
//! to any cell of plate `b` form the "rift strip" and are
//! reassigned to a new plate id `new_pid =
//! kinematics.velocities.len()`. This creates a thin geographic
//! sliver of `new_pid` along the original `(a, b)` boundary,
//! semantically a "nascent rift basin" spawned from `a`'s margin.
//!
//! Edge case — disconnected rift cells: the rift strip may be
//! geometrically discontinuous (multiple components if `a`'s
//! boundary with `b` wraps the periodic domain). The current
//! implementation reassigns ALL plate-a cells touching plate-b
//! regardless of connectivity, so the new plate may itself be
//! disconnected. This is **graceful behaviour, not a bug** — a
//! disconnected new plate represents a real geological pattern
//! ("rift basin in two pieces separated by a transform"). Stage
//! E4 visual review will surface if Stage A needs a connectivity-
//! enforcement refinement.
//!
//! ## Path 3.B — event-driven `age = 0`
//!
//! Track B Phase 2 introduced Path 3.A ridge-aligned `age = 0`
//! init (set at oceanic cells adjacent to divergent boundaries at
//! init time only). Path 3.B is the **event-driven extension**:
//! when a rifting split fires, every cell reassigned to the new
//! plate has its `age` set to `0`, simulating the freshly-formed
//! rift floor with no accumulated thermal age. This preserves the
//! Track B Spearman age-altitude correlation under Track D
//! mutation — without Path 3.B, the rifted-off cells would carry
//! their parent plate's advected age, producing implausible
//! "ancient rift floors".
//!
//! ## Plate-id allocation cap (W6)
//!
//! Splits are gated by `kinematics.velocities.len() <
//! params.plate_id_cap` (default `256`). When the cap is reached,
//! splits are refused gracefully — `stats.splits_count` does not
//! increment, no mutation. Surface architectural finding if Stage A
//! diagnostic shows runs approaching the cap.

use std::collections::{HashMap, HashSet};

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::voronoi::PlateIdField;

use crate::tectonics_c1::kinematics::PlateKinematics;

use super::params::RiftingParams;

const DOT_FLOOR: f64 = 1e-12;

/// Per-pair sustained-divergence counter — symmetric mirror of
/// [`crate::tectonics_c1::closures::accretion::ConvergenceTracker`].
///
/// See module docstring for the canonical ordering, update
/// semantics, and 2-plate-on-torus pathology caveat.
#[derive(Clone, Debug, Default)]
pub struct DivergenceTracker {
    pub divergence_counts: HashMap<(u16, u16), usize>,
}

impl DivergenceTracker {
    pub fn new() -> Self {
        Self { divergence_counts: HashMap::new() }
    }

    /// Re-scan the grid and update per-pair counters.
    ///
    /// Pairs in a net-divergent state this step (more divergent
    /// edges than convergent across their shared boundary) have
    /// their counter incremented by 1. All other tracked pairs are
    /// reset to 0.
    pub fn update(&mut self, plate_id: &PlateIdField, kinematics: &PlateKinematics) {
        let nx = plate_id.nx();
        let ny = plate_id.ny();
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        let mut pair_verdicts: HashMap<(u16, u16), (usize, usize)> = HashMap::new();

        let neighbours: [(i32, i32, f64, f64); 4] = [
            (1, 0, 1.0, 0.0),
            (-1, 0, -1.0, 0.0),
            (0, 1, 0.0, 1.0),
            (0, -1, 0.0, -1.0),
        ];

        for j in 0..ny {
            for i in 0..nx {
                let pid_c = plate_id.get(i, j);
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
                    let pair = if pid_c < pid_n {
                        (pid_c, pid_n)
                    } else {
                        (pid_n, pid_c)
                    };
                    let (vx_a, vy_a) = kinematics.velocities[pair.0 as usize];
                    let (vx_b, vy_b) = kinematics.velocities[pair.1 as usize];
                    let vrel_x = vx_a - vx_b;
                    let vrel_y = vy_a - vy_b;
                    let (n_x, n_y) = if pid_c == pair.0 {
                        (nx_norm, ny_norm)
                    } else {
                        (-nx_norm, -ny_norm)
                    };
                    let dot = vrel_x * n_x + vrel_y * n_y;
                    let entry = pair_verdicts.entry(pair).or_insert((0, 0));
                    if dot > DOT_FLOOR {
                        entry.0 += 1;
                    } else if dot < -DOT_FLOOR {
                        entry.1 += 1;
                    }
                }
            }
        }

        let mut current_divergent: HashSet<(u16, u16)> = HashSet::new();
        for (pair, (conv, div)) in &pair_verdicts {
            if div > conv {
                current_divergent.insert(*pair);
            }
        }

        let mut all_pairs: HashSet<(u16, u16)> =
            self.divergence_counts.keys().copied().collect();
        all_pairs.extend(current_divergent.iter().copied());
        for pair in all_pairs {
            if current_divergent.contains(&pair) {
                *self.divergence_counts.entry(pair).or_insert(0) += 1;
            } else {
                self.divergence_counts.insert(pair, 0);
            }
        }
    }

    pub fn get(&self, a: u16, b: u16) -> usize {
        let pair = if a < b { (a, b) } else { (b, a) };
        self.divergence_counts.get(&pair).copied().unwrap_or(0)
    }
}

/// Per-step rifting-split diagnostics. Returned by
/// [`apply_rifting_split`].
#[derive(Clone, Debug, Default)]
pub struct RiftingSplitStats {
    /// Number of split events that fired this step.
    pub splits_count: usize,
    /// Plate ids allocated for the newly-spawned rift basins
    /// (in fire order — sorted canonically by parent pair).
    pub new_plate_ids_created: Vec<u16>,
    /// Total cells whose `age` was reset to `0` via Path 3.B.
    pub age_zeroed_cells: usize,
}

/// Apply one step of the rifting split closure.
///
/// **Mutates** `plate_id` (rift-strip cells reassigned to new
/// pid), `age` (Path 3.B event-driven reset to `0` on the new
/// plate's cells), and `kinematics.velocities` (pushed with the
/// new plate's velocity = parent + perpendicular offset).
///
/// Returns immediately with `RiftingSplitStats::default()` when
/// `params.enabled == false`.
///
/// **Signature note (spec extension surfaced)**: `plate_type:
/// &PlateTypeField` is added to the user-spec signature. Needed
/// because the rift-strip filter requires continental cells of
/// plate `a` only (oceanic ridge spreading is handled by Track B
/// Path 3.A, not by rifting splits).
pub fn apply_rifting_split(
    plate_id: &mut PlateIdField,
    plate_type: &PlateTypeField,
    age: &mut Field2D,
    s: &Field2D,
    kinematics: &mut PlateKinematics,
    divergence_tracker: &DivergenceTracker,
    params: &RiftingParams,
) -> RiftingSplitStats {
    if !params.enabled {
        return RiftingSplitStats::default();
    }

    // Candidate pairs satisfying the time condition.
    let mut candidate_pairs: Vec<(u16, u16)> = divergence_tracker
        .divergence_counts
        .iter()
        .filter(|&(_pair, &count)| count >= params.split_time_threshold)
        .map(|(&pair, _)| pair)
        .collect();
    candidate_pairs.sort();

    let mut stats = RiftingSplitStats::default();

    for (a, b) in candidate_pairs {
        // Plate-id cap (W6) — refuse split gracefully when full.
        if kinematics.velocities.len() >= params.plate_id_cap {
            continue;
        }

        // Find the rift strip: continental cells of plate `a`
        // 4-adjacent to any cell of plate `b`.
        let rift_strip = find_rift_strip(plate_id, plate_type, a, b);
        if rift_strip.is_empty() {
            // No continental cells of `a` touch `b` — skip
            // (rift would have nothing to split).
            continue;
        }

        // Thickness condition: min S̃ across rift strip <
        // threshold.
        let min_s = rift_strip
            .iter()
            .map(|&(i, j)| s.get(i, j))
            .fold(f64::INFINITY, f64::min);
        if min_s >= params.split_thickness_threshold {
            continue;
        }

        // Both conditions met — fire split.
        let new_pid = kinematics.velocities.len() as u16;

        // Perpendicular velocity offset (Q3.4).
        let (vx_a, vy_a) = kinematics.velocities[a as usize];
        let (vx_b, vy_b) = kinematics.velocities[b as usize];
        let vrel_x = vx_a - vx_b;
        let vrel_y = vy_a - vy_b;
        let vrel_mag = (vrel_x * vrel_x + vrel_y * vrel_y).sqrt();
        // Right-hand-rule perpendicular: rotate v_rel by +90°.
        let (perp_x, perp_y) = if vrel_mag > DOT_FLOOR {
            (-vrel_y / vrel_mag, vrel_x / vrel_mag)
        } else {
            (1.0, 0.0)
        };
        let new_vx = vx_a + perp_x * params.split_velocity_offset;
        let new_vy = vy_a + perp_y * params.split_velocity_offset;
        kinematics.velocities.push((new_vx, new_vy));

        // Reassign strip cells + Path 3.B age = 0 on each.
        for (i, j) in &rift_strip {
            plate_id.set(*i, *j, new_pid);
            age.set(*i, *j, 0.0);
            stats.age_zeroed_cells += 1;
        }

        stats.splits_count += 1;
        stats.new_plate_ids_created.push(new_pid);
    }

    stats
}

/// Cells of plate `a` that are Continental AND 4-adjacent to at
/// least one cell of plate `b`. Used by [`apply_rifting_split`]
/// to delineate the strip that becomes the new rift plate.
fn find_rift_strip(
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
    a: u16,
    b: u16,
) -> Vec<(usize, usize)> {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut strip = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            if plate_id.get(i, j) != a {
                continue;
            }
            if plate_type.get(i, j) != PlateType::Continental {
                continue;
            }
            let neighbours = [
                (idx_x.next(i), j),
                (idx_x.prev(i), j),
                (i, idx_y.next(j)),
                (i, idx_y.prev(j)),
            ];
            for (ni, nj) in neighbours {
                if plate_id.get(ni, nj) == b {
                    strip.push((i, j));
                    break;
                }
            }
        }
    }
    strip
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::tectonics_c1::boundary_classification::classify_boundaries;
    use crate::tectonics_c1::closures::rifting::source_term::apply_rifting_thinning;

    /// Three-plate east-west fixture (mirrors the accretion +
    /// source_term tests' helper).
    fn three_plate_fixture(
        nx: usize,
        ny: usize,
        plate_types: [PlateType; 3],
        velocities: [(f64, f64); 3],
    ) -> (
        Field2D,
        Field2D,
        PlateIdField,
        PlateTypeField,
        PlateKinematics,
    ) {
        assert_eq!(nx % 3, 0);
        let mut s = Field2D::new(nx, ny);
        let mut age = Field2D::new(nx, ny);
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
                age.set(i, j, 5.0); // arbitrary non-zero baseline
            }
        }
        let kinematics = PlateKinematics { velocities: velocities.to_vec() };
        (s, age, plate_id, plate_type, kinematics)
    }

    #[test]
    fn rifting_split_requires_both_time_and_thickness() {
        // Scenario A — time condition met (counter >= 75), but
        // S̃ at rift strip is still above threshold (1.0 > 0.7).
        // → no split.
        let nx = 9;
        let ny = 4;
        let (s, mut age, mut plate_id, plate_type, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        let mut tracker = DivergenceTracker::new();
        tracker.divergence_counts.insert((0, 1), 75);
        let params = RiftingParams::default();
        let stats = apply_rifting_split(
            &mut plate_id,
            &plate_type,
            &mut age,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(
            stats.splits_count, 0,
            "no split — thickness condition not met"
        );

        // Scenario B — thickness condition met (S̃ = 0.5 at strip),
        // but counter below threshold. → no split.
        let nx = 9;
        let ny = 4;
        let (mut s2, mut age2, mut plate_id2, plate_type2, mut kinematics2) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        // Manually thin the strip below threshold.
        s2.set(2, 0, 0.5);
        s2.set(2, 1, 0.5);
        s2.set(2, 2, 0.5);
        s2.set(2, 3, 0.5);
        let mut tracker2 = DivergenceTracker::new();
        tracker2.divergence_counts.insert((0, 1), 10); // below 75
        let stats2 = apply_rifting_split(
            &mut plate_id2,
            &plate_type2,
            &mut age2,
            &s2,
            &mut kinematics2,
            &tracker2,
            &params,
        );
        assert_eq!(
            stats2.splits_count, 0,
            "no split — time condition not met"
        );
    }

    #[test]
    fn rifting_split_triggers_when_both_conditions_met() {
        // Time and thickness both satisfied → split fires.
        let nx = 9;
        let ny = 4;
        let (mut s, mut age, mut plate_id, plate_type, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        // Thin the rift strip on plate 0 (cells at i=2) below
        // threshold 0.7.
        for j in 0..ny {
            s.set(2, j, 0.5);
        }
        let mut tracker = DivergenceTracker::new();
        tracker.divergence_counts.insert((0, 1), 75);
        let params = RiftingParams::default();

        let plate_id_at_2_before: Vec<u16> = (0..ny).map(|j| plate_id.get(2, j)).collect();
        assert!(plate_id_at_2_before.iter().all(|&p| p == 0));

        let stats = apply_rifting_split(
            &mut plate_id,
            &plate_type,
            &mut age,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );

        assert_eq!(stats.splits_count, 1, "split should fire — both conditions met");
        assert_eq!(stats.new_plate_ids_created.len(), 1);
        let new_pid = stats.new_plate_ids_created[0];
        assert_eq!(new_pid, 3, "next pid after 3 base plates");
        // The rift strip cells (plate 0 cells touching plate 1)
        // are at i = 2 (4-adjacent to i = 3 plate 1). Verify
        // they are now the new pid.
        for j in 0..ny {
            assert_eq!(
                plate_id.get(2, j),
                new_pid,
                "rift strip cell (2, {j}) should be reassigned"
            );
        }
        // Plate 1 cells (i = 3..6) unchanged.
        for j in 0..ny {
            assert_eq!(plate_id.get(3, j), 1);
        }
        // kinematics.velocities extended with the new plate's velocity.
        assert_eq!(kinematics.velocities.len(), 4);
    }

    #[test]
    fn rifting_split_velocity_perpendicular_offset() {
        // Verify the perpendicular offset formula:
        //   v_rel = v_a - v_b
        //   perp = (-vrel_y / |vrel|, +vrel_x / |vrel|)  [right-hand-rule]
        //   v_new = v_a + perp · split_velocity_offset
        //
        // With v_a = (-0.01, 0), v_b = (0.01, 0):
        //   vrel = (-0.02, 0), |vrel| = 0.02
        //   perp = (0, -1)
        //   v_new = (-0.01 + 0 · 0.005, 0 + (-1) · 0.005)
        //         = (-0.01, -0.005)
        let nx = 9;
        let ny = 4;
        let (mut s, mut age, mut plate_id, plate_type, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        for j in 0..ny {
            s.set(2, j, 0.5);
        }
        let mut tracker = DivergenceTracker::new();
        tracker.divergence_counts.insert((0, 1), 75);
        let params = RiftingParams::default();

        let _stats = apply_rifting_split(
            &mut plate_id,
            &plate_type,
            &mut age,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );

        let (vx_new, vy_new) = kinematics.velocities[3];
        let expected_vx = -0.01;
        let expected_vy = -0.005;
        assert!(
            (vx_new - expected_vx).abs() < 1e-12,
            "perpendicular offset vx: got {vx_new}, expected {expected_vx}"
        );
        assert!(
            (vy_new - expected_vy).abs() < 1e-12,
            "perpendicular offset vy: got {vy_new}, expected {expected_vy}"
        );
    }

    #[test]
    fn rifting_split_age_zero_on_new_boundary() {
        // Path 3.B verification: every cell reassigned to the new
        // plate has its age reset to 0 (vs the baseline 5.0 set
        // in the fixture).
        let nx = 9;
        let ny = 4;
        let (mut s, mut age, mut plate_id, plate_type, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        for j in 0..ny {
            s.set(2, j, 0.5);
        }
        let mut tracker = DivergenceTracker::new();
        tracker.divergence_counts.insert((0, 1), 75);
        let params = RiftingParams::default();

        // Sanity: baseline age = 5.0 across the rift strip.
        for j in 0..ny {
            assert_eq!(age.get(2, j), 5.0);
        }

        let stats = apply_rifting_split(
            &mut plate_id,
            &plate_type,
            &mut age,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(stats.splits_count, 1);
        assert_eq!(stats.age_zeroed_cells, ny);

        for j in 0..ny {
            assert_eq!(age.get(2, j), 0.0, "Path 3.B age = 0 on rift cell (2, {j})");
        }
        // Non-rift cells retain baseline age 5.0.
        assert_eq!(age.get(0, 0), 5.0, "non-rift cell age must be unchanged");
        assert_eq!(age.get(5, 0), 5.0, "non-rift cell age must be unchanged");
    }

    #[test]
    fn rifting_disabled_no_op() {
        // Both apply_rifting_thinning and apply_rifting_split
        // disabled — bit-identical state across all mutable
        // surfaces.
        let nx = 9;
        let ny = 4;
        let (mut s, mut age, mut plate_id, plate_type, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [PlateType::Continental, PlateType::Continental, PlateType::Continental],
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        // Thin the strip below threshold + counter above
        // threshold → split WOULD fire if enabled.
        for j in 0..ny {
            s.set(2, j, 0.5);
        }
        let mut tracker = DivergenceTracker::new();
        tracker.divergence_counts.insert((0, 1), 100);
        let params = RiftingParams {
            enabled: false,
            ..RiftingParams::default()
        };
        let dt = 0.69;

        let s_before: Vec<f64> = s.data().to_vec();
        let age_before: Vec<f64> = age.data().to_vec();
        let plate_id_before: Vec<u16> = plate_id.data().to_vec();
        let kinematics_before = kinematics.velocities.clone();
        let boundary_info = classify_boundaries(&plate_id, &kinematics);

        let thinning_stats = apply_rifting_thinning(
            &mut s,
            &plate_type,
            &plate_id,
            &boundary_info,
            &kinematics,
            &params,
            dt,
        );
        let split_stats = apply_rifting_split(
            &mut plate_id,
            &plate_type,
            &mut age,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );

        assert_eq!(thinning_stats.cells_thinned, 0);
        assert_eq!(thinning_stats.total_mass_removed, 0.0);
        assert_eq!(split_stats.splits_count, 0);
        assert_eq!(split_stats.age_zeroed_cells, 0);
        assert!(split_stats.new_plate_ids_created.is_empty());

        assert_eq!(s.data(), s_before.as_slice());
        assert_eq!(age.data(), age_before.as_slice());
        assert_eq!(plate_id.data(), plate_id_before.as_slice());
        assert_eq!(kinematics.velocities, kinematics_before);
    }
}
