//! Sustained-convergence plate merge — `ConvergenceTracker` + the
//! `apply_accretion_step` event.
//!
//! ## ConvergenceTracker — sustained-convergence accumulator
//!
//! Per-pair counter keyed by canonical `(a, b)` with `a < b`.
//! [`ConvergenceTracker::update`] re-scans every grid cell each
//! step, aggregates a `(convergent_count, divergent_count)` per
//! adjacent pair, and updates the counter:
//!
//! - **Convergent verdict** (`convergent_count > divergent_count`):
//!   increment the pair's counter by 1.
//! - **Non-convergent verdict** (everything else — divergent,
//!   transform, no shared edges): reset the counter to 0.
//!
//! The reset on non-convergent is **load-bearing**: a plate pair
//! that briefly converges then drifts apart must not accumulate
//! merges over disjoint convergence windows. This is the explicit
//! anti-pattern called out in the Stage E2 spec.
//!
//! ## apply_accretion_step — merge event
//!
//! Once per step, after `ConvergenceTracker::update`. For each
//! pair `(a, b)` whose counter `≥ merge_time_threshold`:
//!
//! 1. Compute per-plate mass `m_a = Σ S̃ over plate-a cells`,
//!    same for `m_b` (W2 — sum of s_field over plate cells).
//! 2. **Lower-index plate wins** (`winner = a` by canonical
//!    ordering). Documented to be deterministic across runs.
//! 3. Compute mass-weighted average velocity (Q2.4) and assign
//!    it to `kinematics.velocities[winner]`.
//! 4. Reassign all `plate_id == loser` cells to `winner`.
//!    `plate_type` is preserved per cell (continental cells stay
//!    continental etc.).
//! 5. **Plate-id gaps are left in place** (W3 recommendation —
//!    compaction deferred to E4 if rifting needs lowest-available
//!    index). The loser id is simply unused after the merge;
//!    `count_distinct_plates(plate_id)` reflects the new total.
//!
//! ## Why the reset doesn't happen inside `apply_accretion_step`
//!
//! The user's E2 spec lists step (f) "Reset convergence counter
//! for merged pair", but the [`ConvergenceTracker`] takes
//! `&ConvergenceTracker` (immutable) by spec. The reset happens
//! **implicitly on the next `tracker.update()` call**: once the
//! loser plate id no longer appears in `plate_id`, no boundary
//! cells will contribute to the merged pair's verdict, so the
//! pair falls into the "no convergent verdict" branch and the
//! counter resets to 0. This keeps the API surface symmetric with
//! the rifting `DivergenceTracker` (Stage E3) and the per-step
//! lifecycle clean.
//!
//! ## Spec extension surfaced
//!
//! The Stage E2 spec listed `apply_accretion_step(plate_id,
//! kinematics, convergence_tracker, params)` — no `s` parameter.
//! But Watchpoint W2 says "Mass calculation per plate: sum of
//! s_field over plate cells". To honour W2 (mass-weighted average,
//! Q2.4), `s: &Field2D` is added to the signature. Read-only — the
//! closure does not mutate `S̃`.

use std::collections::{HashMap, HashSet};

use crate::tectonics_v2::field::{Field2D, PeriodicIndex};
use crate::tectonics_v2::voronoi::PlateIdField;

use crate::tectonics_c1::kinematics::PlateKinematics;

use super::params::{AccretionParams, VelocityMergeMethod};

/// Magnitude floor under which a dot product is treated as zero
/// (no contribution to either convergent or divergent count). Same
/// philosophy as `V_REL_FLOOR` in `boundary_classification.rs`;
/// guards against floating-point edge cases.
const DOT_FLOOR: f64 = 1e-12;

/// Per-pair sustained-convergence counter.
///
/// The map's keys are **canonical ordered** plate-id pairs
/// `(a, b)` with `a < b`. Each value is the number of consecutive
/// steps the pair has spent in a net-convergent state per
/// [`ConvergenceTracker::update`].
///
/// Lifetime: allocated by `run_with_closures` at run start
/// (Stage E4 wiring), dropped at run end. Not part of `C1State`
/// since it is not a save-restorable property of the world (per
/// Issue #132 Q-E1.2 Option (c)).
#[derive(Clone, Debug, Default)]
pub struct ConvergenceTracker {
    pub convergence_counts: HashMap<(u16, u16), usize>,
}

impl ConvergenceTracker {
    /// Empty tracker — no pairs known yet. Memory usage grows on
    /// `update` to `O(num_plate_pairs)` (bounded by `num_plates²`).
    pub fn new() -> Self {
        Self { convergence_counts: HashMap::new() }
    }

    /// Re-scan the grid and update per-pair counters.
    ///
    /// Pairs in a net-convergent state this step (more convergent
    /// edges than divergent across their shared boundary) have
    /// their counter incremented by 1. All other tracked pairs
    /// (formerly convergent, currently divergent / transform /
    /// disconnected) have their counter reset to 0.
    ///
    /// Per-edge classification uses the same outward-normal +
    /// `v_rel · n̂` convention as
    /// [`crate::tectonics_c1::boundary_classification::classify_boundaries`].
    pub fn update(&mut self, plate_id: &PlateIdField, kinematics: &PlateKinematics) {
        let nx = plate_id.nx();
        let ny = plate_id.ny();
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        // Per-pair (convergent_edges, divergent_edges) aggregate.
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
                    // Canonical perspective: compute v_rel = v_a − v_b
                    // and a normal pointing from a's side toward b's
                    // side. If pid_c is a (pair.0), the outward normal
                    // from c to n is from a to b — standard. If pid_c
                    // is b (pair.1), flip the normal sign.
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

        // Derive per-pair net verdict.
        let mut current_convergent: HashSet<(u16, u16)> = HashSet::new();
        for (pair, (conv, div)) in &pair_verdicts {
            if conv > div {
                current_convergent.insert(*pair);
            }
        }

        // Update counters: increment net-convergent pairs, reset
        // every other tracked pair to 0.
        let mut all_pairs: HashSet<(u16, u16)> =
            self.convergence_counts.keys().copied().collect();
        all_pairs.extend(current_convergent.iter().copied());
        for pair in all_pairs {
            if current_convergent.contains(&pair) {
                *self.convergence_counts.entry(pair).or_insert(0) += 1;
            } else {
                self.convergence_counts.insert(pair, 0);
            }
        }
    }

    /// Convenience getter — returns the current count for a pair,
    /// canonicalising the argument order. `0` for unknown pairs.
    pub fn get(&self, a: u16, b: u16) -> usize {
        let pair = if a < b { (a, b) } else { (b, a) };
        self.convergence_counts.get(&pair).copied().unwrap_or(0)
    }
}

/// Per-step accretion diagnostics. Returned by
/// [`apply_accretion_step`].
#[derive(Clone, Copy, Debug, Default)]
pub struct AccretionStats {
    /// Number of plate-pair merges that fired this step.
    pub merges_count: usize,
    /// Number of distinct plate ids remaining in `plate_id` after
    /// this step's merges (post-merge cell count by unique id).
    /// `0` when the closure is disabled — placeholder, not the
    /// true plate count.
    pub plates_remaining: usize,
}

/// Apply one step of the accretion closure.
///
/// **Mutates** `plate_id` (loser cells reassigned to winner) and
/// `kinematics.velocities[winner]` (mass-weighted average). The
/// loser-id slot of `kinematics.velocities` is left unchanged
/// (gap, per W3 — no compaction in E2).
///
/// `s`, `convergence_tracker`, and `params` are read-only inputs.
///
/// Returns immediately with `AccretionStats::default()` (zero
/// merges, zero plates_remaining) when `params.enabled == false`.
pub fn apply_accretion_step(
    plate_id: &mut PlateIdField,
    s: &Field2D,
    kinematics: &mut PlateKinematics,
    convergence_tracker: &ConvergenceTracker,
    params: &AccretionParams,
) -> AccretionStats {
    if !params.enabled {
        return AccretionStats::default();
    }

    // Find pairs ready to merge — counter ≥ threshold.
    let mut merging_pairs: Vec<(u16, u16)> = convergence_tracker
        .convergence_counts
        .iter()
        .filter(|&(_pair, &count)| count >= params.merge_time_threshold)
        .map(|(&pair, _)| pair)
        .collect();
    // Deterministic processing order — sorted by canonical pair.
    merging_pairs.sort();

    let mut stats = AccretionStats::default();

    for (a, b) in merging_pairs {
        let mass_a = compute_plate_mass(plate_id, s, a);
        let mass_b = compute_plate_mass(plate_id, s, b);
        // Skip if either side has been absorbed by an earlier
        // merge in this step (mass == 0 means no cells with that
        // plate id remain).
        if mass_a == 0.0 || mass_b == 0.0 {
            continue;
        }
        let total_mass = mass_a + mass_b;
        if total_mass < DOT_FLOOR {
            // Degenerate — both depleted. Skip rather than divide
            // by zero.
            continue;
        }

        let winner = a;
        let loser = b;

        // Mass-weighted velocity blend per
        // VelocityMergeMethod::MassWeightedAverage (Q2.4).
        let (vx_a, vy_a) = kinematics.velocities[a as usize];
        let (vx_b, vy_b) = kinematics.velocities[b as usize];
        let (avg_vx, avg_vy) = match params.velocity_merge_method {
            VelocityMergeMethod::MassWeightedAverage => (
                (vx_a * mass_a + vx_b * mass_b) / total_mass,
                (vy_a * mass_a + vy_b * mass_b) / total_mass,
            ),
        };
        kinematics.velocities[winner as usize] = (avg_vx, avg_vy);

        // Reassign loser cells to winner. plate_type unchanged.
        for j in 0..plate_id.ny() {
            for i in 0..plate_id.nx() {
                if plate_id.get(i, j) == loser {
                    plate_id.set(i, j, winner);
                }
            }
        }

        stats.merges_count += 1;
    }

    stats.plates_remaining = count_distinct_plates(plate_id);
    stats
}

/// Sum of `s_field` over cells with `plate_id == p`. Returns `0`
/// when no cell carries id `p` (e.g., after the plate has been
/// absorbed by an earlier merge in this step).
fn compute_plate_mass(plate_id: &PlateIdField, s: &Field2D, p: u16) -> f64 {
    let mut total = 0.0_f64;
    for j in 0..plate_id.ny() {
        for i in 0..plate_id.nx() {
            if plate_id.get(i, j) == p {
                total += s.get(i, j);
            }
        }
    }
    total
}

/// Number of distinct plate ids present in `plate_id`.
fn count_distinct_plates(plate_id: &PlateIdField) -> usize {
    let mut seen: HashSet<u16> = HashSet::new();
    for &pid in plate_id.data() {
        seen.insert(pid);
    }
    seen.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::tectonics_v2::boundaries::plate_type::PlateType;
    use crate::tectonics_v2::boundaries::plate_type::PlateTypeField;

    /// Two-plate east-west fixture: plate 0 in left half, plate 1
    /// in right half. Both Continental, S̃ uniform.
    ///
    /// **Caveat — periodic-wrap symmetry.** With only two plates on
    /// a torus, the interior boundary and the wrap-around boundary
    /// share the same pair `(0, 1)` and produce opposite
    /// convergent-vs-divergent verdicts. Counts tie (`conv == div`)
    /// and the `conv > div` net-verdict in
    /// [`ConvergenceTracker::update`] returns "not convergent",
    /// regardless of how `v_left` / `v_right` are oriented. Use
    /// [`three_plate_fixture`] for tests that drive
    /// `tracker.update`. Use this two-plate variant when the
    /// counter is pre-populated manually (the tracker is bypassed).
    fn two_plate_fixture(
        nx: usize,
        ny: usize,
        v_left: (f64, f64),
        v_right: (f64, f64),
    ) -> (Field2D, PlateIdField, PlateTypeField, PlateKinematics) {
        let mut s = Field2D::new(nx, ny);
        let mut plate_id = PlateIdField::new(nx, ny);
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Continental);
        for j in 0..ny {
            for i in 0..nx {
                let pid = if i < nx / 2 { 0_u16 } else { 1_u16 };
                plate_id.set(i, j, pid);
                s.set(i, j, 1.0);
            }
        }
        let kinematics = PlateKinematics { velocities: vec![v_left, v_right] };
        (s, plate_id, plate_type, kinematics)
    }

    /// Three-plate east-west fixture: thirds of the grid carry
    /// `plate_id ∈ {0, 1, 2}`. `nx` must be divisible by 3.
    ///
    /// Used by tests that need `tracker.update` to produce a clear
    /// per-pair verdict — pair `(0, 1)` appears only at the
    /// interior boundary `i = nx/3 - 1 / nx/3` (no periodic wrap
    /// symmetry), so its verdict is purely the interior
    /// convergent/divergent character.
    fn three_plate_fixture(
        nx: usize,
        ny: usize,
        velocities: [(f64, f64); 3],
    ) -> (Field2D, PlateIdField, PlateKinematics) {
        assert_eq!(nx % 3, 0, "three_plate_fixture: nx must be divisible by 3");
        let mut s = Field2D::new(nx, ny);
        let mut plate_id = PlateIdField::new(nx, ny);
        let third = nx / 3;
        for j in 0..ny {
            for i in 0..nx {
                let pid = if i < third {
                    0_u16
                } else if i < 2 * third {
                    1_u16
                } else {
                    2_u16
                };
                plate_id.set(i, j, pid);
                s.set(i, j, 1.0);
            }
        }
        let kinematics = PlateKinematics { velocities: velocities.to_vec() };
        (s, plate_id, kinematics)
    }

    #[test]
    fn accretion_does_not_merge_below_time_threshold() {
        // Convergent pair, but counter at only 1 step (well below
        // default threshold = 50). No merge.
        //
        // Uses the three-plate fixture so pair (0, 1) appears only
        // at the interior boundary (i = 2/3), avoiding the
        // 2-plate-on-torus wrap symmetry that ties conv ↔ div
        // counts to zero net verdict.
        let nx = 9;
        let ny = 4;
        let (s, mut plate_id, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [(0.01, 0.0), (-0.01, 0.0), (0.0, 0.0)],
        );
        let mut tracker = ConvergenceTracker::new();
        tracker.update(&plate_id, &kinematics);
        // Pair (0, 1) is registered as convergent on the interior
        // boundary — counter incremented to 1.
        assert_eq!(
            tracker.get(0, 1),
            1,
            "after one update the convergent (0, 1) pair should have count = 1"
        );
        let params = AccretionParams::default(); // threshold = 50
        let stats = apply_accretion_step(
            &mut plate_id,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(stats.merges_count, 0, "no merge expected below threshold");
        assert_eq!(stats.plates_remaining, 3);
        // Plate 1 cells must still be plate 1.
        assert_eq!(plate_id.get(3, 0), 1);
    }

    #[test]
    fn accretion_merges_at_time_threshold() {
        // Pre-populate the counter at the threshold and verify
        // the merge fires.
        let nx = 8;
        let ny = 4;
        let (s, mut plate_id, _plate_type, mut kinematics) =
            two_plate_fixture(nx, ny, (0.01, 0.0), (-0.01, 0.0));
        let mut tracker = ConvergenceTracker::new();
        tracker.convergence_counts.insert((0, 1), 50);
        let params = AccretionParams::default();
        let stats = apply_accretion_step(
            &mut plate_id,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(stats.merges_count, 1);
        assert_eq!(stats.plates_remaining, 1);
        // Every cell now belongs to the winner plate 0.
        for j in 0..ny {
            for i in 0..nx {
                assert_eq!(plate_id.get(i, j), 0, "cell ({i}, {j}) should be plate 0");
            }
        }
    }

    #[test]
    fn accretion_velocity_mass_weighted_average() {
        // Plate 0 has 4 columns × 4 rows = 16 cells × S̃ = 1.0 → mass 16.
        // Plate 1 has 4 columns × 4 rows = 16 cells × S̃ = 1.0 → mass 16.
        // Equal mass → simple arithmetic mean of velocities.
        let nx = 8;
        let ny = 4;
        let v_a = (0.01, 0.0);
        let v_b = (-0.01, 0.005);
        let (s, mut plate_id, _plate_type, mut kinematics) =
            two_plate_fixture(nx, ny, v_a, v_b);
        let mut tracker = ConvergenceTracker::new();
        tracker.convergence_counts.insert((0, 1), 50);
        let params = AccretionParams::default();
        let _stats = apply_accretion_step(
            &mut plate_id,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        let expected_vx = (v_a.0 + v_b.0) / 2.0;
        let expected_vy = (v_a.1 + v_b.1) / 2.0;
        let (got_vx, got_vy) = kinematics.velocities[0];
        assert!(
            (got_vx - expected_vx).abs() < 1e-12,
            "post-merge vx {got_vx} != mass-weighted expected {expected_vx}"
        );
        assert!(
            (got_vy - expected_vy).abs() < 1e-12,
            "post-merge vy {got_vy} != mass-weighted expected {expected_vy}"
        );
    }

    #[test]
    fn accretion_resets_counter_on_non_convergent() {
        // Pre-populate counter at threshold, then drive a tracker
        // update with DIVERGENT velocities on pair (0, 1). The
        // counter must reset to 0 and the subsequent
        // apply_accretion_step must NOT merge.
        //
        // Three-plate fixture: pair (0, 1) at interior i = 2/3
        // has clear divergent verdict (no wrap symmetry).
        let nx = 9;
        let ny = 4;
        let (s, mut plate_id, mut kinematics) = three_plate_fixture(
            nx,
            ny,
            [(-0.01, 0.0), (0.01, 0.0), (0.0, 0.0)],
        );
        let mut tracker = ConvergenceTracker::new();
        tracker.convergence_counts.insert((0, 1), 50);
        // Tracker update on the divergent kinematics — must reset.
        tracker.update(&plate_id, &kinematics);
        assert_eq!(
            tracker.get(0, 1),
            0,
            "counter for diverging pair must reset to 0"
        );
        let params = AccretionParams::default();
        let stats = apply_accretion_step(
            &mut plate_id,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(stats.merges_count, 0, "no merge — counter was reset");
        assert_eq!(stats.plates_remaining, 3);
    }

    #[test]
    fn accretion_disabled_no_op() {
        // disabled = true, counter above threshold — no merge,
        // bit-identical state.
        let nx = 8;
        let ny = 4;
        let (s, mut plate_id, _plate_type, mut kinematics) =
            two_plate_fixture(nx, ny, (0.01, 0.0), (-0.01, 0.0));
        let mut tracker = ConvergenceTracker::new();
        tracker.convergence_counts.insert((0, 1), 100);
        let params = AccretionParams {
            enabled: false,
            ..AccretionParams::default()
        };
        let plate_id_before: Vec<u16> = plate_id.data().to_vec();
        let kinematics_before = kinematics.velocities.clone();
        let stats = apply_accretion_step(
            &mut plate_id,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(stats.merges_count, 0);
        assert_eq!(stats.plates_remaining, 0); // default-disabled placeholder
        assert_eq!(plate_id.data(), plate_id_before.as_slice());
        assert_eq!(kinematics.velocities, kinematics_before);
    }

    #[test]
    fn accretion_plate_id_compaction_after_merge() {
        // 3 plates. Merge (0, 1) but plate 2 remains independent.
        // Verify: no compaction — surviving ids are {0, 2} (1 is
        // a gap), `plates_remaining == 2`, kinematics.velocities
        // still has 3 entries (gap at index 1).
        let nx = 9;
        let ny = 4;
        let mut s = Field2D::new(nx, ny);
        let mut plate_id = PlateIdField::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let pid = if i < 3 {
                    0
                } else if i < 6 {
                    1
                } else {
                    2
                };
                plate_id.set(i, j, pid);
                s.set(i, j, 1.0);
            }
        }
        let mut kinematics = PlateKinematics {
            velocities: vec![(0.01, 0.0), (-0.01, 0.0), (0.0, 0.01)],
        };
        let mut tracker = ConvergenceTracker::new();
        tracker.convergence_counts.insert((0, 1), 50);
        let params = AccretionParams::default();
        let stats = apply_accretion_step(
            &mut plate_id,
            &s,
            &mut kinematics,
            &tracker,
            &params,
        );
        assert_eq!(stats.merges_count, 1);
        assert_eq!(stats.plates_remaining, 2);
        // Surviving ids: {0, 2}. Plate 1 is a gap.
        let mut seen: HashSet<u16> = HashSet::new();
        for &pid in plate_id.data() {
            seen.insert(pid);
        }
        assert_eq!(seen, HashSet::from([0_u16, 2_u16]));
        // Kinematics vec length unchanged (no compaction).
        assert_eq!(
            kinematics.velocities.len(),
            3,
            "kinematics.velocities should retain the gap entry at index 1 (W3 no-compaction)"
        );
    }
}
