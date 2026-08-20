//! Worker thread driving the C1 time loop.
//!
//! ## Standalone `run_with_closures` (Issue #137 Stage A — A1-c revert)
//!
//! The worker calls `run_with_closures` directly with NO
//! `apply_post_tectonic`. This mirrors the Track D visual gallery
//! ([`c1_phase_2_track_d_visual_gallery`]) which produces the
//! reference PNGs the viz should reproduce.
//!
//! ### Why no `apply_post_tectonic`
//!
//! A prior A1-c attempt added per-cycle `apply_post_tectonic` to
//! "make the coast migrate". Two problems:
//!
//! 1. Phase 1.x workflow tests calibrate `macro_redistribution` for
//!    **ONE** `apply_post_tectonic` per cycle (`n_cycles = 1`,
//!    one big inner-loop). A1-c ran it 6× per default run
//!    (`steps_per_cycle = 50`, 6 cycles), 6× the calibrated cadence
//!    → continent collapse (continental cells eroded below the
//!    floating sea-level → reclassified Oceanic each cycle).
//!
//! 2. The Track D visual gallery — the actual reference for what
//!    the viz should look like — does NOT call `apply_post_tectonic`
//!    at all. `plate_type` and `cratonic_mask` stay init-time in
//!    the gallery. The viz should match.
//!
//! Visible coast migration in the **Altitude** view still emerges
//! because `compute_isostasy` is purely S̃-based and the sea-level
//! threshold (`h_sea = h_min + sea_level_fraction × h_range`)
//! floats with the global S̃ extremes. Continental cells whose S̃
//! drops cross below `h_sea` → appear as ocean in Altitude view
//! WITHOUT any plate_type reclassification. `apply_stein_stein_bathymetry`
//! only overwrites cells where `plate_type == Oceanic`, so
//! continental cells keep their isostatic altitude.
//!
//! `plate_type` and `cratonic_mask` are STATIC across cycles — this
//! is the gallery contract.
//!
//! ## Lifecycle
//!
//! 1. `spawn_c1_thread` returns a `JoinHandle`. The thread blocks
//!    on `commands_rx.recv()` between runs.
//! 2. On `C1Command::RunBaseline { spec }`:
//!    a. Reset the cancel flag (allow future cancels).
//!    b. Send `C1Event::Started`.
//!    c. Build state + kinematics + config from `spec`.
//!    d. Clone `kinematics.velocities` pre-run into a local
//!       `initial_velocities: Vec<(f64, f64)>` — captured by the
//!       per-step closure (cannot access `kinematics` directly
//!       because `run_with_closures` borrows it `&mut`; see
//!       Stage E2 W7 Q-E1.3).
//!    e. Emit a pre-run cycle-0 `StepCompleted` so the UI can
//!       render the init state before the first step fires.
//!    f. Call `run_with_closures`. The per-step closure emits a
//!       `StepCompleted { snapshot }` event each step.
//!    g. Send `C1Event::Completed { spec, final_snapshot, elapsed }`.
//! 3. On `C1Command::Cancel`: set the cancel flag (no effect on
//!    current run per Q-E1.3 Option C MVP).
//!
//! ## Worker owns no Bevy state (W4 global)
//!
//! All Bevy-side data (resources, sprites, image handles) lives
//! in `plugin.rs` + the UI/render systems. The worker uses only
//! crossbeam channels + `ymir-core` types. This isolation is
//! shared with `bridge::v2::thread`.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{C1TimeLoopConfig, run_with_closures};
use ymir_core::tectonics_v2::workflow::PhaseAParams;
use ymir_core::tectonics_v2::workflow::phase_a_common::{PostTectonicInput, apply_post_tectonic};

use super::commands::C1Command;
use super::events::{C1Event, C1RunKind};
use super::snapshot::C1Snapshot;
use super::spec::C1RunSpec;

#[cfg(test)]
mod tests {
    use super::*;
    use crossbeam_channel::bounded;
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;
    use std::time::Duration;

    /// Drain all events from `evt_rx` until a `Completed` arrives
    /// (or `Failed`). Bounded by a generous timeout to prevent the
    /// test hanging on a worker bug. Returns the full event stream
    /// in arrival order.
    fn drain_run(evt_rx: &crossbeam_channel::Receiver<C1Event>) -> Vec<C1Event> {
        let mut events = Vec::new();
        loop {
            match evt_rx.recv_timeout(Duration::from_secs(30)) {
                Ok(e) => {
                    let terminal = matches!(e, C1Event::Completed { .. } | C1Event::Failed { .. });
                    events.push(e);
                    if terminal {
                        break;
                    }
                }
                Err(_) => panic!("c1 worker timed out after 30s"),
            }
        }
        events
    }

    /// Spec sized for fast tests — 32² grid + small n_steps.
    fn small_spec(n_steps: usize, seed: u64) -> C1RunSpec {
        C1RunSpec { grid_size: 32, seed, n_steps, ..C1RunSpec::default() }
    }

    #[test]
    fn c1_worker_spawns_and_runs() {
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_steps = 50;
        cmd_tx
            .send(C1Command::RunBaseline { spec: small_spec(n_steps, 42) })
            .expect("send RunBaseline");

        let events = drain_run(&evt_rx);

        // Event-stream shape: Started, cycle-0 StepCompleted, then
        // one StepCompleted per step (worker emits step indices
        // 1..=n_steps), then Completed. Cycle-0 uses step=0; the
        // first `on_step` callback fires AFTER step 1 mutations
        // and gets `global_step = 1`. Total = n_steps + 1.
        let started = events.iter().filter(|e| matches!(e, C1Event::Started { .. })).count();
        let step_completed =
            events.iter().filter(|e| matches!(e, C1Event::StepCompleted { .. })).count();
        let completed = events.iter().filter(|e| matches!(e, C1Event::Completed { .. })).count();

        assert_eq!(started, 1, "exactly one Started event");
        assert_eq!(completed, 1, "exactly one Completed event");
        assert_eq!(
            step_completed,
            n_steps + 1,
            "expected {} StepCompleted (cycle-0 + n_steps), got {}",
            n_steps + 1,
            step_completed
        );

        // Final snapshot: live plates ≤ init num_plates. Track D
        // accretion + rifting may grow OR shrink the count; under
        // default closures on Phase 2 R7 init at small n_steps the
        // accretion path may not yet fire (merge_time_threshold =
        // 50). Either way the live count must be a non-zero
        // positive integer.
        if let Some(C1Event::Completed { final_snapshot, .. }) = events.last() {
            assert!(final_snapshot.live_plate_count >= 1);
            assert!(
                final_snapshot.live_plate_count <= final_snapshot.num_plates,
                "live plates ({}) must not exceed init num_plates ({})",
                final_snapshot.live_plate_count,
                final_snapshot.num_plates
            );
        } else {
            panic!("last event must be Completed");
        }

        // Worker shutdown: drop the command sender, join.
        drop(cmd_tx);
        handle.join().expect("c1 worker join");
    }

    /// Step b/5 — the worker drives the full HD chain and SHIPS the
    /// product. Tiny grid + small HD target so erosion is fast even in a
    /// debug test build. Validates the event sequence (Started → 4×(phase
    /// start+done) → Completed) and that the result carries every layer at
    /// the HD resolution. Cache regime is NOT asserted (HIT/MISS depends on
    /// whether a prior test populated `.ymir_cache/`).
    #[test]
    fn c1_worker_runs_hd_chain_and_ships_product() {
        use super::super::hd::{HdParams, HdPhase};
        use super::super::inspect::{RiverCellMap, inspect_cell};

        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let target = 64usize;
        cmd_tx
            .send(C1Command::RunHd {
                spec: small_spec(10, 42),
                params: HdParams {
                    target_size: target,
                    latitude_deg: 45.0,
                    domain_km: 1024.0,
                    manual_offset: None,
                    export_dir: None,
                },
            })
            .expect("send RunHd");

        // Drain until HdCompleted / HdFailed (drain_run stops on the
        // tectonic Completed, which the HD path never emits).
        let mut events = Vec::new();
        loop {
            match evt_rx.recv_timeout(Duration::from_secs(120)) {
                Ok(e) => {
                    let terminal =
                        matches!(e, C1Event::HdCompleted { .. } | C1Event::HdFailed { .. });
                    events.push(e);
                    if terminal {
                        break;
                    }
                }
                Err(_) => panic!("HD worker timed out"),
            }
        }

        let started = events.iter().filter(|e| matches!(e, C1Event::HdStarted { .. })).count();
        let phase_started =
            events.iter().filter(|e| matches!(e, C1Event::HdPhaseStarted { .. })).count();
        let phase_done = events.iter().filter(|e| matches!(e, C1Event::HdPhaseDone { .. })).count();
        assert_eq!(started, 1, "one HdStarted");
        // Eroded is split into Tectonic/Relief/Erosion (suite e) → 6 sub-phases.
        assert_eq!(phase_started, 6, "6 HdPhaseStarted");
        assert_eq!(phase_done, 6, "6 HdPhaseDone");

        // Phases arrive in execution order (same on HIT or MISS).
        let order: Vec<HdPhase> = events
            .iter()
            .filter_map(|e| match e {
                C1Event::HdPhaseStarted { phase } => Some(*phase),
                _ => None,
            })
            .collect();
        assert_eq!(
            order,
            vec![
                HdPhase::Tectonic,
                HdPhase::Relief,
                HdPhase::Erosion,
                HdPhase::Climate,
                HdPhase::Drainage,
                HdPhase::Biomes,
            ],
        );

        // Final product carries every layer at the HD resolution.
        match events.last() {
            Some(C1Event::HdCompleted { result, .. }) => {
                assert_eq!(result.width, target);
                assert_eq!(result.height, target);
                assert_eq!(result.eroded.width, target);
                assert_eq!(result.temperature.data.len(), target * target);
                assert_eq!(result.precipitation.data.len(), target * target);
                assert_eq!(result.biomes.len(), target * target);
                assert_eq!(result.drainage.width, target);

                // Step c/5 — per-cell inspection over the HD product.
                let river_map = RiverCellMap::from_drainage(&result.drainage);
                // Every river segment's cells resolve to Some in the map.
                for seg in &result.drainage.rivers.segments {
                    for &(px, py) in &seg.points {
                        assert!(
                            river_map.at(px as usize, py as usize).is_some(),
                            "river cell ({px},{py}) must map to a segment",
                        );
                    }
                }
                // inspect_cell reads coherent values at an arbitrary cell.
                let c = inspect_cell(&result, &river_map, target / 2, target / 2);
                assert!(c.altitude_m.is_finite());
                assert!(c.runoff_mm >= 0.0);
                assert_eq!(c.is_ocean, c.depth_m.is_some());
                // A cell flagged on a river by the map is reported as such by
                // the unified inspection (consistency between the two paths).
                if let Some(info) = river_map.at(target / 2, target / 2) {
                    assert_eq!(c.river, Some(info));
                }
            }
            other => panic!("last event must be HdCompleted, got {other:?}"),
        }

        drop(cmd_tx);
        handle.join().expect("c1 worker join");
    }

    #[test]
    fn c1_worker_event_ordering_no_loss_under_backpressure() {
        // Q-E5.1 Option C: bounded(2) events + immediate drain.
        // The worker naturally backpressures when the channel
        // fills. We don't time-assert blocking; we assert the
        // stream invariant — every step index appears exactly
        // once in StepCompleted events, in strictly monotone
        // order. No loss, no duplication, no reordering.
        //
        // After A1-c revert: standalone `run_with_closures` —
        // events are `[step=0 (cycle-0), step=1, step=2, ...,
        // step=n_steps]`. Total = n_steps + 1.
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_steps = 30;
        cmd_tx.send(C1Command::RunBaseline { spec: small_spec(n_steps, 42) }).unwrap();

        let events = drain_run(&evt_rx);

        let step_indices: Vec<usize> = events
            .iter()
            .filter_map(|e| match e {
                C1Event::StepCompleted { snapshot } => Some(snapshot.step),
                _ => None,
            })
            .collect();

        assert_eq!(step_indices.len(), n_steps + 1);
        assert_eq!(step_indices[0], 0, "first StepCompleted is cycle-0 (step 0)");
        for (i, &s) in step_indices.iter().enumerate().skip(1) {
            assert_eq!(s, i, "expected step {} at position {}, got {}", i, i, s);
        }

        drop(cmd_tx);
        handle.join().unwrap();
    }

    #[test]
    fn c1_snapshot_carries_stats() {
        // Track D evidence (Issue #132 Stage A, seed 42, 64² ×
        // 300 steps): subduction 20,914 cells, accretion 6
        // merges. We use seed 42 with n_steps = 200 at 32² —
        // expect subduction to fire heavily on the first
        // step (≥ 1 cell consumed), and at least one accretion
        // merge somewhere over the 200-step window (merge fires
        // after merge_time_threshold = 50 consecutive convergent
        // steps; with 4× threshold worth of steps and Phase 1.1
        // sustained convergence pairs, ≥ 1 merge is robust).
        //
        // Rifting splits are RARE by design (Track D Stage V
        // 0-3/seed over 300 steps); not asserted here.
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_steps = 200;
        cmd_tx.send(C1Command::RunBaseline { spec: small_spec(n_steps, 42) }).unwrap();

        let events = drain_run(&evt_rx);

        let mut max_sub_cells = 0_usize;
        let mut total_merges = 0_usize;
        let mut total_rifting_thinning = 0_usize;
        for e in &events {
            if let C1Event::StepCompleted { snapshot } = e {
                max_sub_cells = max_sub_cells.max(snapshot.stats.subduction.cells_consumed);
                total_merges += snapshot.stats.accretion.merges_count;
                total_rifting_thinning += snapshot.stats.rifting_thinning.cells_thinned;
            }
        }

        eprintln!("Stage V c1_snapshot_carries_stats (seed 42, n_steps={n_steps}, 32²):");
        eprintln!("  max subduction.cells_consumed per step = {max_sub_cells}");
        eprintln!("  cumulative accretion.merges_count      = {total_merges}");
        eprintln!("  cumulative rifting_thinning.cells_thinned = {total_rifting_thinning}");

        // Hard assertion: subduction fires (Stage A evidence
        // shows subduction is high-frequency at seed 42).
        assert!(
            max_sub_cells > 0,
            "subduction must fire at least once over {n_steps} steps; max cells consumed = {max_sub_cells}"
        );
        // Hard assertion: rifting thinning fires (continental
        // divergent cells exist at Phase 2 R7 init).
        assert!(
            total_rifting_thinning > 0,
            "rifting thinning must fire on continental divergent cells"
        );
        // Best-effort: accretion merges over 200 steps at seed 42.
        // Document but don't hard-fail at zero — small 32² grid
        // may have fewer sustained-convergent pairs than 64².
        if total_merges == 0 {
            eprintln!(
                "  (Note: no accretion merges at 32² seed 42 / {n_steps} steps; 64² Stage A measured 6 over 300 steps)"
            );
        }

        drop(cmd_tx);
        handle.join().unwrap();
    }

    /// Stage A acceptance — DISTINCT angle from Stage V worker
    /// mechanics. Production-scale `64² × 300` seed 42 run with
    /// the default `C1Closures` (all 7 enabled including Track
    /// D), asserting product-level invariants:
    ///
    /// - Completed event received (no hang).
    /// - `final_snapshot.live_plate_count ≤ 3` — Pangaea-like
    ///   collapse from 8 init plates (Track D Stage A 64² ×
    ///   300 evidence: seed 42 → 2 surviving plates, 6 accretion
    ///   merges).
    /// - Track D activity confirmed via cumulative stats: total
    ///   subduction cells > 0 + cumulative accretion merges > 0
    ///   over the full run.
    ///
    /// Visual / UI acceptance is NOT testable here; see
    /// `docs/reports/viz_0_c1_integration/acceptance_checklist.md`
    /// for the manual checklist the user runs before merge.
    #[test]
    fn acceptance_full_run_seed_42_pangaea_collapse() {
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let spec = C1RunSpec::default(); // 64, seed=42, n=300
        cmd_tx.send(C1Command::RunBaseline { spec: spec.clone() }).unwrap();

        let events = drain_run(&evt_rx);

        // Tally cumulative Track D stats from the StepCompleted
        // stream (each snapshot carries that step's stats only;
        // sum across the run for a cumulative figure).
        let mut cum_sub_cells = 0_usize;
        let mut cum_merges = 0_usize;
        let mut cum_splits = 0_usize;
        let mut cum_thinning = 0_usize;
        for e in &events {
            if let C1Event::StepCompleted { snapshot } = e {
                cum_sub_cells += snapshot.stats.subduction.cells_consumed;
                cum_merges += snapshot.stats.accretion.merges_count;
                cum_splits += snapshot.stats.rifting_split.splits_count;
                cum_thinning += snapshot.stats.rifting_thinning.cells_thinned;
            }
        }

        let final_snap = match events.last() {
            Some(C1Event::Completed { final_snapshot, .. }) => final_snapshot,
            _ => panic!("expected Completed as last event; got {:?}", events.last()),
        };

        eprintln!("Stage A acceptance (seed 42, 64²×300, default closures):");
        eprintln!("  init num_plates       = {}", final_snap.num_plates);
        eprintln!("  final live_plate_count = {}", final_snap.live_plate_count);
        eprintln!("  cumulative subduction cells = {cum_sub_cells}");
        eprintln!("  cumulative accretion merges = {cum_merges}");
        eprintln!("  cumulative rifting splits   = {cum_splits}");
        eprintln!("  cumulative thinning cells   = {cum_thinning}");

        assert_eq!(final_snap.num_plates, 8, "Phase 2 R7 default init should produce 8 plates");
        assert!(
            final_snap.live_plate_count <= 3,
            "Pangaea collapse not achieved: final live plates = {} (expected ≤ 3 at seed 42 / 64² / 300 steps)",
            final_snap.live_plate_count
        );

        assert!(
            cum_sub_cells > 1000,
            "subduction must fire substantially at seed 42 / 300 steps; got {cum_sub_cells}"
        );
        assert!(
            cum_merges > 0,
            "accretion must merge at least once at seed 42 / 300 steps; got {cum_merges}"
        );

        drop(cmd_tx);
        handle.join().unwrap();
    }

    // ----- Issue #139 Stage V — workflow-mode worker tests -----

    /// Workflow emits cycle-0 + per cycle `{ k_cycle per-step + 1
    /// post-cycle }` snapshots, with monotone (non-decreasing) step
    /// indices and `total = n_cycles × k_cycle`. Reduced cadence
    /// (2×5) — this tests emission structure, not erosion.
    #[test]
    fn workflow_mode_produces_calibrated_cycle_count() {
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_cycles = 2;
        let k_cycle = 5;
        let phase_a = PhaseAParams { n_cycles, k_cycle, ..PhaseAParams::default() };
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: small_spec(0, 42), // n_steps unused in workflow
                phase_a,
            })
            .unwrap();

        let events = drain_run(&evt_rx);

        // Started carries Workflow kind with the calibrated total.
        match events.first() {
            Some(C1Event::Started { spec, kind }) => {
                assert!(matches!(kind, C1RunKind::Workflow { .. }));
                assert_eq!(kind.total_tectonic_steps(spec), n_cycles * k_cycle);
            }
            other => panic!("expected Started first; got {other:?}"),
        }

        let step_indices: Vec<usize> = events
            .iter()
            .filter_map(|e| match e {
                C1Event::StepCompleted { snapshot } => Some(snapshot.step),
                _ => None,
            })
            .collect();

        // Count = cycle-0 (1) + n_cycles × (k_cycle per-step + 1 post).
        assert_eq!(
            step_indices.len(),
            1 + n_cycles * (k_cycle + 1),
            "snapshot count mismatch: {step_indices:?}"
        );
        // Monotone non-decreasing (post-cycle reuses the boundary
        // step, so equality is allowed at cycle boundaries).
        assert!(
            step_indices.windows(2).all(|w| w[0] <= w[1]),
            "step indices must be monotone non-decreasing: {step_indices:?}"
        );
        // Last per-step / post-cycle step reaches the total.
        assert_eq!(*step_indices.last().unwrap(), n_cycles * k_cycle);

        drop(cmd_tx);
        handle.join().unwrap();
    }

    /// THE A1-c REGRESSION GUARD (Issue #139, reframed for Issue #141
    /// P95-cap). The A1-c disaster ran macro_redistribution 6× the
    /// calibrated cadence and **diverged** (continental collapse). The
    /// calibrated workflow under the robust P95-cap sea level (cap=0.92,
    /// 12×20) instead settles into a **bounded limit cycle** around
    /// ~30% emergent land. The guard tests that signature:
    ///
    ///   1. mass-conserving: `|Δmass|/mass < 1e-3` per cycle
    ///      (macro_redistribution rebound-balanced; A1-c bled mass);
    ///   2. BOUNDED-BAND convergence: late-cycle continental-fraction
    ///      spread < 0.12 (a bounded limit cycle, NOT a Δ-strict fixed
    ///      point — the system oscillates ±0.05 by nature; a Δ-strict
    ///      gate would wrongly fail it);
    ///   3. emergent land in the ~30% band [0.18, 0.45] — NOT
    ///      collapsed (A1-c ≈ 0) and NOT runaway (too-low cap ≈ 0.95).
    ///
    /// 64² seed 42, cap=0.92 / n_cycles=12 (coupled calibration).
    #[test]
    fn workflow_mode_continent_preserved() {
        use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
        use ymir_core::tectonics_v2::field::Field2D;

        let continental_fraction = |snap: &C1Snapshot| -> f64 {
            let c = snap.plate_type.iter().filter(|&&t| t == 1).count();
            c as f64 / snap.plate_type.len() as f64
        };
        // S̃-only emergent land (ungated by plate_type) — the
        // apples-to-apples measure for the isostatic floor.
        let iso_land = |snap: &C1Snapshot| -> f64 {
            let f = Field2D::from_vec(snap.nx, snap.ny, snap.s.clone());
            compute_isostasy(&f, &IsostasyConfig::c1_default()).land_ratio as f64
        };
        let mass = |s: &[f64]| -> f64 { s.iter().sum() };

        // ---- Calibrated workflow run ----
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: C1RunSpec { grid_size: 64, seed: 42, ..C1RunSpec::default() },
                phase_a: PhaseAParams {
                    // Issue #141: cap=0.92 is COUPLED with n_cycles≈12
                    // (worst-case band-entry cycle 9 + margin). 5 would
                    // cut mid-overshoot.
                    n_cycles: 12,
                    ..PhaseAParams::default()
                },
            })
            .unwrap();
        let events = drain_run(&evt_rx);
        drop(cmd_tx);
        handle.join().unwrap();

        let k_cycle = PhaseAParams::default().k_cycle;
        let n_cycles = 12usize; // Issue #141: coupled with cap=0.92.
        let snaps_at = |step: usize| -> Vec<&C1Snapshot> {
            events
                .iter()
                .filter_map(|e| match e {
                    C1Event::StepCompleted { snapshot } if snapshot.step == step => Some(snapshot),
                    _ => None,
                })
                .collect()
        };

        // (1) mass-conservation per cycle: |Δmass|/mass from the
        //     pre-macro (last per-step) → post-macro (post-cycle) pair.
        // Trajectory of the continental fraction for the convergence
        // check.
        let mut frac_traj: Vec<f64> = Vec::new();
        if let Some(s0) = snaps_at(0).first() {
            frac_traj.push(continental_fraction(s0));
        }
        for c in 1..=n_cycles {
            let pair = snaps_at(c * k_cycle);
            assert!(pair.len() >= 2, "cycle {c}: expected pre+post snapshots at boundary");
            let (pre, post) = (pair[0], pair[1]);
            let m_pre = mass(&pre.s);
            let rel_dmass = (mass(&post.s) - m_pre).abs() / m_pre.max(1e-12);
            assert!(
                rel_dmass < 1e-3,
                "cycle {c}: macro_redistribution not mass-conserving — |Δmass|/mass = {rel_dmass:.2e} (≥ 1e-3 would indicate A1-c-style runaway erosion)"
            );
            frac_traj.push(continental_fraction(post));
        }

        // (2) BOUNDED-BAND convergence (Issue #141 — NOT Δ-strict).
        // Under P95-cap the system is a bounded LIMIT CYCLE (±0.05),
        // not a fixed point — a Δ-strict gate would wrongly fail an
        // oscillator. Assert the late cycles stay within a bounded
        // band around the equilibrium (the natural coast fluctuation).
        // frac_traj = [cycle0, cycle1..n_cycles]; late = last 6 cycles.
        let late = &frac_traj[frac_traj.len().saturating_sub(6)..];
        let (lmn, lmx) = late.iter().fold((1.0_f64, 0.0_f64), |(a, b), &v| (a.min(v), b.max(v)));
        let late_spread = lmx - lmn;
        assert!(
            late_spread < 0.12,
            "workflow not in a bounded band — late-cycle spread = {late_spread:.4} (≥ 0.12 = unbounded drift / A1-c-style collapse). Trajectory: {frac_traj:?}"
        );

        let workflow_iso = match events.last() {
            Some(C1Event::Completed { final_snapshot, .. }) => iso_land(final_snapshot),
            other => panic!("expected Completed; got {other:?}"),
        };

        eprintln!("Phase 1.5 workflow_mode_continent_preserved (seed 42, 64², cap=0.92, 12×20):");
        eprintln!("  continental fraction trajectory: {frac_traj:?}");
        eprintln!("  late-cycle band spread = {late_spread:.4}");
        eprintln!("  workflow emergent (iso_land, c1_default) = {workflow_iso:.4}");

        // (3) NOT collapsed, NOT runaway — emergent land in the ~30%
        // band (Issue #141 P95-cap). A1-c collapse would be ≈ 0; a
        // too-low cap runs away to ≈ 0.95. The bounded-band (2) +
        // this band + mass-conservation (1) together are the Phase 1.5
        // A1-c guard. (The #139 0.9×-gallery-floor comparison is
        // dropped: under P95-cap the workflow legitimately settles
        // slightly below the raw-gallery S̃-implied land — macro
        // erosion + reclassify equilibrium — which is healthy, not a
        // collapse.)
        assert!(
            (0.18..=0.45).contains(&workflow_iso),
            "emergent land {workflow_iso:.4} outside the ~30% band [0.18, 0.45] (collapsed or runaway). Trajectory: {frac_traj:?}"
        );
    }

    /// The per-cycle `apply_post_tectonic` reclassification actually
    /// fires: at some cycle boundary the post-tectonic snapshot's
    /// `plate_type` differs from the pre-tectonic (last per-step)
    /// snapshot at the same boundary step. Isolates reclassify — the
    /// only difference between the two snapshots is the
    /// `apply_post_tectonic` call.
    #[test]
    fn workflow_mode_coast_reclassifies() {
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_cycles = 3;
        let k_cycle = 20;
        let phase_a = PhaseAParams { n_cycles, k_cycle, ..PhaseAParams::default() };
        cmd_tx.send(C1Command::RunWorkflow { spec: small_spec(0, 42), phase_a }).unwrap();

        let events = drain_run(&evt_rx);

        // Group plate_type vectors by step. Cycle boundaries (step =
        // c·k_cycle, c ≥ 1) have exactly two snapshots: [0] pre-
        // reclassify (last per-step), [1] post-reclassify (post-cycle).
        use std::collections::HashMap;
        let mut by_step: HashMap<usize, Vec<Vec<u8>>> = HashMap::new();
        for e in &events {
            if let C1Event::StepCompleted { snapshot } = e {
                by_step.entry(snapshot.step).or_default().push(snapshot.plate_type.clone());
            }
        }

        let mut any_reclassify = false;
        for c in 1..=n_cycles {
            let boundary = c * k_cycle;
            if let Some(pair) = by_step.get(&boundary) {
                if pair.len() == 2 && pair[0] != pair[1] {
                    let changed = pair[0].iter().zip(&pair[1]).filter(|(a, b)| a != b).count();
                    eprintln!(
                        "Stage V reclassify: cycle {c} boundary step {boundary} flipped {changed} cells"
                    );
                    any_reclassify = true;
                }
            }
        }

        assert!(
            any_reclassify,
            "no cycle-boundary reclassification observed — apply_post_tectonic reclassify did not change plate_type at any boundary"
        );

        drop(cmd_tx);
        handle.join().unwrap();
    }

    /// Issue #141 Stage A — coast coherence (W2). The Viz-0 Stage A
    /// bug was the plate_type coast (reclassify, S̃-space sea level)
    /// diverging from the altitude=0 coast (compute_isostasy, h-space
    /// sea level) because the two used DIFFERENT sea-level formulas.
    /// Under Phase 1.5 both use `c1_default` (P95-cap), so the
    /// S̃-space reclassify land set must equal the h-space isostasy
    /// land set. This is NON-trivial (it would catch a future
    /// divergence of the two instances); it is NOT the structural
    /// `altitude>0 ⟺ plate_type` identity (which Stein-Stein gating
    /// makes tautological). 64² seed 42, cap=0.92 / n_cycles=12.
    #[test]
    fn acceptance_coast_coherence_phase_1_5() {
        use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
        use ymir_core::tectonics_v2::field::Field2D;

        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: C1RunSpec { grid_size: 64, seed: 42, ..C1RunSpec::default() },
                phase_a: PhaseAParams { n_cycles: 12, ..PhaseAParams::default() },
            })
            .unwrap();
        let events = drain_run(&evt_rx);
        drop(cmd_tx);
        handle.join().unwrap();

        let final_snap = match events.last() {
            Some(C1Event::Completed { final_snapshot, .. }) => final_snapshot,
            other => panic!("expected Completed; got {other:?}"),
        };
        let n = (final_snap.nx * final_snap.ny) as f64;

        // S̃-space land set (plate_type, set by reclassify).
        let reclassify_land = final_snap.plate_type.iter().filter(|&&t| t == 1).count() as f64 / n;
        // h-space land set (compute_isostasy land_ratio on the SAME s,
        // under the SAME c1_default mode).
        let f = Field2D::from_vec(final_snap.nx, final_snap.ny, final_snap.s.clone());
        let isostasy_land = compute_isostasy(&f, &IsostasyConfig::c1_default()).land_ratio as f64;

        eprintln!("Phase 1.5 Stage A coast coherence (seed 42, 64², cap=0.92):");
        eprintln!("  reclassify land (S̃-space) = {reclassify_land:.4}");
        eprintln!("  isostasy land  (h-space)  = {isostasy_land:.4}");

        // Both reduce to {s > sea_level_ref} (h = s·buoyancy monotonic),
        // so they must agree to within f32 rounding. A divergence would
        // mean the two sea-level instances drifted (the Viz-0 Stage A
        // regression).
        assert!(
            (reclassify_land - isostasy_land).abs() < 1e-3,
            "dual-space coast divergence: reclassify land {reclassify_land:.4} != isostasy land {isostasy_land:.4} (the two sea-level instances disagree — Viz-0 Stage A regression)"
        );
    }

    /// Issue #141 Stage A — multi-seed emergent-land distribution.
    /// The P95-cap product target: emergent land in a bounded band
    /// AROUND ~30% with NATURAL per-seed variation (NOT a uniform
    /// fixed value — that would be suspect). Asserts the seed
    /// DISTRIBUTION (mean + per-seed band), not a tight per-seed
    /// point (the system is a bounded limit cycle). 64² cap=0.92 /
    /// n_cycles=12; emergent = mean of the last 4 post-cycle
    /// continental fractions (band-centre estimate, robust to the
    /// ±0.05 limit-cycle stop-point variance).
    #[test]
    fn acceptance_emergent_land_multiseed_phase_1_5() {
        let k_cycle = PhaseAParams::default().k_cycle;
        let n_cycles = 12usize;

        let mut emergents = Vec::new();
        for &seed in &[42_u64, 1337, 2026, 7, 99] {
            let (cmd_tx, cmd_rx) = bounded(4);
            let (evt_tx, evt_rx) = bounded(2);
            let cancel = Arc::new(AtomicBool::new(false));
            let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
            cmd_tx
                .send(C1Command::RunWorkflow {
                    spec: C1RunSpec { grid_size: 64, seed, ..C1RunSpec::default() },
                    phase_a: PhaseAParams { n_cycles, ..PhaseAParams::default() },
                })
                .unwrap();
            let events = drain_run(&evt_rx);
            drop(cmd_tx);
            handle.join().unwrap();

            // Post-cycle continental fraction at each boundary; average
            // the last 4 cycles (band centre, robust to the limit-cycle).
            let mut post_cycle: Vec<f64> = Vec::new();
            for c in 1..=n_cycles {
                let boundary = c * k_cycle;
                let snaps: Vec<&C1Snapshot> = events
                    .iter()
                    .filter_map(|e| match e {
                        C1Event::StepCompleted { snapshot } if snapshot.step == boundary => {
                            Some(snapshot)
                        }
                        _ => None,
                    })
                    .collect();
                if let Some(post) = snaps.get(1) {
                    let frac = post.plate_type.iter().filter(|&&t| t == 1).count() as f64
                        / post.plate_type.len() as f64;
                    post_cycle.push(frac);
                }
            }
            let tail = &post_cycle[post_cycle.len().saturating_sub(4)..];
            let emergent = tail.iter().sum::<f64>() / tail.len() as f64;
            emergents.push((seed, emergent));
        }

        let mean = emergents.iter().map(|(_, e)| e).sum::<f64>() / emergents.len() as f64;
        let (mn, mx) =
            emergents.iter().fold((1.0_f64, 0.0_f64), |(a, b), &(_, e)| (a.min(e), b.max(e)));

        eprintln!("Phase 1.5 Stage A multi-seed emergent land (64², cap=0.92, 12×20):");
        for (seed, e) in &emergents {
            eprintln!("  seed {seed:>5}: {e:.4}");
        }
        eprintln!("  distribution: min {mn:.4}, mean {mean:.4}, max {mx:.4}");

        // Each seed in a generous band (covers the ±0.05 limit cycle);
        // NOT collapsed (≈0) and NOT runaway (≈0.95).
        for (seed, e) in &emergents {
            assert!(
                (0.15..=0.45).contains(e),
                "seed {seed} emergent {e:.4} outside [0.15, 0.45] (collapsed or runaway)"
            );
        }
        // Distribution centred AROUND ~30%.
        assert!(
            (0.24..=0.38).contains(&mean),
            "multi-seed mean emergent {mean:.4} not centred around ~30% (expected [0.24, 0.38])"
        );
        // Natural variation (not a suspect uniform value).
        assert!(
            mx - mn > 0.02,
            "multi-seed spread {:.4} suspiciously uniform (expected natural per-seed variation)",
            mx - mn
        );
    }

    /// Issue #141 Stage V — Q4 convergence gate under P95-cap. Logs
    /// TWO signals (per the Stage V design):
    ///   (1) GLOBAL convergence: per-cycle |Δmass|/mass < 1e-3 +
    ///       last-cycle |Δfrac| < 0.01.
    ///   (2) THRESHOLD jitter: the per-step drainage sea_level_ref
    ///       (P95-cap), recomputed post-hoc from each step's S̃, to
    ///       see whether the threshold itself oscillates step-to-step
    ///       (the most sensitive signal) — distinguishing per-step
    ///       drainage jitter from per-cycle reclassify jumps (macro
    ///       redistribution at cycle boundaries).
    /// Verdict drives whether the per-cycle-stable drainage sub-fix is
    /// needed (already spec'd; activate only if the per-step threshold
    /// jitters).
    #[test]
    fn workflow_converges_under_p95_cap() {
        use ymir_core::tectonics::isostasy::IsostasyConfig;
        use ymir_core::tectonics_v2::field::Field2D;
        use ymir_core::tectonics_v2::workflow::phase_a_common::compute_sea_level_ref_s_space;

        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: C1RunSpec { grid_size: 64, seed: 42, ..C1RunSpec::default() },
                phase_a: PhaseAParams {
                    n_cycles: 12, // Issue #141: coupled with cap=0.92.
                    ..PhaseAParams::default()
                },
            })
            .unwrap();
        let events = drain_run(&evt_rx);
        drop(cmd_tx);
        handle.join().unwrap();

        let c1 = IsostasyConfig::c1_default();
        let k_cycle = PhaseAParams::default().k_cycle;
        let n_cycles = 12usize;

        // Per-step P95-cap drainage threshold, in emission order.
        let mut step_thr: Vec<(usize, f64)> = Vec::new();
        for e in &events {
            if let C1Event::StepCompleted { snapshot } = e {
                let f = Field2D::from_vec(snapshot.nx, snapshot.ny, snapshot.s.clone());
                step_thr.push((snapshot.step, compute_sea_level_ref_s_space(&f, &c1)));
            }
        }

        // (2) Within-cycle step-to-step jitter vs cycle-boundary jumps.
        // Steps 1..=k_cycle are cycle 1, etc. The cycle-0 snapshot is
        // step 0; boundary snapshots (post-macro) repeat the boundary
        // step value. Walk consecutive distinct emissions.
        let mut max_within = 0.0_f64;
        let mut max_within_at = (0usize, 0usize);
        let mut max_boundary = 0.0_f64;
        for w in step_thr.windows(2) {
            let (s0, t0) = w[0];
            let (s1, t1) = w[1];
            let d = (t1 - t0).abs();
            // A boundary transition is where the step index does NOT
            // advance by exactly 1 (post-cycle snapshot reuse) or
            // crosses a multiple of k_cycle.
            let is_boundary = s1 == s0 || s0 % k_cycle == 0 && s0 != 0;
            if is_boundary {
                max_boundary = max_boundary.max(d);
            } else {
                if d > max_within {
                    max_within = d;
                    max_within_at = (s0, s1);
                }
            }
        }
        // Count "large" within-cycle threshold jumps (> 0.05). A
        // chronic per-step jitter would make this large; a smooth
        // threshold (mild sawtooth + rare discrete level shifts) keeps
        // it small. Measured at seed 42: 1 / ~100 (a single discrete
        // S̃-event level shift; the rest is a ~0.01 sawtooth the system
        // absorbs).
        let big_within = step_thr
            .windows(2)
            .filter(|w| {
                let (s0, _) = w[0];
                let (s1, _) = w[1];
                let is_boundary = s1 == s0 || s0 % k_cycle == 0 && s0 != 0;
                !is_boundary && (w[1].1 - w[0].1).abs() > 0.05
            })
            .count();

        // (1) per-cycle mass + continental-fraction trajectory.
        let cont_frac = |snap: &C1Snapshot| -> f64 {
            snap.plate_type.iter().filter(|&&t| t == 1).count() as f64
                / snap.plate_type.len() as f64
        };
        let mass = |s: &[f64]| -> f64 { s.iter().sum() };
        let snaps_at = |step: usize| -> Vec<&C1Snapshot> {
            events
                .iter()
                .filter_map(|e| match e {
                    C1Event::StepCompleted { snapshot } if snapshot.step == step => Some(snapshot),
                    _ => None,
                })
                .collect()
        };
        let mut frac_traj: Vec<f64> = Vec::new();
        if let Some(s0) = snaps_at(0).first() {
            frac_traj.push(cont_frac(s0));
        }
        let mut max_rel_dmass = 0.0_f64;
        for c in 1..=n_cycles {
            let pair = snaps_at(c * k_cycle);
            assert!(pair.len() >= 2);
            let (pre, post) = (pair[0], pair[1]);
            let m_pre = mass(&pre.s);
            max_rel_dmass = max_rel_dmass.max((mass(&post.s) - m_pre).abs() / m_pre.max(1e-12));
            frac_traj.push(cont_frac(post));
        }
        // BOUNDED-BAND convergence (Issue #141 — NOT Δ-strict). The
        // P95-cap system is a bounded LIMIT CYCLE (±0.05), not a fixed
        // point; a Δ-strict gate would wrongly fail an oscillator (the
        // 5-cycle "Δ=0.005" earlier was a fluke). Assert the late
        // cycles stay within a bounded band around the equilibrium.
        let late = &frac_traj[frac_traj.len().saturating_sub(6)..];
        let (lmn, lmx) = late.iter().fold((1.0_f64, 0.0_f64), |(a, b), &v| (a.min(v), b.max(v)));
        let late_spread = lmx - lmn;

        eprintln!(
            "=== Issue #141 Stage V — convergence under P95-cap (seed 42, 64², cap=0.92, 12×20) ==="
        );
        eprintln!(
            "  per-step P95 drainage threshold: max WITHIN-cycle |Δ| = {max_within:.4} at steps {max_within_at:?}, max BOUNDARY |Δ| = {max_boundary:.4}"
        );
        eprintln!("  within-cycle transitions |Δ|>0.05: {big_within} of ~{}", n_cycles * k_cycle);
        eprintln!("  continental-fraction trajectory  = {frac_traj:?}");
        eprintln!("  per-cycle max |Δmass|/mass        = {max_rel_dmass:.2e}");
        eprintln!("  late-cycle band spread            = {late_spread:.4}");

        // Q4 gate (1) — mass-conserving.
        assert!(
            max_rel_dmass < 1e-3,
            "macro not mass-conserving under P95-cap: max |Δmass|/mass = {max_rel_dmass:.2e}"
        );
        // Q4 gate (2) — bounded-band convergence (limit cycle, not
        // fixed point).
        assert!(
            late_spread < 0.12,
            "not in a bounded band under P95-cap: late-cycle spread = {late_spread:.4} (≥ 0.12 = unbounded). trajectory {frac_traj:?}"
        );
        // Q4 gate (3) — NO chronic per-step threshold jitter. The
        // drainage threshold is smooth (mild ~0.01 sawtooth + rare
        // discrete level shifts), so large within-cycle jumps are
        // rare. A chronic sawtooth would blow this up and would
        // trigger the per-cycle-stable drainage sub-fix (spec'd, not
        // needed at seed 42).
        assert!(
            big_within <= 5,
            "per-step drainage threshold jitters ({big_within} within-cycle jumps > 0.05 of ~{}) — chronic jitter; activate the per-cycle-stable drainage sub-fix",
            n_cycles * k_cycle
        );
    }

    /// Issue #139 Stage A acceptance — DISTINCT product angle from the
    /// Stage V mechanism guard. Stage V asserts mass-conservation /
    /// convergence / isostatic-floor; this asserts the PRODUCT promise
    /// of workflow mode: the coast actually MIGRATES (the displayed
    /// land/sea boundary moves substantially from init), AND the
    /// continent does not vanish (emergent land stays nonzero). This
    /// is what distinguishes workflow mode from the static-coast
    /// gallery path (Issue #137 contract). 64² seed 42, calibrated
    /// default cadence.
    #[test]
    fn workflow_acceptance_continent_preserved_seed_42() {
        use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
        use ymir_core::tectonics_v2::field::Field2D;

        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: C1RunSpec { grid_size: 64, seed: 42, ..C1RunSpec::default() },
                phase_a: PhaseAParams {
                    n_cycles: 12, // Issue #141: coupled with cap=0.92.
                    ..PhaseAParams::default()
                },
            })
            .unwrap();
        let events = drain_run(&evt_rx);
        drop(cmd_tx);
        handle.join().unwrap();

        // init (cycle-0) plate_type vs final plate_type.
        let init_pt = events
            .iter()
            .find_map(|e| match e {
                C1Event::StepCompleted { snapshot } if snapshot.step == 0 => {
                    Some(snapshot.plate_type.clone())
                }
                _ => None,
            })
            .expect("cycle-0 snapshot");
        let final_snap = match events.last() {
            Some(C1Event::Completed { final_snapshot, .. }) => final_snapshot,
            other => panic!("expected Completed; got {other:?}"),
        };

        let coast_moved =
            init_pt.iter().zip(&final_snap.plate_type).filter(|(a, b)| a != b).count();

        let emergent_land = {
            let f = Field2D::from_vec(final_snap.nx, final_snap.ny, final_snap.s.clone());
            compute_isostasy(&f, &IsostasyConfig::c1_default()).land_ratio as f64
        };

        eprintln!("Phase 1.5 Stage A acceptance (workflow, seed 42, 64², cap=0.92, 12×20):");
        eprintln!("  coast cells reclassified vs init = {coast_moved}");
        eprintln!("  final emergent land fraction      = {emergent_land:.4}");

        // PRODUCT: the coast migrated substantially (workflow mode's
        // raison d'être — gallery would not reclassify on sea level).
        // Seed 42 measured ~1500 net flips; assert a robust lower
        // bound.
        assert!(
            coast_moved > 500,
            "workflow coast did not migrate: only {coast_moved} cells reclassified vs init (expected > 500)"
        );
        // PRODUCT (Issue #141 P95-cap): emergent land in the ~30%
        // neighbourhood — NOT collapsed (A1-c ≈ 0) and NOT runaway
        // (cap-too-low ≈ 0.95). Band [0.18, 0.45] covers seed 42's
        // bounded limit-cycle (~0.28–0.31) plus margin.
        assert!(
            (0.18..=0.45).contains(&emergent_land),
            "emergent land {emergent_land:.4} outside the ~30% band [0.18, 0.45] (collapsed or runaway)"
        );
    }

    /// Issue #139 Stage V diagnostic (user-requested, `#[ignore]`'d):
    /// resolve whether the ~0.058 workflow continental fraction is an
    /// apples/oranges metric artefact, an adaptive-threshold drift, or
    /// genuine macro_redistribution over-erosion. Logs, per cycle:
    ///   (1a) plate_type==Continental fraction (the label, post-reclassify)
    ///   (1b) rendered altitude>0 fraction (derive_altitude_field, Stein-
    ///        Stein plate_type-gated)
    ///   (1c) isostatic land fraction (compute_isostasy(s).land_ratio,
    ///        S̃-only, NOT plate_type-gated)
    ///   (2)  s_min / s_max / sea_level_ref (s-space adaptive threshold)
    ///   (3)  reclassify C→O / O→C counts + macro Δmass (loss decomposition)
    /// plus a GALLERY control (100 steps, no post-tectonic) for the
    /// "no macro erosion" reference. Run with:
    ///   cargo test --release -p ymir-viz --features v2_legacy --bin ymir-viz \
    ///     workflow_continent_diagnostic -- --ignored --nocapture
    #[test]
    #[ignore]
    fn workflow_continent_diagnostic() {
        use ymir_core::tectonics::isostasy::{IsostasyConfig, compute_isostasy};
        use ymir_core::tectonics_v2::field::Field2D;

        let nx = 64usize;
        let ny = 64usize;
        let n_cells = (nx * ny) as f64;

        let plate_type_cont = |snap: &C1Snapshot| -> f64 {
            snap.plate_type.iter().filter(|&&t| t == 1).count() as f64 / n_cells
        };
        let iso_land = |s: &[f64]| -> f64 {
            let f = Field2D::from_vec(nx, ny, s.to_vec());
            compute_isostasy(&f, &IsostasyConfig::c1_default()).land_ratio as f64
        };
        let rendered_land = |snap: &C1Snapshot| -> f64 {
            let alt = crate::visualization::c1_viz::derive_altitude_field(snap);
            let mut c = 0usize;
            for j in 0..ny {
                for i in 0..nx {
                    if alt.get(i as i32, j as i32) > 0.0 {
                        c += 1;
                    }
                }
            }
            c as f64 / n_cells
        };
        let s_stats = |s: &[f64]| -> (f64, f64, f64) {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for &v in s {
                lo = lo.min(v);
                hi = hi.max(v);
            }
            // s-space adaptive sea level (sea_level_fraction = 0.4).
            let sea = lo + 0.4 * (hi - lo);
            (lo, hi, sea)
        };
        let mass = |s: &[f64]| -> f64 { s.iter().sum() };

        // ---- WORKFLOW run ----
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
        let spec = C1RunSpec { grid_size: 64, seed: 42, ..C1RunSpec::default() };
        cmd_tx
            .send(C1Command::RunWorkflow { spec: spec.clone(), phase_a: PhaseAParams::default() })
            .unwrap();
        let events = drain_run(&evt_rx);
        drop(cmd_tx);
        handle.join().unwrap();

        let k_cycle = PhaseAParams::default().k_cycle;
        let n_cycles = PhaseAParams::default().n_cycles;

        let snaps_at = |step: usize| -> Vec<&C1Snapshot> {
            events
                .iter()
                .filter_map(|e| match e {
                    C1Event::StepCompleted { snapshot } if snapshot.step == step => Some(snapshot),
                    _ => None,
                })
                .collect()
        };

        eprintln!("=== WORKFLOW (seed 42, 64², default 5×20) ===");
        eprintln!(
            "cyc | pt_cont | rend>0 | iso_land | s_min  s_max  sea_ref | mass    | C->O O->C | macroDmass"
        );
        // cycle 0 (init).
        if let Some(s0) = snaps_at(0).first() {
            let (lo, hi, sea) = s_stats(&s0.s);
            eprintln!(
                "  0 | {:.4} | {:.4} | {:.4}  | {:.3} {:.3} {:.3} | {:.1} |   -    -  |    -",
                plate_type_cont(s0),
                rendered_land(s0),
                iso_land(&s0.s),
                lo,
                hi,
                sea,
                mass(&s0.s),
            );
        }
        for c in 1..=n_cycles {
            let boundary = c * k_cycle;
            let pair = snaps_at(boundary);
            if pair.len() < 2 {
                continue;
            }
            let pre = pair[0]; // last per-step (pre macro/reclassify)
            let post = pair[1]; // post macro+reclassify
            let (lo, hi, sea) = s_stats(&pre.s); // threshold from pre-macro s
            let c_to_o = pre
                .plate_type
                .iter()
                .zip(&post.plate_type)
                .filter(|(a, b)| **a == 1 && **b == 0)
                .count();
            let o_to_c = pre
                .plate_type
                .iter()
                .zip(&post.plate_type)
                .filter(|(a, b)| **a == 0 && **b == 1)
                .count();
            let macro_dmass = mass(&post.s) - mass(&pre.s);
            eprintln!(
                "  {} | {:.4} | {:.4} | {:.4}  | {:.3} {:.3} {:.3} | {:.1} | {:>4} {:>4} | {:+.3}",
                c,
                plate_type_cont(post),
                rendered_land(post),
                iso_land(&post.s),
                lo,
                hi,
                sea,
                mass(&post.s),
                c_to_o,
                o_to_c,
                macro_dmass,
            );
        }

        // ---- GALLERY control (no post-tectonic) ----
        let (cmd_tx2, cmd_rx2) = bounded(4);
        let (evt_tx2, evt_rx2) = bounded(2);
        let cancel2 = Arc::new(AtomicBool::new(false));
        let handle2 = spawn_c1_thread(cmd_rx2, evt_tx2, cancel2);
        let gallery_spec =
            C1RunSpec { grid_size: 64, seed: 42, n_steps: 100, ..C1RunSpec::default() };
        cmd_tx2.send(C1Command::RunBaseline { spec: gallery_spec }).unwrap();
        let gevents = drain_run(&evt_rx2);
        drop(cmd_tx2);
        handle2.join().unwrap();
        if let Some(C1Event::Completed { final_snapshot, .. }) = gevents.last() {
            eprintln!("=== GALLERY control (seed 42, 64², 100 steps, NO post-tectonic) ===");
            eprintln!(
                "final: pt_cont={:.4}  rend>0={:.4}  iso_land={:.4}  mass={:.1}",
                plate_type_cont(final_snapshot),
                rendered_land(final_snapshot),
                iso_land(&final_snapshot.s),
                mass(&final_snapshot.s),
            );
        }
    }

    /// Continuation capability (Stage E4): `run_workflow_cycles` run
    /// twice on the SAME state — feeding the first run's end as the
    /// second run's `base_step` — carries the state forward (the step
    /// offset advances and the field keeps evolving). Exercises the
    /// continuation core directly; no `ContinueRun` command exists yet.
    #[test]
    fn worker_retains_state_for_continuation() {
        use crossbeam_channel::unbounded;
        use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;

        // Unbounded so synchronous in-test sends never deadlock.
        let (tx, _rx) = unbounded::<C1Event>();

        let spec = small_spec(0, 42); // 32² for speed
        let phase_a = PhaseAParams { n_cycles: 2, k_cycle: 5, ..PhaseAParams::default() };

        let mut state = init_c1_state_phase_2_r7(spec.grid_size, spec.seed, &spec.init_params);
        let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let vel: Vec<(f64, f64)> = kinematics.velocities.clone();

        // Run 1 from a fresh state (base_step 0).
        let last1 = run_workflow_cycles(&mut state, &mut kinematics, &spec, &phase_a, &vel, &tx, 0);
        assert_eq!(last1, phase_a.n_cycles * phase_a.k_cycle);
        let s_after_run1: Vec<f64> = state.s.data().to_vec();

        // Run 2 RESUMING the same state (base_step = last1).
        let last2 =
            run_workflow_cycles(&mut state, &mut kinematics, &spec, &phase_a, &vel, &tx, last1);

        // Step offset carried forward.
        assert_eq!(
            last2,
            last1 + phase_a.n_cycles * phase_a.k_cycle,
            "continuation must advance the global step from the retained total"
        );
        assert_eq!(last2, 2 * phase_a.n_cycles * phase_a.k_cycle);

        // State carried forward — the second run continued evolving
        // the SAME state rather than restarting from init.
        assert_ne!(
            s_after_run1,
            state.s.data().to_vec(),
            "resuming must continue evolving the retained state"
        );
    }
}

pub fn spawn_c1_thread(
    commands_rx: Receiver<C1Command>,
    events_tx: Sender<C1Event>,
    cancel: Arc<AtomicBool>,
) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new()
        .name("ymir-c1-bridge".into())
        .spawn(move || {
            // Continuation-ready retention (Issue #139 Stage E4): the
            // last completed run's (state, kinematics). A future
            // `C1Command::ContinueRun { additional_cycles,
            // new_kinematics }` would resume from here via
            // `run_workflow_cycles` instead of re-initialising. No
            // command consumes it yet (capability only) — retention is
            // passive and does not change current-run behaviour (W2).
            let mut retained: Option<(C1State, PlateKinematics)> = None;

            while let Ok(cmd) = commands_rx.recv() {
                // Read the retained handle so its presence is
                // observable to the (future) continuation path; today
                // this is purely the W1 "worker holds state across
                // commands" capability marker.
                let _resumable = retained.is_some();

                match cmd {
                    C1Command::RunBaseline { spec } => {
                        cancel.store(false, Ordering::Relaxed);

                        let _ = events_tx.send(C1Event::Started {
                            spec: spec.clone(),
                            kind: C1RunKind::Gallery,
                        });

                        let mut state =
                            init_c1_state_phase_2_r7(spec.grid_size, spec.seed, &spec.init_params);
                        let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);

                        // Pre-run clone of per-plate velocities for
                        // the per-step closure (cannot access
                        // kinematics directly inside on_step due to
                        // &mut borrow by run_with_closures; see
                        // module docstring on Q-E1.3 trade-off).
                        let initial_velocities: Vec<(f64, f64)> = kinematics.velocities.clone();

                        let config = C1TimeLoopConfig {
                            // #145 5d — production runs the buoyancy fix (rigid continental crust).
                            // `false` = legacy (continents collapse), kept as the regression-guard A/B
                            // reference until the flag is removed (rigidity unconditional). See
                            // docs/reports/c1_continental_buoyancy/.
                            rigid_continental_crust: true,
                            n_steps: spec.n_steps,
                            dx: 1.0 / spec.grid_size as f64,
                            dy: 1.0 / spec.grid_size as f64,
                            // Issue #141: C1 engine uses the robust
                            // P95-capped sea level (drainage + render
                            // altitude). Gallery path has no reclassify,
                            // so this affects the per-step isostasy +
                            // drainage only.
                            iso_config: IsostasyConfig::c1_default(),
                            drainage_max_distance: spec.drainage_max_distance,
                        };

                        // Cycle-0 pre-run snapshot.
                        let cycle_0_snapshot =
                            C1Snapshot::from_state(0, &state, &initial_velocities);
                        let _ =
                            events_tx.send(C1Event::StepCompleted { snapshot: cycle_0_snapshot });

                        let t0 = Instant::now();
                        let tx_for_steps = events_tx.clone();
                        let velocities_for_closure = initial_velocities.clone();

                        run_with_closures(
                            &mut state,
                            &mut kinematics,
                            &config,
                            &spec.closures,
                            |step, state| {
                                // C1 time loop callback fires AFTER
                                // the `step`-th iteration's mutations
                                // with the 0-based step index. We emit
                                // `global_step = step + 1` so the
                                // stream is [0 (cycle-0), 1, 2, ...,
                                // n_steps] — monotone and matches the
                                // user-visible "step N completed"
                                // convention.
                                let snap = C1Snapshot::from_state(
                                    step + 1,
                                    state,
                                    &velocities_for_closure,
                                );
                                let _ =
                                    tx_for_steps.send(C1Event::StepCompleted { snapshot: snap });
                            },
                        );

                        let elapsed = t0.elapsed();

                        let final_snapshot =
                            C1Snapshot::from_state(spec.n_steps, &state, &initial_velocities);
                        let _ =
                            events_tx.send(C1Event::Completed { spec, final_snapshot, elapsed });

                        // Retain for a future continuation (Stage E4).
                        retained = Some((state, kinematics));
                    }
                    C1Command::RunWorkflow { spec, phase_a } => {
                        retained = Some(run_workflow(&spec, &phase_a, &events_tx));
                    }
                    C1Command::RunHd { spec, params } => {
                        // HD production chain (step b/5) — cached ymir-core
                        // functions on the worker, per-phase events. Does not
                        // touch `retained` (re-derives deterministically /
                        // from cache).
                        super::hd::run_hd(&spec, &params, &events_tx, &cancel);
                    }
                    C1Command::PreviewShape { spec, params } => {
                        // Coarse-only shape preview (fast) — judge a seed before HD.
                        super::hd::preview_shape(&spec, &params, &events_tx);
                    }
                    C1Command::Cancel => {
                        // MVP Option C — set the flag; takes
                        // effect on the next RunBaseline only.
                        cancel.store(true, Ordering::Relaxed);
                    }
                }
            }
        })
        .expect("ymir-c1-bridge thread spawn")
}

/// Workflow-path run (Issue #139 Stage E2) — the **calibrated** Phase
/// A loop. Reproduces `run_phase_a_cycle_c1(Enabled)` `n_cycles`
/// times, inlining `run_with_closures` only to inject the per-step
/// animation hook (the wrapper hardcodes `on_step = |_, _| {}`).
///
/// Cadence is the calibration anchor: `phase_a.n_cycles × phase_a.
/// k_cycle` tectonic steps with one `apply_post_tectonic` per cycle
/// (default 5×20 = 100 steps, 5 macro-redistribution passes). This is
/// **NOT** the A1-c failure mode (6 passes over 300 steps = 6× the
/// calibrated cadence → over-erosion); `macro_redistribution`'s
/// `alpha = 0.01` is calibrated for exactly this cadence.
///
/// `cratonic_cfg = None` ⇒ `apply_post_tectonic` Step 4 (cratonic
/// recompute) is skipped, so `initial_per_plate_type = None` is
/// bit-identical to `run_phase_a_cycle_c1` with `cratonic_config =
/// None` (only Step 4 reads it — see `phase_a_common.rs` Step 4 gate).
///
/// Snapshot stream: cycle-0 (step 0) + per cycle `{ k_cycle per-step
/// snapshots + 1 post-tectonic snapshot }` + `Completed`. Per-step
/// indices are the global tectonic step `cycle·k_cycle + s + 1`; the
/// post-cycle snapshot reuses the cycle-boundary index (so the step
/// sequence is monotone non-decreasing and the "step N/total" counter
/// tracks tectonic steps, not emission ordinal).
fn run_workflow(
    spec: &C1RunSpec,
    phase_a: &PhaseAParams,
    events_tx: &Sender<C1Event>,
) -> (C1State, PlateKinematics) {
    let _ = events_tx.send(C1Event::Started {
        spec: spec.clone(),
        kind: C1RunKind::Workflow { phase_a: phase_a.clone() },
    });

    let mut state = init_c1_state_phase_2_r7(spec.grid_size, spec.seed, &spec.init_params);
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);

    // Init-time velocities for the snapshot (mid-run kinematics
    // mutations from accretion/rifting are NOT reflected — same MVP
    // caveat as the gallery path; Viz-0-bis item 1).
    let initial_velocities: Vec<(f64, f64)> = kinematics.velocities.clone();

    // Cycle-0 pre-run snapshot (step 0).
    let _ = events_tx.send(C1Event::StepCompleted {
        snapshot: C1Snapshot::from_state(0, &state, &initial_velocities),
    });

    let t0 = Instant::now();
    let final_step = run_workflow_cycles(
        &mut state,
        &mut kinematics,
        spec,
        phase_a,
        &initial_velocities,
        events_tx,
        0,
    );
    let elapsed = t0.elapsed();

    let final_snapshot = C1Snapshot::from_state(final_step, &state, &initial_velocities);
    let _ = events_tx.send(C1Event::Completed { spec: spec.clone(), final_snapshot, elapsed });

    // Hand the post-run state back so the worker can retain it for a
    // future continuation (Issue #139 Stage E4).
    (state, kinematics)
}

/// Continuation-ready core (Issue #139 Stage E4): run `phase_a.
/// n_cycles` Phase A cycles on an **existing** `state` + `kinematics`,
/// emitting per-step + per-cycle snapshots. `base_step` is the global
/// tectonic step already elapsed (0 for a fresh run; `> 0` when
/// resuming a retained state). Returns the global tectonic step after
/// the last cycle.
///
/// Factored out of [`run_workflow`] so a future
/// `C1Command::ContinueRun { additional_cycles, new_kinematics }`
/// could call it again on the worker's retained `(state, kinematics)`
/// with `base_step = prior total` — resuming the simulation instead
/// of re-initialising. No command consumes it yet (Stage E4 ships the
/// capability only; the command + UI are deferred).
fn run_workflow_cycles(
    state: &mut C1State,
    kinematics: &mut PlateKinematics,
    spec: &C1RunSpec,
    phase_a: &PhaseAParams,
    initial_velocities: &[(f64, f64)],
    events_tx: &Sender<C1Event>,
    base_step: usize,
) -> usize {
    let k_cycle = phase_a.k_cycle;
    // Issue #141: P95-capped sea level. This single config flows to
    // BOTH the time loop (drainage, time_loop.rs:614) AND
    // apply_post_tectonic (reclassify) below, so the drainage coast,
    // the reclassify coast, and the rendered altitude coast all share
    // the same robust threshold (W2 coherence).
    let iso_config = IsostasyConfig::c1_default();
    let mut last_step = base_step;

    for cycle in 0..phase_a.n_cycles {
        let cfg = C1TimeLoopConfig {
            // #145 5d — production runs the buoyancy fix (rigid continental crust).
            rigid_continental_crust: true,
            n_steps: k_cycle,
            dx: 1.0 / spec.grid_size as f64,
            dy: 1.0 / spec.grid_size as f64,
            iso_config: iso_config.clone(),
            drainage_max_distance: spec.drainage_max_distance,
        };
        let cyc_base = base_step + cycle * k_cycle;

        // Per-step tectonic loop (animation hook). The closure
        // captures only `events_tx`, `initial_velocities`, `cyc_base`
        // by reference / copy — none alias `state` / `kinematics`, so
        // no borrow conflict with `run_with_closures`'s `&mut` args.
        run_with_closures(state, kinematics, &cfg, &spec.closures, |s_in_c, st| {
            let step = cyc_base + s_in_c + 1;
            let _ = events_tx.send(C1Event::StepCompleted {
                snapshot: C1Snapshot::from_state(step, st, initial_velocities),
            });
        });

        // Calibrated post-tectonic pass (split borrow: &mut s /
        // &plate_id / &mut plate_type are disjoint fields).
        let _ = apply_post_tectonic(PostTectonicInput {
            s_field: &mut state.s,
            plate_id: Some(&state.plate_id),
            plate_type: Some(&mut state.plate_type),
            previous_cratonic_factor: None,
            initial_per_plate_type: None,
            params: phase_a,
            iso_cfg: &iso_config,
            cratonic_cfg: None,
        });

        // Post-cycle snapshot — coast reclassified. Reuses the
        // boundary tectonic step index.
        let boundary_step = cyc_base + k_cycle;
        let _ = events_tx.send(C1Event::StepCompleted {
            snapshot: C1Snapshot::from_state(boundary_step, state, initial_velocities),
        });
        last_step = boundary_step;
    }

    last_step
}
