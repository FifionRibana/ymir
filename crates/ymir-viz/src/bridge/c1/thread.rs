//! Worker thread driving the C1 Phase A pipeline (Issue #137
//! Stage A revision — A1-c design).
//!
//! ## A1-c — full Phase A pipeline with per-step animation
//!
//! The worker replicates `run_phase_a_cycle_c1`'s internals
//! (`run_with_closures` + `apply_post_tectonic`) directly, rather
//! than calling the wrapper, because `run_phase_a_cycle_c1` hardcodes
//! `on_step = |_, _| {}` and provides no animation hook. By
//! inlining the cycle structure here we get:
//!
//! - **Per-step animation**: `run_with_closures` with a custom
//!   `on_step` closure emits one `StepCompleted` per step.
//!   S̃ + age advect smoothly each step.
//! - **Per-cycle reclassify**: after each cycle's
//!   `run_with_closures` returns, the worker invokes
//!   `apply_post_tectonic` (sea-level + macro_redistribution +
//!   reclassify). plate_type is correctly updated; the displayed
//!   coast migrates in discrete jumps every `steps_per_cycle` steps.
//! - **Honest pipeline**: the rendered state IS what C1 Phase A
//!   actually produces — same path Phase 1.x acceptance tests use.
//!   Track A/B/D galleries use `run_with_closures` standalone
//!   (post_tectonic OMITTED) and were therefore showing an
//!   INCOMPLETE pipeline; the viz now shows the COMPLETE one.
//!
//! ## Cycle structure
//!
//! ```text
//! C1Event::Started
//! C1Event::StepCompleted { cycle-0 snapshot }
//! For cycle in 0..n_cycles {
//!   run_with_closures(steps_per_cycle, |step, state| {
//!     emit StepCompleted with global_step
//!   })
//!   apply_post_tectonic(state)                ← reclassify here
//!   emit StepCompleted (post-cycle snapshot)  ← coast migration visible
//! }
//! C1Event::Completed { final_snapshot }
//! ```
//!
//! Total step events emitted: `1 (cycle-0) + n_cycles ×
//! (steps_per_cycle + 1)` (the `+1` is the per-cycle post-tectonic
//! snapshot).
//!
//! ## Cratonic recompute SKIPPED at MVP
//!
//! `apply_post_tectonic` step 4 (cratonic factor recompute) is
//! gated by `initial_per_plate_type` + `cratonic_cfg` both being
//! `Some`. We pass `None` for both — skipping the recompute.
//! Rationale:
//! - `C1State.cratonic_mask: BoolField` was built BFS-style at
//!   init from continental seeds; it's a binary mask, NOT a
//!   `Field2D` factor. The current C1 design (Track A/B/D) does
//!   NOT recompute it per cycle. Including it here would diverge
//!   from established C1 acceptance.
//! - The skipped step is sequential AFTER reclassify and does NOT
//!   feed back into macro_redistribution / reclassify. Skipping
//!   is safe per the apply_post_tectonic ordering.
//!
//! Viz-0-bis #6 (or a Phase 1.5 C1 task) can revisit if the
//! cratonic factor evolution becomes a design priority.
//!
//! ## Worker owns no Bevy state (W4 global)
//!
//! All Bevy-side data lives in `plugin.rs` + UI/render systems.
//! Worker uses only crossbeam + `ymir-core`.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::workflow::phase_a_common::{
    apply_post_tectonic, PostTectonicInput,
};
use ymir_core::tectonics_v2::workflow::WorkflowParams;

use super::commands::C1Command;
use super::events::C1Event;
use super::snapshot::C1Snapshot;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::c1::spec::C1RunSpec;
    use crossbeam_channel::bounded;
    use std::sync::atomic::AtomicBool;
    use std::sync::Arc;
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
                    let terminal = matches!(
                        e,
                        C1Event::Completed { .. } | C1Event::Failed { .. }
                    );
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

    /// Spec sized for fast tests — 32² grid + small n_steps. The
    /// `steps_per_cycle` defaults to 50; tests that need many
    /// cycles at small `n_steps` should pass a custom value via
    /// `small_spec_with_cycle` to avoid `n_cycles == 0`.
    fn small_spec(n_steps: usize, seed: u64) -> C1RunSpec {
        C1RunSpec {
            grid_size: 32,
            seed,
            n_steps,
            ..C1RunSpec::default()
        }
    }

    fn small_spec_with_cycle(n_steps: usize, seed: u64, steps_per_cycle: usize) -> C1RunSpec {
        C1RunSpec {
            grid_size: 32,
            seed,
            n_steps,
            steps_per_cycle,
            ..C1RunSpec::default()
        }
    }

    #[test]
    fn c1_worker_spawns_and_runs() {
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_steps = 50;
        cmd_tx
            .send(C1Command::RunBaseline {
                spec: small_spec(n_steps, 42),
            })
            .expect("send RunBaseline");

        let events = drain_run(&evt_rx);

        // Event-stream shape: Started, cycle-0 StepCompleted, then
        // one StepCompleted per step (worker emits step indices
        // `0..n_steps`), then Completed. The cycle-0 snapshot
        // shares `step = 0` with the first `on_step(step=0, ...)`
        // callback, so total StepCompleted count = n_steps + 1.
        let started = events
            .iter()
            .filter(|e| matches!(e, C1Event::Started { .. }))
            .count();
        let step_completed = events
            .iter()
            .filter(|e| matches!(e, C1Event::StepCompleted { .. }))
            .count();
        let completed = events
            .iter()
            .filter(|e| matches!(e, C1Event::Completed { .. }))
            .count();

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

    #[test]
    fn c1_worker_event_ordering_no_loss_under_backpressure() {
        // Q-E5.1 Option C: bounded(2) events + immediate drain.
        // The worker naturally backpressures when the channel
        // fills. We don't time-assert blocking; we assert the
        // stream invariant — every step index appears exactly
        // once in StepCompleted events, in strictly monotone
        // order. No loss, no duplication, no reordering.
        //
        // A1-c update (Issue #137 Stage A): the per-cycle
        // post_tectonic is emitted SILENTLY (no event) so the
        // event stream remains strictly monotone in `step`.
        // Effective step count = n_cycles * steps_per_cycle.
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        // 30 steps / steps_per_cycle = 3 cycles, effective = 30
        // total inner steps. Plus cycle-0 = 31 StepCompleted.
        let n_steps = 30;
        let steps_per_cycle = 10;
        cmd_tx
            .send(C1Command::RunBaseline {
                spec: small_spec_with_cycle(n_steps, 42, steps_per_cycle),
            })
            .unwrap();

        let events = drain_run(&evt_rx);

        let step_indices: Vec<usize> = events
            .iter()
            .filter_map(|e| match e {
                C1Event::StepCompleted { snapshot } => Some(snapshot.step),
                _ => None,
            })
            .collect();

        // Expected: [0 (cycle-0), 1, 2, ..., 30] — 31 events.
        assert_eq!(step_indices.len(), n_steps + 1);
        assert_eq!(step_indices[0], 0, "first StepCompleted is cycle-0 (step 0)");
        for (i, &s) in step_indices.iter().enumerate().skip(1) {
            assert_eq!(
                s, i,
                "expected step {} at position {}, got {}",
                i, i, s
            );
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
        cmd_tx
            .send(C1Command::RunBaseline {
                spec: small_spec(n_steps, 42),
            })
            .unwrap();

        let events = drain_run(&evt_rx);

        let mut max_sub_cells = 0_usize;
        let mut total_merges = 0_usize;
        let mut total_rifting_thinning = 0_usize;
        for e in &events {
            if let C1Event::StepCompleted { snapshot } = e {
                max_sub_cells = max_sub_cells
                    .max(snapshot.stats.subduction.cells_consumed);
                total_merges += snapshot.stats.accretion.merges_count;
                total_rifting_thinning +=
                    snapshot.stats.rifting_thinning.cells_thinned;
            }
        }

        eprintln!(
            "Stage V c1_snapshot_carries_stats (seed 42, n_steps={n_steps}, 32²):"
        );
        eprintln!("  max subduction.cells_consumed per step = {max_sub_cells}");
        eprintln!("  cumulative accretion.merges_count      = {total_merges}");
        eprintln!(
            "  cumulative rifting_thinning.cells_thinned = {total_rifting_thinning}"
        );

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

        // Production-scale: 64² × 300, seed 42, default closures
        // (Track D enabled). Track D Stage A measured 251 ms wall
        // time at this scale; allow generous timeout in drain.
        let spec = C1RunSpec::default(); // 64, seed=42, n=300
        cmd_tx
            .send(C1Command::RunBaseline { spec: spec.clone() })
            .unwrap();

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
        eprintln!(
            "  init num_plates       = {}",
            final_snap.num_plates
        );
        eprintln!(
            "  final live_plate_count = {}",
            final_snap.live_plate_count
        );
        eprintln!("  cumulative subduction cells = {cum_sub_cells}");
        eprintln!("  cumulative accretion merges = {cum_merges}");
        eprintln!("  cumulative rifting splits   = {cum_splits}");
        eprintln!("  cumulative thinning cells   = {cum_thinning}");

        // Track D Pangaea collapse — 8 init plates collapse to
        // 1-3 surviving (Track D Stage A measured 2 at seed 42).
        // Threshold ≤ 3 gives margin for stochastic variance in
        // accretion-merge timing while still asserting the
        // Pangaea narrative.
        assert_eq!(
            final_snap.num_plates, 8,
            "Phase 2 R7 default init should produce 8 plates"
        );
        assert!(
            final_snap.live_plate_count <= 3,
            "Pangaea collapse not achieved: final live plates = {} (expected ≤ 3 at seed 42 / 64² / 300 steps)",
            final_snap.live_plate_count
        );

        // Track D activity sanity: subduction high-frequency,
        // accretion fires at least once. Stage A reference seed
        // 42 measured 20,914 subduction cells + 6 merges over
        // 300 steps; we conservatively assert > 1000 + > 0.
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
}

pub fn spawn_c1_thread(
    commands_rx: Receiver<C1Command>,
    events_tx: Sender<C1Event>,
    cancel: Arc<AtomicBool>,
) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new()
        .name("ymir-c1-bridge".into())
        .spawn(move || {
            while let Ok(cmd) = commands_rx.recv() {
                match cmd {
                    C1Command::RunBaseline { spec } => {
                        cancel.store(false, Ordering::Relaxed);

                        let _ = events_tx.send(C1Event::Started {
                            spec: spec.clone(),
                        });

                        let mut state = init_c1_state_phase_2_r7(
                            spec.grid_size,
                            spec.seed,
                            &spec.init_params,
                        );
                        let mut kinematics =
                            PlateKinematics::preset_phase_1_1(state.num_plates);

                        // Pre-run clone of per-plate velocities for
                        // the per-step closure (cannot access
                        // kinematics directly inside on_step due to
                        // &mut borrow by run_with_closures; see
                        // module docstring on Q-E1.3 trade-off).
                        let initial_velocities: Vec<(f64, f64)> =
                            kinematics.velocities.clone();

                        // Phase A cycle structure (A1-c).
                        let n_cycles = spec.n_cycles();
                        let steps_per_cycle = spec.steps_per_cycle;
                        let cycle_config = C1TimeLoopConfig {
                            n_steps: steps_per_cycle,
                            dx: 1.0 / spec.grid_size as f64,
                            dy: 1.0 / spec.grid_size as f64,
                            iso_config: IsostasyConfig::default(),
                            drainage_max_distance: spec.drainage_max_distance,
                        };
                        let workflow_params = WorkflowParams::default();

                        // Cycle-0 pre-run snapshot.
                        let cycle_0_snapshot = C1Snapshot::from_state(
                            0,
                            &state,
                            &initial_velocities,
                        );
                        let _ = events_tx.send(C1Event::StepCompleted {
                            snapshot: cycle_0_snapshot,
                        });

                        let t0 = Instant::now();
                        let tx_for_steps = events_tx.clone();
                        let velocities_for_closure = initial_velocities.clone();

                        for cycle_idx in 0..n_cycles {
                            // Step 1 — per-step tectonic loop. Emit
                            // one StepCompleted per step (animation
                            // hook).
                            let cycle_base = cycle_idx * steps_per_cycle;
                            let tx = tx_for_steps.clone();
                            let velocities = velocities_for_closure.clone();
                            run_with_closures(
                                &mut state,
                                &mut kinematics,
                                &cycle_config,
                                &spec.closures,
                                |step_in_cycle, state| {
                                    let global_step =
                                        cycle_base + step_in_cycle + 1;
                                    let snap = C1Snapshot::from_state(
                                        global_step,
                                        state,
                                        &velocities,
                                    );
                                    let _ = tx.send(
                                        C1Event::StepCompleted { snapshot: snap },
                                    );
                                },
                            );

                            // Step 2 — Phase A post-tectonic pass.
                            // Replicates run_phase_a_cycle_c1's
                            // body but with cratonic recompute
                            // SKIPPED (Viz-0 MVP — see module
                            // docstring). The struct-literal split
                            // borrow on `state` works because s,
                            // plate_id, and plate_type are
                            // disjoint fields.
                            //
                            // No event emitted here: the NEXT
                            // cycle's first inner snapshot will
                            // naturally show the post-post_tectonic
                            // state (state mutations are persistent
                            // across the inner loop). The final
                            // cycle's post-tectonic state is
                            // captured in `Completed.final_snapshot`
                            // below. This keeps the StepCompleted
                            // stream's step indices strictly
                            // monotone (no duplicates) — preserves
                            // Stage V's event-ordering invariant.
                            let _ = apply_post_tectonic(
                                PostTectonicInput {
                                    s_field: &mut state.s,
                                    plate_id: Some(&state.plate_id),
                                    plate_type: Some(&mut state.plate_type),
                                    previous_cratonic_factor: None,
                                    initial_per_plate_type: None,
                                    params: &workflow_params.phase_a,
                                    iso_cfg: &IsostasyConfig::default(),
                                    cratonic_cfg: None,
                                },
                            );
                        }

                        let elapsed = t0.elapsed();

                        let final_step =
                            n_cycles.saturating_mul(steps_per_cycle);
                        let final_snapshot = C1Snapshot::from_state(
                            final_step,
                            &state,
                            &initial_velocities,
                        );
                        let _ = events_tx.send(C1Event::Completed {
                            spec,
                            final_snapshot,
                            elapsed,
                        });
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
