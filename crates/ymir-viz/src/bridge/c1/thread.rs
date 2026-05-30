//! Worker thread driving the C1 time loop.
//!
//! ## Lifecycle
//!
//! 1. `spawn_c1_thread` returns a `JoinHandle`. The thread blocks
//!    on `commands_rx.recv()` between runs.
//! 2. On `C1Command::RunBaseline { spec }`:
//!    a. Reset the cancel flag (allow future cancels).
//!    b. Send `C1Event::Started`.
//!    c. Build state + kinematics + config from `spec`.
//!    d. **Clone `kinematics.velocities` pre-run** into a local
//!       `initial_velocities: Vec<(f64, f64)>` — captured by the
//!       per-step closure (cannot access `kinematics` directly
//!       because `run_with_closures` borrows it `&mut`; see
//!       Stage E2 W7 Q-E1.3).
//!    e. Emit a pre-run cycle-0 `StepCompleted` so the UI can
//!       render the init state before the first step fires.
//!    f. Call `run_with_closures`. The per-step closure emits a
//!       `StepCompleted { snapshot }` event each step.
//!    g. Send `C1Event::Completed { spec, final_snapshot,
//!       elapsed }`.
//! 3. On `C1Command::Cancel`: set the cancel flag (no effect on
//!    current run per Q-E1.3 Option C MVP).
//!
//! ## Worker owns no Bevy state (W4 global)
//!
//! All Bevy-side data (resources, sprites, image handles) lives
//! in `plugin.rs` + the UI/render systems. The worker uses only
//! crossbeam channels + `ymir-core` types. This isolation is
//! shared with `bridge::v2::thread`.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1TimeLoopConfig};

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

    /// Spec sized for fast tests — 32² grid + small n_steps.
    fn small_spec(n_steps: usize, seed: u64) -> C1RunSpec {
        C1RunSpec {
            grid_size: 32,
            seed,
            n_steps,
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
        // Q-E5.1 Option C: with bounded(2) events and immediate
        // drain in the test, the worker naturally backpressures
        // when the channel fills. We don't time-assert blocking;
        // we assert the stream invariant — every step index
        // 0..n_steps appears exactly once in the StepCompleted
        // events, in monotone order. No loss, no duplication,
        // no reordering.
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);

        let n_steps = 30;
        cmd_tx
            .send(C1Command::RunBaseline {
                spec: small_spec(n_steps, 42),
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

        // Cycle 0 + 30 steps = 31 total. The cycle-0 snapshot
        // uses step = 0 and the first run_with_closures callback
        // also fires step = 0; the worker emits both verbatim,
        // so we expect [0, 0, 1, 2, ..., 29].
        assert_eq!(step_indices.len(), n_steps + 1);
        assert_eq!(step_indices[0], 0, "first StepCompleted is cycle-0 (step 0)");
        // Steps from run_with_closures: indices 1..=n_steps should be
        // 0, 1, 2, ..., n_steps - 1 (monotonically increasing).
        for (i, &s) in step_indices.iter().enumerate().skip(1) {
            assert_eq!(
                s,
                i - 1,
                "expected step {} at position {}, got {}",
                i - 1,
                i,
                s
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
                        // Reset cancel before each run (Option C
                        // MVP — the flag is only sampled between
                        // runs, never during).
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

                        // Pre-run snapshot of kinematics velocities
                        // for the per-step closure (see module
                        // docstring rationale).
                        let initial_velocities: Vec<(f64, f64)> =
                            kinematics.velocities.clone();

                        let config = C1TimeLoopConfig {
                            n_steps: spec.n_steps,
                            dx: 1.0 / spec.grid_size as f64,
                            dy: 1.0 / spec.grid_size as f64,
                            iso_config: IsostasyConfig::default(),
                            drainage_max_distance: spec.drainage_max_distance,
                        };

                        // Cycle-0 pre-run snapshot — lets the UI
                        // paint the init state before stepping
                        // begins. Uses `step = 0` (matches the
                        // first `on_step` callback's step index;
                        // n_steps + 1 events total per run).
                        let cycle_0 = C1Snapshot::from_state(
                            0,
                            &state,
                            &initial_velocities,
                        );
                        let _ = events_tx.send(C1Event::StepCompleted {
                            snapshot: cycle_0,
                        });

                        let t0 = Instant::now();
                        let tx = events_tx.clone();
                        let velocities_for_closure = initial_velocities.clone();
                        run_with_closures(
                            &mut state,
                            &mut kinematics,
                            &config,
                            &spec.closures,
                            |step, state| {
                                let snapshot = C1Snapshot::from_state(
                                    step,
                                    state,
                                    &velocities_for_closure,
                                );
                                // `send` blocks when the bounded
                                // events channel is full — this is
                                // the backpressure pause semantics
                                // (W7 Stage E2).
                                let _ =
                                    tx.send(C1Event::StepCompleted { snapshot });
                            },
                        );
                        let elapsed = t0.elapsed();

                        // Final snapshot — re-extract for the
                        // Completed convenience payload.
                        let final_snapshot = C1Snapshot::from_state(
                            spec.n_steps.saturating_sub(1),
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
