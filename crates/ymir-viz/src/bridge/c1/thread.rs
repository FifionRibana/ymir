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

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::workflow::phase_a_common::{apply_post_tectonic, PostTectonicInput};
use ymir_core::tectonics_v2::workflow::PhaseAParams;

use super::commands::C1Command;
use super::events::{C1Event, C1RunKind};
use super::snapshot::C1Snapshot;
use super::spec::C1RunSpec;

#[cfg(test)]
mod tests {
    use super::*;
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
        // 1..=n_steps), then Completed. Cycle-0 uses step=0; the
        // first `on_step` callback fires AFTER step 1 mutations
        // and gets `global_step = 1`. Total = n_steps + 1.
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
        // After A1-c revert: standalone `run_with_closures` —
        // events are `[step=0 (cycle-0), step=1, step=2, ...,
        // step=n_steps]`. Total = n_steps + 1.
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

        assert_eq!(
            final_snap.num_plates, 8,
            "Phase 2 R7 default init should produce 8 plates"
        );
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

                        let config = C1TimeLoopConfig {
                            n_steps: spec.n_steps,
                            dx: 1.0 / spec.grid_size as f64,
                            dy: 1.0 / spec.grid_size as f64,
                            iso_config: IsostasyConfig::default(),
                            drainage_max_distance: spec.drainage_max_distance,
                        };

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
                                let _ = tx_for_steps.send(
                                    C1Event::StepCompleted { snapshot: snap },
                                );
                            },
                        );

                        let elapsed = t0.elapsed();

                        let final_snapshot = C1Snapshot::from_state(
                            spec.n_steps,
                            &state,
                            &initial_velocities,
                        );
                        let _ = events_tx.send(C1Event::Completed {
                            spec,
                            final_snapshot,
                            elapsed,
                        });

                        // Retain for a future continuation (Stage E4).
                        retained = Some((state, kinematics));
                    }
                    C1Command::RunWorkflow { spec, phase_a } => {
                        retained =
                            Some(run_workflow(&spec, &phase_a, &events_tx));
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
        kind: C1RunKind::Workflow {
            phase_a: phase_a.clone(),
        },
    });

    let mut state =
        init_c1_state_phase_2_r7(spec.grid_size, spec.seed, &spec.init_params);
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

    let final_snapshot =
        C1Snapshot::from_state(final_step, &state, &initial_velocities);
    let _ = events_tx.send(C1Event::Completed {
        spec: spec.clone(),
        final_snapshot,
        elapsed,
    });

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
    let iso_config = IsostasyConfig::default();
    let mut last_step = base_step;

    for cycle in 0..phase_a.n_cycles {
        let cfg = C1TimeLoopConfig {
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
        run_with_closures(
            state,
            kinematics,
            &cfg,
            &spec.closures,
            |s_in_c, st| {
                let step = cyc_base + s_in_c + 1;
                let _ = events_tx.send(C1Event::StepCompleted {
                    snapshot: C1Snapshot::from_state(step, st, initial_velocities),
                });
            },
        );

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
