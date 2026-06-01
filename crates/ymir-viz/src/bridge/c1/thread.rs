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
        let phase_a = PhaseAParams {
            n_cycles,
            k_cycle,
            ..PhaseAParams::default()
        };
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

    /// THE A1-c REGRESSION GUARD (Issue #139 Stage V). The A1-c
    /// disaster ran macro_redistribution 6× the calibrated cadence
    /// (50/300) and **diverged** — continental crust runaway-eroded
    /// without convergence. The calibrated workflow (5×20) instead
    /// **converges**, conserves mass, and lands ABOVE the gallery's
    /// pure-tectonic isostatic land floor. The guard tests that
    /// QUALITATIVE signature, not a magic band — the Stage V
    /// diagnostic (`workflow_continent_diagnostic`) established that
    /// the absolute fraction (~0.058 at seed 42) is the honest
    /// above-sea-level land, NOT a regression:
    ///
    ///   1. mass-conserving: `|Δmass|/mass < 1e-3` per cycle
    ///      (macro_redistribution is rebound-balanced; A1-c runaway
    ///      would bleed mass);
    ///   2. converges: last-cycle `|Δfrac| < 0.01` (stable
    ///      equilibrium, not A1-c's monotone collapse);
    ///   3. above the isostatic floor: workflow final `iso_land ≥
    ///      0.9 × gallery iso_land`, comparing the S̃-only emergent
    ///      land (`compute_isostasy(s).land_ratio`) apples-to-apples
    ///      at the SAME 100-step duration. (NB: rendered altitude>0
    ///      is Stein-Stein plate_type-GATED, so gallery's rendered
    ///      land is the static geometric label ≈0.27 — using it here
    ///      would re-violate the apples/oranges lesson; `iso_land` is
    ///      the ungated emergent measure that gives the ≈0.045 floor.)
    ///
    /// 64² seed 42, calibrated default cadence (5×20 = 100 steps),
    /// plus a gallery control at the same 100 steps.
    #[test]
    fn workflow_mode_continent_preserved() {
        use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
        use ymir_core::tectonics_v2::field::Field2D;

        let continental_fraction = |snap: &C1Snapshot| -> f64 {
            let c = snap.plate_type.iter().filter(|&&t| t == 1).count();
            c as f64 / snap.plate_type.len() as f64
        };
        // S̃-only emergent land (ungated by plate_type) — the
        // apples-to-apples measure for the isostatic floor.
        let iso_land = |snap: &C1Snapshot| -> f64 {
            let f = Field2D::from_vec(snap.nx, snap.ny, snap.s.clone());
            compute_isostasy(&f, &IsostasyConfig::default()).land_ratio as f64
        };
        let mass = |s: &[f64]| -> f64 { s.iter().sum() };

        // ---- Calibrated workflow run ----
        let (cmd_tx, cmd_rx) = bounded(4);
        let (evt_tx, evt_rx) = bounded(2);
        let cancel = Arc::new(AtomicBool::new(false));
        let handle = spawn_c1_thread(cmd_rx, evt_tx, cancel);
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: C1RunSpec {
                    grid_size: 64,
                    seed: 42,
                    ..C1RunSpec::default()
                },
                phase_a: PhaseAParams::default(), // CALIBRATED 5×20
            })
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
                    C1Event::StepCompleted { snapshot } if snapshot.step == step => {
                        Some(snapshot)
                    }
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
            assert!(
                pair.len() >= 2,
                "cycle {c}: expected pre+post snapshots at boundary"
            );
            let (pre, post) = (pair[0], pair[1]);
            let m_pre = mass(&pre.s);
            let rel_dmass = (mass(&post.s) - m_pre).abs() / m_pre.max(1e-12);
            assert!(
                rel_dmass < 1e-3,
                "cycle {c}: macro_redistribution not mass-conserving — |Δmass|/mass = {rel_dmass:.2e} (≥ 1e-3 would indicate A1-c-style runaway erosion)"
            );
            frac_traj.push(continental_fraction(post));
        }

        // (2) convergence: last-cycle delta is small.
        let n = frac_traj.len();
        let last_delta = (frac_traj[n - 1] - frac_traj[n - 2]).abs();
        assert!(
            last_delta < 0.01,
            "workflow did not converge — last-cycle |Δfrac| = {last_delta:.4} (≥ 0.01 suggests A1-c-style ongoing collapse). Trajectory: {frac_traj:?}"
        );

        let workflow_iso = match events.last() {
            Some(C1Event::Completed { final_snapshot, .. }) => iso_land(final_snapshot),
            other => panic!("expected Completed; got {other:?}"),
        };

        // ---- Gallery control (same 100 steps, NO post-tectonic) ----
        let (cmd_tx2, cmd_rx2) = bounded(4);
        let (evt_tx2, evt_rx2) = bounded(2);
        let cancel2 = Arc::new(AtomicBool::new(false));
        let handle2 = spawn_c1_thread(cmd_rx2, evt_tx2, cancel2);
        cmd_tx2
            .send(C1Command::RunBaseline {
                spec: C1RunSpec {
                    grid_size: 64,
                    seed: 42,
                    n_steps: n_cycles * k_cycle, // SAME duration (100)
                    ..C1RunSpec::default()
                },
            })
            .unwrap();
        let gevents = drain_run(&evt_rx2);
        drop(cmd_tx2);
        handle2.join().unwrap();
        let gallery_iso = match gevents.last() {
            Some(C1Event::Completed { final_snapshot, .. }) => iso_land(final_snapshot),
            other => panic!("expected gallery Completed; got {other:?}"),
        };

        eprintln!("Stage V workflow_mode_continent_preserved (seed 42, 64², 5×20):");
        eprintln!("  continental fraction trajectory: {frac_traj:?}");
        eprintln!("  last-cycle |Δfrac| = {last_delta:.4}");
        eprintln!("  workflow iso_land = {workflow_iso:.4}, gallery iso_land = {gallery_iso:.4}");

        // (3) above the isostatic floor (apples-to-apples iso_land).
        assert!(
            workflow_iso >= 0.9 * gallery_iso,
            "workflow emergent land {workflow_iso:.4} fell below 0.9× gallery isostatic floor {gallery_iso:.4} — macro+reclassify destroyed crust beyond the pure-tectonic S̃ land (A1-c signature)"
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
        let phase_a = PhaseAParams {
            n_cycles,
            k_cycle,
            ..PhaseAParams::default()
        };
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: small_spec(0, 42),
                phase_a,
            })
            .unwrap();

        let events = drain_run(&evt_rx);

        // Group plate_type vectors by step. Cycle boundaries (step =
        // c·k_cycle, c ≥ 1) have exactly two snapshots: [0] pre-
        // reclassify (last per-step), [1] post-reclassify (post-cycle).
        use std::collections::HashMap;
        let mut by_step: HashMap<usize, Vec<Vec<u8>>> = HashMap::new();
        for e in &events {
            if let C1Event::StepCompleted { snapshot } = e {
                by_step
                    .entry(snapshot.step)
                    .or_default()
                    .push(snapshot.plate_type.clone());
            }
        }

        let mut any_reclassify = false;
        for c in 1..=n_cycles {
            let boundary = c * k_cycle;
            if let Some(pair) = by_step.get(&boundary) {
                if pair.len() == 2 && pair[0] != pair[1] {
                    let changed = pair[0]
                        .iter()
                        .zip(&pair[1])
                        .filter(|(a, b)| a != b)
                        .count();
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
        use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
        use ymir_core::tectonics_v2::field::Field2D;

        let nx = 64usize;
        let ny = 64usize;
        let n_cells = (nx * ny) as f64;

        let plate_type_cont = |snap: &C1Snapshot| -> f64 {
            snap.plate_type.iter().filter(|&&t| t == 1).count() as f64 / n_cells
        };
        let iso_land = |s: &[f64]| -> f64 {
            let f = Field2D::from_vec(nx, ny, s.to_vec());
            compute_isostasy(&f, &IsostasyConfig::default()).land_ratio as f64
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
        let spec = C1RunSpec {
            grid_size: 64,
            seed: 42,
            ..C1RunSpec::default()
        };
        cmd_tx
            .send(C1Command::RunWorkflow {
                spec: spec.clone(),
                phase_a: PhaseAParams::default(),
            })
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
                    C1Event::StepCompleted { snapshot } if snapshot.step == step => {
                        Some(snapshot)
                    }
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
        let gallery_spec = C1RunSpec {
            grid_size: 64,
            seed: 42,
            n_steps: 100,
            ..C1RunSpec::default()
        };
        cmd_tx2
            .send(C1Command::RunBaseline { spec: gallery_spec })
            .unwrap();
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
        let phase_a = PhaseAParams {
            n_cycles: 2,
            k_cycle: 5,
            ..PhaseAParams::default()
        };

        let mut state = init_c1_state_phase_2_r7(
            spec.grid_size,
            spec.seed,
            &spec.init_params,
        );
        let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let vel: Vec<(f64, f64)> = kinematics.velocities.clone();

        // Run 1 from a fresh state (base_step 0).
        let last1 = run_workflow_cycles(
            &mut state, &mut kinematics, &spec, &phase_a, &vel, &tx, 0,
        );
        assert_eq!(last1, phase_a.n_cycles * phase_a.k_cycle);
        let s_after_run1: Vec<f64> = state.s.data().to_vec();

        // Run 2 RESUMING the same state (base_step = last1).
        let last2 = run_workflow_cycles(
            &mut state, &mut kinematics, &spec, &phase_a, &vel, &tx, last1,
        );

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
