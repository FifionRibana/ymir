//! Step 8.6 Phase 1 acceptance test #1.
//!
//! Spawns the v2 bridge thread, sends one tiny `RunBaseline` command,
//! verifies it receives `Started` then `Completed` with a populated
//! final-state snapshot, then shuts the thread down cleanly. No Bevy
//! involvement — the bridge is consumed as a plain crossbeam channel
//! pair so this test is fully headless and CI-friendly.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{
    spawn_v2_thread, V2AgeFieldSpec, V2Command, V2CratonicSpec, V2Event, V2ForceKind,
    V2LinearSolverSpec, V2MantleSpec, V2PlateKinematicSpec, V2RunSpec, V2WorkflowSpec,
};

/// Tiny spec — 32² × 5 steps with mantle off (Step 7 shape).
/// Picked to keep wallclock under ~10 s on the development laptop;
/// goal is lifecycle correctness, not numerical fidelity.
fn tiny_spec() -> V2RunSpec {
    V2RunSpec {
        seed: 7,
        grid_nx: 32,
        grid_ny: 32,
        steps: 5,
        num_plates: 4,
        continental_ratio: 0.3,
        bi: 0.15,
        br: 0.05,
        mantle: V2MantleSpec::Off,
        slab_enabled: false,
        cratonic: V2CratonicSpec::Off,
        age_field: V2AgeFieldSpec::On {
            continental_age_init: 7.0,
            oceanic_age_init: 0.5,
        },
        linear_solver: V2LinearSolverSpec::Jacobi,
        force: V2ForceKind::Gpe,
        s_perturbation_amplitude: 0.2,
        total_time_nondim: 0.3,
        cfl_factor: 0.3,
        capture_endpoints: false,
        output_dir: std::env::temp_dir().join("ymir_v2_lifecycle_test"),
        preset_label: "lifecycle_test".to_string(),
        init_mode: ymir_viz::bridge_v2::V2InitModeSpec::default(),
        plate_kinematic: V2PlateKinematicSpec::Zero,
        workflow: V2WorkflowSpec::Off,
    }
}

#[test]
fn v2_bridge_lifecycle() {
    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(16);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel.clone());

    cmd_tx
        .send(V2Command::RunBaseline { spec: tiny_spec() })
        .expect("command send should succeed");

    // Drain events with a generous-but-bounded timeout. The tiny run
    // is ~1-2 s release; the 60 s budget protects against CI
    // pathological scheduling without papering over a hung thread.
    let deadline = Instant::now() + Duration::from_secs(60);
    let mut got_started = false;
    let mut got_completed = false;

    while Instant::now() < deadline && !got_completed {
        match evt_rx.recv_timeout(Duration::from_millis(500)) {
            Ok(V2Event::Started { spec }) => {
                assert_eq!(spec.grid_nx, 32);
                assert_eq!(spec.grid_ny, 32);
                got_started = true;
            }
            Ok(V2Event::Progress { step, total, peek_state }) => {
                // Step 8.6 follow-up — verify the progress stream is
                // wired correctly: step index in [0, total] (the
                // harness emits step=0 with initial state then
                // step=1..total per completed step), matching total,
                // peek snapshot of the right size.
                assert!(step <= total, "step out of range: {}/{}", step, total);
                assert_eq!(peek_state.nx, 32);
                assert_eq!(peek_state.ny, 32);
            }
            Ok(V2Event::Completed { final_state, metrics, elapsed, .. }) => {
                assert!(got_started, "Started must precede Completed");
                assert_eq!(final_state.nx, 32);
                assert_eq!(final_state.ny, 32);
                assert_eq!(final_state.s_field.len(), 32 * 32);
                assert_eq!(final_state.vx.len(), 32 * 32);
                assert_eq!(final_state.vy.len(), 32 * 32);
                assert_eq!(final_state.strain_rate_invariant.len(), 32 * 32);
                assert!(
                    final_state.age_field.is_some(),
                    "age_field requested in spec, must be populated"
                );
                assert!(
                    final_state.cratonic_factor.is_none(),
                    "cratonic disabled in spec, must be None"
                );
                assert!(
                    final_state.plate_id.is_some(),
                    "boundary enabled (Voronoi closed), plate_id must be present"
                );
                // Sanity: S̃ field is finite everywhere (no NaN/Inf).
                assert!(
                    final_state.s_field.iter().all(|v| v.is_finite()),
                    "S̃ contains non-finite values"
                );
                // Sanity: metrics reflect the requested step count.
                assert_eq!(metrics.steps, 5);
                assert!(metrics.cg_iter_mean > 0.0, "CG iters mean must be > 0");
                println!(
                    "v2_bridge_lifecycle: completed in {:.2}s, cg_iter_mean = {:.1}",
                    elapsed.as_secs_f64(),
                    metrics.cg_iter_mean
                );
                got_completed = true;
            }
            Ok(V2Event::WorkflowCycleCompleted { .. })
            | Ok(V2Event::WorkflowPhaseACompleted { .. })
            | Ok(V2Event::WorkflowPhaseBCompleted { .. }) => {}
            Ok(V2Event::Failed { error }) => {
                panic!("bridge reported failure: {}", error);
            }
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("bridge channel disconnected"),
        }
    }

    assert!(got_started, "did not receive V2Event::Started within deadline");
    assert!(got_completed, "did not receive V2Event::Completed within deadline");

    // Clean shutdown.
    cmd_tx.send(V2Command::Shutdown).expect("shutdown send should succeed");
    handle.join().expect("bridge thread should join cleanly");
}
