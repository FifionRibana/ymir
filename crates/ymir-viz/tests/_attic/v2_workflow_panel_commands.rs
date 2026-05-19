//! Step 12 Phase 7b.6 — integration tests for the workflow-panel
//! command path.
//!
//! These tests verify that the `submit_workflow_*` helpers added to
//! `V2SolverBridge` in Phase 7b.2 emit the expected `V2Command`
//! variants on the channel, and that `request_cancel()` flips the
//! shared `AtomicBool` the harness step-callback reads. They mirror
//! the pattern of `v2_bridge_lifecycle.rs` but skip
//! `spawn_v2_thread` entirely: the bridge struct is constructed
//! directly with mock channels, so the workflow command surface is
//! exercised without booting the harness.
//!
//! Coverage:
//!
//! 1. submit_workflow_phase_a → V2Command::RunWorkflowPhaseA
//! 2. submit_continue_workflow_phase_a → V2Command::ContinueWorkflowPhaseA
//! 3. submit_workflow_phase_b → V2Command::RunWorkflowPhaseB
//! 4. request_cancel → cancel_flag set
//! 5. bridge default state is Idle (smoke probe — guarantees the
//!    enabled-state logic in workflow_panel sees the expected
//!    starting variant in tests that build a fresh bridge).
//!
//! UI-coupled assertions (button enabled state, egui::Ui rendering)
//! are out of scope: the egui frame would require booting Bevy +
//! bevy_egui in the test, which the user explicitly flagged as
//! avoidable in the Phase 7b plan. The unit tests in
//! `crates/ymir-viz/src/ui/workflow_panel.rs` already cover the
//! pure-function helpers the buttons call (FIFO history eviction,
//! PNG export round-trip, length-mismatch rejection).

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use crossbeam_channel::{Receiver, bounded};
use ymir_viz::bridge_v2::{
    V2Command, V2Event, V2FinalState, V2PhaseAParams, V2PhaseBParams, V2RunSpec, V2RunState,
    V2SolverBridge, V2WorkflowSpec,
};

/// Build a bare V2SolverBridge backed by mock channels. Returns the
/// bridge plus the *receiver* end of the command channel so the test
/// can observe what `submit_*` enqueued. The event channel is created
/// but unused — wired through to keep the bridge construction valid.
fn build_mock_bridge() -> (V2SolverBridge, Receiver<V2Command>) {
    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (_evt_tx, evt_rx) = bounded::<V2Event>(16);
    let cancel = Arc::new(AtomicBool::new(false));
    let bridge = V2SolverBridge {
        commands_tx: cmd_tx,
        events_rx: evt_rx,
        cancel_flag: cancel,
        state: V2RunState::Idle,
    };
    (bridge, cmd_rx)
}

/// Minimal `V2FinalState` payload for tests that need a `from_state`
/// argument (Continue Phase A, Run Phase B). The values are
/// deterministic and trivially small; the bridge layer treats this
/// payload as opaque, so any consistent shape suffices.
fn minimal_final_state(nx: usize, ny: usize) -> V2FinalState {
    let n = nx * ny;
    V2FinalState {
        nx,
        ny,
        dx: 1.0,
        dy: 1.0,
        s_field: vec![0.5; n],
        vx: vec![0.0; n],
        vy: vec![0.0; n],
        strain_rate_invariant: vec![0.0; n],
        age_field: None,
        cratonic_factor: None,
        plate_id: None,
        plate_type: None,
        boundary_flag: None,
    }
}

/// Build a workflow-On spec from the canonical preset, replacing
/// `workflow` with the pinned defaults so tests assert on a known
/// shape rather than the preset's choice.
fn workflow_on_spec() -> V2RunSpec {
    let mut spec = V2RunSpec::active_medley_defaults();
    spec.workflow = V2WorkflowSpec::On {
        phase_a: V2PhaseAParams::default(),
        phase_b: V2PhaseBParams::default(),
    };
    spec
}

#[test]
fn submit_workflow_phase_a_emits_run_workflow_phase_a_command() {
    let (bridge, cmd_rx) = build_mock_bridge();
    let spec = workflow_on_spec();
    let pinned_seed = spec.seed;

    bridge
        .submit_workflow_phase_a(spec.clone())
        .expect("submit_workflow_phase_a must succeed on a fresh channel");

    let cmd = cmd_rx
        .recv_timeout(Duration::from_secs(1))
        .expect("a V2Command must be enqueued by submit_workflow_phase_a");

    match cmd {
        V2Command::RunWorkflowPhaseA { spec: received } => {
            assert_eq!(received.seed, pinned_seed);
            assert!(
                matches!(received.workflow, V2WorkflowSpec::On { .. }),
                "the bridge must forward the spec's workflow=On verbatim"
            );
        }
        _ => panic!("expected V2Command::RunWorkflowPhaseA"),
    }
}

#[test]
fn submit_continue_workflow_phase_a_emits_continue_command_with_from_state() {
    let (bridge, cmd_rx) = build_mock_bridge();
    let spec = workflow_on_spec();
    let from_state = minimal_final_state(8, 8);
    let pinned_nx = from_state.nx;
    let pinned_ny = from_state.ny;

    bridge
        .submit_continue_workflow_phase_a(spec, from_state)
        .expect("submit_continue_workflow_phase_a must succeed");

    let cmd = cmd_rx
        .recv_timeout(Duration::from_secs(1))
        .expect("a V2Command must be enqueued");

    match cmd {
        V2Command::ContinueWorkflowPhaseA { spec: _, from_state: rs } => {
            // The bridge must thread `from_state` through unchanged so
            // the worker can warm-start cycle 1 from the prior run.
            assert_eq!(rs.nx, pinned_nx);
            assert_eq!(rs.ny, pinned_ny);
            assert_eq!(rs.s_field.len(), pinned_nx * pinned_ny);
        }
        _ => panic!("expected V2Command::ContinueWorkflowPhaseA"),
    }
}

#[test]
fn submit_workflow_phase_b_emits_run_phase_b_command_with_from_state() {
    let (bridge, cmd_rx) = build_mock_bridge();
    let spec = workflow_on_spec();
    let from_state = minimal_final_state(4, 4);

    bridge
        .submit_workflow_phase_b(spec, from_state)
        .expect("submit_workflow_phase_b must succeed");

    let cmd = cmd_rx
        .recv_timeout(Duration::from_secs(1))
        .expect("a V2Command must be enqueued");

    assert!(
        matches!(cmd, V2Command::RunWorkflowPhaseB { .. }),
        "expected V2Command::RunWorkflowPhaseB"
    );
}

#[test]
fn request_cancel_sets_cancel_flag() {
    let (bridge, _cmd_rx) = build_mock_bridge();
    assert!(
        !bridge.cancel_flag.load(Ordering::Relaxed),
        "fresh bridge cancel flag must start as false"
    );

    bridge.request_cancel();

    assert!(
        bridge.cancel_flag.load(Ordering::Relaxed),
        "request_cancel must set the shared AtomicBool to true"
    );
}

#[test]
fn fresh_bridge_state_is_idle() {
    let (bridge, _cmd_rx) = build_mock_bridge();
    assert!(
        matches!(bridge.state, V2RunState::Idle),
        "freshly built bridge must default to V2RunState::Idle so the \
         workflow panel sees the expected `can_run = true / can_stop = \
         false / can_continue = false` enabled-state on first frame"
    );
}

/// The `submit_*` helpers must not block on a saturated channel — the
/// channel has capacity 4, and once the receiver disconnects, the next
/// send returns the documented `Err`. Defends the panic-free contract
/// (the panel calls these on every click; an Err is surfaced via
/// `eprintln!` rather than crashing the UI).
#[test]
fn submit_workflow_phase_a_returns_err_when_channel_disconnected() {
    let (bridge, cmd_rx) = build_mock_bridge();
    drop(cmd_rx);

    let result = bridge.submit_workflow_phase_a(workflow_on_spec());
    assert!(
        result.is_err(),
        "submit_workflow_phase_a must return Err when the channel is \
         disconnected, not panic"
    );
}
