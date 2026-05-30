//! Bevy plugin wrapping the C1 bridge thread.
//!
//! Mirrors `bridge::v2::plugin::V2BridgePlugin` shape but for the
//! C1 worker: spawns the thread, registers a `C1SolverBridge`
//! resource, polls the event channel each frame, **caches the
//! latest snapshot on the resource** so the visualization layer
//! can re-render on field-switch without waiting for a new
//! `StepCompleted` event (W4 view-switch-during-pause contract).
//!
//! ## Cached-snapshot pattern
//!
//! The poll system writes `bridge.state = Running { latest_snapshot,
//! ... }` on each `StepCompleted` (and `Completed`). The Stage E5
//! render system reads `bridge.state.latest_snapshot()` and runs
//! `field_to_rgba(&snap, field)` — pure function of (snapshot,
//! field). When the user toggles `C1Field` (e.g., S̃ → Age), the
//! render system re-fires `field_to_rgba` on the cached snapshot
//! with no channel round-trip.

use std::sync::atomic::AtomicBool;
use std::sync::Arc;
use std::time::{Duration, Instant};

use bevy::prelude::*;
use crossbeam_channel::{bounded, Receiver, Sender};

use super::commands::C1Command;
use super::events::C1Event;
use super::snapshot::C1Snapshot;
use super::spec::C1RunSpec;
use super::thread::spawn_c1_thread;

/// Lifecycle of the most recent C1 run. UI / render systems query
/// this via `C1SolverBridge::state`.
#[derive(Default, Clone)]
pub enum C1RunState {
    #[default]
    Idle,
    Running {
        spec: C1RunSpec,
        step: usize,
        total: usize,
        started_at: Option<Instant>,
        /// Most recent `StepCompleted` snapshot. Cleared back to
        /// `None` on `Started` of a fresh run. Required for the
        /// view-switch-during-pause behaviour.
        latest_snapshot: Option<Box<C1Snapshot>>,
    },
    Completed {
        spec: C1RunSpec,
        elapsed: Duration,
        final_snapshot: Box<C1Snapshot>,
    },
    Failed {
        error: String,
    },
}

impl C1RunState {
    /// Convenience accessor — returns the cached snapshot the
    /// render layer should paint, regardless of state variant.
    /// Returns `None` only when the bridge is `Idle` (nothing
    /// to display).
    pub fn latest_snapshot(&self) -> Option<&C1Snapshot> {
        match self {
            C1RunState::Idle => None,
            C1RunState::Running { latest_snapshot, .. } => {
                latest_snapshot.as_deref()
            }
            C1RunState::Completed { final_snapshot, .. } => {
                Some(final_snapshot.as_ref())
            }
            C1RunState::Failed { .. } => None,
        }
    }
}

/// Runtime resource — command sender + event receiver + cancel
/// flag + latest run state.
#[derive(Resource)]
pub struct C1SolverBridge {
    pub commands_tx: Sender<C1Command>,
    pub events_rx: Receiver<C1Event>,
    pub cancel_flag: Arc<AtomicBool>,
    pub state: C1RunState,
}

impl C1SolverBridge {
    /// Queue a baseline run. Returns `Err` if the channel is
    /// full or disconnected.
    pub fn submit_run(&self, spec: C1RunSpec) -> Result<(), &'static str> {
        self.commands_tx
            .send(C1Command::RunBaseline { spec })
            .map_err(|_| "c1 bridge channel send failed")
    }

    /// MVP Option C cancel — does NOT interrupt the current run
    /// (see `commands.rs` docstring). Sets the AtomicBool; the
    /// worker checks the flag between commands.
    pub fn request_cancel(&self) {
        use std::sync::atomic::Ordering;
        self.cancel_flag.store(true, Ordering::Relaxed);
    }
}

pub struct C1BridgePlugin;

impl Plugin for C1BridgePlugin {
    fn build(&self, app: &mut App) {
        // Per Issue #137 Stage E2 design — commands(4), events(2).
        // The tight events bound enforces backpressure when the
        // UI thread lags.
        let (cmd_tx, cmd_rx) = bounded::<C1Command>(4);
        let (evt_tx, evt_rx) = bounded::<C1Event>(2);
        let cancel = Arc::new(AtomicBool::new(false));

        spawn_c1_thread(cmd_rx, evt_tx, cancel.clone());

        app.insert_resource(C1SolverBridge {
            commands_tx: cmd_tx,
            events_rx: evt_rx,
            cancel_flag: cancel,
            state: C1RunState::Idle,
        });

        app.add_systems(Update, poll_c1_events);
    }
}

fn poll_c1_events(mut bridge: ResMut<C1SolverBridge>) {
    while let Ok(event) = bridge.events_rx.try_recv() {
        match event {
            C1Event::Started { spec } => {
                let total = spec.n_steps;
                bridge.state = C1RunState::Running {
                    spec,
                    step: 0,
                    total,
                    started_at: Some(Instant::now()),
                    latest_snapshot: None,
                };
            }
            C1Event::StepCompleted { snapshot } => {
                // Update the cached snapshot. Preserve spec /
                // started_at / total from the prior Running state.
                let (spec, total, started_at) =
                    match std::mem::take(&mut bridge.state) {
                        C1RunState::Running {
                            spec,
                            total,
                            started_at,
                            ..
                        } => (spec, total, started_at),
                        other => {
                            // Defensive: if StepCompleted arrives
                            // without a preceding Started (unlikely
                            // outside tests), restore the prior
                            // state and skip.
                            bridge.state = other;
                            continue;
                        }
                    };
                let step = snapshot.step;
                bridge.state = C1RunState::Running {
                    spec,
                    step,
                    total,
                    started_at,
                    latest_snapshot: Some(Box::new(snapshot)),
                };
            }
            C1Event::Completed {
                spec,
                final_snapshot,
                elapsed,
            } => {
                bridge.state = C1RunState::Completed {
                    spec,
                    elapsed,
                    final_snapshot: Box::new(final_snapshot),
                };
            }
            C1Event::Failed { error } => {
                bridge.state = C1RunState::Failed { error };
            }
        }
    }
}
