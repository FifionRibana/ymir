//! Step 8.6 — Bevy plugin wrapping the v2 bridge thread.
//!
//! Mirrors `bridge::plugin::TectonicsBridgePlugin` but for the v2
//! solver: spawns a dedicated worker thread, registers a
//! `V2SolverBridge` resource, polls the event channel each frame.
//!
//! Phase 2 scope: thread + resource + event-poll system. UI to
//! actually dispatch commands lives in Phase 3 (parameter_panel_v2);
//! the plugin alone is silent until something pushes a
//! `V2Command::RunBaseline` onto its sender.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use bevy::prelude::*;
use crossbeam_channel::{Receiver, Sender, bounded};

use super::commands::V2Command;
use super::events::{V2Event, V2FinalState};
use super::snapshot::V2ScalarMetrics;
use super::spec::V2RunSpec;
use super::thread::spawn_v2_thread;
use ymir_core::tectonics_v2::diagnostics::metrics::Metrics;

/// Lifecycle of the most recent v2 run, surfaced to the UI through
/// `V2SolverBridge::state`.
///
/// `Running` carries the latest progress payload (step counter, total
/// step count, the bridge thread's start-of-run `Instant` for elapsed
/// wallclock, and the most recent peek-state snapshot for mid-run
/// rendering) so the UI can paint a progress bar / counter / live
/// sprite without having to peek into a separate channel.
#[derive(Default, Clone)]
pub enum V2RunState {
    #[default]
    Idle,
    Running {
        spec: V2RunSpec,
        step: usize,
        total: usize,
        started_at: Option<Instant>,
        peek_state: Option<Box<V2FinalState>>,
    },
    Completed {
        spec: V2RunSpec,
        elapsed: Duration,
        metrics: Box<Metrics>,
        final_state: Box<V2FinalState>,
    },
    /// Step 8.6 Phase 8e — a snapshot loaded from disk (no live solver
    /// run took place this session). Carries only the scalar metrics
    /// the dashboard renders; full `Metrics` is unavailable for
    /// imports because exports drop the heavy histograms /
    /// per-step series. Renders identically to `Completed` from the
    /// sprite's point of view (the `final_state` is the same shape).
    Imported {
        spec: V2RunSpec,
        elapsed: Duration,
        exported_at: String,
        scalar_metrics: V2ScalarMetrics,
        final_state: Box<V2FinalState>,
    },
    Failed {
        error: String,
    },
}

/// Runtime resource — channels into and events back from the v2 worker
/// thread, plus the latest run state for the UI to query.
#[derive(Resource)]
pub struct V2SolverBridge {
    pub commands_tx: Sender<V2Command>,
    pub events_rx: Receiver<V2Event>,
    pub cancel_flag: Arc<AtomicBool>,
    pub state: V2RunState,
}

impl V2SolverBridge {
    /// Convenience: queue a run command. Returns `Err` if the worker
    /// channel is full or disconnected (the latter shouldn't happen
    /// during normal operation).
    pub fn submit_run(&self, spec: V2RunSpec) -> Result<(), &'static str> {
        self.commands_tx
            .send(V2Command::RunBaseline { spec })
            .map_err(|_| "v2 bridge channel send failed")
    }

    /// Convenience: signal the running run (if any) to abort. No-op
    /// during Phase 1-4 since `run_baseline` does not honour the flag
    /// yet (Phase 5 step-callback refactor).
    pub fn request_cancel(&self) {
        let _ = self.commands_tx.send(V2Command::Cancel);
    }
}

pub struct V2BridgePlugin;

impl Plugin for V2BridgePlugin {
    fn build(&self, app: &mut App) {
        let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
        // Step 8.6 follow-up — Progress events fire once per step at
        // simulation cadence (one per ~7-30 s on 64² mantle-on, sub-second
        // on small grids). 256 entries is a generous buffer for
        // small-grid bursts; the Bevy poll system drains at ~60 Hz so
        // the channel never approaches saturation in practice.
        let (evt_tx, evt_rx) = bounded::<V2Event>(256);
        let cancel = Arc::new(AtomicBool::new(false));

        spawn_v2_thread(cmd_rx, evt_tx, cancel.clone());

        app.insert_resource(V2SolverBridge {
            commands_tx: cmd_tx,
            events_rx: evt_rx,
            cancel_flag: cancel,
            state: V2RunState::Idle,
        });

        app.add_systems(Update, poll_v2_events);
    }
}

fn poll_v2_events(mut bridge: ResMut<V2SolverBridge>) {
    while let Ok(event) = bridge.events_rx.try_recv() {
        match event {
            V2Event::Started { spec } => {
                let total = spec.steps;
                bridge.state = V2RunState::Running {
                    spec,
                    step: 0,
                    total,
                    started_at: Some(Instant::now()),
                    peek_state: None,
                };
            }
            V2Event::Progress { step, total, peek_state } => {
                // Preserve the existing `spec` and `started_at` from
                // the prior Running state; only update step/total and
                // the peek snapshot. If we somehow receive a Progress
                // before a Started (shouldn't happen — bridge thread
                // emits Started before invoking the harness callback),
                // fall back to a zero-spec stub.
                let (spec, started_at) = match std::mem::take(&mut bridge.state) {
                    V2RunState::Running { spec, started_at, .. } => (spec, started_at),
                    other => {
                        bridge.state = other;
                        continue;
                    }
                };
                bridge.state = V2RunState::Running {
                    spec,
                    step,
                    total,
                    started_at,
                    peek_state: Some(Box::new(peek_state)),
                };
            }
            V2Event::Completed { spec, final_state, metrics, elapsed } => {
                bridge.state = V2RunState::Completed {
                    spec,
                    elapsed,
                    metrics,
                    final_state: Box::new(final_state),
                };
            }
            V2Event::Failed { error } => {
                bridge.state = V2RunState::Failed { error };
            }
        }
    }
}
