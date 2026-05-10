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
    /// Step 12 Phase 7 — Phase A multi-cycle loop has completed.
    /// Carries the final state for chaining into a subsequent
    /// `V2Command::RunWorkflowPhaseB`. Distinct from `Completed`
    /// because no `Metrics` struct is produced (the workflow's
    /// per-cycle metrics are emitted as `WorkflowCycleCompleted`
    /// events; aggregating them into a single `Metrics` is a Phase
    /// 7b refinement).
    WorkflowPhaseACompleted {
        spec: V2RunSpec,
        cycles_run: usize,
        final_state: Box<V2FinalState>,
        elapsed: Duration,
    },
    /// Step 12 Phase 7 — Phase B HD finalization has completed.
    /// Holds the HD heightmap (flat `Vec<f32>` row-major,
    /// `hd_nx × hd_ny`) for download / visualization.
    /// `grand_scale_deviation_p95` is the formal D5 acceptance
    /// metric (see `PhaseBOutput` for context).
    WorkflowPhaseBCompleted {
        spec: V2RunSpec,
        hd_nx: usize,
        hd_ny: usize,
        hd_heightmap: Vec<f32>,
        grand_scale_deviation_p95: f64,
        elapsed: Duration,
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

    /// Step 8.6 follow-up — queue a "continue from prior final
    /// state" run. The harness uses `from_state` as the initial S̃ /
    /// vx / vy / age / cratonic_factor instead of computing init
    /// from scratch. The voronoi-relevant fields of `spec` (seed,
    /// num_plates, continental_ratio, grid dims) should match the
    /// source run for the override to make physical sense.
    pub fn submit_continue(
        &self,
        spec: V2RunSpec,
        from_state: V2FinalState,
    ) -> Result<(), &'static str> {
        self.commands_tx
            .send(V2Command::ContinueRun { spec, from_state })
            .map_err(|_| "v2 bridge channel send failed")
    }

    /// Signal the running run (if any) to abort. Sets the shared
    /// `AtomicBool` directly rather than queueing a
    /// `V2Command::Cancel`: the bridge thread is blocked inside
    /// `run_baseline_with_progress` while a run is in flight, so a
    /// channel command would not reach the dispatch loop until the
    /// run completed naturally. The harness step-callback reads the
    /// same `AtomicBool` and returns `false` on the next step
    /// boundary, so this lands as a graceful cancellation within
    /// one step (≈ 5–25 s at 32–64² mantle-on).
    pub fn request_cancel(&self) {
        use std::sync::atomic::Ordering;
        self.cancel_flag.store(true, Ordering::Relaxed);
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
            V2Event::WorkflowCycleCompleted {
                cycle_idx,
                n_cycles,
                peek_state,
                ..
            } => {
                // Reuse Running for in-flight workflow Phase A. step
                // tracks cycles, total tracks n_cycles. Per-cycle
                // metric payload (mass_drift, etc.) is dropped here
                // for Phase 7a; the dashboard wires them in 7b.
                let (spec, started_at) = match std::mem::take(&mut bridge.state) {
                    V2RunState::Running { spec, started_at, .. } => (spec, started_at),
                    other => {
                        // First cycle: a Started event populated Running
                        // already. If we land here from another state,
                        // synthesise minimal context so the UI doesn't
                        // crash; the next event will refresh.
                        bridge.state = other;
                        continue;
                    }
                };
                bridge.state = V2RunState::Running {
                    spec,
                    step: cycle_idx + 1,
                    total: n_cycles,
                    started_at,
                    peek_state: Some(Box::new(peek_state)),
                };
            }
            V2Event::WorkflowPhaseACompleted {
                spec,
                cycles_run,
                final_state,
                elapsed,
            } => {
                bridge.state = V2RunState::WorkflowPhaseACompleted {
                    spec,
                    cycles_run,
                    final_state: Box::new(final_state),
                    elapsed,
                };
            }
            V2Event::WorkflowPhaseBCompleted {
                spec,
                hd_nx,
                hd_ny,
                hd_heightmap,
                sediment: _,
                grand_scale_deviation: _,
                grand_scale_deviation_p95,
                elapsed,
            } => {
                bridge.state = V2RunState::WorkflowPhaseBCompleted {
                    spec,
                    hd_nx,
                    hd_ny,
                    hd_heightmap,
                    grand_scale_deviation_p95,
                    elapsed,
                };
            }
        }
    }
}
