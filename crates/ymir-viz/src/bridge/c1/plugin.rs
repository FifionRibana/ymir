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

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use bevy::prelude::*;
use crossbeam_channel::{Receiver, Sender, bounded};

use ymir_core::tectonics_v2::workflow::PhaseAParams;

use super::commands::C1Command;
use super::events::{C1Event, C1RunKind};
use super::hd::{HdParams, HdPhaseRecord, HdState, PreviewShape};
use super::snapshot::C1Snapshot;
use super::spec::C1RunSpec;
use super::thread::spawn_c1_thread;

/// Run-cumulative Track D event totals — accumulated by
/// `poll_c1_events` from every `StepCompleted` snapshot's per-step
/// stats. Reset to all-zero on `Started`. Stage A bug fix (per-step
/// vs cumulative) — the UI panel needs cumulative totals because
/// rare events (accretion merges, rifting splits) have per-step
/// values of 0 almost always and a per-step display reads 0 even
/// when the run has accumulated 6 merges.
#[derive(Default, Clone, Debug)]
pub struct C1CumulativeStats {
    pub subduction_cells: usize,
    pub accretion_merges: usize,
    pub rifting_splits: usize,
    pub thinning_cells: usize,
    /// Plate ids spawned by rifting splits across the full run
    /// (appended each time `apply_rifting_split` fires). Resets
    /// at `Started`.
    pub new_plate_ids: Vec<u16>,
}

/// Lifecycle of the most recent C1 run. UI / render systems query
/// this via `C1SolverBridge::state`.
#[derive(Default, Clone)]
pub enum C1RunState {
    #[default]
    Idle,
    Running {
        spec: C1RunSpec,
        /// Pipeline kind (gallery vs workflow + cadence). Drives the
        /// UI step/cycle counter (Issue #139 Stage E2/E3).
        kind: C1RunKind,
        step: usize,
        total: usize,
        started_at: Option<Instant>,
        /// Most recent `StepCompleted` snapshot. Cleared back to
        /// `None` on `Started` of a fresh run. Required for the
        /// view-switch-during-pause behaviour.
        latest_snapshot: Option<Box<C1Snapshot>>,
        /// Run-cumulative Track D totals (Stage A bug fix).
        cumulative: C1CumulativeStats,
    },
    Completed {
        spec: C1RunSpec,
        /// Pipeline kind of the completed run (Issue #139 Stage E2).
        kind: C1RunKind,
        elapsed: Duration,
        final_snapshot: Box<C1Snapshot>,
        /// Run-cumulative Track D totals at end of run.
        cumulative: C1CumulativeStats,
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
            C1RunState::Running { latest_snapshot, .. } => latest_snapshot.as_deref(),
            C1RunState::Completed { final_snapshot, .. } => Some(final_snapshot.as_ref()),
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
    /// HD pipeline state (step b/5) — independent of the coarse tectonic
    /// `state` so the live gallery view and HD generation coexist.
    pub hd: HdState,
    /// Latest fast tectonic-shape preview (coarse continent + island verdict),
    /// for judging a seed before the HD run. `None` until the first preview.
    pub preview: Option<Arc<PreviewShape>>,
}

impl C1SolverBridge {
    /// Queue a gallery-path baseline run. Returns `Err` if the
    /// channel is full or disconnected.
    pub fn submit_run(&self, spec: C1RunSpec) -> Result<(), &'static str> {
        self.commands_tx
            .send(C1Command::RunBaseline { spec })
            .map_err(|_| "c1 bridge channel send failed")
    }

    /// Queue a workflow-path run (Issue #139 Stage E2). `phase_a`
    /// carries the calibrated Phase A cadence (default
    /// `PhaseAParams::default()`).
    pub fn submit_workflow(
        &self,
        spec: C1RunSpec,
        phase_a: PhaseAParams,
    ) -> Result<(), &'static str> {
        self.commands_tx
            .send(C1Command::RunWorkflow { spec, phase_a })
            .map_err(|_| "c1 bridge channel send failed")
    }

    /// Queue an HD production run (step b/5): the full cached chain
    /// (eroded → climate → drainage → biomes) on the worker thread.
    pub fn submit_hd(&self, spec: C1RunSpec, params: HdParams) -> Result<(), &'static str> {
        self.commands_tx
            .send(C1Command::RunHd { spec, params })
            .map_err(|_| "c1 bridge channel send failed")
    }

    /// Queue a fast tectonic-shape preview (coarse continent only) so a seed
    /// can be judged before the long HD run. Emits `PreviewReady`.
    pub fn submit_preview(&self, spec: C1RunSpec, params: HdParams) -> Result<(), &'static str> {
        self.commands_tx
            .send(C1Command::PreviewShape { spec, params })
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
            hd: HdState::Idle,
            preview: None,
        });

        app.add_systems(Update, poll_c1_events);
    }
}

fn poll_c1_events(mut bridge: ResMut<C1SolverBridge>) {
    while let Ok(event) = bridge.events_rx.try_recv() {
        match event {
            C1Event::Started { spec, kind } => {
                let total = kind.total_tectonic_steps(&spec);
                bridge.state = C1RunState::Running {
                    spec,
                    kind,
                    step: 0,
                    total,
                    started_at: Some(Instant::now()),
                    latest_snapshot: None,
                    cumulative: C1CumulativeStats::default(),
                };
            }
            C1Event::StepCompleted { snapshot } => {
                // Update the cached snapshot + accumulate Track D
                // totals. Stage A bug fix: per-step stats from
                // `snapshot.stats` are this-step-only; for rare
                // events the panel needs the cumulative figure
                // accumulated across the StepCompleted stream.
                let (spec, kind, total, started_at, mut cumulative) =
                    match std::mem::take(&mut bridge.state) {
                        C1RunState::Running {
                            spec, kind, total, started_at, cumulative, ..
                        } => (spec, kind, total, started_at, cumulative),
                        other => {
                            bridge.state = other;
                            continue;
                        }
                    };
                cumulative.subduction_cells += snapshot.stats.subduction.cells_consumed;
                cumulative.accretion_merges += snapshot.stats.accretion.merges_count;
                cumulative.rifting_splits += snapshot.stats.rifting_split.splits_count;
                cumulative.thinning_cells += snapshot.stats.rifting_thinning.cells_thinned;
                cumulative
                    .new_plate_ids
                    .extend(snapshot.stats.rifting_split.new_plate_ids_created.iter().copied());
                let step = snapshot.step;
                bridge.state = C1RunState::Running {
                    spec,
                    kind,
                    step,
                    total,
                    started_at,
                    latest_snapshot: Some(Box::new(snapshot)),
                    cumulative,
                };
            }
            C1Event::Completed { spec, final_snapshot, elapsed } => {
                // Transfer cumulative + kind from Running → Completed.
                let (kind, cumulative) = match std::mem::take(&mut bridge.state) {
                    C1RunState::Running { kind, cumulative, .. } => (kind, cumulative),
                    _ => (C1RunKind::Gallery, C1CumulativeStats::default()),
                };
                bridge.state = C1RunState::Completed {
                    spec,
                    kind,
                    elapsed,
                    final_snapshot: Box::new(final_snapshot),
                    cumulative,
                };
            }
            C1Event::Failed { error } => {
                bridge.state = C1RunState::Failed { error };
            }

            // ── HD pipeline (step b/5 + e progress) — drives `bridge.hd`. ──
            C1Event::HdStarted { params, .. } => {
                bridge.hd =
                    HdState::Running { params, current: None, progress: None, done: Vec::new() };
            }
            C1Event::HdPhaseStarted { phase } => {
                if let HdState::Running { current, progress, .. } = &mut bridge.hd {
                    *current = Some(phase);
                    *progress = None; // reset; a bar appears only once progress arrives
                }
            }
            C1Event::HdPhaseProgress { phase, done, total } => {
                if let HdState::Running { current, progress, .. } = &mut bridge.hd {
                    if *current == Some(phase) {
                        *progress = Some((done, total));
                    }
                }
            }
            C1Event::HdPhaseDone { phase, regime, elapsed } => {
                if let HdState::Running { current, progress, done, .. } = &mut bridge.hd {
                    *current = None;
                    *progress = None;
                    done.push(HdPhaseRecord { phase, regime, elapsed });
                }
            }
            C1Event::HdCompleted { result, elapsed } => {
                let done = match std::mem::take(&mut bridge.hd) {
                    HdState::Running { done, .. } => done,
                    _ => Vec::new(),
                };
                bridge.hd = HdState::Completed { result, total: elapsed, done };
            }
            C1Event::HdFailed { error } => {
                bridge.hd = HdState::Failed { error };
            }
            C1Event::PreviewReady { preview, .. } => {
                bridge.preview = Some(preview);
            }
        }
    }
}
