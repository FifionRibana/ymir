//! Bevy plugin that manages the solver bridge and polls events each frame.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::Duration;

use bevy::prelude::*;
use crossbeam_channel::{Receiver, Sender, bounded};

use ymir_core::tectonics::solver::workspace::StepStats;

use super::commands::SolverCommand;
use super::events::SolverEvent;
use super::thread::spawn_solver_thread;
use crate::visualization::render::TerrainDisplay;

/// Current state of the solver.
#[derive(Default, Debug, Clone)]
pub enum SolverState {
    #[default]
    Idle,
    Running {
        step: usize,
        total_steps: usize,
        stats: Option<StepStats>,
    },
    Completed {
        elapsed: Duration,
    },
    Failed {
        error: String,
    },
}

/// Resource holding the channels to communicate with the solver thread.
#[derive(Resource)]
pub struct SolverBridge {
    pub commands_tx: Sender<SolverCommand>,
    pub events_rx: Receiver<SolverEvent>,
    pub state: SolverState,
    pub cancel_flag: Arc<AtomicBool>,
}

pub struct TectonicsBridgePlugin;

impl Plugin for TectonicsBridgePlugin {
    fn build(&self, app: &mut App) {
        let (cmd_tx, cmd_rx) = bounded::<SolverCommand>(4);
        let (evt_tx, evt_rx) = bounded::<SolverEvent>(256);
        let cancel = Arc::new(AtomicBool::new(false));

        spawn_solver_thread(cmd_rx, evt_tx, cancel.clone());

        app.insert_resource(SolverBridge {
            commands_tx: cmd_tx,
            events_rx: evt_rx,
            state: SolverState::Idle,
            cancel_flag: cancel,
        });

        app.add_systems(Update, poll_solver_events);
    }
}

fn poll_solver_events(
    mut bridge: ResMut<SolverBridge>,
    mut terrain_display: ResMut<TerrainDisplay>,
) {
    while let Ok(event) = bridge.events_rx.try_recv() {
        match event {
            SolverEvent::Progress { step, total_steps, stats } => {
                bridge.state = SolverState::Running { step, total_steps, stats: Some(stats) };
            }
            SolverEvent::Snapshot { s_field, .. } => {
                terrain_display.update_field(s_field);
            }
            SolverEvent::Completed { s_field, elapsed, .. } => {
                terrain_display.update_field(s_field);
                bridge.state = SolverState::Completed { elapsed };
            }
            SolverEvent::Failed { error } => {
                bridge.state = SolverState::Failed { error };
            }
        }
    }
}
