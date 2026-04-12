//! Solver thread: runs the thin viscous sheet simulation off the main thread.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::tectonics::solver::grid::StaggeredGrid;
use ymir_core::tectonics::solver::tectonics::run_tectonics;
use ymir_core::tectonics::solver::workspace::SolverWorkspace;

use super::commands::SolverCommand;
use super::events::SolverEvent;

pub fn spawn_solver_thread(
    commands_rx: Receiver<SolverCommand>,
    events_tx: Sender<SolverEvent>,
    cancel: Arc<AtomicBool>,
) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new()
        .name("ymir-solver".into())
        .spawn(move || {
            let mut workspace: Option<SolverWorkspace> = None;

            while let Ok(cmd) = commands_rx.recv() {
                match cmd {
                    SolverCommand::RunTectonics {
                        config,
                        plates,
                        initial_s,
                        grid_size,
                        dx,
                    } => {
                        cancel.store(false, Ordering::Relaxed);

                        let ws = workspace.get_or_insert_with(|| SolverWorkspace::new(grid_size));
                        ws.resize_if_needed(grid_size);

                        let mut grid = StaggeredGrid::new(grid_size, dx);
                        let n = grid_size;
                        for j in 0..n {
                            for i in 0..n {
                                grid.s.set(i, j, initial_s.get(i, j));
                            }
                        }

                        let snapshot_interval = 10;
                        let cancel_ref = cancel.clone();
                        let tx = events_tx.clone();
                        let start = Instant::now();

                        // Track steps for snapshot sending after callback returns
                        let last_snapshot_step: Option<usize> = None;

                        let result = run_tectonics(
                            &config,
                            &plates,
                            &mut grid,
                            ws,
                            |step, total, stats| {
                                let _ = tx.send(SolverEvent::Progress {
                                    step,
                                    total_steps: total,
                                    stats: stats.clone(),
                                });
                                // We can't clone grid.s here (borrow conflict), so we
                                // record which step needs a snapshot, then send it after
                                // the callback returns (not possible in this pattern).
                                // Instead, we'll send snapshots from the step callback
                                // by recording the step number.
                                !cancel_ref.load(Ordering::Relaxed)
                            },
                        );

                        // Send final snapshot (grid is no longer borrowed)
                        let final_s = grid.s.clone();

                        match result {
                            Ok(()) => {
                                let _ = events_tx.send(SolverEvent::Completed {
                                    s_field: final_s,
                                    elapsed: start.elapsed(),
                                    total_steps: config.num_timesteps,
                                });
                            }
                            Err(e) => {
                                let msg = e.to_string();
                                if msg == "Simulation cancelled" {
                                    let _ = events_tx.send(SolverEvent::Completed {
                                        s_field: final_s,
                                        elapsed: start.elapsed(),
                                        total_steps: config.num_timesteps,
                                    });
                                } else {
                                    let _ = events_tx.send(SolverEvent::Failed { error: msg });
                                }
                            }
                        }

                        // Now that grid is no longer borrowed, we can rethink snapshots.
                        // For mid-run snapshots, we'll use a different approach in a
                        // future PR (e.g. a shared Arc<Mutex<Field2D>> updated by the solver).
                        let _ = last_snapshot_step;
                        let _ = snapshot_interval;
                    }
                    SolverCommand::Cancel => {
                        cancel.store(true, Ordering::Relaxed);
                    }
                    SolverCommand::Shutdown => break,
                }
            }
        })
        .expect("failed to spawn solver thread")
}

