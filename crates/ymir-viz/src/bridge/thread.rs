//! Solver thread: runs the thin viscous sheet simulation off the main thread.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
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
                        mut plate_ctx,
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
                        let dynamic = config.dynamic_boundaries;

                        let result = run_tectonics(
                            &config,
                            &mut plate_ctx,
                            &mut grid,
                            ws,
                            |step, total, stats, s_field, plate_ids| {
                                let _ = tx.send(SolverEvent::Progress {
                                    step,
                                    total_steps: total,
                                    stats: stats.clone(),
                                });
                                if step % snapshot_interval == 0 || step == total - 1 {
                                    let _ = tx.send(SolverEvent::Snapshot {
                                        step,
                                        s_field: s_field.clone(),
                                        plate_ids: plate_ids.map(|ids| ids.to_vec()),
                                        plates: None, // sent only on completion to save bandwidth
                                    });
                                }
                                !cancel_ref.load(Ordering::Relaxed)
                            },
                        );

                        let final_s = grid.s.clone();
                        let final_plate_ids =
                            if dynamic { Some(plate_ctx.ids.clone()) } else { None };
                        let final_plates =
                            if dynamic { Some(plate_ctx.plates.clone()) } else { None };

                        match result {
                            Ok(()) => {
                                let _ = events_tx.send(SolverEvent::Completed {
                                    s_field: final_s,
                                    plate_ids: final_plate_ids,
                                    plates: final_plates,
                                    elapsed: start.elapsed(),
                                    total_steps: config.num_timesteps,
                                });
                            }
                            Err(e) => {
                                let msg = e.to_string();
                                if msg == "Simulation cancelled" {
                                    let _ = events_tx.send(SolverEvent::Completed {
                                        s_field: final_s,
                                        plate_ids: final_plate_ids,
                                        plates: final_plates,
                                        elapsed: start.elapsed(),
                                        total_steps: config.num_timesteps,
                                    });
                                } else {
                                    let _ = events_tx.send(SolverEvent::Failed { error: msg });
                                }
                            }
                        }
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
