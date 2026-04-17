//! Solver thread: runs the thin viscous sheet simulation off the main thread.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::erosion::hydraulic::run_erosion;
use ymir_core::tectonics::solver::grid::StaggeredGrid;
use ymir_core::tectonics::solver::tectonics::run_tectonics;
use ymir_core::tectonics::solver::workspace::SolverWorkspace;
use ymir_core::terrain::flow::compute_flow;
use ymir_core::terrain::upscale::upscale_with_fbm;

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
                        grid_width,
                        grid_height,
                        dx,
                    } => {
                        cancel.store(false, Ordering::Relaxed);

                        let ws = workspace
                            .get_or_insert_with(|| SolverWorkspace::new(grid_width, grid_height));
                        ws.resize_if_needed(grid_width, grid_height);

                        let mut grid = StaggeredGrid::new(grid_width, grid_height, dx);
                        for j in 0..grid_height {
                            for i in 0..grid_width {
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
                            |step, total, stats, snap| {
                                let _ = tx.send(SolverEvent::Progress {
                                    step,
                                    total_steps: total,
                                    stats: stats.clone(),
                                });
                                if step % snapshot_interval == 0 || step == total - 1 {
                                    let _ = tx.send(SolverEvent::Snapshot {
                                        step,
                                        s_field: snap.s_field.clone(),
                                        plate_ids: snap.plate_ids.map(|ids| ids.to_vec()),
                                        plates: snap.plates.map(|p| p.to_vec()),
                                        boundary_types: snap.boundary_types.map(|bt| bt.to_vec()),
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
                        let final_boundary_types = workspace.as_ref().and_then(|w| {
                            w.boundary_field.as_ref().map(|bf| bf.boundary_type.clone())
                        });

                        match result {
                            Ok(()) => {
                                let _ = events_tx.send(SolverEvent::Completed {
                                    s_field: final_s,
                                    plate_ids: final_plate_ids,
                                    plates: final_plates,
                                    boundary_types: final_boundary_types,
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
                                        boundary_types: final_boundary_types,
                                        elapsed: start.elapsed(),
                                        total_steps: config.num_timesteps,
                                    });
                                } else {
                                    let _ = events_tx.send(SolverEvent::Failed { error: msg });
                                }
                            }
                        }
                    }
                    SolverCommand::SingleStep {
                        config,
                        mut plate_ctx,
                        s_field,
                        grid_width,
                        grid_height,
                        dx,
                    } => {
                        cancel.store(false, Ordering::Relaxed);

                        let ws = workspace
                            .get_or_insert_with(|| SolverWorkspace::new(grid_width, grid_height));
                        ws.resize_if_needed(grid_width, grid_height);

                        let mut grid = StaggeredGrid::new(grid_width, grid_height, dx);
                        for j in 0..grid_height {
                            for i in 0..grid_width {
                                grid.s.set(i, j, s_field.get(i, j));
                            }
                        }

                        let mut one_step = config.clone();
                        one_step.num_timesteps = 1;

                        let tx = events_tx.clone();
                        let _ = run_tectonics(
                            &one_step,
                            &mut plate_ctx,
                            &mut grid,
                            ws,
                            |_, _, stats, snap| {
                                let _ = tx.send(SolverEvent::Snapshot {
                                    step: 0,
                                    s_field: snap.s_field.clone(),
                                    plate_ids: snap.plate_ids.map(|i| i.to_vec()),
                                    plates: snap.plates.map(|p| p.to_vec()),
                                    boundary_types: snap.boundary_types.map(|bt| bt.to_vec()),
                                });
                                let _ = tx.send(SolverEvent::Progress {
                                    step: 0,
                                    total_steps: 1,
                                    stats: stats.clone(),
                                });
                                true
                            },
                        );

                        let final_boundary_types = workspace.as_ref().and_then(|w| {
                            w.boundary_field.as_ref().map(|bf| bf.boundary_type.clone())
                        });
                        let _ = events_tx.send(SolverEvent::Completed {
                            s_field: grid.s.clone(),
                            plate_ids: Some(plate_ctx.ids.clone()),
                            plates: Some(plate_ctx.plates.clone()),
                            boundary_types: final_boundary_types,
                            elapsed: std::time::Duration::ZERO,
                            total_steps: 1,
                        });
                    }
                    SolverCommand::RunFbmUpscale { coarse, sea_level, seed, config } => {
                        let start = Instant::now();
                        let result = upscale_with_fbm(&coarse, sea_level, &seed, &config);
                        let _ = events_tx.send(SolverEvent::FbmCompleted {
                            heightmap: result.heightmap,
                            slope: result.slope,
                            elapsed: start.elapsed(),
                        });
                    }
                    SolverCommand::RunErosion { heightmap, config, seed } => {
                        cancel.store(false, Ordering::Relaxed);
                        let start = Instant::now();
                        let tx = events_tx.clone();
                        let cancel_ref = cancel.clone();

                        let snapshot_interval = 3;
                        let mut batch_counter = 0usize;

                        let result =
                            run_erosion(&heightmap, &config, &seed, |completed, total, hmap| {
                                batch_counter += 1;
                                let _ = tx.send(SolverEvent::ErosionProgress { completed, total });
                                if batch_counter % snapshot_interval == 0 {
                                    let _ = tx.send(SolverEvent::ErosionSnapshot {
                                        heightmap: hmap.clone(),
                                        completed,
                                        total,
                                    });
                                }
                                !cancel_ref.load(Ordering::Relaxed)
                            });

                        let _ = events_tx.send(SolverEvent::ErosionCompleted {
                            heightmap: result.heightmap,
                            sediment: result.sediment,
                            stats: result.stats,
                            elapsed: start.elapsed(),
                        });
                    }
                    SolverCommand::RunFlowComputation { heightmap, config } => {
                        let start = Instant::now();
                        let result = compute_flow(&heightmap, &config);
                        let _ = events_tx
                            .send(SolverEvent::FlowCompleted { result, elapsed: start.elapsed() });
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
