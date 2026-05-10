//! Step 8.6 — v2 solver thread.
//!
//! Owns no Bevy state; speaks only over the `crossbeam_channel` pair
//! supplied at spawn. Phase 1 implementation: receives a
//! `V2Command::RunBaseline { spec }`, calls `run_baseline` to
//! completion, ships back a `V2Event::Completed { final_state, metrics }`.
//! Cancellation flag is wired but a no-op until Phase 5 refactors
//! `run_baseline` to accept a step callback.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};
use ymir_core::tectonics_v2::cancel as core_cancel;
use ymir_core::tectonics_v2::diagnostics::harness::{
    run_baseline_with_progress, ContinuationState,
};
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::workflow::{
    final_state_to_continuation, run_phase_a_cycle, run_phase_b, WorkflowConfig,
};

use super::build_config;
use super::commands::V2Command;
use super::events::{V2Event, V2FinalState};

/// Step 12 follow-up — RAII binder for the core thread-local cancel
/// token. Bound at the top of every run-style command branch so the
/// inner CG / Newton / erosion loops on this thread observe the
/// shared `AtomicBool` the UI flips on Stop; dropped at the end of
/// the branch so a subsequent command does not inherit a stale
/// token (and the run-thread sees `is_cancelled() == false` between
/// commands, matching the pre-Step-12 baseline).
struct CancelTokenGuard;

impl CancelTokenGuard {
    fn bind(token: Arc<AtomicBool>) -> Self {
        core_cancel::set(Some(token));
        CancelTokenGuard
    }
}

impl Drop for CancelTokenGuard {
    fn drop(&mut self) {
        core_cancel::clear();
    }
}

pub fn spawn_v2_thread(
    commands_rx: Receiver<V2Command>,
    events_tx: Sender<V2Event>,
    cancel: Arc<AtomicBool>,
) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new()
        .name("ymir-v2-bridge".into())
        .spawn(move || {
            while let Ok(cmd) = commands_rx.recv() {
                match cmd {
                    V2Command::RunBaseline { spec } => {
                        cancel.store(false, Ordering::Relaxed);
                        let _guard = CancelTokenGuard::bind(cancel.clone());
                        let _ = events_tx.send(V2Event::Started { spec: spec.clone() });

                        if spec.capture_endpoints {
                            // Best-effort dir create — harness panics
                            // on missing dirs at PNG write time.
                            let _ = std::fs::create_dir_all(&spec.output_dir);
                        }

                        let cfg = build_config::build(&spec);
                        let t0 = Instant::now();
                        // Step 8.6 follow-up — streaming run. The
                        // callback fires once per completed step,
                        // shipping a Progress event with a peek of
                        // every raster the UI might want to render.
                        // Cancellation is honoured at the next step
                        // boundary by returning `false`.
                        let cancel_ref = cancel.clone();
                        let tx = events_tx.clone();
                        let result = run_baseline_with_progress(&cfg, |progress| {
                            let peek_state = V2FinalState::from_step_progress(progress);
                            let _ = tx.send(V2Event::Progress {
                                step: progress.step,
                                total: progress.total,
                                peek_state,
                            });
                            !cancel_ref.load(Ordering::Relaxed)
                        });
                        let elapsed = t0.elapsed();

                        let final_state = V2FinalState::from_harness(&result.final_state);
                        let _ = events_tx.send(V2Event::Completed {
                            spec,
                            final_state,
                            metrics: Box::new(result.metrics),
                            elapsed,
                        });
                    }
                    V2Command::ContinueRun { spec, from_state } => {
                        cancel.store(false, Ordering::Relaxed);
                        let _guard = CancelTokenGuard::bind(cancel.clone());
                        let _ = events_tx.send(V2Event::Started { spec: spec.clone() });

                        if spec.capture_endpoints {
                            let _ = std::fs::create_dir_all(&spec.output_dir);
                        }

                        let mut cfg = build_config::build(&spec);
                        // Wire the prior run's rasters into the
                        // continuation override. The harness then
                        // skips its `init::init_s_field` /
                        // `AgeFieldState::from_initial_thickness` /
                        // `build_cratonic_factor_field` paths and
                        // uses these instead.
                        let continuation = continuation_from_final_state(&from_state);
                        cfg.continuation = Some(continuation);

                        let t0 = Instant::now();
                        let cancel_ref = cancel.clone();
                        let tx = events_tx.clone();
                        let result = run_baseline_with_progress(&cfg, |progress| {
                            let peek_state = V2FinalState::from_step_progress(progress);
                            let _ = tx.send(V2Event::Progress {
                                step: progress.step,
                                total: progress.total,
                                peek_state,
                            });
                            !cancel_ref.load(Ordering::Relaxed)
                        });
                        let elapsed = t0.elapsed();

                        let final_state = V2FinalState::from_harness(&result.final_state);
                        let _ = events_tx.send(V2Event::Completed {
                            spec,
                            final_state,
                            metrics: Box::new(result.metrics),
                            elapsed,
                        });
                    }
                    V2Command::Cancel => {
                        // Phase 1: cancellation is observed by Phase 5
                        // refactor. For now we just flip the flag so
                        // downstream tests can verify the channel
                        // plumbing.
                        cancel.store(true, Ordering::Relaxed);
                    }
                    V2Command::Shutdown => break,
                    V2Command::RunWorkflowPhaseA { spec } => {
                        cancel.store(false, Ordering::Relaxed);
                        let _guard = CancelTokenGuard::bind(cancel.clone());
                        let _ = events_tx.send(V2Event::Started { spec: spec.clone() });

                        let workflow_cfg = super::build_config::build_workflow(&spec.workflow);
                        let n_cycles = match &workflow_cfg {
                            WorkflowConfig::Enabled(p) => p.phase_a.n_cycles,
                            WorkflowConfig::Disabled => {
                                let _ = events_tx.send(V2Event::Failed {
                                    error: "RunWorkflowPhaseA requires V2WorkflowSpec::On"
                                        .into(),
                                });
                                continue;
                            }
                        };

                        let mut cfg = super::build_config::build(&spec);
                        let t0 = Instant::now();
                        let mut last_final_state: Option<V2FinalState> = None;
                        let mut cycles_run = 0;

                        for cycle_idx in 0..n_cycles {
                            if cancel.load(Ordering::Relaxed) {
                                break;
                            }
                            let cycle_output = run_phase_a_cycle(&cfg, &workflow_cfg);
                            cycles_run += 1;

                            let peek_state =
                                V2FinalState::from_harness(&cycle_output.baseline.final_state);
                            let _ = events_tx.send(V2Event::WorkflowCycleCompleted {
                                cycle_idx,
                                n_cycles,
                                peek_state: peek_state.clone(),
                                erosion_volume_removed: cycle_output.erosion_volume_removed,
                                sea_level_normalized: cycle_output.sea_level_normalized,
                                mass_drift: cycle_output.mass_drift,
                                craton_recomputation_change: cycle_output
                                    .craton_recomputation_change,
                            });

                            if cycle_idx + 1 < n_cycles {
                                cfg.continuation = Some(final_state_to_continuation(
                                    &cycle_output.baseline.final_state,
                                ));
                            }
                            last_final_state = Some(peek_state);
                        }

                        let elapsed = t0.elapsed();
                        let final_state = last_final_state.unwrap_or_else(|| V2FinalState {
                            nx: spec.grid_nx,
                            ny: spec.grid_ny,
                            dx: 0.0,
                            dy: 0.0,
                            s_field: vec![0.0; spec.grid_nx * spec.grid_ny],
                            vx: vec![0.0; spec.grid_nx * spec.grid_ny],
                            vy: vec![0.0; spec.grid_nx * spec.grid_ny],
                            strain_rate_invariant: vec![0.0; spec.grid_nx * spec.grid_ny],
                            age_field: None,
                            cratonic_factor: None,
                            plate_id: None,
                            plate_type: None,
                            boundary_flag: None,
                        });
                        let _ = events_tx.send(V2Event::WorkflowPhaseACompleted {
                            spec,
                            cycles_run,
                            final_state,
                            elapsed,
                        });
                    }
                    V2Command::ContinueWorkflowPhaseA { spec, from_state } => {
                        cancel.store(false, Ordering::Relaxed);
                        let _guard = CancelTokenGuard::bind(cancel.clone());
                        let _ = events_tx.send(V2Event::Started { spec: spec.clone() });

                        let workflow_cfg = super::build_config::build_workflow(&spec.workflow);
                        let n_cycles = match &workflow_cfg {
                            WorkflowConfig::Enabled(p) => p.phase_a.n_cycles,
                            WorkflowConfig::Disabled => {
                                let _ = events_tx.send(V2Event::Failed {
                                    error: "ContinueWorkflowPhaseA requires V2WorkflowSpec::On"
                                        .into(),
                                });
                                continue;
                            }
                        };

                        let mut cfg = super::build_config::build(&spec);
                        // Wire the prior run's rasters into cycle 1's
                        // continuation override; subsequent cycles use
                        // the orchestrator's natural chain.
                        cfg.continuation = Some(continuation_from_final_state(&from_state));

                        let t0 = Instant::now();
                        let mut last_final_state: Option<V2FinalState> = None;
                        let mut cycles_run = 0;

                        for cycle_idx in 0..n_cycles {
                            if cancel.load(Ordering::Relaxed) {
                                break;
                            }
                            let cycle_output = run_phase_a_cycle(&cfg, &workflow_cfg);
                            cycles_run += 1;

                            let peek_state =
                                V2FinalState::from_harness(&cycle_output.baseline.final_state);
                            let _ = events_tx.send(V2Event::WorkflowCycleCompleted {
                                cycle_idx,
                                n_cycles,
                                peek_state: peek_state.clone(),
                                erosion_volume_removed: cycle_output.erosion_volume_removed,
                                sea_level_normalized: cycle_output.sea_level_normalized,
                                mass_drift: cycle_output.mass_drift,
                                craton_recomputation_change: cycle_output
                                    .craton_recomputation_change,
                            });

                            if cycle_idx + 1 < n_cycles {
                                cfg.continuation = Some(final_state_to_continuation(
                                    &cycle_output.baseline.final_state,
                                ));
                            }
                            last_final_state = Some(peek_state);
                        }

                        let elapsed = t0.elapsed();
                        let final_state = last_final_state.unwrap_or(from_state);
                        let _ = events_tx.send(V2Event::WorkflowPhaseACompleted {
                            spec,
                            cycles_run,
                            final_state,
                            elapsed,
                        });
                    }
                    V2Command::RunWorkflowPhaseB { spec, from_state } => {
                        cancel.store(false, Ordering::Relaxed);
                        let _guard = CancelTokenGuard::bind(cancel.clone());
                        let _ = events_tx.send(V2Event::Started { spec: spec.clone() });

                        let workflow_cfg = super::build_config::build_workflow(&spec.workflow);
                        let t0 = Instant::now();

                        let s_field =
                            Field2D::from_vec(from_state.nx, from_state.ny, from_state.s_field);
                        match run_phase_b(&s_field, &workflow_cfg, spec.seed) {
                            Some(output) => {
                                let elapsed = t0.elapsed();
                                let hd_nx = output.heightmap.width;
                                let hd_ny = output.heightmap.height;
                                let _ = events_tx.send(V2Event::WorkflowPhaseBCompleted {
                                    spec,
                                    hd_nx,
                                    hd_ny,
                                    hd_heightmap: output.heightmap.data,
                                    sediment: output.sediment.data,
                                    grand_scale_deviation: output.grand_scale_deviation,
                                    grand_scale_deviation_p95: output
                                        .grand_scale_deviation_p95,
                                    elapsed,
                                });
                            }
                            None => {
                                let _ = events_tx.send(V2Event::Failed {
                                    error: "RunWorkflowPhaseB requires V2WorkflowSpec::On"
                                        .into(),
                                });
                            }
                        }
                    }
                }
            }
        })
        .expect("failed to spawn v2 bridge thread")
}

/// Convert a `V2FinalState` (raw `Vec<f64>` rasters) into a
/// `harness::ContinuationState` (`Field2D` for the scalar fields,
/// `Vec<f64>` for the velocity components — matches the harness's
/// internal storage). Asserts `nx * ny` consistency via
/// `Field2D::from_vec`.
fn continuation_from_final_state(s: &V2FinalState) -> ContinuationState {
    ContinuationState {
        s: Field2D::from_vec(s.nx, s.ny, s.s_field.clone()),
        vx: s.vx.clone(),
        vy: s.vy.clone(),
        age: s
            .age_field
            .as_ref()
            .map(|v| Field2D::from_vec(s.nx, s.ny, v.clone())),
        cratonic_factor: s
            .cratonic_factor
            .as_ref()
            .map(|v| Field2D::from_vec(s.nx, s.ny, v.clone())),
    }
}
