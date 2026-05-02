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
use ymir_core::tectonics_v2::diagnostics::harness::run_baseline_with_progress;

use super::build_config;
use super::commands::V2Command;
use super::events::{V2Event, V2FinalState};

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
                    V2Command::Cancel => {
                        // Phase 1: cancellation is observed by Phase 5
                        // refactor. For now we just flip the flag so
                        // downstream tests can verify the channel
                        // plumbing.
                        cancel.store(true, Ordering::Relaxed);
                    }
                    V2Command::Shutdown => break,
                }
            }
        })
        .expect("failed to spawn v2 bridge thread")
}
