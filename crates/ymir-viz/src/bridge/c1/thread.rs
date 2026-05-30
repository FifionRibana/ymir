//! Worker thread driving the C1 time loop.
//!
//! ## Lifecycle
//!
//! 1. `spawn_c1_thread` returns a `JoinHandle`. The thread blocks
//!    on `commands_rx.recv()` between runs.
//! 2. On `C1Command::RunBaseline { spec }`:
//!    a. Reset the cancel flag (allow future cancels).
//!    b. Send `C1Event::Started`.
//!    c. Build state + kinematics + config from `spec`.
//!    d. **Clone `kinematics.velocities` pre-run** into a local
//!       `initial_velocities: Vec<(f64, f64)>` — captured by the
//!       per-step closure (cannot access `kinematics` directly
//!       because `run_with_closures` borrows it `&mut`; see
//!       Stage E2 W7 Q-E1.3).
//!    e. Emit a pre-run cycle-0 `StepCompleted` so the UI can
//!       render the init state before the first step fires.
//!    f. Call `run_with_closures`. The per-step closure emits a
//!       `StepCompleted { snapshot }` event each step.
//!    g. Send `C1Event::Completed { spec, final_snapshot,
//!       elapsed }`.
//! 3. On `C1Command::Cancel`: set the cancel flag (no effect on
//!    current run per Q-E1.3 Option C MVP).
//!
//! ## Worker owns no Bevy state (W4 global)
//!
//! All Bevy-side data (resources, sprites, image handles) lives
//! in `plugin.rs` + the UI/render systems. The worker uses only
//! crossbeam channels + `ymir-core` types. This isolation is
//! shared with `bridge::v2::thread`.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender};

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1TimeLoopConfig};

use super::commands::C1Command;
use super::events::C1Event;
use super::snapshot::C1Snapshot;

pub fn spawn_c1_thread(
    commands_rx: Receiver<C1Command>,
    events_tx: Sender<C1Event>,
    cancel: Arc<AtomicBool>,
) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new()
        .name("ymir-c1-bridge".into())
        .spawn(move || {
            while let Ok(cmd) = commands_rx.recv() {
                match cmd {
                    C1Command::RunBaseline { spec } => {
                        // Reset cancel before each run (Option C
                        // MVP — the flag is only sampled between
                        // runs, never during).
                        cancel.store(false, Ordering::Relaxed);

                        let _ = events_tx.send(C1Event::Started {
                            spec: spec.clone(),
                        });

                        let mut state = init_c1_state_phase_2_r7(
                            spec.grid_size,
                            spec.seed,
                            &spec.init_params,
                        );
                        let mut kinematics =
                            PlateKinematics::preset_phase_1_1(state.num_plates);

                        // Pre-run snapshot of kinematics velocities
                        // for the per-step closure (see module
                        // docstring rationale).
                        let initial_velocities: Vec<(f64, f64)> =
                            kinematics.velocities.clone();

                        let config = C1TimeLoopConfig {
                            n_steps: spec.n_steps,
                            dx: 1.0 / spec.grid_size as f64,
                            dy: 1.0 / spec.grid_size as f64,
                            iso_config: IsostasyConfig::default(),
                            drainage_max_distance: spec.drainage_max_distance,
                        };

                        // Cycle-0 pre-run snapshot — lets the UI
                        // paint the init state before stepping
                        // begins. Uses `step = 0` (matches the
                        // first `on_step` callback's step index;
                        // n_steps + 1 events total per run).
                        let cycle_0 = C1Snapshot::from_state(
                            0,
                            &state,
                            &initial_velocities,
                        );
                        let _ = events_tx.send(C1Event::StepCompleted {
                            snapshot: cycle_0,
                        });

                        let t0 = Instant::now();
                        let tx = events_tx.clone();
                        let velocities_for_closure = initial_velocities.clone();
                        run_with_closures(
                            &mut state,
                            &mut kinematics,
                            &config,
                            &spec.closures,
                            |step, state| {
                                let snapshot = C1Snapshot::from_state(
                                    step,
                                    state,
                                    &velocities_for_closure,
                                );
                                // `send` blocks when the bounded
                                // events channel is full — this is
                                // the backpressure pause semantics
                                // (W7 Stage E2).
                                let _ =
                                    tx.send(C1Event::StepCompleted { snapshot });
                            },
                        );
                        let elapsed = t0.elapsed();

                        // Final snapshot — re-extract for the
                        // Completed convenience payload.
                        let final_snapshot = C1Snapshot::from_state(
                            spec.n_steps.saturating_sub(1),
                            &state,
                            &initial_velocities,
                        );
                        let _ = events_tx.send(C1Event::Completed {
                            spec,
                            final_snapshot,
                            elapsed,
                        });
                    }
                    C1Command::Cancel => {
                        // MVP Option C — set the flag; takes
                        // effect on the next RunBaseline only.
                        cancel.store(true, Ordering::Relaxed);
                    }
                }
            }
        })
        .expect("ymir-c1-bridge thread spawn")
}
