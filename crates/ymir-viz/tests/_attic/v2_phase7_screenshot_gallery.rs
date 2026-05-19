//! Step 8.6 Phase 7 prep — automated screenshot gallery.
//!
//! Marked `#[ignore]` because it is intentionally slower (runs every
//! preset to completion at 32² × 50 steps with the mantle-on
//! configurations, ~30-60 s per preset on the development laptop).
//! Run explicitly:
//!
//! ```text
//! cargo test --release -p ymir-viz --test v2_phase7_screenshot_gallery -- --ignored --nocapture
//! ```
//!
//! For each registered preset, the test:
//!   1. Loads the preset via `bridge_v2::presets::load`
//!   2. Reduces grid to 32² and steps to 50 to keep wallclock low
//!      (geographic coherence is established within ~50 active steps;
//!      the full 100 steps adds detail but doesn't change the
//!      shape that the reviewer is judging at this checkpoint)
//!   3. Runs the baseline through the bridge
//!   4. Saves all five fields to
//!      `docs/reports/step8_6_phase7_gallery/<preset>/<field>.png`
//!
//! The output directory becomes the artifact the reviewer inspects
//! alongside the live UI to satisfy D7. This is *not* an automated
//! visual-coherence test — D7 explicitly says reviewer judgment.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{presets, spawn_v2_thread, V2Command, V2Event};
use ymir_viz::visualization::v2_viz::{save_field_png, V2Field};

#[test]
#[ignore]
fn v2_phase7_screenshot_gallery() {
    let gallery_root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step8_6_phase7_gallery");
    std::fs::create_dir_all(&gallery_root).expect("create gallery root");

    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(64);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    for name in presets::list() {
        let mut spec = presets::load(name).expect(name);
        // Phase 7 prep — uniform 32² × 30 step shrink so the gallery
        // builds in under ~5 minutes total rather than ~30. Reviewer
        // can re-run any specific preset at full size from the UI.
        // 30 steps is enough to reach a representative visual; the
        // mantle-on presets (~7 s/step at 32²) would otherwise blow
        // the per-preset deadline.
        spec.grid_nx = 32;
        spec.grid_ny = 32;
        spec.steps = 30;
        spec.total_time_nondim = 1.8;
        spec.capture_endpoints = false;
        let preset_dir = gallery_root.join(name);
        std::fs::create_dir_all(&preset_dir).expect("create preset subdir");
        spec.output_dir = preset_dir.clone();

        println!(
            "[gallery] running preset '{}' ({}x{} × {} steps)",
            name, spec.grid_nx, spec.grid_ny, spec.steps
        );
        let t0 = Instant::now();
        cmd_tx
            .send(V2Command::RunBaseline { spec })
            .expect("command send");

        let deadline = Instant::now() + Duration::from_secs(600);
        let mut state = None;
        while Instant::now() < deadline {
            match evt_rx.recv_timeout(Duration::from_millis(500)) {
                Ok(V2Event::Started { .. }) => {}
                Ok(V2Event::Progress { .. }) => {}
                Ok(V2Event::Completed { final_state, .. }) => {
                    state = Some(final_state);
                    break;
                }
                Ok(V2Event::WorkflowCycleCompleted { .. })
                | Ok(V2Event::WorkflowPhaseACompleted { .. })
                | Ok(V2Event::WorkflowPhaseBCompleted { .. }) => {}
                Ok(V2Event::Failed { error }) => panic!("preset '{}' failed: {}", name, error),
                Err(RecvTimeoutError::Timeout) => continue,
                Err(RecvTimeoutError::Disconnected) => panic!("channel disconnected"),
            }
        }
        let final_state =
            state.unwrap_or_else(|| panic!("preset '{}' did not complete in time", name));
        let elapsed = t0.elapsed();

        for &field in V2Field::ALL {
            let tag = match field {
                V2Field::SThickness => "s",
                V2Field::Altitude => "altitude",
                V2Field::Age => "age",
                V2Field::Cratonic => "cratonic",
                V2Field::StrainRate => "strain",
                V2Field::VelocityMagnitude => "vmag",
                V2Field::Slope => "slope",
            };
            let path = preset_dir.join(format!("{}.png", tag));
            save_field_png(&final_state, field, &path).unwrap_or_else(|e| {
                panic!("save_field_png({}, {:?}) failed: {}", name, field, e)
            });
        }
        println!(
            "[gallery] preset '{}' done in {:.1}s, 5 PNGs at {}",
            name,
            elapsed.as_secs_f64(),
            preset_dir.display()
        );
    }

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");

    println!(
        "[gallery] complete — open {} for the reviewer checkpoint",
        gallery_root.display()
    );
}
