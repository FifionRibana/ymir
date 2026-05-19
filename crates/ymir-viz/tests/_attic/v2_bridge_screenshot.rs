//! Step 8.6 Phase 6 acceptance test #4.
//!
//! Runs a tiny baseline through the bridge, then for each of the five
//! field variants writes a PNG via `save_field_png` and verifies the
//! resulting file exists, has the expected dimensions, and decodes
//! cleanly through the `image` crate. Headless — no Bevy.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{
    spawn_v2_thread, V2AgeFieldSpec, V2Command, V2CratonicSpec, V2Event, V2ForceKind,
    V2LinearSolverSpec, V2MantleSpec, V2PlateKinematicSpec, V2RunSpec, V2WorkflowSpec,
};
use ymir_viz::visualization::v2_viz::{save_field_png, screenshot_filename, V2Field};

fn quick_spec() -> V2RunSpec {
    V2RunSpec {
        seed: 7,
        grid_nx: 32,
        grid_ny: 32,
        steps: 5,
        num_plates: 4,
        continental_ratio: 0.3,
        bi: 0.15,
        br: 0.05,
        mantle: V2MantleSpec::Off,
        slab_enabled: false,
        cratonic: V2CratonicSpec::On {
            cr: 0.3,
            k_viscous: 5.0,
            b_factor: 8.0,
            smoothing_width: 0.05,
            plate_area_min: 0.10,
        },
        age_field: V2AgeFieldSpec::On {
            continental_age_init: 7.0,
            oceanic_age_init: 0.5,
        },
        linear_solver: V2LinearSolverSpec::Jacobi,
        force: V2ForceKind::Gpe,
        s_perturbation_amplitude: 0.2,
        total_time_nondim: 0.3,
        cfl_factor: 0.3,
        capture_endpoints: false,
        output_dir: std::env::temp_dir().join("ymir_v2_screenshot_test"),
        preset_label: "screenshot_test".to_string(),
        init_mode: ymir_viz::bridge_v2::V2InitModeSpec::default(),
        plate_kinematic: V2PlateKinematicSpec::Zero,
        workflow: V2WorkflowSpec::Off,
    }
}

#[test]
fn v2_bridge_screenshot() {
    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(16);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    cmd_tx
        .send(V2Command::RunBaseline { spec: quick_spec() })
        .expect("command send");

    let deadline = Instant::now() + Duration::from_secs(60);
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
            Ok(V2Event::Failed { error }) => panic!("bridge failed: {}", error),
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("channel disconnected"),
        }
    }
    let final_state = state.expect("Completed within deadline");

    let scratch = std::env::temp_dir().join(format!(
        "ymir_v2_screenshot_phase6_{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    std::fs::create_dir_all(&scratch).expect("scratch dir");

    for &field in V2Field::ALL {
        let filename = screenshot_filename("test_preset", field);
        // D6 filename hygiene: includes "test_preset" + a field tag
        // + .png.
        assert!(filename.contains("test_preset"), "filename: {}", filename);
        assert!(filename.ends_with(".png"), "filename: {}", filename);

        let path = scratch.join(&filename);
        save_field_png(&final_state, field, &path).unwrap_or_else(|e| {
            panic!("save_field_png({:?}) failed: {}", field, e)
        });

        // Re-decode through the image crate and verify dimensions.
        let img = image::open(&path).unwrap_or_else(|e| {
            panic!("re-decode {} failed: {}", path.display(), e)
        });
        assert_eq!(img.width(), 32, "{:?} png width", field);
        assert_eq!(img.height(), 32, "{:?} png height", field);

        // PNG file is non-empty.
        let meta = std::fs::metadata(&path).expect("metadata");
        assert!(meta.len() > 0, "{:?} png empty", field);
    }

    // Cleanup the scratch dir to keep the temp area tidy.
    let _ = std::fs::remove_dir_all(&scratch);

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");
}
