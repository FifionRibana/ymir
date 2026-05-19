//! Step 8.6 Phase 8e — full export/import roundtrip across the bridge.
//!
//! Spawns a tiny v2 run, exports the resulting `Completed` state to a
//! JSON snapshot on disk, then reloads it as a `V2RunSnapshot` and
//! asserts:
//!
//! 1. Every raster field round-trips byte-identically (s_field, vx,
//!    vy, strain_rate_invariant, age_field, cratonic_factor, plate_id,
//!    plate_type, boundary_flag).
//! 2. `format_version`, `spec.preset_label`, and the snapshot
//!    timestamp survive serialisation.
//! 3. `V2RunSnapshot::load` rejects a snapshot whose `format_version`
//!    differs from `SNAPSHOT_FORMAT_VERSION` (forward-compat probe).

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{
    spawn_v2_thread, V2AgeFieldSpec, V2Command, V2CratonicSpec, V2Event, V2ForceKind,
    V2InitModeSpec, V2LinearSolverSpec, V2MantleSpec, V2PlateKinematicSpec, V2RunSnapshot,
    V2RunSpec, V2WorkflowSpec, SNAPSHOT_FORMAT_VERSION,
};

fn tiny_spec() -> V2RunSpec {
    V2RunSpec {
        seed: 7,
        grid_nx: 16,
        grid_ny: 16,
        steps: 3,
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
        total_time_nondim: 0.18,
        cfl_factor: 0.3,
        capture_endpoints: false,
        output_dir: std::env::temp_dir().join("ymir_v2_export_import_test"),
        preset_label: "export_import_test".to_string(),
        init_mode: V2InitModeSpec::default(),
        plate_kinematic: V2PlateKinematicSpec::Zero,
        workflow: V2WorkflowSpec::Off,
    }
}

#[test]
fn v2_bridge_export_import_roundtrip() {
    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(64);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    let spec = tiny_spec();
    cmd_tx
        .send(V2Command::RunBaseline { spec: spec.clone() })
        .expect("command send");

    let deadline = Instant::now() + Duration::from_secs(120);
    let mut completed: Option<(V2RunSpec, ymir_viz::bridge_v2::V2FinalState, _, Duration)> = None;
    while Instant::now() < deadline {
        match evt_rx.recv_timeout(Duration::from_millis(500)) {
            Ok(V2Event::Started { .. }) | Ok(V2Event::Progress { .. }) => {}
            Ok(V2Event::Completed { spec, final_state, metrics, elapsed }) => {
                completed = Some((spec, final_state, metrics, elapsed));
                break;
            }
            Ok(V2Event::WorkflowCycleCompleted { .. })
            | Ok(V2Event::WorkflowPhaseACompleted { .. })
            | Ok(V2Event::WorkflowPhaseBCompleted { .. }) => {}
            Ok(V2Event::Failed { error }) => panic!("bridge failed: {}", error),
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("disconnected"),
        }
    }
    let (out_spec, final_state, metrics, elapsed) =
        completed.expect("run did not complete within deadline");

    let snap_dir = std::env::temp_dir().join("ymir_v2_export_import_test_snaps");
    let snap_path = snap_dir.join("roundtrip.json");
    let original = V2RunSnapshot::new(out_spec, final_state, &metrics, elapsed);
    original.save(&snap_path).expect("save snapshot");

    let loaded = V2RunSnapshot::load(&snap_path).expect("load snapshot");

    assert_eq!(loaded.format_version, SNAPSHOT_FORMAT_VERSION);
    assert_eq!(loaded.format_version, original.format_version);
    assert_eq!(loaded.exported_at, original.exported_at);
    assert!((loaded.elapsed_seconds - original.elapsed_seconds).abs() < 1e-9);
    assert_eq!(loaded.spec.preset_label, original.spec.preset_label);

    assert_eq!(loaded.final_state.nx, original.final_state.nx);
    assert_eq!(loaded.final_state.ny, original.final_state.ny);

    // f64 round-trip via JSON's `dtoa::shortest` is "round-trip
    // correct" in the IEEE-754 sense for ~all values, but a tiny
    // fraction of f64s produce a 1-ULP drift after parse-serialize
    // (serde_json doesn't bit-pin floats — switching to a binary
    // format like postcard/CBOR would fix that). For the v2 viz
    // snapshot the contract is "same render", not "same bits", and
    // 1 ULP on S̃/v/ε̇_II is invisible at colormap (8-bit) resolution.
    // Compare with a tolerance well below the colormap bin width
    // (1 / 256 ≈ 4e-3) but tight enough to catch real bugs.
    fn vec_close(a: &[f64], b: &[f64], eps: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| (x - y).abs() < eps)
    }
    fn opt_vec_close(a: &Option<Vec<f64>>, b: &Option<Vec<f64>>, eps: f64) -> bool {
        match (a, b) {
            (Some(av), Some(bv)) => vec_close(av, bv, eps),
            (None, None) => true,
            _ => false,
        }
    }
    const EPS: f64 = 1e-12;
    assert!(vec_close(&loaded.final_state.s_field, &original.final_state.s_field, EPS));
    assert!(vec_close(&loaded.final_state.vx, &original.final_state.vx, EPS));
    assert!(vec_close(&loaded.final_state.vy, &original.final_state.vy, EPS));
    assert!(vec_close(
        &loaded.final_state.strain_rate_invariant,
        &original.final_state.strain_rate_invariant,
        EPS,
    ));
    assert!(opt_vec_close(
        &loaded.final_state.age_field,
        &original.final_state.age_field,
        EPS,
    ));
    assert!(opt_vec_close(
        &loaded.final_state.cratonic_factor,
        &original.final_state.cratonic_factor,
        EPS,
    ));
    // u16 / u8 fields are integers — exact equality applies.
    assert_eq!(loaded.final_state.plate_id, original.final_state.plate_id);
    assert_eq!(loaded.final_state.plate_type, original.final_state.plate_type);
    assert_eq!(
        loaded.final_state.boundary_flag,
        original.final_state.boundary_flag
    );

    // Scalar metrics survive (sample probe). Scalars round-trip
    // through dtoa with the same ULP caveat as the rasters.
    assert!(
        (loaded.scalar_metrics.cg_iter_mean - original.scalar_metrics.cg_iter_mean).abs()
            < EPS
    );
    assert!(
        (loaded.scalar_metrics.vmax_peak - original.scalar_metrics.vmax_peak).abs() < EPS
    );

    // Forward-compat probe: load fails on a snapshot with an
    // unknown format_version. Reuse the loaded snapshot, bump the
    // version, save under a new path, expect Err on load.
    let mut bumped = loaded.clone();
    bumped.format_version = SNAPSHOT_FORMAT_VERSION + 1;
    let bumped_path = snap_dir.join("future_version.json");
    let bumped_json = serde_json::to_string(&bumped).expect("serialize");
    std::fs::write(&bumped_path, bumped_json).expect("write");
    let err = V2RunSnapshot::load(&bumped_path).expect_err("must reject");
    assert!(
        err.to_string().contains("format_version"),
        "unexpected error: {}",
        err
    );

    // Cleanup. Failures here are non-fatal — the test passed.
    let _ = std::fs::remove_file(&snap_path);
    let _ = std::fs::remove_file(&bumped_path);

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");
}
