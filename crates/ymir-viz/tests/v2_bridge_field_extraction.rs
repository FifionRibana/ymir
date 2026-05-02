//! Step 8.6 Phase 5 acceptance test #3.
//!
//! Verifies the bridge correctly extracts each of the five raster
//! fields (D5) plus the auxiliary plate / boundary tags. Runs a small
//! end-to-end baseline with cratonic + age field on so every Option
//! lands `Some`. Headless — no Bevy.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{
    spawn_v2_thread, V2AgeFieldSpec, V2Command, V2CratonicSpec, V2Event, V2ForceKind,
    V2LinearSolverSpec, V2MantleSpec, V2RunSpec,
};

fn full_mechanism_spec() -> V2RunSpec {
    V2RunSpec {
        seed: 42,
        grid_nx: 32,
        grid_ny: 32,
        steps: 10,
        num_plates: 6,
        continental_ratio: 0.4,
        bi: 0.15,
        br: 0.05,
        // Mantle on so the run produces non-trivial velocity (else
        // |v| would be near machine zero everywhere and the log
        // colormap test below loses signal).
        mantle: V2MantleSpec::On {
            mf: 1.0,
            coupling: 1.0,
            num_modes: 4,
            seed: 7,
            evolution_rate: 0.0,
        },
        slab_enabled: false,
        // Cratonic + age on so the corresponding Options land Some.
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
        total_time_nondim: 0.6,
        cfl_factor: 0.3,
        capture_endpoints: false,
        output_dir: std::env::temp_dir().join("ymir_v2_field_extraction"),
        preset_label: "field_extraction".to_string(),
        init_mode: ymir_viz::bridge_v2::V2InitModeSpec::default(),
    }
}

#[test]
fn v2_bridge_field_extraction() {
    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(16);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    cmd_tx
        .send(V2Command::RunBaseline { spec: full_mechanism_spec() })
        .expect("command send");

    let deadline = Instant::now() + Duration::from_secs(120);
    let mut completed = None;
    while Instant::now() < deadline {
        match evt_rx.recv_timeout(Duration::from_millis(500)) {
            Ok(V2Event::Started { .. }) => {}
            Ok(V2Event::Progress { .. }) => {}
            Ok(V2Event::Completed { final_state, .. }) => {
                completed = Some(final_state);
                break;
            }
            Ok(V2Event::Failed { error }) => panic!("bridge failed: {}", error),
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("channel disconnected"),
        }
    }
    let final_state = completed.expect("Completed event within deadline");

    let n = 32 * 32;

    // ── Required fields (always populated) ─────────────────────────
    assert_eq!(final_state.s_field.len(), n, "S̃ length");
    assert_eq!(final_state.vx.len(), n, "vx length");
    assert_eq!(final_state.vy.len(), n, "vy length");
    assert_eq!(
        final_state.strain_rate_invariant.len(),
        n,
        "strain_rate_invariant length"
    );
    assert!(
        final_state.s_field.iter().all(|v| v.is_finite()),
        "S̃ has non-finite values"
    );
    assert!(
        final_state.vx.iter().all(|v| v.is_finite()),
        "vx has non-finite values"
    );
    assert!(
        final_state.vy.iter().all(|v| v.is_finite()),
        "vy has non-finite values"
    );
    assert!(
        final_state.strain_rate_invariant.iter().all(|v| v.is_finite() && *v >= 0.0),
        "ε̇_II must be finite and non-negative"
    );

    // ── Optional fields, populated by spec choice ─────────────────
    let age = final_state.age_field.as_ref().expect("age requested but None");
    assert_eq!(age.len(), n);
    assert!(age.iter().all(|v| v.is_finite() && *v >= 0.0), "A must be ≥ 0");

    let crat = final_state
        .cratonic_factor
        .as_ref()
        .expect("cratonic requested but None");
    assert_eq!(crat.len(), n);
    assert!(
        crat.iter().all(|v| v.is_finite() && (0.0..=1.0).contains(v)),
        "cratonic_factor must lie in [0, 1]"
    );

    // ── Auxiliary tags (boundary enabled in spec ⇒ all Some) ──────
    let plate_id = final_state.plate_id.as_ref().expect("plate_id Some");
    let plate_type = final_state.plate_type.as_ref().expect("plate_type Some");
    let boundary_flag = final_state.boundary_flag.as_ref().expect("boundary_flag Some");
    assert_eq!(plate_id.len(), n);
    assert_eq!(plate_type.len(), n);
    assert_eq!(boundary_flag.len(), n);
    // plate_type encodes Oceanic=0 / Continental=1; boundary_flag in [0, 4].
    assert!(plate_type.iter().all(|t| *t <= 1), "plate_type encoded 0/1");
    assert!(boundary_flag.iter().all(|t| *t <= 4), "boundary_flag encoded 0..=4");

    // ── Magnitude check on |v| (compute as in v2_viz) ─────────────
    let v_mag: Vec<f64> = (0..n)
        .map(|k| (final_state.vx[k].powi(2) + final_state.vy[k].powi(2)).sqrt())
        .collect();
    let v_max = v_mag.iter().copied().fold(0.0_f64, f64::max);
    // Mantle is enabled — peak|v| should be well above machine zero,
    // demonstrating the field has signal for the colormap to render.
    assert!(v_max > 1e-6, "|v|_max = {} (mantle on, expected > 1e-6)", v_max);

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");
}
