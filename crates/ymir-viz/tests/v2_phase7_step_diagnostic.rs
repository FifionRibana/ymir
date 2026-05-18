//! Step 8.6 Phase 7 follow-up — per-step diagnostic strip.
//!
//! Marked `#[ignore]`. Runs ONE preset at a small grid and dumps every
//! field to PNG at every (or every Nth) step boundary, so the reviewer
//! can inspect the simulation evolution frame-by-frame from inside
//! Claude (or any image viewer). Distinct from
//! `v2_phase7_screenshot_gallery` which only captures the *final*
//! state per preset.
//!
//! Run explicitly:
//!
//! ```text
//! # default — active_medley, 32² × 30 steps, every step
//! cargo test --release -p ymir-viz --test v2_phase7_step_diagnostic \
//!     -- --ignored --nocapture --jobs 1
//!
//! # custom — single_continent, 50 steps, capture every 5 steps
//! YMIR_DIAG_PRESET=single_continent YMIR_DIAG_STEPS=50 \
//! YMIR_DIAG_INTERVAL=5 cargo test --release -p ymir-viz \
//!     --test v2_phase7_step_diagnostic -- --ignored --nocapture --jobs 1
//! ```
//!
//! Output:
//!   `docs/reports/step8_6_phase7_gallery/diagnostic/<preset>/`
//!   filled with `step_NNNN_<field>.png` (one per
//!   `(captured_step, field)` pair). The directory is gitignored —
//!   regenerate on demand.
//!
//! Acceptance: this is *not* an automated visual-coherence test
//! (Step 8.6 D7 explicitly says reviewer judgment). The test passes
//! as long as the bridge produces a `Completed` event and at least
//! one Progress event with a peek_state of the expected dimensions.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{presets, spawn_v2_thread, V2Command, V2Event, V2FinalState};
use ymir_viz::visualization::v2_viz::{save_field_png, V2Field};

fn env_str(name: &str, default: &str) -> String {
    std::env::var(name).unwrap_or_else(|_| default.to_string())
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn field_stats(buf: &[f64]) -> (f64, f64, f64) {
    if buf.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0_f64;
    for &v in buf {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        sum += v;
    }
    (min, max, sum / buf.len() as f64)
}

fn field_stats_vmag(vx: &[f64], vy: &[f64]) -> (f64, f64, f64) {
    let n = vx.len().min(vy.len());
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0_f64;
    for k in 0..n {
        let m = (vx[k].powi(2) + vy[k].powi(2)).sqrt();
        if m < min {
            min = m;
        }
        if m > max {
            max = m;
        }
        sum += m;
    }
    (min, max, sum / n as f64)
}

fn save_all_fields(state: &V2FinalState, dir: &std::path::Path, step: usize) {
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
        let path = dir.join(format!("step_{:04}_{}.png", step, tag));
        save_field_png(state, field, &path)
            .unwrap_or_else(|e| panic!("save_field_png({:?}, step={}) failed: {}", field, step, e));
    }
}

#[test]
#[ignore]
fn v2_phase7_step_diagnostic() {
    let preset_name = env_str("YMIR_DIAG_PRESET", "active_medley");
    let steps = env_usize("YMIR_DIAG_STEPS", 30);
    let interval = env_usize("YMIR_DIAG_INTERVAL", 1).max(1);
    let grid = env_usize("YMIR_DIAG_GRID", 32);

    let mut spec = presets::load(&preset_name)
        .unwrap_or_else(|e| panic!("preset '{}' load failed: {}", preset_name, e));
    // Diagnostic-friendly overrides — keep the preset's physics knobs
    // (Bi, Br, Mf, Cr, K, B_factor, mantle/cratonic/age toggles) but
    // shrink the run to fit a frame-by-frame inspection budget.
    spec.grid_nx = grid;
    spec.grid_ny = grid;
    spec.steps = steps;
    // Total time scales with steps so per-step `dt` stays close to the
    // preset's original value (preset default = total_time / steps).
    spec.total_time_nondim = (steps as f64) * 0.06;
    spec.capture_endpoints = false;

    let out_root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step8_6_phase7_gallery/diagnostic")
        .join(&preset_name);
    // Wipe the per-preset diagnostic dir on each run so stale frames
    // from a previous step count don't linger. The parent
    // `diagnostic/` is gitignored so this is safe.
    let _ = std::fs::remove_dir_all(&out_root);
    std::fs::create_dir_all(&out_root).expect("create diagnostic dir");

    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(256);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    println!(
        "[diag] preset='{}' grid={}² steps={} interval={} → {}",
        preset_name,
        grid,
        steps,
        interval,
        out_root.display()
    );
    let t0 = Instant::now();
    cmd_tx
        .send(V2Command::RunBaseline { spec })
        .expect("command send");

    // Per-preset deadline — same generous bound as the gallery test
    // (mantle-on presets at 32² take ~7 s/step at full CPU, so 30
    // steps ≈ 210 s; the deadline gives margin for CPU contention).
    let deadline = Instant::now() + Duration::from_secs(900);
    let mut completed = false;
    let mut frames_saved = 0usize;
    while Instant::now() < deadline {
        match evt_rx.recv_timeout(Duration::from_millis(500)) {
            Ok(V2Event::Started { .. }) => {
                println!("[diag] started");
            }
            Ok(V2Event::Progress { step, total, peek_state }) => {
                let is_capture_step = step == 0
                    || step == total
                    || (interval > 0 && step % interval == 0);
                if is_capture_step {
                    save_all_fields(&peek_state, &out_root, step);
                    frames_saved += 1;
                    let (smin, smax, smean) = field_stats(&peek_state.strain_rate_invariant);
                    let (vmin, vmax, vmean) = field_stats_vmag(
                        &peek_state.vx,
                        &peek_state.vy,
                    );
                    println!(
                        "[diag] step {:>3}/{}: strain[min={:.3e} max={:.3e} mean={:.3e}] |v|[min={:.3e} max={:.3e} mean={:.3e}]",
                        step, total, smin, smax, smean, vmin, vmax, vmean
                    );
                }
            }
            Ok(V2Event::Completed { final_state, elapsed, metrics, .. }) => {
                // Final state may differ from the last Progress peek
                // (strain_rate_invariant is computed only at the end);
                // save it under a sentinel `step_final_<field>.png`.
                let final_dir = out_root.join("final");
                std::fs::create_dir_all(&final_dir).expect("create final subdir");
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
                    let path = final_dir.join(format!("{}.png", tag));
                    save_field_png(&final_state, field, &path).unwrap_or_else(|e| {
                        panic!("save_field_png(final, {:?}) failed: {}", field, e)
                    });
                }
                let (fsmin, fsmax, fsmean) = field_stats(&final_state.strain_rate_invariant);
                let yielding_frac = metrics
                    .newton
                    .as_ref()
                    .and_then(|n| n.yielding_cell_fraction_max);
                println!(
                    "[diag] FINAL strain[min={:.3e} max={:.3e} mean={:.3e}] yielding_cell_fraction_max={}",
                    fsmin,
                    fsmax,
                    fsmean,
                    match yielding_frac {
                        Some(f) => format!("{:.4}", f),
                        None => "<none>".into(),
                    }
                );
                println!(
                    "[diag] completed in {:.1}s — peak|v|={:.2e} CG mean={:.0}",
                    elapsed.as_secs_f64(),
                    metrics.vmax_peak,
                    metrics.cg_iter_mean
                );
                completed = true;
                break;
            }
            Ok(V2Event::Failed { error }) => panic!("bridge failed: {}", error),
            // Workflow events never fire in this RunBaseline-only
            // diagnostic; harmlessly ignored if the bridge surfaces
            // them in the future (kept as a wildcard to avoid this
            // test re-breaking on every new V2Event variant).
            Ok(_) => {}
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("channel disconnected"),
        }
    }
    assert!(completed, "diagnostic run did not complete within deadline");
    assert!(frames_saved >= 1, "no Progress frames captured");

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");

    let total_elapsed = t0.elapsed();
    println!(
        "[diag] {} frames saved in {:.1}s — open {}",
        frames_saved,
        total_elapsed.as_secs_f64(),
        out_root.display()
    );
}
