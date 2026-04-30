//! Step 8.6 Phase 8f — equilibrium analysis on `active_medley`.
//!
//! Two back-to-back runs of `active_medley` with the new default
//! `InitMode::Uniform` (Phase 8a) at the same `dt_target`:
//!
//! - **A** — 100 steps, `total_time_nondim = 6.0`.
//! - **B** — 200 steps, `total_time_nondim = 12.0` (i.e. same dt, twice
//!   the simulated time).
//!
//! Compare the two `V2ScalarMetrics` blocks for `peak|v|`, mass drift,
//! yielding-cell fraction, peak yielding in craton, cratonic-cell
//! fraction, and mass-conservation residual. If `|Δ| / |A| < 5%` on
//! every metric the system is at equilibrium by step 100 (the
//! milestone's working assumption); otherwise the report records the
//! actual gaps so step counts in downstream phases can be revised.
//!
//! Marked `#[ignore]`. Wallclock budget on a development laptop:
//! 32² mantle-on ≈ 25 min total (10 + 15 for the two runs);
//! 64² mantle-on ≈ 60–90 min. Override the grid via env var:
//!
//! ```text
//! YMIR_PHASE8F_GRID=64 cargo test --release -p ymir-viz \
//!     --test v2_phase8f_equilibrium \
//!     -- --ignored --nocapture --jobs 1
//! ```
//!
//! Output:
//!   `docs/reports/step8_6_phase8f_equilibrium/active_medley_<grid>sq.md`
//! plus the same summary printed to stdout.

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use crossbeam_channel::{bounded, RecvTimeoutError};
use ymir_viz::bridge_v2::{
    presets, spawn_v2_thread, V2Command, V2Event, V2RunSpec, V2ScalarMetrics,
};

const EQUILIBRIUM_TOLERANCE: f64 = 0.05; // 5 %

struct RunOutcome {
    label: String,
    metrics: V2ScalarMetrics,
    elapsed: Duration,
}

fn env_grid() -> usize {
    std::env::var("YMIR_PHASE8F_GRID")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(32)
}

fn run_one(spec: V2RunSpec, label: &str, deadline_secs: u64) -> RunOutcome {
    let (cmd_tx, cmd_rx) = bounded::<V2Command>(4);
    let (evt_tx, evt_rx) = bounded::<V2Event>(256);
    let cancel = Arc::new(AtomicBool::new(false));
    let handle = spawn_v2_thread(cmd_rx, evt_tx, cancel);

    println!(
        "[phase8f] launching '{}' — {}² × {} steps, t_max={}",
        label, spec.grid_nx, spec.steps, spec.total_time_nondim
    );
    let t0 = Instant::now();
    cmd_tx
        .send(V2Command::RunBaseline { spec })
        .expect("send command");

    let deadline = Instant::now() + Duration::from_secs(deadline_secs);
    let mut outcome: Option<RunOutcome> = None;
    while Instant::now() < deadline {
        match evt_rx.recv_timeout(Duration::from_millis(500)) {
            Ok(V2Event::Started { .. }) => {}
            Ok(V2Event::Progress { step, total, .. }) => {
                if step % 25 == 0 || step == total {
                    println!(
                        "[phase8f]   '{}' step {}/{} ({:.1}s)",
                        label,
                        step,
                        total,
                        t0.elapsed().as_secs_f64()
                    );
                }
            }
            Ok(V2Event::Completed { metrics, elapsed, .. }) => {
                outcome = Some(RunOutcome {
                    label: label.to_string(),
                    metrics: V2ScalarMetrics::from_metrics(&metrics),
                    elapsed,
                });
                break;
            }
            Ok(V2Event::Failed { error }) => panic!("'{}' failed: {}", label, error),
            Err(RecvTimeoutError::Timeout) => continue,
            Err(RecvTimeoutError::Disconnected) => panic!("'{}' channel disconnected", label),
        }
    }

    cmd_tx.send(V2Command::Shutdown).expect("shutdown");
    handle.join().expect("thread join clean");

    outcome.unwrap_or_else(|| panic!("'{}' did not complete within deadline", label))
}

#[derive(Clone, Copy)]
struct CompareRow {
    label: &'static str,
    a: Option<f64>,
    b: Option<f64>,
}

impl CompareRow {
    fn drift_relative(&self) -> Option<f64> {
        match (self.a, self.b) {
            (Some(a), Some(b)) if a.abs() > 1e-15 => Some((b - a).abs() / a.abs()),
            (Some(a), Some(b)) if a.abs() <= 1e-15 && b.abs() <= 1e-15 => Some(0.0),
            _ => None,
        }
    }
    fn within_tolerance(&self, tol: f64) -> Option<bool> {
        self.drift_relative().map(|d| d < tol)
    }
}

fn compare_rows(a: &V2ScalarMetrics, b: &V2ScalarMetrics) -> Vec<CompareRow> {
    vec![
        CompareRow { label: "peak |v|", a: Some(a.vmax_peak), b: Some(b.vmax_peak) },
        CompareRow {
            label: "mass drift |relative|",
            a: Some(a.mass_drift_relative.abs()),
            b: Some(b.mass_drift_relative.abs()),
        },
        CompareRow {
            label: "CG iters mean",
            a: Some(a.cg_iter_mean),
            b: Some(b.cg_iter_mean),
        },
        CompareRow {
            label: "yielding cells max",
            a: a.yielding_cell_fraction_max,
            b: b.yielding_cell_fraction_max,
        },
        CompareRow {
            label: "yielding in craton (peak)",
            a: a.peak_yielding_in_craton,
            b: b.peak_yielding_in_craton,
        },
        CompareRow {
            label: "cratonic cell fraction",
            a: a.cratonic_cell_fraction,
            b: b.cratonic_cell_fraction,
        },
        CompareRow {
            label: "mass conservation residual",
            a: a.mass_conservation_residual,
            b: b.mass_conservation_residual,
        },
    ]
}

fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{:.4e}", x),
        None => "—".to_string(),
    }
}

fn fmt_pct(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{:.2}%", x * 100.0),
        None => "—".to_string(),
    }
}

fn fmt_verdict(v: Option<bool>) -> &'static str {
    match v {
        Some(true) => "yes",
        Some(false) => "**no**",
        None => "—",
    }
}

fn build_report(grid: usize, a: &RunOutcome, b: &RunOutcome) -> String {
    let rows = compare_rows(&a.metrics, &b.metrics);
    let all_pass = rows
        .iter()
        .all(|r| r.within_tolerance(EQUILIBRIUM_TOLERANCE).unwrap_or(true));
    let any_actionable_fail = rows
        .iter()
        .any(|r| r.within_tolerance(EQUILIBRIUM_TOLERANCE) == Some(false));

    let mut s = String::new();
    s.push_str("# Step 8.6 Phase 8f — equilibrium analysis (active_medley)\n\n");
    s.push_str(&format!(
        "Grid: {grid}² · Init mode: `Uniform`. Run A: 100 steps, t_max=6. \
         Run B: 200 steps, t_max=12 (same dt, twice the simulated time).\n\n"
    ));
    s.push_str(&format!(
        "Wallclock A: {:.1}s — B: {:.1}s.\n\n",
        a.elapsed.as_secs_f64(),
        b.elapsed.as_secs_f64()
    ));
    s.push_str("| Metric | A (step 100) | B (step 200) | |Δ| / |A| | <5%? |\n");
    s.push_str("|---|---|---|---|---|\n");
    for r in &rows {
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} |\n",
            r.label,
            fmt_opt(r.a),
            fmt_opt(r.b),
            fmt_pct(r.drift_relative()),
            fmt_verdict(r.within_tolerance(EQUILIBRIUM_TOLERANCE)),
        ));
    }
    s.push_str("\n");

    s.push_str("## Verdict\n\n");
    if all_pass && !any_actionable_fail {
        s.push_str(&format!(
            "All comparable metrics are within {:.0}% relative drift between t=6 and \
             t=12 — `active_medley` reaches equilibrium by **step 100** at this \
             grid. Downstream phases can keep the canonical 100-step budget.\n",
            EQUILIBRIUM_TOLERANCE * 100.0
        ));
    } else {
        s.push_str(&format!(
            "**Equilibrium not reached at step 100.** At least one metric drifts \
             more than {:.0}% between t=6 and t=12 (see the table). The largest \
             drifts:\n\n",
            EQUILIBRIUM_TOLERANCE * 100.0
        ));
        let mut sorted: Vec<_> = rows.iter().collect();
        sorted.sort_by(|x, y| {
            let dx = x.drift_relative().unwrap_or(0.0);
            let dy = y.drift_relative().unwrap_or(0.0);
            dy.partial_cmp(&dx).unwrap_or(std::cmp::Ordering::Equal)
        });
        for r in sorted.iter().take(3) {
            if let (Some(_), Some(false)) = (r.a, r.within_tolerance(EQUILIBRIUM_TOLERANCE)) {
                s.push_str(&format!(
                    "  - **{}**: {} → {} ({} relative drift)\n",
                    r.label,
                    fmt_opt(r.a),
                    fmt_opt(r.b),
                    fmt_pct(r.drift_relative()),
                ));
            }
        }
        s.push_str(&format!(
            "\nDownstream phases (Phase 8g visual revalidation) should run at \
             step ≥ 200 to land in the post-{:.0}% band, or the milestone \
             should re-evaluate the equilibrium definition.\n",
            EQUILIBRIUM_TOLERANCE * 100.0
        ));
    }
    s
}

#[test]
#[ignore]
fn v2_phase8f_equilibrium_active_medley() {
    let grid = env_grid();

    let base_spec = presets::load("active_medley").expect("active_medley preset");

    let mut spec_a = base_spec.clone();
    spec_a.grid_nx = grid;
    spec_a.grid_ny = grid;
    spec_a.steps = 100;
    spec_a.total_time_nondim = 6.0;
    spec_a.preset_label = format!("active_medley_{}sq_100steps", grid);

    let mut spec_b = base_spec.clone();
    spec_b.grid_nx = grid;
    spec_b.grid_ny = grid;
    spec_b.steps = 200;
    // Doubling total_time keeps `dt_target = total_time / steps` at the
    // same value as run A — the per-step physics is identical.
    spec_b.total_time_nondim = 12.0;
    spec_b.preset_label = format!("active_medley_{}sq_200steps", grid);

    // Deadlines scale with grid: 32² mantle-on ≈ 7s/step on the
    // development laptop, 64² ≈ 25s/step. Add 50 % margin for
    // CPU contention.
    let secs_per_step = match grid {
        32 => 7.0,
        64 => 25.0,
        _ => 25.0 * (grid as f64 / 64.0).powi(2),
    };
    let deadline_a = (100.0 * secs_per_step * 1.5).max(120.0) as u64;
    let deadline_b = (200.0 * secs_per_step * 1.5).max(120.0) as u64;

    let a = run_one(spec_a, "100steps", deadline_a);
    let b = run_one(spec_b, "200steps", deadline_b);

    let report = build_report(grid, &a, &b);

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step8_6_phase8f_equilibrium");
    std::fs::create_dir_all(&out_dir).expect("create report dir");
    let out_path = out_dir.join(format!("active_medley_{}sq.md", grid));
    std::fs::write(&out_path, &report).expect("write report");

    println!();
    println!("{}", report);
    println!("[phase8f] report written to {}", out_path.display());

    // The test never fails on the equilibrium criterion — that's a
    // diagnostic verdict, not a regression. We do require that both
    // runs completed (covered by `run_one` panic), produced finite
    // metrics, and that `cratonic_cell_fraction` (D7 static) is
    // identical between A and B (sanity probe — if it shifted,
    // something other than equilibrium is broken).
    assert!(a.metrics.vmax_peak.is_finite() && b.metrics.vmax_peak.is_finite());
    if let (Some(ca), Some(cb)) = (
        a.metrics.cratonic_cell_fraction,
        b.metrics.cratonic_cell_fraction,
    ) {
        assert!(
            (ca - cb).abs() < 1e-9,
            "cratonic_cell_fraction shifted between A and B (D7 should be static): \
             {} vs {}",
            ca,
            cb
        );
    }
}
