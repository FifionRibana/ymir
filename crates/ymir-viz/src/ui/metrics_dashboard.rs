//! Step 8.6 Phase 8c — real-time nondimensional metrics dashboard.
//!
//! Renders a left-side panel that shows:
//!
//! - During a run (`V2RunState::Running`): step / total + progress bar
//!   + ETA + wallclock; live metrics computed on the spot from the
//!   most recent `peek_state` (peak |v|, ⟨S̃⟩, max ε̇_II, cratonic
//!   fraction). These refresh at the bridge's Progress cadence.
//! - After a run (`V2RunState::Completed`): the full
//!   `metrics::Metrics` summary (vmax_peak, mass drift, CG / Newton
//!   stats, yielding fractions, mass-conservation residual, extrap
//!   stats), each with a colour band so the reviewer can scan for
//!   regressions at a glance.
//!
//! The dashboard is read-only — no user inputs flow through it. All
//! data comes from `V2SolverBridge.state`.

use bevy_egui::egui;
use bevy_egui::egui::Color32;

use crate::bridge::v2::{V2FinalState, V2RunState, V2ScalarMetrics, V2SolverBridge};

pub fn draw(ui: &mut egui::Ui, bridge: &V2SolverBridge) {
    ui.heading("Metrics dashboard");
    ui.add_space(4.0);

    match &bridge.state {
        V2RunState::Idle => {
            ui.weak("No run yet — submit a config from the right panel.");
        }
        V2RunState::Failed { error } => {
            ui.colored_label(Color32::RED, format!("Failed: {}", error));
        }
        V2RunState::Running { step, total, started_at, peek_state, spec, .. } => {
            draw_progress(ui, *step, *total, *started_at);
            ui.add_space(6.0);
            ui.separator();
            ui.add_space(4.0);
            ui.label(egui::RichText::new("Live metrics (peek_state)").strong());
            if let Some(peek) = peek_state.as_deref() {
                draw_live_metrics(ui, peek);
            } else {
                ui.weak("(awaiting first peek_state — precompute in flight)");
            }
            ui.add_space(6.0);
            ui.separator();
            ui.add_space(4.0);
            ui.label(egui::RichText::new("Run config").strong());
            ui.label(format!(
                "Grid: {}×{} · {} steps · seed {}",
                spec.grid_nx, spec.grid_ny, spec.steps, spec.seed
            ));
            ui.label(format!(
                "Plates: {} · continental {:.0}%",
                spec.num_plates,
                spec.continental_ratio * 100.0
            ));
        }
        V2RunState::Completed { metrics, elapsed, final_state, spec, .. } => {
            ui.colored_label(
                Color32::LIGHT_GREEN,
                format!("Completed in {:.1}s", elapsed.as_secs_f64()),
            );
            ui.add_space(4.0);
            ui.label(format!(
                "Grid: {}×{} · {} steps · seed {}",
                spec.grid_nx, spec.grid_ny, spec.steps, spec.seed
            ));
            ui.add_space(6.0);
            ui.separator();
            ui.add_space(4.0);
            ui.label(egui::RichText::new("Final metrics").strong());
            draw_final_metrics(ui, metrics);
            ui.add_space(6.0);
            ui.separator();
            ui.add_space(4.0);
            ui.label(egui::RichText::new("Final state").strong());
            draw_live_metrics(ui, final_state);
        }
        V2RunState::Imported {
            scalar_metrics,
            elapsed,
            final_state,
            spec,
            exported_at,
            ..
        } => {
            ui.colored_label(
                Color32::LIGHT_BLUE,
                format!("Imported snapshot · exported {}", exported_at),
            );
            ui.add_space(4.0);
            ui.label(format!(
                "Grid: {}×{} · {} steps · seed {}",
                spec.grid_nx, spec.grid_ny, spec.steps, spec.seed
            ));
            ui.label(format!("Original wallclock: {:.1}s", elapsed.as_secs_f64()));
            ui.add_space(6.0);
            ui.separator();
            ui.add_space(4.0);
            ui.label(egui::RichText::new("Snapshot metrics").strong());
            draw_scalar_metrics(ui, scalar_metrics);
            ui.add_space(6.0);
            ui.separator();
            ui.add_space(4.0);
            ui.label(egui::RichText::new("Final state").strong());
            draw_live_metrics(ui, final_state);
        }
    }
}

/// Render `V2ScalarMetrics` (post-export carrier) using the same
/// colour banding as `draw_final_metrics`. Mirrors that function's
/// layout but reads the `Option<...>` fields directly without the
/// `Metrics → newton` deref hop.
fn draw_scalar_metrics(ui: &mut egui::Ui, m: &V2ScalarMetrics) {
    metric_row(ui, "vmax peak", &format!("{:.3e}", m.vmax_peak), v_color(m.vmax_peak));

    let drift_color = if m.mass_drift_relative.abs() < 1e-6 {
        Color32::LIGHT_GREEN
    } else if m.mass_drift_relative.abs() < 1e-3 {
        Color32::YELLOW
    } else {
        Color32::LIGHT_RED
    };
    metric_row(
        ui,
        "mass drift (relative)",
        &format!("{:.3e}", m.mass_drift_relative),
        drift_color,
    );

    metric_row(
        ui,
        "CG iters mean / max",
        &format!("{:.0} / {}", m.cg_iter_mean, m.cg_iter_max),
        Color32::LIGHT_GRAY,
    );

    if let (Some(conv), Some(total)) = (m.newton_converged, m.newton_total) {
        if total > 0 {
            let pct = conv as f64 / total as f64 * 100.0;
            let color = if pct >= 99.0 {
                Color32::LIGHT_GREEN
            } else if pct >= 95.0 {
                Color32::YELLOW
            } else {
                Color32::LIGHT_RED
            };
            metric_row(
                ui,
                "Newton convergence",
                &format!("{:.1}% ({}/{})", pct, conv, total),
                color,
            );
        }
    }
    if let (Some(mean), Some(max)) = (m.newton_outer_iters_mean, m.newton_outer_iters_max) {
        metric_row(
            ui,
            "Newton outer iters mean / max",
            &format!("{:.1} / {}", mean, max),
            Color32::LIGHT_GRAY,
        );
    }
    if let Some(yf) = m.yielding_cell_fraction_max {
        metric_row(
            ui,
            "yielding cells max",
            &format!("{:.1}%", yf * 100.0),
            Color32::LIGHT_GRAY,
        );
    }
    if let Some(pyc) = m.peak_yielding_in_craton {
        let color = if pyc < 0.01 {
            Color32::LIGHT_GREEN
        } else if pyc < 0.05 {
            Color32::YELLOW
        } else {
            Color32::LIGHT_RED
        };
        metric_row(ui, "yielding in craton (peak)", &format!("{:.4}", pyc), color);
    }
    if let Some(ccf) = m.cratonic_cell_fraction {
        metric_row(
            ui,
            "cratonic cell fraction",
            &format!("{:.1}%", ccf * 100.0),
            Color32::LIGHT_GRAY,
        );
    }
    if let Some(mcr) = m.mass_conservation_residual {
        let color = if mcr < 1e-6 {
            Color32::LIGHT_GREEN
        } else if mcr < 1e-3 {
            Color32::YELLOW
        } else {
            Color32::LIGHT_RED
        };
        metric_row(
            ui,
            "mass conservation residual",
            &format!("{:.3e}", mcr),
            color,
        );
    }
    if let (Some(att), Some(app)) = (m.extrap_attempted, m.extrap_applied) {
        if att > 0 {
            let pct = app as f64 / att as f64 * 100.0;
            metric_row(
                ui,
                "extrapolation applied / attempted",
                &format!("{} / {} ({:.0}%)", app, att, pct),
                Color32::LIGHT_GRAY,
            );
        }
    }
    if let (Some(fbc), Some(att)) = (m.extrap_fallback_count, m.extrap_attempted) {
        if att > 0 {
            let pct = fbc as f64 / att as f64 * 100.0;
            let color = if pct < 5.0 {
                Color32::LIGHT_GREEN
            } else if pct < 20.0 {
                Color32::YELLOW
            } else {
                Color32::LIGHT_RED
            };
            metric_row(
                ui,
                "extrapolation fallback rate",
                &format!("{:.1}%", pct),
                color,
            );
        }
    }
    metric_row(
        ui,
        "wallclock per step (mean)",
        &format!("{:.1} ms", m.wallclock_per_step_mean_s * 1000.0),
        Color32::LIGHT_GRAY,
    );
}

fn draw_progress(ui: &mut egui::Ui, step: usize, total: usize, started_at: Option<std::time::Instant>) {
    let frac = if total > 0 {
        (step as f32 / total as f32).clamp(0.0, 1.0)
    } else {
        0.0
    };
    ui.label(format!("Step {} / {}", step, total));
    ui.add(egui::ProgressBar::new(frac).text(format!("{:.0}%", frac * 100.0)));

    if let Some(start) = started_at {
        let elapsed = start.elapsed().as_secs_f64();
        let eta = if step > 0 && total > step {
            let per_step = elapsed / step as f64;
            (total - step) as f64 * per_step
        } else {
            0.0
        };
        ui.label(format!("Elapsed {:.1}s · ETA {:.1}s", elapsed, eta));
    }
}

/// Live metrics computed from a `V2FinalState` (peek or final). Intentionally
/// scoped to what's directly observable from the raster snapshot — no
/// solver-internal counters, those come from `Metrics`.
fn draw_live_metrics(ui: &mut egui::Ui, state: &V2FinalState) {
    // peak |v| — colour band per Phase 8c spec: green < 1, yellow 1–3,
    // red > 3 (nondim). Active medley typically lands in the yellow
    // band; runaway shows up red on the dashboard before the run
    // completes.
    let peak_v = state
        .vx
        .iter()
        .zip(state.vy.iter())
        .map(|(&vx, &vy)| (vx * vx + vy * vy).sqrt())
        .fold(0.0_f64, f64::max);
    metric_row(ui, "peak |v|", &format!("{:.3e}", peak_v), v_color(peak_v));

    // ⟨S̃⟩ — mean thickness; mass conservation under the closed-mode
    // recycling sees this stay close to the initial mean across the
    // run. Use orange when |Δ⟨S̃⟩ - 0.6| / 0.6 > 5% as a sanity probe.
    let mean_s: f64 = state.s_field.iter().sum::<f64>() / state.s_field.len() as f64;
    metric_row(ui, "⟨S̃⟩", &format!("{:.4}", mean_s), Color32::LIGHT_GRAY);

    // max ε̇_II — second invariant of the strain-rate tensor. The Phase
    // 7 follow-up widened the rendering log range to [1e-3, 1e2]; values
    // above 1e2 indicate localised yielding hotspots.
    let max_eps_ii = state
        .strain_rate_invariant
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    metric_row(ui, "max ε̇_II", &format!("{:.3e}", max_eps_ii), Color32::LIGHT_GRAY);

    // Cratonic cell fraction (cells where the cratonic factor f ≥ 0.5).
    // Only meaningful with cratonic immunity enabled; absent otherwise.
    if let Some(cratonic) = &state.cratonic_factor {
        let n = cratonic.len();
        let count = cratonic.iter().filter(|&&v| v >= 0.5).count();
        let frac = count as f64 / n.max(1) as f64;
        metric_row(
            ui,
            "cratonic cells (f ≥ 0.5)",
            &format!("{:.1}%", frac * 100.0),
            Color32::LIGHT_GRAY,
        );
    }
}

fn draw_final_metrics(ui: &mut egui::Ui, metrics: &ymir_core::tectonics_v2::diagnostics::metrics::Metrics) {
    metric_row(
        ui,
        "vmax peak",
        &format!("{:.3e}", metrics.vmax_peak),
        v_color(metrics.vmax_peak),
    );

    let drift = metrics.mass_drift_relative;
    let drift_color = if drift.abs() < 1e-6 {
        Color32::LIGHT_GREEN
    } else if drift.abs() < 1e-3 {
        Color32::YELLOW
    } else {
        Color32::LIGHT_RED
    };
    metric_row(
        ui,
        "mass drift (relative)",
        &format!("{:.3e}", drift),
        drift_color,
    );

    metric_row(
        ui,
        "CG iters mean / max",
        &format!("{:.0} / {}", metrics.cg_iter_mean, metrics.cg_iter_max),
        Color32::LIGHT_GRAY,
    );

    if let Some(newton) = &metrics.newton {
        let total = newton.converged + newton.stalled + newton.diverged + newton.capped;
        if total > 0 {
            let conv_pct = newton.converged as f64 / total as f64 * 100.0;
            let conv_color = if conv_pct >= 99.0 {
                Color32::LIGHT_GREEN
            } else if conv_pct >= 95.0 {
                Color32::YELLOW
            } else {
                Color32::LIGHT_RED
            };
            metric_row(
                ui,
                "Newton convergence",
                &format!("{:.1}% ({}/{})", conv_pct, newton.converged, total),
                conv_color,
            );
        }
        if !newton.outer_iters.is_empty() {
            let mean_outer: f64 =
                newton.outer_iters.iter().map(|&n| n as f64).sum::<f64>()
                    / newton.outer_iters.len() as f64;
            let max_outer = newton.outer_iters.iter().copied().max().unwrap_or(0);
            metric_row(
                ui,
                "Newton outer iters mean / max",
                &format!("{:.1} / {}", mean_outer, max_outer),
                Color32::LIGHT_GRAY,
            );
        }
        if let Some(yf) = newton.yielding_cell_fraction_max {
            metric_row(
                ui,
                "yielding cells max",
                &format!("{:.1}%", yf * 100.0),
                Color32::LIGHT_GRAY,
            );
        }
        if let Some(pyc) = newton.peak_yielding_in_craton {
            // §9 acceptance: peak_yielding_in_craton should stay below
            // ~0.01 on cratonic-on runs. Above that = immunity broken.
            let craton_color = if pyc < 0.01 {
                Color32::LIGHT_GREEN
            } else if pyc < 0.05 {
                Color32::YELLOW
            } else {
                Color32::LIGHT_RED
            };
            metric_row(
                ui,
                "yielding in craton (peak)",
                &format!("{:.4}", pyc),
                craton_color,
            );
        }
        if let Some(ccf) = newton.cratonic_cell_fraction {
            metric_row(
                ui,
                "cratonic cell fraction",
                &format!("{:.1}%", ccf * 100.0),
                Color32::LIGHT_GRAY,
            );
        }
        if let Some(mcr) = newton.mass_conservation_residual {
            // Step 6+ closed-mode acceptance: residual < 1e-6.
            let mcr_color = if mcr < 1e-6 {
                Color32::LIGHT_GREEN
            } else if mcr < 1e-3 {
                Color32::YELLOW
            } else {
                Color32::LIGHT_RED
            };
            metric_row(
                ui,
                "mass conservation residual",
                &format!("{:.3e}", mcr),
                mcr_color,
            );
        }
    }

    if let Some(extrap) = &metrics.extrapolation {
        if extrap.attempted > 0 {
            let applied_pct = extrap.applied as f64 / extrap.attempted as f64 * 100.0;
            let fallback_pct = extrap.fallback_indices.len() as f64
                / extrap.attempted as f64
                * 100.0;
            metric_row(
                ui,
                "extrapolation applied / attempted",
                &format!(
                    "{} / {} ({:.0}%)",
                    extrap.applied, extrap.attempted, applied_pct
                ),
                Color32::LIGHT_GRAY,
            );
            if !extrap.fallback_indices.is_empty() {
                let fallback_color = if fallback_pct < 5.0 {
                    Color32::LIGHT_GREEN
                } else if fallback_pct < 20.0 {
                    Color32::YELLOW
                } else {
                    Color32::LIGHT_RED
                };
                metric_row(
                    ui,
                    "extrapolation fallback rate",
                    &format!("{:.1}%", fallback_pct),
                    fallback_color,
                );
            }
        }
    }

    metric_row(
        ui,
        "wallclock per step (mean)",
        &format!("{:.1} ms", metrics.wallclock_per_step_mean.as_secs_f64() * 1000.0),
        Color32::LIGHT_GRAY,
    );
}

fn v_color(peak_v: f64) -> Color32 {
    if peak_v < 1.0 {
        Color32::LIGHT_GREEN
    } else if peak_v <= 3.0 {
        Color32::YELLOW
    } else {
        Color32::LIGHT_RED
    }
}

fn metric_row(ui: &mut egui::Ui, label: &str, value: &str, value_color: Color32) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.colored_label(value_color, egui::RichText::new(value).monospace());
        });
    });
}
