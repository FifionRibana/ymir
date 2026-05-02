//! Step 8.6 Phase 8e — export / import of a completed v2 run.
//!
//! Serialises a `(spec, final_state, scalar_metrics)` triple to a
//! single JSON file that round-trips back into a viewable run state
//! without re-running the solver. The format is versioned so future
//! schema changes don't silently corrupt older artefacts:
//!
//! ```json
//! {
//!   "format_version": 1,
//!   "exported_at": "2026-04-30T12:34:56Z",
//!   "elapsed_seconds": 42.3,
//!   "spec": { … V2RunSpec … },
//!   "scalar_metrics": { … V2ScalarMetrics … },
//!   "final_state": { … V2FinalState … }
//! }
//! ```
//!
//! Why a custom `V2ScalarMetrics` instead of the full
//! `ymir_core::…::Metrics`?
//! The harness `Metrics` struct does not carry serde derives and
//! contains heavy histograms / per-step series the import path doesn't
//! need. We pull just the dashboard-relevant scalar fields here.
//!
//! Raster fields (`s_field`, `vx`, `vy`, …) ride along inside the JSON
//! as f64 / u16 / u8 arrays. At 64² × ~9 fields × 8 bytes the file is
//! ~36 KB; at 128² it is ~144 KB — still small enough that base-text
//! JSON is the simplest encoding. If we later move to 256²+ regimes
//! the format can switch to a `data: "snapshot.bin"` sidecar without
//! breaking the schema (bump `format_version`).

use std::path::Path;

use serde::{Deserialize, Serialize};
use ymir_core::tectonics_v2::diagnostics::metrics::Metrics;

use super::events::V2FinalState;
use super::spec::V2RunSpec;

pub const SNAPSHOT_FORMAT_VERSION: u32 = 1;

/// Scalar metrics extracted from `Metrics` for dashboard / export.
/// Each field is `Option<f64>` so disabled mechanisms (cratonic off,
/// boundary off, …) round-trip as `None` rather than 0.0 (which would
/// be a meaningful value for an active mechanism).
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct V2ScalarMetrics {
    pub vmax_peak: f64,
    pub mass_drift_relative: f64,
    pub cg_iter_mean: f64,
    pub cg_iter_max: usize,
    pub wallclock_per_step_mean_s: f64,

    /// Newton outcome counts (None when no Newton was run, e.g. Step 0).
    pub newton_converged: Option<usize>,
    pub newton_total: Option<usize>,
    pub newton_outer_iters_mean: Option<f64>,
    pub newton_outer_iters_max: Option<u32>,

    pub yielding_cell_fraction_max: Option<f64>,
    pub peak_yielding_in_craton: Option<f64>,
    pub cratonic_cell_fraction: Option<f64>,
    pub mass_conservation_residual: Option<f64>,

    pub extrap_attempted: Option<usize>,
    pub extrap_applied: Option<usize>,
    pub extrap_fallback_count: Option<usize>,
}

impl V2ScalarMetrics {
    pub fn from_metrics(m: &Metrics) -> Self {
        let (
            newton_converged,
            newton_total,
            newton_outer_iters_mean,
            newton_outer_iters_max,
            yielding_cf_max,
            peak_yielding_in_craton,
            cratonic_cell_fraction,
            mass_conservation_residual,
        ) = if let Some(n) = &m.newton {
            let total = n.converged + n.stalled + n.diverged + n.capped;
            let mean = if n.outer_iters.is_empty() {
                None
            } else {
                let s: u32 = n.outer_iters.iter().sum();
                Some(s as f64 / n.outer_iters.len() as f64)
            };
            let max = n.outer_iters.iter().copied().max();
            (
                Some(n.converged),
                Some(total),
                mean,
                max,
                n.yielding_cell_fraction_max,
                n.peak_yielding_in_craton,
                n.cratonic_cell_fraction,
                n.mass_conservation_residual,
            )
        } else {
            (None, None, None, None, None, None, None, None)
        };

        let (extrap_attempted, extrap_applied, extrap_fallback_count) =
            if let Some(e) = &m.extrapolation {
                (
                    Some(e.attempted),
                    Some(e.applied),
                    Some(e.fallback_indices.len()),
                )
            } else {
                (None, None, None)
            };

        V2ScalarMetrics {
            vmax_peak: m.vmax_peak,
            mass_drift_relative: m.mass_drift_relative,
            cg_iter_mean: m.cg_iter_mean,
            cg_iter_max: m.cg_iter_max,
            wallclock_per_step_mean_s: m.wallclock_per_step_mean.as_secs_f64(),
            newton_converged,
            newton_total,
            newton_outer_iters_mean,
            newton_outer_iters_max,
            yielding_cell_fraction_max: yielding_cf_max,
            peak_yielding_in_craton,
            cratonic_cell_fraction,
            mass_conservation_residual,
            extrap_attempted,
            extrap_applied,
            extrap_fallback_count,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct V2RunSnapshot {
    pub format_version: u32,
    pub exported_at: String,
    pub elapsed_seconds: f64,
    pub spec: V2RunSpec,
    pub scalar_metrics: V2ScalarMetrics,
    pub final_state: V2FinalState,
}

impl V2RunSnapshot {
    pub fn new(
        spec: V2RunSpec,
        final_state: V2FinalState,
        metrics: &Metrics,
        elapsed: std::time::Duration,
    ) -> Self {
        Self {
            format_version: SNAPSHOT_FORMAT_VERSION,
            exported_at: rfc3339_now(),
            elapsed_seconds: elapsed.as_secs_f64(),
            spec,
            scalar_metrics: V2ScalarMetrics::from_metrics(metrics),
            final_state,
        }
    }

    /// Save as pretty-printed JSON. Creates the parent directory on
    /// demand so a fresh `output_dir/snapshots/` tree works.
    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let json = serde_json::to_string_pretty(self).map_err(io_other)?;
        std::fs::write(path, json)
    }

    pub fn load(path: &Path) -> std::io::Result<Self> {
        let bytes = std::fs::read(path)?;
        let snap: Self = serde_json::from_slice(&bytes).map_err(io_other)?;
        if snap.format_version != SNAPSHOT_FORMAT_VERSION {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "snapshot format_version {} not supported (expected {})",
                    snap.format_version, SNAPSHOT_FORMAT_VERSION
                ),
            ));
        }
        Ok(snap)
    }
}

fn io_other<E: std::fmt::Display>(e: E) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string())
}

/// Produce an RFC3339 UTC timestamp without pulling chrono. Format:
/// `YYYY-MM-DDTHH:MM:SSZ`. Falls back to `unix-epoch-secs N` if the
/// system clock is somehow before 1970.
fn rfc3339_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let dur = match SystemTime::now().duration_since(UNIX_EPOCH) {
        Ok(d) => d,
        Err(_) => return "unix-epoch-secs 0".to_string(),
    };
    let secs = dur.as_secs();
    let (y, m, d, hh, mm, ss) = unix_to_ymdhms(secs);
    format!("{y:04}-{m:02}-{d:02}T{hh:02}:{mm:02}:{ss:02}Z")
}

fn unix_to_ymdhms(t: u64) -> (i64, u32, u32, u32, u32, u32) {
    // Civil-from-days algorithm (Hinnant 2014). Self-contained so we
    // don't add a date dep just for an export timestamp.
    let secs_per_day: u64 = 86_400;
    let days = (t / secs_per_day) as i64;
    let rem = (t % secs_per_day) as u32;
    let hh = rem / 3600;
    let mm = (rem % 3600) / 60;
    let ss = rem % 60;
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d, hh, mm, ss)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::v2::events::V2FinalState;
    use crate::bridge::v2::spec::V2RunSpec;

    fn fake_state(nx: usize, ny: usize) -> V2FinalState {
        V2FinalState {
            nx,
            ny,
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            s_field: vec![0.5; nx * ny],
            vx: vec![0.01; nx * ny],
            vy: vec![-0.02; nx * ny],
            strain_rate_invariant: vec![1.0; nx * ny],
            age_field: Some(vec![3.0; nx * ny]),
            cratonic_factor: Some(vec![0.0; nx * ny]),
            plate_id: Some(vec![0; nx * ny]),
            plate_type: Some(vec![1; nx * ny]),
            boundary_flag: Some(vec![0; nx * ny]),
        }
    }

    #[test]
    fn snapshot_roundtrips_through_disk() {
        let spec = V2RunSpec::active_medley_defaults();
        let state = fake_state(8, 8);
        let scalar = V2ScalarMetrics {
            vmax_peak: 1.234,
            mass_drift_relative: 1.0e-7,
            cg_iter_mean: 12.5,
            cg_iter_max: 30,
            wallclock_per_step_mean_s: 0.5,
            newton_converged: Some(99),
            newton_total: Some(100),
            newton_outer_iters_mean: Some(2.5),
            newton_outer_iters_max: Some(8),
            yielding_cell_fraction_max: Some(0.42),
            peak_yielding_in_craton: Some(0.001),
            cratonic_cell_fraction: Some(0.18),
            mass_conservation_residual: Some(1.0e-9),
            extrap_attempted: Some(50),
            extrap_applied: Some(45),
            extrap_fallback_count: Some(5),
        };
        let snap = V2RunSnapshot {
            format_version: SNAPSHOT_FORMAT_VERSION,
            exported_at: "2026-04-30T12:00:00Z".to_string(),
            elapsed_seconds: 12.34,
            spec,
            scalar_metrics: scalar,
            final_state: state,
        };

        let dir = std::env::temp_dir().join("ymir_snapshot_roundtrip");
        let path = dir.join("snap.json");
        snap.save(&path).expect("save");
        let back = V2RunSnapshot::load(&path).expect("load");

        assert_eq!(back.format_version, snap.format_version);
        assert_eq!(back.exported_at, snap.exported_at);
        assert!((back.elapsed_seconds - snap.elapsed_seconds).abs() < 1e-12);
        assert_eq!(back.scalar_metrics, snap.scalar_metrics);
        assert_eq!(back.final_state.nx, 8);
        assert_eq!(back.final_state.ny, 8);
        assert_eq!(back.final_state.s_field, snap.final_state.s_field);
        assert_eq!(back.final_state.vx, snap.final_state.vx);
        assert_eq!(back.final_state.vy, snap.final_state.vy);

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn load_rejects_unsupported_format_version() {
        // Build a valid snapshot first so the rest of the schema
        // deserialises cleanly; then bump `format_version` to a
        // future value and assert `load` rejects it. Going through a
        // real fixture (rather than a hand-rolled minimal JSON) means
        // the test doesn't accidentally fail on serde's
        // missing-field path before it ever reaches the
        // format_version check.
        let mut snap = V2RunSnapshot {
            format_version: SNAPSHOT_FORMAT_VERSION,
            exported_at: "2026-04-30T12:00:00Z".to_string(),
            elapsed_seconds: 0.0,
            spec: V2RunSpec::active_medley_defaults(),
            scalar_metrics: V2ScalarMetrics::default(),
            final_state: fake_state(4, 4),
        };
        snap.format_version = SNAPSHOT_FORMAT_VERSION + 1;

        let dir = std::env::temp_dir().join("ymir_snapshot_bad_version");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bad.json");
        let json = serde_json::to_string(&snap).expect("serialize");
        std::fs::write(&path, json).unwrap();

        let err = V2RunSnapshot::load(&path).expect_err("must reject");
        assert!(
            err.to_string().contains("format_version"),
            "expected format_version error, got: {}",
            err
        );
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn rfc3339_format_matches_iso_pattern() {
        let s = rfc3339_now();
        assert_eq!(s.len(), 20, "got {}", s);
        assert!(s.ends_with('Z'));
        assert_eq!(&s[4..5], "-");
        assert_eq!(&s[7..8], "-");
        assert_eq!(&s[10..11], "T");
    }
}
