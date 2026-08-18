//! C1 product path WITH the generation cache wired in.
//!
//! Thin wrappers around the (unchanged, pure) generation functions
//! `run_with_closures` + `upscale_from_c1`, adding content-addressed caching of
//! the expensive eroded HD heightmap. The generation functions themselves are
//! NOT touched — the no-cache paths stay byte-identical. See
//! `docs/design/c1_generation_cache.md`.
//!
//! **Key completeness is the load-bearing property.** The cache key must encode
//! EVERYTHING the eroded output depends on, or a changed-but-unencoded input
//! serves a stale terrain (the exact bug content-addressing is meant to make
//! impossible). The key here is built so every input is captured:
//!   - `seed`, `grid` — scalars;
//!   - `init` (`Phase2InitParams`), `run` (`C1TimeLoopConfig`, which embeds the
//!     `IsostasyConfig`), `closures` (`C1Closures`, all 7 nested closure param
//!     structs) — folded via their DERIVED `Debug` (walks every field);
//!   - `ss` (`SteinSteinParams`), `upscale` (`FbmUpscaleConfig`, which embeds
//!     the `Option<ErosionConfig>`) — folded via `Serialize`;
//!   - `ALGO_*` — the source-code version (the one input not in any config).
//!
//! The eroded key is chained onto an upstream "tectonic" key (`derived_from`),
//! so any change to a tectonic input flips the eroded digest too.

use std::path::{Path, PathBuf};

use serde_json::Value;

use crate::cache::{
    ALGO_DRAINAGE, ALGO_TECTONICS, ALGO_UPSCALE_EROSION, CacheKey, RawCodec, cached,
    cached_fallible,
};
use crate::climate::c1_climate_windowed;
use crate::climate::precipitation::PrecipParams;
use crate::export::raw;
use crate::grid::GridF32;
use crate::seed::WorldSeed;
use crate::terrain::flow::FlowResult;
use crate::terrain::upscale::FbmUpscaleConfig;

use super::closures::oceanic_bathymetry::SteinSteinParams;
use super::drainage::C1_SEA_LEVEL_NORM;
use super::drainage::{C1DrainageConfig, C1DrainageResult, DrainageClimate, c1_drainage_windowed};
use super::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use super::kinematics::PlateKinematics;
use super::land_topology::{IslandCriteria, LandTopology, is_island_fit, land_topology};
use super::production_upscale::{
    C1_DOMAIN_KM, EroProgress, c1_coarse_normalized_altitude, c1_coarse_raw_altitude,
    c1_land_centroid_normalized, c1_normalize_coarse, upscale_from_c1_with_progress,
};
use super::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};

/// The cache key for the 64² tectonic build (`init` + `run_with_closures`).
/// Captures every input to that step. Exposed so the eroded/drainage keys can
/// chain onto it (`CacheKey::derived_from`).
pub fn tectonic_key(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
) -> CacheKey {
    CacheKey::root()
        .with("seed", &seed)
        .with("grid", &grid)
        .with_debug("init", init)
        .with_debug("run", run)
        .with_debug("closures", closures)
        .algo(ALGO_TECTONICS)
}

/// The cache key for the eroded HD heightmap — chained onto [`tectonic_key`]
/// (upstream change → this digest changes) plus the upscale-specific inputs.
pub fn eroded_key(
    tectonic: &CacheKey,
    ss: &SteinSteinParams,
    upscale_cfg: &FbmUpscaleConfig,
) -> CacheKey {
    CacheKey::derived_from(tectonic)
        .with("ss", ss)
        .with("upscale", upscale_cfg)
        .algo(ALGO_UPSCALE_EROSION)
}

/// Build the C1 eroded HD heightmap, reusing the cached `.raw` if the inputs are
/// unchanged. On a MISS this runs the full (expensive) build: `init` →
/// `run_with_closures` → `upscale_from_c1`; on a HIT it loads the heightmap and
/// skips all of it (the ~minutes of erosion). Byte-identical to the no-cache
/// path either way (`upscale_from_c1` is untouched).
///
/// `run.iso_config` is the isostasy used both by the time loop and by the
/// upscale (matching the harness convention), so it is not passed separately.
#[allow(clippy::too_many_arguments)]
pub fn cached_c1_eroded(
    cache_dir: &Path,
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    upscale_cfg: &FbmUpscaleConfig,
) -> Result<GridF32, String> {
    cached_c1_eroded_with_progress(
        cache_dir,
        seed,
        grid,
        init,
        run,
        closures,
        ss,
        upscale_cfg,
        &mut |_| {},
        &|| false,
    )
}

/// [`cached_c1_eroded`] with sub-phase progress + mid-build cancel (UI frieze
/// split, suite e). On a MISS the compute emits `EroProgress` (Tectonic step,
/// Relief, Erosion %, Bathymetry) via `progress` and polls `cancel()`; a
/// cancelled build returns `Err` and is NOT written to the cache (via
/// `cached_fallible`). On a HIT nothing runs (no progress) — the caller detects
/// the HIT itself (sidecar) and shows the cached state. Byte-identical payload.
#[allow(clippy::too_many_arguments)]
pub fn cached_c1_eroded_with_progress(
    cache_dir: &Path,
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    upscale_cfg: &FbmUpscaleConfig,
    progress: &mut dyn FnMut(EroProgress),
    cancel: &dyn Fn() -> bool,
) -> Result<GridF32, String> {
    let key = eroded_key(&tectonic_key(seed, grid, init, run, closures), ss, upscale_cfg);
    cached_fallible(cache_dir, "eroded", &key, || {
        let mut state = init_c1_state_phase_2_r7(grid, seed, init);
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let total = run.n_steps;
        run_with_closures(&mut state, &mut kin, run, closures, |step, _| {
            progress(EroProgress::Tectonic { step: step + 1, total });
        });
        if cancel() {
            return Err("cancelled".to_string());
        }
        let up = upscale_from_c1_with_progress(
            &state,
            &run.iso_config,
            ss,
            &WorldSeed::new(seed),
            upscale_cfg,
            progress,
            cancel,
        );
        if cancel() {
            return Err("cancelled".to_string());
        }
        Ok(up.heightmap)
    })
}

/// Land centroid (normalized `[u, v]`) of the coarse C1 continent for these
/// tectonic inputs — used to CENTRE a cropped export window on the continent.
/// Runs the coarse tectonic sim (cheap vs the HD erosion) and reads the SAME
/// normalized altitude the upscale samples. Deterministic in the tectonic key.
///
/// Perf follow-up: NOT cached — on an eroded-cache HIT this still re-runs the
/// coarse tectonics. Could share a cached coarse-altitude pass with
/// `cached_c1_eroded` (the tectonic loop is the cheap phase, so deferred).
pub fn c1_land_centroid(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    target_land_fraction: Option<f32>,
) -> [f64; 2] {
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, run, closures, |_, _| {});
    let coarse = c1_coarse_normalized_altitude(&state, &run.iso_config, ss, target_land_fraction);
    c1_land_centroid_normalized(&coarse)
}

/// Coarse-field land report (M1): the window-origin centroid AND the land-topology
/// diagnostics, from ONE coarse pass (avoids a second tectonic run). Computed on
/// the calibrated coarse field (the same one the upscale renders), so the metrics
/// describe the continent that will actually be exported.
pub struct CoarseLandReport {
    /// Normalized land centroid `[u, v]` (window-origin anchor).
    pub centroid: [f64; 2],
    /// Full-torus land topology (number of masses, largest, wrap flags, bbox).
    pub topology: LandTopology,
}

/// Run the coarse C1 sim and report both the land centroid and the land topology.
pub fn c1_coarse_land_report(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    target_land_fraction: Option<f32>,
) -> CoarseLandReport {
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, run, closures, |_, _| {});
    let coarse = c1_coarse_normalized_altitude(&state, &run.iso_config, ss, target_land_fraction);
    CoarseLandReport {
        centroid: c1_land_centroid_normalized(&coarse),
        topology: land_topology(&coarse, C1_SEA_LEVEL_NORM),
    }
}

/// Sweep the land topology of ONE tectonic config across many target land
/// fractions, running the (expensive) coarse tectonic pass ONCE and evaluating
/// each `tlf` by re-thresholding the same raw altitude field (cheap). Returns
/// `(tlf, topology)` per fraction. Window/margin budgets are a POST-HOC predicate
/// on these metrics — not swept here. Deterministic in the tectonic key.
#[allow(clippy::too_many_arguments)]
pub fn coarse_land_topology_sweep(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    target_land_fractions: &[f32],
) -> Vec<(f32, LandTopology)> {
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, run, closures, |_, _| {});
    let raw = c1_coarse_raw_altitude(&state, &run.iso_config, ss);
    target_land_fractions
        .iter()
        .map(|&f| {
            let norm = c1_normalize_coarse(raw.clone(), Some(f));
            (f, land_topology(&norm, C1_SEA_LEVEL_NORM))
        })
        .collect()
}

/// A tectonic configuration whose largest landmass passes [`is_island_fit`].
#[derive(Debug, Clone, Copy)]
pub struct IslandHit {
    pub seed: u64,
    pub num_plates: usize,
    pub seed_cluster_count: usize,
    pub topology: LandTopology,
    pub centroid: [f64; 2],
}

/// Seed/plate search (M1 #190) — the tool that closes the geometric budget.
/// Scans `plate_counts × cluster_counts × seeds` and returns the FIRST config
/// whose largest landmass does not wrap the torus and fits the export window
/// with an ocean margin (via [`is_island_fit`]). Reuses [`c1_coarse_land_report`]
/// (coarse only — cheap vs HD erosion), so the search runs one tectonic pass per
/// candidate. Deterministic: the scan order is exactly the input order.
#[allow(clippy::too_many_arguments)]
pub fn find_island_config(
    grid: usize,
    base_init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    target_land_fraction: Option<f32>,
    plate_counts: &[usize],
    cluster_counts: &[usize],
    seeds: &[u64],
    criteria: &IslandCriteria,
) -> Option<IslandHit> {
    for &num_plates in plate_counts {
        for &seed_cluster_count in cluster_counts {
            let mut init = base_init.clone();
            init.num_plates = num_plates;
            init.cluster.seed_cluster_count = seed_cluster_count;
            for &seed in seeds {
                let report = c1_coarse_land_report(
                    seed,
                    grid,
                    &init,
                    run,
                    closures,
                    ss,
                    target_land_fraction,
                );
                if is_island_fit(&report.topology, criteria) {
                    return Some(IslandHit {
                        seed,
                        num_plates,
                        seed_cluster_count,
                        topology: report.topology,
                        centroid: report.centroid,
                    });
                }
            }
        }
    }
    None
}

// ── Drainage (composite) cache ─────────────────────────────────────────────

/// Append `_{name}.raw` to a stem path (`dir/drainage_<digest>` → `…_filled.raw`).
fn part(stem: &Path, name: &str) -> PathBuf {
    let mut s = stem.as_os_str().to_owned();
    s.push(format!("_{name}.raw"));
    PathBuf::from(s)
}

/// `C1DrainageResult` is COMPOSITE: the large rasters (filled, accumulation,
/// direction, basins, lake_map) go to one `.raw` each via the shared codec; the
/// structured network (rivers, per-segment km²/navigability, typed lakes) and
/// the scalars (dims, num_basins) ride in the sidecar `shape()` JSON. Every
/// field is persisted AND restored — a forgotten field would be a subtly
/// partial structure (the completeness trap for composites).
impl RawCodec for C1DrainageResult {
    fn shape(&self) -> Value {
        serde_json::json!({
            "width": self.width,
            "height": self.height,
            "num_basins": self.flow.num_basins,
            "rivers": self.rivers,
            "seg_km2": self.segment_drainage_km2,
            "seg_nav": self.segment_navigability,
            "lakes": self.lakes,
        })
    }

    fn write_raw(&self, stem: &Path) -> Result<(), String> {
        self.flow.filled.save_raw(&part(stem, "filled"))?;
        self.flow.accumulation.save_raw(&part(stem, "accum"))?;
        raw::save_u8(&part(stem, "dir"), &self.flow.direction)?;
        raw::save_u32(&part(stem, "basins"), &self.flow.basins)?;
        raw::save_u32(&part(stem, "lakemap"), &self.lake_map)?;
        Ok(())
    }

    fn read_raw(stem: &Path, shape: &Value) -> Result<Self, String> {
        let w = shape["width"].as_u64().ok_or("drainage sidecar: width")? as usize;
        let h = shape["height"].as_u64().ok_or("drainage sidecar: height")? as usize;
        let n = w * h;
        let num_basins = shape["num_basins"].as_u64().ok_or("drainage sidecar: num_basins")? as u32;

        let filled = GridF32::load_raw(&part(stem, "filled"), w, h)?;
        let accumulation = GridF32::load_raw(&part(stem, "accum"), w, h)?;
        let direction = raw::load_u8(&part(stem, "dir"))?;
        if direction.len() != n {
            return Err(format!("drainage direction: expected {n}, got {}", direction.len()));
        }
        let basins = raw::load_u32(&part(stem, "basins"), n)?;
        let lake_map = raw::load_u32(&part(stem, "lakemap"), n)?;

        // NB: separate `from_value` calls (not a shared closure) so each infers
        // its own target type from the struct field it feeds.
        let rivers = serde_json::from_value(shape["rivers"].clone())
            .map_err(|e| format!("drainage sidecar rivers: {e}"))?;
        let segment_drainage_km2 = serde_json::from_value(shape["seg_km2"].clone())
            .map_err(|e| format!("drainage sidecar seg_km2: {e}"))?;
        let segment_navigability = serde_json::from_value(shape["seg_nav"].clone())
            .map_err(|e| format!("drainage sidecar seg_nav: {e}"))?;
        let lakes = serde_json::from_value(shape["lakes"].clone())
            .map_err(|e| format!("drainage sidecar lakes: {e}"))?;

        Ok(C1DrainageResult {
            flow: FlowResult { filled, direction, accumulation, basins, num_basins },
            rivers,
            segment_drainage_km2,
            segment_navigability,
            lakes,
            lake_map,
            width: w,
            height: h,
        })
    }
}

/// The cache key for the drainage product — chained onto the eroded key
/// (upstream change → this digest changes) plus the drainage-specific inputs.
/// `climate` (#drainage water balance) is `Some((latitude_deg, precip_params))`
/// when the hydroclimate layer is on: it is folded into the key so a CLIMATE
/// change re-invalidates the drainage (the new dependency). `None` → the
/// pure-geometry path, byte-identical key to pre-#drainage.
pub fn drainage_key(
    eroded: &CacheKey,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    climate: Option<(f32, &PrecipParams)>,
) -> CacheKey {
    drainage_key_windowed(eroded, cfg, ss, climate, C1_DOMAIN_KM)
}

/// [`drainage_key`] with the window horizontal scale folded in. The window is
/// only added to the digest when it is an ACTUAL crop (`window_km != full
/// torus`), so the full-domain key stays byte-identical to pre-window.
pub fn drainage_key_windowed(
    eroded: &CacheKey,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    climate: Option<(f32, &PrecipParams)>,
    window_km: f32,
) -> CacheKey {
    let mut k = CacheKey::derived_from(eroded).with("drainage_cfg", cfg).with("ss", ss);
    if let Some((lat, pp)) = climate {
        k = k.with("drainage_lat", &lat).with("drainage_precip", pp);
    }
    if window_km != C1_DOMAIN_KM {
        k = k.with("drainage_window_km", &window_km);
    }
    k.algo(ALGO_DRAINAGE)
}

/// Extract the C1 drainage product, reusing the cached payload if `eroded`
/// (the upstream heightmap key) and the drainage inputs are unchanged.
///
/// `climate` = `Some((latitude_deg, &PrecipParams))` activates the #drainage
/// WATER BALANCE: the climate is computed INTERNALLY from the heightmap (so the
/// wrapper owns both the grids it feeds `c1_drainage` AND the key inputs), and
/// folded into the key. `None` → pure-geometry fill-and-spill (byte-identical to
/// pre-#drainage). `c1_drainage` is untouched (pure); this only wraps it.
pub fn cached_c1_drainage(
    cache_dir: &Path,
    eroded: &CacheKey,
    heightmap: &GridF32,
    climate: Option<(f32, &PrecipParams)>,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
) -> Result<C1DrainageResult, String> {
    cached_c1_drainage_windowed(cache_dir, eroded, heightmap, climate, cfg, ss, C1_DOMAIN_KM)
}

/// [`cached_c1_drainage`] for a grid spanning `window_km` (a cropped window):
/// the internal climate + drainage use the window horizontal scale, and the key
/// folds it in. `window_km == C1_DOMAIN_KM` reproduces [`cached_c1_drainage`]
/// exactly (byte-identical key + payload).
#[allow(clippy::too_many_arguments)]
pub fn cached_c1_drainage_windowed(
    cache_dir: &Path,
    eroded: &CacheKey,
    heightmap: &GridF32,
    climate: Option<(f32, &PrecipParams)>,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    window_km: f32,
) -> Result<C1DrainageResult, String> {
    let key = drainage_key_windowed(eroded, cfg, ss, climate, window_km);
    cached(cache_dir, "drainage", &key, || match climate {
        Some((lat, pp)) => {
            let clim = c1_climate_windowed(heightmap, ss, lat, pp, window_km);
            let dc = DrainageClimate {
                precip_internal: &clim.precipitation,
                temperature: &clim.temperature,
            };
            c1_drainage_windowed(heightmap, Some(&dc), cfg, ss, window_km)
        }
        None => c1_drainage_windowed(heightmap, None, cfg, ss, window_km),
    })
}

#[cfg(test)]
mod tests {
    use super::super::drainage::c1_drainage;
    use super::super::production_upscale::upscale_from_c1;
    use super::*;
    use crate::tectonics::isostasy::IsostasyConfig;
    use std::cell::Cell;

    fn tmp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("ymir_cached_product_{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    // A cheap but real pipeline run: tiny grid + tiny target_size + erosion ON
    // (the canonical HD config), so it exercises the true code path fast.
    fn inputs()
    -> (Phase2InitParams, C1TimeLoopConfig, C1Closures, SteinSteinParams, FbmUpscaleConfig) {
        let iso = IsostasyConfig::c1_default();
        let run = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 5,
            dx: 1.0 / 32.0,
            dy: 1.0 / 32.0,
            iso_config: iso,
            drainage_max_distance: 30,
        };
        (
            Phase2InitParams::default(),
            run,
            C1Closures::default(),
            SteinSteinParams::default(),
            FbmUpscaleConfig::c1_hd_production(128),
        )
    }

    /// Island-budget parameter sweep (M1 #190). NOT a unit test — runs a real
    /// coarse-tectonic scan (minutes) and reports DISTRIBUTIONS of land-topology
    /// metrics per target_land_fraction, plus which window budgets close. Cheap
    /// by construction: one tectonic pass per (plates, clusters, seed), all tlf
    /// evaluated by re-thresholding (see `coarse_land_topology_sweep`). Raw rows
    /// dumped to CSV for re-analysis without re-running. Run:
    ///   cargo test -p ymir-core --lib sweep_island_budget -- --ignored --nocapture
    #[test]
    #[ignore]
    fn sweep_island_budget() {
        use std::f32::consts::PI;

        let grid = 64usize;
        let run = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 300,
            dx: 1.0 / grid as f64,
            dy: 1.0 / grid as f64,
            iso_config: IsostasyConfig::c1_default(),
            drainage_max_distance: 30,
        };
        let base = Phase2InitParams::default();
        let clo = C1Closures::default();
        let ss = SteinSteinParams::default();

        // Stage A (this run). Stage B if it does not close: plates {12,16,20},
        // clusters {4,6,8}, seeds 0..20, same tlf list.
        let plate_counts = [16usize];
        let cluster_counts = [6usize];
        let seeds: Vec<u64> = (0..40).collect();
        let tlfs = [0.29f32, 0.20, 0.15, 0.12];
        let margin_km = 25.0f32;
        let windows = [328.0f32, 360.0, 380.0, 400.0];

        struct Row {
            np: usize,
            cc: usize,
            seed: u64,
            tlf: f32,
            masses: usize,
            area: f32,
            frac: f32,
            wx: bool,
            wy: bool,
            bx: f32,
            by: f32,
            traverse: f32,
            disc: f32,
            compact: f32,
        }
        let mut rows: Vec<Row> = Vec::new();
        for &np in &plate_counts {
            for &cc in &cluster_counts {
                let mut init = base.clone();
                init.num_plates = np;
                init.cluster.seed_cluster_count = cc;
                for &seed in &seeds {
                    for (tlf, t) in
                        coarse_land_topology_sweep(seed, grid, &init, &run, &clo, &ss, &tlfs)
                    {
                        let traverse = t.bbox_km.0.max(t.bbox_km.1);
                        let disc = 2.0 * (t.largest_area_km2 / PI).sqrt();
                        let compact = if traverse > 0.0 { disc / traverse } else { 0.0 };
                        rows.push(Row {
                            np,
                            cc,
                            seed,
                            tlf,
                            masses: t.num_landmasses,
                            area: t.largest_area_km2,
                            frac: t.largest_area_frac,
                            wx: t.wraps_x,
                            wy: t.wraps_y,
                            bx: t.bbox_km.0,
                            by: t.bbox_km.1,
                            traverse,
                            disc,
                            compact,
                        });
                    }
                }
            }
        }

        // ── Persist raw rows (CSV) for re-analysis against new budgets. ──
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target");
        let _ = std::fs::create_dir_all(&dir);
        let csv_path = dir.join("island_sweep.csv");
        let mut csv = String::from(
            "num_plates,seed_cluster_count,seed,tlf,num_landmasses,largest_area_km2,\
             largest_area_frac,wraps_x,wraps_y,bbox_x_km,bbox_y_km,traverse_km,\
             equiv_disc_diameter_km,compactness\n",
        );
        for r in &rows {
            csv.push_str(&format!(
                "{},{},{},{:.3},{},{:.1},{:.4},{},{},{:.1},{:.1},{:.1},{:.1},{:.3}\n",
                r.np, r.cc, r.seed, r.tlf, r.masses, r.area, r.frac, r.wx, r.wy, r.bx, r.by,
                r.traverse, r.disc, r.compact
            ));
        }
        std::fs::write(&csv_path, csv).unwrap();
        eprintln!("rows: {} — CSV: {}", rows.len(), csv_path.display());

        // ── Budget table (once). ──
        eprintln!("=== budget table (margin {margin_km:.0} km, HD 8192²) ===");
        for &w in &windows {
            let budget = w - 2.0 * margin_km;
            let max_area = PI * (budget / 2.0).powi(2);
            let cell = w / 8192.0 * 1000.0;
            let flag = if !(30.0..=50.0).contains(&cell) { "  [OUTSIDE 30-50 m band]" } else { "" };
            eprintln!(
                "  window {w:.0} km → traverse budget {budget:.0} km, max compact area \
                 {max_area:.0} km², cell {cell:.1} m{flag}"
            );
        }

        // ── Distributions per target_land_fraction. ──
        let pct = |sorted: &[f32], q: f32| -> f32 {
            if sorted.is_empty() {
                return f32::NAN;
            }
            sorted[(((sorted.len() - 1) as f32) * q).round() as usize]
        };
        let mut fitting: Vec<(f32, f32, &Row)> = Vec::new(); // (tlf, min_window, row)
        for &tlf in &tlfs {
            let rt: Vec<&Row> = rows.iter().filter(|r| (r.tlf - tlf).abs() < 1e-6).collect();
            let total = rt.len();
            let wrap = rt.iter().filter(|r| r.wx || r.wy).count();
            let nonwrap: Vec<&Row> = rt.into_iter().filter(|r| !(r.wx || r.wy)).collect();
            eprintln!(
                "--- tlf {tlf:.2}: {total} configs, {wrap} wrap ({}%), {} non-wrapping",
                if total > 0 { 100 * wrap / total } else { 0 },
                nonwrap.len()
            );
            if nonwrap.is_empty() {
                continue;
            }
            let mut trav: Vec<f32> = nonwrap.iter().map(|r| r.traverse).collect();
            trav.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mut comp: Vec<f32> = nonwrap.iter().map(|r| r.compact).collect();
            comp.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let min_row = nonwrap.iter().min_by(|a, b| a.traverse.partial_cmp(&b.traverse).unwrap()).unwrap();
            eprintln!(
                "    traverse min {:.0} / median {:.0} / p90 {:.0} km; min-traverse area \
                 {:.0} km² ({:.0}%); median compactness {:.2}",
                pct(&trav, 0.0),
                pct(&trav, 0.5),
                pct(&trav, 0.9),
                min_row.area,
                min_row.frac * 100.0,
                pct(&comp, 0.5),
            );
            for &w in &windows {
                let budget = w - 2.0 * margin_km;
                let fit: Vec<&&Row> = nonwrap.iter().filter(|r| r.traverse <= budget).collect();
                let best = fit.iter().min_by(|a, b| a.traverse.partial_cmp(&b.traverse).unwrap());
                match best {
                    Some(b) => eprintln!(
                        "    window {w:.0}: {} fit — best seed {} · {} pl · {} cl → {:.0} km, \
                         {:.0} km² ({:.0}%)",
                        fit.len(),
                        b.seed,
                        b.np,
                        b.cc,
                        b.traverse,
                        b.area,
                        b.frac * 100.0
                    ),
                    None => eprintln!("    window {w:.0}: 0 fit"),
                }
            }
            // Track the highest-tlf fitting rows for the recommendation.
            for r in &nonwrap {
                if let Some(&w) = windows.iter().find(|&&w| r.traverse <= w - 2.0 * margin_km) {
                    fitting.push((tlf, w, r));
                }
            }
        }

        // ── Recommendation: keep the most land (highest tlf), then most compact. ──
        eprintln!("=== recommendation ===");
        if let Some(&(tlf, win, r)) = fitting.iter().max_by(|a, b| {
            a.0.partial_cmp(&b.0).unwrap().then(a.2.compact.partial_cmp(&b.2.compact).unwrap())
        }) {
            eprintln!(
                "RECOMMENDED: tlf {tlf:.2} · {} plates · {} clusters · seed {} → traverse \
                 {:.0} km fits window {win:.0} km (cell {:.1} m); area {:.0} km² ({:.0}%), \
                 {} masses, compactness {:.2}",
                r.np, r.cc, r.seed, r.traverse, win / 8192.0 * 1000.0, r.area, r.frac * 100.0,
                r.masses, r.compact
            );
        } else {
            let best = rows
                .iter()
                .filter(|r| !(r.wx || r.wy))
                .min_by(|a, b| a.traverse.partial_cmp(&b.traverse).unwrap());
            match best {
                Some(b) => {
                    let widest_budget = windows[windows.len() - 1] - 2.0 * margin_km; // 350 km
                    let max_area = PI * (widest_budget / 2.0).powi(2);
                    eprintln!(
                        "NO FIT at any tlf/window. Best non-wrapping: tlf {:.2} · {} pl · {} cl · \
                         seed {} → traverse {:.0} km (misses widest budget {widest_budget:.0} by \
                         {:.0} km); area {:.0} km² (over the {max_area:.0} km² compact budget by \
                         {:.0} km²)",
                        b.tlf, b.np, b.cc, b.seed, b.traverse, b.traverse - widest_budget,
                        b.area, b.area - max_area
                    );
                }
                None => eprintln!("every scanned config wraps the torus at every tlf"),
            }
        }
    }

    #[test]
    fn find_island_config_scans_deterministically() {
        let (init, run, clo, ss, _up) = inputs();
        let crit =
            IslandCriteria { window_km: 328.0, ocean_margin_km: 20.0, min_traverse_km: 50.0 };
        let scan = || {
            find_island_config(
                32,
                &init,
                &run,
                &clo,
                &ss,
                Some(0.29),
                &[8, 12],
                &[1, 3],
                &[42, 7],
                &crit,
            )
            .map(|h| (h.seed, h.num_plates, h.seed_cluster_count))
        };
        // Same inputs → same result (Some or None), and a hit must satisfy the
        // predicate on its own reported topology.
        assert_eq!(scan(), scan(), "the (plates × clusters × seeds) scan is deterministic");
        if let Some(hit) = find_island_config(
            32,
            &init,
            &run,
            &clo,
            &ss,
            Some(0.29),
            &[8, 12],
            &[1, 3],
            &[42, 7],
            &crit,
        ) {
            assert!(is_island_fit(&hit.topology, &crit), "a returned hit must fit the criteria");
        }
    }

    #[test]
    fn pipeline_miss_then_hit_round_trip() {
        let dir = tmp("hit");
        let (init, run, clo, ss, up) = inputs();
        let seed = 42u64;
        let grid = 32usize;

        // MISS — full build.
        let a = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up).unwrap();
        // HIT — must be byte-identical (same field reloaded, downstream identical).
        let b = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up).unwrap();
        assert_eq!((a.width, a.height), (b.width, b.height));
        assert_eq!(a.data, b.data, "HIT must reload byte-identical to MISS");

        // The HIT loads the same file the MISS wrote.
        let key = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &up);
        assert!(dir.join(format!("eroded_{}.raw", key.digest())).exists());
    }

    /// THE critical-completeness test on the REAL pipeline: mutate each input one
    /// at a time and assert the digest moves. A field that does NOT move the
    /// digest is a missing key input — the silent-stale bug.
    #[test]
    fn key_completeness_every_input_moves_the_digest() {
        let (init, run, clo, ss, up) = inputs();
        let seed = 42u64;
        let grid = 32usize;
        let base = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &up).digest();

        // seed
        let d = eroded_key(&tectonic_key(43, grid, &init, &run, &clo), &ss, &up).digest();
        assert_ne!(base, d, "seed must be in the key");

        // grid
        let d = eroded_key(&tectonic_key(seed, 64, &init, &run, &clo), &ss, &up).digest();
        assert_ne!(base, d, "grid must be in the key");

        // init params (R7 sub-component toggled)
        let mut init2 = init;
        init2.r7.enabled = !init.r7.enabled;
        let d = eroded_key(&tectonic_key(seed, grid, &init2, &run, &clo), &ss, &up).digest();
        assert_ne!(base, d, "init params must be in the key");

        // run config — n_steps
        let mut run2 = run.clone();
        run2.n_steps += 1;
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run2, &clo), &ss, &up).digest();
        assert_ne!(base, d, "n_steps must be in the key");

        // run config — rigid flag
        let mut run3 = run.clone();
        run3.rigid_continental_crust = !run.rigid_continental_crust;
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run3, &clo), &ss, &up).digest();
        assert_ne!(base, d, "rigid_continental_crust must be in the key");

        // run config — upstream ISOSTASY (the relief-fix surface; e.g. a future
        // continental_floor_m lives here). Move an existing iso field.
        let mut run4 = run.clone();
        run4.iso_config.max_elevation_m += 100.0;
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run4, &clo), &ss, &up).digest();
        assert_ne!(base, d, "isostasy config (relief fix) must be in the key");

        // closures — a nested closure param (Davis-Suppe O-C boost)
        let mut clo2 = clo;
        clo2.davis_suppe.oc_coupling_boost = Some(9.0);
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo2), &ss, &up).digest();
        assert_ne!(base, d, "nested closure params must be in the key");

        // ss (Stein-Stein, the upscale-level bathymetry params)
        let mut ss2 = ss;
        ss2.depth_scale_m += 1.0;
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss2, &up).digest();
        assert_ne!(base, d, "Stein-Stein params must be in the key");

        // upscale config — target_size
        let mut up2 = up.clone();
        up2.target_size = 256;
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &up2).digest();
        assert_ne!(base, d, "upscale config must be in the key");

        // upscale config — embedded erosion params
        let mut up3 = up.clone();
        if let Some(e) = up3.erosion.as_mut() {
            e.num_droplets += 1;
        }
        let d = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &up3).digest();
        assert_ne!(base, d, "embedded erosion params must be in the key");
    }

    /// Measured acceleration on the REAL 2048² product (the ~minutes of HD
    /// erosion). `#[ignore]` — run explicitly:
    /// `cargo test --release -p ymir-core --lib cached_product::tests::measure -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn measure_real_pipeline_speedup() {
        use std::time::Instant;
        let dir = tmp("measure");
        let iso = IsostasyConfig::c1_default();
        let run = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 300,
            dx: 1.0 / 64.0,
            dy: 1.0 / 64.0,
            iso_config: iso,
            drainage_max_distance: 30,
        };
        let init = Phase2InitParams::default();
        let clo = C1Closures::default();
        let ss = SteinSteinParams::default();
        let up = FbmUpscaleConfig::c1_hd_production(2048);
        let (seed, grid) = (42u64, 64usize);

        let t0 = Instant::now();
        let miss = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up).unwrap();
        let miss_t = t0.elapsed();

        let t1 = Instant::now();
        let hit = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up).unwrap();
        let hit_t = t1.elapsed();

        assert_eq!(miss.data, hit.data, "HIT must be byte-identical to MISS");
        eprintln!(
            "#cache 2048² eroded — MISS (full build) {:.2?}  |  HIT (load .raw) {:.3?}  |  speedup ×{:.0}",
            miss_t,
            hit_t,
            miss_t.as_secs_f64() / hit_t.as_secs_f64().max(1e-9)
        );
    }

    #[test]
    fn no_cache_path_is_byte_identical() {
        // The wrapper's MISS output must equal calling the generation functions
        // directly (the cache does not alter the result).
        let dir = tmp("no_cache");
        let (init, run, clo, ss, up) = inputs();
        let seed = 7u64;
        let grid = 32usize;

        let cached_out = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up).unwrap();

        let mut state = init_c1_state_phase_2_r7(grid, seed, &init);
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        run_with_closures(&mut state, &mut kin, &run, &clo, |_, _| {});
        let direct =
            upscale_from_c1(&state, &run.iso_config, &ss, &WorldSeed::new(seed), &up).heightmap;

        assert_eq!(cached_out.data, direct.data, "cache wrapper must not alter the result");
    }

    /// Whole-structure equivalence: rasters byte-equal AND the structured network
    /// (rivers / per-segment / typed lakes) equal via serde — not field-by-field
    /// happenstance.
    fn assert_drainage_equiv(a: &C1DrainageResult, b: &C1DrainageResult) {
        assert_eq!((a.width, a.height), (b.width, b.height), "dims");
        assert_eq!(a.flow.num_basins, b.flow.num_basins, "num_basins");
        assert_eq!(a.flow.filled.data, b.flow.filled.data, "filled raster");
        assert_eq!(a.flow.accumulation.data, b.flow.accumulation.data, "accumulation raster");
        assert_eq!(a.flow.direction, b.flow.direction, "direction raster");
        assert_eq!(a.flow.basins, b.flow.basins, "basins raster");
        assert_eq!(a.lake_map, b.lake_map, "lake_map raster");
        let structured = |x: &C1DrainageResult| {
            serde_json::json!({
                "rivers": x.rivers,
                "seg_km2": x.segment_drainage_km2,
                "seg_nav": x.segment_navigability,
                "lakes": x.lakes,
            })
        };
        assert_eq!(structured(a), structured(b), "structured network must be equivalent");
    }

    #[test]
    fn drainage_whole_structure_round_trip_and_chaining() {
        let dir = tmp("drainage");
        let (init, run, clo, ss, up) = inputs();
        let (seed, grid) = (42u64, 32usize);
        let h = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up).unwrap();
        let ek = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &up);
        let cfg = C1DrainageConfig::default();

        // MISS then HIT — whole structure must round-trip equivalent.
        let miss = cached_c1_drainage(&dir, &ek, &h, None, &cfg, &ss).unwrap();
        let hit = cached_c1_drainage(&dir, &ek, &h, None, &cfg, &ss).unwrap();
        assert_drainage_equiv(&miss, &hit);

        // And equal to the direct (no-cache) computation (wrapper is pure).
        let direct = c1_drainage(&h, None, &cfg, &ss);
        assert_drainage_equiv(&miss, &direct);

        // Chaining: a different eroded (upstream) key → different drainage digest.
        let ek2 = eroded_key(&tectonic_key(43, grid, &init, &run, &clo), &ss, &up);
        assert_ne!(
            drainage_key(&ek, &cfg, &ss, None).digest(),
            drainage_key(&ek2, &cfg, &ss, None).digest(),
            "upstream (eroded) change must change the drainage digest"
        );

        // Drainage config is in the key.
        let mut cfg2 = cfg.clone();
        cfg2.lake_min_area_km2 += 1.0;
        assert_ne!(
            drainage_key(&ek, &cfg, &ss, None).digest(),
            drainage_key(&ek, &cfg2, &ss, None).digest(),
            "drainage config must be in the key"
        );
        // #drainage — the climate (water-balance dependency) is in the key.
        let pp = PrecipParams::default();
        assert_ne!(
            drainage_key(&ek, &cfg, &ss, None).digest(),
            drainage_key(&ek, &cfg, &ss, Some((45.0, &pp))).digest(),
            "climate (water balance) must be in the drainage key"
        );
    }
}
