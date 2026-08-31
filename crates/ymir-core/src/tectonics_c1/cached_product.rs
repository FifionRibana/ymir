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
    ALGO_BREACH, ALGO_CLIMATE, ALGO_DRAINAGE, ALGO_HD_DRAINAGE, ALGO_TECTONICS,
    ALGO_UPSCALE_EROSION, CacheKey, RawCodec, cached, cached_fallible,
};
use crate::climate::precipitation::PrecipParams;
use crate::climate::{ClimateResult, c1_climate_placed, c1_climate_windowed};
use crate::export::raw;
use crate::grid::GridF32;
use crate::seed::WorldSeed;
use crate::terrain::flow::FlowResult;
use crate::terrain::upscale::FbmUpscaleConfig;

use super::closures::oceanic_bathymetry::SteinSteinParams;
use super::closures::volcanism::{CraterRecord, VolcanismConfig, place_edifices};
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

/// The FULL eroded key including the C-2 volcanism config. The volcanism entry is
/// added ONLY when enabled, so a disabled run keeps the exact pre-C-2 key. This is
/// the SINGLE source of truth for the eroded digest — both `cached_c1_eroded` and
/// the viz bridge's hit-detection / drainage-parent key must use it, or a
/// volcanism-on and a volcanism-off run would silently share the same drainage
/// cache entry (the terrain differs, the drainage would not).
pub fn eroded_key_full(
    tectonic: &CacheKey,
    ss: &SteinSteinParams,
    upscale_cfg: &FbmUpscaleConfig,
    volcanism: &VolcanismConfig,
) -> CacheKey {
    let key = eroded_key(tectonic, ss, upscale_cfg);
    let key = if volcanism.enabled {
        key.with("volcanism", volcanism)
            .with("volc_algo", &super::closures::volcanism::VOLCANISM_ALGO)
    } else {
        key
    };
    // C-3 lithology rides in `upscale_cfg` (skipped when disabled), so only its ALGO
    // version needs appending here, and only when enabled → disabled key unchanged.
    let key = if upscale_cfg.lithology.enabled {
        key.with("litho_algo", &super::closures::lithology::LITHOLOGY_ALGO)
    } else {
        key
    };
    // C-3b fracture — same pattern: config rides in `upscale_cfg`, append the ALGO only
    // when enabled so a disabled key is byte-identical to pre-C-3b.
    if upscale_cfg.fracture.enabled {
        key.with("fracture_algo", &super::closures::fracture::FRACTURE_ALGO)
    } else {
        key
    }
}

/// The cached eroded HD product: the heightmap AND the volcanic craters placed in
/// it, bundled so they are computed and cached TOGETHER. A hit that returned the
/// terrain without its craters would mistype every crater lake silently — the
/// bundle makes that mismatch impossible by construction (C-2, ADR 0001).
pub struct ErodedProduct {
    /// The eroded HD heightmap (post-FBM, post-volcanism, post-erosion, norm).
    pub heightmap: GridF32,
    /// Craters placed by the volcanism closure, in HD pixel coords. Empty when
    /// volcanism is disabled.
    pub craters: Vec<CraterRecord>,
}

impl RawCodec for ErodedProduct {
    fn shape(&self) -> Value {
        serde_json::json!({ "width": self.heightmap.width, "height": self.heightmap.height })
    }
    fn write_raw(&self, stem: &Path) -> Result<(), String> {
        self.heightmap.save_raw(&stem.with_extension("raw"))?;
        // Craters ride in a sibling JSON so they travel with the terrain. Always
        // written (an empty list when volcanism is off) so the reader can require
        // it — under ALGO v2 an old bare-.raw entry has a different digest and is
        // never read here.
        let json = serde_json::to_string(&self.craters)
            .map_err(|e| format!("crater serialize error: {e}"))?;
        std::fs::write(stem.with_extension("craters.json"), json)
            .map_err(|e| format!("crater write error: {e}"))
    }
    fn read_raw(stem: &Path, shape: &Value) -> Result<Self, String> {
        let w = shape["width"].as_u64().ok_or("sidecar: missing width")? as usize;
        let h = shape["height"].as_u64().ok_or("sidecar: missing height")? as usize;
        let heightmap = GridF32::load_raw(&stem.with_extension("raw"), w, h)?;
        let cj = std::fs::read_to_string(stem.with_extension("craters.json"))
            .map_err(|e| format!("crater read error: {e}"))?;
        let craters: Vec<CraterRecord> =
            serde_json::from_str(&cj).map_err(|e| format!("crater parse error: {e}"))?;
        Ok(ErodedProduct { heightmap, craters })
    }
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
    volcanism: &VolcanismConfig,
) -> Result<ErodedProduct, String> {
    cached_c1_eroded_with_progress(
        cache_dir,
        seed,
        grid,
        init,
        run,
        closures,
        ss,
        upscale_cfg,
        volcanism,
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
    volcanism: &VolcanismConfig,
    progress: &mut dyn FnMut(EroProgress),
    cancel: &dyn Fn() -> bool,
) -> Result<ErodedProduct, String> {
    // Volcanism enters the cache key ONLY when enabled, so a disabled run keeps the
    // pre-C-2 terrain key (the ALGO v2 bump already forces the one-time recompute
    // for the payload-format change). Enabling or changing it → new digest →
    // recompute, so a config change can never serve a stale terrain.
    let key =
        eroded_key_full(&tectonic_key(seed, grid, init, run, closures), ss, upscale_cfg, volcanism);
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
        // C-2: derive the edifices from the FINAL tectonic state + kinematics (both
        // live only here, on the miss path — hence bundling craters with the
        // terrain rather than recomputing at typing time). `torus_km = domain_km`
        // (geometric span). Empty when disabled.
        let edifices =
            place_edifices(&state, &kin, &WorldSeed::new(seed), volcanism.domain_km, volcanism);
        if volcanism.enabled {
            use super::closures::volcanism::VolcanoSetting;
            let n = |s: VolcanoSetting| edifices.iter().filter(|e| e.setting == s).count();
            eprintln!(
                "[HD volcanism] edifices: {} arc + {} hotspot + {} rift = {} (domain {:.0} km)",
                n(VolcanoSetting::Arc),
                n(VolcanoSetting::Hotspot),
                n(VolcanoSetting::Rift),
                edifices.len(),
                volcanism.domain_km,
            );
        }
        let (up, craters) = upscale_from_c1_with_progress(
            &state,
            &run.iso_config,
            ss,
            &WorldSeed::new(seed),
            upscale_cfg,
            &edifices,
            volcanism,
            Some(&kin), // C-3b: kin for the fracture field (alive here on the miss path)
            progress,
            cancel,
        );
        if cancel() {
            return Err("cancelled".to_string());
        }
        Ok(ErodedProduct { heightmap: up.heightmap, craters })
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

/// Sweep the normalized coarse field of ONE tectonic config across many target
/// land fractions, running the (expensive) coarse tectonic pass ONCE and
/// producing each `tlf`'s field by re-thresholding the same raw altitude (cheap).
/// Returns `(tlf, normalized_field)` per fraction; the caller derives topology /
/// border-clean metrics. Window/margin budgets are a POST-HOC predicate on those
/// metrics — not swept here. Deterministic in the tectonic key.
#[allow(clippy::too_many_arguments)]
pub fn coarse_normalized_sweep(
    seed: u64,
    grid: usize,
    init: &Phase2InitParams,
    run: &C1TimeLoopConfig,
    closures: &C1Closures,
    ss: &SteinSteinParams,
    target_land_fractions: &[f32],
) -> Vec<(f32, GridF32)> {
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, run, closures, |_, _| {});
    let raw = c1_coarse_raw_altitude(&state, &run.iso_config, ss);
    target_land_fractions.iter().map(|&f| (f, c1_normalize_coarse(raw.clone(), Some(f)))).collect()
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
            "seg_q": self.segment_discharge_m3s,
            "seg_width": self.segment_width_m,
            "seg_profile": self.segment_profile_m,
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
        // seg_q / seg_width / seg_profile: default to empty on older sidecars (forward-compat).
        let segment_discharge_m3s = shape
            .get("seg_q")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .unwrap_or_default();
        let segment_width_m = shape
            .get("seg_width")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .unwrap_or_default();
        let segment_profile_m = shape
            .get("seg_profile")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .unwrap_or_default();
        let lakes = serde_json::from_value(shape["lakes"].clone())
            .map_err(|e| format!("drainage sidecar lakes: {e}"))?;

        Ok(C1DrainageResult {
            flow: FlowResult { filled, direction, accumulation, basins, num_basins },
            rivers,
            segment_drainage_km2,
            segment_navigability,
            segment_discharge_m3s,
            segment_width_m,
            segment_profile_m,
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

// ── Climate cache (#190) ─────────────────────────────────────────────────────
// The climate (temperature + precipitation) is a pure function of the eroded relief,
// the centre latitude and the climatic span. It is recomputed on every HD run today;
// cache it so re-selecting the SAME (latitude, span) is a HIT.

impl RawCodec for ClimateResult {
    fn shape(&self) -> Value {
        serde_json::json!({ "width": self.temperature.width, "height": self.temperature.height })
    }
    fn write_raw(&self, stem: &Path) -> Result<(), String> {
        self.temperature.save_raw(&part(stem, "temp"))?;
        self.precipitation.save_raw(&part(stem, "precip"))?;
        Ok(())
    }
    fn read_raw(stem: &Path, shape: &Value) -> Result<Self, String> {
        let w = shape["width"].as_u64().ok_or("climate sidecar: width")? as usize;
        let h = shape["height"].as_u64().ok_or("climate sidecar: height")? as usize;
        Ok(ClimateResult {
            temperature: GridF32::load_raw(&part(stem, "temp"), w, h)?,
            precipitation: GridF32::load_raw(&part(stem, "precip"), w, h)?,
        })
    }
}

/// Cache key for the climate step — chained onto the eroded key (relief change → new
/// digest) plus the centre latitude, the climatic span (`None` = geographic default),
/// and the precip params. Any of these changing re-invalidates the climate.
pub fn climate_key(
    eroded: &CacheKey,
    latitude_deg: f32,
    span_deg: Option<f32>,
    pp: &PrecipParams,
    window_km: f32,
) -> CacheKey {
    let mut k = CacheKey::derived_from(eroded)
        .with("climate_lat", &latitude_deg)
        .with("climate_span", &span_deg)
        .with("climate_precip", pp);
    if window_km != C1_DOMAIN_KM {
        k = k.with("climate_window_km", &window_km);
    }
    k.algo(ALGO_CLIMATE)
}

/// Compute-or-load the climate. `span` = `Some` uses the placed (explicit-span) path,
/// `None` the windowed path — byte-identical to the direct calls in the HD pipeline.
#[allow(clippy::too_many_arguments)]
pub fn cached_c1_climate(
    cache_dir: &Path,
    eroded: &CacheKey,
    heightmap: &GridF32,
    ss: &SteinSteinParams,
    latitude_deg: f32,
    span_deg: Option<f32>,
    pp: &PrecipParams,
    window_km: f32,
) -> Result<ClimateResult, String> {
    let key = climate_key(eroded, latitude_deg, span_deg, pp, window_km);
    cached(cache_dir, "climate", &key, || match span_deg {
        Some(span) => c1_climate_placed(heightmap, ss, latitude_deg, span, pp, window_km),
        None => c1_climate_windowed(heightmap, ss, latitude_deg, pp, window_km),
    })
}

/// Cache key for the relief-v3 CONDITIONED eroded field (the `breach_monotone` output). Chained
/// onto the PRE-BREACH drainage key (which itself folds the eroded key + drainage config + ss), so
/// any change to the eroded relief OR the pre-breach lakes OR the breach code re-invalidates it.
/// The conditioned field is a pure `GridF32` — cache it with the plain `GridF32` codec.
pub fn conditioned_eroded_key(prebreach_key: &CacheKey) -> CacheKey {
    CacheKey::derived_from(prebreach_key).algo(ALGO_BREACH)
}

// ── HD relief-v3 drainage bundle cache (#190) ────────────────────────────────
// The FINAL drainage (after the below-sea water-balance merge, river clip, spillway
// append and geographic-scale) plus the wetland mask — the climate-dependent,
// per-run-recomputed tail of the HD pipeline. The relief-v3 path bypassed the plain
// drainage cache (it needs the breached field + climate + carried pre-breach lakes and
// the whole post-process), so this bundles the final product under a COMPLETE key.

/// The final HD drainage + wetland mask, cached as one unit.
pub struct HdDrainageBundle {
    pub drainage: C1DrainageResult,
    pub wetland: Vec<u8>,
}

impl RawCodec for HdDrainageBundle {
    fn shape(&self) -> Value {
        // Reuse the drainage shape (dims + structured network) verbatim; the wetland is a
        // plain w*h u8 raster restored from those dims.
        self.drainage.shape()
    }
    fn write_raw(&self, stem: &Path) -> Result<(), String> {
        self.drainage.write_raw(stem)?;
        raw::save_u8(&part(stem, "wetland"), &self.wetland)?;
        Ok(())
    }
    fn read_raw(stem: &Path, shape: &Value) -> Result<Self, String> {
        let drainage = C1DrainageResult::read_raw(stem, shape)?;
        let n = drainage.width * drainage.height;
        let wetland = raw::load_u8(&part(stem, "wetland"))?;
        if wetland.len() != n {
            return Err(format!("hd_drainage wetland: expected {n}, got {}", wetland.len()));
        }
        Ok(HdDrainageBundle { drainage, wetland })
    }
}

/// Cache key for the HD relief-v3 drainage bundle. COMPLETE by construction — it folds
/// every input the final drainage depends on: the eroded key (tectonics + upscale +
/// relief), the drainage config + Stein-Stein, the centre latitude + precip params +
/// climatic span (they set the discharge via climate), the geographic scale ratio (it
/// scales discharge / navigability / width / area), and the window. Missing any one
/// would be the stale-cache bug — see the key-completeness note in `cache.rs`.
#[allow(clippy::too_many_arguments)]
pub fn hd_drainage_key(
    eroded: &CacheKey,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    latitude_deg: f32,
    pp: &PrecipParams,
    span_deg: Option<f32>,
    geo_scale_ratio: f32,
    window_km: f32,
) -> CacheKey {
    CacheKey::derived_from(eroded)
        .with("hd_dr_cfg", cfg)
        .with("hd_dr_ss", ss)
        .with("hd_dr_lat", &latitude_deg)
        .with("hd_dr_precip", pp)
        .with("hd_dr_span", &span_deg)
        .with("hd_dr_geo_scale", &geo_scale_ratio)
        .with("hd_dr_window_km", &window_km)
        .algo(ALGO_HD_DRAINAGE)
}

#[cfg(test)]
mod tests {
    use super::super::drainage::c1_drainage;
    use super::super::land_topology::{evaluate_island, harvest_islands};
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

    /// C-2 PERMANENT GUARD (author review #2): enabling volcanism MUST invalidate
    /// the eroded cache key AND change the terrain — never silently serve the
    /// volcano-free entry. A key that failed to invalidate would report success
    /// while every downstream measurement ran on a terrain with no volcanoes. So
    /// this is tested, not announced: generate WITH and WITHOUT, assert the digests
    /// differ and the two heightmaps differ; the OFF path carries no craters.
    #[test]
    fn volcanism_invalidates_key_and_changes_terrain() {
        let dir = tmp("volcanism_invalidate");
        let (init, run, clo, ss, up) = inputs();
        let (seed, grid) = (42u64, 32usize);
        let off = VolcanismConfig::default(); // disabled
        let on = VolcanismConfig { enabled: true, domain_km: 400.0, ..Default::default() };

        // The full eroded digest must change when volcanism is enabled.
        let tk = tectonic_key(seed, grid, &init, &run, &clo);
        let ka = eroded_key_full(&tk, &ss, &up, &off);
        let kb = eroded_key_full(&tk, &ss, &up, &on);
        assert_ne!(ka.digest(), kb.digest(), "enabled volcanism must change the eroded digest");

        // …and the terrain it addresses must actually differ (edifices injected).
        let a = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up, &off).unwrap();
        let b = cached_c1_eroded(&dir, seed, grid, &init, &run, &clo, &ss, &up, &on).unwrap();
        assert!(a.craters.is_empty(), "disabled volcanism must place no craters");
        let changed = a
            .heightmap
            .data
            .iter()
            .zip(&b.heightmap.data)
            .filter(|(x, y)| (**x - **y).abs() > 1e-6)
            .count();
        assert!(
            changed > 0,
            "volcanism ON must change the terrain; diff 0 would mean the key did not invalidate"
        );
    }

    /// Island-budget sweep with the BORDER-CLEAN predicate (M1 #190 reframe). NOT
    /// a unit test — runs a real coarse-tectonic scan (minutes). The window is
    /// SIZED to the largest landmass (`traverse + 2·margin`); a config is accepted
    /// iff the largest mass does not wrap, the window cell size is in the 30–50 m
    /// band, and NO land (of any mass) lies on the window border ring. Compactness
    /// is reported, never selected on. Cheap: one tectonic pass per (plates,
    /// clusters, seed), all tlf by re-thresholding. Raw rows → CSV. Run:
    ///   cargo test -p ymir-core --lib sweep_island_budget -- --ignored --nocapture
    #[test]
    #[ignore]
    fn sweep_island_budget() {
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

        // Option 3: push tlf lower to shrink the ISOLATED island under the 360 km
        // in-band ceiling — but guard on AREA, not just traverse (a smaller tlf can
        // give a tinier, more tentacular island: barely-reduced window for far less
        // land). Floor at 25 000 km² (seed 8 is 31 232) and report compactness so
        // land-vs-resolution is arbitrated on numbers. Low cc (isolation), same as
        // Stage B; fragmentation is the CAUSE of the satellite-clip, not the cure.
        let plate_counts = [12usize, 16, 20];
        let cluster_counts = [1usize, 2, 3];
        let seeds: Vec<u64> = (0..20).collect();
        let tlfs = [0.10f32, 0.08, 0.06];
        let margin_km = 25.0f32;
        let ring_cells = 1usize;
        let min_area_km2 = 25_000.0f32;
        let km_per_cell = C1_DOMAIN_KM / grid as f32;

        struct Row {
            np: usize,
            cc: usize,
            seed: u64,
            tlf: f32,
            masses: usize,
            area: f32,
            frac: f32,
            traverse: f32,
            window_km: f32,
            m_per_cell: f32,
            res_ok: bool,
            wx: bool,
            wy: bool,
            clean: bool,
            sats: usize,
            compact: f32,
        }
        impl Row {
            fn accepted(&self) -> bool {
                !self.wx && !self.wy && self.res_ok && self.clean
            }
        }

        let mut rows: Vec<Row> = Vec::new();
        for &np in &plate_counts {
            for &cc in &cluster_counts {
                let mut init = base.clone();
                init.num_plates = np;
                init.cluster.seed_cluster_count = cc;
                for &seed in &seeds {
                    for (tlf, field) in
                        coarse_normalized_sweep(seed, grid, &init, &run, &clo, &ss, &tlfs)
                    {
                        let e = evaluate_island(&field, C1_SEA_LEVEL_NORM, margin_km, ring_cells);
                        rows.push(Row {
                            np,
                            cc,
                            seed,
                            tlf,
                            masses: e.topo.num_landmasses,
                            area: e.topo.largest_area_km2,
                            frac: e.topo.largest_area_frac,
                            traverse: e.topo.bbox_km.0.max(e.topo.bbox_km.1),
                            window_km: e.window_km,
                            m_per_cell: e.m_per_cell,
                            res_ok: e.resolution_ok,
                            wx: e.topo.wraps_x,
                            wy: e.topo.wraps_y,
                            clean: e.border_clean,
                            sats: e.satellites_inside,
                            compact: e.compactness,
                        });
                    }
                }
            }
        }

        // ── Persist raw rows (CSV). ──
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target");
        let _ = std::fs::create_dir_all(&dir);
        let csv_path = dir.join("island_sweep.csv");
        let mut csv = String::from(
            "num_plates,seed_cluster_count,seed,tlf,num_landmasses,largest_area_km2,\
             largest_area_frac,traverse_km,window_km,m_per_cell,resolution_ok,wraps_x,wraps_y,\
             border_clean,satellites_inside,compactness\n",
        );
        for r in &rows {
            csv.push_str(&format!(
                "{},{},{},{:.3},{},{:.1},{:.4},{:.1},{:.1},{:.2},{},{},{},{},{},{:.3}\n",
                r.np,
                r.cc,
                r.seed,
                r.tlf,
                r.masses,
                r.area,
                r.frac,
                r.traverse,
                r.window_km,
                r.m_per_cell,
                r.res_ok,
                r.wx,
                r.wy,
                r.clean,
                r.sats,
                r.compact
            ));
        }
        std::fs::write(&csv_path, csv).unwrap();
        eprintln!(
            "rows: {} — CSV: {}\nwindow SIZED to island (traverse + 2·{margin_km:.0} km); ring \
             {ring_cells} cell = {:.0} km; resolution band 30-50 m/cell → window ∈ [246,400] km",
            rows.len(),
            csv_path.display(),
            ring_cells as f32 * km_per_cell
        );

        // ── Distributions per target_land_fraction. ──
        for &tlf in &tlfs {
            let rt: Vec<&Row> = rows.iter().filter(|r| (r.tlf - tlf).abs() < 1e-6).collect();
            let total = rt.len();
            let wrap = rt.iter().filter(|r| r.wx || r.wy).count();
            let nonwrap = total - wrap;
            let clean: Vec<&&Row> = rt.iter().filter(|r| r.accepted()).collect();
            let clean_floor = clean.iter().filter(|r| r.area >= min_area_km2).count();
            eprintln!(
                "--- tlf {tlf:.2}: {total} configs, {wrap} wrap, {nonwrap} non-wrap, \
                 {} BORDER-CLEAN ({} ≥ {:.0}k km²)",
                clean.len(),
                clean_floor,
                min_area_km2 / 1000.0
            );
            if let Some(b) = clean.iter().max_by(|a, b| a.area.partial_cmp(&b.area).unwrap()) {
                eprintln!(
                    "    largest clean: seed {} → area {:.0} km² ({:.0}%), window {:.0} km \
                     ({:.1} m/cell), {} satellites, compactness {:.2}",
                    b.seed,
                    b.area,
                    b.frac * 100.0,
                    b.window_km,
                    b.m_per_cell,
                    b.sats,
                    b.compact
                );
            }
        }

        // ── Recommendation: BORDER-CLEAN first, then LARGEST area, then finest cell. ──
        eprintln!("=== recommendation ===");
        // Border-clean AND in-band AND non-wrap AND area ≥ floor; then LARGEST area.
        let best = rows.iter().filter(|r| r.accepted() && r.area >= min_area_km2).max_by(|a, b| {
            a.area
                .partial_cmp(&b.area)
                .unwrap()
                .then(b.m_per_cell.partial_cmp(&a.m_per_cell).unwrap())
        });
        match best {
            Some(r) => {
                eprintln!(
                    "RECOMMENDED (reproduce): num_plates={} seed_cluster_count={} seed={} \
                     target_land_fraction={:.2} window_km={:.0} margin_km={:.0} ring_cells={ring_cells}",
                    r.np, r.cc, r.seed, r.tlf, r.window_km, margin_km
                );
                eprintln!(
                    "  → largest island {:.0} km² ({:.0}%), traverse {:.0} km, {:.1} m/cell, \
                     {} satellites inside, compactness {:.2} (reported only)",
                    r.area,
                    r.frac * 100.0,
                    r.traverse,
                    r.m_per_cell,
                    r.sats,
                    r.compact
                );

                // ── Step 4: multi-continent harvest on the recommended config. ──
                let mut init = base.clone();
                init.num_plates = r.np;
                init.cluster.seed_cluster_count = r.cc;
                let fields = coarse_normalized_sweep(r.seed, grid, &init, &run, &clo, &ss, &tlfs);
                let field = &fields.iter().find(|(t, _)| (t - r.tlf).abs() < 1e-6).unwrap().1;
                let cands = harvest_islands(field, C1_SEA_LEVEL_NORM, margin_km, ring_cells);
                let placeable: Vec<_> = cands.iter().filter(|c| c.placeable()).collect();
                eprintln!(
                    "=== harvest === {} landmasses, {} independently placeable continents:",
                    cands.len(),
                    placeable.len()
                );
                for c in placeable.iter().take(12) {
                    eprintln!(
                        "  area {:.0} km², traverse {:.0} km, window {:.0} km ({:.1} m/cell)",
                        c.area_km2, c.traverse_km, c.window_km, c.m_per_cell
                    );
                }
            }
            None => {
                // Report the closest miss and WHY (resolution vs border).
                let nonwrap: Vec<&Row> = rows.iter().filter(|r| !(r.wx || r.wy)).collect();
                let dirty = nonwrap.iter().filter(|r| r.res_ok && !r.clean).count();
                let too_fine = nonwrap.iter().filter(|r| r.m_per_cell < 30.0).count();
                let too_coarse = nonwrap.iter().filter(|r| r.m_per_cell > 50.0).count();
                let min_sats =
                    nonwrap.iter().filter(|r| r.res_ok).map(|r| r.sats).min().unwrap_or(0);
                eprintln!(
                    "NO BORDER-CLEAN config. Of {} non-wrapping: {} in-band but border-dirty \
                     (min {} satellites clip the window), {} too fine (<30 m), {} too coarse \
                     (>50 m).",
                    nonwrap.len(),
                    dirty,
                    min_sats,
                    too_fine,
                    too_coarse
                );
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
        let a = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
        // HIT — must be byte-identical (same field reloaded, downstream identical).
        let b = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
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
        let miss = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
        let miss_t = t0.elapsed();

        let t1 = Instant::now();
        let hit = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
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

        let cached_out = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;

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
        let h = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
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

    /// #190 — the climate + HD-drainage-bundle caches: KEY COMPLETENESS (the load-bearing safety
    /// property — a missing key input serves stale terrain) and codec round-trip (a HIT returns the
    /// same payload without recomputing). The turbofish `|| panic!` on the HIT proves the cache does
    /// NOT recompute.
    #[test]
    fn climate_and_hd_drainage_cache_key_and_round_trip() {
        let dir = tmp("hd_cache");
        let (init, run, clo, ss, up) = inputs();
        let (seed, grid) = (42u64, 32usize);
        let h = cached_c1_eroded(
            &dir,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &up,
            &VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
        let ek = eroded_key(&tectonic_key(seed, grid, &init, &run, &clo), &ss, &up);
        let cfg = C1DrainageConfig::default();
        let pp = PrecipParams::default();
        let w = C1_DOMAIN_KM;

        // Climate: MISS then HIT return the same grids.
        let cm = cached_c1_climate(&dir, &ek, &h, &ss, 45.0, Some(10.0), &pp, w).unwrap();
        let cm2 = cached_c1_climate(&dir, &ek, &h, &ss, 45.0, Some(10.0), &pp, w).unwrap();
        assert_eq!(cm.temperature.data, cm2.temperature.data);
        assert_eq!(cm.precipitation.data, cm2.precipitation.data);

        // Climate key completeness: latitude AND span both change the digest.
        assert_ne!(
            climate_key(&ek, 45.0, Some(10.0), &pp, w).digest(),
            climate_key(&ek, 25.0, Some(10.0), &pp, w).digest(),
            "latitude must be in the climate key"
        );
        assert_ne!(
            climate_key(&ek, 45.0, Some(10.0), &pp, w).digest(),
            climate_key(&ek, 45.0, Some(40.0), &pp, w).digest(),
            "climatic span must be in the climate key"
        );

        // HD-drainage-bundle key completeness: span AND geo-scale both change the digest.
        assert_ne!(
            hd_drainage_key(&ek, &cfg, &ss, 45.0, &pp, Some(10.0), 7.5, w).digest(),
            hd_drainage_key(&ek, &cfg, &ss, 45.0, &pp, Some(40.0), 7.5, w).digest(),
            "climatic span must be in the HD-drainage key"
        );
        assert_ne!(
            hd_drainage_key(&ek, &cfg, &ss, 45.0, &pp, Some(10.0), 7.5, w).digest(),
            hd_drainage_key(&ek, &cfg, &ss, 45.0, &pp, Some(10.0), 1.0, w).digest(),
            "geographic scale ratio must be in the HD-drainage key"
        );

        // HD-drainage bundle codec round-trip: MISS writes (drainage rasters + wetland), HIT reads —
        // the closure MUST NOT run on the HIT (it panics if it does).
        let dr = cached_c1_drainage(&dir, &ek, &h, None, &cfg, &ss).unwrap();
        let n = dr.width * dr.height;
        let key = hd_drainage_key(&ek, &cfg, &ss, 45.0, &pp, Some(10.0), 7.5, w);
        let miss = cached::<HdDrainageBundle>(&dir, "hd_drainage", &key, || HdDrainageBundle {
            drainage: cached_c1_drainage(&dir, &ek, &h, None, &cfg, &ss).unwrap(),
            wetland: vec![1u8; n],
        })
        .unwrap();
        let hit = cached::<HdDrainageBundle>(&dir, "hd_drainage", &key, || {
            panic!("HD-drainage HIT must not recompute the bundle")
        })
        .unwrap();
        assert_eq!(miss.drainage.width, hit.drainage.width);
        assert_eq!(miss.drainage.lake_map, hit.drainage.lake_map);
        assert_eq!(hit.wetland, vec![1u8; n], "wetland mask must round-trip");
    }
}
