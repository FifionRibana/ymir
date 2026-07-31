//! HD production pipeline driven on the C1 worker thread (UI rewrite
//! step b/5 — the foundation).
//!
//! The coarse tectonic loop (`RunBaseline` / `RunWorkflow`) renders the
//! live 64² S̃/age/plate/craton fields. This module drives the FULL HD
//! production chain — the same `ymir-core` functions the production /
//! test pipelines use — on the background worker:
//!
//! ```text
//!   eroded  (tectonics → upscale → erosion → bathymetry, CACHED)
//!   climate (temperature + precipitation, computed)
//!   drainage (rivers + lakes + water balance, CACHED)
//!   biomes  (Whittaker classification, computed)
//! ```
//!
//! Each phase emits `HdPhaseStarted` → `HdPhaseDone { regime, elapsed }`
//! events (see [`super::events`]). Per ÉTAPE 0 every phase is an OPAQUE
//! block (no progress callback survives `cached_c1_eroded` / `c1_drainage`),
//! so the UI shows an INDETERMINATE waiter per phase — not an N/total bar.
//! The cache regime (HIT vs MISS) is detected by probing the sidecar file
//! BEFORE the call (exact, not time-inferred): a HIT reloads in ~ms, a MISS
//! computes (erosion at 2048² is the slow one). Re-running the same continent
//! → every cached phase HITs → near-instant.
//!
//! `ymir-core` is CALLED, never re-implemented (the legacy-viz mistake the
//! v2 removal undid).

use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use crossbeam_channel::Sender;

use ymir_core::cache::default_cache_dir;
use ymir_core::climate::biomes::Biome;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::climate::{ClimateResult, c1_biomes, c1_climate};
use ymir_core::export::container::{ContinentMeta, ContinentWriter, Grid};
use ymir_core::export::{hydro, vector};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::{
    cached_c1_drainage, cached_c1_eroded, cached_c1_eroded_with_progress, drainage_key, eroded_key,
    tectonic_key,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, C1DrainageResult};
use ymir_core::tectonics_c1::production_upscale::{C1_DOMAIN_KM, EroProgress};
use ymir_core::tectonics_c1::time_loop::C1TimeLoopConfig;
use ymir_core::terrain::upscale::FbmUpscaleConfig;

use super::events::C1Event;
use super::spec::C1RunSpec;

/// HD-specific run parameters (the coarse tectonic params live in
/// [`C1RunSpec`]). `target_size` is the HD grid edge (2048 = production);
/// `latitude_deg` drives the climate band (45° = the production anchor the
/// drainage/climate tiles use).
#[derive(Clone, Debug, PartialEq)]
pub struct HdParams {
    pub target_size: usize,
    pub latitude_deg: f32,
    /// When `Some(dir)`, after the biome phase the run writes a v1 `.ymir`
    /// delivery container under `dir` (see [`ymir_core::export::container`]).
    /// `None` = no export. Explicit opt-in — the pipeline NEVER auto-exports
    /// (WP-0). The container directory is `dir/<name>.ymir/`.
    pub export_dir: Option<PathBuf>,
}

impl Default for HdParams {
    fn default() -> Self {
        Self { target_size: 2048, latitude_deg: 45.0, export_dir: None }
    }
}

/// The HD pipeline phases, in execution order. The former single `Eroded`
/// phase is split into `Tectonic` → `Relief` → `Erosion` (suite e) so the frieze
/// shows the long erosion progressing; bathymetry is folded into `Erosion`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HdPhase {
    /// Coarse thin-sheet tectonics (`run_with_closures`) — determinate N/steps.
    Tectonic,
    /// HD relief synthesis (bicubic + FBM upscale) — opaque waiter.
    Relief,
    /// HD hydraulic erosion (+ bathymetry) — determinate % (batch callback).
    Erosion,
    /// Temperature + precipitation.
    Climate,
    /// Rivers + lakes + water balance (cached as "drainage").
    Drainage,
    /// Whittaker biome classification.
    Biomes,
}

impl HdPhase {
    pub fn label(self) -> &'static str {
        match self {
            HdPhase::Tectonic => "tectonique",
            HdPhase::Relief => "relief",
            HdPhase::Erosion => "érosion",
            HdPhase::Climate => "climat",
            HdPhase::Drainage => "drainage",
            HdPhase::Biomes => "biomes",
        }
    }
}

/// Outcome of a phase w.r.t. the content-addressed cache.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CacheRegime {
    /// Reloaded from `.ymir_cache/` (the sidecar existed) — near-instant.
    Hit,
    /// Computed and written (no sidecar) — the slow path.
    Miss,
    /// Not a cached phase (climate / biomes are cheap, always computed).
    Computed,
}

impl CacheRegime {
    pub fn label(self) -> &'static str {
        match self {
            CacheRegime::Hit => "depuis le cache",
            CacheRegime::Miss => "calculé",
            CacheRegime::Computed => "calculé",
        }
    }
}

/// The full HD product — every layer the UI inspects. Carried to the UI
/// inside an [`Arc`] (so [`C1Event`] stays `Clone`/`Debug` without
/// requiring `C1DrainageResult: Clone`, which it is not).
pub struct HdResult {
    pub width: usize,
    pub height: usize,
    /// Final relief + bathymetry (post-erosion, normalised altitude).
    pub eroded: GridF32,
    /// Temperature grid (°C).
    pub temperature: GridF32,
    /// Precipitation grid (internal units → `precip_mm_per_year`).
    pub precipitation: GridF32,
    /// Rivers + lakes + per-segment navigability + water balance.
    pub drainage: C1DrainageResult,
    /// Per-cell Whittaker biome (row-major `width × height`).
    pub biomes: Vec<Biome>,
}

impl fmt::Debug for HdResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Manual Debug: the heavy fields (and `C1DrainageResult`) are not
        // `Debug`. Print a compact summary so `C1Event` can derive `Debug`.
        f.debug_struct("HdResult")
            .field("width", &self.width)
            .field("height", &self.height)
            .field("rivers", &self.drainage.rivers.segments.len())
            .field("lakes", &self.drainage.lakes.len())
            .finish()
    }
}

/// A finished phase record (for the UI's per-phase HIT/MISS log).
#[derive(Clone, Copy, Debug)]
pub struct HdPhaseRecord {
    pub phase: HdPhase,
    pub regime: CacheRegime,
    pub elapsed: Duration,
}

/// UI-side state of the HD pipeline, updated by `poll_c1_events` from the
/// HD event stream. Independent of the coarse tectonic `C1RunState` so the
/// live gallery view and the HD generation coexist.
#[derive(Clone, Default)]
pub enum HdState {
    #[default]
    Idle,
    Running {
        params: HdParams,
        /// Phase currently in flight (an indeterminate waiter), if any.
        current: Option<HdPhase>,
        /// `(done, total)` of the current phase when it is DETERMINATE (Tectonic
        /// steps, Erosion batches) → the frieze draws a real bar; `None` → waiter.
        progress: Option<(usize, usize)>,
        /// Phases finished so far, with their HIT/MISS regime + timing.
        done: Vec<HdPhaseRecord>,
    },
    Completed {
        result: Arc<HdResult>,
        total: Duration,
        done: Vec<HdPhaseRecord>,
    },
    Failed {
        error: String,
    },
}

/// Build the production tectonic time-loop config for `spec` (matches the
/// gallery worker + the production/test pipelines so cache keys align).
fn hd_run_config(spec: &C1RunSpec) -> C1TimeLoopConfig {
    C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: spec.n_steps,
        dx: 1.0 / spec.grid_size as f64,
        dy: 1.0 / spec.grid_size as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: spec.drainage_max_distance,
    }
}

/// Drive the HD chain on the worker thread, emitting per-phase events.
/// Cancellable BETWEEN phases (the cached/opaque calls have no interior
/// cancel hook — adding one would reopen `ymir-core`). On any core error,
/// emits `HdFailed` and returns.
pub fn run_hd(spec: &C1RunSpec, params: &HdParams, tx: &Sender<C1Event>, cancel: &Arc<AtomicBool>) {
    cancel.store(false, Ordering::Relaxed);
    let _ = tx.send(C1Event::HdStarted { spec: spec.clone(), params: params.clone() });
    let t_all = Instant::now();

    let run = hd_run_config(spec);
    let ss = SteinSteinParams::default();
    let upscale = FbmUpscaleConfig::c1_hd_production(params.target_size);
    let pp = PrecipParams::default();
    let dcfg = C1DrainageConfig::default();
    let cache_dir = default_cache_dir();
    let lat = params.latitude_deg;

    // Cache keys (for HIT/MISS detection via sidecar existence).
    let tkey = tectonic_key(spec.seed, spec.grid_size, &spec.init_params, &run, &spec.closures);
    let ekey = eroded_key(&tkey, &ss, &upscale);
    let sidecar_exists = |step: &str, digest: String| -> bool {
        cache_dir.join(format!("{step}_{digest}.json")).exists()
    };

    macro_rules! bail_if_cancelled {
        () => {
            if cancel.load(Ordering::Relaxed) {
                let _ = tx.send(C1Event::HdFailed { error: "annulé".to_string() });
                return;
            }
        };
    }

    // ── Phase 1: eroded, split into Tectonic → Relief → Erosion (bathy folded).
    // On a HIT the whole "eroded" .raw reloads at once (no sub-step runs) → mark
    // the three nodes done-from-cache instantly. On a MISS the core progress
    // callback drives the sub-nodes (Tectonic N/steps, Relief waiter, Erosion %),
    // with mid-erosion cancel (a cancelled build is NOT cached).
    bail_if_cancelled!();
    let eroded_hit = sidecar_exists("eroded", ekey.digest());
    let eroded = if eroded_hit {
        for phase in [HdPhase::Tectonic, HdPhase::Relief, HdPhase::Erosion] {
            let _ = tx.send(C1Event::HdPhaseStarted { phase });
            let _ = tx.send(C1Event::HdPhaseDone {
                phase,
                regime: CacheRegime::Hit,
                elapsed: Duration::ZERO,
            });
        }
        match cached_c1_eroded(
            &cache_dir,
            spec.seed,
            spec.grid_size,
            &spec.init_params,
            &run,
            &spec.closures,
            &ss,
            &upscale,
        ) {
            Ok(g) => g,
            Err(e) => {
                let _ = tx.send(C1Event::HdFailed { error: format!("eroded: {e}") });
                return;
            }
        }
    } else {
        // Stateful sub-node emission (shared via Cell so we can finalise the last
        // node after the call returns).
        let cur = std::cell::Cell::new(None::<HdPhase>);
        let started = std::cell::Cell::new(Instant::now());
        let mut progress = |p: EroProgress| {
            let target = match p {
                EroProgress::Tectonic { .. } => HdPhase::Tectonic,
                EroProgress::Relief => HdPhase::Relief,
                EroProgress::Erosion { .. } => HdPhase::Erosion,
                EroProgress::Bathymetry => return, // fold: keep the Erosion node active
            };
            if cur.get() != Some(target) {
                if let Some(prev) = cur.get() {
                    let _ = tx.send(C1Event::HdPhaseDone {
                        phase: prev,
                        regime: CacheRegime::Miss,
                        elapsed: started.get().elapsed(),
                    });
                }
                let _ = tx.send(C1Event::HdPhaseStarted { phase: target });
                cur.set(Some(target));
                started.set(Instant::now());
            }
            match p {
                EroProgress::Tectonic { step, total } => {
                    let _ = tx.send(C1Event::HdPhaseProgress {
                        phase: HdPhase::Tectonic,
                        done: step,
                        total,
                    });
                }
                EroProgress::Erosion { done, total } => {
                    let _ =
                        tx.send(C1Event::HdPhaseProgress { phase: HdPhase::Erosion, done, total });
                }
                _ => {}
            }
        };
        let result = cached_c1_eroded_with_progress(
            &cache_dir,
            spec.seed,
            spec.grid_size,
            &spec.init_params,
            &run,
            &spec.closures,
            &ss,
            &upscale,
            &mut progress,
            &|| cancel.load(Ordering::Relaxed),
        );
        drop(progress);
        // Finalise the last sub-node (unless cancelled → HdFailed below).
        if result.is_ok() {
            if let Some(prev) = cur.get() {
                let _ = tx.send(C1Event::HdPhaseDone {
                    phase: prev,
                    regime: CacheRegime::Miss,
                    elapsed: started.get().elapsed(),
                });
            }
        }
        match result {
            Ok(g) => g,
            Err(_) => {
                let _ = tx.send(C1Event::HdFailed { error: "annulé".to_string() });
                return;
            }
        }
    };

    // ── Phase 2: climate (temperature + precipitation). ──
    bail_if_cancelled!();
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Climate });
    let t = Instant::now();
    let climate = c1_climate(&eroded, &ss, lat, &pp);
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Climate,
        regime: CacheRegime::Computed,
        elapsed: t.elapsed(),
    });

    // ── Phase 3: drainage (rivers + lakes + water balance). ──
    bail_if_cancelled!();
    let dkey = drainage_key(&ekey, &dcfg, &ss, Some((lat, &pp)));
    let drainage_hit = sidecar_exists("drainage", dkey.digest());
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Drainage });
    let t = Instant::now();
    let drainage =
        match cached_c1_drainage(&cache_dir, &ekey, &eroded, Some((lat, &pp)), &dcfg, &ss) {
            Ok(d) => d,
            Err(e) => {
                let _ = tx.send(C1Event::HdFailed { error: format!("drainage: {e}") });
                return;
            }
        };
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Drainage,
        regime: if drainage_hit { CacheRegime::Hit } else { CacheRegime::Miss },
        elapsed: t.elapsed(),
    });

    // ── Phase 4: biomes (Whittaker classification). ──
    bail_if_cancelled!();
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Biomes });
    let t = Instant::now();
    let biomes = c1_biomes(&eroded, &climate);
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Biomes,
        regime: CacheRegime::Computed,
        elapsed: t.elapsed(),
    });

    // ── Optional: write the v1 `.ymir` delivery container. ──
    // Explicit opt-in only (never automatic). Ships height (placeholder) +
    // coastline/cliffs (Y-B) + temperature/precipitation/biome (Y-C).
    if let Some(export_dir) = &params.export_dir {
        if let Err(e) =
            export_ymir_container(spec, &ss, &eroded, &climate, &biomes, &drainage, lat, export_dir)
        {
            // Non-fatal: the product still ships to the UI; surface the reason.
            let _ = tx.send(C1Event::HdFailed { error: format!("export .ymir: {e}") });
        }
    }

    // ── Done — ship the full product. ──
    let result = Arc::new(HdResult {
        width: eroded.width,
        height: eroded.height,
        eroded,
        temperature: climate.temperature,
        precipitation: climate.precipitation,
        drainage,
        biomes,
    });
    let _ = tx.send(C1Event::HdCompleted { result, elapsed: t_all.elapsed() });
}

/// Write a v1 `.ymir` delivery container for `eroded` under `root`
/// (`root/<name>.ymir/`, the destination configured by the caller).
///
/// Emits the manifest + the `height` raster (WP-0 placeholder, see the marker
/// below) + the `coastline` and `cliffs` vector layers (Y-B). The height dump
/// is still a normalized-range linear map, NOT true metres (that is WP-1).
#[allow(clippy::too_many_arguments)]
fn export_ymir_container(
    spec: &C1RunSpec,
    ss: &SteinSteinParams,
    eroded: &GridF32,
    climate: &ClimateResult,
    biomes: &[Biome],
    drainage: &C1DrainageResult,
    lat: f32,
    root: &Path,
) -> Result<(), String> {
    let (w, h) = (eroded.width, eroded.height);

    // WP-1: replace with c1_altitude_norm_to_metres (true metric height).
    // For now: linear-map the eroded normalized range to the full u16 span.
    let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
    for &v in &eroded.data {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    let span = (hi - lo).max(f32::EPSILON);
    let height_u16: Vec<u16> = eroded
        .data
        .iter()
        .map(|&v| (((v - lo) / span) * u16::MAX as f32).round().clamp(0.0, u16::MAX as f32) as u16)
        .collect();

    // WP (later): a real window crop. For now window_km == tectonic_domain_km.
    let meta = ContinentMeta {
        name: format!("seed{}_{}", spec.seed, w),
        seed: spec.seed,
        grid: Grid { width: w, height: h },
        window_km: C1_DOMAIN_KM as f64,
        tectonic_domain_km: C1_DOMAIN_KM as f64,
        window_offset_in_torus: [0.0, 0.0],
        latitude_deg: lat as f64,
        stein_stein: *ss,
        // WP-1: real sea-level / elevation / depth in metres.
        sea_level_m: 0.0,
        max_elevation_m: ss.asymptotic_depth_m,
        max_depth_m: ss.asymptotic_depth_m,
    };

    let dir = root.join(format!("{}.ymir", meta.name));
    let mut writer = ContinentWriter::new(&dir, meta)?;
    writer.add_raster_u16("height", &height_u16)?;

    // ── Y-B vector layers (traced from the same eroded field). ──
    // Coastline: sea-level isoline on the normalized field (sea = 0.5).
    let coastline = vector::coastline_geojson(eroded);
    writer.add_vector_file("coastline", "coastline.geojson", &coastline)?;
    writer.set_level_m("coastline", 0.0)?;

    // Cliffs: slope-threshold isoline (real angle from metric height + km/cell).
    let cell_size_m = (C1_DOMAIN_KM / w as f32) * 1000.0;
    let threshold_deg = vector::DEFAULT_CLIFF_THRESHOLD_DEG;
    let cliffs = vector::cliffs_geojson(eroded, ss, cell_size_m, threshold_deg);
    writer.add_vector_file("cliffs", "cliffs.geojson", &cliffs)?;
    writer.set_slope_threshold_deg("cliffs", threshold_deg as f64)?;

    // ── Y-C climate layers (read c1_climate / c1_biomes outputs; no re-compute). ──
    let n = w * h;
    let temperature: Vec<i16> = climate
        .temperature
        .data
        .iter()
        .map(|&t| (t * 100.0).round().clamp(i16::MIN as f32, i16::MAX as f32) as i16)
        .collect();
    let precipitation: Vec<u16> = climate
        .precipitation
        .data
        .iter()
        .map(|&p| precip_mm_per_year(p).round().clamp(0.0, u16::MAX as f32) as u16)
        .collect();
    let biome: Vec<u8> = biomes.iter().map(|b| b.to_u8()).collect();
    debug_assert_eq!(biome.len(), n);
    writer.add_raster_i16("temperature", &temperature)?;
    writer.add_raster_u16("precipitation", &precipitation)?;
    writer.add_raster_u8("biome", &biome)?;

    // ── Hydro layers (serialize existing drainage outputs; no re-compute). ──
    writer.add_raster_u32("lake_mask", &drainage.lake_map)?;
    writer.add_raster_f32("flow_accumulation", &drainage.flow.accumulation.data)?;
    let rivers = hydro::rivers_json(drainage);
    writer.add_vector_file("rivers", "rivers.json", &rivers)?;
    let lakes = hydro::lakes_json(drainage);
    writer.add_vector_file("lakes", "lakes.json", &lakes)?;

    // water_class: reuse the SAME sea level the coastline trace uses (no second
    // constant) — 0 land, 1 ocean (edge-connected), 2 inland (enclosed below-sea).
    let water = connectivity::water_class(eroded, vector::SEA_LEVEL_NORM);
    writer.add_raster_u8("water_class", &water)?;

    writer.finish()?;
    Ok(())
}
