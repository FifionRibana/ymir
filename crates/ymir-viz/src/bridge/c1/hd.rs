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
use ymir_core::climate::{ClimateResult, c1_biomes_classified_wet, c1_climate_placed, c1_climate_windowed};
use ymir_core::erosion::stream_power::StreamPowerConfig;
use ymir_core::export::container::{ContinentMeta, ContinentWriter, Grid};
use ymir_core::export::{height, hydro, vector};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::{
    c1_coarse_land_report, cached_c1_drainage_windowed, cached_c1_eroded,
    cached_c1_eroded_with_progress, coarse_normalized_sweep, drainage_key_windowed, eroded_key,
    tectonic_key,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, LakeType, apply_geo_scale_ratio,
    below_sea_basin_lakes, c1_drainage_windowed, clip_rivers_to_lakes,
};
use ymir_core::tectonics_c1::land_topology::{
    IslandEval, LandTopology, evaluate_island, land_topology,
};
use ymir_core::tectonics_c1::production_upscale::EroProgress;
use ymir_core::tectonics_c1::time_loop::C1TimeLoopConfig;
use ymir_core::terrain::flow::{RiverSegment, breach_monotone};
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
    /// Physical size (km) the tectonic domain represents — the domain IS the map
    /// (M1 #190, no crop). The whole torus renders at `target_size²`, so
    /// `m_per_px = domain_km · 1000 / target_size` (1024 km @ 2048² → 500 m/px).
    /// A pure relabelling of the 64² tectonic pattern: it stamps the manifest
    /// scale and drives climate/drainage km/cell, but NOT the tectonic sim. The
    /// tectonic cache key does not include it. Default 1024.0.
    pub domain_km: f32,
    /// Framing ROLL to apply at the coarse sampling origin (normalised `[0,1)`).
    /// `None` → `run_hd` centres the largest mass automatically (the default). When
    /// the UI has computed a framing (auto ± manual pan) it passes `Some(offset)`
    /// so the export matches exactly what the preview showed. Determinism: passing
    /// the same value `run_hd` would auto-compute yields a byte-identical export.
    pub manual_offset: Option<[f64; 2]>,
    /// EXPERIMENTAL (ADR 0001): replace droplet erosion with routed stream-power
    /// incision + hillslope regime split for this run. `false` (default) → the
    /// production droplet pass (byte-identical). `true` → droplets OFF, stream-power
    /// ON with the recommended config (K=3, m=0.5, n=1, iters=3, A_c=7.6 km²,
    /// D=0.05, uncoupled). A UI opt-in to eyeball the effect; production default
    /// unchanged until confirmed.
    pub stream_power: bool,
    /// EXPERIMENTAL (ADR 0001, Finding 7): with `stream_power` on, use the `relief_v2`
    /// config — the two bounding CLOSURES (nonlinear hillslope diffusion with a
    /// critical slope → arêtes; hydraulic-geometry lateral widening → valley floors).
    /// `false` → `relief_v1` (v1 slits). A UI opt-in; production default unchanged.
    pub closures: bool,
    /// EXPERIMENTAL (ADR 0001, Finding 9): with `closures` on, apply hillslope diffusion
    /// EVERYWHERE (`diffuse_channels`) at a super-critical D to damp the Smith–Bretherton
    /// parallel-rilling comb. `false` → the regime split (combed at 8192²).
    pub cross_rill: bool,
    /// Cross-rill diffusion coefficient when `cross_rill` is on (2b sweep: 0.25/0.40/0.55).
    pub cross_rill_d: f32,
    /// EXPERIMENTAL (ADR 0001, Finding 11): MFD incision — disperse drainage area to prevent
    /// the comb at its cause (dendritic valleys, no GS solver). Takes precedence over cross_rill.
    pub mfd: bool,
    /// MFD partition exponent when `mfd` is on (lower = more dispersed; p≈2 recommended).
    pub mfd_p: f32,
    /// EXPERIMENTAL: override the FBM `amplitude_base` for this run (`None` = the
    /// production 0.16). Lets the author flip through the striation amplitude ladder
    /// (0.16/0.08/0.04/0.02) with stream-power ON to see the striations shrink.
    pub fbm_amplitude: Option<f64>,
    /// ADR 0001 Finding 24 — GEOGRAPHIC SCALE RATIO (hydrology only). The map DRAWS
    /// `domain_km` but SIGNIFIES `domain_km · ratio` for the EXPORT-DERIVED hydrology
    /// (catchment ×ratio², discharge ×ratio², channel width ×ratio, navigability
    /// re-classed). A pure presentation multiplier — NOTHING that shapes the terrain
    /// (incision, lake balance, climate, biomes) sees it. `1.0` (default) = identity.
    pub geo_scale_ratio: f32,
    /// ADR 0001 Finding 25 — CLIMATIC latitude SPAN in degrees, centred on
    /// `latitude_deg`. Decouples the climate extent from the physical extent so a
    /// small island can span several belts (tundra↔desert). `None` = the geographic
    /// span `domain_km / 111` (identity). REAL physics: it drives temperature, wind
    /// belts, precipitation and hence biomes (recomputed, not cached like the ratio).
    pub latitude_span_deg: Option<f32>,
    /// When `Some(dir)`, after the biome phase the run writes a v1 `.ymir`
    /// delivery container under `dir` (see [`ymir_core::export::container`]).
    /// `None` = no export. Explicit opt-in — the pipeline NEVER auto-exports
    /// (WP-0). The container directory is `dir/<name>.ymir/`.
    pub export_dir: Option<PathBuf>,
}

impl Default for HdParams {
    fn default() -> Self {
        Self {
            target_size: 2048,
            latitude_deg: 45.0,
            domain_km: 1024.0,
            manual_offset: None,
            stream_power: false,
            closures: false,
            cross_rill: false,
            cross_rill_d: 0.40,
            mfd: false,
            mfd_p: 2.0,
            fbm_amplitude: None,
            geo_scale_ratio: 1.0,
            latitude_span_deg: None,
            export_dir: None,
        }
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
    /// M1 land-topology of the full coarse torus (island-continent judgement:
    /// number of masses, largest area, wrap flags). Independent of the window.
    pub land_topology: LandTopology,
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

/// Fast tectonic-shape preview: the calibrated COARSE continent (one ~1-2 s
/// tectonic pass, NO HD upscale/erosion) so a seed can be judged BEFORE the
/// ~20-min HD run. Carries the coarse normalized field (for a heightmap render),
/// the island verdict ([`IslandEval`] — border-clean? area? wraps?), and the
/// export-window rectangle in normalized torus coords (matching `run_hd`).
pub struct PreviewShape {
    /// Coarse normalized altitude (sea = 0.5), calibrated at the HD tlf, in the
    /// UNROLLED torus frame. The UI applies the framing roll (auto ± manual pan) at
    /// render time, so panning is pure relabelling — no tectonic recompute — and
    /// the domain-as-map metrics recompute instantly for any offset.
    pub coarse: GridF32,
    /// Auto-computed framing offset (integer COARSE CELLS) that centres the largest
    /// landmass: `(center_cell − grid/2) mod grid`. The UI seeds the current offset
    /// from this on every new preview and offers a "snap to auto" back to it.
    pub auto_offset_cells: [i64; 2],
    /// Island verdict on the coarse field (largest-mass summary for the Debug log).
    pub eval: IslandEval,
}

impl fmt::Debug for PreviewShape {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PreviewShape")
            .field("grid", &self.coarse.width)
            .field("largest_km2", &self.eval.topo.largest_area_km2)
            .field("border_clean", &self.eval.border_clean)
            .finish()
    }
}

/// Run only the coarse tectonic pass + calibration for `spec`, emitting a
/// [`C1Event::PreviewReady`]. Cheap (no erosion); used by the "Aperçu" button.
pub fn preview_shape(spec: &C1RunSpec, params: &HdParams, tx: &Sender<C1Event>) {
    let t = Instant::now();
    let run = hd_run_config(spec);
    let ss = SteinSteinParams::default();
    // Same sea-level calibration the HD run uses (read it off the HD config so
    // the preview coastline matches what will be exported).
    let tlf = FbmUpscaleConfig::c1_hd_production(params.target_size).target_land_fraction;
    let coarse = coarse_normalized_sweep(
        spec.seed,
        spec.grid_size,
        &spec.init_params,
        &run,
        &spec.closures,
        &ss,
        &[tlf.unwrap_or(0.29)],
    )
    .remove(0)
    .1;
    let eval = evaluate_island(&coarse, vector::SEA_LEVEL_NORM, 25.0, 1);
    // The domain IS the map (M1 #190): no crop. Return the UNROLLED coarse field +
    // the auto framing offset (centre the largest mass). The UI applies the roll
    // (auto ± manual pan) at render time, so panning is instant (no recompute) and
    // what the user frames is exactly what `run_hd` exports.
    let topo = land_topology(&coarse, vector::SEA_LEVEL_NORM);
    let g = spec.grid_size as i64;
    let (cx, cy) = topo.center_cell;
    let auto_offset_cells =
        [((cx as i64) - g / 2).rem_euclid(g), ((cy as i64) - g / 2).rem_euclid(g)];
    let preview = Arc::new(PreviewShape { coarse, auto_offset_cells, eval });
    let _ = tx.send(C1Event::PreviewReady { preview, elapsed: t.elapsed() });
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

    // ── The domain IS the map (M1 #190): NO crop, ever. The whole torus renders at
    // `target_size` (sample_origin [0,0], sample_size 1.0). m_per_px =
    // domain_km·1000/target_size, so the export size is chosen for the cell budget
    // (e.g. 24576² → ~42 m/px at 1024 km). Seed selection — not a window — puts a
    // continent with ocean margin on the map (see the coarse preview / seed scan).
    let window_km = params.domain_km; // domain IS the map: window == domain == domain_km
    let mut upscale = FbmUpscaleConfig::c1_hd_production(params.target_size);
    if let Some(a) = params.fbm_amplitude {
        upscale.amplitude_base = a; // striation-ladder override (default 0.16)
    }
    // EXPERIMENTAL opt-in (ADR 0001): swap the droplet pass for routed stream-power
    // incision + hillslope regime split. A_c is 7.6 km² expressed in cells at THIS
    // resolution (resolution-stable channel head); uncoupled vertical scale.
    if params.stream_power {
        let km_per_cell = window_km / params.target_size as f32; // window_km == domain_km
        let cell_km2 = km_per_cell * km_per_cell;
        let depth = ss.depth_scale_m as f32;
        let mut sp = if params.closures {
            StreamPowerConfig::relief_v2(cell_km2, depth)
        } else {
            StreamPowerConfig::relief_v1(cell_km2, depth)
        };
        // STEP 2b (ADR 0001 Finding 9): cross-rill diffusion — diffuse EVERYWHERE (drop
        // the regime split's channel exclusion) at a super-critical D, to damp the
        // Smith–Bretherton parallel-rilling comb. Experimental, only atop closures.
        if params.closures && params.cross_rill {
            sp.diffuse_channels = true;
            sp.diffusion = params.cross_rill_d;
        }
        // ADR 0001 Findings 11/14: the FINAL relief recipe (relief-v3) — MFD incision
        // (comb prevented at the cause, dendritic valleys) + talus flank grading + light
        // linear hillslope, K×3. No GS solver. The eroded field is then breach-conditioned
        // (below) for monotone rivers + lakes. Takes precedence over cross-rill.
        if params.closures && params.mfd {
            sp = StreamPowerConfig::relief_v3(cell_km2, depth);
            sp.mfd_exponent = Some(params.mfd_p); // let the viz p selector override
        }
        eprintln!(
            "[HD] stream-power incision ON ({}: A_c {:.0} cells = {} km²{}), droplets OFF",
            if params.closures { "relief-v2 + closures" } else { "relief-v1" },
            sp.min_area_cells,
            ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2,
            if params.closures {
                format!(
                    ", S_c=tan(33°), lat {:.0} m/√km²{}",
                    sp.lateral_erosion,
                    if params.mfd {
                        format!(", MFD p={:.1} K×3", params.mfd_p)
                    } else if params.cross_rill {
                        format!(", cross-rill D={:.2}", sp.diffusion)
                    } else {
                        String::new()
                    },
                )
            } else {
                String::new()
            },
        );
        upscale.erosion = None; // droplets off — they collapse the SP valleys
        upscale.stream_power = Some(sp);
    }
    // Land report FIRST (reads only `target_land_fraction`) so the seam-correct
    // torus centre of the largest mass is known before we choose the sampling roll.
    let land = c1_coarse_land_report(
        spec.seed,
        spec.grid_size,
        &spec.init_params,
        &run,
        &spec.closures,
        &ss,
        upscale.target_land_fraction,
    );
    let _ = land.centroid; // superseded by the seam-correct centre below
    let land_topology = land.topology;

    // ── Framing ROLL (M1 #190): shift the coarse SAMPLING ORIGIN so the largest
    // landmass is contiguous and centred in the exported (whole) domain. A mass
    // straddling the torus seam would otherwise export split across the border
    // (Living Landz reads it as two landmasses). This is a ROLL, not a crop:
    // sample_size stays 1.0, the whole torus is exported, no cell is discarded.
    //
    // Offset = integer coarse cells (exact, cache-stable): roll the mass centre to
    // the domain centre — origin = (center_cell − grid/2) mod grid. Because the FBM
    // + coast warp are evaluated in TORUS coords (sx = origin + i·scale), the noise
    // rolls WITH the terrain, so the exported continent is the one the preview shows.
    let g = spec.grid_size as i64;
    let (cx, cy) = land_topology.center_cell;
    let roll_x = ((cx as i64) - g / 2).rem_euclid(g) as usize;
    let roll_y = ((cy as i64) - g / 2).rem_euclid(g) as usize;
    let auto_offset = [roll_x as f64 / g as f64, roll_y as f64 / g as f64];
    // Use the UI framing (auto ± manual pan) when supplied; else auto-centre.
    let window_offset = params.manual_offset.unwrap_or(auto_offset);
    upscale.sample_origin = window_offset;
    upscale.sample_size = 1.0;

    // Telemetry — judge the seed as an island continent. A largest landmass that
    // WRAPS the torus is a band with no coast; a finite mass touching the map edge
    // is split across the border (see the domain-as-map verdict in the preview).
    let t = &land_topology;
    eprintln!(
        "[C1 land] seed {} — {} landmass(es); largest {:.0} km² ({:.1}% of torus), \
         wraps x={} y={}, bbox {:.0}×{:.0} km; coarse emerged {:.1}%; roll ({},{}) cells → offset [{:.3},{:.3}]",
        spec.seed,
        t.num_landmasses,
        t.largest_area_km2,
        t.largest_area_frac * 100.0,
        t.wraps_x,
        t.wraps_y,
        t.bbox_km.0,
        t.bbox_km.1,
        t.emerged_fraction * 100.0,
        roll_x,
        roll_y,
        window_offset[0],
        window_offset[1],
    );

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

    // ── relief-v3 conditioning (ADR 0001 Finding 14): breach the eroded field so exported
    // river long-profiles are monotone, with detected lakes held flat. Done AFTER the cached
    // eroded (in-memory, not re-cached): detect lakes on the pre-breach field, breach
    // (lakes excepted), then compute drainage on the conditioned field for MONOTONE rivers
    // while carrying the PRE-BREACH lakes (the breach flattens them, so post-breach detection
    // would find none — lakes.json must come from the pre-breach set). Uncached (experimental
    // viz path); the cached drainage below is bypassed when this override is present.
    // Returns the conditioned field + the PRE-BREACH drainage (its lakes/lake_map are
    // the real water bodies; the breach flattens them so post-breach detection finds
    // none). The FINAL, discharge-bearing drainage is computed below WITH the climate.
    let (eroded, prebreach_override) = if params.stream_power && params.closures && params.mfd {
        let (gw, gh) = (eroded.width, eroded.height);
        let prebreach = c1_drainage_windowed(&eroded, None, &dcfg, &ss, window_km);
        let conditioned =
            breach_monotone(&eroded, &prebreach.flow.filled, &prebreach.lake_map, 0.5, gw, gh);
        eprintln!(
            "[HD] relief-v3 breach conditioning: monotone rivers + {} lakes held flat",
            prebreach.lakes.len()
        );
        (conditioned, Some(prebreach))
    } else {
        (eroded, None)
    };

    // ── Phase 2: climate (temperature + precipitation). ──
    bail_if_cancelled!();
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Climate });
    let t = Instant::now();
    // Finding 25 — explicit CLIMATIC latitude span (else the geographic domain_km/111).
    // `None` keeps the byte-identical windowed path; `Some` re-derives temperature + the
    // per-belt precipitation over the wider span (real physics → biomes shift).
    let climate = match params.latitude_span_deg {
        Some(span) => c1_climate_placed(&eroded, &ss, lat, span, &pp, window_km),
        None => c1_climate_windowed(&eroded, &ss, lat, &pp, window_km),
    };
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Climate,
        regime: CacheRegime::Computed,
        elapsed: t.elapsed(),
    });

    // ── Phase 3: drainage (rivers + lakes + water balance). ──
    bail_if_cancelled!();
    let dkey = drainage_key_windowed(&ekey, &dcfg, &ss, Some((lat, &pp)), window_km);
    let drainage_hit = sidecar_exists("drainage", dkey.digest());
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Drainage });
    let t = Instant::now();
    let mut drainage = if let Some(prebreach) = prebreach_override {
        // relief-v3: final drainage on the breached field WITH the climate, so per-segment
        // DISCHARGE (→ channel width, Finding 22) is real. The breached field has no pits,
        // so its own lake detection is empty — carry the pre-breach lakes/lake_map instead.
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let mut dr = c1_drainage_windowed(&eroded, Some(&dclim), &dcfg, &ss, window_km);
        dr.lakes = prebreach.lakes;
        dr.lake_map = prebreach.lake_map;
        dr
    } else {
        match cached_c1_drainage_windowed(
            &cache_dir,
            &ekey,
            &eroded,
            Some((lat, &pp)),
            &dcfg,
            &ss,
            window_km,
        ) {
            Ok(d) => d,
            Err(e) => {
                let _ = tx.send(C1Event::HdFailed { error: format!("drainage: {e}") });
                return;
            }
        }
    };
    // ADR 0001 Finding 18: below-sea ENCLOSED basins as typed water bodies. pit_fill/
    // detect_lakes treat them as ocean, so they never enter the lake path; add them here via
    // the water balance (level = min(spill, evaporative)), merging their water cells into the
    // lake_map (only where empty) and their typed lakes into the list → lakes.json + biomes.
    let wetland_mask: Vec<u8> = {
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let bs = below_sea_basin_lakes(&eroded, &dclim, &dcfg, &ss, window_km);
        let (mut endo, mut exo) = (0usize, 0usize);
        for lk in &bs.lakes {
            match lk.lake_type {
                LakeType::Endorheic => endo += 1,
                LakeType::Exorheic => exo += 1,
            }
        }
        for k in 0..bs.lake_map.len() {
            if bs.lake_map[k] != 0 && drainage.lake_map[k] == 0 {
                drainage.lake_map[k] = bs.lake_map[k];
            }
        }
        // Finding 30 — append each exorheic basin's traced SPILLWAY as a watercourse (so an
        // exorheic label always has an outlet; may chain through basins before the sea).
        let (mut to_sea, mut chained) = (0usize, 0usize);
        for sw in &bs.spillways {
            if sw.chained_into.is_some() {
                chained += 1;
            } else {
                to_sea += 1;
            }
            drainage.rivers.segments.push(RiverSegment {
                points: sw.points.clone(),
                strahler_order: 1,
                avg_flow: 0.0,
                max_flow: 0.0,
                basin_id: 0,
                upstream: vec![],
                downstream: None,
            });
            drainage.segment_drainage_km2.push(sw.drainage_km2);
            drainage.segment_navigability.push(sw.navigability);
            drainage.segment_discharge_m3s.push(sw.discharge_m3s);
            drainage.segment_width_m.push(sw.width_m);
            drainage.segment_profile_m.push(sw.profile_m.clone());
        }
        drainage.lakes.extend(bs.lakes);
        eprintln!(
            "[HD] below-sea basins: {} lakes ({exo} exorheic, {endo} endorheic); {} spillways ({to_sea} → sea, {chained} chained)",
            exo + endo,
            bs.spillways.len()
        );
        bs.wetland
    };
    // Finding 24 — geographic scale ratio: a PURE post-process on the export-derived
    // hydrology (catchment/discharge/width/navigability), applied before the clip so the
    // clipped reaches inherit the signified values. Identity at 1.0. Touches nothing
    // physical (the terrain, climate and biomes above ran on the real domain_km).
    apply_geo_scale_ratio(&mut drainage, params.geo_scale_ratio, &dcfg.thresholds);
    if params.geo_scale_ratio != 1.0 {
        eprintln!("[HD] geographic scale ratio {:.2} → hydrology signifies ×{:.1} area", params.geo_scale_ratio, params.geo_scale_ratio * params.geo_scale_ratio);
    }
    // ADR Finding 20 (DEFECT A): the rivers were traced on the breached (monotone)
    // field, so they run straight through lakes and out of endorheic sinks while the
    // final lake_map marks those basins as water. Clip the exported network to that
    // lake_map so every watercourse terminates at its sink (sea or lake); an endorheic
    // basin's phantom outlet is dropped. Topology stays continuous — width is a hint.
    let before_seg = drainage.rivers.segments.len();
    clip_rivers_to_lakes(&mut drainage);
    eprintln!(
        "[HD] rivers clipped to lakes: {before_seg} → {} segments (terminate at sinks)",
        drainage.rivers.segments.len()
    );
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Drainage,
        regime: if drainage_hit { CacheRegime::Hit } else { CacheRegime::Miss },
        elapsed: t.elapsed(),
    });

    // ── Phase 4: biomes (Whittaker classification). ──
    bail_if_cancelled!();
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Biomes });
    let t = Instant::now();
    // ADR Finding 18: biomes from water_class + the drainage lake_map (below-sea enclosed
    // basins → Lake / exposed land, not Ocean), instead of the altitude threshold.
    let biomes = c1_biomes_classified_wet(&eroded, &climate, &drainage.lake_map, &wetland_mask);
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Biomes,
        regime: CacheRegime::Computed,
        elapsed: t.elapsed(),
    });

    // ── Optional: write the v1 `.ymir` delivery container. ──
    // Explicit opt-in only (never automatic). Ships height (placeholder) +
    // coastline/cliffs (Y-B) + temperature/precipitation/biome (Y-C).
    if let Some(export_dir) = &params.export_dir {
        // Resolved climatic span (explicit, else the geographic domain_km/111) so the
        // manifest records the actual gradient behind the climate/biome layers.
        let lat_span = params.latitude_span_deg.unwrap_or(window_km / 111.0);
        if let Err(e) = export_ymir_container(
            spec,
            &ss,
            &eroded,
            &climate,
            &biomes,
            &drainage,
            lat,
            lat_span,
            params.geo_scale_ratio,
            window_km,
            window_offset,
            export_dir,
        ) {
            // Non-fatal: the product still ships to the UI; surface the reason.
            let _ = tx.send(C1Event::HdFailed { error: format!("export .ymir: {e}") });
        }
    }

    // Telemetry: emerged fraction measured AFTER FBM+erosion on the window grid.
    // It differs from the coarse target (resolution-dependent until the FBM band
    // policy is fixed — separate issue); reported, not compensated.
    let post_emerged =
        eroded.data.iter().filter(|&&v| v > 0.5).count() as f32 / eroded.data.len().max(1) as f32;
    eprintln!(
        "[C1 land] post-FBM+erosion emerged fraction (window) = {:.1}%",
        post_emerged * 100.0
    );

    // ── Done — ship the full product. ──
    let result = Arc::new(HdResult {
        width: eroded.width,
        height: eroded.height,
        eroded,
        temperature: climate.temperature,
        precipitation: climate.precipitation,
        drainage,
        biomes,
        land_topology,
    });
    let _ = tx.send(C1Event::HdCompleted { result, elapsed: t_all.elapsed() });
}

/// Write a v1 `.ymir` delivery container for `eroded` under `root`
/// (`root/<name>.ymir/`, the destination configured by the caller).
///
/// Emits the manifest + the `height` raster (true metres via the vertical
/// contract, quantised to u16 over the field's real range — see the body) +
/// the `coastline` and `cliffs` vector layers (Y-B).
#[allow(clippy::too_many_arguments)]
fn export_ymir_container(
    spec: &C1RunSpec,
    ss: &SteinSteinParams,
    eroded: &GridF32,
    climate: &ClimateResult,
    biomes: &[Biome],
    drainage: &C1DrainageResult,
    lat: f32,
    lat_span_deg: f32,
    geo_scale_ratio: f32,
    window_km: f32,
    window_offset: [f64; 2],
    root: &Path,
) -> Result<(), String> {
    let (w, h) = (eroded.width, eroded.height);

    // Metric height (the WP-1 vertical contract). Convert the normalized field
    // to TRUE metres via the single vertical contract, anchored on the SAME
    // sea-level constant the coastline is traced at (so 0 m == the coastline),
    // then quantise linearly to u16 over the field's real [min_m, max_m]:
    //   code = round((m − min_m) / (max_m − min_m) · 65535)  (decode is the
    //   inverse). The min_m/max_m are stamped on the "height" layer below.
    let height = height::metric_height_u16(eroded, ss);

    // The domain IS the map (M1 #190): the HD grid renders the WHOLE torus, so
    // window_km == tectonic_domain_km (no crop). window_offset_in_torus carries the
    // framing ROLL (normalised 0..1) applied at the coarse sampling origin so the
    // continent is contiguous/centred — Living Landz reads it to know where the
    // frame sits on the torus. km_per_cell = window_km / width = domain_km / width.
    let meta = ContinentMeta {
        name: format!("seed{}_{}", spec.seed, w),
        seed: spec.seed,
        grid: Grid { width: w, height: h },
        window_km: window_km as f64,
        tectonic_domain_km: window_km as f64, // window == domain (no crop)
        window_offset_in_torus: window_offset,
        latitude_deg: lat as f64,
        latitude_span_deg: lat_span_deg as f64,
        geographic_scale_ratio: geo_scale_ratio as f64,
        stein_stein: *ss,
        // Honest metric bounds: sea anchored to 0 m; elevation/depth are the
        // field's real metric extrema (depth is the min, i.e. most negative).
        sea_level_m: 0.0,
        max_elevation_m: height.max_m as f64,
        max_depth_m: height.min_m as f64,
    };

    let dir = root.join(format!("{}.ymir", meta.name));
    let mut writer = ContinentWriter::new(&dir, meta)?;
    writer.add_raster_u16("height", &height.codes)?;
    writer.set_metric_range("height", height.min_m as f64, height.max_m as f64)?;

    // ── Y-B vector layers (traced from the same eroded field). ──
    // Coastline: sea-level isoline on the normalized field (sea = 0.5).
    let coastline = vector::coastline_geojson(eroded);
    writer.add_vector_file("coastline", "coastline.geojson", &coastline)?;
    writer.set_level_m("coastline", 0.0)?;

    // Cliffs: slope-threshold isoline (real angle from metric height + km/cell).
    // Window km/cell — the same scale the climate/drainage use.
    let cell_size_m = (window_km / w as f32) * 1000.0;
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
