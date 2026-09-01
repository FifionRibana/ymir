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

use ymir_core::cache::{cached, default_cache_dir};
use ymir_core::climate::biomes::Biome;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::climate::{ClimateResult, c1_biomes_classified_wet};
use ymir_core::erosion::stream_power::StreamPowerConfig;
use ymir_core::export::container::{ContinentMeta, ContinentWriter, Grid};
use ymir_core::export::{height, hydro, vector};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::{
    HdDrainageBundle, c1_coarse_land_report, cached_c1_climate, cached_c1_drainage_windowed,
    cached_c1_eroded, cached_c1_eroded_with_progress, climate_key, coarse_normalized_sweep,
    conditioned_eroded_key, drainage_key_windowed, eroded_key_full, hd_drainage_key, tectonic_key,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, LakeType, apply_geo_scale_ratio,
    apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed_infil,
    clip_rivers_to_lakes, exorheic_lakes_missing_outlet,
};
use ymir_core::tectonics_c1::land_topology::{
    IslandEval, LandTopology, evaluate_island, land_topology,
};
use ymir_core::tectonics_c1::production_upscale::EroProgress;
use ymir_core::tectonics_c1::time_loop::C1TimeLoopConfig;
use ymir_core::terrain::flow::{RiverSegment, breach_monotone_protected};
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
    /// C-2 volcanism (closures roadmap §2). `None` (default) → OFF, byte-identical
    /// to the pre-C-2 terrain. `Some(cfg)` → derive edifices from the tectonic
    /// state and inject them at HD after the FBM. `domain_km` is overwritten with
    /// `params.domain_km` (the geometric span) at use, so the caller need only set
    /// `enabled: true` and any tuning.
    pub volcanism: Option<ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig>,
    /// C-3 lithological heterogeneity (closures roadmap §3). `None` (default) → OFF,
    /// byte-identical (uniform erodibility). `Some(cfg)` with `enabled` → build the
    /// per-cell erodibility multiplier (hard basement ×1.0, rift-soft + volcaniclastic
    /// ABOVE, causal — never noise/geometry) and thread it into the stream-power
    /// incision. Only takes effect with `stream_power` on (that is where K threads in).
    pub lithology: Option<ymir_core::tectonics_c1::closures::lithology::LithologyConfig>,
    /// C-3b inherited structure (closures roadmap §3b). `None` (default) → OFF,
    /// byte-identical. `Some(cfg)` with `enabled` → modulate erodibility ISOTROPICALLY
    /// by fracture density (proximity to plate contacts): intact craton retains relief
    /// (×1 reference), fractured belts erode more. Only bites with `stream_power` on.
    pub fracture: Option<ymir_core::tectonics_c1::closures::fracture::FractureConfig>,
    /// H-1 infiltration (the first subsurface term). `None` (default) → OFF, the pre-H-1
    /// water balance byte-identical. `Some(cfg)` with `enabled` → a causal permeability
    /// field (lithology matrix + fracture density, double porosity) sets the fraction of
    /// the precipitation surplus that infiltrates and never reaches a lake by surface flow.
    pub infiltration: Option<ymir_core::tectonics_c1::closures::infiltration::InfiltrationConfig>,
    /// Debug microscope: derive the coarse tectonic labels (rift/subduction/craton/…)
    /// for the overlay. Costs one extra ~1 s coarse pass; `false` (default) skips it.
    pub emit_tectonic_labels: bool,
}

impl Default for HdParams {
    fn default() -> Self {
        Self {
            target_size: 2048,
            latitude_deg: 45.0,
            domain_km: 1024.0,
            manual_offset: None,
            // relief-v3 is the production default (#190) — stream-power + closures + MFD ON.
            stream_power: true,
            closures: true,
            cross_rill: false,
            cross_rill_d: 0.40,
            mfd: true,
            mfd_p: 2.0,
            fbm_amplitude: None,
            geo_scale_ratio: 1.0,
            latitude_span_deg: None,
            export_dir: None,
            volcanism: None,
            lithology: None,
            fracture: None,
            infiltration: None,
            emit_tectonic_labels: false,
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
    /// C-2 volcanic edifices placed in this render (crater centre in DATA-space
    /// pixels, radius, activity, kind, setting). Empty when volcanism is off. The
    /// UI lists these in the microscope and marks them on the map.
    pub volcanoes: Vec<ymir_core::tectonics_c1::closures::volcanism::CraterRecord>,
    /// Debug microscope (C-2/C-3): coarse tectonic labels (rift, subduction,
    /// collision, craton, lithology class) for the overlay. `None` when not requested
    /// (`HdParams::emit_tectonic_labels`). Nearest-sampled to the HD window via
    /// `sample_origin`/`sample_size`, so every overlay registers with the terrain.
    pub tectonic: Option<ymir_core::tectonics_c1::debug_labels::CoarseTectonicLabels>,
    /// Volcanic edifices (torus-UV centre + basal diameter km) for the lithology
    /// overlay's volcaniclastic footprints. Empty when volcanism is off.
    pub edifices: Vec<ymir_core::tectonics_c1::closures::volcanism::Edifice>,
    /// The window sampling the overlays register against (same as the terrain upscale).
    pub sample_origin: [f64; 2],
    pub sample_size: f64,
    /// Physical km per HD cell (`sample_size · domain_km / width`) — for the basal-disc
    /// radius of the volcaniclastic overlay.
    pub km_per_cell: f32,
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
/// #190 — build the FINAL HD drainage + wetland mask (the climate-dependent tail of the pipeline):
/// base drainage with real discharge, the below-sea water-balance merge + cleanup (Findings 18/39),
/// the river clip, the spillway append and the geographic-scale post-process. A PURE function of
/// exactly the inputs `hd_drainage_key` folds — so caching it under that key is stale-proof.
#[allow(clippy::too_many_arguments)]
fn build_hd_drainage(
    eroded: &GridF32,
    climate: &ClimateResult,
    prebreach: Option<C1DrainageResult>,
    dcfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    window_km: f32,
    geo_scale_ratio: f32,
    // H-1: per-cell infiltrated fraction (None → pre-H-1 balance, byte-identical).
    infiltration: Option<&[f32]>,
) -> HdDrainageBundle {
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    // relief-v3: final drainage on the breached field, carrying the pre-breach lakes; legacy path
    // (prebreach None): plain drainage on the raw eroded. Both with the climate discharge (Finding 22).
    let (mut drainage, carried_geometric_lakes) = if let Some(prebreach) = prebreach {
        let mut dr =
            c1_drainage_windowed_infil(eroded, Some(&dclim), dcfg, ss, window_km, infiltration);
        // The pre-breach lakes are the GEOMETRICALLY correct ones (the breach destroys the
        // depressions) — but they were detected with `climate = None`, so they carry NO
        // water balance. Their geometry is adopted; their classification is fixed below.
        dr.lakes = prebreach.lakes;
        dr.lake_map = prebreach.lake_map;
        (dr, true)
    } else {
        (c1_drainage_windowed_infil(eroded, Some(&dclim), dcfg, ss, window_km, infiltration), false)
    };
    // H-1 — THE MISSING LINK, and a CORRECTION, so it runs BY DEFAULT (not gated): lakes
    // classified without ever seeing the climate are simply wrong. The relief-v3 path
    // carries the climate-free pre-breach lakes, so they were exorheic by pure geometry
    // ("the outlet reaches the sea"). Run the balance on them and adopt the CLASSIFICATION
    // ONLY — geometry, levels and footprints stay untouched (adopting an endorheic
    // equilibrium LEVEL would be H-2 by the back door). Crater types (C-2) are never
    // overwritten. The optional H-1 INFILTRATION rides in `infiltration` (measured a
    // secondary term: 0–3 lakes, against 53–100 % from this reclassification alone).
    if carried_geometric_lakes {
        let (w, h) = (eroded.width, eroded.height);
        let cell_km2 = (window_km / w as f32).powi(2);
        let before_n = drainage.lakes.len();
        let before_km2: f32 = drainage.lakes.iter().map(|l| l.area_km2).sum();
        // H-1c — APPLY the balance: classify AND settle endorheic basins at their
        // evaporative equilibrium (level + footprint), draining the exposed floor from
        // `lake_map`. Runs BEFORE `below_sea_basin_lakes_infil` and BEFORE
        // `clip_rivers_to_lakes` so both see the FINAL footprint — river tracks are clipped
        // to the retreated outline instead of ending in the void (the orphaned-mouth defect
        // already fixed twice: enumerate inlets AFTER the footprint is known).
        let lakes_in = std::mem::take(&mut drainage.lakes);
        drainage.lakes = apply_lake_water_balance(
            eroded,
            &drainage.flow,
            &dclim,
            cell_km2,
            ss,
            &lakes_in,
            &mut drainage.lake_map,
            infiltration,
            w,
            h,
        );
        let endo = drainage.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).count();
        let after_km2: f32 = drainage.lakes.iter().map(|l| l.area_km2).sum();
        eprintln!(
            "[HD] H-1c surface water balance APPLIED: {} → {} lakes | {:.0} → {:.0} km² water ({:.0} km² floor exposed) | {} endorheic",
            before_n,
            drainage.lakes.len(),
            before_km2,
            after_km2,
            (before_km2 - after_km2).max(0.0),
            endo
        );
    }
    let (wetland_mask, below_sea_spillways) = {
        let bs = below_sea_basin_lakes_infil(
            eroded,
            &dclim,
            dcfg,
            ss,
            window_km,
            Some(&drainage.lake_map),
            infiltration,
        );
        let (mut endo, mut exo) = (0usize, 0usize);
        for lk in &bs.lakes {
            match lk.lake_type {
                LakeType::Endorheic => endo += 1,
                LakeType::Exorheic => exo += 1,
                // Crater types are assigned later (post-drainage), never here.
                LakeType::CraterAcidic | LakeType::CraterNeutral => {}
            }
        }
        use std::collections::{HashMap, HashSet};
        let mut det_before: HashMap<u32, usize> = HashMap::new();
        for &id in &drainage.lake_map {
            if id != 0 && id < 1_000_001 {
                *det_before.entry(id).or_default() += 1;
            }
        }
        let mut det_after: HashMap<u32, usize> = HashMap::new();
        for k in 0..bs.lake_map.len() {
            if bs.lake_map[k] != 0 {
                drainage.lake_map[k] = bs.lake_map[k];
            } else if drainage.lake_map[k] != 0 && drainage.lake_map[k] < 1_000_001 {
                *det_after.entry(drainage.lake_map[k]).or_default() += 1;
            }
        }
        let absorbed: HashSet<u32> = det_before
            .keys()
            .copied()
            .filter(|id| det_after.get(id).copied().unwrap_or(0) == 0)
            .collect();
        let n_absorbed = absorbed.len();
        drainage.lakes.retain(|lk| lk.base.id >= 1_000_001 || !absorbed.contains(&lk.base.id));
        let (to_sea, chained) = (
            bs.spillways.iter().filter(|s| s.chained_into.is_none()).count(),
            bs.spillways.iter().filter(|s| s.chained_into.is_some()).count(),
        );
        drainage.lakes.extend(bs.lakes);
        if n_absorbed > 0 {
            eprintln!(
                "[HD] below-sea cleanup: {n_absorbed} detected lake(s) submerged by a filled below-sea lake -> dropped"
            );
        }
        eprintln!(
            "[HD] below-sea basins: {} lakes ({exo} exorheic, {endo} endorheic); {} spillways ({to_sea} -> sea, {chained} chained)",
            exo + endo,
            bs.spillways.len()
        );
        (bs.wetland, bs.spillways)
    };
    let before_seg = drainage.rivers.segments.len();
    clip_rivers_to_lakes(&mut drainage);
    eprintln!(
        "[HD] rivers clipped to lakes: {before_seg} -> {} segments (terminate at sinks)",
        drainage.rivers.segments.len()
    );
    for sw in &below_sea_spillways {
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
    apply_geo_scale_ratio(&mut drainage, geo_scale_ratio, &dcfg.thresholds);
    if geo_scale_ratio != 1.0 {
        eprintln!(
            "[HD] geographic scale ratio {geo_scale_ratio:.2} -> hydrology signifies x{:.1} area",
            geo_scale_ratio * geo_scale_ratio
        );
    }
    HdDrainageBundle { drainage, wetland: wetland_mask }
}

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
    // THE SHIPPED PATH goes through ONE core function — `production_hd_config` — so a bench
    // and production cannot diverge (they did seven times; see ADR "The DEAD KNOB"). The
    // experimental relief-v1/v2 / cross-rill branches below keep the legacy mutation chain,
    // being opt-ins that never ship. `sample_origin` is assigned after the framing roll is
    // computed (it is an OPT of the function, not a drift).
    let ships_relief_v3 =
        params.stream_power && params.closures && params.mfd && !params.cross_rill;
    let mut upscale = if ships_relief_v3 {
        ymir_core::terrain::upscale::production_hd_config(
            &ymir_core::terrain::upscale::ProductionHdOpts {
                target_size: params.target_size,
                domain_km: params.domain_km,
                depth_scale_m: ss.depth_scale_m as f32,
                sample_origin: [0.0, 0.0], // set below, once the framing roll is known
                sample_size: 1.0,
                // ⚠️ INERT on this path (the C-1 relief-budget cap binds everywhere) — kept
                // so the config states the intended level. See ADR "The DEAD KNOB".
                amplitude_base: params.fbm_amplitude.unwrap_or(0.04),
                mfd_p: params.mfd_p,
                lithology: params.lithology.clone().unwrap_or_default(),
                fracture: {
                    let mut f = params.fracture.clone().unwrap_or_default();
                    f.domain_km = params.domain_km;
                    f
                },
            },
        )
    } else {
        let mut u = FbmUpscaleConfig::c1_hd_production(params.target_size);
        if let Some(a) = params.fbm_amplitude {
            u.amplitude_base = a; // striation-ladder override — INERT while conditioning is on
        }
        u
    };
    // EXPERIMENTAL opt-in (ADR 0001): swap the droplet pass for routed stream-power
    // incision + hillslope regime split. A_c is 7.6 km² expressed in cells at THIS
    // resolution (resolution-stable channel head); uncoupled vertical scale.
    if ships_relief_v3 {
        eprintln!(
            "[HD] SHIPPED relief-v3 config from production_hd_config (MFD p={:.1}, C-3 litho {}, C-3b fracture {})",
            params.mfd_p,
            if upscale.lithology.enabled { "ON" } else { "off" },
            if upscale.fracture.enabled { "ON" } else { "off" }
        );
    }
    if params.stream_power && !ships_relief_v3 {
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
        // C-3 lithology — attach the per-cell erodibility multiplier config (the K
        // field is built inside upscale_from_c1). Only meaningful with stream-power on
        // (this branch). None/disabled → uniform K, byte-identical.
        if let Some(litho) = &params.lithology {
            upscale.lithology = litho.clone();
            if litho.enabled {
                eprintln!(
                    "[HD] C-3 lithology ON (soft rift ×{:.0}, volcaniclastic ×{:.0})",
                    litho.soft_multiplier, litho.volcanic_multiplier
                );
            }
        }
        // C-3b inherited structure — fracture-density erodibility (domain_km = the
        // geometric span, like volcanism). None/disabled → uniform, byte-identical.
        if let Some(mut fr) = params.fracture.clone() {
            fr.domain_km = params.domain_km;
            if fr.enabled {
                eprintln!(
                    "[HD] C-3b fracture ON (amplitude ×{:.0}, decay {:.0} km)",
                    fr.amplitude, fr.decay_km
                );
            }
            upscale.fracture = fr;
        }
    }
    // Land report FIRST (reads only `target_land_fraction`) so the seam-correct
    // torus centre of the largest mass is known before we choose the sampling roll.
    let t_land = Instant::now();
    let land = c1_coarse_land_report(
        spec.seed,
        spec.grid_size,
        &spec.init_params,
        &run,
        &spec.closures,
        &ss,
        upscale.target_land_fraction,
    );
    eprintln!("[HD timing] coarse land report {:.1}s", t_land.elapsed().as_secs_f32());
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
    // Finding 37 POINT 2 — extend exported watercourses UPSTREAM from the 20 km² extraction threshold
    // to the fluvial-regime critical area A_c (the channel head), MAIN-STEM only (author's choice:
    // one headwater tail per watercourse; ×33 segments vs the full tree's ×182). The watercourse
    // COUNT is unchanged; only their upstream extent grows, restoring a real source→mouth width range.
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    // H-1: the infiltration TUNABLES ride in the drainage config so BOTH drainage cache
    // keys fold them (the per-cell field is derived from the same seed/tectonic configs
    // already covered by the eroded key). None → key byte-identical to pre-H-1.
    dcfg.infiltration = params.infiltration.clone().filter(|c| c.enabled);
    let cache_dir = default_cache_dir();
    let lat = params.latitude_deg;

    // C-2 volcanism config for this run — geometric domain_km (never geo_scale_ratio,
    // which shapes nothing). None → disabled → byte-identical terrain.
    let volc = {
        let mut v = params
            .volcanism
            .clone()
            .unwrap_or_else(ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig::default);
        v.domain_km = params.domain_km;
        v
    };

    // Cache keys (for HIT/MISS detection via sidecar existence). ekey MUST include
    // volcanism identically to `cached_c1_eroded` (shared `eroded_key_full`), or a
    // volcanism-on/off pair would share the drainage cache derived from ekey.
    let tkey = tectonic_key(spec.seed, spec.grid_size, &spec.init_params, &run, &spec.closures);
    let ekey = eroded_key_full(&tkey, &ss, &upscale, &volc);
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
    let t_ero = Instant::now();
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
            &volc,
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
            &volc,
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
    // Split the bundle: the heightmap flows on as before; the craters (empty when
    // volcanism is off) travel to the lake-typing stage. Bundled with the terrain
    // through the cache, so a HIT carries its craters too (no silent mistyping).
    let craters = eroded.craters;
    let eroded = eroded.heightmap;
    if volc.enabled {
        // Always log when enabled (0 records at a coarse target ≠ "did not run" —
        // per-mechanism edifice counts print from the core miss path above).
        let active = craters.iter().filter(|c| c.active).count();
        eprintln!(
            "[HD volcanism] {} resolved crater records ({} active/acidic, {} extinct)",
            craters.len(),
            active,
            craters.len() - active
        );
    }
    eprintln!(
        "[HD timing] eroded {} ({:.1}s)",
        if eroded_hit { "HIT" } else { "MISS" },
        t_ero.elapsed().as_secs_f32()
    );

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
    // #190 — the hd_drainage bundle key, computed EARLY so we can tell whether the pre-breach drainage
    // is still needed (it is not when BOTH the conditioned eroded and the final drainage are HITs).
    let hd_dr_key = hd_drainage_key(
        &ekey,
        &dcfg,
        &ss,
        lat,
        &pp,
        params.latitude_span_deg,
        params.geo_scale_ratio,
        window_km,
    );
    let bundle_hit = sidecar_exists("hd_drainage", hd_dr_key.digest());
    let (eroded, prebreach_override) = if params.stream_power && params.closures && params.mfd {
        let (gw, gh) = (eroded.width, eroded.height);
        // #190 — the CONDITIONED eroded (the breach output) is CACHED: a pure GridF32 keyed on the
        // pre-breach drainage key + the breach version. On a HIT the ~50 s breach is skipped entirely.
        let prebreach_key = drainage_key_windowed(&ekey, &dcfg, &ss, None, window_km);
        let cond_key = conditioned_eroded_key(&prebreach_key);
        let cond_hit = sidecar_exists("conditioned", cond_key.digest());
        // The pre-breach drainage is climate-independent (cached on the eroded key). It is needed ONLY
        // to RUN the breach (cond MISS) or to carry the pre-breach lakes into build_hd_drainage
        // (bundle MISS); on a full HIT it is not loaded at all.
        let prebreach = if !cond_hit || !bundle_hit {
            let t_pre = Instant::now();
            let p = match cached_c1_drainage_windowed(
                &cache_dir, &ekey, &eroded, None, &dcfg, &ss, window_km,
            ) {
                Ok(p) => p,
                Err(e) => {
                    let _ = tx.send(C1Event::HdFailed { error: format!("prebreach: {e}") });
                    return;
                }
            };
            eprintln!("[HD timing] prebreach drainage {:.1}s", t_pre.elapsed().as_secs_f32());
            Some(p)
        } else {
            None
        };
        // C-2: protect ACTIVE crater bowls from the breach so they survive as closed
        // depressions for the (climate-dependent) crater-lake stage. Climate-INDEPENDENT
        // (crater cells only), so it belongs in this cached, climate-free conditioning.
        let crater_protect: Option<Vec<bool>> = if volc.enabled {
            let mut m = vec![false; gw * gh];
            for c in craters.iter().filter(|c| c.active) {
                let (cx, cy, r) = (c.center_px.0, c.center_px.1, c.radius_px.max(1.0));
                let (i0, i1) =
                    ((cx - r).floor().max(0.0) as usize, ((cx + r).ceil() as usize).min(gw - 1));
                let (j0, j1) =
                    ((cy - r).floor().max(0.0) as usize, ((cy + r).ceil() as usize).min(gh - 1));
                for j in j0..=j1 {
                    for i in i0..=i1 {
                        if ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt()
                            <= r
                        {
                            m[j * gw + i] = true;
                        }
                    }
                }
            }
            Some(m)
        } else {
            None
        };
        let t_breach = Instant::now();
        let conditioned = match cached::<GridF32>(&cache_dir, "conditioned", &cond_key, || {
            let p = prebreach.as_ref().expect("prebreach is present on a conditioned MISS");
            breach_monotone_protected(
                &eroded,
                &p.flow.filled,
                &p.lake_map,
                0.5,
                gw,
                gh,
                crater_protect.as_deref(),
            )
        }) {
            Ok(g) => g,
            Err(e) => {
                let _ = tx.send(C1Event::HdFailed { error: format!("breach: {e}") });
                return;
            }
        };
        eprintln!(
            "[HD] relief-v3 breach conditioning {} ({:.1}s)",
            if cond_hit { "HIT" } else { "MISS" },
            t_breach.elapsed().as_secs_f32()
        );
        (conditioned, prebreach)
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
    // #190 — the climate is cached: keyed on (ekey, lat, span, precip params) it HITs when the same
    // (latitude, span) is re-selected. The conditioned `eroded` it reads is deterministic from ekey.
    let climate_hit = sidecar_exists(
        "climate",
        climate_key(&ekey, lat, params.latitude_span_deg, &pp, window_km).digest(),
    );
    let climate = match cached_c1_climate(
        &cache_dir,
        &ekey,
        &eroded,
        &ss,
        lat,
        params.latitude_span_deg,
        &pp,
        window_km,
    ) {
        Ok(c) => c,
        Err(e) => {
            let _ = tx.send(C1Event::HdFailed { error: format!("climate: {e}") });
            return;
        }
    };
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Climate,
        regime: if climate_hit { CacheRegime::Hit } else { CacheRegime::Miss },
        elapsed: t.elapsed(),
    });
    eprintln!(
        "[HD timing] climate {} ({:.1}s)",
        if climate_hit { "HIT" } else { "MISS" },
        t.elapsed().as_secs_f32()
    );

    // ── Phase 3: drainage (rivers + lakes + water balance). ──
    bail_if_cancelled!();
    // `hd_dr_key` + `bundle_hit` were computed EARLY (before the breach) to gate the pre-breach load.
    let drainage_hit = bundle_hit;
    let _ = tx.send(C1Event::HdPhaseStarted { phase: HdPhase::Drainage });
    let t = Instant::now();
    // #190 — the FINAL drainage (base + below-sea merge + clip + spillways + geo-scale) and the
    // wetland mask are cached as ONE bundle under a COMPLETE key; the same params re-selected HIT.
    // H-1 — build the per-cell INFILTRATED FRACTION from the causal permeability
    // (lithology matrix + fracture density, double porosity; no slope term — the
    // literature does not support one). Costs one cheap coarse tectonic pass. Disabled →
    // `None` → the pre-H-1 balance, byte-identical.
    let infil_field: Option<Vec<f32>> = dcfg.infiltration.as_ref().map(|ic| {
        use ymir_core::tectonics_c1::closures::infiltration::build_hd_infiltration;
        use ymir_core::tectonics_c1::closures::volcanism::place_edifices;
        use ymir_core::tectonics_c1::debug_labels::run_coarse_tectonics;
        let (st, kn) = run_coarse_tectonics(
            spec.seed,
            spec.grid_size,
            &spec.init_params,
            &run,
            &spec.closures,
        );
        let edi = if volc.enabled {
            place_edifices(
                &st,
                &kn,
                &ymir_core::seed::WorldSeed::new(spec.seed),
                volc.domain_km,
                &volc,
            )
        } else {
            Vec::new()
        };
        let km_per_cell = upscale.sample_size as f32 * params.domain_km / eroded.width as f32;
        eprintln!(
            "[HD] H-1 infiltration ON (f_cap {:.2}, k_ref {:.1e} m/day)",
            ic.f_cap, ic.k_ref_m_per_day
        );
        build_hd_infiltration(
            &st,
            &kn,
            &upscale.lithology,
            &upscale.fracture,
            ic,
            &edi,
            eroded.width,
            eroded.height,
            upscale.sample_origin,
            upscale.sample_size,
            km_per_cell,
        )
    });
    let (mut drainage, wetland_mask) = match cached(&cache_dir, "hd_drainage", &hd_dr_key, || {
        build_hd_drainage(
            &eroded,
            &climate,
            prebreach_override,
            &dcfg,
            &ss,
            window_km,
            params.geo_scale_ratio,
            infil_field.as_deref(),
        )
    }) {
        Ok(b) => (b.drainage, b.wetland),
        Err(e) => {
            let _ = tx.send(C1Event::HdFailed { error: format!("drainage: {e}") });
            return;
        }
    };
    // Finding 37 — the exorheic-outlet invariant, over the WHOLE lake population (both provenances),
    // on the PRODUCTION network (not a synthetic-grid subset). A non-empty result means a lake was
    // labelled exorheic without an emitted outflow — the label is wrong. Loud in dev, hard in debug.
    let missing = exorheic_lakes_missing_outlet(&drainage);
    if !missing.is_empty() {
        eprintln!(
            "[HD] WARNING Finding 37: {} exorheic lake(s) with NO traced outlet: {:?}",
            missing.len(),
            &missing[..missing.len().min(12)]
        );
    }
    debug_assert!(missing.is_empty(), "exorheic lakes without a traced outlet: {missing:?}");
    let _ = tx.send(C1Event::HdPhaseDone {
        phase: HdPhase::Drainage,
        regime: if drainage_hit { CacheRegime::Hit } else { CacheRegime::Miss },
        elapsed: t.elapsed(),
    });
    eprintln!(
        "[HD timing] drainage {} ({:.1}s)",
        if drainage_hit { "HIT" } else { "MISS" },
        t.elapsed().as_secs_f32()
    );

    // ── C-2 crater-lake DETECTION + typing. The generic lake detection discards a
    // crater as sub-threshold (~1.2 km² < the 5 km² noise floor), so a dedicated
    // pass fills the ACTIVE craters (reconstructed → closed bowls) with the SAME
    // water balance and adds them as CraterAcidic lakes. Done AFTER the drainage
    // invariant checks and BEFORE export/biomes (which read the updated lake_map).
    if volc.enabled {
        use ymir_core::tectonics_c1::drainage::{DrainageClimate, runoff_accumulation};
        use ymir_core::terrain::flow::{FlowConfig, compute_flow};
        let (gw, gh) = (eroded.width, eroded.height);
        let cell_km2 = (window_km / gw as f32).powi(2);
        let flow = compute_flow(
            &eroded,
            &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false },
        );
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let runoff = runoff_accumulation(&eroded, &flow, &dclim, cell_km2, None, None, gw, gh);
        let (held, dry) = ymir_core::tectonics_c1::closures::volcanism::detect_crater_lakes(
            &eroded,
            &flow.filled,
            &runoff,
            &climate.temperature,
            &craters,
            cell_km2,
            &ss,
            &mut drainage.lakes,
            &mut drainage.lake_map,
        );
        eprintln!("[HD volcanism] crater lakes: {held} acidic | {dry} dry craters");
    }

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
    eprintln!("[HD timing] biomes ({:.1}s)", t.elapsed().as_secs_f32());

    // ── Optional: write the v1 `.ymir` delivery container. ──
    // Explicit opt-in only (never automatic). Ships height (placeholder) +
    // coastline/cliffs (Y-B) + temperature/precipitation/biome (Y-C).
    if let Some(export_dir) = &params.export_dir {
        // Resolved climatic span (explicit, else the geographic domain_km/111) so the
        // manifest records the actual gradient behind the climate/biome layers.
        let lat_span = params.latitude_span_deg.unwrap_or(window_km / 111.0);
        let t_exp = Instant::now();
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
        eprintln!("[HD timing] export .ymir ({:.1}s)", t_exp.elapsed().as_secs_f32());
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

    // ── Debug microscope (C-2/C-3): derive the coarse tectonic labels + volcanic
    // footprints for the overlay. One extra ~1 s coarse pass, opt-in. The labels come
    // from the SAME derivation production erodes (deterministic in seed), so they
    // register with the terrain via `sample_origin`/`sample_size`.
    let km_per_cell = upscale.sample_size as f32 * params.domain_km / eroded.width as f32;
    let (tectonic, edifices) = if params.emit_tectonic_labels {
        use ymir_core::tectonics_c1::closures::volcanism::place_edifices;
        use ymir_core::tectonics_c1::debug_labels::{derive_tectonic_labels, run_coarse_tectonics};
        let (state, kin) = run_coarse_tectonics(
            spec.seed,
            spec.grid_size,
            &spec.init_params,
            &run,
            &spec.closures,
        );
        let labels = derive_tectonic_labels(&state, &kin);
        let edi = if volc.enabled {
            place_edifices(
                &state,
                &kin,
                &ymir_core::seed::WorldSeed::new(spec.seed),
                volc.domain_km,
                &volc,
            )
        } else {
            Vec::new()
        };
        eprintln!("[HD microscope] tectonic labels derived ({} edifices)", edi.len());
        (Some(labels), edi)
    } else {
        (None, Vec::new())
    };

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
        volcanoes: craters,
        tectonic,
        edifices,
        sample_origin: upscale.sample_origin,
        sample_size: upscale.sample_size,
        km_per_cell,
    });
    eprintln!(
        "[HD timing] run_hd TOTAL {:.1}s (render is separate, on the UI thread)",
        t_all.elapsed().as_secs_f32()
    );
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
