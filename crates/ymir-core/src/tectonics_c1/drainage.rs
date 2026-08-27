//! C1-product drainage extraction (#155 complétude).
//!
//! **This does NOT re-implement drainage.** The hydraulic-erosion closure has
//! already carved the river network into the terrain; the geometry here READS
//! that network deterministically. It composes the existing, unit-tested
//! [`crate::terrain::flow`] (priority-flood pit-fill → D8 → flow accumulation →
//! Strahler river extraction) and [`crate::lakes::detection`] (sill-level lakes,
//! `surface_elevation = outlet level` — the physical rule) on the C1 HD product,
//! wired to the C1 **coordinate contract**:
//!
//! - **sea level = 0.5** ([`C1_SEA_LEVEL_NORM`], the Maillon-2 unified scale) —
//!   NOT the legacy 0.1 default the viz hydrology still carries.
//! - **km² network thresholds** (horizontal pin `C1_DOMAIN_KM` →
//!   [`crate::tectonics_c1::production_upscale::c1_cell_area_km2`]) — drainage
//!   area in km², resolution-INdependent, where raw accumulation cell-counts
//!   break across grid sizes.
//! - **lake stats in metres** (vertical pin
//!   [`crate::tectonics_c1::production_upscale::c1_altitude_norm_to_metres`]).
//!
//! The viz hydrology phase runs the same stack on the *v2* interactive erosion
//! cache (a different path); rewiring it to consume this C1-product entry is a
//! DISTINCT maillon (not done here).

use serde::{Deserialize, Serialize};

use crate::climate::precipitation::{e_sat, precip_mm_per_year};
use crate::grid::GridF32;
use crate::lakes::detection::{Lake, LakeConfig, detect_lakes};
use crate::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlatPerturbation, FlowConfig, FlowResult, RiverConfig, RiverNetwork,
    RiverSegment, compute_flow, extract_rivers,
};

use super::closures::oceanic_bathymetry::params::SteinSteinParams;
use super::production_upscale::{C1_DOMAIN_KM, c1_altitude_norm_to_metres};

/// #drainage — climate input for the water balance (the "hydroclimate layer" the
/// geometric fill-and-spill placeholder was waiting for). `Some` couples the
/// balance (basins overflow only if inflow > evaporation; arid basins become
/// terminal/endorheic); `None` → the pure-geometry path (byte-identical pre-fix).
pub struct DrainageClimate<'a> {
    /// Precipitation in `c1_climate` INTERNAL units (→ mm/yr via `precip_mm_per_year`).
    pub precip_internal: &'a GridF32,
    /// Air temperature (°C) — drives the potential-evaporation proxy.
    pub temperature: &'a GridF32,
}

/// #drainage — potential OPEN-WATER evaporation (mm/yr) from temperature: the
/// Clausius-Clapeyron proxy `PE = PE_PER_ESAT · e_sat(T)` (warmer air evaporates
/// more). Anchored on the observed open-water PE range — temperate ~12 °C →
/// ~850 mm/yr, hot desert ~27 °C → ~2200 mm/yr, cold ~0 °C → ~370 mm/yr.
const PE_PER_ESAT: f32 = 61.0;

/// Potential evaporation (mm/yr) at air temperature `t_c` (°C).
pub fn potential_evaporation_mm(t_c: f32) -> f32 {
    PE_PER_ESAT * e_sat(t_c)
}

/// The C1 unified sea level (Maillon 2): continental sea maps to 0.5.
pub const C1_SEA_LEVEL_NORM: f32 = 0.5;

/// River navigability class (TDD §9.2), assigned by **upstream drainage area in
/// km²** — the discharge proxy until climate/runoff couples (a larger basin
/// carries more water). NOT raw cell-counts (those break across resolutions).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Navigability {
    /// Below the small-boat threshold (creeks / minor streams).
    NonNavigable,
    /// Small craft (canoe / skiff).
    SmallBoat,
    /// Barge / riverboat.
    Barge,
    /// Ocean-going ship up-river.
    Ship,
}

/// Drainage classification of a lake by its drainage relationship to the sea.
///
/// **Honest placeholder.** Priority-flood pit-filling routes every depression to
/// an overflow sill that ultimately spills to the ocean, so PURE GEOMETRY yields
/// `Exorheic` for all lakes. A TRUE endorheic (terminal) lake — where evaporation
/// balances inflow so it never overflows — is a CLIMATE phenomenon (precip vs
/// evaporation), not geometry; this field will carry `Endorheic` only once a
/// hydroclimate layer couples. Computed honestly here (trace the outlet
/// downstream): if it reaches the sea → `Exorheic`, else `Endorheic`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LakeType {
    /// Drains to the sea (the geometric default).
    Exorheic,
    /// Terminal — no path to the sea (requires climate to actually occur).
    Endorheic,
}

/// km² drainage-area thresholds for the navigability classes + the minimal
/// mapped channel. **Anchored on real-world river navigability** (drainage area
/// as the discharge proxy), NOT tuned — revisable by visual density review.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DrainageThresholds {
    /// Minimal mapped channel (a cell is part of the river network at/above this
    /// upstream area). ~20 km² ≈ a small perennial stream's catchment.
    pub stream_km2: f32,
    /// Small-boat navigable from this drainage area up. ~500 km².
    pub small_boat_km2: f32,
    /// Barge navigable. ~5 000 km².
    pub barge_km2: f32,
    /// Ship navigable. ~50 000 km².
    pub ship_km2: f32,
    /// ADR 0001 Finding 37 — extend retained watercourses UPSTREAM to this area (the channel head).
    /// `0` = no extension (the first exported point stays at `stream_km2`, byte-identical). The
    /// natural value is the erosion regime-split critical area A_c (`RELIEF_V1_A_C_KM2`, ~0.1 km²).
    #[serde(default)]
    pub head_km2: f32,
    /// With extension on: `true` ramifies the whole upstream tree to `head_km2`; `false` extends the
    /// main stem only. Defaults to `true` (dense tree; render-time filterable by Strahler order).
    #[serde(default = "default_true")]
    pub full_tree: bool,
}

fn default_true() -> bool {
    true
}

impl Default for DrainageThresholds {
    fn default() -> Self {
        Self {
            stream_km2: 20.0,
            small_boat_km2: 500.0,
            barge_km2: 5_000.0,
            ship_km2: 50_000.0,
            head_km2: 0.0,
            full_tree: true,
        }
    }
}

impl DrainageThresholds {
    fn classify(&self, drainage_km2: f32) -> Navigability {
        if drainage_km2 >= self.ship_km2 {
            Navigability::Ship
        } else if drainage_km2 >= self.barge_km2 {
            Navigability::Barge
        } else if drainage_km2 >= self.small_boat_km2 {
            Navigability::SmallBoat
        } else {
            Navigability::NonNavigable
        }
    }
}

/// C1 drainage configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct C1DrainageConfig {
    /// Network / navigability thresholds in km².
    pub thresholds: DrainageThresholds,
    /// Minimum lake depth in METRES (filled − eroded) to classify as a lake.
    pub lake_min_depth_m: f32,
    /// Minimum lake surface area in km² to keep.
    pub lake_min_area_km2: f32,
    /// #drainage fix A — micro-relief restored on the pit-filled flats so rivers
    /// wander instead of running cardinal-straight. `None` → legacy routing
    /// (byte-identical). Folded into the drainage cache key via the config.
    pub flat_perturbation: Option<FlatPerturbation>,
    /// #drainage fix A — D∞ continuous-direction flow on the flats (Tarboton).
    /// `false` → mono D8 (byte-identical). See `FlowConfig::dinf`.
    pub dinf: bool,
}

impl Default for C1DrainageConfig {
    fn default() -> Self {
        // 10 m min lake depth, 5 km² min lake area — plausible mappable lakes.
        Self {
            thresholds: DrainageThresholds::default(),
            lake_min_depth_m: 10.0,
            lake_min_area_km2: 5.0,
            flat_perturbation: Some(FlatPerturbation::default()),
            dinf: false,
        }
    }
}

/// A lake enriched with the C1 coordinate contract (metres + area km² + type).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct C1Lake {
    /// The underlying geometric lake (norm-space surface/depth, outlet, basin).
    pub base: Lake,
    /// Surface (outlet sill) elevation in metres.
    pub level_m: f32,
    /// Maximum depth in metres.
    pub depth_m: f32,
    /// Surface area in km².
    pub area_km2: f32,
    /// Drainage type (Exorheic geometric default; Endorheic needs climate).
    pub lake_type: LakeType,
}

/// Full C1 drainage product: flow field, rivers (+ per-segment navigability /
/// drainage area), and lakes (in metres / km² / typed).
#[derive(Clone)]
pub struct C1DrainageResult {
    pub flow: FlowResult,
    pub rivers: RiverNetwork,
    /// Per-segment upstream drainage area in km² (parallel to `rivers.segments`).
    pub segment_drainage_km2: Vec<f32>,
    /// Per-segment navigability class (parallel to `rivers.segments`).
    pub segment_navigability: Vec<Navigability>,
    /// Per-segment mean DISCHARGE in m³/s (parallel to `rivers.segments`): the runoff
    /// (`precip − PE`) accumulated over the reach's catchment, the water the channel
    /// actually carries. Drives `segment_width_m` (Finding 22).
    pub segment_discharge_m3s: Vec<f32>,
    /// Per-segment bankfull channel WIDTH in metres (parallel to `rivers.segments`),
    /// hydraulic geometry `w = CHANNEL_WIDTH_A · Q^CHANNEL_WIDTH_B` on the DISCHARGE.
    pub segment_width_m: Vec<f32>,
    /// Per-segment LONG PROFILE: bed elevation (metres) at each of the segment's
    /// points, upstream→downstream (parallel to `rivers.segments[i].points`).
    pub segment_profile_m: Vec<Vec<f32>>,
    pub lakes: Vec<C1Lake>,
    pub lake_map: Vec<u32>,
    pub width: usize,
    pub height: usize,
}

/// Extract the drainage network (rivers + lakes) from the C1 HD product
/// `heightmap` (the post-erosion `upscale_from_c1` output), at the C1 unified
/// sea level (0.5) with km² thresholds and metre lake stats.
///
/// `ss` supplies the vertical scale (`depth_scale_m`) for the norm→metre lake
/// stats. The km² scale comes from the heightmap's own resolution
/// (`c1_cell_area_km2(width)`).
pub fn c1_drainage(
    heightmap: &GridF32,
    climate: Option<&DrainageClimate>,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
) -> C1DrainageResult {
    c1_drainage_windowed(heightmap, climate, cfg, ss, C1_DOMAIN_KM)
}

/// [`c1_drainage`] for a grid spanning `window_km` (a cropped playable window).
/// The cell area `(window_km / width)²` drives the km² network / lake-area
/// thresholds, so a zoomed window's rivers/lakes use its OWN metric scale.
/// `window_km == C1_DOMAIN_KM` reproduces [`c1_drainage`] exactly.
pub fn c1_drainage_windowed(
    heightmap: &GridF32,
    climate: Option<&DrainageClimate>,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    window_km: f32,
) -> C1DrainageResult {
    let (w, h) = (heightmap.width, heightmap.height);
    let cell_km2 = {
        let s = window_km / w as f32;
        s * s
    };
    // metres per unit of normalised altitude (the linear vertical scale slope).
    let metres_per_norm = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);

    // 1. Flow at the unified sea level (NOT the legacy 0.1). The flat
    //    perturbation (fix A) restores micro-relief on the pit-filled flats so the
    //    routing wanders; it touches only the routing surface, not the hydrology.
    let flow = compute_flow(
        heightmap,
        &FlowConfig {
            sea_level: C1_SEA_LEVEL_NORM,
            flat_perturbation: cfg.flat_perturbation.clone(),
            dinf: cfg.dinf,
        },
    );

    // 2. Rivers — km² thresholds → accumulation cell-counts for THIS resolution.
    let to_cells = |km2: f32| (km2 / cell_km2).max(1.0);
    let river_cfg = RiverConfig {
        stream_threshold: to_cells(cfg.thresholds.stream_km2),
        // river/major are display tiers in terrain::flow; keep them at the
        // small-boat / barge km² tiers so the existing consumers stay coherent.
        river_threshold: to_cells(cfg.thresholds.small_boat_km2),
        major_river_threshold: to_cells(cfg.thresholds.barge_km2),
        // Finding 37 — extend retained watercourses upstream to A_c (fluvial-regime start), or 0 to
        // keep the first exported point at `stream_km2` (byte-identical). Driven by the config.
        head_threshold: if cfg.thresholds.head_km2 > 0.0 {
            to_cells(cfg.thresholds.head_km2)
        } else {
            0.0
        },
        full_tree: cfg.thresholds.full_tree,
    };
    let rivers = extract_rivers(&flow, &river_cfg, w, h);

    // 3. Lakes — min_depth (m → norm), min_area (km² → cells).
    let lake_cfg = LakeConfig {
        min_depth: cfg.lake_min_depth_m / metres_per_norm,
        min_area: (cfg.lake_min_area_km2 / cell_km2).ceil().max(1.0) as usize,
    };
    let lake_result =
        detect_lakes(heightmap, &flow.filled, &flow.direction, &flow.basins, &lake_cfg);

    // Enrich lakes. With a `climate` (#drainage hydroclimate layer) the WATER
    // BALANCE decides each basin: overflow → exorheic, else shrink to the
    // endorheic equilibrium level (drained cells leave `lake_map`). Without it,
    // the pure-geometry path (every basin fills to its sill → exorheic) is kept
    // byte-identical.
    let mut lake_map = lake_result.lake_map;
    let lakes: Vec<C1Lake> = match climate {
        Some(clim) => water_balance_lakes(
            heightmap,
            &flow,
            clim,
            cell_km2,
            ss,
            &lake_result.lakes,
            &mut lake_map,
            w,
            h,
        ),
        None => lake_result
            .lakes
            .iter()
            .map(|lk| {
                let level_m = c1_altitude_norm_to_metres(lk.surface_elevation, ss);
                let floor_m = c1_altitude_norm_to_metres(lk.surface_elevation - lk.max_depth, ss);
                let lake_type = if outlet_reaches_sea(lk, &flow, heightmap, w, h) {
                    LakeType::Exorheic
                } else {
                    LakeType::Endorheic
                };
                C1Lake {
                    base: lk.clone(),
                    level_m,
                    depth_m: level_m - floor_m,
                    area_km2: lk.area as f32 * cell_km2,
                    lake_type,
                }
            })
            .collect(),
    };

    // 4. Per-segment drainage + navigability. With `climate` (#drainage fix B) the
    // DISCHARGE is the REAL runoff (`max(0, precip−PE)` accumulated downstream),
    // with endorheic basins as SINKS (water dies in the closed lake — no phantom
    // river below it). A segment's discharge = max runoff_accum over its cells;
    // its effective drainage area = discharge / a reference runoff depth, so the
    // existing km² navigability thresholds apply unchanged at that reference.
    // Effect: desert channels (no surplus) and rivers exiting a closed basin → 0
    // discharge → NonNavigable; allochthonous rivers (humid upstream) keep their
    // accumulated discharge across a dry reach (the Nile). Without `climate` →
    // geometric cell-count drainage area (byte-identical).
    // Also carry the per-segment DISCHARGE in m³/s (`runoff·area`, the water actually
    // routed to the reach) — the physical quantity the channel-width law needs (see
    // Finding 22). `dr_km2` stays the effective-area proxy that drives the km²
    // navigability thresholds; discharge drives width.
    let (segment_drainage_km2, segment_navigability, segment_discharge_m3s): (
        Vec<f32>,
        Vec<Navigability>,
        Vec<f32>,
    ) = match climate {
        Some(clim) => {
            let mut endorheic = vec![false; w * h];
            for lk in &lakes {
                if lk.lake_type == LakeType::Endorheic {
                    for k in 0..w * h {
                        if lake_map[k] == lk.base.id {
                            endorheic[k] = true;
                        }
                    }
                }
            }
            let discharge =
                runoff_accumulation(heightmap, &flow, clim, cell_km2, Some(&endorheic), w, h);
            // Per segment: `q` = max runoff_accum over its cells (mm·km²/yr). The reach's
            // downstream-most cell carries the largest accumulation, so `max` = the
            // cumulative catchment discharge at the mouth of the reach — and because
            // `runoff_accumulation` carries flow ACROSS an exorheic lake (flat, routed to
            // the outlet), a lake's outlet reach inherits the whole upstream catchment
            // automatically (Finding 22, TASK 2). An endorheic reach reads 0 (the sink
            // kills propagation), so its channel width is 0.
            let mut dr_km2 = Vec::with_capacity(rivers.segments.len());
            let mut q_m3s = Vec::with_capacity(rivers.segments.len());
            for s in &rivers.segments {
                let q = s
                    .points
                    .iter()
                    .map(|&(x, y)| discharge[y as usize * w + x as usize])
                    .fold(0.0f32, f32::max);
                dr_km2.push(q / REFERENCE_RUNOFF_MM); // discharge → effective km² at reference depth
                q_m3s.push(runoff_km2_to_m3s(q));
            }
            let nav = dr_km2.iter().map(|&km2| cfg.thresholds.classify(km2)).collect();
            (dr_km2, nav, q_m3s)
        }
        None => {
            let dr_km2: Vec<f32> = rivers.segments.iter().map(|s| s.max_flow * cell_km2).collect();
            let nav = dr_km2.iter().map(|&km2| cfg.thresholds.classify(km2)).collect();
            // No climate → assume the reference runoff depth over the geometric area so a
            // discharge-consistent width still results (Q = REFERENCE_RUNOFF · area).
            let q_m3s =
                dr_km2.iter().map(|&km2| runoff_km2_to_m3s(REFERENCE_RUNOFF_MM * km2)).collect();
            (dr_km2, nav, q_m3s)
        }
    };

    // Finding 22 (DEFECT C corrected) — channel WIDTH from DISCHARGE, not area. Hydraulic
    // geometry `w = a·Q^b` (Q in m³/s). A dry/endorheic reach (Q=0) → 0 width (no channel).
    // The long profile is the bed elevation (m) along the segment's own points, up→down.
    let segment_width_m: Vec<f32> = segment_discharge_m3s
        .iter()
        .map(|&q| CHANNEL_WIDTH_A * q.max(0.0).powf(CHANNEL_WIDTH_B))
        .collect();
    let segment_profile_m: Vec<Vec<f32>> = rivers
        .segments
        .iter()
        .map(|s| {
            s.points
                .iter()
                .map(|&(x, y)| {
                    c1_altitude_norm_to_metres(heightmap.data[y as usize * w + x as usize], ss)
                })
                .collect()
        })
        .collect();

    C1DrainageResult {
        flow,
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
    }
}

/// Finding 22 — hydraulic-geometry width law `w = a·Q^b` (Leopold & Maddock),
/// with `Q` the real DISCHARGE in m³/s (NOT the drainage area — the earlier law
/// fed km² into a coefficient calibrated for m³/s, compressing the whole
/// distribution). `b = 0.5` is the classic downstream bankfull-width exponent;
/// `a = 5.0` is the mid-range natural-channel coefficient (a·√Q). Anchored on
/// mean-annual discharge, so a Thames-scale trunk (~65 m³/s mean) reads ~40 m and
/// a headwater (~0.1 m³/s) ~1.5 m — the trunk/headwater RATIO (~50×) is the point.
const CHANNEL_WIDTH_A: f32 = 5.0;
const CHANNEL_WIDTH_B: f32 = 0.5;

/// Seconds per (Julian) year — runoff volume/yr → mean discharge (m³/s).
const SECONDS_PER_YEAR: f32 = 3.155_76e7;

/// Convert an accumulated runoff volume in `mm·km²/yr` (the unit of
/// `runoff_accumulation`) to a mean discharge in m³/s: 1 mm over 1 km² = 1000 m³.
fn runoff_km2_to_m3s(mm_km2_yr: f32) -> f32 {
    mm_km2_yr * 1000.0 / SECONDS_PER_YEAR
}

/// ADR 0001 Finding 20 (DEFECT A) — clip the river network against the FINAL lake
/// surfaces so every watercourse terminates at its SINK instead of running through
/// it. The routing field (breached, monotone) drains every basin, so `extract_rivers`
/// traces straight across lakes and out of endorheic sinks; meanwhile `lake_map`
/// marks those basins as standing water. The visible result is rivers crossing lake
/// polygons and orphan reaches emerging below them.
///
/// Each segment is split into its maximal runs of consecutive NON-lake points; each
/// run of ≥2 points becomes a segment (points + long profile sliced exactly). A run
/// that starts at a lake shore is that lake's OUTLET: kept (with the parent discharge)
/// for an EXORHEIC lake, dropped for an ENDORHEIC one (the water dies in the closed
/// basin — no downstream watercourse). Links are remapped through parent head/tail so
/// a segment's `downstream` is `None` exactly when it ends at a sink (sea or lake).
///
/// Touches ONLY the exported/rendered polylines and their parallel arrays — the
/// routing field, lakes, and `lake_map` are untouched. A display width/threshold is a
/// rendering attribute; topology stays continuous to each sink.
pub fn clip_rivers_to_lakes(dr: &mut C1DrainageResult) {
    let w = dr.width;
    let lake_map = &dr.lake_map;
    let endorheic: std::collections::HashSet<u32> = dr
        .lakes
        .iter()
        .filter(|lk| lk.lake_type == LakeType::Endorheic)
        .map(|lk| lk.base.id)
        .collect();
    let in_lake = |&(x, y): &(u32, u32)| lake_map[y as usize * w + x as usize];

    let src = std::mem::take(&mut dr.rivers.segments);
    let src_km2 = std::mem::take(&mut dr.segment_drainage_km2);
    let src_nav = std::mem::take(&mut dr.segment_navigability);
    let src_q = std::mem::take(&mut dr.segment_discharge_m3s);
    let src_width = std::mem::take(&mut dr.segment_width_m);
    let src_prof = std::mem::take(&mut dr.segment_profile_m);

    // Pass 1 — enumerate each parent's kept runs as (start, end) index ranges into
    // its point list, so head/tail new-indices can be reserved before emission.
    let runs: Vec<Vec<(usize, usize)>> = src
        .iter()
        .map(|s| {
            let mut out = Vec::new();
            let (mut a, n) = (0usize, s.points.len());
            while a < n {
                if in_lake(&s.points[a]) != 0 {
                    a += 1;
                    continue;
                }
                let mut b = a;
                while b < n && in_lake(&s.points[b]) == 0 {
                    b += 1;
                }
                // Drop an outlet run (starts just after a lake cell) when the closed basin
                // has no river outflow: an ENDORHEIC lake, OR any BELOW-SEA basin (id ≥
                // 1_000_001, Finding 31) — the latter's outlet is its traced SPILLWAY (the
                // physically-correct least-sill path), not the arbitrary breached-carve reach,
                // so keeping the carve run here would double-count the outflow.
                let is_outlet = a > 0;
                let src_lake = if is_outlet { in_lake(&s.points[a - 1]) } else { 0 };
                let drop = is_outlet && (endorheic.contains(&src_lake) || src_lake >= 1_000_001);
                if b - a >= 2 && !drop {
                    out.push((a, b));
                }
                a = b + 1;
            }
            out
        })
        .collect();

    // Reserve new indices: parent i's runs occupy a contiguous block.
    let mut offset = vec![0usize; src.len() + 1];
    for i in 0..src.len() {
        offset[i + 1] = offset[i] + runs[i].len();
    }
    let head = |i: usize| runs[i].first().map(|_| offset[i]);
    let tail = |i: usize| runs[i].last().map(|_| offset[i] + runs[i].len() - 1);

    let mut segments = Vec::with_capacity(offset[src.len()]);
    let mut km2 = Vec::with_capacity(offset[src.len()]);
    let mut nav = Vec::with_capacity(offset[src.len()]);
    let mut discharge = Vec::with_capacity(offset[src.len()]);
    let mut width_m = Vec::with_capacity(offset[src.len()]);
    let mut profile_m = Vec::with_capacity(offset[src.len()]);

    for (i, s) in src.iter().enumerate() {
        let n = s.points.len();
        for (ri, &(a, b)) in runs[i].iter().enumerate() {
            let is_first = a == 0; // run holds the parent's true source
            let reaches_end = b == n; // run reaches the parent's downstream end
            // Endorheic outlet runs were dropped; any surviving run that starts at a
            // lake shore is an exorheic outlet → inherit the parent discharge.
            let q_km2 = src_km2.get(i).copied().unwrap_or(0.0);
            let upstream = if is_first {
                s.upstream.iter().filter_map(|&u| tail(u)).collect()
            } else {
                Vec::new()
            };
            let downstream = if reaches_end && ri == runs[i].len() - 1 {
                s.downstream.and_then(head)
            } else {
                None // ends at a lake shore (its sink)
            };
            segments.push(RiverSegment {
                points: s.points[a..b].to_vec(),
                strahler_order: s.strahler_order,
                avg_flow: s.avg_flow,
                max_flow: s.max_flow,
                basin_id: s.basin_id,
                upstream,
                downstream,
            });
            km2.push(q_km2);
            nav.push(src_nav.get(i).copied().unwrap_or(Navigability::NonNavigable));
            // Inherit the parent's DISCHARGE + discharge-based width (Finding 22) — an
            // exorheic outlet reach keeps the full upstream catchment's width; recomputing
            // from the reach's local area would drop it across the lake (the author's bug).
            discharge.push(src_q.get(i).copied().unwrap_or(0.0));
            width_m.push(src_width.get(i).copied().unwrap_or(0.0));
            profile_m.push(src_prof.get(i).map(|p| p[a..b].to_vec()).unwrap_or_default());
        }
    }

    dr.rivers.segments = segments;
    dr.segment_drainage_km2 = km2;
    dr.segment_navigability = nav;
    dr.segment_discharge_m3s = discharge;
    dr.segment_width_m = width_m;
    dr.segment_profile_m = profile_m;
}

/// Finding 37 — the exorheic-outlet invariant over the WHOLE lake population (not a subset). EVERY
/// lake marked `Exorheic`, of EVERY provenance, must have a TRACED OUTLET: a river segment whose
/// source (first point) borders the lake footprint — a detect-lake's overflow reach, OR a below-sea
/// basin's appended spillway (which, after `clip_rivers_to_lakes`, starts just outside the pool).
/// Returns the ids of exorheic lakes with NO such outlet. A non-empty result means a regime was
/// labelled exorheic without an emitted outflow: the label is wrong (the lake should be endorheic).
/// Must be called AFTER spillways are appended and `clip_rivers_to_lakes` has run. The earlier guard
/// (`exorheic_below_sea_basin_has_traced_spillway`) checked ONLY below-sea basins on a synthetic
/// grid; the 21 exorheic below-sea lakes shipped WITHOUT an outlet in the 8192² export slipped
/// through that subset (same blind-spot pattern as Finding 36's config subset — see ADR Finding 37).
pub fn exorheic_lakes_missing_outlet(dr: &C1DrainageResult) -> Vec<u32> {
    let (w, h) = (dr.width, dr.height);
    let sources: Vec<(u32, u32)> =
        dr.rivers.segments.iter().filter_map(|s| s.points.first().copied()).collect();
    let borders = |sx: u32, sy: u32, id: u32| -> bool {
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (sx as i32 + dx, sy as i32 + dy);
                if nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < h
                    && dr.lake_map[ny as usize * w + nx as usize] == id
                {
                    return true;
                }
            }
        }
        false
    };
    dr.lakes
        .iter()
        .filter(|lk| lk.lake_type == LakeType::Exorheic)
        .map(|lk| lk.base.id)
        .filter(|&id| !sources.iter().any(|&(sx, sy)| borders(sx, sy, id)))
        .collect()
}

/// ADR 0001 Finding 24 — GEOGRAPHIC SCALE RATIO (hydrology only). A pure PRESENTATION
/// multiplier: the map DRAWS `real_km` of terrain but SIGNIFIES `real_km · ratio` (a
/// Skyrim-style compression). It is NOT a physical quantity and touches NOTHING that
/// shapes the terrain — it is applied here, on the already-computed drainage, to the
/// EXPORT-DERIVED quantities ONLY:
///   - effective catchment  = real catchment · ratio²   (area scales as length²)
///   - discharge Q          = real discharge · ratio²   (Q = runoff · catchment)
///   - channel width         = a·Q^b  ⇒  scales as ratio (√ of ratio²)
///   - navigability class    re-evaluated on the effective catchment
/// It does NOT touch: the routing field, rivers geometry, lakes, `lake_map`, and —
/// crucially — NOT stream-power incision, the lake water balance, precipitation,
/// temperature, or biomes (those run on the REAL 400 km quantities upstream, so every
/// prior calibration holds). Being a post-process on derived arrays, it is instant and
/// does NOT enter the drainage cache key. `ratio == 1.0` → identity (early return).
pub fn apply_geo_scale_ratio(
    dr: &mut C1DrainageResult,
    ratio: f32,
    thresholds: &DrainageThresholds,
) {
    if ratio == 1.0 {
        return;
    }
    let area_scale = ratio * ratio;
    for i in 0..dr.rivers.segments.len() {
        let km2 = dr.segment_drainage_km2[i] * area_scale;
        let q = dr.segment_discharge_m3s[i] * area_scale;
        dr.segment_drainage_km2[i] = km2;
        dr.segment_discharge_m3s[i] = q;
        dr.segment_width_m[i] = CHANNEL_WIDTH_A * q.max(0.0).powf(CHANNEL_WIDTH_B);
        dr.segment_navigability[i] = thresholds.classify(km2);
    }
}

/// #drainage fix B — reference runoff depth (mm/yr) that maps a river's real
/// DISCHARGE (`runoff_accumulation`, in mm·km²/yr) to an effective drainage area
/// in km², so the existing km² navigability thresholds apply unchanged: a humid
/// river at this runoff depth classifies exactly as the old cell-count area did,
/// a drier one downgrades, a dry/endorheic-below one → NonNavigable. Anchored on
/// a typical humid annual runoff depth (~300 mm).
const REFERENCE_RUNOFF_MM: f32 = 300.0;

/// #drainage — runoff (`max(0, precip − PE)·cell_km2`, mm·km²/yr) accumulated
/// downstream along the D8 flow (decreasing `filled` order). `sinks` cells (e.g.
/// endorheic lake cells) do NOT propagate their accumulation downstream — the
/// water dies there (evaporation), so no phantom discharge continues to the sea.
/// `None` sinks → plain accumulation (used for the lake-inflow classification).
fn runoff_accumulation(
    heightmap: &GridF32,
    flow: &FlowResult,
    climate: &DrainageClimate,
    cell_km2: f32,
    sinks: Option<&[bool]>,
    w: usize,
    h: usize,
) -> Vec<f32> {
    let n = w * h;
    let mut acc = vec![0.0f32; n];
    for k in 0..n {
        if heightmap.data[k] > C1_SEA_LEVEL_NORM {
            let p = precip_mm_per_year(climate.precip_internal.data[k]);
            let pe = potential_evaporation_mm(climate.temperature.data[k]);
            acc[k] = (p - pe).max(0.0) * cell_km2;
        }
    }
    let mut order: Vec<usize> = (0..n).filter(|&k| heightmap.data[k] > C1_SEA_LEVEL_NORM).collect();
    order.sort_unstable_by(|&a, &b| {
        flow.filled.data[b].partial_cmp(&flow.filled.data[a]).unwrap_or(std::cmp::Ordering::Equal)
    });
    for &k in &order {
        if sinks.is_some_and(|s| s[k]) {
            continue; // endorheic sink: water dies here, no downstream discharge
        }
        let d = flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let (i, j) = (k % w, k / w);
        let ni = ((i as i32 + D8_DX[d as usize]).rem_euclid(w as i32)) as usize;
        let nj = ((j as i32 + D8_DY[d as usize]).rem_euclid(h as i32)) as usize;
        acc[nj * w + ni] += acc[k];
    }
    acc
}

/// #drainage — the WATER BALANCE: a basin overflows (exorheic) only if its
/// catchment INFLOW exceeds the lake's EVAPORATION at the sill; otherwise it is
/// ENDORHEIC (terminal), stabilising at the level where evaporation = inflow.
/// Mutates `lake_map` (cells above an endorheic lake's equilibrium are drained →
/// `0`) and returns the reclassified / resized lakes.
///
/// Inflow = runoff (`max(0, precip − PE)·cell_km2`) accumulated downstream to the
/// lake. Endorheic equilibrium SURFACE `A_eq = inflow / PE_lake` (evap balances
/// inflow); the level is read off the lake's own hypsometry (its cells sorted by
/// elevation). The surface is on BOTH sides (evap on surface, surface from level)
/// but `A_eq` is direct (no iteration), then mapped to a level — monotone.
#[allow(clippy::too_many_arguments)]
fn water_balance_lakes(
    heightmap: &GridF32,
    flow: &FlowResult,
    climate: &DrainageClimate,
    cell_km2: f32,
    ss: &SteinSteinParams,
    lakes: &[Lake],
    lake_map: &mut [u32],
    w: usize,
    h: usize,
) -> Vec<C1Lake> {
    let n = w * h;
    // 1-2. Runoff (max(0, precip − PE)·cell_km2) accumulated downstream. No sinks
    //      here: we want the full catchment inflow REACHING each candidate lake.
    let runoff_accum = runoff_accumulation(heightmap, flow, climate, cell_km2, None, w, h);
    // 3. Per lake: inflow vs evaporation → exorheic (overflow) or endorheic level.
    let mut out = Vec::with_capacity(lakes.len());
    for lk in lakes {
        let mut cells: Vec<(usize, f32)> =
            (0..n).filter(|&k| lake_map[k] == lk.id).map(|k| (k, heightmap.data[k])).collect();
        if cells.is_empty() {
            continue;
        }
        // Sort by elevation ascending → the lake's hypsometry (floor first).
        cells.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let inflow = cells.iter().map(|&(k, _)| runoff_accum[k]).fold(0.0f32, f32::max);
        let pe_lake = potential_evaporation_mm(climate.temperature.data[cells[0].0]).max(1.0);
        let a_eq_km2 = inflow / pe_lake; // surface where evaporation = inflow
        let a_sill_km2 = lk.area as f32 * cell_km2;

        if a_eq_km2 >= a_sill_km2 {
            // Inflow fills past the sill → OVERFLOW (exorheic), full sill lake.
            let level_m = c1_altitude_norm_to_metres(lk.surface_elevation, ss);
            let floor_m = c1_altitude_norm_to_metres(lk.surface_elevation - lk.max_depth, ss);
            out.push(C1Lake {
                base: lk.clone(),
                level_m,
                depth_m: level_m - floor_m,
                area_km2: a_sill_km2,
                lake_type: LakeType::Exorheic,
            });
        } else {
            // ENDORHEIC: shrink to A_eq; drain the cells above the equilibrium.
            let n_eq = (a_eq_km2 / cell_km2).floor() as usize;
            for &(k, _) in cells.iter().skip(n_eq.max(1)) {
                lake_map[k] = 0; // drained back to terrain (no longer flooded)
            }
            if n_eq == 0 {
                // Inflow can't sustain even one cell → the basin dries up.
                if let Some(&(k0, _)) = cells.first() {
                    lake_map[k0] = 0;
                }
                continue;
            }
            let level = cells[n_eq - 1].1; // equilibrium surface elevation
            let level_m = c1_altitude_norm_to_metres(level, ss);
            let floor_m = c1_altitude_norm_to_metres(cells[0].1, ss);
            let mut base = lk.clone();
            base.surface_elevation = level;
            base.max_depth = level - cells[0].1;
            base.area = n_eq;
            out.push(C1Lake {
                base,
                level_m,
                depth_m: level_m - floor_m,
                area_km2: n_eq as f32 * cell_km2,
                lake_type: LakeType::Endorheic,
            });
        }
    }
    out
}

/// BELOW-SEA INLAND basins as water bodies (ADR 0001 Finding 18). `pit_fill`/`detect_lakes`
/// treat every below-sea cell as ocean, so an enclosed below-sea basin (a `water_class`
/// class-2 component) never enters the lake path. This finds those basins (≥ `lake_min_area`),
/// runs the SAME water balance as [`water_balance_lakes`] — but reading catchment INFLOW at the
/// basin's LAND INLETS (runoff_accumulation zeroes the below-sea basin cells themselves) — and
/// returns the resulting typed lakes + a lake_map (ids offset ≥ 1_000_001 to avoid colliding
/// with detected-lake ids) marking the WATER cells at the level `min(spill, evaporative)`. Cells
/// above that level stay dry land (possibly below sea) — NOT flooded to 0 m.
/// Depth below which a below-sea basin cell reads as WETLAND rather than open water/lagoon
/// (ADR 0001 Finding 30). Measured: these through-flow basins are 100 % shallow (< ~3 m) with
/// a ~0.1–0.3 % shore slope — the wetland signature; deeper cells are a lagoon/inland sea.
const WETLAND_MAX_DEPTH_M: f32 = 3.0;

/// Finding 33 PART A — minimum below-sea basin size (in CELLS, resolution-independent) to enter
/// the exported inventory (`lakes.json`). A few cells: reject single/double-cell noise, keep every
/// visible lake so no river terminates in a body absent from the export.
const INVENTORY_MIN_CELLS: usize = 4;

/// A traced OVERFLOW path for an exorheic below-sea basin (ADR 0001 Finding 30). Mass balance:
/// a basin receiving more than it evaporates MUST overflow, so an `Exorheic` label REQUIRES an
/// outlet — this is the surplus (`inflow − evaporation`) routed downhill from the basin to its
/// sink, emitted into the river network like any other reach. Ready-to-append parallel data.
#[derive(Clone)]
pub struct Spillway {
    pub lake_id: u32,
    pub points: Vec<(u32, u32)>,
    pub discharge_m3s: f32,
    pub drainage_km2: f32,
    pub width_m: f32,
    pub navigability: Navigability,
    pub profile_m: Vec<f32>,
    /// `None` = the path reaches the SEA; `Some(id)` = it CHAINS into another below-sea basin.
    pub chained_into: Option<u32>,
}

/// Per-basin water-balance summary for EVERY enclosed below-sea basin (Finding 32) — inventory
/// AND sub-threshold — so the invariant can be checked over all of them, not just the exported
/// ones. Real units (the water balance is on the true quantities, before any geographic ratio).
#[derive(Clone, Copy)]
pub struct BasinSummary {
    pub id: u32,
    pub area_km2: f32,
    pub exorheic: bool,
    pub inflow_m3s: f32,
    pub evaporation_m3s: f32,
    /// The two numbers the regime compared (in km²): equilibrium area vs sill area.
    pub a_eq_km2: f32,
    pub a_spill_km2: f32,
    pub spill_level_m: f32,
    pub floor_m: f32,
    /// The RETAINED water level (metres) — must equal the sill for an exorheic basin.
    pub level_m: f32,
    /// Max depth (metres) at the retained level — was ~0 when the level collapsed to the floor.
    pub max_depth_m: f32,
    /// Flooded area (km²) AT THE SILL (`a_spill`) — what an exorheic basin should cover.
    pub area_at_sill_km2: f32,
    /// Regime under the OLD inlet reading (`max` runoff at land neighbours) — for the before/after
    /// that separates the inlet bug (Finding 34) from `a_spill` growth (Finding 33).
    pub exorheic_before_inlet_fix: bool,
    /// Finding 37 POINT 2 — the PREDICTION (would overflow: a sill exists AND `a_eq ≥ a_spill`), before
    /// the trace confirms it. `predicted_exorheic && !exorheic` = a basin the balance said would spill
    /// but for which NO outlet could be traced → correctly demoted to endorheic by the inversion.
    pub predicted_exorheic: bool,
    /// Finding 37 TASK 1 — whether a real LOCAL SILL was found (an escape saddle within the window).
    /// `false` = no escape → the basin is endorheic by ABSENCE (no fabricated sea-level fallback).
    pub has_sill: bool,
}

/// Result of [`below_sea_basin_lakes`]: the typed lakes + their water `lake_map`, the traced
/// spill paths (Finding 30 — every exorheic basin gets one), a wetland mask (shallow through-flow
/// margins → `Biome::Wetland`), and a per-basin balance summary (Finding 32, every basin).
pub struct BelowSeaResult {
    pub lakes: Vec<C1Lake>,
    pub lake_map: Vec<u32>,
    pub spillways: Vec<Spillway>,
    pub wetland: Vec<u8>,
    pub basins: Vec<BasinSummary>,
}

pub fn below_sea_basin_lakes(
    heightmap: &GridF32,
    climate: &DrainageClimate,
    cfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    window_km: f32,
    // Finding 39 — the DETECTED lakes' lake_map (small ids), so a below-sea spillway CHAINS into the
    // FIRST basin it enters — detected OR below-sea — instead of threading UNDER it. `None` = only
    // below-sea lakes are visible to the trace (byte-identical to the pre-fix behaviour).
    detected_lake_map: Option<&[u32]>,
) -> BelowSeaResult {
    use crate::lakes::connectivity::water_class;
    use std::collections::VecDeque;
    let (w, h) = (heightmap.width, heightmap.height);
    let n = w * h;
    let cell_km2 = (window_km / w as f32).powi(2);
    let flow =
        compute_flow(heightmap, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, ..Default::default() });
    let wc = water_class(heightmap, C1_SEA_LEVEL_NORM);
    let runoff = runoff_accumulation(heightmap, &flow, climate, cell_km2, None, w, h);
    // Finding 33 PART A — the below-sea INVENTORY floor is now a few CELLS (reject single-cell
    // noise only), NOT 5 km². At 40 m/cell, 5 km² is 3125 cells — a plainly visible lake that a
    // river can terminate in; excluding it from lakes.json made those terminations land in a body
    // absent from the export (invisible to the consumer). The old fear (erosion-fabricated
    // parasitic pits) is handled by the breach conditioning, so it no longer applies.
    let min_cells = INVENTORY_MIN_CELLS;
    let _ = cfg.lake_min_area_km2; // (kept for the detected-lake path; below-sea uses the cell floor)
    let nb4 = |x: i32, y: i32| -> Option<usize> {
        if x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h {
            Some(y as usize * w + x as usize)
        } else {
            None
        }
    };
    let mut lake_map = vec![0u32; n];
    let mut out = Vec::new();
    let mut spillways: Vec<Spillway> = Vec::new();
    let mut basins: Vec<BasinSummary> = Vec::new();
    let mut wetland = vec![0u8; n];
    let mut seen = vec![false; n];
    let mut next_id = 1_000_001u32;
    let metres_per_norm = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);

    // Finding 37c — the ocean-minimax `spill_receiver` (Finding 32) is GONE: it routed spillways by
    // least-max-elevation to the ocean, which for a large basin threads UNDER other lakes / higher
    // ground and gives a non-monotone profile (#1's path under the 106 m lake #4). Spillways now
    // follow the DOWNHILL flow field from the escape saddle (below), which cannot climb.
    let quant = |e: f32| (e * 1_000_000.0) as i32;

    // Finding 36 — a below-sea LAKE fills to its LOCAL sill (the lowest saddle on the rim of the
    // enclosed hollow), NOT the ocean-minimax barrier `barrier_q`. Using `barrier > height` as the
    // pool test and `min barrier_q` as the level floods everything under the CONTINENTAL pass: for
    // the author's #1000020 the real hollow is 17 km² but the lake rose to the 613 m ocean pass and
    // drowned 267 km² of green (94% of its footprint). The genuine hollow is the connected ENCLOSED
    // below-sea component (`water_class == 2`). A SINGLE priority-flood outward from its floor finds,
    // in one pass, both the local sill (the first rim cell from which water can descend to a DIFFERENT
    // sink) and the bounded bowl (cells ≤ sill connected to the floor, in height order). The lake then
    // fills to min(local_sill, evaporative). `barrier_q`/`spill_receiver` are kept only to TRACE the
    // outlet path (Finding 32), never to set the level.
    for s in 0..n {
        if seen[s] || wc[s] != 2 {
            continue;
        }
        // 1. the PIT = the connected enclosed-below-sea component (8-conn) — the real depression.
        let mut comp = Vec::new();
        let mut q = VecDeque::new();
        q.push_back(s);
        seen[s] = true;
        while let Some(k) = q.pop_front() {
            comp.push(k);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            // 8-connected (Finding 35): diagonally-touching enclosed cells are one water body.
            for (dx, dy) in D8_DX.iter().zip(D8_DY.iter()) {
                if let Some(nk) = nb4(x + dx, y + dy) {
                    if wc[nk] == 2 && !seen[nk] {
                        seen[nk] = true;
                        q.push_back(nk);
                    }
                }
            }
        }
        let floor_k = *comp
            .iter()
            .min_by(|&&a, &&b| heightmap.data[a].partial_cmp(&heightmap.data[b]).unwrap())
            .unwrap();
        if heightmap.data[floor_k] >= C1_SEA_LEVEL_NORM {
            continue; // (defensive — wc == 2 already implies a below-sea floor)
        }
        // Finding 34 — TOTAL INFLOW from the runoff FIELD, SUMMED at the shoreline: `runoff_accumulation`
        // ZEROES below-sea cells, so a tributary's discharge is lost the instant it enters the water.
        // SUM the accumulated runoff of every above-sea cell that drains INTO a hollow cell — each
        // tributary counted once, at the shore, before the zeroing.
        let mut inflow = 0.0f32;
        for &k in &comp {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in D8_DX.iter().zip(D8_DY.iter()) {
                if let Some(nk) = nb4(x + dx, y + dy) {
                    if heightmap.data[nk] > C1_SEA_LEVEL_NORM {
                        let dir = flow.direction[nk];
                        if dir != DIR_NONE {
                            let (tx, ty) = (
                                nk as i32 % w as i32 + D8_DX[dir as usize],
                                nk as i32 / w as i32 + D8_DY[dir as usize],
                            );
                            if nb4(tx, ty) == Some(k) {
                                inflow += runoff[nk]; // this tributary enters the water here
                            }
                        }
                    }
                }
            }
        }
        // 2+3. LOCAL SILL + BOWL in one priority-flood outward from the floor (Barnes 2014). Pop the
        // lowest reachable cell; the FIRST cell with a NOT-YET-FINALISED strictly-lower neighbour is
        // the lowest saddle (water spills there and runs downhill to a different sink). Bowl membership
        // is marked at POP time (finalised), NOT at push — marking at push masked the escape (a saddle's
        // downhill neighbour, pushed earlier by a sibling, read as "seen"), which is exactly what left
        // basins with NO sill: the old `sill_q == MAX` fallback then set the level to sea level (0 m),
        // yielding an "exorheic at 0 m" with no real outlet (Finding 37 POINT 1). `sill_opt` is `None`
        // when no escape exists within the window (a border-truncated basin) → it cannot be exorheic.
        // Also capture the SADDLE and its lowest EXTERIOR escape neighbour — the spillway starts there
        // and runs DOWNHILL from the far side of the divide (Finding 37c), NOT along the ocean-minimax
        // `spill_receiver` (which, for a large basin, threads UNDER other lakes / higher ground and
        // yields a non-monotone profile).
        let (sill_opt, fcells, saddle_opt, escape_opt): (
            Option<f32>,
            Vec<usize>,
            Option<usize>,
            Option<usize>,
        ) = {
            use std::cmp::Reverse;
            use std::collections::{BinaryHeap, HashSet};
            let comp_set: HashSet<usize> = comp.iter().copied().collect();
            let mut in_bowl: HashSet<usize> = HashSet::new();
            let mut queued: HashSet<usize> = HashSet::new();
            let mut heap: BinaryHeap<Reverse<(i32, usize)>> = BinaryHeap::new();
            heap.push(Reverse((quant(heightmap.data[floor_k]), floor_k)));
            queued.insert(floor_k);
            let mut bowl: Vec<usize> = Vec::new();
            let mut sill_q: Option<i32> = None;
            let (mut saddle, mut escape) = (None, None);
            while let Some(Reverse((hq, c))) = heap.pop() {
                let (x, y) = ((c % w) as i32, (c / w) as i32);
                // Escape: the LOWEST neighbour strictly below c that is NOT finalised in the bowl AND is
                // OUTSIDE THIS region's own component (`!comp_set`). Finding 39 refines Finding 38: an
                // internal saddle belongs to the SAME `comp` (another sub-pocket of this region) → it is
                // ABSORBED (the flood fills the whole component before escaping, so no shore sliver is
                // left for an orphan mouth). But a descent to a DIFFERENT below-sea region (a wc==2 cell
                // NOT in this comp) is a real POUR POINT — the lake overflows/chains there at that LOW
                // col, it does NOT climb on to the far continental pass. The old `wc != 2` test absorbed
                // neighbouring regions too, pushing the sill up to the 613 m ocean pass and re-opening
                // the Finding 36 over-flood the instant a humid basin (net_evap 0) filled to it.
                let mut best_e: Option<usize> = None;
                for (dx, dy) in D8_DX.iter().zip(D8_DY.iter()) {
                    if let Some(e) = nb4(x + dx, y + dy) {
                        if quant(heightmap.data[e]) < hq
                            && !in_bowl.contains(&e)
                            && !comp_set.contains(&e)
                            && best_e.map_or(true, |b| heightmap.data[e] < heightmap.data[b])
                        {
                            best_e = Some(e);
                        }
                    }
                }
                if let Some(e) = best_e {
                    sill_q = Some(hq);
                    saddle = Some(c);
                    escape = Some(e);
                    bowl.push(c); // the saddle cell is the pour point (part of the water body)
                    in_bowl.insert(c);
                    break;
                }
                bowl.push(c);
                in_bowl.insert(c);
                for (dx, dy) in D8_DX.iter().zip(D8_DY.iter()) {
                    if let Some(e) = nb4(x + dx, y + dy) {
                        if !queued.contains(&e) {
                            queued.insert(e);
                            heap.push(Reverse((quant(heightmap.data[e]), e)));
                        }
                    }
                }
            }
            (sill_q.map(|q| q as f32 / 1_000_000.0), bowl, saddle, escape)
        };
        // Finding 39 — WATER-BALANCE CLOSURE (replaces the Finding 38b geometric near-sea rule). A
        // below-sea basin OVERFLOWS (exorheic) iff its catchment INFLOW exceeds what the basin can
        // evaporate when filled to its sill; otherwise it is ENDORHEIC, ponding at the level where
        // evaporation balances inflow. DISCHARGE is the discriminant, NOT the sill height: a deep-rimmed
        // basin (#1000022's 471 m col) fed by rivers (961 m³/s signified) in a HUMID climate fills to
        // the col and spills — that high level is now JUSTIFIED BY INFLOW, exactly the case Finding 36's
        // geometric ban wrongly conflated with an unjustified fill. This unifies 36 (arid, low inflow →
        // stays low) and 38b (humid, high inflow → fills & overflows) under one physical law.
        //
        // NO DOUBLE COUNT (precaution 1): `runoff = max(0, precip − PE)` (the inflow source) and
        // `net_evap = max(0, PE − precip)` (the surface loss) are COMPLEMENTARY — for any cell exactly
        // one is nonzero. `runoff_accumulation` sources runoff ONLY from above-sea land and ZEROES the
        // below-sea region, so the lake's own precipitation is never in `inflow`; the surface term
        // charges `net_evap` over the flooded area. A submerged cell therefore contributes to EITHER
        // inflow (humid: net_evap = 0 there) OR surface evaporation (arid: runoff = 0 there), never both.
        let a_spill = fcells.len() as f32 * cell_km2; // area to fill the whole bowl to its sill
        let pe_floor = potential_evaporation_mm(climate.temperature.data[floor_k]);
        let precip_floor = precip_mm_per_year(climate.precip_internal.data[floor_k]);
        let net_evap = (pe_floor - precip_floor).max(0.0); // mm/yr the SURFACE loses net of its own rain
        let pe_lake = pe_floor.max(1.0); // gross PE, kept for the before/after diagnostic below
        // Equilibrium surface where net evaporation balances inflow (km²). net_evap == 0 (a HUMID
        // climate that cannot close the basin) ⇒ infinite equilibrium ⇒ the basin MUST overflow.
        let a_eq = if net_evap > 0.0 { inflow / net_evap } else { f32::INFINITY };
        // Precaution 2 — the reformulated anti-over-flood guard: a HIGH level is admissible ONLY when
        // INFLOW JUSTIFIES it (`a_eq ≥ a_spill`: inflow fills the whole bowl to the sill). The GEOMETRIC
        // invariant is untouched — the footprint is the priority-flood bowl `fcells` (cells ≤ level, all
        // connected to the floor), so `claimed == valid` by construction (no over-flood).
        let fills_to_sill = a_eq >= a_spill;
        let overflow = fills_to_sill; // candidate; confirmed EXORHEIC only if an outlet also traces
        let traced: Option<(Vec<usize>, Option<u32>)> = if overflow {
            // Re-entry into THIS lake's own water body = its flood bowl `fcells` (the loop / case A).
            let in_bowl: std::collections::HashSet<usize> = fcells.iter().copied().collect();
            // Finding 37c — the spillway starts at the SADDLE, crosses to its lowest EXTERIOR escape
            // cell (across the divide), then follows the DOWNHILL flow field to a sink: the OCEAN or a
            // DIFFERENT lake (a chain). Downhill by construction, it never climbs (monotone profile)
            // and never threads under higher ground / other lakes the way the ocean-minimax
            // `spill_receiver` did (that produced #1's path UNDER the 106 m lake #4). It re-enters its
            // own hollow only for a true loop → then INVALID → endorheic.
            match (saddle_opt, escape_opt) {
                (Some(s0), Some(e0)) => {
                    let mut path = vec![s0, e0];
                    let (mut cur, mut chained, mut steps, mut valid) = (e0, None, 0usize, false);
                    loop {
                        if wc[cur] == 1 {
                            valid = true; // reached the OCEAN
                            break;
                        }
                        // THIS lake is not marked yet (marking happens after), so any non-zero
                        // lake_map cell belongs to an ALREADY-PROCESSED, different lake → a valid chain.
                        // Finding 39 — also stop at a DETECTED lake (the `detected_lake_map`): the
                        // spillway must halt at the FIRST basin it reaches, not run UNDER it (river #11
                        // through the 211 m lake #17). Below-sea takes precedence when both are present.
                        let lid = lake_map[cur];
                        let did = detected_lake_map.map_or(0, |d| d[cur]);
                        if lid != 0 || did != 0 {
                            chained = Some(if lid != 0 { lid } else { did });
                            valid = true;
                            break;
                        }
                        if in_bowl.contains(&cur) {
                            break; // re-entered its OWN hollow → a loop → INVALID → endorheic
                        }
                        let d = flow.direction[cur];
                        if d == DIR_NONE {
                            break;
                        }
                        let (cx, cy) = ((cur % w) as i32, (cur / w) as i32);
                        let nx = (cx + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
                        let ny = (cy + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
                        cur = ny * w + nx;
                        path.push(cur);
                        steps += 1;
                        if steps > n {
                            break;
                        }
                    }
                    if valid { Some((path, chained)) } else { None }
                }
                _ => None,
            }
        } else {
            None
        };
        // Finding 39 — regime & surface from the water balance. EXORHEIC iff the inflow fills the bowl
        // to the sill AND a downhill outlet traces to a sink → the surface sits AT the sill (the lake
        // spills there). If it would fill to the sill but no outlet traces (border-truncated / true
        // loop) it is a brim-full ENDORHEIC basin, still at the sill. Otherwise ENDORHEIC BELOW the
        // sill: the surface sits at the level where the bowl area reaches the evaporative equilibrium
        // `a_eq` — read off the priority-flood bowl's HYPSOMETRY (precaution 3): `fcells` sorted by
        // elevation IS the area-vs-level table the flood already swept (floor→sill), so
        // `area(level = sorted[i]) = (i+1)·cell_km2` and the level is `sorted[⌊a_eq/cell⌋−1]`.
        let lake_type = if traced.is_some() { LakeType::Exorheic } else { LakeType::Endorheic };
        let level = if fills_to_sill {
            sill_opt.unwrap_or_else(|| {
                fcells.iter().map(|&k| heightmap.data[k]).fold(f32::MIN, f32::max)
            })
        } else {
            let mut sorted: Vec<f32> = fcells.iter().map(|&k| heightmap.data[k]).collect();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let n_eq = ((a_eq / cell_km2).floor() as usize).clamp(1, sorted.len());
            sorted[n_eq - 1]
        };
        // Finding 39 regression fix — a below-sea basin's SURFACE never reads below sea: the whole
        // enclosed region IS the water body (a would-be-sea depression), so the footprint always covers
        // it up to at least sea level. This restores the Finding 38 guarantee — every below-sea cell is
        // inside the sink, so no river mouth lands on an unmarked shelf (the 3 arid orphan mouths) — and
        // kills the flat-floor depth-0 artefact: an arid basin whose evaporative level collapses onto a
        // flat floor now reads as a sea-level inland sea of depth sea−floor, not a 0-depth sheet marked
        // over the whole floor. A basin whose balance rises ABOVE sea (humid: #1000009 fills to 44.8 m,
        // #1000022 to its col) keeps that higher surface — Finding 39's fill is unchanged there.
        let surface = level.max(C1_SEA_LEVEL_NORM);
        // Invariant 1 (Finding 37 TASK 2) — A LAKE MUST HAVE WATER. A hollow with NO supply (no inflow)
        // AND no computed depth (an endorheic basin whose evaporative level collapses onto its floor)
        // is a DRY DEPRESSION: it belongs in the relief, not in the water bodies. Do not mark it, do
        // not inventory it, do not emit a spillway — its absence of a physical water path is the point.
        // Finding 38 — with `level = sill` the depth is always > 0 (fabricated by the fill), so it can
        // no longer stand in for "has water". A region with NO supply is a DRY salt flat regardless of
        // how deep its rim is: exclude on inflow alone. (This also drops the tiny no-inflow below-sea
        // pockets the author flagged.) A region with any runoff inflow keeps its inland sea.
        let is_dry = inflow <= 0.0;
        if is_dry {
            continue;
        }
        // `spill` for the basin summary: the real local sill if one exists, else the (endorheic) level.
        let spill = sill_opt.unwrap_or(level);
        // 4. mark the WATER cells = the priority-flood BOWL `fcells` up to the SURFACE. Finding 39: the
        // footprint is the whole depression under the surface (the region `comp` PLUS any land ring the
        // surface submerges), so a filled lake COVERS everything below its surface — the montane lake,
        // not a vertical-walled sliver. `surface = level.max(sea)` also guarantees the ENTIRE below-sea
        // region is claimed (Finding 38), so no mouth lands on an unmarked shelf. `fcells` runs from the
        // floor up to the sill, so `fcells ≤ surface` is a connected flood from the lowest point ⟹
        // `claimed == valid` (no over-flood). Claim only free cells (`lake_map == 0`) so two adjacent
        // bowls sharing a ridge stay DISJOINT; the escape at the separating land stops a bowl before it
        // can reach another region, so `comp`s never merge.
        let mark_level = surface;
        let id = next_id;
        next_id += 1;
        let mut water = 0usize;
        let mut floor = f32::MAX;
        let (mut ox, mut oy) = (0u32, 0u32);
        for &k in &fcells {
            if heightmap.data[k] <= mark_level && lake_map[k] == 0 {
                lake_map[k] = id;
                water += 1;
                // Finding 30 — wetland vs lagoon by DEPTH: a shallow (< 3 m) water cell is a
                // through-flow wetland margin; a deeper one is open water (lagoon/inland sea).
                if (surface - heightmap.data[k]) * metres_per_norm < WETLAND_MAX_DEPTH_M {
                    wetland[k] = 1;
                }
                if heightmap.data[k] < floor {
                    floor = heightmap.data[k];
                    ox = (k % w) as u32;
                    oy = (k / w) as u32;
                }
            }
        }
        if water == 0 {
            continue;
        }
        // Finding 39 regression fix — report the SURFACE (never below sea), so a below-sea basin's
        // level/depth describe the inland-sea/pan (depth = surface − floor), not a collapsed 0-depth sheet.
        let level_m = c1_altitude_norm_to_metres(surface, ss);
        let floor_m = c1_altitude_norm_to_metres(floor, ss);
        // Inventory (lakes.json / microscope) keeps the 5 km² threshold — no micro-lakes.
        // The basin is still marked + gets a spillway below regardless (sink validity).
        if water >= min_cells {
            out.push(C1Lake {
                base: Lake {
                    id,
                    surface_elevation: surface,
                    max_depth: surface - floor,
                    area: water,
                    basin_id: flow.basins[fcells[0]],
                    outlet: (ox, oy),
                    shallow: (level_m - floor_m) < 10.0,
                },
                level_m,
                depth_m: level_m - floor_m,
                area_km2: water as f32 * cell_km2,
                lake_type,
            });
        }
        // Per-basin balance summary for EVERY basin (Finding 32), in real units.
        let area_km2_bs = water as f32 * cell_km2;
        // OLD inlet reading (max runoff at the pool's land neighbours) → its regime, for before/after.
        let mut inflow_max = 0.0f32;
        for &k in &comp {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1)] {
                if let Some(nk) = nb4(x + dx, y + dy) {
                    if heightmap.data[nk] > C1_SEA_LEVEL_NORM {
                        inflow_max = inflow_max.max(runoff[nk]);
                    }
                }
            }
        }
        basins.push(BasinSummary {
            id,
            area_km2: area_km2_bs,
            exorheic: lake_type == LakeType::Exorheic,
            inflow_m3s: runoff_km2_to_m3s(inflow),
            evaporation_m3s: runoff_km2_to_m3s(net_evap * area_km2_bs),
            a_eq_km2: if a_eq.is_finite() { a_eq } else { a_spill }, // humid ⇒ ∞; report the sill area it fills
            a_spill_km2: a_spill,
            spill_level_m: c1_altitude_norm_to_metres(spill, ss),
            floor_m,
            level_m,
            max_depth_m: level_m - floor_m,
            area_at_sill_km2: a_spill,
            exorheic_before_inlet_fix: (inflow_max / pe_lake) >= a_spill,
            predicted_exorheic: overflow,
            has_sill: sill_opt.is_some(),
        });
        // Finding 37 POINT 2 — the spillway was ALREADY traced (before the regime decision), and its
        // existence is WHY this basin is exorheic. Emit it from the traced path — no second trace, no
        // way for an exorheic label to stand without one.
        if let Some((path, chained)) = &traced {
            let points: Vec<(u32, u32)> =
                path.iter().map(|&k| ((k % w) as u32, (k / w) as u32)).collect();
            let profile: Vec<f32> =
                path.iter().map(|&k| c1_altitude_norm_to_metres(heightmap.data[k], ss)).collect();
            // discharge = surplus (inflow − net evaporation); net_evap = max(0, PE − precip) · area.
            let area_km2 = water as f32 * cell_km2;
            let surplus = (inflow - net_evap * area_km2).max(0.0); // mm·km²/yr
            let discharge_m3s = runoff_km2_to_m3s(surplus);
            let drainage_km2 = surplus / REFERENCE_RUNOFF_MM;
            spillways.push(Spillway {
                lake_id: id,
                points,
                discharge_m3s,
                drainage_km2,
                width_m: CHANNEL_WIDTH_A * discharge_m3s.max(0.0).powf(CHANNEL_WIDTH_B),
                navigability: cfg.thresholds.classify(drainage_km2),
                profile_m: profile,
                chained_into: *chained,
            });
        }
    }
    BelowSeaResult { lakes: out, lake_map, spillways, wetland, basins }
}

/// Trace a lake's outlet downstream via D8; does it reach the sea (≤ 0.5)?
/// (With priority-flood this is always true — the honest computation that would
/// flag `Endorheic` if a future terminal basin ever existed.)
fn outlet_reaches_sea(
    lake: &Lake,
    flow: &FlowResult,
    heightmap: &GridF32,
    w: usize,
    h: usize,
) -> bool {
    let (mut x, mut y) = (lake.outlet.0 as usize, lake.outlet.1 as usize);
    let max_steps = w * h;
    for _ in 0..=max_steps {
        let k = y * w + x;
        if heightmap.data[k] <= C1_SEA_LEVEL_NORM {
            return true;
        }
        let d = flow.direction[k];
        if d == DIR_NONE {
            return false;
        }
        x = ((x as i32 + D8_DX[d as usize]).rem_euclid(w as i32)) as usize;
        y = ((y as i32 + D8_DY[d as usize]).rem_euclid(h as i32)) as usize;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Finding 30 INVARIANT (permanent, non-ignored) — an EXORHEIC below-sea basin MUST have a
    /// traced spillway reaching a sink. Mass balance: a basin that receives more than it
    /// evaporates overflows; the label is meaningless without an outlet. Same spirit as
    /// `clip_rivers_terminate_at_lake_sinks`. If this ever fails, a regime was labelled without
    /// a network behind it (the class of bug this whole thread chased).
    #[test]
    fn exorheic_below_sea_basin_has_traced_spillway() {
        use crate::lakes::connectivity::water_class;
        // A WIDE grid with the ocean at the far left and a below-sea pit at the far RIGHT — so a
        // bounded per-basin search could miss the coast, but the ocean priority-flood cannot.
        // Plus a TINY (sub-inventory) pit, to cover the blind spot: the invariant is checked over
        // EVERY basin via `r.basins`, not just the inventoried (≥ 5 km²) lakes.
        let (w, h) = (80usize, 12usize);
        let mut hm = GridF32::new(w, h, 0.7);
        for y in 0..h {
            hm.set(0, y, 0.2); // ocean column (class-1)
            hm.set(1, y, 0.2);
        }
        // low NEAR-SEA lane at y=6 across the whole width (0.50015 ≈ +1.7 m) so a below-sea sea can
        // legitimately overflow it to the coast (Finding 38 — a high lane would be an enclosed rim →
        // endorheic, no over-flood).
        for x in 2..w {
            hm.set(x, 6, 0.50015);
        }
        for y in 5..8 {
            for x in 70..73 {
                hm.set(x, y, 0.498); // FAR enclosed below-sea pit (≈ −22 m)
            }
        }
        hm.set(40, 6, 0.498); // a TINY 1-cell below-sea pit mid-grid (sub-inventory)
        let precip = GridF32::new(w, h, 1.0);
        let temp = GridF32::new(w, h, 10.0);
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let r = below_sea_basin_lakes(
            &hm,
            &clim,
            &C1DrainageConfig::default(),
            &SteinSteinParams::default(),
            80.0,
            None,
        );
        let wc = water_class(&hm, C1_SEA_LEVEL_NORM);
        // The INVARIANT over EVERY basin (any area): exorheic ⟹ a traced spillway reaching a sink.
        let exo = r.basins.iter().filter(|b| b.exorheic).count();
        assert!(exo > 0, "at least one below-sea basin must classify EXORHEIC");
        for b in r.basins.iter().filter(|b| b.exorheic) {
            let sw = r.spillways.iter().find(|s| s.lake_id == b.id).unwrap_or_else(|| {
                panic!(
                    "EXORHEIC basin #{} ({:.2} km²) MUST have a traced spillway",
                    b.id, b.area_km2
                )
            });
            let &(lx, ly) = sw.points.last().unwrap();
            let end = ly as usize * w + lx as usize;
            assert!(
                wc[end] == 1 || sw.chained_into.is_some(),
                "the spillway must reach a SINK (ocean or another basin), ended at wc={}",
                wc[end]
            );
        }
        // Finding 35 INVARIANTS (permanent):
        // (a) NO TWO LAKE FOOTPRINTS OVERLAP — each `lake_map` cell carries exactly one id, and a
        //     lake's cell-count equals its id's footprint (no lake was silently overwritten).
        let mut count: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for &id in &r.lake_map {
            if id != 0 {
                *count.entry(id).or_default() += 1;
            }
        }
        for lk in &r.lakes {
            assert_eq!(
                count.get(&lk.base.id).copied().unwrap_or(0),
                lk.base.area,
                "lake #{} footprint overlaps another (marked cells != claimed area)",
                lk.base.id
            );
        }
        // (b) A lake's max depth = level − floor, and the level is NEVER below the floor.
        for b in &r.basins {
            assert!(b.level_m >= b.floor_m - 1e-3, "level below floor for basin #{}", b.id);
            assert!(
                (b.max_depth_m - (b.level_m - b.floor_m)).abs() < 1e-3,
                "depth != level−floor for #{}",
                b.id
            );
        }
    }

    /// Finding 36 INVARIANTS (permanent, non-ignored) — a below-sea lake fills to its LOCAL sill, and
    /// its footprint is a real bowl, not "everything under the continental pass". The prior model set
    /// `level = min(barrier_q)` (the ocean-minimax barrier) and flooded every cell under it: for the
    /// author's #1000020 the true hollow was 17 km² but the lake rose to the 613 m ocean pass and
    /// drowned 267 km² of green (94%). Both existing guards were BLIND to it — "depth == level − floor"
    /// is satisfied by 613 − (−20) = 633 (internal consistency, not plausibility), and the overlap
    /// check compared cell SETS which are genuinely disjoint. These two are the missing plausibility
    /// guards. TASK 2: every lake cell is ≤ level AND connected to the floor through cells ≤ level.
    /// TASK 3: a lake's level never exceeds the ARRIVAL altitude of its inlets (water cannot flow uphill).
    #[test]
    fn below_sea_lake_fills_to_local_sill_not_ocean_barrier() {
        use crate::lakes::connectivity::water_class;
        // A deep pit behind a HIGH ocean ridge (0.9) but with a LOW local saddle (0.52). The ocean
        // minimax barrier is 0.9; the local sill is 0.52. The OLD model filled to 0.9 (drowning the
        // whole shelf); the fix fills to 0.52. A right-hand ramp feeds the pit so it overflows.
        let (w, h) = (30usize, 3usize);
        let mut hm = GridF32::new(w, h, 0.7);
        for y in 0..h {
            hm.set(0, y, 0.2); // ocean column (class-1)
            hm.set(1, y, 0.9); // HIGH ridge — the only path to the ocean (barrier_q = 0.9)
        }
        for x in 2..=18 {
            hm.set(x, 1, 0.51); // above-sea shelf below the saddle (the pit spills onto it)
        }
        hm.set(19, 1, 0.52); // the pit's LOCAL sill (lowest rim saddle)
        for (i, x) in (20..=23).enumerate() {
            hm.set(x, 1, 0.30 + 0.02 * i as f32); // enclosed below-sea pit floor (4 cells → inventoried)
        }
        for (i, x) in (24..=28).enumerate() {
            hm.set(x, 1, 0.55 + 0.01 * i as f32); // ramp draining LEFT into the pit (inflow ⟹ overflow)
        }
        let precip = GridF32::new(w, h, 1.0);
        let temp = GridF32::new(w, h, 10.0);
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let ss = SteinSteinParams::default();
        let r = below_sea_basin_lakes(&hm, &clim, &C1DrainageConfig::default(), &ss, 30.0, None);
        let wc = water_class(&hm, C1_SEA_LEVEL_NORM);
        assert!(!r.lakes.is_empty(), "the below-sea pit must be inventoried");
        // REGRESSION: the pit filled to its LOCAL sill (~0.52), NOT the 0.9 ocean barrier. In metres
        // that is the difference between a 2-cell shore pond and a shelf-drowning flood.
        let a = r.lakes.iter().max_by(|x, y| x.base.area.cmp(&y.base.area)).unwrap();
        assert!(
            a.base.surface_elevation < 0.60,
            "lake filled to {:.3} — it rose toward the 0.9 ocean barrier instead of its 0.52 local sill",
            a.base.surface_elevation
        );
        // Recompute flow to read inlet arrival altitudes (TASK 3).
        let flow =
            compute_flow(&hm, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM, ..Default::default() });
        for lk in &r.lakes {
            let id = lk.base.id;
            let level = lk.base.surface_elevation;
            let claimed: Vec<usize> = (0..w * h).filter(|&k| r.lake_map[k] == id).collect();
            // TASK 2a — every claimed cell is at or below the lake level.
            for &k in &claimed {
                assert!(
                    hm.data[k] <= level + 1e-4,
                    "lake #{id}: cell {k} at {:.3} is ABOVE the lake level {:.3}",
                    hm.data[k],
                    level
                );
            }
            // TASK 2b — every claimed cell is connected to the floor through cells also ≤ level.
            let floor = *claimed
                .iter()
                .min_by(|&&x, &&y| hm.data[x].partial_cmp(&hm.data[y]).unwrap())
                .unwrap();
            let mut inset = vec![false; w * h];
            for &k in &claimed {
                inset[k] = true;
            }
            let mut reach = vec![false; w * h];
            let mut q = std::collections::VecDeque::new();
            q.push_back(floor);
            reach[floor] = true;
            while let Some(k) = q.pop_front() {
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                for (dx, dy) in D8_DX.iter().zip(D8_DY.iter()) {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if inset[nk] && !reach[nk] {
                            reach[nk] = true;
                            q.push_back(nk);
                        }
                    }
                }
            }
            let disconnected = claimed.iter().filter(|&&k| !reach[k]).count();
            assert_eq!(
                disconnected, 0,
                "lake #{id}: {disconnected} cells ≤level but DISCONNECTED from the floor (disjoint puddles)"
            );
            // TASK 3 — the level never exceeds the ARRIVAL altitude of any inlet (no uphill flow).
            let mut min_inlet = f32::MAX;
            for k in 0..w * h {
                if hm.data[k] <= C1_SEA_LEVEL_NORM {
                    continue; // inlets are above-sea land cells
                }
                let dir = flow.direction[k];
                if dir == DIR_NONE {
                    continue;
                }
                let (tx, ty) = (
                    k as i32 % w as i32 + D8_DX[dir as usize],
                    k as i32 / w as i32 + D8_DY[dir as usize],
                );
                if tx >= 0 && ty >= 0 && (tx as usize) < w && (ty as usize) < h {
                    let tk = ty as usize * w + tx as usize;
                    if r.lake_map[tk] == id {
                        min_inlet = min_inlet.min(hm.data[k]);
                    }
                }
            }
            if min_inlet != f32::MAX {
                assert!(
                    level <= min_inlet + 1e-4,
                    "lake #{id}: level {:.3} EXCEEDS its lowest inlet arrival {:.3} — rivers would flow uphill",
                    level,
                    min_inlet
                );
            }
        }
        let _ = &wc;
    }

    /// Finding 37 TASK 2 INVARIANTS (permanent, non-ignored) — the three that close the class:
    /// Finding 39 INVARIANT (permanent) — DISCHARGE IS THE DISCRIMINANT for a deeply-enclosed below-sea
    /// basin (its overflow col far above sea), NOT the col height. The SAME geometry gives OPPOSITE
    /// regimes under two climates, which is the whole point of the water-balance closure and the
    /// resolution of the Finding 36 ↔ 38b conflict:
    ///   • HUMID (precip ≫ PE ⟹ net_evap = 0): inflow cannot be evaporated, so the basin fills to its
    ///     56 m col and OVERFLOWS (exorheic, spillway) — the high level is JUSTIFIED by inflow (this is
    ///     the #1000022 case; Finding 38b wrongly forced it endorheic-at-sea).
    ///   • ARID SINK (a cool humid catchment feeding a hot evaporating floor, PE_floor ≫ precip):
    ///     the surface evaporates the inflow well below the col ⟹ ENDORHEIC, level near the floor (the
    ///     Caspian/Finding 36 case — preserved, because inflow no longer justifies the fill).
    /// Either way the footprint is the priority-flood bowl (every marked cell ≤ level), so the geometric
    /// invariant (claimed == valid, no over-flood) holds regardless of regime.
    fn enclosed_basin_geometry() -> (usize, usize, GridF32) {
        let (w, h) = (24usize, 24usize);
        let mut hm = GridF32::new(w, h, 0.51); // land at ~113 m everywhere
        for y in 0..h {
            hm.set(0, y, 0.2); // ocean column (class-1) at the far left
        }
        // A below-sea pit at the centre (norm 0.498 ≈ −22 m), enclosed by the 0.51 land. Its ONLY rim
        // gap is a HIGH col at (13,11) = 0.505 (≈ 56 m) leading to a slightly lower exterior — the
        // flood's escape (the overflow col) is 56 m up, nothing near sea.
        for x in 11..=13 {
            for y in 11..=13 {
                hm.set(x, y, 0.498);
            }
        }
        hm.set(13, 11, 0.505); // the col (high)
        // A descending lane from the col to the ocean so an exorheic outlet can trace downhill to sea.
        for x in 1..14 {
            hm.set(x, 11, 0.503 - (14 - x) as f32 * 0.0001);
        }
        for x in 11..=13 {
            for y in 12..=13 {
                hm.set(x, y, 0.498); // restore the pit body (lane touched only row 11)
            }
        }
        hm.set(11, 11, 0.498);
        hm.set(12, 11, 0.498);
        hm.set(13, 11, 0.505);
        (w, h, hm)
    }

    #[test]
    fn humid_enclosed_below_sea_fills_to_col_and_overflows() {
        let (w, h, hm) = enclosed_basin_geometry();
        let precip = GridF32::new(w, h, 1.0); // 16500 mm/yr — precip ≫ PE ⟹ net_evap 0
        let temp = GridF32::new(w, h, 10.0);
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let ss = SteinSteinParams::default();
        let r = below_sea_basin_lakes(&hm, &clim, &C1DrainageConfig::default(), &ss, 24.0, None);
        assert!(!r.lakes.is_empty(), "the enclosed pit must be a below-sea lake");
        let lk = &r.lakes[0];
        assert_eq!(
            lk.lake_type,
            LakeType::Exorheic,
            "a HUMID enclosed basin (net_evap 0) must fill to its col and overflow — inflow justifies it"
        );
        assert!(
            lk.level_m > 20.0,
            "the surface must sit at the ~56 m col (justified by inflow), got {:.0} m",
            lk.level_m
        );
        assert!(
            r.spillways.iter().any(|s| s.lake_id == lk.base.id),
            "an exorheic basin MUST emit a traced spillway"
        );
        // Geometric invariant survives the high level — every marked cell ≤ the level.
        for k in 0..w * h {
            if r.lake_map[k] == lk.base.id {
                assert!(
                    hm.data[k] <= lk.base.surface_elevation + 1e-4,
                    "a marked cell is above the lake surface"
                );
            }
        }
    }

    #[test]
    fn arid_sink_enclosed_below_sea_stays_endorheic_below_col() {
        use crate::lakes::connectivity::water_class;
        let (w, h, hm) = enclosed_basin_geometry();
        // Cool humid catchment (generates inflow) feeding a HOT evaporating floor (closes the basin
        // well below the col): PE_floor(30 °C) ≈ 2600 mm/yr ≫ precip 800 mm/yr ⟹ net_evap > 0.
        let precip = GridF32::new(w, h, 800.0 / crate::climate::precipitation::PRECIP_MM_PER_UNIT);
        let mut temp = GridF32::new(w, h, 5.0); // cool catchment → runoff = 800 − PE(5) > 0
        for x in 11..=13 {
            for y in 11..=13 {
                temp.set(x, y, 30.0); // hot pit surface → high PE, evaporates the inflow
            }
        }
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let ss = SteinSteinParams::default();
        let r = below_sea_basin_lakes(&hm, &clim, &C1DrainageConfig::default(), &ss, 24.0, None);
        assert!(!r.lakes.is_empty(), "the fed pit must still be a (terminal) below-sea lake");
        let lk = &r.lakes[0];
        assert_eq!(
            lk.lake_type,
            LakeType::Endorheic,
            "an ARID-sink enclosed basin evaporates its inflow below the col → endorheic"
        );
        assert!(
            lk.level_m < 20.0,
            "the surface must stay near the floor (below the 56 m col), got {:.0} m",
            lk.level_m
        );
        // Finding 39 regression fix — the below-sea region is a would-be-sea depression, so it reads as
        // an inland sea at sea level: depth = sea − floor > 0 (NOT the collapsed 0-depth sheet the
        // evaporative level gave on a flat floor), and the WHOLE below-sea region is marked (no
        // unmarked shelf for a river mouth to strand on). Both would fail before the fix.
        assert!(
            lk.depth_m > 5.0,
            "an arid below-sea basin must read as a sea-level inland sea (depth = sea−floor), not a \
             0-depth sheet; got depth {:.1} m",
            lk.depth_m
        );
        let wc = water_class(&hm, C1_SEA_LEVEL_NORM);
        let unmarked_below_sea =
            (0..w * h).filter(|&k| wc[k] == 2 && r.lake_map[k] == 0).count();
        assert_eq!(
            unmarked_below_sea, 0,
            "every below-sea cell of the region must be inside the sink (no orphan shelf); {unmarked_below_sea} unmarked"
        );
    }

    /// Finding 38 INVARIANT (permanent) — (1) a LAKE MUST HAVE WATER (a dry depression is not
    /// inventoried nor marked); (2) an OUTLET
    /// ARRIVES ELSEWHERE (its sink is the ocean or a DIFFERENT lake, never its own footprint);
    /// (3) an OUTLET MAY NOT RUN THROUGH A RETAINING POCKET (no interior spillway cell below sea
    /// level that is not the ocean). A deep below-sea pit fed by a ramp overflows a low sill into an
    /// arm of the sea — the outlet reaches `water_class == 1` directly.
    #[test]
    fn below_sea_spillway_obeys_invariants() {
        use crate::lakes::connectivity::water_class;
        let (w, h) = (24usize, 3usize);
        let mut hm = GridF32::new(w, h, 0.9); // high walls everywhere
        for y in 0..h {
            hm.set(0, y, 0.2); // ocean column (class-1)
        }
        // An arm of the sea reaching inland to x=3 (all ≤ sea, 8-connected to the border → ocean).
        hm.set(1, 1, 0.40);
        hm.set(2, 1, 0.45);
        hm.set(3, 1, 0.49);
        hm.set(4, 1, 0.5001); // the pit's LOCAL SILL — a NEAR-SEA col (Finding 38: a below-sea sea can
        // only overflow a rim within ~2 m of sea; a high col would make it endorheic instead).
        for x in 5..=8 {
            hm.set(x, 1, 0.498); // enclosed below-sea pit floor (4 cells, ≈ −22 m)
        }
        // A ramp draining LEFT into the pit (inflow ⟹ overflow over the near-sea sill toward the sea arm).
        for (i, x) in (9..=22).enumerate() {
            hm.set(x, 1, 0.60 + 0.02 * i as f32);
        }
        let precip = GridF32::new(w, h, 1.0);
        let temp = GridF32::new(w, h, 10.0);
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let ss = SteinSteinParams::default();
        let r = below_sea_basin_lakes(&hm, &clim, &C1DrainageConfig::default(), &ss, 24.0, None);
        let wc = water_class(&hm, C1_SEA_LEVEL_NORM);
        assert!(!r.lakes.is_empty(), "the fed below-sea pit must be inventoried");
        // Invariant 1 — no inventoried lake is dry: each has inflow > 0 OR a positive depth.
        for lk in &r.lakes {
            let b = r.basins.iter().find(|b| b.id == lk.base.id).unwrap();
            assert!(
                b.inflow_m3s > 0.0 || lk.depth_m > 0.0,
                "lake #{} is a DRY depression (no inflow, no depth) — must not be inventoried",
                lk.base.id
            );
        }
        // Invariants 2 & 3 — every spillway.
        for sw in &r.spillways {
            let &(ex, ey) = sw.points.last().unwrap();
            let end = ey as usize * w + ex as usize;
            // (2) arrives elsewhere: ocean or a DIFFERENT lake, never its own footprint.
            assert_ne!(
                r.lake_map[end], sw.lake_id,
                "spillway #{} loops back into its own lake",
                sw.lake_id
            );
            assert!(
                wc[end] == 1 || (r.lake_map[end] != 0 && r.lake_map[end] != sw.lake_id),
                "spillway #{} does not arrive at a different object (ocean/other lake)",
                sw.lake_id
            );
            // (3, false-positive-free form) no interior point re-enters the source lake's own footprint
            // (the loop catch). A LEVEL/SEA threshold was measured to demote legitimate below-sea chains
            // and sea-overflowing lakes (ADR 37b), so the invariant is "never re-enter one's OWN hollow".
            for &(px, py) in &sw.points[1..sw.points.len().saturating_sub(1)] {
                let k = py as usize * w + px as usize;
                assert_ne!(
                    r.lake_map[k], sw.lake_id,
                    "spillway #{} re-enters its own footprint at ({px},{py}) — a loop",
                    sw.lake_id
                );
            }
        }
        // The pit overflows to the sea arm → EXORHEIC with a spillway.
        assert!(
            r.spillways.iter().any(|s| s.lake_id == r.lakes[0].base.id),
            "the fed pit must have a spillway"
        );
    }

    /// Finding 31 INVARIANT (permanent, non-ignored) — SINK VALIDITY is decoupled from the
    /// inventory threshold, and the SEA is `water_class` alone. A SUB-threshold enclosed
    /// below-sea basin (< 5 km²) must still be MARKED in `lake_map` (so a river ending in it
    /// is not mislabelled "sea"), must NOT appear in the exported inventory, and no marked
    /// basin cell may be ocean (`water_class == 1`). This is the class of bug this thread kept
    /// hitting: a property (the sea label, a sink) read from a proxy (altitude, an area
    /// threshold) instead of the authority (`water_class`, the traced network).
    #[test]
    fn below_sea_sink_decoupled_from_inventory_and_sea_is_water_class() {
        use crate::lakes::connectivity::water_class;
        let (w, h) = (24usize, 12usize);
        let mut hm = GridF32::new(w, h, 0.7);
        for y in 0..h {
            hm.set(0, y, 0.2); // ocean column (class-1)
            hm.set(1, y, 0.2);
        }
        // A TINY enclosed below-sea pit (2 cells → ~2 km² at this scale, below the 5 km²
        // inventory threshold), with a low sill lane toward the ocean.
        hm.set(6, 6, 0.4);
        hm.set(7, 6, 0.4);
        for x in 2..6 {
            hm.set(x, 6, 0.55);
        }
        let precip = GridF32::new(w, h, 1.0);
        let temp = GridF32::new(w, h, 10.0);
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let r = below_sea_basin_lakes(
            &hm,
            &clim,
            &C1DrainageConfig::default(),
            &SteinSteinParams::default(),
            24.0,
            None,
        );
        let wc = water_class(&hm, C1_SEA_LEVEL_NORM);
        // Marked as a sink despite being sub-threshold.
        let marked = r.lake_map.iter().filter(|&&x| x != 0).count();
        assert!(marked > 0, "a sub-threshold below-sea basin MUST still be marked (sink validity)");
        // Inventory floor is now a few CELLS (Finding 33 PART A) — reject single-cell noise but
        // keep visible lakes. The 2-cell pit here is below the 4-cell floor → marked, not listed.
        let cell_km2 = (24.0f32 / w as f32).powi(2);
        assert!(
            r.lakes.iter().all(|l| l.area_km2 >= (INVENTORY_MIN_CELLS as f32 - 0.5) * cell_km2),
            "the inventory must reject sub-{INVENTORY_MIN_CELLS}-cell noise"
        );
        // The SEA authority: a marked basin cell is never ocean.
        for k in 0..r.lake_map.len() {
            if r.lake_map[k] != 0 {
                assert_ne!(
                    wc[k], 1,
                    "a marked below-sea basin cell must NOT be water_class==1 (ocean)"
                );
            }
        }
    }

    /// DEFECT A guard — `clip_rivers_to_lakes` makes every watercourse terminate at
    /// its sink: no segment retains a point inside a lake, an EXORHEIC lake keeps its
    /// outlet reach, an ENDORHEIC lake's phantom outlet is dropped.
    fn one_seg_crossing_lake(lake_type: LakeType) -> C1DrainageResult {
        use crate::terrain::flow::{FlowResult, RiverNetwork};
        let (w, h) = (10usize, 1usize);
        let mut lake_map = vec![0u32; w * h];
        lake_map[4] = 7;
        lake_map[5] = 7; // a two-cell lake spanning x=4..=5
        let seg = RiverSegment {
            points: (0..w as u32).map(|x| (x, 0u32)).collect(),
            strahler_order: 2,
            avg_flow: 10.0,
            max_flow: 20.0,
            basin_id: 1,
            upstream: vec![],
            downstream: None,
        };
        let base = Lake {
            id: 7,
            surface_elevation: 0.5,
            max_depth: 0.1,
            area: 2,
            basin_id: 1,
            outlet: (0, 0),
            shallow: false,
        };
        C1DrainageResult {
            flow: FlowResult {
                filled: GridF32::new(w, h, 0.0),
                direction: vec![0; w * h],
                accumulation: GridF32::new(w, h, 0.0),
                basins: vec![0; w * h],
                num_basins: 1,
            },
            rivers: RiverNetwork { segments: vec![seg] },
            segment_drainage_km2: vec![100.0],
            segment_navigability: vec![Navigability::NonNavigable],
            segment_discharge_m3s: vec![25.0],
            segment_width_m: vec![25.0],
            segment_profile_m: vec![(0..w).map(|x| x as f32).collect()],
            lakes: vec![C1Lake { base, level_m: 0.0, depth_m: 1.0, area_km2: 1.0, lake_type }],
            lake_map,
            width: w,
            height: h,
        }
    }

    #[test]
    fn clip_rivers_terminate_at_lake_sinks() {
        let w = 10usize;
        // Exorheic: inflow reach (x0..3) + outlet reach (x6..9) both survive.
        let mut ex = one_seg_crossing_lake(LakeType::Exorheic);
        clip_rivers_to_lakes(&mut ex);
        assert_eq!(ex.rivers.segments.len(), 2, "exorheic: inflow + outlet reaches");
        for (i, s) in ex.rivers.segments.iter().enumerate() {
            assert!(
                s.points.iter().all(|&(x, y)| ex.lake_map[y as usize * w + x as usize] == 0),
                "no clipped point lies inside a lake"
            );
            assert_eq!(s.points.len(), ex.segment_profile_m[i].len(), "profile parallel to points");
        }
        // Endorheic: the closed basin has no outflow → only the inflow reach survives.
        let mut en = one_seg_crossing_lake(LakeType::Endorheic);
        clip_rivers_to_lakes(&mut en);
        assert_eq!(en.rivers.segments.len(), 1, "endorheic: phantom outlet dropped");
        assert!(en.rivers.segments[0].points.iter().all(|&(x, _)| x < 4));
    }

    /// Finding 37 INVARIANT (permanent, non-ignored) — the exorheic-outlet guard over the WHOLE
    /// lake population, ACROSS PROVENANCES. The prior guard iterated only below-sea `r.basins` on a
    /// synthetic grid, so the 21 exorheic below-sea lakes shipped without an outlet in the 8192²
    /// export slipped through (subset blind spot, third occurrence — see ADR Finding 37). This
    /// checks the checker itself: an exorheic lake with a bordering outlet segment passes; one
    /// without is flagged; an endorheic lake is never required to have one.
    #[test]
    fn every_exorheic_lake_needs_a_traced_outlet() {
        use crate::terrain::flow::{FlowResult, RiverNetwork};
        let (w, h) = (12usize, 1usize);
        let mut lake_map = vec![0u32; w * h];
        lake_map[4] = 7; // lake A (exorheic, WILL have an outlet segment)
        lake_map[5] = 7;
        lake_map[9] = 9; // lake B (exorheic, NO outlet segment)
        lake_map[1] = 11; // lake C (endorheic — never required to have an outlet)
        let mk_lake = |id: u32, lake_type: LakeType| C1Lake {
            base: Lake {
                id,
                surface_elevation: 0.5,
                max_depth: 0.1,
                area: 1,
                basin_id: 1,
                outlet: (0, 0),
                shallow: false,
            },
            level_m: 0.0,
            depth_m: 1.0,
            area_km2: 1.0,
            lake_type,
        };
        let outlet_seg = |sx: u32| RiverSegment {
            points: vec![(sx, 0), (sx + 1, 0)],
            strahler_order: 1,
            avg_flow: 1.0,
            max_flow: 1.0,
            basin_id: 1,
            upstream: vec![],
            downstream: None,
        };
        let mk = |segs: Vec<RiverSegment>| C1DrainageResult {
            flow: FlowResult {
                filled: GridF32::new(w, h, 0.0),
                direction: vec![0; w * h],
                accumulation: GridF32::new(w, h, 0.0),
                basins: vec![0; w * h],
                num_basins: 1,
            },
            segment_drainage_km2: vec![1.0; segs.len()],
            segment_navigability: vec![Navigability::NonNavigable; segs.len()],
            segment_discharge_m3s: vec![1.0; segs.len()],
            segment_width_m: vec![1.0; segs.len()],
            segment_profile_m: segs.iter().map(|s| vec![0.0; s.points.len()]).collect(),
            rivers: RiverNetwork { segments: segs },
            lakes: vec![
                mk_lake(7, LakeType::Exorheic),
                mk_lake(9, LakeType::Exorheic),
                mk_lake(11, LakeType::Endorheic),
            ],
            lake_map: lake_map.clone(),
            width: w,
            height: h,
        };
        // A has an outlet segment (source x=6 borders lake cell 5); B has none; C is endorheic.
        let dr = mk(vec![outlet_seg(6)]);
        assert_eq!(
            exorheic_lakes_missing_outlet(&dr),
            vec![9],
            "only lake B (exorheic, no bordering outlet) must be flagged; C is endorheic"
        );
        // Give B an outlet too (source x=8 borders lake cell 9) → no lake is flagged.
        let dr = mk(vec![outlet_seg(6), outlet_seg(8)]);
        assert!(
            exorheic_lakes_missing_outlet(&dr).is_empty(),
            "every exorheic lake now has a bordering outlet"
        );
    }

    /// A tilted plane draining to an ocean edge: the stack runs at sea=0.5,
    /// produces rivers, and km² drainage areas are finite + monotone with
    /// max_flow.
    #[test]
    fn c1_drainage_runs_on_tilted_plane() {
        let n = 64;
        let mut hm = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                // land on the left (>0.5), ocean on the right (<0.5).
                let t = i as f32 / n as f32;
                hm.set(i, j, 0.85 - t * 0.6);
            }
        }
        let out =
            c1_drainage(&hm, None, &C1DrainageConfig::default(), &SteinSteinParams::default());
        assert_eq!(out.width, n);
        assert_eq!(out.segment_drainage_km2.len(), out.rivers.segments.len());
        assert!(out.segment_drainage_km2.iter().all(|v| v.is_finite() && *v >= 0.0));
        // All lakes (if any) carry finite metre stats.
        for lk in &out.lakes {
            assert!(lk.level_m.is_finite() && lk.depth_m.is_finite() && lk.depth_m >= 0.0);
            assert!(lk.area_km2 > 0.0);
        }
    }

    /// Navigability classification is monotone in drainage area.
    #[test]
    fn navigability_thresholds_monotone() {
        let t = DrainageThresholds::default();
        assert_eq!(t.classify(10.0), Navigability::NonNavigable);
        assert_eq!(t.classify(1_000.0), Navigability::SmallBoat);
        assert_eq!(t.classify(8_000.0), Navigability::Barge);
        assert_eq!(t.classify(80_000.0), Navigability::Ship);
    }
}
