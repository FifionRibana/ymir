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
}

impl Default for DrainageThresholds {
    fn default() -> Self {
        Self { stream_km2: 20.0, small_boat_km2: 500.0, barge_km2: 5_000.0, ship_km2: 50_000.0 }
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
            let q_m3s = dr_km2
                .iter()
                .map(|&km2| runoff_km2_to_m3s(REFERENCE_RUNOFF_MM * km2))
                .collect();
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
pub fn apply_geo_scale_ratio(dr: &mut C1DrainageResult, ratio: f32, thresholds: &DrainageThresholds) {
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
    let metres_per_norm =
        c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);

    // Finding 32 — spillway routing by a SINGLE priority-flood from the ocean (Barnes 2014),
    // computed ONCE. `spill_receiver[k]` is the neighbour toward the ocean along the LEAST-
    // BARRIER (minimax-elevation) path; following it from any cell reaches `water_class == 1`
    // over the lowest sill. This replaces the per-basin bounded Dijkstra, which failed for
    // basins whose ocean was beyond its step budget — leaving exorheic basins with no outlet.
    // O(n log n) once, so a distant/large basin is handled at the same cost as a coastal one.
    // `barrier_q[k]` (Finding 33) is the least MAX-elevation to reach the ocean from `k` — i.e.
    // the true overflow SILL of the basin `k` sits in. A basin FILLS to this sill (not to the
    // shoreline), so its level/area/depth derive from the authority (the flood), not the proxy
    // (the min external neighbour of the below-sea cells, which is only the ~0 m shore).
    let quant = |e: f32| (e * 1_000_000.0) as i32;
    let mut spill_receiver = vec![u32::MAX; n];
    let mut barrier_q = vec![i32::MAX; n];
    if wc.iter().any(|&c| c == 2) {
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;
        let mut done = vec![false; n];
        let mut pq: BinaryHeap<Reverse<(i32, u32)>> = BinaryHeap::new();
        for k in 0..n {
            if wc[k] == 1 {
                done[k] = true;
                barrier_q[k] = quant(heightmap.data[k]);
                pq.push(Reverse((barrier_q[k], k as u32)));
            }
        }
        while let Some(Reverse((c, k))) = pq.pop() {
            let (x, y) = ((k as usize % w) as i32, (k as usize / w) as i32);
            for (dx, dy) in D8_DX.iter().zip(D8_DY.iter()) {
                if let Some(nk) = nb4(x + dx, y + dy) {
                    if !done[nk] {
                        done[nk] = true;
                        spill_receiver[nk] = k; // toward the ocean along the least-barrier path
                        let nb = c.max(quant(heightmap.data[nk]));
                        barrier_q[nk] = nb;
                        pq.push(Reverse((nb, nk as u32)));
                    }
                }
            }
        }
    }

    // Finding 33 — a below-sea LAKE is a connected pool of UNDERWATER cells (a cell is underwater
    // when its escape barrier is above it: `barrier > height`), FILLED to the shared sill. This
    // is the depression itself, so adjacent sub-pockets behind the same sill MERGE automatically
    // (no per-class-2-component double counting). Kept only when the pool's floor is below sea.
    let underwater = |k: usize| barrier_q[k] != i32::MAX && barrier_q[k] > quant(heightmap.data[k]);
    for s in 0..n {
        if seen[s] || !underwater(s) {
            continue;
        }
        // 1. the connected underwater pool (the filled depression, already at its sill).
        let mut comp = Vec::new();
        let mut q = VecDeque::new();
        q.push_back(s);
        seen[s] = true;
        while let Some(k) = q.pop_front() {
            comp.push(k);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1)] {
                if let Some(nk) = nb4(x + dx, y + dy) {
                    if underwater(nk) && !seen[nk] {
                        seen[nk] = true;
                        q.push_back(nk);
                    }
                }
            }
        }
        // Below-sea only: the pool's FLOOR must be below sea (else it is an ordinary above-sea
        // lake handled by `detect_lakes`).
        let floor_norm = comp.iter().map(|&k| heightmap.data[k]).fold(f32::MAX, f32::min);
        if floor_norm >= C1_SEA_LEVEL_NORM {
            continue;
        }
        // Finding 34 — TOTAL INFLOW from the runoff FIELD (the authority), SUMMED at the shoreline:
        // `runoff_accumulation` ZEROES below-sea cells, so a tributary's discharge is lost the
        // instant it enters the water. Reading `max` over the pool therefore saw only the largest
        // single stream, not the total. Instead SUM the accumulated runoff of every above-sea cell
        // that drains INTO a below-sea pool cell — each tributary counted once, at the shore, before
        // the zeroing. Read from the FINAL footprint, so a track ending a cell short is irrelevant.
        let mut inflow = 0.0f32;
        for &k in &comp {
            if heightmap.data[k] > C1_SEA_LEVEL_NORM {
                continue; // want the below-sea WATER cells (where zeroing happens)
            }
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1)] {
                if let Some(nk) = nb4(x + dx, y + dy) {
                    if heightmap.data[nk] > C1_SEA_LEVEL_NORM {
                        let dir = flow.direction[nk];
                        if dir != DIR_NONE {
                            let (tx, ty) =
                                (nk as i32 % w as i32 + D8_DX[dir as usize], nk as i32 / w as i32 + D8_DY[dir as usize]);
                            if nb4(tx, ty) == Some(k) {
                                inflow += runoff[nk]; // this tributary enters the water here
                            }
                        }
                    }
                }
            }
        }
        // 2. the sill (overflow level) = the pool's barrier to the ocean (uniform within the pool).
        let spill = comp.iter().map(|&k| barrier_q[k]).min().unwrap() as f32 / 1_000_000.0;
        let mut fcells: Vec<usize> = comp.clone();
        fcells.sort_by(|&a, &b| heightmap.data[a].partial_cmp(&heightmap.data[b]).unwrap());
        let a_spill = fcells.len() as f32 * cell_km2;
        let pe_lake = potential_evaporation_mm(climate.temperature.data[fcells[0]]).max(1.0);
        let a_eq = inflow / pe_lake;
        let (level, lake_type) = if a_eq >= a_spill {
            (spill, LakeType::Exorheic)
        } else {
            let n_eq = (a_eq / cell_km2).floor().max(1.0) as usize;
            (heightmap.data[fcells[n_eq.min(fcells.len()) - 1]], LakeType::Endorheic)
        };
        // 4. mark the WATER cells (≤ level) and build the typed lake.
        let id = next_id;
        next_id += 1;
        let mut water = 0usize;
        let mut floor = f32::MAX;
        let (mut ox, mut oy) = (0u32, 0u32);
        for &k in &fcells {
            if heightmap.data[k] <= level {
                lake_map[k] = id;
                water += 1;
                // Finding 30 — wetland vs lagoon by DEPTH: a shallow (< 3 m) water cell is a
                // through-flow wetland margin; a deeper one is open water (lagoon/inland sea).
                if (level - heightmap.data[k]) * metres_per_norm < WETLAND_MAX_DEPTH_M {
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
        let level_m = c1_altitude_norm_to_metres(level, ss);
        let floor_m = c1_altitude_norm_to_metres(floor, ss);
        // Inventory (lakes.json / microscope) keeps the 5 km² threshold — no micro-lakes.
        // The basin is still marked + gets a spillway below regardless (sink validity).
        if water >= min_cells {
            out.push(C1Lake {
                base: Lake {
                    id,
                    surface_elevation: level,
                    max_depth: level - floor,
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
                    if !underwater(nk) {
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
            evaporation_m3s: runoff_km2_to_m3s(pe_lake * area_km2_bs),
            a_eq_km2: a_eq,
            a_spill_km2: a_spill,
            spill_level_m: c1_altitude_norm_to_metres(spill, ss),
            floor_m,
            level_m,
            max_depth_m: level_m - floor_m,
            area_at_sill_km2: a_spill,
            exorheic_before_inlet_fix: (inflow_max / pe_lake) >= a_spill,
        });
        // Finding 30 + 32 — an EXORHEIC label REQUIRES a traced outlet (mass balance: a basin
        // that receives more than it evaporates MUST overflow). Follow the ocean priority-flood's
        // `spill_receiver` from the basin's pour point over the lowest sill to a sink — the OCEAN
        // (`water_class == 1`) or, if the path first re-enters another below-sea basin, that basin
        // (a chain). Guaranteed to reach a sink for every basin the flood reached (i.e. every one).
        if lake_type == LakeType::Exorheic {
            // Pour = the basin cell whose receiver leaves the basin (its outflow), lowest sill.
            let mut pour: Option<usize> = None;
            let mut pour_e = f32::MAX;
            for &k in &fcells {
                if lake_map[k] != id {
                    continue;
                }
                let r = spill_receiver[k];
                if r != u32::MAX && lake_map[r as usize] != id && heightmap.data[r as usize] < pour_e
                {
                    pour_e = heightmap.data[r as usize];
                    pour = Some(k);
                }
            }
            if let Some(p0) = pour {
                let mut path = vec![p0];
                let (mut cur, mut chained, mut steps) = (p0, None, 0usize);
                loop {
                    if wc[cur] == 1 {
                        break; // reached the ocean
                    }
                    let lid = lake_map[cur];
                    if lid != 0 && lid != id && cur != p0 {
                        chained = Some(lid); // spilled into another below-sea basin (a chain)
                        break;
                    }
                    let r = spill_receiver[cur];
                    if r == u32::MAX {
                        break;
                    }
                    cur = r as usize;
                    path.push(cur);
                    steps += 1;
                    if steps > n {
                        break;
                    }
                }
                if wc[cur] == 1 || chained.is_some() {
                    let points: Vec<(u32, u32)> = path.iter().map(|&k| ((k % w) as u32, (k / w) as u32)).collect();
                    let profile: Vec<f32> = path.iter().map(|&k| c1_altitude_norm_to_metres(heightmap.data[k], ss)).collect();
                    // discharge = surplus (inflow − evaporation); evaporation = pe_lake · area.
                    let area_km2 = water as f32 * cell_km2;
                    let surplus = (inflow - pe_lake * area_km2).max(0.0); // mm·km²/yr
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
                        chained_into: chained,
                    });
                }
            }
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
        // low lane at y=6 across the whole width (0.55) so a spill route exists to the coast.
        for x in 2..w {
            hm.set(x, 6, 0.55);
        }
        for y in 5..8 {
            for x in 70..73 {
                hm.set(x, y, 0.4); // FAR enclosed below-sea pit (large-ish)
            }
        }
        hm.set(40, 6, 0.4); // a TINY 1-cell below-sea pit mid-grid (sub-inventory)
        let precip = GridF32::new(w, h, 1.0);
        let temp = GridF32::new(w, h, 10.0);
        let clim = DrainageClimate { precip_internal: &precip, temperature: &temp };
        let r = below_sea_basin_lakes(
            &hm,
            &clim,
            &C1DrainageConfig::default(),
            &SteinSteinParams::default(),
            80.0,
        );
        let wc = water_class(&hm, C1_SEA_LEVEL_NORM);
        // The INVARIANT over EVERY basin (any area): exorheic ⟹ a traced spillway reaching a sink.
        let exo = r.basins.iter().filter(|b| b.exorheic).count();
        assert!(exo > 0, "at least one below-sea basin must classify EXORHEIC");
        for b in r.basins.iter().filter(|b| b.exorheic) {
            let sw = r
                .spillways
                .iter()
                .find(|s| s.lake_id == b.id)
                .unwrap_or_else(|| panic!("EXORHEIC basin #{} ({:.2} km²) MUST have a traced spillway", b.id, b.area_km2));
            let &(lx, ly) = sw.points.last().unwrap();
            let end = ly as usize * w + lx as usize;
            assert!(
                wc[end] == 1 || sw.chained_into.is_some(),
                "the spillway must reach a SINK (ocean or another basin), ended at wc={}",
                wc[end]
            );
        }
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
                assert_ne!(wc[k], 1, "a marked below-sea basin cell must NOT be water_class==1 (ocean)");
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
