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
use crate::lakes::detection::{detect_lakes, Lake, LakeConfig};
use crate::terrain::flow::{
    compute_flow, extract_rivers, FlowConfig, FlowResult, RiverConfig, RiverNetwork, D8_DX, D8_DY,
    DIR_NONE,
};

use super::closures::oceanic_bathymetry::params::SteinSteinParams;
use super::production_upscale::{c1_altitude_norm_to_metres, c1_cell_area_km2};

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
}

impl Default for C1DrainageConfig {
    fn default() -> Self {
        // 10 m min lake depth, 5 km² min lake area — plausible mappable lakes.
        Self { thresholds: DrainageThresholds::default(), lake_min_depth_m: 10.0, lake_min_area_km2: 5.0 }
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
pub struct C1DrainageResult {
    pub flow: FlowResult,
    pub rivers: RiverNetwork,
    /// Per-segment upstream drainage area in km² (parallel to `rivers.segments`).
    pub segment_drainage_km2: Vec<f32>,
    /// Per-segment navigability class (parallel to `rivers.segments`).
    pub segment_navigability: Vec<Navigability>,
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
    let (w, h) = (heightmap.width, heightmap.height);
    let cell_km2 = c1_cell_area_km2(w);
    // metres per unit of normalised altitude (the linear vertical scale slope).
    let metres_per_norm = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);

    // 1. Flow at the unified sea level (NOT the legacy 0.1).
    let flow = compute_flow(heightmap, &FlowConfig { sea_level: C1_SEA_LEVEL_NORM });

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

    // Per-segment drainage area (km²) from the segment's max accumulation, and
    // the navigability class.
    let segment_drainage_km2: Vec<f32> =
        rivers.segments.iter().map(|s| s.max_flow * cell_km2).collect();
    let segment_navigability: Vec<Navigability> =
        segment_drainage_km2.iter().map(|&km2| cfg.thresholds.classify(km2)).collect();

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
        Some(clim) => {
            water_balance_lakes(heightmap, &flow, clim, cell_km2, ss, &lake_result.lakes, &mut lake_map, w, h)
        }
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

    C1DrainageResult {
        flow,
        rivers,
        segment_drainage_km2,
        segment_navigability,
        lakes,
        lake_map,
        width: w,
        height: h,
    }
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
    // 1. Runoff per land cell (mm·km²/yr surplus) = max(0, precip − PE)·cell_km2.
    let mut runoff = vec![0.0f32; n];
    for k in 0..n {
        if heightmap.data[k] > C1_SEA_LEVEL_NORM {
            let p = precip_mm_per_year(climate.precip_internal.data[k]);
            let pe = potential_evaporation_mm(climate.temperature.data[k]);
            runoff[k] = (p - pe).max(0.0) * cell_km2;
        }
    }
    // 2. Accumulate runoff downstream (decreasing filled height — the flow
    //    topological order; the flat tiebreak is omitted, negligible for inflow).
    let mut order: Vec<usize> =
        (0..n).filter(|&k| heightmap.data[k] > C1_SEA_LEVEL_NORM).collect();
    order.sort_unstable_by(|&a, &b| {
        flow.filled.data[b].partial_cmp(&flow.filled.data[a]).unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut runoff_accum = runoff;
    for &k in &order {
        let d = flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let (i, j) = (k % w, k / w);
        let ni = ((i as i32 + D8_DX[d as usize]).rem_euclid(w as i32)) as usize;
        let nj = ((j as i32 + D8_DY[d as usize]).rem_euclid(h as i32)) as usize;
        runoff_accum[nj * w + ni] += runoff_accum[k];
    }
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
        let out = c1_drainage(&hm, None, &C1DrainageConfig::default(), &SteinSteinParams::default());
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
