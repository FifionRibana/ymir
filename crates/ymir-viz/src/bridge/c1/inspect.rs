//! Per-cell inspection over the HD product (UI rewrite step c/5).
//!
//! All per-cell quantities are ALREADY readable from `ymir-core` after
//! generation (verified — the inventory was wrong about missing getters):
//! `GridF32::get` (altitude / bathymetry / temperature / precipitation),
//! `Field2D::get` (S̃ / age), `BoolField::get` (craton), `PlateIdField::get`
//! / `PlateTypeField::get` (plate id / type), the `Vec<Biome>` and the
//! drainage arrays (`flow.direction`, `accumulation`, `basins`, `lake_map`)
//! by row-major index. So NO core accessor needed to be added.
//!
//! The ONE non-trivial gap is the **cell → river-segment map**: drainage
//! stores rivers per SEGMENT (points + navigability + drainage area), not
//! per cell, so "is this cell on a river, and which?" needs a reverse map.
//! [`RiverCellMap`] builds it ON DEMAND (the UI memoises it; it is NOT a
//! field on `C1DrainageResult` — only inspection uses it).
//!
//! [`inspect_cell`] is the optional convenience: ONE call → every HD layer
//! for a cell, so step d's inspector has a single access point.

use ymir_core::climate::biomes::Biome;
use ymir_core::climate::precipitation::precip_mm_per_year;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1_SEA_LEVEL_NORM, C1DrainageResult, LakeType, Navigability, potential_evaporation_mm,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;

use super::hd::HdResult;

/// River properties at a cell that sits on a river segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RiverCellInfo {
    /// Index into `C1DrainageResult.rivers.segments`.
    pub segment: usize,
    pub navigability: Navigability,
    pub drainage_km2: f32,
}

/// Reverse map: cell → the river segment passing through it (if any).
/// Built once from a `C1DrainageResult` and memoised by the caller.
pub struct RiverCellMap {
    width: usize,
    height: usize,
    cells: Vec<Option<RiverCellInfo>>,
}

impl RiverCellMap {
    /// Build the per-cell map by walking every segment's points. When two
    /// segments share a cell (junction), the one with the larger drainage
    /// area wins (the main channel).
    pub fn from_drainage(d: &C1DrainageResult) -> Self {
        let (width, height) = (d.width, d.height);
        let mut cells: Vec<Option<RiverCellInfo>> = vec![None; width * height];
        for (si, seg) in d.rivers.segments.iter().enumerate() {
            let info = RiverCellInfo {
                segment: si,
                navigability: d.segment_navigability[si],
                drainage_km2: d.segment_drainage_km2[si],
            };
            for &(px, py) in &seg.points {
                let k = py as usize * width + px as usize;
                let keep = match &cells[k] {
                    Some(prev) => info.drainage_km2 > prev.drainage_km2,
                    None => true,
                };
                if keep {
                    cells[k] = Some(info);
                }
            }
        }
        Self { width, height, cells }
    }

    /// River info at `(x, y)`, or `None` if the cell carries no river.
    pub fn at(&self, x: usize, y: usize) -> Option<RiverCellInfo> {
        if x >= self.width || y >= self.height {
            return None;
        }
        self.cells[y * self.width + x]
    }
}

/// Lake properties at a cell inside a lake.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LakeCellInfo {
    pub lake_type: LakeType,
    pub level_m: f32,
    pub depth_m: f32,
    pub area_km2: f32,
}

/// Every HD quantity at one cell — the single access point for step d's
/// inspector. Coarse tectonic fields (S̃ / plate / craton) live at 64² in
/// the live snapshot, not the HD product; step d maps the HD cell to its
/// coarse cell and reads those directly (their getters already exist).
#[derive(Clone, Copy, Debug)]
pub struct CellInspection {
    pub x: usize,
    pub y: usize,
    pub altitude_norm: f32,
    pub altitude_m: f32,
    pub is_ocean: bool,
    /// Depth below sea level (m), `Some` only for ocean cells.
    pub depth_m: Option<f32>,
    pub temperature_c: f32,
    pub precip_mm: f32,
    /// Water-balance surplus `max(0, precip − PE)` (mm/yr).
    pub runoff_mm: f32,
    pub biome: Biome,
    pub river: Option<RiverCellInfo>,
    pub lake: Option<LakeCellInfo>,
}

/// Gather every HD layer for cell `(x, y)`. `river_map` is the memoised
/// [`RiverCellMap`]; metre conversions use `SteinSteinParams::default()`
/// (the params the HD chain was built with).
pub fn inspect_cell(hd: &HdResult, river_map: &RiverCellMap, x: usize, y: usize) -> CellInspection {
    let ss = SteinSteinParams::default();
    let k = y * hd.width + x;

    let altitude_norm = hd.eroded.get(x as i32, y as i32);
    let altitude_m = c1_altitude_norm_to_metres(altitude_norm, &ss);
    let is_ocean = altitude_norm <= C1_SEA_LEVEL_NORM;
    let depth_m = if is_ocean { Some(-altitude_m) } else { None };

    let temperature_c = hd.temperature.get(x as i32, y as i32);
    let precip_mm = precip_mm_per_year(hd.precipitation.get(x as i32, y as i32));
    let runoff_mm = (precip_mm - potential_evaporation_mm(temperature_c)).max(0.0);
    let biome = hd.biomes[k];

    // Lake membership: lake_map[k] is the lake id (0 = none); resolve to
    // the C1Lake by id (lakes are few — a linear scan is fine).
    let lake = {
        let id = hd.drainage.lake_map[k];
        if id == 0 {
            None
        } else {
            hd.drainage.lakes.iter().find(|l| l.base.id == id).map(|l| LakeCellInfo {
                lake_type: l.lake_type,
                level_m: l.level_m,
                depth_m: l.depth_m,
                area_km2: l.area_km2,
            })
        }
    };

    CellInspection {
        x,
        y,
        altitude_norm,
        altitude_m,
        is_ocean,
        depth_m,
        temperature_c,
        precip_mm,
        runoff_mm,
        biome,
        river: river_map.at(x, y),
        lake,
    }
}
