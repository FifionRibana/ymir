//! GeoJSON vector layers for the v1 `.ymir` container — coastline and cliffs.
//!
//! Emits GeoJSON BY HAND with `serde_json` (no `geojson` dependency). Geometry
//! comes from [`crate::terrain::contour`] marching squares; coordinates are in
//! CELL space, identical to the container's rasters (see that module's header
//! for the orientation invariant). Output bytes are deterministic: same field →
//! same geometry → same JSON.

use serde_json::{Value, json};

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use crate::terrain::contour::{
    COASTLINE_SMOOTH_PASSES, Polyline, SMOOTH_RELAX_BELOW_DEG, marching_squares,
    slope_deg_to_norm_gradient, slope_degree_field, smooth_polylines_on_isoline,
};

/// Sea level in the C1 normalized altitude field. The vertical contract
/// ([`c1_altitude_norm_to_metres`]) pins sea to `0.5` (0 m), so tracing the
/// normalized field at `0.5` is exactly the metric altitude-0 isoline — no
/// metres grid needed for the coastline.
pub const SEA_LEVEL_NORM: f32 = 0.5;

/// Default cliff slope threshold in degrees (WP-3 oracle default).
pub const DEFAULT_CLIFF_THRESHOLD_DEG: f32 = 45.0;

/// Wrap polylines as a GeoJSON `FeatureCollection` with a single
/// `MultiLineString` feature carrying `properties`. Returns pretty JSON bytes.
pub fn polylines_to_geojson(polylines: &[Polyline], properties: Value) -> Vec<u8> {
    let coordinates: Vec<Vec<[f64; 2]>> = polylines
        .iter()
        .map(|pl| pl.iter().map(|&(x, y)| [x as f64, y as f64]).collect())
        .collect();
    let fc = json!({
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": properties,
            "geometry": { "type": "MultiLineString", "coordinates": coordinates }
        }]
    });
    // `to_vec_pretty` is deterministic; `unwrap` is safe (a FeatureCollection of
    // finite numbers always serializes).
    serde_json::to_vec_pretty(&fc).expect("geojson serialization is infallible here")
}

/// Trace the coastline (sea-level isoline) of the normalized eroded heightmap
/// and serialize it as `coastline.geojson` bytes. Documented choice: run on the
/// NORMALIZED field at [`SEA_LEVEL_NORM`] (equivalent to metric altitude 0).
///
/// `smooth` opts into the gradient-weighted contour relaxation (ADR Finding 48).
/// **It ships OFF, and the measurement is why**: on the flattest shores it removes 19 % of the
/// >80° turns, but it ADDS 15 % at 2–5°, and the TOTAL over all classes moves 20.13 % → 20.04 %
/// — a −0.4 % change, i.e. nothing. It redistributes barbs rather than removing them. What it
/// does change is the axial concentration (0.377 → 0.246), so the contour is measurably less
/// axis-aligned; whether that is visible is a judgement for the author, which is exactly why
/// this is a flag and not a default.
pub fn coastline_geojson(eroded_norm: &GridF32, smooth: Option<CoastlineSmoothing>) -> Vec<u8> {
    let mut polylines = marching_squares(eroded_norm, SEA_LEVEL_NORM);
    if let Some(o) = smooth {
        let gate = slope_deg_to_norm_gradient(
            SMOOTH_RELAX_BELOW_DEG,
            o.m_per_cell,
            2.0 * 1.13 * o.depth_scale_m,
        );
        polylines = smooth_polylines_on_isoline(
            eroded_norm,
            SEA_LEVEL_NORM,
            &polylines,
            COASTLINE_SMOOTH_PASSES,
            gate,
        );
    }
    polylines_to_geojson(&polylines, json!({ "id": "coastline", "level_m": 0.0 }))
}

/// The two physical quantities the coastline relaxation needs to turn a slope in DEGREES into a
/// gradient in normalised units per cell. Passing them explicitly keeps the unit conversion at
/// the caller, which is the only place that knows them.
#[derive(Debug, Clone, Copy)]
pub struct CoastlineSmoothing {
    pub m_per_cell: f32,
    pub depth_scale_m: f32,
}

/// Trace cliff edges: build the metric height (via the vertical contract), the
/// per-cell slope-in-degrees field (using `cell_size_m` for a REAL angle), then
/// the `slope_deg >= threshold_deg` isoline. Serialized as `cliffs.geojson`.
///
/// Follow-up: per-segment mean cliff height (Δaltitude across the edge) as a
/// feature property — skipped here (a single MultiLineString feature).
pub fn cliffs_geojson(
    eroded_norm: &GridF32,
    ss: &SteinSteinParams,
    cell_size_m: f32,
    threshold_deg: f32,
) -> Vec<u8> {
    let (w, h) = (eroded_norm.width, eroded_norm.height);
    let metres: Vec<f32> =
        eroded_norm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, ss)).collect();
    let height_m = GridF32::from_vec(w, h, metres);
    let slope_deg = slope_degree_field(&height_m, cell_size_m);
    let polylines = marching_squares(&slope_deg, threshold_deg);
    polylines_to_geojson(
        &polylines,
        json!({ "id": "cliffs", "slope_threshold_deg": threshold_deg as f64 }),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feature_line_count(bytes: &[u8]) -> usize {
        let v: Value = serde_json::from_slice(bytes).unwrap();
        v["features"][0]["geometry"]["coordinates"].as_array().unwrap().len()
    }

    /// A synthetic island: coastline is a non-empty, closed loop, deterministic.
    #[test]
    fn coastline_island_closed_and_deterministic() {
        let (w, h) = (48usize, 48usize);
        let (cx, cy) = (23.5f32, 23.5f32);
        let mut data = vec![0.2f32; w * h]; // sea everywhere (< 0.5)
        for y in 0..h {
            for x in 0..w {
                let d = ((x as f32 - cx).powi(2) + (y as f32 - cy).powi(2)).sqrt();
                // A central landmass rising above sea level (0.5).
                data[y * w + x] = (0.2 + (1.0 - d / 12.0).max(0.0) * 0.6).min(0.9);
            }
        }
        let grid = GridF32::from_vec(w, h, data);
        let bytes = coastline_geojson(&grid, None);

        assert!(feature_line_count(&bytes) > 0, "coastline must have geometry");

        let v: Value = serde_json::from_slice(&bytes).unwrap();
        let ring = &v["features"][0]["geometry"]["coordinates"][0];
        let pts = ring.as_array().unwrap();
        assert!(pts.len() > 4, "island ring must have vertices");
        assert_eq!(pts.first(), pts.last(), "coastline ring must be a closed loop");
        assert_eq!(v["features"][0]["properties"]["level_m"], json!(0.0));

        // Byte-identical across two runs.
        assert_eq!(bytes, coastline_geojson(&grid, None), "coastline bytes must be deterministic");
    }

    /// A vertical step: cliffs appear below the step angle, none above it.
    #[test]
    fn cliffs_present_below_step_angle_absent_above() {
        let ss = SteinSteinParams::default();
        let (w, h) = (16usize, 16usize);
        // Left half sea (norm 0.5 → 0 m), right half raised (norm 0.6 → ~1130 m).
        let mut data = vec![0.5f32; w * h];
        for y in 0..h {
            for x in (w / 2)..w {
                data[y * w + x] = 0.6;
            }
        }
        let grid = GridF32::from_vec(w, h, data);
        // cell_size_m = 100 m → step gradient ~565 m/cell → slope ~80°.
        let below = cliffs_geojson(&grid, &ss, 100.0, 30.0);
        let above = cliffs_geojson(&grid, &ss, 100.0, 85.0);

        assert!(feature_line_count(&below) > 0, "threshold below the step angle → cliff polylines");
        assert_eq!(feature_line_count(&above), 0, "threshold above the step angle → no cliffs");

        let v: Value = serde_json::from_slice(&below).unwrap();
        assert_eq!(v["features"][0]["properties"]["slope_threshold_deg"], json!(30.0));
    }
}
