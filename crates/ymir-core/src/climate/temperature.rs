//! #165 C1 climate — temperature: latitudinal sea-level gradient + adiabatic
//! lapse with altitude. The simple field, computed first (precipitation's
//! Clausius-Clapeyron capacity depends on it).
//!
//! `T(cell) = T_sea(lat) − 6.5 °C/km · altitude_km`. Altitude in metres comes
//! from the vertical coordinate contract (`c1_altitude_norm_to_metres`), so a
//! given thickness yields the same temperature in every seed. Ocean / sub-sea
//! cells take the sea-surface temperature (no lapse under water). Every term is
//! anchored on a real quantity (no free knob): the equator-pole gradient, the
//! environmental lapse rate.

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::production_upscale::{C1_DOMAIN_KM, c1_altitude_norm_to_metres};

/// Environmental (adiabatic) lapse rate — the real value, 6.5 °C per km.
pub const LAPSE_RATE_C_PER_KM: f32 = 6.5;
/// The C1 sea-level normalised value (Maillon 2 unified scale).
pub const SEA_LEVEL_NORM: f32 = 0.5;

/// Sea-level temperature (°C) at a latitude (degrees), anchored on the real
/// annual-mean T-vs-latitude CURVE — equator ~+27 °C, pole ~−25 °C, mid-latitudes
/// convex (flat tropics, steep toward the poles): `T = T_eq − (T_eq − T_pole)·
/// (|lat|/90)^1.8`. The endpoints were always right; the earlier sin²(lat) shape
/// was the bug — it put 45° at the cold MIDPOINT (+1 °C, subarctic) when the real
/// 45° is ~13 °C (temperate). The exponent 1.8 reproduces the observed profile at
/// 0/15/30/45/60/90° (≈27/25/20/12/2/−25 °C), so the slider is right across the
/// whole range, not just the 45° default. (Validated by the biome instrument:
/// this makes the 45° lowlands temperate, not taiga.)
pub fn sea_level_temperature(lat_deg: f32) -> f32 {
    const T_EQ: f32 = 27.0;
    const T_POLE: f32 = -25.0;
    let f = (lat_deg.abs() / 90.0).powf(1.8);
    T_EQ - (T_EQ - T_POLE) * f
}

/// Latitude (degrees) of grid row `j`, full-torus span (`C1_DOMAIN_KM`).
/// Thin wrapper over [`row_latitude_windowed`]; see it for the contract.
pub fn row_latitude(j: usize, ny: usize, lat_deg: f32) -> f32 {
    row_latitude_windowed(j, ny, lat_deg, C1_DOMAIN_KM)
}

/// Latitude (degrees) of grid row `j` when the grid spans `domain_km` north-
/// south: the span is `domain_km / 111` degrees centred on `lat_deg` (the
/// placement parameter). j=0 is the south edge. A cropped window passes its
/// (smaller) `window_km` here so the gradient across the window is correct.
pub fn row_latitude_windowed(j: usize, ny: usize, lat_deg: f32, domain_km: f32) -> f32 {
    row_latitude_span(j, ny, lat_deg, domain_km / 111.0)
}

/// Latitude (degrees) of grid row `j` for an EXPLICIT latitudinal `span_deg`,
/// centred on `centre_deg` (Finding 25 — TASK 2). Decouples the CLIMATIC extent
/// from the physical extent: a 400 km island can be told to span 27° (steep
/// thermal gradient, several wind belts) instead of its geographic ~3.6°. j=0 is
/// the south edge. `span_deg == domain_km / 111` reproduces [`row_latitude_windowed`].
pub fn row_latitude_span(j: usize, ny: usize, centre_deg: f32, span_deg: f32) -> f32 {
    let frac = if ny > 1 { j as f32 / (ny - 1) as f32 } else { 0.5 };
    centre_deg + (frac - 0.5) * span_deg
}

/// Temperature field (°C), full-torus span. Thin wrapper over
/// [`compute_temperature_windowed`] with `domain_km = C1_DOMAIN_KM`.
pub fn compute_temperature(heightmap: &GridF32, ss: &SteinSteinParams, lat_deg: f32) -> GridF32 {
    compute_temperature_windowed(heightmap, ss, lat_deg, C1_DOMAIN_KM)
}

/// Temperature field (°C) for a grid spanning `domain_km` north-south at centre
/// latitude `lat_deg`. `domain_km == C1_DOMAIN_KM` (full torus) reproduces
/// [`compute_temperature`] exactly; a cropped window passes its `window_km`.
pub fn compute_temperature_windowed(
    heightmap: &GridF32,
    ss: &SteinSteinParams,
    lat_deg: f32,
    domain_km: f32,
) -> GridF32 {
    compute_temperature_span(heightmap, ss, lat_deg, domain_km / 111.0)
}

/// Temperature field (°C) for an EXPLICIT latitudinal `span_deg` centred on
/// `centre_deg` (Finding 25 — TASK 2). `span_deg == domain_km / 111` reproduces
/// [`compute_temperature_windowed`] exactly.
pub fn compute_temperature_span(
    heightmap: &GridF32,
    ss: &SteinSteinParams,
    centre_deg: f32,
    span_deg: f32,
) -> GridF32 {
    let (w, h) = (heightmap.width, heightmap.height);
    let mut t = GridF32::new(w, h, 0.0);
    for j in 0..h {
        let t_sea = sea_level_temperature(row_latitude_span(j, h, centre_deg, span_deg));
        for i in 0..w {
            let n = heightmap.get(i as i32, j as i32);
            let cell = if n <= SEA_LEVEL_NORM {
                t_sea // ocean surface — no lapse under water
            } else {
                let alt_km = c1_altitude_norm_to_metres(n, ss).max(0.0) / 1000.0;
                t_sea - LAPSE_RATE_C_PER_KM * alt_km
            };
            t.set(i, j, cell);
        }
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lapse_cools_with_altitude() {
        let ss = SteinSteinParams::default();
        let mut hm = GridF32::new(4, 1, 0.5); // all ocean
        hm.set(1, 0, 0.9); // a high land cell
        let t = compute_temperature(&hm, &ss, 45.0);
        assert!(t.get(1, 0) < t.get(0, 0), "high cell must be colder than sea level");
    }

    #[test]
    fn latitudinal_gradient_present() {
        assert!(sea_level_temperature(0.0) > sea_level_temperature(45.0));
        assert!(sea_level_temperature(45.0) > sea_level_temperature(90.0));
    }

    /// Finding 27 — orientation invariant pinned on the DATA (a display flip cannot
    /// fool this): with `y = 0 = SOUTH`, row 0 is the LOWER latitude, so in the northern
    /// hemisphere row 0 is WARMER than the last (northern, polar) row. If this ever
    /// fails, `row_latitude_*` / the temperature field was inverted — NOT the renderer.
    #[test]
    fn row_zero_is_south_and_warmer_northern_hemisphere() {
        let ss = SteinSteinParams::default();
        let ny = 64usize;
        // row latitude: south edge (row 0) is the lower latitude.
        assert!(
            row_latitude_span(0, ny, 60.0, 40.0) < row_latitude_span(ny - 1, ny, 60.0, 40.0),
            "row 0 must be the SOUTH (lower-latitude) edge"
        );
        // temperature field: flat land → row 0 (south) warmer than the polar last row.
        let land = GridF32::new(8, ny, 0.62);
        let t = compute_temperature_span(&land, &ss, 60.0, 40.0);
        let row_mean = |j: usize| (0..8).map(|i| t.get(i as i32, j as i32)).sum::<f32>() / 8.0;
        assert!(
            row_mean(0) > row_mean(ny - 1),
            "northern hemisphere: south row (0) must be warmer than the north row"
        );
    }
}
