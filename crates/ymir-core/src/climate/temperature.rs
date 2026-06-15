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
use crate::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, C1_DOMAIN_KM};

/// Environmental (adiabatic) lapse rate — the real value, 6.5 °C per km.
pub const LAPSE_RATE_C_PER_KM: f32 = 6.5;
/// The C1 sea-level normalised value (Maillon 2 unified scale).
pub const SEA_LEVEL_NORM: f32 = 0.5;

/// Sea-level temperature (°C) at a latitude (degrees), anchored on the real
/// equator-pole gradient: `T = T_eq − (T_eq − T_pole)·sin²(lat)` (equator ~+27 °C,
/// pole ~−25 °C). At the 45° default this is a temperate sea-level value.
pub fn sea_level_temperature(lat_deg: f32) -> f32 {
    const T_EQ: f32 = 27.0;
    const T_POLE: f32 = -25.0;
    let s = lat_deg.to_radians().sin();
    T_EQ - (T_EQ - T_POLE) * s * s
}

/// Latitude (degrees) of grid row `j`: the domain spans `C1_DOMAIN_KM / 111`
/// ≈ 9.2° centred on `lat_deg` (the placement parameter). j=0 is the south edge.
pub fn row_latitude(j: usize, ny: usize, lat_deg: f32) -> f32 {
    let span = C1_DOMAIN_KM / 111.0;
    let frac = if ny > 1 { j as f32 / (ny - 1) as f32 } else { 0.5 };
    lat_deg + (frac - 0.5) * span
}

/// Temperature field (°C) for the relief heightmap at the given centre latitude.
pub fn compute_temperature(heightmap: &GridF32, ss: &SteinSteinParams, lat_deg: f32) -> GridF32 {
    let (w, h) = (heightmap.width, heightmap.height);
    let mut t = GridF32::new(w, h, 0.0);
    for j in 0..h {
        let t_sea = sea_level_temperature(row_latitude(j, h, lat_deg));
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
}
