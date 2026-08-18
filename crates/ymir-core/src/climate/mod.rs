//! Climate simulation and biome classification.
//!
//! #165: `c1_climate` derives temperature + precipitation from the finished C1
//! relief + a latitude placement parameter — a PURE derived computation (no
//! frozen field, no per-step state), re-runnable whenever the relief changes (a
//! future closure, a different seed). The viz latitude slider will control
//! `latitude_deg` and re-derive; it does not replace the field. Same discipline
//! as the coordinate contracts: one anchored computation, not a tuned output.
//! See `docs/design/c1_climate_design.md`.

pub mod biomes;
pub mod precipitation;
pub mod temperature;

use crate::grid::GridF32;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::production_upscale::{C1_DOMAIN_KM, c1_altitude_norm_to_metres};

use precipitation::PrecipParams;

/// Default continent placement: 45° (westerlies) — the temperate playable-region
/// target (Living Landz §2.2). The slider varies this.
pub const DEFAULT_LATITUDE_DEG: f32 = 45.0;

/// Derived climate of a relief heightmap.
pub struct ClimateResult {
    /// Temperature (°C): latitudinal sea-level gradient − adiabatic lapse.
    pub temperature: GridF32,
    /// Precipitation (internal units; see `precipitation::precip_mm_per_year`).
    pub precipitation: GridF32,
}

/// Whittaker biome map from a climate result (the chain `c1_climate → c1_biomes`).
/// Derived & re-runnable: recompute when the relief/climate change. `heightmap`
/// marks ocean cells. Row-major `Vec<Biome>`.
pub fn c1_biomes(heightmap: &GridF32, climate: &ClimateResult) -> Vec<biomes::Biome> {
    biomes::compute_biomes(heightmap, &climate.temperature, &climate.precipitation)
}

/// Derive (temperature, precipitation) from the C1 relief at a centre latitude.
/// Pure function of `(heightmap, latitude_deg)` (+ the fixed contracts via `ss`)
/// — re-run it when the relief changes. `heightmap` is the upscaled C1 product
/// (normalised, sea = 0.5); `ss` supplies the vertical scale for norm→metres.
pub fn c1_climate(
    heightmap: &GridF32,
    ss: &SteinSteinParams,
    latitude_deg: f32,
    params: &PrecipParams,
) -> ClimateResult {
    c1_climate_windowed(heightmap, ss, latitude_deg, params, C1_DOMAIN_KM)
}

/// [`c1_climate`] for a grid spanning `window_km` (a cropped playable window).
/// The horizontal scale (`km_per_cell = window_km / width`) drives precipitation
/// orographic distance, and `window_km` sets the temperature latitudinal span —
/// so a zoomed window reads its OWN metric scale, not the full torus.
/// `window_km == C1_DOMAIN_KM` reproduces [`c1_climate`] exactly.
pub fn c1_climate_windowed(
    heightmap: &GridF32,
    ss: &SteinSteinParams,
    latitude_deg: f32,
    params: &PrecipParams,
    window_km: f32,
) -> ClimateResult {
    let temperature =
        temperature::compute_temperature_windowed(heightmap, ss, latitude_deg, window_km);
    let km_per_cell = window_km / heightmap.width as f32;
    let precipitation = precipitation::compute_precipitation(
        heightmap,
        &temperature,
        latitude_deg,
        km_per_cell,
        |n| c1_altitude_norm_to_metres(n, ss),
        params,
    );
    ClimateResult { temperature, precipitation }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c1_climate_runs_and_is_finite() {
        let ss = SteinSteinParams::default();
        let mut hm = GridF32::new(32, 32, 0.4); // mostly ocean
        for j in 8..24usize {
            for i in 8..24usize {
                hm.set(i, j, 0.7); // a land block
            }
        }
        let c = c1_climate(&hm, &ss, DEFAULT_LATITUDE_DEG, &PrecipParams::default());
        assert!(c.temperature.data.iter().all(|v| v.is_finite()));
        assert!(c.precipitation.data.iter().all(|v| v.is_finite() && *v >= 0.0));
    }
}
