//! Top-level generation configuration.
//!
//! [`GenerationConfig`] groups all parameters for a full pipeline run.
//! Phase-specific configs are nested within it. The struct derives
//! `Serialize`/`Deserialize` so it can be saved alongside the generated
//! output for full reproducibility: seed + config = identical world.

use crate::seed::WorldSeed;
use serde::{Deserialize, Serialize};

/// Master configuration for a complete continent generation run.
///
/// Saved to `metadata.json` alongside the exported terrain data so that
/// any generated continent can be reproduced exactly from its config file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenerationConfig {
    /// Name of the continent (used for output directory naming).
    pub name: String,

    /// Master seed — deterministic source for all random operations.
    pub seed: WorldSeed,

    /// Meters per pixel on the erosion grid. Controls the tradeoff between
    /// terrain detail and generation speed. Recommended: 35–40 for production,
    /// 70–80 for fast iteration.
    pub meters_per_pixel: f32,

    /// Target continent size in kilometers (approximate side length).
    /// Combined with `meters_per_pixel`, this determines the erosion grid
    /// resolution: `grid_size = (continent_km * 1000) / meters_per_pixel`.
    pub continent_size_km: f32,

    /// Maximum land elevation in meters. Collision zones may reach this value.
    pub max_elevation_m: f32,

    /// Maximum ocean depth in meters (positive value, represents depth below
    /// sea level). Used for bathymetry normalization.
    pub max_depth_m: f32,
}

impl Default for GenerationConfig {
    fn default() -> Self {
        Self {
            name: "unnamed".to_string(),
            seed: WorldSeed::new(42),
            meters_per_pixel: 40.0,
            continent_size_km: 300.0,
            max_elevation_m: 4000.0,
            max_depth_m: 500.0,
        }
    }
}

impl GenerationConfig {
    /// Compute the erosion grid side length from continent size and resolution.
    pub fn erosion_grid_size(&self) -> usize {
        ((self.continent_size_km * 1000.0) / self.meters_per_pixel).ceil() as usize
    }

    /// Compute the tectonic grid size (typically 1/8 to 1/16 of erosion grid).
    pub fn tectonic_grid_size(&self) -> usize {
        let erosion = self.erosion_grid_size();
        // Round to nearest power of 2 for FFT-friendly sizing
        let target = erosion / 16;
        target.next_power_of_two().max(64)
    }
}
