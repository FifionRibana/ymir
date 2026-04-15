//! Output file writers for generated terrain data.
//!
//! Each pipeline phase saves its result as a PNG alongside a metadata.json
//! file that records all parameters for reproducibility.

use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::erosion::hydraulic::{ErosionConfig, ErosionResult, ErosionStats};
use crate::grid::GridF32;
use crate::tectonics::isostasy::{IsostasyConfig, IsostasyResult};
use crate::tectonics::plates::PlateConfig;
use crate::terrain::upscale::FbmUpscaleConfig;

// pub mod raw;     // Native binary format (u16/f32 raw + JSON metadata)
// pub mod png;     // PNG export for compatibility and debugging

// ── Metadata ─────────────────────────────────────────────────────────────

/// Metadata for a complete generation run, saved alongside the output files.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct PipelineMetadata {
    /// Ymir version (for compatibility checking).
    pub version: String,
    /// Master seed.
    pub seed: u64,
    /// Grid resolution.
    pub grid_size: usize,
    /// Meters per pixel.
    pub meters_per_pixel: f32,
    /// Plate generation config.
    pub plates: PlateConfig,
    /// Isostasy config and results.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub isostasy: Option<IsostasyMetadata>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub upscale: Option<FbmUpscaleConfig>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub erosion: Option<ErosionMetadataEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsostasyMetadata {
    pub config: IsostasyConfig,
    pub sea_level_normalized: f32,
    pub peak_altitude_m: f32,
    pub max_depth_m: f32,
    pub land_ratio: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErosionMetadataEntry {
    pub config: ErosionConfig,
    pub stats: ErosionStats,
}

impl Default for PipelineMetadata {
    fn default() -> Self {
        Self {
            version: env!("CARGO_PKG_VERSION").to_string(),
            seed: 0,
            grid_size: 128,
            meters_per_pixel: 40.0,
            plates: PlateConfig::default(),
            isostasy: None,
            upscale: None,
            erosion: None,
        }
    }
}

impl PipelineMetadata {
    pub fn new(seed: u64, grid_size: usize, plates: &PlateConfig) -> Self {
        Self {
            version: env!("CARGO_PKG_VERSION").to_string(),
            seed,
            grid_size,
            meters_per_pixel: 40.0,
            plates: plates.clone(),
            isostasy: None,
            upscale: None,
            erosion: None,
        }
    }
}

// ── Pipeline export ──────────────────────────────────────────────────────

/// Save pipeline state to a directory.
///
/// Creates the directory if needed. Each phase writes its output file and
/// updates metadata.json.
pub struct PipelineExport {
    pub dir: PathBuf,
    pub metadata: PipelineMetadata,
}

impl PipelineExport {
    /// Create a new export directory for the given seed and grid size.
    pub fn new(output_root: &Path, seed: u64, grid_size: usize, plates: &PlateConfig) -> Self {
        let dir_name = format!("seed{}_{}", seed, grid_size);
        let dir = output_root.join(dir_name);
        fs::create_dir_all(&dir).ok();

        Self { dir, metadata: PipelineMetadata::new(seed, grid_size, plates) }
    }

    /// Save the crustal thickness field after tectonics.
    pub fn save_thickness(&self, thickness: &GridF32) -> Result<(), String> {
        thickness.save_png_u16(&self.dir.join("01_thickness.png"))
    }

    /// Save the altitude heightmap after isostasy.
    pub fn save_altitude(
        &mut self,
        result: &IsostasyResult,
        config: &IsostasyConfig,
    ) -> Result<(), String> {
        result.heightmap.save_png_u16(&self.dir.join("02_altitude.png"))?;

        self.metadata.isostasy = Some(IsostasyMetadata {
            config: config.clone(),
            sea_level_normalized: result.sea_level_normalized,
            peak_altitude_m: result.peak_altitude_m,
            max_depth_m: result.max_depth_m,
            land_ratio: result.land_ratio,
        });

        self.save_metadata()
    }

    /// Write metadata.json to disk.
    pub fn save_metadata(&self) -> Result<(), String> {
        let json =
            serde_json::to_string_pretty(&self.metadata).map_err(|e| format!("JSON error: {e}"))?;
        fs::write(self.dir.join("metadata.json"), json).map_err(|e| format!("Write error: {e}"))?;
        Ok(())
    }

    /// Save the upscaled heightmap.
    pub fn save_upscaled(
        &mut self,
        heightmap: &GridF32,
        config: &FbmUpscaleConfig,
    ) -> Result<(), String> {
        heightmap.save_png_u16(&self.dir.join("03_upscaled.png"))?;
        self.metadata.upscale = Some(config.clone());
        self.save_metadata()
    }

    /// Load an existing pipeline export from a directory.
    pub fn load(dir: &Path) -> Result<Self, String> {
        let json = fs::read_to_string(dir.join("metadata.json"))
            .map_err(|e| format!("Read error: {e}"))?;
        let metadata: PipelineMetadata =
            serde_json::from_str(&json).map_err(|e| format!("JSON parse error: {e}"))?;
        Ok(Self { dir: dir.to_path_buf(), metadata })
    }

    /// Load the thickness field from a previously saved export.
    pub fn load_thickness(&self) -> Result<GridF32, String> {
        GridF32::load_png(&self.dir.join("01_thickness.png"))
    }

    /// Load the altitude heightmap from a previously saved export.
    pub fn load_altitude(&self) -> Result<GridF32, String> {
        GridF32::load_png(&self.dir.join("02_altitude.png"))
    }

    /// Load the upscaled heightmap from a previously saved export.
    pub fn load_upscaled(&self) -> Result<GridF32, String> {
        GridF32::load_png(&self.dir.join("03_upscaled.png"))
    }

    /// Save the eroded heightmap and sediment map.
    pub fn save_eroded(
        &mut self,
        result: &ErosionResult,
        config: &ErosionConfig,
    ) -> Result<(), String> {
        result.heightmap.save_png_u16(&self.dir.join("04_eroded.png"))?;
        result.sediment.save_png_u16(&self.dir.join("04_sediment.png"))?;
        self.metadata.erosion =
            Some(ErosionMetadataEntry { config: config.clone(), stats: result.stats.clone() });
        self.save_metadata()
    }

    /// Load the eroded heightmap from a previously saved export.
    pub fn load_eroded(&self) -> Result<GridF32, String> {
        GridF32::load_png(&self.dir.join("04_eroded.png"))
    }

    /// Load the sediment map from a previously saved export.
    pub fn load_sediment(&self) -> Result<GridF32, String> {
        GridF32::load_png(&self.dir.join("04_sediment.png"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deserialize_legacy_metadata_without_upscale() {
        let json = r#"{
            "version": "0.1.0",
            "seed": 42,
            "grid_size": 128,
            "meters_per_pixel": 40.0,
            "plates": {
                "num_plates": 8,
                "continental_ratio": 0.35,
                "velocity_min": 0.5,
                "velocity_max": 2.5,
                "grid_size": 128,
                "boundary_smoothing_sigma": 2.0
            },
            "isostasy": {
                "config": {
                    "rho_crust": 2750.0,
                    "rho_mantle": 3300.0,
                    "rho_water": 1025.0,
                    "max_elevation_m": 4000.0,
                    "max_depth_m": 500.0,
                    "sea_level_fraction": 0.4,
                    "altitude_smoothing_sigma": 2.0
                },
                "sea_level_normalized": 0.11,
                "peak_altitude_m": 2400.0,
                "max_depth_m": 200.0,
                "land_ratio": 0.35
            }
        }"#;

        let meta: PipelineMetadata = serde_json::from_str(json).unwrap();
        assert_eq!(meta.seed, 42);
        assert!(meta.upscale.is_none());
    }

    #[test]
    fn deserialize_upscale_config_without_domain_warp() {
        // An upscale config saved before domain warp fields were added
        let json = r#"{
            "target_size": 1024,
            "octaves": 7,
            "lacunarity": 2.0,
            "persistence": 0.5,
            "amplitude_base": 0.08,
            "amplitude_slope_factor": 3.0,
            "max_anisotropy": 3.0,
            "submarine_damping": 0.3,
            "base_frequency": 1.0
        }"#;

        let config: FbmUpscaleConfig = serde_json::from_str(json).unwrap();
        assert_eq!(config.target_size, 1024);
        // New fields should have their defaults
        assert!((config.domain_warp_strength - 0.0).abs() < 1e-10);
        assert!((config.domain_warp_frequency - 0.5).abs() < 1e-10);
        assert_eq!(config.domain_warp_octaves, 3);
    }
}
