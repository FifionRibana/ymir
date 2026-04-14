//! Output file writers for generated terrain data.
//!
//! Each pipeline phase saves its result as a PNG alongside a metadata.json
//! file that records all parameters for reproducibility.

use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::tectonics::isostasy::{IsostasyConfig, IsostasyResult};
use crate::tectonics::plates::PlateConfig;
use crate::terrain::upscale::FbmUpscaleConfig;

// pub mod raw;     // Native binary format (u16/f32 raw + JSON metadata)
// pub mod png;     // PNG export for compatibility and debugging

// ── Metadata ─────────────────────────────────────────────────────────────

/// Metadata for a complete generation run, saved alongside the output files.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsostasyMetadata {
    pub config: IsostasyConfig,
    pub sea_level_normalized: f32,
    pub peak_altitude_m: f32,
    pub max_depth_m: f32,
    pub land_ratio: f32,
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
}
