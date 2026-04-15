//! Output file writers for generated terrain data.
//!
//! Each pipeline phase saves its result as a PNG alongside a metadata.json
//! file that records all parameters for reproducibility.

use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::erosion::hydraulic::{ErosionConfig, ErosionResult, ErosionStats};
use crate::grid::GridF32;
use crate::lakes::detection::{Lake, LakeConfig, LakeResult};
use crate::tectonics::isostasy::{IsostasyConfig, IsostasyResult};
use crate::tectonics::plates::PlateConfig;
use crate::terrain::flow::{FlowConfig, FlowResult, RiverConfig, RiverNetwork};
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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flow: Option<FlowMetadata>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lakes: Option<LakeMetadata>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub centering_shift: Option<(i32, i32)>,
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowMetadata {
    pub sea_level: f32,
    pub num_basins: u32,
    pub grid_width: usize,
    pub grid_height: usize,
    pub river_config: RiverConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LakeMetadata {
    pub config: LakeConfig,
    pub num_lakes: usize,
    pub total_lake_area: usize,
    pub grid_width: usize,
    pub grid_height: usize,
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
            flow: None,
            lakes: None,
            centering_shift: None,
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
            flow: None,
            lakes: None,
            centering_shift: None,
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
        thickness.save_raw(&self.dir.join("01_thickness.raw"))?;
        thickness.save_png_u16(&self.dir.join("01_thickness.png"))
    }

    /// Save the altitude heightmap after isostasy.
    pub fn save_altitude(
        &mut self,
        result: &IsostasyResult,
        config: &IsostasyConfig,
    ) -> Result<(), String> {
        result.heightmap.save_raw(&self.dir.join("02_altitude.raw"))?;
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
        heightmap.save_raw(&self.dir.join("03_upscaled.raw"))?;
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
        let n = self.metadata.grid_size;
        GridF32::load_raw_or_png(
            &self.dir.join("01_thickness.raw"),
            &self.dir.join("01_thickness.png"),
            n,
            n,
        )
    }

    /// Load the altitude heightmap from a previously saved export.
    pub fn load_altitude(&self) -> Result<GridF32, String> {
        let n = self.metadata.grid_size;
        GridF32::load_raw_or_png(
            &self.dir.join("02_altitude.raw"),
            &self.dir.join("02_altitude.png"),
            n,
            n,
        )
    }

    /// Load the upscaled heightmap from a previously saved export.
    pub fn load_upscaled(&self) -> Result<GridF32, String> {
        let size = self
            .metadata
            .upscale
            .as_ref()
            .map(|u| u.target_size)
            .unwrap_or(self.metadata.grid_size);
        GridF32::load_raw_or_png(
            &self.dir.join("03_upscaled.raw"),
            &self.dir.join("03_upscaled.png"),
            size,
            size,
        )
    }

    /// Save the eroded heightmap and sediment map.
    pub fn save_eroded(
        &mut self,
        result: &ErosionResult,
        config: &ErosionConfig,
    ) -> Result<(), String> {
        result.heightmap.save_raw(&self.dir.join("04_eroded.raw"))?;
        result.heightmap.save_png_u16(&self.dir.join("04_eroded.png"))?;
        result.sediment.save_raw(&self.dir.join("04_sediment.raw"))?;
        result.sediment.save_png_u16(&self.dir.join("04_sediment.png"))?;
        self.metadata.erosion =
            Some(ErosionMetadataEntry { config: config.clone(), stats: result.stats.clone() });
        self.save_metadata()
    }

    /// Load the eroded heightmap from a previously saved export.
    pub fn load_eroded(&self) -> Result<GridF32, String> {
        let size = self
            .metadata
            .upscale
            .as_ref()
            .map(|u| u.target_size)
            .unwrap_or(self.metadata.grid_size);
        GridF32::load_raw_or_png(
            &self.dir.join("04_eroded.raw"),
            &self.dir.join("04_eroded.png"),
            size,
            size,
        )
    }

    /// Load the sediment map from a previously saved export.
    pub fn load_sediment(&self) -> Result<GridF32, String> {
        let size = self
            .metadata
            .upscale
            .as_ref()
            .map(|u| u.target_size)
            .unwrap_or(self.metadata.grid_size);
        GridF32::load_raw_or_png(
            &self.dir.join("04_sediment.raw"),
            &self.dir.join("04_sediment.png"),
            size,
            size,
        )
    }

    /// Save all flow computation results (lossless).
    pub fn save_flow(
        &mut self,
        result: &FlowResult,
        config: &FlowConfig,
        river_config: &RiverConfig,
        rivers: Option<&RiverNetwork>,
    ) -> Result<(), String> {
        let w = result.accumulation.width;
        let h = result.accumulation.height;

        result.filled.save_raw(&self.dir.join("05_filled.raw"))?;
        result.filled.save_png_u16(&self.dir.join("05_filled.png"))?;
        save_raw_f32(&self.dir.join("05_flow_accumulation.raw"), &result.accumulation.data)?;
        save_raw_u8(&self.dir.join("05_flow_direction.raw"), &result.direction)?;
        save_raw_u32(&self.dir.join("05_basins.raw"), &result.basins)?;

        // Visual-only: log-scaled accumulation PNG
        let max_flow = result.accumulation.data.iter().cloned().fold(0.0f32, f32::max);
        if max_flow > 0.0 {
            let log_max = (1.0 + max_flow).ln();
            let viz = GridF32::from_vec(
                w,
                h,
                result.accumulation.data.iter().map(|&v| (1.0 + v).ln() / log_max).collect(),
            );
            let _ = viz.save_png_u16(&self.dir.join("05_flow_accumulation_viz.png"));
        }

        if let Some(network) = rivers {
            let json =
                serde_json::to_string_pretty(network).map_err(|e| format!("JSON error: {e}"))?;
            fs::write(self.dir.join("05_rivers.json"), json)
                .map_err(|e| format!("Write error: {e}"))?;
        }

        self.metadata.flow = Some(FlowMetadata {
            sea_level: config.sea_level,
            num_basins: result.num_basins,
            grid_width: w,
            grid_height: h,
            river_config: river_config.clone(),
        });
        self.save_metadata()
    }

    /// Load flow computation results.
    pub fn load_flow(&self) -> Result<(FlowResult, RiverConfig), String> {
        let meta = self.metadata.flow.as_ref().ok_or("No flow metadata")?;
        let w = meta.grid_width;
        let h = meta.grid_height;
        let n = w * h;

        let filled = GridF32::load_raw_or_png(
            &self.dir.join("05_filled.raw"),
            &self.dir.join("05_filled.png"),
            w,
            h,
        )?;
        let acc_data = load_raw_f32(&self.dir.join("05_flow_accumulation.raw"), n)?;
        let direction = load_raw_u8(&self.dir.join("05_flow_direction.raw"))?;
        let basins = load_raw_u32(&self.dir.join("05_basins.raw"), n)?;

        if direction.len() != n {
            return Err(format!("Direction size mismatch: expected {n}, got {}", direction.len()));
        }
        if basins.len() != n {
            return Err(format!("Basins size mismatch: expected {n}, got {}", basins.len()));
        }

        let accumulation = GridF32::from_vec(w, h, acc_data);
        let result =
            FlowResult { filled, direction, accumulation, basins, num_basins: meta.num_basins };

        Ok((result, meta.river_config.clone()))
    }

    /// Load the pre-extracted river network.
    pub fn load_rivers(&self) -> Result<RiverNetwork, String> {
        let json = fs::read_to_string(self.dir.join("05_rivers.json"))
            .map_err(|e| format!("Read error: {e}"))?;
        serde_json::from_str(&json).map_err(|e| format!("JSON parse error: {e}"))
    }

    /// Save lake detection results.
    pub fn save_lakes(&mut self, result: &LakeResult, config: &LakeConfig) -> Result<(), String> {
        save_raw_u32(&self.dir.join("06_lakes.raw"), &result.lake_map)?;

        let json =
            serde_json::to_string_pretty(&result.lakes).map_err(|e| format!("JSON error: {e}"))?;
        fs::write(self.dir.join("06_lakes.json"), json).map_err(|e| format!("Write error: {e}"))?;

        let total_area: usize = result.lakes.iter().map(|l| l.area).sum();
        self.metadata.lakes = Some(LakeMetadata {
            config: config.clone(),
            num_lakes: result.lakes.len(),
            total_lake_area: total_area,
            grid_width: result.width,
            grid_height: result.height,
        });
        self.save_metadata()
    }

    /// Load lake detection results.
    pub fn load_lakes(&self) -> Result<LakeResult, String> {
        let meta = self.metadata.lakes.as_ref().ok_or("No lake metadata")?;
        let w = meta.grid_width;
        let h = meta.grid_height;
        let n = w * h;

        let lake_map = load_raw_u32(&self.dir.join("06_lakes.raw"), n)?;

        let json = fs::read_to_string(self.dir.join("06_lakes.json"))
            .map_err(|e| format!("Read error: {e}"))?;
        let lakes: Vec<Lake> =
            serde_json::from_str(&json).map_err(|e| format!("JSON parse error: {e}"))?;

        Ok(LakeResult { lake_map, lakes, width: w, height: h })
    }
}

// ── Raw binary helpers ──────────────────────────────────────────────────

fn save_raw_f32(path: &Path, data: &[f32]) -> Result<(), String> {
    let bytes: Vec<u8> = data.iter().flat_map(|f| f.to_le_bytes()).collect();
    fs::write(path, bytes).map_err(|e| format!("Write error: {e}"))
}

fn load_raw_f32(path: &Path, expected_len: usize) -> Result<Vec<f32>, String> {
    let bytes = fs::read(path).map_err(|e| format!("Read error: {e}"))?;
    if bytes.len() != expected_len * 4 {
        return Err(format!(
            "Size mismatch: expected {} bytes, got {}",
            expected_len * 4,
            bytes.len()
        ));
    }
    Ok(bytes.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
}

fn save_raw_u8(path: &Path, data: &[u8]) -> Result<(), String> {
    fs::write(path, data).map_err(|e| format!("Write error: {e}"))
}

fn load_raw_u8(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|e| format!("Read error: {e}"))
}

fn save_raw_u32(path: &Path, data: &[u32]) -> Result<(), String> {
    let bytes: Vec<u8> = data.iter().flat_map(|v| v.to_le_bytes()).collect();
    fs::write(path, bytes).map_err(|e| format!("Write error: {e}"))
}

fn load_raw_u32(path: &Path, expected_len: usize) -> Result<Vec<u32>, String> {
    let bytes = fs::read(path).map_err(|e| format!("Read error: {e}"))?;
    if bytes.len() != expected_len * 4 {
        return Err(format!(
            "Size mismatch: expected {} bytes, got {}",
            expected_len * 4,
            bytes.len()
        ));
    }
    Ok(bytes.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
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
