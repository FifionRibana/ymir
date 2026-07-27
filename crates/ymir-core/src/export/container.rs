//! v1 `.ymir` delivery container — a self-describing continent directory.
//!
//! A continent is exported as a **`<name>.ymir/` directory** holding a single
//! auto-descriptive `manifest.json` plus one file per data layer (binary raster
//! or vector `.json`/`.geojson`). No reader ever has to guess a layout or a
//! dtype — the manifest names the file, `kind`, `dtype`, `endianness` and grid
//! dimensions of every layer, and carries a `present` flag so a missing layer
//! is tolerated rather than fatal.
//!
//! This is the *delivery* format Living Landz consumes; it is **distinct** from
//! the per-step debug export in [`super`] (`PipelineExport`) and from the
//! content-addressed [`crate::cache`]. All three share the ONE binary codec in
//! [`super::raw`] (LE, no header) — this module adds no second serializer.
//!
//! # Raster orientation invariant (documented ONCE, here — the source of flip
//! bugs)
//!
//! Every raster layer is **row-major** with **`y = 0` = the SOUTH edge**: the
//! first row of `width` values is the southernmost, and `x` increases eastward.
//! Writers and readers must both honour this; nothing in the format re-encodes
//! orientation, so a producer that stores north-up MUST flip before calling
//! `add_raster_*`.
//!
//! # Versioning
//!
//! `format_version` is semver. A reader MUST refuse an unknown *major*.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;

use super::raw;

/// The exchange-format version this writer emits.
pub const FORMAT_VERSION: &str = "1.0.0";

// ── Manifest schema (serde, matches docs/WP0_exchange_format.md v1) ────────

/// Top-level `manifest.json` document.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    /// Exchange-format semver (see [`FORMAT_VERSION`]).
    pub format_version: String,
    /// Producing Ymir crate version (`CARGO_PKG_VERSION`).
    pub ymir_version: String,
    /// Physical + provenance description of the exported continent.
    pub continent: Continent,
    /// Multi-continent geodesy — OMITTED in v1 (single continent). Added by
    /// world-edit later; skipped from the JSON when absent.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub geodesy: Option<serde_json::Value>,
    /// Every known layer, each with a `present` flag.
    pub layers: Vec<Layer>,
}

/// Grid dimensions in cells.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Grid {
    pub width: usize,
    pub height: usize,
}

/// `c1_altitude_norm_to_metres` contract — how normalized altitude maps to
/// metres (filled from [`SteinSteinParams`]).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct VerticalScale {
    /// Non-dim altitude half-range (`asymptotic_depth_m / depth_scale_m`).
    pub altitude_norm_half_range: f64,
    /// Depth normalisation `[m / non-dim unit]` (`SteinSteinParams::depth_scale_m`).
    pub depth_scale_m: f64,
    /// Normalized value of sea level (fixed `0.5` by the C1 normalisation).
    pub sea_level_norm: f64,
}

/// Physical + provenance metadata for the continent.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Continent {
    pub name: String,
    pub seed: u64,
    pub grid: Grid,
    /// Physical size of the exported window (km).
    pub window_km: f64,
    /// `= window_km / grid.width` (≈ one hex).
    pub km_per_cell: f64,
    /// Source C1 torus size (km) — context.
    pub tectonic_domain_km: f64,
    /// `0..1` window origin in the torus.
    pub window_offset_in_torus: [f64; 2],
    pub vertical_scale: VerticalScale,
    pub sea_level_m: f64,
    pub max_elevation_m: f64,
    pub max_depth_m: f64,
}

/// One layer entry. A single flat struct covers both raster and vector layers;
/// attributes that do not apply to a given layer are `None` and skipped from
/// the JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Layer {
    pub id: String,
    pub file: String,
    /// `"raster"` or `"vector"`.
    pub kind: String,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub dtype: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub unit: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub encoding: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub min_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub max_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub width: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub height: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub endianness: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub geometry: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub level_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub slope_threshold_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub semantics: Option<String>,
    /// A `false` layer's file may be absent; a reader MUST tolerate it.
    pub present: bool,
}

// ── ContinentMeta ──────────────────────────────────────────────────────────

/// Caller-supplied physical description of a continent, consumed by
/// [`ContinentWriter::new`] to build the manifest `continent` block. The
/// `km_per_cell` and `vertical_scale` fields are *derived* here, not passed.
#[derive(Debug, Clone)]
pub struct ContinentMeta {
    pub name: String,
    pub seed: u64,
    pub grid: Grid,
    /// Physical size of the exported window (km).
    pub window_km: f64,
    /// Source C1 torus size (km).
    pub tectonic_domain_km: f64,
    /// `0..1` window origin in the torus.
    pub window_offset_in_torus: [f64; 2],
    /// Vertical scale is derived from these Stein-Stein params.
    pub stein_stein: SteinSteinParams,
    pub sea_level_m: f64,
    pub max_elevation_m: f64,
    pub max_depth_m: f64,
}

impl ContinentMeta {
    /// Derive the manifest `continent` block, computing
    /// `km_per_cell = window_km / grid.width` and the `vertical_scale`.
    fn to_continent(&self) -> Continent {
        let km_per_cell =
            if self.grid.width > 0 { self.window_km / self.grid.width as f64 } else { 0.0 };
        let ss = &self.stein_stein;
        let altitude_norm_half_range =
            if ss.depth_scale_m != 0.0 { ss.asymptotic_depth_m / ss.depth_scale_m } else { 0.0 };
        Continent {
            name: self.name.clone(),
            seed: self.seed,
            grid: self.grid,
            window_km: self.window_km,
            km_per_cell,
            tectonic_domain_km: self.tectonic_domain_km,
            window_offset_in_torus: self.window_offset_in_torus,
            vertical_scale: VerticalScale {
                altitude_norm_half_range,
                depth_scale_m: ss.depth_scale_m,
                // The C1 normalisation always pins sea level to 0.5.
                sea_level_norm: 0.5,
            },
            sea_level_m: self.sea_level_m,
            max_elevation_m: self.max_elevation_m,
            max_depth_m: self.max_depth_m,
        }
    }
}

// ── ContinentWriter ────────────────────────────────────────────────────────

/// Builds a `<name>.ymir/` directory: writes layer files through
/// [`super::raw`], flips each layer's `present` flag, and finally serialises
/// `manifest.json`.
pub struct ContinentWriter {
    dir: PathBuf,
    manifest: Manifest,
}

impl ContinentWriter {
    /// Create the export directory and stage a manifest with every known layer
    /// declared `present = false`. Raster layers are pre-filled with the grid
    /// dimensions; the `present` flag (and the file on disk) is set by the
    /// matching `add_*` call.
    pub fn new(dir: &Path, meta: ContinentMeta) -> Result<Self, String> {
        std::fs::create_dir_all(dir).map_err(|e| format!("Create dir error: {e}"))?;
        let continent = meta.to_continent();
        let (w, h) = (continent.grid.width, continent.grid.height);
        let manifest = Manifest {
            format_version: FORMAT_VERSION.to_string(),
            ymir_version: env!("CARGO_PKG_VERSION").to_string(),
            continent,
            geodesy: None,
            layers: default_layers(w, h),
        };
        Ok(Self { dir: dir.to_path_buf(), manifest })
    }

    /// Immutable view of the staged manifest (before/after `add_*`).
    pub fn manifest(&self) -> &Manifest {
        &self.manifest
    }

    /// Look up a mutable layer by `id`, erroring if it is not a known layer.
    fn layer_mut(&mut self, id: &str) -> Result<&mut Layer, String> {
        self.manifest
            .layers
            .iter_mut()
            .find(|l| l.id == id)
            .ok_or_else(|| format!("Unknown layer id: {id}"))
    }

    /// Grid cell count the manifest declares (rasters must match it).
    fn cells(&self) -> usize {
        self.manifest.continent.grid.width * self.manifest.continent.grid.height
    }

    /// Common raster bookkeeping: validate length, mark present, fill
    /// dtype/endianness/dims. The file itself is written by the caller.
    fn set_raster_layer(&mut self, id: &str, len: usize, dtype: &str) -> Result<PathBuf, String> {
        let cells = self.cells();
        if len != cells {
            return Err(format!("Layer '{id}': data length {len} != grid cells {cells}"));
        }
        let (w, h) = (self.manifest.continent.grid.width, self.manifest.continent.grid.height);
        let dir = self.dir.clone();
        let layer = self.layer_mut(id)?;
        if layer.kind != "raster" {
            return Err(format!("Layer '{id}' is not a raster layer"));
        }
        layer.present = true;
        layer.dtype = Some(dtype.to_string());
        layer.endianness = Some("le".to_string());
        layer.width = Some(w);
        layer.height = Some(h);
        Ok(dir.join(&layer.file))
    }

    /// Write a `u16` raster (row-major, `y=0`=south). See module invariant.
    pub fn add_raster_u16(&mut self, id: &str, data: &[u16]) -> Result<(), String> {
        let path = self.set_raster_layer(id, data.len(), "u16")?;
        raw::save_u16(&path, data)
    }

    /// Write an `i16` raster (row-major, `y=0`=south).
    pub fn add_raster_i16(&mut self, id: &str, data: &[i16]) -> Result<(), String> {
        let path = self.set_raster_layer(id, data.len(), "i16")?;
        raw::save_i16(&path, data)
    }

    /// Write a `u32` raster (row-major, `y=0`=south).
    pub fn add_raster_u32(&mut self, id: &str, data: &[u32]) -> Result<(), String> {
        let path = self.set_raster_layer(id, data.len(), "u32")?;
        raw::save_u32(&path, data)
    }

    /// Write a `u8` raster (row-major, `y=0`=south).
    pub fn add_raster_u8(&mut self, id: &str, data: &[u8]) -> Result<(), String> {
        let path = self.set_raster_layer(id, data.len(), "u8")?;
        raw::save_u8(&path, data)
    }

    /// Write an `f32` raster (row-major, `y=0`=south).
    pub fn add_raster_f32(&mut self, id: &str, data: &[f32]) -> Result<(), String> {
        let path = self.set_raster_layer(id, data.len(), "f32")?;
        raw::save_f32(&path, data)
    }

    /// Copy/write a pre-serialised vector file (geojson/json) verbatim and mark
    /// the layer present. `filename` overrides the layer's default file name so
    /// the caller controls the on-disk name; no geojson dependency is pulled in.
    pub fn add_vector_file(
        &mut self,
        id: &str,
        filename: &str,
        bytes: &[u8],
    ) -> Result<(), String> {
        let dir = self.dir.clone();
        let layer = self.layer_mut(id)?;
        if layer.kind != "vector" {
            return Err(format!("Layer '{id}' is not a vector layer"));
        }
        layer.file = filename.to_string();
        layer.present = true;
        let path = dir.join(filename);
        std::fs::write(&path, bytes).map_err(|e| format!("Write error: {e}"))
    }

    /// Serialise `manifest.json` (pretty) into the export directory.
    pub fn finish(self) -> Result<PathBuf, String> {
        let json =
            serde_json::to_string_pretty(&self.manifest).map_err(|e| format!("JSON error: {e}"))?;
        let path = self.dir.join("manifest.json");
        std::fs::write(&path, json).map_err(|e| format!("Write error: {e}"))?;
        Ok(path)
    }
}

/// The canonical v1 layer table — every layer the format knows about, with its
/// intrinsic (static) attributes. Raster dims are pre-filled with the grid
/// size; `present` starts `false` for all.
fn default_layers(w: usize, h: usize) -> Vec<Layer> {
    // Helpers keep the table terse and self-documenting.
    let raster = |id: &str, file: &str, dtype: &str| Layer {
        id: id.to_string(),
        file: file.to_string(),
        kind: "raster".to_string(),
        dtype: Some(dtype.to_string()),
        unit: None,
        encoding: None,
        min_m: None,
        max_m: None,
        width: Some(w),
        height: Some(h),
        endianness: Some("le".to_string()),
        geometry: None,
        level_m: None,
        slope_threshold_deg: None,
        semantics: None,
        present: false,
    };
    let vector = |id: &str, file: &str| Layer {
        id: id.to_string(),
        file: file.to_string(),
        kind: "vector".to_string(),
        dtype: None,
        unit: None,
        encoding: None,
        min_m: None,
        max_m: None,
        width: None,
        height: None,
        endianness: None,
        geometry: None,
        level_m: None,
        slope_threshold_deg: None,
        semantics: None,
        present: false,
    };

    let mut height = raster("height", "height.u16", "u16");
    height.unit = Some("meter".to_string());
    height.encoding = Some("linear".to_string());

    let mut coastline = vector("coastline", "coastline.geojson");
    coastline.geometry = Some("multilinestring".to_string());
    coastline.level_m = Some(0.0);

    let mut cliffs = vector("cliffs", "cliffs.geojson");
    cliffs.geometry = Some("multilinestring".to_string());
    cliffs.slope_threshold_deg = Some(45.0);

    let flow_accumulation = raster("flow_accumulation", "flow_accumulation.f32", "f32");

    let rivers = vector("rivers", "rivers.json");
    let lakes = vector("lakes", "lakes.json");

    let lake_mask = raster("lake_mask", "lake_mask.u32", "u32");

    let mut biome = raster("biome", "biome.u8", "u8");
    biome.endianness = None; // single-byte: endianness is meaningless.
    biome.semantics = Some("ymir.WhittakerBiome@v1".to_string());

    let mut temperature = raster("temperature", "temperature.i16", "i16");
    temperature.unit = Some("celsius_x100".to_string());

    let mut precipitation = raster("precipitation", "precipitation.u16", "u16");
    precipitation.unit = Some("mm_per_year".to_string());

    vec![
        height,
        coastline,
        cliffs,
        flow_accumulation,
        rivers,
        lakes,
        lake_mask,
        biome,
        temperature,
        precipitation,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_meta(w: usize, h: usize) -> ContinentMeta {
        ContinentMeta {
            name: "Testonia".to_string(),
            seed: 12345,
            grid: Grid { width: w, height: h },
            window_km: 1024.0,
            tectonic_domain_km: 1024.0,
            window_offset_in_torus: [0.0, 0.0],
            stein_stein: SteinSteinParams::default(),
            sea_level_m: 0.0,
            max_elevation_m: 5650.0,
            max_depth_m: 5650.0,
        }
    }

    #[test]
    fn round_trip_synthetic_continent() {
        let dir = std::env::temp_dir().join("ymir_container_rt");
        // Fresh directory each run.
        let _ = std::fs::remove_dir_all(&dir);
        let (w, h) = (16usize, 16usize);

        let mut writer = ContinentWriter::new(&dir, tiny_meta(w, h)).unwrap();
        let height: Vec<u16> = (0..(w * h) as u16).collect();
        writer.add_raster_u16("height", &height).unwrap();
        let manifest_path = writer.finish().unwrap();

        // Re-read the manifest from disk.
        let json = std::fs::read_to_string(&manifest_path).unwrap();
        let m: Manifest = serde_json::from_str(&json).unwrap();

        // Grid dims + derived km_per_cell.
        assert_eq!(m.continent.grid.width, w);
        assert_eq!(m.continent.grid.height, h);
        assert!((m.continent.km_per_cell - 1024.0 / w as f64).abs() < 1e-9);
        assert!((m.continent.vertical_scale.sea_level_norm - 0.5).abs() < 1e-9);

        // Height file byte length == 16*16*2.
        let height_file = dir.join("height.u16");
        let bytes = std::fs::metadata(&height_file).unwrap().len();
        assert_eq!(bytes, (w * h * 2) as u64);

        // Present flags: only `height` is present.
        for layer in &m.layers {
            if layer.id == "height" {
                assert!(layer.present, "height must be present");
                assert_eq!(layer.dtype.as_deref(), Some("u16"));
                assert_eq!(layer.endianness.as_deref(), Some("le"));
                assert_eq!(layer.width, Some(w));
                assert_eq!(layer.height, Some(h));
            } else {
                assert!(!layer.present, "{} must be absent", layer.id);
            }
        }
    }

    #[test]
    fn wrong_raster_length_is_rejected() {
        let dir = std::env::temp_dir().join("ymir_container_badlen");
        let _ = std::fs::remove_dir_all(&dir);
        let mut writer = ContinentWriter::new(&dir, tiny_meta(4, 4)).unwrap();
        // 15 != 16 cells → error, never a silent truncation.
        assert!(writer.add_raster_u16("height", &vec![0u16; 15]).is_err());
    }
}
