//! # Ymir Core
//!
//! Physically-grounded continent generation library.
//!
//! Ymir produces terrain through a multi-phase pipeline:
//! tectonic plate simulation → isostasy → anisotropic noise → hydraulic erosion →
//! climate modeling → biome classification.
//!
//! Each phase operates on a [`grid::GridF32`] heightmap and is controlled by a
//! deterministic [`seed::WorldSeed`] to ensure reproducibility.

pub mod config;
pub mod grid;
pub mod seed;

// Pipeline phases — module structure follows the TDD.
// Submodules are created as empty placeholders; each will be populated
// by its corresponding milestone issue.
pub mod climate;
pub mod erosion;
pub mod export;
pub mod lakes;
pub mod tectonics;
pub mod tectonics_v2;
pub mod terrain;
