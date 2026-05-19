//! Pipeline phase navigation — restored after Step 8.6 Phase 8h sunset.
//!
//! The post-sunset binary was v2-tectonics-only. This module re-introduces
//! the multi-phase workflow:
//!
//! - **Tectonics** — `tectonics_v2` solver (the v2 sprite + parameter
//!   panel + metrics dashboard from Step 8.6 augmenté).
//! - **Isostasy** — `compute_isostasy` on the v2 final S̃ → altitude
//!   `GridF32`.
//! - **UpscaleFbm** — bicubic interpolation + anisotropic FBM noise on
//!   the isostasy output → high-res `GridF32`.
//! - **Erosion** — particle hydraulic erosion on the upscaled
//!   heightmap → eroded `GridF32` + sediment `GridF32`.
//! - **Hydrology** — flow accumulation + pit-fill + river extraction
//!   + lake detection on the eroded heightmap.
//! - **Climate** — temperature + precipitation (stub in ymir-core).
//! - **Biome** — Whittaker classification (stub in ymir-core).
//!
//! This file ships only the enum + active-phase resource + thin left
//! toolbar. Each phase's panel / cache / dispatch / render systems
//! land in subsequent commits.

use bevy::prelude::*;

/// Pipeline phase — drives which view the central sprite shows and
/// which collapsible section the right parameter panel highlights.
/// `Tectonics` is the v2 view (current behaviour); the other phases
/// chain off the previous phase's output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PipelinePhase {
    Tectonics,
    Isostasy,
    UpscaleFbm,
    Erosion,
    Hydrology,
    Climate,
    Biome,
}

impl PipelinePhase {
    pub const ALL: &'static [PipelinePhase] = &[
        PipelinePhase::Tectonics,
        PipelinePhase::Isostasy,
        PipelinePhase::UpscaleFbm,
        PipelinePhase::Erosion,
        PipelinePhase::Hydrology,
        PipelinePhase::Climate,
        PipelinePhase::Biome,
    ];

    /// 3-letter tag for the left toolbar buttons.
    pub fn short_label(self) -> &'static str {
        match self {
            PipelinePhase::Tectonics => "TEC",
            PipelinePhase::Isostasy => "ISO",
            PipelinePhase::UpscaleFbm => "FBM",
            PipelinePhase::Erosion => "ERO",
            PipelinePhase::Hydrology => "HYD",
            PipelinePhase::Climate => "CLI",
            PipelinePhase::Biome => "BIO",
        }
    }

    /// Long label for hover tooltips + status displays.
    pub fn label(self) -> &'static str {
        match self {
            PipelinePhase::Tectonics => "Tectonics (v2 solver)",
            PipelinePhase::Isostasy => "Isostasy (Airy altitude)",
            PipelinePhase::UpscaleFbm => "Upscale + FBM noise",
            PipelinePhase::Erosion => "Erosion (hydraulic)",
            PipelinePhase::Hydrology => "Hydrology (flow + rivers + lakes)",
            PipelinePhase::Climate => "Climate (temperature + precipitation)",
            PipelinePhase::Biome => "Biome (Whittaker classification)",
        }
    }

    /// Whether this phase has its core logic implemented in
    /// `ymir-core`. `Climate` and `Biome` are still stubs (the
    /// pre-sunset binary also exposed them as placeholders).
    pub fn is_implemented(self) -> bool {
        !matches!(self, PipelinePhase::Climate | PipelinePhase::Biome)
    }
}

/// Resource holding the currently-active pipeline phase. Driven by
/// the left-side toolbar; read by every phase-aware render / panel
/// system.
#[derive(Resource, Debug, Clone, Copy)]
pub struct ActivePhase(pub PipelinePhase);

impl Default for ActivePhase {
    fn default() -> Self {
        ActivePhase(PipelinePhase::Tectonics)
    }
}
