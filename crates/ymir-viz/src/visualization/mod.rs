//! Visualization module — v2-only after Step 8.6 Phase 8h sunset.
//!
//! Pre-sunset modules (`erosion`, `isostasy`, `plugin`, `render`,
//! `rivers`, `upscale`) drove the legacy pipeline phases. Post-sunset
//! the binary renders only the v2 raster sprite (`v2_viz`) plus the
//! shared colormap helpers and the Phase 8b overlay routines.

pub mod c1_plugin;
pub mod c1_viz;
pub mod colormap;
pub mod overlay;
pub mod v2_viz;

pub use c1_plugin::{ActiveEngine, C1VisualizationPlugin, C1VizSprite, C1VizState};
pub use c1_viz::{field_to_rgba, C1Field};
pub use v2_viz::V2VisualizationPlugin;
