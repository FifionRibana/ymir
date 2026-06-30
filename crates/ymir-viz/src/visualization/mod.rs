//! Visualization module — single-engine (C1) after the v2 sunset.
//!
//! Renders the C1 raster sprite (`c1_plugin` + `c1_viz`) plus the shared
//! colormap helpers and overlay routines (Voronoï boundaries, velocity arrows).

pub mod c1_plugin;
pub mod c1_viz;
pub mod colormap;
pub mod overlay;

pub use c1_plugin::C1VisualizationPlugin;
