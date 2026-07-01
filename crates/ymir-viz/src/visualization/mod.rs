//! Visualization module — single-engine (C1) after the v2 sunset.
//!
//! Renders the C1 raster sprite (`c1_plugin` + `c1_viz`) plus the shared
//! colormap helpers and overlay routines (Voronoï boundaries, velocity arrows).

// The coarse-gallery presentation (`c1_plugin` + `c1_viz` + `overlay`) is
// retained but no longer registered — the HD workspace (step d1) supersedes
// it. `colormap` is shared (the workspace reuses the data palettes).
pub mod c1_plugin;
pub mod c1_viz;
pub mod colormap;
pub mod overlay;
