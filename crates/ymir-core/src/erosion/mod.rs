//! Terrain erosion simulation phases.
//!
//! Each erosion type is an independent module that operates on a [`crate::grid::GridF32`]
//! heightmap. Hydraulic erosion (M0) is the foundation; other modes are added in
//! later milestones.

pub mod hydraulic;
// pub mod thermal;     // M5 — rockfall on steep slopes
// pub mod coastal;     // M5 — wave erosion at shorelines
// pub mod aeolian;     // M6 — wind erosion for arid continents
// pub mod glacial;     // M6 — ice carving for cold continents
