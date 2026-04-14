//! Tectonic plate simulation and crustal deformation.
//!
//! Implements the thin viscous sheet model (England & McKenzie 1982) to simulate
//! continental formation from plate dynamics. The output is a crustal thickness
//! field that is converted to altitude via isostasy.

pub mod boundaries; // M1 — subduction, rifting, volcanism source terms
pub mod plates; // M1 — plate initialization, Voronoi partitioning
pub mod solver; // M1 — thin viscous sheet velocity solver
// pub mod advection;   // M1 — crustal thickness advection
pub mod isostasy; // M1 — thickness → altitude conversion
