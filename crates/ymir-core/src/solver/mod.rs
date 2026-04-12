//! Thin viscous sheet solver for tectonic simulation.
//!
//! Implements the England & McKenzie (1982) model: solves for crustal velocity
//! via Stokes equations with Picard linearization and preconditioned CG,
//! then advects crustal thickness with an upwind scheme.

pub mod advection;
pub mod cg;
pub mod config;
pub mod field;
pub mod grid;
pub mod picard;
pub mod plates;
pub mod stokes;
pub mod tectonics;
pub mod workspace;
