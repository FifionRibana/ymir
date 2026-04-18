//! Thin viscous sheet solver for tectonic simulation.
//!
//! Implements the England & McKenzie (1982) model: solves for crustal velocity
//! via Stokes equations with Picard or Newton-Krylov linearization,
//! then advects crustal thickness with an upwind scheme.

pub mod advection;
pub mod config;
pub mod field;
pub mod grid;
pub mod linear_solve;
pub mod newton;
pub mod picard;
pub mod smooth;
pub mod stokes;
pub mod substep;
pub mod tectonics;
pub mod traction;
pub mod workspace;
