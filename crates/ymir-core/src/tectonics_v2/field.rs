//! Re-exports of the two grid utilities whitelisted for tectonics_v2.
//!
//! Audit (Step 0 entry condition #2): `Field2D` and `PeriodicIndex` in
//! the legacy `tectonics/solver/field.rs` have zero external imports,
//! are covered by stride- and wrap-aware unit tests (including coprime
//! and rectangular shapes), and carry no tectonic-specific state. They
//! are safe to re-export as plain grid primitives.
//!
//! Re-exporting (rather than copying) keeps a single source of truth
//! until the legacy `tectonics/` module is retired at the end of the
//! milestone. No other legacy symbol is imported into `tectonics_v2`.

pub use crate::tectonics::solver::field::{Field2D, PeriodicIndex};
