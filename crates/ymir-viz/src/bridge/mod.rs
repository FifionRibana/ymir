//! Bridge module — v2 (Step 8.6) + c1 (Issue #137 Viz-0).
//!
//! The legacy `tectonics::` bridge (commands.rs, events.rs,
//! export_system.rs, plugin.rs, thread.rs) was removed once the
//! Phase 8g visual revalidation passed. The viz binary now wires
//! both `bridge::v2::*` (legacy Stokes solver) and
//! `bridge::c1::*` (lightweight dynamic tectonics) in parallel —
//! neither replaces the other; users select via the engine
//! switcher (Viz-0 Stage E5).

pub mod c1;
pub mod v2;
