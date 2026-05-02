//! Bridge module — v2-only after Step 8.6 Phase 8h sunset.
//!
//! The legacy `tectonics::` bridge (commands.rs, events.rs,
//! export_system.rs, plugin.rs, thread.rs) was removed once the
//! Phase 8g visual revalidation passed. The viz binary is wired to
//! `bridge::v2::*` exclusively.

pub mod v2;
