//! C1 closures — empirical source terms applied per time step.
//!
//! Each closure is an additive contribution to `∂S̃/∂t`, evaluated
//! from the current state. Closures are activatable independently
//! via per-closure `enabled` flags so they can be ablated for
//! diagnostic purposes (W4 watchpoint of Issue #123, §6.3 of the
//! design doc).
//!
//! ## Phase 1.2 (Issue #123) status
//!
//! - [`davis_suppe`] — Davis-Suppe critical taper orogenic profile
//!   (this issue's deliverable). Active on upper-plate cells at
//!   convergent boundaries.
//!
//! ## Forthcoming closures (per design doc §5.1)
//!
//! - Equilibrium height (Molnar-Lyon-Caen) — Phase 1.3
//! - Oceanic bathymetry (Parsons-Sclater) — Phase 2
//! - Macro stream-power erosion (Whipple-Tucker) — Phase 1.4
//! - Airy isostasy — already available via
//!   `crate::tectonics::isostasy::compute_isostasy` (reused, not
//!   re-implemented)

pub mod davis_suppe;
