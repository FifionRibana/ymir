//! C1 closures — empirical source terms applied per time step.
//!
//! Most closures are additive contributions to `∂S̃/∂t` evaluated
//! from the current state. The Phase 2 Track A oceanic-bathymetry
//! closure is an exception (Architecture C — post-isostasy altitude
//! assignment; see [`oceanic_bathymetry`]). All closures are
//! activatable independently via per-closure `enabled` flags so
//! they can be ablated for diagnostic purposes (W4 watchpoint of
//! Issue #123, §6.3 of the design doc).
//!
//! ## Shipped closures
//!
//! - [`davis_suppe`] — Davis-Suppe critical taper orogenic profile
//!   (Phase 1.2, Issue #123). Active on upper-plate cells at
//!   convergent boundaries.
//! - [`equilibrium_height`] — Molnar-Lyon-Caen gravitational
//!   collapse sink (Phase 1.3, Issue #125). Caps `S̃` excess above
//!   `h_eq` via a quadratic relaxation.
//! - [`erosion`] — Whipple-Tucker 1999 stream-power incision sink
//!   (Phase 1.4, Issue #127). Erodes `S̃` based on drainage area
//!   and altitude slope; Lague 2014 calibration discipline applied
//!   to `K`.
//! - [`oceanic_bathymetry`] — Stein-Stein 1992 age-dependent
//!   oceanic depth (Phase 2 Track A, Issue #129). **Architecture C:
//!   modifies altitude on oceanic cells, not `S̃`** — see its
//!   module docstring for the rationale and fallback architectures.
//!
//! ## Forthcoming closures (per design doc §5.1, §5.2)
//!
//! - Subduction arc (Lallemand) — Phase 2+
//! - Rifting / passive margins (McKenzie-Buck) — Phase 2+
//! - Foreland basin flexure (Beaumont) — Phase 2+
//! - Airy isostasy — already available via
//!   `crate::tectonics::isostasy::compute_isostasy` (reused, not
//!   re-implemented)

pub mod davis_suppe;
pub mod equilibrium_height;
pub mod erosion;
pub mod oceanic_bathymetry;
