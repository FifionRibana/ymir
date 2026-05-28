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
//! - [`subduction`] — Lallemand 2005 oceanic-mass consumption +
//!   arc volcanism + floor-triggered `plate_id` reassignment
//!   (Phase 2 Track D, Issue #132). **First C1 closure to mutate
//!   `plate_id` and `plate_type`** — breaks the static-classification
//!   optimisation in `tectonics_c1::time_loop::run_with_closures`
//!   (per-step `classify_boundaries` recompute lands at Stage E4
//!   of Track D).
//!
//! ## Forthcoming closures (per design doc §5.1, §5.2)
//!
//! - Accretion mechanism (sustained-convergence plate_id merge) —
//!   Phase 2 Track D Stage E2 (Issue #132)
//! - Rifting closure + split mechanism (McKenzie-Buck thinning +
//!   "chewing-gum cut" two-condition split) — Phase 2 Track D
//!   Stage E3 (Issue #132)
//! - Foreland basin flexure (Beaumont) — Phase 3
//! - Airy isostasy — already available via
//!   `crate::tectonics::isostasy::compute_isostasy` (reused, not
//!   re-implemented)

pub mod davis_suppe;
pub mod equilibrium_height;
pub mod erosion;
pub mod oceanic_bathymetry;
pub mod subduction;
