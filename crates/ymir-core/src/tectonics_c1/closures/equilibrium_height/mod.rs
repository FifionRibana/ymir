//! Equilibrium height closure — gravitational collapse cap on `S̃`
//! (C1 Phase 1.3, Issue #125).
//!
//! ## Physics in one paragraph
//!
//! Per Molnar & Lyon-Caen 1988, *J. Geophys. Res.* 93(B5), 4885-4923
//! and England & Houseman 1989, *J. Geophys. Res.* 94(B12), 17561-
//! 17579: a thickened continental lithosphere stores gravitational
//! potential energy proportional to the thickness excess above some
//! equilibrium value `h_eq`. Above that threshold, the gravitational
//! body force becomes large enough to drive ductile extension /
//! gravitational collapse, capping the elevated wedge at the
//! equilibrium height. The closure approximates this as a linear
//! one-sided relaxation:
//!
//! ```text
//!     ∂S̃/∂t = − k_collapse · max(0, S̃ − h_eq)
//! ```
//!
//! The `max(0, ·)` makes the closure **asymmetric**: cells above
//! `h_eq` collapse toward it, cells below are untouched (no thin-
//! lithosphere "lift"). The asymmetry is what makes this term a
//! pure sink — never a source.
//!
//! ## Differences from Davis-Suppe
//!
//! - **Applied globally.** No plate-type or boundary skip — every
//!   cell is subject to the same gravitational stability criterion,
//!   regardless of upper/lower-plate identity or boundary
//!   classification. This is the architectural difference that lets
//!   the closure cap Phase 1.2's boundary pile-up
//!   (`global_max ≈ 2297` in the advection-dominated 300-step run).
//! - **First global sink in C1.** Davis-Suppe is source-only;
//!   equilibrium height is the first closure that subtracts mass.
//!   Phase 1.4's stream-power erosion will be the second.
//! - **Bounded contribution.** For `k_collapse · dt ≤ 1` the per-
//!   step decrement is a convex combination of `S̃_old` and `h_eq`
//!   → no undershoot. The [`source_term`] module ships a defensive
//!   clamp that protects against pathological time steps; see
//!   [`source_term::apply_equilibrium_height_step`] for the in-
//!   function comment explaining when the clamp triggers.
//!
//! ## Interaction with Phase 1.2 Davis-Suppe imprint
//!
//! With the default `h_eq = 2.0` (below `DavisSuppeParams::h_max
//! = 2.5`), the Phase 1.2 statistics partition naturally:
//!
//! | Phase 1.2 quantity      | value  | relation to `h_eq = 2.0`     |
//! |-------------------------|--------|------------------------------|
//! | `mean(S̃)`              | 1.574  | below — untouched            |
//! | `wedge_p95`             | 0.376  | below — untouched            |
//! | `mean(d∈0..5)` (near)   | 0.904  | below — untouched            |
//! | `wedge_p99`             | 5.83   | above — collapsed to `h_eq`  |
//! | `global_max` (boundary) | 2297   | above — collapsed to `h_eq`  |
//!
//! The bulk wedge body (where the Davis-Suppe fill-ratio profile
//! sits) is preserved; only the sparse boundary pile-up tail is
//! capped. Stage E3 integration tests will verify both invariants
//! simultaneously.
//!
//! ## Module layout
//!
//! - [`params`] — [`params::EquilibriumHeightParams`] tunables
//!   (`h_eq = 2.0`, `k_collapse = 1.0`, `enabled = true`).
//! - [`source_term`] — [`source_term::apply_equilibrium_height_step`],
//!   the per-step in-place collapse plus 5 unit tests.
//!
//! ## References
//!
//! - Molnar, P. & Lyon-Caen, H. (1988). Some simple physical aspects
//!   of the support, structure, and evolution of mountain belts.
//!   *Geological Society of America Special Papers*, 218, 179-208.
//! - England, P. & Houseman, G. (1989). Extension during continental
//!   convergence, with application to the Tibetan Plateau.
//!   *J. Geophys. Res.*, 94(B12), 17561-17579.

pub mod params;
pub mod source_term;
