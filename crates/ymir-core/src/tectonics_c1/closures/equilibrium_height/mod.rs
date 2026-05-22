//! Equilibrium height closure — gravitational collapse cap on `S̃`
//! (C1 Phase 1.3, Issue #125).
//!
//! ## Physics in one paragraph
//!
//! Per Molnar & Lyon-Caen 1988, *Geological Society of America
//! Special Papers* 218, 179-208 and England & Houseman 1989,
//! *J. Geophys. Res.* 94(B12), 17561-17579: a thickened
//! continental lithosphere stores gravitational potential energy
//! quadratically in the thickness excess above some equilibrium
//! value `h_eq`. Above that threshold, the gravitational body
//! force drives ductile extension / gravitational collapse,
//! capping the elevated wedge. The closure approximates this as
//! a quadratic one-sided relaxation:
//!
//! ```text
//!     ∂S̃/∂t = − k_collapse · max(0, S̃ − h_eq)²
//! ```
//!
//! The `max(0, ·)` makes the closure **asymmetric**: cells above
//! `h_eq` collapse toward it, cells below are untouched (no thin-
//! lithosphere "lift"). The asymmetry is what makes this term a
//! pure sink — never a source. The squared excess implements a
//! soft threshold (see § Formula derivation).
//!
//! ## Formula derivation
//!
//! Per Molnar & Lyon-Caen 1988, GSA Spec. Paper 218, eq. (2):
//!
//! ```text
//!     ΔPE_A = ½ρ_c·g·h² + ρ_c·g·h·H_0 + ½·Δρ·g·ΔH²
//! ```
//!
//! The gravitational potential energy excess is **quadratic** in
//! the excess thickness `ΔH`. Our sink term derives from this:
//!
//! ```text
//!     ∂S̃/∂t |_equilibrium = − k_collapse · max(0, S̃ − h_eq)²
//! ```
//!
//! Quadratic dependence implements threshold behavior:
//!
//! - **Small excess** → weak collapse (preserves the Davis-Suppe
//!   wedge body where excess above `h_eq` is < 1 or so).
//! - **Large excess** → strong collapse (caps boundary outliers,
//!   advection pile-up cells where excess can reach `≈ 2 × 10³`
//!   in the Phase 1.2 regime). One step is enough — the safety
//!   clamp in [`source_term`] catches the overshoot and holds the
//!   cell at `h_eq`.
//!
//! ## Parameter `h_eq = 2.0` (phenomenological)
//!
//! Corresponds to the observed Tibet crustal-thickness ratio
//! (~70 km plateau crust vs ~35 km normal continental crust →
//! ratio ≈ 2). This is **not** a rigorous "equilibrium height"
//! derivation; Molnar & Lyon-Caen 1988 itself defines a critical
//! thickness `S_c ≈ 50-60 km` (ratio 1.43-1.71) above which
//! convective removal of the mantle lithosphere is the active
//! mechanism — a more complex two-layer instability than this
//! single-field linear sink. We adopt `h_eq = 2.0` as a tunable
//! target consistent with observed plateau elevation limits, not
//! as a derived equilibrium value.
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
//! - **Bounded by clamp, not by linearity.** With the quadratic
//!   form, the per-step decrement `k_collapse · excess² · dt` can
//!   exceed `excess` for large `excess`; the safety clamp in
//!   [`source_term`] holds the cell at `h_eq` in that case. This
//!   is the intended threshold behavior — see
//!   [`source_term::apply_equilibrium_height_step`] for the
//!   in-function comment.
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
//! | `wedge_p99`             | 5.83   | above — clamped to `h_eq`    |
//! | `global_max` (boundary) | 2297   | above — clamped to `h_eq`    |
//!
//! The bulk wedge body (where the Davis-Suppe fill-ratio profile
//! sits) is preserved; only the sparse boundary pile-up tail is
//! capped. Stage E3 integration tests will verify both invariants
//! simultaneously.
//!
//! ## Module layout
//!
//! - [`params`] — [`params::EquilibriumHeightParams`] tunables
//!   (`h_eq = 2.0`, `k_collapse = 2.0`, `enabled = true`).
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
