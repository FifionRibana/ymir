//! Stein-Stein 1992 oceanic bathymetry closure — age-dependent
//! depth assignment on oceanic cells (C1 Phase 2 Track A, Issue
//! #129).
//!
//! ## Physics in one paragraph
//!
//! Per Stein, C. A. & Stein, S. (1992). "A model for the global
//! variation in oceanic depth and heat flow with lithospheric
//! age." *Nature* 359, 123-129 — oceanic lithosphere subsides as it
//! cools and thickens away from a mid-ocean ridge. The depth-age
//! relation has two regimes:
//!
//! ```text
//!     d(t) = d_r + b · √t                       for t < t_c
//!     d(t) = d_∞ - C · exp(-α · t)              for t ≥ t_c
//! ```
//!
//! - `d_r = 2600 m` — ridge-axis depth (paper Table 1).
//! - `b = 365 m / √Ma` — young-regime subsidence rate (√t cooling
//!   per half-space conductive model).
//! - `d_∞ = 5651 m` — asymptotic plate-model depth.
//! - `α = 0.0278 Ma⁻¹` — thermal time constant.
//! - `C = 2473 m` — continuity offset (`d_∞ - d_r + 22 m`
//!   accounting for the small jump at `t = t_c`); see
//!   [`source_term::stein_stein_depth`] for the rationale on
//!   hard-coding this constant rather than deriving it.
//! - `t_c = 20 Ma` — crossover age between regimes.
//!
//! ## Architecture C — post-isostasy bathymetry adjustment
//!
//! Unlike the additive source terms shipped Phase 1.2-1.4 (Davis-
//! Suppe, equilibrium height, stream-power erosion), Stein-Stein
//! is **not** an additive contribution to `∂S̃/∂t`. The name
//! `source_term.rs` is preserved for codebase consistency with the
//! other closures, but the operation is **depth assignment** on the
//! altitude field of oceanic cells based on their age.
//!
//! Why this design (Architecture C):
//!
//! - S-S directly publishes a depth value `d(t)`, not a rate
//!   `d/dt`. Inverting it into a source term would require
//!   approximating `∂d/∂t` and integrating back to `d` step by
//!   step — round-trip noise for no benefit.
//! - Altitude (not `S̃`) is what drives downstream consumption
//!   (drainage routing, visual reading). Modifying altitude
//!   directly puts the bathymetric signature exactly where the
//!   downstream pipeline reads it.
//! - The closure does **not** propagate back to `S̃`. `S̃` is left
//!   unchanged on oceanic cells. The next call to
//!   [`crate::tectonics::isostasy::compute_isostasy`] would
//!   recompute altitude from `S̃` and overwrite the S-S adjustment;
//!   this is intentional — S-S is reapplied each step inside the
//!   C1 time loop's stage 4a (see [`super::super::time_loop`]),
//!   immediately after isostasy and before drainage/erosion. Each
//!   step's altitude carries the S-S imprint; `S̃` itself does not.
//!
//! ### Fallback architectures
//!
//! If Stage D visual review finds Architecture C limitations (e.g.,
//! the S-S imprint is overwritten or smoothed out by downstream
//! stages too aggressively), two alternatives are documented for
//! future work:
//!
//! - **Architecture A — `S̃` source term.** Treat S-S as a
//!   prescribed *altitude target* and apply an equilibrium-height-
//!   style relaxation of `S̃` toward the value that, after Airy
//!   isostasy, yields the target altitude. Tighter coupling at the
//!   cost of an extra inversion per step.
//! - **Architecture B — hybrid.** Modify altitude (as C does) AND
//!   push a corresponding `S̃` correction into the closure state so
//!   subsequent isostasy recomputations preserve the S-S
//!   adjustment. Cleaner persistence at the cost of an `S̃`
//!   round-trip per step.
//!
//! These are not implemented in Phase 2 Track A. Stage D's visual
//! gallery is the empirical check on whether Architecture C
//! suffices.
//!
//! ## Parameter choices
//!
//! All defaults read directly from the Stein-Stein 1992 paper
//! (Table 1). The dimensional values (in meters and Ma) are
//! preserved exactly; the conversion to C1 non-dim units happens at
//! application time via two scale factors:
//!
//! - `age_to_ma = 0.667` — maps `1 age step ~ 0.667 Ma`, chosen so
//!   that the canonical `300 steps` Phase 1.x run spans `~200 Ma`,
//!   matching the upper end of typical oceanic-plate lifetimes
//!   from ridge formation to subduction.
//! - `depth_scale_m = 5000` — converts the S-S metric depth range
//!   `[2600, 5651] m` to non-dim altitude offsets `[0.52, 1.13]`,
//!   consistent with the Phase 1.4 isostatic altitude convention
//!   in `docs/c1_lightweight_dynamic_tectonics.md` §11.
//!
//! Both scales are documented in the design doc §11 sub-section
//! ("Phase 2 Track A scales").
//!
//! ## Known limitations
//!
//! 1. **Uniform plate-model parameters.** Stein-Stein's `d_∞ = 5651`
//!    and `α = 0.0278` are global averages. Pacific vs Atlantic
//!    regional variation is documented in the original paper but
//!    not modulated here; C1 uses one global tuning.
//! 2. **No sediment loading.** Real bathymetry is sediment-loaded
//!    in old basins (3-4 km of pelagic sediment on >100 Ma seafloor
//!    raises observed depths by ~1 km). Not modeled; the closure
//!    publishes basement depth, not seafloor depth.
//! 3. **Continental cells unchanged.** Plate-type discrimination is
//!    via `PlateTypeField` from the upstream tessellation. Mixed
//!    cells (e.g., recently-subducted oceanic now under continental
//!    accretion) are treated as their current `PlateType` says, with
//!    no transient handling.
//! 4. **No ridge-position constraint.** S-S assumes age=0 cells sit
//!    on a mid-ocean ridge. Phase 2 Track A does not yet enforce
//!    that the C1 age field initialises to 0 only along ridge-like
//!    geometries; oceanic cells with `age = 0` anywhere get ridge-
//!    depth bathymetry. Phase 2 Track B (R7 init, Issue TBD) will
//!    address the age field initialisation pattern.
//!
//! ## Interaction with Phase 1.4 erosion
//!
//! Architecture C places S-S **inside the C1 time-loop's stage 4a**
//! — immediately after `compute_isostasy` and before drainage
//! routing / erosion (stages 4c-4e). Drainage classification sees
//! the S-S-modulated altitude, so flow accumulation on oceanic
//! cells reflects age-dependent bathymetric gradients (older
//! oceanic plateaus drain toward ridges; younger ridge cells are
//! local highs on the bathymetry). Stream-power erosion then runs
//! on the S-S-modulated altitude and the resulting drainage areas,
//! so the erosion footprint near coastlines emerges from the joint
//! Davis-Suppe (wedge uplift) + isostasy + S-S (oceanic
//! bathymetry) + drainage + erosion stack.
//!
//! See [`super::super::time_loop::run_with_closures`] for the
//! per-step pipeline ordering and the rationale for stage 4a
//! grouping (isostasy + S-S as the "altitude preparation" pair).
//!
//! ## Module layout
//!
//! - [`params`] — [`params::SteinSteinParams`] tunables
//!   (`ridge_depth_m = 2600`, `subsidence_rate = 365`,
//!   `asymptotic_depth_m = 5651`, `time_constant = 0.0278`,
//!   `crossover_age_ma = 20`, `depth_scale_m = 5000`,
//!   `age_to_ma = 0.667`, `enabled = true`).
//! - [`source_term`] — `stein_stein_depth` formula + the per-step
//!   in-place [`source_term::apply_stein_stein_bathymetry`], plus
//!   the 6 unit tests covering the 4 formula regimes and the 2
//!   apply-side gates (`enabled = false`, continental skip).
//!
//! ## References
//!
//! - Stein, C. A. & Stein, S. (1992). A model for the global
//!   variation in oceanic depth and heat flow with lithospheric
//!   age. *Nature* 359, 123-129. doi:10.1038/359123a0
//! - Parsons, B. & Sclater, J. G. (1977). An analysis of the
//!   variation of ocean floor bathymetry and heat flow with age.
//!   *J. Geophys. Res.* 82(5), 803-827. doi:10.1029/JB082i005p00803
//!   (predecessor; superseded by Stein-Stein 1992 for the old-age
//!   regime where the simple half-space √t law overpredicts
//!   subsidence)

pub mod params;
pub mod source_term;

pub use params::SteinSteinParams;
pub use source_term::apply_stein_stein_bathymetry;
