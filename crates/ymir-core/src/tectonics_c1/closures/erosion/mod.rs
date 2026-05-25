//! Stream-power erosion closure — surface incision sink on `S̃`
//! (C1 Phase 1.4, Issue #127).
//!
//! ## Physics in one paragraph
//!
//! Per Whipple & Tucker 1999, *J. Geophys. Res.* 104(B8),
//! 17661-17674, eq. (1) — the canonical Stream-Power Incision
//! Model (SPIM):
//!
//! ```text
//!     ∂h/∂t |_erosion = − K · A^m · S^n
//! ```
//!
//! where:

//! - `A` = drainage area (transitive count of upstream cells
//!   contributing flow through this cell, in cells; computed
//!   from `DrainageMap::target_idx` via path-length-descending
//!   accumulation — see the private `compute_drainage_areas`
//!   helper in [`source_term`]).
//! - `S` = local slope magnitude on the altitude heightmap
//!   (centered-difference gradient with periodic wraparound).
//! - `K` = erosion coefficient — **calibrated for visual
//!   balance, not derived from literature** (see § Parameter
//!   choices below).
//! - `m`, `n` = positive exponents.
//!
//! ## Formula derivation
//!
//! The canonical SPIM is the simplest defensible model of
//! fluvial incision into bedrock. The derivation (W-T § 2):
//! basal shear stress `τ_b ∝ ρ g D S` (where `D` = flow depth
//! ∝ √A under uniform precipitation), incision rate `∝ τ_b - τ_c`
//! (excess shear stress above some threshold), which with
//! `τ_c = 0` and standard hydraulic geometry simplifies to
//! `E ∝ A^m · S^n`.
//!
//! ## Parameter choices
//!
//! Default `(m, n) = (0.5, 1.0)`:
//! - `m / n = 0.5` satisfies the Whipple-Tucker 1999 constraint
//!   (W-T eq. 4 + § "Implications for landscape morphology",
//!   p. 17,665).
//! - `n = 1` is **canonical** in W-T: linear in slope, numerically
//!   stable, time-response to uplift perturbations is
//!   uplift-magnitude-independent (W-T § "Response timescale to
//!   step changes in uplift rate").
//! - Lague 2014 *ESPL* 39(1), 38-61, argues empirically that
//!   `n ≈ 2` better fits field data — leaves a planned upgrade
//!   point analogous to Phase 1.3's `Stage E1.bis` linear →
//!   quadratic refinement on the equilibrium-height closure.
//!
//! Default `K = 0.001`:
//! - Lague 2014 § 5 *explicitly declines* to publish a universal
//!   `K` value, arguing it aggregates lithology + climate +
//!   threshold effects + grain size and is not transferable
//!   between geological contexts. The published-K table is
//!   missing by design.
//! - C1's `K` is calibrated against the joint behaviour of
//!   Davis-Suppe (Phase 1.2) + equilibrium height (Phase 1.3) so
//!   the visual erosion signature is legible without erasing the
//!   wedges. See `docs/c1_lightweight_dynamic_tectonics.md` §11.1
//!   ("Calibration via visual review, not dimensional
//!   derivation") for the discipline applied to this kind of
//!   tunable.
//! - First-pass analytical estimate: `K ≈ Ũ / (Ã_eff^m · S̃^n)
//!   ≈ 0.01 / (50^0.5 · 1) ≈ 0.0014`, rounded down to `0.001`
//!   as the Stage E3 starting value. Stage E3 visual review
//!   may adjust within a 3-iteration budget (see §11.1 of the
//!   design doc).
//!
//! Default `floor = 0.2`:
//! - Matches the oceanic-initialisation S̃ baseline produced by
//!   `init_s_field` (continental ≈ 1.0, oceanic ≈ 0.2 — see
//!   `docs/c1_lightweight_dynamic_tectonics.md` §11 table).
//! - **Defensive clamp.** Continental cells eroded below this
//!   become oceanic via the subsequent reclassification step in
//!   `apply_post_tectonic` (`s > sea_level_ref` → continental,
//!   else oceanic — Phase 3.5 sea-level formula in S̃ units).
//!   The floor prevents pathological `K` calibrations from
//!   driving `S̃` arbitrarily negative.
//!
//! ## Known limitations
//!
//! 1. **No threshold (τ_c = 0).** Per Whipple-Tucker canonical
//!    form. Lague 2014 critiques this as the source of an
//!    at-equilibrium bias (slope responds incorrectly to climate
//!    perturbations without a threshold). Deferred to a Phase 5+
//!    optional enrichment closure.
//! 2. **Uniform `K` everywhere.** Continental and oceanic cells
//!    use the same coefficient. Phase 2+ may refine via
//!    lithology-aware modulation (cratonic immunity already
//!    available via `C1State::cratonic_mask`, but consumed at
//!    the closure layer is a Phase 2+ extension).
//! 3. **Block-uplift assumption violated.** W-T eq. (22) derives
//!    landscape-equilibrium predictions under spatially uniform
//!    block uplift; C1's Davis-Suppe source applies *spatial*
//!    uplift (focused on wedge bodies). The mismatch is tolerated
//!    for visual plausibility per design doc §8.5 ("Departure
//!    from physical fidelity").
//!
//! ## Erosion mutates `S̃` directly — implicit isostatic compensation
//!
//! W-T's `∂h/∂t` formula is on altitude, but C1's state is the
//! crustal-thickness field `S̃`. Erosion in this closure reduces
//! `S̃` directly. Mathematically: removing `Δh` from `S̃` produces
//! altitude reduction `Δh · (1 − ρ_crust / ρ_mantle) ≈ 0.17 · Δh`
//! via Airy isostasy (see `IsostasyConfig::default` for the
//! density ratio). The 0.83 of mass "removed" by surface erosion
//! isostatically rebounds into the same column — implicit,
//! instantaneous compensation.
//!
//! This is the standard C1 simplification: per design doc §8.5
//! ("visual plausibility, not physical fidelity") and §11.1
//! ("calibration via visual review"), the closure is calibrated
//! against the resulting altitude effect, not against W-T's
//! published rates. A future Phase 2+ refinement could decouple
//! erosion-on-altitude from isostatic-rebound-on-S̃ via a separate
//! pass; not in scope here.
//!
//! ## Interaction with Phase 1.3 equilibrium height
//!
//! Both Phase 1.3's equilibrium-height closure and Phase 1.4's
//! stream-power erosion are sinks, with **different physics**:
//!
//! - **Equilibrium height** — gravitational collapse on `S̃`
//!   excess above `h_eq` (vertical mass redistribution).
//!   Asymmetric one-sided sink; clamps at `h_eq`.
//! - **Stream-power erosion** — surface flux process on
//!   altitude + drainage area; scales with slope and upstream
//!   accumulation. No fixed clamp; rate-driven.
//!
//! Per Whipple-Tucker 1999 (citing Molnar-Lyon-Caen 1988): *"the
//! height of mountain ranges is limited by either crustal
//! strength or by a balance between rock uplift and erosion,
//! whichever is more restrictive."* In C1 terms:
//!
//! ```text
//!     h_effective ≈ min(h_collapse, h_erosion)
//! ```
//!
//! The two closures operate independently per step (equilibrium
//! before erosion in the time-loop order — see `time_loop.rs`).
//! Whichever is the more restrictive limit emerges naturally
//! from the joint dynamics; no explicit coupling is required.
//! Stage E3's `K` calibration is what tunes which of the two
//! dominates.
//!
//! ## Module layout
//!
//! - [`params`] — [`params::ErosionParams`] tunables
//!   (`k = 0.001`, `m = 0.5`, `n = 1.0`, `floor = 0.2`,
//!   `enabled = true`).
//! - [`source_term`] — [`source_term::apply_erosion_step`], the
//!   per-step in-place erosion plus the private
//!   `compute_drainage_areas` (transitive accumulation from
//!   `DrainageMap`) and `compute_local_slope` (periodic
//!   centered-difference gradient on `GridF32`). Plus the 7
//!   unit tests.
//!
//! ## References
//!
//! - Whipple, K. X. & Tucker, G. E. (1999). Dynamics of the
//!   stream-power river incision model: Implications for height
//!   limits of mountain ranges, landscape response timescales,
//!   and research needs.
//!   *J. Geophys. Res.* 104(B8), 17661-17674.
//!   doi:10.1029/1999JB900120
//! - Lague, D. (2014). The stream power river incision model:
//!   evidence, theory and beyond.
//!   *Earth Surf. Processes Landforms* 39(1), 38-61.
//!   doi:10.1002/esp.3462
//! - Molnar, P. & Lyon-Caen, H. (1988). Some simple physical
//!   aspects of the support, structure, and evolution of
//!   mountain belts (cited via W-T 1999 for the
//!   `min(h_collapse, h_erosion)` framing).

pub mod params;
pub mod source_term;

pub use params::ErosionParams;
pub use source_term::apply_erosion_step;
