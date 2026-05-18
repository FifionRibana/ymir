//! Step 9 — cratonic immunity (split-mechanism design).
//!
//! Two mechanisms cooperate to differentiate stable plate interiors
//! from mobile belts:
//!
//! - **Primary plastic immunity.** In cratonic cells, the yield
//!   stress is held at the full `Bi` value regardless of accumulated
//!   plastic strain. Today the plastic branch is stateless (no
//!   plastic memory yet — that is post-milestone work), so this hook
//!   is mathematically a no-op in the current code path. The
//!   modulator is wired through `eta_plastic_with_cratonic` so that
//!   when plastic memory is added later, the immunity logic is
//!   already in the right place.
//! - **Secondary viscous contrast `K`.** The effective viscosity
//!   `η_eff[i]` is multiplied by `1 + (K - 1) · cratonic_factor[i]`,
//!   raising η in cratonic cells. K is bounded to [3, 8] (default 5)
//!   to keep `κ(A)` within the budget established in §4.10.
//!
//! The `cratonic_factor` field is computed once at init from the
//! Voronoï partition (see [`factor`] submodule) and stored as part
//! of the simulation state. Static identification matches §4.10's
//! "cells assigned at t = 0" semantics; dynamic craton evolution is
//! out of scope for this milestone.
//!
//! See `step9_issue.md` (D1–D9) for the full design contract.

pub mod factor;

/// Step 9 D1 primary mechanism (operational form) — yield-stress
/// elevation in cratonic cells.
///
/// **Original D1 formula** (literal reading of `step9_issue.md`):
///
/// ```text
///   yield_stress[i] = Bi · (cratonic_factor[i]
///                          + (1 - cratonic_factor[i]) · weakening(plastic_strain[i]))
/// ```
///
/// This formula was discovered to be a trivial no-op in the
/// stateless yielding regime of this milestone: `weakening` is
/// implicitly `1.0` everywhere (no plastic strain field), so the
/// expression collapses to `Bi · 1 = Bi` for any `cratonic_factor`.
/// Diagnostic on the Step 8-shape immunity test (32², mantle on)
/// showed `peak_yielding_in_craton = 0.99` despite the K = 5
/// secondary mechanism — cratons yielded essentially everywhere.
/// Acceptance #6 cannot be met with the original formula in this
/// milestone scope.
///
/// **§4.10 amendment (Step 9, this implementation).** The primary
/// mechanism generalises to "cratons have an *elevated* yield
/// strength" via a new parameter `B_factor ∈ [3, 10]`:
///
/// ```text
///   yield_stress[i] = Bi · (1 + (B_factor - 1) · cratonic_factor[i])
///                       · weakening(plastic_strain[i])
/// ```
///
/// Limits:
/// - `cratonic_factor = 0` (mobile): `yield_stress = Bi · weakening`
///   (= Bi today, plastic-memory-modulated later).
/// - `cratonic_factor = 1` (full cratonic core): `yield_stress =
///   B_factor · Bi · weakening` (= B_factor · Bi today).
///
/// In the current milestone (no plastic memory), `weakening = 1`
/// and this function returns `Bi · (1 + (B_factor - 1) · cratonic_factor)`.
/// When plastic memory is implemented, weakening will modulate
/// mobile belts; cratons' `plastic_strain` stays zero by D1, so
/// `weakening(0) = 1` and `B_factor · Bi` survives unmodified.
#[inline]
pub fn bi_with_cratonic_immunity(bi: f64, cratonic_factor: f64, b_factor: f64) -> f64 {
    // weakening(plastic_strain) = 1.0 (no plastic memory yet).
    let weakening: f64 = 1.0;
    bi * (1.0 + (b_factor - 1.0) * cratonic_factor) * weakening
}

use crate::tectonics_v2::field::Field2D;

/// Concrete parameters for `CratonicConfig::Enabled`.
///
/// All fields are nondimensional. The defaults match §4.10 baselines
/// and the Step 9 issue's D3–D6 design decisions.
#[derive(Clone, Copy, Debug)]
pub struct CratonicConfigEnabled {
    /// Target cratonic fraction within continental plates large
    /// enough to host a craton — `Cr ∈ [0.1, 0.5]`. A construction
    /// local to this design (not a literature number); see §4.10
    /// patch.
    pub cr: f64,
    /// Viscous contrast multiplier in cratonic cells —
    /// `η_eff *= 1 + (K - 1) · cratonic_factor`. `K ∈ [3, 8]`.
    /// `K = 1` is supported (effectively disables the viscous
    /// secondary while keeping the plastic-immunity hook live).
    pub k_viscous: f64,
    /// **Step 9 D1 primary mechanism (operational form).**
    /// Yield-stress (Bi) elevation factor in cratonic cells —
    /// `bi_eff[i] = Bi · (1 + (B_factor - 1) · cratonic_factor[i])`.
    /// `B_factor ∈ [3, 10]`, default `5.0`. `B_factor = 1` reduces
    /// the formula to the no-op identity (matches the
    /// pre-amendment behaviour where the primary mechanism was
    /// trivially `bi · 1 = bi` because `weakening = 1` without
    /// plastic memory).
    ///
    /// Why this exists: §4.10 D1's literal "yield stress maintained
    /// at full Bi" is trivial in stateless yielding (every cell
    /// already has yield = Bi). To make cratons actually resist
    /// viscoplastic yielding `η_p = Bi/(2·(ε̇+ε̇_min))` in active
    /// regimes, Bi itself must be elevated locally — this is the
    /// operational form of the primary mechanism in the absence of
    /// plastic memory. Once plastic memory lands, the formula
    /// `bi_eff · weakening(plastic_strain)` retains B_factor's
    /// effect in cratons (where plastic strain stays zero by D1)
    /// AND lets weakening modulate mobile belts.
    pub b_factor: f64,
    /// Lower bound on plate area (as a fraction of the **domain**) for
    /// a plate to receive a craton at simulation init. Plates below
    /// this threshold get `cratonic_factor = 0` everywhere (they
    /// represent fragments without geological time to consolidate a
    /// cratonic root). Range `[0.05, 0.20]`.
    ///
    /// **Semantics 1 — init-time exclusion.** This parameter is
    /// fraction-of-domain. Used by
    /// [`super::factor::build_cratonic_factor_field`] (Step 9) when the
    /// run starts. **It is *not* the per-cycle retention threshold**;
    /// see [`Self::craton_retention_threshold`] for the post-erosion
    /// recompute counterpart (Step 12).
    ///
    /// The two parameters are deliberately separate (Step 12 Phase 3
    /// finding): a Step 9 init check reads the plate's *cell count*
    /// fraction of the whole domain, while a Step 12 D4 recompute
    /// check reads each plate's *continental cell* fraction *within
    /// the plate*. Numerically equal at the default value (`0.10`)
    /// but conceptually distinct.
    pub plate_area_min: f64,
    /// Step 12 D4 — minimum continental cell fraction **within a
    /// plate** for that plate to retain its cratonic factor across a
    /// Phase A cycle. After erosion drops `continental_count[p] /
    /// total_cells_in_plate[p]` below this value, the plate's
    /// per-plate type flips from `Continental` to oceanic-equivalent
    /// in the recompute, BFS sources this plate's cells, and
    /// `cratonic_factor = 0` for every cell of that plate post-cycle.
    /// Range `[0.05, 0.95]`.
    ///
    /// **Semantics 2 — per-cycle retention.** This parameter is
    /// fraction-within-plate. Used **only** by Step 12's
    /// `tectonics_v2::workflow::phase_a` orchestrator
    /// (`recompute_cratonic_factor_for_cycle`); it has no effect on
    /// any pre-Step-12 code path, so legacy regression tests are
    /// bit-identical regardless of its value.
    ///
    /// **Default `0.10`.** Intentionally aligned with
    /// [`Self::plate_area_min`] default to preserve the
    /// pre-Step-12-Phase-3.5 behaviour bit-for-bit. Empirical
    /// calibration during Phase 8 reports may suggest a different
    /// default (e.g., `0.05` more permissive) once multi-cycle runs
    /// are characterised.
    pub craton_retention_threshold: f64,
    /// Smoothstep transition width, in units of `L_plate`. The
    /// transition runs from `d_mid - smoothing_width / 2` to
    /// `d_mid + smoothing_width / 2`. Default `0.05` (5 % of
    /// `L_plate`).
    pub smoothing_width: f64,
}

impl CratonicConfigEnabled {
    pub const CR_DEFAULT: f64 = 0.3;
    pub const K_VISCOUS_DEFAULT: f64 = 5.0;
    /// `B_factor` default — `8.0`, set to satisfy acceptance #6
    /// (`peak_yielding_in_craton ≤ 0.01`) on the Step 8 shape 32²
    /// immunity test. Derived from analytical threshold
    /// `B > η_v / (2·K·η_p_default) ≈ 6.1` in saturated regimes
    /// (`peak|v| ~ O(1)`, ε̇ large) and validated empirically by
    /// the `B_factor` sweep `tests/v2_step9_physics_and_sweep::
    /// step9_immunity_demo_b_factor_sweep_32sq`. B = 5 produces a
    /// narrow miss (`yc = 0.025`); B = 8 hits zero with margin
    /// (`yc = 0`); B = 10 is the plateau. See `solver-scaling-
    /// step9-patch.md` for the formal §4.10 amendment.
    pub const B_FACTOR_DEFAULT: f64 = 8.0;
    pub const PLATE_AREA_MIN_DEFAULT: f64 = 0.10;
    /// Step 12 D4 default — within-plate continental fraction below
    /// which a plate flips to oceanic-equivalent in the per-cycle
    /// craton recompute. Aligned with `PLATE_AREA_MIN_DEFAULT` for
    /// bit-identical regression behaviour pre-Step-12.
    pub const CRATON_RETENTION_THRESHOLD_DEFAULT: f64 = 0.10;
    pub const SMOOTHING_WIDTH_DEFAULT: f64 = 0.05;
}

impl Default for CratonicConfigEnabled {
    fn default() -> Self {
        Self {
            cr: Self::CR_DEFAULT,
            k_viscous: Self::K_VISCOUS_DEFAULT,
            b_factor: Self::B_FACTOR_DEFAULT,
            plate_area_min: Self::PLATE_AREA_MIN_DEFAULT,
            craton_retention_threshold: Self::CRATON_RETENTION_THRESHOLD_DEFAULT,
            smoothing_width: Self::SMOOTHING_WIDTH_DEFAULT,
        }
    }
}

/// Cratonic-immunity configuration.
///
/// `Disabled` short-circuits the entire pipeline: no factor field is
/// allocated, no eta multiplier is applied, no yield-stress modulator
/// runs. This is the path used by the Step 8 regression to verify
/// bit-identical output with pre-Step-9 code.
///
/// `Enabled(cfg)` activates both the plastic immunity hook and the
/// viscous contrast `K`. The cratonic factor field is computed once
/// from the Voronoï partition.
#[derive(Clone, Copy, Debug)]
pub enum CratonicConfig {
    Disabled,
    Enabled(CratonicConfigEnabled),
}

impl CratonicConfig {
    pub fn label(&self) -> &'static str {
        match self {
            CratonicConfig::Disabled => "disabled",
            CratonicConfig::Enabled(_) => "enabled",
        }
    }

    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "disabled" | "off" => Ok(CratonicConfig::Disabled),
            "enabled" | "on" => Ok(CratonicConfig::Enabled(Default::default())),
            other => Err(format!(
                "unknown --cratonic-config value '{}'; expected disabled|enabled",
                other,
            )),
        }
    }
}

impl Default for CratonicConfig {
    fn default() -> Self {
        CratonicConfig::Disabled
    }
}

/// Per-cell viscous multiplier `m[i] = 1 + (K - 1) · cratonic_factor[i]`.
///
/// Pre-computed once at init from the cratonic factor field and the
/// `K` value. Carried alongside the simulation state and applied to
/// `η_eff` whenever a viscosity field is built. Storing the
/// pre-multiplied form (instead of recomputing from factor + K every
/// call) keeps the inner solver hot path branch-free.
///
/// When `CratonicConfig::Disabled`, this struct is not constructed —
/// callers pass `None` to the eta-build path.
#[derive(Clone)]
pub struct CratonicState {
    /// Continuous factor in `[0, 1]`. `1` = full cratonic, `0` =
    /// fully mobile or non-continental.
    pub factor: Field2D,
    /// Cached secondary-mechanism viscous multiplier
    /// `eta_multiplier[i] = 1 + (K - 1) · factor[i]`. Applied
    /// post-blend to `η_eff` in the Stokes operator (see
    /// `rheology::build_eta_field`). Bounded `[1, K]`.
    pub eta_multiplier: Field2D,
    /// Cached primary-mechanism Bi multiplier
    /// `bi_multiplier[i] = 1 + (B_factor - 1) · factor[i]`. Applied
    /// pre-blend to `η_p = Bi/(2(ε̇+ε̇_min))` so the plastic branch
    /// in cratonic cells uses an elevated yield strength
    /// `B_factor · Bi`. Bounded `[1, B_factor]`. With B_factor = 1
    /// this field is identically 1 and the cratonic state reduces
    /// to "K viscous mult only" (the pre-amendment behaviour).
    pub bi_multiplier: Field2D,
    /// `K` used to build `eta_multiplier` — captured for diagnostics
    /// and reporting.
    pub k_viscous: f64,
    /// `B_factor` used to build `bi_multiplier` — captured for
    /// diagnostics and reporting.
    pub b_factor: f64,
}

impl CratonicState {
    /// Build the eta multiplier and Bi multiplier fields from the
    /// `cratonic_factor` field and the two cratonic mechanism
    /// parameters `K` (viscous, secondary) and `B_factor` (Bi
    /// elevation, primary). Pure function; no dependency on the
    /// BFS / smoothstep pipeline.
    pub fn from_factor(factor: Field2D, k_viscous: f64, b_factor: f64) -> Self {
        let nx = factor.nx();
        let ny = factor.ny();
        let mut eta_multiplier = Field2D::new(nx, ny);
        let mut bi_multiplier = Field2D::new(nx, ny);
        let k_minus_1 = k_viscous - 1.0;
        let b_minus_1 = b_factor - 1.0;
        for j in 0..ny {
            for i in 0..nx {
                let cf = factor.get(i, j);
                eta_multiplier.set(i, j, 1.0 + k_minus_1 * cf);
                bi_multiplier.set(i, j, 1.0 + b_minus_1 * cf);
            }
        }
        Self { factor, eta_multiplier, bi_multiplier, k_viscous, b_factor }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_default_is_disabled() {
        match CratonicConfig::default() {
            CratonicConfig::Disabled => {}
            _ => panic!("default should be Disabled"),
        }
    }

    #[test]
    fn enabled_defaults_match_issue_spec() {
        let cfg = CratonicConfigEnabled::default();
        assert_eq!(cfg.cr, 0.3);
        assert_eq!(cfg.k_viscous, 5.0);
        // B_factor default raised from 5 to 8 after the §4.10
        // amendment validation sweep — see B_FACTOR_DEFAULT
        // docstring for the analytical + empirical justification.
        assert_eq!(cfg.b_factor, 8.0);
        assert_eq!(cfg.plate_area_min, 0.10);
        // Step 12 Phase 3.5 — separate from `plate_area_min`, default
        // intentionally aligned to preserve pre-Step-12 bit-identity.
        assert_eq!(cfg.craton_retention_threshold, 0.10);
        assert_eq!(cfg.smoothing_width, 0.05);
    }

    #[test]
    fn parse_roundtrip() {
        match CratonicConfig::parse("disabled").unwrap() {
            CratonicConfig::Disabled => {}
            _ => panic!(),
        }
        match CratonicConfig::parse("enabled").unwrap() {
            CratonicConfig::Enabled(_) => {}
            _ => panic!(),
        }
        assert!(CratonicConfig::parse("garbage").is_err());
    }

    #[test]
    fn eta_multiplier_at_factor_zero_is_one() {
        let factor = Field2D::filled(4, 4, 0.0);
        let state = CratonicState::from_factor(factor, 5.0, 1.0);
        for v in state.eta_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }

    #[test]
    fn eta_multiplier_at_factor_one_is_k() {
        let factor = Field2D::filled(4, 4, 1.0);
        let state = CratonicState::from_factor(factor, 5.0, 1.0);
        for v in state.eta_multiplier.data() {
            assert_eq!(*v, 5.0);
        }
    }

    #[test]
    fn eta_multiplier_linear_in_factor() {
        let mut factor = Field2D::new(2, 2);
        factor.set(0, 0, 0.0);
        factor.set(1, 0, 0.5);
        factor.set(0, 1, 1.0);
        factor.set(1, 1, 0.25);
        let state = CratonicState::from_factor(factor, 5.0, 1.0);
        // 1 + (K-1)·f = 1 + 4·f
        assert_eq!(state.eta_multiplier.get(0, 0), 1.0);
        assert_eq!(state.eta_multiplier.get(1, 0), 3.0);
        assert_eq!(state.eta_multiplier.get(0, 1), 5.0);
        assert_eq!(state.eta_multiplier.get(1, 1), 2.0);
    }

    #[test]
    fn bi_immunity_b_factor_one_is_identity() {
        // B_factor = 1 ⟹ multiplier = 1 + 0·factor = 1, identity.
        // This recovers the original (pre-amendment) behaviour where
        // the primary mechanism was a no-op in stateless yielding.
        let bi = 0.15;
        for f in [0.0, 0.1, 0.3, 0.5, 0.8, 1.0] {
            let got = bi_with_cratonic_immunity(bi, f, 1.0);
            assert!(
                (got - bi).abs() < 1e-15,
                "bi_with_cratonic_immunity({}, {}, 1.0) = {} ≠ {}",
                bi, f, got, bi
            );
        }
    }

    #[test]
    fn bi_immunity_b_factor_default_elevates_in_cratons() {
        let bi = 0.15;
        let b = 5.0; // default
        // factor = 0 (mobile): bi_eff = bi (unchanged)
        let mobile = bi_with_cratonic_immunity(bi, 0.0, b);
        assert!((mobile - bi).abs() < 1e-15);
        // factor = 1 (full craton): bi_eff = B · bi
        let craton = bi_with_cratonic_immunity(bi, 1.0, b);
        assert!((craton - b * bi).abs() < 1e-15);
        // factor = 0.5 (boundary): bi_eff = bi · (1 + (B-1)·0.5) = 3 · bi
        let mid = bi_with_cratonic_immunity(bi, 0.5, b);
        assert!((mid - bi * (1.0 + (b - 1.0) * 0.5)).abs() < 1e-15);
    }

    #[test]
    fn cratonic_state_bi_multiplier_at_factor_zero_is_one() {
        let factor = Field2D::filled(4, 4, 0.0);
        let state = CratonicState::from_factor(factor, 5.0, 5.0);
        for v in state.bi_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }

    #[test]
    fn cratonic_state_bi_multiplier_at_factor_one_is_b_factor() {
        let factor = Field2D::filled(4, 4, 1.0);
        let state = CratonicState::from_factor(factor, 5.0, 7.0);
        for v in state.bi_multiplier.data() {
            assert_eq!(*v, 7.0);
        }
    }

    #[test]
    fn b_factor_one_makes_bi_multiplier_uniform_one() {
        // B_factor = 1 reduces cratonic state to "K viscous mult only"
        // — bi_multiplier identically 1, no Bi elevation anywhere.
        let mut factor = Field2D::new(2, 2);
        factor.set(0, 0, 0.0);
        factor.set(1, 0, 0.5);
        factor.set(0, 1, 1.0);
        factor.set(1, 1, 0.25);
        let state = CratonicState::from_factor(factor, 5.0, 1.0);
        for v in state.bi_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }

    #[test]
    fn k_one_makes_multiplier_uniform_one() {
        let mut factor = Field2D::new(2, 2);
        factor.set(0, 0, 0.0);
        factor.set(1, 0, 0.7);
        factor.set(0, 1, 1.0);
        factor.set(1, 1, 0.4);
        let state = CratonicState::from_factor(factor, 1.0, 1.0);
        for v in state.eta_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }
}
