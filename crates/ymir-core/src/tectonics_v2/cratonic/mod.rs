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

/// Step 9 D1 primary mechanism — yield-stress modulator hook for
/// plastic memory immunity in cratonic cells.
///
/// Formula (from `step9_issue.md` D1):
///
/// ```text
///   yield_stress[i] = Bi · (cratonic_factor[i]
///                          + (1 - cratonic_factor[i]) · weakening(plastic_strain[i]))
/// ```
///
/// In a cratonic cell (`cratonic_factor → 1`), `yield_stress → Bi`
/// regardless of accumulated plastic strain — the craton is immune
/// to plastic weakening. In a mobile cell (`cratonic_factor → 0`),
/// the formula reduces to `yield_stress = Bi · weakening(plastic_strain)`.
///
/// **Current behaviour — NO-OP.** Plastic memory has not yet landed
/// in the milestone; the `weakening` function is implicitly `1.0`
/// everywhere (no plastic strain accumulates in the stateless Step 0–8
/// yielding model). With `weakening = 1`, the formula collapses to
/// `yield_stress = Bi · 1 = Bi` for any value of `cratonic_factor`,
/// so the function returns the input `bi` unchanged.
///
/// The hook is wired into the codebase now (rather than at the time
/// plastic memory is added) so that the structural integration —
/// where the `cratonic_factor` field is *available* at the yield-
/// stress evaluation site — is already in place. When plastic memory
/// arrives, the only change required will be replacing the literal
/// `1.0` with `weakening(plastic_strain[i])`.
///
/// The active mechanism in current Step 9 baseline is the secondary
/// viscous contrast `K` (see [`CratonicState::eta_multiplier`]); the
/// primary plastic-immunity hook becomes observably effective only
/// once plastic memory is added.
#[inline]
pub fn bi_with_cratonic_immunity(bi: f64, cratonic_factor: f64) -> f64 {
    // weakening(plastic_strain) = 1.0 (no plastic memory yet).
    // → bi_eff = bi · (cratonic_factor + (1 - cratonic_factor) · 1) = bi
    let weakening: f64 = 1.0;
    bi * (cratonic_factor + (1.0 - cratonic_factor) * weakening)
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
    /// Lower bound on plate area (as a fraction of the domain) for a
    /// plate to receive a craton. Plates below this threshold get
    /// `cratonic_factor = 0` everywhere (they represent fragments
    /// without geological time to consolidate a cratonic root).
    /// Range `[0.05, 0.20]`.
    pub plate_area_min: f64,
    /// Smoothstep transition width, in units of `L_plate`. The
    /// transition runs from `d_mid - smoothing_width / 2` to
    /// `d_mid + smoothing_width / 2`. Default `0.05` (5 % of
    /// `L_plate`).
    pub smoothing_width: f64,
}

impl CratonicConfigEnabled {
    pub const CR_DEFAULT: f64 = 0.3;
    pub const K_VISCOUS_DEFAULT: f64 = 5.0;
    pub const PLATE_AREA_MIN_DEFAULT: f64 = 0.10;
    pub const SMOOTHING_WIDTH_DEFAULT: f64 = 0.05;
}

impl Default for CratonicConfigEnabled {
    fn default() -> Self {
        Self {
            cr: Self::CR_DEFAULT,
            k_viscous: Self::K_VISCOUS_DEFAULT,
            plate_area_min: Self::PLATE_AREA_MIN_DEFAULT,
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
    /// Cached `1 + (K - 1) · factor`. Mirrors `factor` cell-for-cell,
    /// kept ready so `build_eta_field` does not branch on a
    /// per-cell K.
    pub eta_multiplier: Field2D,
    /// `K` used to build `eta_multiplier` — captured for diagnostics
    /// and reporting.
    pub k_viscous: f64,
}

impl CratonicState {
    /// Build the viscous multiplier field
    /// `m[i] = 1 + (K - 1) · cratonic_factor[i]` from the factor
    /// field and the configured `K`. Pure function; no dependency
    /// on the BFS / smoothstep pipeline.
    pub fn from_factor(factor: Field2D, k_viscous: f64) -> Self {
        let nx = factor.nx();
        let ny = factor.ny();
        let mut eta_multiplier = Field2D::new(nx, ny);
        let k_minus_1 = k_viscous - 1.0;
        for j in 0..ny {
            for i in 0..nx {
                let cf = factor.get(i, j);
                eta_multiplier.set(i, j, 1.0 + k_minus_1 * cf);
            }
        }
        Self { factor, eta_multiplier, k_viscous }
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
        assert_eq!(cfg.plate_area_min, 0.10);
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
        let state = CratonicState::from_factor(factor, 5.0);
        for v in state.eta_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }

    #[test]
    fn eta_multiplier_at_factor_one_is_k() {
        let factor = Field2D::filled(4, 4, 1.0);
        let state = CratonicState::from_factor(factor, 5.0);
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
        let state = CratonicState::from_factor(factor, 5.0);
        // 1 + (K-1)·f = 1 + 4·f
        assert_eq!(state.eta_multiplier.get(0, 0), 1.0);
        assert_eq!(state.eta_multiplier.get(1, 0), 3.0);
        assert_eq!(state.eta_multiplier.get(0, 1), 5.0);
        assert_eq!(state.eta_multiplier.get(1, 1), 2.0);
    }

    #[test]
    fn bi_immunity_hook_is_no_op_in_current_codepath() {
        // With weakening = 1 (no plastic memory yet), the formula
        // reduces to `bi · 1 = bi` for any cratonic_factor.
        let bi = 0.15;
        for f in [0.0, 0.1, 0.3, 0.5, 0.8, 1.0] {
            let got = bi_with_cratonic_immunity(bi, f);
            assert!(
                (got - bi).abs() < 1e-15,
                "bi_with_cratonic_immunity({}, {}) = {} ≠ {}",
                bi, f, got, bi
            );
        }
    }

    #[test]
    fn k_one_makes_multiplier_uniform_one() {
        let mut factor = Field2D::new(2, 2);
        factor.set(0, 0, 0.0);
        factor.set(1, 0, 0.7);
        factor.set(0, 1, 1.0);
        factor.set(1, 1, 0.4);
        let state = CratonicState::from_factor(factor, 1.0);
        for v in state.eta_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }
}
