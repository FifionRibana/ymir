//! Tunables for the Stein-Stein 1992 oceanic bathymetry closure.
//!
//! See [`super`] for the physics derivation, the Architecture C
//! rationale (post-isostasy altitude assignment rather than `S̃`
//! source term), and the interaction with Phase 1.4 stream-power
//! erosion through the joint stage-4a altitude preparation.

/// Stein-Stein 1992 oceanic bathymetry closure tunables.
///
/// Default values read directly from Stein, C. A. & Stein, S.
/// (1992), *Nature* 359, Table 1 — the canonical "GDH1" plate-
/// model parameters. The two scale factors (`depth_scale_m` and
/// `age_to_ma`) are C1-specific conversions between the paper's
/// SI units and C1's non-dim units; see [`super`] § "Parameter
/// choices" and `docs/c1_lightweight_dynamic_tectonics.md` §11
/// sub-section on Phase 2 Track A scales.
///
/// `enabled` follows the same W4 watchpoint convention as
/// [`crate::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams`]:
/// when `false`, [`super::source_term::apply_stein_stein_bathymetry`]
/// is a no-op and the run reproduces the caller's previous
/// behaviour bit-identically.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct SteinSteinParams {
    /// Master enable/disable. When `false`,
    /// [`super::source_term::apply_stein_stein_bathymetry`] is a
    /// no-op (W4 closure-isolation discipline).
    pub enabled: bool,
    /// Ridge-axis depth `d_r` in meters. S-S 1992 Table 1: `2600 m`.
    pub ridge_depth_m: f64,
    /// Young-regime subsidence rate `b` in `m / √Ma`. S-S 1992
    /// Table 1: `365 m / √Ma`. Square-root cooling per half-space
    /// conductive model; applies for `age_ma < crossover_age_ma`.
    pub subsidence_rate: f64,
    /// Asymptotic plate-model depth `d_∞` in meters. S-S 1992
    /// Table 1: `5651 m`.
    pub asymptotic_depth_m: f64,
    /// Thermal time constant `α` in `Ma⁻¹`. S-S 1992 Table 1:
    /// `0.0278 Ma⁻¹`. Used in the old-regime exponential
    /// `d = d_∞ - C · exp(-α·t)`.
    pub time_constant: f64,
    /// Age threshold `t_c` in Ma separating the young (`√t`) and
    /// old (`exp(-α·t)`) regimes. S-S 1992 § "Plate model
    /// vs half-space cooling": `20 Ma`.
    pub crossover_age_ma: f64,
    /// Depth normalisation `[m / non-dim altitude unit]`. Converts
    /// the S-S metric depth range `[2600, 5651] m` to the C1
    /// non-dim altitude offset range `[0.52, 1.13]`, consistent
    /// with the Phase 1.4 isostatic altitude convention. Documented
    /// in design doc §11 sub-section on Phase 2 Track A scales.
    pub depth_scale_m: f64,
    /// Age normalisation `[Ma / non-dim age step]`. Maps `1 age
    /// step ~ 0.667 Ma` so the canonical `300 steps` Phase 1.x
    /// run spans `~200 Ma`, matching the upper end of typical
    /// oceanic-plate lifetimes from ridge formation to subduction.
    /// Documented in design doc §11 sub-section on Phase 2 Track A
    /// scales.
    pub age_to_ma: f64,
}

impl Default for SteinSteinParams {
    fn default() -> Self {
        Self {
            enabled: true,
            ridge_depth_m: 2600.0,
            subsidence_rate: 365.0,
            asymptotic_depth_m: 5651.0,
            time_constant: 0.0278,
            crossover_age_ma: 20.0,
            depth_scale_m: 5000.0,
            age_to_ma: 0.667,
        }
    }
}
