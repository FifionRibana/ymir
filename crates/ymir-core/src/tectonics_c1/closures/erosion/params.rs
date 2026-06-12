//! Tunables for the stream-power erosion closure.
//!
//! See [`super`] for the physics derivation, the calibration
//! discipline (visual review, not literature), and the
//! interaction with Phase 1.3's equilibrium-height sink.

/// Stream-power erosion closure tunables.
///
/// Defaults selected for the Phase 1.4 64²×300-step demo with
/// both Davis-Suppe (Phase 1.2) and equilibrium-height (Phase
/// 1.3) closures active.
///
/// - `k = 0.001` — first-pass analytical estimate; Stage E3
///   visual review may adjust within a 3-iteration budget per
///   the design doc §11.1 calibration discipline. Range
///   typically `[1e-4, 1e-2]`; below `1e-4` erosion is
///   negligible, above `1e-2` it erases the Davis-Suppe wedges
///   entirely.
/// - `m = 0.5` — Whipple-Tucker 1999 constraint `m / n ≈ 0.5`
///   (W-T eq. 4).
/// - `n = 1.0` — canonical W-T value (linear in slope,
///   uplift-magnitude-independent response timescale). Lague
///   2014 suggests `n ≈ 2` empirically; staged as a future
///   `Stage X.bis` upgrade analogous to Phase 1.3 E1.bis.
/// - `floor = 0.2` — matches the oceanic-init `S̃` baseline.
///   Continental cells eroded below this become oceanic via the
///   subsequent reclassification step in `apply_post_tectonic`
///   (`s > sea_level_ref` → continental, else oceanic).
///   Defensive against pathological `k` calibrations.
///
/// `enabled` follows the same W4 watchpoint convention as
/// [`crate::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams`]
/// and
/// [`crate::tectonics_c1::closures::equilibrium_height::params::EquilibriumHeightParams`]:
/// when `false`, the closure is a no-op and the run reproduces
/// the caller's previous behaviour bit-identically.
#[derive(Clone, Copy, Debug)]
pub struct ErosionParams {
    /// Master enable/disable. When `false`,
    /// [`super::source_term::apply_erosion_step`] is a no-op
    /// (W4 closure-isolation discipline).
    pub enabled: bool,
    /// Erosion coefficient — calibrated for visual balance per
    /// Lague 2014's framework (`K` is not universal across
    /// geological contexts). Units consistent with the
    /// non-dimensional time-loop `dt`; see [`super`] §
    /// "Parameter choices" for the first-pass analytical
    /// estimate.
    pub k: f64,
    /// Drainage-area exponent. W-T 1999 constraint `m / n ≈ 0.5`.
    pub m: f64,
    /// Slope exponent. W-T canonical `n = 1`; Lague 2014
    /// suggests `n ≈ 2` empirically — see [`super`] for the
    /// `Stage X.bis` upgrade marker.
    pub n: f64,
    /// Lower bound on `S̃` after the erosion step. Set to the
    /// oceanic initialisation baseline (`0.2`) so continental
    /// cells eroded below the floor naturally reclassify to
    /// oceanic in the next `apply_post_tectonic` pass.
    /// Defensive against pathological `k · A^m · S^n · dt > S̃`
    /// configurations.
    pub floor: f64,
    /// #155 A′ — craton erosion-resistance. `K` is multiplied by this
    /// factor on cells flagged cratonic, so ancient stable cratons erode
    /// SLOWER (their defining property alongside thick crust). Default
    /// `1.0` = OFF / byte-identical (v2, unit tests, generic erosion). The
    /// canonical C1 config (`C1Closures::default`) sets it < 1 (anchored in
    /// the real craton/non-craton erosion-rate band ~3-10×); without it an
    /// init-thick craton is planed (and can INVERT) by normal erosion — the
    /// measured 2026 craton-S̃ inversion. Pairs with the init thickness
    /// differential (`Phase2InitParams::craton_thickness_ratio`).
    pub craton_resist: f64,
}

impl Default for ErosionParams {
    fn default() -> Self {
        Self {
            enabled: true,
            k: 0.001,
            m: 0.5,
            n: 1.0,
            floor: 0.2,
            craton_resist: 1.0, // OFF by default — byte-identical
        }
    }
}
