//! Tunables for the equilibrium height closure.
//!
//! See [`super`] for the physics motivation, the quadratic
//! formula derivation, and the threshold behavior that the
//! defaults below are calibrated for.

/// Equilibrium height closure tunables.
///
/// Defaults selected for the Phase 1.2 / 1.3 64²×300-step demo:
///
/// - `h_eq = 2.0` sits **below** the Davis-Suppe plateau
///   `h_max = 2.5` and matches the observed Tibet crustal-
///   thickness ratio (~70 km plateau vs ~35 km normal crust).
///   See [`super`] § "Parameter `h_eq = 2.0` (phenomenological)"
///   for why this is a tunable target, not a derived equilibrium
///   value.
///
/// - `k_collapse = 2.0` is calibrated for the quadratic formula
///   derived from Molnar & Lyon-Caen 1988 eq. (2). The squared
///   excess produces a threshold behavior:
///
///   - **Small excess** (wedge body cells, < 1 above `h_eq`):
///     per-step decrement `k · excess² · dt` stays small —
///     wedge cells relax gradually toward `h_eq` without erasing
///     the Davis-Suppe fill-ratio imprint.
///   - **Large excess** (boundary outliers, `excess > 100` in
///     the Phase 1.2 advection-dominated regime): per-step
///     decrement overshoots `h_eq`; the safety clamp in
///     [`super::source_term::apply_equilibrium_height_step`]
///     holds the cell at `h_eq` — effectively a one-step cap on
///     outliers.
///
/// `enabled` follows the same W4 watchpoint convention as
/// [`crate::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams`]:
/// when `false`, the closure is a no-op and the run reproduces
/// the caller's previous behaviour bit-identically.
#[derive(Clone, Copy, Debug)]
pub struct EquilibriumHeightParams {
    /// Master enable/disable. When `false`,
    /// [`super::source_term::apply_equilibrium_height_step`] is a
    /// no-op (W4 closure-isolation discipline).
    pub enabled: bool,
    /// Equilibrium thickness toward which thickened cells
    /// collapse. Set below `DavisSuppeParams::h_max` so the
    /// orogenic wedge plateau is the active cap on the simulated
    /// state.
    pub h_eq: f64,
    /// Relaxation rate for the **quadratic** collapse term
    /// `k_collapse · max(0, S̃ − h_eq)²`. Units `1 / (length ·
    /// time)` consistent with the time-loop `dt`. Calibrated
    /// against Molnar-Lyon-Caen 1988 ΔPE_A; the safety clamp in
    /// [`super::source_term`] catches the large-excess outliers
    /// where one step would overshoot `h_eq`.
    pub k_collapse: f64,
}

impl Default for EquilibriumHeightParams {
    fn default() -> Self {
        Self {
            enabled: true,
            h_eq: 2.0,
            k_collapse: 2.0,
        }
    }
}
