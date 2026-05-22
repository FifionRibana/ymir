//! Tunables for the equilibrium height closure.
//!
//! See [`super`] for the physics motivation and the regime in
//! which the defaults below are well-conditioned.

/// Equilibrium height closure tunables.
///
/// Defaults selected to interact cleanly with the Phase 1.2 Davis-
/// Suppe demo at 64²:
///
/// - `h_eq = 2.0` sits **below** the Davis-Suppe plateau
///   `h_max = 2.5`, so the Phase 1.2 boundary pile-up
///   (`global_max ≈ 2297` from advection-into-boundary cells) and
///   the sparse near-boundary saturation tail (`wedge_p99 ≈ 5.83`)
///   are capped at `h_eq` after a few relaxation time scales. The
///   Phase 1.2 wedge bulk (`wedge_p95 ≈ 0.376`, `mean ≈ 1.574`)
///   is below `h_eq` and stays untouched, so the Davis-Suppe
///   fill-ratio imprint (`fill_near ≈ 0.778`, `fill_far ≈ 0.079`)
///   is preserved.
///
/// - `k_collapse = 1.0` with the Phase 1.1 timestep
///   `dt ≈ 0.69` gives `k_collapse · dt ≈ 0.69` — under 1, so the
///   linear-Euler scheme stays in the safe convex-combination
///   regime (no undershoot of `h_eq`). The [`super::source_term`]
///   module ships a defensive clamp that catches pathological
///   timesteps where `k_collapse · dt > 1`.
///
/// `enabled` follows the same W4 watchpoint convention as
/// [`crate::tectonics_c1::closures::davis_suppe::source_term::DavisSuppeParams`]:
/// when `false`, the closure is a no-op and the run reproduces the
/// caller's previous behaviour bit-identically.
#[derive(Clone, Copy, Debug)]
pub struct EquilibriumHeightParams {
    /// Master enable/disable. When `false`,
    /// [`super::source_term::apply_equilibrium_height_step`] is a
    /// no-op (W4 closure-isolation discipline).
    pub enabled: bool,
    /// Equilibrium thickness toward which thickened cells collapse.
    /// Set below `DavisSuppeParams::h_max` so the orogenic wedge
    /// plateau is the active cap on the simulated state.
    pub h_eq: f64,
    /// Relaxation rate of the one-sided collapse term, in
    /// `1 / time` units consistent with the time-loop `dt`. Keep
    /// `k_collapse · dt < 1` to stay in the well-conditioned
    /// forward-Euler regime; the source-term module clamps to
    /// `h_eq` if this is violated by a pathological caller.
    pub k_collapse: f64,
}

impl Default for EquilibriumHeightParams {
    fn default() -> Self {
        Self {
            enabled: true,
            h_eq: 2.0,
            k_collapse: 1.0,
        }
    }
}
