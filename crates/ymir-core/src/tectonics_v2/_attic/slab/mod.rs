//! Slab-pull body force with regularised slab-mass ODE — Step 7.
//!
//! Installs the first **driving** mechanism in the `tectonics_v2`
//! pipeline: a per-cell body force `f̃_slab = Sp · m̃ · n̂_convergence`
//! (§4.8 of `docs/solver-scaling.md`) driven by a slab-mass field
//! `m̃(x, t)` that accumulates from subduction convergence and
//! decays with a characteristic time `τ̃_slab`. The exponential
//! decay replaces the legacy hard `max_plate_velocity` clamp — the
//! physics self-limits.
//!
//! Layout
//! ------
//!
//! - [`state::SlabState`] owns the `m̃` field; it exposes the ODE
//!   step (`step_ode`) and the upwind advection (`advect`). The
//!   ODE step runs **before** the Stokes solve (so `f_slab` sees
//!   the updated `m`); advection runs **after** (to use the
//!   solved velocity).
//! - [`accumulation::compute_q_sub_conv`] computes the source
//!   term `Q̃_sub_conv = k_slab_accum · max(0, -div(v))` on
//!   oceanic subducting cells. Distinct from `Q_sub` (which drains
//!   `S̃`), `k_slab_accum` is its own rate to keep `Sp` in the
//!   §4.8 target band `[0.5, 3.0]` without entangling the slab
//!   accumulation with the recycling fractions.
//! - [`convergence_direction::compute_convergence_direction`]
//!   produces `n̂_convergence = -∇(div(v)) / |∇(div(v))|` per cell,
//!   with a zero-vector fallback below `epsilon` to avoid
//!   ill-defined unit vectors in quiescent regions.
//! - [`crate::tectonics_v2::forcing::slab_pull::SlabPullForce`] (in the existing
//!   `forcing/` hierarchy) implements the `BodyForce` trait so the
//!   contribution plugs directly into `ForceSum` alongside `GpeForce`.
//!
//! Disabled bypass
//! ---------------
//!
//! `SlabPullConfig::Disabled` structurally short-circuits the whole
//! pipeline: `SlabState` is not allocated, no `Q_sub_conv` or
//! `n̂_convergence` are computed, and `SlabPullForce` is not pushed
//! into the `ForceSum`. The Step 7 regression test exercises this
//! path against the Step 6 physics baseline to enforce the
//! zero-cost-when-disabled invariant.

pub mod accumulation;
pub mod convergence_direction;
pub mod state;

pub use accumulation::{AccumulationConfig, compute_q_sub_conv};
pub use convergence_direction::{ConvergenceDirectionConfig, compute_convergence_direction};
pub use state::SlabState;

/// Default slab-pull coupling. Sits in the middle of the §4.8
/// target band `[0.5, 3.0]`.
pub const SP_DEFAULT: f64 = 1.5;

/// Default slab-mass decay time, nondim. `τ_slab ≈ 0.5 · τ* ≈
/// 15 Myr` with the baseline scales.
pub const TAU_SLAB_DEFAULT: f64 = 0.5;

/// Default slab-accumulation rate. `k_slab_accum = 1.0` is the
/// "absorbed into `Sp`" convention — decoupling the slab mass
/// source from the `Q_sub` drain rate, see D3 in the Step 7 spec.
pub const K_SLAB_ACCUM_DEFAULT: f64 = 1.0;

/// Default gradient epsilon for `n̂_convergence`. Above machine
/// noise on `∇(div v)`, well below the significant-gradient
/// magnitudes observed once slab-pull is active.
pub const EPSILON_DEFAULT: f64 = 1.0e-6;

/// Slab-pull configuration.
///
/// `Disabled` bypasses the whole pipeline (D7). `Enabled` carries
/// the four knobs; sensible defaults live in the `*_DEFAULT`
/// constants above.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SlabPullConfig {
    Disabled,
    Enabled {
        /// §4.8 slab-pull coupling `Sp ∈ [0.5, 3.0]`.
        sp: f64,
        /// Slab-mass decay time `τ_slab` (nondim). Stability
        /// requires `Δt < τ_slab`; the harness enforces
        /// `Δt ≤ 0.1 · τ_slab` defensively.
        tau_slab: f64,
        /// Accumulation rate driving `m_subducted` (D3).
        k_slab_accum: f64,
        /// Fallback threshold on `|∇(div v)|` below which
        /// `n̂_convergence` is set to the zero vector.
        epsilon: f64,
    },
}

impl SlabPullConfig {
    /// Short stable label for reports and logs.
    pub fn label(&self) -> &'static str {
        match self {
            SlabPullConfig::Disabled => "disabled",
            SlabPullConfig::Enabled { .. } => "enabled",
        }
    }

    /// Parse a CLI enable/disable token. Numerical knobs come
    /// through dedicated `--sp`, `--tau-slab`, … flags.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "disabled" | "off" => Ok(SlabPullConfig::Disabled),
            "enabled" | "on" => Ok(SlabPullConfig::Enabled {
                sp: SP_DEFAULT,
                tau_slab: TAU_SLAB_DEFAULT,
                k_slab_accum: K_SLAB_ACCUM_DEFAULT,
                epsilon: EPSILON_DEFAULT,
            }),
            other => Err(format!(
                "unknown --slab-pull-config value '{}'; expected disabled|enabled",
                other,
            )),
        }
    }

    /// Human-readable descriptor threaded into `SolverConfigDump`.
    pub fn describe(&self) -> String {
        match self {
            SlabPullConfig::Disabled => "Disabled".to_string(),
            SlabPullConfig::Enabled { sp, tau_slab, k_slab_accum, epsilon } => format!(
                "Enabled (Sp = {:.3}, τ_slab = {:.3}, k_slab_accum = {:.3}, ε = {:.1e})",
                sp, tau_slab, k_slab_accum, epsilon,
            ),
        }
    }

    /// Convenience: `Sp` when enabled, `0.0` when disabled.
    /// Useful for reporting without match-destructuring at every
    /// call site.
    pub fn sp_or_zero(&self) -> f64 {
        match self {
            SlabPullConfig::Disabled => 0.0,
            SlabPullConfig::Enabled { sp, .. } => *sp,
        }
    }

    /// Convenience: `τ_slab` when enabled, `f64::INFINITY` when
    /// disabled (the disabled case has no decay time).
    pub fn tau_slab(&self) -> f64 {
        match self {
            SlabPullConfig::Disabled => f64::INFINITY,
            SlabPullConfig::Enabled { tau_slab, .. } => *tau_slab,
        }
    }
}

impl Default for SlabPullConfig {
    fn default() -> Self {
        SlabPullConfig::Disabled
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_disabled() {
        assert_eq!(SlabPullConfig::default(), SlabPullConfig::Disabled);
    }

    #[test]
    fn parse_disabled_and_enabled() {
        assert_eq!(SlabPullConfig::parse("disabled").unwrap(), SlabPullConfig::Disabled);
        assert_eq!(SlabPullConfig::parse("off").unwrap(), SlabPullConfig::Disabled);
        match SlabPullConfig::parse("enabled").unwrap() {
            SlabPullConfig::Enabled { sp, tau_slab, k_slab_accum, epsilon } => {
                assert_eq!(sp, SP_DEFAULT);
                assert_eq!(tau_slab, TAU_SLAB_DEFAULT);
                assert_eq!(k_slab_accum, K_SLAB_ACCUM_DEFAULT);
                assert_eq!(epsilon, EPSILON_DEFAULT);
            }
            SlabPullConfig::Disabled => panic!("expected Enabled"),
        }
    }

    #[test]
    fn parse_rejects_unknown() {
        assert!(SlabPullConfig::parse("banana").is_err());
    }

    #[test]
    fn sp_defaults_sit_in_scaling_band() {
        // §4.8: Sp ∈ [0.5, 3.0]. Default must sit in this band.
        assert!(SP_DEFAULT >= 0.5 && SP_DEFAULT <= 3.0);
        // τ_slab ∈ [0.3, 1.0] per §4.8.
        assert!(TAU_SLAB_DEFAULT >= 0.3 && TAU_SLAB_DEFAULT <= 1.0);
    }

    #[test]
    fn sp_or_zero_and_tau_helpers() {
        assert_eq!(SlabPullConfig::Disabled.sp_or_zero(), 0.0);
        assert!(SlabPullConfig::Disabled.tau_slab().is_infinite());
        let en =
            SlabPullConfig::Enabled { sp: 2.0, tau_slab: 0.4, k_slab_accum: 1.0, epsilon: 1e-6 };
        assert_eq!(en.sp_or_zero(), 2.0);
        assert_eq!(en.tau_slab(), 0.4);
    }
}
