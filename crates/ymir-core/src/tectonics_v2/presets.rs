//! Parameter presets for the `tectonics_v2` solver.
//!
//! At Step 1 the only physical mechanism active is power-law
//! rheology, so a preset carries only [`RheologyParams`] and a
//! startup [`ContinuationConfig`]. Fields for the other nondimensional
//! numbers (`Ar`, `Bi`, `De_p`, `Br`, `K`, `Sp`, `Mf`) will be added
//! step by step as the corresponding mechanisms come online.
//!
//! The three presets proposed in `docs/solver-scaling.md` §5.3
//! (`dynamic-accidented`, `stable-shield`, `soft-planet`) are
//! declared here. At Step 1 their rheology parameters are
//! **identical** — the design note pins `n = 1 → 3 via continuation`
//! without a per-preset variation, and the real differentiation
//! between presets rides on numbers introduced later (Ar for
//! spreading, Bi for yielding, `De_p` for plastic memory). Declaring
//! them as distinct names now keeps the CLI stable and avoids a
//! rename when Step 2+ populates their differences.

/// Power-law rheology parameters.
#[derive(Clone, Copy, Debug)]
pub struct RheologyParams {
    /// Power-law exponent `n`. Shear-thinning for `n > 1`.
    pub n: f64,
    /// Prefactor `B̃` in η_newton = B̃ · (ε̇_II + ε̇_min)^(1/n - 1).
    /// Equal to 1 after nondimensionalization.
    pub b_prefactor: f64,
    /// Additive strain-rate floor `ε̇_min` regularising the
    /// divergence of `η` at `ε̇_II → 0`. Range in design note:
    /// `[10⁻³, 10⁻²]`.
    pub strain_rate_floor: f64,
    /// Soft cap `η_max` applied through `smooth_saturate`. Not a hard
    /// clamp — the effective viscosity approaches `η_max` as an
    /// asymptote.
    pub eta_max_cap: f64,
    /// Sharpness exponent of `smooth_saturate`. The documented
    /// default `k = 4` is reused from the legacy implementation.
    pub k_saturation: f64,
}

impl RheologyParams {
    pub fn step1_default() -> Self {
        Self {
            n: 3.0,
            b_prefactor: 1.0,
            strain_rate_floor: 1.0e-3,
            eta_max_cap: 1.0e3,
            k_saturation: 4.0,
        }
    }
}

/// Startup-only continuation schedule on `n`.
///
/// At `t = 0` the nonlinear solver is run once for each value in
/// `n_steps` in order, using the previous solution as initial guess
/// for the next value. After the last entry, `n` is held at
/// `n_steps.last()` for the rest of the run.
#[derive(Clone, Debug)]
pub struct ContinuationConfig {
    pub n_steps: Vec<f64>,
}

impl ContinuationConfig {
    /// Matches the documented `ContinuationConfig::default` of the
    /// legacy solver: `1.0 → 3.0` in 0.5 increments.
    pub fn step1_default() -> Self {
        Self { n_steps: vec![1.0, 1.5, 2.0, 2.5, 3.0] }
    }
}

/// A named preset packages a rheology and a continuation schedule.
/// Future steps will add fields for Ar, Bi, etc.
#[derive(Clone, Debug)]
pub struct Preset {
    pub name: String,
    pub rheology: RheologyParams,
    pub continuation: ContinuationConfig,
}

impl Preset {
    pub fn dynamic_accidented() -> Self {
        Self {
            name: "dynamic-accidented".into(),
            rheology: RheologyParams::step1_default(),
            continuation: ContinuationConfig::step1_default(),
        }
    }

    pub fn stable_shield() -> Self {
        Self {
            name: "stable-shield".into(),
            rheology: RheologyParams::step1_default(),
            continuation: ContinuationConfig::step1_default(),
        }
    }

    pub fn soft_planet() -> Self {
        Self {
            name: "soft-planet".into(),
            rheology: RheologyParams::step1_default(),
            continuation: ContinuationConfig::step1_default(),
        }
    }

    pub fn by_name(name: &str) -> Result<Self, String> {
        match name {
            "dynamic-accidented" => Ok(Self::dynamic_accidented()),
            "stable-shield" => Ok(Self::stable_shield()),
            "soft-planet" => Ok(Self::soft_planet()),
            other => Err(format!(
                "unknown preset '{}'; expected one of: dynamic-accidented, stable-shield, soft-planet",
                other
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_three_presets_are_identical_at_step1() {
        // Explicit guard: prevents the step-1 assumption ("rheology
        // params identical across presets") from silently drifting.
        let a = Preset::dynamic_accidented();
        let b = Preset::stable_shield();
        let c = Preset::soft_planet();
        let same = |p: &RheologyParams, q: &RheologyParams| -> bool {
            p.n == q.n
                && p.b_prefactor == q.b_prefactor
                && p.strain_rate_floor == q.strain_rate_floor
                && p.eta_max_cap == q.eta_max_cap
                && p.k_saturation == q.k_saturation
        };
        assert!(same(&a.rheology, &b.rheology));
        assert!(same(&b.rheology, &c.rheology));
        assert_eq!(a.continuation.n_steps, c.continuation.n_steps);
    }

    #[test]
    fn by_name_is_case_sensitive() {
        assert!(Preset::by_name("dynamic-accidented").is_ok());
        assert!(Preset::by_name("Dynamic-Accidented").is_err());
        assert!(Preset::by_name("no-such-preset").is_err());
    }
}
