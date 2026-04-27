//! Step 10 — geological age field `A(x, t)`.
//!
//! Passive scalar field tracking the geological age (in units of
//! `τ*`) of crustal material at each cell. Advected by plate
//! motion, reset by boundary events (ridge → 0, arc → 0, collision
//! → max-of-contributing-cells). Does **not** feed back into the
//! Stokes operator — `A` is pure state, with no influence on
//! viscosity, yield stress, or any rheological term. Compute cost
//! is one extra advection sweep + one extra event pass per
//! timestep.
//!
//! Enables downstream phases (differential erosion, climate /
//! soil maturity, geological classification export) to
//! distinguish "young Alpine analogue" from "ancient Appalachian
//! analogue" with a single quantitative parameter.
//!
//! See `solver-scaling.md` §4.11 for the full design and the Step
//! 10 issue (`step10.md`) for the milestone-scope decisions.

pub mod advection;
pub mod events;
pub mod init;

use crate::tectonics_v2::field::Field2D;

/// Concrete parameters for `AgeFieldConfig::Enabled`. Both fields
/// are nondimensional ages in units of `τ*` (the characteristic
/// time scale; see `scales::Scales::tau_star`).
#[derive(Clone, Copy, Debug)]
pub struct AgeFieldConfigEnabled {
    /// Initial age for continental cells (`S̃ > 0.5`). Per §4.11
    /// "continental cells start at a large value (representing
    /// pre-simulation cratonic age, e.g. `A_0 = 5–10·τ*`)".
    /// Default: `7.0`. Range: `[3, 15]`.
    pub continental_age_init: f64,
    /// Initial age for oceanic cells (`S̃ ≤ 0.5`). Per §4.11
    /// "oceanic cells at a smaller value reflecting ridge
    /// proximity". Default: `0.5`. Range: `[0, 2]`.
    pub oceanic_age_init: f64,
}

impl AgeFieldConfigEnabled {
    pub const CONTINENTAL_AGE_INIT_DEFAULT: f64 = 7.0;
    pub const OCEANIC_AGE_INIT_DEFAULT: f64 = 0.5;
}

impl Default for AgeFieldConfigEnabled {
    fn default() -> Self {
        Self {
            continental_age_init: Self::CONTINENTAL_AGE_INIT_DEFAULT,
            oceanic_age_init: Self::OCEANIC_AGE_INIT_DEFAULT,
        }
    }
}

/// Geological-age-field configuration (Step 10).
///
/// `Disabled` short-circuits the entire pipeline: no `A` field is
/// allocated, no advection of `A` runs, no event reset fires, no
/// metrics are sampled. This is the path used by the Step 9
/// regression to verify bit-identical output with pre-Step-10
/// code.
///
/// `Enabled(cfg)` activates the full age-field pipeline. The `A`
/// field is allocated at simulation init from the initial `S̃`
/// classification (continental vs oceanic), advected each
/// timestep using the non-conservative upwind scheme of
/// [`advection::step_age_advect`], and reset at boundary cells per
/// [`events::apply_age_events`].
#[derive(Clone, Copy, Debug)]
pub enum AgeFieldConfig {
    Disabled,
    Enabled(AgeFieldConfigEnabled),
}

impl AgeFieldConfig {
    pub fn label(&self) -> &'static str {
        match self {
            AgeFieldConfig::Disabled => "disabled",
            AgeFieldConfig::Enabled(_) => "enabled",
        }
    }

    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "disabled" | "off" => Ok(AgeFieldConfig::Disabled),
            "enabled" | "on" => Ok(AgeFieldConfig::Enabled(Default::default())),
            other => Err(format!(
                "unknown --age-field-config value '{}'; expected disabled|enabled",
                other,
            )),
        }
    }
}

impl Default for AgeFieldConfig {
    fn default() -> Self {
        AgeFieldConfig::Disabled
    }
}

/// Run-local state for the age field. Owned by the harness for the
/// duration of the run; the harness advects `current` each step
/// using `next` as scratch and swaps them at the end of the step.
///
/// Allocated only under `AgeFieldConfig::Enabled`. The `Disabled`
/// path holds `Option<AgeFieldState>` as `None` and pays no cost.
#[derive(Clone)]
pub struct AgeFieldState {
    /// Current age field, indexed `j * nx + i`. Units of `τ*`.
    pub current: Field2D,
    /// Scratch buffer for the next timestep's advection output.
    /// Swapped with `current` after each `step_age_advect` call.
    pub next: Field2D,
    /// `continental_age_init` value used at init. Captured for
    /// diagnostics + the bound-acceptance check
    /// `age_max ≤ age_init_max + simulation_time`.
    pub continental_age_init: f64,
    /// `oceanic_age_init` value used at init.
    pub oceanic_age_init: f64,
}

impl AgeFieldState {
    /// Build the initial `A` field from the initial `S̃` field per
    /// the §4.11 / D2 / D7 spec: continental cells (`S̃ > 0.5`)
    /// start at `continental_age_init`, oceanic cells at
    /// `oceanic_age_init`. Static identification — the
    /// classification is performed once and stored; cells that
    /// later transition continental ↔ oceanic via S̃ advection /
    /// recycling carry their advected `A` value, modulated by
    /// boundary-event resets.
    pub fn from_initial_thickness(s_initial: &Field2D, cfg: &AgeFieldConfigEnabled) -> Self {
        let nx = s_initial.nx();
        let ny = s_initial.ny();
        let mut current = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let s = s_initial.get(i, j);
                let a = if s > 0.5 { cfg.continental_age_init } else { cfg.oceanic_age_init };
                current.set(i, j, a);
            }
        }
        let next = Field2D::new(nx, ny);
        Self {
            current,
            next,
            continental_age_init: cfg.continental_age_init,
            oceanic_age_init: cfg.oceanic_age_init,
        }
    }

    /// Maximum `A` value attained by the initial state — used by
    /// the bound-acceptance check
    /// `age_max ≤ age_init_max + simulation_time` (acceptance #1).
    pub fn age_init_max(&self) -> f64 {
        self.continental_age_init.max(self.oceanic_age_init)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_default_is_disabled() {
        match AgeFieldConfig::default() {
            AgeFieldConfig::Disabled => {}
            _ => panic!("default should be Disabled"),
        }
    }

    #[test]
    fn enabled_defaults_match_issue_spec() {
        let cfg = AgeFieldConfigEnabled::default();
        assert_eq!(cfg.continental_age_init, 7.0);
        assert_eq!(cfg.oceanic_age_init, 0.5);
    }

    #[test]
    fn parse_roundtrip() {
        match AgeFieldConfig::parse("disabled").unwrap() {
            AgeFieldConfig::Disabled => {}
            _ => panic!(),
        }
        match AgeFieldConfig::parse("enabled").unwrap() {
            AgeFieldConfig::Enabled(_) => {}
            _ => panic!(),
        }
        assert!(AgeFieldConfig::parse("garbage").is_err());
    }

    #[test]
    fn from_initial_thickness_classifies_continental_vs_oceanic() {
        // Build a small S̃ field with one continental cell and the
        // rest oceanic. Verify the initial A field uses the right
        // age per cell.
        let mut s = Field2D::filled(4, 4, 0.2); // oceanic
        s.set(2, 2, 1.0); // continental at (2,2)
        let cfg = AgeFieldConfigEnabled::default();
        let state = AgeFieldState::from_initial_thickness(&s, &cfg);
        for j in 0..4 {
            for i in 0..4 {
                let expected = if (i, j) == (2, 2) { 7.0 } else { 0.5 };
                assert_eq!(state.current.get(i, j), expected, "cell ({},{})", i, j);
            }
        }
    }

    #[test]
    fn age_init_max_returns_continental() {
        let cfg = AgeFieldConfigEnabled::default();
        let s = Field2D::filled(2, 2, 1.0);
        let state = AgeFieldState::from_initial_thickness(&s, &cfg);
        assert_eq!(state.age_init_max(), 7.0);
    }
}
