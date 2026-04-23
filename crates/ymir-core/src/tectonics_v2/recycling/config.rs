//! Step 6 recycling configuration — closed-mode budget fractions.
//!
//! In Open mode (Step 5), each source/sink term is driven by its own
//! rate (`k_sub`, `k_arc`, `k_spread`, ...). In Closed mode (Step 6),
//! only `k_sub` survives as a rate. The other terms are **fractions
//! of the recycled budget**:
//!
//! ```text
//!   M_sub_step = Σ_cells |Q_sub_cell| · Δt · dA
//!
//!   M_arc_step    = arc_fraction    · M_sub_step
//!   M_coll_v_step = coll_v_fraction · M_sub_step
//!   M_rift_v_step = rift_v_fraction · M_sub_step
//!   M_spread_step = spread_fraction · M_sub_step   (through delayed buffer)
//!   M_lost_step   = mantle_loss_fraction · M_sub_step (never redistributed)
//! ```
//!
//! Constraint: `arc + coll_v + rift_v + spread + mantle_loss = 1`
//! (validated at config load with tolerance `1e-9`). A strict check
//! `sum == 1.0` would falsely reject configurations that sum to
//! `0.9999999999999999` due to f64 rounding on 5 additions.

/// Fractions of the subduction budget recycled at each step.
///
/// Default values: `(arc, coll_v, rift_v, spread, mantle_loss) =
/// (0.15, 0.03, 0.02, 0.80, 0.00)` — full recycling with spreading
/// dominant, loosely matching the mass-balance targets in
/// `solver-scaling.md` §4.7.
#[derive(Debug, Clone, Copy)]
pub struct RecyclingConfig {
    pub arc_fraction: f64,
    pub coll_v_fraction: f64,
    pub rift_v_fraction: f64,
    pub spread_fraction: f64,
    pub mantle_loss_fraction: f64,
    pub mantle_delay_steps: usize,
}

impl Default for RecyclingConfig {
    fn default() -> Self {
        Self {
            arc_fraction: 0.15,
            coll_v_fraction: 0.03,
            rift_v_fraction: 0.02,
            spread_fraction: 0.80,
            mantle_loss_fraction: 0.0,
            mantle_delay_steps: 20,
        }
    }
}

impl RecyclingConfig {
    /// Validate that the fractions sum to 1 within `1e-9` absolute.
    ///
    /// Rejects both "forgot to set a fraction" errors (sum too far
    /// from 1) and "off-by-one typo" errors (sum > 1.0 or < 1.0),
    /// while tolerating f64 rounding on the 5-term addition.
    pub fn validate(&self) -> Result<(), RecyclingConfigError> {
        let sum = self.arc_fraction
            + self.coll_v_fraction
            + self.rift_v_fraction
            + self.spread_fraction
            + self.mantle_loss_fraction;
        let tol = 1e-9;
        if (sum - 1.0).abs() > tol {
            return Err(RecyclingConfigError::FractionsDoNotSumToOne {
                sum,
                tolerance: tol,
            });
        }
        // Each fraction must be ≥ 0 (negative fractions have no
        // physical meaning and would allow net mass creation).
        for (name, v) in [
            ("arc_fraction", self.arc_fraction),
            ("coll_v_fraction", self.coll_v_fraction),
            ("rift_v_fraction", self.rift_v_fraction),
            ("spread_fraction", self.spread_fraction),
            ("mantle_loss_fraction", self.mantle_loss_fraction),
        ] {
            if v < 0.0 {
                return Err(RecyclingConfigError::NegativeFraction {
                    name,
                    value: v,
                });
            }
        }
        if self.mantle_delay_steps == 0 {
            return Err(RecyclingConfigError::ZeroMantleDelay);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum RecyclingConfigError {
    FractionsDoNotSumToOne { sum: f64, tolerance: f64 },
    NegativeFraction { name: &'static str, value: f64 },
    ZeroMantleDelay,
}

impl std::fmt::Display for RecyclingConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::FractionsDoNotSumToOne { sum, tolerance } => write!(
                f,
                "recycling fractions must sum to 1.0 (±{:.0e}), observed {:.15}",
                tolerance, sum,
            ),
            Self::NegativeFraction { name, value } => write!(
                f,
                "recycling fraction '{}' must be ≥ 0, got {}",
                name, value,
            ),
            Self::ZeroMantleDelay => write!(
                f,
                "mantle_delay_steps must be ≥ 1 (ring buffer of size 0 is undefined)",
            ),
        }
    }
}

impl std::error::Error for RecyclingConfigError {}

/// Per-step accumulators for immediate-distribution fractions that
/// couldn't be distributed (no eligible cell at the step). Same
/// rollover semantics as the delayed buffer: mass held, ready to
/// emerge at the next step where eligibility arises.
#[derive(Debug, Clone, Copy, Default)]
pub struct ImmediateAccumulators {
    pub arc_pending: f64,
    pub coll_v_pending: f64,
    pub rift_v_pending: f64,
}

impl ImmediateAccumulators {
    pub fn sum(&self) -> f64 {
        self.arc_pending + self.coll_v_pending + self.rift_v_pending
    }

    /// Max of the three pending values — the "immediate pending" is
    /// whichever is largest. Normalised externally by
    /// `mean_subducted_per_step` for the diagnostic metric
    /// `immediate_pending_max`.
    pub fn max_pending(&self) -> f64 {
        self.arc_pending.max(self.coll_v_pending).max(self.rift_v_pending)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_validates() {
        RecyclingConfig::default().validate().unwrap();
    }

    #[test]
    fn fractions_summing_exactly_1_pass() {
        let cfg = RecyclingConfig {
            arc_fraction: 0.2,
            coll_v_fraction: 0.1,
            rift_v_fraction: 0.1,
            spread_fraction: 0.5,
            mantle_loss_fraction: 0.1,
            mantle_delay_steps: 10,
        };
        cfg.validate().unwrap();
    }

    #[test]
    fn fractions_with_f64_rounding_pass_within_tolerance() {
        // 0.15 + 0.03 + 0.02 + 0.80 = 1.0 but the floating-point
        // addition may give 0.9999999999999998 or 1.0000000000000002
        // depending on order. Exact equality fails; our 1e-9
        // tolerance should pass.
        let cfg = RecyclingConfig::default();
        cfg.validate().unwrap();
    }

    #[test]
    fn fractions_missing_reject() {
        let cfg = RecyclingConfig {
            arc_fraction: 0.5,
            coll_v_fraction: 0.0,
            rift_v_fraction: 0.0,
            spread_fraction: 0.0,
            mantle_loss_fraction: 0.0,
            mantle_delay_steps: 10,
        };
        assert!(matches!(
            cfg.validate(),
            Err(RecyclingConfigError::FractionsDoNotSumToOne { .. })
        ));
    }

    #[test]
    fn negative_fraction_rejected() {
        let cfg = RecyclingConfig {
            arc_fraction: -0.1,
            coll_v_fraction: 0.1,
            rift_v_fraction: 0.0,
            spread_fraction: 1.0,
            mantle_loss_fraction: 0.0,
            mantle_delay_steps: 10,
        };
        assert!(matches!(
            cfg.validate(),
            Err(RecyclingConfigError::NegativeFraction { .. })
        ));
    }

    #[test]
    fn zero_delay_rejected() {
        let cfg = RecyclingConfig {
            mantle_delay_steps: 0,
            ..RecyclingConfig::default()
        };
        assert!(matches!(
            cfg.validate(),
            Err(RecyclingConfigError::ZeroMantleDelay)
        ));
    }

    #[test]
    fn immediate_accumulators_max_pending_picks_largest() {
        let a = ImmediateAccumulators {
            arc_pending: 0.3,
            coll_v_pending: 0.5,
            rift_v_pending: 0.1,
        };
        assert_eq!(a.max_pending(), 0.5);
        assert!((a.sum() - 0.9).abs() < 1e-14);
    }
}
