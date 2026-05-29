//! Tunables for the accretion mechanism (Track D, Issue #132).
//!
//! ## `merge_time_threshold` — Q3 revision (20 → 50)
//!
//! The original Issue #132 design exploration Q3 set
//! `merge_time_threshold = 20` steps as the sustained-convergence
//! gate before a plate merge fires. Stage E1's W7 analytical pass
//! revised this to **50 steps** for three reasons:
//!
//! 1. **Physical timescale alignment.** At Phase 1.1 `dt ≈ 0.67 Ma/
//!    step`, 50 steps ≈ 33 Ma of sustained convergence — consistent
//!    with collision-zone timescales in the geological record
//!    (Indian-Asian collision matured over ~50 Ma).
//! 2. **Spurious-merge suppression.** Phase 1.1's random-cycled
//!    kinematics produces brief convergent pulses on plate pairs
//!    that subsequently drift apart. A 20-step threshold is
//!    sensitive enough to fire on these transients and produce
//!    spurious merges. 50 steps requires more committed convergence.
//! 3. **Tunable headroom.** If Stage A event-rarity diagnostic
//!    shows fewer merges than the visual narrative needs, the
//!    threshold can be lowered within the [[calibration-via-visual-
//!    review]] 3-iteration budget. Starting conservative leaves
//!    headroom; starting permissive does not.
//!
//! Memory entry [[c1-phase-2-track-d-outcomes]] will record this
//! revision at Stage Final.

/// Strategy for combining the merged plates' velocities into the
/// surviving plate's new velocity.
///
/// Currently the only shipping variant is mass-weighted average
/// (Q2.4 Option A per Issue #132). Future variants (e.g., simple
/// arithmetic average for diagnostic ablation, or momentum-
/// conserving with `Δp = m_a · v_a + m_b · v_b` semantics if
/// physical fidelity becomes the priority) can be added as new
/// enum cases without breaking callers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VelocityMergeMethod {
    /// Post-merge velocity is the mass-weighted average:
    /// `v_new = (m_a · v_a + m_b · v_b) / (m_a + m_b)` where mass
    /// `m_p = Σ S̃ over cells of plate p`. Default — matches the
    /// Issue #132 Q2.4 decision.
    MassWeightedAverage,
}

#[derive(Clone, Copy, Debug)]
pub struct AccretionParams {
    /// Master enable/disable. When `false`, `apply_accretion_step`
    /// is a no-op (W4 closure-isolation discipline).
    pub enabled: bool,

    /// Number of consecutive convergent steps required on a plate
    /// pair before the merge fires. Default `50` revised from the
    /// original Q3 value `20` per the module docstring rationale.
    pub merge_time_threshold: usize,

    /// Strategy for combining merged-plate velocities. Default
    /// [`VelocityMergeMethod::MassWeightedAverage`] per Q2.4.
    pub velocity_merge_method: VelocityMergeMethod,
}

impl Default for AccretionParams {
    fn default() -> Self {
        Self {
            enabled: true,
            merge_time_threshold: 50,
            velocity_merge_method: VelocityMergeMethod::MassWeightedAverage,
        }
    }
}
