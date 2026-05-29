//! Tunables for the rifting closure + split mechanism (Track D,
//! Issue #132).
//!
//! ## Defaults rationale (Stage E2 W7 analytical first-pass)
//!
//! - `thinning_rate = 1.0` (`K_rift`) — half the Davis-Suppe
//!   orogenic coupling. Per-step thinning `Δs ≈ 1.0 × 0.014 × 0.69
//!   ≈ 0.010` at boundary peak (typical Phase 1.1 `|v_rel · n̂|`).
//!   Cumulative over 75 steps `≈ 0.75` — reaches the
//!   `split_thickness_threshold = 0.7` thinning target (continental
//!   baseline `S̃ = 1.0` thinned by 30 %) at around step 43, well
//!   inside the time-binding constraint at step 75. Slower than the
//!   collision orogenic rate is consistent with the geological
//!   record (Atlantic full rifting ~150 Ma vs Tibet collision ~50 Ma).
//! - `split_time_threshold = 75` steps (~50 Ma at Phase 1.1
//!   `dt ≈ 0.67 Ma/step`) — exceeds the accretion
//!   `merge_time_threshold = 50` because rifting needs MORE
//!   sustained extension than collision needs sustained
//!   compression in our timescale framing.
//! - `split_thickness_threshold = 0.7` — McKenzie 1978 β
//!   stretching factor `β = 1.4` ⇔ 30 % crustal thinning before
//!   lithospheric integrity is lost. Continental baseline `S̃ = 1.0`
//!   thinned to `0.7` matches this threshold.
//! - `split_velocity_offset = 0.005` — half of the typical Phase
//!   1.1 plate-speed magnitude `0.01`. Gives the new plate a
//!   noticeable perpendicular separation velocity without being
//!   so fast that it dominates over its parent plate's drift.
//! - `plate_id_cap = 256` — practical limit on the number of
//!   distinct plates a run can carry. `u16` would allow many more
//!   but `256` is enough for the visual narrative (Earth has
//!   ~12-15 major plates, ~50 minor). Beyond this, splits are
//!   refused gracefully (W6 architectural surface for Stage A
//!   if hit).
//!
//! Calibration discipline per [[calibration-via-visual-review]]
//! tier 2 (analytical first-pass + visual review, max 3
//! iterations). 5th C1 first-shot calibration if Stage A
//! confirms the visible-event regime is reached.

#[derive(Clone, Copy, Debug)]
pub struct RiftingParams {
    /// Master enable/disable. When `false`, both
    /// `apply_rifting_thinning` and `apply_rifting_split` are
    /// no-ops (W4 closure-isolation discipline).
    pub enabled: bool,

    /// Thinning-rate coefficient `K_rift`. Per-cell per-step
    /// thinning is `Δs = thinning_rate × |v_rel · n̂| × dt` on
    /// divergent continental cells (mirror of Davis-Suppe
    /// orogenic source on convergent boundaries).
    pub thinning_rate: f64,

    /// Number of consecutive divergent steps on a plate pair
    /// required before the chewing-gum cut split fires
    /// (sustained-extension condition).
    pub split_time_threshold: usize,

    /// Minimum `S̃` value across the rift-strip cells required
    /// before the split fires (sub-threshold-thinness condition).
    /// Default `0.7` = continental baseline `1.0` thinned to
    /// McKenzie 1978's `β = 1.4` stretching factor cap.
    pub split_thickness_threshold: f64,

    /// Magnitude of the perpendicular velocity offset applied to
    /// the newly-spawned plate (Q3.4). Direction is computed at
    /// split-fire time as the perpendicular to the parent's
    /// `v_rel` with the rifted-off plate (right-hand-rule
    /// perpendicular for determinism).
    pub split_velocity_offset: f64,

    /// Cap on the total number of plates a run can carry. When
    /// `kinematics.velocities.len() >= plate_id_cap`, new splits
    /// are refused gracefully (W6 — surface architectural finding
    /// if Stage A approaches this cap).
    pub plate_id_cap: usize,
}

impl Default for RiftingParams {
    fn default() -> Self {
        Self {
            enabled: true,
            thinning_rate: 1.0,
            split_time_threshold: 75,
            split_thickness_threshold: 0.7,
            split_velocity_offset: 0.005,
            plate_id_cap: 256,
        }
    }
}
