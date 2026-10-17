//! Tunables for the subduction closure (Track D, Issue #132).
//!
//! Defaults selected per Phase 2 Track D Stage S analytical first-
//! pass (W6) at Phase 1.1 kinematics scale (`|v_plate| ≈ 0.01`,
//! `dt ≈ 0.69 non-dim/step`, 300-step run):
//!
//! - `consumption_rate = 0.5` (`K_subduction`): cumulative consumption
//!   over 300 steps `≈ 0.5 × 0.014 × 0.69 × 300 ≈ 1.4` per boundary
//!   cell (where `|v_rel · n̂| ≈ 0.014` is the typical convergence
//!   magnitude for plates moving at `±0.01` toward each other).
//!   Comparable to one continental S̃ unit (1.0), so an oceanic cell
//!   at baseline `S̃ = 0.2` is consumed in roughly 40-60 steps — in
//!   the visible-event regime without being so aggressive that
//!   boundaries deplete in a handful of steps.
//! - `arc_efficiency = 0.5`: half of the consumed mass is
//!   redistributed as arc volcanism on continental neighbours; the
//!   rest is lost to "deeper mantle" (out of model).
//! - `arc_distance = 3` (BFS depth in cells): at 64² grid with
//!   typical 16-cell plates, this reaches the first ~24 cells inside
//!   the continental neighbour, mimicking volcanic-arc proximity to
//!   the trench (Lallemand 2005).
//! - `plate_id_reassign_threshold = 0.05`: a quarter of the oceanic
//!   baseline `S̃ = 0.2` — below this, the cell is depleted enough
//!   that calling it "continental" reflects the arc-volcanism
//!   build-up better than oceanic.
//!
//! Calibration discipline per `feedback_calibration_via_visual_review`
//! tier 2 (analytical first-pass + visual review, max 3 iterations
//! in Stage A) — same tier as Phase 1.3 `k_collapse` and Phase 1.4
//! erosion `K`.

#[derive(Clone, Copy, Debug)]
pub struct SubductionParams {
    /// Master enable/disable. When `false`,
    /// `apply_subduction_step` is a no-op (W4 closure-isolation
    /// discipline — must reproduce upstream behaviour bit-identically
    /// when disabled).
    pub enabled: bool,

    /// Consumption-rate coefficient `K_subduction`. Per-cell
    /// consumption per step is `Δs = consumption_rate × |v_rel · n̂|
    /// × dt` where `v_rel · n̂` is the convergence-normal velocity
    /// component at the oceanic-continental edge.
    pub consumption_rate: f64,

    /// Fraction of consumed mass redistributed as arc volcanism on
    /// the continental side. `0.5` matches the Lallemand 2005
    /// observation that a substantial fraction of subducted material
    /// (water, volatiles, sediments) returns to the surface as arc
    /// magmatism while the remainder is recycled into the mantle.
    pub arc_efficiency: f64,

    /// BFS depth (in cells) over which `arc_mass` is distributed
    /// from the consuming oceanic cell into nearby continental
    /// cells. `3` reaches the first ~24 cells of the continental
    /// volume on a 64² grid (typical 16-cell plates).
    pub arc_distance: usize,

    /// Floor on oceanic-cell `S̃` below which the cell is
    /// reassigned to the adjacent continental plate. Below this the
    /// cell is considered "subducted" — its remaining S̃ is small
    /// enough that promoting to continental better reflects the
    /// arc-built-up state than keeping it oceanic.
    pub plate_id_reassign_threshold: f64,
}

impl Default for SubductionParams {
    fn default() -> Self {
        Self {
            enabled: true,
            consumption_rate: 0.5,
            arc_efficiency: 0.5,
            arc_distance: 3,
            plate_id_reassign_threshold: 0.05,
        }
    }
}

/// ADR Finding 175 -- **the ocean-ocean subduction path** (option 5; `docs/reports/relief_method/f174_margin/f174_option5_note.md`),
/// a flag carried by the cascade's continent profile and passed to C1 through `run_with_closures_terms_oo`; never part of
/// [`SubductionParams`] (whose `Debug` is in the tectonic cache key). Its presence also cuts Davis-Suppe's seeds at
/// ocean-ocean boundaries (F175-A7).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OoSubductionParams {
    /// D_init: the cumulative convergence a plate pair's ocean-ocean boundary needs before it subducts, in cells of 64².
    pub d_init_cells: f64,
    /// ADR Finding 178-A2 (a variant, measured only): the O-O arc on the overriding plate's cells at walk depths `lo..=hi` from
    /// the consumed cell (the middle 50 % of arc-trench distances, Syracuse & Abers 2006 p. 10) instead of within
    /// `arc_distance`; `None` is F175's rule.
    pub arc_band: Option<(usize, usize)>,
}

impl OoSubductionParams {
    /// 100 km of convergence (Stern & Gerya 2018 p. 44, after Hall et al. 2003) ÷ 7.5 (the layout's compression) = 13.3 km
    /// on the map = 2.13 cells of 64² (6.25 km).
    pub fn declared() -> Self {
        Self { d_init_cells: 100.0 / 7.5 / 6.25, arc_band: None }
    }
    /// D_init in C1's domain units (the domain is 1; a 64² cell is 1/64).
    pub fn d_init_domain(&self) -> f64 {
        self.d_init_cells / 64.0
    }
}
