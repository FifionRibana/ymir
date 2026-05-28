//! Subduction closure — oceanic-mass consumption + arc volcanism
//! + floor-triggered `plate_id` reassignment (C1 Phase 2 Track D,
//! Issue #132).
//!
//! ## Physics in one paragraph
//!
//! At an oceanic-continental convergent boundary, the denser
//! oceanic plate descends beneath the continental upper plate. The
//! subducted slab dehydrates as it sinks, releasing fluids that
//! lower the melting point of the overlying mantle wedge — magmas
//! rise and form a volcanic arc on the continental side, parallel
//! to and ~100-300 km inland of the trench (Lallemand 2005).
//!
//! C1 implements this with three coupled algorithmic effects per
//! step (per Track D Q1.2-Q1.4 decisions):
//!
//! ```text
//!   Δs_oceanic   = − K_subduction · |v_rel · n̂| · dt
//!   Δarc_total   = − Δs_oceanic · arc_efficiency
//!   distribution = BFS up to arc_distance cells, equal-share on
//!                  continental cells reached
//! ```
//!
//! When the oceanic cell's `S̃` drops below
//! `plate_id_reassign_threshold`, the cell is reassigned to the
//! adjacent continental plate (Q1.4 floor-triggered reassignment
//! — the "subducted" state) and its `plate_type` is promoted to
//! `Continental`.
//!
//! ## Track D position in the closure stack
//!
//! Subduction is the **first C1 closure to mutate `plate_id` and
//! `plate_type`** (Phase 1.1-1.4 + Track A/B treat both as
//! static-after-init). This breaks the static-classification
//! optimisation in `tectonics_c1::time_loop::run_with_closures`:
//! `classify_boundaries` and the Davis-Suppe `wedge_distance` are
//! precomputed ONCE outside the loop. The Track D integration
//! point (Stage E4) will move those into the per-step block under
//! the Track D enabled flag.
//!
//! ## Track D Q1.x decisions (from Issue #132 design pass)
//!
//! - **Q1.2 — rate-based consumption** (vs event-based step
//!   trigger). The closure runs each step on every convergent
//!   oceanic-continental cell; consumption is a continuous rate
//!   proportional to `|v_rel · n̂|`. Rationale: matches the
//!   continuous-process pattern of Phase 1.2 (Davis-Suppe rate
//!   source), Phase 1.3 (equilibrium-height rate sink), Phase 1.4
//!   (stream-power rate sink).
//! - **Q1.3 — local arc distribution** (BFS up to `arc_distance`,
//!   vs global gradient field). Rationale: the volcanic-arc
//!   process is geometrically local (~100-300 km inland of trench);
//!   `arc_distance = 3` cells at 64² maps to ~50-150 km depending
//!   on `dx` interpretation, in-range. A global gradient would
//!   smear arc deposition over the entire continental interior,
//!   destroying the trench-parallel signature.
//! - **Q1.4 — floor-triggered `plate_id` reassignment** (vs
//!   time-window threshold). Rationale: the physical event is
//!   "the oceanic crust has been thinned below sea-level
//!   threshold by accumulated consumption" — a state predicate on
//!   `S̃`, not a memory of how many steps ago consumption started.
//!   Matches the cyclical `apply_post_tectonic` reclassification
//!   pattern (`S̃ vs sea_level_ref` → `plate_type`).
//!
//! ## Calibration discipline
//!
//! Per `feedback_calibration_via_visual_review` tier 2 (analytical
//! first-pass + visual review, max 3 iterations). The Stage S W6
//! first-pass for `consumption_rate = 0.5` is documented in
//! [`params::SubductionParams`]. Stage A visual review may iterate
//! if event frequency falls outside the visible-event regime.
//!
//! ## Module layout
//!
//! - [`params`] — [`params::SubductionParams`] tunables
//!   (`consumption_rate = 0.5`, `arc_efficiency = 0.5`,
//!   `arc_distance = 3`, `plate_id_reassign_threshold = 0.05`,
//!   `enabled = true`).
//! - [`source_term`] — [`source_term::apply_subduction_step`] +
//!   private `distribute_arc_mass` BFS helper + 6 unit tests.
//!   Returns [`source_term::SubductionStats`] for the Stage E4
//!   mass-balance diagnostic.
//!
//! ## References
//!
//! - Lallemand, S., Heuret, A. & Boutelier, D. (2005). On the
//!   relationships between slab dip, back-arc stress, upper plate
//!   absolute motion, and crustal nature in subduction zones.
//!   *Geochem. Geophys. Geosyst.* 6(9), Q09006.
//!   doi:10.1029/2005GC000917
//! - Syracuse, E. M. & Abers, G. A. (2006). Global compilation of
//!   variations in slab depth beneath arc volcanoes and
//!   implications. *Geochem. Geophys. Geosyst.* 7(5), Q05017.
//!   doi:10.1029/2005GC001045

pub mod params;
pub mod source_term;

pub use params::SubductionParams;
pub use source_term::{apply_subduction_step, SubductionStats};
