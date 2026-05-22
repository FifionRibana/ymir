//! Davis-Suppe critical taper orogenic profile closure (C1 Phase 1.2,
//! Issue #123).
//!
//! ## Physics in one paragraph
//!
//! Per Davis, Suppe & Dahlen 1983, *JGR* 88(B2), 1153-1172 and the
//! review by Dahlen 1990, *Annu. Rev. Earth Planet. Sci.* 18, 55-99:
//! a fold-and-thrust wedge at a convergent plate boundary self-
//! organises into a *critical taper* geometry. The surface slope `α`
//! and basal dip `β` satisfy a mechanical balance between basal
//! traction and gravitational driving force; the wedge grows by
//! accretion at the toe and thickens with distance from the
//! boundary, but is bounded by gravitational stability at the
//! critical taper angle. Higher basal friction or lower pore
//! pressure increases the critical taper; submerged wedges (high
//! `ρ_w/ρ` or high `λ`) are systematically lower-tapered than
//! subaerial ones.
//!
//! ## Module layout
//!
//! - [`critical_taper`] — the small-angle Davis-Suppe formula and
//!   its reproduction unit tests (Stage 1 of Issue #123). This is
//!   the **physics validation layer**: if the formula doesn't
//!   reproduce Davis 1983 sandbox values, nothing downstream can
//!   be trusted.
//! - `source_term` (Stage 4) — the per-step `∂S̃/∂t` application on
//!   upper-plate cells at convergent boundaries.
//!
//! ## Phase 1.2 scope reminder
//!
//! Quantitative fidelity to natural orogens (Taiwan, Barbados, …)
//! is **not** a goal (design doc §8.5: Ymir trades physical
//! fidelity for visual plausibility). The critical-taper formula
//! is calibrated against the controlled-input sandbox experiments
//! where μ, λ are measured directly; natural-orogen agreement is
//! tracked informatively in `critical_taper.rs` test output but
//! does not gate.
//!
//! ## Findings during Phase 1.2 (Issue #123)
//!
//! 1. **Architectural skip on Convergent cells.** The intra-plate
//!    Dijkstra (`distance_field::wedge_distance_intra_plate`)
//!    gives `d = 0` to every upper-plate boundary cell. Applying
//!    `h_critical(0) = 0` would produce `driving = h_crit − h <
//!    0` and thin the boundary instead of thickening the
//!    interior. `apply_davis_suppe_step` skips
//!    `BoundaryType::Convergent` cells explicitly; the
//!    architectural lock test `source_term_skips_boundary_cells`
//!    protects this against accidental refactor. Boundary cells
//!    accumulate via advection alone — Phase 1.4 erosion is the
//!    planned mass sink.
//!
//! 2. **Advection-dominated regime.** The Phase 1.1 hand-tuned
//!    kinematics preset produces an advection rate `≈ 32 ×` the
//!    source relaxation rate. The wedge body **bulk-drains**;
//!    only cells immediately adjacent to a Convergent boundary
//!    saturate close to `h_critical(d)`. See the per-bucket
//!    fill-ratio table in `docs/reports/c1_phase_1_2_davis_suppe/
//!    README.md` for the empirical signature.
//!
//! 3. **Acceptance metric: fill ratio.** Absolute-mean tests are
//!    silently regime-dependent. The Stage 5 acceptance test
//!    uses `fill_ratio = mean / h_critical(d_bucket_mid)` —
//!    regime-agnostic — instead of asserting an absolute-mean
//!    direction across distance buckets. The
//!    `mean(near)/mean(far) > 1.5` sub-assertion is **regime-
//!    tagged** (advection-dominated direction); Phase 1.4 / 2
//!    tests must re-evaluate the direction rather than copy.
//!
//! Memory entries — pairs with
//! `feedback_fill_ratio_regime_agnostic_metric` (transferable)
//! and `project_c1_phase_1_2_advection_dominated_regime`
//! (project-specific).

pub mod critical_taper;
pub mod source_term;
