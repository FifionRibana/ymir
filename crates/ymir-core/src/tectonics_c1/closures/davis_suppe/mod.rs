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

pub mod critical_taper;
pub mod source_term;
