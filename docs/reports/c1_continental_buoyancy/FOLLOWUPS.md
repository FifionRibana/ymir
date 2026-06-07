# Issue #145 — registered follow-ups (with acquired diagnosis)

Out of #145 scope (scope kept healthy: core buoyancy fix proven/implemented/flipped). Each carries its acquired diagnosis so it is actionable without re-diagnosis.

## 1. Rigid-boundary refinement (the "curtain" oscillation)

**Diagnosis (acquired):** a 1px mesh-to-mesh oscillation at the continental/oceanic boundary. Cause = **upwind advection on the sharp S̃ contrast (1.0 vs 0.2)** — a known upwind behaviour on a steep discontinuity, made PERMANENT by the rigidity fix (before, collapsing continents never held a sharp persistent boundary). The no-flux wall spreads it (count 648→1063) but **amplitude is BOUNDED on all seeds (capped by equilibrium-height at S̃≈2.1, global_max 2.2 constant)** — cosmetic, NOT a divergent instability, NOT false land. Subduction grid-aligned promotion also contributes a residual (now bounded by the ≥2 rule).
**Lead:** a **LOCAL buoyancy transition at the boundary** (1-2 cells) to soften the net contrast — NOT the global graded mobility (measured + rejected: it fragmented continents). Re-measure oscillation amplitude + spatial bar + visual.

## 2. Oceanic ridge crust-creation (seafloor spreading)

**Diagnosis (acquired):** the raw-S̃ deep oceans come from advection-divergence thinning oceanic crust toward convergence with **nothing recreating crust at divergent boundaries** (ridges). **Invisible in production** — Stein-Stein re-anchors oceanic depth from AGE, absorbing the S̃ swing (measured: production oceanic altitude coherent across the −34%..+26% mass swing). So this is **S̃-fidelity, not a production defect**.
**Lead:** add crust creation at divergent boundaries, tied to the **Stein-Stein / GDH1 relation** (literature), NOT a homemade rate.

## 3. Conservative erosion deposition

**Diagnosis (acquired):** Step 1 (clean non-injecting removal) is a HALF-fix — it removed the floor-clamp mass INJECTION but left erosion a pure SINK (the non-conservation moved from +injection to −deletion). 
**Lead:** route eroded mass downstream (drainage targets exist) → deposit as sediment. Restores mass conservation AND builds plains/lowlands (Living Landz city terrain). Builds on the clean base; possibly with Phase 3.

## 4. Init R7 seeds with zero cratonic cells

**Diagnosis (acquired):** seeds 99 & 1337 init with ZERO cratonic cells (R7 continental-clustering is stochastic). Continents still form (continental clustering), just no cratonic-core label → craton-area metric undefined for those seeds.
**Lead:** if a guaranteed cratonic core is desired, constrain the R7 init to always seed ≥1 craton; else document as expected variability.

## 5. Davis-Suppe `h_critical` is an empirical exponential, not critical-wedge mechanics (Phase 3 Lallemand)

**Diagnosis (acquired):** `h_critical(d) = h_max·(1 − exp(−d/l_taper))` (davis_suppe/source_term.rs) is an empirical exponential rise of plausible shape — NOT the Davis-Suppe-Dahlen critical-wedge mechanics (taper = α_surface + β_décollement). `l_taper` is the exponential's characteristic length, not the run-out at slope tan(α_c). Surfaced during the mesh-invariance scope: anchoring `l_taper` on `h_max/tan(α_c)` would dress a true-physics number onto a non-physics formula (closure-relations trap) — so the invariance fix uses a UNIT conversion of the current 64² value instead, deferring fidelity here.
**Lead:** refound `h_critical` on the true critical-wedge physics (a closure relation, not the homemade exponential) at **Phase 3 Lallemand**, where the code already defers fine fidelity ("Phase 3 Lallemand will refine to the true relative-velocity normal component"). Distinct from invariance: re-founding changes behaviour → reopens closure calibration, which the invariance chantier must not do.

---

**Note — the mesh-invariance chantier itself** (S̃ thickening non-convergent, DS+accretion per-cell length scales) is now its OWN issue/branch (`c1-mesh-invariance`), not a #145 follow-up. See `stage_mesh_convergence.md` + `stage_invariance_scope.md`. Item 5 above is its registered downstream (Objective 2).
