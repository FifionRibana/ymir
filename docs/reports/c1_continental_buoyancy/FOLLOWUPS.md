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

## 6. FOUNDATION — C1 upwind advection is bulk mesh-non-convergent (Issue #147 outcome)

**Diagnosis (acquired, measured):** S̃ does not converge with mesh resolution because the **upwind advection scheme itself** is bulk mesh-dependent (numerical diffusion ∝ dx; init→64 r 0.90 collapses to advection-only r 0.045 over a run). Anchored by an init-convergence CONTROL (init IS convergent → advection destroys it), the contrast counterfactual γ (boundary-contrast smoothing does NOT help → not a boundary issue), and split B (subduction innocent; accretion couples via velocity-averaging, not a parameter). NO closure-PARAMETER fix works: DS physical-width (Fix #1, real units bug, kept, partial) and accretion step→time (degraded 128/256, reverted) are units-hygiene, not the lever. Equilibrium-height MASKS it (full-system r~0.51) but does not cure it. Production S̃ r ~0.51; geography converges (alt r ~0.87).
**Lead:** a higher-order / less-diffusive / flux-limited (or semi-Lagrangian) advection scheme for S̃ + age. **FOUNDATION defect** of the transport scheme — same depth class as continental-buoyancy, surfaced by the same measure-don't-assume method. Out of #147, which REDEFINED the problem from "closure calibration" to "advection scheme".

**GATING ANSWERED → DEFERRABLE (do NOT open the scheme issue).** Measured
(`stage_upscale_robustness.md`): 64²+upscale vs 256²+upscale = SAME world
(upscaled structure r 0.90 ≈ coarse-altitude r 0.88; ~26% land + largest
~0.96 both; visual = same continents/orogen/bathymetry, different fineness).
The upscale orients FBM by the coarse ALTITUDE slope (convergent, r~0.88),
NOT raw S̃ (r~0.51) — the non-convergence is laundered through isostasy +
Stein-Stein, invisible downstream (same as the DS mass swing). So the
advection-scheme milestone is NOT urgent; #6 stays registered but deferred.
**Precondition:** holds because the upscale reads altitude; revisit only if a
consumer reads the raw S̃ gradient. Original reasoning retained below.

**GATING QUESTION (now answered above — kept for the record):** do we actually NEED r→1? The reason to want S̃-field invariance is the upscale (consumes the S̃ gradient). But we have NOT measured whether r~0.51 actually perturbs the upscale downstream. Precedent: the mass swing LOOKED like a problem; Stein-Stein absorbed it (invisible in production). The upscale may be similarly robust to r~0.51 (stochastic FBM merely ORIENTED by the gradient — perhaps no need for a perfectly convergent gradient). **MEASURE FIRST** (next session): 64²+upscale vs 256²+upscale — coherent detail (upscale robust to r~0.51 → scheme milestone DEFERRABLE) or different worlds (upscale diverges → milestone NECESSARY)? Same "is the effect visible downstream?" logic that (correctly) avoided over-fixing DS. The scheme issue's urgency/necessity is decided by that measurement, not assumed.

---

**Note — the mesh-invariance chantier (#147)** delivered Fix #1 (DS physical-width, partial) + the diagnostic that redefined the problem (item 6). See `stage_mesh_convergence.md` + `stage_invariance_scope.md`. Item 5 = Objective 2 (DS fidelity, Phase 3 Lallemand).
