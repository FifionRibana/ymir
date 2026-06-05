# C1 Milestone — Buoyant (non-subducting) continental crust in the advection

**Issue:** #145 · **Branch:** `145-c1-continental-buoyancy`
**Status:** scoped (diagnostic prototyping complete — cause proven, fix proven). Not yet implemented.
**Branch convention:** one branch per issue; integration target is the C1 integration branch, not `main`.
**Roadmap placement:** BEFORE Phase 3 (no point building subduction arcs / boundary evolution on crust that collapses) and BEFORE resuming piste 4 (rendering a correct field). This milestone re-founds C1 crustal transport — it sits *below* the advection scheme itself.

---

## Problem (proven, piste 4 investigation)

C1 advects crustal thickness `S̃` **uniformly, regardless of buoyancy**. Continental lithosphere is physically rigid/buoyant — it does not subduct or funnel into trenches — yet C1 transports it passively like a fluid. On the (globally balanced) Phase-1.1 velocity field, continents are spatially biased into convergence (continental ∇·v = −0.0176, cratonic −0.0265; divergence re-feeds the ocean, +0.0067), so the convergent flow sweeps continental crust into stagnation pile-ups and empties the interiors.

**No advection scheme fixes this** (proven by exhaustion on the diagnostic run — closures off, seed 42, 64², 300):

| scheme | craton area preserved | interior | mass |
|---|---:|---|---:|
| flux-form upwind (current) | 0% (→0.00) | evacuated | −0.0% |
| naive semi-Lagrangian | smeared to 0.60 | diffused | −41% |
| perfect PIC (interp-free) | 1% (4/561) | piled to a point | −0.0% |
| **rigid continent (v=0)** | **92% (517/561)** | **1.00→1.00** | **−0.00%** |

The fix is **making continental/cratonic crust advectively distinct (non-subducting)**, validated on the diagnostic: 92% craton area preserved, interiors untouched, compact dominant landmass (perim/area 0.38, largest 0.92), mass exact. See memory `project_c1_crust_buoyancy_root_cause`.

## Design points to cover

> **Stage S decided (audit + comparative measurement, see `stage_s.md`):** points 1, 3, 4 are resolved; point 2 carries the relocated "does subduction resolve the boundary" measurement.

1. **Rigidity mechanism — DECIDED: BINARY.** `v = 0` on `plate_type == Continental`, via **upstream v-masking** (zero `vx`/`vy` in/after `fill_velocity_field`; mass-conservative — flux-form face cancellation, Δmass −0.00%). Binary beats graded on the SPATIAL bar (RUN 1, closures off): craton area 92% vs 89%, **one compact dominant mass (largest 0.92) vs fragmented (0.43, 31 components)**, drift 0.4, mass exact. Graded's only edge (boundary 1.55 vs 1.99) is moot in production (subduction consumes the inflow, point 2) and it pays by shredding margins (thinned margins become mobile → cascade fragmentation). Hook + mask details in `stage_s.md` S1/S2.

2. **Boundary handling (the one caveat) — carries the production measurement.** Oceanic crust converging on a rigid continental margin accumulates (cont/ocean edge mean 0.60→1.99, localized; closures OFF). Physically subduction/accretion. Route the edge inflow through the existing Track-D subduction/accretion closures. **GATE (relocated from the point-1 deliberation):** Stage S established *in theory* that subduction reads the right inflow, but the prototype measured it **closures OFF**. Point 2 must **MEASURE, closures ON (production)**, that subduction actually drains the 1.99 edge accumulation. This is the point-2 acceptance gate, not a point-1 condition.

   **Point-2 acceptance is THREE-part — never the boundary scalar alone** (E1, by contrast, was verifiable purely numerically: flag-OFF byte-identity + flag-ON 92% craton area = already the spatial bar; the eye adds nothing there). Point 2 is the FIRST time the rigid continent runs WITH the full tectonics (all closures + Track D), so the scalar can lie (the "good number / catastrophic morphology" trap seen 6× this thread — boundary 0.60 says nothing about whether subduction *drains* the ocean cleanly vs *gnaws* the continent inland, nor whether rigidity+subduction interaction creates dentelures / holes / false coast):
   - **(a) scalar** — boundary edge drained (1.99 → ~0.60 with subduction on).
   - **(b) spatial bar** (`morphology.rs`) — craton area preserved, perim/area, n_components, largest-component, centroid drift, on the production run.
   - **(c) VISUAL** — render rigid-continent + closures ON (Viz-0.5 gallery, seed 42 + seed 2) and LOOK: credible continent (sharp coast, compact mass, no interaction artefacts), not merely "edge drained".

3. **Coherence with `plate_type` — SATISFIED BY CONSTRUCTION.** Rigidity = `plate_type==Continental`, recomputed where `plate_type` is rebuilt: once in gallery (static), **per step on the Track-D path** (aligned with the existing per-step boundary/velocity recompute, lines 481–506). Newly-promoted continental crust (subduction Oceanic→Continental) becomes rigid the next step; rifted continental stays rigid. `cratonic_mask` never mutates (Track D doesn't write it), so it is not the mask base.

4. **Determinism — CLEARED.** `step_upwind` + `fill_velocity_field` are serial (no rayon/sort); the mask derives from deterministic fields. v-masking is a serial element-wise scale → no new nondeterministic ordering. Invariant preserved.

   **BLOCKING PREREQUISITE (Stage S discovery) — port `morphology.rs` into #145 (Stage E1).** The spatial gate (`land_morphology`, piste-4 Stage E1, commit `d110cdd`) is the **acceptance invariant** of this milestone, but it lives ONLY on the `piste-4-gallery-production-morphology` branch, NOT on milestone. Cherry-pick / port it into #145 in E1 before any validation stage.

5. **Re-validation matrix** — this change **breaks bit-identity** for every C1 phase AND **invalidates the Phase 1.x / Phase 2 calibrations** (Davis-Suppe rates, erosion K, Stein-Stein anchors — all tuned against uniform advection). Plan the full re-validation: each phase's acceptance re-run and recalibrated where needed. This is the deepest re-validation in the project.

6. **Closure-behaviour measurement post-fix (NOT a formality).** The Phase 1.x/2 closures were calibrated against the *uniform advection that destroyed continents*. Some calibration values may be **compensatory artefacts** of that broken transport — e.g. Davis-Suppe may have been tuned aggressive *because* the convergent boundaries were the only place crust survived the evacuation; on correct transport that same tuning could **over-produce**. So the milestone must include an **explicit measurement of closure behaviour on the new transport**: *do the closures still produce what they were designed for (orogenic uplift magnitude, equilibrium cap, erosion incision, bathymetry anchors), or was their calibration an artefact?* — not merely "do cratons survive". Do NOT assume which; measure before scoping the follow-on. Same discipline as the whole investigation (counterfactual/measure-before-conclude).

   This measurement decides the follow-on scope:
   - **Closures hold** (calibration robust to the transport change) → light re-validation; incremental roadmap (piste 4 / Phase 3 resume quickly).
   - **Calibration was a transport artefact** (closures over/under-produce on correct transport) → heavier Phase 1.x/2 re-foundation.

## Conditional roadmap (sized by the post-fix closure measurement)

1. **Buoyancy milestone** (fix proven) + measure closure behaviour on the new transport.
2. **Per the measurement:** light re-validation OR closure re-foundation.
3. **THEN** piste 4 (rendering a correct field) + Phase 3 (morphology on continents that hold).

## Acceptance criteria

**Spatial bar (the judge — never a scalar/fraction alone):**
- Craton **area** preserved (target ≳ 90% of initial cratonic cells remain emergent).
- Compactness: land-mask perim/area low, largest-component high, n_components bounded.
- Position: craton centroid drift modest.
- Mass conserved (Δ ≈ 0%).
- Determinism test passes (bit-identical across two runs, same seed).
- Track-D closures re-validated; boundary accumulation routed through subduction/accretion (no spurious edge pile-up).
- Multi-seed (not seed 42 only).

## Why no further isolated prototyping

Diagnostic prototyping is **done**: cause proven (buoyancy, by exhaustion) and fix proven (rigid continent, 92%). The remaining boundary-accumulation question is a **design** matter (clean integration with the Track-D closures, their calibration, determinism), to be handled inside this milestone with full re-validation — not a standalone throwaway. Prototyping the caveat alone would be tuning before design.

## Labels

`C1`, `tectonics`, `milestone`, `re-validation`
