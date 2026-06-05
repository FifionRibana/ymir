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

1. **Rigidity mechanism** — how to mark continental/cratonic crust as non-advecting.
   - Binary (`v = 0` on continental cells) works (92%) and is mass-conservative (flux-form face cancellation; no freeze blow-up).
   - Graduated (resistance ∝ buoyancy/thickness) may be more physical — continental crust resists but is not infinitely rigid. Evaluate binary vs graduated.

2. **Boundary handling (the one caveat)** — oceanic crust converging on a rigid continental margin accumulates against it (cont/ocean edge mean 0.60→1.99, localized to convergent edges; p50 unchanged). This is physically subduction (ocean dives under the rigid continent) / accretion (ocean welds onto it). **Route the edge inflow through the existing Track-D subduction/accretion closures** so it is consumed/accreted rather than piled. Consuming the inflow removes the accumulation by construction. This is the integration point with Track D — not a tuning knob to prototype in isolation.

3. **Coherence with `plate_type`** — does rigidity follow `plate_type` (Continental/Cratonic)? `plate_type` evolves (accretion turns Oceanic→Continental; rifting splits). Rigidity must track that evolution: newly-accreted continental crust becomes rigid; rifted continental crust stays rigid. Define the coupling explicitly.

4. **Determinism** — preserve the C1 `Deterministic` invariant (same seed + config → bit-identical). The rigidity mask derives from existing deterministic fields (`plate_type`, `cratonic_mask`); ensure no new nondeterministic ordering.

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
