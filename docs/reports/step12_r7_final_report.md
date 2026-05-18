# Step 12 — Final Report

**Branch:** `112-step-12-interleaved-tectonic-erosion-workflow`
**Sprint duration:** Step 12 R0 → R7 ω.3 (~6 weeks compute + diagnostic)
**Status:** Diagnostic complete, awaiting pivot direction before closure.

---

## TL;DR

Step 12 delivered the interleaved tectonic-erosion workflow (R0–R6 = infrastructure) and tested whether richer S̃ init profiles (R7.A = intra-plate topology) unlock Living Landz dynamics within the existing Phys.A model.

- R7.A.2.4 simulation **falsified** the composite-init hypothesis: Run C (Composite) ≈ Run A (Radial) on 4/5 auto criteria; Run B (Orogenic) produced impulsive cycle-1 dynamics that decayed ×60 over 5 cycles without sustaining.
- ω.3 diagnostic budget 1h15 identified the **structural limitation**: Phys.A is viscous-dissipation-dominated. v amplitude is bounded by η, not by mantle frequency or init mean gradient. Only init distributions with `frac>0.20 ≳ 1 %` (tail-heavy) excite a transient, which then dissipates without regeneration mechanism.
- **Conclusion:** the 2D thin-sheet + Voronoï + Stokes quasi-static + mantle Fourier ψ paradigm produces static + smoothed S̃ fields, not evolving Living Landz. The verdict is paradigmatic, not parametric. Pivot direction TBD with stakeholder.
- All Step 12 acquis (solver D2+D1-ter, macro_redistribution, Phase A workflow, V2 spec layer, test discipline, viz gallery) **transfer** to the next paradigm and should be preserved.

---

## Section 1 — Executive summary

### What Step 12 set out to do

The Step 12 issue framed the work as *Interleaved tectonic + erosion workflow*: replace the linear "tectonic finishes → erosion runs" pipeline with a multi-cycle loop where each cycle alternates between Stokes thin-sheet tectonics and a macro-redistribution pass (drainage + isostatic rebound + deposition). The acceptance contract — multi-dimensional criteria R4.1–R4.6 — required all axes (continents émergés, cratons préservés, bordures + chaînes evolved, mass conservation, drainage active, sustained dynamics) to pass simultaneously.

### What Step 12 delivered

- **R0–R6: infrastructure.** V2Field::Altitude rendering, drainage steepest-descent low-res, conservative macro-redistribution (mass drift < 1e-12 per cycle), Phase A loop integration, multi-dim acceptance scaffolding, solver fix (D2 portage from PR #49 + D1-ter post-macro reinit v=0), and mantle evolution_rate wiring (Phys.A phase-drift, latent bug pre-R6).
- **R7.A: empirical falsification.** R7.A.1 Orogenic-seul produced isolated spikes (insufficient visually); R7.A.2 Composite (dome + ridge with cap) was the alternative hypothesis but ended in FAIL on the R7.A.2.4 A/B/C sweep — Composite dynamics collapsed onto Radial dynamics on every auto metric.
- **R7 ω.3: structural diagnostic.** Four-axis investigation (mantle adiabaticity, Orogenic dynamics nature, Newton stall mechanism, ∇S̃ distribution) identified the viscous-dissipation-dominated regime as the structural cause. Documented as S1–S9.

### What Step 12 did NOT deliver

- Living Landz acceptance gate. R4.3 (visible chains + boundary deformation persistent across cycles) failed on every config tested.
- Phase B HD finalisation in the loop context (out of scope for this step; the Phase B path remains functional as a one-shot from Phase A final state, but no integration test exercises the cycle-driven feed).
- A "tuned" parametric configuration that produces Living Landz. The verdict from ω.3 is that no such tuning exists within Phys.A as currently formulated.

### The lecture-4 verdict

> *Step 12 closes "this approach" — not the project. The 2D thin-sheet + Voronoï + Stokes quasi-static + mantle Fourier ψ paradigm, as instantiated in Phys.A, structurally cannot produce evolving continental landscapes within the operational regime. The pivot is paradigmatic.*

This is the verdict reached by both empirical falsification (R7.A.2.4 FAIL on the last hypothesis available within the model) and theoretical analysis (ω.3 viscous-dissipation argument). It is not a tuning failure or an init choice failure; it is the consequence of the model's force-balance structure.

The user has signalled intent to pivot to a different resolution paradigm. This report does not specify that paradigm; Section 7 is a placeholder to be co-authored with the stakeholder once the direction is shared.

---

## Section 2 — Acquis Step 12 (R0–R6, technical infrastructure)

These are the deliverables of the interleaved workflow build-out, each of which works correctly in the regime where Phys.A operates and is reusable for any future paradigm that retains the underlying 2D field abstraction.

### 2.1 — R0 — V2Field::Altitude (post-isostasy heightmap rendering)

**Commit:** `182168d` — *FEAT : Step 12 R0 — V2Field::Altitude*

Added a new enum variant `V2Field::Altitude` in [`crates/ymir-viz/src/visualization/v2_viz.rs`](crates/ymir-viz/src/visualization/v2_viz.rs) that renders Airy-isostasy altitude derived from S̃ per-frame, with hypsometric colormap remapping so that sea level lands at the green-brown transition in `[0.4, 0.5]`.

Behaviour: `compute_isostasy(&Field2D::from_vec(s_field), IsostasyConfig::default()).heightmap` → piecewise-linear remap → hypsometric ramp. Caveat documented in R7.A.2.3: the underlying `compute_isostasy` uses per-image `h_min` / `h_max`, so `V2Field::Altitude` is **adaptive** by default — see [`feedback_viz_palette_absolute_for_comparison`](../../memory/feedback_viz_palette_absolute_for_comparison.md) for the inter-run pitfall this creates. Mitigation: a fixed-scale variant is used in R4 visual checkpoint and R7.A.2.4 sweep tests (see `save_altitude_fixed_scale` helper in [`crates/ymir-viz/tests/v2_workflow_r4_visual_checkpoint.rs`](../../crates/ymir-viz/tests/v2_workflow_r4_visual_checkpoint.rs)).

### 2.2 — R1 — Drainage steepest-descent low-res

**Commit:** `cd2eeee` — *FEAT : Step 12 R1 — drainage low-res utility*

Implemented a low-res steepest-descent drainage solver in [`crates/ymir-core/src/tectonics_v2/workflow/drainage.rs`](../../crates/ymir-core/src/tectonics_v2/workflow/drainage.rs). Each continental cell with positive altitude points to its lowest neighbour (8-connected, periodic boundaries), producing a target field used by the macro redistribution.

Properties:
- Deterministic order, no Rayon batching needed at this resolution
- Periodic-aware neighbour resolution (wraparound at all four edges)
- Returns `max_path_length` for R4.5 acceptance (drainage active)
- Sink-handling: cells without downhill neighbour stay in-place (later filled by isostatic rebound)

### 2.3 — R2 — macro_redistribution conservative

**Commit:** `dbd5902` — *FEAT : Step 12 R2 — macro mass redistribution (drainage + deposition + isostatic rebound)*

Implemented the macro-scale redistribution kernel in [`crates/ymir-core/src/tectonics_v2/workflow/macro_redistribution.rs`](../../crates/ymir-core/src/tectonics_v2/workflow/macro_redistribution.rs). Each call performs three actions in order:

1. Drainage: source cells contribute a fraction `α` of their excess above sea level into the drainage target.
2. Deposition: target cells receive the contribution with an isostatic rebound ratio applied.
3. Isostatic rebound: shallow cells subject to a global rebound towards equilibrium thickness.

**Conservation:** mass drift per cycle measured at `~ 1e-12` (numerical roundoff only) on every configuration tested across Step 12. This is the strongest invariant of the workflow path.

### 2.4 — R3 — Phase A integration (BREAKING beta drop)

**Commit:** `ea660c3` — *FEAT : Step 12 R3 — integrate macro_redistribution in phase_a, drop beta, add rebound + drainage_distance*

Replaced the legacy `beta` parameter (a hand-tuned reset coefficient) with the macro-redistribution pass invoked between cycles in [`crates/ymir-core/src/tectonics_v2/workflow/phase_a.rs`](../../crates/ymir-core/src/tectonics_v2/workflow/phase_a.rs). BREAKING change to the V2 spec layer: `V2PhaseAParams` lost `beta` and gained `alpha`, `isostatic_rebound_ratio`, `max_drainage_distance`. All presets migrated in the same commit.

The Phase A loop now alternates:

```
for cycle in 0..n_cycles:
    run_baseline(cfg, k_cycle_steps)        # tectonic phase
    s_isostatic = compute_isostasy(s)
    s_post = macro_redistribution(s_isostatic, plate_data, drainage, α, rebound, …)
    plate_type = reclassify_continental_oceanic(s_post, threshold)
    cratonic_factor = recompute_craton(s_post, age, plate_type)
```

### 2.5 — R4 — Multi-dimensional acceptance criteria

**No code commit** — methodology document and test scaffolding.

Replaced the "peak |v| > threshold" single-axis acceptance from earlier steps with the R4.1–R4.6 multi-dimensional gate:

| axis | criterion | source |
|------|-----------|--------|
| R4.1 | continents émergés (peak S̃ > sea level) | scalar metric |
| R4.2 | cratons préservés (retention ≥ 50 %) | per-cycle scalar |
| R4.3 | bordures + chaînes evolved | **visual** review of altitude_fixed PNGs |
| R4.4 | mass conservation (loss < 1 %/cycle) | scalar |
| R4.5 | drainage active (max_path ≥ 5) | scalar |
| R4.6 | dynamique soutenue (peak |v| > 0.1 in ≥ 3/5 cycles) | per-cycle scalar |

The acceptance discipline (`reject if any axis FAIL or unmeasured`) is captured in the [`feedback_multidim_checkpoint_metrics`](../../memory/feedback_multidim_checkpoint_metrics.md) memory and tracks the rule from earlier steps that "peak |v| alone is not validation".

R4 visual checkpoint deliverable: [`docs/reports/step12_r4_visual_checkpoint/`](step12_r4_visual_checkpoint/) — 2-preset × 6-state × 2-view gallery.

### 2.6 — R5b — Solver fix (D2 portage + D1-ter reinit)

The R5/R5b investigation surfaced two solver pathologies introduced in Step 0 (post-clean-rewrite) that hadn't been caught by the existing regression battery:

- **D2 — diagonal/off-diagonal scaling in CG preconditioner.** Ported from upstream PR #49 (a fix that had been merged on a sibling branch but not on the active branch). Restoration of the Jacobi preconditioner's coverage of the cratonic stiffening term; reduced CG iterations per Newton outer step by `~ 85 %` on workflow ON configurations.
- **D1-ter — Newton warm-start reinit `v = 0` post macro-redistribution pass.** Investigation found that warm-starting Newton from the *pre-macro* v field after the field S̃ had been redistributed produced post-macro convergence failures. Setting `v ← 0` before the post-macro Newton solve restored convergence at no performance cost (Newton finds the new equilibrium from rest faster than from a stale prediction).

**R5b mf sweep** confirmed the binary regime: `mf` ∈ {0.5, 0.7, 0.8, 1.0} → preserved / borderline / transition / dissolves. The transition value (~0.8) is sharp; this is the first instance of the same-pattern-across-axes signal (see Section 5.2).

Source: [`crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs`](../../crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs) (post-macro reinit), [`crates/ymir-core/src/tectonics_v2/solver/`](../../crates/ymir-core/src/tectonics_v2/solver/) (D2 portage).

### 2.7 — R6 — mantle.evolution_rate wired

The pre-R6 mantle field was **static after init** despite the `MantleConfig::Enabled.evolution_rate` field being present and serialised. R6 wired the field through the harness step loop and the stream-function builder ([`crates/ymir-core/src/tectonics_v2/mantle/stream_function.rs:135`](../../crates/ymir-core/src/tectonics_v2/mantle/stream_function.rs#L135)).

Formula: `ω = evolution_rate · TAU`; each Fourier mode's `(φx, φy)` shifts by `ω · t_nondim` at sampling time. Wave numbers, amplitudes, and the t=0 normalisation factor are frozen (no per-step renormalisation, which would jitter the argmax position).

**R6.3 sweep** at evo × mf concluded EVO.C (Section 3.2 below), the second same-pattern signal across orthogonal axes. The pre-R6 static-mantle bug is documented retroactively in Section 5.1 — Step 8's "mantle ON" regime was operating on a frozen pattern, which masked the adiabatic-following question for years.

---

## Section 3 — R7 falsification empirique

### 3.1 — R7.A.1 — Orogenic-seul

Implementation: [`crates/ymir-core/src/tectonics_v2/init/orogenic_profile.rs`](../../crates/ymir-core/src/tectonics_v2/init/orogenic_profile.rs).

A linear orogenic ridge per continental plate: PCA principal axis (with periodic-aware circular-mean centroid and min-image unwrap for plates that straddle the domain boundary), Gaussian transverse profile, smoothstep longitudinal modulation. Defaults: `peak_value = 1.20`, `base_continental = 0.85`, `half_length_ratio = 0.40`, `width_sigma_ratio = 0.08` (with a 0.10 variant tested in R7.A.2.4 to widen the chain).

**Visual result (R7.A.1.3 preview at 64²):** the ridge produced **isolated spikes** along the PCA axis rather than continuous chains. At small plates (L_plate ≈ 12 cells, σ ≈ 1 cell), the Gaussian transverse profile is essentially a 1-cell-wide ridge — visually a spike, not a chain. The longitudinal smoothstep produces a few discrete maxima rather than a connected ridge.

**Verdict R7.A.1:** insufficient visually for Living Landz "chains across continents". User proposed pivot to R7.A.2 (Composite: dome + ridge) over R7.B (Voronoï hierarchical) as the cheaper next experiment, with R7.B reserved as fallback.

### 3.2 — R7.A.2 — Composite (dome + ridge with cap)

Implementation: [`crates/ymir-core/src/tectonics_v2/init/composite_profile.rs`](../../crates/ymir-core/src/tectonics_v2/init/composite_profile.rs).

Additive formula with cap:
```
S̃ = clamp(dome + (peak - base) · ridge_amount, 0, cap)
```
where `dome` is a RadialProfile (smoothstep distance-to-boundary), `ridge_amount ∈ [0, 1]` is the Orogenic transverse × longitudinal modulation, and `cap = peak_orogenic` (configurable; the default `UsePeakOrogenic` ties cap to the ridge peak so no point exceeds the orogenic ceiling).

R7.A.2.1 spec doc: [`docs/reports/step12_r7_a_composite_profile/R7_A_2_1_formula_spec.md`](step12_r7_a_composite_profile/R7_A_2_1_formula_spec.md).

**R7.A.2.3 init preview** at 64² showed a visually coherent dome with a superposed ridge for all 8 plates. Numerical statistics confirmed peak S̃ around 1.18 (just under cap = 1.20) and frac>0.95 in the 5–7 % range for continental cells. The init *looked* like Living Landz.

**Then R7.A.2.4 measured the dynamic response.**

### 3.3 — R7.A.2.4 — Sweep A/B/C × 5 cycles × 20 steps

Test fixture: [`crates/ymir-viz/tests/v2_workflow_r4_visual_checkpoint.rs::r7_a_2_4_simulation_abc`](../../crates/ymir-viz/tests/v2_workflow_r4_visual_checkpoint.rs).

Configuration:
- 64² × 8 plates × continental_ratio 0.3 × seed 42
- mf = 1.0, evolution_rate = 0.10, cratonic_amp = 3 (the R6.3 "best available" config)
- workflow ON, n_cycles = 5, k_cycle = 20 → 100 tectonic steps total
- 3 runs: Radial / Orogenic σ=0.10 / Composite (default)

**Results** (data from [`docs/reports/step12_r7_a_2_4_simulation/`](step12_r7_a_2_4_simulation/)):

| Run | Init | Runtime | peak \|v\| c1→c5 | retention | mass loss/c | auto |
|-----|------|---------|-------------------|-----------|--------------|------|
| A Radial   | RadialProfile  | 29 min  | 1.85e-3 → 1.36e-3 | 78.4 % | 0.034 % | **4/5** |
| B Orogenic | Orogenic σ=0.10 | 402 min | 2.53e0 → 4.4e-2  | 25.7 % | 0.587 % | **4/5** |
| C Composite| Composite       | 32 min  | 1.90e-3 → 1.39e-3 | 80.5 % | 0.035 % | **4/5** |

Detailed per-cycle Newton/CG stats in each `<run>/metrics.md`.

**Visual verdict R4.3** (user inspection of `cycle_5_altitude_fixed.png` for each run, shared palette `[0, 1.5]`):
- Run A: smooth dome, no boundary deformation, no persistent chains
- Run B: residual peak at the orogenic spike location, surrounding field smoothed by erosion, no chain
- Run C: visually indistinguishable from Run A (the composite cap had killed the ridge gradient in the cap region; what survived was the dome)

**Numerical verdict:** Run C ≈ Run A on every auto metric. Composite dynamics are Radial dynamics — the addition of an orogenic component, capped at the orogenic peak, contributes nothing to the dynamical response.

### 3.4 — Verdict R7.A.2.4 = FAIL

Per the decision tree the user pre-committed: *"FAIL : Run C ≈ A → Composite n'apporte rien fondamentalement → Pivot R7.B prioritaire → Pas de R7.A.2.bis × N (V6 vigilance)"*.

The "pas de R7.A.2.bis × N" rule is enforced by [`feedback_recursive_tuning_signals_structural`](../../memory/feedback_recursive_tuning_signals_structural.md) — three tuning sweeps had already produced the same pattern (binary regime), the fourth would have been wasted budget. ω.3 was the structural-diagnostic alternative.

---

## Section 4 — ω.3 diagnostic structurel

### 4.1 — Méthode et budget

Budget cible: 1h30, hard ceiling 2h.
Budget consommé: ~1h15.

Four-axis investigation, no long simulation runs (all v field maps for B inferred from existing metrics.md tables; only D required new code, ~270 lines test, ~30s runtime). Test fixture: [`crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs`](../../crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs).

### 4.2 — A — Mantle adiabatic analysis (τ_ψ / τ_relax)

**Phase formula** ([`stream_function.rs:135`](../../crates/ymir-core/src/tectonics_v2/mantle/stream_function.rs#L135)): `ω = evolution_rate · TAU = 2π × 0.10 ≈ 0.628 rad / nondim-time`.

R7.A.2.4 timing: `T_total = 6.0`, `N_steps = 100`, `Δt = 0.06`, `T_cycle = 1.2`.

**Phase drifts:**
- per step: `ω · Δt = 0.038 rad ≈ 2°/step`
- per cycle: `ω · T_cycle = 0.754 rad ≈ 43°/cycle`
- total: `ω · T_total = 3.77 rad ≈ 216°` over 5 cycles (pattern passes through π → sign-flip)

**Characteristic times:**
- `τ_ψ` (half-period at fixed cell) = `π/ω = 5.0 nondim ≈ 83 steps`
- `τ_relax` (Stokes solve ≈ 1 Newton outer) = `1 step = 0.06 nondim`

**Ratio: τ_ψ / τ_relax ≈ 83** → strongly adiabatic. Stokes equilibrates to the mantle pattern instantly at each step; the pattern itself drifts by only 2° per step.

**But the operative bound on dynamics is *not* this ratio.** Advection check:
- Radial: peak |v| 1.4e-3 × T_total 6.0 = 8.4e-3 nondim = **0.54 cell displacement over 5 cycles**
- Composite: same scale → ~0.55 cell
- Orogenic c1: peak |v| 2.5 (impulsive); over T_cycle 1.2 = 3.0 nondim displacement = ~190 cells/cycle (saturating, geometry-dominated)

Conclusion: the system *is* adiabatic, but the cause is **viscous dissipation bounding v amplitude**, not the slowness of ψ rotation. Stokes is quasi-static by construction (`v ~ f · L² / η` instantaneously); increasing ω would only displace the location of the forcing, not its amplitude.

**Phys.A adiabatic regime is structural, not tunable by evolution_rate.**

### 4.3 — B — Orogenic dynamics = transient relaxation

(*Qualifier: not visually confirmed (v field maps absent from R7.A.2.4 capture set); inference from per-cycle metrics.*)

Run B per-cycle signature ([`run_b_orogenic_sigma_10/metrics.md`](step12_r7_a_2_4_simulation/run_b_orogenic_sigma_10/metrics.md)):

| metric            | c1     | c2    | c3    | c4    | c5    | pattern             |
|-------------------|--------|-------|-------|-------|-------|---------------------|
| peak \|v\|        | 2.53e0 | 8.3e-1| 1.2e-1| 9.8e-2| 4.4e-2| ×60 monotone decay  |
| Newton stall+cap  | 69/105 | 1/105 | 0/105 | 0/105 | 0/105 | transient c1 only   |
| frac>0.95         | 0.266  | 0.217 | 0.156 | 0.164 | 0.115 | ridge eroded        |
| mass              | 1987   | 1988  | 1973  | 1960  | 1951  | -1.8 % over 5c      |
| CG iters mean     | ~8700  | ~14600| ~15500| ~15700| ~14700| rise c1→c2, stable  |

**Reading:**
- ×60 monotone decay rules out **sustained convection** (which would maintain a steady-state v amplitude in equilibrium with mantle forcing). Decay tells us no external pump replenishes kinetic energy.
- Newton stall present in c1 only, never repeated → rules out **recurrent numerical instability** (which would oscillate or repeat).
- CG iter counts non-explosive (stable around 15k mid+late) → rules out **mode instability** (which would diverge).
- Mass loss 1.8 % is large but bounded, consistent with **erosion-deposition of an over-steep landscape** — physical, not artefact.

**Verdict B:** Orogenic cycle 1 is an **impulsive relaxation** of an init concentrated-gradient field. Cycle 1 dissipates the energy stored in the gradient concentration; cycles 2–5 then ride the same adiabatic curve as Radial (peak |v| around 1e-1 decaying to 1e-2). Not convection, not Rayleigh-Taylor, not numerical: an impulse response to a sharp init.

### 4.4 — C — Newton stall localised, Orogenic c1

Raw datum: Run B c1 Newton C/S/D/Cap = 36 / 66 / 0 / 3 = **65.7 % stall+cap** vs 100 % converged for A and C.

**Mechanism:** Newton-CG resolves the rheological non-linearity by linearising about the current iterate. The linearisation has a per-cell convergence radius set by the smoothness of the constitutive law. When the local gradient |∇S̃| exceeds that radius in a cell, the linearised step diverges or stagnates → marked "stalled" (reached `max_outer_iters` without converging under tolerance) or "capped" (a step-size clamp applied).

After cycle 1, the macro-redistribution + isostasy + reclassify passes smooth the field; cycle 2 starts with a redistributed gradient close to Radial-like distribution → Newton converges cleanly.

**Important:** this is **not** a solver bug. D2+D1-ter is robust: zero divergences across 100 substeps even in the worst run. It's a pathology of the **coupling Newton-rheology face to a concentrated init**, fully explained by the gradient distribution finding (D below). The 6.7h runtime of Run B is the operational cost of that pathology — every future init with a tail-heavy distribution will incur a similar cost on cycle 1.

### 4.5 — D — ∇S̃ distribution per init mode

Test: `r7_omega_3_d_gradient_diagnostic_three_modes` ([`crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs`](../../crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs)).

Sobel-periodic per-cell magnitude, 64² seed 42 num_plates 8 continental_ratio 0.3 amplitude 0.

| mode               | mean   | p90    | p99    | max    | frac>0.05 | frac>0.10 | **frac>0.20** |
|--------------------|--------|--------|--------|--------|-----------|-----------|---------------|
| radial             | 0.0356 | 0.114  | 0.158  | 0.204  | 32.2 %    | 16.0 %    | **0.22 %**    |
| orogenic σ=0.10    | 0.0357 | 0.115  | 0.363  | 0.363  | 14.1 %    | 13.1 %    | **9.72 %**    |
| composite          | 0.0367 | 0.115  | 0.162  | 0.242  | 33.2 %    | 16.6 %    | **0.29 %**    |

Shared palette: `[0, 0.36]`. Slope PNGs in [`step12_r7_omega3_gradient_diagnostic/{radial,orogenic_sigma_10,composite}/slope.png`](step12_r7_omega3_gradient_diagnostic/).

**Key finding:** mean and p90 are **identical** across the three modes; the distributions only differ in the tail. Composite shifts mean and max very slightly above Radial (additive ridge contribution survives somewhat) but **frac>0.20 stays at the Radial level** (0.29 % vs 0.22 %) — the cap killed the ridge gradient in the cap-active region where the dome was already at peak.

Orogenic's `frac>0.20 = 9.7 %` is **~40 × higher** than Radial/Composite — and that 9.7 % is what drives c1 dynamics + Newton stall.

**Visual confirmation** (user inspection, recorded in conversation):
- Radial: gradient diffus uniforme intra-plate
- Orogenic: noir intérieur + queue lourde aux contours plate + 4 spikes ponctuels au centre
- Composite: Radial + petits points centraux (residual ridge, mostly absorbed by the dome+cap)

Bimodal vs unimodal distribution character empirically and visually validated.

**Characterisation "dynamic" vs "adiabatic" init:**

| criterion       | "dynamic" (Orogenic)    | "adiabatic" (Radial/Composite) |
|-----------------|-------------------------|-------------------------------|
| mean ∇S̃        | indifferent             | indifferent                   |
| p90 ∇S̃         | indifferent             | indifferent                   |
| **frac>0.20**   | > 1 % suffices          | < 0.5 % traps quasi-static    |
| distribution    | bimodal (heavy tail)    | unimodal (light tail)         |

→ **What excites Phys.A is the tail of the distribution, not the bulk.** Composite killed the tail by capping.

### 4.6 — S1–S9 structural limitations

Synthesised from R0–R7 history and ω.3 findings.

#### S1 — Viscous-dissipation-dominated coupling (the principal lock)

`v ~ f_GPE · L² / η` imposes that v is bounded by η for a given f_GPE. At mf=1.0, constant η, and f_GPE controlled by |∇S̃|, the model produces `peak |v| > 1e-2` only if **multiple percent of cells have |∇S̃| > 0.20** (heavy tail).

Implication: smooth init (Radial, Composite-with-cap, Gaussian, Uniform) is structurally adiabatic. The dissipation absorbs the forcing before crust displacement.

#### S2 — Cap + composite annule the concentration

The `cap = peak_orogenic` mechanism destroys the tail where it forms (dome + ridge overlap region). For a Living Landz "dome with emergent chains", the cap would need to permit overshoot, the composition formula would need to be multiplicative or piecewise, or the topology would need an entirely different additive structure. No such fix was tested in R7 (and per the decision tree, would have been "R7.A.2.bis × N" — recursive-tuning territory).

#### S3 — Mantle adiabatic bounded by viscosity, not by frequency

`τ_ψ / τ_relax = 83` → adiabatic following confirmed. But increasing `evolution_rate` does not lift the v-amplitude lock: Stokes is quasi-static by construction, v is set by `f / η`, not by `df/dt`. Even an infinitely slow mantle would produce the same peak |v|. `evolution_rate` is **not a dynamics lever** in this regime; it is a spatial-redistribution lever for the forcing pattern.

#### S4 — Newton converges poorly on a concentrated init

The Newton-CG coupling stalls when ~10 % of cells exceed the linearisation convergence radius (Section 4.4). 6.7h runtime for Orogenic R7.A.2.4 Run B is the budgetary cost of that pathology. Implication: any "tail-heavy" init is expensive on cycle 1. The more dynamic the init, the more expensive the cycle 1.

#### S5 — Cycle 1 is unique relaxation, no sustained mechanism

Orogenic peak |v| decays ×60 across 5 cycles. Cycle 1 is an impulsive release; cycles 2–5 fall onto the Radial adiabatic curve. **Phys.A does not regenerate gradient in flight** — it consumes the tail of |∇S̃| without re-creating it. No active orogenesis, no subduction that lifts edges, no ridge push numerically effective. The model is **dissipative-net**: mantle redistributes spatially but does not create new discontinuities.

A Living Landz "evolving cycle after cycle" requires a gradient-creation mechanism, which is absent from the current model formulation.

#### S6 — Orogenic viscous cap reached in 1 cycle = no margin

Even Orogenic, the only init mode that excites c1, loses its dynamics in 2–3 cycles. There is no durable reservoir of tail. Increasing `cratonic_amp` or `mf` only stretches the relaxation timeline; it does not produce a limit cycle.

#### S7 — Binary regime reproduced across N parametric axes

Sweeps that have produced the same verdict:
- R5b: `mf` ∈ {0.5, 0.7, 0.8, 1.0} → preserved | borderline | transition | dissolves
- R6.2: `evo` ∈ {0, 0.10} → adiabatic in both
- R6.3: `cratonic_amp` ∈ {1, 2, 3} → EVO.C reframe (no config passes auto)
- R7.A.2.4: `init_mode` ∈ {Radial, Orogenic, Composite} → A = C (adiabatic) ≠ B (transient)

Four orthogonal axes, same conclusion: tight binary regime separated by structural thresholds. This is the meta-pattern documented in [`feedback_paradigm_limit_when_N_sweeps_reproduce_pattern`](../../memory/feedback_paradigm_limit_when_N_sweeps_reproduce_pattern.md).

#### S8 — 2D Voronoï geometry ≠ rigid plates

The model treats the lithosphere as a continuous thin-sheet medium, not as an assemblage of rigid floating plates. Consequences:
- No notion of "plate boundary moving as a block"
- Voronoï seeds serve only for **initial classification** (S̃ + craton + plate-type stamp)
- Once running, S̃ evolves cell-by-cell with no internal rigidity
- → **Continental collision** in the geological sense (two rigid blocks meeting) is **inaccessible by construction**

For Living Landz "with collisions", either the model needs to be extended (effective internal rigidity, plate-specific kinematics) or paradigm-shifted (DEM particle plates, distinct elastic substrates).

#### S9 — Robust acquis to preserve

Not everything is a limitation. The reusable Step 12 assets are detailed in Section 6.

---

## Section 5 — Findings rétroactifs et patterns méthodologiques

### 5.1 — Retrospective re-reading of earlier steps

The Step 12 investigation has shed light on three retrospective issues that were either latent or undiagnosed in earlier steps. These should be folded into the project's mental model when reviewing past results.

**Step 0 — clean rewrite lost solver robustness.**
The Step 0 reformulation of the thin-sheet solver replaced an existing implementation with a fresh codebase. R5b's investigation found that the Step 0 version was missing the D2 preconditioner scaling that PR #49 had introduced upstream. Implication: any Step 0–Step 4 benchmark that ran with this gap was operating with a less robust CG; the conclusions about regime boundaries from those steps may have under-estimated the model's actual operating margin.

**Step 8 — mantle calibration ratée for workflow ON regime.**
Step 8's slab + mantle co-calibration (§4.8 and §4.9 bands) was performed in the workflow OFF regime. When workflow ON activated post-Step 12 R3, those bands produced G > 1 runaway as documented in [`project_slab_mantle_cocalibration`](../../memory/project_slab_mantle_cocalibration.md). The lesson: any per-regime calibration must explicitly assert which regime it was performed in, and must be retested when the regime changes.

**Step 8 — evolution_rate latent bug.**
Before R6, `MantleConfig::Enabled.evolution_rate` was a serialised configuration field that **did not flow through to the harness step loop**. The mantle ψ was static after init regardless of the value. This means every "mantle ON" run prior to R6 was operating on a frozen pattern. Step 8's mantle ON regime, in particular, was de facto a "static mantle perturbation", not a time-evolving forcing. The conclusions about mantle's role in dynamics in pre-R6 steps must be qualified accordingly.

**Step 11 — runtime probably victim of the same latent bug.**
Step 11's slow runtime on the mantle-on workflow path may have been driven by a static-mantle pathology that, while not preventing convergence, set up a regime where the Newton solver was repeatedly converging to the same equilibrium under a frozen forcing. The runtime improvement noted in R5b (D2 portage + D1-ter reinit) is partly inherited by Step 11 as well; a re-run of Step 11 with R5b solver would likely complete substantially faster.

### 5.2 — Methodological patterns surfaced

Step 12 has surfaced or sharpened several methodological patterns. The memories listed below have been written or revised in this sprint and form the project's accumulated discipline.

- **[`feedback_recursive_tuning_signals_structural`](../../memory/feedback_recursive_tuning_signals_structural.md)** — when 3+ sweeps each end with "needs another knob", reframe. R5b → R6.3 → R7.A.2.4 each carried this signal forward.
- **[`feedback_paradigm_limit_when_N_sweeps_reproduce_pattern`](../../memory/feedback_paradigm_limit_when_N_sweeps_reproduce_pattern.md)** — sharper rule (this sprint): when N≥3 **orthogonal** axes reproduce the **same pattern**, the model is the limit. R6.3 was the canonical reframe-here point; momentum carried R7 forward despite the signal.
- **[`feedback_smoke_before_long_sweep`](../../memory/feedback_smoke_before_long_sweep.md)** — for sweeps > 1h budget, run one representative config first. R7.A.2.4 Run B's 6.7h could have been signalled by a 5-step smoke before committing the full 5 cycles.
- **[`feedback_multidim_checkpoint_metrics`](../../memory/feedback_multidim_checkpoint_metrics.md)** — list ALL axes (preservation, dynamics, conservation, convergence, visual) before running; reject if any axis FAIL or unmeasured. R4 acceptance discipline.
- **[`feedback_viz_palette_absolute_for_comparison`](../../memory/feedback_viz_palette_absolute_for_comparison.md)** — inter-run side-by-side renders need explicit `[vmin, vmax]`. R7.A.2.3 init preview surfaced this when the user spotted that adaptive palettes were inverting the apparent ordering of Composite vs Radial brightness.
- **[`feedback_init_distribution_tail_drives_dynamics`](../../memory/feedback_init_distribution_tail_drives_dynamics.md)** — tail-aware diagnostic (frac>0.20 ≳ 1 %) is the predictor of Phys.A excitability. Written this sprint based on ω.3 D.
- **[`project_phys_a_viscous_dissipation_dominated`](../../memory/project_phys_a_viscous_dissipation_dominated.md)** — the structural project-specific verdict.
- Warm-start anti-pattern (R5b D1-ter): warm-starting Newton from a stale v field after the underlying scalar field changed is anti-useful; `v ← 0` is faster. Not a general rule, but worth noting for any future iterative-solver coupling.

---

## Section 6 — Acquis réutilisables pour pivot

The pivot will retire the Phys.A formulation but should preserve as much of the surrounding infrastructure as possible. The following components are paradigm-independent (or with minor adaptation) and represent substantial sunk cost worth carrying forward.

### 6.1 — Solver D2+D1-ter

The Newton-CG outer/inner solver at [`crates/ymir-core/src/tectonics_v2/solver/`](../../crates/ymir-core/src/tectonics_v2/solver/) is robust in nominal regime (0 divergences across 100 substeps on every R7 run). The D2 Jacobi preconditioner scaling and the D1-ter post-macro reinit (`v ← 0` before each post-redistribution Newton solve) are documented fixes ported from upstream and bench-validated. For any future paradigm that retains a Stokes-like inner solve, this code is a starting point.

### 6.2 — macro_redistribution

[`crates/ymir-core/src/tectonics_v2/workflow/macro_redistribution.rs`](../../crates/ymir-core/src/tectonics_v2/workflow/macro_redistribution.rs) implements a conservative drainage + isostatic rebound + deposition pass with mass drift `~ 1e-12` per cycle. The algorithm is decoupled from the Phys.A force balance: it only sees S̃ and plate_type. Any future paradigm that produces an S̃-equivalent thickness field can plug into this pipeline directly.

### 6.3 — Workflow Phase A loop

The cycle structure (`tectonic → isostasy → macro → reclassify → craton`) at [`crates/ymir-core/src/tectonics_v2/workflow/phase_a.rs`](../../crates/ymir-core/src/tectonics_v2/workflow/phase_a.rs) is paradigm-agnostic: the "tectonic" step is the only Phys.A-specific component. Replacing it with a different dynamics model leaves the rest of the cycle intact.

### 6.4 — V2 spec layer

[`crates/ymir-viz/src/bridge/v2/`](../../crates/ymir-viz/src/bridge/v2/) provides JSON round-trip for the entire run configuration ([`spec.rs`](../../crates/ymir-viz/src/bridge/v2/spec.rs)), pre-built presets ([`presets/v2/*.json`](../../crates/ymir-viz/presets/v2/)), and a bridge-side `V2FinalState` thread-safe payload ([`events.rs`](../../crates/ymir-viz/src/bridge/v2/events.rs)). Any new dynamics paradigm should keep this surface — it lets the same UI and tests drive comparison between old and new paradigms.

### 6.5 — Test regression discipline

The `Disabled` regression battery (bit-identical Steps 0–10 baselines via opt-in `InitMode::Checkerboard`) and the "no production change without baseline test" rule have caught multiple regressions across the sprint. The pattern at [`crates/ymir-core/tests/`](../../crates/ymir-core/tests/) and [`crates/ymir-viz/tests/`](../../crates/ymir-viz/tests/) is reusable for any future paradigm: keep a frozen-config regression alongside the new development.

### 6.6 — Visualisation gallery

[`crates/ymir-viz/src/visualization/v2_viz.rs`](../../crates/ymir-viz/src/visualization/v2_viz.rs) provides V2Field rendering (S̃, Altitude, Age, Cratonic, StrainRate, VelocityMagnitude, Slope) with hypsometric / log-hot / Sobel-periodic colormaps. The slope rendering at [`v2_viz.rs:751`](../../crates/ymir-viz/src/visualization/v2_viz.rs#L751) was used unchanged for the ω.3 D diagnostic. Adding a new V2Field variant is mechanical.

### 6.7 — Diagnostic tooling

The init-only diagnostic pattern (build state from `init_s_field`, dump stats, render PNGs with shared palette) is captured in [`crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs`](../../crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs). Any future "characterise this init / parameter regime cheaply before committing to a sweep" task can clone this pattern.

---

## Section 7 — Direction pivot (TBD)

*This section is a placeholder. The pivot direction will be co-authored with the stakeholder once the desired resolution paradigm is shared. Candidate directions noted in the conversation include but are not limited to:*

- *Different dynamics formulation (e.g. lower-viscosity continuum, viscoelastic with stress memory, kinematic plate-velocity boundary conditions overriding the GPE-driven Stokes)*
- *Different geometry / topology (e.g. DEM particle plates with rigid-body internal dynamics, hierarchical Voronoï with sub-plate features)*
- *Different solver paradigm (e.g. cellular automata, agent-based plate motion, optimisation-driven landscape evolution)*
- *Different acceptance contract (e.g. statistical Living Landz instead of dynamic Living Landz — accept that the model produces landscapes via post-init redistribution rather than via in-flight tectonic evolution)*

When the direction is chosen, this section will document:
- The new paradigm's force-balance or motion-balance equations
- The Phys.A components retired and the new components introduced
- Migration plan for the Section 6 reusable assets
- New acceptance contract for the equivalent of R4.1–R4.6
- First-iteration milestone definition

---

## Appendix A — Output artefacts produced this sprint

| path | purpose |
|------|---------|
| [`docs/reports/step12_r4_visual_checkpoint/`](step12_r4_visual_checkpoint/) | R4 visual checkpoint gallery (2 presets × 6 states × 2 views + metrics) |
| [`docs/reports/step12_r4b5_mantle_sweep/`](step12_r4b5_mantle_sweep/) | R5b mantle / mf sweep results |
| [`docs/reports/step12_r4b_diagnostic/`](step12_r4b_diagnostic/) | R5b solver diagnostic (D1-ter investigation) |
| [`docs/reports/step12_r7_a_composite_profile/`](step12_r7_a_composite_profile/) | R7.A.2.1 formula spec, R7.A.2.3 init preview stats |
| [`docs/reports/step12_r7_a_2_4_simulation/`](step12_r7_a_2_4_simulation/) | R7.A.2.4 sweep — 3 runs × 5 cycles × {S̃, altitude_fixed} PNGs + per-run metrics |
| [`docs/reports/step12_r7_omega3_gradient_diagnostic/`](step12_r7_omega3_gradient_diagnostic/) | ω.3 D gradient distribution per init mode — shared-palette slope.png + per-mode stats |

## Appendix B — Memory entries written or revised

| memory | type | scope |
|--------|------|-------|
| [`feedback_viz_palette_absolute_for_comparison.md`](../../memory/feedback_viz_palette_absolute_for_comparison.md) | feedback | viz palette discipline for inter-run comparison |
| [`feedback_init_distribution_tail_drives_dynamics.md`](../../memory/feedback_init_distribution_tail_drives_dynamics.md) | feedback | tail-aware gradient diagnostic for viscous-dominated solvers |
| [`feedback_paradigm_limit_when_N_sweeps_reproduce_pattern.md`](../../memory/feedback_paradigm_limit_when_N_sweeps_reproduce_pattern.md) | feedback | meta-pattern: N-axis same-result → paradigm-level reframe |
| [`project_phys_a_viscous_dissipation_dominated.md`](../../memory/project_phys_a_viscous_dissipation_dominated.md) | project | structural verdict — Phys.A is η-bounded, paradigm-level limitation |

## Appendix C — Test fixtures relevant to Step 12

- [`crates/ymir-viz/tests/v2_workflow_r4_visual_checkpoint.rs`](../../crates/ymir-viz/tests/v2_workflow_r4_visual_checkpoint.rs) — R4 gallery + R5b sweeps + R6.3 evo×mf sweep + R7.A.2.4 simulation A/B/C
- [`crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs`](../../crates/ymir-viz/tests/v2_r7_omega_3_gradient_diagnostic.rs) — ω.3 D init-only ∇S̃ diagnostic (~270 lines, 0 production touch, <1s runtime post-build)
- [`crates/ymir-core/src/tectonics_v2/init/composite_profile.rs`](../../crates/ymir-core/src/tectonics_v2/init/composite_profile.rs) — R7.A.2 composite init mode + 6 unit tests
- [`crates/ymir-core/src/tectonics_v2/init/orogenic_profile.rs`](../../crates/ymir-core/src/tectonics_v2/init/orogenic_profile.rs) — R7.A.1 orogenic init mode + PCA + unit tests

## Appendix D — Commits (chronological, branch `112-step-12-interleaved-tectonic-erosion-workflow`)

Pre-R7 commits (R0–R6 + R5b solver fixes + R6 evolution_rate wiring) are listed in `git log` on the branch and indexed by the corresponding section above. R7 + ω.3 commits land in a single bundle alongside this report.

---

*End of Step 12 — Final Report.*
