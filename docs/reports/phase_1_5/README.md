# Issue #141 — Phase 1.5: robust P95-capped sea-level for C1 land/sea classification

The deepest core change since Track D. Replaces the C1 sea-level threshold's `min/max·fraction` formula — fragile to upper-tail outliers — with a percentile-capped formula, raising emergent land from ~5% to ~30% (the extractable-continent target) while keeping v2 byte-identical. Branch off `milestone/c1-lightweight-dynamic-tectonics`.

## The bottleneck (Issue #139 M2 evidence) — and the falsified prediction

Issue #139's Stage V diagnostic established (M2) that the low emergent land (~5%) was **not** macro erosion (mass-conserving) but the **adaptive `sea_level_ref = s_min + 0.4·(s_max − s_min)` drifting up** as `s_max` doubled (1.0→2.18 via Davis-Suppe orographic peaks): the threshold (0.871) sat above ~94% of the crust while the bulk was at P50≈0.28. On the *same* field, a P95-capped threshold gave ~28% land vs ~5.9% — **the formula, not the physics, gates emergent land.**

The Issue #141 ordering prediction was *"Phase 3 relief widens [s_min,s_max] → worsens the drift → Phase 3 would aggravate, so Phase 1.5 must precede Phase 3."* **A relief sweep (M1) FALSIFIED the mechanism**: Davis-Suppe `h_max` does not widen `s_max` (pinned ~2.18 by a soft cap) and *more* relief gives *more* land. The ordering conclusion (Phase 1.5 first) still held, but for the **measured** reason (the min/max formula is fragile to the *existing* ~1% upper tail; a correct sea level is a prerequisite to validate Phase 3 relief at all), not the predicted one. *Measure, don't assume — the prediction's mechanism was wrong; the data found the real cause.*

## The change

- **`SeaLevelMode { MinMaxFraction, PercentileCapped { cap_percentile } }`** on `IsostasyConfig` (`#[serde(default)]` → `MinMaxFraction`, so legacy configs and v2/export are byte-identical).
- **Both** sea-level instances branch on the mode — `compute_isostasy` h_sea (h-space) and `compute_sea_level_ref_s_space` (S̃-space) — coherent because `h = S̃·buoyancy` is monotonic, so both place the boundary at the same S̃ value (W2).
- **`IsostasyConfig::c1_default()`** = `PercentileCapped { 0.92 }`; the C1 engine uses it, v2 + export + gallery PNG generators keep `Default` (`MinMaxFraction`).
- Percentile via O(N) `select_nth_unstable_by` (deterministic).

## Calibration — cap=0.92 / n_cycles=12 (COUPLED), live multi-seed

M2's static post-hoc suggested 0.95→28%; **live, the in-loop feedback contradicted it** (0.95 oscillates persistently ~20%). A live multi-seed sweep settled the coupled calibration:

- **cap=0.92**: emergent distribution mean **30.6%**, range **24.5–36.6%** (natural per-seed variation). 0.95 oscillates+undershoots (~20%); ≤0.85 runs away (~95%).
- **n_cycles=12**: worst-case seed enters the equilibrium band by cycle 9 + margin (n_cycles=5 cuts mid-overshoot). Viz-side workflow default; PhaseAParams core default stays 5 (v2-shared).
- The two are **coupled** — cap=0.92 needs ~12 cycles to settle; documented as one calibration.

## Convergence is a BOUNDED LIMIT CYCLE, not a fixed point

The P95-cap system oscillates ±0.05 around its equilibrium (per-cycle reclassify reacting to the post-macro distribution). The convergence gate was **reframed** from a Δ-strict last-cycle test (which wrongly fails an oscillator — a 5-cycle "Δ=0.005 pass" was a fluke) to a **bounded-band** gate (late-cycle spread < 0.12) + mass-conserving + no chronic per-step jitter. *The criterion must match the system's nature — band for a limit cycle, Δ for a fixed point.* The ±0.05 fluctuation is accepted as natural coast dynamics (transgression/regression).

## Sub-fix NOT applied (measure-before-subfix)

The pre-spec'd Q4 per-cycle drainage sub-fix targets per-STEP drainage jitter — which Stage V measured as smooth (~0.01 sawtooth). The ±0.05 oscillation is per-CYCLE reclassify, a different mechanism the sub-fix would not damp. Not applied. A reclassify-hysteresis damping is a Phase 1.5-bis follow-up (not required; product achieved).

## Re-validation (Stage Cal — physics-first, anti-confirmation-bias)

Judged on the PHYSICS under cap=0.92, not proximity to old numbers.

- **4 dependent — per-step behavior HOLDS** (none re-tuned): Phase 1.4 erosion (continental relief 0.24–2.17, carved not razed); Track A S-S (sensible depths on the reduced ~72% ocean); Track D K_subduction (24,738) + K_rift thinning (4,998).
- **2 independent — byte-identical by MEASUREMENT** (not assumed): Phase 1.3 k_collapse (global_max/wedge identical) and Track B R7 (Spearman ρ −0.5233 identical). No hidden coupling, no init-path leak.
- See [stage_cal.md](stage_cal.md).

## Dual-space coast coherence (W2) — the Viz-0 Stage A divergence does NOT recur

Both the reclassify coast (S̃-space) and the altitude coast (h-space) use `c1_default`, so they agree: seed 42 reclassify land 0.2830 vs isostasy land 0.2832 (within f32/f64 rounding). The Viz-0 Stage A plate_type-vs-altitude divergence (caused by two *different* sea-level formulas) cannot recur. (NB: the structural `altitude>0 ⟺ plate_type` identity is tautological via Stein-Stein gating — the real test is the dual-instance land-set agreement.)

## v2 isolation (the hard guard)

`MinMaxFraction` is the `Default`; v2 + export + gallery PNG generators keep it. Full ymir-core `--features v2_legacy` suite: **586 pass** (only the pre-existing #47 `deserialize_legacy_metadata_without_upscale`, unrelated). The formula change is fully opt-in.

## 9th bit-identical — redefined (contract exact, values redefined)

The two decomposition tests (wrapper == `run_with_closures` + `apply_post_tectonic`) were switched to `c1_default` so they exercise P95-cap. The **contract holds byte-exact** under P95-cap (deterministic percentile → both paths identical); the S̃ reference values are redefined. The STOP condition (decomposition inconsistency) did not trigger.

## Stage S corrected the issue's spec (the audit's job)

The issue predicted production sites `time_loop.rs:811,839` + `phase_a_c1.rs:192`. The Stage S audit found **all three are test fixtures** — the real C1 sea level is **caller-threaded** (`C1TimeLoopConfig.iso_config` / `PhaseACycleInputC1.iso_config`), so the production wiring is in the **viz** (gallery worker, workflow worker, render altitude — the last coherence-critical). The spec was partially wrong; the audit caught it before implementation. This is exactly Stage S's role — recorded, not buried.

## Findings / follow-ups registered

1. **Workflow plate topology is quasi-static → Phase 3 PREREQUISITE.** Accretion merges + rifting splits do NOT fire in workflow mode: the cross-step trackers reset per 20-step `run_with_closures` cycle, below the 50-step threshold (cap-independent — the gallery path at the same cap fires merges). The Pangaea collapse (8→2) is **gallery-only**. Workflow is the Phase 3 production reference → settle before Phase 3.A whether Phase 3 needs an evolving topology (→ run-length / persistent-tracker fix becomes a Phase 3 prerequisite). See [stage_cal.md](stage_cal.md).
2. **Reclassify-hysteresis damping** (Phase 1.5-bis) — to tighten the ±0.05 limit cycle if a fixed point is ever wanted.
3. **Gallery PNG regeneration under (workflow + P95-cap)** — separate visual-reference maintenance; the committed Track A/B/D PNGs stay MinMaxFraction (standalone path) for now.

## Stages

| Stage | Deliverable | Commit |
|---|---|---|
| S | audit (corrected site list; caller-threaded) | `5cc80d2` |
| E1 | SeaLevelMode + dual-space branch | `e1b7c39` |
| E2 | wire 3 viz production sites to c1_default | `ed5d95f` |
| V | convergence + v2 guard; then cap=0.92/n_cycles=12 + bounded-band recalibration | `a2e2c79`, `65888f8` |
| Cal | 6 calibrations re-validated; cross-step → Phase 3 prerequisite | `5ba170a`, `7394a8c` |
| A | acceptance + coast coherence | `fcfbf1a` |
| Final | 9th-bit redefined; this README + memory + PR | `f60a39f`, … |
