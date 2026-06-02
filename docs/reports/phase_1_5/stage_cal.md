# Issue #141 Stage Cal — 6 first-shot calibrations re-validated under cap=0.92

Posture: **judge the PHYSICS under cap=0.92, not proximity to the old numbers.** The sea-level threshold moved (S̃-space 0.871 → ~0.39), so behavior *will* change — "it changed" is NOT a failure; "it produces something physically insensible" is. Re-tune only if physics is broken, documented. Seed 42, 64², workflow context (cap=0.92 / n_cycles=12).

## 2 independent calibrations — confirmed UNCHANGED by measurement

Both have **zero** `sea_level` / `compute_isostasy` / `iso_cfg` references in their code paths (equilibrium-height closure + init), and were re-run with their test `iso_config` flipped to `c1_default` — metrics **byte-identical** to the MinMaxFraction baseline:

| calibration | metric (default) | metric (c1_default) | verdict |
|---|---|---|---|
| Phase 1.3 k_collapse (equilibrium height) | global_max 1079.7/2.024/2.181; wedge_p95 0.376; fill_near 0.207 | **identical** | independent ✓ |
| Track B R7 (init-time age) | Spearman ρ = −0.5233 | **identical** | independent ✓ |

The indirect-coupling worry (equilibrium→S̃→threshold) does **not** bite: these calibration tests use `run_with_closures` with closures that either don't invoke the sea-level-dependent drainage, or measure quantities (equilibrium-capped global_max / wedge; advected age) that are unaffected by the threshold. Independence is real, not assumed.

## 4 dependent calibrations — measured under cap=0.92, all HOLD (no re-tune)

| calibration | measured behavior (cap=0.92) | physics verdict |
|---|---|---|
| **Phase 1.4 K erosion** | continental S̃ relief: min 0.238, mean 0.592, max 2.174 | **HOLDS** — the lower threshold drains/erodes a larger land area (~28% vs ~6%), yet the continent retains real relief (range 0.24–2.17, not razed flat). Erosion carves, doesn't obliterate. K unchanged. |
| **Track A S-S bathymetry** | S-S applies to ~72% oceanic (vs ~94%); depths down to −0.523 non-dim (≈ −2615 m); **coast coherence: 0/4096 mismatches** (altitude>0 ⟺ plate_type==Continental) | **HOLDS** — sensible bathymetry on the reduced ocean; the land/sea transition is exactly where plate_type flips. The Viz-0 Stage A plate_type-vs-altitude divergence does **not** recur (dual-space coherence, W2). depth_scale unchanged. |
| **Track D K_subduction** | cum_sub 24,738 cells over 240 steps | **HOLDS** — subduction fires substantially at convergent boundaries; sensible high-frequency regime under the new ~30/70 geometry. K_subduction unchanged. |
| **Track D K_rift (thinning)** | cum_thinning 4,998 cells | **HOLDS** — rifting thinning fires on continental divergent cells. K_rift unchanged. |

No calibration is physically broken under cap=0.92. None re-tuned.

## Finding (surfaced, NOT smoothed) — workflow suppresses Track-D CROSS-STEP events

`cum_merges = 0` and `cum_splits = 0` in the workflow. **This is NOT a cap effect and NOT a K re-tune.**

- **Mechanism (code-confirmed):** the accretion `ConvergenceTracker` and rifting `DivergenceTracker` are created **fresh at the top of every `run_with_closures` call** ([time_loop.rs:460-466](../../../crates/ymir-core/src/tectonics_c1/time_loop.rs#L460)). The workflow calls `run_with_closures(k_cycle=20)` once per cycle → the tracker resets each cycle → a maximum of 20 consecutive convergent/divergent steps, below the `merge_time_threshold = 50` → the cross-step events **never fire**. The per-STEP events (subduction, thinning) are unaffected and fire normally.
- **Cap-independent (measured):** the trackers use `plate_id` + `kinematics` velocity, never the sea level. The gallery path (`acceptance_full_run_seed_42_pangaea_collapse`, a single 300-step `run_with_closures` under the *same* `c1_default` cap after Stage E2) produces `cum_merges > 0`. So merges depend on **run length per `run_with_closures` call, not the cap** — a cadence × tracker-reset interaction.
- **Latent since #139** (the workflow has always used 20-step cycles; the #139 acceptance that asserts `merges>0` runs the *gallery* path). Surfaced now by the Stage Cal Track-D measurement.
- **Out of Phase 1.5 scope** (orthogonal to the sea-level cap). **Follow-up registered**: "workflow Track-D cross-step events" — thread the tracker across cycles, or lower the threshold for workflow cadence, or document that workflow mode shows subduction/thinning but not accretion-merges/rifting-splits (those need the gallery's long single run). To be scoped separately.

## Summary

- 4 sea-level-dependent calibrations: **all HOLD** under cap=0.92 (physics sensible; none re-tuned).
- 2 independent calibrations: **byte-identical** (confirmed by measurement, not assumption).
- Coast coherence (W2 / Stage A): **perfect** (0/4096) — surfaced early here.
- One cap-INDEPENDENT finding (workflow cross-step Track-D suppression) surfaced + registered as a follow-up, not conflated with the sea-level calibration.
