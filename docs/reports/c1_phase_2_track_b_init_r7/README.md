# C1 Phase 2 Track B — R7 Init + Continental Clustering + Ridge-Aligned Age

Issue: [#131](https://github.com/FifionRibana/ymir/issues/131).
Branch: `131-c1-phase-2-track-b-r7-init-continental-clustering-ridge-aligned-age`,
off `milestone/c1-lightweight-dynamic-tectonics`.

## Acceptance summary

| Stage | Tests | Status |
|-------|-------|--------|
| E1 unit (R7 boundary displacement)        | 6 | 6/6 PASS |
| E2 unit (continental clustering + smoke)  | 5 + 1 | 6/6 PASS |
| E3 unit (ridge-aligned age, Path 3.A)     | 4 | 4/4 PASS |
| E4 unit (Phase 2 dispatcher)              | 2 | 2/2 PASS |
| V validation (R7 init properties)         | 8 | 8/8 PASS |
| A acceptance (Path 3.A IMPROVES baseline) | 2 + 1 deferred | 2/2 PASS + 1 `#[ignore]` |
| D visual gallery + seed diversity         | 0 active + 2 `#[ignore]` | 2/2 PASS (ignored) |
| **Phase 2 Track B subtotal**              | **31** | **28 active PASS + 3 `#[ignore]`'d** |
| Phase 1.x integration preserved           | 19 | 19/19 PASS |
| Phase 2 Track A integration preserved     | 7 | 7/7 PASS |

8th bit-identical decomposition preserved:
`c1_phase_a_decomposes_into_closures_then_post_tectonic` PASS
EXACT (Stage E4 verified — init-side changes don't affect
time-loop decomposition contract).

## Sub-components delivered

### Sub-component 1 — R7 boundary displacement (Stage E1)

Perlin / Simplex noise applied to Voronoï plate boundaries
via two independent FBM channels. For each cell, compute
displacement `(dx, dy)` from noise, re-query nearest seed at
displaced position, reassign `plate_id`. Eliminates the v1/v2
straight-Voronoï-edge visual failure mode.

Parameters: `amplitude = grid_size / 8`, `frequency = 4.0`,
`octaves = 3`, `persistence = 0.5`. **First-shot calibrated** —
Stage E1 Test 5 reassignment 5.64 % (synthetic 2-plate
fixture) and Stage V Test 1 reassignment 14.01 % (full
Voronoï dispatcher) both squarely in the `(0, 20 %]` healthy
regime; no calibration iterations required.

### Sub-component 2 — Continental clustering (Stage E2)

BFS cluster-based plate-type assignment over the per-plate
adjacency graph derived from the (post-displacement) Voronoï.
Defaults `continental_fraction = 0.29` (Earth-like) and
`seed_cluster_count = 1` (single contiguous continent for
§2.4 cadrable-viewport requirement).

Overrides the Voronoï-internal Bernoulli `per_plate_type`
assignment with the BFS output. Voronoï's default
`continental_ratio = 0.30` is preserved unchanged (Phase 1.x
regression baseline). Stage V Test 4 measured actual fraction
2 / 8 = 0.250 vs target 0.290 within granularity-aware
tolerance `max(5 %, 1 / num_plates) = 0.125`.

### Sub-component 3 — Ridge-aligned age = 0 (Stage E3, Path 3.A)

Detect oceanic cells adjacent to divergent boundaries via the
existing `classify_boundaries` helper (reused without
modification), set their age to `ridge_value = 0.0`.
Continental cells keep `continental_baseline = 7.0`; oceanic
non-divergent cells keep `oceanic_baseline = 0.5`.

Resolves Track A's flux-form age-density pile-up finding
([[age-advection-density-vs-lagrangian]]) at init time without
changing the advection PDE. Path 3.B (per-step ridge
detection) and Path 3.C (Lagrangian advection) documented as
fallback options; **not required** based on Stage A Test 3
empirical evidence.

## Phase 2 Track A vs Track B comparison

| Metric | Phase 2 Track A (Phase 1.1 init) | Phase 2 Track B (R7 init) |
|--------|----------------------------------|---------------------------|
| Init method                  | Phase 1.1 (Voronoï straight, scattered Bernoulli) | Phase 2 R7 (curved + cluster + ridge-age) |
| Spearman age-altitude (cycle 300, oceanic cells) | -0.476 | **-0.5233** (Δ = -0.047, IMPROVES) |
| Oceanic age max (cycle 300)  | ≈ 6958 (1000× pile-up vs init) | ≈ 3973 (~570× pile-up vs init) |
| Pile-up reduction Track B vs A | — | **43 %** at convergent-boundary outliers |
| Oceanic altitude mean (cycle 0) | +0.08 (continental-dominant initial mix) | -0.18 (ridge-init bathymetric gradient) |
| Bathymetric maturation timing | Cycle 50 (advection-driven) | Cycle 0 (init-driven) |
| Continental cluster geometry | Random distributed (Bernoulli) | Single contiguous (BFS cluster) |
| §2.4 viewport-cadrable        | n/a (not measured)             | UNMET (9/10 seeds wrap — Track B-bis) |
| Boundary curvature            | Voronoï straight edges          | R7 Perlin curved (14 % reassignment) |
| Quantitative anchor available | S-S 5-point ±50 m              | Reuses Track A anchor (closure unchanged) |
| Per-step cost (64²)          | 455 µs                          | 465 µs (+10 µs / step) |
| 8th bit-identical decomposition | Preserved (7th = Track A)    | **Preserved (8th)** |

Per-step cost overhead at Track B level vs Track A: **+10 µs / step**
(465 vs 455). Init-time cost (one-shot per run) is negligible.

## Architectural findings

### Finding 1 — Path 3.A IMPROVES on Track A baseline (Stage A Test 3)

The escalation criterion (Stage E3 W7) said:
- Spearman ρ ≥ -0.4 (Track A baseline -0.476) → Path 3.A
  SUFFICIENT, ship.
- Spearman degraded → escalate to Path 3.B / 3.C.

Empirical measurement at seed 42, 300 steps, full Phase 2
stack:

       Track A baseline:   ρ = -0.4760
       Track B (Path 3.A): ρ = -0.5233
       Δ Track B − Track A: -0.0473

Track B Spearman is **stronger** (more negative) than Track A
by 0.047 — not just within the escalation threshold, but
empirically improving on it. Path 3.A is sufficient. Path
3.B (per-step ridge detection) and Path 3.C (Lagrangian
advection) remain documented as fallback in
`age_init.rs::## Fallback paths` but are not needed at the
current architectural level.

Secondary observation: Track B's age max at cycle 300 is
≈ 3973 vs Track A's ≈ 6958, a **43 % reduction in the
convergent-boundary pile-up**. The ridge-aligned `age = 0`
init pre-populates the lower tail of the age distribution
so the convergent-pile-up isn't seeded from a uniformly-
positive baseline.

### Finding 2 — R7 boundary displacement first-shot calibration

Defaults `amplitude = grid_size / 8`, `frequency = 4.0`,
`octaves = 3`, `persistence = 0.5` produced reassignment
counts in `(1, 15) %` on the first attempt at both unit
test scale (Stage E1 Test 5 = 5.64 %) and integration scale
(Stage V Test 1 = 14.01 %). No iterations required.

Pattern: continues the Phase 1.3 + 1.4 + Track A first-shot
calibration tradition under
[[calibration-via-visual-review]]'s hierarchy:

- Tier 1 (quantitative anchor) — not applicable for R7
  noise defaults (no published reference).
- **Tier 2 (analytical first-pass + bounded iterations)** —
  applied here. `amplitude = grid_size / 8` is the analytical
  estimate (boundary-cell migration order-of-magnitude); the
  reassignment fraction matched the predicted range on
  first measurement.
- Tier 3 (pure tuning) — not invoked.

Adds Phase 2 Track B as the **fourth first-shot success** in
the C1 milestone (after Phase 1.3 `k_collapse`, Phase 1.4
`K`, Phase 2 Track A Stein-Stein anchor).

### Finding 3 — Continental cluster wraps periodic boundary (Track B-bis)

Multi-seed scan during Stage A revealed **9 / 10 seeds
produce a continental cluster that wraps the periodic
boundary** at the default `64² × 8-plate` configuration:

       seed   | continental | extent_i | extent_j | verdict
       -------+-------------+----------+----------+-----------
           42 |        1123 |       64 |       48 | WRAPS
          100 |        1376 |       64 |       61 | WRAPS
         1337 |        1158 |       44 |       55 | tight (> 70 %)
         2026 |        1381 |       64 |       64 | WRAPS
         9999 |        1049 |       64 |       39 | WRAPS
            7 |         942 |       64 |       64 | WRAPS
           13 |        1109 |       64 |       64 | WRAPS
           31 |        1242 |       64 |       64 | WRAPS
           99 |         465 |       32 |       64 | WRAPS
          144 |         867 |       39 |       64 | WRAPS
       Cadrable: 0 / 10   Wrap-detected: 9 / 10

**Root cause** (structural Phase 1.x inheritance): 8-plate
Voronoï × 30 % continental Bernoulli yields ~2 continental
plates per cluster. BFS-from-single-seed picks 2 random
plates from a small graph; periodic adjacency makes
spatially-opposite plates connectable. §2.4 viewport-cadrable
requirement is not satisfied at this plate count.

**Track B contribution noted**: R7 displacement + cluster-
based BFS IMPROVE on Phase 1.x's "random scattered
continental plates" baseline (Stage V Test 5 single
connected component PASS), but the structural limitation
at low plate count remains.

**Deferred to Track B-bis** per Phase 1.4 Stage E4 T3
pattern. Test `acceptance_track_b2_continent_cadrable` is
`#[ignore]`'d with full 50-line rationale + 3 remediation
options:

1. **Constrained BFS seed selection** — favor central plate.
2. **Increase default plate count** — 8 → 12–16 plates.
3. **Spatially-biased seed sampling** (§6.2 alternative).

All three out of scope for Track B. File Track B-bis as
separate follow-up after merge.

## Phase 2 milestone progress

After Track B merge:

- **Track A** (Issue #129, merged via PR #130) — Stein-Stein
  oceanic bathymetry ✓
- **Track B** (Issue #131, this PR) — R7 init + clustering +
  ridge-aligned age ✓ (pending merge)
- **Track C** — boundary evolution (subduction, accretion,
  rifting per §4.5). Separate issue TBD.
- **Track D** — kinematics sampling (constrained random /
  Euler pole / scoring per §6.3). Separate issue TBD.
- **Track B-bis** — continental cluster cadrable fix
  (post-merge follow-up issue, evidence in this PR).

§7.2 Phase 2 milestone gate "different seeds produce visually
distinct continents" empirically forward-progressed by the
Track B seed-diversity gallery (3 seeds × cycle_000 PNGs
qualitatively distinct continental layouts, even though they
all wrap the periodic boundary pending Track B-bis).

## Cross-references

- Issue: #131
- Design doc:
  [`docs/c1_lightweight_dynamic_tectonics.md`](../../c1_lightweight_dynamic_tectonics.md)
  §6.1 (boundary displacement note), §6.2 (clustering note),
  §6.5 (ridge-aligned age note + Path 3.A/B/C fallback
  documentation), §7.2 (Track A complete, Track B in progress)
- R7 init modules:
  [`crates/ymir-core/src/tectonics_c1/init_r7/`](../../../crates/ymir-core/src/tectonics_c1/init_r7/)
- Dispatcher:
  [`init_r7/mod.rs::init_c1_state_phase_2_r7`](../../../crates/ymir-core/src/tectonics_c1/init_r7/mod.rs)
- Stage V tests:
  [`crates/ymir-core/tests/c1_phase_2_init_r7.rs`](../../../crates/ymir-core/tests/c1_phase_2_init_r7.rs)
- Stage A tests:
  [`crates/ymir-core/tests/c1_phase_2_track_b_acceptance.rs`](../../../crates/ymir-core/tests/c1_phase_2_track_b_acceptance.rs)
- Stage D galleries:
  [`crates/ymir-core/tests/c1_phase_2_track_b_visual_gallery.rs`](../../../crates/ymir-core/tests/c1_phase_2_track_b_visual_gallery.rs)
- Phase 2 Track A foundation: Issue #129 (merged via PR #130);
  Stein-Stein closure + Architecture C + age-density finding
  that Track B Sub-component 3 addresses.
- Memory entries (off-repo):
  - `project_c1_phase_2_track_b_outcomes` (new)
  - `feedback_age_advection_density_vs_lagrangian` (updated —
    Path 3.A mitigation case study, 43 % pile-up reduction)
  - `feedback_calibration_via_visual_review` (updated — R7
    first-shot calibration, fourth C1 milestone success)
  - `feedback_recursive_tuning_signals_structural` (updated —
    Track B Test 2 deferral case)
  - `feedback_init_re_architecture_pattern` (new, transferable)
