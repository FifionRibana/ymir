# Issue #141 Phase 1.5 — acceptance checklist

Robust P95-capped sea level for C1 land/sea classification. cap=0.92 / n_cycles=12 (coupled calibration). Headless acceptance is automated; the manual part is the Viz-0.5 visual review.

## Headless acceptance (automated)

```bash
cargo test --release -p ymir-viz --features v2_legacy --bin ymir-viz bridge::c1::thread -- --nocapture
cargo test --release -p ymir-core --features v2_legacy   # v2 byte-identical guard
```

- **`acceptance_coast_coherence_phase_1_5`** (W2) — the S̃-space reclassify land set == the h-space `compute_isostasy` land set under `c1_default` (seed 42: 0.2830 vs 0.2832, within f32/f64 rounding). The Viz-0 Stage A plate_type-vs-altitude divergence does **not** recur because both instances use the same P95-cap mode.
- **`acceptance_emergent_land_multiseed_phase_1_5`** — emergent land distribution AROUND ~30% with natural per-seed variation: seeds {42, 1337, 2026, 7, 99} → 0.31 / 0.33 / 0.34 / 0.38 / 0.25, **mean 0.3213**, range 0.25–0.38 (each in [0.15, 0.45]; mean in [0.24, 0.38]; spread > 0.02 = not suspiciously uniform). Measured as the mean of the last 4 post-cycle fractions (band centre, robust to the ±0.05 limit-cycle stop-point variance).
- **`workflow_mode_continent_preserved`** (A1-c guard, reframed) — mass-conserving + bounded-band (late-cycle spread < 0.12) + emergent in [0.18, 0.45].
- **`workflow_converges_under_p95_cap`** (Q4 gate) — mass-conserving + bounded-band + no chronic per-step drainage jitter.
- **v2 byte-identical** — full ymir-core `--features v2_legacy` suite passes (MinMaxFraction default preserved); only the pre-existing #47 `deserialize_legacy_metadata_without_upscale` failure (unrelated).

## Manual visual review (Viz-0.5 workflow mode)

```bash
cargo run --release -p ymir-viz --features v2_legacy
```

C1 engine → Pipeline **Workflow** → Run (seed 42, 64², n_cycles defaults to 12).

- [ ] **Emergent land ~30%** — the continent occupies roughly a third of the map (not the ~5% of pre-Phase-1.5, not a near-full-land runaway). Multiple seeds show a *range* (~25–38%), not an identical figure.
- [ ] **Coast coherence (the W2 / Viz-0 Stage A check)** — hover a coast cell in the **Altitude** view: `plate_type` flips Continental↔Oceanic exactly where the non-dim altitude crosses 0. The plate_type contour (PlateType view) and the altitude=0 coastline (Altitude view) coincide — no divergence.
- [ ] **Hover sanity** — land cells show plausible non-dim altitude (> 0) and meters; ocean cells show negative altitude (Stein-Stein depth). Non-dim is the verification value; meters are the cosmetic lens.
- [ ] **Bounded coast fluctuation** — over the 12 cycles the coast breathes within a bounded band (the limit cycle), not a monotone collapse or runaway. This is accepted natural coast dynamics (Stage V finding).
- [ ] **Gallery mode unchanged** — switch Pipeline → Gallery: behaves as Viz-0.5 / Issue #137 (Gallery now also uses c1_default for its drainage + render, but has no reclassify, so its plate_type contour stays init-static).

## Known characteristics (documented, accepted)

- **Bounded limit cycle (±0.05), not a fixed point** — the per-cycle reclassify oscillates around the equilibrium. Accepted as natural coast dynamics (transgression/regression); the acceptance asserts a band, not a tight point. A reclassify-hysteresis damping is a Phase 1.5-bis follow-up (not required; product achieved).
- **Workflow plate topology is quasi-static** — accretion merges + rifting splits do NOT fire in workflow mode (the cross-step trackers reset per 20-step cycle, below the 50-step threshold; cap-independent). The Pangaea collapse (8→2) is **gallery-only**. **This is a registered Phase 3 prerequisite** (see `stage_cal.md`): settle whether Phase 3 needs an evolving topology in workflow mode before Phase 3.A.
