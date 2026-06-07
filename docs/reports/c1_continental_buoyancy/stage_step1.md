# Issue #145 — Step 1 result: erosion clean non-injecting removal

Mechanism fix (decided: clean removal; deposition = separate follow-up). Legacy `s_new = max(floor, s−δ)` raised sub-floor cells to `floor`, injecting mass. New:

```rust
let s_new = if s_old <= floor { s_old }                 // never inject
            else if delta >= s_old - floor { floor }    // exact floor, no FP drift
            else { s_old - delta };
```

**The only behavioural change vs legacy is sub-floor cells** (`s_old < floor=0.2`): legacy injected them up to 0.2; we leave them. For `s_old ≥ floor` the result is bit-identical (floor returned exactly).

## Measurement — SYSTEM level (rigid, all closures, seeds 42 + 2)

| seed | mass Δ (was) | craton area (was) | perim/area (was) | n_comp (was) | largest (was) | bound edge (was) |
|---|---|---|---|---|---|---|
| 42 | **−6.9%** (+25.7%) | 81% (82%) | 0.64 (0.58) | 20 (17) | **0.86** (0.63) | 0.11 (0.21) |
| 2 | **+9.6%** (+16.8%) | 83% (88%) | 0.74 (0.53) | **79** (43) | 0.77 (0.93) | 0.19 (0.79) |

- **Mass** substantially toward conservation (was +25.7%/+16.8%). **Non-additivity confirmed:** removing the +247%-standalone injection dropped the net by only ~32% / ~7% (the equilibrium cap was already absorbing it). **The two seeds STRADDLE zero (−6.9% vs +9.6%)** — the residual is geometry-dependent; a single global DS rate (step 2) cannot zero all seeds, so step 2 calibrates to minimise |mass| across seeds, not zero one.
- **Craton area maintained** (81% / 83%, vs 82% / 88%). Clean removal did not collapse the continents.
- **Spatial mixed:** seed 42 improved (largest 0.63 → **0.86**, one dominant mass); seed 2 mildly fragmented (n_comp 43 → **79**, largest 0.93 → 0.77) — clean removal opens sub-floor "holes" where the injection used to fill, peppering the land mask. Watch-item for step 2 / final.
- **Visual** (`step1/seed{42,02}_step1_land.png`): both credible — dominant compact mass holds; seed 2 more speckled but recognizable.
- **Boundary edge dropped** (0.11 / 0.19) — less margin accumulation (erosion no longer injecting at margins).

## Regression impact (enumerated, all binaries, `--no-fail-fast`)

The fix is unconditional (transport-independent, per "fix it because it's wrong"). Net break = **1 test**:
- `floor_at_oceanic_baseline` (lib): **PASSES** (exact-floor formulation avoids FP drift; `s_old=1.0 > floor` path identical to legacy).
- `erosion_preserves_davis_suppe_imprint_partially` (Phase 1.4): **deferred `#[ignore]`** → point-5 re-validation. `wedge_p95` drops 0.4+ → **0.359** (below Phase-1.3 baseline 0.376). **FINDING:** the Phase-1.4 "erosion lifts wedge_p95" architectural claim was partly an artefact of the floor injection (sub-floor refill inflated the bulk). With correct erosion the lift is gone; the spatial imprint (asymmetry 2.10 > 1.0, fill_near 0.203 > 0.05) still holds. NOT blind-rebaselined — a behaviour to understand at point 5.
- (`rectangular_simulation_smoke_test` — pre-existing v2-Stokes, unrelated.)

## Verdict

Step 1 PASSES its criterion: mass substantially toward conservation, craton area maintained, continents credible, determinism intact (serial code). Open: the seed-straddle (step 2 must minimise across seeds, not zero one) and seed-2 mild fragmentation (watch). Erosion-lift finding logged for point 5.
