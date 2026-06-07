# Issue #145 — Stage 5: finger fix (≥2-neighbour subduction promotion)

The final-visual blocker: a 1px vertical FALSE-LAND finger on seed 2026's production binarymap (subduction grid-aligned promotion). Root pinned: **inter-step cascade** — a promoted cell becomes a continental neighbour for its oceanic neighbour next step, propagating a 1px line into the ocean (1 cell/step). The snapshot-per-step (intra-step) fix was measured INSUFFICIENT (finger 58→59, unchanged — inter-step dominates, as predicted by the 1-cell/step growth).

## Fix — promotion requires ≥2 convergent continental neighbours

A 1px finger tip has exactly ONE convergent continental neighbour (the cell behind it). A real accretion front has ≥2. So: gate the Oceanic→Continental promotion on **≥2 convergent continental neighbours**. Consumption (slab subduction) is unchanged (needs ≥1). Unconditional (not flag-gated) — the cascade is a pre-existing subduction defect the docstring flagged "acceptable for this iteration"; legacy benefits too.

## Validation (rigid, closures ON)

| seed | streaks before (no-flux) | streaks after (≥2) | continental |
|---|---|---|---|
| 2026 | 58 | **3** | 1400 (was 2018 cascade; accretion preserved) |
| 42 | — | 4 | 1135 |
| 1337 | — | 7 | 1158 |
| 2 | — | 5 | 1485 |

- **Finger eliminated** (streaks 58→3; binarymap visually clean on all seeds — user-confirmed credible continents multi-seed).
- **Accretion preserved** (continental still grows via broad fronts; 1400 between subduction-off 1381 and cascade 2018).
- **Regression**: only the subduction unit test `subduction_plate_id_reassignment_below_floor` changed (asserted the old 1-neighbour promotion) → rewritten as `subduction_single_neighbour_consumes_but_does_not_promote` + new positive `subduction_promotes_with_two_convergent_continental_neighbours`. Track D acceptance, mass conservation, boundary events all PASS.

## Status of the two boundary residuals

| residual | status |
|---|---|
| Finger (false land) | **FIXED** (≥2 rule) — was the blocker |
| Curtain (bounded mesh oscillation) | tolerated follow-up (rigid-boundary refinement) — cosmetic, capped by equilibrium, NOT false land |

The blocker (false land) is resolved on all seeds. Core #145 (continents = credible masses) holds. Remaining: 5b re-baseline (on the now-clean state), 5c determinism, 5d flip.
