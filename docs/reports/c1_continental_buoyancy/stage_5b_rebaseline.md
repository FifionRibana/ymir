# Issue #145 — Stage 5b: imprint tests re-baselined on the clean state

Re-measured on the current committed code (post no-flux + finger-fix). The finger-fix does NOT change these tests (they run subduction-OFF), confirmed: values identical to the pre-finger-fix measurement. Migrated the 3 imprint tests to `rigid=true` (production transport; verdict A proved rigid is more faithful to Davis-Suppe than the legacy advection toe-pile).

## Per test

| test | change | new baseline + justification |
|---|---|---|
| 1.2 `davis_suppe_wedge_body_invariants` | config → rigid | PASSES unchanged — its invariants (bucket-count, global_max) aren't wedge_p95-valued. No assertion change. |
| 1.3 `davis_suppe_imprint_preserved_with_equilibrium` | rigid + re-baseline | wedge_p95 [0.3,0.5]→**[1.6,2.0]** (rigid wedge high, below h_eq, bit-identical to 1.2 rigid 1.7985); asymmetry>1.5 → **source taper** fill_near>fill_mid>fill_far (0.747>0.519>0.361) |
| 1.4 `erosion_preserves_davis_suppe_imprint_partially` | rigid + re-baseline + **un-ignore** | wedge_p95 [0.4,1.0]→**[1.3,2.1]** (≈1.64); asymmetry>1 → **source taper** (0.539>0.447>0.281); attribution documented |

## Key re-baseline principle (verdict A, stage_5b_asymmetry.md)

The legacy `asymmetry = near/far mean > 1` tested an **advection toe-pile** that INVERTED the Davis-Suppe critical taper. Under rigid transport the true source signature appears: **fill ratio decreasing with distance** (h_crit grows with distance). The assertion was changed from testing the ARTIFACT (toe-pile) to testing the PHYSICS (source critical taper) — a better assertion, not just a new value.

## Point 5 (curtain in wedge regions) — clear

These tests are subduction-OFF (no finger/stipple). The no-flux curtain (bounded oscillation) does NOT disrupt the wedge: the fill taper is cleanly monotone (0.747/0.519/0.361 and 0.539/0.447/0.281), which it wouldn't be if the curtain dominated the wedge metric.

## wedge_p95 attribution (1.4)

legacy 0.696 → erosion clean-removal (Step 1) 0.359 → rigidity 1.64. Two causes, opposite directions, documented in the test. Finger-fix: no effect (subduction-off).

Full suite green (only pre-existing `rectangular_simulation` v2 failure). Remaining: 5c determinism, 5d flip.
