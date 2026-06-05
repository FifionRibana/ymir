# Issue #145 — Step 2 verdict: SKIP DS recalibration (mass swing invisible in production)

After step 1 (erosion clean removal), the system mass residual is a **wide, geometry-dependent swing** across seeds — but it is **invisible in the production render**. Step 2 (Davis-Suppe recalibration) is therefore **unnecessary**.

## The chain of measurements

1. **6-seed distribution** (rigid, closures ON, DS ×1.0): mass Δ ranges −34% (seed 7/99/1337) to +26% (seed 2026), mean −11.8%, **spread 59.6 pts**. The initial 2 seeds (−6.9 / +9.6) were the tightest — misleading. The residual is **seed-geometry-dependent** → no global DS rate can tighten a 60-pt spread. BUT continents HOLD across all 6 (craton 76–83% where defined, largest 0.55–0.94, visual credible incl. the −33% seed 1337).

2. **Mechanism pinned** (seed 1337, ablation): the deep oceans (ocean S̃ ≈ 0.10, 85% below baseline 0.2) persist WITHOUT erosion (0.112) → **not erosion**. They are **advection-divergence** (oceanic crust advects to convergence, divergent zones empty toward 0 with no ridge crust-creation to refill — conservative). The −33% deletion is **equilibrium-height** capping the convergent piles (removing it → +30.5%). **Deposition (an erosion feature) would NOT fix either** — pinning saved scoping the wrong fix.

3. **Production render decides** (isostasy + Stein-Stein, the real `derive_altitude_field` path — prior renders used `compute_isostasy(S̃)` ONLY = raw field, NOT production): oceanic production altitude is **coherent across the swing**:

   | seed | S̃ mass Δ | RAW ocean alt (p50) | PRODUCTION ocean alt mean/p50/p90 |
   |---|---|---|---|
   | 1337 | −33% | 0.014 | −0.537 / −0.524 / −0.520 |
   | 2026 | +26% | 0.049 | −0.533 / −0.521 / −0.520 |

   Raw oceanic altitude differed ~3.5× (0.014 vs 0.049); **production oceanic altitude is identical within ~0.004**. **Stein-Stein re-anchors oceanic depth from AGE → the S̃ mass swing is absorbed, invisible in production.** Same lesson as the plate_type/altitude honest-render episode: do not judge on the RAW field when a downstream layer (Stein-Stein) re-interprets it.

## Verdict

- **Step 2 (DS recalibration): SKIPPED.** The mass swing it would target is invisible in production; calibrating DS would also worsen the negative-mass seeds. Continents hold.
- **Step 3 (equilibrium-height): confirmed as the load-bearing regulator** (point 6 + the pin: removing it → +30%). It is doing its job (capping convergent over-pile). DO NOT touch.
- The re-foundation reduces to **step 1 alone** (erosion clean removal), measured-justified.

## Registered follow-ups (NOT #145)

1. **Oceanic ridge crust-creation** (seafloor spreading at divergent boundaries) — the real cause of the raw-field deep oceans; a **fidelity** follow-up (invisible in production via Stein-Stein), not a necessary fix.
2. **Conservative erosion deposition** (sediment routing → plains) — richness (Living Landz city terrain), builds on the clean base.
3. **wedge_p95 re-validation** (point 5) — the `#[ignore]`'d Phase-1.4 test; the lift was floor-injection-dependent (now 0.359), understand don't blind-rebaseline.
4. **Init note:** seeds 99 & 1337 have ZERO cratonic cells (R7 clustering is stochastic); continents still form (continental clustering), just no cratonic-core label. The craton-area metric is undefined for those — use continental area.

## Next

Close #145 on the continental goal → **point 5 (re-validation matrix + flip the transitional flag default)**.
