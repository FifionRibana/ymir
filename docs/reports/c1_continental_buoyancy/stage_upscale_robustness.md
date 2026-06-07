# Upscale robustness to S̃ mesh non-convergence (#147 FOLLOWUPS-#6 gating)

**Question (the gate before opening the advection-scheme milestone):** the
S̃ field is mesh-non-convergent in production (r~0.51, #147). The reason to
care is the upscale (it adds the fine terrain). But does r~0.51 actually
perturb the upscale downstream — or is the upscale robust to it (as
Stein-Stein absorbed the mass swing, invisible in production)? **Measure
before opening a foundation milestone.**

**Method (throwaway C1→upscale):** run C1 production (rigid, full closures)
at 64² and 256², same seed 42; build the gallery altitude (isostasy +
Stein-Stein); normalise with the SAME fixed map (sea 0.0→0.5) so both are
comparable; `upscale_with_fbm` both to 1024² anisotropic FBM. Criterion =
**structure convergence, NOT identity** (FBM fills more from 256² by design;
"different paths" accepted). Harness: `c1_closure_morphology::upscale_robustness_64_vs_256`.

## Result

```
 coarse-altitude structure r (64 vs 256, block-meaned →64) = 0.8825
 UPSCALED structure r        (64 vs 256, block-meaned →64) = 0.9027
 upscaled land%  : 25.8 (from 64²)  vs  26.2 (from 256²)
 largest comp    : 0.951            vs  0.969
```

- **Upscaled structure r (0.90) ≈ coarse-altitude r (0.88)** — the upscale
  does NOT amplify the divergence; it PRESERVES the coarse structure
  convergence. The S̃ non-convergence (r~0.51) does NOT propagate.
- **Land morphology matches:** ~26% land both, one dominant landmass
  (largest ~0.96) in both → same continents.
- **Visual (the judge):** both renders show the SAME seed-42 world — same
  continent (upper-left mass + diagonal arm + left-edge island), same brown
  orogenic highland in the same place, same bathymetry/coasts. The 64²-based
  is blockier (coarse stairstep coast); the 256²-based smoother. **Same world,
  different fineness** — structure-convergence criterion MET.

## Why robust (mechanism)

`upscale_with_fbm` orients its anisotropic FBM by the slope of the **coarse
ALTITUDE** it is given (`compute_terrain_analysis(coarse)`), NOT the raw S̃
gradient. We feed it the production altitude (isostasy + Stein-Stein), which
converges at **alt r ~0.88** — the non-convergent S̃ (r~0.51) is "laundered"
through isostasy + Stein-Stein into a convergent altitude BEFORE the upscale
reads it. Same absorption mechanism as the DS mass swing absorbed by
Stein-Stein (invisible in production). The FBM then merely adds detail at the
coarse-cell scale (different fineness per base resolution — the accepted
"different paths"), averaging out under downsampling.

## VERDICT — advection-scheme milestone is DEFERRABLE

The S̃ mesh non-convergence is **invisible downstream**: the upscale produces
the same world at 64² and 256² base. C1 can proceed to upscale wiring /
Phase 3 / piste 4 **without** changing the advection scheme. The foundation
follow-up (FOLLOWUPS #6) stays REGISTERED but is **not urgent / deferrable** —
the gating measurement answered "robust", so the milestone is not opened.

**Caveat (the robustness precondition):** this holds because the upscale
consumes the **altitude** (laundered, convergent), not the raw S̃ gradient.
If a future wiring fed the raw S̃ gradient to the FBM orientation, robustness
would NOT hold (it would read r~0.51 directly). Keep the upscale on the
altitude path; revisit #6 only if that contract changes or a downstream
consumer reads raw S̃.
