# H1/H2/H3 discriminator — why S̃ morphology expresses faintly in altitude (#145 follow-up)

The closure-morphology stage found DS sculpts credible orogens in S̃ but
their imprint on rendered *altitude* at 64² is marginal. Three competing
hypotheses for the faint expression, and the measure that decides between
them (run BEFORE piste 4, per the user — piste 4 targets the render and we
must know whether the problem is the render (H1), the resolution (H2), or
whether C1's render is *meant* to be coarse (H3)):

- **H1 EXPRESSION** — S̃ rich, S̃→altitude conversion poor. → piste 4 fixes
  the render; 64² suffices.
- **H2 RESOLUTION** — 64² (~4096 cells) too coarse to carry detail. →
  raise C1 resolution.
- **H3 ARCHITECTURE** — 64² coarse BY DESIGN; detail is the upscale-to-4096²
  stage's job. → piste 4 renders the coarse structure; the real next
  chantier is the upscale pipeline.

## Measure 1 — resolution sweep (64² / 128² / 256², seed 42, full rigid stack)

Physical time held constant (`n_steps = 300·grid/64`, since CFL `dt ∝ 1/grid`).
Relief amplitude restricted to continental cells (where DS thickening lives):

```
 grid  n_steps  land% |   S̃ std  S̃ p95-5 | alt std alt p95-5 | alt/S̃
  64²     300    27.7 |  0.3097   0.8782  | 0.1078   0.3413   | 0.348
 128²     600    27.0 |  0.3061   0.7921  | 0.1208   0.3560   | 0.394
 256²    1200    28.1 |  0.3330   0.9378  | 0.1245   0.3428   | 0.374
```

**Amplitude is FLAT across resolution.** S̃ relief ~0.31–0.33 std, altitude
relief ~0.11–0.12 std, conversion ratio steady ~0.37. Visual: at 256² the
S̃/altitude fields are the SAME large-scale structure as 64² with SMOOTHER
boundaries — **no finer orogen chains, no richer internal relief**; the
continental interior stays a uniform plateau at every resolution. (Bonus
finding: the bounded curtain oscillation is MORE prominent at 256² —
grid-aligned blue striping radiating from boundaries scales up with res.)

## Measure 2 — architecture (docs + code, see Explore report)

- The upscale stage EXISTS (`terrain/upscale.rs`, `upscale_with_fbm`):
  bilinear interp of the coarse heightmap + anisotropic FBM, 7 octaves, to
  a configurable 1024²–8192² target. It is wired in the v2 workflow.
- It is **downstream by design**: `docs/tdd.md` §3 pipeline = Tectonics
  (128²–512²) → Isostasy → **Upscale + anisotropic FBM (4096²–8192²)** →
  HD erosion → climate → export.
- TDD §6.2 (decisive): the anisotropic FBM **reads the gradient of the
  tectonic thickness field (S̃)** to orient detail ("ridges run parallel
  to the collision front"). So the upscale consumes **S̃ gradient**, not the
  coarse altitude.
- `docs/c1_lightweight_dynamic_tectonics.md`: C1 scope = "Phase 1
  (tectonics) only. Phases 2–6 (isostasy, upscale+FBM, hydraulic erosion,
  …) are preserved and consume C1's output." C1 target resolution 64²–512²
  for runtime (<10 s). The C1→Phase-B(upscale) path is **implemented but
  NOT yet wired into the C1 viz bridge** (a roadmap item).

## Verdict

- **H2 RESOLUTION — REJECTED.** 64²→256² does not produce finer chains or
  richer altitude relief; amplitude is flat, structure merely smoother, and
  the curtain artefact worsens. Detail does not live at the tectonic
  resolution at any tested resolution. Raising C1 resolution is NOT the
  lever (and would re-open the 64²-calibrated cap/n_cycles/closures).
- **H1 EXPRESSION — real but NOT the lever.** Isostasy passes ~37 % of S̃
  relief to altitude (steady, linear-ish, not catastrophic). But there is
  no hidden fine S̃ structure being lost: the S̃ interior is genuinely flat,
  DS thickening is boundary-localised. Fixing the conversion would scale
  amplitude, NOT synthesise the missing fine chains — because they are not
  in S̃ either.
- **H3 ARCHITECTURE — CONFIRMED, decisive.** C1's coarseness is by design.
  Fine morphology is injected by the downstream anisotropic-FBM upscale,
  **oriented by C1's S̃ gradient** — which C1 provides (the validated
  large-scale wedge structure). The faint altitude detail at 64²/256² is
  EXPECTED and irrelevant to the design: the upscale consumes the S̃
  gradient orientation, not altitude detail.

## Consequence for the order of work

The "rich morphology" the closure-morphology caveat flagged appears **after
the upscale**, not in the C1 64² render. Therefore:

1. **The real next chantier is wiring C1 → Phase B (upscale + HD erosion)** —
   the existing-but-unwired path — NOT raising C1 resolution, and NOT a
   piste-4 render fix at 64² (the 64² altitude is already faithful to the
   coarse structure C1 is meant to produce).
2. **Piste 4 at 64² remains valid** for validating the COARSE plate-scale
   structure + the S̃ field that feeds the upscale — but it will not, by
   itself, reveal fine terrain. Set expectations accordingly.
3. **Elevated priority for the curtain follow-up:** the upscale orients FBM
   by the S̃ gradient; the grid-aligned curtain oscillation (worse at higher
   res) would imprint grid-aligned anisotropy into the upscaled detail. The
   curtain should be addressed before (or as part of) the C1→upscale wiring,
   not deferred indefinitely.

Phase 3 (fine morpho closures: arcs, margins, basins) is **not blocked**,
but its added S̃ richness will likewise only express through the upscale —
so wiring + validating the C1→upscale path is the prerequisite that makes
both Phase 3 and piste 4 morphology visible.
