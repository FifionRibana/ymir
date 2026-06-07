# HD export — the eye on the final product (C1 → upscale_from_c1 → 1024²)

Quick export (no UI) through the **contract** function `upscale_from_c1`
(laundered altitude, not raw S̃), seeds 42 / 1988 ("best relief") / 4138
("oceanic world"), 1024². Two upscale configs to test the upscale's own
pattern-breaking knobs: `default` (stock, `domain_warp_strength=0.0`) vs
`warped` (`domain_warp_strength=0.6`, `amplitude_base=0.16`). Cost: ~40 ms
default, ~55 ms warped at 1024² (so 4096² ≈ ~0.6–1 s — cheap). Harness:
`c1_closure_morphology::export_hd_upscaled`.

## What the eye sees (the three questions)

1. **Satisfying? NO — one dominant flaw: the COASTLINE is blocky.** The 64²
   plate polygon is scaled to 1024² as visible stairstep coastlines (+ a thin
   cyan fringe). The upscale adds HEIGHT detail but does NOT break up the
   coast geometry. This is the headline problem on all three seeds.
2. **FBM credible or fake? Credible (not fake) — but weak at default.** At
   `amplitude_base=0.08` (default) the interior is a smooth gradient blob; at
   `0.16` + domain warp (warped) the interior becomes **plausibly rugged**
   mountain texture — NOT obviously fake. The prior "FBM = fake detail"
   concern does not fire here: raised-amplitude FBM holds at the eye. But the
   warp/amplitude knobs improve the INTERIOR only — they do **not** fix the
   coastline.
3. **Flat interior (Lecture A) visible? Yes, but secondary + partially masked.**
   The unsculpted plateau shows as a smooth highland blob at default; stronger
   FBM (warped) masks it acceptably. It is NOT the dominant flaw — the coast is.

(Seed 4138 also shows the F3 ocean age-discontinuity line as a dark dotted
arc — minor, already registered.)

## What the export DICTATES (direction)

The eye points at a target the earlier 3-way framing (Phase 3 / sculpting /
nothing) did NOT list: **coastline refinement is the #1 lever**, and it is
distinct from both:

- **NOT Phase 3 (tectonic morpho: arcs/margins/basins):** adding S̃ structure
  would not break the blocky coast — the coast would stay a 64² stairstep.
- **NOT (primarily) the sculpting chantier (flat interior):** stronger FBM
  already masks the interior acceptably; the interior is secondary.
- **The coastline is the most resolution-starved feature.** Bilinearly
  interpolating the 64² altitude and thresholding at sea level produces a
  64²-resolution coast; `submarine_damping=0.3` further suppresses FBM below
  sea level, so the coastal contour is barely perturbed.

**Next lever to SCOPE (measure, don't assume):** break the coastline geometry
— candidate paths (W7 before any code): perturb the coarse altitude near sea
level BEFORE thresholding; coastal FBM that displaces the sea/land boundary;
reduce `submarine_damping` so coastal amplitude survives; or refine the coast
at higher resolution. The throwaway-config test already shows the in-upscale
knobs (warp, amplitude) do NOT suffice for the coast — so the fix is a
coastline-specific treatment, not a global amplitude bump.

**Positives banked:** the upscale produces a credible same-world product
(robustness held, #147 #6); interior FBM is plausible at raised amplitude;
cost is cheap (≤ ~1 s even at 4096²). Phase 3 and the sculpting chantier
remain valid but are NOT what the HD most needs first — the coast is.
