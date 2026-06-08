# Coastline warp (#151 follow-up) — STEP-1 fix, eye-validated

Cause was pinned to STEP 1 (`stage_hd_export.md`): the sea-level contour
follows the bilinear-interpolated 64² altitude; the domain warp only warps the
NOISE, never the coarse sampling. **Fix (piste 1):** displace the
COARSE-ALTITUDE sampling coords `(sx,sy)` before `sample_bilinear_periodic`,
in coarse-pixel units, so the contour meanders.

## Implementation (minimal, reuse, opt-in)
`FbmUpscaleConfig.coast_warp_strength` (coarse cells) + `coast_warp_frequency`
(cycles/coarse pixel). Dedicated coast-warp FBM (distinct seeds from the
noise/domain-warp). Applied to the sampling position of `base_height` + slope +
direction (coherent). **Default 0.0 = OFF = byte-identical** (the v2 pipeline
is unaffected). The domain warp (noise) is untouched.

## Eye verdict (the judge — credible meander vs procedural-fake)
`export_coast_warp`, seeds 42 + 1988, strength sweep {0, 0.5, 0.8, 1.2} coarse
cells, 1024²:

- **off (0.0):** the blocky 64² stairstep coast (axis-aligned 90° steps).
- **0.8 (sweet spot):** the stairstep is broken into an **irregular,
  credibly-meandering coast** on BOTH seeds — FBM-irregular, NOT regular
  procedural ripples (passes the "not fake" bar; real coasts are irregular).
  No fragmentation; continents stay connected.
- **1.2:** more meander but starts to **over-warp** — a few coastal cells
  detach into pixel-scale fragments. Past the useful range.

**Verdict: piste 1 works.** A coast warp of ~0.8 coarse cells dissolves the
64² stairstep into a credible irregular coastline, eye-confirmed on two seeds,
without fragmentation. The displacement is ~0.8 coarse cells ≈ 13 px at 1024²,
which breaks the ~16 px stairstep. Honest bound: this fixes the PLATE-SCALE
coast geometry; sub-~13 px coastal filigree would need finer warp octaves or a
finer base — but the blocky-coast complaint (the #1 HD flaw) is resolved.

## At the TARGET resolution (4096², seed 1988 — `export_coast_warp_4096`)

~1.7 s per 4096² render. 0.8 coarse cells ≈ 51 px here.

- **off:** even MORE glaring than at 1024² — clean ~64 px axis-aligned
  stairstep blocks. Unacceptable at target res. Confirms the fix is necessary,
  not cosmetic.
- **c08 (+ amplitude 0.16):** credible continent — irregular meandering coast
  (no axis-aligned steps) + rugged mountain interior. Residual: some
  large-scale 64²-polygon straightness remains (the warp wobbles the broad
  shape rather than redrawing it). Plate-scale-acceptable; reads as a real
  coast at this zoom.
- **c12:** stronger meander, breaks the residual large-scale straightness MORE,
  and — unlike at 1024² — shows **no fragmentation** at 4096² (finer pixels
  render the displaced coastal cells smoothly). So warp tolerance is
  **resolution-dependent**: 1.2 fragments at 1024² but is clean + more natural
  at 4096².

## FBM resolution-calibration fix (#151) — the FBM was target-pixel-based

Surfaced by the 4096² render: the SAME config gave a **4× finer** interior FBM
at 4096² than 1024². Cause (pre-existing, NOT the coast warp): the noise was
sampled in TARGET-pixel space — `freq = base/src_w; nx = i·freq` (i = target
pixel) → the FBM cycle count across the image = `(dst/src)·base`, scaling with
`target_size`. The coast warp was already coarse-cell-based (resolution-OK);
only the interior FBM was miscalibrated.

**Fix:** sample ALL noise (base FBM, domain warp, angle perturbation) in
COARSE-CELL space (`sx·nscale`, `nscale = base·NOISE_REF_TARGET/src²`,
`NOISE_REF_TARGET = 1024`). The feature size is now fixed relative to the
terrain, independent of `target_size`. Coefficients reference the prior 1024²
calibration so **1024² is byte-identical to pre-#151** (`sx·nscale ≡ i·freq`
when `target == 1024`; warp off by default) — v2@1024 unaffected; other target
sizes now MATCH the 1024² feature size instead of diverging.

**Verified:** post-fix, c08 at 1024² and 4096² show the same world at the same
relative feature size (broad mottled interior at both); pre-fix 4096² was
visibly finer. Upscale lib tests green (default-off / 1024-identical).
**v2 impact:** v2 renders at a non-1024 target would change (resolution-
corrected) → re-validate v2 visually before relying on it; no v2 test asserts
exact pixels.

## v2 re-validation (the FBM fix touches shared code)

v2 consumers of `upscale_with_fbm`: v2 viz (target = actual heightmap dims,
dropdown) and v2 workflow `phase_b` (`hd_grid_size`, **default 2048**,
workflow/mod.rs:211). So v2's default target is **2048²**, not 1024² — the fix
DOES change v2's default output (the pre-existing resolution-dependence bug
affected v2 too). `export_fbm_2048_isolate`, seed 1988, default config (coast
warp off, to isolate the FBM), rendered before vs after by swapping `upscale.rs`:

- **before (pre-fix) @2048²:** interior FBM slightly finer/granular (2× the
  1024 cycle count — the bug).
- **after (fixed) @2048²:** interior slightly smoother/broader = the
  1024-referenced feature size; **well-formed, well-calibrated continent**
  (amplitude/aniso/gradient read correctly).

**Verdict: v2 IMPROVED / NEUTRAL.** The fix gives resolution-consistency; at
the 2048² default it is a subtle, benign coarsening (2×), NOT a calibration
derangement. v2's calibration was not critically relying on the
resolution-dependent fineness. Residual (non-blocking): this isolation used a
C1 coarse + default params (the FBM feature-size change is source/param-
independent); a v2-NATIVE before/after at 2048² would be the gold-standard
confirmation before relying on v2 at non-1024 targets in production.

## Recommended production setting
`coast_warp_strength ≈ 0.8`, `coast_warp_frequency = 0.5` for the C1 HD export.
Kept OFF by default in `FbmUpscaleConfig` (byte-identical / v2-safe); the C1
production HD path opts in. Wiring the default into the C1 production upscale
config is a product decision (deferred with the rest of the HD/UI wiring).
