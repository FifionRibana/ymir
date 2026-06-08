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

## Recommended production setting
`coast_warp_strength ≈ 0.8`, `coast_warp_frequency = 0.5` for the C1 HD export.
Kept OFF by default in `FbmUpscaleConfig` (byte-identical / v2-safe); the C1
production HD path opts in. Wiring the default into the C1 production upscale
config is a product decision (deferred with the rest of the HD/UI wiring).
