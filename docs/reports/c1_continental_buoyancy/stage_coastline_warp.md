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

## Border-feathering artefact (#151) — cause PINNED

User report: the FBM "deforms the interior near the border" (directional
combing at the coast); macro border + coast warp OK, mountains OK. Native-res
crops (`export_fbm_2048_isolate` + `save_heightmap01_crop`, the artefact is
invisible under Read's downscaling) confirm: a feathered/combed coastline,
worst on the thin near-sea-level bridge/peninsula.

Counterfactual pinning (2048², seed 1988, env-driven):
- `max_anisotropy = 1` (isotropic): **still feathered** → NOT anisotropy.
- `amplitude_slope_factor = 0` (no slope boost): **still feathered** → NOT the
  slope-amplified amplitude.
- `amplitude_base 0.08 → 0.02`: **feathering nearly gone** → the **base FBM
  amplitude crossing the sea-level contour** is the cause.

Mechanism: `final = base_height + amplitude·noise`; near the sea-level contour
`±amplitude·noise` flips cells land↔ocean → a feathered coast. Worst where the
coarse altitude is FLAT and near sea level (the bridge): no relief to anchor
against, so the whole low region is FBM-dominated. Strong-relief interior
(mountains) is barely perturbed → "mountains OK". `submarine_damping` only
damps BELOW sea level, not the near-sea-level band on the land side.

**Fix direction (NOT global amplitude — that flattens mountains too):** a
LOCAL near-sea-level amplitude taper — `amplitude *= smoothstep(|base_height −
sea_level|, 0, coastal_band)` → FBM ≈ 0 at the coast (no feathering), full
inland (mountains preserved). The coast warp supplies the macro coast
irregularity; this removes the micro-feathering. `coastal_band` ~0.05–0.1
altitude.

**FIX DONE + validated** (`coastal_amplitude_band`, opt-in, default 0.0 =
byte-identical / v2-safe; upscale lib tests green). Native-res crops, band
sweep {0, 0.06, 0.12}: **band 0.06 removes the coastal feathering** (coastline
back to a clean contour) while the interior mountain texture is **fully
preserved** (local taper — `smoothstep(|h−sea|,0,band)` = 1 beyond the band, so
inland relief is untouched). band 0.12 also clean (slightly wider flat coastal
strip). Recommend **~0.06** (gentler). Complements the coast warp: the warp
bends the contour (macro irregularity), the taper stops the FBM
micro-feathering at it.

## Recommended production setting
`coast_warp_strength ≈ 1.5`, `coast_warp_frequency = 0.5`,
`coastal_amplitude_band ≈ 0.30`, `amplitude_base ≈ 0.16`,
**`submarine_damping = 0.0`** for the C1 HD export.

(`coast_warp_strength` 0.8 → **1.5**: with FBM removed from the coast (band 0.30),
the warp carries the coastline irregularity ALONE → 0.8 read too light/smooth.
Sweep {0.8, 1.2, 1.6, 2.0}: 1.2–1.6 give a natural irregular coast, 2.0 gets
busy/spiky; **1.5** chosen. Verified across geometries — seeds 1988, 2
(high-relief), 42 (modest), 4138 (islands): irregular natural coasts, NO
fragmentation, smooth coastal plains, rugged mountains preserved. The two coast
levers are now correctly balanced: the warp does the macro irregularity, the
band keeps the plains smooth.)

(`coastal_amplitude_band` 0.20 → **0.30**: a 3rd reported pass showed faint
reddish contour lines on the near-coast land. **Decisive attribution** — a
HILLSHADE of the region shows the actual relief there is GENTLE/smooth (no
sharp structures); the lines are the hypsometric TEST PALETTE's hard green→brown
colour boundary (at normalised height 0.75) being wiggled by FBM that resumes at
the band edge (0.70 for band 0.20) — i.e. a RENDERING artefact of the test
palette, not a terrain defect (the heightmap data is good). Two consequences:
(1) the production renderer/colormap determines whether this shows at all; (2)
for clean test renders AND a sensible model, band 0.30 puts the FBM onset at
height 0.80 — ABOVE the green/brown boundary — so all green lowland is smooth
and FBM lives only in the brown highland. Verified: green coastal margins clean
on seeds 2/42, mountains (seed 2 high-relief) fully preserved. So "FBM ∝
elevation: smooth plains, rugged mountains" with FBM confined to highland.)

(`coastal_amplitude_band` revised 0.06 → **0.20**: a second reported defect —
faint FBM contour ripples on the LOW land NEAR the coast (confirmed real height,
not palette banding, via a grayscale crop). The 0.06 band killed only the
coastline flicker; the near-coast lowland still carried FBM. Widening to ~0.20
suppresses FBM on land within 0.20 altitude of sea (coastal plains → smooth)
while keeping full FBM on the higher inland terrain — verified on seed 2
(high-relief: mountains fully preserved, snow peaks intact) and seed 1988. This
makes FBM amplitude effectively scale with elevation: smooth coastal plains,
rugged mountains — geologically sensible. Pure config value, no code change.)

The last one fixes a separate reported defect: with the default
`submarine_damping = 0.3` the FBM is still visible in the OCEAN (~0.3× amplitude
→ mottled bathymetry / contour ripples — "on voit le FBM, on ne devrait pas").
Counterfactual: `0.3 → 0.0` removes it entirely (smooth bathymetry gradient;
land/mountains unchanged — FBM above sea level is untouched). So set
`submarine_damping = 0.0` for C1 HD (no underwater FBM). Combined production
render (7 seeds: 2/42/99/1337/1988/2026/4138) = credible continents, irregular
coasts, mountains w/ snow peaks, **smooth oceans**. (Residual: F3 ocean
age-discontinuity lines on a few seeds — separate registered follow-up, not the
upscale.)

All these are existing/`FbmUpscaleConfig` fields kept OFF/default
(byte-identical / v2-safe); the C1 production HD path opts in. Wiring these
defaults into the C1 production upscale config is a product decision (deferred
with the rest of the HD/UI wiring).
