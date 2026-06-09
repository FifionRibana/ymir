# F3 — dark dotted lines in the ocean (#151)

Reported residual after the coastline work: thin dark dotted/wandering lines
in the OCEAN on several seeds (2, 99, 4138).

## Cause PINNED (`pin_f3_ocean_lines`)

Rendering the COARSE 64² fields (pre-upscale) for a strong-arc seed (4138):

- The dark arc is **already in the coarse production altitude** (Stein-Stein
  bathymetry) → it originates in C1, NOT the upscale.
- The **age field is ≈0 almost everywhere except a sparse chain of cells with
  age ≈ 447** (vs background ~0), lying exactly on the dark arc.
- That arc **coincides with a convergent plate boundary** (boundary overlay).

So F3 = the **flux-form age-advection pile-up** (the registered
density-vs-Lagrangian artefact: a few cells accumulate ~1000× the background
age at convergent boundaries) → **Stein-Stein turns those age spikes into the
deepest cells** → dark dotted lines. Background ocean is age≈0 (no seafloor
spreading to build an age gradient), so the spikes stand out starkly.

## Fix — despike the age before Stein-Stein (always-on)

A 3×3 **median** on the age field, in `c1_production_altitude`, before
`apply_stein_stein_bathymetry`. Median kills sparse spikes while preserving
smooth structure.

**Verified (before/after age colormap, seeds 2 & 4138):** only the spike chain
disappears; the (faint) background age structure is preserved — no legitimate
age gradient is flattened (the C1 model has none without seafloor spreading).
So always-on is safe (no trade-off → not gated behind a flag, unlike the
coast-warp knob which did have one).

**Scope:** Stein-Stein is the ONLY consumer of age on the altitude path
(`c1_production_altitude` → viz render + upscale), so this fixes F3 **fully**
for both. RENDER-side band-aid: the internal advected age field is unchanged;
the root cause (age advection) is the deferred deep fix
([[feedback_age_advection_density_vs_lagrangian]] — Lagrangian age / ridge
crust-creation). Not byte-identical (despiked age changes Stein-Stein) — but
consistent with the other #151 render changes (FBM resolution, coast warp);
flag in the PR. Regression `upscale_from_c1_structure_converges` still passes
(structure r unchanged; despike only smooths ocean bathymetry).

**Result:** dark ocean lines gone on all affected seeds (2, 99, 4138); ocean
smooth, land/mountains unchanged.
