# #155 drainage — flat-routing artifact diagnostic

The drainage visual validation confirmed geometrically-true rivers on slopes (all
seeds) but revealed an artifact on FLAT interiors at 2048² (invisible at 1024²,
grows with resolution): combs of parallel bars + 45° fans/diamonds in the flat
interiors (cratonic plateaus, basin floors) of 1988 / 2026.

Probe: `c1_closure_morphology::probe_flat_routing_diagnostic` (`#[ignore]`).

## Root — confirmed by counterfactual (not presumed)

`terrain::flow::pit_fill` raises a filled flat by `filled[parent] + eps` (eps=1e-7),
so a flat receives a micro-gradient following the **flood tree** (the BFS order from
the outlet), NOT real geometry. `compute_d8` then routes steepest-descent along that
eps-gradient → parallel bars (flood-tree fronts) + 45° fans (D8 diagonals of the tree).

**eps=0 counterfactual (non-invasive):** steepest-descent on the ORIGINAL (unfilled)
heightmap = eps=0 semantics. Flat cells that get a direction from `compute_flow` but
`DIR_NONE` at eps=0 are purely eps-fill-routed = the artifact. Measured:

| seed @res | flat (% land) | eps-driven (% flat) | filled / native (eps-driven) |
|-----------|---------------|---------------------|------------------------------|
| 1988@1024 | 10.9k (≈9%)   | 1291 (10.9%)        | 1269 / 22                    |
| 1988@2048 | 48156 (5.2%)  | 5256 (10.9%)        | 5056 / 200                   |
| 2026@1024 | 10732 (3.8%)  | 1042 (9.7%)         | 1037 / 5                     |
| 2026@2048 | 43621 (3.9%)  | 4674 (10.7%)        | 4547 / 127                   |

→ **~11% of flat cells are eps-driven** (direction vanishes at eps=0). Root confirmed.

## Family split — the artifact is ~96% on FILLED depressions

Flat cells are ~half filled (pit-filled) / half native (plateau). But the eps-driven
ARTIFACT is **~96% on FILLED depressions** (the pit-filled basins around lakes), only
~4% on native cratonic plateaus. The bars/fans are a **filled-depression phenomenon**.

**Native plateaus** drain via real FBM micro-relief (only ~4% eps-driven); ~13–18% of
native-flat cells exceed the stream km² threshold (minor FBM-driven streams, NOT the
eps artifact). They do NOT need false channels invented — they already drain, and the
threshold keeps most diffuse. **A flat plateau has DIFFUSE drainage by nature** (that's
why it's flat); the fix must not carve fake channels there.

**Resolution growth:** flat fraction ~4–5% of land at both 1024 and 2048 (stable %);
eps-driven ~11% of flat at both. The artifact VISIBILITY grows because absolute pixel
counts scale with grid² — enough pixels at 2048 to see the bar pattern.

## Fix scope (Garbrecht-Martz / Barnes flat resolution — for the next maillon)

Replace the `pit_fill` epsilon increment with a proper, deterministic, knob-free flat
resolution (Garbrecht-Martz 1997, or the Barnes-Lehman-Mulla 2014 priority-flood flat
routing): impose a flat gradient combining **distance-away-from-higher** (the inflow
cells spilling into the flat) + **distance-toward-lower** (the outlet the flat drains
to) → convergent natural drainage to the sill, deterministic, replacing the eps tree.

- **Filled depressions (~96% of artifact):** the main beneficiary — convergent
  drainage to the outlet sill instead of BFS bars.
- **Native plateaus (~4%):** no special handling — already diffuse via FBM; the km²
  threshold gates channels, so GM routing won't carve fake rivers (it makes the routing
  coherent, the threshold decides what's a channel).

**Render-side mitigation already in hand** (NOT the deep fix): the filled depressions
are largely LAKES, so masking rivers under lake cells (the viz convention, applied in
`probe_drainage_overlays`) hides most of the artifact visually. But the underlying D8
directions / accumulation on filled flats are still BFS-tree-wrong — flat resolution is
the correctness fix (matters for accumulation, and for shallow filled flats not deep
enough to register as lakes).

Diagnostic only — the fix follows this verdict.

---

## FIX LANDED — Garbrecht-Martz flat resolution (FEAT 0278119, TEST 7d54afc)

`pit_fill` now fills to the EXACT sill (no `+eps`); `resolve_flats()` imposes the
convergent gradient: `flat_grad = tl·(fhmax+1) + (fhmax − fh)` (toward-lower BFS `tl`
dominant → guaranteed descent, no interior minimum; away-from-higher `fh` breaks the
bar-forming ties). `compute_d8` pass 1 = steepest descent on `filled` (real slopes
unchanged), pass 2 = flat cells route down `flat_grad` to the outlet.

**Self-targeting:** only EXACT-equal flats (pit-filled depressions) are resolved;
native FBM-textured plateaus are never exactly flat → untouched (no invented channels).
`filled` keeps the exact sill → lake levels unchanged.

**Acceptance (probe_flat_routing_fix_render, 2048²):**
1. **Artifact dissolved** — 1988/2026 interiors: parallel-bars/fans → dendritic networks.
2. **Native plateaus intact** — diffuse FBM drainage, no false straight channels.
3. **Slopes unchanged** — 1337 control + margins of 1988/2026 (steepest-descent path
   untouched; eps removal doesn't affect cells with a real lower neighbour).
4. **All resolutions / deterministic** — BFS, no RNG; the resolution-dependent artifact
   is gone at 2048² (the worst case); lib 450/0, all C1 green, only pre-existing v2 red.
5. **Lakes unchanged** — `filled` = exact sill; the fix changes routing, not filling.
6. Residual: convergent fill only in the LARGEST basins (= lakes, hidden by the
   product's lake-masking) — not the systematic artifact.

This CLOSES the drainage maillon.

---

## Quasi-flat residual — measured MARGINAL, documented (probe_quasi_flat_residual)

The visual flagged a faint ladder/parallel residual on QUASI-flats (low gradient but
NOT exactly at the sill → not exact-equal → `resolve_flats` skips them by design). The
local gradient (the instrument the hillshade lacked) quantifies extent + nature:

**Extent (grad <1e-4 norm/cell AND non-exact-flat):** 1988 **0.281 % of land** (3.0 % of
drainage); 2026 **0.141 % of land** (1.4 %). Marginal.

**Nature — no systematic defect:** directional coherence (fraction of drainage
neighbours sharing D8 direction, 5×5) is **~0.50, flat across ALL gradient bins**
(0.40–0.58), including the quasi-flat bins. A parallel/ladder artifact would spike the
low-gradient bins toward the planar baseline (synthetic 0.6 m/km slope → **0.92**); they
don't → drainage is dendritic at every gradient. 88 % of drainage is at franc gradient
(≥2.3 m/km). The synthetic counterfactual also showed planar-slope parallel flow
doesn't accumulate above the stream threshold (so it barely renders as channels).

**Type:** ~50 % of the residual is fringe-of-depression (sill edges just above the
exact-equal flat), ~50 % very-gentle plains (natural D8-planar parallelism).

**Verdict: DOCUMENT, do not fix.** 0.14–0.28 % of land, no gradient-localised defect,
half of it natural parallelism. Extending the criterion (near-sill fringes / a gradient
band) would risk over-correcting correct gentle slopes for a ~0.1 % gain — the
anti-over-correction discipline: the number doesn't justify it. A future flat-routing
refinement (fringe inclusion) is logged but NOT blocking. Drainage maillon stays closed.
