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
