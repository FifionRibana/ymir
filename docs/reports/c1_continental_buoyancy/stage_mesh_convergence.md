# Mesh convergence (resolution invariance) of C1 — measurement (#145 follow-up)

**Chantier:** a well-posed model produces the SAME structure at every mesh
(finer, not different). Measure whether C1's structure CONVERGES or
DIVERGES with resolution — a PREREQUISITE before wiring the upscale (which
reads the S̃ gradient) and before piste 4. Measure + test the proximity
hypothesis BEFORE scoping any fix. Register interior/boundary sculpting as
a DISTINCT downstream chantier (not treated here).

**Measure:** sweep 64²/128²/256²/512², seed 42, full rigid stack, physical
time held constant (`n_steps = 300·grid/64`, since CFL `dt ∝ 1/grid`).
Harness: `c1_closure_morphology::mesh_convergence_sweep` (`#[ignore]`).
Convergence test = block-average the high-res field down to 64² and
correlate (Pearson) with the native 64² field; high r ⇒ large-scale
structure is mesh-stable.

```
 grid  n_steps land% largest | wedge% wedgeN d̄coast | S̃→64 r alt→64 r | time
  64²    300    27.7  0.937  |  2.1    24    0.0000 | 1.0000  1.0000  | 0.28s
 128²    600    27.0  0.963  |  0.5    28    0.0000 | 0.5442  0.9175  | 1.9s
 256²   1200    28.1  0.968  |  0.2    49    0.0000 | 0.4704  0.8799  | 13s
 512²   2400    28.1  0.970  |  0.1    72    0.0000 | 0.4574  0.8718  | 101s
 boundary perim/A:  init≈final at every grid (Δ≤0.003); IoU(final,init)=0.95–0.99
   64²:0.256  128²:0.116  256²:0.057  512²:0.029   (≈ ∝ 1/grid)
```

## Two DISTINCT signatures (they belong to two different chantiers)

### Signature 1 — S̃ thickening is NON-CONVERGENT → THE invariance failure

- **S̃→64² correlation plateaus at ~0.46** (1.0→0.54→0.47→0.46) — it does
  NOT approach 1. The thickness field does not converge.
- **Altitude→64² correlation stays high (~0.87)** — large-scale geography
  IS mesh-stable. The non-convergence is CONFINED to the part of S̃ that
  altitude expresses weakly (the wedge/curtain detail).
- **Wedge area (S̃>1.5) collapses ∝ 1/grid** (2.1→0.5→0.2→0.1 %) with
  **d̄coast = 0.0 at every resolution** and component count growing 24→72.
  This is the signature of a **1-cell-wide margin pile** pinned ON the
  coast: oceanic accretion / Davis-Suppe thickening against the rigid
  margin, whose physical width → 0 as the mesh refines (so its gradient
  blows up and never converges).
- Visual (512² S̃): dense grid-aligned curtain speckle across the interior
  + diagonal striping along the velocity field, **markedly worse than
  64²** — the bounded curtain oscillation is itself mesh-dependent.

**Proximity hypothesis — refined, not as originally stated.** The diverging
feature is NOT a relief chain born from a convergence boundary *merging*
with a nearby coast at low res (the wedges sit exactly ON the coast,
d̄coast=0, at ALL resolutions — there is no separate convergence boundary
moving away as the mesh refines). It is a **coast-pinned, grid-width
accretion/DS pile** whose area ∝ 1/grid. Same root cause family
(low-res-forced numerics), different geometry than the merge story.

**Mesh-dependent constants implicated (step 2 of the chantier):**
1. **DS + accretion deposition is per-CELL, not per-physical-LENGTH** (#1
   culprit). The thickening rim is ~1 cell wide at every grid → physical
   width shrinks with mesh → area ∝ 1/grid, ∂S̃ diverges. Fix = deposit a
   thickness PROFILE of fixed physical width (multiple cells at high res).
2. **No-flux boundary / curtain oscillation** — grid-aligned, worse at
   512². Mesh-dependent by construction (upwind on a 1-cell-sharp rigid
   contrast).
3. Erosion K (stream-power, notoriously mesh-dependent) and upwind
   numerical diffusion (∝ dx) — secondary; the geography converged, so
   these are not the dominant failure here.

### Signature 2 — boundary shape & interior: NO dynamic imprint → SCULPTING (Lecture A)

The `/btw` second signature (smoother boundaries at high res, recognisable
as the init Voronoi polygons). A/B test:

- **IoU(final, init) = 0.95–0.99** — the final continental mask is the INIT
  plate geometry; the dynamics barely move the coast.
- **perim/A final ≈ init at every grid** (Δ ≤ 0.003) and **halves per
  doubling (∝ 1/grid)** — the geometric signature of a FIXED smooth polygon
  sampled finer, NOT a fractal/tectonically-roughened coast (which would
  hold perim/A roughly constant). 512² S̃ shows literal straight-edged
  Voronoi-cell boundaries. (`/btw` tests #2 + #3.)
- **`/btw` test #1 — origin of the 64² edge irregularity** (`boundary_ab_test_64`):
  - **(1a) coast jaggedness UNIFORM:** ocean-neighbours/coast-cell =
    **1.447 at convergent-adjacent coast vs 1.347 elsewhere (ratio 1.07)** —
    the crenellation is essentially uniform along ALL coasts, independent
    of tectonics ⇒ grid pixelation (A).
  - **(1b) coast changes faintly tectonic:** only **12 / 204 coast cells**
    changed (init≠final, matches IoU 0.99), of which **10 (83 %) at
    convergent boundaries vs 64.7 % base → enrichment 1.29×**. A small but
    real edge imprint exists (the ≥2-rule subduction promoting a few
    ocean→land cells at convergent margins).

**Verdict A/B = predominantly Lecture A, faint second-order B residual.**
The dominant signal (uniform jaggedness ratio 1.07, IoU 0.99, perim/A→init
∝1/grid) is decisive: the 64² edge irregularity is grid noise that higher
resolution honestly cleans, and the boundary stays at the un-reworked init
Voronoi geometry — the edge dynamics do not meaningfully shape the coast.
The 1.29× change-enrichment confirms a *faint* real tectonic imprint (12
cells, 6 % of coast) but it is second-order and NOT shown to be
mesh-divergent — not an edge-dynamics convergence failure worth its own
invariance work. So the smooth-boundary signature belongs to the
**SCULPTING chantier** (closures don't rework boundaries OR interior),
**SEPARATE and AFTER invariance** — exactly as pre-registered for the
interior, now extended to the boundary. (Lecture B as stated — real
tectonic shaping LOST at high res — is rejected: there is almost no edge
shaping at 64² to lose.)

## Convergence verdict & fix sizing

- **Large-scale GEOGRAPHY already converges** (alt r ~0.87, land% 27.7→28.1,
  largest 0.94→0.97, IoU≈init) — continents/bathymetry are mesh-stable
  because they are the (convergent) init Voronoi geometry + Stein-Stein.
  Structural-convergence criterion is MET for geography.
- **The invariance FAILURE is localised** to S̃ margin thickening
  (DS+accretion per-cell deposition) + the no-flux curtain. The field the
  upscale would consume (S̃ gradient) is the non-convergent one → wiring the
  upscale now would make 64²+FBM and 512²+FBM different worlds. Confirmed
  prerequisite.
- **Fix is SCOPED, not a wholesale recalibration:** re-express DS +
  accretion deposition in physical width (per-length, fixed-width profile)
  and fix the curtain; then re-run this sweep and require S̃→64 r to climb
  toward ~1 and wedge% to stabilise (not ∝1/grid). The 64² cap/n_cycles
  geography calibration is NOT implicated (geography already converges).

## Registered as DISTINCT (not treated here)

**SCULPTING chantier (boundary + interior), AFTER invariance.** Once the
margin thickening is mesh-convergent, the interior will be revealed as a
flat plateau and the boundary as the un-reworked init polygon (Lecture A).
The dynamics (DS only at margins; Track D barely imprints) do not sculpt.
The fix must be PHYSICAL (erosion that incises, intracontinental rifting,
distributed deformation), NOT FBM pasted on. Causal order: invariance
first (it REVEALS the true flat/un-sculpted state, which at 64² is masked
by the grid-width margin pile overflowing inward); sculpting second.

Anti-patterns honoured: did not aim for bit-identity (criterion is
structural convergence); did not wire the upscale; kept non-convergence
(invariance) and missing-sculpting (separate) distinct; measured + tested
the hypotheses BEFORE proposing the fix.
