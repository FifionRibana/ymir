# C1 mesh-invariance — fix scope (W7 surface, NOT yet implemented)

Issue-in-waiting: **"Resolution invariance (mesh convergence) of C1"**. Built
on the measurement verdict in `stage_mesh_convergence.md`. This document is
the *surface-before-implement* (W7) scoping — it locks the mechanism, the
physical anchor, and the two open points BEFORE any code.

## Problem (from the measurement)

C1's large-scale GEOGRAPHY already converges (alt→64² r ~0.87, IoU≈init,
land% stable). The invariance FAILURE is localised to the **S̃ thickening
field** (S̃→64² r plateaus at 0.46; wedge area ∝ 1/grid; curtain worse at
512²). Because the upscale consumes the **S̃ gradient**, a non-convergent S̃
means 64²+FBM and 512²+FBM are different worlds → invariance is a
prerequisite to wiring the upscale and to piste 4.

## Root cause located in code (the #1 culprit, confirmed)

`crates/ymir-core/src/tectonics_c1/closures/davis_suppe/source_term.rs`
defines the wedge length scales **in CELLS, not physical length**
(verbatim doc, lines 69–77):

```
l_taper: 4.0,       // "Characteristic length … In cells, not non-dim length."
l_decay: 6.0,       // "… In cells."
max_distance: 30.0, // "Outside this distance (cells) … zero."
```

`h_critical(d) = h_max·(1 − exp(−d/l_taper))` with `d` in cells ⇒ the wedge
reaches `h_max` over a FIXED CELL COUNT. Physical wedge width = `l_taper /
grid` → shrinks ∝ 1/grid as the mesh refines. That is precisely the
measured `wedge% ∝ 1/grid`, coast-pinned (`d̄coast=0`), `∂S̃` divergent
signature. **Per-cell, not per-physical-length — root cause, three lines.**

## Fix #1 (principal) — re-express DS length scales in PHYSICAL units

Convert `l_taper`, `l_decay`, `max_distance` from cells to non-dim domain
length (the `[0,1]` unit domain); consume them as cells via `× grid` at the
call site (the wedge-distance field is already in cells, so a single
conversion `l_phys · grid` keeps the kernel unchanged). Result: a wedge of
FIXED PHYSICAL width, sampled by more cells at higher resolution — the
invariance condition.

### Two objectives the original scope conflated — SEPARATE them

A physics check on the code corrects the original Point b. `h_critical(d) =
h_max·(1 − exp(−d/l_taper))` is an **empirical exponential** rise (the doc
calls it "the classical Davis-Suppe wedge profile", but mechanically it is a
plausible-shape exponential, NOT the critical-wedge mechanics of
Davis-Suppe-Dahlen, where the taper = α_surface + β_décollement).
`l_taper` is the exponential's characteristic length, NOT "the run-out at
slope tan(α_c)".

- **Objective 1 = Fix #1 = INVARIANCE (this issue).** The measured bug is a
  UNIT bug: `l_taper/l_decay/max_distance` in CELLS → physical width ∝
  1/grid. Converting cells → non-dim length (× grid at the call site) fixes
  invariance **regardless of whether the formula is physically faithful**.
- **Objective 2 = PHYSICAL FIDELITY of `h_critical` (NOT this issue).**
  Replacing the empirical exponential with the true critical-wedge mechanics
  is a fidelity question, deeper, and the code already defers fine fidelity
  to "Phase 3 Lallemand".

### Point b — REVISED — anchor the physical width on the CURRENT 64² value

Do NOT derive `l_taper_phys` from `h_max / tan(α_c)`: that injects a number
from the TRUE critical-taper relation into a formula that is NOT that
physics → a homemade relation with a parameter dressed as physics (the exact
closure-relations trap). Instead, the conservative, correct anchor:

> take the CURRENT `l_taper = 4 cells` at 64² and convert to non-dim length
> (`4/64 = 0.0625` of the domain); likewise `l_decay = 6/64`, `max_distance
> = 30/64`. This PRESERVES today's already-calibrated 64² behaviour
> (geography converges) while making it mesh-invariant (same physical width,
> cells adapt with resolution).

This is purely a change of UNITS, not of behaviour at 64². The α_c / true
critical-wedge anchor belongs to Objective 2 (Phase 3 Lallemand), where the
code itself says the physics will be refined.

## Fix #2 candidates — to confirm in scope, not presumed

1. **Accretion deposition** — audit `closures/accretion/` for the same
   per-cell length pattern (it piles oceanic crust against the rigid margin,
   contributing to the coast-pinned wedge). Likely the same per-cell→
   per-length conversion. Confirm before lumping into Fix #1.
2. **No-flux curtain** — the bounded grid-aligned oscillation, measured
   WORSE at 512². See Point a.

## Point a (LOCKED as a measurement, not a presumption) — curtain coupled or independent?

Do NOT presume the curtain is a separate defect. It and the margin pile may
share a cause (sharp 1.0/0.2 contrast + upwind at the rigid face, the
oscillation root identified earlier). **Scope order:**
1. Implement Fix #1 (physical-width DS).
2. Re-run `mesh_convergence_sweep`. Inspect whether the curtain ALSO
   diminished (shared cause → one fix) or persists (independent → a second,
   separate treatment of the sharp-contrast no-flux face).
3. The "curtain counterfactual" (disable the no-flux / soften the contrast
   on a throwaway) belongs HERE — after Fix #1 isolates what remains — not
   before scoping.

## Acceptance criterion (the measurement IS the test)

Re-run `c1_closure_morphology::mesh_convergence_sweep`:
- **S̃→64² r must climb toward ~1** (from the 0.46 plateau) across
  64²/128²/256²/512².
- **wedge% must STABILISE** (a fixed physical fraction), not decay ∝ 1/grid.
- alt→64² r and geography metrics (land%, largest, IoU≈init) must NOT
  regress (they already converge).

The structural-convergence criterion (not bit-identity): large formations
same and stable, formations no longer disappear, geography recognisable at
every mesh. Fine detail may still vary (the upscale fills it legitimately).

## Fix #1 RESULT (measured) — correct + byte-identical, but INSUFFICIENT alone

Implemented (`DavisSuppeParams::scaled_to_grid(nx)`, applied at the kernel +
both wedge-distance call sites in `time_loop`). Byte-identical at 64²
(factor = 64/64 = 1; imprint tests 1.2/1.3/1.4 + DS unit tests all green,
unchanged). Re-running `mesh_convergence_sweep`:

```
 grid  S̃→64 r  (pre-#147 → post-Fix#1)   wedge%  (pre → post)
 128²   0.544 → 0.568                      0.5 → 0.6
 256²   0.470 → 0.512                      0.2 → 0.3
 512²   0.457 → 0.473                      0.1 → 0.2
```

**Verdict: Fix #1 is necessary but NOT sufficient.** S̃→64² r still plateaus
~0.47 (target ~1); wedge% still ~∝1/grid. The Davis-Suppe continental wedge
was NOT the dominant non-convergence source — the "evident suspect" (named
#1 from the code reading) is largely innocent for the FIELD-level metric.

**Revised attribution (isolate next, measure before fixing):**
- The **wedge metric (S̃>1.5) is OCEANIC ACCRETION** piled against the rigid
  margin (continental orogens were ~0 above 1.5) — `closures/accretion/`,
  NOT touched by Fix #1. Its per-cell deposition is the likely wedge%∝1/grid
  driver.
- The **curtain** (grid-aligned no-flux oscillation, dense interior speckle
  worse at 512²) is the likely driver of the FIELD decorrelation (S̃ r),
  covering the interior at high spatial frequency.

**Next (Point a, now actionable):** counterfactual sweep — disable
subduction+accretion (isolate the accretion pile) and separately probe the
no-flux curtain — to attribute the residual S̃ r / wedge% BEFORE the next
fix. Fix #1 stays (correct, byte-identical, no regression); one necessary
piece, not the whole.

## Counterfactual ATTRIBUTION (measured — `mesh_convergence_attribution`)

Three variants × {64,128,256}, rigid, post-Fix-#1, each S̃ r vs its OWN 64²:

```
 variant            wedge% (64/128/256)   S̃→64 r (64/128/256)
 A full             2.08 / 0.59 / 0.28    1.0 / 0.568 / 0.512
 B sub+acc OFF      4.76 / 3.71 / 3.78    1.0 / 0.809 / 0.777
 C advection-only   1.46 / 0.71 / 0.36    1.0 / 0.047 / 0.045
```

- **B stabilises wedge%** (3.71→3.78, NOT ∝1/grid) **and lifts S̃ r
  0.51→0.78** → **subduction+accretion carry the wedge%∝1/grid collapse AND
  ~0.27 of the field decorrelation.**
- **C crashes S̃ r to ~0.045** → the **no-flux advection curtain is
  intrinsically mesh-NON-convergent**; equilibrium-height (in A/B, not C)
  is what BOUNDS it to the ~0.5–0.8 range. The curtain is the deep floor.
- **Coupling CONFIRMED:** both at the rigid margin on the sharp contrast;
  removing sub+acc moves BOTH metrics; EH bounds the curtain.

**SPLIT B (measured) — accretion carries it, subduction innocent:**

```
 variant       wedge% (64/128/256)   S̃→64 r (128/256)
 A full        2.08 / 0.59 / 0.28    0.568 / 0.512
 B1 sub_off    2.08 / 0.58 / 0.28    0.564 / 0.508   ≡ A (subduction INNOCENT)
 B2 acc_off    5.03 / 3.84 / 4.17    0.812 / 0.781   ≡ B (accretion CARRIES it)
```

B1 (subduction off) is identical to A → subduction (the per-cell
Oceanic→Continental promotion, the #145 finger/≥2 terrain — a plausible
bet) is INNOCENT. B2 (accretion off) stabilises wedge% + lifts r → the
**oceanic ACCRETION margin pile is the wedge%∝1/grid + ~0.27-r carrier.**
Another evident suspect innocenced by measurement, not bet.

**Next (order = de-risk the hard/critical part first, NOT easy-first):**
1. **Contrast counterfactual (throwaway, BEFORE any fix).** Both accretion
   pile and curtain live on the sharp 1.0/0.2 contrast at the rigid margin
   (coupling: removing accretion moved BOTH metrics). Hypothesis: softening
   the contrast (local buoyancy transition over 2–3 physical cells) treats
   BOTH (upwind stops oscillating on a smoothed contrast; deposition has no
   sharp wall to pile against). Throwaway-soften + re-run the sweep:
   - r→~1 AND wedge% stabilises → single root cause (sharp contrast), one
     elegant fix.
   - r lifts but wedge% doesn't → two fixes (contrast for curtain, physical
     width for accretion margin).
   - even softened contrast doesn't lift r → curtain is a SCHEME floor
     (upwind can't converge on a contrast without changing advection) →
     redefine the criterion (r~0.78 realistic).
2. Then the fix the counterfactual points to — single-root (contrast) or
   two (contrast + accretion physical-width). The curtain is on the critical
   path to r→1 (EH bounds it to ~0.78 but does not converge it); the margin
   is the safe Fix-#1-pattern part. De-risk the curtain/contrast first: it
   tells us whether r→1 is ATTAINABLE at all before investing in the margin.

## Contrast counterfactual (γ) — SCHEME FLOOR (measured, `contrast_counterfactual_gamma`)

Advection-only + rigid no-flux + S̃ band-smoothing (2.5 physical cells),
reusing the real `run_advection_only` one step at a time. S̃→64 r:

```
   λ      64²      128²     256²
 0.00    1.0000   0.0466   0.0454   (= pure curtain, variant C)
 0.25    1.0000   0.0924   0.0158
 0.50    1.0000   0.1385   0.0139
 1.00    1.0000   0.2207   0.0153
```

**Verdict = pre-registered branch 3: SCHEME FLOOR.** Contrast band-smoothing
does NOT lift r toward 1 (256² flat ~0.015 at every λ, even full local
averaging). The single-root "sharp contrast" hypothesis is **rejected**.

**Deeper finding from the r≈0.045 magnitude:** near-zero correlation means
advection-only is **bulk-decorrelated**, not a thin boundary curtain (which
would leave the interior at r~0.7). The curtain speckle is a SYMPTOM; the
disease is **bulk upwind-advection mesh-dependence** — the transport itself
produces a different field at each resolution. A boundary (or even interior)
contrast-smooth cannot fix a bulk-transport mesh-dependence. The full system
reaches r~0.5 only because **equilibrium-height pins S̃ toward h_eq** (the
mesh-stable convergence pattern), NOT because the advection converges.

**Consequence — REDEFINE the criterion (don't chase r→1 via contrast; save
the (α) churn):**
- r→1 is NOT attainable with this upwind+rigid advection (it needs a
  higher-order / less-diffusive advection SCHEME — deep, out of #147).
- The realistic, achievable invariance target is the EH-bounded ceiling:
  fix **accretion** (physical-width margin deposition, the proven Fix-#1
  pattern) to recover wedge% convergence AND lift full-system r from 0.51
  toward the ~0.78 ceiling (the value with accretion's mesh-dependence
  removed). Accept ~0.78 as the realistic structural-convergence r.
- **Register (deep follow-up, out of #147):** "upwind advection is bulk
  mesh-non-convergent (r≈0.045 isolated); r→1 needs an advection scheme
  change (higher-order / less diffusive)." Equilibrium-height masks but does
  not cure it.

**#147 actionable deliverable:** Fix #1 (DS, done) + accretion physical-width
fix → re-run the sweep, expect wedge% to stabilise and full-system r to climb
toward ~0.78. Curtain/contrast is NOT the lever (disproven). (α) full-system
churn is NOT paid (γ already showed the contrast is not the lever).

## Explicitly OUT of scope

- **64² geography calibration (cap / n_cycles, Issue #141)** — NOT
  implicated; geography already converges. This fix does NOT reopen #141.
- **Upscale wiring** — gated behind invariance; not in this chantier.
- **SCULPTING chantier (flat interior + un-reworked init-Voronoi boundary,
  Lecture A)** — a DISTINCT downstream chantier, AFTER invariance. Invariance
  REVEALS the true un-sculpted state (at 64² masked by the grid-width pile
  overflowing inward); the sculpting fix must be PHYSICAL (incising erosion,
  intracontinental rifting, distributed deformation), NOT FBM. Causal order,
  not just sequential.
- **Objective 2 — physical fidelity of `h_critical` (REGISTERED follow-up,
  Phase 3 Lallemand).** `h_critical(d) = h_max·(1−exp(−d/l_taper))` is an
  empirical exponential, NOT the Davis-Suppe-Dahlen critical-wedge mechanics
  (taper = α_surface + β_décollement). Refound it on the true wedge physics
  (a closure relation, not the homemade exponential) at Phase 3 Lallemand,
  where the code already defers fine fidelity. Out of this issue:
  re-founding the formula would change behaviour → reopen the closure
  calibration, which invariance explicitly must NOT do.

## Anti-patterns honoured

Structural-convergence target (not bit-identity); upscale not wired;
invariance vs sculpting kept distinct; **invariance vs physical-fidelity
kept distinct** (Fix #1 is a unit conversion, NOT a formula re-foundation —
a pre-code physics check caught Point b dressing an α_c number into a
non-critical-wedge formula); root cause surfaced BEFORE code; the curtain
coupling left as a post-Fix-#1 measurement rather than a presumption.
