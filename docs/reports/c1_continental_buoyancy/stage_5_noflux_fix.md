# Issue #145 — Stage 5: no-flux rigid boundary (artifact-1 transport fix)

The 1px-band artifact's root (per `stage_5_artifact_band.md`): per-cell `v=0` is a CELL condition, but conservation at the boundary needs a FACE condition (`step_upwind` computes face flux from the NEIGHBOUR's velocity, so a rigid cell leaks crust through its face to a receding mobile neighbour). Fix: **zero the flux on any face touching a rigid cell** — a true no-flux rigid wall (conservative: the face transports nothing on both sides).

## Implementation

- `tectonics_v2::advection::step_upwind_masked` — copy of `step_upwind` that zeroes flux on faces adjacent to a `rigid[]` cell. Shared `step_upwind` signature untouched (other callers unaffected).
- `time_loop`: rigid path builds `continental_rigid_mask` (once gallery / per-step Track D) and calls `step_upwind_masked` for both S̃ and age. Non-rigid path unchanged → **suite byte-identical for rigid=false** (only pre-existing `rectangular_simulation` + the deferred wedge_p95 differ).

## Validation (prototype, closures ON, seeds 42/1337/2026)

- **Transport leak fixed COMPLETELY** — discriminator: with subduction OFF (no-flux rigid), continental-mask streaks drop 58 → **5** (irreducible minimum). The prominent grid-aligned bands (the flagged artifact) are gone; proven to be the transport-leak component.
- **Continents preserved** (craton 75–82%, largest 0.72–0.92), **mass conserved-ish** (no-flux is conservative), **no wall pile-up** (boundary edge 0.14–0.19; subduction drains — the feared subtlety does not materialise).

## Residual — a DISTINCT, pinned issue (subduction, not transport)

Decomposition (seed 2026, no-flux vs legacy continental-mask thin-lines): legacy 26, no-flux 58, **common only 2, ADDED-by-rigid 56, all at the cont/ocean boundary, orientation V=29/H=27 (zero diagonal)**. Subduction discriminator: subduction ON 58 / OFF 5 → **the residual 56 are subduction promotions** (Oceanic→Continental) along the rigid boundary. `plate_type` changes only via Track D, never advection, so these are not transport.

**Nature: GEOMETRIC ARTIFACT, not richness.** The mechanism (margin accretion) is physical, but the geometry is GRID-ALIGNED (a real margin follows the physical continent boundary, not grid axes). The rigid boundary's sharp S̃ discontinuity (1.0/0.2) + convergence makes subduction's floor-trigger fire in grid-aligned lines along the mask. Rigid adds ~32 vs legacy (26→58). Minor visually (boundary stipple, not the prominent bands).

## Decisions (separated)

1. **Commit the no-flux NOW** — it is the correct, proven-complete transport fix; the subduction residual is distinct and does not block the transport fix.
2. **Does the residual block the FLIP (5d)?** — DIFFERENT question (committing the fix ≠ shipping rigid-by-default). Decided by a **legacy-vs-no-flux visual comparison** (is the boundary stipple at legacy level → flip OK, or visibly worse → treat before flip).

## Registered follow-up

**Subduction grid-aligned promotion at the rigid boundary** — geometric artifact (mechanism pinned: floor-trigger on the sharp rigid S̃ discontinuity + convergence, following the cell mask). Fix so accretion follows the physical boundary, not grid axes. +32 streaks vs legacy under rigid.

## Next

- Legacy-vs-no-flux ×8 visual comparison (seed 2026) → flip-safety gate.
- 5b (re-baseline the 3 wedge-imprint tests) can proceed on the committed no-flux base — confirm the boundary stipple doesn't touch the wedge regions of those tests.
- 5c determinism, 5d flip (gated on the comparison).
