# Issue #145 — Stage 5 BLOCKER: the 1px-band artifact (must fix before flip)

Visual review of the production renders (seed 2026) surfaced two artifacts. 5a (test matrix) was reassuring on scalars/contracts but tests don't see spatial structure — the eye caught what the metrics missed (Nth time).

## Artifact 2 (expected) — sharp continent/ocean altitude edge

The binary-mask discontinuity (`v=0` continental / advected oceanic) → a sharp coast. Known caveat (the 1.99 boundary), subduction softens it in production. Not surprising; not investigated further here.

## Artifact 1 (the blocker) — grid-aligned 1px continental streaks

Long, 1-pixel, axis-aligned streaks emanating from the continent toward a domain edge. Diagnosed (seed 2026, rigid):

| Q | finding |
|---|---|
| Q1 at rigidity-mask boundary? | YES — 69% of raw-S̃ thin-line cells (55/80) adjacent to cont/ocean boundary; streaks emanate from the continent |
| Q2 grid-aligned? | YES — clean vertical 1px lines (H/V) |
| Q3 step 0 or builds? | BUILDS — 1 → 80 cells over 300 steps |
| Q4 raw S̃ or altitude? | BOTH — sharp in raw S̃, blurred-but-visible in production altitude (isostasy gaussian softens; the numerical thin-line detector false-negatived altitude=0 due to the blur — **visual caught it**) |
| Q5 one or many? | MANY — systematic pattern (stippled boundary + multiple streaks) |

**Plate_type cross-ref:** the streaks are `plate_type==Continental` cells carrying LOW S̃ (~0.2) = **subduction-promoted** (Oceanic→Continental sets the label, keeps low S̃), then frozen rigid.

## Attribution — rigidity amplifies ×4.5

| | continental cells | continental-mask thin-line streaks |
|---|---|---|
| rigid=false (legacy) | 1486 | 26 |
| rigid=true | 2248 (+51%) | **118 (×4.5)** |

The streaks pre-exist in legacy (26, a minor pre-existing subduction grid-alignment) but rigidity **amplifies them 4.5×**. Introduced/amplified by the flag → blocks the flip.

## Mechanism — rigid × subduction feedback, rooted in boundary flux leakage

1. In `step_upwind`, a rigid cell's right/top **face flux uses the NEIGHBOUR's velocity** — so `v=0` does NOT make the cell a no-flux wall. A rigid continental cell whose ocean neighbour flows *away* **leaks continental crust** out through that face, along grid rows/columns.
2. Leak → grid-aligned low-S̃ features at the boundary.
3. **Subduction** floor-trigger fires along those lines → promotes Oceanic→Continental in 1px trails.
4. New continental cells freeze (rigid) → persist, seed more leak → compounds (1→80→118).

The first hypothesis (pure advection leakage) was half-right: root = the v=0 boundary flux leakage; the *visible streaks* are subduction amplifying it.

## Fix direction — true no-flux rigid boundary

The current hook (`v=0` in `fill_velocity_field`) is **insufficient at the boundary**: zeroing a cell's own velocity leaves the neighbour-driven face flux. The fix: **zero the flux on any face between a rigid and a mobile cell** — a true no-flux rigid boundary, conservative (the face transports nothing on both sides; mass still cancels). Needs a **rigid-aware masked advection** (the per-cell velocity hook cannot express a face condition).

Approaches: (A) extend shared `step_upwind` with an optional rigid mask (zero rigid-adjacent face flux) — clean generalization (internal no-flux boundaries), but touches `tectonics_v2`; (B) C1-specific masked-upwind wrapper — avoids touching shared code but reinvents (W1 watchpoint). DECISION PENDING.

**Blocks 5b (re-baseline) and 5d (flip)** — do not engrave references or ship to production on an un-fixed artifact.
