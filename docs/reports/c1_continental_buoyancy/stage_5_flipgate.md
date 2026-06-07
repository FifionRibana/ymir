# Issue #145 — Stage 5: flip-gate (legacy-vs-no-flux visual comparison)

Question (separated from the commit decision): does the no-flux residual (subduction grid-aligned promotion, +32 vs legacy) block the production flip (5d)? Decided by visual comparison to the legacy reference (seed 2026, ×8).

## Comparison (seed 2026, ×8, production path)

| | plate_type ×8 | production altitude ×8 |
|---|---|---|
| **legacy** (rigid=false) | continental blob + left-edge horizontal stipple + lower-right detached blob + sparse bottom (~26 thin features) | continent w/ internal texture (advected crust) + left-edge stipple/streaks |
| **no-flux** (rigid=true) | continent + lower-center diagonal fan + right-edge stipple (faint vertical hint) | flatter continent (rigid/frozen) + faint left line + diagonal fan |

**Verdict: the no-flux boundary stipple is at LEGACY LEVEL.** Different locations, comparable visual prominence; **neither has the prominent vertical bands** (those were the v=0 transport-leak artifact, now fixed). Legacy already carries ~26 tolerated thin features; the no-flux +32 manifests as a diffuse diagonal fan of comparable prominence, not a new alarming structure.

→ **Flip is visually SAFE** (residual at legacy level). The subduction grid-aligned promotion is a registered follow-up, not a flip blocker.

(Side observation: the legacy continent shows internal altitude texture from advected continental crust; the no-flux continent is flatter — rigid/frozen at thickness ~1.0 — which is precisely the intended buoyancy fix.)

## Registered follow-up (NOT #145)

**Subduction grid-aligned promotion at the rigid boundary** — GEOMETRIC artifact (mechanism: subduction floor-trigger fires on the sharp rigid S̃ discontinuity + convergence, following the cell mask → grid-aligned 1px promotions; +32 vs legacy). Physical mechanism (margin accretion), artificial geometry (grid not physical boundary). Fix so accretion follows the physical continent boundary. Minor (boundary stipple, legacy-level), deferred.

## Remaining for #145

- **5b** — re-baseline the 3 wedge-imprint tests on the committed no-flux base (confirm the boundary stipple doesn't touch their wedge regions).
- **5c** — confirm determinism on the rigid production path (run ×2 byte-identical).
- **5d** — flip production (`RunBaseline`, `phase_a_c1`) to rigid; flip-default vs remove-flag. Flip-gate PASSED.
