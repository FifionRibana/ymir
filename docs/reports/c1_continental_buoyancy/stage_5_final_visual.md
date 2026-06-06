# Issue #145 — Final visual acceptance (complete corrected production state)

Sequencing (per user): the complete production visual MUST be seen BEFORE the flip (5d) and before re-baselining (5b) — each fix was validated locally (streaks, stipple) but the WHOLE corrected render was never looked at = the gap that let filaments through before. Production render (rigid=true, no-flux + erosion-fix committed, closures ON, 300 steps), seeds 42/2/1337/2026, S̃ + land-mask + altitude (`final_visual/`).

## Result

| seed | land mask | verdict |
|---|---|---|
| 42 | compact mass + few islands | ✅ credible |
| 2 | large continent + internal seas | ✅ credible |
| 1337 | single solid continent (largest 0.94) | ✅ credible |
| **2026** | large continent **+ visible 1px vertical finger** | ⚠️ **visible residual artifact** |

**Core #145 goal visually CONFIRMED:** continents are credible MASSES, multi-seed (3/4 clean). The buoyancy fix works.

**BUT the "no residual artifact" criterion FAILS on seed 2026:** the vertical finger (subduction grid-aligned promotion residual) is visible in production (S̃ + land + altitude). The flip-gate (single-seed, "legacy-level" metric) understated it; the complete multi-seed visual catches it — the lesson (visual gates, not metrics).

## Implication

The subduction grid-aligned promotion residual is NOT "invisible/minor" as earlier classified — it is **visible on seed 2026**. Its classification (follow-up vs #145 blocker) is re-opened. **The flip (5d) should not proceed until this is resolved** (addressed, or explicitly accepted), and re-baselining (5b) is paused (don't engrave imprints while a visible residual stands).

## Decision pending

- Address the subduction grid-aligned promotion in #145 (promote from follow-up to blocker), OR
- Accept (1 seed, 1 thin finger, core goal met) and flip with the follow-up registered.
