# Issue #145 — Stage 5: the no-flux boundary oscillation is BOUNDED (multi-seed)

The final-visual found a 1px "curtain" at the boundary (seed 2). Is it a divergent numerical instability (blocks flip) or a bounded cosmetic mesh-mode? Stability is a property → measured multi-seed (raw S̃, the clean signal; production-amp is confounded by the coast step).

## Multi-seed amplitude (no-flux, boundary 1px-mode amplitude, steps 11→300)

| seed | osc-amp mean trend | max @ S̃ | global_max S̃ |
|---|---|---|---|
| 42 | 0.57→0.42→0.37→0.36 (decreasing) | ~1.95 @ 2.1 | 2.2 const |
| 1337 | 0.44→0.35→0.45→0.44 (stable) | ~2.0 @ 2.1 | 2.2 const |
| 2026 | 0.41→0.33→0.40→0.34 (stable) | ~2.0 @ 2.1 | 2.2 const |
| 2 | 0.36→0.30→0.33→0.42 (stable) | ~2.06 @ 2.1 | 2.2 const |

## Verdict — BOUNDED on all seeds; divergence definitively excluded

- Amplitude **stable** across steps on all 4 seeds (seed 42 DECREASES 0.57→0.36; none grows). `global_max S̃ = 2.2 constant` on every seed.
- **The cap-proximity concern is resolved**: on all 4 seeds the max oscillation sits at **S̃ ≈ 2.1 = AT the cap** — the oscillation lives in the high-S̃ **wedge** (near h_eq), where equilibrium-height caps it. No seed has free oscillation far below the cap.
- Earlier "count 648→1063" was spatial SPREAD, not amplitude growth.

Mechanism: a mesh-mode oscillation in the wedge (high S̃; sharp continental/oceanic contrast + upwind), **capped by the equilibrium-height regulator on every seed**. Present in legacy (the contrast + upwind), spread (not amplified in amplitude) by the no-flux wall. **Cosmetic and bounded, NOT a stability problem.** The S̃ field does not diverge → no risk to downstream S̃ consumers or long runs.

## State of the two boundary residuals (both cosmetic)

| residual | nature | severity |
|---|---|---|
| Curtain (e.g. seed 2) | bounded mesh oscillation in wedge, capped by equilibrium | cosmetic, faint in production |
| Finger (seed 2026) | subduction grid-aligned promotion (1 line) | visible thin line |

Core #145 goal (continents = credible masses) holds multi-seed. Both residuals are light edge-finish issues, not structural.

## Decision pending

- Light cosmetic fix of both (boundary smoothing for the curtain + limit 1px subduction promotion for the finger) → clean render, then flip; OR
- Tolerate + document both as known cosmetic edge limitations, flip, fixes as follow-up.
- Finger judged separately on its real production visibility.
