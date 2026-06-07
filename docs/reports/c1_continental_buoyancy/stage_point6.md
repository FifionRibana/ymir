# Issue #145 — Point 6 result: per-closure mass behaviour on rigid transport

The net +25.7% mass gain (point 2) decomposed by closure, on **rigid transport (rigidON)**, seed 42, 64², 300. Add-one (each closure alone vs advection baseline) + leave-one-out (marginal in ALL).

## Add-one — each closure alone (baseline `none` = −0.00%)

| closure | nature | mass Δ alone |
|---|---|---|
| none (advection) | — | −0.00% |
| **davis_suppe** | source | **+237.25%** |
| **equilibrium** | ceiling (h_eq=2.0) | **−31.68%** |
| **erosion** | sink (?) | **+246.81%** |
| stein-stein | render | −0.00% |
| subduction | sink | −1.40% |
| accretion | source | −0.00% |
| rifting | — | −1.61% |
| ALL (point-2 ref) | — | +25.69% |

## Leave-one-out from ALL (marginal)

| remove | mass Δ | marginal |
|---|---|---|
| ds | +7.24% | DS ≈ +18% |
| **eh** | **+201.01%** | **equilibrium ≈ −175% (the regulator)** |
| er | +21.22% | erosion ≈ +4.5% |
| sub | +25.22% | ≈ 0 |
| acc | +21.29% | ≈ +4% |
| rif | +25.77% | ≈ 0 |

## Verdict — mechanisms change ROLE (not a global −20%, not diffuse)

This is the user's **third** case (role change), precisely localised:

1. **equilibrium-height (h_eq=2.0) becomes the LOAD-BEARING regulator.** Not robust/neutral as one might assume — it is doing enormous work: alone it's a −32% sink, and removing it from ALL explodes the system to **+201%**. On the old (broken) transport it rarely triggered (crust was evacuated below h_eq); on rigid transport, preserved+piling crust constantly hits the cap. It is now the single mechanism preventing runaway.

2. **Davis-Suppe massively over-produces (+237% alone).** Consistent with the point-6 hypothesis: tuned aggressive *because* convergent frontiers were the only place crust survived the old evacuation. On rigid transport (interiors preserved) it over-piles. Marginal +18% in ALL (the rest absorbed by the equilibrium cap).

3. **Erosion's floor-clamp flips from sink to mass-INJECTOR (+247% alone).** The known Phase-1.4 floor-clamp non-conservation (`max(floor, …)` with no deposition raises sub-floor cells), negligible on the old transport, becomes dominant when interiors/margins are preserved. This is a **mechanism fix** (the non-conservation), not merely a rate.

4. **subduction / accretion / rifting / stein-stein:** marginal (≤ ±4%). Not re-calibration targets.

The net +25.7% is the **residual leaking past the equilibrium cap** after it absorbs the DS + erosion-floor over-production. The closures interact strongly (non-additive: DS+247, erosion+247 alone, but ALL only +25 because equilibrium caps them).

## seed-2 n_comp=43 islands — explained

Emergent (S̃>0.6) at seed 2: **1296 continental (94%) / 87 oceanic (6%)**. The specks are **continental** (the dominant mass + small continental fragments/peninsulas), **NOT** Davis-Suppe oceanic-frontier peaks. Not a spurious-island artefact; `largest=0.93` confirms one dominant continental mass. No action.

## Re-foundation scope (refines the conditional roadmap)

**TARGETED, medium-weight** — three mechanisms, not the whole Phase 1.x/2 set:
- **Davis-Suppe** source rate — reduce (tuned for broken transport; over-produces on rigid).
- **Erosion floor-clamp** — fix the mass-injection non-conservation (mechanism, not rate).
- **Equilibrium-height** — confirm h_eq=2.0 is the right regulator on rigid transport (it now bears the load).

Subduction/accretion/rifting calibrations look robust to the transport change (marginal mass effect). This is lighter than a full re-tune of every closure, but heavier than a single rate tweak — and erosion is a mechanism fix, not a knob.
