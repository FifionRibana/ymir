# Issue #145 — Closure re-foundation PLAN (point 6 → fix)

Targeted, medium-weight re-foundation of the 3 role-changed mechanisms (point 6 / `stage_point6.md`). **Central constraint: NON-ADDITIVITY.** Add-one DS +237%, erosion +247%, but net ALL +25.7% — equilibrium-height (the regulator) absorbs almost all of it; only the residual leaks. The three mechanisms are **COUPLED through the cap** → they CANNOT be re-calibrated independently and summed. Every measurement is at the **SYSTEM level** (all closures, rigidON), never the isolated mechanism.

## Step 1 — Erosion floor-clamp: MECHANISM fix (not a rate)

- The floor-clamp (`s_new = max(floor, s−δ)`) **injects mass** — a known Phase-1.4 non-conservation. Erosion must NEVER add mass, on any transport. Fix it **because it is wrong** (principled, transport-independent), NOT "because it over-produces" — non-additivity means system mass may not visibly drop (the cap already absorbed the injection).
- **MODEL DECISION (pending):** does erosion *transport* (deposit eroded material downstream → sediment plains, mass-conserving) or *cleanly remove* (no upward injection, but eroded mass still leaves the budget)? → surfaced to the user before coding.
- **Measure AFTER at SYSTEM level** (all closures, rigid): mass + spatial bar + visual. NOT "does erosion still inject" (it won't, by construction) — the effect reads through the coupling.

## Step 2 — Davis-Suppe rate

- Calibration artefact (tuned aggressive for the old transport's frontier-only survival). Reduce the rate.
- **Non-additivity:** halving DS does NOT halve system mass (the cap was already absorbing it). Calibrate against the REAL target (system mass ~+0%, compact continents) WITH the cap present — not against "DS-isolated produces X".
- Analytical first-pass → system measure → visual. **Max 3 iterations** (recursive-tuning = structural limit; do not loop).

## Step 3 — Confirm equilibrium-height (h_eq=2.0) as regulator — DO NOT touch

- After steps 1+2, does it still cap (regulator) or go inactive (over-production gone → system no longer saturates it)?
- Inactive after the two fixes → healthy (the cap was a safety net against over-production, now there's no over-production). Still saturating → the system still over-produces; re-examine.

## Validation criterion (per step AND final)

- **System mass → ~+0%** (conservation: neither destruction nor over-production).
- **Spatial bar maintained** (craton ~82–88%, compact, largest high) — reducing over-production MUST NOT re-fragment (the graded variant did exactly that).
- **VISUAL** (the eye, production, closures ON) at EACH step — mass coming down with fragmented continents = false success.
- **Determinism** preserved.

Method discipline: analytical first-pass → SYSTEM measure (not isolated mechanism) → visual → max 3 iterations → document. (calibration-via-visual-review + recursive-tuning-is-a-structural-limit.) No blind tuning, no looping.
