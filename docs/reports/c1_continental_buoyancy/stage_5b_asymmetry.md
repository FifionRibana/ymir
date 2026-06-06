# Issue #145 — Stage 5b: the wedge-asymmetry change is understood (verdict A)

Re-baselining the 3 imprint tests requires understanding the change first. The alarming one: 1.4 `asymmetry` (near/far wedge S̃) drops 2.12 → 0.91 under rigid. A (advection artifact, OK to lose) vs B (rigidity breaks the orogenic vergence). **Verdict A, proven.**

## Mechanism

- **Davis-Suppe source** (`apply_davis_suppe_step`): drives each wedge cell toward `h_crit(d)`, which **grows from 0 at the boundary toward `h_max` deep in the wedge** → intended profile is **near-THIN / far-THICK** (the critical taper; asymmetry near/far < 1).
- **`asymmetry = bucket_mean[near] / bucket_mean[far]`** (test metric).

## Profile measurement (DS+EH+erosion, seed 42, near/mid/far mean S̃)

| regime | near (0-5) | mid (5-10) | far (10-20) | asymmetry |
|---|---|---|---|---|
| legacy (rigid=false) | 0.236 | 0.176 | 0.112 | **2.10** (monotone ↓ from boundary = toe-pile) |
| rigid (rigid=true) | 0.627 | **0.946** | 0.685 | **0.91** (peak at MID, magnitude 3-5×) |

## Verdict A — proven by three convergent lines

1. **Source intent** is near-thin/far-thick (h_crit grows with distance).
2. **Legacy** is near-thick monotone-decreasing-from-boundary (near 0.236 > far 0.112) = the OPPOSITE of source intent → **advection piling crust at the toe** (boundary). The legacy asymmetry>1 is an advection distortion, not the DS taper.
3. **Rigid** peaks at mid (0.95), 3-5× higher magnitude, source-driven (absolute S̃ = fill × h_crit; fill decreases 0.75→0.36 with distance). No advective toe-pile.

**Losing asymmetry>1 = losing the advection artifact, NOT breaking vergence.** The rigid wedge is MORE faithful to Davis-Suppe (follows the source taper modulated by fill capacity) and stronger (3-5×). The Phase-1.4 `asymmetry > 1` assertion was encoding the advection toe-pile as the expected signature — it tested the artifact.

**B refuted:** rigidity does not break orogenic mechanics; it improves them (source-driven, higher wedge). **Phase 3 will build on STRONGER wedges, not flattened ones.**

## Implication for the re-baseline (now unblocked)

- The rigid wedge is the truer Davis-Suppe wedge → moving the imprint tests to rigid is justified.
- **1.4 `asymmetry > 1`**: re-express. The true signature is the source taper — **fill ratio decreasing with distance** (0.75→0.36) and/or magnitude lift, NOT near/far > 1 (which was the advection toe-pile). Do NOT engrave `asymmetry ≈ 0.91`; assert the source-driven structure.
- **1.4 `wedge_p95`**: rigid 1.640 (attribution: erosion-fix 0.696→0.359 DOWN, rigidity 0.359→1.640 UP). Re-baseline the band with this two-cause attribution documented.
- **1.2 / 1.3 `wedge_p95`**: rigid lifts to ~1.80 (source-driven, no advective dispersion). Re-baseline with the "wedge higher because crust not advected away" justification.

Each new baseline must carry its justification (expected shape + attribution), not just the value.
