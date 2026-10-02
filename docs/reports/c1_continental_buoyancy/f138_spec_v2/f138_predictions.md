# Finding 138 — predictions, written 2026-10-02 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–137, the reviewer's predictions and the eight
review points, and in particular:
- the lakes' levels and lowerings (F136-K7: 2.8–47.6 m);
- their outlet paths, drops and lengths (F136 / F137-Q1, Q4);
- the input lakes' areas (input lake 1 143.8 km² against final 132.8; input lake 2 57.8 against 35.5);
- **S_loi's form, read from the code before writing this**: χ = ∫ (1/A)^0.5 dx (A₀ = 1 km²,
  `valley_construction.rs:573–574`), so `S_loi(A) = k·A^−0.5` with k = 0.07183. At A = 100 km² that is 0.0072
  (0.41°).
- The outlet areas themselves are NOT known: blind there.

## T1 — the gorge by m
- **P-T1**: at m = 10, S_g = 0.72·A^−0.5. That is ≤ 0.07 (≤ 4°) for the lakes' outlets (A ≥ 100 km²), so the
  gorges must be long: 290–860 m of drop needs 4–20 km.
  - **At m = 10, at least HALF of the 23 descents have no room** before the next base. That is more than the
    reviewer's third: same direction, a stronger claim.
  - At m = 3, almost none fit (≤ 3). At m = 30, most fit (≥ 2 / 3).

## T2 — short slopes
- **P-T2**: the max slope over 100 m on today's steps > 50 m is set by the construction's cones (28° walls) and
  by S1's own terrain, which the breach only lowers.
  - In ON extended, ≥ 2 of the 6 > 200 m exceed 28°, but **at most 1 exceeds 40°**, and none 60°.
  - The delivered world is steeper: its incision cuts, so ≥ 2 of its steps > 50 m exceed 40°.
  - This disagrees with the reviewer on ON's > 40°.

## T3 — lakes by age
- **P-T3**: the lowering from r = 0 to r = 1 is small (3–48 m) and the bowls are wide, so the total lake area at
  r = 0 is **1.1–1.4× that at r = 1**, not more than 1.5×. Against the reviewer.
  - The area falls steeply only beyond r = 1, towards the floors. At r = 2 every measured lake is gone, by
    definition.

## T4 — falls
- **P-T4** (with the reviewer): at m = 10 the space-shortage falls are higher in median than the head falls.
  - The head falls are a share φ ≤ 0.5 of the gorge's drop, with zeros.
  - The shortages are the whole remainder of hundreds of metres.
  - No gorge cell passes H_f, by construction of H_f with its margin (checked, not assumed).

## Meta
- Disagreements with the reviewer: P-T2 (ON's > 40°) and P-T3. **At least one of these two is wrong.**
