# F139 — instruments, declared BEFORE the measurements (2026-10-02)

**Everything not stated here is F138's** (`f138_spec_v2/f138_declared.md`, with its two amendments).
- The world: the témoin, ON extended, run_hd tail.
- The 23 descents and their paths, the next base B.
- The outlet area A = max(the accumulation on the first km, the lake's inflow); S_loi = k·A^−0.5.
- The law below = B + k·χ along the path.
- T3's bodies, L_in / L_bed / L_floor, and the hypsometric area grown from the body's lowest cell.
- φ: splitmix64.

## C1 — T3's calibration
For the 22 lakes with a body: T3's area at r = 1 against the lake's final area in the ON world.
- The final area is the cells of `run_hd`'s lake map with its id, times the cell area.
- The ratio and |deviation| > 20 % are given. **If more than a third deviate by > 20 %, it is said first.**

## C2 — T1b, the minimal steepness
- **The required slope** S_req = min over the path's cells i > 0 of (L − law(i)) / s_i (0 if L ≤ law at the col).
- **m_fit** = S_req / S_loi(A).
- If S_req > tan 28°: the remainder beyond the cap = L − tan 28°·X − B, falling at the base.
- The distribution of m_fit: min / p50 / max.

## C3 — the rule "steepen" (m = 10)
- S = max(10·S_loi, S_req), capped at tan 28°. Only the remainder beyond the cap falls.
- Its counts and heights, beside the fixed-slope rule (F138-T1: S = 10·S_loi, the rest falls).

## C4 — T4's variant (m = 10)
- The head fall φ·D_g₀ is taken first: D_g₀ is F138-T4's gorge drop under the fixed-slope rule.
- The gorge then starts at L' = L − φ·D_g₀ and covers the rest.
- "Fit / no room", the required length and the falls are recomputed from L', under both rules (fixed slope and
  steepen).

## R — retreat per lake (14 D8 lakes; a table, no construction)
- **The parameter**: r_lake = min(2, r_world · (A / A_ref)^p), with A_ref = the median outlet area of the 14 D8
  lakes, p ∈ {0, 0.25, 0.5} and r_world ∈ {0, 0.5, 1, 1.5, 2}.
- **The area**: T3's at L(r_lake). Remaining lakes ≥ 1 km² and the total area are given.
- **The emptying r_world** = 2 / (A / A_ref)^p. Beyond 2 the lake is never emptied within the range.
- **The selector age**, by the map (a), piecewise linear: r 0 → ×0.7, r 1 → ×1, r 2 → ×1.4. The emptying age per
  lake, and its spread as a share of the selector's range 0.7–1.4.
- **p's provenance: DECISION.** The link to the celerity is ANCHORED in form (F115: c = K·A^m), not in value.

## C5 — the below-sea lakes: the merged spill
- **L_in_merged** = the median, over the lake's final footprint, of `ocean_flood`'s spill on the breached S1.
  This is F38's level, "fill each enclosed below-sea region as ONE water body to its col".
- Compared with today's level: the count within 1 m, then L(1) = min(L_in_merged, today) against today.
- The 9 spillway lakes (F137-Q4). 1000001 is included, with its footprint.
