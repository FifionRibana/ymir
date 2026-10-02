# F138 — instruments of the tables, declared BEFORE the measurements (2026-10-02)

**World**: the témoin, the viz state C2 /10 col, ON extended, through `run_hd`'s tail. Field = conditioned, in
metres.

**Lakes**: the 14 with a D8 outlet (F136) and the 9 with a spillway path (F137-Q4).
- **The path**:
  - D8 lakes: from the true col down the ON D8;
  - spillway lakes: the spillway segment and its downstream chain, then the D8.
  - It stops at the entry into another lake (**the next base = that lake's level**), the sea (**base 0 m**), the
    edge, or 30 km (base = z at the end, marked "no base").
- **L** = the lake's level today (r = 1).

## T1 — the gorge by m
- **A** = the largest ON accumulation (cells × the cell area) on the path's first 1 km. **S_loi(A) = k·A^−0.5**,
  k = 0.07183 (the construction's own law, `valley_construction.rs:573`).
- **D** = L − B (the drop to the next base). **X** = the path length to the base.
- **The law below**, at a path cell = the extended skeleton's `floor_m` there (the construction's own floor from
  the next base). Where χ is not finite: B + S_loi(A)·(X − s).
- **The gorge**, for m ∈ {3, 10, 30}: S_g = min(m·S_loi(A), tan 28°). z_g(s) = L − S_g·s.
  - **It fits** if z_g(s) ≤ law(s) at some s ≤ X; the required length is the first such s.
  - **Otherwise the shortage fall** = z_g(X) − B, at the base.

## T2 — short slopes
- **The steps** = the drop zones of F136-K6b > 50 m (col → the first graded reach ≥ 2 km), in ON extended and in
  livré.
- **The max slope over 100 m**: max over the zone's cells of (z_i − z_j) / (d_j − d_i), j the first cell ≥ 100 m
  downstream.
- In degrees; counts > 28°, > 40° and > 60°.

## T3 — lakes by age
- **The input body** of each lake is the one with the largest overlap with its final footprint, among:
  - the input lakes (S1, default config);
  - the closed depressions (`basin_base` rule on the breached S1);
  - the inland below-sea components (`water_class` INLAND on S1).
- **L_in** = the lowest cell of its outer ring: S1 for a lake; the breached S1 for a depression or a below-sea
  basin. **L_bed** = min(L_in, today's level). **L_floor** = the body's lowest cell (same field).
- **L(r)** = L_in − r·(L_in − L_bed) for r ≤ 1, then L_bed − (r − 1)·(L_bed − L_floor). Read at
  r ∈ {0, 0.5, 1, 1.5, 2}.
- **The area at L(r)**: the 8-connected cells with z < L(r) grown from the body's lowest cell (the same field),
  times the cell area. Lakes ≥ 1 km² are counted; the total area is summed.
- **The selector age each r would propose**, shown as the two DECISION options:
  - (a) ×0.7 → 0, ×1 → 1, ×1.4 → 2;
  - (b) ×0.7 → 0.5, ×1 → 1, ×1.4 → 1.5.

## T4 — falls (m = 10)
- **φ, the head fall's share** (PROXY of the missing lithology, DECISION form). Per lake, from FNV-1a of
  (seed, lake id): φ = 0 with probability 1/3, otherwise uniform in [0.1, 0.5].
- **The gorge's drop** D_g = the drop from L to the gorge's end (to its meeting with the law if it fits,
  otherwise S_g·X).
- **The head fall** = φ·D_g, at the lip. **The shortage fall** = T1's.
- **H_f(A) = 1.5 × S_g(A) × 100 m, floored at 10 m.** The margin 1.5 covers 2 diagonal cells (138 m). A fall is
  tagged if its height > H_f.
- **Check**: along each lake's path, the designed gorge's drop over any 2 cells (the path's own spacing) stays
  < H_f.
