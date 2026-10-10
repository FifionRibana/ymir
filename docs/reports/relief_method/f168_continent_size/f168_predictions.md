# F168 — predictions, written before any measurement (2026-10-10)

**Not blind**: F167's readings, and the code of A1 (in particular, the land is quantised in whole plates out of 8).

## The reviewer's (from the brief)

- **S1**: the initial fraction for 55–60 % final land is within ±0.1 of 0.58.
- **S2**: the mountain share falls from L0 to L2 on ≥ 3 seeds, but stays above 37.2 % on ≥ 2 seeds in L2 (the 187 km
  wedge still covers half the land).
- **S3**: in L2, plain + plateau exceed 20 % on ≥ 3 seeds.
- **S4**: in L2, continent–continent collision cells appear on ≥ 2 seeds.
- **S5**: in L2, Σ Δs_DS > 0.05 covers ≥ 40 % of the land on ≥ 3 seeds.
- **S6**: in L2, the main continent is NO LONGER circumnavigable on ≥ 1 seed (it wraps the torus).
- **S7**: the facets or the walls still fail in the mountain class on ≥ 3 seeds.
- Meta: at least one is false.

## Mine

- **The trials**: each continental plate adds ~10 % of the map.
  - **L1 = f 0.375 (3 plates, 30–40 %)** and **L2 = f 0.625 (5 plates, 50–62 %)**.
  - The response is monotone (the BFS adds plates).
  - **S1 holds** (0.625 within 0.48–0.68).
  - The target is reached on the témoin, though the quantum (~10 %) may put it just outside 55–60 % (then the nearest
    is taken, as declared).
- **S6 holds**: in L2, with 5 of 8 plates continental and contiguous by BFS, the main mass occupies every row or column
  on ≥ 1 seed, or wraps by winding; the circumnavigation is lost there.
- **S4 holds**: with more continental plates, C-C convergent boundaries appear on ≥ 2 seeds. The accretion still merges
  to ~2 plates at the end.
- **S5 holds**: Σ Δs_DS > 0.05 covers 40–75 % of the L2 land (the wedge reaches 30 cells).
- **S2 holds on its two clauses**: the mountain share falls by 10–25 points from L0 to L2, and stays above 37.2 % on
  ≥ 2 seeds.
- **S3 refuted**: plain + plateau stay under 20 % on ≥ 2 seeds in L2 (the uplift covers most of the continent).
- **S7 holds.**
- **The macro success: not met.** The stop that fires is « the mountains above 37.2 % » (and probably the
  circumnavigation lost on one seed too): decompose, report, stop.
- **The rim**: present wherever land touches the window border after the viz's offset. Absent on the L0 seeds the
  offset centres cleanly; present in L2 where the land wraps.
- **The interior**: the share of the land > 40 km from the coast goes from ~0–5 % (L0) to 15–35 % (L2).
- **The cost**: C1 0.2–0.3 s; the physics level 10–15 s CPU; a frozen chain ~200 s CPU.
- **Meta**: at least two of mine are refuted.
