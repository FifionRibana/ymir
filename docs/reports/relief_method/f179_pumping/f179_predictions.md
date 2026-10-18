# Finding 179 — predictions, written before any measurement

## The reviewer's (hypotheses; at least one is false)

- **P1**: on F176's collision belts, the median pumping factor exceeds 2.
- **P2**: the kinematic bound would be exceeded on more than half of the steps of F176's O-C seeds and of F178's Andean seeds.
- **P3**: the bound would cut the median s_src of F176's collision belts by at least 30 %.
- **P4**: with the seed at the arc, the circumnavigation is lost on at most 1 of 15 seeds.
- **P5**: with the seed at the arc, at least 80 % of the seeds lie on cells emerged at t₀.
- **P6**: the seed at the arc does not halve the Andean seeds' pumping relative to F178, so the bound is needed.

## Mine

- **Part 0**: the `f165` repair is test-only (a layout bug), and production reads C1's own s.
- **P1 held**: the collision belts' median factor is 2–5. Many coastal cells have a net change ≤ 0 (their factor is unbounded),
  as F178 measured on the active band (s 0.59 → 0.24 under F176).
- **P2 held, strongly**: the supply per step is at most ~0.5 cell × s_in / 30 (≈ 0.003 in s for an oceanic s_in of 0.2), while
  a boosted O-C increment near its seed is of order 0.05–0.1. So the bound is exceeded on ≥ 80 % of the steps of the O-C and
  Andean families.
- **P3 held narrowly**: for collision s_in ≈ 1, so the supply is ~5 × larger and the un-boosted increments are smaller. The
  bound cuts the belts' median s_src by 20–50 %.
- **P4 refuted narrowly**: depth 4 from the trench often still falls on the shelf or its emerged edge, and the wedge from d = 1
  still pumps there. The circumnavigation is lost on 1–3 seeds, and the land is 52–62 %.
- **P5 refuted**: 40–70 % of the seeds lie on cells emerged at t₀ (the seed reads no sea level, and the band's continental
  shelf is 2–4 cells wide).
- **P6 held**: the Andean family's pumping factor under F179 is ≥ 0.5 × F178's (the pumping is Davis-Suppe refilling C1's
  erosion, wherever the seed sits).
- **The A1 share**: within ± 3 points of F176 with both rules.
