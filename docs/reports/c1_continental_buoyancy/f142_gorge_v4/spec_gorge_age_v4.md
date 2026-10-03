# Specification v4 — the age sets how far the gorge has retreated (Finding 142, BUILT GATED, OFF by default)

**v4 = v3** (`f140_gorge_build/spec_gorge_age_v3.md`) **with the author's decisions of 2026-10-02 (recorded in F142),
the minimal plain and the disposition of F141's eight differences.** Everything in v3 not restated here stands.
Provenance tags as in v2 / v3: MEASURED, ANCHORED, PROXY, DECISION.

## 0. The decisions (recorded as given, ADR Finding 142)

- **The scope**: « Seulement les 14 lacs D8 (sous la mer inchangés) » for the first version (DECISION).
- **The drained bowl: the minimal plain.** The author asked « Les deux étant possibles, qu'est-ce qui les discerne ? »
  (plain or valley). The reviewer answered that they are two stages of one history:
  - the plain comes from the sediments laid during the lake's life;
  - the valley comes later, when the outlet goes below the old floor, which the range r ≤ 2 does not reach.
  He proposed a minimal plain (the residual hollows filled to their spill), with the valley and the full sediment
  cover for later. The author: « Ok go » (DECISION).
- **F141's rule**: a body and a lake are paired by footprint overlap, never by index.

## 1. The parameters

- **`GorgeRetreat::v4(r_world, p)`** = v3's values (m = 10, A_ref = 418.7 km², a declared constant) plus
  `d8_scope = true` and `plain = true`. Both fields serialise, so v4 moves the cache key against v3.
- `None` stays off, skipped in the serialised key: the guarded digests are unchanged.
- **The viz** builds v4 (the toggle stays HIDDEN until the gates are reported).

## 2. What the construction does (v3 § 2, with the scope and the plain)

1. **The scope (PROXY for "the 14 D8 lakes", deviation 6).** The construction knows no final lake: F132-P4's
   circularity. So `gorge_bodies` keeps only the bodies that are **input lakes standing alone**:
   - an input lake (not a closed-depression component);
   - whose outer ring touches no other body's closed-depression cell (no merged below-sea water body);
   - and no cell under the sea.
   - The other bodies (the spillway lakes and the below-sea basins' fringes, and an input lake merged with one)
     keep the extended lake base unchanged: their χ stops at L_in, no retreat, no rim clamp, no gorge.
   - **The proxy is checked against the 14 by footprint overlap (G-bodies).**
   - **Measured before building (F142-R, declared here): the proxy keeps 18 bodies.** The 14 D8 lakes match 1:1
     (the body's share 60–98 %). The 4 extras are input lakes standing alone:
     - bodies 7 (5.2 km²), 11 (5.3 km²) and 13 (8.2 km²) have no final lake in ON;
     - body 20 (6.1 km², A 1 108 km²) is 20 % under the below-sea lake 1000016, without touching its depression
       cells.
     - **No construction-side criterion separates them without a threshold of convenience.** They stay in the
       scope, declared, and are listed by G-bodies and read by the other gates.
     - **This deviates from the author's « seulement les 14 lacs D8 » by these 4 small bodies.** Named for the
       author.
2. **v3 § 2 points 1–9** (L_in, L_floor, L_bed, the outflow and A, the lake as a base at L(r_lake), the gorge, the
   invariant, the rim, the falls) apply to the kept bodies.
3. **The minimal plain** (`gorge_plain`, `valley_construction.rs`):
   - **Where it runs**: in every kept body with r_lake > 1 (partly or wholly drained), on the construction's output
     (`carve`) and before the light pass (`production_upscale.rs`, right after `carve`).
   - **The fill**: a priority flood restricted to the input footprint, seeded on its ring (every outside neighbour
     at its height). The outlet is the lake left or the gorge's head.
     - A cell whose fill level is the footprint's spill (its lowest ring cell) is the lake's bed, under water: kept.
     - A cell whose fill level stands higher is a residual hollow: raised to it.
   - **It is the construction's FIRST term that raises the terrain: a deposit.** `carve` still only lowers. The
     volume, the deepest fill and the area are published per body (rule 14).
   - In normalised units (the metres map is affine), so no round-trip micro-pit.
   - The valley through the old floor and the full sediment cover are NOT built (later).

## 3. F141's eight differences: their disposition in v4

| # | the difference | v4 | reason |
|---|---|---|---|
| 1 | A = the col's drained area | **maintained** | On the 14 D8 lakes it agrees with F139's instrument within 0.86–1.24, except lake 19 (0.51). Re-tabulated (F142-R) |
| 2 | L_bed = the construction's minimax on the input | **maintained, re-tabulated: NOT the same object** | F142-R: as built, L_bed is within 0–2 m of L_in, and 2.8–46.6 m ABOVE today's level (the bed sill F136-K7 saw is cut by ON's construction, not present in the input). So r ≤ 1 barely shrinks the lakes (1 308.6 → 1 304.5 km² for the 14, against F138-T3's −12 %). At ×1, p = 0.5 the 14 hold 653.4 km² (F139-R: 521.3). Named for the author; the gates read the re-tabulated tables |
| 3 | φ = splitmix64 of the lowest cell | **maintained** | The construction has no seed and no final id. Re-tabulated (T4) |
| 4 | the 3 below-sea lakes without a spillway | **without object** | Out of the scope: the below-sea lakes are unchanged |
| 5 | off-trunk outlets (A < 10 km²) keep L_in | **maintained** | None among the 14 (A ≥ 156 km²); declared behaviour |
| (a) | a body holds land cells only | **without object** | A lone input lake touches no cell under the sea (the scope tests it) |
| (b) | L_floor read on S1 | **maintained** | For input lakes T3 also read S1: no difference within the scope |
| (c) | a merged lake split into several bodies | **corrected by the scope** | A body whose ring touches another body's depression cell is dropped. Checked by G-bodies (no lake split) |
| 6 (new) | the scope is a proxy for "the 14 D8 lakes" | **declared** | The construction knows no final lake. G-bodies checks it 1:1 by footprint |

## 4. The gates (F142; each with its control, rule 15)

- **The worlds**: the témoin at ×0.7 / ×0.85 / ×1 / ×1.2 / ×1.4 (p = 0.5); ×1 (p = 0, 0.25); ×1.4 (p = 0, the worst
  case: every lake drained). OFF and ON at ×1 for the controls. **Each world runs isolated**: a panic is a failed
  gate, attributed by stage, and the others go on.
- **G-bodies**: each of the 14 has exactly one kept body (1:1 by footprint, the body's share > 50 %); no lake is
  split; the kept bodies that match no D8 lake are listed.
- **G-pits**: in every kept body with r_lake > 1, the cells a footprint-restricted fill would raise by > 0.1 m (the
  closed hollows outside the lake). Measured at each stage up to the breach's input:
  - (a) the construction + the plain;
  - (b) + the light pass;
  - (c) + the droplets;
  - (d) + the bathymetry (= the breach's input).
  - Gate: 0 at every stage; F38's invariant never fires.
- **G-sea**: the land cells ≤ sea after the breach are exactly ON's 292 (all under lake 1000016), at every age.
- **G-deposit**: per body, the deposited volume, the deepest fill and the plain's area. Reported with the gorge's
  removal against ON at the same age (no gate without a control).
- **F140's gates, against the re-tabulated tables** (F142-R): G-inert, G-levels, G-rim, G-ring, G-slope100,
  G-drained, G-area, G-tag.
- **Non-regression**: canyons 0, the coast, θ on the carved links (the construction's mask) within OFF's CI.
- **Rule 14** (removal and deposit, km³, where), **rule 18** (OFF / ON / gorge), **the cost**.
