# Specification — the age sets how far the gorge has retreated (Finding 137, NOT BUILT)

**Status**: written for review by the author and the reviewer **before** any construction. No production code
exists for it. Every number carries its provenance (rule 15): MEASURED (this dossier, with its bench), ANCHORED
(a source), PROXY (a calibrated stand-in, labelled), or DECISION (the author's, or proposed for the author).

## 0. The author's criterion (recorded as given)

- « Je veux réduire la taille des lacs avec l'âge, car au départ ils étaient tous remplis à la hauteur
  d'équilibre, ce qui donnait des lacs énormes. » (2026-09-30)
- « Un lac présent est un niveau de base pour les rivières qui s'y jettent ; un lac vidé ne l'est plus. »
  (2026-09-29) « Les lacs des bassins sous la mer sont des lacs présents. » (2026-09-30)
- « Chute ou gorge, les deux sont valides. » « Le cas des chutes lié au type de roche n'est pas encore
  activable […]. Ce sera donc intéressant de faire apparaître des chutes d'eau (et donc de les tagger comme tel)
  là où c'est pertinent. Ça donnera de la diversité au monde. » (2026-10-02)
- « C'est bien d'avoir l'âge qui fixe jusqu'où la gorge a reculé. » (2026-10-02)
  - The reviewer's reading, accepted by the author: **a young continent has a full lake and a gorge downstream;
    an old continent has a notched rim and a lowered lake.**
- **Cost**: a world is produced in minutes, not hours.

## 1. What the measurements give this design (Finding 137 Q / A, Findings 114–136)

| fact | value | provenance |
|---|---|---|
| A step leaves every present lake whose outlet the construction lays from the base below | ON extended: 8 / 6 / 2 steps > 50 / > 200 / > 500 m (14 lakes with a D8 outlet) | MEASURED, F136-K6b, F137-Q1 |
| The delivered stream-power world has the same object | livré 9 / 6 / 3; A1+B2 (F114's world) 5 / 2 / 1 | MEASURED, F137-Q1 (= F114's "family 1") |
| The OFF construction hides it | it cuts through the cols (F135-K1, OFF 11 cols > 100 m), so it has neither the lakes nor their steps | MEASURED |
| The steps are gorges or steep slopes, never walls | ON extended's six > 200 m: mean 8.5–15.0° (3 gorges < 10°, 3 steep); livré: 3.6–28.6°, median 14.5° (n = 9) | MEASURED, F137-Q3 |
| ON's lake lowering is the input lake's own bed, not the notch | 2.8–47.6 m; the notch is 49–941 m deep, 0.4–17 km from the sill | MEASURED, F136-K7 |
| Nothing protects the rim downstream | the outlet valley is laid from the base BELOW and retreats through the rim: an invariant applied on one side | MEASURED, F136-K7 |
| C-3 (`production_k_field`) is smooth | correlation length 12.5–25 km in x, > 100 km in y; it is a fracture ramp, a coarse rift strip and volcanic discs | MEASURED, F137-Q2 |
| Age today | construction: k (`floor_m`), W(k)'s γ (gated), the light pass's S_eq; delivered: S_eq + k_time. No break position and no retreat state exist | READ, F137-A |
| k is not the age | ×1.71 of k buys ×1.079 of relief | MEASURED, F121-C3 |
| K and the duration are not separately observable | the knickpoint celerity at the shipped K is 180 618 m/yr, Courant 3 699 | MEASURED, F44 / F90 / F115 |

## 2. The retreat parameter

**r ∈ [0, 2]**, one scalar per world. DECISION (the form is proposed; the author decides).
- **r = 0, young**: every present lake stands at its input spill level L_in, full. The outlet leaves the col
  through a gorge whose head is AT the col.
- **r = 1**: the gorge head has retreated through the rim to the lake's bed sill. The lake stands at the bed
  sill L_bed. **This is today's ON extended state** (F136-K7: the lowering is the bed).
- **r ∈ (1, 2]**, beyond: the head cuts into the bed towards the floor, and the lake level falls below L_bed.
  - At r = 2 the lake is emptied.
  - **A drained lake is no longer a base** (the author, 2026-09-29): its upstream rivers are re-based on the
    level below, as OFF lays them today.

**Its link to age.** There is no physical anchor for an absolute retreat distance: K and the duration are not
separately observable (F44 / F90), and the celerity is over Courant (F115).
- r is therefore a **DECISION knob tied to the viz's age selector**, not a derived quantity. The proposed map,
  for the author to set: **×0.7 → r = 0, ×1 → r = 1, ×1.4 → r = 2**.
- An alternative, also proposed: r is its own selector {0, 0.5, 1, 1.5, 2}, independent of k (F121-C3: k is not
  the age).
- The two are not equivalent. The author picks one.

## 3. The rim invariant and the lake level

**The lake level as a function of r** (per lake; L_in, L_bed and L_floor are measured on the construction's
input, F136-K7's instrument):
- r ∈ [0, 1]: `L(r) = L_in − r·(L_in − L_bed)`;
- r ∈ [1, 2]: `L(r) = L_bed − (r − 1)·(L_bed − L_floor)`.

The interpolation is linear: a DECISION, to be shown to the author before any other form.

**The invariant: the downstream side gets the protection the upstream side already has.**
- Upstream, χ stops at the shore (F132 / F133).
- Downstream, every cell of the outlet corridor, from the col to the gorge's foot, may not be laid below the
  gorge profile hung from L(r):

  `floor_outlet(s) = max(law_below(s), z_gorge(s))`,

  where `s` is the distance downstream of the gorge head along the outlet path.
- The rim's other cells (the ring, Finding 135-K1) keep `≥ L(r)`.
- **This closes F136-K7's one-sided invariant.**

## 4. The descent's profile

`z_gorge(s) = L(r) − S_g·s`, until it meets `law_below(s)` (the graded reach below). The **mean slope S_g** has
two candidates, both declared:
- **MEASURED (default proposed)**: the median mean slope of the delivered world's outlet descents (F137-Q1,
  livré, n = 9: 3.6, 10.5, 10.6, 12.2, **14.5**, 15.6, 15.9, 16.8, 28.6°): **14.5°** (S_g ≈ 0.26). It is a PROXY from a world measurement, not a law. The gate (§7) checks
  it against the ON and livré distributions.
- **ANCHORED (if a source is supplied)**: a knickpoint or gorge mean slope from the literature in `docs/refs`.
  None is identified yet, so this stays empty until one is.

The gorge's head position: at r = 0 the head is at the col; for r > 0 it lies at the rim cell where the profile
from L(r) meets the input terrain. The gorge's width is the trunk's W(A), the construction's own.

## 5. Waterfalls

- **Criterion**: on the new descents only, a drop ≥ **H_f** over ≤ **2 cells (~100 m)**. H_f is a **DECISION**
  (proposed 20 m), shown with the count it yields before it is set.
- **The role of C-3**: Q2 shows C-3 is not a local hardness map (12.5–25 km in x, > 100 km in y). It may only
  be used **at its contrasts**:
  - the edges of the volcanic discs (4.29 % of land);
  - the edge of the rift class (4.87 %).

  A fall is proposed where a descent crosses such an edge from low K to high K (a caprock). This is an optional,
  gated criterion; the number of outlet descents crossing an edge is to be measured first.
- **The tag** (format, not built): a per-segment list in `rivers.json`, beside `profile_m`:

  `features: [{ "kind": "waterfall" | "gorge_head", "point": i, "x": .., "y": .., "height_m": .., "lake_id": .. | null }]`

  where `point` indexes `points`.
- **Today's steps are an artefact and are NOT tagged**. Only the new rule's falls are.

## 6. The below-sea basins' lakes

- Their outflow is a Spillway over the col, not a D8 outlet (F136: 12 lakes).
- Measured (F137-Q4): 9 of 12 have a spillway path (3 have no spillway segment: 1000002 / 1000008 at 76.3 m,
  1000010 at 0.4 m).
  - Drops > 50 / > 200 / > 500 m: **3 / 2 / 0**, all gorges (3.9–6.2°).
- **The rule applies on the spillway path**: the col is the spill point, the same r and L(r). The descent
  profile is only applied where the measured descent is a step; for most it is already graded.
- The 3 lakes without a spillway segment are listed, not treated, until their outflow is traced.

## 7. Instruments, gates, the viz

**Instruments** (every one has existed since F133–F137):
- the steps, by F136-K6b, classed wall / steep / gorge (F137-Q3);
- the canyons;
- θ on the carved links **with the construction's mask** (the based set, F135-M2);
- the number and area of the lakes as a function of r;
- rule 14 (the removal, in km³, and where);
- rule 18 (OFF / ON on the same world, témoin C2 /10 col);
- the guard (field and lakes).

**Gates, each set from a negative control first:**
1. **r = 1 reproduces today's ON extended bit for bit** when the gorge profile is switched off. This is the
   control that the new code is inert where it should be.
2. **r = 0**: every present lake's level = L_in ± 1 m; 0 cuts > 10 m at the input cols (F135-K1's instrument);
   canyons 0.
3. Steps: at r = 0 no outlet descent is steeper than S_g + its measured spread. No wall is created.
4. θ on the carved links with the construction's mask stays within OFF's CI (0.493 [0.492, 0.495]).
5. Monotone lakes: area(r) non-increasing in r on every lake (the author's criterion).
6. Cost: § 8.

**The viz**: a gated toggle "Recul de la gorge (âge)", off = today's ON. Its values are the r chosen in § 2
(measured values only, no free slider). **Binary questions for the author's look**, at r = 0, 1, 2:
- upstream: does the valley meet each lake at its level?
- downstream: does a gorge, not a cliff, leave each lake?
- does the lake shrink as r grows?
- on the Relief → Ombrage / Différence layers: is the gorge head at the col (r = 0) and in the rim (r ≥ 1)?
- is a tagged waterfall where the profile shows a short drop?

## 8. Cost (against "minutes")

- The construction stage is 120–133 s today (F133-E); a full world is 215–247 s.
- The proposed addition is per present lake (~ 30): one outlet path (≤ 30 km, ~ 600 cells), one profile, one
  `max` over the corridor's cells, plus the ring.
  - **Estimated cost: < 2 s**, i.e. < 1 % of a world. ESTIMATE, to be measured when built.
- Below-sea spillways: the same per path.
- **No iteration, no fixed point**: L_in, L_bed and L_floor are read once on the construction's input. The
  circularity of "present" (F132) is unchanged.

## 9. Open decisions for the author (nothing is built until they are taken)

1. r's link to age: the age selector map (×0.7 / ×1 / ×1.4 → 0 / 1 / 2) or an independent selector.
2. L(r)'s form (linear proposed).
3. S_g: the measured 14.5°, or a source.
4. H_f and whether C-3's contrasts may place falls.
5. Beyond r = 1: do emptied lakes become OFF-like upstream (the 2026-09-29 rule), and from which r?
