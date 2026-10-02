# Specification v2 — the age sets how far the gorge has retreated (Finding 138, NOT BUILT)

**Status**: v2. It replaces v1 (`f137_gorge_age/spec_gorge_age.md`, kept as a draft) and corrects its eight
review points (ADR Finding 137). It is for the author's and the reviewer's review **before** any construction.
No production code exists for it.

Every number carries its provenance (rule 15):
- **MEASURED**: this dossier, bench named;
- **ANCHORED**: a source;
- **PROXY**: a labelled stand-in;
- **DECISION**: the author's, or proposed to the author.

The tables T1–T4 are in `finding_138.md` and `f138_t.txt` (bench `f138_t`).

## 0. The author's criterion and decisions (recorded as given)

**The criterion:**
- « Je veux réduire la taille des lacs avec l'âge, car au départ ils étaient tous remplis à la hauteur
  d'équilibre, ce qui donnait des lacs énormes. » (2026-09-30)
- « Un lac présent est un niveau de base pour les rivières qui s'y jettent ; un lac vidé ne l'est plus. »
  (2026-09-29) « Les lacs des bassins sous la mer sont des lacs présents. » (2026-09-30)
- « Chute ou gorge, les deux sont valides. » « Le cas des chutes lié au type de roche n'est pas encore activable
  […]. Ce sera donc intéressant de faire apparaître des chutes d'eau (et donc de les tagger comme tel) là où
  c'est pertinent. Ça donnera de la diversité au monde. » (2026-10-02)
- « C'est bien d'avoir l'âge qui fixe jusqu'où la gorge a reculé. » (2026-10-02) Young = a full lake and a gorge
  downstream; old = a notched rim and a lowered lake.

**The decisions (2026-10-02):**
- the retreat follows **the existing age selector**;
- the gorge's slope is set **by the river's size**;
- the tagged falls: **« les deux »** (the gorge head and the shortage of room). The head fall is not uniform: its
  per-lake share, zero included, is a PROXY of the missing rock.

**Cost**: a world in minutes.

## 1. The retreat parameter r and its link to age (corrects v1 § 2, points 2 and 4)

**r ∈ [0, 2]** per world, following the existing age selector ("âge k", ×0.7 / ×1 / ×1.4, `workspace.rs:1210`).
**The map age → r is a DECISION**, shown with table T3 and not fixed in advance:
- (a) ×0.7 → 0, ×1 → 1, ×1.4 → 2;
- (b) ×0.7 → 0.5, ×1 → 1, ×1.4 → 1.5.

**What r means for the lake level** (L(r), § 3):
- r = 0: full at the input spill L_in;
- r = 1: at the bed sill L_bed;
- r = 2: down to the floor L_floor, i.e. emptied.

**r = 1 keeps ON's level, not ON's outlet geometry.**
- At r = 1 the lake stands on its bed sill, which the construction already fixes today (F136-K7: on all 14 D8
  lakes, the final sill is a bed cell, untouched by the construction).
- With the invariant (§ 4), the downstream side becomes a gorge hung from L(1), not today's step. So **r = 1 is
  not today's ON state**; it shares only its levels.

## 2. Drained lakes (r → 2) (corrects point 2)

- A lake whose L(r) reaches L_floor is drained, and **is no longer a base** (the author, 2026-09-29).
- **Its upstream is re-based on the level of the notched outlet, L(2) = L_floor, never on χ from the sea.**
  - v1's "re-base like OFF" would recreate F135's kilometre trench (OFF's floor about 1 km under the terrain
    through the bowls).
  - Concretely: the χ walk's lake stop (`valley_construction.rs`, the `lake_base` branch) uses base =
    L_floor + the gorge's floor at the bowl's exit, instead of passing through.
- **Gate** (§ 8, G-drained): no cell upstream of a drained lake lies below L(2) further than the law from that
  base allows.

## 3. The lake level L(r) (unchanged in form, now measured)

- r ∈ [0, 1]: `L(r) = L_in − r·(L_in − L_bed)`;
- r ∈ [1, 2]: `L(r) = L_bed − (r − 1)·(L_bed − L_floor)`.

The form is linear (DECISION).
- L_in: the input body's ring minimum.
- L_bed = min(L_in, today's level).
- L_floor: the body's lowest cell.

All three are MEASURED on the construction's input (S1), F138-T3.

The rim keeps ≥ L(r): § 4.

## 4. Where the invariant acts (corrects point 5)

**On the samples' floor, before the cones.** The construction lays each cell from a polyline sample's cone
(Finding 135-C: the nearest sample lays 98.9 % of the constructed cells). The sample's floor is
`zf = base + age_k·χ` in **`line_samples`** (`valley_construction.rs:864`, the line `let zf = base[c] +
vc.age_k * chi[c];` at :897).
- It is called by `skeleton` (:622) and by `traced_polylines` (:1054, :1179, :1204, :1207).
- `carve_diag`'s cones (`geo`, :1468–1488) and its cross-line minimum (:1598–1624) read that `zf` (`src`,
  :1438–1447).

**The invariant**: for every sample whose cell lies on a present lake's outlet path, from the true col to the
gorge's end (§ 5):

`zf = max(base + age_k·χ, z_gorge(s))`.

- Every cone then hangs from the gorge, so a downstream valley's cone cannot dig beside the corridor below it.
- **The rim**: the input body's outer ring keeps `≥ L(r)`, as a clamp in `carve_diag`'s final `v` (:1626–1650),
  for the ring cells only.

## 5. The gorge's profile by river size (corrects points 3, 6 and 7)

**The gorge's slope**: `S_g(A) = min(m · S_loi(A), tan 28°)`.
- `S_loi(A) = k·A^−0.5` is the construction's own graded law at the outlet area (χ = ∫ (1/A)^0.5 dx, A₀ = 1 km²,
  `valley_construction.rs:573`; k = 0.07183, F121's PROXY).
- **m is a DECISION**, chosen on table T1.
- **The cap is 28°**: the construction's own wall slope (ANCHORED in F121).
- **ANCHORED stays EMPTY** until a source for a gorge's or knickpoint's slope as a function of area is supplied.
- v1's S_g = 14.5° (the delivered steps' median, a circular proxy) is withdrawn.

**What m gives** (T1, 23 descents, MEASURED):

| m | S_g range | fit / no room |
|---|---|---|
| 3 | 0.2–6.6° | 4 / 19 |
| 10 | 0.7–20.9° | 13 / 10 |
| 30 | 2.1–28° | 21 / 2 |

- Today's steps are steeper than any m ≤ 30 on the large outlets: mean 8.5–15°, up to 47° over 100 m (T2). A
  gorge at S_g(A) is gentler and longer.

**The gorge's length** runs from the head to its meeting with the law below (`B + k·χ`, from the next base B
along the path).

**No room** (point 7): when the gorge at S_g(A) does not meet the law before the next base, the remaining drop
falls at that base (§ 6, the shortage fall). At m = 10: 10 of 23, 14–581 m (T1).

**Short slopes** (point 6): v1's "never walls" read mean slopes. The gate (§ 8) is set on the **max slope over
100 m** of every untagged descent.

## 6. The two falls (corrects points 1, 8; the author's « les deux »)

- **The head fall** (at the lip, the gorge's head): a share **φ** of the gorge's drop, concentrated at the lip;
  the gorge below carries the rest at S_g(A).
  - **φ varies per lake, zero included**, as a declared PROXY of the missing lithology. It is not uniform.
  - **The proposed law**: φ = 0 with probability 1/3, otherwise uniform in [0.1, 0.5], drawn from splitmix64 of
    (the world seed, the lake id).
  - Provenance: DECISION (the shape and the bounds are proposed, not anchored), deterministic per world.
- **The shortage fall**: at the next base, the remainder when the gorge has no room (§ 5).
- **The tag**: a fall is tagged if its drop > **H_f** over ≤ 2 cells.
  - `H_f(A) = max(1.5 · S_g(A) · 100 m, 10 m)`. The margin 1.5 covers 2 diagonal cells (138 m), so **no cell of
    a gorge at S_g can pass it** (T4: 0 gorge 2-cell drops over H_f).
  - The 10 m floor is a DECISION (to be shown to the author with the counts).
  - The format is v1's: a per-segment `features` list in `rivers.json`, beside `profile_m`:

    `{ "kind": "waterfall_head" | "waterfall_shortage", "point": i, "x": .., "y": .., "height_m": .., "lake_id": .. }`
  - **Today's steps are an artefact and are not tagged.**
- **What T4 gives** (m = 10, MEASURED on today's levels):
  - 17 head falls > 0 (11 tagged), median 26 m, max 200 m;
  - 10 shortage falls, all tagged, median 129 m, max 581 m.
- **C-3 is EXCLUDED from this version** (point 8): its contrasts are a 6 km grid and perfect circles (F137-Q2),
  and falls placed there would line up.

## 7. The below-sea basins' lakes

- Their outflow is a Spillway (F137-Q4). The rule applies along the spillway path with the same r, L(r), S_g(A)
  and falls.
- The outlet area is the lake's inflow (T1's amended A; their D8 accumulation stops in the basin).
- **T3's instrument does not fit them**: their final level is the below-sea merge's, which can sit above their
  input body's ring (1000011: 459.3 m against L_in 325.8 m), and 1000001 has no input body. Their L(r) needs its
  own definition (the merge's spill) before they are built: **open**.
- The 3 lakes without a spillway segment of their own (F137-Q4) are untreated until traced.

## 8. Instruments and gates (corrects points 3 and 4; rule 15: each gate from an independent control)

**Instruments**:
- the steps (F136-K6b), classed by the **max slope over 100 m** (F138-T2), not by mean slope;
- the canyons;
- θ on the carved links with the construction's mask (F135-M2);
- the lakes' number and area by r (T3);
- rule 14 (the removal in km³, where);
- rule 18 (OFF / ON, the témoin C2 /10 col);
- the guard (field and lakes).

**Gates:**
- **G-inert, "the inertness of the code when off"** (v1's gate 1, renamed). With the gorge switched off, the
  build is bit-identical to today's ON extended (field and lakes hashes). It tests the code's inertness, nothing
  about r = 1's shape.
- **G-levels**: at every r measured, each present lake's level = L(r) ± 1 m.
  - The control is T3's L(r), computed on the input, independently of the built world.
- **G-rim**: no cut > 10 m on any present lake's input ring (F135-K1's instrument) at any r ≤ 1.
  - Its negative control: today's ON fails it (18 cuts > 10 m).
- **G-slope100** (replaces v1's circular gate 3): on every untagged descent, the max slope over 100 m ≤
  max(atan(S_g(A)) × 1.5, 28°).
  - Its negative control: today's ON (2 of 6 > 40°, T2).
- **G-below**: the lake below is unchanged in level and footprint (the lake-by-lake listing of F134, before / after).
- **G-drained**: § 2.
- **G-area**: the total lake area is non-increasing in r on the measured r.
  - The control is T3: 1 308 → 1 224 → 1 146 → 367 → 0 km² on the 14 D8 lakes.
- **G-tag**: no gorge cell passes H_f; every tagged fall lies on a profile break ≤ 2 cells.
- θ on the carved links (construction mask) within OFF's CI (0.493 [0.492, 0.495]).
- Canyons 0.

**The viz**: a gated toggle "Recul de la gorge (âge)", off = today's ON. r is read from the age selector by the
chosen map; the tags are drawn as markers on the Relief layers.

**Binary questions for the author's look** (r at the chosen ages):
- upstream: does the valley meet each lake at its level?
- downstream: does a gorge, not a cliff, leave each lake?
- does each lake shrink as the age grows?
- on Relief → Ombrage / Différence: is the gorge head at the col when young and in the rim when old?
- **is every tagged fall where the profile shows a short break?**
- **is no gorge a wall?**

## 9. Cost

- An outlet path per present lake (~ 30 lakes, ≤ 30 km each, ~ 600 cells), one profile, one `max` per corridor
  sample, one clamp on each ring.
- **Estimated < 2 s per world**, against a world of 215–247 s: an ESTIMATE, to be measured when built.
- No iteration: L_in, L_bed and L_floor are read once on the input.

## 10. Open decisions for the author

1. The map age → r: (a) or (b) (T3).
2. **m** (T1): at m = 10, 10 of 23 descents have no room and fall at the next base. At m = 30, 2 of 23.
3. The φ law (its shape, its share of zeros) and H_f's 10 m floor (T4).
4. The below-sea basins' L(r) (§ 7).
5. Whether the levels at r > 1 may empty lakes at the oldest age (T3: at r = 2 none is left, by definition).
