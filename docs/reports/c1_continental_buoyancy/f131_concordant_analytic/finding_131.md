# Finding 131 — the concordant confluence keeps the law, but not canyon 14, and it flattens the mouths and the sources; LTD beats D8 on every analytic path but makes 894 false confluences on a cone where D8 makes none; the lakes' flats set χ for every catchment above them

**Nothing is promoted; no clause is activated.**
- The benches are in `crates/ymir-core/tests/f126_coast.rs`, Finding 131's section.
- The raw outputs are beside this report.
- mur ↔ mer is ON in every world. bruit ↔ pied is ON where there is a wall profile (B2 + B1 + foot).
- The hydrology is D8 everywhere.
- No number of Orlandini et al. is used as a reference.

**Units.**
- 1 cell = 48.8 m (domain), 8192² over 400 km.
- W is the width of the larger line at the junction (2 × its half-width).
- B works in cell units.

**Reading declaration, unfavourable.** My predictions (`f131_predictions.md`) were written before any grep,
reading or measurement. They are non-blind on Findings 73–130 and on the reviewer's predictions. Every
instrument was declared before its run:
- B's surfaces and metrics before any B run;
- B's gate after the D8 references and before LTD (rule 15);
- A's and P's instruments before their runs;
- P's "where" after P's counts, and before its own run.

## 0 — Finding 130 committed; the guard; the greps

- **Finding 130 was committed** as `595e7b6`. The guard through `run_hd` reads **6 / 6 "= banc"**
  (`guard_six_states.txt`), run with this round's gated additions compiled in.
- **This round's library code** is gated, and bit-neutral when off (see "State"):
  - `ValleyConstruction::trunk_band_concordant` (A);
  - `skeleton_patched` / `SkeletonPatch`, a diagnostic hook on the skeleton's pointers (P).
- **A permanent test with a negative control**:
  `the_concordant_clause_lays_the_tributary_at_the_junction_floor`.
  - Control: without a clause some of these cells dam the trunk, and the full clause lays some above the
    junction's floor.
  - Concordant: every one of them is at the junction's floor, and every other cell is bit-identical to no
    clause.
- **The greps (rule 11/11b/11c, `grep_rule11.txt`)**:
  - **NOTHING FOUND** in the ADR for `Playfair`, `accordant`, `concordant`, `replat`, specific catchment
    area, `r/2` and `Costa-Cabral`; in the code, for Playfair, accordant and concordant.
  - "junction's floor" appears only in Finding 130, as the remedy it named.
  - `confluence` first appears in Finding 3's table.
  - `analytic` first appears in Finding 42 ("analytic continent").
  - The Finding 111 / 113 / 121 hits were read at their first occurrence.
  - `basin_base` is true by default (`ValleyConstruction::new`), which matters for P.

## B — LTD on three analytic surfaces (references we compute)

**The surfaces, declared.**
- **The plane**: 60², z = −(y + x/4), open downslope.
- **The cone**: r ≤ 48, apex up, z = −r. The exact specific area is a = r/2, and there are no confluences.
- **The curve**: z = 0.02 x² + y on 60², outlet at y < 0, closed elsewhere. The exact lines are
  x = x₀ e^{2a(y−y₀)/b}, and the exact a comes in closed form from them.
- **The curve's closed form was checked against a narrow traced flow tube**: to 0.005 at five points.
  My first tube check compared its two lines at different y and disagreed by up to 45 %. It was fixed
  before any LTD run, and the closed form never changed.
- **The metrics**:
  - a = A / (cell width), with the exact a at the cell's downslope edge;
  - the path angle is the chord to the k-th downstream cell against the exact chord at the same progress
    coordinate;
  - false confluences are domain cells with ≥ 2 donors;
  - the cone's asymmetry is the CV of the outer ring's area over 16 sectors.

**The references (D8), then the gate (declared before LTD ran).** B passes iff, on each surface, LTD's
mean angle is smaller than D8's at chords 16 and 32, and LTD makes no more false confluences on the cone
than D8 does (0).

| surface | method | specific area MAE / RMSE | path angle, chords 1 / 8 / 16 / 32 | cone: false confluences · sector CV |
|---|---|---|---|---|
| plane | D8 | 0.354 / 1.800 | 14.04 / 14.04 / **14.04** / **14.04°** | — |
| plane | **LTD** | **0.097 / 0.224** | 18.14 / 0.28 / **0.13** / **0.05°** | — |
| cone | D8 | 0.443 / 0.525 | 11.25 / 11.29 / **11.33** / **11.29°** | **0** · 0.215 |
| cone | **LTD** | **0.720 / 0.809** | 15.94 / 3.47 / **2.49** / **3.02°** | **894** · 0.191 |
| curve | D8 | 0.631 / 1.858 | 9.53 / 7.87 / **7.07** / **7.40°** | — |
| curve | **LTD** | **0.336 / 0.440** | 11.87 / 2.41 / **1.36** / **0.84°** | — |

The exact ring's own CV is 0.024.

- **The angles pass on all three surfaces**, by 5× to 280× at the long chords.
- **B FAILS on the cone's confluences: 894 against D8's 0.** That is Paik 2008's critique, confirmed.
  - On the divergent surface LTD's path memory folds neighbouring rays into each other.
  - Its specific-area errors are **worse** than D8's there (MAE 0.720 against 0.443).
  - The sector asymmetry barely improves (0.191 against 0.215; the exact ring gives 0.024).
- **At chord 1 LTD is worse than D8 everywhere** (it zigzags between two pointers), which is by design.
- **B-plan non-regression**: Finding 128's C0 bench, re-run unchanged, reproduces its run 2 to the digit
  (`b_plane_nonregression_f128_c0.txt`). It also shows LTD confluences on planes (18 / 240 / 74), where D8
  has none.
- **⇒ By the round's rule, C does not run.**

## A — the concordant confluence (Playfair 1802)

**The clause.** Finding 130-A0's domain is the cells on a larger line's floor band, upstream of their owner's
junction with it. They are laid at the larger line's floor **at the junction**; nothing else changes.
- The clause changes 204 781 cells at 30 768 junctions, and **0 cells outside them**. A0's 1 252 were only
  the raised ones.

**On the construction** (`a_law_replat_profiles.txt`, carve on S1):

| | B2 → A_c | B2 (A_c) + B1 + foot |
|---|---|---|
| raised > 1 m: no clause / full / **concordant** | 86 / 1 295 / **50** | 66 / 1 278 / **32** |
| θ laid, no clause [CI] | 0.493 [0.491, 0.496] | 0.492 [0.489, 0.495] |
| θ laid, full | 0.469 [0.464, 0.473] | 0.468 [0.464, 0.472] |
| **θ laid, concordant** | **0.496** [0.493, 0.499] | **0.494** [0.491, 0.498] |
| **replat**: junctions with ≥ 1 flat clause cell, concordant / the same cells without the clause | **58.0 %** / 17.2 % | **53.0 %** / 13.2 % |
| replat length, p50 / p90 / max | 1.20 / 9.54 / 93.5 W | 1.30 / 10.97 / 93.5 W |
| tributaries' floor over the last 2 W, p50: concordant / none | 0.1075 / 0.1186 | 0.0872 / 0.0965 |
| … share ≤ 10⁻⁴: concordant / none | **9.2 %** / 1.2 % | **6.9 %** / 1.4 % |
| **source cones touched**: head → 10th-cell slope, mean / p50 | 1 219 lines: 0.0871 → **0.0071** / 0.0911 → **0.0000** | 1 273: 0.0829 → **0.0116** / 0.0856 → **0.0000** |
| source cones reversed / halved by the clause | 140 / 733 | 147 / 755 |

**On the pipeline's worlds** (`a_worlds.txt`, light pass + A1+B2):

| world | canyons | col 7 | col 14 | 13's dam | canyon 13 | lakes | teeth | R8 terrain / R8 network within 2 W of the junctions |
|---|---|---|---|---|---|---|---|---|
| **témoin C2/10 col** | 0 | 623.6 | 119.6 | 299.7 | — | 25 | 1 440 | 0.0945 / 0.3201 |
| B2 → A_c, no clause | 4 | 787.6 | 493.7 | 546.2 | dammed | 33 | 73 | 0.0308 / 0.2288 |
| B2 → A_c, full | 0 | 714.2 | 231.6 | 362.8 | gone | 25 | 22 | 0.0355 / 0.2301 |
| **B2 → A_c, concordant** | **1** (14: 25.9 km²) | 698.3 | 399.1 | 330.7 | **gone** | 27 | 48 | 0.0312 / **0.2070** |
| B2 + B1 + foot, no clause | 3 | 772.8 | 500.4 | 531.8 | dammed | 35 | 1 084 | 0.0292 / 0.2227 |
| B2 + B1 + foot, full | 0 | 704.7 | 250.7 | 351.6 | gone | 26 | 538 | 0.0370 / 0.2225 |
| **B2 + B1 + foot, concordant** | **1** (14: 38.3 km²) | 689.7 | 436.2 | 320.1 | **gone** | 28 | 810 | **0.0362** / 0.1909 |

Heights are in m. The témoin's R8 values read on each family's union of junction disks (0.0941 / 0.3198 on
B1's).

**The round's reading, applied criterion by criterion (no verdict forced).**

| criterion | B2 → A_c | B2 + B1 + foot |
|---|---|---|
| raised at the floor (86 / 66 ± 20) | **NO by the letter**: 50, i.e. *below* the floor (fewer cells above their law than without any clause) | **NO by the letter**: 32, below |
| θ laid inside the no-clause CI | **YES** (0.496 in [0.491, 0.496]) | **YES** (0.494 in [0.489, 0.495]) |
| canyons 0 | **NO**: 1 (canyon 14) | **NO**: 1 (canyon 14) |
| the replat does not rise | **NO**: 58.0 % against 17.2 % | **NO**: 53.0 % against 13.2 % |
| R8 around the junctions does not rise | terrain **+1.3 %** (0.0312 against 0.0308), network −9.5 % | terrain **+24 %** (0.0362 against 0.0292), network −14 % |

- **The clause is not adoptable by the round's own conditions.** It passes one criterion (θ), fails four,
  and its "raised" count misses the band on the favourable side.
- **What the clause does well.** It keeps the law (θ back inside the CI, where the restricted form could
  not), and it removes canyons 3, 7 and 13 (13's dam falls 546 → 331 m).
- **Where it fails.**
  - **Canyon 14 remains** (its col 399 m against the full clause's 232 m).
  - **The mouths flatten**: 58 % of the junctions get a replat, the p50 1.2 W long.
  - **So do the sources: this is Finding 121's trap.** 1 219 lines whose heads lie in a larger line's band
    are laid, head included, at the junction's floor, and their source slope falls to 0 at the median.
- The network's R8 around the junctions FALLS (−9.5 %, −14 %) while the terrain's rises on B1. The
  reviewer's prediction of a network rise is refuted.
- Named, not built:
  - bound the clause to the cells whose own law is above the junction's floor by less than the dam height;
  - or let the tributary's own law resume beyond ~W/2 from the junction, instead of laying its whole run in
    the band.

## P — do the lakes' flats touch what is constructed?

**The test**: the pointers of the lakes' 684 790 flat cells are replaced, the skeleton rebuilt through
`skeleton_patched`, and the construction compared outside the lake footprints (`p_lakes_flats_patch.txt`).
- The identity patch is the hook's negative control.
- LTD's pointers change 528 863 of those cells; the arbitrary choice changes 505 065.

| world | identity (control) | LTD on the lakes' flats | arbitrary |
|---|---|---|---|
| témoin C2/10 col | **0 cells** outside, 0 inside | **1 173 648** outside (291 694 inside) | 1 164 314 |
| B2 → A_c | **0** | **2 692 389** | 2 681 363 |

- **The flats are NOT without object.** Replacing the pointers inside the lakes moves a million cells of the
  construction outside them. The outside changes lie p50 6.4 km, p90 12.8 km and up to 20 km from a lake.
- **Where (`p_where.txt`, témoin, LTD patch):**
  - **98.2 % of the outside changes lie UPSTREAM of a lake, with a changed skeleton χ**;
  - the base changed nowhere; the drained area changed on 20 cells;
  - |Δz| is p50 6.3 m, p90 14.8 m, max 804 m.
- **The mechanism.** The skeleton's χ is integrated down the D8 path through each lake, so the path's length
  across the flat sets χ, and hence the floor, of every cell above it.
  - The breached lakes are not closed depressions on the breached field, so `basin_base` does not restart χ
    at their col.
  - The flats act through χ, on whole catchments, not near the outlets.

## C — NOT RUN

B failed on the cone's false confluences. Nothing of C was measured.

## The clause table, updated

| unwanted effect | Finding 130's status | what Finding 131 measured | **Finding 131's status** |
|---|---|---|---|
| **B1 spurs, Finding 38** | PROVED, ON | ON everywhere; the base rate held; F38 holds | **PROVED, kept ON** |
| **B1 teeth** | unchanged | concordant: 810 against 1 084 without a clause | unchanged |
| **B2 canyons: the full clause** | the deformation located (A0) | reference | deforms θ (0.469), gated and off |
| **the restricted clause** | does nothing useful | — | gated and off |
| **the concordant clause (A)** | named | θ 0.496 in the CI; canyons 1 / 1 (14 remains); replat 58 % / 53 % of the junctions; source cones flattened (p50 slope → 0); R8 terrain around the junctions +1 % / +24 % | **NOT adoptable** (1 criterion passed of 5). Gated and off |
| **LTD for the skeleton (B)** | the reference not reproducible from [O03] | the angles beat D8 on all three analytic surfaces; **894 false confluences on the cone against 0**; cone MAE worse | **B FAILS on the cone; C not run** |
| **the skeleton's flats (P)** | on the lakes; P1 cardinal | they set χ, so the floor, of 98.2 % of the cells they move, upstream of the lakes: 1.2–2.7 M cells | **they matter for the construction, through χ.** Named, not treated |
| **the construction's uncarved third** | measured | — | unchanged (named; the amplitude-linked k named, not measured) |

## Predictions scored

**Mine** (`f131_predictions.md`):
- **P0 HELD.**
- **A**
  - **P-A1 REFUTED**:
    - raised 86 ± 10: refuted (50, 32);
    - canyons 0 / 0: refuted (1 / 1);
    - the cols within ±15 m of the full clause's: refuted (698 / 399 against 714 / 232);
    - "not back to the témoin's": held.
  - **P-A2 REFUTED**: θ is inside the CI (0.496).
  - **P-A3 HALF**:
    - a replat at 30–60 % of the junctions: HELD (58 %, 53 %);
    - its length 0.2–0.5 W: refuted (1.2 W);
    - the terrain R8 +5–15 %: refuted (+1.3 % and +24 %);
    - the network R8 up: refuted (down).
  - **P-A4 REFUTED**: the source cones are flattened, and the teeth fall 34 % / 25 %.
- **B**
  - **P-B1 HELD** (the plane ≤ 1°, D8 14°).
  - **P-B2 HALF**: D8 0 confluences, LTD some, and the gate failing: held; LTD's MAE 0.6–1.0 of D8's:
    refuted (1.6×).
  - **P-B3 HALF**: LTD ≤ 5° held (1.36°); D8 8–20° refuted (7.07°).
  - **The gate failing: HELD.**
- **P**
  - **P-P HALF**: > 0 cells outside and the témoin's hash changing: held; "within ~2 W of the outlets":
    refuted (upstream, km away).
  - **P-where HELD** (98.2 %).
- **C: not run**, as predicted.
- **Meta ("at least two wrong") HELD.**

**The reviewer's:**
- **A**
  - "canyons 0 in both worlds": REFUTED;
  - "raised at the floor": refuted by the letter (below it);
  - "θ inside the CI": HELD;
  - "cols 7 and 14 at most 10 m above the témoin": REFUTED (+75 m, +280 m);
  - "a replat at 20–40 % of the junctions": REFUTED in degree (58 %; a replat, held);
  - "the network R8 around the junctions +10–20 %": REFUTED (−9.5 %, −14 %);
  - "passes its three criteria, fails on the replat": REFUTED (it fails canyons, replat, R8 terrain and
    the letter of "raised").
- **B**
  - "B-plan passes": HELD;
  - "B-cone: MAE at most a third of D8's": REFUTED (1.6×);
  - "fewer false confluences than D8": REFUTED (894 against 0);
  - "a measurable asymmetry": HELD;
  - "B-curve: LTD ≤ 5°": HELD; "D8 ≥ 15°": REFUTED (7.07°).
- **P** ("zero cells change outside the lakes"): REFUTED (1.2 M).
- **C**: not run.
- **Meta HELD.**

## Limitations, stated

1. **A's "junction" and W** come from `Skeleton::line_parent` and the larger line's half-width at the
   junction sample. The replat's "flat" is "no strictly lower 8-neighbour", the D8 tie surface.
2. **The terrain R8 around the junctions** uses 8 × 8 windows, since the 16 × 16 windows do not fit the 2 W
   disks. It is not comparable in level with the global R8; it is compared only between worlds on the same
   union of disks.
3. **The θ intervals are i.i.d. link bootstraps**, a lower bound (spatially correlated links).
4. **B's cone false confluences** count every domain cell with ≥ 2 donors, the apex's neighbourhood included,
   for both methods alike.
5. **P's patches are pointer replacements on a closed set**: the lakes' flats. The effect is measured on the
   construction (carve on S1), not on the pipeline's delivered worlds.
6. **The end marker of `f131_a_worlds`'s log reads "Finding 130-A2 (worlds)"**, a label copied from the
   generator. It was fixed in the bench after the run, and the numbers are the run's.

## State

**Uncommitted**, awaiting the next round's go-ahead.
- `valley_construction.rs`: `trunk_band_concordant`, `skeleton_patched` / `SkeletonPatch`, and the new
  permanent test.
- `cached_product.rs`: the key test's list.
- `f126_coast.rs`: Finding 131's benches.
- This report, and ADR Finding 131.

**Checks after the last change:**
- The guard reads 6 / 6 "= banc" (`guard_after_this_round.txt`).
- The definition's fresh build hashes to `a8d2d538d692c2f0`, the reference (`identity.txt`).
- `cargo test -p ymir-core --release --lib`: 588 passed, 0 failed.
- `cargo check --workspace --tests --release`: clean.
