# Finding 124 — the canonical framing, the col base by default, and an export that told four lies

Draft for the ADR. Seed 1, 8192², **canonical framing x = 0.09375** (the viz's: "le monde de l'auteur
est l'autorité"). All five parts of the round. Promoted: the EXPORT corrections of part 3 only (the
round's "aucune promotion hors des corrections d'export"); B1, B2 and W(k) are gated, off.
Outputs: `../f124_frame/stages.txt`, `../f124_frame/replay_canonical.txt`, and the `*.txt` beside this
file. Predictions were dated in the scratchpad before each measurement; the scores are in § 5.

**Units.** Every hydrology number says **domain** (what the map draws, 400 km across) or
**signified** (the geographic scale the microscope shows: lengths × 7.5, areas and discharges ×
56.25). The ADR entries (Findings 124-P1 to P5) are the record; this report keeps the tables.

## 1 · The framing (part 1)

**The benches now run at the viz's framing.** `tests/common::CANONICAL_ORIGIN = [0.09375, 0.578125]`
is the default `sample_origin`; Findings 107–123's framing survives as `FINDINGS_107_123_ORIGIN`. The
guard reference (`data/bench_field_hashes.json`) was regenerated at the canonical framing for the six
states the viz can show (livré, A1+B2, C1 nue, C2 /10, C2 /3, C2 /10 base mer (F121)), so the badge
reads "= banc" on every one of them.

**The first non-periodic stage, named** (`f124_frame`, the delivered field built stage by stage at
both framings, the benches' rolled by 768 cells):

| stage | bit-identical | columns ≠ |
|---|---|---|
| S1 bicubic only | 90.669 % | 768 |
| S2 + FBM | 90.669 % | 768 |
| S3–S6 + craters, incision, lithology, fracture | 90.666 % | 784 |
| S7 + bathymetry (= delivered) | 22.304 % | 8 192 |

ADR line: **the product is invariant in resolution (Finding 57) but not in framing (Finding 123),
because the upscale samples its coast warp and its FBM at the UNWRAPPED coordinate
`sx = origin + i·scale`, and the bathymetry normalises by the ocean's mean depth, which is global.**
The first gives the 768 columns past the seam at S1; the second makes every ocean cell differ at S7.
Land cells outside the 768-column strip and the craters are bit-identical.

**Finding 108's window holds at the canonical framing**: delivered Δ excavation / breach 12 / 4;
A1+B2 at s_eq 0.021 / 0.024 / 0.027 → 0 / 3, 0 / 3, 0 / 2.

**The "654 km" entry** still needs the author's resolution and mouth: it is not found at either
framing.

## 2 · The col base by default (part 2)

`ValleyConstruction::new` IS the definition: χ restarts at the col of every closed depression
(`basin_base = true`, the author's decision D). `f121()` keeps Finding 121's sea base, shown in the viz
as "C2 /10 base mer (F121)". k recalibrated at the canonical framing: **0.07183** (0.07186 before).

**Permanent test** (`f122_brush::f124_lake_1000011_is_a_bowl`): lake 1000011's counterpart is
494.3 km², level 459.3 m, **D_L 2.29, 0 constructed floor cells below its col**.

**Finding 121/122's table replayed at the canonical framing**:

| world | coast | Δ e/b | lakes | D_L > 5 | lake % | relief p50 |
|---|---|---|---|---|---|---|
| C1 bare | +5 | 0/0 | 17 | 0 % | 12.33 | 565.8 m |
| C2 /10 | +3 | 0/0 | 18 | 0 % | 12.93 | 515.4 m |
| C2 /3 | +3 | 0/0 | 18 | 6 % | 13.64 | 467.2 m |
| C2 /10 sea base (F121) | +3 | 0/0 | 16 | 25 % | 15.97 | 456.2 m |

C3 (k × 0.7 / × 1.4): relief p50 501.0 / 533.3 m; the k × 0.7 world no longer panics. The skeleton
controls pass for every world (PRE → C p90 ≤ 5 cells, 149–169 of 184 mouths matched).

## 3 · The export tells the truth (part 3)

### 3.1 · What the shipped export said, C2 /10

| | shipped | what it is |
|---|---|---|
| river cells in more than one object | **550 of 383 516** (19 pairs) | mostly a Watercourse × Spillway pair to one mouth (297, 172, 68 cells) |
| Σ catchments / land | watercourses 0.63 + spillways 0.82 = **1.45** | biggest: a spillway, 559 989 km² signified = 9 955 km² (domain) |
| longest trunks | **1 173 / 1 093 km** signified = 156 / 146 km (domain) | under the 566 km (domain) diagonal, but a **10–15 km straight line**, sinuosity 10–14 |
| the rectangle | a block of 2 459 river cells, bbox (4818,2979)–(4876,3026), fill 87 % | § 3.4 |

**Units.** The viz shows SIGNIFIED quantities: lengths × 7.5, areas × 56.25. 1 093 km signified is
146 km (domain). That alone does not make any number true or false; everything below is stated in both.

### 3.2 · Four defects in the network itself (rivers.json and everything built on it)

Measured by `f124_topology` (census) and `f124_heads` / `f124_block` (dissection), fixed, re-measured.
**None changes the terrain**: the terrain's consumers of the drainage (`flint_intercept`, the valley
skeleton) read `flow.accumulation` / `flow.direction` / the detection's `lake_map`, which none of the
fixes touches. The guard reference stays valid.

**T2 — every reach read its RECEIVER's values.** `trace_segment` ends a tributary ON the confluence
cell (`flow.rs`, "Include the junction cell as last point"), and the discharge, the runoff area, the
width and the geometric catchment were read at the last point or as a max over the points. The
confluence's accumulation is the union of every branch meeting there. Before: **57 % (C2 /10) / 73 %
(delivered)** of the reaches with a downstream ended on it; on those, exported area / own area **p50
1.38, p90 339×**, the width > 2× its own on **48 % / 47 %** — the delivered world, which has no comb,
as much as C2 /10. Fix: `drainage::own_end` (one point back when the last point lies on the receiver),
used by `c1_drainage` and by the lake clip. After: **p50 = p90 = 1.00** for the area and the width.

**T1 — the lake clip linked runs the water never reaches.** A head run listed as upstream the tail of
every parent, including a tail dying on the shore of a basin whose outlet run the clip drops; the
downstream link went to the receiver's HEAD run even when the join lay beyond a lake. Before: **791 /
14 770 upstream links asymmetric (5.4 %) on C2 /10, 255 / 5 072 (5.0 %) on delivered**, all naming a
reach whose downstream is None; 22 / 2 downstream links entering the named reach nowhere. Fix: the join
index on the receiver's points, the link to the RUN holding it, the upstream the exact inverse. After:
**0 / 0**, both worlds.

**T3 — `trace_segment` stole junction cells.** Every tributary ending on a junction overwrote
`cell_to_segment` there, taking it from the reach starting there or running through it; the links read
`cell_to_segment` at each reach's receiver, so a tributary was linked to ANOTHER TRIBUTARY and the true
receiver lost it. Found on the delivered world's Finding 93 residual (below): a reach starting with
**1 369 km² (domain) and no upstream**, its two donors ending on its first cell but linked to a third
reach. Fix: claim a junction only if no reach holds it. After: trunk heads with > 10 km² at their first
point **37 → 0**; stolen links 0.

**T3b — a reach ending on a cell its receiver RUNS THROUGH was linked to the next reach.** A reach of
one point does not stop at a junction, so a collector can run through one; a tooth ending there was
linked by the junction's receiver to the collector's NEXT reach, `own_end` did not see the confluence,
and the tooth exported the collector's discharge: in the rectangle, teeth of **0.2–1.1 km² exporting
32–37 m³/s and 28–30 m** (an implied runoff of 16 000–87 000 mm/yr). Fix: a reach ending on a junction
another reach holds enters THAT reach. After: 0 such links held by a watercourse, both worlds; 293 / 54
reaches end on a cell only a SPILLWAY runs through — a spillway is traced over a col outside the D8
network, so this is the Watercourse × Spillway overlap the objects handle (§ 3.3), not a link.

Each fix has a permanent test that fails on the old code (checked by putting the old lines back):
`drainage::the_clip_links_where_the_water_goes_and_reads_the_own_value` (T1, T2),
`flow::a_tributary_links_to_the_reach_its_water_enters` (T3, T3b).
`ALGO_DRAINAGE` 8 → 9 and `ALGO_HD_DRAINAGE` 10 → 11, so no cached drainage is served stale.

### 3.3 · The objects (the microscope's aggregation)

The shipped aggregation grouped by terminal and showed, as "Bassin", `runoff / 300 mm` — a
runoff-EQUIVALENT area. Three gestures, `aggregate_watercourses` (viz), ported in the bench:

- **one reach, one object**: systems whose reaches share a cell are ONE object (union of their
  terminals) — the common downstream is one trunk;
- **the chain crosses every body that overflows**: a reach dying in a water body with a spillway is
  chained to it, as Finding 45 chains an exorheic lake to its outlet; the climb back up is the exact
  inverse of that chaining (a first version rebuilt it separately and climbed into reaches the chaining
  had given to another object — 42–45 km trunks on < 0.4 km²);
- **the basin is an AREA**: `catchment_km2` is the object's share of the GEOMETRIC partition of the
  land — every land cell follows D8 to the first river cell or water body it meets.

| | C2 /10 shipped | C2 /10 F124 | delivered shipped | delivered F124 |
|---|---|---|---|---|
| objects | 636 | 479 | 384 | 362 |
| river cells in > 1 object | 550 | **0** | 266 | **0** |
| Σ "Bassin" / land | 1.40 (runoff-equivalent) | **0.860** (geometric) | 1.34 | **0.844** |
| uncaptured coastal fringe | — | 0.140 | — | 0.156 |
| longest trunk, domain (shipped: points × cell; F124: path length, § 3.4) | 224 / 156 / 146 km, sinuosity 9–14 | 151 km, sinuosity 1.5 | 44 km | 103 km, 1.7 |
| trunks losing a lake crossing vs pre-F123 (Finding 93's control) | — | **0** | — | **0** |

Σ over the land is 1.000 exactly (partition + fringe, counted in cells). The runoff-equivalent sum
stayed above 1 whatever the rule (1.20–1.40): the unit, not the partition, was the lie — **the
runoff area is 1.51× the geometric area on average** on these worlds.

⚠️ **One reading of the author's rule is mine.** "Un tronçon appartient à un seul objet — l'aval commun
est un seul tronc, les affluents des objets séparés": implemented as ONE list entry per system, whose
trunk is the common downstream and whose tributaries are its members (the viz's model since Finding
28). If "des objets séparés" meant a separate list entry per tributary joining the trunk, that is a
change of model (Finding 37's Point 4: a named main stem with each tributary as its own watercourse and
a link to the stem) and it is not done.

### 3.4 · Length and the rectangle

**What the 1 093 km summed.** Points of the trunk's 82 reaches × km per cell × 7.5: the max-area climb
went up a comb of parallel wall columns, 2 633 of its points inside the 59 × 48 block. Two defects made
the climb do that: T2 (every tooth carried its collector's area, so they tied) and T3 (teeth linked to
teeth). T3 alone accounts for most of it: the pre-F123 LONGEST-PATH rule's 224 / 203 / 151 km trunks
fell to 61 km with T3 fixed.

**Length is now the length of the max-area path**: Euclidean steps on the torus (1 or √2), the
confluence cell counted once (it was counted twice: +7.6 % on the 61 km trunk), a crossed water body by
its chord. **Asserted ≤ the domain diagonal** (566 km (domain)) on every corrected trunk in
`f124_objects` — it holds. The longest trunks (`objects_path_length.txt`):

| | longest trunks, domain (signified) | straight line | sinuosity |
|---|---|---|---|
| C2 /10 | **151 km** (1 132 km), 131, 81, 68 | 102, 91, 61, 50 km | 1.5, 1.4, 1.3, 1.3 |
| delivered | **103 km** (771 km), 98, 76, 75 | 61, 53, 64, 47 km | 1.7, 1.8, 1.2, 1.6 |

The two longest C2 /10 systems cross the big below-sea basins by their chord (their straight lines
are 91–102 km); by point count, before the chord, they read 79 and 61 km.

**The rectangle, attributed.** The 2 459-cell block is a collector running east along y = 3026
(reaches of 3 points, own area 64–78 km², an implied runoff ~275 mm/yr, 28–31 m wide) fed by a comb of
parallel teeth running SOUTH (D8 south on 2 557 of the block's 2 632 river cells), 33–51 cells long,
0.2–1.1 km² each, on 100 % carved (constructed) terrain. **The comb is the FIELD's** — Finding 123's
planar walls, the target of part 4's B — **and part of its look was the EXPORT's**: teeth drawn at the
collector's width (T3b). After the fix the block's teeth read p50 1.4 m, p90 4.5 m (29.7 / 30.6 m
before); the collector keeps its 28–31 m.

## 4 · Method

- **An instrument that sums 67 M cells counts in integers.** The first partition summed f32 cell areas
  and read Σ + fringe = 1.149 of the land.
- **An inverse relation is built as the inverse.** The first corrected climb across a body rebuilt its
  own notion of "the inflow of this outlet" instead of inverting the chaining, and walked into another
  object's reaches.
- **A census's holder map keeps every holder.** The first T3b census kept one reach per cell and
  counted spillway overlaps as link defects (392 / 59); with every holder, 0 held by a watercourse.
- Rule 11 greps: `cell_to_segment`, "stale", "asymmetr", "dangling", "junction cell" in the ADR and
  docs — **NOTHING FOUND** for T1 and T3; Finding 47 had PROPOSED the penultimate-point rule for a
  shared MOUTH (L4279) and left it to the author, who called it for the mouth ("a system reports its own
  catchment, never the union"); it had never reached the confluences.

## 5 · Predictions, scored

| | prediction | result |
|---|---|---|
| P1-a | first non-periodic stage = FBM / coast warp, not the mask | **held**: S1 (coast warp, unwrapped `sx`) — FBM adds nothing |
| P1-b | F108 window holds at the canonical framing | **held** |
| P2 | k within ± 1 % | **held** (0.07183 vs 0.07186) |
| P3-a | the 1 093 km / 272 592 km² are partly signified units | **half**: true of the units, but the trunk was also false (T2 / T3) |
| P3-b | shared cells are trunk / shared-mouth cells, 1–5 % | **refuted**: 0.14 %, a Watercourse × Spillway pair |
| P3-c | the rectangle is a flat the router fills | **refuted**: a comb of south-running teeth on carved walls, flat 0 % |
| P3-d | own area straightens the trunk (≤ 60 km, sinuosity < 3) | **partly**, until T1 / T3 were fixed |
| P3-e | union + chaining bring Σ / land ≤ 1 | **refuted as stated** (1.30); its fallback — the unit — held |
| P3-f | 0 trunks lose a lake | **refuted** (3 on delivered) → T3 found, then 0 |
| P3-g | 1–3 % stale links, fewer on delivered | **refuted**: 5.4 % / 5.0 % |
| P3-h | 40–60 % end on the confluence | **half**: 57 % / 73 % |
| P3-i … P3-p | after each fix | **held** (0 asymmetric, ratios 1.00, T3 0, F93 0, heads 37 → 0, T3b 0 on watercourses, objects stable) |

## 6 · Break the plane: B2 and B1 (part 4)

`f124_walls`, on the corrected export (part 3 first: "un export qui ment sur les troncs ne peut pas
juger des tributaires"), the definition (basin base) at the canonical framing, k_time/10 + A1+B2.
Teeth = watercourse reaches with ≥ 80 % of their cells on the world's OWN constructed walls (Finding
123-A); ex-teeth = C2 /10's teeth paired in each world (Finding 123's rule-18 pairing: own mouth ≤ 3
cells, own area ×/÷ 1.5). Output: `../f124_walls/walls.txt`.

| world | teeth | length share | ex-teeth sin p50 | ex-teeth R8 c8 | R8 terrain | σ p50 | Δ e/b | coast | skeleton ≤ 5 at 1 km² | mouths | relief p50 | build |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C2 /10 (reference) | 1 440 | 8.8 % | 1.027 | 0.347 | 0.045 | 5.6 m | 0/0 | +3 | 74.2 % | 509 / 877 | 515 m | 157 s |
| **B2, tributaries ≥ 1 km²** | **262 (−82 %)** | 1.5 % | 1.049 | **0.098** | **0.105** | 9.4 m | 0/0 | +4 | **91.0 %** | 652 / 877 | 446 m | 177 s |
| B2, down to A_c 0.1 km² | 73 (−95 %) | 0.5 % | 1.070 | 0.357 | 0.125 | 12.3 m | **4/0** | +2 | 91.6 % | 671 / 877 | 420 m | 154 s |
| B1, wall profile + detail | **4 719 (+228 %)** | 16.3 % | 1.030 | 0.403 | 0.048 | 5.5 m | 0/0 | **+91** | 74.9 % | 410 / 877 | 500 m | 140 s |
| B2 (A_c) + B1 | ⛔ Finding 38 | | | | | | | | | | 407 m | |

References (the replay, same instrument): delivered R8 terrain 0.092, A1+B2 0.048, PRE 0.013; Finding
123's A1+B2 homologues of the teeth: sinuosity 1.044.

**B2 at 1 km² is the only remedy that makes the teeth fall and keeps every gate**: −82 % teeth, the
paired ex-teeth no longer straight parallel lines (R8 at chord 8 0.347 → 0.098, sinuosity 1.049 —
their A1+B2 homologues' 1.044), canyons 0, coast +4, the skeleton at 1 km² STABLE (91 % within 5
cells — Finding 120 had it only at 100 and 1 000 km²), relief 446 m (within ±10 % of the delivered
424 m for the first time in a construction world), +13 % build time. **But the terrain's R8 doubles**
(0.045 → 0.105, above the delivered 0.092) and σ rises 5.6 → 9.4 m. Not attributed yet; the likely
reading — B2's own tributary valleys are built with the same planar walls, so the anisotropy returns
one scale down — is a hypothesis, and rule 18 applies before any use of it.

**B2 down to A_c fails the canyons gate** (Δ 4/0: four over-dug bodies). **B1 fails twice**: its
detail term (field − trend put back on the wall) roughens the walls so MORE D8 lines run on them
(teeth × 3.3) without bending the old ones (sinuosity 1.030), and it reaches the coast (+91 spurs).
**B2 + B1 violates the production invariant of Finding 38 / 92-B** (`hd_assembly.rs:279`): two
enclosed below-sea components of 1 and 2 cells, floor cells (4257, 2006) and (4613, 2860), carry no
water body — B1's detail digs micro-pits below the sea on B2's walls. The invariant did its job; the
world is not measured beyond its relief (407 m) and its 284 lakes.

The round's question — "laquelle fait tomber les dents sans monter R8 terrain ?" — **neither**. B2 at
1 km² is the candidate, with its R8 terrain to attribute. Not run, a suggestion only: B2 with B1's
PROFILE but no detail (`detail_gain` 0) — the detail is what broke the coast and dug the pits.

**Predictions.** Mine (P4): B2 removes > 70 % — **held** (82 / 95 %); B1 alone < 50 % — **refuted in
sign** (+228 %); skeleton at 1 km² within 5 cells on 70–85 % — **half** (74 % on C2 /10, 91 % on B2).
The reviewer's: B2 > 80 % and ex-teeth sinuosity ~1.04 — held (82 %, 1.049); B1 ~−40 % — refuted;
B2 costs ×2–3 — refuted (×1.13); the skeleton holds at ~70 % — refuted (91 %).

## 7 · W(k) (part 5)

`W = a₀·(k/k₀)^γ·A^0.3`, `k₀ = F121_AGE_K`: the age knob that deepens the floors also widens them.
**γ is a PROXY** — lateral widening is a matter of time, and Clubb 2022's `a` is a snapshot of present
valleys, not a rate. Gated: `ValleyConstruction::width_age_gamma`, `None` = Findings 121–124's width
(permanent test `the_width_widens_with_age_only_when_asked`). Measured with Finding 121's instrument
(`f121_hybrid`, tag `e`, run only when asked), C2 /10, three ages × three γ. Output:
`../f124_width/width.txt`; renders `../f124_width/tile_k*_g*.png` (the comb tile) and
`lake_1000011_k*_g*.png` (the brush lake's crop at the canonical framing), NORTH UP. ⚠ The renders
draw the PRE-CLIP network (every river cell of `c1_drainage_windowed`), not the corrected export.
⚠ **The lake crops do NOT show lake 1000011**: it is a below-sea basin lake (id ≥ 1 000 001), made
by the HD assembly, and the renders read the pre-breach detection's `lake_map`, which has no such
lake — its basin reads as bare ground. The crops show the terrain around it, not the lake. Rendering
it needs the HD drainage per world (~1 min each); not done.

| k | γ | law W p10 / p50 / p90 | relief p50 | lake % | lakes ≥ 1 | Δ e/b | D_L > 5 | coast | R8 terrain |
|---|---|---|---|---|---|---|---|---|---|
| × 0.7 | 0 | 212 / 274 / 524 m | 501.0 m | 12.76 | 17 | 0/0 | 0 % | +4 | 0.0477 |
| × 0.7 | 0.5 | 178 / 229 / 438 m | 503.6 m | 12.77 | 17 | 0/0 | 0 % | +5 | 0.0467 |
| × 0.7 | 1 | 149 / 192 / 367 m | 505.8 m | 12.77 | 17 | 0/0 | 0 % | +5 | 0.0465 |
| × 1 | any | 212 / 274 / 524 m | 515.4 m | 12.93 | 18 | 0/0 | 0 % | +3 | 0.0451 |
| × 1.4 | 0 | 212 / 274 / 524 m | 533.3 m | 13.10 | 19 | 0/0 | 0 % | +3 | 0.0405 |
| × 1.4 | 0.5 | 251 / 325 / 620 m | 530.9 m | 13.10 | 19 | 0/0 | 0 % | +3 | 0.0418 |
| × 1.4 | 1 | 297 / 384 / 733 m | 528.1 m | 13.10 | 19 | 0/0 | 0 % | +3 | 0.0432 |

At k₀ γ changes nothing (k/k₀ = 1, by construction). **The width moves the relief by ±5 m at most**
(+4.8 m at k × 0.7, −5.2 m at k × 1.4, γ 1): the age's own effect (501 → 533 m through the floors)
is six times larger. Nothing else moves — lakes, the Δ class, D_L, the coast. **By eye, at the tile's
50 km, γ 0 and γ 1 at k × 1.4 are hard to tell apart**: slightly broader floors on the big trunks, and
that is all. W(k) is a legible knob with a small footprint; the author's dial is in the viz ("W(k) γ ·
PROXY", EXPERT mode, beside the age), γ 0 keeping the bench's digest.

**Predictions.** P5-a (arithmetic) — as computed. P5-b — **held** (+4.8 m within 2–8; −5.2 m at the
edge of 5–15). P5-c — **held** (Δ 0/0 everywhere, D_L 0 %, lake % within 0.01 point). P5-d — **held**:
γ 1 at k × 1.4 is not visibly "too wide" at the tile scale. The reviewer's E ("γ 1 too wide at k ×
1.4") is not borne out by the eye at this scale; "γ 0.5 reads old without invented" cannot be judged
from renders this similar — the author's eye, not mine.

## 8 · Addendum — the author's answers (2026-09-27)

- **Objects**: one entry per system, the tributaries its members (*"le tronc est une ligne continue,
  tout le reste est secondaire"*). **The "654 km"**: closed, *not reproduced, same cause as the
  1 093 km* — 654 km (signified) = 87 km (domain), a climb through the teeth.
- **#3** (S4, 210 287 km² signified) contains the mouth (3824, 4757) — the image 4–5 pair, fused.
  **#76** (S2, 5 721 km² signified = 101.7 km² domain) is by its area the mouth-(4936, 3051) object
  whose trunk crossed the block; its 17 m³/s (signified) against the collector's 38.3 m³/s is not
  explained (ADR Finding 124-P3).
- **B2's skeleton was taken after T1–T3b**, and it does not read the network at all (ADR Finding
  124-P4); its R8 × 2 is not the error's. The hypothesis for it is written under details there.
