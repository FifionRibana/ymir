# Finding 125 — attribute before remedying: where B1 makes teeth and spurs, where B2 stripes, where B2 → A_c digs

**No remedy this round; no primitive modified.** Every number below comes from
`walls_attribution.txt` (bench `crates/ymir-core/tests/f125_walls.rs`, 11 656 s) or
`guard_six_states.txt`.

**Units.** The domain is 400 km on 8192² cells, so one cell = 48.8 m (domain). "Signified" means ×7.5
for lengths and ×56.25 for areas. Every length and area below says which one it is. Heights are in
metres and are not scaled. The drained-area bands (≥ 10, 3–10, 1–3, < 1 km²) are domain km²: the
skeleton computes `area_km2` with `domain_km`.

**Reading declaration, unfavourable.** My predictions (`f125_predictions.md`, scratchpad) were written
before any grep or measurement. They are non-blind on Findings 73–124 and on the reviewer's
predictions. Two of my own scorings were wrong on first writing and are corrected below: the big
walls' −8 % (I first wrote −5 %), and D's "coastal" (I first scored it held; it was not measured).

## 0 — The prerequisite: the guard through `run_hd`, 6 / 6 "= banc" after one fix

**First run: 2 / 6 (livré, A1+B2). The stop rule fired.** The four construction states read
`NoReference`. Their viz FIELD hashes nevertheless equal the reference's to the bit, so the worlds
ARE the bench's; the KEYS moved.

The cause: `eroded_key_full` folded the construction with `with_debug`, i.e. its Debug string. C3
(1a7d5d3, Finding 124-P4/P5) added two gated fields that print as
`wall_profile: None, width_age_gamma: None`. A gated-off field therefore changed every construction
key.

**The fix, the author's choice ("clé robuste + référence").**
- The construction enters the key by its serde form (`key.with`). `skip_serializing_if` drops a gated
  field left `None`.
- A new permanent test, `a_gated_off_construction_field_does_not_move_the_eroded_key`
  (`cached_product.rs`), checks three things:
  - the gated-off fields are absent from the serialised form;
  - the key equals a hand-rebuilt key (it FAILS on the old `with_debug` form: this is the negative
    control, verified);
  - turning `wall_profile` on DOES move the key.
- The four stored digests were refreshed. The field hashes did not change:

| state | old digest | new digest | field hash (unchanged) |
|---|---|---|---|
| C1 nue | 54545be015b5ff19 | e2a5b05da6eed1fb | c876c200e92e0f60 |
| C2 /10 col (défaut) | 71a1d832683583d1 | cf1d4539c7afe588 | a8d2d538d692c2f0 |
| C2 /3 | 111a7705c9a4131c | 70cf2dc92822681b | c2eb1cb3b24ebccb |
| C2 /10 niveau mer | 01c18a31ec7d0fc0 | 3552097350b238c8 | 35e1bc8123049d8c |

livré (f61ef23f4f4e9d94) and A1+B2 (e05628c2eb1c1c22) are unchanged. The viz guard test
`f123_viz_guard` now covers the six states.

**Second run: 6 / 6 "= banc"** (`guard_six_states.txt`).

## 1 — Worlds, the témoin, and what the bench reproduces

**Worlds.** All use the definition: canonical framing, k_time/10 + A1+B2.
- the témoin, C2/10 col;
- B1a, the profile alone (`detail_gain` 0);
- B1b, the noise alone: a planar wall (`foot_m` 0, `crest_m` 0);
- B1, both (`WallProfile::f124`);
- B2 at 1 km² (domain);
- B2 down to A_c (0.1 km², domain);
- B2 (A_c) + B1.

**Rule 18.** Every population row reads the témoin on the SAME cells, meaning the world's own classes
applied to the témoin's field.

**The bench reproduces Finding 124-P4.**
- Teeth match exactly: 1 440 / 4 719 / 262 / 73.
- Terrain R8 matches within 0.0006: 0.0452 / 0.0490 / 0.1049 / 0.1253, against F124-P4's 0.0451 /
  0.0484 / 0.1049 / 0.1251. The residual comes from the window mask: here a window is dropped if a
  PRE cell is sea.
- F124-P4's "+91 spurs" for B1 counts against the delivered world (the témoin itself read +3 there).
  Against the témoin it is +88 at the breached stage (111 − 23), which is what this bench reads.

**Instruments.**
- **Teeth**: F123-A's instrument. A tooth is a watercourse reach with ≥ 80 % of its points on wall
  cells (carved, not floor), and it is attributed to the class of its mouth.
- **Terrain R8 by population**: 16×16-cell windows (781 m, domain) that are fully land in PRE. Each
  window's `e^{8iθ}` is split by the share of its cells in each class: interfluve, or a carved cell
  labelled by the band of its NEAREST trunk cell (chessboard BFS). "Class R8" = |V_c| / N_c.
  "Contribution" = Re(V_c · conj(V̂)) / N; the contributions add up to the total exactly.
- **Spurs**: Finding 75's detector, at the eroded and at the breached stage. A spur is "new" if it lies
  > 2 km (domain) = 15 km (signified) from every témoin spur of the same stage. A "coastal wall" is a
  cell that is carved, not floor, land, and ≤ 3 cells (146 m, domain) from the sea.
- **Over-dug**: Finding 108's class, with the floor cell. "Cones" is the number of skeleton polylines
  whose planar cone reaches below the PRE field at that cell.
- **F38**: the invariant's own panic message, parsed for its floor cells.
- **Skeleton × planarity**: 8-point chords of the 1 km² skeleton. Planarity is the RMS residual of a
  plane fit to PRE over 21×21 cells (1.03 km, domain), cut into quartiles.

## A — B1 separated: the teeth need both halves; the spurs come from the profile, and they are coastal

| world | teeth (Δ vs témoin) | teeth where (mouth class) | spurs eroded: new / within 2 km (domain) of a coastal wall | spurs breached: new / coastal | lake bodies | terrain R8 (Δ) |
|---|---|---|---|---|---|---|
| **témoin C2/10 col** | **1 440** | trunk ≥ 10: 1 440 | 20 (reference) | 23 (reference) | 25 | 0.0452 |
| B1a profile only | 2 130 (**+48 %**) | trunk ≥ 10: 2 125, interfluve 5 | **42 / 42** | 41 / 40 | 81 | 0.0497 (+0.0045) |
| B1b noise only | 1 681 (**+17 %**) | trunk ≥ 10: 1 681 | **8 / 8** | 8 / 8 | 87 | 0.0452 (−0.0001) |
| B1 both | 4 719 (**+228 %**) | trunk ≥ 10: 4 702, interfluve 17 | **81 / 81** | 83 / 80 | 102 | 0.0490 (+0.0038) |
| B2 ≥ 1 km² (no B1) | 262 (−82 %) | trunk 180, 3–10: 37, 1–3: 44, interfluve 1 | **0** | 1 / 0 | 25 | 0.1049 |
| B2 → A_c (no B1) | 73 (−95 %) | trunk 21, 3–10: 1, 1–3: 3, < 1: 48 | **0** | 0 | 33 | 0.1253 |
| B2 (A_c) + B1 | — (F38) | — | **195 / 195** | 195 / 193 | 284 | 0.1032 |

**The teeth are superadditive.** Profile alone: +690. Noise alone: +241. Both: **+3 279**, 3.5 times
the sum. Every tooth, in every B1 world, has its mouth on a big valley's wall (trunk ≥ 10). **Neither
half carries the teeth; their interaction does.** The noise is the weaker half (+17 %), so the
reviewer's "the teeth come from the noise" and my P-A1 are refuted. WHERE ON THE WALL the teeth sit
(foot, straight section, rim) is not measured: the instrument locates a tooth by its mouth's class
only.

**The spurs come from the PROFILE.** B1a adds 42 new spurs and B1b adds 8. Both together add 81, more
than 42 + 8. They are present from the ERODED stage: the breached counts equal the eroded ones within
±2, so Finding 90's "spurs are born at the breach" does not describe these.

**Where the spurs are, and the stop rule.** All new eroded spurs lie within 2 km (domain) of a coastal
wall: 42 / 42, 8 / 8 and 81 / 81. The first six listed per world sit **1 to 10 cells (49–488 m,
domain) from a coastal wall**. **The stop rule does not fire: the spurs ARE coastal, and they sit on
walls that reach the sea.** "mur ↔ mer" stands.

That test is weaker than it reads, for three reasons:
- B1 only modifies wall cells, so a new B1 spur can only be born near carved cells. "Within 2 km of
  any carved cell" reads the same count (42 / 8 / 81).
- The base rate is missing: the share of the coastline within 2 km (domain) of a coastal wall is not
  measured.
- The discriminating evidence is elsewhere:
  - **B2 alone multiplies the walls, coastal ones included, and makes ZERO new spurs** (0 at 1 km²,
    0 at A_c);
  - **B2 (A_c) + B1 makes 195**, 2.4 times B1's 81.

So walls do not make spurs. The profile on a wall that meets the sea does, and the count scales with
the coastal walls B2 adds. That ties the effect to the wall × sea contact, not to "wall" alone.

**Costs in A.** B1a or B1b alone multiply the lake bodies by 3.2–3.5 (25 → 81 / 87), and both by 4.1
(→ 102). The lake AREA does not move (12.86–12.91 % of land, against 12.93 %), so these are small
bodies. They are not located. Relief p50 does not move (B1 +18.0 %, B1a +17.3 %, B1b +22.1 % against
the témoin's +21.5 %, all against the delivered 424.2 m).

## B — B2's R8: the sub-valleys themselves, not the big walls; planarity is a step at the roughest quartile, not a law

**B2 ≥ 1 km², decomposed by population.** The témoin is read on the same classes.

| population | windows | class R8, world / **témoin** | contribution, world / **témoin** | Δ contribution | share of Δ |
|---|---|---|---|---|---|
| interfluve | 22 758 | 0.0267 / **0.0162** (+65 %) | +0.0135 / **+0.0072** | +0.0063 | 11 % |
| trunk ≥ 10 walls | 3 158 | 0.2161 / **0.2351** (−8 %) | +0.0158 / **+0.0170** | −0.0012 | −2 % |
| sub 3–10 km² | 4 198 | 0.2557 / **0.0857** (× 3.0) | +0.0249 / **+0.0082** | +0.0167 | **28 %** |
| sub 1–3 km² | 12 993 | 0.1689 / **0.0426** (× 4.0) | +0.0508 / **+0.0128** | +0.0379 | **63 %** |
| **total** | 43 107 | — | 0.1049 / **0.0452** | +0.0597 | 100 % |

**91 % of B2's added terrain R8 is in the sub-valleys (1–10 km², domain).** On the same cells their own
class R8 is 3–4 times the témoin's. The big walls do not change (−8 %). The interfluves rise +65 % and
carry 11 %. Those windows are fractional, so a window that is mostly interfluve but also holds sub-
valley cells takes the sub-valley's orientation. That share may be contamination (see the
limitations).

In the témoin's own classes, 90 % of its R8 already lives in the big walls (+0.0405 of 0.0452).
**The terrain's anisotropy has always been the carved walls'. B2 adds walls, so it adds anisotropy.**

The same reading holds for the other two B2 worlds:
- B2 → A_c: sub < 1 carries 83 % of the Δ (+0.0661 of +0.0800), and the sub-valleys 91 % in total.
- B2 (A_c) + B1: sub < 1 carries 89 % (+0.0518 of +0.0580).

**The stripes are in the terrain, not in the exported network.** The network's R8 (chord 8) FALLS
under B2: 0.5613 (témoin) → 0.5060 at 1 km² and 0.5356 at A_c. The routed watercourses are less
axis-aligned than the témoin's while the carved terrain is more.

**The 1 km² skeleton, by band × planarity quartile of PRE.** Values are R8 of 8-point chords. The
residual in parentheses is each quartile's upper bound, in m (vertical).

| band (at the chord's mid cell) | chords | Q1 | Q2 | Q3 | Q4 (roughest) | Q1 / Q4 |
|---|---|---|---|---|---|---|
| **trunk ≥ 10 (≈ the témoin's own skeleton)** | 18 403 | **0.459** (≤ 1.4) | **0.497** (≤ 2.9) | **0.407** (≤ 6.7) | **0.170** (≤ 92.2) | 2.7 |
| sub 3–10 km² | 22 098 | 0.437 (≤ 1.3) | 0.467 (≤ 2.5) | 0.449 (≤ 5.1) | 0.386 (≤ 103.4) | 1.13 |
| sub 1–3 km² | 50 683 | 0.635 (≤ 1.3) | 0.615 (≤ 2.5) | 0.644 (≤ 4.6) | 0.494 (≤ 82.9) | 1.29 |
| sub < 1 (mid cell off the trunk, see limits) | 19 044 | 0.112 | 0.143 | 0.165 | 0.141 | 0.79 |

**The question "is the 1 km² skeleton straight where PRE is planar?" gets a split answer:**
- **Q1–Q3 are FLAT in every band.** The 1–3 km² skeleton reads 0.62–0.64 whether the 1 km (domain)
  plane residual is 1.3 m or 4.6 m.
- The roughest quartile alone lowers it: to 0.49 at 1–3 km², and to 0.17 for the trunks.

This is not a graded planarity law, and it is weakest exactly where B2 adds its R8 (Q1/Q4 1.29 at
1–3 km², 1.13 at 3–10 km²). It is a rough / not-rough step, and 75 % of the chords are on the
"not rough" side.

The 1–3 km² branches are more axis-aligned than the trunks **at every quartile** (0.62–0.64 against
0.41–0.50). The alignment belongs to the small branches' TRACÉ. The terrain they cross does not cause
it.

**Consequence for the future remedy (named, not built).** A hierarchy that stops "where D8 is
straight" would stop on Q1–Q3, i.e. on 75 % of the 1–3 km² skeleton. That is nearly not doing B2. The
numbers favour the reviewer's other candidate: act on the tracé (smooth the skeleton before
construction). No D8 baseline was measured: the R8 of 8-point chords of a D8 trace over isotropic
terrain. So "0.62 is D8's own axis preference" is an inference, not a measurement.

## C — The four canyons of B2 → A_c, located

The témoin has 0 over-dug bodies. So does B2 ≥ 1 km². B2 → A_c has 4:

| body | area, domain / signified | floor cell (x, y) | Δfill | rim p50 | floor class | cones at the floor | in B2 (A_c) + B1? |
|---|---|---|---|---|---|---|---|
| 3 | 1.49 / 83.8 km² | (3892, 2310) | 107.5 m | 31.9° | **sub < 1** | 168 | no |
| 7 | 10.55 / 593.4 km² | (4588, 3240) | 140.6 m | 31.9° | **trunk ≥ 10** | 170 | no |
| 13 | 15.37 / 864.6 km² | (2281, 4563) | 208.3 m | 31.2° | **trunk ≥ 10** | 39 | **yes**: body 14, floor (2283, 4573), 44 cones |
| 14 | 41.97 / 2 360.8 km² | (1835, 4562) | 358.7 m | 32.1° | **trunk ≥ 10** | 34 | **yes**: body 15, same floor, 34 cones |

- **All four are on a constructed FLOOR** (carved and floor).
- **Three of the four are on a big valley's floor** (trunk ≥ 10). One is on a < 1 km² valley's floor.
- **They appear only when the construction goes below 1 km².**
- **The two largest recur when B1 is added**: one on the same floor cell, one 10 cells (488 m, domain) away.

"How many primitives overlap there" does not decide the reviewer's "intersections" hypothesis. The
cone count is 34–170 at every floor. It counts every polyline whose cone reaches below PRE at the
cell, including the valley's own polylines and any deep neighbour, and it has no baseline (the same
count at an ordinary trunk floor of B2 → A_c, or at these cells in B2 ≥ 1 km²).

**Measured:** the canyons need primitives < 1 km², and three of four dig the floor of a BIG valley.
That is where the small valleys meet a trunk, not where two small valleys cross. **Not measured:**
whether a meeting is the mechanism.

## D — B2 (A_c) + B1 and Finding 38: two wall cells, not floors, and not at the canyons

The invariant reports 2 enclosed below-sea components carrying no water body.

| floor cell (x, y) | cells | height | PRE height | class | carved / floor | PRE detail term (500 m box, domain) | cones |
|---|---|---|---|---|---|---|---|
| (4257, 2006) | 1 | **−1.00 m** | 98.7 m | sub < 1 | true / **false** (a wall) | −75.9 m | 5 |
| (4613, 2860) | 2 | **−1.00 m** | 182.6 m | sub < 1 | true / **false** (a wall) | −92.5 m | 11 |

**The same cells in each ingredient's world (the témoin family and B2 alone):**

| cell | B1a | B1b | B1 | B2 ≥ 1 km² | B2 → A_c | B2 (A_c) + B1 |
|---|---|---|---|---|---|---|
| (4257, 2006) | 87.69 m | 68.85 m | 88.40 m | 6.06 m | 4.57 m | **−1.00 m** |
| (4613, 2860) | 1.95 m | 2.73 m | 2.20 m | 5.17 m | 36.05 m | **−1.00 m** |

- **Each ingredient alone leaves both cells above the sea.** The lowest readings are 1.95 m (B1a) and
  4.57 m (B2 → A_c). Only the combination crosses.
- **Both pits are WALL cells of < 1 km² valleys, not floors.** Neither is at a canyon: they lie 23.2
  km and 18.6 km (domain), i.e. 174 km and 139 km (signified), from the nearest canyon floor of C.
- In PRE, both cells sit in local lows at the 500 m (domain) scale (detail term −76 / −93 m). That is
  the sign B1's detail term would put back.
- **Correction:** the instrument's "distance to the sea 0 cells" is TAUTOLOGICAL. The pit is itself ≤
  SEA, so it is its own "sea" cell. **Whether these walls are coastal is NOT measured.**
- Both read exactly −1.00 m. A first grep for a clamp (`1 m below`, `sea_level - 1`, `SEA_M - 1`,
  `-1.0 … sea` in `ymir-core/src`) found NOTHING. Not investigated further.

**D reads as "B1's lowering applied to the walls B2 → A_c adds". It is not "A + C": C's canyons are
floors, and these are walls far from them.**

## E — What it draws: the missing clauses (WRITTEN, NOT MEASURED)

| unwanted effect | where (measured) | the gesture that carries it (measured) | the clause that is missing (written) | status |
|---|---|---|---|---|
| **B1 teeth** +228 % | every mouth on a trunk ≥ 10 wall (4 702 / 4 719); position along the wall not measured | **the interaction** profile × noise (+48 % and +17 % alone, +228 % together) | **bruit ↔ seuil, conditioned on the profile**: the noise must not exceed the local slope where the profile flattens the foot. Not "bruit ↔ seuil" alone, since noise alone gives +17 %. | attributed to an interaction; the clause is a candidate; the "where on the wall" is still missing |
| **B1 spurs** +81 (eroded) | every new spur within 2 km (domain) of a coastal wall, the listed ones at 49–488 m (domain) | **the profile** (42 vs the noise's 8); scales with B2's coastal walls (195); walls alone make 0 | **mur ↔ mer**: F121's "a wall stops where it meets the terrain" has no "and where it meets the sea" | supported (the stop rule does not fire); not proven (no coastline base rate) |
| **B2 terrain R8** × 2.3 | 91 % in the sub-valleys 1–10 km² (1–3: 63 %, 3–10: 28 %); big walls −8 %; interfluves +11 % of the Δ | **the sub-valleys' own walls**: class R8 3–4 × the témoin's on the same cells; the network's R8 falls | **NOT hiérarchie ↔ planéité** (flat over Q1–Q3; a planarity stop would remove 75 % of the 1–3 km² skeleton). **hiérarchie ↔ tracé**: the small branches' axis alignment is the tracé's (0.62–0.64 against the trunks' 0.41–0.50 at the same planarity) | the reviewer's clause not supported; replacement named; no D8 baseline |
| **B2 → A_c canyons** (4) | 3 on trunk ≥ 10 floors, 1 on a < 1 km² floor; absent at 1 km² | **the < 1 km² primitives** | **mur ↔ mur**, in the form "fond ↔ fond": a tributary's construction must not pose its floor under the trunk's where they meet | UNDECIDED: the cone count has no baseline |
| **F38** in B2 (A_c) + B1 | 2 WALL cells of < 1 km² valleys, −1.00 m, 18.6 and 23.2 km (domain) from any canyon | **the combination**: each ingredient alone ≥ 1.95 m above the sea | **bruit ↔ seuil (the sea's)**, a form of mur ↔ mer: a lowered wall must stop above the sea unless it opens onto it. Not mur ↔ mur. | consistent; coastal NOT measured |

The next round's remedy is the clauses, not a fifth gesture. The table rules out one of the reviewer's
four clauses as stated (hiérarchie ↔ planéité). It reshapes two: bruit ↔ seuil needs the profile, and
mur ↔ mur is fond ↔ fond at a trunk. It supports one without proving it (mur ↔ mer).

## Predictions scored

**Mine.**
- **P0 REFUTED**: 2 / 6 on the first run. The cause was neither of my two guesses (a state's params, a
  label); it was the key's Debug form.
- **P-A1 REFUTED**: B1b +17 %, not ≥ 2×; B1a +48 %, not within ±30 %; superadditive.
- **P-A2 REFUTED on the gesture, HELD on the place**: the spurs come from the profile, not the noise.
  They are coastal (100 % within 2 km, domain).
- **P-B1 MOSTLY HELD**: the sub-valleys carry 91 % (≥ 60 %), with 1–3 dominant; big walls −8 % (within
  ±10 %). The interfluves at +65 % are OUTSIDE my ±10 %.
- **P-B2 REFUTED where B2 adds R8** (1–3: Q1/Q4 1.29; 3–10: 1.13); held only for the trunks (2.7).
- **P-C HELD TRIVIALLY**: ≥ 2 cones everywhere, and the metric has no baseline.
- **P-D HALF**: on walls, held; not at the canyons, held; coastal NOT measured.

**The reviewer's.**
- **A: the teeth from the noise, REFUTED** (+17 % against the profile's +48 %; the teeth need both).
  **The spurs from the profile on sea-reaching walls, HELD**; the mechanism ("the convex rim
  overflows") is not measured.
- **B: > 80 % in 1–3 km², REFUTED in degree** (63 %; 91 % over 1–10 km²). **"A strong correlation
  between skeleton R8 and planarity, so the hierarchy threshold is a planarity", REFUTED** in the bands
  that carry the Δ.
- **C: 3 of 4 at intersections, UNDECIDED** (the instrument cannot separate an intersection). The
  measured split is 3 trunk floors + 1 small floor. It matches in number only.
- **D: F38 = A + C, HALF**: A's part is consistent (B1 on walls); C's part is refuted (walls, not the
  canyons' floors, ≥ 18.6 km domain away).
- **Meta ("at least one of the four is wrong"), HELD.**

## Limitations, stated

1. **No base rate for the spurs' proximity.** The share of the coastline within 2 km (domain) of a
   coastal wall is not measured. The evidence for mur ↔ mer is the B2-alone control (0 spurs), not
   the 100 %.
2. **No baseline for the cone count.** It counts polylines whose cone is below PRE, not intersections.
3. **D's "coastal" is not measured.** The distance to the sea is tautological at a below-sea pit. The
   distance to the ocean-connected sea, or to PRE's coastline, is the missing instrument.
4. **Contribution vs class R8.** The contribution projects on each field's OWN total direction. B1a's
   interfluve Δ (+0.0038 of +0.0045) is a projection change with an unchanged class R8 (0.0120 against
   0.0122): B1 does not touch interfluve cells, so it can only come from mixed windows.
5. **Mixed windows.** A window is split fractionally, so "interfluve" R8 includes windows that also
   hold walls. B2's interfluve +65 % may be that.
6. **The skeleton's band is read at the chord's mid cell.** The polyline is sub-cell, so the mid cell
   can be off the trunk. The "sub < 1" row of a 1 km² skeleton (19 044 chords) is that artefact, and
   the other rows carry some of it.
7. **No D8 baseline** for the 8-point chords' R8.
8. **The exact −1.00 m** at both F38 pits is unexplained (first grep: NOTHING FOUND).
9. **The teeth are located by mouth class only**, not by position along the wall (foot / straight
   section / rim). That position is what decides between the reviewer's two B1 mechanisms.
10. **The ×3–4 lake-body count of B1** is not located.

## Uncommitted state (awaiting the author)

- `cached_product.rs`: the serde key plus its permanent test.
- `bench_field_hashes.json`: four refreshed digests.
- `workspace.rs`: the six-state `f123_viz_guard`.
- `tests/f125_walls.rs`: this bench.
- `f125_walls/*`: this report.

**No remedy, no primitive modified, nothing promoted.**
