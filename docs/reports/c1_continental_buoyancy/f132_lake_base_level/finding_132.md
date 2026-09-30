# Finding 132 — a present lake is a base level: the hidden parameter measured, then the author's decision built (gated), its residual located in the closed-depression walk

**The author's decision, recorded as given (2026-09-29):**
> « Un lac présent est un niveau de base pour les rivières qui s'y jettent ; un lac vidé ne l'est plus. »

The base follows the lake's existence, not the state of the breach. `breach_monotone` conditions the drainage; it
does not empty a lake.

**Nothing is promoted.** `ValleyConstruction::lake_base` is gated and `None` by default. The concordant clause
stays gated and off.

- The benches are `f132_*` in `crates/ymir-core/tests/f126_coast.rs`.
- The raw outputs are beside this report.
- mur ↔ mer is ON; bruit ↔ pied is ON where there is a wall profile.
- The hydrology is D8 everywhere.

**Units.** 1 cell = 48.8 m (domain), 8192² over 400 km. Heights are in m, volumes in km³, and W is the larger
line's width at a junction.

**Reading declaration, unfavourable.** My predictions (`f132_predictions.md`) were written before any grep or
measurement. They are non-blind on Findings 73–131 and on the reviewer's predictions. Each block's definitions
were declared before its run:
- 0.5's answer before any measurement;
- P1–P3 before P;
- A before A;
- P4's design before its code;
- the residual's localisation after the remedy's counts, and before its own run.

## 0 — prerequisites

- **Finding 131 was committed** as `39d07dd`. Two lines were added to its ADR entry, as asked:
  - the reviewer's Finding 130 reading ("P without object") is refuted;
  - B's gate failed, LTD is not adopted, C is blocked, and the gate is not rewritten.
- **Before the round** (`checks_before.txt`): the guard reads 6 / 6 "= banc"; the definition's hash is
  `a8d2d538d692c2f0`; the lib passes 588 tests, 0 failed; `cargo check --workspace --tests --release` is clean.
  **After it**, see "State".
- **The greps (`grep_rule11.txt`)**:
  - **"lac présent / lake present": NOTHING FOUND**, in the ADR and in the code.
  - `local base level` / `local_base`: NOTHING FOUND in the code. In the ADR the phrase appears first in
    Finding 4 ("pool to local base levels bounded by contours": the droplet terraces).
  - `basin_base` first appears in Finding 122-B, read. χ counts from the col of a closed depression of an
    open-ocean priority flood on the BREACHED field, so a breached lake is not closed there.
  - `a_eq` first appears in Finding 29; `wc == 2` in Finding 38; `breach` in Finding 12; `vidange` in Finding
    114 ("lakes emptied by one pass").

### 0.5 — the code question, answered before any measurement

- **The stages** in `run_hd` (`ymir-viz/src/bridge/c1/hd.rs`):
  1. `cached_c1_eroded` = `upscale_from_c1`, which holds the construction (`production_upscale.rs:487-497`) and
     the light pass;
  2. the pre-breach drainage of the **delivered** field;
  3. the breach;
  4. the climate;
  5. `build_hd_drainage`.

  **`lakes.json` carries the pre-breach lakes of the delivered field**, typed by the water balance.
- **When χ integrates** (`skeleton`, inside the construction), the lakes known are:
  - the construction INPUT's pre-drainage depressions (`c1_drainage_windowed`'s lake map: topographic, no
    climate);
  - `basin_base`'s closed depressions of the open-ocean flood on the breached input.
- **The circularity, named.** A lake "present" in `lakes.json` exists only after the construction and the
  light pass that it conditions.
  - P4 does not decide it in silence. Its declared choice is `LakeBase::InputLakes`: the input's lakes count
    as present.
  - The fixed point (a second pass on a first pass's `lakes.json`) is named and not built, because it needs
    a lake set from outside the construction's config.

## P — the hidden parameter (measurements, no remedy)

**P1 — its own amplitude** (`p1_p2_p3_lakes_flats.txt`).
- Two different seeded arbitrary acyclic resolutions of the pointers on the lakes' 684 790 flat cells
  differ on 508 804 of them. LTD is the third state.
- The identity patch changes **0** cells in both worlds.

| | témoin C2/10 col | B2 → A_c |
|---|---|---|
| **arbitrary A vs B, outside the lakes** | **1 150 588** cells · \|Δz\| p50 2.64 / **p90 14.45** / max 742.8 m · up to 20.0 km from a lake | **2 645 323** · p50 3.84 / **p90 32.14** / max 727.5 m · 20.0 km |
| LTD vs A / LTD vs B | p90 31.80 / 32.18 m | p90 47.31 / 47.46 m |

- **The stop rule (p90 < 1 m) does not fire.** The effect Finding 131 measured is not LTD's own: two
  arbitrary resolutions differ by the same order.

**P2 — hydrological status**, against each world's delivered lakes (the assembly's pre-breach lakes, with
`classify_lakes_water_balance`).

| | lakes whose flats move the construction | **present** | absent | changed cells attributed to no lake |
|---|---|---|---|---|
| témoin | 23 | **8** | 15 | 20 426 |
| B2 → A_c | 25 | **10** | 15 | 31 776 |

- The per-lake table (area, depth, level, type, a_eq, cells moved) is in the raw output.
- Present lakes include the two below-sea basins' land lakes (pre-lakes 5, 6 → 1000011; 19, 20 → 1000017)
  and lakes 16, 18, 22, 24.
- **Every a_eq is infinite**: the balance is in surplus in this climate.
- **No absent lake is water class 2.** 12 of the 15 absent lakes lie 0–10 % under any final lake: they are
  gone.
- **The stop rule (no lake present) does not fire.**

**P3 — overlap with A** (B2 → A_c). The zone P is the cells whose z differs between A and B: 3 001 152 cells.
- **49.9 % of Finding 131's 30 768 junctions lie in the zone.**
- **col 7 is IN it** (\|Δz\| 3.45 m, path into pre-lake 5), and so is **13's dam cell** (0.68 m, pre-lake 16).
- **col 14 and canyon 14's floor are OUTSIDE** (0.00 m, no lake on their path).
- So A's canyon 14 reading is noise-free. A's junction-level figures (replat, R8 around the junctions) sit on
  a population half inside the zone (see A).

## P4 — "a present lake is a base level", gated

P1 and P2 did not fire, so P4 ran.

**What was built** (gated; `None` is bit-identical; see "State"):
- `ValleyConstruction::lake_base: Option<LakeBase>`, with `LakeBase::InputLakes`.
- **The rule**: in the skeleton's χ walk, a path entering an input lake's footprint stops there. The lake's
  cells take χ = 0 and base = the lake's surface (its col), and upstream cells integrate χ from the lake.
- **The cache key**: skipped in serde when `None`, listed in `a_gated_off_construction_field_does_not_move_the_eroded_key`.
- **The permanent test** `a_present_lake_is_the_base_level_of_its_catchment`: a 150 m bowl in a plane falling
  to the sea.
  - Negative control: without the field, the upstream cell takes the sea's base.
  - With it, that cell takes the lake's surface, and its χ counts from the lake.
  - Downstream of the lake, nothing moves, bit for bit.
- **The viz toggle** "Lac = niveau de base (F132)" is under the valley mode. Its hover text states the
  decision and the circularity, and the guard badge reads "non gardé".

**The remedy test** (`p4_remedy_test.txt`): the two arbitrary resolutions, lake_base ON, compared outside the
lakes.

| | OFF (as P1) | **ON** |
|---|---|---|
| témoin | 1 150 588 cells · p90 14.45 m | **235 602** · p90 **5.14** m (max 631 m) |
| B2 → A_c | 2 645 323 · p90 32.14 m | **481 607** · p90 **11.72** m (max 1 116 m) |

- **Δz does NOT fall to 0 outside the lakes.** ON removes 79.5 % / 81.8 % of the cells and two thirds of the
  p90. The identity control holds at 0.

**The residual, located (not corrected)** (`p4_residual_origin.txt`, témoin, ON). Each residual cell's path
was walked downstream to where χ becomes equal in A and B.
- **No pointer differs outside the lakes.** 508 804 differ inside, as patched.
- **759 717 residual cells stand behind 15 origins that lie INSIDE an input lake**, with χ ≠ 0 there:
  - (4401,3623) in lake 5, χ 346, base 459.3 m;
  - (4709,3609) in lake 6, χ 898, base 459.3 m;
  - lakes 17, 19 and 25, with bases 115.3, 268.8 and 614.9 m.
- 149 032 stand behind 3 "other" origins, with bases 49.3 and 268.8 m.
- 2 363 stand behind 24 origins at a lake's entry.
- **The bases are the levels of the closed depressions of Finding 122-B's ocean flood**: 459.3 m is
  1000011's, and 268.8, 49.3, 614.9 and 115.3 m are the other below-sea basins'.
- **The mechanism.** The `basin_base` branch of the χ walk (`if land[c] && dep(c)`) walks **to the col**,
  pushing every cell on the way, **without the lake test**. A lake lying inside such a closed depression is
  therefore crossed by that walk. Its cells get χ ≠ 0, and the path across its flat re-enters χ upstream.
- **The intra-lake path enters through the closed-depression walk.** It is located and not corrected, as
  asked.

**Consistency (the fixed point, reported, no remedy)**, on the ON world:
- **Bases on a lake absent at the end: 4**, out of 25 input lakes. They are 8, 12, 14 (0 % of their cells
  under a final lake) and 21 (20 %).
- **Final lakes without a base: 19. All are below-sea basins** (ids ≥ 1 000 000): not land lakes of the
  input, and handled by `basin_base`'s col.
- **No land lake is present at the end without a base.**

**Rule 18's pair on the témoin C2/10 col** (`p4_temoin_off_on.txt`):

| | OFF | **ON** |
|---|---|---|
| canyons (over-dug) | 0 | 0 |
| spurs near a coastal wall | 0.0142 /km | 0.0142 /km |
| θ (i) as built [95 %] | 0.338 [0.331, 0.344] | **0.243** [0.234, 0.251] |
| θ (ii) + breach | 0.369 [0.363, 0.376] | **0.297** [0.290, 0.306] |
| θ (iii) the world | 0.351 [0.344, 0.358] | **0.245** [0.236, 0.253] |
| R8 terrain / R8 network | 0.0452 / 0.5613 | 0.0395 / 0.5757 |
| lakes | 25 (3 492 km²) | **33 (4 453 km²)** |
| relief, paired p50 | 515.4 m (median cut 5.7 m) | **551.9 m (+7.1 %)** (median cut 0.6 m) |

**Rule 14 — the erosion ON − OFF, and where.**

| where | cells changed | removed LESS (ON higher) | removed MORE |
|---|---|---|---|
| **upstream of an input lake** | 1 731 586 | **648.5 km³** | 14.5 km³ |
| inside an input lake | 389 882 | 253.2 km³ | 0.02 km³ |
| elsewhere | 407 576 | 5.1 km³ | 2.2 km³ |

- **A lake's base level bounds the incision upstream**: 648.5 km³ less is removed above the lakes. The floors
  there start at the lake's surface instead of integrating down to the sea.
- The lakes gain 8 bodies and 961 km², and the relief p50 rises 7.1 %.
- **θ falls from 0.351 to 0.245.** Where the lake base lifts the law above the terrain, the construction
  leaves the field (Finding 129's uncarved third grows).
- Rule 15: no admissibility gate is declared. The table is read against the témoin, and any later threshold
  will be derived from it with its provenance.

## A — rule 16 on Finding 131's concordant clause (a read-only bench)

**Populations**: (ii) is the junctions whose tributary has a moved cell among its first 4 cells; (i) is the
others.
- B2 → A_c: (i) 26 717, (ii) 4 051 junctions. B2 + B1 + foot: (i) 26 776, (ii) 4 117.
- **On the construction, each population's field is composed exactly** (z0 outside its clause cells, the
  concordant values on them): carve decides every cell alone.

| B2 → A_c | raised > 1 m | θ laid [95 %] | replat: concordant / none (junction share) | replat length p50 | junctions in zone P |
|---|---|---|---|---|---|
| (i) alone | 64 | 0.496 [0.493, 0.499] | **57.3 %** / 16.6 % | 1.14 W | 48.9 % |
| (ii) alone | 72 | 0.493 [0.491, 0.496] | **74.1 %** / 21.3 % | 1.50 W | 56.8 % |

| B2 (A_c) + B1 + foot | raised > 1 m | θ laid [95 %] | replat: concordant / none | length p50 | in zone P |
|---|---|---|---|---|---|
| (i) alone | 46 | 0.494 [0.491, 0.497] | 52.3 % / 12.5 % | 1.24 W | 48.7 % |
| (ii) alone | 52 | 0.492 [0.489, 0.495] | 67.9 % / 17.6 % | 1.68 W | 56.5 % |

| R8 within 2 W of the junctions (pipeline worlds) | (i) terrain / network | (ii) terrain / network |
|---|---|---|
| B2 → A_c, none → concordant | 0.0310 → **0.0310** / 0.2264 → 0.2077 | 0.0415 → **0.0416** / 0.2796 → 0.2323 |
| B2 + B1 + foot, none → concordant | 0.0268 → **0.0342** (+28 %) / 0.2219 → 0.1905 | 0.0312 → **0.0408** (+31 %) / 0.2724 → 0.2169 |

- **The replat's rise does NOT concentrate on (ii).**
  - (ii)'s rate is higher (74 % against 57 %), but (ii) holds only 16.4 % of the excess junctions (2 141 of
    13 035 on B2 → A_c; 16.3 % on B1).
  - The replat is a property of the clause at ordinary junctions (i).
- **The terrain R8 rise does not concentrate on (ii) either.**
  - On B2 → A_c neither population rises.
  - On B2 + B1 + foot both rise alike (+28 %, +31 %).
  - The network R8 falls in all four.
- **Canyon 14, attributed.** In both concordant worlds the témoin's path from canyon 14's floor is DAMMED 28
  cells down at **(1807,4580)**, two cells from Finding 128's col 14. It sits at 384.0 m on B2 → A_c and
  447.1 m on B1, against the témoin's 120.1 m.
  - **It is a clause cell, of population (i)**: a tributary's last cells laid at the junction's floor, which
    here stands above the canyon's water.
  - **It is outside zone P**, so the reading is under no noise from the lakes' flats.
- **Under P1's noise.** About half of each population's junctions lie in zone P (48.9 % / 56.8 %). The
  junction-level differences of A (replat shares, R8 around the junctions) are therefore read on sets half
  exposed to a lake-flat noise of p90 14–32 m. The noise is identical in the none and concordant worlds (the
  same flats, the same D8), so the comparisons stand, but each population's absolute figures do not.
  Canyon 14, col 14 and θ laid are outside that concern (canyon 14 and col 14 out of the zone; θ with CI
  over ~42 000 links).

## Predictions scored

**Mine** (`f132_predictions.md`):
- **0.5 HELD**: the circularity, as named.
- **P**
  - **P-P1 HALF**: p90 8–20 m held on the témoin (14.45; B2 above at 32.1); p50 3–8 m refuted narrowly
    (2.64); the count, the reach, the max and the control held.
  - **P-P2 REFUTED**: 8 of 23 and 10 of 25 are present, not ≥ 2/3.
  - **P-P3 HALF**: col 14 outside held; 15–35 % of the junctions refuted (49.9 %).
- **P4**
  - **P-P4 HALF**:
    - a residual: held in sign, but it comes from the closed-depression walk, not from the "no lake" cells;
      "≥ 98 % vanish" refuted (80 %);
    - ≥ 10 bases on absent lakes: refuted (4);
    - ≥ 1 final lake without a base: held (19, all below-sea basins);
    - canyons ±1: held (0 / 0);
    - relief ±3 %: refuted (+7.1 %);
    - less erosion upstream of the lakes: held (648.5 km³).
- **P-A HALF**:
  - "< 15 % of the replat excess on (ii)": refuted narrowly (16.4 %);
  - "the terrain R8 rise mostly on (i)": refuted (none on B2; alike on B1);
  - "canyon 14's dam on a touched junction": HELD (a population-(i) clause cell).
- **Meta ("at least two wrong") HELD.**

**The reviewer's:**
- **P1** ("p90 5–15 m, reach > 10 km"): HELD on the témoin (14.45 m, 20 km), above it on B2 (32.1 m).
- **P2** ("10–20 present"): REFUTED on the témoin (8); HELD at its edge on B2 (10).
- **P3**: "> 30 % of the junctions in the zone" HELD (49.9 %); "col 14 in the zone" REFUTED.
- **P4**:
  - "Δz = 0 outside the present lakes, to ≤ 1 cell": REFUTED (235 602 cells);
  - "canyons unchanged": HELD;
  - "relief ±3 %": REFUTED (+7.1 %);
  - "≥ 1 base on a lake absent at the end": HELD (4).
- **A**: "> 60 % of the replat excess on (ii)" REFUTED (16.4 %); "canyon 14 dammed outside the touched
  junctions" REFUTED.
- **Meta HELD.**

## Limitations, stated

1. **"Present" at construction time is the input's lakes** (the declared choice). The fixed point on a first
   pass's `lakes.json` is not built. The consistency counts (4 bases on absent lakes; no land lake without a
   base) measure how far this choice is from it.
2. **The residual's mechanism is read from the origin cells and the code's order.** The origins lie in lakes
   inside the ocean flood's closed depressions, with bases equal to those depressions' spill. No patch was
   made to confirm it (the round says: locate, do not correct).
3. **P2's "present"** means ≥ 50 % of a pre-lake's flat cells in one final lake's footprint.
4. **A's population split** uses F131's source-cone census without its ≥ 11-cell cut. The R8 around the
   junctions uses 8 × 8 windows, compared only between worlds on the same union of disks.
5. **Rule 14's volumes** are ON − OFF of the delivered field over land, km³ at 0.00238 km² per cell.
6. **The first A and P4 runs died with the session** (exit 4, 00:46). Both were relaunched unchanged. A's
   B2 → A_c lines are identical in both runs.

## State

**Uncommitted**, awaiting the go-ahead.
- `valley_construction.rs`: `lake_base` / `LakeBase`, the χ walk's lake stop, and the new permanent test.
- `cached_product.rs`: the key test's list.
- `ymir-viz/src/ui/workspace.rs`: the toggle.
- `f126_coast.rs`: the `f132_*` benches.
- This report, and ADR Finding 132 (with the author's decision, dated).

**Checks after the last change** (`checks_after.txt`):
- The guard reads 6 / 6 "= banc", and the hash is `a8d2d538d692c2f0`.
- The lib passes 589 tests, 0 failed.
- `cargo check --workspace --tests --release` is clean.
