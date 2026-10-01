# Finding 133 — the basin lakes are present lakes: the F132 residual is 93 % theirs, the extended base closes 95 % of it; θ falls by mixture AND by the laid law; the new lakes are the input lakes themselves

**The author's decision, recorded as given (2026-09-30):**
> « Les lacs des bassins sous la mer sont des lacs présents. »

The rule of 2026-09-29 therefore applies without exception: a present lake is a base level for what flows into
it, whatever the altitude of its bowl.

**Nothing is promoted.** `LakeBase::InputLakesAndBasins` is a new gated variant. F132's `InputLakes` stays
selectable, and `lake_base` is `None` by default.

- The benches are `f133_*` in `crates/ymir-core/tests/f126_coast.rs`.
- The raw outputs are beside this report.
- mur ↔ mer is ON, and the hydrology is D8.

**Units.** 1 cell = 48.8 m (domain), 0.00238 km², 8192² over 400 km. Heights are in m, volumes in km³.

**Reading declaration, unfavourable.** My predictions (`f133_predictions.md`) were written before any grep or
measurement. They are non-blind on Findings 73–132 and on the reviewer's predictions. Declared before each run:
- 0.5's answer before any measurement;
- R before R;
- F's design after R2 and before its code;
- T / L / E / consistency / cost before the worlds.

## 0 — prerequisites

- **Finding 132 was committed** as `1b2a519`. Three lines were added to its ADR entry, as asked: A refuted,
  "the circularity: a declared PROXY", and "the residual: one invariant, two paths".
- **Before the round** (`checks_before.txt`): the guard reads 6 / 6 "= banc"; the hash is
  `a8d2d538d692c2f0`; the lib passes 589 tests, 0 failed; check clean. **After it**: see "State".
- **The greps (`grep_rule11.txt`)**, earliest hit read:
  - below-sea: Finding 2;
  - `cuvette` / bowl: Finding 33 ("a deep bowl (floor −20 m)"), the below-sea class-2 cells;
  - **"walk to the col": Finding 38**, "fill each enclosed below-sea region as ONE water body … to its col" —
    the same idea, forty findings earlier;
  - `depression_floor`: Finding 96;
  - `basin_base` / F122-B: Finding 122;
  - `lake_base` / `InputLakes`: Finding 132.

### 0.5 — the code question

- **They are not two functions but two branches of ONE function**, the χ walk of `valley_construction::skeleton`:
  - the **main walk** applies the lake test at every cell it visits (under `lake_base`);
  - the **`basin_base` branch** (`if land[c] && dep(c)`) walks to the col, pushing every cell, with **no lake
    test**.
- **Every sibling path that integrates χ or lays a base:**

| path | stops at a present lake's shore? |
|---|---|
| `skeleton`, the main χ walk | **yes** (under `lake_base`) |
| `skeleton`, the `basin_base` col walk | **no** (until this round's extended variant) |
| `line_samples` (zf = base + kχ), `traced_polylines`, `Skeleton::floor_m` | they read the walk's base and χ, and inherit it |
| stream power A1 `depression_floor` (`stream_power.rs:703`, ON in the light pass with A1+B2) | no incision where `filled > field` on its own pit fill: it holds every depression at its spill, by construction, not by a lake test |
| stream power B3 `incision_floor`, Finding 96's χ floor (`stream_power.rs:709`) | χ from the sea, no lake test; `None` in production |
| B2 `slope_floor_uk` | a local slope; lays no base |

## R — reconciliation (measurements only, before any code)

**R1 — one lake per row** (`r_reconciliation.txt`, both worlds).

| totals | témoin | B2 → A_c |
|---|---|---|
| input lakes (all carry a base under `InputLakes`) | 25 | 25 |
| … whose flats move χ under OFF | **23** | 25 |
| … of those, present in the OFF world | **8** | 10 |
| … that lie in a `basin_base` closed depression (all at 100 %) | 7 | 7 |
| bases on a lake absent at the end (ON F132) | **4** | 4 |
| final lakes without a base (ON F132) | **19** (all below-sea basins) | 23 (18 below-sea + 5 land lakes) |

- **The témoin's totals fall back on 23 / 8 / 25 / 4 / 19.**
- The table has two kinds of row: the 25 input lakes, and the ON world's unbased final lakes.
- The 7 input lakes inside a closed depression are 5, 6, 17, 19, 20, 22 and 25. Each lies under a below-sea
  basin's final lake (1000011, 1000016, 1000017, 1000014, 1000019).

**R2 — the overlap.**
- The source of each residual cell was found by walking its D8 path to where χ and base agree again.
- **218 857 of the 235 602 residual cells (92.9 %) have one of the 19 unbased final lakes as their source.**
  Most come through 1000011 (input lakes 5 and 6: 196 938) and 1000017 (19 145).
- On B2 → A_c: **96.6 %**.
- The rest: 16 217 cells are laid by a moved sample with no differing cell on their path; 528 have an
  origin in no lake.
- **The stop rule (< 50 %) does not fire.**

## F — the basin lakes as present lakes (gated)

**What was built.**
- `LakeBase::InputLakesAndBasins`: a `basin_base` closed depression IS a present lake, filled to its col.
- χ stops at its **shore**, the first of its cells a path enters: χ = 0 there, and base = the col. The walk to
  the col, where F132's residual was born, no longer happens.
- **The level is the col**: the open-ocean flood's spill on the breached input, known when χ integrates.
  - A basin lake held below its col by an endorheic balance would be known only after the water balance on
    the delivered field. **That is circular; it is named and not built.**
  - In this climate every a_eq is infinite, so every basin is full.
- **The permanent test** `a_below_sea_basin_lake_stops_chi_at_its_shore`: an inland bowl 165 m below the sea.
  - Negative control: under `InputLakes`, the col walk gives the bowl's land slope χ > 0.
  - Extended: that slope has χ = 0, the upstream cell counts χ from the shore, and the base is the col in both.
  - Outside the bowl's catchment, nothing moves, bit for bit.
- The serde key: the same `Option<LakeBase>`, skipped when `None`; the key test holds. `None` is
  bit-identical.
- **The viz toggle** now takes the extended variant: "Lac = niveau de base (F133)", with both decisions in
  its tooltip. The F132 tooltip had been mangled by a heredoc (its line continuations were eaten); it was
  rewritten.

**The remedy test** (`f_remedy_test.txt`): A vs B outside the present lakes (the input lakes ∪ the closed
depressions); the identity control is 0 everywhere.

| | ON F132 | **ON extended** |
|---|---|---|
| témoin | 235 602 cells · p90 5.14 m | **12 509** · p90 68.8 m (max 631 m) |
| B2 → A_c | 481 606 · p90 11.72 m | **13 899** · p90 58.3 m |

- **It falls by 94.7 % / 97.1 %, not to 0.** The residual is located, not corrected:
  - **12 072 cells (témoin) have NO differing χ or base on their path.** They are laid by a polyline sample
    that moved. The lines crossing a lake follow the patched pointers, so their smoothed positions and their
    nearest-sample laying move near the shores.
  - **437 cells are born at origins outside every lake and depression**, with the same base and a different
    χ (e.g. 6.2 against 218.4). This is the signature of a lake with several spill cells: the two resolutions
    send the outflow through different exits, and the downstream drained area changes.

**Consistency** (the extended ON world):
- Bases on a lake absent at the end: **5**. These are F132's four (inputs 8, 12, 14, 21) and closed
  depression 3 (109 776 cells).
- **Final lakes without a base: 19 → 7.** They are 1000001 (1 km²), 1000006, 1000007 and 1000015 (0 km²,
  degenerate), 1000010 (45 km², 3 % covered), 1000002 (460 km², 38 %) and 1000014 (628 km², 37 %).
- The last two are basin lakes whose final extent exceeds the input's closed depression. They are larger at
  the end than their col predicts on the input.

## T — θ: the law or the mixture? (`t_l_e_table_worlds.txt`)

| | OFF | ON F132 | ON extended |
|---|---|---|---|
| θ (i) as built [95 %] | 0.338 [0.331, 0.344] | 0.243 [0.234, 0.251] | 0.246 [0.238, 0.255] |
| θ (ii) + breach | 0.369 [0.363, 0.376] | 0.297 [0.290, 0.306] | 0.301 [0.293, 0.309] |
| θ (iii) the world | 0.351 [0.344, 0.358] | **0.245** [0.236, 0.253] | **0.243** [0.235, 0.252] |
| **θ on the links LOWERED at both ends** | **0.493** [0.492, 0.495] | **0.436** [0.427, 0.445] | **0.439** [0.431, 0.446] |
| **trunk cells uncarved under their own law** | 32.1 % | **49.6 %** | 45.6 % |
| θ on those links | 0.160 | 0.126 | 0.106 |

**Both.** The fall 0.351 → 0.245 comes from the mixture **and** from the laid law.
- **The mixture**: the share of trunk cells left under their law (at θ ≈ 0.1–0.16) grows by +17.5 points
  (+13.5 extended).
- **The laid law**: θ on the carved links falls to 0.436 [0.427, 0.445], outside OFF's CI.
  - Along a path the law keeps S = kA^−0.5. What changes is WHICH links stay carved (34 786 against 46 885):
    the lake base lifts the floor of the lowland reaches above the terrain, and those leave the carved set.
  - The links that cross a lake's entry, where the base jumps, also enter the fit.
- The split between these two is NOT measured further.

## L — the new lakes (25 → 33)

**11 lakes of the ON F132 world are new against OFF** (8 more lakes net: 3 OFF lakes do not persist). **All 11
are themselves input lakes** (≥ 50 % of their cells), and all appear **at the construction stage**: a lake of
the as-built field's pre-drainage, with the same area there (e.g. 249.1 km² new lake 10, 163.1 km² new lake 11,
160.6 km² new lake 3).

- They are 5.3–249 km², 82–434 m deep; a_eq is infinite except new lake 1 (1 132 km²).
- **They are not "upstream of a based lake". They ARE the based lakes**, which OFF's construction drained. OFF
  cut through their sills, since its χ ran through them to the sea.
- The "meets an input lake at 49 m" of the raw output is the lake's own input footprint, which is larger than
  the final lake. It is one cell from the outlet.
- **ON extended adds none**: the same 11.

## E — rule 14, with denominators

| | removal withheld | share of the témoin's total removal (5 321.2 km³) | share of its removal in / upstream of the present lakes |
|---|---|---|---|
| **ON F132** | **890.0 km³** | **16.7 %** | 50.0 % (of 1 779.3 km³) |
| **ON extended** | **805.2 km³** | **15.1 %** | 28.7 % (of 2 803.5 km³) |

- **Where (ON F132)**: 620.7 km³ upstream of an input lake, 253.1 inside one, 14.2 in or upstream of a closed
  depression, 1.9 elsewhere.
- **Where (ON extended)**: 620.4, 252.8, **−70.0** (MORE removed: χ from the shore is shorter than χ from the
  col, so the floors upstream of the basins are laid lower), and 2.0.
- F132's "648.5 km³" was ON − OFF over the cells ON left higher, upstream of the lakes. Here the removal is
  measured against the construction's input S1. The two agree in order (620.7 km³ upstream).

**Cost** (wall-clock, witness):
- the construction stage alone (skeleton + carve) takes **120–133 s**;
- a full second pass (the world's build + drainage, breach, climate, assembly) takes **215–247 s**.

This is the price of the fixed point that is named, not built: one more such pass per iteration.

## The table against the témoin (rule 18: the same world, OFF / ON)

| | OFF | ON F132 | ON extended |
|---|---|---|---|
| canyons | 0 | 0 | 0 |
| coast (spurs near a coastal wall) | 0.0142 /km | 0.0142 | 0.0142 |
| θ (i) / (ii) / (iii) | 0.338 / 0.369 / 0.351 | 0.243 / 0.297 / 0.245 | 0.246 / 0.301 / 0.243 |
| R8 terrain / network | 0.0452 / 0.5613 | 0.0395 / 0.5757 | 0.0411 / 0.5658 |
| lakes | 25 (3 492 km²) | 33 (4 453 km²) | 33 (4 445 km²) |
| relief, paired p50 | 515.4 m | 551.9 m | 546.6 m |
| removal (S1 − world) | 5 321.2 km³ | 4 431.2 km³ | 4 515.9 km³ |
| Δz between two arbitrary flat resolutions, outside the present lakes (cells · p90) | 1 150 588 · 14.45 m | 235 602 · 5.14 m | **12 509** · 68.8 m |

No admissibility gate is declared (rule 15). The table is reported against the témoin, and nothing is asked of
the author's eye before it.

## Predictions scored

**Mine** (`f133_predictions.md`):
- **0.5 HELD**: two branches of one function, and the siblings as listed.
- **P-R1 HELD** (the totals fall back; two kinds of row). **P-R2 HELD** (92.9 %, 96.6 %).
- **P-F HALF**:
  - the residual < 1 000 cells: refuted (12 509);
  - the unbased final lakes 19 → 0: refuted (19 → 7);
  - the bases on absent lakes 4 → 4–6: held (5).
- **P-T HALF**:
  - θ on the carved links 0.47–0.51: refuted (0.436);
  - the under-law share +≥ 10 points: held (+17.5);
  - "the mixture": half (the mixture AND the law).
- **P-L HELD**: the new lakes are input lakes (11 of 11), at the construction stage; the extended adds 0.
- **P-E REFUTED**: 30–70 % of the total removal (it is 16.7 %).
- **Meta HELD.**

**The reviewer's:**
- **R2** (> 70 %, no stop): HELD.
- **F**:
  - the residual < 5 000 and not 0: refuted (12 509);
  - the unbased lakes 19 → 0: refuted (7);
  - the absent bases stay 4: refuted (5).
- **T**:
  - θ on the carved cells 0.49 ± 0.02: refuted (0.436);
  - the under-law share +≥ 5 points: held (+17.5);
  - "so the mixture, not the law": refuted (both).
- **L**: "≥ 5 of the 8 new lakes upstream of a based lake": refuted (they are the based lakes); "the extended
  adds ≥ 2 more": refuted (0).
- **E** ("< 10 % of the total"): refuted (16.7 %).
- **Meta HELD.**

## Limitations, stated

1. **The residual's two remaining sources** (moved samples; multi-exit lakes) are read from the origin walk and
   the χ values. Neither was patched to confirm (the round says: locate, do not correct).
2. **"Present" is still the input's lakes and the input's closed depressions, filled to their col**: a
   declared PROXY. The fixed point is costed (one full pass, ~215–247 s), not built.
3. **The θ split between mixture and laid law** is shown by two numbers (the under-law share, θ on the carved
   links), not decomposed further.
4. **L's "upstream of a based lake"** walked from the lake's outlet and met the lake's own input footprint. The
   verdict ("they are the based lakes") rests on the ≥ 50 % input coverage of every new lake.
5. **Timings are wall-clock**, on a machine that has suspended during long benches before.

## Viz concordance and the author's look (F133v, added at commit)

- The viz path (`run_hd`, `f133v_viz_path.txt`, `f133v_passage.txt`) lists **26 / 34** lakes against the bench
  assembly's 25 / 33. The difference is the crater lake 2000001 (CraterAcidic, 1.19 km², level 2 104.9 m), added by
  the C-2 crater pass (`hd.rs:1128-1158`) in both states. The lake counts and areas of Findings 131–133 are
  assembly figures (+1 lake, +1.19 km² in the author's world); canyons, θ, R8, relief, volumes and Δz are not
  concerned.
- The author's look (2026-09-30), viz toggle, extended variant:
  - the gate holds on the crops of lakes 1 and 11 (362 m for 11);
  - OFF, lake 11 has another shape and level, and lakes 1 and 2 are absent;
  - no regression seen; the added lakes are acceptable.
  - **Not covered**: lake 2's gate, the eight other crops, Living Landz. **Not promotable as it stands.**
- The author's direction: « Je veux réduire la taille des lacs avec l'âge, car au départ ils étaient tous remplis à
  la hauteur d'équilibre, ce qui donnait des lacs énormes. » It is the entry criterion of H-2 / age (not opened).

## State

**Committed with the go-ahead of the F134 round** (Partie 0.1), together with the `f133v_*` viz benches
(`workspace.rs`, `#[cfg(test)] mod f133v_bench`).
- `valley_construction.rs`: `LakeBase::InputLakesAndBasins`, the shore stop in the `basin_base` branch, the new
  permanent test, and `ocean_flood` made `pub` (visibility only).
- `ymir-viz/src/ui/workspace.rs`: the toggle takes the extended variant, and its tooltip is rewritten.
- `f126_coast.rs`: the `f133_*` benches.
- This report, and ADR Finding 133 (with the author's decision, dated).

**Checks after the last change** (`checks_after.txt`):
- The guard reads 6 / 6 "= banc", and the hash is `a8d2d538d692c2f0`.
- The lib passes 590 tests, 0 failed.
- `cargo check --workspace --tests --release` is clean.
