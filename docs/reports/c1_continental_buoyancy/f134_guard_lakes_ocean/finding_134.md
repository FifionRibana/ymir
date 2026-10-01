# Finding 134 — the guard covers the lakes; the ocean moves by a tenth of a millimetre; lake 2's 1 177 m is OFF's canyon through its col; the new under-law cells are the lake beds

Raw outputs are in this folder:
- `lake_guard_and_identity.txt`: the lake guard bench, the viz guard and the viz listings;
- `bench_lakes_{OFF,ON}.txt` and `viz_lakes_{OFF,ON}.txt`: the lake-by-lake listings;
- `f134_ocean.txt`, `f134_d.txt`, `f134_tu.txt`, `f134_delta.txt`: the four benches;
- `d_lake{2,11,1}_dz.png`: the Δz maps (legend: `maps_legend.md`: sign, scale, saturation, outlines);
- `checks_before.txt` and `checks_after.txt`: the checks around the change.

The predictions were written before any grep or measurement (`f134_predictions.md`). The instruments were
declared before the measurements (`f134_declared.md`).

**Units**: 1 cell = 48.8 m (domain), 0.00238 km².

**Two worlds**:
- **The viz world** is `new(k, Some(0.1))`, the guarded "C2 /10 col". Parts O and D are measured on it.
- **The benches' témoin** is the same plus `wall_sea_floor_m = 0.5` (F126-B's mur ↔ mer). The viz never sets
  it. Parts T, U and Δ are measured on it, where F133's numbers were read.

## Partie 0

**0.1 — F133 committed** (`1a5f3ff`), with the `f133v_*` benches and the three ADR additions.
- Before the commit: guard 6 / 6, check clean, lib 590 passed.

**0.2 — The guard covers the lakes.**

*Code.*
- `bench_guard::check_lakes` is keyed by `hd_drainage_key` (eroded key + drainage config + climate).
- It fingerprints the final lakes **after the C-2 crater pass**: FNV-1a on `lake_map`, then on the list
  (id, type, level, area).
- The reference is `data/bench_lake_hashes.json`. One climate is guarded (45° / 40°); any other reads
  NoReference.
- `HdResult.lake_guard` is computed in `run_hd` after the crater pass. The microscope shows a second badge,
  and a lake Mismatch refuses numbers like a field Mismatch.
- The permanent test `the_lake_guard_sees_one_cell_one_level_one_type_one_area` reads Match, then Mismatch on
  each negative control: one mask cell, one level, one type, one area, one lake fewer. Another digest reads
  NoReference.

*References.* `f134_lake_guard` regenerated them from the bench (six states):

| state | lakes | crater lakes |
|---|---|---|
| livré | 60 | 2 |
| A1+B2 | 38 | 2 |
| C1 nue | 27 | 2 |
| C2 /10 col (défaut) | 26 | 1 |
| C2 /3 | 26 | 1 |
| C2 /10 niveau mer | 24 | 1 |

Before regenerating, each world's field was asserted equal to the field guard's reference.

**0.3 — The bench runs the crater pass.**
- The crater-protect mask and the crater-lake pass were **moved verbatim from `run_hd` to ymir-core**
  (`volcanism::crater_protect_mask`, `crater_lake_pass`). `run_hd` now calls them.
- The bench tail `common::viz_hd_lakes_on` calls the same functions as `run_hd`: the protected breach, the
  placed climate, the H-1 infiltration field, the assembly, then the crater pass.
- **Lake-by-lake identity holds, OFF and ON.** `lake_listing` gives id, type, level and area with their bits,
  the cell count, and the hash of the cell indices.
  - The bench and `run_hd` listings are **identical: 26 / 26 lines OFF, 34 / 34 ON**, the crater lake
    2000001 included.
  - The ON fingerprint is the same on both sides (`d7a888b0a90d795a/93e3bb468a9513a6`, key `64e3c3fddfdf0dd4`).
- No correction was needed beyond the move. The cost is one HD tail per world, 140–260 s.
- Which of the two missing pieces (the protected breach, the crater pass) was needed alone is not tested:
  both went in together.

**0.4 — Checks.**
- Before: guard 6 / 6 field, check clean, lib 590.
- After: guard **6 / 6 field AND lakes**, hash `a8d2d538d692c2f0`, check clean, lib **591** (the new test),
  viz lib 30 passed.

## O — the ocean moves (viz world, OFF / ON extended)

**Grep (rule 11 / 11c).**
- The normalisation is `terrain/bathymetry.rs:188`: `mean_depth` over **every** cell `≤ sea`, inland
  below-sea basins included.
- It was introduced by `563198e` ("submarine bathymetry").
- Its earliest ADR hit is **Finding 124-P1** ("the bathymetry normalises by the ocean's mean depth, which is
  global (every ocean cell differs at S7)"). The "Bound added by Finding 124-P1" note in Finding 57 is a later
  insertion.
- With the shipped `texture = 1`, the re-map is `depth = env(dist) · c / m`. Where `c` is unchanged, the ratio
  is exactly `m_OFF / m_ON`.

**Stages (rule 12).** Each stage is built by the pipeline with a `Knobs` switch; OFF vs ON, bitwise:

| stage | cells differing | OCEAN cells (both sea, edge-connected) | inland sea | sea ↔ land flips |
|---|---|---|---|---|
| (a) construction (+ rims) | 1 940 794 | **0** of 55.4 M | 0 | 0 |
| (b) + light pass | 4 407 766 | **69** | 46 | 3 |
| (c) + droplet erosion | 4 407 766 | 69 | 46 | 3 |
| (d) + bathymetry | 44 067 569 | **39 659 918** | 0 | 3 |

- **The first OCEAN difference is at the light pass**: 69 cells, where the light pass reaches the shore.
- (c) adds nothing: its counts are identical to (b)'s.
- **The mass movement is born at the bathymetry.**
- The 44.3 M of F133v were read on the conditioned field in norm bits. In metres, f32 conversion merges part of
  them: 42 180 409 differ.

**Mean depth.** It is the bathymetry's own, on the stage (c) field.
- OFF: **2 566.6381 m** over 55 780 438 cells. ON: **2 566.6380 m** over 55 780 439 cells.
- m_OFF / m_ON = 1.0000000 at 7 digits.
- The sources: the 69 ocean and 46 inland cells, and one more cell below the sea.

**Ratio z_ON / z_OFF** on the 55 780 437 sea cells:
- median 1.0000001;
- p0.1–p99.9 all at 1.00000;
- **100.000 % within 10⁻³ of the median**.

**Attribution.**
- Of the 54 665 649 unclamped cells whose pre-bathymetry depth `c` is bitwise equal, **99.993 % sit at the
  prediction `m_OFF/m_ON` within 10⁻⁵**.
- No cell has a different `c`. 1 114 788 cells are clamped (shelf_min 1 m or floor 5 650 m).

**⇒ It is the normalisation. Stop rule: reported as a non-local stage, no fix.** The author has decided:
"Toute modification des terres déplace l'océan entier."

The order of the move: the ratio is 1 + O(10⁻⁸), i.e. **about 0.1 mm** at 2 500 m. The 39.7 M ocean cells
differ by the last bits of an f32, through the global mean.

## D — lakes 2, 11, 1 (viz world, reproduced by the bench, lake ids = F133v's)

**The field.** F133v's Δz is read on the **conditioned field** (`HdResult.eroded`). Its pre-breach lakes are at
their surface: the breach sets lake cells to their fill. The pre-breach eroded field (no water) is read
beside it.
- Upstream, the two fields agree to the metre: no OFF lake lies there.
- **Sign: Δz = z_ON − z_OFF > 0** on 96–100 % of the carved cells.
- F133v's populations are reproduced exactly: lake 2 has 88 142 upstream cells, p90 **1 177.4 m**.

**Conditioned field, by lake and population:**

| lake (ON level / depth / col) | population | cells (Δz ≠ 0) | Δz > 0 | p50 / p90 \|Δz\| | Σ\|Δz\|·A |
|---|---|---|---|---|---|
| **2** (1 292.0 / 356.3 / 1 339.6 m) | footprint, carved | 14 735 | 100 % | 942 / 1 108 m | 30.5 km³ |
| | footprint, kept | 154 | 100 % | 236 / 287 m | 0.1 |
| | upstream, carved OFF and ON | 11 882 | 100 % | 1 171 / **1 240 m** | 30.3 |
| | upstream, carved OFF only | 57 719 | 99.4 % | 796 / **1 182 m** | 102.1 |
| | upstream, kept at the terrain | 18 541 of 35 885 | 90.8 % | 90 / 409 m | 6.7 |
| **11** (362.0 / 329.1 / 390.5 m) | footprint, carved | 15 103 | 100 % | 205 / 259 m | 6.7 |
| | footprint, kept | 53 275 | 100 % | 186 / 266 m | 22.5 |
| | upstream, carved both | 69 906 | 99.9 % | 515 / 584 m | 66.3 |
| | upstream, carved OFF only | 67 232 | 98.7 % | 211 / 502 m | 39.7 |
| | upstream, kept | 82 335 of 186 776 | 82.8 % | 10 / 73 m | 5.0 |
| **1** (327.2 / 181.6 / 350.8 m) | footprint, carved | 28 662 | 100 % | 222 / 293 m | 14.3 |
| | footprint, kept | 27 056 | 100 % | 105 / 141 m | 6.3 |
| | upstream, carved both | 154 126 | 99.7 % | 243 / 300 m | 100.3 |
| | upstream, carved OFF only | 63 741 | 96.1 % | 219 / 1 062 m | 55.3 |
| | upstream, kept | 87 828 of 142 410 | 71.7 % | 18 / 116 m | 8.6 |

- No cell is carved by ON only.
- On the eroded field the footprint's "kept" cells fall to p90 2.4 m (lake 2), 68 m (lake 11) and 55 m
  (lake 1): there, the conditioned Δz is the ON lake's water surface.
- Lake 11's footprint lies over an OFF lake on 13 740 cells.

**Reading, lake 2.**
- **The 1 177 m is not in the lake's footprint.** It sits upstream, on cells that **OFF's construction
  carved**: 69 601 of 88 142 carved cells, 132 km³ of the 139.
- In OFF, the footprint's conditioned floor is p50 **352 m**; the col is **1 339.6 m**. OFF's χ, counted from
  the sea through the bowl, laid a valley **about 1 km under the col**: it cut through it. ON lays its base at
  the lake level, 1 292 m.
- Δz (p50 1 171 m on the cells carved in both) is about **three times the lake's depth** (356 m). It is of the
  order of the col above OFF's floor (1 339.6 − 352 ≈ 990 m).
- Lakes 11 and 1 show the same structure on a smaller scale: their upstream carved Δz (515 / 243 m) also
  exceeds their depth (329 / 182 m).
- The map (`d_lake2_dz.png`) shows a saturated red basin (> +1 000 m) north and east of the lake. Its trunks
  are where OFF drained through the col.

## T — θ on the carved links, by population (benches' témoin)

The totals are reproduced: OFF **0.493** [0.492, 0.495]; F132 **0.436** [0.427, 0.445]; extended **0.439**
[0.431, 0.446].

| variant (based set) | upstream of a based lake | downstream | untouched basins |
|---|---|---|---|
| OFF (F132 set) | 0.492 [0.488, 0.496] (13 538) | 0.500 [0.493, 0.507] (1 694) | 0.489 [0.485, 0.492] (19 020) |
| **ON F132** | **0.492** [0.488, 0.496] (13 548) | **0.500** [0.493, 0.507] (1 694) | **0.489** [0.486, 0.493] (19 019) |
| OFF (extended set) | 0.494 [0.492, 0.497] (22 222) | 0.503 [0.485, 0.522] (724) | 0.485 [0.480, 0.490] (11 313) |
| **ON extended** | **0.494** [0.492, 0.497] (24 597) | **0.512** [0.496, 0.530] (871) | **0.485** [0.480, 0.490] (11 470) |

**No population's θ moves.** Upstream of a based lake it is 0.492 / 0.494 in OFF and in ON.

**The fall 0.493 → 0.436 is in none of the three populations.** It comes from:
- the links that leave the carved set: OFF carves 46 885 links, F132 34 786, a loss of 12 099, the lake beds;
- the ~525 ON links outside the three populations (inside the based set, or across a border);
- the pooling of three populations with different intercepts.

The split between these three is not measured: it is named.

**Carved trunk cells more than 1 m off their law:**

| variant | cells | above the law | \|z − law\| p50 / p90 | inside the based set | distance to a based shore p50 / p90 | within 2 km |
|---|---|---|---|---|---|---|
| OFF | 1 354 of 48 128 | 81 | 3.2 / 14.9 m | 344 | 7.23 / 15.53 km (F132 set) | 31.2 % |
| ON F132 | 1 658 of 35 758 | 34 | 11.9 / **324.7 m** | 647 | 3.96 / 15.38 km | 43.8 % |
| ON extended | 1 816 of 38 555 | 35 | 10.0 / **313.7 m** | 697 | 1.42 / 11.62 km | 56.1 % |

- ON adds about 300–460 such cells, half of them **inside the based set**. They are the lakes' carved cells,
  where the law is now the lake level. Their p90 off-law depth rises from 15 m to about 320 m.
- **Upstream** the count barely moves (384 → 385 F132, 514 → 600 extended); **downstream** there are 0–1
  cells.

## U — the 17.5 points of trunk cells under their law

The totals are reproduced: 32.1 % (22 775) → **49.6 %** (35 146) F132, → 45.6 % (32 348) extended.

| | new cells | distance to a based shore p50 / p90 | within 2 km | law − terrain p50 / p90 | where |
|---|---|---|---|---|---|
| ON F132 | **12 495** (+17.6 pt) | **0.00 / 0.00 km** | **99.6 %** | **62.3 / 207.6 m** | inside 12 385, upstream 109, untouched 1 |
| ON extended | 12 486 (+17.6 pt) | 0.00 / 0.00 km | 99.6 % | 62.4 / 207.7 m | inside 12 376, upstream 110 |

- **99.1 % of the new under-law cells are INSIDE the based lakes.** They are the lake beds' trunk cells: ON
  sets their law to the lake level, which lies above the bed, so they are "under their law" and uncarved.
- They are **neither delta-fill places nor long reaches with a law set too high** (F130-E). They are the
  water body, where a floor law does not apply.
- Their run ends at a based shore 99.8 % of the time, after one cell (0.05 km).
- The under-law share of F133-T counts lake beds. **Its +17.5 points are an artefact of the instrument**,
  not a change of the land's profile.

## Δ — the residual of 12 509 cells (témoin, extended; non-blind)

The residual is reproduced: **12 509** cells, |Δz| p50 **6.50** / p90 **68.80** / max **631.2 m**.
- A has 623 polylines (190 283 samples); B has 633 (191 039 samples).

The sample `carve_diag` laid each residual cell from, in A and in B:

| laying sample | F133-F class | cells | p50 / p90 / max \|Δz\| |
|---|---|---|---|
| **MOVED** (another position) | no differing cell on its path | **11 048** | 7.28 / 107.17 / 631.2 m |
| MOVED | an origin on its path | 223 | 1.30 / 13.00 / 204.5 m |
| no sample on one side | no differing cell | 1 003 | 4.63 / 27.11 / 68.3 m |
| no sample on one side | an origin on its path | 213 | 2.68 / 9.12 / 20.9 m |
| the same sample | either | 22 | ≤ 8.8 m (the cross-line minimum with another polyline's cone) |

- **90.1 % are laid by a moved sample.**
- **The function** is `carve_diag` (called by `carve`): the nearest-source vector propagation
  (`valley_construction.rs:1489–1524`) and its cone `geo` (`:1468–1488`). Each cell is laid from the NEAREST
  polyline sample. The samples are built by `line_samples` (`:864`) in `skeleton`, so a different flat pointer
  moves a line and its samples.
- **It is the same code as Finding 127's "laid from the NEAREST sample"**: both blocks are identical, line for
  line, to commit `0ae6243` (Finding 127). `carve` has been a byte-identical wrapper of `carve_diag` since
  Finding 130-A0.

## Predictions

**Mine:**
- **P-0.3**: identity held; whether the protected breach was needed is **not judged** (both built together).
- **P-O**: "> 90 % constant" held. "The first difference at the bathymetry; light pass leaves the sea
  bit-identical" is **refuted**: 69 ocean cells differ at the light pass.
- **P-D**: held. Conditioned field, outside the footprint, z_ON > z_OFF, on OFF-carved cells, of the order of
  the col above OFF's floor, not of the lake's depth. One wording was wrong: "no water" (the conditioned field
  carries the pre-breach lakes' surfaces).
- **P-T**: **refuted**. Upstream ≤ 0.43 against 0.492 / 0.494; extended downstream 0.512 outside 0.47–0.50;
  "within 2 W of a shore" refuted. Untouched 0.49 ± 0.01 held.
- **P-U**: **refuted** (median distance 0, 99.6 % within 2 km). Only "law − terrain > 30 m" held.
- **P-Δ**: held. The numbers were non-blind.
- **Meta**: not held. It was stated about the three predictions that DISAGREED with the reviewer, not
  about all of mine: *"Disagreements: P-0.3 (the protected breach), P-D, P-U. At least two are wrong."*
  P-U is wrong, P-D held and P-0.3 is not judged, so one of the three is wrong. P-T and P-O's first half
  are also refuted, but they were not among the three.

**The reviewer's:**
- **O** held: the normalisation; 100 % within 10⁻³.
- **D** refuted: the 1 177 m is upstream, not in the footprint; upstream p90 1 182–1 240 m, not < 100 m.
- **T** refuted on "upstream < 0.45" (0.492 / 0.494). "Elsewhere 0.49 ± 0.02" held, except extended
  downstream at 0.512.
- **U**: "> 60 % within 2 km" held (99.6 %); "median < 20 m" refuted (62 m).
- **Δ**: "same code as F127" held; "p90 < 2 m" refuted (68.8 m).
- **Identity outside the crater** held.
- **Meta** held.

## Limitations, stated

1. T's populations walk the OFF skeleton. A link is placed by its upstream cell, and upstream wins over
   downstream. The fall of the pooled θ is located only by elimination.
2. D's "carved" is `carve` on S1 (the construction input with the rims applied). The pipeline carves before
   the rims.
3. O's stage (c) shows no droplet-erosion difference: the counts are identical to (b)'s. It is reported as
   measured, not explained.
4. Distances are D8 (Chebyshev) in cells, capped at 48.8 km.

## State

**Uncommitted**, awaiting the go-ahead.
- `bench_guard.rs`: the lake guard and its test.
- `volcanism/mod.rs`: the crater mask and pass, moved from the viz.
- `hd.rs`: calls them; `lake_guard`.
- `workspace.rs`: the lake badge, the guard asserts the lakes, `f134v_lake_listing`.
- `data/bench_lake_hashes.json`.
- `tests/common/mod.rs`: `build_world`, `viz_hd_lakes_on`, `bench_eroded_key`.
- `f123_teeth.rs`: `f134_lake_guard`.
- `f126_coast.rs`: `f134_*`.
- This folder, and ADR Finding 134.
