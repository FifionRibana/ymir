# Finding 122 — the "brush" lakes are constructed valleys drowned by inland below-sea basins; χ from the col removes them, and it is not a local gesture

Draft for the ADR, not yet appended. Seed 1, 8192², 400 km. No promotion:
`ValleyConstruction::basin_base` is `false` by default and Finding 121's path is unchanged (574
library tests pass, including `basin_base_keeps_the_floors_above_the_water_of_a_closed_basin`).

## Grep (rule 11/11b)

`base_level_floor_per_basin` **0** (not an identifier) · `spill_level` 6 (earliest L11548, Finding
85: *"`level = spill` for exorheic"*) · `z_base` 2 (L14355, Finding 96: *"χ = 0 at the sea"*) ·
`Finding 85` 23 · `Finding 39` 45 · `a_eq = ∞` 6 (L2778) · `Finding 121` **0** (not yet in the
ADR) · `drowned` 61 (earliest **L1306**: *"a coastal depression filling to 613 m, which DROWNED #28
(a mountain lake) inside its footprint"*) · `noyé` 10 · `lake_map` 52 · `valley floor` 5 (L28, the
acceptance criterion) · `per-basin` 7 (L1206, the per-basin spillway trace). ⇒ **A χ base per closed
depression was never considered**: Finding 96 integrates χ from the sea only. NOTHING FOUND.

## A — the brush, attributed

**Frames.** The viz frames seed 1 at torus offset x = 0.09375 (the author's export manifest), the
benches at 0.0: a viz cell `x_v` is bench cell `x_v + 768`. The viz shows the rows NORTH UP
(`workspace.rs`, Finding 27: the data is stored y = 0 = south), so "lower right" is large `x_v`,
small data `y`. **The named lake is 1000011** (floor cell 31634531, bench centre (5159, 3549), viz
(4391, 3549)): a family-2 lake, i.e. an inland below-sea basin, 947 km² at 459.3 m.

**The round's criterion does not hold as posed.** The arms (opening by the law's floor width) are
0.3–59 % constructed floor, not ≥ 80 %: what is drowned is the whole constructed VALLEY, walls
included. 65.2 % of 1000011's cells are carved cells.

**The round's control does not read 0.** The delivered world's lakes over the same masks read 1–45 %
floor and up to 100 % carved: lakes form in valleys, and the construction lays its valleys there.
The overlap instrument cannot separate a brush from an ordinary lake on its own.

**What separates them is the same basin across worlds** (same floor cell, or within 3 cells of it:
39767820/39767817 and 38218072/38218071; same level to the decimetre):

| floor cell | C2/10 km² / level / D_L | delivered km² / level / D_L |
| --- | --- | --- |
| 31634531 (**1000011**) | **947 / 459.3 m / 7.89** | 605 / 367.5 m / 7.76 |
| 39767820 | 668 / 115.3 / 8.09 | 757 / 115.3 / 8.58 |
| 38218072 | 656 / 49.3 / 3.38 | 867 / 49.3 / 7.52 |
| 48401756 | 506 / 614.9 / 3.54 | 560 / 614.9 / 5.49 |
| 28970982 | 483 / 76.3 / 3.91 | 661 / 76.3 / 6.73 |

1000011 is the one basin that GREW: +342 km², and its carved cells under water went from ~217 km²
(delivered, same masks) to ~617 km² (C2/10). **The extra area is drowned constructed valley.**

**The mechanism.** The skeleton takes `h ≤ sea` as base level (the drainage chain's convention), so
a trunk ending in an inland below-sea basin gets χ = 0 at sea + 0.5 m, and its valley is laid down
to that height. The basin's water, though, stands at its col (459.3 m) — and that col is never
lowered, because the basin's outflow is a `Spillway` traced over the col OUTSIDE the D8 network
(Finding 121-E): no constructed valley ever reaches it. Valleys dug hundreds of metres below a
water level that does not move are drowned for their length; with flat floors, planar walls and
semicircular cone heads, a drowned valley is a uniform-width stroke with a round end. That is the
brush (`brush_1000011_c2_10_valleys.png`).

## B — χ from the col (`basin_base`, gated)

A priority flood seeded on the OPEN OCEAN only (`water_class` 1) gives every closed depression its
spill level: inland below-sea basins (Finding 85's spill) and the lakes the pre-incision drainage
holds flat. A trunk entering a closed depression takes that spill as its base and counts χ from the
depression's col. Outside closed depressions nothing changes.

⛔ **Local by its rule, not by its effect.** It moves the floor on **50 087 of 70 208 trunk cells
(71 %)**, Δfloor p10 / p50 / p90 **+48.8 / +174.4 / +458.8 m**; carved cells 3.44 M → 2.69 M
(−22 %). On seed 1 most trunks end in an inland basin or cross a held lake.

| world | coast | Δ class e/b | lakes ≥ 1 km² (fam1) | D_L > 5 | lake % | relief p50 | R8 terrain | R8 network | σ p50 | transverse p50 | wall |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C2/10 (Finding 121) | +3 | 0 / 0 | 16 (4) | 19 % | 15.92 | 457.1 m | 0.0608 | 0.3783 | 6.34 m | +0.39 m | 28.0° |
| **C2/10 + basin base** | +5 | 0 / 0 | 18 (6) | **0 %** | 12.92 | **515.9 m** | **0.0461** | 0.4176 | **5.56 m** | +0.37 m | 27.7° |
| ON A1+B2 (Finding 121) | +1 | 0 / 3 | 29 (18) | 10 % | 19.32 | 508.3 m | 0.0484 | 0.3271 | 5.27 m | −0.01 m | 10.6° |

Skeleton control: PRE vs C2/10 + basin base, trunks ≥ 10 km², 90.5 % within 5 cells, p90 5 / 12,
mouths 152 of 184 — passes.

Per lake: every family-2 lake becomes a bowl (D_L 3.4–8.1 → 1.8–4.1, carved share ~0–1 %) at an
UNCHANGED level; **1000011: 947 → 494 km², carved 65.2 → 0.2 %, arms 3 → 1, D_L 7.89 → 2.29**
(`brush_1000011_c2_10_basin_valleys.png`).

**What it costs**: relief p50 +12.9 %, σ p50 −12 %, R8 terrain −24 %, R8 network +10 %, the
lake fraction −3.0 points. The valleys of every basin-draining trunk are shallower by construction.

## C — the tile, for the author's eye

`tile_2048_5120_{livre, a1b2, c1_bare, c2_10, c2_10_basin}.png` and
`brush_1000011_{…}.png` in this folder: hillshade from the NW, lakes blue, `Watercourse` points dark
blue, NORTH UP. The tile is bench (2048, 5120), 1024²; in the viz it sits at x ≈ 1280.

What the renders show that no column above measures (for the author to confirm or dismiss):

- straight parallel river runs (horizontal / vertical dark-blue strokes) filling flats, in every
  world including the delivered one;
- in the SE of the brush crop, polygonal cells where planar walls meet along straight crests, and
  semicircular valley heads at the trunk tips;
- a filled rectangular patch of river points near the SE corner of the brush crop.

The viz gains the state **"C2 /10 col"** (Finding 122-B) beside off / C1 nue / C2 /10 / C2 /3.
