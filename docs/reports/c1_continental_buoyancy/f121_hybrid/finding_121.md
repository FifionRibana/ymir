# Finding 121 — the hybrid at 8192²: the shape constructed along the pre-incision trunks, the texture left to a light incision

Draft for the ADR, not yet appended to it. Seed 1, 8192², 400 km, every column at the product's
resolution. No promotion: `FbmUpscaleConfig::valley_construction` is `None` by default and the
default path is byte-identical (the 573 library tests pass; the rule-13 test of the primitive is
`the_valley_construction_lowers_only_and_lays_its_floor`).

## The three reasons Option A was rejected (C1 design doc §3.1, verbatim)

1. *"cannot produce coherent chronological diversity (young vs old mountain ranges on the same
   continent)"*;
2. *"parameter calibration drifts indefinitely with no physical anchor"*;
3. *"'looks invented' is a real failure mode for users who can identify imposed geometry by
   eye"*.

What this round measured against each:

| reason | what answers it here | status |
| --- | --- | --- |
| (1) chronology | the age `k` is a monotone geometric parameter (C3: relief p50 439.5 → 457.1 → 480.7 m for k × 0.7 / 1 / 1.4) | **necessary, not sufficient**: one `k` per continent says nothing about young and old ranges side by side |
| (2) calibration | every law labelled ANCHORED or PROXY; the one PROXY on a target, `k`, is calibrated on Finding 95's oracle and lands on Finding 96's own value (0.07186 against 0.0719) | answered for `k`; the width coefficient and the skeleton smoothing remain PROXY without a target |
| (3) "looks invented" | the author's eye (block D, pending) — and one statistic that says the bare construction IS detectable: R8 terrain **0.0798**, six times the pre-incision field's 0.0131 | open until D |

## Grep (rule 11/11b)

`Option A` ADR **1** hit (L14013 — a different fork, not C1's) · C1 doc §3.1 (the source) ·
`template` 10 (earliest L15524, the round template) · `gabarit` 1 (L17105) · `invented` 2 (earliest
L274: *"Detail must be either INVENTED (FBM) or DERIVED (erosion physics)"*) · `inventé` **0** ·
`Finding 96 B2` **0** (the heading reads "B2 — the price", 2 hits, L9366) · `Finding 120` **0**
(not yet in the ADR) · `pre_incision network` / `pre-incision network` **0** · `Knobs::no_incision`
1 (L17116: *"Pre-incision = `Knobs::no_incision()` … asserts 20"* spurs) · `light incision` **0** ·
`k_time` 62 (earliest L3738, Finding 90's `SHIPPED_K_TIME = 9000`) · `Finding 88 D1` **0**, `88-D1`
7 (earliest L11615) · `Finding 59` 19 (earliest L6121, a retraction banner) · `Spillway trace` 5
(earliest **L1176**: *"Every enclosed below-sea basin is MARKED … + gets its spillway traced"*) ·
`Finding 94` 2 (L13741: the spillway *"appended with `upstream: vec![]`"*).

## B — the laws, labelled (as coded in `tectonics_c1/valley_construction.rs`)

| law | value | source | label |
| --- | --- | --- | --- |
| floor profile | `z = z_base + k·χ`, `χ = ∫ (A0/A)^0.5 dx`, A0 = 1 km², A clamped at 0.1 km² | Harel et al. 2016 eqs. (4)–(7), PDF p. 15 | ANCHORED |
| age `k` | **0.07186** (floor surface paired p50 = 488.1 m, Finding 95's oracle) | Finding 96 calibrated the same way: 0.0719 | **PROXY** on a target |
| floor width | `W = a·A^0.3`, A km², W m | Clubb et al. 2022 p. 437 (exponent 0.3 ± 0.06) | exponent ANCHORED |
| width coefficient | `a = 100.2` so that W(10 km²) = 200 m ⇒ law W p10/p50/p90 over the trunks **212 / 272 / 525 m** | Clubb Table 2 is unit-inconsistent (Finding 120) | **PROXY** |
| walls | 28°, planar | Whipple & Tucker 1999 Table 1 p. 17,663 | ANCHORED |
| angle of repose 33° | not used (28° sufficed) | code, no page | PROXY, unused |
| base | sea + 0.5 m | Finding 83 | ANCHORED |
| skeleton smoothing | moving average over ±250 m, endpoints pinned | against Finding 96's D8 stamping | **PROXY** |
| trunks | A ≥ 10 km² in-domain = 4 194 cells at 8192² | Finding 120 | DECISION |
| interfluves | `min(field, valley)`: the tectonic + FBM field wherever no wall reaches it | — | construction |
| lakes | the pre-incision hollows no trunk crosses | Finding 117 | decision |

**Two corrections the rule-13 test forced before any measurement.** (i) Taking the minimum of the
cones over every skeleton sample lets a downstream sample, whose floor is lower, flatten the floor
over half a valley width: the long profile is then not the χ law. The cross-section is laid from
the NEAREST sample; the minimum is kept only between samples of different polylines, for
continuity where two valleys meet. (ii) A wall stops where it reaches the terrain, so a valley
never cuts through a ridge into the next basin.

## C1 — the bare construction

**Skeleton control first, declared before the run** (C1 must reproduce PRE at least as well as ON
does, Finding 120: PRE→ON p90 10 cells, mouths 52.7 %): **PASSES**.

| PRE vs | trunks A ≥ 10 km² within 5 cells | p90 there / back | mouths matched |
| --- | --- | --- | --- |
| C1 bare | 90.6 % | 5 / 4 cells | 155 of 184 |
| C2 /10 | 91.4 % | 4 / 9 | 168 of 184 |
| C2 /3 | 91.3 % | 4 / 14 | 165 of 184 |

Cell-exact overlap is only 14.6–18.4 %: the 250 m smoothing moves the valley centreline one or two
cells off the D8 line, and the recomputed D8 follows the valley.

## The columns

| world | coast | Δ class e/b | lakes ≥ 1 km² (fam1) | D_L > 5 | lake % | relief p50 |
| --- | --- | --- | --- | --- | --- | --- |
| PRE | +5 | 0 / 0 | 29 (18) | 0 % | 16.54 | 679.1 m |
| OFF delivered | +115 | 12 / 4 | 54 (43) | 35 % | 24.70 | 424.2 m |
| ON A1+B2 (control) | +1 | 0 / 3 | 29 (18) | 10 % | 19.32 | 508.3 m |
| **C1 bare** | +2 | 0 / 0 | 16 (5) | 19 % | 14.72 | 519.0 m |
| **C2 /10** | +3 | 0 / 0 | 16 (4) | 19 % | 15.92 | 457.1 m |
| **C2 /3** | +3 | 0 / 2 | 16 (4) | 25 % | 16.99 | 405.4 m |
| C3 k×1.4 (/10) | +1 | 0 / 0 | 17 (5) | 12 % | 15.40 | 480.7 m |

| world | R8 terrain | R8 network (≥ 10 km²) | comb tile terrain / network | σ p50 / p90 | in-thalweg (≥ 10 km²) | transverse p50 | interfluve Δ p10 / p50 | wall slope p50 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PRE | 0.0131 | 0.5990 (0.6064) | 0.0400 / 0.6344 | 3.70 / 11.00 m | 7.0 % (22.0 %) | +0.75 m | 0 / 0 | 5.8° |
| OFF | 0.0924 | 0.2056 (0.2014) | 0.1190 / 0.2576 | 6.96 / 36.29 m | 69.3 % (93.2 %) | −5.06 m | −401.2 / −30.5 m | 14.9° |
| ON | 0.0484 | 0.3271 (0.2496) | 0.0599 / 0.3525 | 5.27 / 30.27 m | 50.4 % (90.3 %) | −0.01 m | −277.0 / −0.2 m | 10.6° |
| **C1** | 0.0798 | 0.4875 (0.5691) | 0.0698 / 0.6007 | 5.36 / 21.20 m | 23.9 % (91.4 %) | +0.55 m | 0 / 0 (tautology) | 28.0° |
| **C2 /10** | 0.0608 | 0.3783 (0.5316) | 0.0486 / 0.4612 | 6.34 / 22.60 m | 26.2 % (64.0 %) | +0.39 m | −74.1 / −0.0 m | 28.0° |
| **C2 /3** | 0.0618 | 0.3477 (0.5261) | 0.0331 / 0.3576 | 6.73 / 26.73 m | 32.2 % (66.1 %) | +0.24 m | −209.4 / −0.1 m | 28.0° |
| C3 k×1.4 | 0.0544 | 0.3905 (0.5345) | 0.0422 / 0.4578 | 6.10 / 22.52 m | 26.0 % (67.8 %) | +0.42 m | −71.6 / −0.0 m | 28.0° |

`to_nothing` 0 in every world; unresolved 0–3. Costs (build, one 8192² world): OFF 51 s, ON 50 s,
C1 78 s, C2 105–110 s; the skeleton alone 64.8 s, the carve 1.7 s.

**Readings.**
- **C1 (Finding 96's question): neither too smooth nor combed, but STRIPED at the scale of the
  valleys.** σ p50 5.36 m is ON's (5.27) and −23 % against the delivered; R8 terrain 0.0798 is six
  times PRE's: the walls are planes whose orientation follows trunks that are still D8 lines at
  scales above the 250 m smoothing.
- **C2 lowers the anisotropy instead of raising it**: R8 terrain 0.0798 → 0.0608, network 0.4875 →
  0.3783 (/10) → 0.3477 (/3), the comb tile to 0.0331 at /3, the lowest of every world. σ p50
  climbs to −9 % (/10) and −3 % (/3) of the delivered.
- **Canyons**: the eroded-stage class stays 0 in every constructed world. The breached stage reads
  2 at /3, against ON's own 3 there; whether they are the same bodies is not checked.
- **Entrenchment does NOT appear.** The transverse offset stays positive (+0.24 to +0.55 m: the
  river cell sits above its lower bank) against −5.06 m delivered. C1's 91.4 % "in-thalweg" on the
  ≥ 10 km² reaches is an artefact of a flat floor, where equality with both banks passes the test;
  the light incision breaks the flat and the share falls to 64–66 %.
- **Lakes**: 29 → 16; 13–14 of PRE's 18 family-1 lakes are gone. A trunk crossing them and cutting
  their sill is the expected cause; it is not checked lake by lake.

## C3 — three ages (at k_time/10)

| k | relief p50 | law W p50 | lakes ≥ 1 km² | σ p50 |
| --- | --- | --- | --- | --- |
| × 0.7 | **439.5 m** (Finding 95 criteria, printed before the panic below) | 272 m | — | — |
| × 1 | 457.1 m | 272 m | 16 | 6.34 m |
| × 1.4 | 480.7 m | 272 m | 17 | 6.10 m |

**Monotone.** W does not move: the law has no `k` in it.

⛔ **The k × 0.7 world trips a production invariant in the drainage assembly**, on both runs, at the
same cell: *"ADR Finding 38/92-B: 1 enclosed below-sea component(s) carry no water body. (floor
cell, cells): [(28970982, 120300)]. This is the invariant Finding 38 closed and Finding 86's
`MergedUnionRelevel` re-opened"* (`hd_assembly.rs:279`). The construction cannot create a cell below
the sea (floors and walls stay ≥ 0.5 m above it); it CAN lower the col of an enclosed below-sea
basin. Measured, reproduced, **not attributed**.

## E — the F87-A object

It is basin **1000015**'s spillway: 5 256.41 m³/s, 71 cells, `segment_drainage_km2` 552 932 km²
signified = **9 829.9 km² in-domain**, against a D8 area of **1.3 km²** along its path. The code
says why (`hd_assembly.rs`, the spillway append): *"A spillway's path is traced over a col, OUTSIDE
the accumulation network, so the raster reads ~0 along it (measured). Its contributing area is the
basin's, which `drainage_km2` already carries."* **The discharge is the basin's balance and is
right; the path is the traced col; the D8 raster is not the network that carries it.** It becomes a
control again when read with `segment_drainage_km2`, not the D8 raster.

## D — À VALIDER VISUELLEMENT (the author)

The viz gains "Vallées construites (F121)": `off` / `C1 nue` / `C2 /10` / `C2 /3`, and an age
selector restricted to the three measured values. The four states, seed 1, **8192²**, 400 km:

| state | "Âge du continent" | "Vallées construites" |
| --- | --- | --- |
| livré | off | off |
| A1+B2 | on (0.024) | off |
| C1 nue | any | C1 nue |
| C2 hybride | forced on | C2 /10 |

Same continent crop as Finding 88, same tile (2048, 5120) in 8192² pixels. ⚠️ The viz sets the
same config field the bench sets and runs the same function; the viz run itself was NOT guarded
bit-for-bit against the bench (Finding 109's guard 2 was in-core).
