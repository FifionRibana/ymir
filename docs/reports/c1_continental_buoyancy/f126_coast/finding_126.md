# Finding 126 — the −1.00 m is the bathymetry's clamp, and the construction drowns the coastal walls; mur ↔ mer is proved by its clause; the teeth sit in the concave foot; the sub-valleys' stripes are the D8 lattice's, and no smoothing removes them

**One clause and one skeleton variant were built this round, both gated and off. Nothing is promoted.**
Bench `crates/ymir-core/tests/f126_coast.rs`; its raw outputs are the `.txt` files beside this report.

**Units.** The domain is 400 km on 8192² cells: one cell = 48.8 m (domain). "Signified" = ×7.5 for
lengths, ×56.25 for areas. Every length and area says which one it is. Heights are metres, not scaled.
Drained-area bands are domain km².

**Reading declaration, unfavourable.** My predictions (`f126_predictions.md`, scratchpad) were written
before any grep or measurement. They are non-blind on Findings 73–125 and on the reviewer's
predictions. I declared one fact I already held: the −1.00 m had been read on the field
`build_field_seed` returns.

## 0 — Finding 125 committed; the guard; this round's code is bit-neutral when off

- Finding 125 was committed as `5fb84c9`, one commit. The guard through `run_hd` then read **6 / 6
  "= banc"** (`guard_six_states.txt`).
- **This round's code, off, leaves the definition's world bit-identical.** A fresh, uncached build of
  "C2 /10 col (défaut)" hashes to `a8d2d538d692c2f0`, the reference (`identity.txt`, test
  `f126_identity`). The round added two gated fields and moved the skeleton smoothing into
  `smooth_positions`.
- The construction key's permanent test now also asserts the two new fields absent from the
  serialised form.
- `cargo test -p ymir-core --release --lib`: 582 passed. `cargo check --workspace --tests --release`:
  clean.

## A — the −1.00 m is `shelf_min_depth_m`, and the CONSTRUCTION puts the wall below the sea

**The line.** `crates/ymir-core/src/terrain/bathymetry.rs:190`:
`let min_depth = p.shelf_min_depth_m.max(1.0);`, then `:199`:
`(env * (1.0 + p.texture * dev)).clamp(min_depth, floor_depth)`. `apply_bathymetry_profile` takes every
cell `≤ sea` as ocean (`:173`), whether it is edge-connected or enclosed. For a cell just below the sea,
`cur_depth ≪ mean_depth`, so `1 + dev ≈ 0` and the clamp returns `min_depth` = 1.0 m. It runs LAST in
the upscale (`production_upscale.rs:635-649`), inside `build_field_seed`, so before the bench's breach
and before the HD assembly that asserts Finding 38.

**The dossier already had it.** Finding 80 B2 (ADR L9379): *"the DELIVERED one, whose bathymetry clamp
puts every drowned cell at exactly −1.0 m"*. Finding 125 called the −1.00 m unexplained after grepping
the code with an ASCII minus. Rule 11b (grep the VALUE in the dossier) would have closed it: "−1.0 m"
with U+2212 returns 2 hits, the first at L9379. That was my failure in Finding 125, not a new constant.

**The two cells, stage by stage.** B2 (A_c) + B1 was walked through the pipeline's own stage switches.
The breach is the bench's, as in Finding 125.

| stage | (4257, 2006) | (4613, 2860) | the world's land-in-S1 cells now ≤ sea |
|---|---|---|---|
| S1 input (FBM + C-2 craters) | 98.72 m · land | 182.61 m · land | — |
| **S2 + construction** | **−1.55 m · enclosed** | **−10.52 m · enclosed** | 6 331 (15.1 km², domain): 4 216 edge-connected OCEAN, 2 115 enclosed in 300 components (167 of ≤ 4 cells) |
| S3 + light pass (k/10, A1+B2) | −1.55 m · enclosed | −6.18 m · enclosed | 7 719 (18.4 km²): 5 950 ocean, 1 769 enclosed in 306 components |
| S4 + droplet erosion | −1.55 m | −6.18 m | unchanged |
| **S5 + bathymetry (the world)** | **−1.00 m** | **−1.00 m** | the same 7 719, **all 7 719 at exactly −1.00 m** |
| S6 + breach (the HD assembly's input) | −1.00 m | −1.00 m | 8 165: 6 075 ocean, 2 090 enclosed in 307 components |

**Where they are.** Both cells lie 2–3 cells (0.10–0.15 km, domain) from the open ocean, at S1 and at
S5. **They are coastal.** Finding 125's D left this "NOT measured"; it is measured now, and the
reviewer's "sous-vallées côtières" holds. No lake claims either cell (lake id 0 at S5).

**Why a wall passes below the sea, at the two cells** (switching the terms on the construction's
input; the carve's own output, not a recomposition):

| | (4257, 2006) | (4613, 2860) |
|---|---|---|
| planar wall (B2 → A_c) | 4.71 m | 42.19 m |
| + profile (foot + crest) | 2.60 m (−2.11) | 21.90 m (−20.29) |
| + noise (the detail, −4.16 / −32.43 m) | 0.55 m (−4.16) | 9.76 m (−32.43) |
| **both (B2 (A_c) + B1)** | **−1.55 m** (−6.27 = −2.11 − 4.16) | **−10.52 m** (−52.72 = −20.29 − 32.43) |

The two terms **add to the centimetre**, on walls whose planar value is only 4.7 m and 42.2 m above the
sea. Neither alone crosses at either cell. "The product profile × noise" (the reviewer) is a sum here.

**Where the constant touches without firing.** At the construction stage (direct carves on S1, same
skeletons):

| construction | land → ≤ sea | on walls / floors | edge-connected OCEAN | enclosed (components) |
|---|---|---|---|---|
| **témoin C2/10 col (planar)** | **0** | — | 0 | 0 |
| planar on A_c (B2 alone) | 0 | — | 0 | 0 |
| B1a profile only | 1 385 (3.30 km²) | 1 385 / 0 | 499 | 886 (69) |
| B1b noise only | 365 (0.87 km²) | 365 / 0 | 147 | 218 (72) |
| B1 both | 3 205 (7.64 km²) | 3 205 / 0 | 1 884 | 1 321 (94) |
| profile on A_c | 3 247 (7.74 km²) | 3 247 / 0 | 1 941 | 1 306 (216) |
| noise on A_c | 761 (1.81 km²) | 761 / 0 | 429 | 332 (152) |
| B2 (A_c) + B1 | 6 331 (15.1 km²) | 6 331 / 0 | 4 216 | 2 115 (300) |

- **A planar wall never drowns land** (0 in the témoin, 0 at A_c): `V = zf + tan·u ≥ zf ≥ sea + 0.5 m`.
- **Only WALL cells drown, never a floor**, in every world.
- The halves are superadditive world-wide (3 205 > 1 385 + 365; 6 331 > 3 247 + 761). At a threshold
  (the sea) a sum crosses where neither term alone does.
- **Finding 38 fired in one world out of seven; the construction drowns walls in every B1 world.** B1
  alone lays 94 enclosed components below the sea at the construction stage. Their fate after S2 is
  not measured. F38 did not fire in B1 (Finding 125), so whatever survives to the HD assembly is
  covered by a water body: a flat one at −1 m, Findings 33/34/38's old "−20 m" symptom at its new
  depth. The constant does not only fire; it floors every drowned wall of every B1 world.
- **The edge-connected drowned walls are the coast retreating into the valleys**: 499 (profile), 147
  (noise), 1 884 (both), 4 216 (with A_c). Finding 125's new spurs follow the same order: 42 / 8 / 81 /
  195. B shows they are the same effect: the clause that forbids them removes every new spur.

**Instrument caveat, stated.** The bench's direct carve on S1 is NOT bit-identical to the pipeline's
S2: 203 352 cells differ (0.30 %). S1 is read with the C-2 rim reconstruction already applied (it runs
at the end of every build), while the pipeline constructs before it. That is a hypothesis for the
difference, not verified. The two cells (−1.55 / −10.52 m) and the census (6 331 / 4 216 / 2 115 / 300)
are identical to the pipeline's, so these rows are read. My recomposition from the Euclidean-nearest
sample does NOT reproduce the carve at either cell (11.63 against 4.71 m; 42.64 against 42.19 m). The
carve's value there comes from ANOTHER polyline's cone: either the propagated `who` differs from the
Euclidean nearest, or the continuity minimum across different lines. The two routes are not separated.
**Both cells are junction cells.** Only the switched carves above are read.

**Which term drowns a wall (A2, the author's question "profil concave au pied sous le fond ? bruit
reprojeté négatif ? somme des deux près de la côte ?").** The terms are switched one at a time on S1,
with the same skeletons. "Terrain" is the input terrain at the drowned cells; "to the ocean" is the
distance to the open ocean in S1 (km, domain).

| term(s) | C2/10 skeleton: drowned (ocean / enclosed) | terrain p50 / p90 | to the ocean p50 / p90 | A_c skeleton: drowned (ocean / enclosed) | terrain p50 / p90 | to the ocean p50 / p90 |
|---|---|---|---|---|---|---|
| **planar (témoin)** | **0** | — | — | **0** | — | — |
| foot only | **0** | — | — | **0** | — | — |
| crest only | 173 (18 / 155) | 13.5 / 19.8 m | 0.29 / 0.78 | 385 (171 / 214) | 13.6 / 20.0 m | 0.10 / 0.54 |
| noise only | 365 (147 / 218) | 198.7 / 637.5 m | 0.34 / 1.07 | 761 (429 / 332) | 71.6 / 345.8 m | 0.15 / 0.59 |
| foot + crest (B1a) | 1 385 (499 / 886) | 13.1 / 24.8 m | 0.24 / 0.78 | 3 247 (1 941 / 1 306) | 13.6 / 21.8 m | 0.10 / 0.49 |
| all three (B1) | 3 205 (1 884 / 1 321) | 39.8 / 500.8 m | 0.34 / 1.12 | 6 331 (4 216 / 2 115) | 28.1 / 160.7 m | 0.15 / 0.59 |

- **Not "the foot below the floor".** The foot alone drowns nothing, and it cannot: its rise
  `tan·u²/(2·foot)` is ≥ 0 above a floor laid at ≥ sea + 0.5 m.
- **The crest drowns where the terrain is LOW.** Its smooth minimum takes up to `crest_m·ln 2` = 13.9 m
  off a wall that meets a terrain lying 13–20 m above the sea. Those are low coastal rims.
- **The noise drowns HIGH walls** (p50 199 m on the C2/10 skeleton), wherever its signed detail term
  is strongly negative near the coast.
- **The foot is a MULTIPLIER, not a cause**: foot + crest drowns 8× the crest alone (1 385 against 173).
  My reading, not measured: the lower wall reaches the terrain farther out, so the carve extends into
  more low coastal terrain for the crest to dip.
- **Every drowned wall is coastal**: p90 ≤ 1.12 km (domain) = 8.4 km (signified) from the open ocean.
  Near the coast the floor is laid at sea + 0.5 m (Finding 83's base), so any lowering of the wall
  crosses the sea. Inland the floors are higher (k·χ) and the same lowering does not reach it.
- **At the two cells, it is the foot + the noise**: foot only 2.76 / 21.91 m (−1.95 / −20.28 against
  planar), crest only 4.53 / 42.17 m (−0.18 / −0.02), noise −4.16 / −32.43 m.

## B — mur ↔ mer: the base rate first, then the clause. The clause brings the spurs back to the base rate and Finding 38 to zero, and touches nothing else

**The base rate.** Spurs are counted at the eroded stage with Finding 75's detector (the u16 mask).
Each is placed NEAR a coastal wall (≤ 2 km domain = 15 km signified) or ELSEWHERE, over the coastline
length of each zone (km, domain). The contour's own segments make the length. Each spur is typed
INLET or PENINSULA by the cell halfway from its neck to its farthest sample. Masks and "coastal wall"
are Finding 125's (the world's skeleton on PRE; carved, not floor, land, ≤ 3 cells from the sea).

| world | spurs (inlets / peninsulas) | NEAR a coastal wall: n / km → rate | × témoin | ELSEWHERE: n / km → rate | new vs témoin (inlets / peninsulas) |
|---|---|---|---|---|---|
| **témoin C2/10 col** | 20 (10 / 10) | 3 / 210.8 → **0.0142 /km** | **1** | 17 / 1 426.4 → 0.0119 | — |
| B2 (A_c) | 19 (9 / 10) | 5 / 747.3 → 0.0067 | 0.47 | 14 / 888.3 → 0.0158 | 0 |
| B1a profile | 66 (27 / 39) | 50 / 425.9 → **0.1174** | **8.3** | 16 / 1 308.3 → 0.0122 | 42 (13 / 29) |
| B1b noise | 27 (11 / 16) | 10 / 248.4 → 0.0403 | 2.8 | 17 / 1 424.1 → 0.0119 | 8 (2 / 6) |
| B1 both | 105 (52 / 53) | 89 / 503.2 → **0.1769** | **12.5** | 16 / 1 308.5 → 0.0122 | 81 (41 / 40) |
| B1a + B2 (A_c) | 107 (46 / 61) | 94 / 1 116.9 → 0.0842 | 5.9 | 13 / 723.1 → 0.0180 | 83 (33 / 50) |
| B2 (A_c) + B1 | 223 (112 / 111) | 210 / 1 271.8 → 0.1651 | 11.6 | 13 / 722.5 → 0.0180 | 195 (99 / 96) |
| **B1a + CLAUSE** | 20 (10 / 10) | 4 / 316.5 → **0.0126** | **0.89** | 16 / 1 321.4 → 0.0121 | **0** |
| **B1 + CLAUSE** | 19 (9 / 10) | 3 / 315.0 → **0.0095** | **0.67** | 16 / 1 321.3 → 0.0121 | **0** |
| **B2 (A_c) + B1 + CLAUSE** | 19 (9 / 10) | 6 / 907.6 → **0.0066** | 0.47 (B2's own: 0.99) | 13 / 731.1 → 0.0178 | **0** |

- **The témoin's walls make no spurs**: 0.0142 /km near a wall against 0.0119 elsewhere, ×1.2. **B2
  alone does not either**: it multiplies the coastline next to a wall by 3.5 (211 → 747 km), and the
  rate there FALLS to 0.0067. A wall, coastal or not, does not make a spur.
- **The wall profile does**, at the wall × sea contact: ×8.3 (B1a), ×2.8 (B1b), ×12.5 (both) near a
  coastal wall. **Elsewhere the rate does not move** (0.0119–0.0122 in every 10 km² world), the
  within-world control. The 100 % of Finding 125 now has its denominator.
- The new spurs are inlets AND peninsulas, about half each (B1: 41 / 40). Inlets are drowned valleys
  (Part A: 1 884 wall cells laid below an edge-connected sea by B1's construction). A peninsula would
  be the land left between them. That reading is not measured.

**The clause, gated: `ValleyConstruction::wall_sea_floor_m = Some(ε)`, ε = `base_m` = 0.5 m.** No WALL
cell is laid below sea + ε; the floors' own base (Finding 83), transposed to the walls. It is one
`max` in `carve`, after the crest and the detail, on wall cells only. The clause's other half ("nothing
raises a coastal cell above the original terrain") needs no code: `carve` only lowers,
`out = min(field, V)`, asserted by the tests. Permanent test (rule 13):
`the_coast_clause_keeps_every_wall_above_the_sea`. Its negative control comes first (without the
clause a wall drowns). It also asserts that every cell the unclamped carve left above sea + ε is
bit-identical.

| | B1a (without → with) | B1 (without → with) | B2 (A_c) + B1 (without → with) |
|---|---|---|---|
| cells the clause moves (construction stage, on PRE) | 1 503: 1 397 would drown (503 ocean-connected, 894 enclosed), 106 would stay just above sea | 4 348: 4 211 would drown (2 873 / 1 338), 137 | 24 422: 24 082 would drown (21 575 / 2 507), 340 |
| spurs near a coastal wall, rate | 0.1174 → **0.0126 /km** (témoin 0.0142) | 0.1769 → **0.0095** | 0.1651 → **0.0066** (B2 alone 0.0067) |
| new spurs vs the témoin | 42 → **0** | 81 → **0** | 195 → **0** |
| Finding 38 | holds → holds | holds → holds | **⛔ fires → holds** |
| teeth | 2 130 → 2 145 (**+0.7 %**) | 4 719 → 4 805 (**+1.8 %**) | — (F38) → 1 453 |
| over-dug (canyons) | 0 → 0 | 0 → 0 | 2 → **2** (the same bodies 14 and 15, same floors) |
| R8 terrain | 0.0497 → 0.0498 | 0.0490 → 0.0492 | 0.1032 → 0.1034 |
| R8 network | 0.5513 → 0.5521 | 0.5408 → 0.5407 | — → 0.5054 |
| teeth positions (C) | identical to the percent | identical to the percent | — |

**⇒ "mur ↔ mer" is PROVED at this seed and resolution.** Its clause returns the spurs to the base rate
in all three worlds, B2's own included. It removes every new spur and closes Finding 38. It changes
nothing it should not: teeth +0.7 % and +1.8 %, their positions unchanged, the canyons, both R8s. The
canyons are unchanged because they are FLOORS, which the clause does not touch.

## C — where the teeth sit on their wall: at the concave FOOT under B1, and only under B1

**The instrument.** Every point of every tooth (Finding 123-A's teeth, Finding 125's masks) is placed
on its wall by the planar geometry of the nearest skeleton sample (Euclidean, bucketed): FLOOR
(d ≤ W/2), FOOT (u < 300 m beyond the floor edge, the profile's foot length), CREST (the planar wall
within 60 m = 3·crest_m of the terrain) or STRAIGHT between. The SAME bins serve every world, so B1b's
planar wall is read on B1's foot and crest. The wall's own area (every 13th wall cell) is the base
rate for "uniform".

| world | teeth | mouths: floor / FOOT / straight / crest | heads: floor / FOOT / straight / crest | points: FOOT (× its area share) / straight | the wall's area: foot / straight / crest |
|---|---|---|---|---|---|
| **témoin (planar)** | 1 440 | 67.3 / **28.3** / 4.3 / 0.1 | 0.6 / **8.5** / 83.6 / 7.4 | **21.4 % (×0.96)** / 72.2 | 22.4 / 70.4 / 7.1 |
| B1a profile | 2 130 | 49.7 / **46.4** / 3.8 / 0.1 | 1.1 / **22.2** / 63.6 / 13.1 | **27.8 % (×1.29)** / 63.5 | 21.6 / 62.8 / 15.4 |
| B1b noise | 1 681 | 47.0 / **49.8** / 3.2 / 0.1 | 1.4 / **18.0** / 71.9 / 8.7 | **24.6 % (×1.10)** / 68.8 | 22.4 / 70.5 / 6.9 |
| **B1 both** | 4 719 | 9.2 / **89.6** / 1.2 / 0.0 | 1.4 / **59.3** / 30.1 / 9.2 | **48.9 % (×2.26)** / 44.2 | 21.6 / 62.8 / 15.4 |

- **Under B1 the teeth live IN THE FOOT.** 89.6 % of their mouths, 59.3 % of their heads and 48.9 % of
  their points sit in the concave foot, 2.26 × its share of the wall. They are born there and end
  there, short of the floor (mouths on the floor fall from 67 % to 9 %). **Neither half does this
  alone** (the foot's point share ×1.29 under the profile, ×1.10 under the noise). This is the
  superadditivity of Finding 125's teeth, located.
- **The mechanism the round named is supported where it predicts.** The profile's foot flattens the
  lower wall (its slope ramps from 0 at the floor edge), and the noise puts micro-relief on that flat,
  which the D8 routes ALONG the foot instead of down it. Under B1b the planar wall stays steep and the
  points spread with the wall's area (×1.10).
- **The mouth alone does not discriminate.** In every world ≥ 95 % of mouths sit at the wall's bottom
  (floor + foot), because a tooth ends where the wall meets the floor. The discriminating reads are the
  HEAD and the POINTS.
- The clause changes none of this (B1 + CLAUSE: mouths at the foot 89.7 %, points 49.3 %). The teeth are
  an inland foot defect, and mur ↔ mer is not their clause.

## D — the bent skeleton: the 1–3 km² branches sit AT the D8 floor, and no moving average lifts them off it

**This is not Finding 113's snap.** Finding 113 moved the EXPORTED river polyline toward a thalweg it
already occupied (MFD ≡ D8 at 93 %, Finding 15), and smoothed that drawn line (comb R8 0.466 → 0.452 at
±8 cells). The terrain never changed there, so nothing the terrain's R8 reads could move. Here the
SKELETON's geometry is changed BEFORE the construction digs, so the walls move with it. The two share
one result, measured independently: a moving average does not remove a D8 line's axis alignment.

**The chord instrument.** R8 of the skeleton's chords by drained-area band, the band read at the
chord's middle sample. Two chords are used:
- 8 densified samples (≈ 190 m), Finding 125's chord. It reproduces its population exactly: 50 683
  chords at 1–3 km² on PRE.
- 200 m of arclength: the law in metres.

### D1 — the D8 floor, built before it is used

The synthetic field is a provably isotropic relief (256 plane waves of uniform random direction,
Finding 97's family) on a 1 500 m dome draining to a surrounding sea. It is 2048² at the SAME cell,
48.8 m (a 100 km domain). Its terrain R8 (w 16) is 0.007–0.025, at or below Finding 97's isotropic floor
of 0.0402. The D8 tracé answers to the terrain's ROUGHNESS, so four roughness levels were built.

⚠ The first roughness instrument let the sea into its ±1 km blur and read 0.455 on a field built at
ρ = 0.1. It is not read. The corrected one keeps the window entirely on land (every cell ≥ 21 cells
from the sea).

| field | roughness ratio (±1 km trend) | raw D8, ≥ 10 | 3–10 | **1–3 km²** | with the 250 m smoothing, 1–3 |
|---|---|---|---|---|---|
| synthetic ρ = 0.1 | 0.071 | 0.769 | 0.879 | **0.986** | 0.985 |
| synthetic ρ = 0.3 | 0.209 | 0.734 | 0.749 | **0.796** | 0.778 |
| synthetic ρ = 1 | 0.613 | 0.272 | 0.320 | **0.362** | 0.306 |
| synthetic ρ = 3 | 1.071 | 0.240 | 0.186 | **0.240** | 0.223 |
| **PRE / S1 (the real terrain)** | **0.379** | **0.394** | **0.438** | **0.601** | **0.597** |

- **There is no single D8 floor.** On a smooth field D8 draws straight lines along its eight axes
  (0.986); on a rough one the path wanders (0.240). The floor is a function of the roughness.
- **At the real terrain's roughness (0.379, identical on PRE and S1), the floor at 1–3 km² interpolates
  to 0.56 (log) – 0.61 (linear). The real branches read 0.601. They are AT the floor.** F125's "0.62"
  is the D8 lattice on a terrain this smooth, not a tracé error above it. The reviewer's "~0.45, above →
  reducible" is refuted at this roughness.
- **The real trunks are BELOW their floor** (0.394 against 0.48–0.54 interpolated): real topography
  deflects the trunks and does not deflect the small branches. The small branches run down hillslopes
  that are planar at their scale.
- The interpolation spans a steep interval (0.80 → 0.36 between two synthetic levels), and the synthetic
  spectrum is not PRE's. "At the floor" carries that width.

### D2 — the smoothing length L (the HALF-length, like `smooth_m`)

**Does the existing 250 m smoothing apply to the sub-valleys? YES, to 82 %.** It leaves raw 17.7 % of
the 1–3 km² chords (short lines, and the 5 pinned samples at each end), 6.6 % at 3–10 and 5.7 % at
≥ 10. Those raw chords are LESS aligned (0.458) than the ones it smooths (0.631, read raw). The pinning
is not where the alignment is.

On PRE; S1 reproduces every value within 0.001:

| L | ≥ 10 (8-pt / 200 m) | 3–10 | **1–3 km²** (8-pt / 200 m) | trunks ≥ 10 moved p50 / p90 / max |
|---|---|---|---|---|
| raw D8 | 0.394 / 0.379 | 0.438 / 0.422 | **0.601 / 0.594** | 0 / 0 / 0 m |
| **250 m (the existing, témoin)** | 0.382 / 0.446 | 0.434 / 0.494 | **0.597 / 0.640** | 4.4 / 36.4 / 133 m |
| 0.5 W | 0.359 / 0.431 | 0.392 / 0.460 | **0.572 / 0.625** | 0.0 / 27.9 / 136 m |
| 1 W | 0.395 / 0.458 | 0.409 / 0.480 | **0.568 / 0.634** | 8.6 / 46.5 / 296 m |
| 2 W | 0.450 / 0.507 | 0.452 / 0.514 | **0.584 / 0.651** | 15.6 / 75.5 / 585 m |

- **No L lowers the 1–3 km² alignment by more than 5 %** (best: 1 W, 0.568 on the 8-point chord). On
  the 200 m chord, EVERY smoothing RAISES it. A moving average turns a D8 staircase into its mean
  direction, and on a terrain this smooth the mean direction IS an axis. What is quantised is the
  direction of a long straight run, not a zigzag.
- The trunks become MORE aligned as L grows (0.394 → 0.450 at 2 W) and move farther (max 585 m at 2 W,
  within the 2 W ≈ 1 km of a large trunk). The control ("the trunks move no more than the smoothing
  requires") holds: every maximum lies under its own L.
- **Consequence, named and not built.** The alignment is the D8 LATTICE's at this roughness, and a
  moving average cannot take it off the lattice. What could is a tracé that is not on the lattice: a
  steepest-descent path integrated in continuous space along the terrain's gradient, rather than a chain
  of D8 cells. MFD would not do it (≡ D8 at 93 %, Finding 15).

### D3 — B2 built on the skeleton smoothed at the best L of D2 (L = 1 W): nothing moves

`ValleyConstruction::smooth_w = Some(1.0)`, gated. Permanent test
`the_width_law_smooths_a_line_the_fixed_window_leaves_raw`, whose negative control is the fixed window
leaving a 9-sample line raw. 1 W is the best L of D2 on the 8-point chord (−5 %); none was better. The
témoin and B2 at 250 m were rebuilt in the same run, and both reproduce Finding 124-P4 to the digit.

| | **témoin C2/10 col** | B2 ≥ 1 km², L = 250 m | B2 ≥ 1 km², L = 1 W |
|---|---|---|---|
| skeleton R8, 1–3 km² (8-pt) | — | 0.597 | 0.568 |
| **R8 terrain** (× témoin) | **0.0452** | 0.1049 (**× 2.32**) | 0.1059 (**× 2.34**) |
| contribution, sub 1–3 km² (témoin on the same cells) | — | +0.0508 (+0.0128) | +0.0510 (+0.0129) |
| contribution, sub 3–10 km² | — | +0.0249 (+0.0082) | +0.0252 (+0.0082) |
| class R8, trunk ≥ 10 walls (témoin) | 0.1695 | 0.2161 (0.2351) | 0.2056 (0.2353) |
| class R8, interfluve (témoin) | 0.0126 | 0.0267 (0.0162) | 0.0294 (0.0162) |
| teeth | 1 440 | 262 | 278 |
| R8 network (chord 8) | 0.5613 | 0.5060 | 0.5056 |
| over-dug (canyons) | 0 | 0 | 0 |
| coast (spurs against PRE) | +3 | +4 | +4 |
| skeleton at 1 km² within 5 cells · mouths matched | 74.2 % · 509 / 877 | 91.0 % · 652 / 877 | 90.8 % · 653 / 877 |
| build | 156 s | 149 s | 149 s |

**The terrain's R8 × 2.3 is unchanged (× 2.34 against × 2.32).** The 1–3 km² contribution moves by
+0.4 %, the teeth by +6 %, the skeleton pairing by −0.2 points. The trunks' class R8 falls 5 %, within
the smoothing's own displacement (the control). **"hiérarchie ↔ tracé" in its moving-average form is
REFUTED, twice: on the skeleton (D2) and on the terrain built from it (D3).** The alignment the
sub-valleys engrave is the D8 LATTICE's at this roughness (D1), and a moving average leaves a line on
its lattice direction.

## E — the canyons' cone count is the NORM of a deep floor: "intersection" is refuted as its cause

Finding 125's instrument is reproduced exactly: 168 / 170 / 39 / 34 cones at the four floors, the same
classes. The reference is 200 floor cells drawn at random (with replacement, fixed seed) per
population, counted with the same instrument (skeleton and masks on PRE).

| reference | polylines | p10 | p50 | p90 | p99 | max | the four canyons' ranks |
|---|---|---|---|---|---|---|---|
| **témoin C2/10 col, trunk floors** (the asked one) | 636 | 1 | 1 | 3 | 5 | 6 | at the same four cells, the témoin's own skeleton counts **1 / 2 / 2 / 3** (p0 / p60 / p60 / p82) |
| B2 → A_c, trunk ≥ 10 floors (the canyons' own world) | 103 080 | 2 | 14 | 119 | 232 | 277 | 170 → **p94** · 39 → **p73** · 34 → **p70** |
| B2 → A_c, sub < 1 floors | 103 080 | 3 | 23 | 104 | 209 | 268 | 168 → **p96** |

- Against the témoin, 34–170 lies far beyond its p99 (5), but trivially: B2 → A_c has 162× the
  polylines. At the canyons' own cells the témoin's skeleton counts 1–3 cones, an ordinary floor.
- **Against the canyons' own world, none of the four lies beyond the p99.** The two largest canyons
  (bodies 13 and 14, 15 and 42 km² domain) sit at p70–p73: ordinary trunk floors. "Intersection" as
  this count measures it is REFUTED as the cause; the count reads a deep floor, not a crossing.
- **Where to look instead.** An over-dug body is a CLOSED depression: its floor is lower than its sill.
  The instrument has read the floor twice (F125, F126). What makes the depression closed is the SILL,
  the lowest point of its rim downstream. Why that sill stayed high at A_c and not at 1 km² is the
  unmeasured question: the floor profile along the D8 path downstream of each canyon floor, in B2 → A_c
  against B2 ≥ 1 km², up to the first cell that rises above the floor.

## The clause table, updated: each "supported" line proved or refuted

| unwanted effect | Finding 125's status | what Finding 126 measured, and where | **Finding 126's status** |
|---|---|---|---|
| **B1 spurs** | mur ↔ mer: supported, not proven (no base rate) | base rate ×8.3 (B1a) / ×12.5 (B1) near a coastal wall, ×1.0 elsewhere. The clause (walls ≥ sea + 0.5 m) returns the rate to the témoin's in three worlds, and 0 new spurs | **PROVED** (seed 1, 8192²) |
| **Finding 38** in B2 (A_c) + B1 | "the sea's threshold on a lowered wall": consistent, coastal unmeasured | two COASTAL wall cells (0.10–0.15 km domain from the open ocean), drowned by the CONSTRUCTION (foot −1.95 / −20.28 m + noise −4.16 / −32.43 m, additive), reposed at −1.00 m by `shelf_min_depth_m.max(1.0)` | **PROVED; the SAME clause** (F38 holds with it). Not a separate rule: mur ↔ mer at an enclosed cell |
| **B1 teeth** | "bruit ↔ seuil" conditioned on the profile: a candidate | under B1 ONLY, in the concave FOOT (mouths 89.6 %, heads 59.3 %, points 2.26 × the foot's area share); neither half alone (×1.29, ×1.10); mur ↔ mer leaves them unchanged (+1.8 %) | **LOCATED, clause not built** (one clause per round). The candidate is named: *bruit ↔ pied*, no detail put back on the concave foot, or the foot's slope floored above the detail's |
| **B2 terrain R8** × 2.3 | "hiérarchie ↔ tracé" (replacing hiérarchie ↔ planéité) | the 1–3 km² branches sit AT the D8 floor for the terrain's roughness (0.601 against 0.56–0.61, D1); no moving average lowers them beyond 5 % (D2); B2 on the 1 W-smoothed skeleton keeps × 2.34 (D3) | **"hiérarchie ↔ tracé" as a moving average: REFUTED** (D2, D3). The alignment is the D8 LATTICE's. The candidate is named, not built: a tracé OFF the lattice (steepest descent integrated in continuous space); not MFD (≡ D8 at 93 %, Finding 15) |
| **B2 → A_c canyons** (4) | "fond ↔ fond" at a trunk: undecided | the cone count is p70–p96 of its own population: the norm of a deep floor | **"intersection / mur ↔ mur" REFUTED as the cause** (by this count). Where to look: the SILL, not the floor |

## Predictions scored

**Mine** (written before the round):
- **P0 HELD.**
- **A.** **P-A1 HELD** (`shelf_min_depth_m`, inside the build). **P-A2 HELD** (below the sea at the
  construction). **P-A3 HELD with a correction**: a sum, not a product.
- **B.** **P-B1 HALF**: the témoin's walls ×1.2 held; B1a ×3–6 refuted, it is ×8.3; "B1a + B2 ≈ B1a"
  refuted. **P-B2 HELD except its "not ×1.5"**: the clause reached ×0.67.
- **C.** **P-C HALF**: "the mouth at the foot > 70 % in all three worlds" is refuted as binned (46 / 50
  / 90 %). It holds only in substance: floor + foot ≥ 95 % everywhere, because the mouth is biased by
  construction. "The body moves toward the foot" held.
- **D.** **P-D1 HELD** (at the floor). **P-D2 HALF**: it does apply to the sub-valleys (held), but
  "falls ≤ 0.45 at 1 W" is refuted. **P-D3 HALF**: the R8 ×1.5 and the halving are refuted; teeth,
  skeleton and trunks held.
- **E.** **P-E HELD.**
- **Meta ("at least two wrong") HELD.**

**Reviewer's** (in the prompt):
- **A HELD**: `shelf_min_depth_m`, coastal sub-valley walls, reposed at −1.00, enclosed. "Product" →
  sum.
- **B**: témoin ~0.1 /km **REFUTED** (0.0142); B1a ~×8 **HELD** (×8.3); the clause to the base rate
  within ×1.5 **HELD** (×0.67 / ×0.89); F38 → 0 **HELD**; teeth ±5 % **HELD** (+0.7 / +1.8 %).
- **C**: > 70 % of B1's teeth at the foot **HELD** (89.6 %); B1b uniform **HELD for the points**
  (×1.10), not for the mouths.
- **D1** floor ~0.45 / "the tracé, reducible" **REFUTED** (0.56–0.61, at the floor).
- **D2**: "250 m only on trunks ≥ 10" **REFUTED** (82 % of the 1–3 km² chords); 1 W → ~0.50
  **REFUTED** (0.568).
- **D3**: ×2 → ×1.3 **REFUTED** (×2.34); teeth −80 % **HELD** (−81 %); skeleton 88 % **HELD**
  (90.8 %).
- **E**: 170 beyond the p99 **REFUTED** on the fair reference (p94).
- **Meta HELD.**

## Limitations, stated

1. **A: the bench's direct carve on S1 is not bit-identical to the pipeline's S2** (203 352 cells,
   0.30 %). The cause I suppose is S1 carrying the C-2 rims; it is not verified. Only the two cells and
   the census, which are identical, are read.
2. **A: my Euclidean-nearest recomposition fails at both cells.** The carve's value there comes from
   another polyline, so they are junction cells. Only the switched carves are read.
3. **A: the enclosed components the construction lays below the sea in the non-F38 worlds** are counted
   at the construction stage only.
4. **B/C use Finding 125's instrument** (the skeleton and masks on PRE, not S1). D2 shows PRE's and
   S1's skeletons equal within 0.001; the masks' own agreement is not measured.
5. **B: a spur is typed inlet or peninsula by ONE probe cell.** The coastline length is the contour's
   segment length; "near" is a 2 km (domain) threshold.
6. **C: the positions use the Euclidean-nearest sample**, which item 2 shows can be the wrong one at
   junctions. At A_c the "floor" bin reads 30 % of the wall's area, so C is read on the 10 km² worlds
   only.
7. **D1: the synthetic's spectrum is not PRE's, and one scalar ratio matches them.** The floor is
   interpolated between two levels over a steep interval (0.80 → 0.36). "At the floor" carries that
   width. The first roughness instrument was contaminated by the sea and is not read.
8. **D2's chord band is read at the chord's middle sample.** The "sub < 1" rows of a 1 km² skeleton are
   Finding 125's off-trunk artefact.
9. **E: 200 cells per population, drawn with replacement.** The sill is not measured.
10. **One seed, one resolution (8192²)**, as for every bench of this chantier since Finding 121.

## Uncommitted state (awaiting the author)

- `valley_construction.rs`:
  - two gated fields, `wall_sea_floor_m` (the clause) and `smooth_w` (the skeleton variant), `None`
    by default;
  - `smooth_positions`, extracted byte-identically;
  - two permanent tests with negative controls.
- `cached_product.rs`: the key test lists the two new fields.
- `tests/f126_coast.rs`: this round's bench.
- `docs/reports/c1_continental_buoyancy/f126_coast/`: this report.

**Nothing promoted.**
