# Finding 141 — F140's crash attributed (nothing corrected): body 4 is not lake 4 but half of the merged below-sea lake 1000011; the construction splits that lake into two bodies that retreat apart, and lays the drained half (body 26) at its land fringe's floor, 0.1 m; that floor's cones reach into body 4 (the cell 78 m under its floor), and the breach carves ramps from the drained fringes' pits, anchored on the pit's floor and running to the flood's root, below the sea for up to 416 cells

**Status: diagnosis only. Nothing is corrected, nothing is promoted.** Only point 0.1 is committed (`f6c1caf`).

**Raw outputs, in this folder:**
- `f141_attr_raw.txt`: bench `f141_attr` (Parts A, Br and Z);
- `f141_b4_raw.txt`: bench `f141_b4` (body 4's identity, declared as an amendment after Part A);
- `checks_before.txt`, `checks_after.txt`.

The predictions were written before any measurement (`f141_predictions.md`). The instruments were declared before
the measurements (`f141_declared.md`, including its amendment for body 4).

## Partie 0

1. **The F140 code is committed with the viz toggle HIDDEN** (`f6c1caf`, `GORGE_TOGGLE_VISIBLE = false` in
   `workspace.rs`). The author's world cannot reach the crash from the viz; the gate stays reachable from code and
   benches. The docs and the F140 diagnostic benches went in the same commit.
2. **Added to ADR F140** (in that commit):
   - deviation 1 changed F139-R's tabulated quantity, and the rule "a deviation that touches a tabulated quantity is
     re-tabulated before the construction";
   - the crash will come back at old ages even with A corrected;
   - ON's 292 below-sea land cells, covered by its lakes: F38 holds by covering, not by construction;
   - option 2 is an author's decision, not a fix.
   - **This round corrects the first item's premise** (Part A below): body 4's A is not a deviation-1 artefact.
3. **Checks before and after, the gate off**: the guard 6 / 6 field and lakes (C2 /10 col = `a8d2d538d692c2f0`, the
   reference in `crates/ymir-core/data/bench_field_hashes.json`), lib 593, viz 30, `cargo check` clean
   (`checks_before.txt`, `checks_after.txt`).

## D — the five deviations of spec v3 and the tables they touch (read, not measured)

| # | the deviation | what it changes | why | tabulated quantities it touches | to re-tabulate before building |
|---|---|---|---|---|---|
| 1 | A = the construction's drained area at the col (skeleton D8), not the world's accumulation + inflow | r_lake = min(2, r·(A/A_ref)^p), so the level L(r_lake) and the emptying age. Also S_loi(A), so the gorge's slope m·S_loi, D_g, the head fall φ·D_g and S_req's comparison. And **A_ref = 418.7 km² was measured with the OTHER definition** (the median of (ii) on the 14 D8 lakes) | the construction does not know the final world | F139-R (r_lake, the areas by r, the emptying ages, "0 lakes empty at ×1", −55 % at ×1); A_ref itself; F138-T1 (S_g, fit / no room); F139-C2 (m_fit's denominator), C3 (0 shortage falls), C4 (fit / no room); F138-T4 (head falls, H_f) | R with the chosen A **and A_ref recomputed with the same A**; C2–C4 and T4 with S_loi of that A |
| 2 | L_bed = the construction's minimax from the lowest cell to the outflow, on the input, not "today's level" | the r ≤ 1 branch's end and the r > 1 branch's start of L(r), so every level at 0 < r < 2 | the same object (F136-K7: today's level is the bed sill), computed upstream | F138-T3 (lbed = min(ring, today's level)), F139-C1 (T3 at r = 1 against the final area, 0.982–1.000), F139-R's areas | T3 / C1 with the construction's L_bed: the D8 lakes' L_bed against today's level, lake by lake |
| 3 | φ = splitmix64 of the body's lowest cell, not of (seed, lake id) | which lakes have no head fall (1/3), φ per lake, so the head fall φ·D_g, the start L' of the gorge and the tag counts | the construction has no seed and no final lake id | F138-T4 (17 head falls, 11 tagged, median 26.2 m, max 199.5 m on lake 2); F139-C4 (lakes 2, 3, 9 fit thanks to φ 0.25 / 0.25 / 0.34; lakes 1, 10, 1000016, 1000019 at φ = 0) | T4 / C4 with the construction's φ (the law is the same, the draws are not) |
| 4 | the 3 below-sea lakes without a spillway segment cannot be recognised; their depression components are treated like the others; 1000001 untouched | these 3 get a retreat and a gorge the spec excluded | the construction does not know which final lakes will have a spillway | F139-C5 (8 of 9 within 1 m); F138-T3's "22 lakes with a body" totals. **F139-R tabulated only the 14 D8 lakes**: no below-sea lake's retreat was ever tabulated | R for the 8 spillway lakes and the 3 without a spillway |
| 5 | a lake whose outflow is not on a trunk (A < 10 km²) gets no gorge sample: its col is not cut and it keeps L_in at any r | the level stays at L_in. Its tributaries are still based on L(r_lake) by the χ walk (the lake stop does not test the trunk) | the gorge is laid by the trunk samples only | F139-R and G-area (R let every one of the 14 D8 lakes retreat) | the list of lakes with A(col) < 10 km² under the chosen A, and R without them |

**Three more differences, not among the five declared.** They are differences between F138 / F139's T3 bodies and
the construction's `gorge_bodies`: (a) and (b) were found while reading, (c) by Part A.
- (a) **the construction's bodies hold land cells only** (`land[s0] && dep(s0)`, `valley_construction.rs` in
  `gorge_bodies`). T3's bodies also held the inland below-sea cells (`inland` in `f139_t`). A below-sea lake's body
  is therefore its above-sea fringe, and its L_floor is the fringe's lowest land cell, not the basin's floor.
- (b) **L_floor is read on S1** (`z(low)` on `field`) for every body. T3 read it on the breached input for the
  depression bodies (`zbf`).

- (c) **a merged final lake is several construction bodies, each retreating on its own** (`f141_b4`).
  - 1000011 (merged at 459.3 m) is body 26 (the basin's land fringe) and body 4 (input lake 5), with separate A,
    r_lake and level.
  - The spec's « L_in = the merged spill » for the below-sea lakes (F139-C5) is not what the construction sees:
    body 4's L_in is 325.8 m.

All three touch F138-T3 / F139-R's areas and levels for the depression bodies and the below-sea lakes. They are to
re-tabulate with deviations 2 and 4. **(a) and (c) are, with the breach's ramp, the causes of F140's crash** (below).

## A — A by three definitions (22 lakes with a body; p = 0.5, A_ref = 418.7 km²)

**The set**: ON ×1's final lakes ≥ 1 km² with a D8 outlet or a spillway path (F139's set): 14 D8, 8 spillway.
1000001 has no construction body. Each lake's body is the GORGE skeleton's body with the largest overlap.

| lake | body | (i) col | (ii) F139 | (iii) inflow | (iii)/(ii) | (i)/(ii) | r_lake ×1 (i / ii / iii) | emptying age (i / ii / iii) |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 1 371.7 | 1 295.9 | 2 035.8 | **1.57** | 1.06 | 1.81 / 1.76 / 2.00 | ×1.04 / ×1.05 / ×0.97 |
| 2 | 1 | 442.7 | 357.8 | 491.7 | 1.37 | 1.24 | 1.03 / 0.92 / 1.08 | ×1.38 / never / ×1.34 |
| 3 | 2 | 598.7 | 698.6 | 604.1 | 0.86 | 0.86 | 1.20 / 1.29 / 1.20 | ×1.27 / ×1.22 / ×1.27 |
| **4** | **3** | 300.7 | 298.1 | 305.2 | 1.02 | 1.01 | 0.85 / 0.84 / 0.85 | never (all three) |
| 7 | 6 | 244.7 | 240.4 | 244.4 | 1.02 | 1.02 | 0.76 / 0.76 / 0.76 | never |
| 8 | 8 | 156.4 | 153.7 | 199.2 | 1.30 | 1.02 | 0.61 / 0.61 / 0.69 | never |
| 9 | 10 | 393.3 | 389.0 | 393.4 | 1.01 | 1.01 | 0.97 / 0.96 / 0.97 | never |
| 10 | 9 | 819.4 | 816.7 | 819.5 | 1.00 | 1.00 | 1.40 (all) | ×1.17 (all) |
| 11 | 12 | 947.1 | 932.9 | 1 099.4 | 1.18 | 1.02 | 1.50 / 1.49 / 1.62 | ×1.13 / ×1.14 / ×1.09 |
| 12 | 14 | 425.1 | 428.8 | 425.0 | 0.99 | 0.99 | 1.01 (all) | ×1.39 (all) |
| 13 | 15 | 983.3 | 988.4 | 983.3 | 0.99 | 0.99 | 1.53 / 1.54 / 1.53 | ×1.12 (all) |
| 15 | 17 | 416.6 | 418.7 | 417.9 | 1.00 | 0.99 | 1.00 (all) | never / ×1.40 / never |
| 19 | 22 | 727.2 | 1 437.8 | 1 294.6 | 0.90 | 0.51 | 1.32 / 1.85 / 1.76 | ×1.21 / ×1.03 / ×1.05 |
| 20 | 23 | 227.2 | 240.3 | 323.2 | 1.35 | 0.95 | 0.74 / 0.76 / 0.88 | never |
| 1000011 | 26 | 2 264.3 | 2 409.6 | 4 304.0 | **1.79** | 0.94 | 2.00 (all) | ×0.96 / ×0.95 / ×0.89 |
| 1000012 | 29 | 161.3 | 308.0 | 313.7 | 1.02 | 0.52 | 0.62 / 0.86 / 0.87 | never |
| 1000013 | 28 | 448.8 | 1 263.0 | 1 231.6 | 0.98 | **0.36** | 1.04 / 1.74 / 1.72 | ×1.37 / ×1.06 / ×1.07 |
| 1000014 | 30 | 1 061.3 | 2 986.6 | 3 047.2 | 1.02 | **0.36** | 1.59 / 2.00 / 2.00 | ×1.10 / ×0.92 / ×0.92 |
| 1000016 | 31 | 1 805.1 | 3 608.0 | 3 594.7 | 1.00 | 0.50 | 2.00 (all) | ×0.99 / ×0.90 / ×0.90 |
| 1000017 | 32 | 437.5 | 927.9 | 1 535.8 | **1.66** | **0.47** | 1.02 / 1.49 / 1.92 | ×1.38 / ×1.14 / ×1.02 |
| 1000018 | 33 | 338.3 | 701.1 | 709.2 | 1.01 | **0.48** | 0.90 / 1.29 / 1.30 | never / ×1.22 / ×1.21 |
| 1000019 | 34 | 393.1 | 1 028.8 | 1 169.1 | 1.14 | **0.38** | 0.97 / 1.57 / 1.67 | never / ×1.11 / ×1.08 |

(km²; r_lake at ×0.7 / ×1.2 / ×1.4 and every row in `f141_attr_raw.txt`.)

- **(iii) is within ×1.5 of (ii) for 19 of 22.** The three outside are lake 1 (1.57), 1000011 (1.79) and 1000017
  (1.66).
- **(i) is beyond ×2 of (ii) for 5 of 22, all below-sea spillway lakes, and always SMALLER**: 1000013 0.36, 1000014
  0.36, 1000017 0.47, 1000018 0.48, 1000019 0.38. 1000012 (0.52) and 1000016 (0.50) are at the edge.
- **On the 14 D8 lakes (i) agrees with (ii)**: 0.86–1.24, except lake 19 (0.51).
- **Lakes drained (r_lake = 2) at ×0.7 / ×1 / ×1.2 / ×1.4**:
  - (i): 0 / 2 / 7 / 13;
  - (ii): 0 / 3 / 11 / 15;
  - (iii): 0 / 4 / 11 / 15.
  - **At ×1, every definition drains at least two lakes, all below-sea lakes**: 1000011 and 1000016 under (i);
    1000014 also under (ii); lake 1 also under (iii). F139-R's "no lake empty at ×1" covered only the 14 D8 lakes.

### Lake 4 and body 4: F140 paired them by index

- **F139's lake 4 is body 3.** Its A agrees under all three definitions (300.7 / 298.1 / 305.2 km²), r_lake 0.85 at
  ×1, and it never empties within the range.
- **Body 4 is half of the final lake 1000011** (`f141_b4`):
  - 1000011 is the merged below-sea basin lake, level 459.3 m, 493.2 km²; its footprint holds all of body 4's
    21 991 cells;
  - body 4 is the input lake 5 (L_in 325.8 m), inside the same merged water body;
  - **its A is the same by (i) and (iii): 1 571.0 and 1 571.3 km²**. Its inflow comes entirely from body 26 (the
    largest inflows 654.9, 318.4, 191.3, 78.8 and 18.1 km² all come from body 26). The exit carries 1 568.2 km².
- **F140's "1 571 against F139's 298 for the same lake" was wrong**: the 298 is lake 4 = body 3. Body 4's lake
  (1000011) has (ii) = 2 409.6 km². **Deviation 1 does not make body 4 drain**: option 1 (A = the inflow) gives it
  1 571.3 km², the same r_lake 1.94.

### The gap's stream, in the five lakes beyond ×2

- **No outside stream joins at the col.** In all five, the col's only significant donor is the body's own exit.
- **(i) is one river**: the body's largest single inflow, crossing the body to one exit.
  - 1000014: col 1 061.3 km² = its largest inflow 1 061.3 km² (35 % of (iii)).
  - 1000013: col 448.8 km², largest inflow 441.4 km².
- **The cause**: the five bodies are **the land fringes of below-sea basins**. The construction's bodies hold land
  cells only (difference (a) in D). Their inflow leaves the fringe at many points into the below-sea core, which is
  not in the body. The "exit" (the largest) takes one share.
- **Every one of the five cols lies INSIDE the final lake** (1000013 to 1000019). ON's accumulation there is
  0.0–0.6 km², except 1000018 (342.8 km²). The final outlet is 0.0–15.7 km away (1000014: 15.7 km).
- **The gorge of a below-sea lake is therefore laid from the fringe into the basin's own core, not at the merged
  lake's spill.**

## Br — the breach's below-sea land cells (instrumented copy, bit-identical to production in both worlds: 0 cells differ)

### The code (`terrain/flow.rs`, `breach_monotone_protected`)

- **The bases** (1037, 1049–1054): every cell with `z ≤ sea` (the ocean AND every inland below-sea cell) or a
  pre-breach lake. Lake cells are first set to their fill level (1032–1036).
- **The flood** (1055–1081) pops the lowest cell and records `backlink[nb] = ci` (1066): the tree points back to the
  base the flood came from.
- **The ramp** (1067–1077): when `height[nb] < z[ci]`,
  - `target = height[nb] − EPS` (1070): **the pit's FLOOR**;
  - it walks `ci → backlink → …` (1072–1075), setting `z = target` and lowering the target by EPS (1e-5 norm,
    0.113 m) at each step;
  - it stops only at a base, at a cell already lower than the target, or at the chain's end (1072).
- **The pit's spill (`filled`) is never read by the ramp.** It is read only to flatten the lakes (1034). The fill
  mop-up (1089–1118) cannot raise a cell ≤ sea: it is a base, seeded and never visited.
- **Answer to the central question: the breach does not carve towards each pit's spill.** It carves from the pit's
  floor towards the root of the flood that reached it (an inland below-sea basin, a cell an earlier ramp already
  lowered, or a lake), losing 0.113 m per cell. A long chain therefore crosses the sea before reaching its root.

### ON ×1 (the gorge off): 292 land cells ≤ sea

- **All 292 are laid by a ramp, all inland, 0 ocean-connected.** 17 lie on a construction body (31, 1000016's).
- **Where**: around (2 710–2 960, 5 205–5 270), in 102 ramps.
  - Pits 7.0–12.7 m high, spills 13.5–16.9 m: 91 of the 102 pits have a spill > 0.5 m above their floor.
  - Ramps of 61–140 steps; each crosses the sea at step 61–111.
- **The base**: 101 ramps (291 cells) stop at a cell an earlier ramp had already lowered ≤ sea; 1 ramp stops at an
  inland below-sea cell. **The flood chain's root is the inland below-sea cell (3057, 5288)**, 1000016's basin.
- **Why below the sea although the spill is higher**: the ramp is anchored at the pit's floor (7–13 m), not at its
  spill. At 0.113 m per step it reaches the sea after ~60–110 cells, short of its root.
- **Covered by a final lake: 292 of 292 (100 %), all by lake 1000016**; 0 by a pre-breach lake.

### GORGE ×1 (r 1, p 0.5): 2 821 land cells ≤ sea

- **All laid by a ramp, all inland. 2 670 lie on a construction body.** 0 are covered by a pre-breach lake.
- **143 ramps.** Their pits are mostly **the drained below-sea fringes**, at their L_floor ≈ the sea: floors
  0.0–2.5 m, spill = floor (nothing held).
  - The top twelve ramps lay 1 795 of the 2 821: body 26 (1000011) 1 070 cells in 6 ramps; body 32 (1000017) 321;
    31 (1000016) 185; 28 (1000013) 111; 30 (1000014) 52; no body 56.
  - **The largest ramp**: pit (4732, 3663), floor 0.1 m. It runs 417 steps to the inland below-sea cell
    (4451, 3861), crossing the sea at step 1: **416 cells ≤ sea, the deepest about −47 m.**
- **The base**: 20 ramps (1 391 cells) reach an inland below-sea cell; 120 (1 202 cells) reach a cell an earlier ramp
  lowered; 3 reach a pre-breach lake.
  - Lake 5 at (4162, 3604): 151 cells, from pit (4145, 3675), which is 3 cells from Z.
  - Lake 3 at (4379, 3550): 74 cells, from pit (4299, 3545), floor 5.0 m, spill 178.7 m. **It stops one cell from
    F38's floor cell (4379, 3551).**
  - Lake 20: 3 cells.
- **Coverage**: no final lake map (the tail panics).
  - The assert names 2 uncovered components, 121 + 1 cells. **So at most 122 of the 2 821 are uncovered.** This is
    deduced from F140's panic message, not measured.
- **2 821 against 292**: the gorge makes the below-sea fringes flat at ≈ the sea. Every pit there starts its ramp at
  or under the sea, so ramps cross the sea at step 0–44 instead of 61–111, over chains up to 417 cells.

## Z — the cell (4148, 3674), 78 m under body 4's drained floor

| stage | z (m) | body 4's lowest (m) | body-4 cells > 1 m under L_floor (168.9) |
|---|---|---|---|
| S1 (the input) | 325.6 | 168.9 | 0 |
| `carve_diag(S1)`, ON skeleton | 325.6 | 168.9 | 0 |
| **`carve_diag(S1)`, GORGE skeleton** | **89.4** | **0.1** | **104** |
| (a) the pipeline's construction + rims | 89.4 | 0.1 | 104 |
| (b) + the light pass | 91.1 | 10.5 | 180 |
| (c) + the droplets | 91.1 | 10.5 | 180 |
| (d) + the bathymetry | 91.1 | 10.5 | 180 |
| (e) the protected breach | **−0.2** | −17.0 | 317 |

- **The construction lays it**, and the stage is `carve_diag` (`valley_construction.rs`): the nearest-sample floor
  (`geo`, 1841–1861, read at 1976).
  - The cell's own skeleton base is body 4's L(r) = 178.7 m (χ 0). The rim clamp does not apply (NaN).
- **The laying sample**: line 189, position 172 of 611, at (4148.5, 3675.5), 49 m away, inside the half-width
  (115 m), so on the floor: **zf 89.4 m**.
  - **The sample's cell (4148, 3675) belongs to body 26** (base 0.1 m, χ 0): 1000011's land fringe, at r_lake 2.00
    and level L(2) = its L_floor ≈ 0.1 m.
  - The line runs from zf 13.4 m to zf 178.7 m (the floor RISES downstream: body 26 drains into body 4), with 436 of
    its 611 samples on body-4 cells and a min zf of 0.1 m.
  - **The sample's zf (89.4 m) lies between body 26's floor (0.1) and body 4's (178.7)**: `densify` (1290–1309)
    interpolates zf between consecutive samples.
- **The chain**:
  1. `gorge_bodies` (897–1034) makes body 26's level its land fringe's lowest cell (`l_floor = z(low)`, 969; the
     level at r 2, 1016).
  2. The χ walk bases body 26's cells there (630).
  3. `line_samples` (1266) gives the line's samples zf = 0.1 m on body 26.
  4. `densify` blends them with body 4's 178.7 m.
  5. `carve_diag`'s nearest-sample cone lays body 4's cells at 0.1–89 m.
- **Not the light pass**: it raises the cell by 1.7 m and body 4's lowest from 0.1 to 10.5 m.
- **Not the breach**: it takes the cell from 91.1 to −0.2 m later, through the ramp from pit (4145, 3675).
- **The function that lays it below L_floor: `carve_diag`, from a body-26 sample whose floor `gorge_bodies` set at the
  fringe's lowest cell.**

## The mechanism of F140's crash, attributed

1. **The construction splits one final water body into two bodies.** 1000011 (merged at 459.3 m) is body 26, the
   below-sea basin's land fringe, and body 4, input lake 5. Each retreats with its own A and r_lake. Body 26 drains
   into body 4.
2. **At ×1, body 26 is drained** (r_lake 2.00, all three definitions) to its fringe's lowest cell, ≈ 0.1 m. Body 4
   is at r 1.94 (178.7 m), with the same A under (i) and (iii).
3. **Body 26's floor reaches body 4 through the line that crosses both** (densified samples, nearest-sample cones):
   104 cells of body 4 under its floor after the construction.
4. **The other drained below-sea fringes** (bodies 31, 32, 28, 30) are laid at ≈ the sea the same way.
5. **The breach carves ramps from those pits** (floors 0–5 m, nothing held), anchored on the floor and running to the
   flood's root. 2 821 cells end ≤ sea. Some components reach no water body. One of them, F38's floor cell
   (4379, 3551), sits one cell from the end of a 118-step ramp started in body 26.

**Deviation 1 is not the cause**: (i), (ii) and (iii) all drain body 26 at ×1, and give body 4 the same A. The causes
named here:
- **(α)** the bodies do not follow the final merged lake;
- **(β)** a below-sea lake's body is its land fringe, so "drained" means "at the fringe's lowest cell, ≈ the sea";
- **(γ)** the breach's ramp anchored on the pit's floor, which already makes ON's 292 (covered by 1000016).

## Predictions

**Mine** (`f141_predictions.md`):
- **P-D**: "deviations 1 and 3 touch tabulated quantities, and 2" **held**, but "two others besides 1" **undercounted**:
  2, 3, 4 and 5 all do. Three more differences, not declared, also touch them.
- **P-A**:
  - (iii) ≈ (ii) for ≥ 15: **held** (19);
  - (i) beyond ×2 for ≥ 4: **held** (5);
  - **"lake 4 among them" refuted** (lake 4 = body 3 agrees at 1.01; body 4 is part of 1000011, with (i) = (iii));
  - **the cause (a river the input footprint intercepts) refuted**: there is no outside stream at any col. The gap is
    the fringe's exits.
- **P-Br**:
  - "the ramp starts at the pit's floor, not its spill, and stops at a base or a lower cell": **held**;
  - "ON's pits low (≤ ~10 m) near the below-sea basins": **held in place, the bound loose** (7.0–12.7 m, 1000016's
    basin);
  - "ON ≥ 90 % covered, mostly by the merged below-sea lakes": **held** (100 %, 1000016);
  - "GORGE's from the drained bowls' pits": **half-held**: from the drained below-sea FRINGES (2 670 on bodies,
    body 26 first), not body 4's bowl.
- **P-Z**: "the light pass" **refuted**: the construction lays it. The second choice ("a cone of the outlet line's
  downstream samples") is **refuted in its detail**: it is an UPSTREAM drained body's samples (body 26 drains into
  body 4).
- **Meta** ("at least one of mine is wrong"): **held**.

**The reviewer's:**
- **A**: (iii) ≈ (ii) for ≥ 15 / 22: **held** (19); (i) > ×2 for ≥ 4: **held** (5).
- **Br**: "the ramp targets a base lower than the spill": **held** (the flood's root: an inland below-sea cell, an
  earlier-lowered cell, a lake). "The 292 > 90 % covered": **held** (100 %).
- **Z**: "laid by the construction (a cone of a valley line crossing the cuvette), not the breach": **held** (line 189
  crosses bodies 26 and 4).
- **D**: "at least one deviation other than 1 touches a tabulated quantity": **held** (four do).
- **Meta** ("at least one is false"): **refuted**, none is.

## What this changes in F140's record (stated, not corrected in F140's text)

- **"Body 4: A = 1 571 km² (the F139 bench: 298 km² for the same lake)" is wrong.** The 298 is lake 4 = body 3.
  Body 4 belongs to 1000011, (ii) 2 409.6 km².
- **"The cause is deviation 1" is not supported.** Body 4's A is the same under (iii), and body 26 is drained under
  all three definitions.
- **Option 1 (A = the inflow) would not have prevented the crash**: by (iii), body 26 drains and body 4 keeps r 1.94.
- **The added Partie 0 item stands**: lakes empty by design at old ages, and drained fringes are what the breach
  carves.

## Limitations, stated

1. **The GORGE coverage** is deduced from F140's panic message (≤ 122 uncovered), not measured on a final lake map.
2. **Body 26's own levels** (L_in, L_bed, L_floor) are not printed. Its level ≈ 0.1 m is read from its cell's
   skeleton base.
3. **That the sample's zf comes from `densify`** is read from the code and its value (between 0.1 and 178.7), not
   traced sample by sample. `smooth_positions` (moving a body-26 sample over a body-4 cell) is the other path; both
   are inside `line_samples` → `carve_diag`.
4. **The per-body totals of the 2 821** are given for the top twelve ramps only (1 795 cells). 2 670 are on a body.
5. **ON / GORGE are one world, one seed** (the témoin).

## State

**Uncommitted**, awaiting the go-ahead:
- `f126_coast.rs`: `f141_attr`, `f141_b4`;
- this folder;
- ADR Finding 141.

**No code changed in `src`.**
