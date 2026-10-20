# Finding 186 — the deposition by closures: Yuan's K_eff holds the equilibrium at +0.4 %; the per-cell marine step proposed as the physical reference, to be confirmed by the author (the blocks hold a 60 s budget, since withdrawn, but lose the shelves, stripe the delta and multiply the deep-sea lines); the burial by a level removes the stairs, but where the chain lifts or raises the basement it overshoots by 1.4–3 km; the recycling interpolated; the chain to 8 192² stays ~3 × over the 5-minute target (F187)

**Status: Part 0 measured on mains; A and B built with their controls; tier 1 (the panel, intermediate) and tier 2 (the chain to
2 048², seed 23) measured; tier 3 (the old age) not run: its condition is « tiers 1 and 2 hold their budgets », and the chain's
5-minute budget fails without the deposition already. Nothing committed.** Declared: `f186_declared.md`. Predictions:
`f186_predictions.md`.

## The rules (the author's)

- **The budget**: the complete chain, C1 to 8 192², ≤ 5 minutes of wall time per world (the author's words, verbatim, in the ADR).
- **The measurement**: every time budget is measured on mains power, on the author's reference machine (to be specified by the
  author). A battery measurement is redone before any decision. Part 0 was first measured on battery (≈ 2.5 × slower); it was
  redone on mains, and the battery figures are kept apart (`results/budget_seed23_battery.txt`, `f186_bench_output_battery.txt`).

## Part 0 — the chain's budget on mains (seed 23, intermediate, no deposition)

| stage | what it does | wall |
|---|---|---|
| the selection | 203 C1 runs to their end (`cascade/profile.rs:192`) | 108.8 s |
| C1 | the chosen continent, 300 steps (`cascade/history.rs:359`) | 1.2 s |
| the physics level 256² | 300 steps (`run_loop`) | 24.6 s |
| amplification 512² / 1 024² / 2 048² | iterations, no physical time (`amplify.rs:1398`) | 6.1 / 9.6 / 41.7 s |
| amplification 4 096² | **does not exist** for the continent profile; extrapolated (the 2 048² cost per cell-iteration × cells × `Budget::moyen`'s PROXY row) | ~154 s |
| amplification 8 192² | **does not exist**; extrapolated the same way | ~615 s |
| **the chain to 8 192²** | | **≈ 960 s against 300 s** |

- **W3 held**: without the deposition the chain is already ~3.2 × the target.
- **The most expensive**: the 8 192² amplification (~64 %), then the 4 096² one (~16 %) and the selection (~11 %). They are F187's
  (the author's ruling: neither those levels nor the selection are touched here).
- **0.2**: F184's live combination at 256² took 69.9 s against 24.6 s without the deposition; the marine step was 38.8 s of it
  (65 691 sub-steps).

## Part A — the K_eff closure, and the author's marine fallback

- **A1**: `StreamPowerConfig::keff_g` (`None` in production) divides K by 1 + G (A − 1)/A in the incision's loop, A the
  incision's MFD area in cells. `Deposition::fluvial_keff` sets it, and the pass deposits nothing fluvial: it routes, fills the
  lakes and delivers to the sea.
- **The analytic test with K_eff** (`the_yuan_equilibrium_keff`, without the hillslope): 2.529 against the discrete factor
  2.520, **+0.4 %**; every area bin within ~2 %. **W1 held.**
- **The controls**:
  - K_eff at G = 0 is G = 0 bit for bit;
  - both mass controls < 10⁻⁹;
  - the flags off hash 674074df3c49c6cd (as at 022d3cc and f6ec78e);
  - the block step conserves the solid exactly, and nothing goes above the sea but the stranded part.
- **The author's fallback** (`Deposition::marine_blocks`, declared before tier 1, used because the 60 s budget failed on the
  marine step):
  - the live marine step on 2 × 2 blocks (128²): the restriction (the sources summed, the mean room), F184 A1's rules on the
    blocks, the prolongation in proportion to each cell's room;
  - the stability limit at 128²: 16 538 sub-steps against 65 691;
  - **the justification**: √(K_M d) ≥ 2.7 km from the shortest step (18.7 kyr), against a 3.1 km block.

## Tier 1 — the panel at 256², intermediate (the circumnavigation against the same world without the deposition)

| | A per cell | A with the marine blocks |
|---|---|---|
| the physics level, wall per world | 70.1–81.1 s (**fails 60 s**; the marine step 43–53 s) | **26.9–34.2 s** (the marine step 2.7–3.8 s) |
| the sea kept | 5 / 5 | 5 / 5 |
| the circumnavigation lost | 0 | 0 |
| Δ median altitude (median) | +9 m | +11 m |
| Δ plains / Δ mountain (median) | −1.9 / +2.4 points | −1.4 / +2.9 points |
| E3 (the large lakes' area falls) | 1 / 5 | 2 / 5 |
| the passive shelf (median, per seed) | 2.2 / 6.9 / 8.2 / 37.0 / 3.1 km | 3.1 / 3.1 / 2.2 / 2.2 / 2.2 km |
| the sediment in the lakes / at sea | 0.3–0.7 % / 98.8–99.3 % | 0.2–0.7 % / 99.1–99.4 % |
| the mass controls (worst) | 8.5 × 10⁻¹¹ | 4.7 × 10⁻¹¹ |

- **The time saved by the blocks**: the marine step ×~14 faster, the physics level 70–81 → 27–34 s.
- **What the blocks change**: the shelf. Per cell the deposits build passive shelves of 7–37 km on 3 seeds; with the blocks they
  stay at 2–3 km. The marine deposit above 100 m is about the same (0.3–3.3 %).
- **The slopes in the deposition zones** (blocks):
  - the deposition zone covers 3.7–25 % of the land, and its median slope falls (0.25 → 0.03° on seed 23, 0.54 → 0.10° on
    seed 8) or stays.
  - The plains still fall there on 3 / 5 seeds (−1.8 to −9.8 points) while its slope falls. On seed 8 its median altitude rises
    59 → 367 m: a threshold of the Kapos altitude crossed, not a steeper relief.
  - In the erosion zone the slopes move little (−0.2 to +0.17°).
- **E3 fails** (the large lakes grow on 3 / 5). The reviewer's hypothesis, not measured: K_eff slows the sills' incision (the
  response time × (1 + G)), and A1's load flexure subsides the land behind the loaded coasts. Both make closed depressions last
  longer.
- **(v) and W₁₀ were not measured** in tier 1 (the adoption needs (v)).

## Tier 2 — the chain to 2 048², seed 23, intermediate, with A (blocks) and the transmission

- **The stages**:
  - C1 1.1 s; the physics level with A 34.5 s;
  - the amplifications 4.9 / 8.7 / 39.1 s (without the deposition 5.0 / 9.2 / 42.1 s);
  - **the burial 0.06 s**;
  - with Part 0's selection (108.8 s) and its extrapolation of 4 096² and 8 192² (~145 / ~577 s), **≈ 919 s against 300 s**.
    The deposition adds ~10 s to the chain; the overrun is the chain's own (F187).
- **At 2 048², against the same world without the deposition**:
  - plains 48.2 → 54.9 % (**W6 held**);
  - hills 32.6 → 23.4 %;
  - mountain 19.2 → 21.7 %;
  - the median altitude 73 → 67 m.
- **The cover's volume**: on every cell the fine cover is 134 337 km³ against the coarse 136 600 km³ (−1.7 %). On the land alone
  it is −18.8 % (10 496 against 12 929 km³), because the land / sea line moves between the resolutions.
  - (A first count compared the fine cover on every cell with the coarse cover on land only, +939 %. It was the bench's mistake,
    corrected before reporting.)
- **The plates** (`map2048_off_A_seed23.png`, `crops2048_lake_delta_plateau_off_A_seed23.png`):
  - the sea is shallower and the valleys are filled to a smooth floor (the « plateau » crop);
  - **the burial leaves 8 × 8-pixel stair edges** where the cover meets the basement (the coasts, the filled lakes): the rule
    reads the 256² parent's H > 0, one parent per 8 × 8 fine cells;
  - **the delta crop shows vertical stripes and fingers** in the shallow sea (not diagnosed: the blocks' prolongation, or the
    bicubic surface against the amplified basement).
  - These are for the author's eye.

## The predictions

- **The reviewer's**:
  - W1 held (+0.4 %);
  - W2 refuted per cell (70–81 s), held with the author's fallback (27–34 s);
  - W3 held (~960 s on mains);
  - W4 held (5 / 5, 0 lost);
  - W5 refuted (1 / 5 per cell, 2 / 5 with the blocks);
  - W6 held at 2 048² on seed 23 (+6.7 points of plains).
- **Mine**:
  - Y1 refuted (the test held within 1 %);
  - Y2 refuted (per cell the physics level took 70–81 s; ~40–50 s predicted);
  - Y3 refuted (the relief moved little: +9 to +11 m, not +50 to +100);
  - Y4 held (E3 fails);
  - Y5 held;
  - Y6 refuted on land (−18.8 %), held on every cell (−1.7 %).

## For the author

1. **The marine fallback** (2 × 2 blocks) holds the budget and leaves the relief, the sea and the circumnavigation as per cell, but
   it narrows the passive shelves the per-cell step builds (7–37 km → 2–3 km on 3 seeds). Adopt it, or another way to the 60 s?
2. **The burial's stair edges** (8 × 8 at 2 048²): a continuous cover mask (the bicubic of H, a threshold) would remove them but adds
   a parameter. Your call before tier 3.
3. **E3**: the large lakes grow on 3 / 5. A diagnosis (K_eff's sills against the load's subsidence, by switching A1 of F183 off) is
   one 256² run per seed.
4. **Tier 3** (the old age) was not run: the chain's 5-minute budget fails without the deposition (W3), so its condition « tiers 1
   and 2 hold their budgets » is not met by the chain. The physics level's 60 s holds.
5. **(v)** is still to measure before any adoption.

## The checks (after; the flags off)

- `cargo check --workspace`: clean.
- 672 lib tests pass (F184's 669, plus the K_eff controls, the burial and the block step); 33 viz tests pass.
- Production's C1 hash is unchanged.
- **The guard reads `NoReference` on all six states, but the output is unchanged.**
  - The cache key fingerprints the configs by their derived `Debug` (`cache.rs:134`, `with_debug`), and `StreamPowerConfig` now
    has one more field, `keff_g: None`. So every digest changed.
  - The 12 computed hashes are all present in the stored references: the six fields (`bench_field_hashes.json`: 14953e83…,
    b07eabbe…, c876c200…, a8d2d538…, c2eb1cb3…, 35e1bc81…) and the six lake pairs (`bench_lake_hashes.json`).
  - **Production's output is bit for bit the same.** A run_hd cache entry misses once and recomputes; nothing is stale.
- **For the author**: re-key the two reference files (the values unchanged), or take K_eff out of `StreamPowerConfig` (an incision
  argument instead of a field) so the digests stay as they are.

## The continuation (the author's five steps, in order; on mains; nothing committed)

### Step 1 — the guard re-indexed (F178's rule: the values bit for bit, so a re-recording, not a re-calibration)

- **The identity was checked before writing.** The six field hashes and the six lake pairs that the viz guard computes at the F186
  tree were compared with the stored references; all 12 matched. Only then were the two files' `digest` fields replaced
  (`bench_field_hashes.json`, `bench_lake_hashes.json`). No hash value changed.
- The guard after the round: **6 / 6 Match** (field and lakes; 135 s); production’s C1 hash test passes; 672 lib and 33 viz tests pass; `cargo check --workspace` clean (warnings only).
- **Queued (the author's)**: the cache key fingerprints the configs by their derived `Debug`, so it changes as soon as a field is
  added, even at `None`.

### Step 2 — the marine step: the exact gain first, then the blocks

**The active cells.** At each sub-step the sweep covered only the domain's cells that carry sediment or touch one (the transport
moves only the deposit, F182).
- The control held: the result is bit for bit the full sweep's (the live + K_eff configuration's hash 0c463b9f11e65660, before
  and after).
- **It is slower.** The physics level took a median of 88.3 s (max 94.3 s; the marine step 51–63 s), against 70–81 s for the full
  sweep. The sediment covers most of the sea, so the active set is almost the whole domain, and keeping it costs more than it
  saves. The sweep was reverted to the full one.
- **The rule therefore sends the step to the blocks** (more than 60 s per cell).

**The blocks with the author's prolongation.** The block's thickness change is interpolated bilinearly on the cells (the block
centres at 2b + 1, so weights ¾ and ¼; on the torus the total is kept exactly). A gain is capped by the cell's room under the
sea, and the surplus is shared among the cells of the block and its ring that still have room. The mass holds within 10⁻⁹
(the unit control, on the torus and off it).
- **The losses are not in the author's rule.** Two readings were measured on seed 23 and on the panel:
  - **t1d, the symmetric extension**: a loss capped by the cell's sediment, the excess taken from the ring;
  - **t1e, per block**: the loss taken from the block's own cells in proportion to their sediment.
  - t1d is kept: deficit lines 650 against 848 cells, E3 3/5 against 1/5.

| (the panel, intermediate) | per cell (t1) | blocks, first prolongation (t1b) | blocks, bilinear (t1d, kept) | blocks, bilinear, losses per block (t1e) |
|---|---|---|---|---|
| the physics level, wall per world | 70.1–81.1 s | 26.9–34.2 s | **median 29.9 s, max 33.1 s** | median 30.2 s, max 33.7 s |
| the sea kept / the circumnavigation lost | 5 / 5, 0 | 5 / 5, 0 | 5 / 5, 0 | 5 / 5, 0 |
| Δ median altitude, Δ plains / mountain (medians) | +9 m, −1.9 / +2.4 | +11 m, −1.4 / +2.9 | +7 m, −1.6 / +3.5 | +2 m, −0.4 / +2.9 |
| E3 (the large lakes fall) | 1 / 5 | 2 / 5 | 3 / 5 | 1 / 5 |
| **the passive shelf, median per seed (8, 18, 19, 23, 29)** | 2.2 / 6.9 / 8.2 / 37.0 / 3.1 km | 3.1 / 3.1 / 2.2 / 2.2 / 2.2 km | **2.2 / 2.2 / 2.2 / 2.2 / 2.2 km** | 2.2 / 3.1 / 2.2 / 2.2 / 2.2 km |

- **The shelves.** The bilinear prolongation does not bring back the shelves the per-cell step builds: 2.2 km on 5 / 5 seeds, the
  world-without-deposition value. So the « in proportion to the room » prolongation was not the only cause. What else in the
  block step narrows the shelves (the restriction's mean room, or the step at 128²) is **not measured**.
- **The polygonal lines.** The bilinear blocks leave faint 1-cell deficit lines in the deep sea at 256² (`relief_t1d_seed23.png`).
  - The index: deep-sea cells below 0.8 × their neighbours' median thickness. Seed 23: per cell 201, t1d 650, t1d without the
    recycling 245, t1e 848.
  - **Most of them come from the recycling's steps along the plate edges.** The per-cell diffusion, and the first prolongation
    (which sent the deposit to the deepest cells), fill them; a bilinear prolongation does not resolve a notch narrower than a
    block.

### Step 3 — the burial without steps, without a parameter (replaces B2)

- **`amplify::Transmit::bury`**:
  - on each covered 256² cell, the level L with Σ max(0, L − fine basement) = H × the fine cells' count, exact by sorting
    (`amplify::fill_level`);
  - L interpolated bilinearly on the fine grid, with L = the coarse basement where H = 0;
  - z = max(fine basement, L);
  - a covered lake stays at its flat bed.
- **My reading, to confirm**: the fill applies only on the fine cells whose interpolation stencil holds a covered cell. Without
  that, L = the coarse basement on every uncovered cell would raise every fine hollow of the map to its coarse mean, and the
  control « H = 0 changes nothing » would fail.
- **The controls** (`the_burial_closure`): H = 0 changes nothing; the level is exact (the fill equals H × the count to 10⁻¹²); a
  uniform cover adds exactly its volume.
- **The volume at 2 048²** (seed 23, step 5), fine against coarse:
  - every cell: −0.1 % (B2: −1.7 %);
  - the land alone: −3.4 % (B2: −18.8 %).
- **The previous implementation read the parent cell's H instead of the declared interpolated H.** That is where the 8 × 8 stair
  edges came from; they are gone from the step-5 plates.
- **The delta's vertical stripes are the marine blocks'.**
  - They are already in the bicubic of the 256² surface (`s5_crop2048_delta_bicubic_seed23.png`), so they come from neither the
    burial nor the bicubic.
  - At 256² they appear with all three block prolongations and not per cell (`delta256_cells_t1b_t1d_t1e_seed23.png`).
  - The 2-cell oscillation around the delta, along x / along y, on the sea cells within 20 cells: the surface with A 36.3 / 14.2 m,
    its thickness 53.2 / 22.5 m.
  - The block step's part responsible (the restriction or the step at 128²) is not isolated.

### Step 4 — the measures for the adoption (the panel at 256², the blocks of step 2)

| | without the deposition | with A (blocks) | with A, the load off |
|---|---|---|---|
| (v) at the intermediate age | 5 / 5 | 5 / 5 | 5 / 5 |
| W₁₀ global / active / passive (medians) | 0.150 / 0.354 / 0.090 | 0.149 / 0.230 / 0.093 | — |
| W₁₀ \|Δ\| max (global / active / passive) | | 0.073 / 0.180 / 0.013 | |
| E3 (the large lakes' area falls) | | 3 / 5 (median 735 → 439 km²) | 3 / 5 (554 km²) |

- **W₁₀**: the global index does not move. The active coasts' fall on seeds 23 (0.410 → 0.230) and 29 (0.451 → 0.364): the
  deposit fills their seaward foot.
- **E3 is not the load's alone.** Switching A1 of F183 off leaves the count at 3 / 5, but not on the same seeds:
  - with the load, seed 19's large lakes grow, 654 → 1 941 km² (without the load, 0 km²);
  - without the load, seed 23's grow, 1 033 → 1 340 km² (with the load, 288 km²);
  - seed 8's grow either way (906 → 1 858 / 2 000 km²).
  - Neither K_eff's slow sills nor the subsidence under the load explains the failures alone.
- **The old age** (300 Myr; the circumnavigation against the same world without the deposition):
  - (v) holds 5 / 5 without the deposition and 5 / 5 with it;
  - the sea 47.6 → 53.2 %;
  - circumnavigable 4 / 4, **none lost**;
  - the plains 66.7 → 61.6 %, the mountain 15.6 → 19.6 %;
  - the median altitude 53 → 72 m;
  - the large lakes 1 → 0 (medians).
  - The sea grows with the deposition, at both ages (the load's subsidence drowns the coasts).
- A_100 reproduces t1d bit for bit (the same median altitudes on the 5 seeds): the prolongation in force is t1d's.

### Step 5 — the chain to 2 048², seed 23, with step 3

- **The stages**:
  - C1 1.1 s; the physics level with A 33.9 s;
  - the amplifications 5.0 / 7.9 / 37.3 s;
  - the burial 0.09 s;
  - with Part 0's selection and extrapolation, **≈ 883 s against 300 s** (the chain's own overrun, F187).
- **At 2 048², against the same world without the deposition**: plains 48.2 → 50.3 %, hills 32.6 → 25.7 %, mountain 19.2 →
  23.9 %, median 73 → 70 m.
- **The plates** (for the author's eye):
  - `s5_map2048_off_A_seed23.png`: the whole map;
  - `s5_crops_seed23.png`: lake, delta, the delta's bicubic surface and its basement chain, plateau; each without | with.
  - What they show:
    - the stairs are gone;
    - the deep sea keeps the faint polygonal lines of step 2;
    - the delta's stripes remain (step 3's diagnosis);
    - the filled lake reads as a smooth plain.

### For the author (the continuation)

1. **The shelves.** Neither prolongation of the blocks forms the per-cell shelves (2.2 km on 5 / 5). The blocks hold the budget but
   lose the shelves, the deep-sea notches and the delta's grain. The per-cell marine step forms the shelves but takes 88 s.
2. **The losses' rule** of the bilinear prolongation is mine (t1d); t1e was measured worse.
3. **Step 3's reading**: the fill only on the stencils that hold a cover.
4. **E3** has no single cause between K_eff and the load.

## The decision, the recycling, the diagnosis and the chain per cell

**The decision on the marine step (to be confirmed by the author before it is applied)**, verbatim:

> Décision (à confirmer par l'auteur avant d'appliquer) : le marin en cellules est la référence physique. Les blocs ne
> sont pas adoptés : ils perdent le plateau, créent les stries du delta et multiplient les lignes en mer profonde. Garde
> leur code derrière son drapeau, avec la mesure. Le budget de 60 s pour la couche de physique (A5) était une
> répartition du relecteur : il est retiré. La seule règle est celle de l'auteur, 300 s pour la chaîne entière. Le
> surcoût du marin en cellules (~45 s) sera traité au F187 avec le reste du budget. Inscris-le à l'ADR.
>
> Règle générale, à ajouter à l'ADR : une grandeur grossière ne se lit jamais par sa cellule mère, elle s'interpole
> (cf. les marches de l'enfouissement et les lignes du recyclage).

- What it means here:
  - the per-cell marine step is the physical reference;
  - the blocks are not adopted, and their code stays behind `Deposition::marine_blocks` (off) with its measures;
  - **the 60 s budget for the physics level (A5) is withdrawn**: it was the reviewer's split;
  - the only budget rule is the author's, 300 s for the whole chain; the per-cell marine step's extra cost (~45 s) goes to F187
    with the rest of the budget.
- **The general rule (the author's)**: a coarse quantity is never read through its parent cell; it is interpolated.
  - Its two cases so far: the burial's 8 × 8 stairs (B2 read the parent's H) and the recycling's 1-cell lines at the plate edges
    (A3 read the 64² parent's consumed fraction).

**Before the commit, the author's three steps** (on mains):
1. **The recycling interpolated.**
   - F183 A3's consumed fraction F_j is now interpolated bilinearly from the 64² grid (the cell centres, the torus) in
     `history.rs::recycle_fn`, instead of read by the 64² parent cell.
   - The deep-sea lines' index (seed 23, the per-cell marine step): **201 → 161 cells**.
   - The controls hold (`the_accommodation_controls`): no consumption recycles nothing, bit for bit; with the load, the
     compaction and the recycling on, both mass controls hold at 10⁻⁹.
2. **The diagnosis, no correction** (seed 23, the chain to 2 048², step 5's blocks configuration; `results/diag_white_seed23.txt`,
   `results/diag_bumps_edge_seed23.txt`). The new white bumps are fine cells ≥ 1 200 m with the deposition and < 600 m without;
   there are three families:
   - **The burial on a basement the chain has lifted** (the lake crop's oval spot, the elongated islet top right, the spot at the
     bottom):
     - cells 4980 (map x 231–235, y 35–43), 17360 (x 413–421, y 127–141) and 56482 (x 318–332, y 433–445);
     - interior basins 1.7–2.9 km below the sea at 256², filled to ~sea level by the marine sediment (H 1.8–3.0 km);
     - the chain lifts their basement to the land floor, 0.5 m (`amplify.rs:903`, `lift_interior`), then the burial fills a
       volume H over that lifted basement (`amplify.rs`, `Transmit::bury` / `fill_level`): peaks of 1 761, 2 516 and 3 008 m
       where the basement chain gives 1 m.
   - **The burial on a sea floor the chain has raised** (the open basins near the coasts):
     - cells 12331, 4148, 7216, 9759: basement 1.6–3.5 km below the sea at 256², not interior;
     - the chain's fine basement there is only −0.4 to −1.6 km (its working field puts the ocean at 0 m, `amplify.rs:889`), so
       the fill to H overshoots: 1.4–2.0 km.
   - **The rock raised inland, no cover** (the points at the bottom right):
     - cells 39851, 39603, 41684, 42442, 39102 and others: H = 0, the rock 0.6–1.1 km higher than without the deposition;
     - **neither K_eff nor the load**: the rock rises more without K_eff (+1.5 to +2.1 km) and about as much without the load
       (+0.7 to +1.2 km);
     - **hypothesis, not measured**: the base level moved seaward by the filled basins (the rivers lengthened).
   - **The delta crop's vertical edge** (« avec ») is the basement's own straight coast, identical in the basement chain and in
     the world without the deposition, against which the bay's flat fill now stands. The burial's added height jumps +149 m at
     fine x 247, which is not on a 256² cell edge, so it is not the burial's stencil.
3. **The chain to 2 048² with the per-cell marine step** (seed 23; `results/c3_seed23.txt`):
   - **the stages**: C1 1.2 s; the physics level with A 84.0 s; the amplifications 5.1 / 8.5 / 39.9 s; the burial 0.09 s;
     **≈ 986 s to 8 192²** with Part 0's selection and extrapolation, against 300 s (F187);
   - **at 2 048²**: plains 48.2 → 49.8 %, hills 32.6 → 26.9 %, mountain 19.2 → 23.3 %, median 73 → 70 m;
   - **the cover's volume**, fine against coarse: −0.2 % on every cell, −7.0 % on land;
   - **the delta's 2-cell oscillation along x / y** (m): 18.9 / 15.9 per cell, 36.3 / 14.2 with the blocks;
   - **the plates for the author**: `images/c3_map2048_off_A_seed23.png` (the whole map) and
     `images/c3_crops_lake_delta_plateau_off_A_seed23.png` (lake, delta, plateau; without | with);
   - **what they show**: the shelves along the coasts, no deep-sea lines, no delta stripes; the burial's bumps of step 2 and the
     delta's straight edge remain (diagnosed, not corrected).
- **The adoption** waits for the author's eye on these plates.

- **The order of the measures**: the continuation's t1d, step 4 and step 5 ran before the recycling was interpolated; the lines
  index after it, the diagnosis and the per-cell chain ran after it.
- **The checks at commit** (`checks_f186_commit.txt`): the guard 6 / 6 Match (field and lakes); production's C1 hash unchanged;
  672 lib and 33 viz tests pass; `cargo check --workspace --tests` clean.

**At commit** (the author's words, verbatim): « Puis commit du F186 (drapeaux coupés par défaut, contrôles, rapport). L'adoption
viendra après le regard de l'auteur sur les planches de l'étape 3. »
- Committed: K_eff (`StreamPowerConfig::keff_g`, `Deposition::fluvial_keff`), the live marine step on blocks
  (`Deposition::marine_blocks`), the transmission and the burial by the level (`Deposition::transmit`, `amplify::Transmit`,
  `fill_level`), the recycling interpolated, and F185's `fluvial_dt_max`. Every flag is off by default; production's hash is
  unchanged; the guard 6 / 6 Match after the re-indexing.
