# Finding 159 — the multi-scale cascade prototype, 64² → 128² → 256² → 512² on the témoin

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ».

**Production is unchanged**:
- the guard 6 / 6 (field and lakes) was run before and after (`checks_before.txt`, `checks_after.txt`);
- the only production-file change is a pure extraction, pinned by a hash measured on the code before it.

**The stop rule R fires at the first level, 128².** The cascade holds its trunks and is cheap. But each level raises
the land and sharpens the crests, and the talus makes its own planar band at 33°.

## Files

| file | what |
|---|---|
| `f159_predictions.md` | my predictions, written before any run (the non-blind reading declared) |
| `f159_declared.md` | D1–D9, R and I, declared before any run, with two corrections made before the first run (D3 the sea, D6 the restriction) |
| `crates/ymir-core/src/cascade.rs` | the cascade: `Cascade`, `CascadeConfig::declared`, `LevelRecord` (the additive decomposition), `upsample2` (periodic Catmull-Rom ×2), `uplift_field`, `roll`, `shade_rgba` / `diff_rgba` (the viz's Ombrage formula) |
| `crates/ymir-core/src/erosion/stream_power.rs` | `talus_sweep` and `linear_diffusion` extracted from `incise_with_floor` unchanged, plus the test `the_f159_extraction_changes_no_output` (the hash pinned at 38a353e, with two negative controls: talus off and diffusion off each change it) |
| `crates/ymir-viz/src/ui/cascade.rs` | the « Cascade » window, its worker thread, and the test `the_cascade_worker_runs_level_by_level_then_all` |
| `crates/ymir-core/tests/f159_cascade.rs` | the bench: every D8 instrument, R, the ESTIMATE, the images |
| `f159_bench_output.txt` (run 2), `f159_bench_output_run1.txt` | the bench's printout. The two runs' measurements are **identical** line for line; run 2 differs only by the image amendments below |
| `images/` | the three crops (A range, B plateau, C basin + lake 4), the whole map per level, the level's difference map |

## B — what was built

### The core (`ymir_core::cascade`)

- **The coarse field**: the témoin's 64² C1 altitude, `c1_coarse_normalized_altitude(state, c1_default, ss, None)`,
  rolled (6, 37).
- **Per level**:
  - **agrandissement**: bicubic ×2;
  - **then steps until equilibrium or the cap**, each step:
    1. **soulèvement**: uplift;
    2. **érosion**: `incise_with_floor`, one iteration, relief-v3 at the level's cell, with the talus and the diffusion
       off inside it;
    3. **diffusion**: the extracted `talus_sweep` then `linear_diffusion`, at the declared weight.
- **The sea** is held at sea level during a level and restored at its end.
- **The record**: each level keeps the four sub-steps as an additive decomposition (upscaled + Σ uplift + Σ erosion +
  Σ diffusion = the result), with the time of each.
- **The unit test** (`the_cascade_reaches_equilibrium_and_its_substeps_add_up`), on a synthetic island, checks:
  - equilibrium at each level;
  - the decomposition adds up;
  - no land under the sea;
  - the bathymetry is restored;
  - the calibration within 10 %;
  - the input is untouched.
  - **Negative control**: no uplift, so the mean only falls.

### The viz: the « Cascade » button in the top bar

- **The window** opens on the workspace's world (its seed and framing roll). It never reads or writes the HD world.
- **Its own thread** runs the 64² tectonics and then the cascade:
  - « Niveau suivant » computes one level; « Tout » computes the rest; « Annuler » cancels.
- **The frieze**: 64² tectonique → 128² → 256² → 512², with steps, time and equilibrium (✓, or ⚠ plafond).
- **Per level**: the four sub-steps, with their times.
- **The views**: Ombrage (the workspace's formula), Hypsométrie, Différence; « Contour des dépressions ».
  - The cascade has no lake inventory, so the outlined depressions (filled more than 0.5 m by the priority flood)
    stand in for « Contour des lacs ».
- **A hover readout** in km and m.
- **The test** drives the worker: « Niveau suivant » gives one level, « Tout » the rest, and the coarse is the roll of
  the framed C1.

### Amendments after run 1 (NOT blind, images and viz only; no measurement depends on them)

1. **Crop B**: the declared rule (≥ 90 % land, mean > land p75) finds no window, and run 1 fell back to (0, 0) in
   silence, which is open sea.
   - The amended rule: ≥ 75 % land, mean > land p50.
   - The bench now fails rather than fall back.
2. **The sub-step views**: the cumulated fields are unreadable, because the uplift alone stacks about 30 km over a
   level that the erosion takes back.
   - The images get a 4th row: the level's total change, then Σ uplift, Σ erosion, Σ diffusion (m, saturating at the
     crop's p98, printed).
   - The viz's « Différence » view shows the same: the level's total change on « Agrandissement », then each
     process's own Σ.

## M — the measurements (cascade · témoin block mean at the same level)

Calibration (128²): the trial U₀ = 1e-4 gives a mean land altitude of 71 m against h_iso's 898 m, so **U₀ = 1.27e-3
m/yr per km of h_iso**. 63 steps to equilibrium; the talus inert (0.00 %).

| instrument | 128² | 256² | 512² |
|---|---|---|---|
| steps (all at equilibrium) / time | 36 / 0.1 s | 64 / 0.8 s | 119 / 6.6 s |
| talus active (% of land, last step) | 0.45 | 10.0 | **37.8** |
| planar walls 28° ± 0.5° (%) | **0.150** · 0.000 | **1.036** · 0.255 | **1.541** · 0.796 |
| planar band 33° ± 0.5° (%) — the talus angle | 0.000 · 0.000 | **1.456** · 0.193 | **4.354** · 0.292 |
| crest cells (% land) | 42.1 · 48.9 | 34.8 · 42.7 | 30.3 · 33.9 |
| crest −∇²z p50 / p90 (m/km²) | 39 / 160 · 32 / 136 | 96 / 432 · 74 / 345 | **293 / 1 186** · 122 / 697 |
| crest facets (% crests) | **8.63** · 4.33 | **16.20** · 4.90 | **23.88** · 9.26 |
| facet control: pyramids / paraboloids | 87.5 / 2.3 PASSES | 87.3 / 2.4 PASSES | 87.5 / 2.4 PASSES |
| coherent windows C > 0.7 (%) — the comb stand-in | 0.00 · 0.00 | **22.45** · 4.76 | **27.19** · 22.67 |
| terrain R8 (windows) [noise floor] | 0.139 (13) · 0.201 (20) [0.28 / 0.22, below the floor] | 0.041 · 0.028 [0.10 / 0.09, below the floor] | **0.084** · 0.070 [0.044 / 0.041] |
| valley spacing λ (cells) / Perron band (cells; ×π/2 for transects) | 10.3 · 5.9 / 2.6–5.1 (4.0–8.0) | 12.1 · 7.5 / 2.6–5.1 | 12.6 · 10.7 / 2.6–5.1 |
| relief p50, 12.5 km blocks (m) | **1 626** · 927 | **2 830** · 1 178 | **4 163** · 1 280 |
| mean land / p90 (m) | 893 / 2 300 · 631 / 1 406 | 1 114 / 2 930 · 649 / 1 450 | **1 488 / 3 938** · 659 / 1 469 |
| trunks ≥ 100 km², cascade → témoin, same level: p50 / p90 / max (cells) | 0.00 / 1.00 / 3.16 | 1.00 / 2.83 / 5.83 | 2.00 / 5.66 / 12.73 |
| trunks n+1 → n (cells of n): p50 / p90 / max | — | 0.50 / **1.41** / 6.73 (témoin's own: 0.50 / 1.12 / 5.32) | 0.50 / **0.71** / 6.40 (témoin's own: 0.50 / 1.00 / 5.22) |
| D6: level n+1 restricted − level n, on land: bias / mean \|Δ\| | — | **+227 / 249 m** | **+379 / 384 m** |
| relief p50 change from the previous level | — | **+74 %** | **+47 %** |

**For information, D8.3 at 512² (block means)**: crest facets are 19.25 % on S1, 17.62 % on the construction, 9.26 %
on the témoin and 23.88 % on the cascade. Walls at 28° are 0.160 % on S1 and 2.266 % on the construction.

**Cost**:
- 128² → 512²: **7.8 s**, calibration included;
- the time per step at 512²: 0.055–0.057 s (erosion 4.3 s, diffusion 2.3 s);
- the memory retained per level: 0.5 / 2.2 / 8.7 MB (estimated from the buffers);
- **the ESTIMATE to 8192²: 3 164–3 245 s = 53–54 min**. That is above the 15 min target and under the 1 h stop. It
  assumes 119 steps per finer level and an n log n scaling (D8.9).

## R — the stop rule

**It fires at 128²**:
- **walls 28°**: 0.150 % (4 cells of 2 670) against 0.000 %;
- **crest facets**: 8.63 % against 4.33 %, with the facet control passing.

**It would fire again at 256² and at 512² on every wall and facet criterion.** The trunk criterion never fires (p90
1.41 and 0.71 cells). The cost criterion does not fire (53–54 min < 1 h).

**Read plainly**:
- The 128² wall count is 4 cells and could be dismissed as noise.
- **The facet share cannot**: it doubles with each level (8.6 → 16.2 → 23.9 %), always 2–3× the témoin.

## What the measurements say (the reading; the causes are hypotheses to check)

1. **The trunks hold.** Between levels the trunks ≥ 100 km² move by a p90 of 1.41 cells (128 → 256) and 0.71 cells
   (256 → 512).
   - That is the témoin's own 1.12 and 1.00 between its block means.
   - Bicubic plus a re-equilibrated level keeps the big rivers where they were: **this part of the LOD idea works**.
   - Against the témoin's trunks at the same level they drift more, p90 5.7 cells at 512². The cascade is not the
     témoin's network; it has no coast warp and no FBM.
2. **The relief is not resolution-independent; my D scaling did not make it so.**
   - The calibration holds at 128² (893 m against 898 m).
   - Then every level raises the land: the restricted bias is +227 m, then +379 m.
   - The mean land altitude goes 893 → 1 114 → 1 488 m; the relief p50 goes +74 %, then +47 %.
   - At 512² the relief p50 is 4.2 km, against the témoin's 1.3 km.
   - The difference maps show it: the interfluves rise, the valleys stay.
   - **Hypothesis to check**: every land cell is a channel cell at these levels (A_c = 0.1 km² < one cell). The ridge
     cell's steady slope U / (K √A_cell) then doubles at each halving of the cell. And the diffusion at L_c = 0.4 cell
     (weight 0.032) is too weak to hold the new cell-scale ridges. This is Kwang & Parker's grid dependence, not
     removed by tying D to the cell.
3. **The talus becomes the shape.**
   - The talus moves 0.45 % of the land at 128², 10 % at 256² and **38 % at 512²**.
   - The 33° band reaches 4.35 % of the land (témoin 0.29 %).
   - The crests sharpen: −∇²z p90 is 1 186 against 697 m/km².
   - The talus's straight repose slopes are planar facets by construction. **The cascade has rebuilt the « crêtes en
     facettes » by another route**: their share is 23.9 % against the témoin's 9.3 %.
4. **The comb appears at 256² and 512²**:
   - coherent windows: 22.5 % against 4.8 %, then 27.2 % against 22.7 %;
   - regular straight gullies down the range's flanks are visible on crop A.
   - The terrain R8 is above the témoin's at 512² (0.084 against 0.070, both above the 0.044 floor). At 128² and 256²
     it is below its noise floor (too few windows).
5. **The valley spacing is not Perron's band.**
   - The transects read 10–12.6 cells against 2.6–5.1 (4.0–8.0 with the π/2 transect factor) at every level.
   - So the declared L_c = 0.4 cell does not set the spacing seen. The fluvial term, not D, sets it here (hypothesis to
     check).
6. **The cost is small and the equilibrium is fast.** 36 / 64 / 119 steps, all at equilibrium, 7.8 s in all.
   - The extrapolation (53–54 min) is dominated by the assumption that 4 096² and 8 192² need 119 steps each.
   - It is a weak estimate, not a measurement.
7. **Seen on the images, not measured**:
   - the coarse field's below-sea interior cells become inland seas held at 0 m. The témoin's block means show the
     same patches, but as basins, not open sea.
   - the cascade's coast is the bicubic of the 64² coast, visibly blocky at 128², with no warp. Both are declared
     consequences of D3 / D5.

## The predictions

**The reviewer's** (hypotheses, now checked):

| prediction | verdict |
|---|---|
| at 512², fewer than half the témoin's planar walls | **refuted**: 1.54 % against 0.80 %, about 2× |
| no crest facets beyond the témoin | **refuted**: 23.9 % against 9.3 % |
| trunks > 100 km² move < 1 cell (p90) between 256² and 512² | **held**: 0.71 |
| without scaling, the relief changes > 20 % between levels | **not tested as stated** (D was scaled). **With** the scaling it changes +74 % and +47 %, so the hard point is confirmed in a stronger form |
| 128² → 512² < 60 s | **held**: 7.8 s |
| the extrapolation to 8192² < 15 min | **refuted**: 53–54 min (ESTIMATE) |
| meta: at least one prediction false | **held** |

**Mine** (`f159_predictions.md`):

| prediction | verdict |
|---|---|
| P-M1 walls ≤ the témoin at every level | **refuted** at all three |
| P-M2 facets below the témoin at 512² | **refuted** |
| P-M2 the control holds | **held** (87.5 / 2.4) |
| P-M2 the construction above S1 at 512² | **refuted** (17.6 against 19.3) |
| P-M3 trunks p90 ≤ 1 cell, 128 → 256 | **refuted** (1.41) |
| P-M3 trunks p90 ≤ 1 cell, 256 → 512 | **held** (0.71) |
| P-M3 max > 3 cells | **held** |
| P-M4 λ 2–6 cells near Perron | **refuted** (10–12.6) |
| P-M5 R8 above the témoin at 512² | **held** (0.084 against 0.070) |
| P-M6 relief ±5–20 % per level, within ±25 % of the témoin | **refuted** (+74 / +47 %; 1.75–3.3× the témoin) |
| P-M7 128² in 100–400 steps | **refuted** (36) |
| P-M7 256² / 512² often at the cap | **refuted** (all at equilibrium) |
| P-M8 < 60 s | **held** |
| P-M8 the ESTIMATE 10–30 min | **refuted** (53–54 min) |
| P-R1 the rule does not fire before 512² | **refuted** (it fires at 128²) |
| P-B1 the guard 6 / 6 before and after | **held**: field and lakes = banc on all 6 states, before and after; core lib 621 → 624 tests, viz 32 → 33, all passing |
| meta | **held** |

## I — the questions for the author's look

**The images** (`images/`, north up, 100 km windows, each tile 512 px):

| image | window, 512² origin (= 8192²) |
|---|---|
| `f159_A_chaine.png` | (128, 320) = (2048, 5120) |
| `f159_B_plateau.png` | (176, 192) = (2816, 3072), the amended rule |
| `f159_C_bassin_lac4.png` | (74, 160) = (1184, 2560) |

**The layout of each image**:
- **row 1**: the cascade at 128², 256², 512², then the témoin at 8192² (block mean to 512 px);
- **row 2**: the témoin's block means at 128², 256², 512²;
- **row 3**: the 512² sub-steps, cumulated (agrandissement, + soulèvement, + érosion, + diffusion = the result);
- **row 4**: what each process did at 512²: the level's total change, Σ uplift, Σ erosion, Σ diffusion. The
  saturations are printed in `f159_bench_output.txt`.

`f159_level{128,256,512}.png` is the whole map per level; `_diff` is the level against its agrandissement (±300 m).

**Per level, the brief's four questions**:
1. Des crêtes en facettes, des murs, des promontoires, un « pinceau » ?
2. Des stries, des terrasses droites, des blocs ?
3. Les grands fleuves restent-ils à leur place d'un niveau à l'autre ?
4. Ce niveau te semble-t-il assez bon pour passer au suivant ?

**The same can be walked in the viz**: the « Cascade » button, then « Niveau suivant » or « Tout ». The window computes
the world the workspace frames, so it needs the témoin's seed and framing to match.

## What is open (for the author's decision, nothing chosen)

**The R verdict stands**: the cascade as declared stops at 128².

**What the measurements point at** (hypotheses, not decisions):
- **(a) the relief growth**:
  - a channel head that leaves a resolved hillslope at each level (A_c of a few cells, or L_c of 1–2 cells instead of
    0.4);
  - or keeping the coarse level as a constraint (D6 measured +227 / +379 m of drift; Schott 2024 retargets the coarse
    level).
- **(b) the talus's planar band**: Schott 2024 uses a noisy critical slope. Or a talus only above a level's resolved
  scale.
- **(c) the sea inside the continent**: the coarse field's below-sea interior cells held as open sea.
- **(d) the cost**: the extrapolation depends on the steps at the fine levels. Only a 1 024² level would measure it.

## The author's verdict at commit (2026-10-09), recorded as given
- « C'est pas mal. À tester à des résolutions plus hautes. »
- « Le soulèvement est bien trop fort. À 512², on arrive à des montagnes à +7 000 m sur un continent […] à une
  géométrie type Corse, [qui] ne devrait pas dépasser +3 000 m. »
- « Les structures sont intéressantes. Si on compare à la référence, on a quelque chose d'un peu moins organique. »
- « Pourquoi est-ce qu'on attend l'équilibre, vu que les articles, eux, n'attendent pas ? Pourquoi est-ce qu'on
  reconstruit tout le relief au lieu de l'affiner ? »
- « Les grands fleuves ne bougent que très peu. C'est bien. Par contre, je n'ai rien pour les visualiser sur le viz, ni
  leurs formes et affluents. Pareil pour les lacs. »
- F160 follows: one physics level, then Schott 2024's amplification transposed (no uplift), the final retargeting,
  rivers and lakes in the Cascade window, and a 1024² level.
