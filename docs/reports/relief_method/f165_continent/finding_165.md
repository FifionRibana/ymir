# Finding 165 — the continent's physiography: the physics level follows C1's history

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production is
unchanged**: the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`, `checks_after.txt`).
C1's outputs are bit-identical with the history recorder (a permanent test).

- Declared before any measurement of our worlds: `f165_declared.md`, with two amendments made on the references (the
  ETOPO product, the classifier's one adjustment on Europe).
- Predictions: `f165_predictions.md`.
- Bench: `f165_bench_output.txt`, `f165_europe_output*.txt`, `f165_a2_physics_output.txt`.

**The verdict: the « régime » stop fires.**
- **P1 is refuted on all four seeds.**
  - C1's h_iso above 1 000 m is narrow, but only 31–68 % of it lies near a convergent boundary or a suture. The rest
    is craton shield and arc.
- **The history creates no plateaus**: 0.0–0.1 % of the land, against Europe's 4.7 %.
- **So, as the brief rules, I stop and report, with no second mechanism.**

What the history does show is the finding of the round:
- **C1's land history is mostly SUBSIDENCE.** The continent starts as a uniform high plateau (h_iso(t₀)), and C1's own
  erosion closure thins it.
- **Driven by that history, the physics level loses 23–40 % of its land** to the sea, and 4–11 % of the remaining land
  is interior basins below sea level.
- **Its mountains rise to 55–66 % of the land** (Europe 24.8 %), and its plains fall to 1–4 % (Europe 51.6 %).
- **Hypothesis**: h_iso already contains C1's erosion, so the physics level erodes the same ground twice. The
  tectonic sources alone would be the cleaner drive. That is a second mechanism, for the author.

## Partie 0

1. **F164 committed** (`2d96819`), the author's seed and config edits left out.
   - Its ADR records the reviewer's decisions (a)–(d), the cascade frozen, and the author's six answers, verbatim.
2. **Checks**:
   - before: the guard 6 / 6 (field and lakes = banc on all six states), lib 634, viz 33;
   - after: `checks_after.txt`, lib 634 → 639, viz 33.
   - Five permanent tests are new:
     - the steady physics level, pinned by hash from before the shared loop;
     - the record leaves C1 bit-identical;
     - the history sums its uplift and the rebound lifts;
     - the steady history is the current level, bit for bit (negative control 1);
     - the classifier, the distance transform and the coast on their controls.

## A — the diagnostic: where do the « filins » come from?

### A1 — the inventory (read before writing the declaration; details in `f165_declared.md` A1)

- **C1's time loop**: 300 steps, Δt = 0.5 · dx / max|v| (non-dimensional).
  - The doc reads ~10–30 Ma (or ~200 Ma through the ocean ages). **No pinned physical duration**, hence T.
- **The transport**: Eulerian upwind, one velocity per plate. **With the rigid continental crust, continental cells
  have zero velocity** and a no-flux boundary.
  - So the Eulerian record at a continental cell IS its material history: **following the matter is clean for the
    land**.
  - Measured: 0.1–2.8 % of the final land cells changed type (oceanic → continental) during the run.
- **The continent-type parameter**: `Phase2InitParams`:
  - `num_plates` (8);
  - `cluster.seed_cluster_count` (1);
  - `continental_fraction` (0.29);
  - `craton_shield_fraction` (0.15);
  - `craton_thickness_ratio` (1.25).
- **The gorges' age selector** (`F121_AGE_K`) is « a calibration on a TARGET ». No link to C1's duration.
- **The current physics level**: U = U₀ · max(0, h_iso) / 1 000 m (`cascade.rs::uplift_field`), U₀ calibrated on the
  peak; equilibrium when mean |Δz| / mean U·dt < 0.01 over 5 steps.

### A2 — four seeds

The seeds are 10481999410520546993 (the témoin and the author's seed: the same number), 42, 1 and 9.

| seed | plates (init) | collision | subduction upper plate | craton | rift | sutures | land cells (type changed) |
|---|---|---|---|---|---|---|---|
| témoin | **2** (8) | **0** | 78 | 57 | 28 | 520 | 640 (2.8 %) |
| 42 | **2** (8) | **0** | 78 | 85 | 21 | 448 | 737 (0.1 %) |
| 1 | **2** (8) | **0** | 60 | 58 | 33 | 394 | 897 (0.1 %) |
| 9 | **2** (8) | **0** | 73 | 115 | 25 | 419 | 1 020 (2.3 %) |

- **No seed has a continent–continent collision at the end**, but all four have active margins (subduction upper
  plates), as the brief asks for at least one.
- **Every seed ends with 2 plates**: the accretion merges every pair.

**h_iso > 1 000 m (64²) against the physics level's land > 1 000 m (256², F164's regime)**:

| seed | h_iso width | h_iso near convergent / suture | physics width | ratio |
|---|---|---|---|---|
| témoin | 12.5 km (2 cells), 10 comp. | 50 % | 4.4 km | 2.8 |
| 42 | 12.5 km, 31 comp. | 68 % | 3.1 km | 4.0 |
| 1 | 12.5 km, 30 comp. | 31 % | 3.1 km | 4.0 |
| 9 | 12.5 km, 27 comp. | 41 % | 3.1 km | 4.0 |

- At 500 and 2 000 m the widths are the same: 12.5 km and 3.1–6.2 km. Near the boundaries: 31–52 % at 500 m, 31–100 %
  at 2 000 m.
- **A measurement floor, said plainly**: 2 × the median distance transform reads 2 cells as soon as most cells touch
  the edge. Both h_iso (2 × 6.25 km) and the physics level (2 × 1.5625 km) sit on their floor, so these widths are
  upper bounds.
- **The elongation is not usable as declared**: an area-weighted √(λ₁/λ₂) is dominated by straight one-cell
  components (values of 100–2 300). It is reported in the bench output, not read.

**Verdict A (the declared rule)**:
- **« both »**: h_iso is already narrow (≤ 12.5 km, its grid floor) but not boundary-bound (≥ 70 % fails on all
  four).
- The physics level narrows the land above 1 000 m by ×2.8–4.
- **P1 is refuted** (the near-share clause), on all four seeds.

## B — the references

**Europe**: 3 000 km, LAEA at 50° N 10° E, 1.5625 km cells, from ETOPO 2022 15″.

**The classifier**:
- **Mountain** is Kapos et al. 2000.
- **Hill** is R7 ≥ 150 m, below 500 m.
- **Plateau** is non-mountain at ≥ 500 m.
- The rest is **plain**.
- **The one adjustment on Europe, before any world of ours**: « hills are low », i.e. the hill test applies only below
  the plateau altitude. As first declared, Europe read plateau 1.6 % and the Meseta hill.

| | plain | plateau | hill | mountain | mountain components | width | width / √land |
|---|---|---|---|---|---|---|---|
| **Europe** (5.26 M km²) | **51.6 %** | **4.7 %** | 18.9 % | **24.8 %** | 894 | 15.6 km | **0.0068** |
| the Alps (400 km) | 19.9 % | 7.6 % | 10.4 % | 62.2 % | 11 | 31.2 km | 0.078 |
| Corsica | 1.0 % | 0.0 % | 35.3 % | 63.6 % | 2 | 9.4 km | 0.10 |

- Europe's land deciles: 40 / 81 / 117 / 153 / 195 / 262 / 383 / 580 / 889 m.
- Kapos's mountain share, 24.8 %, is the ~25 % of the literature.
- **The qualitative control**:
  - **holds**: the Paris Basin plain (96 %), the Meseta plateau (61 %), Corsica mountain (61–64 %);
  - **does not hold**: the Massif Central reads mountain (74 %), which is Kapos's definition. It is not adjusted.
- **Meybeck 2001's global figures could not be read** (BioOne serves a page, not the paper). Nothing is quoted.

**Corsica's coast** (B4, 195 m): dimension 1.12, length ratio 1.37, aligned 27.2 %, block index 0.78. At 49 m: 1.07,
1.10, 27.9 %, 0.92.

## C — the mechanism

- **The record**: 31 snapshots (every 10 steps), read-only. C1 is bit-identical (`the_record_leaves_c1_bit_identical`).
- **The historical level**:
  - U = Δh_iso / Δt per interval, the sea not held;
  - the stream power, the talus and the diffusion as F164, MFD p = 6;
  - **the flexural rebound**: (ρc/ρm) · (G ∗ E), Ĝ = 1 / (1 + (kα)⁴/4), α = 64.4 km (Te 25 km), by FFT, never in blocks;
  - 300 steps over T, no equilibrium stop, no deposition.
- **T calibrated once on the témoin: T = 1.155 × 10⁶ yr** (9 trials):
  - 10⁵ yr → 3 890 m; 10⁶ → 3 039 m; 10⁷ → 50 m; then bisection to 2 852 m.
  - **1.2 Myr is short** for C1's nominal 10–30 Ma: with U = Δh_iso / Δt, the uplift and the subsidence are fixed in
    metres, so the erosion's time decides the peak. A longer T erodes the land away (50 m at 10 Myr).
- **The controls**: all as declared. Control 3 is the late 10 %. The steady negative control is a permanent test (bit
  for bit).

## D — the measures

### Macro at 256² (B3)

| seed | run | plain | plateau | hill | mountain | width / √land | land km² | peak |
|---|---|---|---|---|---|---|---|---|
| | **Europe** | 51.6 | 4.7 | 18.9 | 24.8 | 0.0068 | | |
| témoin | F164 | 8.7 | 0.0 | 48.5 | 42.8 | 0.042 | 28 306 | 2 863 |
| | **history** | 1.2 | 0.1 | 39.9 | 58.8 | 0.062 | 20 376 | 2 852 |
| | no rebound | 3.6 | 0.2 | 31.3 | 64.9 | 0.059 | 13 982 | 2 594 |
| | late 10 % | 0.2 | 0.1 | 26.2 | 73.5 | 0.076 | 22 041 | 3 765 |
| 42 | F164 | 30.6 | 0.0 | 57.7 | 11.6 | 0.023 | 36 101 | 1 815 |
| | **history** | 4.1 | 0.0 | 40.9 | 55.0 | 0.042 | 21 748 | 2 455 |
| 1 | F164 | 36.3 | 0.0 | 50.5 | 13.2 | 0.022 | 40 833 | 2 498 |
| | **history** | 2.5 | 0.0 | 36.8 | 60.6 | 0.041 | 28 450 | 3 474 |
| 9 | F164 | 38.9 | 0.0 | 44.0 | 17.1 | 0.029 | 45 879 | 2 545 |
| | **history** | 1.7 | 0.1 | 32.5 | 65.7 | 0.050 | 35 371 | 3 964 |

**What the table says**:
- **Neither F164 nor the history makes plateaus** (≤ 0.3 % in every run, the controls included).
- **The history multiplies the mountain share by ×1.4–5** and drowns the plains (1–4 %).
- **It loses 23–40 % of the land**, and leaves 4–11 % of the rest as interior land at or below 0 m (subsided basins
  with no deposition, as the brief warned).
- **α × 0.5 / × 2** move the land by ±20 % and the peak by ±200–300 m, without changing the classes' picture.
- **The width / √land ratio is 0.02–0.08 against Europe's 0.0068, ×3–11 everywhere, F164 included.**
  - Hypothesis: it is structural. Our ranges keep real-metric widths (Corsica-like slopes), while the layout is
    compressed ×7.5. A Europe-like ratio would need ranges narrower than Corsica's in real km.
  - The ×2 test cannot pass at this scale whatever the drive.

**Snapshots and Σ** (`f165_snapshots.png`, `f165_sums.png`, the témoin):
- h_iso(t₀) is one uniform plateau, the initial continental crust.
- **Σ uplift is negative over almost all the interior**: C1's erosion and equilibrium sinks thin the crust. It is
  positive only on the active eastern margin.
- Σ rebound is a smooth ~α-wide dome.

### P2, P4, P6 (256²)

| seed | early belts: peak F164 → history | young belts | the peak cell | cratons' p50, no rebound against history | drainage density, late 10 % / history |
|---|---|---|---|---|---|
| témoin | 2 863 → 2 852 m (−0 %) | 0 cells | early | −29 % | 1.56× |
| 42 | 1 063 → 2 455 m (+131 %) | 0 | early | +1 % | 1.32× |
| 1 | 2 498 → 3 474 m (+39 %) | 0 | early | −17 % | 1.67× |
| 9 | 2 545 → 3 297 m (+30 %) | 649 | young | −14 % | 1.46× |

- **C1's thickening happens early**: three seeds have no « young » belt at all.
- **P6 is the reverse of the prediction**: the late uplift makes MORE drainage (a young, steep relief), not less.

### Texture at 2 048² (the frozen cascade), against Corsica at 195 m

| seed | F164, whole | F164, mountain class | history, whole | history, mountain class |
|---|---|---|---|---|
| témoin | 17 / 19 (A_dir, F) | 15 / 19 | 16 / 19 (facets, walls, A_dir) | 16 / 19 (facets ×3.4, walls ×2.9, A_dir) |
| 42 | 12 / 19 | 13 / 19 | 15 / 19 | 16 / 19 (facets, walls, A_dir) |
| 1 | 12 / 19 | 12 / 19 | 16 / 19 | 14 / 19 |
| 9 | 12 / 19 | 13 / 19 | 16 / 19 | 16 / 19 (facets, walls, A_dir) |

- **The history fixes the low-relief seeds**:
  - F164's seeds 42, 1 and 9 failed the slopes and the 1.5–12.5 km octaves, their land being low (mean 125–158 m);
  - the history passes them.
- **But it fails the facets and the 28° walls everywhere** (×2–3.4 in the mountain class). **P5 is refuted.**
- **F by class** (history): plain 39–54 %, plateau 0–12 %, hill 7–8 %, **mountain 0.35–0.67 %** (Corsica's mountain
  class 0.6 %). The flat cells are in the plains, as F164's attribution found.

### The coast (B4, 2 048² against Corsica 195 m)

| seed | run | dimension | length ratio | aligned | block index | land lost / gained |
|---|---|---|---|---|---|---|
| Corsica | | 1.12 | 1.37 | 27.2 % | 0.78 | |
| témoin | F164 → history | 1.01 → 1.06 | 0.95 → 1.00 | 28.9 → 27.1 % | 2.8 → 1.2 | 22.4 % / 2.8 % |
| 42 | | 1.01 → 1.08 | 0.98 → 1.06 | 40.5 → 32.4 % | 6.3 → 3.0 | 36.1 / 1.1 |
| 1 | | 1.03 → 1.05 | 0.96 → 1.01 | 38.5 → 33.8 % | 6.0 → 3.6 | 21.9 / 0.8 |
| 9 | | 1.05 → 1.06 | 1.03 → 1.06 | 47.2 → 43.3 % | 8.2 → 6.8 | 19.3 / 5.1 |

- **No coast measure is degraded; all four move toward Corsica.**
  - The subsided margin leaves the 64² block lines: the block index falls by 17–56 %.
  - That is the side effect of losing a fifth to a third of the land.
- **The new depressions at the warp steps** (the clause as reformulated): **7–24 % of the count before**, at every
  level, in the F164 chains as in the history's.
  - F164's net counts (−13 to +19 %) hid them: the warp opens about a tenth of new depressions while closing as many.
  - Reported (not in this round's R).

## R — the stop rules

- **Success: 0 / 4 seeds.** No seed meets the mountain / plateau / plain fractions (×1.5) or the width ratio (×2), and
  P5 fails. The coast clause alone passes.
- **The « C1 » stop: no.** It needs P1 to hold, and P1 is refuted on all four seeds.
- **The « régime » stop: YES.** P1 is refuted, and the history's plateau fraction is under 50 % of Europe's (2.35 %)
  on 4 / 4 seeds. **Stopped and reported. No second mechanism.**

## What C1 lacks, read from the measures (for F166, nothing chosen)

- **Collisions**: none at the end on any seed. Every run ends with 2 plates after the accretion merges. The thickened
  zones are arcs at one active margin, plus craton shields.
- **A history of uplift**: C1's land history is dominated by its own thinning (erosion and equilibrium sinks) from a
  uniform initial plateau. **Driving the physics level by Δh_iso adds the physics erosion to C1's erosion**
  (hypothesis).
  - A drive by the **tectonic sources only**, Davis-Suppe + subduction + rifting (C1's own terms, read-only), is the
    candidate.
  - **It is a second mechanism: the author's to allow.**
- **Wide low-relief high ground**: neither drive makes plateaus. The initial plateau is the only one C1 ever has, and
  the history erodes and drowns it.
- **The scale**: the width / √land test sits ×3–11 above Europe in every run. That follows from the real-metric relief
  in a ×7.5-compressed layout, not from the drive (hypothesis).

## Images (north up; Copernicus notice in `images/NOTICE.md`; ETOPO 2022 public domain)

- `f165_classes_europe.png`, `f165_classes_alps.png`, `f165_classes_corse.png`: the classes of the references (sea
  blue, plain light green, plateau tan, hill olive, mountain brown).
- `f165_seeds_256.png`, `f165_seeds_2048.png`: rows = the four seeds, columns F164 | history, shaded, the same framing.
- `f165_classes_seeds.png`: the same layout, the classes.
- `f165_snapshots.png`: the témoin's h_iso at t₀, mid and final.
- `f165_sums.png`: the témoin's Σ uplift, Σ erosion and Σ rebound (red +, blue −).
- `f165_coast_25km.png`: 25 km coast crops, Corsica | F164 | history.

**Visible artefacts, reported**:
- the drowned margins (light-blue striated shelves) of the history;
- on seed 9, an uneroded rim along the map border, where the land touches it. Hypothesis: border cells route off the
  map and are never eroded.

## Queued (checked only, as the brief asks)

**The river width at export: the 0.5.0 format already carries one** (`export/rivers_ll.rs`):
- `bed_width_m` per vertex, the law a · Q^0.5 (`BED_WIDTH_B` = 0.5, `bed_width_a` in the header's `bed_width_law`);
- `valley_width_m` per vertex;
- `width_m` at each mouth.
- The calibration to 100–150 m at the largest mouth is the queued work; nothing was changed.

## Predictions

**The reviewer's**:
- **P1**: refuted (the near share 31–68 %; the width clause holds at the grid floor).
- **P2**: refuted. The early belts gain (+30 to +131 %) or stay (−0 %); three seeds have no young belt.
- **P3**: held, trivially (plateau 0.0 → 0.0–0.1 %, under half of Europe's).
- **P4**: refuted. The cratons' p50 drops by > 20 % on 1 of 4 seeds.
- **P5**: refuted (the facets and the walls fail in the mountain class).
- **P6**: refuted. The density rises by ×1.3–1.7 under the late uplift.
- **P7**:
  - the alignment part holds on 3 / 4 seeds (≥ 1.5 × Corsica's 27.2 %; the témoin 28.9 %);
  - « the history changes almost nothing » is refuted: the block index falls by 17–56 %, the alignment by 6–20 %.
- **Meta**: held.

**Mine**:
- **Held**:
  - A1's < 5 % changed type;
  - P1 refuted;
  - the active margins present;
  - the Massif Central reads mountain;
  - Meybeck unreadable;
  - P3 held;
  - P6 refuted;
  - P7's alignment;
  - the « régime » stop;
  - the subsided basins ≥ 1 % (4–11 %);
  - success not met;
  - meta.
- **Refuted**:
  - « at least one seed has collision cells » (none);
  - the Europe ranges (plateau 4.7 %, not 8–20 %; plain 51.6 %, above 30–50 %);
  - the mountain width / size 0.05–0.15 (0.0068);
  - T 1–30 Myr (1.16 Myr);
  - P2 « partly »;
  - P4 holds;
  - P5 « 1–3 elements » (2 fail, but they are facets and walls, not λ or an octave);
  - « the history moves the coast < 10 % »;
  - « P3 +2 to +8 points » (+0.0–0.1).

## Cost

- **Wall time is not usable**: the machine suspended during control 3 of seed 9 (33 773 s wall for a ~7 s run;
  [[bench-wallclock-suspend]]).
- From the other runs:
  - the historical physics level takes 4.9–11.5 s against the steady 3.3–12.1 s, so **+2 to +7 s per world** (the
    reference is 249.8 s);
  - a frozen chain to 2 048² takes 27–35 s;
  - the T calibration takes 9 trials, ~50 s.
- The history's recording costs 0.2 s per C1 run.

## At commit (2026-10-10)
- **(a) The drive by C1's tectonic sources only is allowed.** It is not a second mechanism but the same one, shared
  correctly: the tectonics supplies the forcing, and the erosion is done once, by the physics level. It is F166.
- **(b) C1's collisions come after F166.**
- **(c) The « mountain width / √(land area) » test is WITHDRAWN.** It was a design error of the reviewer's: at European
  scale a range would be ~1 km wide on our map, so the test was incompatible with relief in real km in a ×7.5-compressed
  layout.
- **(d) The warp's new depressions**: only how many survive to 2 048² is measured, nothing more.

**The reviewer's errors, recorded**:
- the width test (above);
- **P6**: a late uplift makes the relief younger and steeper, hence MORE drained, not less dissected.

**The author on F165, verbatim**: « Le littoral est intéressant. Trop de terres submergées mais c'est aussi
intéressant d'en avoir un peu. Je pense qu'il faut un équilibre entre les deux. […] Les deux de droite sont à mon avis
mon terrain avec les même structures prévisibles. Mais on n'est pas très loin de la 1ère image. Je pense qu'on peut
s'en rapprocher encore, mais on va gérer déjà ce qui est prioritaire. »

**The author's answers, verbatim**:
1. « Un sommet qui varie avec la tectonique. Un sommet plus haut peut aussi être intéressant. »
2. The sea-level rise: « Oui » (it will be F167).
