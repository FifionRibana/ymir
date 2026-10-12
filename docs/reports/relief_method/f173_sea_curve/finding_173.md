# Finding 173 — a real sea-level curve, the attribution of the drowning, and the origin of the coastal rampart

**Status: built, measured, gated. Production unchanged. Nothing committed without the feu vert.**

**Verdict**:
- **F172's excess came from the amplification levels.** They carry 83–102 % of the length ratio's excess when they erode at
  −134 m, while the physics level's lowstand alone leaves the coast as it was but drowns 2–14 % of the land.
- **The published curve** (Spratt & Lisiecki 2016, at the physics level only, the chain at 0 m) **removes the excess but brings
  no Corsica-like coast.**
  - The dimension stays 0.99–1.02 and the length ratio 0.90–0.99, as without a sea.
  - It still drowns 2.7–7.3 % of the land at the intermediate and old ages: closed lowland basins, not rias.
- **The success rule fails** (0 seeds). Neither the « too drowned » nor the « no effect » stop fires as written.
- **The rampart's origin**: 75–86 % of the coastal prism was added while the coast was a convergent boundary, 65–85 % against an
  oceanic plate.
- **The wall index W₁₀** tracks the 2 048² walls (Pearson 0.96, Spearman 0.97). It fires at young on **33 / 33 screened seeds**
  and at old on none: above F172's 80 % flag, so a probable C1 bias.

**Files**:
- Declared before any measurement: `f173_declared.md`. Predictions: `f173_predictions.md` (written first).
- The bench `crates/ymir-core/tests/f173_cascade.rs`:
  - `f173_measure` → `bench_measure.txt` (A, B4, C1, C2, D), tabulated in `bench_table.txt`;
  - `f173_screen` → `bench_screen.txt` (C3).
- The updated catalogue: `catalogue.tsv`.
- `images/`. Checks: `checks_before.txt`, `checks_after.txt`.

## Part 0

F172 committed (`f6607c6`), with the reviewer's three errors (a constant lowstand, a lowstand applied to levels without duration,
« rampart = convergent perimeter »), decisions (a)–(d) and the author's words, verbatim.

**The panel** is the seeds 8, 18, 19, 23 and 29. Their selections: 16 plates on 8, 23 and 29; 20 on 18 and 19.

## B — the curve (`history::SeaCurve`, `physics_history_sea`)

- **The series**: Spratt & Lisiecki 2016 (Clim. Past 12, 1079–1092), the column `SeaLev_longPC1`.
  - « Scaled first principal component of five sea level reconstructions (0-798 ka) »: 799 points every 1 kyr, in m above
    present.
  - In the paper: the abstract p. 1079, Fig. 2 p. 1085, Table 1 p. 1084.
  - The numbers are not in the PDF. They come from the archive the paper names (Data availability, p. 1090): NOAA's World Data
    Service for Paleoclimatology, study 19982, `spratt2016.txt`, downloaded 2026-10-12, kept in `data/sealevel/` (git-ignored),
    parsed, not digitised.
  - Taken as is: +8.96 m at 0 ka, min −130 m, mean −52.9 m.
- **Before 798 ka: 0 m.**
- **The physics level only**: the base level at each sub-step's mid-age. The steps overlapping the last 798 kyr are cut into
  ⌈d / 2 kyr⌉ sub-steps (43 steps × 10 in linear time), with the uplift spread uniformly. The chain erodes at 0 m.
- **The controls** (the old témoin):
  - s ≡ 0 through the constant path is F171 bit for bit (the existing permanent test).
  - The constant lowstand through the new entry point is F172 bit for bit (a new permanent test,
    `the_sea_curve_entry_point_keeps_f172_and_sub_steps`, which also checks the parsing, the interpolation and the sub-step
    count).
  - **The curve path with s ≡ 0 against F171**: max |Δz| 593 m locally, the classes move 0.03 points. The 2 kyr sub-steps
    re-route some channels; not bit for bit, as declared.
- **The viz**: « Mer : aucune / bas niveau constant / courbe quaternaire » (the curve read from `data/sealevel/spratt2016.txt`).
- `images/f173_sea_curve.png`: s(t) over 0–800 ka, with the sub-steps.

## A — the attribution of F172's drowning (the panel, young and old)

| seed | age | none (dimension / LR) | A-phys (lost) | A-ampl (lost) | A-both = F172 (lost) | the LR excess from the amplifications |
|---|---|---|---|---|---|---|
| 8 | young | 0.99 / 0.92 | 0.99 / 0.90 (2.3 %) | 1.21 / 1.87 (1.0 %) | 1.21 / 1.85 (2.7 %) | 102 % |
| 8 | old | 1.00 / 0.95 | 1.02 / 0.96 (9.5 %) | 1.37 / 2.75 (3.0 %) | 1.41 / 3.09 (10.7 %) | 84 % |
| 18 | young | 0.97 / 0.89 | 1.01 / 0.93 (12.8 %) | 1.26 / 2.15 (1.5 %) | 1.29 / 2.26 (13.5 %) | 92 % |
| 18 | old | 1.02 / 1.00 | 1.01 / 0.94 (14.3 %) | 1.39 / 3.02 (6.4 %) | 1.42 / 3.16 (15.6 %) | 94 % |
| 19 | young | 0.98 / 0.89 | 0.98 / 0.89 (4.3 %) | 1.22 / 1.86 (1.1 %) | 1.23 / 1.94 (5.5 %) | 93 % |
| 19 | old | 1.02 / 0.98 | 1.02 / 0.95 (13.6 %) | 1.39 / 2.91 (4.3 %) | 1.41 / 3.07 (15.1 %) | 92 % |
| 23 | young | 1.00 / 0.92 | 1.00 / 0.93 (4.7 %) | 1.23 / 1.92 (1.1 %) | 1.24 / 2.01 (5.3 %) | 92 % |
| 23 | old | 1.02 / 0.98 | 1.02 / 0.95 (11.3 %) | 1.37 / 2.79 (4.6 %) | 1.39 / 2.83 (12.4 %) | 97 % |
| 29 | young | 1.00 / 0.93 | 1.00 / 0.92 (4.0 %) | 1.25 / 2.03 (1.7 %) | 1.26 / 2.11 (4.4 %) | 94 % |
| 29 | old | 1.01 / 0.95 | 1.02 / 0.97 (13.6 %) | 1.36 / 2.69 (3.9 %) | 1.40 / 3.05 (12.8 %) | 83 % |

The dimension's excess from the amplifications: 89–99 %.

- **The coast's dendritic fringe is the chain's work.** Eroding at −134 m, the amplification levels cut the 0 to −134 m strip
  into 0.4–6 km valleys that the rise floods.
- **The land lost is mostly the physics level's**: 2–14 %, from its 2.58 Myr lowstand grading the lowlands below 0 m.

## D — the curve's effect (the panel × 4 ages)

| seed | none → curve: dimension / LR (bands) | the land lost | the classes' move | the longest drowned valley, Spearman(length, area) | the 0 to −134 m zone (median) |
|---|---|---|---|---|---|
| 8 | B0 0.99 / 0.91 → 1.00 / 0.94 (1); young → 0.99 / 0.91 (2); intermediate → 1.01 / 0.96 (1); old → 1.01 / 0.95 (2) | 0.7 / 1.0 / 2.7 / 3.7 % | 2.1 / 1.6 / 1.9 / 1.8 pt | 49 / 26 / 19 / 23 km; 0.5, 0.5 | 0.2–1.9 km |
| 18 | → 1.02 / 0.99 (0); 1.02 / 0.98 (1); 1.00 / 0.93 (1); 1.00 / 0.95 (2) | 1.5 / 2.3 / 6.1 / 7.3 % | 3.5 / 3.7 / 2.9 / 2.4 pt | 44 / 64 / 28 / 21 km; −0.09, 1.0 | 0.3–2.2 km |
| 19 | → 1.01 / 0.93 (1); 0.99 / 0.90 (2); 1.00 / 0.92 (2); 1.01 / 0.94 (1) | 1.6 / 2.2 / 4.4 / 6.3 % | 3.0 / 2.3 / 2.9 / 2.4 pt | 37 / 31 / 28 / 33 km; 0.95, 0.77 | 0.3–2.1 km |
| 23 | → 1.02 / 0.97 (1); 1.00 / 0.93 (1); 1.02 / 0.95 (2); 1.01 / 0.93 (0) | 1.3 / 1.5 / 3.0 / 5.5 % | 3.8 / 2.3 / 2.5 / 2.4 pt | 54 / 28 / 16 / 25 km; 0.87, 0.63, 1.0, 0.10 | 0.3–2.0 km |
| 29 | → 1.01 / 0.94 (1); 1.01 / 0.95 (1); 1.02 / 0.96 (2); 1.02 / 0.97 (2) | 2.7 / 1.8 / 3.0 / 4.9 % | 3.4 / 2.1 / 2.4 / 2.3 pt | 37 / 26 / 25 / 27 km; 0.10, 0.80, 1.0 | 0.4–2.2 km |

- **The coast moves by −0.02 to +0.05 in dimension and −0.03 to +0.10 in length ratio**, inside 0–2 of the 4 bands. The length ratio
  stays under its band (≥ 1.15) everywhere.
- **The aligned share** stays 0.29–0.34 and the block index 0.63–1.30.
- **The drowned valleys**: 1–7 per world, the longest 16–64 km.
  - Where ≥ 3 valleys exist, Spearman(length, the river's area) is 0.5–1.0 on 4 of 5 seeds, the large rivers' mouths. They are
    fewer and shorter than F172's.
- **The land lost** grows with the age: 0.7–2.7 % at B0 and young, 2.7–7.3 % at intermediate and old.
  - These are closed lowland basins graded below 0 m during the lowstands, then connected to the sea (`f173_coast_25km.png`,
    the deep-blue patches of the curve's column).
- **The texture** (curve): the slope guard holds (mountain-class p50 0.30–0.44); σ_θ at 391 m 23.6–25.6°; R2_g 0.66–0.81.
- **The lakes**: 5–19 per world, 439–6 133 km².
- **The sensitivities** (seed 8, old):
  - the stack's mean (−52.9 m) from 2.58 Ma to 798 ka: 1.02 / 0.96 (2 bands), as the main;
  - the chain at the mean: 1.31 / 2.42. Lowering the amplifications' base again over-cuts.
- **Cost**: the curve doubles the physics level (21–32 → 51–71 s CPU, ~390 more sub-steps) and adds 2–8 % per world (459–508 s
  against 425–483 s).
  - The C3 screening took 20 326 s of wall time for 33 seeds × 2 physics levels; the machine ran slow again.

## C — the rampart

### C1 — the coastal prism's history (today's land within 2 cells of the coast)

| seed | added at a C-O convergent boundary | C-C | no convergent boundary within 2 cells | by the plate count ≥ 10 / 7–9 / ≤ 6 |
|---|---|---|---|---|
| 8 | 73.3 % | 4.8 % | 21.9 % | 41 / 0 / 59 % |
| 18 | 65.2 % | 10.1 % | 24.7 % | 54 / 22 / 25 % |
| 19 | 84.6 % | 1.6 % | 13.8 % | 41 / 17 / 43 % |
| 23 | 79.2 % | 2.6 % | 18.2 % | 37 / 31 / 32 % |
| 29 | 81.1 % | 4.3 % | 14.6 % | 64 / 29 / 8 % |

- **75–86 % of the coastal prism was added at a convergent boundary**, mostly an ocean–continent one.
- Much of it came early (37–64 % while ≥ 10 plates remained).
- So today's « no-boundary » coast is a former active margin, its wedge kept by the hard basement. The reviewer's hypothesis
  holds.

### C2 — the wall index (declared: W₁₀ ≥ 0.30 and W₁₀ ≥ 2 W_in, slope ≥ 8° at 256²)

- **Against the 2 048² coastal walls** (> 28°, within 10 km), over the panel's 20 no-sea worlds: **Pearson 0.96, Spearman 0.97.**
- **W₁₀**: B0 0.76–0.86, young 0.43–0.63, intermediate 0.05–0.17, old 0.02–0.14.
- The rampart fires at B0 and young on 5 / 5 seeds, at intermediate and old on none.

### C3 — the screening (33 seeds; F172's configuration; `catalogue.tsv` gains `W10`, `W_in`, `rampart_W`)

- **rampart_W at young: 33 / 33 (100 %)**; at old: 0 / 33.
- **This is above F172's flag (> 80 % at young): a probable C1 bias, not a variety.** Every young continent has the wall ring.
  With C1: every coast was an active margin early on, and the wedge reaches the whole coast.

## R

| rule | result |
|---|---|
| the series not found | found (NOAA), used |
| the curve's success (≥ 4 seeds × ≥ 2 ages: ≥ 3 bands, < 5 %, < 3 points) | **0 seeds** (no world reaches 3 bands) |
| still too drowned (LR > 1.60 on ≥ 4 seeds) | no: LR ≤ 0.99 everywhere |
| no effect (Δdimension < 0.01 everywhere) | no: up to +0.05 on several worlds (the dimension stays ≈ 1.0) |

No stop fires as written; the curve does not enter the profile.

## The reading (hypotheses to check)

1. **The sea level is not what makes Corsica's coast.**
   - The realistic curve leaves the coast as smooth as without a sea.
   - F172's Corsica-like dimension was an artefact of the chain eroding a −134 m base it never should have had.
2. **What drowns with the curve is closed lowland basins**, cut below 0 m by the physics level during the lowstands (no
   deposition fills them). Deposition (Yuan et al. 2019) is the candidate for both the basins and any rias.
3. **The rampart is C1's early active margins**, a uniform ring at young on every seed. Hence C1 (b): passive margins, and the
   wedge's reach along margins that stopped converging.

## Predictions

### The reviewer's

| held | refuted |
|---|---|
| Z1 (83–102 %, 10 / 10 worlds) | Z2's rises (dimension +0.02 or more on few worlds; LR ≤ +0.10) |
| Z2's « neither exceeds the bands » (trivially) | Z3 at the intermediate and old ages (land lost 5.5–7.3 % on 3 seeds at old) |
| Z3 at young (lost 1.0–2.3 %; the classes < 3 points on 4 / 5) | Z4's « ≤ 50 km » on 2 / 20 worlds (64 and 54 km) |
| Z4 mostly (Spearman ≥ 0.4 on 4 seeds) | |
| Z5 (75–86 %, 5 / 5) | |
| Z6 (100 %) | |
| the meta | |

### Mine

| held | refuted |
|---|---|
| no curve success | Z1 refuted (it held: the amplifications, not the physics level) |
| the « too drowned » stop not firing | Z2's dimension +0.02–0.08 (≈ 0) |
| Z6 (50–80 %; it was 100 %) | Z3's 0.2–3 % (up to 7.3 %) |
| the per-world cost < +10 % | Z4's Spearman refuted (it held) |
| the meta | Z5 refuted (it held) |
| | the physics level's +30–60 % (+100–150 %) |

## Open, for the author (nothing chosen)

- **(a)** The sea in the profile: the curve at the physics level is the physical forcing, but it brings no coast and drowns lowland
  basins without deposition. Deposition next (Yuan et al. 2019), or the sea left out until then.
- **(b)** The rampart: C1 (b), the passive margins (a uniform ring on 33 / 33 seeds at young).
- **(c)** The coast's smoothness has no lever left in the sea level. Candidates: deposition, the structural grain, coastal
  processes (none built).
- **(d)** The wall index (Pearson 0.96) as the rampart's measure.

## At commit (2026-10-12)
- **(a)** The sea does **not** enter the default profile until deposition exists. The code and the curve are kept.
- **(b)** The rampart is a C1 bias (100 % at young). Coastal ranges on passive margins exist on Earth; the defect is that they are
  everywhere.
- **(c)** The sea level is not the coast's lever.
- **(d)** The wall index W₁₀ is adopted.
- **The reviewer's error, recorded**: Corsica taken as the reference of every coast. It stands for steep rocky coasts only; real
  low coasts are often fairly smooth (Landes, the North Sea), their irregularity coming from deposition and waves.
- **To correct (F174)**: Spratt & Lisiecki's series reads +8.96 m at 0 ka. It is re-zeroed to 0 m today: an alignment of
  reference, not a calibration.

**The author's view, verbatim**: « Je ne suis pas convaincu par l'aspect dendritique à la montée des eaux. Ca ne fait pas
"naturel". Ca manque toujours d'irrégularité. »

**The order, decided by the author** (« Ok go »): F174 = the measures, then the passive margin; F175 = deposition (Yuan et al. 2019).
