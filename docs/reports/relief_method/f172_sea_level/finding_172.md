# Finding 172 — the sea-level rise, and a screening of the seeds

**Status: built, measured, gated. Production unchanged. Nothing committed without the feu vert.**

**Verdict**:
- **The rise is a mechanism, and an excessive one.** On all 16 worlds the coast's dimension rises by 0.07–0.39 and its length
  ratio by 0.2–2.0.
  - After the rise they read 1.05–1.39 and 1.11–3.01; Corsica is 1.12 and 1.37.
  - Only one world lands inside the declared bands (the témoin, young). The success rule (« closer to Corsica ») fails on
    0 / 4 seeds.
  - The land lost is 0.7–18 %; the classes move by up to 9 points.
  - **Cause**: without deposition, the lowstand grades the coastal lowlands below 0 m. The rise turns them into a shallow sea
    fringed with dendritic interfluves.
- **The « coastal rampart » is a ring of basement walls, not a ring of mountains.** At B0:
  - the mountain class's walls (> 28°) are 19–26 % within 10 km of the coast, against 5–6 % beyond 20 km;
  - they are 12–18 % on the basement and shields, against ≈ 1 % on the cover.
  - The declared rampart index never fires (R ≤ 1.9; the interior is mountainous too).
- **The screening** (33 seeds, physics 256²):
  - every seed yields a qualified continent (16 plates on 25, 20 on 7, 24 on 1);
  - 15 % of the seeds (5, all at the old age) meet Europe's three bands;
  - the perimeter's convergent share does not predict the rampart index (Spearman −0.08 young, 0.17 old).

**Files**:
- Declared before any measurement: `f172_declared.md`. Predictions: `f172_predictions.md` (written first).
- The bench `crates/ymir-core/tests/f172_cascade.rs`:
  - `f172_rise` → `bench_rise.txt` (A, B, C);
  - `f172_screen` → `bench_screen.txt` (D).
- The catalogue: `catalogue.tsv` and `images/f172_catalogue.png`.
- `images/`. Checks: `checks_before.txt`, `checks_after.txt`.

## Part 0

- F171 committed (`3c0b653`). Its « At commit » block records:
  - (a) K2 as the continent profile's rock;
  - (b) the change of criterion (mechanisms by their effect, the macro layout by a distribution of seeds; A the author's slider,
    its default provisional; ρ_K 10);
  - (c) the coast judged against Corsica in absolute value;
  - (d) the walls, facets and rampart measured here.
- The author's view and position are recorded verbatim.

## The constants (read in the sources; `docs/refs/REFERENCES.md`)

- **Δ = 134 m**: Lambeck et al. 2014, PNAS 111(43), 15296–15303.
  - The abstract: « a slow fall to −134 m from 29 to 21 ka BP ».
  - The Discussion: « greater than the frequently cited −125 m ».
  - Read on PubMed Central (PMC4217469); the PNAS site refused the fetch, HTTP 403. **The article prevails** over the brief's
    « ~120–130 m ».
- **T_q = 2.58 Myr**: ICS, the GSSP list (the Gelasian, Monte San Nicola; the base of the Quaternary, ratified 2009).

## B — the mechanism

- **The physics level** (`HistoryRun::sea_low`, `SeaLow`, `Drive::sea`): during the history's last T_q, the base level is −Δ.
  - The step's remaining time is ≤ T_q; that is the last 137 of 300 steps in linear time.
- **The amplification levels**: the chain runs on the field raised by Δ (its sea at the real −Δ) and is lowered by Δ at the end.
- **At the end**: the sea at 0 m; the final ocean is the largest periodic sea component of the cells ≤ 0.
- **No deposition** (the declared limit).
- **The negative control** (permanent; `the_lowstand_at_zero_is_f171_bit_for_bit`): Δ = 0 is F171's K2 bit for bit, in linear
  time and at an age. Δ = 134 m changes the field.
- **The viz**: « Remontée de la mer » (the history level and the chain, its results shown lowered), and « vue bas niveau » (the
  field shown with the sea at −134 m).

## A — the diagnostics (F171's K2 worlds, no rise)

### A1 — the perimeter (64²)

| seed | convergent | divergent | transform | no boundary |
|---|---|---|---|---|
| témoin | 23.1 % | 4.7 % | 1.2 % | 71.0 % |
| 42 | 14.8 % | 4.1 % | 0.0 % | 81.0 % |
| 1 | 18.5 % | 12.4 % | 1.4 % | 67.8 % |
| 9 | 6.4 % | 5.5 % | 0.6 % | 87.5 % |

- Most of the coast is no plate boundary (68–88 %).
- **Yet Σ Δs_DS peaks at the coast and halves every 1–2 cells inland**:
  - témoin 1.43 → 0.92 → 0.71 → 0.36 …;
  - seed 1: 1.02 → 0.59 → 0.35 → 0.18.
- The Davis-Suppe wedge is the continental margin's own thickening, not only at the convergent stretches.

### A2 — the mountains near the coast; the rampart index (declared: f₂₀ ≥ 0.40 and R = f₂₀ / max(f_in, 0.05) ≥ 2)

| seed | B0: f₁₀ / f₂₀ / interior, R | young | intermediate | old |
|---|---|---|---|---|
| témoin | 90 / 95 / 73 %, 1.30 | 77 / 87 / 67 %, 1.29 | 39 / 58 / 64 %, 0.92 | 31 / 43 / 38 %, 1.12 |
| 42 | 78 / 82 / 59 %, 1.38 | 66 / 75 / 50 %, 1.48 | 37 / 43 / 44 %, 0.97 | 29 / 29 / 17 %, 1.73 |
| 1 | 78 / 82 / 43 %, 1.90 | 64 / 73 / 40 %, 1.84 | 22 / 39 / 24 %, 1.66 | 6 / 9 / 8 %, 1.16 |
| 9 | 73 / 78 / 74 %, 1.05 | 57 / 68 / 65 %, 1.05 | 21 / 35 / 48 %, 0.73 | 11 / 16 / 24 %, 0.69 |

- **The declared index fires on no world**: the coastal band is mountain (f₂₀ 68–95 % at B0 and young), but so is the interior.
- The rampart the author sees is a matter of steepness, not of the mountain class (A3).

### A3 — the walls and facets of the mountain class (2 048²)

| seed | | walls 28°: 0–10 / 10–20 / > 20 km | facets: 0–10 / 10–20 / > 20 km | walls: cover / basement / shield |
|---|---|---|---|---|
| témoin | B0 | 26.0 / 14.4 / 4.7 % | 28.5 / 29.9 / 14.9 | 0.9 / 13.4 / 18.4 % |
| | old | 6.2 / 5.8 / 4.9 % | 18.3 / 16.5 / 12.4 | 1.8 / 5.5 / 4.7 % |
| 42 | B0 | 20.0 / 13.1 / 5.9 % | 29.2 / 28.2 / 15.0 | 1.1 / 13.1 / 14.6 % |
| | old | 4.8 / 5.7 / 6.0 % | 17.5 / 16.5 / 15.2 | — / 5.7 / 5.0 % |
| 1 | B0 | 22.7 / 15.0 / 5.2 % | 29.0 / 27.4 / 14.9 | 1.2 / 16.3 / 12.3 % |
| | old | 4.1 / 1.2 / 0.6 % | 12.1 / 6.4 / 4.1 | — / 1.6 / 0.0 % |
| 9 | B0 | 18.9 / 13.0 / 5.8 % | 27.7 / 25.5 / 16.2 | 1.3 / 12.4 / 11.2 % |
| | old | 9.7 / 3.9 / 0.9 % | 19.8 / 11.2 / 5.2 | — / 2.6 / 0.1 % |

- **At B0, the walls are ×3–5 more frequent near the coast and ×10 more frequent on the basement than on the cover.**
- The facets double near the coast.
- At the old age the contrast fades.
- So the « rampart » is the margin's Davis-Suppe wedge (A1) exposed as hard basement (the cover goes first where the relief is
  highest) and cut into walls.

## C — the effect of the rise (K2; 4 seeds × 4 ages)

### The coast (torus)

Before → after, as dimension / length ratio / aligned share / block index (Corsica 1.12 / 1.37 / 0.27 / 0.78). The bands:
1.06–1.18 / 1.15–1.60 / 0.22–0.32 / 0.58–0.98.

| seed | B0 | young | intermediate | old |
|---|---|---|---|---|
| témoin | 0.98 / 0.90 → 1.05 / 1.11 | 0.99 / 0.90 → **1.16 / 1.60** | 1.00 / 0.94 → 1.29 / 2.28 | 1.01 / 0.95 → 1.32 / 2.42 |
| 42 | 0.98 / 0.89 → 1.22 / 1.87 | 0.98 / 0.89 → 1.24 / 2.00 | 0.99 / 0.93 → 1.33 / 2.51 | 1.00 / 0.96 → 1.38 / 2.83 |
| 1 | 0.98 / 0.90 → 1.21 / 1.81 | 0.99 / 0.93 → 1.22 / 1.89 | 0.99 / 0.95 → 1.33 / 2.55 | 1.00 / 0.96 → 1.39 / 2.94 |
| 9 | 0.99 / 0.92 → 1.21 / 1.81 | 0.99 / 0.94 → 1.27 / 2.16 | 1.00 / 0.98 → 1.36 / 2.72 | 1.01 / 1.01 → 1.39 / 3.01 |

- **Inside the bands after the rise: the témoin's young only** (3 of 4 measures).
- On every other world the dimension and the length ratio overshoot Corsica.
- The aligned share rises to 0.34–0.36 (outside the band) and the block index to 0.9–1.1.
- **Sensitivity Δ = 60 m** (the témoin, old): 1.01 / 0.95 → 1.27 / 2.13, with 23 drowned valleys. Still past Corsica.

### The drowned valleys, the islands, the land

| | B0 | young | intermediate | old |
|---|---|---|---|---|
| drowned valleys ≥ 2 km per world | 3–8 | 2–7 | 1–4 | 1–4 |
| longest (km) | 3–184 | 13–113 | 33–57 | 48–251 |
| the land lost against F171 | 0.7–7.5 % | 1.7–6.2 % | 4.7–9.0 % | 7.2–18 % |
| the classes' largest move | 0.2–7.8 pt | 1.7–5.5 pt | 3.1–5.0 pt | 5.6–9.0 pt |

- **The long « valleys» are flooded lowlands** (`images/f172_coast_25km.png`, `f172_before_after.png`): shallow seas over the
  plains with dendritic interfluve fingers, not rias.
- **Spearman(length, the river's area)**: −0.5 to 1.0 over the 11 worlds with ≥ 3 valleys (≥ 0.4 on 7).
- **Islands ≥ 1 km²**: 1–2 → 1–5.
- **The shelf** (the 0 to −134 m zone at the coast):
  - median 0.4–9.3 km (the témoin at B0 0.4 km; the old ages 4.2–7.0 km), p90 0.8–52 km;
  - it is drowned land, not a marine shelf: under Airy the ocean drops to −5.5 km beyond it (`f172_shelf_profile.png`).
- **The lakes** after the rise: 2–24 per world, 247–6 824 km².
- **Texture and regularity** after the rise:
  - the slope guard holds (p50 0.26–0.45);
  - σ_θ at 391 m 23.7–26.0° (Corsica 29.9°); R2_g 0.49–0.82.
- **The rampart map**: `f172_rampart_map.png` (young; red = mountains within 20 km of the coast).
- **Cost** (CPU per world):
  - the rise adds 0–10 %: 443–475 → 455–501 s on the first two seeds;
  - the last two seeds ran at ~1 050 s for both variants. The same worlds cost 455–479 s at F171, so this is the machine's
    state, not the work.
  - The rise bench took 6 h of wall time.

## D — the screening (33 seeds: 1–32 and the témoin; physics 256²; K2; the rise)

- **The selection**: a continent on 33 / 33 seeds.
  - 16 plates on 25, 20 on 7, 24 on 1.
  - Compact family 22, BFS 11.
- **The distributions** (median, deciles 10 … 90):

  | descriptor | young (20 Myr) | old (300 Myr) |
  |---|---|---|
  | mountain % | 70.7 (60.6 … 80.3) | 28.9 (11.0 … 51.1) |
  | plain + plateau % | 10.7 (5.8 … 14.8) | 26.1 (11.6 … 50.8) |
  | median land (m) | 707 (490 … 1 055) | 222 (93 … 476) |
  | peak (m) | 4 711 (4 246 … 5 636) | 3 957 (1 457 … 5 067) |
  | the rampart index R | 1.3 (1.0 … 1.5) | 0.8 (0.6 … 1.3) |
  | mountain regions > 100 km² | 3 (1 … 5) | 8 (4 … 13) |
  | the perimeter's convergent share | 17 % (0 … 33 %) | (the same, from C1) |

- **Europe's three bands** (mountains 16.5–37.2 %, plain + plateau 37.5–84.5 %, median 98–390 m): **5 / 33 seeds (15 %)**, all at
  the old age: seeds 8, 18, 19, 23, 29.
  - Separately at old: 17 seeds in the mountain band, 10 in the plain + plateau band, 23 in the median band.
  - The old profile's limiting clause is plain + plateau (23 seeds under 37.5 %).
- **The rampart** (the declared index): 0 % at young (no > 80 % bias flag). The index does not capture the walls (A3).
- **Spearman(the perimeter's convergent share, R)**: −0.08 young, 0.17 old.
- **Cost**: 4 494 s CPU in all, 136 s per seed; the selection 57 %.
- **The catalogue**: `catalogue.tsv` (one row per seed and age) and `images/f172_catalogue.png` (one row per seed, in the table's
  order: young shaded | young classes | old shaded | old classes).

## R

| rule | result |
|---|---|
| the rise's success (≥ 3 seeds × ≥ 2 ages closer to Corsica, the classes < 3 points) | **not met**: 1 world closer (the témoin, young); 0 seeds with 2 ages; the classes up to 9 points |
| no effect (Δdimension < 0.01 everywhere) | no: +0.07 to +0.39 everywhere; the screening ran |
| the screening's bias flag (a rampart on > 80 % at young) | no: 0 % with the declared index |

## The reading (hypotheses to check)

1. **The rise works and over-works**: the coast is past Corsica on 15 / 16 worlds.
   - The lowstand grades the lowlands below 0 m.
   - With no deposition and no shelf they all drown, into dendritic fringes.
   - Deposition (Yuan et al. 2019, queued) or a mean Quaternary level would fill or shorten them. The Δ = 60 m sensitivity
     alone does not suffice.
2. **The « rampart » is the margin's wedge in hard basement**:
   - Σ Δs_DS is highest at the coast whatever the boundary there (68–88 % of the coast is no boundary);
   - the cover is stripped first where the relief is highest;
   - the basement cuts into walls (×10 the cover's).
   - It shows at B0 and young and fades by the old age. It is a C1 (b) question: the wedge's reach along passive margins.
3. **The screening reads Europe's layout on 15 % of the seeds at the old age.** The limiting clause is the plains. The catalogue
   is the author's to browse.

## Predictions

### The reviewer's

| held | refuted |
|---|---|
| Y1 (+0.07–0.39 / +0.2–2.0, 16 / 16 worlds) | Y2's « the longest are the large rivers » (mixed, Spearman −0.5 to 1.0) |
| Y2's « > 10 km on ≥ 3 seeds » (4 / 4) | Y3 (the land lost 0.7–18 %, the classes up to 9 points) |
| Y4 (the shelf's median < 5 km on 12 / 16 worlds; but it is drowned land) | Y5 (Spearman −0.08 / 0.17) |
| the meta | Y6 (0 %) |
| | Y7 (15 %) |

### Mine

| held | refuted |
|---|---|
| Y5 refuted | Y1's ranges (far larger) |
| the rampart bias flag not firing | the success |
| A3 (the walls near the coast and on the basement) | Y2's « shorter at B0 / young » |
| the meta | Y3 |
| | Y4's 0.2–1.5 km |
| | Y6 « on the high side » (it was 0 %) |
| | Y7 |
| | 10–25 % of seeds without a continent (0 %) |
| | the selection 60–75 % of the cost (57 %) |
| | the rise's cost < 5 % (0–10 %) |

## Open, for the author (nothing chosen)

- **(a)** The rise as part of the continent profile: it drowns too much without deposition. Options:
  - wait for the deposition (Yuan et al. 2019);
  - a mean Quaternary level (Δ 60 m measured: still past Corsica);
  - a shorter T_q (the last glacial cycle, ~0.1 Myr, would carve less).
  - Each is the author's choice, not a calibration.
- **(b)** The rampart: the wedge along passive margins, a C1 (b) item (the walls ×10 on the basement).
- **(c)** The catalogue: 5 seeds read Europe at the old age (8, 18, 19, 23, 29). Browsing it in the viz is queued.
- **(d)** The rampart index as declared misses the walls; a wall-based index is a candidate.

## At commit (2026-10-12)
- **The reviewer's errors, recorded**:
  1. A constant lowstand at −134 m during 2.58 Myr, whereas the Quaternary sea level oscillates and touches its minimum only
     briefly.
  2. The lowstand applied to the amplification levels, which have no duration.
  3. The hypothesis « rampart = a convergent perimeter », refuted: 68–88 % of the coast is no plate boundary; Spearman −0.08 /
     0.17.
- **(a)** The mechanism is kept, with its forcing corrected at F173.
- **(b)** The rampart goes to C1 (b), after F173's part C diagnosis.
- **(c)** No catalogue browser. The seeds 8, 18, 19, 23 and 29 become the measurement panel. The old témoin serves only the
  bit-for-bit negative controls.
- **(d)** A rampart index based on the walls.

**The author's decision, verbatim**: « non pas besoin d'un navigateur de graines. Juste utiliser ces graines comme témoin pour
la suite. J'ai rajouté spratt2016.pdf dans docs/refs. »
