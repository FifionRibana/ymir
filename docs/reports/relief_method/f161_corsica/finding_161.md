# Finding 161 — the cascade, round 3: Corsica as the reference, N1 (n = 1, the area capped, the smooth retargeting, k calibrated on Corsica)

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ».

**Production is unchanged**: the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`,
`checks_after.txt`). One production file gained an inert field, `StreamPowerConfig::area_cap_cells`: `None`
everywhere, skipped in serialisation. The F159 hash pin still holds.

**The verdict in three sentences.**
- **The smooth retargeting does what it was meant to.** It separates « creuser les vallées » from « abaisser le
  continent »:
  - with it, the drift is −0.28 / −0.07 / −0.01 %;
  - without it, at the same k, the land collapses to a 200 m peak.
- **But N1 cannot reach Corsica's texture at 512² and 1 024².**
  - The 2–4 cell octave is not monotonic in k; it peaks at 40 % and 67 % of Corsica's. Only at 2 048² is it met (14.2
    against 15.3 m).
  - Each level's coarser bands stay those of the level before, so the physics level's deficit at 3–6 km is never
    filled: 47 against 101 m at 2 048².
- **The stop rule, now against Corsica, fires at 512²**, the first amplification level:
  - facets ×7, slope p50 ÷1.9, valley spacing ×2.6;
  - the physics level itself fails the same three.
  - By 2 048² the valley spacing (1.64 against 1.60 km), the facets (1.4×) and the coherence reach Corsica's. The
    slope p50 (÷2.1) and the trunks (p90 1.58) do not.

## Partie 0

1. **F160 committed** (`80f210d`). Its ADR now records:
   - **the stop rule's reference was badly chosen** (the reviewer): the témoin is the rejected production world, and
     it grows more faceted with the resolution itself (5 % at 256², 28 % at 2 048²);
   - **the untested track**: the block retargeting alone.
2. **Checks**:
   - before: the guard 6 / 6, lib 627, viz 33;
   - after: the guard 6 / 6 identical (field and lakes = banc on all six states), lib 627 → 629 (+2), viz 33.
   - **The author's uncommitted edits in the working tree** were present in both runs and are not in any of my
     commits: the default seed (42 → 10481999410520546993) in `config.rs` and `workspace.rs`.
   - The window's `LAST_LEVEL = 8192`, the author's edit, is kept in `ui/cascade.rs`.

## C — the Corsica reference

**The data**:
- Copernicus DEM GLO-30, five tiles N41–N43 / E008–E009, downloaded on 2026-10-09 from the `copernicus-dem-30m` AWS
  bucket (HTTP 200, 52 MB);
- licence and attribution in `docs/refs/REFERENCES.md`;
- the tiles stay outside the repository; the grids live in `data/corsica/` (git-ignored); the script is
  `prep_corse.py`.
- The tool was an isolated venv in the session's scratchpad (numpy, tifffile, imagecodecs). The author's Python is
  untouched.

**The preparation**:
- the tiles are PixelIsPoint, pixel (0, 0) at (lon₀, lat₀ + 1);
- a local metric frame (WGS84 metres per degree at each pixel's latitude);
- 48.83 m cells, then block means;
- **the sea** at −50 m;
- **Corsica alone**: 290 islets set to sea (13 820 cells).

**It checks out**:
- 3.64 M land cells at 48.8 m, i.e. **8 700 km²** (Corsica: about 8 680 km²);
- peak 2 661 m at 49 m cells (Monte Cinto 2 706 m, averaged);
- mean 572 m at every cell size.

**What Corsica measures** (`f161_bench_output.txt`):

| cell | slope p50 / p90 | > 28° | walls 28° | facets | R8 (windows) | coherent | λ | density | relief p50 |
|---|---|---|---|---|---|---|---|---|---|
| 1 563 m | 0.099 / 0.211 | 0 % | 0 % | 3.45 % | 0.288 (36) | 8.3 % | 9.08 km | 0.151 | 1 238 m |
| 781 m | 0.144 / 0.304 | 0.15 % | 0.07 % | 3.77 % | 0.073 (174) | 9.2 % | 4.86 km | 0.339 | 1 431 m |
| 391 m | 0.194 / 0.403 | 2.4 % | 0.70 % | 4.61 % | 0.025 (771) | 9.5 % | 2.74 km | 0.488 | 1 570 m |
| 195 m | 0.247 / 0.501 | 7.7 % | 1.55 % | 6.24 % | 0.024 (3 275) | 13.0 % | 1.60 km | 0.516 | 1 621 m |
| 98 m | 0.292 / 0.578 | 14.0 % | 2.17 % | 8.20 % | 0.007 (13 653) | 19.9 % | 0.98 km | 0.548 | 1 633 m |

- **The spectrum in km is one curve, whatever the cell**:
  - e.g. 3.1–6.3 km reads 96.5 / 100.1 / 100.9 / 100.8 m from 781 to 98 m cells;
  - the exception is each grid's own finest octave, attenuated by 12–16 % by the block mean;
  - Corsica, 98 m cells: 0.2–0.4 km 7.1 m · 0.4–0.8 km 17.3 · 0.8–1.6 km 35.6 · 1.6–3.1 km 63.3 · 3.1–6.3 km 101 ·
    6.3–12.5 km 159 · 12.5–25 km 162 m.
- **The slopes steepen by 15–45 % per halving of the cell** (more at the coarse cells); **λ halves with the cell.**
  - So Corsica's texture is self-similar over the range: valleys at every scale the grid resolves.
- **The scale reservation** (declared): Corsica is a texture reference up to about 20 km. It is not one for the large
  forms, the hypsometry's absolute values, the network's size, or the coast.

## D — the build (gated)

- **`StreamPowerConfig::area_cap_cells`** — inert in production: min(A, cap) in the incision term only. Production's
  `None` path is pinned by `the_f159_extraction_changes_no_output`.
- **`cascade::measure`** — the spectrum:
  - `octave_rms`: a 2D radix-2 FFT, octaves of wavelength, RMS per textured land cell (Parseval);
  - `octave_km`, `upsample_mask2`;
  - test `a_sine_lands_in_its_octave`, with a flat-field negative control.
  - It replaced F160's box-blur bands **before any measurement**: they overlap, and a 3-cell sine read mostly in their
    band 1.
- **`cascade::amplify`**:
  - the variants N1 / N1 sans plafond / N1 sans recalage;
  - `erode_n1`: the implicit solver in blocks of 2 iterations, the area capped, no uplift;
  - `recalage_smooth`: e = R(z) − z_prev, the bicubic of e subtracted from the land; the full-weighting R, centred;
  - the talus angle carried to the cell, tan 33° · (c / 30 m)^−0.2 · [0.8, 1.4] (PROXY);
  - the D5 mask carried by nearest neighbour; the deposition at F160's k (`k_dep`);
  - `Chain::calibrate_next` (D3), `Chain::octave0`;
  - test `n1_keeps_the_coarse_level_and_the_cap_bounds_the_trunks`, with negative controls: without the retargeting
    the restriction drifts, and without the cap the solver carves more.
- **`physics_level`** now returns the D5 mask (`PeakCalibration::lifted_mask`).
- **The viz**:
  - N1 / N1 sans plafond / N1 sans recalage, with **N1 and physics 256² the defaults**;
  - **« k du niveau »**, starting at the bench's calibrated value for this world, and **« Caler sur la Corse »**, which
    runs D3 in the worker;
  - **« Référence Corse »**: Corsica at the displayed level's cell, beside the cascade, half the side (200 against
    400 km), the same km per pixel, the same shading, rivers and lakes.
  - The worker test passes (the physics → a level → « + érosion » → the retargeting).

**Amendments, all before any measurement and stated in `f161_declared.md`**:
- the FFT spectrum, as above;
- **k₀ = 100, not 1e-3.** The first value was 5 orders too small (F159's per-step k·dt was 200); the new N1 test
  caught it (the solver carved nothing). The expansions went from 4 to 6.

## M — N1 against Corsica (and the témoin), per level

**The calibration** (D3, ≤ 8 trials):

| level | Corsica's octave 0 | trials (k → octave 0, m) | kept |
|---|---|---|---|
| 512² | 53.6 m | 1e2 → 14.0 · 1e3 → **21.5** · 1e4 → 8.1 · 1e5 → 7.3 · 1e6 … 1e8 → 7.29 | k = 1e3, **out of reach** (40 %) |
| 1 024² | 30.4 m | 1e2 → 11.6 · 1e3 → **20.4** · 1e4 → 6.3 · 1e5 … 1e8 → 5.05 | k = 1e3, **out of reach** (67 %) |
| 2 048² | 15.3 m | 1e2 → 11.3 · 1e3 → 25.5 · 3.16e2 → 18.4 · **1.78e2 → 14.2** | k = 178, within ±20 % |

- **The octave is not monotonic in k.** It rises, peaks near k = 1e3, then falls to a floor.
  - Hypothesis to check: a strong implicit step grades the channels flat; the retargeting gives back everything
    coarser than 4 cells; what remains of the new band then shrinks.
- The bisection's bracketing assumed monotonicity. At 512² and 1 024² it expanded upward to its cap, then kept the
  closest trial, as declared.

**N1, the judging instruments**:

| | physics 256² | 512² | 1 024² | 2 048² |
|---|---|---|---|---|
| slope p50 / p90 · Corsica | 0.063 / 0.174 · 0.099 / 0.211 | **0.075** / 0.215 · 0.144 / 0.304 | **0.088** / 0.295 · 0.194 / 0.403 | **0.118** / 0.404 · 0.247 / 0.501 |
| facets % · Corsica | **13.0** · 3.45 | **26.4** · 3.77 | **14.1** · 4.61 | 8.8 · 6.24 |
| walls 28° % (cells) · Corsica | 0 · 0 | 0 · 0.07 | 0.02 (29) · 0.70 | 0.63 (4 318) · 1.55 |
| λ (km) · Corsica | **17.6** · 9.1 | **12.6** · 4.9 | **4.24** · 2.74 | 1.64 · 1.60 |
| drainage density (km/km²) · Corsica | 0.17 · 0.15 | 0.52 · 0.34 | 0.72 · 0.49 | 0.65 · 0.52 |
| coherent windows % · Corsica | 18.2 · 8.3 | 40.7 · 9.2 | 33.6 · 9.5 | 15.0 · 13.0 |
| R8 · Corsica | 0.025 · 0.288 | 0.106 · 0.073 | 0.018 · 0.025 | 0.047 · 0.024 |
| relief p50, 12.5 km (m) · Corsica | 1 019 · 1 238 | 1 159 · 1 431 | 1 277 · 1 570 | 1 337 · 1 621 |
| peak / mean (m) | 2 873 / 392 | 3 034 / 394 | 3 080 / 394 | **3 118** / 395 |
| drift (restricted) | — | −0.28 % | −0.07 % | −0.01 % |
| trunks p90 (cells) | — | 0.71 | 1.00 | **1.58** |
| smooth retargeting RMS (m) | — | 348 | 294 | 35 |

**The spectrum, N1 / Corsica / témoin** (RMS m, at 2 048², 195 m cells):

| octave | N1 | Corsica | témoin |
|---|---|---|---|
| 0.4–0.8 km | 14.2 | 15.3 | 10.5 |
| 0.8–1.6 km | 25.5 | 34.4 | 27.7 |
| 1.6–3.1 km | **33.7** | 62.9 | 58.7 |
| 3.1–6.3 km | **46.3** | 100.9 | 138.0 |
| 6.3–12.5 km | 127.9 | 159.8 | 229.6 |
| 12.5–25 km | 187.8 | 162.3 | 260.6 |

- **The physics level sets the 3–6 km octave at 33 against 85 m**, and the retargeting keeps it. Each level adds
  only its own new octave, so **N1's deficit sits at 1.5–6 km**, the octaves of the levels where the calibration
  saturated.
- **The hypsometry** (deciles 9 / 35 / 84 / 146 / 224 / 327 / 470 / 684 / 1 039 m against Corsica's 51 / 117 / 215 /
  331 / 455 / 594 / 763 / 973 / 1 296) is the physics level's, unchanged at every level (`f161_hypsometry.png`).
  - Its lowlands are far lower than Corsica's.
  - **Not a texture measure, declared out of the reference's scope.**
- **The lakes**:
  - the D5 lakes (2 225–2 601 km²) are counted apart;
  - the others: 23 / 34 / 31 km² against Corsica's 75 / 65 / 53 km².
- **The peak rises** from 2 873 (physics) to 3 118 m. The retargeting puts back, on the ridges, the mean lowering the
  solver applied to the block. Inside R's band, outside the 2 850 ± 5 % I predicted.

### The variants (at N1's k; `f161_bench_output.txt`)

| | 512² | 1 024² | 2 048² |
|---|---|---|---|
| **N1 sans recalage**: peak / drift | **853 m / −92 %** | **200 m / −88 %** | 202 m / −12 % |
| **N1 sans plafond**: trunk p90 · facets | 0.71 · 25.9 | 1.00 · 13.6 | 1.58 · 8.7 |
| N1 (with the cap): trunk p90 · facets | 0.71 · 26.4 | 1.00 · 14.1 | 1.58 · 8.8 |
| **N1 talus au quart**: slope p50 / p90 · talus time | 0.075 / 0.215 · 1.1 s (budget unchanged below 1 024²) | 0.088 / 0.295 · **1.5 s** (N1 6.7) | 0.118 / 0.404 · **12.5 s** (N1 47.7) |

- **The retargeting is what holds the continent**: without it, the same k erodes it to sea level.
- **The cap changes almost nothing** at these k (the trunks identical; the facets −0.5 %): the implicit solver
  relaxes every channel onto its receiver anyway.
  - **The trunks' p90 of 1.58 at 2 048² is N1's, with or without the cap.**
- **The talus at a quarter changes the slopes by 0.0 %** and divides its time by 3.8–4.5.
  - The adapted angles (22–38° at 195 m) rarely bind: N1's slopes are gentle.

### Cost (an alert, not a stop)

| | |
|---|---|
| N1, 512² → 2 048² | **70 s** without the calibration |
| with the calibration (its trials at 2 048² are 75 s each) | **396 s** |
| per iteration at 2 048² | erosion 823 ms (the priority flood; 6 iterations = 5 s) · talus 7.9 ms (6 000 = 48 s) · deposition 27 ms |
| **ESTIMATE to 8 192²** | erosion 115 s · **talus 953 s** · deposition 55 s → **19.9 min** |

- The talus is responsible, at its F160 budget. At a quarter it costs nothing in quality here (above).

### R — the stop rule, against Corsica

| chain | where R fires | criteria |
|---|---|---|
| N1 | **512²** | facets 26.4 > 1.5 × 3.77; slope p50 0.075 against 0.144; λ 12.55 against 4.86 km |
| | 1 024² (it would fire again) | facets 14.1; slope p50 0.088 / 0.194; λ 4.24 / 2.74 |
| | 2 048² | slope p50 0.118 / 0.247; trunks p90 1.58 |
| every variant | 512² | the same criteria |
| N1 sans recalage | 512² | adds the drift (−92 %) and the peak (853 m) |

- **The physics level (256², not subject to R) already fails three criteria**: facets 13.0 against 3.45; slope p50
  0.063 against 0.099; λ 17.6 against 9.1 km.
- **Never fired**: the walls (under 1.5 × Corsica or ≤ 20 cells, except 2 048² at 0.63 against 1.55, under); the
  drift; the peak.

## What the measurements say (the causes are hypotheses to check)

1. **The smooth retargeting is the right separation.** Valleys carve, the continent holds:
   - drift −0.01 to −0.28 %;
   - without it the land goes to sea level at the same k.
2. **N1 does make valleys at each level.**
   - λ goes 17.6 → 12.6 → 4.2 → 1.6 km, a new generation per level, and **meets Corsica at 2 048²** (1.64 against
     1.60 km).
   - The facets fall 26 → 14 → 8.8 %, against Corsica's 6.2 % at 2 048².
   - The coherence falls to Corsica's (15 against 13 %). **The comb thins out as the levels add.**
3. **But it does not reach Corsica's amplitude at 512² and 1 024²**, and the retargeting freezes that.
   - The 2–4 cell octave is not monotonic in k and peaks at 40 / 67 % of Corsica's.
   - Since each level only adds its own octave, the **1.5–6 km deficit (about ½ of Corsica's) survives to 2 048²**.
   - So do the slope p50 (÷2) and the gentle lowlands.
   - **The physics level is where the 3–6 km octave is born, and it is born at 0.4 × Corsica.**
4. **The trunks drift at 2 048²** (p90 1.58), with or without the cap. Hypothesis to check: the 2 048² level's finer
   network re-routes the 100 km² trunks across the coarse level's flats (lakes, D5 basins).
5. **The cost is the talus's**, which here does nothing measurable. At a quarter, 2 048² costs 22 s instead of 57 s.

## The predictions

**The reviewer's** (hypotheses, now checked):

| prediction | verdict |
|---|---|
| k calibrates to ±20 % at every level, rising at the fine levels | **refuted**: only 2 048² is reached; k falls (1e3, 1e3, 178) |
| λ falls from level to level | **held** (12.6 → 4.2 → 1.6 km) |
| λ within ×1.5 of Corsica up to 1 024² | **refuted** (×2.6, ×1.55); met at 2 048² |
| the drift < 0.5 % with the retargeting | **held** |
| without it, at the same k, the peak falls by > 20 % | **held** (−72 to −93 %) |
| without the cap, the trunks' p90 > 1.5 at one level at least | **held** (1.58 at 2 048²), but N1 with the cap gives the same 1.58: the cap is not the cause |
| the inherited facets < 1.5 × Corsica from 1 024² | **refuted** (3.1× at 1 024²); held at 2 048² (1.4×) |
| the talus at a quarter: slopes < 10 %, time ÷4 | **held** (0.0 %; ÷3.8–4.5) |
| cost to 2 048² < 3 min | **held** without the calibration (70 s); 6.6 min with it |
| meta | **held** |

**Mine** (`f161_predictions.md`):

| prediction | verdict |
|---|---|
| P-C1: slope p50 / p90 ≈ 0.20 / 0.45 at 391 m | **held** for the p50 (0.194); the p90 0.403 is 11 % lower |
| P-C1: +15–25 % per halving | **refuted** at the coarse cells (+33–45 %); held at the finest halving (+15–18 %) |
| P-C2: the spectrum in km one curve ±10 % | **held**, except each grid's finest octave (−12 to −16 %, the block mean) |
| P-C3: facets 10–20 % | **refuted** (3.5–8.2 %) |
| P-C3: walls ≤ 1 % | held to 391 m, **refuted** at 195 / 98 m (1.5 / 2.2 %) |
| P-C4: R8 at the noise floor | **held** (≤ 1.7 × the floor) |
| P-N1: calibrated at every level, k rising | **refuted** |
| P-N2: λ falls | **held** |
| P-N2: within ×1.5 at 512² / 1 024² | **refuted** |
| P-N2: outside at 2 048² | **refuted** (it matches) |
| P-N3: drift < 0.5 % | **held** |
| P-N3: peak 2 850 ± 5 % | **refuted** (3 034–3 118) |
| P-N4: p90 within ×1.5, p50 below ÷1.5 at 512² | **held** |
| P-N5: facets < 1.5 × Corsica from 1 024² | **refuted** (from 2 048²) |
| P-N6: R fires at 512² on the slope p50 | **held**, with the facets and λ at the same level |
| P-V1: sans recalage, the peak < −20 % and the drift > 2 % at the first level | **held** |
| P-V2: sans plafond, trunks p90 ≤ 1.5 | **refuted** (1.58; N1 too) |
| P-V3: talus ¼, p90 < 5 % and time ÷4 | **held** |
| P-T1: < 3 min without the calibration; the calibration 5–15 min | **held** (70 s; 396 s) |
| P-T2: ESTIMATE > 30 min, erosion largest | **refuted** (19.9 min, the talus) |
| meta | **held** |

## I — the images (`images/`; the Copernicus notice in `images/NOTICE.md`)

| image | content |
|---|---|
| `f161_N1_vs_Corse_{256,512,1024,2048}.png` | N1 (400 km, 1 024 px) beside Corsica (200 km, 512 px), the same km per pixel, north up, rivers of order ≥ 2 and lakes |
| `f161_{A_chaine,B_plateau,C_bassin_lac4}.png` | row 1, N1 at 256² (physics), 512², 1 024², 2 048² on F160's crops (100 km); row 2, Corsica's spine window (100 km) at the same cells; row 3, N1 sans recalage at 512²–2 048² |
| `f161_spectra.png` | RMS (m) against the octave's central wavelength (km), log-log. N1 in colour (red physics, orange 512², green 1 024², blue 2 048²), Corsica black (four cell sizes, superposed), the témoin grey |
| `f161_hypsometry.png` | the land's hypsometric curves: N1 per level in colour, Corsica (195 m) black |

**Questions for the author's look**, level by level (the images, or the viz: « Cascade » → « Calculer la physique »
(256²) → « Niveau suivant », with « Référence Corse » ticked):
1. Crêtes en facettes, murs, promontoires, « pinceau » ? (The facets fall from 26 % to 8.8 %; on crop A at 2 048² the
   flanks carry dense parallel gullies.)
2. Stries, terrasses droites, blocs ? (The coast is still the 64² blocky one.)
3. Les vallées se ramifient-elles à chaque niveau, comme sur la Corse ? (λ: 12.6 → 4.2 → 1.6 km against Corsica's
   4.9 → 2.7 → 1.6.)
4. Les grands fleuves restent-ils en place, et les hauteurs sont-elles justes ? (Trunks p90 1.58 at 2 048²; peak
   3 118 m.)
5. À côté de la Corse, est-ce que ça « fait vrai » ? (The measures say: the fine texture yes at 2 048², the 2–6 km
   ridges no.)

## Open, for the author's decision (nothing chosen)

- **(a) The 1.5–6 km deficit**, born at the physics level (0.4 × Corsica at 3–6 km) and frozen by the retargeting.
  Options:
  - a retargeting that lets each level also amplify the octave below its new one;
  - a physics level with a stronger fluvial term at its cell scale;
  - or calibrating on the two finest octaves.
- **(b) The non-monotonic calibration**: a calibration on the k at the octave's peak, or more iterations at a smaller
  k (the iteration count is not calibrated yet).
- **(c) The talus** costs 80 % of 2 048² and does nothing measurable: at a quarter, or off, at these levels.
- **(d) The trunks at 2 048²** (p90 1.58) and the D5 lakes: the lakes' work.

## At commit (2026-10-09): the author's verdict, recorded as given


- « Déjà beaucoup mieux. On voit apparaître les structures de plus petite échelle tout en conservant l'existant. »
- « Si on peut critiquer, les structures produites sont trop "prévisibles". On sent une régularité importante. On voit
  un schéma se répéter, alors qu'en observant la Corse, on peut voir moins ce détail répétitif. »
- **For a later work, the tectonics (in the queue)**: « Je n'ai vu sur aucune seed pour l'instant […] des terres assez
  larges élevées comme la chaîne de montagne de la Corse. J'ai que des filins montagneux. Mais c'est clairement un
  chantier postérieur. »

**Proposed by the reviewer, not contested by the author** (decisions for F162 on):
- the talus at a quarter by default;
- k calibrated on the octave's peak;
- the trunks at 2 048² and the D5 lakes go to the lakes' work;
- **the texture is judged at the last level**: the intermediate levels are reported, while the drift, the peak and the
  trunks stay checked at every level.
