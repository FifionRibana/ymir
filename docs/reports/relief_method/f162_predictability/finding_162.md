# Finding 162 — the cascade, round 4: the predictability measured against Corsica; ρ, π and R4

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production is
unchanged**: the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`, `checks_after.txt`).

**The verdict in three sentences.**
- **The instruments do not show N1 as regular as the eye does.**
  - Of the six, four pass their control. On those four, F161's N1 is within ×1.5 of Corsica on all.
  - It is more regular on two: confluence angles 22.7° against 27.4°, gullies' parallelism 0.63 against 0.47. It is
    MORE dispersed on the valley spacing: CV 0.96 against 0.68.
- **The variants move them little.**
  - ρ changes nothing.
  - π decorrelates the texture: coherence 3 % against Corsica's 13 %, λ −18 %.
  - **R4 brings the slopes to Corsica's** (p50 0.177 / 0.188 against 0.247, within ×1.5) but closes only 17–22 % of the
    1.5–6 km deficit.
- **Every variant is rejected at 2 048².** The closest, R4+π, fails only on the two 1.5–6 km octaves (39.7 / 58.4
  against 62.9 / 100.9 m). The tight crop shows what the instruments miss: straight D8 trunks and a fishbone of
  parallel tributaries on broad, flat interfluves.

## Partie 0

1. **F161 committed** (`ab9c962`), the author's default-seed edits left out (`config.rs`, `workspace.rs`).
   - Its ADR records:
     - the author's verdict (« trop prévisibles »);
     - the tectonic note in the queue (« des filins montagneux », no Corsica-like wide high land);
     - the reviewer's uncontested proposals: the talus at a quarter by default, k on the octave's peak, the trunks at
       2 048² and the D5 lakes to the lakes' work, the texture judged at the last level.
2. **The scales, declared** (`f162_declared.md` §0):
   - the relief is at the real scale (400 km, F67);
   - the ×7.5 is for the exported hydrology only;
   - every area threshold of the cascade is real.
   - **The test** `the_cascade_reads_no_signified_quantity`: the five cascade sources never name the scale ratio, and
     `hydrology`'s areas are cells × the real cell² exactly.
3. **Corsica at 49 m**, measured for the final level's day (below).
4. **Checks**:
   - before: the guard 6 / 6, lib 629, viz 33;
   - after: the guard 6 / 6 identical (field and lakes = banc on all six states), lib 629 → 631 (+2), viz 33.

## P — the instruments and their controls (`cascade::predict`)

The definitions are in `f162_declared.md` P. **One amendment, at its control and before any measurement on the
cascade or Corsica**: CV_λ is computed along the rows and along the columns apart, then averaged. Mixing the two read
the regular fishbone as bimodal (1.23 against the fractal's 0.61).

| instrument | fishbone (regular) | fractal (dispersed) | the control | Corsica 781 / 391 / 195 / 49 m |
|---|---|---|---|---|
| CV_λ, valley spacing | 0.000 | 0.609 | **passes** | 0.557 / 0.609 / 0.683 / 0.829 |
| CV_L, stream lengths | 0.015 | 1.022 | **passes** | 0.813 / 0.869 / 0.915 / 0.903 |
| CV_Rb, Horton bifurcation | — (no basin of order ≥ 3) | 0.351 | **fails: blind** | 0.203 / 0.215 / 0.283 / 0.265 |
| CV_Rl, Horton length | — | 0.376 | **fails: blind** | 0.409 / 0.477 / 0.449 / 0.472 |
| σ_θ, confluence angles | 4.2° | 30.6° | **passes** | 29.4 / 29.9 / 27.4 / 27.6° |
| P_s, spectral peakiness | 0.060 | 0.174 | **fails: blind** (it reads the fractal as peakier) | 0.045 / 0.073 / 0.109 / 0.168 |
| R2_g, gullies' parallelism | 0.988 | 0.569 | **passes** | — / 0.477 / 0.473 / 0.408 |

- **Horton**: the fishbone's network stops at order 2, so the ratios are undefined. The means are reported anyway:
  Corsica R_b 3.1–3.4, R_l 1.4–2.1.
- **P_s**: the fractal source (`SeededNoise`, lacunarity 2) has steps at its own octaves; the fishbone's V-grooves have
  harmonics that look like a power law. **P_s grows with resolution on Corsica itself** (0.045 → 0.168). It is not a
  regularity measure here.
- **Corsica's confluence mean is 60–75°**, against N1's 49–60°.

## Measured at 1 024² (Corsica 391 m) and 2 048² (Corsica 195 m)

**The four validated instruments**:

| | CV_λ | CV_L | σ_θ (°) | R2_g |
|---|---|---|---|---|
| **Corsica 195 m** | **0.683** | **0.915** | **27.4** | **0.473** |
| témoin 2 048² | 0.962 | 1.019 | 27.6 | 0.636 |
| F161 N1 2 048² | 0.963 | 0.887 | 22.7 | 0.626 |
| N1b | 0.964 | 0.893 | 22.7 | 0.626 |
| ρ | 0.950 | 0.891 | 22.2 | 0.621 |
| π | 0.818 | 0.861 | 22.9 | 0.590 |
| ρ+π | 0.827 | 0.860 | 22.1 | 0.598 |
| R4 | 0.787 | 0.890 | 21.7 | 0.536 |
| **R4+π** | **0.675** | 0.855 | 22.2 | **0.518** |
| (1 024²: Corsica 391 m) | 0.609 | 0.869 | 29.9 | 0.477 |
| (1 024²: N1b) | 0.900 | 0.867 | 23.5 | 0.767 |

- **No variant has a validated instrument outside ×1.5 of Corsica at 2 048²**, N1b included (the ranking: 0
  everywhere; Σ|ln ratio| N1b 0.84, ρ 0.84, ρ+π 0.70, **π 0.64**). So the best predictability variant is **π**.
- **At 1 024² the gullies are far more parallel** (N1b 0.77, ρ 0.78, R4 0.78 against Corsica's 0.48). That is the
  comb of F159–F161 at its level; it thins at 2 048².

**The rest at 2 048²** (Corsica 195 m in brackets):

| | slope p50 / p90 (0.247 / 0.501) | facets (6.24 %) | walls 28° (1.55 %) | λ km (1.60) | R8 (0.024) | coherent (13.0 %) | density (0.52) | peak |
|---|---|---|---|---|---|---|---|---|
| N1b | 0.118 / 0.404 | 8.83 | 0.63 | 1.64 | 0.047 | 15.0 | 0.65 | 3 118 |
| ρ | 0.120 / 0.403 | 8.39 | 0.62 | 1.60 | 0.046 | 13.5 | 0.64 | 3 114 |
| π | 0.127 / 0.408 | 6.91 | 0.65 | 1.35 | 0.046 | **2.98** | 0.61 | 3 130 |
| ρ+π | 0.128 / 0.407 | 6.65 | 0.65 | 1.33 | 0.044 | **2.96** | 0.62 | 3 139 |
| R4 | **0.177 / 0.471** | 9.27 | 1.85 | 1.45 | 0.083 | 5.6 | 0.60 | 3 123 |
| R4+π | **0.188 / 0.473** | 8.54 | 1.88 | 1.32 | 0.042 | **2.38** | 0.59 | 3 121 |

**The spectrum at 2 048²** (m; Corsica in brackets):

| | 0.4–0.8 km (15.3) | 0.8–1.6 km (34.4) | 1.6–3.1 km (62.9) | 3.1–6.3 km (100.9) | 6.3–12.5 km (159.8) |
|---|---|---|---|---|---|
| N1b | 14.2 | 25.5 | 33.7 | 46.3 | 127.9 |
| π | 14.9 | 26.4 | 34.4 | 46.3 | 127.9 |
| R4 | 15.8 | **38.9** | 38.8 | **58.3** | 127.8 |
| R4+π | 16.5 | 39.5 | 39.7 | 58.4 | 127.6 |

- **R4 closes 17 % (1.6–3.1 km) and 22 % (3.1–6.3 km) of N1b's deficit**, and overshoots the 0.8–1.6 km octave (+13 %).
- **The 6–12.5 km octave never moves** (127.6–127.9 m): it is the physics level's.

**The calibration** (D1, ≤ 9 trials):
- 512² and 1 024² stay at the peak, k ≈ 1e3 (R4 at 512²: 562); the octave's peak is 21–23 m at 512² and 19–21 m at
  1 024².
- 2 048² reaches the target: k = 178 (N1b, ρ, π, ρ+π) and 316 (R4, R4+π).

**Per level, everything holds**: drift −0.43 to +0.09 %; peak 3 007–3 143 m. The trunks' p90 is ≤ 1.12 cell before
2 048². At 2 048² it is 1.41–3.20, the lakes' work.

**The π share** of the final 0.4–12.5 km variance: **0.38–0.43 %**, far under the 10 % limit.

**Cost**: a chain with its calibration 212–238 s; N1b at 2 048² alone 22.5 s (erosion 6.0, talus 12.0, deposition
4.1); the bench 1 657 s.

## R — the stop rule at 2 048² (against Corsica 195 m)

| variant | rejected on |
|---|---|
| F161 N1, N1b, ρ, π, ρ+π | slope p50 (0.118–0.128 against 0.247), the octaves 1.6–3.1 km and 3.1–6.3 km |
| **R4, R4+π** | **only the octaves 1.6–3.1 km (38.8 / 39.7 against 62.9) and 3.1–6.3 km (58.3 / 58.4 against 100.9)** |

- **Never fired**: the facets (all within 1.5 × Corsica's 6.24 %); the walls; the slope p90; λ; the four validated
  predictability instruments; the π share.
- **Per level**: nothing fires (drift, peak, trunks before 2 048²).

## What the measurements and the images say (the causes are hypotheses to check)

1. **« Prévisible » is not in the four dispersion instruments at 2 048².**
   - N1 is no more than ×1.3 more regular than Corsica on any of them, and more dispersed on the valley spacing.
   - **The tight crop** (`f162_tight_25km.png`) shows what they miss:
     - long **straight trunks** on the D8 diagonals, over tens of km;
     - a **fishbone** of parallel, equally spaced tributaries;
     - the interfluves as **broad, flat, smooth surfaces**.
   - Corsica at the same detail is tortuous at every scale.
   - Hypothesis: the regularity is (a) the channels' straightness (no sinuosity instrument), and (b) the missing
     1.5–6 km relief, the flat interfluves the eye reads as repetition.
2. **ρ does nothing; π decorrelates.**
   - The hardness noise barely changes the result. Hypothesis: the implicit solver relaxes the channels to their
     receivers whatever k is, and the retargeting keeps the coarse geometry.
   - The perturbation shifts the heads: λ −18 %, coherence from 15 % to 3 % (now below Corsica's 13 %), facets 8.8 →
     6.9 %.
3. **R4 is the one that changes the relief.**
   - Slope p50 +50 %, within ×1.5 of Corsica; λ, R2_g and CV_λ closer.
   - But it closes only a fifth of the 1.5–6 km deficit, and the walls at 28° triple (0.63 → 1.85 %, still under 1.5 ×
     Corsica's).
   - The 6–12.5 km octave comes from the physics level and no amplification touches it.
4. **The comb is a 1 024² feature**: R2_g 0.77 there against Corsica's 0.48, thinning to 0.63 at 2 048².

## The predictions

**The reviewer's** (hypotheses, now checked):

| prediction | verdict |
|---|---|
| N1 more regular than Corsica on ≥ 3 instruments (CV smaller by ≥ ⅓) | **refuted**: two instruments at 17–32 %; CV_λ the opposite; CV_L equal |
| ρ alone reduces the gap on ≥ 2 instruments without exceeding Corsica's facets | **refuted**: no measurable effect; facets 8.4 > 6.2 % |
| π weaker than ρ | **refuted**: π is the stronger (CV_λ, coherence, facets) |
| R4 closes at least half the 1.5–6 km deficit, with a drift < 2 % per level | **refuted** on the deficit (17–22 %); the drift held (≤ 0.43 %) |
| the best combination passes R at 2 048² | **refuted** (R4+π fails on two octaves) |
| meta | **held** |

**Mine** (`f162_predictions.md`):

| prediction | verdict |
|---|---|
| P-P1: every instrument separates ×1.5 except the confluence angles | **refuted**: the angles pass; Horton undefined on the fishbone; P_s reversed |
| P-P2: Corsica ×1.5 more dispersed than N1 on 3–4 instruments | **refuted** |
| P-D1: k near 1e3 at 512² / 1 024², the target only at 2 048² | **held** |
| P-D2: ρ halves the gap on parallelism and CV_λ | **refuted** |
| P-D2: ρ adds < 5 points of facets | **held** (−0.4) |
| P-D3: π < ⅓ of ρ's effect | **refuted** |
| P-D3: π's share < 5 % | **held** (0.4 %) |
| P-D4: R4 closes 30–60 % | **refuted** (17–22 %) |
| P-D4: R4's drift < 1 % | **held** |
| P-D4: R4 raises the facets at 512² | not reported this round |
| P-D5: no variant passes R | **held** |
| P-D5: the slope p50 stays ÷1.5 below Corsica's | **refuted** for R4 / R4+π |
| P-D5: the best combination fails one predictability instrument | **refuted**: it fails none |
| P-C1: a chain 4–8 min | **refuted**, faster (3.5–4 min) |
| P-C1: the bench < 60 min | **held** (28 min) |
| meta | **held** |

## B — the build (gated)

- **`cascade::predict`**: the six instruments; `control_regular` (the fishbone), `control_fractal`.
  - Test `the_fishbone_reads_regular_and_the_fractal_does_not`.
- **`cascade::measure::radial_power`**: the ¼-octave radial periodogram.
- **`cascade::amplify`**:
  - `AmpConfig::{n1_rho, n1_pi_m, n1_recal_depth}`;
  - the N1 hardness as the solver's `k_field`, the π band noise, `recalage_r4`, `Chain::calibrate_peak`;
  - `upscale_full` carrying the level n − 2.
- **`cascade::hydro`** test `the_cascade_reads_no_signified_quantity`.
- **The viz**:
  - « ρ dureté bruitée », « π perturbation » and « R4 recalage 4×4 » can be ticked, with their amplitudes (a, ×
    Corsica's octave);
  - « Caler sur la Corse » uses the peak method;
  - the talus at a quarter by default from 1 024².
- **The bench** `f162_cascade.rs`, with its instruments shared in `tests/common/cascade_bench.rs`.

## I — the images (`images/`; the Copernicus notice in `images/NOTICE.md`)

| image | content |
|---|---|
| `f162_whole_{1024,2048}.png` | N1b, π, R4, R4+π, then Corsica at the same km per pixel |
| `f162_{A_chaine,B_plateau,C_bassin_lac4}.png` | rows 1 024² and 2 048², the same four columns, then Corsica's spine (100 km) |
| `f162_tight_25km.png` | the 25 km crop on the range's southern flank (512² origin (176, 352)) against Corsica's western flank |
| `f162_spectra.png` | the 2 048² octaves: N1b red, π orange, R4 green, R4+π blue, Corsica black |

**The author's four questions**, for each column:
1. Le schéma répétitif a-t-il disparu, ou s'est-il atténué ?
2. Est-ce que ça « fait bruit » ?
3. Les crêtes et les vallées moyennes ont-elles assez de relief ?
4. À côté de la Corse, laquelle fait le plus vrai ?

## Open, for the author's decision (nothing chosen)

- **(a) The regularity the eye sees**:
  - a sinuosity instrument (the channels' length / chord against Corsica's) to name it;
  - the D8 straightness of the trunks (the routing of the implicit solver; an MFD- or D∞-based receiver).
- **(b) The 1.5–6 km deficit**:
  - R4 closes a fifth;
  - the 6–12.5 km octave is fixed by the physics level, so the physics level itself (its D at 256²), or a recalage
    to n − 3, are the options.
- **(c) π or R4+π as the default**:
  - π is the best predictability variant;
  - R4+π passes everything but the two octaves.
- **(d) The tectonics' wide high lands** (the author's note): queued.

## At commit (2026-10-09)
- **R4+π becomes the default setting** of the cascade.
- **Instruments set aside** (blind on their control): Horton's ratios (CV_Rb, CV_Rl) and the spectral peakiness (P_s).
- **The reviewer's reading** (hypothesis to check):
  - the straight trunks at 0° and 45° are D8 paths of the physics level (1.5 km cells), kept by the retargeting;
  - the 6–12.5 km octave and the 1.5–6 km deficit come from that level too.
- **The « recalage n − 3 » track is set aside**: the fine levels cannot create forms of several km.
- F163 follows: redo the physics level (a network less aligned on the grid, more dissected), and name what the author
  sees (directions, sinuosity, flat interfluves).
