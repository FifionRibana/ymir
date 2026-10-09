# Finding 163 — the cascade, round 5: the physics level redone; the directions, the sinuosity and the flat interfluves named

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production is
unchanged**: the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`, `checks_after.txt`).

- Declared before any measurement: `f163_declared.md`, with two amendments made at their controls (below).
- Predictions: `f163_predictions.md`.
- Bench: `f163_bench_output.txt` (27 min); the physics-level crops `f163_crops_output.txt`.

**The verdict in four sentences.**
- **The three instruments pass their controls and name what the author saw.**
  - At 2 048², the cascade's order-≥ 2 rivers are aligned on the grid at 40–48 % (Corsica 20.6 %; isotropic directions
    read 22.2 %).
  - Its flat interfluves are 9–15 % (Corsica 4.3 %).
  - Its trunks are 30–50 % less sinuous in excess (S − 1).
- **Redoing the physics level moves the network little.**
  - The random receiver changes the directions at 256² (34 % against 42 %) but nothing else.
  - The roughness is erased by the equilibrium (1.3 % of the final variance).
  - **The MFD exponent is the lever**: at p = 6 the valley spacing falls from 17.6 to 13.4 km and the facets from
    13 % to 5 %. It is still short of Corsica's 9.1 km.
- **The best variant (roughness + p = 6) through R4+π is the first chain whose every octave between 0.4 and 12.5 km is
  within ×1.5 of Corsica at 2 048²**: 1.6–3.1 km 42.7 against 62.9 m, 3.1–6.3 km 75.6 against 100.9 m. It is still
  rejected:
  - on the directions (44.8 % against 20.6 %);
  - on the sinuosity (S − 1 is 0.11 / 0.15 / 0.19 against 0.17 / 0.24 / 0.28);
  - on the flat interfluves (8.8 % against 4.3 %).
- **The eye agrees with the instruments and adds one thing**: the trunks are still straight lines at 0° and 45° over
  ~10 km in both chains, and p = 6 brings back a comb of parallel gullies on the range's flanks (crop A), which the
  facet instrument does not see.

## Partie 0

1. **F162 committed** (`dca1dcb`), the author's default-seed edits left out (`config.rs`, `workspace.rs`).
   - Its ADR records the reviewer's uncontested decisions:
     - R4+π becomes the default;
     - Horton's ratios and the spectral peakiness are set aside (blind on their control);
     - the reviewer's reading (a hypothesis): the straight trunks are D8 paths of the physics level, kept by the
       retargeting, and the 6–12.5 km octave and the 1.5–6 km deficit come from that level;
     - the « recalage n − 3 » track is set aside.
2. **Checks**:
   - before: the guard 6 / 6 (field and lakes = banc on all six states), lib 631, viz 33;
   - after: the guard 6 / 6 identical (field and lakes = banc on all six states), lib 631 → 633 (+2: the random
     receiver's test and the controls' test), viz 33.

## I — the three instruments (`cascade::planform`) and their controls

The definitions are in `f163_declared.md` I. **Two amendments, both at their controls and before any measurement on
Corsica or the cascade**:
- **I1's « high » control was changed from a funnel to a plane at 10°.** D8 on a cone converges on the 22.5° bisectors,
  zig-zag lines whose 8-step chords read unaligned. This was seen by reasoning, before any run.
- **I3's terraced field: each tread tilts 1 % toward its valley.**
  - As first built it read F = 18 %: the valley floor fell half the plane, so the V widened to ~50 cells.
  - With the floor fixed, it read 38 %: the flats' routing made order-2 streams on the perfectly flat treads.

**I2 is read on the excess S − 1, not on S** (declared before measuring): ×1.5 on S itself cannot reject a straight
line against Corsica. **For the author to confirm or refuse.**

| control | 256² (1 563 m) | 512² (195 m) |
|---|---|---|
| **A_dir**, D8 plane at 10° (high) | 97.2 % | 97.5 % |
| A_dir, fractal (reported) | 23.8 % | **43.1 %** |
| **A_dir, Corsica** | 35.6 % (**18 pieces only**) | **20.6 %** (1 222 pieces) |
| → A_dir | **passes** | **passes** |
| **S − 1**, straight valley at 0° | 0.000 (S_10) | 0.000 / 0.000 / 0.000 |
| S − 1, straight at 22.5° (the D8 floor) | 0.080 | 0.081 / 0.082 / 0.082 |
| S − 1, straight at 30° | 0.070 | 0.073 |
| S − 1, sine valley, read / continuous curve | — | 0.087 / 0.181 / 0.199 against 0.031 / 0.119 / 0.132 |
| **S − 1, Corsica** (S_2 / S_5 / S_10) | — / — / 0.115 | 0.172 / 0.239 / 0.282 |
| → S − 1 | **fails at 10 km: blind** (0.080 > 0.115 / 1.5); 2 and 5 km not read (< 4 cells) | **passes** at 2, 5 and 10 km |
| **F**, terraced field (high) | 88.6 % | 78.8 % |
| F, fractal (reported) | 6.2 % | 9.2 % |
| **F, Corsica** | 4.5 % | 4.3 % |
| → F | **passes** | **passes** |

What the controls teach:
- **Corsica's directions are isotropic at 195 m and 391 m** (20.6 / 20.8 %, against 22.2 %).
- **A synthetic fractal surface is not**: its D8 network reads 43 % aligned at 195 m. A smooth noise lets D8 run
  straight.
- **The D8 floor of the sinuosity is 0.07–0.08** (a straight line read at 22.5° or 30°). The sine valley reads
  0.06–0.07 above its continuous curve: the staircase adds to the true sinuosity.
- **At 1 563 m, Corsica's network is tiny** (18 pieces of order ≥ 2): its A_dir there is an anecdote.

**Corsica at three cells**:

| | 1 563 m | 391 m | 195 m |
|---|---|---|---|
| A_dir | 35.6 % (18) | 20.8 % (406) | 20.6 % (1 222) |
| A_dir, trunks ≥ 100 km² | 57.6 % (2) | 0 % (2) | 18.3 % (201) |
| S_2 / S_5 / S_10 | — / — / 1.115 | 1.112 / 1.161 / 1.194 | 1.172 / 1.239 / 1.282 |
| F | 4.5 % | 2.7 % | 4.3 % |

## D — the physics level at 256², against Corsica at 1 563 m

**The code facts (D0), read before any run**:
- the level already routes its AREA by MFD (p = 2), and its receivers are D8;
- the head threshold (A_c = 0.04 cell) and the diffusion (L_c = 0.4 cell) are both sub-cell, so neither can densify
  the network. The lever left is the MFD exponent.

**Corsica 1 563 m**: λ 9.08 km · facets 3.45 % · octaves 84.9 / 152.4 m · A_dir 35.6 % · F 4.5 %. The roughness's A is
0.5 × 84.9 = 42.5 m. The facet control at 256²: pyramids 87.3 %, paraboloids 2.4 % (it passes).

| variant | s (steps) | peak m | λ km | facets % | octaves 3.1–6.25 / 6.25–12.5 m | A_dir % | S_10 | F % | out of ×1.5 |
|---|---|---|---|---|---|---|---|---|---|
| **P0** | 4.6 (125) | 2 873 | 17.64 | 13.00 | 32.6 / 121.4 | 41.6 | 1.105 | 8.75 | 4 |
| **P-mfd** (random receiver, τ 0.5) | 13.4 (**300, cap**) | 2 719 | 16.62 | 12.26 | 37.2 / 123.0 | 33.5 | 1.102 | 9.20 | 4 |
| **P-bruit** (roughness 42.5 m) | 4.6 (118) | 2 881 | 17.84 | 12.61 | 32.6 / 121.8 | 35.9 | 1.103 | 9.03 | 4 |
| p = 1 | 5.2 | 2 837 | 19.84 | 16.99 | 28.8 / 112.1 | 36.6 | 1.098 | 9.53 | 4 |
| p = 1.5 | 5.8 | 2 929 | 18.83 | 14.85 | 31.1 / 121.8 | 39.1 | 1.100 | 8.81 | 4 |
| p = 3 | 8.7 | 2 751 | 16.18 | 9.26 | 38.7 / 116.8 | 37.1 | 1.101 | 8.89 | 4 |
| p = 4 | 5.7 | 2 797 | 15.00 | 7.38 | 44.3 / 117.0 | 31.8 | 1.106 | 9.03 | 4 |
| **p = 6 (P-dissection)** | 4.3 | 2 890 | **13.41** | **5.12** | 50.8 / 113.9 | 36.5 | 1.110 | 9.14 | 2 |
| D8 area | 2.1 (57) | 2 934 | 14.08 | 3.98 | 40.4 / 108.2 | 34.0 | 1.128 | 9.91 | 3 |
| P-mfd + P-bruit | 12.1 (300, cap) | **2 610: excluded** | 17.00 | 12.72 | 37.5 / 118.5 | 43.8 | 1.100 | 9.79 | 4 |
| P-mfd + P-dis | 10.8 (300, cap) | 2 911 | 13.44 | 6.65 | 57.3 / 127.7 | **21.3** | 1.112 | 8.76 | 3 |
| **P-bruit + P-dis (the best)** | 4.6 (134) | 2 863 | 13.40 | 4.83 | 48.8 / 114.8 | 35.8 | 1.109 | 9.26 | **2** |
| P-mfd + P-bruit + P-dis | 11.4 (300, cap) | 2 934 | 13.34 | 6.53 | 58.6 / 124.9 | 27.9 | 1.108 | 9.03 | 2 |

- The counts are among the measures whose control passed: λ, facets, the two octaves, A_dir and F. **S_10 is blind
  at 1 563 m.**
- **The best is P-bruit + P-dis** (2 out, Σ |ln| 2.29; P0 ranks 10th). It is out on the 3.1–6.25 km octave and on F.
- **P-dissection**: λ by p: 1: 19.84 · 1.5: 18.83 · 2: 17.64 · 3: 16.18 · 4: 15.00 · 6: 13.41 · D8: 14.08. p = 6 is
  kept, **not within ×1.25 of 9.08 km**.
- **The random receiver never reaches equilibrium**: a receiver drawn anew each step keeps |Δz| / U·dt above 0.01.
  - Its runs stop at the 300-step cap.
  - The peak's linear rescale (which assumes a steady state) then misses: 2 719 m, and 2 610 m with the roughness,
    which is excluded.
- **The roughness's share**: (a) A² over the final 3.1–12.5 km band is 9.5–11.7 %; (b) var(z − z_P0) / var(z) is only
  1.3 % alone, 3.2–5.0 % in the combinations. **The equilibrium erases it.**
- **What the physics level does NOT move**: F (8.8–9.9 % in every variant, against 4.5 %) and S_10 (1.10–1.13).

## R — the chains to 2 048², against Corsica at 195 m

Three chains:
- **N1**: F161's, as committed, k fixed; 161 s.
- **R4+π·P0**: F162's default, from P0; 477 s with its calibration.
- **R4+π·best**: from P-bruit + P-dis; 383 s.

Per level, everything holds in all three: drift ≤ 0.44 %, peak 3 003–3 182 m, trunks p90 ≤ 1.12 cells below 2 048².
The facet control at 2 048² passes (87.6 / 2.5 %).

| at 2 048² (Corsica) | N1 | R4+π·P0 | R4+π·best |
|---|---|---|---|
| slope p50 / p90 (0.247 / 0.501) | **0.118** / 0.404 | 0.188 / 0.473 | 0.196 / 0.478 |
| facets % (6.24) | 8.83 | 8.54 | 9.35 |
| walls 28° % (1.547) | 0.628 | 1.875 | 2.071 |
| λ km (1.60) | 1.64 | 1.32 | 1.37 |
| octave 0.39–0.78 km (15.3) | 14.2 | 16.5 | 15.9 |
| octave 0.78–1.56 km (34.4) | 25.5 | 39.5 | 39.3 |
| octave 1.56–3.12 km (62.9) | **33.7** | **39.7** | 42.7 |
| octave 3.12–6.25 km (100.9) | **46.3** | **58.4** | 75.6 |
| octave 6.25–12.5 km (159.8) | 127.9 | 127.6 | 119.2 |
| CV_λ / CV_L / σ_θ / R2_g (0.683 / 0.915 / 27.4° / 0.473) | 0.963 / 0.887 / 22.7° / 0.626 | 0.675 / 0.855 / 22.2° / 0.518 | 0.670 / 0.851 / 22.2° / 0.489 |
| **A_dir** (20.6 %) | **48.0 %** | **40.3 %** | **44.8 %** |
| A_dir, trunks ≥ 100 km² (18.3 %) | 53.3 % | 47.8 % | 48.4 % |
| **S_2 / S_5 / S_10** (1.172 / 1.239 / 1.282) | **1.088 / 1.125 / 1.166** | 1.129 / 1.171 / 1.201 | **1.112 / 1.150 / 1.185** |
| **F** (4.25 %) | **14.95 %** | **9.01 %** | **8.81 %** |
| π share | — | 0.38 % | 0.37 % |
| **R** | rejected (8) | rejected (4) | rejected (5) |

**The reasons**:
- **N1**: slope p50, the 1.6–3.1 and 3.1–6.3 km octaves, A_dir, S − 1 at 2, 5 and 10 km, and F.
- **R4+π·P0**: the 1.6–3.1 and 3.1–6.3 km octaves, A_dir, and F. **Its sinuosity passes** (S − 1 within ×1.5: 0.129
  / 0.171 / 0.201).
- **R4+π·best**: A_dir, S − 1 at 2, 5 and 10 km (0.112 / 0.150 / 0.185), and F. **Every octave passes.**

## The reading (hypotheses to check)

- **The directions at 2 048² are not the physics level's.**
  - At 256² P0's network is barely more aligned than Corsica's (42 % against 36 %, on 60 and 18 pieces).
  - At 2 048² every chain reads 40–48 %, twice Corsica's isotropic 20.6 %, and its trunks ≥ 100 km² 48–53 %.
  - **An 8-step piece at 2 048² is 1.56 km, one physics-level cell.** Each D8 step of a physics-level trunk is a
    straight segment at 0° or 45°, carried by the bicubic and the retargeting. At 2 048² each such segment reads as one
    aligned piece, while at 256² the 8-step chord averages the staircase and reads its true angle.
  - **That would explain why the random receiver cannot help**: it changes the order of the steps, not their
    directions.
  - The amplification levels' erosion (N1, D8 receivers too) carves its own streams along the same steps.
  - **The reviewer's reading (the trunks are physics-level D8 paths) is consistent with this, at the 1.5 km scale.**
- **The sinuosity trades against the octaves.** p = 6 gains the octaves and loses the sinuosity: 0.15 against P0's
  0.17 at 5 km. A more concentrated area makes straighter, deeper valleys.
- **The comb**: p = 6 (as D8) concentrates the area again. Crop A shows the parallel gullies relief_v3's MFD p = 2 was
  chosen to prevent (ADR Finding 10/11). The facets fall but R2_g does not see it at 2 048² (0.489).
- **The flat interfluves come with the physics level** (8.8–9.9 % at 256² in every variant) **and stay** (8.8–9.0 % at
  2 048² for R4+π). A candidate is the hypsometry: the level's land is low (deciles 9 / 35 / 84 m against Corsica's
  57 / 126 / 224 m), so its plains are flat. That is the tectonics' note again (no wide high land).

## I — the images (north up; rivers order ≥ 2; Copernicus notice in `images/NOTICE.md`)

**The images**:
- `f163_physics_256_crops.png`: the physics level alone, the land ×3. Row 1: P0, P-mfd, P-bruit; row 2: p = 6, P-mfd +
  p = 6, P-bruit + p = 6.
  - `f163_physics_256_corse.png`: Corsica at 1 563 m, the same scale.
  - `f163_physics_256.png`: the 13 variants, whole.
- `f163_tight_25km.png`, `f163_A_chaine.png`, `f163_B_plateau.png`, `f163_C_bassin_lac4.png`: at 2 048², R4+π·P0, then
  R4+π·best, then Corsica at 195 m.
- `f163_whole_2048.png`, `f163_spectra.png`.

**For the author's four questions, my reading** (the author's own is the one that counts):
1. **« Les grands cours d'eau sont-ils encore des lignes droites ? »** Yes, in both chains. In the 25 km crop, a 45°
   trunk runs straight for ~10 km and a 0° one for ~6 km. The best variant does not change this.
2. **« Les arêtes de poisson et les interfluves plats ont-ils reculé ? »**
   - The flat interfluves have not (9 % against 4 %).
   - The fishbone has grown on the range's southern flank with p = 6 (crop A, the middle column).
3. **« Le relief entre les vallées (quelques km) est-il assez marqué ? »** Better with the best variant: the 3–6 km
   octave is 75.6 m against P0's 58.4 m and Corsica's 100.9 m, within ×1.5 for the first time.
4. **« À côté de la Corse, est-ce que ça « fait vrai » ? »** Not yet. Corsica's main valleys bend and branch at every
   scale; the cascade's run straight between flat interfluves.

## Predictions

**The reviewer's**:
1. « P0 ≥ 2× more aligned, sinuosity ≥ 20 % lower »: **refuted as written, held in part**.
   - Directions: refuted at 256² (×1.17); at 2 048² ×1.96 for R4+π·P0 (just under) and ×2.33 for N1.
   - Sinuosity on S: refuted (−5.5 %). On S − 1: held (−28 % for R4+π·P0, −48 % for N1).
2. « P-mfd brings the directions within ×1.5 »: **held at 256²**, where P0 already was (×1.17 → ×0.94). Not carried
   to 2 048² (P-mfd was not the best).
3. « P-dissection brings λ and the 3–6 km octave closer and lowers the facets »: **held** (17.6 → 13.4 km;
   32.6 → 50.8 m; 13.0 → 5.1 %).
4. « The best combination misses no octave beyond ×1.5 at 2 048² »: **held**.
5. Meta: held (1 is refuted).

**Mine**:
- **C1**:
  - the plane passes as predicted, but reads 97 %, not 45–70 %;
  - Corsica 25–35 % at every cell: refuted (20.6 % at 195 m);
  - the fractal 22–30 %: refuted at 195 m (43 %).
- **C2**:
  - held: the straight line, Corsica's S_5 − 1 (0.239), blind at 1 563 m and passes at 195 m;
  - refuted: the 22.5° floor, by a hair (0.082 > 0.08).
- **C3**: held, except Corsica 5–15 % at 1 563 m (4.5 %).
- **M1** (P0 at 256²):
  - held: A_dir 40–55 %, and F outside ×1.5;
  - refuted: the ratio (1.17), S_10 − 1 (0.105) and F's range (8.8 %).
- **M2** (2 048²):
  - held: A_dir 40–55 %, F 5–15 %, and Corsica's S − 1 and F;
  - refuted: Corsica's A_dir (20.6 %), the cascade's S_5 − 1 (0.125 / 0.171), and « outside on all three » for
    R4+π·P0 (its sinuosity passes).
- **D1** (P-mfd):
  - held: within ×1.5, λ ±15 %, the cost per step +21 %;
  - refuted: the facets (−6 %).
- **D2** (P-bruit): held on its share (11.3 %); refuted on λ (+1 %) and the octave (0 %).
- **D3**: no p within ×1.25 held; « D8 facets ≥ P0's » refuted (4.0 %).
- **D4**:
  - « the best includes the random receiver »: refuted;
  - « rejected »: held;
  - « on the 1.5–6 km octaves »: refuted (they pass);
  - « on F »: held;
  - « A_dir closes < half at 2 048² »: held (it worsens, 44.8 against 40.3 %).
- **C (cost)**:
  - held: the physics variants (2.1–13.4 s), and the bench under 30 min (27.3);
  - refuted: the chains at ~4 min (6–8 min with the calibration).
- **Meta**: held.

## Cost

- A physics variant takes 2.1–13.4 s; the random receiver costs +21 % per step and runs to the 300-step cap.
- R4+π takes 383–477 s per chain with its calibration; the calibration is 330–400 s of it.
- The whole bench took 1 637 s.

## Open, for the author (nothing chosen)

- (a) **The straight trunks are a 1.5 km step problem, not a routing one** (hypothesis). Bending a trunk inside one
  physics cell needs the amplification. Candidates:
  - a random receiver in the amplification's N1;
  - a course that wanders on the upsampled field (meanders, Schott's multi-scale breach).
- (b) **p = 6 as the physics level's exponent**: it brings the octaves within ×1.5 but loses the sinuosity and brings
  back a comb on the flanks. The author's eye decides.
- (c) **The flat interfluves** come with the physics level and never moved. Is it the hypsometry (the tectonics' wide
  high lands, queued), or the level's own plains?
- (d) **I2 read on S − 1**: to confirm or refuse.

## At commit (2026-10-09)
- **The sinuosity is read on S − 1** (I2's declared reading, confirmed).
- **p = 6 is not adopted as is**: its octaves are within ×1.5, but the sinuosity drops and the comb comes back on the
  flanks (crop A).
- **The reading kept** (hypothesis to check):
  - the straight trunks are the physics level's D8 steps (1.56 km = 8 steps at 2 048²), kept by the upscaling and the
    retargeting;
  - a course does not bend by reordering its steps.
- **A notable control**: a smooth fractal surface reads 43 % aligned at 195 m. Noise alone does not make a natural
  network.
- F164 follows: bend the inherited trunks (a smooth, invertible warp of the upscaled field), arbitrate p, and
  attribute the flat interfluves.
