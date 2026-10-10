# Finding 164 — the cascade, round 6: a smooth invertible warp of the upscaled field, the exponent p arbitrated, the flat interfluves attributed

**Status: built, measured, gated.** Nothing is committed until the author's « feu vert ». **Production is
unchanged**: the guard 6 / 6, field and lakes, was run before and after (`checks_before.txt`, `checks_after.txt`).

- Declared before any measurement: `f164_declared.md`, with one precision written before measuring (how the settings
  are compared).
- Predictions: `f164_predictions.md`.
- Bench: `f164_bench_output.txt` (38 min, 11 chains).

**The verdict in four sentences.**
- **The warp is clean**:
  - invertible everywhere (min J 0.31–0.43);
  - it moves the hypsometry by ≤ 0.46 % of q9 and the coast by ≤ 0.33 % of the land;
  - **but the depression count is noisy in both directions** (−13 % to +19 %), and the 1 % clause fires on 3 of the 4
    warped chains at p = 2 (never at p = 4 or 6).
- **W-tous bends the trunks by a third, not enough.**
  - Their own alignment (pieces on ≥ 100 km²) falls from 47–48 % to 33–36 %.
  - The whole network's A_dir falls from 39–45 % to 33–36 %. It stays above the 30.9 % that ×1.5 of Corsica's isotropic
    20.6 % would allow.
  - W-phys alone does almost nothing at p = 2 (38.9 → 38.8 %).
- **p = 6 with W-tous is the best chain so far**: 2 R reasons only (A_dir 36.4 % and F 9.2 %). Every octave, the slopes,
  λ, the predictability, the sinuosity (S − 1 0.12 / 0.17 / 0.22) and the facets pass.
- **The flat interfluves are mostly low, in the cascade (81.5 % below 100 m) and in Corsica (78.6 %).**
  - About half of the cascade's excess is its hypsometry (too much land below 100 m).
  - A third is a higher flat rate within each band.
  - 37 % of its flat cells lie within 2 km of a D5 lake.

## Partie 0

1. **F163 committed** (`a4831cd`), the author's default-seed edits left out.
   - Its ADR records:
     - S − 1 confirmed;
     - p = 6 not adopted as is;
     - the reading kept as a hypothesis: the trunks are the physics level's D8 steps, and a course does not bend by
       reordering its steps;
     - the notable control: a smooth fractal reads 43 % aligned.
2. **Checks**:
   - before: the guard 6 / 6 (field and lakes = banc on all six states), lib 633, viz 33;
   - after: the guard 6 / 6 identical (field and lakes = banc on all six states), lib 633 → 634 (+1: the warp's
     invertibility test), viz 33.

## D — the warp (`cascade::amplify::{Warp, warp_field, jacobian, warp_grid}`)

**How it works**:
- **u_W(x) = u(x + d(x))**, the Catmull-Rom bicubic, right after the bicubic ×2.
- |d| ≤ A, a one-octave gradient noise of wavelength L, A and L in previous-level cells.
- The D5 mask is warped by nearest cell. The ocean, the D5 lift and the floor are computed on u_W.
- **The retargeting's references are R(u_W) and R(R(u_W))**, so the warp is kept, not pulled back.
- `AmpConfig::warp`, `None` = F159–F163.
- **A permanent test**: W1 and W2 have J > 0 everywhere at 512², 1 024² and 2 048² on two seeds; the negative control
  (A 2, L 2) folds; a zero warp is the identity.

**The warp's controls on the témoin** (before → after, per warp):

| warp | depressions | hypsometry max Δq / q9 | land cells | min J | J < 0.5 |
|---|---|---|---|---|---|
| W1 at 512² (p = 2 / 4 / 6) | 62 → 61 · 61 → 60 · **71 → 62** | 0.34–0.40 % | −0.32 to −0.33 % | 0.348 | 0.008 % |
| W1 at 1 024² (p = 2 / 4 / 6) | 132 → 129 · 150 → 148 · 153 → 152 | 0.12–0.19 % | +0.015 % | 0.425 | 0.003 % |
| W1 at 2 048² (p = 2 / 4 / 6) | **316 → 325 (+2.8 %)** · 346 → 344 · 367 → 360 | 0.02–0.06 % | +0.03 % | 0.339 | 0.002 % |
| W2 at 512² (p = 2) | **62 → 74 (+19 %)** | 0.46 % | −0.31 % | 0.328 | 0.20 % |
| W2 at 1 024² (p = 2) | **131 → 137 (+4.6 %)** | 0.09 % | −0.05 % | 0.306 | 0.18 % |
| W2 at 2 048² (p = 2) | 363 → 342 | 0.03 % | +0.03 % | 0.318 | 0.19 % |

- **The count moves both ways.** A resampled smooth field opens some sills and closes others, so the 1 % clause reads
  noise as much as creation.
  - It fires for p2 · W-tous W1 (at 2 048²), p2 · W-phys W2 and p2 · W-tous W2 (at 512²).
  - It never fires at p = 4 or 6.
- **The coast**: the first warp takes 0.32 % of the land at 512²; the later warps give back 0.02–0.03 %.
  - On the crops, W-phys and W-tous show a speckle of small light-blue cells along some coasts (crop A, the north
    coast). Hypothesis: warped bathymetry leaves small interior below-sea cells that the D5 lift turns into lakes.
  - Reported, not judged (later work).

## The chains (R4+π to 2 048², k calibrated on the peak per level)

**The settings at p = 2**:
- W1: 8 R reasons over its two chains, mean A_dir 37.2 %;
- W2: 10 R reasons (the depression clause twice), mean A_dir 35.5 %.
- **W1 is carried** to p = 4 and 6.
- W2 bends more (A_dir 32.7 % under W-tous, the lowest of the round), but it is the one that opens depressions.

| chain | A_dir (trunks) | S − 1 at 2 / 5 / 10 km | F | octaves 1.56–3.12 / 3.12–6.25 / 6.25–12.5 km | facets | R |
|---|---|---|---|---|---|---|
| **Corsica 195 m** | 20.6 % (18.3) | 0.172 / 0.239 / 0.282 | 4.25 % | 62.9 / 100.9 / 159.8 | 6.24 | |
| p2 · none | 38.9 (46.7) | 0.131 / 0.175 / 0.204 | 8.91 | **39.1** / **58.5** / 127.8 | 8.40 | 4 |
| p2 · W-phys W1 | 38.8 (39.0) | 0.128 / 0.172 / 0.212 | 8.84 | **39.7** / **57.0** / 121.3 | 7.23 | 4 |
| p2 · W-tous W1 | 35.6 (33.4) | 0.126 / 0.172 / 0.214 | 9.42 | 50.8 / **55.7** / 114.6 | 8.95 | 4 (+ depressions) |
| p2 · W-phys W2 | 38.4 (36.6) | 0.127 / 0.173 / 0.214 | 8.78 | **41.0** / **60.4** / 120.1 | 7.25 | 5 (+ depressions) |
| p2 · W-tous W2 | **32.7** (35.4) | 0.131 / 0.183 / 0.228 | 9.08 | 50.4 / **56.9** / 112.9 | 7.68 | 5 (+ depressions) |
| p4 · none | 41.9 (48.0) | 0.122 / 0.162 / 0.195 | 8.94 | **41.1** / 70.6 / 124.0 | 9.04 | 3 |
| p4 · W-phys W1 | 40.0 (42.7) | 0.124 / 0.169 / 0.213 | 8.76 | **41.4** / 68.4 / 117.5 | 7.89 | 3 |
| p4 · W-tous W1 | 36.2 (35.8) | 0.126 / 0.174 / 0.219 | 9.29 | 49.9 / **64.1** / 110.9 | 9.09 | 3 |
| p6 · none | 44.8 (48.4) | **0.112 / 0.150 / 0.185** | 8.81 | 42.7 / 75.6 / 119.2 | 9.35 | 5 |
| p6 · W-phys W1 | 40.2 (42.3) | 0.123 / 0.168 / 0.213 | 8.71 | 43.0 / 73.5 / 112.9 | 8.03 | 2 |
| **p6 · W-tous W1** | 36.4 (33.2) | 0.123 / 0.170 / 0.215 | 9.18 | 50.0 / 67.7 / 106.8 | 8.83 | **2** |

Bold marks a value outside ×1.5 of Corsica. The full reasons, the slopes and λ are in the bench output.

**What the table says**:
- **Every chain fails A_dir and F**:
  - A_dir ×1.6–2.2 of Corsica;
  - F ×2.0–2.2.
- **The warp's effects**:
  - **W-tous cuts the trunks' alignment by about a quarter to a third** (47–48 % → 33–36 %), and the whole network's
    by 3–8 points.
  - **W-tous raises the 1.56–3.12 km octave by 17–30 %** (39–43 → 50–51 m): the bends put relief at the 1.5–3 km
    scale.
  - It lowers the 6.25–12.5 km octave by 10 % (128 → 115 at p = 2).
  - At p = 2 and 4 that moves the failure from one octave to the other: p = 4 W-tous fails 3.12–6.25 km (64.1 against
    67.3).
  - **The sinuosity barely moves at p = 2** (S_5 − 1 0.175 → 0.172), but **at p = 6 the warp brings it back within
    ×1.5** (0.150 → 0.168–0.170). That makes p = 6 · W-phys and W-tous the two chains at 2 reasons.
- **p** sets the octaves (p = 6 passes them all, p = 2 fails 1.5–6 km) and costs sinuosity without the warp.
- Per level everything holds in all 11 chains: drift ≤ 0.61 %, the peak 2 909–3 182 m, the trunks below 2 048² ≤ 1.12
  cells.
- **At 1 024²** A_dir is 41–59 % against Corsica's 20.8 %. The warp's effect is the same there.

## F — the flat interfluves attributed (195 m)

| | p6 · W-tous W1 (2 048²) | Corsica |
|---|---|---|
| F | 9.18 % | 4.25 % |
| land 0–50 / 50–100 / 100–300 / > 300 m | 18.0 / 10.5 / 25.4 / 46.1 % | 9.8 / 7.9 / 19.6 / 62.7 % |
| flat cells in those bands | **66.7 / 14.9** / 14.2 / 4.3 % | **54.6 / 24.1** / 12.2 / 9.1 % |
| flat rate in those bands | 34.0 / 12.9 / 5.1 / 0.85 % | 23.7 / 12.9 / 2.7 / 0.62 % |
| flat cells below 100 m | **81.5 %** | **78.6 %** |
| flat cells 0–2 km from the coast | 39.6 % (rate 32.5 %) | 53.9 % (rate 13.5 %) |
| flat cells > 20 km from the coast | 18.1 % (rate 4.7 %) | 0.8 % (rate 0.3 %) |
| within 2 km of a D5 lake (10 % of the land) | 36.7 % of the flat cells (rate 33.8 %, the rest 6.5 %) | (no D5) |

- **By the declared rule both are « mostly low »** (≥ 60 % below 100 m, and the rate below 100 m 31× / 30× the rate
  above 300 m). **The flat cells are low, not spread out.**
- **The excess, decomposed** (4.9 points of F):
  - with Corsica's hypsometry and its own band rates the cascade would read 5.9 %;
  - with its own hypsometry and Corsica's band rates it would read 6.6 %;
  - so **~47 % of the excess is the hypsometry** (28.5 % of the land below 100 m against 17.7 %);
  - **~33 % is flatter low bands** (0–50 m: 34 % against 24 %; 100–300 m: 5.1 % against 2.7 %);
  - the rest is the interaction.
- **Two cascade traits Corsica has not**:
  - flat cells far inland (18 % of them > 20 km from the coast, against 0.8 %);
  - **a third of them on the D5 lakes' margins** (lifted floors and their flat rims).
- The hypsometry part points to the tectonics (no wide high land, the queued note), and the D5 part to the lakes' work.

## The reading (hypotheses to check)

- **A smooth warp of ≤ 0.5 cell bends the trunks inside the inherited steps but does not erase them.**
  - The fine levels' own D8 erosion and the retargeting keep re-cutting straight segments at 2 048² (A_dir at 1 024²
    stays 45–50 %).
  - A larger A is out of the brief's bound, and W2 already opens depressions.
- **The depression clause as written reads interpolation noise**: the count falls as often as it rises.
  - A count of NEW depressions (those whose bottom cell was not in a depression before), rather than a net count,
    would separate creation from noise.
  - Not done here; open.

## I — the images (north up; rivers order ≥ 2; Copernicus notice in `images/NOTICE.md`)

- `f164_tight_25km.png`, `f164_A_chaine.png`, `f164_B_plateau.png`, `f164_C_bassin_lac4.png`:
  - row 1: p6 · none, p6 · W-phys W1, p6 · W-tous W1, then Corsica;
  - row 2: W-tous W1 at p = 2, 4 and 6, then Corsica.
- `f164_whole_2048.png` (row 1) and `f164_spectra.png`.

**For the author's four questions, my reading** (the author's own is the one that counts):
1. **« Les grands cours d'eau sont-ils encore des lignes droites ? »** Less so with W-tous. In the 25 km crop the 45°
   trunk now wavers along its course, but it still reads as one diagonal over ~10 km. Without W and with W-phys it is
   the same straight line as F163.
2. **« Le peigne sur les flancs : présent, atténué, absent ? »**
   - Present at p = 6 with every placement (crop A, row 1).
   - Attenuated at p = 2 (crop A, row 2, first column).
   - The warp does not remove it.
3. **« Est-ce que la déformation se voit comme telle ? »** Not as twisted shapes at A = 0.5. What shows is the speckle
   of small lakes along some coasts with W-phys and W-tous.
4. **« À côté de la Corse, est-ce que ça « fait vrai » ? »** Closer on the numbers (2 reasons left), not yet on the eye:
   - Corsica's valleys bend at every scale;
   - the cascade's still run in long segments between flat lowlands and lake rims.

## Predictions

**The reviewer's**:
1. « W-tous brings A_dir under 31 % and S − 1 within ×1.5; W-phys less than half »:
   - A_dir **refuted** (35.6 % at W1, 32.7 % at W2);
   - S − 1 **held** (within ×1.5 at all three windows);
   - « W-phys less than half » **held** (0.1 against 3.3 points).
2. « < 1 % new depressions, hypsometry < 1 % »:
   - the hypsometry **held** (≤ 0.46 %);
   - the depressions **refuted** on 3 warps (+2.8 %, +19 %, +4.6 %) and held on the others.
3. « p = 4 with W-tous passes every octave and has less comb than p = 6 »: **refuted** on the octaves (3.12–6.25 km
   64.1 against 67.3). The comb, by the eye: slightly less than p = 6.
4. « Flat interfluves > 60 % below 100 m »: **held** (81.5 %).
5. Meta: held.

**Mine**:
- **The warp's controls**:
  - W1's min J 0.5–0.85: refuted (0.34–0.43);
  - W2's 0.15–0.6: held (0.31–0.33);
  - the negative control folds: held;
  - the depressions within ±1 % at every warp: refuted;
  - the hypsometry < 1 %: held;
  - the coast < 0.5 %: held (0.33 %).
- **Directions and sinuosity at p = 2**:
  - held: none ≈ 40 % (38.9), W-tous 30–37 % (35.6 / 32.7), not within ×1.5, W-phys 36–41 % (38.8 / 38.4) and less
    than half, W2 more than W1;
  - refuted: « S_5 − 1 rises to 0.18–0.24 under W-tous » (W1 0.172; only W2 reaches 0.183).
- **p**:
  - held: p = 6 passes every octave; p = 4 fails 1.56–3.12 km without W and with W-phys;
  - refuted: « with or without W » (with W-tous it passes that one and fails 3.12–6.25), and « the warp changes the
    octaves by ≤ ±10 % » (+17–30 % at 1.5–3 km).
- **R**: held (no chain passes; every chain fails F; all fail A_dir).
- **F**:
  - held: « mostly low »; Corsica > 50 % low (78.6 %); near D5 ≥ 2× the rest (33.8 against 6.5 %);
  - refuted: the 55–75 % range (81.5 %).
- **Cost**: the 70–95 min refuted (38 min: the chains took 3 min each, the calibration faster than at F163).
- **Meta**: held.

## Cost

- A warp is not timed separately: the whole upscale step, warp included, takes 3.1 s at 512² and ~21 s at 2 048² for
  the level.
- Each R4+π chain with its calibration took 182–196 s; the bench took 2 303 s.

## Open, for the author (nothing chosen)

- (a) **p = 6 · W-tous W1 as the cascade's setting** (2 reasons left: A_dir, F), or p = 6 · W-phys (the same 2, with
  straighter trunks).
- (b) **The trunks' remaining alignment**: a stronger bend than the brief's bound, or a fine-level routing that does
  not re-cut D8 segments.
- (c) **The depression clause**: count new depressions rather than the net change.
- (d) **The flat interfluves**: half is the hypsometry (the tectonics' wide high lands, queued), a third the D5 lakes'
  margins (the lakes' work). Neither is the cascade's to fix.

## At commit (2026-10-10)
- **(a) p = 6 · W-tous W1 becomes the cascade's default, provisionally** (with F163's roughness and F162's R4+π).
  Production does not change.
- **(b) The remaining alignment** (36.4 % against 30.9 %) is queued, untreated for now: the network will change with
  the next work.
- **(c) From F165, the depression clause counts NEW depressions**, matched by footprint before and after the step. The
  net count stays in the reports, for information.
- **(d) The flat interfluves (F) go to the tectonics** (F165) **and to the lakes' work.** For F, Corsica stays the
  reference only within the mountain class.
- **The cascade is frozen** at p6 · W-tous W1, with F164's budgets and k. Nothing downstream is recalibrated in F165.

**The author's answers for F165, verbatim**:
1. Type de continent : « Paramètre choisi, mais normalement on a déjà ça. C'est juste qu'on soulève tout au même
   rythme. Souviens toi encore, on est closure-dépendant. Toute la tectonique est basée dessus. Donc il faut un
   mécanisme compatible. Le relief est fait à 64² on devrait avoir toutes ces zones déjà formées. »
2. Échelle : « Techniquement, je veux un équivalent 7.5x la taille. Par contre pour des raisons de dimension, on ne
   peut pas tout scaler et on est obligé de rester à une proportion de type Corse afin d'éviter des pentes trop
   importantes. Mais en soit, les rivières doivent être larges (les plus grosses avec +100-150m rive à rive à
   l'embouchure) ce qui représente 3 à 4 cellule hex de Living Landz. Des chaines de montagnes, des plateaux, etc,
   mais dont la répartition est celle d'un continent, pas de la Corse. C'est dans le but de la jouabilité. On veut
   un continent varié, avec des zones de montagne, des plateaux et plaines permettant de placer les villes. Les
   zones montagneuses pour la récupération de matériaux riches et de valeur. »
3. « La cascade en l'état est bien. Allons donc sur la physio du continent. »
4. « On veut quelque chose qui soit physique pas construit juste pour faire apparaitre tel ou tel phénomène. On s'est
   rendu compte que ce genre de construction va toujours faire plein d'itération de calage et sans réellement
   aboutir car on se retrouve à devoir faire des compromis (facettes, murs, peigne, etc). »
5. Rectification : les hex de Living Landz font 40 à 50 m de **diamètre** en flat-top, pas de rayon. 100 à 150 m de
   rive à rive font donc bien 3 à 4 hex, et une cellule de 8 192² (48,8 m) correspond à environ un hex.
6. « La corse sert aussi pour la granularité du littoral (l'aspect naturel). »
