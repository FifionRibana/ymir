# F161 — declared before any measurement (2026-10-09)

**The third round of the cascade.**
- **A real reference**: the relief of Corsica (Copernicus DEM GLO-30).
- **An amplification that carves**: n = 1, k calibrated per level on Corsica, the area capped.
- **A smooth retargeting** that separates « creuser les vallées » from « abaisser le continent ».

**Production is unchanged**: the chain is gated, and the guard 6 / 6 runs before and after. One production file gains
an inert field: `StreamPowerConfig::area_cap_cells`, `None` everywhere, skipped in serialisation, so the cache keys
are unchanged and the F159 hash pin still holds.

Nothing below has been run on N1 or measured on Corsica. The predictions are in `f161_predictions.md`.

## C — the Corsica reference

1. **The data**:
   - Copernicus DEM GLO-30, five tiles (N41–N43, E008–E009), from `copernicus-dem-30m` (AWS open data), 2026-10-09;
   - the licence and the attribution are in `docs/refs/REFERENCES.md`;
   - the raw tiles stay outside the repository; the grids live in `data/corsica/` (git-ignored).
2. **The preparation** (`prep_corse.py`):
   - each 1″ pixel centre is projected to a local metric frame (lat 42.18°, lon 9.045°; WGS84 metres per degree at the
     pixel's latitude);
   - pixels are averaged into a 200 km square of 48.83 m cells (= 400 km / 8 192), rows south-first;
   - **the sea**: a cell with no land pixel (DEM ≤ 0) is −50 m;
   - **only Corsica's main land mass is kept**: 290 islets (Capraia, La Maddalena…) set to sea, 13 820 cells;
   - block means give 2 048 … 128 cells: cells of 97.7 / 195.3 / 390.6 / 781.3 / 1 562.5 m.
   - **Cascade level N (400 km) ↔ Corsica grid N / 2 (200 km)**: the same cell.
3. **The instruments, on Corsica and on the cascade alike**: F160's `read()`, in metres. The ocean is the cells ≤ 0
   4-connected to the border.
   - **The spectrum in km** (amended before any measurement: the 2D FFT replaces F160's box-blur bands):
     - F160's bands (differences of box blurs) overlap. Their sinc lobes move the power between bands: checked
       analytically, a 3-cell sine reads mostly in F160's band 1. They cannot calibrate on « the 2–4 cell band ».
     - F161 uses the **2D periodogram** (radix-2 FFT; the sea at 0 m; no taper, since every grid is bordered by sea
       and so periodic at 0).
     - It is summed by **octave of wavelength**: octave j = [2^(j+1), 2^(j+2)) cells, octave 0 = 2–4 cells. Each is an
       **RMS in metres per textured land cell**: the variance (Parseval) divided by the share of land minus D5.
     - The D5 lakes stay in the field. They are flat and contribute little, and the normalisation excludes them.
     - Reported in km up to about 25 km. The method follows **Perron, Kirchner & Dietrich 2008, JGR** (« spectral
       signatures », not fetched), without their detrending and windowing.
     - A permanent test checks that a 3-cell and a 12-cell sine land in their octave only, that the octaves sum by
       Parseval, and a flat field (negative control).
   - **The slopes**: p50, p90, and the shares above 28° and 33° (central differences, as F157-B7).
   - **The land hypsometry**: the deciles.
   - **The drainage density**: D8 channel cells (≥ max(1 km², 4 cells), the window's river threshold) × the cell size,
     per km² of land.
   - **The valley spacing λ** in km (transects, prominence ≥ 1 m).
   - **F160's**: crest facets (with the control), walls at 28°, R8, the coherent windows (the comb), and the local relief
     (12.5 km blocks).
4. **The scale reservation**: Corsica is a texture reference up to about 20 km. It is not a reference for the large
   forms, the network's size, the hypsometry's absolute values, or the coast.

## D — the design

1. **The physics level**: 256², F160's P256 (U₀ calibrated on the peak, 2 850 m; D5 lift).
2. **The amplification**, at 512², 1 024² and 2 048²:
   - bicubic ×2;
   - then **3 erosion blocks of 2 iterations**, each block followed by the smooth retargeting;
   - then the talus, then the deposition.
   - **No uplift.**
   - **Erosion N1**: the existing implicit solver `incise_with_floor` (Braun–Willett, n = 1, m = 0.5; relief-v3 at the
     level's cell: MFD p = 2 for the area, the D8 stack on the filled surface, the 0.5 m base level, the lateral
     erosion K_lat = 4, A_c = 0.1 km²), with the talus and the diffusion off inside it.
     - k_L as calibrated (D3), dt = 1.
     - The ocean is held at sea level and the land floor is 0.5 m, as in F160.
   - **The area cap**: min(A, a_max) in the incision term only (the A_c gate and the lateral erosion read the real
     area).
     - **a_max = 250 cells** of the level (Schott's scale selector): 610 / 153 / 38 / 9.5 km² at 256² … 2 048².
   - **The smooth retargeting**, after each erosion block:
     - the error e = R(z_fine) − z_prev, on the working fields (the ocean at 0 m on both);
     - e is upscaled **bicubically** and subtracted from the land cells, with the land floor at 0.5 m.
     - **Never a constant per block.**
     - R is the **centred full weighting** (1/4, 1/2, 1/4). On our cell-corner mapping the brief's plain 2×2 block mean
       is off by half a cell (F159-D6); the full weighting is the centred mean at the 2×2 scale.
     - **The coast**: the sea is the 0 m base level on both fields. The error therefore carries only land change;
       sea cells are never corrected.
     - **What it keeps**: only the wavelengths the previous level cannot carry (< 4 fine cells). Everything coarser
       that the erosion changed is put back.
   - **The talus**: F160's form, h ← h + m(α − β), m = 0.002 cell.
     - **Its angle adapted to the cell**: s₀ = tan 33° · (c / 30 m)^(H − 1) · (0.8 + 0.6 ν), with H = 0.8 and ν F160's
       noise. **PROXY.** A repose slope measured at 30 m, carried to a cell c by a self-affine scaling.
     - So s₀ is 16–29°, 18–31°, 20–35° and 22–38° at 781 / 391 / 195 / 98 m. **Not taken from Corsica.**
   - **The deposition**: F160's (its own k rule, n = 2, m = 0.8, k_c = k_d = 0.1).
   - **The budgets**:
     - T: F160's 1 000 / 2 000 / 6 000 at 512² / 1 024² / 2 048²;
     - D: 700 / 200 / 150.
3. **The calibration of k_L, per level**:
   - **the target**: the spectrum's octave 0 (wavelengths 2–4 cells) of the level equals Corsica's octave 0 at the same
     cell size, within ±20 %;
   - **measured on the whole level's output**, after T and D;
   - **the method**: a trial at **k₀ = 100** (code units, k·dt per iteration), ×10 expansions until bracketed (at most
     6), then bisection in log k. At most 8 trials. The trials are reported.
     - **Amended before any measurement**: the first k₀ = 1e-3 was 5 orders too small (F159's per-step k·dt was 200).
       The permanent N1 test caught it: the solver carved nothing.
     - If the target is not reached (not monotonic, saturated), the closest trial is kept and this is said.
   - **This single measure serves the calibration and does not enter the judgment.**
4. **What judges** (never used for calibrating):
   - the slope distribution;
   - the valley spacing (km);
   - the drainage density;
   - the hypsometry;
   - the facets and the walls;
   - R8 and the comb;
   - **the other bands** of the spectrum.
5. **The variants**, each a full chain from P256 to 2 048²:
   - **(N1)**: as above.
   - **(N1-sans plafond)**: no area cap, at N1's calibrated k_L.
   - **(N1-sans recalage)**: no retargeting, at N1's k_L.
6. **The talus at a quarter**: N1 again with T / 4 at 1 024² and 2 048², at N1's k_L. Its slopes and time are against
   N1's.
7. **The D5 lakes**:
   - **the D5 cells**: the physics level's lifted cells, carried to each level by nearest neighbour, plus each level's
     newly lifted ones;
   - **excluded from every texture instrument** (bands, slopes, facets, walls, R8, comb, λ, density);
   - **the lakes containing a D5 cell are counted apart.** Their treatment goes to the lakes' work.
8. **The physics level's facets**: no fix this round. They are measured at each level.

## V — the Cascade window

- **The variants** N1 / N1 sans plafond / N1 sans recalage, beside F160's S / S+ρ / S+U. **N1 is the default.**
- **« k du niveau »**: the calibrated value displayed and editable; a button calibrates it on Corsica (D3) when the
  data is there.
- **A « Référence Corse » view**: Corsica at the displayed level's cell size, beside the cascade, with the same shading,
  colours and pixels per cell.
- **The rest as in F160.**
- The production world untouched: the guard 6 / 6 and the worker test.

## R — the stop rule (the brief's, verbatim; against Corsica)

« Pour N1, on s'arrête au niveau où l'une de ces conditions est vraie :
- facettes de crête ou murs à 28° au-delà de la Corse à la même taille de cellule de plus de 50 % (relatif), **et**
  plus de 20 cellules ;
- distribution des pentes : p50 ou p90 hors d'un facteur 1,5 de la Corse ;
- espacement des vallées hors d'un facteur 1,5 de la Corse ;
- dérive du relief d'un niveau à l'autre de plus de 2 % ;
- point culminant hors de la cible de plus de 15 % ;
- troncs de plus de 100 km² : p90 de déplacement au-delà de 1,5 cellule. »

The cost stays an alert, not a stop.

**How it reads, declared**:
- **the facets**: the cascade's share > 1.5 × Corsica's, with > 20 facet crests. The same for the walls (> 20 cells).
  The facets count only where their control passes.
- **the slopes**: p50 or p90 outside [Corsica / 1.5, Corsica × 1.5].
- **λ**: outside [Corsica / 1.5, Corsica × 1.5].
- **the drift**: |restricted bias| > 2 % of level n's mean land altitude.
- **the peak**: outside [2 295, 3 450] m (F160's band), the level itself (no retargeting).
- **the trunks**: n+1 → n, p90 > 1.5 cells of level n.
- **The physics level (256²) is measured against Corsica too**, but R applies from 512².

## I — the images

- F160's three crops, plus a whole view per level, north up: **N1 beside Corsica at the same detail**, with rivers and
  lakes;
- the spectra per band in km (cascade, Corsica, témoin) as a figure;
- the hypsometry per level;
- one crop N1 against N1-sans recalage.
- Every figure that shows Corsica carries the Copernicus notice.
