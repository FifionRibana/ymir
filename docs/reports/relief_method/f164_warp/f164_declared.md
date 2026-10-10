# F164 — declared before any measurement (2026-10-09)

**The sixth round of the cascade.** It bends the trunks inherited from the coarse grid with a smooth, invertible warp
of the upscaled field, arbitrates the physics level's MFD exponent p, and attributes the flat interfluves.

**Production is unchanged**: the chain is gated, and the guard 6 / 6 runs before and after. Nothing below has been
measured. The predictions are in `f164_predictions.md`.

**Not blind**:
- F163's readings: the physics-level variants at 256², and the R4+π chains from P0 and from roughness + p = 6 at
  2 048²;
- F163's images: the straight trunks and the flank comb at p = 6.

## 0 — unchanged from F163

- The scales: 400 km, real areas.
- Corsica at the same cell (`corse_{N/2}`).
- The instruments (F161–F163), with I2 read on S − 1 (confirmed at F163's commit).

## D — the design

### D1 — the warp (W)

At an amplification level of n cells, right after the bicubic ×2 of the previous level, the upscaled field u is
resampled:
- u_W(x) = u(x + d(x)), sampled with the Catmull-Rom bicubic (the kernel of `upsample2`), periodic like it.
- **d, a smooth displacement field**: each component is A · (2f − 1) / √2, so |d| ≤ A.
  - f is a one-octave gradient noise (`SeededNoise`), min–max to [0, 1] over the grid, of wavelength L.
  - The seeds are `WorldSeed(PSEED).derive_seed("cascade_warp_x_<n>")` and `"cascade_warp_y_<n>"`.
- **A and L are counted in cells of the previous level** (n / 2). **The two declared settings**:

  | setting | A | L |
  |---|---|---|
  | **W1** | 0.5 | 4 cells |
  | **W2** | 0.5 | 2 cells |

  Both are at the brief's ceiling of A. W2 is twice as tortuous for the same displacement.
- **Invertibility**: J = det(I + ∇d) by central differences on the cell grid. Its minimum is recorded per warp.
  - **A permanent test** checks J > 0 everywhere for W1 and W2 at 512², 1 024² and 2 048², and **J ≤ 0 somewhere** for
    a negative control (A = 2, L = 2 previous cells).
  - **If a declared setting folds (min J ≤ 0) on the témoin, it is not run**; that is reported, with no substitute.
- **Sea and land are warped together**: the ocean mask, the D5 lift and the land floor are computed on u_W. The D5 mask
  carried from the previous level is warped by nearest-cell sampling at x + d(x). The coast's change (land cells,
  coastline cells) is reported, not judged.
- **The retargeting is done on the warped field**:
  - the level's reference (`prev_work`, which the smooth retargeting pulls toward) becomes R(u_W), the full-weighting
    restriction of the warped field, the ocean at 0 m;
  - for R4, the n − 2 reference becomes R(R(u_W)) on a warped level;
  - so the retargeting keeps the warp's geometry instead of pulling the field back to the unwarped one.
- **The warp is a PROXY** for the absent geological heterogeneity (structure, lithology that bends courses). It has no
  physical cause in the model.
- **The placements**:
  - **(W-phys)**: one warp only, on the physics level's output, applied at the first amplification level (512²), at A
    and L in physics-level cells;
  - **(W-tous)**: a warp at each amplification level (512², 1 024², 2 048²), A and L in the previous level's cells, so
    the displacement halves in metres at each level.

### D2 — the warp's controls (per warp, on the upscaled field before and after)

- **The closed depressions**: `cascade::hydro`'s lakes (filled − z > 0.5 m, ≥ 2 cells, 8-connected), counted before
  and after.
- **The hypsometry**: the land deciles d1–d9 before and after. The change is max_d |Δq_d| / q_9 (before). The land
  cell count is reported too.
- **The Jacobian**: min J, and the share of cells with J < 0.5.

### D3 — the physics level's exponent

- p ∈ {2, 4, 6} at 256².
- Each with F163's roughness (0.5 × Corsica's 3.1–6.25 km octave at 1 563 m, from the bench's reading).
- U₀ is calibrated on 2 850 m.

### D4 — the chains (all R4+π to 2 048², k calibrated on the peak per level, F162's method)

**The reduction, declared for the cost.** The full grid {2, 4, 6} × {none, W-phys, W-tous} × {W1, W2} is 15 chains,
~7 min each with the calibration. Instead:
1. **The settings are compared at p = 2 only**: W-phys and W-tous, each with W1 and W2 (4 chains), plus p = 2 without W.
2. **The better setting is carried to p = 4 and p = 6**, with none, W-phys and W-tous (6 chains).
   - It is the one with the fewest R reasons at 2 048², summed over its two chains (W-phys and W-tous); ties go to the
     lower mean A_dir.
   - This precision was written before any measurement.
3. **11 chains in all.** The setting not carried is set aside at p = 4 and 6, because the settings are judged where the
   physics level is F163's P0 geometry.

**The best chain**: the fewest R reasons at 2 048², ties going to the smaller Σ |ln ratio| over R's measures.

## F — the flat interfluves attributed (measured, nothing corrected)

- **On the best chain at 2 048² and on Corsica at 195 m**, with F's flat cells (I3: slope < 0.05, > 3 cells from an
  order-≥ 2 river) among the cells read.
- **By altitude band**: 0–50, 50–100, 100–300, > 300 m. Each band reports:
  - its share of the land read;
  - its share of the flat cells;
  - its own flat rate.
- **By distance to the coast** (Chebyshev, from the ocean cells): 0–2, 2–5, 5–10, 10–20, > 20 km, the same three
  numbers.
- **Near the D5 lakes**: the flat rate of the cells read within 2 km of a D5 cell, against the rest. Corsica has no
  D5 cells.
- **The question**: are the flat cells mostly low, so the hypsometry and the tectonics, or spread out?
  - **« Mostly low » is declared as ≥ 60 % of the flat cells below 100 m AND a flat rate below 100 m ≥ 2× the rate
    above 300 m.**
  - Otherwise they are « spread out ».

## R — the stop rule (the brief's, judged at 2 048² against Corsica at 195 m)

**A chain is rejected** if any holds:
- F163's rule:
  - the facets or the 28° walls beyond 1.5 × Corsica and over 20 cells;
  - the slope p50 or p90 outside ×1.5;
  - λ outside ×1.5;
  - an octave between 0.4 and 12.5 km outside ×1.5;
  - CV_λ, CV_L, σ_θ or R2_g outside ×1.5;
  - π over 10 %;
  - A_dir outside ×1.5;
  - S_W − 1 at 2, 5 or 10 km outside ×1.5;
  - F outside ×1.5;
- **plus: a warp creates closed depressions beyond 1 % of the count before it** (after > 1.01 × before), at any level;
- **and (declared here) a warp with min J ≤ 0**.

**At each level** (unchanged):
- the drift ≤ 2 %;
- the peak within [2 295, 3 450] m;
- the trunks ≥ 100 km² p90 ≤ 1.5 cells below 2 048².

## V — the Cascade window

- **« Déformation »**: none / W-phys / W-tous, with A and L (previous-level cells).
- **The physics level's p**: F163's MFD row, unchanged.
- The roughness is F163's.

## I — the images

At 2 048², north up, beside Corsica, the rivers drawn:
- **the best chain without W, with W-phys and with W-tous**;
- **for the best placement, p = 2, 4 and 6**;
- on the 25 km tight crop and crops A, B, C.

**For the author's four questions**:
1. Les grands cours d'eau sont-ils encore des lignes droites ?
2. Le peigne sur les flancs : présent, atténué, absent ?
3. Est-ce que la déformation se voit comme telle (des formes « tordues » artificiellement) ?
4. À côté de la Corse, est-ce que ça « fait vrai » ?
