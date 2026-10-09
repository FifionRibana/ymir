# F163 — declared before any measurement (2026-10-09)

**The fifth round of the cascade.** It redoes the physics level (256², 1.5 km cells), so its network is less aligned
on the grid and more dissected, and names what the author sees: the directions, the sinuosity and the flat
interfluves.

**Production is unchanged**: the chain is gated, and the guard 6 / 6 runs before and after. Nothing below has been
measured. The predictions are in `f163_predictions.md`.

**Not blind**: F161's reading of the physics level at 256² (P0) is known: λ 17.64 km against Corsica's 9.08 km,
facets 13.0 % against 3.45 %, and the 3.1–6.25 / 6.25–12.5 km octaves 32.6 / 121.4 m against 84.9 / 152.4 m. So are
F162's tight crop (straight trunks, a fishbone, flat interfluves) and the code facts in D0.

## 0 — the scales (F162's, unchanged)

- 400 km; the relief is at the real scale; every area is real (cells × (400 km / N)²).
- The ×7.5 is for the exported hydrology only. `the_cascade_reads_no_signified_quantity` stays.
- The cascade at N cells is compared with `corse_{N/2}`, the same cell:

  | cascade | Corsica | cell |
  |---|---|---|
  | 256² (the physics level) | `corse_128` | 1 563 m |
  | 1 024² | `corse_512` | 391 m |
  | 2 048² | `corse_1024` | 195 m |

## I — the three new instruments (each validated on its controls before any measurement)

**The network** is `cascade::hydro`, F162's: D8 on the filled surface, real km², the threshold max(1 km², 4 cells),
and the Strahler order on it. **The cells read** are the land minus the D5 cells minus a 2-cell border.

**I1 — the directions, A_dir**:
- The rivers of order ≥ 2. An order-ω stream starts at a cell of order ω with no donor of order ω, and follows the D8
  receivers while the order stays ω (F162's CV_L streams).
- Each stream is cut, from its start, into consecutive pieces of 8 D8 steps. A last piece of 4–7 steps is kept; a
  shorter one is dropped. A piece touching an excluded cell is dropped.
- Each piece's chord angle φ, folded modulo 45°. **It is aligned when min(φ, 45° − φ) ≤ 5°**, i.e. within ±5° of 0°,
  45°, 90° or 135°.
- **A_dir = the chord length of the aligned pieces / the chord length of all the pieces.** Isotropic directions read
  10 / 45 = 22.2 %.
- Reported, not in R: the same share on the trunks only (pieces whose first cell drains ≥ 100 km²), the author's
  « grands cours d'eau ».
- **The « high » control: the D8 network of a smooth inclined plane.**
  - The plane falls 2 % toward the west, its gradient at **10° from the x axis**. D8 runs its lines due west
    (cos 10° > cos 35°), while the true gradient is at 10°: the D8 artefact.
  - A pure plane gives lines that never meet, so no order 2. A white hash noise of 0.2 × the per-cell drop is
    therefore added, so a few steps go diagonal and the lines merge.
  - The sea is the 2 western columns and the 2 southern rows.
  - **Set aside before any run: a smooth funnel.** D8 on a cone converges onto the 22.5° bisectors, zig-zag lines
    whose chords read unaligned. For the same reason a smooth trough at a non-grid angle reads its true angle: the
    instrument reads the course's direction over 8 steps, not the D8 steps.
- **The « low » reference**: Corsica at the same cell. F162's fractal control is reported.
- **It passes when the funnel reads ≥ 1.5 × Corsica**, at each cell where it is used.

**I2 — the trunk sinuosity, S_W**:
- The network cells draining ≥ 10 km² (real). From each one, walk down the D8 receivers, adding the path length
  (1 or √2 cells, in km), until it reaches W km.
- S = the path length / the chord between the first and the last cell. A walk that reaches the sea, an outlet or an
  excluded cell first gives no window.
- **S_W = the mean over the windows**, for W = 2, 5 and 10 km. **A window shorter than 4 cells is not read**: at
  1 563 m only W = 10 km is read; at 391 m and 195 m all three are.
- **Amended before any measurement — R reads the excess S_W − 1, not S_W.** S is bounded below by 1, so the brief's
  ×1.5 on S itself cannot reject a straight line (1.0) against a sinuous Corsica (≤ 1.5). A ratio on S − 1 can.
  - This is a reading of the brief's rule, declared here for the author to confirm or refuse.
  - S_W itself is reported beside it.
- **The controls**:
  - **straight**: a V valley (5 % down the axis, 10 % across) along a straight line, read at 0°, 22.5° and 30° to
    the grid. At 0° it should read 1.00. At 22.5° the D8 staircase sets a floor of ~1.08: the D8 floor of the
    instrument, reported;
  - **sinuous**: the same valley along y = a · sin(2πx / P) (P = 32 cells, a = 4 cells), against S_W of the same
    windows walked on the continuous curve, as an accuracy check (reported);
  - **Corsica**.
- **It passes when the straight valley at 0° and at 22.5° reads S − 1 ≤ (Corsica's S − 1) / 1.5**, per window read.

**I3 — the flat interfluves, F**:
- **x = 3 cells; y = 0.05 m/m (2.9°)**.
- F = the share of the textured land with
  - its slope (F161's central differences) < y, **and**
  - its 8-connected (Chebyshev) distance to the nearest order-≥ 2 river > x.
- The same x and y at every cell. The comparison is always with Corsica at the same cell.
- **The « high » control**: a terraced field.
  - A plane (3 % down y) quantised into 100 m treads, so the treads are flat and the risers are steps.
  - V valleys (60 m under the local tread, 8 cells wide) every 64 cells carry the order-≥ 2 rivers to a sea row.
  - **Amended at its control, before any measurement on Corsica or the cascade**: each tread tilts 1 % toward its
    nearest valley (still flat: < y).
    - As first built (perfectly flat treads, the valley floor falling half the plane), it read F = 18 %. The valley
      had widened to ~50 cells, and the flats' routing (Garbrecht–Martz) made order-2 streams on the treads.
    - With the floor fixed at 60 m under the tread, F was 38 %, the treads still carrying order-2 streams.
    - Tilted, the treads drain in parallel order-1 lines, as a real terrace does.
- **It passes when the terraced field reads ≥ 1.5 × Corsica.** F162's fractal control is reported.

**Where the instruments are validated**: at 256² (1 563 m) and at 512² (195 m, F162's control grid). The 391 m reading
uses the 195 m verdict, the same instrument in cells.

**Measured on**:
- Corsica at 1 563 m, 391 m and 195 m;
- F161's N1 (as committed, k fixed);
- F162's R4+π (F162's chain, k calibrated on the peak);
- at 1 563 m that is the physics level P0, shared by both chains; at 391 m (1 024²) and 195 m (2 048²) each chain.

## D — the physics-level variants (256²), each measured at that level against Corsica at 1 563 m

**D0 — the code facts, read before any run** (`CascadeConfig::solver`, `relief_v3`):
- **The level already routes its AREA by MFD** (`mfd_exponent` p = 2). Its RECEIVERS are D8 (`compute_flow`), and
  the implicit solver needs one receiver per cell.
- **The head threshold is sub-cell**: A_c = 0.1 km² = 0.04 cell at 256², so it never gates.
- **So is the diffusion**: D = 0.16 · K_SI · cell², L_c = √0.16 = 0.4 cell.
- So neither the head threshold nor the diffusion can densify the network at this level; raising either can only
  coarsen it. The lever left for P-dissection is how the area concentrates, the MFD exponent.

**The variants**:
1. **(P0)** the current level: F160's `physics_level`, U₀ calibrated on the peak (2 850 m), at most 300 steps.
2. **(P-mfd), as a weighted random receiver** (the brief's fallback: the implicit solver needs one receiver). At each
   step, on the filled surface:
   - the candidates are the in-map neighbours strictly lower than the cell, with a drop slope ≥ τ · s_max,
     **τ = 0.5**;
   - one is drawn with probability ∝ its slope;
   - **it is drawn anew at every step**, from a hash of (the seed, the step, the cell). The seed is
     `WorldSeed(PSEED).derive_seed("cascade_phys_receiver")`;
   - a cell with no strictly lower in-map neighbour (a flat, a sill) keeps the D8 receiver, and so does a cell whose
     D8 receiver leaves the map. Since every draw is strictly lower, no cycle is possible;
   - the MFD area is unchanged.
   - On a plane at 22.5° it draws E and NE half and half; on a plane at 0°, E with 41 % and each forward diagonal with
     29 %. The mean direction follows the gradient instead of the nearest D8 direction.
   - `StreamPowerConfig::random_receiver`, `None` in production, skipped in serialisation; the hash-pinned
     `None` path is unchanged.
3. **(P-bruit)** an initial roughness (not an uplift noise: an uplift noise would stay in the steady state as relief
   with no cause, whereas a roughness only seeds the network, which is the reviewer's hypothesis (2), the smooth
   start):
   - band noise of base wavelength 8 cells, 2 octaves (8 and 4 cells, i.e. the 6.25–12.5 and 3.1–6.25 km octaves),
     persistence 0.5, normalised to unit RMS over the land;
   - amplitude **A = 0.5 × Corsica's 3.1–6.25 km octave at 1 563 m** (taken from the bench's own reading of
     `corse_128`);
   - source `"cascade_phys_rough"`;
   - added to the land of the upsampled field before both the trial run and the calibrated run, never below
     sea + 1 m, the D5 cells untouched.
   - **Its share**, both reported:
     - (a) A² / the variance of the final field's 3.1–12.5 km band (octaves 0 and 1 at 256²);
     - (b) var(z_bruit − z_P0) / var(z_bruit) over the textured land.
4. **(P-dissection)** the MFD exponent:
   - scanned at p ∈ {1, 1.5, 3, 4, 6} and D8 (no MFD), each with its own U₀ calibration;
   - **the kept setting is the one whose λ is nearest 9.08 km in |ln|**. If it is P0's p = 2, P-dissection is « not
     reached » and leaves the combinations;
   - if the nearest is not within ×1.25 of 9.08 km, that is said, and the setting is still kept.
5. **The combinations**: P-mfd + P-bruit, P-mfd + P-dissection, P-bruit + P-dissection, all three.

**Measured at 256², against Corsica at 1 563 m**:
- λ, the facets (with F161's pyramid / paraboloid control), the 3.1–6.25 and 6.25–12.5 km octaves;
- A_dir, S_10 − 1, F;
- the peak (2 700–3 000 m after the calibration); a variant outside it is excluded;
- also reported: the slope p50 / p90, the drainage density, the trunks' share in A_dir, the steps and the seconds.

**The best variant**:
- the fewest of those seven measures outside ×1.5 of Corsica, among the measures whose control passed;
- ties go to the smaller Σ |ln ratio|.
- **The best and P0 then run through F162's R4+π chain** (the talus at ¼ from 1 024², k calibrated on the peak per
  level, R4 retargeting, π) to 2 048².

## R — the stop rule (the brief's, judged at 2 048² against Corsica at 195 m)

**A variant is rejected** if any of these holds:
- F162's rule:
  - the facets or the 28° walls beyond 1.5 × Corsica and over 20 cells;
  - the slope p50 or p90 outside ×1.5;
  - λ outside ×1.5;
  - an octave between 0.4 and 12.5 km outside ×1.5;
  - a predictability instrument that passed its control (CV_λ, CV_L, σ_θ, R2_g) outside ×1.5;
  - π over 10 % of the final variance;
- **plus the three new instruments outside ×1.5 of Corsica, those whose control passed**: A_dir; S_W − 1 at 2, 5 and
  10 km; F.

**At each level** (unchanged):
- the drift beyond 2 %;
- the peak outside [2 295, 3 450] m (the target ±15 %);
- the trunks ≥ 100 km² beyond 1.5 cells (p90), except at 2 048², which goes to the lakes' work.

## V — the Cascade window

The physics level's options can be ticked before « Physique » runs:
- « récepteur aléatoire » (τ);
- « rugosité initiale » (its factor × Corsica's octave);
- « exposant MFD p ».

## I — the images

- At 2 048², north up, beside Corsica at the same detail, **P0 and the best variant with the rivers drawn**:
  - the 25 km tight crop (F162's);
  - crops A, B and C.
- **The physics level alone (256²) for each variant**, the whole map with the rivers, beside Corsica at 1 563 m.
- For the author's four questions:
  1. Les grands cours d'eau sont-ils encore des lignes droites ?
  2. Les arêtes de poisson et les interfluves plats ont-ils reculé ?
  3. Le relief entre les vallées (quelques km) est-il assez marqué ?
  4. À côté de la Corse, est-ce que ça « fait vrai » ?
