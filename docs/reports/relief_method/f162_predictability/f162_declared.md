# F162 — declared before any measurement (2026-10-09)

**The fourth round of the cascade.** It measures the « predictability » (the author: « trop prévisibles … un schéma se
répéter ») against Corsica, tries to remedy it with a noisy hardness and an initial perturbation, and fills the 1.5–6
km deficit with a 4×4 retargeting.

**Production is unchanged**: the chain is gated, and the guard 6 / 6 runs before and after. Nothing below has been
measured. The predictions are in `f162_predictions.md`.

## 0 — the scales and the units, declared

1. **The relief is at the real scale**: 400 km, 48.8 m per cell at 8 192² (F67's decision; the terrain is not
   scaled).
   - The ×7.5 of `geo_scale_ratio` applies only to the exported hydrological quantities (discharges, widths).
   - Every comparison with Corsica is in real km, at the same km per cell.
2. **Every area threshold of the cascade is real** (never the ×56.25 signified area):

   | threshold | unit |
   |---|---|
   | trunks ≥ 100 km² | real km² (cells × (400 km / N)²) |
   | the network ≥ max(1 km², 4 cells) | real |
   | the lakes' area | real |
   | the Strahler order | on that real network |
   | a_max = 250 | cells of the level |
   | A_c = 0.1 km² | real, F1 |

   - F147's lesson (831 junctions announced, 465 real) came from mixing the two units.
   - **A permanent test checks that the cascade reads no signified quantity**:
     - its sources never name `geo_scale_ratio`;
     - `hydrology` gives cells × cell² exactly.
3. **The final level is 8 192² at 48.8 m.** Corsica is measured at 49 m for information (`corse_4096`), for the day the
   cascade gets there. **F162 is judged at 2 048²** (195 m).

## P — the predictability instruments (each validated on two controls first)

**The reviewer's hypotheses to check**: the regularity comes from
- (1) a uniform k;
- (2) the smooth surface each level starts from;
- (3) one valley generation of a single size per level.

A regular spacing is physical (Perron 2009); **what is measured is the dispersion.**

**Where they apply**:
- on the land minus the D5 cells (F161-D7) and minus the cells within 2 of the map border;
- the network is `cascade::hydro` (D8 on the filled surface, real km², the threshold max(1 km², 4 cells)).

**The six instruments**:
1. **CV_λ, the dispersion of the valley spacing**:
   - along every row and column, in land runs ≥ 16 cells, the valley minima (prominence ≥ 1 m, F159);
   - the gaps between successive minima, in km;
   - CV = std / mean of the gaps.
   - **Amended at its control, before any measurement on the cascade or Corsica**: the CV is computed **along the rows
     and along the columns apart**, and the two are averaged. Mixing the two directions read the regular fishbone as
     bimodal: trunks every 64 cells across, gullies every 8 along. Its CV was 1.23, against the fractal's 0.61.
   - Declared: the transects are along the axes, not perpendicular to the crests; a comb at an angle θ is read at
     1 / cos θ, so its CV is unchanged.
2. **CV_L, the tributaries' lengths by Strahler order**:
   - a stream of order ω is the D8 path from a cell with no donor of order ω, down to where the order changes or the
     sea;
   - its length in km (√2 on diagonals);
   - the CV of the lengths per order, ω = 1, 2, 3; **the instrument is their mean.**
3. **Horton's ratios, per basin**:
   - the basins are the D8 outlets with a maximum order ≥ 3;
   - in each, R_b = the geometric mean of N_ω / N_(ω+1), and R_l = the geometric mean of L̄_(ω+1) / L̄_ω;
   - **the instruments are their CVs across basins**, CV_Rb and CV_Rl; the means are reported.
4. **σ_θ, the confluence angles**:
   - at each junction (≥ 2 channel donors), the two largest donors;
   - each one's direction is from the cell 3 steps up its own branch (following its largest donor) to the junction;
   - θ is the angle between the two incoming vectors (degrees). **The instrument is its standard deviation**; the
     mean is reported.
5. **P_s, the spectrum's peakiness**:
   - the radial periodogram (F161's FFT) in ¼-octave bins, between 2 cells and 12.5 km;
   - log₁₀ power against log₁₀ wavelength, a least-squares line;
   - **P_s = the residuals' standard deviation, in dex.** A cascade that puts its power in one narrow peak per level
     reads high; a continuous curve reads low.
6. **R2_g, the gullies' parallelism**:
   - the order-1 streams at least 3 cells long, their orientation the chord from head to end;
   - non-overlapping 5 km windows holding ≥ 5 gully heads;
   - **R2_g = the mean over windows of |⟨e^(2iθ)⟩|**: 1 = all parallel. This is the axis concentration of F84 / F83,
     transposed to the gullies.

**The validation** (rule 13), on 512² grids of 195 m cells:
- **the « regular » control**: a periodic fishbone. Trunks every 64 cells running to a sea row, gullies perpendicular
  every 8 cells; depths 30 m; a hash noise of 0.1 m breaks the flat ties.
- **the « dispersed » control**: a 9-octave fractal (persistence 0.6, 1 500 m) on a tilt toward a sea row.
- **An instrument passes when the regular control reads more regular than the dispersed one by ≥ ×1.5**: CV, σ_θ and
  P_s smaller; R2_g larger.
- **An instrument that fails its control is declared blind** and leaves the stop rule.
- **Corsica** is the second « dispersed » reference, reported (not a pass criterion).
- Measured on: Corsica, F161's N1 and the témoin, at 391 m and 195 m. Corsica also at 49 m, for information.

## D — the design

1. **The base, N1b**: F161's N1, with:
   - the physics level at 256² (n = 1, the area capped at 250 cells, the 2×2 smooth retargeting after each block of 2
     iterations, 3 blocks);
   - **the talus at a quarter** from 1 024² (decided at F161's commit);
   - **k calibrated on the octave's peak.**
     - The method, for a curve that is not monotonic: trials at log₁₀ k = 1.5, 2, 2.5, 3, 3.5.
     - **If the best is ≥ 0.8 × the target**: bisection in log k on the rising branch, between the last trial below
       the target and the first at or above it, until within ±20 % (at most 4 more trials). The smallest such k is
       kept.
     - **Otherwise**: two trials at the best k × 10^(±0.25); the highest octave is kept.
     - At most 9 trials. The trials are reported.
2. **The variants against the predictability**, each alone, then combined:
   - **(ρ)**: k(p) = k · (1 − ρ(p)) / (1 − 0.5), as Schott, with ρ = 0.5 + 0.4 · f.
     - f is a fractal noise, min–max to [−1, 1], from 25 km down to 2 cells of the level (persistence 0.5), source
       `WorldSeed` "cascade_hardness_n1", passed to the implicit solver as its `k_field`.
     - **The reserve** (the brief's): C-3 was rejected for boundaries made of noise. Here the noise modulates the
       erodibility and draws no boundary; the author's eye judges it.
   - **(π)**: an initial perturbation added to the upscaled land before the erosion, so the valley heads shift.
     - It is a band noise at the new octave: noise of wavelength 4 cells (2 octaves, persistence 0.5), normalised to
       unit RMS over the land.
     - Its amplitude is **0.15 × Corsica's octave 0 at the level** (8.0 / 4.6 / 2.3 m at 512² / 1 024² / 2 048²),
       source "cascade_pi".
     - **Its share** = Σ_L A_π,L² / the final field's variance in the octaves 0.4–12.5 km. The levels' perturbations
       are taken as uncorrelated, and the bicubic as variance-preserving. Reported.
   - **(ρ+π)**: both.
3. **(R4), the 1.5–6 km deficit**: the smooth retargeting against the level n − 2 instead of n − 1.
   - e = R(R(z)) − z_(n−2); z ← z − U₂(U₂(e)), smooth (bicubic), land only, the sea at 0 m on both fields. So each
     level can amplify its own octave and the one above (2–8 cells).
   - **z_(n−2)**: at 512² it is the physics level restricted to 128² (R(z_256)); at 1 024² the physics level; at
     2 048² the 512² level.
   - Measured alone, then **combined with the best predictability variant**: among ρ, π and ρ+π, the one with the
     fewest validated instruments outside ×1.5 of Corsica at 2 048², ties going to the smaller Σ |log ratio|. If none
     beats N1b, that is said.
   - **Every chain is calibrated on its own** (D1's method). **The calibration reads only the level's new octave**
     and does not enter the judgment.
4. **What judges**:
   - the P instruments that passed their control;
   - the full spectrum (octaves 0.4–12.5 km);
   - the slope distribution, λ, the facets, the walls, R8 and the comb (coherent windows).

## R — the stop rule (the brief's, verbatim; judged at 2 048²)

« **À 2048²**, contre la Corse à 195 m. Une variante est rejetée si l'une de ces conditions est vraie :
- facettes ou murs à 28° au-delà de 1,5 fois la Corse, et plus de 20 cellules ;
- pente p50 ou p90 hors d'un facteur 1,5 ;
- espacement des vallées hors d'un facteur 1,5 ;
- amplitude d'une octave entre 0,4 et 12,5 km hors d'un facteur 1,5 ;
- un des instruments de prévisibilité hors d'un facteur 1,5 de la Corse ;
- la perturbation (π) représente plus de 10 % de la variance du relief final.

**À chaque niveau** : dérive au-delà de 2 % ; point culminant hors de la cible de plus de 15 % ; troncs de plus de
100 km² au-delà de 1,5 cellule (p90), sauf à 2048², renvoyé au chantier lacs. »

**How it reads, declared**:
- the octaves: the four between 0.4 and 12.5 km at 195 m cells (0.39–0.78, 0.78–1.56, 1.56–3.13, 3.13–6.25 and
  6.25–12.5 km);
- the predictability instruments: only those whose control passed;
- the peak: [2 295, 3 450] m (F160's band);
- the drift: |restricted bias| > 2 % of the level n − 1's mean land altitude.

## V — the Cascade window

- **ρ, π and R4 can be ticked**, with their amplitudes adjustable (ρ's amplitude; π's factor × Corsica's octave), so
  the author can compare by eye beside « Référence Corse ».
- k's calibration (« Caler sur la Corse ») uses D1's peak method.
- The talus at a quarter by default from 1 024².

## I — the images

At 1 024² and 2 048², north up, beside Corsica at the same detail:
- N1b, the best predictability variant, R4, and their combination;
- crops A, B and C;
- **a tight crop (25 km) on a hillslope** to judge the comb.
