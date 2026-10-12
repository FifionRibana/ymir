# F173 — declared before any measurement (2026-10-12)

**A real sea-level curve, the attribution of the drowning, and the origin of the coastal rampart.**

**Production is unchanged**; the guard 6 / 6 is not touched. Frozen: T, Airy, the sources' drive, the torus, K2, 16 plates
with the fallback, F164's cascade.

**Not blind**: F172's readings, the code, and the sea-level series (its numbers read below; no world has run with it).

**Part 0**: F172 committed (`f6607c6`) with the reviewer's three errors, decisions (a)–(d) and the author's words.

**The panel**: the seeds 8, 18, 19, 23 and 29. The old témoin serves only the bit-for-bit negative controls.

## A — the attribution of the drowning (the panel, young and old; K2; a constant lowstand of −134 m)

- **A-phys**: the lowstand at the physics level only; the amplification levels erode with the sea at 0 m.
- **A-ampl**: the physics level with the sea at 0 m; the amplification levels erode raised by Δ (F172's chain).
- **A-both**: F172's mechanism. **None**: no lowstand (F171's K2, recomputed).
- Measured for each: the coast (its 4 measures), the land lost, the classes' move.
- **The share of the excess due to the amplifications** = (LR(A-ampl) − LR(none)) / (LR(A-both) − LR(none)), and the same for the
  dimension, per world.

## B — the sea-level curve

### B1 — the curve

- **The series**: Spratt, R. M., Lisiecki, L. E. (2016). A Late Pleistocene sea level stack. Climate of the Past 12, 1079–1092.
  doi:10.5194/cp-12-1079-2016.
  - **The column `SeaLev_longPC1`**: « Scaled first principal component of five sea level reconstructions (0-798 ka) », in m
    above present, every 1 kyr, 0–798 ka.
  - In the paper: the abstract (p. 1079: « five records from 0 to 798 ka », « mean uncertainty estimates of 9–12 m (1σ) »),
    Fig. 2 (p. 1085, the long and short stacks) and Table 1 (p. 1084, the PCA loadings).
  - The numbers are not in the PDF: « The sea level stack is archived in the Supplement and at the World Data Center for
    Paleoclimatology » (Data availability, p. 1090).
  - **The file used**: NOAA's archive, `https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/spratt2016/spratt2016.txt`
    (study 19982, contributed 2016-04-12), downloaded 2026-10-12. It is kept in `data/sealevel/` (git-ignored) and parsed
    (no digitisation).
- **Taken as published.** At 0 ka the stack reads +8.96 m (within its ±5.7 m 1σ). It is not offset; at the very end the sea is
  at 0 m by definition.
- **The series' figures**: min −130 m; the mean over 0–798 ka is −52.9 m.
- **Before 798 ka: 0 m (the main choice).**
- **The sensitivity**, measured only: the 0–798 ka mean (−52.9 m) applied from 2.58 Ma to 798 ka.

### B2 — where it acts

- **Only at the physics level**: the base level follows s(t), linearly interpolated at each sub-step's mid-age.
- The amplification levels erode with the sea at 0 m.
- **The sensitivity**, measured only: the chain raised by 52.9 m (the stack's mean).

### B3 — the time step

- C1's steps keep their durations (linear time, or F169's age mapping).
- **A step that overlaps the last 798 kyr** (its later end younger than 798 ka) is cut into n = ⌈d / 2 kyr⌉ equal sub-steps; the
  others are not cut.
- **The uplift** of each C1 step is spread uniformly over its sub-steps; the sub-step's dt is d / n.
- **The diffusion's sub-steps** follow F169's rule at the sub-step's weight.
- The cost is reported.
- With d_N = T/300 = 18.7 kyr, the last 798 kyr are ~43 steps of 10 sub-steps: ~390 more physics steps.

### B4 — the controls (the old témoin)

- s ≡ 0 through the constant path (Δ = 0) is F171 bit for bit: the existing permanent test.
- **The constant lowstand is F172 bit for bit**: its path is unchanged. A permanent test compares it with the same call made
  through the new forcing entry point.
- **The curve path with s ≡ 0 is NOT bit for bit with F171**, because of its 2 kyr sub-steps. The difference is measured and
  reported (max |Δz|, the classes).

### B5 — viz

« Mer : aucune / bas niveau constant (F172) / courbe quaternaire » in the Cascade window.

## C — the origin of the rampart, and a wall-based index

### C1 — the coastal prism's history (the panel, 64², one C1 run per seed)

- **The cells**: today's land within 2 cells of the coast.
- **The data**: each one's Σ Δs_DS per step, attributed by its boundary at that step.
  - **C-O** = a convergent cell within Chebyshev 2 has an oceanic plate cell 4-adjacent.
  - **C-C** = a convergent cell within 2 has only continental cells around it.
  - **None** = no convergent cell within 2 (divergent, transform or internal).
  - The boundary comes from C1's `classify_boundaries`, retargeted, at that step.
  - Also attributed by the plate count at that step: ≥ 10, 7–9, ≤ 6.
- **The question**: the share added while the coast was a convergent boundary (C-O + C-C).

### C2 — the wall-based rampart index (computable at 256²)

- **W₁₀** = the share of the land within 10 km of the coast whose 256² slope (`physio::slope_deg`) is ≥ 8°.
- **W_in** = the same share beyond 20 km.
- **A rampart_W** = W₁₀ ≥ 0.30 and W₁₀ ≥ 2 W_in.
- **The validation**: against the 2 048² walls, the share of the land within 10 km of the coast with slope > 28°, over the panel's
  20 no-sea K2 worlds (5 seeds × 4 ages). Pearson and Spearman are reported.

### C3 — the screening

- W₁₀, W_in and rampart_W are recomputed on the 33 screened seeds, young and old, with F172's screening configuration (K2, the
  constant lowstand).
- `catalogue.tsv` gains the columns `W10`, `W_in`, `rampart_W`.

## D — the measures of the curve

- **The grid**: the panel × K2 × B0, young, intermediate, old; physics 256², the chain to 2 048², the torus.
- **The variants**: none, constant (F172), curve.
- **The coast** against F172's bands: dimension 1.06–1.18, length ratio 1.15–1.60, aligned 0.22–0.32, block index 0.58–0.98.
- **The drowned valleys** (F172's detector, with « the lowstand ocean » = the final ocean cells ≤ −134 m, kept for comparability);
  the islands.
- The land lost and the classes' move against « none ».
- The shelf, the lakes, the texture and the regularity.
- **Images**:
  - the whole maps none | constant | curve;
  - the 25 km crops Corsica | none | constant | curve;
  - s(t) with the sub-steps.
- The cost.

## R — the stop rules (the brief's)

- **The series not found**: do A and C, stop and ask. It was found and is used.
- **The curve's success**, on ≥ 4 of 5 seeds and ≥ 2 ages:
  - the coast inside ≥ 3 of the 4 bands;
  - the land lost < 5 %;
  - the classes move < 3 points.
  - Then the curve enters the continent profile.
- **Still too drowned** (length ratio > 1.60 on ≥ 4 seeds, at ≥ 1 age, with the curve): stop. The next candidate is deposition.
- **No effect** (the dimension rises < 0.01 everywhere): report the valleys' width at the mouth and stop.
- **Forbidden**: as F165–F172.
