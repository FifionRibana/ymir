# Finding 174 — the coastline inherited from C1, real reference coasts, the time step, the active perimeter

**Status: Part 0 and Part A measured; Part B not started (its sources could not be read); nothing committed without the
feu vert.** Production unchanged; no cache key touched (no ALGO_* bump).

- Declared before any measurement: `f174_declared.md`. Predictions: `f174_predictions.md`.
- The bench: `crates/ymir-core/tests/f174_cascade.rs` (`f174_measure`), output `f174_bench_output.txt`, 1 593 s wall.
- The real windows: `prep_coasts.py` (output `prep_coasts_output.txt`), grids in `data/coasts/` (git-ignored).
- Images: `images/` (north up, checked; the Copernicus notice in `images/NOTICE.md`).
- The configuration: the panel 8, 18, 19, 23, 29; K2; the torus; no sea; ages B0 and old (300 Myr).

## Part 0 — the Spratt re-zero

- `SeaCurve::from_noaa_spratt2016` now subtracts the youngest row's level, so the curve reads 0 m at 0 ka (it read +8.96 m).
  An alignment of reference: the shape and the amplitudes are unchanged.
- The permanent test `the_sea_curve_entry_point_keeps_f172_and_sub_steps` now checks the re-zero (0 today, every level
  shifted by the same constant). F173's measures are not re-run: the sea is out of the profile (F173 (a)).

## A1 — the coastline is C1's, except where old valleys drown

The share of the final 2 048² coastline within 1 / 2 / 5 cells of 256² of C1's s = 0.915 contour (C1's last Airy snapshot,
upsampled as the cascade does); in brackets the physics level alone at 256².

| seed | B0: ≤ 1 / 2 / 5 cells | B0 physics level | old: ≤ 1 / 2 / 5 cells | old p90 (cells) | old physics level |
|---|---|---|---|---|---|
| 8 | 99.2 / 100.0 / 100.0 % | 99.1 % | 88.8 / 93.2 / 95.0 % | 1.12 | 87.5 % |
| 18 | 97.2 / 99.8 / 100.0 % | 98.1 % | 79.6 / 94.4 / 99.7 % | 1.41 | 65.3 % |
| 19 | 97.5 / 99.9 / 100.0 % | 99.0 % | 70.6 / 77.6 / 79.8 % | 34.12 | 66.9 % |
| 23 | 98.1 / 100.0 / 100.0 % | 99.1 % | 86.7 / 95.3 / 97.7 % | 1.25 | 77.0 % |
| 29 | 99.5 / 100.0 / 100.0 % | 99.3 % | 65.9 / 72.3 / 78.4 % | 23.15 | 64.6 % |

- **At B0, P1 holds on 5 / 5 seeds** (97.2–99.5 % within one cell, the median distance 0.25 cell, nothing beyond 2 cells).
  The coastline is the 0 contour of C1's thickness step, as the reviewer read it.
- **At old, P1 fails on 5 / 5** (65.9–88.8 %). P1 is therefore **mixed**, not refuted at both ages: the stop rule « P1
  refuted » does not fire.
- **What moves the coast at old**: the physics level (its share alone is as low or lower). Where the coast leaves the C1
  contour by more than 5 cells (up to 20–22 % of the coast on seeds 19 and 29), it is an arm of the sea reaching up a
  valley system graded to 0 m over 300 Myr (image `a1_seed19_old.png`: a dendritic inlet ~70 km long). Elsewhere the old
  coast stays within 1–2 cells of C1's.
- So the shape of the coast at every scale above ~1.5 km is C1's 64² mask, smoothed; the only structure the cascade adds
  is drowned dendritic valleys at old (the feature the author found « pas naturel » in F172).

## A2 — real coasts by type

### The real windows (50 km, 195.3 m; the window variant of `physio::coast`)

| window | declared | dimension | length ratio | aligned | block index | coast (km) | rule: rocky share → type |
|---|---|---|---|---|---|---|---|
| Galicia, Rías Baixas | rocky | 1.219 | 1.450 | 0.281 | 1.303 | 266 | 37 % → mixed |
| Dalmatia, Split | rocky | 1.223 | 1.558 | 0.293 | 1.229 | 235 | 39 % → mixed |
| North Brittany | rocky | 1.152 | 1.455 | 0.279 | 0.912 | 136 | 0 % → low |
| Corsica west | rocky | 1.152 | 1.430 | 0.272 | 1.158 | 119 | 99 % → rocky |
| Landes | low | 0.952 | 0.837 | 0.109 | 1.493 | 50 | 0 % → low |
| Vendée | low | 1.002 | 0.977 | 0.226 | 0.804 | 82 | 0 % → low |
| West Jutland | low | 1.354 | 2.344 | 0.345 | 0.869 | 273 | 0 % → low |
| Golfe du Lion | low | 1.250 | 1.823 | 0.321 | 1.017 | 171 | 0 % → low |

- **The control**: on `corse_256` (a closed island off the border) the window variant gives `physio::coast`'s four values
  exactly (1.154 / 1.295 / 0.292 / 1.009); asserted in the bench.
- **The bands** (the declared types, [min, max], not widened):
  - **rocky**: dimension 1.152–1.223, length ratio 1.430–1.558, aligned 0.272–0.293, block index 0.912–1.303;
  - **low**: dimension 0.952–1.354, length ratio 0.837–2.344, aligned 0.109–0.345, block index 0.804–1.493.
- **The low coasts are two families at 195 m**: the straight sandy coasts (Landes 0.952, Vendée 1.002) and the barrier /
  lagoon coasts (Jutland 1.354, Golfe du Lion 1.250), whose irregularity is the lagoon's inner shore (Ringkøbing Fjord is
  connected to the sea at 195 m; the étang de Thau is not, and counts as land). The low band is therefore wider than the
  rocky one and contains it in dimension.
- **P2 is refuted on both ends**: the low coasts span 0.95–1.35 (not 1.02–1.08), the rocky ones 1.15–1.22 (not 1.10–1.20;
  Galicia and Dalmatia exceed 1.20).
- **The rule** (p90 of the altitude within 5 km > 200 m) agrees with the declared type on **5 / 8** windows. It calls Galicia
  and Dalmatia mixed (their islands and ría heads are low) and North Brittany low (a rocky coast under 300 m). The brief's
  example rule separates relief, not coast type. **Flagged for the author; not tuned here.**

### Our worlds (64 windows of 256² cut from the 2 048² field; the rule types each window)

| seed | B0: windows rocky / low / mixed | B0 rocky median (dimension, in-band) | old: windows rocky / low / mixed | old low median (dimension, in-band) | old rocky median (in-band) |
|---|---|---|---|---|---|
| 8 | 24 / 0 / 1 | 0.984 (0 / 4) | 8 / 19 / 1 | 0.998 (3 / 4) | 0.995 (0 / 4) |
| 18 | 18 / 0 / 4 | 0.987 (1 / 4) | 5 / 16 / 2 | 1.016 (2 / 4) | 1.005 (1 / 4) |
| 19 | 26 / 1 / 4 | 0.978 (1 / 4) | 5 / 26 / 5 | 1.024 (3 / 4) | 1.009 (0 / 4) |
| 23 | 28 / 1 / 2 | 0.995 (0 / 4) | 5 / 24 / 3 | 1.017 (3 / 4) | 1.014 (0 / 4) |
| 29 | 22 / 0 / 3 | 0.996 (0 / 4) | 1 / 23 / 7 | 1.007 (4 / 4) | 1.001 (0 / 4) |

- At B0 our coasts are rocky by the rule (the rampart) and **under the rocky band on every seed** (dimension 0.98–1.00,
  length ratio 0.89–0.93, 0–1 of 4 measures in band).
- At old they are mostly low, and their medians fall **in the low band on 2–4 of 4** measures — but only because the low band
  reaches down to the Landes: our low coasts are as smooth as the straightest sandy coast (dimension 1.00–1.02, length ratio
  0.95–1.01), never like the lagoon coasts.
- The whole coast (`coast_torus`, F173's measure) is unchanged in character: dimension 0.97–1.02, length ratio 0.88–1.00.

## A3 — the physics level's time step (300 vs 600 steps)

| seed · age | |Δz| all p50 / p99 / max (m) | land p50 / p99 (m) | classes' largest move (points) | cells changed | order-≥ 3 channels displaced |
|---|---|---|---|---|---|
| 8 B0 | 3.4 / 156 / 1 005 | 6.4 / 202 | 0.26 | 1.62 % | 10.7 % |
| 8 old | 2.4 / 190 / 928 | 3.5 / 313 | 0.43 | 1.47 % | 11.8 % |
| 18 B0 | 1.7 / 123 / 1 301 | 3.1 / 187 | 0.12 | 1.20 % | 9.0 % |
| 18 old | 4.1 / 185 / 1 152 | 3.8 / 332 | 0.20 | 1.47 % | 15.3 % |
| 19 B0 | 3.3 / 102 / 912 | 4.3 / 143 | 0.39 | 1.27 % | 7.5 % |
| 19 old | 3.7 / 129 / 635 | 2.2 / 152 | 0.37 | 1.71 % | 16.4 % |
| 23 B0 | 3.3 / 152 / 726 | 4.3 / 196 | 0.75 | 1.73 % | 8.9 % |
| 23 old | 3.1 / 138 / 1 118 | 3.3 / 206 | 0.06 | 1.33 % | 12.5 % |
| 29 B0 | 2.5 / 106 / 501 | 4.1 / 148 | 0.40 | 1.62 % | 7.7 % |
| 29 old | 3.1 / 78 / 678 | 3.3 / 70 | **1.65** | 2.81 % | 16.5 % |

- The step is **not converged locally**: max |Δz| 501–1 301 m on every world, p99 78–190 m, while the median is 2–4 m. The
  macro is close to converged: the classes move < 1 point on 9 of 10 worlds (1.65 on seed 29 old), 1.2–2.8 % of the cells
  change class.
- The network moves: 7.5–9 % of the order-≥ 3 channels at B0, 12–16.5 % at old.
- **P3**: max |Δz| > 200 m on 5 seeds (held); the classes < 1 point everywhere (refuted on one world). P3 is refuted, narrowly.
- Cost: the physics level at 600 steps is 46–62 s CPU against 22–30 s.

## A4 — why every coast was active

### The measures (the perimeter at 64², every C1 step)

| seed | first 50 steps: C / D / T / none | all steps: C / D / T / none | same-plate share: first 50 / all | plates t₀ → end |
|---|---|---|---|---|
| 8 | 65.9 / 23.6 / 9.7 / 0.8 % | 26.8 / 18.3 / 2.0 / 53.0 % | 0.94 / 55.6 % | 16 → 2 |
| 18 | 48.3 / 37.7 / 13.5 / 0.4 % | 18.9 / 27.3 / 3.7 / 50.0 % | 0.57 / 54.3 % | 16 → 3 |
| 19 | 51.6 / 37.2 / 10.9 / 0.3 % | 30.7 / 13.2 / 2.2 / 53.9 % | 0.49 / 56.0 % | 20 → 2 |
| 23 | 53.4 / 43.0 / 3.1 / 0.5 % | 29.2 / 20.5 / 0.7 / 49.7 % | 0.63 / 53.4 % | 16 → 2 |
| 29 | 57.9 / 35.0 / 6.8 / 0.3 % | 32.0 / 17.4 / 2.6 / 48.0 % | 0.38 / 53.9 % | 20 → 4 |

- Over the first 50 steps the perimeter is **48–66 % convergent (the panel's mean 55.4 %)**, 24–43 % divergent, 3–14 %
  transform, and **0.3–0.8 % without a boundary**. Only 0.4–0.9 % of the perimeter has its ocean on the same plate.
- The shares do not change until the first accretion merges (step ~50); then « none » climbs to 48–78 % as the plates merge
  (16 or 20 → 2–4) and put continent and ocean on one plate. By then the coastal prism is built (F173: 37–64 % of it while
  ≥ 10 plates remained).
- **P4**: > 60 % convergent on 1 / 5 seeds and not on the mean (55.4 %): **refuted** on its number. Its cause: see below.

### The code reading

- **Where the zero velocity is coded**: `tectonics_c1/time_loop.rs:251`, `apply_continental_rigidity` (Issue #145, « the
  buoyancy fix »): after `fill_velocity_field`, vx = vy = 0 on every `plate_type == Continental` cell, with a no-flux
  boundary in the upwind sweep (`continental_rigid_mask`). Called at the advection-only path (l. 161) and in the closures
  loop (l. 627, and l. 720 when Track D re-types cells).
- **Why**: buoyant continental crust must not advect into convergence and subduct or funnel (#145). It is also what makes
  C1's Eulerian record a material history (`cascade/history.rs`'s header).
- **What it implies, and what it does not**:
  - The continental s never moves: the continent is fixed in the grid frame for the whole history.
  - **The boundary types do not read the cell velocities.** `classify_boundaries` (`boundary_classification.rs:199`) uses
    `kinematics.velocities[plate_id]`, the plate's velocity, which is not zeroed. A continental plate still « moves » for
    the classification. **The zero velocity is not what makes the coasts active.**
  - **What does**: (i) the continent profile builds continents out of **whole plates** (`per_plate_type` broadcast to every
    cell of the plate, `init_r7/mod.rs:355`), so at t₀ every coast is a plate boundary (same-plate share < 1 %); (ii)
    `plate_id` is never advected (`time_loop.rs`'s « Static-classification » note), so the boundaries stay put and a
    convergent coast stays convergent until a merge; Davis-Suppe keeps thickening the same coast.
  - The reviewer's « immobile continents » is therefore half right: the plates are immobile (static `plate_id`), not because
    of the zero velocity of the continental cells.
- **Proposals, in order of preference (nothing corrected)**:
  1. **Continents inside their plates** (continent profile, behind a flag): the continental type on a plate's interior,
     with an oceanic rim of the **same** plate between the coast and the plate boundary. A passive margin by construction at
     t₀; it becomes active only where a convergent boundary consumes the rim (the subduction closure exists). It combines
     with Part B's ramp (the rim is the thinned margin). C1 untouched outside the profile.
  2. **Continents that move with their plates** (the queued C1 (b)): advect `plate_id` and the continental s with the plate
     velocity; trailing edges become passive by themselves. A large change: it conflicts with #145's reason and with the
     material-frame assumption of the history record, and needs a new no-subduction rule for continental crust.
  3. **Slow continental plates** in the kinematics (continents are slow in a hotspot frame): fewer convergent metres per
     step on continental margins, but the classification's shares unchanged. Partial; listed for completeness.

## Part B — not started

- Péron-Pinvidic et al. 2013 (Marine and Petroleum Geology 43, 21–47) sits behind HAL's bot protection; Harris et al. 2014
  (Marine Geology 352, 4–24) is an Elsevier page that does not render; Christensen & Mooney 1995 (JGR 100, 9761–9788) was
  not fetched. Neither refusal was circumvented, and none of the three is in `docs/refs`.
- So the margin's profile, w_n and B-hétéro's σ cannot be declared from a read source. **The author is asked to add the
  three articles to `docs/refs`.** Part C is not run; P5–P7 stay open.
- A4's proposal 1 is relevant to B: a ramp on a coast that is a plate boundary stays a plate boundary. If the margin is
  built, building it on the continent's own plate (proposal 1) is what lets it stay passive.

## The stop rules

- P1 refuted at both ages: **no** (held at B0, failed at old) → the rule does not fire.
- Part B stops anyway, on its sources.

## Predictions

- **The reviewer's**:
  - P1 held at B0 (5 / 5), refuted at old (0 / 5): mixed.
  - P2 refuted on both ranges.
  - P3 refuted narrowly (the classes move 1.65 points on one world; the Δz part held).
  - P4 refuted on its number (55.4 %, 1 / 5 seeds > 60 %); its cause half right (static plates, not the zero velocity).
  - P5–P7 not tested.
- **Mine**:
  - held: P1 at B0; P1 failing at old in 60–90 % (65.9–88.8); the Landes ≤ 1.02 and length ratio ≤ 1.05; Jutland and the
    Golfe du Lion > 1.08; Galicia and Dalmatia > 1.20; the bands overlapping; our coasts rocky at B0, low at old, as smooth as
    the smoothest low coast; max |Δz| > 200 m on ≥ 3; the channels displaced 5–20 %; « none » ≈ 0 and the same-plate share
    ≈ 0 early; the cause in the plate geometry, not the zero velocity; A3 doubling the physics level;
  - refuted: ≥ 98 % within 5 cells at old (78.4–99.7); the physics level's share above the 2 048² one (it is equal or
    below); the rule right on 7 / 8 (5 / 8); the length ratio separating the types (the lagoons reach 1.8–2.3); p99 < 50 m
    (78–190); the classes < 1 point everywhere; the A3 cost +15–25 % of a world (+~5 %: the chain dominates);
  - not predicted: the late merges making half of the perimeter passive.

## Cost (CPU seconds per world)

- C1 64² < 1 s; the physics level 22–30 s (46–62 s at 600 steps); the chain to 2 048² 453–512 s.
- A1 and A2 add measures only.

## Open, for the author (nothing chosen)

- (a) The three Part B articles in `docs/refs` (Péron-Pinvidic 2013, Harris 2014, Christensen & Mooney 1995).
- (b) The coast-type rule: the brief's example fails on 3 / 8 real windows. Keep it, change it, or type the coasts another way.
- (c) The low reference: two families (straight sandy, lagoon / barrier); the lagoon family needs deposition (F175).
- (d) A4's proposals: continents inside their plates (1), moving continents (2), slow continental plates (3).
- (e) The physics level's step: locally unconverged (max |Δz| up to 1.3 km, 12–16 % of the large channels displaced at old),
  macro converged; whether to halve it (+25 s per world).

## At commit of Part A (2026-10-12)
- **The reviewer's errors, recorded**:
  - the coast-type rule « p90 of the altitude within 5 km > 200 m » measured the relief, not the coast type;
  - P2's ranges were wrong at both ends;
  - the rampart's cause was half right: the plates are immobile because `plate_id` is static, not because of the zero
    velocity of the continental cells.
- **The coast-type rule, replaced**:
  - on our worlds, a coast segment is typed by the rock at the coastline (basement or shield: rocky; cover: low), with the
    seaside slope reported;
  - the real windows keep their declared type;
  - the references: rocky 1.15–1.22 (dimension) / 1.43–1.56 (length ratio); low sandy 0.95–1.00; low lagoon 1.25–1.35. The
    lagoon family waits for deposition (F175).
- **The time step is not halved.** The criterion is statistical convergence (classes, hypsometry, network statistics), not
  local convergence, which a drainage network does not have.
- **A4**: proposal 1 is retained — continents inside their plates.
- **The reviewer's reading, a hypothesis to check**: the old profile's dendritic sea arms are valley floors worn to a zero
  slope for want of sediment transport. They belong to deposition (F175).
- **The author's instruction for Parts B and C, verbatim**: « attention à ne pas retomber dans le biais d'intégrer un solveur
  dans les mécanismes sauf si réellement nécessaire. Mais c'est justement ce qu'on a essayé de fuir avec C1! »
  - So C1 stays closure-based: Parts B and C use initial conditions (continents inside their plates, the rim w_rim, a margin
    profile prescribed from Péron-Pinvidic et al. 2013, the B-hétéro field) and the existing closures (the rim consumed by
    the subduction closure). Without the author's consent: no advection of `plate_id`, no simulated rifting or stretching, no
    new solver (flexure, creep, PDE) in C1 or in the initial condition. A point that cannot be done without a solver is
    stopped and explained, not implemented.
- The three Part B sources are now in `docs/refs` (`peron-pinvidic2013.pdf`, `harris2014.pdf`, `christensen1995.pdf`).
