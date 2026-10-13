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

# Parts B and C — continents inside their plates and the thinned margin: the selection's stop rule fires

**Status: B built (an initial condition, behind a flag, C1 untouched), the selection stop rule fired on 5 / 5 seeds, Part C not
run; nothing committed without the feu vert.**
- Declared before any measurement: `f174_declared.md` (Parts B and C), predictions `f174_predictions.md` (Parts B and C).
- Code: `crates/ymir-core/src/cascade/margin.rs` (new), `cascade/profile.rs` (`init_state_margin`,
  `select_continent_default_margin`, `select_continent_plates_margin`; the old functions call them with `Margin::off()`).
- The bench: `crates/ymir-core/tests/f174bc_cascade.rs` — `f174bc_measure` (output `f174bc_bench_output.txt`, stopped after the
  selection), `f174bc_diag` and `f174bc_rim` (`f174bc_diagnostics.txt`), `f174bc_fig` (`images/f174bc_partA_vs_margin_c1.png`).

## What was built (no solver)

- **The rim (B1)**: the cells of a continental plate within w_rim = 2.5 cells of an oceanic plate become oceanic on the **same**
  plate, with C1's own oceanic init values (s 0.2, the age by `init_age_field_ridge_aligned`, out of the cratonic mask).
- **The ramp (B2)**: on the continental cells, s follows Péron-Pinvidic et al. 2013's knots (35 → 30 → 10 → 7 km) over the
  proximal / necking / distal widths 200 / 100 / 150 km ÷ 7.5 = 26.7 / 13.3 / 20 km on the map.
- **B-hétéro (B4)**: s_int × (1 + 0.14 n), σ from Christensen & Mooney 1995's shields and platforms (41.5 ± 5.8 km, Table 3).
- **B3's band**: Harris et al. 2014's passive shelf widths per ocean, 34.9–123.0 km (Table 2, p. 11) ÷ 7.5 = 4.65–16.4 km.
  The article gives **no shelf-break depth** (the break is drawn where the slope increases, p. 6); the ~140 m benchmark is not
  supported by it.
- Nothing integrates anything: the initial state is edited once, then C1 runs unchanged. The constraint « no solver » was not
  touched by any point of B.

**The controls (permanent tests, passing)**:
- `the_margin_off_is_the_profile_bit_for_bit`: the margin off, and a margin of zero width, give the profile's initial state
  bit for bit (s, age, types, plate_id, cratonic mask).
- `the_rim_is_the_plates_own_ocean_and_the_ramp_rises`: plate_id unchanged, the rim at s 0.2, the ramp monotonic from 0 to 1
  through the 10 km knot, the coast's ocean ≥ 90 % its own plate's. **Recorded**: where two continental plates meet near the
  coast, one plate's rim can touch the other's continent (3–8 % of the perimeter at t₀, `f174bc_rim`): those few contacts are
  ocean of another plate.
- `the_selection_with_the_margin_off_is_the_default` (ignored, 106 s, run): the margined selection with `Margin::off()` is
  `select_continent_default` bit for bit (seed 8).

## The stop: the selection fails on 5 / 5 seeds

| seed | Part A (off) | margin | w_rim ×0.5 | w_rim ×2 | w_n ×2 | w_n ÷2 | B-hétéro |
|---|---|---|---|---|---|---|---|
| 8 | N 16, k 7, land 50.7 %, r* 3 | none | N 24, k 17, 50.6 %, r* 1 | none | none | none | none |
| 18 | N 16, k 11, 50.4 %, r* 2 | none | none | none | none | none | none |
| 19 | N 20, k 12, 56.2 %, r* 1 | none | none | none | none | none | none |
| 23 | N 16, k 8, 52.3 %, r* 1 | none | none | none | none | none | none |
| 29 | N 20, k 13, 52.3 %, r* 2 | none | none | none | none | none | none |

**The stop rule « the selection fails on ≥ 2 seeds » fired (5 / 5 with the main margin).** Part C was not run; the viz (B6) was
not built.

**Why** (`f174bc_diag`, `f174bc_rim`, `images/f174bc_partA_vs_margin_c1.png`):
- **The margin moves the coast ~9.6 cells of 64² (60 km) inward of the plate boundary**: the rim (2.5 cells) plus the distal,
  necking and part of the proximal domain, under sea level by Airy (the coast at s 0.915 lies 44 km in from the continent's
  edge). On Part A's own continents it halves the land: **50.4–56.2 % → 23.8–26.8 %** at C1's end.
- To get back to 50 % the selection must make 70–75 % of the plates continental, and those continents **wrap the torus**: on
  seed 8 every in-band candidate wraps (16 plates: 6 in band, 0 circumnavigable; 20: 17 and 0; 24: 26 and 0). Part A had 65 in
  band and 2 qualifying at 16 plates.
- **The scale conflict**: at 16–24 plates on 400 km a plate is 80–100 km across; Péron-Pinvidic's margin at ÷ 7.5 is 60 km plus
  the rim. The margin is as wide as a plate. Halving the ramp (w_n ÷ 2) or the rim (w_rim × 0.5) is not enough on 4 / 5 seeds.

**The rim is not consumed** (`f174bc_rim`, the main margin applied to Part A's selected continents, C1 only):
- the perimeter whose ocean belongs to another plate: 3.1–7.9 % at t₀ (the junction contacts) → **0–4.8 % at the end**;
- the active perimeter (a convergent boundary in the 3 × 3): 7.3–20.2 % at t₀ (those same contacts) → **0–10.9 % at the
  end**, never rising (Part A: 48.5–66.4 % → 5.5–18.4 %).
- As read in the code before measuring: `apply_subduction_step` consumes an oceanic convergent cell only next to a
  **continental cell of another plate**; an O-O boundary in front of a rim is « not handled by this closure ». So even with a
  selected continent, the success criterion « 10–60 % active at the end » could not be met by the existing closures.

## Predictions (Parts B and C)

- **The reviewer's**: **Q6 refuted** (the selection finds a continent on 0 / 5, not 5 / 5). Q1–Q5, Q7, Q8 not tested (the stop).
  Q2 is moreover out of reach of the existing subduction closure (the rim check above).
- **Mine**: held: the rim not consumed (0–4.8 % at the end), the ×2 variants failing the selection on ≥ 2 seeds (5 / 5);
  refuted: the selection found on 4–5 seeds (0 / 5), the selection's extra cost (30–120 s; the seven selections of a seed took
  1 500–1 800 s of process CPU together, run in parallel threads, so ~210–260 s each); the rest not tested.

## No solver was needed; what the author would have to decide (nothing chosen)

None of B's points required a solver; the stop is geometric. To go on, one of these (each an initial condition or a design
choice, none a solver):
1. **The margin outward**: annex to the continental plate the band of the neighbouring oceanic plates within the margin's
   width, so the coast stays near Part A's and the land is kept. An initial assignment of `plate_id` (not an advection), but
   it moves the plate boundaries ~10 cells offshore — more than half the width of an oceanic plate at 16 plates.
2. **A narrower margin** than the literature ÷ 7.5 (÷ 2 still fails; ÷ 4 or less): a departure from the declared scale.
3. **A land band for the margin profile** (e.g. 25–62 %): the band is a design target (F168–F169), not physics.
4. **Fewer, larger plates** (a plate ~140 km at 8 plates): the margin becomes relatively narrower; the 16-plate default came
   from F170–F171's selection.
- And independently: **the rim is not consumed by the subduction closure as it is**. Activity at a coast with a rim would need
  the closure to treat an O-O convergent boundary whose upper plate carries continent (consume the rim, then promote) — a
  change of a C1 closure behind the flag, not an initial condition; left to the author.

## Cost

- The selection stage of the bench: 1 382 s wall for 5 seeds × 7 variants (the variants of a seed in parallel threads).
- `f174bc_rim` / `f174bc_fig`: C1 only, < 1 s per run.

## At commit of Part B (2026-10-12)
- **The reviewer's errors, recorded**:
  - the ramp was placed inward (it moved the coast ~60 km inland);
  - the whole proximal domain was included on the land side;
  - the ~140 m shelf-break benchmark is not supported by Harris et al. 2014;
  - the rim's consumption by the subduction closure had not been checked in the code.
- **The author's decisions, verbatim**:
  1. « Oui c'est bien d'avoir un peu de terre dans l'océan. » Option 1 is retained: the margin seaward.
  2. « Ca demanderait de faire de la biblio pour ça + calage sur base des publications et expériences. Mais pourquoi pas. »
     Option 5 (an O-O subduction closure) goes through the bibliography first; it is not coded in F174.
  - Options 2 (a narrower margin), 3 (a land band for the margin profile) and 4 (fewer, larger plates) are discarded.
  - The author's instruction stands: no solver in the mechanisms without a real need.
- The code of Part B is committed as is (`cascade::margin`, off by default, its permanent controls); the selection's stop is
  recorded.

# Parts B′, C′ and D — the margin seaward, and the bibliographic note for option 5

**Status: B′ built (an initial condition, no closure modified); the selection holds (5 / 5) but the success fails on 0 / 5 seeds
(no shelf); one point is stopped because it needs a modified closure; Part D written; nothing committed without the feu vert.**
- Declared before any measurement: `f174_declared.md` (B′ and C′, with two amendments after the smoke run), predictions
  `f174_predictions.md` (B′ and C′).
- Code: `cascade::margin` (`Margin::seaward`, `apply_seaward`, `seaward_conflicts`; Part B's code unchanged); the viz
  (« Marge vers le large », the « Bathymétrie côtière » view, the « Plaques à t₀ (bande) » map in the terms' frieze).
- The bench: `f174bc_cascade.rs` — `f174b2_measure` (`f174b2_bench_output.txt`, 5 642 s), diagnostics `f174b2_ds` and
  `f174b2_profile` (`f174bc_diagnostics.txt`), the supplementary `f174b2_shelf_zone` (`f174b2_shelf_zone_output.txt`).
- Part D: `f174_option5_note.md`.

## B′ — what was built

- The ramp seaward: 32 km of crust (Airy's 0 m) on Part A's continental edge; 35 → 32 km inside (16 km); the rest of the
  proximal (to 30 km at +10.7 km), the necking (to 10 km at +24 km) and the distal domain (to 7 km at +44 km) on a band of the
  ocean **attached to the continental plate** at t₀ — continental while the crust exceeds 10 km (Péron-Pinvidic's taper break,
  p. 25), oceanic beyond. The nearest continent wins where bands meet.
- **The permanent control** `the_seaward_margin_keeps_the_coast_and_annexes_a_band`: a zero-width band is Part A bit for bit;
  plate_id changes only on band cells, each to its nearest continent's plate; the type rule; 32 km on the edge; the first sea
  cell under s = 0.915. The flag off is Part A bit for bit (Part B's permanent controls).
- **At t₀ the profile is the declared one** (`f174b2_profile`, mean Airy altitude by the distance to Part A's edge, 64² cells):
  +250…+460 m at −1…−2 cells, −160 m at +1, −900 m at +2, −2 700 m at +3, −4 500 m at +4, −5 500 m from +5.

## C′ — the measures (the panel, K2, four ages; each variant with its own selection)

**The selection** (R1, R6):

| seed | Part A: N, k, C1 end land | seaward: N, k, C1 end land | \|Δ land\| | B-hétéro: N, k, land |
|---|---|---|---|---|
| 8 | 16, 7, 50.7 % | 16, 8, 52.2 % | 1.5 | 16, 8, 51.1 % |
| 18 | 16, 11, 50.4 % | 16, 8, 56.1 % | 5.7 | 16, 11, 52.0 % |
| 19 | 20, 12, 56.2 % | 16, 7, 55.6 % | 0.6 | 16, 7, 53.5 % |
| 23 | 16, 8, 52.3 % | 16, 8, 50.9 % | 1.4 | 16, 6, 51.1 % |
| 29 | 20, 13, 52.3 % | 16, 9, 54.5 % | 2.2 | 16, 8, 54.2 % |

- Found on 5 / 5 with the seaward margin and with B-hétéro; the same N as Part A on 3 / 5; C1's end land within ±5 points on
  4 / 5. **The chosen continent differs on every seed** (k or start), so Part A and the margin are not the same continent.

**The conflicts (B′2)**, on the seaward continents: (a) and (b) 0 % of the perimeter on every seed (on Part A's own
continents: 0 / 0 / 28.6 / 34.8 / 0 % and 0 / 0 / 12.6 / 16.4 / 0 %); (c) oceanic plates fully annexed 0 / 2 / 1 / 0 / 0, losing
> 50 % 4 / 8 / 8 / 6 / 7, and 64–100 % of the perimeter faces such a plate. C1 ran with the annexed plates without error.

**The margins** (C1, every step): active at t₀ 3.8–10.4 % (the junction contacts: 2.0–3.8 % of the perimeter faces another
plate's ocean), **at the end 0 / 1.9 / 0 / 0 / 0 %** (Part A: 5.5–18.4 %). The margins stay passive.

**The rampart** (W₁₀ at young): Part A 0.591 / 0.532 / 0.552 / 0.629 / 0.434 → seaward 0.695 / 0.051 / 0.202 / 0.111 / 0.097, a
drop of −18 / 90 / 63 / 82 / 78 % (≥ 40 % and ≥ 60 % on 4 / 5); rampart_W still fires on seed 8 only.

**The shelf** (B3's declared measure, at young): **the median width is 0 km on every seed and every variant, Part A included**;
the share of coastline cells with a detected shelf 2 / 28 / 23 / 23 / 23 % (Part A 16–23 %), the break 91–275 m. **R3 refuted;
the success fails on its shelf criterion on 5 / 5.**
- **Supplementary, not declared** (`f174b2_shelf_zone`, F173's zone measure at young): the 0 → −140 m zone is 0.2–1.1 km wide
  (median; Part A 0.4 km), the 0 → −448 m zone 0.8–2.1 km (Part A 1.0–1.1 km); the mean depth 1 km off the coast −560 to
  −2 480 m (Part A −870 to −1 680 m). The seaward margin widens the shallow zone on seeds 18, 19, 29 by ×2–3, but it stays
  ~1–2 km wide, not the declared 10.7 km.

**Why there is no shelf** (`f174b2_profile`, `f174b2_ds`): over C1's history **the band's continental cells are thickened and
emerge**. At the end the Airy altitude is +1.5 to +3.2 km out to +3 cells on seeds 8 and 23 (+0.5 km to +2 cells on seed 18):
the coast moves seaward onto the band, to the taper break (the continental / oceanic type contact at +24 km), where Airy reads
the oceanic cells at s = 0.2 — the step is back, one cell from land to −5 km. Davis-Suppe is what lifts the band: on Part A's
continent its crust falls by 59–73 % with the margin (its coastal share 45–64 % → 17–29 %), while off it — on the band — it
rises (seed 8: 421 → 1 837). The plate boundary is now ocean against ocean 44 km offshore; `retarget_upper_plate_continental`
does not touch O-O boundaries, so Davis-Suppe seeds there with the velocity-based upper plate and its fallback « dome »
geometry, and its wedge reaches the band's continental cells. The selection compensates the larger land by picking smaller
continents (fewer plates), which is why the land stays within ±5 points.

**The coast**:
- A1 at B0 (the share within 1 cell of C1's contour): Part A 97.2–99.5 % → seaward 91.6 / 51.1 / 76.6 / 91.4 / 60.9 % — a drop
  of 7.6 / 46.1 / 20.9 / 6.7 / 38.6 points (≥ 20 on 3 / 5).
- The whole coast's dimension at young: 0.983 / 0.988 / 0.975 / 0.974 / 0.977 against Part A 0.992 / 0.974 / 0.980 / 0.997 / 0.999
  (Δ −0.023 to +0.014): no rise; under the rocky band. **R5 refuted.**
- By the rock at the coastline (young): the seaward coasts are mostly cover (low) windows (9–31 per world against Part A's 0–4),
  their median dimension 0.958–0.991, inside the sandy reference (0.95–1.00); the few rocky windows never reach the rocky band.
- The seaside slope p50 at young: 13.9 / 1.1 / 1.8 / 2.3 / 1.3° against Part A 5.6–7.4° (gentle coasts, except seed 8).
- The dendritic arms: fewer (0 on seeds 8, 19, 23 at every age; 0–4 elsewhere) than Part A (0–8).
- The islands ≥ 1 km² off the main land: 0 on seeds 8, 19, 23; 1–3 on seed 18; 10–14 (~90 km²) on seed 29 at intermediate and
  old. Not more than Part A overall.

**The macro** (no rule): the mountain share at young 89.7 / 29.6 / 58.4 / 62.6 / 48.6 % against Part A 61.1 / 56.4 / 53.0 /
59.7 / 57.4 % (Δ +28.6 / −26.8 / +5.4 / +2.9 / −8.8 points; R7's −10 points on 1 / 5). Seed 8 becomes a high plateau (90 %
mountain, the first decile 759 m) — the band's and the coast's thickening; the others are lower with fewer coastal walls.

**Statistical convergence** (seeds 8 and 23, B0 and old, 300 vs 600 steps): drainage density −1.6 … +0.1 %, λ −2.6 … −0.5 %,
S₁₀ −0.4 … +0.6 % — **R8 holds 4 / 4**; the classes move 0.05–0.73 points except seed 8 old (3.66 points).

**B-hétéro** (young): the dimension 0.975 / 0.984 / 0.982 / 0.978 / 0.983 (Δ against the seaward margin −0.008 … +0.009); it
selects other continents on three seeds, so its macro moves are not separable (the mountain share 83.7 / 38.5 / 58.0 / 37.4 /
70.2 %).

**Lakes**: fewer on seeds 18, 23, 29 (0–454 km²), more on seed 8 (7 092 km²). **Texture**: the 2 048² slope p50 at young
1.3–20.3°. **Cost**: 790–940 s CPU per world with the selection (the selection's CPU is process-wide, its variants run in
parallel threads); 1–3 s of measures.

**Images** (north up, checked): `f174b2_maps_young.png` (rows: seeds; columns: Part A | seaward | B-hétéro — the coastal walls
of Part A, the pale shallow fringe and the low coasts of the margin, seed 8's plateau), `f174b2_plates_t0.png` (the band
around each continent), `f174b2_margins_c1_end.png`, `f174b2_coast_crops_25km.png` (Corsica west | the Landes | Part A |
seaward), `f174b2_transect_seed8.png` and `f174b2_transect_seed18.png` (a smooth sigmoid from 0 to −5 km, no flat).

## The success and the stops

- **The selection**: found on 5 / 5 — the stop rule does not fire.
- **The success**: 0 / 5 — W₁₀ −40 % on 4 / 5 and the land ±5 points on 4 / 5, but **the shelf on 0 / 5**. The seaward margin does
  not enter the profile.
- **A point stopped because it needs a modified closure**: the shelf. With initial conditions only, the band either stays
  continental (Davis-Suppe's O-O fallback lifts it and the coast moves to the taper break) or is made oceanic (Airy reads
  oceanic cells at s = 0.2 and the step returns at the old coast). Getting a shelf needs one of:
  1. Davis-Suppe **not seeding on ocean-ocean convergent boundaries** (they are island-arc settings — option 5's ground); a
     change to a C1 closure;
  2. the band's thinned crust kept continental but **excluded from Davis-Suppe's wedge** beyond the taper break; also a closure
     change;
  3. reading the true s of oceanic cells on continental plates in the Airy mapping — a cascade-side change, but C1's oceanic s
     is advected with the plate (the ramp would be carried away) and the rigidity mask (#145) covers only continental cells; a
     C1 change as well.
  None is implemented.

## Predictions (B′ and C′)

- **The reviewer's**: R1 held (5 / 5, the same N on 3); R2 held (≤ 5 % active on 5 / 5; W₁₀ −60 % on 4 / 5); **R3 refuted** (0 / 5);
  R4 on 3 / 5 (7.6–46.1 points); **R5 refuted** (no rise); R6 held (4 / 5); **R7 refuted** (1 / 5); R8 held (4 / 4).
- **Mine**: held: R1, R6, active ≤ 5 %, R8, the conflicts (a)–(b) small and (c) ≥ 1 plate fully annexed on ≥ 2 seeds; refuted: the
  W₁₀ count (I said ≥ 60 % on 2–3; 4), R3 (the shelf), R4 (I said 3–15 points), R5, R7, more islands, unchanged arms, the success
  ≥ 4 (0).

## D — the bibliographic note for option 5

`f174_option5_note.md`: the four articles were all present and read (Stern 2004, Nikolaeva et al. 2010, Crameri et al. 2020,
Stern & Gerya 2018), cited by page, paragraph, figure or table. In short: the older oceanic side sinks; passive-margin
initiation is rare (no Cenozoic example; 2 of the 13 SZI events near and none along a passive margin; only thin, hot margins in
N10's experiments; Cloetingh's 20 Ma rule via SG18); induced initiation needs ~100 km of convergence at ≥ 1–2 cm/yr and takes
5–10 Myr; a continent arriving at a trench jams it and polarity reversal makes an Andean-type margin. The proposed closure
(O-O consumption of the older side, island arcs on the overriding side, the boundary migrating by reassignment, an initiation
after D_init of convergence, no spontaneous passive-margin collapse) comes with calibration targets and six validation
criteria. **Nothing is implemented.** Note: the B′ result adds a second O-O question to option 5 — Davis-Suppe's fallback at
O-O boundaries.

## Open, for the author (nothing chosen)

- (a) The shelf needs a closure change (1–3 above), or option 5's O-O rule (which would also stop Davis-Suppe's fallback dome at
  O-O boundaries if it replaced it there).
- (b) Option 5 itself (the note).
- (c) Whether the seaward margin, which removes the rampart on 4 / 5 seeds and keeps the margins passive, is worth keeping behind
  its flag in the meantime (the shelf criterion fails).

## At commit of Parts B′, C′ and D (2026-10-13)
- **The reviewer's error, recorded**: a continental shelf was expected from an initial condition alone (B′'s success criterion).
  A flat, shallow shelf is mostly a sedimentary prism (the reviewer's reading, a hypothesis to check); the shelf goes to
  deposition (F176).
- **The seaward margin is off by default** (`Margin::off()`; `Margin::seaward()` only behind its flag); it is activated in F175.
- **A method rule**: an initial-condition variant that changes the selection is compared **also at a fixed continent** (the
  plates chosen for the variant, run with and without it), not only with the selection redone — B′ chose a different continent
  on every seed, which mixed the margin's effect with the continent's.
- **The time-step question is closed**: the physics level's step is not halved; the criterion is statistical convergence
  (density, λ, sinuosity, classes), which holds (F174 A3; B′'s R8 within 2.6 % on 4 / 4).
- **The author's decisions, verbatim**: on option 5, « Ca demanderait de faire de la biblio pour ça + calage sur base des
  publications et expériences. Mais pourquoi pas. »; on the order, « oui ok pour l'ordre » (F175 = option 5, F176 = deposition).
  The instruction stands: no solver; option 5 is a change of an existing closure, behind a flag.
